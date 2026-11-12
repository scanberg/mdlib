#include <md_system.h>
#include "md_attributes_internal.h"
#include "md_run.h"
#include <md_nonbonded.h>
#include "md_system_fields.h"

#include <inttypes.h>
#include <stdio.h>

#include <core/md_log.h>
#include <core/md_format.h>
#include <core/md_simd.h>
#include <core/md_array.h>
#include <core/md_hash.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_file.h>

#ifdef __cplusplus
extern "C" {
#endif

// One counter for every system, so that a version is never handed out twice: a consumer which kept the version of a
// system that has since been freed and loaded anew cannot take the new topology for the one it saw.
static volatile uint64_t topology_version_counter = 0;

uint64_t md_system_topology_changed(md_system_t* sys) {
    MD_ASSERT(sys);
    sys->topology_version = md_atomic_fetch_add64(&topology_version_counter, 1) + 1;
    return sys->topology_version;
}

// ### ASSEMBLIES ###

void md_system_assemblies_free(md_system_t* sys) {
    MD_ASSERT(sys);
    md_allocator_i* alloc = sys->alloc;
    md_assembly_data_t* a = &sys->assembly;
    md_array_free(a->id, alloc);
    md_array_free(a->group_offset, alloc);
    md_array_free(a->group_inst_offset, alloc);
    md_array_free(a->inst, alloc);
    md_array_free(a->group_oper_offset, alloc);
    md_array_free(a->oper, alloc);
    md_array_free(a->oper_id, alloc);
    md_array_free(a->oper_transform, alloc);
    md_memset(a, 0, sizeof(md_assembly_data_t));
}

// ### GROUPS ###

static const char* const group_role_mdp[MD_GROUP_ROLE_COUNT] = {
    "tc-grps", "energygrps", "acc-grps", "freezegrps", "user1-grps", "user2-grps", "comm-grps", "compressed-x-grps", "orire-fitgrp", "qmmm-grps",
};
static const char* const group_role_names[MD_GROUP_ROLE_COUNT] = {
    "tc", "energy", "acc", "freeze", "user1", "user2", "comm", "compressed-x", "orire-fit", "qmmm",
};

const char* md_group_role_mdp_option(md_group_role_t role) {
    return (uint32_t)role < MD_GROUP_ROLE_COUNT ? group_role_mdp[role] : NULL;
}

const char* md_group_role_name(md_group_role_t role) {
    return (uint32_t)role < MD_GROUP_ROLE_COUNT ? group_role_names[role] : NULL;
}

void md_system_groups_free(md_system_t* sys) {
    MD_ASSERT(sys);
    md_allocator_i* alloc = sys->alloc;
    md_group_data_t* g = &sys->group;
    for (size_t i = 0; i < g->count; ++i) {
        if (!md_str_empty(g->name[i])) md_free_str(alloc, g->name[i]);
        md_bitset_free(&g->atoms[i], alloc);
    }
    md_array_free(g->name, alloc);
    md_array_free(g->atoms, alloc);
    md_array_free(g->roles, alloc);
    md_array_free(g->role_group, alloc);
    md_memset(g, 0, sizeof(md_group_data_t));
}

static void groups_copy(md_group_data_t* dst, const md_group_data_t* src, md_allocator_i* alloc) {
    md_memset(dst, 0, sizeof(md_group_data_t));
    if (!src->count) return;
    dst->count = src->count;
    md_array_resize(dst->name, src->count, alloc);
    md_array_resize(dst->atoms, src->count, alloc);
    md_array_push_array(dst->roles, src->roles, src->count, alloc);
    for (size_t i = 0; i < src->count; ++i) {
        dst->name[i]  = md_str_empty(src->name[i]) ? (md_str_t){0} : md_alloc_str(alloc, src->name[i]);
        dst->atoms[i] = md_bitset_copy(src->atoms[i], alloc);
    }
    md_memcpy(dst->role_offset, src->role_offset, sizeof(dst->role_offset));
    if (md_array_size(src->role_group)) md_array_push_array(dst->role_group, src->role_group, md_array_size(src->role_group), alloc);
}

size_t md_system_group_add(md_system_t* sys, md_str_t name, md_bitset_t atoms, uint32_t roles) {
    MD_ASSERT(sys);
    MD_ASSERT(sys->alloc);
    md_group_data_t* g = &sys->group;
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    const uint32_t N = (uint32_t)sys->atom.count;
    if (md_bitset_last_end(atoms) > N) {
        atoms = md_bitset_and(atoms, md_bitset_range(0, N), md_temp_allocator(temp));
    }
    size_t idx = g->count;
    for (size_t i = 0; i < g->count; ++i) {
        if (md_str_eq(g->name[i], name) && md_bitset_equal(g->atoms[i], atoms)) {
            idx = i;
            break;
        }
    }
    if (idx == g->count) {
        md_array_push(g->name, md_str_empty(name) ? (md_str_t){0} : md_alloc_str(sys->alloc, name), sys->alloc);
        md_array_push(g->atoms, md_bitset_copy(atoms, sys->alloc), sys->alloc);
        md_array_push(g->roles, 0, sys->alloc);
        g->count += 1;
    }
    g->roles[idx] |= roles;
    md_temp_end(temp);
    md_system_topology_changed(sys);
    return idx;
}

void md_system_group_set_role(md_system_t* sys, md_group_role_t role, const uint32_t* group_idx, size_t count) {
    MD_ASSERT(sys);
    if ((uint32_t)role >= MD_GROUP_ROLE_COUNT) return;
    md_group_data_t* g = &sys->group;
    const uint32_t beg = g->role_offset[role];
    const uint32_t end = g->role_offset[role + 1];
    // The role's previous groups lose it, unless named again below
    for (uint32_t i = beg; i < end; ++i) g->roles[g->role_group[i]] &= ~(1u << role);
    uint32_t n = 0;
    for (size_t i = 0; i < count; ++i) n += group_idx[i] < g->count;
    // Make room for the new list in place of the old one, the roles after it moved along
    const size_t total = md_array_size(g->role_group);
    const size_t new_total = total - (end - beg) + n;
    if (new_total > total) md_array_resize(g->role_group, new_total, sys->alloc);
    md_memmove(g->role_group + beg + n, g->role_group + end, (total - end) * sizeof(uint32_t));
    md_array_shrink(g->role_group, new_total);
    uint32_t k = beg;
    for (size_t i = 0; i < count; ++i) {
        if (group_idx[i] >= g->count) continue;
        g->role_group[k++] = group_idx[i];
        g->roles[group_idx[i]] |= 1u << role;
    }
    for (uint32_t r = role + 1; r <= MD_GROUP_ROLE_COUNT; ++r) {
        g->role_offset[r] = (uint32_t)(g->role_offset[r] - (end - beg) + n);
    }
    md_system_topology_changed(sys);
}

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

static int compare_i32(const void* a, const void* b) {
    const int32_t x = *(const int32_t*)a;
    const int32_t y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

// The operators of an assembly named by the same rows; those with the same instances are then merged into one group
typedef struct {
    md_array(uint32_t)          rows;   // The rows which name the operators (those of the first, while grouping by rows)
    md_array(uint32_t)          oper;   // In the order they first appear
    md_array(md_instance_idx_t) inst;   // Ascending, unique
} assembly_group_t;

bool md_system_assemblies_init(md_system_t* sys, const md_assembly_gen_t* gen, size_t gen_count, const md_f32x4x4_t* oper_transform, const md_label_t* oper_id, size_t oper_count) {
    MD_ASSERT(sys);
    md_system_assemblies_free(sys);
    if (gen_count == 0 || oper_count == 0) {
        return true;
    }
    if (!gen || !oper_transform) {
        MD_LOG_ERROR("Assemblies: missing rows or operators");
        return false;
    }
    if (oper_count >= UINT32_MAX) {
        MD_LOG_ERROR("Assemblies: too many operators");
        return false;
    }

    md_allocator_i* alloc = sys->alloc;
    md_assembly_data_t* A = &sys->assembly;
    const size_t inst_count = sys->instance.count;

    // The operators, all of them: the transforms of the coordinates. A rotation which is not orthonormal (a scale, a
    // shear, digits lost) would be drawn wrongly, its normals and radii are not transformed as positions are.
    md_array_push_array(A->oper_transform, oper_transform, oper_count, alloc);
    for (size_t i = 0; i < oper_count; ++i) {
        md_label_t id = {0};
        if (oper_id) {
            id = oper_id[i];
        } else {
            char buf[16];
            const int len = snprintf(buf, sizeof(buf), "%zu", i + 1);
            id = md_label_make((md_str_t){buf, (size_t)MD_MAX(len, 0)});
        }
        md_array_push(A->oper_id, id, alloc);

        const md_f32x4x4_t M = oper_transform[i];
        float dev = 0.0f;
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                const float d = M.elem[r][0] * M.elem[c][0] + M.elem[r][1] * M.elem[c][1] + M.elem[r][2] * M.elem[c][2];
                dev = MD_MAX(dev, fabsf(d - (r == c ? 1.0f : 0.0f)));
            }
        }
        if (dev > 1.0e-3f || M.elem[0][3] != 0.0f || M.elem[1][3] != 0.0f || M.elem[2][3] != 0.0f || M.elem[3][3] != 1.0f) {
            MD_LOG_INFO("Assemblies: operator '" MD_STR_FMT "' is not rigid (off by %g)", MD_STR_ARG(MD_LABEL_STR(id)), dev);
        }
    }
    A->oper_count = oper_count;

    md_array_push(A->group_offset, 0, alloc);
    md_array_push(A->group_inst_offset, 0, alloc);
    md_array_push(A->group_oper_offset, 0, alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* ta = md_temp_allocator(temp);

    bool*     row_done = md_temp_alloc_zero(temp, gen_count * sizeof(bool));
    uint32_t* oper_row_beg = md_temp_alloc(temp, oper_count * sizeof(uint32_t));    // The (operator, row) pairs of each operator
    uint32_t* oper_row_end = md_temp_alloc(temp, oper_count * sizeof(uint32_t));

    for (size_t g0 = 0; g0 < gen_count; ++g0) {
        if (row_done[g0]) continue;
        const md_label_t asm_id = gen[g0].assembly_id;

        // The (operator, row) pairs of the rows of the assembly, and the operators in the order they first appear
        md_array(uint64_t) pairs = 0;
        md_array(uint32_t) opers = 0;
        for (size_t g = g0; g < gen_count; ++g) {
            if (row_done[g] || !md_label_eq(gen[g].assembly_id, asm_id)) continue;
            row_done[g] = true;
            for (size_t k = 0; k < gen[g].oper_count; ++k) {
                const uint32_t o = gen[g].oper[k];
                if (o >= oper_count) {
                    MD_LOG_ERROR("Assemblies: assembly '" MD_STR_FMT "' names operator %u of %zu", MD_STR_ARG(MD_LABEL_STR(asm_id)), o, oper_count);
                    continue;
                }
                md_array_push(pairs, ((uint64_t)o << 32) | (uint64_t)g, ta);
                md_array_push(opers, o, ta);
            }
        }
        const size_t pair_count = md_array_size(pairs);
        if (pair_count == 0) continue;
        qsort(pairs, pair_count, sizeof(uint64_t), compare_u64);

        for (size_t i = 0; i < pair_count; ) {
            const uint32_t o = (uint32_t)(pairs[i] >> 32);
            size_t j = i;
            while (j < pair_count && (uint32_t)(pairs[j] >> 32) == o) ++j;
            oper_row_beg[o] = (uint32_t)i;
            oper_row_end[o] = (uint32_t)j;
            i = j;
        }

        // The operators by the rows which name them: the same rows, the same instances
        md_array(assembly_group_t) groups = 0;
        uint8_t*  seen = md_temp_alloc_zero(temp, oper_count);
        uint32_t* rank = md_temp_alloc(temp, oper_count * sizeof(uint32_t));   // Of the first appearance of each operator
        uint32_t  rank_count = 0;
        for (size_t k = 0; k < md_array_size(opers); ++k) {
            const uint32_t o = opers[k];
            if (seen[o]) continue;
            seen[o] = 1;
            rank[o] = rank_count++;
            md_array(uint32_t) rows = 0;
            uint64_t prev = UINT64_MAX;
            for (uint32_t p = oper_row_beg[o]; p < oper_row_end[o]; ++p) {
                if (pairs[p] != prev) md_array_push(rows, (uint32_t)pairs[p], ta);
                prev = pairs[p];
            }
            size_t g = 0;
            for (; g < md_array_size(groups); ++g) {
                if (md_array_size(groups[g].rows) == md_array_size(rows) && md_memcmp(groups[g].rows, rows, md_array_bytes(rows)) == 0) break;
            }
            if (g == md_array_size(groups)) {
                const assembly_group_t grp = {.rows = rows};
                md_array_push(groups, grp, ta);
            }
            md_array_push(groups[g].oper, o, ta);
        }

        // The instances of each group: those of its rows
        for (size_t g = 0; g < md_array_size(groups); ++g) {
            assembly_group_t* grp = groups + g;
            for (size_t r = 0; r < md_array_size(grp->rows); ++r) {
                const md_assembly_gen_t* row = gen + grp->rows[r];
                for (size_t k = 0; k < row->inst_count; ++k) {
                    const md_instance_idx_t inst = row->inst[k];
                    if (inst >= 0 && (size_t)inst < inst_count) md_array_push(grp->inst, inst, ta);
                }
            }
            const size_t n = md_array_size(grp->inst);
            if (n == 0) continue;
            qsort(grp->inst, n, sizeof(md_instance_idx_t), compare_i32);
            size_t u = 1;
            for (size_t k = 1; k < n; ++k) {
                if (grp->inst[k] != grp->inst[u - 1]) grp->inst[u++] = grp->inst[k];
            }
            md_array_shrink(grp->inst, u);
        }

        // Those of the same instances are one group (rows that differ only in what else they name)
        const uint32_t first_group = (uint32_t)A->group_count;
        for (size_t g = 0; g < md_array_size(groups); ++g) {
            const assembly_group_t* grp = groups + g;
            const size_t n = md_array_size(grp->inst);
            if (n == 0) continue;
            bool merged = false;
            for (size_t h = 0; h < g && !merged; ++h) {
                const assembly_group_t* prev_grp = groups + h;
                if (md_array_size(prev_grp->inst) == n && md_memcmp(prev_grp->inst, grp->inst, n * sizeof(md_instance_idx_t)) == 0) {
                    md_array_push_array(groups[h].oper, grp->oper, md_array_size(grp->oper), ta);
                    merged = true;
                }
            }
            if (merged) {
                md_array_shrink(groups[g].inst, 0);
            }
        }
        for (size_t g = 0; g < md_array_size(groups); ++g) {
            assembly_group_t* grp = groups + g;
            if (md_array_size(grp->inst) == 0) continue;
            // The operators of a merged group in the order they first appear
            const size_t n = md_array_size(grp->oper);
            uint64_t* key = md_temp_alloc(temp, n * sizeof(uint64_t));
            for (size_t k = 0; k < n; ++k) key[k] = ((uint64_t)rank[grp->oper[k]] << 32) | grp->oper[k];
            qsort(key, n, sizeof(uint64_t), compare_u64);
            for (size_t k = 0; k < n; ++k) grp->oper[k] = (uint32_t)key[k];

            md_array_push_array(A->inst, grp->inst, md_array_size(grp->inst), alloc);
            md_array_push_array(A->oper, grp->oper, md_array_size(grp->oper), alloc);
            md_array_push(A->group_inst_offset, (uint32_t)md_array_size(A->inst), alloc);
            md_array_push(A->group_oper_offset, (uint32_t)md_array_size(A->oper), alloc);
            A->group_count += 1;
        }
        if (A->group_count == first_group) {
            MD_LOG_INFO("Assemblies: assembly '" MD_STR_FMT "' holds none of the instances of the system and was left out", MD_STR_ARG(MD_LABEL_STR(asm_id)));
            continue;
        }
        md_array_push(A->id, asm_id, alloc);
        md_array_push(A->group_offset, (uint32_t)A->group_count, alloc);
        A->count += 1;
    }

    md_temp_end(temp);
    return true;
}

const char* md_particle_kind_name(md_particle_kind_t kind) {
    switch (kind) {
    case MD_PARTICLE_ATOM:          return "atom";
    case MD_PARTICLE_BEAD:          return "bead";
    case MD_PARTICLE_VIRTUAL_SITE:  return "virtual site";
    default:                        return "";
    }
}

const char* md_residue_kind_name(md_residue_kind_t kind) {
    switch (kind) {
    case MD_RESIDUE_KIND_OTHER:         return "other";
    case MD_RESIDUE_KIND_AMINO_ACID:    return "amino acid";
    case MD_RESIDUE_KIND_NUCLEOTIDE:    return "nucleotide";
    case MD_RESIDUE_KIND_WATER:         return "water";
    case MD_RESIDUE_KIND_ION:           return "ion";
    case MD_RESIDUE_KIND_LIPID:         return "lipid";
    default:                            return "";
    }
}

void md_part_type_data_free(md_part_type_data_t* type, md_allocator_i* alloc) {
    MD_ASSERT(type);
    for (size_t i = 0; i < md_array_size(type->description); ++i) {
        if (!md_str_empty(type->description[i])) md_free_str(alloc, type->description[i]);
    }
    md_array_free(type->id, alloc);
    md_array_free(type->flags, alloc);
    md_array_free(type->description, alloc);
    md_memset(type, 0, sizeof(*type));
}

void md_part_type_data_shrink(md_part_type_data_t* type, size_t count, md_allocator_i* alloc) {
    MD_ASSERT(type);
    if (count >= type->count) return;
    for (size_t i = count; i < md_array_size(type->description); ++i) {
        if (!md_str_empty(type->description[i])) md_free_str(alloc, type->description[i]);
    }
    md_array_shrink(type->id, count);
    md_array_shrink(type->flags, count);
    md_array_shrink(type->description, count);
    type->count = count;
}

md_part_type_idx_t md_part_type_data_push(md_part_type_data_t* type, md_str_t id, md_part_type_flags_t flags, md_str_t description, md_allocator_i* alloc) {
    MD_ASSERT(type);
    md_label_t lbl = md_label_make(id);
    if (md_str_empty(id)) {
        // Its position, or the first number after it which no type has for an id (the ids of a file need not be 1..n)
        for (size_t n = type->count + 1;; ++n) {
            char buf[16];
            const int len = snprintf(buf, sizeof(buf), "%zu", n);
            lbl = md_label_make((md_str_t){buf, (size_t)MD_MAX(len, 0)});
            bool taken = false;
            for (size_t i = 0; i < md_array_size(type->id) && !taken; ++i) taken = md_label_eq(type->id[i], lbl);
            if (!taken) break;
        }
    }
    md_array_push(type->id, lbl, alloc);
    md_array_push(type->flags, flags, alloc);
    md_array_push(type->description, md_alloc_str(alloc, description), alloc);
    return (md_part_type_idx_t)(type->count++);
}

void md_molecule_type_data_free(md_molecule_type_data_t* type, md_allocator_i* alloc) {
    MD_ASSERT(type);
    for (size_t i = 0; i < md_array_size(type->name); ++i) {
        if (!md_str_empty(type->name[i])) md_free_str(alloc, type->name[i]);
    }
    md_array_free(type->name, alloc);
    md_array_free(type->flags, alloc);
    md_memset(type, 0, sizeof(*type));
}

void md_molecule_type_data_shrink(md_molecule_type_data_t* type, size_t count, md_allocator_i* alloc) {
    MD_ASSERT(type);
    if (count >= type->count) return;
    for (size_t i = count; i < md_array_size(type->name); ++i) {
        if (!md_str_empty(type->name[i])) md_free_str(alloc, type->name[i]);
    }
    md_array_shrink(type->name, count);
    md_array_shrink(type->flags, count);
    type->count = count;
}

md_molecule_type_idx_t md_molecule_type_data_push(md_molecule_type_data_t* type, md_str_t name, md_molecule_type_flags_t flags, md_allocator_i* alloc) {
    MD_ASSERT(type);
    md_array_push(type->name, md_alloc_str(alloc, name), alloc);
    md_array_push(type->flags, flags, alloc);
    return (md_molecule_type_idx_t)(type->count++);
}

const char* md_part_kind_name(md_part_kind_t kind) {
    switch (kind) {
    case MD_PART_KIND_UNKNOWN:        return "unknown";
    case MD_PART_KIND_NON_POLYMER:    return "non-polymer";
    case MD_PART_KIND_WATER:          return "water";
    case MD_PART_KIND_BRANCHED:       return "branched";
    case MD_PART_KIND_PEPTIDE:        return "peptide";
    case MD_PART_KIND_DNA:            return "DNA";
    case MD_PART_KIND_RNA:            return "RNA";
    case MD_PART_KIND_NUCLEIC:        return "nucleic acid";
    case MD_PART_KIND_POLYMER:        return "polymer";
    default:                            return "";
    }
}

const char* md_secondary_structure_name(md_secondary_structure_t ss) {
    switch (ss) {
    case MD_SECONDARY_STRUCTURE_UNKNOWN:        return "unknown";
    case MD_SECONDARY_STRUCTURE_COIL:           return "coil";
    case MD_SECONDARY_STRUCTURE_TURN:           return "turn";
    case MD_SECONDARY_STRUCTURE_BEND:           return "bend";
    case MD_SECONDARY_STRUCTURE_HELIX_310:      return "3-10 helix";
    case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:    return "alpha helix";
    case MD_SECONDARY_STRUCTURE_HELIX_PI:       return "pi helix";
    case MD_SECONDARY_STRUCTURE_BETA_SHEET:     return "beta sheet";
    case MD_SECONDARY_STRUCTURE_BETA_BRIDGE:    return "beta bridge";
    default:                                    return "";
    }
}

const char* md_hybridization_name(md_hybridization_t hyb) {
    switch (hyb) {
    case MD_HYBRIDIZATION_UNKNOWN:  return "unknown";
    case MD_HYBRIDIZATION_SP:       return "sp";
    case MD_HYBRIDIZATION_SP2:      return "sp2";
    case MD_HYBRIDIZATION_SP3:      return "sp3";
    default:                        return "";
    }
}

void md_system_free(md_system_t* sys) {
    MD_ASSERT(sys);
    MD_ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;

    if (sys->nonbonded) {
        md_nb_forcefield_free(sys->nonbonded);
        md_free(alloc, sys->nonbonded, sizeof(md_nb_forcefield_t));
        sys->nonbonded = NULL;
    }

    // ATOM
    md_array_free(sys->atom.type_idx, alloc);
    md_array_free(sys->atom.flags, alloc);
    md_array_free(sys->atom.formal_charge, alloc);
    md_array_free(sys->atom.hydrogen_count, alloc);

    // PARTICLE TYPE
    md_array_free(sys->atom.type.name, alloc);
    md_array_free(sys->atom.type.ff_type, alloc);
    md_array_free(sys->atom.type.z, alloc);
    md_array_free(sys->atom.type.mass, alloc);
    md_array_free(sys->atom.type.radius, alloc);
    md_array_free(sys->atom.type.color, alloc);
    md_array_free(sys->atom.type.flags, alloc);

    // RESIDUE
    md_array_free(sys->residue.name, alloc);
    md_array_free(sys->residue.seq_id, alloc);
    md_array_free(sys->residue.atom_offset, alloc);
    md_array_free(sys->residue.flags, alloc);
    md_array_free(sys->residue.type_idx, alloc);

    // RESIDUE TYPE
    md_array_free(sys->residue.type.name, alloc);
    md_array_free(sys->residue.type.kind, alloc);
    md_array_free(sys->residue.type.monomer, alloc);
    md_array_free(sys->residue.type.color, alloc);

    // INSTANCE
    md_array_free(sys->instance.id, alloc);
    md_array_free(sys->instance.auth_id, alloc);
    md_array_free(sys->instance.res_offset, alloc);
    md_array_free(sys->instance.part_type_idx, alloc);
    md_array_free(sys->instance.mol_type_idx, alloc);

    // CHAIN
    md_array_free(sys->chain.id, alloc);
    md_array_free(sys->chain.run_offset, alloc);
    md_array_free(sys->chain.run, alloc);
    md_array_free(sys->chain.res_chain, alloc);

    // PART
    md_array_free(sys->part.res_offset, alloc);
    md_array_free(sys->part.type_idx, alloc);
    md_part_type_data_free(&sys->part.type, alloc);

    // PROTEIN BACKBONE
    md_array_free(sys->protein_backbone.range.offset, alloc);
    md_array_free(sys->protein_backbone.range.part_idx, alloc);
    md_array_free(sys->protein_backbone.segment.atoms, alloc);
    md_array_free(sys->protein_backbone.segment.res_idx, alloc);

    // NUCLEIC BACKBONE
    md_array_free(sys->nucleic_backbone.range.offset, alloc);
    md_array_free(sys->nucleic_backbone.range.part_idx, alloc);
    md_array_free(sys->nucleic_backbone.segment.atoms, alloc);
    md_array_free(sys->nucleic_backbone.segment.res_idx, alloc);

    // LIPIDS
    md_lipid_data_free(&sys->lipid, alloc);

    // BONDS
    md_array_free(sys->bond.pairs, alloc);
    md_array_free(sys->bond.flags, alloc);
    md_array_free(sys->bond.conn.atom_idx, alloc);
    md_array_free(sys->bond.conn.bond_idx, alloc);
    md_array_free(sys->bond.conn.offset, alloc);

    md_ring_data_free(&sys->ring, sys->alloc);

    md_array_free(sys->molecule.offset, alloc);
    md_array_free(sys->molecule.atom_idx, alloc);
    md_array_free(sys->molecule.parent_idx, alloc);
    md_array_free(sys->molecule.atom_slot, alloc);
    md_array_free(sys->molecule.type_idx, alloc);
    md_molecule_type_data_free(&sys->molecule.type, alloc);

    md_fragment_data_free(&sys->fragment, alloc);

    // ASSEMBLY
    md_system_assemblies_free(sys);

    // GROUPS
    md_system_groups_free(sys);

    if (!md_str_empty(sys->description)) {
        md_free_str(alloc, sys->description);
    }

    // REFERENCE STATE
    md_system_state_free(&sys->reference);

    // ATTRIBUTES
    md_attributes_free(&sys->attributes);

    md_memset(sys, 0, sizeof(md_system_t));
}

bool md_system_state_init(md_system_state_t* state, size_t num_atoms) {
    MD_ASSERT(state);

    if (!state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    md_allocator_i* alloc = state->alloc;
    md_system_state_free(state);
    state->alloc = alloc;

    // The table allocates through its own handle, exactly as a system's does, so a producer can
    // publish into a freshly initialised state without a second setup step.
    state->attributes.alloc = alloc;

    // A freshly initialised state did not come from a run. Only md_system_extract_frame and the
    // interpolation which produces a state write a non negative frame. See md_system_state_t.
    state->frame = -1.0;

    if (num_atoms == 0) {
        return true;
    }

    const size_t capacity = MD_ALIGN_TO(num_atoms, 16);

    md_array_resize(state->xyz, capacity, alloc);

    // Zero the padding past num_atoms. The capacity is rounded up so vectorised code may load whole
    // groups of atoms; an uninitialised tail puts garbage floats into those lanes.
    const size_t tail_bytes = (capacity - num_atoms) * sizeof(md_f32x3_t);
    if (tail_bytes > 0) {
        md_memset(state->xyz + num_atoms, 0, tail_bytes);
    }

    state->num_atoms = num_atoms;

    return true;
}

void md_system_state_free(md_system_state_t* state) {
    MD_ASSERT(state);

    // A view owns nothing; zeroing it is the whole job.
    if (state->alloc) {
        md_array_free(state->xyz, state->alloc);
        md_attributes_free(&state->attributes);
    }
    md_allocator_i* alloc = state->alloc;
    md_memset(state, 0, sizeof(md_system_state_t));
    state->alloc = alloc;
}

bool md_system_state_copy(md_system_state_t* dst, const md_system_state_t* src) {
    MD_ASSERT(dst);
    MD_ASSERT(src);

    if (dst == src) {
        return true;
    }
    if (!md_system_state_init(dst, src->num_atoms)) {
        return false;
    }
    if (src->num_atoms > 0 && src->xyz) {
        md_memcpy(dst->xyz, src->xyz, src->num_atoms * sizeof(md_f32x3_t));
    }
    dst->unitcell = src->unitcell;
    dst->frame    = src->frame;

    // The attribute table is deliberately NOT carried over. The one caller is md_system_infer
    // writing sys->reference, and topology inference reads coordinates and the cell - nothing else.
    // Duplicating a frame's velocities into a reference snapshot would put a second copy of them in
    // the system with nothing keeping the two in step, which is the whole class of problem the
    // reference field's own rule exists to prevent.
    return true;
}

void md_lipid_data_free(md_lipid_data_t* lipid, md_allocator_i* alloc) {
    MD_ASSERT(lipid);
    md_array_free(lipid->residue_offset, alloc);
    md_array_free(lipid->residue_idx, alloc);
    md_array_free(lipid->head_offset, alloc);
    md_array_free(lipid->tail_offset, alloc);
    md_array_free(lipid->head.offset, alloc);
    md_array_free(lipid->head.atom_idx, alloc);
    md_array_free(lipid->tail.spine_offset, alloc);
    md_array_free(lipid->tail.spine_idx, alloc);
    md_array_free(lipid->tail.atom_offset, alloc);
    md_array_free(lipid->tail.atom_idx, alloc);
    md_array_free(lipid->tail.anchor_beg, alloc);
    md_array_free(lipid->tail.anchor_end, alloc);
    md_array_free(lipid->tail.head_beg, alloc);
    md_array_free(lipid->tail.head_end, alloc);
    md_array_free(lipid->tail.kind, alloc);
    md_memset(lipid, 0, sizeof(md_lipid_data_t));
}

// ### COPY ###

#define COPY_ARRAY(dst, src, alloc) do { (dst) = NULL; if (md_array_size(src)) md_array_push_array((dst), (src), md_array_size(src), (alloc)); } while (0)

static void copy_strings(md_str_t** dst, const md_str_t* src, md_allocator_i* alloc) {
    *dst = NULL;
    const size_t n = md_array_size(src);
    if (!n) return;
    md_array_resize(*dst, n, alloc);
    for (size_t i = 0; i < n; ++i) (*dst)[i] = md_str_empty(src[i]) ? (md_str_t){0} : md_alloc_str(alloc, src[i]);
}

bool md_system_copy(md_system_t* dst, const md_system_t* src) {
    MD_ASSERT(dst);
    MD_ASSERT(src);
    if (dst == src) return true;
    if (!dst->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    md_allocator_i* alloc = dst->alloc;
    md_system_free(dst);

    // The counts and everything which is not an array, as they are; the arrays below anew
    *dst = *src;
    dst->alloc = alloc;
    dst->nonbonded = NULL;
    md_memset(&dst->attributes, 0, sizeof(dst->attributes));
    dst->attributes.alloc = alloc;
    dst->description = md_str_empty(src->description) ? (md_str_t){0} : md_alloc_str(alloc, src->description);

    // Every array, CORE and DERIVED (md_system_fields.h)
#define COPY_FIELD(f) COPY_ARRAY(dst->f, src->f, alloc);
#define COPY_STRINGS(f) copy_strings(&dst->f, src->f, alloc);
    MD_SYSTEM_CORE_ARRAYS(COPY_FIELD)
    MD_SYSTEM_DERIVED_ARRAYS(COPY_FIELD)
    MD_SYSTEM_DERIVED_STRINGS(COPY_STRINGS)
#undef COPY_FIELD
#undef COPY_STRINGS

    // GROUPS, each set with words of its own
    groups_copy(&dst->group, &src->group, alloc);

    // REFERENCE STATE
    md_memset(&dst->reference, 0, sizeof(dst->reference));
    if (src->reference.xyz || src->reference.num_atoms) {
        dst->reference.alloc = alloc;
        md_system_state_copy(&dst->reference, &src->reference);
    }
    return true;
}

#undef COPY_ARRAY

void md_system_reset(md_system_t* sys) {
    MD_ASSERT(sys);
    md_allocator_i* alloc = sys->alloc;
    md_system_free(sys);
    sys->alloc = alloc;
    // The table allocates through its own handle, so every loader gets a usable table without
    // having to remember to wire this up itself.
    sys->attributes.alloc = alloc;
    // Every loader starts here, so whatever it builds has a version of its own
    md_system_topology_changed(sys);
}

void md_system_bond_build_connectivity(md_system_t* sys) {
    MD_ASSERT(sys);
	md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
}



// ### EXTRACTION ###

typedef enum extract_kind_t {
    EXTRACT_POSITION,   // atom/position into x, y, z
    EXTRACT_CELL,       // unitcell into unitcell
    EXTRACT_OTHER,      // anything else into out->attributes
} extract_kind_t;

typedef struct extract_entry_t {
    md_str_t             path;      // relative to the run, owned
    extract_kind_t    kind;
    md_attribute_id_t id;
    md_attribute_id_t axis_id;
} extract_entry_t;

struct md_system_extract_t {
    const md_system_t*     sys;
    md_allocator_i*        arena;       // everything the context owns
    md_attribute_id_t      run_axis_id;
    size_t                 num_frames;

    extract_entry_t*       entries;
    size_t                 num_entries;

    md_attribute_io_t      io;
    md_thread_id_t         owner;       // the thread that extracts; 0 until the first frame
};

static md_str_t extract_run_path(char* buf, size_t cap, md_str_t run, md_str_t leaf) {
    const int len = md_snprintf(buf, cap, MD_STR_FMT "/" MD_STR_FMT, MD_STR_ARG(run), MD_STR_ARG(leaf));
    return (len > 0 && (size_t)len < cap) ? (md_str_t){buf, (size_t)len} : (md_str_t){0};
}

md_system_extract_t* md_system_extract_begin(const md_system_t* sys, md_str_t run, const md_str_t paths[], size_t num_paths, md_allocator_i* alloc) {
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(paths || num_paths == 0);

    md_allocator_i* arena = md_arena_allocator_create(alloc, MD_KILOBYTES(16));
    md_system_extract_t* ex = md_alloc(arena, sizeof(md_system_extract_t));
    md_memset(ex, 0, sizeof(md_system_extract_t));
    ex->sys      = sys;
    ex->arena    = arena;
    ex->io.alloc = arena;
    ex->entries  = md_alloc(arena, MD_MAX(num_paths, 1) * sizeof(extract_entry_t));

    const md_attributes_t* attributes = &sys->attributes;

    {
        const md_attribute_t* axis = md_str_empty(run) ? NULL : md_attributes_find_in(attributes, run, MD_STR_LIT("time"));
        if (!axis || md_attributes_axis(attributes, axis) != axis) {
            MD_LOG_ERROR("No run '" MD_STR_FMT "' to extract from", MD_STR_ARG(run));
            goto fail;
        }
        ex->run_axis_id = axis->id;
        ex->num_frames  = axis->format.shape[0];
    }

    for (size_t i = 0; i < num_paths; ++i) {
        extract_entry_t* e = &ex->entries[ex->num_entries];
        md_memset(e, 0, sizeof(*e));
        e->path = md_alloc_str(arena, paths[i]);
        e->kind = md_str_eq(paths[i], MD_STR_LIT("atom/position")) ? EXTRACT_POSITION :
                  md_str_eq(paths[i], MD_STR_LIT("unitcell"))      ? EXTRACT_CELL : EXTRACT_OTHER;

        const md_attribute_t* attr = md_attributes_find_in(attributes, run, paths[i]);
        if (!attr) {
            MD_LOG_ERROR("Nothing at '" MD_STR_FMT "' in '" MD_STR_FMT "' to extract", MD_STR_ARG(paths[i]), MD_STR_ARG(run));
            goto fail;
        }

        if (!(attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL)) {
            MD_LOG_ERROR("'" MD_STR_FMT "' does not vary over the run; read it from the system directly", MD_STR_ARG(attr->path));
            goto fail;
        }
        if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {
            MD_LOG_ERROR("'" MD_STR_FMT "' is text and cannot be carried by a state", MD_STR_ARG(attr->path));
            goto fail;
        }
        if (e->kind == EXTRACT_POSITION && !(attr->format.rank == 2 && attr->format.components == 3)) {
            MD_LOG_ERROR("'" MD_STR_FMT "' is not {F,N} with three components", MD_STR_ARG(attr->path));
            goto fail;
        }
        if (e->kind == EXTRACT_CELL && !(attr->format.rank == 3 && attr->format.shape[1] == 3 && attr->format.shape[2] == 3 && attr->format.components == 1)) {
            MD_LOG_ERROR("'" MD_STR_FMT "' is not {F,3,3}", MD_STR_ARG(attr->path));
            goto fail;
        }
        const md_attribute_t* axis = md_attributes_axis(attributes, attr);
        if (!axis) {
            MD_LOG_ERROR("'" MD_STR_FMT "' has no frame axis", MD_STR_ARG(attr->path));
            goto fail;
        }
        e->id      = attr->id;
        e->axis_id = axis->id;
        ex->num_entries += 1;
    }

    return ex;

fail:
    md_arena_allocator_destroy(arena);
    return NULL;
}

void md_system_extract_end(md_system_extract_t* ex) {
    if (!ex) return;
    md_attribute_io_close_all(&ex->io);
    md_arena_allocator_destroy(ex->arena);
}

bool md_system_extract_frame(md_system_extract_t* ex, int64_t frame, md_system_state_t* out) {
    MD_ASSERT(ex);
    MD_ASSERT(out);

    // One thread at a time. Recorded on first use rather than at begin, so a context may be made on
    // one thread and handed to the one that uses it.
    const md_thread_id_t tid = md_thread_id();
    if (ex->owner == 0) {
        ex->owner = tid;
    }
    MD_ASSERT(ex->owner == tid && "an extraction context is used by one thread at a time");

    if (frame < 0 || (size_t)frame >= ex->num_frames) {
        MD_LOG_ERROR("Frame %" PRId64 " is outside the %zu frames being extracted from", frame, ex->num_frames);
        return false;
    }

    const md_attributes_t* attributes = &ex->sys->attributes;
    const bool want_coords = out->xyz != NULL;

    const md_attribute_t* run_axis = md_attributes_get(attributes, ex->run_axis_id);
    if (!run_axis) {
        MD_LOG_ERROR("The run being extracted from is gone; end the context before removing it");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    bool result = true;

    for (size_t i = 0; i < ex->num_entries && result; ++i) {
        const extract_entry_t* e = &ex->entries[i];
        const md_attribute_t* attr = md_attributes_get(attributes, e->id);
        const md_attribute_t* axis = md_attributes_get(attributes, e->axis_id);
        if (!attr || !axis) {
            MD_LOG_ERROR("'" MD_STR_FMT "' is gone from the run being extracted from", MD_STR_ARG(e->path));
            result = false;
            break;
        }

        size_t row = (size_t)frame;
        if (e->axis_id != ex->run_axis_id && !md_attribute_axis_map(&row, run_axis, (size_t)frame, axis)) {
            if (e->kind != EXTRACT_OTHER) {
                MD_LOG_ERROR("'" MD_STR_FMT "' has no value at frame %" PRId64, MD_STR_ARG(e->path), frame);
                result = false;
                break;
            }
            // Sampled at other times than the frames (velocities written every fifth frame): this
            // frame has none, and the state says so by not carrying it - never the value of another
            // frame, which a state reused from the last one would otherwise still hold.
            if (out->attributes.alloc) {
                const md_attribute_t* prev = md_attributes_find(&out->attributes, e->path);
                if (prev) md_attributes_remove(&out->attributes, prev->id);
            }
            continue;
        }
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)row);

        switch (e->kind) {
        case EXTRACT_POSITION: {
            if (!want_coords) break;
            const size_t N = attr->format.shape[1];
            if (out->num_atoms != 0 && out->num_atoms != N) {
                MD_LOG_ERROR("The state holds %zu atoms, '" MD_STR_FMT "' %zu", out->num_atoms, MD_STR_ARG(attr->path), N);
                result = false;
                break;
            }
            if (md_attribute_extract_io_f32((float*)out->xyz, N * 3, attr, slice, md_unit_angstrom(), &ex->io) != N * 3) {
                result = false;
                break;
            }
            out->num_atoms = N;
            break;
        }
        case EXTRACT_CELL: {
            float box[3][3];
            if (md_attribute_extract_io_f32(&box[0][0], 9, attr, slice, md_unit_angstrom(), &ex->io) != 9) {
                result = false;
                break;
            }
            const bool empty = box[0][0] == 0.0f && box[1][1] == 0.0f && box[2][2] == 0.0f;
            out->unitcell = empty ? (md_unitcell_t){0} : md_unitcell_from_matrix_float(box);
            break;
        }
        case EXTRACT_OTHER: {
            if (!out->attributes.alloc) {
                MD_LOG_ERROR("'" MD_STR_FMT "' goes into the state's attributes, and the state has none (no allocator)", MD_STR_ARG(e->path));
                result = false;
                break;
            }
            md_attribute_format_t fmt;
            md_attribute_slice_format(&fmt, attr, slice);
            const size_t bytes = md_attribute_byte_size(&fmt);
            void* data = md_temp_alloc(temp, MD_MAX(bytes, 1));
            if (!data || !md_attribute_read_stored(data, attr, slice, &ex->io)) {
                result = false;
                break;
            }
            const md_attribute_desc_t desc = {
                .path   = e->path,
                .format = fmt,
                .unit   = attr->unit,
                .label  = attr->label,
                .data   = data,
                .byte_size = bytes,
            };
            if (!md_attributes_replace(&out->attributes, &desc)) {
                result = false;
            }
            break;
        }
        }
    }

    md_temp_end(temp);
    if (result) {
        out->frame = (double)frame;
    }
    return result;
}

// ### RUNS ###

static md_str_t run_cache_path(char* buf, size_t cap, md_str_t source_path) {
    const int len = md_snprintf(buf, cap, MD_STR_FMT ".cache", MD_STR_ARG(source_path));
    return (len > 0 && (size_t)len < cap) ? (md_str_t){buf, (size_t)len} : (md_str_t){0};
}

bool md_run_cache_open(md_file_t* out_file, md_run_cache_header_t* out_header, md_str_t source_path, uint64_t magic, uint64_t version) {
    MD_ASSERT(out_file);
    MD_ASSERT(out_header);
    md_memset(out_file, 0, sizeof(*out_file));

    md_file_info_t info = {0};
    if (!md_file_info_extract_from_path(source_path, &info)) {
        return false;
    }
    char buf[4096];
    const md_str_t cache_path = run_cache_path(buf, sizeof(buf), source_path);
    md_file_t file = {0};
    if (md_str_empty(cache_path) || !md_file_open(&file, cache_path, MD_FILE_READ)) {
        return false;
    }

    const char* stale = NULL;
    if (md_file_read(file, out_header, sizeof(*out_header)) != sizeof(*out_header)) {
        stale = "incomplete header";
    } else if (out_header->magic != magic) {
        stale = "not this format's";
    } else if (out_header->version != version) {
        stale = "an older version";
    } else if (out_header->source_size != (uint64_t)info.size) {
        stale = "made from a file of another size";
    } else if (out_header->source_modified != info.modified_time) {
        stale = "made from a file modified at another time";
    } else if (out_header->num_frames == 0 || out_header->num_atoms == 0) {
        stale = "empty";
    }
    if (stale) {
        MD_LOG_INFO("'" MD_STR_FMT "' is %s; it will be made again", MD_STR_ARG(cache_path), stale);
        md_file_close(&file);
        return false;
    }
    *out_file = file;
    return true;
}

bool md_run_cache_create(md_file_t* out_file, md_str_t source_path, const md_file_info_t* scanned, uint64_t magic, uint64_t version, size_t num_atoms, size_t num_frames) {
    MD_ASSERT(out_file);
    MD_ASSERT(scanned);
    md_memset(out_file, 0, sizeof(*out_file));
    const md_file_info_t info = *scanned;
    char buf[4096];
    const md_str_t cache_path = run_cache_path(buf, sizeof(buf), source_path);
    md_file_t file = {0};
    if (md_str_empty(cache_path) || !md_file_open(&file, cache_path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        MD_LOG_INFO("Could not create '" MD_STR_FMT "'", MD_STR_ARG(cache_path));
        return false;
    }
    const md_run_cache_header_t header = {
        .magic           = magic,
        .version         = version,
        .source_size     = (uint64_t)info.size,
        .source_modified = info.modified_time,
        .num_atoms       = num_atoms,
        .num_frames      = num_frames,
    };
    if (md_file_write(file, &header, sizeof(header)) != sizeof(header)) {
        md_file_close(&file);
        return false;
    }
    *out_file = file;
    return true;
}

md_str_t md_run_path(char* buf, size_t cap, md_str_t run, md_str_t leaf) {
    return extract_run_path(buf, cap, run, leaf);
}

bool md_run_publish(md_system_t* sys, md_str_t run, const md_run_desc_t* d) {
    MD_ASSERT(sys);
    MD_ASSERT(d);
    if (md_str_empty(run) || !d->time || !d->source_offset || !d->source_size || !d->position_virt || md_str_empty(d->source_path)) {
        MD_LOG_ERROR("An incomplete description of the run '" MD_STR_FMT "'", MD_STR_ARG(run));
        return false;
    }
    const size_t F = d->num_frames;
    const size_t N = d->num_atoms;
    if (F == 0 || N == 0 || F > UINT32_MAX || N > UINT32_MAX) {
        MD_LOG_ERROR("The run '" MD_STR_FMT "' has %zu frames of %zu atoms, which cannot be published", MD_STR_ARG(run), F, N);
        return false;
    }
    if (sys->atom.count != 0 && sys->atom.count != N) {
        MD_LOG_ERROR("The run '" MD_STR_FMT "' holds %zu atoms, the system %zu", MD_STR_ARG(run), N, sys->atom.count);
        return false;
    }

    md_attributes_t* attributes = &sys->attributes;
    if (!attributes->alloc) {
        attributes->alloc = sys->alloc;
    }

    const md_attribute_format_t series_f64 = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)F } };
    const md_attribute_format_t series_i64 = { .type = MD_ATTRIBUTE_TYPE_I64, .components = 1, .rank = 1, .shape = { (uint32_t)F } };
    const md_attribute_format_t cell_format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 3, .shape = { (uint32_t)F, 3, 3 } };
    char buf[512];

    // The axis first: everything temporal below is checked against it.
    bool ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("time")), .format = series_f64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = d->time_unit, .label = MD_STR_INIT("Time"),
        .data = d->time, .byte_size = F * sizeof(double)});

    if (ok && d->step) {
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("step")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_none(), .label = MD_STR_INIT("Step"),
            .description = MD_STR_INIT("The file's own notion of where a frame sits in the run, not the frame ordinal"),
            .data = d->step, .byte_size = F * sizeof(int64_t)});
    }

    if (ok && (d->unitcell || d->unitcell_virt)) {
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("unitcell")), .format = cell_format,
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_angstrom(), .label = MD_STR_INIT("Unit Cell"),
            .data = d->unitcell_virt ? NULL : d->unitcell, .byte_size = d->unitcell_virt ? 0 : F * 9 * sizeof(float),
            .virt = d->unitcell_virt});
    }

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/path")),
        .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
        .unit = md_unit_none(), .data = &d->source_path, .byte_size = sizeof(md_str_t)});

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/offset")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = d->source_offset, .byte_size = F * sizeof(int64_t)});

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/size")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = d->source_size, .byte_size = F * sizeof(int64_t)});

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("atom/position")),
        .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)F, (uint32_t)N } },
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_angstrom(), .label = MD_STR_INIT("Position"),
        .virt = d->position_virt});

    if (!ok) {
        MD_LOG_ERROR("Failed to publish the run '" MD_STR_FMT "' from '" MD_STR_FMT "'", MD_STR_ARG(run), MD_STR_ARG(d->source_path));
        md_attributes_remove_prefix(attributes, run);
    }
    return ok;
}

bool md_run_source(md_run_source_t* out, const md_attributes_t* attributes, const md_attribute_t* attr, md_str_t leaf) {
    MD_ASSERT(out);
    MD_ASSERT(attributes);
    MD_ASSERT(attr);
    md_memset(out, 0, sizeof(*out));

    const md_str_t p = attr->path;
    if (p.len <= leaf.len + 1 || !md_str_ends_with(p, leaf) || p.ptr[p.len - leaf.len - 1] != '/') {
        return false;
    }
    out->run = md_str_substr(p, 0, p.len - leaf.len - 1);

    const md_attribute_t* path   = md_attributes_find_in(attributes, out->run, MD_STR_LIT("source/path"));
    const md_attribute_t* offset = md_attributes_find_in(attributes, out->run, MD_STR_LIT("source/offset"));
    const md_attribute_t* size   = md_attributes_find_in(attributes, out->run, MD_STR_LIT("source/size"));
    const int64_t* offsets = (const int64_t*)md_attribute_view(offset, MD_ATTRIBUTE_TYPE_I64, 1, 1);
    const int64_t* sizes   = (const int64_t*)md_attribute_view(size,   MD_ATTRIBUTE_TYPE_I64, 1, 1);
    if (!path || path->format.type != MD_ATTRIBUTE_TYPE_STR || !offsets || !sizes ||
        offset->format.shape[0] != size->format.shape[0]) {
        MD_LOG_ERROR("The run '" MD_STR_FMT "' has lost its source attributes", MD_STR_ARG(out->run));
        return false;
    }
    out->path       = md_attribute_str(attributes, path, 0);
    out->offset     = offsets;
    out->size       = sizes;
    out->num_frames = offset->format.shape[0];
    return true;
}

// Lower case letters and digits, anything else one '_' between them
static size_t run_series_slug(char* buf, size_t cap, md_str_t name) {
    size_t len = 0;
    bool pending_sep = false;
    for (size_t i = 0; i < name.len && len + 2 < cap; ++i) {
        char c = name.ptr[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        const bool keep = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
        if (!keep) {
            pending_sep = (len > 0);
            continue;
        }
        if (pending_sep) {
            buf[len++] = '_';
            pending_sep = false;
        }
        buf[len++] = c;
    }
    buf[len] = '\0';
    return len;
}

bool md_run_publish_series(md_system_t* sys, md_str_t run, const md_run_series_desc_t* d) {
    MD_ASSERT(sys);
    MD_ASSERT(d);
    md_attributes_t* attributes = &sys->attributes;
    if (!attributes->alloc) {
        attributes->alloc = sys->alloc;
    }

    char buf[512];
    const md_attribute_t* run_axis = md_str_empty(run) ? NULL : md_attributes_find_in(attributes, run, MD_STR_LIT("time"));
    if (!run_axis || md_attributes_axis(attributes, run_axis) != run_axis) {
        MD_LOG_ERROR("No run '" MD_STR_FMT "' to publish '" MD_STR_FMT "' along", MD_STR_ARG(run), MD_STR_ARG(d->group));
        return false;
    }
    const size_t F = run_axis->format.shape[0];
    const size_t R = d->num_rows;
    if (R == 0 || d->num_columns == 0 || R > UINT32_MAX || !d->names || !d->columns) {
        MD_LOG_ERROR("'" MD_STR_FMT "' holds no values", MD_STR_ARG(d->group));
        return false;
    }

    char group_buf[512];
    const md_str_t group = md_run_path(group_buf, sizeof(group_buf), run, d->group);
    if (md_str_empty(group) || md_str_empty(d->group)) {
        MD_LOG_ERROR("No group to publish the series under");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    md_unit_t time_unit = d->time_unit;
    if (d->time) {
        if (md_unit_is_none(time_unit)) {
            time_unit = run_axis->unit;
        }
        for (size_t i = 1; i < R; ++i) {
            if (d->time[i] < d->time[i - 1]) {
                MD_LOG_ERROR("'" MD_STR_FMT "': time decreases at row %zu", MD_STR_ARG(d->group), i);
                goto done;
            }
        }
        // Whether the rows cover the run is a question about the values, answered against a scratch
        // table before the real one is touched.
        md_attributes_t probe = { .alloc = temp_alloc };
        const md_attribute_id_t probe_id = md_attributes_create(&probe, &(md_attribute_desc_t){
            .path = MD_STR_INIT("time"), .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = time_unit, .data = d->time, .byte_size = R * sizeof(double)});
        const md_attribute_t* axis = probe_id ? md_attributes_get(&probe, probe_id) : NULL;
        if (!axis) goto done;
        for (size_t f = 0; f < F; ++f) {
            size_t row;
            if (!md_attribute_axis_map(&row, run_axis, f, axis)) {
                MD_LOG_ERROR("'" MD_STR_FMT "' has no row at the time of frame %zu of '" MD_STR_FMT "'; it belongs to another run", MD_STR_ARG(d->group), f, MD_STR_ARG(run));
                goto done;
            }
        }
    } else if (R != F) {
        MD_LOG_ERROR("'" MD_STR_FMT "' has %zu rows and no time, and '" MD_STR_FMT "' %zu frames: the rows cannot be the frames", MD_STR_ARG(d->group), R, MD_STR_ARG(run), F);
        goto done;
    }

    md_attributes_remove_prefix(attributes, group);

    const md_attribute_format_t series = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = { (uint32_t)R } };
    bool ok = true;
    if (d->time) {
        ok = md_attributes_create(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, MD_STR_LIT("time")),
            .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = time_unit, .label = MD_STR_INIT("Time"),
            .data = d->time, .byte_size = R * sizeof(double)}) != MD_ATTRIBUTE_INVALID;
    }

    md_array(md_str_t) taken = 0;
    md_array_push(taken, MD_STR_LIT("time"),   temp_alloc);
    md_array_push(taken, MD_STR_LIT("source"), temp_alloc);
    for (size_t k = 0; ok && k < d->num_columns; ++k) {
        char slug[128];
        size_t len = run_series_slug(slug, sizeof(slug) - 8, d->names[k]);
        if (len == 0) {
            len = (size_t)md_snprintf(slug, sizeof(slug), "column%zu", k + 1);
        }
        const size_t base = len;
        for (int n = 2;; ++n) {
            bool clash = false;
            for (size_t t = 0; t < md_array_size(taken); ++t) {
                if (md_str_eq(taken[t], (md_str_t){slug, len})) { clash = true; break; }
            }
            if (!clash) break;
            len = base + (size_t)md_snprintf(slug + base, sizeof(slug) - base, "_%d", n);
        }
        md_array_push(taken, md_alloc_str(temp_alloc, (md_str_t){slug, len}), temp_alloc);

        ok = md_attributes_create(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, (md_str_t){slug, len}),
            .format = series, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = d->units ? d->units[k] : md_unit_none(), .label = d->names[k],
            .data = d->columns[k], .byte_size = R * sizeof(float)}) != MD_ATTRIBUTE_INVALID;
    }
    if (ok && !md_str_empty(d->source_path)) {
        ok = md_attributes_create(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, MD_STR_LIT("source")),
            .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
            .unit = md_unit_none(), .data = &d->source_path, .byte_size = sizeof(md_str_t)}) != MD_ATTRIBUTE_INVALID;
    }
    if (!ok) {
        MD_LOG_ERROR("Failed to publish '" MD_STR_FMT "'", MD_STR_ARG(group));
        md_attributes_remove_prefix(attributes, group);
        goto done;
    }
    result = true;

done:
    md_temp_end(temp);
    return result;
}

#ifdef __cplusplus
}
#endif
