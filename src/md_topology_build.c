#include <md_topology_build.h>

#include <md_system.h>
#include <md_util.h>
#include <md_element.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <math.h>
#include <string.h>

#define TOPO_LOG_INFO(prefix, fmt, ...)  MD_LOG_INFO("%s: " fmt, prefix, ##__VA_ARGS__)
#define TOPO_LOG_ERROR(prefix, fmt, ...) MD_LOG_ERROR("%s: " fmt, prefix, ##__VA_ARGS__)

// The element of an atom from its name and mass, for a topology without ATOMIC_NUMBER. The name first when its
// element is of about the mass (hydrogen mass repartitioning puts some 3 Da on a hydrogen), else the mass alone.
static md_atomic_number_t topology_infer_element(md_str_t name, md_str_t res_name, size_t res_size, double mass) {
    const md_atomic_number_t z = md_atomic_number_infer_from_label(name, res_name, res_size);
    if (z) {
        const double m = md_atomic_number_mass(z);
        // Within what hydrogen mass repartitioning moves (some 2 Da from a heavy atom, onto a hydrogen of up to 4):
        // a name that reads as another element than the mass (SIA, silicon not sulfur) is not trusted
        if (mass <= 0.0 || fabs(mass - m) <= 2.5 || (z == 1 && mass < 4.5)) return z;
    }
    if (!(mass > 0.0)) return 0;
    const md_atomic_number_t exact = md_atomic_number_infer_from_mass((float)mass);
    if (exact) return exact;
    // Force fields round masses their own way (CHARMM's potassium is 39.102): the nearest element within a percent
    md_atomic_number_t best = 0;
    double best_diff = 0.01 * mass;
    for (md_atomic_number_t k = 1; k <= 118; ++k) {
        const double diff = fabs(md_atomic_number_mass(k) - mass);
        if (diff < best_diff) {
            best_diff = diff;
            best = k;
        }
    }
    return best;
}

// Particles without an element and without a mass are virtual sites (extra points); bonds to them are dropped, as
// for a GROMACS topology: a site is tied to its atom by its construction, and kept with its molecule by inference.
static bool topology_is_virtual(const md_topology_desc_t* d, size_t i, md_atomic_number_t z) {
    return z == 0 && (!d->mass || d->mass[i] == 0.0);
}

bool md_topology_system_init(md_system_t* sys, md_system_state_t* state, const md_topology_desc_t* desc) {
    MD_ASSERT(sys);
    MD_ASSERT(state);
    MD_ASSERT(desc);
    if (!sys->alloc || !state->alloc) {
        TOPO_LOG_ERROR(desc->format, "system or state allocator not set");
        return false;
    }
    if (!desc->num_atoms || !desc->num_residues || !desc->atom_name || !desc->res_name || !desc->res_offset) {
        TOPO_LOG_ERROR(desc->format, "no atoms or residues");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;
    const md_topology_desc_t d = *desc;

    md_system_reset(sys);
    md_system_state_init(state, 0);
    md_allocator_i* alloc = sys->alloc;
    const size_t N = d.num_atoms, R = d.num_residues;

    // ## Elements
    md_atomic_number_t* z = md_temp_alloc_array(temp, md_atomic_number_t, N);
    size_t num_inferred = 0;
    for (size_t r = 0; r < R; ++r) {
        const size_t beg = d.res_offset[r], end = d.res_offset[r + 1];
        for (size_t i = beg; i < end; ++i) {
            const int64_t given = d.atomic_number ? d.atomic_number[i] : -1;
            if (given > 0 && given <= 118) {
                z[i] = (md_atomic_number_t)given;
            } else if (given == 0 && d.mass && d.mass[i] == 0.0) {
                z[i] = 0;   // an extra point
            } else {
                z[i] = topology_infer_element(d.atom_name[i], d.res_name[r], end - beg, d.mass ? d.mass[i] : 0.0);
                num_inferred += z[i] != 0;
            }
        }
    }
    if (num_inferred) {
        TOPO_LOG_INFO(d.format, "the elements of %zu atoms were inferred from their names and masses", num_inferred);
    }

    // ## Particle types: one per name, force field type, element, mass and kind
    const size_t capacity = MD_ROUND_UP(N, 16);
    md_array_resize(sys->atom.type_idx, capacity, alloc);
    md_array_resize(sys->atom.flags, capacity, alloc);
    md_memset(sys->atom.type_idx, 0, capacity * sizeof(md_particle_type_idx_t));
    md_memset(sys->atom.flags, 0, capacity * sizeof(md_atom_flags_t));
    sys->atom.type.count = 0;
    md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, alloc);

    md_hashmap32_t type_map = { .allocator = temp_alloc };
    bool names_cut = false;
    for (size_t i = 0; i < N; ++i) {
        const float mass = d.mass ? (float)d.mass[i] : 0.0f;
        const md_str_t ff = d.atom_type ? d.atom_type[i] : (md_str_t){0};
        const md_particle_kind_t kind = topology_is_virtual(&d, i, z[i]) ? MD_PARTICLE_VIRTUAL_SITE : MD_PARTICLE_ATOM;
        const struct { float mass; uint32_t z; uint32_t kind; } key_data = { mass, z[i], (uint32_t)kind };
        const uint64_t key = md_hash64_str(d.atom_name[i], md_hash64_str(ff, md_hash64(&key_data, sizeof(key_data), 0)));
        const uint32_t* cached = md_hashmap_get(&type_map, key);
        if (cached) {
            sys->atom.type_idx[i] = (md_particle_type_idx_t)*cached;
            continue;
        }
        const float radius = kind == MD_PARTICLE_VIRTUAL_SITE ? 0.0f : md_atomic_number_vdw_radius(z[i]);
        const uint32_t color = md_atomic_number_cpk_color(z[i]);
        const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, kind);
        const md_particle_type_idx_t type = md_particle_type_add(&sys->atom.type, d.atom_name[i], ff, z[i], mass, radius, color, flags, alloc);
        names_cut |= d.atom_name[i].len > sizeof(md_label_t) || ff.len > sizeof(md_label16_t);
        md_hashmap_add(&type_map, key, (uint32_t)type);
        sys->atom.type_idx[i] = type;
    }
    sys->atom.count = N;
    if (names_cut) {
        TOPO_LOG_INFO(d.format, "names longer than a label are shown cut short; the types stay apart");
    }

    // ## Residues
    for (size_t r = 0; r < R; ++r) {
        md_array_push(sys->residue.atom_offset, d.res_offset[r], alloc);
        md_array_push(sys->residue.name, md_label_make(d.res_name[r]), alloc);
        md_array_push(sys->residue.seq_id, (md_sequence_id_t)(d.res_number ? d.res_number[r] : (int64_t)r + 1), alloc);
        md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, alloc);
    }
    md_array_push(sys->residue.atom_offset, (uint32_t)N, alloc);
    sys->residue.count = R;

    // ## Bonds
    {
        const md_bond_flags_t flags = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
        md_array_resize(sys->bond.pairs, d.num_bonds, alloc);
        md_array_resize(sys->bond.flags, d.num_bonds, alloc);
        size_t kept = 0;
        for (size_t b = 0; b < d.num_bonds; ++b) {
            const md_atom_pair_t pair = d.bonds[b];
            if (topology_is_virtual(&d, (size_t)pair.idx[0], z[pair.idx[0]]) || topology_is_virtual(&d, (size_t)pair.idx[1], z[pair.idx[1]])) continue;
            sys->bond.pairs[kept] = pair;
            sys->bond.flags[kept] = flags;
            kept += 1;
        }
        md_array_shrink(sys->bond.pairs, kept);
        md_array_shrink(sys->bond.flags, kept);
        sys->bond.count = kept;
        md_bond_build_connectivity(&sys->bond, N, alloc);
    }

    md_util_system_augment_particle_types(sys);
    md_util_system_infer_residue_flags(sys, state);

    // ## Chains, as ParmEd keeps them from the structure the topology was made from
    if (d.res_chain) {
        bool any = false;
        for (size_t r = 0; r < R && !any; ++r) any = !md_str_empty(d.res_chain[r]);
        if (any) md_util_system_instances_from_chain_ids(sys, d.res_chain, NULL);
    }

    // ## Attributes
    if (!sys->attributes.alloc) sys->attributes.alloc = alloc;
    float* column = md_temp_alloc_array(temp, float, N);
    if (d.charge) {
        for (size_t i = 0; i < N; ++i) column[i] = (float)d.charge[i];
        md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, column, N);
    }
    if (d.mass) {
        for (size_t i = 0; i < N; ++i) column[i] = (float)d.mass[i];
        md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/mass"), md_unit_dalton(), 1, column, N);
    }
    if (!md_str_empty(d.title)) {
        sys->description = md_alloc_str(alloc, d.title);
    }

    state->unitcell = d.cell;
    result = true;
    md_temp_end(temp);
    return result;
}

// Names as a structure file has them against the topology's: equal ignoring case, or the structure's a cut short
// form of the topology's.
static bool topology_name_matches(md_str_t structure, md_str_t topology) {
    if (md_str_eq_ignore_case(structure, topology)) return true;
    return structure.len >= 3 && topology.len > structure.len && md_str_begins_with_ignore_case(topology, structure);
}

bool md_topology_system_supplement(md_system_t* sys, const md_topology_desc_t* desc, md_str_t filename) {
    MD_ASSERT(sys);
    MD_ASSERT(desc);
    if (!sys->alloc || sys->atom.count == 0) {
        TOPO_LOG_ERROR(desc->format, "a topology can only supplement a loaded system");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    bool result = false;
    const md_topology_desc_t d = *desc;

    const size_t N = sys->atom.count;
    if (d.num_atoms != N) {
        TOPO_LOG_ERROR(d.format, "'" MD_STR_FMT "' has %zu atoms, the loaded system %zu", MD_STR_ARG(filename), d.num_atoms, N);
        goto done;
    }
    // The atoms in the same order: every name matches. The residues may be cut otherwise in the structure (a PDB file
    // of several chains which reuse residue numbers), so their boundaries are compared only when the counts agree.
    for (size_t i = 0; i < N; ++i) {
        const md_str_t name = md_atom_name(&sys->atom, i);
        if (!topology_name_matches(name, d.atom_name[i])) {
            TOPO_LOG_ERROR(d.format, "atom %zu is '" MD_STR_FMT "' in the loaded system and '" MD_STR_FMT "' in '" MD_STR_FMT "'; the topology is not of these atoms in this order",
                i + 1, MD_STR_ARG(name), MD_STR_ARG(d.atom_name[i]), MD_STR_ARG(filename));
            goto done;
        }
    }
    if (sys->residue.count == d.num_residues) {
        size_t renamed = 0;
        for (size_t r = 0; r < d.num_residues; ++r) {
            const md_urange_t range = md_system_residue_atom_range(sys, r);
            if (range.beg != d.res_offset[r] || range.end != d.res_offset[r + 1]) {
                TOPO_LOG_ERROR(d.format, "residue %zu of '" MD_STR_FMT "' does not have the atoms of the loaded system's", r + 1, MD_STR_ARG(filename));
                goto done;
            }
            renamed += !topology_name_matches(md_system_residue_name(sys, r), d.res_name[r]);
        }
        if (renamed) {
            TOPO_LOG_INFO(d.format, "%zu residues are named otherwise in the topology (HIE for HIS, ...); the loaded names are kept", renamed);
        }
    }

    // ## Bonds: the topology's replace all but the user's, which are kept, as are the non-covalent ones
    {
        const size_t old_count = sys->bond.count;
        const size_t num_noncovalent = sys->bond.noncovalent_count;
        size_t num_user = 0;
        for (size_t i = 0; i < old_count; ++i) num_user += md_bond_origin(sys->bond.flags[i]) == MD_BOND_ORIGIN_USER;

        const size_t cap = d.num_bonds + num_user + num_noncovalent + 1;
        md_atom_pair_t*  pairs = md_temp_alloc_array(temp, md_atom_pair_t, cap);
        md_bond_flags_t* flags = md_temp_alloc_array(temp, md_bond_flags_t, cap);
        size_t count = 0;
        for (size_t b = 0; b < d.num_bonds; ++b) {
            const md_atom_pair_t pair = d.bonds[b];
            if (md_atom_particle_kind(&sys->atom, (size_t)pair.idx[0]) == MD_PARTICLE_VIRTUAL_SITE ||
                md_atom_particle_kind(&sys->atom, (size_t)pair.idx[1]) == MD_PARTICLE_VIRTUAL_SITE) continue;
            pairs[count] = pair;
            flags[count] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
            count += 1;
        }
        const size_t num_topology = count;
        for (size_t i = 0; i < old_count; ++i) {
            if (md_bond_origin(sys->bond.flags[i]) != MD_BOND_ORIGIN_USER) continue;
            pairs[count] = sys->bond.pairs[i];
            flags[count] = sys->bond.flags[i];
            count += 1;
        }
        if (num_noncovalent) {
            md_memcpy(pairs + count, sys->bond.pairs + old_count, num_noncovalent * sizeof(md_atom_pair_t));
            md_memcpy(flags + count, sys->bond.flags + old_count, num_noncovalent * sizeof(md_bond_flags_t));
        }
        md_bond_data_clear(&sys->bond);
        md_array_resize(sys->bond.pairs, count + num_noncovalent, sys->alloc);
        md_array_resize(sys->bond.flags, count + num_noncovalent, sys->alloc);
        if (count + num_noncovalent) {
            md_memcpy(sys->bond.pairs, pairs, (count + num_noncovalent) * sizeof(md_atom_pair_t));
            md_memcpy(sys->bond.flags, flags, (count + num_noncovalent) * sizeof(md_bond_flags_t));
        }
        sys->bond.count = count;
        sys->bond.noncovalent_count = num_noncovalent;
        md_bond_build_connectivity(&sys->bond, N, sys->alloc);
        md_util_system_infer_coordination(sys);
        md_system_topology_changed(sys);
        TOPO_LOG_INFO(d.format, "%zu bonds of the topology replace those inferred", num_topology);
    }

    // ## Charge and mass
    if (!sys->attributes.alloc) sys->attributes.alloc = sys->alloc;
    float* column = md_temp_alloc_array(temp, float, N);
    // Replacing what the structure file published, if anything
    md_attributes_remove_prefix(&sys->attributes, MD_STR_LIT("atom/charge"));
    md_attributes_remove_prefix(&sys->attributes, MD_STR_LIT("atom/mass"));
    if (d.charge) {
        for (size_t i = 0; i < N; ++i) column[i] = (float)d.charge[i];
        md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, column, N);
    }
    if (d.mass) {
        for (size_t i = 0; i < N; ++i) column[i] = (float)d.mass[i];
        md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/mass"), md_unit_dalton(), 1, column, N);
    }

    // ## Derived topology, from the new bonds and the coordinates the system was inferred from
    md_util_system_infer(sys, &sys->reference, MD_UTIL_INFER_PART_BIT | MD_UTIL_INFER_MOLECULE_BIT | MD_UTIL_INFER_BACKBONE_BIT, NULL);
    md_util_system_infer_lipids(sys);
    md_util_system_infer_residue_types(sys);
    result = true;

done:
    md_temp_end(temp);
    return result;
}
