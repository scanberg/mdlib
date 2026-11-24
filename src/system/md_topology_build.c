#include <md_topology_build.h>

#include <md_system.h>
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

// The elements of the atoms: those the topology gives (0 for none), else inferred from names and masses. Returns how
// many were inferred.
static size_t topology_elements(md_atomic_number_t* z, const md_topology_desc_t* d) {
    size_t num_inferred = 0;
    for (size_t r = 0; r < d->num_residues; ++r) {
        const size_t beg = d->res_offset[r], end = d->res_offset[r + 1];
        for (size_t i = beg; i < end; ++i) {
            const int64_t given = d->atomic_number ? d->atomic_number[i] : -1;
            if (given >= 0 && given <= 118) {
                z[i] = (md_atomic_number_t)given;   // 0: none, an extra point when massless
            } else {
                z[i] = topology_infer_element(d->atom_name[i], d->res_name[r], end - beg, d->mass ? d->mass[i] : 0.0);
                num_inferred += z[i] != 0;
            }
        }
    }
    return num_inferred;
}

// 'atom/charge' and 'atom/mass' of the topology, replacing what was published under them
static void topology_publish_quantities(md_system_t* sys, const md_topology_desc_t* d, float* column) {
    if (!sys->attributes.alloc) sys->attributes.alloc = sys->alloc;
    const size_t N = d->num_atoms;
    md_attributes_remove_prefix(&sys->attributes, MD_STR_LIT("atom/charge"));
    md_attributes_remove_prefix(&sys->attributes, MD_STR_LIT("atom/mass"));
    if (d->charge) {
        for (size_t i = 0; i < N; ++i) column[i] = (float)d->charge[i];
        md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, column, N);
    }
    if (d->mass) {
        for (size_t i = 0; i < N; ++i) column[i] = (float)d->mass[i];
        md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/mass"), md_unit_dalton(), 1, column, N);
    }
}

// The system of a topology, onto one without atoms
static void topology_build(md_system_t* sys, md_system_state_t* state, const md_topology_desc_t* desc) {
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    const md_topology_desc_t d = *desc;

    md_system_reset(sys);
    md_system_state_init(state, 0);
    md_allocator_i* alloc = sys->alloc;
    const size_t N = d.num_atoms, R = d.num_residues;

    // ## Elements
    md_atomic_number_t* z = md_temp_alloc_array(temp, md_atomic_number_t, N);
    const size_t num_inferred = topology_elements(z, &d);
    if (num_inferred) {
        TOPO_LOG_INFO(d.format, "the elements of %zu atoms were inferred from their names and masses", num_inferred);
    }

    // ## Particle types: one per name, force field type, element, mass and kind
    const size_t capacity = MD_ROUND_UP(N, 16);
    md_array_resize(sys->atom.type_idx, capacity, alloc);
    md_memset(sys->atom.type_idx, 0, capacity * sizeof(md_particle_type_idx_t));
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
    }
    md_array_push(sys->residue.atom_offset, (uint32_t)N, alloc);
    sys->residue.count = R;

    // ## Bonds
    {
        const md_bond_flags_t flags = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
        md_bond_list_t* given = &sys->given_bond;
        md_array_resize(given->pairs, d.num_bonds, alloc);
        md_array_resize(given->flags, d.num_bonds, alloc);
        size_t kept = 0;
        for (size_t b = 0; b < d.num_bonds; ++b) {
            const md_atom_pair_t pair = d.bonds[b];
            if (topology_is_virtual(&d, (size_t)pair.idx[0], z[pair.idx[0]]) || topology_is_virtual(&d, (size_t)pair.idx[1], z[pair.idx[1]])) continue;
            given->pairs[kept] = pair;
            given->flags[kept] = flags;
            kept += 1;
        }
        md_array_shrink(given->pairs, kept);
        md_array_shrink(given->flags, kept);
        // Of every atom: they have the topology's bonds and no others
        md_system_set_bonds_complete(sys);
    }

    md_system_augment_particle_types(sys);

    // ## Chains, as ParmEd keeps them from the structure the topology was made from
    if (d.res_chain) {
        bool any = false;
        for (size_t r = 0; r < R && !any; ++r) any = !md_str_empty(d.res_chain[r]);
        if (any) md_system_instances_from_chain_ids(sys, d.res_chain, NULL);
    }

    topology_publish_quantities(sys, &d, md_temp_alloc_array(temp, float, N));
    if (!md_str_empty(d.title)) {
        sys->description = md_alloc_str(alloc, d.title);
    }

    state->unitcell = d.cell;
    md_temp_end(temp);
}

md_load_result_t md_topology_system_load(md_system_t* sys, md_system_state_t* state, const md_topology_desc_t* desc, md_str_t filename) {
    MD_ASSERT(sys);
    MD_ASSERT(state);
    MD_ASSERT(desc);
    if (!sys->alloc || !state->alloc) {
        TOPO_LOG_ERROR(desc->format, "system or state allocator not set");
        return MD_LOAD_FAILED;
    }
    if (!desc->num_atoms || !desc->num_residues || !desc->atom_name || !desc->res_name || !desc->res_offset) {
        TOPO_LOG_ERROR(desc->format, "no atoms or residues");
        return MD_LOAD_FAILED;
    }
    if (sys->atom.count) {
        // A topology is a system of its own: joined to one of its atoms by md_system_merge
        TOPO_LOG_ERROR(desc->format, "'" MD_STR_FMT "' builds a system of its own; join it to a loaded one with md_system_merge", MD_STR_ARG(filename));
        return MD_LOAD_FAILED;
    }
    topology_build(sys, state, desc);
    return MD_LOAD_OK;
}
