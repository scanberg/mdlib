#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>

typedef struct md_atom_type_data_t {
    size_t count;

    md_label_t*     name;
    md_atomic_number_t* z;
    float*          mass;
    float*          radius;
    uint32_t*       color;
    uint32_t*       flags;
} md_atom_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates
    float* x;
    float* y;
    float* z;

    md_atom_type_idx_t* type_idx;
    md_flags_t* flags;

    md_atom_type_data_t type;
} md_atom_data_t;

// Component (Residue)
typedef struct md_component_data_t {
    size_t count;
    md_label_t* name;
    md_sequence_id_t* seq_id;
    uint32_t* atom_offset;
    md_flags_t* flags;
} md_component_data_t;

// Instance (e.g. Chains + other)
typedef struct md_instance_data_t {
    size_t count;
    md_label_t* id;
    md_label_t* auth_id;
    uint32_t* comp_offset;
    md_entity_idx_t* entity_idx;
} md_instance_data_t;

// Entities are used to describe the types found within a system
typedef struct md_entity_data_t {
    size_t count;
    md_label_t* id;
    md_flags_t* flags;
    str_t* description;
} md_entity_data_t;

// @TODO: Split this into two or more structucomp,
// One for nucleic backbone and one for protein backbone
typedef struct md_protein_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the segments
        md_instance_idx_t* inst_idx; // Reference to the instance in which the backbone is located
    } range;

    // These fields share the same length 'count'
    struct {
        size_t count;
        md_amino_acid_atoms_t* atoms;
        md_backbone_angles_t* angle;
        md_secondary_structure_t* secondary_structure;
        md_ramachandran_type_t* rama_type;
        md_component_idx_t* comp_idx;                  // Index to the component which contains the backbone
    } segment;
} md_protein_backbone_data_t;

typedef struct md_nucleic_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the backbone fields stored bellow
        md_instance_idx_t* inst_idx; // Reference to the instance in which the backbone is located
    } range;

    // These fields share the same length 'count'
    struct {
        size_t count;
        md_nucleic_acid_atoms_t* atoms;
        md_component_idx_t* comp_idx;                  // Index to the component which contains the backbone segment
    } segment;
} md_nucleic_backbone_data_t;

// This represents symmetries which are instanced, commonly found
// in PDB data. It is up to the renderer to properly render this instanced data.
typedef struct md_assembly_data_t {
    size_t count;
    md_urange_t* atom_range;
    md_label_t* label;
    mat4_t* transform;
} md_assembly_data_t;

// Atom centric representation of bonds
typedef struct md_bond_conn_data_t {
    size_t count;
    md_atom_idx_t* atom_idx; // Indices to the 'other' atoms
    md_bond_idx_t* bond_idx; // Indices to the bonds
    // The offsets into the atom_idx and bond_idx for each atom.
    // Consequently offset_count should be atom count + 1
    size_t offset_count;
    uint32_t* offset;
} md_bond_conn_data_t;

// Bond centric representation
typedef struct md_bond_data_t {
    size_t count;
    md_atom_pair_t* pairs;
    md_flags_t*     flags;
    md_bond_conn_data_t  conn;   // Connectivity
} md_bond_data_t;

typedef struct md_bond_iter_t {
    const md_bond_data_t* data;
	uint32_t i;
    uint32_t end_idx;
} md_bond_iter_t;

typedef struct md_hydrogen_bond_candidates_t {
    struct {
        size_t count;
        md_atom_idx_t* idx;
        int* num_lone_pairs;
    } acceptor;
   
    struct {
        size_t count;
        md_atom_idx_t* d_idx;
        md_atom_idx_t* h_idx;
    } donor;
} md_hydrogen_bond_candidates_t;

typedef struct md_hydrogen_bond_pair_t {
	uint32_t don_idx; // Index into hydrogen bond candidate donors
	uint32_t acc_idx; // Index into hydrogen bond candidate acceptors
} md_hydrogen_bond_pair_t;

typedef struct md_hydrogen_bond_data_t {
    md_hydrogen_bond_candidates_t candidate;

    size_t num_bonds;
    md_hydrogen_bond_pair_t* bonds;
} md_hydrogen_bond_data_t;

// Store atom coordinates in chunks of 4 to strike a good balance between cache locality for fetching individual atoms and SIMD processing
typedef struct md_coordinate_chunk_t {
    ALIGNAS(16) float x[4], y[4], z[4];
} md_coordinate_chunk_t;

// This represents the atom coordinates for a given state (frame) of the system.
typedef struct md_atom_coordinate_data_t {
    size_t num_chunks;
    md_coordinate_chunk_t* chunks;
} md_atom_coordinate_data_t;

// This represents the transient portion of a system which change over time, such as the atom coordinates, bonds, etc.
typedef struct md_system_state_t {
    md_atom_coordinate_data_t coord;
    md_unitcell_t unitcell;
    // Add other dynamic properties here such as e.g. hydrogen bonds, etc. which may change over time.
} md_system_state_t;

// This represents the persistent portion (topology) of a system which does not change over time, such as the atom types, bonds, components, etc.
// It may of course be modified through some special operations, though it is not expected to change frequently.
typedef struct md_system_t {
    md_unitcell_t               unitcell;

    md_atom_data_t              atom;
    md_component_data_t         component;
    md_instance_data_t          instance;
    md_entity_data_t            entity;

    md_protein_backbone_data_t  protein_backbone;
    md_nucleic_backbone_data_t  nucleic_backbone;
    
    md_bond_data_t              bond;               // Persistent covalent bonds
    md_hydrogen_bond_data_t     hydrogen_bond;      // Hydrogen bonds
    
    md_index_data_t             ring;               // Ring structures formed by persistent bonds
    md_index_data_t             structure;          // Isolated structures connected by persistent bonds

    md_assembly_data_t          assembly;           // Instances of the molecule (duplications of ranges with new transforms)

    str_t                       description;
} md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

// Atom coordinate helper functions

static inline vec3_t md_atom_coordinate_get_vec3(const md_atom_coordinate_data_t* coord, size_t atom_idx) {
    ASSERT(coord);
    size_t chunk_idx = atom_idx / 4;
    size_t within_chunk_idx = atom_idx % 4;
    if (chunk_idx < coord->num_chunks) {
        const md_coordinate_chunk_t* chunk = &coord->chunks[chunk_idx];
        return vec3_set(chunk->x[within_chunk_idx], chunk->y[within_chunk_idx], chunk->z[within_chunk_idx]);
    }
    return vec3_zero();
}

static inline vec4_t md_atom_coordinate_get_vec4(const md_atom_coordinate_data_t* coord, size_t atom_idx) {
    ASSERT(coord);
    size_t chunk_idx = atom_idx / 4;
    size_t within_chunk_idx = atom_idx % 4;
    if (chunk_idx < coord->num_chunks) {
        const md_coordinate_chunk_t* chunk = &coord->chunks[chunk_idx];
        return vec4_set(chunk->x[within_chunk_idx], chunk->y[within_chunk_idx], chunk->z[within_chunk_idx], 1.0f);
    }
    return vec4_zero();
}

static inline void md_atom_coordinate_set(md_atom_coordinate_data_t* coord, size_t atom_idx, float x, float y, float z) {
    ASSERT(coord);
    size_t chunk_idx = atom_idx / 4;
    size_t within_chunk_idx = atom_idx % 4;
    if (chunk_idx < coord->num_chunks) {
        md_coordinate_chunk_t* chunk = &coord->chunks[chunk_idx];
        chunk->x[within_chunk_idx] = x;
        chunk->y[within_chunk_idx] = y;
        chunk->z[within_chunk_idx] = z;
    }
}

void md_atom_coordinate_init(md_atom_coordinate_data_t* coord, size_t num_atoms, md_allocator_i* alloc);
void md_atom_coordinate_free(md_atom_coordinate_data_t* coord, md_allocator_i* alloc);

// Atom coordinate iterators (vec4_t, md_128, md_256)
// These are useful for iterating over atom coordinates in a SIMD friendly way

typedef struct md_atom_coordinate_iter_linear_t {
    const md_atom_coordinate_data_t* coord;
    size_t it;
    size_t count;
} md_atom_coordinate_iter_linear_t;

typedef struct md_atom_coordinate_iter_indexed_t {
    const md_atom_coordinate_data_t* coord;
    size_t count;
    const int32_t* idx;
    size_t it;
} md_atom_coordinate_iter_indexed_t;

static inline md_atom_coordinate_iter_linear_t md_atom_coordinate_iter_linear_create(const md_atom_coordinate_data_t* coord, size_t count) {
    ASSERT(coord);
    md_atom_coordinate_iter_linear_t iter = {
        .coord = coord,
        .it = 0,
        .count = count,
    };
    return iter;
}

static inline md_atom_coordinate_iter_indexed_t md_atom_coordinate_iter_indexed_create(const md_atom_coordinate_data_t* coord, const int32_t* idx, size_t count) {
    ASSERT(coord);
    ASSERT(idx);
    md_atom_coordinate_iter_indexed_t iter = {
        .coord = coord,
        .count = count,
        .idx = idx,
        .it = 0,
    };
    return iter;
}

static inline bool md_atom_coordinate_iter_linear_next_vec4(vec4_t* out_coord, float w, md_atom_coordinate_iter_linear_t* iter) {
    ASSERT(iter);
    ASSERT(out_coord);
    if (iter->it < iter->count) {
        size_t idx = iter->it;
        size_t chunk_idx = idx / 4;
        size_t within_chunk_idx = idx % 4;
        ASSERT(chunk_idx < iter->coord->num_chunks);
        const md_coordinate_chunk_t* chunk = &iter->coord->chunks[chunk_idx];
        *out_coord = vec4_set(chunk->x[within_chunk_idx], chunk->y[within_chunk_idx], chunk->z[within_chunk_idx], w);
        iter->it++;
        return true;
    }
    return false;
}

static inline bool md_atom_coordinate_iter_indexed_next_vec4(vec4_t* out_coord, float w, md_atom_coordinate_iter_indexed_t* iter) {
    ASSERT(iter);
    ASSERT(out_coord);
    if (iter->it < iter->count) {
        int idx = iter->idx[iter->it];
        int chunk_idx = idx / 4;
        int within_chunk_idx = idx % 4;
        ASSERT(0 <= chunk_idx && chunk_idx < iter->coord->num_chunks);
        const md_coordinate_chunk_t* chunk = &iter->coord->chunks[chunk_idx];
        *out_coord = vec4_set(chunk->x[within_chunk_idx], chunk->y[within_chunk_idx], chunk->z[within_chunk_idx], w);
        iter->it++;
        return true;
    }
    return false;
}

// Returns number of valid lanes [0..4]. Assumes coordinate data is padded so loads are always safe.
// The caller uses the returned lane count for tail handling.
static inline int md_atom_coordinate_iter_linear_next_chunk_128(md_128* out_x, md_128* out_y, md_128* out_z,
                                                                md_atom_coordinate_iter_linear_t* iter) {
    ASSERT(iter && out_x && out_y && out_z);

    if (iter->it >= iter->count) {
        return 0;
    }

    ASSERT((iter->it & 3) == 0); // chunk aligned (since offsets are disallowed)

    const size_t remaining = iter->count - iter->it;
    const int lanes = MIN(4, (int)remaining);

    const size_t chunk_idx = iter->it >> 2;
    ASSERT(chunk_idx < iter->coord->num_chunks);
    const md_coordinate_chunk_t* c = &iter->coord->chunks[chunk_idx];

    // Always safe due to padding (even for lanes < 4)
    *out_x = md_mm_load_ps(c->x);
    *out_y = md_mm_load_ps(c->y);
    *out_z = md_mm_load_ps(c->z);

    iter->it += 4;
    return lanes;
}

static inline int md_atom_coordinate_iter_linear_next_chunk_256(md_256* out_x, md_256* out_y, md_256* out_z,
                                                                md_atom_coordinate_iter_linear_t* iter) {
    ASSERT(iter && out_x && out_y && out_z);

    if (iter->it >= iter->count) {
        return 0;
    }

    ASSERT((iter->it & 7) == 0); // 8-wide aligned
    const size_t remaining = iter->count - iter->it;
    const int lanes = MIN(8, (int)remaining);

    // Two 128-bit chunks (each chunk is 4 atoms)
    const size_t chunk_idx0 = iter->it >> 2;        // /4
    const size_t chunk_idx1 = chunk_idx0 + 1;

    ASSERT(chunk_idx1 < iter->coord->num_chunks);   // requires padding for the +1 read
    const md_coordinate_chunk_t* c0 = &iter->coord->chunks[chunk_idx0];
    const md_coordinate_chunk_t* c1 = &iter->coord->chunks[chunk_idx1];

    // Equivalent to _mm256_loadu2_m128(high, low):
    // low  128 = c0.{x,y,z}, high 128 = c1.{x,y,z}
    const md_128 x0 = md_mm_load_ps(c0->x);
    const md_128 y0 = md_mm_load_ps(c0->y);
    const md_128 z0 = md_mm_load_ps(c0->z);

    const md_128 x1 = md_mm_load_ps(c1->x);
    const md_128 y1 = md_mm_load_ps(c1->y);
    const md_128 z1 = md_mm_load_ps(c1->z);

    *out_x = simde_mm256_set_m128(x1, x0);
    *out_y = simde_mm256_set_m128(y1, y0);
    *out_z = simde_mm256_set_m128(z1, z0);

    iter->it += 8;
    return lanes;
}

static inline int md_atom_coordinate_iter_indexed_next_chunk_128(md_128* out_x, md_128* out_y, md_128* out_z,
                                                                 md_atom_coordinate_iter_indexed_t* iter) {
    ASSERT(iter && out_x && out_y && out_z);

    if (iter->it >= iter->count) return 0;

    const int remaining = (int)iter->count - (int)iter->it;
    const int lanes = MIN(4, remaining);

    // Base for gathers: treat the entire `chunks` array as a flat float array.
    const float* base = (const float*)&iter->coord->chunks[0].x[0];
    ASSERT(sizeof(md_coordinate_chunk_t) / sizeof(float) == 12); // x[4],y[4],z[4] => 12 floats

    md_128i a;

    if (remaining >= 4) {
        // Fast path: we can safely read 4 indices
        a = md_mm_loadu_si128((const md_128i*)(iter->idx + iter->it));
    } else {
        // Tail: construct a 4-wide index vector using a sentinel atom index in unused lanes.
        const size_t last_idx = iter->count - 1;
        a = md_mm_set_epi32(
            iter->idx[MIN(iter->it + 3, last_idx)],
            iter->idx[MIN(iter->it + 2, last_idx)],
            iter->idx[MIN(iter->it + 1, last_idx)],
            iter->idx[iter->it]
        );
    }

    // base offset is computed using shifts and additions rather than multiplication
    md_128i lane = simde_mm_and_si128(a, simde_mm_set1_epi32(3));

    // three_a = 3*a = a + 2*a
    md_128i two_a   = simde_mm_slli_epi32(a, 1);
    md_128i three_a = simde_mm_add_epi32(two_a, a);

    // base_off = 3*a - 2*lane
    md_128i two_lane = simde_mm_slli_epi32(lane, 1);
    md_128i base_off = simde_mm_sub_epi32(three_a, two_lane);

    md_128i off_x = simde_mm_add_epi32(base_off, simde_mm_set1_epi32(0));
    md_128i off_y = simde_mm_add_epi32(base_off, simde_mm_set1_epi32(4));
    md_128i off_z = simde_mm_add_epi32(base_off, simde_mm_set1_epi32(8));

    *out_x = md_mm_i32gather_ps(base, off_x, 4);
    *out_y = md_mm_i32gather_ps(base, off_y, 4);
    *out_z = md_mm_i32gather_ps(base, off_z, 4);

    iter->it += lanes;
    return lanes;
}

// Atom type table helper functions
static inline size_t md_atom_type_count(const md_atom_type_data_t* atom_type) {
    ASSERT(atom_type);
    return atom_type->count;
}

static inline md_atom_type_idx_t md_atom_type_find(const md_atom_type_data_t* atom_type, str_t name, md_atomic_number_t z) {
    ASSERT(atom_type);
	md_atom_type_idx_t type_idx = 0; // Zero is sentinel for "not found"
    for (size_t i = 0; i < atom_type->count; ++i) {
        str_t atom_type_name = LBL_TO_STR(atom_type->name[i]);
        if (str_eq(atom_type_name, name) && atom_type->z[i] == z) {
            type_idx = (md_atom_type_idx_t)i;
            break;
        }
    }
    return type_idx;
}

static inline md_atom_type_idx_t md_atom_type_find_or_add(md_atom_type_data_t* atom_type, str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_flags_t flags, struct md_allocator_i* alloc) {
    ASSERT(atom_type);
    ASSERT(alloc);
    
    // First try to find existing atom type
    md_atom_type_idx_t type_idx = md_atom_type_find(atom_type, name, z);
    if (type_idx != 0) {
        return type_idx;
    }
    
    // Add new atom type
    md_array_push(atom_type->name, make_label(name), alloc);
    md_array_push(atom_type->z, z, alloc);
    md_array_push(atom_type->mass, mass, alloc);
    md_array_push(atom_type->radius, radius, alloc);
    md_array_push(atom_type->color, color, alloc);
    md_array_push(atom_type->flags, flags, alloc);
    atom_type->count++;
    
    return (md_atom_type_idx_t)(atom_type->count - 1);
}

static inline md_atomic_number_t md_atom_type_atomic_number(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->z[type_idx];
    }
    return 0;
}

static inline float md_atom_type_mass(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->mass[type_idx];
    }
    return 0;
}

static inline md_flags_t md_atom_type_flags(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->flags[type_idx];
    }
    return 0;
}

static inline float md_atom_type_radius(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->radius[type_idx];
    }
    return 0;
}

static inline uint32_t md_atom_type_color(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->color[type_idx];
    }
    return 0;
}

static inline str_t md_atom_type_name(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return LBL_TO_STR(type_data->name[type_idx]);
    }
    return STR_LIT("");
}

// Atom helpers
static inline size_t md_atom_count(const md_atom_data_t* atom_data) {
    ASSERT(atom_data);
    return atom_data->count;
}

static inline vec3_t md_atom_coord(const md_atom_data_t* atom_data, size_t atom_idx) {
    ASSERT(atom_data);
    return vec3_set(atom_data->x[atom_idx], atom_data->y[atom_idx], atom_data->z[atom_idx]);
}

static inline md_atom_type_idx_t md_atom_type_idx(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom->type_idx && atom_idx < atom->count) {
        return atom->type_idx[atom_idx];
	}
    return 0;
}

static inline md_atomic_number_t md_atom_atomic_number(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    
    // Try atom type table first if type_idx is available
    if (atom->type_idx && (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return atom->type.z[atom->type_idx[atom_idx]];
    }
    
    return 0;
}

static inline float md_atom_mass(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    
    if (atom_idx < atom->count) {
		return md_atom_type_mass(&atom->type, atom->type_idx[atom_idx]);
    }
    
    return 0.0f;
}

static inline float md_atom_radius(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    
    if (atom_idx < atom->count) {
		return md_atom_type_radius(&atom->type, atom->type_idx[atom_idx]);
    }
    
    return 0.0f;
}

static inline str_t md_atom_name(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);

    if (atom_idx < atom->count) {
		return md_atom_type_name(&atom->type, atom->type_idx[atom_idx]);
    }

    return STR_LIT("");
}

static inline md_flags_t md_atom_flags(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    if (atom_idx < atom->count && atom->flags) {
        return atom->flags[atom_idx];
	}
    return 0;
}

// Component

static inline size_t md_component_count(const md_component_data_t* comp) {
    ASSERT(comp);
    return comp->count;
}

static inline str_t md_component_name(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    str_t name = STR_LIT("");
    if (comp->name && comp_idx < comp->count) {
        name = LBL_TO_STR(comp->name[comp_idx]);
    }
    return name;
}

static inline md_sequence_id_t md_component_seq_id(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    md_sequence_id_t id = 0;
    if (comp->seq_id && comp_idx < comp->count) {
        id = comp->seq_id[comp_idx];
    }
    return id;
}

static inline md_urange_t md_component_atom_range(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
	md_urange_t range = {0};
	if (comp->atom_offset && comp_idx < comp->count) {
		range.beg = comp->atom_offset[comp_idx];
		range.end = comp->atom_offset[comp_idx + 1];
	}
	return range;
}

static inline md_flags_t md_component_flags(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    md_flags_t flags = 0;
    if (comp->flags && comp_idx < comp->count) {
        flags = comp->flags[comp_idx];
    }
    return flags;
}

static inline md_component_idx_t md_component_find_by_atom_idx(const md_component_data_t* comp, size_t atom_idx) {
    ASSERT(comp);

    md_component_idx_t comp_idx = -1;
    if (comp->atom_offset) {
        for (size_t i = 0; i < comp->count; ++i) {
            size_t comp_beg = comp->atom_offset[i];
            size_t comp_end = comp->atom_offset[i + 1];
            if (comp_beg <= atom_idx && atom_idx < comp_end) {
                comp_idx = (md_component_idx_t)i;
                break;
            }
            if (comp_beg > atom_idx) {
                break;
            }
        }
    }

    return comp_idx;
}

static inline size_t md_component_atom_count(const md_component_data_t* comp, size_t comp_idx) {
    ASSERT(comp);
    size_t count = 0;

    if (comp->atom_offset && comp_idx < comp->count) {
        count = comp->atom_offset[comp_idx + 1] - comp->atom_offset[comp_idx];
    }
    return count;
}

// Instance

static inline size_t md_instance_count(const md_instance_data_t* inst) {
    ASSERT(inst);
    return inst->count;
}

static inline md_urange_t md_instance_component_range(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);

    md_urange_t range = {0};
    if (inst->comp_offset && inst_idx < inst->count) {
        range.beg = inst->comp_offset[inst_idx];
        range.end = inst->comp_offset[inst_idx + 1];
    }
    return range;
}

static inline md_instance_idx_t md_instance_find_by_comp_idx(const md_instance_data_t* inst, size_t comp_idx) {
    ASSERT(inst);

    md_instance_idx_t inst_idx = -1;
    if (inst->comp_offset) {
        for (size_t i = 0; i < inst->count; ++i) {
            size_t inst_beg = inst->comp_offset[i];
            size_t inst_end = inst->comp_offset[i + 1];
            if (inst_beg <= comp_idx && comp_idx < inst_end) {
                inst_idx = (md_instance_idx_t)i;
                break;
            }
            if (inst_beg > i) {
                break;
            }
        }
    }
    return inst_idx;
}

/*
static inline md_instance_idx_t md_inst_find_by_atom_idx(const md_instance_data_t* inst, size_t atom_idx) {
    ASSERT(inst);

    md_instance_idx_t inst_idx = -1;
    if (inst->atom_range) {
        int ai = (int)atom_idx;
        for (size_t i = 0; i < inst->count; ++i) {
            md_range_t range = inst->atom_range[i];
            if (range.beg <= ai && ai < range.end) {
                inst_idx = (md_instance_idx_t)i;
                break;
            }
            if (range.beg > ai) {
                break;
            }
        }
    }
    return inst_idx;
}
*/

static inline size_t md_instance_comp_count(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);

    size_t count = 0;
    if (inst->comp_offset && inst_idx < inst->count) {
        count = inst->comp_offset[inst_idx + 1] - inst->comp_offset[inst_idx];
    }
    return count;
}

/*
static inline md_urange_t md_inst_atom_range(const md_instance_data_t* inst, size_t inst_idx) {
	ASSERT(inst);

    md_range_t range = {0};
    if (inst->comp_range && inst_idx < inst->count) {
        uint32_t cbeg = inst->comp_range[inst_idx].beg;
        uint32_t cend = inst->comp_range[inst_idx].end;
        if (inst->atom_range && cbeg < cend) {
            range.beg = inst->atom_range[cbeg].beg;
            range.end = inst->atom_range[cend - 1].end;
        }
        range = inst->atom_range[inst_idx];
    }
    return range;
}

static inline size_t md_inst_atom_count(const md_instance_data_t* inst, size_t inst_idx) {
    size_t count = 0;
    if (inst->atom_range && inst_idx < inst->count) {
        md_range_t range = inst->atom_range[inst_idx];
        count = range.end - range.beg;
    }
    return count;
}
*/

static inline str_t md_instance_id(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    str_t id = STR_LIT("");
    if (inst->id && inst_idx < inst->count) {
        id = LBL_TO_STR(inst->id[inst_idx]);
    }
    return id;
}

static inline str_t md_instance_auth_id(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    str_t auth_id = STR_LIT("");
    if (inst->auth_id && inst_idx < inst->count) {
        auth_id = LBL_TO_STR(inst->auth_id[inst_idx]);
    }
    return auth_id;
}

static inline md_entity_idx_t md_instance_entity_idx(const md_instance_data_t* inst, size_t inst_idx) {
    ASSERT(inst);
    md_entity_idx_t entity_idx = -1;
    if (inst->entity_idx && inst_idx < inst->count) {
        entity_idx = inst->entity_idx[inst_idx];
    }
    return entity_idx;
}

static inline size_t md_entity_count(const md_entity_data_t* entity) {
    ASSERT(entity);
    return entity->count;
}

static inline md_entity_idx_t md_entity_find_by_id(const md_entity_data_t* entity, str_t id) {
    ASSERT(entity);
    md_entity_idx_t entity_idx = -1;
    if (entity->id) {
        for (size_t i = 0; i < entity->count; ++i) {
            str_t entity_id = LBL_TO_STR(entity->id[i]);
            if (str_eq(entity_id, id)) {
                entity_idx = (md_entity_idx_t)i;
                break;
            }
        }
    }
    return entity_idx;
}

static inline str_t md_entity_id(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    str_t label = STR_LIT("");
    if (entity->id && entity_idx < entity->count) {
        label = LBL_TO_STR(entity->id[entity_idx]);
    }
    return label;
}

static inline str_t md_entity_description(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    str_t desc = STR_LIT("");
    if (entity->description && entity_idx < entity->count) {
        desc = entity->description[entity_idx];
    }
    return desc;
}

static inline md_flags_t md_entity_flags(const md_entity_data_t* entity, size_t entity_idx) {
    ASSERT(entity);
    md_flags_t flags = 0;
    if (entity->flags && entity_idx < entity->count) {
        flags = entity->flags[entity_idx];
    }
    return flags;
}

static inline size_t md_structure_count(const md_index_data_t* structure) {
    ASSERT(structure);
    return md_index_data_num_ranges(structure);
}

static inline size_t md_structure_atom_count(const md_index_data_t* structure, size_t struct_idx) {
    ASSERT(structure);
    return md_index_range_size(structure, struct_idx);
}

static inline const md_atom_idx_t* md_structure_atom_indices(const md_index_data_t* structure, size_t struct_idx) {
    ASSERT(structure);
    return md_index_range_beg(structure, struct_idx);
}


// SYSTEM
// System level convenience accessors

static inline size_t md_system_atom_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_atom_count(&sys->atom);
}

static inline md_flags_t md_system_atom_flags(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    return md_atom_flags(&sys->atom, atom_idx);
}

static inline size_t md_system_atom_type_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_atom_type_count(&sys->atom.type);
}

static inline size_t md_system_component_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_component_count(&sys->component);
}

static inline md_flags_t md_system_component_flags(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_flags(&sys->component, comp_idx);
}

static inline size_t md_system_instance_count(const md_system_t* sys) {
    ASSERT(sys);
    return md_instance_count(&sys->instance);
}

static inline size_t md_system_bond_count(const md_system_t* sys) {
    ASSERT(sys);
    return sys->bond.count;
}

static inline size_t md_system_entity_count(const md_system_t* sys) {
    ASSERT(sys);
    return sys->entity.count;
}

static inline md_flags_t md_system_entity_flags(const md_system_t* sys, size_t ent_idx) {
    ASSERT(sys);
    return md_entity_flags(&sys->entity, ent_idx);
}

static inline md_flags_t md_system_instance_flags(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    if (sys->instance.entity_idx && inst_idx < sys->instance.count) {
        return md_entity_flags(&sys->entity, sys->instance.entity_idx[inst_idx]);
    }
    return 0;
}

static inline str_t md_system_instance_id(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    str_t id = STR_LIT("");
    if (sys->instance.id && inst_idx < sys->instance.count) {
        id = md_instance_id(&sys->instance, inst_idx);
    }
    return id;
}

static inline str_t md_system_instance_auth_id(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    str_t id = STR_LIT("");
    if (sys->instance.id && inst_idx < sys->instance.count) {
        id = md_instance_auth_id(&sys->instance, inst_idx);
    }
    return id;
}

static inline size_t md_system_instance_comp_count(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    return md_instance_comp_count(&sys->instance, inst_idx);
}

static inline md_urange_t md_system_instance_comp_range(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    return md_instance_component_range(&sys->instance, inst_idx);
}

static inline md_urange_t md_system_instance_atom_range(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    md_urange_t range = {0};
    if (inst_idx < sys->instance.count) {
        md_urange_t comp_range = md_instance_component_range(&sys->instance, inst_idx);
        if (comp_range.beg != comp_range.end) {
            range.beg = md_component_atom_range(&sys->component, comp_range.beg).beg;
            range.end = md_component_atom_range(&sys->component, comp_range.end - 1).end;
        }
    }
    return range;
}

static inline size_t md_system_instance_atom_count(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    md_urange_t atom_range = md_system_instance_atom_range(sys, inst_idx);
    return atom_range.end - atom_range.beg;
}

static inline size_t md_system_instance_entity_idx(const md_system_t* sys, size_t inst_idx) {
    ASSERT(sys);
    return md_instance_entity_idx(&sys->instance, inst_idx);
}

static inline md_sequence_id_t md_system_component_seq_id(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_seq_id(&sys->component, comp_idx);
}

static inline md_urange_t md_system_component_atom_range(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_atom_range(&sys->component, comp_idx);
}

static inline size_t md_system_component_atom_count(const md_system_t* sys, size_t comp_idx) {
    ASSERT(sys);
    return md_component_atom_count(&sys->component, comp_idx);
}

static inline md_component_idx_t md_system_component_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    return md_component_find_by_atom_idx(&sys->component, atom_idx);
}

static inline md_instance_idx_t md_system_instance_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    ASSERT(sys);
    md_instance_idx_t inst_idx = -1;
    md_component_idx_t comp_idx = md_system_component_find_by_atom_idx(sys, atom_idx);
    if (comp_idx >= 0) {
        inst_idx = md_instance_find_by_comp_idx(&sys->instance, comp_idx);
    }
    return inst_idx;
}

// Convenience functions to extract atom properties into arrays
static inline void md_atom_extract_radii(float out_radii[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    ASSERT(out_radii);
    ASSERT(atom_data);
    ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_radii[i] = md_atom_radius(atom_data, offset + i);
    }
}

static inline void md_atom_extract_masses(float out_masses[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    ASSERT(out_masses);
    ASSERT(atom_data);
    ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_masses[i] = md_atom_mass(atom_data, offset + i);
    }
}

static inline void md_atom_extract_atomic_numbers(md_atomic_number_t out_z[], size_t offset, size_t length, const md_atom_data_t* atom_data) {
    ASSERT(out_z);
    ASSERT(atom_data);
    ASSERT(offset + length <= atom_data->count);
    
    for (size_t i = 0; i < length; ++i) {
        out_z[i] = md_atom_atomic_number(atom_data, offset + i);
    }
}

static inline md_bond_iter_t md_bond_iter(const md_bond_data_t* bond_data, size_t atom_idx) {
    md_bond_iter_t it = {0};
    if (bond_data && bond_data->conn.offset && atom_idx < bond_data->conn.offset_count) {
        it.data = bond_data;
        it.i = bond_data->conn.offset[atom_idx];
		it.end_idx = bond_data->conn.offset[atom_idx + 1];
    }
    return it;
}

// This is not something which should be done frequently
void md_bond_build_connectivity(md_bond_data_t* in_out_bond, size_t atom_count, md_allocator_i* alloc);
void md_system_bond_build_connectivity(md_system_t* sys, md_allocator_i* alloc);

static inline size_t md_bond_conn_count(const md_bond_data_t* bond_data, size_t atom_idx) {
    ASSERT(bond_data);
    return bond_data->conn.offset[atom_idx + 1] - bond_data->conn.offset[atom_idx];
}

static inline md_atom_idx_t md_bond_conn_atom_idx(const md_bond_data_t* bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    ASSERT(bond_data);
    return bond_data->conn.atom_idx[atom_conn_idx + idx];
}

static inline md_bond_idx_t md_bond_conn_bond_idx(const md_bond_data_t* bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    ASSERT(bond_data);
    return bond_data->conn.bond_idx[atom_conn_idx + idx];
}

static inline bool md_bond_iter_has_next(const md_bond_iter_t* it) {
    ASSERT(it);
    return it->i < it->end_idx;
}

static inline void md_bond_iter_next(md_bond_iter_t* it) {
    ASSERT(it);
    it->i += 1;
}

static inline md_atom_idx_t md_bond_iter_atom_index(const md_bond_iter_t* it) {
    ASSERT(it);
	return it->data->conn.atom_idx[it->i];
}

static inline md_atom_idx_t md_bond_iter_bond_index(const md_bond_iter_t* it) {
    ASSERT(it);
    return it->data->conn.bond_idx[it->i];
}

static inline uint32_t md_bond_iter_bond_flags(const md_bond_iter_t* it) {
    ASSERT(it);
    return it->data->flags[it->data->conn.bond_idx[it->i]];
}

static inline md_bond_idx_t md_bond_find(const md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    ASSERT(bond_data);
	md_bond_iter_t it = md_bond_iter(bond_data, atom_idx_a);
    while (md_bond_iter_has_next(&it)) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(&it);
        if (other_atom_idx == atom_idx_b) {
            return md_bond_iter_bond_index(&it);
        }
        md_bond_iter_next(&it);
	}
	return -1;
}

static inline void md_bond_insert(md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b, md_flags_t flags, md_allocator_i* alloc) {
    ASSERT(bond_data);
    ASSERT(alloc);

    // Ensure that the bond does not exist
	md_bond_idx_t bond_idx = md_bond_find(bond_data, atom_idx_a, atom_idx_b);
    if (bond_idx != -1) {
        return;
    }

    // Add bond
    md_atom_pair_t pair = {atom_idx_a, atom_idx_b};
    md_array_push(bond_data->pairs, pair, alloc);
    md_array_push(bond_data->flags, flags, alloc);
    bond_data->count++;
}

// This requires a rebuild of connectivity to be valid again
static inline void md_bond_remove(md_bond_data_t* bond_data, md_bond_idx_t bond_idx) {
    ASSERT(bond_data);
    ASSERT(bond_idx < (md_bond_idx_t)bond_data->count);

    // Swap and pop
    size_t last_idx = bond_data->count - 1;
    if ((size_t)bond_idx != last_idx) {
        bond_data->pairs[bond_idx] = bond_data->pairs[last_idx];
        bond_data->flags[bond_idx] = bond_data->flags[last_idx];
    }
    bond_data->count--;
    md_array_shrink(bond_data->pairs, bond_data->count);
    md_array_shrink(bond_data->flags, bond_data->count);
}

static inline md_atom_pair_t md_bond_pair(const md_bond_data_t* bond_data, md_bond_idx_t bond_idx) {
    ASSERT(bond_data);
    if (bond_idx < (md_bond_idx_t)bond_data->count) {
        return bond_data->pairs[bond_idx];
    }
    md_atom_pair_t invalid_pair = { -1, -1 };
    return invalid_pair;
}

static inline md_bond_idx_t md_system_bond_find(const md_system_t* sys, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    ASSERT(sys);
    return md_bond_find(&sys->bond, atom_idx_a, atom_idx_b);
}

static inline void md_system_bond_insert(md_system_t* sys, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b, md_flags_t flags, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);
    md_bond_insert(&sys->bond, atom_idx_a,  atom_idx_b, flags, alloc);
}

static inline void md_system_bond_remove(md_system_t* sys, md_bond_idx_t bond_idx) {
    ASSERT(sys);
    md_bond_remove(&sys->bond, bond_idx);
}

static inline void md_bond_conn_clear(md_bond_conn_data_t* conn_data) {
    ASSERT(conn_data);
    conn_data->count = 0;
    md_array_shrink(conn_data->atom_idx, 0);
    md_array_shrink(conn_data->bond_idx, 0);

    conn_data->offset_count = 0;
    md_array_shrink(conn_data->offset, 0);
}

static inline void md_bond_data_clear(md_bond_data_t* bond_data) {
    ASSERT(bond_data);

    bond_data->count = 0;
    md_array_shrink(bond_data->pairs, 0);
    md_array_shrink(bond_data->flags, 0);
    
    bond_data->conn.count = 0;
    md_array_shrink(bond_data->conn.atom_idx, 0);
    md_array_shrink(bond_data->conn.bond_idx, 0);

    md_bond_conn_clear(&bond_data->conn);
}

static inline bool md_atom_is_connected_to_atomic_numbers(const md_atom_data_t* atom_data, const md_bond_data_t* bond_data, size_t atom_idx, const md_atomic_number_t z_list[], size_t z_count) {
    ASSERT(bond_data);
    ASSERT(atom_data);
    ASSERT(atom_idx < atom_data->count);
    bool found = false;
    md_bond_iter_t it = md_bond_iter(bond_data, atom_idx);
    while (md_bond_iter_has_next(&it) && !found) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(&it);
        md_atomic_number_t other_z = md_atom_atomic_number(atom_data, other_atom_idx);
        for (size_t i = 0; i < z_count; ++i) {
            if (other_z == z_list[i]) {
                found = true;
                break;
            }
        }
        md_bond_iter_next(&it);
    }
    return found;
}

static inline md_atom_idx_t md_hydrogen_bond_donor_atom_idx(const md_hydrogen_bond_candidates_t* hbond_cand, size_t don_idx) {
	ASSERT(hbond_cand);
    if (don_idx < hbond_cand->donor.count) {
        return hbond_cand->donor.d_idx[don_idx];
	}
    return -1;
}

static inline md_atom_idx_t md_hydrogen_bond_donor_hydrogen_atom_idx(const md_hydrogen_bond_candidates_t* hbond_cand, size_t don_idx) {
    ASSERT(hbond_cand);
    if (don_idx < hbond_cand->donor.count) {
        return hbond_cand->donor.h_idx[don_idx];
    }
    return -1;
}

static inline md_atom_idx_t md_hydrogen_bond_acceptor_atom_idx(const md_hydrogen_bond_candidates_t* hbond_cand, size_t acc_idx) {
    ASSERT(hbond_cand);
    if (acc_idx < hbond_cand->acceptor.count) {
        return hbond_cand->acceptor.idx[acc_idx];
    }
    return -1;
}

static inline int md_hydrogen_bond_acceptor_num_lone_pairs(const md_hydrogen_bond_candidates_t* hbond_cand, size_t acc_idx) {
    ASSERT(hbond_cand);
    if (acc_idx < hbond_cand->acceptor.count) {
        return hbond_cand->acceptor.num_lone_pairs[acc_idx];
    }
    return 0;
}

/*

The molecule loader interface is just a convenience interface for abstracing the functionality of initializing molecule data

The reason for providing a distinct function for initializing from file is that some molecule files can
also contain their trajectories, such as PDB files. In such case, the whole file would have to be read and passed, but for
molecule data only the first part of the file is used.

*/

typedef struct md_system_loader_i {
    bool (*init_from_str) (md_system_t* sys, str_t string,   const void* arg, struct md_allocator_i* alloc);
    bool (*init_from_file)(md_system_t* sys, str_t filename, const void* arg, struct md_allocator_i* alloc);
} md_system_loader_i;

// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
void md_system_free(md_system_t* sys, struct md_allocator_i* alloc);
void md_system_copy(md_system_t* dst_sys, const md_system_t* src_sys, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
