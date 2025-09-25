#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>

typedef struct md_atom_type_data_t {
    size_t count;

    md_label_t*   name;
    md_atomic_number_t* z;
    float*        mass;
    float*        radius;
} md_atom_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates
    float* x;
    float* y;
    float* z;

    // Atom Type Index (NEW: references into atom_type table)
    md_atom_type_idx_t* type_idx;
    md_flags_t* flags;

    md_atom_type_data_t type;
} md_atom_data_t;

typedef struct md_residue_data_t {
    size_t count;
    md_label_t* name;
    md_residue_id_t* id;
    uint32_t* atom_offset;
    md_flags_t* flags;
} md_residue_data_t;

typedef struct md_chain_data_t {
    size_t count;
    md_label_t* id;
    md_range_t* res_range;
    md_range_t* atom_range;
} md_chain_data_t;

// @TODO: Split this into two or more structures,
// One for nucleic backbone and one for protein backbone
typedef struct md_protein_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the backbone fields stored bellow
        md_chain_idx_t* chain_idx; // Reference to the chain in which the backbone is located
    } range;

    // These fields share the same length 'count'
    size_t count;
    md_protein_backbone_atoms_t* atoms;
    md_backbone_angles_t* angle;
    md_secondary_structure_t* secondary_structure;
    md_ramachandran_type_t* ramachandran_type;
    md_residue_idx_t* residue_idx;                  // Index to the residue which contains the backbone
} md_protein_backbone_data_t;

typedef struct md_nucleic_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the backbone fields stored bellow
        md_chain_idx_t* chain_idx; // Reference to the chain in which the backbone is located
    } range;

    // These fields share the same length 'count'
    size_t count;
    md_nucleic_backbone_atoms_t* atoms;
    md_residue_idx_t* residue_idx;                  // Index to the residue which contains the backbone
} md_nucleic_backbone_data_t;

// This represents symmetries which are instanced, commonly found
// in PDB data. It is up to the renderer to properly render this instanced data.
typedef struct md_instance_data_t {
    size_t count;
    md_range_t* atom_range;
    md_label_t* label;
    mat4_t* transform;
} md_instance_data_t;

// Atom centric representation of bonds
typedef struct md_conn_data_t {
    size_t count;
    md_atom_idx_t* atom_idx; // Indices to the 'other' atoms
    md_bond_idx_t* bond_idx; // Indices to the bonds
    // The offsets into the atom_idx and bond_idx for each atom.
    // Consequently offset_count should be atom count + 1
    size_t offset_count;
    uint32_t* offset;
} md_conn_data_t;

// Bond centric representation
typedef struct md_bond_data_t {
    size_t count;
    md_bond_pair_t* pairs;
    md_order_t*     order;
    md_conn_data_t  conn;   // Connectivity
} md_bond_data_t;

typedef struct md_bond_iter_t {
    const md_bond_data_t* data;
	uint32_t i;
    uint32_t end_idx;
} md_bond_iter_t;

typedef struct md_molecule_t {
    md_unit_cell_t              unit_cell;
    md_atom_data_t              atom;
    md_residue_data_t           residue;
    md_chain_data_t             chain;
    md_protein_backbone_data_t  protein_backbone;
    md_nucleic_backbone_data_t  nucleic_backbone;
    
    md_bond_data_t              bond;               // Persistent covalent bonds
    
    md_index_data_t             ring;               // Ring structures formed by persistent bonds
    md_index_data_t             structure;          // Isolated structures connected by persistent bonds

    md_instance_data_t          instance;           // Instances of the molecule (duplications of ranges with new transforms)
} md_molecule_t;

#ifdef __cplusplus
extern "C" {
#endif

// Atom type table helper functions
static inline md_atom_type_idx_t md_atom_type_find_or_add(md_atom_type_data_t* atom_type, str_t name, md_atomic_number_t z, float mass, float radius, struct md_allocator_i* alloc) {
    ASSERT(atom_type);
    ASSERT(alloc);
    
    // First try to find existing atom type
    for (size_t i = 0; i < atom_type->count; ++i) {
        str_t atom_type_name = LBL_TO_STR(atom_type->name[i]);
        if (str_eq(atom_type_name, name) && 
            atom_type->z[i] == z &&
            atom_type->mass[i] == mass &&
            atom_type->radius[i] == radius) {
            return (md_atom_type_idx_t)i;
        }
    }
    
    // Add new atom type
    md_array_push(atom_type->name, make_label(name), alloc);
    md_array_push(atom_type->z, z, alloc);
    md_array_push(atom_type->mass, mass, alloc);
    md_array_push(atom_type->radius, radius, alloc);
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

static inline float md_atom_type_radius(const md_atom_type_data_t* type_data, size_t type_idx) {
    ASSERT(type_data);
    if (type_idx < type_data->count) {
        return type_data->radius[type_idx];
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

static inline vec3_t md_atom_coord(const md_atom_data_t* atom_data, size_t atom_idx) {
    ASSERT(atom_data);
    return vec3_set(atom_data->x[atom_idx], atom_data->y[atom_idx], atom_data->z[atom_idx]);
}

static inline md_atom_type_idx_t md_atom_type_idx(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    ASSERT(atom_idx < atom->count);
    return atom->type_idx[atom_idx];
}

static inline md_atomic_number_t md_atom_atomic_number(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    ASSERT(atom_idx < atom->count);
    
    // Try atom type table first if type_idx is available
    if (atom->type_idx && atom->type_idx[atom_idx] >= 0 && 
        (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return atom->type.z[atom->type_idx[atom_idx]];
    }
    
    return 0;
}

static inline float md_atom_mass(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    ASSERT(atom_idx < atom->count);
    
    // Try atom type table first if type_idx is available
    if (atom->type_idx && atom->type_idx[atom_idx] >= 0 && 
        (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return atom->type.mass[atom->type_idx[atom_idx]];
    }
    
    return 0.0f;
}

static inline float md_atom_radius(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    ASSERT(atom_idx < atom->count);
    
    if (atom->type_idx && atom->type_idx[atom_idx] >= 0 && 
        (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return atom->type.radius[atom->type_idx[atom_idx]];
    }
    
    return 0.0f;
}

static inline str_t md_atom_name(const md_atom_data_t* atom, size_t atom_idx) {
    ASSERT(atom);
    ASSERT(atom_idx < atom->count);

    if (atom->type_idx && atom->type_idx[atom_idx] >= 0 && 
        (size_t)atom->type_idx[atom_idx] < atom->type.count) {
        return LBL_TO_STR(atom->type.name[atom->type_idx[atom_idx]]);
    }

    return STR_LIT("");
}

static inline str_t md_residue_name(const md_residue_data_t* res, size_t res_idx) {
    ASSERT(res);
    str_t name = STR_LIT("");
    if (res->name && res_idx < res->count) {
        name = LBL_TO_STR(res->name[res_idx]);
    }
    return name;
}

static inline md_residue_id_t md_residue_id(const md_residue_data_t* res, size_t res_idx) {
    ASSERT(res);
    md_residue_id_t id = 0;
    if (res->id && res_idx < res->count) {
        id = res->id[res_idx];
    }
    return id;
}

static inline md_range_t md_residue_atom_range(const md_residue_data_t* res, size_t res_idx) {
    ASSERT(res);
	md_range_t range = {0};
	if (res->atom_offset && res_idx < res->count) {
		range.beg = res->atom_offset[res_idx];
		range.end = res->atom_offset[res_idx + 1];
	}
	return range;
}

static inline md_residue_idx_t md_residue_find_by_atom_idx(const md_residue_data_t* res, size_t atom_idx) {
    ASSERT(res);

    md_residue_idx_t res_idx = -1;
    if (res->atom_offset) {
        for (size_t i = 0; i < res->count; ++i) {
            size_t res_beg = res->atom_offset[i];
            size_t res_end = res->atom_offset[i + 1];
            if (res_beg <= atom_idx && atom_idx < res_end) {
                res_idx = (md_residue_idx_t)i;
                break;
            }
            if (res_beg > atom_idx) {
                break;
            }
        }
    }

    return res_idx;
}

static inline size_t md_residue_atom_count(const md_residue_data_t* res, size_t res_idx) {
    ASSERT(res);
    size_t count = 0;

    if (res->atom_offset && res_idx < res->count) {
        count = res->atom_offset[res_idx + 1] - res->atom_offset[res_idx];
    }
    return count;
}


static inline md_range_t md_chain_residue_range(const md_chain_data_t* chain, size_t chain_idx) {
    ASSERT(chain);

    md_range_t range = {0};
    if (chain->res_range && chain_idx < chain->count) {
        range = chain->res_range[chain_idx];
    }
    return range;
}

static inline md_chain_idx_t md_chain_find_by_residue_idx(const md_chain_data_t* chain, size_t res_idx) {
    ASSERT(chain);

    md_chain_idx_t chain_idx = -1;
    if (chain->res_range) {
        int ri = (int)res_idx;
        for (size_t i = 0; i < chain->count; ++i) {
            md_range_t range = chain->res_range[i];
            if (range.beg <= ri && ri < range.end) {
                chain_idx = (md_chain_idx_t)i;
                break;
            }
            if (range.beg > ri) {
                break;
            }
        }
    }
    return chain_idx;
}

static inline md_chain_idx_t md_chain_find_by_atom_idx(const md_chain_data_t* chain, size_t atom_idx) {
    ASSERT(chain);

    md_chain_idx_t chain_idx = -1;
    if (chain->atom_range) {
        int ai = (int)atom_idx;
        for (size_t i = 0; i < chain->count; ++i) {
            md_range_t range = chain->atom_range[i];
            if (range.beg <= ai && ai < range.end) {
                chain_idx = (md_chain_idx_t)i;
                break;
            }
            if (range.beg > ai) {
                break;
            }
        }
    }
    return chain_idx;
}

static inline size_t md_chain_residue_count(const md_chain_data_t* chain, size_t chain_idx) {
    ASSERT(chain);

    size_t count = 0;
    if (chain->res_range && chain_idx < chain->count) {
        md_range_t range = chain->res_range[chain_idx];
        count = range.end - range.beg;
    }
    return count;
}

static inline md_range_t md_chain_atom_range(const md_chain_data_t* chain, size_t chain_idx) {
	ASSERT(chain);

    md_range_t range = {0};
    if (chain->atom_range && chain_idx < chain->count) {
        range = chain->atom_range[chain_idx];
    }
    return range;
}

static inline size_t md_chain_atom_count(const md_chain_data_t* chain, size_t chain_idx) {
    size_t count = 0;
    if (chain->atom_range && chain_idx < chain->count) {
        md_range_t range = chain->atom_range[chain_idx];
        count = range.end - range.beg;
    }
    return count;
}

static inline str_t md_chain_id(const md_chain_data_t* chain, size_t chain_idx) {
    ASSERT(chain);
    str_t id = STR_LIT("");
    if (chain->id && chain_idx < chain->count) {
        id = LBL_TO_STR(chain->id[chain_idx]);
    }
    return id;
}

static inline size_t md_structure_count(const md_index_data_t* structure) {
    ASSERT(structure);
    return md_index_data_num_ranges(*structure);
}

static inline size_t md_structure_atom_count(const md_index_data_t* structure, size_t struct_idx) {
    ASSERT(structure);
    return md_index_range_size(*structure, struct_idx);
}

static inline const md_atom_idx_t* md_structure_atom_indices(const md_index_data_t* structure, size_t struct_idx) {
    ASSERT(structure);
    return (const md_atom_idx_t*)md_index_range_beg(*structure, struct_idx);
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

#define MD_BOND_FLAG_MASK  0xF0
#define MD_BOND_ORDER_MASK 0x0F

static inline size_t md_bond_conn_count(md_bond_data_t bond_data, size_t atom_idx) {
    return bond_data.conn.offset[atom_idx + 1] - bond_data.conn.offset[atom_idx];
}

static inline void md_bond_order_set(md_bond_data_t* bond_data, md_bond_idx_t bond_idx, uint32_t order) {
    uint32_t o = bond_data->order[bond_idx];
    bond_data->order[bond_idx] = (o & MD_BOND_FLAG_MASK) | (order & MD_BOND_ORDER_MASK);
}

static inline md_atom_idx_t md_bond_conn_atom_idx(md_bond_data_t bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    return bond_data.conn.atom_idx[atom_conn_idx + idx];
}

static inline md_bond_idx_t md_bond_conn_bond_idx(md_bond_data_t bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    return bond_data.conn.bond_idx[atom_conn_idx + idx];
}

static inline bool md_bond_iter_has_next(md_bond_iter_t it) {
    return it.i < it.end_idx;
}

static inline void md_bond_iter_next(md_bond_iter_t* it) {
    ++it->i;
}

static inline md_atom_idx_t md_bond_iter_atom_index(md_bond_iter_t it) {
	return it.data->conn.atom_idx[it.i];
}

static inline md_atom_idx_t md_bond_iter_bond_index(md_bond_iter_t it) {
    return it.data->conn.bond_idx[it.i];
}

static inline uint32_t md_bond_iter_bond_order(md_bond_iter_t it) {
    return it.data->order[it.data->conn.bond_idx[it.i]] & MD_BOND_ORDER_MASK;
}

static inline uint32_t md_bond_iter_bond_flags(md_bond_iter_t it) {
    return it.data->order[it.data->conn.bond_idx[it.i]] & MD_BOND_FLAG_MASK;
}

static inline void md_bond_conn_clear(md_conn_data_t* conn_data) {
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
    md_array_shrink(bond_data->order, 0);
    
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
    while (md_bond_iter_has_next(it) && !found) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(it);
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

/*

The molecule loader interface is just a convenience interface for abstracing the functionality of initializing molecule data

The reason for providing a distinct function for initializing from file is that some molecule files can
also contain their trajectories, such as PDB files. In such case, the whole file would have to be read and passed, but for
molecule data only the first part of the file is used.

*/

typedef struct md_molecule_loader_i {
    bool (*init_from_str) (md_molecule_t* mol, str_t string,   const void* arg, struct md_allocator_i* alloc);
    bool (*init_from_file)(md_molecule_t* mol, str_t filename, const void* arg, struct md_allocator_i* alloc);
} md_molecule_loader_i;

// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
void md_molecule_free(md_molecule_t* mol, struct md_allocator_i* alloc);
void md_molecule_copy(md_molecule_t* dst_mol, const md_molecule_t* src_mol, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
