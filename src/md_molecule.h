#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>

typedef struct md_atom_data_t {
    size_t count;
    // Coordinates
    float* x;
    float* y;
    float* z;
    // Velocities
    float* vx;
    float* vy;
    float* vz;
    // Misc
    float* radius;
    float* mass;
    md_valence_t* valence;
    md_element_t* element;
    md_label_t* type;

    md_residue_id_t* resid;
    md_label_t* resname;
    md_label_t* chainid;

    md_residue_idx_t* res_idx;
    md_chain_idx_t* chain_idx;

    // Stores offset into connectivity data for the atom
    // Length is implicit in the next atom's offset, which means this has length count + 1 if present
    uint32_t* conn_offset;

    md_flags_t* flags;
} md_atom_data_t;

// These are a bit superflous, but just replicate whatever is present in the atom data
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
    uint32_t* res_offset;
    uint32_t* atom_offset;
} md_chain_data_t;

typedef struct md_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    size_t range_count;
    md_range_t* range;

    // These fields share the same length 'count'
    size_t count;
    md_backbone_atoms_t* atoms;
    md_backbone_angles_t* angle;
    md_secondary_structure_t* secondary_structure;
    md_ramachandran_type_t* ramachandran_type;
    md_residue_idx_t* residue_idx;            // Index to the residue which contains the backbone
} md_backbone_data_t;

// This represents symmetries which are instanced, commonly found
// in PDB data. It is up to the renderer to properly render this instanced data.
typedef struct md_instance_data_t {
    size_t count;
    md_range_t* atom_range;
    md_label_t* label;
    mat4_t* transform;
} md_instance_data_t;

// Bond centric representation of bonds
typedef struct md_bond_data_t {
    size_t count;
    md_bond_pair_t* pairs;
    md_order_t*     order;
    md_flags_t*     flags;
} md_bond_data_t;

typedef struct md_hbond_data_t {
    int8_t charge;
    int8_t impl_h;
    int8_t total_h;
    int8_t ideal_geom;
} md_hbond_data_t;

// Atom centric representation of bonds
typedef struct md_conn_data_t {
    size_t count;
    md_atom_idx_t*  index;
    uint8_t*        order;
    md_flags_t*     flags;
} md_conn_data_t;

typedef struct md_conn_iter_t {
	const md_conn_data_t* data;
	int64_t i;
    int64_t beg_idx;
    int64_t end_idx;
} md_conn_iter_t;

typedef struct md_molecule_t {
    md_unit_cell_t          unit_cell;
    md_atom_data_t          atom;
    md_residue_data_t       residue;
    md_chain_data_t         chain;
    md_backbone_data_t      backbone;
    
    md_bond_data_t          bond;       // Bond centric data of (strong persistent bonds)
    md_conn_data_t          conn;       // Atom centric representation of the bond data
    
    // @TODO: move this into some containing structure
    //md_array(md_hbond_data_t)  hbond_data;     

    md_index_data_t         rings;              // Ring structures formed by persistent bonds
    md_index_data_t         structures;         // Isolated structures connected by persistent bonds

    md_instance_data_t      instance;           // Symmetry instances of the molecule (duplications of ranges with new transforms)
    
    // @NOTE(Robin): This should probably move elsewhere.
    // Hydrogen bonds are of interest to be evaluated and analyzed over the trajectory
    // So this should be computed over all frames in the preload
    md_array(md_bond_pair_t)     hydrogen_bonds;
} md_molecule_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline md_range_t md_residue_atom_range(md_residue_data_t res, int64_t res_idx) {
	md_range_t range = {0};
    size_t i = (size_t)res_idx; // Cast to unsigned to only check for positive
	if (res.atom_offset && i < res.count) {
		range.beg = res.atom_offset[i];
		range.end = res.atom_offset[i + 1];
	}
	return range;
}

static inline size_t md_residue_atom_count(md_residue_data_t res, int64_t res_idx) {
    size_t count = 0;
    size_t i = (size_t)res_idx;
    if (res.atom_offset && i < res.count) {
        count = res.atom_offset[i + 1] - res.atom_offset[i];
    }
    return count;
}


static inline md_range_t md_chain_residue_range(md_chain_data_t chain, int64_t chain_idx) {
    md_range_t range = {0};
    size_t i = (size_t)chain_idx; // Cast to unsigned to only check for positive
    if (chain.res_offset && i < chain.count) {
        range.beg = chain.res_offset[i];
        range.end = chain.res_offset[i + 1];
    }
    return range;
}

static inline size_t md_chain_residue_count(md_chain_data_t chain, int64_t chain_idx) {
    size_t count = 0;
    size_t i = (size_t)chain_idx; // Cast to unsigned to only check for positive
    if (chain.res_offset && i < chain.count) {
        count = chain.res_offset[i + 1] - chain.res_offset[i];
    }
    return count;
}

static inline md_range_t md_chain_atom_range(md_chain_data_t chain, int64_t chain_idx) {
    md_range_t range = {0};
    size_t i = (size_t)chain_idx; // Cast to unsigned to only check for positive
    if (chain.atom_offset && i < chain.count) {
        range.beg = chain.atom_offset[i];
        range.end = chain.atom_offset[i + 1];
    }
    return range;
}

static inline size_t md_chain_atom_count(md_chain_data_t chain, int64_t chain_idx) {
    size_t count = 0;
    size_t i = (size_t)chain_idx; // Cast to unsigned to only check for positive
    if (chain.atom_offset && i < chain.count) {
        count = chain.atom_offset[i + 1] - chain.atom_offset[i];
    }
    return count;
}

static inline md_conn_iter_t md_conn_iter(const md_molecule_t* mol, int64_t atom_idx) {
    md_conn_iter_t it = {0};
    if (mol->atom.conn_offset) {
        it.data = &mol->conn;
		it.beg_idx = mol->atom.conn_offset[atom_idx];
		it.end_idx = mol->atom.conn_offset[atom_idx + 1];
        it.i = it.beg_idx;
    }
    return it;
}

static inline bool md_conn_iter_has_next(const md_conn_iter_t* it) {
    return it->i < it->end_idx;
}

static inline void md_conn_iter_next(md_conn_iter_t* it) {
    ++it->i;
}

static inline md_atom_idx_t md_conn_iter_index(const md_conn_iter_t* it) {
	return it->data->index[it->i];
}

static inline md_order_t md_conn_iter_order(const md_conn_iter_t* it) {
	return it->data->order[it->i];
}

static inline md_flags_t md_conn_iter_flags(const md_conn_iter_t* it) {
	return it->data->flags[it->i];
}

/*

The molecule loader api is just a convenience API for abstracing the functionality of initializing molecule data

The reason for providing a distinct function for initializing from file is that some molecule files can
also contain their trajectories, such as PDB files. In such case, the whole file would have to be read and passed, but for
molecule data only the first part of the file is used.

*/

typedef struct md_molecule_loader_i {
    bool (*init_from_str) (md_molecule_t* mol, str_t string,   struct md_allocator_i* alloc);
    bool (*init_from_file)(md_molecule_t* mol, str_t filename, struct md_allocator_i* alloc);
} md_molecule_loader_i;

// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
void md_molecule_free(md_molecule_t* mol, struct md_allocator_i* alloc);

void md_molecule_copy(md_molecule_t* dst_mol, const md_molecule_t* src_mol, struct md_allocator_i* alloc);

// Convenience functions to extract vec3_soa streams from molecule
static inline md_vec3_soa_t md_molecule_soa_coord(md_atom_data_t* atom_data) {
    md_vec3_soa_t soa = {atom_data->x, atom_data->y, atom_data->z, sizeof(float)};
    return soa;
}

static inline md_vec3_soa_t md_molecule_soa_vel(md_atom_data_t* atom_data) {
    md_vec3_soa_t soa = {atom_data->vx, atom_data->vy, atom_data->vz, sizeof(float)};
    return soa;
}

#ifdef __cplusplus
}
#endif
