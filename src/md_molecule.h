#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>

typedef struct md_atom_type_data_t {
    size_t count;

    md_label_t*   name;
    md_element_t* element;
    float*        mass;
    float*        radius;
} md_atom_type_data_t;

typedef struct md_atom_data_t {
    size_t count;

    // Coordinates
    float* x;
    float* y;
    float* z;

    // Atom Type Specific (@TODO: Compress this into a smaller set)
    float* radius;
    float* mass;
    md_element_t* element;
    md_label_t* type;

    // @TODO: Move these into res and chain data structs
    md_residue_id_t* resid;
    md_label_t* resname;
    md_label_t* chainid;

    md_residue_idx_t* res_idx;
    md_chain_idx_t* chain_idx;

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
    md_atom_pair_t* pairs;
    md_order_t*     order;
    md_conn_data_t  conn;   // Connectivity
} md_bond_data_t;

typedef struct md_bond_iter_t {
    const md_bond_data_t* data;
	uint32_t i;
    uint32_t end_idx;
} md_bond_iter_t;

typedef struct md_hydrogen_bond_donor_t {
    md_atom_idx_t d_idx; // donor
    md_atom_idx_t h_idx; // hydrogen
} md_hydrogen_bond_donor_t;

typedef struct md_hydrogen_bond_acceptor_t {
    md_atom_idx_t idx;
    int num_of_lone_pairs;
} md_hydrogen_bond_acceptor_t;

typedef struct md_hydrogen_bond_pair_t {
	uint32_t donor_idx;    // Index into donors array
	uint32_t acceptor_idx; // Index into acceptors array
} md_hydrogen_bond_pair_t;

typedef struct md_hydrogen_bond_data_t {
    size_t num_acceptors;
    md_hydrogen_bond_acceptor_t* acceptors;
   
    size_t num_donors;
    md_hydrogen_bond_donor_t* donors;

    size_t num_bonds;
    md_hydrogen_bond_pair_t* bonds; // index[0] = donor atom idx, index[1] = acceptor atom idx
} md_hydrogen_bond_data_t;

typedef struct md_molecule_t {
    md_unitcell_t               unitcell;
    md_atom_data_t              atom;
    md_residue_data_t           residue;
    md_chain_data_t             chain;
    md_protein_backbone_data_t  protein_backbone;
    md_nucleic_backbone_data_t  nucleic_backbone;
    
    md_bond_data_t              bond;               // Covalent bonds
    md_hydrogen_bond_data_t     hydrogen_bond;      // Hydrogen bonds

    md_index_data_t             ring;               // Ring structures formed by persistent bonds
    md_index_data_t             structure;          // Isolated structures connected by persistent bonds

    md_instance_data_t          instance;           // Instances of the molecule (duplications of ranges with new transforms)
} md_molecule_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline vec3_t md_atom_coord(md_atom_data_t atom_data, size_t atom_idx) {
    return vec3_set(atom_data.x[atom_idx], atom_data.y[atom_idx], atom_data.z[atom_idx]);
}

static inline md_range_t md_residue_atom_range(md_residue_data_t res, size_t res_idx) {
	md_range_t range = {0};
	if (res.atom_offset && res_idx < res.count) {
		range.beg = res.atom_offset[res_idx];
		range.end = res.atom_offset[res_idx + 1];
	}
	return range;
}

static inline size_t md_residue_atom_count(md_residue_data_t res, size_t res_idx) {
    size_t count = 0;
    if (res.atom_offset && res_idx < res.count) {
        count = res.atom_offset[res_idx + 1] - res.atom_offset[res_idx];
    }
    return count;
}


static inline md_range_t md_chain_residue_range(md_chain_data_t chain, size_t chain_idx) {
    md_range_t range = {0};
    if (chain.res_range && chain_idx < chain.count) {
        range = chain.res_range[chain_idx];
    }
    return range;
}

static inline size_t md_chain_residue_count(md_chain_data_t chain, size_t chain_idx) {
    size_t count = 0;
    if (chain.res_range && chain_idx < chain.count) {
        md_range_t range = chain.res_range[chain_idx];
        count = range.end - range.beg;
    }
    return count;
}

static inline md_range_t md_chain_atom_range(md_chain_data_t chain, size_t chain_idx) {
    md_range_t range = {0};
    if (chain.atom_range && chain_idx < chain.count) {
        range = chain.atom_range[chain_idx];
    }
    return range;
}

static inline size_t md_chain_atom_count(md_chain_data_t chain, size_t chain_idx) {
    size_t count = 0;
    if (chain.atom_range && chain_idx < chain.count) {
        md_range_t range = chain.atom_range[chain_idx];
        count = range.end - range.beg;
    }
    return count;
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

static inline md_atom_idx_t md_hydrogen_bond_donor_atom_idx(const md_hydrogen_bond_data_t* hbond_data, size_t donor_idx) {
	ASSERT(hbond_data);
    if (donor_idx < hbond_data->num_donors) {
        return hbond_data->donors[donor_idx].d_idx;
	}
    return -1;
}

static inline md_atom_idx_t md_hydrogen_bond_donor_hydrogen_atom_idx(const md_hydrogen_bond_data_t* hbond_data, size_t donor_idx) {
    ASSERT(hbond_data);
    if (donor_idx < hbond_data->num_donors) {
        return hbond_data->donors[donor_idx].h_idx;
    }
    return -1;
}

static inline md_atom_idx_t md_hydrogen_bond_acceptor_atom_idx(const md_hydrogen_bond_data_t* hbond_data, size_t acceptor_idx) {
    ASSERT(hbond_data);
    if (acceptor_idx < hbond_data->num_acceptors) {
        return hbond_data->acceptors[acceptor_idx].idx;
    }
    return -1;
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
