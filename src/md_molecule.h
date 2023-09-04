#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>

typedef struct md_atom_data_t {
    int64_t count;
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
    md_label_t* name;
    md_flags_t* flags;                          // Auxillary bit buffer for flagging individual atoms
    md_residue_idx_t* residue_idx;
    md_chain_idx_t* chain_idx;
} md_atom_data_t;

typedef struct md_residue_data_t {
    int64_t count;
    md_label_t* name;
    md_residue_id_t* id;
    md_range_t* atom_range;
} md_residue_data_t;

typedef struct md_chain_data_t {
    int64_t count;
    md_label_t* id;
    md_range_t* residue_range;
    md_range_t* atom_range;
} md_chain_data_t;

typedef struct md_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    int64_t range_count;
    md_range_t* range;

    // These fields share the same length 'count'
    int64_t count;
    md_backbone_atoms_t* atoms;
    md_backbone_angles_t* angle;
    md_secondary_structure_t* secondary_structure;
    md_ramachandran_type_t* ramachandran_type;
    md_residue_idx_t* residue_idx;            // Index to the residue which contains the backbone
} md_backbone_data_t;

// This represents symmetries which are instanced, commonly found
// in PDB data. It is up to the renderer to properly render this instanced data.
typedef struct md_instance_data_t {
    int64_t count;
    md_range_t* atom_range;
    md_label_t* label;
    mat4_t* transform;
} md_instance_data_t;

/*
typedef struct md_bond_data_t {
    // Should have the same length
    md_array(md_bond_t)     bonds;
    md_array(md_order_t)    order;

    md_array(md_range_t)    ranges; // should have length atom.count
    md_array(md_bond_idx_t) atom_bonds;
    md_array(md_atom_idx_t) atom_connections;
} md_bond_data_t;
*/

typedef struct md_hbond_data_t {
    int8_t charge;
    int8_t impl_h;
    int8_t total_h;
    int8_t ideal_geom;
} md_hbond_data_t;

typedef struct md_molecule_t {
    md_unit_cell_t          unit_cell;
    md_atom_data_t          atom;
    md_residue_data_t       residue;
    md_chain_data_t         chain;
    md_backbone_data_t      backbone;
    
    md_array(md_bond_t)     bonds;              // Strong bonds which persist throughout the simulation
    
    // @TODO: move this into some containing structure
    md_array(md_hbond_data_t)  hbond_data;     

    md_index_data_t         connectivity;       // Connectivity graph of the atoms within the molecule
    md_index_data_t         rings;              // Ring structures formed by persistent bonds
    md_index_data_t         structures;         // Isolated structures connected by persistent bonds

    md_instance_data_t      instance;           // Symmetry instances of the molecule (duplications of ranges with new transforms)
    
    // @NOTE(Robin): This should probably move elsewhere.
    // Hydrogen bonds are of interest to be evaluated and analyzed over the trajectory
    // So this should be computed over all frames in the preload
    md_array(md_bond_t)     hydrogen_bonds;
} md_molecule_t;

/*

The molecule api is just a convenience API for abstracing the functionality of initializing molecule data

The reason for providing a distinct function for initializing from file is that some molecule files can
also contain their trajectories, such as PDB files. In such case, the whole file would have to be read and passed, but for
molecule data only the first part of the file is used.

*/

#ifdef __cplusplus
extern "C" {
#endif

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
