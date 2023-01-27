#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>

// Forward declarations
struct md_allocator_i;

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_residue_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_residue_id_t;
typedef int32_t     md_chain_idx_t;
typedef int32_t     md_molecule_idx_t;
typedef uint32_t    md_secondary_structure_t;
typedef uint32_t    md_flags_t;
typedef uint8_t     md_element_t;
typedef uint8_t     md_ramachandran_type_t;
typedef uint8_t     md_valence_t;
typedef mat3_t      md_coordinate_frame_t;

// We are sneaky, we encode the secondary structure as a uint8x4 unorm where the the components encode the fraction of each secondary structure type
// This allows us later on to interpolate between them
enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MD_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MD_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000
};

enum {
    MD_RAMACHANDRAN_TYPE_UNKNOWN    = 0,
    MD_RAMACHANDRAN_TYPE_GENERAL    = 1,
    MD_RAMACHANDRAN_TYPE_GLYCINE    = 2,
    MD_RAMACHANDRAN_TYPE_PROLINE    = 3,
    MD_RAMACHANDRAN_TYPE_PREPROL    = 4,
};

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_range_t {
    int32_t beg;
    int32_t end;
} md_range_t;

// Single bond between two entities represented by the indices
typedef struct md_bond_t {
    md_atom_idx_t idx[2];
} md_bond_t;

typedef struct md_backbone_atoms_t {
    md_atom_idx_t n;
    md_atom_idx_t ca;
    md_atom_idx_t c;
    md_atom_idx_t o;
} md_backbone_atoms_t;

typedef struct md_backbone_angles_t {
    float phi;
    float psi;
} md_backbone_angles_t;

// Miniature string buffer with length
typedef struct md_label_t {
    char    buf[7];
    uint8_t len;

#ifdef __cplusplus
    constexpr operator str_t() { return {buf, len}; }
    constexpr operator const char*() { return buf; }
#endif
} md_label_t;

// Structure Of Array layout version of vec3_t
// This is to simplify the interfaces a bit when dealing with multiple coordinate streams
typedef struct md_vec3_soa_t {
    float* x;
    float* y;
    float* z;
    // The stride signifies the byte stride between each entry of the pointers.
    // If the arrays are packed and not part of a larger struct, this value would be sizeof(float) == 4
    uint64_t stride;
} md_vec3_soa_t;

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

// Container structure for ranges of indices.
// It stores a set of indices for every stored entity. (Think of it as an Array of Arrays of integers)
// The indices are packed into a single array and we store explicit ranges referencing spans within this index array.
// It is used to represent connected structures, rings and atom connectivity etc.
typedef struct md_index_data_t {
    md_array(md_range_t) ranges;
    md_array(int32_t)    indices;
} md_index_data_t;

// BOND DATA
typedef struct md_bond_data_t {
    int64_t count;
    md_array(md_bond_t) bond;
    md_index_data_t connectivity;   // Connectivity of atoms
    md_index_data_t structures;     // Isolated structures formed by bonds
    md_index_data_t rings;          // Rings structures formed by bonds
} md_bond_data_t;

// This represents symmetries which are instanced, commonly found
// in PDB data. It is up to the renderer to properly render this instanced data.
typedef struct md_instance_data_t {
    int64_t count;
    md_range_t* atom_range;
    md_label_t* label;
    mat4_t* transform;
} md_instance_data_t;

typedef struct md_molecule_t {
    md_coordinate_frame_t   coord_frame;
    md_atom_data_t          atom;
    md_residue_data_t       residue;
    md_chain_data_t         chain;
    md_backbone_data_t      backbone;
    md_bond_data_t          covalent;
    md_bond_data_t          hydrogen;
    md_instance_data_t      instance;
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

typedef struct md_molecule_api {
    bool (*init_from_str) (md_molecule_t* mol, str_t string,   struct md_allocator_i* alloc);
    bool (*init_from_file)(md_molecule_t* mol, str_t filename, struct md_allocator_i* alloc);
} md_molecule_api;

// @NOTE(Robin): This is just to be thorough,
// I would recommend using an explicit arena allocator for the molecule and just clearing that in one go instead of calling this.
void md_molecule_free(md_molecule_t* mol, struct md_allocator_i* alloc);

// Append the internal data of one molecule to another
void md_molecule_append(md_molecule_t* dst_mol, const md_molecule_t* src_mol, struct md_allocator_i* alloc);

void md_molecule_copy(md_molecule_t* dst_mol, const md_molecule_t* src_mol, struct md_allocator_i* alloc);

// macro concatenate trick to assert that the input is a valid compile time C-string
#define MAKE_LABEL(cstr) {cstr"", sizeof(cstr)-1}

#ifdef __cplusplus
#define LBL_TO_STR(lbl) {lbl.buf, lbl.len}
#else
#define LBL_TO_STR(lbl) (str_t){lbl.buf, lbl.len}
#endif // __cplusplus

static inline md_label_t make_label(str_t str) {
    md_label_t lbl = {0};
    if (str.ptr) {
        const int64_t len = MIN(str.len, (int64_t)sizeof(lbl.buf) - 1);
        for (int64_t i = 0; i < len; ++i) {
            lbl.buf[i] = str.ptr[i];
        }
        lbl.len = (uint8_t)len;
    }
    return lbl;
}

// Convenience functions to extract vec3_soa streams from molecule
static inline md_vec3_soa_t md_molecule_soa_coord(md_atom_data_t* atom_data) {
    md_vec3_soa_t soa = {atom_data->x, atom_data->y, atom_data->z, sizeof(float)};
    return soa;
}

static inline md_vec3_soa_t md_molecule_soa_vel(md_atom_data_t* atom_data) {
    md_vec3_soa_t soa = {atom_data->vx, atom_data->vy, atom_data->vz, sizeof(float)};
    return soa;
}

static inline md_vec3_soa_t md_vec3_soa_from_vec3(vec3_t* ptr) {
    md_vec3_soa_t soa = {&ptr->x, &ptr->y, &ptr->z, sizeof(vec3_t)};
    return soa;
}

static inline vec3_t md_vec3_soa_get(md_vec3_soa_t soa, int64_t idx) {
    ASSERT(soa.stride != 0);
    vec3_t v = {*(float*)((char*)soa.x + idx * soa.stride), *(float*)((char*)soa.y + idx * soa.stride), *(float*)((char*)soa.z + idx * soa.stride)};
    return v;
}

// Access to substructure data
static inline void md_index_data_free (md_index_data_t* data, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    md_array_free(data->ranges,  alloc);
    md_array_free(data->indices, alloc);
}

static inline int64_t md_index_data_push (md_index_data_t* data, int32_t* index_data, int64_t index_count, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    ASSERT(index_count >= 0);

    const int64_t range_idx = md_array_size(data->ranges);
    const md_range_t range = {(int)md_array_size(data->indices), (int)md_array_size(data->indices) + (int)index_count};
    md_array_push(data->ranges, range, alloc);
    if (index_count > 0) {
        ASSERT(index_data);
        md_array_push_array(data->indices, index_data, index_count, alloc);
    }

    return range_idx;
}

static inline void md_index_data_data_clear(md_index_data_t* data) {
    ASSERT(data);
    md_array_shrink(data->ranges,  0);
    md_array_shrink(data->indices, 0);
}

FORCE_INLINE int64_t md_index_data_count(const md_index_data_t* data) {
    ASSERT(data);
    return md_array_size(data->ranges);
}

// Access to individual substructures
static inline int32_t* md_index_range_beg(const md_index_data_t* data, int64_t idx) {
    ASSERT(data);
    ASSERT(idx >= 0 && idx < md_array_size(data->ranges));
    return data->indices + data->ranges[idx].beg;
}

static inline int32_t* md_index_range_end(const md_index_data_t* data, int64_t idx) {
    ASSERT(data);
    ASSERT(idx >= 0 && idx < md_array_size(data->ranges));
    return data->indices + data->ranges[idx].end;
}

static inline int64_t md_index_range_size(const md_index_data_t* data, int64_t idx) {
    ASSERT(data);
    ASSERT(idx >= 0 && idx < md_array_size(data->ranges));
    return data->ranges[idx].end - data->ranges[idx].beg;
}

static inline void md_bond_data_clear(md_bond_data_t* data) {
    ASSERT(data);
    md_array_shrink(data->bond, 0);
    md_index_data_data_clear(&data->connectivity);
}

static inline int64_t md_bond_count(const md_bond_data_t* bond_data) {
    ASSERT(bond_data);
    return md_array_size(bond_data->bond);
}

#ifdef __cplusplus
}
#endif
