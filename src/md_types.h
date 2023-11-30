#pragma once


#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>

#include <stdint.h>

// We are sneaky, we encode the secondary structure as a uint8x4 unorm where the the components encode the fraction of each secondary structure type
// This allows us later on to interpolate between them
enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN  = 0,
    MD_SECONDARY_STRUCTURE_COIL     = 0x000000FF,
    MD_SECONDARY_STRUCTURE_HELIX    = 0x0000FF00,
    MD_SECONDARY_STRUCTURE_SHEET    = 0x00FF0000
};

enum {
    MD_RAMACHANDRAN_TYPE_UNKNOWN    = 0,
    MD_RAMACHANDRAN_TYPE_GENERAL    = 1,
    MD_RAMACHANDRAN_TYPE_GLYCINE    = 2,
    MD_RAMACHANDRAN_TYPE_PROLINE    = 3,
    MD_RAMACHANDRAN_TYPE_PREPROL    = 4,
};

// Cell information
enum {
    MD_CELL_NONE                    = 0,
    MD_CELL_ORTHOGONAL              = 1,
    MD_CELL_TRICLINIC               = 2,
    MD_CELL_X                       = 4,
    MD_CELL_Y                       = 8,
    MD_CELL_Z                       = 16,
};

// These flags are not specific to any distinct subtype, but can appear in both atoms, residues, bonds and whatnot.
// Where ever they make sense, they can appear. This makes it easy to propagate the flags upwards and downwards
enum {
    MD_FLAG_RES_BEG 		    = 1,
    MD_FLAG_RES_END 		    = 2,
    MD_FLAG_CHAIN_BEG 		    = 4,
    MD_FLAG_CHAIN_END 		    = 8,
    MD_FLAG_HETATM              = 16,
    MD_FLAG_AMINO_ACID		    = 32,
    MD_FLAG_NUCLEOBASE          = 64,
    MD_FLAG_NUCLEOTIDE	        = 128,
    MD_FLAG_WATER			    = 256,
    MD_FLAG_ION			        = 512,

    MD_FLAG_AROMATIC            = 1024,
    MD_FLAG_INTER_BOND          = 2048,
};

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_residue_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_residue_id_t;
typedef int32_t     md_chain_idx_t;
typedef int32_t     md_bond_idx_t;
typedef uint32_t    md_secondary_structure_t;
typedef uint32_t    md_flags_t;
typedef uint8_t     md_element_t;
typedef uint8_t     md_ramachandran_type_t;
typedef uint8_t     md_valence_t;
typedef uint8_t     md_order_t;

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

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_range_t {
    int32_t beg;
    int32_t end;
} md_range_t;

typedef struct md_unit_cell_t {
    mat3_t basis;
    mat3_t inv_basis;
    uint32_t flags;
} md_unit_cell_t;

/*
// Single bond between two entities represented by atom indices
typedef struct md_bond_t {
    md_atom_idx_t idx[2];
    uint16_t order;
    uint16_t flags;
} md_bond_t;
*/

typedef struct md_bond_pair_t {
    md_atom_idx_t idx[2];
} md_bond_pair_t;

typedef struct md_backbone_atoms_t {
    md_atom_idx_t n;
    md_atom_idx_t ca;
    md_atom_idx_t c;
    md_atom_idx_t o;
} md_backbone_atoms_t;

// Backbone angles
// φ (phi) is the angle in the chain C' − N − Cα − C'
// ψ (psi) is the angle in the chain N − Cα − C' − N
// https://en.wikipedia.org/wiki/Dihedral_angle#/media/File:Protein_backbone_PhiPsiOmega_drawing.svg)
typedef struct md_backbone_angles_t {
    float phi;
    float psi;
} md_backbone_angles_t;

// Miniature string buffer with explicit length
// It can store up to 6 characters + null terminator which makes it compatible as a C-string
// This is to store the atom symbol and other short strings common in molecular data
// The motiviation is that it saves an indirection
typedef struct md_label_t {
    char    buf[7];
    uint8_t len;

#ifdef __cplusplus
    constexpr operator str_t() { return {buf, len}; }
    constexpr operator const char*() { return buf; }
#endif
} md_label_t;

// Container structure for ranges of indices.
// It stores a set of indices for every stored entity. (Think of it as an Array of Arrays of integers)
// The indices are packed into a single array and we store explicit offsets to represent ranges within this index array.
// It is used to represent connected structures, rings and atom connectivity etc.

typedef struct md_index_data_t {
    md_array(uint32_t) offsets;
    md_array(int32_t)  indices;
} md_index_data_t;

// OPERATIONS ON THE TYPES

// macro concatenate trick to assert that the input is a valid compile time C-string
#define MAKE_LABEL(cstr) {cstr"", sizeof(cstr)-1}

#ifdef __cplusplus
#define LBL_TO_STR(lbl) {lbl.buf, lbl.len}
inline bool operator==(const md_label_t& a, const md_label_t& b) {
    return MEMCMP(&a, &b, sizeof(md_label_t)) == 0;
}
inline bool operator!=(const md_label_t& a, const md_label_t& b) {
    return MEMCMP(&a, &b, sizeof(md_label_t)) != 0;
}
#else
#define LBL_TO_STR(lbl) (str_t){lbl.buf, lbl.len}
#endif // __cplusplus

static inline bool label_empty(md_label_t lbl) {
    uint64_t ref = 0;
    return MEMCMP(&lbl, &ref, sizeof(md_label_t)) == 0;
}

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
    if (data->offsets) md_array_free(data->offsets,  alloc);
    if (data->indices) md_array_free(data->indices, alloc);
    MEMSET(data, 0, sizeof(md_index_data_t));
}

static inline int64_t md_index_data_push_arr (md_index_data_t* data, const int32_t* index_data, int64_t index_count, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    ASSERT(index_count >= 0);

    if (md_array_size(data->offsets) == 0) {
        md_array_push(data->offsets, 0, alloc);
    }

    int64_t offset = *md_array_last(data->offsets);
    if (index_count > 0) {
        ASSERT(index_data);
        md_array_push_array(data->indices, index_data, index_count, alloc);
        md_array_push(data->offsets, (int32_t)md_array_size(data->indices), alloc);
    }

    return offset;
}

static inline void md_index_data_data_clear(md_index_data_t* data) {
    ASSERT(data);
    md_array_shrink(data->offsets,  0);
    md_array_shrink(data->indices, 0);
}

static inline int64_t md_index_data_count(md_index_data_t data) {
    return MAX(md_array_size(data.offsets) - 1, 0);
}

// Access to individual substructures
static inline int32_t* md_index_range_beg(md_index_data_t data, int64_t idx) {
    ASSERT(idx >= 0 && idx < md_array_size(data.offsets) - 1);
    return data.indices + data.offsets[idx];
}

static inline int32_t* md_index_range_end(md_index_data_t data, int64_t idx) {
    ASSERT(idx >= 0 && idx < md_array_size(data.offsets) - 1);
    return data.indices + data.offsets[idx+1];
}

static inline int64_t md_index_range_size(md_index_data_t data, int64_t idx) {
    ASSERT(idx >= 0 && idx < md_array_size(data.offsets) - 1);
    return data.offsets[idx+1] - data.offsets[idx];
}
