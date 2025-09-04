#pragma once


#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>

#include <stdint.h>

#if DEBUG
#include <core/md_log.h>
#endif

enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL,
    MD_SECONDARY_STRUCTURE_TURN,
    MD_SECONDARY_STRUCTURE_BEND,
    MD_SECONDARY_STRUCTURE_HELIX_310,
    MD_SECONDARY_STRUCTURE_HELIX_ALPHA,
    MD_SECONDARY_STRUCTURE_HELIX_PI,
    MD_SECONDARY_STRUCTURE_BETA_SHEET,
    MD_SECONDARY_STRUCTURE_BETA_BRIDGE,
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
    MD_UNITCELL_NONE          = 0,
    MD_UNITCELL_ORTHO         = 1,
    MD_UNITCELL_TRICLINIC     = 2,
    MD_UNITCELL_PBC_X         = 4,
    MD_UNITCELL_PBC_Y         = 8,
    MD_UNITCELL_PBC_Z         = 16,
    MD_UNITCELL_PBC_ALL       = 4 | 8 | 16,
};

// These flags are not specific to any distinct subtype, but can appear in both atoms, residues, bonds and whatnot.
// Where ever they make sense, they can appear. This makes it easy to propagate the flags upwards and downwards between structures
enum {
    MD_FLAG_SEQ_TERM	        = 0x1,

    MD_FLAG_HETATM              = 0x40,
    MD_FLAG_AMINO_ACID		    = 0x80,
    MD_FLAG_SIDE_CHAIN          = 0x100,
    MD_FLAG_NUCLEOTIDE	        = 0x200,
    MD_FLAG_NUCLEOBASE          = 0x400,
    MD_FLAG_NUCLEOSIDE          = 0x800,
    MD_FLAG_WATER			    = 0x1000,
    MD_FLAG_ION			        = 0x2000,
    MD_FLAG_BACKBONE            = 0x4000,

    // Experimental
    MD_FLAG_SP                  = 0x10000,
    MD_FLAG_SP2                 = 0x20000,
    MD_FLAG_SP3                 = 0x40000,
    MD_FLAG_AROMATIC            = 0x80000,
    MD_FLAG_HBOND_DONOR         = 0x100000,
    MD_FLAG_HBOND_ACCEPTOR      = 0x200000,
};

// In bonds, the order and flags are merged where the lower 4 bits encode the order and the upper 4 bits encode flags.
enum {
    MD_BOND_FLAG_AROMATIC       = 0x10,
    MD_BOND_FLAG_INTER          = 0x20,
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

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_range_t {
    int32_t beg;
    int32_t end;
} md_range_t;

typedef struct md_unitcell_t {
    double x, xy, xz;
    double y, yz;
    double z;
    uint32_t flags;
} md_unitcell_t;

// constructors for unitcell

static inline md_unitcell_t md_unitcell_from_basis_parameters(double x, double y, double z, double xy, double xz, double yz) {
    uint32_t flags = 0;
    flags |= (xy == 0.0 && xz == 0.0 && yz == 0.0) ? MD_UNITCELL_ORTHO : MD_UNITCELL_TRICLINIC;
    if (x != 0.0) flags |= MD_UNITCELL_PBC_X;
    if (y != 0.0) flags |= MD_UNITCELL_PBC_Y;
    if (z != 0.0) flags |= MD_UNITCELL_PBC_Z;

    md_unitcell_t cell = {.x = x, .xy = xy, .xz = xz, .y = y, .yz = yz, .z = z, .flags = flags};
    return cell;
}

static inline md_unitcell_t md_unitcell_from_extent(double x, double y, double z) {
    return md_unitcell_from_basis_parameters(x, y, z, 0, 0, 0);    
}

// From here: https://www.arianarab.com/post/crystal-structure-software
// Assumes that all input angles (alpha, beta, gamma) are given in degrees
static inline md_unitcell_t md_unitcell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma) {
    if (alpha == 90.0 && beta == 90.0 && gamma == 90.0) {
        return md_unitcell_from_basis_parameters(a, b, c, 0, 0, 0);
    }

    alpha = DEG_TO_RAD(alpha);
    beta  = DEG_TO_RAD(beta);
    gamma = DEG_TO_RAD(gamma);

    const double cb = cos(beta);
    const double cg = cos(gamma);
    const double sg = sin(gamma);
    const double x = (cos(alpha) - cb * cg) / sg;

    return md_unitcell_from_basis_parameters(a, b * sg, sqrt(1.0 - cb * cb - x * x), b * cg, c * cb, c * x);
}

static inline md_unitcell_t md_unitcell_from_matrix_float(const float A[3][3]) {
    return md_unitcell_from_basis_parameters(A[0][0], A[1][1], A[2][2], A[0][1], A[0][2], A[1][2]);
}

static inline md_unitcell_t md_unitcell_from_matrix_double(const double A[3][3]) {
    return md_unitcell_from_basis_parameters(A[0][0], A[1][1], A[2][2], A[0][1], A[0][2], A[1][2]);
}

// Getters and helper functionality
static inline uint32_t md_unitcell_flags(const md_unitcell_t* cell) {
    if (cell) return cell->flags;
    return 0;
}

static inline bool md_unitcell_is_triclinic(const md_unitcell_t* cell) { return cell->flags & MD_UNITCELL_TRICLINIC; }
static inline bool md_unitcell_is_orthorhombic(const md_unitcell_t* cell) { return cell->flags & MD_UNITCELL_ORTHO; }

static inline vec3_t md_unitcell_diag_vec3(const md_unitcell_t* cell) {
    if (cell) return vec3_set(cell->x, cell->y, cell->z);
    return vec3_zero();
}

static inline vec4_t md_unitcell_diag_vec4(const md_unitcell_t* cell) {
    if (cell) return vec4_set(cell->x, cell->y, cell->z, 0);
    return vec4_zero();
}

// returns the unit_cell basis matrix (A)
static inline mat3_t md_unitcell_basis_mat3(const md_unitcell_t* cell) {
    if (cell) {
        mat3_t A = {cell->x, 0, 0, cell->xy, cell->y, 0, cell->xz, cell->yz, cell->z};
        return A;
    }
    return mat3_ident();
}

static inline void md_unitcell_basis_extract(double out_A[3][3], const md_unitcell_t* cell) {
    if (cell) {
        out_A[0][0] = cell->x;
        out_A[0][1] = 0;
        out_A[0][2] = 0;
        out_A[1][0] = cell->xy;
        out_A[1][1] = cell->y;
        out_A[1][2] = 0;
        out_A[2][0] = cell->xz;
        out_A[2][1] = cell->yz;
        out_A[2][2] = cell->z;
    }
}

// returns the unit_cell basis matrix (A)
static inline mat3_t md_unitcell_inv_basis_mat3(const md_unitcell_t* cell) {
    if (cell) {
        const double i11 = 1.0 / cell->x;
        const double i22 = 1.0 / cell->y;
        const double i33 = 1.0 / cell->z;
        const double i12 = -cell->xy / (cell->x * cell->y);
        const double i13 = (cell->xy * cell->yz - cell->xz * cell->y) / (cell->x * cell->y * cell->z);
        const double i23 = -cell->yz / (cell->y * cell->z);
        mat3_t Ai = {i11, 0, 0, i12, i22, 0, i13, i23, i33};
        return Ai;
    }
    return mat3_ident();
}

static inline void md_unitcell_inv_basis_extract(double out_Ai[3][3], const md_unitcell_t* cell) {
    if (cell) {
        out_Ai[0][0] = 1.0 / cell->x;
        out_Ai[0][1] = 0;
        out_Ai[0][2] = 0;
        out_Ai[1][0] = -cell->xy / (cell->x * cell->y);
        out_Ai[1][1] = 1.0 / cell->y;
        out_Ai[1][2] = 0;
        out_Ai[2][0] = (cell->xy * cell->yz - cell->xz * cell->y) / (cell->x * cell->y * cell->z);
        out_Ai[2][1] = -cell->yz / (cell->y * cell->z);
        out_Ai[2][2] = 1.0 / cell->z;
    }
}

// returns the unitcell metric tensor G=(A^T)A
static inline mat3_t md_unitcell_G_mat3(const md_unitcell_t* cell) {
    if (cell) {
        const double g11 = cell->x * cell->x;
        const double g22 = cell->xy * cell->xy + cell->y * cell->y;
        const double g33 = cell->xz * cell->xz + cell->yz * cell->yz + cell->z * cell->z;
        const double g12 = cell->x * cell->xy;
        const double g13 = cell->x * cell->xz;
        const double g23 = cell->xy * cell->xz + cell->y * cell->yz;
        mat3_t G = {g11, g12, g13, g12, g22, g23, g13, g23, g33};
        return G;
    }
    return mat3_ident();
}

static inline void md_unitcell_extract_cell_parameters(double* a, double* b, double* c, double* alpha, double* beta, double* gamma, const md_unitcell_t* cell) {
    ASSERT(cell);

    // lengths
    *a = fabs(cell->x);
    *b = sqrt(cell->xy * cell->xy + cell->y * cell->y);
    *c = sqrt(cell->xz * cell->xz + cell->yz * cell->yz + cell->z * cell->z);

    // angles (radians)
    *alpha = acos((cell->xy * cell->xz + cell->y * cell->yz) / (*b * *c));
    *beta  = acos((cell->x * cell->xz) / (*a * *c));
    *gamma = acos((cell->x * cell->xy) / (*a * *b));
}

// Create a vec4 mask which represents the periodic dimensions from a unit cell.
// I.e. [0,1,1,0] -> periodic in y and z, but not x
static inline vec4_t md_unitcell_pbc_mask_vec4(const md_unitcell_t* unit_cell) {
    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    return vec4_set((unit_cell->flags & MD_UNITCELL_PBC_X) ? val : 0, (unit_cell->flags & MD_UNITCELL_PBC_Y) ? val : 0,
                    (unit_cell->flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);
}

static inline uint32_t md_unit_cell_flags(const md_unitcell_t* unit_cell) { return unit_cell->flags; }

typedef struct md_atom_pair_t {
    md_atom_idx_t idx[2];
} md_atom_pair_t;

typedef struct md_protein_backbone_atoms_t {
    md_atom_idx_t n;
    md_atom_idx_t ca;
    md_atom_idx_t c;
    md_atom_idx_t o;
    md_atom_idx_t hn;
} md_protein_backbone_atoms_t;

// Backbone angles
// φ (phi) is the angle in the chain C' − N − Cα − C'
// ψ (psi) is the angle in the chain N − Cα − C' − N
// https://en.wikipedia.org/wiki/Dihedral_angle#/media/File:Protein_backbone_PhiPsiOmega_drawing.svg)
typedef struct md_backbone_angles_t {
    float phi;
    float psi;
} md_backbone_angles_t;

typedef struct md_nucleic_backbone_atoms_t {
    md_atom_idx_t c5;
    md_atom_idx_t c4;
    md_atom_idx_t c3;
    md_atom_idx_t o3;
    md_atom_idx_t p;
    md_atom_idx_t o5;
} md_nucleic_backbone_atoms_t;

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
    struct md_allocator_i* alloc;
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
        const size_t len = MIN(str.len, sizeof(lbl.buf) - 1);
        MEMCPY(lbl.buf, str.ptr, len);
        lbl.len = (uint8_t)len;
    }
    return lbl;
}

// Access to substructure data
static inline void md_index_data_free (md_index_data_t* data) {
    ASSERT(data);
    if (data->offsets) {
        ASSERT(data->alloc);
        md_array_free(data->offsets,  data->alloc);
    }
    if (data->indices) {
        ASSERT(data->alloc);
        md_array_free(data->indices,  data->alloc);
    }
#if DEBUG
    MEMSET(data, 0, sizeof(md_index_data_t));
#endif
}

static inline size_t md_index_data_push_arr (md_index_data_t* data, const int32_t* index_data, size_t index_count) {
    ASSERT(data);
    ASSERT(data->alloc);
    ASSERT(index_count > 0);

    if (md_array_size(data->offsets) == 0) {
        md_array_push(data->offsets, 0, data->alloc);
    }

    size_t offset = *md_array_last(data->offsets);
    if (index_count > 0) {
        if (index_data) {
            md_array_push_array(data->indices, index_data, index_count, data->alloc);
        } else {
            md_array_grow(data->indices, md_array_size(data->indices) + index_count, data->alloc);
        }
        offset = md_array_size(data->indices);
        md_array_push(data->offsets, (uint32_t)offset, data->alloc);
    }

    return offset;
}

static inline void md_index_data_clear(md_index_data_t* data) {
    ASSERT(data);
    md_array_shrink(data->offsets, 0);
    md_array_shrink(data->indices, 0);
}

static inline size_t md_index_data_num_ranges(md_index_data_t data) {
    return data.offsets ? md_array_size(data.offsets) - 1 : 0;
}

// Access to individual substructures
static inline int32_t* md_index_range_beg(md_index_data_t data, size_t range_idx) {
    ASSERT(data.offsets && range_idx < md_array_size(data.offsets) - 1);
    return data.indices + data.offsets[range_idx];
}

static inline int32_t* md_index_range_end(md_index_data_t data, size_t range_idx) {
    ASSERT(data.offsets && range_idx < md_array_size(data.offsets) - 1);
    return data.indices + data.offsets[range_idx+1];
}

static inline size_t md_index_range_size(md_index_data_t data, size_t range_idx) {
    ASSERT(data.offsets && range_idx < md_array_size(data.offsets) - 1);
    return data.offsets[range_idx+1] - data.offsets[range_idx];
}
