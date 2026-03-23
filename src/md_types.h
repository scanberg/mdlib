#pragma once

#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>

#include <stdint.h>

#if DEBUG
#include <core/md_log.h>
#endif

// Secondary structure types (DSSP-like)
typedef enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL,
    MD_SECONDARY_STRUCTURE_TURN,
    MD_SECONDARY_STRUCTURE_BEND,
    MD_SECONDARY_STRUCTURE_HELIX_310,
    MD_SECONDARY_STRUCTURE_HELIX_ALPHA,
    MD_SECONDARY_STRUCTURE_HELIX_PI,
    MD_SECONDARY_STRUCTURE_BETA_SHEET,
    MD_SECONDARY_STRUCTURE_BETA_BRIDGE,
} md_secondary_structure_t;

enum {
    MD_RAMACHANDRAN_TYPE_UNKNOWN = 0,
    MD_RAMACHANDRAN_TYPE_GENERAL,
    MD_RAMACHANDRAN_TYPE_GLYCINE,
    MD_RAMACHANDRAN_TYPE_PROLINE,
    MD_RAMACHANDRAN_TYPE_PREPROL,
};

// Cell information
typedef enum {
    MD_UNITCELL_NONE          = 0,
    MD_UNITCELL_ORTHO         = 1,
    MD_UNITCELL_TRICLINIC     = 2,
    MD_UNITCELL_PBC_X         = 4,
    MD_UNITCELL_PBC_Y         = 8,
    MD_UNITCELL_PBC_Z         = 16,
    MD_UNITCELL_PBC_ALL       = 4 | 8 | 16,
} md_unitcell_flags_t;

ENUM_FLAGS(md_unitcell_flags_t)

// These flags are not specific to any distinct subtype, but can appear in both atoms, residues and whatnot.
// Where ever they make sense, they can appear. This makes it easy to propagate the flags upwards and downwards between structures
typedef enum {
    MD_FLAG_NONE                = 0,
    MD_FLAG_COARSE_GRAINED      = 0x1,      // Coarse grained representation

    MD_FLAG_POLYMER             = 0x2,      // Flag for connected polymers
    MD_FLAG_BACKBONE            = 0x4,      // Backbone atoms
    MD_FLAG_TERMINAL_BEG        = 0x10,     // Terminal atoms (N and C terminus in proteins, 5' and 3' in nucleic acids)
    MD_FLAG_TERMINAL_END        = 0x20,     // Terminal atoms (N and C terminus in proteins, 5' and 3' in nucleic acids)

    // Proteins, Nucleic acids and their components
    MD_FLAG_POLYPEPTIDE         = 0x100,    // Top level for connected polypeptide chains (i.e. proteins)
    MD_FLAG_AMINO_ACID		    = 0x200,    // For expressing an amino acid monomer in the polypeptide chain
    MD_FLAG_SIDE_CHAIN          = 0x400,    // Amino acid side chain atoms

    MD_FLAG_NUCLEIC_ACID        = 0x800,    // Top level for connected nucleic acid chains (i.e. DNA, RNA)
    MD_FLAG_NUCLEOTIDE          = 0x1000,   // For expressing a nucleotide monomer in the nucleic acid chain
    MD_FLAG_NUCLEOSIDE          = 0x2000,   // Nucleoside component of a nucleotide
    MD_FLAG_NUCLEOBASE          = 0x4000,   // Nucleobase component of a nucleotide

    // HETERO types
    MD_FLAG_HETERO              = 0x10000,
    MD_FLAG_WATER			    = 0x20000,
    MD_FLAG_ION			        = 0x40000,

    // Chirality
    MD_FLAG_ISOMER_L            = 0x100000,
    MD_FLAG_ISOMER_D            = 0x200000,

    // Experimental
    MD_FLAG_SP                  = 0x1000000,
    MD_FLAG_SP2                 = 0x2000000,
    MD_FLAG_SP3                 = 0x4000000,
    MD_FLAG_AROMATIC            = 0x8000000,

//    MD_FLAG_HBOND_DONOR         = 0x1000000,
//    MD_FLAG_HBOND_ACCEPTOR      = 0x2000000,
} md_flags_t;

ENUM_FLAGS(md_flags_t)

// In bonds, the order and flags are merged where the lower 4 bits encode the order and the upper 4 bits encode flags.
typedef enum {
    MD_BOND_FLAG_COVALENT       = 0x1,
    MD_BOND_FLAG_DOUBLE         = 0x2,
    MD_BOND_FLAG_TRIPLE         = 0x4,
    MD_BOND_FLAG_QUADRUPLE      = 0x8,
    MD_BOND_FLAG_AROMATIC       = 0x10,
    MD_BOND_FLAG_COORDINATE     = 0x20, // Coordinate / Dative
    MD_BOND_FLAG_METAL          = 0x40, // Involves a metal atom
	MD_BOND_FLAG_USER_DEFINED   = 0x80, // User defined bond
} md_bond_flags_t;

ENUM_FLAGS(md_bond_flags_t)

// Atomic number constants for all elements (Z values)
enum {
    MD_Z_X  = 0,   // Unknown
    MD_Z_H  = 1,   // Hydrogen
    MD_Z_He = 2,   // Helium
    MD_Z_Li = 3,   // Lithium
    MD_Z_Be = 4,   // Beryllium
    MD_Z_B  = 5,   // Boron
    MD_Z_C  = 6,   // Carbon
    MD_Z_N  = 7,   // Nitrogen
    MD_Z_O  = 8,   // Oxygen
    MD_Z_F  = 9,   // Fluorine
    MD_Z_Ne = 10,  // Neon
    MD_Z_Na = 11,  // Sodium
    MD_Z_Mg = 12,  // Magnesium
    MD_Z_Al = 13,  // Aluminium
    MD_Z_Si = 14,  // Silicon
    MD_Z_P  = 15,  // Phosphorus
    MD_Z_S  = 16,  // Sulfur
    MD_Z_Cl = 17,  // Chlorine
    MD_Z_Ar = 18,  // Argon
    MD_Z_K  = 19,  // Potassium
    MD_Z_Ca = 20,  // Calcium
    MD_Z_Sc = 21,  // Scandium
    MD_Z_Ti = 22,  // Titanium
    MD_Z_V  = 23,  // Vanadium
    MD_Z_Cr = 24,  // Chromium
    MD_Z_Mn = 25,  // Manganese
    MD_Z_Fe = 26,  // Iron
    MD_Z_Co = 27,  // Cobalt
    MD_Z_Ni = 28,  // Nickel
    MD_Z_Cu = 29,  // Copper
    MD_Z_Zn = 30,  // Zinc
    MD_Z_Ga = 31,  // Gallium
    MD_Z_Ge = 32,  // Germanium
    MD_Z_As = 33,  // Arsenic
    MD_Z_Se = 34,  // Selenium
    MD_Z_Br = 35,  // Bromine
    MD_Z_Kr = 36,  // Krypton
    MD_Z_Rb = 37,  // Rubidium
    MD_Z_Sr = 38,  // Strontium
    MD_Z_Y  = 39,  // Yttrium
    MD_Z_Zr = 40,  // Zirconium
    MD_Z_Nb = 41,  // Niobium
    MD_Z_Mo = 42,  // Molybdenum
    MD_Z_Tc = 43,  // Technetium
    MD_Z_Ru = 44,  // Ruthenium
    MD_Z_Rh = 45,  // Rhodium
    MD_Z_Pd = 46,  // Palladium
    MD_Z_Ag = 47,  // Silver
    MD_Z_Cd = 48,  // Cadmium
    MD_Z_In = 49,  // Indium
    MD_Z_Sn = 50,  // Tin
    MD_Z_Sb = 51,  // Antimony
    MD_Z_Te = 52,  // Tellurium
    MD_Z_I  = 53,  // Iodine
    MD_Z_Xe = 54,  // Xenon
    MD_Z_Cs = 55,  // Caesium
    MD_Z_Ba = 56,  // Barium
    MD_Z_La = 57,  // Lanthanum
    MD_Z_Ce = 58,  // Cerium
    MD_Z_Pr = 59,  // Praseodymium
    MD_Z_Nd = 60,  // Neodymium
    MD_Z_Pm = 61,  // Promethium
    MD_Z_Sm = 62,  // Samarium
    MD_Z_Eu = 63,  // Europium
    MD_Z_Gd = 64,  // Gadolinium
    MD_Z_Tb = 65,  // Terbium
    MD_Z_Dy = 66,  // Dysprosium
    MD_Z_Ho = 67,  // Holmium
    MD_Z_Er = 68,  // Erbium
    MD_Z_Tm = 69,  // Thulium
    MD_Z_Yb = 70,  // Ytterbium
    MD_Z_Lu = 71,  // Lutetium
    MD_Z_Hf = 72,  // Hafnium
    MD_Z_Ta = 73,  // Tantalum
    MD_Z_W  = 74,  // Tungsten
    MD_Z_Re = 75,  // Rhenium
    MD_Z_Os = 76,  // Osmium
    MD_Z_Ir = 77,  // Iridium
    MD_Z_Pt = 78,  // Platinum
    MD_Z_Au = 79,  // Gold
    MD_Z_Hg = 80,  // Mercury
    MD_Z_Tl = 81,  // Thallium
    MD_Z_Pb = 82,  // Lead
    MD_Z_Bi = 83,  // Bismuth
    MD_Z_Po = 84,  // Polonium
    MD_Z_At = 85,  // Astatine
    MD_Z_Rn = 86,  // Radon
    MD_Z_Fr = 87,  // Francium
    MD_Z_Ra = 88,  // Radium
    MD_Z_Ac = 89,  // Actinium
    MD_Z_Th = 90,  // Thorium
    MD_Z_Pa = 91,  // Protactinium
    MD_Z_U  = 92,  // Uranium
    MD_Z_Np = 93,  // Neptunium
    MD_Z_Pu = 94,  // Plutonium
    MD_Z_Am = 95,  // Americium
    MD_Z_Cm = 96,  // Curium
    MD_Z_Bk = 97,  // Berkelium
    MD_Z_Cf = 98,  // Californium
    MD_Z_Es = 99,  // Einsteinium
    MD_Z_Fm = 100, // Fermium
    MD_Z_Md = 101, // Mendelevium
    MD_Z_No = 102, // Nobelium
    MD_Z_Lr = 103, // Lawrencium
    MD_Z_Rf = 104, // Rutherfordium
    MD_Z_Db = 105, // Dubnium
    MD_Z_Sg = 106, // Seaborgium
    MD_Z_Bh = 107, // Bohrium
    MD_Z_Hs = 108, // Hassium
    MD_Z_Mt = 109, // Meitnerium
    MD_Z_Ds = 110, // Darmstadtium
    MD_Z_Rg = 111, // Roentgenium
    MD_Z_Cn = 112, // Copernicium
    MD_Z_Nh = 113, // Nihonium
    MD_Z_Fl = 114, // Flerovium
    MD_Z_Mc = 115, // Moscovium
    MD_Z_Lv = 116, // Livermorium
    MD_Z_Ts = 117, // Tennessine
    MD_Z_Og = 118, // Oganesson
    MD_Z_Count
};

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_component_idx_t;
typedef int32_t     md_instance_idx_t;
typedef int32_t     md_entity_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_sequence_id_t;
typedef int32_t     md_bond_idx_t;
typedef uint16_t    md_atom_type_idx_t;
typedef uint8_t     md_atomic_number_t;
typedef uint8_t     md_ramachandran_type_t;
typedef uint8_t     md_order_t;

// For backwards compatibility
typedef uint8_t md_element_t;  // Atomic number (1-118), 0 is unknown

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_irange_t {
    int32_t beg;
    int32_t end;
} md_irange_t;

typedef struct md_urange_t {
    uint32_t beg;
    uint32_t end;
} md_urange_t;

typedef struct md_unitcell_t {
    double x, xy, xz;
    double y, yz;
    double z;
    md_unitcell_flags_t flags;
} md_unitcell_t;

typedef struct md_atom_pair_t {
    md_atom_idx_t idx[2];
} md_atom_pair_t;

typedef struct md_amino_acid_atoms_t {
    // Backbone
    md_atom_idx_t n;
    md_atom_idx_t ca;
    md_atom_idx_t c;
	// Side chain and appendages
    md_atom_idx_t o;
    md_atom_idx_t cb;
    md_atom_idx_t hn;
} md_amino_acid_atoms_t;

// Backbone angles
// φ (phi) is the angle in the chain C' − N − Cα − C'
// ψ (psi) is the angle in the chain N − Cα − C' − N
// https://en.wikipedia.org/wiki/Dihedral_angle#/media/File:Protein_backbone_PhiPsiOmega_drawing.svg)
typedef struct md_backbone_angles_t {
    float phi;
    float psi;
} md_backbone_angles_t;

typedef struct md_nucleic_acid_atoms_t {
    // Backbone
    md_atom_idx_t p;
    md_atom_idx_t o5;
    md_atom_idx_t c5;
    md_atom_idx_t c4;
    md_atom_idx_t c3;
    md_atom_idx_t o3;
    // Nucleoside
    md_atom_idx_t c1;
    md_atom_idx_t c2;
    md_atom_idx_t o4;
} md_nucleic_acid_atoms_t;

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

#ifdef __cplusplus
extern "C" {
#endif

// Element property functions
str_t md_atomic_number_name(md_atomic_number_t z);
str_t md_atomic_number_symbol(md_atomic_number_t z);
float md_atomic_number_mass(md_atomic_number_t z);
float md_atomic_number_vdw_radius(md_atomic_number_t z);
float md_atomic_number_covalent_radius(md_atomic_number_t z);
int   md_atomic_number_max_valence(md_atomic_number_t z);
uint32_t md_atomic_number_cpk_color(md_atomic_number_t z);

// Element symbol and name lookup functions
md_atomic_number_t md_atomic_number_from_symbol(str_t sym, bool ignore_case);

// Infer atomic number from other properties, this is a heuristic and may fail. 0 is returned on failure
// Per-atom inference from labels (atom name + residue)
// res_name is the name of the component of the atom [optional]
// res_size is the size of the component of the atom [optional]
md_atomic_number_t md_atomic_number_infer_from_label(str_t atom_name, str_t res_name, size_t res_size);
md_atomic_number_t md_atomic_number_infer_from_mass(float mass);

// Batch form wired to molecule structure
// Returns the number of successfully inferred atomic numbers
// All of the supplied arrays must have at least 'count' elements
size_t md_atomic_number_infer_from_mass_batch(md_atomic_number_t out_z[], const float masses[], size_t count);

#ifdef __cplusplus
}
#endif

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

static inline size_t md_index_data_num_ranges(const md_index_data_t* data) {
    ASSERT(data);
    size_t num_ranges = 0;
    if (data->offsets) {
        num_ranges = md_array_size(data->offsets) - 1;
    }
    return num_ranges;
}

// Access to individual substructures
static inline const int32_t* md_index_range_beg(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    const int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx < md_array_size(data->offsets) - 1) {
        ptr = data->indices + data->offsets[range_idx];
    }
    return ptr;
}

static inline const int32_t* md_index_range_end(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    const int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx < md_array_size(data->offsets) - 1) {
        ptr = data->indices + data->offsets[range_idx + 1];
    }
    return ptr;
}

// This is a mutable version of md_index_range_beg
static inline int32_t* md_index_range_ptr(md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx < md_array_size(data->offsets) - 1) {
        ptr = data->indices + data->offsets[range_idx];
    }
    return ptr;
}

static inline size_t md_index_range_size(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    size_t size = 0;
    if (data->offsets && range_idx < md_array_size(data->offsets) - 1) {
        size = data->offsets[range_idx+1] - data->offsets[range_idx];
    }
    return size;
}

static inline void md_index_data_merge(md_index_data_t* dest, const md_index_data_t* src) {
    ASSERT(dest);
    ASSERT(src);
    ASSERT(dest->alloc == src->alloc);
    size_t dest_num_indices = md_array_size(dest->indices);
    size_t src_num_ranges = md_index_data_num_ranges(src);
    for (size_t i = 0; i < src_num_ranges; ++i) {
        uint32_t beg_offset = src->offsets[i];
        uint32_t end_offset = src->offsets[i + 1];
        size_t range_size = end_offset - beg_offset;
        // Push adjusted offsets
        uint32_t new_offset = (uint32_t)(dest_num_indices + range_size);
        md_array_push(dest->offsets, new_offset, dest->alloc);
        // Push indices
        for (uint32_t j = beg_offset; j < end_offset; ++j) {
            int32_t idx = src->indices[j];
            md_array_push(dest->indices, idx, dest->alloc);
        }
        dest_num_indices += range_size;
    }
}


// Atom coordinate helper functions

#define MD_COORD_CHUNK_SIZE 4

// Store atom coordinates in chunks of 4 to strike a good balance between cache locality for fetching individual atoms and SIMD processing
typedef struct md_coord_chunk_t {
    ALIGNAS(16) float x[MD_COORD_CHUNK_SIZE], y[MD_COORD_CHUNK_SIZE], z[MD_COORD_CHUNK_SIZE];
} md_coord_chunk_t;

// This represents the atom coordinates for a given state (frame) of the system.
typedef struct md_atom_coord_data_t {
    size_t num_chunks;
    md_coord_chunk_t* chunks;
} md_atom_coord_data_t;

// This represents the transient portion of a system which change over time, such as the atom coordinates, bonds, etc.
typedef struct md_system_state_t {
    md_atom_coord_data_t coord;
    md_unitcell_t unitcell;
    // Add other dynamic properties here such as e.g. hydrogen bonds, etc. which may change over time.
} md_system_state_t;

static inline void md_atom_coord_init(md_atom_coord_data_t* coord, size_t num_atoms, md_allocator_i* alloc) {
    ASSERT(coord);
    ASSERT(alloc);
    coord->num_chunks = DIV_UP(num_atoms, MD_COORD_CHUNK_SIZE);
    md_array_resize(coord->chunks, coord->num_chunks, alloc);
}

static inline void md_atom_coord_free(md_atom_coord_data_t* coord, md_allocator_i* alloc) {
    ASSERT(coord);
    ASSERT(alloc);
    if (coord->chunks) {
        md_array_free(coord->chunks, alloc);
        coord->chunks = NULL;
        coord->num_chunks = 0;
    }
}

static inline vec3_t md_atom_coord_vec3(const md_atom_coord_data_t* coord, size_t atom_idx) {
    ASSERT(coord);
    STATIC_ASSERT(MD_COORD_CHUNK_SIZE == 4, "Unexpected chunk size");

    size_t chunk_idx = atom_idx >> 2;
    size_t within_chunk_idx = atom_idx & 3;
    if (chunk_idx < coord->num_chunks) {
        const md_coord_chunk_t* chunk = &coord->chunks[chunk_idx];
        return vec3_set(chunk->x[within_chunk_idx], chunk->y[within_chunk_idx], chunk->z[within_chunk_idx]);
    }
    return vec3_zero();
}

static inline vec4_t md_atom_coord_vec4(const md_atom_coord_data_t* coord, float w, size_t atom_idx) {
    ASSERT(coord);
    STATIC_ASSERT(MD_COORD_CHUNK_SIZE == 4, "Unexpected chunk size");

    size_t chunk_idx = atom_idx >> 2;
    size_t within_chunk_idx = atom_idx & 3;
    if (chunk_idx < coord->num_chunks) {
        const md_coord_chunk_t* chunk = &coord->chunks[chunk_idx];
        return vec4_set(chunk->x[within_chunk_idx], chunk->y[within_chunk_idx], chunk->z[within_chunk_idx], w);
    }
    return vec4_zero();
}

static inline void md_atom_coord_set(md_atom_coord_data_t* coord, size_t atom_idx, float x, float y, float z) {
    ASSERT(coord);
    STATIC_ASSERT(MD_COORD_CHUNK_SIZE == 4, "Unexpected chunk size");

    size_t chunk_idx = atom_idx >> 2;
    size_t within_chunk_idx = atom_idx & 3;
    if (chunk_idx < coord->num_chunks) {
        md_coord_chunk_t* chunk = &coord->chunks[chunk_idx];
        chunk->x[within_chunk_idx] = x;
        chunk->y[within_chunk_idx] = y;
        chunk->z[within_chunk_idx] = z;
    }
}

// Atom coordinate iterators (vec4_t, md_128, md_256)
// These are useful for iterating over atom coordinates in a SIMD friendly way

typedef struct md_atom_coord_iter_linear_t {
    const md_atom_coord_data_t* coord;
    size_t it;
    size_t count;
} md_atom_coord_iter_linear_t;

typedef struct md_atom_coord_iter_indexed_t {
    const md_atom_coord_data_t* coord;
    size_t count;
    const int32_t* idx;
    size_t it;
} md_atom_coord_iter_indexed_t;

static inline md_atom_coord_iter_linear_t md_atom_coord_iter_linear_create(const md_atom_coord_data_t* coord, size_t count) {
    ASSERT(coord);
    md_atom_coord_iter_linear_t iter = {
        .coord = coord,
        .it = 0,
        .count = count,
    };
    return iter;
}

static inline md_atom_coord_iter_indexed_t md_atom_coord_iter_indexed_create(const md_atom_coord_data_t* coord, const int32_t* idx, size_t count) {
    ASSERT(coord);
    ASSERT(idx);
    md_atom_coord_iter_indexed_t iter = {
        .coord = coord,
        .count = count,
        .idx = idx,
        .it = 0,
    };
    return iter;
}

static inline bool md_atom_coord_iter_linear_next_xyz(float* out_xyz, md_atom_coord_iter_linear_t* iter) {
    ASSERT(iter);
    ASSERT(out_xyz);
    if (iter->it < iter->count) {
        size_t idx = iter->it;
        size_t chunk_idx = idx >> 2;
        size_t within_chunk_idx = idx & 3;
        ASSERT(chunk_idx < iter->coord->num_chunks);
        const md_coord_chunk_t* chunk = &iter->coord->chunks[chunk_idx];
        out_xyz[0] = chunk->x[within_chunk_idx];
        out_xyz[1] = chunk->y[within_chunk_idx];
        out_xyz[2] = chunk->z[within_chunk_idx];
        iter->it++;
        return true;
    }
    return false;
}

static inline bool md_atom_coord_iter_indexed_next_xyz(float* out_xyz, md_atom_coord_iter_indexed_t* iter) {
    ASSERT(iter);
    ASSERT(out_xyz);
    if (iter->it < iter->count) {
        int idx = iter->idx[iter->it];
        int chunk_idx = idx / 4;
        int within_chunk_idx = idx % 4;
        ASSERT(0 <= chunk_idx && chunk_idx < iter->coord->num_chunks);
        const md_coord_chunk_t* chunk = &iter->coord->chunks[chunk_idx];
        out_xyz[0] = chunk->x[within_chunk_idx];
        out_xyz[1] = chunk->y[within_chunk_idx];
        out_xyz[2] = chunk->z[within_chunk_idx];
        iter->it++;
        return true;
    }
    return false;
}

// Returns number of valid lanes [0..4]. Assumes coordinate data is padded so loads are always safe.
// The caller uses the returned lane count for tail handling.
static inline int md_atom_coord_iter_linear_next_chunk_128(md_128* out_x, md_128* out_y, md_128* out_z,
    md_atom_coord_iter_linear_t* iter) {
    ASSERT(iter && out_x && out_y && out_z);

    if (iter->it >= iter->count) {
        return 0;
    }

    ASSERT((iter->it & 3) == 0); // chunk aligned (since offsets are disallowed)

    const size_t remaining = iter->count - iter->it;
    const int lanes = MIN(4, (int)remaining);

    const size_t chunk_idx = iter->it >> 2;
    ASSERT(chunk_idx < iter->coord->num_chunks);
    const md_coord_chunk_t* c = &iter->coord->chunks[chunk_idx];

    // Always safe due to padding (even for lanes < 4)
    *out_x = md_mm_load_ps(c->x);
    *out_y = md_mm_load_ps(c->y);
    *out_z = md_mm_load_ps(c->z);

    iter->it += 4;
    return lanes;
}

static inline int md_atom_coord_iter_linear_next_chunk_256(md_256* out_x, md_256* out_y, md_256* out_z,
    md_atom_coord_iter_linear_t* iter) {
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
    const md_coord_chunk_t* c0 = &iter->coord->chunks[chunk_idx0];
    const md_coord_chunk_t* c1 = &iter->coord->chunks[chunk_idx1];

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

static inline int md_atom_coord_iter_indexed_next_chunk_128(md_128* out_x, md_128* out_y, md_128* out_z,
    md_atom_coord_iter_indexed_t* iter) {
    ASSERT(iter && out_x && out_y && out_z);

    if (iter->it >= iter->count) return 0;

    const int remaining = (int)iter->count - (int)iter->it;
    const int lanes = MIN(4, remaining);

    // Base for gathers: treat the entire `chunks` array as a flat float array.
    const float* base = (const float*)&iter->coord->chunks[0].x[0];
    ASSERT(sizeof(md_coord_chunk_t) / sizeof(float) == 12); // x[4],y[4],z[4] => 12 floats

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

#include <md_unitcell.inl>
