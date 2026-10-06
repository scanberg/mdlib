#pragma once

#include <core/md_str.h>
#include <core/md_array.h>

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

// These flags are not specific to any distinct subtype, but can appear in both atoms, residues and whatnot.
// Where ever they make sense, they can appear. This makes it easy to propagate the flags upwards and downwards between structures
typedef enum {
    MD_FLAG_NONE                = 0,
    MD_FLAG_COARSE_GRAINED      = 0x1,      // Coarse grained

    MD_FLAG_POLYMER             = 0x2,      // Flag for connected polymers
    MD_FLAG_BACKBONE            = 0x4,      // Backbone atoms
    MD_FLAG_DERIVED             = 0x8,      // Not defined by the source but derived when loading, e.g. entities (and thereby instances) from md_util_system_infer_entity_and_instance
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

    // Chemistry of an atom, see md_chem.h. Hybridization of heavy atoms: SP2 includes the N and O whose lone pair
    // takes part in a pi system (amide and aniline N, pyrrole N, furan O). AROMATIC: in an aromatic ring.
    MD_FLAG_SP                  = 0x1000000,
    MD_FLAG_SP2                 = 0x2000000,
    MD_FLAG_SP3                 = 0x4000000,
    MD_FLAG_AROMATIC            = 0x8000000,

    // Hydrogen bond roles, see md_hbond.h (md_hbond_infer_atom_flags)
    MD_FLAG_HBOND_DONOR         = 0x10000000,
    MD_FLAG_HBOND_ACCEPTOR      = 0x20000000,

    // A massless interaction site without an element: the M site of 4 site water models (TIP4P, OPC), the lone pairs
    // of 5 site models (TIP5P). It belongs to its molecule but forms no covalent bonds. Set on the atom type.
    MD_FLAG_VIRTUAL_SITE        = 0x40000000,
} md_flags_t;

ENUM_FLAGS(md_flags_t)

// Bond flags, on two axes: the low byte describes the bond chemically, the bits above it where the bond and its order
// came from. Only the low byte is part of a bond's identity when comparing structures. Each field is either one bit or
// one value, so every combination describes a case (AROMATIC with DELOCALIZED aside, which never occur together).
//
// CHEMISTRY
//   ORDER        A value in bits 0-2, not a set of bits: read it with md_bond_order and write it with
//                md_bond_flags_set_order. 0 is unknown (not perceived, and not given by the file), which is not the
//                same as single.
//   AROMATIC     In an aromatic ring.
//   DELOCALIZED  Resonance outside of rings: carboxylate, guanidinium, nitro, phosphate, ... AROMATIC and DELOCALIZED
//                come on top of an order which holds one Kekule structure: a benzene bond is AROMATIC with order 1 or
//                2, the two C-O of a carboxylate are DELOCALIZED with orders 2 and 1.
//   COORDINATE   A coordination (dative) bond, which counts towards the valence of neither atom: every bond between a
//                metal and a non-metal, whatever its origin (md_util_system_infer_coordination), and any other bond a
//                file says is one. A bond without it is covalent, or metallic between two metals: whether a bond
//                involves a metal is a question for the elements of its atoms.
//
// PROVENANCE
//   ORIGIN           A value in bits 8-9 (md_bond_origin_t): read it with md_bond_origin and write it with
//                    md_bond_flags_set_origin.
//   ORDER_PERCEIVED  The order, AROMATIC and DELOCALIZED were perceived (md_chem_perceive), not given. A known order
//                    without it was given, by the file or by hand, and stays as it is when the chemistry is perceived
//                    again; one with it is perceived anew.
typedef enum {
    MD_BOND_FLAG_NONE            = 0,
    MD_BOND_FLAG_ORDER_MASK      = 0x7,     // See md_bond_order
    MD_BOND_FLAG_AROMATIC        = 0x8,     // In an aromatic ring
    MD_BOND_FLAG_DELOCALIZED     = 0x10,    // Resonance outside of rings
    MD_BOND_FLAG_COORDINATE      = 0x20,    // Coordination (dative), always between a metal and a non-metal

    MD_BOND_FLAG_ORIGIN_MASK     = 0x300,   // See md_bond_origin
    MD_BOND_FLAG_ORDER_PERCEIVED = 0x400,   // The order (and AROMATIC, DELOCALIZED) was perceived, not given
} md_bond_flags_t;

ENUM_FLAGS(md_bond_flags_t)

enum {
    MD_BOND_ORDER_UNKNOWN   = 0,
    MD_BOND_ORDER_SINGLE    = 1,
    MD_BOND_ORDER_DOUBLE    = 2,
    MD_BOND_ORDER_TRIPLE    = 3,
    MD_BOND_ORDER_QUADRUPLE = 4,
};

// Where a bond came from. FILE is 0, so a bond added without saying where it came from counts as given.
typedef enum md_bond_origin_t {
    MD_BOND_ORIGIN_FILE     = 0,    // Given by the file, for some of its atoms (CONECT records, a component's bonds)
    MD_BOND_ORIGIN_TOPOLOGY = 1,    // A force field topology: all of the bonds of the atoms it covers
    MD_BOND_ORIGIN_USER     = 2,    // Added by hand
    MD_BOND_ORIGIN_INFERRED = 3,    // Inferred from the geometry, and replaced when inferred again
} md_bond_origin_t;

#define MD_BOND_ORDER_SHIFT  0
#define MD_BOND_ORIGIN_SHIFT 8

static inline int md_bond_order(md_bond_flags_t flags) {
    return (int)((flags & MD_BOND_FLAG_ORDER_MASK) >> MD_BOND_ORDER_SHIFT);
}

static inline md_bond_flags_t md_bond_flags_set_order(md_bond_flags_t flags, int order) {
    const unsigned int value = (order < 0 ? 0u : (order > 4 ? 4u : (unsigned int)order)) << MD_BOND_ORDER_SHIFT;
    return (md_bond_flags_t)(((unsigned int)flags & ~(unsigned int)MD_BOND_FLAG_ORDER_MASK) | value);
}

static inline md_bond_origin_t md_bond_origin(md_bond_flags_t flags) {
    return (md_bond_origin_t)(((unsigned int)flags & (unsigned int)MD_BOND_FLAG_ORIGIN_MASK) >> MD_BOND_ORIGIN_SHIFT);
}

static inline md_bond_flags_t md_bond_flags_set_origin(md_bond_flags_t flags, md_bond_origin_t origin) {
    const unsigned int value = ((unsigned int)origin << MD_BOND_ORIGIN_SHIFT) & (unsigned int)MD_BOND_FLAG_ORIGIN_MASK;
    return (md_bond_flags_t)(((unsigned int)flags & ~(unsigned int)MD_BOND_FLAG_ORIGIN_MASK) | value);
}

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
