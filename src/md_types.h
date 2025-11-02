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
    MD_FLAG_COARSE_GRAINED      = 0x1,  // Coarse grained representation

    // Polymers
    MD_FLAG_POLYMER             = 0x10,     // Flag for connected polymers
    MD_FLAG_BACKBONE            = 0x20,     // Backbone atoms

    // Proteins
    MD_FLAG_AMINO_ACID		    = 0x40,
    MD_FLAG_SIDE_CHAIN          = 0x80,

    // RNA DNA
    MD_FLAG_NUCLEOTIDE	        = 0x100,
    MD_FLAG_NUCLEOBASE          = 0x200,
    MD_FLAG_NUCLEOSIDE          = 0x400,

    // HETERO types
    MD_FLAG_HETERO              = 0x1000,
    MD_FLAG_WATER			    = 0x2000,
    MD_FLAG_ION			        = 0x4000,

    // Chirality
    MD_FLAG_ISOMER_L            = 0x10000,
    MD_FLAG_ISOMER_D            = 0x20000,

    // Experimental
    MD_FLAG_SP                  = 0x100000,
    MD_FLAG_SP2                 = 0x200000,
    MD_FLAG_SP3                 = 0x400000,
    MD_FLAG_AROMATIC            = 0x800000,

//    MD_FLAG_HBOND_DONOR         = 0x1000000,
//    MD_FLAG_HBOND_ACCEPTOR      = 0x2000000,
};

// In bonds, the order and flags are merged where the lower 4 bits encode the order and the upper 4 bits encode flags.
enum {
    MD_BOND_FLAG_COVALENT       = 0x1,
    MD_BOND_FLAG_DOUBLE         = 0x2,
    MD_BOND_FLAG_TRIPLE         = 0x4,
    MD_BOND_FLAG_QUADRUPLE      = 0x8,
    MD_BOND_FLAG_AROMATIC       = 0x10,
    MD_BOND_FLAG_INTER          = 0x20, // Inter residue / component
    MD_BOND_FLAG_COORDINATE     = 0x40, // Coordinate / Dative
    MD_BOND_FLAG_METAL          = 0x80, // Involves a metal atom
};

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
};

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_comp_idx_t;
typedef int32_t     md_inst_idx_t;
typedef int32_t     md_entity_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_seq_id_t;
typedef int32_t     md_bond_idx_t;
typedef uint16_t    md_atom_type_idx_t;
typedef uint32_t    md_secondary_structure_t;
typedef uint32_t    md_flags_t;
typedef uint8_t     md_atomic_number_t;
typedef uint8_t     md_ramachandran_type_t;
typedef uint8_t     md_valence_t;
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
    uint32_t flags;
} md_unitcell_t;

static inline size_t md_unitcell_print(char* out_buf, size_t buf_cap, const md_unitcell_t* cell) {
    int len = snprintf(out_buf, buf_cap, "x: %f, y: %f, z: %f \nxy: %f, xz: %f, yz: %f \nflags: %i", cell->x, cell->y, cell->z, cell->xy, cell->xz, cell->yz, cell->flags);
    return len;
}

// constructors for unitcell

static inline md_unitcell_t md_unitcell_none(void) {
    md_unitcell_t cell = {0};
    return cell;
}

static inline md_unitcell_t md_unitcell_from_basis_parameters(double x, double y, double z, double xy, double xz, double yz) {
    uint32_t flags = 0;

    if (xy == 0.0 && xz == 0.0 && yz == 0.0) {
        if (!(x == 0.0 && y == 0.0 && z == 0.0) && !(x == 1.0 && y == 1.0 && z == 1.0)) {
            flags |= MD_UNITCELL_ORTHO;
        }
    } else {
        flags |= MD_UNITCELL_TRICLINIC;
    }

    if (flags) {
        if (x != 0.0) flags |= MD_UNITCELL_PBC_X;
        if (y != 0.0) flags |= MD_UNITCELL_PBC_Y;
        if (z != 0.0) flags |= MD_UNITCELL_PBC_Z;
    }

    md_unitcell_t cell = {.x = x, .xy = xy, .xz = xz, .y = y, .yz = yz, .z = z, .flags = flags};
    return cell;
}

static inline md_unitcell_t md_unitcell_from_extent(double x, double y, double z) {
    return md_unitcell_from_basis_parameters(x, y, z, 0, 0, 0);    
}

// From here: https://www.arianarab.com/post/crystal-structure-software
// Assumes that all input angles (alpha, beta, gamma) are given in degrees
static inline md_unitcell_t md_unitcell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma) {
    if (a == 0 && b == 0 && c == 0 && alpha == 0 && beta == 0 && gamma == 0) {
        return md_unitcell_none();
    }

    if (alpha == 90.0 && beta == 90.0 && gamma == 90.0) {
        return md_unitcell_from_basis_parameters(a, b, c, 0, 0, 0);
    }

    alpha = DEG_TO_RAD(alpha);
    beta  = DEG_TO_RAD(beta);
    gamma = DEG_TO_RAD(gamma);

    // https://docs.lammps.org/Howto_triclinic.html
    double x = a;
    double xy = b * cos(gamma);
    double xz = c * cos(beta);
    double y = b * sin(gamma);
    double yz = (b * c * cos(alpha) - xy * xz) / y;
    double z = sqrt(c * c - xz * xz - yz * yz);

    return md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
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
    if (cell) return vec3_set((float)cell->x, (float)cell->y, (float)cell->z);
    return vec3_zero();
}

static inline vec4_t md_unitcell_diag_vec4(const md_unitcell_t* cell) {
    if (cell) return vec4_set((float)cell->x, (float)cell->y, (float)cell->z, 0);
    return vec4_zero();
}

// returns the unit_cell basis matrix (A)
static inline mat3_t md_unitcell_basis_mat3(const md_unitcell_t* cell) {
    if (cell) {
        mat3_t A = {(float)cell->x, 0, 0, (float)cell->xy, (float)cell->y, 0, (float)cell->xz, (float)cell->yz, (float)cell->z};
        return A;
    }
    return mat3_ident();
}

// returns the unit_cell basis matrix (A)
static inline mat4_t md_unitcell_basis_mat4(const md_unitcell_t* cell) {
    if (cell) {
        mat4_t A = {
            (float)cell->x, 0, 0, 0,
            (float)cell->xy, (float)cell->y, 0, 0,
            (float)cell->xz, (float)cell->yz, (float)cell->z, 0,
            0, 0, 0, 0
        };
        return A;
    }
    return mat4_ident();
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
        if (!cell->flags) {
            mat3_t zero = {0};
            return zero;
        }
        const double i11 = 1.0 / cell->x;
        const double i22 = 1.0 / cell->y;
        const double i33 = 1.0 / cell->z;
        const double i12 = -cell->xy / (cell->x * cell->y);
        const double i13 = (cell->xy * cell->yz - cell->xz * cell->y) / (cell->x * cell->y * cell->z);
        const double i23 = -cell->yz / (cell->y * cell->z);
        mat3_t Ai = {(float)i11, 0, 0, (float)i12, (float)i22, 0, (float)i13, (float)i23, (float)i33};
        return Ai;
    }
    return mat3_ident();
}

static inline void md_unitcell_inv_basis_extract(double out_Ai[3][3], const md_unitcell_t* cell) {
    if (cell) {
        if (!cell->flags) {
            MEMSET(out_Ai, 0, sizeof(double) * 9);
            return;
        }
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
        mat3_t G = {(float)g11, (float)g12, (float)g13, (float)g12, (float)g22, (float)g23, (float)g13, (float)g23, (float)g33};
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
    md_atom_idx_t p;
    md_atom_idx_t o5;
    md_atom_idx_t c5;
    md_atom_idx_t c4;
    md_atom_idx_t c3;
    md_atom_idx_t o3;
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
static inline int32_t* md_index_range_beg(md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx < md_array_size(data->offsets) - 1) {
        ptr = data->indices + data->offsets[range_idx];
    }
    return ptr;
}

static inline int32_t* md_index_range_end(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx < md_array_size(data->offsets) - 1) {
        ptr = data->indices + data->offsets[range_idx + 1];
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
