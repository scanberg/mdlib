#include "md_util.h"

#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_log.h>

#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// Tables
enum {
    Unknown = 0,
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V,
    Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru,
    Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb,
    Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh,
    Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og, Num_Elements
};

static const char* element_symbols[] = {
    "Xx", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

static const char* element_names[] = {
    "Unknown",     "Hydrogen",     "Helium",       "Lithium",     "Beryllium",   "Boron",         "Carbon",     "Nitrogen",   "Oxygen",
    "Fluorine",    "Neon",         "Sodium",       "Magnesium",   "Aluminium",   "Silicon",       "Phosphorus", "Sulfur",     "Chlorine",
    "Argon",       "Potassium",    "Calcium",      "Scandium",    "Titanium",    "Vanadium",      "Chromium",   "Manganese",  "Iron",
    "Cobalt",      "Nickel",       "Copper",       "Zinc",        "Gallium",     "Germanium",     "Arsenic",    "Selenium",   "Bromine",
    "Krypton",     "Rubidium",     "Strontium",    "Yttrium",     "Zirconium",   "Niobium",       "Molybdenum", "Technetium", "Ruthenium",
    "Rhodium",     "Palladium",    "Silver",       "Cadmium",     "Indium",      "Tin",           "Antimony",   "Tellurium",  "Iodine",
    "Xenon",       "Caesium",      "Barium",       "Lanthanum",   "Cerium",      "Praseodymium",  "Neodymium",  "Promethium", "Samarium",
    "Europium",    "Gadolinium",   "Terbium",      "Dysprosium",  "Holmium",     "Erbium",        "Thulium",    "Ytterbium",  "Lutetium",
    "Hafnium",     "Tantalum",     "Tungsten",     "Rhenium",     "Osmium",      "Iridium",       "Platinum",   "Gold",       "Mercury",
    "Thallium",    "Lead",         "Bismuth",      "Polonium",    "Astatine",    "Radon",         "Francium",   "Radium",     "Actinium",
    "Thorium",     "Protactinium", "Uranium",      "Neptunium",   "Plutonium",   "Americium",     "Curium",     "Berkelium",  "Californium",
    "Einsteinium", "Fermium",      "Mendelevium",  "Nobelium",    "Lawrencium",  "Rutherfordium", "Dubnium",    "Seaborgium", "Bohrium",
    "Hassium",     "Meitnerium",   "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium",      "Flerovium",  "Moscovium",  "Livermorium",
    "Tennessine",  "Oganesson"};

// http://dx.doi.org/10.1039/b801115j
static float element_covalent_radii[] = {
    0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53,
    1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45,
    1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
    1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80,
    1.69, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60};

// https://dx.doi.org/10.1021/jp8111556
static float element_vdw_radii[] = {
    1.00, 1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2.30, 2.15, 2.05,
    2.05, 2.05, 2.05, 2.00, 2.00, 2.00, 2.10, 1.87, 2.11, 1.85, 1.90, 1.83, 2.02, 3.03, 2.49, 2.40, 2.30, 2.15, 2.10, 2.05, 2.05, 2.00, 2.05, 2.10,
    2.20, 2.20, 1.93, 2.17, 2.06, 1.98, 2.16, 3.43, 2.68, 2.50, 2.48, 2.47, 2.45, 2.43, 2.42, 2.40, 2.38, 2.37, 2.35, 2.33, 2.32, 2.30, 2.28, 2.27,
    2.25, 2.20, 2.10, 2.05, 2.00, 2.00, 2.05, 2.10, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, 2.00, 2.30, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00};

// http://chemistry.wikia.com/wiki/List_of_elements_by_atomic_mass
// last element padded with 295
static float element_atomic_mass[] = {
    0,       1.00794,    4.002602, 6.941,     9.012182, 10.811,    12.0107,   14.0067,     15.9994, 18.9984032, 20.1797, 22.98976928,
    24.3050, 26.9815386, 28.0855,  30.973762, 32.065,   35.453,    39.0983,   39.948,      40.078,  44.955912,  47.867,  50.9415,
    51.9961, 54.938045,  55.845,   58.933195, 58.6934,  63.546,    65.409,    69.723,      72.64,   74.92160,   78.96,   79.904,
    83.798,  85.4678,    87.62,    88.90585,  91.224,   92.906,    95.94,     98,          101.07,  102.905,    106.42,  107.8682,
    112.411, 114.818,    118.710,  121.760,   127.60,   126.904,   131.293,   132.9054519, 137.327, 138.90547,  140.116, 140.90765,
    144.242, 145,        150.36,   151.964,   157.25,   158.92535, 162.500,   164.930,     167.259, 168.93421,  173.04,  174.967,
    178.49,  180.94788,  183.84,   186.207,   190.23,   192.217,   195.084,   196.966569,  200.59,  204.3833,   207.2,   208.98040,
    210,     210,        220,      223,       226,      227,       231.03588, 232.03806,   237,     238.02891,  243,     244,
    247,     247,        251,      252,       257,      258,       259,       262,         261,     262,        266,     264,
    277,     268,        271,      272,       285,      284,       289,       288,         292,     294,        295};

// http://jmol.sourceforge.net/jscolors/
static uint32_t element_cpk_colors[] = {
    0xFFFF00FF, 0xFFFFFFFF, 0xFFFFFFD9, 0xFF2222B2, 0xFF00FFC2, 0xFFB5B5FF, 0xFFB0B0B0, 0xFFFF8F8F, 0xFF0000F0, 0xFF50E090, 0xFFF5E3B3, 0xFFF25CAB,
    0xFF00FF8A, 0xFF908080, 0xFFA0C8F0, 0xFF00A5FF, 0xFF32C8FF, 0xFF1FF01F, 0xFFE3D180, 0xFFD4408F, 0xFF908080, 0xFFE6E6E6, 0xFF908080, 0xFFABA6A6,
    0xFF908080, 0xFF908080, 0xFF00A5FF, 0xFFA090F0, 0xFF2A2AA5, 0xFF2A2AA5, 0xFF2A2AA5, 0xFF8F8FC2, 0xFF8F8F66, 0xFFE380BD, 0xFF00A1FF, 0xFF2A2AA5,
    0xFFD1B85C, 0xFFB02E70, 0xFF00FF00, 0xFFFFFF94, 0xFFE0E094, 0xFFC9C273, 0xFFB5B554, 0xFF9E9E3B, 0xFF8F8F24, 0xFF8C7D0A, 0xFF856900, 0xFF908080,
    0xFF8FD9FF, 0xFF7375A6, 0xFF808066, 0xFFB5639E, 0xFF007AD4, 0xFF940094, 0xFFB09E42, 0xFF8F1757, 0xFF00A5FF, 0xFFFFD470, 0xFFC7FFFF, 0xFFC7FFD9,
    0xFFC7FFC7, 0xFFC7FFA3, 0xFFC7FF8F, 0xFFC7FF61, 0xFFC7FF45, 0xFFC7FF30, 0xFFC7FF1F, 0xFF9CFF00, 0xFF75E600, 0xFF52D400, 0xFF38BF00, 0xFF24AB00,
    0xFFFFC24D, 0xFFFFA64D, 0xFFD69421, 0xFFAB7D26, 0xFF966626, 0xFF875417, 0xFFE0D0D0, 0xFF23D1FF, 0xFFD0B8B8, 0xFF4D54A6, 0xFF615957, 0xFFB54F9E,
    0xFF005CAB, 0xFF454F75, 0xFF968242, 0xFF660042, 0xFF007D00, 0xFFFAAB70, 0xFFFFBA00, 0xFFFFA100, 0xFFFF8F00, 0xFFFF8000, 0xFFFF6B00, 0xFFF25C54,
    0xFFE35C78, 0xFFE34F8A, 0xFFD436A1, 0xFFD41FB3, 0xFFBA1FB3, 0xFFA60DB3, 0xFF870DBD, 0xFF6600C7, 0xFF5900CC, 0xFF4F00D1, 0xFF4500D9, 0xFF3800E0,
    0xFF2E00E6, 0xFF2600EB, 0xFF2200F0, 0xFF2000F6, 0xFF1E00F8, 0xFF1C00FA, 0xFF1A00FC, 0xFF1800FD, 0xFF1600FE, 0xFF1400FF, 0xFF1200FF};

static const char* amino_acids[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE",
};

static const char* dna[] = {"DA", "DA3", "DA5", "DC", "DC3", "DC5", "DG", "DG3", "DG5", "DT", "DT3", "DT5"};
static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

md_element md_util_lookup_element(str_t str) {
    if (str.len == 1 || str.len == 2) {
        for (uint8_t i = 0; i < ARRAY_SIZE(element_symbols); ++i) {
            const char* sym = element_symbols[i];
            if ((str.ptr[0] != sym[0]) || (str.len == 2 && str.ptr[1] != sym[1])) continue;
            return i;
        }
    }
    return 0;
}

md_element md_util_decode_element(str_t atom_name, str_t res_name) {
    if (!atom_name.ptr || !atom_name.len) return 0;
    if (!res_name.ptr  || !res_name.len) return 0;

    const char* beg = atom_name.ptr;
    const char* end = atom_name.ptr + atom_name.len;

    // TRIM
    while (beg < end && *beg && (is_digit(*beg) || is_whitespace(*beg))) ++beg;
    while (beg < end && (is_digit(end[-1]) || is_whitespace(*beg))) --end;

    // @TODO: Implement

    return 0;
}

const char* md_util_symbol(md_element element) {
    ASSERT(element < Num_Elements);
    return element_symbols[element];
}

const char* md_util_name(md_element element) {
    ASSERT(element < Num_Elements);
    return element_names[element];
}

float md_util_vdw_radius(md_element element) {
    ASSERT(element < Num_Elements);
    return element_vdw_radii[element];
}

float md_util_covalent_radius(md_element element) {
    ASSERT(element < Num_Elements);
    return element_covalent_radii[element];
}

float md_util_atomic_mass(md_element element) {
    ASSERT(element < Num_Elements);
    return element_atomic_mass[element];
}

uint32_t md_util_cpk_color(md_element element) {
    ASSERT(element < Num_Elements);
    return element_cpk_colors[element];
}


static inline float distance(const float v[3], const float u[3]) {
    const float d[3] = {v[0] - u[0], v[1] - u[1], v[2] - u[2]};
    return sqrtf(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

static inline bool cmp1(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == '\0';
}

static inline bool cmp2(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == ref[1] && str[2] == '\0';
}

static inline bool md_util_extract_backbone_atoms(md_backbone_atoms* backbone_atoms, const char** atom_names, md_range atom_range) {
    uint32_t bits = 0;
    md_backbone_atoms bb = {0};
    for (uint32_t i = atom_range.beg; i < atom_range.end; ++i) {
        if (!(bits & 1) && cmp1(atom_names[i], "N"))  { bb.n  = i; bits |= 1;  continue; }
        if (!(bits & 2) && cmp2(atom_names[i], "CA")) { bb.ca = i; bits |= 2;  continue; }
        if (!(bits & 4) && cmp1(atom_names[i], "C"))  { bb.c  = i; bits |= 4;  continue; }
        if (!(bits & 8) && cmp1(atom_names[i], "O"))  { bb.o  = i; bits |= 8;  continue; }
        if (!(bits & 8) && i == (atom_range.end - 1) && atom_names[i][0] == 'O') {
            bb.o = i; bits |= 8; continue;
        }
    }

    // If we have CA, C and O, we have enough for computing the backbone
    if (bits & (2|4|8)) {
        if (backbone_atoms) *backbone_atoms = bb;
        return true;
    }
    return false;
}


bool md_util_extract_backbone(md_backbone_atoms* backbone_atoms, const md_util_backbone_args_t* args) {
    ASSERT(backbone_atoms);
    ASSERT(args);
    memset(backbone_atoms, 0, args->residue.count * sizeof(md_backbone_atoms));
    bool result = true;
    for (uint32_t i = 0; i < args->residue.count; ++i) {
        if (!md_util_extract_backbone_atoms(&backbone_atoms[i], args->atom.name, args->residue.atom_range[i])) {
            md_print(MD_LOG_TYPE_INFO, "Could not extract backbone of residue[%i], possible that the residue is not an amino acid");
            result = false;
        }
    }
    return result;
}

static inline bool zhang_skolnick_ss(const md_util_secondary_structure_args_t* args, md_range res_range, int i, const float distances[3], float delta) {
    for (int j = MAX((int)res_range.beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= (int)res_range.end) continue;
            const int ca_j = args->residue.backbone_atoms[j].ca;
            const int ca_k = args->residue.backbone_atoms[j + k].ca;
            const float pos_j[3] = {args->atom.x[ca_j], args->atom.y[ca_j], args->atom.z[ca_j]};
            const float pos_k[3] = {args->atom.x[ca_k], args->atom.y[ca_k], args->atom.z[ca_k]};
            const float d = distance(pos_j, pos_k);
            if (ABS(d - distances[k - 2]) > delta) {
                return false;
            }
        }
    }
    return true;
}

static inline bool is_helical(const md_util_secondary_structure_args_t* args, md_range res_range, int i) {
    const float distances[] = { 5.45f, 5.18f, 6.37f };
    const float delta = 2.1f;
    return zhang_skolnick_ss(args, res_range, i, distances, delta);
}

static inline bool is_sheet(const md_util_secondary_structure_args_t* args, md_range res_range, int i) {
    const float distances[] = { 6.1f, 10.4f, 13.0f };
    const float delta = 1.42f;
    return zhang_skolnick_ss(args, res_range, i, distances, delta);
}

// TM-align: a protein structure alignment algorithm based on the Tm-score
// doi:10.1093/nar/gki524
bool md_util_compute_secondary_structure(md_secondary_structure* secondary_structure, const md_util_secondary_structure_args_t* args) {
    ASSERT(args);
    ASSERT(args->atom.x);
    ASSERT(args->atom.y);
    ASSERT(args->atom.z);
    ASSERT(args->residue.backbone_atoms);
    ASSERT(args->chain.residue_range);

    for (uint32_t chain_idx = 0; chain_idx < args->chain.count; ++chain_idx) {
        const md_range range = args->chain.residue_range[chain_idx];
        if (range.end - range.beg < 4) continue;

        // Classify residues
        for (uint32_t i = range.beg; i < range.end; ++i) {
            md_secondary_structure ss = MD_SECONDARY_STRUCTURE_COIL;
            if (is_sheet(args, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_SHEET;
            }
            else if (is_helical(args, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_HELIX;
            }

            secondary_structure[i] = ss;
        }

        // Set squished isolated structures to the surrounding (only for sheets and helices)
        md_secondary_structure* ss = secondary_structure;
        for (uint32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i-1] != MD_SECONDARY_STRUCTURE_COIL && ss[i] != ss[i-1] && ss[i-1] == ss[i+1]) ss[i] = ss[i-1];
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MD_SECONDARY_STRUCTURE_COIL;
        for (uint32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i] != ss[i-1] && ss[i] != ss[i+1]) ss[i] = MD_SECONDARY_STRUCTURE_COIL;
        }
        if (ss[range.end - 1] != ss[range.end - 2]) ss[range.end - 1] = MD_SECONDARY_STRUCTURE_COIL;
    }

    return true;
}

static inline bool covelent_bond_heuristic(float x0, float y0, float z0, md_element e0, float x1, float y1, float z1, md_element e1) {
    ASSERT(e0 < Num_Elements);
    ASSERT(e1 < Num_Elements);
    const float d = element_covalent_radii[e0] + element_covalent_radii[e1];
    const float d_min = d - 0.5f;
    const float d_max = d + 0.3f;
    const float d2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    return (d_min * d_min) < d2 && d2 < (d_max * d_max);
}

md_bond* md_util_extract_covalent_bonds(const md_util_covalent_bond_args_t* args, struct md_allocator_i* alloc) {
    md_bond* bonds = 0;

    for (int64_t ri = 0; ri < args->residue.count; ++ri) {
        const int64_t pre_internal_bond_count = md_array_size(bonds);

        // Compute internal bonds (within residue)
        for (int64_t i = args->residue.atom_range[ri].beg; i < args->residue.atom_range[ri].end - 1; ++i) {
            for (int64_t j = i + 1; j < args->residue.atom_range[ri].end; ++j) {
                if (covelent_bond_heuristic(args->atom.x[i], args->atom.y[i], args->atom.z[i], args->atom.element[i],
                                            args->atom.x[j], args->atom.y[j], args->atom.z[j], args->atom.element[j])) {
                    md_bond bond = {i, j};
                    md_array_push(bonds, bond, alloc);
                }
            }
        }

        const int64_t post_internal_bond_count = md_array_size(bonds);

        if (args->residue.internal_bond_range) {
            args->residue.internal_bond_range[ri].beg = pre_internal_bond_count;
            args->residue.internal_bond_range[ri].end = post_internal_bond_count;
        }

        const int64_t rj = ri + 1;
        if (rj < args->residue.count) {
            // Compute extarnal bonds (between residues)
            for (int64_t i = args->residue.atom_range[ri].beg; i < args->residue.atom_range[ri].end; ++i) {
                for (int64_t j = args->residue.atom_range[rj].beg; j < args->residue.atom_range[rj].end; ++j) {
                    if (covelent_bond_heuristic(args->atom.x[i], args->atom.y[i], args->atom.z[i], args->atom.element[i],
                        args->atom.x[j], args->atom.y[j], args->atom.z[j], args->atom.element[j])) {
                        md_bond bond = {i, j};
                        md_array_push(bonds, bond, alloc);
                    }
                }
            }

            const int64_t post_external_bond_count = md_array_size(bonds);

            if (args->residue.complete_bond_range) {
                args->residue.complete_bond_range[ri].end = post_external_bond_count;
                args->residue.complete_bond_range[rj].beg = post_internal_bond_count;
            }
        }
    }

    if (args->residue.complete_bond_range) {
        args->residue.complete_bond_range[0].beg = 0;
        md_array_last(args->residue.complete_bond_range)->end = md_array_size(bonds);
    }

    return bonds;
}

#ifdef __cplusplus
}
#endif