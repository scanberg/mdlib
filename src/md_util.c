#include "md_util.h"

#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_ring_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_simd.h>
#include <core/md_spatial_hash.h>
#include <core/md_bitfield.h>

#include <md_trajectory.h>
#include <md_molecule.h>

#include <math.h>
#include <string.h>
#include <float.h>

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

#define bake(str) {str, sizeof(str)-1}

static const str_t element_symbols[] = {
    bake("Xx"), bake("H"),  bake("He"), bake("Li"), bake("Be"), bake("B"),  bake("C"),  bake("N"),  bake("O"),  bake("F"),  bake("Ne"), bake("Na"), bake("Mg"), bake("Al"), bake("Si"), bake("P"),  bake("S"),  bake("Cl"), bake("Ar"), bake("K"),  bake("Ca"), bake("Sc"), bake("Ti"), bake("V"),
    bake("Cr"), bake("Mn"), bake("Fe"), bake("Co"), bake("Ni"), bake("Cu"), bake("Zn"), bake("Ga"), bake("Ge"), bake("As"), bake("Se"), bake("Br"), bake("Kr"), bake("Rb"), bake("Sr"), bake("Y"),  bake("Zr"), bake("Nb"), bake("Mo"), bake("Tc"), bake("Ru"), bake("Rh"), bake("Pd"), bake("Ag"),
    bake("Cd"), bake("In"), bake("Sn"), bake("Sb"), bake("Te"), bake("I"),  bake("Xe"), bake("Cs"), bake("Ba"), bake("La"), bake("Ce"), bake("Pr"), bake("Nd"), bake("Pm"), bake("Sm"), bake("Eu"), bake("Gd"), bake("Tb"), bake("Dy"), bake("Ho"), bake("Er"), bake("Tm"), bake("Yb"), bake("Lu"),
    bake("Hf"), bake("Ta"), bake("W"),  bake("Re"), bake("Os"), bake("Ir"), bake("Pt"), bake("Au"), bake("Hg"), bake("Tl"), bake("Pb"), bake("Bi"), bake("Po"), bake("At"), bake("Rn"), bake("Fr"), bake("Ra"), bake("Ac"), bake("Th"), bake("Pa"), bake("U"),  bake("Np"), bake("Pu"), bake("Am"),
    bake("Cm"), bake("Bk"), bake("Cf"), bake("Es"), bake("Fm"), bake("Md"), bake("No"), bake("Lr"), bake("Rf"), bake("Db"), bake("Sg"), bake("Bh"), bake("Hs"), bake("Mt"), bake("Ds"), bake("Rg"), bake("Cn"), bake("Nh"), bake("Fl"), bake("Mc"), bake("Lv"), bake("Ts"), bake("Og"), };


static const str_t element_names[] = {
    bake("Unknown"),     bake("Hydrogen"),     bake("Helium"),       bake("Lithium"),     bake("Beryllium"),   bake("Boron"),         bake("Carbon"),     bake("Nitrogen"),   bake("Oxygen"),
    bake("Fluorine"),    bake("Neon"),         bake("Sodium"),       bake("Magnesium"),   bake("Aluminium"),   bake("Silicon"),       bake("Phosphorus"), bake("Sulfur"),     bake("Chlorine"),
    bake("Argon"),       bake("Potassium"),    bake("Calcium"),      bake("Scandium"),    bake("Titanium"),    bake("Vanadium"),      bake("Chromium"),   bake("Manganese"),  bake("Iron"),
    bake("Cobalt"),      bake("Nickel"),       bake("Copper"),       bake("Zinc"),        bake("Gallium"),     bake("Germanium"),     bake("Arsenic"),    bake("Selenium"),   bake("Bromine"),
    bake("Krypton"),     bake("Rubidium"),     bake("Strontium"),    bake("Yttrium"),     bake("Zirconium"),   bake("Niobium"),       bake("Molybdenum"), bake("Technetium"), bake("Ruthenium"),
    bake("Rhodium"),     bake("Palladium"),    bake("Silver"),       bake("Cadmium"),     bake("Indium"),      bake("Tin"),           bake("Antimony"),   bake("Tellurium"),  bake("Iodine"),
    bake("Xenon"),       bake("Caesium"),      bake("Barium"),       bake("Lanthanum"),   bake("Cerium"),      bake("Praseodymium"),  bake("Neodymium"),  bake("Promethium"), bake("Samarium"),
    bake("Europium"),    bake("Gadolinium"),   bake("Terbium"),      bake("Dysprosium"),  bake("Holmium"),     bake("Erbium"),        bake("Thulium"),    bake("Ytterbium"),  bake("Lutetium"),
    bake("Hafnium"),     bake("Tantalum"),     bake("Tungsten"),     bake("Rhenium"),     bake("Osmium"),      bake("Iridium"),       bake("Platinum"),   bake("Gold"),       bake("Mercury"),
    bake("Thallium"),    bake("Lead"),         bake("Bismuth"),      bake("Polonium"),    bake("Astatine"),    bake("Radon"),         bake("Francium"),   bake("Radium"),     bake("Actinium"),
    bake("Thorium"),     bake("Protactinium"), bake("Uranium"),      bake("Neptunium"),   bake("Plutonium"),   bake("Americium"),     bake("Curium"),     bake("Berkelium"),  bake("Californium"),
    bake("Einsteinium"), bake("Fermium"),      bake("Mendelevium"),  bake("Nobelium"),    bake("Lawrencium"),  bake("Rutherfordium"), bake("Dubnium"),    bake("Seaborgium"), bake("Bohrium"),
    bake("Hassium"),     bake("Meitnerium"),   bake("Darmstadtium"), bake("Roentgenium"), bake("Copernicium"), bake("Nihonium"),      bake("Flerovium"),  bake("Moscovium"),  bake("Livermorium"),
    bake("Tennessine"),  bake("Oganesson"),                          
};

#undef bake

// http://dx.doi.org/10.1039/b801115j
static float element_covalent_radii[] = {
    0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53,
    1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45,
    1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
    1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80,
    1.69, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60
};

// https://dx.doi.org/10.1021/jp8111556
static float element_vdw_radii[] = {
    1.00, 1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2.30, 2.15, 2.05,
    2.05, 2.05, 2.05, 2.00, 2.00, 2.00, 2.10, 1.87, 2.11, 1.85, 1.90, 1.83, 2.02, 3.03, 2.49, 2.40, 2.30, 2.15, 2.10, 2.05, 2.05, 2.00, 2.05, 2.10,
    2.20, 2.20, 1.93, 2.17, 2.06, 1.98, 2.16, 3.43, 2.68, 2.50, 2.48, 2.47, 2.45, 2.43, 2.42, 2.40, 2.38, 2.37, 2.35, 2.33, 2.32, 2.30, 2.28, 2.27,
    2.25, 2.20, 2.10, 2.05, 2.00, 2.00, 2.05, 2.10, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, 2.00, 2.30, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00
};

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
    277,     268,        271,      272,       285,      284,       289,       288,         292,     294,        295
};

// https://en.wikipedia.org/wiki/Valence_(chemistry)
// Elements with unknown max valence is assigned with 0
static uint8_t element_max_valence[] = {
    0,
    1, 0,
    1, 2, 3, 4, 5, 2, 1, 0,
    1, 2, 3, 4, 5, 6, 7, 0,
    1, 2, 3, 4, 5, 6, 7, 7, 5, 4, 4, 2, 3, 4, 5, 6, 7, 2,
    1, 2, 3, 4, 5, 6, 7, 8, 6, 4, 3, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 4, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 6, 5, 2, 3, 4, 5, 6, 7, 6,
    1, 2, 3, 4, 5, 6, 7, 7, 7, 6, 5, 5, 4, 3, 3, 3, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

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
    0xFF2E00E6, 0xFF2600EB, 0xFF2200F0, 0xFF2000F6, 0xFF1E00F8, 0xFF1C00FA, 0xFF1A00FC, 0xFF1800FD, 0xFF1600FE, 0xFF1400FF, 0xFF1200FF
};

static const char* amino_acids[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE", "HISE"
};

static const char* dna[] = {"DA", "DA3", "DA5", "DC", "DC3", "DC5", "DG", "DG3", "DG5", "DT", "DT3", "DT5"};
static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

static inline int64_t match_str_in_array(str_t str, const char* arr[], int64_t arr_len) {
    for (int64_t i = 0; i < arr_len; ++i) {
        if (str_equal_cstr(str, arr[i])) return i;
    }
    return -1;
}

static inline md_simd_f32_t simd_deperiodize(md_simd_f32_t x, md_simd_f32_t r, md_simd_f32_t period) {
    md_simd_f32_t d = md_simd_sub(x, r);
    md_simd_f32_t dx = md_simd_div(d, period);
    dx = md_simd_sub(dx, md_simd_round(dx));
    return md_simd_add(r, md_simd_mul(dx, period));
}
    
bool md_util_resname_dna(str_t str) {
    return match_str_in_array(str, dna, ARRAY_SIZE(dna)) != -1;
}
    
bool md_util_resname_acidic(str_t str) {
    return match_str_in_array(str, acidic, ARRAY_SIZE(acidic)) != -1;
}

bool md_util_resname_basic(str_t str) {
    return match_str_in_array(str, basic, ARRAY_SIZE(basic)) != -1;
}
    
bool md_util_resname_neutral(str_t str) {
    return match_str_in_array(str, neutral, ARRAY_SIZE(neutral)) != -1;
}
    
bool md_util_resname_water(str_t str) {
    return match_str_in_array(str, water, ARRAY_SIZE(water)) != -1;
}
    
bool md_util_resname_hydrophobic(str_t str) {
    return match_str_in_array(str, hydrophobic, ARRAY_SIZE(hydrophobic)) != -1;
}
    
bool md_util_resname_amino_acid(str_t str) {
    return match_str_in_array(str, amino_acids, ARRAY_SIZE(amino_acids)) != -1;
}

md_element_t md_util_element_lookup(str_t str) {
    if (str.len == 1 || str.len == 2) {
        for (md_element_t i = 0; i < ARRAY_SIZE(element_symbols); ++i) {
            str_t sym = element_symbols[i];
            if (str_equal(str, sym)) return i;
        }
    }
    return 0;
}

md_element_t md_util_element_lookup_ignore_case(str_t str) {
    if (str.len == 1 || str.len == 2) {
        for (md_element_t i = 0; i < ARRAY_SIZE(element_symbols); ++i) {
            str_t sym = element_symbols[i];
            if (str.len == sym.len) {
                if (str.len == 1) {
                    if (to_upper(str.ptr[0]) == to_upper(sym.ptr[0])) return i;
                } else {
                    if (to_upper(str.ptr[0]) == to_upper(sym.ptr[0]) && to_upper(str.ptr[1]) == to_upper(sym.ptr[1])) return i;
                }
            }
        }
    }
    return 0;
}

// Trim whitespace, digits and 'X's
str_t trim_type(str_t type) {
    const char* beg = str_beg(type);
    const char* end = str_end(type);
    const char* c = beg;
    while (c < end && *c && (is_digit(*c) || is_whitespace(*c) || *c == 'x' || *c == 'X')) ++c;
    beg = c;

    c = beg;
    while (c < end && is_alpha(*c) && !(*c == 'x' || *c == 'X')) ++c;
    end = c;

    return (str_t) { .ptr = beg, .len = end-beg };
}

static bool amino_acid_heuristic(const md_label_t labels[], int size) {
    // This is the minimal set of types which needs to be present in the case of Glycine and excluding hydrogen
#define BIT_N  1
#define BIT_CA 2
#define BIT_C  4
#define BIT_O  8

    int count = 0;
    int bits  = 0;
    for (int i = 0; i < size; ++i) {
        str_t lbl = { labels[i].buf, labels[i].len };
        lbl = trim_type(lbl);
        if (lbl.len && lbl.ptr[0] != 'H') {
            if (str_equal_cstr(lbl, "N")) bits |= BIT_N;
            else if (str_equal_cstr(lbl, "CA")) bits |= BIT_CA;
            else if (str_equal_cstr(lbl, "C")) bits |= BIT_C;
            else if (str_equal_cstr(lbl, "O")) bits |= BIT_O;
            count += 1;
        }
    }
    return 4 <= count && count <= 14 && bits == (BIT_N | BIT_CA | BIT_C | BIT_O);
#undef BIT_N
#undef BIT_CA
#undef BIT_C
#undef BIT_O
}

bool md_util_element_decode(md_element_t element[], int64_t capacity, const struct md_molecule_t* mol) {
    ASSERT(capacity >= 0);
    ASSERT(mol);
    ASSERT(mol->atom.count >= 0);

    // @TODO, PERF: Iterate over residues and check if the entire residue is an amino acid
    const int64_t count = MIN(capacity, mol->atom.count);
    for (int64_t i = 0; i < count; ++i) {
        if (element[i] != 0) continue;

        str_t original = {mol->atom.name[i].buf, mol->atom.name[i].len};

        // Trim whitespace, digits and 'X's
        str_t name = trim_type(original);

        if (name.len > 0) {

            md_element_t elem = 0;
            if ((elem = md_util_element_lookup(name)) != 0) goto done;

            // If amino acid, try to deduce the element from that
            if (mol->atom.residue_idx && mol->residue.name) {
                str_t resname = str_trim(LBL_TO_STR(mol->residue.name[mol->atom.residue_idx[i]]));
                md_range_t res_range = mol->residue.atom_range[mol->atom.residue_idx[i]];
                if (md_util_resname_amino_acid(resname) ||
                    amino_acid_heuristic(mol->atom.name + res_range.beg, res_range.end - res_range.beg))
                {
                    // EASY-PEASY, we just try to match against the first character
                    name.len = 1;
                    elem = md_util_element_lookup_ignore_case(name);
                    goto done;
                }
            }

            // Heuristic cases

            // CA -> Carbon, not Calcium (if part of a residue > 2 atoms)


            // 2 letters + 1 digit (e.g. HO[0-9]) usually means just look at the first letter
            if (original.len == 3 && is_digit(original.ptr[2])) {
                elem = md_util_element_lookup_ignore_case(str_substr(original, 0, 1));
                goto done;
            }

            // Try to match against several characters but ignore the case
            if (name.len > 1) {
                name.len = 2;
                elem = md_util_element_lookup_ignore_case(name);
            }

            // Last resort, try to match against single first character
            if (elem == 0) {
                name.len = 1;
                elem = md_util_element_lookup_ignore_case(name);
            }

        done:
            element[i] = elem;
        }
    }

    return true;
}

bool md_util_element_from_mass(md_element_t element[], const float mass[], int64_t count) {
    if (!element) {
        MD_LOG_ERROR("element is null");
        return false;
    }

    if (!mass) {
        MD_LOG_ERROR("element is null");
        return false;
    }

    const float eps = 1.0e-2;
    for (int64_t i = 0; i < count; ++i) {
        md_element_t elem = 0;
        const float m = mass[i];
        for (uint8_t j = 1; j < (uint8_t)ARRAY_SIZE(element_atomic_mass); ++j) {
            if (fabs(m - element_atomic_mass[j]) < eps) {
                elem = j;
                break;
            }
        }

        element[i] = elem;
    }

    return false;
}

str_t md_util_element_symbol(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return element_symbols[element];
}

str_t md_util_element_name(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return element_names[element];
}

float md_util_element_vdw_radius(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return element_vdw_radii[element];
}

float md_util_element_covalent_radius(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return element_covalent_radii[element];
}

float md_util_element_atomic_mass(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return element_atomic_mass[element];
}

int md_util_element_max_valence(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return (int)element_max_valence[element];
}

uint32_t md_util_element_cpk_color(md_element_t element) {
    element = element < Num_Elements ? element : 0;
    return element_cpk_colors[element];
}

static inline bool cmp1(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == '\0';
}

static inline bool cmp2(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == ref[1] && str[2] == '\0';
}

static inline bool extract_backbone_atoms(md_backbone_atoms_t* backbone_atoms, const md_label_t* atom_names, md_range_t atom_range) {
    uint32_t bits = 0;
    md_backbone_atoms_t bb = {0};
    for (int32_t i = atom_range.beg; i < atom_range.end; ++i) {
        if (!(bits & 1) && cmp1(atom_names[i].buf, "N"))  { bb.n  = i; bits |= 1;  continue; }
        if (!(bits & 2) && cmp2(atom_names[i].buf, "CA")) { bb.ca = i; bits |= 2;  continue; }
        if (!(bits & 4) && cmp1(atom_names[i].buf, "C"))  { bb.c  = i; bits |= 4;  continue; }
        if (!(bits & 8) && cmp1(atom_names[i].buf, "O"))  { bb.o  = i; bits |= 8;  continue; }
        if (!(bits & 8) && i == (atom_range.end - 1) && atom_names[i].buf[0] == 'O') {
            bb.o = i; bits |= 8; continue;
        }
    }

    // If we have CA, C and O, we have enough for computing the backbone
    if ((bits & (2|4|8)) == (2|4|8)) {
        if (backbone_atoms) *backbone_atoms = bb;
        return true;
    }
    return false;
}

bool md_util_backbone_atoms_extract_from_residue_idx(md_backbone_atoms_t* backbone_atoms, md_residue_idx_t res_idx, const md_molecule_t* mol) {
    ASSERT(backbone_atoms);
    ASSERT(mol);
    if (res_idx < 0 || mol->residue.count <= res_idx) return false;
    return extract_backbone_atoms(backbone_atoms, mol->atom.name, mol->residue.atom_range[res_idx]);
}

static inline bool zhang_skolnick_ss(const md_molecule_t* args, md_range_t bb_range, int i, const float distances[3], float delta) {
    for (int j = MAX((int)bb_range.beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= (int)bb_range.end) continue;
            const int ca_j = args->backbone.atoms[j].ca;
            const int ca_k = args->backbone.atoms[j + k].ca;
            const vec3_t pos_j = {args->atom.x[ca_j], args->atom.y[ca_j], args->atom.z[ca_j]};
            const vec3_t pos_k = {args->atom.x[ca_k], args->atom.y[ca_k], args->atom.z[ca_k]};
            const float d = vec3_distance(pos_j, pos_k);
            if (ABS(d - distances[k - 2]) > delta) {
                return false;
            }
        }
    }
    return true;
}

static inline bool is_helical(const md_molecule_t* mol, md_range_t bb_range, int i) {
    const float distances[] = { 5.45f, 5.18f, 6.37f };
    const float delta = 2.1f;
    return zhang_skolnick_ss(mol, bb_range, i, distances, delta);
}

static inline bool is_sheet(const md_molecule_t* mol, md_range_t bb_range, int i) {
    const float distances[] = { 6.1f, 10.4f, 13.0f };
    const float delta = 1.42f;
    return zhang_skolnick_ss(mol, bb_range, i, distances, delta);
}

// TM-align: a protein structure alignment algorithm based on the Tm-score
// doi:10.1093/nar/gki524
bool md_util_backbone_secondary_structure_compute(md_secondary_structure_t secondary_structure[], int64_t capacity, const struct md_molecule_t* mol) {
    if (!secondary_structure) return false;
    if (capacity < 0) return false;

    MEMSET(secondary_structure, 0, capacity * sizeof(md_secondary_structure_t));

    if (!mol) return false;
    if (!mol->atom.x) return false;
    if (!mol->atom.y) return false;
    if (!mol->atom.z) return false;
    if (!mol->backbone.atoms) return false;
    if (!mol->backbone.range) return false;

    for (int64_t bb_idx = 0; bb_idx < mol->backbone.range_count; ++bb_idx) {
        const md_range_t range = mol->backbone.range[bb_idx];
        ASSERT(range.end <= capacity);

        if (range.end - range.beg < 4) {
            continue;
        }

        // Classify residues
        for (int32_t i = range.beg; i < range.end; ++i) {
            md_secondary_structure_t ss = MD_SECONDARY_STRUCTURE_COIL;

            if (is_sheet(mol, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_SHEET;
            }
            else if (is_helical(mol, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_HELIX;
            }
            secondary_structure[i] = ss;
        }

        // Set squished isolated structures to the surrounding (only for sheets and helices)
        md_secondary_structure_t* ss = secondary_structure;
        for (int64_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i-1] != MD_SECONDARY_STRUCTURE_COIL && ss[i] != ss[i-1] && ss[i-1] == ss[i+1]) ss[i] = ss[i-1];
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MD_SECONDARY_STRUCTURE_COIL;
        for (int64_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i] != ss[i-1] && ss[i] != ss[i+1]) ss[i] = MD_SECONDARY_STRUCTURE_COIL;
        }
        if (ss[range.end - 1] != ss[range.end - 2]) ss[range.end - 1] = MD_SECONDARY_STRUCTURE_COIL;
    }

    return true;
}

static inline float dihedral_angle(vec3_t p0, vec3_t p1, vec3_t p2, vec3_t p3) {
    const vec3_t b1 = vec3_normalize(vec3_sub(p1, p0));
    const vec3_t b2 = vec3_normalize(vec3_sub(p2, p1));
    const vec3_t b3 = vec3_normalize(vec3_sub(p3, p2));
    const vec3_t c1 = vec3_cross(b1, b2);
    const vec3_t c2 = vec3_cross(b2, b3);
    return atan2f(vec3_dot(vec3_cross(c1, c2), b2), vec3_dot(c1, c2));
}

bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], int64_t capacity, const md_molecule_t* mol) {
    if (!backbone_angles) return false;
    if (capacity < 0) return false;

    MEMSET(backbone_angles, 0, sizeof(md_backbone_angles_t) * capacity);

    if (!mol) return false;
    if (!mol->atom.x) return false;
    if (!mol->atom.y) return false;
    if (!mol->atom.z) return false;
    if (!mol->backbone.atoms) return false;

    for (int64_t bb_idx = 0; bb_idx < mol->backbone.range_count; ++bb_idx) {
        const md_range_t range = mol->backbone.range[bb_idx];
        ASSERT(range.end <= capacity);

        if (range.end - range.beg < 4) {
            continue;
        }

        for (int64_t i = range.beg + 1; i < range.end - 1; ++i) {
            const vec3_t c_prev = { mol->atom.x[mol->backbone.atoms[i-1].c], mol->atom.y[mol->backbone.atoms[i-1].c], mol->atom.z[mol->backbone.atoms[i-1].c] };
            const vec3_t n      = { mol->atom.x[mol->backbone.atoms[i].n]  , mol->atom.y[mol->backbone.atoms[i].n]  , mol->atom.z[mol->backbone.atoms[i].n]   };
            const vec3_t ca     = { mol->atom.x[mol->backbone.atoms[i].ca] , mol->atom.y[mol->backbone.atoms[i].ca] , mol->atom.z[mol->backbone.atoms[i].ca]  };
            const vec3_t c      = { mol->atom.x[mol->backbone.atoms[i].c]  , mol->atom.y[mol->backbone.atoms[i].c]  , mol->atom.z[mol->backbone.atoms[i].c]   };
            const vec3_t n_next = { mol->atom.x[mol->backbone.atoms[i+1].n], mol->atom.y[mol->backbone.atoms[i+1].n], mol->atom.z[mol->backbone.atoms[i+1].n] };
            backbone_angles[i].phi = dihedral_angle(c_prev, n, ca, c);
            backbone_angles[i].psi = dihedral_angle(n, ca, c, n_next);
        }
    }

    return true;
}

bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], int64_t capacity, const md_molecule_t* mol) {
    ASSERT(ramachandran_types);
    MEMSET(ramachandran_types, MD_RAMACHANDRAN_TYPE_UNKNOWN, sizeof(md_ramachandran_type_t) * capacity);

    if (capacity < 0) return false;
    if (mol->backbone.count == 0) return false;
    if (mol->residue.count == 0) return false;

    ASSERT(mol->residue.name);
    ASSERT(mol->backbone.residue_idx);

    int64_t size = MIN(capacity, mol->backbone.count);

    for (int64_t i = 0; i < size; ++i) {
        int64_t res_idx = mol->backbone.residue_idx[i];
        ASSERT(res_idx < mol->residue.count);

        str_t resname = LBL_TO_STR(mol->residue.name[res_idx]);
        if (str_equal_cstr_n(resname, "GLY", 3)) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GLYCINE;
        } else if (str_equal_cstr_n(resname, "PRO", 3)) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PROLINE;
            ramachandran_types[i - 1] = MD_RAMACHANDRAN_TYPE_PREPROL;
        } else {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GENERAL;
        }
    }

    return true;
}

static md_index_data_t compute_connectivity_data(int64_t atom_count, const md_bond_t bonds[], int64_t bond_count, md_allocator_i* alloc) {
    md_index_data_t data = { 0 };

    md_array_resize(data.ranges, atom_count, alloc);
    MEMSET(data.ranges, 0, md_array_bytes(data.ranges));

    // This have length of 2 * bond_count (2 == one for direction)
    md_array_resize(data.indices, 2 * bond_count, alloc);

    uint32_t* local_offset = 0;
    md_array_resize(local_offset, bond_count, alloc);

    // Two packed 16-bit local offsets for each of the bond idx
    // Accumulate the length of each range in the .end component
    for (int64_t i = 0; i < bond_count; ++i) {
        const uint32_t off0 = (uint32_t)data.ranges[bonds[i].idx[0]].end++;
        const uint32_t off1 = (uint32_t)data.ranges[bonds[i].idx[1]].end++;
        local_offset[i] = (off1 << 16) | off0;
    }

    // Compute complete edge ranges (compute beg and add to end in order to convert the stored length to a proper absolute offset)
    for (int64_t i = 1; i < atom_count; ++i) {
        data.ranges[i].beg = data.ranges[i - 1].end;
        data.ranges[i].end += data.ranges[i].beg;
    }

    // Write edge indices to correct location
    for (int64_t i = 0; i < bond_count; ++i) {
        const md_bond_t bond = bonds[i];
        const int atom_a = bond.idx[0];
        const int atom_b = bond.idx[1];
        const int local_a = (int)(local_offset[i] & 0xFFFF);
        const int local_b = (int)(local_offset[i] >> 16);
        const md_range_t range_a = data.ranges[atom_a];
        const md_range_t range_b = data.ranges[atom_b];

        const int idx_a = range_a.beg + local_a;
        const int idx_b = range_b.beg + local_b;
        ASSERT(idx_a < range_a.end);
        ASSERT(idx_b < range_b.end);

        ASSERT(idx_a < md_array_size(data.indices));
        ASSERT(idx_b < md_array_size(data.indices));

        // Store the cross references to the 'other' atom index signified by the bond in the correct location
        data.indices[idx_a] = atom_b;
        data.indices[idx_b] = atom_a;
    }

    md_array_free(local_offset, alloc);

    return data;
}

static inline bool covelent_bond_heuristic(float dist_squared, md_element_t elem_a, md_element_t elem_b) {
    const float d = element_covalent_radii[elem_a] + element_covalent_radii[elem_b];
    const float d_min = d - 0.5f;
    const float d_max = d + 0.3f;
    return (d_min * d_min) < dist_squared && dist_squared < (d_max * d_max);
}

bool md_util_compute_covalent_bonds(md_bond_data_t* bond_data, const float* x, const float* y, const float* z, const md_element_t* element, const md_residue_idx_t* res_idx, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    if (!bond_data) {
        MD_LOG_ERROR("missing parameter bond_data");
        return false;
    }

    if (!x) {
        MD_LOG_ERROR("missing parameter atom_x");
        return false;
    }

    if (!y) {
        MD_LOG_ERROR("missing parameter atom_y");
        return false;
    }

    if (!z) {
        MD_LOG_ERROR("missing parameter atom_z");
        return false;
    }

    if (!element) {
        MD_LOG_ERROR("missing parameter atom_element");
        return false;
    }

    const vec4_t pbc_ext4 = vec4_from_vec3(pbc_ext, 0);

    if (res_idx) {
        // atom residue indices are given,
        // First find connections first within the residue, then to the next residue

        md_range_t range = {0, 0};
        md_range_t prev_range = {0, 0};
        for (int i = 0; i < (int)count; ++i) {
            while (i < (int)count && res_idx[i] == res_idx[range.beg]) ++i;
            range.end = i;

            for (int j = prev_range.beg; j < prev_range.end; ++j) {
                for (int k = range.beg; k < range.end; ++k) {
                    const vec4_t a = {x[j], y[j], z[j], 0};
                    const vec4_t b = {x[k], y[k], z[k], 0};
                    const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext4);
                    if (covelent_bond_heuristic(d2, element[j], element[k])) {
                        md_array_push(bond_data->bond, ((md_bond_t){j, k}), alloc);
                    }
                }
            }

            for (int j = range.beg; j < range.end - 1; ++j) {
                for (int k = j + 1; k < range.end; ++k) {
                    const vec4_t a = {x[j], y[j], z[j], 0};
                    const vec4_t b = {x[k], y[k], z[k], 0};
                    const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext4);
                    if (covelent_bond_heuristic(d2, element[j], element[k])) {
                        md_array_push(bond_data->bond, ((md_bond_t){j, k}), alloc);
                    }
                }
            }

            prev_range = range;
            range.beg  = range.end;
        }
        
        bond_data->connectivity = compute_connectivity_data(count, bond_data->bond, md_array_size(bond_data->bond), alloc);
    }
    else {
        md_spatial_hash_t sh = {0};
        md_spatial_hash_init_soa(&sh, x, y, z, count, pbc_ext, default_allocator);

        const float cutoff = 3.0f;

        md_ring_allocator_t* ring_alloc = get_thread_ring_allocator();
        md_allocator_i temp_alloc = md_ring_allocator_create_interface(ring_alloc);
        uint64_t reset_pos = md_ring_allocator_get_pos(ring_alloc);

        for (int i = 0; i < (int)count; ++i) {
            const vec3_t pos = {x[i], y[i], z[i]};
            
            // Reset ring allocator pos
            md_ring_allocator_set_pos(ring_alloc, reset_pos);
            md_array(uint32_t) indices = md_spatial_hash_query_idx(&sh, pos, cutoff, &temp_alloc);
            md_array(int32_t) atom_connectivity = 0;
        
            const int64_t num_indices = md_array_size(indices);
            for (int64_t iter = 0; iter < num_indices; ++iter) {
                const int j = indices[iter];

                if (res_idx && abs(res_idx[j] - res_idx[i]) > 2) {
                    continue;
                }

                const vec4_t pos_a = vec4_from_vec3(pos, 0);
                const vec4_t pos_b = {x[j], y[j], z[j], 0};
                const float d2 = vec4_periodic_distance_squared(pos_a, pos_b, pbc_ext4);
                if (covelent_bond_heuristic(d2, element[i], element[j])) {
                    if (i < j) {
                        // Only store monotonic bonds connections
                        md_array_push(bond_data->bond, ((md_bond_t){i, j}), alloc);
                    }
                    // Store all connectivity to atom
                    md_array_push(atom_connectivity, j, &temp_alloc);
                }
            
                md_index_data_push(&bond_data->connectivity, atom_connectivity, md_array_size(atom_connectivity), alloc);
            }
        }

        md_spatial_hash_free(&sh);
    }

    bond_data->count = md_array_size(bond_data->bond);
    return true;
}

bool md_util_compute_chain_data(md_chain_data_t* chain_data, const md_residue_idx_t* res_idx, int64_t atom_count, const md_bond_t* bonds, int64_t bond_count, md_allocator_i* alloc) {
    if (!chain_data) {
        MD_LOG_ERROR("chain data is missing");
        return false;
    }
    
    if (!res_idx) {
        MD_LOG_ERROR("residue index is null");
        return false;
    }

    if (!bonds) {
        MD_LOG_ERROR("bonds is null");
        return false;
    }

    md_allocator_i* temp_alloc = default_temp_allocator;
    
    md_array(md_range_t) res_atom_range = 0;
    
    int res_count = 0;
    for (int i = 0; i < (int)atom_count; ++i) {
        res_count = MAX(res_count, res_idx[i]);
        if (res_count >= md_array_size(res_atom_range)) {
            md_range_t atom_range = {i, i};
            md_array_push(res_atom_range, atom_range, temp_alloc);
        }
        res_atom_range[res_idx[i]].end += 1;
    }
    
    md_array(bool) res_bond_to_next = 0;
    md_array_resize(res_bond_to_next, res_count, temp_alloc);
    MEMSET(res_bond_to_next, 0, md_array_bytes(res_bond_to_next));

    // Iterate over bonds and determine bonds between residues (res_bonds)
    for (int64_t i = 0; i < bond_count; ++i) {
        const md_bond_t bond = bonds[i];
        const md_residue_idx_t res_a = res_idx[bond.idx[0]];
        const md_residue_idx_t res_b = res_idx[bond.idx[1]];
        
        if (res_b == res_a + 1) {
            res_bond_to_next[res_a] = true;
        }
    }

    md_array_shrink(chain_data->id, 0);
    md_array_shrink(chain_data->atom_range, 0);
    md_array_shrink(chain_data->residue_range, 0);
    chain_data->count = 0;

    char next_char = 0;
    int beg_idx = 0;
    for (int i = 0; i < res_count; ++i) {
        if (res_bond_to_next[i] == false || i == res_count - 1) {
            int end_idx = i + 1;
            if (end_idx - beg_idx > 1) {
                const md_label_t id = {.buf = {'A' + next_char}, .len = 1};
                const md_range_t atom_range = {res_atom_range[beg_idx].beg, res_atom_range[end_idx - 1].end};
                const md_range_t res_range  = {beg_idx, end_idx};
                
                md_array_push(chain_data->id, id, alloc);
                md_array_push(chain_data->atom_range, atom_range, alloc);
                md_array_push(chain_data->residue_range, res_range, alloc);
                
                next_char = (next_char + 1) % 26;
                chain_data->count += 1;
            }
            beg_idx = i + 1;
        }
    }

    return chain_data;
}

bool md_util_compute_atom_valence(md_valence_t atom_valence[], int64_t atom_count, const md_bond_t bonds[], int64_t bond_count) {
    if (!atom_valence) {
        MD_LOG_ERROR("Missing input: atom_valence");
        return false;
    }

    if (atom_count < 0) {
        MD_LOG_ERROR("Invalid input: atom_count");
        return false;
    }

    if (!bonds) {
        MD_LOG_ERROR("Missing input: bonds");
        return false;
    }

    if (bond_count < 0) {
        MD_LOG_ERROR("Invalid input: bond_count");
        return false;
    }

    MEMSET(atom_valence, 0, sizeof(md_valence_t) * atom_count);

    for (int64_t i = 0; i < bond_count; ++i) {
        const md_bond_t bond = bonds[i];
        atom_valence[bond.idx[0]] += 1;
        atom_valence[bond.idx[1]] += 1;
    }

    return true;
}


typedef struct fifo_t {
    md_allocator_i* alloc;
    int* data;
    unsigned int head;
    unsigned int tail;
    unsigned int cap;
} fifo_t;

static bool fifo_empty(fifo_t* fifo) { return fifo->head == fifo->tail; }
static bool fifo_full(fifo_t* fifo)  { return ((fifo->head + 1) & (fifo->cap - 1)) == fifo->tail; }

static void fifo_grow(fifo_t* fifo, int64_t new_capacity) {
    uint32_t new_cap = next_power_of_two32((uint32_t)new_capacity);
    fifo->data = md_realloc(fifo->alloc, fifo->data, fifo->cap, new_cap);
    fifo->cap = new_cap;
}

static void fifo_init(fifo_t* fifo, int64_t capacity, md_allocator_i* alloc) {
    ASSERT(fifo);
    ASSERT(0 <= capacity && capacity < UINT32_MAX);
    fifo->alloc = alloc;
    fifo->head = 0;
    fifo->tail = 0;
    fifo->cap = next_power_of_two32((uint32_t)capacity);
    fifo->data = md_alloc(alloc, sizeof(int) * fifo->cap);
    MEMSET(fifo->data, 0, sizeof(int) * fifo->cap);
}

static inline fifo_t fifo_create(int64_t capacity, md_allocator_i* alloc) {
    fifo_t fifo;
    fifo_init(&fifo, capacity, alloc);
    return fifo;
}

static void fifo_free(fifo_t* fifo) {
    ASSERT(fifo);
    md_free(fifo->alloc, fifo->data, sizeof(int) * fifo->cap);
    MEMSET(fifo, 0, sizeof(fifo_t));
}

static void fifo_clear(fifo_t* fifo) {
    fifo->head = fifo->tail = 0;
}

static void fifo_push(fifo_t* fifo, int value) {
    if (fifo_full(fifo)) {
        fifo_grow(fifo, fifo->cap * 2);
    }
    fifo->data[fifo->head] = value;
    fifo->head = (fifo->head + 1) & (fifo->cap - 1);
}

static int fifo_pop(fifo_t* fifo) {
    ASSERT(!fifo_empty(fifo));
    int val = fifo->data[fifo->tail];
    fifo->tail = (fifo->tail + 1) & (fifo->cap - 1);
    return val;
}

static inline bool compare_ring(const int* ring_a, int ring_a_size, const int* ring_b, int ring_b_size) {
    if (ring_a_size != ring_b_size) {
        return false;
    }
    for (int i = 0; i < ring_a_size; ++i) {
        if (ring_a[i] != ring_b[i]) {
            return false;
        }
    }
    return true;
}

static bool has_ring(const md_index_data_t* ring_data, const int* ring, int ring_size) {
    const int64_t num_rings = md_index_data_count(ring_data);
    for (int64_t i = 0; i < num_rings; ++i) {
        const int* ring_i     = md_index_range_beg(ring_data, i);
        const int ring_i_size = (int)md_index_range_size(ring_data, i);
        if (compare_ring(ring_i, ring_i_size, ring, ring_size)) return true;
    }
    return false;
}

// Simplistic inplace bubble sort for small arrays
static void sort_arr(int* arr, int n) {
    bool swapped = true;
    while (swapped) {
        swapped = false;
        for (int i = 0; i < n - 1; ++i) {
            if (arr[i] > arr[i + 1]) {
                int tmp = arr[i];
                arr[i] = arr[i + 1];
                arr[i + 1] = tmp;
                swapped = true;
            }
        }
    }
}

#define MIN_RING_SIZE 3
#define MAX_RING_SIZE 6

bool md_util_compute_rings(md_index_data_t* ring_data, int64_t atom_count, const md_bond_t bonds[], int64_t bond_count, md_allocator_i* alloc) {
    ASSERT(alloc);
    
    if (!ring_data) {
        MD_LOG_ERROR("ring data is missing");
        return false;
    }

    if (!bonds) {
        MD_LOG_ERROR("bonds are missing");
        return false;
    }

    if (bond_count < 0) {
        MD_LOG_ERROR("invalid bond count");
        return false;
    }
    
    md_index_data_data_clear(ring_data);
    
    md_vm_arena_t arena = {0};
    md_vm_arena_init(&arena, GIGABYTES(1));
    md_allocator_i arena_alloc = md_vm_arena_create_interface(&arena);

    md_index_data_t edge_data = compute_connectivity_data(atom_count, bonds, bond_count, &arena_alloc);
    
    int* depth = md_vm_arena_push(&arena, atom_count * sizeof(int));
    MEMSET(depth, 0, atom_count * sizeof(int));
    
    int* pred = md_vm_arena_push(&arena, atom_count * sizeof(int));
    for (int64_t i = 0; i < atom_count; ++i) {
        pred[i] = -1;
    }

    // The capacity is arbitrary here, but will be resized if needed.
    fifo_t queue = fifo_create(64, &arena_alloc);

    for (int i = 0; i < atom_count; ++i) {
        // Skip any atom which has alread been touched
        if (depth[i]) continue;

        fifo_clear(&queue);
        fifo_push(&queue, i);
        depth[i] = 1;
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);

            const int* eb = md_index_range_beg(&edge_data, idx);
            const int* ee = md_index_range_end(&edge_data, idx);
            for (const int* it = eb; it != ee; ++it) {
                int next = *it;
                if (next == pred[idx]) continue;  // avoid adding parent to search queue

                if (depth[next]) {
                    // We found a junction point where the graph connects
                    // Now we need to traverse both branches of this graph back up
                    // In order to find the other junction where the graph branches
                    // We follow both branches up towards the root.

                    int l = idx;
                    int r = next;

                    // Only process one of the two cases
                    if (r < l) continue;

                    int n = 0;
                    int ring[MAX_RING_SIZE + 3];

                    int d = MAX(depth[l], depth[r]);
                    while (d-- && n < MAX_RING_SIZE) {
                        if (depth[l] >= d) {
                            ring[n++] = l;
                            l = pred[l];
                        }

                        if (depth[r] >= d) {
                            ring[n++] = r;
                            r = pred[r];
                        }

                        if (l == r) {
                            ring[n++] = l;
                            break;
                        }
                    }

                    if (MIN_RING_SIZE <= n && n <= MAX_RING_SIZE) {
                        sort_arr(ring, n);
                        if (!has_ring(ring_data, ring, n)) {
                            md_index_data_push(ring_data, ring, n, alloc);
                        }
                    }

                } else {
                    depth[next] = depth[idx] + 1;
                    pred[next]  = idx;
                    fifo_push(&queue, next);
                }
            }
        }
    }

    md_vm_arena_free(&arena);

    return true;
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE

bool md_util_compute_structures(md_index_data_t* structure_data, int64_t atom_count, const md_bond_t bonds[], int64_t bond_count, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    if (!structure_data) {
        MD_LOG_ERROR("structure data is missing");
        return false;
    }

    if (!bonds) {
        MD_LOG_ERROR("bonds are missing");
        return false;
    }

    if (bond_count < 0) {
        MD_LOG_ERROR("invalid bond count");
        return false;
    }

    md_index_data_data_clear(structure_data);

    md_vm_arena_t arena = { 0 };
    md_vm_arena_init(&arena, GIGABYTES(1));
    md_allocator_i arena_alloc = md_vm_arena_create_interface(&arena);

    md_index_data_t edge_data = compute_connectivity_data(atom_count, bonds, bond_count, &arena_alloc);

    int* depth = md_vm_arena_push(&arena, atom_count * sizeof(int));
    MEMSET(depth, 0, atom_count * sizeof(int));

    // The capacity is arbitrary here, but will be resized if needed.
    fifo_t queue = fifo_create(64, &arena_alloc);

    md_array(int) indices = 0;

    for (int i = 0; i < atom_count; ++i) {
        // Skip any atom which has already been touched
        if (depth[i]) continue;

        md_array_shrink(indices, 0);
        
        fifo_clear(&queue);
        fifo_push(&queue, i);
        depth[i] = 1;
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);
            md_array_push(indices, idx, &arena_alloc);
            
            const int* eb = md_index_range_beg(&edge_data, idx);
            const int* ee = md_index_range_end(&edge_data, idx);
            for (const int* it = eb; it != ee; ++it) {
                int next = *it;
                if (depth[next] == 0) {
                    depth[next] = depth[idx] + 1;
                    fifo_push(&queue, next);
                }
            }
        }

        sort_arr(indices, (int)md_array_size(indices));
        
        // Here we should have exhausted every atom that is connected to index i.
        md_index_data_push(structure_data, indices, md_array_size(indices), alloc);
    }
    
    md_vm_arena_free(&arena);

    return true;
}

void md_util_compute_aabb_xyz(vec3_t* aabb_min, vec3_t* aabb_max, const float* in_x, const float* in_y, const float* in_z, int64_t count) {
    md_simd_f32_t vx_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vy_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vz_min = md_simd_set1_f32(+FLT_MAX);

    md_simd_f32_t vx_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vy_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vz_max = md_simd_set1_f32(-FLT_MAX);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, md_simd_f32_width);
    for (; i < simd_count; i += md_simd_f32_width) {
        md_simd_f32_t x = md_simd_load_f32(in_x + i);
        md_simd_f32_t y = md_simd_load_f32(in_y + i);
        md_simd_f32_t z = md_simd_load_f32(in_z + i);

        vx_min = md_simd_min(vx_min, x);
        vy_min = md_simd_min(vy_min, y);
        vz_min = md_simd_min(vz_min, z);

        vx_max = md_simd_max(vx_max, x);
        vy_max = md_simd_max(vy_max, y);
        vz_max = md_simd_max(vz_max, z);
    }

    float x_min = md_simd_hmin(vx_min);
    float y_min = md_simd_hmin(vy_min);
    float z_min = md_simd_hmin(vz_min);

    float x_max = md_simd_hmax(vx_max);
    float y_max = md_simd_hmax(vy_max);
    float z_max = md_simd_hmax(vz_max);

    for (; i < count; ++i) {
        x_min = MIN(x_min, in_x[i]);
        y_min = MIN(y_min, in_y[i]);
        z_min = MIN(z_min, in_z[i]);

        x_max = MAX(x_max, in_x[i]);
        y_max = MAX(y_max, in_y[i]);
        z_max = MAX(z_max, in_z[i]);
    }

    *aabb_min = (vec3_t){x_min, y_min, z_min};
    *aabb_max = (vec3_t){x_max, y_max, z_max};
}

void md_util_compute_aabb_xyzr(vec3_t* aabb_min, vec3_t* aabb_max, const float* in_x, const float* in_y, const float* in_z, const float* in_r, int64_t count) {
    md_simd_f32_t vx_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vy_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vz_min = md_simd_set1_f32(+FLT_MAX);

    md_simd_f32_t vx_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vy_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vz_max = md_simd_set1_f32(-FLT_MAX);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, md_simd_f32_width);
    for (; i < simd_count; i += md_simd_f32_width) {
        md_simd_f32_t x = md_simd_load_f32(in_x + i);
        md_simd_f32_t y = md_simd_load_f32(in_y + i);
        md_simd_f32_t z = md_simd_load_f32(in_z + i);
        md_simd_f32_t r = md_simd_load_f32(in_r + i);

        vx_min = md_simd_min(vx_min, md_simd_sub(x, r));
        vy_min = md_simd_min(vy_min, md_simd_sub(y, r));
        vz_min = md_simd_min(vz_min, md_simd_sub(z, r));

        vx_max = md_simd_max(vx_max, md_simd_add(x, r));
        vy_max = md_simd_max(vy_max, md_simd_add(y, r));
        vz_max = md_simd_max(vz_max, md_simd_add(z, r));
    }

    float x_min = md_simd_hmin(vx_min);
    float y_min = md_simd_hmin(vy_min);
    float z_min = md_simd_hmin(vz_min);

    float x_max = md_simd_hmax(vx_max);
    float y_max = md_simd_hmax(vy_max);
    float z_max = md_simd_hmax(vz_max);

    for (; i < count; ++i) {
        float r = in_r[i];
        x_min = MIN(x_min, in_x[i] - r);
        y_min = MIN(y_min, in_y[i] - r);
        z_min = MIN(z_min, in_z[i] - r);

        x_max = MAX(x_max, in_x[i] + r);
        y_max = MAX(y_max, in_y[i] + r);
        z_max = MAX(z_max, in_z[i] + r);
    }

    *aabb_min = (vec3_t){x_min, y_min, z_min};
    *aabb_max = (vec3_t){x_max, y_max, z_max};
}

void md_util_compute_aabb_ortho_xyz(vec3_t* out_aabb_min, vec3_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, int64_t count, uint64_t stride, vec3_t pbc_ext) {
    md_simd_f32_t vx_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vy_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vz_min = md_simd_set1_f32(+FLT_MAX);

    md_simd_f32_t vx_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vy_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vz_max = md_simd_set1_f32(-FLT_MAX);

    vec4_t ext = vec4_from_vec3(pbc_ext, 0);
    vec4_t ref = vec4_mul_f(ext, 0.5f);

    int64_t i = 0;
    const int64_t simd_count = (count / md_simd_f32_width) * md_simd_f32_width;
    for (; i < simd_count; i += md_simd_f32_width) {
        md_simd_f32_t x = md_simd_load_f32((const float*)((const char*)in_x + i * stride));
        md_simd_f32_t y = md_simd_load_f32((const float*)((const char*)in_y + i * stride));
        md_simd_f32_t z = md_simd_load_f32((const float*)((const char*)in_z + i * stride));

        x = md_simd_deperiodize(x, md_simd_set1_f32(ref.x), md_simd_set1_f32(ext.x));
        y = md_simd_deperiodize(y, md_simd_set1_f32(ref.y), md_simd_set1_f32(ext.y));
        z = md_simd_deperiodize(z, md_simd_set1_f32(ref.z), md_simd_set1_f32(ext.z));

        vx_min = md_simd_min(vx_min, x);
        vy_min = md_simd_min(vy_min, y);
        vz_min = md_simd_min(vz_min, z);

        vx_max = md_simd_max(vx_max, x);
        vy_max = md_simd_max(vy_max, y);
        vz_max = md_simd_max(vz_max, z);
    }

    vec4_t aabb_min = {
        md_simd_hmin(vx_min),
        md_simd_hmin(vy_min),
        md_simd_hmin(vz_min),
        0
    };

    vec4_t aabb_max = {
        md_simd_hmax(vx_max),
        md_simd_hmax(vy_max),
        md_simd_hmax(vz_max),
        0
    };

    // Handle remainder
    for (; i < count; ++i) {
        vec4_t c = {
            *(const float*)((const char*)in_x + i * stride),
            *(const float*)((const char*)in_y + i * stride),
            *(const float*)((const char*)in_z + i * stride),
            0
        };
        c = vec4_deperiodize(c, ref, ext);
        aabb_min = vec4_min(aabb_min, c);
        aabb_max = vec4_max(aabb_max, c);
    }

    *out_aabb_min = vec3_from_vec4(aabb_min);
    *out_aabb_max = vec3_from_vec4(aabb_max);
}

vec3_t md_util_compute_com(const vec3_t* xyz, const float* w, int64_t count) {
    ASSERT(xyz);
    ASSERT(count >= 0);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    // Use vec4 here so we can utilize SSE vectorization if applicable
    // @TODO: Vectorize with full register width
    vec4_t sum_xyzw = {0,0,0,0};
    if (w) {
        for (int64_t i = 0; i < count; ++i) {
            vec4_t xyz1 = vec4_from_vec3(xyz[i], 1.0f);
            sum_xyzw = vec4_add(sum_xyzw, vec4_mul_f(xyz1, w[i]));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            vec4_t xyz0 = vec4_from_vec3(xyz[i], 0);
            sum_xyzw = vec4_add(sum_xyzw, xyz0);
        }
        sum_xyzw.w = (float)count;
    }

    return vec3_div_f(vec3_from_vec4(sum_xyzw), sum_xyzw.w);
}

vec3_t md_util_compute_com_soa(const float* x, const float* y, const float* z, const float* w, int64_t count) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(count >= 0);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    // Use vec4 here so we can utilize SSE vectorization if applicable
    // @TODO: Vectorize with full register width
    vec4_t sum_xyzw = {0,0,0,0};
    if (w) {
        for (int64_t i = 0; i < count; ++i) {
            const vec4_t xyz1 = vec4_set(x[i], y[i], z[i], 1.0f);
            sum_xyzw = vec4_add(sum_xyzw, vec4_mul_f(xyz1, w[i]));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            const vec4_t xyz0 = vec4_set(x[i], y[i], z[i], 0);
            sum_xyzw = vec4_add(sum_xyzw, xyz0);
        }
        sum_xyzw.w = (float)count;
    }

    return vec3_div_f(vec3_from_vec4(sum_xyzw), sum_xyzw.w);
}

vec3_t md_util_compute_com_indexed_soa(const float *x, const float* y, const float* z, const float* w, const int32_t* indices, int64_t index_count) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(indices);
    ASSERT(index_count >= 0);

    if (index_count == 0) {
        return (vec3_t) {0,0,0};
    }

    // Use vec4 here so we can utilize SSE vectorization if applicable
    // @TODO: Vectorize with full register width
    vec4_t sum_xyzw = {0,0,0,0};
    if (w) {
        for (int64_t i = 0; i < index_count; ++i) {
            const int idx = indices[i];
            const vec4_t xyz1 = vec4_set(x[idx], y[idx], z[idx], 1.0f);
            sum_xyzw = vec4_add(sum_xyzw, vec4_mul_f(xyz1, w[idx]));
        }
    } else {
        for (int64_t i = 0; i < index_count; ++i) {
            const int idx = indices[i];
            const vec4_t xyz0 = vec4_set(x[idx], y[idx], z[idx], 0);
            sum_xyzw = vec4_add(sum_xyzw, xyz0);
        }
        sum_xyzw.w = (float)index_count;
    }

    return vec3_div_f(vec3_from_vec4(sum_xyzw), sum_xyzw.w);
}

static inline double compute_com_periodic(const float* in_x, const float* in_w, int64_t count, float pbc_ext, size_t stride_x) {
    double com;

    double acc_c = 0;
    double acc_s = 0;
    double acc_w = 0;

    const double scl = TWO_PI / pbc_ext;
    stride_x = MAX(1, stride_x);

    for (int64_t i = 0; i < count; ++i) {
        double theta = in_x[stride_x * i] * scl;
        double w = in_w ? in_w[i] : 1.0;
        acc_c += w * cos(theta);
        acc_s += w * sin(theta);
        acc_w += w;
    }

    acc_w = (acc_w == 0) ? count : acc_w;
    const double theta_prim = atan2(-acc_s / acc_w, -acc_c / acc_w) + PI;
    com = (theta_prim / TWO_PI) * pbc_ext;

    return com;
}

static inline double compute_com(const float* in_x, const float* in_w, int64_t count, size_t stride_x) {
    double com = 0;
    double acc_x = 0;
    double acc_w = 0;

    stride_x = MAX(1, stride_x);

    for (int64_t i = 0; i < count; ++i) {
        double x = in_x[stride_x * i];
        double w = in_w ? in_w[i] : 1.0;
        acc_x += x * w;
        acc_w += w;
    }

    com = acc_x / acc_w;
    return com;
}

// Love the elegance of using trigonometric functions, unsure of the performance...
// @TODO: sin, cos and atan2 can and should of course be vectorized.
vec3_t md_util_compute_com_ortho(const vec3_t* in_xyz, const float* in_w, int64_t count, vec3_t pbc_ext) {
    ASSERT(in_xyz);
    ASSERT(count >= 0);
    ASSERT(pbc_ext.x >= 0);
    ASSERT(pbc_ext.y >= 0);
    ASSERT(pbc_ext.z >= 0);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    // We need to pick the version based on each component of pdc_ext, since one or more may be zero
    double x = pbc_ext.x > 0 ? compute_com_periodic(in_xyz->elem + 0, in_w, count, pbc_ext.x, 3) : compute_com(in_xyz->elem + 0, in_w, count, 3);
    double y = pbc_ext.y > 0 ? compute_com_periodic(in_xyz->elem + 1, in_w, count, pbc_ext.y, 3) : compute_com(in_xyz->elem + 1, in_w, count, 3);
    double z = pbc_ext.z > 0 ? compute_com_periodic(in_xyz->elem + 2, in_w, count, pbc_ext.z, 3) : compute_com(in_xyz->elem + 2, in_w, count, 3);

    return (vec3_t) {(float)x, (float)y, (float)z};
}

vec3_t md_util_compute_com_ortho_soa(const float* in_x, const float* in_y, const float* in_z, const float* in_w, int64_t count, vec3_t pbc_ext) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(count >= 0);
    ASSERT(pbc_ext.x >= 0);
    ASSERT(pbc_ext.y >= 0);
    ASSERT(pbc_ext.z >= 0);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    // We need to pick the version based on each component of pdc_ext, since one or more may be zero
    double x = pbc_ext.x > 0 ? compute_com_periodic(in_x, in_w, count, pbc_ext.x, 1) : compute_com(in_x, in_w, count, 1);
    double y = pbc_ext.y > 0 ? compute_com_periodic(in_y, in_w, count, pbc_ext.y, 1) : compute_com(in_y, in_w, count, 1);
    double z = pbc_ext.z > 0 ? compute_com_periodic(in_z, in_w, count, pbc_ext.z, 1) : compute_com(in_z, in_w, count, 1);

    return (vec3_t) {(float)x, (float)y, (float)z};
}

static inline vec4_t minimum_image(vec4_t x, vec4_t ext, vec4_t inv_ext) {
    const vec4_t s = vec4_mul(x, inv_ext);
    const vec4_t dx = vec4_mul(ext, vec4_sub(s, vec4_round(s)));
    return vec4_blend(dx, x, vec4_cmp_eq(ext, vec4_zero()));
}

vec3_t md_util_compute_com_ortho_indexed_soa(const float *in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* indices, int64_t count, vec3_t box) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(indices);
    ASSERT(count >= 0);

    if (count == 0) {
        return (vec3_t){0, 0, 0};
    }
    
    const vec4_t ext     = vec4_from_vec3(box, 0);
    const vec4_t inv_ext = vec4_set(box.x ? 1.0f / box.x : 0, box.y ? 1.0f / box.y : 0, box.z ? 1.0f / box.z : 0, 0);

    vec4_t xyzw_acc = {in_x[0], in_y[0], in_z[0], in_w ? in_w[0] : 1};

    if (in_w) {
        for (int64_t i = 1; i < count; ++i) {
            const int idx = indices[i];
            const vec4_t com     = vec4_div_f(xyzw_acc, xyzw_acc.w);
            const vec4_t xyzw    = vec4_deperiodize(vec4_set(in_x[idx], in_y[idx], in_z[idx], 1.0f), com, ext);
            xyzw_acc = vec4_add(xyzw_acc, vec4_mul_f(xyzw, in_w[idx]));
        }
    } else {
        for (int64_t i = 1; i < count; ++i) {
            const int idx = indices[i];
            const vec4_t com     = vec4_div_f(xyzw_acc, xyzw_acc.w);
            const vec4_t xyzw    = vec4_deperiodize(vec4_set(in_x[idx], in_y[idx], in_z[idx], 1.0f), com, ext);
            xyzw_acc = vec4_add(xyzw_acc, xyzw);
        }
    }


    return vec3_from_vec4(vec4_div_f(xyzw_acc, xyzw_acc.w));
}

/*
vec3_t md_util_compute_com_ortho(const vec3_t* in_xyz, const float* in_w, int64_t count, vec3_t pbc_ext) {
    ASSERT(in_xyz);

    float w = in_w ? in_w[0] : 1.0f;
    float sum_x = in_xyz[0].x * w;
    float sum_y = in_xyz[0].y * w;
    float sum_z = in_xyz[0].z * w;
    float sum_w = w;

    for (int64_t i = 1; i < count; ++i) {
        float x = pbc_ext.x > 0 ? deperiodizef(in_xyz[i].x, sum_x / sum_w, pbc_ext.x) : in_xyz[i].x;
        float y = pbc_ext.y > 0 ? deperiodizef(in_xyz[i].y, sum_y / sum_w, pbc_ext.y) : in_xyz[i].y;
        float z = pbc_ext.z > 0 ? deperiodizef(in_xyz[i].z, sum_z / sum_w, pbc_ext.z) : in_xyz[i].z;

        w = in_w ? in_w[i] : 1.0f;
        sum_x += x * w;
        sum_y += y * w;
        sum_z += z * w;
        sum_w += w;
    }

    // Ensure the resulting com is within the period of the box
    float com_x = sum_x / sum_w;
    float com_y = sum_y / sum_w;
    float com_z = sum_z / sum_w;

    return (vec3_t){com_x, com_y, com_z};
}

vec3_t md_util_compute_com_ortho_soa(const float* in_x, const float* in_y, const float* in_z, const float* in_w, int64_t count, vec3_t pbc_ext) {
    float w = in_w ? in_w[0] : 1.0f;
    float sum_x = in_x[0] * w;
    float sum_y = in_y[0] * w;
    float sum_z = in_z[0] * w;
    float sum_w = w;

    for (int64_t i = 1; i < count; ++i) {
        float x = pbc_ext.x > 0 ? deperiodizef(in_x[i], sum_x / sum_w, pbc_ext.x) : in_x[i];
        float y = pbc_ext.y > 0 ? deperiodizef(in_y[i], sum_y / sum_w, pbc_ext.y) : in_y[i];
        float z = pbc_ext.z > 0 ? deperiodizef(in_z[i], sum_z / sum_w, pbc_ext.z) : in_z[i];

        w = in_w ? in_w[i] : 1.0f;
        sum_x += x * w;
        sum_y += y * w;
        sum_z += z * w;
        sum_w += w;
    }

    float com_x = sum_x / sum_w;
    float com_y = sum_y / sum_w;
    float com_z = sum_z / sum_w;

    return (vec3_t){com_x, com_y, com_z};
}
*/

// From here: https://www.arianarab.com/post/crystal-structure-software
mat3_t md_util_compute_unit_cell_basis(double a, double b, double c, double alpha, double beta, double gamma) {
    alpha = DEG_TO_RAD(alpha);
    beta  = DEG_TO_RAD(beta);
    gamma = DEG_TO_RAD(gamma);

    const double cb = cos(beta);
    const double x = (cos(alpha) - cb * cos(gamma)) / sin(gamma);
    mat3_t M = {
        .col = {
            {(float)a, 0, 0},
            {(float)(b * cos(gamma)), (float)(b * sin(gamma)), 0},
            {(float)(c * cb), (float)(c * x), (float)(c * sqrt(1 - cb * cb - x * x))},
        },
    };
    return M;
}

// This is a bit silly, but I want to make sure we are conservative and not always compute the length of the vectors
// In the odd case that we don't exactly get back what we had due to precision errors sqrtf(^2).
vec3_t md_util_compute_unit_cell_extent(mat3_t M) {
    vec3_t ext = {0};

    if (M.elem[0][1] == 0 && M.elem[0][2] == 0)
        ext.x = M.elem[0][0];
    else
        ext.x = vec3_length(M.col[0]);

    if (M.elem[1][0] == 0 && M.elem[1][2] == 0)
        ext.y = M.elem[1][1];
    else
        ext.y = vec3_length(M.col[1]);

    if (M.elem[2][0] == 0 && M.elem[2][1] == 0)
        ext.z = M.elem[2][2];
    else
        ext.z = vec3_length(M.col[2]);

    return ext;
}

bool md_util_pbc_ortho(float* x, float* y, float* z, int64_t count, vec3_t box_ext) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (count < 0) {
        MD_LOG_ERROR("Invalid count");
        return false;
    }

    if (box_ext.x < 0 || box_ext.y < 0 || box_ext.z < 0) {
        MD_LOG_ERROR("Invalid PBC extent: One or more components were negative.");
        return false;
    }

    if (vec3_equal(box_ext, vec3_zero())) {
        MD_LOG_ERROR("PBC extent is zero, unable to apply it!");
        return false;
    }

    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    const vec4_t ref = vec4_mul_f(ext, 0.5f);
    for (int64_t i = 0; i < count; ++i) {
        vec4_t pos = {x[i], y[i], z[i], 0};
        pos = vec4_deperiodize(pos, ref, ext);
        x[i] = pos.x;
        y[i] = pos.y;
        z[i] = pos.z;
    }

    return true;
}

bool md_util_unwrap_ortho(float* x, float* y, float* z, const md_index_data_t* structures, vec3_t box_ext) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (!structures) {
        MD_LOG_ERROR("Missing required input: structures");
        return false;
    }

    if (box_ext.x < 0 || box_ext.y < 0 || box_ext.z < 0) {
        MD_LOG_ERROR("Invalid box extent: One or more components were negative.");
        return false;
    }

    if (vec3_equal(box_ext, vec3_zero())) {
        MD_LOG_ERROR("Box extent is zero.");
        return false;
    }

    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    
    const int64_t num_structures = md_index_data_count(structures);
    for (int64_t s_idx = 0; s_idx < num_structures; ++s_idx) {
        const int* indices = md_index_range_beg(structures, s_idx);
        const int count = (int)md_index_range_size(structures, s_idx);
        
        // Indices are sorted with respect to covalent bond 'depth' from the first atom
        // This means we can use the inherent order of indices to iterate and unwrap.
        int idx = indices[0];
        vec4_t ref_pos = {x[idx], y[idx], z[idx], 0};
        for (int i = 1; i < count; ++i) {
            idx = indices[i];
            const vec4_t pos = vec4_deperiodize((vec4_t){x[idx], y[idx], z[idx], 0}, ref_pos, ext);
            x[idx] = pos.x;
            y[idx] = pos.y;
            z[idx] = pos.z;
            ref_pos = pos;
        }
    }
    return true;
}

bool md_util_deperiodize_system_ortho(float* x, float* y, float* z, const float* w, int64_t count, const md_index_data_t* structures, vec3_t box) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    
    md_util_pbc_ortho(x,y,z, count, box);
    if (structures) {
        md_util_unwrap_ortho(x,y,z, structures, box);
    
        // Place any defined covalent structure such that its center of mass is within the correct periodic image
        const vec4_t ext = vec4_from_vec3(box, 0);
        const vec4_t ref = vec4_mul_f(ext, 0.5f);

        const int64_t num_structures = md_index_data_count(structures);
        for (int64_t s_idx = 0; s_idx < num_structures; ++s_idx) {
            const int32_t* indices = md_index_range_beg(structures, s_idx);
            const int64_t num_indices = md_index_range_size(structures, s_idx);

            const vec4_t com = vec4_from_vec3(md_util_compute_com_indexed_soa(x, y, z, w, indices, num_indices), 0);
            const vec4_t pbc_com = vec4_deperiodize(com, ref, ext);
            const vec4_t delta = vec4_sub(pbc_com, com);
            if (vec4_dot(delta, delta) > 0.001f) {
                for (int64_t i = 0; i < num_indices; ++i) {
                    const int32_t j = indices[i];
                    x[j] += delta.x;
                    y[j] += delta.y;
                    z[j] += delta.z;
                }
            }
        }
    }

    return true;
}

#if 0
bool md_util_apply_pbc_preserve_covalent(float* x, float* y, float* z, const md_bond_data_t* covalent, vec3_t pbc_ext) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (!covalent) {
        MD_LOG_ERROR("Missing required input: covalent_data");
        return false;
    }

    if (pbc_ext.x < 0 || pbc_ext.y < 0 || pbc_ext.z < 0) {
        MD_LOG_ERROR("Invalid PBC extent: One or more components were negative.");
        return false;
    }

    if (vec3_equal(pbc_ext, vec3_zero())) {
        MD_LOG_ERROR("PBC extent is zero, unable to apply it!");
        return false;
    }

    const vec4_t ext = vec4_from_vec3(pbc_ext, 0);

    fifo_t queue = fifo_create(64, default_temp_allocator);
    md_array(int)  pred = 0;

    md_array(vec3_t) xyz = 0;
    const int64_t num_structures = md_index_data_count(&covalent->structures);
    for (int64_t s_idx = 0; s_idx < num_structures; ++s_idx) {
        md_array_shrink(xyz, 0);
        md_array_shrink(pred, 0);
        
        const md_atom_idx_t* beg_idx = md_index_range_beg(&covalent->structures, s_idx);
        const md_atom_idx_t* end_idx = md_index_range_end(&covalent->structures, s_idx);
        
        for (const md_atom_idx_t* it = beg_idx; it != end_idx; ++it) {
            int i = *it;
            const vec3_t pos = {x[i], y[i], z[i]};
            md_array_push(xyz, pos, default_temp_allocator);
            md_array_push(pred, -1, default_temp_allocator);
        }
        
        //const vec4_t com = vec4_from_vec3(md_util_compute_com_ortho(xyz, NULL, md_array_size(xyz), pbc_ext), 0);
        const vec4_t com = vec4_mul_f(ext, 0.5f);
        
        /*
        for (int64_t i = 0; i < md_array_size(xyz); ++i) {
            const vec4_t pos = vec4_from_vec3(xyz[i], 0);
            const vec4_t dp  = vec4_deperiodize(pos, ref, vec4_from_vec3(pbc_ext, 0));
            const int idx = beg_idx[i];
            x[idx] = dp.x;
            y[idx] = dp.y;
            z[idx] = dp.z;
        }
        */

        // NOTE: This is some black box magic, dark voodoo shit.
        // The problem is: Even if we have a valid center of mass of our structure, and we would like to deperiodize the atoms within the structure
        // such that every atom reside within the same period as the center of mass, it is not possible if the structure has an extent which is more than
        // half a period. In such case, the _deperiodize function will fail and place the atom in the next period, since that was closer to the reference (com).
        // Since we know that every atom within this structure is connected by covalent bonds, we can pick the atom which is the closest to the center of mass.
        // Then perform a breadth-first traversal which propagates outwards from this reference atom. When iterating, we deperiodize the atoms with respect to its
        // Predecessor. This should enforce that the entire structure resides within the same period (Even if it was split into periods before).

        // Find index of atom which is closest to com
        int ref_idx = -1;
        float ref_dist = FLT_MAX;
        for (int i = 0; i < md_array_size(xyz); ++i) {
            const float d2 = vec4_distance_squared(com, vec4_from_vec3(xyz[i], 0));
            if (d2 < ref_dist) {
                ref_dist = d2;
                ref_idx = i;
            }
        }

        fifo_clear(&queue);
        fifo_push(&queue, ref_idx);
        pred[ref_idx] = ref_idx;
        
        while (!fifo_empty(&queue)) {
            const int idx = fifo_pop(&queue);

            const vec4_t ref = vec4_from_vec3(xyz[idx], 0);

            const int* n_beg = md_index_range_beg(&covalent->connectivity, idx);
            const int* n_end = md_index_range_end(&covalent->connectivity, idx);
            for (const int* it = n_beg; it < n_end; ++it) {
                const int neighbor_idx = *it;
                if (pred[neighbor_idx] != -1) continue;
                pred[neighbor_idx] = idx;

                const vec4_t p = vec4_from_vec3(xyz[neighbor_idx], 0);
                const vec4_t dp = vec4_deperiodize(p, ref, ext);
                const int i = beg_idx[idx];
                x[i] = dp.x;
                y[i] = dp.y;
                z[i] = dp.z;
                
                fifo_push(&queue, neighbor_idx);
            }
        }        
    }
    
    return true;
}
#endif

/*
bool md_util_pbc_ortho(md_molecule_t* mol, vec3_t pbc_ext) {
    if (!mol) return false;

    vec3_t* residue_com = 0;
    if (mol->chain.residue_range && mol->chain.count) {
        md_array_resize(residue_com, mol->residue.count, default_allocator);
    }
    
    for (int64_t i = 0; i < mol->residue.count; ++i) {
        md_range_t atom_range = mol->residue.atom_range[i];

        vec3_t com = md_util_compute_com_ortho_soa(
            mol->atom.x + atom_range.beg,
            mol->atom.y + atom_range.beg,
            mol->atom.z + atom_range.beg,
            mol->atom.mass ? mol->atom.mass + atom_range.beg : NULL,
            atom_range.end - atom_range.beg,
            pbc_ext);

        com.x = deperiodizef(com.x, pbc_ext.x * 0.5f, pbc_ext.x);
        com.y = deperiodizef(com.y, pbc_ext.y * 0.5f, pbc_ext.y);
        com.z = deperiodizef(com.z, pbc_ext.z * 0.5f, pbc_ext.z);

        for (int64_t j = atom_range.beg; j < atom_range.end; ++j) {
            mol->atom.x[j] = deperiodizef(mol->atom.x[j], com.x, pbc_ext.x);
            mol->atom.y[j] = deperiodizef(mol->atom.y[j], com.y, pbc_ext.y);
            mol->atom.z[j] = deperiodizef(mol->atom.z[j], com.z, pbc_ext.z);
        }

        if (residue_com) {
            residue_com[i] = com;
        }
    }

    // Chains are a bit problematic and can span more than half the extent of a simulation box
    // The question is then what is the correct placement of it with respect to the box?
    // Ideally we don't want to break up chains since they will certainly hold secondary structures
    // Which span multiple residues

    for (int64_t j = 0; j < mol->chain.count; ++j) {
        md_range_t res_range = mol->chain.residue_range[j];
        vec3_t chain_com = md_util_compute_com_ortho(residue_com + res_range.beg, 0, res_range.end - res_range.beg, pbc_ext);

        {
            vec3_t dp;
            dp.x = deperiodizef(chain_com.x, pbc_ext.x * 0.5f, pbc_ext.x);
            dp.y = deperiodizef(chain_com.y, pbc_ext.y * 0.5f, pbc_ext.y);
            dp.z = deperiodizef(chain_com.z, pbc_ext.z * 0.5f, pbc_ext.z);

            vec3_t d = vec3_sub(dp, chain_com);
            if (vec3_dot(d, d) > 0) {
                chain_com = dp;
                for (int64_t k = res_range.beg; k < res_range.end; ++k) {
                    residue_com[k] = vec3_add(residue_com[k], d);
                }
            }
        }

        // Pick a a residue reference which is closest to the COM
        // Then we deperiodize from this residue forward and backward each residue with respect to the previous.
        int ref_idx = res_range.beg;
        float ref_d2 = vec3_distance_squared(residue_com[ref_idx], chain_com);
        for (int i = res_range.beg + 1; i < res_range.end; ++i) {
            float d2 = vec3_distance_squared(residue_com[i], chain_com);
            if (d2 < ref_d2) {
                ref_d2 = d2;
                ref_idx = i;
            }
        }

        // First we deperiodize the reference to ensure it is within the bounds of the simulation box
        // Ensure that the residue is within the period of the chain
        {
            vec3_t res_com = residue_com[ref_idx];
            res_com.x = deperiodizef(res_com.x, chain_com.x, pbc_ext.x);
            res_com.y = deperiodizef(res_com.y, chain_com.y, pbc_ext.y);
            res_com.z = deperiodizef(res_com.z, chain_com.z, pbc_ext.z);

            vec3_t d = vec3_sub(res_com, residue_com[ref_idx]);

            if (vec3_dot(d,d) > 0) {
                residue_com[ref_idx] = res_com;
                md_range_t atom_range = mol->residue.atom_range[ref_idx];
                for (int64_t k = atom_range.beg; k < atom_range.end; ++k) {
                    mol->atom.x[k] += d.x;
                    mol->atom.y[k] += d.y;
                    mol->atom.z[k] += d.z;
                }
            }
        }

        // Propagate backwards
        for (int i = ref_idx - 1; i >= res_range.beg; --i) {
            vec3_t res_com = residue_com[i];
            vec3_t ref_com = residue_com[i+1];

            // Ensure that the residue is within the period of the chain
            res_com.x = deperiodizef(res_com.x, ref_com.x, pbc_ext.x);
            res_com.y = deperiodizef(res_com.y, ref_com.y, pbc_ext.y);
            res_com.z = deperiodizef(res_com.z, ref_com.z, pbc_ext.z);

            vec3_t d = vec3_sub(res_com, residue_com[i]);

            if (vec3_dot(d,d) > 0) {
                residue_com[i] = res_com;
                md_range_t atom_range = mol->residue.atom_range[i];
                for (int64_t k = atom_range.beg; k < atom_range.end; ++k) {
                    mol->atom.x[k] += d.x;
                    mol->atom.y[k] += d.y;
                    mol->atom.z[k] += d.z;
                }
            }    
        }
        // Propagate forwards
        for (int i = ref_idx + 1; i < res_range.end; ++i) {
            vec3_t res_com = residue_com[i];
            vec3_t ref_com = residue_com[i-1];

            // Ensure that the residue is within the period of the chain
            res_com.x = deperiodizef(res_com.x, ref_com.x, pbc_ext.x);
            res_com.y = deperiodizef(res_com.y, ref_com.y, pbc_ext.y);
            res_com.z = deperiodizef(res_com.z, ref_com.z, pbc_ext.z);

            vec3_t d = vec3_sub(res_com, residue_com[i]);

            if (vec3_dot(d,d) > 0) {
                residue_com[i] = res_com;
                md_range_t atom_range = mol->residue.atom_range[i];
                for (int64_t k = atom_range.beg; k < atom_range.end; ++k) {
                    mol->atom.x[k] += d.x;
                    mol->atom.y[k] += d.y;
                    mol->atom.z[k] += d.z;
                }
            }    
        }
    }

    if (residue_com) {
        md_array_free(residue_com, default_allocator);
    }
    return true;
}
*/

mat3_t md_util_compute_optimal_rotation(const md_vec3_soa_t coord[2], const vec3_t com[2], const float* w, int64_t count) {
    if (count < 1) {
        return mat3_ident();
    }

    mat3_t cov_mat = {0};
    if (w) {
        cov_mat = mat3_weighted_cross_covariance_matrix(coord[0].x, coord[0].y, coord[0].z, coord[1].x, coord[1].y, coord[1].z, w, com[0], com[1], count);
    } else {
        cov_mat = mat3_cross_covariance_matrix(coord[0].x, coord[0].y, coord[0].z, coord[1].x, coord[1].y, coord[1].z, com[0], com[1], count);
    }

    return mat3_extract_rotation(cov_mat);
}

double md_util_compute_rmsd(const md_vec3_soa_t coord[2], const vec3_t com[2], const float* w, int64_t count) {
    const mat3_t R = md_util_compute_optimal_rotation(coord, com, w, count);
    double d_sum = 0;
    double w_sum = 0;
    for (int64_t i = 0; i < count; ++i) {
        vec3_t u = {coord[0].x[i] - com[0].x, coord[0].y[i] - com[0].y, coord[0].z[i] - com[0].z};
        vec3_t v = {coord[1].x[i] - com[1].x, coord[1].y[i] - com[1].y, coord[1].z[i] - com[1].z};
        vec3_t vp = mat3_mul_vec3(R, v);
        vec3_t d = vec3_sub(u, vp);
        float weight = w ? w[i] : 1.0f;
        d_sum += weight * vec3_dot(d, d);
        w_sum += weight;
    }

    return sqrt(d_sum / w_sum);
}

bool md_util_linear_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[2], int64_t count, vec3_t pbc_ext, float t) {
    const bool use_pbc = !vec3_equal(pbc_ext, (vec3_t){0,0,0});
    t = CLAMP(t, 0.0f, 1.0f);

    if (use_pbc) {
        md_simd_f32_t box_ext_x = md_simd_set1_f32(pbc_ext.x);
        md_simd_f32_t box_ext_y = md_simd_set1_f32(pbc_ext.y);
        md_simd_f32_t box_ext_z = md_simd_set1_f32(pbc_ext.z);

        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_f32_width) * md_simd_f32_width;
        for (; i < simd_count; i += md_simd_f32_width) {
            md_simd_f32_t x0 = md_simd_load_f32(src_coord[0].x + i);
            md_simd_f32_t y0 = md_simd_load_f32(src_coord[0].y + i);
            md_simd_f32_t z0 = md_simd_load_f32(src_coord[0].z + i);

            md_simd_f32_t x1 = md_simd_load_f32(src_coord[1].x + i);
            md_simd_f32_t y1 = md_simd_load_f32(src_coord[1].y + i);
            md_simd_f32_t z1 = md_simd_load_f32(src_coord[1].z + i);

            x1 = simd_deperiodize(x1, x0, box_ext_x);
            y1 = simd_deperiodize(y1, y0, box_ext_y);
            z1 = simd_deperiodize(z1, z0, box_ext_z);

            md_simd_f32_t x = md_simd_lerp(x0, x1, t);
            md_simd_f32_t y = md_simd_lerp(y0, y1, t);
            md_simd_f32_t z = md_simd_lerp(z0, z1, t);

            md_simd_store(dst_coord.x + i, x);
            md_simd_store(dst_coord.y + i, y);
            md_simd_store(dst_coord.z + i, z);
        }

        // Do the rest
        const vec4_t pbc_ext4 = vec4_from_vec3(pbc_ext, 0.0f);
        for (; i < count; ++i) {
            vec4_t src[2] = {
                {src_coord[0].x[i], src_coord[0].y[i], src_coord[0].z[i], 1},
                {src_coord[1].x[i], src_coord[1].y[i], src_coord[1].z[i], 1},
            };

            src[1] = vec4_deperiodize(src[1], src[0], pbc_ext4);
            const vec4_t coord = vec4_lerp(src[0], src[1], t);

            dst_coord.x[i] = coord.x;
            dst_coord.y[i] = coord.y;
            dst_coord.z[i] = coord.z;
        }
    } else {
        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_f32_width) * md_simd_f32_width;
        for (; i < simd_count; i += md_simd_f32_width) {
            md_simd_f32_t x0 = md_simd_load_f32(src_coord[0].x + i);
            md_simd_f32_t y0 = md_simd_load_f32(src_coord[0].y + i);
            md_simd_f32_t z0 = md_simd_load_f32(src_coord[0].z + i);

            md_simd_f32_t x1 = md_simd_load_f32(src_coord[1].x + i);
            md_simd_f32_t y1 = md_simd_load_f32(src_coord[1].y + i);
            md_simd_f32_t z1 = md_simd_load_f32(src_coord[1].z + i);

            md_simd_f32_t x = md_simd_lerp(x0, x1, t);
            md_simd_f32_t y = md_simd_lerp(y0, y1, t);
            md_simd_f32_t z = md_simd_lerp(z0, z1, t);

            md_simd_store(dst_coord.x + i, x);
            md_simd_store(dst_coord.y + i, y);
            md_simd_store(dst_coord.z + i, z);
        }

        // Do the rest
        for (; i < count; ++i) {
            vec4_t src[2] = {
                {src_coord[0].x[i], src_coord[0].y[i], src_coord[0].z[i], 1},
                {src_coord[1].x[i], src_coord[1].y[i], src_coord[1].z[i], 1},
            };

            const vec4_t coord = vec4_lerp(src[0], src[1], t);

            dst_coord.x[i] = coord.x;
            dst_coord.y[i] = coord.y;
            dst_coord.z[i] = coord.z;
        }
    }

    return true;
}

bool md_util_cubic_spline_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[4], int64_t count, vec3_t pbc_ext, float t, float s) {
    const bool use_pbc = !vec3_equal(pbc_ext, (vec3_t){0,0,0});
    t = CLAMP(t, 0.0f, 1.0f);
    s = CLAMP(s, 0.0f, 1.0f);

    if (use_pbc) {
        md_simd_f32_t box_ext_x = md_simd_set1_f32(pbc_ext.x);
        md_simd_f32_t box_ext_y = md_simd_set1_f32(pbc_ext.y);
        md_simd_f32_t box_ext_z = md_simd_set1_f32(pbc_ext.z);

        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_f32_width) * md_simd_f32_width;
        for (; i < simd_count; i += md_simd_f32_width) {
            md_simd_f32_t x0 = md_simd_load_f32(src_coord[0].x + i);
            md_simd_f32_t y0 = md_simd_load_f32(src_coord[0].y + i);
            md_simd_f32_t z0 = md_simd_load_f32(src_coord[0].z + i);

            md_simd_f32_t x1 = md_simd_load_f32(src_coord[1].x + i);
            md_simd_f32_t y1 = md_simd_load_f32(src_coord[1].y + i);
            md_simd_f32_t z1 = md_simd_load_f32(src_coord[1].z + i);

            md_simd_f32_t x2 = md_simd_load_f32(src_coord[2].x + i);
            md_simd_f32_t y2 = md_simd_load_f32(src_coord[2].y + i);
            md_simd_f32_t z2 = md_simd_load_f32(src_coord[2].z + i);

            md_simd_f32_t x3 = md_simd_load_f32(src_coord[3].x + i);
            md_simd_f32_t y3 = md_simd_load_f32(src_coord[3].y + i);
            md_simd_f32_t z3 = md_simd_load_f32(src_coord[3].z + i);

            x0 = simd_deperiodize(x0, x1, box_ext_x);
            x2 = simd_deperiodize(x2, x1, box_ext_x);
            x3 = simd_deperiodize(x3, x2, box_ext_x);

            y0 = simd_deperiodize(y0, y1, box_ext_y);
            y2 = simd_deperiodize(y2, y1, box_ext_y);
            y3 = simd_deperiodize(y3, y2, box_ext_y);

            z0 = simd_deperiodize(z0, z1, box_ext_z);
            z2 = simd_deperiodize(z2, z1, box_ext_z);
            z3 = simd_deperiodize(z3, z2, box_ext_z);

            md_simd_f32_t x = md_simd_cubic_spline(x0, x1, x2, x3, md_simd_set1_f32(t), md_simd_set1_f32(s));
            md_simd_f32_t y = md_simd_cubic_spline(y0, y1, y2, y3, md_simd_set1_f32(t), md_simd_set1_f32(s));
            md_simd_f32_t z = md_simd_cubic_spline(z0, z1, z2, z3, md_simd_set1_f32(t), md_simd_set1_f32(s));

            md_simd_store(dst_coord.x + i, x);
            md_simd_store(dst_coord.y + i, y);
            md_simd_store(dst_coord.z + i, z);
        }

        // Do the rest
        const vec4_t pbc_ext4 = vec4_from_vec3(pbc_ext, 0.0f);
        for (; i < count; ++i) {
            vec4_t src[4] = {
                {src_coord[0].x[i], src_coord[0].y[i], src_coord[0].z[i], 1},
                {src_coord[1].x[i], src_coord[1].y[i], src_coord[1].z[i], 1},
                {src_coord[2].x[i], src_coord[2].y[i], src_coord[2].z[i], 1},
                {src_coord[3].x[i], src_coord[3].y[i], src_coord[3].z[i], 1},
            };

            src[0] = vec4_deperiodize(src[0], src[1], pbc_ext4);
            src[2] = vec4_deperiodize(src[2], src[1], pbc_ext4);
            src[3] = vec4_deperiodize(src[3], src[2], pbc_ext4);

            const vec4_t coord = vec4_cubic_spline(src[0], src[1], src[2], src[3], t, s);

            dst_coord.x[i] = coord.x;
            dst_coord.y[i] = coord.y;
            dst_coord.z[i] = coord.z;
        }
    } else {
        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_f32_width) * md_simd_f32_width;
        for (; i < simd_count; i += md_simd_f32_width) {
            md_simd_f32_t x0 = md_simd_load_f32(src_coord[0].x + i);
            md_simd_f32_t y0 = md_simd_load_f32(src_coord[0].y + i);
            md_simd_f32_t z0 = md_simd_load_f32(src_coord[0].z + i);

            md_simd_f32_t x1 = md_simd_load_f32(src_coord[1].x + i);
            md_simd_f32_t y1 = md_simd_load_f32(src_coord[1].y + i);
            md_simd_f32_t z1 = md_simd_load_f32(src_coord[1].z + i);

            md_simd_f32_t x2 = md_simd_load_f32(src_coord[2].x + i);
            md_simd_f32_t y2 = md_simd_load_f32(src_coord[2].y + i);
            md_simd_f32_t z2 = md_simd_load_f32(src_coord[2].z + i);

            md_simd_f32_t x3 = md_simd_load_f32(src_coord[3].x + i);
            md_simd_f32_t y3 = md_simd_load_f32(src_coord[3].y + i);
            md_simd_f32_t z3 = md_simd_load_f32(src_coord[3].z + i);

            md_simd_f32_t x = md_simd_cubic_spline(x0, x1, x2, x3, md_simd_set1_f32(t), md_simd_set1_f32(s));
            md_simd_f32_t y = md_simd_cubic_spline(y0, y1, y2, y3, md_simd_set1_f32(t), md_simd_set1_f32(s));
            md_simd_f32_t z = md_simd_cubic_spline(z0, z1, z2, z3, md_simd_set1_f32(t), md_simd_set1_f32(s));

            md_simd_store_f32(dst_coord.x + i, x);
            md_simd_store_f32(dst_coord.y + i, y);
            md_simd_store_f32(dst_coord.z + i, z);
        }

        // Do the rest
        for (; i < count; ++i) {
            vec4_t src[4] = {
                {src_coord[0].x[i], src_coord[0].y[i], src_coord[0].z[i], 1},
                {src_coord[1].x[i], src_coord[1].y[i], src_coord[1].z[i], 1},
                {src_coord[2].x[i], src_coord[2].y[i], src_coord[2].z[i], 1},
                {src_coord[3].x[i], src_coord[3].y[i], src_coord[3].z[i], 1},
            };

            vec4_t coord = vec4_cubic_spline(src[0], src[1], src[2], src[3], t, s);

            dst_coord.x[i] = coord.x;
            dst_coord.y[i] = coord.y;
            dst_coord.z[i] = coord.z;
        }
    }

    return true;
}

static inline bool ranges_overlap(md_range_t a, md_range_t b) {
    return (a.beg < b.end && b.beg < a.end);
}

static inline void commit_backbone(md_backbone_atoms_t* backbone, md_range_t res_range, md_molecule_t* mol, md_allocator_i* alloc) {
    md_range_t  bb_range = {(int32_t)md_array_size(mol->backbone.atoms), (int32_t)(md_array_size(mol->backbone.atoms) + md_array_size(backbone))};
    md_array_push_array(mol->backbone.atoms, backbone, md_array_size(backbone), alloc);
    md_array_push(mol->backbone.range, bb_range, alloc);

    for (int32_t i = res_range.beg; i < res_range.end; ++i) {
        md_array_push(mol->backbone.residue_idx, i, alloc);
    }
}

// Try to fill in missing fields for molecule struct
// (Labels)   -> Elements
// (Elements) -> Mass & Radius
// (Coordinates & Elements) -> Covalent Bonds
bool md_util_postprocess_molecule(struct md_molecule_t* mol, struct md_allocator_i* alloc, md_util_postprocess_flags_t flags) {
    ASSERT(mol);
    ASSERT(alloc);

    if (mol->atom.vx == 0) {
        md_array_resize(mol->atom.vx, mol->atom.count, alloc);
        MEMSET(mol->atom.vx, 0, md_array_size(mol->atom.vx) * sizeof(*mol->atom.vx));
    }
    if (mol->atom.vy == 0) {
        md_array_resize(mol->atom.vy, mol->atom.count, alloc);
        MEMSET(mol->atom.vy, 0, md_array_size(mol->atom.vy) * sizeof(*mol->atom.vy));
    }
    if (mol->atom.vz == 0) {
        md_array_resize(mol->atom.vz, mol->atom.count, alloc);
        MEMSET(mol->atom.vz, 0, md_array_size(mol->atom.vz) * sizeof(*mol->atom.vz));
    }

    if (mol->atom.flags == 0) {
        md_array_resize(mol->atom.flags, mol->atom.count, alloc);
        MEMSET(mol->atom.flags, 0, md_array_size(mol->atom.flags) * sizeof(*mol->atom.flags));
    }

    if (flags & MD_UTIL_POSTPROCESS_ELEMENT_BIT) {
        if (mol->atom.element == 0) {
            md_array_resize(mol->atom.element, mol->atom.count, alloc);
            MEMSET(mol->atom.element, 0, mol->atom.count * sizeof(*mol->atom.element));
        }
        md_util_element_decode(mol->atom.element, mol->atom.count, mol);
    }

    if (flags & MD_UTIL_POSTPROCESS_RADIUS_BIT) {
        if (mol->atom.radius == 0)  md_array_resize(mol->atom.radius, mol->atom.count, alloc);
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            mol->atom.radius[i] = mol->atom.element ? md_util_element_vdw_radius(mol->atom.element[i]) : 1.0f;
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_MASS_BIT) {
        if (mol->atom.mass == 0) md_array_resize(mol->atom.mass, mol->atom.count, alloc);
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            mol->atom.mass[i] = mol->atom.element ? md_util_element_atomic_mass(mol->atom.element[i]) : 1.0f;
        }
    }
   
    if (flags & MD_UTIL_POSTPROCESS_COVALENT_BIT) {
        if (mol->covalent.count == 0) {
            const vec3_t pbc_ext = md_util_compute_unit_cell_extent(mol->coord_frame);
            md_util_compute_covalent_bonds(&mol->covalent, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.element, mol->atom.residue_idx, mol->atom.count, pbc_ext, alloc);
            md_util_compute_structures(&mol->covalent.structures, mol->atom.count, mol->covalent.bond, mol->covalent.count, alloc);
            md_util_compute_rings(&mol->covalent.rings, mol->atom.count, mol->covalent.bond, mol->covalent.count, alloc);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_VALENCE_BIT) {
        if (mol->atom.valence == 0 && mol->covalent.count > 0) {
            md_array_resize(mol->atom.valence, mol->atom.count, alloc);
            md_util_compute_atom_valence(mol->atom.valence, mol->atom.count, mol->covalent.bond, mol->covalent.count);    
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_CHAINS_BIT) {
        if (mol->chain.count == 0 && mol->residue.count > 0 && mol->covalent.count > 0) {
            md_util_compute_chain_data(&mol->chain, mol->atom.residue_idx, mol->atom.count, mol->covalent.bond, mol->covalent.count, alloc);

            if (mol->chain.count) {
                md_array_resize(mol->atom.chain_idx, mol->atom.count, alloc);
                for (int i = 0; i < (int)mol->atom.count; ++i) {
                    mol->atom.chain_idx[i] = -1;
                }

                // Update references from atoms and residues to chains
                for (md_chain_idx_t c = 0; c < (md_chain_idx_t)mol->chain.count; ++c) {
                    md_range_t atom_range = mol->chain.atom_range[c];
                    for (int i = atom_range.beg; i < atom_range.end; ++i) {
                        mol->atom.chain_idx[i] = c;
                    }
                }
            }
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_BACKBONE_BIT) {
        if (mol->chain.count) {
            // Compute backbone data
            // 
            // @NOTE: We should only attempt to compute backbone data for valid residues (e.g. amino acids / dna)
            // Backbones are not directly tied to chains and therefore we cannot use chains as a 1:1 mapping for the backbones.
            // We look within the chains and see if we can find consecutive ranges which form backbones.

            static const int64_t MIN_BACKBONE_LENGTH = 3;
            md_backbone_atoms_t* backbone = 0;
        
            for (int64_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
                for (int64_t res_idx = mol->chain.residue_range[chain_idx].beg; res_idx < mol->chain.residue_range[chain_idx].end; ++res_idx) {
                    md_backbone_atoms_t atoms;
                    if (md_util_backbone_atoms_extract_from_residue_idx(&atoms, (int32_t)res_idx, mol)) {
                        md_array_push(backbone, atoms, default_allocator);
                    } else {
                        if (md_array_size(backbone) >= MIN_BACKBONE_LENGTH) {
                            // Commit the backbone
                            md_range_t res_range = {(int32_t)(res_idx - md_array_size(backbone)), (int32_t)res_idx};
                            commit_backbone(backbone, res_range, mol, alloc);
                        }
                        md_array_shrink(backbone, 0);
                    }
                }
                // Possibly commit remainder of the chain
                if (md_array_size(backbone) >= MIN_BACKBONE_LENGTH) {
                    int64_t res_idx = mol->chain.residue_range[chain_idx].end;
                    md_range_t res_range = {(int32_t)(res_idx - md_array_size(backbone)), (int32_t)res_idx};
                    commit_backbone(backbone, res_range, mol, alloc);
                }
                md_array_shrink(backbone, 0);
            }

            md_array_free(backbone, default_allocator);

            mol->backbone.range_count = md_array_size(mol->backbone.range);
            mol->backbone.count = md_array_size(mol->backbone.atoms);

            md_array_resize(mol->backbone.angle, mol->backbone.count, alloc);
            md_array_resize(mol->backbone.secondary_structure, mol->backbone.count, alloc);
            md_array_resize(mol->backbone.ramachandran_type, mol->backbone.count, alloc);

            if (mol->backbone.count > 0) {
                md_util_backbone_angles_compute(mol->backbone.angle, mol->backbone.count, mol);
                md_util_backbone_secondary_structure_compute(mol->backbone.secondary_structure, mol->backbone.count, mol);
                md_util_backbone_ramachandran_classify(mol->backbone.ramachandran_type, mol->backbone.count, mol);
            }
        }
    }

    return true;
}

/*
This is blatantly stolen and modified from Arseny Kapoulkines Mesh optimizer
https://github.com/zeux/meshoptimizer/

MIT License
Copyright (c) 2016-2022 Arseny Kapoulkine
*/

// "Insert" two 0 bits after each of the 10 low bits of x
static inline uint32_t part1By2(uint32_t x) {
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x << 8)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x << 4)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x << 2)) & 0x09249249;  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return x;
}

static void computeOrder(uint32_t* result, const float* vx, const float* vy, const float* vz, size_t count, size_t byte_stride) {
    float minv[3] = {+FLT_MAX, +FLT_MAX, +FLT_MAX};
    float maxv[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    byte_stride = MAX(byte_stride, sizeof(float));
    for (size_t i = 0; i < count; ++i) {
        float x = *(float*)((const char*)vx + byte_stride * i);
        float y = *(float*)((const char*)vy + byte_stride * i);
        float z = *(float*)((const char*)vz + byte_stride * i);

        minv[0] = MIN(minv[0], x);
        maxv[0] = MAX(maxv[0], x);

        minv[1] = MIN(minv[1], y);
        maxv[1] = MAX(maxv[1], y);

        minv[2] = MIN(minv[2], z);
        maxv[2] = MAX(maxv[2], z);
    }

    float extent = 0.f;

    extent = (maxv[0] - minv[0]) < extent ? extent : (maxv[0] - minv[0]);
    extent = (maxv[1] - minv[1]) < extent ? extent : (maxv[1] - minv[1]);
    extent = (maxv[2] - minv[2]) < extent ? extent : (maxv[2] - minv[2]);

    float scale = extent == 0 ? 0.f : 1.f / extent;

    // generate Morton order based on the position inside a unit cube
    for (size_t i = 0; i < count; ++i) {
        int x = (int)((vx[i] - minv[0]) * scale * 1023.f + 0.5f);
        int y = (int)((vy[i] - minv[1]) * scale * 1023.f + 0.5f);
        int z = (int)((vz[i] - minv[2]) * scale * 1023.f + 0.5f);

        result[i] = part1By2(x) | (part1By2(y) << 1) | (part1By2(z) << 2);
    }
}

static void computeHistogram(uint32_t hist[1024][3], const uint32_t* data, size_t count) {
    // compute 3 10-bit histograms in parallel
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = data[i];

        hist[(id >> 0)  & 1023][0]++;
        hist[(id >> 10) & 1023][1]++;
        hist[(id >> 20) & 1023][2]++;
    }

    uint32_t sumx = 0, sumy = 0, sumz = 0;

    // replace histogram data with prefix histogram sums in-place
    for (int i = 0; i < 1024; ++i) {
        uint32_t hx = hist[i][0], hy = hist[i][1], hz = hist[i][2];

        hist[i][0] = sumx;
        hist[i][1] = sumy;
        hist[i][2] = sumz;

        sumx += hx;
        sumy += hy;
        sumz += hz;
    }

    ASSERT(sumx == count && sumy == count && sumz == count);
}

static void radixPass(uint32_t* destination, const uint32_t* source, const uint32_t* keys, size_t count, uint32_t hist[1024][3], int pass) {
    int bitoff = pass * 10;
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = (keys[source[i]] >> bitoff) & 1023;
        destination[hist[id][pass]++] = source[i];
    }
}

void md_util_spatial_sort_soa(uint32_t* source, const float* x, const float* y, const float* z, int64_t count) {
    if (!source || !z || !y || !z || count <= 0) return;

    md_allocator_i* alloc = default_allocator;

    uint32_t* keys = md_alloc(alloc, sizeof(uint32_t) * count);
    computeOrder(keys, x, y, z, count, 0);

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram(hist, keys, count);

    uint32_t* scratch = md_alloc(alloc, sizeof(uint32_t) * count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass radix sort computes the resulting order into scratch
    radixPass(source, scratch, keys, count, hist, 0);
    radixPass(scratch, source, keys, count, hist, 1);
    radixPass(source, scratch, keys, count, hist, 2);

    md_free(alloc, scratch, sizeof(uint32_t) * count);
    md_free(alloc, keys,    sizeof(uint32_t) * count);
}

void md_util_spatial_sort(uint32_t* source, const vec3_t* xyz, int64_t count) {
    if (!source || !xyz || count <= 0) return;

    md_allocator_i* alloc = default_allocator;

    uint32_t* keys = md_alloc(alloc, sizeof(uint32_t) * count);
    const float* base = (const float*)xyz;
    const float* x = base + 0;
    const float* y = base + 1;
    const float* z = base + 2;
    computeOrder(keys, x, y, z, count, sizeof(vec3_t));

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram(hist, keys, count);

    uint32_t* scratch = md_alloc(alloc, sizeof(uint32_t) * count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass radix sort computes the resulting order into scratch
    radixPass(source, scratch, keys, count, hist, 0);
    radixPass(scratch, source, keys, count, hist, 1);
    radixPass(source, scratch, keys, count, hist, 2);

    md_free(alloc, scratch, sizeof(uint32_t) * count);
    md_free(alloc, keys,    sizeof(uint32_t) * count);
}

#ifdef __cplusplus
}
#endif
