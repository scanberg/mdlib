#include <md_util.h>

#include <md_trajectory.h>
#include <md_molecule.h>

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

#include <math.h>
#include <string.h>
#include <float.h>

#define bake(str) {str, sizeof(str)-1}

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

#define RGBA(r,g,b,a) ( ((a & 255) << 24) | ((b & 255) << 16) | ((g & 255) << 8) | (r & 255) )
#define RGB(r,g,b) RGBA(r,g,b,255)

// Based on this with some values modified
// Some stem from Jmol, some from Rasmol
// http://jmol.sourceforge.net/jscolors/
static uint32_t element_cpk_colors[] = {
    0xFFFF00FF, 0xFFFFFFFF, 0xFFFFFFD9, 0xFF2222B2, 0xFF00FFC2, 0xFFB5B5FF, 0xFFB0B0B0, 0xFFFF8F8F, 0xFF0000F0, 0xFF50E090, 0xFFF5E3B3, 0xFFF25CAB,
    0xFF00FF8A, 0xFF908080, 0xFFA0C8F0, 0xFF00A5FF, 0xFF32C8FF, 0xFF1FF01F, 0xFFE3D180, 0xFFD4408F, 0xFF908080, 0xFFE6E6E6, 0xFF908080, 0xFFABA6A6,
    0xFF908080, 0xFF908080, 0xFF3366E0, 0xFFA090F0, 0xFF2A2AA5, 0xFF2A2AA5, 0xFF2A2AA5, 0xFF8F8FC2, 0xFF8F8F66, 0xFFE380BD, 0xFF00A1FF, 0xFF2A2AA5,
    0xFFD1B85C, 0xFFB02E70, 0xFF00FF00, 0xFFFFFF94, 0xFFE0E094, 0xFFC9C273, 0xFFB5B554, 0xFF9E9E3B, 0xFF8F8F24, 0xFF8C7D0A, 0xFF856900, 0xFF908080,
    0xFF8FD9FF, 0xFF7375A6, 0xFF808066, 0xFFB5639E, 0xFF007AD4, 0xFF940094, 0xFFB09E42, 0xFF8F1757, 0xFF00A5FF, 0xFFFFD470, 0xFFC7FFFF, 0xFFC7FFD9,
    0xFFC7FFC7, 0xFFC7FFA3, 0xFFC7FF8F, 0xFFC7FF61, 0xFFC7FF45, 0xFFC7FF30, 0xFFC7FF1F, 0xFF9CFF00, 0xFF75E600, 0xFF52D400, 0xFF38BF00, 0xFF24AB00,
    0xFFFFC24D, 0xFFFFA64D, 0xFFD69421, 0xFFAB7D26, 0xFF966626, 0xFF875417, 0xFFE0D0D0, 0xFF23D1FF, 0xFFD0B8B8, 0xFF4D54A6, 0xFF615957, 0xFFB54F9E,
    0xFF005CAB, 0xFF454F75, 0xFF968242, 0xFF660042, 0xFF007D00, 0xFFFAAB70, 0xFFFFBA00, 0xFFFFA100, 0xFFFF8F00, 0xFFFF8000, 0xFFFF6B00, 0xFFF25C54,
    0xFFE35C78, 0xFFE34F8A, 0xFFD436A1, 0xFFD41FB3, 0xFFBA1FB3, 0xFFA60DB3, 0xFF870DBD, 0xFF6600C7, 0xFF5900CC, 0xFF4F00D1, 0xFF4500D9, 0xFF3800E0,
    0xFF2E00E6, 0xFF2600EB, 0xFF2200F0, 0xFF2000F6, 0xFF1E00F8, 0xFF1C00FA, 0xFF1A00FC, 0xFF1800FD, 0xFF1600FE, 0xFF1400FF, 0xFF1200FF
};

// This has been filled with some entries found in molstar (github.com/molstar)
static const char* amino_acids[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU",
    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASC", "GLX", "XLE", "HISE",
    
    // CCD
    "UNK",
    "MSE", "SEP", "TPO", "PTR", "PCA",

    // Charmm
    "HSD", "HSP", "LSN", "ASPP", "GLUP",

    // Amber
    "HID", "HIE", "HIP", "LYN", "ASH", "GLH",
};

static const char* rna[] = {"A", "C", "T", "G", "I", "U", "N"};
static const char* dna[] = {"DA", "DC", "DG", "DT"};

static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4", "W", "DOD", "D30" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

static inline int64_t find_str_in_array(str_t str, const char* arr[], int64_t arr_len) {
    for (int64_t i = 0; i < arr_len; ++i) {
        if (str_equal_cstr(str, arr[i])) return i;
    }
    return -1;
}


// Trim whitespace, digits and 'X's
static inline str_t trim_label(str_t lbl) {
    const char* beg = str_beg(lbl);
    const char* end = str_end(lbl);
    const char* c = beg;
    while (c < end && *c && (is_digit(*c) || is_whitespace(*c) || *c == 'x' || *c == 'X')) ++c;
    beg = c;
    while (c < end && is_alpha(*c) && !(*c == 'x' || *c == 'X')) ++c;
    end = c;

    return (str_t) { .ptr = beg, .len = end-beg };
}

static inline md_simd_f32_t simd_deperiodize(md_simd_f32_t x, md_simd_f32_t r, md_simd_f32_t period) {
    md_simd_f32_t d = md_simd_sub(x, r);
    md_simd_f32_t dx = md_simd_div(d, period);
    dx = md_simd_sub(dx, md_simd_round(dx));
    return md_simd_add(r, md_simd_mul(dx, period));
}

bool md_util_resname_rna(str_t str) {
    str = trim_label(str);
    return find_str_in_array(str, rna, ARRAY_SIZE(rna)) != -1;
}
    
bool md_util_resname_dna(str_t str) {
    str = trim_label(str);
    return find_str_in_array(str, dna, ARRAY_SIZE(dna)) != -1;
}

bool md_util_resname_nucleic_acid(str_t str) {
    str = trim_label(str);
    return find_str_in_array(str, rna, ARRAY_SIZE(rna)) > -1 || find_str_in_array(str, dna, ARRAY_SIZE(dna)) > -1;
}
    
bool md_util_resname_acidic(str_t str) {
    return find_str_in_array(str, acidic, ARRAY_SIZE(acidic)) != -1;
}

bool md_util_resname_basic(str_t str) {
    return find_str_in_array(str, basic, ARRAY_SIZE(basic)) != -1;
}
    
bool md_util_resname_neutral(str_t str) {
    return find_str_in_array(str, neutral, ARRAY_SIZE(neutral)) != -1;
}
    
bool md_util_resname_water(str_t str) {
    return find_str_in_array(str, water, ARRAY_SIZE(water)) != -1;
}
    
bool md_util_resname_hydrophobic(str_t str) {
    return find_str_in_array(str, hydrophobic, ARRAY_SIZE(hydrophobic)) != -1;
}
    
bool md_util_resname_amino_acid(str_t str) {
    return find_str_in_array(str, amino_acids, ARRAY_SIZE(amino_acids)) != -1;
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
        lbl = trim_label(lbl);
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

#define MIN_RES_LEN 4
#define MAX_RES_LEN 25

bool md_util_element_guess(md_element_t element[], int64_t capacity, const struct md_molecule_t* mol) {
    ASSERT(capacity >= 0);
    ASSERT(mol);
    ASSERT(mol->atom.count >= 0);

    // @TODO, PERF: Iterate over residues and check if the entire residue is an amino acid
    const int64_t count = MIN(capacity, mol->atom.count);
    for (int64_t i = 0; i < count; ++i) {
        if (element[i] != 0) continue;

        str_t original = {mol->atom.name[i].buf, mol->atom.name[i].len};

        // Trim whitespace, digits and 'X's
        str_t name = trim_label(original);

        if (name.len > 0) {

            md_element_t elem = 0;
            if ((elem = md_util_element_lookup(name)) != 0) goto done;

            // If amino acid, try to deduce the element from that
            if (mol->atom.residue_idx) {
                const md_range_t res_range = mol->residue.atom_range[mol->atom.residue_idx[i]];
                const int res_len = res_range.end - res_range.beg;
                
                if (MIN_RES_LEN < res_len && res_len < MAX_RES_LEN && mol->residue.name) {
                    str_t resname = LBL_TO_STR(mol->residue.name[mol->atom.residue_idx[i]]);
                    if (md_util_resname_amino_acid(resname) ||
                        amino_acid_heuristic(mol->atom.name + res_range.beg, res_range.end - res_range.beg) ||
                        md_util_resname_dna(resname))
                    {
                        // EASY-PEASY, we just try to match against the first character
                        name.len = 1;
                        elem = md_util_element_lookup_ignore_case(name);
                        goto done;
                    }
                }
            }

            // Heuristic cases

            int num_alpha = 0;
            while (num_alpha < original.len && is_alpha(original.ptr[num_alpha])) ++num_alpha;
            
            int num_digits = 0;
            str_t digits = str_substr(original, num_alpha, -1);
            while (num_digits < digits.len && is_digit(digits.ptr[num_digits])) ++num_digits;

            if (str_equal_cstr(name, "HOH")) {
                elem = H;
                goto done;
            }

            // 2-3 letters + 1-2 digit (e.g. HO(H)[0-99]) usually means just look at the first letter
            if ((num_alpha == 2 || num_alpha == 3) && (num_digits == 1 || num_digits == 2)) {
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

    const float eps = 1.0e-2f;
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
    const uint32_t bb_bits = 2 | 4 | 8;
    if ((bits & bb_bits) == bb_bits) {
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

static inline bool zhang_skolnick_ss(const md_molecule_t* mol, md_range_t bb_range, int i, const float distances[3], float delta) {
    ASSERT(mol);
    ASSERT(mol->atom.x);
    ASSERT(mol->atom.y);
    ASSERT(mol->atom.z);
    ASSERT(mol->backbone.atoms);
    for (int j = MAX((int)bb_range.beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= (int)bb_range.end) continue;
            const int ca_j = mol->backbone.atoms[j].ca;
            const int ca_k = mol->backbone.atoms[j + k].ca;
            const vec3_t pos_j = {mol->atom.x[ca_j], mol->atom.y[ca_j], mol->atom.z[ca_j]};
            const vec3_t pos_k = {mol->atom.x[ca_k], mol->atom.y[ca_k], mol->atom.z[ca_k]};
            const float d = vec3_distance(pos_j, pos_k);
            if (fabsf(d - distances[k - 2]) > delta) {
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
        for (int32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i-1] != MD_SECONDARY_STRUCTURE_COIL && ss[i] != ss[i-1] && ss[i-1] == ss[i+1]) ss[i] = ss[i-1];
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MD_SECONDARY_STRUCTURE_COIL;
        for (int32_t i = range.beg + 1; i < range.end - 1; ++i) {
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

static md_index_data_t md_compute_atom_bond_indices(const md_bond_t bonds[], int64_t bond_count, int64_t atom_count, md_allocator_i* alloc) {
    md_index_data_t data = {0};
    md_array_resize(data.ranges, atom_count, alloc);
    MEMSET(data.ranges, 0, md_array_bytes(data.ranges));

    // This have length of 2 * bond_count (2 = one for each atom involved)
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
        data.ranges[i].beg =  data.ranges[i - 1].end;
        data.ranges[i].end += data.ranges[i].beg;
    }

    // Write edge indices to correct location
    for (int i = 0; i < (int)bond_count; ++i) {
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
        data.indices[idx_a] = i;
        data.indices[idx_b] = i;
    }

    md_array_free(local_offset, alloc);

    return data;
}

static md_index_data_t md_compute_connectivity(const md_bond_t bonds[], int64_t bond_count, int64_t atom_count, md_allocator_i* alloc) {
    md_index_data_t data = {0};
    md_array_resize(data.ranges, atom_count, alloc);
    MEMSET(data.ranges, 0, md_array_bytes(data.ranges));

    // This have length of 2 * bond_count (2 = one for direction)
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
        data.ranges[i].beg =  data.ranges[i - 1].end;
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

// If an entry exists in this table, it has order 2
static const char* intra_bond_order_table[] = {
    "HIS|CD2|CG",
    "HIS|CE1|ND1",
    "ARG|CZ|NH2",
    "PHE|CE1|CZ",
    "PHE|CD2|CE2",
    "PHE|CD1|CG",
    "TRP|CD1|CG",
    "TRP|CD2|CE2",
    "TRP|CE3|CZ3",
    "TRP|CH2|CZ2",
    "ASN|CG|OD1",
    "GLN|CD|OE1",
    "TYR|CD1|CG",
    "TYR|CD2|CE2",
    "TYR|CE1|CZ",
    "ASP|CG|OD1",
    "GLU|CD|OE1",

    "G|C8|N7",
    "G|C4|C5",
    "G|C2|N3",
    "G|C6|O6",
    "C|C4|N3",
    "C|C5|C6",
    "C|C2|O2",
    "A|C2|N3",
    "A|C6|N1",
    "A|C4|C5",
    "A|C8|N7",
    "U|C5|C6",
    "U|C2|O2",
    "U|C4|O4",

    "DG|C8|N7",
    "DG|C4|C5",
    "DG|C2|N3",
    "DG|C6|O6",
    "DC|C4|N3",
    "DC|C5|C6",
    "DC|C2|O2",
    "DA|C2|N3",
    "DA|C6|N1",
    "DA|C4|C5",
    "DA|C8|N7",
    "DT|C5|C6",
    "DT|C2|O2",
    "DT|C4|O4",
};

static inline int build_key(char* buf, str_t res, str_t a, str_t b) {
    int len = 0;
    while (res.len) {
        buf[len++] = *res.ptr++;
        res.len--;
    }
    buf[len++] = '|';
    while (a.len) {
        buf[len++] = *a.ptr++;
        a.len--;
    }
    buf[len++] = '|';
    while (b.len) {
        buf[len++] = *b.ptr++;
        b.len--;
    }
    buf[len] = '\0';
    return len;
}

// This is a c implementation of the tabular data found in Molstar (github.com/molstar)
md_array(md_order_t) md_util_compute_covalent_bond_order(const md_bond_t* bonds, int64_t bond_count, const md_label_t* atom_label, const md_residue_idx_t* atom_res_idx, const md_label_t* res_name, md_allocator_i* alloc) {
    ASSERT(bonds);
    ASSERT(alloc);

    md_array(md_order_t) order = 0;
    md_array_resize(order, bond_count, alloc);

    if (!atom_label || !atom_res_idx || !res_name) {
        MD_LOG_DEBUG("No atom label or atom residue indices or residue names were given, will default to covalent order of 1 for all bonds.");
        STATIC_ASSERT(sizeof(md_order_t) == 1, "Incorrect size of md_order_t");
        MEMSET(order, 1, md_array_bytes(order));
        return order;
    }
    
    for (int64_t i = 0; i < bond_count; ++i) {
        const int atom_a = bonds[i].idx[0];
        const int atom_b = bonds[i].idx[1];
        str_t atomname_a = LBL_TO_STR(atom_label[atom_a]);
        str_t atomname_b = LBL_TO_STR(atom_label[atom_b]);
        const int res_a = atom_res_idx[atom_a];
        const int res_b = atom_res_idx[atom_b];
        str_t resname_a = LBL_TO_STR(res_name[res_a]);
        str_t resname_b = LBL_TO_STR(res_name[res_b]);

        if (str_equal(resname_a, resname_b)) {
            if (strcmp(atomname_a.ptr, atomname_b.ptr) > 0) {
                str_t temp = atomname_a;
                atomname_a = atomname_b;
                atomname_b = temp;
            }
            // Intra
            if (md_util_resname_amino_acid(resname_a) && str_equal_cstr(atomname_a, "C") && str_equal_cstr(atomname_a, "O")) {
                order[i] = 2;
                continue;
            }
            
            char buf[32];
            int len = build_key(buf, resname_a, atomname_a, atomname_b);
            if (find_str_in_array((str_t){buf,len}, intra_bond_order_table, ARRAY_SIZE(intra_bond_order_table)) >= 0) {
                order[i] = 2;
            } else {
                order[i] = 1;
            }
        } else {
            // Inter
            if ( (str_equal_cstr(resname_a, "LYS") && str_equal_cstr(atomname_a, "CZ") && str_equal_cstr(resname_b, "RET") && str_equal_cstr(atomname_b, "C15")) ||
                 (str_equal_cstr(resname_b, "LYS") && str_equal_cstr(atomname_b, "CZ") && str_equal_cstr(resname_a, "RET") && str_equal_cstr(atomname_a, "C15")) ){
                order[i] = 2;
            } else {
                order[i] = 1;
            }
        }
    }
    
    return order;
}

static inline bool covalent_bond_heuristic(float dist_squared, md_element_t elem_a, md_element_t elem_b) {
    const float d = element_covalent_radii[elem_a] + element_covalent_radii[elem_b];
    const float d_min = d - 0.5f;
    const float d_max = d + 0.3f;
    return (d_min * d_min) < dist_squared && dist_squared < (d_max * d_max);
}

int64_t md_util_compute_covalent_bounds_upper_bound(const md_element_t* element, int64_t count) {
    int64_t result = 0;
    for (int64_t i = 0; i < count; ++i) {
        result += md_util_element_max_valence(element[i]);
    }
    return result;
}

md_array(md_bond_t) md_util_compute_covalent_bonds(const md_atom_data_t* atom, const md_unit_cell_t* cell, md_allocator_i* alloc) {
    ASSERT(atom);
    ASSERT(alloc);

    md_array(md_bond_t) bonds = 0;

    if (atom->count < 0) {
        MD_LOG_ERROR("Incorrect number of atoms: %i", (int)atom->count);
        return bonds;
    }
    
    if (!atom->x || !atom->y || !atom->z) {
        MD_LOG_ERROR("Missing atom field (x/y/z)");
        return bonds;
    }

    if (!atom->element) {
        MD_LOG_ERROR("Missing atom field element");
        return bonds;
    }

    if (cell->flags & MD_CELL_TRICLINIC) {
        MD_LOG_ERROR("Triclinic cells are not supported yet! Sorry!");
        return bonds;
    }
    
    const vec4_t pbc_ext = cell ? vec4_from_vec3(mat3_diag(cell->basis), 0) : vec4_zero();
    const int atom_count = (int)atom->count;
    
    if (atom->residue_idx) {
        // atom residue indices are given,
        // First find connections first within the residue, then to the next residue

        md_range_t range = {0, 0};
        md_range_t prev_range = {0, 0};
        for (int i = 0; i < atom_count; ++i) {
            while (i < atom_count && atom->residue_idx[i] == atom->residue_idx[range.beg]) ++i;
            range.end = i;

            for (int j = prev_range.beg; j < prev_range.end; ++j) {
                for (int k = range.beg; k < range.end; ++k) {
                    const vec4_t a = {atom->x[j], atom->y[j], atom->z[j], 0};
                    const vec4_t b = {atom->x[k], atom->y[k], atom->z[k], 0};
                    const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext);
                    if (covalent_bond_heuristic(d2, atom->element[j], atom->element[k])) {
                        md_array_push(bonds, ((md_bond_t){j, k}), alloc);
                    }
                }
            }

            for (int j = range.beg; j < range.end - 1; ++j) {
                for (int k = j + 1; k < range.end; ++k) {
                    const vec4_t a = {atom->x[j], atom->y[j], atom->z[j], 0};
                    const vec4_t b = {atom->x[k], atom->y[k], atom->z[k], 0};
                    const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext);
                    if (covalent_bond_heuristic(d2, atom->element[j], atom->element[k])) {
                        md_array_push(bonds, ((md_bond_t){j, k}), alloc);
                    }
                }
            }

            prev_range = range;
            range.beg  = range.end;
        }
    }
    else {
        if (atom_count < 100) {
            // If the system is small, just brute-force it
            for (int i = 0; i < (int)atom_count - 1; ++i) {
                for (int j = i + 1; j < (int)atom_count; ++j) {
                    const vec4_t a = {atom->x[i], atom->y[i], atom->z[i], 0};
                    const vec4_t b = {atom->x[j], atom->y[j], atom->z[j], 0};
                    const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext);
                    if (covalent_bond_heuristic(d2, atom->element[i], atom->element[j])) {
                        md_array_push(bonds, ((md_bond_t){i, j}), alloc);
                    }
                }
            }
        } else {
            // Resort to spatial acceleration structure
            md_spatial_hash_t* sh = md_spatial_hash_create_soa(atom->x, atom->y, atom->z, NULL, atom->count, cell, md_heap_allocator);

            const float cutoff = 3.0f;
            
            int indices[128];

            for (int i = 0; i < (int)atom_count; ++i) {
                const vec3_t pos = {atom->x[i], atom->y[i], atom->z[i]};

                const int64_t num_indices = md_spatial_hash_query_idx(indices, ARRAY_SIZE(indices), sh, pos, cutoff);
                for (int64_t iter = 0; iter < num_indices; ++iter) {
                    const int j = indices[iter];

                    if (atom->residue_idx && abs(atom->residue_idx[j] - atom->residue_idx[i]) > 2) {
                        continue;
                    }

                    const vec4_t pos_a = vec4_from_vec3(pos, 0);
                    const vec4_t pos_b = {atom->x[j], atom->y[j], atom->z[j], 0};
                    const float d2 = vec4_periodic_distance_squared(pos_a, pos_b, pbc_ext);
                    if (covalent_bond_heuristic(d2, atom->element[i], atom->element[j])) {
                        if (i < j) {
                            // Only store monotonic bonds connections
                            md_array_push(bonds, ((md_bond_t){i, j}), alloc);
                        }
                    }
                }
            }

            md_spatial_hash_free(sh);
        }
    }

    return bonds;
}

// Compute the prerequisite fields to enable hydrogen bond determination
// Heavily inspired by molstar
// https://github.com/molstar/molstar/blob/master/src/mol-model-props/computed/chemistry/valence-model.ts

/*
md_array(md_hbond_data_t) md_util_compute_hbond_data(const md_molecule_t* mol, md_index_data_t connectivity, md_allocator_i* alloc) {
    md_array(md_valence_t) atom_valence = md_array_create(md_valence_t, mol->atom.count, md_temp_allocator);
    memset(atom_valence, 0, md_array_bytes(atom_valence));

    md_array(int8_t) atom_hcount = md_array_create(int8_t, mol->atom.count, md_temp_allocator);
    memset(atom_hcount, 0, md_array_bytes(atom_hcount));

    md_array(md_order_t) order = md_util_compute_covalent_bond_order(mol->bonds, md_array_size(mol->bonds), mol->atom.name, mol->atom.residue_idx, mol->residue.name, md_temp_allocator);

    for (int64_t i = 0; i < md_array_size(mol->bonds); ++i) {
        const md_atom_idx_t *idx = mol->bonds[i].idx;
        atom_valence[idx[0]] += order[i];
        atom_valence[idx[1]] += order[i];

        if (mol->atom.element[idx[0]] == H) atom_hcount[idx[1]] += 1;
        if (mol->atom.element[idx[1]] == H) atom_hcount[idx[0]] += 1;
    }

    for (int64_t i = 0; i < mol->atom.count; ++i) {
        md_element_t elem = mol->atom.element[i];
        int charge  = 0;
        int implicit_hcount = 0;
        int degree  = (int)md_index_range_size(connectivity, i);
        int valence = atom_valence[i];
        int geom    = 0;

        bool multi_bond = (valence - degree) > 0;

        switch (elem)
        {
        case H:
            if (degree == 0) {
                charge = 1;
                geom = 0;
            } else if (degree == 1) {
                charge = 0;
                geom = 1;
            }
            break;
        case C:
            charge = 0;
            implicit_hcount = MAX(0, 4 - valence + MAX(0, -charge));
            geom = degree + implicit_hcount + MAX(0, -charge);
            break;
        case N:
            // @TODO: Complete this
            break;
        default:
            break;
        }
    }
}

md_array(md_bond_t) md_util_compute_hydrogen_bonds(const md_molecule_t* mol, md_allocator_i* alloc) {
    for (int64_t i = 0; i < mol->atom.count; ++i) {

    }
}
*/

bool md_util_compute_chain_data(md_chain_data_t* chain_data, const md_residue_idx_t res_idx[], int64_t atom_count, const md_bond_t bonds[], int64_t bond_count, md_allocator_i* alloc) {
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

    md_allocator_i* temp_alloc = md_temp_allocator;
    
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

    res_count += 1;

    if (res_count == 0) {
        MD_LOG_DEBUG("Dataset only contains single residue, no artificial chains can be created");
        return false;
    }
    
    md_array(bool) res_bond_to_next = md_array_create(bool, res_count, temp_alloc);
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
    const uint32_t new_cap = next_power_of_two32((uint32_t)new_capacity);
    md_array_grow(fifo->data, new_cap, fifo->alloc);
    fifo->cap = new_cap;
}

static void fifo_init(fifo_t* fifo, int64_t capacity, md_allocator_i* alloc) {
    ASSERT(fifo);
    ASSERT(0 <= capacity && capacity < UINT32_MAX);
    fifo->alloc = alloc;
    fifo->head = 0;
    fifo->tail = 0;
    fifo->cap = next_power_of_two32((uint32_t)capacity);
    fifo->data = md_array_create(int, fifo->cap, alloc);
#if DEBUG
    // Clear memory to make debugging easier
    MEMSET(fifo->data, 0, sizeof(int) * fifo->cap);
#endif
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
    fifo->head = 0;
    fifo->tail = 0;
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

static inline bool has_ring(md_index_data_t ring_data, const int* ring, int ring_size) {
    const int64_t num_rings = md_index_data_count(ring_data);
    for (int64_t i = 0; i < num_rings; ++i) {
        const int* ring_i     = md_index_range_beg(ring_data, i);
        const int ring_i_size = (int)md_index_range_size(ring_data, i);
        if (compare_ring(ring_i, ring_i_size, ring, ring_size)) return true;
    }
    return false;
}

#include <stdlib.h>

static int compare_int(void const* a, void const* b) {
    return ( *(const int*)a - *(const int*)b );
}

// Simplistic inplace bubble sort for small arrays
// In practice, this is only used to sort the small rings
static void sort_arr(int* arr, int n) {
    if (n < 16) {
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
    } else {
        qsort(arr, n, sizeof(int), compare_int);
    }
}

#define MIN_RING_SIZE 3
#define MAX_RING_SIZE 6

#define STB_DS_IMPLEMENTATION
#include <stb_ds.h>

md_index_data_t md_util_compute_rings(md_index_data_t atom_connectivity, md_allocator_i* alloc) {
    ASSERT(alloc);
    
    md_index_data_t ring_data = {0};
    const int64_t atom_count = md_index_data_count(atom_connectivity);    

    if (atom_count == 0) {
        return ring_data;
    }
    
    // We don't know how large datasets we operate on here, so its better to go for the big dog (vm arena) directly
    // Instead of risking not being able to fit into the temp ring buffer
    md_vm_arena_t arena = {0};
    md_vm_arena_init(&arena, GIGABYTES(1));
    md_allocator_i arena_alloc = md_vm_arena_create_interface(&arena);
    
    int* depth = md_vm_arena_push(&arena, atom_count * sizeof(int));
    MEMSET(depth, 0, atom_count * sizeof(int));
    
    int* pred = md_vm_arena_push(&arena, atom_count * sizeof(int));
    for (int64_t i = 0; i < atom_count; ++i) {
        pred[i] = -1;
    }

    // The capacity is arbitrary here, but will be resized if needed.
    fifo_t queue = fifo_create(64, &arena_alloc);

    typedef struct T {
        uint64_t key;
    } T;
    
    T* hm = NULL;
    size_t seed = 12;

    for (int i = 0; i < atom_count; ++i) {
        // Skip any atom which has alread been touched
        if (depth[i]) continue;

        fifo_clear(&queue);
        fifo_push(&queue, i);
        depth[i] = 1;
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);

            const int* eb = md_index_range_beg(atom_connectivity, idx);
            const int* ee = md_index_range_end(atom_connectivity, idx);
            ASSERT(eb);
            ASSERT(ee);
            ASSERT(eb <= ee);
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
                        if (l > -1 && depth[l] >= d) {
                            ring[n++] = l;
                            l = pred[l];
                        }

                        if (r > -1 && depth[r] >= d) {
                            ring[n++] = r;
                            r = pred[r];
                        }

                        if (l == r) {
                            ring[n++] = l;
                            if (l == -1)
                                n = 0;
                            break;
                        }
                    }

                    if (MIN_RING_SIZE <= n && n <= MAX_RING_SIZE) {
                        sort_arr(ring, n);
                        size_t key = stbds_hash_bytes(ring, n * sizeof(int), seed);
                        if (hmgeti(hm, key) == -1) {
                            hmputs(hm, (T){key});
                            md_index_data_push(&ring_data, ring, n, alloc);
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
    hmfree(hm);
    md_vm_arena_free(&arena);

    return ring_data;
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE

static inline bool test_bit(uint64_t* bits, int idx) {
	return (bits[idx >> 6] & (1ULL << (idx & 63)));
}

static inline void set_bit(uint64_t* bits, int idx) {
    bits[idx >> 6] |= (1ULL << (idx & 63));
}

md_index_data_t md_util_compute_structures(md_index_data_t atom_connectivity, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    md_index_data_t structures = {0};
    const int64_t num_atoms = md_index_data_count(atom_connectivity);

    md_vm_arena_t arena = { 0 };
    md_vm_arena_init(&arena, GIGABYTES(1));
    md_allocator_i arena_alloc = md_vm_arena_create_interface(&arena);

    // Create a bitfield to keep track of which atoms have been visited
    const uint64_t bytes = DIV_UP(num_atoms, 64) * sizeof(uint64_t);
    uint64_t* visited = md_vm_arena_push(&arena, bytes);

    // Clear the visited bitfield, this might be redundant due to virtual memory
    MEMSET(visited, 0, bytes);

    // The capacity is arbitrary here, and will be resized if needed.
    fifo_t queue = fifo_create(128, &arena_alloc);

    md_array(int) indices = 0;
    md_array_ensure(indices, 256, &arena_alloc);

    for (int i = 0; i < num_atoms; ++i) {
        // Skip any atom which has already been touched
        if (test_bit(visited, i)) continue;

        fifo_clear(&queue);
        fifo_push(&queue, i);
        while (!fifo_empty(&queue)) {
            int cur = fifo_pop(&queue);
            if (test_bit(visited, cur)) continue;
            set_bit(visited, cur);

            md_array_push(indices, cur, &arena_alloc);
            
            const int* eb = md_index_range_beg(atom_connectivity, cur);
            const int* ee = md_index_range_end(atom_connectivity, cur);
            for (const int* it = eb; it != ee; ++it) {
                int next = *it;
                if (!test_bit(visited, next)) {
                    fifo_push(&queue, next);
                }
            }
        }

        // Sort the indices within the structure for more coherent memory access
        sort_arr(indices, (int)md_array_size(indices));
        
        // Here we should have exhausted every atom that is connected to index i.
        md_index_data_push(&structures, indices, md_array_size(indices), alloc);
        md_array_shrink(indices, 0);
    }
    
    md_vm_arena_free(&arena);

    return structures;
}

void md_util_grow_mask_by_bonds(md_bitfield_t* mask, const struct md_molecule_t* mol, int extent, const md_bitfield_t* viable_mask) {
    ASSERT(mask);
    ASSERT(mol);

    if (extent <= 0) return;
    if (extent >= 255) {
        MD_LOG_DEBUG("Maximum supported growth extent is 255, the extent will be clamped to 255");
        extent = 255;
    }
    // The initial depth is 1, so we need to add one.
    extent += 1;
    
    const int64_t mask_size = md_bitfield_popcount(mask);
    if (!mask_size) return;

    md_vm_arena_t arena = { 0 };
    md_vm_arena_init(&arena, GIGABYTES(1));
    md_allocator_i arena_alloc = md_vm_arena_create_interface(&arena);

    int* indices = md_vm_arena_push(&arena, mask_size * sizeof(int));
    const int64_t num_indices = md_bitfield_extract_indices(indices, mask_size, mask);
    ASSERT(num_indices == mask_size);

    fifo_t queue = fifo_create(64, &arena_alloc);

    int64_t num_atoms = mol->atom.count;
    uint8_t* depth = md_vm_arena_push(&arena, num_atoms * sizeof(uint8_t));
    MEMSET(depth, 0, num_atoms * sizeof(uint8_t));

    {
        md_bitfield_iter_t it = md_bitfield_iter_create(mask);
        while (md_bitfield_iter_next(&it)) {
            int idx = (int)md_bitfield_iter_idx(&it);
            depth[idx] = 1;
        }
    }
    
    for (int j = 0; j < num_indices; ++j) {
        int i = indices[j];

        fifo_clear(&queue);
        fifo_push(&queue, i);
        
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);
            md_bitfield_set_bit(mask, idx);

            const int* eb = md_index_range_beg(mol->connectivity, idx);
            const int* ee = md_index_range_end(mol->connectivity, idx);
            for (const int* it = eb; it != ee; ++it) {
                int next = *it;
                if (viable_mask && md_bitfield_test_bit(viable_mask, next) == false) {
                    continue;
                }
                // Only explore the path if it is not already marked and it is within the extent.
                // Or if the current path is 'cheaper'
                const uint8_t dc = depth[idx] + 1;
                const uint8_t dn = depth[next];
                if (dn == 0 || dc < dn) {
                    depth[next] = dc;
                    if (dc <= extent) {
                        fifo_push(&queue, next);
                    }
                }
            }
        }
    }

    md_vm_arena_free(&arena);
}

void md_util_grow_mask_by_radius(md_bitfield_t* mask, const struct md_molecule_t* mol, float radius, const md_bitfield_t* viable_mask) {
    ASSERT(mask);
    ASSERT(mol);
    
    if (radius <= 0.0f) return;
    

    md_vm_arena_t arena = { 0 };
    md_vm_arena_init(&arena, GIGABYTES(1));
    md_allocator_i arena_alloc = md_vm_arena_create_interface(&arena);
    
    int32_t* indices = 0;
    int64_t count = mol->atom.count;
        
    if (viable_mask) {
        md_bitfield_t tmp_bf = md_bitfield_create(&arena_alloc);
        md_bitfield_andnot(&tmp_bf, viable_mask, mask);
            
        const int64_t num_atoms = md_bitfield_popcount(&tmp_bf);
        indices = md_vm_arena_push(&arena, num_atoms * sizeof(int32_t));
        count = md_bitfield_extract_indices(indices, num_atoms, &tmp_bf);
    }

    md_spatial_hash_t* ctx = md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, indices, count, &mol->unit_cell, &arena_alloc);
    
    md_bitfield_t old_mask = md_bitfield_create(&arena_alloc);
    md_bitfield_copy(&old_mask, mask);

    md_bitfield_iter_t it = md_bitfield_iter_create(&old_mask);
    while (md_bitfield_iter_next(&it)) {
        int idx = (int)md_bitfield_iter_idx(&it);
        const vec3_t pos = {mol->atom.x[idx], mol->atom.y[idx], mol->atom.z[idx]};
        md_spatial_hash_query_bits(mask, ctx, pos, radius);
    }

    md_vm_arena_free(&arena);
}

static void compute_aabb(vec3_t* out_aabb_min, vec3_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, const float* in_r, int64_t count, uint64_t xyz_stride) {
    md_simd_f32_t vx_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vy_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vz_min = md_simd_set1_f32(+FLT_MAX);

    md_simd_f32_t vx_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vy_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vz_max = md_simd_set1_f32(-FLT_MAX);

    int64_t i = 0;
    const int64_t simd_count = (count / md_simd_width_f32) * md_simd_width_f32;

    if (in_r) {
        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t x = md_simd_load_f32((const float*)((const char*)in_x + i * xyz_stride));
            md_simd_f32_t y = md_simd_load_f32((const float*)((const char*)in_y + i * xyz_stride));
            md_simd_f32_t z = md_simd_load_f32((const float*)((const char*)in_z + i * xyz_stride));
            md_simd_f32_t r = md_simd_load_f32(in_r + i);            

            vx_min = md_simd_min(vx_min, md_simd_sub(x, r));
            vy_min = md_simd_min(vy_min, md_simd_sub(y, r));
            vz_min = md_simd_min(vz_min, md_simd_sub(z, r));

            vx_max = md_simd_max(vx_max, md_simd_add(x, r));
            vy_max = md_simd_max(vy_max, md_simd_add(y, r));
            vz_max = md_simd_max(vz_max, md_simd_add(z, r));
        }
    } else {
        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t x = md_simd_load_f32((const float*)((const char*)in_x + i * xyz_stride));
            md_simd_f32_t y = md_simd_load_f32((const float*)((const char*)in_y + i * xyz_stride));
            md_simd_f32_t z = md_simd_load_f32((const float*)((const char*)in_z + i * xyz_stride));

            vx_min = md_simd_min(vx_min, x);
            vy_min = md_simd_min(vy_min, y);
            vz_min = md_simd_min(vz_min, z);

            vx_max = md_simd_max(vx_max, x);
            vy_max = md_simd_max(vy_max, y);
            vz_max = md_simd_max(vz_max, z);
        }
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
    if (in_r) {
        for (; i < count; ++i) {
            for (; i < count; ++i) {
                vec4_t c = {
                    *(const float*)((const char*)in_x + i * xyz_stride),
                    *(const float*)((const char*)in_y + i * xyz_stride),
                    *(const float*)((const char*)in_z + i * xyz_stride),
                    0
                };
                vec4_t r = vec4_set1(in_r[i]);
                
                aabb_min = vec4_min(aabb_min, vec4_sub(c, r));
                aabb_max = vec4_max(aabb_max, vec4_add(c, r));
            }
        }
    } else {
        for (; i < count; ++i) {
            vec4_t c = {
                *(const float*)((const char*)in_x + i * xyz_stride),
                *(const float*)((const char*)in_y + i * xyz_stride),
                *(const float*)((const char*)in_z + i * xyz_stride),
                0
            };
            
            aabb_min = vec4_min(aabb_min, c);
            aabb_max = vec4_max(aabb_max, c);
        }
    }

    *out_aabb_min = vec3_from_vec4(aabb_min);
    *out_aabb_max = vec3_from_vec4(aabb_max);
}

static void compute_aabb_indexed(vec3_t* out_aabb_min, vec3_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int* indices, int64_t count, uint64_t xyz_stride) {   
    if (count == 0) return;
    
    vec4_t aabb_min = vec4_set1( FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);

    // @PERF(Robin), we are only using 4 lanes here since we are using idirect loads from supplied indices.
    // This could boil down to a gather instruction. However, the gather instruction
    // does not allow arbitrary strides, but perhaps there could be a work around where the indices are scaled directly
    // HOWEVER this only applies for AVX2 and above, so its not a big deal for now.

    if (in_r) {
        for (int64_t i = 0; i < count; ++i) {
            int idx = indices[i];
            vec4_t c = {
                *(const float*)((const char*)in_x + idx * xyz_stride),
                *(const float*)((const char*)in_y + idx * xyz_stride),
                *(const float*)((const char*)in_z + idx * xyz_stride),
                0
            };
            vec4_t r = vec4_set1(in_r[idx]);

            aabb_min = vec4_min(aabb_min, vec4_sub(c, r));
            aabb_max = vec4_max(aabb_max, vec4_add(c, r));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            int idx = indices[i];
            vec4_t c = {
                *(const float*)((const char*)in_x + idx * xyz_stride),
                *(const float*)((const char*)in_y + idx * xyz_stride),
                *(const float*)((const char*)in_z + idx * xyz_stride),
                0
            };

            aabb_min = vec4_min(aabb_min, c);
            aabb_max = vec4_max(aabb_max, c);
        }
    }

    *out_aabb_min = vec3_from_vec4(aabb_min);
    *out_aabb_max = vec3_from_vec4(aabb_max);
}

void md_util_compute_aabb(vec3_t* aabb_min, vec3_t* aabb_max, const vec3_t* xyz, const float* r, int64_t count) {
    compute_aabb(aabb_min, aabb_max, &xyz->x, &xyz->y, &xyz->z, r, count, sizeof(vec3_t));
}
void md_util_compute_aabb_soa(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, const float* r, int64_t count) {
    compute_aabb(aabb_min, aabb_max, x, y, z, r, count, sizeof(float));
}

// Computes the minimum axis aligned bounding box for a set of points with a given radius (radius is optional), indices are used to select a subset of points
void md_util_compute_aabb_indexed(vec3_t* aabb_min, vec3_t* aabb_max, const vec3_t* xyz, const float* r, const int32_t* indices, int64_t index_count) {
    compute_aabb_indexed(aabb_min, aabb_max, &xyz->x, &xyz->y, &xyz->z, r, indices, index_count, sizeof(vec3_t));
}

void md_util_compute_aabb_indexed_soa(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, const float* r, const int32_t* indices, int64_t index_count) {
    compute_aabb_indexed(aabb_min, aabb_max, x, y, z, r, indices, index_count, sizeof(float));
}

static vec3_t compute_com_periodic_trig_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, int64_t count, vec3_t ext_max) {
    const vec4_t scl = vec4_div(vec4_set1(TWO_PI), vec4_from_vec3(ext_max, 1.0f));
    vec4_t acc_c = {0};
    vec4_t acc_s = {0};
    vec4_t acc_xyzw = {0};

    for (int64_t i = 0; i < count; ++i) {
        int64_t idx  = in_idx ? in_idx[i] : i;
        vec4_t xyzw  = in_xyzw[idx];
        vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(1,0,0,0));
        vec4_t theta = vec4_mul(xyzw, scl);
        vec4_t s,c;
        vec4_sincos(theta, &s, &c);
        acc_s = vec4_add(acc_s, vec4_mul(s, www1));
        acc_c = vec4_add(acc_c, vec4_mul(c, www1));
        acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
    }

    vec3_t com;
    const float w = acc_xyzw.w;
    for (int i = 0; i < 3; ++i) {
        if (ext_max.elem[i] > 0) {
            const double y = acc_s.elem[i] / w;
            const double x = acc_c.elem[i] / w;
            const double r2 = x*x + y*y;
            double theta_prim = PI;
            if (r2 > 1.0e-15) {
                theta_prim += atan2(-y, -x);
            }
            com.elem[i] = (float)((theta_prim / TWO_PI) * ext_max.elem[i]);
        } else {
            com.elem[i] = acc_xyzw.elem[i] / w;
        }
    }

    return com;
}

// Regular version, deperiodization is done with respect to the previous element
static vec3_t compute_com_periodic_reg_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, int64_t count, vec3_t ext) {
    const vec4_t period = vec4_from_vec3(ext, 0.0f);

    int32_t idx0     = in_idx ? in_idx[0] : 0;
    vec4_t acc_xyzw  = in_xyzw[idx0];
    vec4_t prev_xyzw = in_xyzw[idx0];

    for (int64_t i = 1; i < count; ++i) {
        int64_t idx = in_idx ? in_idx[i] : i;
        vec4_t xyzw = vec4_deperiodize(in_xyzw[idx], prev_xyzw, period);
        vec4_t www1 = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(1,0,0,0));
        acc_xyzw    = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
        prev_xyzw   = xyzw;
    }

    vec3_t com = vec3_from_vec4(vec4_div(acc_xyzw, vec4_splat_w(acc_xyzw)));
    return com;
}

static float compute_com_periodic_trig(const float* in_x, const float* in_w, const int32_t* in_idx, int64_t count, float x_max) {
    double acc_c = 0;
    double acc_s = 0;
    double acc_w = 0;

    const double scl = TWO_PI / x_max;
    int64_t i = 0;

#if 0
    const float fscl = (float)(TWO_PI / x_max);
    md_simd_f32_t v_acc_c = md_simd_zero_f32();
    md_simd_f32_t v_acc_s = md_simd_zero_f32();
    md_simd_f32_t v_acc_w = md_simd_zero_f32();

    const int64_t simd_count = ALIGN_TO(count, md_simd_width_f32) - md_simd_width_f32;
    if (in_idx) {
        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t v_x = md_simd_gather_f32(in_x, in_idx + i);
            md_simd_f32_t v_w = in_w ? md_simd_gather_f32(in_w, in_idx + i) : md_simd_set1_f32(1.0f);
            md_simd_f32_t v_theta = md_simd_mul(v_x, md_simd_set1_f32(fscl));
            md_simd_f32_t v_c, v_s;
            md_simd_sincos(v_theta, &v_s, &v_c);
            v_acc_s = md_simd_add(v_acc_s, md_simd_mul(v_s, v_w));
            v_acc_c = md_simd_add(v_acc_c, md_simd_mul(v_c, v_w));
            v_acc_w = md_simd_add(v_acc_w, v_w);
        }
    } else {
        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t v_x = md_simd_load_f32(in_x + i);
            md_simd_f32_t v_w = in_w ? md_simd_load_f32(in_w + i) : md_simd_set1_f32(1.0f);
            md_simd_f32_t v_theta = md_simd_mul(v_x, md_simd_set1_f32(fscl));
            md_simd_f32_t v_c, v_s;
            md_simd_sincos(v_theta, &v_s, &v_c);
            v_acc_s = md_simd_add(v_acc_s, md_simd_mul(v_s, v_w));
            v_acc_c = md_simd_add(v_acc_c, md_simd_mul(v_c, v_w));
            v_acc_w = md_simd_add(v_acc_w, v_w);
        }
    }
    acc_s = md_simd_hsum(v_acc_s);
    acc_c = md_simd_hsum(v_acc_c);
    acc_w = md_simd_hsum(v_acc_w);
#endif

    for (; i < count; ++i) {
        int64_t idx = in_idx ? in_idx[i] : i;
        double theta = in_x[idx] * scl;
        double w = in_w ? in_w[idx] : 1.0;
        acc_c += w * cos(theta);
        acc_s += w * sin(theta);
        acc_w += w;
    }

    acc_w = (acc_w == 0) ? count : acc_w;
    const double y = acc_s / acc_w;
    const double x = acc_c / acc_w;
    const double r2 = x*x + y*y;

    double theta_prim = PI;
    if (r2 > 1.0e-15) {
        theta_prim += atan2(-y, -x);
    }

    return (float)((theta_prim / TWO_PI) * x_max);
}

// Regular version, deperiodization is done with respect to the previous element
static float compute_com_periodic_reg(const float* in_x, const float* in_w, const int32_t* in_idx, int64_t count, float x_max) {
    if (count <= 0) {
        return 0;
    }

    int32_t idx0 = in_idx ? in_idx[0] : 0;
    double acc_x = in_x[idx0];
    double acc_w = in_w ? in_w[idx0] : 1.0;
    double prev_x = in_x[idx0];
    double d_x_max = x_max;

    for (int64_t i = 1; i < count; ++i) {
        int64_t idx = in_idx ? in_idx[i] : i;
        double x = in_x[idx];
        double dx = deperiodize(x, prev_x, d_x_max);
        double w = in_w ? in_w[idx] : 1.0;
        acc_x += dx * w;
        acc_w += w;
        prev_x = dx;
    }

    return (float)(acc_x / acc_w);
}

static float compute_com(const float* in_x, const float* in_w, const int32_t* in_idx, int64_t count) {
    ASSERT(in_x);

    if (count == 0)
        return 0.0f;

    float acc_x = 0;
    float acc_w = 0;
    int64_t i = 0;

    md_simd_f32_t v_acc_x = md_simd_zero_f32();
    md_simd_f32_t v_acc_w = md_simd_zero_f32();
    const int64_t simd_count = ALIGN_TO(count, md_simd_width_f32) - md_simd_width_f32;
    if (in_idx) {
        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t v_x = md_simd_gather_f32(in_x, in_idx + i);
            md_simd_f32_t v_w = in_w ? md_simd_gather_f32(in_w, in_idx) : md_simd_set1_f32(1.0f);
            v_acc_x = md_simd_add(v_acc_x, md_simd_mul(v_x, v_w));
            v_acc_w = md_simd_add(v_acc_w, v_w);
        }
    } else {
        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t v_x = md_simd_load_f32(in_x + i);
            md_simd_f32_t v_w = in_w ? md_simd_load_f32(in_w) : md_simd_set1_f32(1.0f);
            v_acc_x = md_simd_add(v_acc_x, md_simd_mul(v_x, v_w));
            v_acc_w = md_simd_add(v_acc_w, v_w);
        }
    }
    acc_x = md_simd_hsum(v_acc_x);
    acc_w = md_simd_hsum(v_acc_w);

    for (; i < count; ++i) {
        int64_t idx = in_idx ? in_idx[i] : i;
        float x = in_x[idx];
        float w = in_w ? in_w[idx] : 1.0f;
        acc_x += x * w;
        acc_w += w;
    }

    float com = acc_x / acc_w;
    return com;
}

vec3_t md_util_compute_com_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, int64_t count) {
    ASSERT(in_xyzw);

    if (count <= 0) {
        return (vec3_t) {0,0,0};
    }

    // Use vec4 here so we can utilize SSE vectorization if applicable
    // @TODO: Vectorize with full register width
    vec4_t sum_xyzw = {0,0,0,0};
    for (int64_t i = 0; i < count; ++i) {
        int64_t idx = in_idx ? in_idx[i] : i;
        vec4_t xyzw = in_xyzw[idx];
        vec4_t www1 = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(1,0,0,0));
        sum_xyzw = vec4_add(sum_xyzw, vec4_mul(xyzw, www1));
    }

    return vec3_div_f(vec3_from_vec4(sum_xyzw), sum_xyzw.w);
}

vec3_t md_util_compute_com(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* indices, int64_t count) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    if (count <= 0) {
        return (vec3_t) {0,0,0};
    }

    // Use vec4 here so we can utilize SSE vectorization if applicable
    // @TODO: Vectorize with full register width
    vec4_t sum_xyzw = {0,0,0,0};
    if (in_w) {
        for (int64_t i = 0; i < count; ++i) {
            int64_t idx = indices ? indices[i] : i;
            const vec4_t xyz1 = vec4_set(in_x[idx], in_y[idx], in_z[idx], 1.0f);
            sum_xyzw = vec4_add(sum_xyzw, vec4_mul_f(xyz1, in_w[idx]));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            int64_t idx = indices ? indices[i] : i;
            const vec4_t xyz0 = vec4_set(in_x[idx], in_y[idx], in_z[idx], 0);
            sum_xyzw = vec4_add(sum_xyzw, xyz0);
        }
        sum_xyzw.w = (float)count;
    }

    return vec3_div_f(vec3_from_vec4(sum_xyzw), sum_xyzw.w);
}

#ifndef MD_UTIL_COMPUTE_COM_USE_TRIG
#define MD_UTIL_COMPUTE_COM_USE_TRIG 1
#endif

// Love the elegance of using trigonometric functions, unsure of the performance...
// @TODO: sin, cos and atan2 can and should of course be vectorized.
vec3_t md_util_compute_com_vec4_ortho(const vec4_t* in_xyzw, const int32_t* indices, int64_t count, vec3_t pbc_ext) {
    ASSERT(in_xyzw);

    if (count <= 0) {
        return (vec3_t) {0,0,0};
    }

#if MD_UTIL_COMPUTE_COM_USE_TRIG
    return compute_com_periodic_trig_vec4(in_xyzw, indices, count, pbc_ext);
#else
    return compute_com_periodic_reg_vec4(in_xyzw, indices, count, pbc_ext);
#endif
}

vec3_t md_util_compute_com_ortho(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* indices, int64_t count, vec3_t pbc_ext) {
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

#if MD_UTIL_COMPUTE_COM_USE_TRIG
    // We need to pick the version based on each component of pdc_ext, since one or more may be zero
    float x = pbc_ext.x > 0 ? compute_com_periodic_trig(in_x, in_w, indices, count, pbc_ext.x) : compute_com(in_x, in_w, indices, count);
    float y = pbc_ext.y > 0 ? compute_com_periodic_trig(in_y, in_w, indices, count, pbc_ext.y) : compute_com(in_y, in_w, indices, count);
    float z = pbc_ext.z > 0 ? compute_com_periodic_trig(in_z, in_w, indices, count, pbc_ext.z) : compute_com(in_z, in_w, indices, count);
#else
    // We need to pick the version based on each component of pdc_ext, since one or more may be zero
    float x = pbc_ext.x > 0 ? compute_com_periodic_reg(in_x, in_w, indices, count, pbc_ext.x) : compute_com(in_x, in_w, indices, count);
    float y = pbc_ext.y > 0 ? compute_com_periodic_reg(in_y, in_w, indices, count, pbc_ext.y) : compute_com(in_y, in_w, indices, count);
    float z = pbc_ext.z > 0 ? compute_com_periodic_reg(in_z, in_w, indices, count, pbc_ext.z) : compute_com(in_z, in_w, indices, count);
#endif

    return (vec3_t) {x, y, z};
}

/*
vec3_t md_util_compute_com_ortho(const float *in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* indices, int64_t count, vec3_t box) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    if (count <= 0) {
        return (vec3_t){0, 0, 0};
    }
    const int idx0 = indices ? indices[0] : 0;

    if (count == 0) {
        return (vec3_t){in_x[idx0], in_y[idx0], in_z[idx0]};
    }
    
    const vec4_t ext = vec4_from_vec3(box, 0);

    vec4_t ref = vec4_set(in_x[idx0], in_y[idx0], in_z[idx0], 1.0f);
    vec4_t acc = ref;

    if (in_w) {
        acc = vec4_mul_f(acc, in_w[idx0]);
        for (int64_t i = 1; i < count; ++i) {
            const int idx = indices ? indices[i] : i;
            const vec4_t xyzw = vec4_deperiodize(vec4_set(in_x[idx], in_y[idx], in_z[idx], 1.0f), ref, ext);
            acc = vec4_add(acc, vec4_mul_f(xyzw, in_w[idx]));
            ref = xyzw;
        }
    } else {
        for (int64_t i = 1; i < count; ++i) {
            const int idx = indices ? indices[i] : i;
            const vec4_t xyzw = vec4_deperiodize(vec4_set(in_x[idx], in_y[idx], in_z[idx], 1.0f), ref, ext);
            acc = vec4_add(acc, xyzw);
            ref = xyzw;
        }
    }

    return vec3_from_vec4(vec4_div_f(acc, acc.w));
}

vec3_t md_util_compute_com_vec3_ortho(const vec3_t* in_xyz, const float* in_w, int64_t count, vec3_t pbc_ext) {
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

#if MD_COMPILER_MSVC
#   pragma warning( push )
#   pragma warning( disable : 4244 )
#endif

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

md_unit_cell_t md_util_unit_cell_from_extent(double x, double y, double z) {
    if (x == 0.0 && y == 0.0 && z == 0.0) {
        return (md_unit_cell_t) {0};
    }

    md_unit_cell_t cell = {
        .basis = {
            (float)x, 0, 0,
            0, (float)y, 0,
            0, 0, (float)z,
        },
        .inv_basis = {
            x != 0.0 ? 1.0/x : 0, 0, 0,
            0, y != 0.0 ? 1.0/y : 0, 0,
            0, 0, z != 0.0 ? 1.0/z : 0,
        },
        .flags = MD_CELL_ORTHOGONAL | (x != 0.0 ? MD_CELL_X : 0) | (y != 0.0 ? MD_CELL_Y : 0) | (z != 0.0 ? MD_CELL_Z : 0),
    };
    return cell;
}

md_unit_cell_t md_util_unit_cell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma) {
    if (a == 0 && b == 0 && c == 0) {
        return (md_unit_cell_t) {0};
    }

    if (alpha == 90.0 && beta == 90.0 && gamma == 90.0) {
        return md_util_unit_cell_from_extent(a,b,c);
    }

    alpha = DEG_TO_RAD(alpha);
    beta  = DEG_TO_RAD(beta);
    gamma = DEG_TO_RAD(gamma);

    const double cosa = cos(alpha);
    const double cosb = cos(beta);
    const double cosg = cos(gamma);
    const double sing = sin(gamma);
    const double x = (cosa - cosb * cosg) / sing;

    // Values for basis. Keep in double precision for inverse
    const double M[9] = {
               a,        0,                             0,
        b * cosg, b * sing,                             0,
        c * cosb,    c * x, c * sqrt(1 - cosb*cosb - x*x),
    };

    // Inverse of M
    const double I[9] = {
                                        1.0/M[0],                 0,        0,
                               -M[3]/(M[0]*M[4]),          1.0/M[4],        0,
        (M[3]*M[7]/(M[0]*M[4]) - M[6]/M[0])/M[8], -M[7]/(M[4]*M[8]), 1.0/M[8],
    };
    
    md_unit_cell_t cell = {
        .basis = {
            M[0], M[1], M[2], M[3], M[4], M[5], M[6], M[7], M[8],
        },
        .inv_basis = {
            I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8],
        },
        .flags = MD_CELL_TRICLINIC | (a != 0.0 ? MD_CELL_X : 0) | (b != 0.0 ? MD_CELL_Y : 0) | (c != 0.0 ? MD_CELL_Z : 0),
    };

    return cell;
}

md_unit_cell_t md_util_unit_cell_from_matrix(mat3_t M) {
    if (M.elem[0][1] == 0 && M.elem[0][2] == 0 && M.elem[1][0] == 0 && M.elem[1][2] == 0 && M.elem[2][0] == 0 && M.elem[2][1] == 0) {
        return md_util_unit_cell_from_extent(M.elem[0][0], M.elem[1][1], M.elem[2][2]);
    } else {
        const float* N = &M.elem[0][0];
        
        // Inverse of M
        const double I[9] = {
            1.0/(double)N[0], 0, 0,
            -(double)N[3]/((double)N[0]*(double)N[4]), 1.0/(double)N[4], 0,
            ((double)N[3]*(double)N[7]/((double)N[0]*(double)N[4]) - (double)N[6]/(double)N[0])/(double)N[8], -(double)N[7]/((double)N[4]*(double)N[8]), 1.0/(double)N[8],
        };
        
        const float a = vec3_length_squared(M.col[0]);
        const float b = vec3_length_squared(M.col[1]);
        const float c = vec3_length_squared(M.col[2]);
        
        md_unit_cell_t cell = {
            .basis = M,
            .inv_basis = {
                I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8],
            },
            .flags = MD_CELL_TRICLINIC | (a != 0 ? MD_CELL_X : 0) | (b != 0 ? MD_CELL_Y : 0) | (c != 0 ? MD_CELL_Z : 0),
        };
        return cell;
    }
}

// Blatantly stolen from MDAnalysis project
// https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/include/calc_distances.h
static void minimum_image_triclinic(float dx[3], const float box[3][3]) {
    /*
    * Minimum image convention for triclinic systems, modelled after domain.cpp
    * in LAMMPS.
    * Assumes that there is a maximum separation of 1 box length (enforced in
    * dist functions by moving all particles to inside the box before
    * calculating separations).
    * Assumes box having zero values for box[1], box[2] and box[5]:
    *   /  a_x   0    0   \                 /  0    1    2  \
    *   |  b_x  b_y   0   |       indices:  |  3    4    5  |
    *   \  c_x  c_y  c_z  /                 \  6    7    8  /
    */
    double dx_min[3] = {0.0, 0.0, 0.0};
    double dsq_min = FLT_MAX;
    double dsq;
    double rx;
    double ry[2];
    double rz[3];
    int ix, iy, iz;
    for (ix = -1; ix < 2; ++ix) {
        rx = dx[0] + box[0][0] * ix;
        for (iy = -1; iy < 2; ++iy) {
            ry[0] = rx + box[1][0] * iy;
            ry[1] = dx[1] + box[1][1] * iy;
            for (iz = -1; iz < 2; ++iz) {
                rz[0] = ry[0] + box[2][0] * iz;
                rz[1] = ry[1] + box[2][1] * iz;
                rz[2] = dx[2] + box[2][2] * iz;
                dsq = rz[0] * rz[0] + rz[1] * rz[1] + rz[2] * rz[2];
                if (dsq < dsq_min) {
                    dsq_min = dsq;
                    dx_min[0] = rz[0];
                    dx_min[1] = rz[1];
                    dx_min[2] = rz[2];
                }
            }
        }
    }
    dx[0] = (float)dx_min[0];
    dx[1] = (float)dx_min[1];
    dx[2] = (float)dx_min[2];
}

void md_util_unit_cell_distance_array(float* out_dist, const vec3_t* coord_a, int64_t num_a, const vec3_t* coord_b, int64_t num_b, const md_unit_cell_t* cell) {
    if (cell->flags == 0) {
        for (int64_t i = 0; i < num_a; ++i) {
        	for (int64_t j = 0; j < num_b; ++j) {
                out_dist[i * num_b + j] = vec3_distance(coord_a[i], coord_b[j]);
            }
        }
    }
    else if (cell->flags & MD_CELL_ORTHOGONAL) {
        const vec4_t box = {cell->basis.elem[0][0], cell->basis.elem[1][1], cell->basis.elem[2][2], 0};
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                vec4_t a = vec4_from_vec3(coord_a[i], 0);
                vec4_t b = vec4_from_vec3(coord_b[j], 0);
                out_dist[i * num_b + j] = vec4_periodic_distance(a, b, box);
            }
        }
    } else if (cell->flags & MD_CELL_TRICLINIC) {
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, cell->basis.elem);
                out_dist[i * num_b + j] = vec3_length(dx);
            }
        }
    }
}

float md_util_unit_cell_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, int64_t num_a, const vec3_t* coord_b, int64_t num_b, const md_unit_cell_t* cell) {
    int64_t min_i = 0;
    int64_t min_j = 0;
    float min_dist = FLT_MAX;

    if (cell->flags == 0) {
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                const float d = vec3_distance(coord_a[i], coord_b[j]);
                if (d < min_dist) {
                    min_dist = d;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }
    else if (cell->flags & MD_CELL_ORTHOGONAL) {
        const vec4_t box = {cell->basis.elem[0][0], cell->basis.elem[1][1], cell->basis.elem[2][2], 0};
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                vec4_t a = vec4_from_vec3(coord_a[i], 0);
                vec4_t b = vec4_from_vec3(coord_b[j], 0);
                const float d = vec4_periodic_distance(a, b, box);
                if (d < min_dist) {
                    min_dist = d;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    } else if (cell->flags & MD_CELL_TRICLINIC) {
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, cell->basis.elem);
                const float d = vec3_length(dx);
                if (d < min_dist) {
                    min_dist = d;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }
    if (min_dist < FLT_MAX) {
        *out_idx_a = min_i;
        *out_idx_b = min_j;
    }
    return min_dist;
}

float md_util_unit_cell_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, int64_t num_a, const vec3_t* coord_b, int64_t num_b, const md_unit_cell_t* cell) {
    int64_t max_i = 0;
    int64_t max_j = 0;
    float max_dist = 0;

    if (cell->flags == 0) {
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                const float d = vec3_distance(coord_a[i], coord_b[j]);
                if (d > max_dist) {
                    max_dist = d;
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }
    else if (cell->flags & MD_CELL_ORTHOGONAL) {
        const vec4_t box = {cell->basis.elem[0][0], cell->basis.elem[1][1], cell->basis.elem[2][2], 0};
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                vec4_t a = vec4_from_vec3(coord_a[i], 0);
                vec4_t b = vec4_from_vec3(coord_b[j], 0);
                const float d = vec4_periodic_distance(a, b, box);
                if (d > max_dist) {
                    max_dist = d;
                    max_i = i;
                    max_j = j;
                }
            }
        }
    } else if (cell->flags & MD_CELL_TRICLINIC) {
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < num_a; ++i) {
            for (int64_t j = 0; j < num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, cell->basis.elem);
                const float d = vec3_length(dx);
                if (d > max_dist) {
                    max_dist = d;
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }

    if (max_dist > 0) {
        if (out_idx_a) *out_idx_a = max_i;
        if (out_idx_b) *out_idx_b = max_j;
    }

    return max_dist;
}

#if MD_COMPILER_MSVC
#   pragma warning( pop )
#endif

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

bool md_util_unwrap_ortho(float* x, float* y, float* z, md_index_data_t structures, vec3_t box_ext) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
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

bool md_util_deperiodize_system(float* x, float* y, float* z, const md_unit_cell_t* cell, const md_molecule_t* mol) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(cell);

    const float* w = mol->atom.mass;
    const int64_t num_atoms = mol->atom.count;
    const int64_t num_structures = md_index_data_count(mol->structures);
    const md_index_data_t structures = mol->structures;

    if (cell->flags & MD_CELL_TRICLINIC) {
        MD_LOG_ERROR("Triclinic cells are not supported!");
        return false;
    }

    const vec3_t box = mat3_diag(cell->basis);
    
    md_util_pbc_ortho(x, y, z, num_atoms, box);
    if (num_structures > 0) {
        md_util_unwrap_ortho(x, y, z, structures, box);
    
        // Place any defined covalent structure such that its center of mass is within the correct periodic image
        const vec4_t ext = vec4_from_vec3(box, 0);
        const vec4_t ref = vec4_mul_f(ext, 0.5f);

        for (int64_t s_idx = 0; s_idx < num_structures; ++s_idx) {
            const int32_t* indices = md_index_range_beg(structures, s_idx);
            const int64_t num_indices = md_index_range_size(structures, s_idx);

            const vec4_t com = vec4_from_vec3(md_util_compute_com(x, y, z, w, indices, num_indices), 0);
            const vec4_t pbc_com = vec4_deperiodize(com, ref, ext);
            const vec4_t delta = vec4_sub(pbc_com, com);
            const vec4_t abs_delta = vec4_abs(delta);
            
            if (vec4_dot(delta, delta) > 0) {
                for (int64_t i = 0; i < num_indices; ++i) {
                    const int32_t j = indices[i];
                    if (abs_delta.x > 0) x[j] += delta.x;
                    if (abs_delta.y > 0) y[j] += delta.y;
                    if (abs_delta.z > 0) z[j] += delta.z;
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

    fifo_t queue = fifo_create(64, md_temp_allocator);
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
            md_array_push(xyz, pos, md_temp_allocator);
            md_array_push(pred, -1, md_temp_allocator);
        }
        
        //const vec4_t com = vec4_from_vec3(md_util_compute_com_vec3_ortho(xyz, NULL, md_array_size(xyz), pbc_ext), 0);
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

            const int* n_beg = md_index_range_beg(&covalent->atom_connectivity, idx);
            const int* n_end = md_index_range_end(&covalent->atom_connectivity, idx);
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

double md_util_compute_rmsd(const md_vec3_soa_t coord[2], const vec3_t com[2], const float* w, int64_t count) {
    const mat3_t R = mat3_optimal_rotation(coord[0].x, coord[0].y, coord[0].z, coord[1].x, coord[1].y, coord[1].z, w, com[0], com[1], count);
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

vec3_t md_util_shape_weights(const mat3_t* covariance_matrix) {
    ASSERT(covariance_matrix);
	mat3_eigen_t eigen = mat3_eigen(*covariance_matrix);
	float scl = 1.0f / (eigen.values.x + eigen.values.y + eigen.values.z);
	vec3_t weights = {(eigen.values.x - eigen.values.y) * scl, 2.0f * (eigen.values.y - eigen.values.z) * scl, 3.0f * eigen.values.z * scl};
	return weights;
}

bool md_util_linear_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[2], int64_t count, vec3_t pbc_ext, float t) {
    t = CLAMP(t, 0.0f, 1.0f);

    int64_t i = 0;
    const int64_t simd_count = (count / md_simd_width_f32) * md_simd_width_f32;
    for (; i < simd_count; i += md_simd_width_f32) {
        md_simd_f32_t x0 = md_simd_load_f32(src_coord[0].x + i);
        md_simd_f32_t y0 = md_simd_load_f32(src_coord[0].y + i);
        md_simd_f32_t z0 = md_simd_load_f32(src_coord[0].z + i);

        md_simd_f32_t x1 = md_simd_load_f32(src_coord[1].x + i);
        md_simd_f32_t y1 = md_simd_load_f32(src_coord[1].y + i);
        md_simd_f32_t z1 = md_simd_load_f32(src_coord[1].z + i);

        if (pbc_ext.x > 0.0f) {
            x1 = simd_deperiodize(x1, x0, md_simd_set1_f32(pbc_ext.x));
        }
        if (pbc_ext.y > 0.0f) {
            y1 = simd_deperiodize(y1, y0, md_simd_set1_f32(pbc_ext.y));
        }
        if (pbc_ext.z > 0.0f) {
            z1 = simd_deperiodize(z1, z0, md_simd_set1_f32(pbc_ext.z));
        }

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

    return true;
}

bool md_util_cubic_spline_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[4], int64_t count, vec3_t pbc_ext, float t, float s) {
    t = CLAMP(t, 0.0f, 1.0f);
    s = CLAMP(s, 0.0f, 1.0f);

    int64_t i = 0;
    const int64_t simd_count = (count / md_simd_width_f32) * md_simd_width_f32;
    for (; i < simd_count; i += md_simd_width_f32) {
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

        if (pbc_ext.x > 0.0f) {
            const md_simd_f32_t box_ext_x = md_simd_set1_f32(pbc_ext.x);
            x0 = simd_deperiodize(x0, x1, box_ext_x);
            x2 = simd_deperiodize(x2, x1, box_ext_x);
            x3 = simd_deperiodize(x3, x2, box_ext_x);
        }

        if (pbc_ext.y > 0.0f) {
            const md_simd_f32_t box_ext_y = md_simd_set1_f32(pbc_ext.y);
            y0 = simd_deperiodize(y0, y1, box_ext_y);
            y2 = simd_deperiodize(y2, y1, box_ext_y);
            y3 = simd_deperiodize(y3, y2, box_ext_y);
        }

        if (pbc_ext.z > 0.0f) {
            const md_simd_f32_t box_ext_z = md_simd_set1_f32(pbc_ext.z);
            z0 = simd_deperiodize(z0, z1, box_ext_z);
            z2 = simd_deperiodize(z2, z1, box_ext_z);
            z3 = simd_deperiodize(z3, z2, box_ext_z);
        }

        const md_simd_f32_t x = md_simd_cubic_spline(x0, x1, x2, x3, md_simd_set1_f32(t), md_simd_set1_f32(s));
        const md_simd_f32_t y = md_simd_cubic_spline(y0, y1, y2, y3, md_simd_set1_f32(t), md_simd_set1_f32(s));
        const md_simd_f32_t z = md_simd_cubic_spline(z0, z1, z2, z3, md_simd_set1_f32(t), md_simd_set1_f32(s));

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
        md_util_element_guess(mol->atom.element, mol->atom.count, mol);
    }

    if (flags & MD_UTIL_POSTPROCESS_RADIUS_BIT) {
        if (mol->atom.radius == 0)  md_array_resize(mol->atom.radius, mol->atom.count, alloc);
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            mol->atom.radius[i] = mol->atom.element ? md_util_element_vdw_radius(mol->atom.element[i]) : 1.0f;
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_MASS_BIT) {
        if (!mol->atom.mass) {
            md_array_resize(mol->atom.mass, mol->atom.count, alloc);
            for (int64_t i = 0; i < mol->atom.count; ++i) {
                mol->atom.mass[i] = mol->atom.element ? md_util_element_atomic_mass(mol->atom.element[i]) : 1.0f;
            }
        }
    }
   
    if (flags & MD_UTIL_POSTPROCESS_BOND_BIT) {
        if (!mol->bonds) {
            mol->bonds = md_util_compute_covalent_bonds(&mol->atom, &mol->unit_cell, alloc);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT) {
        if (mol->bonds) {
            mol->connectivity   = md_compute_connectivity(mol->bonds, md_array_size(mol->bonds), mol->atom.count, alloc);
            mol->structures     = md_util_compute_structures(mol->connectivity, alloc);
            mol->rings          = md_util_compute_rings(mol->connectivity, alloc);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_CHAINS_BIT) {
        if (mol->chain.count == 0 && mol->residue.count > 0 && mol->bonds) {
            md_util_compute_chain_data(&mol->chain, mol->atom.residue_idx, mol->atom.count, mol->bonds, md_array_size(mol->bonds), alloc);

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
                        md_array_push(backbone, atoms, md_heap_allocator);
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

            md_array_free(backbone, md_heap_allocator);

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

    md_allocator_i* alloc = md_heap_allocator;

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

    md_allocator_i* alloc = md_heap_allocator;

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
