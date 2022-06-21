#include "md_util.h"

#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_simd.h>

#include <ext/svd3/svd3.h>

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

#define bake(str) {str, ARRAY_SIZE(str)-1}

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
        if (compare_str_cstr(str, arr[i])) return i;
    }
    return -1;
}

static inline md_simd_typef simd_deperiodize(md_simd_typef x, md_simd_typef r, md_simd_typef period) {
    md_simd_typef d = md_simd_subf(x, r);
    md_simd_typef dx = md_simd_divf(d, period);
    dx = md_simd_subf(dx, md_simd_roundf(dx));
    return md_simd_addf(r, md_simd_mulf(dx, period));
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
            if (compare_str(str, sym)) return i;
        }
    }
    return 0;
}

static inline md_element_t lookup_element_ignore_case(str_t str) {
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

bool md_util_element_decode(md_element_t element[], int64_t capacity, const struct md_molecule_t* mol) {
    ASSERT(capacity >= 0);
    ASSERT(mol);
    ASSERT(mol->atom.count >= 0);

    const int64_t count = MIN(capacity, mol->atom.count);
    for (int64_t i = 0; i < count; ++i) {
        if (element[i] != 0) continue;

        const char* beg = mol->atom.name[i].buf;
        const char* end = mol->atom.name[i].buf + mol->atom.name[i].len;

        // Trim whitespace, digits and 'X's
        const char* c = beg;
        while (c < end && *c && (is_digit(*c) || is_whitespace(*c) || *c == 'x' || *c == 'X')) ++c;
        beg = c;

        c = beg;
        while (c < end && is_alpha(*c) && !(*c == 'x' || *c == 'X')) ++c;
        end = c;

        str_t name = {
            .ptr = beg,
            .len = end-beg
        };

        if (name.len > 0) {

            md_element_t elem = 0;
            if ((elem = md_util_element_lookup(name)) != 0) goto done;

            // If amino acid, try to deduce the element from that
            if (mol->atom.residue_idx && mol->residue.name) {
                str_t resname = trim_whitespace(label_to_str(&mol->residue.name[mol->atom.residue_idx[i]]));
                if (md_util_resname_amino_acid(resname)) {
                    // EASY-PEASY, we just try to match against the first character
                    name.len = 1;
                    elem = lookup_element_ignore_case(name);
                    goto done;
                }
            }

            // Try to match against several characters but ignore the case
            if (name.len > 1) {
                name.len = 2;
                elem = lookup_element_ignore_case(name);
            }

            // Last resort, try to match against single first character
            if (elem == 0) {
                name.len = 1;
                elem = lookup_element_ignore_case(name);
            }

        done:
            element[i] = elem;
        }
    }

    return true;
}

str_t md_util_element_symbol(md_element_t element) {
    ASSERT(element < Num_Elements);
    return element_symbols[element];
}

str_t md_util_element_name(md_element_t element) {
    ASSERT(element < Num_Elements);
    return element_names[element];
}

float md_util_element_vdw_radius(md_element_t element) {
    ASSERT(element < Num_Elements);
    return element_vdw_radii[element];
}

float md_util_element_covalent_radius(md_element_t element) {
    ASSERT(element < Num_Elements);
    return element_covalent_radii[element];
}

float md_util_element_atomic_mass(md_element_t element) {
    ASSERT(element < Num_Elements);
    return element_atomic_mass[element];
}

int md_util_element_max_valence(md_element_t element) {
    ASSERT(element < Num_Elements);
    return (int)element_max_valence[element];
}

uint32_t md_util_element_cpk_color(md_element_t element) {
    ASSERT(element < Num_Elements);
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
    if (bits & (2|4|8)) {
        if (backbone_atoms) *backbone_atoms = bb;
        return true;
    }
    return false;
}

bool md_util_backbone_atoms_extract(md_backbone_atoms_t backbone_atoms[], int64_t capacity, const md_molecule_t* mol) {
    if (!backbone_atoms) return false;
    if (capacity < 0) return false;
    memset(backbone_atoms, 0, capacity * sizeof(md_backbone_atoms_t));

    if (mol) {
        int64_t size = MIN(capacity, mol->backbone.count);
        for (int64_t bb_idx = 0; bb_idx < size; ++bb_idx) {
            const md_residue_idx_t res_idx = mol->backbone.residue_idx[bb_idx];
            if (!extract_backbone_atoms(&backbone_atoms[bb_idx], mol->atom.name, mol->residue.atom_range[res_idx])) {
                md_printf(MD_LOG_TYPE_INFO, "Failed to extract backbone atoms of residue[%i], possible that the residue is not an amino acid", res_idx);
            }
        }
        return true;
    }
    return false;
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

    memset(secondary_structure, 0, capacity * sizeof(md_secondary_structure_t));

    if (!mol) return false;
    if (!mol->atom.x) return false;
    if (!mol->atom.y) return false;
    if (!mol->atom.z) return false;
    if (!mol->backbone.atoms) return false;
    if (!mol->chain.backbone_range) return false;

    for (int64_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
        const md_range_t range = mol->chain.backbone_range[chain_idx];
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

    memset(backbone_angles, 0, sizeof(md_backbone_angles_t) * capacity);

    if (!mol) return false;
    if (!mol->atom.x) return false;
    if (!mol->atom.y) return false;
    if (!mol->atom.z) return false;
    if (!mol->backbone.atoms) return false;
    

    for (int64_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
        const md_range_t range = mol->chain.backbone_range[chain_idx];
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
    memset(ramachandran_types, MD_RAMACHANDRAN_TYPE_UNKNOWN, sizeof(md_ramachandran_type_t) * capacity);

    if (capacity < 0) return false;
    if (mol->backbone.count == 0) return false;
    if (mol->residue.count == 0) return false;

    ASSERT(mol->residue.name);
    ASSERT(mol->backbone.residue_idx);

    int64_t size = MIN(capacity, mol->backbone.count);

    for (int64_t i = 0; i < size; ++i) {
        int64_t res_idx = mol->backbone.residue_idx[i];
        ASSERT(res_idx < mol->residue.count);

        str_t resname = label_to_str(&mol->residue.name[res_idx]);
        if (compare_str_cstr_n(resname, "GLY", 3)) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GLYCINE;
        } else if (compare_str_cstr_n(resname, "PRO", 3)) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PROLINE;
            ramachandran_types[i - 1] = MD_RAMACHANDRAN_TYPE_PREPROL;
        } else {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GENERAL;
        }
    }

    return true;
}

static inline bool covelent_bond_heuristic(float dist_squared, const md_element_t elem[2]) {
    const float d = element_covalent_radii[elem[0]] + element_covalent_radii[elem[1]];
    const float d_min = d - 0.5f;
    const float d_max = d + 0.3f;
    return (d_min * d_min) < dist_squared && dist_squared < (d_max * d_max);
}

bool md_util_extract_covalent_bonds(md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    if (mol->covalent_bond.bond) {
        md_array_shrink(mol->covalent_bond.bond, 0);
    }

    md_array_resize(mol->atom.valence, mol->atom.count, alloc);
    memset(mol->atom.valence, 0, md_array_bytes(mol->atom.valence));

    const vec4_t period = vec4_from_vec3(md_util_compute_unit_cell_extent(mol->coord_frame), 0);

    if (mol->residue.count > 0) {
        // Residues are defined,
        // First find connections within residues, then between residues

        for (int64_t ri = 0; ri < mol->residue.count; ++ri) {
            const int64_t pre_internal_bond_count = md_array_size(mol->covalent_bond.bond);

            // Compute internal bonds (within residue)
            for (int64_t i = mol->residue.atom_range[ri].beg; i < mol->residue.atom_range[ri].end - 1; ++i) {
                const vec4_t a = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};
                for (int64_t j = i + 1; j < mol->residue.atom_range[ri].end; ++j) {
                    const vec4_t b = {mol->atom.x[j], mol->atom.y[j], mol->atom.z[j], 0};
                    const float d2 = vec4_periodic_distance_squared(a, b, period);
                    const md_element_t elem[2] = { mol->atom.element[i], mol->atom.element[j] };
                    if (covelent_bond_heuristic(d2, elem)) {
                        md_bond_t bond = {(int32_t)i, (int32_t)j};
                        md_array_push(mol->covalent_bond.bond, bond, alloc);
                        mol->atom.valence[i] += 1;
                        mol->atom.valence[j] += 1;
                    }
                }
            }

            const int64_t post_internal_bond_count = md_array_size(mol->covalent_bond.bond);

            if (mol->residue.internal_covalent_bond_range) {
                mol->residue.internal_covalent_bond_range[ri].beg = (uint32_t)pre_internal_bond_count;
                mol->residue.internal_covalent_bond_range[ri].end = (uint32_t)post_internal_bond_count;
            }

            const int64_t rj = ri + 1;
            if (rj < mol->residue.count) {
                // Compute external bonds (between residues)
                for (int64_t i = mol->residue.atom_range[ri].beg; i < mol->residue.atom_range[ri].end; ++i) {
                    const vec4_t a = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};
                    for (int64_t j = mol->residue.atom_range[rj].beg; j < mol->residue.atom_range[rj].end; ++j) {
                        const vec4_t b = {mol->atom.x[j], mol->atom.y[j], mol->atom.z[j], 0};
                        const float d2 = vec4_periodic_distance_squared(a, b, period);
                        const md_element_t elem[2] = { mol->atom.element[i], mol->atom.element[j] };
                        if (covelent_bond_heuristic(d2, elem)) {
                            md_bond_t bond = {(int32_t)i, (int32_t)j};
                            md_array_push(mol->covalent_bond.bond, bond, alloc);
                            mol->atom.valence[i] += 1;
                            mol->atom.valence[j] += 1;
                        }
                    }
                }

                const int64_t post_external_bond_count = md_array_size(mol->covalent_bond.bond);

                if (mol->residue.complete_covalent_bond_range) {
                    mol->residue.complete_covalent_bond_range[ri].end = (uint32_t)post_external_bond_count;
                    mol->residue.complete_covalent_bond_range[rj].beg = (uint32_t)post_internal_bond_count;
                }
            }
        }

        if (mol->residue.complete_covalent_bond_range) {
            ASSERT(mol->residue.count > 0);
            mol->residue.complete_covalent_bond_range[0].beg = 0;
            mol->residue.complete_covalent_bond_range[mol->residue.count - 1].end = (uint32_t)md_array_size(mol->covalent_bond.bond);
        }
    }
    else {
        // No residues present, test all against all
        // @TODO: This should be accelerated by spatial hashing
        for (int64_t i = 0; i < mol->atom.count - 1; ++i) {
            const vec4_t a = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};
            for (int64_t j = i + 1; j < mol->atom.count; ++j) {
                const vec4_t b = {mol->atom.x[j], mol->atom.y[j], mol->atom.z[j], 0};
                const float d2 = vec4_periodic_distance_squared(a, b, period);
                const md_element_t elem[2] = { mol->atom.element[i], mol->atom.element[j] };
                if (covelent_bond_heuristic(d2, elem)) {
                    md_bond_t bond = {(int32_t)i, (int32_t)j};
                    md_array_push(mol->covalent_bond.bond, bond, alloc);
                    mol->atom.valence[i] += 1;
                    mol->atom.valence[j] += 1;
                }
            }
        }
    }

    mol->covalent_bond.count = md_array_size(mol->covalent_bond.bond);

    return true;
}

vec3_t md_util_compute_com(const float* x, const float* y, const float* z, const float* w, int64_t count) {
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
            vec4_t xyz1 = {x[i], y[i], z[i], 1.0f};
            sum_xyzw = vec4_add(sum_xyzw, vec4_mul_f(xyz1, w[i]));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            vec4_t xyz0 = {x[i], y[i], z[i], 0};
            sum_xyzw = vec4_add(sum_xyzw, xyz0);
        }
        sum_xyzw.w = (float)count;
    }

    return vec3_div_f(vec3_from_vec4(sum_xyzw), sum_xyzw.w);
}

// We need to support cases where pbc_ext is zero, therefore we need to pick our algorithm based on this
// Single component version
static inline double compute_com_periodic(const float* in_x, const float* in_w, int64_t count, float pbc_ext) {
    double com = 0;

    if (pbc_ext > 0) {
        double acc_c = 0;
        double acc_s = 0;
        double acc_w = 0;

        const double scl = TWO_PI / (double)pbc_ext;

        for (int64_t i = 0; i < count; ++i) {
            double theta = in_x[i] * scl;
            double w = in_w ? in_w[i] : 1.0;
            acc_c += w * cos(theta);
            acc_s += w * sin(theta);
            acc_w += w;
        }

        const double theta_prim = atan2(-acc_s / acc_w, -acc_c / acc_w) + PI;
        com = (theta_prim / TWO_PI) * (double)pbc_ext;
    } else {
        double acc_x = 0;
        double acc_w = 0;
        for (int64_t i = 0; i < count; ++i) {
            double x = in_x[i];
            double w = in_w ? in_w[i] : 1.0;
            acc_x += x * w;
            acc_w += w;
        }
        com = acc_x / acc_w;
    }

    return com;
}

// Love the elegance of using trigonometric functions, unsure of the performance...
// @TODO: sin, cos and atan2 can and should of course be vectorized.
vec3_t md_util_compute_com_periodic(const float* in_x, const float* in_y, const float* in_z, const float* in_w, int64_t count, vec3_t pbc_ext) {
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

    double x = compute_com_periodic(in_x, in_w, count, pbc_ext.x);
    double y = compute_com_periodic(in_y, in_w, count, pbc_ext.y);
    double z = compute_com_periodic(in_z, in_w, count, pbc_ext.z);

    return (vec3_t) {(float)x, (float)y, (float)z};
}

// From here: https://www.arianarab.com/post/crystal-structure-software
mat3_t md_util_compute_unit_cell_basis(double a, double b, double c, double alpha, double beta, double gamma) {
    alpha = DEG_TO_RAD(alpha);
    beta  = DEG_TO_RAD(beta);
    gamma = DEG_TO_RAD(gamma);

    const double cb = cos(beta);
    const double x = (cos(alpha) - cb * cos(gamma)) / sin(gamma);
    mat3_t M = {
        .col = {
            {a, 0, 0},
            {b * cos(gamma), b * sin(gamma), 0},
            {c * cb, c * x, c * sqrt(1 - cb * cb - x * x)},
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

bool md_util_apply_pbc(md_molecule_t* mol, vec3_t pbc_ext) {
    if (!mol) return false;

    float* residue_com_x = 0;
    float* residue_com_y = 0;
    float* residue_com_z = 0;
    float* residue_aabb_min_x = 0;
    float* residue_aabb_max_x = 0;
    float* residue_aabb_min_y = 0;
    float* residue_aabb_max_y = 0;
    float* residue_aabb_min_z = 0;
    float* residue_aabb_max_z = 0;

    if (mol->chain.residue_range && mol->chain.count) {
        int64_t stride = ROUND_UP(mol->residue.count, md_simd_widthf);
        float* mem = md_alloc(default_temp_allocator, stride * sizeof(float) * 9);
        ASSERT(mem);
        residue_com_x       = mem + stride * 0;
        residue_com_y       = mem + stride * 1;
        residue_com_z       = mem + stride * 2;
        residue_aabb_min_x  = mem + stride * 3;
        residue_aabb_max_x  = mem + stride * 4;
        residue_aabb_min_y  = mem + stride * 5;
        residue_aabb_max_y  = mem + stride * 6;
        residue_aabb_min_z  = mem + stride * 7;
        residue_aabb_max_z  = mem + stride * 8;
    }
    
    for (int64_t i = 0; i < mol->residue.count; ++i) {
        md_range_t atom_range = mol->residue.atom_range[i];
        vec3_t com = md_util_compute_com_periodic(
            mol->atom.x + atom_range.beg,
            mol->atom.y + atom_range.beg,
            mol->atom.z + atom_range.beg,
            mol->atom.mass + atom_range.beg,
            atom_range.end - atom_range.beg,
            pbc_ext);

        //if (atom_range.end - atom_range.beg > 3) {
            for (int64_t j = atom_range.beg; j < atom_range.end; ++j) {
                mol->atom.x[j] = deperiodizef(mol->atom.x[j], com.x, pbc_ext.x);
                mol->atom.y[j] = deperiodizef(mol->atom.y[j], com.y, pbc_ext.y);
                mol->atom.z[j] = deperiodizef(mol->atom.z[j], com.z, pbc_ext.z);
            }
        //}

        if (residue_com_x) {
            ASSERT(residue_com_y);
            ASSERT(residue_com_z);
            ASSERT(residue_aabb_min_x);
            ASSERT(residue_aabb_max_x);
            ASSERT(residue_aabb_min_y);
            ASSERT(residue_aabb_max_y);
            ASSERT(residue_aabb_min_z);
            ASSERT(residue_aabb_max_z);
            residue_com_x[i] = com.x;
            residue_com_y[i] = com.y;
            residue_com_z[i] = com.z;
            residue_aabb_min_x[i] = +FLT_MAX;
            residue_aabb_max_x[i] = -FLT_MAX;
            residue_aabb_min_y[i] = +FLT_MAX;
            residue_aabb_max_y[i] = -FLT_MAX;
            residue_aabb_min_z[i] = +FLT_MAX;
            residue_aabb_max_z[i] = -FLT_MAX;

            for (int64_t j = atom_range.beg; j < atom_range.end; ++j) {
                residue_aabb_min_x[i] = MIN(residue_aabb_min_x[i], mol->atom.x[j]);
                residue_aabb_max_x[i] = MAX(residue_aabb_max_x[i], mol->atom.x[j]);
                residue_aabb_min_y[i] = MIN(residue_aabb_min_y[i], mol->atom.y[j]);
                residue_aabb_max_y[i] = MAX(residue_aabb_max_y[i], mol->atom.y[j]);
                residue_aabb_min_z[i] = MIN(residue_aabb_min_z[i], mol->atom.z[j]);
                residue_aabb_max_z[i] = MAX(residue_aabb_max_z[i], mol->atom.z[j]);
            }

            if (residue_aabb_max_x[i] - residue_aabb_min_x[i] > pbc_ext.x * 0.5f) {
                while(0) {};
            }
        }
    }

    for (int64_t i = 0; i < mol->chain.count; ++i) {
        md_range_t res_range = mol->chain.residue_range[i];
        vec3_t chain_com = md_util_compute_com_periodic(
            residue_com_x + res_range.beg,
            residue_com_y + res_range.beg,
            residue_com_z + res_range.beg,
            NULL,
            res_range.end - res_range.beg,
            pbc_ext);

        // Compute the aabb of the chain
        vec3_t chain_aabb_min = {+FLT_MAX, +FLT_MAX, +FLT_MAX};
        vec3_t chain_aabb_max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (int64_t j = res_range.beg; j < res_range.end; ++j) {
            vec3_t res_com = {residue_com_x[j], residue_com_y[j], residue_com_z[j]};

            // Ensure that the residue is within the period of the chain
            res_com.x = deperiodizef(res_com.x, chain_com.x, pbc_ext.x);
            res_com.y = deperiodizef(res_com.y, chain_com.y, pbc_ext.y);
            res_com.z = deperiodizef(res_com.z, chain_com.z, pbc_ext.z);

            vec3_t d = vec3_sub(res_com, (vec3_t){residue_com_x[j], residue_com_y[j], residue_com_z[j]});
            const float d2 = vec3_dot(d,d);

            if (d2 > 0) {
                residue_aabb_min_x[j] += d.x;
                residue_aabb_max_x[j] += d.x;
                residue_aabb_min_y[j] += d.y;
                residue_aabb_max_y[j] += d.y;
                residue_aabb_min_z[j] += d.z;
                residue_aabb_max_z[j] += d.z;
            }

            chain_aabb_min.x = MIN(chain_aabb_min.x, residue_aabb_min_x[j]);
            chain_aabb_max.x = MAX(chain_aabb_max.x, residue_aabb_max_x[j]);
            chain_aabb_min.y = MIN(chain_aabb_min.y, residue_aabb_min_y[j]);
            chain_aabb_max.y = MAX(chain_aabb_max.y, residue_aabb_max_y[j]);
            chain_aabb_min.z = MIN(chain_aabb_min.z, residue_aabb_min_z[j]);
            chain_aabb_max.z = MAX(chain_aabb_max.z, residue_aabb_max_z[j]);
        }

        // Create mask on which axis to apply the deperiodization
        // Don't apply deperiodization if the chain spans more than half the simulation extent in that dimension
        vec3_t chain_dp_mask = vec3_less_than(vec3_sub(chain_aabb_max, chain_aabb_min), vec3_mul_f(pbc_ext, 0.5f));

        if (!vec3_equal(chain_dp_mask, (vec3_t){0,0,0})) {
            for (int64_t j = res_range.beg; j < res_range.end; ++j) {
                vec3_t res_com = {residue_com_x[j], residue_com_y[j], residue_com_z[j]};

                // Ensure that the residue is within the period of the chain
                res_com.x = deperiodizef(res_com.x, chain_com.x, pbc_ext.x);
                res_com.y = deperiodizef(res_com.y, chain_com.y, pbc_ext.y);
                res_com.z = deperiodizef(res_com.z, chain_com.z, pbc_ext.z);

                vec3_t d = vec3_sub(res_com, (vec3_t){residue_com_x[j], residue_com_y[j], residue_com_z[j]});
                d = vec3_mul(d, chain_dp_mask);
                const float d2 = vec3_dot(d,d);

                if (d2 > 0) {
                    md_range_t atom_range = mol->residue.atom_range[j];
                    for (int64_t k = atom_range.beg; k < atom_range.end; ++k) {

                        if (strncmp(mol->atom.name->buf, "HW1", 3) == 0) {
                            while(0){};
                        }
                        mol->atom.x[k] += d.x;
                        mol->atom.y[k] += d.y;
                        mol->atom.z[k] += d.z;
                    }
                }            
            }
        }
    }

    return true;
}

/*
vec3_t md_util_compute_periodic_com(const float* in_x, const float* in_y, const float* in_z, const float* in_w, int64_t count, vec3_t pbc_ext) {
    float w = in_w ? in_w[0] : 1.0f;
    float sum_x = in_x[0] * w;
    float sum_y = in_y[0] * w;
    float sum_z = in_z[0] * w;
    float sum_w = w;

    for (int64_t i = 1; i < count; ++i) {
        float com_x = sum_x / sum_w;
        float com_y = sum_y / sum_w;
        float com_z = sum_z / sum_w;

        float x = deperiodizef(in_x[i], com_x, pbc_ext.x);
        float y = deperiodizef(in_y[i], com_y, pbc_ext.y);
        float z = deperiodizef(in_z[i], com_z, pbc_ext.z);

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
        md_simd_typef box_ext_x = md_simd_set1f(pbc_ext.x);
        md_simd_typef box_ext_y = md_simd_set1f(pbc_ext.y);
        md_simd_typef box_ext_z = md_simd_set1f(pbc_ext.z);

        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_widthf) * md_simd_widthf;
        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef x0 = md_simd_loadf(src_coord[0].x + i);
            md_simd_typef y0 = md_simd_loadf(src_coord[0].y + i);
            md_simd_typef z0 = md_simd_loadf(src_coord[0].z + i);

            md_simd_typef x1 = md_simd_loadf(src_coord[1].x + i);
            md_simd_typef y1 = md_simd_loadf(src_coord[1].y + i);
            md_simd_typef z1 = md_simd_loadf(src_coord[1].z + i);

            x1 = simd_deperiodize(x1, x0, box_ext_x);
            y1 = simd_deperiodize(y1, y0, box_ext_y);
            z1 = simd_deperiodize(z1, z0, box_ext_z);

            md_simd_typef x = md_simd_lerpf(x0, x1, t);
            md_simd_typef y = md_simd_lerpf(y0, y1, t);
            md_simd_typef z = md_simd_lerpf(z0, z1, t);

            md_simd_storef(dst_coord.x + i, x);
            md_simd_storef(dst_coord.y + i, y);
            md_simd_storef(dst_coord.z + i, z);
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
        const int64_t simd_count = (count / md_simd_widthf) * md_simd_widthf;
        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef x0 = md_simd_loadf(src_coord[0].x + i);
            md_simd_typef y0 = md_simd_loadf(src_coord[0].y + i);
            md_simd_typef z0 = md_simd_loadf(src_coord[0].z + i);

            md_simd_typef x1 = md_simd_loadf(src_coord[1].x + i);
            md_simd_typef y1 = md_simd_loadf(src_coord[1].y + i);
            md_simd_typef z1 = md_simd_loadf(src_coord[1].z + i);

            md_simd_typef x = md_simd_lerpf(x0, x1, t);
            md_simd_typef y = md_simd_lerpf(y0, y1, t);
            md_simd_typef z = md_simd_lerpf(z0, z1, t);

            md_simd_storef(dst_coord.x + i, x);
            md_simd_storef(dst_coord.y + i, y);
            md_simd_storef(dst_coord.z + i, z);
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

bool md_util_cubic_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[4], int64_t count, vec3_t pbc_ext, float t, float tension) {
    const bool use_pbc = !vec3_equal(pbc_ext, (vec3_t){0,0,0});
    t = CLAMP(t, 0.0f, 1.0f);
    tension = CLAMP(tension, 0.0f, 1.0f);

    if (use_pbc) {
        md_simd_typef box_ext_x = md_simd_set1f(pbc_ext.x);
        md_simd_typef box_ext_y = md_simd_set1f(pbc_ext.y);
        md_simd_typef box_ext_z = md_simd_set1f(pbc_ext.z);

        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_widthf) * md_simd_widthf;
        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef x0 = md_simd_loadf(src_coord[0].x + i);
            md_simd_typef y0 = md_simd_loadf(src_coord[0].y + i);
            md_simd_typef z0 = md_simd_loadf(src_coord[0].z + i);

            md_simd_typef x1 = md_simd_loadf(src_coord[1].x + i);
            md_simd_typef y1 = md_simd_loadf(src_coord[1].y + i);
            md_simd_typef z1 = md_simd_loadf(src_coord[1].z + i);

            md_simd_typef x2 = md_simd_loadf(src_coord[2].x + i);
            md_simd_typef y2 = md_simd_loadf(src_coord[2].y + i);
            md_simd_typef z2 = md_simd_loadf(src_coord[2].z + i);

            md_simd_typef x3 = md_simd_loadf(src_coord[3].x + i);
            md_simd_typef y3 = md_simd_loadf(src_coord[3].y + i);
            md_simd_typef z3 = md_simd_loadf(src_coord[3].z + i);

            x0 = simd_deperiodize(x0, x1, box_ext_x);
            x2 = simd_deperiodize(x2, x1, box_ext_x);
            x3 = simd_deperiodize(x3, x2, box_ext_x);

            y0 = simd_deperiodize(y0, y1, box_ext_y);
            y2 = simd_deperiodize(y2, y1, box_ext_y);
            y3 = simd_deperiodize(y3, y2, box_ext_y);

            z0 = simd_deperiodize(z0, z1, box_ext_z);
            z2 = simd_deperiodize(z2, z1, box_ext_z);
            z3 = simd_deperiodize(z3, z2, box_ext_z);

            md_simd_typef x = md_simd_cubic_splinef(x0, x1, x2, x3, t, tension);
            md_simd_typef y = md_simd_cubic_splinef(y0, y1, y2, y3, t, tension);
            md_simd_typef z = md_simd_cubic_splinef(z0, z1, z2, z3, t, tension);

            md_simd_storef(dst_coord.x + i, x);
            md_simd_storef(dst_coord.y + i, y);
            md_simd_storef(dst_coord.z + i, z);
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

            const vec4_t coord = vec4_cubic_spline(src[0], src[1], src[2], src[3], t, tension);

            dst_coord.x[i] = coord.x;
            dst_coord.y[i] = coord.y;
            dst_coord.z[i] = coord.z;
        }
    } else {
        int64_t i = 0;
        const int64_t simd_count = (count / md_simd_widthf) * md_simd_widthf;
        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef x0 = md_simd_loadf(src_coord[0].x + i);
            md_simd_typef y0 = md_simd_loadf(src_coord[0].y + i);
            md_simd_typef z0 = md_simd_loadf(src_coord[0].z + i);

            md_simd_typef x1 = md_simd_loadf(src_coord[1].x + i);
            md_simd_typef y1 = md_simd_loadf(src_coord[1].y + i);
            md_simd_typef z1 = md_simd_loadf(src_coord[1].z + i);

            md_simd_typef x2 = md_simd_loadf(src_coord[2].x + i);
            md_simd_typef y2 = md_simd_loadf(src_coord[2].y + i);
            md_simd_typef z2 = md_simd_loadf(src_coord[2].z + i);

            md_simd_typef x3 = md_simd_loadf(src_coord[3].x + i);
            md_simd_typef y3 = md_simd_loadf(src_coord[3].y + i);
            md_simd_typef z3 = md_simd_loadf(src_coord[3].z + i);

            md_simd_typef x = md_simd_cubic_splinef(x0, x1, x2, x3, t, tension);
            md_simd_typef y = md_simd_cubic_splinef(y0, y1, y2, y3, t, tension);
            md_simd_typef z = md_simd_cubic_splinef(z0, z1, z2, z3, t, tension);

            md_simd_storef(dst_coord.x + i, x);
            md_simd_storef(dst_coord.y + i, y);
            md_simd_storef(dst_coord.z + i, z);
        }

        // Do the rest
        for (; i < count; ++i) {
            vec4_t src[4] = {
                {src_coord[0].x[i], src_coord[0].y[i], src_coord[0].z[i], 1},
                {src_coord[1].x[i], src_coord[1].y[i], src_coord[1].z[i], 1},
                {src_coord[2].x[i], src_coord[2].y[i], src_coord[2].z[i], 1},
                {src_coord[3].x[i], src_coord[3].y[i], src_coord[3].z[i], 1},
            };

            vec4_t coord = vec4_cubic_spline(src[0], src[1], src[2], src[3], t, tension);

            dst_coord.x[i] = coord.x;
            dst_coord.y[i] = coord.y;
            dst_coord.z[i] = coord.z;
        }
    }

    return true;
}

static inline bool ranges_overlap(md_range_t a, md_range_t b) {
    return (a.beg < b.end&& b.beg < a.end);
}

// Try to fill in missing fields for molecule struct
// Element -> Mass / Radius
// Bonds
bool md_util_postprocess_molecule(struct md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    if (mol->atom.element == 0) {
        md_array_resize(mol->atom.element, mol->atom.count, alloc);
        memset(mol->atom.element, 0, mol->atom.count * sizeof(*mol->atom.element));
    }
    md_util_element_decode(mol->atom.element, mol->atom.count, mol);

    if (mol->atom.radius == 0)  md_array_resize(mol->atom.radius, mol->atom.count, alloc);
    if (mol->atom.mass == 0)    md_array_resize(mol->atom.mass, mol->atom.count, alloc);
    if (mol->atom.valence == 0) {
        md_array_resize(mol->atom.valence, mol->atom.count, alloc);
        memset(mol->atom.valence, 0, md_array_size(mol->atom.valence) * sizeof(*mol->atom.valence));
    }
    if (mol->atom.flags == 0) {
        md_array_resize(mol->atom.flags, mol->atom.count, alloc);
        memset(mol->atom.flags, 0, md_array_size(mol->atom.flags) * sizeof(*mol->atom.flags));
    }
    if (mol->atom.vx == 0) {
        md_array_resize(mol->atom.vx, mol->atom.count, alloc);
        memset(mol->atom.vx, 0, md_array_size(mol->atom.vx) * sizeof(*mol->atom.vx));
    }
    if (mol->atom.vy == 0) {
        md_array_resize(mol->atom.vy, mol->atom.count, alloc);
        memset(mol->atom.vy, 0, md_array_size(mol->atom.vy) * sizeof(*mol->atom.vy));
    }
    if (mol->atom.vz == 0) {
        md_array_resize(mol->atom.vz, mol->atom.count, alloc);
        memset(mol->atom.vz, 0, md_array_size(mol->atom.vz) * sizeof(*mol->atom.vz));
    }
   
    for (int64_t i = 0; i < mol->atom.count; ++i) {
        mol->atom.radius[i] = md_util_element_vdw_radius(mol->atom.element[i]);
        mol->atom.mass[i]   = md_util_element_atomic_mass(mol->atom.element[i]);
    }

    if (mol->covalent_bond.count == 0) {
        md_array_ensure(mol->residue.internal_covalent_bond_range, mol->residue.count, alloc);
        md_array_ensure(mol->residue.complete_covalent_bond_range, mol->residue.count, alloc);

        // Use heuristical method of finding covalent bonds
        md_util_extract_covalent_bonds(mol, alloc);
    }

    if (mol->chain.count == 0 && mol->residue.count > 0) {
        // Compute artificial chains and label them A - Z
        if (mol->residue.complete_covalent_bond_range) {
            // Identify connected residues (through covalent bonds) if more than one, we store it as a chain
            md_range_t res_range = { 0, 1 };
            for (int64_t res_idx = 0; res_idx < mol->residue.count - 1; ++res_idx) {
                if (ranges_overlap(mol->residue.complete_covalent_bond_range[res_idx],
                    mol->residue.complete_covalent_bond_range[res_idx + 1])) {
                    res_range.end += 1;
                }
                else {
                    if (res_range.end - res_range.beg > 1) {
                        md_array_push(mol->chain.residue_range, res_range, alloc);
                    }
                    res_range.beg = (int32_t)res_idx + 1;
                    res_range.end = (int32_t)res_idx + 2;
                }
            }

            if (res_range.beg == 0 && res_range.end > 1) {
                md_array_push(mol->chain.residue_range, res_range, alloc);
            }

            mol->chain.count = md_array_size(mol->chain.residue_range);
            if (mol->chain.count) {
                const char* id_arr[] = { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z" };

                // Add labels
                md_array_resize(mol->chain.id, mol->chain.count, alloc);
                for (int64_t i = 0; i < mol->chain.count; ++i) {
                    str_t id = { id_arr[i % ARRAY_SIZE(id_arr)], 1 };
                    mol->chain.id[i] = make_label(id);
                }

                // Add atom ranges
                md_array_resize(mol->chain.atom_range, mol->chain.count, alloc);
                for (int64_t i = 0; i < mol->chain.count; ++i) {
                    md_range_t range = mol->chain.residue_range[i];
                    mol->chain.atom_range[i].beg = mol->residue.atom_range[range.beg].beg;
                    mol->chain.atom_range[i].end = mol->residue.atom_range[range.end - 1].end;
                }

                // Set atom chain indices
                md_array_resize(mol->atom.chain_idx, mol->atom.count, alloc);
                for (int64_t i = 0; i < mol->atom.count; ++i) {
                    mol->atom.chain_idx[i] = -1;
                }

                for (int64_t ci = 0; ci < mol->chain.count; ++ci) {
                    for (int64_t ai = (int64_t)mol->chain.atom_range[ci].beg; ai < (int64_t)mol->chain.atom_range[ci].end; ++ai) {
                        mol->atom.chain_idx[ai] = (md_chain_idx_t)ci;
                    }
                }
            }
        }
    }

    if (mol->chain.count) {
        // Compute backbone data
        int64_t bb_offset = 0;
        for (int64_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
            const int64_t res_count = mol->chain.residue_range[chain_idx].end - mol->chain.residue_range[chain_idx].beg;
            md_range_t bb_range = { (int32_t)bb_offset, (int32_t)(bb_offset + res_count) };
            md_array_push(mol->chain.backbone_range, bb_range, alloc);
            bb_offset += res_count;
        }

        mol->backbone.count = bb_offset;
        md_array_resize(mol->backbone.atoms, mol->backbone.count, alloc);
        md_array_resize(mol->backbone.angle, mol->backbone.count, alloc);
        md_array_resize(mol->backbone.secondary_structure, mol->backbone.count, alloc);
        md_array_resize(mol->backbone.ramachandran_type, mol->backbone.count, alloc);
        md_array_resize(mol->backbone.residue_idx, mol->backbone.count, alloc);

        // Set the residue indices for the backbone
        bb_offset = 0;
        for (int64_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
            int64_t res_count = (int64_t)mol->chain.residue_range[chain_idx].end - (int64_t)mol->chain.residue_range[chain_idx].beg;
            for (int64_t j = 0; j < res_count; ++j) {
                mol->backbone.residue_idx[bb_offset + j] = (md_residue_idx_t)((int64_t)mol->chain.residue_range[chain_idx].beg + j);
            }
            bb_offset += res_count;
        }

        md_util_backbone_atoms_extract(mol->backbone.atoms, mol->backbone.count, mol);
        md_util_backbone_secondary_structure_compute(mol->backbone.secondary_structure, mol->backbone.count, mol);
        md_util_backbone_ramachandran_classify(mol->backbone.ramachandran_type, mol->backbone.count, mol);
    }

    return true;
}

#ifdef __cplusplus
}
#endif
