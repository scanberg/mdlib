#include <md_util.h>

#include <md_molecule.h>
#include <md_smiles.h>

#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_intrinsics.h>
#include <core/md_simd.h>
#include <core/md_spatial_hash.h>
#include <core/md_bitfield.h>
#include <core/md_str_builder.h>
#include <core/md_os.h>
#include <core/md_hash.h>

#include <math.h>
#include <string.h>
#include <float.h>

//#define PROFILE

#define BAKE(str) {str, sizeof(str)-1}

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
    BAKE("Xx"), BAKE("H"),  BAKE("He"), BAKE("Li"), BAKE("Be"), BAKE("B"),  BAKE("C"),  BAKE("N"),  BAKE("O"),  BAKE("F"),  BAKE("Ne"), BAKE("Na"), BAKE("Mg"), BAKE("Al"), BAKE("Si"), BAKE("P"),  BAKE("S"),  BAKE("Cl"), BAKE("Ar"), BAKE("K"),  BAKE("Ca"), BAKE("Sc"), BAKE("Ti"), BAKE("V"),
    BAKE("Cr"), BAKE("Mn"), BAKE("Fe"), BAKE("Co"), BAKE("Ni"), BAKE("Cu"), BAKE("Zn"), BAKE("Ga"), BAKE("Ge"), BAKE("As"), BAKE("Se"), BAKE("Br"), BAKE("Kr"), BAKE("Rb"), BAKE("Sr"), BAKE("Y"),  BAKE("Zr"), BAKE("Nb"), BAKE("Mo"), BAKE("Tc"), BAKE("Ru"), BAKE("Rh"), BAKE("Pd"), BAKE("Ag"),
    BAKE("Cd"), BAKE("In"), BAKE("Sn"), BAKE("Sb"), BAKE("Te"), BAKE("I"),  BAKE("Xe"), BAKE("Cs"), BAKE("Ba"), BAKE("La"), BAKE("Ce"), BAKE("Pr"), BAKE("Nd"), BAKE("Pm"), BAKE("Sm"), BAKE("Eu"), BAKE("Gd"), BAKE("Tb"), BAKE("Dy"), BAKE("Ho"), BAKE("Er"), BAKE("Tm"), BAKE("Yb"), BAKE("Lu"),
    BAKE("Hf"), BAKE("Ta"), BAKE("W"),  BAKE("Re"), BAKE("Os"), BAKE("Ir"), BAKE("Pt"), BAKE("Au"), BAKE("Hg"), BAKE("Tl"), BAKE("Pb"), BAKE("Bi"), BAKE("Po"), BAKE("At"), BAKE("Rn"), BAKE("Fr"), BAKE("Ra"), BAKE("Ac"), BAKE("Th"), BAKE("Pa"), BAKE("U"),  BAKE("Np"), BAKE("Pu"), BAKE("Am"),
    BAKE("Cm"), BAKE("Bk"), BAKE("Cf"), BAKE("Es"), BAKE("Fm"), BAKE("Md"), BAKE("No"), BAKE("Lr"), BAKE("Rf"), BAKE("Db"), BAKE("Sg"), BAKE("Bh"), BAKE("Hs"), BAKE("Mt"), BAKE("Ds"), BAKE("Rg"), BAKE("Cn"), BAKE("Nh"), BAKE("Fl"), BAKE("Mc"), BAKE("Lv"), BAKE("Ts"), BAKE("Og"),
};


static const str_t element_names[] = {
    BAKE("Unknown"),     BAKE("Hydrogen"),     BAKE("Helium"),       BAKE("Lithium"),     BAKE("Beryllium"),   BAKE("Boron"),         BAKE("Carbon"),     BAKE("Nitrogen"),   BAKE("Oxygen"),
    BAKE("Fluorine"),    BAKE("Neon"),         BAKE("Sodium"),       BAKE("Magnesium"),   BAKE("Aluminium"),   BAKE("Silicon"),       BAKE("Phosphorus"), BAKE("Sulfur"),     BAKE("Chlorine"),
    BAKE("Argon"),       BAKE("Potassium"),    BAKE("Calcium"),      BAKE("Scandium"),    BAKE("Titanium"),    BAKE("Vanadium"),      BAKE("Chromium"),   BAKE("Manganese"),  BAKE("Iron"),
    BAKE("Cobalt"),      BAKE("Nickel"),       BAKE("Copper"),       BAKE("Zinc"),        BAKE("Gallium"),     BAKE("Germanium"),     BAKE("Arsenic"),    BAKE("Selenium"),   BAKE("Bromine"),
    BAKE("Krypton"),     BAKE("Rubidium"),     BAKE("Strontium"),    BAKE("Yttrium"),     BAKE("Zirconium"),   BAKE("Niobium"),       BAKE("Molybdenum"), BAKE("Technetium"), BAKE("Ruthenium"),
    BAKE("Rhodium"),     BAKE("Palladium"),    BAKE("Silver"),       BAKE("Cadmium"),     BAKE("Indium"),      BAKE("Tin"),           BAKE("Antimony"),   BAKE("Tellurium"),  BAKE("Iodine"),
    BAKE("Xenon"),       BAKE("Caesium"),      BAKE("Barium"),       BAKE("Lanthanum"),   BAKE("Cerium"),      BAKE("Praseodymium"),  BAKE("Neodymium"),  BAKE("Promethium"), BAKE("Samarium"),
    BAKE("Europium"),    BAKE("Gadolinium"),   BAKE("Terbium"),      BAKE("Dysprosium"),  BAKE("Holmium"),     BAKE("Erbium"),        BAKE("Thulium"),    BAKE("Ytterbium"),  BAKE("Lutetium"),
    BAKE("Hafnium"),     BAKE("Tantalum"),     BAKE("Tungsten"),     BAKE("Rhenium"),     BAKE("Osmium"),      BAKE("Iridium"),       BAKE("Platinum"),   BAKE("Gold"),       BAKE("Mercury"),
    BAKE("Thallium"),    BAKE("Lead"),         BAKE("Bismuth"),      BAKE("Polonium"),    BAKE("Astatine"),    BAKE("Radon"),         BAKE("Francium"),   BAKE("Radium"),     BAKE("Actinium"),
    BAKE("Thorium"),     BAKE("Protactinium"), BAKE("Uranium"),      BAKE("Neptunium"),   BAKE("Plutonium"),   BAKE("Americium"),     BAKE("Curium"),     BAKE("Berkelium"),  BAKE("Californium"),
    BAKE("Einsteinium"), BAKE("Fermium"),      BAKE("Mendelevium"),  BAKE("Nobelium"),    BAKE("Lawrencium"),  BAKE("Rutherfordium"), BAKE("Dubnium"),    BAKE("Seaborgium"), BAKE("Bohrium"),
    BAKE("Hassium"),     BAKE("Meitnerium"),   BAKE("Darmstadtium"), BAKE("Roentgenium"), BAKE("Copernicium"), BAKE("Nihonium"),      BAKE("Flerovium"),  BAKE("Moscovium"),  BAKE("Livermorium"),
    BAKE("Tennessine"),  BAKE("Oganesson"),                          
};

// http://dx.doi.org/10.1039/b801115j
static const uint8_t element_covalent_radii_u8[] = {
      0,  31,  28, 128,  96,  84,  76,  71,  66,  57,  58, 166, 141, 121, 111, 107, 105, 102, 106, 203, 176, 170, 160, 153,
    139, 139, 132, 126, 124, 132, 122, 122, 120, 119, 120, 120, 116, 220, 195, 190, 175, 164, 154, 147, 146, 142, 139, 145,
    144, 142, 139, 139, 138, 139, 140, 244, 215, 207, 204, 203, 201, 199, 198, 198, 196, 194, 192, 192, 189, 190, 187, 187,
    175, 170, 162, 151, 144, 141, 136, 136, 132, 145, 146, 148, 140, 150, 150, 255, 221, 215, 206, 200, 196, 190, 187, 180,
    169, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160
};

static const float element_covalent_radii_f32[] = {
       0, 0.31, 0.28, 1.28,  .96,  .84,  .76,  .71,  .66,  .57,  .58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53,
    1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45,
    1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87,
    1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80,
    1.69, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60
};

// Covalent radii for elements (single, double, triple) bonds in pikometers
// https://en.wikipedia.org/wiki/Covalent_radius
static const uint8_t element_covalent_radii3[][3] = {
    {  0,    0,   0 }, // Unknown
    {  32,   0,   0 }, // H
    {  46,   0,   0 }, // He
    { 133, 124,   0 }, // Li
    { 102,  90,  85 }, // Be
    {  85,  78,  73 }, // B
    {  75,  67,  60 }, // C
    {  71,  60,  54 }, // N
    {  63,  57,  53 }, // O
    {  64,  59,  53 }, // F
    {  67,  96,   0 }, // Ne
    { 155, 160,   0 }, // Na
    { 139, 132, 127 }, // Mg
    { 126, 113, 111 }, // Al
    { 116, 107, 102 }, // Si
    { 111, 102,  94 }, // P
    { 103,  94,  95 }, // S
    {  99,  95,  93 }, // Cl
    {  96, 107,  96 }, // Ar

    {  19, 196, 193 }, // K
    { 171, 147, 133 }, // Ca
    { 148, 116, 114 }, // Sc
    { 136, 117, 108 }, // Ti
    { 134, 112, 106 }, // V
    { 122, 111, 103 }, // Cr
    { 119, 105, 103 }, // Mn
    { 116, 109, 102 }, // Fe
    { 111, 103,  96 }, // Co
    { 110, 101, 101 }, // Ni
    { 112, 115, 120 }, // Cu
    { 118, 120,   0 }, // Zn
    { 124, 117, 121 }, // Ga
    { 121, 111, 114 }, // Ge
    { 121, 114, 106 }, // As
    { 116, 107, 107 }, // Se
    { 114, 109, 110 }, // Br
    { 117, 121, 108 }, // Kr

    { 210, 202,   0 }, // Rb
    { 185, 157, 139 }, // Sr
    { 163, 130, 124 }, // Y
    { 154, 127, 121 }, // Zr
    { 147, 125, 116 }, // Nb
    { 138, 121, 113 }, // Mo
    { 128, 120, 110 }, // Tc
    { 125, 114, 103 }, // Ru
    { 125, 110, 106 }, // Rh
    { 120, 117, 112 }, // Pd
    { 128, 139, 137 }, // Ag
    { 136, 144,   0 }, // Cd
    { 142, 136, 146 }, // In
    { 140, 130, 132 }, // Sn
    { 140, 133, 127 }, // Sb
    { 136, 128, 121 }, // Te
    { 133, 129, 121 }, // I
    { 131, 135, 122 }, // Xe
    { 232, 209,   0 }, // Cs
    { 196, 161, 149 }, // Ba
    { 180, 139, 139 }, // La
    { 163, 137, 131 }, // Ce
    { 176, 138, 128 }, // Pr
    { 174, 137,   0 }, // Nd
    { 173, 135,   0 }, // Pm
    { 172, 134,   0 }, // Sm
    { 168, 134,   0 }, // Eu
    { 169, 135, 132 }, // Gd
    { 168, 135,   0 }, // Tb
    { 167, 133,   0 }, // Dy
    { 166, 133,   0 }, // Ho
    { 165, 133,   0 }, // Er
    { 164, 131,   0 }, // Tm
    { 170, 129,   0 }, // Yb

    { 162, 131, 131 }, // Lu
    { 152, 128, 122 }, // Hf
    { 146, 126, 119 }, // Ta
    { 137, 120, 115 }, // W
    { 131, 119, 110 }, // Re
    { 129, 116, 109 }, // Os
    { 122, 115, 107 }, // Ir
    { 123, 112, 110 }, // Pt
    { 124, 121, 123 }, // Au
    { 133, 142,   0 }, // Hg
    { 144, 142, 150 }, // Tl
    { 144, 135, 137 }, // Pb
    { 151, 141, 135 }, // Bi
    { 145, 135, 129 }, // Po
    { 147, 138, 138 }, // At
    { 142, 145, 133 }, // Rn
    { 223, 218,   0 }, // Fr
    { 201, 173, 159 }, // Ra

    { 186, 153, 140 }, // Ac
    { 175, 143, 136 }, // Th
    { 169, 138, 129 }, // Pa
    { 170, 134, 118 }, // U
    { 171, 136, 116 }, // Np
    { 172, 135,   0 }, // Pu
    { 166, 135,   0 }, // Am
    { 166, 136,   0 }, // Cm
    { 168, 139,   0 }, // Bk
    { 168, 140,   0 }, // Cf
    { 165, 140,   0 }, // Es
    { 167,   0,   0 }, // Fm
    { 173, 139,   0 }, // Md
    { 176,   0,   0 }, // No

    { 161, 141,   0 }, // Lr
    { 157, 140, 131 }, // Rf
    { 149, 136, 126 }, // Db
    { 143, 128, 121 }, // Sg
    { 141, 128, 119 }, // Bh
    { 134, 125, 118 }, // Hs
    { 129, 125, 113 }, // Mt
    { 128, 116, 112 }, // Ds
    { 121, 116, 112 }, // Rg
    { 122, 137, 130 }, // Cn
    { 136,   0,   0 }, // Nh
    { 143,   0,   0 }, // Fl
    { 162,   0,   0 }, // Mc
    { 175,   0,   0 }, // Lv
    { 165,   0,   0 }, // Ts
    { 157,   0,   0 }, // Og
};

// https://dx.doi.org/10.1021/jp8111556
static const float element_vdw_radii[] = {
    1.00, 1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2.30, 2.15, 2.05,
    2.05, 2.05, 2.05, 2.00, 2.00, 2.00, 2.10, 1.87, 2.11, 1.85, 1.90, 1.83, 2.02, 3.03, 2.49, 2.40, 2.30, 2.15, 2.10, 2.05, 2.05, 2.00, 2.05, 2.10,
    2.20, 2.20, 1.93, 2.17, 2.06, 1.98, 2.16, 3.43, 2.68, 2.50, 2.48, 2.47, 2.45, 2.43, 2.42, 2.40, 2.38, 2.37, 2.35, 2.33, 2.32, 2.30, 2.28, 2.27,
    2.25, 2.20, 2.10, 2.05, 2.00, 2.00, 2.05, 2.10, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, 2.00, 2.30, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00
};

// http://chemistry.wikia.com/wiki/List_of_elements_by_atomic_mass
// last element padded with 295
static const float element_atomic_mass[] = {
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
static const uint8_t element_max_valence[] = {
    0,
    1, 0,
    1, 2, 3, 4, 5, 2, 1, 0,
    1, 2, 3, 4, 5, 6, 7, 0,
    1, 2, 3, 4, 5, 6, 7, 7, 5, 4, 4, 2, 3, 4, 5, 6, 7, 2,
    1, 2, 3, 4, 5, 6, 7, 8, 6, 4, 3, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 4, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 6, 5, 2, 3, 4, 5, 6, 7, 6,
    1, 2, 3, 4, 5, 6, 7, 7, 7, 6, 5, 5, 4, 3, 3, 3, 3, 4, 5, 6, 7, 8, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
};

#define RGBA(r,g,b,a) ( ((a & 255) << 24) | ((b & 255) << 16) | ((g & 255) << 8) | (r & 255) )
#define RGB(r,g,b) RGBA(r,g,b,255)

// Based on this with some values modified
// Some stem from Jmol, some from Rasmol
// http://jmol.sourceforge.net/jscolors/
static const uint32_t element_cpk_colors[] = {
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

// This is bitmasks expressed as combinations of two uint64_t that encodes properties of atomic elements
static const uint64_t element_metal_mask[2] = {
    0xff87ffe0fff83818,
    0x00001fffff9fffff,
};

static const uint64_t element_halogen_mask[2] = {
    0x0020000800020200,
    0x0000000000200000,
};

static const uint64_t element_alkali_mask[2] = {
    0x0080002000080808,
    0x0000000000800000,
};

static const uint64_t element_alkaline_earth[2] = {
    0x0100004000101010,
    0x0000000001000000,
};

static const uint64_t element_aromatic_ring[2] = {
    0x000c00030001c1e0,
    0x0000000000080000,
};

// This is a macro to generate test functions of the bitmasks above
// This seems to generate decent assembly across Clang, GCC and MSVC
// If this is a static inline function, it fails to represent the raw values of the mask in assembly
#define GENERATE_MASK_FUNC(name, mask) \
static inline bool name(uint32_t elem) { \
    uint64_t off = 0; \
    uint64_t pat = mask[0]; \
    if (elem >= 64) { \
        off = 64; \
        pat = mask[1]; \
    } \
    uint64_t msk = (elem < 119) ? 0xFFFFFFFFULL : 0; \
    uint64_t idx = (elem - off) & msk; \
    return pat & (1ULL << idx); \
}

GENERATE_MASK_FUNC(is_metal, element_metal_mask)

GENERATE_MASK_FUNC(is_halogen, element_halogen_mask)

GENERATE_MASK_FUNC(is_aromatic, element_aromatic_ring)

// Some from here https://daylight.com/meetings/mug01/Sayle/m4xbondage.html
// Some from molstar (github.com/molstar)
static const char* amino_acids[] = {
    "ACE", "ALA", "ARG", "ASC", "ASN", "ASP", "ASX",
    "CYS", "CYX",
    "FOR",
    "GLN", "GLU", "GLX", "GLY",
    "HIS", "HYP", "HISE",
    "ILE", "LEU",
    "LYS",
    "MET",
    "PHE", "PRO", "PCA", "PYL",
    "SER",
    "THR", "TRP", "TYR",
    "VAL", "SEC",
    "XLE",
    
    // CCD
    "UNK",
    "MSE", "SEP", "TPO", "PTR",

    // Charmm
    "HSD", "HSP", "LSN", "ASPP", "GLUP",

    // Amber
    "HID", "HIE", "HIP", "LYN", "ASH", "GLH",
};

static const char* nucleic_acids[] = {
    "A", "C", "G", "T", "U", "+U", "YG", "1MA", "1MG", "2MG", "5MC", "5MU", "7MG", "H2U", "M2G", "OMC", "OMG", "PSU",
};

/*static const char* peptides[] = {"APN", "CPN", "TPN", "GPN"};*/
static const char* rna[] = {"A", "C", "T", "G", "I", "U", "N"};
static const char* dna[] = {"DA", "DC", "DT", "DG", "DI", "DU", "DN"};

static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4", "W", "DOD", "D30" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

// Taken from here
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/model/properties/atomic/bonds.ts
/*
* Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
* @author Alexander Rose <alexander.rose@weirdbyte.de>
*/

// If an entry exists in this table, it has order 2
static str_t intra_bond_order_table[] = {
    BAKE("HIS|CD2|CG"),
    BAKE("HIS|CE1|ND1"),
    BAKE("ARG|CZ|NH2"),
    BAKE("PHE|CE1|CZ"),
    BAKE("PHE|CD2|CE2"),
    BAKE("PHE|CD1|CG"),
    BAKE("TRP|CD1|CG"),
    BAKE("TRP|CD2|CE2"),
    BAKE("TRP|CE3|CZ3"),
    BAKE("TRP|CH2|CZ2"),
    BAKE("ASN|CG|OD1"),
    BAKE("GLN|CD|OE1"),
    BAKE("TYR|CD1|CG"),
    BAKE("TYR|CD2|CE2"),
    BAKE("TYR|CE1|CZ"),
    BAKE("ASP|CG|OD1"),
    BAKE("GLU|CD|OE1"),

    BAKE("G|C8|N7"),
    BAKE("G|C4|C5"),
    BAKE("G|C2|N3"),
    BAKE("G|C6|O6"),
    BAKE("C|C4|N3"),
    BAKE("C|C5|C6"),
    BAKE("C|C2|O2"),
    BAKE("A|C2|N3"),
    BAKE("A|C6|N1"),
    BAKE("A|C4|C5"),
    BAKE("A|C8|N7"),
    BAKE("U|C5|C6"),
    BAKE("U|C2|O2"),
    BAKE("U|C4|O4"),

    BAKE("DG|C8|N7"),
    BAKE("DG|C4|C5"),
    BAKE("DG|C2|N3"),
    BAKE("DG|C6|O6"),
    BAKE("DC|C4|N3"),
    BAKE("DC|C5|C6"),
    BAKE("DC|C2|O2"),
    BAKE("DA|C2|N3"),
    BAKE("DA|C6|N1"),
    BAKE("DA|C4|C5"),
    BAKE("DA|C8|N7"),
    BAKE("DT|C5|C6"),
    BAKE("DT|C2|O2"),
    BAKE("DT|C4|O4"),
};

#include <stdlib.h>

static int compare_int(void const* a, void const* b) {
    return ( *(const int*)a - *(const int*)b );
}

// Simplistic inplace bubble sort for small arrays
// In practice, this is only used to sort the small rings
// Fallback is qsort for data larger than N
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

static void sort_arr_masked(int* arr, int n, int mask) {
    bool swapped = true;
    while (swapped) {
        swapped = false;
        for (int i = 0; i < n - 1; ++i) {
            if ((arr[i] & mask) > (arr[i + 1] & mask)) {
                int tmp = arr[i];
                arr[i] = arr[i + 1];
                arr[i + 1] = tmp;
                swapped = true;
            }
        }
    }
}

static void sort_arr_masked64(int64_t* arr, size_t n, int64_t mask) {
    bool swapped = true;
    while (swapped) {
        swapped = false;
        for (size_t i = 0; i < n - 1; ++i) {
            if ((arr[i] & mask) > (arr[i + 1] & mask)) {
                int64_t tmp = arr[i];
                arr[i] = arr[i + 1];
                arr[i + 1] = tmp;
                swapped = true;
            }
        }
    }
}

static inline void radix_pass_8(uint32_t* dst, const uint32_t* src, size_t count, uint32_t hist[256][4], int pass) {
    int bitoff = pass * 8;

    for (size_t i = 0; i < count; ++i) {
        uint32_t id = (src[i] >> bitoff) & 255;
        dst[hist[id][pass]++] = src[i];
    }
}

static void sort_radix_inplace_uint32(uint32_t* data, size_t count, md_allocator_i* temp_arena) {
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_arena);
    uint32_t* scratch = md_vm_arena_push_array(temp_arena,uint32_t, count);

    typedef union {
        uint32_t u32[256][4];
        md_128i  vec[256];
    } hist_t;

    hist_t hist = {0};
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = data[i];
        hist.u32[(id >> 0)  & 255][0]++;
        hist.u32[(id >> 8)  & 255][1]++;
        hist.u32[(id >> 16) & 255][2]++;
        hist.u32[(id >> 24) & 255][3]++;
    }

    // Prefix sum
    //uint32_t sum[4] = {0};
    md_128i sum = md_mm_setzero_si128();
    for (size_t i = 0; i < ARRAY_SIZE(hist.u32); ++i) {
        md_128i val = hist.vec[i];
        hist.vec[i] = sum;
        sum = md_mm_add_epi32(sum, val);
    }

    // 4-pass 8-bit radix sort computes the resulting order into scratch
    radix_pass_8(scratch, data, count, hist.u32, 0);
    radix_pass_8(data, scratch, count, hist.u32, 1);
    radix_pass_8(scratch, data, count, hist.u32, 2);
    radix_pass_8(data, scratch, count, hist.u32, 3);

    md_vm_arena_temp_end(temp);
}


static inline void radix_pass_indices_8(uint32_t* dst_indices, const uint32_t* keys, const uint32_t* src_indices, size_t count, uint32_t hist[256],
                                        int pass) {
    int bitoff = pass * 8;
    uint32_t pos[256];
    memcpy(pos, hist, sizeof(uint32_t) * 256);

    for (size_t i = 0; i < count; ++i) {
        uint32_t key = (keys[src_indices[i]] >> bitoff) & 0xFF;
        dst_indices[pos[key]++] = src_indices[i];
    }
}

void radix_sort_indices_uint32(const uint32_t* keys,    // array of keys (e.g., min_x of AABBs)
                               uint32_t* indices,       // output: sorted indices
                               uint32_t* temp_indices,  // scratch array of size 'count'
                               size_t count) {
    uint32_t hist[256];

    // Initialize indices
    for (size_t i = 0; i < count; ++i) indices[i] = (uint32_t)i;

    uint32_t* src = indices;
    uint32_t* dst = temp_indices;

    // 4 passes of 8-bit radix
    for (int pass = 0; pass < 4; ++pass) {
        MEMSET(hist, 0, sizeof(hist));

        // Compute histogram for this pass
        for (size_t i = 0; i < count; ++i) {
            uint32_t key = (keys[src[i]] >> (pass * 8)) & 0xFF;
            hist[key]++;
        }

        // Prefix sum
        uint32_t sum = 0;
        for (int i = 0; i < 256; ++i) {
            uint32_t tmp = hist[i];
            hist[i] = sum;
            sum += tmp;
        }

        // Sort pass
        radix_pass_indices_8(dst, keys, src, count, hist, pass);

        // Swap src and dst for next pass
        uint32_t* tmp = src;
        src = dst;
        dst = tmp;
    }

    // After even number of passes, result is in 'indices' already
    // If odd, copy back from temp
    if (src != indices) {
        MEMCPY(indices, src, count * sizeof(uint32_t));
    }
}

typedef struct fifo_t {
    int* data;
    uint32_t head;
    uint32_t tail;
    uint32_t size;
    uint32_t capacity;
    md_allocator_i* alloc;
} fifo_t;

static inline bool fifo_empty(fifo_t* fifo) { return fifo->size == 0; }
static inline bool fifo_full(fifo_t* fifo)  { return fifo->size == fifo->capacity; }

static inline fifo_t fifo_create(size_t capacity, md_allocator_i* alloc) {
    ASSERT(alloc);
    fifo_t fifo = {
        .data = md_alloc(alloc, sizeof(int) * next_power_of_two64(MAX(16, (uint64_t)capacity))),
        .head = 0,
        .tail = 0,
        .size = 0,
        .capacity = (uint32_t)capacity,
        .alloc = alloc,
    };
#if DEBUG
    // Clear memory to make debugging easier
    MEMSET(fifo.data, 0, sizeof(int) * fifo.capacity);
#endif
    return fifo;
}

static inline void fifo_free(fifo_t* fifo) {
    ASSERT(fifo);
    md_free(fifo->alloc, fifo->data, sizeof(int) * fifo->capacity);
    MEMSET(fifo, 0, sizeof(fifo_t));
}

static inline void fifo_clear(fifo_t* fifo) {
    fifo->head = 0;
    fifo->tail = 0;
    fifo->size = 0;
#if DEBUG
    // Clear memory to make debugging easier
    MEMSET(fifo->data, 0, sizeof(int) * fifo->capacity);
#endif
}

static inline void fifo_push(fifo_t* fifo, int value) {
    ASSERT(fifo);
    ASSERT(fifo->data);
    if (fifo_full(fifo)) {
        uint32_t new_capacity = next_power_of_two32(fifo->capacity * 2);
        fifo->data = md_realloc(fifo->alloc, fifo->data, sizeof(int) * fifo->capacity, sizeof(int) * new_capacity);
        fifo->capacity = new_capacity;
    }
    fifo->data[fifo->head] = value;
    fifo->head = (fifo->head + 1) & (fifo->capacity - 1);
    fifo->size += 1;
}

static inline int fifo_pop(fifo_t* fifo) {
    ASSERT(!fifo_empty(fifo));
    int val = fifo->data[fifo->tail];
    fifo->tail = (fifo->tail + 1) & (fifo->capacity - 1);
    fifo->size -= 1;
    return val;
}

static md_array(uint64_t) make_bitfield(size_t num_bits, md_allocator_i* alloc) {
    size_t num_elem = DIV_UP(num_bits, 64);
    md_array(uint64_t) bits = md_array_create(uint64_t, num_elem, alloc);
    MEMSET(bits, 0, md_array_bytes(bits));
    return bits;
}

static inline md_array(uint64_t) bitfield_copy(const md_array(uint64_t) src, md_allocator_i* alloc) {
    md_array(uint64_t) copy = md_array_create(uint64_t, md_array_size(src), alloc);
    MEMCPY(copy, src, md_array_bytes(src));
    return copy;
}

static inline void bitfield_copy_inplace(md_array(uint64_t) dst, const md_array(uint64_t) src) {
    MEMCPY(dst, src, md_array_bytes(src));
}

static inline void bitfield_clear_all(uint64_t* bits, size_t num_bits) {
    size_t num_bytes = DIV_UP(num_bits, 8);
    MEMSET(bits, 0, num_bytes);
}

static void bitfield_set_all(uint64_t* bits, size_t num_bits) {
    size_t num_bytes = DIV_UP(num_bits, 8);
    MEMSET(bits, 0xFF, num_bytes);
}

static bool bitfield_find_first_bit_set(int* out_idx, const uint64_t* bits, size_t num_bits) {
    ASSERT(out_idx);
    size_t num_elem = DIV_UP(num_bits, 64);
    for (size_t i = 0; i < num_elem; ++i) {
        if (bits[i]) {
            *out_idx = (int)ctz64(bits[i]);
            return true;
        }
    }
    return false;
}

static inline bool bitfield_test_bit(const uint64_t* bits, int64_t idx) {
    return (bits[idx >> 6] & (1ULL << (idx & 63)));
}

static inline void bitfield_set_bit(uint64_t* bits, int64_t idx) {
    bits[idx >> 6] |= (1ULL << (idx & 63));
}

static inline void bitfield_clear_bit(uint64_t* bits, int64_t idx) {
    bits[idx >> 6] &= ~(1ULL << (idx & 63));
}

static inline size_t bitfield_popcount(const uint64_t* bits, size_t num_bits) {
    ASSERT(bits);
    ASSERT(num_bits > 0);
    const uint64_t* it = bits;
    const uint64_t* end = bits + DIV_UP(num_bits, 64) - 1;
    size_t count = 0;
    while (it != end) {
        count += (size_t)popcnt64(*it);
        it++;
    }
    const uint64_t mask = (1ULL << (num_bits & 63)) - 1;
    count += (size_t)popcnt64(*end & mask);
    return count;
}

typedef struct graph_t {
    size_t    vertex_count;
    uint8_t*  vertex_type;
    uint32_t* edge_offset;  // offset, length is implicitly encoded by the next offset, last offset is the total number of edges and therefore length is count + 1
    uint32_t* edge_data;    // packed 32-bit data consiting of (from hi to low) type : 8, index : 24      

    // These map the graphs internal indices to the provided indices from which the graph was constructed
    md_atom_idx_t* atom_idx_map;
    md_bond_idx_t* bond_idx_map;
} graph_t;

static inline size_t graph_vertex_count(const graph_t* g) {
    return g->vertex_count;
}

static inline int graph_vertex_type(const graph_t* g, int64_t vidx) {
    return g->vertex_type[vidx];
}

static inline size_t graph_vertex_edge_count(const graph_t* g, int64_t vidx) {
    return g->edge_offset[vidx + 1] - g->edge_offset[vidx];
}

static inline int graph_edge_type(const graph_t* g, int64_t eidx) {
    return ((g->edge_data[eidx]) >> 24) & 0xFF;
}

static inline int graph_edge_vertex_idx(const graph_t* g, int64_t eidx) {
    return g->edge_data[eidx] & 0x00FFFFFF;
}

typedef struct graph_edge_iter_t {
    uint32_t* cur;
    uint32_t* end;
} graph_edge_iter_t;

static inline graph_edge_iter_t graph_edge_iter(const graph_t* g, int64_t vidx) {
    graph_edge_iter_t it = {
        .cur = g->edge_data + g->edge_offset[vidx],
        .end = g->edge_data + g->edge_offset[vidx + 1]
    };
    return it;
}

static inline bool graph_edge_iter_has_next(graph_edge_iter_t it) {
    return it.cur != it.end;
}

static inline bool graph_edge_iter_valid(graph_edge_iter_t it) {
    return it.cur != 0;
}

static inline void graph_edge_iter_next(graph_edge_iter_t* it) {
    ++it->cur;
}

static inline int graph_edge_iter_type(graph_edge_iter_t it) {
    return ((*it.cur) >> 24) & 0xFF;
}

static inline int graph_edge_iter_vidx(graph_edge_iter_t it) {
    return (*it.cur) & 0x00FFFFFF;
}

static inline bool graph_vertex_is_connected_to(const graph_t* g, int vidx, int other_vidx) {
    graph_edge_iter_t it = graph_edge_iter(g, vidx);
    while (graph_edge_iter_has_next(it)) {
        if (graph_edge_iter_vidx(it) == other_vidx) return true;
        graph_edge_iter_next(&it);
    }
    return false;
}

static inline bool graph_vertex_has_connection(const graph_t* g, int vidx, int other_vidx, int other_type) {
    graph_edge_iter_t it = graph_edge_iter(g, vidx);
    while (graph_edge_iter_has_next(it)) {
        int evidx = graph_edge_iter_vidx(it);
        int etype = graph_edge_iter_type(it);
        if (evidx == other_vidx && etype == other_type) {
            return true;
        }
        graph_edge_iter_next(&it);
    }
    return false;
}

static bool graph_equivalent(const graph_t* a, const graph_t* b) {
    if (a->vertex_count != b->vertex_count) return false;
    for (int64_t i = 0; i < (int64_t)a->vertex_count; ++i) {
        if (graph_vertex_type(a, i) != graph_vertex_type(b, i)) return false;
        if (graph_vertex_edge_count(a, i) != graph_vertex_edge_count(b, i)) return false;

        graph_edge_iter_t a_it = graph_edge_iter(a, i);
        graph_edge_iter_t b_it = graph_edge_iter(b, i);
        while (graph_edge_iter_has_next(a_it)) {
            bool found = false;
            while (graph_edge_iter_has_next(b_it)) {
                if (graph_edge_iter_vidx(a_it) == graph_edge_iter_vidx(b_it) &&
                    graph_edge_iter_type(a_it) == graph_edge_iter_type(b_it)) {
                    found = true;
                    break;
                }
                graph_edge_iter_next(&b_it);
            }
            if (!found) {
                return false;
            }
            graph_edge_iter_next(&a_it);
        }
    }
    return true;
}

// Confusing function
// Extracts a graph from an atom index range with supplied atom types

static graph_t make_graph(const md_bond_data_t* bond, const uint8_t atom_types[], const int indices[], size_t count, md_allocator_i* vm_arena) {   
    ASSERT(bond);

    // This is just an upper estimate of the number of edges that could potentially exist
    const size_t edge_data_cap = (count * 4);
    size_t edge_data_len = 0;

    graph_t graph = {
        .vertex_count = count,
        .vertex_type  = md_vm_arena_push     (vm_arena, sizeof(uint8_t)  * count),
        .edge_offset  = md_vm_arena_push_zero(vm_arena, sizeof(uint32_t) * (count + 1)),
        .edge_data    = md_vm_arena_push     (vm_arena, sizeof(uint32_t) * edge_data_cap),
        .atom_idx_map = md_vm_arena_push     (vm_arena, sizeof(uint32_t) * count),
        .bond_idx_map = md_vm_arena_push     (vm_arena, sizeof(uint32_t) * edge_data_cap),
    };

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(vm_arena);

    // Map from global indices (which the connectivity info is given in) to local (graph) indices
    md_hashmap32_t map = { .allocator = vm_arena };
    md_hashmap_reserve(&map, count);

    for (int i = 0; i < (int)count; ++i) {
        int idx = indices[i];
        md_hashmap_add(&map, (uint64_t)idx, i);
        graph.vertex_type[i] = atom_types[idx];
        graph.atom_idx_map[i] = idx;
    }

    // Only store edges which point to vertices within the graph as this will be used later as a traversal template
    for (size_t i = 0; i < count; ++i) {
        int idx = indices[i];
        // Translate the global atom indices to local structure indices
        uint64_t edge_data_arr[8];
        uint32_t length = 0;

        md_bond_iter_t it = md_bond_iter(bond, idx);
        while (md_bond_iter_has_next(it)) {
            uint32_t bond_idx = md_bond_iter_bond_index(it);
            uint32_t atom_idx = md_bond_iter_atom_index(it);
            uint32_t order    = md_bond_iter_bond_order(it);
            uint32_t* local_idx = md_hashmap_get(&map, atom_idx);
            if (local_idx) {
                // Only commit the edge if it is referring to a local index within the graph
                edge_data_arr[length++] = ((uint64_t)bond_idx << 32) | ((uint64_t)order << 24) | (uint32_t)(*local_idx);
            }
            md_bond_iter_next(&it);
        }

        if (length > 0) {
            // Sort on indices
            sort_arr_masked64((int64_t*)edge_data_arr, length, 0x00FFFFFF);

            // Note that this is not the true offset yet, only the number of local edges
            graph.edge_offset[i] = length;
            ASSERT(edge_data_len + length < edge_data_cap);
            for (uint32_t j = 0; j < length; ++j) {
                uint32_t edge_data = (uint32_t)(edge_data_arr[j]);
                uint32_t bond_idx  = (uint32_t)(edge_data_arr[j] >> 32);
                graph.edge_data   [edge_data_len] = edge_data;
                graph.bond_idx_map[edge_data_len] = bond_idx;
                edge_data_len += 1;
            }
        }
    }

    // Perform exclusive scan to convert local edge count to global edge offsets
    uint32_t offset = 0;
    for (size_t i = 0; i < count + 1; ++i) {
        uint32_t len = graph.edge_offset[i];
        graph.edge_offset[i] = offset;
        offset += len;
    }

    md_vm_arena_temp_end(temp);
    return graph;
}

typedef bool (*solution_callback)(const int map[], size_t length, void* user);

typedef struct state_t {
    bool abort;
    uint16_t flags;

    solution_callback callback;
    void* user_data;

    md_array(int) map;

    const graph_t* n_graph;
    const graph_t* h_graph;

    md_array(int) n_path;
    md_array(int) h_path;

    // Terminal sets
    md_array(uint64_t) n_path_bits;
    md_array(uint64_t) h_path_bits;
    md_array(uint32_t) n_depths;
    md_array(uint32_t) h_depths;
} state_t;

static void state_reset(state_t* state) {
    state->abort = false;
    md_array_shrink(state->n_path, 0);
    md_array_shrink(state->h_path, 0);
    MEMSET(state->map, -1, md_array_bytes(state->map));
    bitfield_clear_all(state->n_path_bits, state->n_graph->vertex_count);
    bitfield_clear_all(state->h_path_bits, state->h_graph->vertex_count);
    MEMSET(state->n_depths, 0, md_array_bytes(state->n_depths));
    MEMSET(state->h_depths, 0, md_array_bytes(state->h_depths));
}

static void state_init(state_t* state, const graph_t* n_graph, const graph_t* h_graph, md_allocator_i* alloc) {
    state->n_graph = n_graph;
    state->h_graph = h_graph;
    state->map = md_array_create(int, n_graph->vertex_count, alloc);
    state->n_path_bits = make_bitfield(n_graph->vertex_count, alloc);
    state->h_path_bits = make_bitfield(h_graph->vertex_count, alloc);
    state->n_depths = md_array_create(uint32_t, n_graph->vertex_count, alloc);
    state->h_depths = md_array_create(uint32_t, h_graph->vertex_count, alloc);
    md_array_ensure(state->n_path, n_graph->vertex_count, alloc);
    md_array_ensure(state->h_path, h_graph->vertex_count, alloc);

    state_reset(state);
}

static void state_free(state_t* state, md_allocator_i* alloc) {
    md_array_free(state->map, alloc);
    md_array_free(state->n_path, alloc);
    md_array_free(state->h_path, alloc);
    md_array_free(state->n_path_bits, alloc);
    md_array_free(state->h_path_bits, alloc);
    md_array_free(state->n_depths, alloc);
    md_array_free(state->h_depths, alloc);
}

static inline bool find_str_in_str_arr(size_t* out_loc, str_t str, const str_t str_arr[], size_t arr_len) {
    for (size_t i = 0; i < arr_len; ++i) {
        if (str_eq(str, str_arr[i])) {
            if (out_loc) *out_loc = i;
            return true;
        }
    }
    return false;
}

static inline bool find_str_in_cstr_arr(size_t* loc, str_t str, const char* arr[], size_t arr_len) {
    for (size_t i = 0; i < arr_len; ++i) {
        if (str_eq_cstr(str, arr[i])) {
            if (loc) *loc = i;
            return true;
        }
    }
    return false;
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

// Modeled after gromacs implementation which assumes a separation of less than half smallest box dimension (diagonal element)
// https://github.com/gromacs/gromacs/blob/4f423711fabc081226d43beef1ab80c33c25edf8/src/gromacs/mdlib/ns.cpp#L1263
// Mostly used as a reference for correctness
static float min_image_tric(float dx[3], const float box[3][3], const float inv_box_ext[3]) {
    int   t[3];

     t[2]  = (int)(dx[2]*inv_box_ext[2] + 2.5f) - 2;
    dx[2] -= t[2]*box[2][2];
    dx[1] -= t[2]*box[2][1];
    dx[0] -= t[2]*box[2][0];

    t[1]   = (int)(dx[1]*inv_box_ext[1] + 2.5f) - 2;
    dx[1] -= t[1]*box[1][1];
    dx[0] -= t[1]*box[1][0];

    t[0]   = (int)(dx[0]*inv_box_ext[0] + 2.5f) - 2;
    dx[0] -= t[0]*box[0][0];

    float r2 = (dx[0] * dx[0]) + (dx[1] * dx[1]) + (dx[2] * dx[2]);
    return r2;
}

static md_256 simd_min_image_tric(md_256 dx, md_256 dy, md_256 dz, const float box[3][3], const float inv_box_ext[3]) {
    md_256 tz = md_mm256_round_ps(md_mm256_mul_ps(dz, md_mm256_set1_ps(inv_box_ext[2])));
    dz = md_mm256_sub_ps(dz, md_mm256_mul_ps(tz, md_mm256_set1_ps(box[2][2])));
    dy = md_mm256_sub_ps(dy, md_mm256_mul_ps(tz, md_mm256_set1_ps(box[2][1])));
    dx = md_mm256_sub_ps(dx, md_mm256_mul_ps(tz, md_mm256_set1_ps(box[2][0])));

    md_256 ty = md_mm256_round_ps(md_mm256_mul_ps(dy, md_mm256_set1_ps(inv_box_ext[1])));
    dy = md_mm256_sub_ps(dy, md_mm256_mul_ps(ty, md_mm256_set1_ps(box[1][1])));
    dx = md_mm256_sub_ps(dx, md_mm256_mul_ps(ty, md_mm256_set1_ps(box[1][0])));

    md_256 tx = md_mm256_round_ps(md_mm256_mul_ps(dx, md_mm256_set1_ps(inv_box_ext[0])));
    dx = md_mm256_sub_ps(dx, md_mm256_mul_ps(tx, md_mm256_set1_ps(box[0][0])));

    md_256 r2 = md_mm256_fmadd_ps(dz, dz, md_mm256_fmadd_ps(dx, dx, md_mm256_mul_ps(dy, dy)));
    return r2;
}

static inline md_256 simd_deperiodize_ortho(md_256 x, md_256 r, md_256 period, md_256 inv_period) {
    md_256 d = md_mm256_sub_ps(x, r);
    md_256 dx = md_mm256_mul_ps(d, inv_period);
    dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
    md_256 x_prim = md_mm256_add_ps(r, md_mm256_mul_ps(dx, period));
    return md_mm256_blendv_ps(x_prim, x, md_mm256_cmpeq_ps(period, md_mm256_setzero_ps()));
}

static inline int simd_xyz_mask(float ext[3]) {
    int mask = 0;
    if (ext[0] > 0.0f) mask |= 1;
    if (ext[1] > 0.0f) mask |= 2;
    if (ext[2] > 0.0f) mask |= 4;
    return mask;
}

bool md_util_resname_rna(str_t str) {
    str = trim_label(str);
    return find_str_in_cstr_arr(NULL, str, rna, ARRAY_SIZE(rna));
}
    
bool md_util_resname_dna(str_t str) {
    str = trim_label(str);
    return find_str_in_cstr_arr(NULL, str, dna, ARRAY_SIZE(dna));
}

bool md_util_resname_nucleic_acid(str_t str) {
    str = trim_label(str);
    return find_str_in_cstr_arr(NULL, str, rna, ARRAY_SIZE(rna)) || find_str_in_cstr_arr(NULL, str, dna, ARRAY_SIZE(dna));
}
    
bool md_util_resname_acidic(str_t str) {
    return find_str_in_cstr_arr(NULL, str, acidic, ARRAY_SIZE(acidic));
}

bool md_util_resname_basic(str_t str) {
    return find_str_in_cstr_arr(NULL, str, basic, ARRAY_SIZE(basic));
}
    
bool md_util_resname_neutral(str_t str) {
    return find_str_in_cstr_arr(NULL, str, neutral, ARRAY_SIZE(neutral));
}
    
bool md_util_resname_water(str_t str) {
    return find_str_in_cstr_arr(NULL, str, water, ARRAY_SIZE(water));
}
    
bool md_util_resname_hydrophobic(str_t str) {
    return find_str_in_cstr_arr(NULL, str, hydrophobic, ARRAY_SIZE(hydrophobic));
}
    
bool md_util_resname_amino_acid(str_t str) {
    return find_str_in_cstr_arr(NULL, str, amino_acids, ARRAY_SIZE(amino_acids));
}

static inline bool aromatic_ring_element(md_element_t elem) {
    switch (elem) {
    case B:
    case C:
    case N:
    case O:
    case Si:
    case P:
    case S:
    case Ge:
    case As:
    case Sn:
    case Sb:
    case Bi:
        return true;
    default:
        return false;
    }
}

static inline bool metal_element(md_element_t elem) {
    switch(elem) {
    case Li:
    case Be:
    case Mg:
    case Al:
    case Ar:
    case V:
    case Cr:
    case Mn:
    case Co:
    case Cu:
    case Zn:
    case Ga:
    case As:
    case Kr:
    case Rb:
    case Sr:
    case Y:
    case Mo:
    case Ag:
    case Cd:
    case In:
    case Sb:
    case Te:
    case Xe:
    case Cs:
    case Ba:
    case La:
    case Ce:
    case Sm:
    case Eu:
    case Gd:
    case Tb:
    case Ho:
    case Yb:
    case Lu:
    case W:
    case Re:
    case Os:
    case Ir:
    case Pt:
    case Au:
    case Hg:
    case Tl:
    case Pb:
    case U:
        return true;
    default:
        return false;
    }
}

static inline bool halogen_element(md_element_t elem) {
    switch (elem) {
    case F:
    case Cl:
    case Br:
    case I:
    case At:
        return true;
    default:
        return false;
    }
}

md_element_t md_util_element_lookup(str_t str) {
    if (str.len == 1 || str.len == 2) {
        for (md_element_t i = 0; i < ARRAY_SIZE(element_symbols); ++i) {
            str_t sym = element_symbols[i];
            if (str_eq(str, sym)) return i;
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

static bool water_heuristic(const md_element_t elements[], size_t size) {
    if (size == 3) {
        size_t h_count = 0;
        size_t o_count = 0;
        for (size_t i = 0; i < 3; ++i) {
            h_count += (elements[i] == H) ? 1 : 0;
            o_count += (elements[i] == O) ? 1 : 0;
        }
        return o_count == 1 && h_count == 2;
    }
    return false;
}

static inline bool cmp1(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == '\0';
}

static inline bool cmp2(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == ref[1] && str[2] == '\0';
}

static inline bool cmp3(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == ref[1] && str[2] == ref[2] && str[3] == '\0';
}

static bool md_util_protein_backbone_atoms_extract(md_protein_backbone_atoms_t* backbone_atoms, const md_label_t* atom_types, size_t count, int32_t atom_offset) {
    if (count == 0) return false;

    const uint32_t all_bits = 1 | 2 | 4 | 8 | 16;
    uint32_t bits = 0;
    md_protein_backbone_atoms_t bb = {-1, -1, -1, -1, -1};
    for (int i = 0; i < (int)count; ++i) {
        if (!(bits & 1)  && cmp1(atom_types[i].buf, "N"))  { bb.n  = atom_offset + i; bits |= 1;  continue; }
        if (!(bits & 16) && cmp2(atom_types[i].buf, "HN")) { bb.hn = atom_offset + i; bits |= 16; continue; }
        if (!(bits & 2)  && cmp2(atom_types[i].buf, "CA")) { bb.ca = atom_offset + i; bits |= 2;  continue; }
        if (!(bits & 4)  && cmp1(atom_types[i].buf, "C"))  { bb.c  = atom_offset + i; bits |= 4;  continue; }
        if (!(bits & 8)  && cmp1(atom_types[i].buf, "O"))  { bb.o  = atom_offset + i; bits |= 8;  continue; }

        // Check if done
        if (bits == all_bits) break;
    }
    // Explicitly check and assign last atom if Oxgen
    if (!(bits & 8) && atom_types[count - 1].buf[0] == 'O') {
        bb.o = atom_offset + (int)count - 1;
        bits |= 8;
    }

    // HN is not
    const uint32_t req_bits = 1 | 2 | 4 | 8;
    if ((bits & req_bits) == req_bits) {
        if (backbone_atoms) *backbone_atoms = bb;
        return true;
    }
    return false;
}

static bool md_util_nucleic_backbone_atoms_extract(md_nucleic_backbone_atoms_t* backbone_atoms, const md_label_t* atom_types, size_t count, int32_t atom_offset) {
    if (count == 0) return false;

    const uint32_t all_bits = 1 | 2 | 4 | 8 | 16 | 32;
    uint32_t bits = 0;
    md_nucleic_backbone_atoms_t bb = {0};
    for (int i = 0; i < (int)count; ++i) {
        if (!(bits & 1)  && cmp3(atom_types[i].buf, "C5'")) { bb.c5 = atom_offset + i; bits |= 1;  continue; }
        if (!(bits & 2)  && cmp3(atom_types[i].buf, "C4'")) { bb.c4 = atom_offset + i; bits |= 2;  continue; }
        if (!(bits & 4)  && cmp3(atom_types[i].buf, "C3'")) { bb.c3 = atom_offset + i; bits |= 4;  continue; }
        if (!(bits & 8)  && cmp3(atom_types[i].buf, "O3'")) { bb.o3 = atom_offset + i; bits |= 8;  continue; }
        if (!(bits & 16) && cmp1(atom_types[i].buf, "P"))   { bb.p  = atom_offset + i; bits |= 16; continue; }
        if (!(bits & 32) && cmp3(atom_types[i].buf, "O5'")) { bb.o5 = atom_offset + i; bits |= 32; continue; }

        // Check if done
        if (bits == all_bits) break;
    }

    if (bits == all_bits) {
        if (backbone_atoms) *backbone_atoms = bb;
        return true;
    }
    return false;
}

static inline bool is_organic(char c) {
    return c == 'C' || c == 'N' || c == 'O' || c == 'S' || c == 'P';
}

bool md_util_element_guess(md_element_t element[], size_t capacity, const struct md_molecule_t* mol) {
    ASSERT(capacity > 0);
    ASSERT(mol);
    ASSERT(mol->atom.count > 0);

    md_hashmap32_t map = { .allocator = md_get_temp_allocator() };
    md_hashmap_reserve(&map, 256);

    // Just for pure elements which have not been salted with resname
    md_hashmap32_t elem_map = { .allocator = md_get_temp_allocator() };
    md_hashmap_reserve(&elem_map, 256);

    typedef struct {
        str_t name;
        md_element_t elem;
    } entry_t;
    
    // Extra table for predefined atom types
    entry_t entries[] = {
        {STR_LIT("SOD"), Na},
        {STR_LIT("OW"),  O},
        {STR_LIT("HW"),  H},
    };

    for (size_t i = 0; i < ARRAY_SIZE(entries); ++i) {
        md_hashmap_add(&elem_map, md_hash64(entries[i].name.ptr, entries[i].name.len, 0), entries[i].elem);
    }

    const size_t count = MIN(capacity, mol->atom.count);
    for (size_t i = 0; i < count; ++i) {
        if (element[i] != 0) continue;

        str_t original = LBL_TO_STR(mol->atom.type[i]);

        // Trim whitespace, digits and 'X's
        str_t name = trim_label(original);

        if (name.len > 0) {
            md_element_t elem = 0;

            str_t resname = STR_LIT("");
            uint64_t res_key = 0;
            if (mol->atom.resname) {
                resname = LBL_TO_STR(mol->atom.resname[i]);
                res_key = md_hash64(resname.ptr, resname.len, 0);
            }
            uint64_t key = md_hash64(name.ptr, name.len, res_key);
            uint32_t* ptr = md_hashmap_get(&map, key);
            if (ptr) {
                element[i] = (md_element_t)*ptr;
                continue;
            } else {
                uint64_t elem_key = md_hash64(name.ptr, name.len, 0);
                ptr = md_hashmap_get(&elem_map, elem_key);
                if (ptr) {
                    elem = (md_element_t)*ptr;
                    goto done;
                }
            }

            if ((elem = md_util_element_lookup(name)) != 0) goto done;

            // If amino acid, try to deduce the element from that
            if (mol->atom.flags) {
                if (mol->atom.flags[i] & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE)) {
                    // Try to match against the first character
                    name.len = 1;
                    elem = md_util_element_lookup_ignore_case(name);
                    goto done;
                }
            }

            // This is the same logic as above but more general, for the natural organic elements
            if (is_organic(name.ptr[0]) && name.len > 1) {
                if (name.ptr[1] - 'A' < 5) {
                    if (mol->residue.count > 0 && mol->atom.res_idx) {
                        int32_t res_idx = mol->atom.res_idx[i];
                        uint32_t res_beg = mol->residue.atom_offset[res_idx];
                        uint32_t res_end = mol->residue.atom_offset[res_idx+1];
                        uint32_t res_len = res_end - res_beg;
                        if (res_len > 3) {
                            name.len = 1;
                            elem = md_util_element_lookup_ignore_case(name);
                            goto done;
                        }
                    }
                }
            }

            // Heuristic cases

            // This can be fishy...
            if (str_eq_cstr(name, "HOH")) {
                elem = H;
                goto done;
            }
            if (str_eq_cstr(name, "HS")) {
                elem = H;
                goto done;
            }

            size_t num_alpha = 0;
            while (num_alpha < original.len && is_alpha(original.ptr[num_alpha])) ++num_alpha;
            
            size_t num_digits = 0;
            str_t digits = str_substr(original, num_alpha, SIZE_MAX);
            while (num_digits < digits.len && is_digit(digits.ptr[num_digits])) ++num_digits;

            // 2-3 letters + 1-2 digit (e.g. HO(H)[0-99]) usually means just look at the first letter
            if ((num_alpha == 2 || num_alpha == 3) && (num_digits == 1 || num_digits == 2)) {
                name.len = 1;
                elem = md_util_element_lookup_ignore_case(name);
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
            if (elem != 0) {
                md_hashmap_add(&map, key, elem);
            }
        }
    }

    md_hashmap_free(&map);

    return true;
}

bool md_util_element_from_mass(md_element_t element[], const float mass[], size_t count) {
    if (!element) {
        MD_LOG_ERROR("element is null");
        return false;
    }

    if (!mass) {
        MD_LOG_ERROR("mass is null");
        return false;
    }

    size_t failed_matches = 0;

    const float eps = 1.0e-3f;
    for (size_t i = 0; i < count; ++i) {
        md_element_t elem = 0;
        const float m = mass[i];

        if (0.0f < m && m != 1.0f) {
            // Linear search for matching atomic mass
            for (uint8_t j = 1; j < (uint8_t)ARRAY_SIZE(element_atomic_mass); ++j) {
                if (fabs(m - element_atomic_mass[j]) < eps) {
                    elem = j;
                    break;
                }
            }
        }

        element[i] = elem;
        if (elem == 0) {
            //Found no matching element for mass
            failed_matches++;
        }
    }
    //Returns true if all masses found a matching element
    //Elements with a value of 0 indicates no match
    if (failed_matches == 0) {
        return true;
    }
    else {
        MD_LOG_ERROR("%zu masses had no matching element", failed_matches);
        return false;
    }
}

const str_t* md_util_element_symbols(void) {
    return element_symbols;
}

const str_t* md_util_element_names(void) {
    return element_names;
}

const float* md_util_element_vdw_radii(void) {
    return element_vdw_radii;
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

static inline float element_covalent_radius(md_element_t element) {
    return element < Num_Elements ? element_covalent_radii_u8[element] * 0.01f : 0;
}

float md_util_element_covalent_radius(md_element_t element) {
    return element_covalent_radius(element);
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

// Blatantly stolen from MDAnalysis project
// MDAnalysis --- https://www.mdanalysis.org
// Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
// https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/include/calc_distances.h

static inline void minimum_image_triclinic(float dx[3], const float box[3][3]) {
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

static inline void simd_minimum_image_triclinic(md_256 dx[3], const float box[3][3]) {
    md_256 dx_min[3] = {md_mm256_set1_ps(0), md_mm256_set1_ps(0), md_mm256_set1_ps(0)};
    md_256 dsq_min   = md_mm256_set1_ps(FLT_MAX);
    md_256 dsq;
    md_256 rx;
    md_256 ry[2];
    md_256 rz[3];
    int ix, iy, iz;
    for (ix = -1; ix < 2; ++ix) {
        md_256 fix = md_mm256_set1_ps((float)ix);
        rx = md_mm256_add_ps(dx[0], md_mm256_mul_ps(md_mm256_set1_ps(box[0][0]), fix));
        for (iy = -1; iy < 2; ++iy) {
            md_256 fiy = md_mm256_set1_ps((float)iy);
            ry[0] = md_mm256_add_ps(rx,    md_mm256_mul_ps(md_mm256_set1_ps(box[1][0]), fiy));
            ry[1] = md_mm256_add_ps(dx[1], md_mm256_mul_ps(md_mm256_set1_ps(box[1][1]), fiy));
            for (iz = -1; iz < 2; ++iz) {
                md_256 fiz = md_mm256_set1_ps((float)iz);
                rz[0] = md_mm256_add_ps(ry[0], md_mm256_mul_ps(md_mm256_set1_ps(box[2][0]), fiz));
                rz[1] = md_mm256_add_ps(ry[1], md_mm256_mul_ps(md_mm256_set1_ps(box[2][1]), fiz));
                rz[2] = md_mm256_add_ps(dx[2], md_mm256_mul_ps(md_mm256_set1_ps(box[2][2]), fiz));
                dsq = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(rz[0], rz[0]), md_mm256_mul_ps(rz[1], rz[1])), md_mm256_mul_ps(rz[2], rz[2]));
                md_256 mask = md_mm256_cmplt_ps(dsq, dsq_min);
                dsq_min = md_mm256_blendv_ps(dsq_min, dsq, mask);
                dx_min[0] = md_mm256_blendv_ps(dx_min[0], rz[0], mask);
                dx_min[1] = md_mm256_blendv_ps(dx_min[1], rz[1], mask);
                dx_min[2] = md_mm256_blendv_ps(dx_min[2], rz[2], mask);
            }
        }
    }
    dx[0] = dx_min[0];
    dx[1] = dx_min[1];
    dx[2] = dx_min[2];
}

static inline void deperiodize_triclinic(float x[3], const float rx[3], const float box[3][3]) {
    float dx[3];

    dx[0] = x[0] - rx[0];
    dx[1] = x[1] - rx[1];
    dx[2] = x[2] - rx[2];

    minimum_image_triclinic(dx, box);

    x[0] = rx[0] + dx[0];
    x[1] = rx[1] + dx[1];
    x[2] = rx[2] + dx[2];
}

static inline void simd_deperiodize_triclinic(md_256 x[3], const md_256 rx[3], const float box[3][3]) {
    md_256 dx[3];

    dx[0] = md_mm256_sub_ps(x[0], rx[0]);
    dx[1] = md_mm256_sub_ps(x[1], rx[1]);
    dx[2] = md_mm256_sub_ps(x[2], rx[2]);

    simd_minimum_image_triclinic(dx, box);

    x[0] = md_mm256_add_ps(rx[0], dx[0]);
    x[1] = md_mm256_add_ps(rx[1], dx[1]);
    x[2] = md_mm256_add_ps(rx[2], dx[2]);
}

static inline bool zhang_skolnick_ss(const md_molecule_t* mol, md_range_t bb_range, int i, const float distances[3], float delta) {
    ASSERT(mol);
    ASSERT(mol->atom.x);
    ASSERT(mol->atom.y);
    ASSERT(mol->atom.z);
    ASSERT(mol->protein_backbone.atoms);
    for (int j = MAX((int)bb_range.beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= (int)bb_range.end) continue;
            const int ca_j = mol->protein_backbone.atoms[j].ca;
            const int ca_k = mol->protein_backbone.atoms[j + k].ca;
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

bool tm_align(md_secondary_structure_t secondary_structure[], size_t capacity, const md_molecule_t* mol) {
    if (!secondary_structure) return false;
    if (capacity == 0) return false;

    MEMSET(secondary_structure, 0, capacity * sizeof(md_secondary_structure_t));

    if (!mol) return false;
    if (!mol->atom.x) return false;
    if (!mol->atom.y) return false;
    if (!mol->atom.z) return false;
    if (!mol->protein_backbone.atoms) return false;
    if (!mol->protein_backbone.range.offset) return false;

    for (size_t bb_idx = 0; bb_idx < mol->protein_backbone.range.count; ++bb_idx) {
        const md_range_t range = {mol->protein_backbone.range.offset[bb_idx], mol->protein_backbone.range.offset[bb_idx + 1]};
        ASSERT(range.end <= (int)capacity);

        if (range.end - range.beg < 4) {
            continue;
        }

        // Classify residues
        for (int32_t i = range.beg; i < range.end; ++i) {
            md_secondary_structure_t ss = MD_SECONDARY_STRUCTURE_COIL;

            if (is_sheet(mol, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_BETA_SHEET;
            }
            else if (is_helical(mol, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
            }
            secondary_structure[i] = ss;
        }

        // Set squished isolated structures to the surrounding (only for sheets and helices)
        md_secondary_structure_t* ss = secondary_structure;
        for (int32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i-1] == ss[i+1] &&
                ss[i-1] != MD_SECONDARY_STRUCTURE_COIL &&
                ss[i+1] != MD_SECONDARY_STRUCTURE_COIL)
            {
                ss[i] = ss[i-1];
            }
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MD_SECONDARY_STRUCTURE_COIL;
        for (int32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i] != ss[i-1] &&
                ss[i] != ss[i+1])
            {
                ss[i] = MD_SECONDARY_STRUCTURE_COIL;
            }
        }
        if (ss[range.end - 1] != ss[range.end - 2]) ss[range.end - 1] = MD_SECONDARY_STRUCTURE_COIL;
    }

    return true;
}

typedef struct dssp_hbond_t {
    float energy;
    uint32_t res_idx;
} dssp_hbond_t;

typedef struct dssp_res_coords_t {
    vec4_t N;
    vec4_t O;
    vec4_t C;
    vec4_t H;
    vec4_t CA;
} dssp_res_coords_t;

typedef struct dssp_res_hbonds_t {
    dssp_hbond_t acc[2];
    dssp_hbond_t don[2];
} dssp_res_hbonds_t;

static inline float calc_hbond_energy(dssp_res_hbonds_t res_hbonds[], const dssp_res_coords_t res_coords[], size_t don_idx, size_t acc_idx) {
    /*
    float d_NO = vec3_distance(donor->N, acceptor->O);
    float d_HC = vec3_distance(donor->H, acceptor->C);
    float d_NC = vec3_distance(donor->N, acceptor->C);
    float d_OH = vec3_distance(donor->O, acceptor->H);
    double e_ref = 0.084 * (1.0 / d_NO + 1.0 / d_HC - 1.0 / d_NC - 1.0 / d_OH) * 332.0;
    */

    const dssp_res_coords_t* don = &res_coords[don_idx];
    const dssp_res_coords_t* acc = &res_coords[acc_idx];

    float result = 0.0f;

    const float min_distance_sq = 0.5f * 0.5f;
    const float min_bond_energy = -9.9f;

    const vec4_t d2 = {
        vec4_distance_squared(don->N, acc->O),
        vec4_distance_squared(don->H, acc->C),
        vec4_distance_squared(don->N, acc->C),
        vec4_distance_squared(don->O, acc->H),
    };
    
    const int mask_min = vec4_move_mask(vec4_cmp_lt(d2, vec4_set1(min_distance_sq)));
    if (mask_min > 0) {
        result = min_bond_energy;
    } else {
        const vec4_t v = vec4_rsqrt(d2);
        result = (0.084f * 332.0f) * vec4_reduce_add(vec4_mul(v, vec4_set(1.0f, 1.0f, -1.0f, -1.0f)));
    }

    result = MAX(result, min_bond_energy);

    dssp_res_hbonds_t* don_hbonds = &res_hbonds[don_idx];
    dssp_res_hbonds_t* acc_hbonds = &res_hbonds[acc_idx];

    if (result < don_hbonds->acc[0].energy) {
        don_hbonds->acc[1] = don_hbonds->acc[0];
        don_hbonds->acc[0].res_idx = (uint32_t)acc_idx;
        don_hbonds->acc[0].energy  = result;
    } else if (result < don_hbonds->acc[1].energy) {
        don_hbonds->acc[1].res_idx = (uint32_t)acc_idx;
        don_hbonds->acc[1].energy  = result;
    }

    if (result < acc_hbonds->don[0].energy) {
        acc_hbonds->don[1] = acc_hbonds->don[0];
        acc_hbonds->don[0].res_idx = (uint32_t)don_idx;
        acc_hbonds->don[0].energy  = result;
    } else if (result < acc_hbonds->don[1].energy) {
        acc_hbonds->don[1].res_idx = (uint32_t)don_idx;
        acc_hbonds->don[1].energy  = result;
    }

    return result;
}

static inline vec4_t estimate_HN(vec4_t N, vec4_t CA, vec4_t C_prev) {
    const float NH_dist = 1.01f;
    vec4_t U = vec4_normalize(vec4_sub(N, CA));
    vec4_t V = vec4_normalize(vec4_sub(N, C_prev));
    return vec4_add(N, vec4_mul_f(vec4_normalize(vec4_add(U, V)), NH_dist));
}

static inline vec4_t estimate_terminal_HN(vec4_t N, vec4_t CA, vec4_t C) {
    // i == range.beg (N-terminus): estimate previous C by extrapolation
    // C_prev_est = N + (CA - C)
    vec4_t CA_minus_C = vec4_sub(CA, C);
    vec4_t Cprev_est  = vec4_add(N, CA_minus_C);

    // Try to construct H using the standard estimator
    return estimate_HN(N, CA, Cprev_est);
}

static inline bool dssp_test_bond(const dssp_res_hbonds_t res_hbonds[], size_t res_a, size_t res_b) {
    const float max_hbond_energy = -0.5f;
    return (res_hbonds[res_a].acc[0].res_idx == res_b && res_hbonds[res_a].acc[0].energy < max_hbond_energy) ||
           (res_hbonds[res_a].acc[1].res_idx == res_b && res_hbonds[res_a].acc[1].energy < max_hbond_energy);
}

enum {
    SS_FLAG_NONE = 0,
    SS_FLAG_TURN = 1,
    SS_FLAG_BEND = 2,
    SS_FLAG_HELIX_310 = 4,
    SS_FLAG_HELIX_ALPHA = 8,
    SS_FLAG_HELIX_PI = 16,
    SS_FLAG_SHEET = 32,
    SS_FLAG_BRIDGE = 64,
};

void dssp(md_secondary_structure_t out_secondary_structure[], const float* x, const float* y, const float* z,
    const struct md_protein_backbone_atoms_t* backbone_atoms, size_t backbone_count, const uint32_t* backbone_range_offsets, size_t backbone_range_count) {
    ASSERT(out_secondary_structure);

    md_allocator_i* temp_alloc = md_vm_arena_create(GIGABYTES(1));

    dssp_res_coords_t* res_coords = md_vm_arena_push(temp_alloc, sizeof(dssp_res_coords_t) * backbone_count);
    dssp_res_hbonds_t* res_hbonds = md_vm_arena_push(temp_alloc, sizeof(dssp_res_hbonds_t) * backbone_count);
    uint32_t* ss_flags = md_vm_arena_push_zero(temp_alloc, sizeof(uint32_t) * backbone_count);

    for (size_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_range_t range = {backbone_range_offsets[range_idx], backbone_range_offsets[range_idx + 1]};

        for (int i = range.beg; i < range.end; ++i) {
            md_atom_idx_t ca_idx = backbone_atoms[i].ca;
            md_atom_idx_t n_idx  = backbone_atoms[i].n;
            md_atom_idx_t c_idx  = backbone_atoms[i].c;
            md_atom_idx_t o_idx  = backbone_atoms[i].o;
            md_atom_idx_t hn_idx = backbone_atoms[i].hn;

            res_coords[i].CA = vec4_set(x[ca_idx], y[ca_idx], z[ca_idx], 0);
            res_coords[i].N  = vec4_set(x[n_idx],  y[n_idx],  z[n_idx],  0);
            res_coords[i].C  = vec4_set(x[c_idx],  y[c_idx],  z[c_idx],  0);
            res_coords[i].O  = vec4_set(x[o_idx],  y[o_idx],  z[o_idx],  0);

            if (hn_idx > 0) {
                res_coords[i].H = vec4_set(x[hn_idx], y[hn_idx], z[hn_idx], 0);
            } else {
                if (i > 0)
                    res_coords[i].H = estimate_HN(res_coords[i].N, res_coords[i].CA, res_coords[i - 1].C);
                else
                    res_coords[i].H = estimate_terminal_HN(res_coords[i].N, res_coords[i].CA, res_coords[i].C);
            }
        }
    }

    // Do N^2 bond energy calculations
    MEMSET(res_hbonds, 0, sizeof(dssp_res_hbonds_t) * backbone_count);
    const float min_ca_dist2 = (9.0f * 9.0f);
    for (size_t i = 0; i < backbone_count - 1; ++i) {
        for (size_t j = i + 1; j < backbone_count; ++j) {
            if (vec4_distance_squared(res_coords[i].CA, res_coords[j].CA) > min_ca_dist2) {
                continue;
            }

            calc_hbond_energy(res_hbonds, res_coords, i, j);
            if (j != i + 1) {
                calc_hbond_energy(res_hbonds, res_coords, j, i);
            }
        }
    }

    // Classify secondary structure
    MEMSET(out_secondary_structure, MD_SECONDARY_STRUCTURE_COIL, sizeof(md_secondary_structure_t) * backbone_count);

    // ---- Per-range passes: helices, turns, bends (set flags only) ----
    for (size_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_range_t range = {backbone_range_offsets[range_idx], backbone_range_offsets[range_idx + 1]};

        // Helices (i -> i+3, i+4, i+5)
        for (size_t i = range.beg; i < range.end; ++i) {
            if (i + 3 < range.end && dssp_test_bond(res_hbonds, i, i + 3)) {
                ss_flags[i] |= SS_FLAG_HELIX_310;
            }
            if (i + 4 < range.end && dssp_test_bond(res_hbonds, i, i + 4)) {
                ss_flags[i] |= SS_FLAG_HELIX_ALPHA;
            }
            if (i + 5 < range.end && dssp_test_bond(res_hbonds, i, i + 5)) {
                ss_flags[i] |= SS_FLAG_HELIX_PI;
            }
        }

        // Turns: 4-residue H-bonded turn i->i+3 (only set if coil-ish)
        for (size_t i = range.beg; i + 3 < range.end; ++i) {
            if (dssp_test_bond(res_hbonds, i, i + 3)) {
                // mark the 4 residues as turn candidates (don't override helix flags)
                for (size_t k = 0; k < 4; ++k) {
                    ss_flags[i + k] |= SS_FLAG_TURN;
                }
            }
        }

        // Bends: geometric Cα pseudo-angle around i (i-2 .. i+2)
        for (size_t i = range.beg + 2; i + 2 < range.end; ++i) {
            vec4_t v1 = vec4_sub(res_coords[i - 2].CA, res_coords[i].CA);
            vec4_t v2 = vec4_sub(res_coords[i + 2].CA, res_coords[i].CA);
            float angle = vec4_angle(v1, v2);  // radians
            if (RAD_TO_DEG(angle) < 70.0f) {
                ss_flags[i] |= SS_FLAG_BEND;
            }
        }
    }

    // ---- Global pass: sheets / bridges (non-local, can cross ranges) ----
    const float min_ca_dist2_sheet = 5.5f * 5.5f;  // slightly relaxed CA cutoff
    for (size_t i = 0; i < backbone_count; ++i) {
        // We will track how many partner hbonds found and whether ladders exist.
        int partner_count = 0;

        for (size_t j = 0; j < backbone_count; ++j) {
            if (i == j) continue;
            if (abs((int)i - (int)j) <= 2) continue;  // skip sequence neighbors
            if (vec4_distance_squared(res_coords[i].CA, res_coords[j].CA) > min_ca_dist2_sheet) continue;

            bool ij = dssp_test_bond(res_hbonds, i, j);
            bool ji = dssp_test_bond(res_hbonds, j, i);
            if (!(ij || ji)) continue;

            // Check ladder-style partners to prefer sheets over isolated bridges:
            bool ladder = false;
            // antiparallel: i bonds j  AND (i-1 binds j+1 OR i+1 binds j-1)
            if (i > 0 && j + 1 < backbone_count) {
                if (dssp_test_bond(res_hbonds, i - 1, j + 1) || dssp_test_bond(res_hbonds, j + 1, i - 1)) ladder = true;
            }
            if (i + 1 < backbone_count && j > 0) {
                if (dssp_test_bond(res_hbonds, i + 1, j - 1) || dssp_test_bond(res_hbonds, j - 1, i + 1)) ladder = true;
            }
            // parallel ladder: i+1 binds j+1 or i-1 binds j-1 also OK
            if (i + 1 < backbone_count && j + 1 < backbone_count) {
                if (dssp_test_bond(res_hbonds, i + 1, j + 1) || dssp_test_bond(res_hbonds, j + 1, i + 1)) ladder = true;
            }
            if (i > 0 && j > 0) {
                if (dssp_test_bond(res_hbonds, i - 1, j - 1) || dssp_test_bond(res_hbonds, j - 1, i - 1)) ladder = true;
            }

            if (ladder) {
                ss_flags[i] |= SS_FLAG_SHEET;
                ss_flags[j] |= SS_FLAG_SHEET;
            } else {
                // isolated single H-bond — candidate bridge
                ss_flags[i] |= SS_FLAG_BRIDGE;
                ss_flags[j] |= SS_FLAG_BRIDGE;
            }
            partner_count++;
        }
    }

    // ---- Final resolution: convert flags to enum with priority ----
    // Priority: HELIX (alpha > 310 > pi) > SHEET > BRIDGE > TURN > BEND > COIL
    for (size_t i = 0; i < backbone_count; ++i) {
        uint8_t f = ss_flags[i];

        // prefer alpha helix, but allow other helix flags if alpha not present
        if (f & SS_FLAG_HELIX_ALPHA) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
        } else if (f & SS_FLAG_HELIX_310) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_310;
        } else if (f & SS_FLAG_HELIX_PI) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_PI;
        } else if (f & SS_FLAG_SHEET) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_BETA_SHEET;
        } else if (f & SS_FLAG_BRIDGE) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_BETA_BRIDGE;
        } else if (f & SS_FLAG_TURN) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_TURN;
        } else if (f & SS_FLAG_BEND) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_BEND;
        } else {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_COIL;
        }
    }

    md_vm_arena_destroy(temp_alloc);
}

bool md_util_backbone_secondary_structure_compute(md_secondary_structure_t secondary_structure[], size_t capacity, const struct md_molecule_t* mol) {
    return tm_align(secondary_structure, capacity, mol);
    
    //dssp(secondary_structure, mol->atom.x, mol->atom.y, mol->atom.z, mol->protein_backbone.atoms, mol->protein_backbone.count, mol->protein_backbone.range.offset, mol->protein_backbone.range.count);
}

bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const md_molecule_t* mol) {
    if (!backbone_angles) return false;
    if (capacity == 0) return false;

    MEMSET(backbone_angles, 0, sizeof(md_backbone_angles_t) * capacity);

    if (!mol) return false;
    if (!mol->atom.x) return false;
    if (!mol->atom.y) return false;
    if (!mol->atom.z) return false;
    if (!mol->protein_backbone.atoms) return false;

    for (size_t bb_idx = 0; bb_idx < mol->protein_backbone.range.count; ++bb_idx) {
        const md_range_t range = {mol->protein_backbone.range.offset[bb_idx], mol->protein_backbone.range.offset[bb_idx + 1]};
        ASSERT(range.end <= (int)capacity);

        if (range.end - range.beg < 4) {
            continue;
        }

        for (int64_t i = range.beg + 1; i < range.end - 1; ++i) {
            const vec3_t c_prev = { mol->atom.x[mol->protein_backbone.atoms[i-1].c], mol->atom.y[mol->protein_backbone.atoms[i-1].c], mol->atom.z[mol->protein_backbone.atoms[i-1].c] };
            const vec3_t n      = { mol->atom.x[mol->protein_backbone.atoms[i].n]  , mol->atom.y[mol->protein_backbone.atoms[i].n]  , mol->atom.z[mol->protein_backbone.atoms[i].n]   };
            const vec3_t ca     = { mol->atom.x[mol->protein_backbone.atoms[i].ca] , mol->atom.y[mol->protein_backbone.atoms[i].ca] , mol->atom.z[mol->protein_backbone.atoms[i].ca]  };
            const vec3_t c      = { mol->atom.x[mol->protein_backbone.atoms[i].c]  , mol->atom.y[mol->protein_backbone.atoms[i].c]  , mol->atom.z[mol->protein_backbone.atoms[i].c]   };
            const vec3_t n_next = { mol->atom.x[mol->protein_backbone.atoms[i+1].n], mol->atom.y[mol->protein_backbone.atoms[i+1].n], mol->atom.z[mol->protein_backbone.atoms[i+1].n] };

            vec3_t d[4] = {
                vec3_sub(n, c_prev),
                vec3_sub(ca, n),
                vec3_sub(c, ca),
                vec3_sub(n_next, c),
            };

            md_util_min_image_vec3(d, ARRAY_SIZE(d), &mol->unit_cell);

            backbone_angles[i].phi = vec3_dihedral_angle(d[0], d[1], d[2]);
            backbone_angles[i].psi = vec3_dihedral_angle(d[1], d[2], d[3]);
        }
    }

    return true;
}

bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const md_molecule_t* mol) {
    ASSERT(ramachandran_types);
    MEMSET(ramachandran_types, MD_RAMACHANDRAN_TYPE_UNKNOWN, sizeof(md_ramachandran_type_t) * capacity);

    if (capacity == 0) return false;
    if (mol->protein_backbone.count == 0) return false;
    if (mol->residue.count == 0) return false;

    ASSERT(mol->residue.name);
    ASSERT(mol->protein_backbone.residue_idx);

    size_t size = MIN(capacity, mol->protein_backbone.count);

    for (size_t i = 0; i < size; ++i) {
        size_t res_idx = mol->protein_backbone.residue_idx[i];
        ASSERT(res_idx < mol->residue.count);

        str_t resname = LBL_TO_STR(mol->residue.name[res_idx]);
        if (str_eq_cstr_n(resname, "GLY", 3)) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GLYCINE;
        } else if (str_eq_cstr_n(resname, "PRO", 3)) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PROLINE;
            ramachandran_types[i - 1] = MD_RAMACHANDRAN_TYPE_PREPROL;
        } else {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GENERAL;
        }
    }

    return true;
}

static void compute_connectivity(md_conn_data_t* out_conn, const md_atom_pair_t* bond_pairs, size_t bond_count, size_t atom_count, md_allocator_i* alloc, md_allocator_i* temp_arena) {    
    ASSERT(out_conn);
    ASSERT(alloc);

    if (bond_pairs == NULL) return;
    if (bond_count == 0) return;
    if (atom_count == 0) return;

    out_conn->offset_count = atom_count + 1;
    out_conn->offset = md_array_create(uint32_t, out_conn->offset_count, alloc);
    MEMSET(out_conn->offset, 0, md_array_bytes(out_conn->offset));

    // This have length of 2 * bond_count (one for each direction of the bond)
    out_conn->count = 2 * bond_count;
    out_conn->atom_idx = md_array_create(md_atom_idx_t, out_conn->count, alloc);
    out_conn->bond_idx = md_array_create(md_bond_idx_t, out_conn->count, alloc);

    typedef struct {
        uint16_t off[2];
    } offset_t;

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_arena);
    offset_t* local_offset = md_vm_arena_push_zero_array(temp_arena, offset_t, bond_count);

    // Two packed 16-bit local offsets for each of the bond idx
    // Use offsets as accumulators for length
    for (size_t i = 0; i < bond_count; ++i) {
        local_offset[i].off[0] = (uint16_t)out_conn->offset[bond_pairs[i].idx[0]]++;
        local_offset[i].off[1] = (uint16_t)out_conn->offset[bond_pairs[i].idx[1]]++;
    }

    // Compute complete edge offsets (exclusive scan)
    uint32_t off = 0;
    for (size_t i = 0; i < out_conn->offset_count; ++i) {
        const uint32_t len = out_conn->offset[i];
        out_conn->offset[i] = off;
        off += len;
    }

    // Write edge indices to correct location
    for (size_t i = 0; i < bond_count; ++i) {
        const md_atom_pair_t p = bond_pairs[i];
        const int atom_a = p.idx[0];
        const int atom_b = p.idx[1];
        const int local_a = (int)local_offset[i].off[0];
        const int local_b = (int)local_offset[i].off[1];
        const int off_a = out_conn->offset[atom_a];
        const int off_b = out_conn->offset[atom_b];

        const int idx_a = off_a + local_a;
        const int idx_b = off_b + local_b;

        ASSERT(idx_a < (int)out_conn->count);
        ASSERT(idx_b < (int)out_conn->count);

        // Store the cross references to the 'other' atom index signified by the bond in the correct location
        out_conn->atom_idx[idx_a] = atom_b;
        out_conn->atom_idx[idx_b] = atom_a;

        out_conn->bond_idx[idx_a] = (md_bond_idx_t)i;
        out_conn->bond_idx[idx_b] = (md_bond_idx_t)i;
    }

    md_vm_arena_temp_end(temp);
}

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

static const int max_valence_bonds[] = {
    0,  // Z = 0 (placeholder, no element with Z = 0)
    1,  // H
    0,  // He
    1,  // Li
    2,  // Be
    3,  // B
    4,  // C
    4,  // N
    3,  // O
    1,  // F
    0,  // Ne
    1,  // Na
    2,  // Mg
    4,  // Al
    4,  // Si
    5,  // P
    6,  // S
    7,  // Cl
    4,  // Ar
    1,  // K
    2,  // Ca
    6,  // Sc
    6,  // Ti
    6,  // V
    6,  // Cr
    6,  // Mn
    6,  // Fe
    6,  // Co
    6,  // Ni
    4,  // Cu
    4,  // Zn
    3,  // Ga
    4,  // Ge
    5,  // As
    6,  // Se
    7,  // Br
    4,  // Kr
    1,  // Rb
    2,  // Sr
    6,  // Y
    6,  // Zr
    6,  // Nb
    6,  // Mo
    6,  // Tc
    6,  // Ru
    6,  // Rh
    6,  // Pd
    4,  // Ag
    4,  // Cd
    3,  // In
    4,  // Sn
    5,  // Sb
    6,  // Te
    7,  // I
    4,  // Xe
    1,  // Cs
    2,  // Ba
    6,  // La
    6,  // Ce
    6,  // Pr
    6,  // Nd
    6,  // Pm
    6,  // Sm
    6,  // Eu
    6,  // Gd
    6,  // Tb
    6,  // Dy
    6,  // Ho
    6,  // Er
    6,  // Tm
    6,  // Yb
    6,  // Lu
    6,  // Hf
    6,  // Ta
    6,  // W
    6,  // Re
    6,  // Os
    6,  // Ir
    6,  // Pt
    4,  // Au
    4,  // Hg
    3,  // Tl
    4,  // Pb
    5,  // Bi
    6,  // Po
    7,  // At
    4,  // Rn
    1,  // Fr
    2,  // Ra
    6,  // Ac
    6,  // Th
    6,  // Pa
    6,  // U
    6,  // Np
    6,  // Pu
    6,  // Am
    6,  // Cm
    6,  // Bk
    6,  // Cf
    6,  // Es
    6,  // Fm
    6,  // Md
    6,  // No
    6,  // Lr
    6,  // Rf
    6,  // Db
    6,  // Sg
    6,  // Bh
    6,  // Hs
    6,  // Mt
    6,  // Ds
    6,  // Rg
    6,  // Cn
    6,  // Nh
    6,  // Fl
    6,  // Mc
    6,  // Lv
    6,  // Ts
    6   // Og
};

static inline int max_neighbors_element(md_element_t elem) {
    const int idx = MIN(elem, ARRAY_SIZE(max_valence_bonds)-1);
    return max_valence_bonds[idx];
}

static inline float angle(vec3_t a, vec3_t b, vec3_t c) {
    vec3_t v0 = vec3_normalize(vec3_sub(a, b));
    vec3_t v1 = vec3_normalize(vec3_sub(c, b));
    return (float)RAD_TO_DEG(acosf(vec3_dot(v0, v1)));
}

static inline size_t extract_neighborhood(md_atom_idx_t out_idx[], size_t out_cap, const md_bond_data_t* bond, size_t atom_idx, size_t depth) {
    ASSERT(out_idx);
    ASSERT(out_cap > 0);

#if DEBUG
    MEMSET(out_idx, -1, sizeof(md_atom_idx_t) * out_cap);
#endif

    size_t out_len = 0;
    out_idx[out_len++] = (md_atom_idx_t)atom_idx;
    if (depth == 0) return out_len;

    size_t beg = 0;
    size_t end = 1;
    md_atom_idx_t parent = -1;

    while (depth-- > 0) {
        size_t pre = out_len;
        for (size_t i = beg; i < end; ++i) {
            md_atom_idx_t cur_idx = out_idx[i];
            uint32_t conn_len = (uint32_t)md_bond_conn_count(*bond, cur_idx);
            if (conn_len == 1) continue;
            uint32_t conn_i   = bond->conn.offset[cur_idx];
            for (uint32_t j = 0; j < conn_len; ++j) {
                md_atom_idx_t idx = md_bond_conn_atom_idx(*bond, conn_i, j);
                if (idx != parent) {
                    out_idx[out_len++] = idx;
                    if (out_len == out_cap) {
                        MD_LOG_DEBUG("Maximum capacity was reached, but perhaps not all of the neighborhood was written");
                        return out_len;
                    }
                }
            }
        }
        parent = out_idx[0];
        beg = pre;
        end = out_len;
    }
    return out_len;
}

static size_t find_isomorphisms(md_index_data_t*, const graph_t* ,const graph_t* ,md_util_match_mode_t ,int , md_allocator_i*, md_allocator_i*);
static void find_isomorphisms_callback(const graph_t*, const graph_t*, int, state_t*);

static graph_t smiles_to_graph(str_t smiles_str, md_util_match_flags_t flags, md_allocator_i* arena, md_allocator_i* temp_arena) {
    typedef struct vertex_t {
        int type;
        int edge_count;
        int edge_idx[8];
        int edge_type[8];
    } vertex_t;

    size_t node_cap = str_len(smiles_str);
    md_smiles_node_t* nodes = md_vm_arena_push_array(temp_arena, md_smiles_node_t, node_cap);
    const size_t num_nodes = md_smiles_parse(nodes, node_cap, str_ptr(smiles_str), str_len(smiles_str));
    vertex_t* verts = md_vm_arena_push(temp_arena, sizeof(vertex_t) * num_nodes * 2);
    size_t num_verts = 0;

    graph_t graph = {0};

    if (nodes) {
        int bridge_idx[128];
        int stack[128];
        int stack_size = 0;
        int hub = -1;
        int order = 0;

        MEMSET(bridge_idx, -1, sizeof(bridge_idx));

        for (size_t i = 0; i < num_nodes; ++i) {
            const md_smiles_node_t* node = &nodes[i];
            if (node->type == MD_SMILES_NODE_ATOM) {
                md_element_t elem = node->atom.element;

                if ((flags & MD_UTIL_MATCH_FLAGS_NO_H) && elem == H) {
                    continue;
                }

                int cur = (int)num_verts;
                verts[num_verts++] = (vertex_t){.type = elem};

                if (hub != -1) {
                    verts[cur].edge_idx[verts[cur].edge_count] = hub;
                    verts[cur].edge_type[verts[cur].edge_count] = order;

                    verts[hub].edge_idx[verts[hub].edge_count] = cur;
                    verts[hub].edge_type[verts[hub].edge_count] = order;

                    verts[cur].edge_count++;
                    verts[hub].edge_count++;

                    // reset order state
                    order = 0;
                }
                hub = cur;

                if ( (flags & MD_UTIL_MATCH_FLAGS_NO_H) ||
                    ((flags & MD_UTIL_MATCH_FLAGS_NO_CH) && elem == C)) {
                    continue;
                }

                for (int j = 0; j < node->atom.hydrogen_count; ++j) {
                    const int h_idx = (int)num_verts;
                    const int edge_type = 1;
                    verts[num_verts++] = (vertex_t){.type = H};

                    verts[cur].edge_idx[verts[cur].edge_count] = h_idx;
                    verts[cur].edge_type[verts[cur].edge_count] = edge_type;

                    verts[h_idx].edge_idx[verts[h_idx].edge_count] = cur;
                    verts[h_idx].edge_type[verts[h_idx].edge_count] = edge_type;

                    verts[cur].edge_count++;
                    verts[h_idx].edge_count++;
                }
            }
            else if (node->type == MD_SMILES_NODE_BOND) {
                switch (node->bond.symbol) {
                case '-':
                case '/':
                case '\\':
                    order = 1; break;
                case '=': order = 2; break;
                case '#': order = 3; break;
                case '$': order = 4; break;
                case ':': order = MD_BOND_FLAG_AROMATIC; break;
                default: ASSERT(false); break;
                }
            }
            else if (node->type == MD_SMILES_NODE_BRANCH_OPEN) {
                if (stack_size < (int)ARRAY_SIZE(stack)) {
                    stack[stack_size++] = hub;
                } else {
                    MD_LOG_ERROR("Branch stack overflow in SMILES string");
                    goto done;
                }
            }
            else if (node->type == MD_SMILES_NODE_BRANCH_CLOSE) {
                if (stack_size > 0) {
                    hub = stack[--stack_size];
                } else {
                    MD_LOG_ERROR("Unbalanced branch close in SMILES string");
                    goto done;
                }
            }
            else if (node->type == MD_SMILES_NODE_BRIDGE) {
                ASSERT(node->bridge.index < ARRAY_SIZE(bridge_idx));
                int pi = bridge_idx[node->bridge.index];
                int ci = hub;
                if (pi != -1) {
                    verts[ci].edge_idx[verts[ci].edge_count] = pi;
                    verts[ci].edge_type[verts[ci].edge_count] = 0;

                    verts[pi].edge_idx[verts[pi].edge_count] = ci;
                    verts[pi].edge_type[verts[pi].edge_count] = 0;

                    verts[ci].edge_count++;
                    verts[pi].edge_count++;
                    bridge_idx[node->bridge.index] = -1;
                } else {
                    bridge_idx[node->bridge.index] = ci;
                }
            }
        }
    }

    if (num_verts > 0) {
        size_t tot_edge_count = 0;
        for (size_t i = 0; i < num_verts; ++i) {
            tot_edge_count += verts[i].edge_count;
        }
        graph.vertex_count = num_verts;
        graph.vertex_type = md_vm_arena_push_array     (arena, uint8_t,  num_verts);
        graph.edge_offset = md_vm_arena_push_zero_array(arena, uint32_t, num_verts + 1);
        graph.edge_data   = md_vm_arena_push_array     (arena, uint32_t, tot_edge_count);

        size_t edge_idx = 0;
        for (size_t i = 0; i < num_verts; ++i) {
            graph.vertex_type[i] = (uint8_t)verts[i].type;
            if (verts[i].edge_count > 0) {
                const uint32_t len = verts[i].edge_count;
                graph.edge_offset[i] = len;
                for (uint32_t j = 0; j < len; ++j) {
                    const uint32_t idx  = verts[i].edge_idx[j];
                    const uint32_t type = verts[i].edge_type[j] & 0x7F;
                    graph.edge_data[edge_idx++] = (type << 24) | idx;
                }
            }
        }

        // exclusive scan to set the offsets
        uint32_t offset = 0;
        for (size_t i = 0; i < num_verts + 1; ++i) {
            uint32_t len = graph.edge_offset[i];
            graph.edge_offset[i] = offset;
            offset += len;
        }
    }

done:
    return graph;
}

static inline bool element_is_aromatic(md_element_t elem) {
    switch (elem) {
    case C:
    case N:
    case O:
    case S:
    case Se:
        return true;
    default:
        return false;
    }
}

static inline uint64_t compute_key(const md_bond_data_t* bond, const md_element_t* atom_element, uint32_t conn_beg, uint32_t conn_end) {
    uint64_t key = 0;
    for (uint32_t i = conn_beg; i < conn_end; ++i) {
        md_element_t elem = atom_element[bond->conn.atom_idx[i]];
        md_order_t  order = bond->order[bond->conn.bond_idx[i]];
        key += elem * order;
    }
    return key;
}

typedef struct {
    md_element_t node_type;
    md_order_t   ring_order_sum;
    md_order_t   ext_order;
    md_element_t ext_type;
} ring_pattern_t;

static inline bool compare_ring_pattern(const ring_pattern_t* pat, const ring_pattern_t* node) {
    if (pat->node_type != node->node_type) return false;
    if (pat->ring_order_sum != node->ring_order_sum) return false;
    if (pat->ext_order != node->ext_order) return false;
    if (pat->ext_type && pat->ext_type != node->ext_type) return false;
    return true;
}

typedef struct {
    const graph_t* neighborhood;
    const graph_t* pattern;
    md_bond_data_t* bond;
} pattern_match_payload_t;

static bool pattern_match_callback(const int map[], size_t length, void* user) {
    pattern_match_payload_t* data = (pattern_match_payload_t*)user;

    uint32_t p_edge_beg = data->pattern->edge_offset[0];
    uint32_t p_edge_end = data->pattern->edge_offset[1];

    //md_atom_idx_t na = map[0];
    for (uint32_t i = p_edge_beg; i < p_edge_end; ++i) {
        md_order_t order = (md_order_t)graph_edge_type(data->pattern, i);
        if (order > 1) {
            md_atom_idx_t nb = map[graph_edge_vertex_idx(data->pattern, i)];
            uint32_t n_edge_beg = data->neighborhood->edge_offset[0];
            uint32_t n_edge_end = data->neighborhood->edge_offset[1];
            for (uint32_t j = n_edge_beg; j < n_edge_end; ++j) {
                if (graph_edge_vertex_idx(data->neighborhood, j) == nb) {
                    md_bond_idx_t bidx = data->neighborhood->bond_idx_map[j];
                    md_bond_order_set(data->bond, bidx, order);
                    break;
                }
            }
        }
    }

    return true;
}

static int graph_depth(const graph_t* graph, int start_idx, md_allocator_i* temp_arena) {
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_arena);

    int* depth = md_vm_arena_push_zero_array(temp_arena, int, graph->vertex_count);
    fifo_t fifo = fifo_create(32, temp_arena);

    fifo_push(&fifo, start_idx);
    while (!fifo_empty(&fifo)) {
        int idx = fifo_pop(&fifo);
        int d = depth[idx] + 1;
        graph_edge_iter_t it = graph_edge_iter(graph, idx);
        while (graph_edge_iter_has_next(it)) {
            int other_idx = graph_edge_iter_vidx(it);
            // @NOTE: We do not care about the minimum depth, for each node
            if (other_idx != start_idx && depth[other_idx] == 0) {
                depth[other_idx] = d;
                fifo_push(&fifo, other_idx);
            }
            graph_edge_iter_next(&it);
        }
    }

    int max_depth = 0;
    for (size_t i = 0; i < graph->vertex_count; ++i) {
        max_depth = MAX(max_depth, depth[i]);
    }

    md_vm_arena_temp_end(temp);

    return max_depth;
}

static bool compute_covalent_bond_order(md_bond_data_t* bond, const md_atom_data_t* atom, const md_index_data_t* rings) {
    if (!bond || bond->count == 0 || !atom || atom->count == 0) {
        return false;
    }

#if 0
    if (!type || !resname) {
        MD_LOG_DEBUG("No atom type or atom resname were given, will default to covalent order of 1 for all bonds.");
        STATIC_ASSERT(sizeof(md_order_t) == 1, "Incorrect size of md_order_t");
        MEMSET(bond_order, 1, md_array_bytes(bond_order));
        return true;
    }
    
    for (int64_t i = 0; i < bond_count; ++i) {
        const int atom_a = bond_pairs[i].idx[0];
        const int atom_b = bond_pairs[i].idx[1];
        str_t type_a = LBL_TO_STR(type[atom_a]);
        str_t type_b = LBL_TO_STR(type[atom_b]);
        str_t comp_a = LBL_TO_STR(resname[atom_a]);
        str_t comp_b = LBL_TO_STR(resname[atom_b]);

        md_order_t order = 1;
        if (str_eq(comp_a, comp_b)) {
            if (str_cmp_lex(type_a, type_b) > 0) {
                str_swap(type_a, type_b);
            }
            // Intra
            if (str_eq_cstr(type_a, "C") && str_eq_cstr(type_b, "O") && md_util_resname_amino_acid(comp_a)) {
                order = 2;
            } else {
                char buf[32];
                int len = build_key(buf, comp_a, type_a, type_b);
                if (find_str_in_str_arr(NULL, (str_t){buf,len}, intra_bond_order_table, ARRAY_SIZE(intra_bond_order_table))) {
                    bond_order[i] = 2;
                }
            }
        } else {
            // Inter
            if ( (str_eq_cstr(comp_a, "LYS") && str_eq_cstr(type_a, "CZ") && str_eq_cstr(comp_b, "RET") && str_eq_cstr(type_b, "C15")) ||
                 (str_eq_cstr(comp_b, "LYS") && str_eq_cstr(type_b, "CZ") && str_eq_cstr(comp_a, "RET") && str_eq_cstr(type_a, "C15")) ){
                bond_order[i] = 2;
            }
        }
        bond_order[i] = order;
    }
#else
    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));
    uint8_t* type = md_vm_arena_push_zero_array(temp_arena, uint8_t, atom->count);

    // Identify Sp Sp2 and Sp3 types
    for (size_t i = 0; i < atom->count; ++i) {
        uint32_t   conn_len = (uint32_t)md_bond_conn_count(*bond, i);
        uint32_t   conn_i   = bond->conn.offset[i];

        if (conn_len < 2) continue;
        vec3_t xi = md_atom_coord(*atom, i);

        switch (conn_len) {
        case 2: {
            uint32_t j = md_bond_conn_atom_idx(*bond, conn_i, 0);
            uint32_t k = md_bond_conn_atom_idx(*bond, conn_i, 1);
            vec3_t xj = md_atom_coord(*atom, j);
            vec3_t xk = md_atom_coord(*atom, k);
            float theta = angle(xj, xi, xk);
            if (theta > 155.f) {
                type[i] = 1;
            } else if (theta > 100.f) {
                type[i] = 2;
            }
            break;
        }
        case 3:
        case 4: {
            float theta = 0.0f;
            for (uint32_t idx = 0; idx < conn_len; ++idx) {
                uint32_t j = md_bond_conn_atom_idx(*bond, conn_i, idx);
                uint32_t k = md_bond_conn_atom_idx(*bond, conn_i, (idx + 1) % conn_len);
                vec3_t xj = md_atom_coord(*atom, j);
                vec3_t xk = md_atom_coord(*atom, k);
                theta += angle(xj, xi, xk);
            }
            theta /= (float)conn_len;
            if (theta > 115.f) {
                type[i] = 2;
            } else {
                type[i] = 3;
            }
            break;
        }
        default:
            continue;
        }
    }

    // 6b
    if (rings) {
        size_t num_rings = md_index_data_num_ranges(*rings);
        for (size_t ring_idx = 0; ring_idx < num_rings; ++ring_idx) {
            size_t ring_size = md_index_range_size(*rings, ring_idx);
            md_atom_idx_t* atom_idx = md_index_range_beg(*rings, ring_idx);
            bool set_c2_to_sp2 = false;
            if (ring_size == 5) {
                vec3_t c[5];
                for (size_t i = 0; i < 5; ++i) {
                    c[i] = md_atom_coord(*atom, atom_idx[i]);
                }

                vec3_t d[5] = {
                    vec3_sub(c[1], c[0]),
                    vec3_sub(c[2], c[1]),
                    vec3_sub(c[3], c[2]),
                    vec3_sub(c[4], c[3]),
                    vec3_sub(c[0], c[4]),
                };

                // @TODO: Perform min image on d
                
                double avg_angle = 0;
                avg_angle += fabs(vec3_dihedral_angle(d[0], d[1], d[2]));
                avg_angle += fabs(vec3_dihedral_angle(d[1], d[2], d[3]));
                avg_angle += fabs(vec3_dihedral_angle(d[2], d[3], d[4]));
                avg_angle += fabs(vec3_dihedral_angle(d[3], d[4], d[0]));
                avg_angle += fabs(vec3_dihedral_angle(d[4], d[0], d[1]));
                avg_angle = RAD_TO_DEG(avg_angle) / 5.0;

                set_c2_to_sp2 = (avg_angle < 7.5);
            } else if (ring_size == 6) {
                vec3_t c[6];
                for (size_t i = 0; i < 6; ++i) {
                    c[i] = md_atom_coord(*atom, atom_idx[i]);
                }

                vec3_t d[6] = {
                    vec3_sub(c[1], c[0]),
                    vec3_sub(c[2], c[1]),
                    vec3_sub(c[3], c[2]),
                    vec3_sub(c[4], c[3]),
                    vec3_sub(c[0], c[4]),
                };

                // @TODO: Perform min image on d

                double avg_angle = 0;
                avg_angle += fabs(vec3_dihedral_angle(d[0], d[1], d[2]));
                avg_angle += fabs(vec3_dihedral_angle(d[1], d[2], d[3]));
                avg_angle += fabs(vec3_dihedral_angle(d[2], d[3], d[4]));
                avg_angle += fabs(vec3_dihedral_angle(d[3], d[4], d[5]));
                avg_angle += fabs(vec3_dihedral_angle(d[4], d[5], d[0]));
                avg_angle += fabs(vec3_dihedral_angle(d[5], d[0], d[1]));
                avg_angle = RAD_TO_DEG(avg_angle) / 6.0;

                set_c2_to_sp2 = (avg_angle < 12.0);
            }

            if (set_c2_to_sp2) {
                for (size_t i = 0; i < ring_size; ++i) {
                    if (md_bond_conn_count(*bond, atom_idx[i]) == 2) {
                        type[atom_idx[i]] = 2;
                    }
                }
            }
        }
    }

    // 6c
    for (size_t i = 0; i < atom->count; ++i) {
        if (type[i] == 1) {
            md_bond_iter_t it = md_bond_iter(bond, i);
            while (md_bond_iter_has_next(it)) {
                // @TODO: Check valence of terminal atom
                if (type[md_bond_iter_atom_index(it)] == 1 ||
                    type[md_bond_iter_atom_index(it)] == 0) {
                    goto next_atom;
                }
                md_bond_iter_next(&it);
            }
            type[i] = 2;
        } else if (type[i] == 2) {
            md_bond_iter_t it = md_bond_iter(bond, i);
            while (md_bond_iter_has_next(it)) {
                // @TODO: Check valence of terminal atom
                if (type[md_bond_iter_atom_index(it)] == 2 ||
                    type[md_bond_iter_atom_index(it)] == 0) {
                    goto next_atom;
                }
                md_bond_iter_next(&it);
            }
            type[i] = 3;
        }
        next_atom:;
    }


    md_timestamp_t t0 = md_time_current();

    // 7
    // Patterns of functional groups to look for as graphs
    // The patterns are not strictly matched on the bond order, only on the element (if present) and if its terminal
    // The pattern encodes bond order which is applied if the pattern can be matched

    str_t functional_group_patterns[] = {
        STR_LIT("[CH](=O)N(*)(*)"),
        STR_LIT("C(=O)([OH])(*)"),
        STR_LIT("C(=O)(O*)(*)"),
        STR_LIT("C(=O)(S*)(*)"),
        STR_LIT("C(=O)([NH]*)(*)"),
        STR_LIT("C(=S)(S*)(*)"),
        STR_LIT("C(=S)([NH]*)(*)"),
        STR_LIT("C(=N*)([NH]*)([NH]*)"),
        STR_LIT("C(=[NH])([NH2])(*)"),
        STR_LIT("[N+](=[N-])(=N*)"),
        STR_LIT("P(=O)(*)(*)(*)"),
        STR_LIT("S(=O)(=O)(*)(*)"),
        STR_LIT("N(=O)(=O)(*)"),
        STR_LIT("[Se](=O)([OH])(*)"),
    };

    typedef struct {
        graph_t graph;
        int depth;
        uint64_t mask;
        uint8_t bond_count[128];
    } func_group_t;

    func_group_t func_group[ARRAY_SIZE(functional_group_patterns)] = {0};
    uint8_t  element_min_connectivity[128];
    uint8_t  element_max_connectivity[128];

    MEMSET(element_min_connectivity, 255, sizeof(element_min_connectivity));
    MEMSET(element_max_connectivity,   0, sizeof(element_max_connectivity));

    for (size_t i = 0; i < ARRAY_SIZE(functional_group_patterns); ++i) {
        graph_t graph = smiles_to_graph(functional_group_patterns[i], 0, temp_arena, temp_arena);
        func_group[i].graph = graph;
        func_group[i].depth = graph_depth(&graph, 0, temp_arena);
        md_element_t elem = graph.vertex_type[0];
        element_min_connectivity[elem] = MIN(element_min_connectivity[elem], (uint8_t)graph_vertex_edge_count(&graph, 0));
        element_max_connectivity[elem] = MAX(element_max_connectivity[elem], (uint8_t)graph_vertex_edge_count(&graph, 0));

        for (size_t j = 0; j < graph.vertex_count; ++j) {
            func_group[i].mask |= ((uint64_t)1 << (elem & 63));
        }
        graph_edge_iter_t it = graph_edge_iter(&graph, 0);
        while (graph_edge_iter_has_next(it)) {
            int et = graph_edge_iter_type(it);
            func_group[i].bond_count[et] += 1;
            graph_edge_iter_next(&it);
        }
    }

    for (size_t i = 0; i < atom->count; ++i) {
        md_element_t elem_i = atom->element[i];
        uint32_t conn_beg = bond->conn.offset[i];
        uint32_t conn_end = bond->conn.offset[i+1];
        uint32_t conn_len = conn_end - conn_beg;

        if (element_min_connectivity[elem_i] == 0 || conn_len < element_min_connectivity[elem_i] || element_max_connectivity[elem_i] < conn_len) {
            continue;
        }

        uint64_t mask = 0;
        for (uint32_t j = conn_beg; j < conn_end; ++j) {
            md_atom_idx_t idx = bond->conn.atom_idx[j];
            uint64_t elem = atom->element[idx];
            mask |= ((uint64_t)1 << (elem & 63));
        }


        int patterns_to_test[ARRAY_SIZE(func_group)];
        int pattern_count = 0;
        int max_depth = 0;

        for (int j = 0; j < (int)ARRAY_SIZE(func_group); ++j) {
            int    vert_type  = graph_vertex_type(&func_group[j].graph, 0);
            size_t edge_count = graph_vertex_edge_count(&func_group[j].graph, 0);
            if (elem_i == vert_type &&
                conn_len == edge_count &&
                (mask & func_group[j].mask) == func_group[j].mask)
            {
                patterns_to_test[pattern_count++] = j;
                max_depth = MAX(max_depth, func_group[j].depth);
            }
        }

        if (pattern_count == 0) continue;

        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_arena);

        // Extract neighborhood of atom i
        md_atom_idx_t n_idx[32];
        size_t n_len = extract_neighborhood(n_idx, ARRAY_SIZE(n_idx), bond, i, max_depth);
        ASSERT(n_len > 0);

        graph_t neighborhood = make_graph(bond, atom->element, n_idx, n_len, temp_arena);

        uint64_t n_mask = 0;
        for (size_t j = 0; j < neighborhood.vertex_count ; ++j) {
            uint64_t elem = graph_vertex_type(&neighborhood, j);
            n_mask |= ((uint64_t)1 << (elem & 63));
        }

        pattern_match_payload_t payload = {
            .neighborhood = & neighborhood,
            .pattern = 0,
            .bond = bond,
        };

        for (int j = 0; j < pattern_count; ++j) {
            int fi = patterns_to_test[j];
            const graph_t* pattern = &func_group[fi].graph;

            if ((n_mask & func_group[j].mask) != func_group[j].mask) continue;

            payload.pattern = pattern;

            state_t state = {0};
            state_init(&state, pattern, &neighborhood, temp_arena);
            state.callback = pattern_match_callback;
            state.user_data = &payload;

            find_isomorphisms_callback(pattern, &neighborhood, pattern->vertex_type[0], &state);
        }

        md_vm_arena_temp_end(temp);
    }

    md_timestamp_t t1 = md_time_current();
    printf("Time for functional group perception: %.3f ms\n", md_time_as_milliseconds(t1-t0));

    // 8. Aromatic Ring Perception
    if (rings) {
        const ring_pattern_t ring_patterns[] = {
            {6, 2, 1, 0},
            {6, 3, 1, 0},
            {6, 2, 2, 0},
            {6, 2, 0, 0},
            {6, 3, 0, 0},

            {6, 2, 1, 8},
            {6, 2, 2, 8},
            {6, 3, 1, 8},

            {7, 3, 0, 0},
            {7, 2, 2, 8},
            {7, 3, 2, 0},
            {7, 2, 1, 0},
            {7, 3, 1, 0},

            {7, 2, 0, 0},
            {7, 2, 1, 0},

            {8, 2, 0, 0},
            {16, 2, 0, 0},
            {34, 2, 0, 0},
        };

        const uint8_t pattern_electron_count[20] = {
            1,1,1,1,1,
            1,0,1,
            1,1,1,1,1,
            1,2,
            2,2,2,
        };

        const bool pattern_incident_multi_bond[20] = {
            false, true, true, false, true,
            false, true, true,
            true, true, true, false, true,
            false, false,
            false, false, false,
        };

        t0 = md_time_current();
        size_t num_rings = md_index_data_num_ranges(*rings);
        for (size_t ring_idx = 0; ring_idx < num_rings; ++ring_idx) {
            md_atom_idx_t* atom_beg = md_index_range_beg(*rings, ring_idx);
            md_atom_idx_t* atom_end = md_index_range_end(*rings, ring_idx);
            size_t ring_size        = md_index_range_size(*rings, ring_idx);

            if (ring_size < 5 && 6 < ring_size) continue;

            // First iterate to see that all atoms are sp2 and of correct element
            for (md_atom_idx_t *it = atom_beg; it != atom_end; ++it) {
                if (type[*it] != 2 || !is_aromatic(atom->element[*it])) {
                    goto next_ring;
                }
            }

            uint8_t pidx[8];
            int electron_count = 0;
            for (int *it = atom_beg, j = 0; it != atom_end; ++it, ++j) {
                int i = *it;
                md_element_t elem_i = atom->element[i];

                ring_pattern_t p = {0};
                p.node_type = atom->element[i];
                
                uint32_t conn_beg = bond->conn.offset[i];
                uint32_t conn_end = bond->conn.offset[i+1];
                uint32_t conn_len = conn_end - conn_beg;

                if (conn_len < 2 || 3 < conn_len) continue;

                for (uint32_t k = conn_beg; k < conn_end; ++k) {
                    md_atom_idx_t idx = bond->conn.atom_idx[k];
                    bool ring_bond = false;
                    if (i == atom_beg[0]) {
                        ring_bond = idx == atom_end[-1] || idx == it[1];
                    } else if (i == atom_end[-1]) {
                        ring_bond = idx == atom_beg[0]  || idx == it[-1];
                    } else {
                        ring_bond = idx == it[-1] || idx == it[1];
                    }
                    md_element_t elem = atom->element[bond->conn.atom_idx[k]];
                    md_order_t  order = bond->order  [bond->conn.bond_idx[k]];
                    if (ring_bond) {
                        p.ring_order_sum += order;
                    } else {
                        p.ext_type  = elem;
                        p.ext_order = order;
                    }
                }

                for (size_t k = 0; k < ARRAY_SIZE(ring_patterns); ++k) {
                    if (compare_ring_pattern(&ring_patterns[k], &p)) {
                        pidx[j] = (uint8_t)k;
                        goto next;
                    }
                }
                // No match found
                goto next_ring;
                next:;
            }

            int ambibous_c[8];
            int ambigous_c_count = 0;

            int ambigous_n[8];
            int ambigous_n_count = 0;

            int dbl_bond_targets[4];
            int dbl_bond_target_count = 0;

            if (ring_idx == 15) {
                while(0) {};
            }

            // 8a: Try to resolve ambigous cases, by looking at neighbors
            for (int *it = atom_beg, j = 0; it != atom_end; ++it, ++j) {
                if (*it == 801) {
                    while(0) {};
                }
                int pi = pidx[j];
                if (pi == 5 || pi == 13) {
                    int prev_pat_idx = (j == 0) ? pidx[ring_size-1] : pidx[j-1];
                    int next_pat_idx = (j == ring_size-1) ? pidx[0] : pidx[j+1];

                    int prev_atom_idx = it == atom_beg   ? atom_end[-1] : it[-1];
                    int next_atom_idx = it == atom_end-1 ? atom_beg[ 0] : it[ 1];

                    // @TODO: test valance of neighbors
                    if (pattern_incident_multi_bond[prev_pat_idx] && pattern_incident_multi_bond[next_pat_idx]) {
                        if (pi == 5) {
                            electron_count += 1;
                            dbl_bond_targets[dbl_bond_target_count++] = *it;
                        } else {
                            electron_count += 2;
                        }
                    } else {
                        if (pi == 5) {
                            ambibous_c[ambigous_c_count++] = *it;
                            electron_count += 1;
                        } else {
                            ambigous_n[ambigous_n_count++] = *it;
                            electron_count += 1;
                        }
                    }
                } else {
                    electron_count += pattern_electron_count[pi];
                }
            }
            
            // 8b
            if ((electron_count & 3) == 1 && ambigous_n_count > 0) {
                // convert *-N-* to *-[NH]-*
                electron_count += 1;
                ambigous_n_count -= 1;
            }

            // 8c
            if ((electron_count & 3) == 3 && ambigous_c_count > 0) {
                // convert *-[C](O)-* to *-C(=O)-*
                electron_count -= 1;
                dbl_bond_targets[dbl_bond_target_count++] = ambibous_c[0];
            }
            
            // 8d
            if ((electron_count & 3) == 3 && ambigous_n_count > 0) {
                // convert *-[NH]-* to *-[NH+]-*
                electron_count -= 1;
            }

            // 8e
            if ((electron_count & 3) == 2) {
                for (md_atom_idx_t *it = atom_beg; it != atom_end; ++it) {
                    atom->flags[*it] |= MD_FLAG_AROMATIC;
                }
            }

            // Set double bonds
            for (int i = 0; i < dbl_bond_target_count; ++i) {
                md_atom_idx_t idx = dbl_bond_targets[i];
                uint32_t conn_beg = bond->conn.offset[i];
                uint32_t conn_end = bond->conn.offset[i+1];
                for (uint32_t j = conn_beg; j < conn_end; ++j) {
                    if (atom->element[bond->conn.atom_idx[j]] == O) {
                        bond->order[bond->conn.bond_idx[j]] = 2;
                        break;
                    }
                }
            }

        next_ring:;
        }
        t1 = md_time_current();
        printf("Time for aromatic ring perception: %.3f ms\n", md_time_as_milliseconds(t1-t0));
    }

    for (size_t i = 0; i < atom->count; ++i) {
        if (type[i] == 1) {
            atom->flags[i] |= MD_FLAG_SP;
        } else if (type[i] == 2) {
            atom->flags[i] |= MD_FLAG_SP2;
        } else if (type[i] == 3) {
            atom->flags[i] |= MD_FLAG_SP3;
        }
    }

    for (size_t i = 0; i < bond->count; ++i) {
        if (atom->flags[bond->pairs[i].idx[0]] & MD_FLAG_AROMATIC &&
            atom->flags[bond->pairs[i].idx[1]] & MD_FLAG_AROMATIC) {
            bond->order[i] |= MD_BOND_FLAG_AROMATIC;
        }
    }

    md_vm_arena_destroy(temp_arena);
#endif
    
    return true;
}

// Returns the bond-order (0 if no bond is present)
static inline uint8_t covalent_bond_heuristic2(float d2, md_element_t a, md_element_t b) {
    const uint8_t va = element_max_valence[a];
    const uint8_t vb = element_max_valence[b];

    const int v[3] = {
        (int)(va > 0 && vb > 0),
        (int)(va > 1 && vb > 1),
        (int)(va > 2 && vb > 2),
    };

    const float r[3] = {
        (element_covalent_radii3[a][0] * 0.01f + element_covalent_radii3[b][0] * 0.01f),
        (element_covalent_radii3[a][1] * 0.01f + element_covalent_radii3[b][1] * 0.01f),
        (element_covalent_radii3[a][2] * 0.01f + element_covalent_radii3[b][2] * 0.01f),
    };

    const float x = sqrtf(d2);

    const float k[3] = {-10, -20, -40}; 
    const float z[3] = {
        (x - r[0]),
        (x - r[1]),
        (x - r[2]),
    };
    const float y[3] = {
        v[0] * (k[0] * z[0] * z[0] + 1),
        v[1] * (k[1] * z[1] * z[1] + 1),
        v[2] * (k[2] * z[2] * z[2] + 1),
    };

    uint8_t val = 0;
    float max_val = 0;
    for (int i = 0; i < 3; ++i) {
        if (y[i] > max_val) {
            max_val = y[i];
            val = (uint8_t)(i + 1);
        }
    }
    return val;
}

#define R_MIN 0.8f
#define R_MAX 0.45f // 0.3f ???

static inline bool bond_heuristic(float d2, float r_min, float r_max) {
    return (r_min * r_min) < d2 && d2 < (r_max * r_max);
}

static inline int simd_bond_heuristic(md_256 d2, md_256 r_min, md_256 r_max) {
    const md_256 mask_d2  = md_mm256_and_ps(md_mm256_cmpgt_ps(d2, md_mm256_mul_ps(r_min, r_min)), md_mm256_cmplt_ps(d2, md_mm256_mul_ps(r_max, r_max)));
    return md_mm256_movemask_ps(mask_d2);
}

static inline bool covalent_bond_heuristic(float d2, float ra, float rb) {
    const float r_sum = ra + rb;
    const float r_min = R_MIN;
    const float r_max = r_sum + R_MAX;
    return r_min * r_min < d2 && d2 < (r_max * r_max);
}

static inline int simd_covalent_bond_heuristic(md_256 d2, md_256 cov_rad_a, md_256 cov_rad_b) {
    const md_256 r_sum = md_mm256_add_ps(cov_rad_a, cov_rad_b);
    const md_256 r_min = md_mm256_set1_ps(R_MIN);
    const md_256 r_max = md_mm256_add_ps(r_sum, md_mm256_set1_ps(R_MAX));
    const md_256 mask_d2  = md_mm256_and_ps(md_mm256_cmpgt_ps(d2, md_mm256_mul_ps(r_min, r_min)), md_mm256_cmplt_ps(d2, md_mm256_mul_ps(r_max, r_max)));
    const md_256 mask_rad = md_mm256_and_ps(md_mm256_cmpgt_ps(cov_rad_a, md_mm256_setzero_ps()),  md_mm256_cmpgt_ps(cov_rad_b, md_mm256_setzero_ps()));
    return md_mm256_movemask_ps(md_mm256_and_ps(mask_d2, mask_rad));
}

static inline float distance_squared(vec4_t dx, const md_unit_cell_t* cell) {
    if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        const vec4_t p  = vec4_from_vec3(mat3_diag(cell->basis), 0);
        const vec4_t rp = vec4_from_vec3(mat3_diag(cell->inv_basis), 0);
        const vec4_t d  = vec4_sub(dx, vec4_mul(vec4_round(vec4_mul(dx, rp)), p));
        dx = vec4_blend(d, dx, vec4_cmp_eq(p, vec4_zero()));
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        minimum_image_triclinic(dx.elem, cell->basis.elem);
    }
    return vec4_length_squared(dx);
}

static inline md_256 simd_distance_squared(const md_256 dx[3], const md_unit_cell_t* cell) {
    md_256 d[3] = {
        dx[0],
        dx[1],
        dx[2],
    };
    if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        md_256 p[3] = {
            md_mm256_set1_ps(cell->basis.elem[0][0]),
            md_mm256_set1_ps(cell->basis.elem[1][1]),
            md_mm256_set1_ps(cell->basis.elem[2][2]),
        };
        md_256 rp[3] = {
            md_mm256_set1_ps(cell->inv_basis.elem[0][0]),
            md_mm256_set1_ps(cell->inv_basis.elem[1][1]),
            md_mm256_set1_ps(cell->inv_basis.elem[2][2]),
        };
        d[0] = md_mm256_minimage_ps(d[0], p[0], rp[0]);
        d[1] = md_mm256_minimage_ps(d[1], p[1], rp[1]);
        d[2] = md_mm256_minimage_ps(d[2], p[2], rp[2]);
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        simd_minimum_image_triclinic(d, cell->basis.elem);
    }
    return md_mm256_fmadd_ps(d[0], d[0], md_mm256_fmadd_ps(d[1], d[1], md_mm256_mul_ps(d[2], d[2])));
}

// Compile-time bond_allowed table
// Indexed as bond_allowed[NUM_ELEMENTS][NUM_ELEMENTS]
// Assumes md_element_t enum matches periodic order (H=1, He=2, ...)

// ---------------------- Covalent Bonds ----------------------
static const float cov_r_min[Num_Elements][Num_Elements] = {
    [0] = {0},  // dummy

    // Main group elements
    [H]  = { [H]=0.33f, [B]=0.33f, [C]=0.33f, [N]=0.33f, [O]=0.33f, [F]=0.33f, [Si]=0.33f, [P]=0.33f, [S]=0.33f, [Cl]=0.33f, [Br]=0.33f, [I]=0.33f },
    [Li] = { [H]=1.68f, [C]=1.68f, [O]=1.68f, [N]=1.68f },
    [Be] = { [O]=1.42f, [N]=1.42f, [C]=1.42f },
    [B]  = { [H]=0.33f, [C]=0.99f, [N]=0.99f, [O]=0.99f, [F]=0.99f, [Si]=1.21f, [P]=1.21f, [S]=1.32f },
    [C]  = { [H]=0.33f, [C]=0.77f, [N]=0.77f, [O]=0.77f, [F]=0.77f, [Si]=1.16f, [P]=1.37f, [S]=0.77f, [Cl]=1.79f, [Br]=1.94f, [I]=2.14f, [B]=0.99f },
    [N]  = { [H]=0.33f, [C]=0.77f, [N]=0.77f, [O]=0.77f, [F]=0.77f, [P]=1.37f, [S]=1.58f, [B]=0.99f },
    [O]  = { [H]=0.33f, [C]=0.77f, [N]=0.77f, [O]=0.77f, [S]=1.40f, [P]=1.37f, [B]=0.99f },
    [F]  = { [H]=0.33f, [C]=0.77f, [N]=0.77f, [O]=0.77f, [F]=0.77f },
    [Na] = { [H]=1.89f, [O]=1.89f, [N]=1.89f },
    [Mg] = { [O]=1.99f, [N]=1.99f, [S]=2.04f },
    [Al] = { [H]=1.26f, [C]=1.58f, [N]=1.47f, [O]=1.47f, [Si]=1.47f, [P]=1.68f, [S]=1.84f },
    [Si] = { [H]=0.33f, [C]=1.16f, [O]=1.47f, [Si]=1.22f, [P]=1.63f, [S]=1.84f, [B]=1.21f },
    [P]  = { [H]=0.33f, [C]=1.37f, [N]=1.37f, [O]=1.37f, [S]=1.79f },
    [S]  = { [H]=0.33f, [C]=0.77f, [O]=1.40f, [P]=1.79f, [S]=1.84f },
    [Cl] = { [H]=0.33f, [C]=1.79f, [N]=1.63f, [O]=1.53f, [Cl]=0.95f },
    [K]  = { [H]=2.24f, [O]=2.24f, [N]=2.24f },
    [Ca] = { [O]=2.04f, [S]=2.14f },

    // Transition metals
    [Sc] = { [O]=1.68f, [N]=1.68f, [C]=1.78f },
    [Ti] = { [O]=1.53f, [N]=1.68f, [C]=1.78f },
    [V]  = { [O]=1.63f, [N]=1.73f, [C]=1.84f },
    [Cr] = { [O]=1.63f, [N]=1.73f, [C]=1.84f, [S]=1.84f },
    [Mn] = { [O]=1.63f, [N]=1.73f, [C]=1.84f, [S]=1.84f },
    [Fe] = { [O]=1.68f, [N]=1.68f, [C]=1.78f, [S]=1.73f },
    [Co] = { [O]=1.68f, [N]=1.68f, [C]=1.78f },
    [Ni] = { [O]=1.68f, [N]=1.68f, [C]=1.78f },
    [Cu] = { [O]=1.43f, [N]=1.43f, [S]=1.53f },
    [Zn] = { [O]=1.43f, [N]=1.43f, [S]=1.53f },
    [Ga] = { [O]=1.42f, [N]=1.42f, [C]=1.58f, [Ga]=1.22f },
    [Ge] = { [O]=1.43f, [C]=1.58f, [Ge]=1.22f },
    [As] = { [H]=0.33f, [C]=1.37f, [O]=1.37f, [S]=1.63f, [P]=1.53f },
    [Se] = { [H]=0.33f, [C]=1.58f, [O]=1.58f, [S]=1.73f },
    [Br] = { [H]=0.33f, [C]=1.99f, [Cl]=1.99f, [Br]=1.27f },
    [Rb] = { [O]=2.24f },
    [Sr] = { [O]=2.14f },
    [Y]  = { [O]=1.68f },
    [Zr] = { [O]=1.68f },
    [Nb] = { [O]=1.68f },
    [Mo] = { [O]=1.58f },
    [Tc] = { [O]=1.58f },
    [Ru] = { [O]=1.58f },
    [Rh] = { [O]=1.58f },
    [Pd] = { [O]=1.58f },
    [Ag] = { [O]=1.53f },
    [Cd] = { [O]=1.53f },
    [In] = { [O]=1.43f },
    [Sn] = { [O]=1.43f },
    [Sb] = { [O]=1.53f },
    [Te] = { [O]=1.53f },
    [I]  = { [H]=0.33f, [C]=2.14f, [Cl]=2.14f, [I]=1.40f },

    // Lanthanides
    [La] = { [O]=2.04f, [N]=2.04f },
    [Ce] = { [O]=2.04f, [N]=2.04f },
    [Pr] = { [O]=2.04f, [N]=2.04f },
    [Nd] = { [O]=2.04f, [N]=2.04f },
    [Sm] = { [O]=2.04f },
    [Eu] = { [O]=2.04f },
    [Gd] = { [O]=2.04f },
    [Tb] = { [O]=2.04f },
    [Dy] = { [O]=2.04f },
    [Ho] = { [O]=2.04f },
    [Er] = { [O]=2.04f },
    [Tm] = { [O]=2.04f },
    [Yb] = { [O]=2.04f },
    [Lu] = { [O]=2.04f },

    // Actinides
    [Th] = { [O]=2.04f },
    [Pa] = { [O]=2.04f },
    [U]  = { [O]=2.04f },
    [Np] = { [O]=2.04f },
    [Pu] = { [O]=2.04f },
};

static const float cov_r_max[Num_Elements][Num_Elements] = {
    [0] = {0},  // dummy

    // Main group elements
    [H]  = { [H]=0.951f, [B]=1.251f, [C]=1.251f, [N]=1.251f, [O]=1.251f, [F]=1.151f, [Si]=1.651f, [P]=1.651f, [S]=1.651f, [Cl]=1.651f, [Br]=1.951f, [I]=2.151f },
    [Li] = { [H]=2.151f, [C]=2.251f, [O]=2.151f, [N]=2.151f },
    [Be] = { [O]=1.951f, [N]=1.951f, [C]=1.951f },
    [B]  = { [H]=1.251f, [C]=1.651f, [N]=1.651f, [O]=1.551f, [F]=1.551f, [Si]=1.851f, [P]=1.851f, [S]=2.051f },
    [C]  = { [H]=1.251f, [C]=1.651f, [N]=1.651f, [O]=1.551f, [F]=1.551f, [Si]=1.951f, [P]=2.051f, [S]=1.951f, [Cl]=2.151f, [Br]=2.351f, [I]=2.651f, [B]=1.651f },
    [N]  = { [H]=1.251f, [C]=1.651f, [N]=1.651f, [O]=1.551f, [F]=1.551f, [P]=1.951f, [S]=2.051f, [B]=1.651f },
    [O]  = { [H]=1.251f, [C]=1.551f, [N]=1.551f, [O]=1.551f, [S]=2.051f, [P]=1.951f, [B]=1.551f },
    [F]  = { [H]=1.151f, [C]=1.551f, [N]=1.551f, [O]=1.551f, [F]=1.451f },
    [Na] = { [H]=2.451f, [O]=2.451f, [N]=2.451f },
    [Mg] = { [O]=2.251f, [N]=2.251f, [S]=2.351f },
    [Al] = { [H]=1.651f, [C]=2.151f, [N]=2.051f, [O]=2.051f, [Si]=2.151f, [P]=2.451f, [S]=2.551f },
    [Si] = { [H]=1.651f, [C]=1.951f, [O]=2.151f, [Si]=2.151f, [P]=2.451f, [S]=2.551f, [B]=1.851f },
    [P]  = { [H]=1.651f, [C]=2.051f, [N]=2.051f, [O]=1.951f, [S]=2.551f },
    [S]  = { [H]=1.651f, [C]=1.951f, [O]=2.051f, [P]=2.551f, [S]=2.651f },
    [Cl] = { [H]=1.651f, [C]=2.151f, [N]=1.951f, [O]=1.951f, [Cl]=1.951f },
    [K]  = { [O]=2.851f, [H]=2.851f, [N]=2.851f },
    [Ca] = { [O]=2.451f, [S]=2.551f },

    // Transition metals
    [Sc] = { [O]=2.151f, [N]=2.151f, [C]=2.251f },
    [Ti] = { [O]=2.051f, [N]=2.151f, [C]=2.251f },
    [V]  = { [O]=2.051f, [N]=2.151f, [C]=2.251f },
    [Cr] = { [O]=2.051f, [N]=2.151f, [C]=2.251f, [S]=2.351f },
    [Mn] = { [O]=2.051f, [N]=2.151f, [C]=2.251f, [S]=2.351f },
    [Fe] = { [O]=2.151f, [N]=2.151f, [C]=2.251f, [S]=2.251f },
    [Co] = { [O]=2.151f, [N]=2.151f, [C]=2.251f },
    [Ni] = { [O]=2.151f, [N]=2.151f, [C]=2.251f },
    [Cu] = { [O]=2.051f, [N]=2.051f, [S]=2.151f },
    [Zn] = { [O]=2.051f, [N]=2.051f, [S]=2.151f },
    [Ga] = { [O]=1.951f, [N]=1.951f, [C]=2.151f, [Ga]=2.051f },
    [Ge] = { [O]=2.051f, [C]=2.151f, [Ge]=2.151f },
    [As] = { [H]=1.651f, [C]=2.051f, [O]=2.051f, [S]=2.351f, [P]=2.151f },
    [Se] = { [H]=1.651f, [C]=2.151f, [O]=2.151f, [S]=2.451f },
    [Br] = { [H]=1.951f, [C]=2.351f, [Cl]=2.351f, [Br]=2.051f },
    [Rb] = { [O]=2.851f },
    [Sr] = { [O]=2.551f },
    [Y]  = { [O]=2.151f },
    [Zr] = { [O]=2.151f },
    [Nb] = { [O]=2.151f },
    [Mo] = { [O]=2.051f },
    [Tc] = { [O]=2.051f },
    [Ru] = { [O]=2.051f },
    [Rh] = { [O]=2.051f },
    [Pd] = { [O]=2.051f },
    [Ag] = { [O]=2.151f },
    [Cd] = { [O]=2.151f },
    [In] = { [O]=2.051f },
    [Sn] = { [O]=2.051f },
    [Sb] = { [O]=2.151f },
    [Te] = { [O]=2.151f },
    [I]  = { [H]=2.151f, [C]=2.651f, [Cl]=2.651f, [I]=2.251f },

    // Lanthanides
    [La] = { [O]=2.451f, [N]=2.451f },
    [Ce] = { [O]=2.451f, [N]=2.451f },
    [Pr] = { [O]=2.451f, [N]=2.451f },
    [Nd] = { [O]=2.451f, [N]=2.451f },
    [Sm] = { [O]=2.451f },
    [Eu] = { [O]=2.451f },
    [Gd] = { [O]=2.451f },
    [Tb] = { [O]=2.451f },
    [Dy] = { [O]=2.451f },
    [Ho] = { [O]=2.451f },
    [Er] = { [O]=2.451f },
    [Tm] = { [O]=2.451f },
    [Yb] = { [O]=2.451f },
    [Lu] = { [O]=2.451f },

    // Actinides
    [Ac] = { [O]=2.451f },
    [Th] = { [O]=2.451f },
    [Pa] = { [O]=2.451f },
    [U]  = { [O]=2.451f },
};

// ---------------------- Coordination Bonds ----------------------
static const float coord_r_min[Num_Elements][Num_Elements] = {
    [0] = {0},  // dummy

    // Common donor atoms
    [C]  = { [H]=0.90f, [C]=1.20f, [N]=1.10f, [O]=1.10f, [S]=1.50f, [P]=1.50f },
    [N]  = { [H]=0.90f, [C]=1.10f, [N]=1.20f, [O]=1.20f, [P]=1.50f, [S]=1.50f },
    [O]  = { [H]=0.90f, [C]=1.10f, [N]=1.20f, [O]=1.20f, [S]=1.50f, [P]=1.50f },
    [S]  = { [H]=1.00f, [C]=1.50f, [N]=1.50f, [O]=1.50f, [S]=1.60f, [P]=1.60f },
    [P]  = { [H]=1.10f, [C]=1.50f, [N]=1.50f, [O]=1.50f, [S]=1.60f, [P]=1.60f },

    // Halogens as donors are uncommon but possible
    [F]  = { [H]=0.90f, [C]=1.20f, [N]=1.20f, [O]=1.20f },
    [Cl] = { [H]=1.20f, [C]=1.80f, [N]=1.80f, [O]=1.80f },
    [Br] = { [H]=1.20f, [C]=2.00f, [N]=2.00f, [O]=2.00f },
    [I]  = { [H]=1.20f, [C]=2.20f, [N]=2.20f, [O]=2.20f },

    // Metals (common acceptors in coordination)
    [Li] = { [O]=1.80f, [N]=1.80f, [S]=2.00f },
    [Be] = { [O]=1.60f, [N]=1.60f, [C]=1.70f },
    [Na] = { [O]=2.00f, [N]=2.00f },
    [Mg] = { [O]=1.90f, [N]=1.90f },
    [Al] = { [O]=1.80f, [N]=1.80f },
    [K]  = { [O]=2.20f, [N]=2.20f },
    [Ca] = { [O]=2.00f, [N]=2.00f },

    // Transition metals (typical coordination distances, minimal)
    [Ti] = { [O]=1.80f, [N]=1.85f, [S]=1.90f },
    [V]  = { [O]=1.80f, [N]=1.85f },
    [Cr] = { [O]=1.85f, [N]=1.90f },
    [Mn] = { [O]=1.85f, [N]=1.90f },
    [Fe] = { [O]=1.85f, [N]=1.90f },
    [Co] = { [O]=1.85f, [N]=1.90f },
    [Ni] = { [O]=1.85f, [N]=1.90f },
    [Cu] = { [O]=1.85f, [N]=1.90f, [S]=2.00f },
    [Zn] = { [O]=1.85f, [N]=1.90f, [S]=2.00f },

    // Post-transition metals
    [Ga] = { [O]=1.80f, [N]=1.80f },
    [In] = { [O]=1.90f, [N]=1.90f },
    [Sn] = { [O]=1.90f, [N]=1.90f },

    // Selected lanthanides (typical coordination distances)
    [La] = { [O]=2.20f, [N]=2.20f },
    [Ce] = { [O]=2.20f, [N]=2.20f },
    [Nd] = { [O]=2.20f, [N]=2.20f },
    [Eu] = { [O]=2.20f, [N]=2.20f },
    [Yb] = { [O]=2.20f, [N]=2.20f },

    // Actinides
    [Th] = { [O]=2.30f, [N]=2.30f },
    [U]  = { [O]=2.30f, [N]=2.30f },
};

// Coordination max distances
static const float coord_r_max[Num_Elements][Num_Elements] = {
    [0] = {0},  // dummy

    // Common donor atoms
    [C]  = { [H]=1.40f, [C]=1.90f, [N]=1.80f, [O]=1.80f, [S]=2.20f, [P]=2.20f },
    [N]  = { [H]=1.40f, [C]=1.80f, [N]=1.90f, [O]=1.90f, [P]=2.20f, [S]=2.20f },
    [O]  = { [H]=1.40f, [C]=1.80f, [N]=1.90f, [O]=1.90f, [S]=2.20f, [P]=2.20f },
    [S]  = { [H]=1.50f, [C]=2.20f, [N]=2.20f, [O]=2.20f, [S]=2.30f, [P]=2.30f },
    [P]  = { [H]=1.60f, [C]=2.20f, [N]=2.20f, [O]=2.20f, [S]=2.30f, [P]=2.30f },

    // Halogens as donors
    [F]  = { [H]=1.40f, [C]=1.90f, [N]=1.90f, [O]=1.90f },
    [Cl] = { [H]=1.60f, [C]=2.50f, [N]=2.50f, [O]=2.50f },
    [Br] = { [H]=1.60f, [C]=2.70f, [N]=2.70f, [O]=2.70f },
    [I]  = { [H]=1.60f, [C]=3.00f, [N]=3.00f, [O]=3.00f },

    // Metals (typical acceptors in coordination)
    [Li] = { [O]=2.50f, [N]=2.50f, [S]=2.80f },
    [Be] = { [O]=2.10f, [N]=2.10f, [C]=2.20f },
    [Na] = { [O]=2.80f, [N]=2.80f },
    [Mg] = { [O]=2.50f, [N]=2.50f },
    [Al] = { [O]=2.30f, [N]=2.30f },
    [K]  = { [O]=3.00f, [N]=3.00f },
    [Ca] = { [O]=2.80f, [N]=2.80f },

    // Transition metals
    [Ti] = { [O]=2.30f, [N]=2.40f, [S]=2.50f },
    [V]  = { [O]=2.30f, [N]=2.40f },
    [Cr] = { [O]=2.40f, [N]=2.50f },
    [Mn] = { [O]=2.40f, [N]=2.50f },
    [Fe] = { [O]=2.40f, [N]=2.50f },
    [Co] = { [O]=2.40f, [N]=2.50f },
    [Ni] = { [O]=2.40f, [N]=2.50f },
    [Cu] = { [O]=2.40f, [N]=2.50f, [S]=2.60f },
    [Zn] = { [O]=2.40f, [N]=2.50f, [S]=2.60f },

    // Post-transition metals
    [Ga] = { [O]=2.30f, [N]=2.30f },
    [In] = { [O]=2.50f, [N]=2.50f },
    [Sn] = { [O]=2.50f, [N]=2.50f },

    // Lanthanides
    [La] = { [O]=2.80f, [N]=2.80f },
    [Ce] = { [O]=2.80f, [N]=2.80f },
    [Nd] = { [O]=2.80f, [N]=2.80f },
    [Eu] = { [O]=2.80f, [N]=2.80f },
    [Yb] = { [O]=2.80f, [N]=2.80f },

    // Actinides
    [Th] = { [O]=3.00f, [N]=3.00f },
    [U]  = { [O]=3.00f, [N]=3.00f },
};

static size_t find_bonds_in_ranges(md_bond_data_t* bond, const float* x, const float* y, const float* z, const md_element_t* element, const md_unit_cell_t* cell, md_range_t range_a, md_range_t range_b, md_allocator_i* alloc, md_allocator_i* temp_arena) {
    ASSERT(bond);
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(cell);

    size_t pre_count = bond->count;

    const int N = (range_a.end - range_a.beg) * (range_b.end - range_b.beg);
    // Swap ranges if required
    if (range_b.beg < range_a.beg) {
        md_range_t tmp = range_a;
        range_a = range_b;
        range_b = tmp;
    }
    if (N < 20000) {
        // Brute force
        for (int i = range_a.beg; i < range_a.end; ++i) {
            const md_256  ci[3] = {
                md_mm256_set1_ps(x[i]),
                md_mm256_set1_ps(y[i]),
                md_mm256_set1_ps(z[i]),
            };

            const float* table_cov_r_min   = cov_r_min[element[i]];
            const float* table_cov_r_max   = cov_r_max[element[i]];
            //const float* table_coord_r_min = coord_r_min[element[i]];
            //const float* table_coord_r_max = coord_r_max[element[i]];

            for (int j = MAX(i+1, range_b.beg); j < range_b.end; j += 8) {
                const md_256  cj[3] = {
                    md_mm256_loadu_ps(x + j),
                    md_mm256_loadu_ps(y + j),
                    md_mm256_loadu_ps(z + j),
                };
                const md_256i vj = md_mm256_add_epi32(md_mm256_set1_epi32(j), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
                const md_256i nj = md_mm256_cmpgt_epi32(md_mm256_set1_epi32(range_b.end), vj);
                const md_256i ej = md_mm256_and_epi32(md_mm256_cvtepu8_epi32(md_mm_loadu_epi32(element + j)), nj);

                const md_256  r_min_cov   = md_mm256_i32gather_ps(table_cov_r_min,   ej, 4);
                const md_256  r_max_cov   = md_mm256_i32gather_ps(table_cov_r_max,   ej, 4);
                //const md_256  r_min_coord = md_mm256_i32gather_ps(table_coord_r_min, ej, 4);
                //const md_256  r_max_coord = md_mm256_i32gather_ps(table_coord_r_max, ej, 4);

                const md_256 dx[3] = {
                    md_mm256_sub_ps(ci[0], cj[0]),
                    md_mm256_sub_ps(ci[1], cj[1]),
                    md_mm256_sub_ps(ci[2], cj[2]),
                };
                const md_256 d2 = simd_distance_squared(dx, cell);
                int cov_mask    = simd_bond_heuristic(d2, r_min_cov,   r_max_cov);
                //int coord_mask  = simd_bond_heuristic(d2, r_min_coord, r_max_coord);

                int mask = cov_mask;

                md_array_ensure(bond->pairs, md_array_size(bond->pairs) + popcnt32(mask), alloc);
                while (mask) {
                    const int idx = ctz32(mask);
                    mask = mask & ~(1 << idx);
                    md_array_push_no_grow(bond->pairs, ((md_atom_pair_t){i, j + idx}));
                    bond->count += 1;
                }
            }
        }
    }
    else {
        size_t temp_pos = md_vm_arena_get_pos(temp_arena);
        int capacity = range_b.end - range_b.beg;
        int* indices = md_vm_arena_push(temp_arena, sizeof(int) * capacity);

        int size = 0;
        for (int i = range_b.beg; i < range_b.end; ++i) {
            if (element[i] != 0) {
                indices[size++] = i;
            }
        }

        // Spatial acceleration structure
        md_spatial_hash_t* sh = md_spatial_hash_create_soa(x, y, z, indices, size, cell, temp_arena);

        const float cutoff = 3.0f;

        for (int i = range_a.beg; i < range_a.end; ++i) {
            if (element[i] != 0) {
                const vec4_t ci = {x[i], y[i], z[i], 0};
                const float* table_cov_r_min   = cov_r_min[element[i]];
                const float* table_cov_r_max   = cov_r_max[element[i]];
                const float* table_coord_r_min = coord_r_min[element[i]];
                const float* table_coord_r_max = coord_r_max[element[i]];

                const size_t num_indices = md_spatial_hash_query_idx(indices, capacity, sh, vec3_from_vec4(ci), cutoff);
                for (size_t idx = 0; idx < num_indices; ++idx) {
                    const int j = indices[idx];
                    // Only store monotonic bond connections
                    if (j < i) continue;

                    const float r_min_cov   = table_cov_r_min[element[j]];
                    const float r_max_cov   = table_cov_r_max[element[j]];
                    const float r_min_coord = table_coord_r_min[element[j]];
                    const float r_max_coord = table_coord_r_max[element[j]];

                    const vec4_t cj = {x[j], y[j], z[j], 0};
                    const vec4_t dx = vec4_sub(ci, cj);
                    const float d2 = distance_squared(dx, cell);

                    if (bond_heuristic(d2, r_min_cov,   r_max_cov))
                    {
                        md_array_push(bond->pairs, ((md_atom_pair_t){i, j}), alloc);
                        bond->count += 1;
                    }
                }
            }
        }

        md_vm_arena_set_pos_back(temp_arena, temp_pos);
    }

    return bond->count - pre_count;
}

static inline void find_bonds(md_bond_data_t* bond, const float* x, const float* y, const float* z, const md_element_t* element, size_t count, const md_unit_cell_t* cell) {

}

typedef struct aabb_t {
    vec3_t min_box;
    vec3_t max_box;
} aabb_t;

static inline bool aabb_overlap(aabb_t a, aabb_t b) {
    return
        (a.min_box.x <= b.max_box.x && a.max_box.x >= b.min_box.x) &&
        (a.min_box.y <= b.max_box.y && a.max_box.y >= b.min_box.y) &&
        (a.min_box.z <= b.max_box.z && a.max_box.z >= b.min_box.z);
}

void md_util_covalent_bonds_compute_exp(md_bond_data_t* bond, const float* x, const float* y, const float* z, const md_element_t* elem, size_t atom_count, const md_residue_data_t* res, const md_unit_cell_t* cell, md_allocator_i* alloc) {
    ASSERT(bond);
    ASSERT(alloc);

    md_bond_data_clear(bond);

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));
    
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing atom field (x/y/z)");
        goto done;
    }

    if (!elem) {
        MD_LOG_ERROR("Missing atom field element");
        goto done;
    }

    int* atom_res_idx = 0;
       
    if (res && res->count > 0) {
        atom_res_idx = md_vm_arena_push(temp_arena, atom_count * sizeof(int));

        // The padding applied to residue AABBs in order to determine potential overlap
        const float aabb_pad = 1.5f;
        size_t prev_idx = 0;
        md_range_t prev_range = md_residue_atom_range(*res, prev_idx);
        md_flags_t prev_flags = res->flags[prev_idx];

        for (int i = prev_range.beg; i < prev_range.end; ++i) {
            atom_res_idx[i] = prev_idx;
        }

        aabb_t prev_aabb = {0};
        md_util_aabb_compute(prev_aabb.min_box.elem, prev_aabb.max_box.elem, x + prev_range.beg, y + prev_range.beg, z + prev_range.beg, 0, 0, prev_range.end - prev_range.beg);
        prev_aabb.min_box = vec3_sub_f(prev_aabb.min_box, aabb_pad);
        prev_aabb.max_box = vec3_add_f(prev_aabb.max_box, aabb_pad);

        find_bonds_in_ranges(bond, x, y, z, elem, cell, prev_range, prev_range, alloc, temp_arena);
        for (size_t curr_idx = 1; curr_idx < res->count; ++curr_idx) {
            md_range_t curr_range = md_residue_atom_range(*res, curr_idx);
            md_flags_t curr_flags = res->flags[curr_idx];

            for (int i = curr_range.beg; i < curr_range.end; ++i) {
                atom_res_idx[i] = curr_idx;
            }

            aabb_t curr_aabb = {0};
            if (!(curr_flags & (MD_FLAG_WATER | MD_FLAG_ION))) {
                md_util_aabb_compute(curr_aabb.min_box.elem, curr_aabb.max_box.elem, x + curr_range.beg, y + curr_range.beg, z + curr_range.beg, 0, 0, curr_range.end - curr_range.beg);
                curr_aabb.min_box = vec3_sub_f(curr_aabb.min_box, aabb_pad);
                curr_aabb.max_box = vec3_add_f(curr_aabb.max_box, aabb_pad);

                // @NOTE: Interresidual bonds
                if (aabb_overlap(prev_aabb, curr_aabb)) {
                    find_bonds_in_ranges(bond, x, y, z, elem, cell, prev_range, curr_range, alloc, temp_arena);
                    // We want to flag these bonds with INTER flag to signify that they connect residues (which are used to identify chains)
                }
            }
            find_bonds_in_ranges(bond, x, y, z, elem, cell, curr_range, curr_range, alloc, temp_arena);
            prev_range = curr_range;
            prev_flags = curr_flags;
            prev_idx   = curr_idx;
            prev_aabb  = curr_aabb;
        }
    }
    else {
        md_range_t range = {0, (int)atom_count};
        find_bonds_in_ranges(bond, x, y, z, elem, cell, range, range, alloc, temp_arena);
    }

    // Compute connectivity count
    uint8_t* c_count = md_vm_arena_push_zero_array(temp_arena, uint8_t, atom_count);
    for (size_t i = 0; i < bond->count; ++i) {
        c_count[bond->pairs[i].idx[0]] += 1;
        c_count[bond->pairs[i].idx[1]] += 1;
    }

    // This is allocated in temp until we know that all connections are final
    md_conn_data_t conn = {0};
    compute_connectivity(&conn, bond->pairs, bond->count, atom_count, temp_arena, temp_arena);

    if (conn.count == 0) {
        goto done;
    }

    // We cannot remove bonds inplace as they will mess up the indexing
    md_array(uint32_t) bond_indices_to_remove = 0;

    typedef struct {
        float len2;
        int   idx;
    } bond_t;

    md_array(bond_t) bond_buf = 0;
    md_array(uint64_t) checked_bonds = make_bitfield(bond->count, temp_arena);

    // Prune over-connected atoms by removing longest bonds
    for (size_t i = 0; i < atom_count; ++i) {
        md_element_t e = elem[i];
        const int max_con = max_neighbors_element(e);

        if ((int)c_count[i] > max_con) {
            const vec3_t xi = {x[i], y[i], z[i]};
            const uint32_t conn_beg = conn.offset[i];
            const uint32_t conn_end = conn.offset[i+1];
            md_array_shrink(bond_buf, 0);

            for (uint32_t conn_idx = conn_beg; conn_idx < conn_end; ++conn_idx) {
                md_bond_idx_t bij = conn.bond_idx[conn_idx];
                if (bitfield_test_bit(checked_bonds, bij)) continue;

                bitfield_set_bit(checked_bonds, bij);
                md_atom_idx_t j = conn.atom_idx[conn_idx];
                const vec3_t xj = {x[j], y[j], z[j]};
                bond_t b = {
                    .len2 = vec3_distance_squared(xi, xj),
                    .idx  = bij,
                };
                md_array_push(bond_buf, b, temp_arena);
            }

            while ((int)md_array_size(bond_buf) > max_con) {
                size_t max_k = 0;
                for (size_t k = 1; k < md_array_size(bond_buf); ++k) {
                    if (bond_buf[k].len2 > bond_buf[max_k].len2) {
                        max_k = k;
                    }
                }
                md_array_push(bond_indices_to_remove, bond_buf[max_k].idx, temp_arena);
                md_array_swap_back_and_pop(bond_buf, max_k);
            }
        }
    }

    if (bond_indices_to_remove) {
        size_t count = md_array_size(bond_indices_to_remove);
        for (size_t i = 0; i < count; ++i) {
            uint32_t bond_idx = md_array_back(bond_indices_to_remove);
            md_array_pop(bond_indices_to_remove);
            md_array_swap_back_and_pop(bond->pairs, bond_idx);
        }
        bond->count -= count;

        // Recompute connectivity information (since we removed bonds)
        compute_connectivity(&bond->conn, bond->pairs, bond->count, atom_count, alloc, temp_arena);
    } else {
        // Commit the connectivity to molecule
        md_array_push_array(bond->conn.offset,   conn.offset,   conn.offset_count, alloc);
        md_array_push_array(bond->conn.atom_idx, conn.atom_idx, conn.count, alloc);
        md_array_push_array(bond->conn.bond_idx, conn.bond_idx, conn.count, alloc);
        bond->conn.count = conn.count;
        bond->conn.offset_count = conn.offset_count;
    }

    md_array_resize(bond->order, bond->count, alloc);

    if (res && res->count > 0) {
        // Mark interresidual bonds
        for (size_t i = 0; i < bond->count; ++i) {
            md_residue_idx_t res_idx[2] = {
                atom_res_idx[bond->pairs[i].idx[0]],
                atom_res_idx[bond->pairs[i].idx[1]],
            };
            bond->order[i] = (res_idx[0] == res_idx[1]) ? 1 : MD_BOND_FLAG_INTER;
        }
    } else {
        MEMSET(bond->order, 1, md_array_bytes(bond->order));
    }
    
done:
    md_vm_arena_destroy(temp_arena);
}

// Compute the prerequisite fields to enable hydrogen bond determination
// ported from molstar
// https://github.com/molstar/molstar/blob/master/src/mol-model-props/computed/chemistry/valence-model.ts

/*
md_array(md_hbond_data_t) md_util_compute_hbond_data(const md_molecule_t* mol, md_index_data_t connectivity, md_allocator_i* alloc) {
    md_array(md_valence_t) atom_valence = md_array_create(md_valence_t, mol->atom.count, md_get_temp_allocator());
    memset(atom_valence, 0, md_array_bytes(atom_valence));

    md_array(int8_t) atom_hcount = md_array_create(int8_t, mol->atom.count, md_get_temp_allocator());
    memset(atom_hcount, 0, md_array_bytes(atom_hcount));

    md_array(md_order_t) order = md_util_compute_covalent_bond_order(mol->bonds, md_array_size(mol->bonds), mol->atom.name, mol->atom.residue_idx, mol->residue.name, md_get_temp_allocator());

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

#define MIN_RES_LEN 4
#define MAX_RES_LEN 25

#define MIN_NUC_LEN 6
#define MAX_NUC_LEN 35

#define MIN_RESIDUE_CHAIN_LEN 2

bool md_util_compute_residue_data(md_residue_data_t* res, md_atom_data_t* atom, md_allocator_i* alloc) {
    ASSERT(res);
    ASSERT(atom);
    ASSERT(alloc);

    md_residue_id_t prev_resid = -1;
    str_t prev_resstr = STR_LIT("");

    for (size_t i = 0; i < atom->count; ++i) {
        const md_residue_id_t resid = atom->resid[i];
        const md_flags_t flags = atom->flags[i];
        const str_t resstr = atom->resname ? LBL_TO_STR(atom->resname[i]) : STR_LIT("");

        if (resid != prev_resid || !str_eq(resstr, prev_resstr) || flags & (MD_FLAG_SEQ_TERM)) {
            md_array_push(res->id, resid, alloc);
            md_array_push(res->name, make_label(resstr), alloc);
            md_array_push(res->atom_offset, (uint32_t)i, alloc);
            md_array_push(res->flags, 0, alloc);
            res->count += 1;
        }

        prev_resstr = resstr;
        prev_resid = resid;
    }

    md_array_push(res->atom_offset, (uint32_t)atom->count, alloc);

    atom->res_idx = md_array_create(md_residue_idx_t, atom->count, alloc);
    MEMSET(atom->res_idx, -1, md_array_bytes(atom->res_idx));

    for (size_t i = 0; i < res->count; ++i) {
        str_t resname = LBL_TO_STR(res->name[i]);
        md_range_t range = md_residue_atom_range(*res, i);
        size_t len = (size_t)(range.end - range.beg);

        if (i == 130) {
            while(0) {};
        }

        md_protein_backbone_atoms_t prot_atoms = {0};
        md_nucleic_backbone_atoms_t nucl_atoms = {0};
        if (MIN_RES_LEN <= len && len <= MAX_RES_LEN && (atom->type && md_util_protein_backbone_atoms_extract(&prot_atoms, atom->type + range.beg, len, range.beg)))
        {
            res->flags[i] |= MD_FLAG_AMINO_ACID;
            atom->flags[prot_atoms.n]  |= MD_FLAG_BACKBONE;
            atom->flags[prot_atoms.ca] |= MD_FLAG_BACKBONE;
            atom->flags[prot_atoms.c]  |= MD_FLAG_BACKBONE;
            //atom->flags[prot_atoms.o]  |= MD_FLAG_BACKBONE;
        } else if (MIN_NUC_LEN <= len && len <= MAX_NUC_LEN && (atom->type && md_util_nucleic_backbone_atoms_extract(&nucl_atoms, atom->type + range.beg, len, range.beg))) {
            res->flags[i] |= MD_FLAG_NUCLEOTIDE;
            atom->flags[nucl_atoms.c5] |= MD_FLAG_BACKBONE;
            atom->flags[nucl_atoms.c4] |= MD_FLAG_BACKBONE;
            atom->flags[nucl_atoms.c3] |= MD_FLAG_BACKBONE;
            atom->flags[nucl_atoms.o3] |= MD_FLAG_BACKBONE;
            atom->flags[nucl_atoms.p]  |= MD_FLAG_BACKBONE;
            atom->flags[nucl_atoms.o5] |= MD_FLAG_BACKBONE;
        } else if ((len == 1 || len == 3) && (md_util_resname_water(resname) || (atom->type && len == 1 && md_util_resname_water(LBL_TO_STR(atom->type[range.beg]))))) {
            res->flags[i] |= MD_FLAG_WATER;
        } else if (md_util_resname_amino_acid(resname)) {
            res->flags[i] |= MD_FLAG_AMINO_ACID;
        }
#if 0
        else if (len == 1) {
            res->flags[i] |= MD_FLAG_ION;
        }
#endif

        for (int j = range.beg; j < range.end; ++j) {
            atom->res_idx[j] = (md_residue_idx_t)i;

            // Propagate flags to atoms
            atom->flags[j] |= res->flags[i];

            res->flags[i] |= atom->flags[j] & MD_FLAG_SEQ_TERM;
        }
    }

    return true;
}

// @TODO: Convert to a table
static bool monatomic_ion_element(md_element_t elem) {
    switch (elem) {
    case 1:   // H+ / H-
    case 3:   // Li+
    case 9:   // F-
    case 11:  // Na+
    case 17:  // Cl-
    case 19:  // K+
    case 35:  // Br-
    case 37:  // Rb+
    case 47:  // Ag+
    case 53:  // I-
    case 55:  // Cs+
    case 4:   // Be2+
    case 8:   // O2-
    case 12:  // Mg2+
    case 20:  // Ca2+
    case 16:  // S2-
    case 30:  // Zn2+
    case 34:  // Se2-
    case 38:  // Sr2+
    case 48:  // Cd2+
    case 56:  // Ba2+
    case 66:  // Hg2+
    case 7:   // N3-
    case 13:  // Al3+
    case 15:  // P3-
    case 29:  // Cu1+ / Cu2+
    case 24:  // Cr2+ / Cr3+
    case 25:  // Mn2+ / Mn3+
    case 26:  // Fe2+ / Fe3+
    case 27:  // Co2+ / Co3+
    case 50:  // Sn2+ / Sn4+
    case 68:  // Pb2+ / Pb4+
        return true;
    default:
        return false;
    }
}

bool md_util_identify_ions(md_atom_data_t* atom, const md_bond_data_t* bond) {
    ASSERT(atom);
    ASSERT(bond);
    if (!atom->element || !bond->conn.offset) {
        return false;
    }
    for (size_t i = 0; i < atom->count; ++i) {
        // Check if it has no bonds
        if (md_bond_conn_count(*bond, i) == 0 && monatomic_ion_element(atom->element[i]) && !(atom->flags[i] & MD_FLAG_WATER)) {
            atom->flags[i] |= MD_FLAG_ION;
        }
    }
    return true;
}

void md_util_init_hydrogen_bond_data(md_hydrogen_bond_data_t* hbond_data, md_atom_data_t* atom_data, const md_bond_data_t* bond_data, md_allocator_i* alloc) {
    ASSERT(hbond_data);
    ASSERT(atom_data);
    ASSERT(bond_data);

    if (!atom_data->element || !bond_data->conn.offset) {
        return;
    }

    hbond_data->num_donors = 0;
    md_array_shrink(hbond_data->donors, 0);

    hbond_data->num_acceptors = 0;
    md_array_shrink(hbond_data->acceptors, 0);

    // Identify donors and acceptors
    for (size_t i = 0; i < atom_data->count; ++i) {
        int max_conn = 0;
        int num_of_lone_pairs = 2;
        switch (atom_data->element[i]) {
        case N:
        case S:
            max_conn = 3;
            break;
        case O:
            max_conn = 2;
            break;
        default:
            continue;
        }

        md_flags_t flags = 0;
        
        md_bond_iter_t it = md_bond_iter(bond_data, i);
        while (md_bond_iter_has_next(it)) {
            md_atom_idx_t j = md_bond_iter_atom_index(it);
            if (atom_data->element[j] == H) {
                flags |= MD_FLAG_H_DONOR;
                md_hydrogen_bond_donor_t donor = {(md_atom_idx_t)i, j};
                md_array_push(hbond_data->donors, donor, alloc);
            }
            md_bond_iter_next(&it);
        }
        
        size_t num_conn = md_bond_conn_count(*bond_data, i);
        if (num_conn <= max_conn) {
            flags |= MD_FLAG_H_ACCEPTOR;

            if (atom_data->element[i] == S) {
                num_of_lone_pairs = 4 - num_conn;
            }

            md_hydrogen_bond_acceptor_t acceptor = {(md_atom_idx_t)i, num_of_lone_pairs};
            md_array_push(hbond_data->acceptors, acceptor, alloc);
            hbond_data->num_acceptors += 1;
        }

        atom_data->flags[i] |= flags;
    }

    md_array_ensure(hbond_data->bonds, 2 * hbond_data->num_donors, alloc);
}

// @NOTE(Robin): This could certainly be improved to incorporate more characters
// Perhaps first A-Z, then [A-Z]0-9, then AA-ZZ etc.
static inline md_label_t generate_chain_id_from_index(size_t idx) {
    char c = 'A' + (idx % 26);
    str_t str = {&c, 1};
    return make_label(str);
}

bool md_util_compute_chain_data(md_chain_data_t* chain, md_atom_data_t* atom, const md_residue_data_t* res, const md_bond_data_t* bond, md_allocator_i* alloc) {
    ASSERT(chain);
    ASSERT(alloc);
    
    if (!atom) {
        MD_LOG_ERROR("atom data is null");
        return false;
    }

    if (!bond) {
        MD_LOG_ERROR("bond data is null");
        return false;
    }

    if (res->count == 0) {
        MD_LOG_DEBUG("Dataset only contains single residue, no artificial chains can be created");
        return false;
    }

    MEMSET(chain, 0, sizeof(md_chain_data_t));

    if (res->count <= 1) {
        return true;
    }

    md_array(uint64_t) res_bond_to_prev = make_bitfield(res->count + 1, md_get_temp_allocator());
    md_residue_idx_t max_bonded_residue_idx = 0;
    for (size_t i = 0; i < bond->count; ++i) {
        if (bond->order[i] & MD_BOND_FLAG_INTER) {
            const md_residue_idx_t res_a = atom->res_idx[bond->pairs[i].idx[0]];
            const md_residue_idx_t res_b = atom->res_idx[bond->pairs[i].idx[1]];
            const md_residue_idx_t res_max = MAX(res_a, res_b);
            if (abs(res_b - res_a) == 1) {
                bitfield_set_bit(res_bond_to_prev, res_max);
                max_bonded_residue_idx = MAX(max_bonded_residue_idx, res_max);
            }
        }
    }

    int beg_idx = 0;
    str_t prev_id = {0};
    md_flags_t prev_flags = 0;
    int chain_end_idx = 0;

    for (int i = 0; i <= max_bonded_residue_idx + 1; ++i) {
        str_t id = {0};
        md_flags_t flags = 0;
        if (i < res->count) {
            const md_range_t atom_range = md_residue_atom_range(*res, i);
            id = atom->chainid ? LBL_TO_STR(atom->chainid[atom_range.beg]) : (str_t){0};
            flags = res->flags[i];
        }

        if (i < max_bonded_residue_idx && !(flags & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE))) {
            beg_idx += 1;
            goto next;
        }

        bool different_id = !str_empty(id) && !str_eq(id, prev_id);
        bool terminal_res = flags & MD_FLAG_SEQ_TERM;
        bool disconnected = !bitfield_test_bit(res_bond_to_prev, i);

        if (different_id || terminal_res || disconnected) {
            int end_idx = i;
            if (end_idx - beg_idx >= MIN_RESIDUE_CHAIN_LEN) {
                md_label_t lbl = str_empty(prev_id) ? generate_chain_id_from_index(chain->count) : make_label(prev_id);
                md_range_t res_rng = {beg_idx, end_idx};
                md_range_t atom_rng = {res->atom_offset[beg_idx], res->atom_offset[end_idx]};

                md_array_push(chain->id, lbl, alloc);
                md_array_push(chain->atom_range, atom_rng, alloc);
                md_array_push(chain->res_range,   res_rng, alloc);

                chain_end_idx = i;
                chain->count += 1;
            }
            beg_idx = i;
        }

next:
        prev_id = id;
        prev_flags = flags;
    }

    md_array_resize(atom->chain_idx, atom->count, alloc);
    MEMSET(atom->chain_idx, -1, md_array_bytes(atom->chain_idx));

    for (size_t i = 0; i < chain->count; ++i) {
        const md_range_t atom_range = md_chain_atom_range(*chain, i);
        for (int j = atom_range.beg; j < atom_range.end; ++j) {
            atom->chain_idx[j] = (md_chain_idx_t)i;
        }
    }

    return true;
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
    const int64_t num_rings = md_index_data_num_ranges(ring_data);
    for (int64_t i = 0; i < num_rings; ++i) {
        const int* ring_i     = md_index_range_beg(ring_data, i);
        const int ring_i_size = (int)md_index_range_size(ring_data, i);
        if (compare_ring(ring_i, ring_i_size, ring, ring_size)) return true;
    }
    return false;
}

#define SORT_RING_PRINT 0

// This is a specialized deterministic sort function for reordering the indices of a ring in a particular way such that it can be subjected to MEMCMP
static void sort_ring(int* arr, int n) {
    ASSERT(n >= 0);
    if (n == 0) return;

#if SORT_RING_PRINT
    printf("pre:  [");
    for (int i = 0; i < n; ++i) {
        char delim = i < (n-1) ? ' ' : ']';
        printf("%i%c", arr[i], delim);
    }
    printf("\n");
#endif

    int lowest_idx = 0;
    for (int i = 1; i < n; ++i) {
        if (arr[i] < arr[lowest_idx]) {
            lowest_idx = i;
        }
    }

    // Rotate ring such that the lowest index in the first slot
    if (lowest_idx != 0) {
        int* buf = md_temp_push(sizeof(*arr) * n);
        int j = 0;
        for (int i = lowest_idx; i < n; ++i) {
            buf[j++] = arr[i];
        }
        for (int i = 0; i < lowest_idx; ++i) {
            buf[j++] = arr[i];
        }
        MEMCPY(arr, buf, sizeof(*arr) * n);
    }

#if SORT_RING_PRINT
    printf("mid:  [");
    for (int i = 0; i < n; ++i) {
        char delim = i < (n-1) ? ' ' : ']';
        printf("%i%c", arr[i], delim);
    }
    printf("\n");
#endif

    // And we want to store the lowest neighboring index in the second slot (right in this cyclic array)
    if (arr[n - 1] < arr[1]) {
        // Need to swap
        for (int i = 1, j = n - 1; i < j; ++i, --j) {
            int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }
#if SORT_RING_PRINT
    printf("post: [");
    for (int i = 0; i < n; ++i) {
        char delim = i < (n-1) ? ' ' : ']';
        printf("%i%c", arr[i], delim);
    }
    printf("\n\n");
#endif
}

#define MIN_RING_SIZE 3
#define MAX_RING_SIZE 6
#define MAX_DEPTH (MAX_RING_SIZE/2 + 1)

// Inspired by molstars implementation
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/structure/unit/rings/compute.ts#L249
// Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
// @author David Sehnal <david.sehnal@gmail.com>

size_t md_util_compute_rings(md_index_data_t* out_rings, const md_atom_data_t* atom, const md_bond_data_t* bond) {
    ASSERT(out_rings);
    ASSERT(out_rings->alloc);

    ASSERT(atom);
    ASSERT(bond);

    if (atom->count == 0) {
        return 0;
    }
    if (bond->count == 0) {
        return 0;
    }

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    // Some typedefs to enable laboration of types
    typedef uint16_t color_t;
    typedef uint8_t  depth_t;
    typedef uint16_t mark_t;
    typedef  int32_t pred_t;
        
    color_t* color = md_vm_arena_push_zero_array(temp_arena, color_t, atom->count);
    depth_t* depth = md_vm_arena_push_zero_array(temp_arena, depth_t, atom->count);
    mark_t*  mark  = md_vm_arena_push_zero_array(temp_arena, mark_t,  atom->count);
    pred_t*  pred  = md_vm_arena_push_array     (temp_arena, pred_t,  atom->count);
    MEMSET(pred, -1, atom->count * sizeof(pred_t));    // We can do memset as the representation of -1 under two's complement is 0xFFFFFFFF

    // The capacity is arbitrary here, but will be resized if needed.
    fifo_t queue = fifo_create(64, temp_arena);

    md_hashset_t ring_set = {.allocator = temp_arena};
    
    color_t current_color = 1;
    mark_t  current_mark  = 1;

#if DEBUG
    size_t processed_ring_elements = 0;
#endif

    size_t num_rings = 0;

    for (int atom_idx = 0; atom_idx < (int)atom->count; ++atom_idx) {
        if (atom->flags[atom_idx] & (MD_FLAG_WATER | MD_FLAG_ION)) continue;

        // Skip any atom that has already been colored in the previous search
        if (color[atom_idx] == current_color) continue;

        // Set i as our root
        depth[atom_idx] = 1;
        pred [atom_idx] = -1;
        mark [atom_idx] = current_mark++;

        fifo_clear(&queue);
        fifo_push(&queue, atom_idx);
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);

#if DEBUG
            processed_ring_elements += 1;
#endif

            md_bond_iter_t it = md_bond_iter(bond, idx);
            while (md_bond_iter_has_next(it)) {
                int next = md_bond_iter_atom_index(it);
                md_bond_iter_next(&it);

                //if (next > idx) continue;
                if (next == pred[idx]) continue;  // avoid adding parent to search queue

                if (mark[next] == current_mark && depth[idx] >= MIN_RING_SIZE) {
                    // We found a junction point where the graph connects
                    // Now we need to traverse both branches of this graph back up
                    // In order to find the other junction where the graph branches
                    // We follow both branches up towards the root.

                    int l = idx;
                    int r = next;

                    // Only process one of the two branches/cases
                    if (r > l) continue;

                    int len = 0;
                    int ring[MAX_DEPTH * 2];

                    color_t col = current_color++;
                    int cur;

                    cur = l;
                    for (int d = 0; d < MAX_DEPTH; ++d) {
                        color[cur] = col;
                        cur = pred[cur];
                        if (cur < 0) break;
                    }

                    int target = -1;
                    cur = r;
                    for (int d = 0; d < MAX_DEPTH; ++d) {
                        if (color[cur] == col) {
                            target = cur;
                            break;
                        }
                        ring[len++] = cur;
                        cur = pred[cur];
                        if (cur == -1) break;
                    }

                    // No junction point found
                    if (target == -1) continue;

                    cur = l;
                    int left_beg = len;
                    for (int d = 0; d < MAX_DEPTH; ++d) {
                        ring[len++] = cur;
                        if (target == cur) break;
                        cur = pred[cur];
                        if (cur == -1) break;
                    }

                    // Otherwise we made a big whoopsie
                    ASSERT(len < (int)ARRAY_SIZE(ring));

                    if (MIN_RING_SIZE <= len && len <= MAX_RING_SIZE) {
                        // Ther order of the left branch is reversed so we need to swap that
                        // Such that the indices of the ring sequentially iterates the ring
                        for (int i = left_beg, j = len - 1; i < j; ++i, --j) {
                            int tmp = ring[i];
                            ring[i] = ring[j];
                            ring[j] = tmp;
                        }

                        // @TODO: Implement a custom sort for rings such that the indices within a ring
                        // Are predictable and can be easily MEMCMP or hashed
                        // But does not violate the sequential indexing such that it can be iterated
                        sort_ring(ring, len);
                        uint64_t key = md_hash64(ring, len * sizeof(*ring), 0);
                        if (!md_hashset_get(&ring_set, key)) {
                            md_hashset_add(&ring_set, key);
                            md_index_data_push_arr(out_rings, ring, len);
                            num_rings += 1;
                        }
                    }
                } else {
                    // Avoid expanding too far from the root as we are only interested in rings up to a certain size
                    // Otherwise, we will may have alot of potential large cycles such as in the case of C60 or graphene
                    depth_t d = depth[idx] + 1;
                    if (d > MAX_DEPTH) continue;

                    depth[next] = d;
                    pred[next]  = idx;
                    mark[next]  = current_mark;
                    fifo_push(&queue, next);
                }
            }
        }
    }
    
    md_vm_arena_destroy(temp_arena);

#if 0
    MD_LOG_DEBUG("Processed ring elements: %llu\n", processed_ring_elements);
#endif

    return num_rings;
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE
#undef MAX_DEPTH

// Identifies isolated 'structures' defined by covalent bonds. Any set of atoms connected by covalent bonds are considered a structure
size_t md_util_compute_structures(md_index_data_t* out_structures, const md_bond_data_t* bond) {
    ASSERT(out_structures);
    ASSERT(out_structures->alloc);

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    const size_t atom_count = bond->conn.offset_count - 1;

    // Create a bitfield to keep track of which atoms have been visited
    uint64_t* visited = make_bitfield(atom_count, temp_arena);

    // The capacity is arbitrary here, and will be resized if needed.
    fifo_t queue = fifo_create(1024, temp_arena);

    md_array(int) indices = md_array_create(int, 1024, temp_arena);

    size_t num_structures = 0;

    for (int i = 0; i < (int)atom_count; ++i) {
        md_array_shrink(indices, 0);
        // Skip any atom which has already been touched
        if (bitfield_test_bit(visited, i)) continue;

        fifo_clear(&queue);
        fifo_push(&queue, i);
        while (!fifo_empty(&queue)) {
            int cur = fifo_pop(&queue);
            if (bitfield_test_bit(visited, cur)) continue;
            bitfield_set_bit(visited, cur);

            md_array_push(indices, cur, temp_arena);
            
            md_bond_iter_t it = md_bond_iter(bond, cur);
            while (md_bond_iter_has_next(it)) {
                int next = md_bond_iter_atom_index(it);
                if (!bitfield_test_bit(visited, next)) {
                    fifo_push(&queue, next);
                }
                md_bond_iter_next(&it);
            }
        }

        // Sort the indices within the structure for more coherent memory access
        size_t size = md_array_size(indices);
        if (size < 128) {
            sort_arr(indices, (int)md_array_size(indices));
        } else {
            sort_radix_inplace_uint32((uint32_t*)indices, md_array_size(indices), temp_arena);
        }
        
        // Here we should have exhausted every atom that is connected to index i.
        md_index_data_push_arr(out_structures, indices, md_array_size(indices));
        num_structures += 1;
    }
    
    md_vm_arena_destroy(temp_arena);

    return num_structures;
}

void md_util_mask_grow_by_bonds(md_bitfield_t* mask, const md_molecule_t* mol, size_t extent, const md_bitfield_t* viable_mask) {
    ASSERT(mask);
    ASSERT(mol);

    if (extent >= 255) {
        MD_LOG_DEBUG("Maximum supported growth extent is 255, the extent will be clamped to 255");
        extent = 255;
    }
    // The initial depth is 1, so we need to add one.
    extent += 1;
    
    const size_t mask_size = md_bitfield_popcount(mask);
    if (!mask_size) return;

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    int* indices = md_alloc(temp_alloc, mask_size * sizeof(int));
    const size_t num_indices = md_bitfield_iter_extract_indices(indices, mask_size, md_bitfield_iter_create(mask));
    ASSERT(num_indices == mask_size);

    fifo_t queue = fifo_create(64, temp_alloc);

    uint8_t* depth = md_alloc(temp_alloc, mol->atom.count * sizeof(uint8_t));
    MEMSET(depth, 0, mol->atom.count * sizeof(uint8_t));

    {
        md_bitfield_iter_t it = md_bitfield_iter_create(mask);
        while (md_bitfield_iter_next(&it)) {
            int idx = (int)md_bitfield_iter_idx(&it);
            depth[idx] = 1;
        }
    }
    
    for (size_t j = 0; j < num_indices; ++j) {
        int i = indices[j];

        fifo_clear(&queue);
        fifo_push(&queue, i);
        
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);
            md_bitfield_set_bit(mask, idx);

            md_bond_iter_t it = md_bond_iter(&mol->bond, idx);
            while (md_bond_iter_has_next(it)) {
                int next = md_bond_iter_atom_index(it);
                md_bond_iter_next(&it);
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

    md_arena_allocator_destroy(temp_alloc);
}

void md_util_mask_grow_by_radius(md_bitfield_t* mask, const struct md_molecule_t* mol, double radius, const md_bitfield_t* viable_mask) {
    ASSERT(mask);
    ASSERT(mol);
    
    if (radius <= 0.0) return;
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    
    int32_t* indices = 0;
    size_t count = mol->atom.count;
        
    if (viable_mask) {
        md_bitfield_t tmp_bf = md_bitfield_create(arena);
        md_bitfield_andnot(&tmp_bf, viable_mask, mask);
            
        const size_t num_atoms = md_bitfield_popcount(&tmp_bf);
        indices = md_vm_arena_push(arena, num_atoms * sizeof(int32_t));
        count = md_bitfield_iter_extract_indices(indices, num_atoms, md_bitfield_iter_create(&tmp_bf));
    }

    md_spatial_hash_t* ctx = md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, indices, count, &mol->unit_cell, arena);
    
    md_bitfield_t old_mask = md_bitfield_create(arena);
    md_bitfield_copy(&old_mask, mask);

    md_bitfield_iter_t it = md_bitfield_iter_create(&old_mask);
    const float rad = (float)radius;
    while (md_bitfield_iter_next(&it)) {
        int idx = (int)md_bitfield_iter_idx(&it);
        const vec3_t pos = {mol->atom.x[idx], mol->atom.y[idx], mol->atom.z[idx]};
        md_spatial_hash_query_bits(mask, ctx, pos, rad);
    }

    md_vm_arena_destroy(arena);
}

void md_util_grow_mask_by_residue(md_bitfield_t* mask, const md_molecule_t* mol) {

}

void md_util_grow_mask_by_structure(struct md_bitfield_t* mask, const uint32_t* structure_offsets, const int32_t* structure_indices, size_t num_structures) {

}

void md_util_aabb_compute(float out_aabb_min[3], float out_aabb_max[3], const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int32_t* in_idx, size_t count) {
    ASSERT(out_aabb_min);
    ASSERT(out_aabb_max);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    md_256 vx_min = md_mm256_set1_ps(+FLT_MAX);
    md_256 vy_min = md_mm256_set1_ps(+FLT_MAX);
    md_256 vz_min = md_mm256_set1_ps(+FLT_MAX);

    md_256 vx_max = md_mm256_set1_ps(-FLT_MAX);
    md_256 vy_max = md_mm256_set1_ps(-FLT_MAX);
    md_256 vz_max = md_mm256_set1_ps(-FLT_MAX);

    size_t i = 0;
    const size_t simd_elem = 8;
    const size_t simd_count = ROUND_DOWN(count, simd_elem);

    if (in_idx) {
        if (in_r) {
            for (; i < simd_count; i += simd_elem) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);

                md_256 x = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 y = md_mm256_i32gather_ps(in_y, idx, 4);
                md_256 z = md_mm256_i32gather_ps(in_z, idx, 4);
                md_256 r = md_mm256_i32gather_ps(in_r, idx, 4);

                vx_min = md_mm256_min_ps(vx_min, md_mm256_sub_ps(x, r));
                vy_min = md_mm256_min_ps(vy_min, md_mm256_sub_ps(y, r));
                vz_min = md_mm256_min_ps(vz_min, md_mm256_sub_ps(z, r));

                vx_max = md_mm256_max_ps(vx_max, md_mm256_add_ps(x, r));
                vy_max = md_mm256_max_ps(vy_max, md_mm256_add_ps(y, r));
                vz_max = md_mm256_max_ps(vz_max, md_mm256_add_ps(z, r));
            }
        } else {
            for (; i < simd_count; i += simd_elem) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);

                md_256 x = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 y = md_mm256_i32gather_ps(in_y, idx, 4);
                md_256 z = md_mm256_i32gather_ps(in_z, idx, 4);

                vx_min = md_mm256_min_ps(vx_min, x);
                vy_min = md_mm256_min_ps(vy_min, y);
                vz_min = md_mm256_min_ps(vz_min, z);

                vx_max = md_mm256_max_ps(vx_max, x);
                vy_max = md_mm256_max_ps(vy_max, y);
                vz_max = md_mm256_max_ps(vz_max, z);
            }
        }
    } else {
        if (in_r) {
            for (; i < simd_count; i += simd_elem) {
                md_256 x = md_mm256_loadu_ps(in_x + i);
                md_256 y = md_mm256_loadu_ps(in_y + i);
                md_256 z = md_mm256_loadu_ps(in_z + i);
                md_256 r = md_mm256_loadu_ps(in_r + i);

                vx_min = md_mm256_min_ps(vx_min, md_mm256_sub_ps(x, r));
                vy_min = md_mm256_min_ps(vy_min, md_mm256_sub_ps(y, r));
                vz_min = md_mm256_min_ps(vz_min, md_mm256_sub_ps(z, r));

                vx_max = md_mm256_max_ps(vx_max, md_mm256_add_ps(x, r));
                vy_max = md_mm256_max_ps(vy_max, md_mm256_add_ps(y, r));
                vz_max = md_mm256_max_ps(vz_max, md_mm256_add_ps(z, r));
            }
        } else {
            for (; i < simd_count; i += simd_elem) {
                md_256 x = md_mm256_loadu_ps(in_x + i);
                md_256 y = md_mm256_loadu_ps(in_y + i);
                md_256 z = md_mm256_loadu_ps(in_z + i);

                vx_min = md_mm256_min_ps(vx_min, x);
                vy_min = md_mm256_min_ps(vy_min, y);
                vz_min = md_mm256_min_ps(vz_min, z);

                vx_max = md_mm256_max_ps(vx_max, x);
                vy_max = md_mm256_max_ps(vy_max, y);
                vz_max = md_mm256_max_ps(vz_max, z);
            }
        }
    }

    vec4_t aabb_min = (vec4_t) { md_mm256_reduce_min_ps(vx_min), md_mm256_reduce_min_ps(vy_min), md_mm256_reduce_min_ps(vz_min) };
    vec4_t aabb_max = (vec4_t) { md_mm256_reduce_max_ps(vx_max), md_mm256_reduce_max_ps(vy_max), md_mm256_reduce_max_ps(vz_max) };

    // Handle remainder
    if (in_idx) {
        if (in_r) {
            for (; i < count; ++i) {
                int32_t idx = in_idx[i];
                const vec4_t c = { in_x[idx], in_y[idx], in_z[idx] };
                const vec4_t r = vec4_set1(in_r[idx]);
                aabb_min = vec4_min(aabb_min, vec4_sub(c, r));
                aabb_max = vec4_max(aabb_max, vec4_add(c, r));
            }
        } else {
            for (; i < count; ++i) {
                int32_t idx = in_idx[i];
                const vec4_t c = { in_x[idx], in_y[idx], in_z[idx] };
                aabb_min = vec4_min(aabb_min, c);
                aabb_max = vec4_max(aabb_max, c);
            }
        }
    } else {
        if (in_r) {
            for (; i < count; ++i) {
                const vec4_t c = { in_x[i], in_y[i], in_z[i] };
                const vec4_t r = vec4_set1(in_r[i]);
                aabb_min = vec4_min(aabb_min, vec4_sub(c, r));
                aabb_max = vec4_max(aabb_max, vec4_add(c, r));
            }
        } else {
            for (; i < count; ++i) {
                const vec4_t c = { in_x[i], in_y[i], in_z[i] };
                aabb_min = vec4_min(aabb_min, c);
                aabb_max = vec4_max(aabb_max, c);
            }
        }
    }

    MEMCPY(out_aabb_min, &aabb_min, sizeof(float) * 3);
    MEMCPY(out_aabb_max, &aabb_max, sizeof(float) * 3);
}

// Computes the minimum axis aligned bounding box for a set of points with a given radius (radius is optional), indices are used to select a subset of points
void md_util_aabb_compute_vec4(float out_aabb_min[3], float out_aabb_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count) {
    vec4_t aabb_min = vec4_set1( FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);
    if (in_idx) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx = in_idx[i];
            vec4_t xyzr = in_xyzr[idx];
            vec4_t r = vec4_splat_w(xyzr);
            aabb_min = vec4_min(aabb_min, vec4_sub(xyzr, r));
            aabb_max = vec4_max(aabb_max, vec4_add(xyzr, r));
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            vec4_t xyzr = in_xyzr[i];
            vec4_t r = vec4_splat_w(xyzr);
            aabb_min = vec4_min(aabb_min, vec4_sub(xyzr, r));
            aabb_max = vec4_max(aabb_max, vec4_add(xyzr, r));
        }
    }
    MEMCPY(out_aabb_min, &aabb_min, sizeof(float) * 3);
    MEMCPY(out_aabb_max, &aabb_max, sizeof(float) * 3);
}


void md_util_oobb_compute(float out_basis[3][3], float out_ext_min[3], float out_ext_max[3], const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int32_t* in_idx, size_t count, const md_unit_cell_t* cell) {
    ASSERT(out_basis);
    ASSERT(out_ext_min);
    ASSERT(out_ext_max);

    if (count == 0) return;

    vec3_t com = md_util_com_compute(in_x, in_y, in_z, NULL, in_idx, count, cell);
    double cov[3][3] = {0};
    if (cell) {
        vec4_t ref = vec4_from_vec3(com, 0);
        if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
            vec4_t ext = vec4_from_vec3(mat3_diag(cell->basis), 0);

            for (size_t i = 0; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                const vec4_t c = vec4_deperiodize_ortho(vec4_set(in_x[idx], in_y[idx], in_z[idx], 0.0f), ref, ext);

                cov[0][0] += c.x * c.x;
                cov[0][1] += c.x * c.y;
                cov[0][2] += c.x * c.z;
                cov[1][0] += c.y * c.x;
                cov[1][1] += c.y * c.y;
                cov[1][2] += c.y * c.z;
                cov[2][0] += c.z * c.x;
                cov[2][1] += c.z * c.y;
                cov[2][2] += c.z * c.z;
            }
        } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
            for (size_t i = 1; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                vec4_t c = vec4_set(in_x[idx], in_y[idx], in_z[idx], 0.0f);
                deperiodize_triclinic(c.elem, ref.elem, cell->basis.elem);

                cov[0][0] += c.x * c.x;
                cov[0][1] += c.x * c.y;
                cov[0][2] += c.x * c.z;
                cov[1][0] += c.y * c.x;
                cov[1][1] += c.y * c.y;
                cov[1][2] += c.y * c.z;
                cov[2][0] += c.z * c.x;
                cov[2][1] += c.z * c.y;
                cov[2][2] += c.z * c.z;
            }
        }
    }

    mat3_t cov_mat;
    const double scl = 1.0 / count;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            cov_mat.elem[i][j] = (float)(cov[i][j] * scl);
        }
    }

    mat3_eigen_t eigen = mat3_eigen(cov_mat);
    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    mat4_t Ri  = mat4_from_mat3(PCA);

    // Compute min and maximum extent along the PCA axes
    vec4_t min_ext = vec4_set1( FLT_MAX);
    vec4_t max_ext = vec4_set1(-FLT_MAX);
    mat3_t basis = mat3_transpose(PCA);

    // Transform the gto (x,y,z,radius) into the PCA frame to find the min and max extend within it
    for (size_t i = 0; i < count; ++i) {
        int32_t idx = in_idx ? in_idx[i] : i;
        float  r = in_r ? in_r[idx] : 0.0f;
        vec4_t c = { in_x[idx], in_y[idx], in_z[idx], 1.0f };

        vec4_t p = mat4_mul_vec4(Ri, c);
        min_ext = vec4_min(min_ext, vec4_sub_f(p, r));
        max_ext = vec4_max(max_ext, vec4_add_f(p, r));
    }

    MEMCPY(out_basis,   &basis,   sizeof(mat3_t));
    MEMCPY(out_ext_min, &min_ext, sizeof(vec3_t));
    MEMCPY(out_ext_max, &max_ext, sizeof(vec3_t));
}

void md_util_oobb_compute_vec4(float out_basis[3][3], float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count, const md_unit_cell_t* cell) {
    ASSERT(out_basis);
    ASSERT(out_ext_min);
    ASSERT(out_ext_max);

    if (count == 0) return;

    vec3_t com = md_util_com_compute_vec4(in_xyzr, in_idx, count, cell);
    double cov[3][3] = {0};
    if (cell) {
        vec4_t ref = vec4_from_vec3(com, 0);
        if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
            vec4_t ext = vec4_from_vec3(mat3_diag(cell->basis), 0);

            for (size_t i = 0; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                const vec4_t c = vec4_deperiodize_ortho(in_xyzr[idx], ref, ext);

                cov[0][0] += c.x * c.x;
                cov[0][1] += c.x * c.y;
                cov[0][2] += c.x * c.z;
                cov[1][0] += c.y * c.x;
                cov[1][1] += c.y * c.y;
                cov[1][2] += c.y * c.z;
                cov[2][0] += c.z * c.x;
                cov[2][1] += c.z * c.y;
                cov[2][2] += c.z * c.z;
            }
        } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
            for (size_t i = 1; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                vec4_t c = in_xyzr[idx];
                deperiodize_triclinic(c.elem, ref.elem, cell->basis.elem);

                cov[0][0] += c.x * c.x;
                cov[0][1] += c.x * c.y;
                cov[0][2] += c.x * c.z;
                cov[1][0] += c.y * c.x;
                cov[1][1] += c.y * c.y;
                cov[1][2] += c.y * c.z;
                cov[2][0] += c.z * c.x;
                cov[2][1] += c.z * c.y;
                cov[2][2] += c.z * c.z;
            }
        }
    }

    mat3_t cov_mat;
    const double scl = 1.0 / count;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            cov_mat.elem[i][j] = (float)(cov[i][j] * scl);
        }
    }

    mat3_eigen_t eigen = mat3_eigen(cov_mat);
    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    mat4_t Ri  = mat4_from_mat3(PCA);

    // Compute min and maximum extent along the PCA axes
    vec4_t min_ext = vec4_set1( FLT_MAX);
    vec4_t max_ext = vec4_set1(-FLT_MAX);
    mat3_t basis = mat3_transpose(PCA);

    // Transform the gto (x,y,z,radius) into the PCA frame to find the min and max extend within it
    for (size_t i = 0; i < count; ++i) {
        int32_t idx = in_idx ? in_idx[i] : i;
        vec4_t xyzr = in_xyzr[idx];
        vec4_t c = vec4_blend_mask(xyzr, vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
        vec4_t r = vec4_splat_w(xyzr);

        vec4_t p = mat4_mul_vec4(Ri, c);
        min_ext = vec4_min(min_ext, vec4_sub(p, r));
        max_ext = vec4_max(max_ext, vec4_add(p, r));
    }

    MEMCPY(out_basis,   &basis,   sizeof(mat3_t));
    MEMCPY(out_ext_min, &min_ext, sizeof(vec3_t));
    MEMCPY(out_ext_max, &max_ext, sizeof(vec3_t));
}

#define TRIG_ATAN2_R2_THRESHOLD 1.0e-8

static vec3_t compute_com_periodic_trig_xyz(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, vec3_t xyz_max) {
    double acc_cx = 0;
    double acc_sx = 0;
    double acc_cy = 0;
    double acc_sy = 0;
    double acc_cz = 0;
    double acc_sz = 0;
    double acc_w = 0;

    size_t i = 0;

#if defined(__AVX512F__)
    const __m512 v_scl_x = _mm512_set1_ps((float)(TWO_PI / xyz_max.x));
    const __m512 v_scl_y = _mm512_set1_ps((float)(TWO_PI / xyz_max.y));
    const __m512 v_scl_z = _mm512_set1_ps((float)(TWO_PI / xyz_max.z));
    __m512 v_acc_cx = _mm512_setzero_ps();
    __m512 v_acc_sx = _mm512_setzero_ps();
    __m512 v_acc_cy = _mm512_setzero_ps();
    __m512 v_acc_sy = _mm512_setzero_ps();
    __m512 v_acc_cz = _mm512_setzero_ps();
    __m512 v_acc_sz = _mm512_setzero_ps();
    __m512 v_acc_w = _mm512_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 16);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 v_x = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 v_y = _mm512_i32gather_ps(idx, in_y, 4);
                __m512 v_z = _mm512_i32gather_ps(idx, in_z, 4);
                __m512 v_w = _mm512_i32gather_ps(idx, in_w, 4);

                __m512 v_theta_x = _mm512_mul_ps(v_x, v_scl_x);
                __m512 v_theta_y = _mm512_mul_ps(v_y, v_scl_y);
                __m512 v_theta_z = _mm512_mul_ps(v_z, v_scl_z);

                __m512 v_cx, v_sx;
                __m512 v_cy, v_sy;
                __m512 v_cz, v_sz;
                md_mm512_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm512_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm512_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = _mm512_add_ps(v_acc_sx, _mm512_mul_ps(v_sx, v_w));
                v_acc_cx = _mm512_add_ps(v_acc_cx, _mm512_mul_ps(v_cx, v_w));

                v_acc_sy = _mm512_add_ps(v_acc_sy, _mm512_mul_ps(v_sy, v_w));
                v_acc_cy = _mm512_add_ps(v_acc_cy, _mm512_mul_ps(v_cy, v_w));

                v_acc_sz = _mm512_add_ps(v_acc_sz, _mm512_mul_ps(v_sz, v_w));
                v_acc_cz = _mm512_add_ps(v_acc_cz, _mm512_mul_ps(v_cz, v_w));

                v_acc_w = _mm512_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 v_x = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 v_y = _mm512_i32gather_ps(idx, in_y, 4);
                __m512 v_z = _mm512_i32gather_ps(idx, in_z, 4);

                __m512 v_theta_x = _mm512_mul_ps(v_x, v_scl_x);
                __m512 v_theta_y = _mm512_mul_ps(v_y, v_scl_y);
                __m512 v_theta_z = _mm512_mul_ps(v_z, v_scl_z);

                __m512 v_cx, v_sx;
                __m512 v_cy, v_sy;
                __m512 v_cz, v_sz;
                md_mm512_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm512_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm512_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = _mm512_add_ps(v_acc_sx, v_sx);
                v_acc_cx = _mm512_add_ps(v_acc_cx, v_cx);

                v_acc_sy = _mm512_add_ps(v_acc_sy, v_sy);
                v_acc_cy = _mm512_add_ps(v_acc_cy, v_cy);

                v_acc_sz = _mm512_add_ps(v_acc_sz, v_sz);
                v_acc_cz = _mm512_add_ps(v_acc_cz, v_cz);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512 v_x = _mm512_loadu_ps(in_x + i);
                __m512 v_y = _mm512_loadu_ps(in_y + i);
                __m512 v_z = _mm512_loadu_ps(in_z + i);
                __m512 v_w = _mm512_loadu_ps(in_w + i);

                __m512 v_theta_x = _mm512_mul_ps(v_x, v_scl_x);
                __m512 v_theta_y = _mm512_mul_ps(v_y, v_scl_y);
                __m512 v_theta_z = _mm512_mul_ps(v_z, v_scl_z);

                __m512 v_cx, v_sx;
                __m512 v_cy, v_sy;
                __m512 v_cz, v_sz;
                md_mm512_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm512_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm512_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = _mm512_add_ps(v_acc_sx, _mm512_mul_ps(v_sx, v_w));
                v_acc_cx = _mm512_add_ps(v_acc_cx, _mm512_mul_ps(v_cx, v_w));

                v_acc_sy = _mm512_add_ps(v_acc_sy, _mm512_mul_ps(v_sy, v_w));
                v_acc_cy = _mm512_add_ps(v_acc_cy, _mm512_mul_ps(v_cy, v_w));

                v_acc_sz = _mm512_add_ps(v_acc_sz, _mm512_mul_ps(v_sz, v_w));
                v_acc_cz = _mm512_add_ps(v_acc_cz, _mm512_mul_ps(v_cz, v_w));

                v_acc_w = _mm512_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512 v_x = _mm512_loadu_ps(in_x + i);
                __m512 v_y = _mm512_loadu_ps(in_y + i);
                __m512 v_z = _mm512_loadu_ps(in_z + i);

                __m512 v_theta_x = _mm512_mul_ps(v_x, v_scl_x);
                __m512 v_theta_y = _mm512_mul_ps(v_y, v_scl_y);
                __m512 v_theta_z = _mm512_mul_ps(v_z, v_scl_z);

                __m512 v_cx, v_sx;
                __m512 v_cy, v_sy;
                __m512 v_cz, v_sz;
                md_mm512_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm512_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm512_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = _mm512_add_ps(v_acc_sx, v_sx);
                v_acc_cx = _mm512_add_ps(v_acc_cx, v_cx);

                v_acc_sy = _mm512_add_ps(v_acc_sy, v_sy);
                v_acc_cy = _mm512_add_ps(v_acc_cy, v_cy);

                v_acc_sz = _mm512_add_ps(v_acc_sz, v_sz);
                v_acc_cz = _mm512_add_ps(v_acc_cz, v_cz);
            }
        }
    }

    acc_sx = _mm512_reduce_add_ps(v_acc_sx);
    acc_cx = _mm512_reduce_add_ps(v_acc_cx);

    acc_sy = _mm512_reduce_add_ps(v_acc_sy);
    acc_cy = _mm512_reduce_add_ps(v_acc_cy);

    acc_sz = _mm512_reduce_add_ps(v_acc_sz);
    acc_cz = _mm512_reduce_add_ps(v_acc_cz);

    acc_w  = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    const md_256 v_scl_x = md_mm256_set1_ps((float)(TWO_PI / xyz_max.x));
    const md_256 v_scl_y = md_mm256_set1_ps((float)(TWO_PI / xyz_max.y));
    const md_256 v_scl_z = md_mm256_set1_ps((float)(TWO_PI / xyz_max.z));
    md_256 v_acc_cx = md_mm256_setzero_ps();
    md_256 v_acc_sx = md_mm256_setzero_ps();
    md_256 v_acc_cy = md_mm256_setzero_ps();
    md_256 v_acc_sy = md_mm256_setzero_ps();
    md_256 v_acc_cz = md_mm256_setzero_ps();
    md_256 v_acc_sz = md_mm256_setzero_ps();
    md_256 v_acc_w = md_mm256_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 8);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 v_x = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 v_y = md_mm256_i32gather_ps(in_y, idx, 4);
                md_256 v_z = md_mm256_i32gather_ps(in_z, idx, 4);
                md_256 v_w = md_mm256_i32gather_ps(in_w, idx, 4);

                md_256 v_theta_x = md_mm256_mul_ps(v_x, v_scl_x);
                md_256 v_theta_y = md_mm256_mul_ps(v_y, v_scl_y);
                md_256 v_theta_z = md_mm256_mul_ps(v_z, v_scl_z);

                md_256 v_cx, v_sx;
                md_256 v_cy, v_sy;
                md_256 v_cz, v_sz;
                md_mm256_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm256_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm256_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm256_add_ps(v_acc_sx, md_mm256_mul_ps(v_sx, v_w));
                v_acc_cx = md_mm256_add_ps(v_acc_cx, md_mm256_mul_ps(v_cx, v_w));

                v_acc_sy = md_mm256_add_ps(v_acc_sy, md_mm256_mul_ps(v_sy, v_w));
                v_acc_cy = md_mm256_add_ps(v_acc_cy, md_mm256_mul_ps(v_cy, v_w));

                v_acc_sz = md_mm256_add_ps(v_acc_sz, md_mm256_mul_ps(v_sz, v_w));
                v_acc_cz = md_mm256_add_ps(v_acc_cz, md_mm256_mul_ps(v_cz, v_w));

                v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 v_x = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 v_y = md_mm256_i32gather_ps(in_y, idx, 4);
                md_256 v_z = md_mm256_i32gather_ps(in_z, idx, 4);

                md_256 v_theta_x = md_mm256_mul_ps(v_x, v_scl_x);
                md_256 v_theta_y = md_mm256_mul_ps(v_y, v_scl_y);
                md_256 v_theta_z = md_mm256_mul_ps(v_z, v_scl_z);

                md_256 v_cx, v_sx;
                md_256 v_cy, v_sy;
                md_256 v_cz, v_sz;
                md_mm256_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm256_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm256_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm256_add_ps(v_acc_sx, v_sx);
                v_acc_cx = md_mm256_add_ps(v_acc_cx, v_cx);

                v_acc_sy = md_mm256_add_ps(v_acc_sy, v_sy);
                v_acc_cy = md_mm256_add_ps(v_acc_cy, v_cy);

                v_acc_sz = md_mm256_add_ps(v_acc_sz, v_sz);
                v_acc_cz = md_mm256_add_ps(v_acc_cz, v_cz);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256 v_x = md_mm256_loadu_ps(in_x + i);
                md_256 v_y = md_mm256_loadu_ps(in_y + i);
                md_256 v_z = md_mm256_loadu_ps(in_z + i);
                md_256 v_w = md_mm256_loadu_ps(in_w + i);

                md_256 v_theta_x = md_mm256_mul_ps(v_x, v_scl_x);
                md_256 v_theta_y = md_mm256_mul_ps(v_y, v_scl_y);
                md_256 v_theta_z = md_mm256_mul_ps(v_z, v_scl_z);

                md_256 v_cx, v_sx;
                md_256 v_cy, v_sy;
                md_256 v_cz, v_sz;
                md_mm256_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm256_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm256_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm256_add_ps(v_acc_sx, md_mm256_mul_ps(v_sx, v_w));
                v_acc_cx = md_mm256_add_ps(v_acc_cx, md_mm256_mul_ps(v_cx, v_w));

                v_acc_sy = md_mm256_add_ps(v_acc_sy, md_mm256_mul_ps(v_sy, v_w));
                v_acc_cy = md_mm256_add_ps(v_acc_cy, md_mm256_mul_ps(v_cy, v_w));

                v_acc_sz = md_mm256_add_ps(v_acc_sz, md_mm256_mul_ps(v_sz, v_w));
                v_acc_cz = md_mm256_add_ps(v_acc_cz, md_mm256_mul_ps(v_cz, v_w));

                v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256 v_x = md_mm256_loadu_ps(in_x + i);
                md_256 v_y = md_mm256_loadu_ps(in_y + i);
                md_256 v_z = md_mm256_loadu_ps(in_z + i);

                md_256 v_theta_x = md_mm256_mul_ps(v_x, v_scl_x);
                md_256 v_theta_y = md_mm256_mul_ps(v_y, v_scl_y);
                md_256 v_theta_z = md_mm256_mul_ps(v_z, v_scl_z);

                md_256 v_cx, v_sx;
                md_256 v_cy, v_sy;
                md_256 v_cz, v_sz;
                md_mm256_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm256_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm256_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm256_add_ps(v_acc_sx, v_sx);
                v_acc_cx = md_mm256_add_ps(v_acc_cx, v_cx);

                v_acc_sy = md_mm256_add_ps(v_acc_sy, v_sy);
                v_acc_cy = md_mm256_add_ps(v_acc_cy, v_cy);

                v_acc_sz = md_mm256_add_ps(v_acc_sz, v_sz);
                v_acc_cz = md_mm256_add_ps(v_acc_cz, v_cz);
            }
        }
    }

    acc_sx = md_mm256_reduce_add_ps(v_acc_sx);
    acc_cx = md_mm256_reduce_add_ps(v_acc_cx);

    acc_sy = md_mm256_reduce_add_ps(v_acc_sy);
    acc_cy = md_mm256_reduce_add_ps(v_acc_cy);

    acc_sz = md_mm256_reduce_add_ps(v_acc_sz);
    acc_cz = md_mm256_reduce_add_ps(v_acc_cz);

    acc_w  = md_mm256_reduce_add_ps(v_acc_w);

#elif defined(__SSE2__)
    const md_128 v_scl_x = md_mm_set1_ps((float)(TWO_PI / xyz_max.x));
    const md_128 v_scl_y = md_mm_set1_ps((float)(TWO_PI / xyz_max.y));
    const md_128 v_scl_z = md_mm_set1_ps((float)(TWO_PI / xyz_max.z));
    md_128 v_acc_cx = md_mm_setzero_ps();
    md_128 v_acc_sx = md_mm_setzero_ps();
    md_128 v_acc_cy = md_mm_setzero_ps();
    md_128 v_acc_sy = md_mm_setzero_ps();
    md_128 v_acc_cz = md_mm_setzero_ps();
    md_128 v_acc_sz = md_mm_setzero_ps();
    md_128 v_acc_w = md_mm_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 4);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 v_x = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 v_y = md_mm_i32gather_ps(in_y, idx, 4);
                md_128 v_z = md_mm_i32gather_ps(in_z, idx, 4);
                md_128 v_w = md_mm_i32gather_ps(in_w, idx, 4);

                md_128 v_theta_x = md_mm_mul_ps(v_x, v_scl_x);
                md_128 v_theta_y = md_mm_mul_ps(v_y, v_scl_y);
                md_128 v_theta_z = md_mm_mul_ps(v_z, v_scl_z);

                md_128 v_cx, v_sx;
                md_128 v_cy, v_sy;
                md_128 v_cz, v_sz;
                md_mm_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm_add_ps(v_acc_sx, md_mm_mul_ps(v_sx, v_w));
                v_acc_cx = md_mm_add_ps(v_acc_cx, md_mm_mul_ps(v_cx, v_w));

                v_acc_sy = md_mm_add_ps(v_acc_sy, md_mm_mul_ps(v_sy, v_w));
                v_acc_cy = md_mm_add_ps(v_acc_cy, md_mm_mul_ps(v_cy, v_w));

                v_acc_sz = md_mm_add_ps(v_acc_sz, md_mm_mul_ps(v_sz, v_w));
                v_acc_cz = md_mm_add_ps(v_acc_cz, md_mm_mul_ps(v_cz, v_w));

                v_acc_w = md_mm_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 v_x = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 v_y = md_mm_i32gather_ps(in_y, idx, 4);
                md_128 v_z = md_mm_i32gather_ps(in_z, idx, 4);

                md_128 v_theta_x = md_mm_mul_ps(v_x, v_scl_x);
                md_128 v_theta_y = md_mm_mul_ps(v_y, v_scl_y);
                md_128 v_theta_z = md_mm_mul_ps(v_z, v_scl_z);

                md_128 v_cx, v_sx;
                md_128 v_cy, v_sy;
                md_128 v_cz, v_sz;
                md_mm_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm_add_ps(v_acc_sx, v_sx);
                v_acc_cx = md_mm_add_ps(v_acc_cx, v_cx);

                v_acc_sy = md_mm_add_ps(v_acc_sy, v_sy);
                v_acc_cy = md_mm_add_ps(v_acc_cy, v_cy);

                v_acc_sz = md_mm_add_ps(v_acc_sz, v_sz);
                v_acc_cz = md_mm_add_ps(v_acc_cz, v_cz);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128 v_x = md_mm_loadu_ps(in_x + i);
                md_128 v_y = md_mm_loadu_ps(in_y + i);
                md_128 v_z = md_mm_loadu_ps(in_z + i);
                md_128 v_w = md_mm_loadu_ps(in_w + i);

                md_128 v_theta_x = md_mm_mul_ps(v_x, v_scl_x);
                md_128 v_theta_y = md_mm_mul_ps(v_y, v_scl_y);
                md_128 v_theta_z = md_mm_mul_ps(v_z, v_scl_z);

                md_128 v_cx, v_sx;
                md_128 v_cy, v_sy;
                md_128 v_cz, v_sz;
                md_mm_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm_add_ps(v_acc_sx, md_mm_mul_ps(v_sx, v_w));
                v_acc_cx = md_mm_add_ps(v_acc_cx, md_mm_mul_ps(v_cx, v_w));

                v_acc_sy = md_mm_add_ps(v_acc_sy, md_mm_mul_ps(v_sy, v_w));
                v_acc_cy = md_mm_add_ps(v_acc_cy, md_mm_mul_ps(v_cy, v_w));

                v_acc_sz = md_mm_add_ps(v_acc_sz, md_mm_mul_ps(v_sz, v_w));
                v_acc_cz = md_mm_add_ps(v_acc_cz, md_mm_mul_ps(v_cz, v_w));

                v_acc_w = md_mm_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128 v_x = md_mm_loadu_ps(in_x + i);
                md_128 v_y = md_mm_loadu_ps(in_y + i);
                md_128 v_z = md_mm_loadu_ps(in_z + i);

                md_128 v_theta_x = md_mm_mul_ps(v_x, v_scl_x);
                md_128 v_theta_y = md_mm_mul_ps(v_y, v_scl_y);
                md_128 v_theta_z = md_mm_mul_ps(v_z, v_scl_z);

                md_128 v_cx, v_sx;
                md_128 v_cy, v_sy;
                md_128 v_cz, v_sz;
                md_mm_sincos_ps(v_theta_x, &v_sx, &v_cx);
                md_mm_sincos_ps(v_theta_y, &v_sy, &v_cy);
                md_mm_sincos_ps(v_theta_z, &v_sz, &v_cz);

                v_acc_sx = md_mm_add_ps(v_acc_sx, v_sx);
                v_acc_cx = md_mm_add_ps(v_acc_cx, v_cx);

                v_acc_sy = md_mm_add_ps(v_acc_sy, v_sy);
                v_acc_cy = md_mm_add_ps(v_acc_cy, v_cy);

                v_acc_sz = md_mm_add_ps(v_acc_sz, v_sz);
                v_acc_cz = md_mm_add_ps(v_acc_cz, v_cz);
            }
        }
    }

    acc_sx = md_mm_reduce_add_ps(v_acc_sx);
    acc_cx = md_mm_reduce_add_ps(v_acc_cx);

    acc_sy = md_mm_reduce_add_ps(v_acc_sy);
    acc_cy = md_mm_reduce_add_ps(v_acc_cy);

    acc_sz = md_mm_reduce_add_ps(v_acc_sz);
    acc_cz = md_mm_reduce_add_ps(v_acc_cz);

    acc_w  = md_mm_reduce_add_ps(v_acc_w);
#endif

    const double scl_x = TWO_PI / xyz_max.x;
    const double scl_y = TWO_PI / xyz_max.y;
    const double scl_z = TWO_PI / xyz_max.z;

    if (in_idx) {
        if (in_w) {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                double theta_x = in_x[idx] * scl_x;
                double theta_y = in_y[idx] * scl_y;
                double theta_z = in_z[idx] * scl_z;
                double w = in_w[idx];
                acc_cx += w * cos(theta_x);
                acc_sx += w * sin(theta_x);
                acc_cy += w * cos(theta_y);
                acc_sy += w * sin(theta_y);
                acc_cz += w * cos(theta_z);
                acc_sz += w * sin(theta_z);
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                double theta_x = in_x[idx] * scl_x;
                double theta_y = in_y[idx] * scl_y;
                double theta_z = in_z[idx] * scl_z;
                acc_cx += cos(theta_x);
                acc_sx += sin(theta_x);
                acc_cy += cos(theta_y);
                acc_sy += sin(theta_y);
                acc_cz += cos(theta_z);
                acc_sz += sin(theta_z);
            }
            acc_w += (double)count;
        }
    } else {
        if (in_w) {
            for (; i < count; ++i) {
                double theta_x = in_x[i] * scl_x;
                double theta_y = in_y[i] * scl_y;
                double theta_z = in_z[i] * scl_z;
                double w = in_w[i];
                acc_cx += w * cos(theta_x);
                acc_sx += w * sin(theta_x);
                acc_cy += w * cos(theta_y);
                acc_sy += w * sin(theta_y);
                acc_cz += w * cos(theta_z);
                acc_sz += w * sin(theta_z);
                acc_w  += w;
            }
        } else {
            for (; i < count; ++i) {
                double theta_x = in_x[i] * scl_x;
                double theta_y = in_y[i] * scl_y;
                double theta_z = in_z[i] * scl_z;
                acc_cx += cos(theta_x);
                acc_sx += sin(theta_x);
                acc_cy += cos(theta_y);
                acc_sy += sin(theta_y);
                acc_cz += cos(theta_z);
                acc_sz += sin(theta_z);
            }
            acc_w += (double)count;
        }
    }

    double theta_prim_x = PI;
    double theta_prim_y = PI;
    double theta_prim_z = PI;
    {
        const double y = acc_sx / acc_w;
        const double x = acc_cx / acc_w;
        const double r2 = x*x + y*y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta_prim_x += atan2(-y, -x);
        }
    }
    {
        const double y = acc_sy / acc_w;
        const double x = acc_cy / acc_w;
        const double r2 = x*x + y*y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta_prim_y += atan2(-y, -x);
        }
    }
    {
        const double y = acc_sz / acc_w;
        const double x = acc_cz / acc_w;
        const double r2 = x*x + y*y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta_prim_z += atan2(-y, -x);
        }
    }

    return vec3_set(
        (float)((theta_prim_x / TWO_PI) * xyz_max.x),
        (float)((theta_prim_y / TWO_PI) * xyz_max.y),
        (float)((theta_prim_z / TWO_PI) * xyz_max.z));
}

static float compute_com_periodic_trig(const float* in_x, const float* in_w, const int32_t* in_idx, size_t count, float x_max) {
    double acc_c = 0;
    double acc_s = 0;
    double acc_w = 0;

    const double scl = TWO_PI / x_max;
    size_t i = 0;

#if defined(__AVX512F__)
    const __m512 v_scl = _mm512_set1_ps((float)(TWO_PI / x_max));
    __m512 v_acc_c = _mm512_setzero_ps();
    __m512 v_acc_s = _mm512_setzero_ps();
    __m512 v_acc_w = _mm512_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 16);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 v_x = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 v_w = _mm512_i32gather_ps(idx, in_w, 4);
                __m512 v_theta = _mm512_mul_ps(v_x, v_scl);
                __m512 v_c, v_s;
                md_mm512_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = _mm512_add_ps(v_acc_s, _mm512_mul_ps(v_s, v_w));
                v_acc_c = _mm512_add_ps(v_acc_c, _mm512_mul_ps(v_c, v_w));
                v_acc_w = _mm512_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 v_x = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 v_theta = _mm512_mul_ps(v_x, v_scl);
                __m512 v_c, v_s;
                md_mm512_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = _mm512_add_ps(v_acc_s, v_s);
                v_acc_c = _mm512_add_ps(v_acc_c, v_c);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512 v_x = _mm512_loadu_ps(in_x + i);
                __m512 v_w = _mm512_loadu_ps(in_w + i);
                __m512 v_theta = _mm512_mul_ps(v_x, v_scl);
                __m512 v_c, v_s;
                md_mm512_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = _mm512_add_ps(v_acc_s, _mm512_mul_ps(v_s, v_w));
                v_acc_c = _mm512_add_ps(v_acc_c, _mm512_mul_ps(v_c, v_w));
                v_acc_w = _mm512_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512 v_x = _mm512_loadu_ps(in_x + i);
                __m512 v_theta = _mm512_mul_ps(v_x, v_scl);
                __m512 v_c, v_s;
                md_mm512_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = _mm512_add_ps(v_acc_s, v_s);
                v_acc_c = _mm512_add_ps(v_acc_c, v_c);
            }
        }
    }

    acc_s = _mm512_reduce_add_ps(v_acc_s);
    acc_c = _mm512_reduce_add_ps(v_acc_c);
    acc_w = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    const md_256 v_scl = md_mm256_set1_ps((float)(TWO_PI / x_max));
    md_256 v_acc_c = md_mm256_setzero_ps();
    md_256 v_acc_s = md_mm256_setzero_ps();
    md_256 v_acc_w = md_mm256_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 8);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 v_x = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 v_w = md_mm256_i32gather_ps(in_w, idx, 4);
                md_256 v_theta = md_mm256_mul_ps(v_x, v_scl);
                md_256 v_c, v_s;
                md_mm256_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm256_add_ps(v_acc_s, md_mm256_mul_ps(v_s, v_w));
                v_acc_c = md_mm256_add_ps(v_acc_c, md_mm256_mul_ps(v_c, v_w));
                v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 v_x = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 v_theta = md_mm256_mul_ps(v_x, v_scl);
                md_256 v_c, v_s;
                md_mm256_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm256_add_ps(v_acc_s, v_s);
                v_acc_c = md_mm256_add_ps(v_acc_c, v_c);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256 v_x = md_mm256_loadu_ps(in_x + i);
                md_256 v_w = md_mm256_loadu_ps(in_w + i);
                md_256 v_theta = md_mm256_mul_ps(v_x, v_scl);
                md_256 v_c, v_s;
                md_mm256_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm256_add_ps(v_acc_s, md_mm256_mul_ps(v_s, v_w));
                v_acc_c = md_mm256_add_ps(v_acc_c, md_mm256_mul_ps(v_c, v_w));
                v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256 v_x = md_mm256_loadu_ps(in_x + i);
                md_256 v_theta = md_mm256_mul_ps(v_x, v_scl);
                md_256 v_c, v_s;
                md_mm256_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm256_add_ps(v_acc_s, v_s);
                v_acc_c = md_mm256_add_ps(v_acc_c, v_c);
            }
        }
    }

    acc_s = md_mm256_reduce_add_ps(v_acc_s);
    acc_c = md_mm256_reduce_add_ps(v_acc_c);
    acc_w = md_mm256_reduce_add_ps(v_acc_w);
#elif defined(__SSE2__)
    const md_128 v_scl = md_mm_set1_ps((float)(TWO_PI / x_max));
    md_128 v_acc_c = md_mm_setzero_ps();
    md_128 v_acc_s = md_mm_setzero_ps();
    md_128 v_acc_w = md_mm_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 4);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 v_x = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 v_w = md_mm_i32gather_ps(in_w, idx, 4);
                md_128 v_theta = md_mm_mul_ps(v_x, v_scl);
                md_128 v_c, v_s;
                md_mm_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm_add_ps(v_acc_s, md_mm_mul_ps(v_s, v_w));
                v_acc_c = md_mm_add_ps(v_acc_c, md_mm_mul_ps(v_c, v_w));
                v_acc_w = md_mm_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 v_x = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 v_theta = md_mm_mul_ps(v_x, v_scl);
                md_128 v_c, v_s;
                md_mm_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm_add_ps(v_acc_s, v_s);
                v_acc_c = md_mm_add_ps(v_acc_c, v_c);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128 v_x = md_mm_loadu_ps(in_x + i);
                md_128 v_w = md_mm_loadu_ps(in_w + i);
                md_128 v_theta = md_mm_mul_ps(v_x, v_scl);
                md_128 v_c, v_s;
                md_mm_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm_add_ps(v_acc_s, md_mm_mul_ps(v_s, v_w));
                v_acc_c = md_mm_add_ps(v_acc_c, md_mm_mul_ps(v_c, v_w));
                v_acc_w = md_mm_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128 v_x = md_mm_loadu_ps(in_x + i);
                md_128 v_theta = md_mm_mul_ps(v_x, v_scl);
                md_128 v_c, v_s;
                md_mm_sincos_ps(v_theta, &v_s, &v_c);
                v_acc_s = md_mm_add_ps(v_acc_s, v_s);
                v_acc_c = md_mm_add_ps(v_acc_c, v_c);
            }
        }
    }

    acc_s = md_mm_reduce_add_ps(v_acc_s);
    acc_c = md_mm_reduce_add_ps(v_acc_c);
    acc_w = md_mm_reduce_add_ps(v_acc_w);
#endif

    if (in_idx) {
        if (in_w) {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                double theta = in_x[idx] * scl;
                double w = in_w[idx];
                acc_c += w * cos(theta);
                acc_s += w * sin(theta);
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                double theta = in_x[idx] * scl;
                acc_c += cos(theta);
                acc_s += sin(theta);
            }
            acc_w += (double)count;
        }
    } else {
        if (in_w) {
            for (; i < count; ++i) {
                double theta = in_x[i] * scl;
                double w = in_w[i];
                acc_c += w * cos(theta);
                acc_s += w * sin(theta);
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                double theta = in_x[i] * scl;
                acc_c += cos(theta);
                acc_s += sin(theta);
            }
            acc_w += (double)count;
        }
    }

    const double y = acc_s / acc_w;
    const double x = acc_c / acc_w;
    const double r2 = x*x + y*y;

    double theta_prim = PI;
    if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
        theta_prim += atan2(-y, -x);
    }

    return (float)((theta_prim / TWO_PI) * x_max);
}

static float compute_com_periodic_reg(const float* in_x, const float* in_w, const int32_t* in_idx, size_t count, float x_max) {
    if (count <= 0) {
        return 0;
    }

    double acc_x = 0;
    double acc_w = 0;

    size_t i = 0;
#if defined(__AVX512F__)
    const size_t simd_count = ROUND_DOWN(count, 16);
    if (simd_count > 0) {
        const __m512 v_ext = _mm512_set1_ps(x_max);
        __m512 v_acc_x = _mm512_setzero_ps();
        __m512 v_acc_w = _mm512_setzero_ps();
        i += 8;
        if (in_idx) {
            __m512i idx = _mm512_loadu_si512(in_idx);
            if (in_w) {
                v_acc_x = _mm512_i32gather_ps(idx, in_x, 4);
                v_acc_w = _mm512_i32gather_ps(idx, in_w, 4);
                for (; i < simd_count; i += 16) {
                    idx = _mm512_loadu_si512(in_idx + i);
                    __m512 r = _mm512_div_ps(v_acc_x, v_acc_w);
                    __m512 x = _mm512_i32gather_ps(idx, in_x, 4);
                    __m512 w = _mm512_i32gather_ps(idx, in_w, 4);
                    x = md_mm512_deperiodize_ps(x, r, v_ext);
                    v_acc_x = _mm512_add_ps(v_acc_x, _mm512_mul_ps(x, w));
                    v_acc_w = _mm512_add_ps(v_acc_w, w);
                }
            } else {
                v_acc_x = _mm512_i32gather_ps(idx, in_x, 4);
                v_acc_w = _mm512_set1_ps(16);
                for (; i < simd_count; i += 16) {
                    idx = _mm512_loadu_si512(in_idx + i);
                    __m512 r = _mm512_div_ps(v_acc_x, v_acc_w);
                    __m512 x = _mm512_i32gather_ps(idx, in_x, 4);
                    x = md_mm512_deperiodize_ps(x, r, v_ext);
                    v_acc_x = _mm512_add_ps(v_acc_x, x);
                    v_acc_w = _mm512_add_ps(v_acc_w, _mm512_set1_ps(8));
                }
            }
        } else {
            if (in_w) {
                v_acc_x = _mm512_loadu_ps(in_x);
                v_acc_w = _mm512_loadu_ps(in_w);
                for (; i < simd_count; i += 16) {
                    __m512 r = _mm512_div_ps(v_acc_x, v_acc_w);
                    __m512 x = _mm512_loadu_ps(in_x + i);
                    __m512 w = _mm512_loadu_ps(in_w + i);
                    x = md_mm512_deperiodize_ps(x, r, v_ext);
                    v_acc_x = _mm512_add_ps(v_acc_x, _mm512_mul_ps(x, w));
                    v_acc_w = _mm512_add_ps(v_acc_w, w);
                }
            } else {
                v_acc_x = _mm512_loadu_ps(in_x);
                v_acc_w = _mm512_set1_ps(16);
                for (; i < simd_count; i += 16) {
                    __m512 r = _mm512_div_ps(v_acc_x, v_acc_w);
                    __m512 x = _mm512_loadu_ps(in_x + i);
                    x = md_mm512_deperiodize_ps(x, r, v_ext);
                    v_acc_x = _mm512_add_ps(v_acc_x, x);
                    v_acc_w = _mm512_add_ps(v_acc_w, _mm512_set1_ps(16));
                }
            }
        }
        acc_x = _mm512_reduce_add_ps(v_acc_x);
        acc_w = _mm512_reduce_add_ps(v_acc_w);
    }
#elif defined(__AVX2__)
    const size_t simd_count = ROUND_DOWN(count, 8);
    if (simd_count > 0) {
       const md_256 v_ext = md_mm256_set1_ps(x_max);
        md_256 v_acc_x = md_mm256_setzero_ps();
        md_256 v_acc_w = md_mm256_setzero_ps();
        i += 8;
        if (in_idx) {
            md_256i idx = md_mm256_loadu_si256(in_idx);
            if (in_w) {
                v_acc_x = md_mm256_i32gather_ps(in_x, idx, 4);
                v_acc_w = md_mm256_i32gather_ps(in_w, idx, 4);
                for (; i < simd_count; i += 8) {
                    idx = md_mm256_loadu_si256(in_idx + i);
                    md_256 r = md_mm256_div_ps(v_acc_x, v_acc_w);
                    md_256 x = md_mm256_i32gather_ps(in_x, idx, 4);
                    md_256 w = md_mm256_i32gather_ps(in_w, idx, 4);
                    x = md_mm256_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm256_add_ps(v_acc_x, md_mm256_mul_ps(x, w));
                    v_acc_w = md_mm256_add_ps(v_acc_w, w);
                }
            } else {
                v_acc_x = md_mm256_i32gather_ps(in_x, idx, 4);
                v_acc_w = md_mm256_set1_ps(8);
                for (; i < simd_count; i += 8) {
                    idx = md_mm256_loadu_si256(in_idx + i);
                    md_256 r = md_mm256_div_ps(v_acc_x, v_acc_w);
                    md_256 x = md_mm256_i32gather_ps(in_x, idx, 4);
                    x = md_mm256_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm256_add_ps(v_acc_x, x);
                    v_acc_w = md_mm256_add_ps(v_acc_w, md_mm256_set1_ps(8));
                }
            }
        } else {
            if (in_w) {
                v_acc_x = md_mm256_loadu_ps(in_x);
                v_acc_w = md_mm256_loadu_ps(in_w);
                for (; i < simd_count; i += 8) {
                    md_256 r = md_mm256_div_ps(v_acc_x, v_acc_w);
                    md_256 x = md_mm256_loadu_ps(in_x + i);
                    md_256 w = md_mm256_loadu_ps(in_w + i);
                    x = md_mm256_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm256_add_ps(v_acc_x, md_mm256_mul_ps(x, w));
                    v_acc_w = md_mm256_add_ps(v_acc_w, w);
                }
            } else {
                v_acc_x = md_mm256_loadu_ps(in_x);
                v_acc_w = md_mm256_set1_ps(8);
                for (; i < simd_count; i += 8) {
                    md_256 r = md_mm256_div_ps(v_acc_x, v_acc_w);
                    md_256 x = md_mm256_loadu_ps(in_x + i);
                    x = md_mm256_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm256_add_ps(v_acc_x, x);
                    v_acc_w = md_mm256_add_ps(v_acc_w, md_mm256_set1_ps(8));
                }
            }
        }
        acc_x = md_mm256_reduce_add_ps(v_acc_x);
        acc_w = md_mm256_reduce_add_ps(v_acc_w);
    }
#elif defined(__SSE2__)
    const size_t simd_count = ROUND_DOWN(count, 4);
    if (simd_count > 0) {
        const md_128 v_ext = md_mm_set1_ps(x_max);
        md_128 v_acc_x = md_mm_setzero_ps();
        md_128 v_acc_w = md_mm_setzero_ps();
        i += 8;
        if (in_idx) {
            md_128i idx = md_mm_loadu_si128(in_idx);
            if (in_w) {
                v_acc_x = md_mm_i32gather_ps(in_x, idx, 4);
                v_acc_w = md_mm_i32gather_ps(in_w, idx, 4);
                for (; i < simd_count; i += 4) {
                    idx = md_mm_loadu_si128(in_idx + i);
                    md_128 r = md_mm_div_ps(v_acc_x, v_acc_w);
                    md_128 x = md_mm_i32gather_ps(in_x, idx, 4);
                    md_128 w = md_mm_i32gather_ps(in_w, idx, 4);
                    x = md_mm_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm_add_ps(v_acc_x, md_mm_mul_ps(x, w));
                    v_acc_w = md_mm_add_ps(v_acc_w, w);
                }
            } else {
                v_acc_x = md_mm_i32gather_ps(in_x, idx, 4);
                v_acc_w = md_mm_set1_ps(4);
                for (; i < simd_count; i += 4) {
                    idx = md_mm_loadu_si128(in_idx + i);
                    md_128 r = md_mm_div_ps(v_acc_x, v_acc_w);
                    md_128 x = md_mm_i32gather_ps(in_x, idx, 4);
                    x = md_mm_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm_add_ps(v_acc_x, x);
                    v_acc_w = md_mm_add_ps(v_acc_w, md_mm_set1_ps(4));
                }
            }
        } else {
            if (in_w) {
                v_acc_x = md_mm_loadu_ps(in_x);
                v_acc_w = md_mm_loadu_ps(in_w);
                for (; i < simd_count; i += 4) {
                    md_128 r = md_mm_div_ps(v_acc_x, v_acc_w);
                    md_128 x = md_mm_loadu_ps(in_x + i);
                    md_128 w = md_mm_loadu_ps(in_w + i);
                    x = md_mm_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm_add_ps(v_acc_x, md_mm_mul_ps(x, w));
                    v_acc_w = md_mm_add_ps(v_acc_w, w);
                }
            } else {
                v_acc_x = md_mm_loadu_ps(in_x);
                v_acc_w = md_mm_set1_ps(4);
                for (; i < simd_count; i += 4) {
                    md_128 r = md_mm_div_ps(v_acc_x, v_acc_w);
                    md_128 x = md_mm_loadu_ps(in_x + i);
                    x = md_mm_deperiodize_ps(x, r, v_ext);
                    v_acc_x = md_mm_add_ps(v_acc_x, x);
                    v_acc_w = md_mm_add_ps(v_acc_w, md_mm_set1_ps(4));
                }
            }
        }
        acc_x = md_mm_reduce_add_ps(v_acc_x);
        acc_w = md_mm_reduce_add_ps(v_acc_w);
    }
#endif

    if (i == 0) {
        int64_t idx = in_idx ? in_idx[i] : (int64_t)i;
        acc_x = in_x[idx];
        acc_w = in_w ? in_w[idx] : 1.0;
    }
    const double d_x_max = x_max;

    if (in_idx) {
        if (in_w) {
            int64_t idx = in_idx[i];
            double r = acc_x / acc_w;
            double x = in_x[idx];
            double w = in_w[idx];
            double dx = deperiodize_ortho(x, r, d_x_max);
            acc_x += dx * w;
            acc_w += w;
        } else {
            int64_t idx = in_idx[i];
            double r = acc_x / acc_w;
            double x = in_x[idx];
            double dx = deperiodize_ortho(x, r, d_x_max);
            acc_x += dx;
            acc_w += 1.0;
        }
    } else {
        if (in_w) {
            double r = acc_x / acc_w;
            double x = in_x[i];
            double w = in_w[i];
            double dx = deperiodize_ortho(x, r, d_x_max);
            acc_x += dx * w;
            acc_w += w;
        } else {
            double r = acc_x / acc_w;
            double x = in_x[i];
            double dx = deperiodize_ortho(x, r, d_x_max);
            acc_x += dx;
            acc_w += 1.0;
        }
    }

    return (float)(acc_x / acc_w);
}

static float compute_com(const float* in_x, const float* in_w, const int32_t* in_idx, size_t count) {
    ASSERT(in_x);

    if (count <= 0)
        return 0.0f;

    size_t i = 0;
    double acc_x = 0;
    double acc_w = 0;

#if defined(__AVX512F__)
    __m512 v_acc_x = _mm512_setzero_ps();
    __m512 v_acc_w = _mm512_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 16);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 v_x  = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 v_w  = _mm512_i32gather_ps(idx, in_w, 4);
                v_acc_x = _mm512_add_ps(v_acc_x, _mm512_mul_ps(v_x, v_w));
                v_acc_w = _mm512_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 v_x  = _mm512_i32gather_ps(idx, in_x, 4);
                v_acc_x = _mm512_add_ps(v_acc_x, v_x);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512 v_x = _mm512_loadu_ps(in_x + i);
                __m512 v_w = _mm512_loadu_ps(in_w + i);
                v_acc_x = _mm512_add_ps(v_acc_x, _mm512_mul_ps(v_x, v_w));
                v_acc_w = _mm512_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512 v_x = _mm512_loadu_ps(in_x + i);
                v_acc_x = _mm512_add_ps(v_acc_x, v_x);
            }
        }
    }

    acc_x = _mm512_reduce_add_ps(v_acc_x);
    acc_w = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    md_256 v_acc_x = md_mm256_setzero_ps();
    md_256 v_acc_w = md_mm256_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 8);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 v_x  = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 v_w  = md_mm256_i32gather_ps(in_w, idx, 4);
                v_acc_x = md_mm256_add_ps(v_acc_x, md_mm256_mul_ps(v_x, v_w));
                v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 v_x  = md_mm256_i32gather_ps(in_x, idx, 4);
                v_acc_x = md_mm256_add_ps(v_acc_x, v_x);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256 v_x = md_mm256_loadu_ps(in_x + i);
                md_256 v_w = md_mm256_loadu_ps(in_w + i);
                v_acc_x = md_mm256_add_ps(v_acc_x, md_mm256_mul_ps(v_x, v_w));
                v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256 v_x = md_mm256_loadu_ps(in_x + i);
                v_acc_x = md_mm256_add_ps(v_acc_x, v_x);
            }
        }
    }

    acc_x = md_mm256_reduce_add_ps(v_acc_x);
    acc_w = md_mm256_reduce_add_ps(v_acc_w);
#elif defined(__SSE2__)
    md_128 v_acc_x = md_mm_setzero_ps();
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 4);
    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 v_x  = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 v_w  = md_mm_i32gather_ps(in_w, idx, 4);
                v_acc_x = md_mm_add_ps(v_acc_x, md_mm_mul_ps(v_x, v_w));
                v_acc_w = md_mm_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 v_x  = md_mm_i32gather_ps(in_x, idx, 4);
                v_acc_x = md_mm_add_ps(v_acc_x, v_x);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128 v_x = md_mm_loadu_ps(in_x + i);
                md_128 v_w = md_mm_loadu_ps(in_w + i);
                v_acc_x = md_mm_add_ps(v_acc_x, md_mm_mul_ps(v_x, v_w));
                v_acc_w = md_mm_add_ps(v_acc_w, v_w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128 v_x = md_mm_loadu_ps(in_x + i);
                v_acc_x = md_mm_add_ps(v_acc_x, v_x);
            }
        }
    }

    acc_x = md_mm_reduce_add_ps(v_acc_x);
    acc_w = md_mm_reduce_add_ps(v_acc_w);
#endif

    if (in_idx) {
        if (in_w) {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                float w = in_w[idx];
                acc_x += in_x[idx] * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                acc_x += in_x[idx];
            }
            acc_w += count;
        }
    } else {
        if (in_w) {
            for (; i < count; ++i) {
                float w = in_w[i];
                acc_x += in_x[i] * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                acc_x += in_x[i];
            }
            acc_w += count;
        }
    }

    float com = (float)(acc_x / acc_w);
    return com;
}

static void com(float* out_com, const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count) {
    size_t i = 0;
    double acc_x = 0;
    double acc_y = 0;
    double acc_z = 0;
    double acc_w = 0;
#if defined(__AVX512F__)
    __m512 vx = _mm512_setzero_ps();
    __m512 vy = _mm512_setzero_ps();
    __m512 vz = _mm512_setzero_ps();
    __m512 vw = _mm512_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 16);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 x    = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 y    = _mm512_i32gather_ps(idx, in_y, 4);
                __m512 z    = _mm512_i32gather_ps(idx, in_z, 4);
                __m512 w    = _mm512_i32gather_ps(idx, in_w, 4);

                vx = _mm512_add_ps(vx, _mm512_mul_ps(x, w));
                vy = _mm512_add_ps(vy, _mm512_mul_ps(y, w));
                vz = _mm512_add_ps(vz, _mm512_mul_ps(z, w));
                vw = _mm512_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 x    = _mm512_i32gather_ps(idx, in_x, 4);
                __m512 y    = _mm512_i32gather_ps(idx, in_y, 4);
                __m512 z    = _mm512_i32gather_ps(idx, in_z, 4);

                vx = _mm512_add_ps(vx, x);
                vy = _mm512_add_ps(vy, y);
                vz = _mm512_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512 x = _mm512_loadu_ps(in_x + i);
                __m512 y = _mm512_loadu_ps(in_y + i);
                __m512 z = _mm512_loadu_ps(in_z + i);
                __m512 w = _mm512_loadu_ps(in_w + i);

                vx = _mm512_add_ps(vx, _mm512_mul_ps(x, w));
                vy = _mm512_add_ps(vy, _mm512_mul_ps(y, w));
                vz = _mm512_add_ps(vz, _mm512_mul_ps(z, w));
                vw = _mm512_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512 x = _mm512_loadu_ps(in_x + i);
                __m512 y = _mm512_loadu_ps(in_y + i);
                __m512 z = _mm512_loadu_ps(in_z + i);

                vx = _mm512_add_ps(vx, x);
                vy = _mm512_add_ps(vy, y);
                vz = _mm512_add_ps(vz, z);
            }
        }
    }

    acc_x = _mm512_reduce_add_ps(vx);
    acc_y = _mm512_reduce_add_ps(vy);
    acc_z = _mm512_reduce_add_ps(vz);
    acc_w = _mm512_reduce_add_ps(vw);
#elif defined (__AVX__)
    md_256 vx = md_mm256_setzero_ps();
    md_256 vy = md_mm256_setzero_ps();
    md_256 vz = md_mm256_setzero_ps();
    md_256 vw = md_mm256_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 8);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 x    = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 y    = md_mm256_i32gather_ps(in_y, idx, 4);
                md_256 z    = md_mm256_i32gather_ps(in_z, idx, 4);
                md_256 w    = md_mm256_i32gather_ps(in_w, idx, 4);

                vx = md_mm256_add_ps(vx, md_mm256_mul_ps(x, w));
                vy = md_mm256_add_ps(vy, md_mm256_mul_ps(y, w));
                vz = md_mm256_add_ps(vz, md_mm256_mul_ps(z, w));
                vw = md_mm256_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 x    = md_mm256_i32gather_ps(in_x, idx, 4);
                md_256 y    = md_mm256_i32gather_ps(in_y, idx, 4);
                md_256 z    = md_mm256_i32gather_ps(in_z, idx, 4);

                vx = md_mm256_add_ps(vx, x);
                vy = md_mm256_add_ps(vy, y);
                vz = md_mm256_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256 x = md_mm256_loadu_ps(in_x + i);
                md_256 y = md_mm256_loadu_ps(in_y + i);
                md_256 z = md_mm256_loadu_ps(in_z + i);
                md_256 w = md_mm256_loadu_ps(in_w + i);

                vx = md_mm256_add_ps(vx, md_mm256_mul_ps(x, w));
                vy = md_mm256_add_ps(vy, md_mm256_mul_ps(y, w));
                vz = md_mm256_add_ps(vz, md_mm256_mul_ps(z, w));
                vw = md_mm256_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256 x = md_mm256_loadu_ps(in_x + i);
                md_256 y = md_mm256_loadu_ps(in_y + i);
                md_256 z = md_mm256_loadu_ps(in_z + i);

                vx = md_mm256_add_ps(vx, x);
                vy = md_mm256_add_ps(vy, y);
                vz = md_mm256_add_ps(vz, z);
            }
        }
    }

    acc_x = md_mm256_reduce_add_ps(vx);
    acc_y = md_mm256_reduce_add_ps(vy);
    acc_z = md_mm256_reduce_add_ps(vz);
    acc_w = md_mm256_reduce_add_ps(vw);

#elif defined(__SSE2__)
    md_128 vx = md_mm_setzero_ps();
    md_128 vy = md_mm_setzero_ps();
    md_128 vz = md_mm_setzero_ps();
    md_128 vw = md_mm_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 4);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 x    = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 y    = md_mm_i32gather_ps(in_y, idx, 4);
                md_128 z    = md_mm_i32gather_ps(in_z, idx, 4);
                md_128 w    = md_mm_i32gather_ps(in_w, idx, 4);

                vx = md_mm_add_ps(vx, md_mm_mul_ps(x, w));
                vy = md_mm_add_ps(vy, md_mm_mul_ps(y, w));
                vz = md_mm_add_ps(vz, md_mm_mul_ps(z, w));
                vw = md_mm_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 x    = md_mm_i32gather_ps(in_x, idx, 4);
                md_128 y    = md_mm_i32gather_ps(in_y, idx, 4);
                md_128 z    = md_mm_i32gather_ps(in_z, idx, 4);

                vx = md_mm_add_ps(vx, x);
                vy = md_mm_add_ps(vy, y);
                vz = md_mm_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128 x = md_mm_loadu_ps(in_x + i);
                md_128 y = md_mm_loadu_ps(in_y + i);
                md_128 z = md_mm_loadu_ps(in_z + i);
                md_128 w = md_mm_loadu_ps(in_w + i);

                vx = _mm_add_ps(vx, md_mm_mul_ps(x, w));
                vy = _mm_add_ps(vy, md_mm_mul_ps(y, w));
                vz = _mm_add_ps(vz, md_mm_mul_ps(z, w));
                vw = _mm_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128 x = md_mm_loadu_ps(in_x + i);
                md_128 y = md_mm_loadu_ps(in_y + i);
                md_128 z = md_mm_loadu_ps(in_z + i);

                vx = md_mm_add_ps(vx, x);
                vy = md_mm_add_ps(vy, y);
                vz = md_mm_add_ps(vz, z);
            }
        }
    }

    acc_x = md_mm_reduce_add_ps(vx);
    acc_y = md_mm_reduce_add_ps(vy);
    acc_z = md_mm_reduce_add_ps(vz);
    acc_w = md_mm_reduce_add_ps(vw);
#endif

    if (in_idx) {
        if (in_w) {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                float w = in_w[idx];
                acc_x += in_x[idx] * w;
                acc_y += in_y[idx] * w;
                acc_z += in_z[idx] * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                acc_x += in_x[idx];
                acc_y += in_y[idx];
                acc_z += in_z[idx];
            }
            acc_w = (double)count;
        }
    } else {
        if (in_w) {
            for (; i < count; ++i) {
                float w = in_w[i];
                acc_x += in_x[i] * w;
                acc_y += in_y[i] * w;
                acc_z += in_z[i] * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                acc_x += in_x[i];
                acc_y += in_y[i];
                acc_z += in_z[i];
            }
            acc_w = (double)count;
        }
    }

    out_com[0] = (float)(acc_x / acc_w);
    out_com[1] = (float)(acc_y / acc_w);
    out_com[2] = (float)(acc_z / acc_w);
}

// Internal versions for COM computation supporting triclinic and ortho PBC
static void _com_pbc_w(float out_com[3], const float* in_x, const float* in_y, const float* in_z, const float* in_w, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(in_w);
    ASSERT(M);
    ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};
    double acc_w = 0;

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    md_512 v_acc_w = _mm512_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512 v_x = _mm512_loadu_ps(in_x + i);
        md_512 v_y = _mm512_loadu_ps(in_y + i);
        md_512 v_z = _mm512_loadu_ps(in_z + i);
        md_512 v_w = _mm512_loadu_ps(in_w + i);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_mul_ps(v_x, _mm512_set1_ps(M[0][2]))));
        md_512 v_ty = _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_y, _mm512_set1_ps(M[1][2]))));
        md_512 v_tz = _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][0]), _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = _mm512_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = _mm512_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = _mm512_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = _mm512_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = _mm512_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = _mm512_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
    acc_w    = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    md_256 v_acc_w = md_mm256_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256 v_x = md_mm256_loadu_ps(in_x + i);
        md_256 v_y = md_mm256_loadu_ps(in_y + i);
        md_256 v_z = md_mm256_loadu_ps(in_z + i);
        md_256 v_w = md_mm256_loadu_ps(in_w + i);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_mul_ps(v_x, md_mm256_set1_ps(M[0][2]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_y, md_mm256_set1_ps(M[1][2]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][0]), md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm256_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm256_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm256_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm256_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm256_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm256_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm256_reduce_add_ps(v_acc_w);
#elif defined(__SSE2__)
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128 v_x = md_mm_loadu_ps(in_x + i);
        md_128 v_y = md_mm_loadu_ps(in_y + i);
        md_128 v_z = md_mm_loadu_ps(in_z + i);
        md_128 v_w = md_mm_loadu_ps(in_w + i);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_mul_ps(v_x, md_mm_set1_ps(M[0][2]))));
        md_128 v_ty = md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_y, md_mm_set1_ps(M[1][2]))));
        md_128 v_tz = md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][0]), md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm_reduce_add_ps(v_acc_w);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        double x = in_x[i];
        double y = in_y[i];
        double z = in_z[i];
        double w = in_w[i];
        double tx = x * M[0][0] + x * M[0][1] + x * M[0][2];
        double ty = y * M[1][0] + y * M[1][1] + y * M[1][2];
        double tz = z * M[2][0] + z * M[2][1] + z * M[2][2];
        acc_c[0] += w * cos(tx);
        acc_s[0] += w * sin(tx);
        acc_c[1] += w * cos(ty);
        acc_s[1] += w * sin(ty);
        acc_c[2] += w * cos(tz);
        acc_s[2] += w * sin(tz);
        acc_w    += w;
    }

    const double inv_w = 1.0 / acc_w;
    for (int j = 0; j < 3; ++j) {
        double theta = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta += atan2(-y, -x);
        }
        // This is essentially a matrix vector multiplication, but for the single row
        out_com[j] = (float)(theta * I[j][0] + theta * I[j][1] + theta * I[j][2]);
    }
}

// Internal versions for COM computation supporting triclinic and ortho PBC
static void _com_pbc(float out_com[3], const float* in_x, const float* in_y, const float* in_z, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(M);
    ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512 v_x = _mm512_loadu_ps(in_x + i);
        md_512 v_y = _mm512_loadu_ps(in_y + i);
        md_512 v_z = _mm512_loadu_ps(in_z + i);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_mul_ps(v_x, _mm512_set1_ps(M[0][2]))));
        md_512 v_ty = _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_y, _mm512_set1_ps(M[1][2]))));
        md_512 v_tz = _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][0]), _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = _mm512_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = _mm512_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = _mm512_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = _mm512_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = _mm512_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256 v_x = md_mm256_loadu_ps(in_x + i);
        md_256 v_y = md_mm256_loadu_ps(in_y + i);
        md_256 v_z = md_mm256_loadu_ps(in_z + i);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_mul_ps(v_x, md_mm256_set1_ps(M[0][2]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_y, md_mm256_set1_ps(M[1][2]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][0]), md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm256_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm256_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm256_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm256_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm256_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
#elif defined(__SSE2__)
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128 v_x = md_mm_loadu_ps(in_x + i);
        md_128 v_y = md_mm_loadu_ps(in_y + i);
        md_128 v_z = md_mm_loadu_ps(in_z + i);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_mul_ps(v_x, md_mm_set1_ps(M[0][2]))));
        md_128 v_ty = md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_y, md_mm_set1_ps(M[1][2]))));
        md_128 v_tz = md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][0]), md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        double x = in_x[i];
        double y = in_y[i];
        double z = in_z[i];
        double tx = x * M[0][0] + x * M[0][1] + x * M[0][2];
        double ty = y * M[1][0] + y * M[1][1] + y * M[1][2];
        double tz = z * M[2][0] + z * M[2][1] + z * M[2][2];
        acc_c[0] += cos(tx);
        acc_s[0] += sin(tx);
        acc_c[1] += cos(ty);
        acc_s[1] += sin(ty);
        acc_c[2] += cos(tz);
        acc_s[2] += sin(tz);
    }

    const double inv_w = 1.0 / (double)count;
    for (int j = 0; j < 3; ++j) {
        double theta = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta += atan2(-y, -x);
        }
        // This is essentially a matrix vector multiplication, but for the single row
        out_com[j] = (float)(theta * I[j][0] + theta * I[j][1] + theta * I[j][2]);
    }
}


static void _com_pbc_i(float out_com[3], const float* in_x, const float* in_y, const float* in_z, const int32_t* in_idx, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(M);
    ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512i idx = _mm512_loadu_si512(in_idx + i);
        md_512 v_x  = _mm512_i32gather_ps(idx, in_x, 4);
        md_512 v_y  = _mm512_i32gather_ps(idx, in_y, 4);
        md_512 v_z  = _mm512_i32gather_ps(idx, in_z, 4);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_mul_ps(v_x, _mm512_set1_ps(M[0][2]))));
        md_512 v_ty = _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_y, _mm512_set1_ps(M[1][2]))));
        md_512 v_tz = _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][0]), _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = _mm512_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = _mm512_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = _mm512_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = _mm512_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = _mm512_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256i idx = md_mm256_loadu_si256(in_idx + i);
        md_256 v_x  = md_mm256_i32gather_ps(in_x, idx, 4);
        md_256 v_y  = md_mm256_i32gather_ps(in_y, idx, 4);
        md_256 v_z  = md_mm256_i32gather_ps(in_z, idx, 4);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_mul_ps(v_x, md_mm256_set1_ps(M[0][2]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_y, md_mm256_set1_ps(M[1][2]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][0]), md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm256_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm256_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm256_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm256_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm256_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
#elif defined(__SSE2__)
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128i idx = md_mm_loadu_si128(in_idx + i);
        md_128 v_x  = md_mm_i32gather_ps(in_x, idx, 4);
        md_128 v_y  = md_mm_i32gather_ps(in_y, idx, 4);
        md_128 v_z  = md_mm_i32gather_ps(in_z, idx, 4);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_mul_ps(v_x, md_mm_set1_ps(M[0][2]))));
        md_128 v_ty = md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_y, md_mm_set1_ps(M[1][2]))));
        md_128 v_tz = md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][0]), md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        int32_t idx = in_idx[i];
        double x = in_x[idx];
        double y = in_y[idx];
        double z = in_z[idx];
        double tx = x * M[0][0] + x * M[0][1] + x * M[0][2];
        double ty = y * M[1][0] + y * M[1][1] + y * M[1][2];
        double tz = z * M[2][0] + z * M[2][1] + z * M[2][2];
        acc_c[0] += cos(tx);
        acc_s[0] += sin(tx);
        acc_c[1] += cos(ty);
        acc_s[1] += sin(ty);
        acc_c[2] += cos(tz);
        acc_s[2] += sin(tz);
    }

    const double inv_w = 1.0 / (double)count;
    for (int j = 0; j < 3; ++j) {
        double theta = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta += atan2(-y, -x);
        }
        // This is essentially a matrix vector multiplication, but for the single row
        out_com[j] = (float)(theta * I[j][0] + theta * I[j][1] + theta * I[j][2]);
    }
}

static void _com_pbc_iw(float out_com[3], const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(in_w);
    ASSERT(M);
    ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};
    double acc_w = 0;

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    md_512 v_acc_w = _mm512_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512i idx = _mm512_loadu_si512(in_idx + i);
        md_512 v_x  = _mm512_i32gather_ps(idx, in_x, 4);
        md_512 v_y  = _mm512_i32gather_ps(idx, in_y, 4);
        md_512 v_z  = _mm512_i32gather_ps(idx, in_z, 4);
        md_512 v_w  = _mm512_i32gather_ps(idx, in_w, 4);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_mul_ps(v_x, _mm512_set1_ps(M[0][2]))));
        md_512 v_ty = _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_y, _mm512_set1_ps(M[1][2]))));
        md_512 v_tz = _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][0]), _mm512_fmadd_ps(v_z, _mm512_set1_ps(M[2][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = _mm512_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = _mm512_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = _mm512_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = _mm512_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = _mm512_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = _mm512_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
    acc_w    = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    md_256 v_acc_w = md_mm256_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256i idx = md_mm256_loadu_si256(in_idx + i);
        md_256 v_x  = md_mm256_i32gather_ps(in_x, idx, 4);
        md_256 v_y  = md_mm256_i32gather_ps(in_y, idx, 4);
        md_256 v_z  = md_mm256_i32gather_ps(in_z, idx, 4);
        md_256 v_w  = md_mm256_i32gather_ps(in_w, idx, 4);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_mul_ps(v_x, md_mm256_set1_ps(M[0][2]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_y, md_mm256_set1_ps(M[1][2]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][0]), md_mm256_fmadd_ps(v_z, md_mm256_set1_ps(M[2][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm256_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm256_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm256_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm256_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm256_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm256_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm256_reduce_add_ps(v_acc_w);
#elif defined(__SSE2__)
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128i idx = md_mm_loadu_si128(in_idx + i);
        md_128 v_x  = md_mm_i32gather_ps(in_x, idx, 4);
        md_128 v_y  = md_mm_i32gather_ps(in_y, idx, 4);
        md_128 v_z  = md_mm_i32gather_ps(in_z, idx, 4);
        md_128 v_w  = md_mm_i32gather_ps(in_w, idx, 4);
        // Compute thetas 
        // In the non orthogonal case, this corresponds to 'multiplying' by the inverse of the Matrix
        // This is achieved by performing a dot product with the corresponding row of the inverse matrix
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_mul_ps(v_x, md_mm_set1_ps(M[0][2]))));
        md_128 v_ty = md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_y, md_mm_set1_ps(M[1][2]))));
        md_128 v_tz = md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][0]), md_mm_fmadd_ps(v_z, md_mm_set1_ps(M[2][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm_reduce_add_ps(v_acc_w);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        int32_t idx = in_idx[i];
        double x = in_x[idx];
        double y = in_y[idx];
        double z = in_z[idx];
        double w = in_w[idx];
        double tx = x * M[0][0] + x * M[0][1] + x * M[0][2];
        double ty = y * M[1][0] + y * M[1][1] + y * M[1][2];
        double tz = z * M[2][0] + z * M[2][1] + z * M[2][2];
        acc_c[0] += w * cos(tx);
        acc_s[0] += w * sin(tx);
        acc_c[1] += w * cos(ty);
        acc_s[1] += w * sin(ty);
        acc_c[2] += w * cos(tz);
        acc_s[2] += w * sin(tz);
        acc_w    += w;
    }

    const double inv_w = 1.0 / acc_w;
    for (int j = 0; j < 3; ++j) {
        double theta = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            theta += atan2(-y, -x);
        }
        // This is essentially a matrix vector multiplication, but for the single row
        out_com[j] = (float)(theta * I[j][0] + theta * I[j][1] + theta * I[j][2]);
    }
}

// This is uses the trigonometric com algorithm presented in: INSERT REF HERE
static void com_pbc(float* out_com, const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell) {

    // Should transform cartesian coordinates to fractional coordinates and scale by 2 pi
    mat3_t M = mat3_mul(mat3_scale(TWO_PI, TWO_PI, TWO_PI), unit_cell->inv_basis);
    // Inverse transform
    mat3_t I = mat3_mul(unit_cell->basis, mat3_scale(1.0/TWO_PI, 1.0/TWO_PI, 1.0/TWO_PI));

    // Here we select the correct internal version based on the available inputs
    if (in_w) {
        if (in_idx) {
            _com_pbc_iw(out_com, in_x, in_y, in_z, in_w, in_idx, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        } else {
            _com_pbc_w(out_com, in_x, in_y, in_z, in_w, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        }
    } else {
        if (in_idx) {
            _com_pbc_i(out_com, in_x, in_y, in_z, in_idx, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        } else {
            _com_pbc(out_com, in_x, in_y, in_z, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        }
    }
}

static void com_vec4(float* out_com, const vec4_t* in_xyzw, const int32_t* in_idx, size_t count) {
    // Use vec4 here so we can utilize SSE vectorization if applicable
    // @TODO: Vectorize with full register width
    vec4_t acc = {0,0,0,0};
    for (size_t i = 0; i < count; ++i) {
        int64_t idx = in_idx ? in_idx[i] : (int64_t)i;
        vec4_t xyzw = in_xyzw[idx];
        vec4_t www1 = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
        acc = vec4_add(acc, vec4_mul(xyzw, www1));
    }

    acc = vec4_div(acc, vec4_splat_w(acc));
    MEMCPY(out_com, &acc, sizeof(float) * 3);
}

static void com_pbc_vec4(float* out_com, const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell) {
    ASSERT(out_com);
    ASSERT(in_xyzw);
    ASSERT(unit_cell);

    vec4_t acc_c = {0};
    vec4_t acc_s = {0};
    vec4_t acc_xyzw = {0};

    if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        const vec3_t ext = mat3_diag(unit_cell->basis);
        const vec4_t scl = vec4_div(vec4_set1(TWO_PI), vec4_from_vec3(ext, TWO_PI));

        if (in_idx) {
            for (size_t i = 0; i < count; ++i) {
                int32_t idx  = in_idx[i];
                vec4_t xyzw  = in_xyzw[idx];
                vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                vec4_t theta = vec4_mul(xyzw, scl);
                vec4_t s,c;
                vec4_sincos(theta, &s, &c);
                acc_s = vec4_add(acc_s, vec4_mul(s, www1));
                acc_c = vec4_add(acc_c, vec4_mul(c, www1));
                acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
            }
        } else {
            for (size_t i = 0; i < count; ++i) {
                vec4_t xyzw  = in_xyzw[i];
                vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                vec4_t theta = vec4_mul(xyzw, scl);
                vec4_t s,c;
                vec4_sincos(theta, &s, &c);
                acc_s = vec4_add(acc_s, vec4_mul(s, www1));
                acc_c = vec4_add(acc_c, vec4_mul(c, www1));
                acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
            }
        }

        const float w = acc_xyzw.w;
        for (int i = 0; i < 3; ++i) {
            const double y = acc_s.elem[i] / w;
            const double x = acc_c.elem[i] / w;
            const double r2 = x*x + y*y;
            double theta_prim = PI;
            if (r2 > 1.0e-15) {
                theta_prim += atan2(-y, -x);
            }
            out_com[i] = (float)((theta_prim / TWO_PI) * ext.elem[i]);
        }
    } else {
        ASSERT(unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC);
        const mat4_t M = mat4_mul(mat4_scale(TWO_PI, TWO_PI, TWO_PI), mat4_from_mat3(unit_cell->inv_basis));
        const mat4_t I = mat4_mul(mat4_scale(1.0f / TWO_PI, 1.0f / TWO_PI, 1.0f / TWO_PI), mat4_from_mat3(unit_cell->basis));

        if (in_idx) {
            for (size_t i = 0; i < count; ++i) {
                int32_t idx  = in_idx[i];
                vec4_t xyzw  = in_xyzw[idx];
                vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                vec4_t theta = mat4_mul_vec4(M, xyzw);
                vec4_t s,c;
                vec4_sincos(theta, &s, &c);
                acc_s = vec4_add(acc_s, vec4_mul(s, www1));
                acc_c = vec4_add(acc_c, vec4_mul(c, www1));
                acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
            }
        } else {
            for (size_t i = 0; i < count; ++i) {
                vec4_t xyzw  = in_xyzw[i];
                vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                vec4_t theta = mat4_mul_vec4(I, xyzw);
                vec4_t s,c;
                vec4_sincos(theta, &s, &c);
                acc_s = vec4_add(acc_s, vec4_mul(s, www1));
                acc_c = vec4_add(acc_c, vec4_mul(c, www1));
                acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
            }
        }

        for (int i = 0; i < 3; ++i) {
            const double y = acc_s.elem[i] / acc_xyzw.w;
            const double x = acc_c.elem[i] / acc_xyzw.w;
            const double r2 = x*x + y*y;
            double theta_prim = PI;
            if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
                theta_prim += atan2(-y, -x);
            }
            out_com[i] = (float)(theta_prim * I.elem[i][0] + theta_prim * I.elem[i][1] + theta_prim * I.elem[i][2]);
        }
    }
}

vec3_t md_util_com_compute(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    vec3_t xyz = {0};
    if (!unit_cell || unit_cell->flags == MD_UNIT_CELL_FLAG_NONE) {
        com(xyz.elem, in_x, in_y, in_z, in_w, in_idx, count);
        return xyz;
    }

    if (unit_cell->flags & (MD_UNIT_CELL_FLAG_ORTHO | MD_UNIT_CELL_FLAG_TRICLINIC)) {
        com_pbc(xyz.elem, in_x, in_y, in_z, in_w, in_idx, count, unit_cell);
        return xyz;
    }

    // Error
    MD_LOG_ERROR("Invalid unit cell flags: %d", unit_cell->flags);
    return xyz;
}

vec3_t md_util_com_compute_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell) {
    ASSERT(in_xyzw);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    vec3_t xyz = {0};
    if (unit_cell && (unit_cell->flags & (MD_UNIT_CELL_FLAG_ORTHO | MD_UNIT_CELL_FLAG_TRICLINIC))) {
        com_pbc_vec4(xyz.elem, in_xyzw, in_idx, count, unit_cell);
    } else {
        com_vec4(xyz.elem, in_xyzw, in_idx, count);
    }

    return xyz;
}

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

//Triclinic
mat3_t md_util_compute_triclinic_unit_cell_basis(double x, double y, double z, double xy, double xz, double yz) {
    mat3_t M = {
        .col = {
            {x, 0, 0},
            {xy, y, 0},
            {xz, yz, z},
        },
    };
    return M;
}

md_unit_cell_t md_util_unit_cell_from_triclinic(double x, double y, double z, double xy, double xz, double yz) {
    mat3_t matrix = md_util_compute_triclinic_unit_cell_basis(x, y, z, xy, xz, yz);
    return md_util_unit_cell_from_matrix(matrix.elem);
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
        .flags = MD_UNIT_CELL_FLAG_ORTHO | (x != 0.0 ? MD_UNIT_CELL_FLAG_PBC_X : 0) | (y != 0.0 ? MD_UNIT_CELL_FLAG_PBC_Y : 0) | (z != 0.0 ? MD_UNIT_CELL_FLAG_PBC_Z : 0),
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
        .flags = MD_UNIT_CELL_FLAG_TRICLINIC | (a != 0.0 ? MD_UNIT_CELL_FLAG_PBC_X : 0) | (b != 0.0 ? MD_UNIT_CELL_FLAG_PBC_Y : 0) | (c != 0.0 ? MD_UNIT_CELL_FLAG_PBC_Z : 0),
    };

    return cell;
}

md_unit_cell_t md_util_unit_cell_from_matrix(float M[3][3]) {
    if (M[0][1] == 0 && M[0][2] == 0 && M[1][0] == 0 && M[1][2] == 0 && M[2][0] == 0 && M[2][1] == 0) {
        return md_util_unit_cell_from_extent(M[0][0], M[1][1], M[2][2]);
    } else {
        const float* N = &M[0][0];
        
        // Inverse of M
        const double I[9] = {
            1.0/(double)N[0], 0, 0,
            -(double)N[3]/((double)N[0]*(double)N[4]), 1.0/(double)N[4], 0,
            ((double)N[3]*(double)N[7]/((double)N[0]*(double)N[4]) - (double)N[6]/(double)N[0])/(double)N[8], -(double)N[7]/((double)N[4]*(double)N[8]), 1.0/(double)N[8],
        };
        

        const float a = M[0][0]*M[0][0] + M[0][1]*M[0][1] + M[0][2]*M[0][2];
        const float b = M[1][0]*M[1][0] + M[1][1]*M[1][1] + M[1][2]*M[1][2];
        const float c = M[2][0]*M[2][0] + M[2][1]*M[2][1] + M[2][2]*M[2][2];
        
        
        mat3_t basis;
        MEMCPY(basis.elem, M, sizeof(basis));
        md_unit_cell_t cell = {
            .basis = basis,
            .inv_basis = {
                I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8],
            },
            .flags = MD_UNIT_CELL_FLAG_TRICLINIC | (a != 0 ? MD_UNIT_CELL_FLAG_PBC_X : 0) | (b != 0 ? MD_UNIT_CELL_FLAG_PBC_Y : 0) | (c != 0 ? MD_UNIT_CELL_FLAG_PBC_Z : 0),
        };
        return cell;
    }
}

void md_util_unit_cell_distance_array(float* out_dist, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unit_cell_t* cell) {
    if (cell->flags == 0) {
        for (size_t i = 0; i < num_a; ++i) {
            for (size_t j = 0; j < num_b; ++j) {
                out_dist[i * num_b + j] = vec3_distance(coord_a[i], coord_b[j]);
            }
        }
    }
    else if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        const vec4_t box = {cell->basis.elem[0][0], cell->basis.elem[1][1], cell->basis.elem[2][2], 0};
        for (size_t i = 0; i < num_a; ++i) {
            for (size_t j = 0; j < num_b; ++j) {
                vec4_t a = vec4_from_vec3(coord_a[i], 0);
                vec4_t b = vec4_from_vec3(coord_b[j], 0);
                out_dist[i * num_b + j] = vec4_periodic_distance(a, b, box);
            }
        }
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (size_t i = 0; i < num_a; ++i) {
            for (size_t j = 0; j < num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, cell->basis.elem);
                out_dist[i * num_b + j] = vec3_length(dx);
            }
        }
    }
}

float md_util_unit_cell_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unit_cell_t* cell) {
    int64_t min_i = 0;
    int64_t min_j = 0;
    float min_dist = FLT_MAX;

    if (cell->flags == 0) {
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
                const float d = vec3_distance(coord_a[i], coord_b[j]);
                if (d < min_dist) {
                    min_dist = d;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }
    else if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        const vec4_t box = {cell->basis.elem[0][0], cell->basis.elem[1][1], cell->basis.elem[2][2], 0};
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
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
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
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

float md_util_unit_cell_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unit_cell_t* cell) {
    int64_t max_i = 0;
    int64_t max_j = 0;
    float max_dist = 0;

    if (cell->flags == 0) {
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
                const float d = vec3_distance(coord_a[i], coord_b[j]);
                if (d > max_dist) {
                    max_dist = d;
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }
    else if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        const vec4_t box = {cell->basis.elem[0][0], cell->basis.elem[1][1], cell->basis.elem[2][2], 0};
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
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
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
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

// Blatantly stolen from MDAnalysis project
// MDAnalysis --- https://www.mdanalysis.org
// Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
// https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/include/calc_distances.h

static inline void min_image_triclinic(float dx[3], const float box[3][3], const float half_diag[3]) {
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

    // Ensure that dx is within the required "zone"
    for (int i = 2; i >= 0; i--) {
        if (half_diag[i] > 0.0f) {
            while (dx[i] > half_diag[i]) {
                for (int j = i; j >= 0; j--) {
                    dx[j] -= box[i][j];
                }
            }
            while (dx[i] <= -half_diag[i]) {
                for (int j = i; j >= 0; j--) {
                    dx[j] += box[i][j];
                }
            }
        }
    }

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

static inline void min_image_ortho(float dx[3], float ext[3], float half_ext[3]) {
    for (int i = 0; i < 3; i++) {
        if (ext[i] > 0.0f) {
            while (dx[i] > half_ext[i]) {
                dx[i] -= ext[i];
            }
            while (dx[i] <= -half_ext[i]) {
                dx[i] += ext[i];
            }
        }
    }
}

void md_util_min_image_vec3(vec3_t dx[], size_t count, const md_unit_cell_t* unit_cell) {
    if (unit_cell) {
        vec3_t diag = mat3_diag(unit_cell->basis);
        vec3_t half_diag = vec3_mul_f(diag, 0.5f);
        if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
            for (size_t i = 0; i < count; ++i) {
                min_image_ortho(dx[i].elem, diag.elem, half_diag.elem);
            }
        } else if (unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
            for (size_t i = 0; i < count; ++i) {
                min_image_triclinic(dx[i].elem, unit_cell->basis.elem, half_diag.elem);
            }
        }
    }
}

void md_util_min_image_vec4(vec4_t dx[], size_t count, const md_unit_cell_t* unit_cell) {
    if (unit_cell) {
        vec3_t diag = mat3_diag(unit_cell->basis);
        vec3_t half_diag = vec3_mul_f(diag, 0.5f);
        if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
            for (size_t i = 0; i < count; ++i) {
                min_image_ortho(dx[i].elem, diag.elem, half_diag.elem);
            }
        } else if (unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
            for (size_t i = 0; i < count; ++i) {
                min_image_triclinic(dx[i].elem, unit_cell->basis.elem, half_diag.elem);
            }
        }
    }
}

static void pbc_ortho(float* x, float* y, float* z, const int32_t* indices, size_t count, vec3_t box_ext) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);

    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    if (indices) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx = indices[i];
            vec4_t pos = {x[idx], y[idx], z[idx], 0};
            pos = vec4_deperiodize_ortho(pos, ref, ext);
            x[i] = pos.x;
            y[i] = pos.y;
            z[i] = pos.z;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            vec4_t pos = {x[i], y[i], z[i], 0};
            pos = vec4_deperiodize_ortho(pos, ref, ext);
            x[i] = pos.x;
            y[i] = pos.y;
            z[i] = pos.z;
        }
    }
}

static void pbc_ortho_vec4(vec4_t* xyzw, size_t count, vec3_t box_ext) {
    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    const vec4_t ref = vec4_mul_f(ext, 0.5f);
    for (size_t i = 0; i < count; ++i) {
        xyzw[i] = vec4_deperiodize_ortho(xyzw[i], ref, ext);
    }
}

static void pbc_triclinic(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unit_cell_t* cell) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);

    const vec4_t mask = md_unit_cell_pbc_mask(cell);
    mat4x3_t I = mat4x3_from_mat3(cell->inv_basis);
    mat4x3_t M = mat4x3_from_mat3(cell->basis);

    if (indices) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx = indices[i];
            vec4_t original = {x[idx], y[idx], z[idx], 0};
            vec4_t coord = mat4x3_mul_vec4(I, original);
            coord = vec4_fract(coord);
            coord = mat4x3_mul_vec4(M, coord);
            coord = vec4_blend(original, coord, mask);

            x[idx] = coord.x;
            y[idx] = coord.y;
            z[idx] = coord.z;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            vec4_t original = {x[i], y[i], z[i], 0};
            vec4_t coord = mat4x3_mul_vec4(I, original);
            coord = vec4_fract(coord);
            coord = mat4x3_mul_vec4(M, coord);
            coord = vec4_blend(original, coord, mask);

            x[i] = coord.x;
            y[i] = coord.y;
            z[i] = coord.z;
        }
    }
}

static void pbc_triclinic_vec4(vec4_t* xyzw, size_t count, const md_unit_cell_t* cell) {
    const vec4_t mask = md_unit_cell_pbc_mask(cell);
    mat4x3_t I = mat4x3_from_mat3(cell->inv_basis);
    mat4x3_t M = mat4x3_from_mat3(cell->basis);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = xyzw[i];
        c = mat4x3_mul_vec4(I, c);
        c = vec4_fract(c);
        c = mat4x3_mul_vec4(M, c);
        xyzw[i] = vec4_blend(xyzw[i], c, mask);
    }
}

bool md_util_pbc(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unit_cell_t* unit_cell) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (!unit_cell) {
        MD_LOG_ERROR("Missing unit cell");
        return false;
    }

    if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        vec3_t ext = mat3_diag(unit_cell->basis);
        pbc_ortho(x, y, z, indices, count, ext);
        return true;
    } else if (unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        pbc_triclinic(x, y, z, indices, count, unit_cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    MD_LOG_ERROR("Unrecognized unit_cell type");
    return false;
}

bool md_util_pbc_vec4(vec4_t* in_out_xyzw, size_t count, const md_unit_cell_t* unit_cell) {
    if (!in_out_xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (!unit_cell) {
        MD_LOG_ERROR("Missing unit cell");
        return false;
    }

    if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        vec3_t ext = mat3_diag(unit_cell->basis);
        pbc_ortho_vec4(in_out_xyzw, count, ext);
        return true;
    } else if (unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        pbc_triclinic_vec4(in_out_xyzw, count, unit_cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    MD_LOG_ERROR("Unrecognized unit_cell type");
    return false;
}

static void unwrap_ortho(float* x, float* y, float* z, const int32_t* indices, size_t count, vec3_t box_ext) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);

    const vec4_t ext = vec4_from_vec3(box_ext, 0);

    if (indices) {
        int idx = indices[0];
        vec4_t ref_pos = vec4_set(x[idx], y[idx], z[idx], 0);
        for (size_t i = 1; i < count; ++i) {
            idx = indices[i];
            const vec4_t pos = vec4_deperiodize_ortho((vec4_t){x[idx], y[idx], z[idx], 0}, ref_pos, ext);
            x[idx] = pos.x;
            y[idx] = pos.y;
            z[idx] = pos.z;
            ref_pos = pos;
        }
    } else {
        vec4_t ref_pos = {x[0], y[0], z[0], 0};
        for (size_t i = 1; i < count; ++i) {
            const vec4_t pos = vec4_deperiodize_ortho((vec4_t){x[i], y[i], z[i], 0}, ref_pos, ext);
            x[i] = pos.x;
            y[i] = pos.y;
            z[i] = pos.z;
            ref_pos = pos;
        }
    }
}

static void unwrap_ortho_vec4(vec4_t* xyzw, size_t count, vec3_t box_ext) {
    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    vec4_t ref_pos = xyzw[0];
    for (size_t i = 1; i < count; ++i) {
        const vec4_t pos = vec4_deperiodize_ortho(xyzw[i], ref_pos, ext);
        xyzw[i] = pos;
        ref_pos = pos;
    }
}

static void unwrap_triclinic(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unit_cell_t* cell) {
    if (indices) {
        int idx = indices[0];
        vec3_t ref_pos = {x[idx], y[idx], z[idx]};
        for (size_t i = 1; i < count; ++i) {
            idx = indices[i];
            vec3_t pos = {x[idx], y[idx], z[idx]};
            deperiodize_triclinic(pos.elem, ref_pos.elem, cell->basis.elem);
            x[idx] = pos.x;
            y[idx] = pos.y;
            z[idx] = pos.z;
            ref_pos = pos;
        }
    } else {
        vec3_t ref_pos = {x[0], y[0], z[0]};
        for (size_t i = 1; i < count; ++i) {
            vec3_t pos = {x[i], y[i], z[i]};
            deperiodize_triclinic(pos.elem, ref_pos.elem, cell->basis.elem);
            x[i] = pos.x;
            y[i] = pos.y;
            z[i] = pos.z;
            ref_pos = pos;
        }
    }
}

static void unwrap_triclinic_vec4(vec4_t* xyzw, size_t count, const md_unit_cell_t* cell) {
    vec4_t ref_pos = xyzw[0];
    for (size_t i = 1; i < count; ++i) {
        vec4_t pos = xyzw[i];
        deperiodize_triclinic(pos.elem, ref_pos.elem, cell->basis.elem);
        xyzw[i] = pos;
        ref_pos = pos;
    }

}

bool md_util_unwrap(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unit_cell_t* cell) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing cell");
        return false;
    }

    if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        vec3_t ext = mat3_diag(cell->basis);
        unwrap_ortho(x, y, z, indices, count, ext);
        return true;
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        unwrap_triclinic(x, y, z, indices, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    return false;
}

bool md_util_unwrap_vec4(vec4_t* xyzw, size_t count, const md_unit_cell_t* cell) {
    if (!xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing cell");
        return false;
    }

    if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        vec3_t ext = mat3_diag(cell->basis);
        unwrap_ortho_vec4(xyzw, count, ext);
        return true;
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        unwrap_triclinic_vec4(xyzw, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    return false;
}

bool md_util_deperiodize_vec4(vec4_t* xyzw, size_t count, vec3_t ref_xyz, const md_unit_cell_t* cell) {
    if (!xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing cell");
        return false;
    }

    if (cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
        vec4_t ref_pos = vec4_from_vec3(ref_xyz, 0);
        vec4_t ext = vec4_from_vec3(mat3_diag(cell->basis), 0);
        for (size_t i = 0; i < count; ++i) {
            const vec4_t pos = vec4_deperiodize_ortho(xyzw[i], ref_pos, ext);
            xyzw[i] = pos;
        }
        return true;
    } else if (cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
        for (size_t i = 1; i < count; ++i) {
            vec4_t pos = xyzw[i];
            deperiodize_triclinic(pos.elem, ref_xyz.elem, cell->basis.elem);
            xyzw[i] = pos;
        }
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    return false;
}

double md_util_rmsd_compute(const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], const float* const in_w[2], const int32_t* const in_idx[2], const size_t count, const vec3_t in_com[2]) {
    const mat3_t R = mat3_optimal_rotation(in_x, in_y, in_z, in_w, in_idx, count, in_com);
    double d_sum = 0;
    double w_sum = 0;
    if (in_idx) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx0 = in_idx[0][i];
            int32_t idx1 = in_idx[1][i];
            vec3_t u = {in_x[0][idx0] - in_com[0].x, in_y[0][idx0] - in_com[0].y, in_z[0][idx0] - in_com[0].z};
            vec3_t v = {in_x[1][idx1] - in_com[1].x, in_y[1][idx1] - in_com[1].y, in_z[1][idx1] - in_com[1].z};
            vec3_t vp = mat3_mul_vec3(R, v);
            vec3_t d = vec3_sub(u, vp);
            float weight = in_w ? (in_w[0][idx0] + in_w[1][idx1]) * 0.5f : 1.0f;
            d_sum += weight * vec3_dot(d, d);
            w_sum += weight;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            vec3_t u = {in_x[0][i] - in_com[0].x, in_y[0][i] - in_com[0].y, in_z[0][i] - in_com[0].z};
            vec3_t v = {in_x[1][i] - in_com[1].x, in_y[1][i] - in_com[1].y, in_z[1][i] - in_com[1].z};
            vec3_t vp = mat3_mul_vec3(R, v);
            vec3_t d = vec3_sub(u, vp);
            float weight = in_w ? (in_w[0][i] + in_w[1][i]) * 0.5f : 1.0f;
            d_sum += weight * vec3_dot(d, d);
            w_sum += weight;
        }
    }

    return sqrt(d_sum / w_sum);
}

double md_util_rmsd_compute_vec4(const vec4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]) {
    const mat3_t R = mat3_optimal_rotation_vec4(in_xyzw, in_idx, count, in_com);
    double d_sum = 0;
    double w_sum = 0;
    if (in_idx) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx0 = in_idx[0][i];
            int32_t idx1 = in_idx[1][i];
            vec3_t u = {in_xyzw[0][idx0].x - in_com[0].x, in_xyzw[0][idx0].y - in_com[0].y, in_xyzw[0][idx0].z - in_com[0].z};
            vec3_t v = {in_xyzw[1][idx1].x - in_com[1].x, in_xyzw[1][idx1].y - in_com[1].y, in_xyzw[1][idx1].z - in_com[1].z};
            vec3_t vp = mat3_mul_vec3(R, v);
            vec3_t d = vec3_sub(u, vp);
            float w = (in_xyzw[0][idx0].w + in_xyzw[1][idx1].w) * 0.5f;
            d_sum += w * vec3_dot(d, d);
            w_sum += w;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            vec3_t u = {in_xyzw[0][i].x - in_com[0].x, in_xyzw[0][i].y - in_com[0].y, in_xyzw[0][i].z - in_com[0].z};
            vec3_t v = {in_xyzw[1][i].x - in_com[1].x, in_xyzw[1][i].y - in_com[1].y, in_xyzw[1][i].z - in_com[1].z};
            vec3_t vp = mat3_mul_vec3(R, v);
            vec3_t d = vec3_sub(u, vp);
            float w = (in_xyzw[0][i].w + in_xyzw[1][i].w) * 0.5f;
            d_sum += w * vec3_dot(d, d);
            w_sum += w;
        }
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

bool md_util_interpolate_linear(float* out_x, float* out_y, float* out_z, const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], size_t count, const md_unit_cell_t* unit_cell, float t) {
    ASSERT(out_x);
    ASSERT(out_y);
    ASSERT(out_z);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    t = CLAMP(t, 0.0f, 1.0f);

    if (unit_cell && unit_cell->flags & (MD_UNIT_CELL_FLAG_ORTHO | MD_UNIT_CELL_FLAG_TRICLINIC)) {
        for (size_t i = 0; i < count; i += 8) {
            md_256 x0[3] = {
                md_mm256_loadu_ps(in_x[0] + i),
                md_mm256_loadu_ps(in_y[0] + i),
                md_mm256_loadu_ps(in_z[0] + i)
            };

            md_256 x1[3] = {
                md_mm256_loadu_ps(in_x[1] + i),
                md_mm256_loadu_ps(in_y[1] + i),
                md_mm256_loadu_ps(in_z[1] + i)
            };

            if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
                simd_deperiodize_ortho(x1[0], x0[0], md_mm256_set1_ps(unit_cell->basis.elem[0][0]), md_mm256_set1_ps(unit_cell->inv_basis.elem[0][0]));
                simd_deperiodize_ortho(x1[1], x0[1], md_mm256_set1_ps(unit_cell->basis.elem[1][1]), md_mm256_set1_ps(unit_cell->inv_basis.elem[1][1]));
                simd_deperiodize_ortho(x1[2], x0[2], md_mm256_set1_ps(unit_cell->basis.elem[2][2]), md_mm256_set1_ps(unit_cell->inv_basis.elem[2][2]));
            } else if (unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
                simd_deperiodize_triclinic(x1, x0, unit_cell->basis.elem);
            }

            md_256 x = md_mm256_lerp_ps(x0[0], x1[0], t);
            md_256 y = md_mm256_lerp_ps(x0[1], x1[1], t);
            md_256 z = md_mm256_lerp_ps(x0[2], x1[2], t);

            md_mm256_storeu_ps(out_x + i, x);
            md_mm256_storeu_ps(out_y + i, y);
            md_mm256_storeu_ps(out_z + i, z);
        }
    } else {
        // Straight forward lerp
        for (size_t i = 0; i < count; i += 8) {
            md_256 x0 = md_mm256_loadu_ps(in_x[0] + i);
            md_256 y0 = md_mm256_loadu_ps(in_y[0] + i);
            md_256 z0 = md_mm256_loadu_ps(in_z[0] + i);

            md_256 x1 = md_mm256_loadu_ps(in_x[1] + i);
            md_256 y1 = md_mm256_loadu_ps(in_y[1] + i);
            md_256 z1 = md_mm256_loadu_ps(in_z[1] + i);

            md_256 x = md_mm256_lerp_ps(x0, x1, t);
            md_256 y = md_mm256_lerp_ps(y0, y1, t);
            md_256 z = md_mm256_lerp_ps(z0, z1, t);

            md_mm256_storeu_ps(out_x + i, x);
            md_mm256_storeu_ps(out_y + i, y);
            md_mm256_storeu_ps(out_z + i, z);
        }
    }

    return true;
}

bool md_util_interpolate_cubic_spline(float* out_x, float* out_y, float* out_z, const float* const in_x[4], const float* const in_y[4], const float* const in_z[4], size_t count, const md_unit_cell_t* unit_cell, float t, float s) {
    ASSERT(out_x);
    ASSERT(out_y);
    ASSERT(out_z);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    
    t = CLAMP(t, 0.0f, 1.0f);
    s = CLAMP(s, 0.0f, 1.0f);

    for (size_t i = 0; i < count; i += 8) {
        md_256 x0[3] = {
            md_mm256_loadu_ps(in_x[0] + i),
            md_mm256_loadu_ps(in_y[0] + i),
            md_mm256_loadu_ps(in_z[0] + i),
        };

        md_256 x1[3] = {
            md_mm256_loadu_ps(in_x[1] + i),
            md_mm256_loadu_ps(in_y[1] + i),
            md_mm256_loadu_ps(in_z[1] + i),
        };

        md_256 x2[3] = {
            md_mm256_loadu_ps(in_x[2] + i),
            md_mm256_loadu_ps(in_y[2] + i),
            md_mm256_loadu_ps(in_z[2] + i),
        };

        md_256 x3[3] = {
            md_mm256_loadu_ps(in_x[3] + i),
            md_mm256_loadu_ps(in_y[3] + i),
            md_mm256_loadu_ps(in_z[3] + i),
        };

        if (unit_cell->flags & MD_UNIT_CELL_FLAG_ORTHO) {
            x0[0] = simd_deperiodize_ortho(x0[0], x1[0], md_mm256_set1_ps(unit_cell->basis.elem[0][0]), md_mm256_set1_ps(unit_cell->inv_basis.elem[0][0]));
            x2[0] = simd_deperiodize_ortho(x2[0], x1[0], md_mm256_set1_ps(unit_cell->basis.elem[0][0]), md_mm256_set1_ps(unit_cell->inv_basis.elem[0][0]));
            x3[0] = simd_deperiodize_ortho(x3[0], x2[0], md_mm256_set1_ps(unit_cell->basis.elem[0][0]), md_mm256_set1_ps(unit_cell->inv_basis.elem[0][0]));
            x0[1] = simd_deperiodize_ortho(x0[1], x1[1], md_mm256_set1_ps(unit_cell->basis.elem[1][1]), md_mm256_set1_ps(unit_cell->inv_basis.elem[1][1]));
            x2[1] = simd_deperiodize_ortho(x2[1], x1[1], md_mm256_set1_ps(unit_cell->basis.elem[1][1]), md_mm256_set1_ps(unit_cell->inv_basis.elem[1][1]));
            x3[1] = simd_deperiodize_ortho(x3[1], x2[1], md_mm256_set1_ps(unit_cell->basis.elem[1][1]), md_mm256_set1_ps(unit_cell->inv_basis.elem[1][1]));
            x0[2] = simd_deperiodize_ortho(x0[2], x1[2], md_mm256_set1_ps(unit_cell->basis.elem[2][2]), md_mm256_set1_ps(unit_cell->inv_basis.elem[2][2]));
            x2[2] = simd_deperiodize_ortho(x2[2], x1[2], md_mm256_set1_ps(unit_cell->basis.elem[2][2]), md_mm256_set1_ps(unit_cell->inv_basis.elem[2][2]));
            x3[2] = simd_deperiodize_ortho(x3[2], x2[2], md_mm256_set1_ps(unit_cell->basis.elem[2][2]), md_mm256_set1_ps(unit_cell->inv_basis.elem[2][2]));
        } else if (unit_cell->flags & MD_UNIT_CELL_FLAG_TRICLINIC) {
            simd_deperiodize_triclinic(x0, x1, unit_cell->basis.elem);
            simd_deperiodize_triclinic(x2, x1, unit_cell->basis.elem);
            simd_deperiodize_triclinic(x3, x2, unit_cell->basis.elem);
        }

        const md_256 x = md_mm256_cubic_spline_ps(x0[0], x1[0], x2[0], x3[0], md_mm256_set1_ps(t), md_mm256_set1_ps(s));
        const md_256 y = md_mm256_cubic_spline_ps(x0[1], x1[1], x2[1], x3[1], md_mm256_set1_ps(t), md_mm256_set1_ps(s));
        const md_256 z = md_mm256_cubic_spline_ps(x0[2], x1[2], x2[2], x3[2], md_mm256_set1_ps(t), md_mm256_set1_ps(s));

        md_mm256_storeu_ps(out_x + i, x);
        md_mm256_storeu_ps(out_y + i, y);
        md_mm256_storeu_ps(out_z + i, z);
    }

    return true;
}

static inline bool ranges_overlap(md_range_t a, md_range_t b) {
    return (a.beg < b.end && b.beg < a.end);
}

static inline void commit_protein_backbone(md_protein_backbone_atoms_t* bb_atoms, size_t bb_length, size_t res_offset, md_chain_idx_t chain_idx, md_molecule_t* mol, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(mol->protein_backbone.atoms);
    md_array_push_array(mol->protein_backbone.atoms, bb_atoms, bb_length, alloc);
    md_array_push(mol->protein_backbone.range.offset, offset, alloc);
    md_array_push(mol->protein_backbone.range.chain_idx, chain_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        int32_t res_idx = (int32_t)(res_offset + i);
        md_array_push(mol->protein_backbone.residue_idx, res_idx, alloc);
    }
}

static inline void commit_nucleic_backbone(md_nucleic_backbone_atoms_t* bb_atoms, size_t bb_length, size_t res_offset, md_chain_idx_t chain_idx, md_molecule_t* mol, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(mol->nucleic_backbone.atoms);
    md_array_push_array(mol->nucleic_backbone.atoms, bb_atoms, bb_length, alloc);
    md_array_push(mol->nucleic_backbone.range.offset, offset, alloc);
    md_array_push(mol->nucleic_backbone.range.chain_idx, chain_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        int32_t res_idx = (int32_t)(res_offset + i);
        md_array_push(mol->nucleic_backbone.residue_idx, res_idx, alloc);
    }
}

#if 0
static void md_util_compute_backbone_data(md_protein_backbone_data_t* protein_backbone, const md_molecule_t* mol, md_allocator_i* alloc) {
    if (mol->chain.count) {
        // Compute backbone data
        // 
        // @NOTE: We should only attempt to compute backbone data for valid residues (e.g. amino acids / dna)
        // Backbones are not directly tied to chains and therefore we cannot use chains as a 1:1 mapping for the backbones.
        // We look within the chains and see if we can find consecutive ranges which form backbones.

        static const int64_t MIN_BACKBONE_LENGTH = 3;
        md_protein_backbone_atoms_t* bb_atoms = 0;
        md_allocator_i* temp_alloc = md_get_temp_allocator();

        for (int64_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
            for (int64_t res_idx = mol->chain.residue_range[chain_idx].beg; res_idx < mol->chain.residue_range[chain_idx].end; ++res_idx) {
                md_protein_backbone_atoms_t atoms;
                if (md_util_backbone_atoms_extract_from_residue_idx(&atoms, (int32_t)res_idx, mol)) {
                    md_array_push(bb_atoms, atoms, temp_alloc);
                } else {
                    if (md_array_size(protein_backbone) >= MIN_BACKBONE_LENGTH) {
                        // Commit the backbone
                        md_range_t res_range = {(int32_t)(res_idx - md_array_size(bb_atoms)), (int32_t)res_idx};
                        commit_backbone(protein_backbone, bb_atoms, md_array_size(bb_atoms), res_range, alloc);
                    }
                    md_array_shrink(bb_atoms, 0);
                }
            }
            // Possibly commit remainder of the chain
            if (md_array_size(bb_atoms) >= MIN_BACKBONE_LENGTH) {
                int64_t res_idx = mol->chain.residue_range[chain_idx].end;
                md_range_t res_range = {(int32_t)(res_idx - md_array_size(bb_atoms)), (int32_t)res_idx};
                commit_backbone(protein_backbone, bb_atoms, md_array_size(bb_atoms), res_range, alloc);
            }
            md_array_shrink(bb_atoms, 0);
        }

        md_array_free(bb_atoms, temp_alloc);

        protein_backbone->range_count = md_array_size(protein_backbone->range);
        protein_backbone->count = md_array_size(protein_backbone->atoms);

        md_array_resize(protein_backbone->angle, protein_backbone->count, alloc);
        md_array_resize(protein_backbone->secondary_structure, protein_backbone->count, alloc);
        md_array_resize(protein_backbone->ramachandran_type, protein_backbone->count, alloc);

        if (mol->protein_backbone.count > 0) {
            md_util_backbone_angles_compute(protein_backbone->angle, protein_backbone->count, mol);
            md_util_secondary_structure_compute(protein_backbone->secondary_structure, protein_backbone->count, mol);
            md_util_backbone_ramachandran_classify(protein_backbone->ramachandran_type, protein_backbone->count, mol);
        }
    }
}
#endif

// Try to fill in missing fields for molecule struct
// (resid)                  -> Residues
// (Types & resname)        -> Elements
// (Elements)               -> Mass & Radius
// (Coordinates & Elements) -> Covalent Bonds
// (resid & Bonds)          -> Chains
// (resname & types)        -> Backbone
bool md_util_molecule_postprocess(md_molecule_t* mol, md_allocator_i* alloc, md_util_postprocess_flags_t flags) {
    ASSERT(mol);
    ASSERT(alloc);

    //MD_LOG_DEBUG("Starting Postprocessing Molecule");

    if (mol->atom.count == 0) return false;

    {
#ifdef PROFILE
        md_timestamp_t t0 = md_time_current();
#endif
        if (!mol->atom.flags) {
            md_array_resize(mol->atom.flags, mol->atom.count, alloc);
            MEMSET(mol->atom.flags, 0, md_array_bytes(mol->atom.flags));
        }
        if (!mol->atom.type) {
            md_array_resize(mol->atom.type, mol->atom.count, alloc);
            MEMSET(mol->atom.type, 0, md_array_bytes(mol->atom.type));
        }
#ifdef PROFILE
        md_timestamp_t t1 = md_time_current();
        MD_LOG_DEBUG("Postprocess: allocate missing fields %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
    }

    if (flags & MD_UTIL_POSTPROCESS_RESIDUE_BIT) {
        // Create residues from resids
        if (mol->residue.count == 0 && mol->atom.resid) {
#ifdef PROFILE
            md_timestamp_t t0 = md_time_current();
#endif
            md_util_compute_residue_data(&mol->residue, &mol->atom, alloc);
#ifdef PROFILE
            md_timestamp_t t1 = md_time_current();
            MD_LOG_DEBUG("Postprocess: compute residues %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_ELEMENT_BIT) {
#ifdef PROFILE
        md_timestamp_t t0 = md_time_current();
#endif
        if (!mol->atom.element) {
            md_array_resize(mol->atom.element, mol->atom.count, alloc);
            MEMSET(mol->atom.element, 0, md_array_bytes(mol->atom.element));
        }
        md_util_element_guess(mol->atom.element, mol->atom.count, mol);
#ifdef PROFILE
        md_timestamp_t t1 = md_time_current();
        MD_LOG_DEBUG("Postprocess: guess elements %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
    }

    if (flags & MD_UTIL_POSTPROCESS_RADIUS_BIT) {
#ifdef PROFILE
        md_timestamp_t t0 = md_time_current();
#endif
        if (!mol->atom.radius) {
            md_array_resize(mol->atom.radius, mol->atom.count, alloc);
            for (size_t i = 0; i < mol->atom.count; ++i) {
                mol->atom.radius[i] = mol->atom.element ? md_util_element_vdw_radius(mol->atom.element[i]) : 1.0f;
            }
        }
#ifdef PROFILE
        md_timestamp_t t1 = md_time_current();
        MD_LOG_DEBUG("Postprocess: compute radii %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
    }

    if (flags & MD_UTIL_POSTPROCESS_MASS_BIT) {
        if (!mol->atom.mass) {
#ifdef PROFILE
            md_timestamp_t t0 = md_time_current();
#endif
            md_array_resize(mol->atom.mass, mol->atom.count, alloc);
            for (size_t i = 0; i < mol->atom.count; ++i) {
                mol->atom.mass[i] = mol->atom.element ? md_util_element_atomic_mass(mol->atom.element[i]) : 1.0f;
            }
#ifdef PROFILE
            md_timestamp_t t1 = md_time_current();
            MD_LOG_DEBUG("Postprocess: compute masses %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
        }
    }
   
    if (flags & MD_UTIL_POSTPROCESS_BOND_BIT) {
#ifdef PROFILE
        md_timestamp_t t0 = md_time_current();
#endif
        md_util_covalent_bonds_compute(&mol->bond, mol, alloc);
#ifdef PROFILE
        md_timestamp_t t1 = md_time_current();
        MD_LOG_DEBUG("Postprocess: compute bonds %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
    }

    if (flags & MD_UTIL_POSTPROCESS_STRUCTURE_BIT) {
        if (mol->bond.count) {
#ifdef PROFILE
            md_timestamp_t t0 = md_time_current();
#endif
            if (!mol->structure.alloc) {
                mol->structure.alloc = alloc;
            }
            md_index_data_clear(&mol->structure);
            md_util_compute_structures(&mol->structure, &mol->bond);
#ifdef PROFILE
            md_timestamp_t t1 = md_time_current();
#endif
            if (!mol->ring.alloc) {
                mol->ring.alloc = alloc;
            }
            md_index_data_clear(&mol->ring);
            md_util_compute_rings(&mol->ring, &mol->atom, &mol->bond);
#ifdef PROFILE
            md_timestamp_t t2 = md_time_current();
            MD_LOG_DEBUG("Postprocess: compute structures %.3f ms\n", md_time_as_milliseconds(t1-t0));
            MD_LOG_DEBUG("Postprocess: compute rings %.3f ms\n", md_time_as_milliseconds(t2-t1));
#endif
        }
    }
    
#if 0
    if (flags & MD_UTIL_POSTPROCESS_ORDER_BIT) {
#ifdef PROFILE
        md_timestamp_t t0 = md_time_current();
#endif
        compute_covalent_bond_order(&mol->bond, &mol->atom, &mol->ring);
#ifdef PROFILE
        md_timestamp_t t1 = md_time_current();
        MD_LOG_DEBUG("Postprocess: compute bond order %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
    }
#endif

    if (flags & MD_UTIL_POSTPROCESS_ION_BIT) {
        if (mol->atom.element) {
#ifdef PROFILE
            md_timestamp_t t0 = md_time_current();
#endif
            md_util_identify_ions(&mol->atom, &mol->bond);
#ifdef PROFILE
            md_timestamp_t t1 = md_time_current();
            MD_LOG_DEBUG("Postprocess: identify ions %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_HBOND_BIT) {
        if (mol->atom.element && mol->atom.count > 0 && mol->bond.count > 0) {
            
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_CHAINS_BIT) {
        if (mol->chain.count == 0 && mol->residue.count > 0 && mol->bond.pairs) {
#ifdef PROFILE
            md_timestamp_t t0 = md_time_current();
#endif
            md_util_compute_chain_data(&mol->chain, &mol->atom, &mol->residue, &mol->bond, alloc);
#ifdef PROFILE
            md_timestamp_t t1 = md_time_current();
            MD_LOG_DEBUG("Postprocess: compute chains %.3f ms\n", md_time_as_milliseconds(t1-t0));
#endif
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_BACKBONE_BIT) {
        if (mol->chain.count && mol->atom.type) {
            // Compute backbone data
            // 
            // @NOTE: We should only attempt to compute backbone data for valid residues (e.g. amino acids / dna)
            // Backbones are not directly tied to chains and therefore we cannot use chains as a 1:1 mapping for the backbones.
            // We look within the chains and see if we can find consecutive ranges which form backbones.

            static const size_t MIN_BACKBONE_LENGTH = 3;
            static const size_t MAX_BACKBONE_LENGTH = 1<<15;

            size_t temp_pos = md_temp_get_pos();

            {
                md_protein_backbone_atoms_t* backbone_atoms = md_temp_push(MAX_BACKBONE_LENGTH * sizeof(md_protein_backbone_atoms_t));
                size_t backbone_length = 0;
                md_residue_idx_t res_base = -1;
        
                for (size_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
                    md_range_t range = md_chain_residue_range(mol->chain, chain_idx);
                    for (md_residue_idx_t res_idx = range.beg; res_idx < range.end; ++res_idx) {

                        if (res_idx == 130) {
                            while(0) {};
                        }

                        md_protein_backbone_atoms_t atoms = {-1, -1, -1, -1, -1};
                        md_range_t atom_range = md_residue_atom_range(mol->residue, res_idx);
                        if (md_util_protein_backbone_atoms_extract(&atoms, mol->atom.type + atom_range.beg, atom_range.end - atom_range.beg, atom_range.beg)) {
                            backbone_atoms[backbone_length++] = atoms;
                            if (res_base == -1) {
                                res_base = res_idx;
                            }
                        } else {
                            if (backbone_length >= MIN_BACKBONE_LENGTH && res_base != -1) {
                                // Commit the backbone
                                commit_protein_backbone(backbone_atoms, backbone_length, res_base, (md_chain_idx_t)chain_idx, mol, alloc);
                            }
                            backbone_length = 0;
                            res_base = -1;
                        }
                    }
                    // Possibly commit remainder of the chain
                    if (backbone_length >= MIN_BACKBONE_LENGTH && res_base != -1) {
                        commit_protein_backbone(backbone_atoms, backbone_length, res_base, (md_chain_idx_t)chain_idx, mol, alloc);
                    }
                    backbone_length = 0;
                    res_base = -1;
                }

                mol->protein_backbone.range.count = md_array_size(mol->protein_backbone.range.offset);
                if (mol->protein_backbone.range.count) {
                    // Add end offset
                    md_array_push(mol->protein_backbone.range.offset, (uint32_t)md_array_size(mol->protein_backbone.atoms), alloc);
                }
                mol->protein_backbone.count = md_array_size(mol->protein_backbone.atoms);

                md_array_resize(mol->protein_backbone.angle, mol->protein_backbone.count, alloc);
                md_array_resize(mol->protein_backbone.secondary_structure, mol->protein_backbone.count, alloc);
                md_array_resize(mol->protein_backbone.ramachandran_type, mol->protein_backbone.count, alloc);

                if (mol->protein_backbone.count > 0) {
                    md_util_backbone_angles_compute(mol->protein_backbone.angle, mol->protein_backbone.count, mol);
                    md_util_backbone_secondary_structure_compute(mol->protein_backbone.secondary_structure, mol->protein_backbone.count, mol);
                    md_util_backbone_ramachandran_classify(mol->protein_backbone.ramachandran_type, mol->protein_backbone.count, mol);
                }
            }
            md_temp_set_pos_back(temp_pos);

#if 1
            {
                md_nucleic_backbone_atoms_t* backbone_atoms = md_temp_push(MAX_BACKBONE_LENGTH * sizeof(md_nucleic_backbone_atoms_t));
                size_t backbone_length = 0;
                md_residue_idx_t res_base = -1;

                for (size_t chain_idx = 0; chain_idx < mol->chain.count; ++chain_idx) {
                    md_range_t range = md_chain_residue_range(mol->chain, chain_idx);
                    for (md_residue_idx_t res_idx = range.beg; res_idx < range.end; ++res_idx) {
                        md_nucleic_backbone_atoms_t atoms;
                        md_range_t atom_range = md_residue_atom_range(mol->residue, res_idx);
                        const md_label_t* atom_labels = mol->atom.type + atom_range.beg;
                        if (md_util_nucleic_backbone_atoms_extract(&atoms, atom_labels, atom_range.end - atom_range.beg, atom_range.beg)) {
                            backbone_atoms[backbone_length++] = atoms;
                            if (res_base == -1) {
                                res_base = res_idx;
                            }
                        } else {
                            if (backbone_length >= MIN_BACKBONE_LENGTH && res_base != -1) {
                                // Commit the backbone
                                commit_nucleic_backbone(backbone_atoms, backbone_length, res_base, (md_chain_idx_t)chain_idx, mol, alloc);
                            }
                            backbone_length = 0;
                            res_base = -1;
                        }
                    }
                    // Possibly commit remainder of the chain
                    if (backbone_length >= MIN_BACKBONE_LENGTH && res_base != -1) {
                        int32_t res_idx = mol->chain.res_range[chain_idx].end - (int32_t)backbone_length;
                        commit_nucleic_backbone(backbone_atoms, backbone_length, res_base, (md_chain_idx_t)chain_idx, mol, alloc);
                    }
                    backbone_length = 0;
                    res_base = -1;
                }

                mol->nucleic_backbone.range.count = md_array_size(mol->protein_backbone.range.offset);
                if (mol->nucleic_backbone.range.count) {
                    // Add end offset
                    md_array_push(mol->nucleic_backbone.range.offset, (uint32_t)md_array_size(mol->protein_backbone.atoms), alloc);
                }
                mol->nucleic_backbone.count = md_array_size(mol->nucleic_backbone.atoms);
            }
            md_temp_set_pos_back(temp_pos);
#endif
        }
    }

    if ((mol->unit_cell.flags != MD_UNIT_CELL_FLAG_NONE) &&
        mol->protein_backbone.count)
    {
        for (size_t i = 0; i < mol->protein_backbone.range.count; ++i) {
            md_chain_idx_t chain_idx = mol->protein_backbone.range.chain_idx[i];
            size_t beg = md_chain_atom_range(mol->chain, chain_idx).beg;
            size_t len = md_chain_atom_count(mol->chain, chain_idx);
            md_util_unwrap(mol->atom.x + beg, mol->atom.y + beg, mol->atom.z + beg, 0, len, &mol->unit_cell);
        }
    }

    //MD_LOG_DEBUG("Finished Postprocessing Molecule");

    return true;
}

/*
This is blatantly stolen and modified from 'Arseny Kapoulkines' goldnugget 'Mesh optimizer'
https://github.com/zeux/meshoptimizer/

MIT License
Copyright (c) 2016-2022 Arseny Kapoulkine
*/

// "Insert" two 0 bits after each of the 10 low bits of x
static inline uint32_t part1By2(uint32_t x) {
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x << 8))  & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x << 4))  & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x << 2))  & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
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

static void computeHistogram10(uint32_t hist[1024][3], const uint32_t* data, size_t count) {
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

static void radixPass10(uint32_t* destination, const uint32_t* source, const uint32_t* keys, size_t count, uint32_t hist[1024][3], int pass) {
    int bitoff = pass * 10;
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = (keys[source[i]] >> bitoff) & 1023;
        destination[hist[id][pass]++] = source[i];
    }
}

void md_util_sort_spatial(uint32_t* source, const float* x, const float* y, const float* z, size_t count) {
    if (!source || !z || !y || !z || count <= 0) return;

    md_allocator_i* alloc = md_get_heap_allocator();

    uint32_t* keys = md_alloc(alloc, sizeof(uint32_t) * count);
    computeOrder(keys, x, y, z, count, 0);

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram10(hist, keys, count);

    uint32_t* scratch = md_alloc(alloc, sizeof(uint32_t) * count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass 10-bit radix sort computes the resulting order into scratch
    radixPass10(source, scratch, keys, count, hist, 0);
    radixPass10(scratch, source, keys, count, hist, 1);
    radixPass10(source, scratch, keys, count, hist, 2);

    md_free(alloc, scratch, sizeof(uint32_t) * count);
    md_free(alloc, keys,    sizeof(uint32_t) * count);
}

void md_util_sort_spatial_vec3(uint32_t* source, const vec3_t* xyz, size_t count) {
    if (!source || !xyz || count <= 0) return;

    md_allocator_i* alloc = md_get_heap_allocator();

    uint32_t* keys = md_alloc(alloc, sizeof(uint32_t) * count);
    const float* base = (const float*)xyz;
    const float* x = base + 0;
    const float* y = base + 1;
    const float* z = base + 2;
    computeOrder(keys, x, y, z, count, sizeof(vec3_t));

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram10(hist, keys, count);

    uint32_t* scratch = md_alloc(alloc, sizeof(uint32_t) * count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass radix sort computes the resulting order into scratch
    radixPass10(source, scratch, keys, count, hist, 0);
    radixPass10(scratch, source, keys, count, hist, 1);
    radixPass10(source, scratch, keys, count, hist, 2);

    md_free(alloc, scratch, sizeof(uint32_t) * count);
    md_free(alloc, keys,    sizeof(uint32_t) * count);
}

void md_util_sort_radix_inplace_uint32(uint32_t* data, size_t count) {
    if (!data || count <= 0) return;

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));
    sort_radix_inplace_uint32(data, count, temp_arena);
    md_vm_arena_destroy(temp_arena);
}

// Non inplace version of radix sort which should fill in the source indices which represents the sorted order
//void md_util_sort_radix_uint32(uint32_t* source_indices, const uint32_t* data, size_t count) {
//}

#define DEBUG_PRINT 0

static void backtrack(state_t* state) {
    uint32_t depth = (uint32_t)md_array_size(state->n_path);

    if (md_array_size(state->n_path)) {
        int n_idx = md_array_back(state->n_path);
        md_array_pop(state->n_path);
        state->map[n_idx] = -1;
        bitfield_clear_bit(state->n_path_bits, n_idx);
    }
    if (md_array_size(state->h_path)) {
        int h_idx = md_array_back(state->h_path);
        md_array_pop(state->h_path);
        bitfield_clear_bit(state->h_path_bits, h_idx);
    }

    for (size_t i = 0; i < md_array_size(state->n_depths); ++i) {
        //state->n_depths[i] = (state->n_depths[i] == depth) ? 0 : state->n_depths[i];
        if (state->n_depths[i] == depth) {
            state->n_depths[i] = 0;
        }
    }
    for (size_t i = 0; i < md_array_size(state->h_depths); ++i) {
        //state->h_depths[i] = (state->h_depths[i] == depth) ? 0 : state->h_depths[i];
        if (state->h_depths[i] == depth) {
            state->h_depths[i] = 0;
        }
    }
}

static inline bool is_in_terminal_set(const uint32_t depths[], const uint64_t path_bits[], int64_t i) {
    return !bitfield_test_bit(path_bits, i) && depths[i];
    //if (bitfield_test_bit(path_bits, i)) return false;
    //if (!depths[i]) return false;
    //return true;
}

static inline bool check_bonds(state_t* state, int n_idx) {
    graph_edge_iter_t it = graph_edge_iter(state->n_graph, n_idx);
    int h_idx = state->map[n_idx];
    while (graph_edge_iter_has_next(it)) {
        int other_n_idx = graph_edge_iter_vidx(it);
        int other_h_idx = state->map[other_n_idx];

        if (other_h_idx != -1) {
            int type = graph_edge_iter_type(it);
            if (type == 0) {
                if (!graph_vertex_is_connected_to(state->h_graph, h_idx, other_h_idx)) {
                    return false;
                }
            } else {
                if (!graph_vertex_has_connection(state->h_graph, h_idx, other_h_idx, type)) {
                    return false;
                }
            }
        }

        graph_edge_iter_next(&it);
    }
    return true;
}

// Check for solution
static inline bool check_map(state_t* state) {
    if (md_array_size(state->n_path) == state->n_graph->vertex_count) {
        return state->callback(state->map, md_array_size(state->map), state->user_data);
    }
    return false;
}

static inline bool test_edge_constraints(state_t* state, int n_idx, int h_idx) {
    if (state->flags & MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT && graph_vertex_type(state->n_graph, n_idx) > 0) {
        size_t cn = graph_vertex_edge_count(state->n_graph, n_idx);
        size_t ch = graph_vertex_edge_count(state->h_graph, h_idx);
        if (cn != ch) return false;
    }

    graph_edge_iter_t it = graph_edge_iter(state->n_graph, n_idx);
    while (graph_edge_iter_has_next(it)) {
        int other_n_idx = graph_edge_iter_vidx(it);
        int other_h_idx = state->map[other_n_idx];

        if (other_h_idx != -1) {
            int edge_type = graph_edge_iter_type(it);
            if (!(state->flags & MD_UTIL_MATCH_FLAGS_STRICT_EDGE_TYPE)) {
                edge_type = 0;
            }
            if (edge_type == 0) {
                if (!graph_vertex_is_connected_to(state->h_graph, h_idx, other_h_idx)) {
                    return false;
                }
            } else {
                if (!graph_vertex_has_connection(state->h_graph, h_idx, other_h_idx, edge_type)) {
                    return false;
                }
            }
        }

        graph_edge_iter_next(&it);
    }
    return true;
}

static bool match_candidate(state_t* state, int n_idx, int h_idx) {
    int vn = graph_vertex_type(state->n_graph, n_idx);
    int vh = graph_vertex_type(state->h_graph, h_idx);
    if (vn && vh && vn != vh) {
        return false;
    }

    state->map[n_idx] = h_idx;

    if (!test_edge_constraints(state, n_idx, h_idx)) {
        state->map[n_idx] = -1;
        return false;
    }

    md_array_push_no_grow(state->n_path, n_idx);
    md_array_push_no_grow(state->h_path, h_idx);

    bitfield_set_bit(state->n_path_bits, n_idx);
    bitfield_set_bit(state->h_path_bits, h_idx);

    // Update n_depths
    {
        const uint32_t d = (uint32_t)md_array_size(state->n_path);
        state->n_depths[n_idx] = state->n_depths[n_idx] == 0 ? d : state->n_depths[n_idx];

        graph_edge_iter_t it = graph_edge_iter(state->n_graph, n_idx);
        while (graph_edge_iter_has_next(it)) {
            int vidx = graph_edge_iter_vidx(it);
            state->n_depths[vidx] = state->n_depths[vidx] == 0 ? d : state->n_depths[vidx];
            graph_edge_iter_next(&it);
        }
    }

    // Update h_depths
    {
        const uint32_t d = (uint32_t)md_array_size(state->h_path);
        graph_edge_iter_t it = graph_edge_iter(state->h_graph, h_idx);
        while (graph_edge_iter_has_next(it)) {
            int vidx = graph_edge_iter_vidx(it);
            if (!state->h_depths[vidx]) {
                state->h_depths[vidx] = d;
            }
            graph_edge_iter_next(&it);
        }
    }

    //if (!check_bonds(state, n_idx)) {
    //    backtrack(state);
    //	return false;
    //}

    // VF2 Feasability check
    const size_t N1 = state->n_graph->vertex_count;
    const size_t N2 = state->h_graph->vertex_count;

    size_t T1 = 0;
    for (size_t i = 0; i < N1; ++i) {
        if (is_in_terminal_set(state->n_depths, state->n_path_bits, i)) {
            ++T1;
        }
    }

    size_t T2 = 0;
    for (size_t i = 0; i < N2; ++i) {
        if (is_in_terminal_set(state->h_depths, state->h_path_bits, i)) {
            ++T2;
        }
    }

    const size_t M1 = md_array_size(state->n_path);
    const size_t M2 = md_array_size(state->h_path);

    if (T1 > T2 || (N1 - M1 - T1) > (N2 - M2 - T2)) {
        backtrack(state);
        return false;
    }

    state->abort = check_map(state);
    return true;
}

static inline int terminal_size(const uint32_t depths[], size_t len) {
    int count = 0;
    for (size_t i = 0; i < len; ++i) {
        count += depths[i] ? 1 : 0;
    }
    return count;
}

static bool next_candidate(const state_t* state, int* n_idx, int* h_idx) {
    int last_n_idx = *n_idx == -1 ? 0 : *n_idx;
    int last_h_idx = *h_idx == -1 ? 0 : *h_idx + 1;

    int n_size = (int)state->n_graph->vertex_count;
    int h_size = (int)state->h_graph->vertex_count;

    int n_term_size = terminal_size(state->n_depths, md_array_size(state->n_depths));
    int h_term_size = terminal_size(state->h_depths, md_array_size(state->h_depths));

    int map_size = (int)md_array_size(state->n_path);

    if (n_term_size > map_size && h_term_size > map_size) {
        ASSERT(state->n_depths);
        while (last_n_idx < n_size && (bitfield_test_bit(state->n_path_bits, last_n_idx) || !state->n_depths[last_n_idx])) {
            ++last_n_idx;
            last_h_idx = 0;
        }
    } else {
        while (last_n_idx < n_size && bitfield_test_bit(state->n_path_bits, last_n_idx)) {
            ++last_n_idx;
            last_h_idx = 0;
        }
    }

    if (n_term_size > map_size && h_term_size > map_size) {
        ASSERT(state->h_depths);
        while (last_h_idx < h_size && (bitfield_test_bit(state->h_path_bits, last_h_idx) || !state->h_depths[last_h_idx])) {
            last_h_idx++;
        }
    } else {
        while (last_h_idx < h_size && bitfield_test_bit(state->h_path_bits, last_h_idx)) {
            last_h_idx++;
        }
    }

    if (last_n_idx < n_size && last_h_idx < h_size) {
        *n_idx = last_n_idx;
        *h_idx = last_h_idx;
        return true;
    }

    return false;
}

static void map_next(state_t* state) {
    // @TODO: Timeout check?

    if (state->abort) {
        return;
    }

    int n_idx = -1;
    int h_idx = -1;
    while (!state->abort) {
        if (!next_candidate(state, &n_idx, &h_idx)) {
            return;
        }

#if DEBUG_PRINT
        const int depth = (int)md_array_size(state->n_path);
        printf("%*s [%d, %d]\n", depth, "", n_idx, h_idx);
#endif

        if (match_candidate(state, n_idx, h_idx)) {
            map_next(state);
            backtrack(state);
        }
    }
}

typedef struct {
    uint32_t value;
    uint32_t index;
} result_value_t;

typedef struct MD_HASHMAP_T(result_value_t) result_map_t;

typedef struct store_data_t {
    size_t* count;
    md_index_data_t* result;
    md_allocator_i*  alloc;
    // Only for unique
    result_map_t map;
} store_data_t;

static bool store_unique_callback(const int map[], size_t length, void* user) {
    store_data_t* data = (store_data_t*)user;

    uint64_t key = 0;
    uint32_t fit = 0;
    for (size_t i = 0; i < length; ++i) {
        key += map[i];
        // The fit value is just a arbitrary metric that tells us how good of a 'fit' this mapping is
        // It is the sum of the difference in index, this means a value of zero would correspond to an identical fit in terms of indices
        fit += (uint32_t)map[i] * (uint32_t)i;
    }
    result_value_t* entry = md_hashmap_get(&data->map, key);
    if (entry) {
        if (entry->value < fit) {
            // Update result if the new value is lower (more strictly sorted)
            entry->value = fit;
            int* ptr = md_index_range_beg(*data->result, entry->index);
            MEMCPY(ptr, map, sizeof(int) * length);
        }
    } else {
        // Store new entry if unique solution
        uint32_t result_idx = (uint32_t)md_index_data_push_arr(data->result, map, length);
        *data->count += 1;
        result_value_t val = {fit, result_idx};
        md_hashmap_add(&data->map, key, val);
    }

#if DEBUG_PRINT
    const int depth = (int)length;
    printf("%*sSolution found!\n", depth, "");
#endif

    return false;
}

static bool store_first_callback(const int map[], size_t length, void* user) {
    store_data_t* data = (store_data_t*)user;

    md_index_data_push_arr(data->result, map, length);
    *data->count += 1;
#if DEBUG_PRINT
    const int depth = (int)length;
    printf("%*sSolution found!\n", depth, "");
#endif

    return true;
}

static bool store_all_callback(const int map[], size_t length, void* user) {
    store_data_t* data = (store_data_t*)user;

    md_index_data_push_arr(data->result, map, length);
    *data->count += 1;
#if DEBUG_PRINT
    const int depth = (int)length;
    printf("%*sSolution found!\n", depth, "");
#endif

    return false;
}

static bool count_first_callback(const int map[], size_t length, void* user) {
    (void)map;
    (void)length;
    size_t* count = (size_t*)user;
    *count += 1;
    return true;
}

static bool count_all_callback(const int map[], size_t length, void* user) {
    (void)map;
    (void)length;
    size_t* count = (size_t*)user;
    *count += 1;
    return false;
}

// Attempt to find subgraphs in a larger graph (haystack) which match a reference graph (needle)
// Returns an array of graphs which match the topologically match the reference
// start_type is a hint of the most unusual type in the graphs and serve as good starting points
static void find_isomorphisms_callback(const graph_t* needle, const graph_t* haystack, int start_type, state_t* state) {
    ASSERT(needle);
    ASSERT(haystack);
    ASSERT(state);

    // Impossible case
    if (needle->vertex_count > haystack->vertex_count) {
        return;
    }

    // Check for equivalence, if we have a 1:1 mapping
    if (needle->vertex_count == haystack->vertex_count) {
        if (graph_equivalent(needle, haystack)) {
            size_t length = haystack->vertex_count;
            for (int i = 0; i < (int)length; ++i) {
                state->map[i] = i;
            }
            if (state->callback(state->map, length, state->user_data)) {
                return;
            }
        }
    }

    // The problematic case (subgraph isomorphism)

    // Create list of starting candidate pairs
    for (int h_idx = 0; h_idx < (int)haystack->vertex_count; ++h_idx) {
        if (graph_vertex_type(haystack, h_idx) != start_type) continue;
        for (int n_idx = 0; n_idx < (int)needle->vertex_count; ++n_idx) {
            if (graph_vertex_type(needle, n_idx) != start_type) continue;
            // Set initial state
            if (match_candidate(state, n_idx, h_idx)) {
                map_next(state);
            }

            if (state->abort) {
                return;
            }

            // Reset state
            state_reset(state);
        }
    }
}

// Attempt to find subgraphs in a larger graph (haystack) which match a reference graph (needle)
// Returns an array of graphs which match the topologically match the reference
// start_type is a hint of the most unusual type in the graphs and serve as good starting points
static size_t find_isomorphisms(md_index_data_t* mappings, const graph_t* needle, const graph_t* haystack, md_util_match_mode_t mode, int start_type, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    size_t count = 0;
    
    // Impossible case
    if (needle->vertex_count > haystack->vertex_count) {
        return 0;
    }

    // Check for equivalence (graph isomorphism)
    if (mode == MD_UTIL_MATCH_MODE_FIRST &&
        needle->vertex_count == haystack->vertex_count) {
        if (graph_equivalent(needle, haystack)) {
            // This should be a 1:1 mapping
            if (mappings) {
                size_t length = haystack->vertex_count;
                int* indices = md_vm_arena_push_array(temp_alloc, int, length);
                for (int i = 0; i < (int)length; ++i) {
                    indices[i] = i;
                }
                md_index_data_push_arr(mappings, indices, length);
            }
            return 1;
        }
    }

    // The problematic case (subgraph isomorphism)
    typedef struct {
        int n_idx;
        int h_idx;
    } pair_t;

    // Create list of starting candidate pairs
    md_array(pair_t) start_candidates = 0;
    for (int i = 0; i < (int)needle->vertex_count; ++i) {
        if (graph_vertex_type(needle, i) != start_type) continue;
        for (int j = 0; j < (int)haystack->vertex_count; ++j) {
            if (graph_vertex_type(haystack, j) != start_type) continue;
            md_array_push(start_candidates, ((pair_t){i,j}), temp_alloc);
        }
    }

    if (start_candidates == 0) goto done;

    state_t state = {0};
    state_init(&state, needle, haystack, temp_alloc);
    store_data_t store_data = {
        .count = &count,
        .map = { .allocator = temp_alloc },
    };

    if (!mappings) {
        state.user_data = &count;
        switch (mode) {
        case MD_UTIL_MATCH_MODE_UNIQUE:
            MD_LOG_ERROR("Cannot count unique occurrences without supplying mappings");
            goto done;
            break;
        case MD_UTIL_MATCH_MODE_FIRST:
            state.callback = count_first_callback;
            break;
        case MD_UTIL_MATCH_MODE_ALL:
            state.callback = count_all_callback;
            break;
        default:
            ASSERT(false);
            break;
        }
    } else {
        store_data.result = mappings;
        store_data.alloc = alloc;
        state.user_data = &store_data;

        switch (mode) {
        case MD_UTIL_MATCH_MODE_UNIQUE:
            state.callback = store_unique_callback;
            break;
        case MD_UTIL_MATCH_MODE_FIRST:
            state.callback = store_first_callback;
            break;
        case MD_UTIL_MATCH_MODE_ALL:
            state.callback = store_all_callback;
            break;
        default:
            ASSERT(false);
            break;
        }
    }

    const size_t num_candidates = md_array_size(start_candidates);
    for (size_t i = 0; i < num_candidates; ++i) {
        int n_idx = start_candidates[i].n_idx;
        int h_idx = start_candidates[i].h_idx;

#if DEBUG_PRINT
        printf("STARTING ATTEMPT TO MATCH %d -> %d\n", n_idx, h_idx);
#endif
        // Reset state
        state_reset(&state);

        // Set initial state
        if (match_candidate(&state, n_idx, h_idx)) {
            map_next(&state);
        }

        if (mode == MD_UTIL_MATCH_MODE_FIRST && count > 0) {
            goto done;
        }
    }
    
done:    
    return count;
}

enum {
    VERTEX_TYPE_MAPPING_UNDEFINED = 0,
    VERTEX_TYPE_MAPPING_TYPE = 1,
    VERTEX_TYPE_MAPPING_ELEMENT = 2,
    VERTEX_TYPE_MAPPING_MASS = 3
};

typedef int vertex_type_mapping_t;

// Create a new reference structure which is pruned of certain atoms (Hydrogen) and loosely connected subcomponents
// There are simply too many permutations to cover and the result will explode.
static md_array(int) filter_structure_connectivity(const int* indices, size_t count, const md_index_data_t* connectivity, int min_val, md_allocator_i* alloc) {
    md_array(int) filt = 0;
    for (size_t i = 0; i < count; ++i) {
        int idx = indices[i];
        size_t order = md_index_range_size(*connectivity, idx);
        if (order >= min_val) {
            md_array_push(filt, idx, alloc);
        }
    }
    return filt;
}

static md_array(int) filter_structure_type(const int* indices, int64_t count, const uint8_t* type, uint8_t type_to_remove, md_allocator_i* alloc) {
    md_array(int) filt = 0;
    for (int64_t i = 0; i < count; ++i) {
        int idx = indices[i];
        if (type[idx] != type_to_remove) {
            md_array_push(filt, idx, alloc);
        }
    }
    return filt;
}

static inline bool filter_atom(const md_molecule_t* mol, int atom_i, md_util_match_flags_t filter) {
    if (filter & MD_UTIL_MATCH_FLAGS_NO_H) {
        return mol->atom.element[atom_i] != H;
    } else if (filter & MD_UTIL_MATCH_FLAGS_NO_CH) {
        if (mol->atom.element[atom_i] == H) {
            uint32_t   conn_len = (uint32_t)md_bond_conn_count(mol->bond, atom_i);
            uint32_t   conn_i   = mol->bond.conn.offset[atom_i];
            for (uint32_t i = 0; i < conn_len; ++i) {
                int atom_j = md_bond_conn_atom_idx(mol->bond, conn_i, i);
                if (mol->atom.element[atom_j] == C) {
                    return false;
                }
            }
        }
    }
    return true;
}

typedef struct idx_range_t {
    md_index_data_t* idx_data;
    size_t idx_offset;
} idx_range_t;

static inline idx_range_t idx_range_create(md_index_data_t* idx_data) {
    idx_range_t range = {idx_data, md_array_size(idx_data->indices)};
    return range;
}

static inline void idx_range_commit(idx_range_t range) {
    size_t cur_offset = md_array_size(range.idx_data->indices);

    if (cur_offset != range.idx_offset) {
        if (md_array_size(range.idx_data->offsets) == 0) {
            md_array_push(range.idx_data->offsets, 0, range.idx_data->alloc);
        }

        // First insertion in this range, create a new offset pair corresponding to the new range
        md_array_push(range.idx_data->offsets, (uint32_t)cur_offset, range.idx_data->alloc);
    }
}

static void idx_range_push(idx_range_t range, int32_t idx) {
    md_array_push(range.idx_data->indices, idx, range.idx_data->alloc);
}

static size_t extract_structures(md_index_data_t* out_structures, const md_molecule_t* mol, md_util_match_level_t level, md_util_match_flags_t filter, size_t min_size) {
    ASSERT(out_structures);
    size_t pre_offset = md_index_data_num_ranges(*out_structures);

    switch (level) {
    case MD_UTIL_MATCH_LEVEL_STRUCTURE:
        if (filter) {
            for (size_t s_idx = 0; s_idx < md_index_data_num_ranges(mol->structure); ++s_idx) {
                if (md_index_range_size(mol->structure, s_idx) < min_size) {
                    continue;
                }
                idx_range_t range = idx_range_create(out_structures);
                for (int* it = md_index_range_beg(mol->structure, s_idx); it < md_index_range_end(mol->structure, s_idx); ++it) {
                    int i = *it;
                    if (filter_atom(mol, i, filter)) {
                        idx_range_push(range, i);
                    }
                }
                idx_range_commit(range);
            }
        } else {
            for (size_t s_idx = 0; s_idx < md_index_data_num_ranges(mol->structure); ++s_idx) {
                if (md_index_range_size(mol->structure, s_idx) < min_size) {
                    continue;
                }
                md_index_data_push_arr(out_structures, md_index_range_beg(mol->structure, s_idx), md_index_range_size(mol->structure, s_idx));
            }
        }
        break;
    case MD_UTIL_MATCH_LEVEL_RESIDUE:
        for (size_t r_idx = 0; r_idx < mol->residue.count; ++r_idx) {
            if (md_residue_atom_count(mol->residue, r_idx) < min_size) {
                continue;
            }
            idx_range_t range = idx_range_create(out_structures);
            md_range_t res_range = md_residue_atom_range(mol->residue, r_idx);
            for (int i = res_range.beg; i < res_range.end; ++i) {
                if (filter_atom(mol, i, filter)) {
                    idx_range_push(range, i);
                }
            }
            idx_range_commit(range);
        }
        break;
    case MD_UTIL_MATCH_LEVEL_CHAIN:
        for (size_t c_idx = 0; c_idx < mol->chain.count; ++c_idx) {
            if (md_chain_atom_count(mol->chain, c_idx) < min_size) {
                continue;
            }
            idx_range_t range = idx_range_create(out_structures);
            md_range_t chain_range = md_chain_atom_range(mol->chain, c_idx);
            for (int i = chain_range.beg; i < chain_range.end; ++i) {
                if (filter_atom(mol, i, filter)) {
                    idx_range_push(range, i);
                }
            }
            idx_range_commit(range);
        }
        break;
    default:
        ASSERT(false);
    }

    size_t post_offset = md_index_data_num_ranges(*out_structures);

    return post_offset - pre_offset;
}

#define MAX_TYPES 256
md_index_data_t match_structure(const int* ref_idx, size_t ref_len, md_util_match_mode_t mode, md_util_match_level_t level, vertex_type_mapping_t mapping, const md_molecule_t* mol, md_allocator_i* alloc) {
    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    graph_t ref_graph = {0};
    md_index_data_t result = {.alloc = alloc};

    md_hashmap32_t map_table = { .allocator = temp_arena };

    md_array(uint8_t) atom_type = md_vm_arena_push_array(temp_arena, uint8_t, mol->atom.count);
    md_array(int) structure_idx = md_vm_arena_push_array(temp_arena, int,     mol->atom.count);

    if (mapping == VERTEX_TYPE_MAPPING_ELEMENT) {
        for (size_t i = 0; i < mol->atom.count; ++i) {
            atom_type[i] = mol->atom.element[i];
        }
    } else {
        uint32_t type_count = 0;
        for (size_t i = 0; i < mol->atom.count; ++i) {
            uint64_t key = 0;
            switch (mapping) {
            case VERTEX_TYPE_MAPPING_TYPE:
                key = md_hash64(mol->atom.type[i].buf, mol->atom.type[i].len, 0);
                break;
            case VERTEX_TYPE_MAPPING_MASS:
                key = md_hash64(&mol->atom.mass[i], sizeof(mol->atom.mass[i]), 0);
                break;
            default: ASSERT(false); break;
            }
            uint8_t type = 0;
            uint32_t* val = md_hashmap_get(&map_table, key);
            if (!val) {
                type = (uint8_t)type_count++;
                md_hashmap_add(&map_table, key, 1);
            } else {
                type = (uint8_t)*val;
                *val += 1;
            }
            atom_type[i] = type;
        }
    }

    for (size_t i = 0; i < md_index_data_num_ranges(mol->structure); ++i) {
        int* beg = md_index_range_beg(mol->structure, i);
        int* end = md_index_range_end(mol->structure, i);
        for (int* it = beg; it != end; ++it) {
            structure_idx[*it] = (int)i;
        }
    }

    // Ensure that the reference indices all belong to the same structure and are not disjoint
    int ref_structure_idx = structure_idx[ref_idx[0]];
    for (size_t i = 1; i < ref_len; ++i) {
        if (structure_idx[ref_idx[i]] != ref_structure_idx) {
            MD_LOG_ERROR("Reference indices are not part of the same structure, they are disconnected");
            goto done;
        }
    }

    bool ref_hydro_present = false;
    for (size_t i = 0; i < ref_len; ++i) {
        int idx = ref_idx[i];
        if (mol->atom.element[idx] == H) {
            ref_hydro_present = true;
            break;
        }
    }

    const md_util_match_flags_t flags = ref_hydro_present ? 0 : MD_UTIL_MATCH_FLAGS_NO_H;

    ref_graph = make_graph(&mol->bond, atom_type, ref_idx, ref_len, temp_arena);

    int ref_type_count[MAX_TYPES] = {0};
    for (int i = 0; i < ref_len; ++i) {
        uint8_t type = (uint8_t)graph_vertex_type(&ref_graph, i);
        ref_type_count[type]++;
    }

    md_index_data_t structures = { .alloc = temp_arena };
    const size_t num_structures = extract_structures(&structures, mol, level, flags, ref_len);

    int starting_type = -1;
    int min_freq = INT_MAX;
    for (int i = 0; i < MAX_TYPES; ++i) {
        int freq = ref_type_count[i];
        if (freq > 0 && freq < min_freq) {
            min_freq = freq;
            starting_type = i;
        }
    }

    size_t result_count = 0;
    store_data_t data = {
        .alloc = alloc,
        .result = &result,
        .count = &result_count,
        .map = {.allocator = temp_arena },
    };


    for (size_t i = 0; i < num_structures; ++i) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_arena);

        const int*   s_idx = md_index_range_beg(structures, i);
        const size_t s_len = md_index_range_size(structures, i);

        md_util_match_mode_t s_mode = mode;
        if (s_len == ref_len && s_mode == MD_UTIL_MATCH_MODE_UNIQUE) {
            s_mode = MD_UTIL_MATCH_MODE_FIRST;
        }

        solution_callback cb = 0;
        switch (s_mode) {
            case MD_UTIL_MATCH_MODE_UNIQUE:
                cb = store_unique_callback;
                break;
            case MD_UTIL_MATCH_MODE_FIRST: 
                cb = store_first_callback;
                break;
            case MD_UTIL_MATCH_MODE_ALL:
                cb = store_all_callback;
                break;
            default:
                ASSERT(false);
        }
        
        int s_type_count[MAX_TYPES] = {0};
        for (size_t j = 0; j < s_len; ++j) {
            int idx = s_idx[j];
            s_type_count[atom_type[idx]]++;
        }

        // Sanity check
        // The structure needs to have atleast the same amount of types as the reference
        for (int j = 0; j < MAX_TYPES; ++j) {
            if (ref_type_count[j] > s_type_count[j]) {
                goto next;
            }
        }

        graph_t graph = make_graph(&mol->bond, atom_type, s_idx, s_len, temp_arena);
        size_t pre_count = md_index_data_num_ranges(result);

        state_t state = {0};
        state_init(&state, &ref_graph, &graph, temp_arena);
        state.callback = cb;
        state.user_data = &data;

        MEMSET(&data.map, 0, sizeof(result_map_t));
        data.map.allocator = temp_arena;

        find_isomorphisms_callback(&ref_graph, &graph, starting_type, &state);
        size_t post_count = md_index_data_num_ranges(result);

        // Remap indices to global indices in result
        for (size_t j = pre_count; j < post_count; ++j) {
            int* beg = md_index_range_beg(result, j);
            int* end = md_index_range_end(result, j);
            for (int* it = beg; it != end; ++it) {
                *it = s_idx[*it];
            }
        }
    next:
        md_vm_arena_temp_end(temp);
    }

done:
    md_vm_arena_destroy(temp_arena);
    return result;
}

md_index_data_t md_util_match_by_type(const int ref_indices[], size_t ref_count, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc) {
    return match_structure(ref_indices, ref_count, mode, level, VERTEX_TYPE_MAPPING_TYPE, mol, alloc);
}

md_index_data_t md_util_match_by_element(const int ref_indices[], size_t ref_count, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc) {
    return match_structure(ref_indices, ref_count, mode, level, VERTEX_TYPE_MAPPING_ELEMENT, mol, alloc);
}

size_t md_util_match_smiles(md_index_data_t* idx_data, str_t smiles, md_util_match_mode_t mode, md_util_match_level_t level, md_util_match_flags_t flags, const md_molecule_t* mol, md_allocator_i* alloc) {
    ASSERT(mol);

    if (idx_data && !idx_data->alloc) {
        MD_LOG_ERROR("Incomplete idx_data structure supplied");
        return 0;
    }

    if (!mol->atom.element) {
        MD_LOG_ERROR("Molecule does not have any atom element field");
        return 0;
    }

    md_allocator_i* temp_alloc = md_vm_arena_create(GIGABYTES(4));
    graph_t ref_graph = smiles_to_graph(smiles, flags, temp_alloc, temp_alloc);

    // Histogram of types present in reference structure
    int ref_type_count[256] = {0};
    for (size_t i = 0; i < ref_graph.vertex_count; ++i) {
        uint8_t type = ref_graph.vertex_type[i];
        ref_type_count[type]++;
    }

    // If there is no hydrogen present in the reference pattern, add that to filter to limit the search scope.
    if (ref_type_count[1] == 0) {
        flags |= MD_UTIL_MATCH_FLAGS_NO_H;
    }

    md_index_data_t structures = { .alloc = temp_alloc };
    const size_t num_structures = extract_structures(&structures, mol, level, flags, ref_graph.vertex_count);

    size_t match_count = 0;

    if (num_structures == 0) {
        goto done;
    }

    // Find most uncommon type in dataset. This will be our starting point(s)
    int starting_type = -1;
    int min_freq = INT_MAX;
    // Exclude hydrogen
    for (int i = H + 1; i < Num_Elements; ++i) {
        int freq = ref_type_count[i];
        if (freq > 0 && freq < min_freq) {
            starting_type = i;
            min_freq = freq;
        }
    }

    for (size_t i = 0; i < num_structures; ++i) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_alloc);
        const size_t s_size = md_index_range_size(structures, i);
        const int*   s_idx  = md_index_range_beg(structures, i);

        int s_type_count[256] = {0};
        for (size_t j = 0; j < s_size; ++j) {
            int idx = s_idx[j];
            uint8_t type = mol->atom.element[idx];
            s_type_count[type]++;
        }

        // Sanity check
        for (size_t j = 0; j < Num_Elements; ++j) {
            if (ref_type_count[j] > s_type_count[j]) goto next;
        }

        graph_t s_graph = make_graph(&mol->bond, mol->atom.element, s_idx, s_size, temp_alloc);

        if (flags & MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT) {
            if (s_graph.vertex_count != ref_graph.vertex_count) {
                continue;
            }
        }
        
        size_t pre_count = md_index_data_num_ranges(*idx_data);
        size_t count = find_isomorphisms(idx_data, &ref_graph, &s_graph, mode, starting_type, alloc, temp_alloc);
        size_t post_count = md_index_data_num_ranges(*idx_data);

        // Remap indices to global indices in result
        if (idx_data && count > 0) {
            size_t num_ranges = md_index_data_num_ranges(*idx_data);
            ASSERT(post_count - pre_count == count);
            for (size_t j = pre_count; j < post_count; ++j) {
                int* beg = md_index_range_beg(*idx_data, j);
                int* end = md_index_range_end(*idx_data, j);
                for (int* it = beg; it != end; ++it) {
                    *it = s_idx[*it];
                }
            }
        }

        match_count += count;
    next:
        md_vm_arena_temp_end(temp);
    }

done:
    md_vm_arena_destroy(temp_alloc);

    return match_count;
}

#ifdef __cplusplus
}
#endif
