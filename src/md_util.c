#include <md_util.h>

#include <md_system.h>
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
#include <core/md_spatial_acc.h>
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
// Covalent radii (in Å ×100) — Cordero et al., Dalton Trans., 2008, 2832–2838
static const uint8_t element_covalent_radii_u8[] = {
      0,  31,  28, 128,  96,  84,  76,  71,  66,  57,  58, 166, 141, 121, 111, 107, 105, 102, 106, 203, 176, 170, 160, 153,
    139, 139, 132, 126, 124, 132, 122, 122, 120, 119, 120, 120, 116, 220, 195, 190, 175, 164, 154, 147, 146, 142, 139, 145,
    144, 142, 139, 139, 138, 139, 140, 244, 215, 207, 204, 203, 201, 199, 198, 198, 196, 194, 192, 192, 189, 190, 187, 187,
    175, 170, 162, 151, 144, 141, 136, 136, 132, 145, 146, 148, 140, 150, 150, 255, 221, 215, 206, 200, 196, 190, 187, 180,
    169, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160
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

#define RGBA(r,g,b,a) (uint32_t)( ((a & 255) << 24) | ((b & 255) << 16) | ((g & 255) << 8) | (r & 255) )
#define RGB(r,g,b) RGBA(r,g,b,255)

// Based on this with some values modified
// Some stem from Jmol, some from Rasmol
// http://jmol.sourceforge.net/jscolors/
static const uint32_t element_cpk_colors[] = {
    RGB(255,   0, 255), // Unknown (0)  rgb(255, 0, 255)
    RGB(255, 255, 255), // H (1)        rgb(255, 255, 255)
    RGB(217, 255, 255), // He (2)       rgb(217, 255, 255)
    RGB(204, 128, 255), // Li (3)       rgb(204, 128, 255)
    RGB(194, 255,   0), // Be (4)       rgb(194, 255, 0)
    RGB(255, 181, 181), // B (5)        rgb(255, 181, 181)
    RGB(131, 131, 131), // C (6)        rgb(131, 131, 131)
    RGB( 70,  70, 240), // N (7)        rgb(70, 70, 240)
    RGB(240,   0,   0), // O (8)        rgb(240, 0, 0)
    RGB(144, 224,  80), // F (9)        rgb(144, 224, 80)
    RGB(179, 227, 245), // Ne (10)      rgb(179, 227, 245)
    RGB(171,  92, 242), // Na (11)      rgb(171, 92, 242)
    RGB(138, 255,   0), // Mg (12)      rgb(138, 255, 0)
    RGB(128, 128, 144), // Al (13)      rgb(128, 128, 144)
    RGB(240, 200, 160), // Si (14)      rgb(240, 200, 160)
    RGB(255, 165,   0), // P (15)       rgb(255, 165, 0)
    RGB(240, 240,  10), // S (16)       rgb(240, 240, 10)
    RGB( 31, 240,  31), // Cl (17)      rgb(31, 240, 31)
    RGB(128, 209, 227), // Ar (18)      rgb(128, 209, 227)
    RGB(143,  64, 212), // K (19)       rgb(143, 64, 212)
    RGB(128, 128, 144), // Ca (20)      rgb(128, 128, 144)
    RGB(230, 230, 230), // Sc (21)      rgb(230, 230, 230)
    RGB(194, 194, 199), // Ti (22)      rgb(194, 194, 199)
    RGB(166, 166, 171), // V (23)       rgb(166, 166, 171)
    RGB(128, 128, 144), // Cr (24)      rgb(128, 128, 144)
    RGB(128, 128, 144), // Mn (25)      rgb(128, 128, 144)
    RGB(224, 102,  51), // Fe (26)      rgb(224, 102, 51)
    RGB(240, 144, 160), // Co (27)      rgb(240, 144, 160)
    RGB(80,  208,  80), // Ni (28)      rgb(80, 208, 80)
    RGB(200, 128,  51), // Cu (29)      rgb(165, 42, 42)
    RGB(125, 128, 176), // Zn (30)      rgb(125, 128, 176)
    RGB(194, 143, 143), // Ga (31)      rgb(194, 143, 143)
    RGB(102, 143, 143), // Ge (32)      rgb(102, 143, 143)
    RGB(189, 128, 227), // As (33)      rgb(189, 128, 227)
    RGB(255, 161,   0), // Se (34)      rgb(255, 161, 0)
    RGB(165,  42,  42), // Br (35)      rgb(165, 42, 42)
    RGB( 92, 184, 209), // Kr (36)      rgb(92, 184, 209)
    RGB(112,  46, 176), // Rb (37)      rgb(112, 46, 176)
    RGB(  0, 255,   0), // Sr (38)      rgb(0, 255, 0)
    RGB(148, 255, 255), // Y (39)       rgb(148, 255, 255)
    RGB(148, 224, 224), // Zr (40)      rgb(148, 224, 224)
    RGB(115, 194, 201), // Nb (41)      rgb(115, 194, 201)
    RGB( 84, 181, 181), // Mo (42)      rgb(84, 181, 181)
    RGB( 59, 158, 158), // Tc (43)      rgb(59, 158, 158)
    RGB( 36, 143, 143), // Ru (44)      rgb(36, 143, 143)
    RGB( 10, 125, 140), // Rh (45)      rgb(10, 125, 140)
    RGB(  0, 105, 133), // Pd (46)      rgb(0, 105, 133)
    RGB(128, 128, 144), // Ag (47)      rgb(128, 128, 144)
    RGB(255, 217, 143), // Cd (48)      rgb(255, 217, 143)
    RGB(166, 117, 115), // In (49)      rgb(166, 117, 115)
    RGB(102, 128, 128), // Sn (50)      rgb(102, 128, 128)
    RGB(158,  99, 181), // Sb (51)      rgb(158, 99, 181)
    RGB(212, 122,   0), // Te (52)      rgb(212, 122, 0)
    RGB(148,   0, 148), // I (53)       rgb(148, 0, 148)
    RGB( 66, 158, 176), // Xe (54)      rgb(66, 158, 176)
    RGB( 87,  23, 143), // Cs (55)      rgb(87, 23, 143)
    RGB(255, 165,   0), // Ba (56)      rgb(255, 165, 0)
    RGB(112, 212, 255), // La (57)      rgb(112, 212, 255)
    RGB(255, 255, 199), // Ce (58)      rgb(255, 255, 199)
    RGB(217, 255, 199), // Pr (59)      rgb(217, 255, 199)
    RGB(199, 255, 199), // Nd (60)      rgb(199, 255, 199)
    RGB(163, 255, 199), // Pm (61)      rgb(163, 255, 199)
    RGB(143, 255, 199), // Sm (62)      rgb(143, 255, 199)
    RGB( 97, 255, 199), // Eu (63)      rgb(97, 255, 199)
    RGB( 69, 255, 199), // Gd (64)      rgb(69, 255, 199)
    RGB( 48, 255, 199), // Tb (65)      rgb(48, 255, 199)
    RGB( 31, 255, 199), // Dy (66)      rgb(31, 255, 199)
    RGB(  0, 255, 156), // Ho (67)      rgb(0, 255, 156)
    RGB(  0, 230, 117), // Er (68)      rgb(0, 230, 117)
    RGB(  0, 212,  82), // Tm (69)      rgb(0, 212, 82)
    RGB(  0, 191,  56), // Yb (70)      rgb(0, 191, 56)
    RGB(  0, 171,  36), // Lu (71)      rgb(0, 171, 36)
    RGB( 77, 194, 255), // Hf (72)      rgb(77, 194, 255)
    RGB( 77, 166, 255), // Ta (73)      rgb(77, 166, 255)
    RGB( 33, 148, 214), // W (74)       rgb(33, 148, 214)
    RGB( 38, 125, 171), // Re (75)      rgb(38, 125, 171)
    RGB( 38, 102, 150), // Os (76)      rgb(38, 102, 150)
    RGB( 23,  84, 135), // Ir (77)      rgb(23, 84, 135)
    RGB(208, 208, 224), // Pt (78)      rgb(208, 208, 224)
    RGB(255, 209,  35), // Au (79)      rgb(255, 209, 35)
    RGB(184, 184, 208), // Hg (80)      rgb(184, 184, 208)
    RGB(166,  84,  77), // Tl (81)      rgb(166, 84, 77)
    RGB( 87,  89,  97), // Pb (82)      rgb(87, 89, 97)
    RGB(158,  79, 181), // Bi (83)      rgb(158, 79, 181)
    RGB(171,  92,   0), // Po (84)      rgb(171, 92, 0)
    RGB(117,  79,  69), // At (85)      rgb(117, 79, 69)
    RGB( 66, 130, 150), // Rn (86)      rgb(66, 130, 150)
    RGB( 66,   0, 102), // Fr (87)      rgb(66, 0, 102)
    RGB(  0, 125,   0), // Ra (88)      rgb(0, 125, 0)
    RGB(112, 171, 250), // Ac (89)      rgb(112, 171, 250)
    RGB(  0, 186, 255), // Th (90)      rgb(0, 186, 255)
    RGB(  0, 161, 255), // Pa (91)      rgb(0, 161, 255)
    RGB(  0, 143, 255), // U (92)       rgb(0, 143, 255)
    RGB(  0, 128, 255), // Np (93)      rgb(0, 128, 255)
    RGB(  0, 107, 255), // Pu (94)      rgb(0, 107, 255)
    RGB( 84,  92, 242), // Am (95)      rgb(84, 92, 242)
    RGB(120,  92, 227), // Cm (96)      rgb(120, 92, 227)
    RGB(138,  79, 227), // Bk (97)      rgb(138, 79, 227)
    RGB(161,  54, 212), // Cf (98)      rgb(161, 54, 212)
    RGB(179,  31, 212), // Es (99)      rgb(179, 31, 212)
    RGB(179,  31, 186), // Fm (100)     rgb(179, 31, 186)
    RGB(179,  13, 166), // Md (101)     rgb(179, 13, 166)
    RGB(189,  13, 135), // No (102)     rgb(189, 13, 135)
    RGB(199,   0, 102), // Lr (103)     rgb(199, 0, 102)
    RGB(204,   0,  89), // Rf (104)     rgb(204, 0, 89)
    RGB(209,   0,  79), // Db (105)     rgb(209, 0, 79)
    RGB(217,   0,  69), // Sg (106)     rgb(217, 0, 69)
    RGB(224,   0,  56), // Bh (107)     rgb(224, 0, 56)
    RGB(230,   0,  46), // Hs (108)     rgb(230, 0, 46)
    RGB(235,   0,  38), // Mt (109)     rgb(235, 0, 38)
    RGB(240,   0,  34), // Ds (110)     rgb(240, 0, 34)
    RGB(246,   0,  32), // Rg (111)     rgb(246, 0, 32)
    RGB(248,   0,  30), // Cn (112)     rgb(248, 0, 30)
    RGB(250,   0,  28), // Nh (113)     rgb(250, 0, 28)
    RGB(252,   0,  26), // Fl (114)     rgb(252, 0, 26)
    RGB(253,   0,  24), // Mc (115)     rgb(253, 0, 24)
    RGB(254,   0,  22), // Lv (116)     rgb(254, 0, 22)
    RGB(255,   0,  20), // Ts (117)     rgb(255, 0, 20)
    RGB(255,   0,  18), // Og (118)     rgb(255, 0, 18)
};

// Structure to hold string and color to represent predefined color mappings for e.g. amino acids
typedef struct {
    const char* name;
    uint32_t color;
} named_color_t;

// Table for predefined amino acid color mappings
named_color_t amino_acid_colors[] = {
    { "ALA", RGB(200, 200, 200) },  // Alanine                      rgb(200, 200, 200)
    { "ARG", RGB(20, 90, 255) },    // Arginine                     rgb(20, 90, 255)
    { "ASN", RGB(220, 220, 220) },  // Asparagine                   rgb(220, 220, 220)
    { "ASP", RGB(230, 10, 10) },    // Aspartic acid                rgb(230, 10, 10)
    { "CYS", RGB(255, 255, 0) },    // Cysteine                     rgb(255, 255, 0)
    { "GLN", RGB(220, 220, 220) },  // Glutamine                    rgb(220, 220, 220)
    { "GLU", RGB(230, 10, 10) },    // Glutamic acid                rgb(230, 10, 10)
    { "GLY", RGB(235, 235, 235) },  // Glycine                      rgb(235, 235, 235)
    { "HIS", RGB(130, 130, 210) },  // Histidine                    rgb(130, 130, 210)
    { "ILE", RGB(200, 150, 100) },  // Isoleucine                   rgb(200, 150, 100)
    { "LEU", RGB(200, 150, 100) },  // Leucine                      rgb(200, 150, 100)
    { "LYS", RGB(20, 90, 255) },    // Lysine                       rgb(20, 90, 255)
    { "MET", RGB(255, 150, 50) },   // Methionine                   rgb(255, 150, 50)
    { "PHE", RGB(50, 50, 170) },    // Phenylalanine                rgb(50, 50, 170)
    { "PRO", RGB(220, 150, 130) },  // Proline                      rgb(220, 150, 130)
    { "SER", RGB(250, 150, 0) },    // Serine                       rgb(250, 150, 0)
    { "THR", RGB(250, 150, 0) },    // Threonine                    rgb(250, 150, 0)
    { "TRP", RGB(180, 90, 180) },   // Tryptophan                   rgb(180, 90, 180)
    { "TYR", RGB(50, 50, 170) },    // Tyrosine                     rgb(50, 50, 170)
    { "VAL", RGB(200, 150, 100) },  // Valine                       rgb(200, 150, 100)
    { "SEC", RGB(255, 255, 0) },    // Selenocysteine               rgb(255, 255, 0)

    // CCD
    { "MSE", RGB(255, 150, 50) },   // Selenomethionine             rgb(255, 150, 50)
    { "SEP", RGB(250, 150, 0) },    // Phosphoserine                rgb(250, 150, 0)
    { "TPO", RGB(250, 150, 0) },    // Phosphothreonine             rgb(250, 150, 0)
    { "PTR", RGB(250, 150, 0) },    // Phosphotyrosine              rgb(250, 150, 0)

    // Charmm
    { "HSD", RGB(130, 130, 210) },  // Histidine (delta)            rgb(130, 130, 210)
    { "HSP", RGB(130, 130, 210) },  // Histidine (epsilon)          rgb(130, 130, 210)
    { "LSN", RGB(220, 220, 220) },  // Asparagine (protonated)      rgb(220, 220, 220)
    { "ASPP", RGB(230, 10, 10) },   // Aspartic acid (protonated)   rgb(230, 10, 10)
    { "GLUP", RGB(230, 10, 10) },   // Glutamic acid (protonated)   rgb(230, 10, 10)

    // Amber
    { "HID", RGB(130, 130, 210) },  // Histidine (delta)            rgb(130, 130, 210)
    { "HIE", RGB(130, 130, 210) },  // Histidine (epsilon)          rgb(130, 130, 210)
    { "HIP", RGB(130, 130, 210) },  // Histidine (protonated)       rgb(130, 130, 210)
    { "LYN", RGB(20, 90, 255) },    // Lysine (deprotonated)        rgb(20, 90, 255)
    { "ASH", RGB(230, 10, 10) },    // Aspartic acid (protonated)   rgb(230, 10, 10)
    { "GLH", RGB(230, 10, 10) },    // Glutamic acid (protonated)   rgb(230, 10, 10)
};

// Table for predefined colors for nucleic acids

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

static const uint64_t element_covalent_mask[2] = {
    0x003E0000001FC7E6ULL, // atomic numbers 1–63
    0x0000000000000000ULL, // atomic numbers 64–127 (unused)
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

static inline bool can_form_covalent_bond(size_t atomic_nr) {
    // 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54
    switch (atomic_nr) {
        // Nonmetals
        case H:
        case C:
        case N:
        case O:
        case F:
        case P:
        case S:
        case Cl:
        case Se:
        case Br:
        case I:

        // Noble gases
        case He:
        case Ne:
        case Ar:
        case Kr:
        case Xe:

        // Metalloids
        case B:
        case Si:
        case Ge:
        case As:
        case Sb:
        case Te:

        // Transition metals
        case Ti:
        case V:
        case Cr:
        case Mn:
        case Fe:
        case Co:
        case Ni:
        case Cu:
        case Zn:
        case Nb:
        case Mo:
        case Tc:
        case Ru:
        case Rh:
        case Pd:
        case Ag:
        case Ta:
        case W:
        case Re:
        case Os:
        case Ir:
        case Pt:
        case Au:

            return true;
        default:
            return false;
    }
}

static inline bool is_metal_donor(size_t atomic_nr) {
    switch(atomic_nr) {
    case O:
    case N:
    case S:
    case Cl:
    case Br:
    case F:
        return true;
    default:
        return false;
    }
}

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
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4", "TIP5", "W", "DOD", "D30", "SPC" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

static const char* common_ions[] = { "NA", "K", "CA", "MG", "ZN", "CL", "F", "MN", "FE", "CU", "CO", "NI", "CD", "BR", "I", "CS", "SR"};

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

static void sort_radix_inplace_uint32(uint32_t* data, size_t count, uint32_t* temp) {
    // Histograms for each pass
    uint32_t hist[256][4] = {0};

    // Populate histogram
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = data[i];
        hist[(id >> 0)  & 255][0]++;
        hist[(id >> 8)  & 255][1]++;
        hist[(id >> 16) & 255][2]++;
        hist[(id >> 24) & 255][3]++;
    }

    // Prefix sum
    md_128i sum = md_mm_setzero_si128();
    for (size_t i = 0; i < 256; ++i) {
        md_128i val = md_mm_loadu_si128(hist[i]);
        md_mm_storeu_si128(hist[i], sum);
        sum = md_mm_add_epi32(sum, val);
    }

    // 4-pass 8-bit radix sort
    radix_pass_8(temp, data, count, hist, 0);
    radix_pass_8(data, temp, count, hist, 1);
    radix_pass_8(temp, data, count, hist, 2);
    radix_pass_8(data, temp, count, hist, 3);
}

static inline void radix_pass_idx_8(uint32_t* dst, const uint32_t* src, const uint32_t* keys, size_t count, const uint32_t hist[256], int pass) {
    int bitoff = pass * 8;

    // Make a working copy of the counts for prefix sums
    uint32_t offsets[256];
    MEMCPY(offsets, hist, sizeof(uint32_t) * 256);

    // Prefix sum on offsets
    uint32_t sum = 0;
    for (int i = 0; i < 256; ++i) {
        uint32_t v = offsets[i];
        offsets[i] = sum;
        sum += v;
    }

    // Scatter indices to dst using offsets
    for (size_t i = 0; i < count; ++i) {
        uint32_t idx = src[i];
        uint32_t key = keys[idx];
        uint32_t id  = (key >> bitoff) & 255;
        dst[offsets[id]++] = idx;
    }
}

static void sort_radix_uint32(uint32_t* out_indices, const uint32_t* keys, size_t count, uint32_t* tmp_indices) {
    // Initialize indices
    for (size_t i = 0; i < count; ++i) {
        out_indices[i] = (uint32_t)i;
    }

    // Precompute histograms for each pass (counts only)
    uint32_t hist[4][256] = {0};
    for (size_t i = 0; i < count; ++i) {
        uint32_t key = keys[i];
        hist[0][(key >> 0)  & 255]++;
        hist[1][(key >> 8)  & 255]++;
        hist[2][(key >> 16) & 255]++;
        hist[3][(key >> 24) & 255]++;
    }

    // 4-pass radix sort using precomputed counts
    radix_pass_idx_8(tmp_indices, out_indices, keys, count, hist[0], 0);
    radix_pass_idx_8(out_indices, tmp_indices, keys, count, hist[1], 1);
    radix_pass_idx_8(tmp_indices, out_indices, keys, count, hist[2], 2);
    radix_pass_idx_8(out_indices, tmp_indices, keys, count, hist[3], 3);
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

// Confusing procedure name
// Extracts a graph from an atom index range with supplied atom types

typedef enum {
    VERTEX_TYPE_MAPPING_MODE_ATOMIC_NUMBER = 0,
    VERTEX_TYPE_MAPPING_MODE_ATOM_TYPE_IDX = 1,
} vertex_type_mapping_mode_t;

static inline uint8_t vertex_type_from_atom_idx(const md_atom_data_t* atom, size_t atom_idx, vertex_type_mapping_mode_t vertex_type_mapping) {
    switch (vertex_type_mapping) {
    case VERTEX_TYPE_MAPPING_MODE_ATOMIC_NUMBER:
        return md_atom_atomic_number(atom, atom_idx);
    case VERTEX_TYPE_MAPPING_MODE_ATOM_TYPE_IDX:
        return (uint8_t)md_atom_type_idx(atom, atom_idx);
    default:
        ASSERT(false);
    }
    return 0;
}

static graph_t extract_graph(const md_system_t* sys, const int indices[], size_t count, vertex_type_mapping_mode_t vertex_mapping, md_allocator_i* vm_arena) {
    ASSERT(sys);
    ASSERT(indices);
    ASSERT(vm_arena);

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
    md_hashmap32_t global_to_local = { .allocator = vm_arena };
    md_hashmap_reserve(&global_to_local, count);

    for (int i = 0; i < (int)count; ++i) {
        int idx = indices[i];
        md_hashmap_add(&global_to_local, (uint64_t)idx, i);
        graph.vertex_type[i] = vertex_type_from_atom_idx(&sys->atom, idx, vertex_mapping);
        graph.atom_idx_map[i] = idx;
    }

    // Only store edges which point to vertices within the graph as this will be used later as a traversal template
    for (size_t i = 0; i < count; ++i) {
        int idx = indices[i];
        // Translate the global atom indices to local structure indices
        uint64_t edge_data_arr[8];
        uint32_t length = 0;

        md_bond_iter_t it = md_bond_iter(&sys->bond, idx);
        while (md_bond_iter_has_next(&it)) {
            uint32_t bond_idx = md_bond_iter_bond_index(&it);
            uint32_t atom_idx = md_bond_iter_atom_index(&it);
            uint32_t flags    = md_bond_iter_bond_flags(&it);
            uint32_t* local_idx = md_hashmap_get(&global_to_local, atom_idx);
            if (local_idx) {
                // Only commit the edge if it is referring to a local index within the graph
                edge_data_arr[length++] = ((uint64_t)bond_idx << 32) | ((uint64_t)flags << 24) | (uint32_t)(*local_idx);
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

bool md_util_resname_nucleotide(str_t str) {
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

bool md_util_resname_ion(str_t str) {
    return find_str_in_cstr_arr(NULL, str, common_ions, ARRAY_SIZE(common_ions));
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

md_element_t md_util_element_lookup(str_t str, bool ignore_case) {
    if (str.len == 1 || str.len == 2) {
        for (md_element_t i = 0; i < ARRAY_SIZE(element_symbols); ++i) {
            str_t sym = element_symbols[i];
            if (str_len(sym) != str_len(str)) continue;
            if (ignore_case) {
                if (str_eq_ignore_case(str, sym)) return i;   
            } else {
                if (str_eq(str, sym)) return i;
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

static inline bool cmp(str_t str, const char* ref, int len) {
    // We compare against len + 1 to include the zero terminator character
    return MEMCMP(str.ptr, ref, len + 1) == 0;
}

static bool md_util_protein_backbone_atoms_extract(md_protein_backbone_atoms_t* out_backbone_atoms, const md_atom_data_t* atom_data, md_urange_t atom_range) {
    ASSERT(atom_data);
    int len = atom_range.end - atom_range.beg;
    if (len < 4) return false;

    static const uint32_t all_bits = 1 | 2 | 4 | 8 | 16;
    uint32_t bits = 0;
    md_protein_backbone_atoms_t bb;
    MEMSET(&bb, -1, sizeof(bb));
    for (uint32_t i = atom_range.beg; i < atom_range.end; ++i) {
        str_t id = md_atom_name(atom_data, i);
        if (str_empty(id)) continue;

        if (!(bits & 1)  && cmp(id, "N",  1)) { bb.n  = (md_atom_idx_t)i; bits |= 1;  continue; }
        if (!(bits & 16) && cmp(id, "HN", 2)) { bb.hn = (md_atom_idx_t)i; bits |= 16; continue; }
        if (!(bits & 2)  && cmp(id, "CA", 2)) { bb.ca = (md_atom_idx_t)i; bits |= 2;  continue; }
        if (!(bits & 4)  && cmp(id, "C",  1)) { bb.c  = (md_atom_idx_t)i; bits |= 4;  continue; }
        if (!(bits & 8)) {
            if (cmp(id, "O", 1)     ||
                cmp(id, "O1", 2)    ||
                cmp(id, "OT1", 3)   ||
                cmp(id, "OC1", 3))
            {
                bb.o = (md_atom_idx_t)i; bits |= 8; continue;
            }
        }

        // Check if done
        if (bits == all_bits) break;
    }

    // HN is optional
    static const uint32_t req_bits = 1 | 2 | 4 | 8;

    if (out_backbone_atoms) *out_backbone_atoms = bb;
    return (bits & req_bits) == req_bits;
}

// There are several conventions to label Nucleic backbone atoms (e.g. C3', C3* or just C3)
static bool cmp_nuc_atom(str_t id, const char* base, size_t len) {
    if (!str_eq_cstr_n(id, base, len)) return false;
    char next = str_ptr(id)[len];
    return next == 0 || next == '\'' || next == '*';
}

static bool md_util_nucleic_backbone_atoms_extract(md_nucleic_backbone_atoms_t* out_backbone_atoms, const md_atom_data_t* atom_data, md_urange_t atom_range) {
    ASSERT(atom_data);
    int len = atom_range.end - atom_range.beg;
    if (len < 6) return false;

    static const uint32_t all_bits = 1 | 2 | 4 | 8 | 16 | 32;

    uint32_t bits = 0;
    md_nucleic_backbone_atoms_t bb = {0};
    for (uint32_t i = atom_range.beg; i < atom_range.end; ++i) {
        str_t id = md_atom_name(atom_data, i);
        if (str_empty(id)) continue;

        if (!(bits & 1)  && cmp(id, "P",   1))               { bb.p  = i; bits |= 1;  continue; }
        if (!(bits & 2)  && cmp_nuc_atom(id, "O5", 2))       { bb.o5 = i; bits |= 2;  continue; }
        if (!(bits & 4)  && cmp_nuc_atom(id, "C5", 2))       { bb.c5 = i; bits |= 4;  continue; }
        if (!(bits & 8)  && cmp_nuc_atom(id, "C4", 2))       { bb.c4 = i; bits |= 8;  continue; }
        if (!(bits & 16) && cmp_nuc_atom(id, "C3", 2))       { bb.c3 = i; bits |= 16; continue; }
        if (!(bits & 32) && (cmp_nuc_atom(id, "O3", 2) || cmp(id, "O3P", 3))) {
            bb.o3 = i; bits |= 32; continue;
        }

        // Check if done
        if (bits == all_bits) break;
    }

    static const uint32_t req_bits = 2 | 4 | 8 | 16 | 32; // No phosphor in first terminal component (starts with O5)

    if (out_backbone_atoms) *out_backbone_atoms = bb;
    return (bits & req_bits) == req_bits;
}

size_t md_util_element_from_mass(md_element_t element[], const float mass[], size_t count) {
    if (!element) {
        MD_LOG_ERROR("element is null");
        return false;
    }

    if (!mass) {
        MD_LOG_ERROR("mass is null");
        return false;
    }

    size_t successful_matches = 0;

    for (size_t i = 0; i < count; ++i) {
        md_element_t elem = 0;
        const float m = mass[i];

        // Use an epsilon which scales with the mass
        const float eps = 0.001f * fmaxf(1.0f, m * 0.01f);

        if (0.0f < m && m != 1.0f) {
            // Linear search for matching atomic mass
            for (size_t j = 1; j < ARRAY_SIZE(element_atomic_mass); ++j) {
                if (fabs(m - element_atomic_mass[j]) < eps) {
                    elem = (md_element_t)j;
                    break;
                } else if (m < element_atomic_mass[j]) {
                    //Masses are sorted, no need to continue searching
                    break;
                }
            }
        }

        element[i] = elem;
        successful_matches += (size_t)(elem != 0);
    }

    return successful_matches;
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

static inline bool zhang_skolnick_ss(const vec3_t* ca_pos, int seq_beg, int seq_end, int i, const float distances[3], float delta) {
    for (int j = MAX(seq_beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= seq_end) continue;
            const vec3_t pos_j = ca_pos[j];
            const vec3_t pos_k = ca_pos[j + k];
            const float d = vec3_distance(pos_j, pos_k);
            if (fabsf(d - distances[k - 2]) > delta) {
                return false;
            }
        }
    }
    return true;
}

static inline bool is_helical(const vec3_t* ca_pos, int seq_beg, int seq_end, int i) {
    const float distances[] = { 5.45f, 5.18f, 6.37f };
    const float delta = 2.1f;
    return zhang_skolnick_ss(ca_pos, seq_beg, seq_end, i, distances, delta);
}

static inline bool is_sheet(const vec3_t* ca_pos, int seq_beg, int seq_end, int i) {
    const float distances[] = { 6.1f, 10.4f, 13.0f };
    const float delta = 1.42f;
    return zhang_skolnick_ss(ca_pos, seq_beg, seq_end, i, distances, delta);
}

// TM-align: a protein structure alignment algorithm based on the Tm-score
// doi:10.1093/nar/gki524

enum {
    TM_ALIGN_VALID = 1,
	TM_ALIGN_HELIX = 2,
	TM_ALIGN_SHEET = 4,
};

bool tm_align(md_secondary_structure_t secondary_structure[], size_t capacity, const float* x, const float* y, const float* z, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {
    if (!secondary_structure) return false;
    if (capacity == 0) return false;

    MEMSET(secondary_structure, 0, capacity * sizeof(md_secondary_structure_t));

    if (!backbone) return false;
    if (!x) return false;
    if (!y) return false;
    if (!z) return false;
    if (!backbone->segment.atoms) return false;
    if (!backbone->range.offset)  return false;

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));
    int* type_mask = md_vm_arena_push(temp_arena, sizeof(int)    * backbone->segment.count);
    vec3_t* ca_pos = md_vm_arena_push(temp_arena, sizeof(vec3_t) * backbone->segment.count);

    for (size_t i = 0; i < backbone->segment.count; ++i) {
        const md_atom_idx_t ca_idx = backbone->segment.atoms[i].ca;
        if (ca_idx >= 0) {
			ca_pos[i] = vec3_set(x[ca_idx], y[ca_idx], z[ca_idx]);
        }
		type_mask[i] = ca_idx >= 0 ? TM_ALIGN_VALID : 0;
	}

    const float helix_distances[] = { 5.45f, 5.18f, 6.37f };
    const float helix_delta = 2.1f;

    const float sheet_distances[] = { 6.1f, 10.4f, 13.0f };
    const float sheet_delta = 1.42f;

    for (size_t bb_idx = 0; bb_idx < backbone->range.count; ++bb_idx) {
        const md_urange_t range = {backbone->range.offset[bb_idx], backbone->range.offset[bb_idx + 1]};
        ASSERT(range.end <= (int)capacity);

        if (range.end - range.beg < 4) {
            continue;
        }

        // Classify residues
        for (uint32_t i = range.beg; i < range.end; ++i) {
            md_secondary_structure_t ss = MD_SECONDARY_STRUCTURE_COIL;

            if (is_sheet(ca_pos, range.beg, range.end, i)) {
                ss = MD_SECONDARY_STRUCTURE_BETA_SHEET;
            }
            else if (is_helical(ca_pos, range.beg, range.end, i)) {
                ss = MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
            }
            secondary_structure[i] = ss;
        }

        // Set squished isolated structures to the surrounding (only for sheets and helices)
        md_secondary_structure_t* ss = secondary_structure;
        for (int i = (int)range.beg + 1; i < (int)range.end - 1; ++i) {
            if (ss[i-1] == ss[i+1] &&
                ss[i-1] != MD_SECONDARY_STRUCTURE_COIL &&
                ss[i+1] != MD_SECONDARY_STRUCTURE_COIL)
            {
                ss[i] = ss[i-1];
            }
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MD_SECONDARY_STRUCTURE_COIL;
        for (int i = (int)range.beg + 1; i < (int)range.end - 1; ++i) {
            if (ss[i] != ss[i-1] &&
                ss[i] != ss[i+1])
            {
                ss[i] = MD_SECONDARY_STRUCTURE_COIL;
            }
        }
        if (ss[range.end - 1] != ss[range.end - 2]) ss[range.end - 1] = MD_SECONDARY_STRUCTURE_COIL;
    }

    md_vm_arena_destroy(temp_arena);

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

    float energy = 0.0f;

    const float min_distance_sq = 0.5f * 0.5f;
    const float min_bond_energy = -9.9f;
	const float energy_cutoff   = -0.5f;

    const vec4_t d2 = {
        vec4_distance_squared(don->N, acc->O),
        vec4_distance_squared(don->H, acc->C),
        vec4_distance_squared(don->N, acc->C),
        vec4_distance_squared(don->H, acc->O),
    };
    
    const int mask_min = vec4_move_mask(vec4_cmp_lt(d2, vec4_set1(min_distance_sq)));
    if (mask_min > 0) {
        energy = min_bond_energy;
    } else {
        const vec4_t inv_d = vec4_rsqrt(d2);
		const vec4_t w = vec4_set(1.0f, 1.0f, -1.0f, -1.0f);
        const float Q  = 27.888f; // 0.084f * 332.0f;
        energy = Q * vec4_reduce_add(vec4_mul(inv_d, w));
    }

    if (energy > energy_cutoff)
		return 0.0f;

    energy = MAX(energy, min_bond_energy);

    dssp_res_hbonds_t* don_hbonds = &res_hbonds[don_idx];
    dssp_res_hbonds_t* acc_hbonds = &res_hbonds[acc_idx];

    if (energy < don_hbonds->acc[0].energy) {
        don_hbonds->acc[1] = don_hbonds->acc[0];
        don_hbonds->acc[0].res_idx = (uint32_t)acc_idx;
        don_hbonds->acc[0].energy  = energy;
    } else if (energy < don_hbonds->acc[1].energy) {
        don_hbonds->acc[1].res_idx = (uint32_t)acc_idx;
        don_hbonds->acc[1].energy  = energy;
    }

    if (energy < acc_hbonds->don[0].energy) {
        acc_hbonds->don[1] = acc_hbonds->don[0];
        acc_hbonds->don[0].res_idx = (uint32_t)don_idx;
        acc_hbonds->don[0].energy  = energy;
    } else if (energy < acc_hbonds->don[1].energy) {
        acc_hbonds->don[1].res_idx = (uint32_t)don_idx;
        acc_hbonds->don[1].energy  = energy;
    }

    return energy;
}

static inline vec4_t estimate_HN(vec4_t N, vec4_t CA, vec4_t C_prev) {
    const float NH_dist = 1.01f;
    vec4_t U = vec4_normalize(vec4_sub(N, CA));
    vec4_t V = vec4_normalize(vec4_sub(N, C_prev));
    return vec4_add(N, vec4_mul_f(vec4_normalize(vec4_add(U, V)), NH_dist));
}

static inline vec4_t estimate_HN_term(vec4_t N, vec4_t CA, vec4_t C) {
    // i == range.beg (N-terminus): estimate previous C by extrapolation
    // C_prev_est = N + (CA - C)
    vec4_t CA_minus_C = vec4_sub(CA, C);
    vec4_t Cprev_est  = vec4_add(N, CA_minus_C);

    // Try to construct H using the standard estimator
    return estimate_HN(N, CA, Cprev_est);
}

static inline bool dssp_test_bond(const dssp_res_hbonds_t res_hbonds[], size_t acceptor, size_t donor) {
    const float max_hbond_energy = -0.5f;
    return (res_hbonds[donor].acc[0].res_idx == acceptor && res_hbonds[donor].acc[0].energy < max_hbond_energy) ||
        (res_hbonds[donor].acc[1].res_idx == acceptor && res_hbonds[donor].acc[1].energy < max_hbond_energy);
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

typedef struct residue_pair_t {
    uint32_t i;
    uint32_t j;
} residue_pair_t;

typedef struct dssp_hbond_energy_user_param_t {
	const uint32_t* res_range_id;
    const dssp_res_coords_t* res_coords;
    dssp_res_hbonds_t* res_hbonds;
	md_array(residue_pair_t)* sheet_candidates;
	md_allocator_i* alloc;
} dssp_hbond_energy_user_param_t;

void dssp_hbond_energy_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    dssp_hbond_energy_user_param_t* param = (dssp_hbond_energy_user_param_t*)user_param;
    for (size_t idx = 0; idx < num_pairs; ++idx) {
        uint32_t i = MIN(i_idx[idx], j_idx[idx]);
        uint32_t j = MAX(i_idx[idx], j_idx[idx]);

		uint32_t ri = param->res_range_id[i];
		uint32_t rj = param->res_range_id[j];

        int d = (int)j - (int)i;

        if (ri == rj && d <= 1) continue;

        calc_hbond_energy(param->res_hbonds, param->res_coords, i, j);
        calc_hbond_energy(param->res_hbonds, param->res_coords, j, i);

        // DSSP excludes very local sequence contacts on the same chain
        if (ri == rj && d <= 2) continue;

		// Ensure that pairs are sorted (i < j) to avoid duplicates
		residue_pair_t pair = { i, j };
		md_array_push(*param->sheet_candidates, pair, param->alloc);
    }
}

void dssp(md_secondary_structure_t out_secondary_structure[], size_t capacity, const float* x, const float* y, const float* z, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {
    ASSERT(out_secondary_structure);
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(cell);
    ASSERT(backbone);

    md_allocator_i* temp_alloc = md_vm_arena_create(GIGABYTES(1));

    size_t backbone_segment_count = backbone->segment.count;
    size_t backbone_range_count   = backbone->range.count;
    const uint32_t* backbone_range_offsets = backbone->range.offset;
    const md_protein_backbone_atoms_t* backbone_atoms = backbone->segment.atoms;

    dssp_res_coords_t* res_coords = md_vm_arena_push(temp_alloc, sizeof(dssp_res_coords_t) * backbone_segment_count);
    dssp_res_hbonds_t* res_hbonds = md_vm_arena_push(temp_alloc, sizeof(dssp_res_hbonds_t) * backbone_segment_count);
    uint32_t* res_range_id = md_vm_arena_push_zero(temp_alloc, sizeof(uint32_t) * backbone_segment_count);
    uint8_t* turn_mask = md_vm_arena_push_zero(temp_alloc, sizeof(uint8_t) * backbone_segment_count);
    float* res_ca_x = md_vm_arena_push(temp_alloc, sizeof(float) * backbone_segment_count);
    float* res_ca_y = md_vm_arena_push(temp_alloc, sizeof(float) * backbone_segment_count);
    float* res_ca_z = md_vm_arena_push(temp_alloc, sizeof(float) * backbone_segment_count);

    md_array(residue_pair_t) sheet_candidates = NULL;
	md_array_ensure(sheet_candidates, 1024, temp_alloc);

    for (uint32_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_urange_t range = { backbone_range_offsets[range_idx], backbone_range_offsets[range_idx + 1] };
        for (uint32_t i = range.beg; i < range.end; ++i) {
            res_range_id[i] = range_idx;
        }
    }
    uint32_t* ss_flags = md_vm_arena_push_zero(temp_alloc, sizeof(uint32_t) * backbone_segment_count);

    for (size_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_irange_t range = {(int)backbone_range_offsets[range_idx], (int)backbone_range_offsets[range_idx + 1]};

        for (int i = range.beg; i < range.end; ++i) {
            md_atom_idx_t ca_idx = backbone_atoms[i].ca;
            md_atom_idx_t n_idx  = backbone_atoms[i].n;
            md_atom_idx_t c_idx  = backbone_atoms[i].c;
            md_atom_idx_t o_idx  = backbone_atoms[i].o;
            md_atom_idx_t hn_idx = backbone_atoms[i].hn;

			res_ca_x[i] = x[ca_idx];
			res_ca_y[i] = y[ca_idx];
			res_ca_z[i] = z[ca_idx];

            res_coords[i].CA = vec4_set(x[ca_idx], y[ca_idx], z[ca_idx], 0);
            res_coords[i].N  = vec4_set(x[n_idx],  y[n_idx],  z[n_idx],  0);
            res_coords[i].C  = vec4_set(x[c_idx],  y[c_idx],  z[c_idx],  0);
            res_coords[i].O  = vec4_set(x[o_idx],  y[o_idx],  z[o_idx],  0);

            if (hn_idx >= 0) {
                res_coords[i].H = vec4_set(x[hn_idx], y[hn_idx], z[hn_idx], 0);
            } else {
                if (i > range.beg)
                    res_coords[i].H = estimate_HN(res_coords[i].N, res_coords[i].CA, res_coords[i - 1].C);
                else
                    res_coords[i].H = estimate_HN_term(res_coords[i].N, res_coords[i].CA, res_coords[i].C);
            }
        }
    }

	// Safe to set to zero, as we will only be comparing energies (which are negative)
	MEMSET(res_hbonds, 0, sizeof(dssp_res_hbonds_t) * backbone_segment_count);

    md_spatial_acc_t acc = { .alloc = temp_alloc };
    md_spatial_acc_init(&acc, res_ca_x, res_ca_y, res_ca_z, NULL, backbone_segment_count, 9.0, cell);

    dssp_hbond_energy_user_param_t user_param = {
		.res_range_id = res_range_id,
        .res_coords = res_coords,
        .res_hbonds = res_hbonds,
		.sheet_candidates = &sheet_candidates,
		.alloc = temp_alloc,
	};

    /*
    // Do N^2 bond energy calculations
    const float min_ca_dist2 = (9.0f * 9.0f);
    for (size_t i = 0; i + 1 < backbone_segment_count; ++i) {
    for (size_t j = i + 1; j < backbone_segment_count; ++j) {
    // Only apply local-sequence exclusion within the same backbone range (chain)
    if (res_range_id[i] == res_range_id[j] && (j - i) <= 1) {
    continue;
    }
    if (vec4_distance_squared(res_coords[i].CA, res_coords[j].CA) > min_ca_dist2) {
    continue;
    }

    calc_hbond_energy(res_hbonds, res_coords, i, j);
    calc_hbond_energy(res_hbonds, res_coords, j, i);
    }
    }
    */
    md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, 9.0, dssp_hbond_energy_cb, &user_param);

    // ---- Per-range passes: helices, turns, bends (set flags only) ----
    // Ground-truth DSSP helix assignment:
    // Two consecutive n-turns define a helix, and then residues i+1..i+n get helix flag.
    for (size_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_urange_t range = {backbone_range_offsets[range_idx], backbone_range_offsets[range_idx + 1]};

        // Bits in turn_mask (local indices): 1 = 3-turn, 2 = 4-turn, 4 = 5-turn
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (i + 3 < range.end && dssp_test_bond(res_hbonds, i, i + 3)) {
                turn_mask[i] |= 1;
            }
            if (i + 4 < range.end && dssp_test_bond(res_hbonds, i, i + 4)) {
                turn_mask[i] |= 2;
            }
            if (i + 5 < range.end && dssp_test_bond(res_hbonds, i, i + 5)) {
                turn_mask[i] |= 4;
            }
        }

        // DSSP helix assignment from consecutive turns:
        // If turn_n(i) && turn_n(i+1) then mark residues (i+1 ... i+n) as helix_n.
        // n=3 => 3_10 helix, n=4 => alpha helix, n=5 => pi helix.
        for (uint32_t i = range.beg; i + 1 < range.end; ++i) {
            const uint8_t a = turn_mask[i];
            const uint8_t b = turn_mask[i + 1];

            // 3_10 helix (G)
            if ((a & 1) && (b & 1)) {
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + 3, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_HELIX_310;
                }
            }

            // alpha helix (H)
            if ((a & 2) && (b & 2)) {
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + 4, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_HELIX_ALPHA;
                }
            }

            // pi helix (I)
            if ((a & 4) && (b & 4)) {
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + 5, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_HELIX_PI;
                }
            }

            // Just turn
            if (a) {
				const int len = (a & 4) ? 5 : (a & 2) ? 4 : 3;
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + len, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_TURN;
                }
			}
        }

        // Bends: unchanged
        for (uint32_t i = range.beg + 2; i + 2 < range.end; ++i) {
            vec4_t v1 = vec4_sub(res_coords[i - 2].CA, res_coords[i].CA);
            vec4_t v2 = vec4_sub(res_coords[i + 2].CA, res_coords[i].CA);
            double angle = vec4_angle(v1, v2);
            if (angle < DEG_TO_RAD(70.0)) {
                for (int k = -2; k <= 2; ++k) {
                    ss_flags[i + k] |= SS_FLAG_BEND;
				}
            }
        }
    }

    // ---- Global pass: sheets / bridges ----

    typedef struct {
        uint32_t i;
        uint32_t j;
		uint64_t chain_pair_id_and_type; // Unique id for the chain_pair from (res_range_id[i], res_range_id[j]) and the bridge type (parallel/antiparallel)
    } dssp_bridge_t;

    // Worst-case: could be large; allocate generously but not insanely.
    dssp_bridge_t* bridges = md_vm_arena_push(temp_alloc, sizeof(dssp_bridge_t) * backbone_segment_count * 8);
    size_t num_bridges = 0;

    // Helper lambdas are not allowed in C; keep as local inline-ish macros
#define SAME_RANGE(a,b) (res_range_id[(a)] == res_range_id[(b)])

#if 1
	md_hashmap32_t parallel_bridge_map     = { .allocator = temp_alloc };
	md_hashmap32_t antiparallel_bridge_map = { .allocator = temp_alloc };

	size_t num_sheet_candidates = md_array_size(sheet_candidates);
    for (size_t c = 0; c < num_sheet_candidates; ++c) {
        uint32_t i = sheet_candidates[c].i;
        uint32_t j = sheet_candidates[c].j;

        // Neighbor existence within each chain/range
        const bool i_has_im1 = (i > 0) && SAME_RANGE(i - 1, i);
        const bool i_has_ip1 = (i + 1 < (uint32_t)backbone_segment_count) && SAME_RANGE(i + 1, i);
        const bool j_has_jm1 = (j > 0) && SAME_RANGE(j - 1, j);
        const bool j_has_jp1 = (j + 1 < (uint32_t)backbone_segment_count) && SAME_RANGE(j + 1, j);

        // Evaluate DSSP bridge patterns (directed H-bonds)
        bool antiparallel = false;
        bool parallel     = false;

        // Antiparallel:
        // (i -> j && j -> i) OR (i-1 -> j+1 && j-1 -> i+1)
        if (dssp_test_bond(res_hbonds, i, j) && dssp_test_bond(res_hbonds, j, i)) {
            antiparallel = true;
        } else if (i_has_im1 && i_has_ip1 && j_has_jm1 && j_has_jp1) {
            if (dssp_test_bond(res_hbonds, i - 1, j + 1) && dssp_test_bond(res_hbonds, j - 1, i + 1)) {
                antiparallel = true;
            }
        }

        // Parallel:
        // (i-1 -> j && j -> i+1) OR (j-1 -> i && i -> j+1)
        if (!antiparallel) {
            if (i_has_im1 && i_has_ip1) {
                if (dssp_test_bond(res_hbonds, i - 1, j) && dssp_test_bond(res_hbonds, j, i + 1)) {
                    parallel = true;
                }
            }
            if (!parallel && j_has_jm1 && j_has_jp1) {
                if (dssp_test_bond(res_hbonds, j - 1, i) && dssp_test_bond(res_hbonds, i, j + 1)) {
                    parallel = true;
                }
            }
        }

        if (!(antiparallel || parallel)) {
            continue;
        }

		uint32_t ri = res_range_id[i];
		uint32_t rj = res_range_id[j];

		uint64_t res_range_id = (uint64_t)MIN(ri, rj) << 32 | MAX(ri, rj);

        uint32_t idx = (uint32_t)(num_bridges++);

        // Record bridge (normalize ordering i<j already true)
        bridges[idx] = (dssp_bridge_t){
            .i = i,
            .j = j,
            .chain_pair_id_and_type = (uint64_t)(parallel ? 1 : 0) << 63 | res_range_id,
        };

        // DSSP marks both residues as "bridge" even before ladder/sheet conversion
        ss_flags[i] |= SS_FLAG_BRIDGE;
        ss_flags[j] |= SS_FLAG_BRIDGE;

		uint64_t key = md_hash64(&bridges[idx], sizeof(dssp_bridge_t), 0);
		if (parallel) {
			md_hashmap_add(&parallel_bridge_map, key, idx);
		} else {
			md_hashmap_add(&antiparallel_bridge_map, key, idx);
		}
    }

#else
    const float min_ca_dist2_sheet = 9.0f * 9.0f;
    for (uint32_t i = 0; i < (uint32_t)backbone_segment_count; ++i) {
        const uint32_t ri = res_range_id[i];

        for (uint32_t j = i + 1; j < (uint32_t)backbone_segment_count; ++j) {
            const uint32_t rj = res_range_id[j];

            // DSSP excludes very local sequence contacts on the same chain
            if (ri == rj && (int)j - (int)i <= 2) continue;

            if (vec4_distance_squared(res_coords[i].CA, res_coords[j].CA) > min_ca_dist2_sheet) continue;

            // Neighbor existence within each chain/range
            const bool i_has_im1 = (i > 0) && SAME_RANGE(i - 1, i);
            const bool i_has_ip1 = (i + 1 < (uint32_t)backbone_segment_count) && SAME_RANGE(i + 1, i);
            const bool j_has_jm1 = (j > 0) && SAME_RANGE(j - 1, j);
            const bool j_has_jp1 = (j + 1 < (uint32_t)backbone_segment_count) && SAME_RANGE(j + 1, j);

            // Evaluate DSSP bridge patterns (directed H-bonds)
            bool antiparallel = false;
            bool parallel     = false;

            // Antiparallel:
            // (i -> j && j -> i) OR (i-1 -> j+1 && j-1 -> i+1)
            if (dssp_test_bond(res_hbonds, i, j) && dssp_test_bond(res_hbonds, j, i)) {
                antiparallel = true;
            } else if (i_has_im1 && i_has_ip1 && j_has_jm1 && j_has_jp1) {
                if (dssp_test_bond(res_hbonds, i - 1, j + 1) && dssp_test_bond(res_hbonds, j - 1, i + 1)) {
                    antiparallel = true;
                }
            }

            // Parallel:
            // (i-1 -> j && j -> i+1) OR (j-1 -> i && i -> j+1)
            if (!antiparallel) {
                if (i_has_im1 && i_has_ip1) {
                    if (dssp_test_bond(res_hbonds, i - 1, j) && dssp_test_bond(res_hbonds, j, i + 1)) {
                        parallel = true;
                    }
                }
                if (!parallel && j_has_jm1 && j_has_jp1) {
                    if (dssp_test_bond(res_hbonds, j - 1, i) && dssp_test_bond(res_hbonds, i, j + 1)) {
                        parallel = true;
                    }
                }
            }

            if (!(antiparallel || parallel)) {
                continue;
            }

		    uint64_t res_range_id = (uint64_t)MIN(ri, rj) << 32 | MAX(ri, rj);

            // Record bridge (normalize ordering i<j already true)
            bridges[num_bridges++] = (dssp_bridge_t){
                .i = i,
                .j = j,
                .chain_pair_id_and_type = (uint64_t)(parallel ? 1 : 0) << 63 | res_range_id,
            };

            // DSSP marks both residues as "bridge" even before ladder/sheet conversion
            ss_flags[i] |= SS_FLAG_BRIDGE;
            ss_flags[j] |= SS_FLAG_BRIDGE;
        }
    }
#endif
#undef SAME_RANGE

#if 1
	// ---- Build ladders by directly looking up bridges in hashmaps ----
    // We mark as SHEET if a ladder has >= 2 bridges (i.e. spans >= 3 residues on at least one strand).
    bool* used = md_vm_arena_push_zero(temp_alloc, sizeof(bool) * num_bridges);
    for (size_t b = 0; b < num_bridges; ++b) {
        if (used[b]) continue;
        const dssp_bridge_t seed = bridges[b];
        used[b] = true;
        // Ladder endpoints
        uint32_t i_min = seed.i, i_max = seed.i;
        uint32_t j_min = seed.j, j_max = seed.j;
        // Chain pair and type
        const uint64_t chain_pair_id_and_type = seed.chain_pair_id_and_type;
        const bool is_parallel = (chain_pair_id_and_type >> 63) != 0;
        const uint64_t chain_pair_id = chain_pair_id_and_type & 0x7FFFFFFFFFFFFFFF;
        md_hashmap32_t* bridge_map = is_parallel ? &parallel_bridge_map : &antiparallel_bridge_map;

        bool grew = true;
        while (grew) {
            grew = false;
            // Look for the next bridge in the ladder by checking the expected next positions based on the ladder type
            uint32_t next_i_min, next_i_max, next_j_min, next_j_max;
            if (is_parallel) {
                // Parallel ladder: (i+1, j+1)
                next_i_min = i_max + 1;
                next_j_min = j_max + 1;
                next_i_max = i_min - 1;
                next_j_max = j_min - 1;
            } else {
                // Antiparallel ladder: (i+1, j-1)
                next_i_min = i_max + 1;
                next_j_min = j_min - 1;
                next_i_max = i_min - 1;
                next_j_max = j_max + 1;
            }
            // Check both possible extensions of the ladder
            dssp_bridge_t query_bridge_1 = { .i = next_i_min, .j = next_j_min, .chain_pair_id_and_type = chain_pair_id_and_type };
            dssp_bridge_t query_bridge_2 = { .i = next_i_max, .j = next_j_max, .chain_pair_id_and_type = chain_pair_id_and_type };
			uint64_t key_1 = md_hash64(&query_bridge_1, sizeof(dssp_bridge_t), 0);
			uint64_t key_2 = md_hash64(&query_bridge_2, sizeof(dssp_bridge_t), 0);

            uint32_t* found_idx_1 = md_hashmap_get(bridge_map, key_1);
            uint32_t* found_idx_2 = md_hashmap_get(bridge_map, key_2);
            if (found_idx_1) {
                const dssp_bridge_t next_bridge = bridges[*found_idx_1];
                used[*found_idx_1] = true;
                i_min = MIN(i_min, next_bridge.i);
                i_max = MAX(i_max, next_bridge.i);
                j_min = MIN(j_min, next_bridge.j);
                j_max = MAX(j_max, next_bridge.j);
                grew = true;
            } else if (found_idx_2) {
                const dssp_bridge_t next_bridge = bridges[*found_idx_2];
                used[*found_idx_2] = true;
                i_min = MIN(i_min, next_bridge.i);
                i_max = MAX(i_max, next_bridge.i);
                j_min = MIN(j_min, next_bridge.j);
                j_max = MAX(j_max, next_bridge.j);
                grew = true;
            }
        }
        const uint32_t len_i = i_max - i_min + 1;
        const uint32_t len_j = j_max - j_min + 1;
        if (len_i >= 3 || len_j >= 3) {
            for (uint32_t k = i_min; k <= i_max; ++k) ss_flags[k] |= SS_FLAG_SHEET;
            for (uint32_t k = j_min; k <= j_max; ++k) ss_flags[k] |= SS_FLAG_SHEET;
        }
	}
#else
    // ---- Build ladders by chaining bridges ----
    // We mark as SHEET if a ladder has >= 2 bridges (i.e. spans >= 3 residues on at least one strand).
    bool* used = md_vm_arena_push_zero(temp_alloc, sizeof(bool) * num_bridges);

    for (size_t b = 0; b < num_bridges; ++b) {
        if (used[b]) continue;

        const dssp_bridge_t seed = bridges[b];
		int seed_type = (seed.chain_pair_id_and_type >> 63) != 0 ? 1 : 0;
        used[b] = true;

        // Ladder endpoints
        uint32_t i_min = seed.i, i_max = seed.i;
        uint32_t j_min = seed.j, j_max = seed.j;

        // Chain constraints: ladders should stay within the same chain pair
        const uint32_t ri = res_range_id[seed.i];
        const uint32_t rj = res_range_id[seed.j];

        bool grew = true;
        while (grew) {
            grew = false;

            for (size_t c = 0; c < num_bridges; ++c) {
                if (used[c]) continue;

                const dssp_bridge_t br = bridges[c];
				int br_type = (br.chain_pair_id_and_type >> 63) != 0 ? 1 : 0;
                if (br_type != seed_type) continue;

                // Must be between the same chain pair (order-insensitive)
                const uint32_t bri = res_range_id[br.i];
                const uint32_t brj = res_range_id[br.j];
                if (!((bri == ri && brj == rj) || (bri == rj && brj == ri))) continue;

                bool connects = false;

                if (seed_type == 1) {
                    // Parallel ladder: (i+1, j+1)
                    connects = (br.i == i_max + 1 && br.j == j_max + 1) ||
                        (br.i + 1 == i_min && br.j + 1 == j_min);
                } else {
                    // Antiparallel ladder: (i+1, j-1)
                    connects = (br.i == i_max + 1 && br.j + 1 == j_min) ||
                        (br.i + 1 == i_min && br.j == j_max + 1);
                }

                if (connects) {
                    used[c] = true;
                    i_min = MIN(i_min, br.i);
                    i_max = MAX(i_max, br.i);
                    j_min = MIN(j_min, br.j);
                    j_max = MAX(j_max, br.j);
                    grew = true;
                }
            }
        }

        const uint32_t len_i = i_max - i_min + 1;
        const uint32_t len_j = j_max - j_min + 1;

        if (len_i >= 3 || len_j >= 3) {
            for (uint32_t k = i_min; k <= i_max; ++k) ss_flags[k] |= SS_FLAG_SHEET;
            for (uint32_t k = j_min; k <= j_max; ++k) ss_flags[k] |= SS_FLAG_SHEET;
        }
    }
#endif

    // ---- Final resolution: convert flags to enum with priority ----
    // Priority: HELIX (PI > 310 > ALPHA) > SHEET > BRIDGE > TURN > BEND > COIL
    MEMSET(out_secondary_structure, 0, sizeof(md_secondary_structure_t) * backbone_segment_count);

    for (size_t i = 0; i < backbone_segment_count; ++i) {
        uint32_t f = ss_flags[i];

        if (f & SS_FLAG_HELIX_PI) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_PI;
        } else if (f & SS_FLAG_HELIX_ALPHA) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
        } else if (f & SS_FLAG_HELIX_310) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_310;
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

bool md_util_backbone_secondary_structure_infer(md_secondary_structure_t secondary_structure[], size_t capacity, const float* x, const float* y, const float* z, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {

    if (!backbone) return false;
    if (capacity < backbone->segment.count) {
        MD_LOG_ERROR("Secondary structure capacity is not sufficient to hold the number of segments in the protein backbone");
        return false;
    }

    // TODO: Use tm_align when coarse grained
#if 0
    return tm_align(secondary_structure, capacity, sys);
#else
    if (secondary_structure && x && y && z && backbone) {
        dssp(secondary_structure, capacity, x, y, z, cell, backbone);
        return true;
    }
    return false;
#endif
}

bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const float* x, const float* y, const float* z, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {
    if (!backbone_angles) return false;
    if (capacity == 0) return false;

    MEMSET(backbone_angles, 0, sizeof(md_backbone_angles_t) * capacity);

    if (!backbone) return false;
    if (!x) return false;
    if (!y) return false;
    if (!z) return false;
    if (!backbone->segment.atoms) return false;

    for (size_t bb_idx = 0; bb_idx < backbone->range.count; ++bb_idx) {
        const md_urange_t range = {backbone->range.offset[bb_idx], backbone->range.offset[bb_idx + 1]};
        ASSERT(range.end <= (int)capacity);

        if (range.end - range.beg < 4) {
            continue;
        }

        for (size_t i = range.beg + 1; i + 1 < range.end; ++i) {
			int c_prev_idx = backbone->segment.atoms[i - 1].c;
			int ca_idx     = backbone->segment.atoms[i].ca;
			int  c_idx     = backbone->segment.atoms[i].c;
			int  n_idx     = backbone->segment.atoms[i].n;
			int n_next_idx = backbone->segment.atoms[i + 1].n;

			vec3_t c_prev = vec3_set(x[c_prev_idx], y[c_prev_idx], z[c_prev_idx]);
            vec3_t n      = vec3_set(x[n_idx], y[n_idx], z[n_idx]);
            vec3_t ca     = vec3_set(x[ca_idx], y[ca_idx], z[ca_idx]);
            vec3_t c      = vec3_set(x[c_idx], y[c_idx], z[c_idx]);
            vec3_t n_next = vec3_set(x[n_next_idx], y[n_next_idx], z[n_next_idx]);

            vec3_t d[4] = {
                vec3_sub(n, c_prev),
                vec3_sub(ca, n),
                vec3_sub(c, ca),
                vec3_sub(n_next, c),
            };

            md_util_min_image_vec3(d, ARRAY_SIZE(d), cell);

            backbone_angles[i].phi = vec3_dihedral_angle(d[0], d[1], d[2]);
            backbone_angles[i].psi = vec3_dihedral_angle(d[1], d[2], d[3]);
        }
    }

    return true;
}

bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const md_system_t* sys) {
    ASSERT(ramachandran_types);
    MEMSET(ramachandran_types, MD_RAMACHANDRAN_TYPE_UNKNOWN, sizeof(md_ramachandran_type_t) * capacity);

    if (capacity == 0) return false;
    if (sys->protein_backbone.segment.count == 0) return false;
    if (sys->comp.count == 0) return false;

    ASSERT(sys->comp.name);
    ASSERT(sys->protein_backbone.segment.comp_idx);

    size_t size = MIN(capacity, sys->protein_backbone.segment.count);

    for (size_t i = 0; i < size; ++i) {
        size_t comp_idx = sys->protein_backbone.segment.comp_idx[i];

        str_t name = md_comp_name(&sys->comp, comp_idx);
        if (str_eq(name, STR_LIT("GLY"))) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GLYCINE;
        } else if (str_eq(name, STR_LIT("PRO"))) {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PROLINE;
            ramachandran_types[i - 1] = MD_RAMACHANDRAN_TYPE_PREPROL;
        } else {
            ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GENERAL;
        }
    }

    return true;
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
    if (elem < Num_Elements)
        return max_valence_bonds[elem];
    return 0;
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
            uint32_t conn_len = (uint32_t)md_bond_conn_count(bond, cur_idx);
            if (conn_len == 1) continue;
            uint32_t conn_i   = bond->conn.offset[cur_idx];
            for (uint32_t j = 0; j < conn_len; ++j) {
                md_atom_idx_t idx = md_bond_conn_atom_idx(bond, conn_i, j);
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
    vertex_t* verts = md_vm_arena_push(temp_arena, sizeof(vertex_t) * num_nodes * 4);
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
        uint64_t elem  = atom_element[bond->conn.atom_idx[i]];
        uint64_t flags = bond->flags[bond->conn.bond_idx[i]];
        key += elem * flags;
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
        md_flags_t flags = (md_flags_t)graph_edge_type(data->pattern, i);
        if (flags) {
            md_atom_idx_t nb = map[graph_edge_vertex_idx(data->pattern, i)];
            uint32_t n_edge_beg = data->neighborhood->edge_offset[0];
            uint32_t n_edge_end = data->neighborhood->edge_offset[1];
            for (uint32_t j = n_edge_beg; j < n_edge_end; ++j) {
                if (graph_edge_vertex_idx(data->neighborhood, j) == nb) {
                    md_bond_idx_t b_idx = data->neighborhood->bond_idx_map[j];
                    data->bond->flags[b_idx] = flags;
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

#if 0
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
    uint8_t* elements = md_vm_arena_push_array(temp_arena, uint8_t, atom->count);

    md_atom_extract_atomic_numbers(elements, 0, atom->count, atom);

    // Identify Sp Sp2 and Sp3 types
    for (size_t i = 0; i < atom->count; ++i) {
        uint32_t   conn_len = (uint32_t)md_bond_conn_count(*bond, i);
        uint32_t   conn_i   = bond->conn.offset[i];

        if (conn_len < 2) continue;
        vec3_t xi = md_atom_coord(atom, i);

        switch (conn_len) {
        case 2: {
            uint32_t j = md_bond_conn_atom_idx(*bond, conn_i, 0);
            uint32_t k = md_bond_conn_atom_idx(*bond, conn_i, 1);
            vec3_t xj = md_atom_coord(atom, j);
            vec3_t xk = md_atom_coord(atom, k);
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
                vec3_t xj = md_atom_coord(atom, j);
                vec3_t xk = md_atom_coord(atom, k);
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
            md_atom_idx_t* ring_atom_idx = md_index_range_beg(*rings, ring_idx);
            bool set_c2_to_sp2 = false;
            if (ring_size == 5) {
                vec3_t c[5];
                for (size_t i = 0; i < 5; ++i) {
                    c[i] = md_atom_coord(atom, ring_atom_idx[i]);
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
                    c[i] = md_atom_coord(atom, ring_atom_idx[i]);
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
                    md_atom_idx_t atom_idx = ring_atom_idx[i];
                    if (md_bond_conn_count(*bond, atom_idx) == 2) {
                        type[atom_idx] = 2;
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
        md_atomic_number_t elem = graph.vertex_type[0];
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
		md_element_t elem_i = md_atom_atomic_number(atom, i);
        uint32_t conn_beg = bond->conn.offset[i];
        uint32_t conn_end = bond->conn.offset[i+1];
        uint32_t conn_len = conn_end - conn_beg;

        if (element_min_connectivity[elem_i] == 0 || conn_len < element_min_connectivity[elem_i] || element_max_connectivity[elem_i] < conn_len) {
            continue;
        }

        uint64_t mask = 0;
        for (uint32_t j = conn_beg; j < conn_end; ++j) {
            md_atom_idx_t idx = bond->conn.atom_idx[j];
			uint64_t elem_j = md_atom_atomic_number(atom, idx);
            mask |= ((uint64_t)1 << (elem_j & 63));
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

        graph_t neighborhood = extract_graph(atom, bond, n_idx, n_len, VERTEX_TYPE_MAPPING_MODE_ATOMIC_NUMBER, temp_arena);

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
                if (type[*it] != 2 || !is_aromatic(elements[*it])) {
                    goto next_ring;
                }
            }

            uint8_t pidx[8];
            int electron_count = 0;
            for (int *it = atom_beg, j = 0; it != atom_end; ++it, ++j) {
                int i = *it;

                ring_pattern_t p = {0};
                p.node_type = elements[i];
                
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
                    md_element_t elem = elements[bond->conn.atom_idx[k]];
                    md_order_t  order = bond->order[bond->conn.bond_idx[k]];
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

            // 8a: Try to resolve ambigous cases, by looking at neighbors
            for (int *it = atom_beg, j = 0; it != atom_end; ++it, ++j) {
                int pi = pidx[j];
                if (pi == 5 || pi == 13) {
                    int prev_pat_idx = (j == 0) ? pidx[ring_size-1] : pidx[j-1];
                    int next_pat_idx = (j == ring_size-1) ? pidx[0] : pidx[j+1];

                    //int prev_atom_idx = it == atom_beg   ? atom_end[-1] : it[-1];
                    //int next_atom_idx = it == atom_end-1 ? atom_beg[ 0] : it[ 1];

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
                    if (elements[bond->conn.atom_idx[j]] == O) {
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
#endif

typedef struct {
    int atom_i;
    int atom_j;
    float dist;
} bond_pair_t;

typedef struct {
	const float* radius; // Bonding radius for each atom
    float k_min;
    float k_max;
    md_array(bond_pair_t)* candidates;
    md_allocator_i* alloc;
} cov_bond_callback_param_t;

static void test_cov_bond_pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    cov_bond_callback_param_t* data = (cov_bond_callback_param_t*)user_param;

    const md_256  v0   = md_mm256_setzero_ps();
    const md_256i v8   = md_mm256_set1_epi32(8);
    const md_256i vlen = md_mm256_set1_epi32((int)num_pairs);

    const md_256 vk_min = md_mm256_set1_ps(data->k_min);
    const md_256 vk_max = md_mm256_set1_ps(data->k_max);

    md_256i vi = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    for (size_t i = 0; i < num_pairs; i += 8) {
        md_256 vmask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(vi, vlen));

        md_256i vi_idx = md_mm256_loadu_epi32(i_idx + i);
        md_256i vj_idx = md_mm256_loadu_epi32(j_idx + i);
        md_256 v_dist2 = md_mm256_loadu_ps(ij_dist2 + i);
        //md_256 v_dist  = md_mm256_sqrt_ps(v_dist2);  // There are alot of gathers here, so we can afford latency of sqrt

        md_256 rad_i   = md_mm256_mask_i32gather_ps(v0, data->radius, vi_idx, vmask, 4);
        md_256 rad_j   = md_mm256_mask_i32gather_ps(v0, data->radius, vj_idx, vmask, 4);
        md_256 rad_sum = md_mm256_add_ps(rad_i, rad_j);

        md_256 r_min_cov = md_mm256_mul_ps(vk_min, rad_sum);
        md_256 r_max_cov = md_mm256_mul_ps(vk_max, rad_sum);

        md_256 r2_min_cov = md_mm256_mul_ps(r_min_cov, r_min_cov);
        md_256 r2_max_cov = md_mm256_mul_ps(r_max_cov, r_max_cov);

        md_256 cmp_min = md_mm256_cmpgt_ps(v_dist2, r2_min_cov);
        md_256 cmp_max = md_mm256_cmplt_ps(v_dist2, r2_max_cov);

        int mask = md_mm256_movemask_ps(md_mm256_and_ps(cmp_min, cmp_max));
        if (mask) {
            size_t new_size = md_array_size(*data->candidates) + popcnt32(mask);
            md_array_ensure(*data->candidates, new_size, data->alloc);

            while (mask) {
                const int bit_idx = ctz32(mask);
                mask = mask & ~(1 << bit_idx);
                const int idx = i + bit_idx;
                bond_pair_t cp = {
                    .atom_i = i_idx[idx],
                    .atom_j = j_idx[idx],
                    .dist   = sqrtf(ij_dist2[idx]),
                };
                md_array_push_no_grow(*data->candidates, cp);
            }
        }
        vi = md_mm256_add_epi32(vi, v8);
    }
}

static inline void test_bb_pair(int atom_i, int atom_j, float cutoff, const float* x, const float* y, const float* z, const md_unitcell_t* cell, md_array(bond_pair_t)* candidates, md_allocator_i* alloc) {
    vec3_t d = vec3_set(
        x[atom_i] - x[atom_j],
        y[atom_i] - y[atom_j],
        z[atom_i] - z[atom_j]
    );
    md_util_min_image_vec3(&d, 1, cell);
    float dist = vec3_length(d);
    if (dist < cutoff) {
        bond_pair_t cp = {
            .atom_i = atom_i,
            .atom_j = atom_j,
            .dist   = dist,
        };
        md_array_push(*candidates, cp, alloc);
    }
}

void md_util_infer_covalent_bonds(md_bond_data_t* bond, const float* x, const float* y, const float* z, const md_unitcell_t* cell, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(bond);
    ASSERT(sys);
    ASSERT(alloc);

    md_bond_data_clear(bond);

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));
    
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing atom field (x/y/z)");
        goto done;
    }

    // Check if system is coarse grained
	bool is_coarse_grained = false;
    for (size_t i = 0; i < sys->atom.type.count; ++i) {
		md_flags_t atype_flags = md_atom_type_flags(&sys->atom.type, i);
		if (atype_flags & MD_FLAG_COARSE_GRAINED) {
			is_coarse_grained = true;
            break;
        }
	}

    size_t num_atoms = sys->atom.count;
    md_array(bond_pair_t) candidates = 0;
    md_array_ensure(candidates, num_atoms, temp_arena);

    if (is_coarse_grained) {
		// In coarse-grained systems, use a simple distance cutoff and only check the backbone atoms in component sequence
        int bb_prev = -1;
		md_flags_t comp_i_flags = md_comp_flags(&sys->comp, 0);
        if ((comp_i_flags & MD_FLAG_AMINO_ACID)) {
			// Get the backbone atom index of the first component
			md_urange_t atom_range = md_comp_atom_range(&sys->comp, 0);
			for (size_t i = atom_range.beg; i < atom_range.end; ++i) {
				md_flags_t atom_flags = md_atom_flags(&sys->atom, i);
				if (atom_flags & MD_FLAG_BACKBONE) {
					bb_prev = (int)i;
                    break;
				}
            }
        }
        // Test consecutive components
        for (size_t ci = 1; ci < sys->comp.count; ++ci) {
			int bb_i = -1;
			md_flags_t comp_j_flags = md_comp_flags(&sys->comp, ci);
            if ((comp_j_flags & MD_FLAG_AMINO_ACID)) {
			    // Get the backbone atom index of the first component
			    md_urange_t atom_range = md_comp_atom_range(&sys->comp, ci);
			    for (size_t i = atom_range.beg; i < atom_range.end; ++i) {
				    md_flags_t atom_i_flags = md_atom_flags(&sys->atom, i);
				    if (atom_i_flags & MD_FLAG_BACKBONE) {
					    bb_i = (int)i;
                        break;
				    }
                }
                // Test internal components
                for (size_t i = atom_range.beg; i + 1 < atom_range.end; ++i) {
                    for (size_t j = atom_range.beg + 1; j < atom_range.end; ++j) {
                        test_bb_pair((int)i, (int)j, 4.0f, x, y, z, cell, &candidates, temp_arena);
                    }
				}
            }
            if (bb_prev >= 0 && bb_i >= 0) {
				test_bb_pair(bb_prev, bb_i, 4.0f, x, y, z, cell, &candidates, temp_arena);
			}
			bb_prev = bb_i;
        }

        size_t num_candidates = md_array_size(candidates);
		md_array_ensure(bond->pairs, num_candidates, alloc);
        md_array_ensure(bond->flags, num_candidates, alloc);

        for (size_t i = 0; i < num_candidates; ++i) {
            md_atom_pair_t pair = {
                candidates[i].atom_i,
                candidates[i].atom_j,
            };
            md_array_push_no_grow(bond->pairs, pair);
            md_array_push_no_grow(bond->flags, MD_BOND_FLAG_COVALENT);
            bond->count += 1;
        }
	} else {
        // Covalent radius sum factors
        static const float k_min   = 0.65f;
        static const float k_cov   = 1.15f;
        static const float k_tight = 1.05f;
        static const float k_coord = 1.35f;
        static const float k_metal = 0.90f;

        double max_cov_rad = 0.0;

        float* cov_radius = md_vm_arena_push_array(temp_arena, float, num_atoms);
        md_atomic_number_t* atomic_nr = md_vm_arena_push_array(temp_arena, md_atomic_number_t, num_atoms);
        for (size_t i = 0; i < num_atoms; ++i) {
            atomic_nr[i]  = md_atom_atomic_number(&sys->atom, i);
            cov_radius[i] = element_covalent_radius(atomic_nr[i]);
            max_cov_rad = MAX(max_cov_rad, cov_radius[i]);
        }

        {
            // Find connections for all atoms
            cov_bond_callback_param_t param = {
                .radius = cov_radius,
                .k_min = k_min,
                .k_max = k_coord,
                .candidates = &candidates,
                .alloc      = temp_arena,
            };

            // Compute a cell size based on the max cov radius within the set
            const double cell_ext = 2.0 * max_cov_rad * k_coord;

            // Build candidate list
            md_timestamp_t ts_start = md_time_current();
            md_spatial_acc_t acc = {.alloc = temp_arena};
            md_spatial_acc_init(&acc, x, y, z, NULL, num_atoms, cell_ext, cell);
            //md_spatial_acc_for_each_pair_in_neighboring_cells(&acc, test_cov_bond_pair_callback, &param);
            md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cell_ext, test_cov_bond_pair_callback, &param);
            md_timestamp_t ts_end = md_time_current();

            double dt_ms = md_time_as_milliseconds(ts_end - ts_start);
            //MD_LOG_DEBUG("Constructed candidate bond list with cell size of %f in %f ms", cell_ext, dt_ms);
        }

        size_t num_candidates = md_array_size(candidates);

        //MD_LOG_DEBUG("Found %zu candidates", num_candidates);

        bond_pair_t* temp_bond_pairs = 0;
        md_flags_t*  temp_bond_flags = 0;
        md_array_ensure(temp_bond_pairs, num_candidates, temp_arena);
        md_array_ensure(temp_bond_flags, num_candidates, temp_arena);
        size_t temp_bond_count = 0;

        for (size_t i = 0; i < num_candidates; ++i) {
            int ai = candidates[i].atom_i;
            int aj = candidates[i].atom_j;
            float d = candidates[i].dist;
            int ei = atomic_nr[ai];
            int ej = atomic_nr[aj];
            int mi = is_metal(ei);
            int mj = is_metal(ej);
            float sum_r = cov_radius[ai] + cov_radius[aj];

            if (!mi && !mj) {
                // Both non-metals, check against k_cov
                if (d < k_cov * sum_r) {
                    md_array_push(temp_bond_pairs, candidates[i], temp_arena);
                    md_array_push(temp_bond_flags, MD_BOND_FLAG_COVALENT, temp_arena);
                    temp_bond_count += 1;
                }
            } else if (mi ^ mj) { // XOR here (either is metal, but not both)
                int non_metal_e = mi ? ej : ei;
                if (is_metal_donor(non_metal_e) && d < k_coord * sum_r) {
                    md_array_push(temp_bond_pairs, candidates[i], temp_arena);
                    md_array_push(temp_bond_flags, MD_BOND_FLAG_METAL | MD_BOND_FLAG_COORDINATE, temp_arena);
                    temp_bond_count += 1;
                    // Potentially strong / partially covalent if d < ~1.2
                    // Coordination bond (@TODO validate by geometrical matching)
                }
            } else {
                // Both metals
                if (d < k_metal * sum_r) {
                    md_array_push(temp_bond_pairs, candidates[i], temp_arena);
                    md_array_push(temp_bond_flags, MD_BOND_FLAG_METAL, temp_arena);
                    temp_bond_count += 1;
                }
            }
        }

            // Populate real bond data
        md_array_ensure(bond->pairs, temp_bond_count, alloc);
        md_array_ensure(bond->flags, temp_bond_count, alloc);
        for (size_t i = 0; i < temp_bond_count; ++i) {
            if (temp_bond_flags[i] == 0) continue;

            md_atom_pair_t pair = {
                temp_bond_pairs[i].atom_i,
                temp_bond_pairs[i].atom_j,
            };
            md_array_push_no_grow(bond->pairs, pair);
            md_array_push_no_grow(bond->flags, temp_bond_flags[i]);
            bond->count += 1;
        }
    }

    size_t num_comp = md_system_comp_count(sys);
    if (num_comp > 0) {
        // Create map from atom to component index
        md_comp_idx_t* atom_comp_idx = md_vm_arena_push(temp_arena, sizeof(md_comp_idx_t) * num_atoms);
        for (size_t i = 0; i < num_comp; ++i) {
            md_urange_t range = md_system_comp_atom_range(sys, i);
            for (uint32_t j = range.beg; j < range.end; ++j) {
                atom_comp_idx[j] = (md_comp_idx_t)i;
            }
        }

        // Mark interresidual bonds
        for (size_t i = 0; i < bond->count; ++i) {
            md_comp_idx_t comp_idx[2] = {
                atom_comp_idx[bond->pairs[i].idx[0]],
                atom_comp_idx[bond->pairs[i].idx[1]],
            };
            if (comp_idx[0] != comp_idx[1]) {
                bond->flags[i] |= MD_BOND_FLAG_INTER;
            }
        }
    }
    
done:
    md_vm_arena_destroy(temp_arena);
}

#define MIN_RES_LEN 4
#define MAX_RES_LEN 25

#define MIN_NUC_LEN 6
#define MAX_NUC_LEN 35

#define MIN_POLYMER_COMP_LEN 2

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

bool md_util_system_infer_comp_flags(md_system_t* sys) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }

    size_t temp_pos = md_temp_get_pos();
    md_allocator_i* temp_alloc = md_get_temp_allocator();
    md_array(int) ambigous_amino_acid_indices = 0;
    md_array(int) ambigous_nucleotide_indices = 0;

    for (size_t i = 0; i < sys->comp.count; ++i) {
        str_t comp_name = md_comp_name(&sys->comp, i);
		md_flags_t comp_flags = md_comp_flags(&sys->comp, i);
        md_urange_t range = md_comp_atom_range(&sys->comp, i);
        size_t len = (size_t)(range.end - range.beg);

        if (comp_flags & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE | MD_FLAG_WATER | MD_FLAG_ION)) {
            // Already assigned
            continue;
		}

        md_protein_backbone_atoms_t prot_atoms = {0};
        md_nucleic_backbone_atoms_t nucl_atoms = {0};
        if (!(sys->comp.flags[i] & MD_FLAG_HETERO)) {
            if (MIN_RES_LEN <= len && len <= MAX_RES_LEN && md_util_protein_backbone_atoms_extract(&prot_atoms, &sys->atom, range)) {
                sys->comp.flags[i] |= MD_FLAG_AMINO_ACID;
                sys->atom.flags[prot_atoms.n]  |= MD_FLAG_BACKBONE;
                sys->atom.flags[prot_atoms.ca] |= MD_FLAG_BACKBONE;
                sys->atom.flags[prot_atoms.c]  |= MD_FLAG_BACKBONE;
                sys->atom.flags[prot_atoms.o]  |= MD_FLAG_BACKBONE;
                if (prot_atoms.hn != -1) {
                    sys->atom.flags[prot_atoms.hn]  |= MD_FLAG_BACKBONE;
                }
                goto done;
            } else if (md_util_resname_amino_acid(comp_name)) {
                sys->comp.flags[i] |= MD_FLAG_AMINO_ACID;
                md_array_push(ambigous_amino_acid_indices, (int)i, temp_alloc);
                goto done;
            }

            if (MIN_NUC_LEN <= len && len <= MAX_NUC_LEN && md_util_nucleic_backbone_atoms_extract(&nucl_atoms, &sys->atom, range)) {
                sys->comp.flags[i] |= MD_FLAG_NUCLEOTIDE;
                if (sys->atom.flags) {
                    if (nucl_atoms.p != -1) {
                        sys->atom.flags[nucl_atoms.p]  |= MD_FLAG_BACKBONE;
                    }
                    sys->atom.flags[nucl_atoms.o5] |= MD_FLAG_BACKBONE;
                    sys->atom.flags[nucl_atoms.c5] |= MD_FLAG_BACKBONE;
                    sys->atom.flags[nucl_atoms.c4] |= MD_FLAG_BACKBONE;
                    sys->atom.flags[nucl_atoms.c3] |= MD_FLAG_BACKBONE;
                    sys->atom.flags[nucl_atoms.o3] |= MD_FLAG_BACKBONE;
                }
                goto done;
            } else if (md_util_resname_nucleotide(comp_name)) {
                sys->comp.flags[i] |= MD_FLAG_NUCLEOTIDE;
                md_array_push(ambigous_nucleotide_indices, (int)i, temp_alloc);
                goto done;
            }
        }

        if (((len == 1 || len == 3) && md_util_resname_water(comp_name)) ||
           (len == 1 && md_util_resname_water(md_atom_name(&sys->atom, range.beg))))
        {
            sys->comp.flags[i] |= MD_FLAG_WATER;
        } else if (len == 1 && (md_util_resname_ion(comp_name) || monatomic_ion_element(md_atom_atomic_number(&sys->atom, range.beg)))) {
            sys->comp.flags[i] |= MD_FLAG_ION;
        }

done:
        // Propagate flags to atoms
        for (unsigned j = range.beg; j < range.end; ++j) {
            sys->atom.flags[j] |= sys->comp.flags[i];
        }
    }

    const size_t max_warnings = 10;
    
    // Print warnings about unsafe assignments
    size_t num_ambigous_amino_acids = md_array_size(ambigous_amino_acid_indices);
    if (num_ambigous_amino_acids > 0) {
        MD_LOG_DEBUG("Structure contains %zu amino acid residue names, but their backbone cannot be mapped from atom labels:", num_ambigous_amino_acids);
        for (size_t i = 0; i < MIN(num_ambigous_amino_acids, max_warnings); ++i) {
            int res_idx = ambigous_amino_acid_indices[i];
            MD_LOG_DEBUG(" - [%d] " STR_FMT " (%d)", res_idx + 1, STR_ARG(md_comp_name(&sys->comp, res_idx)), md_comp_seq_id(&sys->comp, res_idx));
        }
        if (num_ambigous_amino_acids > max_warnings) {
            MD_LOG_DEBUG(" - And %zu more...", num_ambigous_amino_acids - max_warnings);
        }
    }

    size_t num_ambigous_nucleotides = md_array_size(ambigous_nucleotide_indices);
    if (num_ambigous_nucleotides > 0) {
        MD_LOG_DEBUG("Structure contains %zu nucleotide residue names, but their backbone cannot be mapped from atom labels:", num_ambigous_nucleotides);
        for (size_t i = 0; i < MIN(num_ambigous_nucleotides, max_warnings); ++i) {
            int res_idx = ambigous_nucleotide_indices[i];
            MD_LOG_DEBUG(" - [%d] " STR_FMT " (%d)", res_idx + 1, STR_ARG(md_comp_name(&sys->comp, res_idx)), md_comp_seq_id(&sys->comp, res_idx));
        }
        if (num_ambigous_nucleotides > max_warnings) {
            MD_LOG_DEBUG(" - And %zu more...", num_ambigous_nucleotides - max_warnings);
        }
    }
    
    md_temp_set_pos_back(temp_pos);
    return true;
}

void md_util_hydrogen_bond_init(md_hydrogen_bond_data_t* hbond_data, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(hbond_data);
    ASSERT(sys);
    ASSERT(alloc);

    if (sys->bond.count == 0) {
        return;
    }

    hbond_data->candidate.num_donors = 0;
    md_array_shrink(hbond_data->candidate.donors, 0);

    hbond_data->candidate.num_acceptors = 0;
    md_array_shrink(hbond_data->candidate.acceptors, 0);

    size_t num_atoms = md_system_atom_count(sys);

    // Identify donors and acceptors
    for (size_t i = 0; i < num_atoms; ++i) {
        int max_conn = 0;
        int num_of_lone_pairs = 2;
        md_atomic_number_t z_i = md_atom_atomic_number(&sys->atom, i);
        switch (z_i) {
        case MD_Z_N:
        case MD_Z_S:
            max_conn = 3;
            break;
        case MD_Z_O:
            max_conn = 2;
            break;
        default:
            continue;
        }
        
        md_bond_iter_t it = md_bond_iter(&sys->bond, i);
        while (md_bond_iter_has_next(&it)) {
            md_atom_idx_t j = md_bond_iter_atom_index(&it);
            md_atomic_number_t z_j = md_atom_atomic_number(&sys->atom, j);
            if (z_j== MD_Z_H) {
                md_hydrogen_bond_donor_t donor = {(md_atom_idx_t)i, j};
                md_array_push(hbond_data->candidate.donors, donor, alloc);
                hbond_data->candidate.num_donors += 1;
            }
            md_bond_iter_next(&it);
        }
        
        size_t num_conn = md_bond_conn_count(&sys->bond, i);
        if (num_conn <= max_conn) {
            if (z_i == MD_Z_S) {
                num_of_lone_pairs = 4 - num_conn;
            }

            md_hydrogen_bond_acceptor_t acceptor = {(md_atom_idx_t)i, num_of_lone_pairs};
            md_array_push(hbond_data->candidate.acceptors, acceptor, alloc);
            hbond_data->candidate.num_acceptors += 1;
        }
    }

    md_array_ensure(hbond_data->bonds, 2 * hbond_data->candidate.num_donors, alloc);
}

typedef struct hbond_candidate_t {
    float score;
    int acc_idx;
} hbond_candidate_t;

typedef struct hbond_payload_t {
    vec4_t don;
    vec4_t hyd;
    hbond_candidate_t candidates[16];
    size_t num_candidates;
    float min_angle_in_radians;
} hbond_payload_t;

static bool hbond_iter(const md_spatial_hash_elem_t* elem, void* param) {
    hbond_payload_t* data = (hbond_payload_t*)param;
    if (data->num_candidates == ARRAY_SIZE(data->candidates)) {
        return false;
    }

    vec4_t acc = vec4_from_vec3(elem->xyz, 0);
    float angle = vec4_angle(vec4_sub(data->don, data->hyd), vec4_sub(acc, data->hyd));
    float dist  = vec4_length(vec4_sub(data->don, acc));
    if (angle > data->min_angle_in_radians) {
        float h_dist = vec4_length(vec4_sub(data->hyd, acc));
        float score = cosf(PI - angle) / h_dist;
        hbond_candidate_t candidate = {
            .score = score,
            .acc_idx = elem->idx,
        };
        data->candidates[data->num_candidates++] = candidate;
    }

    return true;
}

void md_util_hydrogen_bond_infer(md_hydrogen_bond_data_t* hbond_data, const float* atom_x, const float* atom_y, const float* atom_z,
                                 const md_unitcell_t* unitcell, double max_dist, double min_angle) {
    ASSERT(hbond_data);
    ASSERT(atom_x);
    ASSERT(atom_y);
    ASSERT(atom_z);

    hbond_data->num_bonds = 0;
    md_array_shrink(hbond_data->bonds, 0);

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));

    // Clamp to some form of reasonable value range
    max_dist = CLAMP(max_dist, 1.0, 7.0);
    min_angle = CLAMP(min_angle, 100.0, 180.0);

    vec3_t*  acc_pos = md_vm_arena_push_array(arena, vec3_t,  hbond_data->candidate.num_acceptors);
    uint8_t* acc_cap = md_vm_arena_push_array(arena, uint8_t, hbond_data->candidate.num_acceptors);

    for (size_t i = 0; i < hbond_data->candidate.num_acceptors; ++i) {
        int atom_idx = hbond_data->candidate.acceptors[i].idx;
        acc_pos[i] = vec3_set(atom_x[atom_idx], atom_y[atom_idx], atom_z[atom_idx]);
        acc_cap[i] = (uint8_t)hbond_data->candidate.acceptors[i].num_of_lone_pairs;
    }

    hbond_payload_t payload = {0};
    payload.min_angle_in_radians = DEG_TO_RAD(min_angle);
    const float radius = max_dist;
    
    md_spatial_hash_t* sh = md_spatial_hash_create_vec3(acc_pos, NULL, hbond_data->candidate.num_acceptors, unitcell, arena);

    for (size_t i = 0; i < hbond_data->candidate.num_donors; ++i) {
        int d_idx = hbond_data->candidate.donors[i].d_idx;
        int h_idx = hbond_data->candidate.donors[i].h_idx;
        payload.don = vec4_set(atom_x[d_idx], atom_y[d_idx], atom_z[d_idx], 0);
        payload.hyd = vec4_set(atom_x[h_idx], atom_y[h_idx], atom_z[h_idx], 0);
        payload.num_candidates = 0;

        md_spatial_hash_query(sh, vec3_from_vec4(payload.don), radius, hbond_iter, &payload);

        // Pick top candidates
        if (payload.num_candidates > 0) {
            hbond_candidate_t* arr = payload.candidates;
            // Sort them (insertion sort)
            int n = payload.num_candidates;
            for (int j = 1; j < n; j++) {
                int k = j - 1;

                // descending order
                while (k >= 0 && arr[k].score < arr[j].score) {
                    arr[k + 1] = arr[k];
                    k--;
                }
                arr[k + 1] = arr[j];
            }

            // Try to form bonds
            for (int j = 0; j < n ; ++j) {
                uint32_t a_idx = (uint32_t)arr[j].acc_idx;
                if (acc_cap[a_idx] > 0) {
                    md_hydrogen_bond_pair_t bond = { (uint32_t)i, a_idx };
                    md_array_push_no_grow(hbond_data->bonds, bond);
                    hbond_data->num_bonds += 1;
                    acc_cap[a_idx] -= 1;
                    break;
                }
            }
        }
    }

    md_vm_arena_destroy(arena);
}

// Excel-style chain-id generator: A..Z, AA..ZZ, AAA..
// Notes:
// - Uses only upper-case A-Z to avoid format-specific surprises.
// - kMaxLen limits growth to a safe small size for md_label_t (tune if needed).
// - If overflow occurs, it clamps at the maximum representable length.
static inline md_label_t md_util_chain_id_from_index(size_t idx) {
    // Convert 0-based idx to 1-based Excel-style base-26
    // 0 -> A, 25 -> Z, 26 -> AA, ...
    char buf[8];
    int  len = 0;

    size_t n = idx;
    while (len < (int)ARRAY_SIZE(buf)) {
        size_t q = n / 26;
        size_t r = n % 26;
        buf[len++] = (char)('A' + (int)r);
        if (q == 0) break;
        n = q - 1;
    }

    // Reverse to get most significant first
    for (int i = 0, j = len - 1; i < j; ++i, --j) {
        char t = buf[i]; buf[i] = buf[j]; buf[j] = t;
    }

    return make_label((str_t){ buf, (size_t)len });
}

// @NOTE(Robin): This could certainly be improved to incorporate more characters
// Perhaps first A-Z, then [A-Z]0-9, then AA-ZZ etc.
static inline md_label_t generate_chain_id_from_index(size_t idx) {
    char c = 'A' + (idx % 26);
    str_t str = {&c, 1};
    return make_label(str);
}

static inline md_label_t md_util_next_inst_id(str_t last) {

    // Fallback for empty/invalid: start at "A"
    if (str_empty(last)) return make_label(STR_LIT("A"));

    // Copy and normalize to upper-case A-Z; if anything else, fallback to "A"
    char buf[8];
    int len = (int)MIN(ARRAY_SIZE(buf), str_len(last));
    for (int i = 0; i < len; ++i) {
        char c = last.ptr[i];
        if ('a' <= c && c <= 'z') c = (char)(c - 'a' + 'A');
        if (c < 'A' || c > 'Z') return make_label(STR_LIT("A"));
        buf[i] = c;
    }

    // Increment with carry: Z->A and carry to the left, prepend 'A' on overflow
    int i = len - 1;
    while (i >= 0) {
        if (buf[i] < 'Z') {
            buf[i] = (char)(buf[i] + 1);
            return make_label((str_t){ buf, (size_t)len });
        }
        buf[i] = 'A';
        --i;
    }

    // Overflow: prepend 'A' if there's room, otherwise clamp to all 'Z'
    if (len < (int)ARRAY_SIZE(buf)) {
        // shift right and set first to 'A'
        for (int j = len; j > 0; --j) buf[j] = buf[j - 1];
        buf[0] = 'A';
        len += 1;
        return make_label((str_t){ buf, (size_t)len });
    } else {
        for (int j = 0; j < len; ++j) buf[j] = 'Z';
        return make_label((str_t){ buf, (size_t)len });
    }
}

// Define what size of components we group into same instances
#define MAX_GROUPED_COMP_SIZE 4

bool md_util_system_infer_entity_and_instance(md_system_t* sys, const str_t comp_auth_asym_id[], struct md_allocator_i* alloc) {
    if (!sys || sys->comp.count == 0) {
        MD_LOG_ERROR("Missing system or components");
        return false;
    }
    ASSERT(alloc);

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    bool* connected_to_prev = (bool*)md_vm_arena_push_zero(temp_arena, sys->comp.count * sizeof(bool));
    if (sys->bond.count > 0) {
        // Pass 1: map atoms to components
        md_comp_idx_t* atom_comp_idx = md_vm_arena_push_array(temp_arena, md_comp_idx_t, sys->atom.count);
        MEMSET(atom_comp_idx, 0, sizeof(md_comp_idx_t) * sys->atom.count);
        for (size_t i = 0; i < sys->comp.count; ++i) {
            md_urange_t atom_range = md_comp_atom_range(&sys->comp, i);
            for (uint32_t j = atom_range.beg; j < atom_range.end; ++j) {
                atom_comp_idx[j] = (md_comp_idx_t)i;
            }
        }

        // Pass 2: find connected residues
        for (size_t i = 0; i < sys->bond.count; ++i) {
            md_atom_pair_t pair = sys->bond.pairs[i];
            md_comp_idx_t res_a = atom_comp_idx[pair.idx[0]];
            md_comp_idx_t res_b = atom_comp_idx[pair.idx[1]];
            if (abs(res_a - res_b) == 1) {
                int res_max = MAX(res_a, res_b);
                connected_to_prev[res_max] = true;
            }
        }
    }

    md_array(uint64_t) entity_keys = 0;
    md_hashset_t inst_id_set = { .allocator = temp_arena };

    // Pass 3: Construct instances (sequential ranges of components) either from connected components (== Polymer?) or from just sequential components with the same name
    {
        size_t i = 0;
        while (i < sys->comp.count) {
            size_t j = i + 1;

            str_t  comp_name = md_comp_name(&sys->comp, i);
            md_seq_id_t comp_seq_id = md_comp_seq_id(&sys->comp, i);
            size_t comp_size = md_comp_atom_count(&sys->comp, i);
            md_flags_t comp_flags = md_comp_flags(&sys->comp, i);
            uint64_t entity_key = md_hash64_str(comp_name, 0);
            str_t  comp_auth_id = comp_auth_asym_id ? str_trim(comp_auth_asym_id[i]) : STR_LIT("");

            md_flags_t entity_flags = comp_flags;

            bool is_amino_or_nucleotide = (comp_flags & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE)) != 0;

            bool test_name      = (comp_flags & (MD_FLAG_WATER | MD_FLAG_ION)) && comp_size <= MAX_GROUPED_COMP_SIZE;
            bool test_seq_id    = !is_amino_or_nucleotide && comp_size > MAX_GROUPED_COMP_SIZE;
            bool test_auth_id   = comp_auth_asym_id && !str_empty(comp_auth_id);
            bool test_bond      = !test_auth_id && is_amino_or_nucleotide && connected_to_prev[j];

#if 0
            MD_LOG_DEBUG("Identifying new instance");
            MD_LOG_DEBUG("\t test_name: %i",    (int)test_name);
            MD_LOG_DEBUG("\t test_seq_id: %i",  (int)test_seq_id);
            MD_LOG_DEBUG("\t test_auth_id: %i", (int)test_auth_id);
            MD_LOG_DEBUG("\t test_bond: %i",    (int)test_bond);
#endif

            if (is_amino_or_nucleotide) {
                entity_flags |= MD_FLAG_POLYMER;
            }

            if (test_name || test_seq_id || test_bond || test_auth_id) {
                while (j < sys->comp.count) {
                    str_t comp_name_j = md_comp_name(&sys->comp, j);
                    if (test_name && !str_eq(comp_name, comp_name_j)) break;
                    if (test_seq_id && comp_seq_id != md_comp_seq_id(&sys->comp, j)) break;
                    if (test_bond && !connected_to_prev[j]) break;
                    if (test_auth_id && !str_eq(comp_auth_id, comp_auth_asym_id[j])) break;

                    md_flags_t flags = md_comp_flags(&sys->comp, j);
                    if ((flags ^ comp_flags) & (MD_FLAG_HETERO | MD_FLAG_WATER | MD_FLAG_ION | MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE)) break;
                    if (!test_name) {
                        entity_key = md_hash64_str(comp_name_j, entity_key);
                    }
                    ++j;
                }
            }

            entity_key = md_hash64(&entity_flags, sizeof(md_flags_t), entity_key);

            // See if this is a unique entity type or not
            md_entity_idx_t entity_idx = -1;
            for (size_t k = 0; k < md_array_size(entity_keys); ++k) {
                if (entity_keys[k] == entity_key) {
                    entity_idx = (md_entity_idx_t)k;
                    break;
                }
            }
            if (entity_idx == -1) {
                entity_idx = (md_entity_idx_t)sys->entity.count;
                md_label_t entity_id;
                entity_id.len = snprintf(entity_id.buf, sizeof(entity_id.buf), "%i", entity_idx + 1);

                const char* entity_type_str = "";
                if (entity_flags & MD_FLAG_POLYMER) {
                    if (entity_flags & MD_FLAG_AMINO_ACID) {
                        entity_type_str = "polypeptide";
                    } else if (entity_flags & MD_FLAG_NUCLEOTIDE) {
                        entity_type_str = "nucleic acid";
                    } else {
                        entity_type_str = "polymer";
                    }
                } else if (entity_flags & MD_FLAG_WATER) {
                    entity_type_str = "water";
                } else {
                    entity_type_str = comp_name.ptr;
                }

                char buf[256];
                snprintf(buf, sizeof(buf), "%s (*)", entity_type_str);

                // Create new entity
                md_array_push(sys->entity.id,    entity_id, alloc);
                md_array_push(sys->entity.flags, entity_flags, alloc);
                md_array_push(sys->entity.description, str_copy_cstr(buf, alloc), alloc);
                md_array_push(entity_keys, entity_key, temp_arena);
                sys->entity.count += 1;
#if 0
                MD_LOG_DEBUG("New entity: %s, %s", entity_id.buf, buf);
#endif
            }

            // Commit range (i,j) as an instance

            md_label_t inst_id = {0};
			// Create unique instance id
            {
				str_t last = sys->inst.count > 0 ? md_system_inst_id(sys, sys->inst.count - 1) : STR_LIT("");
                inst_id = md_util_next_inst_id(last);
                str_t id_str = LBL_TO_STR(inst_id);
				uint64_t key = md_hash64_str(id_str, 0);

                while (md_hashset_get(&inst_id_set, key)) {
                    inst_id = md_util_next_inst_id(id_str);
                    id_str = LBL_TO_STR(inst_id);
					key = md_hash64_str(id_str, 0);
                }
            }

            md_array_push(sys->inst.id, inst_id, alloc);
            md_array_push(sys->inst.auth_id, make_label(str_trim(comp_auth_id)), alloc);  // No auth id info as its generated
            md_array_push(sys->inst.comp_offset, (uint32_t)i, alloc);
            md_array_push(sys->inst.entity_idx, entity_idx, alloc);
            sys->inst.count += 1;

			md_hashset_add(&inst_id_set, md_hash64_str(LBL_TO_STR(inst_id), 0));
#if 0
            MD_LOG_DEBUG("New instance: %s (" STR_FMT "), %zu", inst_id.buf, STR_ARG(comp_auth_id), i);
#endif
            i = j;
        }
        md_array_push(sys->inst.comp_offset, (uint32_t)i, alloc);
    }

    md_vm_arena_destroy(temp_arena);
    return true;
}

bool md_util_set_ion_flags(md_flags_t* out_flags, const md_atom_data_t* atom, const md_bond_data_t* bond) {
    ASSERT(atom);
    ASSERT(bond);

    for (size_t i = 0; i < atom->count; ++i) {
        // Check if it has no bonds
        md_atomic_number_t z = md_atom_atomic_number(atom, i);

        if (md_bond_conn_count(bond, i) == 0 && monatomic_ion_element(z) && !(atom->flags[i] & MD_FLAG_WATER)) {
            atom->flags[i] |= MD_FLAG_ION;
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

static inline bool has_ring(const md_index_data_t* ring_data, const int* ring, int ring_size) {
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

bool md_util_system_infer_rings(md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);

    size_t num_atoms = md_system_atom_count(sys);
    size_t num_bonds = md_system_bond_count(sys);

    if (num_atoms == 0) {
        return false;
    }
    if (num_bonds == 0) {
        return false;
    }

    if (sys->ring.alloc && sys->ring.alloc != alloc) {
        md_index_data_free(&sys->ring);
    }
    sys->ring.alloc = alloc;
    md_index_data_clear(&sys->ring);

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    // Some typedefs to enable laboration of types
    typedef uint16_t color_t;
    typedef uint8_t  depth_t;
    typedef uint16_t mark_t;
    typedef  int32_t pred_t;
        
    color_t* color = md_vm_arena_push_zero_array(temp_arena, color_t, num_atoms);
    depth_t* depth = md_vm_arena_push_zero_array(temp_arena, depth_t, num_atoms);
    mark_t*  mark  = md_vm_arena_push_zero_array(temp_arena, mark_t,  num_atoms);
    pred_t*  pred  = md_vm_arena_push_array     (temp_arena, pred_t,  num_atoms);
    MEMSET(pred, -1, num_atoms * sizeof(pred_t));    // We can do memset as the representation of -1 under two's complement is 0xFFFFFFFF

    // The capacity is arbitrary here, but will be resized if needed.
    fifo_t queue = fifo_create(64, temp_arena);

    md_hashset_t ring_set = {.allocator = temp_arena};
    
    color_t current_color = 1;
    mark_t  current_mark  = 1;

#if DEBUG
    size_t processed_ring_elements = 0;
#endif

    size_t num_rings = 0;

    for (int atom_idx = 0; atom_idx < (int)num_atoms; ++atom_idx) {
        if (sys->atom.flags[atom_idx] & (MD_FLAG_WATER | MD_FLAG_ION)) continue;

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

            md_bond_iter_t it = md_bond_iter(&sys->bond, idx);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
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
                            md_index_data_push_arr(&sys->ring, ring, len);
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

    return true;
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE
#undef MAX_DEPTH

// Identifies isolated 'structures' defined by covalent bonds. Any set of atoms connected by covalent bonds are considered a structure
bool md_util_system_infer_structures(md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);

    if (sys->structure.alloc && sys->structure.alloc != alloc) {
        md_index_data_free(&sys->structure);
    }
    sys->structure.alloc = alloc;
    md_index_data_clear(&sys->structure);

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    const size_t atom_count = md_system_atom_count(sys);

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
            
            md_bond_iter_t it = md_bond_iter(&sys->bond, cur);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
                if (!bitfield_test_bit(visited, next)) {
                    fifo_push(&queue, next);
                }
                md_bond_iter_next(&it);
            }
        }

        // Sort the indices within the structure for more coherent memory access
        size_t size = md_array_size(indices);
        if (size < 128) {
            sort_arr(indices, (int)size);
        } else {
            uint32_t* temp = md_vm_arena_push(temp_arena, size * sizeof(uint32_t));
            sort_radix_inplace_uint32((uint32_t*)indices, size, temp);
        }
        
        // Here we should have exhausted every atom that is connected to index i.
        md_index_data_push_arr(&sys->structure, indices, md_array_size(indices));
        num_structures += 1;
    }
    
    md_vm_arena_destroy(temp_arena);

    return true;
}

typedef struct {
    const char* comp;
    const char* atom;
    int   atomic_nr;
    float mass;
    float radius;
    md_flags_t flags;
} atom_type_t;

// Estimated radius based on (C12/C6)^(1/6) parameters given for martini particle types
// P5 4.7
// P4 4.7
// P3 4.7
// P1 4.7
// C5 4.7          CYS
// C3 4.7          PRO
// C2 4.7          VAL
// C1              LEU
// Nda             ASN
// N0
// D
// Qa              ASP
// Qd
// SC5
// SC4
// SP1
// SNd
// SQd


// Predefined atom types (This include coarse grained types)
static const atom_type_t predefined_atom_types[] = {
    // Martini CG types (BB) + (SC*)
    {"ALA", "BB",  0,   89.09f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"ARG", "BB",  0,  174.20f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"ARG", "SC1", 0,  101.19f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"ARG", "SC2", 0,   70.09f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"ASN", "BB",  0,  132.12f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"ASN", "SC1", 0,   87.09f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"ASP", "BB",  0,  133.10f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"ASP", "SC1", 0,   96.06f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"CYS", "BB",  0,  121.16f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"CYS", "SC1", 0,  122.17f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"GLN", "BB",  0,  146.15f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"GLN", "SC1", 0,  111.14f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"GLU", "BB",  0,  147.13f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"GLU", "SC1", 0,  109.12f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"GLY", "BB",  0,   75.07f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"HIS", "BB",  0,  155.16f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"HIS", "SC1", 0,  110.14f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"HIS", "SC2", 0,   82.11f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"HIS", "SC3", 0,   40.04f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"ILE", "BB",  0,  131.18f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"ILE", "SC1", 0,  113.16f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"LEU", "BB",  0,  131.18f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"LEU", "SC1", 0,  113.16f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"LYS", "BB",  0,  146.19f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"LYS", "SC1", 0,  128.17f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"LYS", "SC2", 0,   84.11f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"MET", "BB",  0,  149.21f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"MET", "SC1", 0,  149.21f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"PHE", "BB",  0,  165.19f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"PHE", "SC1", 0,  135.18f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"PHE", "SC2", 0,   77.15f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"PHE", "SC3", 0,   39.04f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"PRO", "BB",  0,  115.13f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"SER", "BB",  0,  105.09f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"SER", "SC1", 0,   73.06f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"THR", "BB",  0,  119.12f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"THR", "SC1", 0,   87.09f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"TRP", "BB",  0,  204.23f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"TRP", "SC1", 0,  162.20f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"TRP", "SC2", 0,   77.15f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"TRP", "SC3", 0,   44.07f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"TRP", "SC4", 0,   15.04f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"TYR", "BB",  0,  181.19f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
    {"TYR", "SC1", 0,  136.17f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"TYR", "SC2", 0,   91.11f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
	{"TYR", "SC3", 0,   33.04f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},
    {"VAL", "BB",  0,  117.15f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_BACKBONE},
	{"VAL", "SC1", 0,   99.13f,    4.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_AMINO_ACID | MD_FLAG_SIDE_CHAIN},

	{"POPE", "PO4", 0,  95.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "GL1", 0,  55.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "GL2", 0,  55.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C1A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "D2A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C3A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C4A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C1B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C2B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C3B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
	{"POPE", "C4B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},

    {"POPG", "GL0", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "PO4", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "GL1", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "GL2", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C1A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "D2A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C3A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C4A", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C1B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C2B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C3B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},
    {"POPG", "C4B", 0,  72.00f,    4.1f, MD_FLAG_COARSE_GRAINED},

    {"RAMP", "PO1", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GM1", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GM2", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GM3", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GM4", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GM5", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GM6", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "PO2", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL1", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL2", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C1A", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C2A", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C3A", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C1B", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C2B", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C3B", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL3", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL4", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C1C", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C2C", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C3C", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C1D", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C2D", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C3D", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL5", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL6", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C1E", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C2E", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL7", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "GL8", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C1F", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "C2F", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S01", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S02", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S03", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S04", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S05", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S06", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S07", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S08", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S09", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S10", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S11", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S12", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S13", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S14", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S15", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S16", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S17", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S18", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S19", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S20", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S21", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S22", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S23", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S24", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S25", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S26", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S27", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S28", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S29", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S30", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S31", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S32", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S33", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S34", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S35", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S36", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S37", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S38", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},
    {"RAMP", "S39", 0, 50.f, 4.0f, MD_FLAG_COARSE_GRAINED},

	// Depending on the Martini version, water can be represented as a single bead or as 4-to-1 mapping
	{ "W", "W",     0,  18.015f,     2.3f, MD_FLAG_COARSE_GRAINED | MD_FLAG_WATER },

    { "ION", "CL",  17, 35.45f,     1.8f, MD_FLAG_ION },
    { "ION", "NA",  11, 22.99f,     2.3f, MD_FLAG_ION },
    { "ION", "K",   19, 39.10f,     2.7f, MD_FLAG_ION },
    { "ION", "CA",  20, 40.08f,     2.7f, MD_FLAG_ION },
	{ "ION", "MG",  12, 24.31f,     2.0f, MD_FLAG_ION },

    {"*", "IC", 0,   3240.0f,   19.5f, MD_FLAG_COARSE_GRAINED},
    {"*", "OC", 0,   3240.0f,   13.5f, MD_FLAG_COARSE_GRAINED},
    {"*", "CC", 0,   3900.0f,   10.0f, MD_FLAG_COARSE_GRAINED},

    {"P10", "X1", 0, 3.31E6f,   50.0f,  MD_FLAG_COARSE_GRAINED},
    {"P20", "X2", 0, 2.65E7f,   100.0f, MD_FLAG_COARSE_GRAINED},
    {"P30", "X3", 0, 8.93E7f,   150.0f, MD_FLAG_COARSE_GRAINED},
    {"P40", "X4", 0, 2.12E8f,   200.0f, MD_FLAG_COARSE_GRAINED},

    {"C##", "C#?", 0, 1727.7f,   13.5f, MD_FLAG_COARSE_GRAINED},
};

static inline uint64_t gen_key_from_names(str_t comp_name, str_t atom_name) {
    char buf[64];
    int len = snprintf(buf, sizeof(buf), STR_FMT "_" STR_FMT, STR_ARG(comp_name), STR_ARG(atom_name));
    str_t str = {buf, (size_t)len};
    return md_hash64_str(str, 0);
}

// Searches the predefined atom types for a matching atom name within a component
// We support simple pattern matching here like wildcards and numerical substitutes
// * is supported to match any sequence of characters
// ? is supported to match any single character
// # is supported to match any single digit

bool pattern_match(const char* pattern, const char* str) {
    while (*pattern && *str) {
        if (*pattern == '*') {
            pattern++;
            if (!*pattern) return true;
            while (*str) {
                if (pattern_match(pattern, str)) return true;
                str++;
            }
            return false;
        } else if (*pattern == '?') {
            pattern++;
            str++;
        } else if (*pattern == '#') {
            if (!is_digit((unsigned char)*str)) return false;
            pattern++;
            str++;
        } else {
            if (*pattern != *str) return false;
            pattern++;
            str++;
        }
    }
    while (*pattern == '*') pattern++;
    return !*pattern && !*str;
}

static atom_type_t* find_predefined_atom_type(str_t comp_name, str_t atom_name) {
    for (size_t i = 0; i < ARRAY_SIZE(predefined_atom_types); ++i) {
        atom_type_t type = predefined_atom_types[i];
        if (pattern_match(type.comp, comp_name.ptr) && pattern_match(type.atom, atom_name.ptr)) {
            return (atom_type_t*)&predefined_atom_types[i];
        }
    }
    return NULL;
}

void md_util_system_infer_atom_types(md_system_t* sys, const str_t atom_labels[], md_allocator_i* alloc) {
    if (!sys) {
        MD_LOG_ERROR("system was null");
        return;
    }
    if (!atom_labels) {
        MD_LOG_ERROR("Must provide atom_labels");
        return;
    }

    size_t temp_pos = md_temp_get_pos();

    md_hashmap32_t atom_type_cache = {.allocator = md_get_temp_allocator() };
    md_hashmap_reserve(&atom_type_cache, 256);
    
    size_t num_success = 0;

    if (sys->comp.count > 0) {
        for (size_t comp_idx = 0; comp_idx < sys->comp.count; ++comp_idx) {
            str_t comp_name        = md_comp_name(&sys->comp, comp_idx);
            md_urange_t comp_range = md_comp_atom_range(&sys->comp, comp_idx);
            size_t comp_size       = comp_range.end - comp_range.beg;
            uint64_t comp_key      = md_hash64_str(comp_name, comp_size);

            // Flags to propagate to the component from atom types
            md_flags_t comp_flags = 0;
            for (size_t i = comp_range.beg; i < comp_range.end; ++i) {
                if (sys->atom.type_idx[i] != 0) continue;

                uint64_t key = md_hash64_str(atom_labels[i], comp_key);
                uint32_t* cached_type = md_hashmap_get(&atom_type_cache, key);
                if (cached_type) {
                    sys->atom.type_idx[i] = *cached_type;
                } else {
                    str_t atom_name = atom_labels[i];
                    md_atomic_number_t z = 0;
                    float mass = 0;
                    float radius = 0;
                    uint32_t color = 0;
                    md_flags_t flags = 0;

                    // Try to find in predefined set
                    atom_type_t* predef_type = find_predefined_atom_type(comp_name, atom_name);
                    if (predef_type) {
                        z       = predef_type->atomic_nr;
                        mass    = predef_type->mass;
                        radius  = predef_type->radius;
                        color   = z == 0 ? 0 : md_atomic_number_cpk_color(z);
                        flags   = predef_type->flags;
                    } else {
                        z       = md_atomic_number_infer_from_label(atom_name, comp_name, comp_size);
                        mass    = md_atomic_number_mass(z);
                        radius  = md_atomic_number_vdw_radius(z);
                        color   = md_atomic_number_cpk_color(z);
                        flags   = 0;
                    }

                    comp_flags |= flags;

                    md_atom_type_idx_t type = md_atom_type_find_or_add(&sys->atom.type, atom_name, z, mass, radius, color, flags, alloc);
                    sys->atom.type_idx[i] = type;
                    md_hashmap_add(&atom_type_cache, key, (uint32_t)type);
                }
                sys->atom.flags[i] |= sys->atom.type.flags[sys->atom.type_idx[i]];
            }
            sys->comp.flags[comp_idx] |= comp_flags;
        }
    } else {
        for (size_t i = 0; i < sys->atom.count; ++i) {
            if (sys->atom.type_idx[i] != 0) continue;

            uint64_t key = md_hash64_str(atom_labels[i], 0);
            uint32_t* cached_type = md_hashmap_get(&atom_type_cache, key);
            if (cached_type) {
                sys->atom.type_idx[i] = *cached_type;
            } else {
                md_atomic_number_t z = md_atomic_number_infer_from_label(atom_labels[i], (str_t){0}, 0);
                float mass   = md_atomic_number_mass(z);
                float radius = md_atomic_number_vdw_radius(z);
                uint32_t color = md_atomic_number_cpk_color(z);
                md_flags_t flags = 0;
                md_atom_type_idx_t type = md_atom_type_find_or_add(&sys->atom.type, atom_labels[i], z, mass, radius, color, flags, alloc);
                sys->atom.type_idx[i] = type;
                md_hashmap_add(&atom_type_cache, key, (uint32_t)type);
            }
        }
    }

    md_temp_set_pos_back(temp_pos);
}

void md_util_mask_grow_by_bonds(md_bitfield_t* mask, const md_system_t* sys, size_t extent, const md_bitfield_t* viable_mask) {
    ASSERT(mask);
    ASSERT(sys);

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

    uint8_t* depth = md_alloc(temp_alloc, sys->atom.count * sizeof(uint8_t));
    MEMSET(depth, 0, sys->atom.count * sizeof(uint8_t));

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

            md_bond_iter_t it = md_bond_iter(&sys->bond, idx);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
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

void md_util_mask_grow_by_radius(md_bitfield_t* mask, const md_system_t* sys, double radius, const md_bitfield_t* viable_mask) {
    ASSERT(mask);
    ASSERT(sys);
    
    if (radius <= 0.0) return;
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    
    int32_t* indices = 0;
    size_t count = sys->atom.count;
        
    if (viable_mask) {
        md_bitfield_t tmp_bf = md_bitfield_create(arena);
        md_bitfield_andnot(&tmp_bf, viable_mask, mask);
            
        const size_t num_atoms = md_bitfield_popcount(&tmp_bf);
        indices = md_vm_arena_push(arena, num_atoms * sizeof(int32_t));
        count = md_bitfield_iter_extract_indices(indices, num_atoms, md_bitfield_iter_create(&tmp_bf));
    }

    md_spatial_hash_t* ctx = md_spatial_hash_create_soa(sys->atom.x, sys->atom.y, sys->atom.z, indices, count, &sys->unitcell, arena);

    md_bitfield_t old_mask = md_bitfield_create(arena);
    md_bitfield_copy(&old_mask, mask);

    md_bitfield_iter_t it = md_bitfield_iter_create(&old_mask);
    const float rad = (float)radius;
    while (md_bitfield_iter_next(&it)) {
        int idx = (int)md_bitfield_iter_idx(&it);
        const vec3_t pos = md_atom_coord(&sys->atom, idx);
        md_spatial_hash_query_bits(mask, ctx, pos, rad);
    }

    md_vm_arena_destroy(arena);
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


void md_util_oobb_compute(float out_basis[3][3], float out_ext_min[3], float out_ext_max[3], const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {
    ASSERT(out_basis);
    ASSERT(out_ext_min);
    ASSERT(out_ext_max);

    if (count == 0) return;

    vec3_t com = md_util_com_compute(in_x, in_y, in_z, NULL, in_idx, count, cell);
    double cov[3][3] = {0};
    if (cell) {
        vec4_t ref = vec4_from_vec3(com, 0);
        if (cell->flags & MD_UNITCELL_ORTHO) {
            vec4_t ext = md_unitcell_diag_vec4(cell);

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
        } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
            mat3_t basis = md_unitcell_basis_mat3(cell);
            for (size_t i = 1; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                vec4_t c = vec4_set(in_x[idx], in_y[idx], in_z[idx], 0.0f);
                deperiodize_triclinic(c.elem, ref.elem, basis.elem);

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
        int32_t idx = in_idx ? in_idx[i] : (int)i;
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

void md_util_oobb_compute_vec4(float out_basis[3][3], float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {
    ASSERT(out_basis);
    ASSERT(out_ext_min);
    ASSERT(out_ext_max);

    if (count == 0) return;

    vec3_t com = md_util_com_compute_vec4(in_xyzr, in_idx, count, cell);
    double cov[3][3] = {0};
    if (cell) {
        vec4_t ref = vec4_from_vec3(com, 0);
        if (cell->flags & MD_UNITCELL_ORTHO) {
            vec4_t ext = md_unitcell_diag_vec4(cell);

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
        } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
            mat3_t basis = md_unitcell_basis_mat3(cell);
            for (size_t i = 1; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                vec4_t c = in_xyzr[idx];
                deperiodize_triclinic(c.elem, ref.elem, basis.elem);

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
        int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
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
static void com_pbc(float* out_com, const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {

    mat3_t A  = md_unitcell_basis_mat3(cell);
    mat3_t Ai = md_unitcell_inv_basis_mat3(cell);

    // Should transform cartesian coordinates to fractional coordinates and scale by 2 pi
    mat3_t M = mat3_mul(mat3_scale(TWO_PI, TWO_PI, TWO_PI), Ai);
    // Inverse transform
    mat3_t I = mat3_mul(A, mat3_scale(1.0/TWO_PI, 1.0/TWO_PI, 1.0/TWO_PI));

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

static void com_pbc_vec4(float* out_com, const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {
    ASSERT(out_com);
    ASSERT(in_xyzw);
    ASSERT(cell);

    vec4_t acc_c = {0};
    vec4_t acc_s = {0};
    vec4_t acc_xyzw = {0};

    uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        const vec3_t ext = md_unitcell_diag_vec3(cell);
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
        ASSERT(flags & MD_UNITCELL_TRICLINIC);
        const mat3_t A  = md_unitcell_inv_basis_mat3(cell);
        const mat3_t Ai = md_unitcell_basis_mat3(cell);

        const mat4x3_t M = mat4x3_from_mat3(mat3_mul(mat3_scale(TWO_PI, TWO_PI, TWO_PI), Ai));
        const mat4x3_t I = mat4x3_from_mat3(mat3_mul(mat3_scale(1.0f / TWO_PI, 1.0f / TWO_PI, 1.0f / TWO_PI), A));

        if (in_idx) {
            for (size_t i = 0; i < count; ++i) {
                int32_t idx  = in_idx[i];
                vec4_t xyzw  = in_xyzw[idx];
                vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                vec4_t theta = mat4x3_mul_vec4(M, xyzw);
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
                vec4_t theta = mat4x3_mul_vec4(I, xyzw);
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

vec3_t md_util_com_compute(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* unit_cell) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    vec3_t xyz = {0};
    if (!unit_cell || unit_cell->flags == MD_UNITCELL_NONE) {
        com(xyz.elem, in_x, in_y, in_z, in_w, in_idx, count);
        return xyz;
    }

    if (unit_cell->flags & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) {
        com_pbc(xyz.elem, in_x, in_y, in_z, in_w, in_idx, count, unit_cell);
        return xyz;
    }

    // Error
    MD_LOG_ERROR("Invalid unit cell flags: %d", unit_cell->flags);
    return xyz;
}

vec3_t md_util_com_compute_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unitcell_t* unit_cell) {
    ASSERT(in_xyzw);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    vec3_t xyz = {0};
    if (unit_cell && (unit_cell->flags & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC))) {
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

void md_util_distance_array(float* out_dist, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    if (cell->flags == 0) {
        for (size_t i = 0; i < num_a; ++i) {
            for (size_t j = 0; j < num_b; ++j) {
                out_dist[i * num_b + j] = vec3_distance(coord_a[i], coord_b[j]);
            }
        }
    }
    else if (cell->flags & MD_UNITCELL_ORTHO) {
        const vec4_t box = md_unitcell_diag_vec4(cell);
        for (size_t i = 0; i < num_a; ++i) {
            for (size_t j = 0; j < num_b; ++j) {
                vec4_t a = vec4_from_vec3(coord_a[i], 0);
                vec4_t b = vec4_from_vec3(coord_b[j], 0);
                out_dist[i * num_b + j] = vec4_periodic_distance(a, b, box);
            }
        }
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        mat3_t A = md_unitcell_basis_mat3(cell);
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (size_t i = 0; i < num_a; ++i) {
            for (size_t j = 0; j < num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, A.elem);
                out_dist[i * num_b + j] = vec3_length(dx);
            }
        }
    }
}

float md_util_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
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
    else if (cell->flags & MD_UNITCELL_ORTHO) {
        const vec4_t box = md_unitcell_diag_vec4(cell);
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
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        const mat3_t A = md_unitcell_basis_mat3(cell);
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, A.elem);
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

float md_util_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
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
    else if (cell->flags & MD_UNITCELL_ORTHO) {
        const vec4_t box = md_unitcell_diag_vec4(cell);
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
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        const mat3_t A = md_unitcell_basis_mat3(cell);
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (int64_t i = 0; i < (int64_t)num_a; ++i) {
            for (int64_t j = 0; j < (int64_t)num_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, A.elem);
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

void md_util_min_image_vec3(vec3_t dx[], size_t count, const md_unitcell_t* cell) {
    if (cell) {
        vec3_t diag = md_unitcell_diag_vec3(cell);
        vec3_t half_diag = vec3_mul_f(diag, 0.5f);
        uint32_t flags = md_unitcell_flags(cell);
        if (flags & MD_UNITCELL_ORTHO) {
            for (size_t i = 0; i < count; ++i) {
                min_image_ortho(dx[i].elem, diag.elem, half_diag.elem);
            }
        } else if (flags & MD_UNITCELL_TRICLINIC) {
            mat3_t A = md_unitcell_basis_mat3(cell);
            for (size_t i = 0; i < count; ++i) {
                min_image_triclinic(dx[i].elem, A.elem, half_diag.elem);
            }
        }
    }
}

void md_util_min_image_vec4(vec4_t dx[], size_t count, const md_unitcell_t* cell) {
    if (cell) {
        vec3_t diag = md_unitcell_diag_vec3(cell);
        vec3_t half_diag = vec3_mul_f(diag, 0.5f);
        uint32_t flags = md_unitcell_flags(cell);
        if (flags & MD_UNITCELL_ORTHO) {
            for (size_t i = 0; i < count; ++i) {
                min_image_ortho(dx[i].elem, diag.elem, half_diag.elem);
            }
        } else if (flags & MD_UNITCELL_TRICLINIC) {
            mat3_t A = md_unitcell_basis_mat3(cell);
            for (size_t i = 0; i < count; ++i) {
                min_image_triclinic(dx[i].elem, A.elem, half_diag.elem);
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

static void pbc_triclinic(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    
    const vec4_t mask = md_unitcell_pbc_mask_vec4(cell);
    mat4x3_t I = mat4x3_from_mat3(md_unitcell_inv_basis_mat3(cell));
    mat4x3_t M = mat4x3_from_mat3(md_unitcell_basis_mat3(cell));

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

static void pbc_triclinic_vec4(vec4_t* xyzw, size_t count, const md_unitcell_t* cell) {
    const vec4_t mask = md_unitcell_pbc_mask_vec4(cell);
    mat4x3_t I = mat4x3_from_mat3(md_unitcell_inv_basis_mat3(cell));
    mat4x3_t M = mat4x3_from_mat3(md_unitcell_basis_mat3(cell));

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = xyzw[i];
        c = mat4x3_mul_vec4(I, c);
        c = vec4_fract(c);
        c = mat4x3_mul_vec4(M, c);
        xyzw[i] = vec4_blend(xyzw[i], c, mask);
    }
}

bool md_util_pbc(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing unitcell");
        return false;
    }

    uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        vec3_t ext = md_unitcell_diag_vec3(cell);
        pbc_ortho(x, y, z, indices, count, ext);
        return true;
    } else if (flags & MD_UNITCELL_TRICLINIC) {
        pbc_triclinic(x, y, z, indices, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    MD_LOG_ERROR("Unrecognized unit_cell type");
    return false;
}

bool md_util_pbc_vec4(vec4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell) {
    if (!in_out_xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing unit cell");
        return false;
    }

    uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        vec3_t ext = md_unitcell_diag_vec3(cell);
        pbc_ortho_vec4(in_out_xyzw, count, ext);
        return true;
    } else if (flags & MD_UNITCELL_TRICLINIC) {
        pbc_triclinic_vec4(in_out_xyzw, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    MD_LOG_ERROR("Unrecognized unit_cell type");
    return false;
}

bool md_util_system_pbc(md_system_t* sys) {
    ASSERT(sys);

	md_flags_t cell_flags = sys->unitcell.flags;

    if (md_unitcell_is_orthorhombic(&sys->unitcell)) {
        vec3_t ext = md_unitcell_diag_vec3(&sys->unitcell);
        pbc_ortho(sys->atom.x, sys->atom.y, sys->atom.z, 0, sys->atom.count, ext);
    } else if (md_unitcell_is_triclinic(&sys->unitcell)) {
		pbc_triclinic(sys->atom.x, sys->atom.y, sys->atom.z, 0, sys->atom.count, &sys->unitcell);
    }

    return true;
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

static void unwrap_triclinic(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    mat3_t A = md_unitcell_basis_mat3(cell);
    if (indices) {
        int idx = indices[0];
        vec3_t ref_pos = {x[idx], y[idx], z[idx]};
        for (size_t i = 1; i < count; ++i) {
            idx = indices[i];
            vec3_t pos = {x[idx], y[idx], z[idx]};
            deperiodize_triclinic(pos.elem, ref_pos.elem, A.elem);
            x[idx] = pos.x;
            y[idx] = pos.y;
            z[idx] = pos.z;
            ref_pos = pos;
        }
    } else {
        vec3_t ref_pos = {x[0], y[0], z[0]};
        for (size_t i = 1; i < count; ++i) {
            vec3_t pos = {x[i], y[i], z[i]};
            deperiodize_triclinic(pos.elem, ref_pos.elem, A.elem);
            x[i] = pos.x;
            y[i] = pos.y;
            z[i] = pos.z;
            ref_pos = pos;
        }
    }
}

static void unwrap_triclinic_vec4(vec4_t* xyzw, size_t count, const md_unitcell_t* cell) {
    mat3_t A = md_unitcell_basis_mat3(cell);
    vec4_t ref_pos = xyzw[0];
    for (size_t i = 1; i < count; ++i) {
        vec4_t pos = xyzw[i];
        deperiodize_triclinic(pos.elem, ref_pos.elem, A.elem);
        xyzw[i] = pos;
        ref_pos = pos;
    }

}

bool md_util_unwrap(float* x, float* y, float* z, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    if (!x || !y || !z) {
        MD_LOG_ERROR("Missing required input: x,y or z");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing cell");
        return false;
    }

    if (cell->flags & MD_UNITCELL_ORTHO) {
        vec3_t ext = md_unitcell_diag_vec3(cell);
        unwrap_ortho(x, y, z, indices, count, ext);
        return true;
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        unwrap_triclinic(x, y, z, indices, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    return false;
}

bool md_util_unwrap_vec4(vec4_t* xyzw, size_t count, const md_unitcell_t* cell) {
    if (!xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing cell");
        return false;
    }

    if (md_unitcell_is_orthorhombic(cell)) {
        vec3_t ext = md_unitcell_diag_vec3(cell);
        unwrap_ortho_vec4(xyzw, count, ext);
        return true;
    } else if (md_unitcell_is_triclinic(cell)) {
        unwrap_triclinic_vec4(xyzw, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    return false;
}

bool md_util_system_unwrap(md_system_t* sys) {
    ASSERT(sys);
    size_t num_structures = md_index_data_num_ranges(&sys->structure);
	float* x = sys->atom.x;
	float* y = sys->atom.y;
	float* z = sys->atom.z;

    md_flags_t cell_flags = sys->unitcell.flags;
    if (cell_flags & MD_UNITCELL_ORTHO) {
        vec3_t ext = md_unitcell_diag_vec3(&sys->unitcell);
        for (size_t i = 0; i < num_structures; ++i) {
            const int32_t* s_idx = md_index_range_beg(&sys->structure, i);
            const size_t   s_len = md_index_range_size(&sys->structure, i);
            unwrap_ortho(x, y, z, s_idx, s_len, ext);
        }
    } else if (cell_flags & MD_UNITCELL_TRICLINIC) {
        for (size_t i = 0; i < num_structures; ++i) {
            const int32_t* s_idx = md_index_range_beg(&sys->structure, i);
            const size_t   s_len = md_index_range_size(&sys->structure, i);
            unwrap_triclinic(x, y, z, s_idx, s_len, &sys->unitcell);
        }
    }
    
    return true;
}

bool md_util_deperiodize_vec4(vec4_t* xyzw, size_t count, vec3_t ref_xyz, const md_unitcell_t* cell) {
    if (!xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing cell");
        return false;
    }

    if (cell->flags & MD_UNITCELL_ORTHO) {
        vec4_t ref_pos = vec4_from_vec3(ref_xyz, 0);
        vec4_t ext = md_unitcell_diag_vec4(cell);
        for (size_t i = 0; i < count; ++i) {
            const vec4_t pos = vec4_deperiodize_ortho(xyzw[i], ref_pos, ext);
            xyzw[i] = pos;
        }
        return true;
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        mat3_t A = md_unitcell_basis_mat3(cell);
        for (size_t i = 1; i < count; ++i) {
            vec4_t pos = xyzw[i];
            deperiodize_triclinic(pos.elem, ref_xyz.elem, A.elem);
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

bool md_util_interpolate_linear(float* out_x, float* out_y, float* out_z, const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], size_t count, const md_unitcell_t* cell, float t) {
    ASSERT(out_x);
    ASSERT(out_y);
    ASSERT(out_z);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    t = CLAMP(t, 0.0f, 1.0f);

    bool ortho = md_unitcell_is_orthorhombic(cell);
    bool tricl = md_unitcell_is_triclinic(cell);
    
    if (ortho | tricl) {
        mat3_t A = md_unitcell_basis_mat3(cell);
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

            if (ortho) {
                simd_deperiodize_ortho(x1[0], x0[0], md_mm256_set1_ps(A.elem[0][0]), md_mm256_set1_ps(A.elem[0][0]));
                simd_deperiodize_ortho(x1[1], x0[1], md_mm256_set1_ps(A.elem[1][1]), md_mm256_set1_ps(A.elem[1][1]));
                simd_deperiodize_ortho(x1[2], x0[2], md_mm256_set1_ps(A.elem[2][2]), md_mm256_set1_ps(A.elem[2][2]));
            } else if (tricl) {
                simd_deperiodize_triclinic(x1, x0, A.elem);
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

bool md_util_interpolate_cubic_spline(float* out_x, float* out_y, float* out_z, const float* const in_x[4], const float* const in_y[4], const float* const in_z[4], size_t count, const md_unitcell_t* cell, float t, float s) {
    ASSERT(out_x);
    ASSERT(out_y);
    ASSERT(out_z);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    
    t = CLAMP(t, 0.0f, 1.0f);
    s = CLAMP(s, 0.0f, 1.0f);

    bool ortho = md_unitcell_is_orthorhombic(cell);
    bool tricl = md_unitcell_is_triclinic(cell);
    mat3_t A   = md_unitcell_basis_mat3(cell);
    mat3_t Ai  = md_unitcell_inv_basis_mat3(cell);

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

        if (ortho) {
            x0[0] = simd_deperiodize_ortho(x0[0], x1[0], md_mm256_set1_ps(A.elem[0][0]), md_mm256_set1_ps(Ai.elem[0][0]));
            x2[0] = simd_deperiodize_ortho(x2[0], x1[0], md_mm256_set1_ps(A.elem[0][0]), md_mm256_set1_ps(Ai.elem[0][0]));
            x3[0] = simd_deperiodize_ortho(x3[0], x2[0], md_mm256_set1_ps(A.elem[0][0]), md_mm256_set1_ps(Ai.elem[0][0]));
            x0[1] = simd_deperiodize_ortho(x0[1], x1[1], md_mm256_set1_ps(A.elem[1][1]), md_mm256_set1_ps(Ai.elem[1][1]));
            x2[1] = simd_deperiodize_ortho(x2[1], x1[1], md_mm256_set1_ps(A.elem[1][1]), md_mm256_set1_ps(Ai.elem[1][1]));
            x3[1] = simd_deperiodize_ortho(x3[1], x2[1], md_mm256_set1_ps(A.elem[1][1]), md_mm256_set1_ps(Ai.elem[1][1]));
            x0[2] = simd_deperiodize_ortho(x0[2], x1[2], md_mm256_set1_ps(A.elem[2][2]), md_mm256_set1_ps(Ai.elem[2][2]));
            x2[2] = simd_deperiodize_ortho(x2[2], x1[2], md_mm256_set1_ps(A.elem[2][2]), md_mm256_set1_ps(Ai.elem[2][2]));
            x3[2] = simd_deperiodize_ortho(x3[2], x2[2], md_mm256_set1_ps(A.elem[2][2]), md_mm256_set1_ps(Ai.elem[2][2]));
        } else if (tricl) {
            simd_deperiodize_triclinic(x0, x1, A.elem);
            simd_deperiodize_triclinic(x2, x1, A.elem);
            simd_deperiodize_triclinic(x3, x2, A.elem);
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

static inline void commit_protein_backbone(md_protein_backbone_atoms_t* bb_atoms, size_t bb_length, size_t comp_base_idx, size_t inst_idx, md_protein_backbone_data_t* bb_data, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(bb_data->segment.atoms);
    md_array_push_array(bb_data->segment.atoms, bb_atoms, bb_length, alloc);
    md_array_push(bb_data->range.offset, offset, alloc);
    md_array_push(bb_data->range.inst_idx, (md_inst_idx_t)inst_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        int32_t comp_idx = (int32_t)(comp_base_idx + i);
        md_array_push(bb_data->segment.comp_idx, comp_idx, alloc);
    }
}

static inline void flush_protein_backbone_segment(md_protein_backbone_atoms_t* backbone_atoms, md_comp_idx_t* comp_base, size_t inst_idx, md_protein_backbone_data_t* bb_data, md_allocator_i* alloc, size_t min_backbone_length) {
    size_t backbone_length = md_array_size(backbone_atoms);
    if (backbone_length >= min_backbone_length && *comp_base != -1) {
        commit_protein_backbone(backbone_atoms, backbone_length, *comp_base, inst_idx, bb_data, alloc);
    }
    md_array_shrink(backbone_atoms, 0);
    *comp_base = -1;
}

static inline void commit_nucleic_backbone(md_nucleic_backbone_atoms_t* bb_atoms, size_t bb_length, size_t comp_base_idx, size_t inst_idx, md_nucleic_backbone_data_t* bb_data, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(bb_data->segment.atoms);
    md_array_push_array(bb_data->segment.atoms, bb_atoms, bb_length, alloc);
    md_array_push(bb_data->range.offset, offset, alloc);
    md_array_push(bb_data->range.inst_idx, (md_inst_idx_t)inst_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        int32_t comp_idx = (int32_t)(comp_base_idx + i);
        md_array_push(bb_data->segment.comp_idx, comp_idx, alloc);
    }
}

static inline void flush_nucleic_backbone_segment(md_nucleic_backbone_atoms_t* backbone_atoms, md_comp_idx_t* comp_base, size_t inst_idx, md_nucleic_backbone_data_t* bb_data, md_allocator_i* alloc, size_t min_backbone_length) {
    size_t backbone_length = md_array_size(backbone_atoms);
    if (backbone_length >= min_backbone_length && *comp_base != -1) {
        commit_nucleic_backbone(backbone_atoms, backbone_length, *comp_base, inst_idx, bb_data, alloc);
    }
    md_array_shrink(backbone_atoms, 0);
    *comp_base = -1;
}

static inline vec3_t hcl_to_rgb(float h, float c, float l) {
    float r = 0, g = 0, b = 0;

    if (c == 0) {
        r = g = b = l;
    } else {
        float h_prime = h * 6.0f;
        float x = c * (1.0f - fabsf(fmodf(h_prime, 2.0f) - 1.0f));
        float m = l - c * 0.5f;

        if (0.0f <= h_prime && h_prime < 1.0f) {
            r = c; g = x; b = 0;
        } else if (1.0f <= h_prime && h_prime < 2.0f) {
            r = x; g = c; b = 0;
        } else if (2.0f <= h_prime && h_prime < 3.0f) {
            r = 0; g = c; b = x;
        } else if (3.0f <= h_prime && h_prime < 4.0f) {
            r = 0; g = x; b = c;
        } else if (4.0f <= h_prime && h_prime < 5.0f) {
            r = x; g = 0; b = c;
        } else if (5.0f <= h_prime && h_prime < 6.0f) {
            r = c; g = 0; b = x;
        }

        r += m;
        g += m;
        b += m;
    }

    return (vec3_t){r, g, b};
}

// Try to fill in missing fields for molecule struct
// (Coordinates & Elements) -> Covalent Bonds
// (residues & Bonds)       -> Chains
// (Chains)                 -> Backbone
bool md_util_molecule_postprocess(md_system_t* sys, md_allocator_i* alloc, md_util_postprocess_flags_t flags) {
    ASSERT(sys);
    ASSERT(alloc);

    if (sys->atom.count == 0) return false;

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    bool cg = false;
    for (size_t i = 0; i < sys->atom.type.count; ++i) {
        if (sys->atom.type.flags[i] & MD_FLAG_COARSE_GRAINED) {
            cg = true;
            break;
        }
    }

    if (cg) {
        flags &= ~MD_UTIL_POSTPROCESS_BOND_BIT;
        flags &= ~MD_UTIL_POSTPROCESS_HBOND_BIT;
    }

    if (flags & MD_UTIL_POSTPROCESS_COLOR_BIT) {
        if (sys->atom.type.color) {
            for (size_t i = 1; i < sys->atom.type.count; ++i) {
                uint32_t color = sys->atom.type.color[i];
                md_atomic_number_t z = sys->atom.type.z[i];
                if (color == 0) {
                    // Try to assign a color based on element if available
                    if (z != 0) {
                        color = md_util_element_cpk_color(z);
                        sys->atom.type.color[i] = color;
                    } else {
                        // Assign a random color for the type
                        float hue = (float)(i * 97 % 360) / 360.0f;
                        float chroma = 0.8f;
                        float light  = 1.0f;
                        float alpha  = 1.0f;
                        vec3_t rgb = hcl_to_rgb(hue, chroma, light);
                        // Pack into 0xAARRGGBB
                        color = (uint32_t)((uint8_t)(alpha * 255) << 24 | (uint8_t)(rgb.x * 255) << 16 | (uint8_t)(rgb.y * 255) << 8 | (uint8_t)(rgb.z * 255));
                        sys->atom.type.color[i] = color;
                    }
                }
            }
        }
    }

    if (!sys->atom.flags) {
        md_array_resize(sys->atom.flags, sys->atom.count, alloc);
        MEMSET(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    }
   
    if (flags & MD_UTIL_POSTPROCESS_BOND_BIT) {
        if (sys->bond.count == 0) {
            md_util_system_infer_covalent_bonds(sys, alloc);
        }
        if (sys->bond.conn.count == 0) {
			md_bond_build_connectivity(&sys->bond, sys->atom.count, alloc);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_STRUCTURE_BIT) {
        if (sys->bond.count) {
            md_util_system_infer_structures(sys, alloc);
            md_util_system_infer_rings(sys, alloc);
        }
    }
    
#if 0
    if (flags & MD_UTIL_POSTPROCESS_ORDER_BIT) {
        compute_covalent_bond_order(&sys->bond, &sys->atom, &sys->ring);
    }
#endif

    if (flags & MD_UTIL_POSTPROCESS_INSTANCE_BIT) {
        if (sys->inst.count == 0 && sys->comp.count > 0 && sys->bond.pairs) {
            md_util_system_infer_entity_and_instance(sys, NULL, alloc);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_HBOND_BIT) {
        if (sys->atom.count > 0 && sys->bond.count > 0) {
            md_util_hydrogen_bond_init (&sys->hydrogen_bond, sys, alloc);
            md_util_hydrogen_bond_infer(&sys->hydrogen_bond, sys->atom.x, sys->atom.y, sys->atom.z, &sys->unitcell, 3.0, 150.0);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_BACKBONE_BIT) {
        if (sys->inst.count && sys->atom.type_idx) {
            // Compute backbone data
            // 
            // @NOTE: We should only attempt to compute backbone data for valid components (e.g. amino acids / dna / rna)
            // Backbones are not directly tied to instances, but a polymer instance has the potential to have a backbone.
            // We look within these instances and see if we can find consecutive ranges which form backbones.

            static const size_t MIN_BACKBONE_LENGTH = 3;
            static const size_t MAX_BACKBONE_LENGTH = 1 << 15;

            size_t temp_pos = md_vm_arena_get_pos(temp_arena);

            {
                md_flags_t req_flags = MD_FLAG_POLYMER | MD_FLAG_AMINO_ACID;
                md_array(md_protein_backbone_atoms_t) backbone_atoms = 0;
                md_array_ensure(backbone_atoms, MAX_BACKBONE_LENGTH, temp_arena);
                md_comp_idx_t comp_base = -1;
                md_seq_id_t prev_seq_id = -1;

                for (size_t inst_idx = 0; inst_idx < sys->inst.count; ++inst_idx) {
                    md_flags_t inst_flags = md_system_inst_flags(sys, inst_idx);
                    
                    // Check for polymer and amino acid otherwise skip
                    if ((inst_flags & req_flags) != req_flags) {
                        continue;
                    }

                    md_urange_t range = md_inst_comp_range(&sys->inst, inst_idx);
                    for (md_comp_idx_t res_idx = range.beg; res_idx < range.end; ++res_idx) {
                        md_protein_backbone_atoms_t atoms;
                        md_urange_t atom_range = md_comp_atom_range(&sys->comp, res_idx);
                        bool has_backbone_atoms = md_util_protein_backbone_atoms_extract(&atoms, &sys->atom, atom_range);

                        if (has_backbone_atoms) {
                            md_seq_id_t seq_id = md_system_comp_seq_id(sys, res_idx);
                            size_t backbone_length = md_array_size(backbone_atoms);

                            if (backbone_length > 0 && seq_id != prev_seq_id + 1) {
                                flush_protein_backbone_segment(backbone_atoms, &comp_base, inst_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                            }

                            if (comp_base == -1) {
                                comp_base = res_idx;
                            }
                            md_array_push(backbone_atoms, atoms, temp_arena);
                            prev_seq_id = seq_id;
                            continue;
                        }

                        flush_protein_backbone_segment(backbone_atoms, &comp_base, inst_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                        prev_seq_id = -1;
                    }
                    // Possibly commit remainder of the chain
                    flush_protein_backbone_segment(backbone_atoms, &comp_base, inst_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                    prev_seq_id = -1;
                }

                sys->protein_backbone.range.count = md_array_size(sys->protein_backbone.range.offset);
                if (sys->protein_backbone.range.count) {
                    // Add end offset
                    md_array_push(sys->protein_backbone.range.offset, (uint32_t)md_array_size(sys->protein_backbone.segment.atoms), alloc);
                }
                sys->protein_backbone.segment.count = md_array_size(sys->protein_backbone.segment.atoms);

                md_array_resize(sys->protein_backbone.segment.angle, sys->protein_backbone.segment.count, alloc);
                md_array_resize(sys->protein_backbone.segment.secondary_structure, sys->protein_backbone.segment.count, alloc);
                md_array_resize(sys->protein_backbone.segment.rama_type, sys->protein_backbone.segment.count, alloc);

                if (sys->protein_backbone.segment.count > 0) {
                    md_util_backbone_angles_compute(sys->protein_backbone.segment.angle, sys->protein_backbone.segment.count, sys->atom.x, sys->atom.y, sys->atom.z, &sys->unitcell, &sys->protein_backbone);
                    md_util_backbone_secondary_structure_infer(sys->protein_backbone.segment.secondary_structure, sys->protein_backbone.segment.count, sys->atom.x, sys->atom.y, sys->atom.z, &sys->unitcell, &sys->protein_backbone);
                    md_util_backbone_ramachandran_classify(sys->protein_backbone.segment.rama_type, sys->protein_backbone.segment.count, sys);
                }
            }
            md_vm_arena_set_pos_back(temp_arena, temp_pos);

            {
                md_flags_t req_flags = MD_FLAG_POLYMER | MD_FLAG_NUCLEOTIDE;
                md_array(md_nucleic_backbone_atoms_t) backbone_atoms = 0;
                md_array_ensure(backbone_atoms, MAX_BACKBONE_LENGTH, temp_arena);
                md_comp_idx_t comp_base = -1;
                md_seq_id_t prev_seq_id = -1;

                for (size_t inst_idx = 0; inst_idx < sys->inst.count; ++inst_idx) {
                    md_flags_t inst_flags = md_system_inst_flags(sys, inst_idx);

                    // Check for polymer and nucleotide otherwise skip
                    if ((inst_flags & req_flags) != req_flags) {
                        continue;
                    }

                    md_urange_t range = md_inst_comp_range(&sys->inst, inst_idx);
                    for (md_comp_idx_t comp_idx = range.beg; comp_idx < range.end; ++comp_idx) {
                        md_nucleic_backbone_atoms_t atoms;
                        md_urange_t atom_range = md_comp_atom_range(&sys->comp, comp_idx);
                        bool has_backbone_atoms = md_util_nucleic_backbone_atoms_extract(&atoms, &sys->atom, atom_range);

                        if (has_backbone_atoms) {
                            md_seq_id_t seq_id = md_system_comp_seq_id(sys, comp_idx);
                            size_t backbone_length = md_array_size(backbone_atoms);

                            if (backbone_length > 0 && seq_id != prev_seq_id + 1) {
                                flush_nucleic_backbone_segment(backbone_atoms, &comp_base, inst_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                                backbone_length = 0;
                            }

                            bool can_extend_segment = (backbone_length == 0) || (atoms.p != -1);

                            if (can_extend_segment) {
                                if (comp_base == -1) {
                                    comp_base = comp_idx;
                                }
                                md_array_push(backbone_atoms, atoms, temp_arena);
                                prev_seq_id = seq_id;
                                continue;
                            }
                        }

                        flush_nucleic_backbone_segment(backbone_atoms, &comp_base, inst_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                        prev_seq_id = -1;
                    }
                    // Possibly commit remainder of the chain
                    flush_nucleic_backbone_segment(backbone_atoms, &comp_base, inst_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                    prev_seq_id = -1;
                }

                sys->nucleic_backbone.range.count = md_array_size(sys->protein_backbone.range.offset);
                if (sys->nucleic_backbone.range.count) {
                    // Add end offset
                    md_array_push(sys->nucleic_backbone.range.offset, (uint32_t)md_array_size(sys->protein_backbone.segment.atoms), alloc);
                }
                sys->nucleic_backbone.segment.count = md_array_size(sys->nucleic_backbone.segment.atoms);
            }
            md_vm_arena_set_pos_back(temp_arena, temp_pos);
        }
    }

    if ((flags & MD_UTIL_POSTPROCESS_UNWRAP_STRUCTURE_BIT) && 
        (sys->unitcell.flags != 0) && md_structure_count(&sys->structure) > 0)
    {
        md_util_system_unwrap(sys);
    }

	md_vm_arena_destroy(temp_arena);
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

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    uint32_t* temp = md_vm_arena_push(arena, count * sizeof(uint32_t));
    sort_radix_inplace_uint32(data, count, temp);
    md_vm_arena_destroy(arena);
}

void md_util_sort_radix_uint32(uint32_t* out_indices, const uint32_t* in_keys, size_t count, uint32_t* tmp_indices) {
    ASSERT(in_keys);
    ASSERT(out_indices);
    ASSERT(tmp_indices);
    sort_radix_uint32(out_indices, in_keys, count, tmp_indices);
}

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
            int* ptr = md_index_range_beg(data->result, entry->index);
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

// Create a new reference structure which is pruned of certain atoms (Hydrogen) and loosely connected subcomponents
// There are simply too many permutations to cover and the result will explode.
static md_array(int) filter_structure_connectivity(const int* indices, size_t count, const md_index_data_t* connectivity, int min_val, md_allocator_i* alloc) {
    md_array(int) filt = 0;
    for (size_t i = 0; i < count; ++i) {
        int idx = indices[i];
        int order = (int)md_index_range_size(connectivity, idx);
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

static inline bool filter_atom(const md_system_t* sys, int atom_i, md_util_match_flags_t filter) {
    md_atomic_number_t z_i = md_atom_atomic_number(&sys->atom, atom_i);
    if (filter & MD_UTIL_MATCH_FLAGS_NO_H) {
        return z_i != MD_Z_H;
    } else if (filter & MD_UTIL_MATCH_FLAGS_NO_CH) {
        if (z_i == MD_Z_H) {
            uint32_t   conn_len = (uint32_t)md_bond_conn_count(&sys->bond, atom_i);
            uint32_t   conn_i   = sys->bond.conn.offset[atom_i];
            for (uint32_t i = 0; i < conn_len; ++i) {
                int atom_j = md_bond_conn_atom_idx(&sys->bond, conn_i, i);
                if (md_atom_atomic_number(&sys->atom, atom_j) == MD_Z_C) {
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

static size_t extract_structures(md_index_data_t* out_structures, const md_system_t* sys, md_util_match_level_t level, md_util_match_flags_t filter, size_t min_size) {
    ASSERT(out_structures);
    size_t pre_offset = md_index_data_num_ranges(out_structures);

    switch (level) {
    case MD_UTIL_MATCH_LEVEL_STRUCTURE:
        if (filter) {
            for (size_t s_idx = 0; s_idx < md_index_data_num_ranges(&sys->structure); ++s_idx) {
                if (md_index_range_size(&sys->structure, s_idx) < min_size) {
                    continue;
                }
                idx_range_t range = idx_range_create(out_structures);
                for (int* it = md_index_range_beg(&sys->structure, s_idx); it < md_index_range_end(&sys->structure, s_idx); ++it) {
                    int i = *it;
                    if (filter_atom(sys, i, filter)) {
                        idx_range_push(range, i);
                    }
                }
                idx_range_commit(range);
            }
        } else {
            for (size_t s_idx = 0; s_idx < md_index_data_num_ranges(&sys->structure); ++s_idx) {
                if (md_index_range_size(&sys->structure, s_idx) < min_size) {
                    continue;
                }
                md_index_data_push_arr(out_structures, md_index_range_beg(&sys->structure, s_idx), md_index_range_size(&sys->structure, s_idx));
            }
        }
        break;
    case MD_UTIL_MATCH_LEVEL_COMPONENT:
        for (size_t r_idx = 0; r_idx < sys->comp.count; ++r_idx) {
            if (md_comp_atom_count(&sys->comp, r_idx) < min_size) {
                continue;
            }
            idx_range_t range = idx_range_create(out_structures);
            md_urange_t comp_range = md_comp_atom_range(&sys->comp, r_idx);
            for (int i = comp_range.beg; i < comp_range.end; ++i) {
                if (filter_atom(sys, i, filter)) {
                    idx_range_push(range, i);
                }
            }
            idx_range_commit(range);
        }
        break;
    case MD_UTIL_MATCH_LEVEL_INSTANCE:
        for (size_t i_idx = 0; i_idx < sys->inst.count; ++i_idx) {
            md_urange_t inst_range = md_system_inst_atom_range(sys, i_idx);
            if (inst_range.end - inst_range.beg < min_size) {
                continue;
            }
            idx_range_t range = idx_range_create(out_structures);
            for (int i = inst_range.beg; i < inst_range.end; ++i) {
                if (filter_atom(sys, i, filter)) {
                    idx_range_push(range, i);
                }
            }
            idx_range_commit(range);
        }
        break;
    default:
        ASSERT(false);
    }

    size_t post_offset = md_index_data_num_ranges(out_structures);

    return post_offset - pre_offset;
}

#define MAX_TYPES 256
md_index_data_t match_structure(const int* ref_idx, size_t ref_len, md_util_match_mode_t mode, md_util_match_level_t level, vertex_type_mapping_mode_t vertex_mapping, const md_system_t* sys, md_allocator_i* alloc) {
    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    graph_t ref_graph = {0};
    md_index_data_t result = {.alloc = alloc};

    size_t atom_count = md_system_atom_count(sys);
    md_array(int) structure_idx = md_vm_arena_push_array(temp_arena, int, atom_count);

    for (size_t i = 0; i < md_structure_count(&sys->structure); ++i) {
        int* beg = md_index_range_beg(&sys->structure, i);
        int* end = md_index_range_end(&sys->structure, i);
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
        if (md_atom_atomic_number(&sys->atom, idx) == MD_Z_H) {
            ref_hydro_present = true;
            break;
        }
    }

    const md_util_match_flags_t flags = ref_hydro_present ? 0 : MD_UTIL_MATCH_FLAGS_NO_H;

    ref_graph = extract_graph(sys, ref_idx, ref_len, vertex_mapping, temp_arena);

    int ref_type_count[MAX_TYPES] = {0};
    for (size_t i = 0; i < ref_len; ++i) {
        uint8_t type = (uint8_t)graph_vertex_type(&ref_graph, i);
        ref_type_count[type]++;
    }

    md_index_data_t structures = { .alloc = temp_arena };
    const size_t num_structures = extract_structures(&structures, sys, level, flags, ref_len);

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

        const int*   s_idx = md_index_range_beg(&structures, i);
        const size_t s_len = md_index_range_size(&structures, i);

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
            uint8_t type = vertex_type_from_atom_idx(&sys->atom, idx, vertex_mapping);
            s_type_count[type]++;
        }

        // Sanity check
        // The structure needs to have atleast the same amount of types as the reference
        for (int j = 0; j < MAX_TYPES; ++j) {
            if (ref_type_count[j] > s_type_count[j]) {
                goto next;
            }
        }

        graph_t graph = extract_graph(sys, s_idx, s_len, vertex_mapping, temp_arena);
        size_t pre_count = md_index_data_num_ranges(&result);

        state_t state = {0};
        state_init(&state, &ref_graph, &graph, temp_arena);
        state.callback = cb;
        state.user_data = &data;

        MEMSET(&data.map, 0, sizeof(result_map_t));
        data.map.allocator = temp_arena;

        find_isomorphisms_callback(&ref_graph, &graph, starting_type, &state);
        size_t post_count = md_index_data_num_ranges(&result);

        // Remap indices to global indices in result
        for (size_t j = pre_count; j < post_count; ++j) {
            int* beg = md_index_range_beg(&result, j);
            int* end = md_index_range_end(&result, j);
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

md_index_data_t md_util_match_by_type(const int ref_indices[], size_t ref_count, md_util_match_mode_t mode, md_util_match_level_t level, const md_system_t* sys, md_allocator_i* alloc) {
    return match_structure(ref_indices, ref_count, mode, level, VERTEX_TYPE_MAPPING_MODE_ATOM_TYPE_IDX, sys, alloc);
}

md_index_data_t md_util_match_by_element(const int ref_indices[], size_t ref_count, md_util_match_mode_t mode, md_util_match_level_t level, const md_system_t* sys, md_allocator_i* alloc) {
    return match_structure(ref_indices, ref_count, mode, level, VERTEX_TYPE_MAPPING_MODE_ATOMIC_NUMBER, sys, alloc);
}

size_t md_util_match_smiles(md_index_data_t* idx_data, str_t smiles, md_util_match_mode_t mode, md_util_match_level_t level, md_util_match_flags_t flags, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(sys);

    if (idx_data && !idx_data->alloc) {
        MD_LOG_ERROR("Incomplete idx_data structure supplied");
        return 0;
    }

    if (!sys->atom.type_idx) {
        MD_LOG_ERROR("Missing required atom type field");
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
    const size_t num_structures = extract_structures(&structures, sys, level, flags, ref_graph.vertex_count);

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

    size_t atom_count = md_system_atom_count(sys);
    uint8_t* atom_types = md_vm_arena_push(temp_alloc, atom_count);
    md_atom_extract_atomic_numbers(atom_types, 0, atom_count, &sys->atom);

    for (size_t i = 0; i < num_structures; ++i) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_alloc);
        const size_t s_size = md_index_range_size(&structures, i);
        const int*   s_idx  = md_index_range_beg(&structures, i);

        int s_type_count[256] = {0};
        for (size_t j = 0; j < s_size; ++j) {
            int idx = s_idx[j];
            md_atomic_number_t type = md_atom_atomic_number(&sys->atom, idx);
            s_type_count[type]++;
        }

        // Sanity check
        for (size_t j = 0; j < Num_Elements; ++j) {
            if (ref_type_count[j] > s_type_count[j]) goto next;
        }

        graph_t s_graph = extract_graph(sys, s_idx, s_size, VERTEX_TYPE_MAPPING_MODE_ATOMIC_NUMBER, temp_alloc);

        if (flags & MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT) {
            if (s_graph.vertex_count != ref_graph.vertex_count) {
                continue;
            }
        }
        
        size_t pre_count = md_index_data_num_ranges(idx_data);
        size_t count = find_isomorphisms(idx_data, &ref_graph, &s_graph, mode, starting_type, alloc, temp_alloc);
        size_t post_count = md_index_data_num_ranges(idx_data);

        // Remap indices to global indices in result
        if (idx_data && count > 0) {
            ASSERT(post_count - pre_count == count);
            for (size_t j = pre_count; j < post_count; ++j) {
                int* beg = md_index_range_beg(idx_data, j);
                int* end = md_index_range_end(idx_data, j);
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
