#include <md_util.h>

#include <md_system.h>
#include <md_chem.h>
#include <md_dssp.h>
#include <md_match.h>

#include <core/md_coord_stream.h>
#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_intrinsics.h>
#include <core/md_simd.h>
#include <core/md_spatial_acc.h>
#include <core/md_bitset.h>
#include <core/md_os.h>
#include <core/md_hash.h>
#include <core/md_fifo.h>

#include <math.h>
#include <string.h>
#include <float.h>

//#define PROFILE


#ifdef __cplusplus
extern "C" {
#endif

// Tables
// The zero slot is 'Xx', not 'Unknown': winioctl.h declares an enumerator by the latter name at
// file scope, so on any build whose include chain reaches windows.h the two collide. 'Xx' is what
// element_symbols[0] holds anyway.
enum {
    Xx = 0,
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V,
    Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru,
    Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb,
    Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh,
    Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og, Num_Elements
};

static const str_t element_symbols[] = {
    STR_INIT("Xx"), STR_INIT("H"),  STR_INIT("He"), STR_INIT("Li"), STR_INIT("Be"), STR_INIT("B"),  STR_INIT("C"),  STR_INIT("N"),  STR_INIT("O"),  STR_INIT("F"),  STR_INIT("Ne"), STR_INIT("Na"), STR_INIT("Mg"), STR_INIT("Al"), STR_INIT("Si"), STR_INIT("P"),  STR_INIT("S"),  STR_INIT("Cl"), STR_INIT("Ar"), STR_INIT("K"),  STR_INIT("Ca"), STR_INIT("Sc"), STR_INIT("Ti"), STR_INIT("V"),
    STR_INIT("Cr"), STR_INIT("Mn"), STR_INIT("Fe"), STR_INIT("Co"), STR_INIT("Ni"), STR_INIT("Cu"), STR_INIT("Zn"), STR_INIT("Ga"), STR_INIT("Ge"), STR_INIT("As"), STR_INIT("Se"), STR_INIT("Br"), STR_INIT("Kr"), STR_INIT("Rb"), STR_INIT("Sr"), STR_INIT("Y"),  STR_INIT("Zr"), STR_INIT("Nb"), STR_INIT("Mo"), STR_INIT("Tc"), STR_INIT("Ru"), STR_INIT("Rh"), STR_INIT("Pd"), STR_INIT("Ag"),
    STR_INIT("Cd"), STR_INIT("In"), STR_INIT("Sn"), STR_INIT("Sb"), STR_INIT("Te"), STR_INIT("I"),  STR_INIT("Xe"), STR_INIT("Cs"), STR_INIT("Ba"), STR_INIT("La"), STR_INIT("Ce"), STR_INIT("Pr"), STR_INIT("Nd"), STR_INIT("Pm"), STR_INIT("Sm"), STR_INIT("Eu"), STR_INIT("Gd"), STR_INIT("Tb"), STR_INIT("Dy"), STR_INIT("Ho"), STR_INIT("Er"), STR_INIT("Tm"), STR_INIT("Yb"), STR_INIT("Lu"),
    STR_INIT("Hf"), STR_INIT("Ta"), STR_INIT("W"),  STR_INIT("Re"), STR_INIT("Os"), STR_INIT("Ir"), STR_INIT("Pt"), STR_INIT("Au"), STR_INIT("Hg"), STR_INIT("Tl"), STR_INIT("Pb"), STR_INIT("Bi"), STR_INIT("Po"), STR_INIT("At"), STR_INIT("Rn"), STR_INIT("Fr"), STR_INIT("Ra"), STR_INIT("Ac"), STR_INIT("Th"), STR_INIT("Pa"), STR_INIT("U"),  STR_INIT("Np"), STR_INIT("Pu"), STR_INIT("Am"),
    STR_INIT("Cm"), STR_INIT("Bk"), STR_INIT("Cf"), STR_INIT("Es"), STR_INIT("Fm"), STR_INIT("Md"), STR_INIT("No"), STR_INIT("Lr"), STR_INIT("Rf"), STR_INIT("Db"), STR_INIT("Sg"), STR_INIT("Bh"), STR_INIT("Hs"), STR_INIT("Mt"), STR_INIT("Ds"), STR_INIT("Rg"), STR_INIT("Cn"), STR_INIT("Nh"), STR_INIT("Fl"), STR_INIT("Mc"), STR_INIT("Lv"), STR_INIT("Ts"), STR_INIT("Og"),
};


static const str_t element_names[] = {
    STR_INIT("Unknown"),     STR_INIT("Hydrogen"),     STR_INIT("Helium"),       STR_INIT("Lithium"),     STR_INIT("Beryllium"),   STR_INIT("Boron"),         STR_INIT("Carbon"),     STR_INIT("Nitrogen"),   STR_INIT("Oxygen"),
    STR_INIT("Fluorine"),    STR_INIT("Neon"),         STR_INIT("Sodium"),       STR_INIT("Magnesium"),   STR_INIT("Aluminium"),   STR_INIT("Silicon"),       STR_INIT("Phosphorus"), STR_INIT("Sulfur"),     STR_INIT("Chlorine"),
    STR_INIT("Argon"),       STR_INIT("Potassium"),    STR_INIT("Calcium"),      STR_INIT("Scandium"),    STR_INIT("Titanium"),    STR_INIT("Vanadium"),      STR_INIT("Chromium"),   STR_INIT("Manganese"),  STR_INIT("Iron"),
    STR_INIT("Cobalt"),      STR_INIT("Nickel"),       STR_INIT("Copper"),       STR_INIT("Zinc"),        STR_INIT("Gallium"),     STR_INIT("Germanium"),     STR_INIT("Arsenic"),    STR_INIT("Selenium"),   STR_INIT("Bromine"),
    STR_INIT("Krypton"),     STR_INIT("Rubidium"),     STR_INIT("Strontium"),    STR_INIT("Yttrium"),     STR_INIT("Zirconium"),   STR_INIT("Niobium"),       STR_INIT("Molybdenum"), STR_INIT("Technetium"), STR_INIT("Ruthenium"),
    STR_INIT("Rhodium"),     STR_INIT("Palladium"),    STR_INIT("Silver"),       STR_INIT("Cadmium"),     STR_INIT("Indium"),      STR_INIT("Tin"),           STR_INIT("Antimony"),   STR_INIT("Tellurium"),  STR_INIT("Iodine"),
    STR_INIT("Xenon"),       STR_INIT("Caesium"),      STR_INIT("Barium"),       STR_INIT("Lanthanum"),   STR_INIT("Cerium"),      STR_INIT("Praseodymium"),  STR_INIT("Neodymium"),  STR_INIT("Promethium"), STR_INIT("Samarium"),
    STR_INIT("Europium"),    STR_INIT("Gadolinium"),   STR_INIT("Terbium"),      STR_INIT("Dysprosium"),  STR_INIT("Holmium"),     STR_INIT("Erbium"),        STR_INIT("Thulium"),    STR_INIT("Ytterbium"),  STR_INIT("Lutetium"),
    STR_INIT("Hafnium"),     STR_INIT("Tantalum"),     STR_INIT("Tungsten"),     STR_INIT("Rhenium"),     STR_INIT("Osmium"),      STR_INIT("Iridium"),       STR_INIT("Platinum"),   STR_INIT("Gold"),       STR_INIT("Mercury"),
    STR_INIT("Thallium"),    STR_INIT("Lead"),         STR_INIT("Bismuth"),      STR_INIT("Polonium"),    STR_INIT("Astatine"),    STR_INIT("Radon"),         STR_INIT("Francium"),   STR_INIT("Radium"),     STR_INIT("Actinium"),
    STR_INIT("Thorium"),     STR_INIT("Protactinium"), STR_INIT("Uranium"),      STR_INIT("Neptunium"),   STR_INIT("Plutonium"),   STR_INIT("Americium"),     STR_INIT("Curium"),     STR_INIT("Berkelium"),  STR_INIT("Californium"),
    STR_INIT("Einsteinium"), STR_INIT("Fermium"),      STR_INIT("Mendelevium"),  STR_INIT("Nobelium"),    STR_INIT("Lawrencium"),  STR_INIT("Rutherfordium"), STR_INIT("Dubnium"),    STR_INIT("Seaborgium"), STR_INIT("Bohrium"),
    STR_INIT("Hassium"),     STR_INIT("Meitnerium"),   STR_INIT("Darmstadtium"), STR_INIT("Roentgenium"), STR_INIT("Copernicium"), STR_INIT("Nihonium"),      STR_INIT("Flerovium"),  STR_INIT("Moscovium"),  STR_INIT("Livermorium"),
    STR_INIT("Tennessine"),  STR_INIT("Oganesson"),                          
};

// http://dx.doi.org/10.1039/b801115j
// Covalent radii (in � �100) � Cordero et al., Dalton Trans., 2008, 2832�2838
static const uint8_t element_covalent_radii_u8[] = {
      0,  31,  28, 128,  96,  84,  76,  71,  66,  57,  58, 166, 141, 121, 111, 107, 105, 102, 106, 203, 176, 170, 160, 153,
    139, 139, 132, 126, 124, 132, 122, 122, 120, 119, 120, 120, 116, 220, 195, 190, 175, 164, 154, 147, 146, 142, 139, 145,
    144, 142, 139, 139, 138, 139, 140, 244, 215, 207, 204, 203, 201, 199, 198, 198, 196, 194, 192, 192, 189, 190, 187, 187,
    175, 170, 162, 151, 144, 141, 136, 136, 132, 145, 146, 148, 140, 150, 150, 255, 221, 215, 206, 200, 196, 190, 187, 180,
    169, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160
};

// Approximate ionic radii (� �100), mostly Shannon radii (CN�6, common oxidation state)
static const uint8_t element_ionic_radii_u8[] = {
    0,
    0,   0,  76,  45,  27,  16, 146, 140, 133,   0,
    102,  72,  54,  40, 212, 184, 181,   0, 138, 100,
    75,  61,  54,  62,  83,  78,  74,  69,  73,  74,
    62,  53, 222, 198, 196,   0, 152, 118,  90,  72,
    64,  59,  56,  68,  67,  86, 115,  95,  80,  69,
    245, 221, 220,   0, 167, 135, 103, 101,  99,  98,
    97,  96,  95,  94,  92,  91,  90,  89,  88,  87,
    86,  71,  64,  60,  56,  77,  77,  80, 137, 102,
    150, 119, 103, 230, 227,   0, 180, 148, 112,  94,
    89, 100,  98,  96, 109, 108, 107, 106, 105, 104,
    103, 110, 103,  78,  74,  69,  67,  66,   0,   0,
    0, 100,   0, 120,   0, 240, 230,   0
};

#if 0
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
#endif

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

// Prefixed, because wingdi.h defines its own RGB() with the opposite byte order. Colliding with it
// is not just a warning: whichever definition wins decides whether these tables are RGBA or BGR.
#define CPK_RGBA(r,g,b,a) ((((uint32_t)a & 255) << 24) | (((uint32_t)b & 255) << 16) | (((uint32_t)g & 255) << 8) | ((uint32_t)r & 255) )
#define CPK_RGB(r,g,b) CPK_RGBA(r,g,b,255)

// Based on this with some values modified
// Some stem from Jmol, some from Rasmol
// http://jmol.sourceforge.net/jscolors/

static const uint32_t element_cpk_colors[] = {
    CPK_RGB(255,   0, 255), // Unknown (0)  rgb(255, 0, 255)
    CPK_RGB(255, 255, 255), // H (1)        rgb(255, 255, 255)
    CPK_RGB(217, 255, 255), // He (2)       rgb(217, 255, 255)
    CPK_RGB(226, 152, 255), // Li (3)       rgb(226, 152, 255)
    CPK_RGB(194, 255,   0), // Be (4)       rgb(194, 255, 0)
    CPK_RGB(254, 162, 176), // B (5)        rgb(254, 162, 176)
    CPK_RGB(131, 131, 131), // C (6)        rgb(131, 131, 131)
    CPK_RGB( 45,  95, 253), // N (7)        rgb(45, 95, 253)
    CPK_RGB(240,   0,   0), // O (8)        rgb(240, 0, 0)
    CPK_RGB(182, 217, 124), // F (9)        rgb(182, 217, 124)
    CPK_RGB(179, 227, 245), // Ne (10)      rgb(179, 227, 245)
    CPK_RGB(177, 110, 255), // Na (11)      rgb(177, 110, 255)
    CPK_RGB(135, 249,  58), // Mg (12)      rgb(135, 249, 58)
    CPK_RGB(173, 160, 159), // Al (13)      rgb(173, 160, 159)
    CPK_RGB(255, 203, 168), // Si (14)      rgb(255, 203, 168)
    CPK_RGB(255, 157,  41), // P (15)       rgb(255, 157, 41)
    CPK_RGB(253, 237,   1), // S (16)       rgb(253, 237, 1)
    CPK_RGB( 56, 204,   2), // Cl (17)      rgb(56, 204, 2)
    CPK_RGB(128, 209, 227), // Ar (18)      rgb(128, 209, 227)
    CPK_RGB(144,  86, 180), // K (19)       rgb(144, 86, 180)
    CPK_RGB(  2, 194, 199), // Ca (20)      rgb(2, 194, 199)
    CPK_RGB(230, 230, 230), // Sc (21)      rgb(230, 230, 230)
    CPK_RGB(194, 199, 220), // Ti (22)      rgb(194, 199, 220)
    CPK_RGB(166, 166, 171), // V (23)       rgb(166, 166, 171)
    CPK_RGB(138, 153, 199), // Cr (24)      rgb(138, 153, 199)
    CPK_RGB(158, 141, 209), // Mn (25)      rgb(158, 141, 209)
    CPK_RGB(236, 104,  36), // Fe (26)      rgb(236, 104, 36)
    CPK_RGB(247, 140, 144), // Co (27)      rgb(247, 140, 144)
    CPK_RGB( 87, 171, 102), // Ni (28)      rgb(87, 171, 102)
    CPK_RGB(173, 122,  63), // Cu (29)      rgb(173, 122, 63)
    CPK_RGB( 96, 126, 194), // Zn (30)      rgb(96, 126, 194)
    CPK_RGB(188, 125, 136), // Ga (31)      rgb(188, 125, 136)
    CPK_RGB( 68, 124, 132), // Ge (32)      rgb(68, 124, 132)
    CPK_RGB(189, 128, 227), // As (33)      rgb(189, 128, 227)
    CPK_RGB(184, 157,  71), // Se (34)      rgb(184, 157, 71)
    CPK_RGB(144,  18,  13), // Br (35)      rgb(144, 18, 13)
    CPK_RGB( 92, 184, 209), // Kr (36)      rgb(92, 184, 209)
    CPK_RGB( 96,  47, 175), // Rb (37)      rgb(96, 47, 175)
    CPK_RGB(  3, 161,  81), // Sr (38)      rgb(3, 161, 81)
    CPK_RGB(148, 255, 255), // Y (39)       rgb(148, 255, 255)
    CPK_RGB(143, 235, 223), // Zr (40)      rgb(143, 235, 223)
    CPK_RGB(115, 194, 201), // Nb (41)      rgb(115, 194, 201)
    CPK_RGB( 35, 155, 164), // Mo (42)      rgb(35, 155, 164)
    CPK_RGB( 59, 158, 158), // Tc (43)      rgb(59, 158, 158)
    CPK_RGB( 36, 143, 143), // Ru (44)      rgb(36, 143, 143)
    CPK_RGB( 10, 125, 140), // Rh (45)      rgb(10, 125, 140)
    CPK_RGB(  0, 101, 141), // Pd (46)      rgb(0, 101, 141)
    CPK_RGB(182, 182, 182), // Ag (47)      rgb(182, 182, 182)
    CPK_RGB(255, 229, 142), // Cd (48)      rgb(255, 229, 142)
    CPK_RGB(158,  99,  87), // In (49)      rgb(158, 99, 87)
    CPK_RGB( 62, 107, 103), // Sn (50)      rgb(62, 107, 103)
    CPK_RGB(158,  99, 181), // Sb (51)      rgb(158, 99, 181)
    CPK_RGB(206, 113,   1), // Te (52)      rgb(206, 113, 1)
    CPK_RGB(150,   0, 140), // I (53)       rgb(150, 0, 140)
    CPK_RGB( 66, 158, 176), // Xe (54)      rgb(66, 158, 176)
    CPK_RGB( 80,   0, 113), // Cs (55)      rgb(80, 0, 113)
    CPK_RGB(  0, 112,  78), // Ba (56)      rgb(0, 112, 78)
    CPK_RGB(112, 212, 255), // La (57)      rgb(112, 212, 255)
    CPK_RGB(255, 255, 199), // Ce (58)      rgb(255, 255, 199)
    CPK_RGB(217, 255, 199), // Pr (59)      rgb(217, 255, 199)
    CPK_RGB(199, 255, 199), // Nd (60)      rgb(199, 255, 199)
    CPK_RGB(163, 255, 199), // Pm (61)      rgb(163, 255, 199)
    CPK_RGB(143, 255, 199), // Sm (62)      rgb(143, 255, 199)
    CPK_RGB( 97, 255, 199), // Eu (63)      rgb(97, 255, 199)
    CPK_RGB( 69, 255, 199), // Gd (64)      rgb(69, 255, 199)
    CPK_RGB( 48, 255, 199), // Tb (65)      rgb(48, 255, 199)
    CPK_RGB( 31, 255, 199), // Dy (66)      rgb(31, 255, 199)
    CPK_RGB(  0, 255, 156), // Ho (67)      rgb(0, 255, 156)
    CPK_RGB(  0, 230, 117), // Er (68)      rgb(0, 230, 117)
    CPK_RGB(  0, 212,  82), // Tm (69)      rgb(0, 212, 82)
    CPK_RGB(  0, 191,  56), // Yb (70)      rgb(0, 191, 56)
    CPK_RGB(  0, 171,  36), // Lu (71)      rgb(0, 171, 36)
    CPK_RGB( 77, 194, 255), // Hf (72)      rgb(77, 194, 255)
    CPK_RGB( 77, 166, 255), // Ta (73)      rgb(77, 166, 255)
    CPK_RGB( 18, 150, 225), // W (74)       rgb(18, 150, 225)
    CPK_RGB( 38, 125, 171), // Re (75)      rgb(38, 125, 171)
    CPK_RGB( 38, 102, 150), // Os (76)      rgb(38, 102, 150)
    CPK_RGB( 23,  84, 135), // Ir (77)      rgb(23, 84, 135)
    CPK_RGB(226, 221, 247), // Pt (78)      rgb(226, 221, 247)
    CPK_RGB(255, 209,  35), // Au (79)      rgb(255, 209, 35)
    CPK_RGB(184, 184, 208), // Hg (80)      rgb(184, 184, 208)
    CPK_RGB(166,  84,  77), // Tl (81)      rgb(166, 84, 77)
    CPK_RGB( 73,  69,  86), // Pb (82)      rgb(73, 69, 86)
    CPK_RGB(158,  79, 181), // Bi (83)      rgb(158, 79, 181)
    CPK_RGB(171,  92,   0), // Po (84)      rgb(171, 92, 0)
    CPK_RGB(117,  79,  69), // At (85)      rgb(117, 79, 69)
    CPK_RGB( 66, 130, 150), // Rn (86)      rgb(66, 130, 150)
    CPK_RGB( 66,   0, 102), // Fr (87)      rgb(66, 0, 102)
    CPK_RGB(  0, 125,   0), // Ra (88)      rgb(0, 125, 0)
    CPK_RGB(112, 171, 250), // Ac (89)      rgb(112, 171, 250)
    CPK_RGB(  0, 186, 255), // Th (90)      rgb(0, 186, 255)
    CPK_RGB(  0, 161, 255), // Pa (91)      rgb(0, 161, 255)
    CPK_RGB(  0, 143, 255), // U (92)       rgb(0, 143, 255)
    CPK_RGB(  0, 128, 255), // Np (93)      rgb(0, 128, 255)
    CPK_RGB(  0, 107, 255), // Pu (94)      rgb(0, 107, 255)
    CPK_RGB( 84,  92, 242), // Am (95)      rgb(84, 92, 242)
    CPK_RGB(120,  92, 227), // Cm (96)      rgb(120, 92, 227)
    CPK_RGB(138,  79, 227), // Bk (97)      rgb(138, 79, 227)
    CPK_RGB(161,  54, 212), // Cf (98)      rgb(161, 54, 212)
    CPK_RGB(179,  31, 212), // Es (99)      rgb(179, 31, 212)
    CPK_RGB(179,  31, 186), // Fm (100)     rgb(179, 31, 186)
    CPK_RGB(179,  13, 166), // Md (101)     rgb(179, 13, 166)
    CPK_RGB(189,  13, 135), // No (102)     rgb(189, 13, 135)
    CPK_RGB(199,   0, 102), // Lr (103)     rgb(199, 0, 102)
    CPK_RGB(204,   0,  89), // Rf (104)     rgb(204, 0, 89)
    CPK_RGB(209,   0,  79), // Db (105)     rgb(209, 0, 79)
    CPK_RGB(217,   0,  69), // Sg (106)     rgb(217, 0, 69)
    CPK_RGB(224,   0,  56), // Bh (107)     rgb(224, 0, 56)
    CPK_RGB(230,   0,  46), // Hs (108)     rgb(230, 0, 46)
    CPK_RGB(235,   0,  38), // Mt (109)     rgb(235, 0, 38)
    CPK_RGB(240,   0,  34), // Ds (110)     rgb(240, 0, 34)
    CPK_RGB(246,   0,  32), // Rg (111)     rgb(246, 0, 32)
    CPK_RGB(248,   0,  30), // Cn (112)     rgb(248, 0, 30)
    CPK_RGB(250,   0,  28), // Nh (113)     rgb(250, 0, 28)
    CPK_RGB(252,   0,  26), // Fl (114)     rgb(252, 0, 26)
    CPK_RGB(253,   0,  24), // Mc (115)     rgb(253, 0, 24)
    CPK_RGB(254,   0,  22), // Lv (116)     rgb(254, 0, 22)
    CPK_RGB(255,   0,  20), // Ts (117)     rgb(255, 0, 20)
    CPK_RGB(255,   0,  18), // Og (118)     rgb(255, 0, 18)
};

// Structure to hold string and color to represent predefined color mappings for e.g. amino acids
typedef struct {
    const char* name;
    uint32_t color;
} named_color_t;

// Table for predefined amino acid color mappings
named_color_t amino_acid_colors[] = {
    { "ALA", CPK_RGB(200, 200, 200) },  // Alanine                      rgb(200, 200, 200)
    { "ARG", CPK_RGB(20, 90, 255) },    // Arginine                     rgb(20, 90, 255)
    { "ASN", CPK_RGB(220, 220, 220) },  // Asparagine                   rgb(220, 220, 220)
    { "ASP", CPK_RGB(230, 10, 10) },    // Aspartic acid                rgb(230, 10, 10)
    { "CYS", CPK_RGB(255, 255, 0) },    // Cysteine                     rgb(255, 255, 0)
    { "GLN", CPK_RGB(220, 220, 220) },  // Glutamine                    rgb(220, 220, 220)
    { "GLU", CPK_RGB(230, 10, 10) },    // Glutamic acid                rgb(230, 10, 10)
    { "GLY", CPK_RGB(235, 235, 235) },  // Glycine                      rgb(235, 235, 235)
    { "HIS", CPK_RGB(130, 130, 210) },  // Histidine                    rgb(130, 130, 210)
    { "ILE", CPK_RGB(200, 150, 100) },  // Isoleucine                   rgb(200, 150, 100)
    { "LEU", CPK_RGB(200, 150, 100) },  // Leucine                      rgb(200, 150, 100)
    { "LYS", CPK_RGB(20, 90, 255) },    // Lysine                       rgb(20, 90, 255)
    { "MET", CPK_RGB(255, 150, 50) },   // Methionine                   rgb(255, 150, 50)
    { "PHE", CPK_RGB(50, 50, 170) },    // Phenylalanine                rgb(50, 50, 170)
    { "PRO", CPK_RGB(220, 150, 130) },  // Proline                      rgb(220, 150, 130)
    { "SER", CPK_RGB(250, 150, 0) },    // Serine                       rgb(250, 150, 0)
    { "THR", CPK_RGB(250, 150, 0) },    // Threonine                    rgb(250, 150, 0)
    { "TRP", CPK_RGB(180, 90, 180) },   // Tryptophan                   rgb(180, 90, 180)
    { "TYR", CPK_RGB(50, 50, 170) },    // Tyrosine                     rgb(50, 50, 170)
    { "VAL", CPK_RGB(200, 150, 100) },  // Valine                       rgb(200, 150, 100)
    { "SEC", CPK_RGB(255, 255, 0) },    // Selenocysteine               rgb(255, 255, 0)

    // CCD
    { "MSE", CPK_RGB(255, 150, 50) },   // Selenomethionine             rgb(255, 150, 50)
    { "SEP", CPK_RGB(250, 150, 0) },    // Phosphoserine                rgb(250, 150, 0)
    { "TPO", CPK_RGB(250, 150, 0) },    // Phosphothreonine             rgb(250, 150, 0)
    { "PTR", CPK_RGB(250, 150, 0) },    // Phosphotyrosine              rgb(250, 150, 0)

    // Charmm
    { "HSD", CPK_RGB(130, 130, 210) },  // Histidine (delta)            rgb(130, 130, 210)
    { "HSP", CPK_RGB(130, 130, 210) },  // Histidine (epsilon)          rgb(130, 130, 210)
    { "LSN", CPK_RGB(220, 220, 220) },  // Asparagine (protonated)      rgb(220, 220, 220)
    { "ASPP", CPK_RGB(230, 10, 10) },   // Aspartic acid (protonated)   rgb(230, 10, 10)
    { "GLUP", CPK_RGB(230, 10, 10) },   // Glutamic acid (protonated)   rgb(230, 10, 10)

    // Amber
    { "HID", CPK_RGB(130, 130, 210) },  // Histidine (delta)            rgb(130, 130, 210)
    { "HIE", CPK_RGB(130, 130, 210) },  // Histidine (epsilon)          rgb(130, 130, 210)
    { "HIP", CPK_RGB(130, 130, 210) },  // Histidine (protonated)       rgb(130, 130, 210)
    { "LYN", CPK_RGB(20, 90, 255) },    // Lysine (deprotonated)        rgb(20, 90, 255)
    { "ASH", CPK_RGB(230, 10, 10) },    // Aspartic acid (protonated)   rgb(230, 10, 10)
    { "GLH", CPK_RGB(230, 10, 10) },    // Glutamic acid (protonated)   rgb(230, 10, 10)
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

// Groups 3-12: Sc-Zn, Y-Cd, Hf-Hg, Rf-Cn
static const uint64_t element_transition_metal_mask[2] = {
    0x0001FF807FE00000ULL,
    0x0001FF000001FF00ULL,
};

#if 0
static const uint64_t element_alkaline_earth[2] = {
    0x0100004000101010,
    0x0000000001000000,
};

static const uint64_t element_covalent_mask[2] = {
    0x003E0000001FC7E6ULL, // atomic numbers 1�63
    0x0000000000000000ULL, // atomic numbers 64�127 (unused)
};
#endif

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

GENERATE_MASK_FUNC(is_alkali, element_alkali_mask)

GENERATE_MASK_FUNC(is_transition_metal, element_transition_metal_mask)

GENERATE_MASK_FUNC(is_halogen, element_halogen_mask)

GENERATE_MASK_FUNC(is_aromatic, element_aromatic_ring)

static inline bool is_metal_donor(size_t atomic_nr) {
    switch(atomic_nr) {
    case O:
    case N:
    case S:
    case P:
    case F:
    case Cl:
    case Br:
    case I:
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
    "HSD", "HSE", "HSP", "LSN", "ASPP", "GLUP",

    // Amber
    "HID", "HIE", "HIP", "LYN", "ASH", "GLH",

    // Martini: the neutral Arg, Asp, Glu, Lys and the charged His of Martini 2 (HISH is GROMOS's too) and of Martini 3
    // beta (HIH), and the lipidated residues of Martini 3: farnesyl, geranylgeranyl and palmitoyl Cys, myristoyl Gly
    "ARG0", "ASP0", "GLU0", "LYS0", "HISH", "HIH",
    "CYSF", "CYSG", "CYSP", "GLYM",
};

#if 0
static const char* nucleic_acids[] = {
    "A", "C", "G", "T", "U", "+U", "YG", "1MA", "1MG", "2MG", "5MC", "5MU", "7MG", "H2U", "M2G", "OMC", "OMG", "PSU",
};
#endif

/*static const char* peptides[] = {"APN", "CPN", "TPN", "GPN"};*/
static const char* rna[] = {"A", "C", "T", "G", "I", "U", "N"};
static const char* dna[] = {"DA", "DC", "DT", "DG", "DI", "DU", "DN"};

static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4", "TIP5", "TIP4P", "TIP5P", "TP3", "TP4", "TP5", "T3P", "T4P", "T5P", "T4E", "OPC", "W", "DOD", "D30", "SPC" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

static const char* common_ions[] = { "NA", "K", "CA", "MG", "ZN", "CL", "F", "MN", "FE", "CU", "CO", "NI", "CD", "BR", "I", "CS", "SR"};

#if 0
// Taken from here
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/model/properties/atomic/bonds.ts
/*
* Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
* @author Alexander Rose <alexander.rose@weirdbyte.de>
*/

// If an entry exists in this table, it has order 2
static str_t intra_bond_order_table[] = {
    STR_INIT("HIS|CD2|CG"),
    STR_INIT("HIS|CE1|ND1"),
    STR_INIT("ARG|CZ|NH2"),
    STR_INIT("PHE|CE1|CZ"),
    STR_INIT("PHE|CD2|CE2"),
    STR_INIT("PHE|CD1|CG"),
    STR_INIT("TRP|CD1|CG"),
    STR_INIT("TRP|CD2|CE2"),
    STR_INIT("TRP|CE3|CZ3"),
    STR_INIT("TRP|CH2|CZ2"),
    STR_INIT("ASN|CG|OD1"),
    STR_INIT("GLN|CD|OE1"),
    STR_INIT("TYR|CD1|CG"),
    STR_INIT("TYR|CD2|CE2"),
    STR_INIT("TYR|CE1|CZ"),
    STR_INIT("ASP|CG|OD1"),
    STR_INIT("GLU|CD|OE1"),

    STR_INIT("G|C8|N7"),
    STR_INIT("G|C4|C5"),
    STR_INIT("G|C2|N3"),
    STR_INIT("G|C6|O6"),
    STR_INIT("C|C4|N3"),
    STR_INIT("C|C5|C6"),
    STR_INIT("C|C2|O2"),
    STR_INIT("A|C2|N3"),
    STR_INIT("A|C6|N1"),
    STR_INIT("A|C4|C5"),
    STR_INIT("A|C8|N7"),
    STR_INIT("U|C5|C6"),
    STR_INIT("U|C2|O2"),
    STR_INIT("U|C4|O4"),

    STR_INIT("DG|C8|N7"),
    STR_INIT("DG|C4|C5"),
    STR_INIT("DG|C2|N3"),
    STR_INIT("DG|C6|O6"),
    STR_INIT("DC|C4|N3"),
    STR_INIT("DC|C5|C6"),
    STR_INIT("DC|C2|O2"),
    STR_INIT("DA|C2|N3"),
    STR_INIT("DA|C6|N1"),
    STR_INIT("DA|C4|C5"),
    STR_INIT("DA|C8|N7"),
    STR_INIT("DT|C5|C6"),
    STR_INIT("DT|C2|O2"),
    STR_INIT("DT|C4|O4"),
};
#endif

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

// ### RADIX SORT ###

// Least significant digit first, 8 bits a digit. The histograms of every digit come from one read of the keys, and a pass
// whose digit all keys share would leave the order as it is, so it is skipped: keys with constant upper bytes (indices,
// pairs of them packed) take only the passes they need. Up to RADIX_SMALL_COUNT keys an insertion sort is faster than
// the histograms, and keys already in order are left as they are.
//
// The histograms have a row per bucket and a column per digit, so the prefix sums of every digit are one vector sum per
// row. The argsorts move each key together with its index, packed in one value, so a pass reads them in order and writes
// one stream per bucket, rather than looking up the key of each index.

#define RADIX_SMALL_COUNT 64

static void insertion_sort_u32(uint32_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const uint32_t k = data[i];
        size_t j = i;
        while (j > 0 && data[j - 1] > k) {
            data[j] = data[j - 1];
            --j;
        }
        data[j] = k;
    }
}

static void insertion_sort_u64(uint64_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const uint64_t k = data[i];
        size_t j = i;
        while (j > 0 && data[j - 1] > k) {
            data[j] = data[j - 1];
            --j;
        }
        data[j] = k;
    }
}

// The key above the index: the values are unique, in the order of a stable sort
static void insertion_argsort_u32(uint32_t* out_indices, const uint32_t* keys, size_t count) {
    ASSERT(count <= RADIX_SMALL_COUNT);
    uint64_t v[RADIX_SMALL_COUNT];
    for (size_t i = 0; i < count; ++i) {
        v[i] = ((uint64_t)keys[i] << 32) | i;
    }
    insertion_sort_u64(v, count);
    for (size_t i = 0; i < count; ++i) {
        out_indices[i] = (uint32_t)v[i];
    }
}

static void insertion_argsort_u64(uint32_t* out_indices, const uint64_t* keys, size_t count) {
    ASSERT(count <= RADIX_SMALL_COUNT);
    uint64_t k[RADIX_SMALL_COUNT];
    for (size_t i = 0; i < count; ++i) {
        const uint64_t key = keys[i];
        size_t j = i;
        // Only past greater keys: equal keys keep their order
        while (j > 0 && k[j - 1] > key) {
            k[j] = k[j - 1];
            out_indices[j] = out_indices[j - 1];
            --j;
        }
        k[j] = key;
        out_indices[j] = (uint32_t)i;
    }
}

static bool is_sorted_u32(const uint32_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        if (data[i] < data[i - 1]) return false;
    }
    return true;
}

static bool is_sorted_u64(const uint64_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        if (data[i] < data[i - 1]) return false;
    }
    return true;
}

static void radix_hist_u32(uint32_t hist[256][4], const uint32_t* keys, size_t count) {
    MEMSET(hist, 0, sizeof(uint32_t) * 256 * 4);
    for (size_t i = 0; i < count; ++i) {
        const uint32_t k = keys[i];
        hist[(k >>  0) & 255][0]++;
        hist[(k >>  8) & 255][1]++;
        hist[(k >> 16) & 255][2]++;
        hist[(k >> 24) & 255][3]++;
    }
}

static void radix_hist_u64(uint32_t hist[256][8], const uint64_t* keys, size_t count) {
    MEMSET(hist, 0, sizeof(uint32_t) * 256 * 8);
    for (size_t i = 0; i < count; ++i) {
        const uint64_t k = keys[i];
        hist[(k >>  0) & 255][0]++;
        hist[(k >>  8) & 255][1]++;
        hist[(k >> 16) & 255][2]++;
        hist[(k >> 24) & 255][3]++;
        hist[(k >> 32) & 255][4]++;
        hist[(k >> 40) & 255][5]++;
        hist[(k >> 48) & 255][6]++;
        hist[(k >> 56) & 255][7]++;
    }
}

// The digits whose pass reorders anything, then the counts turned into the offsets of the buckets. Any key tells whether
// a digit is shared by all: then its bucket holds every key.
static int radix_offsets_u32(int passes[4], uint32_t hist[256][4], uint32_t any_key, size_t count) {
    int num_passes = 0;
    for (int p = 0; p < 4; ++p) {
        if (hist[(any_key >> (p * 8)) & 255][p] != count) passes[num_passes++] = p;
    }
    md_128i sum = md_mm_setzero_si128();
    for (size_t i = 0; i < 256; ++i) {
        const md_128i c = md_mm_loadu_si128(hist[i]);
        md_mm_storeu_si128(hist[i], sum);
        sum = md_mm_add_epi32(sum, c);
    }
    return num_passes;
}

static int radix_offsets_u64(int passes[8], uint32_t hist[256][8], uint64_t any_key, size_t count) {
    int num_passes = 0;
    for (int p = 0; p < 8; ++p) {
        if (hist[(any_key >> (p * 8)) & 255][p] != count) passes[num_passes++] = p;
    }
    md_128i sum_lo = md_mm_setzero_si128();
    md_128i sum_hi = md_mm_setzero_si128();
    for (size_t i = 0; i < 256; ++i) {
        const md_128i lo = md_mm_loadu_si128(hist[i] + 0);
        const md_128i hi = md_mm_loadu_si128(hist[i] + 4);
        md_mm_storeu_si128(hist[i] + 0, sum_lo);
        md_mm_storeu_si128(hist[i] + 4, sum_hi);
        sum_lo = md_mm_add_epi32(sum_lo, lo);
        sum_hi = md_mm_add_epi32(sum_hi, hi);
    }
    return num_passes;
}

// The offsets of one digit, out of its column into a table of their own, for the pass to keep close
static inline void radix_column(uint32_t offset[256], const uint32_t* hist, int num_digits, int digit) {
    for (int b = 0; b < 256; ++b) {
        offset[b] = hist[b * num_digits + digit];
    }
}

static inline void radix_pass_u32(uint32_t* dst, const uint32_t* src, size_t count, uint32_t hist[256][4], int digit) {
    const int shift = digit * 8;
    for (size_t i = 0; i < count; ++i) {
        const uint32_t k = src[i];
        dst[hist[(k >> shift) & 255][digit]++] = k;
    }
}

static void sort_radix_u32(uint32_t* data, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][4];
    radix_hist_u32(hist, data, count);
    int passes[4];
    const int num_passes = radix_offsets_u32(passes, hist, data[0], count);
    if (num_passes == 0) return;

    uint32_t* src = data;
    uint32_t* dst = md_temp_alloc_array(temp, uint32_t, count);
    for (int p = 0; p < num_passes; ++p) {
        // The digit a constant in each: these loops are the time the sort takes
        switch (passes[p]) {
        case 0: radix_pass_u32(dst, src, count, hist, 0); break;
        case 1: radix_pass_u32(dst, src, count, hist, 1); break;
        case 2: radix_pass_u32(dst, src, count, hist, 2); break;
        case 3: radix_pass_u32(dst, src, count, hist, 3); break;
        default: ASSERT(false); break;
        }
        uint32_t* tmp = src;
        src = dst;
        dst = tmp;
    }

    // An odd number of passes leaves the keys in the scratch
    if (src != data) {
        MEMCPY(data, src, count * sizeof(uint32_t));
    }
}

static void sort_radix_u64(uint64_t* data, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][8];
    radix_hist_u64(hist, data, count);
    int passes[8];
    const int num_passes = radix_offsets_u64(passes, hist, data[0], count);
    if (num_passes == 0) return;

    uint64_t* src = data;
    uint64_t* dst = md_temp_alloc_array(temp, uint64_t, count);
    uint32_t offset[256];
    for (int p = 0; p < num_passes; ++p) {
        const int shift = passes[p] * 8;
        radix_column(offset, hist[0], 8, passes[p]);
        for (size_t i = 0; i < count; ++i) {
            const uint64_t k = src[i];
            dst[offset[(k >> shift) & 255]++] = k;
        }
        uint64_t* tmp = src;
        src = dst;
        dst = tmp;
    }

    if (src != data) {
        MEMCPY(data, src, count * sizeof(uint64_t));
    }
}

// Each key with its index below it in one 64-bit value: the passes sort on the upper half, and the last keeps the index
static void argsort_radix_u32(uint32_t* out_indices, const uint32_t* keys, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][4];
    radix_hist_u32(hist, keys, count);
    int passes[4];
    const int num_passes = radix_offsets_u32(passes, hist, keys[0], count);
    uint32_t offset[256];

    if (num_passes == 0) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }
    if (num_passes == 1) {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 4, passes[0]);
        for (size_t i = 0; i < count; ++i) {
            out_indices[offset[(keys[i] >> shift) & 255]++] = (uint32_t)i;
        }
        return;
    }

    uint64_t* buf[2] = {
        md_temp_alloc_array(temp, uint64_t, count),
        num_passes > 2 ? md_temp_alloc_array(temp, uint64_t, count) : NULL,
    };

    {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 4, passes[0]);
        uint64_t* dst = buf[0];
        for (size_t i = 0; i < count; ++i) {
            const uint32_t k = keys[i];
            dst[offset[(k >> shift) & 255]++] = ((uint64_t)k << 32) | i;
        }
    }
    for (int p = 1; p < num_passes - 1; ++p) {
        const int shift = 32 + passes[p] * 8;
        radix_column(offset, hist[0], 4, passes[p]);
        const uint64_t* src = buf[(p - 1) & 1];
        uint64_t* dst = buf[p & 1];
        for (size_t i = 0; i < count; ++i) {
            const uint64_t v = src[i];
            dst[offset[(v >> shift) & 255]++] = v;
        }
    }
    {
        const int shift = 32 + passes[num_passes - 1] * 8;
        radix_column(offset, hist[0], 4, passes[num_passes - 1]);
        const uint64_t* src = buf[(num_passes - 2) & 1];
        for (size_t i = 0; i < count; ++i) {
            const uint64_t v = src[i];
            out_indices[offset[(v >> shift) & 255]++] = (uint32_t)v;
        }
    }
}

// A 64-bit key leaves no room for its index: the two go together as a record of 16 bytes
typedef struct radix_record_u64_t {
    uint64_t key;
    uint64_t idx;
} radix_record_u64_t;

static void argsort_radix_u64(uint32_t* out_indices, const uint64_t* keys, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][8];
    radix_hist_u64(hist, keys, count);
    int passes[8];
    const int num_passes = radix_offsets_u64(passes, hist, keys[0], count);
    uint32_t offset[256];

    if (num_passes == 0) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }
    if (num_passes == 1) {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 8, passes[0]);
        for (size_t i = 0; i < count; ++i) {
            out_indices[offset[(keys[i] >> shift) & 255]++] = (uint32_t)i;
        }
        return;
    }

    radix_record_u64_t* buf[2] = {
        md_temp_alloc_array(temp, radix_record_u64_t, count),
        num_passes > 2 ? md_temp_alloc_array(temp, radix_record_u64_t, count) : NULL,
    };

    {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 8, passes[0]);
        radix_record_u64_t* dst = buf[0];
        for (size_t i = 0; i < count; ++i) {
            const uint64_t k = keys[i];
            dst[offset[(k >> shift) & 255]++] = (radix_record_u64_t){ k, i };
        }
    }
    for (int p = 1; p < num_passes - 1; ++p) {
        const int shift = passes[p] * 8;
        radix_column(offset, hist[0], 8, passes[p]);
        const radix_record_u64_t* src = buf[(p - 1) & 1];
        radix_record_u64_t* dst = buf[p & 1];
        for (size_t i = 0; i < count; ++i) {
            const radix_record_u64_t r = src[i];
            dst[offset[(r.key >> shift) & 255]++] = r;
        }
    }
    {
        const int shift = passes[num_passes - 1] * 8;
        radix_column(offset, hist[0], 8, passes[num_passes - 1]);
        const radix_record_u64_t* src = buf[(num_passes - 2) & 1];
        for (size_t i = 0; i < count; ++i) {
            out_indices[offset[(src[i].key >> shift) & 255]++] = (uint32_t)src[i].idx;
        }
    }
}

static md_array(uint64_t) make_bitfield(size_t num_bits, md_allocator_i* alloc) {
    size_t num_elem = DIV_UP(num_bits, 64);
    if (num_elem == 0) {
        return NULL;
    }
    md_array(uint64_t) bits = md_array_create(uint64_t, num_elem, alloc);
    MEMSET(bits, 0, md_array_bytes(bits));
    return bits;
}

static void bitfield_set_all(uint64_t* bits, size_t num_bits) {
    // No bits is a NULL bitfield (make_bitfield), which memset must not be handed even for zero bytes
    if (num_bits == 0) return;
    size_t num_bytes = DIV_UP(num_bits, 8);
    MEMSET(bits, 0xFF, num_bytes);
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
    if (bits == NULL || num_bits == 0) {
        return 0;
    }
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
    // As md_mm256_deperiodize2_ps: d from a fused multiply-subtract, see the deperiodize comment in md_simd.h
    md_256 dx = md_mm256_fmsub_ps(x, inv_period, md_mm256_mul_ps(r, inv_period));
    dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
    md_256 x_prim = md_mm256_fmadd_ps(dx, period, r);
    return md_mm256_select_ps(x_prim, x, md_mm256_cmpeq_ps(period, md_mm256_setzero_ps()));
}

// PACKED XYZ, COMPONENT WISE
//
// Eight packed atoms are 24 floats, and float p of them is component p % 3. An operation that treats
// x, y and z alike but apart - a min image in an orthorhombic cell, a blend, a sum - works on those
// floats as they are loaded, with its per axis constants laid out in the same pattern, and never
// splits them into x, y and z. The 24 floats are XYZ_GROUP_VECS native vectors (md_xv): three of
// eight lanes with AVX, six of four on NEON.
#define XYZ_GROUP_VECS (24 / MD_XV_WIDTH)

// Per axis constants in the lane pattern of a group
static inline void xyz_lane_pattern(md_xv out[XYZ_GROUP_VECS], float x, float y, float z) {
    const float v[3] = { x, y, z };
    float p[24];
    for (int k = 0; k < 24; ++k) p[k] = v[k % 3];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) out[j] = md_xv_loadu_ps(p + j * MD_XV_WIDTH);
}

// Eight per atom values (weights) spread to the lanes of their atoms' three components
static inline void xyz_lane_spread(md_xv out[XYZ_GROUP_VECS], const float* w8) {
#if MD_XV_WIDTH == 8
    const md_256 w = md_mm256_loadu_ps(w8);
    out[0] = md_mm256_permutevar8x32_ps(w, simde_mm256_setr_epi32(0, 0, 0, 1, 1, 1, 2, 2));
    out[1] = md_mm256_permutevar8x32_ps(w, simde_mm256_setr_epi32(2, 3, 3, 3, 4, 4, 4, 5));
    out[2] = md_mm256_permutevar8x32_ps(w, simde_mm256_setr_epi32(5, 5, 6, 6, 6, 7, 7, 7));
#else
    for (int h = 0; h < 2; ++h) {
        const md_128 w = md_mm_loadu_ps(w8 + 4 * h);
        out[3 * h + 0] = md_mm_shuffle_ps(w, w, SIMDE_MM_SHUFFLE(1, 0, 0, 0));   // w0 w0 w0 w1
        out[3 * h + 1] = md_mm_shuffle_ps(w, w, SIMDE_MM_SHUFFLE(2, 2, 1, 1));   // w1 w1 w2 w2
        out[3 * h + 2] = md_mm_shuffle_ps(w, w, SIMDE_MM_SHUFFLE(3, 3, 3, 2));   // w2 w3 w3 w3
    }
#endif
}

// A group's lanes summed per component
static inline void xyz_lane_reduce_add(double out[3], const md_xv v[XYZ_GROUP_VECS]) {
    float t[24];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) md_xv_storeu_ps(t + j * MD_XV_WIDTH, v[j]);
    for (int p = 0; p < 24; ++p) out[p % 3] += t[p];
}

// x brought to the image nearest r, lane by lane; a lane with no period (zero) is left alone. As
// md_mm256_deperiodize2_ps: d from a fused multiply-subtract, see the deperiodize comment in md_simd.h.
static inline md_xv xv_deperiodize_ortho(md_xv x, md_xv r, md_xv period, md_xv inv_period) {
    md_xv d = md_xv_fmsub_ps(x, inv_period, md_xv_mul_ps(r, inv_period));
    d = md_xv_sub_ps(d, md_xv_round_ps(d));
    const md_xv p = md_xv_fmadd_ps(d, period, r);
    return md_xv_select_ps(p, x, md_xv_cmpeq_ps(period, md_xv_setzero_ps()));
}

// Atoms to step over before packed coordinates at p are aligned to the given boundary, so every load
// of a group stays within its cache line; from a 16 byte boundary, which is where an md_array's
// elements begin, half of the 32 byte loads would straddle two. Zero if p is not float aligned.
static inline size_t xyz_align_peel_to(const void* p, size_t count, size_t alignment) {
    const size_t a = (size_t)((uintptr_t)p & (alignment - 1));
    if (a & 3) return 0;
    for (size_t k = 0; k < alignment / 4; ++k) {
        if (((a + 12 * k) & (alignment - 1)) == 0) return MIN(k, count);
    }
    return 0;
}

static inline size_t xyz_align_peel(const void* p, size_t count) {
    return xyz_align_peel_to(p, count, 32);
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

void md_util_system_extract_xyzw_from_mask(vec4_t* out_xyzw, const md_bitset_t* mask, const md_system_t* sys, const md_system_state_t* state) {
    ASSERT(out_xyzw);
    ASSERT(mask);
    ASSERT(sys);

    md_bitset_iter_t it = md_bitset_iter(*mask);
    size_t count = 0;
    while (md_bitset_iter_next(&it)) {
        size_t i = it.idx;
        vec3_t xyz = md_state_coord(state, i);
        float mass = md_atom_mass(&sys->atom, i);
        out_xyzw[count] = vec4_from_vec3(xyz, mass);
        count++;
    }
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

static bool md_util_amino_acid_atoms_extract(md_amino_acid_atoms_t* out_atoms, const md_atom_data_t* atom_data, md_urange_t atom_range) {
    ASSERT(atom_data);
    int len = atom_range.end - atom_range.beg;
    if (len < 4) return false;

    static const uint32_t all_bits = 1 | 2 | 4 | 8 | 16 | 32;
    uint32_t bits = 0;
    md_amino_acid_atoms_t aa;
	MEMSET(&aa, -1, sizeof(md_amino_acid_atoms_t));
    for (uint32_t i = atom_range.beg; i < atom_range.end; ++i) {
        str_t id = md_atom_name(atom_data, i);
        if (str_empty(id)) continue;

        if (!(bits & 1)  && cmp(id, "N",  1)) { aa.n  = (md_atom_idx_t)i; bits |= 1;  continue; }
        if (!(bits & 16) && cmp(id, "HN", 2)) { aa.hn = (md_atom_idx_t)i; bits |= 16; continue; }
        if (!(bits & 2)  && cmp(id, "CA", 2)) { aa.ca = (md_atom_idx_t)i; bits |= 2;  continue; }
        if (!(bits & 4)  && cmp(id, "C",  1)) { aa.c  = (md_atom_idx_t)i; bits |= 4;  continue; }
        if (!(bits & 32) && cmp(id, "CB", 2)) { aa.cb = (md_atom_idx_t)i; bits |= 32; continue; }
        if (!(bits & 8)) {
            if (cmp(id, "O", 1)     ||
                cmp(id, "O1", 2)    ||
                cmp(id, "OT1", 3)   ||
                cmp(id, "OC1", 3))
            {
                aa.o = (md_atom_idx_t)i; bits |= 8; continue;
            }
        }

        // Check if done
        if (bits == all_bits) break;
    }

    // HN is optional
    static const uint32_t backbone_bits = 1 | 2 | 4 | 8;

    if (out_atoms) *out_atoms = aa;
    return (bits & backbone_bits) == backbone_bits;
}

// There are several conventions to label Nucleic backbone atoms (e.g. C3', C3* or just C3). How the name matches: 0 not
// at all, 1 the bare name (C3), 2 the primed name (C3' or C3*).
static int cmp_nuc_atom(str_t id, const char* base, size_t len) {
    if (id.len < len || !str_eq_cstr_n(id, base, len)) return 0;
    if (id.len == len) return 1;
    const char next = str_ptr(id)[len];
    return (id.len == len + 1 && (next == '\'' || next == '*')) ? 2 : 0;
}

// Whether the atom whose name matches (cmp_nuc_atom) takes the slot (bit) of the sugar or backbone: a free slot takes any
// match, and a primed name takes the slot from a bare one. A residue which primes the atoms of its sugar has base atoms
// of the bare names (C2, C4, C5, O4 of a pyrimidine), and some list the base first: OMU and PSU of the PDB (3G71).
static inline bool nuc_atom_take_slot(uint32_t* bits, uint32_t* primed, uint32_t bit, int match) {
    if (!match || (*primed & bit) || (match == 1 && (*bits & bit))) return false;
    *bits |= bit;
    if (match == 2) *primed |= bit;
    return true;
}

static bool md_util_nucleic_acid_atoms_extract(md_nucleic_acid_atoms_t* out_atoms, const md_atom_data_t* atom_data, md_urange_t atom_range) {
    ASSERT(atom_data);
    int len = atom_range.end - atom_range.beg;
    if (len < 6) return false;

    static const uint32_t all_bits = 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 256;

    uint32_t bits = 0;
    uint32_t primed = 0; // The slots taken by a primed name, which no other atom takes back
    md_nucleic_acid_atoms_t na;
	MEMSET(&na, -1, sizeof(md_nucleic_acid_atoms_t));
    for (uint32_t i = atom_range.beg; i < atom_range.end; ++i) {
        str_t id = md_atom_name(atom_data, i);
        if (str_empty(id)) continue;

        if (!(bits & 1)  && cmp(id, "P",   1))               { na.p  = i; bits |= 1;  continue; }
        if (nuc_atom_take_slot(&bits, &primed, 2,  cmp_nuc_atom(id, "O5", 2))) { na.o5 = i; continue; }
        if (nuc_atom_take_slot(&bits, &primed, 4,  cmp_nuc_atom(id, "C5", 2))) { na.c5 = i; continue; }
        if (nuc_atom_take_slot(&bits, &primed, 8,  cmp_nuc_atom(id, "C4", 2))) { na.c4 = i; continue; }
        if (nuc_atom_take_slot(&bits, &primed, 16, cmp_nuc_atom(id, "C3", 2))) { na.c3 = i; continue; }
        // O3P (exact) counts as a bare O3
        if (nuc_atom_take_slot(&bits, &primed, 32, cmp(id, "O3P", 3) ? 1 : cmp_nuc_atom(id, "O3", 2))) { na.o3 = i; continue; }

        if (nuc_atom_take_slot(&bits, &primed, 64,  cmp_nuc_atom(id, "C1", 2))) { na.c1 = i; continue; }
        if (nuc_atom_take_slot(&bits, &primed, 128, cmp_nuc_atom(id, "C2", 2))) { na.c2 = i; continue; }
        if (nuc_atom_take_slot(&bits, &primed, 256, cmp_nuc_atom(id, "O4", 2))) { na.o4 = i; continue; }

        // Done when the P is found and every other slot holds a primed name: a bare name could still be taken back
        if ((bits & 1) && primed == (all_bits & ~1u)) break;
    }

    static const uint32_t req_bits = 2 | 4 | 8 | 16 | 32; // No phosphor in first terminal residue (starts with O5)

    if (out_atoms) *out_atoms = na;
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

static inline float element_ionic_radius(md_element_t element) {
    return element < Num_Elements ? element_ionic_radii_u8[element] * 0.01f : 0;
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

    // The search below only looks one lattice vector in each direction, so first bring dx into the
    // zone it covers. Without this a separation of more than one box - any point of an unwrapped
    // trajectory which has diffused out of the cell - silently gets a wrong image.
    // Reduce along c, then b, then a: the box is lower triangular, so each step leaves the
    // components already reduced untouched.
    for (int i = 2; i >= 0; --i) {
        if (box[i][i] != 0.0f) {
            const float n = roundf(dx[i] / box[i][i]);
            if (n != 0.0f) {
                for (int j = i; j >= 0; --j) {
                    dx[j] -= n * box[i][j];
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

static inline void simd_minimum_image_triclinic(md_256 dx[3], const float box[3][3]) {
    // Into the zone the one lattice vector search covers first - see minimum_image_triclinic
    for (int i = 2; i >= 0; --i) {
        if (box[i][i] != 0.0f) {
            const md_256 n = md_mm256_round_ps(md_mm256_mul_ps(dx[i], md_mm256_set1_ps(1.0f / box[i][i])));
            for (int j = i; j >= 0; --j) {
                dx[j] = md_mm256_sub_ps(dx[j], md_mm256_mul_ps(n, md_mm256_set1_ps(box[i][j])));
            }
        }
    }

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
                dsq_min = md_mm256_select_ps(dsq_min, dsq, mask);
                dx_min[0] = md_mm256_select_ps(dx_min[0], rz[0], mask);
                dx_min[1] = md_mm256_select_ps(dx_min[1], rz[1], mask);
                dx_min[2] = md_mm256_select_ps(dx_min[2], rz[2], mask);
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

bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const vec3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {
    if (!backbone_angles) return false;
    if (capacity == 0) return false;

    MEMSET(backbone_angles, 0, sizeof(md_backbone_angles_t) * capacity);

    if (!backbone) return false;
    if (!xyz) return false;
    if (!backbone->segment.atoms) return false;

    for (size_t bb_idx = 0; bb_idx < backbone->range.count; ++bb_idx) {
        const md_urange_t range = {backbone->range.offset[bb_idx], backbone->range.offset[bb_idx + 1]};
        ASSERT(range.end <= capacity);

        // A trace has one bead per residue, and so no phi and psi
        if (range.end - range.beg < 4 || md_protein_backbone_range_is_trace(backbone, bb_idx)) {
            continue;
        }

        for (size_t i = range.beg + 1; i + 1 < range.end; ++i) {
			int c_prev_idx = backbone->segment.atoms[i - 1].c;
			int ca_idx     = backbone->segment.atoms[i].ca;
			int  c_idx     = backbone->segment.atoms[i].c;
			int  n_idx     = backbone->segment.atoms[i].n;
			int n_next_idx = backbone->segment.atoms[i + 1].n;

			vec3_t c_prev = xyz[c_prev_idx];
            vec3_t n      = xyz[n_idx];
            vec3_t ca     = xyz[ca_idx];
            vec3_t c      = xyz[c_idx];
            vec3_t n_next = xyz[n_next_idx];

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

// ### THE BACKBONE OF A STATE ###

STATIC_ASSERT(sizeof(md_secondary_structure_t) == sizeof(int32_t), "The secondary structure is stored as I32");
STATIC_ASSERT(sizeof(md_backbone_angles_t) == 2 * sizeof(float), "The backbone angles are stored as F32 x 2");

static const md_attribute_t* state_backbone_attr(const md_system_state_t* state, str_t path, md_attribute_type_t type, uint32_t components, size_t num_segments) {
    if (!state || num_segments == 0) return NULL;
    const md_attribute_t* attr = md_attributes_find(&state->attributes, path);
    if (!md_attribute_view(attr, type, components, 1) || attr->format.shape[0] != num_segments) return NULL;
    return attr;
}

static void* state_backbone_write(md_system_state_t* state, const md_system_t* sys, str_t path, md_attribute_type_t type, uint32_t components, md_unit_t unit, str_t label) {
    if (!state || !sys) return NULL;
    const size_t num_segments = sys->protein_backbone.segment.count;
    if (num_segments == 0) return NULL;
    if (!state->attributes.alloc) {
        // The table of a state allocates with the state's allocator (see md_system_state_t); a view has neither
        if (!state->alloc) return NULL;
        state->attributes.alloc = state->alloc;
    }

    const md_attribute_t* attr = state_backbone_attr(state, path, type, components, num_segments);
    md_attribute_id_t id = attr ? attr->id : MD_ATTRIBUTE_INVALID;
    if (id == MD_ATTRIBUTE_INVALID) {
        const md_attribute_desc_t desc = {
            .path   = path,
            .format = { .type = type, .components = components, .rank = 1, .shape = { (uint32_t)num_segments } },
            .unit   = unit,
            .label  = label,
        };
        id = md_attributes_replace(&state->attributes, &desc);
        if (id == MD_ATTRIBUTE_INVALID) return NULL;
    }
    md_attributes_touch(&state->attributes, id);
    return md_attributes_data(&state->attributes, id, type);
}

const md_backbone_angles_t* md_util_state_backbone_angles(const md_system_state_t* state, const md_system_t* sys) {
    if (!sys) return NULL;
    const md_attribute_t* attr = state_backbone_attr(state, STR_LIT(MD_BACKBONE_ANGLE_PATH), MD_ATTRIBUTE_TYPE_F32, 2, sys->protein_backbone.segment.count);
    return attr ? (const md_backbone_angles_t*)attr->data : NULL;
}

const md_secondary_structure_t* md_util_state_secondary_structure(const md_system_state_t* state, const md_system_t* sys) {
    if (!sys) return NULL;
    const md_attribute_t* attr = state_backbone_attr(state, STR_LIT(MD_BACKBONE_SECONDARY_STRUCTURE_PATH), MD_ATTRIBUTE_TYPE_I32, 1, sys->protein_backbone.segment.count);
    return attr ? (const md_secondary_structure_t*)attr->data : NULL;
}

md_backbone_angles_t* md_util_state_backbone_angles_write(md_system_state_t* state, const md_system_t* sys) {
    return (md_backbone_angles_t*)state_backbone_write(state, sys, STR_LIT(MD_BACKBONE_ANGLE_PATH), MD_ATTRIBUTE_TYPE_F32, 2, md_unit_radian(), STR_LIT("Backbone Angles"));
}

md_secondary_structure_t* md_util_state_secondary_structure_write(md_system_state_t* state, const md_system_t* sys) {
    return (md_secondary_structure_t*)state_backbone_write(state, sys, STR_LIT(MD_BACKBONE_SECONDARY_STRUCTURE_PATH), MD_ATTRIBUTE_TYPE_I32, 1, md_unit_none(), STR_LIT("Secondary Structure"));
}

bool md_util_state_backbone_compute(md_system_state_t* state, const md_system_t* sys) {
    if (!state || !sys || !md_system_state_has_coords(state)) return false;
    const size_t num_segments = sys->protein_backbone.segment.count;
    if (num_segments == 0) return false;
    md_backbone_angles_t*     angle = md_util_state_backbone_angles_write(state, sys);
    md_secondary_structure_t* ss    = md_util_state_secondary_structure_write(state, sys);
    if (!angle || !ss) return false;
    const bool ok_angle = md_util_backbone_angles_compute(angle, num_segments, state->xyz, &state->unitcell, &sys->protein_backbone);
    const bool ok_ss    = md_dssp_compute(ss, num_segments, state->xyz, &state->unitcell, &sys->protein_backbone);
    return ok_angle && ok_ss;
}

// Of segment i of an atomistic backbone: 'G' a glycine, 'P' a proline, 0 any other residue. By the name of its
// residue, or for a segment without one (BACKBONES FROM THE BONDS) by its atoms: a glycine has no CB, and the N of a
// proline closes the ring N-CA-CB-CG-CD: it is bonded to a CD (besides its CA and the C of the segment before it, if
// it has_prev) which is bonded to a CG which is bonded to the CB. The carbonyl of a cap (ACE) on an N terminus closes
// none.
static char ramachandran_residue(const md_system_t* sys, const md_protein_backbone_data_t* backbone, size_t i, bool has_prev) {
    const md_residue_idx_t res_idx = backbone->segment.res_idx[i];
    if (res_idx >= 0) {
        const str_t name = md_system_residue_name(sys, res_idx);
        return str_eq(name, STR_LIT("GLY")) ? 'G' : str_eq(name, STR_LIT("PRO")) ? 'P' : 0;
    }
    const md_amino_acid_atoms_t* atoms = &backbone->segment.atoms[i];
    const md_atom_idx_t prev_c = has_prev ? backbone->segment.atoms[i - 1].c : -1;
    if (atoms->cb < 0) return 'G';
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)atoms->n);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const md_atom_idx_t cd = md_bond_iter_atom_index(&it);
        if (cd == atoms->ca || cd == prev_c || md_atom_atomic_number(&sys->atom, (size_t)cd) != MD_Z_C) continue;
        md_bond_iter_t it_cd = md_bond_iter(&sys->bond, (size_t)cd);
        for (; md_bond_iter_has_next(&it_cd); md_bond_iter_next(&it_cd)) {
            const md_atom_idx_t cg = md_bond_iter_atom_index(&it_cd);
            if (cg != atoms->n && md_bond_find(&sys->bond, cg, atoms->cb) >= 0) return 'P';
        }
    }
    return 0;
}

bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const md_system_t* sys) {
    ASSERT(ramachandran_types);
    MEMSET(ramachandran_types, MD_RAMACHANDRAN_TYPE_UNKNOWN, sizeof(md_ramachandran_type_t) * capacity);

    if (capacity == 0) return false;

    const md_protein_backbone_data_t* backbone = &sys->protein_backbone;
    if (backbone->segment.count == 0) return false;
    if (backbone->range.count == 0 || !backbone->range.offset) return false;

    ASSERT(backbone->segment.res_idx);

    // The types are MolProbity's (mmtbx.validation.ramalyze), decided in its order: glycine, proline, then pre-proline
    // for any other residue followed by a proline. A glycine or a proline before a proline keeps its own type.
    //
    // Walked range by range, because a residue is pre-proline only when the proline follows it in the same chain. The
    // last segment of a range has no successor; segment i + 1 is then the first residue of another chain.
    for (size_t range_idx = 0; range_idx < backbone->range.count; ++range_idx) {
        // A trace has no phi and psi to plot: its segments stay UNKNOWN
        if (md_protein_backbone_range_is_trace(backbone, range_idx)) continue;
        const size_t range_beg = backbone->range.offset[range_idx];
        const size_t range_end = backbone->range.offset[range_idx + 1];
        // Written up to the capacity, but the successor is read from the system, so the last written type is still
        // decided by the segment after it.
        const size_t end = MIN(range_end, capacity);

        for (size_t i = range_beg; i < end; ++i) {
            const char res = ramachandran_residue(sys, backbone, i, i > range_beg);
            if (res == 'G') {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GLYCINE;
            } else if (res == 'P') {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PROLINE;
            } else if (i + 1 < range_end && ramachandran_residue(sys, backbone, i + 1, true) == 'P') {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PREPROL;
            } else {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GENERAL;
            }
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

// Functional group and bond order perception, disabled. It was written against the graph matcher which md_match.c
// replaced, and md_chem.c does its job now. Kept for reference.
#if 0
static size_t find_isomorphisms(md_index_data_t*, const graph_t* ,const graph_t* ,md_util_match_mode_t ,int , md_allocator_i*);
static void find_isomorphisms_callback(const graph_t*, const graph_t*, int, state_t*);

static graph_t smiles_to_graph(str_t smiles_str, md_util_match_flags_t flags, md_allocator_i* arena) {
    typedef struct vertex_t {
        int type;
        int edge_count;
        int edge_idx[8];
        int edge_type[8];
    } vertex_t;

    md_temp_scope_t temp = md_temp_begin_avoid(arena);

    size_t node_cap = str_len(smiles_str);
    md_smiles_node_t* nodes = md_temp_alloc_array(temp, md_smiles_node_t, node_cap);
    const size_t num_nodes = md_smiles_parse(nodes, node_cap, str_ptr(smiles_str), str_len(smiles_str));
    vertex_t* verts = md_temp_alloc_array(temp, vertex_t, num_nodes * 4);
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
    md_temp_end(temp);
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
        md_bond_flags_t flags = (md_bond_flags_t)graph_edge_type(data->pattern, i);
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
    md_temp_scope_t temp = md_temp_begin_in(temp_arena);

    int* depth = md_temp_alloc_array(temp, int, graph->vertex_count);
    md_fifo_t fifo = md_fifo_create(32, temp_arena);

    md_fifo_push(&fifo, start_idx);
    while (!md_fifo_empty(&fifo)) {
        int idx = md_fifo_pop(&fifo);
        int d = depth[idx] + 1;
        graph_edge_iter_t it = graph_edge_iter(graph, idx);
        while (graph_edge_iter_has_next(it)) {
            int other_idx = graph_edge_iter_vidx(it);
            // @NOTE: We do not care about the minimum depth, for each node
            if (other_idx != start_idx && depth[other_idx] == 0) {
                depth[other_idx] = d;
                md_fifo_push(&fifo, other_idx);
            }
            graph_edge_iter_next(&it);
        }
    }

    int max_depth = 0;
    for (size_t i = 0; i < graph->vertex_count; ++i) {
        max_depth = MAX(max_depth, depth[i]);
    }

    md_temp_end(temp);

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
        str_t res_a = LBL_TO_STR(resname[atom_a]);
        str_t res_b = LBL_TO_STR(resname[atom_b]);

        md_order_t order = 1;
        if (str_eq(res_a, res_b)) {
            if (str_cmp_lex(type_a, type_b) > 0) {
                str_swap(&type_a, &type_b);
            }
            // Intra
            if (str_eq_cstr(type_a, "C") && str_eq_cstr(type_b, "O") && md_util_resname_amino_acid(res_a)) {
                order = 2;
            } else {
                char buf[32];
                int len = build_key(buf, res_a, type_a, type_b);
                if (find_str_in_str_arr(NULL, (str_t){buf,len}, intra_bond_order_table, ARRAY_SIZE(intra_bond_order_table))) {
                    bond_order[i] = 2;
                }
            }
        } else {
            // Inter
            if ( (str_eq_cstr(res_a, "LYS") && str_eq_cstr(type_a, "CZ") && str_eq_cstr(res_b, "RET") && str_eq_cstr(type_b, "C15")) ||
                 (str_eq_cstr(res_b, "LYS") && str_eq_cstr(type_b, "CZ") && str_eq_cstr(res_a, "RET") && str_eq_cstr(type_a, "C15")) ){
                bond_order[i] = 2;
            }
        }
        bond_order[i] = order;
    }
#else
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    uint8_t* type = md_temp_alloc_zero_array(temp, uint8_t, atom->count);
    uint8_t* elements = md_temp_alloc_array(temp, uint8_t, atom->count);

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


    md_tick_t t0 = md_tick_now();

    // 7
    // Patterns of functional groups to look for as graphs
    // The patterns are not strictly matched on the bond order, only on the element (if present) and if its terminal
    // The pattern encodes bond order which is applied if the pattern can be matched

    str_t functional_group_patterns[] = {
        STR_INIT("[CH](=O)N(*)(*)"),
        STR_INIT("C(=O)([OH])(*)"),
        STR_INIT("C(=O)(O*)(*)"),
        STR_INIT("C(=O)(S*)(*)"),
        STR_INIT("C(=O)([NH]*)(*)"),
        STR_INIT("C(=S)(S*)(*)"),
        STR_INIT("C(=S)([NH]*)(*)"),
        STR_INIT("C(=N*)([NH]*)([NH]*)"),
        STR_INIT("C(=[NH])([NH2])(*)"),
        STR_INIT("[N+](=[N-])(=N*)"),
        STR_INIT("P(=O)(*)(*)(*)"),
        STR_INIT("S(=O)(=O)(*)(*)"),
        STR_INIT("N(=O)(=O)(*)"),
        STR_INIT("[Se](=O)([OH])(*)"),
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
        graph_t graph = smiles_to_graph(functional_group_patterns[i], 0, temp_arena);
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

        md_temp_scope_t temp = md_temp_begin_in(temp_arena);

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

        md_temp_end(temp);
    }

    md_tick_t t1 = md_tick_now();
    printf("Time for functional group perception: %.3f ms\n", md_tick_to_milliseconds(t1-t0));

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

        t0 = md_tick_now();
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
                    atom->flags[*it] |= MD_ATOM_FLAG_AROMATIC;
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
        t1 = md_tick_now();
        printf("Time for aromatic ring perception: %.3f ms\n", md_tick_to_milliseconds(t1-t0));
    }

    for (size_t i = 0; i < atom->count; ++i) {
        atom->flags[i] = md_atom_flags_set_hybridization(atom->flags[i], (md_hybridization_t)type[i]);
    }

    for (size_t i = 0; i < bond->count; ++i) {
        if (atom->flags[bond->pairs[i].idx[0]] & MD_ATOM_FLAG_AROMATIC &&
            atom->flags[bond->pairs[i].idx[1]] & MD_ATOM_FLAG_AROMATIC) {
            bond->order[i] |= MD_BOND_FLAG_AROMATIC;
        }
    }

    md_temp_end(temp_scope);
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
        md_256 v_dist  = md_mm256_sqrt_ps(v_dist2);  // There are alot of gathers here, so we can afford latency of sqrt

        float dist[8];
        md_mm256_storeu_ps(dist, v_dist);

        md_256 rad_i   = md_mm256_mask_i32gather_ps(v0, data->radius, vi_idx, vmask, 4);
        md_256 rad_j   = md_mm256_mask_i32gather_ps(v0, data->radius, vj_idx, vmask, 4);
        md_256 rad_sum = md_mm256_add_ps(rad_i, rad_j);

        md_256 r_min_cov = md_mm256_mul_ps(vk_min, rad_sum);
        md_256 r_max_cov = md_mm256_mul_ps(vk_max, rad_sum);

        md_256 cmp_min = md_mm256_cmpgt_ps(v_dist, r_min_cov);
        md_256 cmp_max = md_mm256_cmplt_ps(v_dist, r_max_cov);

        int mask = md_mm256_movemask_ps(md_mm256_and_ps(cmp_min, cmp_max));
        if (mask) {
            size_t new_size = md_array_size(*data->candidates) + popcnt32(mask);
            md_array_ensure(*data->candidates, new_size, data->alloc);

            while (mask) {
                int bit_idx = ctz32(mask);
                mask = mask & ~(1 << bit_idx);
                int idx = (int)i + bit_idx;
                bond_pair_t cp = {
                    .atom_i = i_idx[idx],
                    .atom_j = j_idx[idx],
                    .dist   = dist[bit_idx],
                };
                md_array_push_no_grow(*data->candidates, cp);
            }
        }
        vi = md_mm256_add_epi32(vi, v8);
    }
}

// A water molecule whose atoms beyond its O and two H are virtual sites (4 and 5 site models)
// One oxygen, two hydrogens and possibly virtual sites (4 and 5 site models): a water molecule, whatever it is called
static bool atoms_are_water(const md_system_t* sys, md_urange_t range) {
    uint32_t num_o = 0, num_h = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
        if (z == MD_Z_O) num_o += 1;
        else if (z == MD_Z_H) num_h += 1;
        else if (md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_VIRTUAL_SITE) return false;
    }
    return num_o == 1 && num_h == 2;
}

static int compare_atom_pair(const void* a, const void* b) {
    const md_atom_pair_t* pa = (const md_atom_pair_t*)a;
    const md_atom_pair_t* pb = (const md_atom_pair_t*)b;
    if (pa->idx[0] != pb->idx[0]) return pa->idx[0] < pb->idx[0] ? -1 : 1;
    if (pa->idx[1] != pb->idx[1]) return pa->idx[1] < pb->idx[1] ? -1 : 1;
    return 0;
}

// The most covalent bonds an atom of the element takes, 0 for no limit. Generous: the limit is there to catch
// atoms that clash (a hydrogen within bonding distance of two atoms), not to judge chemistry.
static inline int max_covalent_bonds(md_atomic_number_t z) {
    switch (z) {
    case H: case F:                     return 1;
    case O:                             return 3;   // Oxonium, within one molecule (see below)
    case B: case C: case N: case Si:    return 4;
    case Cl: case Br: case I:           return 4;   // Perchlorate and its kind
    case P: case S: case Se: case As:   return 6;
    default:                            return 0;
    }
}

// The most coordination bonds a metal takes
static inline int max_coordination(md_atomic_number_t z) {
    if (z == Be || z == Mg || z == Al || z == Zn || z == Ga) return 6;
    if ((z >= 57 && z <= 71) || (z >= 89 && z <= 103)) return 9;   // Lanthanides, actinides
    return 8;   // Seven coordinate Ru, eight coordinate Ca, ...
}

// A bond in excess of what its atoms can take. Removed in order: bonds between residues (residues, molecules)
// before bonds within one, and within each the bonds whose length is furthest from the sum of the radii (a clash
// as much as a stretch).
typedef struct bond_excess_t {
    uint32_t inter;
    float    deviation;
    uint32_t idx;
} bond_excess_t;

static int compare_bond_excess(const void* a, const void* b) {
    const bond_excess_t* x = (const bond_excess_t*)a;
    const bond_excess_t* y = (const bond_excess_t*)b;
    if (x->inter != y->inter) return x->inter > y->inter ? -1 : 1;
    if (x->deviation != y->deviation) return x->deviation > y->deviation ? -1 : 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

// Factors on the sum of the covalent radii (Cordero et al. 2008, which covers metal-ligand distances too)
static const float bond_k_min    = 0.65f;   // Shorter than this is a clash, not a bond
static const float bond_k_cov    = 1.15f;   // Covalent bond between non-metals
static const float bond_k_coord  = 1.30f;   // Upper end of the candidate search
static const float bond_k_ligand = 1.25f;   // Metal to a donor (N, O, S, P, halogen): 1.0-1.1 in well resolved structures
static const float bond_k_carbon = 1.15f;   // Transition metal to carbon (carbenes, carbonyls, cyclopentadienyl)
static const float bond_k_metal  = 1.00f;   // Metal-metal

// How two atoms are bonded from their distance: within factor times the sum of their covalent radii (0: never, whatever
// the distance), and whether the bond is a coordination bond. same_res: the atoms are of one residue, where a bond
// may be a little longer (large flexible ligands). water: one of them is of a water molecule.
typedef struct bond_rule_t {
    float factor;
    bool  coordinate;
} bond_rule_t;

static bond_rule_t bond_rule(int ei, int ej, bool same_res, bool water) {
    const bond_rule_t none = {0};
    const bool mi = is_metal(ei);
    const bool mj = is_metal(ej);
    if (!mi && !mj) {
        const bond_rule_t r = { same_res ? bond_k_cov * 1.05f : bond_k_cov, false };
        return r;
    }
    if (mi != mj) {
        const int metal_e = mi ? ei : ej;
        const int non_metal_e = mi ? ej : ei;
        // Alkali ions are free ions in solution: a coordination frozen from the reference frame would be
        // wrong in every other. Water around metals is left out for the same reason.
        if (is_alkali(metal_e) || water) return none;
        if (is_metal_donor(non_metal_e)) {
            const bond_rule_t r = { bond_k_ligand, true };
            return r;
        }
        if (non_metal_e == C && is_transition_metal(metal_e)) {
            const bond_rule_t r = { bond_k_carbon, true };
            return r;
        }
        return none;
    }
    if (is_transition_metal(ei) && is_transition_metal(ej)) {
        const bond_rule_t r = { bond_k_metal, false };
        return r;
    }
    return none;
}

// A halogen bonded to a non-metal (an organohalogen, is_halogen) is no ligand of a metal: its contact with one, a
// halogen lying on a metal surface, is no coordination bond. A halide is.

// ### METALLIC PHASES ###
//
// A metal in bulk, at a surface or in a particle is not bonded as a molecule is: each atom has as many neighbours as its
// lattice gives it, and a cutoff on the sum of covalent radii falls within their shell (2.72 Å for gold, against a
// first shell of 2.6-3.1 Å), which bonds a part of it that looks random and changes from frame to frame. The atoms of
// a metallic phase (MD_ATOM_FLAG_METALLIC) are instead bonded to every metal atom of their first coordination shell,
// by METALLIC bonds, which are not bonds of the graph (md_types.h), and to nothing else: what touches a surface moves
// across it, so its contacts are those of a frame (md_util_bond_delta) rather than of the topology.
//
// THE FIRST SHELL ends at the first minimum of the distances between the metal atoms, each over the sum of the
// covalent radii of its atoms: the gap between two shells of a crystal, or the first minimum of the radial distribution
// of a liquid or a strained lattice. It is the lowest point of their histogram beyond its first peak and up to 1.5
// times it, which holds the gap after the first shell in the lattices of the metals: their second shell is at 1.41 times
// the first in a close packed lattice (fcc, hcp), and at 1.15 in a body centred cubic one, whose two shells (8 and 6
// neighbours) merge into one when warm. Told from the system rather than given per element, as the covalent radii put
// the first shell of the metals at 0.94 (Fe) to 1.18 (Al) times their sum, and strain and temperature move and widen it.
// METALLIC_SHELL without enough pairs to tell.
//
// A METALLIC ATOM is a metal with at least METALLIC_MIN_NEIGHBORS metal atoms in its first shell, and at least three
// times as many as the non-metals within reach of a bond to it (a donor, N, O, S, P or a halogen, within bond_k_ligand;
// carbon within bond_k_carbon; whether the pair would be bonded or not, as an alkali to its anion is not): a cation of
// an ionic crystal or an oxide is not (NaCl and MgO: 12 metal neighbours to 6 anions, CaF2: 12 to 8, TiO2: 10 to 6),
// nor a metal of a molecule or a cluster (Fe4S4, FeMoco, Mn4CaO5: 3 or 4 metal neighbours at most), whereas an atom at
// the surface of a metal has 7 to 9 metal neighbours and an adsorbate or two. Then so is a metal atom in the first
// shell of a metallic atom with at least METALLIC_MIN_EDGE_NEIGHBORS metal neighbours, again at least three times its
// non-metals: the edges and corners of a slab or a particle, an adatom in a hollow site.
#define METALLIC_MIN_NEIGHBORS      6
#define METALLIC_MIN_EDGE_NEIGHBORS 3
#define METALLIC_SHELL              1.25f       // On the sum of the covalent radii, without enough pairs to tell
#define METALLIC_SHELL_MIN          1.10f
#define METALLIC_SHELL_MAX          1.60f
#define METALLIC_MIN_PAIRS          64          // Pairs within the shell for its end to be told from them
#define METALLIC_HIST_LO            0.80f       // The histogram of the pairs, on the sum of the radii
#define METALLIC_HIST_HI            1.90f
#define METALLIC_HIST_BIN           0.02f
#define METALLIC_HIST_BINS          55          // (HI - LO) / BIN

typedef struct metal_pair_search_t {
    const float* radius;
    float        q_max;                         // On the sum of the radii
    uint32_t*    hist;                          // [METALLIC_HIST_BINS], or NULL
    md_array(bond_pair_t)* pairs;               // Within q_max, or NULL
    md_allocator_i* alloc;
} metal_pair_search_t;

static void metal_pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    metal_pair_search_t* s = (metal_pair_search_t*)user_param;
    for (size_t k = 0; k < num_pairs; ++k) {
        const uint32_t i = i_idx[k];
        const uint32_t j = j_idx[k];
        const float sum_r = s->radius[i] + s->radius[j];
        if (i == j || sum_r <= 0.0f) continue;
        const float d = sqrtf(ij_dist2[k]);
        const float q = d / sum_r;
        if (s->hist && q >= METALLIC_HIST_LO && q < METALLIC_HIST_HI) {
            s->hist[(int)((q - METALLIC_HIST_LO) / METALLIC_HIST_BIN)] += 1;
        }
        if (s->pairs && q > bond_k_min && q <= s->q_max) {
            const bond_pair_t p = { (int)MIN(i, j), (int)MAX(i, j), d };
            md_array_push(*s->pairs, p, s->alloc);
        }
    }
}

// The end of the first shell, on the sum of the radii, from the histogram of the metal-metal pairs (THE FIRST SHELL)
static float metallic_shell(const uint32_t hist[METALLIC_HIST_BINS]) {
    float smooth[METALLIC_HIST_BINS];
    for (int b = 0; b < METALLIC_HIST_BINS; ++b) {
        const uint32_t l = b > 0 ? hist[b - 1] : 0;
        const uint32_t r = b + 1 < METALLIC_HIST_BINS ? hist[b + 1] : 0;
        smooth[b] = (float)(l + 2 * hist[b] + r) * 0.25f;
    }
    const float q_peak_max = 1.35f;
    int peak = 0;
    for (int b = 1; b < METALLIC_HIST_BINS && METALLIC_HIST_LO + (b + 0.5f) * METALLIC_HIST_BIN < q_peak_max; ++b) {
        if (smooth[b] > smooth[peak]) peak = b;
    }
    uint32_t in_shell = 0;
    for (int b = 0; b <= peak; ++b) in_shell += hist[b];
    if (smooth[peak] <= 0.0f || in_shell < METALLIC_MIN_PAIRS) return METALLIC_SHELL;

    // The lowest point beyond the peak up to 1.5 times it, the middle of the run of bins as low as it
    const float q_peak = METALLIC_HIST_LO + (peak + 0.5f) * METALLIC_HIST_BIN;
    const int last = MIN(METALLIC_HIST_BINS - 1, (int)((1.5f * q_peak - METALLIC_HIST_LO) / METALLIC_HIST_BIN));
    int lo = -1;
    for (int b = peak + 1; b <= last; ++b) {
        if (lo < 0 || smooth[b] < smooth[lo]) lo = b;
    }
    if (lo < 0) return METALLIC_SHELL;
    int hi = lo;
    while (hi + 1 <= last && smooth[hi + 1] == smooth[lo]) hi += 1;
    const float q = METALLIC_HIST_LO + (0.5f * (lo + hi) + 0.5f) * METALLIC_HIST_BIN;
    return CLAMP(q, METALLIC_SHELL_MIN, METALLIC_SHELL_MAX);
}

// The metallic atoms (out_metallic, one per atom) and the pairs of them within the first shell (out_pairs), see
// METALLIC PHASES. contacts: the number of non-metals within reach of each metal atom.
static void infer_metallic_phases(uint8_t* out_metallic, md_array(bond_pair_t)* out_pairs, const md_system_state_t* state, const md_system_t* sys,
                                  const md_atomic_number_t* atomic_nr, const float* radius, const uint32_t* contacts, md_allocator_i* temp_arena) {
    const size_t num_atoms = state->num_atoms;
    MEMSET(out_metallic, 0, num_atoms);

    md_array(int) metals = 0;
    float max_rad = 0.0f;
    for (size_t i = 0; i < num_atoms; ++i) {
        if (!is_metal(atomic_nr[i]) || md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_VIRTUAL_SITE) continue;
        md_array_push(metals, (int)i, temp_arena);
        max_rad = MAX(max_rad, radius[i]);
    }
    const size_t num_metals = md_array_size(metals);
    if (num_metals < METALLIC_MIN_NEIGHBORS + 1 || max_rad <= 0.0f) return;

    // The histogram of the pairs, then the pairs within the first shell
    uint32_t hist[METALLIC_HIST_BINS] = {0};
    md_coord_stream_t coords = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), metals, num_metals);
    {
        const double cutoff = 2.0 * max_rad * METALLIC_HIST_HI;
        md_spatial_acc_t acc = { .alloc = temp_arena };
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &coords, .cutoff = cutoff, .unitcell = &state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX });
        metal_pair_search_t search = { .radius = radius, .hist = hist };
        md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, metal_pair_callback, &search);
    }
    const float shell = metallic_shell(hist);

    md_array(bond_pair_t) pairs = 0;
    {
        const double cutoff = 2.0 * max_rad * shell;
        md_spatial_acc_t acc = { .alloc = temp_arena };
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &coords, .cutoff = cutoff, .unitcell = &state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX });
        metal_pair_search_t search = { .radius = radius, .q_max = shell, .pairs = &pairs, .alloc = temp_arena };
        md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, metal_pair_callback, &search);
    }
    const size_t num_pairs = md_array_size(pairs);

    uint32_t* neighbors = md_alloc(temp_arena, sizeof(uint32_t) * num_atoms);
    MEMSET(neighbors, 0, sizeof(uint32_t) * num_atoms);
    for (size_t k = 0; k < num_pairs; ++k) {
        neighbors[pairs[k].atom_i] += 1;
        neighbors[pairs[k].atom_j] += 1;
    }
    for (size_t m = 0; m < num_metals; ++m) {
        const int i = metals[m];
        out_metallic[i] = neighbors[i] >= METALLIC_MIN_NEIGHBORS && neighbors[i] >= 3 * contacts[i];
    }
    // The edges, corners and adatoms: in the shell of a metallic atom of the core
    uint8_t* edge = md_alloc(temp_arena, num_atoms);
    MEMSET(edge, 0, num_atoms);
    for (size_t k = 0; k < num_pairs; ++k) {
        const int i = pairs[k].atom_i;
        const int j = pairs[k].atom_j;
        if (out_metallic[i] && !out_metallic[j]) edge[j] = 1;
        if (out_metallic[j] && !out_metallic[i]) edge[i] = 1;
    }
    for (size_t m = 0; m < num_metals; ++m) {
        const int i = metals[m];
        if (edge[i] && neighbors[i] >= METALLIC_MIN_EDGE_NEIGHBORS && neighbors[i] >= 3 * contacts[i]) out_metallic[i] = 1;
    }

    for (size_t k = 0; k < num_pairs; ++k) {
        if (out_metallic[pairs[k].atom_i] && out_metallic[pairs[k].atom_j]) md_array_push(*out_pairs, pairs[k], temp_arena);
    }
}

// ### MARTINI ###

// The versions of Martini, which give beads of the same names different masses and sizes, and some residues different
// beads. Which one a system is from its residue and bead names tell (martini_system).
enum {
    MARTINI_2   = 1,    // martini22, elnedyn21, elnedyn22
    MARTINI_2P  = 2,    // Polarizable Martini 2: martini22p, elnedyn22p
    MARTINI_3   = 4,    // martini3001
    MARTINI_ANY = MARTINI_2 | MARTINI_2P | MARTINI_3,
};

// The radius of a regular Martini bead (sigma 4.7 A), as md_tpr gives it: half the distance of the minimum of the
// Lennard-Jones potential, 2^(1/6) sigma / 2
#define MARTINI_R   2.638f

typedef struct {
    const char*         name;
    float               mass;
    float               radius;     // From the sigma of the bead type with itself, as md_tpr gives it
    md_particle_kind_t  particle;   // A massless bead is a virtual site; one without Lennard-Jones a charge site
} martini_db_bead_t;

typedef struct {
    const char*         name;
    uint32_t            versions;
    md_residue_kind_t kind;
    uint16_t            bead_offset;    // Into martini_db_residue_beads
    uint8_t             num_beads;
    uint16_t            bond_offset;    // Into martini_db_residue_bonds, pairs of beads of the residue
    uint8_t             num_bonds;
} martini_db_residue_t;

typedef struct {
    const char* name;       // A residue of a polymer (polyply), and its bead bonded to
    const char* bead;
    const char* next_name;  // the bead of the residue after it
    const char* next_bead;
    uint32_t    versions;
} martini_db_link_t;

// The residues of the Martini force fields: proteins (vermouth), nucleic acids, the lipidomes, sterols, solvents, ions,
// small molecules, sugars and polymers (polyply) of Martini 2 and 3, read from their topology files by
// tools/gen_martini_db.py, and the links between the residues of the polymers.
#include "md_martini_db.inl"

static inline const martini_db_bead_t* martini_db_bead(const martini_db_residue_t* res, size_t k) {
    return &martini_db_beads[martini_db_residue_beads[res->bead_offset + k]];
}

// An amino acid or a nucleotide, of a chain whose residues vary (termini, modifications, missing beads)
static inline bool martini_db_biopolymer(const martini_db_residue_t* res) {
    return res->kind == MD_RESIDUE_KIND_AMINO_ACID || res->kind == MD_RESIDUE_KIND_NUCLEOTIDE;
}

// The residues of the name, [*beg, *end)
static void martini_db_find(str_t name, size_t* beg, size_t* end) {
    size_t lo = 0, hi = ARRAY_SIZE(martini_db_residues);
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (str_cmp_lex(str_from_cstr(martini_db_residues[mid].name), name) < 0) lo = mid + 1; else hi = mid;
    }
    size_t e = lo;
    while (e < ARRAY_SIZE(martini_db_residues) && str_eq_cstr(name, martini_db_residues[e].name)) ++e;
    *beg = lo;
    *end = e;
}

// Whether the labels are the beads of the residue, each once (in any order). map[i] is then the bead of label i.
static bool martini_db_match(const martini_db_residue_t* res, const str_t labels[], size_t num_labels, uint8_t map[]) {
    if (num_labels != res->num_beads) return false;
    bool used[256] = {0};
    for (size_t i = 0; i < num_labels; ++i) {
        size_t k = 0;
        while (k < res->num_beads && (used[k] || !str_eq_cstr(labels[i], martini_db_bead(res, k)->name))) ++k;
        if (k == res->num_beads) return false;
        used[k] = true;
        map[i] = (uint8_t)k;
    }
    return true;
}

// The molecule (not an amino acid or a nucleotide) a residue of the labels is, of the version if there is one of it,
// else of any version. NULL when it is none.
static const martini_db_residue_t* martini_db_match_molecule(uint32_t version, str_t res_name, const str_t labels[], size_t num_labels, uint8_t map[]) {
    if (num_labels == 0 || num_labels > 255) return NULL;
    size_t beg, end;
    martini_db_find(res_name, &beg, &end);
    const martini_db_residue_t* other = NULL;
    for (size_t k = beg; k < end; ++k) {
        const martini_db_residue_t* res = &martini_db_residues[k];
        if (martini_db_biopolymer(res) || !martini_db_match(res, labels, num_labels, map)) continue;
        if (res->versions & version) return res;
        if (!other) other = res;
    }
    if (other) martini_db_match(other, labels, num_labels, map);
    return other;
}

// The molecule (not an amino acid or a nucleotide) whose beads the labels are, whatever its name: a residue renamed, as the
// three letters of a PDB file make CHOL CHL, PSM SPM and PCER CER. Of the version where there is one of it, and with the
// masses (optional) of its beads. NULL for a single bead, which says too little, and where molecules of different kinds
// have those beads. The names of a molecule of the Martini force fields are those of its chemistry, so molecules of the
// same beads are the same molecule as far as anything here goes (its kind, its bonds, its anatomy as a lipid).
static const martini_db_residue_t* martini_db_match_beads(uint32_t version, const str_t labels[], const float masses[], size_t num_labels, uint8_t map[]) {
    if (num_labels < 2 || num_labels > 255) return NULL;
    const martini_db_residue_t* match = NULL;   // Of the version
    const martini_db_residue_t* other = NULL;   // Of any
    for (size_t r = 0; r < ARRAY_SIZE(martini_db_residues); ++r) {
        const martini_db_residue_t* res = &martini_db_residues[r];
        if (res->num_beads != num_labels || martini_db_biopolymer(res) || !martini_db_match(res, labels, num_labels, map)) continue;
        bool same = true;
        for (size_t i = 0; i < num_labels && masses && same; ++i) same = masses[i] == martini_db_bead(res, map[i])->mass;
        if (!same) continue;
        if (other && other->kind != res->kind) return NULL;
        if (!other) other = res;
        if (!match && (res->versions & version)) match = res;
    }
    const martini_db_residue_t* res = match ? match : other;
    if (res) martini_db_match(res, labels, num_labels, map);
    return res;
}

// A molecule of the force field (not an amino acid or a nucleotide) of the name whose beads are those of none of its
// topologies, as most of the labels (more than half) are beads of the name: a variant of the molecule whose topology is
// not in the database, or a residue edited. NULL for a name of no such molecule, for a single bead (which says too
// little), and where most labels are not beads of it (atoms named as an atomistic force field names them, of which C1 and
// C2 are beads of a Martini cholesterol too). Of several molecules of the name, the first.
static const martini_db_residue_t* martini_db_unknown_variant(str_t res_name, const str_t labels[], size_t num_labels) {
    if (num_labels < 2) return NULL;
    size_t beg, end;
    martini_db_find(res_name, &beg, &end);
    const martini_db_residue_t* first = NULL;
    for (size_t k = beg; k < end && !first; ++k) {
        if (!martini_db_biopolymer(&martini_db_residues[k])) first = &martini_db_residues[k];
    }
    if (!first) return NULL;
    size_t known = 0;
    for (size_t i = 0; i < num_labels; ++i) {
        bool found = false;
        for (size_t k = beg; k < end && !found; ++k) {
            const martini_db_residue_t* res = &martini_db_residues[k];
            if (martini_db_biopolymer(res)) continue;
            for (size_t b = 0; b < res->num_beads && !found; ++b) found = str_eq_cstr(labels[i], martini_db_bead(res, b)->name);
        }
        known += found;
    }
    return 2 * known > num_labels ? first : NULL;
}

// Whether the residue is a monomer of a polymer: a link names it
static bool martini_db_linked(str_t name) {
    for (size_t k = 0; k < ARRAY_SIZE(martini_db_links); ++k) {
        if (str_eq_cstr(name, martini_db_links[k].name) || str_eq_cstr(name, martini_db_links[k].next_name)) return true;
    }
    return false;
}

// Bead idx of an amino acid or a nucleotide of the labels: of the residue of the version whose beads they are, else the
// first bead of its name in a residue of the version (the residue may lack beads), else in one of any version
static const martini_db_bead_t* martini_db_biopolymer_bead(uint32_t version, str_t res_name, const str_t labels[], size_t num_labels, size_t idx) {
    size_t beg, end;
    martini_db_find(res_name, &beg, &end);
    if (num_labels <= 255) {
        uint8_t map[256];
        for (size_t r = beg; r < end; ++r) {
            const martini_db_residue_t* res = &martini_db_residues[r];
            if (martini_db_biopolymer(res) && (res->versions & version) && martini_db_match(res, labels, num_labels, map)) return martini_db_bead(res, map[idx]);
        }
    }
    const str_t label = labels[idx];
    const martini_db_bead_t* other = NULL;
    for (size_t r = beg; r < end; ++r) {
        const martini_db_residue_t* res = &martini_db_residues[r];
        if (!martini_db_biopolymer(res)) continue;
        for (size_t k = 0; k < res->num_beads; ++k) {
            const martini_db_bead_t* bead = martini_db_bead(res, k);
            if (!str_eq_cstr(label, bead->name)) continue;
            if (res->versions & version) return bead;
            if (!other) other = bead;
        }
    }
    return other;
}

typedef struct {
    uint32_t version;
    bool     martini;   // A Martini system: residues of several beads which are Martini's, Martini water or polymers
} martini_system_t;

// Which version of Martini a system is from. Each residue whose beads are those of a residue of the database in some
// versions only votes for those: Ala with SC1, Tyr with SC4 and Trp with SC5 for Martini 3, the charge dummies and the
// polarizable water for polarizable Martini 2, cholesterol of 8 beads for Martini 2, ... Martini 3 when nothing tells,
// the current one and martinize2's default, and Martini 2 over polarizable Martini 2 on a tie.
static martini_system_t martini_system(const md_system_t* sys, str_t (*label)(const void*, size_t), const void* data) {
    size_t votes[3] = {0};  // MARTINI_2, MARTINI_2P, MARTINI_3
    bool martini = false;
    str_t labels[256];
    uint8_t map[256];
    for (size_t res_idx = 0; res_idx < sys->residue.count; ++res_idx) {
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        const size_t n = range.end - range.beg;
        if (n == 0 || n > 255) continue;
        size_t beg, end;
        martini_db_find(md_system_residue_name(sys, res_idx), &beg, &end);
        if (beg == end) continue;
        for (size_t i = 0; i < n; ++i) labels[i] = label(data, range.beg + i);
        uint32_t versions = 0;
        for (size_t k = beg; k < end; ++k) {
            const martini_db_residue_t* res = &martini_db_residues[k];
            if (!martini_db_match(res, labels, n, map)) continue;
            versions |= res->versions;
            // No atomistic residue has the name and the atom names of a Martini residue of several beads, of Martini water
            // or of a monomer of a Martini polymer
            if (n > 1 || res->kind == MD_RESIDUE_KIND_WATER || martini_db_linked(str_from_cstr(res->name))) martini = true;
        }
        if (versions == 0 || versions == MARTINI_ANY) continue;
        if (versions & MARTINI_2)  votes[0] += 1;
        if (versions & MARTINI_2P) votes[1] += 1;
        if (versions & MARTINI_3)  votes[2] += 1;
    }
    martini_system_t result = { .version = MARTINI_3, .martini = martini };
    if (votes[0] > votes[2] && votes[0] >= votes[1]) result.version = MARTINI_2;
    else if (votes[1] > votes[2] && votes[1] > votes[0]) result.version = MARTINI_2P;
    return result;
}

// Lipids by name: those of the Martini lipidomes (whose names the atomistic force fields mostly share: POPC, DPPC, CHOL,
// ...) and the cholesterol of CHARMM
bool md_util_resname_lipid(str_t str) {
    size_t beg, end;
    martini_db_find(str, &beg, &end);
    for (size_t k = beg; k < end; ++k) {
        if (martini_db_residues[k].kind == MD_RESIDUE_KIND_LIPID) return true;
    }
    return str_eq_cstr(str, "CHL1");
}

static str_t label_of_array(const void* data, size_t i) {
    return ((const str_t*)data)[i];
}

static str_t label_of_atom(const void* data, size_t i) {
    return md_atom_name((const md_atom_data_t*)data, i);
}

// The type of the name, element, mass, size and particle kind, added when there is none. md_particle_type_find_or_add goes by
// name and element alone, which gives every bead of a name the mass and size of the first: a Martini SC1 is 36 in Ala
// and 54 in Leu (Martini 3), and Trp's SC3 is a massless virtual site.
static md_particle_type_idx_t find_or_add_particle_type(md_particle_type_data_t* type, str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_particle_type_flags_t flags, md_allocator_i* alloc) {
    for (size_t t = 1; t < type->count; ++t) {
        if (type->z[t] == z && type->mass[t] == mass && type->radius[t] == radius && md_particle_type_flags(type, t) == flags && str_eq(md_particle_type_name(type, t), name)) {
            return (md_particle_type_idx_t)t;
        }
    }
    return md_particle_type_add(type, name, (str_t){0}, z, mass, radius, color, flags, alloc);
}

// Coarse grained (Martini) bond inference. A residue whose beads are those of a residue of the Martini force fields
// (martini_db), with their masses, has the bonds of its topology: those within the residue for an amino acid or a
// nucleotide, all of them for a molecule (a lipid, a sterol, a solvent, a sugar, ...). Amino acids and nucleotides which
// are not (missing beads, modified residues, other force fields) are bonded by their geometry:
//
// Within a residue
// - A minimum spanning tree of its beads, of bonds up to 6 A: the beads hang together the way their topologies have
//   them, and a stray bead stays unbonded.
// - The rings: side chain beads closer than 3.6 A (base beads closer than 3.0 A) besides those the tree bonds. Longer
//   ring bonds (Trp SC2-SC4 of Martini 3, 4.7 A; purine SC2-SC4 of Martini 2 DNA) are left out, as beads of the rings
//   that are not bonded come as close as 3.3 A.
// - A massless bead is a virtual site, which the topologies do not bond (Trp SC3 of Martini 3), and stays out, unless
//   it is all of the side chain there is (the polar SC1 of polarizable Martini 2, bonded to BB in its topology).
// - The charge dummies of polarizable Martini 2 (SCP, SCN) are bonded to the bead they sit by, and to each other.
// Between residues (of any amino acid or nucleotide)
// - The last backbone bead of a residue to the first of the next (BB to BB; BB3 to BB1 in DNA and RNA, and BB2 to BB1
//   as well in Martini 3 RNA), up to 5.5 A when the next is numbered next (bonds stretch in a simulation), else up to
//   4.3 A for amino acids (bonded pairs of structures reach 4.05 A) and 5.5 A for nucleotides (4.4 A), unless the
//   numbering starts over: a new chain. The numbering may jump forward over a bond (an engineered insert), so a gap is
//   told by the distance alone.
// - The SC1 beads of two cysteines closer than 3 A: a disulfide bridge (2.3-2.5 A; unbonded pairs are 4.5 A apart).
//
// The residues of a polymer (polyply's library: PEO, PS, PE, PMMA, dextran, ..., a PEGylated lipid) are bonded to the
// residue before them by the links of the force field.
// A virtual site is bonded to the first bead it is constructed from in the topologies, as md_tpr bonds it.
// Validated against the topologies martinize2 (vermouth 0.15) writes for 16 structures (dipeptides to 300 residues,
// several chains, disulfides, a chain with gaps) in Martini 3, 2 and polarizable 2, against Martini 2 DNA (polyply)
// mapped from a nucleosome, and against simulated systems of proteins, DNA, RNA, membranes and solvents in Martini 2
// and 3 (their .tpr): no bond which the topologies do not have, besides elastic networks and the restraints of
// secondary structure, which are not chemical bonds.
static inline float cg_bead_distance(const md_system_state_t* state, int i, int j) {
    vec3_t d = vec3_sub(state->xyz[i], state->xyz[j]);
    md_util_min_image_vec3(&d, 1, &state->unitcell);
    return vec3_length(d);
}

static inline void cg_push_bond(md_array(bond_pair_t)* candidates, int i, int j, float dist, md_allocator_i* alloc) {
    bond_pair_t cp = { .atom_i = i, .atom_j = j, .dist = dist };
    md_array_push(*candidates, cp, alloc);
}

static inline bool martini_charge_dummy(str_t atom_name) {
    return str_eq_cstr(atom_name, "SCP") || str_eq_cstr(atom_name, "SCN");
}

// The residue of the database a residue is: the one of its name whose beads are its beads, each with the mass of its
// atom (its type was given by it), and of the version where there are several. map[i] is then the bead of atom i.
static const martini_db_residue_t* martini_db_match_residue(const md_system_t* sys, size_t res_idx, uint32_t version, str_t labels[], uint8_t map[]) {
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    const size_t n = range.end - range.beg;
    if (n == 0 || n > 255) return NULL;
    const md_residue_kind_t kind = md_system_residue_kind(sys, res_idx);
    const bool biopolymer = kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
    for (size_t i = 0; i < n; ++i) labels[i] = md_atom_name(&sys->atom, range.beg + i);
    size_t beg, end;
    martini_db_find(md_system_residue_name(sys, res_idx), &beg, &end);
    const martini_db_residue_t* other = NULL;
    for (size_t r = beg; r < end; ++r) {
        const martini_db_residue_t* res = &martini_db_residues[r];
        if (martini_db_biopolymer(res) != biopolymer || !martini_db_match(res, labels, n, map)) continue;
        bool same = true;
        for (size_t i = 0; i < n && same; ++i) same = md_atom_mass(&sys->atom, range.beg + i) == martini_db_bead(res, map[i])->mass;
        if (!same) continue;
        if (res->versions & version) return res;
        if (!other) other = res;
    }
    if (other) {
        martini_db_match(other, labels, n, map);
        return other;
    }
    // A molecule renamed: by its beads and their masses
    if (biopolymer) return NULL;
    float masses[256];
    for (size_t i = 0; i < n; ++i) masses[i] = md_atom_mass(&sys->atom, range.beg + i);
    return martini_db_match_beads(version, labels, masses, n, map);
}

// The bonds of a residue of a polymer (polyply) to the one before it, by the links of the force field, up to 6 A apart
// (the longest are 4.7 A), unless the numbering starts over: a new chain
static void cg_polymer_links(md_array(bond_pair_t)* candidates, const md_system_t* sys, const md_system_state_t* state, size_t ci, md_allocator_i* alloc) {
    const float polymer_cap = 6.0f;
    const str_t prev_name = md_system_residue_name(sys, ci - 1);
    const str_t name      = md_system_residue_name(sys, ci);
    const md_sequence_id_t prev_seq = md_system_residue_seq_id(sys, ci - 1);
    const md_sequence_id_t seq      = md_system_residue_seq_id(sys, ci);
    if (seq < prev_seq && !(seq <= 1 && prev_seq >= 9999)) return;
    const md_urange_t prev = md_system_residue_atom_range(sys, ci - 1);
    const md_urange_t cur  = md_system_residue_atom_range(sys, ci);
    for (size_t k = 0; k < ARRAY_SIZE(martini_db_links); ++k) {
        const martini_db_link_t* link = &martini_db_links[k];
        if (!str_eq_cstr(prev_name, link->name) || !str_eq_cstr(name, link->next_name)) continue;
        int a = -1, b = -1;
        for (uint32_t i = prev.beg; i < prev.end && a == -1; ++i) if (str_eq_cstr(md_atom_name(&sys->atom, i), link->bead)) a = (int)i;
        for (uint32_t i = cur.beg;  i < cur.end  && b == -1; ++i) if (str_eq_cstr(md_atom_name(&sys->atom, i), link->next_bead)) b = (int)i;
        if (a == -1 || b == -1) continue;
        const float d = cg_bead_distance(state, a, b);
        if (d < polymer_cap) cg_push_bond(candidates, a, b, d, alloc);
    }
}

static void infer_coarse_grained_bonds(md_array(bond_pair_t)* candidates, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    const float tree_cap      = 6.0f;
    const float dummy_cap     = 3.0f;   // A charge dummy sits 1.1-1.4 A from its bead
    const float dummy_pair    = 3.2f;   // SCP-SCN are constrained at 2.8 A
    const float disulfide_cap = 3.0f;

    size_t max_len = 0;
    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_urange_t r = md_system_residue_atom_range(sys, ci);
        max_len = MAX(max_len, (size_t)(r.end - r.beg));
    }
    int*   tree   = md_alloc(alloc, sizeof(int)   * (max_len + 1));
    int*   parent = md_alloc(alloc, sizeof(int)   * (max_len + 1));
    float* best   = md_alloc(alloc, sizeof(float) * (max_len + 1));
    bool*  done   = md_alloc(alloc, sizeof(bool)  * (max_len + 1));
    md_array(int) cys_sc1 = 0;

    int prev_last_bb = -1, prev_also_bb = -1;
    size_t prev_ci = 0;
    md_residue_kind_t prev_kind = MD_RESIDUE_KIND_OTHER;

    str_t* labels = md_alloc(alloc, sizeof(str_t) * (max_len + 1));
    uint8_t map[256];
    const uint32_t version = martini_system(sys, label_of_atom, &sys->atom).version;

    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        const bool biopolymer = kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
        const md_urange_t range = md_system_residue_atom_range(sys, ci);

        if (ci > 0) cg_polymer_links(candidates, sys, state, ci, alloc);

        // The bonds of its topology. A charge site (a virtual site without Lennard-Jones) has none, as md_tpr gives it.
        const martini_db_residue_t* res = martini_db_match_residue(sys, ci, version, labels, map);
        if (res) {
            uint8_t local[256];
            for (size_t i = 0; i < range.end - range.beg; ++i) local[map[i]] = (uint8_t)i;
            for (size_t b = 0; b < res->num_bonds; ++b) {
                const int x = (int)range.beg + local[martini_db_residue_bonds[res->bond_offset + b][0]];
                const int y = (int)range.beg + local[martini_db_residue_bonds[res->bond_offset + b][1]];
                if (md_atom_particle_kind(&sys->atom, x) == MD_PARTICLE_VIRTUAL_SITE || md_atom_particle_kind(&sys->atom, y) == MD_PARTICLE_VIRTUAL_SITE) continue;
                cg_push_bond(candidates, x, y, cg_bead_distance(state, x, y), alloc);
            }
        }
        if (!biopolymer) {
            prev_last_bb = prev_also_bb = -1;
            continue;
        }
        const bool  amino_acid = kind == MD_RESIDUE_KIND_AMINO_ACID;
        const float ring_cap   = amino_acid ? 3.6f : 3.0f;
        const float link_cap      = amino_acid ? 4.3f : 5.5f;
        const float link_cap_next = 5.5f;

        // Without a topology, by the geometry. Which beads make up the tree: all but the charge dummies and, where the
        // side chain has beads with a mass, the massless ones
        bool side_with_mass = false;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) && !martini_charge_dummy(md_atom_name(&sys->atom, i)) && md_atom_mass(&sys->atom, i) > 0.0f) {
                side_with_mass = true;
            }
        }
        // The next residue is linked to the last bead of the backbone, and, in a nucleotide of Martini 3, to BB2 as well
        const bool m3_nucleotide = !amino_acid && ((res ? res->versions : version) & (MARTINI_2 | MARTINI_2P)) == 0;
        int n = 0;
        int first_bb = -1, last_bb = -1, also_bb = -1;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const bool backbone = md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE;
            if (backbone) {
                if (first_bb == -1) first_bb = (int)i;
                last_bb = (int)i;
                if (m3_nucleotide && str_eq_cstr(md_atom_name(&sys->atom, i), "BB2")) also_bb = (int)i;
            }
            if (res || martini_charge_dummy(md_atom_name(&sys->atom, i))) continue;
            if (!backbone && side_with_mass && md_atom_mass(&sys->atom, i) == 0.0f) continue;
            tree[n++] = (int)i;
        }

        // Prim's minimum spanning tree
        for (int k = 0; k < n; ++k) {
            best[k] = FLT_MAX;
            parent[k] = -1;
            done[k] = false;
        }
        if (n > 0) best[0] = 0.0f;
        for (int iter = 0; iter < n; ++iter) {
            int k = -1;
            for (int m = 0; m < n; ++m) {
                if (!done[m] && (k == -1 || best[m] < best[k])) k = m;
            }
            done[k] = true;
            if (parent[k] >= 0 && best[k] <= tree_cap) {
                cg_push_bond(candidates, tree[parent[k]], tree[k], best[k], alloc);
            } else {
                parent[k] = -1;
            }
            for (int m = 0; m < n; ++m) {
                if (done[m]) continue;
                const float d = cg_bead_distance(state, tree[k], tree[m]);
                if (d < best[m]) {
                    best[m] = d;
                    parent[m] = k;
                }
            }
        }

        // Ring bonds between side chain beads, besides the tree's
        for (int a = 0; a < n; ++a) {
            if (md_atom_flags(&sys->atom, tree[a]) & MD_ATOM_FLAG_BACKBONE) continue;
            for (int b = a + 1; b < n; ++b) {
                if (md_atom_flags(&sys->atom, tree[b]) & MD_ATOM_FLAG_BACKBONE) continue;
                if (parent[a] == b || parent[b] == a) continue;
                const float d = cg_bead_distance(state, tree[a], tree[b]);
                if (d < ring_cap) cg_push_bond(candidates, tree[a], tree[b], d, alloc);
            }
        }

        // The charge dummies
        for (uint32_t i = range.beg; i < range.end && !res; ++i) {
            if (!martini_charge_dummy(md_atom_name(&sys->atom, i))) continue;
            int nearest = -1;
            float nearest_d = dummy_cap;
            for (uint32_t j = range.beg; j < range.end; ++j) {
                if (j == i) continue;
                const float d = cg_bead_distance(state, (int)i, (int)j);
                if (martini_charge_dummy(md_atom_name(&sys->atom, j))) {
                    if (j > i && d < dummy_pair) cg_push_bond(candidates, (int)i, (int)j, d, alloc);
                } else if (d < nearest_d) {
                    nearest = (int)j;
                    nearest_d = d;
                }
            }
            if (nearest != -1) cg_push_bond(candidates, nearest, (int)i, nearest_d, alloc);
        }

        // To the previous residue of the chain. The next number allows for the stretch of a bond in a simulation
        // (Martini 2's backbone bonds, the softest, pass 4.3 A in one of 200 frames and 5.5 A in one of 10^5).
        if (prev_last_bb != -1 && first_bb != -1 && prev_ci + 1 == ci && prev_kind == kind) {
            const md_sequence_id_t prev_seq = md_system_residue_seq_id(sys, prev_ci);
            const md_sequence_id_t seq      = md_system_residue_seq_id(sys, ci);
            const bool restart = seq < prev_seq && !(seq <= 1 && prev_seq >= 9999);  // GRO and PDB wrap at 99999 and 9999
            const bool next    = seq == prev_seq + 1 || (seq <= 1 && prev_seq >= 9999);
            if (!restart) {
                const float cap = next ? link_cap_next : link_cap;
                const float d = cg_bead_distance(state, prev_last_bb, first_bb);
                if (d < cap) {
                    cg_push_bond(candidates, prev_last_bb, first_bb, d, alloc);
                    if (prev_also_bb != -1 && prev_also_bb != prev_last_bb) {
                        const float d2 = cg_bead_distance(state, prev_also_bb, first_bb);
                        if (d2 < cap) cg_push_bond(candidates, prev_also_bb, first_bb, d2, alloc);
                    }
                }
            }
        }
        prev_last_bb = last_bb;
        prev_also_bb = also_bb;
        prev_ci = ci;
        prev_kind = kind;

        const str_t res_name = md_system_residue_name(sys, ci);
        if (amino_acid && (str_eq_cstr(res_name, "CYS") || str_eq_cstr(res_name, "CYX"))) {
            for (uint32_t i = range.beg; i < range.end; ++i) {
                if (str_eq_cstr(md_atom_name(&sys->atom, i), "SC1")) md_array_push(cys_sc1, (int)i, alloc);
            }
        }
    }

    // Disulfide bridges
    const size_t num_cys = md_array_size(cys_sc1);
    for (size_t a = 0; a < num_cys; ++a) {
        for (size_t b = a + 1; b < num_cys; ++b) {
            const float d = cg_bead_distance(state, cys_sc1[a], cys_sc1[b]);
            if (d < disulfide_cap) cg_push_bond(candidates, cys_sc1[a], cys_sc1[b], d, alloc);
        }
    }
}

void md_util_infer_covalent_bonds(md_bond_data_t* bond, const md_system_state_t* state, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(bond);
    ASSERT(state);
    ASSERT(sys);
    ASSERT(alloc);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    // Inference only replaces what inference produced. Every bond of another origin (user defined, topology, or read
    // from the file) is kept, wherever it sits in the arrays, and put back afterwards.
    md_array(md_atom_pair_t)  bond_pairs = 0;
    md_array(md_bond_flags_t) bond_flags = 0;
    const md_bond_flags_t inferred = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_INFERRED);

    for (size_t i = 0; i < bond->count; ++i) {
        if (md_bond_origin(bond->flags[i]) != MD_BOND_ORIGIN_INFERRED) {
            md_array_push(bond_pairs, bond->pairs[i], temp_arena);
            md_array_push(bond_flags, bond->flags[i], temp_arena);
        }
    }

    // The non-covalent bonds are not inferred: they are kept as they are, after the graph (md_bond_data_t)
    const size_t num_noncovalent = bond->noncovalent_count;
    md_atom_pair_t*  noncovalent_pairs = md_temp_alloc_array(temp_scope, md_atom_pair_t,  MAX(num_noncovalent, 1));
    md_bond_flags_t* noncovalent_flags = md_temp_alloc_array(temp_scope, md_bond_flags_t, MAX(num_noncovalent, 1));
    if (num_noncovalent) {
        MEMCPY(noncovalent_pairs, bond->pairs + bond->count, sizeof(md_atom_pair_t)  * num_noncovalent);
        MEMCPY(noncovalent_flags, bond->flags + bond->count, sizeof(md_bond_flags_t) * num_noncovalent);
    }

    md_bond_data_clear(bond);
    
    if (!state->xyz) {
        MD_LOG_ERROR("Missing atom coordinates");
        goto done;
    }

    const bool is_coarse_grained = md_system_is_coarse_grained(sys);

    size_t num_atoms = state->num_atoms;
    md_array(bond_pair_t) candidates = 0;
    md_array_ensure(candidates, num_atoms, temp_arena);

    if (is_coarse_grained) {
        infer_coarse_grained_bonds(&candidates, sys, state, temp_arena);

        size_t num_candidates = md_array_size(candidates);
		md_array_ensure(bond->pairs, num_candidates, alloc);
        md_array_ensure(bond->flags, num_candidates, alloc);

        for (size_t i = 0; i < num_candidates; ++i) {
            md_atom_pair_t pair = {
                candidates[i].atom_i,
                candidates[i].atom_j,
            };
            md_array_push_no_grow(bond->pairs, pair);
            md_array_push_no_grow(bond->flags, inferred);
            bond->count += 1;
        }
	} else {
        double max_atom_rad = 0.0;

        float* atom_radius = md_temp_alloc_array(temp_scope, float, num_atoms);
        md_atomic_number_t* atomic_nr = md_temp_alloc_array(temp_scope, md_atomic_number_t, num_atoms);
        for (size_t i = 0; i < num_atoms; ++i) {
            atomic_nr[i]  = md_atom_atomic_number(&sys->atom, i);
            // The covalent radius for metals too: the ionic radius of the metal plus the covalent radius of the
            // ligand mixes two scales and comes out 0.4-0.6 Å short of real metal-ligand distances.
            atom_radius[i] = element_covalent_radius(atomic_nr[i]);
            max_atom_rad = MAX(max_atom_rad, atom_radius[i]);
        }

        int* comp_id = md_temp_alloc_array(temp_scope, int, num_atoms);
        MEMSET(comp_id, -1, sizeof(int) * num_atoms);
        size_t num_res = md_system_residue_count(sys);

        for (size_t i = 0; i < num_res; ++i) {
            md_urange_t range = md_system_residue_atom_range(sys, i);
            for (uint32_t j = range.beg; j < range.end; ++j) {
                comp_id[j] = (int)i;
            }
        }

        // The residues are usually not classified yet (that needs the bonds), so whether a residue is a water
        // molecule is told here from its atoms. Bond perception needs it around metals.
        bool* res_water = md_temp_alloc_array(temp_scope, bool, num_res);
        for (size_t i = 0; i < num_res; ++i) {
            res_water[i] = md_system_residue_kind(sys, i) == MD_RESIDUE_KIND_WATER || atoms_are_water(sys, md_system_residue_atom_range(sys, i));
        }

        {
            // Find connections for all atoms
            cov_bond_callback_param_t param = {
                .radius = atom_radius,
                .k_min = bond_k_min,
                .k_max = bond_k_coord,
                .candidates = &candidates,
                .alloc      = temp_arena,
            };

            // Candidates within reach of the largest covalent radius of the set: about 2.7 Å for organic matter
            const double cutoff = MAX(1.0, 2.0 * max_atom_rad * bond_k_coord);

            // Build candidate list
            md_coord_stream_t coords = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), NULL, state->num_atoms);
            md_spatial_acc_t acc = {.alloc = temp_arena};
            md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &coords, .cutoff = cutoff, .unitcell = &state->unitcell });
            md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, test_cov_bond_pair_callback, &param);
        }

        size_t num_candidates = md_array_size(candidates);

        //MD_LOG_DEBUG("Found %zu candidates", num_candidates);

        // The rule of each candidate, 0 for a pair which is not bonded: virtual sites sit within bonding distance of
        // their parent atom (the lone pairs of TIP5P at 0.7 Å) but are not bonded to anything
        bond_rule_t* rules = md_temp_alloc_zero_array(temp_scope, bond_rule_t, MAX(num_candidates, 1));
        uint8_t* organohalogen = md_temp_alloc_zero_array(temp_scope, uint8_t, num_atoms);
        for (size_t i = 0; i < num_candidates; ++i) {
            const int ai = candidates[i].atom_i;
            const int aj = candidates[i].atom_j;
            if (md_atom_particle_kind(&sys->atom, ai) == MD_PARTICLE_VIRTUAL_SITE || md_atom_particle_kind(&sys->atom, aj) == MD_PARTICLE_VIRTUAL_SITE) continue;
            const int ci = comp_id[ai];
            const int cj = comp_id[aj];
            const bool is_water = (ci >= 0 && res_water[ci]) || (cj >= 0 && res_water[cj]);
            const bond_rule_t rule = bond_rule(atomic_nr[ai], atomic_nr[aj], ci == cj, is_water);
            if (rule.factor <= 0.0f || candidates[i].dist >= rule.factor * (atom_radius[ai] + atom_radius[aj])) continue;
            rules[i] = rule;
            if (!rule.coordinate && !is_metal(atomic_nr[ai]) && !is_metal(atomic_nr[aj])) {
                if (is_halogen(atomic_nr[ai])) organohalogen[ai] = 1;
                if (is_halogen(atomic_nr[aj])) organohalogen[aj] = 1;
            }
        }
        // An organohalogen coordinates nothing
        for (size_t i = 0; i < num_candidates; ++i) {
            if (rules[i].coordinate && (organohalogen[candidates[i].atom_i] || organohalogen[candidates[i].atom_j])) {
                rules[i].factor = 0.0f;
            }
        }

        // The metallic phases, from the metal atoms and the non-metals within reach of each (METALLIC PHASES)
        uint32_t* contacts = md_temp_alloc_zero_array(temp_scope, uint32_t, num_atoms);
        for (size_t i = 0; i < num_candidates; ++i) {
            const int ai = candidates[i].atom_i;
            const int aj = candidates[i].atom_j;
            if (is_metal(atomic_nr[ai]) == is_metal(atomic_nr[aj])) continue;
            const int metal = is_metal(atomic_nr[ai]) ? ai : aj;
            const int other = metal == ai ? aj : ai;
            if (md_atom_particle_kind(&sys->atom, other) == MD_PARTICLE_VIRTUAL_SITE) continue;
            const int z = atomic_nr[other];
            const float factor = z == C ? bond_k_carbon : is_metal_donor(z) ? bond_k_ligand : 0.0f;
            if (candidates[i].dist < factor * (atom_radius[ai] + atom_radius[aj])) contacts[metal] += 1;
        }
        uint8_t* metallic = md_temp_alloc_array(temp_scope, uint8_t, num_atoms);
        md_array(bond_pair_t) metallic_pairs = 0;
        infer_metallic_phases(metallic, &metallic_pairs, state, sys, atomic_nr, atom_radius, contacts, temp_arena);

        md_bond_data_t temp_bond = { 0 };
        // We are lazy and propagate the correct periodic distance for each bond, so we can do more informed pruning later (e.g. for coordination bonds)
        md_array(float) temp_bond_dist = 0;

        for (size_t i = 0; i < num_candidates; ++i) {
            if (rules[i].factor <= 0.0f) continue;
            const int ai = candidates[i].atom_i;
            const int aj = candidates[i].atom_j;
            // An atom of a metallic phase has its metallic bonds, below, and no other: what touches a metal moves
            // across it, and has its contacts in a frame (md_util_bond_delta)
            if (metallic[ai] || metallic[aj]) continue;
            md_bond_insert(&temp_bond, ai, aj, rules[i].coordinate ? md_bond_flags_set_kind(inferred, MD_BOND_KIND_COORDINATE) : inferred, temp_arena);
            md_array_push(temp_bond_dist, candidates[i].dist, temp_arena);
        }
        {
            const md_bond_flags_t flags = md_bond_flags_set_kind(inferred, MD_BOND_KIND_METALLIC);
            for (size_t k = 0; k < md_array_size(metallic_pairs); ++k) {
                md_bond_insert(&temp_bond, metallic_pairs[k].atom_i, metallic_pairs[k].atom_j, flags, temp_arena);
                md_array_push(temp_bond_dist, metallic_pairs[k].dist, temp_arena);
            }
        }

        // Create bitfield to keep track of valid bonds during pruning (we do not remove bonds to ensure index stability)
        md_array(uint64_t) temp_bond_valid = make_bitfield(temp_bond.count, temp_arena);
        bitfield_set_all(temp_bond_valid, temp_bond.count); // Set all bits to 1 (valid)

        // ### PRUNE ###
        // Build temporary connectivity data for pruning
        md_bond_build_connectivity(&temp_bond, sys->atom.count, temp_arena);

        // Atoms with more bonds than they can take: atoms that clash (overlapping molecules in a packed starting
        // structure, unresolved alternate conformations), most often a hydrogen within reach of two atoms. Bonds go
        // in the order of compare_bond_excess until every atom is within its limit. Covalent bonds and coordination
        // bonds are limited separately.
        {
            const size_t num_temp = temp_bond.count;
            int* num_cov   = md_temp_alloc_zero_array(temp_scope, int, num_atoms);
            int* num_coord = md_temp_alloc_zero_array(temp_scope, int, num_atoms);
            int* max_cov   = md_temp_alloc_zero_array(temp_scope, int, num_atoms);
            for (size_t b = 0; b < num_temp; ++b) {
                // A metallic bond counts towards nothing: its atoms have as many as their lattice gives them
                if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_METALLIC) continue;
                const int x = temp_bond.pairs[b].idx[0];
                const int y = temp_bond.pairs[b].idx[1];
                int* n = md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_COORDINATE ? num_coord : num_cov;
                n[x] += 1;
                n[y] += 1;
                // A third bond on an oxygen only within its own molecule (oxonium); one reaching into another
                // molecule is a clash
                if (md_bond_kind(temp_bond.flags[b]) != MD_BOND_KIND_COORDINATE && comp_id[x] != comp_id[y]) {
                    if (atomic_nr[x] == O) max_cov[x] = 2;
                    if (atomic_nr[y] == O) max_cov[y] = 2;
                }
            }
            for (size_t i = 0; i < num_atoms; ++i) {
                if (!max_cov[i]) max_cov[i] = max_covalent_bonds(atomic_nr[i]);
            }

            md_array(bond_excess_t) over = 0;
            for (size_t b = 0; b < num_temp; ++b) {
                if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_METALLIC) continue;
                const int x = temp_bond.pairs[b].idx[0];
                const int y = temp_bond.pairs[b].idx[1];
                bool exceeds = false;
                if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_COORDINATE) {
                    const int metal = is_metal(atomic_nr[x]) ? x : y;
                    exceeds = num_coord[metal] > max_coordination(atomic_nr[metal]);
                } else {
                    const int mx = max_cov[x];
                    const int my = max_cov[y];
                    exceeds = (mx && num_cov[x] > mx) || (my && num_cov[y] > my);
                }
                if (exceeds) {
                    const float sum_r = atom_radius[x] + atom_radius[y];
                    bond_excess_t e = {
                        .inter = comp_id[x] != comp_id[y],
                        .deviation = sum_r > 0 ? fabsf(temp_bond_dist[b] / sum_r - 1.0f) : 0.0f,
                        .idx = (uint32_t)b,
                    };
                    md_array_push(over, e, temp_arena);
                }
            }

            const size_t num_over = md_array_size(over);
            if (num_over) {
                qsort(over, num_over, sizeof(bond_excess_t), compare_bond_excess);
                for (size_t k = 0; k < num_over; ++k) {
                    const uint32_t b = over[k].idx;
                    const int x = temp_bond.pairs[b].idx[0];
                    const int y = temp_bond.pairs[b].idx[1];
                    if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_COORDINATE) {
                        const int metal = is_metal(atomic_nr[x]) ? x : y;
                        if (num_coord[metal] <= max_coordination(atomic_nr[metal])) continue;
                        num_coord[x] -= 1;
                        num_coord[y] -= 1;
                    } else {
                        const int mx = max_cov[x];
                        const int my = max_cov[y];
                        if (!((mx && num_cov[x] > mx) || (my && num_cov[y] > my))) continue;
                        num_cov[x] -= 1;
                        num_cov[y] -= 1;
                    }
                    bitfield_clear_bit(temp_bond_valid, b);
                }
            }
        }

# if 0
        // @NOTE(Robin): If this is to be re-enabled, it needs to be verified as it did not produce exactly the expected results
        // without side effects in some dataset
        // 
        // Remove *triangles* i.e. if A-B and B-C are covalent, then A-C should not be covalent, (unless all carbon)
        // In such case remove the longest bond
        for (int bij = 0; bij < (int)temp_bond.count; ++bij) {
            if (!bitfield_test_bit(temp_bond_valid, bij)) {
                continue; // Skip invalidated bonds
            }

            int i = temp_bond.pairs[bij].idx[0];
            int j = temp_bond.pairs[bij].idx[1];
            float dij = temp_bond_dist[bij];

            // Check if atoms i and j have a common neighbor k
            md_bond_iter_t iter_i = md_bond_iter(&temp_bond, i);
            while (md_bond_iter_has_next(&iter_i)) {
                int k = md_bond_iter_atom_index(&iter_i);
                bool all_carbon = (atomic_nr[i] == 6 && atomic_nr[j] == 6 && atomic_nr[k] == 6);

                if (k != j && !all_carbon) {
                    int bik = md_bond_iter_bond_index(&iter_i);
                    float dik = temp_bond_dist[bik];

                    // Find other path between i and j, through k
                    int bjk = md_bond_find(&temp_bond, j, k);
                    if (bjk != -1) {
                        float djk = temp_bond_dist[bjk];
                    
                        // Compare all three distances and flag the longest bond for removal
                        if (dij >= dik && dij >= djk) {
                            bitfield_clear_bit(temp_bond_valid, bij);
                        } else if (dik >= dij && dik >= djk) {
                            bitfield_clear_bit(temp_bond_valid, bik);
                        } else {
                            bitfield_clear_bit(temp_bond_valid, bjk);
                        }
                    }
                }
                md_bond_iter_next(&iter_i);
            }
        }
#endif

        size_t num_valid_bonds = bitfield_popcount(temp_bond_valid, temp_bond.count);

        // Populate real bond data
        md_array_ensure(bond->pairs, num_valid_bonds, alloc);
        md_array_ensure(bond->flags, num_valid_bonds, alloc);

        for (size_t i = 0; i < temp_bond.count; ++i) {
            if (bitfield_test_bit(temp_bond_valid, i)) {
                md_bond_insert(bond, temp_bond.pairs[i].idx[0], temp_bond.pairs[i].idx[1], temp_bond.flags[i], alloc);
            }
        }
    }

    // Reconcile the fresh inferred bonds with the kept ones:
    // - an inferred bond naming a pair that is already kept is a duplicate, the kept one wins
    // - where a topology (md_itp) defines the connectivity it is authoritative, so inferred bonds between two
    //   atoms it covers are dropped, or re-inferring would add the distance guesses back on top of it
    // - a pair with a non-covalent bond keeps it: an atom pair has one bond, and that one was given
    {
        const size_t num_kept  = md_array_size(bond_pairs);
        const size_t num_saved = num_kept + num_noncovalent;
        if (num_saved > 0 && bond->count > 0) {
            const size_t num_atoms = state->num_atoms;
            uint64_t* covered = make_bitfield(num_atoms, temp_arena);

            // An elastic bond is no bond of the molecule: it neither stands in for an inferred bond nor says that a
            // topology covers its atoms
            md_atom_pair_t* saved_sorted = md_temp_alloc_array(temp_scope, md_atom_pair_t, num_saved);
            size_t num_sorted = 0;
            for (size_t i = 0; i < num_saved; ++i) {
                const md_bond_flags_t f = i < num_kept ? bond_flags[i] : noncovalent_flags[i - num_kept];
                if (md_bond_kind(f) == MD_BOND_KIND_ELASTIC) continue;
                md_atom_pair_t pair = i < num_kept ? bond_pairs[i] : noncovalent_pairs[i - num_kept];
                if (pair.idx[0] > pair.idx[1]) { md_atom_idx_t t = pair.idx[0]; pair.idx[0] = pair.idx[1]; pair.idx[1] = t; }
                saved_sorted[num_sorted++] = pair;
                if (i < num_kept && md_bond_origin(f) == MD_BOND_ORIGIN_TOPOLOGY) {
                    for (int k = 0; k < 2; ++k) {
                        if (pair.idx[k] >= 0 && (size_t)pair.idx[k] < num_atoms) bitfield_set_bit(covered, pair.idx[k]);
                    }
                }
            }
            qsort(saved_sorted, num_sorted, sizeof(md_atom_pair_t), compare_atom_pair);

            size_t w = 0;
            for (size_t i = 0; i < bond->count; ++i) {
                const md_atom_pair_t pair = bond->pairs[i];
                if (bitfield_test_bit(covered, pair.idx[0]) && bitfield_test_bit(covered, pair.idx[1])) continue;
                md_atom_pair_t key = pair;
                if (key.idx[0] > key.idx[1]) { md_atom_idx_t t = key.idx[0]; key.idx[0] = key.idx[1]; key.idx[1] = t; }
                if (num_sorted && bsearch(&key, saved_sorted, num_sorted, sizeof(md_atom_pair_t), compare_atom_pair)) continue;
                bond->pairs[w] = pair;
                bond->flags[w] = bond->flags[i];
                w += 1;
            }
            bond->count = w;
            md_array_shrink(bond->pairs, w);
            md_array_shrink(bond->flags, w);
        }
    }

done:
    // Add back the kept bonds, after the inferred ones and in their original order. Also on the early exit above:
    // missing coordinates mean nothing could be inferred, not that the kept bonds stop existing.
    // @NOTE: bond->count used to be left out here, so the re-added bonds sat past the end of the valid range.
    for (size_t i = 0; i < md_array_size(bond_pairs); ++i) {
        md_array_push(bond->pairs, bond_pairs[i], alloc);
        md_array_push(bond->flags, bond_flags[i], alloc);
        bond->count += 1;
    }
    for (size_t i = 0; i < num_noncovalent; ++i) {
        md_array_push(bond->pairs, noncovalent_pairs[i], alloc);
        md_array_push(bond->flags, noncovalent_flags[i], alloc);
    }
    bond->noncovalent_count = num_noncovalent;

    md_temp_end(temp_scope);
}

void md_util_system_infer_covalent_bonds(md_system_t* sys, const md_system_state_t* state) {
    md_system_topology_changed(sys);
    // Overwrite semantics: discard whatever was there so a repeat call cannot accumulate, except the bonds a user
    // added by hand and those the file gives (CONECT records), which md_util_infer_covalent_bonds keeps and reconciles
    // with the inferred ones: an inferred bond never stands in for one of them. The elastic bonds are kept as well,
    // whatever their origin: they are no bonds of the molecule, which is what is inferred, and nothing infers them.
    size_t num_user = 0;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const md_bond_origin_t origin = md_bond_origin(sys->bond.flags[i]);
        if (origin == MD_BOND_ORIGIN_USER || origin == MD_BOND_ORIGIN_FILE || md_bond_kind(sys->bond.flags[i]) == MD_BOND_KIND_ELASTIC) {
            sys->bond.pairs[num_user] = sys->bond.pairs[i];
            sys->bond.flags[num_user] = sys->bond.flags[i];
            num_user += 1;
        }
    }
    // The non-covalent bonds stay, right after the graph (md_bond_data_t)
    const size_t num_noncovalent = sys->bond.noncovalent_count;
    for (size_t i = 0; i < num_noncovalent; ++i) {
        sys->bond.pairs[num_user + i] = sys->bond.pairs[sys->bond.count + i];
        sys->bond.flags[num_user + i] = sys->bond.flags[sys->bond.count + i];
    }
    md_bond_conn_clear(&sys->bond.conn);
    sys->bond.count = num_user;
    md_array_shrink(sys->bond.pairs, num_user + num_noncovalent);
    md_array_shrink(sys->bond.flags, num_user + num_noncovalent);
    md_util_infer_covalent_bonds(&sys->bond, state, sys, sys->alloc);
    md_bond_build_connectivity(&sys->bond, state->num_atoms, sys->alloc);

    // The atoms of the metallic phases are those of the metallic bonds (METALLIC PHASES)
    if (md_array_size(sys->atom.flags) < sys->atom.count) {
        const size_t size = md_array_size(sys->atom.flags);
        md_array_resize(sys->atom.flags, sys->atom.count, sys->alloc);
        MEMSET(sys->atom.flags + size, 0, (sys->atom.count - size) * sizeof(md_atom_flags_t));
    }
    for (size_t i = 0; i < sys->atom.count; ++i) {
        sys->atom.flags[i] = (md_atom_flags_t)((uint32_t)sys->atom.flags[i] & ~(uint32_t)MD_ATOM_FLAG_METALLIC);
    }
    for (size_t i = 0; i < sys->bond.count; ++i) {
        if (md_bond_kind(sys->bond.flags[i]) != MD_BOND_KIND_METALLIC) continue;
        for (int k = 0; k < 2; ++k) {
            const md_atom_idx_t a = sys->bond.pairs[i].idx[k];
            if (a >= 0 && (size_t)a < sys->atom.count) sys->atom.flags[a] |= MD_ATOM_FLAG_METALLIC;
        }
    }
}

size_t md_util_system_infer_coordination(md_system_t* sys) {
    ASSERT(sys);
    size_t count = 0;
    if (!sys->bond.flags) return 0;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const md_atom_pair_t p = sys->bond.pairs[i];
        if (p.idx[0] < 0 || p.idx[1] < 0 || (size_t)p.idx[0] >= sys->atom.count || (size_t)p.idx[1] >= sys->atom.count) continue;
        if (!md_bond_is_graph(sys->bond.flags[i])) continue;   // A spring, of no chemistry
        const bool metal_a = is_metal(md_atom_atomic_number(&sys->atom, p.idx[0]));
        const bool metal_b = is_metal(md_atom_atomic_number(&sys->atom, p.idx[1]));
        // A kind which was given stays as it was given: a metal-ligand bond a user says is covalent
        const md_bond_flags_t f = sys->bond.flags[i];
        if (metal_a != metal_b && md_bond_kind(f) == MD_BOND_KIND_COVALENT && !(f & MD_BOND_FLAG_KIND_GIVEN)) {
            sys->bond.flags[i] = md_bond_flags_set_kind(f, MD_BOND_KIND_COORDINATE);
            count += 1;
        }
    }
    if (count) md_system_topology_changed(sys);
    return count;
}

// ### BONDS FORMED AND BROKEN ###

typedef struct delta_search_t {
    const md_bitset_t* scope;   // NULL when every atom is looked at, and the pairs come once each
    md_array(bond_pair_t) pairs;
    md_allocator_i* alloc;
} delta_search_t;

static void delta_pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    delta_search_t* s = (delta_search_t*)user_param;
    for (size_t k = 0; k < num_pairs; ++k) {
        const uint32_t i = i_idx[k];
        const uint32_t j = j_idx[k];
        if (i == j) continue;
        // i is of the scope: a pair with both atoms in it comes twice, and is taken as the one with the lower atom first
        if (s->scope && j < i && md_bitset_test(*s->scope, j)) continue;
        const bond_pair_t p = { (int)MIN(i, j), (int)MAX(i, j), sqrtf(ij_dist2[k]) };
        md_array_push(s->pairs, p, s->alloc);
    }
}

// What the delta needs to know of an atom, looked up as atoms come up rather than for every atom of the system
typedef struct delta_atoms_t {
    const md_system_t* sys;
    int8_t* res_water;          // [residue count] -1 not looked at yet, 0, 1
} delta_atoms_t;

static inline int delta_res(const delta_atoms_t* a, uint32_t atom) {
    return a->sys->residue.count ? (int)md_system_residue_find_by_atom_idx(a->sys, atom) : -1;
}

static inline bool delta_water(delta_atoms_t* a, int res) {
    if (res < 0) return false;
    if (a->res_water[res] < 0) {
        a->res_water[res] = (int8_t)(md_system_residue_kind(a->sys, res) == MD_RESIDUE_KIND_WATER || atoms_are_water(a->sys, md_system_residue_atom_range(a->sys, res)));
    }
    return a->res_water[res] != 0;
}

// The cutoff of a pair (factor times the sum of the radii, 0 when it is never bonded) and whether it is coordination
static inline bond_rule_t delta_rule(delta_atoms_t* a, uint32_t i, uint32_t j, float* out_sum_r) {
    const md_atomic_number_t zi = md_atom_atomic_number(&a->sys->atom, i);
    const md_atomic_number_t zj = md_atom_atomic_number(&a->sys->atom, j);
    *out_sum_r = element_covalent_radius(zi) + element_covalent_radius(zj);
    const bond_rule_t none = {0};
    if (zi == 0 || zj == 0) return none;
    if (md_atom_particle_kind(&a->sys->atom, i) == MD_PARTICLE_VIRTUAL_SITE || md_atom_particle_kind(&a->sys->atom, j) == MD_PARTICLE_VIRTUAL_SITE) return none;
    // The bonds within a metallic phase are its metallic bonds, which do not change (METALLIC PHASES)
    if ((md_atom_flags(&a->sys->atom, i) & MD_ATOM_FLAG_METALLIC) && (md_atom_flags(&a->sys->atom, j) & MD_ATOM_FLAG_METALLIC)) return none;
    const int ci = delta_res(a, i);
    const int cj = delta_res(a, j);
    return bond_rule(zi, zj, ci == cj, delta_water(a, ci) || delta_water(a, cj));
}

static inline float delta_distance(const md_system_state_t* state, uint32_t i, uint32_t j) {
    vec3_t dx = vec3_sub(state->xyz[j], state->xyz[i]);
    md_util_min_image_vec3(&dx, 1, &state->unitcell);
    return vec3_length(dx);
}

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

static inline uint64_t delta_pair_key(md_atom_pair_t p) {
    return ((uint64_t)(uint32_t)MIN(p.idx[0], p.idx[1]) << 32) | (uint32_t)MAX(p.idx[0], p.idx[1]);
}

// The position of atom in the sorted atoms, which holds it
static inline size_t delta_atom_slot(const uint32_t* atoms, size_t count, uint32_t atom) {
    size_t lo = 0, hi = count;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (atoms[mid] < atom) lo = mid + 1;
        else hi = mid;
    }
    ASSERT(lo < count && atoms[lo] == atom);
    return lo;
}

static int compare_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

void md_bond_delta_free(md_bond_delta_t* delta, md_allocator_i* alloc) {
    ASSERT(delta);
    ASSERT(alloc);
    md_bitset_free(&delta->broken, alloc);
    if (delta->formed_pair)  md_free(alloc, delta->formed_pair,  sizeof(md_atom_pair_t)  * delta->num_formed);
    if (delta->formed_flags) md_free(alloc, delta->formed_flags, sizeof(md_bond_flags_t) * delta->num_formed);
    MEMSET(delta, 0, sizeof(md_bond_delta_t));
}

bool md_util_bond_delta(md_bond_delta_t* out, const md_system_t* sys, const md_system_state_t* state, const md_bond_delta_desc_t* desc, md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(sys);
    ASSERT(state);
    ASSERT(desc);
    ASSERT(alloc);
    MEMSET(out, 0, sizeof(md_bond_delta_t));

    const size_t num_atoms = sys->atom.count;
    if (!state->xyz || state->num_atoms != num_atoms || num_atoms == 0) return false;
    if (md_system_is_coarse_grained(sys)) return false;
    if (sys->bond.count > 0 && !sys->bond.conn.offset) {
        MD_LOG_ERROR("The bonds of the system have no connectivity (md_bond_build_connectivity)");
        return false;
    }
    const float tol = CLAMP(desc->tolerance, 0.0f, 0.49f);
    const md_bitset_t* scope = desc->scope;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    delta_atoms_t atoms = { .sys = sys };
    atoms.res_water = md_temp_alloc_array(temp, int8_t, MAX(sys->residue.count, 1));
    MEMSET(atoms.res_water, -1, MAX(sys->residue.count, 1));

    // ## The bonds of the graph which touch the scope: broken beyond the band, and stretched beyond the cutoff
    // (they no longer count against what their atoms can take)
    const size_t num_bonds = sys->bond.count;
    md_bitmask_t broken    = md_bitmask_create(0, (uint32_t)num_bonds, temp_alloc);
    md_bitmask_t stretched = md_bitmask_create(0, (uint32_t)num_bonds, temp_alloc);
    {
        md_bitset_iter_t it = scope ? md_bitset_iter(*scope) : md_bitset_iter(md_bitset_range(0, (uint32_t)num_atoms));
        while (md_bitset_iter_next(&it)) {
            const uint32_t a = (uint32_t)it.idx;
            if (a >= num_atoms) break;
            md_bond_iter_t bit = md_bond_iter(&sys->bond, a);
            while (md_bond_iter_has_next(&bit)) {
                const uint32_t o = (uint32_t)md_bond_iter_atom_index(&bit);
                const uint32_t b = (uint32_t)md_bond_iter_bond_index(&bit);
                md_bond_iter_next(&bit);
                // Each bond once: from its lower atom when both are in the scope
                const bool o_in = scope ? md_bitset_test(*scope, o) : true;
                if (o_in && o < a) continue;
                float sum_r = 0.0f;
                const bond_rule_t rule = delta_rule(&atoms, a, o, &sum_r);
                if (sum_r <= 0.0f) continue;
                const md_atomic_number_t za = md_atom_atomic_number(&sys->atom, a);
                const md_atomic_number_t zo = md_atom_atomic_number(&sys->atom, o);
                if (za == 0 || zo == 0) continue;
                // A pair inference would not bond (an alkali ion and its ligand from a file, ...) breaks beyond the reach
                // of the search for candidates
                const float cutoff = (rule.factor > 0.0f ? rule.factor : bond_k_coord) * sum_r;
                const float d = delta_distance(state, a, o);
                if (d > cutoff) md_bitmask_set_bit(&stretched, b);
                const uint32_t origin_bit = 1u << md_bond_origin(sys->bond.flags[b]);
                if ((desc->breakable & origin_bit) && d > (1.0f + tol) * cutoff) {
                    md_bitmask_set_bit(&broken, b);
                }
            }
        }
    }

    // ## Candidates for the formed bonds: the pairs with an atom of the scope within reach
    md_array(bond_pair_t) formed = 0;
    {
        // The reach of the largest rule, between the largest atoms of the scope and of the system
        float max_rad_scope = 0.0f, max_rad_all = 0.0f;
        for (size_t i = 0; i < num_atoms; ++i) {
            const float r = element_covalent_radius(md_atom_atomic_number(&sys->atom, i));
            max_rad_all = MAX(max_rad_all, r);
            if (!scope || md_bitset_test(*scope, (uint32_t)i)) max_rad_scope = MAX(max_rad_scope, r);
        }
        const double cutoff = MAX(1.0, (double)(max_rad_scope + max_rad_all) * bond_k_coord);

        delta_search_t search = { .scope = scope, .alloc = temp_alloc };
        md_coord_stream_t coords = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), NULL, num_atoms);
        md_spatial_acc_t acc = { .alloc = temp_alloc };
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &coords, .cutoff = cutoff, .unitcell = &state->unitcell });
        if (!scope) {
            md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, delta_pair_callback, &search);
        } else {
            const size_t num_scope = md_bitset_count(*scope);
            if (num_scope > 0) {
                int* idx = md_temp_alloc_array(temp, int, num_scope);
                size_t n = 0;
                md_bitset_iter_t it = md_bitset_iter(*scope);
                while (md_bitset_iter_next(&it) && n < num_scope) {
                    if (it.idx < num_atoms) idx[n++] = (int)it.idx;
                }
                md_coord_stream_t ext = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), idx, n);
                md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext, cutoff, delta_pair_callback, &search, MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX);
            }
        }

        // The non-covalent bonds of the system, whose pairs form nothing
        const size_t num_nc = sys->bond.noncovalent_count;
        uint64_t* nc_keys = md_temp_alloc_array(temp, uint64_t, MAX(num_nc, 1));
        for (size_t k = 0; k < num_nc; ++k) nc_keys[k] = delta_pair_key(sys->bond.pairs[num_bonds + k]);
        qsort(nc_keys, num_nc, sizeof(uint64_t), compare_u64);

        for (size_t k = 0; k < md_array_size(search.pairs); ++k) {
            const bond_pair_t p = search.pairs[k];
            float sum_r = 0.0f;
            const bond_rule_t rule = delta_rule(&atoms, (uint32_t)p.atom_i, (uint32_t)p.atom_j, &sum_r);
            if (rule.factor <= 0.0f) continue;
            if (p.dist <= bond_k_min * sum_r || p.dist >= (1.0f - tol) * rule.factor * sum_r) continue;
            if (md_bond_find(&sys->bond, p.atom_i, p.atom_j) != -1) continue;
            const md_atom_pair_t pair = {{ p.atom_i, p.atom_j }};
            const uint64_t key = delta_pair_key(pair);
            if (num_nc && bsearch(&key, nc_keys, num_nc, sizeof(uint64_t), compare_u64)) continue;
            bond_pair_t f = p;
            // The kind rides in the sign of the distance until the flags are written
            if (rule.coordinate) f.dist = -f.dist;
            md_array_push(formed, f, temp_alloc);
        }
    }

    // ## An organohalogen coordinates nothing: a halogen bonded to a non-metal in the state, by a bond of the graph which
    // is not broken in it or by one formed in it (is_halogen)
    if (md_array_size(formed)) {
        md_array(uint32_t) halogens = 0;    // With a covalent bond formed
        for (size_t k = 0; k < md_array_size(formed); ++k) {
            if (formed[k].dist < 0.0f) continue;
            const uint32_t i = (uint32_t)formed[k].atom_i;
            const uint32_t j = (uint32_t)formed[k].atom_j;
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
            const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, j);
            if (is_halogen(zi) && !is_metal(zj)) md_array_push(halogens, i, temp_alloc);
            if (is_halogen(zj) && !is_metal(zi)) md_array_push(halogens, j, temp_alloc);
        }
        const size_t num_halogens = md_array_size(halogens);
        if (num_halogens) qsort(halogens, num_halogens, sizeof(uint32_t), compare_u32);

        size_t w = 0;
        for (size_t k = 0; k < md_array_size(formed); ++k) {
            if (formed[k].dist < 0.0f) {
                const uint32_t i = (uint32_t)formed[k].atom_i;
                const uint32_t j = (uint32_t)formed[k].atom_j;
                const uint32_t h = is_halogen(md_atom_atomic_number(&sys->atom, i)) ? i : j;
                if (is_halogen(md_atom_atomic_number(&sys->atom, h))) {
                    bool organo = num_halogens && bsearch(&h, halogens, num_halogens, sizeof(uint32_t), compare_u32);
                    md_bond_iter_t it = md_bond_iter(&sys->bond, h);
                    for (; md_bond_iter_has_next(&it) && !organo; md_bond_iter_next(&it)) {
                        const uint32_t o = (uint32_t)md_bond_iter_atom_index(&it);
                        const uint32_t b = (uint32_t)md_bond_iter_bond_index(&it);
                        organo = md_bond_kind((md_bond_flags_t)md_bond_iter_bond_flags(&it)) == MD_BOND_KIND_COVALENT &&
                                 !is_metal(md_atom_atomic_number(&sys->atom, o)) && !md_bitmask_test_bit(&broken, b);
                    }
                    if (organo) continue;
                }
            }
            formed[w++] = formed[k];
        }
        md_array_shrink(formed, w);
    }

    // ## What the atoms of the formed bonds can take. The bonds of the graph count unless broken or stretched, and
    // only a formed bond is ever removed.
    const size_t num_candidates = md_array_size(formed);
    md_array(bool) keep = 0;
    if (num_candidates) {
        md_array_resize(keep, num_candidates, temp_alloc);
        for (size_t k = 0; k < num_candidates; ++k) keep[k] = true;

        uint32_t* ends = md_temp_alloc_array(temp, uint32_t, 2 * num_candidates);
        for (size_t k = 0; k < num_candidates; ++k) {
            ends[2 * k + 0] = (uint32_t)formed[k].atom_i;
            ends[2 * k + 1] = (uint32_t)formed[k].atom_j;
        }
        qsort(ends, 2 * num_candidates, sizeof(uint32_t), compare_u32);
        size_t num_ends = 0;
        for (size_t k = 0; k < 2 * num_candidates; ++k) {
            if (num_ends == 0 || ends[num_ends - 1] != ends[k]) ends[num_ends++] = ends[k];
        }
        int* num_cov   = md_temp_alloc_array(temp, int, num_ends);
        int* num_coord = md_temp_alloc_array(temp, int, num_ends);
        MEMSET(num_cov,   0, sizeof(int) * num_ends);
        MEMSET(num_coord, 0, sizeof(int) * num_ends);

        for (size_t e = 0; e < num_ends; ++e) {
            const uint32_t a = ends[e];
            md_bond_iter_t bit = md_bond_iter(&sys->bond, a);
            while (md_bond_iter_has_next(&bit)) {
                const uint32_t o = (uint32_t)md_bond_iter_atom_index(&bit);
                const uint32_t b = (uint32_t)md_bond_iter_bond_index(&bit);
                const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&bit);
                md_bond_iter_next(&bit);
                if (md_bitmask_test_bit(&broken, b)) continue;
                // A bond of an atom outside of the scope has not been measured
                bool stretch = md_bitmask_test_bit(&stretched, b);
                if (!stretch && scope && !md_bitset_test(*scope, a) && !md_bitset_test(*scope, o)) {
                    float sum_r = 0.0f;
                    const bond_rule_t rule = delta_rule(&atoms, a, o, &sum_r);
                    stretch = sum_r > 0.0f && delta_distance(state, a, o) > (rule.factor > 0.0f ? rule.factor : bond_k_coord) * sum_r;
                }
                if (stretch) continue;
                if (md_bond_kind(f) == MD_BOND_KIND_COORDINATE) num_coord[e] += 1;
                else                             num_cov[e]   += 1;
            }
        }
        for (size_t k = 0; k < num_candidates; ++k) {
            const size_t ei = delta_atom_slot(ends, num_ends, (uint32_t)formed[k].atom_i);
            const size_t ej = delta_atom_slot(ends, num_ends, (uint32_t)formed[k].atom_j);
            int* n = formed[k].dist < 0.0f ? num_coord : num_cov;
            n[ei] += 1;
            n[ej] += 1;
        }

        // In the order of md_util_infer_covalent_bonds: bonds between residues first, then the furthest from the radii
        md_array(bond_excess_t) over = 0;
        for (size_t k = 0; k < num_candidates; ++k) {
            const uint32_t i = (uint32_t)formed[k].atom_i;
            const uint32_t j = (uint32_t)formed[k].atom_j;
            const size_t ei = delta_atom_slot(ends, num_ends, i);
            const size_t ej = delta_atom_slot(ends, num_ends, j);
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
            const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, j);
            bool exceeds = false;
            if (formed[k].dist < 0.0f) {
                const size_t m = is_metal(zi) ? ei : ej;
                exceeds = num_coord[m] > max_coordination(is_metal(zi) ? zi : zj);
            } else {
                const int mi = max_covalent_bonds(zi);
                const int mj = max_covalent_bonds(zj);
                exceeds = (mi && num_cov[ei] > mi) || (mj && num_cov[ej] > mj);
            }
            if (exceeds) {
                const float sum_r = element_covalent_radius(zi) + element_covalent_radius(zj);
                bond_excess_t x = {
                    .inter = delta_res(&atoms, i) != delta_res(&atoms, j),
                    .deviation = sum_r > 0 ? fabsf(fabsf(formed[k].dist) / sum_r - 1.0f) : 0.0f,
                    .idx = (uint32_t)k,
                };
                md_array_push(over, x, temp_alloc);
            }
        }
        const size_t num_over = md_array_size(over);
        if (num_over) {
            qsort(over, num_over, sizeof(bond_excess_t), compare_bond_excess);
            for (size_t o = 0; o < num_over; ++o) {
                const uint32_t k = over[o].idx;
                const uint32_t i = (uint32_t)formed[k].atom_i;
                const uint32_t j = (uint32_t)formed[k].atom_j;
                const size_t ei = delta_atom_slot(ends, num_ends, i);
                const size_t ej = delta_atom_slot(ends, num_ends, j);
                const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
                const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, j);
                if (formed[k].dist < 0.0f) {
                    const size_t m = is_metal(zi) ? ei : ej;
                    if (num_coord[m] <= max_coordination(is_metal(zi) ? zi : zj)) continue;
                    num_coord[ei] -= 1;
                    num_coord[ej] -= 1;
                } else {
                    const int mi = max_covalent_bonds(zi);
                    const int mj = max_covalent_bonds(zj);
                    if (!((mi && num_cov[ei] > mi) || (mj && num_cov[ej] > mj))) continue;
                    num_cov[ei] -= 1;
                    num_cov[ej] -= 1;
                }
                keep[k] = false;
            }
        }
    }

    // ## The delta, in the caller's allocator
    out->broken = md_bitmask_extract(&broken, alloc);
    size_t num_formed = 0;
    for (size_t k = 0; k < num_candidates; ++k) num_formed += keep[k];
    if (num_formed) {
        // Sorted by pair: the candidates came in the order of the search
        uint64_t* keys = md_temp_alloc_array(temp, uint64_t, num_formed);
        size_t n = 0;
        for (size_t k = 0; k < num_candidates; ++k) {
            if (!keep[k]) continue;
            const md_atom_pair_t pair = {{ formed[k].atom_i, formed[k].atom_j }};
            keys[n++] = delta_pair_key(pair);
        }
        qsort(keys, num_formed, sizeof(uint64_t), compare_u64);
        out->formed_pair  = md_alloc(alloc, sizeof(md_atom_pair_t)  * num_formed);
        out->formed_flags = md_alloc(alloc, sizeof(md_bond_flags_t) * num_formed);
        const md_bond_flags_t inferred = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_INFERRED);
        for (size_t k = 0; k < num_formed; ++k) {
            const md_atom_pair_t pair = {{ (md_atom_idx_t)(keys[k] >> 32), (md_atom_idx_t)(keys[k] & 0xFFFFFFFFu) }};
            out->formed_pair[k] = pair;
            // A rule gives a coordination bond for a metal and a non-metal, and only for them
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, pair.idx[0]);
            const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, pair.idx[1]);
            out->formed_flags[k] = (is_metal(zi) != is_metal(zj)) ? md_bond_flags_set_kind(inferred, MD_BOND_KIND_COORDINATE) : inferred;
        }
        out->num_formed = num_formed;
    }

    md_bitmask_free(&broken);
    md_bitmask_free(&stretched);
    md_temp_end(temp);
    return true;
}

// No upper bound for an amino acid either: a modified one is larger than any standard one (a palmitoylated cysteine is
// 58 atoms with its hydrogens), and classify_amino_acid tells ligands named like amino acids by their bonds
#define MIN_RES_LEN 4

// Nor for a nucleotide: a hypermodified one is larger than any standard one (wybutosine, 39 heavy atoms), as
// is a modified one with its hydrogens, and classify_nucleotide tells nucleotide ligands (NAD, FAD, CoA, ...) by their
// bonds rather than their size. The nucleic backbone (MD_UTIL_INFER_BACKBONE_BIT) takes any size, and the two agree.
#define MIN_NUC_LEN 6

#define MIN_POLYMER_RES_LEN 2

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

// helper function to get atom indices within a 'tree branch' of the connectivity graph of a system.
static inline bool get_branch_atom_indices(md_array(int)* out_indices, const md_bond_data_t* bond, int start_atom, int exclude_atom, md_allocator_i* alloc) {
    md_array_push(*out_indices, start_atom, alloc);
    md_bond_iter_t iter = md_bond_iter(bond, start_atom);
    while (md_bond_iter_has_next(&iter)) {
        int neighbor = md_bond_iter_atom_index(&iter);
        if (neighbor != exclude_atom) {
            if (!get_branch_atom_indices(out_indices, bond, neighbor, start_atom, alloc)) {
                return false;
            }
        }
        md_bond_iter_next(&iter);
    }
    return true;
}

// Gives the flag to start_atom and to what its bonds reach from it within the residue (range): not past exclude_atom,
// nor into an atom with any of stop_flags (the backbone, for a side chain). A bond leaving the residue is not followed:
// the peptide bond (which a proline side chain reaches through its N), a disulfide, a metal coordinated to it.
static inline void topo_floodfill_flag(md_system_t* sys, int start_atom, int exclude_atom, md_atom_flags_t flag, md_urange_t range, md_atom_flags_t stop_flags) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_array(int) stack = 0;
    md_array_push(stack, start_atom, alloc);
    while (md_array_size(stack) > 0) {
        int atom_idx = stack[md_array_size(stack) - 1];
        md_array_pop(stack);
        if ((sys->atom.flags[atom_idx] & flag) == 0) {
            sys->atom.flags[atom_idx] |= flag;
            md_bond_iter_t iter = md_bond_iter(&sys->bond, atom_idx);
            while (md_bond_iter_has_next(&iter)) {
                int nbr = md_bond_iter_atom_index(&iter);
                if (nbr != exclude_atom && range.beg <= (uint32_t)nbr && (uint32_t)nbr < range.end &&
                    (sys->atom.flags[nbr] & flag) == 0 && (sys->atom.flags[nbr] & stop_flags) == 0) {
                    md_array_push(stack, nbr, alloc);
                }
                md_bond_iter_next(&iter);
            }
        }
    }
    md_temp_end(temp);
}

static inline void set_residue_kind(md_system_t* sys, size_t res_idx, md_residue_kind_t kind) {
    sys->residue.flags[res_idx] = md_residue_flags_set_kind(sys->residue.flags[res_idx], kind);
}

// Whether the atoms of the range are all coarse grained beads (or their virtual sites)
static bool atoms_are_beads(const md_system_t* sys, md_urange_t range) {
    if (range.beg >= range.end) return false;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if (md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_ATOM || md_atom_atomic_number(&sys->atom, i) != 0) return false;
    }
    return true;
}

// True when every neighbor of the atom other than 'except' is of element allowed_z: the atom ends its chain
static bool atom_is_chain_end(const md_system_t* sys, int atom, int except, md_atomic_number_t allowed_z) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, atom);
    while (md_bond_iter_has_next(&it)) {
        const int nbr = md_bond_iter_atom_index(&it);
        if (nbr != except && md_atom_atomic_number(&sys->atom, nbr) != allowed_z) {
            return false;
        }
        md_bond_iter_next(&it);
    }
    return true;
}

// Whether the atom is bonded to an atom of element z outside of the range (another residue)
static bool atom_bonded_outside(const md_system_t* sys, int atom, md_urange_t range, md_atomic_number_t z) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, atom);
    while (md_bond_iter_has_next(&it)) {
        const int nbr = md_bond_iter_atom_index(&it);
        if ((nbr < (int)range.beg || nbr >= (int)range.end) && md_atom_atomic_number(&sys->atom, nbr) == z) {
            return true;
        }
        md_bond_iter_next(&it);
    }
    return false;
}

// Amino acid whose backbone atoms were found by name: verify them by their bonds and give the atoms their roles.
// It must also be a standard residue by name or be linked into a chain by a peptide bond: plenty of ligands name
// their atoms N, CA, C and O without being amino acids.
static bool classify_amino_acid(md_system_t* sys, size_t res_idx, const md_amino_acid_atoms_t* aa) {
    // Ensure that there are bonds between the backbone atoms, otherwise it might be a false posititive
    if (md_bond_find(&sys->bond, aa->n,  aa->ca) == -1 ||
        md_bond_find(&sys->bond, aa->ca, aa->c)  == -1 ||
        md_bond_find(&sys->bond, aa->c,  aa->o)  == -1) {
        return false;
    }
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    if (!md_util_resname_amino_acid(md_system_residue_name(sys, res_idx)) &&
        !atom_bonded_outside(sys, aa->n, range, MD_Z_C) && !atom_bonded_outside(sys, aa->c, range, MD_Z_N)) {
        return false;
    }

    set_residue_kind(sys, res_idx, MD_RESIDUE_KIND_AMINO_ACID);
    sys->residue.flags[res_idx] |= MD_RESIDUE_FLAG_RESOLVED;

    // O and H are not part of the backbone proper
    sys->atom.flags[aa->n]  |= MD_ATOM_FLAG_BACKBONE;
    sys->atom.flags[aa->ca] |= MD_ATOM_FLAG_BACKBONE;
    sys->atom.flags[aa->c]  |= MD_ATOM_FLAG_BACKBONE;

    if (aa->cb != -1) {
        // Flood fill sidechain from CB, this will mark all connected atoms as sidechain
        topo_floodfill_flag(sys, aa->cb, aa->ca, MD_ATOM_FLAG_SIDE_CHAIN, range, MD_ATOM_FLAG_BACKBONE);
    }

    // A free amino group (only hydrogens besides CA) begins the chain, a free carboxyl group (only oxygens besides CA) ends it
    if (atom_is_chain_end(sys, aa->n, aa->ca, MD_Z_H)) {
        topo_floodfill_flag(sys, aa->n, aa->ca, MD_ATOM_FLAG_TERMINAL_BEG, range, MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_SIDE_CHAIN);
        sys->residue.flags[res_idx] |= MD_RESIDUE_FLAG_TERMINAL_BEG;
    }
    if (atom_is_chain_end(sys, aa->c, aa->ca, MD_Z_O)) {
        topo_floodfill_flag(sys, aa->c, aa->ca, MD_ATOM_FLAG_TERMINAL_END, range, MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_SIDE_CHAIN);
        sys->residue.flags[res_idx] |= MD_RESIDUE_FLAG_TERMINAL_END;
    }
    return true;
}

// Nucleotide whose backbone atoms were found by name: verify them by their bonds and give the atoms their roles.
// It must also be a standard residue by name or be linked into a chain by a phosphodiester bond: the nucleotide
// ligands (GTP, ATP, NAD, FAD, SAM, ...) are named like nucleotides atom by atom.
static bool classify_nucleotide(md_system_t* sys, size_t res_idx, const md_nucleic_acid_atoms_t* na) {
    // O3 - C3 - C4 - C5 - O5  (and P for nucleotides if not terminal)
    if (md_bond_find(&sys->bond, na->o3, na->c3) == -1 ||
        md_bond_find(&sys->bond, na->c3, na->c4) == -1 ||
        md_bond_find(&sys->bond, na->c4, na->c5) == -1 ||
        md_bond_find(&sys->bond, na->c5, na->o5) == -1 ||
        (na->p != -1 && md_bond_find(&sys->bond, na->o5, na->p) == -1)) {
        return false;
    }
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    if (!md_util_resname_nucleotide(md_system_residue_name(sys, res_idx)) &&
        !(na->p != -1 && atom_bonded_outside(sys, na->p, range, MD_Z_O)) && !atom_bonded_outside(sys, na->o3, range, MD_Z_P)) {
        return false;
    }

    set_residue_kind(sys, res_idx, MD_RESIDUE_KIND_NUCLEOTIDE);
    sys->residue.flags[res_idx] |= MD_RESIDUE_FLAG_RESOLVED;

    if (na->p != -1) {
        sys->atom.flags[na->p] |= MD_ATOM_FLAG_BACKBONE;
    }
    sys->atom.flags[na->o5] |= MD_ATOM_FLAG_BACKBONE;
    sys->atom.flags[na->c5] |= MD_ATOM_FLAG_BACKBONE;
    sys->atom.flags[na->c4] |= MD_ATOM_FLAG_BACKBONE;
    sys->atom.flags[na->c3] |= MD_ATOM_FLAG_BACKBONE;
    sys->atom.flags[na->o3] |= MD_ATOM_FLAG_BACKBONE;

    if (na->c1 != -1 && na->c2 != -1 && na->o4 != -1) {
        sys->atom.flags[na->c1] |= MD_ATOM_FLAG_NUCLEOSIDE;
        sys->atom.flags[na->c2] |= MD_ATOM_FLAG_NUCLEOSIDE;
        sys->atom.flags[na->c3] |= MD_ATOM_FLAG_NUCLEOSIDE;
        sys->atom.flags[na->c4] |= MD_ATOM_FLAG_NUCLEOSIDE;
        sys->atom.flags[na->o4] |= MD_ATOM_FLAG_NUCLEOSIDE;

        // The atom of the base bonded to the sugar (C1'), from which the base is flood filled
        int base_start_atom = -1;
        md_bond_iter_t it = md_bond_iter(&sys->bond, na->c1);
        while (md_bond_iter_has_next(&it)) {
            int n = md_bond_iter_atom_index(&it);
            if (n != na->c2 && n != na->o4 && md_atom_atomic_number(&sys->atom, n) != MD_Z_H) {
                base_start_atom = n;
                break;
            }
            md_bond_iter_next(&it);
        }
        if (base_start_atom != -1) {
            topo_floodfill_flag(sys, base_start_atom, na->c1, MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_NUCLEOSIDE, range, MD_ATOM_FLAG_BACKBONE);
        }

        // A free 3' hydroxyl ends the chain, a free 5' hydroxyl begins it
        if (atom_is_chain_end(sys, na->o3, na->c3, MD_Z_H)) {
            topo_floodfill_flag(sys, na->o3, na->c3, MD_ATOM_FLAG_TERMINAL_END, range, 0);
            sys->residue.flags[res_idx] |= MD_RESIDUE_FLAG_TERMINAL_END;
        }
        if (atom_is_chain_end(sys, na->o5, na->c5, MD_Z_H)) {
            topo_floodfill_flag(sys, na->o5, na->c5, MD_ATOM_FLAG_TERMINAL_BEG, range, 0);
            sys->residue.flags[res_idx] |= MD_RESIDUE_FLAG_TERMINAL_BEG;
        }
    }
    return true;
}

// The ends of coarse grained chains: an amino acid or a nucleotide of beads whose backbone is not bonded to that of the
// residue before (after) it begins (ends) its chain. Martini 3 gives the backbone bead at either end of a protein chain
// a regular type, charged (Q5) or neutral (P6, P1, P2) alike, where Gly, Ala, Val and Pro otherwise have a small one: a
// small backbone bead (54) at an end gets the mass and size of a regular one.
static void infer_coarse_grained_chain_ends(md_system_t* sys) {
    if (!md_system_is_coarse_grained(sys) || !sys->residue.flags || sys->bond.conn.offset == NULL) return;
    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        if (kind != MD_RESIDUE_KIND_AMINO_ACID && kind != MD_RESIDUE_KIND_NUCLEOTIDE) continue;
        if (sys->residue.flags[ci] & MD_RESIDUE_FLAG_RESOLVED) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, ci);
        if (range.beg >= range.end || md_atom_particle_kind(&sys->atom, range.beg) != MD_PARTICLE_BEAD) continue;
        const md_urange_t prev = ci > 0 ? md_system_residue_atom_range(sys, ci - 1) : (md_urange_t){0, 0};
        const md_urange_t next = ci + 1 < sys->residue.count ? md_system_residue_atom_range(sys, ci + 1) : (md_urange_t){0, 0};

        bool to_prev = false, to_next = false;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE)) continue;
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            while (md_bond_iter_has_next(&it)) {
                const uint32_t j = (uint32_t)md_bond_iter_atom_index(&it);
                if (md_atom_flags(&sys->atom, j) & MD_ATOM_FLAG_BACKBONE) {
                    if (prev.beg <= j && j < prev.end) to_prev = true;
                    if (next.beg <= j && j < next.end) to_next = true;
                }
                md_bond_iter_next(&it);
            }
        }
        if (!to_prev) sys->residue.flags[ci] |= MD_RESIDUE_FLAG_TERMINAL_BEG;
        if (!to_next) sys->residue.flags[ci] |= MD_RESIDUE_FLAG_TERMINAL_END;

        if (kind == MD_RESIDUE_KIND_AMINO_ACID && (!to_prev || !to_next)) {
            for (uint32_t i = range.beg; i < range.end; ++i) {
                const md_particle_type_idx_t t = sys->atom.type_idx[i];
                if (md_particle_type_mass(&sys->atom.type, t) != 54.0f || !str_eq_cstr(md_particle_type_name(&sys->atom.type, t), "BB")) continue;
                sys->atom.type_idx[i] = find_or_add_particle_type(&sys->atom.type, md_particle_type_name(&sys->atom.type, t), 0, 72.0f, MARTINI_R,
                    md_particle_type_color(&sys->atom.type, t), md_particle_type_flags(&sys->atom.type, t), sys->alloc);
            }
        }
    }
}

// Classifies the residues (amino acid, nucleotide, water, ion, lipid) from their names, atoms and bonds, and gives the
// atoms of the amino acids and nucleotides whose backbone could be verified their roles.
bool md_util_system_infer_residue_flags(md_system_t* sys) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }
    if (sys->residue.count == 0) {
        return true;
    }
    ASSERT(sys->alloc);
    md_system_topology_changed(sys);

    if (md_array_size(sys->residue.flags) < sys->residue.count) {
        const size_t size = md_array_size(sys->residue.flags);
        md_array_resize(sys->residue.flags, sys->residue.count, sys->alloc);
        MEMSET(sys->residue.flags + size, 0, (sys->residue.count - size) * sizeof(md_residue_flags_t));
    }
    if (md_array_size(sys->atom.flags) < sys->atom.count) {
        const size_t size = md_array_size(sys->atom.flags);
        md_array_resize(sys->atom.flags, sys->atom.count, sys->alloc);
        MEMSET(sys->atom.flags + size, 0, (sys->atom.count - size) * sizeof(md_atom_flags_t));
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(int) ambigous_amino_acid_indices = 0;
    md_array(int) ambigous_nucleotide_indices = 0;

    for (size_t residue_idx = 0; residue_idx < sys->residue.count; ++residue_idx) {
        str_t res_name = md_system_residue_name(sys, residue_idx);
        md_urange_t res_range = md_system_residue_atom_range(sys, residue_idx);
        size_t len = (size_t)(res_range.end - res_range.beg);

        const md_residue_kind_t kind = md_system_residue_kind(sys, residue_idx);
        if (kind == MD_RESIDUE_KIND_WATER || kind == MD_RESIDUE_KIND_ION || kind == MD_RESIDUE_KIND_LIPID) {
            // Already classified (by a predefined coarse grained type)
            continue;
		}

        // A molecule of the Martini force fields (a lipid, water, an ion, a sugar, ...): the kind of its topology. Its
        // name alone may say something else: PCA, a lipid of Martini 2, is pyroglutamic acid in the PDB.
        if (len <= 255 && atoms_are_beads(sys, res_range)) {
            str_t labels[256];
            uint8_t map[256];
            for (size_t i = 0; i < len; ++i) labels[i] = md_atom_name(&sys->atom, res_range.beg + i);
            const martini_db_residue_t* mol = martini_db_match_molecule(MARTINI_ANY, res_name, labels, len, map);
            // Renamed (CHL for CHOL): by its beads, unless it is named as a residue of a chain
            if (!mol && !md_util_resname_amino_acid(res_name) && !md_util_resname_nucleotide(res_name)) {
                mol = martini_db_match_beads(MARTINI_ANY, labels, NULL, len, map);
            }
            if (mol) {
                set_residue_kind(sys, residue_idx, mol->kind);
                continue;
            }
        }

        // Amino acids and nucleotides are classified whatever their kind, as a predefined coarse grained type may
        // have given the kind but the atoms of an atomistic one still want their roles
        md_amino_acid_atoms_t prot_atoms = {0};
        if (MIN_RES_LEN <= len && md_util_amino_acid_atoms_extract(&prot_atoms, &sys->atom, res_range)) {
            // Atoms named like an amino acid's backbone but not bonded like one, or a ligand which is neither a standard
            // residue nor linked into a chain: not an amino acid
            classify_amino_acid(sys, residue_idx, &prot_atoms);
            continue;
        }
        if (md_util_resname_amino_acid(res_name)) {
            set_residue_kind(sys, residue_idx, MD_RESIDUE_KIND_AMINO_ACID);
            md_array_push(ambigous_amino_acid_indices, (int)residue_idx, temp_alloc);
            continue;
        }

        md_nucleic_acid_atoms_t nucl_atoms = {0};
        if (MIN_NUC_LEN <= len && md_util_nucleic_acid_atoms_extract(&nucl_atoms, &sys->atom, res_range)) {
            classify_nucleotide(sys, residue_idx, &nucl_atoms);
            continue;
        }
        if (md_util_resname_nucleotide(res_name)) {
            set_residue_kind(sys, residue_idx, MD_RESIDUE_KIND_NUCLEOTIDE);
            md_array_push(ambigous_nucleotide_indices, (int)residue_idx, temp_alloc);
            continue;
        }

        // Water by its atoms, or by name when its hydrogens are left out (united atom, coarse grained)
        if (atoms_are_water(sys, res_range) ||
           ((len == 1 || len == 3) && md_util_resname_water(res_name)) ||
           (len == 1 && md_util_resname_water(md_atom_name(&sys->atom, res_range.beg))))
        {
            set_residue_kind(sys, residue_idx, MD_RESIDUE_KIND_WATER);
        } else if (len == 1 && (md_util_resname_ion(res_name) || monatomic_ion_element(md_atom_atomic_number(&sys->atom, res_range.beg)))) {
            set_residue_kind(sys, residue_idx, MD_RESIDUE_KIND_ION);
        } else if (md_util_resname_lipid(res_name)) {
            set_residue_kind(sys, residue_idx, MD_RESIDUE_KIND_LIPID);
        }
    }

    infer_coarse_grained_chain_ends(sys);

    const size_t max_warnings = 10;
    
    // Print warnings about unsafe assignments
    size_t num_ambigous_amino_acids = md_array_size(ambigous_amino_acid_indices);
    if (num_ambigous_amino_acids > 0) {
        MD_LOG_DEBUG("Structure contains %zu amino acid residue names, but their backbone cannot be mapped from atom labels:", num_ambigous_amino_acids);
        for (size_t i = 0; i < MIN(num_ambigous_amino_acids, max_warnings); ++i) {
            int res_idx = ambigous_amino_acid_indices[i];
            MD_LOG_DEBUG(" - [%d] " STR_FMT " (%d)", res_idx + 1, STR_ARG(md_system_residue_name(sys, res_idx)), md_system_residue_seq_id(sys, res_idx));
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
            MD_LOG_DEBUG(" - [%d] " STR_FMT " (%d)", res_idx + 1, STR_ARG(md_system_residue_name(sys, res_idx)), md_system_residue_seq_id(sys, res_idx));
        }
        if (num_ambigous_nucleotides > max_warnings) {
            MD_LOG_DEBUG(" - And %zu more...", num_ambigous_nucleotides - max_warnings);
        }
    }
    
    md_temp_end(temp);
    return true;
}

static void clear_chains(md_system_t* sys, md_allocator_i* alloc) {
    md_array_free(sys->chain.id, alloc);
    md_array_free(sys->chain.run_offset, alloc);
    md_array_free(sys->chain.run, alloc);
    md_array_free(sys->chain.res_chain, alloc);
    MEMSET(&sys->chain, 0, sizeof(sys->chain));
}

static bool polymer_chain_ends_between(const md_system_t* sys, size_t a, size_t b);

static void clear_instances(md_system_t* sys, md_allocator_i* alloc) {
    md_array_free(sys->instance.id, alloc);
    md_array_free(sys->instance.auth_id, alloc);
    md_array_free(sys->instance.res_offset, alloc);
    md_array_free(sys->instance.part_type_idx, alloc);
    md_array_free(sys->instance.mol_type_idx, alloc);
    MEMSET(&sys->instance, 0, sizeof(sys->instance));
}

bool md_util_system_instances_from_chain_ids(md_system_t* sys, const str_t res_chain_id[]) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;
    clear_instances(sys, alloc);
    md_system_topology_changed(sys);

    const size_t num_res = sys->residue.count;
    if (!res_chain_id || num_res == 0) return true;

    // No ids at all is no units: the data gives none
    bool any_id = false;
    for (size_t c = 0; c < num_res && !any_id; ++c) any_id = !str_empty(str_trim(res_chain_id[c]));
    if (!any_id) return true;

    // The runs of one id. Every residue is in one, as the instances are back to back; a run without an id is an
    // instance the data does not name (md_instance_data_t)
    size_t i = 0;
    while (i < num_res) {
        const str_t id = str_trim(res_chain_id[i]);
        size_t j = i + 1;
        while (j < num_res && str_eq(id, str_trim(res_chain_id[j]))) ++j;
        md_array_push(sys->instance.id, (md_label_t){0}, alloc);
        md_array_push(sys->instance.auth_id, make_label(id), alloc);
        md_array_push(sys->instance.res_offset, (uint32_t)i, alloc);
        md_array_push(sys->instance.part_type_idx, -1, alloc);
        md_array_push(sys->instance.mol_type_idx, -1, alloc);
        sys->instance.count += 1;
        i = j;
    }
    md_array_push(sys->instance.res_offset, (uint32_t)num_res, alloc);
    return true;
}

bool md_util_system_infer_chains(md_system_t* sys) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;
    clear_chains(sys, alloc);
    md_system_topology_changed(sys);

    const size_t num_inst = sys->instance.count;
    const size_t num_res = sys->residue.count;
    if (num_inst == 0 || num_res == 0) return true;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The chain of each id, by the hash of the id
    md_hashmap32_t chain_of_id = { .allocator = temp_arena };
    md_array_resize(sys->chain.res_chain, num_res, alloc);
    for (size_t c = 0; c < num_res; ++c) sys->chain.res_chain[c] = -1;

    // The runs, in the order they are found, each with its chain. A chain is a label: only what carries an id is in one.
    md_array(md_urange_t) runs = 0;
    md_array(md_chain_idx_t) run_chain = 0;

    for (size_t i = 0; i < num_inst; ++i) {
        const str_t id = str_trim(md_system_instance_auth_id(sys, i));
        const md_urange_t residues = md_system_instance_residue_range(sys, i);
        if (str_empty(id) || residues.beg >= residues.end || residues.end > num_res) continue;
        const uint64_t key = md_hash64_str(id, 0);
        const uint32_t* found = md_hashmap_get(&chain_of_id, key);
        md_chain_idx_t chain;
        if (found) {
            chain = (md_chain_idx_t)*found;
        } else {
            chain = (md_chain_idx_t)sys->chain.count++;
            md_array_push(sys->chain.id, make_label(id), alloc);
            md_hashmap_add(&chain_of_id, key, (uint32_t)chain);
        }
        md_array_push(runs, residues, temp_arena);
        md_array_push(run_chain, chain, temp_arena);
    }

    // The runs of each chain, in the order of the chains, the consecutive ones of a chain joined
    if (sys->chain.count) {
        const size_t num_runs = md_array_size(runs);
        uint32_t* count = md_temp_alloc_zero_array(temp, uint32_t, sys->chain.count + 1);
        for (size_t r = 0; r < num_runs; ++r) count[run_chain[r] + 1] += 1;
        for (size_t c = 0; c < sys->chain.count; ++c) count[c + 1] += count[c];
        md_urange_t* sorted = md_temp_alloc_array(temp, md_urange_t, MAX(1, num_runs));
        uint32_t* fill = md_temp_alloc_array(temp, uint32_t, sys->chain.count);
        MEMCPY(fill, count, sizeof(uint32_t) * sys->chain.count);
        for (size_t r = 0; r < num_runs; ++r) sorted[fill[run_chain[r]]++] = runs[r];

        md_array_resize(sys->chain.run_offset, sys->chain.count + 1, alloc);
        sys->chain.run_offset[0] = 0;
        for (size_t c = 0; c < sys->chain.count; ++c) {
            for (uint32_t r = count[c]; r < count[c + 1]; ++r) {
                const md_urange_t run = sorted[r];
                const size_t n = md_array_size(sys->chain.run);
                if (r > count[c] && sys->chain.run[n - 1].end == run.beg) {
                    sys->chain.run[n - 1].end = run.end;
                } else {
                    md_array_push(sys->chain.run, run, alloc);
                }
                for (uint32_t k = run.beg; k < run.end; ++k) sys->chain.res_chain[k] = (md_chain_idx_t)c;
            }
            sys->chain.run_offset[c + 1] = (uint32_t)md_array_size(sys->chain.run);
        }
    }

    md_temp_end(temp);
    return true;
}

md_particle_type_idx_t md_util_particle_type_find_or_add(md_system_t* sys, md_hashmap32_t* map, str_t name, md_atomic_number_t z) {
    ASSERT(sys && map);
    // The name as the type keeps it, its first 8 characters
    const uint64_t kv[2] = { md_label_key(make_label(name)), (uint64_t)z };
    const uint64_t key = md_hash64(kv, sizeof(kv), 0);
    const uint32_t* found = md_hashmap_get(map, key);
    if (found) return (md_particle_type_idx_t)*found;
    const md_particle_type_idx_t type = md_particle_type_find_or_add(&sys->atom.type, name, z, md_atomic_number_mass(z), md_atomic_number_vdw_radius(z),
        md_atomic_number_cpk_color(z), 0, sys->alloc);
    md_hashmap_add(map, key, (uint32_t)type);
    return type;
}

static inline bool residue_kind_is_monomer(md_residue_kind_t kind) {
    return kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
}

// Whether a nucleotide is DNA (1) or RNA (2), 0 when it cannot be told. The name decides when it is one of the
// standard names (DA, DC, ... and A, C, ...). Otherwise (ADE, CYT, ... of CHARMM, which are both) the sugar does: a
// resolved nucleotide with a 2' oxygen is RNA, one without DNA.
static int nucleotide_type(const md_system_t* sys, size_t res_idx) {
    const str_t name = md_system_residue_name(sys, res_idx);
    if (md_util_resname_dna(name)) return 1;
    if (md_util_resname_rna(name)) return 2;
    if (md_system_residue_flags(sys, res_idx) & MD_RESIDUE_FLAG_RESOLVED) {
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const str_t atom = md_atom_name(&sys->atom, i);
            if (str_eq(atom, STR_LIT("O2'")) || str_eq(atom, STR_LIT("O2*"))) return 2;
        }
        return 1;
    }
    return 0;
}

// What a part made of the residues in res_range is, from what its residues are. Two or more amino acids
// (nucleotides) which make up at least half of the residues are a peptide (nucleic acid), which lets the caps of a
// chain (ACE, NME) and modified residues the names do not cover come along.
static md_part_kind_t classify_part(const md_system_t* sys, md_urange_t res_range) {
    const size_t n = res_range.end - res_range.beg;
    if (n == 0) return MD_PART_KIND_UNKNOWN;
    if (n == 1) {
        return md_system_residue_kind(sys, res_range.beg) == MD_RESIDUE_KIND_WATER ? MD_PART_KIND_WATER : MD_PART_KIND_NON_POLYMER;
    }

    size_t num_amino = 0, num_nucl = 0, num_dna = 0, num_rna = 0;
    for (uint32_t c = res_range.beg; c < res_range.end; ++c) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, c);
        if (kind == MD_RESIDUE_KIND_AMINO_ACID) {
            num_amino += 1;
        } else if (kind == MD_RESIDUE_KIND_NUCLEOTIDE) {
            num_nucl += 1;
            const int type = nucleotide_type(sys, c);
            num_dna += type == 1;
            num_rna += type == 2;
        }
    }

    if (num_amino >= 2 && 2 * num_amino >= n) {
        return MD_PART_KIND_PEPTIDE;
    }
    if (num_nucl >= 2 && 2 * num_nucl >= n) {
        if (num_dna == num_nucl) return MD_PART_KIND_DNA;
        if (num_rna == num_nucl) return MD_PART_KIND_RNA;
        return MD_PART_KIND_NUCLEIC;
    }
    return MD_PART_KIND_NON_POLYMER;
}

// The given part types are kept, the inferred ones (which come after them) dropped: the parts are about to infer them anew
static void drop_inferred_part_types(md_part_type_data_t* type, md_allocator_i* alloc) {
    size_t n = type->count;
    while (n > 0 && (type->flags[n - 1] & MD_PART_TYPE_FLAG_INFERRED)) --n;
    md_part_type_data_shrink(type, n, alloc);
}

// The type of a part which is not all of an instance of a given type: one per kind and sequence of residue names. That
// the type was inferred is carried by its flags, the description only says what it is: the kind of a polymer, the
// residue names of anything else.
static md_part_type_idx_t inferred_part_type(md_system_t* sys, md_hashmap32_t* type_of_key, md_urange_t residues) {
    const md_part_kind_t kind = classify_part(sys, residues);
    const bool polymer = md_part_kind_is_polymer(kind);

    uint64_t key = md_hash64(&kind, sizeof(kind), 0);
    for (uint32_t c = residues.beg; c < residues.end; ++c) {
        key = md_hash64_str(md_system_residue_name(sys, c), key);
    }
    const uint32_t* found = md_hashmap_get(type_of_key, key);
    if (found) return (md_part_type_idx_t)*found;

    char buf[64];
    str_t desc = {0};
    if (polymer || kind == MD_PART_KIND_WATER) {
        desc = str_from_cstr(md_part_kind_name(kind));
    } else if (residues.end - residues.beg <= 4) {
        int len = 0;
        for (uint32_t c = residues.beg; c < residues.end; ++c) {
            const str_t name = md_system_residue_name(sys, c);
            len += md_snprintf(buf + len, sizeof(buf) - len, "%s" STR_FMT, c > residues.beg ? "-" : "", STR_ARG(name));
            len = MIN(len, (int)sizeof(buf) - 1);
        }
        desc = (str_t){ buf, (size_t)len };
    } else {
        desc = md_system_residue_name(sys, residues.beg);
    }
    const md_part_type_idx_t type = md_part_type_data_push(&sys->part.type, STR_LIT(""), md_part_type_flags_set_kind(MD_PART_TYPE_FLAG_INFERRED, kind), desc, sys->alloc);
    md_hashmap_add(type_of_key, key, (uint32_t)type);
    return type;
}

bool md_util_system_infer_parts(md_system_t* sys) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;
    md_array_free(sys->part.res_offset, alloc);
    md_array_free(sys->part.type_idx, alloc);
    sys->part.res_offset = NULL;
    sys->part.type_idx = NULL;
    sys->part.count = 0;
    drop_inferred_part_types(&sys->part.type, alloc);
    md_system_topology_changed(sys);

    const size_t num_res = sys->residue.count;
    if (num_res == 0) return true;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The bonds between residues, but those of a metal coordinated by residues: the residue bonded to the one before
    // it (by an atom of a backbone: a cap, a modified residue in a strand, rather than a glycan on a side chain), and
    // the last residue before each which it is bonded to
    uint64_t* bonded_to_prev = make_bitfield(num_res, temp_arena);
    uint64_t* backbone_to_prev = make_bitfield(num_res, temp_arena);
    int32_t* last_bonded_before = md_temp_alloc_array(temp, int32_t, num_res);
    for (size_t c = 0; c < num_res; ++c) last_bonded_before[c] = -1;
    if (sys->bond.count > 0) {
        md_residue_idx_t* atom_res = md_temp_alloc_array(temp, md_residue_idx_t, MAX(1, sys->atom.count));
        MEMSET(atom_res, -1, sizeof(md_residue_idx_t) * sys->atom.count);
        for (size_t c = 0; c < num_res; ++c) {
            const md_urange_t r = md_system_residue_atom_range(sys, c);
            for (uint32_t a = r.beg; a < r.end && a < sys->atom.count; ++a) atom_res[a] = (md_residue_idx_t)c;
        }
        for (size_t b = 0; b < sys->bond.count; ++b) {
            if (sys->bond.flags && (!md_bond_is_graph(sys->bond.flags[b]) || md_bond_kind(sys->bond.flags[b]) == MD_BOND_KIND_COORDINATE)) continue;
            const md_atom_pair_t pair = sys->bond.pairs[b];
            if ((size_t)pair.idx[0] >= sys->atom.count || (size_t)pair.idx[1] >= sys->atom.count) continue;
            const md_residue_idx_t x = atom_res[pair.idx[0]];
            const md_residue_idx_t y = atom_res[pair.idx[1]];
            if (x < 0 || y < 0 || x == y) continue;
            const md_residue_idx_t lo = MIN(x, y), hi = MAX(x, y);
            if (hi - lo == 1) {
                bitfield_set_bit(bonded_to_prev, hi);
                if ((md_atom_flags(&sys->atom, pair.idx[0]) | md_atom_flags(&sys->atom, pair.idx[1])) & MD_ATOM_FLAG_BACKBONE) {
                    bitfield_set_bit(backbone_to_prev, hi);
                }
            }
            last_bonded_before[hi] = MAX(last_bonded_before[hi], lo);
        }
    }

    // The spans the parts are found in: the instances, and the residues before the first and after the last, which no
    // instance holds (all of them in a system without instances). A span which is an instance the data names (an id or
    // an author id) lets a strand go on over a gap while the numbering does.
    typedef struct {
        md_urange_t residues;
        md_instance_idx_t inst;     // -1 for residues in no instance
        bool named;
    } span_t;
    md_array(span_t) spans = 0;
    {
        uint32_t covered = 0;
        for (size_t i = 0; i < sys->instance.count; ++i) {
            md_urange_t r = md_system_instance_residue_range(sys, i);
            r.end = MIN(r.end, (uint32_t)num_res);
            if (r.beg < covered) r.beg = covered;
            if (r.beg >= r.end) continue;
            if (covered < r.beg) {
                const span_t none = {{covered, r.beg}, -1, false};
                md_array_push(spans, none, temp_arena);
            }
            const bool named = !str_empty(md_system_instance_id(sys, i)) || !str_empty(md_system_instance_auth_id(sys, i));
            const span_t span = {r, (md_instance_idx_t)i, named};
            md_array_push(spans, span, temp_arena);
            covered = r.end;
        }
        if (covered < num_res) {
            const span_t none = {{covered, (uint32_t)num_res}, -1, false};
            md_array_push(spans, none, temp_arena);
        }
    }

    md_hashmap32_t type_of_key = { .allocator = temp_arena };

    md_array_push(sys->part.res_offset, 0, alloc);
    for (size_t s = 0; s < md_array_size(spans); ++s) {
        const span_t span = spans[s];
        const md_urange_t residues = span.residues;
        const size_t first_part = md_array_size(sys->part.res_offset) - 1;

        uint32_t beg = residues.beg;
        bool run_polymer = residue_kind_is_monomer(md_system_residue_kind(sys, beg));
        for (uint32_t j = residues.beg + 1; j <= residues.end; ++j) {
            bool goes_on = false;
            if (j < residues.end) {
                const md_residue_kind_t kj = md_system_residue_kind(sys, j);
                const md_residue_kind_t kp = md_system_residue_kind(sys, j - 1);
                const bool bonded = bitfield_test_bit(bonded_to_prev, j);
                if (kj == MD_RESIDUE_KIND_WATER || kj == MD_RESIDUE_KIND_ION || kp == MD_RESIDUE_KIND_WATER || kp == MD_RESIDUE_KIND_ION) {
                    goes_on = false;
                } else if (residue_kind_is_monomer(kj) && kj == kp) {
                    // A strand, over a gap unless the chain ends there: within an instance the data names while the
                    // numbering goes on, elsewhere only in a system without bonds while it runs on by one
                    if (bonded) {
                        goes_on = true;
                    } else if (!polymer_chain_ends_between(sys, j - 1, j)) {
                        const md_sequence_id_t seq_prev = md_system_residue_seq_id(sys, j - 1);
                        const md_sequence_id_t seq_j    = md_system_residue_seq_id(sys, j);
                        goes_on = span.named ? seq_j > seq_prev : (sys->bond.count == 0 && seq_j == seq_prev + 1);
                    }
                } else if (!residue_kind_is_monomer(kj) && !run_polymer) {
                    // Of no polymer: on while bonded to a residue of the run (the branches of a glycan)
                    goes_on = bonded || last_bonded_before[j] >= (int32_t)beg;
                } else {
                    // A cap or a modified residue bonded into a strand by its backbone, not a glycan on a side chain
                    goes_on = bitfield_test_bit(backbone_to_prev, j);
                }
            }
            if (goes_on) {
                run_polymer |= residue_kind_is_monomer(md_system_residue_kind(sys, j));
                continue;
            }
            md_array_push(sys->part.res_offset, j, alloc);
            beg = j;
            if (j < residues.end) run_polymer = residue_kind_is_monomer(md_system_residue_kind(sys, j));
        }

        // The types: that of the instance where the part is all of it, else an inferred one
        const size_t end_part = md_array_size(sys->part.res_offset) - 1;
        const md_part_type_idx_t given = span.inst >= 0 ? md_system_instance_part_type_idx(sys, (size_t)span.inst) : -1;
        for (size_t p = first_part; p < end_part; ++p) {
            const md_urange_t r = { sys->part.res_offset[p], sys->part.res_offset[p + 1] };
            md_part_type_idx_t type = -1;
            if (given >= 0 && (size_t)given < sys->part.type.count && end_part - first_part == 1) {
                type = given;
                if (md_system_part_type_kind(sys, (size_t)type) == MD_PART_KIND_UNKNOWN) {
                    // A given type which does not say what it is (a polymer entity without _entity_poly)
                    sys->part.type.flags[type] = md_part_type_flags_set_kind(sys->part.type.flags[type], classify_part(sys, r));
                }
            } else {
                type = inferred_part_type(sys, &type_of_key, r);
            }
            md_array_push(sys->part.type_idx, type, alloc);
        }
    }
    sys->part.count = md_array_size(sys->part.res_offset) - 1;
    ASSERT(md_array_size(sys->part.type_idx) == sys->part.count);

    md_temp_end(temp);
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

static inline bool has_ring(const md_ring_data_t* ring_data, const int* ring, int ring_size) {
    const int64_t num_rings = md_ring_data_count(ring_data);
    for (int64_t i = 0; i < num_rings; ++i) {
        const int* ring_i     = md_ring_atoms(ring_data, i);
        const int ring_i_size = (int)md_ring_size(ring_data, i);
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
        md_temp_scope_t temp = md_temp_begin();
        int* buf = md_temp_alloc_array(temp, int, n);
        int j = 0;
        for (int i = lowest_idx; i < n; ++i) {
            buf[j++] = arr[i];
        }
        for (int i = 0; i < lowest_idx; ++i) {
            buf[j++] = arr[i];
        }
        MEMCPY(arr, buf, sizeof(*arr) * n);
        md_temp_end(temp);
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

// Rings of up to 7 atoms: the 7 membered aromatic rings (tropolone, azulene, tropylium) count, larger ones are loops
#define MIN_RING_SIZE MD_RING_MIN_SIZE
#define MAX_RING_SIZE MD_RING_MAX_SIZE
#define MAX_DEPTH (MAX_RING_SIZE/2 + 1)

// Removes the rings which are sums of smaller ones: rings whose bonds are those of strictly smaller rings, each bond
// counted an odd number of times (the 7 membered ring around a cyclopropane fused to a cyclohexane, the 6 membered one
// around the two 5 membered rings of norbornane, the 4 membered one around the two 3 membered rings of
// bicyclo[1.1.0]butane, the 6 membered rings around two or three faces of cubane). What remains are the relevant rings
// of the graph up to the size limit: the smallest set of smallest rings, made symmetric (each of the three 6 membered
// rings of bicyclo[2.2.2]octane is the sum of the other two, which are not smaller, and all three are kept).
// The smaller rings considered are those sharing an atom with the ring, which every ring of such a sum does when the
// rings are this small.
#define ENVELOPE_MAX_CAND 64
#define ENVELOPE_MAX_BONDS ((ENVELOPE_MAX_CAND + 1) * MAX_RING_SIZE)
#define ENVELOPE_WORDS ((ENVELOPE_MAX_BONDS + 63) / 64)
static void remove_envelope_rings(md_ring_data_t* rings, const md_bond_data_t* bond, size_t num_atoms, md_allocator_i* temp_arena) {
    const size_t num_rings = md_ring_data_count(rings);
    if (num_rings < 2) return;

    // The rings of each atom
    uint32_t* off = md_alloc(temp_arena, (num_atoms + 1) * sizeof(uint32_t));
    MEMSET(off, 0, (num_atoms + 1) * sizeof(uint32_t));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        for (size_t k = 0; k < n; ++k) off[ring[k] + 1] += 1;
    }
    for (size_t i = 0; i < num_atoms; ++i) off[i + 1] += off[i];
    uint32_t* idx = md_alloc(temp_arena, MAX(off[num_atoms], 1) * sizeof(uint32_t));
    uint32_t* fill = md_alloc(temp_arena, num_atoms * sizeof(uint32_t));
    MEMCPY(fill, off, num_atoms * sizeof(uint32_t));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        for (size_t k = 0; k < n; ++k) idx[fill[ring[k]]++] = (uint32_t)r;
    }

    uint8_t* removed = md_alloc(temp_arena, num_rings);
    MEMSET(removed, 0, num_rings);
    size_t removed_count = 0;

    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const int n = (int)md_ring_size(rings, r);
        if (n < 4) continue;    // Two rings of 3 make one of 4 at the least

        // Smaller rings which share an atom with this one
        uint32_t cand[ENVELOPE_MAX_CAND];
        int num_cand = 0;
        for (int k = 0; k < n; ++k) {
            for (uint32_t t = off[ring[k]]; t < off[ring[k] + 1]; ++t) {
                const uint32_t s = idx[t];
                if (s == r || (int)md_ring_size(rings, s) >= n) continue;
                bool dup = false;
                for (int c = 0; c < num_cand; ++c) dup |= cand[c] == s;
                if (!dup && num_cand < ENVELOPE_MAX_CAND) cand[num_cand++] = s;
            }
        }
        if (num_cand < 2) continue;

        // The bonds, numbered locally, as bit vectors over GF(2): the ring is a sum of the candidates when it reduces to
        // nothing by them
        int32_t bonds[ENVELOPE_MAX_BONDS];
        int num_bonds = 0;
        uint64_t basis[ENVELOPE_MAX_BONDS][ENVELOPE_WORDS];   // basis[p]: the row whose lowest bit is p
        uint8_t  has_row[ENVELOPE_MAX_BONDS];
        MEMSET(has_row, 0, sizeof(has_row));
        bool valid = true;
        bool envelope = false;
        for (int c = 0; c <= num_cand && valid; ++c) {
            const bool self = c == num_cand;
            const uint32_t s = self ? (uint32_t)r : cand[c];
            const int32_t* rs = md_ring_atoms(rings, s);
            const int ns = (int)md_ring_size(rings, s);
            uint64_t v[ENVELOPE_WORDS] = {0};
            for (int k = 0; k < ns; ++k) {
                const md_bond_idx_t b = md_bond_find(bond, rs[k], rs[(k + 1) % ns]);
                if (b < 0) { valid = false; break; }
                int j = 0;
                while (j < num_bonds && bonds[j] != b) ++j;
                if (j == num_bonds) bonds[num_bonds++] = b;
                v[j >> 6] ^= 1ull << (j & 63);
            }
            if (!valid) break;
            // Reduce by the rows so far
            int p = -1;
            for (;;) {
                p = -1;
                for (int w = 0; w < ENVELOPE_WORDS; ++w) {
                    if (v[w]) { p = w * 64 + (int)ctz64(v[w]); break; }
                }
                if (p < 0 || !has_row[p]) break;
                for (int w = 0; w < ENVELOPE_WORDS; ++w) v[w] ^= basis[p][w];
            }
            if (self) {
                envelope = p < 0;
            } else if (p >= 0) {
                MEMCPY(basis[p], v, sizeof(v));
                has_row[p] = 1;
            }
        }
        if (valid && envelope) {
            removed[r] = 1;
            removed_count += 1;
        }
    }
    if (!removed_count) return;

    size_t kept = 0;
    for (size_t r = 0; r < num_rings; ++r) {
        if (!removed[r]) rings->rings[kept++] = rings->rings[r];
    }
    md_array_shrink(rings->rings, kept);
}
#undef ENVELOPE_MAX_CAND
#undef ENVELOPE_MAX_BONDS
#undef ENVELOPE_WORDS

// Inspired by molstars implementation
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/structure/unit/rings/compute.ts#L249
// Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
// @author David Sehnal <david.sehnal@gmail.com>

// The rings of the bond graph (bond, with its connectivity), pushed to out: see md_util_system_infer_rings. The atoms
// set in skip (optional) are neither roots nor passed through.
static void find_rings(md_ring_data_t* out, const md_bond_data_t* bond, size_t num_atoms, const uint64_t* skip, md_allocator_i* alloc, md_temp_scope_t temp) {
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // Some typedefs to enable laboration of types
    typedef uint16_t color_t;
    typedef uint8_t  depth_t;
    typedef uint16_t mark_t;
    typedef  int32_t pred_t;
        
    color_t* color = md_temp_alloc_zero_array(temp, color_t, num_atoms);
    depth_t* depth = md_temp_alloc_zero_array(temp, depth_t, num_atoms);
    mark_t*  mark  = md_temp_alloc_zero_array(temp, mark_t, num_atoms);
    pred_t*  pred  = md_temp_alloc_zero_array(temp, pred_t, num_atoms);
    MEMSET(pred, -1, num_atoms * sizeof(pred_t));    // We can do memset as the representation of -1 under two's complement is 0xFFFFFFFF

    // The capacity is arbitrary here, but will be resized if needed.
    md_fifo_t queue = md_fifo_create(64, temp_arena);

    md_hashset_t ring_set = {.allocator = temp_arena};
    
    color_t current_color = 1;
    mark_t  current_mark  = 1;

    for (int atom_idx = 0; atom_idx < (int)num_atoms; ++atom_idx) {
        if (skip && bitfield_test_bit(skip, atom_idx)) continue;

        // Skip any atom that has already been colored in the previous search
        if (color[atom_idx] == current_color) continue;

        // Set i as our root
        depth[atom_idx] = 1;
        pred [atom_idx] = -1;
        mark [atom_idx] = current_mark++;

        md_fifo_clear(&queue);
        md_fifo_push(&queue, atom_idx);
        while (!md_fifo_empty(&queue)) {
            int idx = md_fifo_pop(&queue);

            md_bond_iter_t it = md_bond_iter(bond, idx);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);

                //if (next > idx) continue;
                if (next == pred[idx]) continue;  // avoid adding parent to search queue
                if (skip && bitfield_test_bit(skip, next)) continue;

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
                            md_ring_data_push(out, ring, len, MD_RING_FLAG_NONE, alloc);
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
                    md_fifo_push(&queue, next);
                }
            }
        }
    }

    remove_envelope_rings(out, bond, num_atoms, temp_arena);
}

void md_util_infer_rings(md_ring_data_t* out_rings, const md_bond_data_t* bond, size_t num_atoms, md_allocator_i* alloc) {
    ASSERT(out_rings && alloc);
    ASSERT(bond);
    md_ring_data_clear(out_rings);
    if (num_atoms == 0 || bond->count == 0 || !bond->conn.offset) return;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    find_rings(out_rings, bond, num_atoms, NULL, alloc, temp);
    md_temp_end(temp);
}

bool md_util_system_infer_rings(md_system_t* sys) {
    ASSERT(sys);
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;

    size_t num_atoms = md_system_atom_count(sys);
    size_t num_bonds = md_system_bond_count(sys);

    if (num_atoms == 0) {
        return false;
    }
    if (num_bonds == 0) {
        return false;
    }
    md_system_topology_changed(sys);

    md_ring_data_clear(&sys->ring);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // Water and ions have no rings, and neither do coarse grained beads or virtual sites: a ring is made of atoms.
    // Martini closes a few beads into a rigid triangle where the atoms had a ring (one bead per two or three ring
    // atoms, fused rings a triangulated patch), and the elastic network and secondary structure bonds of a coarse
    // grained protein close thousands of cycles which are no rings at all. Skipped atoms are neither roots nor
    // passed through, so a ring of atoms in a molecule which also holds beads is still found.
    uint64_t* skip = make_bitfield(num_atoms, temp_arena);
    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        if (kind == MD_RESIDUE_KIND_WATER || kind == MD_RESIDUE_KIND_ION) {
            const md_urange_t range = md_system_residue_atom_range(sys, ci);
            for (uint32_t i = range.beg; i < range.end && i < num_atoms; ++i) bitfield_set_bit(skip, i);
        }
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        if (md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_ATOM) bitfield_set_bit(skip, i);
    }

    find_rings(&sys->ring, &sys->bond, num_atoms, skip, alloc, temp);

    md_temp_end(temp);

    md_util_system_infer_ring_flags(sys);

    return true;
}

// Is ring (n atoms in order) the ring of a monosaccharide (MD_RING_FLAG_SUGAR). Returns the position in the ring of
// its anomeric C, -1 if it is not one, and the position of the ring O in out_ring_o. The anomeric C is the neighbor of
// the ring O with an O or N outside of the ring (C1 of an aldose, C2 of a ketose); with both neighbors carrying one, the
// one with more, and the first of the two in the ring when they carry as many.
static int sugar_anomeric(const md_system_t* sys, const int32_t* ring, int n, int* out_ring_o) {
    if (n != 5 && n != 6) return -1;
    int ring_o = -1;
    for (int k = 0; k < n; ++k) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, ring[k]);
        if (z == MD_Z_O) {
            if (ring_o >= 0) return -1;
            ring_o = k;
        } else if (z != MD_Z_C) {
            return -1;
        }
    }
    if (ring_o < 0) return -1;

    int substituted = 0;
    int sub[MD_RING_MAX_SIZE] = {0};    // O and N outside of the ring, on each ring C
    for (int k = 0; k < n; ++k) {
        if (k == ring_o) continue;
        md_bond_iter_t it = md_bond_iter(&sys->bond, ring[k]);
        while (md_bond_iter_has_next(&it)) {
            const md_atom_idx_t x = md_bond_iter_atom_index(&it);
            const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            md_bond_iter_next(&it);
            // No double or aromatic bond on the ring atoms, within the ring or out of it
            if (md_bond_order(f) >= 2 || (f & MD_BOND_FLAG_AROMATIC)) return -1;
            bool in_ring = false;
            for (int j = 0; j < n; ++j) in_ring |= ring[j] == x;
            if (in_ring) continue;
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, x);
            sub[k] += z == MD_Z_O || z == MD_Z_N;
        }
        substituted += sub[k] > 0;
    }
    const int next = (ring_o + 1) % n;
    const int prev = (ring_o + n - 1) % n;
    if (substituted < 2 || (!sub[next] && !sub[prev])) return -1;
    *out_ring_o = ring_o;
    if (sub[next] != sub[prev]) return sub[next] > sub[prev] ? next : prev;
    return MIN(next, prev);
}

// Puts the atoms of a sugar ring in their canonical order: the ring O first, then the anomeric C and on around the
// ring away from the O (O5 C1 C2 C3 C4 C5 of a pyranose, O4' C1' C2' C3' C4' of the furanose of a nucleotide)
static void sugar_canonical_order(int32_t* ring, int n, int ring_o, int anomeric) {
    int32_t tmp[MD_RING_MAX_SIZE];
    const int step = anomeric == (ring_o + 1) % n ? 1 : n - 1;
    for (int k = 0; k < n; ++k) tmp[k] = ring[(ring_o + k * step) % n];
    MEMCPY(ring, tmp, n * sizeof(int32_t));
}

void md_util_system_infer_ring_flags(md_system_t* sys) {
    ASSERT(sys);
    const size_t num_rings = md_ring_data_count(&sys->ring);
    if (num_rings == 0) return;
    for (size_t r = 0; r < num_rings; ++r) sys->ring.rings[r].flags = 0;
    if (!sys->bond.conn.offset || !sys->bond.flags) return;

    // Bonds in more than one ring make those rings fused
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    uint8_t* bond_rings = md_temp_alloc(temp, MAX(sys->bond.count, 1));
    MEMSET(bond_rings, 0, MAX(sys->bond.count, 1));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        for (size_t k = 0; k < n; ++k) {
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b >= 0 && bond_rings[b] < 255) bond_rings[b] += 1;
        }
    }

    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        md_ring_flags_t flags = MD_RING_FLAG_NONE;
        bool aromatic = n > 0;
        for (size_t k = 0; k < n; ++k) {
            if (md_atom_atomic_number(&sys->atom, ring[k]) != MD_Z_C) flags |= MD_RING_FLAG_HETERO;
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b < 0) {
                aromatic = false;
                continue;
            }
            aromatic &= (sys->bond.flags[b] & MD_BOND_FLAG_AROMATIC) != 0;
            if (bond_rings[b] > 1) flags |= MD_RING_FLAG_FUSED;
        }
        if (aromatic) flags |= MD_RING_FLAG_AROMATIC;
        int ring_o = -1;
        const int anomeric = sugar_anomeric(sys, ring, (int)n, &ring_o);
        if (anomeric >= 0) {
            flags |= MD_RING_FLAG_SUGAR;
            sugar_canonical_order(md_ring_atoms_mut(&sys->ring, r), (int)n, ring_o, anomeric);
        }
        sys->ring.rings[r].flags = (uint16_t)flags;
    }

    // Bonds of aromatic rings fuse the other rings they are in to an aromatic ring: an aromatic ring which shares a
    // bond with another is of a fused aromatic system (naphthalene), one fused only to rings which are not aromatic
    // keeps a sextet of its own (tetralin)
    uint8_t* bond_aromatic = md_temp_alloc(temp, MAX(sys->bond.count, 1));
    MEMSET(bond_aromatic, 0, MAX(sys->bond.count, 1));
    for (size_t r = 0; r < num_rings; ++r) {
        if (!(md_ring_flags(&sys->ring, r) & MD_RING_FLAG_AROMATIC)) continue;
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        for (size_t k = 0; k < n; ++k) {
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b >= 0 && bond_aromatic[b] < 255) bond_aromatic[b] += 1;
        }
    }
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        // Aromatic rings other than this one
        const uint8_t self = (md_ring_flags(&sys->ring, r) & MD_RING_FLAG_AROMATIC) ? 1 : 0;
        for (size_t k = 0; k < n; ++k) {
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b >= 0 && bond_aromatic[b] > self) {
                sys->ring.rings[r].flags |= MD_RING_FLAG_FUSED_AROMATIC;
                break;
            }
        }
    }
    md_temp_end(temp);
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE
#undef MAX_DEPTH

// The molecules (md_molecule_data_t) are the connected components of the bond graph, joined where the bonds do not
// hold a molecule together (see there).
// Breadth first traversal of the connected component containing 'seed'.
//
// Marks every atom it reaches in 'visited' and never clears it; the caller supplies a bitfield which
// is clear over the component. When out_atoms is supplied the atoms are written in BFS order and
// out_count receives how many there were. When out_pred is supplied it records, per global atom
// index, the atom each one was reached from, with the seed recording itself.
//
// Returns the last atom popped from the queue. BFS pops in non decreasing distance order, so that
// atom is at maximum edge distance from the seed - which is what the double sweep below needs.
static int bfs_traverse(int32_t* out_atoms, size_t* out_count, int32_t* out_pred, uint64_t* visited, md_fifo_t* queue, const md_bond_data_t* bond, int seed) {
    md_fifo_clear(queue);

    bitfield_set_bit(visited, seed);
    if (out_pred) out_pred[seed] = seed;
    md_fifo_push(queue, seed);

    size_t count = 0;
    int last = seed;

    while (!md_fifo_empty(queue)) {
        const int cur = md_fifo_pop(queue);
        last = cur;
        if (out_atoms) out_atoms[count] = cur;
        count += 1;

        md_bond_iter_t it = md_bond_iter(bond, cur);
        while (md_bond_iter_has_next(&it)) {
            const int next = md_bond_iter_atom_index(&it);
            md_bond_iter_next(&it);

            if (bitfield_test_bit(visited, next)) continue;

            bitfield_set_bit(visited, next);
            if (out_pred) out_pred[next] = cur;
            md_fifo_push(queue, next);
        }
    }

    if (out_count) *out_count = count;
    return last;
}

static int32_t uf_find(int32_t* parent, int32_t x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

// Returns true if a and b were in different sets
static bool uf_union(int32_t* parent, int32_t a, int32_t b) {
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a == b) return false;
    if (a < b) parent[b] = a; else parent[a] = b;
    return true;
}

// ### ELASTIC BONDS ###
//
// The chemistry of a coarse grained chain is that of infer_coarse_grained_bonds: between two of its residues there is the
// link from the last backbone bead of one to the first of the next, and a disulfide between the SC1 of two cysteines. A
// chain goes on while the residues are of one kind and their numbering does not start over.
static inline bool cg_chain_restarts(md_sequence_id_t prev_seq, md_sequence_id_t seq) {
    return seq < prev_seq && !(seq <= 1 && prev_seq >= 9999);  // GRO and PDB wrap at 99999 and 9999
}

static inline bool cg_cysteine_sc1(const md_system_t* sys, int32_t res_idx, int32_t atom_idx) {
    const str_t res_name = md_system_residue_name(sys, res_idx);
    return (str_eq_cstr(res_name, "CYS") || str_eq_cstr(res_name, "CYX")) && str_eq_cstr(md_atom_name(&sys->atom, atom_idx), "SC1");
}

size_t md_util_system_infer_elastic_bonds(md_system_t* sys) {
    ASSERT(sys);
    const size_t N = sys->atom.count;
    const size_t C = sys->residue.count;
    if (sys->bond.count == 0 || !sys->bond.flags || N == 0 || C == 0 || !sys->residue.atom_offset) return 0;
    if (!md_system_is_coarse_grained(sys)) return 0;

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* arena = md_temp_allocator(temp);

    int32_t* res_of = md_alloc(arena, sizeof(int32_t) * N);
    for (size_t i = 0; i < N; ++i) res_of[i] = -1;
    // The beads a residue of a chain links by: the first and the last bead of its backbone, and in a nucleotide of
    // Martini 3 its BB2, which is linked to the next residue as well
    int32_t* first_bb = md_alloc(arena, sizeof(int32_t) * C);
    int32_t* last_bb  = md_alloc(arena, sizeof(int32_t) * C);
    int32_t* also_bb  = md_alloc(arena, sizeof(int32_t) * C);
    const uint32_t version = martini_system(sys, label_of_atom, &sys->atom).version;
    const bool m3 = (version & (MARTINI_2 | MARTINI_2P)) == 0;
    for (size_t ci = 0; ci < C; ++ci) {
        first_bb[ci] = last_bb[ci] = also_bb[ci] = -1;
        const md_urange_t range = md_system_residue_atom_range(sys, ci);
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        const bool biopolymer = kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
        for (uint32_t i = range.beg; i < range.end && i < N; ++i) {
            res_of[i] = (int32_t)ci;
            if (!biopolymer || md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_BEAD) continue;
            if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE)) continue;
            if (first_bb[ci] == -1) first_bb[ci] = (int32_t)i;
            last_bb[ci] = (int32_t)i;
            if (m3 && kind == MD_RESIDUE_KIND_NUCLEOTIDE && str_eq_cstr(md_atom_name(&sys->atom, i), "BB2")) also_bb[ci] = (int32_t)i;
        }
    }

    // The links of the chains first: they say where a chain begins and ends, which the bond closing a cyclic one needs
    uint8_t* has_prev = md_alloc(arena, C);
    uint8_t* has_next = md_alloc(arena, C);
    MEMSET(has_prev, 0, C);
    MEMSET(has_next, 0, C);
    int32_t* chain = md_alloc(arena, sizeof(int32_t) * C);
    for (size_t ci = 0; ci < C; ++ci) chain[ci] = (int32_t)ci;
    md_array(uint32_t) candidates = 0;

    for (size_t b = 0; b < sys->bond.count; ++b) {
        if (!md_bond_is_graph(sys->bond.flags[b])) continue;
        int32_t x = sys->bond.pairs[b].idx[0];
        int32_t y = sys->bond.pairs[b].idx[1];
        if (x < 0 || y < 0 || (size_t)x >= N || (size_t)y >= N) continue;
        if (res_of[x] > res_of[y]) { const int32_t t = x; x = y; y = t; }
        const int32_t ca = res_of[x];
        const int32_t cb = res_of[y];
        if (ca < 0 || ca == cb) continue;  // Within a residue
        if (md_atom_particle_kind(&sys->atom, x) != MD_PARTICLE_BEAD || md_atom_particle_kind(&sys->atom, y) != MD_PARTICLE_BEAD) continue;
        const md_residue_kind_t ka = md_system_residue_kind(sys, ca);
        const md_residue_kind_t kb = md_system_residue_kind(sys, cb);
        if (!(ka == MD_RESIDUE_KIND_AMINO_ACID || ka == MD_RESIDUE_KIND_NUCLEOTIDE) || ka != kb) continue;

        const bool link = cb == ca + 1 && (x == last_bb[ca] || x == also_bb[ca]) && y == first_bb[cb] &&
                          !cg_chain_restarts(md_system_residue_seq_id(sys, ca), md_system_residue_seq_id(sys, cb));
        if (link) {
            has_next[ca] = 1;
            has_prev[cb] = 1;
            uf_union(chain, ca, cb);
        } else if (ka == MD_RESIDUE_KIND_AMINO_ACID && cg_cysteine_sc1(sys, ca, x) && cg_cysteine_sc1(sys, cb, y)) {
            // A disulfide
        } else {
            md_array_push(candidates, (uint32_t)b, arena);
        }
    }

    size_t count = 0;
    for (size_t k = 0; k < md_array_size(candidates); ++k) {
        const uint32_t b = candidates[k];
        int32_t x = sys->bond.pairs[b].idx[0];
        int32_t y = sys->bond.pairs[b].idx[1];
        if (res_of[x] > res_of[y]) { const int32_t t = x; x = y; y = t; }
        const int32_t ca = res_of[x];
        const int32_t cb = res_of[y];
        // The bond which closes a cyclic chain: from the first bead of the backbone at its start to the last at its end
        const bool closure = cb > ca + 1 && x == first_bb[ca] && y == last_bb[cb] && !has_prev[ca] && !has_next[cb] && uf_find(chain, ca) == uf_find(chain, cb);
        if (closure) continue;
        // A spring has neither an order nor resonance
        sys->bond.flags[b] = md_bond_flags_set_kind(sys->bond.flags[b] & ~(MD_BOND_FLAG_ORDER_MASK | MD_BOND_FLAG_RESONANCE_MASK), MD_BOND_KIND_ELASTIC);
        count += 1;
    }

    md_temp_end(temp);
    if (count) {
        md_bond_build_connectivity(&sys->bond, N, sys->alloc);
        md_system_topology_changed(sys);
    }
    return count;
}

// Whether a bond holds a molecule together (md_molecule_data_t): a bond of the graph, but not a coordination bond
// between two residues, a metal and what surrounds it (the Mg of a GTP with its phosphates and the waters of its shell,
// the Zn of a zinc finger with its cysteines), which are molecules of their own. Unless a topology gives it: a bond
// within one of its molecules, as GROMACS has it. A coordination bond within a residue (the Fe of a heme) holds it.
static bool bond_holds_molecule(const md_system_t* sys, size_t b) {
    if (!sys->bond.flags) return true;
    const md_bond_flags_t flags = sys->bond.flags[b];
    if (!md_bond_is_graph(flags)) return false;
    if (md_bond_kind(flags) != MD_BOND_KIND_COORDINATE || md_bond_origin(flags) == MD_BOND_ORIGIN_TOPOLOGY) return true;
    if (!sys->residue.count || !sys->residue.atom_offset) return false;
    const md_atom_pair_t pair = sys->bond.pairs[b];
    return md_system_residue_find_by_atom_idx(sys, (size_t)pair.idx[0]) == md_system_residue_find_by_atom_idx(sys, (size_t)pair.idx[1]);
}

// The molecules of a topology are whole (md_instance_data_t::mol_type_idx): the atoms of each of its instances, joined
// to the first where nothing else joins them (the two strands of a duplex which martinize-dna writes as one molecule
// type, held together by an elastic network which is no bond of the graph)
static void link_topology_molecules(const md_system_t* sys, int32_t* set, md_array(md_atom_pair_t)* links, md_allocator_i* alloc) {
    if (!sys->instance.mol_type_idx) return;
    const size_t atom_count = md_system_atom_count(sys);
    for (size_t inst = 0; inst < sys->instance.count; ++inst) {
        if (md_system_instance_molecule_type_idx(sys, inst) < 0) continue;
        md_urange_t r = md_system_instance_atom_range(sys, inst);
        if (r.end > atom_count) r.end = (uint32_t)atom_count;
        for (uint32_t a = r.beg + 1; a < r.end; ++a) {
            if (uf_union(set, (int32_t)r.beg, (int32_t)a)) {
                md_atom_pair_t pair = {{ (int32_t)r.beg, (int32_t)a }};
                md_array_push(*links, pair, alloc);
            }
        }
    }
}

// Hierarchy links for coarse grained systems, used where bonds leave the topology disconnected.
//
// md_molecule_data_t is a spanning forest and unwrap only needs each atom's parent to stay within half a cell
// of it, not to be bonded to it. Coarse grained systems frequently have no bonds (md_util_system_infer skips
// bond inference for them) or have a distance guess that misses stretched ones. So the missing links are
// synthesized from the residue hierarchy instead, WITHOUT looking at coordinates:
//   - within a residue, every bead hangs off the residue's anchor (its backbone bead, else its first bead)
//   - the anchors of consecutive polymer residues (amino acid or nucleotide, same kind, consecutive sequence
//     ids, or of one part whatever their ids when the parts are known) are linked along the chain
// A link is only emitted between atoms that are not already connected, by bonds or by an earlier link, so a
// fully bonded system gets none and the result stays a forest over bonds + links.
//
// Deliberately NOT done: distance based linking. In a dense coarse grained membrane or solvent that percolates
// into one box spanning structure, and unwrapping that is meaningless. Non polymer residues (lipids, sterols)
// are not chained together either - a span of residues may hold hundreds of them in file order, not spatial order.
//
// The links are returned as pairs and never written to sys->bond: they are not bonds, must not be drawn and
// must not feed ring detection.
static md_array(md_atom_pair_t) synthesize_hierarchy_links(const md_system_t* sys, md_allocator_i* alloc) {
    md_array(md_atom_pair_t) links = 0;

    const size_t atom_count = md_system_atom_count(sys);
    const size_t res_count = sys->residue.count;
    if (atom_count == 0 || res_count == 0) {
        return links;
    }

    int32_t* set = md_alloc(alloc, sizeof(int32_t) * atom_count);
    for (size_t i = 0; i < atom_count; ++i) set[i] = (int32_t)i;

    for (size_t i = 0; i < sys->bond.count; ++i) {
        if (!bond_holds_molecule(sys, i)) continue;
        const md_atom_pair_t pair = sys->bond.pairs[i];
        if ((size_t)pair.idx[0] < atom_count && (size_t)pair.idx[1] < atom_count) {
            uf_union(set, pair.idx[0], pair.idx[1]);
        }
    }

    // Part of each residue, -1 when the parts are unknown (md_util_system_infer infers them before the molecules)
    int32_t* res_part = md_alloc(alloc, sizeof(int32_t) * res_count);
    for (size_t i = 0; i < res_count; ++i) res_part[i] = -1;
    for (size_t p = 0; p < sys->part.count; ++p) {
        md_urange_t range = md_system_part_residue_range(sys, p);
        for (uint32_t ci = range.beg; ci < range.end && ci < res_count; ++ci) {
            res_part[ci] = (int32_t)p;
        }
    }

    int32_t             prev_anchor = -1;
    size_t              prev_res    = 0;
    md_residue_kind_t prev_kind     = MD_RESIDUE_KIND_OTHER;

    for (size_t ci = 0; ci < res_count; ++ci) {
        md_urange_t range = md_system_residue_atom_range(sys, ci);
        if (range.end > atom_count) range.end = (uint32_t)atom_count;
        if (range.beg >= range.end) {
            prev_anchor = -1;
            continue;
        }

        int32_t anchor = (int32_t)range.beg;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) {
                anchor = (int32_t)i;
                break;
            }
        }

        for (uint32_t i = range.beg; i < range.end; ++i) {
            if ((int32_t)i != anchor && uf_union(set, anchor, (int32_t)i)) {
                md_atom_pair_t pair = {{ anchor, (int32_t)i }};
                md_array_push(links, pair, alloc);
            }
        }

        const md_residue_kind_t ck = md_system_residue_kind(sys, ci);
        const md_residue_kind_t kind = residue_kind_is_monomer(ck) ? ck : MD_RESIDUE_KIND_OTHER;
        // Within a part the residues are of one strand whatever their numbering (a strand goes on over a gap,
        // md_util_system_infer_parts), as a molecule holds its parts whole; without parts by consecutive numbers
        const bool same_part = res_part[prev_res] >= 0 && res_part[prev_res] == res_part[ci];
        const bool consecutive = !sys->residue.seq_id || md_system_residue_seq_id(sys, ci) == md_system_residue_seq_id(sys, prev_res) + 1;
        if (kind && prev_anchor >= 0 && prev_res + 1 == ci && prev_kind == kind &&
            res_part[prev_res] == res_part[ci] && (same_part || consecutive))
        {
            if (uf_union(set, prev_anchor, anchor)) {
                md_atom_pair_t pair = {{ prev_anchor, anchor }};
                md_array_push(links, pair, alloc);
            }
        }

        prev_anchor = kind ? anchor : -1;
        prev_res    = ci;
        prev_kind   = kind;
    }

    link_topology_molecules(sys, set, &links, alloc);
    return links;
}

// The atoms which join two residues of a polymer chain over a gap: the C and N of a peptide bond, the O3' and P of a
// phosphodiester bond (the O5' where the second has no phosphate: the 5' end of a strand, a phosphate not modelled),
// else the last atom of the first and the first of the second. Never an atom the residues do not have.
static md_atom_pair_t chain_gap_link(const md_system_t* sys, size_t res_idx) {
    const md_urange_t range_a = md_system_residue_atom_range(sys, res_idx);
    const md_urange_t range_b = md_system_residue_atom_range(sys, res_idx + 1);
    md_atom_pair_t pair = {{ (int32_t)range_a.end - 1, (int32_t)range_b.beg }};
    md_amino_acid_atoms_t aa_a, aa_b;
    md_nucleic_acid_atoms_t na_a, na_b;
    if (md_util_amino_acid_atoms_extract(&aa_a, &sys->atom, range_a) && md_util_amino_acid_atoms_extract(&aa_b, &sys->atom, range_b)) {
        if (aa_a.c >= 0) pair.idx[0] = aa_a.c;
        if (aa_b.n >= 0) pair.idx[1] = aa_b.n;
    } else if (md_util_nucleic_acid_atoms_extract(&na_a, &sys->atom, range_a) && md_util_nucleic_acid_atoms_extract(&na_b, &sys->atom, range_b)) {
        if (na_a.o3 >= 0) pair.idx[0] = na_a.o3;
        if (na_b.p >= 0) pair.idx[1] = na_b.p;
        else if (na_b.o5 >= 0) pair.idx[1] = na_b.o5;
    }
    return pair;
}

// The links which make a molecule whole where its bonds leave it in pieces (see md_molecule_data_t). The atoms of a
// residue are of one molecule: an atom which no bond reaches (the virtual site of a TIP4P water, an ion a file names
// with the residue it sits in) is linked to the first heavy atom of its residue. The residues of a polymer chain are of
// one molecule: over a gap (residues missing from a crystal structure, a chain of a topology whose bonds are not
// given) its two ends are linked. Between the pieces only: a residue or a chain its bonds join gets no link.
static md_array(md_atom_pair_t) synthesize_molecule_links(const md_system_t* sys, md_allocator_i* alloc) {
    md_array(md_atom_pair_t) links = 0;
    const size_t atom_count = md_system_atom_count(sys);
    if (atom_count == 0 || sys->residue.count == 0) return links;

    int32_t* set = md_alloc(alloc, sizeof(int32_t) * atom_count);
    for (size_t i = 0; i < atom_count; ++i) set[i] = (int32_t)i;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        if (!bond_holds_molecule(sys, i)) continue;
        const md_atom_pair_t pair = sys->bond.pairs[i];
        if ((size_t)pair.idx[0] < atom_count && (size_t)pair.idx[1] < atom_count) {
            uf_union(set, pair.idx[0], pair.idx[1]);
        }
    }

    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        md_urange_t range = md_system_residue_atom_range(sys, ci);
        if (range.end > atom_count) range.end = (uint32_t)atom_count;
        if (range.beg >= range.end) continue;
        int32_t anchor = -1;
        for (uint32_t i = range.beg; i < range.end && anchor < 0; ++i) {
            if (md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_VIRTUAL_SITE && md_atom_atomic_number(&sys->atom, i) > MD_Z_H) anchor = (int32_t)i;
        }
        if (anchor < 0) anchor = (int32_t)range.beg;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (uf_union(set, anchor, (int32_t)i)) {
                md_atom_pair_t pair = {{ anchor, (int32_t)i }};
                md_array_push(links, pair, alloc);
            }
        }
    }

    for (size_t p = 0; p < sys->part.count; ++p) {
        if (!md_part_kind_is_polymer(md_system_part_kind(sys, p))) continue;
        const md_urange_t residues = md_system_part_residue_range(sys, p);
        for (uint32_t ci = residues.beg; ci + 1 < residues.end; ++ci) {
            const md_urange_t ra = md_system_residue_atom_range(sys, ci);
            const md_urange_t rb = md_system_residue_atom_range(sys, ci + 1);
            if (ra.beg >= ra.end || rb.beg >= rb.end || rb.end > atom_count) continue;
            if (uf_find(set, (int32_t)ra.beg) == uf_find(set, (int32_t)rb.beg)) continue;
            const md_atom_pair_t pair = chain_gap_link(sys, ci);
            ASSERT(pair.idx[0] >= 0 && pair.idx[1] >= 0);
            uf_union(set, pair.idx[0], pair.idx[1]);
            md_array_push(links, pair, alloc);
        }
    }
    link_topology_molecules(sys, set, &links, alloc);
    return links;
}

// The name of an inferred molecule type: what its parts are ("peptide + NAG-NAG"), or its formula where it has no parts
// (an XYZ file), C and H first and the other elements by symbol
static void molecule_type_name(char* buf, size_t cap, const md_system_t* sys, const md_part_idx_t* parts, size_t num_parts, const md_molecule_t* mol) {
    ASSERT(cap > 0);
    buf[0] = '\0';
    int len = 0;
    if (num_parts > 0) {
        const size_t MAX_SHOWN = 3;
        for (size_t i = 0; i < num_parts && i < MAX_SHOWN; ++i) {
            const str_t desc = md_system_part_type_description(sys, (size_t)md_system_part_type_idx(sys, (size_t)parts[i]));
            len += md_snprintf(buf + len, cap - len, "%s" STR_FMT, i > 0 ? " + " : "", STR_ARG(desc));
            len = MIN(len, (int)cap - 1);
        }
        if (num_parts > MAX_SHOWN) {
            len += md_snprintf(buf + len, cap - len, " + %zu more", num_parts - MAX_SHOWN);
        }
        return;
    }

    uint32_t count[MD_Z_Count] = {0};
    for (size_t i = 0; i < mol->count; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, (size_t)mol->atom_idx[i]);
        if (z < MD_Z_Count) count[z] += 1;
    }
    const md_atomic_number_t first[2] = {MD_Z_C, MD_Z_H};
    for (int k = 0; k < 2; ++k) {
        if (count[MD_Z_C] == 0 && first[k] == MD_Z_H) break;    // Without carbon, H goes in order (Hill)
        const md_atomic_number_t z = first[k];
        if (!count[z]) continue;
        const str_t sym = md_util_element_symbol(z);
        len += count[z] > 1 ? md_snprintf(buf + len, cap - len, STR_FMT "%u", STR_ARG(sym), count[z]) : md_snprintf(buf + len, cap - len, STR_FMT, STR_ARG(sym));
        len = MIN(len, (int)cap - 1);
        count[z] = 0;
    }
    // The rest by symbol: few elements, a selection each time is cheap
    for (;;) {
        md_atomic_number_t best = 0;
        for (int z = 1; z < MD_Z_Count; ++z) {
            if (count[z] && (!best || str_cmp_lex(md_util_element_symbol((md_atomic_number_t)z), md_util_element_symbol(best)) < 0)) best = (md_atomic_number_t)z;
        }
        if (!best) break;
        const str_t sym = md_util_element_symbol(best);
        len += count[best] > 1 ? md_snprintf(buf + len, cap - len, STR_FMT "%u", STR_ARG(sym), count[best]) : md_snprintf(buf + len, cap - len, STR_FMT, STR_ARG(sym));
        len = MIN(len, (int)cap - 1);
        count[best] = 0;
    }
    if (len == 0) {
        md_snprintf(buf, cap, "%zu particles", mol->count);
    }
}

// The types of the molecules (md_molecule_type_data_t). The molecule which is all of an instance of a given type (a
// molecule of a topology) is of that type. The others are of inferred types, one for the molecules made alike: of the
// same part types in the same order, as many atoms of each particle type, and as many bonds between each pair of
// particle types. The inferred types of the last inference are dropped first.
static void infer_molecule_types(md_system_t* sys) {
    md_allocator_i* alloc = sys->alloc;
    md_molecule_data_t* mol = &sys->molecule;

    {
        md_molecule_type_data_t* type = &mol->type;
        size_t n = type->count;
        while (n > 0 && (type->flags[n - 1] & MD_MOLECULE_TYPE_FLAG_INFERRED)) --n;
        md_molecule_type_data_shrink(type, n, alloc);
    }
    md_array_resize(mol->type_idx, mol->count, alloc);
    if (mol->count == 0) return;

    const size_t atom_count = sys->atom.count;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    int32_t* mol_of_atom = md_temp_alloc_array(temp, int32_t, MAX(1, atom_count));
    for (size_t m = 0; m < mol->count; ++m) {
        for (uint32_t slot = mol->offset[m]; slot < mol->offset[m + 1]; ++slot) {
            mol_of_atom[mol->atom_idx[slot]] = (int32_t)m;
        }
    }

    // What each molecule is made of: its part types in order (hashed in sequence), and its particle types and bonds
    // (summed, as the order of its atoms and bonds is the file's)
    uint64_t* seq_key  = md_temp_alloc_zero_array(temp, uint64_t, mol->count);
    uint64_t* atom_sum = md_temp_alloc_zero_array(temp, uint64_t, mol->count);
    uint64_t* bond_sum = md_temp_alloc_zero_array(temp, uint64_t, mol->count);
    uint32_t* first_part = md_temp_alloc_array(temp, uint32_t, mol->count);
    uint32_t* num_parts  = md_temp_alloc_zero_array(temp, uint32_t, mol->count);
    for (size_t m = 0; m < mol->count; ++m) first_part[m] = UINT32_MAX;

    for (size_t p = 0; p < sys->part.count; ++p) {
        const md_urange_t atoms = md_system_part_atom_range(sys, p);
        if (atoms.beg >= atoms.end || atoms.beg >= atom_count) continue;
        const int32_t m = mol_of_atom[atoms.beg];
        const md_part_type_idx_t t = md_system_part_type_idx(sys, p);
        seq_key[m] = md_hash64(&t, sizeof(t), seq_key[m] + 1);
        if (first_part[m] == UINT32_MAX) first_part[m] = (uint32_t)p;
        num_parts[m] += 1;
    }
    for (size_t a = 0; a < atom_count; ++a) {
        const md_particle_type_idx_t t = sys->atom.type_idx ? sys->atom.type_idx[a] : 0;
        atom_sum[mol_of_atom[a]] += md_hash64(&t, sizeof(t), 0x9E3779B97F4A7C15ull);
    }
    for (size_t b = 0; b < sys->bond.count; ++b) {
        if (!bond_holds_molecule(sys, b)) continue;
        const md_atom_pair_t pair = sys->bond.pairs[b];
        if ((size_t)pair.idx[0] >= atom_count || (size_t)pair.idx[1] >= atom_count) continue;
        const md_particle_type_idx_t ta = sys->atom.type_idx ? sys->atom.type_idx[pair.idx[0]] : 0;
        const md_particle_type_idx_t tb = sys->atom.type_idx ? sys->atom.type_idx[pair.idx[1]] : 0;
        const uint32_t key = ((uint32_t)MIN(ta, tb) << 16) | (uint32_t)MAX(ta, tb);
        bond_sum[mol_of_atom[pair.idx[0]]] += md_hash64(&key, sizeof(key), 0xC2B2AE3D27D4EB4Full);
    }

    md_hashmap32_t type_of_key = { .allocator = temp_arena };
    md_array(md_part_idx_t) parts = 0;
    char name[128];

    for (size_t m = 0; m < mol->count; ++m) {
        const md_molecule_t molecule = md_system_molecule(sys, m);

        // All of one instance of a given type
        md_molecule_type_idx_t type = -1;
        if (sys->instance.mol_type_idx && molecule.count > 0) {
            int32_t lo = INT32_MAX, hi = -1;
            for (size_t i = 0; i < molecule.count; ++i) {
                lo = MIN(lo, molecule.atom_idx[i]);
                hi = MAX(hi, molecule.atom_idx[i]);
            }
            const md_instance_idx_t inst = md_system_instance_find_by_atom_idx(sys, (size_t)lo);
            if (inst >= 0) {
                const md_urange_t r = md_system_instance_atom_range(sys, (size_t)inst);
                const md_molecule_type_idx_t given = md_system_instance_molecule_type_idx(sys, (size_t)inst);
                if (given >= 0 && (size_t)given < mol->type.count && r.beg == (uint32_t)lo && r.end == (uint32_t)hi + 1 && r.end - r.beg == molecule.count) {
                    type = given;
                }
            }
        }

        if (type < 0) {
            const size_t n = molecule.count;
            const uint64_t key = md_hash64_combine(md_hash64(&n, sizeof(n), seq_key[m]), md_hash64_combine(atom_sum[m], bond_sum[m]));
            const uint32_t* found = md_hashmap_get(&type_of_key, key);
            if (found) {
                type = (md_molecule_type_idx_t)*found;
            } else {
                md_array_shrink(parts, 0);
                for (uint32_t p = first_part[m]; p < sys->part.count && md_array_size(parts) < num_parts[m]; ++p) {
                    const md_urange_t atoms = md_system_part_atom_range(sys, p);
                    if (atoms.beg < atoms.end && atoms.beg < atom_count && mol_of_atom[atoms.beg] == (int32_t)m) {
                        md_array_push(parts, (md_part_idx_t)p, temp_arena);
                    }
                }
                molecule_type_name(name, sizeof(name), sys, parts, md_array_size(parts), &molecule);
                type = md_molecule_type_data_push(&mol->type, str_from_cstr(name), MD_MOLECULE_TYPE_FLAG_INFERRED, alloc);
                md_hashmap_add(&type_of_key, key, (uint32_t)type);
            }
        }
        mol->type_idx[m] = type;
    }

    md_temp_end(temp);
}

bool md_util_system_infer_molecules(md_system_t* sys) {
    ASSERT(sys);
    ASSERT(sys->alloc);

    md_allocator_i* alloc = sys->alloc;
    md_system_topology_changed(sys);

    md_array_shrink(sys->molecule.offset, 0);
    md_array_shrink(sys->molecule.atom_idx, 0);
    md_array_shrink(sys->molecule.parent_idx, 0);
    md_array_shrink(sys->molecule.atom_slot, 0);
    sys->molecule.count = 0;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    const size_t atom_count = md_system_atom_count(sys);

    // Reverse map from atom to slot, one entry per atom in the system
    md_array_resize(sys->molecule.atom_slot, atom_count, alloc);

    // Which atoms have been assigned to a molecule. Persists across components.
    uint64_t* visited = make_bitfield(atom_count, temp_arena);

    // Reused by the sweeps below. Cleared over the atoms of the current component only, so the cost
    // stays proportional to the component rather than to the system.
    uint64_t* sweep_visited = make_bitfield(atom_count, temp_arena);

    // Atoms of the residue currently being processed, and the predecessor of each atom as recorded
    // by the second sweep. Only the entries belonging to the current component are meaningful.
    int32_t* component = md_temp_alloc_zero_array(temp, int32_t, atom_count);
    int32_t* pred      = md_temp_alloc_zero_array(temp, int32_t, atom_count);

    // Parallel queues:
    // - atom_queue   : global atom index to visit
    // - parent_queue : global index of the atom it was reached from (the root is its own parent)
    md_fifo_t atom_queue   = md_fifo_create(1024, temp_arena);
    md_fifo_t parent_queue = md_fifo_create(1024, temp_arena);

    // The traversal runs over the bonds which hold molecules together (bond_holds_molecule: not the coordination of a
    // metal by other residues), plus hierarchy links for coarse grained systems where bonds leave beads disconnected
    // (see synthesize_hierarchy_links), and otherwise the links which keep a residue and a polymer chain in one molecule
    // where no bond does (see synthesize_molecule_links). The links only shape the forest; sys->bond is untouched.
    const md_bond_data_t* graph = &sys->bond;
    md_bond_data_t link_graph = {0};
    {
        md_array(md_atom_pair_t) links = md_system_is_coarse_grained(sys) ? synthesize_hierarchy_links(sys, temp_arena) : synthesize_molecule_links(sys, temp_arena);
        const size_t num_links = md_array_size(links);
        // The bonds which hold molecules together (bond_holds_molecule): the coordination bonds between residues left out
        size_t num_held = 0;
        for (size_t b = 0; b < sys->bond.count; ++b) num_held += bond_holds_molecule(sys, b);
        if (num_links > 0 || num_held < sys->bond.count) {
            link_graph.count = num_held + num_links;
            md_array_resize(link_graph.pairs, link_graph.count, temp_arena);
            md_array_resize(link_graph.flags, link_graph.count, temp_arena);
            MEMSET(link_graph.flags, 0, sizeof(md_bond_flags_t) * link_graph.count);
            size_t n = 0;
            for (size_t b = 0; b < sys->bond.count; ++b) {
                if (!bond_holds_molecule(sys, b)) continue;
                link_graph.pairs[n] = sys->bond.pairs[b];
                if (sys->bond.flags) link_graph.flags[n] = sys->bond.flags[b];
                n += 1;
            }
            if (num_links) MEMCPY(link_graph.pairs + n, links, sizeof(md_atom_pair_t) * num_links);
            md_bond_build_connectivity(&link_graph, atom_count, temp_arena);
            graph = &link_graph;
        }
    }

    // Sentinel start offset
    md_array_push(sys->molecule.offset, 0, alloc);

    for (int i = 0; i < (int)atom_count; ++i) {
        if (bitfield_test_bit(visited, i)) continue;

        // Sweep 1: discover the component from an arbitrary seed and note an atom at maximum distance
        // from it. This is the sweep that claims the atoms in 'visited'; the atom list it produces is
        // what lets the following sweeps clear their bits without touching the rest of the system.
        size_t component_size = 0;
        const int far_a = bfs_traverse(component, &component_size, NULL, visited, &atom_queue, graph, i);

        // Locate the graph center: the atom whose greatest distance to any other atom in the
        // component is smallest. Standard double sweep - the midpoint of a longest path. Exact for
        // acyclic structures, a close approximation in the presence of rings.
        // Components of two atoms or fewer are their own center.
        int root = i;
        if (component_size > 2) {
            // Sweep 2: from far_a, find the opposite end of the path and the chain leading back
            for (size_t k = 0; k < component_size; ++k) bitfield_clear_bit(sweep_visited, component[k]);
            const int far_b = bfs_traverse(NULL, NULL, pred, sweep_visited, &atom_queue, graph, far_a);

            int length = 0;
            for (int a = far_b; a != pred[a]; a = pred[a]) length += 1;

            root = far_b;
            for (int k = 0; k < length / 2; ++k) root = pred[root];
        }

        // Sweep 3: emit the component in BFS order from the root, which establishes the
        // slot(parent) < slot(child) invariant md_molecule_data_t documents.
        for (size_t k = 0; k < component_size; ++k) bitfield_clear_bit(sweep_visited, component[k]);

        md_fifo_clear(&atom_queue);
        md_fifo_clear(&parent_queue);

        bitfield_set_bit(sweep_visited, root);
        md_fifo_push(&atom_queue, root);
        md_fifo_push(&parent_queue, root);   // the root is its own parent

        while (!md_fifo_empty(&atom_queue)) {
            const int cur        = md_fifo_pop(&atom_queue);
            const int parent_idx = md_fifo_pop(&parent_queue);

            sys->molecule.atom_slot[cur] = (int32_t)md_array_size(sys->molecule.atom_idx);
            md_array_push(sys->molecule.atom_idx, cur, alloc);
            md_array_push(sys->molecule.parent_idx, parent_idx, alloc);

            md_bond_iter_t it = md_bond_iter(graph, cur);
            while (md_bond_iter_has_next(&it)) {
                const int next = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);

                if (bitfield_test_bit(sweep_visited, next)) continue;

                bitfield_set_bit(sweep_visited, next);
                md_fifo_push(&atom_queue, next);
                md_fifo_push(&parent_queue, cur);
            }
        }

        sys->molecule.count += 1;
        md_array_push(sys->molecule.offset, (uint32_t)md_array_size(sys->molecule.atom_idx), alloc);
    }
    
    md_temp_end(temp);

    infer_molecule_types(sys);
    return true;
}

typedef struct {
    const char* res;                // NULL for any component named as its res_kind: an amino acid or a nucleotide
    const char* atom;
    int   atomic_nr;
    float mass;
    float radius;
    md_particle_kind_t  particle;   // What the particle is
    md_atom_flags_t     role;       // Its role in its residue (backbone, side chain)
    md_residue_kind_t res_kind;     // What its residue is, when the particle tells (OTHER when it does not)
} particle_type_t;


// Predefined particle types (This include coarse grained types)
static const particle_type_t predefined_particle_types[] = {
    // Martini proteins: the backbone bead and the side chain beads of any residue named as an amino acid, which covers
    // the names martinize2 gives protonation states and lipidated residues. No atomistic force field names an atom BB
    // or SC#. The mass and radius are those of a regular bead; the residues of martini_db have those of each version.
    // Martini molecules other than amino acids and nucleotides (lipids, water, ions, ...) are told by their residue
    // name and all their bead names at once (martini_db_match_molecule), not here.
    {NULL, "BB",      0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_BACKBONE,   MD_RESIDUE_KIND_AMINO_ACID},
    {NULL, "SC#",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_SIDE_CHAIN, MD_RESIDUE_KIND_AMINO_ACID},
    // The charge dummies of the polar side chains of polarizable Martini 2 (martini22p, elnedyn22p): no Lennard-Jones,
    // so the radius of element 0, as md_tpr gives a particle without
    {NULL, "SCP|SCN", 0, 36.0f, 1.0f,      MD_PARTICLE_BEAD, MD_ATOM_FLAG_SIDE_CHAIN, MD_RESIDUE_KIND_AMINO_ACID},
    // The lipid tails of the lipidated residues of Martini 3: farnesyl, geranylgeranyl and palmitoyl Cys, myristoyl Gly
    {"CYSF|CYSG|CYSP|GLYM", "C#", 0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_SIDE_CHAIN, MD_RESIDUE_KIND_AMINO_ACID},

    // Martini nucleotides (Martini 2 DNA, Uusitalo et al. 2015): the phosphate (BB1) and the sugar (BB2, BB3) of the
    // backbone, and the base (SC1-SC4), of any residue named as a nucleotide
    {NULL, "BB1",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_BACKBONE, MD_RESIDUE_KIND_NUCLEOTIDE},
    {NULL, "BB#",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_NUCLEOSIDE, MD_RESIDUE_KIND_NUCLEOTIDE},
    {NULL, "SC#",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_NUCLEOSIDE, MD_RESIDUE_KIND_NUCLEOTIDE},

    {"RAMP", "PO1", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM1", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM2", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM3", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM4", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM5", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM6", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "PO2", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL1", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL2", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1A", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2A", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3A", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1B", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2B", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3B", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL3", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL4", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1C", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2C", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3C", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1D", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2D", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3D", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL5", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL6", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1E", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2E", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL7", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL8", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1F", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2F", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S01", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S02", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S03", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S04", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S05", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S06", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S07", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S08", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S09", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S10", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S11", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S12", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S13", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S14", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S15", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S16", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S17", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S18", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S19", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S20", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S21", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S22", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S23", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S24", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S25", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S26", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S27", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S28", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S29", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S30", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S31", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S32", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S33", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S34", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S35", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S36", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S37", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S38", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S39", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},

    // The massless charge sites of 4 and 5 site water models: M of TIP4P and its variants and OPC (MW in GROMACS, OM
    // in CHARMM, EPW in AMBER, M in OpenMM) and the two lone pairs of TIP5P (LP1 LP2, EP1 EP2, M1 M2). These are all
    // atom models of a single molecule, NOT coarse grained: the coarse grained flag would turn off covalent bond and
    // hydrogen bond inference for the whole system. No element, no mass, no radius, and no bonds.
    { "SOL|WAT|HOH|H2O|TIP4|TIP5|TIP4P|TIP5P|TP4|TP5|T4P|T5P|T4E|OPC", "MW|M|OM|EPW|EP|EP#|LP#|M#", 0, 0.0f, 0.0f, MD_PARTICLE_VIRTUAL_SITE, 0, MD_RESIDUE_KIND_OTHER },

    { "ION", "CL",  17, 35.45f,     1.8f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
    { "ION", "NA",  11, 22.99f,     2.3f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
    { "ION", "K",   19, 39.10f,     2.7f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
    { "ION", "CA",  20, 40.08f,     2.7f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
	{ "ION", "MG",  12, 24.31f,     2.0f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },

    {"*", "IC", 0,   3240.0f,   12.5f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"*", "OC", 0,   3240.0f,   12.5f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"*", "CC", 0,   3900.0f,   14.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},

    {"P10", "X1", 0, 3.31E6f,   50.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"P20", "X2", 0, 2.65E7f,   100.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"P30", "X3", 0, 8.93E7f,   150.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"P40", "X4", 0, 2.12E8f,   200.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},

    {"C##", "C#?", 0, 1727.7f,   13.5f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
};

static inline uint64_t gen_key_from_names(str_t res_name, str_t atom_name) {
    char buf[64];
    int len = md_snprintf(buf, sizeof(buf), STR_FMT "_" STR_FMT, STR_ARG(res_name), STR_ARG(atom_name));
    str_t str = {buf, (size_t)len};
    return md_hash64_str(str, 0);
}

// Searches the predefined particle types for a matching atom name within a residue
// We support simple pattern matching here like wildcards and numerical substitutes
// * is supported to match any sequence of characters
// ? is supported to match any single character
// # is supported to match any single digit
// | separates alternatives

static inline bool pattern_match(const char* pattern, const char* str) {
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

// Alternatives separated by '|': "SOL|WAT" matches either
static bool pattern_match_any(const char* patterns, const char* str) {
    char buf[128];
    const char* beg = patterns;
    for (;;) {
        const char* end = beg;
        while (*end && *end != '|') ++end;
        const size_t len = (size_t)(end - beg);
        if (len < sizeof(buf)) {
            MEMCPY(buf, beg, len);
            buf[len] = '\0';
            if (pattern_match(buf, str)) return true;
        }
        if (!*end) return false;
        beg = end + 1;
    }
}

// Whether a residue name is that of the kind: what a predefined type without a residue name matches
static bool resname_of_kind(str_t res_name, md_residue_kind_t kind) {
    switch (kind) {
    case MD_RESIDUE_KIND_AMINO_ACID: return md_util_resname_amino_acid(res_name);
    case MD_RESIDUE_KIND_NUCLEOTIDE: return md_util_resname_nucleotide(res_name);
    default: return false;
    }
}

static particle_type_t* find_predefined_particle_type(str_t res_name, str_t atom_name) {
    for (size_t i = 0; i < ARRAY_SIZE(predefined_particle_types); ++i) {
        particle_type_t type = predefined_particle_types[i];
        if (!pattern_match_any(type.atom, atom_name.ptr)) continue;
        const bool res_match = type.res ? pattern_match_any(type.res, res_name.ptr) : resname_of_kind(res_name, type.res_kind);
        if (res_match) {
            return (particle_type_t*)&predefined_particle_types[i];
        }
    }
    return NULL;
}

// The residue kind a predefined type tells, applied when the residue has none yet
static inline void apply_predefined_res_kind(md_system_t* sys, size_t res_idx, md_residue_kind_t kind) {
    if (kind != MD_RESIDUE_KIND_OTHER && sys->residue.flags && md_system_residue_kind(sys, res_idx) == MD_RESIDUE_KIND_OTHER) {
        sys->residue.flags[res_idx] = md_residue_flags_set_kind(sys->residue.flags[res_idx], kind);
    }
}

void md_util_system_infer_particle_types(md_system_t* sys, const str_t atom_labels[]) {
    if (!sys) {
        MD_LOG_ERROR("system was null");
        return;
    }
    if (!atom_labels) {
        MD_LOG_ERROR("Must provide atom_labels");
        return;
    }
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator was null");
        return;
    }
    md_allocator_i* alloc = sys->alloc;
    md_system_topology_changed(sys);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    // Maps (residue name and size, atom label) to the type in the low 16 bits and the predefined type it came from
    // (its index + 1, 0 for none) in the high 16 bits: the predefined type also gives the atom its role.
    md_hashmap32_t particle_type_cache = {.allocator = temp_alloc };
    md_hashmap_reserve(&particle_type_cache, 256);
    STATIC_ASSERT(sizeof(md_particle_type_idx_t) == 2, "The cache packs a type index in 16 bits");
    STATIC_ASSERT(ARRAY_SIZE(predefined_particle_types) < 0xFFFF, "The cache packs a predefined type in 16 bits");

    if (sys->residue.count > 0) {
        // The masses and sizes of Martini beads depend on the version of Martini
        const martini_system_t martini = martini_system(sys, label_of_array, atom_labels);
        // The molecules of the force field of beads of none of its topologies, by name and size: the first and how many
        typedef struct { str_t name; size_t size; size_t first; size_t count; } unknown_variant_t;
        md_array(unknown_variant_t) unknown_variants = 0;
        for (size_t res_idx = 0; res_idx < sys->residue.count; ++res_idx) {
            str_t res_name         = md_system_residue_name(sys, res_idx);
            md_urange_t res_range = md_system_residue_atom_range(sys, res_idx);
            size_t res_size        = res_range.end - res_range.beg;

            // A Martini molecule (a lipid, water, an ion, ...): its residue name and all its bead names tell it, and give
            // every bead its mass, size and particle kind. A residue of a single bead is one only in a Martini system:
            // an ion of Martini 3 (ION NA) is named like an atomistic one.
            uint8_t map[256];
            const martini_db_residue_t* mol = martini_db_match_molecule(martini.version, res_name, atom_labels + res_range.beg, res_size, map);
            if (mol && res_size == 1 && !martini.martini && mol->kind != MD_RESIDUE_KIND_WATER) mol = NULL;
            // A molecule renamed (CHL for CHOL, the three letters of a PDB file) is told by its beads, in a Martini system
            // only: atomistic residues are named as their atoms, which may be the names of beads
            if (!mol && martini.martini && !md_util_resname_amino_acid(res_name) && !md_util_resname_nucleotide(res_name)) {
                mol = martini_db_match_beads(martini.version, atom_labels + res_range.beg, NULL, res_size, map);
            }
            if (mol) {
                for (size_t i = res_range.beg; i < res_range.end; ++i) {
                    if (sys->atom.type_idx[i] != 0) continue;
                    const martini_db_bead_t* bead = martini_db_bead(mol, map[i - res_range.beg]);
                    const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, bead->particle);
                    sys->atom.type_idx[i] = find_or_add_particle_type(&sys->atom.type, atom_labels[i], 0, bead->mass, bead->radius, 0, flags, alloc);
                }
                apply_predefined_res_kind(sys, res_idx, mol->kind);
                continue;
            }
            // A molecule of the force field whose beads are those of none of its topologies (a variant of it, a residue
            // edited), in a Martini system: beads still, rather than the elements their names would read as (the dummy V1
            // of a sterol vanadium, its C1 carbon). Regular beads, as nothing tells their masses and sizes, of the kind of
            // the molecule, and without a topology of it no bonds. Not of a polymer, whose residues are bonded by name.
            if (martini.martini && !md_util_resname_amino_acid(res_name) && !md_util_resname_nucleotide(res_name) && !martini_db_linked(res_name)) {
                const martini_db_residue_t* like = martini_db_unknown_variant(res_name, atom_labels + res_range.beg, res_size);
                if (like) {
                    const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, MD_PARTICLE_BEAD);
                    for (size_t i = res_range.beg; i < res_range.end; ++i) {
                        if (sys->atom.type_idx[i] != 0) continue;
                        sys->atom.type_idx[i] = find_or_add_particle_type(&sys->atom.type, atom_labels[i], 0, 72.0f, MARTINI_R, 0, flags, alloc);
                    }
                    apply_predefined_res_kind(sys, res_idx, like->kind);
                    // Told once per name and size
                    size_t u = 0;
                    while (u < md_array_size(unknown_variants) && !(str_eq(unknown_variants[u].name, res_name) && unknown_variants[u].size == res_size)) ++u;
                    if (u == md_array_size(unknown_variants)) {
                        unknown_variant_t v = { res_name, res_size, res_idx, 0 };
                        md_array_push(unknown_variants, v, temp_alloc);
                    }
                    unknown_variants[u].count += 1;
                    continue;
                }
            }

            // The labels of the residue are part of the key, as which of them a name stands for can depend on the others
            uint64_t res_key = md_hash64_str(res_name, res_size);
            for (size_t i = res_range.beg; i < res_range.end; ++i) res_key = md_hash64_str(atom_labels[i], res_key);

            // What the predefined types of its atoms say the residue is
            md_residue_kind_t res_kind = MD_RESIDUE_KIND_OTHER;
            for (size_t i = res_range.beg; i < res_range.end; ++i) {
                if (sys->atom.type_idx[i] != 0) continue;

                const particle_type_t* predef = NULL;
                uint64_t key = md_hash64_str(atom_labels[i], res_key);
                uint32_t* cached = md_hashmap_get(&particle_type_cache, key);
                if (cached) {
                    sys->atom.type_idx[i] = (md_particle_type_idx_t)(*cached & 0xFFFF);
                    const uint32_t entry = *cached >> 16;
                    predef = entry ? &predefined_particle_types[entry - 1] : NULL;
                } else {
                    str_t atom_name = atom_labels[i];
                    md_atomic_number_t z = 0;
                    float mass = 0;
                    float radius = 0;
                    uint32_t color = 0;
                    md_particle_type_flags_t flags = MD_PARTICLE_TYPE_FLAG_NONE;

                    // Try to find in predefined set
                    predef = find_predefined_particle_type(res_name, atom_name);
                    if (predef) {
                        z       = (md_atomic_number_t)predef->atomic_nr;
                        mass    = predef->mass;
                        radius  = predef->radius;
                        color   = z == 0 ? 0 : md_atomic_number_cpk_color(z);
                        flags   = md_particle_type_flags_set_particle_kind(flags, predef->particle);
                        if (predef->res_kind == MD_RESIDUE_KIND_AMINO_ACID || predef->res_kind == MD_RESIDUE_KIND_NUCLEOTIDE) {
                            // A Martini bead of an amino acid or a nucleotide: its version tells its mass and size
                            const martini_db_bead_t* bead = martini_db_biopolymer_bead(martini.version, res_name, atom_labels + res_range.beg, res_size, i - res_range.beg);
                            if (bead) {
                                mass   = bead->mass;
                                radius = bead->radius;
                                flags  = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, bead->particle);
                            }
                        }
                    } else {
                        z       = md_atomic_number_infer_from_label(atom_name, res_name, res_size);
                        mass    = md_atomic_number_mass(z);
                        radius  = md_atomic_number_vdw_radius(z);
                        color   = md_atomic_number_cpk_color(z);
                    }

                    md_particle_type_idx_t type = find_or_add_particle_type(&sys->atom.type, atom_name, z, mass, radius, color, flags, alloc);
                    sys->atom.type_idx[i] = type;
                    const uint32_t entry = predef ? (uint32_t)(predef - predefined_particle_types) + 1 : 0;
                    md_hashmap_add(&particle_type_cache, key, (uint32_t)type | (entry << 16));
                }
                if (predef) {
                    if (sys->atom.flags) sys->atom.flags[i] |= predef->role;
                    if (predef->res_kind != MD_RESIDUE_KIND_OTHER) res_kind = predef->res_kind;
                }
            }
            apply_predefined_res_kind(sys, res_idx, res_kind);
        }
        for (size_t u = 0; u < md_array_size(unknown_variants); ++u) {
            const unknown_variant_t* v = &unknown_variants[u];
            const md_urange_t range = md_system_residue_atom_range(sys, v->first);
            char beads[256];
            size_t len = 0;
            beads[0] = '\0';
            for (uint32_t i = range.beg; i < range.end && len + 8 < sizeof(beads); ++i) {
                len += (size_t)snprintf(beads + len, sizeof(beads) - len, "%s" STR_FMT, len ? " " : "", STR_ARG(atom_labels[i]));
            }
            MD_LOG_INFO(STR_FMT " of %zu beads (%s), %zu of them: no Martini topology of the name has these beads; shown as regular beads, without bonds",
                STR_ARG(v->name), v->size, beads, v->count);
        }
    } else {
        for (size_t i = 0; i < sys->atom.count; ++i) {
            if (sys->atom.type_idx[i] != 0) continue;

            uint64_t key = md_hash64_str(atom_labels[i], 0);
            uint32_t* cached_type = md_hashmap_get(&particle_type_cache, key);
            if (cached_type) {
                sys->atom.type_idx[i] = (md_particle_type_idx_t)*cached_type;
            } else {
                md_atomic_number_t z = md_atomic_number_infer_from_label(atom_labels[i], (str_t){0}, 0);
                float mass   = md_atomic_number_mass(z);
                float radius = md_atomic_number_vdw_radius(z);
                uint32_t color = md_atomic_number_cpk_color(z);
                md_particle_type_idx_t type = md_particle_type_find_or_add(&sys->atom.type, atom_labels[i], z, mass, radius, color, MD_PARTICLE_TYPE_FLAG_NONE, alloc);
                sys->atom.type_idx[i] = type;
                md_hashmap_add(&particle_type_cache, key, (uint32_t)type);
            }
        }
    }

    md_temp_end(temp);
}

void md_util_system_augment_particle_types(md_system_t* sys) {
    if (!sys || !sys->atom.type_idx || sys->residue.count == 0) {
        return;
    }
    md_system_topology_changed(sys);

    const size_t num_types = sys->atom.type.count;
    for (size_t res_idx = 0; res_idx < sys->residue.count; ++res_idx) {
        const str_t res_name = md_system_residue_name(sys, res_idx);
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        md_residue_kind_t res_kind = MD_RESIDUE_KIND_OTHER;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const md_particle_type_idx_t t = sys->atom.type_idx[i];
            // Only types without an element: an atom with one is not what these tables describe
            if (t <= 0 || (size_t)t >= num_types || sys->atom.type.z[t] != 0) continue;

            const particle_type_t* entry = find_predefined_particle_type(res_name, md_particle_type_name(&sys->atom.type, t));
            if (!entry) continue;

            // The loader's word on what the particle is (a tpr knows its virtual sites) stands over the table's
            if (sys->atom.type.flags && md_particle_type_particle_kind(&sys->atom.type, t) == MD_PARTICLE_ATOM) {
                sys->atom.type.flags[t] = md_particle_type_flags_set_particle_kind(sys->atom.type.flags[t], entry->particle);
            }
            if (sys->atom.flags) sys->atom.flags[i] |= entry->role;
            if (entry->res_kind != MD_RESIDUE_KIND_OTHER) res_kind = entry->res_kind;
        }
        apply_predefined_res_kind(sys, res_idx, res_kind);
    }
}

void md_util_mask_grow_by_bonds(md_bitmask_t* mask, const md_system_t* sys, size_t extent, const md_bitset_t* viable_mask) {
    ASSERT(mask);
    ASSERT(sys);

    if (extent >= 255) {
        MD_LOG_DEBUG("Maximum supported growth extent is 255, the extent will be clamped to 255");
        extent = 255;
    }
    // The initial depth is 1, so we need to add one.
    extent += 1;
    
    const size_t mask_size = md_bitmask_count(mask);
    if (!mask_size) return;

    md_temp_scope_t temp_scope = md_temp_begin_avoid(mask->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    int* indices = md_alloc(temp_alloc, mask_size * sizeof(int));
    const size_t num_indices = md_bitset_iter_extract_indices(indices, mask_size, md_bitmask_iter(mask));
    ASSERT(num_indices == mask_size);

    md_fifo_t queue = md_fifo_create(64, temp_alloc);

    // Zero means 'not reached'. Temp memory is not cleared, so this has to be, or the walk depends on what the
    // allocator handed out last: fresh pages happen to be zero, reused ones are not.
    uint8_t* depth = md_temp_alloc_array(temp_scope, uint8_t, sys->atom.count);
    MEMSET(depth, 0, sys->atom.count * sizeof(uint8_t));

    {
        md_bitset_iter_t it = md_bitmask_iter(mask);
        while (md_bitset_iter_next(&it)) {
            int idx = (int)it.idx;
            depth[idx] = 1;
        }
    }
    
    for (size_t j = 0; j < num_indices; ++j) {
        int i = indices[j];

        md_fifo_clear(&queue);
        md_fifo_push(&queue, i);
        
        while (!md_fifo_empty(&queue)) {
            int idx = md_fifo_pop(&queue);
            md_bitmask_set_bit(mask, (uint32_t)idx);

            md_bond_iter_t it = md_bond_iter(&sys->bond, idx);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (viable_mask && md_bitset_test(*viable_mask, (uint32_t)next) == false) {
                    continue;
                }
                // Only explore the path if it is not already marked and it is within the extent.
                // Or if the current path is 'cheaper'
                const uint8_t dc = depth[idx] + 1;
                const uint8_t dn = depth[next];
                if (dn == 0 || dc < dn) {
                    depth[next] = dc;
                    if (dc <= extent) {
                        md_fifo_push(&queue, next);
                    }
                }
            }
        }
    }

    md_temp_end(temp_scope);
}

// typedef void (*md_spatial_acc_pair_callback_t)(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param);
static inline void spatial_acc_pair_set_bits_callback(const uint32_t* idx_i, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)idx_i;
    (void)ij_dist2;
    md_bitmask_t* mask = (md_bitmask_t*)user_param;
    md_bitmask_set_indices(mask, j_idx, num_pairs);
}

void md_util_mask_grow_by_radius(md_bitmask_t* mask, const md_system_state_t* state, double radius, const md_bitset_t* viable_mask) {
    ASSERT(mask);
    ASSERT(state);
    
    if (radius <= 0.0) return;
    md_temp_scope_t temp = md_temp_begin_avoid(mask->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    size_t num_indices = md_bitmask_count(mask);
    if (num_indices > 0) {
        int32_t* indices = md_temp_alloc(temp, num_indices * sizeof(int32_t));
        num_indices = md_bitset_iter_extract_indices(indices, num_indices, md_bitmask_iter(mask));

        int32_t* viable_indices = 0;
        size_t viable_count = state->num_atoms;
            
        if (viable_mask) {
            md_bitset_t temp_mask = md_bitset_andnot(*viable_mask, md_bitmask_view(mask), temp_arena);
            const size_t count = md_bitset_count(temp_mask);
            viable_indices = md_temp_alloc(temp, count * sizeof(int32_t));
            viable_count = md_bitset_iter_extract_indices(viable_indices, count, md_bitset_iter(temp_mask));
        }

        if (viable_count > 0) {
            md_spatial_acc_t acc = {.alloc = temp_arena};
            md_coord_stream_t stream = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), viable_indices, viable_count);
            md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &stream, .cutoff = radius, .unitcell = &state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX });

            md_coord_stream_t ext_stream = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), indices, num_indices);
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, radius, spatial_acc_pair_set_bits_callback, mask, 0);
        }
    }

    md_temp_end(temp);
}

void md_util_aabb_compute(float out_aabb_min[3], float out_aabb_max[3], const vec3_t* in_xyz, const float* in_r, const int32_t* in_idx, size_t count) {
    ASSERT(out_aabb_min);
    ASSERT(out_aabb_max);
    ASSERT(in_xyz);

    const float* xyz = (const float*)in_xyz;
    vec4_t aabb_min = vec4_set1( FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);
    size_t i = 0;

    if (!in_idx && !in_r) {
        // A minimum and a maximum need no split: float p of a group of eight atoms is component p % 3
        // Scalar up to 32 byte alignment
        const size_t head = xyz_align_peel(xyz, count);
        for (; i < head; ++i) {
            const vec4_t c = vec4_from_vec3(in_xyz[i], 0);
            aabb_min = vec4_min(aabb_min, c);
            aabb_max = vec4_max(aabb_max, c);
        }
        md_xv v_min[XYZ_GROUP_VECS], v_max[XYZ_GROUP_VECS];
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            v_min[j] = md_xv_set1_ps(+FLT_MAX);
            v_max[j] = md_xv_set1_ps(-FLT_MAX);
        }
        for (; i + 8 <= count; i += 8) {
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                const md_xv v = md_xv_loadu_ps(xyz + i * 3 + j * MD_XV_WIDTH);
                v_min[j] = md_xv_min_ps(v_min[j], v);
                v_max[j] = md_xv_max_ps(v_max[j], v);
            }
        }
        float lo[24], hi[24];
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            md_xv_storeu_ps(lo + j * MD_XV_WIDTH, v_min[j]);
            md_xv_storeu_ps(hi + j * MD_XV_WIDTH, v_max[j]);
        }
        for (int p = 0; p < 24; ++p) {
            aabb_min.elem[p % 3] = MIN(aabb_min.elem[p % 3], lo[p]);
            aabb_max.elem[p % 3] = MAX(aabb_max.elem[p % 3], hi[p]);
        }
    } else {
        md_256 vx_min = md_mm256_set1_ps(+FLT_MAX);
        md_256 vy_min = md_mm256_set1_ps(+FLT_MAX);
        md_256 vz_min = md_mm256_set1_ps(+FLT_MAX);
        md_256 vx_max = md_mm256_set1_ps(-FLT_MAX);
        md_256 vy_max = md_mm256_set1_ps(-FLT_MAX);
        md_256 vz_max = md_mm256_set1_ps(-FLT_MAX);

        for (; i + 8 <= count; i += 8) {
            md_256 x, y, z, r;
            if (in_idx) {
                const md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_mm256_i32gather_xyz_ps(&x, &y, &z, xyz, idx);
                r = in_r ? md_mm256_i32gather_ps(in_r, idx, 4) : md_mm256_setzero_ps();
            } else {
                md_mm256_load_xyz_packed_ps(&x, &y, &z, xyz + i * 3);
                r = md_mm256_loadu_ps(in_r + i);
            }
            vx_min = md_mm256_min_ps(vx_min, md_mm256_sub_ps(x, r));
            vy_min = md_mm256_min_ps(vy_min, md_mm256_sub_ps(y, r));
            vz_min = md_mm256_min_ps(vz_min, md_mm256_sub_ps(z, r));
            vx_max = md_mm256_max_ps(vx_max, md_mm256_add_ps(x, r));
            vy_max = md_mm256_max_ps(vy_max, md_mm256_add_ps(y, r));
            vz_max = md_mm256_max_ps(vz_max, md_mm256_add_ps(z, r));
        }
        aabb_min = vec4_min(aabb_min, vec4_set(md_mm256_reduce_min_ps(vx_min), md_mm256_reduce_min_ps(vy_min), md_mm256_reduce_min_ps(vz_min), 0));
        aabb_max = vec4_max(aabb_max, vec4_set(md_mm256_reduce_max_ps(vx_max), md_mm256_reduce_max_ps(vy_max), md_mm256_reduce_max_ps(vz_max), 0));
    }

    // Remainder
    for (; i < count; ++i) {
        const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
        const vec4_t c = vec4_from_vec3(in_xyz[idx], 0);
        const vec4_t r = vec4_set1(in_r ? in_r[idx] : 0.0f);
        aabb_min = vec4_min(aabb_min, vec4_sub(c, r));
        aabb_max = vec4_max(aabb_max, vec4_add(c, r));
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


void md_util_oobb_compute(float out_basis[3][3], float out_ext_min[3], float out_ext_max[3], const vec3_t* in_xyz, const float* in_r, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {
    ASSERT(out_basis);
    ASSERT(out_ext_min);
    ASSERT(out_ext_max);

    if (count == 0) return;

    vec3_t com = md_util_com_compute(in_xyz, NULL, in_idx, count, cell);
    double cov[3][3] = {0};
    if (cell) {
        vec4_t ref = vec4_from_vec3(com, 0);
        if (cell->flags & MD_UNITCELL_ORTHO) {
            vec4_t ext = { 0 };
            md_unitcell_diag_extract_float(ext.elem, cell);

            for (size_t i = 0; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                const vec4_t c = vec4_deperiodize_ortho(vec4_from_vec3(in_xyz[idx], 0.0f), ref, ext);

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
            float A[3][3] = {0};
            md_unitcell_A_extract_float(A, cell);
            for (size_t i = 1; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                vec4_t c = vec4_from_vec3(in_xyz[idx], 0.0f);
                deperiodize_triclinic(c.elem, ref.elem, MD_AS_CONST_MAT3(A));

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
        vec4_t c = vec4_from_vec3(in_xyz[idx], 1.0f);

        vec4_t p = mat4_mul_vec4(Ri, c);
        min_ext = vec4_min(min_ext, vec4_sub1(p, r));
        max_ext = vec4_max(max_ext, vec4_add1(p, r));
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
            vec4_t ext = { 0 };
            md_unitcell_diag_extract_float(ext.elem, cell);

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
            float A[3][3] = {0};
            md_unitcell_A_extract_float(A, cell);
            for (size_t i = 1; i < count; ++i) {
                const int32_t idx = in_idx ? in_idx[i] : (int32_t)i;
                vec4_t c = in_xyzr[idx];
                deperiodize_triclinic(c.elem, ref.elem, MD_AS_CONST_MAT3(A));

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

// Contiguous atoms without an index: the sums run on the packed floats, the weights spread to the
// lanes of their atoms' three components. Groups of eight atoms from the first, and the last few as
// one more group padded with zeros, so how the sums group follows the atoms' indices and count alone,
// never where in memory they lie: the same positions give the same bits from any buffer.
static void com_packed(float out_com[3], const vec3_t* in_xyz, const float* in_w, size_t count) {
    const float* f = (const float*)in_xyz;
    double acc[3] = {0};
    double acc_w = 0;

    md_xv s[XYZ_GROUP_VECS];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) s[j] = md_xv_setzero_ps();
    md_xv sw = md_xv_setzero_ps();
    size_t i = 0;
    // The last count % 8 atoms padded to a group: zero positions and, weighted, zero weights
    float tail_f[24] = {0};
    float tail_w[8]  = {0};
    const size_t body = ROUND_DOWN(count, 8);
    if (body < count) {
        MEMCPY(tail_f, f + body * 3, (count - body) * 3 * sizeof(float));
        if (in_w) MEMCPY(tail_w, in_w + body, (count - body) * sizeof(float));
    }
    for (; i < count; i += 8) {
        const float* g = i < body ? f + i * 3 : tail_f;
        if (in_w) {
            const float* gw = i < body ? in_w + i : tail_w;
            md_xv w[XYZ_GROUP_VECS];
            xyz_lane_spread(w, gw);
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                s[j] = md_xv_fmadd_ps(md_xv_loadu_ps(g + j * MD_XV_WIDTH), w[j], s[j]);
            }
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(gw));
#if MD_XV_WIDTH == 4
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(gw + 4));
#endif
        } else {
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                s[j] = md_xv_add_ps(s[j], md_xv_loadu_ps(g + j * MD_XV_WIDTH));
            }
        }
    }
    acc_w = in_w ? md_xv_reduce_add_ps(sw) : (double)count;
    xyz_lane_reduce_add(acc, s);

    out_com[0] = (float)(acc[0] / acc_w);
    out_com[1] = (float)(acc[1] / acc_w);
    out_com[2] = (float)(acc[2] / acc_w);
}

// The periodic mean in an orthorhombic cell, contiguous atoms without an index. The angle of each
// coordinate is its own axis scaled (M is diagonal), so the sines and cosines are taken on the
// packed floats as they lie and summed per lane; the same angles and the same epilogue as _com_pbc.
static void com_pbc_ortho_packed(float out_com[3], const vec3_t* in_xyz, const float* in_w, size_t count, const float M[3][3], const float I[3][3]) {
    const float* f = (const float*)in_xyz;
    double acc_c[3] = {0};
    double acc_s[3] = {0};
    double acc_w = 0;

    // From the first atom, sixteen or eight at a time, and the last few as one more group of eight padded with zero
    // weights: every atom's sine and cosine come from the same vector sincos, and how the sums group follows the
    // atoms' indices and count alone, never where in memory they lie
    size_t i = 0;

#if defined(__AVX512F__) && defined(__AVX512DQ__)
    {
        // Sixteen atoms, 48 floats: lane k of the j:th vector holds component (16j + k) % 3
        float   pat[48];
        int32_t spread[48];
        for (int p = 0; p < 48; ++p) {
            pat[p]    = M[p % 3][p % 3];
            spread[p] = p / 3;
        }
        __m512  m5[3], c5[3], s5[3];
        __m512i w5_idx[3];
        for (int j = 0; j < 3; ++j) {
            m5[j] = _mm512_loadu_ps(pat + 16 * j);
            w5_idx[j] = _mm512_loadu_si512(spread + 16 * j);
            c5[j] = _mm512_setzero_ps();
            s5[j] = _mm512_setzero_ps();
        }
        __m512 sw5 = _mm512_setzero_ps();
        for (; i + 16 <= count; i += 16) {
            __m512 w16 = _mm512_setzero_ps();
            if (in_w) {
                w16 = _mm512_loadu_ps(in_w + i);
                sw5 = _mm512_add_ps(sw5, w16);
            }
            for (int j = 0; j < 3; ++j) {
                const __m512 t = _mm512_mul_ps(_mm512_loadu_ps(f + i * 3 + j * 16), m5[j]);
                __m512 s, c;
                md_mm512_sincos_ps(t, &s, &c);
                if (in_w) {
                    const __m512 w = _mm512_permutexvar_ps(w5_idx[j], w16);
                    c5[j] = _mm512_fmadd_ps(c, w, c5[j]);
                    s5[j] = _mm512_fmadd_ps(s, w, s5[j]);
                } else {
                    c5[j] = _mm512_add_ps(c5[j], c);
                    s5[j] = _mm512_add_ps(s5[j], s);
                }
            }
        }
        float tc[48], ts[48];
        for (int j = 0; j < 3; ++j) {
            _mm512_storeu_ps(tc + 16 * j, c5[j]);
            _mm512_storeu_ps(ts + 16 * j, s5[j]);
        }
        for (int p = 0; p < 48; ++p) {
            acc_c[p % 3] += tc[p];
            acc_s[p % 3] += ts[p];
        }
        if (in_w) acc_w += _mm512_reduce_add_ps(sw5);
    }
#endif

    md_xv m[XYZ_GROUP_VECS];
    xyz_lane_pattern(m, M[0][0], M[1][1], M[2][2]);
    md_xv vc[XYZ_GROUP_VECS];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) vc[j] = md_xv_setzero_ps();
    md_xv vs[XYZ_GROUP_VECS];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) vs[j] = md_xv_setzero_ps();
    md_xv sw = md_xv_setzero_ps();
    for (; i + 8 <= count; i += 8) {
        md_xv w[XYZ_GROUP_VECS];
        if (in_w) {
            xyz_lane_spread(w, in_w + i);
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(in_w + i));
#if MD_XV_WIDTH == 4
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(in_w + i + 4));
#endif
        }
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            const md_xv t = md_xv_mul_ps(md_xv_loadu_ps(f + i * 3 + j * MD_XV_WIDTH), m[j]);
            md_xv s, c;
            md_xv_sincos_ps(t, &s, &c);
            if (in_w) {
                vc[j] = md_xv_fmadd_ps(c, w[j], vc[j]);
                vs[j] = md_xv_fmadd_ps(s, w[j], vs[j]);
            } else {
                vc[j] = md_xv_add_ps(vc[j], c);
                vs[j] = md_xv_add_ps(vs[j], s);
            }
        }
    }
    if (i < count) {
        // The rest: a weight of 1 (or its own) for each atom and 0 for the padding, so the padding's cos 0 = 1
        // adds nothing; c * 1 + acc is c + acc, the same as the unweighted lanes
        const size_t rest = count - i;
        float tf[24] = {0};
        float tw[8]  = {0};
        MEMCPY(tf, f + i * 3, rest * 3 * sizeof(float));
        for (size_t q = 0; q < rest; ++q) tw[q] = in_w ? in_w[i + q] : 1.0f;
        md_xv w[XYZ_GROUP_VECS];
        xyz_lane_spread(w, tw);
        if (in_w) {
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(tw));
#if MD_XV_WIDTH == 4
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(tw + 4));
#endif
        }
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            const md_xv t = md_xv_mul_ps(md_xv_loadu_ps(tf + j * MD_XV_WIDTH), m[j]);
            md_xv s, c;
            md_xv_sincos_ps(t, &s, &c);
            vc[j] = md_xv_fmadd_ps(c, w[j], vc[j]);
            vs[j] = md_xv_fmadd_ps(s, w[j], vs[j]);
        }
        i = count;
    }
    acc_w += in_w ? md_xv_reduce_add_ps(sw) : (double)count;
    xyz_lane_reduce_add(acc_c, vc);
    xyz_lane_reduce_add(acc_s, vs);

    const double inv_w = 1.0 / acc_w;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = PI;
        const double x = acc_c[j] * inv_w;
        const double y = acc_s[j] * inv_w;
        if (x * x + y * y > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

static void com(float* out_com, const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count) {
    if (!in_idx) {
        com_packed(out_com, in_xyz, in_w, count);
        return;
    }
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
                __m512 x, y, z;
                md_mm512_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);
                __m512 w    = _mm512_i32gather_ps(idx, in_w, 4);

                vx = _mm512_add_ps(vx, _mm512_mul_ps(x, w));
                vy = _mm512_add_ps(vy, _mm512_mul_ps(y, w));
                vz = _mm512_add_ps(vz, _mm512_mul_ps(z, w));
                vw = _mm512_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 x, y, z;
                md_mm512_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);

                vx = _mm512_add_ps(vx, x);
                vy = _mm512_add_ps(vy, y);
                vz = _mm512_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512 x, y, z;
                md_mm512_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));
                __m512 w = _mm512_loadu_ps(in_w + i);

                vx = _mm512_add_ps(vx, _mm512_mul_ps(x, w));
                vy = _mm512_add_ps(vy, _mm512_mul_ps(y, w));
                vz = _mm512_add_ps(vz, _mm512_mul_ps(z, w));
                vw = _mm512_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512 x, y, z;
                md_mm512_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

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
                md_256 x, y, z;
                md_mm256_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);
                md_256 w    = md_mm256_i32gather_ps(in_w, idx, 4);

                vx = md_mm256_add_ps(vx, md_mm256_mul_ps(x, w));
                vy = md_mm256_add_ps(vy, md_mm256_mul_ps(y, w));
                vz = md_mm256_add_ps(vz, md_mm256_mul_ps(z, w));
                vw = md_mm256_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 x, y, z;
                md_mm256_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);

                vx = md_mm256_add_ps(vx, x);
                vy = md_mm256_add_ps(vy, y);
                vz = md_mm256_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256 x, y, z;
                md_mm256_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));
                md_256 w = md_mm256_loadu_ps(in_w + i);

                vx = md_mm256_add_ps(vx, md_mm256_mul_ps(x, w));
                vy = md_mm256_add_ps(vy, md_mm256_mul_ps(y, w));
                vz = md_mm256_add_ps(vz, md_mm256_mul_ps(z, w));
                vw = md_mm256_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256 x, y, z;
                md_mm256_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

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

#else
    md_128 vx = md_mm_setzero_ps();
    md_128 vy = md_mm_setzero_ps();
    md_128 vz = md_mm_setzero_ps();
    md_128 vw = md_mm_setzero_ps();

    const size_t simd_count = ROUND_DOWN(count, 4);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 x, y, z;
                md_mm_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);
                md_128 w    = md_mm_i32gather_ps(in_w, idx, 4);

                vx = md_mm_add_ps(vx, md_mm_mul_ps(x, w));
                vy = md_mm_add_ps(vy, md_mm_mul_ps(y, w));
                vz = md_mm_add_ps(vz, md_mm_mul_ps(z, w));
                vw = md_mm_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 x, y, z;
                md_mm_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);

                vx = md_mm_add_ps(vx, x);
                vy = md_mm_add_ps(vy, y);
                vz = md_mm_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128 x, y, z;
                md_mm_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));
                md_128 w = md_mm_loadu_ps(in_w + i);

                vx = md_mm_add_ps(vx, md_mm_mul_ps(x, w));
                vy = md_mm_add_ps(vy, md_mm_mul_ps(y, w));
                vz = md_mm_add_ps(vz, md_mm_mul_ps(z, w));
                vw = md_mm_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128 x, y, z;
                md_mm_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

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
                acc_x += in_xyz[idx].x * w;
                acc_y += in_xyz[idx].y * w;
                acc_z += in_xyz[idx].z * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                acc_x += in_xyz[idx].x;
                acc_y += in_xyz[idx].y;
                acc_z += in_xyz[idx].z;
            }
            acc_w = (double)count;
        }
    } else {
        if (in_w) {
            for (; i < count; ++i) {
                float w = in_w[i];
                acc_x += in_xyz[i].x * w;
                acc_y += in_xyz[i].y * w;
                acc_z += in_xyz[i].z * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                acc_x += in_xyz[i].x;
                acc_y += in_xyz[i].y;
                acc_z += in_xyz[i].z;
            }
            acc_w = (double)count;
        }
    }

    out_com[0] = (float)(acc_x / acc_w);
    out_com[1] = (float)(acc_y / acc_w);
    out_com[2] = (float)(acc_z / acc_w);
}

// Internal versions for COM computation supporting triclinic and ortho PBC
static void _com_pbc_w(float out_com[3], const vec3_t* in_xyz, const float* in_w, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_xyz);
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
        md_512 v_x, v_y, v_z;
        md_mm512_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        md_512 v_w = _mm512_loadu_ps(in_w + i);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
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
        md_256 v_x, v_y, v_z;
        md_mm256_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        md_256 v_w = md_mm256_loadu_ps(in_w + i);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
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
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128 v_x, v_y, v_z;
        md_mm_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        md_128 v_w = md_mm_loadu_ps(in_w + i);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
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
        double x = in_xyz[i].x;
        double y = in_xyz[i].y;
        double z = in_xyz[i].z;
        double w = in_w[i];
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += w * cos(tx);
        acc_s[0] += w * sin(tx);
        acc_c[1] += w * cos(ty);
        acc_s[1] += w * sin(ty);
        acc_c[2] += w * cos(tz);
        acc_s[2] += w * sin(tz);
        acc_w    += w;
    }

    const double inv_w = 1.0 / acc_w;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

// Internal versions for COM computation supporting triclinic and ortho PBC
static void _com_pbc(float out_com[3], const vec3_t* in_xyz, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_xyz);
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
        md_512 v_x, v_y, v_z;
        md_mm512_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
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
        md_256 v_x, v_y, v_z;
        md_mm256_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
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
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128 v_x, v_y, v_z;
        md_mm_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
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
        double x = in_xyz[i].x;
        double y = in_xyz[i].y;
        double z = in_xyz[i].z;
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += cos(tx);
        acc_s[0] += sin(tx);
        acc_c[1] += cos(ty);
        acc_s[1] += sin(ty);
        acc_c[2] += cos(tz);
        acc_s[2] += sin(tz);
    }

    const double inv_w = 1.0 / (double)count;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}


static void _com_pbc_i(float out_com[3], const vec3_t* in_xyz, const int32_t* in_idx, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_xyz);
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
        md_512 v_x, v_y, v_z;
        md_mm512_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
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
        md_256 v_x, v_y, v_z;
        md_mm256_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
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
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128i idx = md_mm_loadu_si128(in_idx + i);
        md_128 v_x, v_y, v_z;
        md_mm_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
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
        double x = in_xyz[idx].x;
        double y = in_xyz[idx].y;
        double z = in_xyz[idx].z;
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += cos(tx);
        acc_s[0] += sin(tx);
        acc_c[1] += cos(ty);
        acc_s[1] += sin(ty);
        acc_c[2] += cos(tz);
        acc_s[2] += sin(tz);
    }

    const double inv_w = 1.0 / (double)count;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

static void _com_pbc_iw(float out_com[3], const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const float M[3][3], const float I[3][3]) {
    ASSERT(out_com);
    ASSERT(in_xyz);
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
        md_512 v_x, v_y, v_z;
        md_mm512_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        md_512 v_w  = _mm512_i32gather_ps(idx, in_w, 4);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
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
        md_256 v_x, v_y, v_z;
        md_mm256_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        md_256 v_w  = md_mm256_i32gather_ps(in_w, idx, 4);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
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
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128i idx = md_mm_loadu_si128(in_idx + i);
        md_128 v_x, v_y, v_z;
        md_mm_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        md_128 v_w  = md_mm_i32gather_ps(in_w, idx, 4);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
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
        double x = in_xyz[idx].x;
        double y = in_xyz[idx].y;
        double z = in_xyz[idx].z;
        double w = in_w[idx];
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += w * cos(tx);
        acc_s[0] += w * sin(tx);
        acc_c[1] += w * cos(ty);
        acc_s[1] += w * sin(ty);
        acc_c[2] += w * cos(tz);
        acc_s[2] += w * sin(tz);
        acc_w    += w;
    }

    const double inv_w = 1.0 / acc_w;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

// This is uses the trigonometric com algorithm presented in: INSERT REF HERE
static void com_pbc(float out_com[3], const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {

    mat3_t A = { 0 };
    md_unitcell_A_extract_float(A.elem, cell);

    mat3_t Ai = { 0 };
    md_unitcell_I_extract_float(Ai.elem, cell);

    // Should transform cartesian coordinates to fractional coordinates and scale by 2 pi
    mat3_t M = mat3_mul(mat3_scale(TWO_PI, TWO_PI, TWO_PI), Ai);
    // Inverse transform
    mat3_t I = mat3_mul(A, mat3_scale(1.0/TWO_PI, 1.0/TWO_PI, 1.0/TWO_PI));

    if (!in_idx && md_unitcell_is_orthorhombic(cell)) {
        com_pbc_ortho_packed(out_com, in_xyz, in_w, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        return;
    }

    // Here we select the correct internal version based on the available inputs
    if (in_w) {
        if (in_idx) {
            _com_pbc_iw(out_com, in_xyz, in_w, in_idx, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        } else {
            _com_pbc_w(out_com, in_xyz, in_w, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        }
    } else {
        if (in_idx) {
            _com_pbc_i(out_com, in_xyz, in_idx, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        } else {
            _com_pbc(out_com, in_xyz, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        }
    }
}

static void com_vec4(float out_com[3], const vec4_t* in_xyzw, const int32_t* in_idx, size_t count) {
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
        vec3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
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
        mat3_t Ai = { 0 };
        md_unitcell_I_extract_float(Ai.elem, cell);

        mat3_t A = { 0 };
        md_unitcell_A_extract_float(A.elem, cell);

        // Cartesian -> fractional, scaled into angle space, and the inverse which takes the mean
        // angle back out to a Cartesian position. A triclinic basis mixes the components, so both
        // directions have to be full matrix products: resolving an axis on its own is only valid
        // when the basis is diagonal, which is why the orthorhombic branch above can get away with
        // per component scalars and this one cannot.
        const mat4x3_t M = mat4x3_from_mat3(mat3_mul(mat3_scale(TWO_PI, TWO_PI, TWO_PI), Ai));
        const mat3_t   I = mat3_mul(A, mat3_scale(1.0f / TWO_PI, 1.0f / TWO_PI, 1.0f / TWO_PI));

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
                vec4_t theta = mat4x3_mul_vec4(M, xyzw);
                vec4_t s,c;
                vec4_sincos(theta, &s, &c);
                acc_s = vec4_add(acc_s, vec4_mul(s, www1));
                acc_c = vec4_add(acc_c, vec4_mul(c, www1));
                acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
            }
        }

        vec3_t theta_vec = { 0 };
        for (int i = 0; i < 3; ++i) {
            const double y = acc_s.elem[i] / acc_xyzw.w;
            const double x = acc_c.elem[i] / acc_xyzw.w;
            const double r2 = x*x + y*y;
            double theta_prim = PI;
            if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
                theta_prim += atan2(-y, -x);
            }
            theta_vec.elem[i] = (float)theta_prim;
        }

        // The three mean angles ARE a fractional coordinate (times two pi). Carry it back to
        // Cartesian as one vector - component i of the result draws on all three angles.
        const vec3_t com = mat3_mul_vec3(I, theta_vec);
        out_com[0] = com.x;
        out_com[1] = com.y;
        out_com[2] = com.z;
    }
}

vec3_t md_util_com_compute(const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* unit_cell) {
    ASSERT(in_xyz);

    if (count == 0) {
        return (vec3_t) {0,0,0};
    }

    vec3_t xyz = {0};
    if (!unit_cell || unit_cell->flags == MD_UNITCELL_NONE) {
        com(xyz.elem, in_xyz, in_w, in_idx, count);
        return xyz;
    }

    if (unit_cell->flags & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC)) {
        com_pbc(xyz.elem, in_xyz, in_w, in_idx, count, unit_cell);
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

void md_util_distance_array(float* out_dist, const vec3_t* coord_a, size_t len_a, const vec3_t* coord_b, size_t len_b, const md_unitcell_t* cell) {
    if (cell->flags == 0) {
        for (size_t i = 0; i < len_a; ++i) {
            for (size_t j = 0; j < len_b; ++j) {
                out_dist[i * len_b + j] = vec3_distance(coord_a[i], coord_b[j]);
            }
        }
    }
    else if (cell->flags & MD_UNITCELL_ORTHO) {
		vec4_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
        for (size_t i = 0; i < len_a; ++i) {
            vec4_t a = vec4_from_vec3(coord_a[i], 0);
            for (size_t j = 0; j < len_b; ++j) {
                vec4_t b = vec4_from_vec3(coord_b[j], 0);
                out_dist[i * len_b + j] = vec4_periodic_distance(a, b, ext);
            }
        }
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        float A[3][3] = {0};
        md_unitcell_A_extract_float(A, cell);
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (size_t i = 0; i < len_a; ++i) {
            for (size_t j = 0; j < len_b; ++j) {
                vec3_t dx = vec3_sub(coord_a[i], coord_b[j]);
                minimum_image_triclinic(dx.elem, MD_AS_CONST_MAT3(A));
                out_dist[i * len_b + j] = vec3_length(dx);
            }
        }
    }
}

// ### MINIMUM AND MAXIMUM DISTANCE BETWEEN SETS ###
//
// For each group of points of a, the nearest (farthest) point of b under the minimum image convention.
//
// There is no spatial structure. The points of a group are taken eight at a time, a cluster in one SIMD register, and
// b is streamed once per cluster, each point of b tested against the sphere which bounds the cluster: a point of b can
// only improve on the best pair so far if it lies within best + R of the centre (beyond best - R for the farthest
// pair). All but a handful are rejected by that one test, and those are tested against every point of the cluster.
// The best pair is shared by the clusters of a group and starts from the pair of the previous group, which is usually
// close, so the test is tight from the first point. A grid over b would make the cost of a cluster sublinear in b,
// but building one costs more than the whole stream until a has a few hundred clusters, and the common case is a
// small a (a ligand, a residue, a QM region) against a large b.
//
// Minimum image. A displacement is first reduced into the brick |d_x| <= x/2, |d_y| <= y/2, |d_z| <= z/2: along c,
// then b, then a, which the lower triangular basis allows without disturbing the components already reduced. In an
// orthorhombic cell that is the minimum image. In a triclinic one it is the minimum image whenever that is shorter
// than R_SAFE = min(x, y, z) / 2: the brick is a fundamental domain of the lattice, the reduction returns the image
// inside it, and a vector shorter than R_SAFE lies inside it. So a reduced length below R_SAFE is exact, and one at or
// above it only says that the true length lies between R_SAFE and the reduced one. Where that does not decide the
// comparison at hand, the minimum image is searched for (dist_search_d2). That takes a group about as far from b as
// half the cell, or the farthest pair.
//
// The search is bounded rather than assumed: md_util_min_image_vec3 looks one lattice vector either way of the reduced
// image, which is the minimum image for a cell of a sensible shape but not for every cell. The bound is two at most for
// the usual cells (a reduced triclinic box, a rhombic dodecahedron, a truncated octahedron), and only a degenerate cell
// would reach DIST_MAX_SEARCH.

enum {
    DIST_PBC_NONE      = 0,
    DIST_PBC_ORTHO     = 1,
    DIST_PBC_TRICLINIC = 2,
};

// The search goes no further than this many lattice vectors along an axis. A cell which needs more is degenerate.
#define DIST_MAX_SEARCH 8

typedef struct dist_cell_t {
    // The basis, lower triangular: a = (ax, 0, 0), b = (bx, by, 0), c = (cx, cy, cz). Zero along an axis without a
    // period, as is the inverse of its diagonal (ix, iy, iz), which leaves such an axis alone in the reduction.
    md_256 ax, bx, by, cx, cy, cz;
    md_256 ix, iy, iz;
    float  basis[3][3];     // [vector][residue]
    float  dual[3];         // Lengths of the rows of the inverse basis, zero along an axis without a period
    float  rsafe2;          // A reduced length squared below this is the minimum image
    float  eps;             // Slack for the rounding of the coordinates in the sphere test
    int    mode;
} dist_cell_t;

static void dist_cell_init(dist_cell_t* c, const md_unitcell_t* cell, float coord_scale) {
    MEMSET(c, 0, sizeof(*c));
    const uint32_t flags = cell ? (uint32_t)cell->flags : 0;
    const bool pbc[3] = {
        (flags & MD_UNITCELL_PBC_X) != 0,
        (flags & MD_UNITCELL_PBC_Y) != 0,
        (flags & MD_UNITCELL_PBC_Z) != 0,
    };

    double A[3][3] = {{0}};
    if (cell) md_unitcell_A_extract_double(A, cell);
    for (int i = 0; i < 3; ++i) {
        if (!pbc[i] || A[i][i] == 0.0) {
            A[i][0] = A[i][1] = A[i][2] = 0.0;
        }
    }

    const bool any = A[0][0] != 0.0 || A[1][1] != 0.0 || A[2][2] != 0.0;
    c->mode = !any ? DIST_PBC_NONE : (flags & MD_UNITCELL_TRICLINIC) ? DIST_PBC_TRICLINIC : DIST_PBC_ORTHO;
    if (c->mode == DIST_PBC_ORTHO) {
        A[1][0] = A[2][0] = A[2][1] = 0.0;
    }

    double extent = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            c->basis[i][j] = (float)A[i][j];
            extent += fabs(A[i][j]);
        }
    }

    c->ax = md_mm256_set1_ps((float)A[0][0]);
    c->bx = md_mm256_set1_ps((float)A[1][0]);
    c->by = md_mm256_set1_ps((float)A[1][1]);
    c->cx = md_mm256_set1_ps((float)A[2][0]);
    c->cy = md_mm256_set1_ps((float)A[2][1]);
    c->cz = md_mm256_set1_ps((float)A[2][2]);
    c->ix = md_mm256_set1_ps(A[0][0] != 0.0 ? (float)(1.0 / A[0][0]) : 0.0f);
    c->iy = md_mm256_set1_ps(A[1][1] != 0.0 ? (float)(1.0 / A[1][1]) : 0.0f);
    c->iz = md_mm256_set1_ps(A[2][2] != 0.0 ? (float)(1.0 / A[2][2]) : 0.0f);

    // Several ulps of the largest magnitude a difference is taken between
    c->eps = (float)((coord_scale + extent) * 4.0e-6 + 1.0e-5);

    // The rows of the inverse of the basis, completed by a unit vector along each axis without a period. A lattice
    // vector A n has n_i = r_i . A n, so |n_i| <= |r_i| |A n|.
    double M[3][3];
    MEMCPY(M, A, sizeof(M));
    for (int i = 0; i < 3; ++i) {
        if (M[i][i] == 0.0) {
            M[i][0] = M[i][1] = M[i][2] = 0.0;
            M[i][i] = 1.0;
        }
    }
    const double x = M[0][0], xy = M[1][0], y = M[1][1], xz = M[2][0], yz = M[2][1], z = M[2][2];
    const double r0[3] = { 1.0 / x, -xy / (x * y), (xy * yz - xz * y) / (x * y * z) };
    const double r1[3] = { 0.0, 1.0 / y, -yz / (y * z) };
    const double r2[3] = { 0.0, 0.0, 1.0 / z };
    c->dual[0] = A[0][0] != 0.0 ? (float)sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]) : 0.0f;
    c->dual[1] = A[1][1] != 0.0 ? (float)sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]) : 0.0f;
    c->dual[2] = A[2][2] != 0.0 ? (float)sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]) : 0.0f;

    if (c->mode == DIST_PBC_TRICLINIC) {
        double rsafe = DBL_MAX;
        for (int i = 0; i < 3; ++i) {
            if (A[i][i] != 0.0) rsafe = MIN(rsafe, 0.5 * fabs(A[i][i]));
        }
        rsafe = MAX(0.0, rsafe - c->eps);
        c->rsafe2 = (float)(rsafe * rsafe);
    } else {
        c->rsafe2 = FLT_MAX;
    }
}

// Into the brick, along c, then b, then a
static FORCE_INLINE void dist_reduce(md_256* dx, md_256* dy, md_256* dz, const dist_cell_t* c, const int mode) {
    if (mode == DIST_PBC_ORTHO) {
        *dx = md_mm256_fnmadd_ps(md_mm256_round_ps(md_mm256_mul_ps(*dx, c->ix)), c->ax, *dx);
        *dy = md_mm256_fnmadd_ps(md_mm256_round_ps(md_mm256_mul_ps(*dy, c->iy)), c->by, *dy);
        *dz = md_mm256_fnmadd_ps(md_mm256_round_ps(md_mm256_mul_ps(*dz, c->iz)), c->cz, *dz);
    } else if (mode == DIST_PBC_TRICLINIC) {
        const md_256 nz = md_mm256_round_ps(md_mm256_mul_ps(*dz, c->iz));
        *dx = md_mm256_fnmadd_ps(nz, c->cx, *dx);
        *dy = md_mm256_fnmadd_ps(nz, c->cy, *dy);
        *dz = md_mm256_fnmadd_ps(nz, c->cz, *dz);
        const md_256 ny = md_mm256_round_ps(md_mm256_mul_ps(*dy, c->iy));
        *dx = md_mm256_fnmadd_ps(ny, c->bx, *dx);
        *dy = md_mm256_fnmadd_ps(ny, c->by, *dy);
        const md_256 nx = md_mm256_round_ps(md_mm256_mul_ps(*dx, c->ix));
        *dx = md_mm256_fnmadd_ps(nx, c->ax, *dx);
    }
}

static FORCE_INLINE md_256 dist_length2(md_256 dx, md_256 dy, md_256 dz) {
    return md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
}

// The minimum image length squared of reduced displacements, by search. A lattice vector A n which shortens a reduced
// d has |A n| <= |d| + |d + A n| <= 2 |d|, so |n_i| <= 2 |d| |r_i|. In a cell of a sensible shape that is one vector
// either way along each axis, or two along a short one.
static md_256 dist_search_d2(md_256 dx, md_256 dy, md_256 dz, const dist_cell_t* c) {
    md_256 best = dist_length2(dx, dy, dz);
    const float len = sqrtf(md_mm256_reduce_max_ps(best));
    int k[3];
    for (int i = 0; i < 3; ++i) {
        const float n = floorf(2.0f * len * c->dual[i] * 1.0001f);
        k[i] = (int)CLAMP(n, 0.0f, (float)DIST_MAX_SEARCH);
    }
    for (int nz = -k[2]; nz <= k[2]; ++nz) {
        for (int ny = -k[1]; ny <= k[1]; ++ny) {
            for (int nx = -k[0]; nx <= k[0]; ++nx) {
                if (nx == 0 && ny == 0 && nz == 0) continue;
                const float sx = nx * c->basis[0][0] + ny * c->basis[1][0] + nz * c->basis[2][0];
                const float sy = ny * c->basis[1][1] + nz * c->basis[2][1];
                const float sz = nz * c->basis[2][2];
                const md_256 vx = md_mm256_add_ps(dx, md_mm256_set1_ps(sx));
                const md_256 vy = md_mm256_add_ps(dy, md_mm256_set1_ps(sy));
                const md_256 vz = md_mm256_add_ps(dz, md_mm256_set1_ps(sz));
                best = md_mm256_min_ps(best, dist_length2(vx, vy, vz));
            }
        }
    }
    return best;
}

// Up to eight points of a group, as one image of them: each moved next to the first. Only the size of the bounding
// sphere depends on which image, the distances do not. Lanes past the end repeat the first point.
typedef struct dist_cluster_t {
    md_256  x, y, z;
    float   cen[3];
    float   rad;
    int64_t idx[8];
} dist_cluster_t;

static FORCE_INLINE void dist_cluster_init(dist_cluster_t* cl, const vec3_t* a, size_t beg, size_t n, const dist_cell_t* c, const int mode) {
    const vec3_t p0 = a[beg];
    float x[8], y[8], z[8];
    for (size_t i = 0; i < 8; ++i) {
        const size_t k = beg + (i < n ? i : 0);
        x[i] = a[k].x - p0.x;
        y[i] = a[k].y - p0.y;
        z[i] = a[k].z - p0.z;
        cl->idx[i] = (int64_t)k;
    }
    md_256 vx = md_mm256_loadu_ps(x);
    md_256 vy = md_mm256_loadu_ps(y);
    md_256 vz = md_mm256_loadu_ps(z);
    dist_reduce(&vx, &vy, &vz, c, mode);
    cl->x = md_mm256_add_ps(vx, md_mm256_set1_ps(p0.x));
    cl->y = md_mm256_add_ps(vy, md_mm256_set1_ps(p0.y));
    cl->z = md_mm256_add_ps(vz, md_mm256_set1_ps(p0.z));
    md_mm256_storeu_ps(x, cl->x);
    md_mm256_storeu_ps(y, cl->y);
    md_mm256_storeu_ps(z, cl->z);

    float lo[3] = { x[0], y[0], z[0] };
    float hi[3] = { x[0], y[0], z[0] };
    for (int i = 1; i < 8; ++i) {
        lo[0] = MIN(lo[0], x[i]); hi[0] = MAX(hi[0], x[i]);
        lo[1] = MIN(lo[1], y[i]); hi[1] = MAX(hi[1], y[i]);
        lo[2] = MIN(lo[2], z[i]); hi[2] = MAX(hi[2], z[i]);
    }
    cl->cen[0] = 0.5f * (lo[0] + hi[0]);
    cl->cen[1] = 0.5f * (lo[1] + hi[1]);
    cl->cen[2] = 0.5f * (lo[2] + hi[2]);
    float r2 = 0.0f;
    for (int i = 0; i < 8; ++i) {
        const float dx = x[i] - cl->cen[0], dy = y[i] - cl->cen[1], dz = z[i] - cl->cen[2];
        r2 = MAX(r2, dx * dx + dy * dy + dz * dz);
    }
    cl->rad = sqrtf(r2);
}

// The points of b, one array per axis, padded to a multiple of eight
typedef struct dist_set_t {
    const float* x;
    const float* y;
    const float* z;
    size_t n;
} dist_set_t;

// The squared distances from the points of the cluster to point j of b, exact wherever that can change the
// comparison with best2. If one improves on best2 it becomes the best pair and true is returned.
static FORCE_INLINE bool dist_cluster_test(float* best2, int64_t* ia, int64_t* ib, const dist_cluster_t* cl, const dist_set_t* b, size_t j, const dist_cell_t* c, const int mode, const bool largest) {
    md_256 dx = md_mm256_sub_ps(cl->x, md_mm256_set1_ps(b->x[j]));
    md_256 dy = md_mm256_sub_ps(cl->y, md_mm256_set1_ps(b->y[j]));
    md_256 dz = md_mm256_sub_ps(cl->z, md_mm256_set1_ps(b->z[j]));
    dist_reduce(&dx, &dy, &dz, c, mode);
    md_256 d2 = dist_length2(dx, dy, dz);
    if (mode == DIST_PBC_TRICLINIC) {
        // At or beyond R_SAFE a reduced length is an upper bound only: the true one lies in [R_SAFE, reduced]
        const md_256 far = md_mm256_cmpge_ps(d2, md_mm256_set1_ps(c->rsafe2));
        const int far_mask = md_mm256_movemask_ps(far);
        if (far_mask) {
            const bool matters = largest ? (md_mm256_movemask_ps(md_mm256_and_ps(far, md_mm256_cmpgt_ps(d2, md_mm256_set1_ps(*best2)))) != 0) : (c->rsafe2 < *best2);
            if (matters) {
                d2 = dist_search_d2(dx, dy, dz, c);
            }
        }
    }
    const float m = largest ? md_mm256_reduce_max_ps(d2) : md_mm256_reduce_min_ps(d2);
    if (largest ? !(m > *best2) : !(m < *best2)) return false;
    const int lane = ctz32((uint32_t)md_mm256_movemask_ps(md_mm256_cmpeq_ps(d2, md_mm256_set1_ps(m))));
    *best2 = m;
    *ia = cl->idx[lane];
    *ib = (int64_t)j;
    return true;
}

// The bound of the sphere test, squared: a point of b at a reduced distance from the centre of the cluster at or beyond
// it (for the farthest pair: at or within it) cannot improve on best2. Negative for the farthest pair while every
// point can.
static FORCE_INLINE float dist_threshold2(float best2, const dist_cluster_t* cl, const dist_cell_t* c, const bool largest) {
    if (largest) {
        if (best2 < 0.0f) return -1.0f;
        const float t = sqrtf(best2) - cl->rad - c->eps;
        return t > 0.0f ? t * t : -1.0f;
    }
    if (best2 >= FLT_MAX) return FLT_MAX;
    const float t = sqrtf(best2) + cl->rad + c->eps;
    return t * t;
}

static FORCE_INLINE void dist_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* a, const size_t* a_offsets, size_t num_groups,
                                     const dist_set_t* b, const dist_cell_t* c, const int mode, const bool largest) {
    int64_t prev = -1;
    for (size_t g = 0; g < num_groups; ++g) {
        const size_t beg = a_offsets[g];
        const size_t end = a_offsets[g + 1];
        float   best2 = largest ? -1.0f : FLT_MAX;
        int64_t ia = -1;
        int64_t ib = -1;
        bool done = false;

        for (size_t i = beg; i < end && !done; i += 8) {
            dist_cluster_t cl;
            dist_cluster_init(&cl, a, i, MIN(8, end - i), c, mode);

            // The partner of the previous group: a real pair, usually close to the best one
            if (prev >= 0) {
                dist_cluster_test(&best2, &ia, &ib, &cl, b, (size_t)prev, c, mode, largest);
            }

            const md_256 cx = md_mm256_set1_ps(cl.cen[0]);
            const md_256 cy = md_mm256_set1_ps(cl.cen[1]);
            const md_256 cz = md_mm256_set1_ps(cl.cen[2]);
            float  thr2 = dist_threshold2(best2, &cl, c, largest);
            md_256 vthr = md_mm256_set1_ps(thr2);
            const md_256 vrsafe2 = md_mm256_set1_ps(c->rsafe2);

            for (size_t j = 0; j < b->n && !done; j += 8) {
                md_256 dx = md_mm256_sub_ps(cx, md_mm256_loadu_ps(b->x + j));
                md_256 dy = md_mm256_sub_ps(cy, md_mm256_loadu_ps(b->y + j));
                md_256 dz = md_mm256_sub_ps(cz, md_mm256_loadu_ps(b->z + j));
                dist_reduce(&dx, &dy, &dz, c, mode);
                const md_256 d2 = dist_length2(dx, dy, dz);
                int mask = md_mm256_movemask_ps(largest ? md_mm256_cmpgt_ps(d2, vthr) : md_mm256_cmplt_ps(d2, vthr));
                if (mode == DIST_PBC_TRICLINIC && !largest && thr2 > c->rsafe2) {
                    // The sphere reaches past R_SAFE: a reduced distance there may overstate the true one
                    const int far = md_mm256_movemask_ps(md_mm256_cmpge_ps(d2, vrsafe2)) & ~mask;
                    if (far) {
                        mask |= md_mm256_movemask_ps(md_mm256_cmplt_ps(dist_search_d2(dx, dy, dz, c), vthr));
                    }
                }
                if (j + 8 > b->n) {
                    mask &= (1 << (b->n - j)) - 1;
                }
                while (mask) {
                    const size_t k = j + ctz32((uint32_t)mask);
                    mask &= mask - 1;
                    if (dist_cluster_test(&best2, &ia, &ib, &cl, b, k, c, mode, largest)) {
                        if (!largest && best2 <= 0.0f) {
                            // Nothing is nearer than a shared point
                            done = true;
                            break;
                        }
                        thr2 = dist_threshold2(best2, &cl, c, largest);
                        vthr = md_mm256_set1_ps(thr2);
                    }
                }
            }
        }

        out_dist[g] = ib >= 0 ? sqrtf(MAX(best2, 0.0f)) : 0.0f;
        if (out_idx_a) out_idx_a[g] = ia;
        if (out_idx_b) out_idx_b[g] = ib;
        if (ib >= 0) prev = ib;
    }
}

static void dist_groups_none_min (float* d, int64_t* ia, int64_t* ib, const vec3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_NONE,      false); }
static void dist_groups_ortho_min(float* d, int64_t* ia, int64_t* ib, const vec3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_ORTHO,     false); }
static void dist_groups_tri_min  (float* d, int64_t* ia, int64_t* ib, const vec3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_TRICLINIC, false); }
static void dist_groups_none_max (float* d, int64_t* ia, int64_t* ib, const vec3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_NONE,      true);  }
static void dist_groups_ortho_max(float* d, int64_t* ia, int64_t* ib, const vec3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_ORTHO,     true);  }
static void dist_groups_tri_max  (float* d, int64_t* ia, int64_t* ib, const vec3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_TRICLINIC, true);  }

static void distance_extent_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell, bool largest) {
    ASSERT(out_dist);
    ASSERT(a_offsets || num_groups == 0);
    if (num_groups == 0) return;

    if (num_b == 0) {
        for (size_t g = 0; g < num_groups; ++g) {
            out_dist[g] = 0.0f;
            if (out_idx_a) out_idx_a[g] = -1;
            if (out_idx_b) out_idx_b[g] = -1;
        }
        return;
    }
    ASSERT(coord_b);

    md_temp_scope_t temp = md_temp_begin();

    // b, one array per axis. The lanes past the end repeat the last point and are masked off.
    const size_t cap = ALIGN_TO(num_b, 8);
    float* bx = md_temp_alloc_array(temp, float, cap);
    float* by = md_temp_alloc_array(temp, float, cap);
    float* bz = md_temp_alloc_array(temp, float, cap);
    float scale = 0.0f;
    for (size_t i = 0; i < cap; ++i) {
        const vec3_t p = coord_b[MIN(i, num_b - 1)];
        bx[i] = p.x;
        by[i] = p.y;
        bz[i] = p.z;
        scale = MAX(scale, MAX(fabsf(p.x), MAX(fabsf(p.y), fabsf(p.z))));
    }
    for (size_t g = 0; g < num_groups; ++g) {
        for (size_t i = a_offsets[g]; i < a_offsets[g + 1]; ++i) {
            const vec3_t p = coord_a[i];
            scale = MAX(scale, MAX(fabsf(p.x), MAX(fabsf(p.y), fabsf(p.z))));
        }
    }

    dist_cell_t c;
    dist_cell_init(&c, cell, scale);
    const dist_set_t b = { bx, by, bz, num_b };

    switch (c.mode) {
    case DIST_PBC_NONE:      (largest ? dist_groups_none_max  : dist_groups_none_min) (out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, &b, &c); break;
    case DIST_PBC_ORTHO:     (largest ? dist_groups_ortho_max : dist_groups_ortho_min)(out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, &b, &c); break;
    case DIST_PBC_TRICLINIC: (largest ? dist_groups_tri_max   : dist_groups_tri_min)  (out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, &b, &c); break;
    default: ASSERT(false); break;
    }

    md_temp_end(temp);
}

void md_util_min_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    distance_extent_groups(out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, coord_b, num_b, cell, false);
}

void md_util_max_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    distance_extent_groups(out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, coord_b, num_b, cell, true);
}

float md_util_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    if (num_a == 0 || num_b == 0) return FLT_MAX;
    const size_t offsets[2] = { 0, num_a };
    float   dist = 0.0f;
    int64_t ia = -1;
    int64_t ib = -1;
    distance_extent_groups(&dist, &ia, &ib, coord_a, offsets, 1, coord_b, num_b, cell, false);
    if (out_idx_a) *out_idx_a = ia;
    if (out_idx_b) *out_idx_b = ib;
    return dist;
}

float md_util_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    if (num_a == 0 || num_b == 0) return 0.0f;
    const size_t offsets[2] = { 0, num_a };
    float   dist = 0.0f;
    int64_t ia = -1;
    int64_t ib = -1;
    distance_extent_groups(&dist, &ia, &ib, coord_a, offsets, 1, coord_b, num_b, cell, true);
    if (dist > 0.0f) {
        if (out_idx_a) *out_idx_a = ia;
        if (out_idx_b) *out_idx_b = ib;
    }
    return dist;
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

static inline void min_image_ortho(float dx[3], const float ext[3], const float inv_ext[3]) {
    dx[0] -= ext[0] * nearbyintf(dx[0] * inv_ext[0]);
    dx[1] -= ext[1] * nearbyintf(dx[1] * inv_ext[1]);
    dx[2] -= ext[2] * nearbyintf(dx[2] * inv_ext[2]);
}

// Reciprocal extent for min_image_ortho.
// @NOTE: this is NOT the half extent min_image_triclinic wants - the two take different quantities
// and the same local must not be handed to both.
// A zero extent means that axis is not periodic; a zero reciprocal makes the minimum image a no-op
// along it rather than producing a NaN.
static inline void min_image_inv_ext(float out_inv_ext[3], const float ext[3]) {
    out_inv_ext[0] = ext[0] != 0.0f ? 1.0f / ext[0] : 0.0f;
    out_inv_ext[1] = ext[1] != 0.0f ? 1.0f / ext[1] : 0.0f;
    out_inv_ext[2] = ext[2] != 0.0f ? 1.0f / ext[2] : 0.0f;
}

void md_util_min_image_vec3(vec3_t dx[], size_t count, const md_unitcell_t* cell) {
    if (cell) {
        vec3_t diag = { 0 };
        md_unitcell_diag_extract_float(diag.elem, cell);
        vec3_t half_diag = vec3_mul1(diag, 0.5f);
        uint32_t flags = md_unitcell_flags(cell);
        if (flags & MD_UNITCELL_ORTHO) {
            vec3_t inv_diag = { 0 };
            min_image_inv_ext(inv_diag.elem, diag.elem);
            for (size_t i = 0; i < count; ++i) {
                min_image_ortho(dx[i].elem, diag.elem, inv_diag.elem);
            }
        } else if (flags & MD_UNITCELL_TRICLINIC) {
            mat3_t A = { 0 };
            md_unitcell_A_extract_float(A.elem, cell);
            for (size_t i = 0; i < count; ++i) {
                min_image_triclinic(dx[i].elem, MD_AS_CONST_MAT3(A.elem), half_diag.elem);
            }
        }
    }
}

void md_util_min_image_vec4(vec4_t dx[], size_t count, const md_unitcell_t* cell) {
    if (cell) {
        vec3_t diag = { 0 };
        md_unitcell_diag_extract_float(diag.elem, cell);
        vec3_t half_diag = vec3_mul1(diag, 0.5f);
        uint32_t flags = md_unitcell_flags(cell);
        if (flags & MD_UNITCELL_ORTHO) {
            vec3_t inv_diag = { 0 };
            min_image_inv_ext(inv_diag.elem, diag.elem);
            for (size_t i = 0; i < count; ++i) {
                min_image_ortho(dx[i].elem, diag.elem, inv_diag.elem);
            }
        } else if (flags & MD_UNITCELL_TRICLINIC) {
            mat3_t A = { 0 };
            md_unitcell_A_extract_float(A.elem, cell);
            for (size_t i = 0; i < count; ++i) {
                min_image_triclinic(dx[i].elem, MD_AS_CONST_MAT3(A.elem), half_diag.elem);
            }
        }
    }
}

// Wraps into [0, ext) along each axis with a period; an axis without one (ext zero) is left alone.
// Without an index the axes are independent and it runs on the packed floats directly.
static void pbc_ortho(vec3_t* xyz, const int32_t* indices, size_t count, vec3_t box_ext) {
    ASSERT(xyz);

    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    const vec4_t ref = vec4_mul1(ext, 0.5f);

    size_t i = 0;
    if (!indices) {
        md_xv v_ext[XYZ_GROUP_VECS], v_inv[XYZ_GROUP_VECS], v_ref[XYZ_GROUP_VECS], v_refi[XYZ_GROUP_VECS], v_keep[XYZ_GROUP_VECS];
        xyz_lane_pattern(v_ext, ext.x, ext.y, ext.z);
        xyz_lane_pattern(v_inv, ext.x > 0 ? 1.0f / ext.x : 0.0f, ext.y > 0 ? 1.0f / ext.y : 0.0f, ext.z > 0 ? 1.0f / ext.z : 0.0f);
        xyz_lane_pattern(v_ref, ref.x, ref.y, ref.z);
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            v_keep[j] = md_xv_cmpeq_ps(v_ext[j], md_xv_setzero_ps());
            v_refi[j] = md_xv_mul_ps(v_ref[j], v_inv[j]);
        }

        // Scalar up to alignment, then eight atoms at a time
        const size_t head = xyz_align_peel(xyz, count);
        for (; i < head; ++i) {
            xyz[i] = vec3_from_vec4(vec4_deperiodize_ortho(vec4_from_vec3(xyz[i], 0), ref, ext));
        }
        float* f = (float*)xyz;
        for (; i + 8 <= count; i += 8) {
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                const md_xv v = md_xv_loadu_ps(f + i * 3 + j * MD_XV_WIDTH);
                // Into [0, ext]: d from a fused multiply-subtract, see the deperiodize comment in md_simd.h
                md_xv d = md_xv_fmsub_ps(v, v_inv[j], v_refi[j]);
                d = md_xv_sub_ps(d, md_xv_round_ps(d));
                const md_xv p = md_xv_fmadd_ps(d, v_ext[j], v_ref[j]);
                md_xv_storeu_ps(f + i * 3 + j * MD_XV_WIDTH, md_xv_select_ps(p, v, v_keep[j]));
            }
        }
    }
    for (; i < count; ++i) {
        const int32_t idx = indices ? indices[i] : (int32_t)i;
        const vec4_t pos = vec4_deperiodize_ortho(vec4_from_vec3(xyz[idx], 0), ref, ext);
        xyz[idx] = vec3_from_vec4(pos);
    }
}

static void pbc_ortho_vec4(vec4_t* xyzw, size_t count, vec3_t box_ext) {
    const vec4_t ext = vec4_from_vec3(box_ext, 0);
    const vec4_t ref = vec4_mul1(ext, 0.5f);
    for (size_t i = 0; i < count; ++i) {
        xyzw[i] = vec4_deperiodize_ortho(xyzw[i], ref, ext);
    }
}

// Into the cell: to fractional coordinates, the fraction of each, and back.
static void pbc_triclinic(vec3_t* xyz, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    ASSERT(xyz);

    mat3_t A = { 0 };
    md_unitcell_A_extract_float(A.elem, cell);

    mat3_t I = { 0 };
    md_unitcell_I_extract_float(I.elem, cell);

    size_t i = 0;
    if (!indices) {
        md_256 vA[3][3], vI[3][3];
        for (int c = 0; c < 3; ++c) {
            for (int r = 0; r < 3; ++r) {
                vA[c][r] = md_mm256_set1_ps(A.elem[c][r]);
                vI[c][r] = md_mm256_set1_ps(I.elem[c][r]);
            }
        }
        for (; i + 8 <= count; i += 8) {
            md_256 x, y, z;
            md_mm256_load_xyz_packed_ps(&x, &y, &z, (const float*)(xyz + i));
            md_256 f[3];
            for (int r = 0; r < 3; ++r) {
                f[r] = md_mm256_fmadd_ps(vI[0][r], x, md_mm256_fmadd_ps(vI[1][r], y, md_mm256_mul_ps(vI[2][r], z)));
                f[r] = md_mm256_fract_ps(f[r]);
            }
            md_256 c[3];
            for (int r = 0; r < 3; ++r) {
                c[r] = md_mm256_fmadd_ps(vA[0][r], f[0], md_mm256_fmadd_ps(vA[1][r], f[1], md_mm256_mul_ps(vA[2][r], f[2])));
            }
            md_mm256_store_xyz_packed_ps((float*)(xyz + i), c[0], c[1], c[2]);
        }
    }

    const mat4x3_t I4x3 = mat4x3_from_mat3(I);
    const mat4x3_t A4x3 = mat4x3_from_mat3(A);
    for (; i < count; ++i) {
        const int32_t idx = indices ? indices[i] : (int32_t)i;
        vec4_t c = mat4x3_mul_vec4(I4x3, vec4_from_vec3(xyz[idx], 0));
        c = vec4_fract(c);
        c = mat4x3_mul_vec4(A4x3, c);
        xyz[idx] = vec3_from_vec4(c);
    }
}

static void pbc_triclinic_vec4(vec4_t* xyzw, size_t count, const md_unitcell_t* cell) {
    vec4_t mask = { 0 };
    md_unitcell_pbc_mask_extract((int*)mask.elem, cell);

	mat3_t A = { 0 };
	md_unitcell_A_extract_float(A.elem, cell);

	mat3_t I = { 0 };
	md_unitcell_I_extract_float(I.elem, cell);

    mat4x3_t A4x3 = mat4x3_from_mat3(A);
    mat4x3_t I4x3 = mat4x3_from_mat3(I);

    for (size_t i = 0; i < count; ++i) {
        vec4_t c = xyzw[i];
        c = mat4x3_mul_vec4(I4x3, c);
        c = vec4_fract(c);
        c = mat4x3_mul_vec4(A4x3, c);
        xyzw[i] = vec4_blend(xyzw[i], c, mask);
    }
}

bool md_util_pbc(vec3_t* xyz, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    if (!xyz) {
        MD_LOG_ERROR("Missing required input: xyz");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing unitcell");
        return false;
    }

    uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        vec3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
        pbc_ortho(xyz, indices, count, ext);
        return true;
    } else if (flags & MD_UNITCELL_TRICLINIC) {
        pbc_triclinic(xyz, indices, count, cell);
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
        vec3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
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

bool md_util_system_pbc(md_system_state_t* state) {
    ASSERT(state);

    if (md_unitcell_is_orthorhombic(&state->unitcell)) {
        vec3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, &state->unitcell);
        pbc_ortho(state->xyz, 0, state->num_atoms, ext);
    } else if (md_unitcell_is_triclinic(&state->unitcell)) {
		pbc_triclinic(state->xyz, 0, state->num_atoms, &state->unitcell);
    }

    return true;
}

void md_util_unwrap_molecule(md_system_state_t* state, const md_molecule_t* mol) {
    ASSERT(state);
    ASSERT(mol);

    if (mol->count == 0) {
        return;
    }

    const bool is_ortho     = md_unitcell_is_orthorhombic(&state->unitcell);
    const bool is_triclinic = md_unitcell_is_triclinic(&state->unitcell);
    if (!is_ortho && !is_triclinic) {
        // Nothing to unwrap against
        return;
    }

    float ext[3] = { 0 };
    float inv_ext[3] = { 0 };
	float half_ext[3] = { 0 };
    float box[3][3] = { 0 };

    if (is_ortho) {
        md_unitcell_diag_extract_float(ext, &state->unitcell);
        min_image_inv_ext(inv_ext, ext);
    } else {
        md_unitcell_A_extract_float(box, &state->unitcell);
		half_ext[0] = 0.5f * vec3_length((vec3_t) { box[0][0], box[1][0], box[2][0] });
		half_ext[1] = 0.5f * vec3_length((vec3_t) { box[0][1], box[1][1], box[2][1] });
		half_ext[2] = 0.5f * vec3_length((vec3_t) { box[0][2], box[1][2], box[2][2] });
    }

    // md_molecule_data_t stores each molecule in BFS order, so a single forward pass places every
    // atom after the one it was reached from. No queue, no visited set and no bond lookups: the
    // traversal was already computed by md_util_system_infer_molecules.
    // The seed is its own parent, so it deperiodizes against itself and stays put - no branch.
    for (size_t i = 0; i < mol->count; ++i) {
        const int32_t a = mol->atom_idx[i];
        const int32_t p = mol->parent_idx[i];

        const vec3_t x = state->xyz[a];
        const vec3_t r = state->xyz[p];
	    vec3_t dx = {x.x - r.x, x.y - r.y, x.z - r.z};

        if (is_ortho) {
            min_image_ortho(dx.elem, ext, inv_ext);
        } else {
            min_image_triclinic(dx.elem, MD_AS_CONST_MAT3(box), half_ext);
        }

        // The atom is placed relative to its PARENT, not relative to where it happened to sit:
        // dx is already the minimum image of (x - r), so the unwrapped position is r + dx.
        // Adding dx to x instead yields 2x - r, which is wrong even when no image correction applies.
        state->xyz[a] = vec3_add(r, dx);
    }
}

void md_util_unwrap_system(md_system_state_t* state, const md_system_t* sys) {
    ASSERT(state);
    ASSERT(sys);

    for (size_t i = 0; i < sys->molecule.count; ++i) {
        const md_molecule_t mol = md_system_molecule(sys, i);
        md_util_unwrap_molecule(state, &mol);
    }
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
        vec4_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
        for (size_t i = 0; i < count; ++i) {
            const vec4_t pos = vec4_deperiodize_ortho(xyzw[i], ref_pos, ext);
            xyzw[i] = pos;
        }
        return true;
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        mat3_t A = { 0 };
        md_unitcell_A_extract_float(A.elem, cell);

        for (size_t i = 0; i < count; ++i) {
            vec4_t pos = xyzw[i];
            deperiodize_triclinic(pos.elem, ref_xyz.elem, MD_AS_CONST_MAT3(A.elem));
            xyzw[i] = pos;
        }
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    return false;
}

// PERIODIC IMAGE SELECTION
//
// Everything below shares one idea: the estimator picks the periodic image of each point, choosing
// whichever image minimises the estimator's own objective. That is strictly better than deciding the
// images up front and estimating afterwards, because a preprocessing pass has to commit using some
// criterion unrelated to what is about to be measured.
//
// The payoff is stability. Every image choice has a margin - how far the deciding quantity sits from
// the half cell flip boundary - and the margin is what determines whether a point flips images
// between two nearly identical frames. Deciding by distance to a centre gives an atom at radius r a
// margin of only (L/2 - r), and it is worst for the outer atoms which dominate a rotational fit.
// Deciding by fit residual gives a margin of nearly L/2, because minimising the residual is the
// whole point. The method that is most accurate here is also the one that is hardest to flip.

// Decoded cell, hoisted out of the inner loops.
typedef struct pbc_image_t {
    bool   periodic;
    bool   ortho;
    vec4_t ext;   // orthorhombic extent, w = 0 so the weight lane is carried through untouched
    mat3_t A;     // triclinic basis
    mat3_t I;     // triclinic inverse basis, Cartesian -> fractional
} pbc_image_t;

static pbc_image_t pbc_image_init(const md_unitcell_t* cell) {
    pbc_image_t im = {0};
    if (!cell) {
        return im;
    }
    if (md_unitcell_is_orthorhombic(cell)) {
        vec3_t ext = {0};
        md_unitcell_diag_extract_float(ext.elem, cell);
        im.periodic = true;
        im.ortho    = true;
        im.ext      = vec4_from_vec3(ext, 0);
    } else if (md_unitcell_is_triclinic(cell)) {
        md_unitcell_A_extract_float(im.A.elem, cell);
        md_unitcell_I_extract_float(im.I.elem, cell);
        im.periodic = true;
    }
    return im;
}

// The image of 'pos' nearest to 'ref'. The weight in w is carried through untouched.
static inline vec4_t pbc_image_nearest(const pbc_image_t* im, vec4_t pos, vec3_t ref) {
    if (!im->periodic) {
        return pos;
    }
    if (im->ortho) {
        return vec4_deperiodize_ortho(pos, vec4_from_vec3(ref, 0), im->ext);
    }
    vec4_t out = pos;
    deperiodize_triclinic(out.elem, ref.elem, im->A.elem);
    return out;
}

// The lattice vector nearest to 'delta'. An axis with no extent is not periodic and contributes
// nothing.
static inline vec3_t pbc_image_lattice_round(const pbc_image_t* im, vec3_t delta) {
    if (!im->periodic) {
        return vec3_zero();
    }
    if (im->ortho) {
        const vec3_t ext = vec3_from_vec4(im->ext);
        return vec3_set(ext.x != 0.0f ? roundf(delta.x / ext.x) * ext.x : 0.0f,
                        ext.y != 0.0f ? roundf(delta.y / ext.y) * ext.y : 0.0f,
                        ext.z != 0.0f ? roundf(delta.z / ext.z) * ext.z : 0.0f);
    }
    const vec3_t f = mat3_mul_vec3(im->I, delta);
    return mat3_mul_vec3(im->A, vec3_set(roundf(f.x), roundf(f.y), roundf(f.z)));
}

// KEEPING THE SET IN ITS OWN IMAGE
//
// The alternations below settle on images that are consistent RELATIVE to a centre, and that centre
// is seeded from the circular mean - whose angle is taken in [0, 2pi), so it always lands in the
// reference cell no matter where the input lives. Left at that, the whole set is quietly translated
// into image (0,0,0) as a side effect of being made whole.
//
// That is invisible to a caller that only reads relative geometry (an rmsd, a covariance, a
// rotation) and wrong for one that mixes the result back in with coordinates it did not hand over:
// a transform applied to the full system, an AABB queried against raw positions, a plane offset
// drawn over the atoms it was fitted to. Those callers get a result a whole cell vector off, and
// with a rotation in play, off by R times a cell vector, which is not even an image.
//
// So the contract is: make the set mutually consistent and otherwise leave it where it came in.
// The shift back is the lattice vector closest to the displacement of the weighted mean, i.e. the
// one that moves the input least. A set that genuinely straddles a boundary has no image to
// preserve - the two candidates are one lattice vector apart and equally good - and a caller that
// needs the set inside the cell wraps afterwards, explicitly.
static void pbc_image_restore(const pbc_image_t* im, vec4_t* xyzw, size_t count, vec3_t input_com, vec3_t* in_out_com) {
    const vec3_t off = pbc_image_lattice_round(im, vec3_sub(input_com, *in_out_com));
    if (vec3_dot(off, off) == 0.0f) {
        return;
    }
    const vec4_t off4 = vec4_from_vec3(off, 0);
    for (size_t i = 0; i < count; ++i) {
        xyzw[i] = vec4_add(xyzw[i], off4);
    }
    *in_out_com = vec3_add(*in_out_com, off);
}

// Squared distance below which two placements count as the same image. Generous: a genuine image
// change moves a point by a whole cell vector, so anything smaller is float noise.
#define PBC_IMAGE_EPS_SQ 1.0e-6f

// Alternation is a fixed point iteration and settles almost immediately - the assignment is usually
// stable on the second pass. The cap only bounds pathological input.
#define PBC_MAX_ITER 8

bool md_util_deperiodize_self_vec4(vec4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell, vec3_t* out_com) {
    if (!in_out_xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (count == 0) {
        if (out_com) *out_com = vec3_zero();
        return true;
    }

    const pbc_image_t im = pbc_image_init(cell);

    if (!im.periodic) {
        if (out_com) *out_com = md_util_com_compute_vec4(in_out_xyzw, NULL, count, NULL);
        return true;
    }

    // Where the input sits, before anything is moved. Only the image it implies is used - see
    // pbc_image_restore.
    const vec3_t input_com = md_util_com_compute_vec4(in_out_xyzw, NULL, count, NULL);

    // Seed with the circular mean. It needs no image assignment to compute, so unlike an arbitrary
    // reference point it cannot be biased by whichever images the coordinates happened to arrive in.
    vec3_t com = md_util_com_compute_vec4(in_out_xyzw, NULL, count, cell);

    for (int iter = 0; iter < PBC_MAX_ITER; ++iter) {
        bool stable = true;
        for (size_t i = 0; i < count; ++i) {
            const vec4_t before = in_out_xyzw[i];
            const vec4_t after  = pbc_image_nearest(&im, before, com);
            if (vec3_distance_squared(vec3_from_vec4(before), vec3_from_vec4(after)) > PBC_IMAGE_EPS_SQ) {
                stable = false;
            }
            in_out_xyzw[i] = after;
        }

        // Recompute the centre from the placed points. Both halves decrease the same objective.
        com = md_util_com_compute_vec4(in_out_xyzw, NULL, count, NULL);

        if (stable) {
            break;
        }
    }

    // Whole, and back in the image it arrived in.
    pbc_image_restore(&im, in_out_xyzw, count, input_com, &com);

    if (out_com) *out_com = com;
    return true;
}

float md_util_optimal_rotation_pbc_vec4_iter(mat3_t* out_rot, vec3_t* out_com, vec4_t* out_trg_xyzw,
                                             const vec4_t* ref_xyzw, vec3_t ref_com,
                                             const vec4_t* trg_xyzw, size_t count,
                                             const md_unitcell_t* cell, int max_iter, float tol) {
    ASSERT(ref_xyzw);
    ASSERT(trg_xyzw);

    if (out_rot) *out_rot = mat3_ident();
    if (out_com) *out_com = vec3_zero();

    if (count == 0) {
        return 0.0f;
    }

    const pbc_image_t im = pbc_image_init(cell);

    md_temp_scope_t temp = md_temp_begin();

    // The caller can have the placed points if it wants them (rmsd does); otherwise they are scratch.
    vec4_t* cur = out_trg_xyzw ? out_trg_xyzw : md_temp_alloc_array(temp, vec4_t, count);

    // Where the target sits before anything is moved. Read now: 'cur' is allowed to alias trg_xyzw
    // and the loop below writes through it.
    const vec3_t input_com = md_util_com_compute_vec4(trg_xyzw, NULL, count, NULL);

    // Image free seed, safe to compute before any assignment exists.
    vec3_t com = md_util_com_compute_vec4(trg_xyzw, NULL, count, im.periodic ? cell : NULL);
    mat3_t rot = mat3_ident();

    // Alternating minimisation of  sum_k w_k |R (q_k + A n_k - c) - p_k|^2  over R, c and n_k.
    // Each block has a closed form optimum given the other, so the objective decreases monotonically
    // and the discrete part reaches a genuine fixed point rather than a tolerance.
    bool stable = false;
    for (int iter = 0; iter < max_iter && !stable; ++iter) {
        // rot maps the target frame onto the reference frame, so its transpose carries the reference
        // template back out into the target's frame, which is where the prediction has to live.
        const mat3_t rot_inv = mat3_transpose(rot);

        stable = (iter > 0);
        for (size_t i = 0; i < count; ++i) {
            const vec3_t p    = vec3_sub(vec3_from_vec4(ref_xyzw[i]), ref_com);
            const vec3_t pred = vec3_add(com, mat3_mul_vec3(rot_inv, p));
            const vec4_t placed = pbc_image_nearest(&im, trg_xyzw[i], pred);

            if (iter > 0 && vec3_distance_squared(vec3_from_vec4(cur[i]), vec3_from_vec4(placed)) > tol) {
                stable = false;
            }
            cur[i] = placed;
        }

        // Refit. Both are the exact optimum for this assignment.
        com = md_util_com_compute_vec4(cur, NULL, count, NULL);

        const vec4_t* const sets[2] = { ref_xyzw, cur };
        const vec3_t        coms[2] = { ref_com,  com };
        rot = mat3_extract_rotation(mat3_cross_covariance_matrix_vec4(sets, NULL, count, coms));
    }

    // Back in the image the target arrived in. The rotation and the residual below are measured
    // relative to 'com', so both are untouched by this.
    pbc_image_restore(&im, cur, count, input_com, &com);

    // The largest per point residual is the distance to an ambiguous image choice. Small compared to
    // half a cell means the assignment above is not close to flipping.
    float max_dist_sq = 0.0f;
    for (size_t i = 0; i < count; ++i) {
        const vec3_t p = vec3_sub(vec3_from_vec4(ref_xyzw[i]), ref_com);
        const vec3_t q = mat3_mul_vec3(rot, vec3_sub(vec3_from_vec4(cur[i]), com));
        const float  d = vec3_distance_squared(p, q);
        max_dist_sq = MAX(max_dist_sq, d);
    }

    md_temp_end(temp);

    if (out_rot) *out_rot = rot;
    if (out_com) *out_com = com;

    return sqrtf(max_dist_sq);
}

void md_util_convert_to_relative_coordinates_vec4(vec4_t in_out_xyzw[], vec3_t ref_xyz, size_t count, const md_unitcell_t* cell) {
	ASSERT(in_out_xyzw);
	vec4_t com = vec4_from_vec3(ref_xyz, 0);
	for (size_t i = 0; i < count; ++i) {
		in_out_xyzw[i] = vec4_sub(in_out_xyzw[i], com);
	}
	if (cell) {
		md_util_min_image_vec4(in_out_xyzw, count, cell);
	}
}

mat3_t md_util_optimal_rotation_rel_vec4(const vec4_t* ref_rel_xyzw, const vec4_t* trg_rel_xyzw, size_t count) {
    ASSERT(ref_rel_xyzw);
    ASSERT(trg_rel_xyzw);

    if (count == 0) {
        return mat3_ident();
    }

    double A[3][3] = { 0 };
    double w_sum = 0.0;

    for (size_t i = 0; i < count; ++i) {
		const vec4_t p = ref_rel_xyzw[i];
		const vec4_t q = trg_rel_xyzw[i];

        // The question here is how to combine the weights.
        // For now we just take the average.
        // This should be equivalent to the other case where the weights for both sets are equal.
        const float w = (p.w + q.w) * 0.5f;

        A[0][0] += w * p.x * q.x;
        A[0][1] += w * p.x * q.y;
        A[0][2] += w * p.x * q.z;
        A[1][0] += w * p.y * q.x;
        A[1][1] += w * p.y * q.y;
        A[1][2] += w * p.y * q.z;
        A[2][0] += w * p.z * q.x;
        A[2][1] += w * p.z * q.y;
        A[2][2] += w * p.z * q.z;
        w_sum += w;
    }

	if (w_sum == 0.0) {
		return mat3_ident();
	}

    mat3_t C;
	C.elem[0][0] = (float)(A[0][0] / w_sum);
    C.elem[0][1] = (float)(A[0][1] / w_sum);
	C.elem[0][2] = (float)(A[0][2] / w_sum);
	C.elem[1][0] = (float)(A[1][0] / w_sum);
	C.elem[1][1] = (float)(A[1][1] / w_sum);
	C.elem[1][2] = (float)(A[1][2] / w_sum);
	C.elem[2][0] = (float)(A[2][0] / w_sum);
	C.elem[2][1] = (float)(A[2][1] / w_sum);
	C.elem[2][2] = (float)(A[2][2] / w_sum);

    return mat3_extract_rotation(C);
}

double md_util_rmsd_compute(const vec3_t* const in_xyz[2], const float* const in_w[2], const int32_t* const in_idx[2], const size_t count, const vec3_t in_com[2]) {
    const mat3_t R = mat3_optimal_rotation(in_xyz, in_w, in_idx, count, in_com);
    double d_sum = 0;
    double w_sum = 0;
    if (in_idx) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx0 = in_idx[0][i];
            int32_t idx1 = in_idx[1][i];
            vec3_t u = {in_xyz[0][idx0].x - in_com[0].x, in_xyz[0][idx0].y - in_com[0].y, in_xyz[0][idx0].z - in_com[0].z};
            vec3_t v = {in_xyz[1][idx1].x - in_com[1].x, in_xyz[1][idx1].y - in_com[1].y, in_xyz[1][idx1].z - in_com[1].z};
            vec3_t vp = mat3_mul_vec3(R, v);
            vec3_t d = vec3_sub(u, vp);
            float weight = in_w ? (in_w[0][idx0] + in_w[1][idx1]) * 0.5f : 1.0f;
            d_sum += weight * vec3_dot(d, d);
            w_sum += weight;
        }
    } else {
        // A native vector of atoms at a time: the deviation of each from the rotated second set
        size_t i = 0;
        md_xv vR[3][3];
        for (int c = 0; c < 3; ++c) for (int r = 0; r < 3; ++r) vR[c][r] = md_xv_set1_ps(R.elem[c][r]);
        const md_xv c0[3] = { md_xv_set1_ps(in_com[0].x), md_xv_set1_ps(in_com[0].y), md_xv_set1_ps(in_com[0].z) };
        const md_xv c1[3] = { md_xv_set1_ps(in_com[1].x), md_xv_set1_ps(in_com[1].y), md_xv_set1_ps(in_com[1].z) };
        md_xv vd = md_xv_setzero_ps();
        md_xv vw = md_xv_setzero_ps();
        for (; i + MD_XV_WIDTH <= count; i += MD_XV_WIDTH) {
            md_xv u[3], v[3];
            md_xv_load_xyz_packed_ps(&u[0], &u[1], &u[2], (const float*)(in_xyz[0] + i));
            md_xv_load_xyz_packed_ps(&v[0], &v[1], &v[2], (const float*)(in_xyz[1] + i));
            for (int k = 0; k < 3; ++k) {
                u[k] = md_xv_sub_ps(u[k], c0[k]);
                v[k] = md_xv_sub_ps(v[k], c1[k]);
            }
            md_xv d2 = md_xv_setzero_ps();
            for (int r = 0; r < 3; ++r) {
                // (R v)_r = R[0][r] v.x + R[1][r] v.y + R[2][r] v.z, R column major
                const md_xv vp = md_xv_fmadd_ps(vR[0][r], v[0], md_xv_fmadd_ps(vR[1][r], v[1], md_xv_mul_ps(vR[2][r], v[2])));
                const md_xv d  = md_xv_sub_ps(u[r], vp);
                d2 = md_xv_fmadd_ps(d, d, d2);
            }
            md_xv w = md_xv_set1_ps(1.0f);
            if (in_w) {
                w = md_xv_mul_ps(md_xv_add_ps(md_xv_loadu_ps(in_w[0] + i), md_xv_loadu_ps(in_w[1] + i)), md_xv_set1_ps(0.5f));
            }
            vd = md_xv_fmadd_ps(w, d2, vd);
            vw = md_xv_add_ps(vw, w);
        }
        d_sum += md_xv_reduce_add_ps(vd);
        w_sum += md_xv_reduce_add_ps(vw);

        for (; i < count; ++i) {
            vec3_t u = {in_xyz[0][i].x - in_com[0].x, in_xyz[0][i].y - in_com[0].y, in_xyz[0][i].z - in_com[0].z};
            vec3_t v = {in_xyz[1][i].x - in_com[1].x, in_xyz[1][i].y - in_com[1].y, in_xyz[1][i].z - in_com[1].z};
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

// One group of eight atoms, packed floats as they lie. The second frame is brought to the minimum
// image of the first when the cell is periodic.
static inline void lerp_xyz8(float* dst, const float* a, const float* b, const md_xv ext[XYZ_GROUP_VECS], const md_xv inv[XYZ_GROUP_VECS], bool periodic, float t) {
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
        const md_xv x0 = md_xv_loadu_ps(a + j * MD_XV_WIDTH);
        md_xv x1 = md_xv_loadu_ps(b + j * MD_XV_WIDTH);
        if (periodic) x1 = xv_deperiodize_ortho(x1, x0, ext[j], inv[j]);
        md_xv_storeu_ps(dst + j * MD_XV_WIDTH, md_xv_fmadd_ps(md_xv_sub_ps(x1, x0), md_xv_set1_ps(t), x0));
    }
}

// One group of eight atoms over four frames, weighted (see md_util_interpolate_cubic_weighted). The weights are applied
// in order, so a weight of exactly one with the others zero gives that frame back exactly: an interpolating scheme
// lands on its frames.
static inline void cubic_xyz8(float* dst, const float* const src[4], const md_xv ext[XYZ_GROUP_VECS], const md_xv inv[XYZ_GROUP_VECS], bool periodic, const md_xv w[4]) {
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
        md_xv x0 = md_xv_loadu_ps(src[0] + j * MD_XV_WIDTH);
        const md_xv x1 = md_xv_loadu_ps(src[1] + j * MD_XV_WIDTH);
        md_xv x2 = md_xv_loadu_ps(src[2] + j * MD_XV_WIDTH);
        md_xv x3 = md_xv_loadu_ps(src[3] + j * MD_XV_WIDTH);
        if (periodic) {
            x0 = xv_deperiodize_ortho(x0, x1, ext[j], inv[j]);
            x2 = xv_deperiodize_ortho(x2, x1, ext[j], inv[j]);
            x3 = xv_deperiodize_ortho(x3, x2, ext[j], inv[j]);
        }
        md_xv r = md_xv_mul_ps(w[0], x0);
        r = md_xv_fmadd_ps(w[1], x1, r);
        r = md_xv_fmadd_ps(w[2], x2, r);
        r = md_xv_fmadd_ps(w[3], x3, r);
        md_xv_storeu_ps(dst + j * MD_XV_WIDTH, r);
    }
}

// Linear interpolation of packed coordinates. Without a cell or in an orthorhombic one the axes are
// independent and it runs on the packed floats directly, eight atoms at a time, with a few leading
// atoms first so the loads are aligned and a partial group last, each through a small staging copy;
// any count works. A triclinic cell couples the axes: that path splits eight atoms into x, y and z
// and relies on the arrays being padded to a multiple of eight atoms.
bool md_util_interpolate_linear(vec3_t* out_xyz, const vec3_t* const in_xyz[2], size_t count, const md_unitcell_t* cell, float t) {
    ASSERT(out_xyz);
    ASSERT(in_xyz);
    ASSERT(in_xyz[0]);
    ASSERT(in_xyz[1]);

    t = CLAMP(t, 0.0f, 1.0f);

    const bool ortho = cell && md_unitcell_is_orthorhombic(cell);
    const bool tricl = cell && md_unitcell_is_triclinic(cell);

    mat3_t A = { 0 };
    mat3_t I = { 0 };
    if (ortho | tricl) {
        md_unitcell_A_extract_float(A.elem, cell);
        md_unitcell_I_extract_float(I.elem, cell);
    }

    const float* src0 = (const float*)in_xyz[0];
    const float* src1 = (const float*)in_xyz[1];
    float* dst = (float*)out_xyz;

    if (!tricl) {
        md_xv ext[XYZ_GROUP_VECS], inv[XYZ_GROUP_VECS];
        xyz_lane_pattern(ext, A.elem[0][0], A.elem[1][1], A.elem[2][2]);
        xyz_lane_pattern(inv, I.elem[0][0], I.elem[1][1], I.elem[2][2]);

        size_t i = 0;
        const size_t head = xyz_align_peel(src0, count);
        while (i < count) {
            const size_t n = (i == 0 && head) ? head : MIN(8, count - i);
            if (n == 8) {
                lerp_xyz8(dst + i * 3, src0 + i * 3, src1 + i * 3, ext, inv, ortho, t);
            } else {
                float a[24] = {0}, b[24] = {0}, r[24];
                MEMCPY(a, src0 + i * 3, n * sizeof(vec3_t));
                MEMCPY(b, src1 + i * 3, n * sizeof(vec3_t));
                lerp_xyz8(r, a, b, ext, inv, ortho, t);
                MEMCPY(dst + i * 3, r, n * sizeof(vec3_t));
            }
            i += n;
        }
        return true;
    }

    for (size_t i = 0; i < count; i += 8) {
        md_256 x0[3], x1[3];
        md_mm256_load_xyz_packed_ps(&x0[0], &x0[1], &x0[2], src0 + i * 3);
        md_mm256_load_xyz_packed_ps(&x1[0], &x1[1], &x1[2], src1 + i * 3);

        simd_deperiodize_triclinic(x1, x0, MD_AS_CONST_MAT3(A.elem));

        const md_256 x = md_mm256_lerp_ps(x0[0], x1[0], t);
        const md_256 y = md_mm256_lerp_ps(x0[1], x1[1], t);
        const md_256 z = md_mm256_lerp_ps(x0[2], x1[2], t);

        md_mm256_store_xyz_packed_ps(dst + i * 3, x, y, z);
    }

    return true;
}

void md_util_cubic_spline_weights(float w[4], float t, float s) {
    ASSERT(w);
    // The Hermite form, p(t) = h00 p1 + h10 m1 + h01 p2 + h11 m2, with the tangents m1 = s (p2 - p0) and m2 = s (p3 - p1)
    const float t2 = t * t;
    const float t3 = t2 * t;
    const float h00 =  2.0f * t3 - 3.0f * t2 + 1.0f;
    const float h10 =         t3 - 2.0f * t2 + t;
    const float h01 = -2.0f * t3 + 3.0f * t2;
    const float h11 =         t3 -        t2;
    w[0] = -s * h10;
    w[1] = h00 - s * h11;
    w[2] = h01 + s * h10;
    w[3] = s * h11;
}

void md_util_cubic_bspline_weights(float w[4], float t) {
    ASSERT(w);
    const float t2 = t * t;
    const float t3 = t2 * t;
    const float u  = 1.0f - t;
    w[0] = u * u * u * (1.0f / 6.0f);
    w[1] = ( 3.0f * t3 - 6.0f * t2 + 4.0f) * (1.0f / 6.0f);
    w[2] = (-3.0f * t3 + 3.0f * t2 + 3.0f * t + 1.0f) * (1.0f / 6.0f);
    w[3] = t3 * (1.0f / 6.0f);
}

bool md_util_interpolate_cubic_spline(vec3_t* out_xyz, const vec3_t* const in_xyz[4], size_t count, const md_unitcell_t* cell, float t, float s) {
    float w[4];
    md_util_cubic_spline_weights(w, CLAMP(t, 0.0f, 1.0f), CLAMP(s, 0.0f, 1.0f));
    return md_util_interpolate_cubic_weighted(out_xyz, in_xyz, count, cell, w);
}

// As md_util_interpolate_linear: packed floats directly without a cell or in an orthorhombic one,
// split into x, y and z (and padded) in a triclinic one.
bool md_util_interpolate_cubic_weighted(vec3_t* out_xyz, const vec3_t* const in_xyz[4], size_t count, const md_unitcell_t* cell, const float w[4]) {
    ASSERT(out_xyz);
    ASSERT(in_xyz);
    ASSERT(in_xyz[0]);
    ASSERT(in_xyz[1]);
    ASSERT(in_xyz[2]);
    ASSERT(in_xyz[3]);
    ASSERT(w);

    const bool ortho = cell && md_unitcell_is_orthorhombic(cell);
    const bool tricl = cell && md_unitcell_is_triclinic(cell);
    mat3_t A = { 0 };
    mat3_t I = { 0 };
    if (ortho | tricl) {
        md_unitcell_A_extract_float(A.elem, cell);
        md_unitcell_I_extract_float(I.elem, cell);
    }

    const md_256 vw[4] = { md_mm256_set1_ps(w[0]), md_mm256_set1_ps(w[1]), md_mm256_set1_ps(w[2]), md_mm256_set1_ps(w[3]) };
    const md_xv  xw[4] = { md_xv_set1_ps(w[0]),    md_xv_set1_ps(w[1]),    md_xv_set1_ps(w[2]),    md_xv_set1_ps(w[3]) };

    const float* src[4] = { (const float*)in_xyz[0], (const float*)in_xyz[1], (const float*)in_xyz[2], (const float*)in_xyz[3] };
    float* dst = (float*)out_xyz;

    if (!tricl) {
        md_xv ext[XYZ_GROUP_VECS], inv[XYZ_GROUP_VECS];
        xyz_lane_pattern(ext, A.elem[0][0], A.elem[1][1], A.elem[2][2]);
        xyz_lane_pattern(inv, I.elem[0][0], I.elem[1][1], I.elem[2][2]);

        size_t i = 0;
        const size_t head = xyz_align_peel(src[1], count);
        while (i < count) {
            const size_t n = (i == 0 && head) ? head : MIN(8, count - i);
            if (n == 8) {
                const float* const g[4] = { src[0] + i * 3, src[1] + i * 3, src[2] + i * 3, src[3] + i * 3 };
                cubic_xyz8(dst + i * 3, g, ext, inv, ortho, xw);
            } else {
                float stage[4][24] = {0}, r[24];
                for (int k = 0; k < 4; ++k) {
                    MEMCPY(stage[k], src[k] + i * 3, n * sizeof(vec3_t));
                }
                const float* const g[4] = { stage[0], stage[1], stage[2], stage[3] };
                cubic_xyz8(r, g, ext, inv, ortho, xw);
                MEMCPY(dst + i * 3, r, n * sizeof(vec3_t));
            }
            i += n;
        }
        return true;
    }

    for (size_t i = 0; i < count; i += 8) {
        md_256 x0[3], x1[3], x2[3], x3[3];
        md_mm256_load_xyz_packed_ps(&x0[0], &x0[1], &x0[2], src[0] + i * 3);
        md_mm256_load_xyz_packed_ps(&x1[0], &x1[1], &x1[2], src[1] + i * 3);
        md_mm256_load_xyz_packed_ps(&x2[0], &x2[1], &x2[2], src[2] + i * 3);
        md_mm256_load_xyz_packed_ps(&x3[0], &x3[1], &x3[2], src[3] + i * 3);

        simd_deperiodize_triclinic(x0, x1, MD_AS_CONST_MAT3(A.elem));
        simd_deperiodize_triclinic(x2, x1, MD_AS_CONST_MAT3(A.elem));
        simd_deperiodize_triclinic(x3, x2, MD_AS_CONST_MAT3(A.elem));

        md_256 r[3];
        for (int k = 0; k < 3; ++k) {
            r[k] = md_mm256_mul_ps(vw[0], x0[k]);
            r[k] = md_mm256_fmadd_ps(vw[1], x1[k], r[k]);
            r[k] = md_mm256_fmadd_ps(vw[2], x2[k], r[k]);
            r[k] = md_mm256_fmadd_ps(vw[3], x3[k], r[k]);
        }

        md_mm256_store_xyz_packed_ps(dst + i * 3, r[0], r[1], r[2]);
    }

    return true;
}

// The segments of one range, of consecutive residues from res_base, or of none (res_base -1, BACKBONES FROM THE
// BONDS)
static inline void commit_protein_backbone(md_amino_acid_atoms_t* bb_atoms, size_t bb_length, md_residue_idx_t res_base, md_part_idx_t part_idx, md_protein_backbone_data_t* bb_data, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(bb_data->segment.atoms);
    md_array_push_array(bb_data->segment.atoms, bb_atoms, bb_length, alloc);
    md_array_push(bb_data->range.offset, offset, alloc);
    md_array_push(bb_data->range.part_idx, part_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        const md_residue_idx_t res_idx = res_base < 0 ? -1 : res_base + (md_residue_idx_t)i;
        md_array_push(bb_data->segment.res_idx, res_idx, alloc);
    }
}

// The bead which stands for the backbone of an amino acid of a coarse grained model (the BB bead of Martini): its one
// bead with the backbone role. -1 for a residue which is no amino acid, or has none or several.
static md_atom_idx_t amino_acid_backbone_bead(const md_system_t* sys, size_t res_idx) {
    if (md_system_residue_kind(sys, res_idx) != MD_RESIDUE_KIND_AMINO_ACID) return -1;
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    md_atom_idx_t bead = -1;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if ((md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) && md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_BEAD) {
            if (bead >= 0) return -1;
            bead = (md_atom_idx_t)i;
        }
    }
    return bead;
}

// The backbone beads of a nucleotide of a coarse grained model, in the fields of the atoms they hold (Martini 2 DNA and
// RNA, Uusitalo et al. 2015, 2017; Martini 3 RNA of reForge): the phosphate bead BB1 in p, the sugar bead BB2 which
// holds the C4' in c4 and the sugar bead BB3 which holds the C1' in c1, the other fields -1. A 5' end has no BB1. False
// for a residue which is no nucleotide, or lacks BB2 or BB3. Their o5 is -1, which that of atoms never is
// (md_util_nucleic_acid_atoms_extract requires it): that tells beads from atoms (nucleic_atoms_are_beads).
static bool nucleotide_backbone_beads(md_nucleic_acid_atoms_t* out, const md_system_t* sys, size_t res_idx) {
    if (md_system_residue_kind(sys, res_idx) != MD_RESIDUE_KIND_NUCLEOTIDE) return false;
    md_nucleic_acid_atoms_t na;
    MEMSET(&na, -1, sizeof(na));
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) || md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_BEAD) continue;
        const str_t name = md_atom_name(&sys->atom, i);
        if      (str_eq_cstr(name, "BB1")) na.p  = (md_atom_idx_t)i;
        else if (str_eq_cstr(name, "BB2")) na.c4 = (md_atom_idx_t)i;
        else if (str_eq_cstr(name, "BB3")) na.c1 = (md_atom_idx_t)i;
    }
    if (out) *out = na;
    return na.c4 >= 0 && na.c1 >= 0;
}

static inline bool nucleic_atoms_are_beads(const md_nucleic_acid_atoms_t* na) {
    return na->o5 < 0;
}

// ### CHAIN BREAKS ### (see md_util.h)

// Ångström within which the atoms of a link are linked when they have no bonds: a peptide bond (C-N) is 1.33 long, a
// phosphodiester bond (O3'-P) 1.61, and the shortest distance over a missing residue well above 3
#define CHAIN_LINK_MAX_DIST 2.0f

// Whether the atoms a and b of consecutive residues of a chain link their backbones. The bonds decide when the system has
// bonds for both. Atoms without bonds (a structure file without bond inference) are linked within max_dist of each
// other, the shortest way round a periodic cell, when max_dist > 0 and the state has coordinates; otherwise when their
// sequence numbers are consecutive (seq_next). An atom which is not there (-1) links nothing.
static bool chain_atoms_linked(const md_system_t* sys, const md_system_state_t* state, md_atom_idx_t a, md_atom_idx_t b, float max_dist, bool seq_next) {
    if (a < 0 || b < 0) return false;
    const md_bond_iter_t it_a = md_bond_iter(&sys->bond, (size_t)a);
    const md_bond_iter_t it_b = md_bond_iter(&sys->bond, (size_t)b);
    if (md_bond_iter_has_next(&it_a) && md_bond_iter_has_next(&it_b)) {
        return md_bond_find(&sys->bond, a, b) >= 0;
    }
    if (max_dist > 0.0f && md_system_state_has_coords(state) && (size_t)a < state->num_atoms && (size_t)b < state->num_atoms) {
        vec3_t d = vec3_sub(md_state_coord(state, (size_t)b), md_state_coord(state, (size_t)a));
        md_util_min_image_vec3(&d, 1, &state->unitcell);
        return vec3_dot(d, d) <= max_dist * max_dist;
    }
    return seq_next;
}

// Whether the consecutive residues of a protein chain with the backbone atoms prev and next, both atomistic or both
// traces, are linked: by a peptide bond from the C of prev to the N of next, or for the beads of a trace (coarse grained)
// by a bond between them. Beads without bonds are linked by consecutive sequence numbers alone, not by their distance:
// beads are bonded farther apart than atoms are.
static bool protein_backbone_linked(const md_system_t* sys, const md_system_state_t* state, const md_amino_acid_atoms_t* prev, const md_amino_acid_atoms_t* next, bool seq_next) {
    if (prev->n < 0) {
        return chain_atoms_linked(sys, state, prev->ca, next->ca, 0.0f, seq_next);
    }
    return chain_atoms_linked(sys, state, prev->c, next->n, CHAIN_LINK_MAX_DIST, seq_next);
}

// Whether the consecutive residues of a nucleic acid chain with the backbone atoms prev and next, both atomistic or
// both beads, are linked: by a phosphodiester bond from the O3' of prev to the P of next, or for beads by a bond from a
// sugar bead of prev to the phosphate bead of next (BB3 to BB1 in Martini 2, BB2 and BB3 to BB1 in Martini 3), by
// consecutive sequence numbers alone without bonds. Without its P or phosphate bead (a 5' end) next is linked to
// nothing.
static bool nucleic_backbone_linked(const md_system_t* sys, const md_system_state_t* state, const md_nucleic_acid_atoms_t* prev, const md_nucleic_acid_atoms_t* next, bool seq_next) {
    if (nucleic_atoms_are_beads(prev)) {
        return chain_atoms_linked(sys, state, prev->c1, next->p, 0.0f, seq_next) || chain_atoms_linked(sys, state, prev->c4, next->p, 0.0f, seq_next);
    }
    return chain_atoms_linked(sys, state, prev->o3, next->p, CHAIN_LINK_MAX_DIST, seq_next);
}

// The atoms of element z bonded to the atom, and whether it has bonds at all
static size_t bonded_count_of_element(const md_system_t* sys, md_atom_idx_t atom, md_atomic_number_t z, bool* has_bonds) {
    size_t count = 0;
    *has_bonds = false;
    if (atom < 0) return 0;
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)atom);
    while (md_bond_iter_has_next(&it)) {
        *has_bonds = true;
        count += md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(&it)) == z;
        md_bond_iter_next(&it);
    }
    return count;
}

// Whether the carboxyl group of an amino acid (its C) is free: its C bonded to two oxygens, O and OXT. Without bonds, by
// the name of the second oxygen.
static bool amino_acid_carboxyl_free(const md_system_t* sys, size_t res_idx, md_atom_idx_t c) {
    bool has_bonds;
    const size_t num_o = bonded_count_of_element(sys, c, MD_Z_O, &has_bonds);
    if (has_bonds) return num_o >= 2;
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    for (uint32_t i = range.beg; i < range.end; ++i) {
        const str_t name = md_atom_name(&sys->atom, i);
        if (str_eq_cstr(name, "OXT") || str_eq_cstr(name, "OT2") || str_eq_cstr(name, "OC2") || str_eq_cstr(name, "O2")) return true;
    }
    return false;
}

// Whether the chain ends between the consecutive residues a and b of a protein chain whose backbones are not linked,
// rather than going on over a break. Only what is unambiguous counts, since a crystal structure rarely has the
// hydrogens and often leaves out the OXT: a free carboxyl group on a (O and OXT), a free amino group on b (two or
// three hydrogens on its N), or b a lone amino acid with a free carboxyl group of its own: a free amino acid (the
// arginine bound to 1LAF), which a PDB file can give the chain id of the chain and the next number. The beads of a
// trace (coarse grained) tell none of this.
static bool protein_chain_ends_between(const md_system_t* sys, size_t a, size_t b, const md_amino_acid_atoms_t* atoms_a, const md_amino_acid_atoms_t* atoms_b) {
    if (atoms_a->n < 0) return false;
    if (amino_acid_carboxyl_free(sys, a, atoms_a->c)) return true;
    if (amino_acid_carboxyl_free(sys, b, atoms_b->c)) return true;
    bool has_bonds;
    return bonded_count_of_element(sys, atoms_b->n, MD_Z_H, &has_bonds) >= 2;
}

// Whether the chain ends between the consecutive residues of a nucleic acid chain whose backbones are not linked: b
// has no P (a 5' end, which a crystal structure gives the first residue of a strand), or a hydrogen on the O3' of a or
// the O5' of b makes a free 3' or 5' hydroxyl. Beads tell only the first: a b without its phosphate bead.
static bool nucleic_chain_ends_between(const md_system_t* sys, const md_nucleic_acid_atoms_t* atoms_a, const md_nucleic_acid_atoms_t* atoms_b) {
    if (nucleic_atoms_are_beads(atoms_a)) return atoms_b->p < 0;
    if (atoms_a->o3 < 0 || atoms_b->p < 0) return true;
    bool has_bonds;
    return bonded_count_of_element(sys, atoms_a->o3, MD_Z_H, &has_bonds) > 0 || bonded_count_of_element(sys, atoms_b->o5, MD_Z_H, &has_bonds) > 0;
}

// Whether a polymer chain ends between the consecutive residues a and b = a + 1, rather than going on over a break, by
// what their atoms or beads tell (CHAIN BREAKS in md_util.h): the amino acids by protein_chain_ends_between, the
// nucleotides of atoms or of beads by nucleic_chain_ends_between. Whether their backbones are linked is not asked.
static bool polymer_chain_ends_between(const md_system_t* sys, size_t a, size_t b) {
    const md_urange_t range_a = md_system_residue_atom_range(sys, a);
    const md_urange_t range_b = md_system_residue_atom_range(sys, b);
    md_amino_acid_atoms_t aa_a, aa_b;
    md_nucleic_acid_atoms_t na_a, na_b;
    if (md_util_amino_acid_atoms_extract(&aa_a, &sys->atom, range_a) && md_util_amino_acid_atoms_extract(&aa_b, &sys->atom, range_b)) {
        return protein_chain_ends_between(sys, a, b, &aa_a, &aa_b);
    }
    if (md_util_nucleic_acid_atoms_extract(&na_a, &sys->atom, range_a) && md_util_nucleic_acid_atoms_extract(&na_b, &sys->atom, range_b)) {
        return nucleic_chain_ends_between(sys, &na_a, &na_b);
    }
    if (nucleotide_backbone_beads(&na_a, sys, a) && nucleotide_backbone_beads(&na_b, sys, b)) {
        return nucleic_chain_ends_between(sys, &na_a, &na_b);
    }
    return false;
}

// Marks the break between the consecutive residues a and b = a + 1 of a chain. The chain goes on past the break, so
// neither residue ends it: what the classification took for a terminus there (the carboxyl group of a residue whose
// successor is missing, the amino group of one whose predecessor is, a hydroxyl on either side of a nucleic acid
// break) is not one, and loses its terminal flags.
static void mark_chain_break(md_system_t* sys, size_t a, size_t b) {
    ASSERT(b == a + 1 && b < sys->residue.count);
    sys->residue.flags[a] = (md_residue_flags_t)(((uint32_t)sys->residue.flags[a] | MD_RESIDUE_FLAG_BREAK_AFTER)          & ~(uint32_t)MD_RESIDUE_FLAG_TERMINAL_END);
    sys->residue.flags[b] = (md_residue_flags_t)(((uint32_t)sys->residue.flags[b] | MD_RESIDUE_FLAG_BREAK_BEFORE) & ~(uint32_t)MD_RESIDUE_FLAG_TERMINAL_BEG);
    if (sys->atom.flags) {
        const md_urange_t range_a = md_system_residue_atom_range(sys, a);
        const md_urange_t range_b = md_system_residue_atom_range(sys, b);
        for (uint32_t i = range_a.beg; i < range_a.end; ++i) {
            sys->atom.flags[i] = (md_atom_flags_t)((uint32_t)sys->atom.flags[i] & ~(uint32_t)MD_ATOM_FLAG_TERMINAL_END);
        }
        for (uint32_t i = range_b.beg; i < range_b.end; ++i) {
            sys->atom.flags[i] = (md_atom_flags_t)((uint32_t)sys->atom.flags[i] & ~(uint32_t)MD_ATOM_FLAG_TERMINAL_BEG);
        }
    }
}

bool md_util_chain_break_atoms(md_atom_idx_t out_atoms[2], const md_system_t* sys, size_t res_idx) {
    ASSERT(out_atoms);
    out_atoms[0] = -1;
    out_atoms[1] = -1;
    if (!sys || res_idx + 1 >= sys->residue.count) return false;
    if (!(md_system_residue_flags(sys, res_idx) & MD_RESIDUE_FLAG_BREAK_AFTER)) return false;

    // The atoms the backbone was found by (md_util_system_infer): an amino acid by its atom names, else by its backbone
    // bead, a nucleotide by its atom names, else by its backbone beads
    const size_t next = res_idx + 1;
    const md_urange_t range_a = md_system_residue_atom_range(sys, res_idx);
    const md_urange_t range_b = md_system_residue_atom_range(sys, next);
    md_atom_idx_t a = -1;
    md_atom_idx_t b = -1;
    md_amino_acid_atoms_t aa_a, aa_b;
    md_nucleic_acid_atoms_t na_a, na_b;
    if (md_util_amino_acid_atoms_extract(&aa_a, &sys->atom, range_a) && md_util_amino_acid_atoms_extract(&aa_b, &sys->atom, range_b)) {
        a = aa_a.c;
        b = aa_b.n;
    } else if (amino_acid_backbone_bead(sys, res_idx) >= 0 && amino_acid_backbone_bead(sys, next) >= 0) {
        a = amino_acid_backbone_bead(sys, res_idx);
        b = amino_acid_backbone_bead(sys, next);
    } else if (md_util_nucleic_acid_atoms_extract(&na_a, &sys->atom, range_a) && md_util_nucleic_acid_atoms_extract(&na_b, &sys->atom, range_b)) {
        a = na_a.o3;
        b = na_b.p;
    } else if (nucleotide_backbone_beads(&na_a, sys, res_idx) && nucleotide_backbone_beads(&na_b, sys, next)) {
        a = na_a.c1;
        b = na_b.p;
    }
    if (a < 0 || b < 0) return false;
    out_atoms[0] = a;
    out_atoms[1] = b;
    return true;
}

static inline void flush_protein_backbone_segment(md_amino_acid_atoms_t* backbone_atoms, md_residue_idx_t* res_base, size_t part_idx, md_protein_backbone_data_t* bb_data, md_allocator_i* alloc, size_t min_backbone_length) {
    size_t backbone_length = md_array_size(backbone_atoms);
    if (backbone_length >= min_backbone_length && *res_base != -1) {
        commit_protein_backbone(backbone_atoms, backbone_length, *res_base, (md_part_idx_t)part_idx, bb_data, alloc);
    }
    md_array_shrink(backbone_atoms, 0);
    *res_base = -1;
}

// As commit_protein_backbone
static inline void commit_nucleic_backbone(md_nucleic_acid_atoms_t* bb_atoms, size_t bb_length, md_residue_idx_t res_base, md_part_idx_t part_idx, md_nucleic_backbone_data_t* bb_data, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(bb_data->segment.atoms);
    md_array_push_array(bb_data->segment.atoms, bb_atoms, bb_length, alloc);
    md_array_push(bb_data->range.offset, offset, alloc);
    md_array_push(bb_data->range.part_idx, part_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        const md_residue_idx_t res_idx = res_base < 0 ? -1 : res_base + (md_residue_idx_t)i;
        md_array_push(bb_data->segment.res_idx, res_idx, alloc);
    }
}

static inline void flush_nucleic_backbone_segment(md_nucleic_acid_atoms_t* backbone_atoms, md_residue_idx_t* res_base, size_t part_idx, md_nucleic_backbone_data_t* bb_data, md_allocator_i* alloc, size_t min_backbone_length) {
    size_t backbone_length = md_array_size(backbone_atoms);
    if (backbone_length >= min_backbone_length && *res_base != -1) {
        commit_nucleic_backbone(backbone_atoms, backbone_length, *res_base, (md_part_idx_t)part_idx, bb_data, alloc);
    }
    md_array_shrink(backbone_atoms, 0);
    *res_base = -1;
}

// ### BACKBONES FROM THE BONDS ###
//
// A system without residues (an XYZ file, the geometry of a QM calculation) has no residues whose atoms could be found
// by their names. Its backbones are found from the covalent bonds alone (md_match within each structure), which gives
// what a cartoon and the secondary structure need: the atoms of each segment, in order along each chain.
//
//   AMINO ACIDS  N-CA-C-O, the bond orders not tested: they are md_chem's perception of bonds which were themselves
//                inferred from distances, and one wrong bond (a clash between two chains taken for a bond) moves the
//                double bond of a carbonyl onto its C-N and loses the residue. Without orders N-CA-CB-OG of a serine or a
//                threonine matches too, as does N-CD-CG-OD1 of a hydroxyproline. The matches are taken best first, each
//                where its N, CA and C are not taken yet: one whose C is bonded to the N of another match (a peptide
//                bond), then one whose C has two oxygens (the carboxylate of a C terminus), then one whose C-O is as short
//                as a double bond (which tells the carbonyl of a C terminus without its OXT from a hydroxyl). CB is the
//                heavy atom on CA besides N and C (none for a glycine), HN the hydrogen of an N which has one.
//   NUCLEOTIDES  The sugar with its backbone, O5'-C5'-C4'(-O4'-C1'-C2'-C3'-)C3'-O3', and the P bonded to the O5'. A
//                sugar matches in one way only.
//   CHAINS       A residue follows the one whose C (O3') is bonded to its N (P). A chain runs from a residue which follows
//                none; what is left are rings (a cyclic peptide), cut before their first residue. A chain shorter than
//                min_length is no backbone: that leaves out the free amino acids and nucleotides (ATP, GTP) and the
//                groups which match without being residues (N7-C5-C6=O6 of a guanine, the ethanolamine of a lipid).
//
// The segments have no residue (res_idx -1) and the ranges no part (part_idx -1). The chains come in the order
// of their first residue, the residues by their CA (C4').

// Ångström below which a bond from C to O is double: C=O is 1.23 long, C-O 1.43
#define CARBONYL_MAX_DIST 1.33f

static inline bool atoms_closer_than(const md_system_state_t* state, md_atom_idx_t a, md_atom_idx_t b, float max_dist) {
    if (!md_system_state_has_coords(state) || (size_t)a >= state->num_atoms || (size_t)b >= state->num_atoms) return false;
    vec3_t d = vec3_sub(md_state_coord(state, (size_t)b), md_state_coord(state, (size_t)a));
    md_util_min_image_vec3(&d, 1, &state->unitcell);
    return vec3_dot(d, d) < max_dist * max_dist;
}

// The first neighbour of atom a of element z, which is neither x nor y (z 0: any element but hydrogen)
static md_atom_idx_t bonded_atom_other(const md_system_t* sys, md_atom_idx_t a, md_atomic_number_t z, md_atom_idx_t x, md_atom_idx_t y) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)a);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const md_atom_idx_t b = md_bond_iter_atom_index(&it);
        if (b == x || b == y) continue;
        const md_atomic_number_t zb = md_atom_atomic_number(&sys->atom, (size_t)b);
        if (z ? zb == z : zb != MD_Z_H) return b;
    }
    return -1;
}

// The matches of a query of the backbone atoms in the system, at the structure level, bonds of any order
static bool backbone_matches(md_match_result_t* out, str_t smiles, const md_system_t* sys, md_allocator_i* alloc) {
    md_match_query_t query = {0};
    if (!md_match_query_init_smiles(&query, smiles, alloc, NULL)) return false;
    const md_match_desc_t desc = {
        .query       = &query,
        .level       = MD_MATCH_LEVEL_MOLECULE,
        .mode        = MD_MATCH_MODE_UNIQUE,
        .hydrogens   = MD_MATCH_RESOLVE_NEVER,
        .bond_orders = MD_MATCH_RESOLVE_NEVER,
        .charges     = MD_MATCH_RESOLVE_NEVER,
    };
    return md_match_find(out, &desc, sys, alloc) && out->count > 0;
}

static int cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// The order of the chains of num residues, from the residue each follows (prev, -1 for none) and is followed by (next):
// out_order lists the residues chain by chain, a chain beginning at each residue of out_begin (num + 1 entries each).
// Returns the number of chains.
static size_t backbone_chains(uint32_t* out_order, uint32_t* out_begin, const int32_t* prev, const int32_t* next, size_t num, md_temp_scope_t temp) {
    uint8_t* done = md_temp_alloc_zero_array(temp, uint8_t, num);
    size_t num_order = 0;
    size_t num_chains = 0;
    // Open chains from their first residue, then the rings that are left
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t u = 0; u < num; ++u) {
            if (done[u] || (pass == 0 && prev[u] >= 0)) continue;
            out_begin[num_chains++] = (uint32_t)num_order;
            for (int32_t v = (int32_t)u; v >= 0 && !done[v]; v = next[v]) {
                done[v] = 1;
                out_order[num_order++] = (uint32_t)v;
            }
        }
    }
    out_begin[num_chains] = (uint32_t)num_order;
    return num_chains;
}

static void protein_backbones_from_bonds(md_system_t* sys, const md_system_state_t* state, size_t min_length, md_temp_scope_t temp) {
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_match_result_t res = {0};
    // Columns: N CA C O
    if (!backbone_matches(&res, STR_LIT("NCCO"), sys, temp_alloc)) return;
    ASSERT(res.width == 4);

    const size_t num_atoms = sys->atom.count;
    const size_t num_matches = res.count;

    // The score of each match, its key: score descending, then the order found
    uint8_t* is_n = md_temp_alloc_zero_array(temp, uint8_t, num_atoms);
    for (size_t i = 0; i < num_matches; ++i) {
        is_n[res.atom_idx[i * 4 + 0]] = 1;
    }
    uint64_t* key = md_temp_alloc_array(temp, uint64_t, num_matches);
    for (size_t i = 0; i < num_matches; ++i) {
        const int32_t* m = res.atom_idx + i * 4;
        uint32_t score = 0;
        uint32_t num_o = 0;
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)m[2]);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_atom_idx_t b = md_bond_iter_atom_index(&it);
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, (size_t)b);
            if (z == MD_Z_N && b != m[0] && is_n[b]) score |= 4;
            if (z == MD_Z_O) num_o += 1;
        }
        if (num_o >= 2) score |= 2;
        if (atoms_closer_than(state, m[2], m[3], CARBONYL_MAX_DIST)) score |= 1;
        key[i] = (uint64_t)(7 - score) << 32 | i;
    }
    qsort(key, num_matches, sizeof(uint64_t), cmp_u64);

    // Taken best first where N, CA and C are free; the residue whose N, CA or C an atom is, -1 for none
    int32_t* residue_of = md_temp_alloc_array(temp, int32_t, num_atoms);
    MEMSET(residue_of, -1, sizeof(int32_t) * num_atoms);
    uint64_t* by_ca = md_temp_alloc_array(temp, uint64_t, num_matches);
    size_t num = 0;
    for (size_t k = 0; k < num_matches; ++k) {
        const uint32_t i = (uint32_t)key[k];
        const int32_t* m = res.atom_idx + i * 4;
        if (residue_of[m[0]] >= 0 || residue_of[m[1]] >= 0 || residue_of[m[2]] >= 0) continue;
        residue_of[m[0]] = residue_of[m[1]] = residue_of[m[2]] = 0;
        by_ca[num++] = (uint64_t)m[1] << 32 | i;
    }
    qsort(by_ca, num, sizeof(uint64_t), cmp_u64);

    md_amino_acid_atoms_t* res_atoms = md_temp_alloc_array(temp, md_amino_acid_atoms_t, num);
    for (size_t u = 0; u < num; ++u) {
        const int32_t* m = res.atom_idx + (uint32_t)by_ca[u] * 4;
        md_amino_acid_atoms_t* a = &res_atoms[u];
        a->n  = m[0];
        a->ca = m[1];
        a->c  = m[2];
        a->o  = m[3];
        a->cb = bonded_atom_other(sys, a->ca, 0, a->n, a->c);
        a->hn = bonded_atom_other(sys, a->n, MD_Z_H, -1, -1);
        if (a->hn >= 0 && bonded_atom_other(sys, a->n, MD_Z_H, a->hn, -1) >= 0) a->hn = -1;
        residue_of[a->n] = residue_of[a->ca] = residue_of[a->c] = (int32_t)u;
    }

    // Residue v follows u when the C of u is bonded to the N of v
    int32_t* prev = md_temp_alloc_array(temp, int32_t, num);
    int32_t* next = md_temp_alloc_array(temp, int32_t, num);
    MEMSET(prev, -1, sizeof(int32_t) * num);
    MEMSET(next, -1, sizeof(int32_t) * num);
    for (size_t u = 0; u < num; ++u) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)res_atoms[u].c);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_atom_idx_t b = md_bond_iter_atom_index(&it);
            const int32_t v = residue_of[b];
            if (v >= 0 && v != (int32_t)u && res_atoms[v].n == b && prev[v] < 0) {
                next[u] = v;
                prev[v] = (int32_t)u;
                break;
            }
        }
    }

    uint32_t* order = md_temp_alloc_array(temp, uint32_t, num);
    uint32_t* begin = md_temp_alloc_array(temp, uint32_t, num + 1);
    const size_t num_chains = backbone_chains(order, begin, prev, next, num, temp);

    md_amino_acid_atoms_t* chain = md_temp_alloc_array(temp, md_amino_acid_atoms_t, num);
    for (size_t c = 0; c < num_chains; ++c) {
        const size_t len = begin[c + 1] - begin[c];
        if (len < min_length) continue;
        for (size_t j = 0; j < len; ++j) {
            chain[j] = res_atoms[order[begin[c] + j]];
        }
        commit_protein_backbone(chain, len, -1, -1, &sys->protein_backbone, sys->alloc);
    }
}

static void nucleic_backbones_from_bonds(md_system_t* sys, size_t min_length, md_temp_scope_t temp) {
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_match_result_t res = {0};
    // Columns: O5' C5' C4' O4' C1' C2' C3' O3'
    if (!backbone_matches(&res, STR_LIT("OCC1OCCC1O"), sys, temp_alloc)) return;
    ASSERT(res.width == 8);

    const size_t num_atoms = sys->atom.count;

    // A sugar matches in one way only; the sugars by their C4'
    uint64_t* by_c4 = md_temp_alloc_array(temp, uint64_t, res.count);
    for (size_t i = 0; i < res.count; ++i) {
        by_c4[i] = (uint64_t)res.atom_idx[i * 8 + 2] << 32 | i;
    }
    qsort(by_c4, res.count, sizeof(uint64_t), cmp_u64);

    // The residue whose P an atom is, -1 for none
    int32_t* residue_of_p = md_temp_alloc_array(temp, int32_t, num_atoms);
    MEMSET(residue_of_p, -1, sizeof(int32_t) * num_atoms);
    md_nucleic_acid_atoms_t* res_atoms = md_temp_alloc_array(temp, md_nucleic_acid_atoms_t, res.count);
    size_t num = 0;
    for (size_t k = 0; k < res.count; ++k) {
        if (k > 0 && by_c4[k] >> 32 == by_c4[k - 1] >> 32) continue;
        const int32_t* m = res.atom_idx + (uint32_t)by_c4[k] * 8;
        md_nucleic_acid_atoms_t* a = &res_atoms[num];
        a->o5 = m[0];
        a->c5 = m[1];
        a->c4 = m[2];
        a->o4 = m[3];
        a->c1 = m[4];
        a->c2 = m[5];
        a->c3 = m[6];
        a->o3 = m[7];
        a->p  = bonded_atom_other(sys, a->o5, MD_Z_P, -1, -1);
        if (a->p >= 0 && residue_of_p[a->p] < 0) residue_of_p[a->p] = (int32_t)num;
        num += 1;
    }

    // Residue v follows u when the O3' of u is bonded to the P of v
    int32_t* prev = md_temp_alloc_array(temp, int32_t, num);
    int32_t* next = md_temp_alloc_array(temp, int32_t, num);
    MEMSET(prev, -1, sizeof(int32_t) * num);
    MEMSET(next, -1, sizeof(int32_t) * num);
    for (size_t u = 0; u < num; ++u) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)res_atoms[u].o3);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const int32_t v = residue_of_p[md_bond_iter_atom_index(&it)];
            if (v >= 0 && v != (int32_t)u && prev[v] < 0) {
                next[u] = v;
                prev[v] = (int32_t)u;
                break;
            }
        }
    }

    uint32_t* order = md_temp_alloc_array(temp, uint32_t, num);
    uint32_t* begin = md_temp_alloc_array(temp, uint32_t, num + 1);
    const size_t num_chains = backbone_chains(order, begin, prev, next, num, temp);

    md_nucleic_acid_atoms_t* chain = md_temp_alloc_array(temp, md_nucleic_acid_atoms_t, num);
    for (size_t c = 0; c < num_chains; ++c) {
        const size_t len = begin[c + 1] - begin[c];
        if (len < min_length) continue;
        for (size_t j = 0; j < len; ++j) {
            chain[j] = res_atoms[order[begin[c] + j]];
        }
        commit_nucleic_backbone(chain, len, -1, -1, &sys->nucleic_backbone, sys->alloc);
    }
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

// ### LIPIDS ###
//
// md_util_system_infer_lipids: the heads and tails of the lipids (md_lipid_data_t). A molecule is a set of residues
// which are not amino acids, nucleotides, water or ions, bonded to each other, wherever they are in the system. Its apolar atoms
// are the carbons bonded to carbons (and hydrogens) only, or, for coarse grained beads, the beads named as the beads of
// the tails of the Martini lipidomes (lipid_bead_kind). Each connected part of them which hangs off the rest of the
// molecule at one atom (its anchor) is a tail: a chain, or the skeleton of a sterol when it holds fused rings (or the
// ring beads of a coarse grained sterol). A chain which hangs between two anchors spans them (a bolalipid). The rest of the molecule, in connected parts, is its heads.
//
// A molecule whose residues are lipids by name is a lipid with the tails found. One of no lipid is a lipid by its
// structure: a polar head (an atom other than carbon outside of the tails) and a chain without rings of at least
// LIPID_MIN_CHAIN_UNNAMED carbons, or the skeleton of a sterol. Coarse grained molecules are lipids by the force field
// only: by their names, by their beads where the name is not Martini's (md_util_system_infer_residue_flags), or by beads
// named as those of the tails of the Martini lipidomes where the molecule is none of the force field's.
//
// The lipid of a molecule of several residues is the part of it which holds the tails, with what is bonded to that
// part (lipid_extent): a sugar chain or a polymer attached to a lipid, each sugar or monomer a residue of its own, is no
// part of the lipid.

#define LIPID_MIN_CHAIN_NAMED       4   // Carbons on the spine of a chain of a lipid by name (hexanoyl: 5)
#define LIPID_MIN_CHAIN_UNNAMED     8   // Of a molecule which is a lipid by its structure alone (octanoyl: 7, decanoyl: 9)
#define LIPID_MIN_CHAIN_BEADS       1   // Beads on the spine of a chain of a coarse grained lipid (the C1B of a lyso PC)
#define LIPID_MIN_CHAIN_LETTERLESS  2   // Of a chain of beads named without the letter of their tail (COO C1 C2 C3 C4)
#define LIPID_MIN_SKELETON_ATOMS    15  // Carbons of the skeleton of a sterol (cholesterol: 26 besides C3)
#define LIPID_MIN_SKELETON_BEADS    4   // Beads of it (Martini cholesterol: 7 and 8)
#define LIPID_MAX_ATOMS             16384   // Larger molecules (polymers) are no lipids

// What a bead of a lipid of the Martini lipidomes is to its tails, by its name:
//   CHAIN       A bead of a tail: C (saturated), D (unsaturated) or T (trans), the number of the bead and the letter of
//               its tail, with the number of the acyl chain after it in cardiolipins (C1A, D2A, C4B, T1A, C1A1); the ring
//               skeleton of a sterol R and a number (R1 .. R6), and in a sterol its side chain C and a number (C1, C2)
//               and its dummies V (VOH V1 V2 V3): the optimal LINCS cholesterol of Martini 2 carries the mass of its
//               rings on dummies which do not interact, its ring beads virtual sites of them, so its skeleton holds
//               together through them, and its head is the hydroxyl bead ROH alone, as that of the cholesterol of 8 beads.
//   LETTERLESS  C, D or T and a number alone elsewhere: the chain of a fatty acid or a detergent of Martini 2 (COO C1 C2
//               C3 C4, SO3 C1 C2 C4, the B1 B2 B3 C1 C2 of an octyl glucoside), but also the inositol of a
//               phosphatidylinositol (C1 .. C4) and the sugars of glycolipids. Such beads are a tail only where they
//               make a chain: an unbranched run of two beads or more, without rings, which the inositol and sugars are.
//   POLAR       Any other: the beads of the heads (NC3, NH3, CNO, PO4, GL1, GL2, AM1, ROH, COO, ...).
typedef enum {
    LIPID_BEAD_POLAR = 0,
    LIPID_BEAD_CHAIN,
    LIPID_BEAD_LETTERLESS,
} lipid_bead_t;

static lipid_bead_t lipid_bead_kind(str_t name, bool sterol) {
    if (name.len < 2) return LIPID_BEAD_POLAR;
    const char c = name.ptr[0];
    if (c == 'V' && sterol) return LIPID_BEAD_CHAIN;
    if (c != 'C' && c != 'D' && c != 'T' && c != 'R') return LIPID_BEAD_POLAR;
    size_t i = 1;
    if (!is_digit(name.ptr[i])) return LIPID_BEAD_POLAR;
    while (i < name.len && is_digit(name.ptr[i])) ++i;
    if (i == name.len) {
        if (c == 'R' || (c == 'C' && sterol)) return LIPID_BEAD_CHAIN;
        return LIPID_BEAD_LETTERLESS;
    }
    if (name.ptr[i] < 'A' || name.ptr[i] > 'H') return LIPID_BEAD_POLAR;
    ++i;
    while (i < name.len && is_digit(name.ptr[i])) ++i;
    return i == name.len ? LIPID_BEAD_CHAIN : LIPID_BEAD_POLAR;
}

static inline bool lipid_atom_is_bead(const md_system_t* sys, size_t i) {
    return md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_ATOM && md_atom_atomic_number(&sys->atom, i) == 0;
}

// A heavy atom, or a bead: what a head or the spine of a tail is made of
static inline bool lipid_atom_is_heavy(const md_system_t* sys, size_t i) {
    return lipid_atom_is_bead(sys, i) || md_atom_atomic_number(&sys->atom, i) > MD_Z_H;
}

// A ring bead of a coarse grained sterol (R1 .. R6)
static inline bool lipid_bead_sterol_ring(str_t name) {
    return name.len >= 2 && name.ptr[0] == 'R' && is_digit(name.ptr[1]);
}

// What an atom is to the tails: a bead by its name, an atom a CHAIN when it is a carbon bonded to carbons and hydrogens only
static lipid_bead_t lipid_atom_kind(const md_system_t* sys, size_t i, bool sterol) {
    if (lipid_atom_is_bead(sys, i)) {
        return lipid_bead_kind(md_atom_name(&sys->atom, i), sterol);
    }
    if (md_atom_atomic_number(&sys->atom, i) != MD_Z_C) return LIPID_BEAD_POLAR;
    md_bond_iter_t it = md_bond_iter(&sys->bond, i);
    while (md_bond_iter_has_next(&it)) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(&it));
        if (z != MD_Z_C && z != MD_Z_H) return LIPID_BEAD_POLAR;
        md_bond_iter_next(&it);
    }
    return LIPID_BEAD_CHAIN;
}

static inline bool lipid_candidate_res(const md_system_t* sys, size_t res_idx) {
    const md_residue_kind_t kind = md_system_residue_kind(sys, res_idx);
    return kind == MD_RESIDUE_KIND_OTHER || kind == MD_RESIDUE_KIND_LIPID;
}

// A bead named as one of a tail with the letter of its tail (C1A, D2B, T1A, C1A1), as only the lipids of the Martini
// lipidomes have them
static inline bool lipid_bead_lettered(str_t name) {
    return name.len >= 3 && name.ptr[0] != 'R' && lipid_bead_kind(name, false) == LIPID_BEAD_CHAIN;
}

// The residues of the molecules are joined by the bonds between them (union-find, each set rooted at its first
// residue)
static int32_t lipid_find(int32_t* root, int32_t c) {
    while (root[c] != c) {
        root[c] = root[root[c]];
        c = root[c];
    }
    return c;
}

static void lipid_union(int32_t* root, int32_t a, int32_t b) {
    a = lipid_find(root, a);
    b = lipid_find(root, b);
    if (a < b) root[b] = a;
    else if (b < a) root[a] = b;
}

// A molecule: residues bonded to each other, ascending, and their atoms. Its residues need not be consecutive (the
// lipid A of a PDB entry, whose acyl chains are ligands of their own listed after its sugars).
typedef struct {
    const int32_t* residues; // [num_residues] ascending
    int32_t  num_residues;
    int32_t* atoms;         // [n] ascending: the atoms of its residues one after the other
    int32_t* member;        // [n] the residue (index into residues) of each atom
    int32_t  n;
    bool     contiguous;    // atoms[k] == atoms[0] + k
} lipid_mol_t;

// The index within the molecule of an atom, -1 for one not of it
static inline int32_t lipid_mol_local(const lipid_mol_t* mol, int32_t atom) {
    if (mol->contiguous) {
        const int32_t k = atom - mol->atoms[0];
        return (0 <= k && k < mol->n) ? k : -1;
    }
    int32_t lo = 0, hi = mol->n;
    while (lo < hi) {
        const int32_t mid = (lo + hi) / 2;
        if (mol->atoms[mid] < atom) lo = mid + 1;
        else hi = mid;
    }
    return (lo < mol->n && mol->atoms[lo] == atom) ? lo : -1;
}

typedef struct {
    int32_t  beg, end;      // its atoms (local indices) in part_atoms
    int32_t  anchor[2];     // atoms (global) of the rest of the molecule it hangs off
    int32_t  num_anchors;   // distinct anchors, more than 2 counts as 3
    bool     outside;       // hangs off an atom outside of the molecule
    bool     skeleton;      // holds an atom of a fused ring, or a ring bead of a sterol (R1 .. R6)
    bool     letterless;    // all of its beads are named without the letter of a tail (LIPID_BEAD_LETTERLESS)
    bool     lettered;      // all of its atoms are beads named with the letter of a tail (lipid_bead_lettered)
    bool     branched;      // an atom of it is bonded to more than two others of it
    int32_t  bonds;         // bonds within it, each counted from both ends
} lipid_part_t;

static int lipid_cmp_i32(const void* a, const void* b) {
    const int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

static int lipid_cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// Appends a range of indices to a ragged array of offsets (the first of them 0) and indices
static void push_index_range(uint32_t** offset, int32_t** idx, const int32_t* src, size_t n, md_allocator_i* alloc) {
    if (md_array_size(*offset) == 0) md_array_push(*offset, 0, alloc);
    if (n > 0) md_array_push_array(*idx, src, n, alloc);
    md_array_push(*offset, (uint32_t)md_array_size(*idx), alloc);
}

// The empty lipid data, whose offsets start at 0
static void init_lipid_data(md_lipid_data_t* lipid, md_allocator_i* alloc) {
    md_array_push(lipid->residue_offset, 0, alloc);
    md_array_push(lipid->head_offset, 0, alloc);
    md_array_push(lipid->tail_offset, 0, alloc);
    md_array_push(lipid->head.offset, 0, alloc);
    md_array_push(lipid->tail.spine_offset, 0, alloc);
    md_array_push(lipid->tail.atom_offset, 0, alloc);
}

// Records a lipid of the residues residues (ascending), with the heads and tails pushed for it since the last one (none
// when its anatomy is not known), and its residues as lipids
static void record_lipid(md_system_t* sys, const int32_t* residues, size_t num_residues) {
    md_lipid_data_t* lipid = &sys->lipid;
    md_allocator_i* alloc = sys->alloc;
    push_index_range(&lipid->residue_offset, &lipid->residue_idx, residues, num_residues, alloc);
    md_array_push(lipid->head_offset, (uint32_t)lipid->head.count, alloc);
    md_array_push(lipid->tail_offset, (uint32_t)lipid->tail.count, alloc);
    lipid->count += 1;
    for (size_t i = 0; i < num_residues; ++i) {
        set_residue_kind(sys, (size_t)residues[i], MD_RESIDUE_KIND_LIPID);
    }
}

// The lipid of a molecule of several residues: the residues which hold its tails and those between them (the core),
// and the residues bonded to the core: a head group which is a residue of its own (the PC of a lyso PC of the AMBER
// lipids), the first sugar of a glycolipid, the first Kdo of a lipid A. What hangs off that (the further sugars of a
// glycolipid, the core and O antigen of a lipopolysaccharide, the polymer of a PEG lipid of several residues) is attached
// to the lipid but not of it. The core is what is left of the residues when those without tails are taken off the ends
// one after the other. in_lipid[c] for each residue of the molecule.
static void lipid_extent(uint8_t* in_lipid, const md_system_t* sys, const lipid_mol_t* mol, const uint8_t* has_tail, md_allocator_i* temp_arena) {
    const int32_t nc = mol->num_residues;
    MEMSET(in_lipid, 1, (size_t)nc);
    if (nc < 2) return;

    // The bonds between the residues, each pair once
    md_array(uint64_t) edges = 0;
    for (int32_t k = 0; k < mol->n; ++k) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, mol->atoms[k]);
        while (md_bond_iter_has_next(&it)) {
            const int32_t l = lipid_mol_local(mol, md_bond_iter_atom_index(&it));
            md_bond_iter_next(&it);
            if (l < 0) continue;
            const uint32_t a = (uint32_t)mol->member[k], b = (uint32_t)mol->member[l];
            if (a < b) md_array_push(edges, ((uint64_t)a << 32) | b, temp_arena);
        }
    }
    size_t num_edges = md_array_size(edges);
    if (num_edges == 0) return;
    qsort(edges, num_edges, sizeof(uint64_t), lipid_cmp_u64);
    size_t u = 0;
    for (size_t e = 0; e < num_edges; ++e) {
        if (u == 0 || edges[e] != edges[u - 1]) edges[u++] = edges[e];
    }
    num_edges = u;

    int32_t* deg = md_alloc(temp_arena, (size_t)nc * sizeof(int32_t));
    int32_t* off = md_alloc(temp_arena, ((size_t)nc + 1) * sizeof(int32_t));
    int32_t* adj = md_alloc(temp_arena, 2 * num_edges * sizeof(int32_t));
    int32_t* queue = md_alloc(temp_arena, (size_t)nc * sizeof(int32_t));
    uint8_t* alive = md_alloc(temp_arena, (size_t)nc);
    MEMSET(deg, 0, (size_t)nc * sizeof(int32_t));
    for (size_t e = 0; e < num_edges; ++e) {
        deg[edges[e] >> 32] += 1;
        deg[edges[e] & 0xFFFFFFFF] += 1;
    }
    off[0] = 0;
    for (int32_t c = 0; c < nc; ++c) off[c + 1] = off[c] + deg[c];
    MEMSET(deg, 0, (size_t)nc * sizeof(int32_t));
    for (size_t e = 0; e < num_edges; ++e) {
        const int32_t a = (int32_t)(edges[e] >> 32), b = (int32_t)(edges[e] & 0xFFFFFFFF);
        adj[off[a] + deg[a]++] = b;
        adj[off[b] + deg[b]++] = a;
    }

    // The core: the ends without tails taken off until none is left
    int32_t qe = 0;
    for (int32_t c = 0; c < nc; ++c) {
        alive[c] = 1;
        if (!has_tail[c] && deg[c] <= 1) {
            alive[c] = 0;
            queue[qe++] = c;
        }
    }
    for (int32_t qb = 0; qb < qe; ++qb) {
        const int32_t c = queue[qb];
        for (int32_t e = off[c]; e < off[c + 1]; ++e) {
            const int32_t d = adj[e];
            if (!alive[d]) continue;
            deg[d] -= 1;
            if (!has_tail[d] && deg[d] <= 1) {
                alive[d] = 0;
                queue[qe++] = d;
            }
        }
    }

    // And what is bonded to it
    for (int32_t c = 0; c < nc; ++c) {
        bool in = alive[c];
        for (int32_t e = off[c]; e < off[c + 1] && !in; ++e) in = alive[adj[e]];
        in_lipid[c] = (uint8_t)in;
    }
}

// The lipid of a molecule, if it is or holds one
static void infer_lipid_molecule(md_system_t* sys, const lipid_mol_t* mol, const uint8_t* in_fused_ring, md_allocator_i* temp_arena) {
    md_lipid_data_t* lipid = &sys->lipid;
    md_allocator_i* alloc = sys->alloc;

    const int32_t n = mol->n;
    const int32_t* atoms = mol->atoms;
    const int32_t nc = mol->num_residues;
    if (n <= 0) return;

    // Coarse grained when it holds a bead (a molecule too large to look into: when the first atom of a residue is one)
    bool cg = false;
    bool sterol_beads = false;
    if (atoms) {
        for (int32_t k = 0; k < n; ++k) {
            if (lipid_atom_is_bead(sys, atoms[k])) {
                cg = true;
                sterol_beads |= lipid_bead_sterol_ring(md_atom_name(&sys->atom, atoms[k]));
            }
        }
    } else {
        for (int32_t c = 0; c < nc && !cg; ++c) {
            const md_urange_t r = md_system_residue_atom_range(sys, (size_t)mol->residues[c]);
            cg = r.beg < r.end && lipid_atom_is_bead(sys, r.beg);
        }
    }
    // Lipids by name. Those which are lipids by their structure alone are classified anew, so that inferring again gives
    // what inferring once did. A coarse grained lipid is one of the force field, by its name or, where that is not
    // Martini's, by its beads (the kind md_util_system_infer_residue_flags gave it), or one whose beads are named as those
    // of the lipids of the force field (below): never by its structure alone.
    bool named = false;
    for (int32_t c = 0; c < nc; ++c) {
        const size_t res = (size_t)mol->residues[c];
        if (md_util_resname_lipid(md_system_residue_name(sys, res)) || (cg && md_system_residue_kind(sys, res) == MD_RESIDUE_KIND_LIPID)) {
            named = true;
        } else if (md_system_residue_kind(sys, res) == MD_RESIDUE_KIND_LIPID) {
            set_residue_kind(sys, res, MD_RESIDUE_KIND_OTHER);
        }
    }
    // A lipid by name too small or too large to look into is one without heads and tails
    if (n <= 1 || n > LIPID_MAX_ATOMS || !atoms) {
        if (named) record_lipid(sys, mol->residues, (size_t)nc);
        return;
    }

    const size_t pos = md_vm_arena_get_pos(temp_arena);

    uint8_t* apolar     = md_alloc(temp_arena, n);
    uint8_t* letterless = md_alloc(temp_arena, n);                      // apolar, if it makes a chain (LIPID_BEAD_LETTERLESS)
    int32_t* part       = md_alloc(temp_arena, n * sizeof(int32_t));   // the apolar part of an atom, -1
    int32_t* part_atoms = md_alloc(temp_arena, n * sizeof(int32_t));   // the atoms of each part, contiguous
    int32_t* dist       = md_alloc(temp_arena, n * sizeof(int32_t));
    int32_t* prev       = md_alloc(temp_arena, n * sizeof(int32_t));
    int32_t* queue      = md_alloc(temp_arena, n * sizeof(int32_t));
    int32_t* head       = md_alloc(temp_arena, n * sizeof(int32_t));   // the head of a heavy atom outside of the tails, -1
    uint8_t* tail       = md_alloc(temp_arena, n);                      // in a tail
    uint8_t* has_tail   = md_alloc(temp_arena, nc);                     // a residue holds an atom of a tail
    uint8_t* in_lipid   = md_alloc(temp_arena, nc);                     // a residue is of the lipid (lipid_extent)
    int32_t* lipid_residues = md_alloc(temp_arena, nc * sizeof(int32_t));
    md_array(lipid_part_t) parts = 0;

    for (int32_t k = 0; k < n; ++k) {
        const lipid_bead_t kind = lipid_atom_is_heavy(sys, atoms[k]) ? lipid_atom_kind(sys, atoms[k], sterol_beads) : LIPID_BEAD_POLAR;
        apolar[k] = (uint8_t)(kind != LIPID_BEAD_POLAR);
        letterless[k] = (uint8_t)(kind == LIPID_BEAD_LETTERLESS);
        part[k] = -1;
        head[k] = -1;
        tail[k] = 0;
    }
    MEMSET(has_tail, 0, nc);

    // The connected apolar parts and what they hang off. The atoms of a part are part_atoms[beg, end), in the order
    // they were reached.
    int32_t num_part_atoms = 0;
    for (int32_t s = 0; s < n; ++s) {
        if (!apolar[s] || part[s] != -1) continue;
        const int32_t id = (int32_t)md_array_size(parts);
        lipid_part_t p = { .beg = num_part_atoms, .anchor = {-1, -1}, .letterless = true, .lettered = true };
        int32_t qb = num_part_atoms;
        part_atoms[num_part_atoms++] = s;
        part[s] = id;
        while (qb < num_part_atoms) {
            const int32_t k = part_atoms[qb++];
            int32_t degree = 0;
            p.letterless &= letterless[k] != 0;
            const bool bead = lipid_atom_is_bead(sys, atoms[k]);
            p.lettered &= bead && lipid_bead_lettered(md_atom_name(&sys->atom, atoms[k]));
            if (in_fused_ring && in_fused_ring[atoms[k]]) p.skeleton = true;
            // The ring beads of a sterol, whose bonds (a star of virtual sites in Martini 3) need not form rings
            if (bead && lipid_bead_sterol_ring(md_atom_name(&sys->atom, atoms[k]))) p.skeleton = true;
            md_bond_iter_t it = md_bond_iter(&sys->bond, atoms[k]);
            while (md_bond_iter_has_next(&it)) {
                const int32_t j = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (!lipid_atom_is_heavy(sys, j)) continue;
                const int32_t l = lipid_mol_local(mol, j);
                if (l < 0) {
                    p.outside = true;
                    continue;
                }
                if (apolar[l]) {
                    p.bonds += 1;
                    degree += 1;
                    if (part[l] == -1) {
                        part[l] = id;
                        part_atoms[num_part_atoms++] = l;
                    }
                } else if (j != p.anchor[0] && j != p.anchor[1]) {
                    if (p.num_anchors < 2) p.anchor[p.num_anchors] = j;
                    p.num_anchors = MIN(p.num_anchors + 1, 3);
                }
            }
            p.branched |= degree > 2;
        }
        p.end = num_part_atoms;
        md_array_push(parts, p, temp_arena);
    }

    // A coarse grained molecule the force field does not know, from a topology of its own (a ganglioside, a sphingomyelin
    // named otherwise): a lipid when a chain of beads named as those of the tails of the Martini lipidomes (C1A C2A C3A,
    // the letter of the tail after the number) hangs off the rest of the molecule, which has a bead of its own. Then it is
    // looked into as a lipid by name, so that inferring again (the kind lipid now its name) gives the same.
    if (cg && !named) {
        bool lettered_tail = false;
        const size_t num_parts = md_array_size(parts);
        for (size_t pi = 0; pi < num_parts && !lettered_tail; ++pi) {
            const lipid_part_t* p = &parts[pi];
            const int32_t size = p->end - p->beg;
            lettered_tail = p->lettered && !p->skeleton && !p->outside && !p->branched && p->num_anchors >= 1 && p->num_anchors <= 2 &&
                            size >= 2 && p->bonds / 2 == size - 1;
        }
        bool polar = false;
        for (int32_t k = 0; k < n && !polar; ++k) polar = !apolar[k] && lipid_atom_is_heavy(sys, atoms[k]);
        if (!lettered_tail || !polar) goto done;
        named = true;
    }

    // The tails among the parts, and their spines. The spine is found by a breadth first search within the part: from
    // every atom bonded to the anchor out to the farthest atom, or for a spanning chain to the nearest atom bonded to its
    // other anchor.
    typedef struct { int32_t beg, end; } spine_range_t;
    md_array(int32_t) spines = 0;
    md_array(spine_range_t) spine_ranges = 0;
    md_array(int32_t) tail_part = 0;
    md_array(md_lipid_tail_kind_t) tail_kind = 0;
    md_array(int32_t) tail_anchor_beg = 0;
    md_array(int32_t) tail_anchor_end = 0;

    const size_t num_parts = md_array_size(parts);
    for (size_t pi = 0; pi < num_parts; ++pi) {
        const lipid_part_t* p = &parts[pi];
        if (p->outside || p->num_anchors == 0 || p->num_anchors > 2) continue;
        const int32_t size = p->end - p->beg;
        const bool skeleton = p->skeleton;
        if (skeleton) {
            // The skeleton of a sterol hangs off its hydroxyl (or ester) alone
            if (p->num_anchors != 1 || size < (cg ? LIPID_MIN_SKELETON_BEADS : LIPID_MIN_SKELETON_ATOMS)) continue;
        } else if (p->letterless && (p->bonds / 2 != size - 1 || p->branched || size < LIPID_MIN_CHAIN_LETTERLESS)) {
            // Beads named without the letter of a tail are one only as a chain: the inositol of a phosphatidylinositol and
            // the sugars of a glycolipid, named so too, are rings
            continue;
        } else if (!named && p->bonds / 2 != size - 1) {
            // The chains of a lipid of no name hold no rings: a ring with alkyl groups (a drug, a dye) is no tail. Those of
            // lipids by name may (the cyclopropane of a fatty acid, the cyclopentanes of the chains of archaeal lipids).
            continue;
        }

        // Anchors in order of their index, so that a spanning chain runs from the head of the lower one
        int32_t anchor_beg = p->anchor[0];
        int32_t anchor_end = p->num_anchors == 2 ? p->anchor[1] : -1;
        if (anchor_end != -1 && anchor_end < anchor_beg) {
            const int32_t t = anchor_beg; anchor_beg = anchor_end; anchor_end = t;
        }

        int32_t qb = 0, qe = 0;
        for (int32_t m = p->beg; m < p->end; ++m) {
            const int32_t k = part_atoms[m];
            prev[k] = -1;
            dist[k] = -1;
            if (md_bond_find(&sys->bond, atoms[k], anchor_beg) >= 0) {
                dist[k] = 0;
                queue[qe++] = k;
            }
        }
        int32_t tip = -1;
        while (qb < qe) {
            const int32_t k = queue[qb++];
            if (anchor_end != -1) {
                if (md_bond_find(&sys->bond, atoms[k], anchor_end) >= 0) {
                    tip = k;
                    break;
                }
            } else if (tip == -1 || dist[k] > dist[tip] || (dist[k] == dist[tip] && k > tip)) {
                // The farthest atom, the last of them on a tie (the end of the side chain of a Martini cholesterol)
                tip = k;
            }
            md_bond_iter_t it = md_bond_iter(&sys->bond, atoms[k]);
            while (md_bond_iter_has_next(&it)) {
                const int32_t l = lipid_mol_local(mol, md_bond_iter_atom_index(&it));
                md_bond_iter_next(&it);
                if (l >= 0 && part[l] == (int32_t)pi && dist[l] == -1) {
                    dist[l] = dist[k] + 1;
                    prev[l] = k;
                    queue[qe++] = l;
                }
            }
        }
        if (tip == -1) continue;

        const int32_t len = dist[tip] + 1;
        const int32_t min_len = skeleton ? 1 : (cg ? (p->letterless ? LIPID_MIN_CHAIN_LETTERLESS : LIPID_MIN_CHAIN_BEADS) : (named ? LIPID_MIN_CHAIN_NAMED : LIPID_MIN_CHAIN_UNNAMED));
        if (len < min_len) continue;

        const int32_t base = (int32_t)md_array_size(spines);
        spine_range_t r = { base, base + len };
        md_array_resize(spines, (size_t)r.end, temp_arena);
        int32_t at = r.end - 1;
        for (int32_t k = tip; k != -1; k = prev[k]) {
            spines[at--] = atoms[k];
        }
        ASSERT(at == base - 1);
        md_array_push(spine_ranges, r, temp_arena);
        md_array_push(tail_part, (int32_t)pi, temp_arena);
        md_array_push(tail_kind, skeleton ? MD_LIPID_TAIL_SKELETON : MD_LIPID_TAIL_CHAIN, temp_arena);
        md_array_push(tail_anchor_beg, anchor_beg, temp_arena);
        md_array_push(tail_anchor_end, anchor_end, temp_arena);
        for (int32_t m = p->beg; m < p->end; ++m) {
            tail[part_atoms[m]] = 1;
            has_tail[mol->member[part_atoms[m]]] = 1;
        }
    }

    // A lipid by name in which no tails are found (no bonds, a naming the beads do not follow) is one without heads and
    // tails: it is a lipid molecule still
    const size_t num_tails = md_array_size(spine_ranges);
    if (num_tails == 0) {
        if (named) record_lipid(sys, mol->residues, (size_t)nc);
        goto done;
    }

    // The residues of the lipid, and of what is attached to it
    lipid_extent(in_lipid, sys, mol, has_tail, temp_arena);

    // A molecule of no lipid needs a polar head: an atom other than carbon outside of its tails. Its head is of the
    // elements of lipids (C, N, O, P, S) and is not aromatic: a porphyrin with a phytyl or farnesyl tail (chlorophyll,
    // pheophytin, heme A) and the like are cofactors to be shown as such.
    if (!named) {
        bool polar = false;
        for (int32_t k = 0; k < n; ++k) {
            if (tail[k] || !in_lipid[mol->member[k]]) continue;
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, atoms[k]);
            if (z <= MD_Z_H) continue;
            if (z != MD_Z_C && z != MD_Z_N && z != MD_Z_O && z != MD_Z_P && z != MD_Z_S) goto done;
            if (md_atom_flags(&sys->atom, atoms[k]) & MD_ATOM_FLAG_AROMATIC) goto done;
            polar |= z != MD_Z_C;
        }
        if (!polar) goto done;
    }

    // The heads: the connected parts of the heavy atoms of the lipid outside of the tails, each in ascending order
    const uint32_t head_base = (uint32_t)lipid->head.count;
    uint32_t num_heads = 0;
    for (int32_t s = 0; s < n; ++s) {
        if (tail[s] || head[s] != -1 || !in_lipid[mol->member[s]] || !lipid_atom_is_heavy(sys, atoms[s])) continue;
        int32_t qb = 0, qe = 0;
        queue[qe++] = s;
        head[s] = (int32_t)num_heads;
        while (qb < qe) {
            const int32_t k = queue[qb++];
            md_bond_iter_t it = md_bond_iter(&sys->bond, atoms[k]);
            while (md_bond_iter_has_next(&it)) {
                const int32_t j = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                const int32_t l = lipid_mol_local(mol, j);
                if (l < 0) continue;
                if (!tail[l] && head[l] == -1 && in_lipid[mol->member[l]] && lipid_atom_is_heavy(sys, j)) {
                    head[l] = (int32_t)num_heads;
                    queue[qe++] = l;
                }
            }
        }
        for (int32_t m = 0; m < qe; ++m) queue[m] = atoms[queue[m]];
        qsort(queue, (size_t)qe, sizeof(int32_t), lipid_cmp_i32);
        push_index_range(&lipid->head.offset, &lipid->head.atom_idx, queue, (size_t)qe, alloc);
        num_heads += 1;
    }

    for (size_t t = 0; t < num_tails; ++t) {
        const spine_range_t r = spine_ranges[t];
        push_index_range(&lipid->tail.spine_offset, &lipid->tail.spine_idx, spines + r.beg, (size_t)(r.end - r.beg), alloc);

        // Its atoms and their hydrogens, ascending
        const lipid_part_t* p = &parts[tail_part[t]];
        int32_t m = 0;
        for (int32_t q = p->beg; q < p->end; ++q) {
            const int32_t i = atoms[part_atoms[q]];
            queue[m++] = i;
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            while (md_bond_iter_has_next(&it)) {
                const int32_t j = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (lipid_mol_local(mol, j) >= 0 && md_atom_atomic_number(&sys->atom, j) == MD_Z_H && m < n) {
                    queue[m++] = j;
                }
            }
        }
        qsort(queue, (size_t)m, sizeof(int32_t), lipid_cmp_i32);
        {
            // A hydrogen bonded to two atoms of the tail (an inferred bond too many) is in it once
            int32_t u = 0;
            for (int32_t q = 0; q < m; ++q) {
                if (u == 0 || queue[q] != queue[u - 1]) queue[u++] = queue[q];
            }
            m = u;
        }
        push_index_range(&lipid->tail.atom_offset, &lipid->tail.atom_idx, queue, (size_t)m, alloc);
        for (int32_t q = 0; q < m; ++q) {
            sys->atom.flags[queue[q]] |= MD_ATOM_FLAG_LIPID_TAIL;
        }

        const int32_t ab = tail_anchor_beg[t];
        const int32_t ae = tail_anchor_end[t];
        const uint32_t hb = head_base + (uint32_t)head[lipid_mol_local(mol, ab)];
        const uint32_t he = ae != -1 ? head_base + (uint32_t)head[lipid_mol_local(mol, ae)] : UINT32_MAX;
        md_array_push(lipid->tail.anchor_beg, ab, alloc);
        md_array_push(lipid->tail.anchor_end, ae, alloc);
        md_array_push(lipid->tail.head_beg, hb, alloc);
        md_array_push(lipid->tail.head_end, he, alloc);
        md_array_push(lipid->tail.kind, tail_kind[t], alloc);
    }

    lipid->head.count += num_heads;
    lipid->tail.count += num_tails;

    {
        // What is attached to the lipid is no lipid, whatever its name
        int32_t num_lipid_residues = 0;
        for (int32_t c = 0; c < nc; ++c) {
            if (in_lipid[c]) {
                lipid_residues[num_lipid_residues++] = mol->residues[c];
            } else if (md_system_residue_kind(sys, (size_t)mol->residues[c]) == MD_RESIDUE_KIND_LIPID) {
                set_residue_kind(sys, (size_t)mol->residues[c], MD_RESIDUE_KIND_OTHER);
            }
        }
        record_lipid(sys, lipid_residues, (size_t)num_lipid_residues);
    }

done:
    md_vm_arena_set_pos_back(temp_arena, pos);
}

bool md_util_system_infer_lipids(md_system_t* sys) {
    ASSERT(sys);
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    md_allocator_i* alloc = sys->alloc;
    md_system_topology_changed(sys);

    // Derived from scratch
    md_lipid_data_free(&sys->lipid, alloc);
    md_lipid_data_t* lipid = &sys->lipid;
    init_lipid_data(lipid, alloc);

    if (md_array_size(sys->atom.flags) < sys->atom.count) {
        const size_t size = md_array_size(sys->atom.flags);
        md_array_resize(sys->atom.flags, sys->atom.count, alloc);
        MEMSET(sys->atom.flags + size, 0, (sys->atom.count - size) * sizeof(md_atom_flags_t));
    }
    for (size_t i = 0; i < sys->atom.count; ++i) sys->atom.flags[i] &= ~MD_ATOM_FLAG_LIPID_TAIL;
    // Without bonds every lipid by name is a molecule of its own residue, without heads and tails
    if (sys->atom.count == 0 || sys->residue.count == 0) {
        return true;
    }
    if (md_array_size(sys->residue.flags) < sys->residue.count) {
        const size_t size = md_array_size(sys->residue.flags);
        md_array_resize(sys->residue.flags, sys->residue.count, alloc);
        MEMSET(sys->residue.flags + size, 0, (sys->residue.count - size) * sizeof(md_residue_flags_t));
    }

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The atoms of fused rings: a skeleton (a sterol) rather than a chain
    uint8_t* in_fused_ring = NULL;
    const size_t num_rings = md_ring_data_count(&sys->ring);
    for (size_t r = 0; r < num_rings; ++r) {
        if (!(md_system_ring_flags(sys, r) & MD_RING_FLAG_FUSED)) continue;
        if (!in_fused_ring) {
            in_fused_ring = md_alloc(temp_arena, sys->atom.count);
            MEMSET(in_fused_ring, 0, sys->atom.count);
        }
        const int32_t* atoms = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        for (size_t k = 0; k < n; ++k) in_fused_ring[atoms[k]] = 1;
    }

    // The molecules: the residues which may be of a lipid (not amino acids, nucleotides, water or ions), joined by the
    // bonds between them, wherever they are in the system
    const size_t num_residues = sys->residue.count;
    int32_t* root = md_alloc(temp_arena, num_residues * sizeof(int32_t));
    for (size_t c = 0; c < num_residues; ++c) {
        root[c] = lipid_candidate_res(sys, c) ? (int32_t)c : -1;
    }
    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] < 0) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        for (uint32_t i = range.beg; i < range.end; ++i) {
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            while (md_bond_iter_has_next(&it)) {
                const uint32_t j = (uint32_t)md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (range.beg <= j && j < range.end) continue;
                const md_residue_idx_t d = md_system_residue_find_by_atom_idx(sys, j);
                if (d < 0 || root[d] < 0) continue;
                lipid_union(root, (int32_t)c, (int32_t)d);
            }
        }
    }

    // The residues of each molecule, ascending, after one another
    int32_t* count = md_alloc(temp_arena, num_residues * sizeof(int32_t));
    int32_t* start = md_alloc(temp_arena, num_residues * sizeof(int32_t));
    MEMSET(count, 0, num_residues * sizeof(int32_t));
    size_t num_members = 0;
    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] < 0) continue;
        root[c] = lipid_find(root, (int32_t)c);
        count[root[c]] += 1;
        num_members += 1;
    }
    int32_t* members = md_alloc(temp_arena, MAX(num_members, 1) * sizeof(int32_t));
    int32_t pos = 0;
    for (size_t c = 0; c < num_residues; ++c) {
        start[c] = pos;
        if (root[c] == (int32_t)c) pos += count[c];
    }
    MEMSET(count, 0, num_residues * sizeof(int32_t));
    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] < 0) continue;
        const int32_t r = root[c];
        members[start[r] + count[r]++] = (int32_t)c;
    }

    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] != (int32_t)c) continue;
        const int32_t* residues = members + start[c];
        const int32_t nc = count[c];

        size_t n = 0;
        for (int32_t m = 0; m < nc; ++m) {
            const md_urange_t r = md_system_residue_atom_range(sys, (size_t)residues[m]);
            n += r.end - r.beg;
        }
        lipid_mol_t mol = { .residues = residues, .num_residues = nc, .n = (int32_t)MIN(n, (size_t)INT32_MAX) };
        const size_t mark = md_vm_arena_get_pos(temp_arena);
        if (n <= LIPID_MAX_ATOMS) {
            mol.atoms  = md_alloc(temp_arena, MAX(n, 1) * sizeof(int32_t));
            mol.member = md_alloc(temp_arena, MAX(n, 1) * sizeof(int32_t));
            int32_t k = 0;
            for (int32_t m = 0; m < nc; ++m) {
                const md_urange_t r = md_system_residue_atom_range(sys, (size_t)residues[m]);
                for (uint32_t i = r.beg; i < r.end; ++i) {
                    mol.atoms[k] = (int32_t)i;
                    mol.member[k] = m;
                    k += 1;
                }
            }
            mol.contiguous = n > 0 && mol.atoms[n - 1] - mol.atoms[0] == (int32_t)n - 1;
        }
        infer_lipid_molecule(sys, &mol, in_fused_ring, temp_arena);
        md_vm_arena_set_pos_back(temp_arena, mark);
    }

    md_temp_end(temp);
    return true;
}

// ### RESIDUE TYPES ###

// The color a new component type starts out with: that of its standard monomer, that of the element of a monatomic
// ion, and otherwise one of its name, the same in every system
static uint32_t residue_type_default_color(const md_system_t* sys, size_t res_idx, str_t name, md_residue_kind_t kind, md_monomer_t monomer) {
    if (monomer != MD_MONOMER_NONE) {
        return md_monomer_color(monomer);
    }
    if (kind == MD_RESIDUE_KIND_ION) {
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        const md_atomic_number_t z = range.end > range.beg ? md_atom_atomic_number(&sys->atom, range.beg) : 0;
        if (z) return md_atomic_number_cpk_color(z);
    }
    // The golden ratio spreads the hues of names which hash alike
    const uint32_t h = (uint32_t)md_hash64_str(name, 0) * 0x9E3779B9u;
    const float hue = (float)(h >> 8) * (1.0f / 16777216.0f);
    const vec3_t rgb = hcl_to_rgb(hue, 0.5f, 0.7f);    // Channels within [0.45, 0.95]
    const uint32_t r = (uint32_t)(CLAMP(rgb.x, 0.0f, 1.0f) * 255.0f + 0.5f);
    const uint32_t g = (uint32_t)(CLAMP(rgb.y, 0.0f, 1.0f) * 255.0f + 0.5f);
    const uint32_t b = (uint32_t)(CLAMP(rgb.z, 0.0f, 1.0f) * 255.0f + 0.5f);
    return (0xFFu << 24) | (b << 16) | (g << 8) | r;
}

static inline uint64_t residue_type_key(str_t name, md_residue_kind_t kind) {
    return md_hash64_str(name, (uint64_t)kind + 1);
}

bool md_util_system_infer_residue_types(md_system_t* sys) {
    ASSERT(sys);
    md_residue_data_t* res = &sys->residue;
    if (res->count == 0) return false;
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    md_allocator_i* alloc = sys->alloc;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The types the system has already, so that inferring again keeps them and what was made of them (their colors).
    // A type no residue is of any more (a residue which was reclassified) stays, used by none.
    md_hashmap32_t lookup = { .allocator = temp_arena };
    for (size_t t = 0; t < res->type.count; ++t) {
        md_hashmap_add(&lookup, residue_type_key(LBL_TO_STR(res->type.name[t]), res->type.kind[t]), (uint32_t)t);
    }

    md_array_resize(res->type_idx, res->count, alloc);
    for (size_t c = 0; c < res->count; ++c) {
        const str_t name = md_system_residue_name(sys, c);
        const md_residue_kind_t kind = md_system_residue_kind(sys, c);
        const uint64_t key = residue_type_key(name, kind);

        md_residue_type_idx_t type_idx = -1;
        const uint32_t* found = md_hashmap_get(&lookup, key);
        if (found && res->type.kind[*found] == kind && str_eq(LBL_TO_STR(res->type.name[*found]), name)) {
            type_idx = (md_residue_type_idx_t)*found;
        } else if (found) {
            type_idx = md_residue_type_find(&sys->residue.type, name, kind);       // Two keys which collide
        }
        if (type_idx < 0) {
            const md_monomer_t monomer = md_monomer_from_residue_name(name, kind);
            const uint32_t color = residue_type_default_color(sys, c, name, kind, monomer);
            type_idx = md_residue_type_add(&sys->residue.type, name, kind, monomer, color, alloc);
            if (!found) md_hashmap_add(&lookup, key, (uint32_t)type_idx);
        }
        res->type_idx[c] = type_idx;
    }

    md_temp_end(temp);
    return true;
}

// Try to fill in missing fields for molecule struct
// (Coordinates & Elements) -> Covalent Bonds
// (residues & Bonds)       -> Chains
// (Chains)                 -> Backbone
bool md_util_system_infer(md_system_t* sys, const md_system_state_t* state, md_infer_flags_t flags) {
    ASSERT(sys);

    if (sys->atom.count == 0) return false;
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    md_system_topology_changed(sys);

    md_allocator_i* alloc = sys->alloc;

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    const bool cg = md_system_is_coarse_grained(sys);
    size_t num_particle_types = md_system_particle_type_count(sys);

    if (cg) {
        flags &= ~MD_UTIL_INFER_BOND_BIT;
    }

    if (flags & MD_UTIL_INFER_COLOR_BIT) {
        if (sys->atom.type.color) {
            // @Note: The 0'th type is a sentinel type for representing invalid or unset atoms.
            // Hence it's skipped here
            for (size_t i = 1; i < num_particle_types; ++i) {
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
                        color = (uint32_t)(uint8_t)(alpha * 255) << 24 | (uint32_t)(uint8_t)(rgb.x * 255) << 16 | (uint32_t)(uint8_t)(rgb.y * 255) << 8 | (uint32_t)(uint8_t)(rgb.z * 255);
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
   
    if (flags & MD_UTIL_INFER_BOND_BIT) {
        // Setting the bit means "derive this", unconditionally. A format which supplies its own
        // bonds (LAMMPS, PSF) is respected by the caller not setting the bit, not by a guard here.
        md_util_system_infer_covalent_bonds(sys, state);
    }

    // Whatever the bonds came from (inferred, a topology, a file, by hand): those between a metal and a non-metal are
    // coordination. Cheap, and every reader of the flags relies on it, so not behind a bit.
    md_util_system_infer_coordination(sys);

    // The parts before the molecules: a molecule holds its polymer strands whole (synthesize_molecule_links). They are
    // found from the instances the data gives, the bonds and the classified residues, which the molecules do not change.
    if (flags & MD_UTIL_INFER_PART_BIT) {
        md_util_system_infer_chains(sys);
        md_util_system_infer_parts(sys);
    }

    if (flags & MD_UTIL_INFER_MOLECULE_BIT) {
        // Coarse grained systems get a hierarchy even without bonds (synthesize_hierarchy_links), and a residue is a
        // molecule whatever its bonds (synthesize_molecule_links): ions and waters given without bonds are molecules
        if (sys->bond.count || cg || sys->residue.count) {
            md_util_system_infer_molecules(sys);
        }
        if (sys->bond.count) {
            md_util_system_infer_rings(sys);
        }
    }
    
    if (flags & MD_UTIL_INFER_CHEMISTRY_BIT) {
        // After the rings, before the hydrogen bond roles which read the orders. Protonation at pH 7 applies to
        // residues without hydrogens only.
        if (sys->bond.count) {
            md_chem_perceive(sys, state, MD_CHEM_FLAG_PROTONATE_PH7);
        }
    }

    if (flags & MD_UTIL_INFER_BACKBONE_BIT) {
        static const size_t MIN_BACKBONE_LENGTH = 3;

        // Derived from scratch: a repeated inference must not append to the backbones of the last one
        md_array_shrink(sys->protein_backbone.range.offset, 0);
        md_array_shrink(sys->protein_backbone.range.part_idx, 0);
        md_array_shrink(sys->protein_backbone.segment.atoms, 0);
        md_array_shrink(sys->protein_backbone.segment.res_idx, 0);
        sys->protein_backbone.range.count = 0;
        sys->protein_backbone.segment.count = 0;
        md_array_shrink(sys->nucleic_backbone.range.offset, 0);
        md_array_shrink(sys->nucleic_backbone.range.part_idx, 0);
        md_array_shrink(sys->nucleic_backbone.segment.atoms, 0);
        md_array_shrink(sys->nucleic_backbone.segment.res_idx, 0);
        sys->nucleic_backbone.range.count = 0;
        sys->nucleic_backbone.segment.count = 0;

        if (sys->part.count && sys->atom.type_idx) {
            // The backbones are found within the strands: the parts of a polymer, in which consecutive ranges of
            // residues form backbones

            static const size_t MAX_BACKBONE_LENGTH = 1 << 15;

            // The chain breaks too (mark_chain_break), found anew with the backbones
            if (md_array_size(sys->residue.flags) < sys->residue.count) {
                const size_t size = md_array_size(sys->residue.flags);
                md_array_resize(sys->residue.flags, sys->residue.count, alloc);
                MEMSET(sys->residue.flags + size, 0, (sys->residue.count - size) * sizeof(md_residue_flags_t));
            }
            for (size_t i = 0; i < sys->residue.count; ++i) {
                sys->residue.flags[i] = (md_residue_flags_t)((uint32_t)sys->residue.flags[i] & ~(uint32_t)(MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
            }

            size_t temp_pos = md_vm_arena_get_pos(temp_arena);

            {
                md_array(md_amino_acid_atoms_t) backbone_atoms = 0;
                md_array_ensure(backbone_atoms, MAX_BACKBONE_LENGTH, temp_arena);
                md_residue_idx_t res_base = -1;
                md_sequence_id_t prev_seq_id = -1;

                for (size_t part_idx = 0; part_idx < sys->part.count; ++part_idx) {
                    if (md_system_part_kind(sys, part_idx) != MD_PART_KIND_PEPTIDE) {
                        continue;
                    }

                    md_urange_t range = md_system_part_residue_range(sys, part_idx);
                    for (size_t res_idx = range.beg; res_idx < range.end; ++res_idx) {
                        md_amino_acid_atoms_t atoms;
                        md_urange_t atom_range = md_system_residue_atom_range(sys, res_idx);
                        bool has_backbone_atoms = md_util_amino_acid_atoms_extract(&atoms, &sys->atom, atom_range);
                        if (!has_backbone_atoms) {
                            // A coarse grained residue: its backbone bead, a segment of a trace
                            const md_atom_idx_t bead = amino_acid_backbone_bead(sys, res_idx);
                            if (bead >= 0) {
                                MEMSET(&atoms, -1, sizeof(atoms));
                                atoms.ca = bead;
                                has_backbone_atoms = true;
                            }
                        }

                        if (has_backbone_atoms) {
                            md_sequence_id_t seq_id = md_system_residue_seq_id(sys, res_idx);
                            size_t backbone_length = md_array_size(backbone_atoms);

                            if (backbone_length > 0) {
                                // The residue before it in the chain. A range is atomistic or a trace throughout, and
                                // ends where the chain breaks: where the backbone is not linked, whatever the sequence
                                // numbers say (CHAIN BREAKS in md_util.h)
                                const md_amino_acid_atoms_t prev = *md_array_last(backbone_atoms);
                                const bool trace = atoms.n < 0;
                                if (trace != (prev.n < 0)) {
                                    flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                                } else if (!protein_backbone_linked(sys, state, &prev, &atoms, seq_id == prev_seq_id + 1)) {
                                    if (!protein_chain_ends_between(sys, res_idx - 1, res_idx, &prev, &atoms)) {
                                        mark_chain_break(sys, res_idx - 1, res_idx);
                                    }
                                    flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                                }
                            }

                            if (res_base == -1) {
                                res_base = (md_residue_idx_t)res_idx;
                            }
                            md_array_push(backbone_atoms, atoms, temp_arena);
                            prev_seq_id = seq_id;
                            continue;
                        }

                        flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                        prev_seq_id = -1;
                    }
                    // Possibly commit remainder of the chain
                    flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->protein_backbone, alloc, MIN_BACKBONE_LENGTH);
                    prev_seq_id = -1;
                }

                sys->protein_backbone.range.count = md_array_size(sys->protein_backbone.range.offset);
                if (sys->protein_backbone.range.count) {
                    // Add end offset
                    md_array_push(sys->protein_backbone.range.offset, (uint32_t)md_array_size(sys->protein_backbone.segment.atoms), alloc);
                }
                sys->protein_backbone.segment.count = md_array_size(sys->protein_backbone.segment.atoms);
            }
            md_vm_arena_set_pos_back(temp_arena, temp_pos);

            {
                md_array(md_nucleic_acid_atoms_t) backbone_atoms = 0;
                md_array_ensure(backbone_atoms, MAX_BACKBONE_LENGTH, temp_arena);
                md_residue_idx_t res_base = -1;
                md_sequence_id_t prev_seq_id = -1;

                for (size_t part_idx = 0; part_idx < sys->part.count; ++part_idx) {
                    if (!md_part_kind_is_nucleic_acid(md_system_part_kind(sys, part_idx))) {
                        continue;
                    }

                    md_urange_t range = md_system_part_residue_range(sys, part_idx);
                    for (md_residue_idx_t res_idx = (int)range.beg; res_idx < (int)range.end; ++res_idx) {
                        md_nucleic_acid_atoms_t atoms;
                        md_urange_t atom_range = md_system_residue_atom_range(sys, res_idx);
                        bool has_backbone_atoms = md_util_nucleic_acid_atoms_extract(&atoms, &sys->atom, atom_range);
                        if (!has_backbone_atoms) {
                            // A coarse grained nucleotide: its backbone beads
                            has_backbone_atoms = nucleotide_backbone_beads(&atoms, sys, (size_t)res_idx);
                        }

                        if (has_backbone_atoms) {
                            md_sequence_id_t seq_id = md_system_residue_seq_id(sys, res_idx);
                            size_t backbone_length = md_array_size(backbone_atoms);

                            if (backbone_length > 0) {
                                // The residue before it in the chain: a range is of atoms or of beads throughout, and
                                // ends where the chain breaks (CHAIN BREAKS in md_util.h); a residue without its P (or
                                // phosphate bead), which is linked to nothing before it, begins one
                                const md_nucleic_acid_atoms_t prev = *md_array_last(backbone_atoms);
                                if (nucleic_atoms_are_beads(&atoms) != nucleic_atoms_are_beads(&prev)) {
                                    flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                                } else if (!nucleic_backbone_linked(sys, state, &prev, &atoms, seq_id == prev_seq_id + 1)) {
                                    if (!nucleic_chain_ends_between(sys, &prev, &atoms)) {
                                        mark_chain_break(sys, (size_t)res_idx - 1, (size_t)res_idx);
                                    }
                                    flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                                }
                            }

                            if (res_base == -1) {
                                res_base = res_idx;
                            }
                            md_array_push(backbone_atoms, atoms, temp_arena);
                            prev_seq_id = seq_id;
                            continue;
                        }

                        flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                        prev_seq_id = -1;
                    }
                    // Possibly commit remainder of the chain
                    flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &sys->nucleic_backbone, alloc, MIN_BACKBONE_LENGTH);
                    prev_seq_id = -1;
                }

                sys->nucleic_backbone.range.count = md_array_size(sys->nucleic_backbone.range.offset);
                if (sys->nucleic_backbone.range.count) {
                    // Add end offset
                    md_array_push(sys->nucleic_backbone.range.offset, (uint32_t)md_array_size(sys->nucleic_backbone.segment.atoms), alloc);
                }
                sys->nucleic_backbone.segment.count = md_array_size(sys->nucleic_backbone.segment.atoms);
            }
            md_vm_arena_set_pos_back(temp_arena, temp_pos);
        } else if (sys->residue.count == 0 && sys->bond.count && sys->atom.type_idx) {
            // No residues to find them by: from the bonds alone (BACKBONES FROM THE BONDS)
            const size_t temp_pos = md_vm_arena_get_pos(temp_arena);
            protein_backbones_from_bonds(sys, state, MIN_BACKBONE_LENGTH, temp_scope);
            nucleic_backbones_from_bonds(sys, MIN_BACKBONE_LENGTH, temp_scope);
            md_vm_arena_set_pos_back(temp_arena, temp_pos);

            sys->protein_backbone.range.count = md_array_size(sys->protein_backbone.range.offset);
            if (sys->protein_backbone.range.count) {
                md_array_push(sys->protein_backbone.range.offset, (uint32_t)md_array_size(sys->protein_backbone.segment.atoms), alloc);
            }
            sys->protein_backbone.segment.count = md_array_size(sys->protein_backbone.segment.atoms);

            sys->nucleic_backbone.range.count = md_array_size(sys->nucleic_backbone.range.offset);
            if (sys->nucleic_backbone.range.count) {
                md_array_push(sys->nucleic_backbone.range.offset, (uint32_t)md_array_size(sys->nucleic_backbone.segment.atoms), alloc);
            }
            sys->nucleic_backbone.segment.count = md_array_size(sys->nucleic_backbone.segment.atoms);
        }
    }

    if (flags & MD_UTIL_INFER_LIPID_BIT) {
        // After the bonds, the rings (a sterol is told by its fused rings) and the classification of the residues
        md_util_system_infer_lipids(sys);
    }

    // Last of what depends on the classification of the residues, which the lipids are the last to change. Cheap,
    // and what a residue type is follows from the residues alone, so not behind a bit.
    md_util_system_infer_residue_types(sys);

#if 0
    if ((flags & MD_UTIL_INFER_UNWRAP_MOLECULE_BIT) && 
        (state->unitcell.flags != 0) && md_system_molecule_count(sys) > 0)
    {
        md_util_system_unwrap(sys);
    }
#endif

    // Record the state the topology above was inferred from. Written here, as part of performing
    // the inference, so it cannot disagree with what it produced.
    if (state != &sys->reference) {
        sys->reference.alloc = sys->alloc;
        md_system_state_copy(&sys->reference, state);
    }

    md_temp_end(temp_scope);
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
        const float px = *(const float*)((const char*)vx + byte_stride * i);
        const float py = *(const float*)((const char*)vy + byte_stride * i);
        const float pz = *(const float*)((const char*)vz + byte_stride * i);
        int x = (int)((px - minv[0]) * scale * 1023.f + 0.5f);
        int y = (int)((py - minv[1]) * scale * 1023.f + 0.5f);
        int z = (int)((pz - minv[2]) * scale * 1023.f + 0.5f);

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

// Packed coordinates: the bounds come from the packed AABB, which reads them as they lie, and the
// keys from each atom in turn.
void md_util_sort_spatial(uint32_t* source, const vec3_t* xyz, size_t count) {
    if (!source || !xyz || count == 0) return;

    md_temp_scope_t temp_scope = md_temp_begin();

    float minv[3], maxv[3];
    md_util_aabb_compute(minv, maxv, xyz, NULL, NULL, count);
    const float extent = MAX(maxv[0] - minv[0], MAX(maxv[1] - minv[1], maxv[2] - minv[2]));
    const float scale = extent == 0 ? 0.f : 1023.f / extent;

    // Morton order based on the position inside a unit cube
    uint32_t* keys = md_temp_alloc_array(temp_scope, uint32_t, count);
    for (size_t i = 0; i < count; ++i) {
        const int x = (int)((xyz[i].x - minv[0]) * scale + 0.5f);
        const int y = (int)((xyz[i].y - minv[1]) * scale + 0.5f);
        const int z = (int)((xyz[i].z - minv[2]) * scale + 0.5f);
        keys[i] = part1By2(x) | (part1By2(y) << 1) | (part1By2(z) << 2);
    }

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram10(hist, keys, count);

    uint32_t* scratch = md_temp_alloc_array(temp_scope, uint32_t, count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass 10-bit radix sort computes the resulting order into scratch
    radixPass10(source, scratch, keys, count, hist, 0);
    radixPass10(scratch, source, keys, count, hist, 1);
    radixPass10(source, scratch, keys, count, hist, 2);

    md_temp_end(temp_scope);
}

void md_util_sort_spatial_xyz(uint32_t* source, const float* xyz, size_t stride_in_bytes, size_t count) {
    if (!source || !xyz || count <= 0) return;

    md_temp_scope_t temp_scope = md_temp_begin();

    uint32_t* keys = md_temp_alloc_array(temp_scope, uint32_t, count);
    const float* base = (const float*)xyz;
    const float* x = base + 0;
    const float* y = base + 1;
    const float* z = base + 2;
    computeOrder(keys, x, y, z, count, stride_in_bytes);

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram10(hist, keys, count);

    uint32_t* scratch = md_temp_alloc_array(temp_scope, uint32_t, count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass radix sort computes the resulting order into scratch
    radixPass10(source, scratch, keys, count, hist, 0);
    radixPass10(scratch, source, keys, count, hist, 1);
    radixPass10(source, scratch, keys, count, hist, 2);

    md_temp_end(temp_scope);
}

void md_util_sort_radix_u32(uint32_t* data, size_t count) {
    if (!data || count <= 1) return;
    ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_sort_u32(data, count);
        return;
    }
    if (is_sorted_u32(data, count)) return;

    md_temp_scope_t temp_scope = md_temp_begin();
    sort_radix_u32(data, count, temp_scope);
    md_temp_end(temp_scope);
}

void md_util_sort_radix_u64(uint64_t* data, size_t count) {
    if (!data || count <= 1) return;
    ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_sort_u64(data, count);
        return;
    }
    if (is_sorted_u64(data, count)) return;

    md_temp_scope_t temp_scope = md_temp_begin();
    sort_radix_u64(data, count, temp_scope);
    md_temp_end(temp_scope);
}

void md_util_argsort_radix_u32(uint32_t* out_indices, const uint32_t* keys, size_t count) {
    if (count == 0) return;
    ASSERT(out_indices);
    ASSERT(keys);
    ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_argsort_u32(out_indices, keys, count);
        return;
    }
    if (is_sorted_u32(keys, count)) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    argsort_radix_u32(out_indices, keys, count, temp_scope);
    md_temp_end(temp_scope);
}

void md_util_argsort_radix_u64(uint32_t* out_indices, const uint64_t* keys, size_t count) {
    if (count == 0) return;
    ASSERT(out_indices);
    ASSERT(keys);
    ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_argsort_u64(out_indices, keys, count);
        return;
    }
    if (is_sorted_u64(keys, count)) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    argsort_radix_u64(out_indices, keys, count, temp_scope);
    md_temp_end(temp_scope);
}

#ifdef __cplusplus
}
#endif
