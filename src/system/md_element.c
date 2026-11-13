#include <md_element.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_format.h>
#include <md_monomer.h>

#include <math.h>

// The element tables, indexed by atomic number (md_element.h)
enum { Num_Elements = MD_Z_Og + 1 };

static const md_str_t element_symbols[] = {
    MD_STR_INIT("Xx"), MD_STR_INIT("H"),  MD_STR_INIT("He"), MD_STR_INIT("Li"), MD_STR_INIT("Be"), MD_STR_INIT("B"),  MD_STR_INIT("C"),  MD_STR_INIT("N"),  MD_STR_INIT("O"),  MD_STR_INIT("F"),  MD_STR_INIT("Ne"), MD_STR_INIT("Na"), MD_STR_INIT("Mg"), MD_STR_INIT("Al"), MD_STR_INIT("Si"), MD_STR_INIT("P"),  MD_STR_INIT("S"),  MD_STR_INIT("Cl"), MD_STR_INIT("Ar"), MD_STR_INIT("K"),  MD_STR_INIT("Ca"), MD_STR_INIT("Sc"), MD_STR_INIT("Ti"), MD_STR_INIT("V"),
    MD_STR_INIT("Cr"), MD_STR_INIT("Mn"), MD_STR_INIT("Fe"), MD_STR_INIT("Co"), MD_STR_INIT("Ni"), MD_STR_INIT("Cu"), MD_STR_INIT("Zn"), MD_STR_INIT("Ga"), MD_STR_INIT("Ge"), MD_STR_INIT("As"), MD_STR_INIT("Se"), MD_STR_INIT("Br"), MD_STR_INIT("Kr"), MD_STR_INIT("Rb"), MD_STR_INIT("Sr"), MD_STR_INIT("Y"),  MD_STR_INIT("Zr"), MD_STR_INIT("Nb"), MD_STR_INIT("Mo"), MD_STR_INIT("Tc"), MD_STR_INIT("Ru"), MD_STR_INIT("Rh"), MD_STR_INIT("Pd"), MD_STR_INIT("Ag"),
    MD_STR_INIT("Cd"), MD_STR_INIT("In"), MD_STR_INIT("Sn"), MD_STR_INIT("Sb"), MD_STR_INIT("Te"), MD_STR_INIT("I"),  MD_STR_INIT("Xe"), MD_STR_INIT("Cs"), MD_STR_INIT("Ba"), MD_STR_INIT("La"), MD_STR_INIT("Ce"), MD_STR_INIT("Pr"), MD_STR_INIT("Nd"), MD_STR_INIT("Pm"), MD_STR_INIT("Sm"), MD_STR_INIT("Eu"), MD_STR_INIT("Gd"), MD_STR_INIT("Tb"), MD_STR_INIT("Dy"), MD_STR_INIT("Ho"), MD_STR_INIT("Er"), MD_STR_INIT("Tm"), MD_STR_INIT("Yb"), MD_STR_INIT("Lu"),
    MD_STR_INIT("Hf"), MD_STR_INIT("Ta"), MD_STR_INIT("W"),  MD_STR_INIT("Re"), MD_STR_INIT("Os"), MD_STR_INIT("Ir"), MD_STR_INIT("Pt"), MD_STR_INIT("Au"), MD_STR_INIT("Hg"), MD_STR_INIT("Tl"), MD_STR_INIT("Pb"), MD_STR_INIT("Bi"), MD_STR_INIT("Po"), MD_STR_INIT("At"), MD_STR_INIT("Rn"), MD_STR_INIT("Fr"), MD_STR_INIT("Ra"), MD_STR_INIT("Ac"), MD_STR_INIT("Th"), MD_STR_INIT("Pa"), MD_STR_INIT("U"),  MD_STR_INIT("Np"), MD_STR_INIT("Pu"), MD_STR_INIT("Am"),
    MD_STR_INIT("Cm"), MD_STR_INIT("Bk"), MD_STR_INIT("Cf"), MD_STR_INIT("Es"), MD_STR_INIT("Fm"), MD_STR_INIT("Md"), MD_STR_INIT("No"), MD_STR_INIT("Lr"), MD_STR_INIT("Rf"), MD_STR_INIT("Db"), MD_STR_INIT("Sg"), MD_STR_INIT("Bh"), MD_STR_INIT("Hs"), MD_STR_INIT("Mt"), MD_STR_INIT("Ds"), MD_STR_INIT("Rg"), MD_STR_INIT("Cn"), MD_STR_INIT("Nh"), MD_STR_INIT("Fl"), MD_STR_INIT("Mc"), MD_STR_INIT("Lv"), MD_STR_INIT("Ts"), MD_STR_INIT("Og"),
};


static const md_str_t element_names[] = {
    MD_STR_INIT("Unknown"),     MD_STR_INIT("Hydrogen"),     MD_STR_INIT("Helium"),       MD_STR_INIT("Lithium"),     MD_STR_INIT("Beryllium"),   MD_STR_INIT("Boron"),         MD_STR_INIT("Carbon"),     MD_STR_INIT("Nitrogen"),   MD_STR_INIT("Oxygen"),
    MD_STR_INIT("Fluorine"),    MD_STR_INIT("Neon"),         MD_STR_INIT("Sodium"),       MD_STR_INIT("Magnesium"),   MD_STR_INIT("Aluminium"),   MD_STR_INIT("Silicon"),       MD_STR_INIT("Phosphorus"), MD_STR_INIT("Sulfur"),     MD_STR_INIT("Chlorine"),
    MD_STR_INIT("Argon"),       MD_STR_INIT("Potassium"),    MD_STR_INIT("Calcium"),      MD_STR_INIT("Scandium"),    MD_STR_INIT("Titanium"),    MD_STR_INIT("Vanadium"),      MD_STR_INIT("Chromium"),   MD_STR_INIT("Manganese"),  MD_STR_INIT("Iron"),
    MD_STR_INIT("Cobalt"),      MD_STR_INIT("Nickel"),       MD_STR_INIT("Copper"),       MD_STR_INIT("Zinc"),        MD_STR_INIT("Gallium"),     MD_STR_INIT("Germanium"),     MD_STR_INIT("Arsenic"),    MD_STR_INIT("Selenium"),   MD_STR_INIT("Bromine"),
    MD_STR_INIT("Krypton"),     MD_STR_INIT("Rubidium"),     MD_STR_INIT("Strontium"),    MD_STR_INIT("Yttrium"),     MD_STR_INIT("Zirconium"),   MD_STR_INIT("Niobium"),       MD_STR_INIT("Molybdenum"), MD_STR_INIT("Technetium"), MD_STR_INIT("Ruthenium"),
    MD_STR_INIT("Rhodium"),     MD_STR_INIT("Palladium"),    MD_STR_INIT("Silver"),       MD_STR_INIT("Cadmium"),     MD_STR_INIT("Indium"),      MD_STR_INIT("Tin"),           MD_STR_INIT("Antimony"),   MD_STR_INIT("Tellurium"),  MD_STR_INIT("Iodine"),
    MD_STR_INIT("Xenon"),       MD_STR_INIT("Caesium"),      MD_STR_INIT("Barium"),       MD_STR_INIT("Lanthanum"),   MD_STR_INIT("Cerium"),      MD_STR_INIT("Praseodymium"),  MD_STR_INIT("Neodymium"),  MD_STR_INIT("Promethium"), MD_STR_INIT("Samarium"),
    MD_STR_INIT("Europium"),    MD_STR_INIT("Gadolinium"),   MD_STR_INIT("Terbium"),      MD_STR_INIT("Dysprosium"),  MD_STR_INIT("Holmium"),     MD_STR_INIT("Erbium"),        MD_STR_INIT("Thulium"),    MD_STR_INIT("Ytterbium"),  MD_STR_INIT("Lutetium"),
    MD_STR_INIT("Hafnium"),     MD_STR_INIT("Tantalum"),     MD_STR_INIT("Tungsten"),     MD_STR_INIT("Rhenium"),     MD_STR_INIT("Osmium"),      MD_STR_INIT("Iridium"),       MD_STR_INIT("Platinum"),   MD_STR_INIT("Gold"),       MD_STR_INIT("Mercury"),
    MD_STR_INIT("Thallium"),    MD_STR_INIT("Lead"),         MD_STR_INIT("Bismuth"),      MD_STR_INIT("Polonium"),    MD_STR_INIT("Astatine"),    MD_STR_INIT("Radon"),         MD_STR_INIT("Francium"),   MD_STR_INIT("Radium"),     MD_STR_INIT("Actinium"),
    MD_STR_INIT("Thorium"),     MD_STR_INIT("Protactinium"), MD_STR_INIT("Uranium"),      MD_STR_INIT("Neptunium"),   MD_STR_INIT("Plutonium"),   MD_STR_INIT("Americium"),     MD_STR_INIT("Curium"),     MD_STR_INIT("Berkelium"),  MD_STR_INIT("Californium"),
    MD_STR_INIT("Einsteinium"), MD_STR_INIT("Fermium"),      MD_STR_INIT("Mendelevium"),  MD_STR_INIT("Nobelium"),    MD_STR_INIT("Lawrencium"),  MD_STR_INIT("Rutherfordium"), MD_STR_INIT("Dubnium"),    MD_STR_INIT("Seaborgium"), MD_STR_INIT("Bohrium"),
    MD_STR_INIT("Hassium"),     MD_STR_INIT("Meitnerium"),   MD_STR_INIT("Darmstadtium"), MD_STR_INIT("Roentgenium"), MD_STR_INIT("Copernicium"), MD_STR_INIT("Nihonium"),      MD_STR_INIT("Flerovium"),  MD_STR_INIT("Moscovium"),  MD_STR_INIT("Livermorium"),
    MD_STR_INIT("Tennessine"),  MD_STR_INIT("Oganesson"),                          
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

// https://dx.doi.org/10.1021/jp8111556 (Mantina et al. 2009) for the main group elements; the others are not from it
static const float element_vdw_radii[] = {
    1.00, 1.10, 1.40, 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, 2.30, 2.15, 2.05,
    2.05, 2.05, 2.05, 2.00, 2.00, 2.00, 2.10, 1.87, 2.11, 1.85, 1.90, 1.83, 2.02, 3.03, 2.49, 2.40, 2.30, 2.15, 2.10, 2.05, 2.05, 2.00, 2.05, 2.10,
    2.20, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16, 3.43, 2.68, 2.50, 2.48, 2.47, 2.45, 2.43, 2.42, 2.40, 2.38, 2.37, 2.35, 2.33, 2.32, 2.30, 2.28, 2.27,
    2.25, 2.20, 2.10, 2.05, 2.00, 2.00, 2.05, 2.10, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, 2.00, 2.30, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00
};

// http://chemistry.wikia.com/wiki/List_of_elements_by_atomic_mass
// last element padded with 295
static const float element_atomic_mass[] = {
    0,       1.00794,    4.002602, 6.941,     9.012182, 10.811,    12.0107,   14.0067,     15.9994, 18.9984032, 20.1797, 22.98976928,
    24.3050, 26.9815386, 28.0855,  30.973762, 32.065,   35.453,    39.948,    39.0983,     40.078,  44.955912,  47.867,  50.9415,
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

// The element properties of md_element.h, from the tables above
md_atomic_number_t md_atomic_number_from_symbol(md_str_t str, bool ignore_case) {
    if (str.len == 1 || str.len == 2) {
        for (md_atomic_number_t i = 0; i < MD_COUNTOF(element_symbols); ++i) {
            md_str_t sym = element_symbols[i];
            if (md_str_len(sym) != md_str_len(str)) continue;
            if (ignore_case) {
                if (md_str_eq_ignore_case(str, sym)) return i;   
            } else {
                if (md_str_eq(str, sym)) return i;
            }
        }
    }
    return 0;
}


size_t md_atomic_number_infer_from_mass_batch(md_atomic_number_t element[], const float mass[], size_t count) {
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
        md_atomic_number_t elem = 0;
        const float m = mass[i];

        // Use an epsilon which scales with the mass
        const float eps = 0.001f * fmaxf(1.0f, m * 0.01f);

        if (0.0f < m && m != 1.0f) {
            // The nearest within eps. The masses are not in the order of the elements throughout (Ar is heavier
            // than K, Co than Ni, Te than I), so all are looked at.
            float best = eps;
            for (size_t j = 1; j < MD_COUNTOF(element_atomic_mass); ++j) {
                const float d = fabsf(m - element_atomic_mass[j]);
                if (d < best) {
                    best = d;
                    elem = (md_atomic_number_t)j;
                }
            }
        }

        element[i] = elem;
        successful_matches += (size_t)(elem != 0);
    }

    return successful_matches;
}

md_str_t md_atomic_number_symbol(md_atomic_number_t element) {
    element = element < Num_Elements ? element : 0;
    return element_symbols[element];
}

md_str_t md_atomic_number_name(md_atomic_number_t element) {
    element = element < Num_Elements ? element : 0;
    return element_names[element];
}

float md_atomic_number_vdw_radius(md_atomic_number_t element) {
    element = element < Num_Elements ? element : 0;
    return element_vdw_radii[element];
}

static inline float element_covalent_radius(md_atomic_number_t element) {
    return element < Num_Elements ? element_covalent_radii_u8[element] * 0.01f : 0;
}

float md_atomic_number_covalent_radius(md_atomic_number_t element) {
    return element_covalent_radius(element);
}

float md_atomic_number_mass(md_atomic_number_t element) {
    element = element < Num_Elements ? element : 0;
    return element_atomic_mass[element];
}

int md_atomic_number_max_valence(md_atomic_number_t element) {
    element = element < Num_Elements ? element : 0;
    return (int)element_max_valence[element];
}

uint32_t md_atomic_number_cpk_color(md_atomic_number_t element) {
    element = element < Num_Elements ? element : 0;
    return element_cpk_colors[element];
}

// INFERENCE OF THE ELEMENT: from the names of an atom and its residue, or its mass

typedef struct {
    md_str_t str;
    md_atomic_number_t z;
} Mapping;


// Predefined mappings (applies to a combination from res_name:atom_name)
static const Mapping predefined_mappings[] = {
    // Common ions
    {MD_STR_INIT("NA:NA"), MD_Z_Na},
    {MD_STR_INIT("SOD:SOD"), MD_Z_Na},
    {MD_STR_INIT("SOD:NA"),MD_Z_Na},
    {MD_STR_INIT("K:K"),   MD_Z_K},
    {MD_STR_INIT("POT:K"), MD_Z_K},
    {MD_STR_INIT("POT:POT"), MD_Z_K},  // CHARMM names the ion as its residue
    {MD_STR_INIT("CL:CL"), MD_Z_Cl},
    {MD_STR_INIT("CLA:CL"),MD_Z_Cl},
    {MD_STR_INIT("CLA:CLA"),MD_Z_Cl},
    {MD_STR_INIT("CA:CA"), MD_Z_Ca},
    {MD_STR_INIT("CAL:CA"),MD_Z_Ca},
    {MD_STR_INIT("MG:MG"), MD_Z_Mg},
    {MD_STR_INIT("MAG:MG"),MD_Z_Mg},
    {MD_STR_INIT("ZN:ZN"), MD_Z_Zn},
    {MD_STR_INIT("FE:FE"), MD_Z_Fe},
    {MD_STR_INIT("CU:CU"), MD_Z_Cu},
    {MD_STR_INIT("MN:MN"), MD_Z_Mn},
    {MD_STR_INIT("F: F"),  MD_Z_F},
    {MD_STR_INIT("BR:BR"), MD_Z_Br},
    {MD_STR_INIT("I: I"),  MD_Z_I},
    {MD_STR_INIT("IOD:I"), MD_Z_I},

    // Other
    {MD_STR_INIT("POPC:HS"), MD_Z_H},
    {MD_STR_INIT("HOH:O"), MD_Z_O},
    {MD_STR_INIT("WAT:O"), MD_Z_O},
    {MD_STR_INIT("H2O:O"), MD_Z_O},
    {MD_STR_INIT("DOD:O"), MD_Z_O}, // Heavy water
    {MD_STR_INIT("CO3:C"), MD_Z_C}, // Carbonate
    {MD_STR_INIT("CO3:O"), MD_Z_O},
    {MD_STR_INIT("CO2:C"), MD_Z_C}, // Carbon dioxide
    {MD_STR_INIT("CO2:O"), MD_Z_O},
    {MD_STR_INIT("NH4:N"), MD_Z_N}, // Ammonium
    {MD_STR_INIT("NH4:H"), MD_Z_H},
    {MD_STR_INIT("NO3:N"), MD_Z_N}, // Nitrate
    {MD_STR_INIT("NO3:O"), MD_Z_O},
    {MD_STR_INIT("SO4:S"), MD_Z_S}, // Sulfate
    {MD_STR_INIT("SO4:O"), MD_Z_O},

    // Hem
    {MD_STR_INIT("HEM:FE"), MD_Z_Fe},
};

static const md_str_t resname_fallback_to_single_character[] = {
    MD_STR_INIT("ALA"), MD_STR_INIT("ARG"),
    MD_STR_INIT("ASN"), MD_STR_INIT("ASP"), MD_STR_INIT("CYS"), MD_STR_INIT("GLN"),
    MD_STR_INIT("GLU"), MD_STR_INIT("GLY"), MD_STR_INIT("HIS"), MD_STR_INIT("ILE"), MD_STR_INIT("LEU"), MD_STR_INIT("LYS"), MD_STR_INIT("MET"), MD_STR_INIT("PHE"), MD_STR_INIT("PRO"), MD_STR_INIT("SER"), MD_STR_INIT("THR"), MD_STR_INIT("TRP"), MD_STR_INIT("TYR"), MD_STR_INIT("VAL"),
    MD_STR_INIT("DA"),  MD_STR_INIT("DC"),  MD_STR_INIT("DG"),  MD_STR_INIT("DT"),  MD_STR_INIT("A"),   MD_STR_INIT("C"),   MD_STR_INIT("G"),   MD_STR_INIT("U"),
    MD_STR_INIT("HOH"), MD_STR_INIT("H2O"), MD_STR_INIT("DOD"), MD_STR_INIT("WAT"),
    
    MD_STR_INIT("HEM"), MD_STR_INIT("HEME"),
    MD_STR_INIT("ADP"), MD_STR_INIT("ATP"), MD_STR_INIT("GDP"), MD_STR_INIT("GTP"),

    MD_STR_INIT("POPC"), MD_STR_INIT("DPPC"), MD_STR_INIT("DMPC"), MD_STR_INIT("DOPC"), MD_STR_INIT("POPE"), MD_STR_INIT("DOPE"),

    MD_STR_INIT("CHOL"),  // Cholesterol
};

static inline md_atomic_number_t find_in_mappings(md_str_t key, const Mapping* mappings, size_t num_mappings) {
    for (size_t i = 0; i < num_mappings; ++i) {
        if (md_str_eq_ignore_case(key, mappings[i].str)) {
            return mappings[i].z;
        }
    }
    return MD_Z_X;
}

static inline bool resname_in_fallback_list(md_str_t res_name) {
    for (size_t i = 0; i < MD_COUNTOF(resname_fallback_to_single_character); ++i) {
        if (md_str_eq_ignore_case(res_name, resname_fallback_to_single_character[i])) {
            return true;
        }
    }
    return false;
}

// Inference functions
md_atomic_number_t md_atomic_number_infer_from_label(md_str_t atom_name, md_str_t res_name, size_t res_size) {    

    // Special case: if atom name is empty but residue name is an element (ion case)
    if (atom_name.len == 0 && res_name.len > 0) {
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            return res_element;
        }
        return MD_Z_X;
    }
    
    if (atom_name.len == 0) return MD_Z_X;

    if (md_ascii_is_digit(atom_name.ptr[0])) {
        atom_name = md_str_substr(atom_name, 1, SIZE_MAX); // Skip leading digits
    }
    
    // First try residue+atom combination
    if (res_name.len > 0) {
        // Special case: if water, amino acid or nucleotide: Match against first character only
        if (md_resname_water(res_name) || md_resname_amino_acid(res_name) || md_resname_nucleotide(res_name)) {
            md_str_t first_char = md_str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }

		char combined_key_buf[16];
		int len = md_snprintf(combined_key_buf, sizeof(combined_key_buf), MD_STR_FMT ":" MD_STR_FMT, MD_STR_ARG(res_name), MD_STR_ARG(atom_name));
		md_str_t combined_key = { combined_key_buf, (size_t)len };

		// Check predefined resname : atomlabel mappings first
        {
        md_atomic_number_t z = find_in_mappings(combined_key, predefined_mappings, MD_COUNTOF(predefined_mappings));
            if (z != MD_Z_X) return z;
        }

        // Check if we should fallback to single-character atom name
        if (resname_in_fallback_list(res_name)) {
            md_str_t first_char = md_str_substr(atom_name, 0, 1);
            md_atomic_number_t z = md_atomic_number_from_symbol(first_char, true);
            if (z != MD_Z_X) return z;
        }
        
        // If residue name itself is an element (ions for example)
        md_atomic_number_t res_element = md_atomic_number_from_symbol(res_name, true);
        if (res_element != MD_Z_X) {
            // If atom name is empty or equals residue, return that element
            if (md_str_eq_ignore_case(atom_name, res_name)) {
                return res_element;
            }
        }
    }

    // User has provided a residue size and its likely not an ion
    if (res_size > 2) {
        md_str_t two_letter_str = md_str_substr(atom_name, 0, 2);
        md_atomic_number_t two_z = md_atomic_number_from_symbol(two_letter_str, false);
        if (two_z != MD_Z_X) {
            return two_z;
        }

        // Try single-letter element heuristic (e.g., C1 => C, N2 => N)
        md_str_t first_letter_str = md_str_substr(atom_name, 0, 1);
        md_atomic_number_t first_z = md_atomic_number_from_symbol(first_letter_str, true);
        if (first_z != MD_Z_X) {
            return first_z;
        }
    }
    
    // Try two-letter element heuristic (e.g., CL12 => Cl, BR1 => Br)
    if (md_str_len(atom_name) >= 2) {
        md_str_t two_letter_str = md_str_substr(atom_name, 0, 2);
        md_atomic_number_t two_z = md_atomic_number_from_symbol(two_letter_str, true);
        if (two_z != MD_Z_X) {
            return two_z;
        }
    }
    
    // Final fallback: first-letter element mapping
    md_str_t first_letter_str = md_str_substr(atom_name, 0, 1);
    return md_atomic_number_from_symbol(first_letter_str, true);
}

md_atomic_number_t md_atomic_number_infer_from_mass(float mass) {
    md_atomic_number_t res = MD_Z_X;
    md_atomic_number_infer_from_mass_batch(&res, &mass, 1);
    return res;
}

