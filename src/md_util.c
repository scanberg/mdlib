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

#include <math.h>
#include <string.h>
#include <float.h>

#define STB_DS_IMPLEMENTATION
#include <stb_ds.h>

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
static const float element_covalent_radii[] = {
    0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53,
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

static const char* peptides[] = {"APN", "CPN", "TPN", "GPN"};
static const char* rna[] = {"A", "C", "T", "G", "I", "U", "N"};
static const char* dna[] = {"DA", "DC", "DT", "DG", "DI", "DU", "DN"};

static const char* acidic[] = { "ASP", "GLU" };
static const char* basic[] = { "ARG", "HIS", "LYS" };

static const char* neutral[] = { "VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP" };
static const char* water[] = { "H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4", "W", "DOD", "D30" };
static const char* hydrophobic[] = { "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "CYX" };

static const uint8_t aromatic_ring_elements[] = {
    B, C, N, O, Si, P, S, Ge, As, Sn, Sb, Bi
};

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

static inline md_256 simd_deperiodize(md_256 x, md_256 r, md_256 period, md_256 inv_period) {
    md_256 d = md_mm256_sub_ps(x, r);
    md_256 dx = md_mm256_mul_ps(d, inv_period);
    dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
    return md_mm256_add_ps(r, md_mm256_mul_ps(dx, period));
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
    return find_str_in_array(str, rna, ARRAY_SIZE(rna)) != -1;
}
    
bool md_util_resname_dna(str_t str) {
    str = trim_label(str);
    return find_str_in_array(str, dna, ARRAY_SIZE(dna)) != -1;
}

bool md_util_resname_nucleic_acid(str_t str) {
    str = trim_label(str);
    return find_str_in_array(str, rna, ARRAY_SIZE(rna)) != -1 || find_str_in_array(str, dna, ARRAY_SIZE(dna)) != -1;
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
            if (str_equal_cstr(lbl, "N"))       bits |= BIT_N;
            else if (str_equal_cstr(lbl, "CA")) bits |= BIT_CA;
            else if (str_equal_cstr(lbl, "C"))  bits |= BIT_C;
            else if (str_equal_cstr(lbl, "O"))  bits |= BIT_O;
            count += 1;
        }
    }
    return 5 <= count && count <= 15 && bits == (BIT_N | BIT_CA | BIT_C | BIT_O);
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

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    // @TODO, PERF: Iterate over residues and check if the entire residue is an amino acid
    md_array(uint8_t) res_flags = 0;
    if (mol->residue.count > 0) {
        res_flags = md_array_create(uint8_t, mol->residue.count, temp_alloc);
        MEMSET(res_flags, 0, md_array_bytes(res_flags));
    }

    enum {
        FLAG_AMINO_OR_NUCLEIC = 1,
    };

    for (int64_t i = 0; i < mol->residue.count; ++i) {
        const md_range_t res_range = mol->residue.atom_range[i];
        const int res_len = res_range.end - res_range.beg;

        if (MIN_RES_LEN < res_len && res_len < MAX_RES_LEN && mol->residue.name) {
            str_t resname = LBL_TO_STR(mol->residue.name[i]);
            if (md_util_resname_amino_acid(resname) || md_util_resname_nucleic_acid(resname) ||
                amino_acid_heuristic(mol->atom.type + res_range.beg, res_range.end - res_range.beg))
            {
                res_flags[i] |= FLAG_AMINO_OR_NUCLEIC;
            }
        }
    }

    const int64_t count = MIN(capacity, mol->atom.count);
    for (int64_t i = 0; i < count; ++i) {
        if (element[i] != 0) continue;

        str_t original = LBL_TO_STR(mol->atom.type[i]);

        // Trim whitespace, digits and 'X's
        str_t name = trim_label(original);

        if (name.len > 0) {

            md_element_t elem = 0;
            if ((elem = md_util_element_lookup(name)) != 0) goto done;

            // If amino acid, try to deduce the element from that
            if (mol->atom.res_idx) {
                const int64_t res_idx = mol->atom.res_idx[i];
                if (res_flags[res_idx] & FLAG_AMINO_OR_NUCLEIC) {
                    // EASY-PEASY, we just try to match against the first character
                    name.len = 1;
                    elem = md_util_element_lookup_ignore_case(name);
                    goto done;
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
        MD_LOG_ERROR("mass is null");
        return false;
    }

    int64_t failed_matches = 0;

    const float eps = 1.0e-2f;
    for (int64_t i = 0; i < count; ++i) {
        md_element_t elem = 0;
        const float m = mass[i];
        //Loop through the mass options, break when match is found
        for (uint8_t j = 1; j < (uint8_t)ARRAY_SIZE(element_atomic_mass); ++j) {
            if (fabs(m - element_atomic_mass[j]) < eps) {
                elem = j;
                break;
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
        MD_LOG_ERROR("%i masses had no matching element", (int)failed_matches);
        return false;
    }
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
    return extract_backbone_atoms(backbone_atoms, mol->atom.type, mol->residue.atom_range[res_idx]);
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

static bool md_compute_connectivity(md_conn_data_t* conn, md_atom_data_t* atom, const md_bond_data_t* bond, md_allocator_i* alloc) {
    ASSERT(conn);
    ASSERT(atom);
    ASSERT(bond);
    ASSERT(alloc);

    md_array_resize(atom->conn_offset, atom->count + 1, alloc);
    MEMSET(atom->conn_offset, 0, md_array_bytes(atom->conn_offset));

    // This have length of 2 * bond_count (one for each direction of the bond)
    conn->count = 2 * bond->count;
    md_array_resize(conn->index, conn->count, alloc);

    if (bond->order) {
        md_array_resize(conn->order, conn->count, alloc);
    }

    if (bond->flags) {
        md_array_resize(conn->flags, conn->count, alloc);
    }

    md_array(uint32_t) local_offset = md_array_create(uint32_t, bond->count, md_heap_allocator);

    // Two packed 16-bit local offsets for each of the bond idx
    // Use offsets as accumulators for length
    for (int64_t i = 0; i < bond->count; ++i) {
        const uint32_t off0 = atom->conn_offset[bond->pairs[i].idx[0]]++;
        const uint32_t off1 = atom->conn_offset[bond->pairs[i].idx[1]]++;
        local_offset[i] = (off1 << 16) | off0;
    }

    // Compute complete edge offsets (exclusive scan)
    uint32_t off = 0;
    for (int64_t i = 0; i < atom->count + 1; ++i) {
        const uint32_t len = atom->conn_offset[i];
        atom->conn_offset[i] = off;
        off += len;
    }

    // Write edge indices to correct location
    for (int64_t i = 0; i < bond->count; ++i) {
        const md_bond_pair_t p = bond->pairs[i];
        const int atom_a = p.idx[0];
        const int atom_b = p.idx[1];
        const int local_a = (int)(local_offset[i] & 0xFFFF);
        const int local_b = (int)(local_offset[i] >> 16);
        const int off_a = atom->conn_offset[atom_a];
        const int off_b = atom->conn_offset[atom_b];

        const int idx_a = off_a + local_a;
        const int idx_b = off_b + local_b;

        ASSERT(idx_a < conn->count);
        ASSERT(idx_b < conn->count);

        // Store the cross references to the 'other' atom index signified by the bond in the correct location
        conn->index[idx_a] = atom_b;
        conn->index[idx_b] = atom_a;

        if (bond->order) {
            conn->order[idx_a] = bond->order[i];
            conn->order[idx_b] = bond->order[i];
        }

        if (bond->flags) {
            conn->flags[idx_a] = bond->flags[i];
            conn->flags[idx_b] = bond->flags[i];
        }
    }

    md_array_free(local_offset, md_heap_allocator);
    return true;
}

// Taken from here
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/model/properties/atomic/bonds.ts
/*
* Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
* @author Alexander Rose <alexander.rose@weirdbyte.de>
*/

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
bool md_util_compute_covalent_bond_order(md_order_t* bond_order, const md_bond_pair_t* bond_pairs, int64_t bond_count, const md_label_t* type, const md_label_t* resname) {
    ASSERT(bond_order);

    if (!bond_pairs || bond_count == 0) {
        return false;
    }

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
        if (str_equal(comp_a, comp_b)) {
            if (str_compare_lex(type_a, type_b) > 0) {
                str_swap(type_a, type_b);
            }
            // Intra
            if (md_util_resname_amino_acid(comp_a) && str_equal_cstr(type_a, "C") && str_equal_cstr(type_b, "O")) {
                order = 2;
            } else {
                char buf[32];
                int len = build_key(buf, comp_a, type_a, type_b);
                if (find_str_in_array((str_t){buf,len}, intra_bond_order_table, ARRAY_SIZE(intra_bond_order_table)) >= 0) {
                    bond_order[i] = 2;
                }
            }
        } else {
            // Inter
            if ( (str_equal_cstr(comp_a, "LYS") && str_equal_cstr(type_a, "CZ") && str_equal_cstr(comp_b, "RET") && str_equal_cstr(type_b, "C15")) ||
                 (str_equal_cstr(comp_b, "LYS") && str_equal_cstr(type_b, "CZ") && str_equal_cstr(comp_a, "RET") && str_equal_cstr(type_a, "C15")) ){
                bond_order[i] = 2;
            }
        }
        bond_order[i] = order;
    }
    
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

    if (a == C && b == C) {
        while(0);
    }

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

static inline uint8_t covalent_bond_heuristic(float d2, md_element_t elem_a, md_element_t elem_b) {
    const float ra = element_covalent_radii[elem_a];
    const float rb = element_covalent_radii[elem_b];
    const float r_sum = ra + rb;
    const float r_min = r_sum - 0.5f;
    const float r_max = r_sum + 0.3f;

    return (r_min * r_min) < d2 && d2 < (r_max * r_max) ? 1 : 0;
}

int64_t md_util_compute_covalent_bounds_upper_bound(const md_element_t* element, int64_t count) {
    int64_t result = 0;
    for (int64_t i = 0; i < count; ++i) {
        result += md_util_element_max_valence(element[i]);
    }
    return result;
}

void find_intra_bonds_in_range(md_bond_data_t* bond, const md_atom_data_t* atom, const md_unit_cell_t* cell, md_range_t range, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    ASSERT(bond);
    ASSERT(atom);
    ASSERT(cell);

    const vec4_t pbc_ext = cell ? vec4_from_vec3(mat3_diag(cell->basis), 0) : vec4_zero();

    if (range.end - range.beg < 100) {
        // Brute force
        for (int i = range.beg; i < range.end - 1; ++i) {
            for (int j = i + 1; j < range.end; ++j) {
                const vec4_t a = {atom->x[i], atom->y[i], atom->z[i], 0};
                const vec4_t b = {atom->x[j], atom->y[j], atom->z[j], 0};
                const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext);
                const uint8_t order = covalent_bond_heuristic(d2, atom->element[i], atom->element[j]);
                if (order) {
                    md_array_push(bond->pairs, ((md_bond_pair_t){i, j}), alloc);
                    md_array_push(bond->order, order, alloc);
                    md_array_push(bond->flags, 0, alloc);
                    bond->count += 1;
                }
            }
        }
    }
    else {
		// Spatial acceleration structure
        md_spatial_hash_t* sh = md_spatial_hash_create_soa(atom->x, atom->y, atom->z, NULL, atom->count, cell, temp_alloc);

        const float cutoff = 3.0f;
        int indices[128];

        for (int i = 0; i < (int)atom->count; ++i) {
            const vec3_t pos = {atom->x[i], atom->y[i], atom->z[i]};

            const int64_t num_indices = md_spatial_hash_query_idx(indices, ARRAY_SIZE(indices), sh, pos, cutoff);
            for (int64_t iter = 0; iter < num_indices; ++iter) {
                const int j = indices[iter];
                // Only store monotonic bond connections
                if (j < i) continue;

                const vec4_t pos_a = vec4_from_vec3(pos, 0);
                const vec4_t pos_b = {atom->x[j], atom->y[j], atom->z[j], 0};
                const float d2 = vec4_periodic_distance_squared(pos_a, pos_b, pbc_ext);
                const uint8_t order = covalent_bond_heuristic(d2, atom->element[i], atom->element[j]);
                if (order) {
                    md_array_push(bond->pairs, ((md_bond_pair_t){i, j}), alloc);
                    md_array_push(bond->order, order, alloc);
                    md_array_push(bond->flags, 0, alloc);
                    bond->count += 1;
                }
            }
        }
	}
}

void find_inter_bonds_in_ranges(md_bond_data_t* bond, const md_atom_data_t* atom, const md_unit_cell_t* cell, md_range_t range_a, md_range_t range_b, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    ASSERT(bond);
    ASSERT(atom);
    ASSERT(cell);

    const vec4_t pbc_ext = cell ? vec4_from_vec3(mat3_diag(cell->basis), 0) : vec4_zero();

    const int N = (range_a.end - range_a.beg) * (range_b.end - range_b.beg);
    if (N < 1000) {
        // Brute force
        for (int i = range_a.beg; i < range_a.end; ++i) {
            for (int j = range_b.beg; j < range_b.end; ++j) {
                const vec4_t a = {atom->x[i], atom->y[i], atom->z[i], 0};
                const vec4_t b = {atom->x[j], atom->y[j], atom->z[j], 0};
                const float d2 = vec4_periodic_distance_squared(a, b, pbc_ext);
                const uint8_t order = covalent_bond_heuristic(d2, atom->element[i], atom->element[j]);
                if (order) {
                    md_array_push(bond->pairs, ((md_bond_pair_t){i, j}), alloc);
                    md_array_push(bond->order, order, alloc);
                    md_array_push(bond->flags, 0, alloc);
                    bond->count += 1;
                }
            }
        }
    }
    else {
        // Spatial acceleration structure
        int count_b = range_b.end - range_a.beg;
        md_spatial_hash_t* sh = md_spatial_hash_create_soa(atom->x + range_b.beg, atom->y + range_b.beg, atom->z + range_b.beg, NULL, count_b, cell, temp_alloc);

        const float cutoff = 3.0f;
        int indices[128];

        for (int i = range_a.beg; i < range_a.end; ++i) {
            const vec3_t pos = {atom->x[i], atom->y[i], atom->z[i]};

            const int64_t num_indices = md_spatial_hash_query_idx(indices, ARRAY_SIZE(indices), sh, pos, cutoff);
            for (int64_t iter = 0; iter < num_indices; ++iter) {
                const int j = indices[iter];
                // Only store monotonic bond connections
                if (j < i) continue;

                const vec4_t pos_a = vec4_from_vec3(pos, 0);
                const vec4_t pos_b = {atom->x[j], atom->y[j], atom->z[j], 0};
                const float d2 = vec4_periodic_distance_squared(pos_a, pos_b, pbc_ext);
                const uint8_t order = covalent_bond_heuristic(d2, atom->element[i], atom->element[j]);
                if (order) {
                    md_array_push(bond->pairs, ((md_bond_pair_t){i, j}), alloc);
                    md_array_push(bond->order, order, alloc);
                    md_array_push(bond->flags, 0, alloc);
                    bond->count += 1;
                }
            }
        }

        md_spatial_hash_free(sh);
    }
}

md_bond_data_t md_util_compute_covalent_bonds(const md_atom_data_t* atom, const md_unit_cell_t* cell, md_allocator_i* alloc) {
    ASSERT(atom);
    ASSERT(alloc);

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    md_bond_data_t bond = {0};

    if (atom->count < 0) {
        MD_LOG_ERROR("Incorrect number of atoms: %i", (int)atom->count);
        goto done;
    }
    
    if (!atom->x || !atom->y || !atom->z) {
        MD_LOG_ERROR("Missing atom field (x/y/z)");
        goto done;
    }

    if (!atom->element) {
        MD_LOG_ERROR("Missing atom field element");
        goto done;
    }

    if (cell->flags & MD_CELL_TRICLINIC) {
        MD_LOG_ERROR("Triclinic cells are not supported yet! Sorry!");
        goto done;
    }
    const int atom_count = (int)atom->count;
       
    if (atom->resid) {
        md_range_t range = {0, 0};
        md_range_t prev_range = {0, 0};
        for (int i = 0; i < atom->count; ++i) {
            while (i < atom_count && atom->resid[i] == atom->resid[range.beg]) ++i;
            range.end = i;

			find_inter_bonds_in_ranges(&bond, atom, cell, prev_range, range, alloc, temp_alloc);
            find_intra_bonds_in_range(&bond, atom, cell, range, alloc, temp_alloc);

            prev_range = range;
            range.beg  = range.end;
        }
    }
    else {
        find_intra_bonds_in_range(&bond, atom, cell, (md_range_t){0, atom_count}, alloc, temp_alloc);
    }

    md_util_compute_covalent_bond_order(bond.order, bond.pairs, bond.count, atom->type, atom->resname);

done:
    md_arena_allocator_destroy(temp_alloc);
    return bond;
}

// Compute the prerequisite fields to enable hydrogen bond determination
// ported from molstar
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

bool md_util_compute_residue_data(md_residue_data_t* res_data, md_atom_data_t* atom, md_allocator_i* alloc) {
    ASSERT(res_data);
    ASSERT(atom);
    ASSERT(alloc);

    md_residue_id_t prev_resid = -1;
    md_flags_t prev_flags = 0;
    for (int64_t i = 0; i < atom->count; ++i) {
        const md_residue_id_t resid = atom->resid[i];
        const md_flags_t flags = atom->flags[i];
        const md_label_t resname = atom->resname[i];        

        if (resid != prev_resid || flags & MD_FLAG_RES_BEG || prev_flags & MD_FLAG_RES_END) {
            md_range_t range = {(int)i, (int)i};
            md_array_push(res_data->id, resid, alloc);
            md_array_push(res_data->name, resname, alloc);
            md_array_push(res_data->atom_range, range, alloc);
            md_array_push(res_data->flags, 0, alloc);
            prev_resid = resid;
            res_data->count += 1;
        }

		md_array_last(res_data->atom_range)->end += 1;
        prev_flags = flags;
    }

    atom->res_idx = md_array_create(md_residue_idx_t, atom->count, alloc);

    for (int64_t i = 0; i < res_data->count; ++i) {
        str_t resname = LBL_TO_STR(res_data->name[i]);
        md_range_t range = res_data->atom_range[i];
        if (md_util_resname_amino_acid(resname) || amino_acid_heuristic(atom->type + range.beg, range.end - range.beg)) {
			res_data->flags[i] |= MD_FLAG_AMINO_ACID;
		} else if (md_util_resname_nucleic_acid(resname)) {
            res_data->flags[i] |= MD_FLAG_NUCLEIC_ACID;
        } else if (md_util_resname_water(resname)) {
            res_data->flags[i] |= MD_FLAG_WATER;
        }

        for (int32_t j = range.beg; j < range.end; ++j) {
            atom->res_idx[j] = (md_residue_idx_t)i;
        }
    }

    return true;
}

bool md_util_compute_chain_data(md_chain_data_t* chain_data, md_atom_data_t* atom, const md_residue_data_t* res, const md_bond_data_t* bond, md_allocator_i* alloc) {
    ASSERT(chain_data);
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

    md_allocator_i* temp_alloc = md_temp_allocator;
    
    md_array(bool) res_bond_to_next = md_array_create(bool, res->count, temp_alloc);
    MEMSET(res_bond_to_next, 0, md_array_bytes(res_bond_to_next));

    // Iterate over bonds and determine bonds between residues (res_bonds)
    for (int64_t i = 0; i < bond->count; ++i) {
        const md_bond_pair_t pair = bond->pairs[i];
        const md_residue_idx_t res_a = atom->res_idx[pair.idx[0]];
        const md_residue_idx_t res_b = atom->res_idx[pair.idx[1]];
        
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
    for (int i = 0; i < res->count; ++i) {
        if (res_bond_to_next[i] == false || i == res->count - 1) {
            int end_idx = i + 1;
            if (end_idx - beg_idx > 1) {
                const md_label_t id = {.buf = {'A' + next_char}, .len = 1};
                const md_range_t atom_range = {res->atom_range[beg_idx].beg, res->atom_range[end_idx - 1].end};
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

    atom->chain_idx = md_array_create(md_chain_idx_t, atom->count, alloc);

    for (int64_t i = 0; i < chain_data->count; ++i) {
        const md_range_t atom_range = chain_data->atom_range[i];
        for (int64_t j = atom_range.beg; j < atom_range.end; ++j) {
            atom->chain_idx[j] = (md_chain_idx_t)i;
        }
    }

    return chain_data;
}

bool md_util_compute_atom_valence(md_valence_t atom_valence[], int64_t atom_count, const md_bond_pair_t bond_pairs[], int64_t bond_count) {
    if (!atom_valence) {
        MD_LOG_ERROR("Missing input: atom valence");
        return false;
    }

    if (atom_count < 0) {
        MD_LOG_ERROR("Invalid input: atom count");
        return false;
    }

    if (!bond_pairs) {
        MD_LOG_ERROR("Missing input: bond pairs");
        return false;
    }

    if (bond_count < 0) {
        MD_LOG_ERROR("Invalid input: bond count");
        return false;
    }

    MEMSET(atom_valence, 0, sizeof(md_valence_t) * atom_count);

    for (int64_t i = 0; i < bond_count; ++i) {
        const md_bond_pair_t bond = bond_pairs[i];
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
    MEMSET(fifo->data, 0, md_array_bytes(fifo->data));
#endif
}

static inline fifo_t fifo_create(int64_t capacity, md_allocator_i* alloc) {
    fifo_t fifo;
    fifo_init(&fifo, capacity, alloc);
    return fifo;
}

static void fifo_free(fifo_t* fifo) {
    ASSERT(fifo);
    md_free(fifo->alloc, fifo->data, sizeof(int64_t) * fifo->cap);
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

#define MIN_RING_SIZE 3
#define MAX_RING_SIZE 6
#define MAX_DEPTH (MAX_RING_SIZE/2 + 1)

// Inspired by molstars implementation
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/structure/unit/rings/compute.ts#L249
//
// Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
// @author David Sehnal <david.sehnal@gmail.com>

md_index_data_t md_util_compute_rings(const md_molecule_t* mol, md_allocator_i* alloc) {
    ASSERT(alloc);

    md_index_data_t ring_data = {0};

    if (mol->atom.count == 0) {
        return ring_data;
    }

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
        
    int* color = md_alloc(temp_alloc, mol->atom.count * sizeof(int));
    int* depth = md_alloc(temp_alloc, mol->atom.count * sizeof(int));
    int* mark  = md_alloc(temp_alloc, mol->atom.count * sizeof(int));
    int* pred  = md_alloc(temp_alloc, mol->atom.count * sizeof(int));

    MEMSET(color, 0, mol->atom.count * sizeof(int));
#if DEBUG
    MEMSET(depth, 0, mol->atom.count * sizeof(int));
    MEMSET(mark,  0, mol->atom.count * sizeof(int));
    MEMSET(pred, -1, mol->atom.count * sizeof(int));    // We can do memset as the representation of -1 under two's complement is 0xFFFFFFFF
#endif

    // The capacity is arbitrary here, but will be resized if needed.
    fifo_t queue = fifo_create(64, temp_alloc);

    // Hashmap element struct
    typedef struct T {
        uint64_t key;
    } T;
    T *hm = NULL;
    
    const size_t seed = 12;
    int current_color = 1;
    int current_mark  = 1;

#ifdef DEBUG
    uint64_t processed_ring_elements = 0;
#endif

    for (int i = 0; i < mol->atom.count; ++i) {
        // Skip any atom that has already been colored in the previous search
        if (color[i] == current_color) continue;

        current_mark++;

        // Set i as our root
        depth[i] = 1;
        pred[i]  = -1;
        mark[i] = current_mark;

        fifo_clear(&queue);
        fifo_push(&queue, i);
        while (!fifo_empty(&queue)) {
            int idx = fifo_pop(&queue);

#ifdef DEBUG
            processed_ring_elements += 1;
#endif

            md_conn_iter_t it = md_conn_iter(mol, idx);
            while (md_conn_iter_has_next(&it)) {
                int next = md_conn_iter_index(&it);
                md_conn_iter_next(&it);
                if (next == pred[idx]) continue;  // avoid adding parent to search queue

                if (mark[next] == current_mark) {
                    // We found a junction point where the graph connects
                    // Now we need to traverse both branches of this graph back up
                    // In order to find the other junction where the graph branches
                    // We follow both branches up towards the root.

                    int l = idx;
                    int r = next;

                    // Only process one of the two branches/cases
                    if (r < l) {
                        int len = 0;
                        int ring[MAX_DEPTH * 2];

                        int col = current_color++;
                        int cur;

                        cur = l;
                        for (int d = 0; d < MAX_DEPTH; ++d) {
							color[cur] = col;
							cur = pred[cur];
                            if (cur < 0) break;
                        }

                        bool found = false;
                        int target = 0;
                        cur = r;
                        for (int d = 0; d < MAX_DEPTH; ++d) {
                            if (color[cur] == col) {
                                target = cur;
                                found = true;
								break;
							}
                            ring[len++] = cur;
                            cur = pred[cur];
                            if (cur < 0) break;
                        }

                        if (found) {
                            cur = l;
                            for (int d = 0; d < MAX_DEPTH; ++d) {
                                ring[len++] = cur;
                                if (target == cur) break;
                                cur = pred[cur];
                                if (cur < 0) break;
                            }

                            if (MIN_RING_SIZE <= len && len <= MAX_RING_SIZE) {
                                sort_arr(ring, len);
                                size_t key = stbds_hash_bytes(ring, len * sizeof(int), seed);
                                if (stbds_hmgeti(hm, key) == -1) {
                                    stbds_hmputs(hm, (T){key});
                                    md_index_data_push(&ring_data, ring, len, alloc);
                                    goto next;
                                }
                            }
                        }
                    }
                } else {
                    // Avoid expanding too far from the root as we are only interested in rings up to a certain size
                    // Otherwise, we will may have alot of potential large cycles such as in the case of C60 or graphene
                    int d = depth[idx] + 1;
                    if (d > MAX_DEPTH) continue;

                    depth[next] = d;
                    pred[next]  = idx;
                    mark[next]  = current_mark;
                    fifo_push(&queue, next);
                }
            }
        }
        next:;
    }
    stbds_hmfree(hm);
    
    md_arena_allocator_destroy(temp_alloc);

#if 0
    MD_LOG_DEBUG("Processed ring elements: %llu\n", processed_ring_elements);
#endif

    return ring_data;
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE
#undef MAX_DEPTH

static md_array(uint64_t) copy_bitfield(const md_array(uint64_t) src, md_allocator_i* alloc) {
	md_array(uint64_t) copy = md_array_create(uint64_t, md_array_size(src), alloc);
	MEMCPY(copy, src, md_array_bytes(src));
	return copy;
}

static void set_bitfield(md_array(uint64_t) dst, const md_array(uint64_t) src) {
	MEMCPY(dst, src, md_array_bytes(src));
}

static void clear_all_bitfield(md_array(uint64_t) bits) {
	MEMSET(bits, 0, md_array_bytes(bits));
}

static void set_all_bitfield(md_array(uint64_t) bits, int64_t num_bits) {
    MEMSET(bits, 0xffffffff, md_array_bytes(bits));
    uint64_t mask = 1ULL << (num_bits & 63);
    *md_array_last(bits) &= mask - 1;
}

static bool find_first_bit_set_bitfield(int* idx, const md_array(uint64_t) bits) {
    for (int64_t i = 0; i < md_array_size(bits); ++i) {
		if (bits[i]) {
			*idx = (int)ctz64(bits[i]);
            return true;
		}
	}
    return false;
}

static md_array(uint64_t) make_bitfield(uint64_t num_bits, md_allocator_i* alloc) {
	md_array(uint64_t) bits = md_array_create(uint64_t, DIV_UP(num_bits, 64), alloc);
	MEMSET(bits, 0, md_array_bytes(bits));
	return bits;
}

static inline bool test_bit(const uint64_t* bits, int64_t idx) {
	return (bits[idx >> 6] & (1ULL << (idx & 63)));
}

static inline void set_bit(uint64_t* bits, int64_t idx) {
    bits[idx >> 6] |= (1ULL << (idx & 63));
}

static inline void clear_bit(uint64_t* bits, int64_t idx) {
    bits[idx >> 6] &= ~(1ULL << (idx & 63));
}

static inline uint64_t popcount_bits(const uint64_t* beg, const uint64_t* end) {
    const uint64_t* it = beg;
    uint64_t count = 0;
    while (it != end) {
        count += popcnt64(*it);
        it++;
    }
    return count;
}

static inline uint64_t popcount_bitfield(const md_array(uint64_t) bits) {
    return popcount_bits(bits, bits + md_array_size(bits));
}

md_index_data_t md_util_compute_structures(const md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    md_index_data_t structures = {0};

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    // Create a bitfield to keep track of which atoms have been visited
    const uint64_t bytes = DIV_UP(mol->atom.count, 64) * sizeof(uint64_t);
    uint64_t* visited = md_alloc(temp_alloc, bytes);
    MEMSET(visited, 0, bytes);

    // The capacity is arbitrary here, and will be resized if needed.
    fifo_t queue = fifo_create(128, temp_alloc);

    md_array(int) indices = 0;
    md_array_ensure(indices, 256, temp_alloc);

    for (int i = 0; i < mol->atom.count; ++i) {
        // Skip any atom which has already been touched
        if (test_bit(visited, i)) continue;

        fifo_clear(&queue);
        fifo_push(&queue, i);
        while (!fifo_empty(&queue)) {
            int cur = fifo_pop(&queue);
            if (test_bit(visited, cur)) continue;
            set_bit(visited, cur);

            md_array_push(indices, cur, temp_alloc);
            
            md_conn_iter_t it = md_conn_iter(mol, cur);
            while (md_conn_iter_has_next(&it)) {
				int next = md_conn_iter_index(&it);
				if (!test_bit(visited, next)) {
					fifo_push(&queue, next);
				}
                md_conn_iter_next(&it);
			}
        }

        // Sort the indices within the structure for more coherent memory access
        sort_arr(indices, (int)md_array_size(indices));
        
        // Here we should have exhausted every atom that is connected to index i.
        md_index_data_push(&structures, indices, md_array_size(indices), alloc);
        md_array_shrink(indices, 0);
    }
    
    md_arena_allocator_destroy(temp_alloc);

    return structures;
}

void md_util_grow_mask_by_bonds(md_bitfield_t* mask, const md_molecule_t* mol, int extent, const md_bitfield_t* viable_mask) {
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

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    int* indices = md_alloc(temp_alloc, mask_size * sizeof(int));
    const int64_t num_indices = md_bitfield_extract_indices(indices, mask_size, mask);
    ASSERT(num_indices == mask_size);

    fifo_t queue = fifo_create(64, temp_alloc);

    int64_t num_atoms = mol->atom.count;
    uint8_t* depth = md_alloc(temp_alloc, num_atoms * sizeof(uint8_t));
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

            md_conn_iter_t it = md_conn_iter(mol, idx);
            while (md_conn_iter_has_next(&it)) {
                int next = md_conn_iter_index(&it);
                md_conn_iter_next(&it);
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

// Computes the minimum axis aligned bounding box for a set of points with a given radius (radius is optional), indices are used to select a subset of points
void md_util_compute_aabb_vec4(vec3_t* out_aabb_min, vec3_t* out_aabb_max, const vec4_t* in_xyzr, const int32_t* in_idx, int64_t count) {
    vec4_t aabb_min = vec4_set1( FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);
    if (in_idx) {
        for (int64_t i = 0; i < count; ++i) {
            int32_t idx = in_idx[i];
            vec4_t xyzr = in_xyzr[idx];
            vec4_t r = vec4_splat_w(xyzr);
            aabb_min = vec4_min(aabb_min, vec4_sub(xyzr, r));
            aabb_max = vec4_max(aabb_max, vec4_add(xyzr, r));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            vec4_t xyzr = in_xyzr[i];
            vec4_t r = vec4_splat_w(xyzr);
            aabb_min = vec4_min(aabb_min, vec4_sub(xyzr, r));
            aabb_max = vec4_max(aabb_max, vec4_add(xyzr, r));
        }
    }
    *out_aabb_min = vec3_from_vec4(aabb_min);
    *out_aabb_max = vec3_from_vec4(aabb_max);
}

void md_util_compute_aabb(vec3_t* out_aabb_min, vec3_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int32_t* in_idx, int64_t count) {
    md_256 vx_min = md_mm256_set1_ps(+FLT_MAX);
    md_256 vy_min = md_mm256_set1_ps(+FLT_MAX);
    md_256 vz_min = md_mm256_set1_ps(+FLT_MAX);

    md_256 vx_max = md_mm256_set1_ps(-FLT_MAX);
    md_256 vy_max = md_mm256_set1_ps(-FLT_MAX);
    md_256 vz_max = md_mm256_set1_ps(-FLT_MAX);

    int64_t i = 0;
    const int64_t  simd_elem = 8;
    const int64_t simd_count = ROUND_DOWN(count, simd_elem);

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

    *out_aabb_min = vec3_from_vec4(aabb_min);
    *out_aabb_max = vec3_from_vec4(aabb_max);
}

static vec3_t compute_com_periodic_trig_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, int64_t count, vec3_t ext_max) {
    const vec4_t scl = vec4_div(vec4_set1(TWO_PI), vec4_from_vec3(ext_max, TWO_PI));
    vec4_t acc_c = {0};
    vec4_t acc_s = {0};
    vec4_t acc_xyzw = {0};

    if (in_idx) {
        for (int64_t i = 0; i < count; ++i) {
            int64_t idx  = in_idx[i];
            vec4_t xyzw  = in_xyzw[idx];
            vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(1,0,0,0));
            vec4_t theta = vec4_mul(xyzw, scl);
            vec4_t s,c;
            vec4_sincos(theta, &s, &c);
            acc_s = vec4_add(acc_s, vec4_mul(s, www1));
            acc_c = vec4_add(acc_c, vec4_mul(c, www1));
            acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
        }
    } else {
        for (int64_t i = 0; i < count; ++i) {
            vec4_t xyzw  = in_xyzw[i];
            vec4_t www1  = vec4_blend_mask(vec4_splat_w(xyzw), vec4_set1(1.0f), MD_SIMD_BLEND_MASK(1,0,0,0));
            vec4_t theta = vec4_mul(xyzw, scl);
            vec4_t s,c;
            vec4_sincos(theta, &s, &c);
            acc_s = vec4_add(acc_s, vec4_mul(s, www1));
            acc_c = vec4_add(acc_c, vec4_mul(c, www1));
            acc_xyzw = vec4_add(acc_xyzw, vec4_mul(xyzw, www1));
        }
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

#define TRIG_ATAN2_R2_THRESHOLD 1.0e-8

static vec3_t compute_com_periodic_trig_xyz(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, int64_t count, vec3_t xyz_max) {
    double acc_cx = 0;
    double acc_sx = 0;
    double acc_cy = 0;
    double acc_sy = 0;
    double acc_cz = 0;
    double acc_sz = 0;
    double acc_w = 0;

    int64_t i = 0;

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

    const int64_t simd_count = ROUND_DOWN(count, 16);
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

    const int64_t simd_count = ROUND_DOWN(count, 8);
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

    const int64_t simd_count = ROUND_DOWN(count, 4);
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

static float compute_com_periodic_trig(const float* in_x, const float* in_w, const int32_t* in_idx, int64_t count, float x_max) {
    double acc_c = 0;
    double acc_s = 0;
    double acc_w = 0;

    const double scl = TWO_PI / x_max;
    int64_t i = 0;

#if defined(__AVX512F__)
    const __m512 v_scl = _mm512_set1_ps((float)(TWO_PI / x_max));
    __m512 v_acc_c = _mm512_setzero_ps();
    __m512 v_acc_s = _mm512_setzero_ps();
    __m512 v_acc_w = _mm512_setzero_ps();

    const int64_t simd_count = ROUND_DOWN(count, 16);
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

    const int64_t simd_count = ROUND_DOWN(count, 8);
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

    const int64_t simd_count = ROUND_DOWN(count, 4);
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

static float compute_com_periodic_reg(const float* in_x, const float* in_w, const int32_t* in_idx, int64_t count, float x_max) {
    if (count <= 0) {
        return 0;
    }

    double acc_x = 0;
    double acc_w = 0;

    int64_t i = 0;
#if defined(__AVX512F__)
    const int64_t simd_count = ROUND_DOWN(count, 16);
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
    const int64_t simd_count = ROUND_DOWN(count, 8);
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
    const int64_t simd_count = ROUND_DOWN(count, 4);
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
        int64_t idx = in_idx ? in_idx[i] : i;
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
            double dx = deperiodize(x, r, d_x_max);
            acc_x += dx * w;
            acc_w += w;
        } else {
            int64_t idx = in_idx[i];
            double r = acc_x / acc_w;
            double x = in_x[idx];
            double dx = deperiodize(x, r, d_x_max);
            acc_x += dx;
            acc_w += 1.0;
        }
    } else {
        if (in_w) {
            double r = acc_x / acc_w;
            double x = in_x[i];
            double w = in_w[i];
            double dx = deperiodize(x, r, d_x_max);
            acc_x += dx * w;
            acc_w += w;
        } else {
            double r = acc_x / acc_w;
            double x = in_x[i];
            double dx = deperiodize(x, r, d_x_max);
            acc_x += dx;
            acc_w += 1.0;
        }
    }

    return (float)(acc_x / acc_w);
}

static float compute_com(const float* in_x, const float* in_w, const int32_t* in_idx, int64_t count) {
    ASSERT(in_x);

    if (count <= 0)
        return 0.0f;

    int64_t i = 0;
    double acc_x = 0;
    double acc_w = 0;

#if defined(__AVX512F__)
    __m512 v_acc_x = _mm512_setzero_ps();
    __m512 v_acc_w = _mm512_setzero_ps();
    const int64_t simd_count = ROUND_DOWN(count, 16);
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
    const int64_t simd_count = ROUND_DOWN(count, 8);
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
    const int64_t simd_count = ROUND_DOWN(count, 4);
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

vec3_t md_util_compute_com(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, int64_t count) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    if (count <= 0) {
        return (vec3_t) {0,0,0};
    }

    double acc_x = 0;
    double acc_y = 0;
    double acc_z = 0;
    double acc_w = 0;

    int64_t i = 0;
#if defined(__AVX512F__)
    __m512 vx = _mm512_setzero_ps();
    __m512 vy = _mm512_setzero_ps();
    __m512 vz = _mm512_setzero_ps();
    __m512 vw = _mm512_setzero_ps();

    const int64_t simd_count = ROUND_DOWN(count, 16);

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

    const int64_t simd_count = ROUND_DOWN(count, 8);

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

    const int64_t simd_count = ROUND_DOWN(count, 4);

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
        }
    }
    if (!in_w) acc_w = (double)count;
    double scl = 1.0 / acc_w;
    return (vec3_t) {(float)(acc_x * scl), (float)(acc_y * scl), (float)(acc_z * scl)};
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

vec3_t md_util_compute_com_ortho(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, int64_t count, vec3_t pbc_ext) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);

    if (count <= 0) {
        return (vec3_t) {0,0,0};
    }

#if MD_UTIL_COMPUTE_COM_USE_TRIG
    int mask = (pbc_ext.x > 0 ? 1 : 0) | (pbc_ext.y > 0 ? 2 : 0) | (pbc_ext.z > 0 ? 4 : 0);
    if (mask == 7) {
		return compute_com_periodic_trig_xyz(in_x, in_y, in_z, in_w, in_idx, count, pbc_ext);
	} else if (mask == 0) {
        return md_util_compute_com(in_x, in_y, in_z, in_w, in_idx, count);
    }

    // We need to pick the version based on each component of pdc_ext, since one or more may be zero
    float x = (mask & 1) ? compute_com_periodic_trig(in_x, in_w, in_idx, count, pbc_ext.x) : compute_com(in_x, in_w, in_idx, count);
    float y = (mask & 2) ? compute_com_periodic_trig(in_y, in_w, in_idx, count, pbc_ext.y) : compute_com(in_y, in_w, in_idx, count);
    float z = (mask & 3) ? compute_com_periodic_trig(in_z, in_w, in_idx, count, pbc_ext.z) : compute_com(in_z, in_w, in_idx, count);
#else
    // We need to pick the version based on each component of pdc_ext, since one or more may be zero
    float x = (mask & 1) ? compute_com_periodic_reg(in_x, in_w, in_idx, count, pbc_ext.x) : compute_com(in_x, in_w, in_idx, count);
    float y = (mask & 2) ? compute_com_periodic_reg(in_y, in_w, in_idx, count, pbc_ext.y) : compute_com(in_y, in_w, in_idx, count);
    float z = (mask & 3) ? compute_com_periodic_reg(in_z, in_w, in_idx, count, pbc_ext.z) : compute_com(in_z, in_w, in_idx, count);
#endif

    return (vec3_t) {x, y, z};
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
    return md_util_unit_cell_from_matrix(matrix);
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

double md_util_compute_rmsd_vec4(const vec4_t* xyzw[2], const vec3_t com[2], int64_t count) {
    const mat3_t R = mat3_optimal_rotation_vec4(xyzw[0], xyzw[1], com[0], com[1], count);
    double d_sum = 0;
    double w_sum = 0;
    for (int64_t i = 0; i < count; ++i) {
        vec3_t u = {xyzw[0][i].x - com[0].x, xyzw[0][i].y - com[0].y, xyzw[0][i].z - com[0].z};
		vec3_t v = {xyzw[1][i].x - com[1].x, xyzw[1][i].y - com[1].y, xyzw[1][i].z - com[1].z};
		vec3_t vp = mat3_mul_vec3(R, v);
		vec3_t d = vec3_sub(u, vp);
		float w = (xyzw[0][i].w + xyzw[1][i].w) * 0.5f;
		d_sum += w * vec3_dot(d, d);
		w_sum += w;
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

    vec3_t ext = pbc_ext;
    vec3_t inv_ext = vec3_div(vec3_set1(1.0f), pbc_ext);

    int pbc_mask = simd_xyz_mask(pbc_ext.elem);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        md_256 x0 = md_mm256_loadu_ps(src_coord[0].x + i);
        md_256 y0 = md_mm256_loadu_ps(src_coord[0].y + i);
        md_256 z0 = md_mm256_loadu_ps(src_coord[0].z + i);

        md_256 x1 = md_mm256_loadu_ps(src_coord[1].x + i);
        md_256 y1 = md_mm256_loadu_ps(src_coord[1].y + i);
        md_256 z1 = md_mm256_loadu_ps(src_coord[1].z + i);

        switch (pbc_mask)
        {
        case 1:
            x1 = simd_deperiodize(x1, x0, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            break;
        case 2:
            y1 = simd_deperiodize(y1, y0, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            break;
        case 3:
            x1 = simd_deperiodize(x1, x0, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            y1 = simd_deperiodize(y1, y0, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            break;
        case 4:
            z1 = simd_deperiodize(z1, z0, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        case 5:
            x1 = simd_deperiodize(x1, x0, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            z1 = simd_deperiodize(z1, z0, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        case 6:
            y1 = simd_deperiodize(y1, y0, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            z1 = simd_deperiodize(z1, z0, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        case 7:
            x1 = simd_deperiodize(x1, x0, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            y1 = simd_deperiodize(y1, y0, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            z1 = simd_deperiodize(z1, z0, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        default:
            break;
        }

        md_256 x = md_mm256_lerp_ps(x0, x1, t);
        md_256 y = md_mm256_lerp_ps(y0, y1, t);
        md_256 z = md_mm256_lerp_ps(z0, z1, t);

        md_mm256_storeu_ps(dst_coord.x + i, x);
        md_mm256_storeu_ps(dst_coord.y + i, y);
        md_mm256_storeu_ps(dst_coord.z + i, z);
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

    vec3_t ext = pbc_ext;
    vec3_t inv_ext = vec3_div(vec3_set1(1.0f), pbc_ext);

    int pbc_mask = simd_xyz_mask(pbc_ext.elem);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        md_256 x0 = md_mm256_loadu_ps(src_coord[0].x + i);
        md_256 y0 = md_mm256_loadu_ps(src_coord[0].y + i);
        md_256 z0 = md_mm256_loadu_ps(src_coord[0].z + i);

        md_256 x1 = md_mm256_loadu_ps(src_coord[1].x + i);
        md_256 y1 = md_mm256_loadu_ps(src_coord[1].y + i);
        md_256 z1 = md_mm256_loadu_ps(src_coord[1].z + i);

        md_256 x2 = md_mm256_loadu_ps(src_coord[2].x + i);
        md_256 y2 = md_mm256_loadu_ps(src_coord[2].y + i);
        md_256 z2 = md_mm256_loadu_ps(src_coord[2].z + i);

        md_256 x3 = md_mm256_loadu_ps(src_coord[3].x + i);
        md_256 y3 = md_mm256_loadu_ps(src_coord[3].y + i);
        md_256 z3 = md_mm256_loadu_ps(src_coord[3].z + i);

        switch (pbc_mask) {
        case 1:
            x0 = simd_deperiodize(x0, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x2 = simd_deperiodize(x2, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x3 = simd_deperiodize(x3, x2, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            break;
        case 2:
            y0 = simd_deperiodize(y0, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y2 = simd_deperiodize(y2, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y3 = simd_deperiodize(y3, y2, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            break;
        case 3:
            x0 = simd_deperiodize(x0, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x2 = simd_deperiodize(x2, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x3 = simd_deperiodize(x3, x2, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            y0 = simd_deperiodize(y0, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y2 = simd_deperiodize(y2, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y3 = simd_deperiodize(y3, y2, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            break;
        case 4:
            z0 = simd_deperiodize(z0, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z2 = simd_deperiodize(z2, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z3 = simd_deperiodize(z3, z2, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        case 5:
            x0 = simd_deperiodize(x0, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x2 = simd_deperiodize(x2, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x3 = simd_deperiodize(x3, x2, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            z0 = simd_deperiodize(z0, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z2 = simd_deperiodize(z2, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z3 = simd_deperiodize(z3, z2, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        case 6:
            y0 = simd_deperiodize(y0, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y2 = simd_deperiodize(y2, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y3 = simd_deperiodize(y3, y2, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            z0 = simd_deperiodize(z0, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z2 = simd_deperiodize(z2, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z3 = simd_deperiodize(z3, z2, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        case 7:
            x0 = simd_deperiodize(x0, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x2 = simd_deperiodize(x2, x1, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            x3 = simd_deperiodize(x3, x2, md_mm256_set1_ps(ext.x), md_mm256_set1_ps(inv_ext.x));
            y0 = simd_deperiodize(y0, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y2 = simd_deperiodize(y2, y1, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            y3 = simd_deperiodize(y3, y2, md_mm256_set1_ps(ext.y), md_mm256_set1_ps(inv_ext.y));
            z0 = simd_deperiodize(z0, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z2 = simd_deperiodize(z2, z1, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            z3 = simd_deperiodize(z3, z2, md_mm256_set1_ps(ext.z), md_mm256_set1_ps(inv_ext.z));
            break;
        default:
            break;
        }

        const md_256 x = md_mm256_cubic_spline_ps(x0, x1, x2, x3, md_mm256_set1_ps(t), md_mm256_set1_ps(s));
        const md_256 y = md_mm256_cubic_spline_ps(y0, y1, y2, y3, md_mm256_set1_ps(t), md_mm256_set1_ps(s));
        const md_256 z = md_mm256_cubic_spline_ps(z0, z1, z2, z3, md_mm256_set1_ps(t), md_mm256_set1_ps(s));

        md_mm256_storeu_ps(dst_coord.x + i, x);
        md_mm256_storeu_ps(dst_coord.y + i, y);
        md_mm256_storeu_ps(dst_coord.z + i, z);
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
// (resid)                  -> Residues
// (Types & resname)        -> Elements
// (Elements)               -> Mass & Radius
// (Coordinates & Elements) -> Covalent Bonds
// (resid & Bonds)          -> Chains
// (resname & types)        -> Backbone
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

    if (flags & MD_UTIL_POSTPROCESS_RESIDUE_BIT) {
        // Create residues from resids
        if (mol->residue.count == 0 && mol->atom.resid) {
			md_util_compute_residue_data(&mol->residue, &mol->atom, alloc);
		}
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
        mol->bond = md_util_compute_covalent_bonds(&mol->atom, &mol->unit_cell, alloc);
    }

    if (flags & MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT) {
        if (mol->bond.pairs) {
            md_compute_connectivity(&mol->conn, &mol->atom, &mol->bond, alloc);
        }
    }

    if (flags & MD_UTIL_POSTPROCESS_STRUCTURE_BIT) {
        if (mol->conn.count) {
            mol->structures = md_util_compute_structures(mol, alloc);
            mol->rings      = md_util_compute_rings(mol, alloc);
		}
    }

    if (flags & MD_UTIL_POSTPROCESS_CHAINS_BIT) {
        if (mol->chain.count == 0 && mol->residue.count > 0 && mol->bond.pairs) {
            md_util_compute_chain_data(&mol->chain, &mol->atom, &mol->residue, &mol->bond, alloc);
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
This is blatantly stolen and modified from 'Arseny Kapoulkines' goldnugget 'Mesh optimizer'
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

typedef struct graph_t {
    int64_t   vertex_count;
	uint8_t*  vertex_type;
    uint32_t* edge_offset;  // offset, length is implicitly encoded by the next offset, last offset is the total number of edges and therefore length is count + 1
    uint32_t* edge_data;    // packed 32-bit data consiting of (from hi to low) grow : 1, type : 7, index : 24      
} graph_t;

static int64_t graph_vertex_count(const graph_t* g) {
	return g->vertex_count;
}

static int graph_vertex_type(const graph_t* g, int64_t vidx) {
	return g->vertex_type[vidx];
}

static int64_t graph_vertex_edge_count(const graph_t* g, int64_t vidx) {
	return g->edge_offset[vidx + 1] - g->edge_offset[vidx];
}

typedef struct graph_edge_iter_t {
    uint32_t* cur;
    uint32_t* end;
} graph_edge_iter_t;

static graph_edge_iter_t graph_edge_iter(const graph_t* g, int64_t vidx) {
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

static inline bool graph_edge_iter_grow(graph_edge_iter_t it) {
    return (*it.cur) >> 31;
}

static inline int graph_edge_iter_type(graph_edge_iter_t it) {
	return ((*it.cur) >> 24) & 0xff;
}

static inline int graph_edge_iter_vidx(graph_edge_iter_t it) {
	return (*it.cur) & 0xffffff;
}

static bool graph_vertex_is_connected_to(const graph_t* g, int vidx, int other_vidx) {
    graph_edge_iter_t it = graph_edge_iter(g, vidx);
    while (graph_edge_iter_has_next(it)) {
        if (graph_edge_iter_vidx(it) == other_vidx) return true;
        graph_edge_iter_next(&it);
    }
    return false;
}

static bool graph_vertex_has_connection(const graph_t* g, int vidx, int other_vidx, int other_type) {
    graph_edge_iter_t it = graph_edge_iter(g, vidx);
    while (graph_edge_iter_has_next(it)) {
		if (graph_edge_iter_vidx(it) == other_vidx &&
            graph_edge_iter_type(it) == other_type) return true;
        graph_edge_iter_next(&it);
    }
	return false;
}

static bool graph_equivalent(const graph_t* a, const graph_t* b) {
	if (a->vertex_count != b->vertex_count) return false;
	for (int i = 0; i < a->vertex_count; ++i) {
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

// VF2 Isomorphism algorithm
typedef struct candidate_t {
    int n_idx;
    int h_idx;
} candidate_t;

static bool candidate_equal(candidate_t a, candidate_t b) {
    return MEMCMP(&a, &b, sizeof(candidate_t)) == 0;
}

typedef bool (*solution_callback)(const int map[], int64_t length, void* user);

typedef struct state_t {
    bool abort;

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

    md_allocator_i* alloc;
} state_t;

static void state_init(state_t* state, const graph_t* n_graph, const graph_t* h_graph, md_allocator_i* alloc) {
	state->abort = false;
	state->n_graph = n_graph;
	state->h_graph = h_graph;
	state->n_path = 0;
	state->h_path = 0;
	state->map = md_array_create(int, n_graph->vertex_count, alloc);
	state->n_path_bits = make_bitfield(n_graph->vertex_count, alloc);
	state->h_path_bits = make_bitfield(h_graph->vertex_count, alloc);
	state->n_depths = 0;
	state->h_depths = 0;

    md_array_ensure(state->n_path, n_graph->vertex_count, alloc);
    md_array_ensure(state->h_path, h_graph->vertex_count, alloc);
    state->n_depths = md_array_create(uint32_t, n_graph->vertex_count, alloc);
    state->h_depths = md_array_create(uint32_t, h_graph->vertex_count, alloc);

    MEMSET(state->n_depths, 0, md_array_bytes(state->n_depths));
    MEMSET(state->h_depths, 0, md_array_bytes(state->h_depths));

    state->alloc = alloc;
}

static void state_reset(state_t* state) {
    state->abort = false;
    md_array_shrink(state->n_path, 0);
    md_array_shrink(state->h_path, 0);
    MEMSET(state->map, -1, md_array_bytes(state->map));
    clear_all_bitfield(state->n_path_bits);
    clear_all_bitfield(state->h_path_bits);
    MEMSET(state->n_depths, 0, md_array_bytes(state->n_depths));
    MEMSET(state->h_depths, 0, md_array_bytes(state->h_depths));
}

static void state_free(state_t* state) {
    md_array_free(state->map, state->alloc);
	md_array_free(state->n_path, state->alloc);
	md_array_free(state->h_path, state->alloc);
	md_array_free(state->n_path_bits, state->alloc);
	md_array_free(state->h_path_bits, state->alloc);
	md_array_free(state->n_depths, state->alloc);
	md_array_free(state->h_depths, state->alloc);
}

static void backtrack(state_t* state) {
    if (md_array_size(state->n_path)) {
        int n_idx = md_array_pop(state->n_path);
        state->map[n_idx] = -1;
        clear_bit(state->n_path_bits, n_idx);
    }
    if (md_array_size(state->h_path)) {
        int h_idx = md_array_pop(state->h_path);
        clear_bit(state->h_path_bits, h_idx);
    }

    uint32_t depth = (uint32_t)md_array_size(state->n_path) + 1;
    for (int64_t i = 0; i < md_array_size(state->n_depths); ++i) {
        //state->n_depths[i] = (state->n_depths[i] == depth) ? 0 : state->n_depths[i];
        if (state->n_depths[i] == depth) {
			state->n_depths[i] = 0;
		}
    }
    for (int64_t i = 0; i < md_array_size(state->h_depths); ++i) {
        //state->h_depths[i] = (state->h_depths[i] == depth) ? 0 : state->h_depths[i];
        if (state->h_depths[i] == depth) {
            state->h_depths[i] = 0;
        }
    }
}

static bool is_in_terminal_set(const uint32_t depths[], const uint64_t path_bits[], int64_t i) {
    if (test_bit(path_bits, i)) return false;
    if (!depths[i]) return false;
    return true;
}

static bool check_bonds(state_t* state, int n_idx) {
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
static bool check_map(state_t* state) {
    if (md_array_size(state->n_path) != state->n_graph->vertex_count) {
        return false;
    }

    return state->callback(state->map, md_array_size(state->map), state->user_data);
}

static bool test_edge_constraints(state_t* state, int n_idx, int h_idx) {
    graph_edge_iter_t it = graph_edge_iter(state->n_graph, n_idx);
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

static bool match_candidate(state_t* state, int n_idx, int h_idx) {
	if (graph_vertex_type(state->n_graph, n_idx) != graph_vertex_type(state->h_graph, h_idx)) {
		return false;
	}

    state->map[n_idx] = h_idx;

    if (!test_edge_constraints(state, n_idx, h_idx)) {
        state->map[n_idx] = -1;
        return false;
    }

    md_array_push(state->n_path, n_idx, state->alloc);
    md_array_push(state->h_path, h_idx, state->alloc);

    set_bit(state->n_path_bits, n_idx);
    set_bit(state->h_path_bits, h_idx);

    // Update n_depths
    {
        const uint32_t d = (uint32_t)md_array_size(state->n_path);
        if (!state->n_depths[n_idx]) {
            state->n_depths[n_idx] = d;
        }

        graph_edge_iter_t it = graph_edge_iter(state->n_graph, n_idx);
        while (graph_edge_iter_has_next(it)) {
    	    int vidx = graph_edge_iter_vidx(it);
            if (!state->n_depths[vidx]) {
			    state->n_depths[vidx] = d;
		    }
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
    const int64_t N1 = state->n_graph->vertex_count;
    const int64_t N2 = state->h_graph->vertex_count;

    int64_t T1 = 0;
    for (int64_t i = 0; i < N1; ++i) {
        if (is_in_terminal_set(state->n_depths, state->n_path_bits, i)) {
			++T1;
		}
    }

    int64_t T2 = 0;
    for (int64_t i = 0; i < N2; ++i) {
		if (is_in_terminal_set(state->h_depths, state->h_path_bits, i)) {
            ++T2;
        }
    }

    if (T1 > T2) {
        backtrack(state);
        return false;
    }

    const int64_t M1 = md_array_size(state->n_path);
    const int64_t M2 = md_array_size(state->h_path);

    if ((N1 - M1 - T1) > (N2 - M2 - T2)) {
		backtrack(state);
		return false;
	}

    state->abort = check_map(state);

	return true;
}

static int terminal_size(const uint32_t depths[], int64_t len) {
    int count = 0;
    for (int64_t i = 0; i < len; ++i) {
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
        while (last_n_idx < n_size && (test_bit(state->n_path_bits, last_n_idx) || !state->n_depths[last_n_idx])) {
			++last_n_idx;
            last_h_idx = 0;
		}
    } else {
        while (last_n_idx < n_size && test_bit(state->n_path_bits, last_n_idx)) {
            ++last_n_idx;
			last_h_idx = 0;
        }
    }

    if (n_term_size > map_size && h_term_size > map_size) {
        while (last_h_idx < h_size && (test_bit(state->h_path_bits, last_h_idx) || !state->h_depths[last_h_idx])) {
            last_h_idx++;
        }
    } else {
        while (last_h_idx < h_size && test_bit(state->h_path_bits, last_h_idx)) {
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

#if DEBUG
        const int depth = (int)md_array_size(state->n_path);
        printf("%*s [%d, %d]\n", depth, "", n_idx, h_idx);
#endif

        if (match_candidate(state, n_idx, h_idx)) {
			map_next(state);
			backtrack(state);
        }
    }
}

// This is used to map whatever mapping type we have to a unique integer
typedef struct type_map_entry_t {
    size_t key;
    size_t value;
} type_map_entry_t;

typedef struct result_entry_t {
    size_t key;
    size_t value;
    int result_idx;
} result_entry_t;

typedef struct subgraph_context_t {
	const graph_t* n_graph;
	const graph_t* h_graph;
	
	md_array(int) h_idx_map;
    md_array(uint64_t) h_occupied;
    md_array(uint64_t) n_available;

    result_entry_t* result_map; // This is to govern the permutations and to only store a single permutation for a given set of indices
    md_index_data_t* result;
    md_allocator_i*  alloc;
    md_allocator_i*  temp_alloc;
} subgraph_context_t;

static void print_solution(const subgraph_context_t* ctx) {
    md_strb_t sb = md_strb_create(ctx->temp_alloc);
    for (int64_t i = 0; i < ctx->n_graph->vertex_count; ++i) {
        md_strb_fmt(&sb, "%d", ctx->h_idx_map[i]);
        if (i < ctx->n_graph->vertex_count - 1) {
            md_strb_push_cstr(&sb, ", ");
        }
    }
    MD_LOG_DEBUG("Found solution: {%s}", md_strb_to_cstr(&sb));
    md_strb_free(&sb);
}

static bool check_constraints(int n_idx, int h_idx, const subgraph_context_t* ctx) {
    int n_type = graph_vertex_type(ctx->n_graph, n_idx);
    int h_type = graph_vertex_type(ctx->h_graph, h_idx);
    
    // Vertex type constraint
    if (n_type != h_type) return false;

    // Edge constraints
    graph_edge_iter_t n_it = graph_edge_iter(ctx->n_graph, n_idx);
    while (graph_edge_iter_has_next(n_it)) {
        int type = graph_edge_iter_type(n_it);
        int n_other_idx = graph_edge_iter_vidx(n_it);
        int h_other_idx = ctx->h_idx_map[n_other_idx];
        if (h_other_idx != -1) {
            if (type == 0) {
                if (!graph_vertex_is_connected_to(ctx->h_graph, h_idx, h_other_idx)) {
					return false;
				}
            } else {
                if (!graph_vertex_has_connection(ctx->h_graph, h_idx, h_other_idx, type)) {
                    return false;
                }
            }
        }
        graph_edge_iter_next(&n_it);
    }

    return true;
}

static bool get_next_pair(int* n_idx, int* h_idx, const subgraph_context_t* ctx) {
    for (int i = 0; i < (int)ctx->n_graph->vertex_count; ++i) {
        if (ctx->h_idx_map[i] > -1) {
            int n_par = i;
            int h_par = ctx->h_idx_map[i];
            graph_edge_iter_t n_it = graph_edge_iter(ctx->n_graph, n_par);
            while (graph_edge_iter_has_next(n_it)) {
                int n = graph_edge_iter_vidx(n_it);
                if (ctx->h_idx_map[n] == -1) {
                    // Unconnected edge in n_par
                    // Must match with an unconnected edge in h_par
                    int n_type = ctx->n_graph->vertex_type[n];
                    int n_edge_type = graph_edge_iter_type(n_it);

                    graph_edge_iter_t h_it = graph_edge_iter(ctx->h_graph, h_par);
                    while (graph_edge_iter_has_next(h_it)) {
                        int h = graph_edge_iter_vidx(h_it);
                        if (!test_bit(ctx->h_occupied, h)) {
                            int h_type = ctx->h_graph->vertex_type[h];
                            int h_edge_type = graph_edge_iter_type(h_it);
                            if (n_type == h_type) {
                                if (n_type == 0 || n_edge_type == h_edge_type) {
                                    *n_idx = n;
                                    *h_idx = h;
                                    return true;
                                }
                            }
                        }
                        graph_edge_iter_next(&h_it);
                    }
                }
                graph_edge_iter_next(&n_it);
            }
        }
    }
    return false;
}

static int get_next_n_idx(const subgraph_context_t* ctx) {
    for (int i = 0; i < (int)ctx->n_graph->vertex_count; ++i) {
        if (ctx->h_idx_map[i] > -1) {
            graph_edge_iter_t it = graph_edge_iter(ctx->n_graph, i);
            while (graph_edge_iter_has_next(it)) {
                int n_idx = graph_edge_iter_vidx(it);
                if (ctx->h_idx_map[n_idx] == -1) {
					return n_idx;
				}
                graph_edge_iter_next(&it);
            }
        }
    }

    return -1;

#if 0
    // Disjoint set if we end up here
	for (int i = 0; i < (int)ctx->n_graph->count; ++i) {
		if (ctx->h_idx_map[i] == -1) {
			return i;
		}
	}
	return -1;
#endif
}

static bool next_n_idx(int* n_idx, const subgraph_context_t* ctx) {
    int i;
    if (!find_first_bit_set_bitfield(&i, ctx->n_available)) {
        return false;
    }
    for (i = 0; i < (int)ctx->n_graph->vertex_count; ++i) {
        if (ctx->h_idx_map[i] > -1) {
            graph_edge_iter_t it = graph_edge_iter(ctx->n_graph, i);
            while (graph_edge_iter_has_next(it)) {
                int n = graph_edge_iter_vidx(it);
                if (ctx->h_idx_map[n] == -1) {
                    *n_idx = n;
                    return true;
                }
                graph_edge_iter_next(&it);
            }
        }
    }
    return false;
}

static int get_next_h_idx(int n_idx, int h_idx, const subgraph_context_t* ctx) {
	for (int i = h_idx + 1; i < (int)ctx->h_graph->vertex_count; ++i) {
		if (test_bit(ctx->h_occupied, i)) continue;
        if (check_constraints(n_idx, i, ctx)) {
			return i;
		}
	}
	return -1;
}

static bool next_h_idx(int* h_idx, int n_idx, const subgraph_context_t* ctx) {
    for (int i = *h_idx + 1; i < (int)ctx->h_graph->vertex_count; ++i) {
        if (test_bit(ctx->h_occupied, i)) continue;
        if (check_constraints(n_idx, i, ctx)) {
            *h_idx = i;
            return true;
        }
    }
    return false;
}

static void record_solution(subgraph_context_t* ctx) {
    size_t key = 0;
    size_t value = 0;
    for (size_t i = 0; i < (size_t)md_array_size(ctx->h_idx_map); ++i) {
        key += ctx->h_idx_map[i];
        value = ctx->h_idx_map[i] * (i + 1);
    }
    result_entry_t* entry = stbds_hmgetp_null(ctx->result_map, key);
    if (entry) {
        if (entry->value < value) {
            entry->value = value;
            int* ptr = md_index_range_beg(*ctx->result, entry->result_idx);
            MEMCPY(ptr, ctx->h_idx_map, sizeof(int) * ctx->n_graph->vertex_count);
        }
    } else {
		int result_idx = (int)md_index_data_push(ctx->result, ctx->h_idx_map, ctx->n_graph->vertex_count, ctx->alloc);
        result_entry_t e = {key, value, result_idx};
		stbds_hmputs(ctx->result_map, e);
	}
    
}

static void mark_slot(int n_idx, int h_idx, subgraph_context_t* ctx) {
    ctx->h_idx_map[n_idx] = h_idx;
    set_bit(ctx->h_occupied, h_idx);
    clear_bit(ctx->n_available, n_idx);
}

static void unmark_slot(int n_idx, int h_idx, subgraph_context_t* ctx) {
    set_bit(ctx->n_available, n_idx);
    clear_bit(ctx->h_occupied, h_idx);
    ctx->h_idx_map[n_idx] = -1;
}

static bool is_solution(subgraph_context_t* ctx) {
    for (int64_t i = 0; i < (int)ctx->n_graph->vertex_count; ++i) {
		if (ctx->h_idx_map[i] == -1) return false;
	}
    return true;
}

#define DEBUG_PRINT DEBUG

static void brute_force(subgraph_context_t* ctx, int depth) {
    if (is_solution(ctx)) {
		record_solution(ctx);
#if DEBUG_PRINT
        print_solution(ctx);
#endif
        return;
    }

#if 0
    int n_idx;
    int h_idx;
    if (!get_next_pair(&n_idx, &h_idx, ctx)) {
        return;
    }

    // Find potential mapping for n_idx
    do {
        mark_slot(n_idx, h_idx, ctx);

#if DEBUG_PRINT
        printf("%*s [%d, %d]\n", depth, " ", n_idx, h_idx);
#endif
        brute_force(ctx, depth + 1);

        unmark_slot(n_idx, h_idx, ctx);
    } while (next_h_idx(&h_idx, n_idx, ctx));
#else 

    int n_idx;
    next_n_idx(&n_idx, ctx);

    int h_idx = -1;
    while (next_h_idx(&h_idx, n_idx, ctx)) {
		mark_slot(n_idx, h_idx, ctx);

#if DEBUG_PRINT
		printf("%*s [%d, %d]\n", depth, "", n_idx, h_idx);
#endif
        brute_force(ctx, depth + 1);

		unmark_slot(n_idx, h_idx, ctx);
	}

#endif
}

typedef struct pair_t{
    int n_idx;
    int h_idx;
} pair_t;

static void brute_force_non_recursive(subgraph_context_t* ctx) {
    md_array(pair_t) stack = 0;

    {
        int n_idx = get_next_n_idx(ctx);
        md_array_push(stack, ((pair_t){n_idx,-1}), ctx->temp_alloc);
    }

    while (md_array_size(stack) > 0) {
    	pair_t* p = md_array_last(stack);
        if (p->h_idx != -1) {
            unmark_slot(p->n_idx, p->h_idx, ctx);
        }

        p->h_idx = get_next_h_idx(p->n_idx, p->h_idx, ctx);

        if (p->h_idx == -1) {
            int h_idx = ctx->h_idx_map[p->n_idx];
            unmark_slot(p->n_idx, h_idx, ctx);
            md_array_pop(stack);
        } else {
            mark_slot(p->n_idx, p->h_idx, ctx);

			int n_idx = get_next_n_idx(ctx);
            if (n_idx == -1) {
                record_solution(ctx);
                unmark_slot(p->n_idx, p->h_idx, ctx);
                md_array_pop(stack);
            } else {
			    md_array_push(stack, ((pair_t){n_idx,-1}), ctx->temp_alloc);
            }
		}
    }
}

typedef struct store_data_t {
    md_index_data_t* result;
    md_allocator_i*  alloc;
    // Only for unique
    result_entry_t* map;
} store_data_t;

static bool store_unique_callback(const int map[], int64_t length, void* user) {
    store_data_t* data = (store_data_t*)user;

    size_t key = 0;
    size_t value = 0;
    for (size_t i = 0; i < (size_t)length; ++i) {
        key += map[i];
        value = map[i] * (i + 1);
    }
    result_entry_t* entry = stbds_hmgetp_null(data->map, key);
    if (entry) {
        if (entry->value < value) {
            // Update result if the new value is lower (more strictly sorted)
            entry->value = value;
            int* ptr = md_index_range_beg(*data->result, entry->result_idx);
            MEMCPY(ptr, map, sizeof(int) * length);
        }
    } else {
        // Store new entry if unique solution
        int result_idx = (int)md_index_data_push(data->result, map, length, data->alloc);
        result_entry_t e = {key, value, result_idx};
        stbds_hmputs(data->map, e);
    }

#if DEBUG
    const int depth = (int)length;
    printf("%*sSolution found!\n", depth, "");
#endif

    return false;
}

static bool store_first_callback(const int map[], int64_t length, void* user) {
    store_data_t* data = (store_data_t*)user;

    md_index_data_push(data->result, map, length, data->alloc);
#if DEBUG
    const int depth = (int)length;
    printf("%*sSolution found!\n", depth, "");
#endif

    return true;
}

static bool store_all_callback(const int map[], int64_t length, void* user) {
    store_data_t* data = (store_data_t*)user;

    md_index_data_push(data->result, map, length, data->alloc);
#if DEBUG
    const int depth = (int)length;
    printf("%*sSolution found!\n", depth, "");
#endif

    return false;
}

static bool store_count_callback(const int map[], int64_t length, void* user) {
    (void)map;
    (void)length;
	int64_t* count = (int64_t*)user;
	++(*count);
	return false;
}

// Attempt to find subgraphs in a larger graph (haystack) which match a reference graph (needle)
// Returns an array of graphs which match the topologically match the reference
// start_type is a hint of the most unusual type in the graphs and serve as good starting points
static md_index_data_t find_isomorphisms(const graph_t* needle, const graph_t* haystack, md_util_match_mode_t mode, int start_type, md_allocator_i* alloc) {
    md_index_data_t result = {0};

    // Impossible case
    if (needle->vertex_count > haystack->vertex_count) {
        return result;
    }

    // The not so problematic case
#if 0
    if (needle->count == haystack->count) {
        if (graph_equivalent(needle, haystack)) {
            // This should be a 1:1 mapping
            md_array_resize(result.indices, haystack->count, alloc);
            for (int i = 0; i < (int)haystack->count; ++i) {
                result.indices[i] = i;
            }
            md_array_resize(result.ranges, 1, alloc);
            result.ranges[0] = (md_range_t){0, (int)haystack->count};
		}
        return result;
    }
#endif

    // The problematic case
    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    // Create list of starting candidate pairs
    md_array(pair_t) start_candidates = 0;
    for (int i = 0; i < needle->vertex_count; ++i) {
        if (graph_vertex_type(needle, i) != start_type) continue;
        for (int j = 0; j < haystack->vertex_count; ++j) {
            if (graph_vertex_type(haystack, j) != start_type) continue;
            md_array_push(start_candidates, ((pair_t){i,j}), temp_alloc);
        }
    }

    if (start_candidates == 0) goto done;
    
    // Corresponding index array which needle indices correspond to the indices within the haystack;
    md_array(int) corr_idx = md_array_create(int, needle->vertex_count, temp_alloc);
    md_array(uint64_t) h_occupied  = make_bitfield(haystack->vertex_count, temp_alloc);
    md_array(uint64_t) n_available = make_bitfield(needle->vertex_count, temp_alloc);

    subgraph_context_t ctx = {
		.n_graph = needle,
		.h_graph = haystack,
		.h_idx_map = corr_idx,
        .h_occupied = h_occupied,
        .n_available = n_available,
		.result = &result,
		.alloc = alloc,
		.temp_alloc = temp_alloc,
	};

    state_t state = {0};
    state_init(&state, needle, haystack, temp_alloc);

    store_data_t store_data = {&result, alloc};
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

    const int64_t num_candidates = md_array_size(start_candidates);
    for (int64_t i = 0; i < num_candidates; ++i) {
        int n_idx = start_candidates[i].n_idx;
        int h_idx = start_candidates[i].h_idx;

#if 0
        printf("STARTING ATTEMPT TO MATCH %d -> %d\n", n_idx, h_idx);
#endif

#if 0
        // Reset state
        for (int64_t k = 0; k < needle->vertex_count; ++k) {
            corr_idx[k] = -1;
        }
        clear_all_bitfield(ctx.h_occupied);
        set_all_bitfield(ctx.n_available, ctx.n_graph->vertex_count);

        // Set initial state
        mark_slot(n_idx, h_idx, &ctx);
        brute_force(&ctx, 0);

        //brute_force_non_recursive(&ctx);

#else

        // Reset state
		state_reset(&state);

		// Set initial state
		if (match_candidate(&state, n_idx, h_idx)) {
			map_next(&state);
		}
#endif

        if (mode == MD_UTIL_MATCH_MODE_FIRST && md_index_data_count(result) > 0) {
            goto done;
        }
    }
    
done:    
    md_arena_allocator_destroy(temp_alloc);
    return result;
}

enum {
    VERTEX_TYPE_MAPPING_UNDEFINED = 0,
    VERTEX_TYPE_MAPPING_TYPE = 1,
    VERTEX_TYPE_MAPPING_ELEMENT = 2,
    VERTEX_TYPE_MAPPING_MASS = 3
};

typedef int vertex_type_mapping_t;

// Confusing function
// Extracts a graph from an atom index range with supplied atom types

static graph_t make_graph(const md_molecule_t* mol, const uint8_t atom_types[], const int indices[], int64_t count, md_allocator_i* alloc) {   
    graph_t graph = {0};

    // We create a map from the global atom indices to new indices local to the graph
    struct { int key; int value; } *map = NULL;
    stbds_hmdefault(map, -1);

    graph.vertex_count = count;
    graph.vertex_type = md_array_create(uint8_t, count, alloc);

    for (int64_t i = 0; i < count; ++i) {
        int idx = indices[i];
        stbds_hmput(map, idx, (int)i);
        graph.vertex_type[i] = atom_types[idx];
    }

    graph.edge_offset = md_array_create(uint32_t, count + 1, alloc);
    MEMSET(graph.edge_offset, 0, md_array_bytes(graph.edge_offset));

    // Only store edges which point to vertexes within the graph as this will be used later as a traversal template
    for (int64_t i = 0; i < count; ++i) {
        int idx = indices[i];
        // Translate the global atom indices to local structure indices
        uint32_t edge_data_arr[8];
        uint32_t length = 0;

        if (i == 238) {
            while(0);
        }

        md_conn_iter_t it = md_conn_iter(mol, idx);
        while (md_conn_iter_has_next(&it)) {
            int index = md_conn_iter_index(&it);
            int order = md_conn_iter_order(&it);
            int local_idx = stbds_hmget(map, index);
            if (local_idx != -1) {
                edge_data_arr[length++] = (uint32_t)(order << 24) | (uint32_t)local_idx;
            }
            md_conn_iter_next(&it);
        }

        // Sort on indices
        sort_arr_masked((int*)edge_data_arr, length, 0x00FFFFFF);

        if (length > 0) {
            graph.edge_offset[i] = length;
            md_array_push_array(graph.edge_data, edge_data_arr, length, alloc);
        }
    }

    // exclusive scan to set the offsets
    uint32_t offset = 0;
    for (int64_t i = 0; i < count + 1; ++i) {
        uint32_t len = graph.edge_offset[i];
        graph.edge_offset[i] = offset;
        offset += len;
	}

    stbds_hmfree(map);
    return graph;
}

// Create a new reference structure which is pruned of certain atoms (Hydrogen) and loosely connected subcomponents
// There are simply too many permutations to cover and the result will explode.
static md_array(int) filter_structure_connectivity(const int* indices, int64_t count, const md_index_data_t* connectivity, int min_val, md_allocator_i* alloc) {
	md_array(int) filt = 0;
	for (int64_t i = 0; i < count; ++i) {
		int idx = indices[i];
		int64_t order = md_index_range_size(*connectivity, idx);
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

md_index_data_t get_structures(const md_molecule_t* mol, md_util_match_level_t level, bool filter_hydrogen, md_allocator_i* alloc) {
    md_index_data_t result = {0};

    switch (level) {
    case MD_UTIL_MATCH_LEVEL_STRUCTURE:
        if (filter_hydrogen) {
            for (int64_t s_idx = 0; s_idx < md_index_data_count(mol->structures); ++s_idx) {
                md_array(int) structure = 0;
                for (int* it = md_index_range_beg(mol->structures, s_idx); it < md_index_range_end(mol->structures, s_idx); ++it) {
                    int i = *it;
                    if (filter_hydrogen && mol->atom.element[i] == H) {
                        continue;
                    }
                    md_array_push(structure, i, alloc);
                }
                md_index_data_push(&result, structure, md_array_size(structure), alloc);
            }
        } else {
            result = mol->structures;
        }
        break;
    case MD_UTIL_MATCH_LEVEL_RESIDUE:
        for (int64_t r_idx = 0; r_idx < mol->residue.count; ++r_idx) {
            md_array(int) structure = 0;
            for (int i = mol->residue.atom_range[r_idx].beg; i < mol->residue.atom_range[r_idx].end; ++i) {
                if (filter_hydrogen && mol->atom.element[i] == H) {
                    continue;
                }
                md_array_push(structure, i, alloc);
            }
            md_index_data_push(&result, structure, md_array_size(structure), alloc);
        }
        break;
    case MD_UTIL_MATCH_LEVEL_CHAIN:
        for (int64_t c_idx = 0; c_idx < mol->chain.count; ++c_idx) {
            md_array(int) structure = 0;
            for (int i = mol->chain.atom_range[c_idx].beg; i < mol->chain.atom_range[c_idx].end; ++i) {
                if (filter_hydrogen && mol->atom.element[i] == H) {
                    continue;
                }
                md_array_push(structure, i, alloc);
            }
            md_index_data_push(&result, structure, md_array_size(structure), alloc);
        }
        break;
    default:
        ASSERT(false);
    }
    return result;
}

md_index_data_t match_structure(const int* ref_idx, int64_t ref_len, md_util_match_mode_t mode, md_util_match_level_t level, vertex_type_mapping_t mapping, const md_molecule_t* mol, md_allocator_i* alloc) {
    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    graph_t ref_graph = {0};
    md_index_data_t result = {0};

    type_map_entry_t* map_table = 0;

    md_array(uint8_t) atom_type = md_array_create(uint8_t, mol->atom.count, temp_alloc);
    md_array(int) structure_idx = md_array_create(int, mol->atom.count, temp_alloc);

    if (mapping == VERTEX_TYPE_MAPPING_ELEMENT) {
        for (int i = 0; i < Num_Elements; ++i) {
            stbds_hmput(map_table, i, i);
        }
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            atom_type[i] = mol->atom.element[i];
        }
    } else {
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            size_t key;
            switch (mapping) {
            case VERTEX_TYPE_MAPPING_TYPE:
                key = stbds_hash_string(mol->atom.type[i].buf, 0);
                break;
            case VERTEX_TYPE_MAPPING_MASS:
                key = stbds_hash_bytes(&mol->atom.mass[i], sizeof(mol->atom.mass[i]), 0);
                break;
            default: ASSERT(false); break;
            }
            uint8_t type = 0;
            int64_t idx = stbds_hmgeti(map_table, key);
            if (idx == -1) {
                type = (uint8_t)stbds_hmlen(map_table);
                stbds_hmput(map_table, key, 1);
            } else {
                type = (uint8_t)idx;
                map_table[idx].value += 1;
            }
            ASSERT(type >= 0);
            atom_type[i] = type;
        }
    }

    for (int64_t i = 0; i < md_index_data_count(mol->structures); ++i) {
        int* beg = md_index_range_beg(mol->structures, i);
        int* end = md_index_range_end(mol->structures, i);
        for (int* it = beg; it != end; ++it) {
            structure_idx[*it] = (int)i;
        }
    }

    // Ensure that the reference indices all belong to the same structure and are not disjoint
    int ref_structure_idx = structure_idx[ref_idx[0]];
    for (int64_t i = 1; i < ref_len; ++i) {
        if (structure_idx[ref_idx[i]] != ref_structure_idx) {
            MD_LOG_ERROR("Reference indices are not part of the same structure, they are disconnected");
            goto done;
        }
    }

    bool ref_hydro_present = false;
    for (int64_t i = 0; i < ref_len; ++i) {
        int idx = ref_idx[i];
		if (mol->atom.element[idx] == H) {
			ref_hydro_present = true;
			break;
		}
    }

    ref_graph = make_graph(mol, atom_type, ref_idx, ref_len, temp_alloc);

    const int max_types = (int)stbds_hmlenu(map_table);

    int ref_type_count[256] = {0};
    for (int i = 0; i < ref_len; ++i) {
        uint8_t type = (uint8_t)graph_vertex_type(&ref_graph, i);
        ref_type_count[type]++;
    }

    md_index_data_t structures = get_structures(mol, level, !ref_hydro_present, temp_alloc);
    const int64_t num_structures = md_index_data_count(structures);

    int starting_type = -1;
    int min_freq = INT_MAX;
    for (int i = 0; i < max_types; ++i) {
        int freq = ref_type_count[i];
        if (freq > 0 && freq < min_freq) {
            min_freq = freq;
            starting_type = i;
        }
    }

    for (int64_t i = 0; i < num_structures; ++i) {

        md_range_t s_range = structures.ranges[i];
        const int* s_idx = structures.indices + s_range.beg;
        const int64_t s_len = s_range.end - s_range.beg;
        if (s_len < ref_len) continue;

        md_util_match_mode_t s_mode = mode;
        if (s_len == ref_len && s_mode == MD_UTIL_MATCH_MODE_UNIQUE) {
			s_mode = MD_UTIL_MATCH_MODE_FIRST;
		}
        
        int s_type_count[256] = {0};
        for (int64_t j = 0; j < s_len; ++j) {
            int idx = s_idx[j];
            uint8_t type = atom_type[idx];
            s_type_count[type]++;
        }

        // Sanity check
        for (int j = 0; j < max_types; ++j) {
            if (ref_type_count[j] > s_type_count[j]) {
                goto next;
            }
        }

        graph_t graph = make_graph(mol, atom_type, s_idx, s_len, temp_alloc);
        md_index_data_t mappings = find_isomorphisms(&ref_graph, &graph, s_mode, starting_type, temp_alloc);

        // Remap indices to global indices in result
        const int64_t num_mappings = md_index_data_count(mappings);
        for (int64_t j = 0; j < num_mappings; ++j) {
            int* beg = md_index_range_beg(mappings, j);
            int* end = md_index_range_end(mappings, j);
            int r_beg = (int)md_array_size(result.indices);
            int r_end = r_beg + (int)(end - beg);
            md_range_t range = {r_beg, r_end};
            md_array_push(result.ranges, range, alloc);

            for (int* it = beg; it != end; ++it) {
                int idx = s_idx[*it];
                md_array_push(result.indices, idx, alloc);
            }
        }
    next:;
    }

done:
    md_arena_allocator_destroy(temp_alloc);
    return result;
}

md_index_data_t md_util_match_by_type(const int ref_indices[], int64_t ref_count, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc) {
	return match_structure(ref_indices, ref_count, mode, level, VERTEX_TYPE_MAPPING_TYPE, mol, alloc);
}

md_index_data_t md_util_match_by_element(const int ref_indices[], int64_t ref_count, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc) {
    return match_structure(ref_indices, ref_count, mode, level, VERTEX_TYPE_MAPPING_ELEMENT, mol, alloc);
}

graph_t smiles_to_graph(str_t smiles_str, md_allocator_i* alloc) {
    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    typedef struct vertex_t {
        int type;
        int edge_count;
        int edge_idx[8];
        int edge_type[8];
    } vertex_t;

    md_array(vertex_t) verts = 0;
    md_array(md_smiles_node_t) nodes = md_array_create(md_smiles_node_t, str_len(smiles_str), temp_alloc);
    const int64_t num_nodes = md_smiles_parse(nodes, md_array_size(nodes), str_ptr(smiles_str), str_len(smiles_str));
    md_array_shrink(nodes, num_nodes);

    graph_t graph = {0};

    if (nodes) {
        struct {
           int key;
           int value;
        } *map = 0;
        stbds_hmdefault(map, -1);

        int stack[128];
        int stack_size = 0;
        int hub = -1;
        int order = 0;
        
        for (int64_t i = 0; i < num_nodes; ++i) {
            const md_smiles_node_t* node = &nodes[i];
            if (node->type == MD_SMILES_NODE_ATOM) {
                md_element_t elem = node->atom.element;
                md_array_push(verts, (vertex_t){.type = elem}, temp_alloc);
                int cur = (int)md_array_size(verts) - 1;

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

                for (int j = 0; j < node->atom.hydrogen_count; ++j) {
                    md_array_push(verts, (vertex_t){.type = H}, temp_alloc);
					const int h_idx = (int)md_array_size(verts) - 1;
                    const int edge_type = 1;

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
                case ':': order = 0; break;
                default: ASSERT(false); break;
                }
            }
            else if (node->type == MD_SMILES_NODE_BRANCH_OPEN) {
                if (stack_size < ARRAY_SIZE(stack)) {
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
                int idx = (int)stbds_hmgeti(map, node->bridge.index);
                int ci = (int)md_array_size(verts) - 1;
                if (idx != -1) {
					int pi = map[idx].value;
                    verts[ci].edge_idx[verts[ci].edge_count] = pi;
                    verts[ci].edge_type[verts[ci].edge_count] = 0;

                    verts[pi].edge_idx[verts[pi].edge_count] = ci;
                    verts[pi].edge_type[verts[pi].edge_count] = 0;

                    verts[ci].edge_count++;
                    verts[pi].edge_count++;

                    stbds_hmdel(map, idx);
                } else {
					stbds_hmput(map, node->bridge.index, ci);
				}
            }
        }
    }

    const int64_t num_verts = md_array_size(verts);
    if (num_verts > 0) {
		graph.vertex_count = num_verts;
        graph.vertex_type = md_array_create(uint8_t, num_verts, alloc);
        graph.edge_offset = md_array_create(uint32_t, num_verts + 1, alloc);
        MEMSET(graph.edge_offset, 0, md_array_bytes(graph.edge_offset));

        for (int64_t i = 0; i < num_verts; ++i) {
        	graph.vertex_type[i] = (uint8_t)verts[i].type;
            if (verts[i].edge_count > 0) {
                const uint32_t len = verts[i].edge_count;
                graph.edge_offset[i] = len;
				for (uint32_t j = 0; j < len; ++j) {
					const uint32_t idx  = verts[i].edge_idx[j];
					const uint32_t type = verts[i].edge_type[j] & 0x7F;
                    const uint32_t data = (type << 24) | idx;
					md_array_push(graph.edge_data, data, alloc);
				}
			}
        }

        // exclusive scan to set the offsets
        uint32_t offset = 0;
        for (int64_t i = 0; i < num_verts + 1; ++i) {
            uint32_t len = graph.edge_offset[i];
            graph.edge_offset[i] = offset;
            offset += len;
        }
    }

done:
    md_arena_allocator_destroy(temp_alloc);
    return graph;
}

md_index_data_t md_util_match_smiles(str_t smiles, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc) {
    ASSERT(mol);
    md_index_data_t result = {0};

    if (!mol->atom.element) {
        MD_LOG_ERROR("Molecule does not have any atom element field");
        return result;
    }

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    graph_t ref_graph = smiles_to_graph(smiles, temp_alloc);

    int ref_type_count[256] = {0};
    for (int64_t i = 0; i < ref_graph.vertex_count; ++i) {
		uint8_t type = ref_graph.vertex_type[i];
		ref_type_count[type]++;
	}
    const bool exclude_H = ref_type_count[1] == 0;

    md_index_data_t structures = get_structures(mol, level, exclude_H, temp_alloc);
    const int64_t num_structures = md_index_data_count(structures);

    if (num_structures == 0) {
        MD_LOG_ERROR("Molecule does not have any structures");
        return result;
    }

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

    for (int64_t i = 0; i < num_structures; ++i) {
        md_range_t s_range = structures.ranges[i];
        int32_t s_size = s_range.end - s_range.beg;
        const int* s_indices = structures.indices + s_range.beg;
        if (s_size < ref_graph.vertex_count) continue;

        int s_type_count[256] = {0};
        for (int64_t j = 0; j < s_size; ++j) {
            int idx = s_indices[j];
            uint8_t type = mol->atom.element[idx];
            s_type_count[type]++;
        }

        // Sanity check
        for (int j = 0; j < Num_Elements; ++j) {
            if (ref_type_count[j] > s_type_count[j]) goto next;
        }
        
        graph_t s_graph = make_graph(mol, mol->atom.element, s_indices, s_size, temp_alloc);

        md_index_data_t mappings = find_isomorphisms(&ref_graph, &s_graph, mode, starting_type, temp_alloc);
        // Remap indices to global indices in result
        for (int64_t j = 0; j < md_index_data_count(mappings); ++j) {
            int* beg = md_index_range_beg(mappings, j);
            int* end = md_index_range_end(mappings, j);
            int r_beg = (int)md_array_size(result.indices);
            int r_end = r_beg + (int)(end - beg);
            md_range_t range = {r_beg, r_end};
            md_array_push(result.ranges, range, alloc);

            for (int* it = beg; it != end; ++it) {
                int idx = s_indices[*it];
                md_array_push(result.indices, idx, alloc);
            }
        }
next:;
    }

    md_arena_allocator_destroy(temp_alloc);

    return result;
}

#ifdef __cplusplus
}
#endif
