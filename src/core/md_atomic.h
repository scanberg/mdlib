#pragma once

#include <core/md_str.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Atomic number type: 0 = unknown, 1-118 = element atomic numbers
typedef uint8_t md_atomic_number_t;

// Legacy alias for compatibility
typedef md_atomic_number_t md_element_t;

// Forward declaration of molecule type
struct md_molecule_t;

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

// New preferred API names

// Element symbol and name lookup functions
md_atomic_number_t md_atomic_number_from_symbol(str_t sym);
md_atomic_number_t md_atomic_number_from_symbol_icase(str_t sym);
str_t md_symbol_from_atomic_number(md_atomic_number_t z);
str_t md_name_from_atomic_number(md_atomic_number_t z);

// Element property functions
float md_atomic_mass(md_atomic_number_t z);
float md_vdw_radius(md_atomic_number_t z);
float md_covalent_radius(md_atomic_number_t z);
int   md_max_valence(md_atomic_number_t z);
uint32_t md_cpk_color(md_atomic_number_t z);

// Per-atom inference from labels (atom name + residue)
md_atomic_number_t md_atom_infer_atomic_number(str_t atom_name, str_t res_name);

// Batch form wired to molecule structure
bool md_atoms_infer_atomic_numbers(md_atomic_number_t out[], size_t n, const struct md_molecule_t* mol);

#ifdef __cplusplus
}
#endif