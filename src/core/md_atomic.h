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
    MD_Z_HE = 2,   // Helium
    MD_Z_LI = 3,   // Lithium
    MD_Z_BE = 4,   // Beryllium
    MD_Z_B  = 5,   // Boron
    MD_Z_C  = 6,   // Carbon
    MD_Z_N  = 7,   // Nitrogen
    MD_Z_O  = 8,   // Oxygen
    MD_Z_F  = 9,   // Fluorine
    MD_Z_NE = 10,  // Neon
    MD_Z_NA = 11,  // Sodium
    MD_Z_MG = 12,  // Magnesium
    MD_Z_AL = 13,  // Aluminium
    MD_Z_SI = 14,  // Silicon
    MD_Z_P  = 15,  // Phosphorus
    MD_Z_S  = 16,  // Sulfur
    MD_Z_CL = 17,  // Chlorine
    MD_Z_AR = 18,  // Argon
    MD_Z_K  = 19,  // Potassium
    MD_Z_CA = 20,  // Calcium
    MD_Z_SC = 21,  // Scandium
    MD_Z_TI = 22,  // Titanium
    MD_Z_V  = 23,  // Vanadium
    MD_Z_CR = 24,  // Chromium
    MD_Z_MN = 25,  // Manganese
    MD_Z_FE = 26,  // Iron
    MD_Z_CO = 27,  // Cobalt
    MD_Z_NI = 28,  // Nickel
    MD_Z_CU = 29,  // Copper
    MD_Z_ZN = 30,  // Zinc
    MD_Z_GA = 31,  // Gallium
    MD_Z_GE = 32,  // Germanium
    MD_Z_AS = 33,  // Arsenic
    MD_Z_SE = 34,  // Selenium
    MD_Z_BR = 35,  // Bromine
    MD_Z_KR = 36,  // Krypton
    MD_Z_RB = 37,  // Rubidium
    MD_Z_SR = 38,  // Strontium
    MD_Z_Y  = 39,  // Yttrium
    MD_Z_ZR = 40,  // Zirconium
    MD_Z_NB = 41,  // Niobium
    MD_Z_MO = 42,  // Molybdenum
    MD_Z_TC = 43,  // Technetium
    MD_Z_RU = 44,  // Ruthenium
    MD_Z_RH = 45,  // Rhodium
    MD_Z_PD = 46,  // Palladium
    MD_Z_AG = 47,  // Silver
    MD_Z_CD = 48,  // Cadmium
    MD_Z_IN = 49,  // Indium
    MD_Z_SN = 50,  // Tin
    MD_Z_SB = 51,  // Antimony
    MD_Z_TE = 52,  // Tellurium
    MD_Z_I  = 53,  // Iodine
    MD_Z_XE = 54,  // Xenon
    MD_Z_CS = 55,  // Caesium
    MD_Z_BA = 56,  // Barium
    MD_Z_LA = 57,  // Lanthanum
    MD_Z_CE = 58,  // Cerium
    MD_Z_PR = 59,  // Praseodymium
    MD_Z_ND = 60,  // Neodymium
    MD_Z_PM = 61,  // Promethium
    MD_Z_SM = 62,  // Samarium
    MD_Z_EU = 63,  // Europium
    MD_Z_GD = 64,  // Gadolinium
    MD_Z_TB = 65,  // Terbium
    MD_Z_DY = 66,  // Dysprosium
    MD_Z_HO = 67,  // Holmium
    MD_Z_ER = 68,  // Erbium
    MD_Z_TM = 69,  // Thulium
    MD_Z_YB = 70,  // Ytterbium
    MD_Z_LU = 71,  // Lutetium
    MD_Z_HF = 72,  // Hafnium
    MD_Z_TA = 73,  // Tantalum
    MD_Z_W  = 74,  // Tungsten
    MD_Z_RE = 75,  // Rhenium
    MD_Z_OS = 76,  // Osmium
    MD_Z_IR = 77,  // Iridium
    MD_Z_PT = 78,  // Platinum
    MD_Z_AU = 79,  // Gold
    MD_Z_HG = 80,  // Mercury
    MD_Z_TL = 81,  // Thallium
    MD_Z_PB = 82,  // Lead
    MD_Z_BI = 83,  // Bismuth
    MD_Z_PO = 84,  // Polonium
    MD_Z_AT = 85,  // Astatine
    MD_Z_RN = 86,  // Radon
    MD_Z_FR = 87,  // Francium
    MD_Z_RA = 88,  // Radium
    MD_Z_AC = 89,  // Actinium
    MD_Z_TH = 90,  // Thorium
    MD_Z_PA = 91,  // Protactinium
    MD_Z_U  = 92,  // Uranium
    MD_Z_NP = 93,  // Neptunium
    MD_Z_PU = 94,  // Plutonium
    MD_Z_AM = 95,  // Americium
    MD_Z_CM = 96,  // Curium
    MD_Z_BK = 97,  // Berkelium
    MD_Z_CF = 98,  // Californium
    MD_Z_ES = 99,  // Einsteinium
    MD_Z_FM = 100, // Fermium
    MD_Z_MD = 101, // Mendelevium
    MD_Z_NO = 102, // Nobelium
    MD_Z_LR = 103, // Lawrencium
    MD_Z_RF = 104, // Rutherfordium
    MD_Z_DB = 105, // Dubnium
    MD_Z_SG = 106, // Seaborgium
    MD_Z_BH = 107, // Bohrium
    MD_Z_HS = 108, // Hassium
    MD_Z_MT = 109, // Meitnerium
    MD_Z_DS = 110, // Darmstadtium
    MD_Z_RG = 111, // Roentgenium
    MD_Z_CN = 112, // Copernicium
    MD_Z_NH = 113, // Nihonium
    MD_Z_FL = 114, // Flerovium
    MD_Z_MC = 115, // Moscovium
    MD_Z_LV = 116, // Livermorium
    MD_Z_TS = 117, // Tennessine
    MD_Z_OG = 118, // Oganesson
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