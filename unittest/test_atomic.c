#include "utest.h"
#include <md_util.h>
#include <md_types.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

UTEST(atomic, enum_constants) {
    // Test that enum constants are correct
    EXPECT_EQ(MD_Z_X, 0);   // Unknown
    EXPECT_EQ(MD_Z_H, 1);   // Hydrogen
    EXPECT_EQ(MD_Z_He, 2);  // Helium
    EXPECT_EQ(MD_Z_C, 6);   // Carbon
    EXPECT_EQ(MD_Z_N, 7);   // Nitrogen
    EXPECT_EQ(MD_Z_O, 8);   // Oxygen
    EXPECT_EQ(MD_Z_P, 15);  // Phosphorus
    EXPECT_EQ(MD_Z_S, 16);  // Sulfur
    EXPECT_EQ(MD_Z_Ca, 20); // Calcium
    EXPECT_EQ(MD_Z_Cl, 17); // Chlorine
    EXPECT_EQ(MD_Z_Br, 35); // Bromine
    EXPECT_EQ(MD_Z_Na, 11); // Sodium
    EXPECT_EQ(MD_Z_Fe, 26); // Iron
    EXPECT_EQ(MD_Z_Og, 118); // Oganesson
}

UTEST(atomic, symbol_lookup) {
    // Test symbol lookup functions
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("H"), false), MD_Z_H);
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("C"), false), MD_Z_C);
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("He"), false), MD_Z_He);
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("Ca"), false), MD_Z_Ca);
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("Unknown"), false), MD_Z_X);
    
    // Test case insensitive lookup
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("h"), true), MD_Z_H);
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("ca"), true), MD_Z_Ca);
    EXPECT_EQ(md_atomic_number_from_symbol(STR_LIT("HE"), true), MD_Z_He);
}

UTEST(atomic, symbol_from_number) {
    // Test reverse lookup
    str_t h_symbol = md_atomic_number_symbol(MD_Z_H);
    EXPECT_TRUE(str_eq_cstr(h_symbol, "H"));
    
    str_t c_symbol = md_atomic_number_symbol(MD_Z_C);
    EXPECT_TRUE(str_eq_cstr(c_symbol, "C"));
    
    str_t ca_symbol = md_atomic_number_symbol(MD_Z_Ca);
    EXPECT_TRUE(str_eq_cstr(ca_symbol, "Ca"));
}

UTEST(atomic, element_properties) {
    // Test that we can get basic properties
    EXPECT_GT(md_atomic_number_mass(MD_Z_H), 0.0f);
    EXPECT_GT(md_atomic_number_mass(MD_Z_C), 0.0f);
    EXPECT_GT(md_atomic_number_vdw_radius(MD_Z_H), 0.0f);
    EXPECT_GT(md_atomic_number_covalent_radius(MD_Z_C), 0.0f);
    EXPECT_GT(md_atomic_number_max_valence(MD_Z_C), 0);
    EXPECT_GT(md_atomic_number_cpk_color(MD_Z_C), 0);
}

UTEST(atomic, inference_water) {
    // Test water atom inference
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("O"), STR_LIT("HOH"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("OW"), STR_LIT("HOH"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("OH2"), STR_LIT("WAT"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("H"), STR_LIT("HOH"), 0), MD_Z_H);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("H1"), STR_LIT("TIP3"), 0), MD_Z_H);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("HW"), STR_LIT("SPC"), 0), MD_Z_H);
}

UTEST(atomic, inference_amino_acid) {
    // Test amino acid atom inference
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("CA"), STR_LIT("ALA"), 0), MD_Z_C);  // Alpha carbon, not calcium
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("N"), STR_LIT("GLY"), 0), MD_Z_N);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("C"), STR_LIT("SER"), 0), MD_Z_C);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("O"), STR_LIT("TRP"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("OXT"), STR_LIT("PHE"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("OG"), STR_LIT("SER"), 0), MD_Z_O);  // Serine hydroxyl
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("SG"), STR_LIT("CYS"), 0), MD_Z_S);  // Cysteine sulfur
}

UTEST(atomic, inference_nucleic_acid) {
    // Test nucleic acid atom inference
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("P"), STR_LIT("DA"), 0), MD_Z_P);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("OP1"), STR_LIT("DG"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("O2P"), STR_LIT("A"), 0), MD_Z_O);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("P"), STR_LIT("U"), 0), MD_Z_P);
}

UTEST(atomic, inference_ions) {
    // Test ion inference (residue name is element)
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("NA"), STR_LIT("NA"), 0), MD_Z_Na);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT(""), STR_LIT("K"), 0), MD_Z_K);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("MG"), STR_LIT("MG"), 0), MD_Z_Mg);
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("CL"), STR_LIT("CL"), 0), MD_Z_Cl);
}

UTEST(atomic, inference_fallbacks) {
    // Test fallback mechanisms
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("CL12"), STR_LIT(""), 0), MD_Z_Cl);  // Two-letter heuristic
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("BR1"), STR_LIT(""), 0), MD_Z_Br);   // Two-letter heuristic
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("H123"), STR_LIT(""), 0), MD_Z_H);   // First letter fallback
    EXPECT_EQ(md_atomic_number_infer_from_label(STR_LIT("C99"), STR_LIT(""), 0), MD_Z_C);    // First letter fallback
}

UTEST(atomic, backward_compatibility) {
    // Test that old API still works through wrappers
    EXPECT_EQ(md_util_element_lookup(STR_LIT("H"), true), MD_Z_H);
    EXPECT_EQ(md_util_element_lookup(STR_LIT("ca"), true), MD_Z_Ca);
    
    str_t symbol = md_util_element_symbol(MD_Z_C);
    EXPECT_TRUE(str_eq_cstr(symbol, "C"));
    
    str_t name = md_util_element_name(MD_Z_O);
    EXPECT_TRUE(str_eq_cstr(name, "Oxygen"));
    
    EXPECT_GT(md_util_element_atomic_mass(MD_Z_C), 0.0f);
    EXPECT_GT(md_util_element_vdw_radius(MD_Z_H), 0.0f);
    EXPECT_GT(md_util_element_covalent_radius(MD_Z_N), 0.0f);
    EXPECT_GT(md_util_element_max_valence(MD_Z_C), 0);
    EXPECT_GT(md_util_element_cpk_color(MD_Z_O), 0);
}