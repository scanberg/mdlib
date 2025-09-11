// Test to verify API equivalence between old and new systems
#include "utest.h"
#include <core/md_atomic.h>
#include <md_util.h>

UTEST(api_equivalence, symbol_lookup_consistency) {
    // Test that old and new APIs return the same results
    for (int i = 1; i <= 118; ++i) {
        str_t symbol_old = md_util_element_symbol(i);
        str_t symbol_new = md_symbol_from_atomic_number(i);
        
        EXPECT_TRUE(str_eq(symbol_old, symbol_new));
        
        // Test reverse lookup
        md_atomic_number_t old_lookup = md_util_element_lookup(symbol_old);
        md_atomic_number_t new_lookup = md_atomic_number_from_symbol(symbol_old);
        
        EXPECT_EQ(old_lookup, new_lookup);
        EXPECT_EQ(old_lookup, i);
    }
}

UTEST(api_equivalence, property_consistency) {
    // Test a few key elements for property consistency
    md_atomic_number_t test_elements[] = {MD_Z_H, MD_Z_C, MD_Z_N, MD_Z_O, MD_Z_CA, MD_Z_FE};
    
    for (size_t i = 0; i < ARRAY_SIZE(test_elements); ++i) {
        md_atomic_number_t z = test_elements[i];
        
        // Test masses
        float mass_old = md_util_element_atomic_mass(z);
        float mass_new = md_atomic_mass(z);
        EXPECT_EQ(mass_old, mass_new);
        
        // Test radii
        float vdw_old = md_util_element_vdw_radius(z);
        float vdw_new = md_vdw_radius(z);
        EXPECT_EQ(vdw_old, vdw_new);
        
        float cov_old = md_util_element_covalent_radius(z);
        float cov_new = md_covalent_radius(z);
        EXPECT_EQ(cov_old, cov_new);
        
        // Test valence and color
        int val_old = md_util_element_max_valence(z);
        int val_new = md_max_valence(z);
        EXPECT_EQ(val_old, val_new);
        
        uint32_t color_old = md_util_element_cpk_color(z);
        uint32_t color_new = md_cpk_color(z);
        EXPECT_EQ(color_old, color_new);
    }
}