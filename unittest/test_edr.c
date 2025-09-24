#include "utest.h"

#include <md_edr.h>
#include <core/md_allocator.h>

UTEST(edr, pullout) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, orires) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/orires.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, ener) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/ener.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, dhdl) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR_LIT(MD_UNITTEST_DATA_DIR "/dhdl.edr"), md_get_heap_allocator());
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, comprehensive_validation) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test all available EDR files systematically
    str_t edr_files[] = {
        STR_LIT(MD_UNITTEST_DATA_DIR "/ener.edr"),
        STR_LIT(MD_UNITTEST_DATA_DIR "/dhdl.edr"),
        STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"),
        STR_LIT(MD_UNITTEST_DATA_DIR "/orires.edr")
    };
    
    for (int f = 0; f < 4; ++f) {
        md_edr_energies_t energies = {0};
        bool result = md_edr_energies_parse_file(&energies, edr_files[f], alloc);
        ASSERT_TRUE(result);
        
        EXPECT_GT(energies.num_energies, 0);
        EXPECT_GT(energies.num_frames, 0);
        
        // Validate that frame times are reasonable
        for (int64_t i = 0; i < MIN(100, energies.num_frames); ++i) {
            EXPECT_FALSE(isnan(energies.frame_time[i]));
            EXPECT_FALSE(isinf(energies.frame_time[i]));
            EXPECT_GE(energies.frame_time[i], 0.0);
        }
        
        // Validate that energy data is accessible
        for (int64_t j = 0; j < energies.num_energies; ++j) {
            if (energies.energy[j].values) {
                for (int64_t i = 0; i < MIN(10, energies.num_frames); ++i) {
                    EXPECT_FALSE(isnan(energies.energy[j].values[i]));
                    EXPECT_FALSE(isinf(energies.energy[j].values[i]));
                }
            }
        }
        
        md_edr_energies_free(&energies);
    }
}

UTEST(edr, nonexistent_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.edr");
    
    md_edr_energies_t energies = {0};
    bool result = md_edr_energies_parse_file(&energies, path, alloc);
    EXPECT_FALSE(result);
    
    md_edr_energies_free(&energies);
}
