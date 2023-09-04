#include "utest.h"

#include <md_edr.h>
#include <core/md_allocator.h>

UTEST(edr, pullout) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"), md_heap_allocator);
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, orires) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR(MD_UNITTEST_DATA_DIR "/orires.edr"), md_heap_allocator);
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, ener) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR(MD_UNITTEST_DATA_DIR "/ener.edr"), md_heap_allocator);
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}

UTEST(edr, dhdl) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_parse_file(&energies, STR(MD_UNITTEST_DATA_DIR "/dhdl.edr"), md_heap_allocator);
	EXPECT_TRUE(result);

	//EXPECT_TRUE(md_unit_equal(energies.units[0], unit_ki))
	md_edr_energies_free(&energies);
}