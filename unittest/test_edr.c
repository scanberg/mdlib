#include "utest.h"

#include <md_edr.h>
#include <core/md_allocator.h>

UTEST(edr, read_file) {
	md_edr_energies_t energies = {0};
	bool result = md_edr_energies_read_file(&energies, STR(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"), default_allocator);
	EXPECT_TRUE(result);
	md_edr_energies_free(&energies);
}
