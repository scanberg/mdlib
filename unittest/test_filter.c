#include "utest.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <md_molecule.h>
#include <md_gro.h>
#include <md_util.h>

#include <md_filter.h>

#define TEST(str) md_filter(&bf, STR(str), &mol, NULL, &is_dynamic, err, sizeof(err))

UTEST(filter, centered) {
    md_molecule_t mol = {0};
    const str_t gro_file = STR(MD_UNITTEST_DATA_DIR "/centered.gro");

    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, gro_file, md_heap_allocator));
    ASSERT_TRUE(md_util_postprocess_molecule(&mol, md_heap_allocator, MD_UTIL_POSTPROCESS_ALL));
    
    md_bitfield_t bf = md_bitfield_create(md_heap_allocator);
    char err[256];
    bool is_dynamic = false;
    
    EXPECT_TRUE(TEST("resname('ALA')"));
    EXPECT_TRUE(TEST("within(10, residue(1))"));

    md_bitfield_free(&bf);
    md_molecule_free(&mol, md_heap_allocator);
}
