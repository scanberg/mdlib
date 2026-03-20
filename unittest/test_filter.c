#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <md_system.h>
#include <md_gro.h>
#include <md_util.h>

#include <md_filter.h>

#define TEST(str) md_filter(&bf, STR_LIT(str), &sys, sys.atom.x, sys.atom.y, sys.atom.z, NULL, &is_dynamic, err, sizeof(err))

UTEST(filter, centered) {
    const str_t gro_file = STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, gro_file));
    ASSERT_TRUE(md_util_system_postprocess(&sys, MD_UTIL_POSTPROCESS_ALL));
    
    md_bitfield_t bf = md_bitfield_create(alloc);
    char err[256];
    bool is_dynamic = false;
    
    EXPECT_TRUE(TEST("resname('ALA')"));
    EXPECT_TRUE(TEST("within(10, residue(1))"));

    md_arena_allocator_destroy(alloc);
}
