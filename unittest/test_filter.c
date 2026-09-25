#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <md_system.h>
#include <md_gro.h>
#include <md_util.h>

#include <md_filter.h>

#define TEST(str) md_filter(&bf, STR_LIT(str), &sys, &state, NULL, &is_dynamic, err, sizeof(err))

UTEST(filter, centered) {
    const str_t gro_file = STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    
    md_bitfield_t bf = md_bitfield_create(alloc);
    char err[256];
    bool is_dynamic = false;
    
    EXPECT_TRUE(TEST("resname('ALA')"));
    EXPECT_TRUE(TEST("within(10, residue(1))"));

    md_arena_allocator_destroy(alloc);
}

// A filter whose number of bitfields is only known once evaluated: the residues within some distance
UTEST(filter, variable_length) {
    const str_t gro_file = STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    char err[256];
    bool is_dynamic = false;

    md_bitfield_t near = md_bitfield_create(alloc);
    ASSERT_TRUE(md_filter(&near, STR_LIT("within(5, residue(1))"), &sys, &state, NULL, &is_dynamic, err, sizeof(err)));

    // Expected: every residue with an atom within the distance, in full
    md_bitfield_t expected = md_bitfield_create(alloc);
    size_t num_residues = 0;
    for (size_t i = 0; i < sys.component.count; ++i) {
        const md_urange_t range = md_component_atom_range(&sys.component, i);
        if (md_bitfield_popcount_range(&near, range.beg, range.end) > 0) {
            md_bitfield_set_range(&expected, range.beg, range.end);
            num_residues += 1;
        }
    }
    ASSERT_GT(num_residues, (size_t)1);

    md_bitfield_t bf = md_bitfield_create(alloc);
    is_dynamic = false;
    EXPECT_TRUE(md_filter(&bf, STR_LIT("residue(within(5, residue(1)))"), &sys, &state, NULL, &is_dynamic, err, sizeof(err)));
    EXPECT_TRUE(is_dynamic);
    EXPECT_EQ(md_bitfield_popcount(&expected), md_bitfield_popcount(&bf));
    md_bitfield_t both = md_bitfield_create(alloc);
    md_bitfield_and(&both, &bf, &expected);
    EXPECT_EQ(md_bitfield_popcount(&expected), md_bitfield_popcount(&both));

    // md_filter_evaluate keeps the residues apart
    md_array(md_bitfield_t) arr = 0;
    EXPECT_TRUE(md_filter_evaluate(&arr, STR_LIT("residue(within(5, residue(1)))"), &sys, &state, NULL, &is_dynamic, err, sizeof(err), alloc));
    EXPECT_EQ(num_residues, md_array_size(arr));

    md_arena_allocator_destroy(alloc);
}
