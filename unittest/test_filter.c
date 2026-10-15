#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_bitset.h>
#include <md_system.h>
#include <md_gro.h>
#include <md_util.h>

#include <md_filter.h>

#include <string.h>

// Whether the filter picks atoms, with why not in res.error
#define TEST(str) (res = md_filter(&bf, STR_LIT(str), &fctx, alloc)).success

UTEST(filter, centered) {
    const str_t gro_file = STR_INIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    
    md_bitset_t bf = {0};
    const md_filter_context_t fctx = {&sys, &state, NULL};
    md_filter_result_t res;

    EXPECT_TRUE(TEST("resname('ALA')"));
    EXPECT_TRUE(TEST("within(10, residue(1))"));

    md_arena_allocator_destroy(alloc);
}

// A filter whose number of bitfields is only known once evaluated: the residues within some distance
UTEST(filter, variable_length) {
    const str_t gro_file = STR_INIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    const md_filter_context_t fctx = {&sys, &state, NULL};

    md_bitset_t near = {0};
    ASSERT_TRUE(md_filter(&near, STR_LIT("within(5, residue(1))"), &fctx, alloc).success);

    // Expected: every residue with an atom within the distance, in full
    md_bitmask_t expected_b = md_bitmask_create(0, (uint32_t)sys.atom.count, alloc);
    size_t num_residues = 0;
    for (size_t i = 0; i < sys.component.count; ++i) {
        const md_urange_t range = md_component_atom_range(&sys.component, i);
        if (md_bitset_count_range(near, range.beg, range.end) > 0) {
            md_bitmask_set_range(&expected_b, range.beg, range.end);
            num_residues += 1;
        }
    }
    ASSERT_GT(num_residues, (size_t)1);
    const md_bitset_t expected = md_bitmask_extract(&expected_b, alloc);

    md_bitset_t bf = {0};
    const md_filter_result_t res = md_filter(&bf, STR_LIT("residue(within(5, residue(1)))"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_TRUE(res.is_dynamic);
    EXPECT_EQ(md_bitset_count(expected), md_bitset_count(bf));
    EXPECT_TRUE(md_bitset_equal(expected, bf));

    // md_filter_sets keeps the residues apart
    md_array(md_bitset_t) arr = 0;
    EXPECT_TRUE(md_filter_sets(&arr, STR_LIT("residue(within(5, residue(1)))"), &fctx, alloc).success);
    EXPECT_EQ(num_residues, md_array_size(arr));

    md_arena_allocator_destroy(alloc);
}

// The QM region and its environment are selectable only in a system loaded from a quantum chemistry
// calculation: in any other, both are compile errors rather than empty selections
UTEST(filter, qm_regions_need_a_qm_system) {
    const str_t gro_file = STR_INIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    md_bitset_t bf = {0};
    const md_filter_context_t fctx = {&sys, &state, NULL};
    md_filter_result_t res;

    EXPECT_FALSE(TEST("qm"));
    EXPECT_TRUE(strstr(res.error, "no QM region") != NULL);
    EXPECT_FALSE(TEST("environment"));
    EXPECT_TRUE(strstr(res.error, "no QM region") != NULL);
    EXPECT_FALSE(TEST("within(5, qm)"));

    md_arena_allocator_destroy(alloc);
}

// What goes wrong is said, and where in the expression
UTEST(filter, errors_say_where) {
    const str_t gro_file = STR_INIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    md_bitset_t bf = {0};
    const md_filter_context_t fctx = {&sys, &state, NULL};

    // A name nothing goes by: the error is about it
    const char* expr = "resname('ALA') and nonsense";
    md_filter_result_t res = md_filter(&bf, str_from_cstr(expr), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_NE('\0', res.error[0]);
    EXPECT_EQ(19, res.error_range.beg);
    EXPECT_EQ(27, res.error_range.end);
    EXPECT_EQ(0u, md_bitset_count(bf));

    // Something that does not pick atoms
    res = md_filter(&bf, STR_LIT("distance(1, 2)"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "does not pick atoms") != NULL);

    // Nothing at all picks every atom, as one set
    res = md_filter(&bf, STR_LIT("   "), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_FALSE(res.is_dynamic);
    EXPECT_EQ('\0', res.error[0]);
    EXPECT_EQ(sys.atom.count, md_bitset_count(bf));
    md_array(md_bitset_t) all_sets = 0;
    EXPECT_TRUE(md_filter_sets(&all_sets, STR_LIT(""), &fctx, alloc).success);
    ASSERT_EQ(1, md_array_size(all_sets));
    EXPECT_EQ(sys.atom.count, md_bitset_count(all_sets[0]));

    // Without a state, in the reference one; static or not
    const md_filter_context_t reference = {&sys, NULL, NULL};
    res = md_filter(&bf, STR_LIT("within(10, residue(1))"), &reference, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_TRUE(res.is_dynamic);
    EXPECT_EQ('\0', res.error[0]);
    res = md_filter(&bf, STR_LIT("resname('ALA')"), &reference, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_FALSE(res.is_dynamic);

    md_arena_allocator_destroy(alloc);
}
