#include "utest.h"

#include <md_ndx.h>
#include <md_tpr.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_str.h>

#define TPR_DIR MD_UNITTEST_DATA_DIR "/tpr/"

UTEST(ndx, parse) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    const char* text =
        "; written by hand\n"
        "[ System ]\n"
        "   1    2    3    4    5\r\n"
        "   6    7    8\n"
        "\n"
        "[ Protein-H ]   ; a comment\n"
        "1 3 5\n"
        "[ empty ]\n"
        "[ Water_and_ions ]\n"
        "8 6 7 6\n";
    md_ndx_data_t data = {0};
    ASSERT_TRUE(md_ndx_data_parse_str(&data, md_str_from_cstr(text), arena));
    ASSERT_EQ(4u, data.num_groups);
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("System"), data.groups[0].name));
    EXPECT_EQ(8u, data.groups[0].count);
    EXPECT_EQ(7u, data.groups[0].atom_idx[7]);
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("Protein-H"), data.groups[1].name));
    EXPECT_EQ(3u, data.groups[1].count);
    EXPECT_EQ(4u, data.groups[1].atom_idx[2]);
    EXPECT_EQ(0u, data.groups[2].count);
    EXPECT_EQ(4u, data.groups[3].count);
    EXPECT_EQ(8u, data.max_atom);
    md_ndx_data_free(&data, arena);

    EXPECT_FALSE(md_ndx_data_parse_str(&data, MD_STR_LIT("1 2 3\n[ A ]\n4\n"), arena));     // Atoms before a header
    EXPECT_FALSE(md_ndx_data_parse_str(&data, MD_STR_LIT("[ A ]\n1 0 2\n"), arena));        // 1 based
    EXPECT_FALSE(md_ndx_data_parse_str(&data, MD_STR_LIT("[ A ]\n1 x\n"), arena));
    EXPECT_FALSE(md_ndx_data_parse_str(&data, MD_STR_LIT("[ A \n1\n"), arena));
    EXPECT_FALSE(md_ndx_data_parse_str(&data, MD_STR_LIT("; nothing\n"), arena));
    md_arena_allocator_destroy(arena);
}

// The index file grompp was given for the nb tprs: the same groups as the tpr keeps, which they merge with
UTEST(ndx, supplement_tpr_system) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &state, MD_STR_LIT(TPR_DIR "nb_fswitch_rf.tpr")));
    ASSERT_EQ(2u, sys.group.count);
    const uint64_t version = sys.topology_version;

    md_ndx_data_t data = {0};
    ASSERT_TRUE(md_ndx_data_parse_file(&data, MD_STR_LIT(TPR_DIR "nb_index.ndx"), arena));
    ASSERT_EQ(2u, data.num_groups);
    ASSERT_TRUE(md_ndx_system_supplement(&sys, &data));
    EXPECT_EQ(2u, sys.group.count);                                     // Merged: same names, same atoms
    EXPECT_EQ(1u << MD_GROUP_ROLE_ENERGY, sys.group.roles[0]);
    EXPECT_NE(version, sys.topology_version);

    // Another index file of the system adds what is new to it
    md_ndx_data_t more = {0};
    ASSERT_TRUE(md_ndx_data_parse_str(&more, MD_STR_LIT("[ CHN ]\n1 2 3\n[ first_and_last ]\n1 19\n"), arena));
    ASSERT_TRUE(md_ndx_system_supplement(&sys, &more));
    ASSERT_EQ(4u, sys.group.count);
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("CHN"), sys.group.name[2]));       // A name twice, of other atoms
    EXPECT_EQ(0u, sys.group.roles[2]);
    EXPECT_EQ(2u, md_bitset_count(sys.group.atoms[3]));
    EXPECT_TRUE(md_bitset_test(sys.group.atoms[3], 18));

    // Of another system
    md_ndx_data_t other = {0};
    ASSERT_TRUE(md_ndx_data_parse_str(&other, MD_STR_LIT("[ big ]\n1 20\n"), arena));
    EXPECT_FALSE(md_ndx_system_supplement(&sys, &other));
    EXPECT_EQ(4u, sys.group.count);

    md_arena_allocator_destroy(arena);
}
