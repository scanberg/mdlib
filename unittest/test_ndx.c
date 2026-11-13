#include "utest.h"
#include <md_io.h>

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
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &state, MD_STR_LIT(TPR_DIR "nb_fswitch_rf.tpr"), NULL));
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

#include <md_gro.h>
#include <md_edr.h>

static uint64_t default_count(const md_system_t* sys, const char* name, md_allocator_i* alloc) {
    md_bitset_t s;
    if (!md_ndx_default_group(&s, sys, md_str_from_cstr(name), alloc)) return UINT64_MAX;
    return md_bitset_count(s);
}

// GROMACS' default groups of inside-md-pullout: a protein of 1189 atoms (119 CA, 275 hydrogens), 49302 of water, the
// ligand AIN (20), a calcium (CA) and three sodium (NA)
UTEST(ndx, default_groups) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.gro"), NULL));
    ASSERT_EQ(50515u, sys.atom.count);

    EXPECT_EQ(50515u, default_count(&sys, "System", alloc));
    EXPECT_EQ(50515u, default_count(&sys, "system", alloc));         // As grompp matches, whatever the case
    EXPECT_EQ(1189u, default_count(&sys, "Protein", alloc));
    EXPECT_EQ(1189u - 275u, default_count(&sys, "Protein-H", alloc));
    EXPECT_EQ(119u, default_count(&sys, "C-alpha", alloc));
    EXPECT_EQ(119u * 3u, default_count(&sys, "Backbone", alloc));
    EXPECT_EQ(50515u - 1189u, default_count(&sys, "non-Protein", alloc));
    EXPECT_EQ(50515u - 1189u, default_count(&sys, "Non-Protein", alloc));
    EXPECT_EQ(49302u, default_count(&sys, "Water", alloc));
    EXPECT_EQ(49302u, default_count(&sys, "SOL", alloc));
    EXPECT_EQ(50515u - 49302u, default_count(&sys, "non-Water", alloc));
    EXPECT_EQ(4u, default_count(&sys, "Ion", alloc));
    EXPECT_EQ(49302u + 4u, default_count(&sys, "Water_and_ions", alloc));
    EXPECT_EQ(20u, default_count(&sys, "Other", alloc));
    EXPECT_EQ(20u, default_count(&sys, "AIN", alloc));                // A group per residue name of what is not protein or water
    EXPECT_EQ(1u, default_count(&sys, "CA", alloc));                  // The calcium, not the C-alphas
    EXPECT_EQ(3u, default_count(&sys, "NA", alloc));
    EXPECT_EQ(UINT64_MAX, default_count(&sys, "DNA", alloc));          // None of it, so no group
    EXPECT_EQ(UINT64_MAX, default_count(&sys, "ALA", alloc));          // Protein residues get no group of their name
    EXPECT_EQ(UINT64_MAX, default_count(&sys, "Lipids", alloc));
    // The main chain and side chain split the protein's heavy atoms and hydrogens between them
    EXPECT_EQ(1189u, default_count(&sys, "MainChain+H", alloc) + default_count(&sys, "SideChain", alloc));
    md_arena_allocator_destroy(alloc);
}

// Only the structure and the energy file: the energy groups CA Protein AIN SOL rest and the T-coupling group System are
// GROMACS' own of the topology, so they are found without the run's tpr or index file; rest is the three NA
UTEST(ndx, energy_file_names_default_groups) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.gro"), NULL));
    ASSERT_EQ(0u, sys.group.count);
    ASSERT_TRUE(md_edr_system_supplement_from_file(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"), &(md_load_desc_t){ .run = MD_STR_LIT("run/md") }));

    const md_urange_t er = md_system_group_role_range(&sys.group, MD_GROUP_ROLE_ENERGY);
    ASSERT_EQ(5u, er.end - er.beg);
    const char* order[] = { "CA", "Protein", "AIN", "SOL", "rest" };
    const uint64_t counts[] = { 1, 1189, 20, 49302, 3 };
    for (uint32_t i = 0; i < 5; ++i) {
        const uint32_t g = sys.group.role_group[er.beg + i];
        EXPECT_TRUE(md_str_eq(md_str_from_cstr(order[i]), sys.group.name[g]));
        EXPECT_EQ(counts[i], md_bitset_count(sys.group.atoms[g]));
    }
    const md_urange_t tr = md_system_group_role_range(&sys.group, MD_GROUP_ROLE_TC);
    ASSERT_EQ(1u, tr.end - tr.beg);
    const uint32_t system = sys.group.role_group[tr.beg];
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("System"), sys.group.name[system]));
    EXPECT_EQ(50515u, md_bitset_count(sys.group.atoms[system]));

    md_edr_term_groups_t tg;
    ASSERT_TRUE(md_edr_term_groups(&tg, &sys, MD_STR_LIT("T-System")));
    EXPECT_EQ(system, tg.group[0]);
    ASSERT_TRUE(md_edr_term_groups(&tg, &sys, MD_STR_LIT("Lamb-System")));
    ASSERT_TRUE(md_edr_term_groups(&tg, &sys, MD_STR_LIT("Coul-SR:AIN-rest")));
    EXPECT_EQ(3u, md_bitset_count(sys.group.atoms[tg.group[1]]));
    EXPECT_TRUE(md_edr_term_is_of_groups(MD_STR_LIT("Coul-SR:AIN-rest")));
    EXPECT_TRUE(md_edr_term_is_of_groups(MD_STR_LIT("T-System")));
    EXPECT_FALSE(md_edr_term_is_of_groups(MD_STR_LIT("Potential")));
    EXPECT_FALSE(md_edr_term_is_of_groups(MD_STR_LIT("Pres-DC")));
    md_arena_allocator_destroy(alloc);
}
