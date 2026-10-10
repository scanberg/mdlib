#include "utest.h"

// Units given to structure matching (md_match_desc_t.units), the rings above all, the any bond '~' of a SMILES
// query, and ring() in the script: by pattern, kind and size. The systems are single components of unittest/ccd_chem.txt (all atoms given),
// which need no test data.

#include <md_system.h>
#include <md_chem.h>
#include <md_match.h>
#include <md_smiles.h>
#include <md_filter.h>
#include <md_util.h>

#include "component_util.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_bitset.h>

#include <stdio.h>
#include <string.h>

typedef struct ring_sys_t {
    md_allocator_i*   alloc;
    md_system_t       sys;
    md_system_state_t st;
} ring_sys_t;

struct ring {
    md_allocator_i* alloc;
    md_array(ref_comp_t) comps;
    ring_sys_t trp;     // Indole: a benzene and a pyrrole ring sharing a bond
    ring_sys_t ala;     // No ring, a carboxylic acid
    ring_sys_t pro;     // Pyrrolidine
    ring_sys_t mes;     // Morpholine
    ring_sys_t ade;     // Adenosine monophosphate: purine and ribose
};

static void ring_sys_init(ring_sys_t* r, const ref_comp_t* comps, const char* id) {
    const ref_comp_t* c = find_component(comps, id);
    ASSERT(c);
    r->alloc = md_vm_arena_create(GIGABYTES(1));
    r->sys = (md_system_t){ .alloc = r->alloc };
    r->st  = (md_system_state_t){ .alloc = r->alloc };
    int* map = md_alloc(r->alloc, sizeof(int) * (md_array_size(c->atom) + 1));
    build_component(&r->sys, &r->st, map, c, false, false);
    md_chem_perceive(&r->sys, &r->st, MD_CHEM_FLAG_NONE);
}

UTEST_F_SETUP(ring) {
    utest_fixture->alloc = md_vm_arena_create(GIGABYTES(1));
    utest_fixture->comps = read_components(STR_LIT(MD_UNITTEST_SOURCE_DIR "/ccd_chem.txt"), utest_fixture->alloc);
    ASSERT_GT((int)md_array_size(utest_fixture->comps), 100);
    ring_sys_init(&utest_fixture->trp, utest_fixture->comps, "TRP");
    ring_sys_init(&utest_fixture->ala, utest_fixture->comps, "ALA");
    ring_sys_init(&utest_fixture->pro, utest_fixture->comps, "PRO");
    ring_sys_init(&utest_fixture->mes, utest_fixture->comps, "MES");
    ring_sys_init(&utest_fixture->ade, utest_fixture->comps, "A");
}

UTEST_F_TEARDOWN(ring) {
    md_vm_arena_destroy(utest_fixture->trp.alloc);
    md_vm_arena_destroy(utest_fixture->ala.alloc);
    md_vm_arena_destroy(utest_fixture->pro.alloc);
    md_vm_arena_destroy(utest_fixture->mes.alloc);
    md_vm_arena_destroy(utest_fixture->ade.alloc);
    md_vm_arena_destroy(utest_fixture->alloc);
}

static size_t ring_count_matches(const char* smiles, md_match_level_t level, md_match_mode_t mode, const md_system_t* sys) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_match_query_t query = {0};
    md_smiles_error_t err = {0};
    if (!md_match_query_init_smiles(&query, str_from_cstr(smiles), alloc, &err)) {
        printf("SMILES error in '%s' at %zu: %s\n", smiles, err.offset, err.message);
        return SIZE_MAX;
    }
    const md_match_desc_t desc = { .query = &query, .level = level, .mode = mode };
    size_t count = 0;
    const bool ok = md_match_for_each(&count, &desc, sys, NULL, NULL);
    md_match_query_free(&query);
    return ok ? count : SIZE_MAX;
}

// Matches with the rings of the system as the units
static size_t ring_count_in_rings(const char* smiles, md_match_mode_t mode, const md_system_t* sys) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    const size_t num_rings = md_index_data_num_ranges(&sys->ring);
    md_bitset_t* units = md_alloc(alloc, sizeof(md_bitset_t) * (num_rings + 1));
    for (size_t r = 0; r < num_rings; ++r) {
        units[r] = md_bitset_from_indices((const uint32_t*)md_index_range_beg(&sys->ring, r), md_index_range_size(&sys->ring, r), alloc);
    }
    md_match_query_t query = {0};
    size_t count = SIZE_MAX;
    if (md_match_query_init_smiles(&query, str_from_cstr(smiles), alloc, NULL)) {
        const md_match_desc_t desc = { .query = &query, .mode = mode, .units = units, .num_units = num_rings };
        if (!md_match_for_each(&count, &desc, sys, NULL, NULL)) count = SIZE_MAX;
        md_match_query_free(&query);
    }
    md_arena_allocator_destroy(alloc);
    return count;
}

// The number of selections ring() and the like give, SIZE_MAX on an error
static size_t ring_filter_count(const char* expr, const ring_sys_t* r, size_t* out_atoms) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_array(md_bitfield_t) arr = 0;
    char err[256] = {0};
    bool is_dynamic = false;
    size_t result = SIZE_MAX;
    if (md_filter_evaluate(&arr, str_from_cstr(expr), &r->sys, &r->st, NULL, &is_dynamic, err, sizeof(err), alloc)) {
        result = md_array_size(arr);
        if (out_atoms) {
            md_bitfield_t all = md_bitfield_create(alloc);
            for (size_t i = 0; i < result; ++i) md_bitfield_or_inplace(&all, &arr[i]);
            *out_atoms = md_bitfield_popcount(&all);
        }
    } else {
        printf("'%s': %s\n", expr, err);
    }
    md_arena_allocator_destroy(alloc);
    return result;
}

static bool ring_filter_fails_with(const char* expr, const char* expected, const ring_sys_t* r) {
    md_bitfield_t bf = md_bitfield_create(md_get_heap_allocator());
    char err[256] = {0};
    bool is_dynamic = false;
    const bool ok = md_filter(&bf, str_from_cstr(expr), &r->sys, &r->st, NULL, &is_dynamic, err, sizeof(err));
    md_bitfield_free(&bf);
    if (ok) {
        printf("'%s' was expected to fail\n", expr);
        return false;
    }
    if (!strstr(err, expected)) {
        printf("'%s' failed with '%s', expected '%s'\n", expr, err, expected);
        return false;
    }
    return true;
}

// ### '~' ###

UTEST_F(ring, any_bond) {
    md_smiles_t g = {0};
    md_smiles_error_t err = {0};
    ASSERT_TRUE(md_smiles_parse(&g, STR_LIT("C~N"), md_get_heap_allocator(), &err));
    EXPECT_EQ(MD_SMILES_BOND_ANY, g.bonds[0].order);
    md_smiles_free(&g);
    EXPECT_FALSE(md_smiles_parse(&g, STR_LIT("C~=N"), md_get_heap_allocator(), &err));

    // Alanine with all of its hydrogens: a carboxylic acid, its C=O and its C-OH
    const md_system_t* ala = &utest_fixture->ala.sys;
    EXPECT_EQ(1, (int)ring_count_matches("C=O", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(1, (int)ring_count_matches("C-O", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(1, (int)ring_count_matches("CO",  MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
    EXPECT_EQ(2, (int)ring_count_matches("C~O", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, ala));
}

// ### The rings as units ###

UTEST_F(ring, match_level) {
    const md_system_t* trp = &utest_fixture->trp.sys;
    ASSERT_EQ(2, (int)md_index_data_num_ranges(&trp->ring));

    // The C-C bonds: 6 in the benzene ring and 3 in the pyrrole ring, the bond they share in both. Within the
    // molecule it is one bond.
    EXPECT_EQ(2, (int)ring_count_in_rings("cc", MD_MATCH_MODE_ONE_PER_UNIT, trp));
    EXPECT_EQ(9, (int)ring_count_in_rings("cc", MD_MATCH_MODE_UNIQUE, trp));
    EXPECT_EQ(8, (int)ring_count_matches("cc", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, trp));
    // Three C-C without a shared atom in the benzene ring and two in the pyrrole ring: DISJOINT is per ring as well,
    // the shared bond being free in each
    EXPECT_EQ(5, (int)ring_count_in_rings("cc", MD_MATCH_MODE_DISJOINT, trp));

    // Along the ring, not into what is bonded to it: the CB of the side chain is not in a ring
    EXPECT_EQ(1, (int)ring_count_matches("Cc", MD_MATCH_LEVEL_STRUCTURE, MD_MATCH_MODE_UNIQUE, trp));
    EXPECT_EQ(0, (int)ring_count_in_rings("Cc", MD_MATCH_MODE_UNIQUE, trp));

    // A closed chain is the whole ring
    EXPECT_EQ(1, (int)ring_count_in_rings("c1ccccc1",   MD_MATCH_MODE_UNIQUE, trp));
    EXPECT_EQ(1, (int)ring_count_in_rings("c1cc[nH]c1", MD_MATCH_MODE_UNIQUE, trp));
    EXPECT_EQ(0, (int)ring_count_in_rings("c1ccccc1c",  MD_MATCH_MODE_UNIQUE, trp));
    // An atom has two neighbors along a ring: c(c)c is a chain of three, c(c)(c)c a branch which never matches
    EXPECT_EQ(8, (int)ring_count_in_rings("c(c)c",      MD_MATCH_MODE_UNIQUE, trp));
    EXPECT_EQ(0, (int)ring_count_in_rings("c(c)(c)c",   MD_MATCH_MODE_UNIQUE, trp));
    // Parts anywhere in one ring, on atoms of their own
    EXPECT_EQ(1, (int)ring_count_in_rings("n.c.c.c.c", MD_MATCH_MODE_ONE_PER_UNIT, trp));
    EXPECT_EQ(0, (int)ring_count_in_rings("n.n",       MD_MATCH_MODE_ONE_PER_UNIT, trp));

    // The unit is the index of the ring
    md_allocator_i* alloc = md_get_heap_allocator();
    md_bitset_t units[3];
    for (size_t r = 0; r < 2; ++r) {
        units[r] = md_bitset_from_indices((const uint32_t*)md_index_range_beg(&trp->ring, r), md_index_range_size(&trp->ring, r), alloc);
    }
    md_match_query_t query = {0};
    ASSERT_TRUE(md_match_query_init_smiles(&query, STR_LIT("[nH]"), alloc, NULL));
    md_match_desc_t desc = { .query = &query, .units = units, .num_units = 2 };
    md_match_result_t res = {0};
    ASSERT_TRUE(md_match_find(&res, &desc, trp, md_get_heap_allocator()));
    ASSERT_EQ(1, (int)res.count);
    EXPECT_EQ(5, (int)md_index_range_size(&trp->ring, res.unit[0]));
    md_match_result_free(&res);

    // Any sets of atoms: the whole molecule as a third unit, overlapping both rings, finds the N-H again
    units[2] = md_bitset_range(0, (uint32_t)trp->atom.count);
    desc.num_units = 3;
    size_t count = 0;
    ASSERT_TRUE(md_match_for_each(&count, &desc, trp, NULL, NULL));
    EXPECT_EQ(2, (int)count);
    // A unit beyond the atoms of the system is an error
    units[2] = md_bitset_range(0, (uint32_t)trp->atom.count + 1);
    EXPECT_FALSE(md_match_for_each(&count, &desc, trp, NULL, NULL));

    md_match_query_free(&query);
    for (size_t r = 0; r < 2; ++r) md_bitset_free(&units[r], alloc);

    // A system without rings has no units
    EXPECT_EQ(0, (int)ring_count_in_rings("C", MD_MATCH_MODE_UNIQUE, &utest_fixture->ala.sys));
}

// ### ring() ###

UTEST_F(ring, script_ring) {
    const ring_sys_t* trp = &utest_fixture->trp;
    const ring_sys_t* ade = &utest_fixture->ade;
    size_t atoms = 0;

    // All, and by size
    EXPECT_EQ(2, (int)ring_filter_count("ring()", trp, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring(size=5)", trp, NULL));
    EXPECT_EQ(2, (int)ring_filter_count("ring(size=5:6)", trp, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring(size=6:)", trp, NULL));
    EXPECT_EQ(0, (int)ring_filter_count("ring(size=7)", trp, NULL));

    // By kind: a list is any of them
    EXPECT_EQ(2, (int)ring_filter_count("ring(kind='aromatic')", trp, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring(kind='heteroaromatic')", trp, NULL));
    EXPECT_EQ(0, (int)ring_filter_count("ring(kind='aliphatic')", trp, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring(kind='sugar')", ade, NULL));
    EXPECT_EQ(2, (int)ring_filter_count("ring(kind='heteroaromatic')", ade, NULL));
    EXPECT_EQ(3, (int)ring_filter_count("ring(kind={'sugar', 'heteroaromatic'})", ade, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring(kind='aliphatic')", ade, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring(kind='hetero')", &utest_fixture->pro, NULL));

    // By pattern: along the ring, in either direction, from any atom
    EXPECT_EQ(1, (int)ring_filter_count("ring('c[nH]c')", trp, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring('c1ccccc1')", trp, &atoms));
    EXPECT_EQ(6, (int)atoms);
    EXPECT_EQ(2, (int)ring_filter_count("ring('cc')", trp, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring('CCNCC')", &utest_fixture->pro, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring('CNCCC')", &utest_fixture->pro, NULL));
    EXPECT_EQ(0, (int)ring_filter_count("ring('CNNCC')", &utest_fixture->pro, NULL));
    // Anywhere in the ring: '.'
    EXPECT_EQ(1, (int)ring_filter_count("ring('N.O')", &utest_fixture->mes, NULL));
    EXPECT_EQ(0, (int)ring_filter_count("ring('N.O')", trp, NULL));
    // A list of patterns is any of them
    EXPECT_EQ(2, (int)ring_filter_count("ring({'c1ccccc1', 'c[nH]c'})", trp, NULL));

    // The filters together: all of them
    EXPECT_EQ(1, (int)ring_filter_count("ring('cc', kind='aromatic', size=6)", trp, NULL));
    EXPECT_EQ(2, (int)ring_filter_count("ring('n', kind='heteroaromatic')", ade, NULL));
    EXPECT_EQ(1, (int)ring_filter_count("ring('n', size=5)", ade, NULL));

    // Rings fully within the context
    EXPECT_EQ(1, (int)ring_filter_count("ring() in element('C')", trp, NULL));
    EXPECT_EQ(0, (int)ring_filter_count("ring('[nH]') in element('C')", trp, NULL));

    // As one selection: the atoms of all of them
    EXPECT_EQ(1, (int)ring_filter_count("flatten(ring(kind='aromatic'))", trp, &atoms));
    EXPECT_EQ(9, (int)atoms);

    // Errors, found before evaluating
    EXPECT_TRUE(ring_filter_fails_with("ring('c(C)(C)c')", "use smiles() for substituents", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring('CCCCCCCC')", "more than the largest ring", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring('C1CC')", "Invalid SMILES", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(kind='flat')", "Unknown ring kind 'flat'", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(size=6:5)", "ascending", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring('cc', 'aromatic')", "can only be given by name", trp));
}
