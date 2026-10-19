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
#include <core/md_bitset.h>
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
    const bool ok = md_match_for_each(&count, NULL, &desc, sys, NULL, NULL);
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
        if (!md_match_for_each(&count, NULL, &desc, sys, NULL, NULL)) count = SIZE_MAX;
        md_match_query_free(&query);
    }
    md_arena_allocator_destroy(alloc);
    return count;
}

// The number of selections ring() and the like give, SIZE_MAX on an error
static size_t ring_filter_count(const char* expr, const ring_sys_t* r, size_t* out_atoms) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_array(md_bitset_t) arr = 0;
    size_t result = SIZE_MAX;
    const md_filter_context_t ctx = {&r->sys, &r->st, NULL};
    const md_filter_result_t res = md_filter_sets(&arr, str_from_cstr(expr), &ctx, alloc);
    if (res.success) {
        result = md_array_size(arr);
        if (out_atoms) {
            md_bitmask_t all = md_bitmask_create(0, (uint32_t)r->sys.atom.count, alloc);
            for (size_t i = 0; i < result; ++i) md_bitmask_or(&all, arr[i]);
            *out_atoms = md_bitmask_count(&all);
        }
    } else {
        printf("'%s': %s\n", expr, res.error);
    }
    md_arena_allocator_destroy(alloc);
    return result;
}

static bool ring_filter_fails_with(const char* expr, const char* expected, const ring_sys_t* r) {
    md_bitset_t bf = {0};
    const md_filter_context_t ctx = {&r->sys, &r->st, NULL};
    const md_filter_result_t res = md_filter(&bf, str_from_cstr(expr), &ctx, md_get_heap_allocator());
    md_bitset_free(&bf, md_get_heap_allocator());
    if (res.success) {
        printf("'%s' was expected to fail\n", expr);
        return false;
    }
    if (!strstr(res.error, expected)) {
        printf("'%s' failed with '%s', expected '%s'\n", expr, res.error, expected);
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
    ASSERT_TRUE(md_match_for_each(&count, NULL, &desc, trp, NULL, NULL));
    EXPECT_EQ(2, (int)count);
    // A unit beyond the atoms of the system is an error
    units[2] = md_bitset_range(0, (uint32_t)trp->atom.count + 1);
    EXPECT_FALSE(md_match_for_each(&count, NULL, &desc, trp, NULL, NULL));

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
    EXPECT_EQ(1, (int)ring_filter_count("union(ring(kind='aromatic'))", trp, &atoms));
    EXPECT_EQ(9, (int)atoms);

    // Errors, found before evaluating
    EXPECT_TRUE(ring_filter_fails_with("ring('c(C)(C)c')", "use smiles() for substituents", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring('CCCCCCCC')", "more than the largest ring", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring('C1CC')", "Invalid SMILES", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(kind='flat')", "Unknown ring kind 'flat'", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(size=6:5)", "ascending", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring('cc', 'aromatic')", "can only be given by name", trp));
}

// The rings ring() and the like give, as indices into md_system_t::ring (-1 for a selection which is no ring), in
// their order. Returns how many, -1 on an error.
static int ring_filter_indices(int* out, int cap, const char* expr, const ring_sys_t* r) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_array(md_bitset_t) arr = 0;
    int result = -1;
    const md_filter_context_t ctx = {&r->sys, &r->st, NULL};
    const md_filter_result_t res = md_filter_sets(&arr, str_from_cstr(expr), &ctx, alloc);
    if (res.success) {
        result = (int)md_array_size(arr);
        for (int i = 0; i < result && i < cap; ++i) {
            out[i] = -1;
            for (size_t k = 0; k < md_index_data_num_ranges(&r->sys.ring); ++k) {
                const int32_t* ring = md_index_range_beg(&r->sys.ring, k);
                const size_t n = md_index_range_size(&r->sys.ring, k);
                bool same = md_bitset_count(arr[i]) == n;
                for (size_t j = 0; j < n && same; ++j) same = md_bitset_test(arr[i], (uint32_t)ring[j]);
                if (same) { out[i] = (int)k; break; }
            }
        }
    } else {
        printf("'%s': %s\n", expr, res.error);
    }
    md_arena_allocator_destroy(alloc);
    return result;
}

// ring(index): by position among the rings the filters let through, within the context. Without either, ring(k) is
// the ring k - 1 of the system, which is what the circle of an aromatic ring picks as (md_gl) and the tooltip names.
UTEST_F(ring, script_ring_index) {
    const ring_sys_t* trp = &utest_fixture->trp;
    const ring_sys_t* ade = &utest_fixture->ade;
    ASSERT_EQ(2, (int)md_index_data_num_ranges(&trp->sys.ring));
    ASSERT_EQ(3, (int)md_index_data_num_ranges(&ade->sys.ring));
    int idx[8];

    // The rings of the system, in its order
    EXPECT_EQ(1, ring_filter_indices(idx, 8, "ring(1)", trp));
    EXPECT_EQ(0, idx[0]);
    EXPECT_EQ(1, ring_filter_indices(idx, 8, "ring(2)", trp));
    EXPECT_EQ(1, idx[0]);
    EXPECT_EQ(2, ring_filter_indices(idx, 8, "ring(1:2)", trp));
    EXPECT_EQ(0, idx[0]);
    EXPECT_EQ(1, idx[1]);
    EXPECT_EQ(3, ring_filter_indices(idx, 8, "ring(:)", ade));
    EXPECT_EQ(2, ring_filter_indices(idx, 8, "ring(2:)", ade));
    EXPECT_EQ(1, idx[0]);
    EXPECT_EQ(2, idx[1]);
    // In the order given, as residue() takes them
    EXPECT_EQ(2, ring_filter_indices(idx, 8, "ring({3, 1})", ade));
    EXPECT_EQ(2, idx[0]);
    EXPECT_EQ(0, idx[1]);

    // Among the rings the filters let through: the 6 membered ring of indole is ring(1, size=6) whatever its index
    const int six = md_index_range_size(&trp->sys.ring, 0) == 6 ? 0 : 1;
    EXPECT_EQ(1, ring_filter_indices(idx, 8, "ring(1, size=6)", trp));
    EXPECT_EQ(six, idx[0]);
    EXPECT_EQ(1, ring_filter_indices(idx, 8, "ring(1, kind='heteroaromatic')", trp));
    EXPECT_EQ(1 - six, idx[0]);
    EXPECT_EQ(2, ring_filter_indices(idx, 8, "ring(:, kind='heteroaromatic')", ade));

    // Within a context: the rings within it. The benzene ring is the only one of indole made of C alone.
    EXPECT_EQ(1, ring_filter_indices(idx, 8, "ring(1) in element('C')", trp));
    EXPECT_EQ(six, idx[0]);
    EXPECT_EQ(0, ring_filter_indices(idx, 8, "ring(2) in element('C')", trp));    // Past the last: nothing
    EXPECT_EQ(1, ring_filter_indices(idx, 8, "ring(1) in residue(:)", trp));

    // The same selection as the patterns give
    size_t atoms = 0;
    EXPECT_EQ(1, (int)ring_filter_count("union(ring(1:2))", trp, &atoms));
    EXPECT_EQ(9, (int)atoms);

    // Errors, found before evaluating
    EXPECT_TRUE(ring_filter_fails_with("ring(3)", "past the last ring (2)", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(1:5)", "past the last ring (2)", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(2, kind='heteroaromatic')", "past the last ring (1)", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(1)", "there are none", &utest_fixture->ala));
    EXPECT_TRUE(ring_filter_fails_with("ring(0)", "before the first ring", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(0) in element('C')", "before the first ring", trp));
    EXPECT_TRUE(ring_filter_fails_with("ring(1, 'aromatic')", "can only be given by name", trp));
}
