#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_bitset.h>
#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_mmcif.h>
#include <md_util.h>

#include <md_filter.h>
#include <md_script.h>

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

    // An operator without its left operand: an error at it, not a crash
    res = md_filter(&bf, STR_LIT("resname('ALA') or *2"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_EQ(18, res.error_range.beg);
    EXPECT_EQ(19, res.error_range.end);

    // An argument list which failed to parse is not taken for an empty one
    res = md_filter(&bf, STR_LIT("all(<2)"), &fctx, alloc);
    EXPECT_FALSE(res.success);

    // Something that picks neither atoms nor bonds
    res = md_filter(&bf, STR_LIT("distance(1, 2)"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "picks neither atoms nor bonds") != NULL);

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

// Within a context the expression is the left hand side of 'in' with the context as one structure: the same as
// (expression) in union(context) wherever 'in' can be written, and also where it cannot, for a context which
// depends on the coordinates
UTEST(filter, in_context) {
    const str_t gro_file = STR_INIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    const md_filter_context_t whole = {&sys, &state, NULL, NULL};

    static const char* cases[][2] = {
        // context              expression
        {"chain(2)",            "residue(1)"},      // relative to the context: the first residue of the chain
        {"chain(2)",            "residue(1:3)"},
        {"chain(2)",            "all"},
        {"chain(2)",            "element('C')"},
        {"protein",             "residue(1)"},      // one context, however many structures it is made of
        {"resname('PFT')",      "residue(1)"},
        {"chain(1:2)",          "chain(2)"},
        {"chain(1)",            "out resname('PFT')"},
    };

    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_bitset_t context = {0};
        ASSERT_TRUE(md_filter(&context, str_from_cstr(cases[i][0]), &whole, alloc).success);

        const md_filter_context_t within = {&sys, &state, NULL, &context};
        md_bitset_t in_context = {0};
        EXPECT_TRUE(md_filter(&in_context, str_from_cstr(cases[i][1]), &within, alloc).success);

        char buf[256];
        snprintf(buf, sizeof(buf), "(%s) in union(%s)", cases[i][1], cases[i][0]);
        md_bitset_t in_union = {0};
        EXPECT_TRUE(md_filter(&in_union, str_from_cstr(buf), &whole, alloc).success);

        EXPECT_TRUE(md_bitset_count(in_context) > 0);
        EXPECT_TRUE(md_bitset_equal(in_union, in_context));
    }

    md_bitset_t context = {0};
    md_bitset_t bf = {0};
    md_filter_result_t res;

    // The context is what the expression sees of the system, its validation included
    ASSERT_TRUE(md_filter(&context, STR_LIT("chain(1)"), &whole, alloc).success);
    const md_filter_context_t chain = {&sys, &state, NULL, &context};
    EXPECT_TRUE(md_filter(&bf, STR_LIT("all"), &chain, alloc).success);
    EXPECT_TRUE(md_bitset_equal(context, bf));
    EXPECT_TRUE(md_filter(&bf, STR_LIT("  "), &chain, alloc).success);     // Nothing written is 'all', so the context
    EXPECT_TRUE(md_bitset_equal(context, bf));
    res = md_filter(&bf, STR_LIT("resname('PFT')"), &chain, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "within the context") != NULL);
    EXPECT_LT(res.error_range.beg, res.error_range.end);     // and where, in the expression as it was written

    // A context which depends on the coordinates: 'in' cannot have it, as it evaluates its context when compiling
    md_bitset_t near = {0};
    res = md_filter(&near, STR_LIT("within(5, resname('PFT'))"), &whole, alloc);
    ASSERT_TRUE(res.success);
    EXPECT_TRUE(res.is_dynamic);
    EXPECT_FALSE(md_filter(&bf, STR_LIT("(all) in (within(5, resname('PFT')))"), &whole, alloc).success);

    const md_filter_context_t around = {&sys, &state, NULL, &near};
    res = md_filter(&bf, STR_LIT("all"), &around, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_TRUE(md_bitset_equal(near, bf));
    EXPECT_FALSE(res.is_dynamic);   // Of the expression alone: the context is the caller's to know about

    md_bitset_t carbons = {0};
    EXPECT_TRUE(md_filter(&carbons, STR_LIT("element('C') and within(5, resname('PFT'))"), &whole, alloc).success);
    EXPECT_TRUE(md_filter(&bf, STR_LIT("element('C')"), &around, alloc).success);
    EXPECT_TRUE(md_bitset_count(bf) > 0);
    EXPECT_TRUE(md_bitset_equal(carbons, bf));

    md_arena_allocator_destroy(alloc);
}

static bool completes_with(md_script_completions_t c, const char* text) {
    for (size_t i = 0; i < c.count; ++i) {
        if (c.items[i].kind == MD_SCRIPT_COMPLETION_VALUE && str_eq(c.items[i].text, str_from_cstr(text))) return true;
    }
    return false;
}

// The values a selector is completed with are those within the context: one that is not there would not match
UTEST(filter, completion_within_context) {
    const str_t gro_file = STR_INIT(MD_UNITTEST_DATA_DIR "/centered.gro");
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, gro_file));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    const md_filter_context_t whole = {&sys, &state, NULL, NULL};
    md_bitset_t context = {0};
    ASSERT_TRUE(md_filter(&context, STR_LIT("chain(1)"), &whole, alloc).success);
    const md_filter_context_t chain = {&sys, &state, NULL, &context};

    const str_t expr = STR_LIT("resname()");
    const int cursor = 8;

    md_script_completions_t c = md_filter_complete(expr, cursor, &whole, alloc);
    EXPECT_TRUE(completes_with(c, "\"PFT\""));
    EXPECT_TRUE(completes_with(c, "\"ALA\""));

    c = md_filter_complete(expr, cursor, &chain, alloc);
    EXPECT_FALSE(completes_with(c, "\"PFT\""));
    EXPECT_TRUE(completes_with(c, "\"ALA\""));

    md_arena_allocator_destroy(alloc);
}

// An acetate with a magnesium on one of its oxygens, a water and a copper pair, and a hydrogen bond a user marks from
// a water hydrogen to the other oxygen
static const char filter_bonds_pdb[] =
    "HETATM    1  C1  ACT A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    2  C2  ACT A   1       1.520   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3  O1  ACT A   1       2.150   1.080   0.000  1.00  0.00           O\n"
    "HETATM    4  O2  ACT A   1       2.150  -1.080   0.000  1.00  0.00           O\n"
    "HETATM    5  MG  MG  A   2       3.198   2.877   0.000  1.00  0.00          MG\n"
    "HETATM    6  O   HOH A   3       0.000   6.000   0.000  1.00  0.00           O\n"
    "HETATM    7  H1  HOH A   3       0.760   6.590   0.000  1.00  0.00           H\n"
    "HETATM    8  H2  HOH A   3      -0.760   6.590   0.000  1.00  0.00           H\n"
    "HETATM    9  CU  CU  A   4      20.000   0.000   0.000  1.00  0.00          CU\n"
    "HETATM   10  CU  CU  A   5      22.450   0.000   0.000  1.00  0.00          CU\n"
    "END\n";

static bool has_pair(const md_atom_pair_t* pairs, int a, int b) {
    for (size_t i = 0; i < md_array_size(pairs); ++i) {
        if ((pairs[i].idx[0] == a && pairs[i].idx[1] == b) || (pairs[i].idx[0] == b && pairs[i].idx[1] == a)) return true;
    }
    return false;
}

static bool offers(md_script_completions_t c, const char* label) {
    for (size_t i = 0; i < c.count; ++i) {
        if (c.items[i].kind == MD_SCRIPT_COMPLETION_PROCEDURE && str_eq(c.items[i].label, str_from_cstr(label))) return true;
    }
    return false;
}

// A filter of bonds picks those of a bond set, or the persistent bonds among the atoms of a set of atoms
UTEST(filter, bonds) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(filter_bonds_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    // C1-C2, C2-O1, C2-O2, Mg-O1, O-H1, O-H2, Cu-Cu
    ASSERT_EQ(7, (int)sys.bond.count);
    md_system_bond_insert(&sys, 6, 3, md_bond_flags_set_origin(md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_KIND_HYDROGEN), MD_BOND_ORIGIN_USER));
    ASSERT_EQ(1, (int)sys.bond.noncovalent_count);

    const md_filter_context_t fctx = {&sys, &state, NULL, NULL};
    md_array(md_atom_pair_t) pairs = 0;
    md_filter_result_t res;

    // A bond set
    res = md_filter_bonds(&pairs, STR_LIT("bond(kind='hydrogen')"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_FALSE(res.is_dynamic);
    EXPECT_EQ(1, (int)md_array_size(pairs));
    EXPECT_TRUE(has_pair(pairs, 6, 3));

    res = md_filter_bonds(&pairs, STR_LIT("bond(kind='coordinate') or bond(kind='hydrogen')"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(2, (int)md_array_size(pairs));
    EXPECT_TRUE(has_pair(pairs, 4, 2));

    // A set of atoms: the bonds among them, of both ranges. Cleared first, not appended to.
    res = md_filter_bonds(&pairs, STR_LIT("resname('HOH')"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(2, (int)md_array_size(pairs));
    res = md_filter_bonds(&pairs, STR_LIT("resname('HOH') or resname('ACT')"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(6, (int)md_array_size(pairs));
    EXPECT_TRUE(has_pair(pairs, 6, 3));
    EXPECT_FALSE(has_pair(pairs, 4, 2));

    // Nothing written is every atom, and so every persistent bond
    res = md_filter_bonds(&pairs, STR_LIT(""), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(8, (int)md_array_size(pairs));

    // Within a context, both atoms of every bond are in it
    md_bitset_t water = {0};
    ASSERT_TRUE(md_filter(&water, STR_LIT("resname('HOH')"), &fctx, alloc).success);
    const md_filter_context_t in_water = {&sys, &state, NULL, &water};
    res = md_filter_bonds(&pairs, STR_LIT("bond()"), &in_water, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(2, (int)md_array_size(pairs));
    // also those of a bond_delta(), whose scope (the context) takes in the bonds with one atom in it: the hydrogen bond
    // from the water to the acetate is not within the water
    res = md_filter_bonds(&pairs, STR_LIT("bond(bond_delta())"), &in_water, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(2, (int)md_array_size(pairs));
    EXPECT_FALSE(has_pair(pairs, 6, 3));

    // What is neither fails, and says so
    res = md_filter_bonds(&pairs, STR_LIT("distance(1, 2)"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_EQ(0, (int)md_array_size(pairs));
    EXPECT_NE('\0', res.error[0]);

    // Both filters are offered a bond set where a selection goes: a filter of atoms takes the atoms of the bonds
    md_script_completions_t c = md_filter_complete(STR_LIT("bon"), 3, &fctx, alloc);
    EXPECT_TRUE(offers(c, "bond"));

    // A spring of the force field (an elastic bond) and the bond of a metal (a metallic bond), here the Cu-Cu one, are
    // kept with the bonds but are no bonds of a molecule. Every atom is still every bond, as bond() is: nothing written
    // leaves none out. A bond_delta() leaves them out: they neither break nor form.
    const md_bond_idx_t cu = md_system_bond_find(&sys, 8, 9);
    ASSERT_GE(cu, 0);
    const struct { md_bond_kind_t kind; const char* expr; } outside[] = {
        { MD_BOND_KIND_ELASTIC, "bond(kind='elastic')" }, { MD_BOND_KIND_METALLIC, "bond(kind='metallic')" },
    };
    for (size_t k = 0; k < ARRAY_SIZE(outside); ++k) {
        sys.bond.flags[cu] = md_bond_flags_set_kind(sys.bond.flags[cu], outside[k].kind);
        md_system_bond_build_connectivity(&sys);
        res = md_filter_bonds(&pairs, STR_LIT(""), &fctx, alloc);
        EXPECT_TRUE(res.success);
        EXPECT_EQ(8, (int)md_array_size(pairs));
        EXPECT_TRUE(has_pair(pairs, 8, 9));
        res = md_filter_bonds(&pairs, STR_LIT("bond()"), &fctx, alloc);
        EXPECT_TRUE(res.success);
        EXPECT_EQ(8, (int)md_array_size(pairs));
        res = md_filter_bonds(&pairs, STR_LIT("resname('CU')"), &fctx, alloc);
        EXPECT_TRUE(res.success);
        EXPECT_EQ(1, (int)md_array_size(pairs));
        EXPECT_TRUE(has_pair(pairs, 8, 9));
        res = md_filter_bonds(&pairs, str_from_cstr(outside[k].expr), &fctx, alloc);
        EXPECT_TRUE(res.success);
        EXPECT_EQ(1, (int)md_array_size(pairs));
        EXPECT_TRUE(has_pair(pairs, 8, 9));
        res = md_filter_bonds(&pairs, STR_LIT("bond(kind='covalent') or bond(bond_delta())"), &fctx, alloc);
        EXPECT_TRUE(res.success);
        EXPECT_FALSE(has_pair(pairs, 8, 9));
    }

    md_arena_allocator_destroy(alloc);
}

// A filter of atoms which picks bonds picks the atoms of those bonds; within an expression, bonds and atoms do not mix
UTEST(filter, bonds_as_atoms) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(filter_bonds_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    md_system_bond_insert(&sys, 6, 3, md_bond_flags_set_origin(md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_KIND_HYDROGEN), MD_BOND_ORIGIN_USER));

    const md_filter_context_t fctx = {&sys, &state, NULL, NULL};
    md_bitset_t atoms = {0};
    md_filter_result_t res;

    // The hydrogen bond a user marked: H1 of the water and O2 of the acetate, nothing else of their residues
    res = md_filter(&atoms, STR_LIT("bond(origin='user')"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(2, (int)md_bitset_count(atoms));
    EXPECT_TRUE(md_bitset_test(atoms, 6));
    EXPECT_TRUE(md_bitset_test(atoms, 3));

    // The complement of the atoms of the hydrogen bonds is not that of the bonds: every other bond has those atoms too
    res = md_filter(&atoms, STR_LIT("not atom(bond(kind='hydrogen'))"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(8, (int)md_bitset_count(atoms));
    EXPECT_FALSE(md_bitset_test(atoms, 6));
    // here every atom: each has a bond which is not a hydrogen bond
    res = md_filter(&atoms, STR_LIT("not bond(kind='hydrogen')"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(10, (int)md_bitset_count(atoms));

    // One set of atoms per set of bonds: the bonds within each residue, which the Mg and the Cu have none of
    md_array(md_bitset_t) sets = 0;
    res = md_filter_sets(&sets, STR_LIT("bond(residue(:))"), &fctx, alloc);
    EXPECT_TRUE(res.success);
    ASSERT_EQ(sys.component.count, md_array_size(sets));
    EXPECT_EQ(4, (int)md_bitset_count(sets[0]));
    EXPECT_EQ(0, (int)md_bitset_count(sets[1]));
    EXPECT_EQ(3, (int)md_bitset_count(sets[2]));
    EXPECT_EQ(0, (int)md_bitset_count(sets[3]));

    // Within a context: of a set of bonds made without it (a variable of a script), only the atoms in the context. Of
    // the Mg-O1 bond, O1 is in the acetate and the Mg is not.
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("coord_b = bond(kind='coordinate');"), &sys, NULL));
    md_bitset_t acetate = {0};
    ASSERT_TRUE(md_filter(&acetate, STR_LIT("resname('ACT')"), &fctx, alloc).success);
    const md_filter_context_t in_acetate = {&sys, &state, ir, &acetate};
    res = md_filter(&atoms, STR_LIT("coord_b"), &in_acetate, alloc);
    EXPECT_TRUE(res.success);
    EXPECT_EQ(1, (int)md_bitset_count(atoms));
    EXPECT_TRUE(md_bitset_test(atoms, 2));

    // Within an expression they do not mix, and the error says how to write it
    res = md_filter(&atoms, STR_LIT("resname('HOH') and bond()"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "atom(b)") != NULL);
    res = md_filter(&atoms, STR_LIT("within(3.0, bond(kind='hydrogen'))"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "atom(b)") != NULL);
    res = md_filter(&atoms, STR_LIT("within(3.0, atom(bond(kind='hydrogen')))"), &fctx, alloc);
    EXPECT_TRUE(res.success);

    md_arena_allocator_destroy(alloc);
}

// chain("A") is every atom which carries the author's chain id A, as one set: the polymer and the ligands, ions and
// waters filed under it, wherever they are in the file. chain(i) is the i-th chain, one set per chain.
UTEST(filter, chains_by_author_id) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    const md_filter_context_t fctx = {&sys, &state, NULL, NULL};

    // The atoms whose instance the author filed under A
    md_bitmask_t want = md_bitmask_create(0, (uint32_t)sys.atom.count, alloc);
    for (size_t i = 0; i < md_system_instance_count(&sys); ++i) {
        if (str_eq(md_system_instance_auth_id(&sys, i), STR_LIT("A"))) {
            const md_urange_t r = md_system_instance_atom_range(&sys, i);
            md_bitmask_set_range(&want, r.beg, r.end);
        }
    }
    ASSERT_GT(md_bitmask_count(&want), md_system_instance_atom_count(&sys, 0));   // More than the polymer

    md_array(md_bitset_t) sets = 0;
    for (size_t k = 0; k < 3; ++k) {
        static const char* exprs[] = {"chain('A')", "chain_id('A')", "auth_id('A')"};
        md_array_shrink(sets, 0);
        ASSERT_TRUE(md_filter_sets(&sets, str_from_cstr(exprs[k]), &fctx, alloc).success);
        ASSERT_EQ(1, (int)md_array_size(sets));
        EXPECT_TRUE(md_bitset_equal(sets[0], md_bitmask_view(&want)));
    }
    // The same atoms as the first chain by position, which is one set per chain
    md_array_shrink(sets, 0);
    ASSERT_TRUE(md_filter_sets(&sets, STR_LIT("chain(1:2)"), &fctx, alloc).success);
    ASSERT_EQ(2, (int)md_array_size(sets));
    EXPECT_TRUE(md_bitset_equal(sets[0], md_bitmask_view(&want)));
    EXPECT_EQ(sys.atom.count, md_bitset_count(sets[0]) + md_bitset_count(sets[1]));

    // Several ids, or a pattern, are their atoms together
    md_bitset_t bf = {0};
    ASSERT_TRUE(md_filter(&bf, STR_LIT("chain({'A', 'B'})"), &fctx, alloc).success);
    EXPECT_EQ(sys.atom.count, md_bitset_count(bf));
    md_array_shrink(sets, 0);
    ASSERT_TRUE(md_filter_sets(&sets, STR_LIT("chain('*')"), &fctx, alloc).success);
    EXPECT_EQ(1, (int)md_array_size(sets));

    // The chains a selection touches, whole: a water of chain A gives all of chain A
    md_array_shrink(sets, 0);
    ASSERT_TRUE(md_filter_sets(&sets, STR_LIT("chain(water())"), &fctx, alloc).success);
    ASSERT_EQ(2, (int)md_array_size(sets));
    EXPECT_TRUE(md_bitset_equal(sets[0], md_bitmask_view(&want)));

    // Within a context, what of it carries the id. Every water is not of chain A, and 'chain('A') in water()' is an
    // error as 'name("CB") in residue(:)' is where a residue has none: the waters are taken as one here.
    md_bitset_t waters = {0};
    ASSERT_TRUE(md_filter(&waters, STR_LIT("water()"), &fctx, alloc).success);
    const md_filter_context_t in_water = {&sys, &state, NULL, &waters};
    ASSERT_TRUE(md_filter(&bf, STR_LIT("chain('A')"), &in_water, alloc).success);
    md_bitset_t water_a = {0};
    ASSERT_TRUE(md_filter(&water_a, STR_LIT("water() and chain('A')"), &fctx, alloc).success);
    EXPECT_GT(md_bitset_count(water_a), 0);
    EXPECT_TRUE(md_bitset_equal(bf, water_a));

    // An id there is not: the ids there are
    md_filter_result_t res = md_filter(&bf, STR_LIT("chain('C')"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "'A', 'B'") != NULL);

    // Instances are not a word of the script
    EXPECT_FALSE(md_filter(&bf, STR_LIT("instance('E')"), &fctx, alloc).success);

    // Structure matching within a chain: the first amino acid of each of the two
    md_array(md_bitset_t) per_chain = 0;
    ASSERT_TRUE(md_filter_sets(&per_chain, STR_LIT("smiles('NCC=O', level='chain', mode='one_per_unit')"), &fctx, alloc).success);
    EXPECT_EQ(2, (int)md_array_size(per_chain));

    md_arena_allocator_destroy(alloc);
}

// Where the data gives no chain ids none are made up: the polymers are chains found by position
UTEST(filter, chains_without_ids) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    const md_filter_context_t fctx = {&sys, &state, NULL, NULL};
    ASSERT_GT(md_system_chain_count(&sys), 1);

    md_bitset_t bf = {0};
    md_filter_result_t res = md_filter(&bf, STR_LIT("chain('A')"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "no chain ids") != NULL);

    md_array(md_bitset_t) sets = 0;
    ASSERT_TRUE(md_filter_sets(&sets, STR_LIT("chain(:)"), &fctx, alloc).success);
    EXPECT_EQ(md_system_chain_count(&sys), md_array_size(sets));
    ASSERT_TRUE(md_filter(&bf, STR_LIT("chain(2)"), &fctx, alloc).success);
    EXPECT_EQ(md_system_chain_atom_count(&sys, 1), md_bitset_count(bf));

    md_arena_allocator_destroy(alloc);
}

// A molecule is atoms joined by bonds, closed over its residues and its polymer chains: molecule(i) by position, one set
// each, molecule(sel) grown to whole molecules
UTEST(filter, molecules) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    const md_filter_context_t fctx = {&sys, &state, NULL, NULL};
    const size_t num = md_structure_count(&sys.structure);
    ASSERT_GT(num, 1);

    md_array(md_bitset_t) sets = 0;
    ASSERT_TRUE(md_filter_sets(&sets, STR_LIT("molecule(:)"), &fctx, alloc).success);
    ASSERT_EQ(num, md_array_size(sets));
    size_t total = 0;
    for (size_t i = 0; i < num; ++i) total += md_bitset_count(sets[i]);
    EXPECT_EQ(sys.atom.count, total);

    // The first chain is a molecule of its own: growing an atom of it gives all of it
    md_bitset_t chain = {0};
    ASSERT_TRUE(md_filter(&chain, STR_LIT("chain(1)"), &fctx, alloc).success);
    md_array_shrink(sets, 0);
    ASSERT_TRUE(md_filter_sets(&sets, STR_LIT("molecule(atom(1))"), &fctx, alloc).success);
    ASSERT_EQ(1, (int)md_array_size(sets));
    EXPECT_TRUE(md_bitset_equal(sets[0], chain));

    // Counted in molecules, not in instances
    md_script_ir_t* ir = md_script_ir_create(alloc);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("a = count(all(), 'molecule');"), &sys, NULL));
    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("b = count(all(), 'instance');"), &sys, NULL));
    md_script_ir_free(ir);

    md_bitset_t bf = {0};
    md_filter_result_t res = md_filter(&bf, str_from_cstr("molecule(100000)"), &fctx, alloc);
    EXPECT_FALSE(res.success);
    EXPECT_TRUE(strstr(res.error, "past the last molecule") != NULL);

    md_arena_allocator_destroy(alloc);
}
