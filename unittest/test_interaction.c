#include "utest.h"

#include <md_interaction.h>
#include <md_vial.h>
#include <md_bond.h>
#include <md_pdb.h>
#include <md_system.h>
#include <md_infer.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_geom.h>
#include <core/md_unitcell.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>

#include "vial_test_systems.h"

#include <math.h>

// Acetate 1 and methylammonium 2 in a salt bridge, N...O 2.78 Å; acetate 3 and methylammonium 4 alone. Chloromethane 5
// towards the O of formaldehyde 6, Cl...O 3.00 Å along the C-Cl axis; 7 and 8 at 90 degrees; 9 and 10 at 3.62 Å, beyond
// the 3.27 Å of the van der Waals radii. No hydrogens: the charges are those at pH 7.
static const char interaction__pdb[] =
    "HETATM    1  C1  ACT A   1      -1.520   0.000   0.000  1.00  0.00           C\n"
    "HETATM    2  C2  ACT A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3  O1  ACT A   1       0.625   1.083   0.000  1.00  0.00           O\n"
    "HETATM    4  O2  ACT A   1       0.625  -1.083   0.000  1.00  0.00           O\n"
    "HETATM    5  N1  MAM A   2       3.400   1.300   0.000  1.00  0.00           N\n"
    "HETATM    6  C1  MAM A   2       4.870   1.300   0.000  1.00  0.00           C\n"
    "HETATM    7  C1  ACT A   3      28.480   0.000   0.000  1.00  0.00           C\n"
    "HETATM    8  C2  ACT A   3      30.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    9  O1  ACT A   3      30.625   1.083   0.000  1.00  0.00           O\n"
    "HETATM   10  O2  ACT A   3      30.625  -1.083   0.000  1.00  0.00           O\n"
    "HETATM   11  N1  MAM A   4      60.000   0.000   0.000  1.00  0.00           N\n"
    "HETATM   12  C1  MAM A   4      61.470   0.000   0.000  1.00  0.00           C\n"
    "HETATM   13  C1  CLM A   5       0.000   0.000  10.000  1.00  0.00           C\n"
    "HETATM   14 CL1  CLM A   5       1.780   0.000  10.000  1.00  0.00          CL\n"
    "HETATM   15  O1  FOR A   6       4.780   0.000  10.000  1.00  0.00           O\n"
    "HETATM   16  C1  FOR A   6       6.000   0.000  10.000  1.00  0.00           C\n"
    "HETATM   17  C1  CLM A   7       0.000   0.000  20.000  1.00  0.00           C\n"
    "HETATM   18 CL1  CLM A   7       1.780   0.000  20.000  1.00  0.00          CL\n"
    "HETATM   19  O1  FOR A   8       1.780   3.000  20.000  1.00  0.00           O\n"
    "HETATM   20  C1  FOR A   8       1.780   4.220  20.000  1.00  0.00           C\n"
    "HETATM   21  C1  CLM A   9       0.000   0.000  30.000  1.00  0.00           C\n"
    "HETATM   22 CL1  CLM A   9       1.780   0.000  30.000  1.00  0.00          CL\n"
    "HETATM   23  O1  FOR A  10       5.400   0.000  30.000  1.00  0.00           O\n"
    "HETATM   24  C1  FOR A  10       6.620   0.000  30.000  1.00  0.00           C\n"
    "END\n";

enum { ACT1_O1 = 2, MAM2_N1 = 4, CLM5_CL = 13, FOR6_O = 14 };

static bool interaction__load(md_system_t* sys, md_system_state_t* st, md_allocator_i* alloc) {
    *sys = (md_system_t){.alloc = alloc};
    *st = (md_system_state_t){.alloc = alloc};
    if ((md_pdb_system_load_from_str(sys, st, (md_str_t){interaction__pdb, sizeof(interaction__pdb) - 1}, NULL) != MD_LOAD_OK)) return false;
    md_system_infer(sys, st, NULL);
    return true;
}

UTEST(interaction, ionic_and_halogen) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(interaction__load(&sys, &st, alloc));
    ASSERT_EQ(24u, (uint32_t)sys.atom.count);
    EXPECT_EQ(-1, md_atom_formal_charge(&sys.atom, ACT1_O1) + md_atom_formal_charge(&sys.atom, ACT1_O1 + 1));
    EXPECT_EQ(1, md_atom_formal_charge(&sys.atom, MAM2_N1));

    md_interaction_query_t q;
    ASSERT_TRUE(md_interaction_query_init(&q, NULL, &sys, alloc));
    EXPECT_EQ(2u, q.num_cations);
    EXPECT_EQ(2u, q.num_anions);
    EXPECT_EQ(3u, (uint32_t)q.num_halogens);
    EXPECT_TRUE(md_interaction_query_can_find(&q, MD_INTERACTION_IONIC));
    EXPECT_TRUE(md_interaction_query_can_find(&q, MD_INTERACTION_HALOGEN));
    md_interaction_set_t s;
    ASSERT_TRUE(md_interaction_query_eval(&s, &q, &st, alloc));
    ASSERT_EQ(2u, (uint32_t)s.count);
    // Ionic first, the cation's atom first, of the closest two: O1, not O2
    EXPECT_EQ((uint8_t)MD_INTERACTION_IONIC, s.kind[0]);
    EXPECT_EQ((uint32_t)MAM2_N1, s.a[0]);
    EXPECT_EQ((uint32_t)ACT1_O1, s.b[0]);
    EXPECT_NEAR(2.783f, s.dist[0], 1e-3f);
    EXPECT_EQ((uint8_t)MD_INTERACTION_HALOGEN, s.kind[1]);
    EXPECT_EQ((uint32_t)CLM5_CL, s.a[1]);
    EXPECT_EQ((uint32_t)FOR6_O, s.b[1]);
    EXPECT_NEAR(3.0f, s.dist[1], 1e-3f);
    md_interaction_query_free(&q);

    // One kind; a gate let wider: 9 and 10 at 3.62 Å, past the radii by 0.35 Å
    md_interaction_params_t p = md_interaction_params_default();
    p.kinds = MD_INTERACTION_HALOGEN;
    p.halogen_tolerance = 0.4f;
    ASSERT_TRUE(md_interaction_query_init(&q, &(md_interaction_desc_t){.params = &p}, &sys, alloc));
    EXPECT_FALSE(md_interaction_query_can_find(&q, MD_INTERACTION_IONIC));
    ASSERT_TRUE(md_interaction_query_eval(&s, &q, &st, alloc));
    EXPECT_EQ(2u, (uint32_t)s.count);
    md_interaction_query_free(&q);

    // The angle gate let down to 80 degrees: 7 and 8 at 90 too
    p = md_interaction_params_default();
    p.halogen_min_cxa = 80.0f;
    ASSERT_TRUE(md_interaction_query_init(&q, &(md_interaction_desc_t){.params = &p}, &sys, alloc));
    ASSERT_TRUE(md_interaction_query_eval(&s, &q, &st, alloc));
    EXPECT_EQ(3u, (uint32_t)s.count);
    md_interaction_query_free(&q);

    // Within a set: the halogen's pair alone
    const uint32_t sel[] = {12, 13, 14, 15};
    md_bitset_t set = md_bitset_from_indices(sel, 4, alloc);
    ASSERT_TRUE(md_interaction_query_init(&q, &(md_interaction_desc_t){.set = &set}, &sys, alloc));
    ASSERT_TRUE(md_interaction_query_eval(&s, &q, &st, alloc));
    ASSERT_EQ(1u, (uint32_t)s.count);
    EXPECT_EQ((uint8_t)MD_INTERACTION_HALOGEN, s.kind[0]);
    md_interaction_query_free(&q);

    // Across the cell: the methylammonium moved by a period, the bridge is the same
    md_system_state_t moved = st;
    moved.xyz = md_alloc(alloc, st.num_atoms * sizeof(md_float3_t));
    for (size_t i = 0; i < st.num_atoms; ++i) moved.xyz[i] = st.xyz[i];
    moved.unitcell = md_unitcell_from_extent(100.0, 100.0, 100.0);
    moved.xyz[MAM2_N1].x += 100.0f;
    moved.xyz[MAM2_N1 + 1].x += 100.0f;
    ASSERT_TRUE(md_interaction_query_init(&q, NULL, &sys, alloc));
    ASSERT_TRUE(md_interaction_query_eval(&s, &q, &moved, alloc));
    ASSERT_EQ(2u, (uint32_t)s.count);
    EXPECT_NEAR(2.783f, s.dist[0], 1e-3f);
    md_interaction_query_free(&q);

    md_arena_allocator_destroy(alloc);
}

// Through VIAL, as viamd's representations filter: bond(kind=) finds them per frame, with their kinds in the flags of
// the pairs, of origin inferred
UTEST(interaction, vial_bond) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(interaction__load(&sys, &st, alloc));
    static const struct { const char* expr; uint32_t n; int kind; } cases[] = {   // kind -1: any
        {"bond(kind='ionic')", 1, MD_BOND_KIND_IONIC},
        {"bond(kind='halogen')", 1, MD_BOND_KIND_HALOGEN},
        {"bond(kind='interaction' and not 'hydrogen')", 2, -1},
        {"bond(particle(resid=5..6), kind='interaction')", 1, MD_BOND_KIND_HALOGEN},
        {"bond(particle(resid=[1, 3]), kind='ionic')", 0, -1},
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        md_vial_filter_t f;
        const bool ok = md_vial_filter(&f, md_str_from_cstr(cases[i].expr), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc);
        if (!ok) printf("%s: %.*s %.*s\n", cases[i].expr, 5, f.error.code, (int)f.error.message.len, f.error.message.ptr);
        ASSERT_TRUE(ok);
        EXPECT_TRUE(f.ok);
        EXPECT_TRUE(f.dynamic);
        EXPECT_EQ(2 * cases[i].n, (uint32_t)f.num_sets);
        for (uint32_t k = 0; k < cases[i].n && f.pair_flags; ++k) {
            const md_bond_flags_t fl = (md_bond_flags_t)f.pair_flags[k];
            EXPECT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin(fl));
            if (cases[i].kind >= 0) EXPECT_EQ((md_bond_kind_t)cases[i].kind, md_bond_kind(fl));
        }
    }
    md_arena_allocator_destroy(alloc);
}

// On T, the salt bridges as a brute force over every pair of charged atoms finds them: per pair of groups of opposite
// charges with two atoms within 4 Å, more than 3 bonds apart, one bond, of the closest two
UTEST(interaction, ionic_brute_force) {
    const md_system_t* sys = vial_test_system('T');
    const md_system_state_t* st = vial_test_state('T');
    ASSERT_TRUE(sys && st);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    md_interaction_params_t p = md_interaction_params_default();
    p.kinds = MD_INTERACTION_IONIC;
    md_interaction_query_t q;
    ASSERT_TRUE(md_interaction_query_init(&q, &(md_interaction_desc_t){.params = &p}, sys, alloc));
    EXPECT_GT(q.num_cations, 0u);
    EXPECT_GT(q.num_anions, 0u);
    md_interaction_set_t s;
    ASSERT_TRUE(md_interaction_query_eval(&s, &q, st, alloc));
    EXPECT_GT(s.count, 0u);

    // By brute force: the closest pair of each pair of groups
    const size_t n = q.num_charged;
    typedef struct { uint32_t ga, gb, a, b; float d; } best_t;
    md_array(best_t) best = 0;
    for (size_t i = 0; i < n; ++i) {
        if (q.charged_sign[i] != 1) continue;
        for (size_t j = 0; j < n; ++j) {
            if (q.charged_sign[j] != -1) continue;
            md_float3_t d = md_float3_sub(st->xyz[q.charged_atom[j]], st->xyz[q.charged_atom[i]]);
            md_geom_min_image(d.elem, 0, 1, &st->unitcell);
            const float l = md_float3_length(d);
            if (l > 4.0f) continue;
            bool excluded = false;
            for (uint32_t k = q.excl_off[i]; k < q.excl_off[i + 1]; ++k) excluded |= q.excl_atom[k] == q.charged_atom[j];
            if (excluded) continue;
            size_t k = 0;
            while (k < md_array_size(best) && !(best[k].ga == q.charged_group[i] && best[k].gb == q.charged_group[j])) ++k;
            if (k == md_array_size(best)) {
                const best_t b = {q.charged_group[i], q.charged_group[j], q.charged_atom[i], q.charged_atom[j], l};
                md_array_push(best, b, alloc);
            } else if (l < best[k].d) {
                best[k].a = q.charged_atom[i];
                best[k].b = q.charged_atom[j];
                best[k].d = l;
            }
        }
    }
    ASSERT_EQ(md_array_size(best), s.count);
    for (size_t k = 0; k < md_array_size(best); ++k) {
        bool found = false;
        for (size_t m = 0; m < s.count && !found; ++m) found = s.a[m] == best[k].a && s.b[m] == best[k].b && fabsf(s.dist[m] - best[k].d) < 1e-4f;
        EXPECT_TRUE(found);
    }
    printf("T: %u cationic and %u anionic groups, %zu salt bridges\n", q.num_cations, q.num_anions, s.count);
    md_interaction_query_free(&q);
    md_arena_allocator_destroy(alloc);
}
