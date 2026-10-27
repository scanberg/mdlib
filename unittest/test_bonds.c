#include "utest.h"

// Covalent and coordination bond inference: against force field topologies, on structures that clash, at metal sites

#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_mmcif.h>
#include <md_tpr.h>
#include <md_lammps.h>
#include <md_util.h>
#include <md_xyz.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_os.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

static inline uint64_t pair_key(md_atom_pair_t p) {
    const uint32_t a = (uint32_t)MIN(p.idx[0], p.idx[1]);
    const uint32_t b = (uint32_t)MAX(p.idx[0], p.idx[1]);
    return ((uint64_t)a << 32) | b;
}

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// Bonds of the system (from its topology) against what inference finds from the coordinates alone
static void expect_inference_matches_topology(int* utest_result, const md_system_t* sys, const md_system_state_t* st, md_allocator_i* alloc) {
    const size_t n = sys->bond.count;
    uint64_t* ref = md_alloc(alloc, sizeof(uint64_t) * (n + 1));
    for (size_t i = 0; i < n; ++i) ref[i] = pair_key(sys->bond.pairs[i]);
    qsort(ref, n, sizeof(uint64_t), compare_u64);

    md_bond_data_t inferred = { 0 };
    md_util_infer_covalent_bonds(&inferred, st, sys, alloc);
    ASSERT_EQ(n, inferred.count);
    uint64_t* inf = md_alloc(alloc, sizeof(uint64_t) * (inferred.count + 1));
    for (size_t i = 0; i < inferred.count; ++i) inf[i] = pair_key(inferred.pairs[i]);
    qsort(inf, inferred.count, sizeof(uint64_t), compare_u64);
    EXPECT_EQ(0, memcmp(ref, inf, sizeof(uint64_t) * n));
}

UTEST(bonds, topology_tpr) {
    const char* files[] = { MD_UNITTEST_DATA_DIR "/tpr/peptide_tip3p.tpr", MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr" };
    for (size_t f = 0; f < ARRAY_SIZE(files); ++f) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &st, str_from_cstr(files[f])));
        expect_inference_matches_topology(utest_result, &sys, &st, alloc);
        md_vm_arena_destroy(alloc);
    }
}

// The initial configuration of this system has hundreds of atoms of different molecules closer than 1 Å. Without a
// limit on the bonds per atom, inference joins them (900 bonds too many, hydrogens with two or three bonds).
UTEST(bonds, clashing_molecules) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"), NULL));
    EXPECT_TRUE(sys.bond.count > 0 && md_bond_origin(sys.bond.flags[0]) == MD_BOND_ORIGIN_TOPOLOGY);
    expect_inference_matches_topology(utest_result, &sys, &st, alloc);
    md_vm_arena_destroy(alloc);
}

// Acetate with Mg at 2.08 Å from one oxygen and Na at 2.40 Å from the other, a water 2.10 Å from the Mg, and two
// Cu 2.45 Å apart
static const char metal_site_pdb[] =
    "HETATM    1  C1  ACT A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    2  C2  ACT A   1       1.520   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3  O1  ACT A   1       2.150   1.080   0.000  1.00  0.00           O\n"
    "HETATM    4  O2  ACT A   1       2.150  -1.080   0.000  1.00  0.00           O\n"
    "HETATM    5  MG  MG  A   2       3.198   2.877   0.000  1.00  0.00          MG\n"
    "HETATM    6  NA  NA  A   3       3.359  -3.153   0.000  1.00  0.00          NA\n"
    "HETATM    7  O   HOH A   4       3.198   4.977   0.000  1.00  0.00           O\n"
    "HETATM    8  H1  HOH A   4       3.958   5.567   0.000  1.00  0.00           H\n"
    "HETATM    9  H2  HOH A   4       2.438   5.567   0.000  1.00  0.00           H\n"
    "HETATM   10  CU  CU  A   5      20.000   0.000   0.000  1.00  0.00          CU\n"
    "HETATM   11  CU  CU  A   6      22.450   0.000   0.000  1.00  0.00          CU\n"
    "END\n";

static bool has_bond(const md_system_t* sys, int a, int b, md_bond_flags_t* out_flags) {
    const md_bond_idx_t bi = md_system_bond_find(sys, a, b);
    if (bi < 0) return false;
    if (out_flags) *out_flags = sys->bond.flags[bi];
    return true;
}

UTEST(bonds, metal_coordination) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, str_from_cstr(metal_site_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    md_bond_flags_t f = 0;

    // Acetate
    EXPECT_TRUE(has_bond(&sys, 0, 1, NULL));
    EXPECT_TRUE(has_bond(&sys, 1, 2, NULL));
    EXPECT_TRUE(has_bond(&sys, 1, 3, NULL));
    // Mg to the carboxylate: a coordination bond
    ASSERT_TRUE(has_bond(&sys, 2, 4, &f));
    EXPECT_EQ(MD_BOND_KIND_COORDINATE, md_bond_kind(f));
    EXPECT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin(f));
    // Na is a free ion: nothing
    EXPECT_EQ(0, (int)md_bond_conn_count(&sys.bond, 5));
    // Water is never bonded to a metal
    EXPECT_FALSE(has_bond(&sys, 4, 6, NULL));
    // Cu-Cu: a metal-metal bond, which is not coordination
    ASSERT_TRUE(has_bond(&sys, 9, 10, &f));
    EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind(f));

    // A bond between a metal and a non-metal is coordination whatever its origin: the Mg-O bond as a topology gives it
    for (size_t b = 0; b < sys.bond.count; ++b) {
        sys.bond.flags[b] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
    }
    EXPECT_EQ(1, (int)md_util_system_infer_coordination(&sys));
    ASSERT_TRUE(has_bond(&sys, 2, 4, &f));
    EXPECT_EQ(MD_BOND_KIND_COORDINATE, md_bond_kind(f));
    EXPECT_EQ(MD_BOND_ORIGIN_TOPOLOGY, md_bond_origin(f));
    ASSERT_TRUE(has_bond(&sys, 9, 10, &f));
    EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind(f));
    ASSERT_TRUE(has_bond(&sys, 1, 2, &f));
    EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind(f));
    EXPECT_EQ(0, (int)md_util_system_infer_coordination(&sys));
    md_vm_arena_destroy(alloc);
}

// Metal sites of deposited structures: the four cysteines of the zinc fingers (Zn-S 2.4-2.7 Å in this cryo-EM
// model), Mg to the GTP phosphates
UTEST(bonds, metal_sites_deposited) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(2));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    int num_zn = 0;
    for (size_t i = 0; i < sys.atom.count; ++i) {
        if (md_atom_atomic_number(&sys.atom, i) != 30) continue;
        num_zn += 1;
        int num_s = 0;
        md_bond_iter_t it = md_bond_iter(&sys.bond, i);
        while (md_bond_iter_has_next(&it)) {
            if (md_bond_kind(md_bond_iter_bond_flags(&it)) == MD_BOND_KIND_COORDINATE && md_atom_atomic_number(&sys.atom, md_bond_iter_atom_index(&it)) == 16) num_s += 1;
            md_bond_iter_next(&it);
        }
        EXPECT_GE(num_s, 2);
    }
    EXPECT_EQ(2, num_zn);

    md_system_t tub = { .alloc = alloc };
    md_system_state_t tub_st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&tub, &tub_st, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&tub, &tub_st, MD_UTIL_INFER_ALL));
    for (size_t i = 0; i < tub.atom.count; ++i) {
        if (md_atom_atomic_number(&tub.atom, i) == 12) {
            EXPECT_GE((int)md_bond_conn_count(&tub.bond, i), 2);
        }
    }
    md_vm_arena_destroy(alloc);
}

// Inferring the bonds again keeps the bonds a user added by hand
UTEST(bonds, reinference_keeps_user_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const size_t inferred = sys.bond.count;
    const int a = 0, b = (int)sys.atom.count - 1;
    ASSERT_FALSE(has_bond(&sys, a, b, NULL));
    md_system_bond_insert(&sys, a, b, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER));
    md_util_system_infer_covalent_bonds(&sys, &st);
    md_bond_flags_t f = 0;
    EXPECT_TRUE(has_bond(&sys, a, b, &f));
    EXPECT_EQ(MD_BOND_ORIGIN_USER, md_bond_origin(f));
    EXPECT_EQ(inferred + 1, sys.bond.count);
    md_vm_arena_destroy(alloc);
}

// The persistent non-covalent bonds sit in a range of their own after the bond graph (md_bond_data_t): the graph and
// everything derived from it never sees them, and inserting, removing and inferring bonds keeps both ranges
UTEST(bonds, noncovalent_range) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const size_t graph = sys.bond.count;
    const int n = (int)sys.atom.count;
    ASSERT_GE(n, 6);
    EXPECT_EQ(0, (int)sys.bond.noncovalent_count);

    const md_bond_flags_t user = md_bond_flags_set_origin(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_ORIGIN_USER);
    const md_bond_flags_t hbond = md_bond_flags_set_kind(user, MD_BOND_KIND_HYDROGEN);
    const md_bond_flags_t inter = md_bond_flags_set_kind(user, MD_BOND_KIND_INTERACTION);
    ASSERT_FALSE(has_bond(&sys, 0, n - 1, NULL));
    ASSERT_FALSE(has_bond(&sys, 1, n - 2, NULL));
    md_system_bond_insert(&sys, 0, n - 1, hbond);
    md_system_bond_insert(&sys, 1, n - 2, inter);
    EXPECT_EQ(graph, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);
    EXPECT_EQ(graph + 2, md_array_size(sys.bond.pairs));
    EXPECT_EQ(graph + 2, md_array_size(sys.bond.flags));

    // An atom pair has one bond: neither range takes a second one
    md_system_bond_insert(&sys, n - 1, 0, inter);
    md_system_bond_insert(&sys, sys.bond.pairs[0].idx[0], sys.bond.pairs[0].idx[1], hbond);
    EXPECT_EQ(graph, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);

    // A bond of the graph inserted goes before them
    md_system_bond_insert(&sys, 2, n - 3, user);
    EXPECT_EQ(graph + 1, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);
    for (size_t i = 0; i < sys.bond.count; ++i) EXPECT_FALSE(md_bond_is_noncovalent(sys.bond.flags[i]));
    for (size_t i = sys.bond.count; i < sys.bond.count + sys.bond.noncovalent_count; ++i) EXPECT_TRUE(md_bond_is_noncovalent(sys.bond.flags[i]));
    md_system_bond_build_connectivity(&sys);

    // The graph does not have them, the second range does
    EXPECT_FALSE(has_bond(&sys, 0, n - 1, NULL));
    EXPECT_TRUE(has_bond(&sys, 2, n - 3, NULL));
    const md_bond_idx_t hb = md_bond_find_noncovalent(&sys.bond, n - 1, 0);
    ASSERT_GE(hb, (md_bond_idx_t)sys.bond.count);
    EXPECT_EQ(hbond, md_system_bond_flags(&sys, hb));
    EXPECT_EQ(-1, md_bond_find_noncovalent(&sys.bond, 2, n - 3));

    // Inferring the bonds again keeps them, and the user's bond of the graph
    md_util_system_infer_covalent_bonds(&sys, &st);
    EXPECT_EQ(graph + 1, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);
    EXPECT_GE(md_bond_find_noncovalent(&sys.bond, 0, n - 1), 0);
    EXPECT_GE(md_bond_find_noncovalent(&sys.bond, 1, n - 2), 0);
    EXPECT_TRUE(has_bond(&sys, 2, n - 3, NULL));

    // A pair the user marks non-covalent is not given a covalent bond by inference as well
    const md_atom_pair_t covalent = sys.bond.pairs[0];
    ASSERT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin(sys.bond.flags[0]));
    md_system_bond_remove(&sys, 0);
    md_system_bond_build_connectivity(&sys);
    md_system_bond_insert(&sys, covalent.idx[0], covalent.idx[1], hbond);
    EXPECT_EQ(3, (int)sys.bond.noncovalent_count);
    md_util_system_infer_covalent_bonds(&sys, &st);
    EXPECT_FALSE(has_bond(&sys, covalent.idx[0], covalent.idx[1], NULL));
    EXPECT_GE(md_bond_find_noncovalent(&sys.bond, covalent.idx[0], covalent.idx[1]), 0);
    EXPECT_EQ(graph, sys.bond.count);
    EXPECT_EQ(3, (int)sys.bond.noncovalent_count);

    // Removing keeps the ranges: a bond of the graph (the last non-covalent one moves down), then a non-covalent one
    md_system_bond_remove(&sys, 0);
    EXPECT_EQ(graph - 1, sys.bond.count);
    EXPECT_EQ(3, (int)sys.bond.noncovalent_count);
    EXPECT_EQ(graph + 2, md_array_size(sys.bond.pairs));
    for (size_t i = 0; i < sys.bond.count; ++i) EXPECT_FALSE(md_bond_is_noncovalent(sys.bond.flags[i]));
    for (size_t i = sys.bond.count; i < sys.bond.count + sys.bond.noncovalent_count; ++i) EXPECT_TRUE(md_bond_is_noncovalent(sys.bond.flags[i]));
    md_system_bond_remove(&sys, md_bond_find_noncovalent(&sys.bond, 1, n - 2));
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);
    EXPECT_EQ(-1, md_bond_find_noncovalent(&sys.bond, 1, n - 2));
    EXPECT_GE(md_bond_find_noncovalent(&sys.bond, 0, n - 1), 0);
    md_vm_arena_destroy(alloc);
}

// A kind which was given stays: a metal-ligand bond marked covalent by hand is not made coordination
UTEST(bonds, kind_given_stays) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, str_from_cstr(metal_site_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const md_bond_idx_t b = md_system_bond_find(&sys, 2, 4);
    ASSERT_GE(b, 0);
    sys.bond.flags[b] = md_bond_flags_set_origin(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_ORIGIN_USER);
    EXPECT_EQ(0, (int)md_util_system_infer_coordination(&sys));
    EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind(sys.bond.flags[b]));
    // Without it, the kind is perceived
    sys.bond.flags[b] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER);
    EXPECT_EQ(1, (int)md_util_system_infer_coordination(&sys));
    EXPECT_EQ(MD_BOND_KIND_COORDINATE, md_bond_kind(sys.bond.flags[b]));
    md_vm_arena_destroy(alloc);
}

// Every flag combination of a kind is of the kinds the name table says, and of no other
UTEST(bonds, kind_names) {
    size_t count = 0;
    const md_bond_kind_name_t* names = md_bond_kind_names(&count);
    ASSERT_EQ(7, (int)count);
    const md_bond_flags_t given = md_bond_flags_set_origin(MD_BOND_FLAG_KIND_GIVEN | MD_BOND_FLAG_ORDER_PERCEIVED, MD_BOND_ORIGIN_INFERRED);
    const struct { md_bond_flags_t flags; const char* kinds; } cases[] = {
        { MD_BOND_FLAG_NONE,                                                    "covalent" },
        { MD_BOND_FLAG_AROMATIC | 2,                                            "covalent" },
        { (md_bond_flags_t)MD_BOND_KIND_COORDINATE,                             "coordinate" },
        // A spring is of no other kind
        { (md_bond_flags_t)MD_BOND_KIND_ELASTIC,                                "elastic" },
        { (md_bond_flags_t)MD_BOND_KIND_METALLIC,                               "metallic" },
        { (md_bond_flags_t)MD_BOND_KIND_INTERACTION,                            "noncovalent interaction" },
        { (md_bond_flags_t)MD_BOND_KIND_HYDROGEN,                               "noncovalent hydrogen" },
        // The fields around the kind take no part in it
        { md_bond_flags_set_kind(given | MD_BOND_FLAG_DELOCALIZED | 4, MD_BOND_KIND_COVALENT), "covalent" },
        { md_bond_flags_set_kind(given, MD_BOND_KIND_ELASTIC),                  "elastic" },
        { md_bond_flags_set_kind(given, MD_BOND_KIND_HYDROGEN),                 "noncovalent hydrogen" },
    };
    for (size_t c = 0; c < ARRAY_SIZE(cases); ++c) {
        char buf[128] = "";
        for (size_t i = 0; i < count; ++i) {
            if ((cases[c].flags & names[i].mask) == names[i].value) {
                if (buf[0]) strcat(buf, " ");
                strcat(buf, names[i].name);
            }
        }
        EXPECT_STREQ(cases[c].kinds, buf);
    }
}

// The kind is one value of its own: setting it keeps every other field, and its bits say where the bond belongs
UTEST(bonds, kind_code) {
    const md_bond_kind_t kinds[] = { MD_BOND_KIND_COVALENT, MD_BOND_KIND_COORDINATE, MD_BOND_KIND_ELASTIC, MD_BOND_KIND_METALLIC, MD_BOND_KIND_INTERACTION, MD_BOND_KIND_HYDROGEN };
    const bool graph[]       = { true,  true,  false, false, false, false };
    const bool noncovalent[] = { false, false, false, false, true,  true  };
    const md_bond_flags_t rest = md_bond_flags_set_origin(md_bond_flags_set_order(MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_ORDER_PERCEIVED | MD_BOND_FLAG_KIND_GIVEN, 2), MD_BOND_ORIGIN_USER);
    for (size_t i = 0; i < ARRAY_SIZE(kinds); ++i) {
        // Over a bond of another kind, so every bit of the old one goes
        const md_bond_flags_t f = md_bond_flags_set_kind(md_bond_flags_set_kind(rest, kinds[(i + 1) % ARRAY_SIZE(kinds)]), kinds[i]);
        EXPECT_EQ(kinds[i], md_bond_kind(f));
        EXPECT_EQ(rest, f & ~MD_BOND_FLAG_KIND_MASK);
        EXPECT_EQ(graph[i], md_bond_is_graph(f));
        EXPECT_EQ(noncovalent[i], md_bond_is_noncovalent(f));
    }
}

// Two waters, the hydrogen H1a of the first pointing at the oxygen of the second (2.8 Å apart)
static const char water_pair_pdb[] =
    "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
    "HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H\n"
    "HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H\n"
    "HETATM    4  O   HOH A   2       2.800   0.000   0.000  1.00  0.00           O\n"
    "HETATM    5  H1  HOH A   2       3.040   0.930   0.000  1.00  0.00           H\n"
    "HETATM    6  H2  HOH A   2       3.040  -0.930   0.000  1.00  0.00           H\n"
    "END\n";

typedef struct delta_fixture_t {
    md_allocator_i* alloc;
    md_system_t sys;
    md_system_state_t ref;
    md_system_state_t st;   // A copy of ref to move atoms in
} delta_fixture_t;

static bool delta_fixture_init(delta_fixture_t* f) {
    f->alloc = md_vm_arena_create(GIGABYTES(1));
    f->sys = (md_system_t){ .alloc = f->alloc };
    f->ref = (md_system_state_t){ .alloc = f->alloc };
    if (!md_pdb_system_init_from_str(&f->sys, &f->ref, str_from_cstr(water_pair_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE)) return false;
    if (!md_util_system_infer(&f->sys, &f->ref, MD_UTIL_INFER_ALL)) return false;
    f->st = f->ref;
    f->st.xyz = md_alloc(f->alloc, sizeof(vec3_t) * f->ref.num_atoms);
    MEMCPY(f->st.xyz, f->ref.xyz, sizeof(vec3_t) * f->ref.num_atoms);
    return true;
}

static bool delta_has_formed(const md_bond_delta_t* d, int a, int b) {
    for (size_t k = 0; k < d->num_formed; ++k) {
        if (d->formed_pair[k].idx[0] == MIN(a, b) && d->formed_pair[k].idx[1] == MAX(a, b)) return true;
    }
    return false;
}

static const md_bond_delta_desc_t delta_default = { .tolerance = MD_BOND_DELTA_TOLERANCE, .breakable = MD_BOND_DELTA_BREAKABLE };

// The bonds of the reference state are the graph's: nothing formed or broken
UTEST(bonds, delta_reference_is_empty) {
    delta_fixture_t f = {0};
    ASSERT_TRUE(delta_fixture_init(&f));
    ASSERT_EQ(4, (int)f.sys.bond.count);
    md_bond_delta_t d = {0};
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.ref, &delta_default, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    EXPECT_EQ(0, (int)d.num_formed);

    md_system_t tub = { .alloc = f.alloc };
    md_system_state_t tub_st = { .alloc = f.alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&tub, &tub_st, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&tub, &tub_st, MD_UTIL_INFER_ALL));
    md_bond_delta_t td = {0};
    ASSERT_TRUE(md_util_bond_delta(&td, &tub, &tub_st, &delta_default, f.alloc));
    EXPECT_EQ(0, (int)md_bitset_count(td.broken));
    EXPECT_EQ(0, (int)td.num_formed);
    md_vm_arena_destroy(f.alloc);
}

// A proton moves to the other water: its bond to the first breaks, one to the second forms
UTEST(bonds, delta_proton_transfer) {
    delta_fixture_t f = {0};
    ASSERT_TRUE(delta_fixture_init(&f));
    const md_bond_idx_t oh = md_system_bond_find(&f.sys, 0, 1);
    ASSERT_GE(oh, 0);
    f.st.xyz[1] = vec3_set(1.80f, 0.0f, 0.0f);

    md_bond_delta_t d = {0};
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_EQ(1, (int)md_bitset_count(d.broken));
    EXPECT_TRUE(md_bitset_test(d.broken, (uint32_t)oh));
    ASSERT_EQ(1, (int)d.num_formed);
    EXPECT_TRUE(delta_has_formed(&d, 1, 3));
    EXPECT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin(d.formed_flags[0]));
    EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind(d.formed_flags[0]));

    // Only within the scope: the second water alone still sees its bond form (an atom of it is in the scope), the
    // first oxygen alone sees its bond break, an atom far from both sees nothing
    md_bitset_t water2 = md_bitset_range(3, 6);
    md_bond_delta_desc_t desc = delta_default;
    desc.scope = &water2;
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &desc, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    EXPECT_TRUE(delta_has_formed(&d, 1, 3));
    md_bitset_t o1 = md_bitset_range(0, 1);
    desc.scope = &o1;
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &desc, f.alloc));
    EXPECT_TRUE(md_bitset_test(d.broken, (uint32_t)oh));
    EXPECT_EQ(0, (int)d.num_formed);
    md_bitset_t h2 = md_bitset_range(2, 3);
    desc.scope = &h2;
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &desc, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    EXPECT_EQ(0, (int)d.num_formed);

    // A bond of an origin which is not breakable stays
    f.sys.bond.flags[oh] = md_bond_flags_set_origin(f.sys.bond.flags[oh], MD_BOND_ORIGIN_TOPOLOGY);
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    EXPECT_TRUE(delta_has_formed(&d, 1, 3));
    md_vm_arena_destroy(f.alloc);
}

// Within the band around the cutoffs nothing changes: a bond is neither broken nor formed until it is clearly past them
UTEST(bonds, delta_hysteresis) {
    delta_fixture_t f = {0};
    ASSERT_TRUE(delta_fixture_init(&f));
    // O1-H at 1.25 Å: past the cutoff of an O-H of one molecule (1.05 * 1.15 * 0.97 = 1.17) but within the band (1.29)
    f.st.xyz[1] = vec3_set(1.25f, 0.0f, 0.0f);
    md_bond_delta_t d = {0};
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    md_bond_delta_desc_t desc = delta_default;
    desc.tolerance = 0.0f;
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &desc, f.alloc));
    EXPECT_EQ(1, (int)md_bitset_count(d.broken));

    // H to O2 at 1.05 Å: within the cutoff between molecules (1.15 * 0.97 = 1.12) but not past the band (1.00)
    f.st.xyz[1] = vec3_set(1.75f, 0.0f, 0.0f);
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_FALSE(delta_has_formed(&d, 1, 3));
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &desc, f.alloc));
    EXPECT_TRUE(delta_has_formed(&d, 1, 3));
    md_vm_arena_destroy(f.alloc);
}

// The transition: the old bond stretched past its cutoff but not broken, the new one formed. Both show.
UTEST(bonds, delta_transition_shows_both) {
    delta_fixture_t f = {0};
    ASSERT_TRUE(delta_fixture_init(&f));
    // The second water 0.6 Å closer, the proton 1.22 Å from O1 and 0.98 Å from O2
    for (int i = 3; i < 6; ++i) f.st.xyz[i].x -= 0.6f;
    f.st.xyz[1] = vec3_set(1.22f, 0.0f, 0.0f);
    md_bond_delta_t d = {0};
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    EXPECT_TRUE(delta_has_formed(&d, 1, 3));

    // With the old bond within its cutoff, the hydrogen has all it can take: nothing forms
    f.st.xyz[1] = vec3_set(1.10f, 0.0f, 0.0f);
    md_bond_delta_desc_t desc = delta_default;
    desc.tolerance = 0.0f;
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &desc, f.alloc));
    EXPECT_TRUE(md_bitset_empty(d.broken));
    EXPECT_FALSE(delta_has_formed(&d, 1, 3));
    md_vm_arena_destroy(f.alloc);
}

// A pair with a non-covalent bond in the system forms no bond
UTEST(bonds, delta_keeps_noncovalent) {
    delta_fixture_t f = {0};
    ASSERT_TRUE(delta_fixture_init(&f));
    md_system_bond_insert(&f.sys, 1, 3, md_bond_flags_set_origin(md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_KIND_HYDROGEN), MD_BOND_ORIGIN_USER));
    f.st.xyz[1] = vec3_set(1.80f, 0.0f, 0.0f);
    md_bond_delta_t d = {0};
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_EQ(1, (int)md_bitset_count(d.broken));
    EXPECT_EQ(0, (int)d.num_formed);
    md_vm_arena_destroy(f.alloc);
}

// ### METALLIC PHASES ###

// A crystal of n x n x n cubic cells of side a as XYZ text: the atoms of a cell (in fractions of it) with their elements,
// periodic (a Lattice of the crystal) or not
typedef struct crystal_site_t {
    const char* element;
    float x, y, z;
} crystal_site_t;

static str_t crystal_xyz(const crystal_site_t* sites, int num_sites, int n, float a, bool periodic, md_allocator_i* alloc) {
    const int num = num_sites * n * n * n;
    char* buf = md_alloc(alloc, 160 + 64 * (size_t)num);
    int len = 0;
    if (periodic) {
        len += sprintf(buf + len, "%d\nLattice=\"%.4f 0 0 0 %.4f 0 0 0 %.4f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n", num, n * a, n * a, n * a);
    } else {
        len += sprintf(buf + len, "%d\ncrystal\n", num);
    }
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {
        for (int s = 0; s < num_sites; ++s) {
            len += sprintf(buf + len, "%s %.4f %.4f %.4f\n", sites[s].element, (i + sites[s].x) * a, (j + sites[s].y) * a, (k + sites[s].z) * a);
        }
    }
    return (str_t){ buf, (size_t)len };
}

static const crystal_site_t fcc_cu[] = { {"Cu", 0, 0, 0}, {"Cu", 0.5f, 0.5f, 0}, {"Cu", 0.5f, 0, 0.5f}, {"Cu", 0, 0.5f, 0.5f} };

static size_t count_kind(const md_system_t* sys, md_bond_kind_t kind) {
    size_t n = 0;
    for (size_t i = 0; i < sys->bond.count; ++i) n += md_bond_kind(sys->bond.flags[i]) == kind;
    return n;
}

static size_t count_flag(const md_system_t* sys, md_atom_flags_t flag) {
    size_t n = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) n += (md_atom_flags(&sys->atom, i) & flag) != 0;
    return n;
}

// A copper crystal: every atom is metallic, bonded to its 12 nearest neighbours (fewer at the faces, edges and corners
// of a crystal which is not periodic: a corner has 3) and to nothing else, and its bonds are of no ring nor structure
UTEST(bonds, metallic_crystal) {
    for (int periodic = 0; periodic < 2; ++periodic) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        const int n = 4;
        const float a = 3.615f;
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &st, crystal_xyz(fcc_cu, 4, n, a, periodic, alloc), MD_XYZ_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
        const size_t num = sys.atom.count;
        EXPECT_EQ(num, count_flag(&sys, MD_ATOM_FLAG_METALLIC));
        EXPECT_EQ(0, (int)count_kind(&sys, MD_BOND_KIND_COVALENT));

        // The nearest neighbour pairs, a / sqrt(2) apart
        const float nn = a / sqrtf(2.0f);
        size_t expected = 0;
        for (size_t i = 0; i < num; ++i) {
            for (size_t j = i + 1; j < num; ++j) {
                vec3_t d = vec3_sub(st.xyz[j], st.xyz[i]);
                md_util_min_image_vec3(&d, 1, &st.unitcell);
                expected += vec3_length(d) < 1.1f * nn;
            }
        }
        EXPECT_EQ(periodic ? num * 12 / 2 : expected, count_kind(&sys, MD_BOND_KIND_METALLIC));
        EXPECT_EQ(0, (int)md_ring_data_count(&sys.ring));
        md_vm_arena_destroy(alloc);
    }
}

// Ionic crystals are no metal, though each cation has 12 cations as close as the atoms of a metal are: NaCl and MgO
// (rock salt), with 6 anions to each cation
UTEST(bonds, metallic_not_ionic) {
    const struct { const char* cation; const char* anion; float a; } salts[] = { {"Na", "Cl", 5.640f}, {"Mg", "O", 4.212f} };
    for (size_t s = 0; s < ARRAY_SIZE(salts); ++s) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        const crystal_site_t rock_salt[] = {
            {salts[s].cation, 0, 0, 0}, {salts[s].cation, 0.5f, 0.5f, 0}, {salts[s].cation, 0.5f, 0, 0.5f}, {salts[s].cation, 0, 0.5f, 0.5f},
            {salts[s].anion, 0.5f, 0, 0}, {salts[s].anion, 0, 0.5f, 0}, {salts[s].anion, 0, 0, 0.5f}, {salts[s].anion, 0.5f, 0.5f, 0.5f},
        };
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &st, crystal_xyz(rock_salt, 8, 3, salts[s].a, true, alloc), MD_XYZ_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
        EXPECT_EQ(0, (int)count_flag(&sys, MD_ATOM_FLAG_METALLIC));
        EXPECT_EQ(0, (int)count_kind(&sys, MD_BOND_KIND_METALLIC));
        md_vm_arena_destroy(alloc);
    }
}

// Deterministic normal deviates (Box-Muller over a 64-bit LCG)
static float gauss(uint64_t* state) {
    float u[2];
    for (int k = 0; k < 2; ++k) {
        *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
        u[k] = ((float)(*state >> 40) + 0.5f) / (float)(1 << 24);
    }
    return sqrtf(-2.0f * logf(u[0])) * cosf(6.2831853f * u[1]);
}

// A slab of two (111) layers of gold, periodic in the plane, its atoms displaced as by heat (0.1 Å in each direction):
// the distances of the first shell spread across the cutoff of the covalent radii (2.72 Å). Over it a 4-bromophenyl,
// its free carbon 2.0 Å over a gold atom, and a bromomethane, its bromine 2.7 Å over another: within reach of a
// coordination each. The gold is metallic and has its 9 neighbours as metallic bonds, the molecules keep their bonds
// and have none to the gold; per frame the carbon forms its coordination, the bromine (of an organohalide) none.
UTEST(bonds, metallic_slab) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const int nx = 8, ny = 5;
    const float d = 2.884f;
    const float h = d * 0.8660254f;
    const float lz = 40.0f;
    const float z0 = 10.0f;
    const float dz = d * 0.8164966f;
    const int num_au = 2 * 2 * nx * ny;
    const int num = num_au + 11 + 5;
    char* buf = md_alloc(alloc, 256 + 64 * (size_t)num);
    int len = sprintf(buf, "%d\nLattice=\"%.4f 0 0 0 %.4f 0 0 0 %.4f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n", num, nx * d, ny * 2 * h, lz);

    uint64_t rng = 12345;
    vec3_t au[2 * 2 * 8 * 5];
    int n = 0;
    for (int layer = 0; layer < 2; ++layer) {
        for (int i = 0; i < nx; ++i) for (int j = 0; j < ny; ++j) for (int b = 0; b < 2; ++b) {
            const float x = i * d + b * 0.5f * d + layer * 0.5f * d;
            const float y = j * 2 * h + b * h + layer * h / 3.0f;
            au[n++] = vec3_set(x + 0.1f * gauss(&rng), y + 0.1f * gauss(&rng), z0 + layer * dz + 0.1f * gauss(&rng));
        }
    }
    for (int i = 0; i < num_au; ++i) len += sprintf(buf + len, "Au %.4f %.4f %.4f\n", au[i].x, au[i].y, au[i].z);

    // 4-bromophenyl upright in the xz-plane, C1 (the free carbon) 2.0 Å over the gold atom au_c of the top layer
    const int au_c = num_au / 2 + 2 * (1 * ny + 1);
    const vec3_t c1 = vec3_add(au[au_c], vec3_set(0, 0, 2.0f));
    const float r = 1.39f;
    const vec3_t ring[6] = {
        c1,
        vec3_add(c1, vec3_set( r * 0.8660254f, 0, 0.5f * r)), vec3_add(c1, vec3_set( r * 0.8660254f, 0, 1.5f * r)),
        vec3_add(c1, vec3_set(0, 0, 2 * r)),
        vec3_add(c1, vec3_set(-r * 0.8660254f, 0, 1.5f * r)), vec3_add(c1, vec3_set(-r * 0.8660254f, 0, 0.5f * r)),
    };
    const vec3_t center = vec3_add(c1, vec3_set(0, 0, r));
    for (int k = 0; k < 6; ++k) len += sprintf(buf + len, "C %.4f %.4f %.4f\n", ring[k].x, ring[k].y, ring[k].z);
    for (int k = 1; k < 6; ++k) {
        const vec3_t out = vec3_normalize(vec3_sub(ring[k], center));
        const vec3_t p = vec3_add(ring[k], vec3_mul1(out, k == 3 ? 1.90f : 1.08f));
        len += sprintf(buf + len, "%s %.4f %.4f %.4f\n", k == 3 ? "Br" : "H", p.x, p.y, p.z);
    }

    // Bromomethane, its bromine 2.7 Å over the gold atom au_br (far from au_c) and its carbon above
    const int au_br = num_au / 2 + 2 * (5 * ny + 3);
    const vec3_t br = vec3_add(au[au_br], vec3_set(0, 0, 2.7f));
    const vec3_t c = vec3_add(br, vec3_set(0, 0, 1.94f));
    len += sprintf(buf + len, "Br %.4f %.4f %.4f\nC %.4f %.4f %.4f\n", br.x, br.y, br.z, c.x, c.y, c.z);
    for (int k = 0; k < 3; ++k) {
        const float phi = k * 2.0943951f;
        const vec3_t p = vec3_add(c, vec3_mul1(vec3_set(0.943f * cosf(phi), 0.943f * sinf(phi), 0.333f), 1.09f));
        len += sprintf(buf + len, "H %.4f %.4f %.4f\n", p.x, p.y, p.z);
    }

    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &st, (str_t){ buf, (size_t)len }, MD_XYZ_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    ASSERT_EQ(num, (int)sys.atom.count);

    for (int i = 0; i < num; ++i) {
        EXPECT_EQ(i < num_au, (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_METALLIC) != 0);
    }

    // The metallic bonds are the pairs of the first shell, 9 to each atom; the rest are the bonds of the molecules
    size_t shell = 0;
    for (int i = 0; i < num_au; ++i) {
        for (int j = i + 1; j < num_au; ++j) {
            vec3_t v = vec3_sub(st.xyz[j], st.xyz[i]);
            md_util_min_image_vec3(&v, 1, &st.unitcell);
            shell += vec3_length(v) < 3.5f;
        }
    }
    EXPECT_EQ((size_t)num_au * 9 / 2, shell);
    EXPECT_EQ(shell, count_kind(&sys, MD_BOND_KIND_METALLIC));
    EXPECT_EQ(15, (int)(sys.bond.count - count_kind(&sys, MD_BOND_KIND_METALLIC)));
    size_t c_br = 0;
    for (size_t b = 0; b < sys.bond.count; ++b) {
        if (md_bond_kind(sys.bond.flags[b]) == MD_BOND_KIND_METALLIC) continue;
        const int za = md_atom_atomic_number(&sys.atom, sys.bond.pairs[b].idx[0]);
        const int zb = md_atom_atomic_number(&sys.atom, sys.bond.pairs[b].idx[1]);
        EXPECT_TRUE(za != 79 && zb != 79);
        c_br += (za == 6 && zb == 35) || (za == 35 && zb == 6);
    }
    EXPECT_EQ(2, (int)c_br);
    EXPECT_EQ(1, (int)md_ring_data_count(&sys.ring));

    // The frame's own contacts: the free carbon on its gold atom, nothing else
    const md_bond_delta_desc_t desc = { .tolerance = MD_BOND_DELTA_TOLERANCE, .breakable = MD_BOND_DELTA_BREAKABLE };
    md_bond_delta_t delta = {0};
    ASSERT_TRUE(md_util_bond_delta(&delta, &sys, &st, &desc, alloc));
    EXPECT_TRUE(md_bitset_empty(delta.broken));
    ASSERT_EQ(1, (int)delta.num_formed);
    const int fa = (int)delta.formed_pair[0].idx[0];
    const int fb = (int)delta.formed_pair[0].idx[1];
    EXPECT_TRUE((fa == au_c && fb == num_au) || (fb == au_c && fa == num_au));
    EXPECT_EQ(MD_BOND_KIND_COORDINATE, md_bond_kind(delta.formed_flags[0]));
    md_bond_delta_free(&delta, alloc);
    md_vm_arena_destroy(alloc);
}

// Two DBBA molecules (C28H16Br2) moving over a slab of two layers of gold (448 atoms), whose first shell (2.3-3.3 Å) a
// cutoff on the covalent radii (2.72 Å) would cut through: every gold atom is metallic, bonded to the gold atoms of
// its shell alone, the same in every frame, and the molecules have no bonds to the gold but those formed in a frame.
// Their bromine stays on its carbon throughout, so it coordinates no gold. The trajectory is not of the test data
// repository: the test runs where it has been put there.
UTEST(bonds, metallic_slab_trajectory) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
    if (!md_path_is_valid(path)) UTEST_SKIP("40-40-2-ddba-dyna.xmol not in the test data");
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_xyz_data_t data = {0};
    ASSERT_TRUE(md_xyz_data_parse_file(&data, path, alloc));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_data(&sys, &st, &data, MD_XYZ_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const size_t num = sys.atom.count;
    ASSERT_EQ(540, (int)num);

    size_t num_gold = 0;
    for (size_t i = 0; i < num; ++i) {
        const bool gold = md_atom_atomic_number(&sys.atom, i) == 79;
        num_gold += gold;
        EXPECT_EQ(gold, (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_METALLIC) != 0);
    }
    ASSERT_EQ(448, (int)num_gold);

    // No bond of the graph touches the gold; the molecules have theirs, the C-Br among them
    size_t c_br = 0;
    for (size_t b = 0; b < sys.bond.count; ++b) {
        const int za = md_atom_atomic_number(&sys.atom, sys.bond.pairs[b].idx[0]);
        const int zb = md_atom_atomic_number(&sys.atom, sys.bond.pairs[b].idx[1]);
        if (md_bond_kind(sys.bond.flags[b]) == MD_BOND_KIND_METALLIC) {
            EXPECT_TRUE(za == 79 && zb == 79);
            continue;
        }
        EXPECT_TRUE(za != 79 && zb != 79);
        c_br += (za == 6 && zb == 35) || (za == 35 && zb == 6);
    }
    EXPECT_EQ(4, (int)c_br);
    EXPECT_EQ(102, (int)(sys.bond.count - count_kind(&sys, MD_BOND_KIND_METALLIC)));

    // The metallic bonds are the pairs of the first shell: all of them within 3.25 Å, none beyond the gap at 3.3-3.4 Å
    size_t within = 0;
    for (size_t i = 0; i < num; ++i) {
        for (size_t j = i + 1; j < num; ++j) {
            if (md_atom_atomic_number(&sys.atom, i) != 79 || md_atom_atomic_number(&sys.atom, j) != 79) continue;
            within += vec3_distance(st.xyz[i], st.xyz[j]) < 3.25f;
        }
    }
    const size_t num_metallic = count_kind(&sys, MD_BOND_KIND_METALLIC);
    EXPECT_GE(num_metallic, within);
    for (size_t b = 0; b < sys.bond.count; ++b) {
        if (md_bond_kind(sys.bond.flags[b]) != MD_BOND_KIND_METALLIC) continue;
        EXPECT_LT(vec3_distance(st.xyz[sys.bond.pairs[b].idx[0]], st.xyz[sys.bond.pairs[b].idx[1]]), 3.45f);
    }

    // No ring of gold, and the molecules are two structures of their own
    for (size_t r = 0; r < md_ring_data_count(&sys.ring); ++r) {
        EXPECT_NE(79, (int)md_atom_atomic_number(&sys.atom, md_ring_atoms(&sys.ring, r)[0]));
    }

    // In every frame: the same metallic bonds inferred (within a percent), and the bonds formed to the gold are the
    // contacts of the molecules with it, never a bond between two gold atoms nor the bromine on the gold
    size_t frames_with_contacts = 0;
    for (size_t f = 0; f < data.num_models; ++f) {
        md_system_state_t fs = { .alloc = alloc };
        md_system_state_init(&fs, num);
        for (size_t i = 0; i < num; ++i) {
            const md_xyz_coordinate_t c = data.coordinates[data.models[f].beg_coord_index + i];
            fs.xyz[i] = vec3_set(c.x, c.y, c.z);
        }
        md_bond_data_t bonds = {0};
        md_util_infer_covalent_bonds(&bonds, &fs, &sys, alloc);
        size_t metallic = 0;
        for (size_t b = 0; b < bonds.count; ++b) metallic += md_bond_kind(bonds.flags[b]) == MD_BOND_KIND_METALLIC;
        EXPECT_LT(fabs((double)metallic - (double)num_metallic), 0.01 * (double)num_metallic);

        const md_bond_delta_desc_t desc = { .tolerance = MD_BOND_DELTA_TOLERANCE, .breakable = MD_BOND_DELTA_BREAKABLE };
        md_bond_delta_t d = {0};
        ASSERT_TRUE(md_util_bond_delta(&d, &sys, &fs, &desc, alloc));
        EXPECT_TRUE(md_bitset_empty(d.broken));
        size_t contacts = 0;
        for (size_t k = 0; k < d.num_formed; ++k) {
            const int za = md_atom_atomic_number(&sys.atom, d.formed_pair[k].idx[0]);
            const int zb = md_atom_atomic_number(&sys.atom, d.formed_pair[k].idx[1]);
            EXPECT_FALSE(za == 79 && zb == 79);
            EXPECT_FALSE(za == 35 || zb == 35);
            if (za == 79 || zb == 79) {
                EXPECT_EQ(MD_BOND_KIND_COORDINATE, md_bond_kind(d.formed_flags[k]));
                contacts += 1;
            }
        }
        frames_with_contacts += contacts > 0;
        md_bond_delta_free(&d, alloc);
    }
    EXPECT_GT(frames_with_contacts, 0);
    md_vm_arena_destroy(alloc);
}
