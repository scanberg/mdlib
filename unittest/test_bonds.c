#include "utest.h"

// Covalent and coordination bond inference: against force field topologies, on structures that clash, at metal sites

#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_mmcif.h>
#include <md_tpr.h>
#include <md_lammps.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>

#include <stdlib.h>
#include <string.h>

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
    EXPECT_TRUE(f & MD_BOND_FLAG_COORDINATE);
    EXPECT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin(f));
    // Na is a free ion: nothing
    EXPECT_EQ(0, (int)md_bond_conn_count(&sys.bond, 5));
    // Water is never bonded to a metal
    EXPECT_FALSE(has_bond(&sys, 4, 6, NULL));
    // Cu-Cu: a metal-metal bond, which is not coordination
    ASSERT_TRUE(has_bond(&sys, 9, 10, &f));
    EXPECT_FALSE(f & MD_BOND_FLAG_COORDINATE);

    // A bond between a metal and a non-metal is coordination whatever its origin: the Mg-O bond as a topology gives it
    for (size_t b = 0; b < sys.bond.count; ++b) {
        sys.bond.flags[b] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
    }
    EXPECT_EQ(1, (int)md_util_system_infer_coordination(&sys));
    ASSERT_TRUE(has_bond(&sys, 2, 4, &f));
    EXPECT_TRUE(f & MD_BOND_FLAG_COORDINATE);
    EXPECT_EQ(MD_BOND_ORIGIN_TOPOLOGY, md_bond_origin(f));
    ASSERT_TRUE(has_bond(&sys, 9, 10, &f));
    EXPECT_FALSE(f & MD_BOND_FLAG_COORDINATE);
    ASSERT_TRUE(has_bond(&sys, 1, 2, &f));
    EXPECT_FALSE(f & MD_BOND_FLAG_COORDINATE);
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
            if ((md_bond_iter_bond_flags(&it) & MD_BOND_FLAG_COORDINATE) && md_atom_atomic_number(&sys.atom, md_bond_iter_atom_index(&it)) == 16) num_s += 1;
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
    const md_bond_flags_t hbond = user | MD_BOND_FLAG_NONCOVALENT | MD_BOND_FLAG_HYDROGEN;
    const md_bond_flags_t ionic = user | MD_BOND_FLAG_NONCOVALENT | MD_BOND_FLAG_IONIC;
    ASSERT_FALSE(has_bond(&sys, 0, n - 1, NULL));
    ASSERT_FALSE(has_bond(&sys, 1, n - 2, NULL));
    md_system_bond_insert(&sys, 0, n - 1, hbond);
    md_system_bond_insert(&sys, 1, n - 2, ionic);
    EXPECT_EQ(graph, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);
    EXPECT_EQ(graph + 2, md_array_size(sys.bond.pairs));
    EXPECT_EQ(graph + 2, md_array_size(sys.bond.flags));

    // An atom pair has one bond: neither range takes a second one
    md_system_bond_insert(&sys, n - 1, 0, ionic);
    md_system_bond_insert(&sys, sys.bond.pairs[0].idx[0], sys.bond.pairs[0].idx[1], hbond);
    EXPECT_EQ(graph, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);

    // A bond of the graph inserted goes before them
    md_system_bond_insert(&sys, 2, n - 3, user);
    EXPECT_EQ(graph + 1, sys.bond.count);
    EXPECT_EQ(2, (int)sys.bond.noncovalent_count);
    for (size_t i = 0; i < sys.bond.count; ++i) EXPECT_FALSE(sys.bond.flags[i] & MD_BOND_FLAG_NONCOVALENT);
    for (size_t i = sys.bond.count; i < sys.bond.count + sys.bond.noncovalent_count; ++i) EXPECT_TRUE(sys.bond.flags[i] & MD_BOND_FLAG_NONCOVALENT);
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
    for (size_t i = 0; i < sys.bond.count; ++i) EXPECT_FALSE(sys.bond.flags[i] & MD_BOND_FLAG_NONCOVALENT);
    for (size_t i = sys.bond.count; i < sys.bond.count + sys.bond.noncovalent_count; ++i) EXPECT_TRUE(sys.bond.flags[i] & MD_BOND_FLAG_NONCOVALENT);
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
    EXPECT_FALSE(sys.bond.flags[b] & MD_BOND_FLAG_COORDINATE);
    // Without it, the kind is perceived
    sys.bond.flags[b] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER);
    EXPECT_EQ(1, (int)md_util_system_infer_coordination(&sys));
    EXPECT_TRUE(sys.bond.flags[b] & MD_BOND_FLAG_COORDINATE);
    md_vm_arena_destroy(alloc);
}

// Every flag combination of a kind is of the kinds the name table says, and of no other
UTEST(bonds, kind_names) {
    size_t count = 0;
    const md_bond_kind_name_t* names = md_bond_kind_names(&count);
    ASSERT_EQ(7, (int)count);
    const struct { md_bond_flags_t flags; const char* kinds; } cases[] = {
        { MD_BOND_FLAG_NONE,                                                    "covalent" },
        { MD_BOND_FLAG_AROMATIC | 2,                                            "covalent" },
        { MD_BOND_FLAG_COORDINATE,                                              "coordinate" },
        // A spring is of no other kind
        { MD_BOND_FLAG_ELASTIC,                                                 "elastic" },
        { MD_BOND_FLAG_NONCOVALENT,                                             "noncovalent" },
        { MD_BOND_FLAG_NONCOVALENT | MD_BOND_FLAG_HYDROGEN,                     "noncovalent hydrogen" },
        { MD_BOND_FLAG_NONCOVALENT | MD_BOND_FLAG_HYDROGEN | MD_BOND_FLAG_IONIC, "noncovalent hydrogen ionic" },
        { MD_BOND_FLAG_NONCOVALENT | MD_BOND_FLAG_SIGMA_HOLE,                   "noncovalent sigma_hole" },
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
    EXPECT_FALSE(d.formed_flags[0] & MD_BOND_FLAG_KIND_MASK);

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
    md_system_bond_insert(&f.sys, 1, 3, md_bond_flags_set_origin(MD_BOND_FLAG_NONCOVALENT | MD_BOND_FLAG_HYDROGEN | MD_BOND_FLAG_KIND_GIVEN, MD_BOND_ORIGIN_USER));
    f.st.xyz[1] = vec3_set(1.80f, 0.0f, 0.0f);
    md_bond_delta_t d = {0};
    ASSERT_TRUE(md_util_bond_delta(&d, &f.sys, &f.st, &delta_default, f.alloc));
    EXPECT_EQ(1, (int)md_bitset_count(d.broken));
    EXPECT_EQ(0, (int)d.num_formed);
    md_vm_arena_destroy(f.alloc);
}
