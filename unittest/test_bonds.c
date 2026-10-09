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
