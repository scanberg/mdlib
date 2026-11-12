#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <stdio.h>
#include <math.h>
#include <md_system.h>
#include <md_run.h>
#include <md_infer.h>
#include <md_attributes.h>
#include <core/md_arena_allocator.h>
#include <core/md_file.h>
#include <core/md_vec_math.h>

#include "run_check.h"

// The file's atoms, its cell and the bonds of its CONECT records, as the system has them
UTEST(pdb, ordinary) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE, NULL));
    EXPECT_EQ((size_t)9084, sys.atom.count);
    EXPECT_EQ((size_t)9084, st.num_atoms);
    // Its CRYST1 is 1 x 1 x 1 at right angles, what writers put there for a structure of no cell
    EXPECT_EQ(0u, st.unitcell.flags);
    EXPECT_EQ((size_t)18, sys.bond.count);
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, tryptophan) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), MD_PDB_OPTION_NONE, NULL));
    EXPECT_EQ((size_t)28, sys.atom.count);
    md_vm_arena_destroy(alloc);
}

// A MODEL record without its ENDMDL: all of the atoms after it
UTEST(pdb, unmatched_model_entry) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"), MD_PDB_OPTION_NONE, NULL));
    EXPECT_EQ((size_t)14738, sys.atom.count);
    EXPECT_EQ((size_t)0, sys.bond.count);
    EXPECT_NE(0u, st.unitcell.flags);
    md_vm_arena_destroy(alloc);
}

// The models of a file of several: the first is the system, and the run's frames begin at the MODEL records
UTEST(pdb, models) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t path = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, path, MD_PDB_OPTION_NONE, NULL));
    EXPECT_EQ((size_t)153, sys.atom.count);
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, path, MD_STR_LIT("run/models"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* pos = md_attributes_find(&sys.attributes, MD_STR_LIT("run/models/atom/position"));
    ASSERT_TRUE(pos != NULL);
    md_run_source_t src;
    ASSERT_TRUE(md_run_source(&src, &sys.attributes, pos, MD_STR_LIT("atom/position")));
    ASSERT_EQ(38u, (uint32_t)pos->format.shape[0]);

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (uint32_t i = 0; i < 38; ++i) {
        char data[6] = {0};
        md_file_seek(file, src.offset[i], MD_FILE_BEG);
        md_file_read(file, data, 5);
        EXPECT_EQ(strncmp(data, "MODEL", 5), 0);
    }
    md_file_close(&file);
    md_vm_arena_destroy(alloc);
}

// REMARK 350 of 1k4r: biomolecule 1 is its three chains by 60 transforms, the shell of the virion
UTEST(pdb, assemblies) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE, NULL));
    const md_assembly_data_t* a = &sys.assembly;
    ASSERT_EQ(1, a->count);
    EXPECT_STREQ("1", md_str_ptr(md_assembly_id(a, 0)));
    ASSERT_EQ(1, a->group_count);
    EXPECT_EQ(60, a->oper_count);
    const md_urange_t ir = md_assembly_group_inst_range(a, 0);
    EXPECT_EQ(3, ir.end - ir.beg);
    EXPECT_EQ(60, md_assembly_copy_count(a, 0));
    const md_float4x4_t I = md_float4x4_ident();
    EXPECT_EQ(0, md_memcmp(&I, &a->oper_transform[a->oper[0]], sizeof(md_float4x4_t)));
    md_vm_arena_destroy(alloc);
}

// Two groups of chains in one biomolecule, each with transforms of its own: a chain is copied by the transforms of its
// group only, and the identity both have is one copy
UTEST(pdb, assembly_groups) {
    static const char pdb[] =
        "REMARK 350 BIOMOLECULE: 1\n"
        "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n"
        "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n"
        "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\n"
        "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n"
        "REMARK 350   BIOMT1   2  0.000000 -1.000000  0.000000       10.00000\n"
        "REMARK 350   BIOMT2   2  1.000000  0.000000  0.000000        0.00000\n"
        "REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000\n"
        "REMARK 350 APPLY THE FOLLOWING TO CHAINS: B,\n"
        "REMARK 350                    AND CHAINS: C\n"
        "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n"
        "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\n"
        "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  N   ALA B   1       0.000  20.000   0.000  1.00  0.00           N\n"
        "ATOM      4  CA  ALA B   1       1.458  20.000   0.000  1.00  0.00           C\n"
        "ATOM      5  N   ALA C   1       0.000  40.000   0.000  1.00  0.00           N\n"
        "ATOM      6  CA  ALA C   1       1.458  40.000   0.000  1.00  0.00           C\n"
        "END\n";
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, (md_str_t){pdb, sizeof(pdb) - 1}, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ(3, sys.instance.count);
    const md_assembly_data_t* a = &sys.assembly;
    ASSERT_EQ(1, a->count);
    // The identity copies A, B and C; the other only A
    ASSERT_EQ(2, a->oper_count);
    ASSERT_EQ(2, a->group_count);
    md_urange_t ir = md_assembly_group_inst_range(a, 0);
    md_urange_t orr = md_assembly_group_oper_range(a, 0);
    EXPECT_EQ(3, ir.end - ir.beg);
    EXPECT_EQ(1, orr.end - orr.beg);
    EXPECT_EQ(0, a->oper[orr.beg]);
    ir = md_assembly_group_inst_range(a, 1);
    orr = md_assembly_group_oper_range(a, 1);
    EXPECT_EQ(1, ir.end - ir.beg);
    EXPECT_EQ(1, orr.end - orr.beg);
    EXPECT_EQ(1, a->oper[orr.beg]);
    // A quarter turn about z, by rows as the records give it: x' = -y + 10, y' = x
    const md_float4_t q = md_float4x4_mul_float4(a->oper_transform[1], (md_float4_t){1, 2, 3, 1});
    EXPECT_NEAR(8.0f, q.x, 1e-6f);
    EXPECT_NEAR(1.0f, q.y, 1e-6f);
    EXPECT_NEAR(3.0f, q.z, 1e-6f);
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, create_system) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = { .alloc = arena };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)9084, sys.atom.count);

    EXPECT_EQ(1185, sys.residue.count);
    // The chains of the file are its units, and their strands its parts: three of one sequence
    EXPECT_EQ(3, sys.instance.count);
    // The parts are found by the bonds and the kinds of the residues, which the load leaves to inference
    ASSERT_TRUE(md_system_infer(&sys, &sys_state, MD_SYSTEM_INFER_BOND_BIT | MD_SYSTEM_INFER_RESIDUE_BIT | MD_SYSTEM_INFER_PART_BIT, NULL));
    EXPECT_EQ(3, sys.part.count);
    EXPECT_EQ(1, md_system_part_type_count(&sys));
    md_vm_arena_destroy(arena);
}

static size_t count_bonds_of_origin(const md_system_t* sys, md_bond_origin_t origin) {
    size_t count = 0;
    for (size_t i = 0; i < sys->bond.count; ++i) count += md_bond_origin(sys->bond.flags[i]) == origin;
    return count;
}

// The CONECT records are bonds of origin MD_BOND_ORIGIN_FILE, each once though both its atoms list it. Inference adds the
// bonds they leave out, never one of theirs again, and keeps them when it runs again. 1k4r lists its 18 disulfides.
// The cell as some tools write it: CRYST with a blank for the 1 and all angles 0 (the Conformational Analysis Tools of
// CHARMM-GUI). A rectangular box: the bonds of a residue which lies across its boundary are found through it.
UTEST(pdb, cryst_record_without_its_1) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    static const char pdb[] =
        "CRYST   246.108  246.108  398.534   0.00   0.00   0.00\n"
        "ATOM      1  N   ALA A   1     245.000  10.000  10.000  0.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       0.300  10.300  10.000  0.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       1.500  11.200  10.000  0.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.500  12.400  10.000  0.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       0.300   9.000  11.000  0.00  0.00           C\n"
        "END\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(pdb), MD_PDB_OPTION_NONE, NULL));
    EXPECT_NEAR(246.108, state.unitcell.x, 1e-3);
    EXPECT_NEAR(246.108, state.unitcell.y, 1e-3);
    EXPECT_NEAR(398.534, state.unitcell.z, 1e-3);
    EXPECT_NEAR(0.0, state.unitcell.xy, 1e-9);
    EXPECT_NEAR(0.0, state.unitcell.xz, 1e-9);
    EXPECT_NEAR(0.0, state.unitcell.yz, 1e-9);
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_GE(md_bond_find(&sys.bond, 0, 1), 0);    // N-CA, across the boundary in x
    EXPECT_EQ((int)MD_RESIDUE_KIND_AMINO_ACID, (int)md_system_residue_kind(&sys, 0));
    md_arena_allocator_destroy(alloc);
}

// Past 99999 atoms VMD writes the serials in hexadecimal (186a0 for 100000), and '*****' past FFFFF: not hybrid-36,
// which begins a number past 99999 with a letter. They are not numbers of the format: the atoms are read in their
// order, and no CONECT record can name them. Hexadecimal without a letter (18700) is decimal; A0000 is 100000.
UTEST(pdb, serials_which_are_not_numbers) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    static const char pdb[] =
        "ATOM  99999  OW  SOL X3333       0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  186a0  OW  SOL X3334       3.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  18700  OW  SOL X3335       6.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  *****  OW  SOL X3336       9.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  A0000  OW  SOL X3337      12.000   0.000   0.000  1.00  0.00           O\n"
        "CONECT99999A0000\n"
        "CONECT1870099999\n"
        "CONECT186a099999\n"
        "CONECT*****99999\n"
        "END\n";
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, md_str_from_cstr(pdb), MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)5, sys.atom.count);
    EXPECT_NEAR(9.0f, st.xyz[3].x, 1e-6f);
    EXPECT_EQ((size_t)2, sys.bond.count);
    EXPECT_NE(-1, md_bond_find(&sys.bond, 0, 4));   // 99999 and A0000
    EXPECT_NE(-1, md_bond_find(&sys.bond, 0, 2));   // 18700 and 99999
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, conect_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE, NULL));
    // The load gives the file's bonds, and those alone
    EXPECT_EQ((size_t)18, sys.bond.count);
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_BOND_BIT, NULL));
    EXPECT_EQ((size_t)18, count_bonds_of_origin(&sys, MD_BOND_ORIGIN_FILE));
    for (size_t i = 0; i < sys.bond.count; ++i) {
        if (md_bond_origin(sys.bond.flags[i]) != MD_BOND_ORIGIN_FILE) continue;
        EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, sys.bond.pairs[i].idx[0]), MD_STR_LIT("SG")));
        EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, sys.bond.pairs[i].idx[1]), MD_STR_LIT("SG")));
    }
    // Atoms 23 and 226 (Cys 3 and Cys 30): one bond, from the file
    const md_bond_idx_t b = md_bond_find(&sys.bond, 22, 225);
    ASSERT_NE(-1, b);
    EXPECT_EQ(MD_BOND_ORIGIN_FILE, md_bond_origin(sys.bond.flags[b]));
    EXPECT_GT(count_bonds_of_origin(&sys, MD_BOND_ORIGIN_INFERRED), (size_t)9000);

    const size_t num_bonds = sys.bond.count;
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    EXPECT_EQ((size_t)18, count_bonds_of_origin(&sys, MD_BOND_ORIGIN_FILE));
    EXPECT_EQ(num_bonds, sys.bond.count);
    for (size_t i = 0; i < sys.bond.count; ++i) {
        const md_atom_pair_t p = sys.bond.pairs[i];
        if (md_bond_origin(sys.bond.flags[i]) == MD_BOND_ORIGIN_FILE) {
            EXPECT_EQ((md_bond_idx_t)i, md_bond_find(&sys.bond, p.idx[0], p.idx[1]));
        }
    }
    md_vm_arena_destroy(alloc);
}

// A record of more than four bonded atoms goes on over several CONECT lines, and the serials of atoms past 99999 are
// hybrid-36 (A0000 is 100000)
UTEST(pdb, conect_continued_and_hybrid36) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00          FE\n"
        "HETATM    2  N1  HEM A   1       2.000   0.000   0.000  1.00  0.00           N\n"
        "HETATM    3  N2  HEM A   1      -2.000   0.000   0.000  1.00  0.00           N\n"
        "HETATM    4  N3  HEM A   1       0.000   2.000   0.000  1.00  0.00           N\n"
        "HETATM    5  N4  HEM A   1       0.000  -2.000   0.000  1.00  0.00           N\n"
        "HETATM    6  O1  HEM A   1       0.000   0.000   2.000  1.00  0.00           O\n"
        "HETATM    7  O2  HEM A   1       0.000   0.000  -2.000  1.00  0.00           O\n"
        "HETATM    8  O3  HEM A   1       1.500   1.500   0.000  1.00  0.00           O\n"
        "HETATM    9  O4  HEM A   1      -1.500   1.500   0.000  1.00  0.00           O\n"
        "HETATM99999  C1  LIG B   1      10.000   0.000   0.000  1.00  0.00           C\n"
        "HETATMA0000  C2  LIG B   1      11.500   0.000   0.000  1.00  0.00           C\n"
        "CONECT    1    2    3    4    5\n"
        "CONECT    1    6    7    8\n"
        "CONECT    1    9\n"
        "CONECT99999A0000\n");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)11, sys.atom.count);
    EXPECT_EQ((size_t)9, sys.bond.count);
    for (int k = 1; k < 9; ++k) {
        EXPECT_NE(-1, md_bond_find(&sys.bond, 0, k));
    }
    EXPECT_NE(-1, md_bond_find(&sys.bond, 9, 10));
    md_vm_arena_destroy(alloc);
}

// CONECT records refer to atoms by their serials however the file numbers them (not counting up from the first, the
// same serial twice: the first atom of it), wherever in the file they are, and a bond given twice (each of its atoms
// names the other, or a writer repeats it for its order) is one
UTEST(pdb, conect_serials_out_of_order) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "HETATM   30  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM   10  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "CONECT   30   10\n"
        "HETATM   20  C3  LIG A   1       3.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM   20  C4  LIG A   1       4.500   0.000   0.000  1.00  0.00           C\n"
        "CONECT   10   30   20\n"
        "CONECT   20   10   10\n"
        "CONECT   30   99\n");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)4, sys.atom.count);
    ASSERT_EQ((size_t)2, sys.bond.count);
    // In the order the records first give them
    EXPECT_EQ(0, (int)sys.bond.pairs[0].idx[0]);
    EXPECT_EQ(1, (int)sys.bond.pairs[0].idx[1]);
    EXPECT_EQ(1, (int)sys.bond.pairs[1].idx[0]);
    EXPECT_EQ(2, (int)sys.bond.pairs[1].idx[1]);
    md_vm_arena_destroy(alloc);
}

// The bonds of a coarse grained structure (martinize2 writes CONECT for every bead) are between neighbouring atoms,
// each given by both of them: as many as there are, in time linear in them (a map of pairs keyed by their atoms took
// time quadratic in them, as neighbouring atoms hashed alike)
UTEST(pdb, conect_dense) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const size_t n = 50000;
    const size_t cap = n * 82 + (n - 1) * 34 + 1;
    char* buf = md_alloc(alloc, cap);
    size_t len = 0;
    for (size_t i = 0; i < n; ++i) {
        len += (size_t)snprintf(buf + len, cap - len, "ATOM  %5d  BB  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
            (int)(i + 1), (int)(i % 9999) + 1, (float)(i % 1000), 0.0f, 0.0f);
    }
    for (size_t i = 1; i < n; ++i) {
        len += (size_t)snprintf(buf + len, cap - len, "CONECT%5d%5d\nCONECT%5d%5d\n", (int)i, (int)(i + 1), (int)(i + 1), (int)i);
    }
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, (md_str_t){buf, len}, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ(n, sys.atom.count);
    ASSERT_EQ(n - 1, sys.bond.count);
    size_t wrong = 0;
    for (size_t i = 0; i + 1 < n; ++i) {
        wrong += sys.bond.pairs[i].idx[0] != (md_atom_idx_t)i || sys.bond.pairs[i].idx[1] != (md_atom_idx_t)(i + 1);
    }
    EXPECT_EQ((size_t)0, wrong);
    md_vm_arena_destroy(alloc);
}

// An atom name of four characters fills columns 13 to 16 and is kept whole: HN4' and HN4A are two atoms, not two HN4
UTEST(pdb, four_character_atom_names) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    const md_str_t text = MD_STR_LIT(
        "HETATM    1  N4' TD6 B   1      43.654  64.837  38.628  1.00  0.00           N\n"
        "HETATM    2 HN4' TD6 B   1      44.050  64.155  39.051  1.00  0.00           H\n"
        "HETATM    3 HN4A TD6 B   1      44.059  65.650  38.601  1.00  0.00           H\n"
        "HETATM    4  CA  TD6 B   1      45.000  65.000  38.000  1.00  0.00           C\n");
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)4, sys.atom.count);
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 0), MD_STR_LIT("N4'")));
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 1), MD_STR_LIT("HN4'")));
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 2), MD_STR_LIT("HN4A")));
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 3), MD_STR_LIT("CA")));
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, parse_nonexistent_file) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    EXPECT_FALSE(md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.pdb"), MD_PDB_OPTION_NONE, NULL));
    EXPECT_EQ((size_t)0, sys.atom.count);
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, parse_empty_path) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    EXPECT_FALSE(md_pdb_system_init_from_file(&sys, &st, (md_str_t){0}, MD_PDB_OPTION_NONE, NULL));
    md_vm_arena_destroy(alloc);
}
// ### RUN ###

#define PDB_RUN MD_STR_LIT("run/1ala")

// Recorded from the trajectory reader this replaced
static const run_ref_t pdb_refs[] = {
    { 0,  {3568.34002, 7394.96399, 3699.66098}, {23.7040005, 23, 21.3549995},         {25.0020008, 75.5550003, 23.5319996}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
    { 19, {3297.43299, 6833.70401, 3418.89901}, {22.4950008, 43.8180008, 11.6870003}, {7.09700012, 53.9360008, 23.4740009}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
    { 37, {3298.315, 6835.46697, 3419.76601},   {28.4710007, 57.6910019, 25.2150002}, {8.61200047, 35.6669998, 14.9329996}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
};

// One frame per model, time as model ordinals, and the file's one CRYST1 cell at every frame.
UTEST(pdb, run_matches_reference) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = {.alloc = arena};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &sys_state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL));
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, path, PDB_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, PDB_RUN, 38, 153, pdb_refs, MD_COUNTOF(pdb_refs));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, MD_STR_LIT("run/1ala/time"));
    const md_attribute_t* pos  = md_attributes_find(&sys.attributes, MD_STR_LIT("run/1ala/atom/position"));
    ASSERT_TRUE(time && pos);
    EXPECT_TRUE(md_unit_is_none(time->unit));
    EXPECT_EQ(9.0, ((const double*)time->data)[9]);

    // The first model is the structure's own coordinates.
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, sys.atom.count);
    ASSERT_TRUE(run_extract_one(&got, &sys, PDB_RUN, 0));
    EXPECT_EQ(0, md_memcmp(sys_state.xyz, got.xyz, sys.atom.count * sizeof(md_float3_t)));

    // One atom of one model: the text is read up to that atom and no further.
    float xyz[3];
    ASSERT_TRUE(run_extract_one(&got, &sys, PDB_RUN, 9));
    md_attribute_slice_t one = md_attribute_slice_2(9, 100);
    ASSERT_EQ(3u, md_attribute_extract_f32(xyz, 3, pos, one, md_unit_none()));
    EXPECT_EQ(got.xyz[100].x, xyz[0]);
    EXPECT_EQ(got.xyz[100].z, xyz[2]);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// The atoms of a model are those the system has: an alternate location other than the first is in neither, so the
// frames are of the system's atoms (they were not, when the models were read whole)
UTEST(pdb, run_alternate_locations) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
        "MODEL        1\n"
        "ATOM      1  N   SER A   1       1.000   1.000   1.000  1.00  0.00           N\n"
        "ATOM      2  CA ASER A   1       2.000   1.000   1.000  0.60  0.00           C\n"
        "ATOM      3  CA BSER A   1       2.100   1.100   1.100  0.40  0.00           C\n"
        "ATOM      4  C   SER A   1       3.000   1.000   1.000  1.00  0.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      1  N   SER A   1       1.500   1.000   1.000  1.00  0.00           N\n"
        "ATOM      2  CA ASER A   1       2.500   1.000   1.000  0.60  0.00           C\n"
        "ATOM      3  CA BSER A   1       2.600   1.100   1.100  0.40  0.00           C\n"
        "ATOM      4  C   SER A   1       3.500   1.000   1.000  1.00  0.00           C\n"
        "ENDMDL\n"
        "END\n");
    const char* path = "md_unittest_altloc_models.pdb";
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, md_str_from_cstr(path), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    ASSERT_EQ(text.len, md_file_write(file, text.ptr, text.len));
    md_file_close(&file);

    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, md_str_from_cstr(path), MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)3, sys.atom.count);
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, md_str_from_cstr(path), MD_STR_LIT("run/alt"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));

    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, sys.atom.count);
    ASSERT_TRUE(run_extract_one(&got, &sys, MD_STR_LIT("run/alt"), 0));
    EXPECT_EQ(0, md_memcmp(st.xyz, got.xyz, sys.atom.count * sizeof(md_float3_t)));
    ASSERT_TRUE(run_extract_one(&got, &sys, MD_STR_LIT("run/alt"), 1));
    EXPECT_EQ(1.5f, got.xyz[0].x);
    EXPECT_EQ(2.5f, got.xyz[1].x);     // the A location of the CA
    EXPECT_EQ(3.5f, got.xyz[2].x);     // the C, not the B location of the CA
    EXPECT_NE(0u, got.unitcell.flags);

    remove(path);
    md_vm_arena_destroy(arena);
}

// A file of one model is a structure, not a run.
UTEST(pdb, run_needs_several_models) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    EXPECT_FALSE(md_pdb_system_publish_run(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), PDB_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, MD_STR_LIT("run/1ala/time")) == NULL);
    md_vm_arena_destroy(arena);
}

// The columns read without a scan give what parse_float gives on them: every coordinate, occupancy
// and temperature factor of the reference files
#include <core/md_parse.h>

static size_t pdb_test_compare_with_columns(md_str_t path) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    if (!md_pdb_system_init_from_file(&sys, &st, path, MD_PDB_OPTION_NONE, NULL)) {
        md_vm_arena_destroy(alloc);
        return SIZE_MAX;
    }
    // A column of one value throughout is not published (md_attributes_publish_atom_column): NULL, and not compared
    const float* occupancy = (const float*)md_attribute_view(md_attributes_find(&sys.attributes, MD_STR_LIT("atom/occupancy")), MD_ATTRIBUTE_TYPE_F32, 1, 1);
    const float* b_factor  = (const float*)md_attribute_view(md_attributes_find(&sys.attributes, MD_STR_LIT("atom/b_factor")),  MD_ATTRIBUTE_TYPE_F32, 1, 1);
    md_str_t text = md_load_textfile(path, alloc);
    md_str_t rest = text, line;
    size_t i = 0, mismatches = 0;
    static const size_t cols[5][2] = { {30, 8}, {38, 8}, {46, 8}, {54, 6}, {60, 6} };
    // The atoms of the first model, of the first alternate location
    while (md_parse_next_line(&line, &rest) && i < sys.atom.count) {
        if (!(md_str_begins_with(line, MD_STR_LIT("ATOM")) || md_str_begins_with(line, MD_STR_LIT("HETATM")))) continue;
        if (line.len > 16 && line.ptr[16] != ' ' && line.ptr[16] != 'A') continue;
        const float got[5] = { st.xyz[i].x, st.xyz[i].y, st.xyz[i].z, occupancy ? occupancy[i] : 0.0f, b_factor ? b_factor[i] : 0.0f };
        i += 1;
        for (int k = 0; k < 5; ++k) {
            if (line.len < cols[k][0] + cols[k][1]) continue;
            if ((k == 3 && !occupancy) || (k == 4 && !b_factor)) continue;
            const float want = (float)parse_float(md_str_substr(line, cols[k][0], cols[k][1]));
            if (memcmp(&want, &got[k], sizeof(float)) != 0) mismatches += 1;
        }
    }
    if (i != sys.atom.count) mismatches += 1;
    md_vm_arena_destroy(alloc);
    return mismatches;
}

UTEST(pdb, columns_match_parse_float) {
    EXPECT_EQ(0, pdb_test_compare_with_columns(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
    EXPECT_EQ(0, pdb_test_compare_with_columns(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb")));
    EXPECT_EQ(0, pdb_test_compare_with_columns(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb")));
}

// ### CHAIN BREAKS AND INSERTION CODES ###
//
// The backbone (N, CA, C, O) of the first ten residues of chain A of 1k4r: a stretch of a real chain, to cut and number
// as a test needs
typedef struct {
    const char* name;
    float xyz[4][3];
} peptide_residue_t;

static const peptide_residue_t k4r_backbone[10] = {
    {"SER", {{-42.732f, -19.665f, 235.272f}, {-41.832f, -20.531f, 234.522f}, {-40.476f, -19.848f, 234.349f}, {-40.419f, -18.708f, 233.902f}}},
    {"ARG", {{-39.391f, -20.558f, 234.642f}, {-38.057f, -19.971f, 234.546f}, {-37.710f, -19.352f, 233.192f}, {-37.158f, -18.261f, 233.145f}}},
    {"CYS", {{-38.028f, -20.029f, 232.098f}, {-37.719f, -19.493f, 230.781f}, {-38.347f, -18.131f, 230.486f}, {-37.844f, -17.403f, 229.623f}}},
    {"THR", {{-39.415f, -17.754f, 231.195f}, {-40.030f, -16.449f, 230.952f}, {-39.151f, -15.306f, 231.495f}, {-39.421f, -14.127f, 231.234f}}},
    {"HIS", {{-38.116f, -15.655f, 232.265f}, {-37.208f, -14.673f, 232.855f}, {-36.184f, -14.111f, 231.867f}, {-35.542f, -13.086f, 232.156f}}},
    {"LEU", {{-35.940f, -14.814f, 230.765f}, {-34.949f, -14.342f, 229.809f}, {-35.477f, -14.072f, 228.420f}, {-36.570f, -14.498f, 228.063f}}},
    {"GLU", {{-34.677f, -13.353f, 227.641f}, {-35.030f, -12.991f, 226.273f}, {-34.728f, -14.158f, 225.333f}, {-35.482f, -14.422f, 224.388f}}},
    {"ASN", {{-33.645f, -14.875f, 225.612f}, {-33.249f, -16.007f, 224.792f}, {-34.011f, -17.276f, 225.201f}, {-33.574f, -18.011f, 226.085f}}},
    {"ARG", {{-35.120f, -17.555f, 224.519f}, {-35.940f, -18.715f, 224.825f}, {-36.367f, -19.442f, 223.563f}, {-36.713f, -18.805f, 222.576f}}},
    {"ASP", {{-36.336f, -20.772f, 223.589f}, {-36.761f, -21.574f, 222.451f}, {-37.908f, -22.442f, 222.896f}, {-38.035f, -22.753f, 224.083f}}},
};

// A residue of a test peptide: one of k4r_backbone, its number, insertion code and name (NULL for its own)
typedef struct {
    int res;
    int seq;
    char icode;
    const char* name;
} peptide_pick_t;

// The peptide as the ATOM records of chain A of a PDB file
static md_str_t peptide_pdb(md_allocator_i* alloc, const peptide_pick_t* pick, size_t count) {
    static const char* atom_name[4] = {" N  ", " CA ", " C  ", " O  "};
    static const char* element[4]   = {"N", "C", "C", "O"};
    const size_t cap = count * 4 * 82 + 1;
    char* buf = md_alloc(alloc, cap);
    size_t len = 0;
    int serial = 1;
    for (size_t i = 0; i < count; ++i) {
        const peptide_residue_t* r = &k4r_backbone[pick[i].res];
        for (int a = 0; a < 4; ++a) {
            len += (size_t)snprintf(buf + len, cap - len, "ATOM  %5d %4s %3s A%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                serial++, atom_name[a], pick[i].name ? pick[i].name : r->name, pick[i].seq, pick[i].icode ? pick[i].icode : ' ',
                r->xyz[a][0], r->xyz[a][1], r->xyz[a][2], 1.0, 0.0, element[a]);
        }
    }
    return (md_str_t){buf, len};
}

static bool load_peptide(md_system_t* sys, md_system_state_t* state, const peptide_pick_t* pick, size_t count, md_allocator_i* alloc) {
    if (!md_pdb_system_init_from_str(sys, state, peptide_pdb(alloc, pick, count), MD_PDB_OPTION_NONE, NULL)) return false;
    return md_system_infer(sys, state, MD_SYSTEM_INFER_ALL, NULL);
}

static size_t count_residues_with(const md_system_t* sys, md_residue_flags_t flag) {
    size_t count = 0;
    for (size_t i = 0; i < sys->residue.count; ++i) count += (md_system_residue_flags(sys, i) & flag) != 0;
    return count;
}

// 52, 52A and 52B are three residues of the same number, also when two of them share a name, and a chain whose backbone is
// bonded through them is not broken by their numbers
UTEST(pdb, insertion_codes) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 50, 0, NULL}, {1, 51, 0, NULL}, {2, 52, 0, "CYS"}, {3, 52, 'A', "CYS"}, {4, 52, 'B', NULL}, {5, 53, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, MD_COUNTOF(pick), alloc));

    ASSERT_EQ((size_t)6, sys.residue.count);
    for (size_t i = 0; i < sys.residue.count; ++i) {
        EXPECT_EQ(pick[i].seq, md_system_residue_seq_id(&sys, i));
        EXPECT_EQ((size_t)4, md_system_residue_atom_count(&sys, i));
    }
    EXPECT_EQ((size_t)1, sys.instance.count);
    EXPECT_EQ((size_t)0, count_residues_with(&sys, MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
    ASSERT_EQ((size_t)1, sys.protein_backbone.range.count);
    EXPECT_EQ((size_t)6, sys.protein_backbone.segment.count);

    md_vm_arena_destroy(alloc);
}

// Numbering that skips (1, 2, 3, 14, 15, 16) over a bonded backbone is no break
UTEST(pdb, numbering_jump_is_no_break) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 1, 0, NULL}, {1, 2, 0, NULL}, {2, 3, 0, NULL}, {3, 14, 0, NULL}, {4, 15, 0, NULL}, {5, 16, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, MD_COUNTOF(pick), alloc));

    ASSERT_EQ((size_t)6, sys.residue.count);
    EXPECT_EQ((size_t)0, count_residues_with(&sys, MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
    ASSERT_EQ((size_t)1, sys.protein_backbone.range.count);
    EXPECT_EQ((size_t)6, sys.protein_backbone.segment.count);

    md_atom_idx_t atoms[2];
    for (size_t i = 0; i < sys.residue.count; ++i) {
        EXPECT_FALSE(md_system_chain_break_atoms(atoms, &sys, i));
        EXPECT_EQ(-1, atoms[0]);
        EXPECT_EQ(-1, atoms[1]);
    }

    md_vm_arena_destroy(alloc);
}

// Two residues cut out of the chain: a break after the fourth residue, whether the numbering skips them or runs on. The
// chain is one instance with two backbone ranges, and neither side of the break is a terminus.
static void expect_break_after_fourth(int* utest_result, const md_system_t* sys) {
    ASSERT_EQ((size_t)8, sys->residue.count);
    EXPECT_EQ((size_t)1, sys->instance.count);
    for (size_t i = 0; i < sys->residue.count; ++i) {
        const md_residue_flags_t flags = md_system_residue_flags(sys, i);
        EXPECT_EQ(i == 3, (flags & MD_RESIDUE_FLAG_BREAK_AFTER)    != 0);
        EXPECT_EQ(i == 4, (flags & MD_RESIDUE_FLAG_BREAK_BEFORE) != 0);
    }
    EXPECT_FALSE(md_system_residue_flags(sys, 3) & MD_RESIDUE_FLAG_TERMINAL_END);
    EXPECT_FALSE(md_system_residue_flags(sys, 4) & MD_RESIDUE_FLAG_TERMINAL_BEG);

    ASSERT_EQ((size_t)2, sys->protein_backbone.range.count);
    EXPECT_EQ(0u, sys->protein_backbone.range.offset[0]);
    EXPECT_EQ(4u, sys->protein_backbone.range.offset[1]);
    EXPECT_EQ(8u, sys->protein_backbone.range.offset[2]);

    // The C of the fourth residue and the N of the fifth
    md_atom_idx_t atoms[2];
    ASSERT_TRUE(md_system_chain_break_atoms(atoms, sys, 3));
    EXPECT_EQ(3 * 4 + 2, atoms[0]);
    EXPECT_EQ(4 * 4 + 0, atoms[1]);
    EXPECT_FALSE(md_system_chain_break_atoms(atoms, sys, 4));
}

UTEST(pdb, chain_break_numbering_runs_on) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 1, 0, NULL}, {1, 2, 0, NULL}, {2, 3, 0, NULL}, {3, 4, 0, NULL}, {6, 5, 0, NULL}, {7, 6, 0, NULL}, {8, 7, 0, NULL}, {9, 8, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, MD_COUNTOF(pick), alloc));
    expect_break_after_fourth(utest_result, &sys);

    // The backbones of atoms are found from their bonds alone: without bonds there are none, and no breaks
    sys.bond.count = 0;
    md_array_shrink(sys.bond.pairs, 0);
    md_array_shrink(sys.bond.flags, 0);
    md_bond_build_connectivity(&sys.bond, sys.atom.count, alloc);
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_BACKBONE_BIT, NULL));
    EXPECT_EQ(0u, sys.protein_backbone.range.count);
    EXPECT_EQ(0u, sys.protein_backbone.segment.count);
    for (size_t i = 0; i < sys.residue.count; ++i) {
        EXPECT_EQ(0u, (unsigned)(md_system_residue_flags(&sys, i) & (MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER)));
    }

    md_vm_arena_destroy(alloc);
}

UTEST(pdb, chain_break_numbering_skips) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 1, 0, NULL}, {1, 2, 0, NULL}, {2, 3, 0, NULL}, {3, 4, 0, NULL}, {6, 7, 0, NULL}, {7, 8, 0, NULL}, {8, 9, 0, NULL}, {9, 10, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, MD_COUNTOF(pick), alloc));
    expect_break_after_fourth(utest_result, &sys);

    // Inferring again finds the same breaks, not more
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    expect_break_after_fourth(utest_result, &sys);

    md_vm_arena_destroy(alloc);
}

// 1LAF ends chain E at Asp 238, whose carboxyl group is whole (OXT): the arginine numbered 240 after it is a free
// amino acid bound to the protein, not the chain going on. 1k4r and 1a64 are whole (1a64 misses the first residues of its
// chains only, which leaves no break)
UTEST(pdb, chain_breaks_of_files) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb"), MD_PDB_OPTION_NONE, NULL));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ((size_t)0, count_residues_with(&sys, MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
        size_t asp = SIZE_MAX;
        for (size_t i = 0; i + 1 < sys.residue.count; ++i) {
            if (md_system_residue_seq_id(&sys, i) == 238 && md_system_residue_seq_id(&sys, i + 1) == 240) asp = i;
        }
        ASSERT_NE(SIZE_MAX, asp);
        EXPECT_TRUE(md_str_eq(md_system_residue_name(&sys, asp), MD_STR_LIT("ASP")));
        EXPECT_TRUE(md_system_residue_flags(&sys, asp) & MD_RESIDUE_FLAG_TERMINAL_END);
    }
    const md_str_t whole[] = {MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1a64.pdb")};
    for (size_t f = 0; f < MD_COUNTOF(whole); ++f) {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, whole[f], MD_PDB_OPTION_NONE, NULL));
        ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ((size_t)0, count_residues_with(&sys, MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
    }
    md_vm_arena_destroy(alloc);
}

// The writer ends each part with a TER before the next part begins, never within a residue: two waters of a GRO file,
// which has no chains, read back as two residues of three atoms
UTEST(pdb, write_ter_between_parts) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const char* gro =
        "two waters\n6\n"
        "    1SOL     OW    1   0.000   0.000   0.000\n"
        "    1SOL    HW1    2   0.096   0.000   0.000\n"
        "    1SOL    HW2    3  -0.024   0.093   0.000\n"
        "    2SOL     OW    4   0.500   0.500   0.500\n"
        "    2SOL    HW1    5   0.596   0.500   0.500\n"
        "    2SOL    HW2    6   0.476   0.593   0.500\n"
        "   2.00000   2.00000   2.00000\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, md_str_from_cstr(gro), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(2u, md_system_part_count(&sys));

    const md_str_t path = MD_STR_LIT("md_unittest_write_ter.pdb");
    md_file_t f = {0};
    ASSERT_TRUE(md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    ASSERT_TRUE(md_pdb_system_write_state_to_file(f, &sys, &st, NULL, sys.atom.count, 0));
    md_file_close(&f);

    md_system_t back = {.alloc = alloc};
    md_system_state_t back_st = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&back, &back_st, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL));
    remove(path.ptr);
    ASSERT_EQ(2u, back.residue.count);
    EXPECT_EQ(3u, md_system_residue_atom_range(&back, 0).end - md_system_residue_atom_range(&back, 0).beg);
    EXPECT_EQ(3u, md_system_residue_atom_range(&back, 1).end - md_system_residue_atom_range(&back, 1).beg);
    md_vm_arena_destroy(alloc);
}

// ---------------------------------------------------------------------------
// Straight from the text, in blocks and pieces
//
// A file of more than one block (md_reader.h): chains ended by TER, residues straddling the blocks and the pieces,
// insertion codes, alternate locations, atoms without an element, CONECT records at the end and a second model which is
// not loaded. The same system read on the calling thread, on threads, and with the pieces out of order; from the file
// and from a string.
// ---------------------------------------------------------------------------

#include "test_executor.h"
#include "test_load_blocks.h"

static void pdb_blocks_atom(md_array(char)* text, md_allocator_i* alloc, int serial, const char* name, const char* res, char chain, int seq, char icode, char alt, float x, float y, float z, const char* element) {
    char line[128];
    const int n = snprintf(line, sizeof(line), "ATOM  %5d %-4s%c%-3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                           serial % 100000, name, alt, res, chain, seq % 10000, icode, x, y, z, 1.0, 0.5, element);
    md_array_push_array(*text, line, (size_t)n, alloc);
}

UTEST(pdb, blocks_and_pieces) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(4));
    md_array(char) text = 0;
    const char* header = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\nMODEL        1\n";
    md_array_push_array(text, header, strlen(header), alloc);

    static const char* names[4] = { "N", "CA", "C", "O" };
    static const char* elements[4] = { "N", "C", "C", "O" };
    int serial = 1;
    const int num_chains = 4, num_res = 30000;  // 480000 atom lines and some: two blocks
    for (int c = 0; c < num_chains; ++c) {
        for (int r = 0; r < num_res; ++r) {
            // Every 50th residue an insertion after it, of the same number; every 70th residue without elements
            const int copies = r % 50 == 0 ? 2 : 1;
            for (int k = 0; k < copies; ++k) {
                for (int a = 0; a < 4; ++a) {
                    const float x = (float)(serial % 701) * 0.1f, y = (float)(serial % 709) * 0.1f, z = (float)(serial % 719) * 0.1f;
                    const char* el = r % 70 == 0 ? "" : elements[a];
                    pdb_blocks_atom(&text, alloc, serial++, names[a], "ALA", (char)('A' + c), r + 1, k ? 'A' : ' ', ' ', x, y, z, el);
                    // An alternate location B of every 30th residue's CA, which is not loaded
                    if (a == 1 && r % 30 == 0) pdb_blocks_atom(&text, alloc, serial++, names[a], "ALA", (char)('A' + c), r + 1, k ? 'A' : ' ', 'B', x + 1, y, z, el);
                }
            }
        }
        const char* ter = "TER\n";  // As GROMACS writes it: the three letters alone
        md_array_push_array(text, ter, strlen(ter), alloc);
        char line[64];
        const int n = snprintf(line, sizeof(line), "TER   %5d      ALA %c%4d\n", serial++, 'A' + c, num_res % 10000);
        md_array_push_array(text, line, (size_t)n, alloc);
    }
    // Two residues of the same name and number across a TER: two residues, and two units of the same chain id
    pdb_blocks_atom(&text, alloc, serial++, "ZN", "ZN", 'Z', 1, ' ', ' ', 1, 2, 3, "ZN");
    md_array_push_array(text, "TER\n", 4, alloc);
    md_array_push_array(text, "TER   \n", 7, alloc);
    pdb_blocks_atom(&text, alloc, serial++, "ZN", "ZN", 'Z', 1, ' ', ' ', 4, 5, 6, "ZN");
    md_array_push_array(text, "ENDMDL\nMODEL        2\n", 22, alloc);
    pdb_blocks_atom(&text, alloc, serial++, "XX", "XXX", 'X', 1, ' ', ' ', 0, 0, 0, "C");
    md_array_push_array(text, "ENDMDL\n", 7, alloc);
    // A disulfide between two atoms far apart, in the CONECT records at the end
    const char* conect = "CONECT    2  400\nCONECT  400    2\nEND\n";
    md_array_push_array(text, conect, strlen(conect), alloc);
    const md_str_t str = { text, md_array_size(text) };
    ASSERT_GT(str.len, (size_t)MD_MEGABYTES(33));

    const char* path = "md_unittest_blocks.pdb";
    ASSERT_TRUE(load_blocks_write(path, str));

    md_system_t ref = { .alloc = alloc };
    md_system_state_t ref_st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&ref, &ref_st, str, MD_PDB_OPTION_NONE, NULL));
    EXPECT_EQ((size_t)(num_chains * (num_res + num_res / 50) * 4 + 2), ref.atom.count);
    EXPECT_EQ((size_t)(num_chains * (num_res + num_res / 50) + 2), ref.residue.count);
    EXPECT_EQ((size_t)1, ref.bond.count);
    EXPECT_EQ((size_t)(num_chains + 2), ref.instance.count);

    test_executor_t threads  = test_executor_threads(4);
    test_executor_t shuffled = test_executor_shuffled(3);
    const md_executor_t* executors[] = { NULL, test_executor(&threads), test_executor(&shuffled) };
    for (size_t e = 0; e < MD_COUNTOF(executors); ++e) {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, md_str_from_cstr(path), MD_PDB_OPTION_NONE, executors[e]));
        EXPECT_EQ((size_t)0, load_blocks_diff(&ref, &ref_st, &sys, &st));
        md_system_t sys2 = { .alloc = alloc };
        md_system_state_t st2 = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys2, &st2, str, MD_PDB_OPTION_NONE, executors[e]));
        EXPECT_EQ((size_t)0, load_blocks_diff(&ref, &ref_st, &sys2, &st2));
    }

    // All models: the atoms of the second as well
    {
        md_system_t all = { .alloc = alloc }, sys = { .alloc = alloc };
        md_system_state_t all_st = { .alloc = alloc }, st = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_str(&all, &all_st, str, MD_PDB_OPTION_CONCAT_MODELS, NULL));
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, md_str_from_cstr(path), MD_PDB_OPTION_CONCAT_MODELS, test_executor(&threads)));
        EXPECT_EQ(ref.atom.count + 1, all.atom.count);
        EXPECT_EQ((size_t)0, load_blocks_diff(&all, &all_st, &sys, &st));
    }

    remove(path);
    md_vm_arena_destroy(alloc);
}

// A TER record ends a chain, also the three letters alone and also between two chains of the same id: the segments of
// an EM structure (T4SS, chain J of segment Z and of segment ZA) are two units, with no break of a chain between them
UTEST(pdb, ter_ends_the_unit) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "ATOM      1  N   GLY J  48      10.000  10.000  10.000  1.00  0.00      Z    N\n"
        "ATOM      2  CA  GLY J  48      11.458  10.000  10.000  1.00  0.00      Z    C\n"
        "ATOM      3  C   GLY J  48      12.009  11.420  10.000  1.00  0.00      Z    C\n"
        "ATOM      4  O   GLY J  48      11.251  12.390  10.000  1.00  0.00      Z    O\n"
        "ATOM      5  N   GLY J  49      13.332  11.540  10.000  1.00  0.00      Z    N\n"
        "ATOM      6  CA  GLY J  49      13.970  12.850  10.000  1.00  0.00      Z    C\n"
        "ATOM      7  C   GLY J  49      15.480  12.730  10.000  1.00  0.00      Z    C\n"
        "ATOM      8  O   GLY J  49      16.050  11.640  10.000  1.00  0.00      Z    O\n"
        "TER\n"
        "ATOM      9  N   GLY J 208      60.000  60.000  60.000  1.00  0.00      ZA   N\n"
        "ATOM     10  CA  GLY J 208      61.458  60.000  60.000  1.00  0.00      ZA   C\n"
        "ATOM     11  C   GLY J 208      62.009  61.420  60.000  1.00  0.00      ZA   C\n"
        "ATOM     12  O   GLY J 208      61.251  62.390  60.000  1.00  0.00      ZA   O\n"
        "ATOM     13  N   GLY J 209      63.332  61.540  60.000  1.00  0.00      ZA   N\n"
        "ATOM     14  CA  GLY J 209      63.970  62.850  60.000  1.00  0.00      ZA   C\n"
        "ATOM     15  C   GLY J 209      65.480  62.730  60.000  1.00  0.00      ZA   C\n"
        "ATOM     16  O   GLY J 209      66.050  61.640  60.000  1.00  0.00      ZA   O\n"
        "END\n");
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t st = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
        ASSERT_EQ((size_t)4, sys.residue.count);
        ASSERT_EQ((size_t)2, sys.instance.count);
        EXPECT_EQ(2u, sys.instance.res_offset[1]);
        EXPECT_TRUE(md_str_eq(md_system_instance_auth_id(&sys, 0), MD_STR_LIT("J")));
        EXPECT_TRUE(md_str_eq(md_system_instance_auth_id(&sys, 1), MD_STR_LIT("J")));
        ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ((size_t)0, count_residues_with(&sys, MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
        md_atom_idx_t atoms[2];
        EXPECT_FALSE(md_system_chain_break_atoms(atoms, &sys, 1));
    }
    md_vm_arena_destroy(alloc);
}

// The charge column is the formal charge the file gives (atom/formal_charge), which chemistry perception starts from:
// the glutamate of T4SS whose OE2 is "O1-" keeps its charge there, not on the OE1 perception would choose
UTEST(pdb, formal_charge_column) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "ATOM  99999  N   GLU J  48     278.601 255.559 389.405  1.00 44.59      Z    N  \n"
        "ATOM  99999  CA  GLU J  48     278.674 256.483 390.532  1.00 44.59      Z    C  \n"
        "ATOM  99999  C   GLU J  48     278.403 257.917 390.089  1.00 44.59      Z    C  \n"
        "ATOM  99999  O   GLU J  48     279.067 258.429 389.186  1.00 44.59      Z    O  \n"
        "ATOM  99999  CB  GLU J  48     280.040 256.392 391.216  1.00 44.59      Z    C  \n"
        "ATOM  99999  CG  GLU J  48     279.973 256.457 392.735  1.00 44.59      Z    C  \n"
        "ATOM  99999  CD  GLU J  48     281.267 256.030 393.411  1.00 44.59      Z    C  \n"
        "ATOM  99999  OE1 GLU J  48     281.325 256.024 394.663  1.00 44.59      Z    O  \n"
        "ATOM  99999  OE2 GLU J  48     282.229 255.699 392.685  1.00 44.59      Z    O1-\n");
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t st = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
        ASSERT_TRUE(md_attributes_find(&sys.attributes, MD_STR_LIT("atom/formal_charge")) != NULL);
        ASSERT_TRUE(md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL));
        EXPECT_EQ(-1, (int)md_atom_formal_charge(&sys.atom, 8));
        EXPECT_EQ(0,  (int)md_atom_formal_charge(&sys.atom, 7));
    }
    md_vm_arena_destroy(alloc);
}

// The chain id is whatever character column 22 holds, a digit as well as a letter, kept as a string (a blank is no
// chain): '1' and 'A' are two chains, and '1' is not the number 1
UTEST(pdb, chain_ids_as_given) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "ATOM      1  CA  GLY 1   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY 1   2       3.800   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CA  GLY A   1      20.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  CA  GLY a   1      40.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      5  CA  GLY #   1      60.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      6  CA  GLY     1      80.000   0.000   0.000  1.00  0.00           C\n");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)6, sys.atom.count);
    ASSERT_EQ((size_t)5, sys.instance.count);
    const char* ids[] = { "1", "A", "a", "#", "" };
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_TRUE(md_str_eq(md_system_instance_auth_id(&sys, i), md_str_from_cstr(ids[i])));
    }
    md_vm_arena_destroy(alloc);
}

// A residue number which is not a number of the format ('****' past 9999) is the one after the residue before
UTEST(pdb, residue_numbers_which_are_not_numbers) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "ATOM      1  OW  SOL X9998       0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  OW  SOL X9999       3.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      3  OW  SOL X****       6.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      4  OW  SOL X####       9.000   0.000   0.000  1.00  0.00           O\n");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL));
    ASSERT_EQ((size_t)4, sys.residue.count);
    EXPECT_EQ(9998, sys.residue.seq_id[0]);
    EXPECT_EQ(9999, sys.residue.seq_id[1]);
    EXPECT_EQ(10000, sys.residue.seq_id[2]);
    EXPECT_EQ(10001, sys.residue.seq_id[3]);
    md_vm_arena_destroy(alloc);
}

// A file of many models is read to the end of the first, and then from the last MODEL or ENDMDL record on, for the
// CONECT records which follow the models: of a file larger than what is looked through at the end at first, and of one
// whose models have no ENDMDL. And the same where the rest is a single model's CONECT records, which are read through.
UTEST(pdb, models_then_conect) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    test_executor_t threads = test_executor_threads(3);
    const md_executor_t* executors[] = { NULL, test_executor(&threads) };
    for (int variant = 0; variant < 3; ++variant) {
        md_array(char) text = 0;
        // More than a pass (md_text_pass_size): the end of the first model is not in the pass which has the CONECT records
        const int num_models = variant == 2 ? 1 : 7000;
        const int num_atoms = variant == 2 ? 70000 : 10;
        char line[128];
        for (int m = 0; m < num_models; ++m) {
            int n = snprintf(line, sizeof(line), "MODEL     %4d\n", m + 1);
            md_array_push_array(text, line, (size_t)n, alloc);
            for (int i = 0; i < num_atoms; ++i) {
                n = snprintf(line, sizeof(line), "ATOM  %5d  C%-2d LIG A   1    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
                    i + 1, i % 100, (float)i * 1.5f, (float)m * 0.001f, 0.0f);
                md_array_push_array(text, line, (size_t)n, alloc);
            }
            if (variant != 1) md_array_push_array(text, "ENDMDL\n", 7, alloc);
        }
        // The bonds of the first model's atoms, numbered as they are in every model
        for (int i = 1; i < num_atoms; ++i) {
            const int n = snprintf(line, sizeof(line), "CONECT%5d%5d\n", i, i + 1);
            md_array_push_array(text, line, (size_t)n, alloc);
        }
        md_array_push_array(text, "END\n", 4, alloc);
        const md_str_t str = { text, md_array_size(text) };
        ASSERT_GT(str.len, (size_t)MD_MEGABYTES(5));

        const char* path = "md_unittest_models_conect.pdb";
        ASSERT_TRUE(load_blocks_write(path, str));
        for (size_t e = 0; e < MD_COUNTOF(executors); ++e) {
            md_system_t sys = { .alloc = alloc }, sys2 = { .alloc = alloc };
            md_system_state_t st = { .alloc = alloc }, st2 = { .alloc = alloc };
            ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, md_str_from_cstr(path), MD_PDB_OPTION_NONE, executors[e]));
            ASSERT_TRUE(md_pdb_system_init_from_str(&sys2, &st2, str, MD_PDB_OPTION_NONE, executors[e]));
            EXPECT_EQ((size_t)num_atoms, sys.atom.count);
            EXPECT_EQ((size_t)(num_atoms - 1), sys.bond.count);
            EXPECT_EQ((size_t)0, load_blocks_diff(&sys, &st, &sys2, &st2));
            // Of the first model
            EXPECT_EQ(0.0f, st.xyz[num_atoms - 1].y);
        }
        remove(path);
    }
    md_vm_arena_destroy(alloc);
}

// CRYST1 records cut short, or of no volume. Columns past the end of the line were read from beyond it; now a short
// record gives what it has (lengths without angles are a box of right angles), and a record that describes no cell
// - the PDB's 1 1 1 placeholder, a zero length, an angle of 0 - gives none rather than a cell of infinities.
UTEST(pdb, cryst1_short_or_empty) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const char* atom = "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n";
    static const struct { const char* cryst1; bool cell; float x; } cases[] = {
        { "CRYST1   40.000   50.000   60.000\n",                                         true,  40.0f },
        { "CRYST1   40.000   50.000   60.000  90.00  90.00  90.00 P 1           1\n",   true,  40.0f },
        { "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n",   false, 0.0f },
        { "CRYST1    0.000    0.000    0.000   0.00   0.00   0.00 P 1           1\n",   false, 0.0f },
        { "CRYST1   40.000   50.000   60.000   0.00  90.00  90.00 P 1           1\n",   false, 0.0f },
        { "CRYST1   40.000\n",                                                          false, 0.0f },
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        char text[512];
        snprintf(text, sizeof(text), "%s%s", cases[i].cryst1, atom);
        md_system_t sys = {.alloc = alloc};
        md_system_state_t st = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, md_str_from_cstr(text), MD_PDB_OPTION_NONE, NULL));
        const bool has_cell = st.unitcell.flags != 0;
        EXPECT_EQ(cases[i].cell, has_cell);
        if (cases[i].cell) {
            EXPECT_NEAR(cases[i].x, (float)st.unitcell.x, 1e-4f);
            EXPECT_TRUE(isfinite(st.unitcell.y) && isfinite(st.unitcell.z));
        }
    }
    md_vm_arena_destroy(alloc);
}
