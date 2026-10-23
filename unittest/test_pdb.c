#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include "run_check.h"

UTEST(pdb, parse_ordinary) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/1k4r.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 9084);
    EXPECT_EQ(pdb_data.num_connections, 36);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_helices, 24);
    EXPECT_EQ(pdb_data.num_sheets, 102);

    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, tryptophan) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/tryptophan.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 28);

    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    EXPECT_TRUE(md_pdb_system_init_from_data(&sys, &sys_state, &pdb_data, MD_PDB_OPTION_NONE));

    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, unmatched_model_entry) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/dppc64.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 14738);
    EXPECT_EQ(pdb_data.num_connections, 0);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_helices, 0);
    EXPECT_EQ(pdb_data.num_sheets, 0);

    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, parse_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 38);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 5814);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_connections, 0);
    EXPECT_EQ(pdb_data.num_helices, 0);
    EXPECT_EQ(pdb_data.num_sheets, 0);

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (int64_t i = 0; i < pdb_data.num_models; ++i) {
        char data[6] = {0};
        md_file_seek(file, pdb_data.models[i].byte_offset, MD_FILE_BEG);
        md_file_read(file, data, 5);
        EXPECT_EQ(strncmp(data, "MODEL", 5), 0);
    }
    md_file_close(&file);

    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

// REMARK 350 of 1k4r: biomolecule 1 is its three chains by 60 transforms, the shell of the virion
UTEST(pdb, assemblies) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE));
    const md_assembly_data_t* a = &sys.assembly;
    ASSERT_EQ(1, a->count);
    EXPECT_STREQ("1", str_ptr(md_assembly_id(a, 0)));
    ASSERT_EQ(1, a->part_count);
    EXPECT_EQ(60, a->oper_count);
    const md_urange_t ir = md_assembly_part_inst_range(a, 0);
    EXPECT_EQ(3, ir.end - ir.beg);
    EXPECT_EQ(60, md_assembly_copy_count(a, 0));
    const mat4_t I = mat4_ident();
    EXPECT_EQ(0, MEMCMP(&I, &a->oper_transform[a->oper[0]], sizeof(mat4_t)));
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
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, (str_t){pdb, sizeof(pdb) - 1}, MD_PDB_OPTION_NONE));
    ASSERT_EQ(3, sys.instance.count);
    const md_assembly_data_t* a = &sys.assembly;
    ASSERT_EQ(1, a->count);
    // The identity copies A, B and C; the other only A
    ASSERT_EQ(2, a->oper_count);
    ASSERT_EQ(2, a->part_count);
    md_urange_t ir = md_assembly_part_inst_range(a, 0);
    md_urange_t orr = md_assembly_part_oper_range(a, 0);
    EXPECT_EQ(3, ir.end - ir.beg);
    EXPECT_EQ(1, orr.end - orr.beg);
    EXPECT_EQ(0, a->oper[orr.beg]);
    ir = md_assembly_part_inst_range(a, 1);
    orr = md_assembly_part_oper_range(a, 1);
    EXPECT_EQ(1, ir.end - ir.beg);
    EXPECT_EQ(1, orr.end - orr.beg);
    EXPECT_EQ(1, a->oper[orr.beg]);
    // A quarter turn about z, by rows as the records give it: x' = -y + 10, y' = x
    const vec4_t q = mat4_mul_vec4(a->oper_transform[1], (vec4_t){1, 2, 3, 1});
    EXPECT_NEAR(8.0f, q.x, 1e-6f);
    EXPECT_NEAR(1.0f, q.y, 1e-6f);
    EXPECT_NEAR(3.0f, q.z, 1e-6f);
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, create_system) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, path, arena));
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = { .alloc = arena };
    EXPECT_TRUE(md_pdb_system_init_from_data(&sys, &sys_state, &pdb_data, MD_PDB_OPTION_NONE));
    ASSERT_EQ(sys.atom.count, pdb_data.num_atom_coordinates);

    EXPECT_EQ(1185, sys.component.count);
    EXPECT_EQ(3, sys.instance.count);
    EXPECT_EQ(1, sys.entity.count);

    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, pdb_data.atom_coordinates[i].x);
        EXPECT_EQ(sys_state.xyz[i].y, pdb_data.atom_coordinates[i].y);
        EXPECT_EQ(sys_state.xyz[i].z, pdb_data.atom_coordinates[i].z);
    }

    md_system_free(&sys);
    
    md_pdb_data_free(&pdb_data, arena);
    md_arena_allocator_destroy(arena);
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
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
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(pdb), MD_PDB_OPTION_NONE));
    EXPECT_NEAR(246.108, state.unitcell.x, 1e-3);
    EXPECT_NEAR(246.108, state.unitcell.y, 1e-3);
    EXPECT_NEAR(398.534, state.unitcell.z, 1e-3);
    EXPECT_NEAR(0.0, state.unitcell.xy, 1e-9);
    EXPECT_NEAR(0.0, state.unitcell.xz, 1e-9);
    EXPECT_NEAR(0.0, state.unitcell.yz, 1e-9);
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    EXPECT_GE(md_bond_find(&sys.bond, 0, 1), 0);    // N-CA, across the boundary in x
    EXPECT_EQ((int)MD_COMPONENT_KIND_AMINO_ACID, (int)md_component_kind(&sys.component, 0));
    md_arena_allocator_destroy(alloc);
}

// Past 99999 atoms VMD writes the serials in hexadecimal (186a0 for 100000), and '*****' past FFFFF: not hybrid-36,
// which begins a number past 99999 with a letter. They are not numbers of the format: 0, counted once rather than
// reported line by line (millions of lines in a large system), and the atoms are read in their order.
UTEST(pdb, serials_which_are_not_numbers) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    static const char pdb[] =
        "ATOM  99999  OW  SOL X3333       0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  186a0  OW  SOL X3334       3.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  18700  OW  SOL X3335       6.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  *****  OW  SOL X3336       9.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM  A0000  OW  SOL X3337      12.000   0.000   0.000  1.00  0.00           O\n"
        "END\n";
    md_pdb_data_t data = {0};
    ASSERT_TRUE(md_pdb_data_parse_str(&data, str_from_cstr(pdb), alloc));
    ASSERT_EQ(5, (int)data.num_atom_coordinates);
    EXPECT_EQ(99999,  data.atom_coordinates[0].atom_serial);
    EXPECT_EQ(0,      data.atom_coordinates[1].atom_serial);    // hexadecimal: not a number of the format
    EXPECT_EQ(18700,  data.atom_coordinates[2].atom_serial);    // hexadecimal without a letter is decimal
    EXPECT_EQ(0,      data.atom_coordinates[3].atom_serial);
    EXPECT_EQ(100000, data.atom_coordinates[4].atom_serial);    // hybrid-36
    EXPECT_NEAR(9.0f, data.atom_coordinates[3].x, 1e-6f);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, conect_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE));
    EXPECT_EQ((size_t)18, count_bonds_of_origin(&sys, MD_BOND_ORIGIN_FILE));
    for (size_t i = 0; i < sys.bond.count; ++i) {
        if (md_bond_origin(sys.bond.flags[i]) != MD_BOND_ORIGIN_FILE) continue;
        EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, sys.bond.pairs[i].idx[0]), STR_LIT("SG")));
        EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, sys.bond.pairs[i].idx[1]), STR_LIT("SG")));
    }
    // Atoms 23 and 226 (Cys 3 and Cys 30): one bond, from the file
    const md_bond_idx_t b = md_bond_find(&sys.bond, 22, 225);
    ASSERT_NE(-1, b);
    EXPECT_EQ(MD_BOND_ORIGIN_FILE, md_bond_origin(sys.bond.flags[b]));
    EXPECT_GT(count_bonds_of_origin(&sys, MD_BOND_ORIGIN_INFERRED), (size_t)9000);

    const size_t num_bonds = sys.bond.count;
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
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
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_pdb_data_t data = {0};
    const str_t text = STR_LIT(
        "HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00  0.00          FE\n"
        "CONECT    1    2    3    4    5\n"
        "CONECT    1    6    7    8\n"
        "CONECT    1    9\n"
        "CONECT99999A0000\n");
    ASSERT_TRUE(md_pdb_data_parse_str(&data, text, alloc));
    ASSERT_EQ((size_t)2, data.num_connections);
    for (int k = 1; k < 8; ++k) {
        EXPECT_EQ(k + 1, data.connections[0].atom_serial[k]);
    }
    EXPECT_EQ(99999, data.connections[1].atom_serial[0]);
    EXPECT_EQ(100000, data.connections[1].atom_serial[1]);
    md_vm_arena_destroy(alloc);
}

// An atom name of four characters fills columns 13 to 16 and is kept whole: HN4' and HN4A are two atoms, not two HN4
UTEST(pdb, four_character_atom_names) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    const str_t text = STR_LIT(
        "HETATM    1  N4' TD6 B   1      43.654  64.837  38.628  1.00  0.00           N\n"
        "HETATM    2 HN4' TD6 B   1      44.050  64.155  39.051  1.00  0.00           H\n"
        "HETATM    3 HN4A TD6 B   1      44.059  65.650  38.601  1.00  0.00           H\n"
        "HETATM    4  CA  TD6 B   1      45.000  65.000  38.000  1.00  0.00           C\n");
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE));
    ASSERT_EQ((size_t)4, sys.atom.count);
    EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, 0), STR_LIT("N4'")));
    EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, 1), STR_LIT("HN4'")));
    EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, 2), STR_LIT("HN4A")));
    EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, 3), STR_LIT("CA")));
    md_vm_arena_destroy(alloc);
}

UTEST(pdb, parse_nonexistent_file) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/nonexistent.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_FALSE(result);
    
    // Should be safe to free even when parsing failed
    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, parse_empty_path) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = {0};  // Empty string
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, md_get_heap_allocator());
    EXPECT_FALSE(result);
    
    md_pdb_data_free(&pdb_data, md_get_heap_allocator());
    md_arena_allocator_destroy(alloc);
}
// ### RUN ###

#define PDB_RUN STR_LIT("run/1ala")

// Recorded from the trajectory reader this replaced
static const run_ref_t pdb_refs[] = {
    { 0,  {3568.34002, 7394.96399, 3699.66098}, {23.7040005, 23, 21.3549995},         {25.0020008, 75.5550003, 23.5319996}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
    { 19, {3297.43299, 6833.70401, 3418.89901}, {22.4950008, 43.8180008, 11.6870003}, {7.09700012, 53.9360008, 23.4740009}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
    { 37, {3298.315, 6835.46697, 3419.76601},   {28.4710007, 57.6910019, 25.2150002}, {8.61200047, 35.6669998, 14.9329996}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
};

// One frame per model, time as model ordinals, and the file's one CRYST1 cell at every frame.
UTEST(pdb, run_matches_reference) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = {.alloc = arena};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &sys_state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, path, PDB_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, PDB_RUN, 38, 153, pdb_refs, ARRAY_SIZE(pdb_refs));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/1ala/time"));
    const md_attribute_t* pos  = md_attributes_find(&sys.attributes, STR_LIT("run/1ala/atom/position"));
    ASSERT_TRUE(time && pos);
    EXPECT_TRUE(md_unit_is_none(time->unit));
    EXPECT_EQ(9.0, ((const double*)time->data)[9]);

    // The first model is the structure's own coordinates.
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, sys.atom.count);
    ASSERT_TRUE(run_extract_one(&got, &sys, PDB_RUN, 0));
    EXPECT_EQ(0, MEMCMP(sys_state.xyz, got.xyz, sys.atom.count * sizeof(vec3_t)));

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

// A file of one model is a structure, not a run.
UTEST(pdb, run_needs_several_models) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    EXPECT_FALSE(md_pdb_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), PDB_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/1ala/time")) == NULL);
    md_vm_arena_destroy(arena);
}

// The columns read without a scan give what parse_float gives on them: every coordinate, occupancy
// and temperature factor of the reference files
#include <core/md_parse.h>

static size_t pdb_test_compare_with_columns(str_t path) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_pdb_data_t pdb = {0};
    if (!md_pdb_data_parse_file(&pdb, path, alloc)) return SIZE_MAX;
    str_t text = load_textfile(path, alloc);
    str_t rest = text, line;
    size_t i = 0, mismatches = 0;
    static const size_t cols[5][2] = { {30, 8}, {38, 8}, {46, 8}, {54, 6}, {60, 6} };
    while (str_extract_line(&line, &rest) && i < md_array_size(pdb.atom_coordinates)) {
        if (!(str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6))) continue;
        const md_pdb_coordinate_t* c = &pdb.atom_coordinates[i++];
        const float got[5] = { c->x, c->y, c->z, c->occupancy, c->temp_factor };
        for (int k = 0; k < 5; ++k) {
            if (line.len < cols[k][0] + cols[k][1]) continue;
            const float want = (float)parse_float(str_substr(line, cols[k][0], cols[k][1]));
            if (memcmp(&want, &got[k], sizeof(float)) != 0) mismatches += 1;
        }
    }
    if (i != md_array_size(pdb.atom_coordinates)) mismatches += 1;
    str_free(text, alloc);
    md_pdb_data_free(&pdb, alloc);
    return mismatches;
}

UTEST(pdb, columns_match_parse_float) {
    EXPECT_EQ(0, pdb_test_compare_with_columns(STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
    EXPECT_EQ(0, pdb_test_compare_with_columns(STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb")));
    EXPECT_EQ(0, pdb_test_compare_with_columns(STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb")));
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
static str_t peptide_pdb(md_allocator_i* alloc, const peptide_pick_t* pick, size_t count) {
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
    return (str_t){buf, len};
}

static bool load_peptide(md_system_t* sys, md_system_state_t* state, const peptide_pick_t* pick, size_t count, md_allocator_i* alloc) {
    md_pdb_data_t data = {0};
    if (!md_pdb_data_parse_str(&data, peptide_pdb(alloc, pick, count), alloc)) return false;
    if (!md_pdb_system_init_from_data(sys, state, &data, MD_PDB_OPTION_NONE)) return false;
    return md_util_system_infer(sys, state, MD_UTIL_INFER_ALL);
}

static size_t count_components_with(const md_system_t* sys, md_component_flags_t flag) {
    size_t count = 0;
    for (size_t i = 0; i < sys->component.count; ++i) count += (md_component_flags(&sys->component, i) & flag) != 0;
    return count;
}

// 52, 52A and 52B are three residues of the same number, also when two of them share a name, and a chain whose backbone is
// bonded through them is not broken by their numbers
UTEST(pdb, insertion_codes) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 50, 0, NULL}, {1, 51, 0, NULL}, {2, 52, 0, "CYS"}, {3, 52, 'A', "CYS"}, {4, 52, 'B', NULL}, {5, 53, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, ARRAY_SIZE(pick), alloc));

    ASSERT_EQ((size_t)6, sys.component.count);
    for (size_t i = 0; i < sys.component.count; ++i) {
        EXPECT_EQ(pick[i].seq, md_component_seq_id(&sys.component, i));
        EXPECT_EQ((size_t)4, md_component_atom_count(&sys.component, i));
    }
    EXPECT_EQ((size_t)1, sys.instance.count);
    EXPECT_EQ((size_t)0, count_components_with(&sys, MD_COMPONENT_FLAG_BREAK_BEFORE | MD_COMPONENT_FLAG_BREAK_AFTER));
    ASSERT_EQ((size_t)1, sys.protein_backbone.range.count);
    EXPECT_EQ((size_t)6, sys.protein_backbone.segment.count);

    md_vm_arena_destroy(alloc);
}

// Numbering that skips (1, 2, 3, 14, 15, 16) over a bonded backbone is no break
UTEST(pdb, numbering_jump_is_no_break) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 1, 0, NULL}, {1, 2, 0, NULL}, {2, 3, 0, NULL}, {3, 14, 0, NULL}, {4, 15, 0, NULL}, {5, 16, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, ARRAY_SIZE(pick), alloc));

    ASSERT_EQ((size_t)6, sys.component.count);
    EXPECT_EQ((size_t)0, count_components_with(&sys, MD_COMPONENT_FLAG_BREAK_BEFORE | MD_COMPONENT_FLAG_BREAK_AFTER));
    ASSERT_EQ((size_t)1, sys.protein_backbone.range.count);
    EXPECT_EQ((size_t)6, sys.protein_backbone.segment.count);

    md_atom_idx_t atoms[2];
    for (size_t i = 0; i < sys.component.count; ++i) {
        EXPECT_FALSE(md_util_chain_break_atoms(atoms, &sys, i));
        EXPECT_EQ(-1, atoms[0]);
        EXPECT_EQ(-1, atoms[1]);
    }

    md_vm_arena_destroy(alloc);
}

// Two residues cut out of the chain: a break after the fourth residue, whether the numbering skips them or runs on. The
// chain is one instance with two backbone ranges, and neither side of the break is a terminus.
static void expect_break_after_fourth(int* utest_result, const md_system_t* sys) {
    ASSERT_EQ((size_t)8, sys->component.count);
    EXPECT_EQ((size_t)1, sys->instance.count);
    for (size_t i = 0; i < sys->component.count; ++i) {
        const md_component_flags_t flags = md_component_flags(&sys->component, i);
        EXPECT_EQ(i == 3, (flags & MD_COMPONENT_FLAG_BREAK_AFTER)  != 0);
        EXPECT_EQ(i == 4, (flags & MD_COMPONENT_FLAG_BREAK_BEFORE) != 0);
    }
    EXPECT_FALSE(md_component_flags(&sys->component, 3) & MD_COMPONENT_FLAG_TERMINAL_END);
    EXPECT_FALSE(md_component_flags(&sys->component, 4) & MD_COMPONENT_FLAG_TERMINAL_BEG);

    ASSERT_EQ((size_t)2, sys->protein_backbone.range.count);
    EXPECT_EQ(0u, sys->protein_backbone.range.offset[0]);
    EXPECT_EQ(4u, sys->protein_backbone.range.offset[1]);
    EXPECT_EQ(8u, sys->protein_backbone.range.offset[2]);

    // The C of the fourth residue and the N of the fifth
    md_atom_idx_t atoms[2];
    ASSERT_TRUE(md_util_chain_break_atoms(atoms, sys, 3));
    EXPECT_EQ(3 * 4 + 2, atoms[0]);
    EXPECT_EQ(4 * 4 + 0, atoms[1]);
    EXPECT_FALSE(md_util_chain_break_atoms(atoms, sys, 4));
}

UTEST(pdb, chain_break_numbering_runs_on) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 1, 0, NULL}, {1, 2, 0, NULL}, {2, 3, 0, NULL}, {3, 4, 0, NULL}, {6, 5, 0, NULL}, {7, 6, 0, NULL}, {8, 7, 0, NULL}, {9, 8, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, ARRAY_SIZE(pick), alloc));
    expect_break_after_fourth(utest_result, &sys);

    // Without bonds the distance of the C and the N tells the same
    sys.bond.count = 0;
    md_array_shrink(sys.bond.pairs, 0);
    md_array_shrink(sys.bond.flags, 0);
    md_bond_build_connectivity(&sys.bond, sys.atom.count, alloc);
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_BACKBONE_BIT));
    expect_break_after_fourth(utest_result, &sys);

    md_vm_arena_destroy(alloc);
}

UTEST(pdb, chain_break_numbering_skips) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const peptide_pick_t pick[] = {
        {0, 1, 0, NULL}, {1, 2, 0, NULL}, {2, 3, 0, NULL}, {3, 4, 0, NULL}, {6, 7, 0, NULL}, {7, 8, 0, NULL}, {8, 9, 0, NULL}, {9, 10, 0, NULL},
    };
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(load_peptide(&sys, &state, pick, ARRAY_SIZE(pick), alloc));
    expect_break_after_fourth(utest_result, &sys);

    // Inferring again finds the same breaks, not more
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    expect_break_after_fourth(utest_result, &sys);

    md_vm_arena_destroy(alloc);
}

// 1LAF ends chain E at Asp 238, whose carboxyl group is whole (OXT): the arginine numbered 240 after it is a free
// amino acid bound to the protein, not the chain going on. 1k4r and 1a64 are whole (1a64 misses the first residues of its
// chains only, which leaves no break)
UTEST(pdb, chain_breaks_of_files) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb"), MD_PDB_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        EXPECT_EQ((size_t)0, count_components_with(&sys, MD_COMPONENT_FLAG_BREAK_BEFORE | MD_COMPONENT_FLAG_BREAK_AFTER));
        size_t asp = SIZE_MAX;
        for (size_t i = 0; i + 1 < sys.component.count; ++i) {
            if (md_component_seq_id(&sys.component, i) == 238 && md_component_seq_id(&sys.component, i + 1) == 240) asp = i;
        }
        ASSERT_NE(SIZE_MAX, asp);
        EXPECT_TRUE(str_eq(md_component_name(&sys.component, asp), STR_LIT("ASP")));
        EXPECT_TRUE(md_component_flags(&sys.component, asp) & MD_COMPONENT_FLAG_TERMINAL_END);
    }
    const str_t whole[] = {STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), STR_LIT(MD_UNITTEST_DATA_DIR "/1a64.pdb")};
    for (size_t f = 0; f < ARRAY_SIZE(whole); ++f) {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, whole[f], MD_PDB_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        EXPECT_EQ((size_t)0, count_components_with(&sys, MD_COMPONENT_FLAG_BREAK_BEFORE | MD_COMPONENT_FLAG_BREAK_AFTER));
    }
    md_vm_arena_destroy(alloc);
}
