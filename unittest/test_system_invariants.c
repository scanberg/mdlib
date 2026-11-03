#include "system_invariants.h"

#include <md_pdb.h>
#include <md_gro.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <string.h>

typedef bool (*sys_load_fn)(md_system_t*, md_system_state_t*, str_t);

// Load a file, run full inference, and assert the structural invariants hold.
static void check_file(int* utest_result, sys_load_fn load, str_t path) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };

    if (!load(&sys, &state, path)) {
        printf("failed to load '" STR_FMT "'\n", STR_ARG(path));
        *utest_result = UTEST_TEST_FAILURE;
        md_temp_end(temp);
        return;
    }

    md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL, NULL);
    expect_valid_system(utest_result, __FILE__, __LINE__, &sys);

    md_system_free(&sys);
    md_temp_end(temp);
}

static bool load_pdb(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_pdb_system_init_from_file(sys, state, path, MD_PDB_OPTION_NONE, NULL);
}

static bool load_gro(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_gro_system_init_from_file(sys, state, path, NULL);
}

static bool load_lammps(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_lammps_system_init_from_file(sys, state, path, NULL);
}

UTEST(system_invariants, pdb) {
    check_file(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"));
    check_file(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_file(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"));
}

UTEST(system_invariants, mmcif) {
    check_file(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"));
    check_file(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    check_file(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
}

UTEST(system_invariants, gro) {
    check_file(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    check_file(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_file(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"));
}

UTEST(system_invariants, lammps) {
    check_file(utest_result, load_lammps, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
}

// A copy is the same system in another allocator: valid, and inferred on it gives what inferred on the original gave
UTEST(system_invariants, copy) {
    md_allocator_i* a1 = md_vm_arena_create(GIGABYTES(1));
    md_allocator_i* a2 = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = a1 };
    md_system_state_t state = { .alloc = a1 };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE, NULL));

    // The core, copied and inferred on, as viamd does in the background
    md_system_t core = { .alloc = a2 };
    ASSERT_TRUE(md_system_copy(&core, &sys));
    EXPECT_EQ(sys.topology_version, core.topology_version);
    EXPECT_EQ(sys.atom.count, core.atom.count);
    EXPECT_NE((void*)sys.atom.type_idx, (void*)core.atom.type_idx);
    md_system_state_t core_state = { .alloc = a2 };
    md_system_state_copy(&core_state, &state);
    ASSERT_TRUE(md_util_system_infer(&core, &core_state, MD_UTIL_INFER_ALL, NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL, NULL));
    expect_valid_system(utest_result, __FILE__, __LINE__, &core);
    EXPECT_EQ(sys.bond.count, core.bond.count);
    EXPECT_EQ(0, memcmp(sys.bond.pairs, core.bond.pairs, sys.bond.count * sizeof(md_atom_pair_t)));
    EXPECT_EQ(0, memcmp(sys.residue.flags, core.residue.flags, sys.residue.count * sizeof(md_residue_flags_t)));
    EXPECT_EQ(sys.part.count, core.part.count);
    EXPECT_EQ(md_system_molecule_count(&sys), md_system_molecule_count(&core));

    // And the inferred system copied: everything there, nothing shared
    md_allocator_i* a3 = md_vm_arena_create(GIGABYTES(1));
    md_system_t full = { .alloc = a3 };
    ASSERT_TRUE(md_system_copy(&full, &sys));
    expect_valid_system(utest_result, __FILE__, __LINE__, &full);
    EXPECT_EQ(sys.protein_backbone.segment.count, full.protein_backbone.segment.count);
    EXPECT_EQ(md_array_size(sys.ring.rings), md_array_size(full.ring.rings));
    EXPECT_EQ(sys.reference.num_atoms, full.reference.num_atoms);
    EXPECT_NE((void*)sys.reference.xyz, (void*)full.reference.xyz);
    for (size_t t = 0; t < sys.part.type.count; ++t) {
        EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, t), md_system_part_type_description(&full, t)));
    }

    md_vm_arena_destroy(a1);
    md_vm_arena_destroy(a2);
    md_vm_arena_destroy(a3);
}
