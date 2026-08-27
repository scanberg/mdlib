#include "system_invariants.h"

#include <md_pdb.h>
#include <md_gro.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_util.h>

#include <core/md_allocator.h>

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

    md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL);
    expect_valid_system(utest_result, __FILE__, __LINE__, &sys);

    md_system_free(&sys);
    md_temp_end(temp);
}

static bool load_pdb(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_pdb_system_init_from_file(sys, state, path, MD_PDB_OPTION_NONE);
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
    check_file(utest_result, md_gro_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    check_file(utest_result, md_gro_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_file(utest_result, md_gro_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"));
}

UTEST(system_invariants, lammps) {
    check_file(utest_result, load_lammps, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
}
