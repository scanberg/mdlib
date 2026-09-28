#include "utest.h"

#include <md_ase_traj.h>
#include <md_system.h>
#include <core/md_allocator.h>

#include "run_check.h"

#define ASE_RUN STR_LIT("run/ase")
#define ASE_FILE(name) STR_LIT(MD_ASE_TRAJ_TEST_DATA_DIR "/" name)

UTEST(ase_traj, fixed_atoms_and_time) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_system_t sys = {.alloc = alloc};
    md_system_state_t initial = {.alloc = alloc};
    const str_t path = ASE_FILE("ase_fixed.traj");

    ASSERT_TRUE(md_ase_traj_system_init_from_file(&sys, &initial, path));
    ASSERT_EQ(3u, sys.atom.count);
    EXPECT_EQ(1, md_atom_atomic_number(&sys.atom, 0));
    EXPECT_EQ(8, md_atom_atomic_number(&sys.atom, 2));
    EXPECT_NEAR(0.9, initial.xyz[1].x, 1.0e-6);
    EXPECT_NEAR(10.0, initial.unitcell.x, 1.0e-6);

    ASSERT_TRUE(md_ase_traj_system_publish_run(&sys, path, ASE_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    ASSERT_EQ(3u, run_num_frames(&sys, ASE_RUN));
    ASSERT_EQ(3u, run_num_atoms(&sys, ASE_RUN));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/ase/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    EXPECT_NEAR(0.25, ((const double*)time->data)[1], 1.0e-6);
    EXPECT_NEAR(0.50, ((const double*)time->data)[2], 1.0e-6);

    md_system_state_t frame = {.alloc = alloc};
    ASSERT_TRUE(md_system_state_init(&frame, sys.atom.count));
    ASSERT_TRUE(run_extract_one(&frame, &sys, ASE_RUN, 2));
    EXPECT_NEAR(0.3, frame.xyz[0].x, 1.0e-6);
    EXPECT_NEAR(1.2, frame.xyz[1].x, 1.0e-6);
    EXPECT_NEAR(0.8, frame.xyz[2].y, 1.0e-6);
    EXPECT_NEAR(12.0, frame.unitcell.x, 1.0e-6);
    EXPECT_NEAR(10.0, frame.unitcell.y, 1.0e-6);
    EXPECT_NEAR(1.5, frame.unitcell.xy, 1.0e-6);
    EXPECT_NEAR(-0.5, frame.unitcell.xz, 1.0e-6);
    EXPECT_NEAR(0.75, frame.unitcell.yz, 1.0e-6);
    EXPECT_FALSE(run_extract_one(&frame, &sys, ASE_RUN, 3));

    md_system_state_free(&frame);
    md_system_state_free(&initial);
    md_system_free(&sys);
}

UTEST(ase_traj, float32_positions) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_system_t sys = {.alloc = alloc};
    md_system_state_t initial = {.alloc = alloc};
    const str_t path = ASE_FILE("ase_float32.traj");

    ASSERT_TRUE(md_ase_traj_system_init_from_file(&sys, &initial, path));
    ASSERT_EQ(2u, sys.atom.count);
    EXPECT_EQ(0, md_atom_atomic_number(&sys.atom, 0));
    EXPECT_EQ(8, md_atom_atomic_number(&sys.atom, 1));
    EXPECT_NEAR(4.0, initial.xyz[1].x, 1.0e-6);
    EXPECT_EQ(0, initial.unitcell.flags & MD_UNITCELL_PBC_ALL);

    ASSERT_TRUE(md_ase_traj_system_publish_run(&sys, path, ASE_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_EQ(1u, run_num_frames(&sys, ASE_RUN));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/ase/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_TRUE(md_unit_is_none(time->unit));
    EXPECT_EQ(0.0, ((const double*)time->data)[0]);

    md_system_state_t frame = {.alloc = alloc};
    ASSERT_TRUE(md_system_state_init(&frame, sys.atom.count));
    ASSERT_TRUE(run_extract_one(&frame, &sys, ASE_RUN, 0));
    EXPECT_NEAR(4.0, frame.xyz[1].x, 1.0e-6);
    EXPECT_NEAR(6.0, frame.xyz[1].z, 1.0e-6);

    md_system_state_free(&frame);
    md_system_state_free(&initial);
    md_system_free(&sys);
}

UTEST(ase_traj, rejects_unsupported_input) {
    const str_t invalid[] = {
        ASE_FILE("ase_rotated_cell.traj"),
        ASE_FILE("ase_partial_pbc.traj"),
        ASE_FILE("make_fixtures.py"),
    };
    for (size_t i = 0; i < ARRAY_SIZE(invalid); ++i) {
        md_system_t sys = {.alloc = md_get_heap_allocator()};
        md_system_state_t state = {.alloc = md_get_heap_allocator()};
        EXPECT_FALSE(md_ase_traj_system_init_from_file(&sys, &state, invalid[i]));
        md_system_state_free(&state);
        md_system_free(&sys);
    }
}

UTEST(ase_traj, changing_atoms_rejected_when_published) {
    md_system_t sys = {.alloc = md_get_heap_allocator()};
    md_system_state_t state = {.alloc = md_get_heap_allocator()};
    const str_t path = ASE_FILE("ase_changing_atoms.traj");
    // Initialization reads only frame 0; run publication validates every frame.
    ASSERT_TRUE(md_ase_traj_system_init_from_file(&sys, &state, path));
    EXPECT_FALSE(md_ase_traj_system_publish_run(&sys, path, ASE_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_EQ(0u, run_num_frames(&sys, ASE_RUN));
    md_system_state_free(&state);
    md_system_free(&sys);
}
