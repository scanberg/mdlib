#include "utest.h"
#include <string.h>

#include <md_xtc.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_file.h>

UTEST(xtc, trajectory_i) {
    const str_t path = make_cstr(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_trajectory_i traj = {0};

    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(16));
    ASSERT_TRUE(md_xtc_trajectory_open(&traj, path, alloc));

    EXPECT_EQ(md_trajectory_num_atoms(&traj), 1336);
    EXPECT_EQ(md_trajectory_num_frames(&traj), 501);

    const int64_t mem_size = md_trajectory_num_atoms(&traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(alloc, mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(&traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(&traj) * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < md_trajectory_num_frames(&traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(&traj, i, &header, x, y, z));
    }

    md_xtc_trajectory_close(&traj);

    md_arena_allocator_destroy(alloc);
}
