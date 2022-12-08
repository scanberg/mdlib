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
    const str_t path = STR(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_trajectory_i* traj = md_xtc_trajectory_create(path, default_allocator);
    ASSERT_TRUE(traj);

    EXPECT_EQ(md_trajectory_num_atoms(traj), 1336);
    EXPECT_EQ(md_trajectory_num_frames(traj), 501);

    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(default_allocator, mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
    }

    md_free(default_allocator, mem_ptr, mem_size);
    md_xtc_trajectory_free(traj);
}
