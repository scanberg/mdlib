#include "utest.h"

#include <md_trr.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UTEST(trr, trajectory_i) {
    md_trajectory_i* traj = md_trr_trajectory_create(STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), md_get_heap_allocator(), MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    EXPECT_EQ(md_trajectory_num_atoms(traj), 6495);
    EXPECT_EQ(md_trajectory_num_frames(traj), 101);

    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_get_heap_allocator(), mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
    }

    md_free(md_get_heap_allocator(), mem_ptr, mem_size);
    md_trr_trajectory_free(traj);
}
