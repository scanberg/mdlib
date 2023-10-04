#include "ubench.h"

#include <md_xtc.h>
#include <md_trajectory.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

UBENCH_EX(xtc, load) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    md_trajectory_i* traj = md_xtc_trajectory_create(STR(MD_BENCHMARK_DATA_DIR "/catalyst.xtc"), alloc);

    if (!traj) {
        printf("Failed to load trajectory\n");
        return;
    }
    
    const int64_t num_atoms  = md_trajectory_num_atoms(traj);
    const int64_t num_frames = md_trajectory_num_frames(traj);

    float* x = md_alloc(alloc, sizeof(float) * num_atoms);
    float* y = md_alloc(alloc, sizeof(float) * num_atoms);
    float* z = md_alloc(alloc, sizeof(float) * num_atoms);
    
    md_trajectory_frame_header_t frame_header;
    UBENCH_DO_BENCHMARK() {
        for (int64_t i = 0; i < num_frames; ++i) {
            md_trajectory_load_frame(traj, i, &frame_header, x, y, z);
        }
    }

    md_arena_allocator_destroy(alloc);
}
