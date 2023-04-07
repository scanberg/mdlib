#include "ubench.h"

#include <md_xtc.h>
#include <md_trajectory.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

UBENCH_EX(xtc, load) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, MEGABYTES(1));
    md_trajectory_i* traj = md_xtc_trajectory_create(STR(MD_BENCHMARK_DATA_DIR "/catalyst.xtc"), alloc);
    md_trajectory_header_t header = {0};
    
    md_trajectory_get_header(traj, &header);

    float* x = md_alloc(alloc, sizeof(float) * header.num_atoms);
    float* y = md_alloc(alloc, sizeof(float) * header.num_atoms);
    float* z = md_alloc(alloc, sizeof(float) * header.num_atoms);
    
    md_trajectory_frame_header_t frame_header;
    UBENCH_DO_BENCHMARK() {
        const int64_t i = rand() % header.num_frames;
        md_trajectory_load_frame(traj, i, &frame_header, x, y, z);
    }

    md_arena_allocator_destroy(alloc);
}