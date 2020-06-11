#include "mold_trajectory.h"

#include <stdlib.h> // malloc

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#endif

#include <core/file.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef uint32_t mold_trajectory_type;

enum {
    TRAJECTORY_TYPE_UNDEFINED,
    TRAJECTORY_TYPE_PDB,
    TRAJECTORY_TYPE_XTC,
    TRAJECTORY_TYPE_DCD,
};

struct mold_trajectory {
    uint32_t magic;
    mold_trajectory_type type;
    uint32_t num_atoms;
    uint32_t num_frames;
    float start_time;       // Probably almost always zero?
    float frame_delta_time;
    mold_file* file;
    uint64_t* frame_offsets;
};

mold_error mold_trajectory_init(mold_trajectory* traj, const char* filename, const mold_trajectory_cache_desc* cache_desc) {

    // Determine trajectory type based on file ending

    // Read or Generate frame index cache

}

mold_error mold_trajectory_free(mold_trajectory* traj) {

}

mold_error mold_trajectory_read_header(const mold_trajectory* traj, mold_trajectory_header* header) {
    ASSERT(traj);
    ASSERT(header);
    header->num_atoms = traj->num_atoms;
    header->num_frames = traj->num_frames;
    header->start_time = traj->start_time;
    header->frame_delta_time = traj->frame_delta_time;
    return MOLD_TRAJECTORY_SUCCESS;
}

// Loads frame data
// The target descriptor (target_desc) is optional and if set to NULL the function pre-fetchs the data and stores it in its internal cache (if supplied during initialization).
mold_error mold_trajectory_load_frame(mold_trajectory* traj, uint32_t frame_idx, mold_trajectory_load_frame_target_desc* target_desc) {

}



#ifdef __cplusplus
}
#endif