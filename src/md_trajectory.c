#include "md_trajectory.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_file.h>
#include <core/md_str.h>

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

bool md_trajectory_load_frame(const md_trajectory_i* traj, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z, int64_t num_coord) {
    ASSERT(traj);
    ASSERT(traj->inst);

    bool result = true;
    const int64_t frame_size = traj->extract_frame_data(traj->inst, frame_idx, NULL);
    if (frame_size > 0) {
        void* frame_data = md_alloc(default_allocator, frame_size);
        const int64_t read_size = traj->extract_frame_data(traj->inst, frame_idx, frame_data);
        ASSERT(read_size == frame_size);
        bool read_header = header != NULL;
        bool read_coords = x != NULL || y != NULL || z != NULL;
        result |= (read_header ? traj->decode_frame_header(traj->inst, frame_data, frame_size, header) : true) &&
                 (read_coords ? traj->decode_frame_coords(traj->inst, frame_data, frame_size, x, y, z, num_coord) : true);
        md_free(default_allocator, frame_data, frame_size);
    }

    return result;
}

#ifdef __cplusplus
}
#endif