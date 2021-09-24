#include "md_trajectory.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_file.h>
#include <core/md_str.h>

#include <string.h>

bool md_trajectory_default_load_frame(md_trajectory_i* traj, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(traj);

    bool result = true;
    const int64_t frame_size = traj->extract_frame_data(traj->inst, frame_idx, NULL);
    if (frame_size > 0) {
        // This is a borderline case if one should use the default_temp_allocator as the raw frame size could potentially be several megabytes...
        void* frame_data = md_alloc(default_allocator, frame_size);
        const int64_t read_size = traj->extract_frame_data(traj->inst, frame_idx, frame_data);
        ASSERT(read_size == frame_size);

        result = traj->decode_frame_data(traj->inst, frame_data, frame_size, header, x, y, z);

        md_free(default_allocator, frame_data, frame_size);
    }

    return result;
}
