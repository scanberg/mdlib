#include "md_trajectory.h"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_file.h>
#include <core/md_str.h>

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

bool md_trajectory_load_frame(const struct md_trajectory* traj, int64_t frame_idx, md_trajectory_field_flags_t requested_fields, md_trajectory_data_t* write_target) {
    ASSERT(traj);
    ASSERT(traj->inst);

    bool result = false;
    const int64_t frame_size = traj->extract_frame(traj->inst, frame_idx, NULL);
    if (frame_size > 0) {
        void* frame_data = md_alloc(default_allocator, frame_size);
        const int64_t read_size = traj->extract_frame(traj->inst, frame_idx, frame_data);
        ASSERT(read_size == frame_size);
        result = traj->decode_frame(traj->inst, frame_data, frame_size, requested_fields, write_target);
        md_free(default_allocator, frame_data, frame_size);
    }

    return result;
}

#ifdef __cplusplus
}
#endif