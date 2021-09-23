#ifndef _MD_TRAJ_H_
#define _MD_TRAJ_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_trajectory_o;

typedef struct md_trajectory_frame_header {
	int32_t num_atoms;
	int32_t step;
	double  timestamp;
	float	box[3][3];
} md_trajectory_frame_header_t;

typedef struct md_trajectory_i {
	struct md_trajectory_o* inst; // Opaque trajectory data

	int64_t num_atoms;
	int64_t num_frames;
	int64_t max_frame_data_size;

	// Returns size in bytes of frame, frame_data_ptr is optional and if supplied, the frame data will be written to it.
	int64_t (*extract_frame_data)(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr);

	// Decodes the raw frame data
	bool (*decode_frame_header)(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_frame_header_t* header);
	bool (*decode_frame_coords)(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, float* x, float* y, float* z, int64_t num_coord);
} md_trajectory_i;

// This is the easymode version which reads and decodes the frame in one function call.
bool md_trajectory_load_frame(const struct md_trajectory_i* traj, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z, int64_t num_coord);

#ifdef __cplusplus
}
#endif

#endif
