#pragma once

#include <stdint.h>
#include <stdbool.h>



struct md_trajectory_o;

typedef struct md_trajectory_header_t {
	int64_t num_frames;
	int64_t num_atoms;
	int64_t max_frame_data_size; // This represents the maximum size of any frame which is extracted using extract_frame_data.
} md_trajectory_header_t;

typedef struct md_trajectory_frame_header_t {
	int64_t num_atoms;
	int64_t step;
	double  timestamp;
	float	box[3][3];
} md_trajectory_frame_header_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_trajectory_i {
	struct md_trajectory_o* inst; // Opaque trajectory data

	bool (*get_header)(struct md_trajectory_o* inst, md_trajectory_header_t* header);

	// --- EASY MODE ---
	// This the parameters header, x, y, z are optional and if you provide those pointers, it is assumed that they have enough space to hold the data
	bool (*load_frame)(const struct md_trajectory_i* traj, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z);

	// --- LOW LEVEL ---
	// Returns size in bytes of frame, frame_data_ptr is optional and if supplied, the frame data will be written to it.
	int64_t (*extract_frame_data)(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr);

	// Decodes the raw frame data
	bool (*decode_frame_data)(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z);
} md_trajectory_i;

// A convenient default option for load_frame if one does not want to do anything fancy when loading frame data.
bool md_trajectory_default_load_frame(const struct md_trajectory_i* traj, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z);

#ifdef __cplusplus
}
#endif

static inline bool md_trajectory_load_frame(const md_trajectory_i* traj, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
	return traj->load_frame(traj, frame_idx, header, x, y, z);
}

// Easy mode accessors
static inline int64_t md_trajectory_num_frames(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.num_frames;
	}
	return 0;
}

static inline int64_t md_trajectory_num_atoms(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.num_atoms;
	}
	return 0;
}

static inline int64_t md_trajectory_max_frame_data_size(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.max_frame_data_size;
	}
	return 0;
}


