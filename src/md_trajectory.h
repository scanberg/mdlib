#ifndef _MD_TRAJ_H_
#define _MD_TRAJ_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator;
struct md_trajectory_o;

enum {
	MD_TRAJ_FIELD_XYZ		= 1,
	MD_TRAJ_FIELD_BOX		= 2,
	MD_TRAJ_FIELD_TIMESTAMP = 4,
};

typedef uint32_t md_trajectory_field_flags_t;

typedef struct md_trajectory_box {
	float basis[3][3];
} md_trajectory_box_t;

typedef struct md_trajectory_data {
	int64_t num_atoms;
	float *x;
	float *y;
	float *z;

	md_trajectory_box_t *box;
	double *timestamp;
	md_trajectory_field_flags_t written_fields;
} md_trajectory_data_t;

typedef struct md_trajectory {
	struct md_trajectory_o* inst; // Opaque trajectory data

	int64_t num_atoms;
	int64_t num_frames;

	// Returns size in bytes of frame, frame_data_ptr is optional and is the destination to write the frame data to.
	int64_t (*extract_frame)(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr);

	// Decodes the raw frame data
	bool (*decode_frame)(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_field_flags_t requested_fields, md_trajectory_data_t* write_target);
} md_trajectory_i;

// This is the easymode version which reads and decodes the frame in one function call.
bool md_trajectory_load_frame(const struct md_trajectory* traj, int64_t frame_idx, md_trajectory_field_flags_t requested_fields, md_trajectory_data_t* write_target);



#ifdef __cplusplus
}
#endif

#endif