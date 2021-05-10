#ifndef _MD_TRAJ_H_
#define _MD_TRAJ_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_trajectory_i md_trajectory_i;
typedef struct md_trajectory_header_t md_trajectory_header_t;
typedef struct md_trajectory_data_t md_trajectory_data_t;

typedef enum md_trajectory_field_flags_t {
	MD_TRAJ_FIELD_XYZ		= 1,
	MD_TRAJ_FIELD_BOX		= 2,
	MD_TRAJ_FIELD_TIMESTAMP = 4,
} md_trajectory_field_flags_t;

struct md_allocator_i;
struct md_trajectory_o;

struct md_trajectory_header_t {
	int64_t num_atoms;
	int64_t num_frames;
};

struct md_trajectory_data_t {
	int64_t num_atoms;
	float *x;
	float *y;
	float *z;

	float *box[3][3];
	double *timestamp;
	md_trajectory_field_flags_t written_fields;
};

struct md_trajectory_i {
	struct md_trajectory_o* inst; // Opaque trajectory data

	// For extracting the trajectory header
	bool (*extract_header)(struct md_trajectory_o* inst, struct md_trajectory_header_t* header);

	// Returns size in bytes of frame, frame_data_ptr is optional and is the destination to write the frame data to.
	int64_t (*extract_frame)(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr);

	// Decodes the raw frame data
	bool (*decode_frame)(const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_field_flags_t requested_fields, md_trajectory_data_t* write_target);
};

bool md_trajectory_extract_header(const struct md_trajectory_i* traj, struct md_trajectory_header_t* header);

// This is the easymode version which reads and decodes the frame in one function call.
bool md_trajectory_load_frame(const struct md_trajectory_i* traj, int64_t frame_idx, md_trajectory_field_flags_t requested_fields, md_trajectory_data_t* write_target);



#ifdef __cplusplus
}
#endif

#endif