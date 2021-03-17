#ifndef _MD_TRAJ_H_
#define _MD_TRAJ_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_trajectory_header_t {
	uint64_t num_atoms;
	uint64_t num_frames;
	float start_time;
	float frame_delta_time;
};

struct md_trajectory_o;

struct md_trajectory_i {
	struct md_trajectory_o* inst;
	bool (*read_header)(struct md_trajectory_i* traj, struct md_trajectory_header* header);
	bool (*load_frame_data)(struct md_trajectory_i* traj, uint64_t frame_idx, float* x, float* y, float* z, float* box[3][3], double* time);
};

#ifdef __cplusplus
}
#endif

#endif