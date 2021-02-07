#ifndef _MD_TRAJ_H_
#define _MD_TRAJ_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

	/*
// This is used to define the target storage location for loading frame data
struct md_trajectory_frame_data {
	float* x;	 // target storage for num_atoms floats,	x coordinates of atoms
	float* y;	 // target storage for num_atoms floats,	y coordinates of atoms
	float* z;	 // target storage for num_atoms floats,	z coordinates of atoms
	float time;
	float box[3][3];
};

struct md_trajectory_header {
	uint32_t num_atoms;
	uint32_t num_frames;
	float start_time;
	float frame_delta_time;
};

struct md_trajectory_o;

struct md_trajectory_i {
	struct md_trajectory_o* inst;
	bool (*init_from_file)(struct md_trajectory* traj, const char* filename);
	bool (*init_from_mem)(struct md_trajectory* traj)
};

bool md_trajectory_init_from_file(struct md_trajectory* traj, const char* filename, struct md_trajectory_config* config);
bool md_trajectory_init_from_mem(struct md_trajectory* traj, const void* memory, struct md_trajectory_config* config);
void md_trajectory_free(struct md_trajectory* traj);

// Loads frame data
// target is optional and if set to NULL the function pre-fetchs the data and stores it in its internal cache.
bool md_trajectory_load_frame(struct md_trajectory* traj, uint32_t frame_idx, struct md_trajectory_frame_data* target);
*/
#ifdef __cplusplus
}
#endif

#endif