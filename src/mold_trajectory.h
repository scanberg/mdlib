#ifndef _MOLD_TRAJ_H_
#define _MOLD_TRAJ_H_

#include "mold.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MOLD_TRAJECTORY_SUCCESS			 0
#define MOLD_TRAJECTORY_INVALID_FILE	-1
#define MOLD_TRAJECTORY_CORRUPT			-2

typedef struct mold_trajectory_cache_desc mold_trajectory_cache_desc;
typedef struct mold_trajectory_load_frame_target_desc mold_trajectory_load_frame_target_desc;
typedef struct mold_trajectory_header mold_trajectory_header;
typedef struct mold_trajectory mold_trajectory;

struct mold_trajectory_load_frame_target_desc {
	float* x;	 // storage for num_atoms floats,	x coordinates of atoms
	float* y;	 // storage for num_atoms floats,	y coordinates of atoms
	float* z;	 // storage for num_atoms floats,	z coordinates of atoms
	float* time; // storage for 1 float,			time stamp of frame
	float* box;	 // storage for 9 floats,			3x3 column major matrix defining the simulation box / unitcell
};

struct mold_trajectory_header {
	uint32_t num_atoms;
	uint32_t num_frames;
	float start_time;
	float frame_delta_time;
};

struct mold_trajectory;

// Initializes a trajectory given a path (filename) to the trajectory file.
// The cache descriptor (cache_desc) is optional but recommended to give the trajectory some internal cache space to manage frames.
mold_error mold_trajectory_init(mold_trajectory* traj, const char* filename, uint64_t cache_byte_size);
mold_error mold_trajectory_free(mold_trajectory* traj);

mold_error mold_trajectory_read_header(const mold_trajectory* traj, mold_trajectory_header* header);

// Loads frame data
// The target descriptor (target_desc) is optional and if set to NULL the function pre-fetchs the data and stores it in its internal cache (if supplied during initialization).
mold_error mold_trajectory_load_frame(mold_trajectory* traj, uint32_t frame_idx, mold_trajectory_load_frame_target_desc* target_desc);

#ifdef __cplusplus
}
#endif

#endif