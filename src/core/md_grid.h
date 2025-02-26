#pragma once

#include <stddef.h>

// The grid data is assumed to be given in Z,Y,X order (e.g. data[Z][Y][X])
typedef struct md_grid_t {
	float* data;
	int   dim[3];
	float origin[3];
	float step_x[3];
	float step_y[3];
	float step_z[3];
} md_grid_t;

#ifdef __cplusplus
extern "C" {
#endif

// Procedures for extracting 4x4 column major transformation matrices from a grid
void md_grid_extract_model_to_world(float* out_mat, const md_grid_t* grid);
void md_grid_extract_world_to_model(float* out_mat, const md_grid_t* grid);
void md_grid_extract_index_to_world(float* out_mat, const md_grid_t* grid);

static inline size_t md_grid_num_points(const md_grid_t* grid) {
	if (grid) {
		size_t dim = (size_t)grid->dim[0] * (size_t)grid->dim[1] * (size_t)grid->dim[2];
		return dim;
	}
	return 0;
}

void md_grid_extract_points(float* out_xyz, const md_grid_t* grid);

#ifdef __cplusplus
}
#endif
