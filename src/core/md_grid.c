#include <core/md_grid.h>
#include <core/md_common.h>

// Procedures for extracting 4x4 column major transformation matrices from a grid
void md_grid_extract_model_to_world(float* out_mat, const md_grid_t* grid) {
	ASSERT(grid);
	ASSERT(out_mat);
}

void md_grid_extract_world_to_model(float* out_mat, const md_grid_t* grid) {
	ASSERT(grid);
	ASSERT(out_mat);
}

void md_grid_extract_index_to_world(float* out_mat, const md_grid_t* grid) {
	ASSERT(grid);
	ASSERT(out_mat);
}

void md_grid_extract_points(float* out_xyz, const md_grid_t* grid) {
	ASSERT(grid);
	ASSERT(out_xyz);
}
