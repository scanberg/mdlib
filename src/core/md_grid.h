#pragma once

#include <core/md_vec_math.h>
#include <stddef.h>

typedef struct md_grid_t {
    mat3_t orientation;     // Orientation of the grid in world space
    vec3_t origin;          // Origin of the grid in world units
    vec3_t spacing;         // Voxel spacing of the grid in world units
    int    dim[3];			// Dimensions of the grid in number of voxels
} md_grid_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline size_t md_grid_num_points(const md_grid_t* grid) {
    if (grid) {
        size_t dim = (size_t)grid->dim[0] * (size_t)grid->dim[1] * (size_t)grid->dim[2];
        return dim;
    }
    return 0;
}

// Extracts the points which make up the grid
// - out_xyz: The output array to store the points in world space coordinates (length = grid->dim[0] * grid->dim[1] * grid->dim[2] * 3)
// - grid: The grid to extract the points from
void md_grid_extract_points(float* out_xyz, const md_grid_t* grid);

// The convention used is that the model space does not include any scaling, only rotation and translation
mat4_t md_grid_model_to_world(const md_grid_t* grid);
mat4_t md_grid_world_to_model(const md_grid_t* grid);

// These includes scaling as well (by spacing)
mat4_t md_grid_index_to_world(const md_grid_t* grid);
mat4_t md_grid_world_to_index(const md_grid_t* grid);

static inline vec3_t md_grid_origin(const md_grid_t* grid) {
    ASSERT(grid);
    return grid->origin;
}

static inline mat3_t md_grid_orientation(const md_grid_t* grid) {
    ASSERT(grid);
    return grid->orientation;
}

static inline vec3_t md_grid_extent(const md_grid_t* grid) {
    ASSERT(grid);
    vec3_t extent = vec3_mul(grid->spacing, vec3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
    return extent;
}

static inline vec3_t md_grid_center(const md_grid_t* grid) {
    ASSERT(grid);
    vec3_t center = vec3_add(grid->origin, vec3_mul_f(md_grid_extent(grid), 0.5f));
    return center;
}

#ifdef __cplusplus
}
#endif
