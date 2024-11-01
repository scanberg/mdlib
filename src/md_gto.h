#pragma once
#include <stdint.h>
#include <stddef.h>

// Gaussian Type Orbital
// Are evaluated as f(x',y',z') = (x'-x)^i (y'-y)^j (z'-z)^k c exp(-a ((x'-x)^2 + (y'-y)^2 + (z'-z)^2))
// Where x' y' and z' are the observer coordinates we evaluate the function at
typedef struct md_gto_t {
	float x;
	float y;
	float z;
	float coeff;		// Baked coefficient (should include cartesian normalization factors)
	float alpha;		// Exponent alpha
	float cutoff;		// Radial cutoff
	uint8_t i, j, k, l;
	uint32_t _pad;
} md_gto_t;

// The grid data is assumed to be given in Z,Y,X order (e.g. data[Z][Y][X])
typedef struct md_grid_t {
	float* data;
	int   dim[3];
	float origin[3];
	float step_x[3];
	float step_y[3];
	float step_z[3];
} md_grid_t;

typedef enum {
	MD_GTO_EVAL_MODE_PSI = 0,
	MD_GTO_EVAL_MODE_PSI_SQUARED = 1,
} md_gto_eval_mode_t;

#ifdef __cplusplus
extern "C" {
#endif

// Evaluates GTOs over a grid
// - grid: The grid to evaluate a subportion of
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
void md_gto_grid_evaluate(md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode);

// Evaluates GTOs over a grid on the GPU and stores the result into a supplied volume
// - vol_tex: The texture handle to the volume
// - vol_dim: The dimensions of the volume
// - vol_step: The voxel spacing in world space length units
// - world_to_model: float[4][4] (col-major) transformation matrix to transform a point in world_space coordinates into the volumes model space (note not texture space, but a space which is rotated and translated such that the axes align with the volume and its origin is placed at (0,0,0))
// - index_to_world: float[4][4] (col-major) transformation matrix to transform a point in the volumes index coordinates [0, dim[ into world space coordinates
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
void md_gto_grid_evaluate_GPU(uint32_t vol_tex, const int vol_dim[3], const float vol_step[3], const float* world_to_model, const float* index_to_world, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode);

// This is malplaced at the moment, but this is for the moment, the best match in where to place the functionality
// Performs voronoi segmentation of the supplied volume to points with a supplied radius and accumulates the value of each voxel into the corresponding group of the closest point
// - out_group_values: Destination array holding the group values that are written to
// - cap_groups: Capacity of group array
// - vol_tex: The texture handle to the volume
// - vol_dim: The dimensions of the volume
// - vol_step: The voxel spacing in world space length units
// - world_to_model: float[4][4] (col-major) transformation matrix to transform a point in world_space coordinates into the volumes model space (note not texture space, but a space which is rotated and translated such that the axes align with the volume and its origin is placed at (0,0,0))
// - index_to_world: float[4][4] (col-major) transformation matrix to transform a point in the volumes index coordinates [0, dim[ into world space coordinates
// - point_xyzr: Point coordinates + radius, packed xyzrxyzrxyzr
// - point_group_idx: Point group index [0, num_groups-1]
// - num_points: Number of points
void md_gto_segment_and_attribute_to_groups_GPU(float* out_group_values, size_t cap_groups, uint32_t vol_tex, const int vol_dim[3], const float vol_step[3], const float* world_to_model, const float* index_to_world, const float* point_xyzr, const uint32_t* point_group_idx, size_t num_points);

// Evaluate GTOs over subportion of a grid
// - grid: The grid to evaluate a subportion of
// - grid_idx_off: Index offset for x,y,z
// - grid_idx_len: Index length for x,y,z
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
// @NOTE: It is strongly recommended that the evaluation occurs over a 8x8x8 blocks as this will get the fastpath
void md_gto_grid_evaluate_sub(md_grid_t* grid, const int grid_idx_off[3], const int grid_idx_len[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t eval_mode);

// Evaluates GTOs over a set of given XYZ coordinates
// - out_psi: Array of values to write evaluated values to, should have length 'num_xyz'
// - xyz: Pointer to base addr of xyz packed coordinates to be evaluated, should have length 'num_xyz'
// - num_xyz: Number of coordinates to evaluate
// - stride_xyz [OPTIONAL]: Stride in bytes between the given xyz coordinates. A value of zero implies fully packed XYZXYZ... -> (12 bytes)
// - gtos: input gtos to be evaluated for the supplied coordinates
// - num_gtos: Number of gtos
// - eval_mode: GTO evaluation mode
void md_gto_xyz_evaluate(float* out_psi, const float* xyz, size_t num_xyz, size_t stride_xyz, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t eval_mode);

// Evaluates GTOs over space 
//void md_gto_xyz_voronoi_evaluate(float* out_val, const float* xyz, size_t num_xyz, size_t stride_xyz, const md_gto_t* gtos, size_t)

// Compute the cutoff parameter within the supplied GTOs based on the given value
// Typically this could be somewhere around 1.0e-6
void md_gto_cutoff_compute(md_gto_t* gtos, size_t num_gtos, double value);

// Extracts a subset of gtos from an input array which overlap a given aabb with its radii of influence
size_t md_gto_aabb_test(md_gto_t* out_gtos, const float aabb_min[3], const float aabb_max[3], const md_gto_t* in_gtos, size_t num_gtos);

#ifdef __cplusplus
}
#endif
