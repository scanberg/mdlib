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
	float coeff;		// Baked coefficient (should include normalization factors)
	float alpha;		// Exponent alpha
	float cutoff;		// Radial cutoff
	// The integer type here is arbitrary as we only need to store values 0-4 in reality.
	// uint16_t was chosen to pad the struct to 32 bytes in size
	uint16_t i, j, k, l;
} md_gto_t;

// The grid data is assumed to be given in Z,Y,X order (e.g. data[Z][Y][X])
typedef struct md_grid_t {
	float* data;
	int   dim[3];
	float origin[3];
	float stepsize[3];
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
