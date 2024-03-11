#pragma once
#include <stddef.h>

struct md_allocator_i;

typedef struct md_gto_data_t {
	size_t count;
	float* x;
	float* y;
	float* z;
	float* neg_alpha;	// Negated Alpha term
	float* coeff;		// Coefficient: Should include normalization factors
	float* cutoff;		// Cutoff radial distance of the gto
	int* i;
	int* j;
	int* k;
	int* l;
} md_gto_data_t;

typedef struct md_grid_t {
	float* data;
	int   dim[3];
	float origin[3];
	float stepsize[3];
} md_grid_t;

#ifdef __cplusplus
extern "C" {
#endif

// Allocates gto data for given count
void md_gto_data_init(md_gto_data_t* gto, size_t count, struct md_allocator_i* alloc);
void md_gto_data_free(md_gto_data_t* gto, struct md_allocator_i* alloc);

// Evaluates a GTOs over a grid
// The grid data is assumed to be given in Z,Y,X order (e.g. data[Z][Y][X])
void md_gto_grid_evaluate    (md_grid_t* grid, const md_gto_data_t* gto);

// Evaluate over subportion of a grid
// - grid: The grid to evaluate a subportion of
// - grid_idx_off: Index offset for x,y,z
// - grid_idx_len: Index length for x,y,z
// - gto: The gto data to evaluate
void md_gto_grid_evaluate_sub(md_grid_t* grid, const int grid_idx_off[3], const int grid_idx_len[3], const md_gto_data_t* gto);


// Evaluates GTOs over a set of given XYZ coordinates
// 
// - out_psi: Array of values to write evaluated values to, should have length 'count'
// - xyz: Pointer to base addr of xyz packed coordinates to be evaluated, should have length 'count'
// - count: Number of coordinates to evaluate
// - stride: Stride in bytes between the given xyz coordinates. [OPTIONAL], a value of zero implies fully packed XYZXYZ... -> (12 bytes)
// - gto: input gtos to be evaluated for the supplied coordinates
void md_gto_xyz_evaluate(float* out_psi, const float* xyz, size_t count, size_t stride, const md_gto_data_t* gto);

// Compute the cutoff parameter within the gto data based on the given value
// Typically this could be somewhere around 1.0e-6
void md_gto_cutoff_compute(md_gto_data_t* gto, double value);

// Prune GTO data by creating a subset of GTOs which contribute to a spatial region (box)
// This is determined by the GTOs xyz + cutoff
void md_gto_spatial_prune(md_gto_data_t* out_gto, const md_gto_data_t* in_gto, float min_box[3], float max_box[3]);

#ifdef __cplusplus
}
#endif
