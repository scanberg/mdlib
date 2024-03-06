#pragma once
#include <stddef.h>

typedef struct md_gto_data_t {
	size_t count;
	float* x;
	float* y;
	float* z;
	float* neg_alpha;
	float* coeff;
	float* cutoff;
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

void md_gto_grid_evaluate    (md_grid_t* grid, const md_gto_data_t* gto);
void md_gto_grid_evaluate_sub(md_grid_t* grid, const int grid_idx_min[3], const int grid_idx_max[3], const md_gto_data_t* gto);

// Compute the cutoff parameter within the gto data based on the given value
// Typically this could be somewhere around 1.0e-6
void md_gto_compute_cutoff(md_gto_data_t* gto, float value);

#ifdef __cplusplus
}
#endif
