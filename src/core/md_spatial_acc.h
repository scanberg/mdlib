#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_unitcell_t;

typedef struct md_spatial_acc_t {
    size_t num_elems;
    float* elem_x;
    float* elem_y;
    float* elem_z;
    uint32_t* elem_idx;

    size_t num_cells;
    uint32_t* cell_off;
    uint64_t  cell_mask[3][16];
    uint32_t  cell_dim[3];
    float     cell_ext[3];

    float G00, G11, G22;
    float H01, H02, H12;

    float A[3][3];
    float I[3][3];
    float origin[3];

    uint32_t flags;

    struct md_allocator_i* alloc;
} md_spatial_acc_t;

#ifdef __cplusplus
extern "C" {
#endif

// Callback signatures

// Callback for single point query.
typedef void (*md_spatial_acc_point_callback_t)(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user_param);

// Callback for pairwise interactions
typedef void (*md_spatial_acc_pair_callback_t)(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param);

// It is recommended to use a cell extent that is equal to the maximum search radius
void md_spatial_acc_init(md_spatial_acc_t* acc, const float* in_x, const float* in_y, const float* in_z, const int32_t* in_idx, size_t count, double cell_ext, const struct md_unitcell_t* unitcell);
void md_spatial_acc_free(md_spatial_acc_t* acc);

// --- INTERNAL PAIR TESTS ---

// Perform full N^2 test of points within the spatial acceleration structure for a supplied radius
// It is recommended that the radius <~ cell_ext, then a tight neighbor search is performed
void md_spatial_acc_for_each_internal_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param);

// Iterate over each point within the spatial acceleration structure within a 1-cell neighborhood (periodic if applicable)
void md_spatial_acc_for_each_internal_pair_in_neighboring_cells(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param);

// --- EXTERNAL PAIR TESTS ---

// Test external points against internal points within the spatial acceleration structure for a supplied cutoff
// The external points are not part of the spatial acceleration structure and will be represented in the callback as the 'i' indices and the internal points are the 'j' indices in the callback
void md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, const int32_t* ext_idx, size_t ext_count, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param);

// Perform a spatial query for points within the spatial acceleration structure within a bounding box defined by center and extent
void md_spatial_acc_for_each_point_in_aabb(const md_spatial_acc_t* acc, const double aabb_min[3], const double aabb_max[3], md_spatial_acc_point_callback_t callback, void* user_param);

void md_spatial_acc_for_each_point_in_sphere(const md_spatial_acc_t* acc, double center[3], double radius, md_spatial_acc_point_callback_t callback, void* user_param);

#if 0
// Iterate over external points against points within the spatial acceleration structure in neighboring cells (1-cell neighborhood)
bool md_spatial_acc_for_each_external_point_in_neighboring_cells(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, size_t ext_count, md_spatial_acc_pair_callback_t callback, void* user_param);
#endif

// --- HELPER FUNCTIONS ---

// Helper functions for partial functionality
static inline void md_spatial_acc_cell_range(uint32_t out_cell_range[2], const md_spatial_acc_t* acc, size_t cell_idx) {
	out_cell_range[0] = 0;
	out_cell_range[1] = 0;
    if (cell_idx < acc->num_cells) {
        out_cell_range[0] = acc->cell_off[cell_idx];
        out_cell_range[1] = acc->cell_off[cell_idx + 1];
    }
}

#ifdef __cplusplus
}
#endif
