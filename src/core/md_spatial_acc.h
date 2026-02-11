#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_simd.h>

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
    uint32_t  cell_min[3];
    uint32_t  cell_max[3];
    uint32_t  cell_dim[3];
    float     inv_cell_ext[3];

    float G00, G11, G22;
    float H01, H02, H12;

    float I[3][3];
    // Origin offset (added to input coords before multiplying with I) to maintain same fractional frame as build
    float origin[3];

    uint32_t flags;

    struct md_allocator_i* alloc;
} md_spatial_acc_t;

#ifdef __cplusplus
extern "C" {
#endif

// Callback for cells
// The coordinates are fractional coordinates
typedef void (*md_spatial_acc_cell_callback_t)(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_elems, void* user_param);

// Callback for pairwise interactions
typedef void (*md_spatial_acc_pair_callback_t)(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param);

// It is recommended to use a cell extent that is equal to the maximum search radius
void md_spatial_acc_init(md_spatial_acc_t* acc, const float* in_x, const float* in_y, const float* in_z, const int32_t* in_idx, size_t count, double cell_ext, const struct md_unitcell_t* unitcell);
void md_spatial_acc_free(md_spatial_acc_t* acc);

// --- INTERNAL PAIR TESTS ---

// Perform full N^2 test of points within the spatial acceleration structure for a supplied radius
// It is recommended that the radius <~ cell_ext, then a tight neighbor search is performed
bool md_spatial_acc_for_each_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param);

// Iterate over each point within the spatial acceleration structure within a 1-cell neighborhood (periodic if applicable)
void md_spatial_acc_for_each_pair_in_neighboring_cells(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param);

// --- EXTERNAL PAIR TESTS ---

// Iterate over external points against points within the spatial acceleration structure for a supplied cutoff
// The external points are not part of the spatial acceleration structure and will be represented in the callback as the 'i' indices and the internal points are the 'j' indices
bool md_spatial_acc_for_each_external_point_within_cutoff(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, const int32_t* ext_idx, size_t ext_count, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param);

#if 0
// Iterate over external points against points within the spatial acceleration structure in neighboring cells (1-cell neighborhood)
bool md_spatial_acc_for_each_external_point_in_neighboring_cells(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, size_t ext_count, md_spatial_acc_pair_callback_t callback, void* user_param);
#endif

// --- HELPER FUNCTIONS ---

// Helper functions for partial functionality
static inline bool md_spatial_acc_cell_range(uint32_t* cell_beg, uint32_t* cell_end, const md_spatial_acc_t* acc, size_t cell_idx) {
    if (cell_idx >= acc->num_cells)
        return false;
    *cell_beg = acc->cell_off[cell_idx];
    *cell_end = acc->cell_off[cell_idx + 1];
    return true;
}

#ifdef __cplusplus
}
#endif
