#pragma once

#include <stdint.h>
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
    uint32_t  cell_min[3];
    uint32_t  cell_max[3];
    uint32_t  cell_dim[3];

    float G00, G11, G22;
    float H01, H02, H12;

    uint32_t flags;

    struct md_allocator_i* alloc;
} md_spatial_acc_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*md_spatial_acc_query_callback)(uint32_t i_idx, const uint32_t* j_idx, size_t j_len, void* user_param);

// It is recommended to use a cell extent that is equal to the maximum search radius
void md_spatial_acc_init(md_spatial_acc_t* acc, const float* in_x, const float* in_y, const float* in_z, size_t count, double cell_ext, const struct md_unitcell_t* unitcell);
void md_spatial_acc_free(md_spatial_acc_t* acc);

// Perform full N^2 test of points within the spatial acceleration structure for a supplied radius
// It is recommended that the radius < cell_ext, then a fast-path will be chosen
void md_spatial_acc_for_each_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_query_callback callback);

// Query all points within the spatial acceleration structure against a supplied list of points + radius
void md_spatial_acc_query(const md_spatial_acc_t* acc, const float* xyz, size_t xyz_stride, size_t count, double cutoff, md_spatial_acc_query_callback callback);

#ifdef __cplusplus
}
#endif
