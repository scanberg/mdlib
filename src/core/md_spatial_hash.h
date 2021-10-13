#pragma once

#include <stdint.h>
#include <stdbool.h>
#include "md_vec_math.h"

struct md_allocator_i;
struct md_spatial_hash_cell_t;
struct md_spatial_hash_coord_t;

typedef struct md_spatial_hash_t {
    struct md_spatial_hash_cell_t*  cells;
    struct md_spatial_hash_coord_t* coords;
    uint32_t*                        index;      // original given index
    vec3_t   cell_min;
    float    cell_ext;
    uint32_t cell_dim[3];
    uint64_t magic;
} md_spatial_hash_t;

typedef bool (*md_spatial_hash_iterator_fn)(uint32_t idx, vec3_t coord, void* user_param);

#ifdef __cplusplus
extern "C" {
#endif

bool md_spatial_hash_init(md_spatial_hash_t* spatial_hash, const float* x, const float* y, const float* z, int64_t count, float cell_ext, struct md_allocator_i* alloc);
bool md_spatial_hash_init_vec3(md_spatial_hash_t* spatial_hash, const vec3_t* xyz, int64_t count, float cell_ext, struct md_allocator_i* alloc);

bool md_spatial_hash_free(md_spatial_hash_t* spatial_hash, struct md_allocator_i* alloc);

// Perform a spatial query for a given position + radius.
// Will iterate over all points within the given position + radius and call the supplied function iter.
bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec4_t pos_rad, md_spatial_hash_iterator_fn iter, void* user_param);
//bool md_spatial_hash_query_n(const md_spatial_hash_t* spatial_hash, vec4_t* pos_rad, int64_t count, md_spatial_hash_iterator_fn iter, void* user_param);

#ifdef __cplusplus
}
#endif