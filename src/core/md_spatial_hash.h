#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_vec_math.h>

struct md_allocator_i;
struct md_spatial_hash_cell_t;

typedef bool (*md_spatial_hash_iterator_fn)(uint32_t idx, vec3_t coord, void* user_param);

typedef struct md_spatial_hash_t {
    struct md_spatial_hash_cell_t*  cells;
    float*                          coords;
    uint32_t*                       indices;    // original given index
    struct md_allocator_i*          alloc;
    vec3_t   cell_min;
    float    cell_ext;
    uint32_t cell_dim[3];
    uint32_t magic;
} md_spatial_hash_t;

// Initialization arguments
typedef struct md_spatial_hash_args_t {
    struct {
        int64_t count;
        int64_t stride; // The byte stride for the pointers. Set to zero if the data is packed
        const float* x;
        const float* y;
        const float* z;
    } coords;

    // Optional, set to zero if the implementation should decide the extent
    // Otherwise, the ideal value is roughly the intended search diameter / 3
    float cell_ext; 

    struct md_allocator_i* alloc;
    struct md_allocator_i* temp_alloc; // Optional, default is default_temp_allocator
} md_spatial_hash_args_t;

#ifdef __cplusplus
extern "C" {
#endif

bool md_spatial_hash_init(md_spatial_hash_t* spatial_hash, const md_spatial_hash_args_t* args);
bool md_spatial_hash_free(md_spatial_hash_t* spatial_hash);

// Perform a spatial query for a given position + radius.
// Will iterate over all points within the given position + radius and call the supplied function iter.
bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param);

// Perform a spatial query for a given position + radius in a periodic domain given by pbc_min and pbc_max.
// This will also include periodic occurrences of points accross the periodic boundries.
bool md_spatial_hash_query_periodic(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, vec3_t pbc_min, vec3_t pbc_max, md_spatial_hash_iterator_fn iter, void* user_param);

#ifdef __cplusplus
}
#endif
