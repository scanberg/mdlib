#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_vec_math.h>

struct md_allocator_i;
struct md_spatial_hash_cell_t;

typedef bool (*md_spatial_hash_iterator_fn)(uint32_t idx, vec3_t coord, void* user_param);

typedef struct md_spatial_hash_t {
    struct md_spatial_hash_cell_t*  cells;
    uint32_t*                       coords;     // Compressed coordinates 10 bits precision relative to cell
    uint32_t*                       indices;    // original given index
    struct md_allocator_i*          alloc;
    int32_t cell_min[3];
    uint32_t coord_count;
    int32_t cell_dim[3];
    uint32_t magic;
    vec4_t pbc_ext;
} md_spatial_hash_t;

#ifdef __cplusplus
extern "C" {
#endif

// Initialize a spatial hash structure with a set of coordinates given by array of vec3_t (xyz) or separate arrays (x,y,z).
// pbc_ext is optional and supplies the periodic extent for each axis. To disable periodicity, supply (0,0,0).
bool md_spatial_hash_init(md_spatial_hash_t* spatial_hash, const vec3_t* xyz, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);
bool md_spatial_hash_init_soa(md_spatial_hash_t* spatial_hash, const float* x, const float* y, const float* z, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);

bool md_spatial_hash_free(md_spatial_hash_t* spatial_hash);


// Perform a spatial query for a given position + radius.
// Will iterate over all points within the given position + radius and call the supplied function iter.
//bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param);

// Perform a spatial query for a given position + radius in a periodic domain given by pbc_min and pbc_max.
// This will also include periodic occurrences of points accross the periodic boundries.
bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param);

#ifdef __cplusplus
}
#endif
