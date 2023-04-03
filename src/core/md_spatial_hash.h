#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_array.h>
#include <core/md_vec_math.h>

/*
 *  This is an attempt of a spatial hashing structure intended for accelerated spatial queries.
 *  It makes some assumptions as it is not intended to be a general purpose spatial hash:
 *  - Only 3-dimensions (You could surely use it for 2D and 1D, but it would still do some redundant work)
 *  - Each cell can only hold up to 1023 entries, this is because the offset and length is packed into a single 32-bit integer.
 *      This implies the input data has a certain spatial distribution and the cell size is not too small.
 *  - The internal coordinates are compressed as 3x10-bit integers, which means there is some loss in precision. (Could change over time)
 *  - The cell size is currently fixed for simplicity. Could certainly change over time.
 *  - Internally works for periodic data as well as non-periodic
 * 
 *  Things to improve:
 *  Simlify implementation. The current version of using fixed cell size makes it a bit more complicated in the periodic case than it needs to be.
 *  Query should be handled with an external state object which can be used to iterate. The current callback is just cumersome and fugly.
 *  There should be an init function which should accept x,y,z + indices and xyz + indices
 * 
 */

struct md_allocator_i;
struct md_bitfield_t;

// Callback function for when a point is found within the search space of the query
// The return value (bool) signifies if it should continue its search (true) or if it should early exit (false).
typedef bool (*md_spatial_hash_iterator_fn)(uint32_t idx, vec3_t coord, void* user_param);

typedef struct md_spatial_hash_elem_t {
    vec3_t coord;
    uint32_t idx;
} md_spatial_hash_elem_t;

typedef struct md_spatial_hash_t {
    uint32_t*               cells;
    uint32_t*               coords;     // Compressed coordinates 10 bits precision relative to cell
    uint32_t*               indices;    // original given index
    struct md_allocator_i*  alloc;
    int32_t cell_min[3];
    uint32_t coord_count;
    int32_t cell_dim[3];
    uint32_t magic;
    vec4_t pbc_ext;
} md_spatial_hash_t;

typedef struct md_spatial_acc_t md_spatial_acc_t;

#ifdef __cplusplus
extern "C" {
#endif
    
md_spatial_acc_t* md_spatial_acc_create_vec3(const vec3_t* xyz, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);
md_spatial_acc_t* md_spatial_acc_create_vec3_indexed(const vec3_t* xyz, const int32_t* indices, int64_t index_count, vec3_t pbc_ext, struct md_allocator_i* alloc);
md_spatial_acc_t* md_spatial_acc_create_soa(const float* x, const float* y, const float* z, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);
md_spatial_acc_t* md_spatial_acc_create_soa_indexed(const float* x, const float* y, const float* z, const int32_t* indices, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);

void              md_spatial_acc_free(md_spatial_acc_t* acc);

// Initialize a spatial hash structure with a set of coordinates given by array of vec3_t (xyz) or separate arrays (x,y,z).
// pbc_ext is optional and supplies the periodic extent for each axis. To disable periodicity, supply (0,0,0).
bool md_spatial_hash_init(md_spatial_hash_t* spatial_hash, const vec3_t* xyz, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);
bool md_spatial_hash_init_indexed(md_spatial_hash_t* spatial_hash, const vec3_t* xyz, const int32_t* indices, int64_t index_count, vec3_t pbc_ext, struct md_allocator_i* alloc);

bool md_spatial_hash_init_soa(md_spatial_hash_t* spatial_hash, const float* x, const float* y, const float* z, int64_t count, vec3_t pbc_ext, struct md_allocator_i* alloc);
bool md_spatial_hash_init_indexed_soa(md_spatial_hash_t* spatial_hash, const float* x, const float* y, const float* z, const int32_t* indices, int64_t index_count, vec3_t pbc_ext, struct md_allocator_i* alloc);

bool md_spatial_hash_free(md_spatial_hash_t* spatial_hash);

// Perform a spatial query for a given position + radius in a periodic domain given by pbc_min and pbc_max.
// This will also include periodic occurrences of points accross the periodic boundries.
bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param);

// Get a list of indices which fall within the search space (pos + radius)
// Writes directly to the supplied buffer and will return the number of indices written.
int64_t md_spatial_hash_query_idx(int32_t* buf_idx, int64_t buf_cap, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius);
void    md_spatial_hash_query_bits(struct md_bitfield_t* bf, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius);

#ifdef __cplusplus
}
#endif
