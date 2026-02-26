#pragma once

#include <core/md_vec_math.h>

#include <stdint.h>
#include <stdbool.h>

/*
 *  This is an attempt of a spatial hashing structure intended for accelerated spatial queries.
 *  It makes some assumptions as it is not intended to be a general purpose spatial hash:
 *  - Only 3-dimensions (You could surely use it for 2D and 1D, but it would still do some redundant work)
 *  - The cell size is currently fixed for simplicity. Could certainly change over time.
 *  - Internally works for periodic data as well as non-periodic
 * 
 *  Things to improve:
 *  Simplify implementation. The current version of using fixed cell size makes it a bit more complicated in the periodic case than it needs to be.
 *  Query should be handled with an external state object which can be used to iterate. The current callback is just cumersome and fugly.
 *  There should be an init function which should accept x,y,z + indices and xyz + indices
 * 
 */

struct md_allocator_i;
struct md_bitfield_t;
struct md_unitcell_t;

typedef struct md_spatial_hash_t md_spatial_hash_t;

typedef struct md_spatial_hash_elem_t {
    ALIGNAS(16) vec3_t xyz;
    uint32_t idx;
} md_spatial_hash_elem_t;

// Callback function for when a point is found within the search space of the query
// The return value (bool) signifies if it should continue its search (true) or if it should early exit (false).
typedef bool (*md_spatial_hash_iterator_fn)(const md_spatial_hash_elem_t* elem, void* user_param);

// This is a bit quirky, but it better maps to how the elements are processed internally (in batches).
// elem_mask contains the set bits of elements which passed the test for the query.
// To iterate over the elem_array, just find the set bits and clear them. Or just do a popcount for example if counting occurrences.
typedef bool (*md_spatial_hash_batch_iter_fn)(const md_spatial_hash_elem_t* elem_ptr, int elem_mask, void* user_param);
typedef bool (*md_spatial_hash_n2_batch_iter_fn)(const md_spatial_hash_elem_t* elem_ptr, md_256 d2, int elem_mask, size_t i, void* user_param);

#ifdef __cplusplus
extern "C" {
#endif
    
#if 0
typedef struct md_spatial_acc_t md_spatial_acc_t;
// Create spatial acceleration structures from a set of coordinates as input. Either supplied as vec3 or separate x, y, z streams with optional indices.
// If indices are supplied, then the count represents the number of indices which are supplied, otherwise count represents the number of coordinates to read directly
// as packed array elements.
// unit_cell is optional and represents the periodic bounds of the system (if applicable).
md_spatial_acc_t* md_spatial_acc_create_vec3(const vec3_t* xyz, const int32_t* indices, int64_t count, const struct md_unit_cell_t* unit_cell, struct md_allocator_i* alloc);
md_spatial_acc_t* md_spatial_acc_create_soa (const float* x, const float* y, const float* z, const int32_t* indices, int64_t count, const struct md_unit_cell_t* unit_cell, struct md_allocator_i* alloc);

// Free a spatial acceleration structure
void md_spatial_acc_free(md_spatial_acc_t* acc);

typedef struct md_spatial_acc_iter_t {
    const md_spatial_acc_t* acc;

    uint32_t cell_idx_end;
    uint32_t cell_idx;

    int32_t cc[3];
    int32_t cc_beg[3];
    int32_t cc_ext[3];
    int32_t cc_pbc[3];

    vec3_t   elem_pos;
    uint32_t elem_idx;

    vec3_t   query_pos;
    float    query_rad2;
} md_spatial_acc_iter_t;

// Query the spatial acceleration structure for all potential elements within the search space of the query
md_spatial_acc_iter_t md_spatial_acc_query(const md_spatial_acc_t* acc, vec3_t pos, float radius);

// Iterate over all potential elements within the search space of the query
bool md_spatial_acc_iter_next(md_spatial_acc_iter_t* iter);

// Get the index of the current element
static inline uint32_t md_spatial_acc_iter_idx(const md_spatial_acc_iter_t* iter) {
    return iter->elem_idx;
}

static inline vec3_t md_spatial_acc_iter_pos(const md_spatial_acc_iter_t* iter) {
    return iter->elem_pos;
}

#endif

// Initialize a spatial hash structure with a set of coordinates given by array of vec3_t (xyz) or separate arrays (x,y,z).
// pbc_ext is optional and supplies the periodic extent for each axis. To disable periodicity, supply (0,0,0).
md_spatial_hash_t* md_spatial_hash_create_vec3(const vec3_t in_xyz[], const int32_t in_idx[], size_t count, const struct md_unitcell_t* unit_cell, struct md_allocator_i* alloc);
md_spatial_hash_t* md_spatial_hash_create_soa (const float in_x[], const float in_y[], const float in_z[], const int32_t in_idx[], size_t count, const struct md_unitcell_t* unit_cell, struct md_allocator_i* alloc);

void md_spatial_hash_free(md_spatial_hash_t* spatial_hash);

// Perform a spatial query for a given position + radius in a periodic domain given by pbc_min and pbc_max.
// This will also include periodic occurrences of points accross the periodic boundries.
void md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param);

void md_spatial_hash_query_multi(const md_spatial_hash_t* spatial_hash, const vec3_t pos[], size_t count, float radius, md_spatial_hash_iterator_fn iter, void* user_param);

void md_spatial_hash_query_batch(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_batch_iter_fn iter, void* user_param);

void md_spatial_hash_query_multi_batch(const md_spatial_hash_t* spatial_hash, const vec3_t pos[], size_t count, float radius, md_spatial_hash_batch_iter_fn iter, void* user_param);

void md_spatial_hash_query_n2_batch(const md_spatial_hash_t* hash, float radius, md_spatial_hash_n2_batch_iter_fn iter, void* user_param);

// Get a list of indices which fall within the search space (pos + radius)
// Writes directly to the supplied buffer and will return the number of indices written.
size_t  md_spatial_hash_query_idx(int32_t* buf_idx, size_t buf_cap, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius);
void    md_spatial_hash_query_bits(struct md_bitfield_t* bf, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius);

#ifdef __cplusplus
}
#endif
