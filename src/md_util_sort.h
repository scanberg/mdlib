#pragma once

// Sorting: spatial (Morton order) and radix sorts of unsigned integer keys.

#include <stdint.h>
#include <stddef.h>

#include <core/md_vec.h>

#ifdef __cplusplus
extern "C" {
#endif

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial(uint32_t* source_indices, const md_f32x3_t* xyz, size_t count);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial_xyz(uint32_t* source_indices, const float* xyz, size_t stride_in_bytes, size_t count);

// Radix sorts of unsigned integer keys, in place (sort) or as the order of the keys (argsort): out_indices[0] is the index
// of the smallest key. The argsorts are stable: indices of equal keys come in increasing order. They take temporary memory
// from the temp arena of the calling thread, for the keys and indices between passes.
void md_util_sort_radix_u32(uint32_t* data, size_t count);
void md_util_sort_radix_u64(uint64_t* data, size_t count);
void md_util_argsort_radix_u32(uint32_t* out_indices, const uint32_t* keys, size_t count);
void md_util_argsort_radix_u64(uint32_t* out_indices, const uint64_t* keys, size_t count);

#ifdef __cplusplus
}
#endif
