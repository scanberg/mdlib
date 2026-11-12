#pragma once

// Sorting: spatial (Morton order) and radix sorts of unsigned integer keys.

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// The order of the points along a Morton (Z-order) curve through their bounding box: out_order[k] is the index of the
// k-th point, of count, so that points close in the order are close in space and consecutive ranges make spatially
// coherent clusters (with a jump now and then between octants, which a cluster made of a range should allow for).
// xyz points at the x of the first point, its y and z following; stride_in_bytes is the distance from one point to the
// next (sizeof(md_f32x3_t) for packed xyz, sizeof(md_f32x4_t) for xyzw), 0 for packed xyz. out_order has room for count.
void md_util_sort_spatial(uint32_t* out_order, const float* xyz, size_t stride_in_bytes, size_t count);

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
