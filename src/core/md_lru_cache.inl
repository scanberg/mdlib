// We do not want any include guard:
// This is only meant to be used directly inside translation units

#include <stdint.h>
#include <core/md_intrinsics.h>

// 8-way fully associative cache using Least Recently Used (LRU) eviction policy
// Based on the '8x8 bit matrix technique' presented in Hacker's Delight
// The bit matrix is initialized to a upper triangular matrix (only 4x4 bits are shown here for simplicity)
// 0  1  1  1
// 0  0  1  1
// 0  0  0  1
// 0  0  0  0	<- LRU
//
// To set the Most Recently Used index, we set the row of the index and then clear the column of the index:
// set_mru(1) ->
// 0  0  1  1
// 1  0  1  1
// 0  0  0  1
// 0  0  0  0	<- LRU
//
// set_mru(3) ->
// 0  0  1  0
// 1  0  1  0
// 0  0  0  0	<- LRU
// 1  1  1  0
//
// To get the Least Recently Used index, we simply find the first zero row = byte (shown as <- LRU in the diagrams)
// This we can find first zero byte in a uint64 building on intrinsic operations which are supported on all modern CPUs
//
// To validate the matrix, we know that at any point in time, the number of bits set should be equal
// To that of an upper triangular matrix, this count should remain constant through out the operations
// 7 + 6 + 5 + 4 + 3 + 2 + 1 = 28
//
// @TODO: Provide example of how to use this as a building-block for creating a LRU controlled cache

typedef uint64_t lru_cache_t;

// Initializes the cache matrix to an upper triangular bit matrix
#define md_lru_cache_init(cache) (cache = 0x0080c0e0f0f8fcfeLLU)

// Validates the integrity of the matrix
#define md_lru_cache_validate(cache) (popcnt64(cache) == 28LLU)

// Get the Least Recently Used index
#define md_lru_cache_get_lru(cache) (ASSERT(md_lru_cache_validate(cache)), find_first_zero_byte64(cache))

// Set the Most Recently Used index
#define md_lru_cache_set_mru(cache, idx) (ASSERT(md_lru_cache_validate(cache)), ASSERT(0 <= idx && idx < 8), cache = (cache | (0xFFLLU << (8 * idx)) & ~(0x0101010101010101LLU << idx)))
