#pragma once

#include <core/md_intrinsics.h>
#include <stdint.h>

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

typedef uint64_t md_lru_cache8_t;

// 8-way LRU helpers (static inline -> type-checked, no multiple-eval, easy to test)
static inline void md_lru_cache8_init(md_lru_cache8_t* cache) {
	*cache = 0x0080C0E0F0F8FCFELLU;
}

static inline int md_lru_cache8_validate(md_lru_cache8_t cache) {
	return popcnt64(cache) == 28LLU;
}

static inline int md_lru_cache8_get_lru(md_lru_cache8_t cache) {
	ASSERT(md_lru_cache8_validate(cache));
	return (int)find_first_zero_byte64(cache);
}

// Set `idx` as MRU
static inline void md_lru_cache8_set_mru(md_lru_cache8_t* cache, int idx) {
	ASSERT(md_lru_cache8_validate(*cache));
	ASSERT(0 <= idx && idx < 8);
	*cache |= 0xFFULL << (8 * idx);
	*cache &= ~(0x0101010101010101ULL << idx);
}

// -------------------------------------------------------------------------
// 4-way (4x4) LRU cache helpers
// Uses the same technique but stored in 16 bits (4x4 bit matrix)
typedef uint16_t md_lru_cache4_t;

static inline md_lru_cache4_t md_lru_cache4_init(md_lru_cache4_t* cache) {
    // 0  1  1  1
    // 0  0  1  1
    // 0  0  0  1
    // 0  0  0  0	<- LRU
	return *cache = 0x8CEU;
}

static inline int md_lru_cache4_validate(md_lru_cache4_t cache) {
	// number of set bits in upper-triangular 4x4 = 4*3/2 = 6
	return popcnt32(cache) == 6LLU;
}

static inline int md_lru_cache4_get_lru(md_lru_cache4_t cache) {
	ASSERT(md_lru_cache4_validate(cache));
	return (int)find_first_zero_nibble16(cache);
}

static inline void md_lru_cache4_set_mru(md_lru_cache4_t* cache, int idx) {
	ASSERT(md_lru_cache4_validate(*cache));
	ASSERT(0 <= idx && idx < 4);
	// Promote idx to a byte-row mask then clear the column
	*cache |= 0xFU << (4 * idx); // Set row
    *cache &= ~(0x1111U << idx); // Clear column
}
