#ifndef _LRU_CACHE_H_
#define _LRU_CACHE_H_

#include <stdint.h>
#include "intrinsics.h"

// 8-way fully associative cache using Least Recently Used (LRU) eviction policy
// Based on the '8x8 bit matrix technique' presented in Hacker's Delight
// The bit matrix is initialized to a upper triangular matrix (only 4x4 bits are shown here for simplicity)
// 0  1  1  1
// 0  0  1  1
// 0  0  0  1
// 0  0  0  0	<- LRU
//
// To set the most recently used index, we set the row of the index and then clear the column of the index:
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
// To get the least recently used index, we simply find the first zero row = byte (shown as <- LRU in the diagrams)
// Key feature that we resort to is that we can find first zero byte in a uint64 in a fast way (intrisic) which is supported on all modern CPUs
//
// The cache only stores pointers to the actual data.
// This is meant to be used as a building block for constructing a cache together with some small fixed size hash_table which holds keys that map to the cache indices.

typedef struct lru_cache_8 lru_cache_8;

struct lru_cache_8 {
	uint64_t mat;
	void*    line[8];
};

inline void lru_cache_init(lru_cache_8* cache) {
	cache->mat = 0x0080c0e0f0f8fcfe; // Upper triangular bit matrix
	memset(cache->line, 0, sizeof(void*) * 8);
}

inline int lru_cache_get_lru_idx(lru_cache_8* cache) {
	ASSERT(cache);
	ASSERT(cache->mat != 0 && "Cache matrix is 0, lru_cache is not initialized");
	return find_first_zero_byte64(cache->mat);
}

inline void lru_cache_set_mru_idx(lru_cache_8* cache, int idx) {
	ASSERT(cache);
	ASSERT(cache->mat != 0 && "Cache matrix is 0, lru_cache is not initialized");
	cache->mat = cache->mat | (0xFFLLU << (8 * idx));
	cache->mat = cache->mat & ~(0x0101010101010101LLU << idx);
}

inline void* lru_cache_get(lru_cache_8* cache, int idx) {
	ASSERT(cache);
	ASSERT(0 <= idx && idx < 8);
	lru_cache_set_mru_idx(cache, idx);
	return cache->line[idx];
}

inline void lru_cache_put(lru_cache_8* cache, void* data) {
	ASSERT(cache);
	int idx = lru_cache_get_lru_idx(cache);
	ASSERT(0 <= idx && idx <= 8);
	lru_cache_set_mru_idx(cache, idx);
	cache->line[idx] = data;
}

/*
typedef struct lru_cache_8xN lru_cache_8xN;
struct {
	uint32_t	 num_blocks;
	lru_cache_8* cache_blocks;
};
*/

#endif