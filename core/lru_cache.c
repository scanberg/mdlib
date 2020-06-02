#include "lru_cache.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifndef ASSERT
#define ASSERT assert
#endif

typedef struct node_t node_t;
struct node_t {
	int next;
	int prev;
};

typedef struct internal_lru_cache_t internal_lru_cache_t;
struct internal_lru_cache_t {
	int capacity;
	int size;
	int item_size;
	int node_front_idx;
	int node_back_idx;

	lru_cache_key_t* keys;
	void* 			 items;
	node_t* 		 nodes;	
};

static const int cache_size = sizeof(lru_cache_t);

void lru_cache_init(lru_cache_t* ext_cache, int capacity, int item_size) {
	internal_lru_cache_t* cache = (internal_lru_cache_t*)ext_cache;

	ASSERT(cache);
	ASSERT(capacity > 0);
	ASSERT(item_size > 0);

	int mem_size = capacity * (sizeof(lru_cache_key_t) + item_size + sizeof(node_t));
	void* mem = calloc(mem_size, 1);
	ASSERT(mem);

	cache->capacity = capacity;
	cache->item_size = item_size;
	cache->keys  	  = (char*)mem;
	cache->items 	  = (char*)mem + capacity * sizeof(lru_cache_key_t);
	cache->nodes 	  = (char*)mem + capacity * (item_size + sizeof(lru_cache_key_t));

	lru_cache_clear(ext_cache);
}

void lru_cache_free(lru_cache_t* ext_cache) {
	ASSERT(ext_cache);
	internal_lru_cache_t* cache = (internal_lru_cache_t*)ext_cache;
	free(cache->keys);
}

void lru_cache_clear(lru_cache_t* ext_cache) {
	ASSERT(ext_cache);
	internal_lru_cache_t* cache = (internal_lru_cache_t*)ext_cache;
	cache->size = 0;

	for (int i = 0; i < cache->capacity; ++i) {
		cache->nodes[i].prev = (i == 0) ? 0 : i-1;
		cache->nodes[i].next = (i == cache->capacity - 1) ? 0 : i+1;
	}
	cache->node_front_idx = 0;
	cache->node_back_idx = cache->capacity - 1;
}

static inline int lru_cache_find(internal_lru_cache_t* cache, lru_cache_key_t key) {
	for (int i = 0; i < cache->size; ++i) {
		if (cache->keys[i] == key) return i;
	}
	return -1;
}

static inline int lru_cache_pop_back(internal_lru_cache_t* cache) {
	int idx = cache->node_back_idx;
	cache->node_back_idx = cache->nodes[idx].prev;
	return idx;
}

static inline void lru_cache_pop(internal_lru_cache_t* cache, int idx) {
	int prev = cache->nodes[idx].prev;
	int next = cache->nodes[idx].next;
	cache->nodes[prev].next = next;
	cache->nodes[next].prev = prev;
}

static inline void lru_cache_push_front(internal_lru_cache_t* cache, int idx) {
	int old = cache->node_front_idx;
	cache->nodes[old].prev = idx;
	cache->nodes[idx].next = old;
	cache->node_front_idx = idx;
}

void lru_cache_put(lru_cache_t* ext_cache, lru_cache_key_t key, const void* item) {
	internal_lru_cache_t* cache = (internal_lru_cache_t*)ext_cache;
	ASSERT(lru_cache_find(cache, key) == -1);

	int idx = lru_cache_pop_back(cache);
	lru_cache_push_front(cache, idx);
	cache->keys[idx] = key;
	memcpy((char*)cache->items + cache->item_size * idx, item, cache->item_size);
	cache->size = cache->size + 1 < cache->capacity ? cache->size + 1 : cache->capacity;
}

void* lru_cache_get(lru_cache_t* ext_cache, lru_cache_key_t key) {
	internal_lru_cache_t* cache = (internal_lru_cache_t*)ext_cache;
	int idx = lru_cache_find(cache, key);
	if (idx != -1) {
		lru_cache_pop(cache, idx);
		lru_cache_push_front(cache, idx);
		return (char*)cache->items + cache->item_size * idx;
	}
	return NULL;
}

#if 1
int main() {
	lru_cache_t cache = {0};

	lru_cache_init(&cache, 8, sizeof(int));

	for (uint32_t i = 0; i < 8; ++i) {
		lru_cache_put(&cache, i, &i);
	}

	for (uint32_t i = 0; i < 10000; ++i) {
		lru_cache_get(&cache, i % 8);
	}

	lru_cache_free(&cache);

	return 0;
}
#endif