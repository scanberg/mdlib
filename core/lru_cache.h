#ifndef _LRU_CACHE_H_
#define _LRU_CACHE_H_

#include <stdint.h>

typedef uint64_t lru_cache_key_t;

struct lru_cache_t {
    uint64_t _mem[6];
};
typedef struct lru_cache_t lru_cache_t;

void lru_cache_init(lru_cache_t* cache, int item_size, int capacity);
void lru_cache_free(lru_cache_t* cache);

void lru_cache_clear(lru_cache_t* cache);

void  lru_cache_put(lru_cache_t* cache, lru_cache_key_t key, const void* item);
void* lru_cache_get(lru_cache_t* cache, lru_cache_key_t key);

#endif