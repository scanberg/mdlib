#ifndef _LRU_CACHE_H_
#define _LRU_CACHE_H_

#include <stdint.h>


// 8-way fully associative cache using Least Recently Used (LRU) eviction policy
typedef struct lru_cache_8 lru_cache_8;
struct lru_cache_8;

void  lru_cache_init(lru_cache_8*);
void* lru_cache_get(lru_cache_8* cache, int64_t idx);
void  lru_cache_set(lru_cache_8* cache, void* ptr);


#endif