#pragma once

#include <stdint.h>
#include <core/md_common.h>

#define MD_HANDLE_INVALID_ID (0)
#define MD_HANDLE_SLOT_SHIFT (16)
#define MD_HANDLE_SLOT_MASK ((1 << MD_HANDLE_SLOT_SHIFT) - 1)

struct md_allocator_i;

typedef struct {
	int size;
	int queue_top;
	int* free_queue;
	uint32_t* gen_counters;
} md_handle_pool_t;

void md_handle_pool_init(md_handle_pool_t* pool, int count, struct md_allocator_i* alloc);
void md_handle_pool_free(md_handle_pool_t* pool, struct md_allocator_i* alloc);
uint32_t md_handle_pool_alloc_slot(md_handle_pool_t* pool);
void	 md_handle_pool_free_slot (md_handle_pool_t* pool, uint32_t id);

static inline int md_handle_index(uint32_t id) {
	int slot_index = (int) (id & MD_HANDLE_SLOT_MASK);
	ASSERT(slot_index != 0);
	return slot_index;
}
