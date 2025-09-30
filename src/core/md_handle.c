#include <core/md_handle.h>
#include <core/md_common.h>
#include <core/md_allocator.h>

// This is a great introduction and argument for using handles instead of pointers
// https://floooh.github.io/2018/06/17/handles-vs-pointers.html
// This is massively inspired by SOKOL libraries (Fantastic set of libraries and design)
// https://github.com/floooh/sokol

void md_handle_pool_init(md_handle_pool_t* pool, int count, md_allocator_i* alloc) {
	ASSERT(pool && count >= 1);
	pool->size = count + 1;
	pool->queue_top = 0;
	size_t gen_counters_size = sizeof(uint32_t) * (size_t)pool->size;
	pool->gen_counters = (uint32_t*)md_alloc(alloc, gen_counters_size);
	MEMSET(pool->gen_counters, 0, sizeof(gen_counters_size));
	size_t free_queue_size = sizeof(int) * (size_t)count;
	pool->free_queue = (int*)md_alloc(alloc, free_queue_size);
	MEMSET(pool->free_queue, 0, free_queue_size);
	for (int i = pool->size-1; i >= 1; i--) {
		pool->free_queue[pool->queue_top++] = i;
	}
}

void md_handle_pool_free(md_handle_pool_t* pool, md_allocator_i* alloc) {
	ASSERT(pool);
	ASSERT(pool->free_queue);
	ASSERT(pool->gen_counters);
	md_free(alloc, pool->free_queue, sizeof(int) * (pool->size-1));
	md_free(alloc, pool->gen_counters, sizeof(uint32_t) * pool->size);
	MEMSET(pool, 0, sizeof(md_handle_pool_t));
}

uint32_t md_handle_pool_alloc_slot(md_handle_pool_t* pool) {
	ASSERT(pool);
	ASSERT(pool->free_queue);
	if (pool->queue_top > 0) {
		int slot_index = pool->free_queue[--pool->queue_top];
		ASSERT((slot_index > 0) && (slot_index < pool->size));
		uint32_t gen = ++pool->gen_counters[slot_index];
		return (gen<<MD_HANDLE_SLOT_SHIFT) | (slot_index & MD_HANDLE_SLOT_MASK);
	} else {
		// pool exhausted
		return MD_HANDLE_INVALID_ID;
	}
}

void md_handle_pool_free_slot(md_handle_pool_t* pool, uint32_t id) {
	ASSERT(pool);
	int slot_index = (int) (id & MD_HANDLE_SLOT_MASK);
	ASSERT((slot_index > 0) && (slot_index < pool->size));
	ASSERT(pool->free_queue);
	ASSERT(pool->queue_top < pool->size);
#ifdef DEBUG
	// debug check against double-free
	for (int i = 0; i < pool->queue_top; i++) {
		ASSERT(pool->free_queue[i] != slot_index);
	}
#endif
	pool->free_queue[pool->queue_top++] = slot_index;
	ASSERT(pool->queue_top <= (pool->size-1));
}
