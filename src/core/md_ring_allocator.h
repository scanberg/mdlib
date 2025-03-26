#pragma once

#include <stddef.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// This is a small, efficient linear allocator which just pushes a pointer offset within a ring buffer and wraps around when full
// For simplicity, all allocations > 2 bytes are aligned to 16 bytes

struct md_allocator_i* md_ring_allocator_create(void* backing_buffer, size_t buffer_capacity);

void* md_ring_allocator_push_aligned(struct md_allocator_i* alloc, size_t size, size_t align);
void* md_ring_allocator_push(struct md_allocator_i* alloc, size_t size);

void md_ring_allocator_pop(struct md_allocator_i* alloc, size_t size);
void md_ring_allocator_reset(struct md_allocator_i* alloc);

void md_ring_allocator_set_pos_back(struct md_allocator_i* alloc, size_t pos);
size_t md_ring_allocator_get_pos(struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
