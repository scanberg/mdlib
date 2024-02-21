#pragma once

#include <stddef.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// This is a small, efficient linear allocator which just pushes a pointer offset within a buffer
// For simplicity, all allocations > 2 bytes are aligned to 16 bytes, which is default on x64 systems anyways
// The buffer is fixed in size and cannot grow beyond its capacity
// If there is a need for a growable linear allocator 

// This creates an allocator interface to the linear allocator in the beginning of the buffer
// No need to explicitly destroy it since it is contained within the supplied buffer
struct md_allocator_i* md_linear_allocator_create(void* backing_buffer, size_t buffer_capacity);

size_t md_linear_allocator_avail_bytes(struct md_allocator_i* alloc);

void* md_linear_allocator_push(struct md_allocator_i* alloc, size_t size);
void* md_linear_allocator_push_aligned(struct md_allocator_i* alloc, size_t size, size_t align);

void  md_linear_allocator_pop(struct md_allocator_i* alloc, size_t size);
void  md_linear_allocator_reset(struct md_allocator_i* alloc);

// Modify the position within the buffer
void   md_linear_allocator_set_pos_back(struct md_allocator_i* alloc, size_t pos);
size_t md_linear_allocator_get_pos(struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
