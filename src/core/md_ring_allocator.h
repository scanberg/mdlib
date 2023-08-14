#pragma once

#include <core/md_allocator.h>
#include <core/md_common.h>

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_RING_ALLOCATOR_DEFAULT_ALIGNMENT 16
#define MD_RING_ALLOCATOR_MAGIC 0x8a7ccbba82745852

// This is a small, efficient linear allocator which just pushes a pointer offset within a buffer and wraps around when full
// For simplicity, all allocations > 2 bytes are aligned to 16 bytes
typedef struct md_ring_allocator_t {
    uint64_t pos;
    uint64_t cap;
    uint64_t magic;
    void*    ptr;
} md_ring_allocator_t;

static inline void md_ring_allocator_init(md_ring_allocator_t* ring, void* backing_buffer, uint64_t buffer_capacity) {
    ASSERT(ring);
    ASSERT(backing_buffer);
    ASSERT(buffer_capacity > 0);

    ring->pos = 0;
    ring->cap = buffer_capacity;
    ring->magic = MD_RING_ALLOCATOR_MAGIC;
    ring->ptr = backing_buffer;
}

static inline void* md_ring_allocator_push_aligned(md_ring_allocator_t* ring, uint64_t size, uint64_t align) {
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);
    ASSERT(IS_POW2(align));

    uint64_t pos_addr = (uint64_t)ring->ptr + ring->pos;
    uint64_t alignment_size = ALIGN_TO(pos_addr, align) - pos_addr;
    uint64_t pos = (ring->pos + alignment_size + size <= ring->cap) ? ring->pos + alignment_size : ALIGN_TO((uint64_t)ring->ptr, align) - (uint64_t)ring->ptr;
    ring->pos = pos + size;

    return (char*)ring->ptr + pos;
}

static inline void* md_ring_allocator_push(md_ring_allocator_t* ring, uint64_t size) {
    uint64_t alignment = size <= 2 ? size : MD_RING_ALLOCATOR_DEFAULT_ALIGNMENT;
    return md_ring_allocator_push_aligned(ring, size, alignment);
}

static inline void md_ring_allocator_pop(md_ring_allocator_t* ring, uint64_t size) {
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);
    ASSERT(size <= ring->pos);
    ring->pos -= size;
}

static inline void md_ring_allocator_reset(md_ring_allocator_t* ring) {
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);
    ring->pos = 0;
}

static inline void md_ring_allocator_set_pos(md_ring_allocator_t* ring, uint64_t offset) {
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);
    ASSERT(offset < ring->cap);
    ring->pos = offset;
}

static inline uint64_t md_ring_allocator_get_pos(md_ring_allocator_t* ring) {
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);
    return ring->pos;
}

md_allocator_i md_ring_allocator_create_interface(md_ring_allocator_t* ring_alloc);


#ifdef __cplusplus
}
#endif
