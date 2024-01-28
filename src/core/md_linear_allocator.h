#pragma once

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_LINEAR_ALLOCATOR_DEFAULT_ALIGNMENT 16
#define MD_LINEAR_ALLOCATOR_MAGIC 0x8a7ccbba81980234

// This is a small, efficient linear allocator which just pushes a pointer offset within a buffer
// For simplicity, all allocations > 2 bytes are aligned to 16 bytes, which is default on x64 systems anyways
typedef struct md_linear_allocator_t {
    size_t pos;
    size_t cap;
    uint64_t magic;
    void* ptr;
} md_linear_allocator_t;

static inline void md_linear_allocator_init(md_linear_allocator_t* linear, void* backing_buffer, size_t buffer_capacity) {
    ASSERT(linear);
    ASSERT(backing_buffer);
    ASSERT(buffer_capacity > 0);

    linear->pos = 0;
    linear->cap = buffer_capacity;
    linear->magic = MD_LINEAR_ALLOCATOR_MAGIC;
    linear->ptr = backing_buffer;
}

static inline void* md_linear_allocator_push_aligned(md_linear_allocator_t* linear, size_t size, size_t align) {
    ASSERT(linear && linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    ASSERT(IS_POW2(align));

    void* mem = 0;
    size_t pos_addr = (size_t)linear->ptr + linear->pos;
    size_t alignment_size = ALIGN_TO(pos_addr, align) - pos_addr;

    if (linear->pos + alignment_size + size <= linear->cap) {
        mem = (char*)linear->ptr + linear->pos + alignment_size;
        linear->pos += alignment_size + size;
    }

    return mem;
}

static inline size_t md_linear_allocator_avail_bytes(md_linear_allocator_t* linear) {
    ASSERT(linear && linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    return (size_t)MAX(0, (int64_t)linear->cap - (int64_t)linear->pos);
}

static inline void* md_linear_allocator_push(md_linear_allocator_t* linear, size_t size) {
    size_t alignment = size <= 2 ? size : MD_LINEAR_ALLOCATOR_DEFAULT_ALIGNMENT;
    return md_linear_allocator_push_aligned(linear, size, alignment);
}

static inline void md_linear_allocator_pop(md_linear_allocator_t* linear, size_t size) {
    ASSERT(linear && linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    ASSERT(size <= linear->pos);
    linear->pos -= size;
}

static inline void md_linear_allocator_reset(md_linear_allocator_t* linear) {
    ASSERT(linear && linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    linear->pos = 0;
}

static inline void md_linear_allocator_set_pos(md_linear_allocator_t* linear, size_t pos) {
    ASSERT(linear && linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    ASSERT(pos < linear->cap);
    linear->pos = pos;
}

static inline uint64_t md_linear_allocator_get_pos(md_linear_allocator_t* linear) {
    ASSERT(linear && linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    return linear->pos;
}

md_allocator_i md_linear_allocator_create_interface(md_linear_allocator_t* linear_alloc);


#ifdef __cplusplus
}
#endif
