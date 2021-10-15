#pragma once

#include <stdint.h>
#include "core/md_allocator.h"
#include "core/md_common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MD_STACK_ALLOCATOR_ALIGNMENT 16
#define MD_STACK_ALLOCATOR_MAGIC 0x8a7ccbba81980234

// This is a small, efficient linear allocator which just pushes a pointer offset within a buffer
// For simplicity, all allocations are aligned to 16 bytes
typedef struct md_stack_allocator_t {
    uint8_t* buf;
    uint64_t cap;
    uint64_t cur;
    uint64_t prv;
    uint64_t magic;
} md_stack_allocator_t;

static inline void md_stack_allocator_init(md_stack_allocator_t* stack_alloc, void* backing_buffer, int64_t buffer_capacity) {
    ASSERT(stack_alloc);
    ASSERT(backing_buffer);
    ASSERT(buffer_capacity > 0);

    stack_alloc->buf = (uint8_t*)backing_buffer;
    stack_alloc->cap = buffer_capacity;
    stack_alloc->cur = 0;
    stack_alloc->prv = 0;
    stack_alloc->magic = MD_STACK_ALLOCATOR_MAGIC;
}

static inline void* md_stack_allocator_alloc(md_stack_allocator_t* stack_alloc, int64_t size) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);

    const int64_t aligned_size = (size + (MD_STACK_ALLOCATOR_ALIGNMENT-1)) & -MD_STACK_ALLOCATOR_ALIGNMENT;
    ASSERT(stack_alloc->cur + aligned_size < stack_alloc->cap);

    stack_alloc->prv = stack_alloc->cur;
    stack_alloc->cur += aligned_size;
    return stack_alloc->buf + stack_alloc->prv;
}

static inline void md_stack_allocator_free(md_stack_allocator_t* stack_alloc, void* mem, uint64_t size) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    (void)size;

    if (mem == stack_alloc->buf + stack_alloc->prv) {
        // compare against aligned size
        ASSERT(((size + (MD_STACK_ALLOCATOR_ALIGNMENT-1)) & -MD_STACK_ALLOCATOR_ALIGNMENT) == stack_alloc->cur - stack_alloc->prv);
        stack_alloc->cur = stack_alloc->prv;
    }
}

static inline void md_stack_allocator_reset(md_stack_allocator_t* stack_alloc) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    stack_alloc->cur = 0;
    stack_alloc->prv = 0;
}

static inline void md_stack_allocator_set_offset(md_stack_allocator_t* stack_alloc, uint64_t offset) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    ASSERT(0 <= offset && offset < stack_alloc->cap);
    stack_alloc->cur = offset;
    stack_alloc->prv = offset;
}

static inline uint64_t md_stack_allocator_get_offset(md_stack_allocator_t* stack_alloc) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    return stack_alloc->cur;
}

md_allocator_i md_stack_allocator_create_interface(md_stack_allocator_t* stack_alloc);


#ifdef __cplusplus
}
#endif
