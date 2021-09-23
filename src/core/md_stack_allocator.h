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
    char* buf;
    uint64_t capacity;
    uint64_t curr;
    uint64_t prev;
    uint64_t magic;
} md_stack_allocator_t;

static inline void md_stack_allocator_init(md_stack_allocator_t* stack_alloc, void* backing_buffer, int64_t buffer_capacity) {
    ASSERT(stack_alloc);
    ASSERT(backing_buffer);
    ASSERT(buffer_capacity > 0);

    stack_alloc->buf = (char*)backing_buffer;
    stack_alloc->capacity = buffer_capacity;
    stack_alloc->curr = 0;
    stack_alloc->prev = 0;
    stack_alloc->magic = MD_STACK_ALLOCATOR_MAGIC;
}

static inline void* md_stack_allocator_alloc(md_stack_allocator_t* stack_alloc, int64_t size) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);

    const int64_t aligned_size = (size + (MD_STACK_ALLOCATOR_ALIGNMENT-1)) & -MD_STACK_ALLOCATOR_ALIGNMENT;
    ASSERT(stack_alloc->curr + aligned_size < stack_alloc->capacity);

    stack_alloc->prev = stack_alloc->curr;
    stack_alloc->curr += aligned_size;
    return stack_alloc->buf + stack_alloc->prev;
}

static inline void md_stack_allocator_free(md_stack_allocator_t* stack_alloc, void* mem, uint64_t size) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    (void)size;

    if (mem == stack_alloc->buf + stack_alloc->prev) {
        // compare against aligned size
        ASSERT(((size + (MD_STACK_ALLOCATOR_ALIGNMENT-1)) & -MD_STACK_ALLOCATOR_ALIGNMENT) == stack_alloc->curr - stack_alloc->prev);
        stack_alloc->curr = stack_alloc->prev;
    }
}

static inline void md_stack_allocator_reset(md_stack_allocator_t* stack_alloc) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    stack_alloc->curr = 0;
    stack_alloc->prev = 0;
}

static inline void md_stack_allocator_set_offset(md_stack_allocator_t* stack_alloc, uint64_t offset) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    ASSERT(0 <= offset && offset < stack_alloc->capacity);
    stack_alloc->curr = offset;
    stack_alloc->prev = offset;
}

static inline uint64_t md_stack_allocator_get_offset(md_stack_allocator_t* stack_alloc) {
    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);
    return stack_alloc->curr;
}

md_allocator_i md_stack_allocator_create_interface(md_stack_allocator_t* stack_alloc);


#ifdef __cplusplus
}
#endif
