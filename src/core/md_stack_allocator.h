#pragma once

#include "core/md_allocator.h"
#include "core/md_common.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_STACK_ALLOCATOR_DEFAULT_ALIGNMENT 16
#define MD_STACK_ALLOCATOR_MAGIC 0x8a7ccbba81980234

// This is a small, efficient linear allocator which just pushes a pointer offset within a buffer
// For simplicity, all allocations > 2 bytes are aligned to 16 bytes, which is default on x64 systems anyways
typedef struct md_stack_allocator_t {
    uint64_t pos;
    uint64_t cap;
    uint64_t magic;
    void* ptr;
} md_stack_allocator_t;

static inline void md_stack_allocator_init(md_stack_allocator_t* stack, void* backing_buffer, int64_t buffer_capacity) {
    ASSERT(stack);
    ASSERT(backing_buffer);
    ASSERT(buffer_capacity > 0);

    stack->pos = 0;
    stack->cap = buffer_capacity;
    stack->magic = MD_STACK_ALLOCATOR_MAGIC;
    stack->ptr = backing_buffer;
}

static inline void* md_stack_allocator_push_aligned(md_stack_allocator_t* stack, uint64_t size, uint64_t align) {
    ASSERT(stack && stack->magic == MD_STACK_ALLOCATOR_MAGIC);
    ASSERT(IS_POW2(align));

    void* mem = 0;
    uint64_t pos_addr = (uint64_t)stack->ptr + stack->pos;
    uint64_t alignment_size = ALIGN_TO(pos_addr, align) - pos_addr;

    if (stack->pos + alignment_size <= stack->cap) {
        mem = (char*)stack->ptr + stack->pos + alignment_size;
        stack->pos += alignment_size + size;
    }

    return mem;
}

static inline void* md_stack_allocator_push(md_stack_allocator_t* stack, uint64_t size) {
    uint64_t alignment = size <= 2 ? size : MD_STACK_ALLOCATOR_DEFAULT_ALIGNMENT;
    return md_stack_allocator_push_aligned(stack, size, alignment);
}

static inline void md_stack_allocator_pop(md_stack_allocator_t* stack, uint64_t size) {
    ASSERT(stack && stack->magic == MD_STACK_ALLOCATOR_MAGIC);
    ASSERT(size <= stack->pos);
    stack->pos -= size;
}

static inline void md_stack_allocator_reset(md_stack_allocator_t* stack) {
    ASSERT(stack && stack->magic == MD_STACK_ALLOCATOR_MAGIC);
    stack->pos = 0;
}

static inline void md_stack_allocator_set_pos(md_stack_allocator_t* stack, uint64_t offset) {
    ASSERT(stack && stack->magic == MD_STACK_ALLOCATOR_MAGIC);
    ASSERT(offset < stack->cap);
    stack->pos = offset;
}

static inline uint64_t md_stack_allocator_get_pos(md_stack_allocator_t* stack) {
    ASSERT(stack && stack->magic == MD_STACK_ALLOCATOR_MAGIC);
    return stack->pos;
}

md_allocator_i md_stack_allocator_create_interface(md_stack_allocator_t* stack_alloc);


#ifdef __cplusplus
}
#endif
