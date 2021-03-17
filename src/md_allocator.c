#include "md_allocator.h"
#include "core/common.h"
#include <stdlib.h>

#define THREAD_LOCAL_RING_BUFFER_SIZE MEGABYTES(1)

// @NOTE: Perhaps going through thread local storage for this type of default temporary allocator is bat-shit insane.
// This should be properly profiled and tested accross platforms and compilers to see what the performance implications are.

typedef struct ring_buffer {
    uint64_t curr;
    uint64_t prev;
    char mem[THREAD_LOCAL_RING_BUFFER_SIZE];
} ring_buffer;

THREAD_LOCAL ring_buffer ring;

internal void* realloc_internal(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)inst;
    (void)old_size;
    (void)file;
    (void)line;
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    return realloc(ptr, (size_t)new_size);
}

internal void* ring_alloc_internal(struct md_allocator_o *inst, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)inst;
    (void)ptr;
    (void)old_size;
    (void)file;
    (void)line;
    ASSERT(new_size < sizeof(ring.mem));
    
    if (ring.curr + new_size > sizeof(ring.mem)) {
        ring.curr = new_size;
        ring.prev = 0;
        return &ring.mem[0];
    }
    else {
        ring.prev = ring.curr;
        ring.curr += new_size;
        return &ring.mem[ring.prev];
    }
}

static struct md_allocator_i _default_allocator = {
    NULL,
    realloc_internal
};

static struct md_allocator_i _default_temp_allocator = {
    NULL,
    ring_alloc_internal
};

struct md_allocator_i* default_allocator = &_default_allocator;
struct md_allocator_i* default_temp_allocator = &_default_temp_allocator;
