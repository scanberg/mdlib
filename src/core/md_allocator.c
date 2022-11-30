#include "md_allocator.h"
#include <core/md_ring_allocator.h>
#include <core/md_common.h>
#include <stdlib.h>
#include <string.h>

#ifndef MD_TEMP_ALLOC_SIZE
#define MD_TEMP_ALLOC_SIZE MEGABYTES(32)
#endif 

#define MAX_TEMP_ALLOCATION_SIZE (MD_TEMP_ALLOC_SIZE / 2)
#define DEFAULT_ALIGNMENT 16

uint64_t default_temp_allocator_max_allocation_size() {
    return MAX_TEMP_ALLOCATION_SIZE;
}

// @NOTE: Perhaps going through thread local storage for this type of default temporary allocator is bat-shit insane.
// This should be properly profiled and tested accross platforms and compilers to see what the performance implications are.

// This is a fixed buffer hacked
typedef struct ring_alloc_t {
    md_ring_allocator_t alloc;
    char buf[MD_TEMP_ALLOC_SIZE];
} ring_alloc_t;

THREAD_LOCAL static ring_alloc_t _ring = {
    .alloc = {
        .magic = MD_RING_ALLOCATOR_MAGIC,
        .pos = 0,
        .cap = MD_TEMP_ALLOC_SIZE,
        .ptr = 0,
    }
};

static void* realloc_internal(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
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

extern void* ring_realloc(struct md_allocator_o* inst, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line);

static void* ring_realloc_internal(struct md_allocator_o* inst, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)inst;
    return ring_realloc((md_allocator_o*)get_thread_ring_allocator(), ptr, old_size, new_size, file, line);
}

struct md_ring_allocator_t* get_thread_ring_allocator() {
    // This is an Ugly hack since we cannot initialize the thread local ring buffer object with the address of the other thread local buffer
    // This can probably be resolved in some fashion,
    // The constraint is that initialization fields must be constant-expressions when initializeing thread local objects.

    _ring.alloc.ptr = _ring.buf;
    return &_ring.alloc;
}

/*
static inline void* ring_alloc_internal(uint64_t size) {
    // We want to make sure that we maintain some type of alignment for allocations
    // By default on x64, the alignment should be 16 bytes when using malloc or stack allocations (alloca)
    // For single byte and double byte allocations, we could skip this default alignment, but for 4 byte and above
    // we should enforce this as x64 calling convention uses xmm registers for passing floats, which could result in a
    // performance penalty if the data is not 16 byte aligned.

    const uint64_t alignment = (size <= 2) ? size : DEFAULT_ALIGNMENT;

    ring.prev = ring.curr;

    if (ring.curr + size + alignment > sizeof(ring.mem)) {
        ring.curr = 0;
    }

    const uint64_t addr = (uint64_t)ring.mem + ring.curr;
    const uint64_t mod = addr & (alignment - 1);      // modulo: addr % alignment, but fast for power of 2
    const uint64_t curr = mod ? ring.curr + alignment - mod : ring.curr;

    ring.curr = curr + size;
    return ring.mem + curr;
}

static void* ring_realloc_internal(struct md_allocator_o* inst, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)inst;    
    (void)file;
    (void)line;
    ASSERT(new_size < MAX_TEMP_ALLOCATION_SIZE);

    if (new_size == 0) {
        // free
        // If this is the last allocation, then we move the pointer back
        if (ptr == (&ring.mem[ring.prev])) {
            ring.curr = ring.prev;
        }

        return NULL;
    }

    if (ptr && old_size) {
        // realloc
        int64_t diff = (int64_t)new_size - (int64_t)old_size;
        if (ptr == (char*)ring.mem + ring.curr) {
            // This is the last allocation that occured, then we can shrink or grow that sucker
            int64_t new_curr = (int64_t)ring.curr + diff;
            ASSERT(new_curr > 0);
            if (new_curr < MD_TEMP_ALLOC_SIZE) {
                ring.curr = new_curr;
                return ptr;
            }
        }

        // ptr is not the last allocation or the new size did not fit into the buffer which is left
        void* new_ptr = ring_alloc_internal(new_size);
        memcpy(new_ptr, ptr, old_size);
        return new_ptr;
    }
    
    // alloc
    return ring_alloc_internal(new_size);
}
*/

static struct md_allocator_i _default_allocator = {
    NULL,
    realloc_internal,
};

static struct md_allocator_i _default_temp_allocator = {
    NULL,
    ring_realloc_internal,
};

struct md_allocator_i* default_allocator = &_default_allocator;
struct md_allocator_i* default_temp_allocator = &_default_temp_allocator;
