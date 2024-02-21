#pragma once

#include <stddef.h>

// THIS IS MASSIVELY INSPIRED BY THE API WHICH IS USED AT OURMACHINERY (https://ourmachinery.com/) 

#define md_alloc(a, sz) (a)->realloc((a)->inst, 0, 0, sz, __FILE__, __LINE__)
#define md_realloc(a, ptr, old_sz, new_sz) (a)->realloc((a)->inst, ptr, old_sz, new_sz, __FILE__, __LINE__)
#define md_alloc_at(a, sz, file, line) (a)->realloc((a)->inst, 0, 0, sz, file, line)
#define md_free(a, p, sz) (a)->realloc((a)->inst, p, sz, 0, __FILE__, __LINE__)

struct md_ring_allocator_t;
typedef struct md_allocator_o md_allocator_o;

typedef struct md_allocator_i {
    struct md_allocator_o* inst; // Opaque data associated with allocator

                                 // One function to rule them all
                                 // To allocate something new, then use 0 as ptr and 0 as old_size.
                                 // To reallocate something existing, then its quite straight forward.
                                 // To free something, call with a new_size of 0. ptr and old_size are still required for freeing the old data.
    void *(*realloc)(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line);
} md_allocator_i;

static inline void* md_aligned_alloc(struct md_allocator_i* alloc, size_t size, size_t align) {
    size_t offset = align - 1 + sizeof(void*) + sizeof(size_t);
    void* p1 = md_alloc(alloc, size + offset);
    if (!p1) return 0;
    void** p2 = (void**)(((size_t)(p1) + offset) & ~(align - 1));
    p2[-1] = p1;
    p2[-2] = (void*)(size + offset); // store real allocation size

    return p2;
}

static inline void md_aligned_free(struct md_allocator_i* alloc, void* ptr, size_t size) {
    (void)size;
    void** p2 = (void**)ptr;
    void* p1 = p2[-1];
    size_t alloc_size = (size_t)p2[-2];
    md_free(alloc, p1, alloc_size);
}

#ifdef __cplusplus
extern "C" {
#endif

size_t md_temp_allocator_max_allocation_size(void);

// General allocator interface to heap (malloc)
extern struct md_allocator_i* md_heap_allocator;

// General allocator interface to thread local ring buffer
extern struct md_allocator_i* md_temp_allocator;

// Simple interface to the thread local ring buffer
void*  md_temp_push(size_t bytes);
void   md_temp_pop (size_t bytes);

size_t md_temp_get_pos();
void   md_temp_set_pos_back(size_t pos);

#ifdef __cplusplus
}
#endif



