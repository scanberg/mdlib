#ifndef _MD_ALLOCATOR_H_
#define _MD_ALLOCATOR_H_

#include <stdint.h>

// THIS IS MASSIVELY INSPIRED BY THE API WHICH IS USED AT OURMACHINERY (https://ourmachinery.com/) 

#define md_alloc(a, sz) a->realloc(a->inst, NULL, 0, sz, __FILE__, __LINE__)
#define md_realloc(a, ptr, old_sz, new_sz) a->realloc(a->inst, ptr, old_sz, new_sz, __FILE__, __LINE__)
#define md_alloc_at(a, sz, file, line) a->realloc(a->inst, NULL, 0, sz, file, line)
#define md_free(a, p, sz) a->realloc(a->inst, p, sz, 0, __FILE__, __LINE__)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_allocator_o md_allocator_o;

typedef struct md_allocator_i {
    struct md_allocator_o* inst; // Opaque data associated with allocator

    // One function to rule them all
    // To allocate something new, then use 0 as ptr and 0 as old_size.
    // To reallocate something existing, then its quite straight forward.
    // To free something, call with a new_size of 0. ptr and old_size are still required for freeing the old data.
    void *(*realloc)(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line);
} md_allocator_i;

extern struct md_allocator_i* default_allocator;
extern struct md_allocator_i* default_temp_allocator;

#ifdef __cplusplus
}
#endif

#endif