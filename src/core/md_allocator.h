#pragma once

#include <core/md_common.h>
#include <stddef.h>

// THIS IS MASSIVELY INSPIRED BY THE API WHICH IS USED AT OURMACHINERY (https://ourmachinery.com/) 

#define md_alloc(a, sz) (a)->realloc((a)->inst, 0, 0, sz, __FILE__, __LINE__)
#define md_realloc(a, ptr, old_sz, new_sz) (a)->realloc((a)->inst, ptr, old_sz, new_sz, __FILE__, __LINE__)
#define md_alloc_at(a, sz, file, line) (a)->realloc((a)->inst, 0, 0, sz, file, line)
#define md_free(a, p, sz) (a)->realloc((a)->inst, p, sz, 0, __FILE__, __LINE__)

typedef struct md_allocator_o md_allocator_o;

typedef struct md_allocator_i {
    struct md_allocator_o* inst; // Opaque data associated with allocator

                                 // One function to rule them all
                                 // To allocate something new, then use 0 as ptr and 0 as old_size.
                                 // To reallocate something existing, then its quite straight forward.
                                 // To free something, call with a new_size of 0. ptr and old_size are still required for freeing the old data.
    void *(*realloc)(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line);
} md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

size_t md_temp_allocator_max_allocation_size(void);

// Get general allocator interface to heap (malloc)
struct md_allocator_i* md_get_heap_allocator(void);

// Get general allocator interface to thread local ring buffer
struct md_allocator_i* md_get_temp_allocator(void);

// Simple interface to the thread local ring buffer
void*  md_temp_push(size_t bytes);
void*  md_temp_push_zero(size_t bytes);
void*  md_temp_push_aligned(size_t bytes, size_t alignment);
void*  md_temp_push_zero_aligned(size_t bytes, size_t alignment);

void   md_temp_pop (size_t bytes);

size_t md_temp_get_pos(void);
void   md_temp_set_pos_back(size_t pos);

#ifdef __cplusplus
}

struct ScopedTemp {
    size_t pos;
    ScopedTemp() {
        pos = md_temp_get_pos();
    }
    ~ScopedTemp() {
        md_temp_set_pos_back(pos);
    }
};

#endif
