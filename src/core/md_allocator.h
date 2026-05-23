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

                                 // One procedure to rule them all
                                 // To allocate something new, then use 0 as ptr and 0 as old_size.
                                 // To reallocate something existing, then its quite straight forward.
                                 // To free something, call with a new_size of 0. ptr and old_size are still required for freeing the old data.
    void *(*realloc)(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line);
} md_allocator_i;

typedef struct md_temp_t {
    struct md_allocator_i* arena;
    size_t pos;
    struct md_allocator_i* prev_arena;
    size_t depth;
} md_temp_t;

#ifdef __cplusplus
extern "C" {
#endif

size_t md_temp_arena_reservation_size(void);

// Get general allocator interface to heap (malloc)
struct md_allocator_i* md_get_heap_allocator(void);

// Get current thread local temporary allocator
struct md_allocator_i* md_get_temp_arena(void);

// Get a thread local temporary allocator which does not alias any supplied allocator
struct md_allocator_i* md_get_temp_arena_avoid(struct md_allocator_i* const* conflicts, size_t conflict_count);

md_temp_t md_temp_begin(void);
md_temp_t md_temp_begin_avoid(md_allocator_i* const* conflicts, size_t conflict_count);
md_temp_t md_temp_begin_arena(md_allocator_i* arena);
void md_temp_end(md_temp_t temp);

md_allocator_i* md_temp_allocator(md_temp_t temp);

void* md_temp_push(size_t size);
void* md_temp_push_zero(size_t size);
void* md_temp_push_aligned(size_t size, size_t alignment);

#define md_temp_push_array(type, count) ((type*)md_temp_push(sizeof(type) * (count)))
#define md_temp_push_zero_array(type, count) ((type*)md_temp_push_zero(sizeof(type) * (count)))

#ifdef __cplusplus
}

struct ScopedTemp {
    md_temp_t temp;

    ScopedTemp() {
        temp = md_temp_begin();
    }

    explicit ScopedTemp(md_allocator_i* conflict) {
        temp = md_temp_begin_avoid(&conflict, 1);
    }

    ScopedTemp(md_allocator_i* const* conflicts, size_t conflict_count) {
        temp = md_temp_begin_avoid(conflicts, conflict_count);
    }

    ~ScopedTemp() {
        md_temp_end(temp);
    }

    md_allocator_i* allocator() const {
        return md_temp_allocator(temp);
    }
};

#endif
