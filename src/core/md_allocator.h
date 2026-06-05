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

typedef struct md_temp_scope_t {
    struct md_allocator_i* arena;
    size_t pos;
} md_temp_scope_t;

typedef md_temp_scope_t md_temp_scope_t;

#ifdef __cplusplus
extern "C" {
#endif

size_t md_temp_arena_reservation_size(void);

// Legacy explicit initialization entry point. Temp arenas are also initialized lazily on first use.
void md_temp_arena_system_init(void);

// Get general allocator interface to heap (malloc)
struct md_allocator_i* md_get_heap_allocator(void);

// Begin a temporary allocation scope using the default temp arena for this thread.
md_temp_scope_t md_temp_begin(void);

// Begin a temporary allocation scope using a temp arena different from the supplied allocator.
md_temp_scope_t md_temp_begin_avoid(const md_allocator_i* avoid);

// Begin a temporary allocation scope using an explicitly supplied arena.
md_temp_scope_t md_temp_begin_in(struct md_allocator_i* arena);

// End a temporary allocation scope and rewind the arena.
void md_temp_end(md_temp_scope_t temp);

// Access the underlying allocator associated with the scope.
struct md_allocator_i* md_temp_allocator(md_temp_scope_t temp);

// Allocate from the scope's arena.
void* md_temp_alloc(md_temp_scope_t temp, size_t size);
void* md_temp_alloc_zero(md_temp_scope_t temp, size_t size);
void* md_temp_alloc_aligned(md_temp_scope_t temp, size_t size, size_t alignment);

#define md_temp_alloc_array(temp, type, count) ((type*)md_temp_alloc(temp, sizeof(type) * (count)))
#define md_temp_alloc_zero_array(temp, type, count) ((type*)md_temp_alloc_zero(temp, sizeof(type) * (count)))

#ifdef __cplusplus
}
#endif
