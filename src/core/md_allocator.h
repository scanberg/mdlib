#pragma once

#include <core/md_common.h>
#include <core/md_str.h>
#include <stddef.h>
#include <stdarg.h>

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

// A point of a scope to come back to: md_temp_reset gives back what was allocated after it, and the scope stays open.
// For a phase whose scratch the next phase can use again (the spatial structure of bond inference, then the bonds),
// without a scope of its own around it. Only the innermost scope (on its thread) is reset, to a mark of its own.
typedef struct md_temp_mark_t {
    size_t pos;
} md_temp_mark_t;

md_temp_mark_t md_temp_mark(md_temp_scope_t temp);
void md_temp_reset(md_temp_scope_t temp, md_temp_mark_t mark);

// Gives back to the OS what the calling thread's temp arenas hold committed above what is in use. An arena
// keeps what it once grew to (ending a scope frees nothing), so a thread that ran something large once - a
// pool worker given a protein's orbitals - holds that much until it exits. For such a thread after such work.
void md_temp_trim(void);

// Access the underlying allocator associated with the scope.
struct md_allocator_i* md_temp_allocator(md_temp_scope_t temp);

// Allocate from the scope's arena.
void* md_temp_alloc(md_temp_scope_t temp, size_t size);
void* md_temp_alloc_zero(md_temp_scope_t temp, size_t size);
void* md_temp_alloc_aligned(md_temp_scope_t temp, size_t size, size_t alignment);

#define md_temp_alloc_array(temp, type, count) ((type*)md_temp_alloc(temp, sizeof(type) * (count)))
#define md_temp_alloc_zero_array(temp, type, count) ((type*)md_temp_alloc_zero(temp, sizeof(type) * (count)))

// ### STRINGS ###
// md_str_t is a view and never owns what it points at: these copy into memory from alloc, which owns the
// result. There is no freeing one string on its own; pass an arena (or a temp scope's allocator) that lives as
// long as the strings do, and free that as a whole.

// A zero terminated copy of str (the terminator is not counted in len). An empty str gives an empty result and
// allocates nothing.
md_str_t md_alloc_str(md_allocator_i* alloc, md_str_t str);

// Formats with md_vsnprintf (the "C" locale) into memory from alloc. Empty, and nothing allocated, when nothing
// is written.
md_str_t md_alloc_printf (md_allocator_i* alloc, const char* format, ...) MD_PRINTF_FORMAT(2, 3);

// TRANSITIONAL: frees a string from the functions above, made with an allocator that frees (the heap). It is
// only here for the owners which still free their strings one by one; it goes once they have an arena.
void md_free_str(md_allocator_i* alloc, md_str_t str);

#ifdef __cplusplus
}
#endif
