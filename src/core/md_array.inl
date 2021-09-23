#ifndef __MD_ARRAY_INL__
#define __MD_ARRAY_INL__

// THIS IS INSIRED BY OUR MACHINERY'S IMPLEMENTATION OF MD ARRAY AT A PERVERTED LEVEL. 
// (https://ourmachinery.com)

#include "md_allocator.h"
#include "md_common.h"
#include "md_compiler.h"

#if MD_COMPILER_MSVC
#pragma warning( disable : 6387 6011 )
#endif

#include <stdint.h>
#include <string.h>

typedef struct md_array_header_t {
    int64_t capacity;
    int64_t size;
} md_array_header_t;

#define md_array_header(a)          ((md_array_header_t*)((uint8_t*)(a) - sizeof(md_array_header_t)))
#define md_array_size(a)            ((a) ? md_array_header(a)->size : 0)
#define md_array_bytes(a)           (md_array_size(a) * sizeof(*(a)))
#define md_array_end(a)             ((a) ? (a + md_array_size(a)) : NULL)
#define md_array_last(a)            ((a) ? md_array_end(a) - 1 : NULL)
#define md_array_capacity(a)        ((a) ? md_array_header(a)->capacity : 0)
#define md_array_needs_to_grow(a,n) (md_array_capacity(a) < n)
#define md_array_pop(a)             ((a)[--md_array_header(a)->size])
#define md_array_shrink(a,n)        ((a) ? md_array_header(a)->size = n : 0)

// Allocator time!
#define md_array_grow(a, n, alloc)      ((*(void **)&(a)) = md_array_grow_internal((void *)a, n, sizeof(*a), alloc, __FILE__, __LINE__))
#define md_set_capacity(a, n, alloc)    ((*(void **)&(a)) = md_array_set_capacity_internal((void *)a, n, sizeof(*a), alloc, __FILE__, __LINE__))
#define md_array_resize(a, n, alloc)    ((md_array_needs_to_grow(a, n) ? md_set_capacity(a, n, alloc) : 0), (a) ? md_array_header(a)->size = n : 0)
#define md_array_ensure(a, n, alloc)    (md_array_needs_to_grow(a, n) ? md_array_grow(a, n, alloc) : 0)
#define md_array_push(a, item, alloc)   (md_array_ensure(a, md_array_size(a) + 1, alloc), (a)[md_array_header(a)->size++] = (item), (a) + md_array_header(a)->size - 1)
#define md_array_push_array(a, items, n, alloc) ((n) ? ((md_array_ensure(a, md_array_size(a) + n, alloc), memcpy(a + md_array_size(a), items, n * sizeof(*(a))), md_array_header(a)->size += n), 0) : 0)
#define md_array_free(a, alloc)         ((*(void **)&(a)) = md_array_set_capacity_internal((void *)a, 0, sizeof(*(a)), alloc, __FILE__, __LINE__))

static inline void* md_array_set_capacity_internal(void* arr, int64_t new_cap, int64_t item_size, struct md_allocator_i* alloc, const char* file, uint32_t line) {
    ASSERT(alloc);
    uint8_t* p = arr ? (uint8_t*)md_array_header(arr) : 0;
    const int64_t extra = sizeof(md_array_header_t);
    const int64_t size = md_array_size(arr);
    const int64_t old_cap = md_array_capacity(arr);
    const int64_t bytes_before = arr ? item_size * old_cap + extra : 0;
    const int64_t bytes_after = new_cap ? item_size * new_cap + extra : 0;
    if (p && old_cap == 0) {
        // This is to deal with the fact that an array can be statically allocated from the beginning and
        // would like to grow that sucker anyways. A statically allocated array will have capacity 0
        uint8_t* old_p = p;
        p = (uint8_t*)alloc->realloc(alloc->inst, 0, 0, bytes_after, file, line);
        const int64_t static_bytes = item_size * size + extra;
        const int64_t bytes_to_copy = static_bytes > bytes_after ? bytes_after : static_bytes;
        memcpy(p, old_p, bytes_to_copy);
    } else {
        p = (uint8_t*)alloc->realloc(alloc->inst, p, bytes_before, bytes_after, file, line);
    }
    void *new_arr = p ? p + extra : p;
    if (new_arr) {
        md_array_header(new_arr)->size = size;
        md_array_header(new_arr)->capacity = new_cap;
    }
    return new_arr;
}

static inline void* md_array_grow_internal(void* arr, int64_t n, int64_t item_size, struct md_allocator_i* alloc, const char* file, uint32_t line) {
    const int64_t cap = arr ? md_array_capacity(arr) : 0;
    if (cap >= n) {
        // No need for growth
        return arr;
    }
    const int64_t min_cap = cap ? cap * 2 : 16;
    const int64_t new_cap = min_cap > n ? min_cap : n;
    return md_array_set_capacity_internal(arr, new_cap, item_size, alloc, file, line);
}

#endif
