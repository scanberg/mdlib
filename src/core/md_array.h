#pragma once

// THIS IS INSIRED BY OUR MACHINERY'S IMPLEMENTATION OF ARRAY AT A PERVERTED LEVEL. 
// (https://ourmachinery.com)

#include <core/md_common.h>
#include <core/md_allocator.h>

#include <stdint.h>

#if MD_COMPILER_GCC
//#   pragma GCC diagnostic ignored "-Wunused-value"      // Ignore warning of unused value when returning last pushed element in md_array_push
#elif MD_COMPILER_MSVC
//#   pragma warning( disable : 6011 6387 )               // Dereferencing NULL pointer, Could be 0
#endif

typedef struct md_array_header_t {
    void*  ptr;
    size_t capacity;
    size_t size;
    size_t __unused;
} md_array_header_t;

// Silly typedef to semantically differentiate between a pointer and an array
#define md_array(type) type*

#define md_array_header(a)              ((md_array_header_t*)((uint8_t*)(a) - sizeof(md_array_header_t)))
#define md_array_size(a)                ((a) ? md_array_header(a)->size : 0)
#define md_array_bytes(a)               (md_array_size(a) * sizeof(*(a)))
#define md_array_end(a)                 ((a) ? (a + md_array_size(a)) : NULL)
#define md_array_last(a)                ((a) ? md_array_end(a) - 1 : NULL)
#define md_array_back(a)                ((a)[md_array_size(a) - 1])
#define md_array_capacity(a)            ((a) ? md_array_header(a)->capacity : 0)
#define md_array_needs_to_grow(a,n)     (md_array_capacity(a) < n)
#define md_array_pop(a)                 (--md_array_header(a)->size)
#define md_array_shrink(a,n)            ((a) ? md_array_header(a)->size = n : 0)
#define md_array_swap_back_and_pop(a,i) ((a)[i] = (a)[--md_array_header(a)->size])

// Allocator time!
#define md_array_create(type, n, alloc) ((type*)md_array_create_internal( (n), sizeof(type), alloc, __FILE__, __LINE__))
#define md_array_grow(a, n, alloc)      ((*(void **)&(a)) = md_array_capacity(a) >= (n) ? (a) : md_array_set_capacity_internal((void*)(a), md_array_grow_cap((void*)(a), (n)), sizeof(*(a)), alloc, __FILE__, __LINE__))
#define md_array_set_capacity(a, n, alloc)    ((*(void **)&(a)) = md_array_set_capacity_internal((void*)(a), (n), sizeof(*(a)), alloc, __FILE__, __LINE__))
#define md_array_resize(a, n, alloc)    ((md_array_needs_to_grow((a), (n)) ? md_array_set_capacity((a), (n), alloc) : 0), (a) ? md_array_header(a)->size = (n) : 0)
#define md_array_ensure(a, n, alloc)    (md_array_needs_to_grow((a), (n)) ? md_array_grow((a), (n), alloc) : 0)
#if MD_COMPILER_MSVC
// Suppress incorrect warnings for macros in MSVC
#define md_array_push(a, item, alloc) \
    __pragma(warning(suppress:6011 6387)) \
    (md_array_ensure((a), md_array_size(a) + 1, alloc), (a)[(md_array_header(a)->size)++] = item)

#define md_array_push_no_grow(a, item) \
    __pragma(warning(suppress:6011 6387)) \
    (a)[(md_array_header(a)->size)++] = (item)

#define md_array_push_array(a, items, n, alloc) \
    __pragma(warning(suppress:6011 6387)) \
    ((n) ? (md_array_ensure((a), md_array_size(a) + (n), alloc), MEMCPY((a) + md_array_size(a), items, (n) * sizeof(*(a))), md_array_header(a)->size += (n)) : 0)

#else
#define md_array_push(a, item, alloc) (md_array_ensure((a), md_array_size(a) + 1, alloc), (a)[md_array_header(a)->size++] = item)
#define md_array_push_no_grow(a, item)  (a)[md_array_header(a)->size++] = (item)
#define md_array_push_array(a, items, n, alloc) ((n) ? (md_array_ensure((a), md_array_size(a) + (n), alloc), MEMCPY((a) + md_array_size(a), items, (n) * sizeof(*(a))), md_array_header(a)->size += (n)) : 0)
#endif
#define md_array_free(a, alloc)         ((*(void **)&(a)) = md_array_set_capacity_internal((void *)(a), 0, sizeof(*(a)), alloc, __FILE__, __LINE__))


#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline void* md_array_set_capacity_internal(void* arr, size_t new_cap, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(alloc);
    uint8_t* p = arr ? (uint8_t*)md_array_header(arr) : 0;
    const size_t extra = sizeof(md_array_header_t);
    const size_t size = md_array_size(arr);
    const size_t old_cap = md_array_capacity(arr);
    const size_t bytes_before = arr ? item_size * old_cap + extra : 0;
    const size_t bytes_after = new_cap ? item_size * new_cap + extra : 0;

    uint8_t* new_p = (uint8_t*)alloc->realloc(alloc->inst, p, bytes_before, bytes_after, file, line);
    //void*  new_arr = new_p ? (void*)ALIGN_TO((uintptr_t)new_p + sizeof(md_array_header_t), item_alignment) : 0;
    void*  new_arr = new_p ? new_p + extra : 0;

    if (new_arr) {
        md_array_header(new_arr)->ptr = new_p;
        md_array_header(new_arr)->size = size;
        md_array_header(new_arr)->capacity = new_cap;
    }

    return new_arr;
}

static inline void* md_array_create_internal(size_t size, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    void* arr = md_array_set_capacity_internal(NULL, size, item_size, alloc, file, line);
    md_array_header(arr)->size = size;
    return arr;
}

static inline size_t md_array_grow_cap(void* arr, size_t n) {
    const size_t cap = md_array_capacity(arr);
    const size_t min_cap = cap ? cap * 2 : 8;
    return min_cap > n ? min_cap : n;
}

#ifdef __cplusplus
}
#endif
