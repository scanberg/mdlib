#pragma once

// THIS IS INSIRED BY OUR MACHINERY'S IMPLEMENTATION OF ARRAY AT A PERVERTED LEVEL. 
// (https://ourmachinery.com)

#include <core/md_common.h>
#include <core/md_allocator.h>

#include <stdint.h>
#include <stdbool.h>

#if MD_COMPILER_GCC
//#   pragma GCC diagnostic ignored "-Wunused-value"      // Ignore warning of unused value when returning last pushed element in md_array_push
#elif MD_COMPILER_MSVC
//#   pragma warning( disable : 6011 6387 )               // Dereferencing NULL pointer, Could be 0
#endif

typedef struct md_array_header_t {
    size_t capacity;
    size_t size;
} md_array_header_t;

// Silly typedef to semantically differentiate between a pointer and an array
#define md_array(type) type*

#define md_array_header(a)              ((md_array_header_t*)((uint8_t*)(a) - sizeof(md_array_header_t)))
#define md_array_size(a)                ((a) ? md_array_header(a)->size : 0)
#define md_array_bytes(a)               (md_array_size(a) * sizeof(*(a)))
#define md_array_beg(a)                 (a)
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
#define md_array_extend(a, n, alloc)    ((n) ? (md_array_ensure((a), md_array_size(a) + (n), alloc), md_array_header(a)->size += (n), &((a)[md_array_header(a)->size - (n)])) : 0)

#if MD_COMPILER_MSVC
// Suppress incorrect warnings for macros in MSVC
#define md_array_push(a, item, alloc) \
    __pragma(warning(suppress:6011 6387)) \
    (md_array_ensure((a), md_array_size(a) + 1, alloc), (a)[(md_array_header(a)->size)++] = (item))

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
#define md_array_free(a, alloc)         ((a) ? (*(void **)&(a)) = md_array_set_capacity_internal((void *)(a), 0, sizeof(*(a)), alloc, __FILE__, __LINE__) : 0)

// -----------------------------------------------------------------------------
// ORTHOGONAL TYPED RESIZABLE ARRAY API (RAW + STRUCT WRAPPERS)
// -----------------------------------------------------------------------------
// This API is independent from the hidden-header pointer arrays above.
// It is intended for explicit, user-defined array structs and gradual migration.
//
// RAW tier:
//   Operates directly on lvalues at the callsite:
//     capacity, size, data, alloc
//   Example lvalues: s.capacity, s.size, s.data, s.alloc
//
// STRUCT tier:
//   Convenience wrappers over RAW for structs with fields:
//     capacity, size, data, alloc
//
// Notes:
// - The RAW operations are strongly typed through 'data' (T*).
// - Use TRY variants when allocation failure should be handled gracefully.
// - Non-TRY variants ASSERT on allocation failure.

// Define a standalone typed array struct.
// Example:
//   MD_ARRAY_RAW_DEFINE_TYPE(int_array_t, int);
//   int_array_t arr;
//   MD_ARRAY_INIT(arr, alloc);
//   MD_ARRAY_PUSH(arr, 42);
//   MD_ARRAY_FREE(arr);
#define MD_ARRAY_RAW_DEFINE_TYPE(name, type) \
    typedef struct name {                    \
        size_t capacity;                     \
        size_t size;                         \
        type* data;                          \
        struct md_allocator_i* alloc;        \
    } name

// Convenience form for inline declarations:
//   MD_ARRAY_RAW_TYPE(int) tmp = {0};
#define MD_ARRAY_RAW_TYPE(type) struct { size_t capacity; size_t size; type* data; struct md_allocator_i* alloc; }

// RAW lifecycle & queries
#define MD_ARRAY_RAW_INIT(capacity, size, data) \
    do {                                        \
        (capacity) = 0;                         \
        (size) = 0;                             \
        (data) = NULL;                          \
    } while (0)

#define MD_ARRAY_RAW_INIT_WITH_ALLOC(capacity, size, data, alloc, allocator) \
    do {                                                                       \
        MD_ARRAY_RAW_INIT((capacity), (size), (data));                         \
        (alloc) = (allocator);                                                 \
    } while (0)

#define MD_ARRAY_RAW_SIZE(size)      (size)
#define MD_ARRAY_RAW_CAPACITY(cap)   (cap)
#define MD_ARRAY_RAW_DATA(data)      (data)
#define MD_ARRAY_RAW_CLEAR(size)     ((size) = 0)

// RAW allocation primitives
#define MD_ARRAY_RAW_TRY_SET_CAPACITY(capacity, size, data, new_capacity, alloc) \
    md_array_raw_set_capacity_impl(&(data), &(capacity), &(size), (size_t)(new_capacity), sizeof(*(data)), (alloc), __FILE__, __LINE__)

#define MD_ARRAY_RAW_SET_CAPACITY(capacity, size, data, new_capacity, alloc) \
    do {                                                                       \
        ASSERT(MD_ARRAY_RAW_TRY_SET_CAPACITY((capacity), (size), (data), (new_capacity), (alloc))); \
    } while (0)

#define MD_ARRAY_RAW_TRY_ENSURE(capacity, size, data, min_capacity, alloc) \
    md_array_raw_ensure_impl(&(data), &(capacity), (size), (size_t)(min_capacity), sizeof(*(data)), (alloc), __FILE__, __LINE__)

#define MD_ARRAY_RAW_ENSURE(capacity, size, data, min_capacity, alloc) \
    do {                                                                 \
        ASSERT(MD_ARRAY_RAW_TRY_ENSURE((capacity), (size), (data), (min_capacity), (alloc))); \
    } while (0)

// RAW element operations
#define MD_ARRAY_RAW_TRY_PUSH(capacity, size, data, item, alloc) \
    md_array_raw_push_impl(&(data), &(capacity), &(size), sizeof(*(data)), &(item), (alloc), __FILE__, __LINE__)

#define MD_ARRAY_RAW_PUSH(capacity, size, data, item, alloc) \
    do {                                                       \
        ASSERT(MD_ARRAY_RAW_TRY_PUSH((capacity), (size), (data), (item), (alloc))); \
    } while (0)

#define MD_ARRAY_RAW_POP(size) (--(size))

#define MD_ARRAY_RAW_BACK(data, size) ((data)[(size) - 1])

#define MD_ARRAY_RAW_TRY_EXTEND(capacity, size, data, count, alloc) \
    md_array_raw_extend_impl(&(data), &(capacity), &(size), sizeof(*(data)), (size_t)(count), (alloc), __FILE__, __LINE__)

#define MD_ARRAY_RAW_EXTEND(capacity, size, data, count, alloc) \
    md_array_raw_extend_assert_impl(&(data), &(capacity), &(size), sizeof(*(data)), (size_t)(count), (alloc), __FILE__, __LINE__)

#define MD_ARRAY_RAW_PUSH_ARRAY(capacity, size, data, items, count, alloc) \
    do {                                                                     \
        void* _md_dst = MD_ARRAY_RAW_EXTEND((capacity), (size), (data), (count), (alloc)); \
        if (_md_dst) MEMCPY(_md_dst, (items), (count) * sizeof(*(data)));    \
    } while (0)

#define MD_ARRAY_RAW_SWAP_BACK_AND_POP(data, size, idx) \
    do {                                                 \
        (data)[(idx)] = (data)[--(size)];               \
    } while (0)

#define MD_ARRAY_RAW_TRY_FREE(capacity, size, data, alloc) \
    md_array_raw_set_capacity_impl(&(data), &(capacity), &(size), 0, sizeof(*(data)), (alloc), __FILE__, __LINE__)

#define MD_ARRAY_RAW_FREE(capacity, size, data, alloc) \
    do {                                                \
        ASSERT(MD_ARRAY_RAW_TRY_FREE((capacity), (size), (data), (alloc))); \
    } while (0)

// STRUCT convenience wrappers (expects fields: capacity, size, data, alloc)
#define MD_ARRAY_INIT(arr, allocator) \
    MD_ARRAY_RAW_INIT_WITH_ALLOC((arr).capacity, (arr).size, (arr).data, (arr).alloc, (allocator))

#define MD_ARRAY_SIZE(arr)      MD_ARRAY_RAW_SIZE((arr).size)
#define MD_ARRAY_CAPACITY(arr)  MD_ARRAY_RAW_CAPACITY((arr).capacity)
#define MD_ARRAY_DATA(arr)      MD_ARRAY_RAW_DATA((arr).data)
#define MD_ARRAY_CLEAR(arr)     MD_ARRAY_RAW_CLEAR((arr).size)

#define MD_ARRAY_TRY_SET_CAPACITY(arr, new_capacity) \
    MD_ARRAY_RAW_TRY_SET_CAPACITY((arr).capacity, (arr).size, (arr).data, (new_capacity), (arr).alloc)

#define MD_ARRAY_SET_CAPACITY(arr, new_capacity) \
    MD_ARRAY_RAW_SET_CAPACITY((arr).capacity, (arr).size, (arr).data, (new_capacity), (arr).alloc)

#define MD_ARRAY_TRY_ENSURE(arr, min_capacity) \
    MD_ARRAY_RAW_TRY_ENSURE((arr).capacity, (arr).size, (arr).data, (min_capacity), (arr).alloc)

#define MD_ARRAY_ENSURE(arr, min_capacity) \
    MD_ARRAY_RAW_ENSURE((arr).capacity, (arr).size, (arr).data, (min_capacity), (arr).alloc)

#define MD_ARRAY_TRY_PUSH(arr, item) \
    MD_ARRAY_RAW_TRY_PUSH((arr).capacity, (arr).size, (arr).data, (item), (arr).alloc)

#define MD_ARRAY_PUSH(arr, item) \
    MD_ARRAY_RAW_PUSH((arr).capacity, (arr).size, (arr).data, (item), (arr).alloc)

#define MD_ARRAY_POP(arr) \
    MD_ARRAY_RAW_POP((arr).size)

#define MD_ARRAY_BACK(arr) \
    MD_ARRAY_RAW_BACK((arr).data, (arr).size)

#define MD_ARRAY_TRY_EXTEND(arr, count) \
    MD_ARRAY_RAW_TRY_EXTEND((arr).capacity, (arr).size, (arr).data, (count), (arr).alloc)

#define MD_ARRAY_EXTEND(arr, count) \
    MD_ARRAY_RAW_EXTEND((arr).capacity, (arr).size, (arr).data, (count), (arr).alloc)

#define MD_ARRAY_PUSH_ARRAY(arr, items, count) \
    MD_ARRAY_RAW_PUSH_ARRAY((arr).capacity, (arr).size, (arr).data, (items), (count), (arr).alloc)

#define MD_ARRAY_SWAP_BACK_AND_POP(arr, idx) \
    MD_ARRAY_RAW_SWAP_BACK_AND_POP((arr).data, (arr).size, (idx))

#define MD_ARRAY_TRY_FREE(arr) \
    MD_ARRAY_RAW_TRY_FREE((arr).capacity, (arr).size, (arr).data, (arr).alloc)

#define MD_ARRAY_FREE(arr) \
    MD_ARRAY_RAW_FREE((arr).capacity, (arr).size, (arr).data, (arr).alloc)

// -----------------------------------------------------------------------------
// Usage examples (for quick copy/paste)
// -----------------------------------------------------------------------------
// 1) Define a typed array and use struct wrappers:
//
//   MD_ARRAY_RAW_DEFINE_TYPE(float_array_t, float);
//   float_array_t arr;
//   MD_ARRAY_INIT(arr, alloc);
//   MD_ARRAY_PUSH(arr, 1.0f);
//   MD_ARRAY_PUSH(arr, 2.0f);
//   float* dst = MD_ARRAY_EXTEND(arr, 3);
//   if (dst) {
//       dst[0] = 3.0f;
//       dst[1] = 4.0f;
//       dst[2] = 5.0f;
//   }
//   MD_ARRAY_FREE(arr);
//
// 2) Use RAW ops directly on custom field paths:
//
//   MD_ARRAY_RAW_TRY_PUSH(state.buf.capacity,
//                         state.buf.size,
//                         state.buf.data,
//                         value,
//                         shared_alloc);
//
// 3) If allocation must be handled, use TRY variants:
//
//   if (!MD_ARRAY_TRY_PUSH(arr, value)) {
//       // Handle out-of-memory.
//   }

#ifdef __cplusplus
extern "C" {
#endif

static inline size_t md_array_raw_grow_cap_internal(size_t cap, size_t n) {
    const size_t min_cap = cap ? cap * 2 : 8;
    return min_cap > n ? min_cap : n;
}

static inline bool md_array_raw_set_capacity_internal(void** data_ptr, size_t* cap_ptr, size_t size, size_t new_cap, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(data_ptr);
    ASSERT(cap_ptr);
    ASSERT(alloc);
    ASSERT(size <= *cap_ptr);

    const size_t old_cap = *cap_ptr;
    size_t bytes_before = 0;
    size_t bytes_after = 0;

    if (old_cap && item_size > SIZE_MAX / old_cap) {
        return false;
    }
    bytes_before = old_cap * item_size;

    if (new_cap && item_size > SIZE_MAX / new_cap) {
        return false;
    }
    bytes_after = new_cap * item_size;

    void* new_data = alloc->realloc(alloc->inst, *data_ptr, bytes_before, bytes_after, file, line);
    if (!new_data && bytes_after > 0) {
        return false;
    }

    *data_ptr = new_data;
    *cap_ptr = new_cap;
    return true;
}

static inline void* md_array_set_capacity_internal(void* arr, size_t new_cap, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(alloc);
    uint8_t* p = arr ? (uint8_t*)md_array_header(arr) : 0;
    const size_t extra = sizeof(md_array_header_t);
    const size_t size = md_array_size(arr);
    const size_t old_cap = md_array_capacity(arr);
    const size_t bytes_before = arr ? item_size * old_cap + extra : 0;
    const size_t bytes_after = new_cap ? item_size * new_cap + extra : 0;

    uint8_t* new_p = (uint8_t*)alloc->realloc(alloc->inst, p, bytes_before, bytes_after, file, line);
    void*  new_arr = new_p ? new_p + extra : 0;

    if (new_arr) {
        md_array_header(new_arr)->size = size;
        md_array_header(new_arr)->capacity = new_cap;
    }

    return new_arr;
}

static inline bool md_array_raw_set_capacity_impl(void* data_addr, size_t* capacity, size_t* size, size_t new_capacity, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(data_addr);
    ASSERT(capacity);
    ASSERT(size);
    ASSERT(alloc);

    if (item_size == 0) {
        return false;
    }

    void** data = (void**)data_addr;
    const size_t old_capacity = *capacity;

    if (old_capacity > SIZE_MAX / item_size || new_capacity > SIZE_MAX / item_size) {
        return false;
    }

    const size_t bytes_before = old_capacity * item_size;
    const size_t bytes_after = new_capacity * item_size;

    void* new_data = alloc->realloc(alloc->inst, *data, bytes_before, bytes_after, file, line);
    if (!new_data && bytes_after) {
        return false;
    }

    *data = new_data;
    *capacity = new_capacity;
    if (*size > new_capacity) {
        *size = new_capacity;
    }

    return true;
}

static inline bool md_array_raw_count_add(size_t a, size_t b, size_t* out_sum) {
    ASSERT(out_sum);
    if (a > SIZE_MAX - b) {
        return false;
    }
    *out_sum = a + b;
    return true;
}

static inline size_t md_array_raw_grow_cap(size_t capacity, size_t min_capacity) {
    size_t new_capacity = capacity ? capacity : 8;
    while (new_capacity < min_capacity) {
        if (new_capacity > SIZE_MAX / 2) {
            return min_capacity;
        }
        new_capacity *= 2;
    }
    return new_capacity;
}

static inline bool md_array_raw_ensure_impl(void* data_addr, size_t* capacity, size_t size, size_t min_capacity, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(capacity);
    ASSERT(size <= *capacity);

    if (*capacity >= min_capacity) {
        return true;
    }

    const size_t new_capacity = md_array_raw_grow_cap(*capacity, min_capacity);
    size_t tmp_size = size;
    return md_array_raw_set_capacity_impl(data_addr, capacity, &tmp_size, new_capacity, item_size, alloc, file, line);
}

static inline bool md_array_raw_push_impl(void* data_addr, size_t* capacity, size_t* size, size_t item_size, const void* item, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(size);
    ASSERT(item);

    size_t needed = 0;
    if (!md_array_raw_count_add(*size, 1, &needed)) {
        return false;
    }

    if (!md_array_raw_ensure_impl(data_addr, capacity, *size, needed, item_size, alloc, file, line)) {
        return false;
    }

    void** data = (void**)data_addr;
    MEMCPY((uint8_t*)(*data) + (*size) * item_size, item, item_size);
    *size += 1;
    return true;
}

static inline void* md_array_raw_extend_impl(void* data_addr, size_t* capacity, size_t* size, size_t item_size, size_t count, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(size);

    size_t old_size = *size;
    size_t needed = 0;
    if (!md_array_raw_count_add(old_size, count, &needed)) {
        return NULL;
    }

    if (!md_array_raw_ensure_impl(data_addr, capacity, old_size, needed, item_size, alloc, file, line)) {
        return NULL;
    }

    void** data = (void**)data_addr;
    *size = needed;
    return (uint8_t*)(*data) + old_size * item_size;
}

static inline void* md_array_raw_extend_assert_impl(void* data_addr, size_t* capacity, size_t* size, size_t item_size, size_t count, struct md_allocator_i* alloc, const char* file, size_t line) {
    void* ptr = md_array_raw_extend_impl(data_addr, capacity, size, item_size, count, alloc, file, line);
    ASSERT(ptr || count == 0);
    return ptr;
}

static inline void* md_array_create_internal(size_t size, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    void* arr = md_array_set_capacity_internal(NULL, size, item_size, alloc, file, line);
    if (arr) {
        md_array_header(arr)->size = size;
    }
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
