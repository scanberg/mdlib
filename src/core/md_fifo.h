#pragma once

#include <core/md_common.h>
#include <core/md_allocator.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

typedef struct md_fifo_t {
    size_t capacity;
    size_t size;
    size_t head;
    size_t tail;
    int* data;
    md_allocator_i* alloc;
} md_fifo_t;

#ifdef __cplusplus
extern "C" {
#endif

bool md_fifo_empty(const md_fifo_t* fifo);
bool md_fifo_full(const md_fifo_t* fifo);

md_fifo_t md_fifo_create(size_t capacity, md_allocator_i* alloc);
void md_fifo_free(md_fifo_t* fifo);
void md_fifo_clear(md_fifo_t* fifo);

void md_fifo_push(md_fifo_t* fifo, int value);
int  md_fifo_pop(md_fifo_t* fifo);

#ifdef __cplusplus
}
#endif

// -----------------------------------------------------------------------------
// ORTHOGONAL TYPED CIRCULAR FIFO API (RAW + STRUCT WRAPPERS)
// -----------------------------------------------------------------------------
// Generalization of the fixed int-based md_fifo_t above to arbitrary POD types,
// mirroring the RAW/STRUCT split used by md_array.h.
//
// RAW tier:
//   Operates directly on lvalues at the callsite:
//     capacity, size, head, tail, data, alloc
//
// STRUCT tier:
//   Convenience wrappers over RAW for structs with fields:
//     capacity, size, head, tail, data, alloc
//
// Notes:
// - The RAW operations are strongly typed through 'data' (T*).
// - Use TRY variants when allocation failure should be handled gracefully.
// - Non-TRY variants ASSERT on allocation failure.
//
// Example:
//   MD_FIFO_RAW_DEFINE_TYPE(int_fifo_t, int);
//   int_fifo_t q;
//   MD_FIFO_INIT(q, alloc);
//   MD_FIFO_PUSH(q, 42);
//   int v = MD_FIFO_FRONT(q);
//   MD_FIFO_POP(q);
//   MD_FIFO_FREE(q);

// Define a standalone typed fifo struct.
#define MD_FIFO_RAW_DEFINE_TYPE(name, type) \
    typedef struct name {                   \
        size_t capacity;                    \
        size_t size;                        \
        size_t head;                        \
        size_t tail;                        \
        type* data;                         \
        struct md_allocator_i* alloc;       \
    } name

// Convenience form for inline declarations.
#define MD_FIFO_RAW_TYPE(type) struct { size_t capacity; size_t size; size_t head; size_t tail; type* data; struct md_allocator_i* alloc; }

// RAW lifecycle & queries
#define MD_FIFO_RAW_INIT(capacity, size, head, tail, data) \
    do {                                                    \
        (capacity) = 0;                                     \
        (size) = 0;                                          \
        (head) = 0;                                          \
        (tail) = 0;                                          \
        (data) = NULL;                                       \
    } while (0)

#define MD_FIFO_RAW_INIT_WITH_ALLOC(capacity, size, head, tail, data, alloc, allocator) \
    do {                                                                                \
        MD_FIFO_RAW_INIT((capacity), (size), (head), (tail), (data));                    \
        (alloc) = (allocator);                                                           \
    } while (0)

#define MD_FIFO_RAW_SIZE(size)          (size)
#define MD_FIFO_RAW_CAPACITY(capacity)  (capacity)
#define MD_FIFO_RAW_DATA(data)          (data)
#define MD_FIFO_RAW_EMPTY(size)         ((size) == 0)
#define MD_FIFO_RAW_FULL(size, capacity) ((size) == (capacity))
#define MD_FIFO_RAW_CLEAR(size, head, tail) \
    do {                                     \
        (size) = 0;                          \
        (head) = 0;                          \
        (tail) = 0;                          \
    } while (0)

// RAW allocation primitives
#define MD_FIFO_RAW_TRY_ENSURE(capacity, size, head, tail, data, min_capacity, alloc) \
    md_fifo_raw_ensure_impl(&(data), &(capacity), &(head), &(tail), (size), (size_t)(min_capacity), sizeof(*(data)), (alloc), __FILE__, __LINE__)

#define MD_FIFO_RAW_ENSURE(capacity, size, head, tail, data, min_capacity, alloc) \
    do {                                                                            \
        const bool _md_fifo_ok = MD_FIFO_RAW_TRY_ENSURE((capacity), (size), (head), (tail), (data), (min_capacity), (alloc)); \
        ASSERT(_md_fifo_ok); (void)_md_fifo_ok;                                      \
    } while (0)

// RAW element operations
#define MD_FIFO_RAW_TRY_PUSH(capacity, size, head, tail, data, item, alloc) \
    md_fifo_raw_push_impl(&(data), &(capacity), &(head), &(tail), &(size), sizeof(*(data)), &(item), (alloc), __FILE__, __LINE__)

#define MD_FIFO_RAW_PUSH(capacity, size, head, tail, data, item, alloc) \
    do {                                                                  \
        const bool _md_fifo_ok = MD_FIFO_RAW_TRY_PUSH((capacity), (size), (head), (tail), (data), (item), (alloc)); \
        ASSERT(_md_fifo_ok); (void)_md_fifo_ok;                            \
    } while (0)

#define MD_FIFO_RAW_FRONT(data, tail) ((data)[(tail)])

#define MD_FIFO_RAW_POP(capacity, size, head, tail) \
    do {                                              \
        (tail) = ((tail) + 1) & ((capacity) - 1);      \
        (size) -= 1;                                   \
    } while (0)

#define MD_FIFO_RAW_TRY_FREE(capacity, size, head, tail, data, alloc) \
    md_fifo_raw_set_capacity_impl(&(data), &(capacity), &(head), &(tail), &(size), 0, sizeof(*(data)), (alloc), __FILE__, __LINE__)

#define MD_FIFO_RAW_FREE(capacity, size, head, tail, data, alloc) \
    do {                                                            \
        const bool _md_fifo_ok = MD_FIFO_RAW_TRY_FREE((capacity), (size), (head), (tail), (data), (alloc)); \
        ASSERT(_md_fifo_ok); (void)_md_fifo_ok;                      \
    } while (0)

// STRUCT convenience wrappers (expects fields: capacity, size, head, tail, data, alloc)
#define MD_FIFO_INIT(q, allocator) \
    MD_FIFO_RAW_INIT_WITH_ALLOC((q).capacity, (q).size, (q).head, (q).tail, (q).data, (q).alloc, (allocator))

#define MD_FIFO_SIZE(q)     MD_FIFO_RAW_SIZE((q).size)
#define MD_FIFO_CAPACITY(q) MD_FIFO_RAW_CAPACITY((q).capacity)
#define MD_FIFO_DATA(q)     MD_FIFO_RAW_DATA((q).data)
#define MD_FIFO_EMPTY(q)    MD_FIFO_RAW_EMPTY((q).size)
#define MD_FIFO_FULL(q)     MD_FIFO_RAW_FULL((q).size, (q).capacity)
#define MD_FIFO_CLEAR(q)    MD_FIFO_RAW_CLEAR((q).size, (q).head, (q).tail)

#define MD_FIFO_TRY_ENSURE(q, min_capacity) \
    MD_FIFO_RAW_TRY_ENSURE((q).capacity, (q).size, (q).head, (q).tail, (q).data, (min_capacity), (q).alloc)

#define MD_FIFO_ENSURE(q, min_capacity) \
    MD_FIFO_RAW_ENSURE((q).capacity, (q).size, (q).head, (q).tail, (q).data, (min_capacity), (q).alloc)

#define MD_FIFO_TRY_PUSH(q, item) \
    MD_FIFO_RAW_TRY_PUSH((q).capacity, (q).size, (q).head, (q).tail, (q).data, (item), (q).alloc)

#define MD_FIFO_PUSH(q, item) \
    MD_FIFO_RAW_PUSH((q).capacity, (q).size, (q).head, (q).tail, (q).data, (item), (q).alloc)

#define MD_FIFO_FRONT(q) MD_FIFO_RAW_FRONT((q).data, (q).tail)

#define MD_FIFO_POP(q) \
    MD_FIFO_RAW_POP((q).capacity, (q).size, (q).head, (q).tail)

#define MD_FIFO_TRY_FREE(q) \
    MD_FIFO_RAW_TRY_FREE((q).capacity, (q).size, (q).head, (q).tail, (q).data, (q).alloc)

#define MD_FIFO_FREE(q) \
    MD_FIFO_RAW_FREE((q).capacity, (q).size, (q).head, (q).tail, (q).data, (q).alloc)

#ifdef __cplusplus
extern "C" {
#endif

// Rounds up to the next power of two, with a minimum of 16 (so head/tail wrap-around masking works).
static inline size_t md_fifo_raw_grow_cap(size_t capacity, size_t min_capacity) {
    size_t new_capacity = capacity ? capacity : 16;
    while (new_capacity < min_capacity) {
        if (new_capacity > SIZE_MAX / 2) {
            return min_capacity;
        }
        new_capacity *= 2;
    }
    return new_capacity;
}

// Reallocates the backing buffer to new_capacity, linearizing any wrapped content so that
// tail becomes 0 and head becomes the (unwrapped) size. Shrinking below size is not supported.
static inline bool md_fifo_raw_set_capacity_impl(void* data_addr, size_t* capacity, size_t* head, size_t* tail, size_t* size, size_t new_capacity, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(data_addr);
    ASSERT(capacity);
    ASSERT(head);
    ASSERT(tail);
    ASSERT(size);
    ASSERT(alloc);

    if (item_size == 0) {
        return false;
    }
    if (new_capacity < *size) {
        return false;
    }

    void** data = (void**)data_addr;
    void* new_data = NULL;

    if (new_capacity > 0) {
        if (new_capacity > SIZE_MAX / item_size) {
            return false;
        }
        new_data = alloc->realloc(alloc->inst, NULL, 0, new_capacity * item_size, file, line);
        if (!new_data) {
            return false;
        }

        if (*size) {
            if (*tail < *head) {
                MEMCPY(new_data, (uint8_t*)(*data) + (*tail) * item_size, (*size) * item_size);
            } else {
                const size_t old_capacity = *capacity;
                const size_t right_count = old_capacity - *tail;
                MEMCPY(new_data, (uint8_t*)(*data) + (*tail) * item_size, right_count * item_size);
                MEMCPY((uint8_t*)new_data + right_count * item_size, *data, (*head) * item_size);
            }
        }
    }

    if (*data) {
        alloc->realloc(alloc->inst, *data, (*capacity) * item_size, 0, file, line);
    }

    *data = new_data;
    *capacity = new_capacity;
    *tail = 0;
    *head = (new_capacity > 0) ? (*size % new_capacity) : 0;

    return true;
}

static inline bool md_fifo_raw_ensure_impl(void* data_addr, size_t* capacity, size_t* head, size_t* tail, size_t size, size_t min_capacity, size_t item_size, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(capacity);
    ASSERT(size <= *capacity);

    if (*capacity >= min_capacity) {
        return true;
    }

    const size_t new_capacity = md_fifo_raw_grow_cap(*capacity, min_capacity);
    size_t tmp_size = size;
    return md_fifo_raw_set_capacity_impl(data_addr, capacity, head, tail, &tmp_size, new_capacity, item_size, alloc, file, line);
}

static inline bool md_fifo_raw_push_impl(void* data_addr, size_t* capacity, size_t* head, size_t* tail, size_t* size, size_t item_size, const void* item, struct md_allocator_i* alloc, const char* file, size_t line) {
    ASSERT(size);
    ASSERT(item);

    if (*size >= *capacity) {
        if (!md_fifo_raw_ensure_impl(data_addr, capacity, head, tail, *size, *size + 1, item_size, alloc, file, line)) {
            return false;
        }
    }

    void** data = (void**)data_addr;
    MEMCPY((uint8_t*)(*data) + (*head) * item_size, item, item_size);
    *head = (*head + 1) & (*capacity - 1);
    *size += 1;
    return true;
}

#ifdef __cplusplus
}
#endif
