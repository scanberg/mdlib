#ifndef _ARENA_ALLOC_H_
#define _ARENA_ALLOC_H_

// Implementation based on ginger bills

#ifndef DEFAULT_ALIGNMENT
#define DEFAULT_ALIGNMENT (2*sizeof(void *))
#endif

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#endif

inline uintptr_t align_forward(uintptr_t ptr, size_t align) {
    ASSERT((align & (align - 1)) == 0); // Check that alignment is a power of two (otherwise the mod trick bellow will not work)
    uintptr_t mod = ptr & ((uintptr_t)align - 1);
    return mod ? ptr + (uintptr_t)align - mod : ptr;
}

typedef struct mem_arena_t {
    unsigned char* buf;
    size_t         buf_len;
    size_t         offset;
} mem_arena_t;

inline void arena_init(mem_arena_t* a, void* backing_buffer, size_t backing_buffer_length) {
    a->buf = (unsigned char*)backing_buffer;
    a->buf_len = backing_buffer_length;
    a->offset = 0;
}

inline void* arena_alloc_aligned(mem_arena_t* a, size_t size, size_t align) {
    const bool is_power_of_two = (align & (align - 1)) == 0;
    if (!is_power_of_two) {

        return NULL;
    }

    // Align 'curr_offset' forward to the specified alignment
    uintptr_t curr_ptr = (uintptr_t)a->buf + (uintptr_t)a->offset;
    uintptr_t offset = align_forward(curr_ptr, align);
    offset -= (uintptr_t)a->buf; // Change to relative offset

    // Check to see if the backing memory has space left
    if (offset + size <= a->buf_len) {
        void* ptr = &a->buf[offset];
        a->offset = offset + size;

        // Zero new memory by default
        memset(ptr, 0, size);
        return ptr;
    }
    // Return NULL if the arena is out of memory (or handle differently)
    return NULL;
}

inline void* arena_alloc(mem_arena_t* a, size_t size) {
    return arena_alloc_aligned(a, size, DEFAULT_ALIGNMENT);
}

inline void arena_free(mem_arena_t* a, void* mem) {
    (void)a;
    (void)mem;
    // Nothing!
}

#endif