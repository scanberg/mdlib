#include "md_allocator.h"
#include "core/compiler.h"
#include "core/common.h"
#include <stdlib.h>

typedef struct ring_buffer {
    uint64_t curr;
    uint64_t prev;
    char mem[MEGABYTES(1)];
} ring_buffer;

#ifdef MD_COMPILER_MSVC
__declspec(thread) ring_buffer ring;
#else
_Thread_local ring_buffer ring;
#endif

static void* _realloc(struct md_allocator_i *a, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)a;
    (void)old_size;
    (void)file;
    (void)line;
    return realloc(ptr, (size_t)new_size);
}

static void* _ring_alloc(struct md_allocator_i *a, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)a;
    (void)ptr;
    (void)old_size;
    (void)file;
    (void)line;
    if (new_size > sizeof(ring.mem)) {
        return NULL;
    }
    else if (new_size > (sizeof(ring.mem) - ring.curr)) {
        ring.curr = new_size;
        ring.prev = 0;
        return &ring.mem[0];
    }
    else {
        ring.prev = ring.curr;
        ring.curr += new_size;
        return &ring.mem[ring.prev];
    }
}

static struct md_allocator_i _default_allocator = {
    NULL,
    _realloc
};

static struct md_allocator_i _default_temp_allocator = {
    NULL,
    _ring_alloc
};

struct md_allocator_i* default_allocator = &_default_allocator;
struct md_allocator_i* default_temp_allocator = &_default_temp_allocator;
