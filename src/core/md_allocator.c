#include <core/md_allocator.h>

#include <core/md_ring_allocator.h>
#include <core/md_common.h>
#include <core/md_os.h>

#include <stdlib.h>
#include <string.h>

#ifndef MD_TEMP_ALLOC_SIZE
#define MD_TEMP_ALLOC_SIZE MEGABYTES(16)
#endif 

#define MAX_TEMP_ALLOCATION_SIZE (MD_TEMP_ALLOC_SIZE / 2)

uint64_t default_temp_allocator_max_allocation_size() {
    return MAX_TEMP_ALLOCATION_SIZE;
}

THREAD_LOCAL md_ring_allocator_t _ring_alloc = {
    .pos = 0,
    .cap = MD_TEMP_ALLOC_SIZE,
    .magic = MD_RING_ALLOCATOR_MAGIC,
    .ptr = 0,
};


static void* realloc_internal(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)inst;
    (void)old_size;
    (void)file;
    (void)line;
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    return realloc(ptr, (size_t)new_size);
}

extern void* ring_realloc(struct md_allocator_o* inst, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line);

static void* ring_realloc_internal(struct md_allocator_o* inst, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)inst;
    return ring_realloc((md_allocator_o*)get_thread_ring_allocator(), ptr, old_size, new_size, file, line);
}

static void release_ring_buffer(void* data) {
    (void)data;
    ASSERT(_ring_alloc.ptr);
    if (_ring_alloc.ptr) {
        free(_ring_alloc.ptr);
    }
}

struct md_ring_allocator_t* get_thread_ring_allocator() {
    if (!_ring_alloc.ptr) {
        _ring_alloc.ptr = malloc(MD_TEMP_ALLOC_SIZE);
        ASSERT(_ring_alloc.ptr);
        md_thread_on_exit(release_ring_buffer);
    }
    return &_ring_alloc;
}

static struct md_allocator_i _default_allocator = {
    NULL,
    realloc_internal,
};

static struct md_allocator_i _default_temp_allocator = {
    NULL,
    ring_realloc_internal,
};

struct md_allocator_i* default_allocator = &_default_allocator;
struct md_allocator_i* default_temp_allocator = &_default_temp_allocator;
