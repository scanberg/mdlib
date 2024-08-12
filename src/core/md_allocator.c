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

size_t md_temp_allocator_max_allocation_size(void) {
    return MAX_TEMP_ALLOCATION_SIZE;
}

THREAD_LOCAL md_allocator_i* _ring_alloc = 0;
THREAD_LOCAL void* _ring_buf = 0;

static void release_ring_buffer(void* data) {
    (void)data;
    if (_ring_buf) {
        free(_ring_buf);
        _ring_buf = 0;
        _ring_alloc = 0;
    }
}

static inline md_allocator_i* md_thread_ring_allocator(void) {
    if (!_ring_alloc) {
        ASSERT(!_ring_buf);
        _ring_buf = malloc(MD_TEMP_ALLOC_SIZE);
        _ring_alloc = md_ring_allocator_create(_ring_buf, MD_TEMP_ALLOC_SIZE);
        md_thread_on_exit(release_ring_buffer);
    }
    ASSERT(_ring_alloc);
    ASSERT(_ring_alloc->inst);
    return _ring_alloc;
}

static void* realloc_internal(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
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

// Ugly^2 indirection here just to keep the temp API somewhat homogeneous
static void* ring_realloc_internal(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)inst;
    md_allocator_i* thread_ring = md_thread_ring_allocator();
    return thread_ring->realloc(thread_ring->inst, ptr, old_size, new_size, file, line);
}

static struct md_allocator_i _heap_allocator = {
    NULL,
    realloc_internal,
};

// Get general allocator interface to heap (malloc)
md_allocator_i* md_get_heap_allocator(void) {
    return &_heap_allocator;
}

// Get general allocator interface to thread local ring buffer
md_allocator_i* md_get_temp_allocator(void) {
    return md_thread_ring_allocator();
}

void* md_temp_push(size_t bytes) {
    return md_ring_allocator_push(md_thread_ring_allocator(), bytes);
}

void* md_temp_push_aligned(size_t bytes, size_t alignment) {
    return md_ring_allocator_push_aligned(md_thread_ring_allocator(), bytes, alignment);
}

void md_temp_pop (size_t bytes) {
    md_ring_allocator_pop(md_thread_ring_allocator(), bytes);
}

size_t md_temp_get_pos(void) {
    return md_ring_allocator_get_pos(md_thread_ring_allocator());
}

void md_temp_set_pos_back(size_t pos) {
    md_ring_allocator_set_pos_back(md_thread_ring_allocator(), pos);
}
