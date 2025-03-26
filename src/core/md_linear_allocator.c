#include <core/md_linear_allocator.h>
#include <core/md_allocator.h>
#include <core/md_common.h>

#include <stdint.h>

#define DEFAULT_ALIGNMENT 16
#define LINEAR_MAGIC 0x8a7ccbba81980234
#define HEADER_SIZE (sizeof(linear_alloc_t) + sizeof(md_allocator_i))
#define VALIDATE_AND_EXTRACT_LINEAR(alloc) \
    ASSERT(alloc); \
    linear_alloc_t* linear = (linear_alloc_t*)alloc->inst; \
    ASSERT(linear && linear->magic == LINEAR_MAGIC);

typedef struct linear_alloc_t {
    size_t pos;
    size_t cap;
    uint64_t magic;
    void* ptr;
} linear_alloc_t;

static inline void* linear_push(linear_alloc_t* linear, size_t size, size_t alignment) {
    void* mem = 0;
    size_t pos_addr = (size_t)linear->ptr + linear->pos;
    size_t alignment_size = ALIGN_TO(pos_addr, alignment) - pos_addr;

    if (linear->pos + alignment_size + size <= linear->cap) {
        mem = (char*)linear->ptr + linear->pos + alignment_size;
        linear->pos += alignment_size + size;
    }
    return mem;
}

static inline void linear_set_pos(linear_alloc_t* linear, size_t pos) {
    ASSERT(HEADER_SIZE <= pos && pos <= linear->pos);
    linear->pos = pos;
}

static void* linear_realloc(struct md_allocator_o* alloc, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file;
    (void)line;
    linear_alloc_t* linear = (linear_alloc_t*)alloc;
    ASSERT(linear && linear->magic == LINEAR_MAGIC);

    // Realloc or Free
    if (old_size) {
        ASSERT(ptr);
        if ((char*)linear->ptr + linear->pos == (char*)ptr + old_size) {
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            const int64_t new_cur = linear->pos + diff;
            ASSERT(0 <= new_cur && new_cur < (int64_t)linear->cap);
            linear->pos = new_cur;
            return ptr;
        }

        if (new_size == 0) {
            // Free
            return NULL;
        }

        // Realloc
        size_t alignment = new_size <= 2 ? new_size : DEFAULT_ALIGNMENT;
        void* new_ptr = linear_push(linear, new_size, alignment);
        MEMCPY(new_ptr, ptr, old_size);
        return new_ptr;
    }

    // Alloc
    size_t alignment = new_size <= 2 ? new_size : DEFAULT_ALIGNMENT;
    return linear_push(linear, new_size, alignment);
}

md_allocator_i* md_linear_allocator_create(void* backing_buffer, size_t buffer_capacity) {
    ASSERT(buffer_capacity > HEADER_SIZE);
    linear_alloc_t* linear = backing_buffer;
    linear->pos = HEADER_SIZE;
    linear->cap = buffer_capacity;
    linear->magic = LINEAR_MAGIC;
    linear->ptr = backing_buffer;

    md_allocator_i* alloc = (md_allocator_i*)((char*)backing_buffer + sizeof(linear_alloc_t));
    alloc->inst = (md_allocator_o*)linear;
    alloc->realloc = linear_realloc;

    return alloc;
}

size_t md_linear_allocator_avail_bytes(md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc)
        return (size_t)MAX(0, (int64_t)linear->cap - (int64_t)linear->pos);
}

void* md_linear_allocator_push_aligned(md_allocator_i* alloc, size_t size, size_t align) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc)
    ASSERT(IS_POW2(align));
    return linear_push(linear, size, align);
}

void* md_linear_allocator_push(md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc);
    size_t alignment = size <= 2 ? size : DEFAULT_ALIGNMENT;
    return linear_push(linear, size, alignment);
}

void md_linear_allocator_pop(md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc);
    linear_set_pos(linear, linear->pos - size);
}

void md_linear_allocator_reset(md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc);
    linear_set_pos(linear, HEADER_SIZE);
}

void md_linear_allocator_set_pos_back(md_allocator_i* alloc, size_t pos) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc);
    linear_set_pos(linear, pos);
}

size_t md_linear_allocator_get_pos(md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_LINEAR(alloc);
    return linear->pos;
}
