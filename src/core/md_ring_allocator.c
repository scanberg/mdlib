#include <core/md_ring_allocator.h>

#include <core/md_allocator.h>
#include <core/md_common.h>

#include <stdint.h>

#define DEFAULT_ALIGNMENT 16
#define RING_MAGIC 0x8a7ccbba82745852
#define HEADER_SIZE (sizeof(ring_alloc_t) + sizeof(md_allocator_i))
#define VALIDATE_AND_EXTRACT_RING(alloc) \
    ASSERT(alloc); \
    ring_alloc_t* ring = (ring_alloc_t*)alloc->inst; \
    ASSERT(ring && ring->magic == RING_MAGIC); \

typedef struct ring_alloc_t {
    size_t   pos;
    size_t   cap;
    uint64_t magic;
    void*    ptr;
} ring_alloc_t;

static inline void* ring_push(ring_alloc_t* ring, size_t size, size_t align) {
    ASSERT(IS_POW2(align));
    ASSERT(HEADER_SIZE + size + align < ring->cap);

    size_t pos_addr = (size_t)ring->ptr + ring->pos;
    size_t alignment_size = ALIGN_TO(pos_addr, align) - pos_addr;
    size_t pos;
    if (ring->pos + alignment_size + size <= ring->cap) {
        pos = ring->pos + alignment_size;
    } else {
        // Start from beginning and align
        pos_addr = (size_t)ring->ptr + HEADER_SIZE;
        pos = ALIGN_TO(pos_addr, align) - pos_addr;
    }
    ring->pos = pos + size;
    return (char*)ring->ptr + pos;
}

static inline void ring_set_pos(ring_alloc_t* ring, size_t pos) {
    ASSERT(HEADER_SIZE <= pos && pos <= ring->cap);
    ring->pos = pos;
}

void* ring_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)inst;    
    (void)file;
    (void)line;
    ring_alloc_t* ring = (ring_alloc_t*)inst;
    ASSERT(ring && ring->magic == RING_MAGIC);

    if (new_size == 0) {
        // free
        // If this is the last allocation, then we move the pointer back
        if ((char*)ring->ptr + ring->pos == (char*)ptr + old_size) {
            ring_set_pos(ring, ring->pos - old_size);
        }
        return NULL;
    }

    if (old_size) {
        ASSERT(ptr);
        // realloc
        if ((char*)ring->ptr + ring->pos == (char*)ptr + old_size) {
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            const int64_t new_cur = ring->pos + diff;
            ASSERT(0 <= new_cur && new_cur < (int64_t)ring->cap);
            ring->pos = new_cur;
            return ptr;
        }
        size_t alignment = new_size <= 2 ? new_size : DEFAULT_ALIGNMENT;
        void* new_ptr = ring_push(ring, new_size, alignment);
        MEMCPY(new_ptr, ptr, old_size);
        return new_ptr;
    }

    // alloc
    size_t alignment = new_size <= 2 ? new_size : DEFAULT_ALIGNMENT;
    return ring_push(ring, new_size, alignment);
}

md_allocator_i* md_ring_allocator_create(void* backing_buffer, size_t buffer_capacity) {
    ASSERT(backing_buffer);
    ASSERT(buffer_capacity > HEADER_SIZE);

    ring_alloc_t* ring = backing_buffer;
    ring->pos = HEADER_SIZE;
    ring->cap = buffer_capacity;
    ring->magic = RING_MAGIC;
    ring->ptr = backing_buffer;

    md_allocator_i* alloc = (md_allocator_i*)((char*)backing_buffer + sizeof(ring_alloc_t));
    alloc->inst = (md_allocator_o*)ring;
    alloc->realloc = ring_realloc;

    return alloc;
}

void* md_ring_allocator_push_aligned(md_allocator_i* alloc, size_t size, size_t align) {
    VALIDATE_AND_EXTRACT_RING(alloc)
    return ring_push(ring, size, align);
}

void* md_ring_allocator_push(md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_RING(alloc)
    size_t align = size <= 2 ? size : DEFAULT_ALIGNMENT;
    return ring_push(ring, size, align);
}

void md_ring_allocator_pop(md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_RING(alloc)
    ring_set_pos(ring, ring->pos - size);
}

void md_ring_allocator_reset(md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_RING(alloc)
    ring_set_pos(ring, HEADER_SIZE);
}

void md_ring_allocator_set_pos_back(md_allocator_i* alloc, size_t pos) {
    VALIDATE_AND_EXTRACT_RING(alloc)
    ring_set_pos(ring, pos);
}

uint64_t md_ring_allocator_get_pos(md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_RING(alloc)
    return ring->pos;
}
