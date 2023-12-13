#include <core/md_ring_allocator.h>

void* ring_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)inst;    
    (void)file;
    (void)line;
    md_ring_allocator_t* ring = (md_ring_allocator_t*)inst;
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);

    if (new_size == 0) {
        // free
        // If this is the last allocation, then we move the pointer back
        if ((char*)ring->ptr + ring->pos == (char*)ptr + old_size) {
            md_ring_allocator_pop(ring, old_size);
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
        void* new_ptr = md_ring_allocator_push(ring, new_size);
        MEMCPY(new_ptr, ptr, old_size);
        return new_ptr;
    }

    // alloc
    return md_ring_allocator_push(ring, new_size);
}

md_allocator_i md_ring_allocator_create_interface(md_ring_allocator_t* ring) {
    ASSERT(ring && ring->magic == MD_RING_ALLOCATOR_MAGIC);
    md_allocator_i alloc = {
        .inst = (md_allocator_o*)ring,
        .realloc = ring_realloc
    };
    return alloc;
}
