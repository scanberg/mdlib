#include "md_linear_allocator.h"
#include <string.h>

static void* linear_realloc(struct md_allocator_o* alloc, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    md_linear_allocator_t* linear = (md_linear_allocator_t*)alloc;

    ASSERT(linear);
    ASSERT(linear->magic == MD_LINEAR_ALLOCATOR_MAGIC);

    if (new_size == 0) {
        if ((char*)linear->ptr + linear->pos == (char*)ptr + old_size) {
            md_linear_allocator_pop(linear, old_size);
        }
        return NULL;
    }

    if (old_size) {
        ASSERT(ptr);
        if ((char*)linear->ptr + linear->pos == (char*)ptr + old_size) {
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            const int64_t new_cur = linear->pos + diff;
            ASSERT(0 <= new_cur && new_cur < (int64_t)linear->cap);
            linear->pos = new_cur;
            return ptr;
        }
        void* new_ptr = md_linear_allocator_push(linear, new_size);
        MEMCPY(new_ptr, ptr, old_size);
        return new_ptr;
    }

    return md_linear_allocator_push(linear, new_size);
}

md_allocator_i md_linear_allocator_create_interface(md_linear_allocator_t* linear_alloc) {
    ASSERT(linear_alloc->magic == MD_LINEAR_ALLOCATOR_MAGIC);
    md_allocator_i alloc = {
        .inst = (md_allocator_o*)linear_alloc,
        .realloc = linear_realloc
    };
    return alloc;
}
