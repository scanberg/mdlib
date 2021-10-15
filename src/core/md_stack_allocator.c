#include "md_stack_allocator.h"
#include <string.h>

void* stack_realloc(struct md_allocator_o* alloc, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    md_stack_allocator_t* stack_alloc = (md_stack_allocator_t*)alloc;

    ASSERT(stack_alloc);
    ASSERT(stack_alloc->magic == MD_STACK_ALLOCATOR_MAGIC);

    if (new_size == 0) {
        md_stack_allocator_free(stack_alloc, ptr, old_size);
        return NULL;
    }

    if (old_size) {
        ASSERT(ptr);
        if (ptr == stack_alloc->buf + stack_alloc->prv) {
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            const int64_t new_cur = stack_alloc->cur + diff;
            ASSERT(0 <= new_cur && new_cur < (int64_t)stack_alloc->cap);
            stack_alloc->cur = new_cur;
            return ptr;
        }
        void* new_ptr = md_stack_allocator_alloc(stack_alloc, new_size);
        memcpy(new_ptr, ptr, old_size);
        return new_ptr;
    }

    return md_stack_allocator_alloc(stack_alloc, new_size);
}

md_allocator_i md_stack_allocator_create_interface(md_stack_allocator_t* stack_alloc) {
    md_allocator_i alloc = {
        .inst = (md_allocator_o*)stack_alloc,
        .realloc = stack_realloc
    };
    return alloc;
}
