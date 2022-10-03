#include "md_stack_allocator.h"
#include <string.h>

static void* stack_realloc(struct md_allocator_o* alloc, void* ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    md_stack_allocator_t* stack = (md_stack_allocator_t*)alloc;

    ASSERT(stack);
    ASSERT(stack->magic == MD_STACK_ALLOCATOR_MAGIC);

    if (new_size == 0) {
        if ((char*)stack->ptr + stack->pos == (char*)ptr + old_size) {
            md_stack_allocator_pop(stack, old_size);
        }
        return NULL;
    }

    if (old_size) {
        ASSERT(ptr);
        if ((char*)stack->ptr + stack->pos == (char*)ptr + old_size) {
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            const int64_t new_cur = stack->pos + diff;
            ASSERT(0 <= new_cur && new_cur < (int64_t)stack->cap);
            stack->pos = new_cur;
            return ptr;
        }
        void* new_ptr = md_stack_allocator_push(stack, new_size);
        memcpy(new_ptr, ptr, old_size);
        return new_ptr;
    }

    return md_stack_allocator_push(stack, new_size);
}

md_allocator_i md_stack_allocator_create_interface(md_stack_allocator_t* stack_alloc) {
    md_allocator_i alloc = {
        .inst = (md_allocator_o*)stack_alloc,
        .realloc = stack_realloc
    };
    return alloc;
}
