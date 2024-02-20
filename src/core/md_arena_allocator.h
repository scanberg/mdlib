#pragma once

#include <stddef.h>

struct md_allocator_i;

#define MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE (64*1024)

typedef struct md_vm_arena_temp_t {
    struct md_allocator_i* vm_arena;
    size_t pos;
} md_vm_arena_temp_t;

#ifdef __cplusplus
extern "C" {
#endif

// This is a simple allocator which does not free any memory until reset or destroy is called, and then all memory is freed at once.
// Internally it allocates memory in pages from a backing allocator to keep up with the demand. All pages are freed upon reset or destroy.
// This composes nicely with other allocators as it uses them as backing allocators

struct md_allocator_i* md_arena_allocator_create(struct md_allocator_i* backing, size_t page_size);
void md_arena_allocator_destroy(struct md_allocator_i* arena);

void* md_arena_allocator_push(struct md_allocator_i* arena, size_t size);
void* md_arena_allocator_push_zero(struct md_allocator_i* arena, size_t size);
void* md_arena_allocator_push_aligned(struct md_allocator_i* arena, size_t size, size_t alignment);

// Pop does not take into consideration the alignment of the allocation
// So there may be some bytes that are not fully recovered
void md_arena_allocator_pop(struct md_allocator_i* arena, size_t size);

void md_arena_allocator_reset(struct md_allocator_i* arena);

// Specialized arena allocator which uses resorts to the os virtual memory allocator
// Should be prefered whenever possible over the very generic version above

struct md_allocator_i* md_vm_arena_create(size_t reservation_size);
void md_vm_arena_destroy(struct md_allocator_i* vm_arena);

//void md_vm_arena_init(md_vm_arena_t* vm_arena, size_t reservation_size);
//void md_vm_arena_free(md_vm_arena_t* vm_arena);

void* md_vm_arena_push(struct md_allocator_i* vm_arena, size_t size);
void* md_vm_arena_push_zero(struct md_allocator_i* vm_arena, size_t size);
void* md_vm_arena_push_aligned(struct md_allocator_i* vm_arena, size_t size, size_t align);
void  md_vm_arena_pop(struct md_allocator_i* vm_arena, size_t size);

void  md_vm_arena_reset(struct md_allocator_i* vm_arena);

void   md_vm_arena_set_pos(struct md_allocator_i* vm_arena, size_t pos);
size_t md_vm_arena_get_pos(struct md_allocator_i* vm_arena);

//struct md_allocator_i md_vm_arena_create_interface(md_vm_arena_t* arena);

// Convenience functions for storing the position upon begin and resetting to that position upon end
md_vm_arena_temp_t md_vm_arena_temp_begin(struct md_allocator_i* vm_arena);
void md_vm_arena_temp_end(md_vm_arena_temp_t temp);

#ifdef __cplusplus
}
#endif
