#pragma once

#include <stdint.h>

struct md_allocator_i;

#define MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE (64*1024)

typedef struct md_vm_arena_t {
    void* base;     // Base adress of reserved memory space
    uint64_t size;  // Reservation size
    uint64_t pos;
    uint64_t commit_pos;
    uint64_t align;
    uint64_t magic;
} md_vm_arena_t;

typedef struct md_vm_arena_temp_t {
    md_vm_arena_t* arena;
    uint64_t pos;
} md_vm_arena_temp_t;

#ifdef __cplusplus
extern "C" {
#endif

// This is a simple allocator which does not free any memory until reset or destroy is called, and then all memory is freed at once.
// Internally it allocates memory in pages from a backing allocator to keep up with the demand. All pages are of course freed upon reset or destroy.

struct md_allocator_i* md_arena_allocator_create(struct md_allocator_i* backing, uint64_t page_size);
void md_arena_allocator_reset(struct md_allocator_i* arena);
void md_arena_allocator_destroy(struct md_allocator_i* arena);


// Specialized arena allocator which uses resorts to the os virtual memory allocator
// Should be prefered whenever possible over the very generic version above

void md_vm_arena_init(md_vm_arena_t* vm_arena, uint64_t reservation_size);
void md_vm_arena_free(md_vm_arena_t* vm_arena);

void* md_vm_arena_push(md_vm_arena_t* vm_arena, uint64_t size);
void* md_vm_arena_push_zero(md_vm_arena_t* vm_arena, uint64_t size);
void* md_vm_arena_push_aligned(md_vm_arena_t* vm_arena, uint64_t size, uint64_t align);
void  md_vm_arena_pop(md_vm_arena_t* vm_arena, uint64_t size);

void  md_vm_arena_reset(md_vm_arena_t* vm_arena);

void	 md_vm_arena_set_pos(md_vm_arena_t* vm_arena, uint64_t pos);
uint64_t md_vm_arena_get_pos(md_vm_arena_t* vm_arena);

struct md_allocator_i md_vm_arena_create_interface(md_vm_arena_t* arena);

// Convenience functions for storing the position upon begin and resetting to that position upon end
md_vm_arena_temp_t md_vm_arena_temp_begin(md_vm_arena_t* arena);
void md_vm_arena_temp_end(md_vm_arena_temp_t temp);

#ifdef __cplusplus
}
#endif
