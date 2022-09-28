#pragma once

#include <stdint.h>

struct md_allocator_i;

#define MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE (64*1024)
#define MD_ARENA_ALLOCATOR_DEFAULT_CAPACITY (4 * 1024ULL * 1024ULL)

typedef struct md_arena_allocator_t md_arena_allocator_t;

typedef struct md_arena_temp_t md_arena_temp_t;
struct md_arena_temp_t {
	md_arena_allocator_t* arena;
	uint64_t pos;
};

#ifdef __cplusplus
extern "C" {
#endif

// This is a simple allocator which does not free any memory until reset or destroy is called, and then all memory is freed at once.
// Internally it allocates memory in pages from a backing allocator to keep up with the demand. All pages are of course freed upon reset or destroy.

struct md_allocator_i* md_arena_allocator_create(struct md_allocator_i* backing, uint64_t page_size);
void md_arena_allocator_reset(struct md_allocator_i* arena);
void md_arena_allocator_destroy(struct md_allocator_i* arena);

//struct md_allocator_i* md_arena_allocator_create_interface(md_arena_allocator_t* arena);

/*
md_arena_allocator_t* md_arena_allocator_create(uint64_t capacity);
void md_arena_allocator_destroy(md_arena_allocator_t* arena);

void* md_arena_push(md_arena_allocator_t* arena, uint64_t size);
void* md_arena_push_zero(md_arena_allocator_t* arena, uint64_t size);
void* md_arena_push_aligned(md_arena_allocator_t* arena, uint64_t size, uint64_t align);
void  md_arena_pop(md_arena_allocator_t* arena, uint64_t size);

void  md_arena_clear(md_arena_allocator_t* arena);

void	 md_arena_set_pos(struct md_arena_allocator_t* arena, uint64_t pos);
uint64_t md_arena_get_pos(struct md_arena_allocator_t* arena);

// Convenience functions for storing a position upon begin and resetting to that position upon end
md_arena_temp_t md_arena_temp_begin(md_arena_allocator_t* arena);
void md_arena_temp_end(md_arena_temp_t temp);
*/

#ifdef __cplusplus
}
#endif
