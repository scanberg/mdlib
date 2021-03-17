#ifndef _MD_ARENA_ALLOCATOR_H_
#define _MD_ARENA_ALLOCATOR_H_

#include <stdint.h>

struct md_allocator_i;

#define MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE 4096

// This is a simple allocator which does not free any memory until reset or destroy is called, and then all memory is freed in one go.
// Internally it allocates memory in pages from a backing allocator to keep up with the demand. All pages are of course freed when reset or destroy is called.

struct md_allocator_i* md_arena_allocator_create(struct md_allocator_i* backing, uint64_t page_size);
void md_arena_allocator_reset(struct md_allocator_i* a);
void md_arena_allocator_destroy(struct md_allocator_i* a);

#endif