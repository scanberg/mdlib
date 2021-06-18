#ifndef _MD_ARENA_ALLOCATOR_H_
#define _MD_ARENA_ALLOCATOR_H_

#include <stdint.h>

struct md_allocator;

#define MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE 4096

#ifdef __cplusplus
extern "C" {
#endif

// This is a simple allocator which does not free any memory until reset or destroy is called, and then all memory is freed in one go.
// Internally it allocates memory in pages from a backing allocator to keep up with the demand. All pages are of course freed when reset or destroy is called.
struct md_allocator* md_arena_allocator_create(struct md_allocator* backing, uint64_t page_size);
void md_arena_allocator_reset(struct md_allocator* arena);
void md_arena_allocator_destroy(struct md_allocator* arena);

#ifdef __cplusplus
}
#endif

#endif