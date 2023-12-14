#pragma once

#include <stddef.h>

typedef struct md_pool_allocator_t md_pool_allocator_t;
struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// The pool allocator is ideal for allocating and freeing objects of the same size (or less than the defined object_size)
// It internally allocates larger pages to restrict the strain on the backing allocator.
// NOTE: YOU CAN ONLY ALLOCATE OBJECT_SIZE OR LESS FROM THIS ALLOCATOR, OTHERWISE ASSERTIONS WILL FIRE!
// NOTE: This allocator is not well designed, don't use it.

struct md_allocator_i* md_pool_allocator_create(struct md_allocator_i* backing, size_t slot_size_bytes);
void md_pool_allocator_destroy(struct md_allocator_i* a);

#ifdef __cplusplus
}
#endif
