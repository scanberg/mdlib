#ifndef _MD_POOL_ALLOCATOR_H_
#define _MD_POOL_ALLOCATOR_H_

#include <stdint.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// The pool allocator is ideal for allocating and freeing objects of the same size (or less than the defined object_size)
// It internally allocates larger pages to restrict the strain on the backing allocator.
// NOTE: YOU CAN ONLY ALLOCATE OBJECT_SIZE OR LESS FROM THIS ALLOCATOR, OTHERWISE ASSERTIONS WILL FIRE!

struct md_allocator_i* md_create_pool_allocator(struct md_allocator_i* backing, uint32_t slot_size);
void md_destroy_pool_allocator(struct md_allocator_i* a);

#ifdef __cplusplus
}
#endif

#endif