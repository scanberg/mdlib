#ifndef _MD_POOL_ALLOCATOR_H_
#define _MD_POOL_ALLOCATOR_H_

#include <core/md_allocator.h>

typedef struct md_pool_allocator_t md_pool_allocator_t;

#ifdef __cplusplus
extern "C" {
#endif

// The pool allocator is ideal for allocating and freeing objects of the same size (or less than the defined object_size)
// It internally allocates larger pages to restrict the strain on the backing allocator.
// NOTE: YOU CAN ONLY ALLOCATE OBJECT_SIZE OR LESS FROM THIS ALLOCATOR, OTHERWISE ASSERTIONS WILL FIRE!

struct md_allocator_i* md_pool_allocator_create(struct md_allocator_i* backing, uint64_t slot_size_bytes);
void md_pool_allocator_destroy(struct md_allocator_i* a);

#ifdef __cplusplus
}
#endif

#endif
