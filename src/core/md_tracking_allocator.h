#ifndef _MD_TRACKING_ALLOCATOR_H_
#define _MD_TRACKING_ALLOCATOR_H_

#include <stdint.h>

struct md_allocator;

#ifdef __cplusplus
extern "C" {
#endif

// This is a tracking allocator which tracks every allocation and free to make sure that the correct amount of memory is freed and that non of the memory is leaked upon destruction.
struct md_allocator* md_tracking_allocator_create(struct md_allocator* backing);
void md_tracking_allocator_destroy(struct md_allocator* tracking_alloc);

#ifdef __cplusplus
}
#endif

#endif