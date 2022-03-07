#pragma once

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// This is a tracking allocator which tracks every allocation and free to make sure that the correct amount of memory is freed and that non of the memory is leaked upon destruction.
struct md_allocator_i* md_tracking_allocator_create(struct md_allocator_i* backing);
void md_tracking_allocator_destroy(struct md_allocator_i* tracking_alloc);

#ifdef __cplusplus
}
#endif
