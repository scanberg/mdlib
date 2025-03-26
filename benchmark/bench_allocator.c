#include "ubench.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_pool_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_virtual_allocator.h>
#include <core/md_linear_allocator.h>
#include <core/md_ring_allocator.h>

#define COMMON_ALLOCATOR_TEST_BODY \
    void* mem = 0; \
    uint32_t size[] = {16, 7238, 1, 2, 7, 3, 2, 4, 11, 12, 13, 14, 15, 17}; \
    for (uint64_t i = 0; i < ARRAY_SIZE(size); ++i) { \
        mem = md_alloc(alloc, size[i]); \
        md_free(alloc, mem, size[i]); \
    }

UBENCH_EX(allocator, default) {
    md_allocator_i* alloc = md_get_heap_allocator();
    UBENCH_DO_BENCHMARK() {
        COMMON_ALLOCATOR_TEST_BODY
    }
}
