#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_arena_allocator.h>
#include <core/md_virtual_allocator.h>
#include <core/md_linear_allocator.h>
#include <core/md_ring_allocator.h>

#define COMMON_ALLOCATOR_TEST_BODY \
    void* mem = 0; \
    mem = md_alloc(alloc, 16); \
    EXPECT_NE(mem, NULL); \
    md_free(alloc, mem, 16); \
    \
    int64_t* arr = NULL; \
    for (int64_t i = 0; i < 1000; ++i) { \
        md_array_push(arr, i, alloc); \
        md_array_push(arr, i, alloc); \
        md_array_pop(arr); \
    } \
    \
    for (int64_t i = 0; i < 1000; ++i) { \
        ASSERT_EQ(arr[i], i); \
    } \
    md_array_free(arr, alloc); \
    \
    uint32_t size[] = {16, KILOBYTES(1), 1, 2, 7, 3, 2, 4, 11, 12, 13, 14, KILOBYTES(2), 17, 1, 2, 256, 1, 16, 7, 2, KILOBYTES(3), 23, 24, 25, 26, 27, 1, 29, 30, 4, 8}; \
    for (int64_t i = 0; i < ARRAY_SIZE(size); ++i) { \
        uint64_t expected_alignment = size[i] > 2 ? 16 : size[i]; \
        mem = md_alloc(alloc, size[i]); \
        EXPECT_EQ((uint64_t)mem % expected_alignment, 0ULL); \
        md_free(alloc, mem, size[i]); \
    }

// Common buffer for allocators to use
static char buf[KILOBYTES(16)];

// COMMON TESTS
UTEST(allocator, default) {
    md_allocator_i* alloc = md_get_heap_allocator();
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, default_temp) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    COMMON_ALLOCATOR_TEST_BODY
    md_temp_end(temp);
}

UTEST(allocator, arena) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);
    COMMON_ALLOCATOR_TEST_BODY
    md_arena_allocator_destroy(alloc);
}

UTEST(allocator, linear_generic) {
    md_allocator_i* alloc = md_linear_allocator_create(buf, sizeof(buf));
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, ring_generic) {
    md_allocator_i* alloc = md_ring_allocator_create(buf, sizeof(buf));
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, arena_extended) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);

    for (int j = 0; j < 1000; ++j) {
        EXPECT_EQ((uint64_t)md_alloc(arena, 16) % 16, 0ULL); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4)  % 16, 0ULL); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4)  % 16, 0ULL); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4)  % 16, 0ULL); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4)  % 16, 0ULL); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4)  % 16, 0ULL); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 1)  % 1,  0ULL); // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1)  % 1,  0ULL); // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1)  % 1,  0ULL); // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1)  % 1,  0ULL); // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1)  % 1,  0ULL); // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 2)  % 2,  0ULL); // Expect to be aligned to 2 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 1)  % 1,  0ULL); // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 2)  % 2,  0ULL); // Expect to be aligned to 2 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4)  % 16, 0ULL); // Expect to be aligned to 16 bytes

        EXPECT_NE(md_alloc(arena, 9000), NULL); // Exceeds page size, should still be good

        // Make sure we get some internal pages going.
        for (int i = 0; i < 20; ++i) {
            EXPECT_EQ((uint64_t)md_alloc(arena, 1024) % 16, 0ULL);
        }

        md_arena_allocator_reset(arena);
    }

    md_arena_allocator_destroy(arena);
}

UTEST(allocator, vm_arena) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    ASSERT_TRUE(arena);

    for (uint32_t i = 0; i < 1000; ++i) {
        size_t size = rand() % 63 + 1;
        void* mem = md_vm_arena_push(arena, size);
        MEMSET(mem, 0, size);   // Important to touch memory in order to validate that it is writable
    }

    md_vm_arena_reset(arena);

    for (uint32_t i = 0; i < 1000; ++i) {
        size_t size = rand() % 127 + 1;
        void* mem = md_vm_arena_push_aligned(arena, size, 32);
        MEMSET(mem, 0, size); // Important to touch memory in order to validate that it is writable
        md_vm_arena_pop(arena, size);
    }

    md_vm_arena_destroy(arena);
}

UTEST(allocator, vm_arena_generic) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(4));
    ASSERT_TRUE(alloc);
    COMMON_ALLOCATOR_TEST_BODY
    md_vm_arena_destroy(alloc);
}


UTEST(allocator, vm) {
    md_vm_allocator_t vm;
    md_vm_allocator_init(&vm, GIGABYTES(64)); // Reserve ridiculous amounts of space
    md_allocator_i alloc = md_vm_allocator_create_interface(&vm);

    // Commit to smaller portion by allocating
    for (uint32_t i = 0; i < 1000; ++i) {
        md_alloc(&alloc, rand() % 63 + 1);
    }

    // Release reserved space
    md_vm_allocator_free(&vm);
}

UTEST(allocator, linear) {
    md_allocator_i* linear = md_linear_allocator_create(buf, sizeof(buf));

    for (uint32_t i = 0; i < 1000; ++i) {
        md_linear_allocator_push(linear, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 500; ++i) {
        md_linear_allocator_pop(linear, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 1000; ++i) {
        md_linear_allocator_push_aligned(linear, sizeof(uint64_t), 32);
    }
}

UTEST(allocator, ring) {
    md_allocator_i* ring = md_ring_allocator_create(buf, sizeof(buf));

    for (uint32_t i = 0; i < 100000; ++i) {
        uint32_t size = i % 20 + 1;
        md_ring_allocator_push(ring, size);
    }

    for (uint32_t i = 0; i < 500; ++i) {
        md_ring_allocator_pop(ring, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 1000; ++i) {
        md_ring_allocator_push_aligned(ring, sizeof(uint64_t), 32);
    }
}

UTEST(allocator, temp) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);

    for (uint32_t i = 0; i < 100000; ++i) {
        uint32_t size = i % 20 + 1;
        md_temp_alloc(temp, size);
    }

    for (uint32_t i = 0; i < 1000; ++i) {
        md_temp_alloc_aligned(temp, sizeof(uint64_t), 32);
    }

    md_temp_end(temp);
}

UTEST(allocator, temp_scope) {
    md_temp_scope_t outer = md_temp_begin();
    md_allocator_i* outer_alloc = md_temp_allocator(outer);

    uint64_t* outer_data = md_temp_alloc_zero_array(outer, uint64_t, 4);
    outer_data[0] = 0x1234;

    md_allocator_i* conflicts[] = { outer_alloc };
    md_temp_scope_t inner = md_temp_begin_avoid(outer_alloc);
    EXPECT_NE(md_temp_allocator(inner), outer_alloc);

    void* aligned = md_temp_alloc_aligned(inner, 64, 64);
    EXPECT_EQ((uint64_t)aligned % 64, 0ULL);

    uint32_t* zero_data = md_temp_alloc_zero_array(inner, uint32_t, 8);
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(zero_data[i], 0U);
    }

    md_temp_end(inner);

    EXPECT_EQ(outer_data[0], 0x1234ULL);

    md_temp_end(outer);
}
