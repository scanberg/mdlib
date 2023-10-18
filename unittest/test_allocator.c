#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_pool_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_virtual_allocator.h>
#include <core/md_linear_allocator.h>
#include <core/md_ring_allocator.h>

#define COMMON_ALLOCATOR_TEST_BODY \
    void* mem = 0; \
    mem = md_alloc(alloc, 16); \
    EXPECT_NE(mem, NULL); \
    EXPECT_EQ(md_free(alloc, mem, 16), NULL); \
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

struct allocator {
    void*  buf;
    size_t cap;
};

UTEST_F_SETUP(allocator) {
    utest_fixture->buf = md_alloc(md_heap_allocator, KILOBYTES(32));
	utest_fixture->cap = KILOBYTES(32);
}

UTEST_F_TEARDOWN(allocator) {
	md_free(md_heap_allocator, utest_fixture->buf, utest_fixture->cap);
}

// COMMON TESTS
UTEST(allocator, default) {
    md_allocator_i* alloc = md_heap_allocator;
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, default_temp) {
    md_allocator_i* alloc = md_temp_allocator;
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, arena) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);
    COMMON_ALLOCATOR_TEST_BODY
    md_arena_allocator_destroy(alloc);
}

UTEST_F(allocator, linear_generic) {
    md_linear_allocator_t linear;
    md_linear_allocator_init(&linear, utest_fixture->buf, utest_fixture->cap);
    md_allocator_i linear_alloc = md_linear_allocator_create_interface(&linear);
    md_allocator_i* alloc = &linear_alloc;
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST_F(allocator, ring_generic) {
    md_ring_allocator_t ring;
    md_ring_allocator_init(&ring, utest_fixture->buf, utest_fixture->cap);
    md_allocator_i ring_alloc = md_ring_allocator_create_interface(&ring);
    md_allocator_i* alloc = &ring_alloc;
    COMMON_ALLOCATOR_TEST_BODY
}


// @NOTE: Pool is an outlier here, since it is meant for allocations of a fixed size, thus cannot be tested with the common allocator test

UTEST(allocator, pool) {
    md_allocator_i* pool = md_pool_allocator_create(md_heap_allocator, sizeof(uint64_t));

    uint64_t **items = {0};

    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            uint64_t *item = *md_array_push(items, md_alloc(pool, sizeof(uint64_t)), md_heap_allocator);
            *item = i;
        }

        for (int i = 100; i < 1000; ++i) {
            md_free(pool, items[i], sizeof(uint64_t));
        }
        md_array_shrink(items, 100);

        int indices[10] = {1, 2, 5, 6, 70, 90, 18, 16, 12, 10};

        for (int i = 0; i < (int)ARRAY_SIZE(indices); ++i) {
            const int idx = indices[i];
            md_free(pool, items[idx], sizeof(uint64_t));
        }

        for (int i = 0; i < (int)ARRAY_SIZE(indices); ++i) {
            const int idx = indices[i];
            items[idx] = md_alloc(pool, sizeof(uint64_t));
        }
    }

    md_array_free(items, md_heap_allocator);
}

UTEST(allocator, arena_extended) {
    md_allocator_i* arena = md_arena_allocator_create(md_heap_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);

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
    md_vm_arena_t arena;
    md_vm_arena_init(&arena, GIGABYTES(4));

    for (uint32_t i = 0; i < 1000; ++i) {
        size_t size = rand() % 63 + 1;
        void* mem = md_vm_arena_push(&arena, size);
        memset(mem, 0, size);   // Important to touch memory in order to validate that it is writable
    }

    md_vm_arena_reset(&arena);

    for (uint32_t i = 0; i < 1000; ++i) {
        size_t size = rand() % 127 + 1;
        void* mem = md_vm_arena_push_aligned(&arena, size, 32);
        memset(mem, 0, size); // Important to touch memory in order to validate that it is writable
    }

    md_vm_arena_free(&arena);
}

UTEST(allocator, vm_arena_generic) {
    md_vm_arena_t arena;
    md_vm_arena_init(&arena, GIGABYTES(4));
    md_allocator_i vm_alloc = md_vm_arena_create_interface(&arena);
    md_allocator_i* alloc = &vm_alloc;

    COMMON_ALLOCATOR_TEST_BODY

    md_vm_arena_free(&arena);
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

UTEST_F(allocator, linear) {
    md_linear_allocator_t linear;
    md_linear_allocator_init(&linear, utest_fixture->buf, utest_fixture->cap);

    for (uint32_t i = 0; i < 1000; ++i) {
        md_linear_allocator_push(&linear, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 500; ++i) {
        md_linear_allocator_pop(&linear, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 1000; ++i) {
        md_linear_allocator_push_aligned(&linear, sizeof(uint64_t), 32);
    }
}

UTEST_F(allocator, ring) {
    md_ring_allocator_t ring;
    md_ring_allocator_init(&ring, utest_fixture->buf, utest_fixture->cap);

    for (uint32_t i = 0; i < 1000; ++i) {
        md_ring_allocator_push(&ring, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 500; ++i) {
        md_ring_allocator_pop(&ring, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 1000; ++i) {
        md_ring_allocator_push_aligned(&ring, sizeof(uint64_t), 32);
    }
}