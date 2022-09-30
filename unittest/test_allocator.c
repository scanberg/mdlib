#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_pool_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_virtual_allocator.h>
#include <core/md_stack_allocator.h>

#define COMMON_ALLOCATOR_TEST_BODY \
    void* mem = 0; \
    mem = md_alloc(alloc, 16); \
    EXPECT_NE(mem, NULL); \
    EXPECT_EQ(md_free(alloc, mem, 16), NULL); \
    \
    int64_t* arr = NULL; \
    for (int64_t i = 0; i < 1000; ++i) { \
        md_array_push(arr, i, alloc); \
    } \
    \
    for (int64_t i = 0; i < 1000; ++i) { \
        ASSERT_EQ(arr[i], i); \
    } \
    md_array_free(arr, alloc); \
    \
    int64_t size[] = {16, 7238, 1, 2, 7, 3, 2, 4}; \
    for (int64_t i = 0; i < ARRAY_SIZE(size); ++i) { \
        uint64_t expected_alignment = size[i] > 2 ? 16 : size[i]; \
        mem = md_alloc(alloc, size[i]); \
        EXPECT_EQ((uint64_t)mem % expected_alignment, 0ULL); \
        md_free(alloc, mem, size[i]); \
    }

// COMMON TESTS
UTEST(allocator, default) {
    md_allocator_i* alloc = default_allocator;
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, default_temp) {
    md_allocator_i* alloc = default_temp_allocator;
    COMMON_ALLOCATOR_TEST_BODY
}

UTEST(allocator, arena) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);
    COMMON_ALLOCATOR_TEST_BODY
    md_arena_allocator_destroy(alloc);
}

UTEST(allocator, stack_generic) {
    void* buf = md_alloc(default_allocator, MEGABYTES(32));
    md_stack_allocator_t stack;
    md_stack_allocator_init(&stack, buf, MEGABYTES(32));
    md_allocator_i stack_alloc = md_stack_allocator_create_interface(&stack);
    md_allocator_i* alloc = &stack_alloc;
    COMMON_ALLOCATOR_TEST_BODY
    md_free(default_allocator, buf, MEGABYTES(32));
}

// @NOTE: Pool is an outlier here, since it is meant for allocations of a fixed size, thus cannot be tested with the common allocator test

UTEST(allocator, pool) {
    md_allocator_i* pool = md_pool_allocator_create(default_allocator, sizeof(uint64_t));

    uint64_t **items = {0};

    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            uint64_t *item = *md_array_push(items, md_alloc(pool, sizeof(uint64_t)), default_allocator);
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

    md_array_free(items, default_allocator);
}

UTEST(allocator, arena_extended) {
    md_allocator_i* arena = md_arena_allocator_create(default_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);

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
        void* mem = md_vm_arena_push(&arena, sizeof(int64_t));
    }

    md_vm_arena_reset(&arena);

    for (uint32_t i = 0; i < 1000; ++i) {
        void* mem = md_vm_arena_push_aligned(&arena, 128, 32);
    }

    md_vm_arena_free(&arena);
}


UTEST(allocator, vm_allocator) {
    md_vm_allocator_t vm;
    md_vm_allocator_init(&vm, GIGABYTES(64)); // Reserve ridiculous amounts of space
    md_allocator_i alloc = md_vm_allocator_create_interface(&vm);

    // Commit to smaller portion by allocating
    for (uint32_t i = 0; i < 1000; ++i) {
        md_alloc(&alloc, MEGABYTES(1));
    }

    // Release reserved space
    md_vm_allocator_free(&vm);
}

UTEST(allocator, stack_allocator) {
    void* buf = md_alloc(default_allocator, MEGABYTES(64));
    md_stack_allocator_t stack;
    md_stack_allocator_init(&stack, buf, MEGABYTES(64));

    for (uint32_t i = 0; i < 1000; ++i) {
        md_stack_allocator_push(&stack, sizeof(uint64_t));
    }

    for (uint32_t i = 0; i < 1000; ++i) {
        md_stack_allocator_push_aligned(&stack, sizeof(uint64_t), 32);
    }

    md_free(default_allocator, buf, MEGABYTES(64));
}