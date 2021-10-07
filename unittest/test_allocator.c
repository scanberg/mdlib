#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_array.inl>
#include <core/md_pool_allocator.h>
#include <core/md_arena_allocator.h>

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
