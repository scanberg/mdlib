#include "utest.h"

#include <md_allocator.h>
#include <core/array.inl>
#include <core/pool_allocator.h>
#include <core/arena_allocator.h>

// @TODO: Implement this test case for all allocators using Fixtures?
// It's vital that we test realloc properly since array is implemented using it.
UTEST(allocator, default) {
    md_allocator_i* alloc = default_allocator;

    void* mem = md_alloc(alloc, 16);
    EXPECT_NE(mem, NULL);
    EXPECT_EQ(md_free(alloc, mem, 16), NULL);
    
    // This will use realloc, when it grows internally
    uint64_t* arr = NULL;
    for (uint64_t i = 0; i < 1000; ++i) {
        md_array_push(arr, i, alloc);
    }

    // We now know what we expect to find within arr, this will catch if realloc is not implemented properly
    for (uint64_t i = 0; i < 1000; ++i) {
        ASSERT_EQ(arr[i], i);
    }

    md_array_free(arr, alloc);
}

UTEST(allocator, default_temp) {
    md_allocator_i* alloc = default_temp_allocator;

    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
           md_alloc(alloc, i);
        }
    }
}

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

        for (int i = 0; i < ARRAY_SIZE(indices); ++i) {
            const int idx = indices[i];
            md_free(pool, items[idx], sizeof(uint64_t));
        }

        for (int i = 0; i < ARRAY_SIZE(indices); ++i) {
            const int idx = indices[i];
            items[idx] = md_alloc(pool, sizeof(uint64_t));
        }
    }

    md_array_free(items, default_allocator);
}

UTEST(allocator, arena) {
    md_allocator_i* arena = md_arena_allocator_create(default_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);

    for (int j = 0; j < 1000; ++j) {
        EXPECT_EQ((uint64_t)md_alloc(arena, 16) % 16, 0); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);   // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);   // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);   // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);   // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);   // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 2) % 2, 0);   // Expect to be aligned to 2 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);   // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 2) % 2, 0);   // Expect to be aligned to 2 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes

        EXPECT_NE(md_alloc(arena, 9000), NULL); // Exceeds page size, should still be good

        // Make sure we get some internal pages going.
        for (int i = 0; i < 20; ++i) {
            EXPECT_EQ((uint64_t)md_alloc(arena, 1024) % 16, 0);
        }

        md_arena_allocator_reset(arena);
    }

    md_arena_allocator_destroy(arena);
}