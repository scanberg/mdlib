#include "utest.h"

#include <core/md_handle.h>
#include <core/md_allocator.h>

#include <stdlib.h>
#include <string.h>

/* The pool hands out ids built from a slot index and a generation counter. The index is what makes
 * a lookup cheap; the generation is the entire reason to prefer a handle over a pointer, because it
 * is what lets a stale handle be recognised after its slot has been reused. So these tests are as
 * interested in the generation as in the index.
 *
 * An allocator that poisons what it hands back, rather than the heap allocator, because a pool's
 * counters must start at a known value whatever memory it is given. With fresh pages from the OS
 * the distinction is invisible - they arrive zeroed - and a pool built on recycled memory would
 * behave differently from one built on new memory, which is the sort of difference that only shows
 * up later and somewhere else. */
static void* poison_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)inst; (void)file; (void)line;
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    void* mem = realloc(ptr, new_size);
    if (mem && new_size > old_size) {
        memset((char*)mem + old_size, 0xCD, new_size - old_size);
    }
    return mem;
}

static md_allocator_i poison_allocator(void) {
    md_allocator_i alloc = {NULL, poison_realloc};
    return alloc;
}

#define POOL_COUNT 8

UTEST(handle, every_slot_starts_at_generation_one) {
    md_allocator_i alloc = poison_allocator();
    md_handle_pool_t pool = {0};
    md_handle_pool_init(&pool, POOL_COUNT, &alloc);

    /* The first handle drawn from a slot has been through the counter exactly once, so its
     * generation is 1 - for every slot, not just the first couple. */
    for (int i = 0; i < POOL_COUNT; ++i) {
        const uint32_t id = md_handle_pool_alloc_slot(&pool);
        ASSERT_NE(MD_HANDLE_INVALID_ID, id);
        EXPECT_EQ(1u, id >> MD_HANDLE_SLOT_SHIFT);
    }

    md_handle_pool_free(&pool, &alloc);
}

UTEST(handle, slots_are_distinct_and_in_range) {
    md_allocator_i alloc = poison_allocator();
    md_handle_pool_t pool = {0};
    md_handle_pool_init(&pool, POOL_COUNT, &alloc);

    int seen[POOL_COUNT + 1] = {0};
    for (int i = 0; i < POOL_COUNT; ++i) {
        const uint32_t id = md_handle_pool_alloc_slot(&pool);
        ASSERT_NE(MD_HANDLE_INVALID_ID, id);
        const int index = md_handle_index(id);
        /* Slot 0 is reserved so that an all-zero id can mean "no handle". */
        ASSERT_GT(index, 0);
        ASSERT_LE(index, POOL_COUNT);
        EXPECT_EQ(0, seen[index]);
        seen[index] = 1;
    }

    md_handle_pool_free(&pool, &alloc);
}

UTEST(handle, an_exhausted_pool_returns_the_invalid_id) {
    md_allocator_i alloc = poison_allocator();
    md_handle_pool_t pool = {0};
    md_handle_pool_init(&pool, POOL_COUNT, &alloc);

    for (int i = 0; i < POOL_COUNT; ++i) {
        ASSERT_NE(MD_HANDLE_INVALID_ID, md_handle_pool_alloc_slot(&pool));
    }
    /* Exhaustion is reported, not asserted on, and it keeps being reported. */
    EXPECT_EQ((uint32_t)MD_HANDLE_INVALID_ID, md_handle_pool_alloc_slot(&pool));
    EXPECT_EQ((uint32_t)MD_HANDLE_INVALID_ID, md_handle_pool_alloc_slot(&pool));

    md_handle_pool_free(&pool, &alloc);
}

UTEST(handle, a_reused_slot_gets_a_new_id) {
    md_allocator_i alloc = poison_allocator();
    md_handle_pool_t pool = {0};
    md_handle_pool_init(&pool, POOL_COUNT, &alloc);

    const uint32_t first = md_handle_pool_alloc_slot(&pool);
    ASSERT_NE(MD_HANDLE_INVALID_ID, first);
    md_handle_pool_free_slot(&pool, first);

    const uint32_t second = md_handle_pool_alloc_slot(&pool);
    ASSERT_NE(MD_HANDLE_INVALID_ID, second);

    /* Same storage, different identity. This is the property the whole scheme exists for: holding
     * on to 'first' after freeing it must not silently address whatever now lives in that slot. */
    EXPECT_EQ(md_handle_index(first), md_handle_index(second));
    EXPECT_NE(first, second);
    EXPECT_EQ((first >> MD_HANDLE_SLOT_SHIFT) + 1u, second >> MD_HANDLE_SLOT_SHIFT);

    md_handle_pool_free(&pool, &alloc);
}

UTEST(handle, freeing_returns_capacity) {
    md_allocator_i alloc = poison_allocator();
    md_handle_pool_t pool = {0};
    md_handle_pool_init(&pool, POOL_COUNT, &alloc);

    uint32_t ids[POOL_COUNT];
    for (int i = 0; i < POOL_COUNT; ++i) {
        ids[i] = md_handle_pool_alloc_slot(&pool);
        ASSERT_NE(MD_HANDLE_INVALID_ID, ids[i]);
    }
    ASSERT_EQ((uint32_t)MD_HANDLE_INVALID_ID, md_handle_pool_alloc_slot(&pool));

    for (int i = 0; i < POOL_COUNT; ++i) {
        md_handle_pool_free_slot(&pool, ids[i]);
    }
    /* Everything handed back is available again - the pool neither leaks slots nor overcounts. */
    for (int i = 0; i < POOL_COUNT; ++i) {
        EXPECT_NE(MD_HANDLE_INVALID_ID, md_handle_pool_alloc_slot(&pool));
    }
    EXPECT_EQ((uint32_t)MD_HANDLE_INVALID_ID, md_handle_pool_alloc_slot(&pool));

    md_handle_pool_free(&pool, &alloc);
}

UTEST(handle, free_clears_the_pool) {
    md_allocator_i alloc = poison_allocator();
    md_handle_pool_t pool = {0};
    md_handle_pool_init(&pool, POOL_COUNT, &alloc);
    md_handle_pool_free(&pool, &alloc);

    EXPECT_EQ(0, pool.size);
    EXPECT_EQ(0, pool.queue_top);
    EXPECT_TRUE(pool.free_queue == NULL);
    EXPECT_TRUE(pool.gen_counters == NULL);
}

#undef POOL_COUNT
