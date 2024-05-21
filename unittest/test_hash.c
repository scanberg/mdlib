#include "utest.h"

#include <core/md_hash.h>
#include <core/md_allocator.h>

UTEST(hash, u32) {
    md_hashmap32_t map = {.allocator = md_get_heap_allocator() };
    
    for (uint32_t i = 0; i < 100; ++i) {
        md_hashmap_add(&map, i, i * 2);
    }

    for (uint32_t i = 0; i < 100; ++i) {
        uint32_t *val = md_hashmap_get(&map, i);
        EXPECT_TRUE(val);
        if (val) {
            EXPECT_EQ(2 * i, *val);
        }
    }
    md_hashmap_free(&map);
}


UTEST(hash, u64) {
    md_hashmap64_t map = {.allocator = md_get_heap_allocator() };

    for (uint64_t i = 0; i < 100; ++i) {
        md_hashmap_add(&map, i, i * 2);
    }

    for (uint64_t i = 0; i < 100; ++i) {
        uint64_t *val = md_hashmap_get(&map, i);
        EXPECT_TRUE(val);
        if (val) {
            EXPECT_EQ(2 * i, *val);
        }
    }
    md_hashmap_free(&map);
}

UTEST(hash, custom_type) {
    typedef struct {
        uint32_t a;
        uint32_t b;
    } value_t;

    typedef struct MD_HASHMAP_T(value_t) map_t;
    map_t map = {.allocator = md_get_heap_allocator() };

    for (uint64_t i = 0; i < 100; ++i) {
        value_t val = {1 * i, 2 * i};
        md_hashmap_add(&map, i, val);
    }

    for (uint64_t i = 0; i < 100; ++i) {
        value_t *val = md_hashmap_get(&map, i);
        EXPECT_TRUE(val);
        if (val) {
            EXPECT_EQ(1 * i, val->a);
            EXPECT_EQ(2 * i, val->b);
        }
    }
    md_hashmap_free(&map);
}
