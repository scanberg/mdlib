#include "utest.h"
#include <core/md_allocator.h>
#include <core/md_bitfield.h>
#include <core/md_tracking_allocator.h>

#include <string.h>
#include <stdint.h>

UTEST(bitfield, test) {
    md_allocator_i* alloc = md_tracking_allocator_create(default_allocator);

    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);
    md_bitfield_set_range(&bf, 10000, 10100);
    EXPECT_EQ(md_bitfield_popcount(&bf), 100);

    md_bitfield_set_range(&bf, 500, 600);
    EXPECT_EQ(md_bitfield_popcount(&bf), 200);
    md_bitfield_clear(&bf);

    md_bitfield_set_bit(&bf, 1000);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 1000));
    EXPECT_EQ(md_bitfield_popcount(&bf), 1);

    md_bitfield_set_bit(&bf, 2000);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2000));
    EXPECT_EQ(md_bitfield_popcount(&bf), 2);

    md_bitfield_set_bit(&bf, 2001);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2001));

    md_bitfield_set_bit(&bf, 2005);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2005));

    md_bitfield_set_bit(&bf, 2018);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2018));

    md_bitfield_set_bit(&bf, 2022);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2022));

    md_bitfield_set_bit(&bf, 2089);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2089));

    md_bitfield_set_bit(&bf, 2100);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2100));

    md_bitfield_set_bit(&bf, 2101);
    EXPECT_TRUE(md_bitfield_test_bit(&bf, 2101));

    EXPECT_EQ(md_bitfield_popcount(&bf), 9);

    int64_t beg_idx = 0;
    int64_t end_idx = 3000;

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 1000);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2000);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2001);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2005);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2018);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2022);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2089);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2100);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx-1, 2101);

    beg_idx = md_bitfield_scan(&bf, beg_idx, end_idx);
    EXPECT_EQ(beg_idx, 0);

    md_bitfield_set_range(&bf, 2500, 2600);
    EXPECT_EQ(md_bitfield_popcount(&bf), 109);

    md_bitfield_clear_range(&bf, 2550, 2560);
    EXPECT_EQ(md_bitfield_popcount(&bf), 99);

    md_bitfield_clear(&bf);
    EXPECT_EQ(md_bitfield_popcount(&bf), 0);

    md_bitfield_t mask = {0};
    md_bitfield_init(&mask, default_allocator);
    md_bitfield_set_range(&mask, 100, 200);
    EXPECT_EQ(md_bitfield_popcount(&mask), 100);

    md_bitfield_or_inplace(&bf, &mask);
    EXPECT_EQ(md_bitfield_popcount(&bf), 100);

    for (int64_t i = 0; i < 100; ++i) {
        EXPECT_FALSE(md_bitfield_test_bit(&bf, i));
    }
    for (int64_t i = 100; i < 200; ++i) {
        EXPECT_TRUE(md_bitfield_test_bit(&bf, i));
    }
    for (int64_t i = 200; i < 1000; ++i) {
        EXPECT_FALSE(md_bitfield_test_bit(&bf, i));
    }

    md_bitfield_not(&bf, &mask, 0, 1200);
    EXPECT_EQ(md_bitfield_popcount(&bf), 1100);

    for (int64_t i = 0; i < 100; ++i) {
        EXPECT_TRUE(md_bitfield_test_bit(&bf, i));
    }
    for (int64_t i = 100; i < 200; ++i) {
        EXPECT_FALSE(md_bitfield_test_bit(&bf, i));
    }
    for (int64_t i = 200; i < 1200; ++i) {
        EXPECT_TRUE(md_bitfield_test_bit(&bf, i));
    }

    md_bitfield_clear(&bf);
    md_bitfield_set_range(&bf, 10000, 10001);
    md_bitfield_not_inplace(&bf, 0, 12000);
    int64_t count = md_bitfield_popcount(&bf);
    EXPECT_EQ(count, 12000-1);
    
    md_bitfield_free(&mask);
    md_bitfield_free(&bf);

    md_tracking_allocator_destroy(alloc);
}

UTEST(bitfield, serialization) {
    md_allocator_i* alloc = md_tracking_allocator_create(default_allocator);

    md_bitfield_t a = {0};
    md_bitfield_init(&a, alloc);

    md_bitfield_set_bit(&a, 100);
    md_bitfield_set_bit(&a, 110);
    md_bitfield_set_bit(&a, 510);
    md_bitfield_set_bit(&a, 600);
    md_bitfield_set_bit(&a, 700);
    md_bitfield_set_bit(&a, 11992);
    md_bitfield_set_range(&a, 1000, 2000);
    md_bitfield_set_bit(&a, 1 << 16);

    int64_t est_bytes = md_bitfield_serialize_size_in_bytes(&a);
    void* mem = md_alloc(alloc, est_bytes);
    int64_t real_bytes = md_bitfield_serialize(mem, &a);

    md_bitfield_t b = {0};
    md_bitfield_init(&b, alloc);

    EXPECT_TRUE(md_bitfield_deserialize(&b, mem, real_bytes));

    for (int64_t i = 0; i < 70000; ++i) {
        bool bit_a = md_bitfield_test_bit(&a, i);
        bool bit_b = md_bitfield_test_bit(&b, i);
        ASSERT_TRUE(bit_a == bit_b);
    }

    md_bitfield_t c = {0};
    md_bitfield_init(&c, alloc);
    md_bitfield_set_bit(&c, 0);

    real_bytes = md_bitfield_serialize(mem, &c);

    md_bitfield_clear(&b);
    EXPECT_TRUE(md_bitfield_deserialize(&b, mem, real_bytes));
    for (int64_t i = 0; i < 100; ++i) {
        bool bit_c = md_bitfield_test_bit(&c, i);
        bool bit_b = md_bitfield_test_bit(&b, i);
        ASSERT_TRUE(bit_c == bit_b);
    }

    md_free(alloc, mem, est_bytes);
    md_bitfield_free(&a);
    md_bitfield_free(&b);
    md_bitfield_free(&c);
}