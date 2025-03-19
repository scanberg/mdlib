#include "utest.h"
#include <core/md_allocator.h>
#include <core/md_bitfield.h>
#include <core/md_tracking_allocator.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_base64.h>

#include <string.h>
#include <stdint.h>

static void set_bitfield(md_bitfield_t* bf, const char* str) {
    md_bitfield_clear(bf);
    const size_t len = strlen(str);
    for (size_t i = 0; i < len; ++i) {
        if (str[i] == '1') {
            md_bitfield_set_bit(bf, i);
        }
    }
}

static bool cmp_bitfield(const md_bitfield_t* bf, const char* str) {
    const uint32_t beg_idx = 0;
    const uint32_t end_idx = bf->end_bit;
    const size_t len = strlen(str);

    for (uint32_t i = beg_idx; i < end_idx; ++i) {
        const bool  bf_bit = md_bitfield_test_bit(bf, i);
        const bool ref_bit = i < len ? str[i] == '1' : false;
        if (bf_bit != ref_bit) {
            return false;
        }
    }
    return true;
}

UTEST(bitfield, bit_op) {
    md_bitfield_t a = {0};
    md_bitfield_t b = {0};
    md_bitfield_t c = {0};

    md_bitfield_init(&a, md_get_heap_allocator());
    md_bitfield_init(&b, md_get_heap_allocator());
    md_bitfield_init(&c, md_get_heap_allocator());

    set_bitfield(&a, "10001110101");
    set_bitfield(&b, "01101011101");

    md_bitfield_xor(&c, &a, &a);
    EXPECT_EQ(0, md_bitfield_popcount(&c));

    md_bitfield_xor(&c, &b, &b);
    EXPECT_EQ(0, md_bitfield_popcount(&c));

    md_bitfield_not(&c, &a, 0, a.end_bit);
    EXPECT_TRUE(cmp_bitfield(&c, "01110001010"));

    md_bitfield_not(&c, &b, 0, b.end_bit);
    EXPECT_TRUE(cmp_bitfield(&c, "10010100010"));

    md_bitfield_andnot(&c, &a, &b);
    EXPECT_TRUE(cmp_bitfield(&c, "10000100000"));
    
    md_bitfield_or(&c, &a, &b);
    EXPECT_TRUE(cmp_bitfield(&c, "11101111101"));

    md_bitfield_and(&c, &a, &b);
    EXPECT_TRUE(cmp_bitfield(&c, "00001010101"));

    md_bitfield_free(&a);
    md_bitfield_free(&b);
    md_bitfield_free(&c);
}

UTEST(bitfield, realloc) {
    md_allocator_i* alloc = md_tracking_allocator_create(md_get_heap_allocator());

    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);
    md_bitfield_set_bit(&bf, 0);
    EXPECT_EQ(md_bitfield_popcount(&bf), 1);

    md_bitfield_set_bit(&bf, 10001);
    EXPECT_EQ(md_bitfield_popcount(&bf), 2);

    md_bitfield_free(&bf);

    md_tracking_allocator_destroy(alloc);
}

UTEST(bitfield, general) {
    md_allocator_i* alloc = md_tracking_allocator_create(md_get_heap_allocator());

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
    md_bitfield_init(&mask, md_get_heap_allocator());
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
    {
        size_t count = md_bitfield_popcount(&bf);
        EXPECT_EQ(0, count);
    }

    md_bitfield_set_range(&bf, 10000, 10001);
    {
        size_t count = md_bitfield_popcount(&bf);
        EXPECT_EQ(1, count);
    }

    md_bitfield_not_inplace(&bf, 0, 12000);
    {
        size_t count = md_bitfield_popcount(&bf);
        EXPECT_EQ(count, 12000-1);
    }
    
    md_bitfield_free(&mask);
    md_bitfield_free(&bf);

    md_tracking_allocator_destroy(alloc);
}

UTEST(bitfield, serialization) {
    md_allocator_i* alloc = md_tracking_allocator_create(md_get_heap_allocator());

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

    EXPECT_EQ(1007, md_bitfield_popcount(&a));

    size_t est_bytes = md_bitfield_serialize_size_in_bytes(&a);
    void* mem = md_alloc(alloc, est_bytes);
    size_t real_bytes = md_bitfield_serialize(mem, &a);

    MD_LOG_INFO("Estimated serialization bytes for a: %i, actual bytes: %i", (int)est_bytes, (int)real_bytes);

    md_bitfield_t b = {0};
    md_bitfield_init(&b, alloc);

    bool result = md_bitfield_deserialize(&b, mem, real_bytes);
    MD_LOG_INFO("Deserialization of b: %s", result ? "true" : "false");
    EXPECT_TRUE(result);

    EXPECT_EQ(a.beg_bit, b.beg_bit);
    EXPECT_EQ(a.end_bit, b.end_bit);

    for (size_t i = 0; i < 70000; ++i) {
        bool bit_a = md_bitfield_test_bit(&a, i);
        bool bit_b = md_bitfield_test_bit(&b, i);
        if (bit_a != bit_b) {
            while(0){};
        }
        ASSERT_TRUE(bit_a == bit_b);
    }

    md_bitfield_t c = {0};
    md_bitfield_init(&c, alloc);
    md_bitfield_set_bit(&c, 0);

    real_bytes = md_bitfield_serialize(mem, &c);
    MD_LOG_INFO("Serialized bytes bytes for c: %i", (int)est_bytes, (int)real_bytes);


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

UTEST(bitfield, deserialize_base64) {
    str_t str = STR_LIT("P38AhACFAIcAiACJAIoAiwCMAI0AjgCPgJCAkYCSgJOAH5SAlYCWgJeAmICZgJqAm4CcgJ2AnoCfgKCAoYCigKOAH6SApYCmgKeAqICpgKqAq4CsgK2AroCvgLCAsYCygLOAH7SAtYC2gLeAuIC5gLqAu4C8gL2AvoC/gMCAwYDCgMOAH8SAxYDGgMeAyIDJgMqAy4DMgM2AzoDPgNCA0YDSgNOAH9SA1YDWgNeA2IDZgNqA24DcgN2A3oDfgOCA4YDigOOAH+SA5YDmgOeA6IDpgOqA64DsgO2A7oDvgPCA8YDygPOAH/SA9YD2gPeA+ID5gPqA+4D8gP0ANgE3ATgBOQE6ATsBAADgMwAIJzsO/v/vzgYC4DNE4AcABQSACIBz/+ADAAR/NzMcBOAGJgOAACII4AQSAYA/YA4EwP4RwQNgCeAQAAuhEQQg9z8m9+tvNPQgZgUfDzPyiGDgEC0KAABEgwAAue+7+w8gLOAQAAMbAiKCICsBAO+gIQLP//5ACQD3QATgBgAJv+GOIAiIeBEBBuACwQGg++ALKsAAINvhBmMG8AcAyD31F0BkAOPAJeBHACDLA79B9EDgG1YAOyCICNeBs7v8ffR/0CAwAL8gA+A7AAAfIVLgHABgbADfgs4AwGAMoAAAD2A+AICADYIZAOCADIAAAAdgHwNAF++lYAjgFQCAMoAAIZNAACD4IADgCmWgACAcIACAZeAEAADwIBYi2wUAAPwBAODgBBngEQABPAwg26HvgMvgBAAA+CBQ4AaYgAAgZeABAKD+gBDgBDKAAOAAmEAAAPxgQOEMhgQAAAAAAA==");
    md_allocator_i* alloc = md_tracking_allocator_create(md_get_heap_allocator());

    const int cap = md_base64_decode_size_in_bytes((int)str.len);
    void* mem = md_alloc(alloc, cap);
    int   len = md_base64_decode(mem, str.ptr, (int)str.len);

    md_bitfield_t bf = md_bitfield_create(alloc);
    EXPECT_TRUE(md_bitfield_deserialize(&bf, mem, len));
    md_bitfield_free(&bf);

    md_free(alloc, mem, cap);

    md_tracking_allocator_destroy(alloc);
    
}

UTEST(bitfield, iterator) {
    md_allocator_i* alloc = md_tracking_allocator_create(md_get_heap_allocator());

    md_bitfield_t bf = md_bitfield_create(alloc);

    md_bitfield_set_bit(&bf, 60);
    md_bitfield_set_bit(&bf, 64);
    md_bitfield_set_bit(&bf, 1000);
    
    md_bitfield_iter_t it = md_bitfield_iter_create(&bf);

    size_t pop = md_bitfield_popcount(&bf);
    EXPECT_EQ(3, pop);
    
    EXPECT_TRUE(md_bitfield_iter_next(&it));
    EXPECT_EQ(md_bitfield_iter_idx(&it), 60);

    EXPECT_TRUE(md_bitfield_iter_next(&it));
    EXPECT_EQ(md_bitfield_iter_idx(&it), 64);

    EXPECT_TRUE(md_bitfield_iter_next(&it));
    EXPECT_EQ(md_bitfield_iter_idx(&it), 1000);

    EXPECT_FALSE(md_bitfield_iter_next(&it));

    md_bitfield_free(&bf);

    md_tracking_allocator_destroy(alloc);
}