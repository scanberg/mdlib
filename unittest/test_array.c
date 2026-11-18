#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_array.h>

UTEST(array, test) {
    uint64_t* arr = 0;
    for (uint64_t i = 0; i < 10000; ++i) {
        md_array_push(arr, i, md_get_heap_allocator());
    }

    EXPECT_EQ(md_array_size(arr), 10000);

    for (uint64_t i = 0; i < 10000; ++i) {
        ASSERT_EQ(arr[i], i);
    }

    for (uint64_t i = 0; i < 10000; ++i) {
        ASSERT_EQ(md_array_back(arr), 9999-i);
        md_array_pop(arr);
    }

    md_array_free(arr, md_get_heap_allocator());

    for (uint64_t i = 0; i < 10000; ++i) {
        md_array_push(arr, i, md_get_heap_allocator());
    }

    EXPECT_EQ(md_array_size(arr), 10000);

    for (uint64_t i = 0; i < 10000; ++i) {
        ASSERT_EQ(arr[i], i);
    }

    for (uint64_t i = 0; i < 10000; ++i) {
        ASSERT_EQ(md_array_back(arr), 9999-i);
        md_array_pop(arr);
    }

    md_array_free(arr, md_get_heap_allocator());
}

MD_ARRAY_RAW_DEFINE_TYPE(test_u32_array_t, uint32_t);

// The asserting operations must still do their work where MD_ASSERT expands to nothing (release builds).
UTEST(array, typed_push_free) {
    test_u32_array_t arr;
    MD_ARRAY_INIT(arr, md_get_heap_allocator());
    for (uint32_t i = 0; i < 1000; ++i) MD_ARRAY_PUSH(arr, i);
    ASSERT_EQ(1000, (int)MD_ARRAY_SIZE(arr));
    for (uint32_t i = 0; i < 1000; ++i) ASSERT_EQ(i, arr.data[i]);
    MD_ARRAY_ENSURE(arr, 5000);
    EXPECT_GE(MD_ARRAY_CAPACITY(arr), (size_t)5000);
    MD_ARRAY_FREE(arr);
    EXPECT_TRUE(arr.data == NULL);
    EXPECT_EQ(0, (int)arr.size);
}

UTEST(array, typed_insert_remove_index_of) {
    test_u32_array_t arr;
    MD_ARRAY_INIT(arr, md_get_heap_allocator());
    for (uint32_t i = 0; i < 5; ++i) ASSERT_TRUE(MD_ARRAY_TRY_PUSH(arr, i));    // 0 1 2 3 4

    uint32_t v = 10;
    ASSERT_TRUE(MD_ARRAY_TRY_INSERT(arr, 0, v));                                  // 10 0 1 2 3 4
    v = 11;
    ASSERT_TRUE(MD_ARRAY_TRY_INSERT(arr, 3, v));                                  // 10 0 1 11 2 3 4
    v = 12;
    ASSERT_TRUE(MD_ARRAY_TRY_INSERT(arr, arr.size, v));                           // 10 0 1 11 2 3 4 12
    const uint32_t after_insert[] = {10, 0, 1, 11, 2, 3, 4, 12};
    ASSERT_EQ(8, (int)arr.size);
    for (int i = 0; i < 8; ++i) EXPECT_EQ(after_insert[i], arr.data[i]);

    v = 11;
    EXPECT_EQ(3, (int)MD_ARRAY_INDEX_OF(arr, v));
    v = 99;
    EXPECT_EQ((int)arr.size, (int)MD_ARRAY_INDEX_OF(arr, v));
    EXPECT_FALSE(MD_ARRAY_REMOVE_ITEM(arr, v));

    MD_ARRAY_REMOVE(arr, 0);                                                      // 0 1 11 2 3 4 12
    v = 11;
    EXPECT_TRUE(MD_ARRAY_REMOVE_ITEM(arr, v));                                    // 0 1 2 3 4 12
    MD_ARRAY_REMOVE(arr, arr.size - 1);                                           // 0 1 2 3 4
    ASSERT_EQ(5, (int)arr.size);
    for (uint32_t i = 0; i < 5; ++i) EXPECT_EQ(i, arr.data[i]);

    MD_ARRAY_FREE(arr);
}
