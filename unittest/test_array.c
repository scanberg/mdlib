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
