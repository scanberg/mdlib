#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_array.h>

UTEST(array, test) {
    uint64_t* arr = 0;
    for (uint64_t i = 0; i < 1000; ++i) {
        md_array_push(arr, i, md_heap_allocator);
    }

    EXPECT_EQ(md_array_size(arr), 1000);

    for (uint64_t i = 0; i < 1000; ++i) {
        ASSERT_EQ(arr[i], i);
    }

    for (uint64_t i = 0; i < 1000; ++i) {
        ASSERT_EQ(md_array_pop(arr), 999-i);
    }
}