#include "utest.h"

#include <core/md_intrinsics.h>

UTEST(intrinsics, test) {
    {
        static const uint32_t mask[5] = {0, 0xFFU, 0xFFFFU, 0xFFFFFFU, 0xFFFFFFFFU};

        EXPECT_EQ(clz(mask[1]), 24);
        EXPECT_EQ(clz(mask[2]), 16);
        EXPECT_EQ(clz(mask[3]), 8);
        EXPECT_EQ(clz(mask[4]), 0);

        EXPECT_EQ(popcnt(mask[0]), 0);
        EXPECT_EQ(popcnt(mask[1]), 8);
        EXPECT_EQ(popcnt(mask[2]), 16);
        EXPECT_EQ(popcnt(mask[3]), 24);
        EXPECT_EQ(popcnt(mask[4]), 32);
    }

    {
        static const uint32_t mask[5] = {0, (1U << 8), (1U << 16), (1U << 24), (1U << 31)};

        EXPECT_EQ(bit_scan_forward(mask[0]), 0);
        EXPECT_EQ(bit_scan_forward(mask[1]), 9);
        EXPECT_EQ(bit_scan_forward(mask[2]), 17);
        EXPECT_EQ(bit_scan_forward(mask[3]), 25);
        EXPECT_EQ(bit_scan_forward(mask[4]), 32);
    }

    EXPECT_EQ(next_power_of_two(3), 4);
    EXPECT_EQ(next_power_of_two(31), 32);
    EXPECT_EQ(next_power_of_two(32), 32);
    EXPECT_EQ(next_power_of_two(5), 8);
    EXPECT_EQ(next_power_of_two(8000), 8192);

    EXPECT_EQ(next_power_of_two64(sizeof(void*)), sizeof(void*));

    // @TODO: Implement more tests for example find_first_zero_byte
}