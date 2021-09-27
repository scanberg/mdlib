#include "utest.h"

#include <core/md_intrinsics.h>

UTEST(intrinsics, test) {
    {
        static const uint32_t mask[5] = {0, 0xFFU, 0xFFFFU, 0xFFFFFFU, 0xFFFFFFFFU};

        EXPECT_EQ(clz32(mask[1]), 24);
        EXPECT_EQ(clz32(mask[2]), 16);
        EXPECT_EQ(clz32(mask[3]), 8);
        EXPECT_EQ(clz32(mask[4]), 0);

        EXPECT_EQ(popcnt32(mask[0]), 0);
        EXPECT_EQ(popcnt32(mask[1]), 8);
        EXPECT_EQ(popcnt32(mask[2]), 16);
        EXPECT_EQ(popcnt32(mask[3]), 24);
        EXPECT_EQ(popcnt32(mask[4]), 32);
    }

    {
        static const uint32_t mask[5] = {0, (1U << 8), (1U << 16), (1U << 24), (1U << 31)};

        EXPECT_EQ(bit_scan_forward32(mask[0]), 0);
        EXPECT_EQ(bit_scan_forward32(mask[1]), 9);
        EXPECT_EQ(bit_scan_forward32(mask[2]), 17);
        EXPECT_EQ(bit_scan_forward32(mask[3]), 25);
        EXPECT_EQ(bit_scan_forward32(mask[4]), 32);
    }

    EXPECT_EQ(next_power_of_two32(3), 4);
    EXPECT_EQ(next_power_of_two32(31), 32);
    EXPECT_EQ(next_power_of_two32(32), 32);
    EXPECT_EQ(next_power_of_two32(5), 8);
    EXPECT_EQ(next_power_of_two32(8000), 8192);

    EXPECT_EQ(next_power_of_two64(sizeof(void*)), sizeof(void*));



    uint32_t res = 0;
    res = find_first_zero_byte64(65280);
    res = find_first_zero_byte64(255);
    res = find_first_zero_byte64(65535);


    EXPECT_EQ(find_first_zero_byte32(65280),      0); // 1111111100000000
    EXPECT_EQ(find_first_zero_byte32(255),        1); // 0000000011111111
    EXPECT_EQ(find_first_zero_byte32(65535),      2); // 000000001111111111111111

}
