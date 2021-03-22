#include "utest.h"
#include <md_allocator.h>
#include <core/bitop.h>

#include <string.h>
#include <stdint.h>

static uint64_t set_bits(uint64_t* bits, const char* bit_str) {    
    const uint64_t num_bits = strlen(bit_str);
    for (uint64_t i = 0; i < num_bits; i++) {
        const uint64_t blk_idx = i / 64;
        const uint64_t bit_mask = (1LLU << i);
        if (bit_str[i] != '0')
            bits[blk_idx] |= bit_mask;
        else
            bits[blk_idx] &= ~bit_mask;
    }
    return num_bits;
}

static void print_bits(uint64_t* bits, uint64_t num_bits) {
    for (uint64_t i = 0; i < num_bits; ++i) {
        const uint64_t blk_idx = i / 64;
        const uint64_t bit_mask = (1LLU << i);
        printf("%i", bits[blk_idx] & bit_mask ? 1 : 0);
    }
}

UTEST(bitop, test) {
    uint64_t bf[5][2] = {0};
    uint64_t num_bits = 0;
    bool res = false;
    
    num_bits = set_bits(bf[0], "00100000000000000000000001000000000000000000100000000000000000000001");
               set_bits(bf[1], "00000000000100000000010001000000000000000000000000000000010000000001");

    EXPECT_TRUE (bit_cmp(bf[0], bf[0], 0, num_bits));
    EXPECT_FALSE(bit_cmp(bf[0], bf[1], 0, num_bits));

    EXPECT_EQ(bit_count(bf[0], 0, num_bits), 4);
    EXPECT_EQ(bit_count(bf[1], 0, num_bits), 5);

    // SET
    bit_set(bf[2], 0, num_bits);
    set_bits(bf[3], "11111111111111111111111111111111111111111111111111111111111111111111");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));

    // CLEAR
    bit_clear(bf[2], 0, num_bits);
    set_bits(bf[3], "00000000000000000000000000000000000000000000000000000000000000000000");
    res = bit_cmp(bf[2], bf[3], 0, num_bits);
    EXPECT_TRUE(res);
    if (!res) {
        printf("Got:\n");
        print_bits(bf[2], num_bits);
        printf("\nExpected:\n");
        print_bits(bf[3], num_bits);
        printf("\n");
    }

    // NOT
    bit_not(bf[2], bf[0], 0, num_bits);
    set_bits(bf[3], "11011111111111111111111110111111111111111111011111111111111111111110");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));

    // OR
    bit_or(bf[2], bf[0], bf[1], 0, num_bits);
    set_bits(bf[3], "00100000000100000000010001000000000000000000100000000000010000000001");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));

    // OR NOT
    bit_or_not(bf[2], bf[0], bf[1], 0, num_bits);
    set_bits(bf[3], "11111111111011111111101111111111111111111111111111111111101111111111");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));
    
    // AND
    bit_and(bf[2], bf[0], bf[1], 0, num_bits);
    set_bits(bf[3], "00000000000000000000000001000000000000000000000000000000000000000001");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));

    // AND NOT
    bit_and_not(bf[2], bf[0], bf[1], 0, num_bits);
    set_bits(bf[3], "00100000000000000000000000000000000000000000100000000000000000000000");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));

    // XOR
    bit_xor(bf[2], bf[0], bf[1], 0, num_bits);
    set_bits(bf[3], "00100000000100000000010000000000000000000000100000000000010000000000");
    EXPECT_TRUE(bit_cmp(bf[2], bf[3], 0, num_bits));

    // SCAN
    EXPECT_EQ(bit_scan(bf[0], 0, num_bits), 3);
    EXPECT_EQ(bit_scan(bf[1], 0, num_bits), 12);
}