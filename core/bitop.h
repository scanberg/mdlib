#ifndef _BIT_OP_H_
#define _BIT_OP_H_

#include <stdint.h>
#include <stdbool.h>

#if _MSC_VER && !__INTEL_COMPILER
#include <intrin.h>
#pragma intrinsic(_BitScanForward)
#pragma intrinsic(_BitScanForward64)
#endif

inline uint64_t bit_count_mask(uint64_t mask) {
#if 0
    x = x - ((x >> 1) & 0x5555555555555555);
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
    x = ((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F);
    return (x * (0x0101010101010101)) >> 56;
#endif
    return __popcnt64(mask);
}

inline bool bit_scan_forward_mask(uint64_t* bit_idx, uint64_t mask) {
    unsigned long i;
    bool result = (bool)_BitScanForward64(&i, mask);
    *bit_idx = i;
    return result;
}

void bit_set    (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);
void bit_clear  (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);
void bit_invert (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

void bit_or     (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_or_not (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_and    (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_and_not(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_xor    (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_not    (uint64_t* dst, const uint64_t* src, uint64_t bit_offset, uint64_t bit_count);

// Counts the number of bits set
uint64_t bit_count(const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

bool find_next_bit_set(uint64_t* bit_idx, const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

#endif