#ifndef _BIT_OP_H_
#define _BIT_OP_H_

#include <stdint.h>
#include <stdbool.h>

#if _MSC_VER && !__INTEL_COMPILER && !__clang__
#include <intrin.h>
#pragma intrinsic(_popcnt64)
#pragma intrinsic(_BitScanForward64)

inline uint64_t bit_count_mask(uint64_t mask) {
    return __popcnt64(mask);
}

inline bool bit_scan_forward_mask(uint64_t* bit_idx, uint64_t mask) {
    unsigned long i;
    bool result = (bool)_BitScanForward64(&i, mask);
    *bit_idx = i;
    return result;
}

#elif defined(__GNUC__)

inline uint64_t bit_count_mask(uint64_t mask) {
	return __builtin_popcountll(mask);
}

inline bool bit_scan_forward_mask(uint64_t* bit_idx, uint64_t mask) {
    uint64_t idx = __builtin_ffsll(mask);
    if (idx) {
        *bit_idx = idx - 1;
        return true;
    }
    return false;
}

#endif

void bit_set    (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);
void bit_clear  (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);
void bit_invert (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

void bit_or     (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_or_not (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_and    (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_and_not(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_xor    (uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);
void bit_not    (uint64_t* dst, const uint64_t* src, uint64_t bit_offset, uint64_t bit_count);

//void bit_copy   (uint64_t* dst, uint64_t dst_offset, const uint64_t* src, uint64_t src_offset, uint64_t bit_count);

// Counts the number of bits set
uint64_t bit_count(const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

// Test if single bit in field is set
bool bit_test(const uint64_t* bits, uint64_t idx);

bool find_next_bit_set(uint64_t* bit_idx, const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

#endif