#ifndef _BIT_OP_H_
#define _BIT_OP_H_

#include <stdint.h>
#include <stdbool.h>

void bit_set    (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);
void bit_clear  (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

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

// Test if bitfields are equivalent
bool bit_cmp(const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count);

// Finds the next bit set within the bitarray, if no bit is found -1 is returned (i.e. UINT64_MAX = 0xffffffffffffffff)
//uint64_t find_next_bit_set(const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count);

#endif