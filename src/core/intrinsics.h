#ifndef _MD_INTRINSICS_H_
#define _MD_INTRINSICS_H_

#include <stdint.h>
#include <stdbool.h>
#include "common.h"

#if (MD_COMPILER_CLANG || MD_COMPILER_GCC)
// count leading zeros
static inline uint32_t clz(uint32_t v) {
    ASSERT(v != 0); // 0 is undefined behaviour
    return __builtin_clz(v);
}

// count leading zeros 64-bit
static inline uint64_t clz64(uint64_t v) {
    ASSERT(v != 0); // 0 is undefined behaviour
    return __builtin_clzll(v);
}

// Population count: counts the number of bits set
static inline uint32_t popcnt(uint32_t v) {
    return __builtin_popcount(v);
}

// Population count: counts the number of bits set
static inline uint64_t popcnt64(uint64_t v) {
    return __builtin_popcountll(v);
}

// Scans for the first bit set, from least significant to most significant bit,
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bit_scan_forward(uint32_t v) {
    return __builtin_ffs(v);
}

// Scans for the first bit set, from least significant to most significant bit,
// indexing starts at 1, returns 0 if no bit is set 
static inline uint64_t bit_scan_forward64(uint64_t v) {
    return __builtin_ffsll(v);
}

#elif MD_COMPILER_MSVC
#include <intrin.h>

// count leading zeros
static inline uint32_t clz(uint32_t v) {
    ASSERT(v != 0); // 0 is undefined behaviour
    return __lzcnt(v);
}

// count leading zeros 64-bit
static inline uint64_t clz64(uint64_t v) {
    return __lzcnt64(v);
}

// count bits
static inline uint32_t popcnt(uint32_t v) {
    return __popcnt(v);
}

// count bits 64-bit
static inline uint64_t popcnt64(uint64_t v) {
    return __popcnt64(v);
}

// Scans for the first bit set from least significant bit (LSB) to most significant bit (MSB)
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bit_scan_forward(uint32_t mask) {
    if (mask == 0) return 0;
    unsigned long idx;
    _BitScanForward(&idx, mask);
    return idx + 1;
}

// Searches for the first bit set from most significant bit (MSB) to least significant bit (LSB)
// index starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bit_scan_reverse(uint32_t mask) {
    if (mask == 0) return 0;
    unsigned long idx;
    _BitScanReverse(&idx, mask);
    return idx + 1;
}

// Searches for the first bit set from most significant bit (MSB) to least significant bit (LSB)
// index starts at 1, returns 0 if no bit is set
static inline uint64_t bit_scan_reverse64(uint64_t mask) {
    if (mask == 0) return 0;
    unsigned long idx;
    _BitScanReverse64(&idx, mask);
    return (uint64_t)idx + 1;
} 

// Scans for the first bit set, from least significant to most significant bit,
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint64_t bit_scan_forward64(uint64_t mask) {
    if (mask == 0) return 0;
    unsigned long idx;
    _BitScanForward64(&idx, mask);
    return (uint64_t)idx + 1;
}

#endif

// find first zero Byte
static inline uint32_t find_first_zero_byte(uint32_t x) {
    uint32_t y = (x & 0x7F7F7F7F) + 0x7F7F7F7F;
    y = ~(y | x | 0x7F7F7F7F);
    return (31U - clz(y)) >> 3;
}

// find first zero Byte 64-bit
static inline uint64_t find_first_zero_byte64(uint64_t x) {
    //uint64_t y = (x & 0x7F7F7F7F7F7F7F7F) + 0x7F7F7F7F7F7F7F7F;
    //y = ~(y | x | 0x7F7F7F7F7F7F7F7F);
    //return (63LLU - clz(y)) >> 3LLU;
    
    uint64_t y = (x - 0x0101010101010101) & ~x & 0x8080808080808080;
    return clz64(y) >> 3;
}

static inline uint32_t next_power_of_two(uint32_t x) {
    if (x == 0) return 0;
    if (x == 1) return 1; // avoid clz(0)
    return 1UL << (sizeof(uint32_t) * 8 - clz(x-1));
}

static inline uint64_t next_power_of_two64 (uint64_t x) {
    if (x == 0) return 0;
    if (x == 1) return 1; // avoid clz(0)
    return 1ULL << (sizeof(uint64_t) * 8 - clz64(x-1));
}

#endif

