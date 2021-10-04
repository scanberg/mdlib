#ifndef _MD_INTRINSICS_H_
#define _MD_INTRINSICS_H_

#include "md_common.h"
#include "md_compiler.h"

#include <stdint.h>
#include <stdbool.h>

#if (MD_COMPILER_CLANG || MD_COMPILER_GCC)
#include <x86intrin.h>

// Population count: counts the number of bits set
static inline uint32_t popcnt32(uint32_t x) {
    return __builtin_popcount(x);
}

// Population count: counts the number of bits set
static inline uint64_t popcnt64(uint64_t x) {
    return __builtin_popcountll(x);
}

// count leading zeros
static inline uint32_t clz32(uint32_t x) {
    ASSERT(x != 0); // 0 is undefined behaviour
    return __builtin_clz(x);
}

// count leading zeros 64-bit
static inline uint64_t clz64(uint64_t x) {
    ASSERT(x != 0); // 0 is undefined behaviour
    return __builtin_clzll(x);
}

static inline uint32_t ctz32(uint32_t x) {
    return __builtin_ctz(x);
}

static inline uint64_t ctz64(uint64_t x) {
    return __builtin_ctzll(x);
}

// Scans for the first bit set, from least significant to most significant bit,
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bit_scan_forward32(uint32_t x) {
    return __builtin_ffs(x);
}

// Scans for the first bit set, from least significant to most significant bit,
// indexing starts at 1, returns 0 if no bit is set 
static inline uint64_t bit_scan_forward64(uint64_t x) {
    return __builtin_ffsll(x);
}

#elif MD_COMPILER_MSVC
#include <intrin.h>

// count bits
static inline uint32_t popcnt32(uint32_t x) {
    return __popcnt(x);
}

// count bits 64-bit
static inline uint64_t popcnt64(uint64_t x) {
    return __popcnt64(x);
}

// count leading zeros
static inline uint32_t clz32(uint32_t x) {
    ASSERT(x != 0); // 0 is undefined behaviour
    return __lzcnt(x);
}

// count leading zeros 64-bit
static inline uint64_t clz64(uint64_t x) {
    return __lzcnt64(x);
}

static inline uint32_t ctz32(uint32_t x) {
    return popcnt32(~x & (x-1U));
}

static inline uint64_t ctz64(uint64_t x) {
    return popcnt64(~x & (x-1LLU));
}

// Scans for the first bit set from least significant bit (LSB) to most significant bit (MSB)
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bit_scan_forward32(uint32_t mask) {
    if (mask == 0) return 0;
    unsigned long idx;
    _BitScanForward(&idx, mask);
    return idx + 1;
}

// Searches for the first bit set from most significant bit (MSB) to least significant bit (LSB)
// index starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bit_scan_reverse32(uint32_t mask) {
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
static inline uint32_t find_first_zero_byte32(uint32_t x) {
    uint32_t y = (x - 0x01010101) & ~x & 0x80808080;
    return ctz32(y) >> 3;
}

// find first zero Byte 64-bit
static inline uint64_t find_first_zero_byte64(uint64_t x) {
    uint64_t y = (x - 0x0101010101010101) & ~x & 0x8080808080808080;
    return ctz64(y) >> 3;
}

static inline uint32_t next_power_of_two32(uint32_t x) {
    if (x < 2) return x;  // avoid clz(0)
    return 1UL << (sizeof(uint32_t) * 8 - clz32(x-1));
}

static inline uint64_t next_power_of_two64 (uint64_t x) {
    if (x < 2) return x;  // avoid clz(0)
    return 1ULL << (sizeof(uint64_t) * 8 - clz64(x-1));
}

#endif

