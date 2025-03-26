#pragma once

#include <core/md_common.h>

#include <stdint.h>

#if (MD_COMPILER_CLANG || MD_COMPILER_GCC)
#include <x86intrin.h>

#define popcnt32 __builtin_popcount
#define popcnt64 __builtin_popcountll

#define ctz32 __builtin_ctz
#define ctz64 __builtin_ctzll

#define bsf32 __builtin_ffs
#define bsf64 __builtin_ffsll

// count leading zeros
static inline uint32_t clz32(uint32_t x) {
    return x ? __builtin_clz(x) : 32;
}

// count leading zeros 64-bit
static inline uint64_t clz64(uint64_t x) {
    return x ? __builtin_clzll(x) : 32;
}

static inline uint32_t bsr32(uint32_t x) {
    return x ? 32 - __builtin_clz(x) : 0;
}

static inline uint64_t bsr64(uint64_t x) {
    return x ? 64 - __builtin_clzll(x) : 0;
}

#elif MD_COMPILER_MSVC
#include <intrin.h>

#define popcnt32 __popcnt
#define popcnt64 __popcnt64
#define clz32 __lzcnt
#define clz64 __lzcnt64
#define ctz32 _tzcnt_u32
#define ctz64 _tzcnt_u64

// Scans for the first bit set from least significant bit (LSB) to most significant bit (MSB)
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bsf32(uint32_t mask) {
    unsigned long idx;
    return _BitScanForward(&idx, mask) ? idx + 1 : 0;
}

// Searches for the first bit set from most significant bit (MSB) to least significant bit (LSB)
// index starts at 1, returns 0 if no bit is set (posix convention)
static inline uint32_t bsr32(uint32_t mask) {
    unsigned long idx;
    return _BitScanReverse(&idx, mask) ? idx + 1 : 0;
}

// Scans for the first bit set, from least significant to most significant bit,
// indexing starts at 1, returns 0 if no bit is set (posix convention)
static inline uint64_t bsf64(uint64_t mask) {
    unsigned long idx;
    return _BitScanForward64(&idx, mask) ? (uint64_t)idx + 1 : 0;
}

// Searches for the first bit set from most significant bit (MSB) to least significant bit (LSB)
// index starts at 1, returns 0 if no bit is set (posix convention)
static inline uint64_t bsr64(uint64_t mask) {
    unsigned long idx;
    return _BitScanReverse64(&idx, mask) ? (uint64_t)idx + 1 : 0;
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
    // avoid clz(0)
    return (x < 2) ? x : 1U << (32 - clz32(x-1));
}

static inline uint64_t next_power_of_two64 (uint64_t x) {
    // avoid clz(0)
    return (x < 2) ? x : 1ULL << (64 - clz64(x-1));
}

static inline uint32_t swap_endian32(uint32_t x) {
	return (x >> 24) | ((x >> 8) & 0x0000FF00) | ((x << 8) & 0x00FF0000) | (x << 24);
}

static inline uint64_t swap_endian64(uint64_t x) {
	return (x >> 56) | ((x >> 40) & 0x000000000000FF00) | ((x >> 24) & 0x0000000000FF0000) | ((x >> 8) & 0x00000000FF000000) |
	       ((x << 8) & 0x000000FF00000000) | ((x << 24) & 0x0000FF0000000000) | ((x << 40) & 0x00FF000000000000) | (x << 56);
}
