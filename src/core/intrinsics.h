#ifndef _MD_INTRINSICS_H_
#define _MD_INTRINSICS_H_

#include <stdint.h>
#include <stdbool.h>
#include "compiler.h"

#if (MD_COMPILER_CLANG || MD_COMPILER_GCC)
// count leading zeros
inline uint32_t clz(uint32_t v) {
    return __builtin_clz(v);
}

// count leading zeros 64-bit
inline uint64_t clz64(uint64_t v) {
    return __builtin_clzll(v);
}

inline uint64_t popcnt(uint64_t v) {
    return __builtin_popcount(v);
}

inline uint64_t popcnt64(uint64_t v) {
    return __builtin_popcountll(v);
}

#elif MD_COMPILER_MSVC
#include <intrin.h>

// count leading zeros
inline uint32_t clz(uint32_t v) {
    return __lzcnt(v);
}

// count leading zeros 64-bit
inline uint64_t clz64(uint64_t v) {
    return __lzcnt64(v);
}

// count bits
inline uint64_t popcnt(uint64_t v) {
    return __popcnt(v);
}

// count bits 64-bit
inline uint64_t popcnt64(uint64_t v) {
    return __popcnt64(v);
}

// bit scan forward
inline bool bsf(uint64_t* bit_idx, uint64_t v) {
    unsigned long i;
    bool result = (bool)_BitScanForward64(&i, v);
    *bit_idx = i;
    return result;
}

#endif

// find first zero Byte
inline uint32_t ffzB(uint32_t x) {
    uint32_t y = (x & 0x7F7F7F7F) + 0x7F7F7F7F;
    y = ~(y | x | 0x7F7F7F7F);
    return (31U - clz(y)) >> 3;
}

// find first zero Byte 64-bit
inline uint64_t ffzB64(uint64_t x) {
    //uint64_t y = (x & 0x7F7F7F7F7F7F7F7F) + 0x7F7F7F7F7F7F7F7F;
    //y = ~(y | x | 0x7F7F7F7F7F7F7F7F);
    //return (63LLU - clz(y)) >> 3LLU;
    
    uint64_t y = (x - 0x0101010101010101) & ~x & 0x8080808080808080;
    return clz64(y) >> 3;
}

#endif