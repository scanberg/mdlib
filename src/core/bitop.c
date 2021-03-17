#include "bitop.h"
#include "intrinsics.h"
#include "common.h"
#include "compiler.h"
#include <string.h> // memset

#if MD_COMPILER_MSVC
#define _CRT_SECURE_NO_WARNINGS
#pragma warning( disable : 6011 )
#endif

#define BITS_PER_BLOCK (sizeof(uint64_t) * 8)
#define MASKED_OP(mask, op, a, b) ((mask) & ((a) op (b)))

static inline uint64_t fast_mod(uint64_t x, uint64_t y) {
    return x & (y - 1);
}

static inline uint64_t bit_pattern(uint64_t bit_idx) {
    return 1LLU << (bit_idx & 63);
}

static inline uint64_t block_idx(uint64_t bit_idx) {
    return bit_idx / BITS_PER_BLOCK;
}

void bit_set(uint64_t* bits, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        bits[beg_idx] |= mask;
        return;
    }
    bits[beg_idx] |= beg_mask;
    memset(bits + beg_idx + 1, 0xFF, (end_idx - beg_idx - 1) * sizeof(uint64_t));
    bits[end_idx] |= end_mask;
}

void bit_clear(uint64_t* bits, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        bits[beg_idx] &= ~(beg_mask & end_mask);
        return;
    }
    
    bits[beg_idx] &= ~beg_mask;
    memset(bits + beg_idx + 1, 0x00, (end_idx - beg_idx - 1) * sizeof(uint64_t));
    bits[end_idx] &= ~end_mask;
}

void bit_invert (uint64_t* bits, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);

    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        bits[beg_idx] = (~mask & bits[beg_idx]) | (mask & ~bits[beg_idx]);
        return;
    }
    bits[beg_idx] = (~beg_mask & bits[beg_idx]) | (beg_mask & ~bits[beg_idx]);
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        bits[i] = ~bits[i];
    }
    bits[end_idx] = (~end_mask & bits[end_idx]) | (end_mask & ~bits[end_idx]);
}

void bit_or(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        dst[beg_idx] = (~mask & dst[beg_idx]) | (mask & (src_a[beg_idx] | src_b[beg_idx]));
        return;
    }
    
    dst[beg_idx] = (~beg_mask & dst[beg_idx]) | (beg_mask & (src_a[beg_idx] | src_b[beg_idx]));
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        dst[i] = (src_a[i] | src_b[i]);
    }
    dst[end_idx] = (~end_mask & dst[end_idx]) | (end_mask & (src_a[end_idx] | src_b[end_idx]));
}

void bit_and(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        dst[beg_idx] = (~mask & dst[beg_idx]) | (mask & (src_a[beg_idx] & src_b[beg_idx]));
        return;
    }

    dst[beg_idx] = (~beg_mask & dst[beg_idx]) | (beg_mask & (src_a[beg_idx] & src_b[beg_idx]));
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        dst[i] = (src_a[i] & src_b[i]);
    }
    dst[end_idx] = (~end_mask & dst[end_idx]) | (end_mask & (src_a[end_idx] & src_b[end_idx]));
}

void bit_xor(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        dst[beg_idx] = (~mask & dst[beg_idx]) | (mask & (src_a[beg_idx] ^ src_b[beg_idx]));
        return;
    }

    dst[beg_idx] = (~beg_mask & dst[beg_idx]) | (beg_mask & (src_a[beg_idx] ^ src_b[beg_idx]));
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        dst[i] = (src_a[i] ^ src_b[i]);
    }
    dst[end_idx] = (~end_mask & dst[end_idx]) | (end_mask & (src_a[end_idx] ^ src_b[end_idx]));
}

void bit_not(uint64_t* dst, const uint64_t* src, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        dst[beg_idx] = (~mask & dst[beg_idx]) | (mask & ~src[beg_idx]);
        return;
    }

    dst[beg_idx] = (~beg_mask & dst[beg_idx]) | (beg_mask & ~src[beg_idx]);
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        dst[i] = ~src[i];
    }
    dst[end_idx] = (~end_mask & dst[end_idx]) | (end_mask & ~src[end_idx]);
}

void bit_or_not(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        dst[beg_idx] = (~mask & dst[beg_idx]) | (mask & (src_a[beg_idx] | ~src_b[beg_idx]));
        return;
    }

    dst[beg_idx] = (~beg_mask & dst[beg_idx]) | (beg_mask & (src_a[beg_idx] | ~src_b[beg_idx]));
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        dst[i] = (src_a[i] | ~src_b[i]);
    }
    dst[end_idx] = (~end_mask & dst[end_idx]) | (end_mask & (src_a[end_idx] | ~src_b[end_idx]));
}

void bit_and_not(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        dst[beg_idx] = (~mask & dst[beg_idx]) | (mask & (src_a[beg_idx] & ~src_b[beg_idx]));
        return;
    }

    dst[beg_idx] = (~beg_mask & dst[beg_idx]) | (beg_mask & (src_a[beg_idx] & ~src_b[beg_idx]));
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        dst[i] = (src_a[i] & ~src_b[i]);
    }
    dst[end_idx] = (~end_mask & dst[end_idx]) | (end_mask & (src_a[end_idx] & ~src_b[end_idx]));
}

/*
void bit_copy(uint64_t* dst, uint64_t dst_offset, const uint64_t* src, uint64_t src_offset, uint64_t bit_count) {
    const uint64_t dst_beg_idx = block_idx(dst_offset);
    const uint64_t dst_end_idx = block_idx(dst_offset + bit_count);
    const uint64_t src_beg_idx = block_idx(src_offset);
    const uint64_t src_end_idx = block_idx(src_offset + bit_count);
}
*/

uint64_t bit_count(const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count) {
    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);
    uint64_t count = 0;

    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        return popcnt64(bits[beg_idx] & mask);
    }
    
    count += popcnt64(bits[beg_idx] & beg_mask);
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        count += popcnt64(bits[i]);
    }
    count += popcnt64(bits[end_idx] & end_mask);

    return count;
}

bool bit_test(const uint64_t* bits, uint64_t idx) {
    const uint64_t block = bits[block_idx(idx)];
    return block & (1ULL << fast_mod(idx, BITS_PER_BLOCK));
}

/*
uint64_t find_next_bit_set(uint64_t* bit_idx, const uint64_t* bits, uint64_t bit_offset, uint64_t bit_count) {
#ifdef ASSERT
    ASSERT(bit_idx);
#endif

    const uint64_t beg_idx = block_idx(bit_offset);
    const uint64_t end_idx = block_idx(bit_offset + bit_count);
    const uint64_t beg_mask = ~(bit_pattern(bit_offset) - 1);
    const uint64_t end_mask = (bit_pattern(bit_offset + bit_count) - 1);

    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        return bit_scan_forward(bit_idx, bits[beg_idx] & mask);
    }

    if (bit_scan_forward(bit_idx, beg_mask & bits[beg_idx])) {
        *bit_idx += bit_offset;
        return true;
    }
    bit_offset += BITS_PER_BLOCK;
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i, bit_offset += BITS_PER_BLOCK) {
        if (bit_scan_forward(bit_idx, bits[i])) {
            *bit_idx += bit_offset;
            return true;
        }
    }
    if (bit_scan_forward(bit_idx, end_mask & bits[end_idx])) {
        *bit_idx += bit_offset;
        return true;
    }
    return false;
}
*/