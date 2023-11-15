#include <core/md_intrinsics.h>
#include <core/md_common.h>
#include <core/md_compiler.h>

#include <stdint.h>
#include <stdbool.h>

#if MD_COMPILER_MSVC
#pragma warning( disable : 6011 )
#endif

static inline uint64_t bit_pattern(uint64_t bit_idx) {
    return 1LLU << (bit_idx & 63);
}

#define BITOP_PREAMBLE                                     \
    const uint64_t beg_idx = beg_bit >> 6;                 \
    const uint64_t end_idx = end_bit >> 6;                 \
    const uint64_t beg_mask = ~(bit_pattern(beg_bit) - 1); \
    const uint64_t end_mask =  (bit_pattern(end_bit) - 1)

static inline void bit_set(uint64_t* bits, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        bits[beg_idx] |= mask;
        return;
    }
    bits[beg_idx] |= beg_mask;
    MEMSET(bits + beg_idx + 1, 0xFF, (end_idx - beg_idx - 1) * sizeof(uint64_t));
    bits[end_idx] |= end_mask;
}

static inline void bit_set_idx(uint64_t* bits, uint64_t bit_idx) {
    const uint64_t blk_idx = bit_idx >> 6;
    bits[blk_idx] |= bit_pattern(bit_idx);
}

static inline void bit_clear(uint64_t* bits, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
    if (beg_idx == end_idx) {
        bits[beg_idx] &= ~(beg_mask & end_mask);
        return;
    }
    
    bits[beg_idx] &= ~beg_mask;
    MEMSET(bits + beg_idx + 1, 0x00, (end_idx - beg_idx - 1) * sizeof(uint64_t));
    bits[end_idx] &= ~end_mask;
}

static inline void bit_clear_idx(uint64_t* bits, uint64_t bit_idx) {
    const uint64_t blk_idx = bit_idx >> 6;    
    bits[blk_idx] &= ~bit_pattern(bit_idx);
}

static inline void bit_or(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
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

static inline void bit_and(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
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

static inline void bit_xor(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
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

static inline void bit_not(uint64_t* dst, const uint64_t* src, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;

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

static inline void bit_or_not(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
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

static inline void bit_and_not(uint64_t* dst, const uint64_t* src_a, const uint64_t* src_b, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
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


static inline uint64_t bit_count(const uint64_t* bits, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;
    
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

static inline bool bit_test(const uint64_t* bits, uint64_t idx) {
    const uint64_t block = bits[idx >> 6];
    return block & (1ULL << (idx & 63));
}

static inline bool bit_cmp(const uint64_t* src_a, const uint64_t* src_b, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;

    if (beg_idx == end_idx) {
        const uint64_t mask = beg_mask & end_mask;
        return (mask & src_a[beg_idx]) == (mask & src_b[beg_idx]);
    }

    if ((beg_mask & src_a[beg_idx]) != (beg_mask & src_b[beg_idx])) return false;
    if ((end_idx - beg_idx) > 2 && (MEMCMP(src_a + beg_idx + 1, src_b + beg_idx + 1, ((end_idx - beg_idx) - 2) * sizeof(uint64_t)) != 0)) return false;
    if ((end_mask & src_a[end_idx]) != (end_mask & src_b[end_idx])) return false;
    return true;
}

static inline uint64_t bit_scan_forward(const uint64_t* bits, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;

    uint64_t val;
    if (beg_idx == end_idx) {
        val = bits[beg_idx] & (beg_mask & end_mask);
        if (val) return (beg_idx << 6) + ctz64(val) + 1;
        return 0;
    }

    val = bits[beg_idx] & beg_mask;
    if (val) return (beg_idx << 6) + ctz64(val) + 1;
    for (uint64_t i = beg_idx + 1; i < end_idx; ++i) {
        val = bits[i];
        if (val) return (i << 6) + ctz64(val) + 1;
    }
    val = bits[end_idx] & end_mask;
    if (val) return (end_idx << 6) + ctz64(val) + 1;
    
    return 0;
}

static inline uint64_t bit_scan_reverse(const uint64_t* bits, uint64_t beg_bit, uint64_t end_bit) {
    BITOP_PREAMBLE;

    uint64_t val;
    if (beg_idx == end_idx) {
        val = bits[beg_idx] & (beg_mask & end_mask);
        if (val) return ((beg_idx + 1) << 6) - clz64(val);
        return 0;
    }

    val = bits[end_idx] & end_mask;
    if (val) return ((end_idx + 1) << 6) - clz64(val);

    for (uint64_t i = end_idx - 1; i > beg_idx; --i) {
        val = bits[i];
        if (val) return ((i + 1) << 6) - clz64(val);
    }

    val = bits[beg_idx] & beg_mask;
    if (val) return ((beg_idx + 1) << 6) - clz64(val);

    return 0;
}

#undef BITOP_PREAMBLE
