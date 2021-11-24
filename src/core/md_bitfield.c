#include "md_bitfield.h"
#include "md_allocator.h"
#include "md_log.h"
#include "md_array.inl"

#include "md_common.h"
#include "md_bitop.h"
#include "md_intrinsics.h"
#include "md_simd.h"

#include <string.h>
#include <ext/fastlz/fastlz.h>

#define MAGIC 0xcad86278
#define ALIGNMENT sizeof(block_t)

/*

TODO:
Fix beg_bit and end_bit to mark the real bits for the range, and not the first bits within the blocks.
This is to help out when sampling the bitfield to easily return the correct value when sampling outside of the range.
We don't really have to care about the block values outside of the range and can just use full block operations.

*/

typedef union block_t {
    uint64_t u64[8];
    uint32_t u32[16];
    uint16_t u16[32];
    uint8_t   u8[64];
#if md_simd_widthi == 8
    md_simd_typei mm256[2];
#else
    md_simd_typei mm128[4];
#endif
} block_t;

#define BITS_PER_BLOCK (sizeof(block_t) * 8)

static inline int64_t block_idx(int64_t bit_idx) {
    return bit_idx / BITS_PER_BLOCK;
}

static inline int64_t block_bit(int64_t bit_idx) {
    return block_idx(bit_idx) * BITS_PER_BLOCK;
}

static inline int64_t num_blocks(int64_t beg_bit, int64_t end_bit) {
    if (beg_bit == 0 && end_bit == 0) return 0;
    int64_t beg_blk = block_idx(beg_bit);
    int64_t end_blk = block_idx(end_bit);
    return end_blk - beg_blk + 1;
}

static inline block_t block_and(block_t a, block_t b) {
    block_t res;
#ifdef __AVX__
    res.mm256[0] = md_simd_andi(a.mm256[0], b.mm256[0]);
    res.mm256[1] = md_simd_andi(a.mm256[1], b.mm256[1]);
#else
    res.mm128[0] = md_simd_andi(a.mm128[0], b.mm128[0]);
    res.mm128[1] = md_simd_andi(a.mm128[1], b.mm128[1]);
    res.mm128[2] = md_simd_andi(a.mm128[2], b.mm128[2]);
    res.mm128[3] = md_simd_andi(a.mm128[3], b.mm128[3]);
#endif
    /*
    block_t res = {
        a.u64[0] & b.u64[0],
        a.u64[1] & b.u64[1],
        a.u64[2] & b.u64[2],
        a.u64[3] & b.u64[3],
        a.u64[4] & b.u64[4],
        a.u64[5] & b.u64[5],
        a.u64[6] & b.u64[6],
        a.u64[7] & b.u64[7],
    };
    */
    return res;
}

static inline block_t block_andnot(block_t a, block_t b) {
    block_t res;
#ifdef __AVX__
    res.mm256[0] = md_simd_andnoti(a.mm256[0], b.mm256[0]);
    res.mm256[1] = md_simd_andnoti(a.mm256[1], b.mm256[1]);
#else
    res.mm128[0] = md_simd_andnoti(a.mm128[0], b.mm128[0]);
    res.mm128[1] = md_simd_andnoti(a.mm128[1], b.mm128[1]);
    res.mm128[2] = md_simd_andnoti(a.mm128[2], b.mm128[2]);
    res.mm128[3] = md_simd_andnoti(a.mm128[3], b.mm128[3]);
#endif
    /*
    block_t res = {
        a.u64[0] & ~b.u64[0],
        a.u64[1] & ~b.u64[1],
        a.u64[2] & ~b.u64[2],
        a.u64[3] & ~b.u64[3],
        a.u64[4] & ~b.u64[4],
        a.u64[5] & ~b.u64[5],
        a.u64[6] & ~b.u64[6],
        a.u64[7] & ~b.u64[7],
    };
    */
    return res;
}

static inline block_t block_or(block_t a, block_t b) {
    block_t res;
#ifdef __AVX__
    res.mm256[0] = md_simd_ori(a.mm256[0], b.mm256[0]);
    res.mm256[1] = md_simd_ori(a.mm256[1], b.mm256[1]);
#else
    res.mm128[0] = md_simd_ori(a.mm128[0], b.mm128[0]);
    res.mm128[1] = md_simd_ori(a.mm128[1], b.mm128[1]);
    res.mm128[2] = md_simd_ori(a.mm128[2], b.mm128[2]);
    res.mm128[3] = md_simd_ori(a.mm128[3], b.mm128[3]);
#endif
    /*
    block_t res = {
        a.u64[0] | b.u64[0],
        a.u64[1] | b.u64[1],
        a.u64[2] | b.u64[2],
        a.u64[3] | b.u64[3],
        a.u64[4] | b.u64[4],
        a.u64[5] | b.u64[5],
        a.u64[6] | b.u64[6],
        a.u64[7] | b.u64[7],
    };
    */
    return res;
}

static inline block_t block_xor(block_t a, block_t b) {
    block_t res;
#ifdef __AVX__
    res.mm256[0] = md_simd_xori(a.mm256[0], b.mm256[0]);
    res.mm256[1] = md_simd_xori(a.mm256[1], b.mm256[1]);
#else
    res.mm128[0] = md_simd_xori(a.mm128[0], b.mm128[0]);
    res.mm128[1] = md_simd_xori(a.mm128[1], b.mm128[1]);
    res.mm128[2] = md_simd_xori(a.mm128[2], b.mm128[2]);
    res.mm128[3] = md_simd_xori(a.mm128[3], b.mm128[3]);
#endif
    /*
    block_t res = {
        a.u64[0] ^ b.u64[0],
        a.u64[1] ^ b.u64[1],
        a.u64[2] ^ b.u64[2],
        a.u64[3] ^ b.u64[3],
        a.u64[4] ^ b.u64[4],
        a.u64[5] ^ b.u64[5],
        a.u64[6] ^ b.u64[6],
        a.u64[7] ^ b.u64[7],
    };
    */
    return res;
}

static inline block_t block_not(block_t blk) {
#if 0
    block_t res;
#ifdef __AVX__
    res.mm256[0] = md_simd_noti(blk.mm256[0]);
    res.mm256[1] = md_simd_noti(blk.mm256[1]);
#elif (defined(_M_AMD64) || defined(_M_X64))
    res.mm128[0] = md_simd_noti(blk.mm128[0]);
    res.mm128[1] = md_simd_noti(blk.mm128[1]);
    res.mm128[2] = md_simd_noti(blk.mm128[2]);
    res.mm128[3] = md_simd_noti(blk.mm128[3]);
#endif
#endif
    block_t res = {
        ~blk.u64[0],
        ~blk.u64[1],
        ~blk.u64[2],
        ~blk.u64[3],
        ~blk.u64[4],
        ~blk.u64[5],
        ~blk.u64[6],
        ~blk.u64[7],
    };
    return res;
}

// Creates a mask where all the bits less than index (0-512) is set
block_t block_mask_lo(uint32_t idx) {
    block_t res;
#if __AVX2__
    __m256i eq_idx = _mm256_set1_epi32(idx / 32);
    __m256i eq_bits = _mm256_set1_epi32((1UL << (idx & 31)) - 1);

    __m256i lo_idx = _mm256_set_epi32( 7,  6,  5,  4,  3,  2,  1,  0);
    __m256i hi_idx = _mm256_set_epi32(15, 14, 13, 12, 11, 10,  9,  8);

    res.mm256[0] = _mm256_blendv_epi8(_mm256_cmpgt_epi32(eq_idx, lo_idx), eq_bits, _mm256_cmpeq_epi32(lo_idx, eq_idx));
    res.mm256[1] = _mm256_blendv_epi8(_mm256_cmpgt_epi32(eq_idx, hi_idx), eq_bits, _mm256_cmpeq_epi32(hi_idx, eq_idx));
#else
    memset(&res, 0, sizeof(block_t));
    for (uint64_t i = 0; i < idx / 64; ++i) res.u64[i] = ~0ULL;
    res.u64[idx / 64] = ((uint64_t)1 << (idx & 63)) - 1;
#endif
    return res;
}

// Creates a mask where all the bits greater or equal to the index (0-511) is set
block_t block_mask_hi(uint32_t idx) {
    block_t res;

#if __AVX2__
    __m256i eq_idx = _mm256_set1_epi32(idx / 32);
    __m256i eq_bit = _mm256_set1_epi32(~((1UL << (idx & 31)) - 1));

    __m256i lo_idx = _mm256_set_epi32( 7,  6,  5,  4,  3,  2,  1,  0);
    __m256i hi_idx = _mm256_set_epi32(15, 14, 13, 12, 11, 10,  9,  8);

    res.mm256[0] = _mm256_blendv_epi8(_mm256_cmpgt_epi32(lo_idx, eq_idx), eq_bit, _mm256_cmpeq_epi32(lo_idx, eq_idx));
    res.mm256[1] = _mm256_blendv_epi8(_mm256_cmpgt_epi32(hi_idx, eq_idx), eq_bit, _mm256_cmpeq_epi32(hi_idx, eq_idx));
#else
    memset(&res, 0, sizeof(block_t));
    res.u64[idx / 64] = ~(((uint64_t)1 << (idx & 63)) - 1);
    for (uint64_t i = idx / 64 + 1; i < 8; ++i) res.u64[i] = ~0ULL;
#endif
    return res;
}

// According to this source: http://0x80.pl/articles/sse-popcount.html
// There is very little benefit to move beyond the popcnt intrinsic
// Only AVX512 have a popcnt intrinsic for vector registers.
static inline int64_t block_count_bits(block_t blk) {
    int64_t count =
        popcnt64(blk.u64[0]) +
        popcnt64(blk.u64[1]) +
        popcnt64(blk.u64[2]) +
        popcnt64(blk.u64[3]) +
        popcnt64(blk.u64[4]) +
        popcnt64(blk.u64[5]) +
        popcnt64(blk.u64[6]) +
        popcnt64(blk.u64[7]);
    return count;
}

static inline void block_set_bit(block_t* blk, int64_t bit_idx) {
    ASSERT(0 <= bit_idx && bit_idx < BITS_PER_BLOCK);
    blk->u64[bit_idx / 64] |= (1ULL << (bit_idx & 63));
}

static inline bool block_test_bit(block_t blk, int64_t bit_idx) {
    ASSERT(0 <= bit_idx && bit_idx < BITS_PER_BLOCK);
    return blk.u64[bit_idx / 64] & (1ULL << (bit_idx & 63));
}

// Posix convention, returns 0 if no bit is set
// otherwise it returns the index of the first set bit in the interval (1-512)
static inline int64_t block_scan_forward(block_t blk) {
    for (int64_t i = 0; i < ARRAY_SIZE(blk.u64); ++i) {
        uint64_t res = bit_scan_forward64(blk.u64[i]);
        if (res) return res + i * 64;
    }
    return 0;
}

// Posix convention, returns 0 if no bit is set
// otherwise it returns the index of the last set bit in the interval (1-512)
static inline int64_t block_scan_reverse(block_t blk) {
    for (int64_t i = ARRAY_SIZE(blk.u64) - 1; i >= 0; --i) {
        uint64_t res = bit_scan_reverse64(blk.u64[i]);
        if (res) return res + i * 64;
    }
    return 0;
}

static inline void set_block(md_exp_bitfield_t* bf, int64_t idx, block_t blk) {
    ASSERT(bf->bits);
    const int64_t beg_blk = block_idx(bf->beg_bit);
    const int64_t end_blk = block_idx(bf->end_bit);
    ASSERT(beg_blk <= idx && idx <= end_blk);
    ((block_t*)bf->bits)[idx - beg_blk] = blk;
}

static inline block_t get_block(const md_exp_bitfield_t* bf, int64_t blk_idx) {
    const int64_t beg_blk = block_idx(bf->beg_bit);
    const int64_t end_blk = block_idx(bf->end_bit);
    if (bf->bits && beg_blk <= blk_idx && blk_idx <= end_blk)
        return ((block_t*)bf->bits)[blk_idx - beg_blk];
    return (block_t) {0};
}

static inline block_t* get_block_ptr(md_exp_bitfield_t* bf, int64_t blk_idx) {
    const int64_t beg_blk = block_idx(bf->beg_bit);
    const int64_t end_blk = block_idx(bf->end_bit);
    ASSERT(bf->bits && beg_blk <= blk_idx && blk_idx <= end_blk);
    return &((block_t*)bf->bits)[blk_idx - beg_blk];
}

static inline void free_blocks(md_exp_bitfield_t* bf) {
    ASSERT(bf);
    ASSERT(bf->bits != NULL);
    md_aligned_free(bf->alloc, bf->bits, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));
    bf->bits = NULL;
}

static inline void allocate_blocks(md_exp_bitfield_t* bf, int64_t num_blocks) {
    ASSERT(bf);
    ASSERT(bf->bits == NULL);
    ASSERT(bf->alloc);
    if (num_blocks > 0) {
        bf->bits = md_aligned_alloc(bf->alloc, num_blocks * sizeof(block_t), ALIGNMENT);
    }
}

// This fits the range of the bitfield to a bit interval.
// This does not clear the memory so make sure you write to it later
static inline void fit_to_range(md_exp_bitfield_t* bf, int64_t beg_bit, int64_t end_bit) {
    ASSERT(bf);
    ASSERT(beg_bit >= 0);
    ASSERT(beg_bit <= end_bit);

    if (beg_bit == 0 && end_bit == 0) {
        
    }

    const int64_t beg_blk = block_idx(beg_bit);
    const int64_t end_blk = block_idx(end_bit);
    
    if (bf->bits) {
        if (block_idx(bf->beg_bit) == beg_blk && block_idx(bf->end_bit) == end_blk) {
            // We already have the correct blocks, just need to update the bit range
            bf->beg_bit = (uint32_t)beg_bit;
            bf->end_bit = (uint32_t)end_bit;
            return;
        }
        free_blocks(bf);
    }

    bf->beg_bit = (uint32_t)beg_bit;
    bf->end_bit = (uint32_t)end_bit;
    allocate_blocks(bf, num_blocks(bf->beg_bit, bf->end_bit));
}

// This ensures that the bits are represented by the bitfield and will grow it if necessary.
// In case of growth, the bits are preserved.
static inline void ensure_range(md_exp_bitfield_t* bf, int64_t beg_bit, int64_t end_bit) {
    ASSERT(bf);
    ASSERT(beg_bit >= 0);
    ASSERT(beg_bit <= end_bit);

    if (bf->beg_bit <= beg_bit && end_bit <= bf->end_bit) {
        return;
    }

    // If the bitfield is empty, we need to set the bits for min max to work
    if (bf->beg_bit == 0 && bf->end_bit == 0) {
        if (bf->bits) {
            free_blocks(bf);
        }
        bf->beg_bit = (uint32_t)beg_bit;
        bf->end_bit = (uint32_t)end_bit;
        int64_t num_blk = num_blocks(beg_bit, end_bit);
        allocate_blocks(bf, num_blk);
        memset(bf->bits, 0, num_blk * sizeof(block_t));
        return;
    }

    const int64_t new_beg_bit = MIN(bf->beg_bit, beg_bit);
    const int64_t new_end_bit = MAX(bf->end_bit, end_bit);
    const int64_t new_beg_blk = block_idx(new_beg_bit);
    const int64_t new_end_blk = block_idx(new_end_bit);
    const int64_t cur_beg_blk = block_idx(bf->beg_bit);
    const int64_t cur_end_blk = block_idx(bf->end_bit);

    if (new_beg_blk < cur_beg_blk || cur_end_blk < new_end_blk) {
        // Grow!
        const int64_t new_num_blk = num_blocks(new_beg_bit, new_end_bit);
        block_t* new_bits = md_aligned_alloc(bf->alloc, new_num_blk * sizeof(block_t), ALIGNMENT);
        ASSERT(new_bits);
        memset(new_bits, 0, new_num_blk * sizeof(block_t));

        if (bf->bits) {
            // Copy old data blocks
            const int64_t cur_num_blk = num_blocks(bf->beg_bit, bf->end_bit);
            const int64_t offset_old_in_new = cur_beg_blk - new_beg_blk;
            ASSERT(offset_old_in_new >= 0);
            memcpy(new_bits + offset_old_in_new, bf->bits, cur_num_blk * sizeof(block_t));
            md_aligned_free(bf->alloc, bf->bits, cur_num_blk * sizeof(block_t));
        }

        bf->bits = new_bits;
    }

    bf->beg_bit = (uint32_t)new_beg_bit;
    bf->end_bit = (uint32_t)new_end_bit;
}

static inline bool validate_bitfield(const md_exp_bitfield_t* bf) {
    ASSERT(bf);
    ASSERT(bf->magic == MAGIC);
    return true;
}

void md_bitfield_init(md_exp_bitfield_t* bf, struct md_allocator_i* alloc) {
    ASSERT(bf);
    ASSERT(alloc);
    if (bf->magic == MAGIC && bf->bits) {
        md_bitfield_free(bf);
    }
    memset(bf, 0, sizeof(md_exp_bitfield_t));
    bf->magic = MAGIC;
    bf->alloc = alloc;
}

bool md_bitfield_free(md_exp_bitfield_t* bf) {
    validate_bitfield(bf);

    if (bf->bits) md_aligned_free(bf->alloc, bf->bits, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));
    memset(bf, 0, sizeof(md_exp_bitfield_t));
    return true;
}

bool md_bitfield_empty(const md_exp_bitfield_t* bf) {
    return bf->bits == NULL || (bf->beg_bit == bf->end_bit);
}

void md_bitfield_set_range(md_exp_bitfield_t* bf, int64_t beg, int64_t end) {
    validate_bitfield(bf);

    ensure_range(bf, beg, end);
    bit_set((uint64_t*)bf->bits, beg - block_bit(bf->beg_bit), end - beg);
}

void md_bitfield_set_bit(md_exp_bitfield_t* bf, int64_t bit_idx) {
    validate_bitfield(bf);

    ensure_range(bf, bit_idx, bit_idx+1);
    block_set_bit(get_block_ptr(bf, block_idx(bit_idx)), bit_idx - block_bit(bit_idx));
    //bit_set((uint64_t*)bf->bits, bit_idx - block_bit(bf->beg_bit), 1);
}

void md_bitfield_set_indices_u32(md_exp_bitfield_t* bf, uint32_t* indices, int64_t num_indices) {
    validate_bitfield(bf);
    ASSERT(indices);
    ASSERT(num_indices >= 0);

    for (int64_t i = 0; i < num_indices; ++i) {
        int64_t bit_idx = (int64_t)indices[i];
        ensure_range(bf, bit_idx, bit_idx+1);
        bit_set((uint64_t*)bf->bits, bit_idx - block_bit(bf->beg_bit), 1);
    }
}

void md_bitfield_clear(md_exp_bitfield_t* bf) {
    validate_bitfield(bf);

    memset(bf->bits, 0, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));
}

void md_bitfield_clear_range(md_exp_bitfield_t* bf, int64_t beg, int64_t end) {
    validate_bitfield(bf);
        
    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    if (end-beg == 0) return;
    bit_clear((uint64_t*)bf->bits, beg - block_bit(bf->beg_bit), end - beg);
}

void md_bitfield_clear_bit(md_exp_bitfield_t* bf, int64_t bit_idx) {
    validate_bitfield(bf);

    if (bf->beg_bit <= bit_idx && bit_idx < bf->end_bit) {
        bit_clear_idx((uint64_t*)bf->bits, bit_idx - block_bit(bf->beg_bit));
    }
}

void md_bitfield_or(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b) {
    validate_bitfield(dst);
    validate_bitfield(src_a);
    validate_bitfield(src_b);

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    int64_t beg_bit = MIN(src_a->beg_bit, src_b->beg_bit);
    int64_t end_bit = MAX(src_a->end_bit, src_b->end_bit);

    fit_to_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_or(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_or_inplace(md_exp_bitfield_t* a, const md_exp_bitfield_t* b) {
    validate_bitfield(a);
    validate_bitfield(b);

    if (a == b) return;
    if (a->bits == 0) {
        md_bitfield_copy(a, b);
        return;
    }

    int64_t beg_bit = MIN(a->beg_bit, b->beg_bit);
    int64_t end_bit = MAX(a->end_bit, b->end_bit);
    ensure_range(a, beg_bit, end_bit);

    if (a->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_or(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_and(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b) {
    validate_bitfield(dst);
    validate_bitfield(src_a);
    validate_bitfield(src_b);

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    int64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    int64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(dst);
        return;
    }

    fit_to_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_and(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_and_inplace(md_exp_bitfield_t* a, const md_exp_bitfield_t* b) {
    validate_bitfield(a);
    validate_bitfield(b);

    if (a == b) return;

    int64_t beg_bit = MAX(a->beg_bit, b->beg_bit);
    int64_t end_bit = MIN(a->end_bit, b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(a);
        return;
    }

    if (a->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_and(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_andnot(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b) {
    validate_bitfield(dst);
    validate_bitfield(src_a);
    validate_bitfield(src_b);

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    if (src_a == src_b) {
        md_bitfield_clear(dst);
        return;
    }

    int64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    int64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);
    
    md_bitfield_copy(dst, src_a);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap, done.
        return;
    }

    ensure_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_andnot(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_andnot_inplace(md_exp_bitfield_t* a, const md_exp_bitfield_t* b) {
    validate_bitfield(a);
    validate_bitfield(b);

    if (a == b) {
        md_bitfield_clear(a);
        return;
    }

    int64_t beg_bit = MAX(a->beg_bit, b->beg_bit);
    int64_t end_bit = MIN(a->end_bit, b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap, done.
        return;
    }

    // The only place we need to modify is where the two ranges overlap.
    // ensure_range(a, beg_bit, end_bit);

    if (a->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_andnot(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_xor(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b) {
    validate_bitfield(dst);
    validate_bitfield(src_a);
    validate_bitfield(src_b);

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    int64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    int64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(dst);
        return;
    }

    fit_to_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_xor(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_xor_inplace(md_exp_bitfield_t* a, const md_exp_bitfield_t* b) {
    validate_bitfield(a);
    validate_bitfield(b);

    if (a == b) return;

    int64_t beg_bit = MAX(a->beg_bit, b->beg_bit);
    int64_t end_bit = MIN(a->end_bit, b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(a);
        return;
    }

    ensure_range(a, beg_bit, end_bit);
    if (a->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_xor(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_not(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src, int64_t beg, int64_t end) {
    validate_bitfield(dst);
    validate_bitfield(src);

    ASSERT((dst != src) && "dst is same as src, use inplace version instead!");

    fit_to_range(dst, beg, end);

    if (dst->bits) {
        for (int64_t i = block_idx(beg); i <= block_idx(end); ++i) {
            set_block(dst, i, block_not(get_block(src, i)));
        }
    }
}

void md_bitfield_not_inplace(md_exp_bitfield_t* bf, int64_t beg, int64_t end) {
    validate_bitfield(bf);

    ensure_range(bf, beg, end);

    if (bf->bits) {
        for (int64_t i = block_idx(beg); i <= block_idx(end); ++i) {
            set_block(bf, i, block_not(get_block(bf, i)));
        }
    }
}

void md_bitfield_copy(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src) {
    validate_bitfield(dst);
    validate_bitfield(src);
    if (dst == src) return;

    fit_to_range(dst, src->beg_bit, src->end_bit);
    if (dst->bits) {
        memcpy(dst->bits, src->bits, num_blocks(dst->beg_bit, dst->end_bit) * sizeof(block_t));
    }
    dst->flags = src->flags;
}

// Counts the number of bits set
int64_t md_bitfield_popcount(const md_exp_bitfield_t* bf) {
    validate_bitfield(bf);

    const int64_t count = (int64_t)bf->end_bit - (int64_t)bf->beg_bit;
    if (count == 0) return 0;

    // If we can find a way to efficiently generating a mask for a block we should use the block_popcount proc instead...

    return bit_count((uint64_t*)bf->bits, bf->beg_bit - block_bit(bf->beg_bit), count);
}

int64_t md_bitfield_popcount_range(const md_exp_bitfield_t* bf, int64_t beg, int64_t end) {
    ASSERT(beg <= end);
    validate_bitfield(bf);

    if (bf->bits == NULL) return 0;
    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    const int64_t count = end - beg;
    if (count == 0) return 0;
    return bit_count((uint64_t*)bf->bits, beg - block_bit(bf->beg_bit), count);
}

// Test if single bit in field is set
bool md_bitfield_test_bit(const md_exp_bitfield_t* bf, int64_t idx) {
    validate_bitfield(bf);

    if (idx < bf->beg_bit || bf->end_bit <= idx) return false;
    return block_test_bit(get_block(bf, block_idx(idx)), idx - block_bit(idx));
    //return bit_test((uint64_t*)bf->bits, idx - block_bit(bf->beg_bit));
}

// Test if bitfields are equivalent
//bool md_bitfield_cmp        (const uint64_t* src_a, const uint64_t* src_b);
//bool md_bitfield_cmp_range  (const uint64_t* src_a, const uint64_t* src_b, int64_t beg, int64_t end);

// Bit scan forward, finds the first bit set from the given offset.
// Returns 0 if no bit is found, otherwise it returns the offset to that bit (indexing starts at 1, posix convention)
int64_t md_bitfield_scan(const md_exp_bitfield_t* bf, int64_t beg, int64_t end) {
    validate_bitfield(bf);

    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    if (beg == end)
        return 0;
    int64_t result = bit_scan((uint64_t*)bf->bits, beg - block_bit(bf->beg_bit), end - beg);
    if (!result)
        return 0;
    return result + block_bit(bf->beg_bit);
}

static inline uint64_t bit_pattern(uint64_t bit_idx) {
    return 1LLU << (bit_idx & 63);
}

static inline uint64_t get_u64(const md_exp_bitfield_t* bf, int64_t idx_u64) {
    if (!bf->bits) return 0;

    const int64_t beg_u64 = bf->beg_bit / 64;
    const int64_t end_u64 = bf->end_bit / 64;
    
    if (beg_u64 < idx_u64 && idx_u64 < end_u64) {
        return ((uint64_t*)bf->bits)[idx_u64 - beg_u64];
    }

    uint64_t result = 0;
    if (idx_u64 == beg_u64) {
        const uint64_t beg_mask = ~(bit_pattern(bf->beg_bit) - 1);
        result |= beg_mask & ((uint64_t*)bf->bits)[beg_u64];
    }

    if (idx_u64 == end_u64) {
        const uint64_t end_mask = (bit_pattern(bf->end_bit) - 1);
        result |= end_mask & ((uint64_t*)bf->bits)[end_u64];
    }

    return result;
}

bool md_bitfield_extract_u64(uint64_t* dst_ptr, int64_t num_bits, const md_exp_bitfield_t* src) {
    ASSERT(dst_ptr);
    ASSERT(num_bits >= 0);
    validate_bitfield(src);
    
    if (num_bits == 0) return true;

    const int64_t count = (num_bits + 63) / 64;
    for (int64_t i = 0; i < count; ++i) {
        dst_ptr[i] = get_u64(src, i);
    }
    
    return true;
}

uint32_t* md_bitfield_extract_indices_u32(const md_exp_bitfield_t* bf, md_allocator_i* alloc) {
    ASSERT(bf);
    ASSERT(alloc);
    validate_bitfield(bf);

    uint32_t* indices = 0;

    int64_t beg_bit = bf->beg_bit;
    int64_t end_bit = bf->end_bit;
    while ((beg_bit = md_bitfield_scan(bf, beg_bit, end_bit)) != 0) {
        uint32_t bit_idx = (uint32_t)(beg_bit - 1);
        md_array_push(indices, bit_idx, alloc);
    }

    return indices;
}

uint32_t* md_bitfield_extract_bits_u32(const md_exp_bitfield_t* bf, struct md_allocator_i* alloc) {
    ASSERT(bf);
    ASSERT(alloc);
    validate_bitfield(bf);

    uint32_t* bits = 0;

    int64_t block_count = num_blocks(bf->beg_bit, bf->end_bit);
    if (block_count) {
        for (int64_t i = block_idx(bf->beg_bit); i <= block_idx(bf->end_bit); ++i) {
            block_t block = get_block(bf, i);
            for (int64_t j = 0; j < ARRAY_SIZE(block.u32); ++j) {
                md_array_push(bits, block.u32[j], alloc);
            }
        }
    }

    return bits;
}

#define BLOCK_IDX_FLAG_ALL_SET (0x8000)

uint32_t* get_non_empty_block_indices(const md_exp_bitfield_t* bf, md_allocator_i* alloc) {
    uint32_t* indices = 0;

    int64_t beg_blk_idx = block_idx(bf->beg_bit);
    int64_t end_blk_idx = block_idx(bf->end_bit);
    for (int64_t i = beg_blk_idx; i <= end_blk_idx; ++i) {
        block_t blk = get_block(bf, i);
        if (i == beg_blk_idx) block_and(blk, block_mask_hi(bf->beg_bit & 511));
        if (i == beg_blk_idx) block_and(blk, block_mask_lo(bf->end_bit & 511));

        bool empty = true;
        for (int j = 0; j < ARRAY_SIZE(blk.u64); ++j) {
            if (blk.u64[j] != 0) empty = false;
        }
        if (empty) continue;
        md_array_push(indices, (uint32_t)i, alloc);
    }

    return indices;
}

uint16_t* get_serialization_block_indices(const md_exp_bitfield_t* bf, md_allocator_i* alloc) {
    uint16_t* indices = 0;

    int64_t beg_blk_idx = block_idx(bf->beg_bit);
    int64_t end_blk_idx = block_idx(bf->end_bit);
    for (int64_t i = beg_blk_idx; i <= end_blk_idx; ++i) {
        block_t blk = get_block(bf, i);
        if (i == beg_blk_idx) block_and(blk, block_mask_hi(bf->beg_bit & 511));
        if (i == beg_blk_idx) block_and(blk, block_mask_lo(bf->end_bit & 511));

        bool empty = true;
        bool all_set = true;
        for (int j = 0; j < ARRAY_SIZE(blk.u64); ++j) {
            if (blk.u64[j] != 0) empty = false;
            if (blk.u64[j] != 0xFFFFFFFFFFFFFFFF) all_set = false;
        }

        if (empty) continue;

        // Store non empty block indices with additional flag to see if it is empty or not
        md_array_push(indices, (uint16_t)i | (all_set ? BLOCK_IDX_FLAG_ALL_SET : 0), alloc);
    }
    return indices;
}

// Returns the maximum serialization size in bytes of a bitfield
int64_t md_bitfield_serialize_size_in_bytes(const md_exp_bitfield_t* bf) {
    uint64_t size = sizeof(uint16_t);
    uint16_t* indices = get_serialization_block_indices(bf, default_temp_allocator);
    for (int64_t i = 0; i < md_array_size(indices); ++i) {
        if (indices[i] & BLOCK_IDX_FLAG_ALL_SET) continue;
        size += sizeof(indices[i]) + sizeof(block_t);
    }

    return MAX(66, size);
}

typedef union block_idx_t {
    uint16_t u16;
    struct {
        uint8_t occupancy;
        uint8_t offset_to_next;
    };
} block_idx_t;

// Serializes a bitfield into a destination buffer
// It is expected that the supplied buffer has the size_in_bytes supplied by bitfield_serialize_size_in_bytes()
int64_t md_bitfield_serialize(void* dst, const md_exp_bitfield_t* bf) {
    md_allocator_i* alloc = default_temp_allocator;

    uint16_t* indices = get_serialization_block_indices(bf, alloc);
    uint16_t* data = 0;

    uint16_t block_count = (uint16_t)md_array_size(indices);
    int64_t beg_blk = block_idx(bf->beg_bit);
    int64_t end_blk = block_idx(bf->end_bit);
    
    md_array_push(data, block_count, alloc);
    md_array_push_array(data, indices, block_count, alloc);

    // Add actual block data from the blocks (if not all bits set within block)
    for (int64_t i = 0; i < block_count; ++i) {
        int64_t blk_idx = data[1 + i];
        bool    all_set = data[1 + i] & BLOCK_IDX_FLAG_ALL_SET;

        if (!all_set) {
            block_t blk = get_block(bf, blk_idx);
            if (blk_idx == beg_blk) block_and(blk, block_mask_hi(bf->beg_bit & 511));
            if (blk_idx == end_blk) block_and(blk, block_mask_lo(bf->end_bit & 511));
            md_array_push_array(data, blk.u16, ARRAY_SIZE(blk.u16), alloc);
        }
    }

    int lz_bytes = fastlz_compress_level(2, data, (int)(md_array_size(data) * sizeof(uint16_t)), dst);
    md_printf(MD_LOG_TYPE_INFO, "LZ:\t%lli number of bits compressed into %i bytes\n", md_bitfield_popcount(bf), lz_bytes);

    md_array_free(data, alloc);
    md_array_free(indices, alloc);

    return lz_bytes;
}

bool md_bitfield_deserialize(md_exp_bitfield_t* bf, const void* src, int64_t num_bytes) {
    ASSERT(bf);
    validate_bitfield(bf);
    ASSERT(src);

    // Temporary buffer to decompress into.
    // The current compression is no where near 100:1 ratio.
    const int64_t mem_bytes = num_bytes * 100;
    void* mem = md_alloc(default_allocator, mem_bytes);

    int size = fastlz_decompress(src, (int)num_bytes, mem, (int)mem_bytes);

    if (size == 0) {
        md_printf(MD_LOG_TYPE_ERROR, "Failed, to decompress bitfield.");
        return false;
    }

    const uint16_t* data = (const uint16_t*)mem;
    uint16_t block_count = data[0];
    const uint16_t* block_indices = data + 1;
    const block_t* block_data = (const block_t*)(block_indices + block_count);

    if (block_count == 0) {
        md_printf(MD_LOG_TYPE_ERROR, "Block count was zero");
        return false;
    }

    uint16_t beg_blk_idx = block_indices[0];
    uint16_t end_blk_idx = block_indices[block_count - 1];

    // Allocate the blocks
    fit_to_range(bf, beg_blk_idx * 512, end_blk_idx * 512 + 511);
    memset(bf->bits, 0, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));

    // Fetch block_data and store
    block_t* dst_block = (block_t*)bf->bits;
    int64_t src_offset = 0;
    for (int64_t i = 0; i < block_count; ++i) {
        uint16_t blk_idx = block_indices[i];
        if (blk_idx & BLOCK_IDX_FLAG_ALL_SET) {
            blk_idx &= ~BLOCK_IDX_FLAG_ALL_SET;
            memset(dst_block + blk_idx, 0xFFFFFFFF, sizeof(block_t));
        } else {
            memcpy(dst_block + blk_idx, block_data + src_offset, sizeof(block_t));
            src_offset += 1;
        }
    }

    // Now we have the data, compute the true beg_bit and end_bit
    int64_t beg_bit = block_scan_forward(get_block(bf, beg_blk_idx));
    int64_t end_bit = block_scan_reverse(get_block(bf, end_blk_idx));

    bf->beg_bit = (uint32_t)(beg_blk_idx * 512 + (beg_bit ? beg_bit - 1 : 0));
    bf->end_bit = (uint32_t)(end_blk_idx * 512 + (end_bit ? end_bit : 0));

    md_free(default_allocator, mem, mem_bytes);

    return true;
}