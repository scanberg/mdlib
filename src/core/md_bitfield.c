#include "md_bitfield.h"
#include "md_allocator.h"
#include "md_log.h"
#include "md_array.inl"

#include "md_common.h"
#include "md_bitop.h"
#include "md_intrinsics.h"
#include "md_simd.h"

#include <string.h>

#define MAGIC 0xcad86278
#define ALIGNMENT sizeof(block_t)

/*

TODO:
Fix beg_bit and end_bit to mark the real bits for the range, and not the first bits within the blocks.
This is to help out when sampling the bitfield to easily return the correct value when sampling outside of the range.
We don't really have to care about the block values outside of the range and can just use full block operations.

*/

typedef union block {
    uint64_t u64[8];
    uint32_t u32[16];
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

// Generate mask for all bits lower than the specified index
static inline block_t block_mask_lo(uint64_t idx) {
    ASSERT(idx < BITS_PER_BLOCK);
    block_t res;
    uint64_t j = idx / 64;
    for (uint64_t i = 0; i < j; ++i) res.u64[0] = (uint64_t)~0;
    res.u64[j] = (uint64_t)1 << (idx & 63);
    for (uint64_t i = j + 1; i < ARRAY_SIZE(res.u64); ++i) res.u64[0] = 0;

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


static inline void set_block(md_exp_bitfield_t* bf, int64_t idx, block_t blk) {
    ASSERT(bf->bits);
    const int64_t beg_blk = block_idx(bf->beg_bit);
    const int64_t end_blk = block_idx(bf->end_bit);
    ASSERT(beg_blk <= idx && idx <= end_blk);
    ((block_t*)bf->bits)[idx - beg_blk] = blk;
}

static inline block_t get_block(const md_exp_bitfield_t* bf, int64_t idx) {
    const int64_t beg_blk = block_idx(bf->beg_bit);
    const int64_t end_blk = block_idx(bf->end_bit);
    if (bf->bits && beg_blk <= idx && idx <= end_blk)
        return ((block_t*)bf->bits)[idx - beg_blk];
    return (block_t) {0};
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

    ensure_range(bf, bit_idx, bit_idx + 1);
    bit_set((uint64_t*)bf->bits, bit_idx - block_bit(bf->beg_bit), 1);
}

void md_bitfield_set_indices_u32(md_exp_bitfield_t* bf, uint32_t* indices, int64_t num_indices) {
    validate_bitfield(bf);
    ASSERT(indices);
    ASSERT(num_indices >= 0);

    for (int64_t i = 0; i < num_indices; ++i) {
        int64_t bit_idx = (int64_t)indices[i];
        ensure_range(bf, bit_idx, bit_idx + 1);
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
    return bit_test((uint64_t*)bf->bits, idx - block_bit(bf->beg_bit));
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