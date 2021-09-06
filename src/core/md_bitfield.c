#include "md_bitfield.h"
#include "md_allocator.h"
#include "md_log.h"
#include "md_array.inl"

#include "md_common.h"
#include "md_bitop.h"
#include "md_intrinsics.h"

#include <string.h>

#define DIV_UP(x,y) ((x + (y-1)) / y)

#define MAGIC 0xcad86278

/*

TODO:
Fix beg_bit and end_bit to mark the real bits for the range, and not the first bits within the blocks.
This is to help out when sampling the bitfield to easily return the correct value when sampling outside of the range.
We don't really have to care about the block values outside of the range and can just use full block operations.

*/

enum {
    BITFIELD_FLAG_INVERTED = 1
};

typedef struct block {
    uint64_t data[8];
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
    block_t res = {
        a.data[0] & b.data[0],
        a.data[1] & b.data[1],
        a.data[2] & b.data[2],
        a.data[3] & b.data[3],
        a.data[4] & b.data[4],
        a.data[5] & b.data[5],
        a.data[6] & b.data[6],
        a.data[7] & b.data[7],
    };
    return res;
}

static inline block_t block_or(block_t a, block_t b) {
    block_t res = {
        a.data[0] | b.data[0],
        a.data[1] | b.data[1],
        a.data[2] | b.data[2],
        a.data[3] | b.data[3],
        a.data[4] | b.data[4],
        a.data[5] | b.data[5],
        a.data[6] | b.data[6],
        a.data[7] | b.data[7],
    };
    return res;
}

static inline block_t block_not(block_t blk) {
    block_t res = {
        ~blk.data[0],
        ~blk.data[1],
        ~blk.data[2],
        ~blk.data[3],
        ~blk.data[4],
        ~blk.data[5],
        ~blk.data[6],
        ~blk.data[7],
    };
    return res;
}

static inline int64_t block_count_bits(block_t blk) {
    int64_t count =
        popcnt64(blk.data[0]) +
        popcnt64(blk.data[1]) +
        popcnt64(blk.data[2]) +
        popcnt64(blk.data[3]) +
        popcnt64(blk.data[4]) +
        popcnt64(blk.data[5]) +
        popcnt64(blk.data[6]) +
        popcnt64(blk.data[7]);
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
    block_t empty_block = {0};
    return empty_block;
}

static inline void free_blocks(md_exp_bitfield_t* bf) {
    ASSERT(bf);
    ASSERT(bf->bits != NULL);
    md_free(bf->alloc, bf->bits, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));
    bf->bits = NULL;
}

static inline void allocate_blocks(md_exp_bitfield_t* bf, int64_t num_blocks) {
    ASSERT(bf);
    ASSERT(bf->bits == NULL);
    ASSERT(bf->alloc);
    if (num_blocks > 0) {
        bf->bits = md_alloc(bf->alloc, num_blocks * sizeof(block_t));
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
        ASSERT(bf->bits == NULL);
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
        block_t* new_bits = md_alloc(bf->alloc, new_num_blk * sizeof(block_t));
        ASSERT(new_bits);
        memset(new_bits, 0, new_num_blk * sizeof(block_t));

        if (bf->bits) {
            // Copy old data blocks
            const int64_t cur_num_blk = num_blocks(bf->beg_bit, bf->end_bit);
            const int64_t offset_old_in_new = cur_beg_blk - new_beg_blk;
            ASSERT(offset_old_in_new >= 0);
            memcpy(new_bits + offset_old_in_new, bf->bits, cur_num_blk * sizeof(block_t));
            md_free(bf->alloc, bf->bits, cur_num_blk * sizeof(block_t));
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

    if (bf->bits) md_free(bf->alloc, bf->bits, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));
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

/*
void md_bitfield_or(md_exp_bitfield_t* target, const md_exp_bitfield_t* mask) {
    ASSERT(target);
    ASSERT(mask);

    ensure_range(target, mask->beg_bit, mask->end_bit);

    uint64_t* dst = (block_t*)target->bits + block_idx(mask->beg_bit);
    bit_or(dst, dst, (const uint64_t*)mask->bits, 0, mask->end_bit - mask->beg_bit);
}

void md_bitfield_and(md_exp_bitfield_t* target, const md_exp_bitfield_t* mask) {
    ASSERT(target);
    ASSERT(mask);

    ensure_range(target, mask->beg_bit, mask->end_bit);

    md_bitfield_clear_range(target, target->beg_bit, mask->beg_bit);
    md_bitfield_clear_range(target, mask->end_bit, target->end_bit);
    uint64_t* dst = (block_t*)target->bits + block_idx(mask->beg_bit);
    bit_and(dst, dst, (const uint64_t*)mask->bits, 0, mask->end_bit - mask->beg_bit);
}
*/

void md_bitfield_or(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b) {
    validate_bitfield(dst);
    validate_bitfield(src_a);
    validate_bitfield(src_b);

    int64_t beg_bit = MIN(src_a->beg_bit, src_b->beg_bit);
    int64_t end_bit = MAX(src_a->end_bit, src_b->end_bit);

    if (dst == src_a || dst == src_b) {
        ensure_range(dst, beg_bit, end_bit);
    } else {
        fit_to_range(dst, beg_bit, end_bit);
    }
    if (dst->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_or(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_and(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src_a, const md_exp_bitfield_t* src_b) {
    validate_bitfield(dst);
    validate_bitfield(src_a);
    validate_bitfield(src_b);

    int64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    int64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(dst);
        return;
    }

    if (dst == src_a || dst == src_b) {
        ensure_range(dst, beg_bit, end_bit);
    } else {
        fit_to_range(dst, beg_bit, end_bit);
    }
    if (dst->bits) {
        for (int64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_and(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_not(md_exp_bitfield_t* dst, const md_exp_bitfield_t* src, int64_t beg, int64_t end) {
    validate_bitfield(dst);
    validate_bitfield(src);

    if (dst == src) {
        ensure_range(dst, beg, end);
    } else {
        fit_to_range(dst, beg, end);
    }

    if (dst->bits) {
        for (int64_t i = block_idx(beg); i <= block_idx(end); ++i) {
            set_block(dst, i, block_not(get_block(src, i)));
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

/*
int64_t md_bitfield_count_range(const md_exp_bitfield_t* bf, int64_t beg, int64_t end) {
    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    const int64_t count = end - beg;
    if (count == 0) return 0;
    return bit_count((uint64_t*)bf->bits, beg - bf->beg_bit, count);
}
*/

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