#include <core/md_bitfield.h>

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_intrinsics.h>
#include <core/md_simd.h>
#include <core/md_bitop.inl>
#include <core/md_hash.h>

#include <fastlz.h>

#define MAGIC 0xcad86278
#define ALIGNMENT sizeof(block_t)

#if MD_COMPILER_MSVC
#pragma warning(disable : 4701)  // potentially uninitialized local variable used
#endif

typedef union block_t {
    md_256i    v[2];
    uint64_t u64[8];
    uint32_t u32[16];
    uint16_t u16[32];
    uint8_t   u8[64];
} block_t;

#define BITS_PER_BLOCK (sizeof(block_t) * 8)

// This is the number of bits we shift the bit index to get the block index (should equal log2(BITS_PER_BLOCK))
#define BASE_SHIFT 6
#define BASE_MASK  ~0x7LLU

static inline const uint64_t* u64_base(const md_bitfield_t* bf) {
    ASSERT(bf);
    const uint64_t offset = ((uint64_t)bf->beg_bit >> BASE_SHIFT) & BASE_MASK;
    return (const uint64_t*)bf->bits - offset;
}

static inline uint64_t* u64_base_mut(md_bitfield_t* bf) {
    ASSERT(bf);
    const uint64_t offset = ((uint64_t)bf->beg_bit >> BASE_SHIFT) & BASE_MASK;
    return (uint64_t*)bf->bits - offset;
}

static inline uint64_t block_idx(uint64_t bit_idx) {
    return bit_idx / BITS_PER_BLOCK;
}

static inline uint64_t block_bit(uint64_t bit_idx) {
    return block_idx(bit_idx) * BITS_PER_BLOCK;
}

static inline uint64_t num_blocks(uint64_t beg_bit, uint64_t end_bit) {
    if (beg_bit == 0 && end_bit == 0) return 0;
    uint64_t beg_blk = block_idx(beg_bit);
    uint64_t end_blk = block_idx(end_bit);
    return end_blk - beg_blk + 1;
}

static inline block_t block_and(block_t a, block_t b) {
    block_t res;
    for (size_t i = 0; i < ARRAY_SIZE(a.v); ++i) {
        res.v[i] = md_mm256_and_si256(a.v[i], b.v[i]);
    };
    return res;
}

static inline block_t block_andnot(block_t a, block_t b) {
    block_t res;
    // The arguments are flipped here because it is the first operand in the intrinsic that is bitwise NOT
    for (int i = 0; i < (int)ARRAY_SIZE(a.v); ++i) {
        res.v[i] = md_mm256_andnot_si256(b.v[i], a.v[i]);
    }
    return res;
}

static inline block_t block_or(block_t a, block_t b) {
    block_t res;
    for (int i = 0; i < (int)ARRAY_SIZE(a.v); ++i) {
        res.v[i] = md_mm256_or_si256(a.v[i], b.v[i]);
    }
    return res;
}

static inline block_t block_xor(block_t a, block_t b) {
    block_t res;
    for (int i = 0; i < (int)ARRAY_SIZE(a.v); ++i) {
        res.v[i] = md_mm256_xor_si256(a.v[i], b.v[i]);
    }
    return res;
}

static inline block_t block_not(block_t blk) {
    block_t res;
    for (int i = 0; i < (int)ARRAY_SIZE(blk.v); ++i) {
        res.v[i] = md_mm256_not_si256(blk.v[i]);
    }
    return res;
}

// Creates a mask where all the bits less than index (0-512) is set
block_t block_mask_lo(uint32_t idx) {
    block_t res;

    __m256i eq_idx = md_mm256_set1_epi64(idx / 64);
    __m256i eq_bit = md_mm256_set1_epi64((1ULL << (idx & 63)) - 1);

    __m256i lo_idx = md_mm256_set_epi64(3, 2, 1, 0);
    __m256i hi_idx = md_mm256_set_epi64(7, 6, 5, 4);

    res.v[0] = md_mm256_blendv_epi8(md_mm256_cmpgt_epi64(eq_idx, lo_idx), eq_bit, md_mm256_cmpeq_epi64(lo_idx, eq_idx));
    res.v[1] = md_mm256_blendv_epi8(md_mm256_cmpgt_epi64(eq_idx, hi_idx), eq_bit, md_mm256_cmpeq_epi64(hi_idx, eq_idx));
#if 0
    MEMSET(&res, 0, sizeof(block_t));
    for (uint64_t i = 0; i < idx / 64; ++i) res.u64[i] = ~0ULL;
    res.u64[idx / 64] = ((uint64_t)1 << (idx & 63)) - 1;
#endif
    return res;
}

// Creates a mask where all the bits greater or equal to the index (0-511) is set
block_t block_mask_hi(uint32_t idx) {
    block_t res;

    __m256i eq_idx = md_mm256_set1_epi64(idx / 32);
    __m256i eq_bit = md_mm256_set1_epi64(~((1UL << (idx & 31)) - 1));

    __m256i lo_idx = md_mm256_set_epi64(3, 2, 1, 0);
    __m256i hi_idx = md_mm256_set_epi64(7, 6, 5, 4);

    res.v[0] = md_mm256_blendv_epi8(md_mm256_cmpgt_epi64(lo_idx, eq_idx), eq_bit, md_mm256_cmpeq_epi64(lo_idx, eq_idx));
    res.v[1] = md_mm256_blendv_epi8(md_mm256_cmpgt_epi64(hi_idx, eq_idx), eq_bit, md_mm256_cmpeq_epi64(hi_idx, eq_idx));
#if 0
    MEMSET(&res, 0, sizeof(block_t));
    res.u64[idx / 64] = ~(((uint64_t)1 << (idx & 63)) - 1);
    for (uint64_t i = idx / 64 + 1; i < 8; ++i) res.u64[i] = ~0ULL;
#endif
    return res;
}

// According to this source: http://0x80.pl/articles/sse-popcount.html
// There is very little benefit to move beyond the popcnt intrinsic
// Only AVX512 have a popcnt intrinsic for vector registers.
static inline uint64_t block_count_bits(block_t blk) {
    uint64_t count = 0;
    for (size_t i = 0; i < ARRAY_SIZE(blk.u64); ++i) {
        count += popcnt64(blk.u64[i]);
    }
    return count;
}

static inline void block_set_bit(block_t* blk, uint64_t bit_idx) {
    ASSERT(bit_idx < BITS_PER_BLOCK);
    blk->u64[bit_idx >> 6] |= (1ULL << (bit_idx & 63));
}

static inline bool block_test_bit(block_t blk, uint64_t bit_idx) {
    ASSERT(bit_idx < BITS_PER_BLOCK);
    return blk.u64[bit_idx >> 6] & (1ULL << (bit_idx & 63));
}

// Posix convention, returns 0 if no bit is set
// otherwise it returns the index of the first set bit in the interval (1-512)
static inline uint64_t block_scan_forward(block_t blk) {
    for (uint64_t i = 0; i < ARRAY_SIZE(blk.u64); ++i) {
        uint64_t res = bsf64(blk.u64[i]);
        if (res) return res + i * 64;
    }
    return 0;
}

// Posix convention, returns 0 if no bit is set
// otherwise it returns the index of the last set bit in the interval (1-512)
static inline uint64_t block_scan_reverse(block_t blk) {
    for (int64_t i = (int64_t)ARRAY_SIZE(blk.u64) - 1; i >= 0; --i) {
        uint64_t res = bsr64(blk.u64[i]);
        if (res) return res + i * 64;
    }
    return 0;
}

static inline void set_block(md_bitfield_t* bf, uint64_t idx, block_t blk) {
    ASSERT(bf->bits);
    const uint64_t beg_blk = block_idx(bf->beg_bit);
    ASSERT(beg_blk <= idx && idx <= block_idx(bf->end_bit));
    MEMCPY(((block_t*)bf->bits) + (idx - beg_blk), &blk, sizeof(block_t));
}

static inline block_t get_block(const md_bitfield_t* bf, uint64_t idx) {
    const uint64_t beg_blk = block_idx(bf->beg_bit);
    const uint64_t end_blk = block_idx(bf->end_bit);
    if (bf->bits && beg_blk <= idx && idx <= end_blk) {
        block_t blk;
        MEMCPY(&blk, ((block_t*)bf->bits) + (idx - beg_blk), sizeof(block_t));
        return blk;
    }
    return (block_t){0};
}

#define get_blk_ptr(bf, bit_idx) ((block_t*)bf->bits + (bit_idx - block_idx(bf->beg_bit)))

static inline void free_blocks(md_allocator_i* alloc, void* ptr, size_t num_blocks) {
    if (ptr) {
        md_free(alloc, ptr, num_blocks * sizeof(block_t));
        //md_aligned_free(alloc, ptr);
    }
}

static inline void* alloc_blocks(md_allocator_i* alloc, size_t num_blocks) {
    void* mem = 0;
    if (num_blocks > 0) {
        mem = md_alloc(alloc, num_blocks * sizeof(block_t));
        //mem = md_aligned_alloc(alloc, num_blocks * sizeof(block_t), sizeof(block_t));
    }
    return mem;
}

static inline void* realloc_blocks(md_allocator_i* alloc, void* ptr, size_t cur_blocks, size_t new_blocks) {
    return md_realloc(alloc, ptr, cur_blocks * sizeof(block_t), new_blocks * sizeof(block_t));
    //return md_aligned_realloc(alloc, ptr, num_blocks * sizeof(block_t), sizeof(block_t));
}

// This fits the range of the bitfield to a bit interval.
// This does not clear the memory so make sure you write to it later
static inline void fit_to_range(md_bitfield_t* bf, uint64_t beg_bit, uint64_t end_bit) {
    ASSERT(bf);
    ASSERT(beg_bit <= end_bit);

    const uint64_t beg_blk = block_idx(beg_bit);
    const uint64_t end_blk = block_idx(end_bit);
    const uint64_t cur_beg_blk = block_idx(bf->beg_bit);
    const uint64_t cur_end_blk = block_idx(bf->end_bit);

    uint64_t cur_blocks = num_blocks(bf->beg_bit, bf->end_bit);
    uint64_t new_blocks = num_blocks(beg_bit, end_bit);

    bf->beg_bit = (uint32_t)beg_bit;
    bf->end_bit = (uint32_t)end_bit;
    
    if (bf->bits && cur_beg_blk == beg_blk && cur_end_blk == end_blk) {
        return;
    }

    bf->bits = realloc_blocks(bf->alloc, bf->bits, cur_blocks, new_blocks);
}

// This ensures that the bits are represented by the bitfield and will grow it if necessary.
// In case of growth, the bits are preserved.
static inline void ensure_range(md_bitfield_t* bf, uint64_t beg_bit, uint64_t end_bit) {
    ASSERT(bf);
    ASSERT(beg_bit <= end_bit);

    if (bf->beg_bit <= beg_bit && end_bit <= bf->end_bit) {
        return;
    }

    // If the bitfield is empty, we need to set the bits for min max to work
    if (bf->beg_bit == 0 && bf->end_bit == 0) {
        if (bf->bits) {
            ASSERT(false);
            free_blocks(bf->alloc, bf->bits, 0);
            bf->bits = 0;
        }
        bf->beg_bit = (uint32_t)beg_bit;
        bf->end_bit = (uint32_t)end_bit;
        uint64_t num_blk = num_blocks(beg_bit, end_bit);
        if (num_blk > 0) {
            bf->bits = alloc_blocks(bf->alloc, num_blk);
            ASSERT(bf->bits);
            MEMSET(bf->bits, 0, num_blk * sizeof(block_t));
        }
        return;
    }

    const uint64_t new_beg_bit = MIN(bf->beg_bit, beg_bit);
    const uint64_t new_end_bit = MAX(bf->end_bit, end_bit);
    const uint64_t new_beg_blk = block_idx(new_beg_bit);
    const uint64_t new_end_blk = block_idx(new_end_bit);
    const uint64_t cur_beg_blk = block_idx(bf->beg_bit);
    const uint64_t cur_end_blk = block_idx(bf->end_bit);

    if (new_beg_blk < cur_beg_blk || cur_end_blk < new_end_blk) {
        // Grow!
        // @NOTE: This is not great. Its not even good...
        // The entire internal workings of the bitfield should change
        // To a sparse layout of blocks that are indexed by a hashmap

        const uint64_t cur_num_blk = num_blocks(bf->beg_bit, bf->end_bit);
        const uint64_t new_num_blk = num_blocks(new_beg_bit, new_end_bit);

#if 1
        if (new_beg_blk == cur_beg_blk) {
            // The starting block has not changed, meaning we can use realloc
            block_t* new_bits = realloc_blocks(bf->alloc, bf->bits, cur_num_blk, new_num_blk);

            // Only clear newly allocated region
            MEMSET(new_bits + cur_num_blk, 0, (new_num_blk - cur_num_blk) * sizeof(block_t));

            bf->bits = new_bits;
        } else {
#endif
            block_t* new_bits = alloc_blocks(bf->alloc, new_num_blk);
            MEMSET(new_bits, 0, new_num_blk * sizeof(block_t));

            if (bf->bits) {
                size_t blk_diff = cur_beg_blk - new_beg_blk;
                MEMCPY(new_bits + blk_diff, bf->bits, cur_num_blk * sizeof(block_t));
                free_blocks(bf->alloc, bf->bits, cur_num_blk);
            }

            bf->bits = new_bits;
#if 1
        }
#endif
    }

    bf->beg_bit = (uint32_t)new_beg_bit;
    bf->end_bit = (uint32_t)new_end_bit;
}

md_bitfield_t md_bitfield_create(md_allocator_i* alloc) {
    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);
    return bf;
}

void md_bitfield_init(md_bitfield_t* bf, md_allocator_i* alloc) {
    ASSERT(bf);
    ASSERT(alloc);
    if (bf->magic == MAGIC && bf->bits) {
        md_bitfield_free(bf);
    }
    *bf = (md_bitfield_t){
        .magic = MAGIC,
        .alloc = alloc,
    };
}

bool md_bitfield_free(md_bitfield_t* bf) {
    ASSERT(md_bitfield_validate(bf));

    if (bf->bits) {
        free_blocks(bf->alloc, bf->bits, num_blocks(bf->beg_bit, bf->end_bit));
    }
    MEMSET(bf, 0, sizeof(md_bitfield_t));
    return true;
}

bool md_bitfield_validate(const md_bitfield_t* bf) {
    return bf && bf->magic == MAGIC;
}

void md_bitfield_reserve_range(md_bitfield_t* bf, uint64_t beg_bit, uint64_t end_bit) {
    ASSERT(md_bitfield_validate(bf));
    ensure_range(bf, beg_bit, end_bit);
}

bool md_bitfield_empty(const md_bitfield_t* bf) {
    return bf->bits == NULL || (bf->beg_bit == bf->end_bit);
}

uint64_t md_bitfield_beg_bit(const md_bitfield_t* bf) {
    return bf->bits ? bf->beg_bit : 0;
}

uint64_t md_bitfield_end_bit(const md_bitfield_t* bf) {
    return bf->bits ? bf->end_bit : 0;
}

void md_bitfield_set_range(md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(bf));
    if (end <= beg) return;

    ensure_range(bf, beg, end);
    bit_set(u64_base_mut(bf), beg, end);
}

void md_bitfield_set_bit(md_bitfield_t* bf, uint64_t bit_idx) {
    ASSERT(md_bitfield_validate(bf));

    ensure_range(bf, bit_idx, bit_idx+1);
    block_set_bit(get_blk_ptr(bf, block_idx(bit_idx)), bit_idx - block_bit(bit_idx));
}

void md_bitfield_set_indices_u32(md_bitfield_t* bf, const uint32_t* indices, size_t num_indices) {
    ASSERT(md_bitfield_validate(bf));
    ASSERT(indices);    

    for (uint64_t i = 0; i < num_indices; ++i) {
        uint64_t bit_idx = indices[i];
        ensure_range(bf, bit_idx, bit_idx+1);
        bit_set((uint64_t*)bf->bits, bit_idx - block_bit(bf->beg_bit), 1);
    }
}

void md_bitfield_reset(md_bitfield_t* bf) {
    ASSERT(md_bitfield_validate(bf));
    if (bf->bits) {
        free_blocks(bf->alloc, bf->bits, num_blocks(bf->beg_bit, bf->end_bit));
        bf->bits = NULL;
        bf->beg_bit = 0;
        bf->end_bit = 0;
    }
}

void md_bitfield_clear(md_bitfield_t* bf) {
    ASSERT(md_bitfield_validate(bf));
    if (bf->bits) {
        MEMSET(bf->bits, 0, num_blocks(bf->beg_bit, bf->end_bit) * sizeof(block_t));
    }
}

void md_bitfield_clear_range(md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(bf));
        
    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    if (end-beg == 0) return;
    bit_clear(u64_base_mut(bf), beg, end);
}

void md_bitfield_clear_bit(md_bitfield_t* bf, uint64_t bit_idx) {
    ASSERT(md_bitfield_validate(bf));

    if (bf->beg_bit <= bit_idx && bit_idx < bf->end_bit) {
        bit_clear_idx(u64_base_mut(bf), bit_idx);
    }
}

void md_bitfield_or(md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b) {
    ASSERT(md_bitfield_validate(dst));
    ASSERT(md_bitfield_validate(src_a));
    ASSERT(md_bitfield_validate(src_b));

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    uint64_t beg_bit = MIN(src_a->beg_bit, src_b->beg_bit);
    uint64_t end_bit = MAX(src_a->end_bit, src_b->end_bit);

    fit_to_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_or(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_or_inplace(md_bitfield_t* a, const md_bitfield_t* b) {
    ASSERT(md_bitfield_validate(a));
    ASSERT(md_bitfield_validate(b));

    if (a == b) return;
    if (a->bits == 0) {
        md_bitfield_copy(a, b);
        return;
    }

    uint64_t beg_bit = MIN(a->beg_bit, b->beg_bit);
    uint64_t end_bit = MAX(a->end_bit, b->end_bit);
    ensure_range(a, beg_bit, end_bit);

    if (a->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_or(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_and(md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b) {
    ASSERT(md_bitfield_validate(dst));
    ASSERT(md_bitfield_validate(src_a));
    ASSERT(md_bitfield_validate(src_b));

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    uint64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    uint64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(dst);
        return;
    }

    fit_to_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_and(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_and_inplace(md_bitfield_t* a, const md_bitfield_t* b) {
    ASSERT(md_bitfield_validate(a));
    ASSERT(md_bitfield_validate(b));

    if (a == b) return;

    uint64_t beg_bit = MAX(a->beg_bit, b->beg_bit);
    uint64_t end_bit = MIN(a->end_bit, b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(a);
        return;
    }

    if (a->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_and(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_andnot(md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b) {
    ASSERT(md_bitfield_validate(dst));
    ASSERT(md_bitfield_validate(src_a));
    ASSERT(md_bitfield_validate(src_b));

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    if (src_a == src_b) {
        md_bitfield_clear(dst);
        return;
    }

    uint64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    uint64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);
    
    md_bitfield_copy(dst, src_a);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap, done.
        return;
    }

    ensure_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_andnot(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_andnot_inplace(md_bitfield_t* a, const md_bitfield_t* b) {
    ASSERT(md_bitfield_validate(a));
    ASSERT(md_bitfield_validate(b));

    if (a == b) {
        md_bitfield_clear(a);
        return;
    }

    uint64_t beg_bit = MAX(a->beg_bit, b->beg_bit);
    uint64_t end_bit = MIN(a->end_bit, b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap, done.
        return;
    }

    // The only place we need to modify is where the two ranges overlap.
    // ensure_range(a, beg_bit, end_bit);

    if (a->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_andnot(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_xor(md_bitfield_t* dst, const md_bitfield_t* src_a, const md_bitfield_t* src_b) {
    ASSERT(md_bitfield_validate(dst));
    ASSERT(md_bitfield_validate(src_a));
    ASSERT(md_bitfield_validate(src_b));

    ASSERT((dst != src_a && dst != src_b) && "dst is same as src_a or src_b, use inplace version instead!");

    uint64_t beg_bit = MAX(src_a->beg_bit, src_b->beg_bit);
    uint64_t end_bit = MIN(src_a->end_bit, src_b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(dst);
        return;
    }

    fit_to_range(dst, beg_bit, end_bit);
    if (dst->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(dst, i, block_xor(get_block(src_a, i), get_block(src_b, i)));
        }
    }
}

void md_bitfield_xor_inplace(md_bitfield_t* a, const md_bitfield_t* b) {
    ASSERT(md_bitfield_validate(a));
    ASSERT(md_bitfield_validate(b));

    if (a == b) return;

    uint64_t beg_bit = MAX(a->beg_bit, b->beg_bit);
    uint64_t end_bit = MIN(a->end_bit, b->end_bit);

    if (end_bit <= beg_bit) {
        // Ranges do not overlap
        md_bitfield_clear(a);
        return;
    }

    ensure_range(a, beg_bit, end_bit);
    if (a->bits) {
        for (uint64_t i = block_idx(beg_bit); i <= block_idx(end_bit); ++i) {
            set_block(a, i, block_xor(get_block(a, i), get_block(b, i)));
        }
    }
}

void md_bitfield_not(md_bitfield_t* dst, const md_bitfield_t* src, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(dst));
    ASSERT(md_bitfield_validate(src));

    ASSERT((dst != src) && "dst is same as src, use inplace version instead!");

    fit_to_range(dst, beg, end);

    if (dst->bits) {
        for (uint64_t i = block_idx(beg); i <= block_idx(end); ++i) {
            set_block(dst, i, block_not(get_block(src, i)));
        }
    }
}

void md_bitfield_not_inplace(md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(bf));

    ensure_range(bf, beg, end);

    if (bf->bits) {
        for (uint64_t i = block_idx(beg); i <= block_idx(end); ++i) {
            set_block(bf, i, block_not(get_block(bf, i)));
        }
    }
}

void md_bitfield_copy(md_bitfield_t* dst, const md_bitfield_t* src) {
    ASSERT(md_bitfield_validate(dst));
    ASSERT(md_bitfield_validate(src));
    if (dst == src) return;

    fit_to_range(dst, src->beg_bit, src->end_bit);
    if (dst->bits) {
        MEMCPY(dst->bits, src->bits, num_blocks(dst->beg_bit, dst->end_bit) * sizeof(block_t));
    }
    dst->flags = src->flags;
}

// Counts the number of bits set
size_t md_bitfield_popcount(const md_bitfield_t* bf) {
    ASSERT(md_bitfield_validate(bf));

    if (bf->bits == NULL || bf->beg_bit == bf->end_bit) return 0;
    return bit_count(u64_base(bf), bf->beg_bit, bf->end_bit);
}

size_t md_bitfield_popcount_range(const md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(beg <= end);
    ASSERT(md_bitfield_validate(bf));

    if (bf->bits == NULL) return 0;
    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    if (beg == end) return 0;

    return bit_count(u64_base(bf), beg, end);
}

// Test if single bit in field is set
bool md_bitfield_test_bit(const md_bitfield_t* bf, uint64_t idx) {
    ASSERT(md_bitfield_validate(bf));

    if (idx < bf->beg_bit || bf->end_bit <= idx) return false;
    return block_test_bit(get_block(bf, block_idx(idx)), idx - block_bit(idx));
   
}

bool md_bitfield_test_all (const md_bitfield_t* bf) {
    return md_bitfield_test_all_range(bf, bf->beg_bit, bf->end_bit);
}

bool md_bitfield_test_all_range (const md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(bf));

    if (end < bf->beg_bit || bf->end_bit < beg) return false;
    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    if (end <= beg) return false;
    
    uint64_t count = md_bitfield_popcount_range(bf, beg, end);
    return count == (end - beg);
}

bool md_bitfield_test_any (const md_bitfield_t* bf) {
	return md_bitfield_test_any_range(bf, bf->beg_bit, bf->end_bit);
}

bool md_bitfield_test_any_range (const md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(bf));
    return md_bitfield_popcount_range(bf, beg, end) > 0;
}

// Test if bitfields are equivalent
//bool md_bitfield_cmp        (const uint64_t* src_a, const uint64_t* src_b);
//bool md_bitfield_cmp_range  (const uint64_t* src_a, const uint64_t* src_b, int64_t beg, int64_t end);

// Bit scan forward, finds the first bit set from the given offset.
// Returns 0 if no bit is found, otherwise it returns the offset to that bit (indexing starts at 1, posix convention)
uint64_t md_bitfield_scan(const md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    ASSERT(md_bitfield_validate(bf));

    beg = CLAMP(beg, bf->beg_bit, bf->end_bit);
    end = CLAMP(end, bf->beg_bit, bf->end_bit);
    if (beg == end)
        return 0;
    
    return bit_scan_forward(u64_base(bf), beg, end);
}

md_bitfield_iter_t md_bitfield_iter_create(const md_bitfield_t* bf) {
    md_bitfield_iter_t iter = {
        .bf = bf,
        .idx = bf->beg_bit,
        .beg_bit = bf->beg_bit,
        .end_bit = bf->end_bit,
    };
    return iter;
}

md_bitfield_iter_t md_bitfield_iter_range_create(const md_bitfield_t* bf, uint64_t beg, uint64_t end) {
    md_bitfield_iter_t iter = {
        .bf = bf,
        .idx = bf->beg_bit,
        .beg_bit = MAX((uint32_t)beg, bf->beg_bit),
        .end_bit = MIN((uint32_t)end, bf->end_bit),
    };
    return iter;
}

bool md_bitfield_iter_next(md_bitfield_iter_t* it) {
    ASSERT(it);
    ASSERT(md_bitfield_validate(it->bf));

    if (it->idx >= it->end_bit) return false;
    
    const uint64_t res = bit_scan_forward(u64_base(it->bf), it->idx, it->end_bit);
    if (res) {
        it->idx = res;
        return true;
    }
    
    it->idx = it->end_bit;
    return false;
}

void md_bitfield_iter_skip_to_idx(md_bitfield_iter_t* it, uint64_t idx) {
    ASSERT(it);
    it->idx = idx;
}

bool md_bitfield_get_range(uint64_t* first_idx, uint64_t* last_idx, const md_bitfield_t* bf) {
    ASSERT(first_idx);
    ASSERT(last_idx);
    const uint64_t first = bit_scan_forward(u64_base(bf), bf->beg_bit, bf->end_bit);
    if (first) {
        const uint64_t last = bit_scan_reverse(u64_base(bf), bf->beg_bit, bf->end_bit);
        *first_idx = first - 1;
        *last_idx  = last - 1;
        return true;
    }
    return false;
}

static inline uint64_t get_u64(const md_bitfield_t* bf, uint64_t idx_u64) {
    if (!bf->bits) return 0;

    const uint64_t beg_u64 = bf->beg_bit / 64;
    const uint64_t end_u64 = bf->end_bit / 64;

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

bool md_bitfield_extract_u64(uint64_t* dst_ptr, uint64_t num_bits, const md_bitfield_t* src) {
    ASSERT(dst_ptr);
    ASSERT(md_bitfield_validate(src));
    
    if (num_bits == 0) return true;

    const int64_t count = (num_bits + 63) / 64;
    for (int64_t i = 0; i < count; ++i) {
        dst_ptr[i] = get_u64(src, i);
    }
    
    return true;
}

size_t md_bitfield_iter_extract_indices(int32_t* buf, size_t cap, md_bitfield_iter_t it) {
    ASSERT(md_bitfield_validate(it.bf));
    ASSERT(cap >= 0);

    if (!buf || cap == 0) return 0;

    size_t len = 0;
    while (len < cap && md_bitfield_iter_next(&it)) {
        buf[len++] = (int32_t)md_bitfield_iter_idx(&it);
    }

    return len;
}

/*
uint32_t* md_bitfield_extract_bits_u32(const md_bitfield_t* bf, struct md_allocator_i* alloc) {
    ASSERT(bf);
    ASSERT(alloc);
    ASSERT(md_bitfield_validate(bf));

    uint32_t* bits = 0;

    uint64_t block_count = num_blocks(bf->beg_bit, bf->end_bit);
    if (block_count) {
        for (uint64_t i = block_idx(bf->beg_bit); i <= block_idx(bf->end_bit); ++i) {
            block_t block = get_block(bf, i);
            for (uint64_t j = 0; j < ARRAY_SIZE(block.u32); ++j) {
                md_array_push(bits, block.u32[j], alloc);
            }
        }
    }

    return bits;
}
*/

#define BLOCK_IDX_FLAG_ALL_SET (0x8000)

uint32_t* get_non_empty_block_indices(const md_bitfield_t* bf, md_allocator_i* alloc) {
    uint32_t* indices = 0;

    uint64_t beg_blk_idx = block_idx(bf->beg_bit);
    uint64_t end_blk_idx = block_idx(bf->end_bit);
    for (uint64_t i = beg_blk_idx; i <= end_blk_idx; ++i) {
        block_t blk = get_block(bf, i);
        if (i == beg_blk_idx) block_and(blk, block_mask_hi(bf->beg_bit & (BITS_PER_BLOCK-1)));
        if (i == beg_blk_idx) block_and(blk, block_mask_lo(bf->end_bit & (BITS_PER_BLOCK-1)));

        bool empty = true;
        for (size_t j = 0; j < ARRAY_SIZE(blk.u64); ++j) {
            if (blk.u64[j] != 0) empty = false;
        }
        if (empty) continue;
        md_array_push(indices, (uint32_t)i, alloc);
    }

    return indices;
}

uint16_t* get_serialization_block_indices(const md_bitfield_t* bf, md_allocator_i* alloc) {
    uint16_t* indices = 0;

    uint64_t beg_blk_idx = block_idx(bf->beg_bit);
    uint64_t end_blk_idx = block_idx(bf->end_bit);
    for (uint64_t i = beg_blk_idx; i <= end_blk_idx; ++i) {
        block_t blk = get_block(bf, i);
        if (i == beg_blk_idx) block_and(blk, block_mask_hi(bf->beg_bit & (BITS_PER_BLOCK-1)));
        if (i == beg_blk_idx) block_and(blk, block_mask_lo(bf->end_bit & (BITS_PER_BLOCK-1)));

        bool empty = true;
        bool all_set = true;
        for (size_t j = 0; j < ARRAY_SIZE(blk.u64); ++j) {
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
size_t md_bitfield_serialize_size_in_bytes(const md_bitfield_t* bf) {
    size_t temp_pos = md_temp_get_pos();
    size_t size = sizeof(uint16_t);
    uint16_t* indices = get_serialization_block_indices(bf, md_get_temp_allocator());
    for (size_t i = 0; i < md_array_size(indices); ++i) {
        if (indices[i] & BLOCK_IDX_FLAG_ALL_SET) continue;
        size += sizeof(indices[i]) + sizeof(block_t);
    }

    md_temp_set_pos_back(temp_pos);
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
size_t md_bitfield_serialize(void* dst, const md_bitfield_t* bf) {
    size_t temp_pos = md_temp_get_pos();
    md_allocator_i* alloc = md_get_temp_allocator();

    uint16_t* indices = get_serialization_block_indices(bf, alloc);
    uint16_t* data = 0;

    uint16_t block_count = (uint16_t)md_array_size(indices);
    uint64_t beg_blk = block_idx(bf->beg_bit);
    uint64_t end_blk = block_idx(bf->end_bit);
    
    md_array_push(data, block_count, alloc);
    md_array_push_array(data, indices, block_count, alloc);

    // Add actual block data from the blocks (if not all bits set within block)
    for (uint64_t i = 0; i < block_count; ++i) {
        uint64_t blk_idx = data[1 + i];
        bool    all_set = data[1 + i] & BLOCK_IDX_FLAG_ALL_SET;

        if (!all_set) {
            block_t blk = get_block(bf, blk_idx);
            if (blk_idx == beg_blk) block_and(blk, block_mask_hi(bf->beg_bit & (BITS_PER_BLOCK-1)));
            if (blk_idx == end_blk) block_and(blk, block_mask_lo(bf->end_bit & (BITS_PER_BLOCK-1)));
            md_array_push_array(data, blk.u16, ARRAY_SIZE(blk.u16), alloc);
        }
    }

    int lz_bytes = fastlz_compress_level(2, data, (int)md_array_bytes(data), dst);

    md_temp_set_pos_back(temp_pos);

    return (size_t)lz_bytes;
}

bool md_bitfield_deserialize(md_bitfield_t* bf, const void* src, size_t num_bytes) {
    ASSERT(bf);
    ASSERT(md_bitfield_validate(bf));
    ASSERT(src);

    // Temporary buffer to decompress into.
    // The current compression is no where near 100:1 ratio.
    const size_t mem_bytes = num_bytes * 100;
    void* mem = md_alloc(md_get_heap_allocator(), mem_bytes);

    if (!mem) {
        MD_LOG_ERROR("Failed to allocate temporary data for deserialization");
        return false;
    }

    bool result = false;

    int size = fastlz_decompress(src, (int)num_bytes, mem, (int)mem_bytes);    
    if (size == 0) {
        MD_LOG_ERROR("Failed, to decompress bitfield.");
        goto done;
    }

    const uint16_t* data = (const uint16_t*)mem;
    size_t block_count = data[0];
    const uint16_t* block_indices = data + 1;
    const uint8_t*  block_data = (const uint8_t*)(block_indices + block_count);

    if (block_count == 0) {
        MD_LOG_ERROR("Block count was zero");
        goto done;
    }

    size_t beg_blk_idx = block_indices[0];
    size_t end_blk_idx = block_indices[block_count - 1];

    // Allocate the blocks
    
    fit_to_range(bf, beg_blk_idx * BITS_PER_BLOCK, end_blk_idx * BITS_PER_BLOCK + (BITS_PER_BLOCK-1));

    size_t nblocks = num_blocks(bf->beg_bit, bf->end_bit);
    MEMSET(bf->bits, 0, nblocks * sizeof(block_t));

    // Fetch block_data and store
    size_t src_offset = 0;
    for (size_t i = 0; i < block_count; ++i) {
        uint16_t blk_idx = block_indices[i];
        if (blk_idx & BLOCK_IDX_FLAG_ALL_SET) {
            blk_idx &= ~BLOCK_IDX_FLAG_ALL_SET;
            block_t* blk_ptr = get_blk_ptr(bf, blk_idx);
            MEMSET(blk_ptr, 0xFFFFFFFF, sizeof(block_t));
        } else {
            block_t* blk_ptr = get_blk_ptr(bf, blk_idx);
            MEMCPY(blk_ptr, block_data + src_offset, sizeof(block_t));
            src_offset += sizeof(block_t);
        }
    }

    // Now we have the data, compute the true beg_bit and end_bit
    uint64_t beg_bit = block_scan_forward(get_block(bf, beg_blk_idx));
    uint64_t end_bit = block_scan_reverse(get_block(bf, end_blk_idx));

    bf->beg_bit = (uint32_t)(beg_blk_idx * BITS_PER_BLOCK + (beg_bit ? beg_bit - 1 : 0));
    bf->end_bit = (uint32_t)(end_blk_idx * BITS_PER_BLOCK + (end_bit ? end_bit : 0));

    result = true;
done:
    if (mem) {
        md_free(md_get_heap_allocator(), mem, mem_bytes);
    }

    return result;
}

uint64_t md_bitfield_hash64(const md_bitfield_t* bf, uint64_t seed) {
    ASSERT(bf);
    if (bf->beg_bit == bf->end_bit) return seed;
    const block_t* ptr = (const block_t*)bf->bits;
    const size_t num_blk = num_blocks(bf->beg_bit, bf->end_bit);
    return md_hash64(ptr, sizeof(block_t) * num_blk, seed);
}
