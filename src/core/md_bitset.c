#include <core/md_bitset.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_hash.h>

#include <core/md_arena_allocator.h>
#include <core/md_log.h>

#include <fastlz.h>
#include <string.h>

STATIC_ASSERT(sizeof(md_bitset_gpu_t) == 16, "GPU descriptor must be 16 bytes");

#define MAX_INDEX MD_BITSET_MAX_INDEX

static inline uint32_t min_u32(uint32_t a, uint32_t b) { return a < b ? a : b; }
static inline uint32_t max_u32(uint32_t a, uint32_t b) { return a > b ? a : b; }

// Bits [lo, hi) of a 64-bit word, 0 <= lo < hi <= 64
static inline uint64_t bits_between(uint32_t lo, uint32_t hi) {
    const uint64_t a = ~0ull << lo;
    const uint64_t b = hi >= 64 ? ~0ull : ~(~0ull << hi);
    return a & b;
}

// The bits of word g (indices [64 g, 64 g + 64)) which lie within [beg, end), beg < end
static inline uint64_t word_range_mask(uint32_t g, uint32_t beg, uint32_t end) {
    const uint64_t wb = (uint64_t)g << 6;
    if (wb + 64 <= beg || wb >= end) return 0;
    const uint32_t lo = beg > wb ? (uint32_t)(beg - wb) : 0;
    const uint32_t hi = (uint64_t)end < wb + 64 ? (uint32_t)(end - wb) : 64;
    return bits_between(lo, hi);
}

static inline uint64_t* alloc_words(md_allocator_i* alloc, size_t num_words) {
    ASSERT(alloc);
    uint64_t* p = (uint64_t*)md_alloc(alloc, num_words * sizeof(uint64_t));
    ASSERT(p);
    ASSERT(((uintptr_t)p & 7) == 0);
    return p;
}

// ---------------------------------------------------------------------------------------------------------
// Word level evaluation of a binary operation

typedef enum {
    OP_AND,
    OP_OR,
    OP_ANDNOT,
    OP_XOR,
} op_t;

typedef struct {
    uint64_t count;
    uint64_t first_g;   // global word index (index >> 6) of the first non-zero result word
    uint64_t last_g;
    uint64_t first_w;
    uint64_t last_w;
} stats_t;

// Global word indices where an operand changes behaviour. Between two consecutive breakpoints each operand
// is either a constant word (0 or ~0) or a contiguous slice of its words; RUN edge words, which are
// partial, are isolated into single word intervals.
static inline int push_breakpoints(uint64_t* bp, int n, md_bitset_t s) {
    if (s.beg >= s.end) return n;
    const uint64_t gb = s.beg >> 6;
    const uint64_t ge = ((uint64_t)s.end + 63) >> 6;
    if (s.words) {
        bp[n++] = gb;
        bp[n++] = ge;
    } else {
        const uint64_t gl = ((uint64_t)s.end - 1) >> 6;
        bp[n++] = gb;
        bp[n++] = gb + 1;
        bp[n++] = gl;
        bp[n++] = gl + 1;
    }
    return n;
}

// Pointer and stride for reading operand s over the word interval [g, g + n). A constant word is served
// with stride 0 from *scratch.
static inline const uint64_t* operand_slice(size_t* stride, uint64_t* scratch, md_bitset_t s, uint64_t g, uint64_t n) {
    if (s.words && s.beg < s.end) {
        const uint64_t g0 = s.beg >> 6;
        const uint64_t g1 = ((uint64_t)s.end + 63) >> 6;
        if (g0 <= g && g + n <= g1) {
            *stride = 1;
            return s.words + (g - g0);
        }
    }
    *scratch = md_bitset_word(s, (uint32_t)(g << 6));
    *stride = 0;
    return scratch;
}

#define OP_LOOP(EXPR)                                               \
    for (uint64_t k = 0; k < n; ++k) {                              \
        const uint64_t x = pa[k * sa];                              \
        const uint64_t y = pb[k * sb];                              \
        const uint64_t w = (EXPR);                                  \
        if (out) out[g - out_g + k] = w;                            \
        if (w) {                                                    \
            if (st->count == 0) { st->first_g = g + k; st->first_w = w; } \
            st->last_g = g + k; st->last_w = w;                     \
            st->count += popcnt64(w);                               \
        }                                                           \
    }

// Evaluate op over the global words [G0, G1). Accumulates stats; if out is given, writes word g to out[g - out_g].
static void eval_op(stats_t* st, uint64_t* out, uint64_t out_g, md_bitset_t a, md_bitset_t b, op_t op, uint64_t G0, uint64_t G1) {
    uint64_t bp[10];
    int nbp = 0;
    nbp = push_breakpoints(bp, nbp, a);
    nbp = push_breakpoints(bp, nbp, b);

    // Insertion sort, at most 8 entries
    for (int i = 1; i < nbp; ++i) {
        uint64_t v = bp[i];
        int j = i - 1;
        while (j >= 0 && bp[j] > v) { bp[j + 1] = bp[j]; --j; }
        bp[j + 1] = v;
    }

    uint64_t g = G0;
    int bi = 0;
    while (g < G1) {
        while (bi < nbp && bp[bi] <= g) ++bi;
        const uint64_t e = (bi < nbp && bp[bi] < G1) ? bp[bi] : G1;
        const uint64_t n = e - g;

        uint64_t ca, cb;
        size_t sa, sb;
        const uint64_t* pa = operand_slice(&sa, &ca, a, g, n);
        const uint64_t* pb = operand_slice(&sb, &cb, b, g, n);

        switch (op) {
        case OP_AND:    OP_LOOP(x & y);  break;
        case OP_OR:     OP_LOOP(x | y);  break;
        case OP_ANDNOT: OP_LOOP(x & ~y); break;
        case OP_XOR:    OP_LOOP(x ^ y);  break;
        default: ASSERT(false); break;
        }
        g = e;
    }
}

#undef OP_LOOP

static md_bitset_t binary_op(md_bitset_t a, md_bitset_t b, op_t op, uint32_t dom_beg, uint32_t dom_end, md_allocator_i* alloc) {
    if (dom_beg >= dom_end) return md_bitset_empty_set();

    const uint64_t G0 = dom_beg >> 6;
    const uint64_t G1 = ((uint64_t)dom_end + 63) >> 6;

    stats_t st = {0};
    eval_op(&st, NULL, 0, a, b, op, G0, G1);
    if (st.count == 0) return md_bitset_empty_set();

    md_bitset_t r;
    r.beg = (uint32_t)((st.first_g << 6) + ctz64(st.first_w));
    r.end = (uint32_t)((st.last_g << 6) + 64 - clz64(st.last_w));
    r.words = NULL;
    if (st.count == (uint64_t)(r.end - r.beg)) return r;  // dense: RUN

    const uint64_t nw = st.last_g - st.first_g + 1;
    uint64_t* words = alloc_words(alloc, nw);
    stats_t st2 = {0};
    eval_op(&st2, words, st.first_g, a, b, op, st.first_g, st.last_g + 1);
    ASSERT(st2.count == st.count);
    r.words = words;
    return r;
}

// ---------------------------------------------------------------------------------------------------------
// Producers

md_bitset_t md_bitset_from_words(const uint64_t* words, size_t num_words, uint32_t base, md_allocator_i* alloc) {
    ASSERT((base & 63) == 0);
    // Every index the words can hold must be below MAX_INDEX, or end (one past the last member) can wrap to 0
    ASSERT(((uint64_t)base + (uint64_t)num_words * 64) <= (uint64_t)MAX_INDEX);
    if (!words || !num_words) return md_bitset_empty_set();

    size_t first = SIZE_MAX, last = 0;
    uint64_t count = 0;
    for (size_t k = 0; k < num_words; ++k) {
        if (words[k]) {
            if (first == SIZE_MAX) first = k;
            last = k;
            count += popcnt64(words[k]);
        }
    }
    if (!count) return md_bitset_empty_set();

    md_bitset_t r;
    r.beg = base + (uint32_t)(first * 64 + ctz64(words[first]));
    r.end = base + (uint32_t)(last  * 64 + 64 - clz64(words[last]));
    r.words = NULL;
    if (count == (uint64_t)(r.end - r.beg)) return r;

    const size_t nw = last - first + 1;
    uint64_t* dst = alloc_words(alloc, nw);
    MEMCPY(dst, words + first, nw * sizeof(uint64_t));
    r.words = dst;
    return r;
}

md_bitset_t md_bitset_from_indices(const uint32_t* indices, size_t num_indices, md_allocator_i* alloc) {
    if (!indices || !num_indices) return md_bitset_empty_set();

    uint32_t lo = UINT32_MAX, hi = 0;
    for (size_t i = 0; i < num_indices; ++i) {
        lo = min_u32(lo, indices[i]);
        hi = max_u32(hi, indices[i]);
    }
    ASSERT(hi < MAX_INDEX);

    const uint32_t base = lo & ~63u;
    const size_t nw = ((size_t)hi >> 6) - ((size_t)lo >> 6) + 1;

    // Scratch from a temp arena other than alloc: if alloc is itself a temp arena, the result must not land
    // inside the scope we rewind.
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    uint64_t* scratch = (uint64_t*)md_temp_alloc_zero(temp, nw * sizeof(uint64_t));
    for (size_t i = 0; i < num_indices; ++i) {
        const uint32_t k = indices[i] - base;
        scratch[k >> 6] |= 1ull << (k & 63);
    }
    md_bitset_t r = md_bitset_from_words(scratch, nw, base, alloc);
    md_temp_end(temp);
    return r;
}

md_bitset_t md_bitset_copy(md_bitset_t s, md_allocator_i* alloc) {
    if (s.beg >= s.end) return md_bitset_empty_set();
    if (!s.words) return s;
    // Through from_words: a copy of a view of a mask comes out tidy
    return md_bitset_from_words(s.words, md_bitset_num_words(s), s.beg & ~63u, alloc);
}

md_bitset_t md_bitset_and(md_bitset_t a, md_bitset_t b, md_allocator_i* alloc) {
    if (md_bitset_empty(a) || md_bitset_empty(b)) return md_bitset_empty_set();
    const uint32_t beg = max_u32(a.beg, b.beg);
    const uint32_t end = min_u32(a.end, b.end);
    if (beg >= end) return md_bitset_empty_set();
    if (!a.words && !b.words) return md_bitset_range(beg, end);
    return binary_op(a, b, OP_AND, beg, end, alloc);
}

md_bitset_t md_bitset_or(md_bitset_t a, md_bitset_t b, md_allocator_i* alloc) {
    if (md_bitset_empty(a)) return md_bitset_copy(b, alloc);
    if (md_bitset_empty(b)) return md_bitset_copy(a, alloc);
    if (!a.words && !b.words && a.beg <= b.end && b.beg <= a.end) {
        return md_bitset_range(min_u32(a.beg, b.beg), max_u32(a.end, b.end));
    }
    return binary_op(a, b, OP_OR, min_u32(a.beg, b.beg), max_u32(a.end, b.end), alloc);
}

md_bitset_t md_bitset_andnot(md_bitset_t a, md_bitset_t b, md_allocator_i* alloc) {
    if (md_bitset_empty(a)) return md_bitset_empty_set();
    if (md_bitset_empty(b) || b.end <= a.beg || a.end <= b.beg) return md_bitset_copy(a, alloc);
    if (!b.words) {
        if (b.beg <= a.beg && a.end <= b.end) return md_bitset_empty_set();
        if (!a.words) {
            if (b.beg <= a.beg) return md_bitset_range(b.end, a.end);
            if (a.end <= b.end) return md_bitset_range(a.beg, b.beg);
            // b punches a hole in a: needs words
        }
    }
    return binary_op(a, b, OP_ANDNOT, a.beg, a.end, alloc);
}

md_bitset_t md_bitset_xor(md_bitset_t a, md_bitset_t b, md_allocator_i* alloc) {
    if (md_bitset_empty(a)) return md_bitset_copy(b, alloc);
    if (md_bitset_empty(b)) return md_bitset_copy(a, alloc);
    return binary_op(a, b, OP_XOR, min_u32(a.beg, b.beg), max_u32(a.end, b.end), alloc);
}

md_bitset_t md_bitset_not(md_bitset_t s, uint32_t beg, uint32_t end, md_allocator_i* alloc) {
    ASSERT(end <= MAX_INDEX);
    return md_bitset_andnot(md_bitset_range(beg, end), s, alloc);
}

void md_bitset_free(md_bitset_t* s, md_allocator_i* alloc) {
    ASSERT(s);
    if (s->words && s->beg < s->end) {
        ASSERT(alloc);
        md_free(alloc, (void*)s->words, md_bitset_num_words(*s) * sizeof(uint64_t));
    }
    *s = md_bitset_empty_set();
}

void md_bitset_replace(md_bitset_t* dst, md_bitset_t value, md_allocator_i* alloc) {
    ASSERT(dst);
    if (dst->words && dst->words == value.words) {  // same words: nothing to release
        *dst = value;
        return;
    }
    md_bitset_free(dst, alloc);
    *dst = value;
}

void md_bitset_assign(md_bitset_t* dst, md_bitset_t src, md_allocator_i* alloc) {
    ASSERT(dst);
    if (dst->words && dst->words == src.words && dst->beg == src.beg && dst->end == src.end) return;  // already this value
    md_bitset_replace(dst, md_bitset_copy(src, alloc), alloc);
}

// ---------------------------------------------------------------------------------------------------------
// Readers

uint64_t md_bitset_count(md_bitset_t s) {
    if (s.beg >= s.end) return 0;
    if (!s.words) return s.end - s.beg;
    const size_t nw = md_bitset_num_words(s);
    uint64_t count = 0;
    for (size_t k = 0; k < nw; ++k) count += popcnt64(s.words[k]);
    return count;
}

uint64_t md_bitset_count_range(md_bitset_t s, uint32_t beg, uint32_t end) {
    beg = max_u32(beg, s.beg);
    end = min_u32(end, s.end);
    if (beg >= end) return 0;
    if (!s.words) return end - beg;

    const uint32_t wb = s.beg & ~63u;
    const uint32_t k0 = (beg - wb) >> 6;
    const uint32_t k1 = (end - 1 - wb) >> 6;
    const uint64_t lo = ~0ull << ((beg - wb) & 63);
    const uint64_t hi = ~0ull >> (63 - ((end - 1 - wb) & 63));
    if (k0 == k1) return popcnt64(s.words[k0] & lo & hi);

    uint64_t count = popcnt64(s.words[k0] & lo) + popcnt64(s.words[k1] & hi);
    for (uint32_t k = k0 + 1; k < k1; ++k) count += popcnt64(s.words[k]);
    return count;
}

bool md_bitset_test_any_range(md_bitset_t s, uint32_t beg, uint32_t end) {
    beg = max_u32(beg, s.beg);
    end = min_u32(end, s.end);
    if (beg >= end) return false;
    if (!s.words) return true;
    // beg and end - 1 of a set are members, so a range covering either is a hit without reading words
    if (beg == s.beg || end == s.end) return true;

    const uint32_t wb = s.beg & ~63u;
    const uint32_t k0 = (beg - wb) >> 6;
    const uint32_t k1 = (end - 1 - wb) >> 6;
    const uint64_t lo = ~0ull << ((beg - wb) & 63);
    const uint64_t hi = ~0ull >> (63 - ((end - 1 - wb) & 63));
    if (k0 == k1) return (s.words[k0] & lo & hi) != 0;
    if (s.words[k0] & lo) return true;
    if (s.words[k1] & hi) return true;
    for (uint32_t k = k0 + 1; k < k1; ++k) {
        if (s.words[k]) return true;
    }
    return false;
}

bool md_bitset_test_all_range(md_bitset_t s, uint32_t beg, uint32_t end) {
    if (beg >= end) return true;
    if (beg < s.beg || s.end < end) return false;
    if (!s.words) return true;
    return md_bitset_count_range(s, beg, end) == (uint64_t)(end - beg);
}

bool md_bitset_equal(md_bitset_t a, md_bitset_t b) {
    const bool ea = a.beg >= a.end, eb = b.beg >= b.end;
    if (ea || eb) return ea == eb;
    if (a.beg != b.beg || a.end != b.end) return false;
    if (!a.words && !b.words) return true;
    if (!a.words || !b.words) {
        // A RUN equals a BITS set exactly when every index of the span is a member of it
        const md_bitset_t bits = a.words ? a : b;
        return md_bitset_count(bits) == (uint64_t)(bits.end - bits.beg);
    }
    if (a.words == b.words) return true;
    return memcmp(a.words, b.words, md_bitset_num_words(a) * sizeof(uint64_t)) == 0;
}

uint64_t md_bitset_hash64(md_bitset_t s, uint64_t seed) {
    if (s.beg >= s.end) s = md_bitset_empty_set();
    const uint32_t range[2] = {s.beg, s.end};
    uint64_t h = md_hash64(range, sizeof(range), seed);
    // A full span hashes as its range only, whatever its form
    if (s.words && md_bitset_count(s) != (uint64_t)(s.end - s.beg)) {
        h = md_hash64(s.words, md_bitset_num_words(s) * sizeof(uint64_t), h);
    }
    return h;
}

size_t md_bitset_extract_indices(uint32_t* out, size_t cap, md_bitset_t s) {
    size_t n = 0;
    md_bitset_iter_t it = md_bitset_iter(s);
    while (n < cap && md_bitset_iter_next(&it)) out[n++] = it.idx;
    return n;
}

size_t md_bitset_iter_extract_indices(int32_t* out, size_t cap, md_bitset_iter_t it) {
    size_t n = 0;
    while (n < cap && md_bitset_iter_next(&it)) out[n++] = (int32_t)it.idx;
    return n;
}

bool md_bitset_validate(md_bitset_t s, bool tidy) {
    if (s.beg >= s.end) return s.beg == s.end && (!tidy || (s.beg == 0 && s.words == NULL));
    if (s.end > MAX_INDEX) return false;
    if (!s.words) return true;

    const size_t nw = md_bitset_num_words(s);
    if (tidy) {
        uint64_t count = 0;
        for (size_t k = 0; k < nw; ++k) count += popcnt64(s.words[k]);
        if (count == (uint64_t)(s.end - s.beg)) return false;  // should have been a RUN
    }

    const uint32_t lo = s.beg & 63;
    const uint32_t hi = (s.end - 1) & 63;
    if (((s.words[0] >> lo) & 1) == 0) return false;                       // beg is a member
    if (lo && (s.words[0] & ~(~0ull << lo))) return false;                 // nothing below beg
    if (((s.words[nw - 1] >> hi) & 1) == 0) return false;                  // end - 1 is a member
    if (hi != 63 && (s.words[nw - 1] & (~0ull << (hi + 1)))) return false; // nothing at or above end
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// GPU packing

size_t md_bitset_pack_num_words(const md_bitset_t* sets, size_t num_sets) {
    size_t n = 0;
    for (size_t i = 0; i < num_sets; ++i) n += md_bitset_num_words(sets[i]);
    return n;
}

void md_bitset_pack(md_bitset_gpu_t* out_sets, uint64_t* out_words, const md_bitset_t* sets, size_t num_sets) {
    ASSERT(out_sets || !num_sets);
    size_t off = 0;
    for (size_t i = 0; i < num_sets; ++i) {
        md_bitset_t s = sets[i];
        if (s.beg >= s.end) s = md_bitset_empty_set();
        md_bitset_gpu_t g = {s.beg, s.end, MD_BITSET_RUN, (uint32_t)md_bitset_count(s)};
        if (s.words) {
            const size_t nw = md_bitset_num_words(s);
            ASSERT(off + nw < MD_BITSET_RUN);
            MEMCPY(out_words + off, s.words, nw * sizeof(uint64_t));
            g.word_off = (uint32_t)off;
            off += nw;
        }
        out_sets[i] = g;
    }
}

// ---------------------------------------------------------------------------------------------------------
// Serialization. Two formats, both of 512-bit blocks, where empty blocks are omitted so that a sparse set stays
// small however far apart its members are:
//
// V1, byte compatible with the format of the former md_bitfield, written whenever the set fits it (every index
// below 2^24), so files stay readable by versions which know only it:
//   FastLZ( u16 block_count, u16 block_index[block_count], u8 block_data[64][#blocks not all set] )
// An index with the top bit set is a block with every bit set, which has no data. The 15 bits of the index are
// what limit it to 2^24 indices.
//
// V2, for the sets which do not fit V1:
//   u8 magic[4] = "MBS2", u32 raw_size, FastLZ( u32 entry_count, entries, u8 block_data[64][#blocks not all set] )
// An entry is a u32 block index; with the top bit set it is followed by a u32 count, of consecutive blocks with
// every bit set from that one on (so a range of 100M indices is one entry, not 200k).
// A FastLZ stream starts with a byte below 0x40 (its level in the top 3 bits, 1 or 2), so 'M' (0x4D) tells the
// formats apart.

#define BLOCK_BITS             512u
#define BLOCK_WORDS            8u
#define SERIAL_MIN_INPUT_SIZE  16u

#define V1_FLAG_ALL_SET        0x8000u
#define V1_IDX_MASK            0x7FFFu
#define V1_MAX_RAW_SIZE        (sizeof(uint16_t) + ((size_t)V1_IDX_MASK + 1) * (sizeof(uint16_t) + BLOCK_WORDS * sizeof(uint64_t)))

#define V2_FLAG_ALL_SET        0x80000000u
#define V2_IDX_MASK            0x7FFFFFFFu
#define V2_MAGIC_SIZE          4u
#define V2_HEADER_SIZE         (V2_MAGIC_SIZE + sizeof(uint32_t))
static const uint8_t v2_magic[V2_MAGIC_SIZE] = {'M', 'B', 'S', '2'};

static inline size_t fastlz_bound(size_t input_size) {
    size_t b = input_size + (input_size + 19) / 20 + 16;
    return b > 66 ? b : 66;
}

// The layout of the blocks in V1
typedef struct {
    size_t   block_count;
    size_t   block_data_bytes;
    size_t   input_size;
    uint32_t beg_blk;
    uint32_t end_blk;  // inclusive
} serial_info_t;

static inline void get_block(uint64_t out[BLOCK_WORDS], md_bitset_t s, uint32_t blk) {
    for (uint32_t j = 0; j < BLOCK_WORDS; ++j) out[j] = md_bitset_word(s, blk * BLOCK_BITS + j * 64);
}

static inline bool block_empty(const uint64_t blk[BLOCK_WORDS]) {
    uint64_t x = 0;
    for (uint32_t j = 0; j < BLOCK_WORDS; ++j) x |= blk[j];
    return x == 0;
}

static inline bool block_full(const uint64_t blk[BLOCK_WORDS]) {
    uint64_t x = ~0ull;
    for (uint32_t j = 0; j < BLOCK_WORDS; ++j) x &= blk[j];
    return x == ~0ull;
}

// Whether the set can be written in V1
static inline bool fits_v1(md_bitset_t s) {
    return md_bitset_empty(s) || (s.end - 1) / BLOCK_BITS <= V1_IDX_MASK;
}

static void serial_info(serial_info_t* info, md_bitset_t s) {
    MEMSET(info, 0, sizeof(*info));
    if (md_bitset_empty(s)) return;
    info->beg_blk = s.beg / BLOCK_BITS;
    info->end_blk = (s.end - 1) / BLOCK_BITS;
    for (uint32_t blk = info->beg_blk; blk <= info->end_blk; ++blk) {
        uint64_t w[BLOCK_WORDS];
        get_block(w, s, blk);
        if (block_empty(w)) continue;
        info->block_count += 1;
        if (!block_full(w)) info->block_data_bytes += sizeof(w);
    }
    const size_t raw = sizeof(uint16_t) + info->block_count * sizeof(uint16_t) + info->block_data_bytes;
    info->input_size = raw > SERIAL_MIN_INPUT_SIZE ? raw : SERIAL_MIN_INPUT_SIZE;
}

// V1: writes the blocks of s, uncompressed, to raw (info->input_size bytes, zeroed)
static void write_blocks(uint8_t* raw, const serial_info_t* info, md_bitset_t s) {
    const uint16_t count = (uint16_t)info->block_count;
    MEMCPY(raw, &count, sizeof(count));
    uint8_t* block_indices = raw + sizeof(uint16_t);
    uint8_t* block_data = block_indices + info->block_count * sizeof(uint16_t);

    size_t index_pos = 0, data_pos = 0;
    for (uint32_t blk = info->beg_blk; blk <= info->end_blk && info->block_count; ++blk) {
        uint64_t w[BLOCK_WORDS];
        get_block(w, s, blk);
        if (block_empty(w)) continue;
        uint16_t idx = (uint16_t)blk;
        if (block_full(w)) {
            idx |= V1_FLAG_ALL_SET;
        } else {
            MEMCPY(block_data + data_pos, w, sizeof(w));
            data_pos += sizeof(w);
        }
        MEMCPY(block_indices + index_pos * sizeof(uint16_t), &idx, sizeof(idx));
        index_pos += 1;
    }
    ASSERT(index_pos == info->block_count);
    ASSERT(data_pos == info->block_data_bytes);
}

static inline uint16_t get_u16(const uint8_t* src) { uint16_t v; MEMCPY(&v, src, sizeof(v)); return v; }

// V1: reads the blocks written by write_blocks
static bool read_blocks(md_bitset_t* dst, const uint8_t* data, size_t size, md_allocator_i* alloc) {
    if (size < sizeof(uint16_t)) return false;
    const size_t block_count = get_u16(data);
    if (block_count == 0) {
        *dst = md_bitset_empty_set();
        return true;
    }
    const size_t header_bytes = sizeof(uint16_t) + block_count * sizeof(uint16_t);
    if (size < header_bytes) return false;

    const uint8_t* idx_ptr = data + sizeof(uint16_t);
    size_t data_bytes = 0;
    for (size_t i = 0; i < block_count; ++i) {
        const uint16_t v = get_u16(idx_ptr + i * sizeof(uint16_t));
        if (i > 0 && (v & V1_IDX_MASK) <= (get_u16(idx_ptr + (i - 1) * sizeof(uint16_t)) & V1_IDX_MASK)) return false;
        if (!(v & V1_FLAG_ALL_SET)) data_bytes += BLOCK_WORDS * sizeof(uint64_t);
    }
    if (size < header_bytes + data_bytes) return false;

    const uint32_t beg_blk = get_u16(idx_ptr) & V1_IDX_MASK;
    const uint32_t end_blk = get_u16(idx_ptr + (block_count - 1) * sizeof(uint16_t)) & V1_IDX_MASK;

    // Through a mask: it holds the words of the blocks there are, not those of the gaps between them
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_bitmask_t m = md_bitmask_create(beg_blk * BLOCK_BITS, (end_blk + 1) * BLOCK_BITS, md_temp_allocator(temp));
    const uint8_t* blk_data = data + header_bytes;
    for (size_t i = 0; i < block_count; ++i) {
        const uint16_t v = get_u16(idx_ptr + i * sizeof(uint16_t));
        const uint32_t base = (v & V1_IDX_MASK) * BLOCK_BITS;
        if (v & V1_FLAG_ALL_SET) {
            md_bitmask_set_range(&m, base, base + BLOCK_BITS);
        } else {
            uint64_t w[BLOCK_WORDS];
            MEMCPY(w, blk_data, sizeof(w));
            blk_data += sizeof(w);
            // Not a tidy set (its span is the block's), which the operations of a mask do not need
            const md_bitset_t block = {base, base + BLOCK_BITS, w};
            md_bitmask_or(&m, block);
        }
    }
    *dst = md_bitmask_extract(&m, alloc);
    md_temp_end(temp);
    return !md_bitset_empty(*dst);
}

// V2: the entries and the data of the blocks of s, written to raw when it is given (input_size bytes, zeroed).
// Returns the size of the uncompressed stream.
static size_t v2_blocks(uint8_t* raw, md_bitset_t s) {
    uint32_t num_entries = 0;
    size_t entry_bytes = 0, data_bytes = 0;
    // Two passes when writing: the first sizes the entries, the data follows them
    for (int pass = raw ? 0 : 1; pass < 2; ++pass) {
        uint8_t* entries = raw ? raw + sizeof(uint32_t) : NULL;
        uint8_t* data = (raw && pass == 1) ? raw + sizeof(uint32_t) + entry_bytes : NULL;
        num_entries = 0;
        size_t e_pos = 0, d_pos = 0;
        size_t run_pos = SIZE_MAX;      // where the count of the open run of full blocks is
        uint32_t run_next = 0;          // the block which would extend it
        if (!md_bitset_empty(s)) {
            const uint32_t beg_blk = s.beg / BLOCK_BITS;
            const uint32_t end_blk = (s.end - 1) / BLOCK_BITS;
            for (uint32_t blk = beg_blk; blk <= end_blk; ++blk) {
                uint64_t w[BLOCK_WORDS];
                get_block(w, s, blk);
                if (block_empty(w)) { run_pos = SIZE_MAX; continue; }
                if (block_full(w)) {
                    if (run_pos != SIZE_MAX && run_next == blk) {
                        if (pass == 1 && raw) { uint32_t c; MEMCPY(&c, entries + run_pos, sizeof(c)); c += 1; MEMCPY(entries + run_pos, &c, sizeof(c)); }
                        run_next = blk + 1;
                        continue;
                    }
                    if (pass == 1 && raw) {
                        const uint32_t v = blk | V2_FLAG_ALL_SET, c = 1;
                        MEMCPY(entries + e_pos, &v, sizeof(v));
                        MEMCPY(entries + e_pos + sizeof(v), &c, sizeof(c));
                    }
                    run_pos = e_pos + sizeof(uint32_t);
                    run_next = blk + 1;
                    e_pos += 2 * sizeof(uint32_t);
                } else {
                    if (pass == 1 && raw) {
                        MEMCPY(entries + e_pos, &blk, sizeof(blk));
                        MEMCPY(data + d_pos, w, sizeof(w));
                    }
                    run_pos = SIZE_MAX;
                    e_pos += sizeof(uint32_t);
                    d_pos += sizeof(w);
                }
                num_entries += 1;
            }
        }
        entry_bytes = e_pos;
        data_bytes = d_pos;
    }
    if (raw) MEMCPY(raw, &num_entries, sizeof(num_entries));
    const size_t size = sizeof(uint32_t) + entry_bytes + data_bytes;
    return num_entries ? (size > SERIAL_MIN_INPUT_SIZE ? size : SERIAL_MIN_INPUT_SIZE) : 0;
}

static bool v2_read(md_bitset_t* dst, const uint8_t* data, size_t size, md_allocator_i* alloc) {
    if (size < sizeof(uint32_t)) return false;
    uint32_t num_entries;
    MEMCPY(&num_entries, data, sizeof(num_entries));
    if (num_entries == 0) return false;

    // The entries first: their extent, and where the data starts
    size_t pos = sizeof(uint32_t), data_bytes = 0;
    uint64_t next = 0;      // the first block an entry may name
    uint32_t beg_blk = 0, end_blk = 0;
    for (uint32_t i = 0; i < num_entries; ++i) {
        uint32_t v;
        if (pos + sizeof(v) > size) return false;
        MEMCPY(&v, data + pos, sizeof(v));
        pos += sizeof(v);
        const uint32_t blk = v & V2_IDX_MASK;
        uint64_t count = 1;
        if (v & V2_FLAG_ALL_SET) {
            uint32_t c;
            if (pos + sizeof(c) > size) return false;
            MEMCPY(&c, data + pos, sizeof(c));
            pos += sizeof(c);
            if (c == 0) return false;
            count = c;
        } else {
            data_bytes += BLOCK_WORDS * sizeof(uint64_t);
        }
        if (blk < next || ((uint64_t)blk + count - 1) * BLOCK_BITS >= MAX_INDEX) return false;
        if (i == 0) beg_blk = blk;
        end_blk = (uint32_t)(blk + count - 1);
        next = (uint64_t)blk + count;
    }
    if (pos + data_bytes > size) return false;

    // Through a mask: it holds the words of the blocks there are, not those of the gaps between them. The last
    // block may reach past the indices a set can hold: the domain drops what lies beyond.
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    const uint32_t dom_end = (uint32_t)MIN((uint64_t)(end_blk + 1) * BLOCK_BITS, (uint64_t)MAX_INDEX);
    md_bitmask_t m = md_bitmask_create(beg_blk * BLOCK_BITS, dom_end, md_temp_allocator(temp));
    const uint8_t* blk_data = data + pos;
    pos = sizeof(uint32_t);
    for (uint32_t i = 0; i < num_entries; ++i) {
        uint32_t v;
        MEMCPY(&v, data + pos, sizeof(v));
        pos += sizeof(v);
        const uint32_t base = (v & V2_IDX_MASK) * BLOCK_BITS;
        if (v & V2_FLAG_ALL_SET) {
            uint32_t c;
            MEMCPY(&c, data + pos, sizeof(c));
            pos += sizeof(c);
            md_bitmask_set_range(&m, base, (uint32_t)MIN((uint64_t)base + (uint64_t)c * BLOCK_BITS, (uint64_t)MAX_INDEX));
        } else {
            uint64_t w[BLOCK_WORDS];
            MEMCPY(w, blk_data, sizeof(w));
            blk_data += sizeof(w);
            // Not a tidy set (its span is the block's), which the operations of a mask do not need
            const md_bitset_t block = {base, (uint32_t)MIN((uint64_t)base + BLOCK_BITS, (uint64_t)MAX_INDEX), w};
            md_bitmask_or(&m, block);
        }
    }
    *dst = md_bitmask_extract(&m, alloc);
    md_temp_end(temp);
    return !md_bitset_empty(*dst);
}

size_t md_bitset_serialize_size_in_bytes(md_bitset_t s) {
    if (md_bitset_empty(s)) return 0;
    if (fits_v1(s)) {
        serial_info_t info;
        serial_info(&info, s);
        return info.block_count ? fastlz_bound(info.input_size) : 0;
    }
    const size_t raw_size = v2_blocks(NULL, s);
    return raw_size ? V2_HEADER_SIZE + fastlz_bound(raw_size) : 0;
}

size_t md_bitset_serialize(void* dst, md_bitset_t s) {
    if (!dst || md_bitset_empty(s)) return 0;
    const bool v1 = fits_v1(s);
    serial_info_t info = {0};
    size_t raw_size = 0;
    if (v1) {
        serial_info(&info, s);
        raw_size = info.block_count ? info.input_size : 0;
    } else {
        raw_size = v2_blocks(NULL, s);
    }
    if (raw_size == 0) return 0;
    if (raw_size > INT32_MAX) {
        MD_LOG_ERROR("Bitset is too large to serialize");
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    size_t result = 0;
    uint8_t* raw = (uint8_t*)md_temp_alloc_zero(temp, raw_size);
    if (raw) {
        uint8_t* out = (uint8_t*)dst;
        if (v1) {
            write_blocks(raw, &info, s);
        } else {
            v2_blocks(raw, s);
            const uint32_t raw_size32 = (uint32_t)raw_size;
            MEMCPY(out, v2_magic, V2_MAGIC_SIZE);
            MEMCPY(out + V2_MAGIC_SIZE, &raw_size32, sizeof(raw_size32));
            out += V2_HEADER_SIZE;
        }
        const int lz_bytes = fastlz_compress_level(2, raw, (int)raw_size, out);
        result = lz_bytes > 0 ? (v1 ? 0 : V2_HEADER_SIZE) + (size_t)lz_bytes : 0;
    }
    md_temp_end(temp);
    return result;
}

static bool deserialize_v2(md_bitset_t* dst, const uint8_t* src, size_t num_bytes, md_allocator_i* alloc) {
    if (num_bytes <= V2_HEADER_SIZE) return false;
    uint32_t raw_size;
    MEMCPY(&raw_size, src + V2_MAGIC_SIZE, sizeof(raw_size));
    if (raw_size < SERIAL_MIN_INPUT_SIZE || raw_size > INT32_MAX) return false;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    bool result = false;
    uint8_t* raw = (uint8_t*)md_temp_alloc(temp, raw_size);
    if (raw && fastlz_decompress(src + V2_HEADER_SIZE, (int)(num_bytes - V2_HEADER_SIZE), raw, (int)raw_size) == (int)raw_size) {
        result = v2_read(dst, raw, raw_size, alloc);
    }
    md_temp_end(temp);
    return result;
}

static bool deserialize(md_bitset_t* dst, const void* src, size_t num_bytes, md_allocator_i* alloc) {
    if (num_bytes >= V2_MAGIC_SIZE && memcmp(src, v2_magic, V2_MAGIC_SIZE) == 0) {
        return deserialize_v2(dst, (const uint8_t*)src, num_bytes, alloc);
    }

    // V1 does not store its raw size: decompress into a buffer which grows until it fits
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;
    size_t cap = num_bytes * 4;
    if (cap < SERIAL_MIN_INPUT_SIZE) cap = SERIAL_MIN_INPUT_SIZE;
    if (cap > V1_MAX_RAW_SIZE)       cap = V1_MAX_RAW_SIZE;

    while (cap <= V1_MAX_RAW_SIZE) {
        md_vm_arena_set_pos_back(temp_alloc, temp.pos);
        uint8_t* mem = (uint8_t*)md_temp_alloc(temp, cap);
        if (!mem) break;
        int size = fastlz_decompress(src, (int)num_bytes, mem, (int)cap);
        if (size > 0) {
            result = read_blocks(dst, mem, (size_t)size, alloc);
            break;
        }
        if (cap == V1_MAX_RAW_SIZE) break;
        cap = cap * 2 < V1_MAX_RAW_SIZE ? cap * 2 : V1_MAX_RAW_SIZE;
    }
    md_temp_end(temp);
    return result;
}

bool md_bitset_deserialize(md_bitset_t* dst, const void* src, size_t num_bytes, md_allocator_i* alloc) {
    ASSERT(dst);
    if (!dst || !src || num_bytes == 0 || num_bytes > INT32_MAX) return false;

    // FastLZ may read a byte past the end of a malformed stream: decompress from a copy which has zeroed padding
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    bool result = false;
    uint8_t* in = (uint8_t*)md_temp_alloc_zero(temp, num_bytes + 8);
    if (in) {
        MEMCPY(in, src, num_bytes);
        result = deserialize(dst, in, num_bytes, alloc);
    }
    md_temp_end(temp);

    if (!result) {
        MD_LOG_ERROR("Failed to deserialize bitset");
    }
    return result;
}

// ---------------------------------------------------------------------------------------------------------
// md_bitmask_t
//
// The words of the allocation are [g0, g0 + cap) (word indices). Of those, [lo, hi) are written; the others are
// undefined and never read. Writing a word outside [lo, hi) first zeroes the words between it and [lo, hi), so
// the written words stay contiguous and clearing is just lo = hi.

#define W(m, g) ((m)->words[(g) - (m)->g0])

static inline uint32_t dom_w0(const md_bitmask_t* m) { return m->beg >> 6; }
static inline uint32_t dom_w1(const md_bitmask_t* m) { return (uint32_t)(((uint64_t)m->end + 63) >> 6); }

void md_bitmask_init(md_bitmask_t* m, uint32_t beg, uint32_t end, md_allocator_i* alloc) {
    ASSERT(m);
    ASSERT(alloc);
    ASSERT(end <= MAX_INDEX);
    MEMSET(m, 0, sizeof(*m));
    m->alloc = alloc;
    if (beg < end) {
        m->beg = beg;
        m->end = end;
    }
}

md_bitmask_t md_bitmask_create(uint32_t beg, uint32_t end, md_allocator_i* alloc) {
    md_bitmask_t m;
    md_bitmask_init(&m, beg, end, alloc);
    return m;
}

void md_bitmask_free(md_bitmask_t* m) {
    ASSERT(m);
    if (m->words) {
        ASSERT(m->alloc);
        md_free(m->alloc, m->words, (size_t)m->cap * sizeof(uint64_t));
    }
    MEMSET(m, 0, sizeof(*m));
}

// Make the allocation cover the words [w0, w1) (within the domain), keeping the written words
static void mask_alloc(md_bitmask_t* m, uint32_t w0, uint32_t w1) {
    ASSERT(w0 < w1);
    if (m->cap && m->g0 <= w0 && w1 <= m->g0 + m->cap) return;
    ASSERT(m->alloc);

    const uint32_t d0 = dom_w0(m), d1 = dom_w1(m);
    ASSERT(d0 <= w0 && w1 <= d1);

    uint32_t lo = w0, hi = w1;
    if (m->lo < m->hi) {
        lo = min_u32(lo, m->lo);
        hi = max_u32(hi, m->hi);
    }
    // Grow geometrically, the slack on the side that grew, within the domain
    uint32_t cap = max_u32(max_u32(hi - lo, m->cap * 2), 8);
    cap = min_u32(cap, d1 - d0);
    uint32_t g0 = (m->cap && w0 < m->g0) ? (hi > d0 + cap ? hi - cap : d0) : lo;
    if (g0 + cap > d1) g0 = d1 - cap;
    ASSERT(d0 <= g0 && g0 <= lo && hi <= g0 + cap);

    uint64_t* words = (uint64_t*)md_alloc(m->alloc, (size_t)cap * sizeof(uint64_t));
    ASSERT(words);
    if (m->lo < m->hi) {
        MEMCPY(words + (m->lo - g0), m->words + (m->lo - m->g0), (size_t)(m->hi - m->lo) * sizeof(uint64_t));
    }
    if (m->words) {
        md_free(m->alloc, m->words, (size_t)m->cap * sizeof(uint64_t));
    }
    m->words = words;
    m->g0 = g0;
    m->cap = cap;
}

// Make the words [w0, w1) written (zero where they were not), w0 < w1, within the domain
static void mask_touch(md_bitmask_t* m, uint32_t w0, uint32_t w1) {
    if (m->lo < m->hi && m->lo <= w0 && w1 <= m->hi) return;
    mask_alloc(m, w0, w1);
    if (m->lo >= m->hi) {
        MEMSET(&W(m, w0), 0, (size_t)(w1 - w0) * sizeof(uint64_t));
        m->lo = w0;
        m->hi = w1;
        return;
    }
    if (w0 < m->lo) {
        MEMSET(&W(m, w0), 0, (size_t)(m->lo - w0) * sizeof(uint64_t));
        m->lo = w0;
    }
    if (w1 > m->hi) {
        MEMSET(&W(m, m->hi), 0, (size_t)(w1 - m->hi) * sizeof(uint64_t));
        m->hi = w1;
    }
}

// Clamp [*beg, *end) to the domain; false if nothing is left
static inline bool clamp_to_domain(const md_bitmask_t* m, uint32_t* beg, uint32_t* end) {
    *beg = max_u32(*beg, m->beg);
    *end = min_u32(*end, m->end);
    return *beg < *end;
}

void md_bitmask_reserve(md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (!clamp_to_domain(m, &beg, &end)) return;
    mask_alloc(m, beg >> 6, ((end - 1) >> 6) + 1);
}

void md_bitmask_set_domain(md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    ASSERT(end <= MAX_INDEX);
    if (beg >= end) beg = end = 0;
    m->beg = beg;
    m->end = end;
    if (m->lo >= m->hi) return;
    if (beg >= end) {
        m->lo = m->hi = 0;
        return;
    }
    // Keep the written words within the new domain, without the bits outside it
    const uint32_t lo = max_u32(m->lo, dom_w0(m));
    const uint32_t hi = min_u32(m->hi, dom_w1(m));
    if (lo >= hi) {
        m->lo = m->hi = 0;
        return;
    }
    m->lo = lo;
    m->hi = hi;
    W(m, lo)     &= word_range_mask(lo,     beg, end);
    W(m, hi - 1) &= word_range_mask(hi - 1, beg, end);
    // An allocation reaching outside the new domain is fine: nothing is written there
}

void md_bitmask_clear(md_bitmask_t* m) {
    ASSERT(m);
    m->lo = m->hi = 0;
}

void md_bitmask_set_bit(md_bitmask_t* m, uint32_t idx) {
    ASSERT(m);
    if (idx - m->beg >= m->end - m->beg) return;  // outside the domain
    const uint32_t g = idx >> 6;
    if (g - m->lo >= m->hi - m->lo) mask_touch(m, g, g + 1);
    W(m, g) |= 1ull << (idx & 63);
}

void md_bitmask_clear_bit(md_bitmask_t* m, uint32_t idx) {
    ASSERT(m);
    const uint32_t g = idx >> 6;
    if (g - m->lo < m->hi - m->lo) W(m, g) &= ~(1ull << (idx & 63));
}

void md_bitmask_set_range(md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (!clamp_to_domain(m, &beg, &end)) return;
    const uint32_t g0 = beg >> 6, g1 = ((end - 1) >> 6) + 1;
    mask_touch(m, g0, g1);
    for (uint32_t g = g0; g < g1; ++g) W(m, g) |= word_range_mask(g, beg, end);
}

void md_bitmask_clear_range(md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (m->lo >= m->hi) return;
    beg = max_u32(beg, m->lo << 6);
    end = (uint32_t)MIN((uint64_t)end, (uint64_t)m->hi << 6);
    if (beg >= end) return;
    const uint32_t g0 = beg >> 6, g1 = ((end - 1) >> 6) + 1;
    for (uint32_t g = g0; g < g1; ++g) W(m, g) &= ~word_range_mask(g, beg, end);
}

void md_bitmask_invert_range(md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (!clamp_to_domain(m, &beg, &end)) return;
    const uint32_t g0 = beg >> 6, g1 = ((end - 1) >> 6) + 1;
    mask_touch(m, g0, g1);
    for (uint32_t g = g0; g < g1; ++g) W(m, g) ^= word_range_mask(g, beg, end);
}

void md_bitmask_set_indices(md_bitmask_t* m, const uint32_t* indices, size_t num_indices) {
    ASSERT(m);
    if (!indices || !num_indices) return;
    uint32_t lo = UINT32_MAX, hi = 0;
    for (size_t i = 0; i < num_indices; ++i) {
        const uint32_t idx = indices[i];
        if (idx - m->beg >= m->end - m->beg) continue;
        lo = min_u32(lo, idx);
        hi = max_u32(hi, idx);
    }
    if (lo > hi) return;
    mask_touch(m, lo >> 6, (hi >> 6) + 1);
    for (size_t i = 0; i < num_indices; ++i) {
        const uint32_t idx = indices[i];
        if (idx - m->beg >= m->end - m->beg) continue;
        W(m, idx >> 6) |= 1ull << (idx & 63);
    }
}

// The words of s within the domain: [*w0, *w1). False if s has no members within the domain.
static inline bool operand_words(const md_bitmask_t* m, md_bitset_t s, uint32_t* w0, uint32_t* w1) {
    if (s.beg >= s.end) return false;
    uint32_t beg = s.beg, end = s.end;
    if (!clamp_to_domain(m, &beg, &end)) return false;
    *w0 = beg >> 6;
    *w1 = ((end - 1) >> 6) + 1;
    return true;
}

static inline uint64_t domain_word_mask(const md_bitmask_t* m, uint32_t g) {
    return word_range_mask(g, m->beg, m->end);
}

void md_bitmask_assign(md_bitmask_t* m, md_bitset_t s) {
    ASSERT(m);
    uint32_t w0, w1;
    if (!operand_words(m, s, &w0, &w1)) {
        md_bitmask_clear(m);
        return;
    }
    // A view of m lies within [lo, hi), so touching moves nothing; each word is read before it is written
    mask_touch(m, w0, w1);
    for (uint32_t g = w0; g < w1; ++g) W(m, g) = md_bitset_word(s, g << 6) & domain_word_mask(m, g);
    m->lo = w0;
    m->hi = w1;
}

void md_bitmask_or(md_bitmask_t* m, md_bitset_t s) {
    ASSERT(m);
    uint32_t w0, w1;
    if (!operand_words(m, s, &w0, &w1)) return;
    mask_touch(m, w0, w1);
    for (uint32_t g = w0; g < w1; ++g) W(m, g) |= md_bitset_word(s, g << 6) & domain_word_mask(m, g);
}

void md_bitmask_xor(md_bitmask_t* m, md_bitset_t s) {
    ASSERT(m);
    uint32_t w0, w1;
    if (!operand_words(m, s, &w0, &w1)) return;
    mask_touch(m, w0, w1);
    for (uint32_t g = w0; g < w1; ++g) W(m, g) ^= md_bitset_word(s, g << 6) & domain_word_mask(m, g);
}

void md_bitmask_and(md_bitmask_t* m, md_bitset_t s) {
    ASSERT(m);
    for (uint32_t g = m->lo; g < m->hi; ++g) W(m, g) &= md_bitset_word(s, g << 6);
}

void md_bitmask_andnot(md_bitmask_t* m, md_bitset_t s) {
    ASSERT(m);
    if (s.beg >= s.end || m->lo >= m->hi) return;
    const uint32_t w0 = max_u32(m->lo, s.beg >> 6);
    const uint32_t w1 = min_u32(m->hi, (uint32_t)(((uint64_t)s.end + 63) >> 6));
    for (uint32_t g = w0; g < w1; ++g) W(m, g) &= ~md_bitset_word(s, g << 6);
}

bool md_bitmask_empty(const md_bitmask_t* m) {
    ASSERT(m);
    for (uint32_t g = m->lo; g < m->hi; ++g) if (W(m, g)) return false;
    return true;
}

uint64_t md_bitmask_count(const md_bitmask_t* m) {
    ASSERT(m);
    uint64_t c = 0;
    for (uint32_t g = m->lo; g < m->hi; ++g) c += popcnt64(W(m, g));
    return c;
}

// Clamp [*beg, *end) to the written words; false if nothing is left
static inline bool clamp_to_written(const md_bitmask_t* m, uint32_t* beg, uint32_t* end) {
    if (m->lo >= m->hi) return false;
    *beg = max_u32(*beg, m->lo << 6);
    *end = (uint32_t)MIN((uint64_t)*end, (uint64_t)m->hi << 6);
    return *beg < *end;
}

uint64_t md_bitmask_count_range(const md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (!clamp_to_written(m, &beg, &end)) return 0;
    const uint32_t g0 = beg >> 6, g1 = ((end - 1) >> 6) + 1;
    uint64_t count = 0;
    for (uint32_t g = g0; g < g1; ++g) count += popcnt64(W(m, g) & word_range_mask(g, beg, end));
    return count;
}

bool md_bitmask_test_any_range(const md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (!clamp_to_written(m, &beg, &end)) return false;
    const uint32_t g0 = beg >> 6, g1 = ((end - 1) >> 6) + 1;
    for (uint32_t g = g0; g < g1; ++g) {
        if (W(m, g) & word_range_mask(g, beg, end)) return true;
    }
    return false;
}

bool md_bitmask_test_all_range(const md_bitmask_t* m, uint32_t beg, uint32_t end) {
    ASSERT(m);
    if (beg >= end) return true;
    return md_bitmask_count_range(m, beg, end) == (uint64_t)(end - beg);
}

md_bitset_t md_bitmask_view(const md_bitmask_t* m) {
    ASSERT(m);
    uint32_t gf = m->lo, gl = m->hi;
    while (gf < gl && W(m, gf) == 0) ++gf;
    if (gf == gl) return md_bitset_empty_set();
    while (W(m, gl - 1) == 0) --gl;
    md_bitset_t s;
    s.beg = (gf << 6) + (uint32_t)ctz64(W(m, gf));
    s.end = ((gl - 1) << 6) + 64 - (uint32_t)clz64(W(m, gl - 1));
    s.words = &W(m, gf);
    return s;
}

md_bitset_t md_bitmask_extract(const md_bitmask_t* m, md_allocator_i* alloc) {
    return md_bitset_copy(md_bitmask_view(m), alloc);
}

#undef W
