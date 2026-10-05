#include <core/md_bitset.h>

#include <core/md_allocator.h>
#include <core/md_bitfield.h>
#include <core/md_common.h>
#include <core/md_hash.h>

#include <string.h>

STATIC_ASSERT(sizeof(md_bitset_gpu_t) == 16, "GPU descriptor must be 16 bytes");

#define MAX_INDEX 0xFFFFFFC0u  // sets hold indices < 2^32 - 64, so a window base + 64 never wraps

// ---------------------------------------------------------------------------------------------------------
// Payload

static inline size_t payload_bytes(size_t num_words) {
    return (1 + num_words) * sizeof(uint64_t);
}

static inline uint64_t* alloc_payload(md_allocator_i* alloc, size_t num_words, uint64_t count) {
    ASSERT(alloc);
    uint64_t* p = (uint64_t*)md_alloc(alloc, payload_bytes(num_words));
    ASSERT(p);
    ASSERT(((uintptr_t)p & 7) == 0);
    p[0] = count;
    return p + 1;
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
// is either a constant word (0 or ~0) or a contiguous slice of its payload; RUN edge words, which are
// partial, are isolated into single word intervals.
static inline int push_breakpoints(uint64_t* bp, int n, md_bitset_t s) {
    if (s.beg == s.end) return n;
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
    if (s.words) {
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
    uint64_t* words = alloc_payload(alloc, nw, st.count);
    stats_t st2 = {0};
    eval_op(&st2, words, st.first_g, a, b, op, st.first_g, st.last_g + 1);
    ASSERT(st2.count == st.count);
    r.words = words;
    return r;
}

static inline uint32_t min_u32(uint32_t a, uint32_t b) { return a < b ? a : b; }
static inline uint32_t max_u32(uint32_t a, uint32_t b) { return a > b ? a : b; }

// ---------------------------------------------------------------------------------------------------------
// Producers

md_bitset_t md_bitset_copy(md_bitset_t s, md_allocator_i* alloc) {
    if (!s.words) return s;
    const size_t nw = md_bitset_num_words(s);
    uint64_t* words = alloc_payload(alloc, nw, s.words[-1]);
    MEMCPY(words, s.words, nw * sizeof(uint64_t));
    s.words = words;
    return s;
}

void md_bitset_free(md_bitset_t* s, md_allocator_i* alloc) {
    ASSERT(s);
    if (s->words) {
        ASSERT(alloc);
        md_free(alloc, (void*)(s->words - 1), payload_bytes(md_bitset_num_words(*s)));
    }
    *s = md_bitset_empty_set();
}

md_bitset_t md_bitset_from_words(const uint64_t* words, size_t num_words, uint32_t base, md_allocator_i* alloc) {
    ASSERT((base & 63) == 0);
    ASSERT(((uint64_t)base + (uint64_t)num_words * 64) <= (uint64_t)MAX_INDEX + 64);
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
    uint64_t* dst = alloc_payload(alloc, nw, count);
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
            // b punches a hole in a: needs a payload
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

// ---------------------------------------------------------------------------------------------------------
// Readers

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
    // beg/end of a BITS set are members, so a range covering either is a hit without reading words
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
    if (a.beg != b.beg || a.end != b.end) return false;
    if (!a.words || !b.words) return a.words == b.words || a.beg == a.end;
    if (a.words == b.words) return true;
    return memcmp(a.words - 1, b.words - 1, payload_bytes(md_bitset_num_words(a))) == 0;
}

uint64_t md_bitset_hash64(md_bitset_t s, uint64_t seed) {
    const uint32_t range[2] = {s.beg, s.end};
    uint64_t h = md_hash64(range, sizeof(range), seed);
    if (s.words) h = md_hash64(s.words, md_bitset_num_words(s) * sizeof(uint64_t), h);
    return h;
}

size_t md_bitset_extract_indices(uint32_t* out, size_t cap, md_bitset_t s) {
    size_t n = 0;
    md_bitset_iter_t it = md_bitset_iter(s);
    while (n < cap && md_bitset_iter_next(&it)) out[n++] = it.idx;
    return n;
}

bool md_bitset_validate(md_bitset_t s) {
    if (s.beg == s.end) return s.beg == 0 && s.words == NULL;
    if (s.beg > s.end || s.end > MAX_INDEX) return false;
    if (!s.words) return true;

    const size_t nw = md_bitset_num_words(s);
    uint64_t count = 0;
    for (size_t k = 0; k < nw; ++k) count += popcnt64(s.words[k]);
    if (count != s.words[-1]) return false;
    if (count == (uint64_t)(s.end - s.beg)) return false;  // should have been a RUN

    const uint32_t lo = s.beg & 63;
    const uint32_t hi = (s.end - 1) & 63;
    if (((s.words[0] >> lo) & 1) == 0) return false;                       // beg is a member
    if (lo && (s.words[0] & ~(~0ull << lo))) return false;                 // nothing below beg
    if (((s.words[nw - 1] >> hi) & 1) == 0) return false;                  // end - 1 is a member
    if (hi != 63 && (s.words[nw - 1] & (~0ull << (hi + 1)))) return false; // nothing at or above end
    return true;
}

// ---------------------------------------------------------------------------------------------------------
// md_bitfield interop

md_bitset_t md_bitset_from_bitfield(const md_bitfield_t* bf, md_allocator_i* alloc) {
    ASSERT(bf);
    uint64_t first, last;
    if (!md_bitfield_get_range(&first, &last, bf)) return md_bitset_empty_set();
    ASSERT(last < MAX_INDEX);

    const uint32_t base = (uint32_t)first & ~63u;
    const size_t nw = (size_t)(last >> 6) - (size_t)(first >> 6) + 1;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    uint64_t* scratch = (uint64_t*)md_temp_alloc_zero(temp, nw * sizeof(uint64_t));
    md_bitfield_iter_t it = md_bitfield_iter_range_create(bf, first, last + 1);
    while (md_bitfield_iter_next(&it)) {
        const uint32_t k = (uint32_t)md_bitfield_iter_idx(&it) - base;
        scratch[k >> 6] |= 1ull << (k & 63);
    }
    md_bitset_t r = md_bitset_from_words(scratch, nw, base, alloc);
    md_temp_end(temp);
    return r;
}

void md_bitset_to_bitfield(md_bitfield_t* dst, md_bitset_t s) {
    ASSERT(dst);
    md_bitfield_clear(dst);
    if (md_bitset_empty(s)) return;
    if (!s.words) {
        md_bitfield_set_range(dst, s.beg, s.end);
        return;
    }
    md_bitfield_reserve_range(dst, s.beg, s.end);
    md_bitset_iter_t it = md_bitset_iter(s);
    while (md_bitset_iter_next(&it)) md_bitfield_set_bit(dst, it.idx);
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
        const md_bitset_t s = sets[i];
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
