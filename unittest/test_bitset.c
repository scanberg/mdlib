#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_tracking_allocator.h>
#include <core/md_bitset.h>
#include <core/md_base64.h>

#include <stdlib.h>
#include <string.h>

// Differential tests against a plain bool-per-index reference. The reference universe is small, but the
// sets are placed at varying offsets (including large ones and ones straddling word boundaries) so every
// edge case of the span arithmetic is hit.


// Counting allocator: tests assert that everything allocated is freed, with the size it was allocated with
typedef struct { md_allocator_i iface; int64_t live_bytes; int64_t live_allocs; } count_alloc_t;
static void* count_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file; (void)line;
    count_alloc_t* c = (count_alloc_t*)inst;
    if (ptr) { c->live_bytes -= (int64_t)old_size; c->live_allocs -= 1; }
    if (new_size == 0) { free(ptr); return NULL; }
    c->live_bytes += (int64_t)new_size; c->live_allocs += 1;
    return realloc(ptr, new_size);
}
static md_allocator_i* count_alloc_init(count_alloc_t* c) {
    c->live_bytes = 0; c->live_allocs = 0;
    c->iface.inst = (struct md_allocator_o*)c;
    c->iface.realloc = count_realloc;
    return &c->iface;
}

#define REF_N 1536

typedef struct {
    uint32_t base;       // index of ref[0]
    uint8_t  ref[REF_N];
} ref_t;

static uint64_t rng_state = 0x9E3779B97F4A7C15ull;
static inline uint32_t rnd(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return (uint32_t)(rng_state >> 16);
}
static inline uint32_t rnd_range(uint32_t lo, uint32_t hi) { return lo + rnd() % (hi - lo); }  // [lo, hi)

// Shapes that matter for a single-span representation
enum { SHAPE_EMPTY, SHAPE_SINGLE, SHAPE_RUN, SHAPE_RUN_EDGE, SHAPE_SPARSE, SHAPE_DENSE_HOLES, SHAPE_TWO_CLUSTERS, SHAPE_COUNT };

static void gen(uint8_t* ref, int shape) {
    memset(ref, 0, REF_N);
    switch (shape) {
    case SHAPE_EMPTY: break;
    case SHAPE_SINGLE: ref[rnd_range(0, REF_N)] = 1; break;
    case SHAPE_RUN: {
        uint32_t b = rnd_range(0, REF_N - 1), e = rnd_range(b + 1, REF_N);
        memset(ref + b, 1, e - b);
    } break;
    case SHAPE_RUN_EDGE: {  // runs starting/ending exactly at, or one off, a word boundary
        uint32_t b = 64 * rnd_range(1, 8) + (int)rnd_range(0, 3) - 1;
        uint32_t e = b + 64 * rnd_range(1, 6) + (int)rnd_range(0, 3) - 1;
        memset(ref + b, 1, e - b);
    } break;
    case SHAPE_SPARSE:
        for (int i = 0; i < REF_N; ++i) ref[i] = (rnd() % 37) == 0;
        break;
    case SHAPE_DENSE_HOLES: {
        uint32_t b = rnd_range(0, REF_N / 2), e = rnd_range(b + 2, REF_N);
        memset(ref + b, 1, e - b);
        int holes = 1 + rnd() % 5;
        for (int h = 0; h < holes; ++h) ref[rnd_range(b + 1, e - 1 > b + 1 ? e - 1 : b + 2)] = 0;
    } break;
    case SHAPE_TWO_CLUSTERS: {
        for (int i = 0; i < 40; ++i) ref[rnd_range(0, 100)] = 1;
        for (int i = 0; i < 40; ++i) ref[rnd_range(REF_N - 100, REF_N)] = 1;
    } break;
    }
}

static md_bitset_t make(const ref_t* r, md_allocator_i* alloc) {
    uint32_t idx[REF_N];
    size_t n = 0;
    for (uint32_t i = 0; i < REF_N; ++i) if (r->ref[i]) idx[n++] = r->base + i;
    // Shuffle: from_indices accepts any order
    for (size_t i = n; i > 1; --i) { size_t j = rnd() % i; uint32_t t = idx[i-1]; idx[i-1] = idx[j]; idx[j] = t; }
    return md_bitset_from_indices(idx, n, alloc);
}

static bool ref_get(const ref_t* r, uint64_t i) {
    return i >= r->base && i < (uint64_t)r->base + REF_N && r->ref[i - r->base];
}

// Checks every observable property of s against r
static bool matches_form(md_bitset_t s, const ref_t* r, bool tidy) {
    if (!md_bitset_validate(s, tidy)) return false;

    const uint32_t lo = r->base >= 128 ? r->base - 128 : 0;
    const uint32_t hi = r->base + REF_N + 128;
    uint64_t count = 0;
    for (uint32_t i = lo; i < hi; ++i) {
        const bool e = ref_get(r, i);
        if (md_bitset_test(s, i) != e) return false;
        count += e;
    }
    if (md_bitset_count(s) != count) return false;

    // Window view
    for (uint32_t base = lo & ~63u; base < hi; base += 64) {
        uint64_t m = 0;
        for (uint32_t k = 0; k < 64; ++k) m |= (uint64_t)ref_get(r, base + k) << k;
        if (md_bitset_word(s, base) != m) return false;
    }

    // Iteration visits exactly the members, ascending
    md_bitset_iter_t it = md_bitset_iter(s);
    uint64_t seen = 0;
    int64_t prev = -1;
    while (md_bitset_iter_next(&it)) {
        if ((int64_t)it.idx <= prev || !ref_get(r, it.idx)) return false;
        prev = it.idx;
        ++seen;
    }
    return seen == count;
}

// A set a producer returned: tidy as well
static bool matches(md_bitset_t s, const ref_t* r) {
    return matches_form(s, r, true);
}

static const int NUM_ROUNDS = 400;

static uint32_t pick_base(void) {
    static const uint32_t bases[] = {0, 64, 1, 63, 65, 4096 - 5, 1u << 20, (1u << 24) + 31};
    return bases[rnd() % (sizeof(bases) / sizeof(bases[0]))];
}

UTEST(bitset, construct) {
    count_alloc_t counter; md_allocator_i* alloc = count_alloc_init(&counter);
    static ref_t r;
    for (int round = 0; round < NUM_ROUNDS; ++round) {
        r.base = pick_base();
        gen(r.ref, round % SHAPE_COUNT);
        md_bitset_t s = make(&r, alloc);
        EXPECT_TRUE(matches(s, &r));

        // Same set via dense words must be identical, payload included
        const uint32_t wb = r.base & ~63u;
        const size_t nw = (r.base - wb + REF_N + 63) / 64;
        uint64_t* words = (uint64_t*)calloc(nw, sizeof(uint64_t));
        for (uint32_t i = 0; i < REF_N; ++i) if (r.ref[i]) { uint32_t k = r.base + i - wb; words[k >> 6] |= 1ull << (k & 63); }
        md_bitset_t s2 = md_bitset_from_words(words, nw, wb, alloc);
        free(words);
        EXPECT_TRUE(md_bitset_equal(s, s2));
        EXPECT_EQ(md_bitset_hash64(s, 7), md_bitset_hash64(s2, 7));

        // Deep copy
        md_bitset_t s3 = md_bitset_copy(s, alloc);
        EXPECT_TRUE(md_bitset_equal(s, s3));
        EXPECT_TRUE(s.words == NULL || s3.words != s.words);

        md_bitset_free(&s, alloc);
        md_bitset_free(&s2, alloc);
        md_bitset_free(&s3, alloc);
        EXPECT_TRUE(md_bitset_empty(s));
    }
    EXPECT_EQ(0, counter.live_bytes); EXPECT_EQ(0, counter.live_allocs);
}

UTEST(bitset, canonical_forms) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_bitset_t e = md_bitset_from_indices(NULL, 0, alloc);
    EXPECT_TRUE(md_bitset_empty(e));
    EXPECT_TRUE(md_bitset_validate(e, true));
    EXPECT_FALSE(md_bitset_test(e, 0));

    // A dense set built from indices is a RUN, with no payload
    uint32_t idx[100];
    for (int i = 0; i < 100; ++i) idx[i] = 1000 + i;
    md_bitset_t r = md_bitset_from_indices(idx, 100, alloc);
    EXPECT_TRUE(md_bitset_is_run(r));
    EXPECT_EQ(1000u, r.beg);
    EXPECT_EQ(1100u, r.end);
    EXPECT_TRUE(md_bitset_equal(r, md_bitset_range(1000, 1100)));

    // Punching a hole gives BITS; filling it again gives back the RUN
    md_bitset_t hole = md_bitset_range(1050, 1051);
    md_bitset_t a = md_bitset_andnot(r, hole, alloc);
    EXPECT_FALSE(md_bitset_is_run(a));
    EXPECT_EQ(99u, md_bitset_count(a));
    md_bitset_t b = md_bitset_or(a, hole, alloc);
    EXPECT_TRUE(md_bitset_is_run(b));
    EXPECT_TRUE(md_bitset_equal(b, r));

    // x ^ x is EMPTY in canonical form
    md_bitset_t z = md_bitset_xor(a, a, alloc);
    EXPECT_TRUE(md_bitset_validate(z, true));
    EXPECT_TRUE(md_bitset_empty(z));

    md_bitset_free(&a, alloc);
    md_bitset_free(&b, alloc);
}

typedef md_bitset_t (*binop_fn)(md_bitset_t, md_bitset_t, md_allocator_i*);

UTEST(bitset, binary_ops) {
    count_alloc_t counter; md_allocator_i* alloc = count_alloc_init(&counter);
    static ref_t ra, rb, rr;
    const binop_fn fns[4] = {md_bitset_and, md_bitset_or, md_bitset_andnot, md_bitset_xor};
    const char* names[4] = {"and", "or", "andnot", "xor"};

    for (int round = 0; round < NUM_ROUNDS * 4; ++round) {
        // Same base for both, offset by a random shift, so spans overlap partially, fully or not at all
        const uint32_t base = pick_base() + 2 * REF_N;
        ra.base = base;
        rb.base = base + (int)rnd_range(0, 2 * REF_N) - REF_N;
        gen(ra.ref, rnd() % SHAPE_COUNT);
        gen(rb.ref, rnd() % SHAPE_COUNT);
        md_bitset_t a = make(&ra, alloc);
        md_bitset_t b = make(&rb, alloc);

        for (int op = 0; op < 4; ++op) {
            md_bitset_t c = fns[op](a, b, alloc);
            // Reference over a window covering both
            const uint32_t lo = ra.base < rb.base ? ra.base : rb.base;
            const uint32_t hi = (ra.base > rb.base ? ra.base : rb.base) + REF_N;
            bool ok = md_bitset_validate(c, true);
            uint64_t count = 0;
            for (uint32_t i = lo - 64; ok && i < hi + 64; ++i) {
                const bool x = ref_get(&ra, i), y = ref_get(&rb, i);
                bool e = false;
                switch (op) {
                case 0: e = x && y; break;
                case 1: e = x || y; break;
                case 2: e = x && !y; break;
                case 3: e = x != y; break;
                }
                ok = md_bitset_test(c, i) == e;
                count += e;
            }
            ok = ok && md_bitset_count(c) == count;
            if (!ok) {
                printf("op %s failed (round %d): a=[%u,%u) %s  b=[%u,%u) %s\n", names[op], round,
                    a.beg, a.end, a.words ? "bits" : "run", b.beg, b.end, b.words ? "bits" : "run");
            }
            EXPECT_TRUE(ok);
            md_bitset_free(&c, alloc);
        }

        // not within a domain
        {
            const uint32_t dbeg = ra.base + rnd_range(0, REF_N / 2) - REF_N / 4;
            const uint32_t dend = dbeg + rnd_range(1, REF_N);
            md_bitset_t c = md_bitset_not(a, dbeg, dend, alloc);
            bool ok = md_bitset_validate(c, true);
            for (uint32_t i = dbeg - 64; ok && i < dend + 64; ++i) {
                const bool e = i >= dbeg && i < dend && !ref_get(&ra, i);
                ok = md_bitset_test(c, i) == e;
            }
            EXPECT_TRUE(ok);
            md_bitset_free(&c, alloc);
        }

        md_bitset_free(&a, alloc);
        md_bitset_free(&b, alloc);
    }
    EXPECT_EQ(0, counter.live_bytes); EXPECT_EQ(0, counter.live_allocs);
}

UTEST(bitset, range_queries) {
    md_allocator_i* alloc = md_get_heap_allocator();
    static ref_t r;
    for (int round = 0; round < NUM_ROUNDS; ++round) {
        r.base = pick_base() + 256;
        gen(r.ref, round % SHAPE_COUNT);
        md_bitset_t s = make(&r, alloc);
        for (int q = 0; q < 50; ++q) {
            const uint32_t beg = r.base - 128 + rnd_range(0, REF_N + 256);
            const uint32_t end = beg + rnd_range(0, 300);
            uint64_t count = 0;
            for (uint32_t i = beg; i < end; ++i) count += ref_get(&r, i);
            EXPECT_EQ(count, md_bitset_count_range(s, beg, end));
            EXPECT_EQ(count > 0, md_bitset_test_any_range(s, beg, end));
            EXPECT_EQ(count == (uint64_t)(end - beg), md_bitset_test_all_range(s, beg, end));

            // Iterating a sub range visits exactly the members in it
            md_bitset_iter_t it = md_bitset_iter_range(s, beg, end);
            uint64_t seen = 0;
            bool ok = true;
            while (md_bitset_iter_next(&it)) { ok = ok && it.idx >= beg && it.idx < end && ref_get(&r, it.idx); ++seen; }
            EXPECT_TRUE(ok);
            EXPECT_EQ(count, seen);
        }
        md_bitset_free(&s, alloc);
    }
}

UTEST(bitset, gpu_pack) {
    md_allocator_i* alloc = md_get_heap_allocator();
    enum { NUM = 64 };
    static ref_t refs[NUM];
    md_bitset_t sets[NUM];
    for (int i = 0; i < NUM; ++i) {
        refs[i].base = pick_base();
        gen(refs[i].ref, i % SHAPE_COUNT);
        sets[i] = make(&refs[i], alloc);
    }
    const size_t nw = md_bitset_pack_num_words(sets, NUM);
    md_bitset_gpu_t gsets[NUM];
    uint64_t* words = (uint64_t*)malloc((nw ? nw : 1) * sizeof(uint64_t));
    md_bitset_pack(gsets, words, sets, NUM);

    for (int s = 0; s < NUM; ++s) {
        EXPECT_EQ(md_bitset_count(sets[s]), (uint64_t)gsets[s].count);
        bool ok = true;
        for (uint32_t i = (refs[s].base > 64 ? refs[s].base - 64 : 0); ok && i < refs[s].base + REF_N + 64; ++i) {
            ok = md_bitset_gpu_test(gsets[s], words, i) == ref_get(&refs[s], i);
        }
        EXPECT_TRUE(ok);
        md_bitset_free(&sets[s], alloc);
    }
    free(words);
}

// Producing into a temp arena: the scratch used internally must come from a different arena, otherwise
// rewinding it would release the result.
UTEST(bitset, temp_arena_result) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* talloc = md_temp_allocator(temp);
    static ref_t r;
    r.base = 4096;
    gen(r.ref, SHAPE_SPARSE);
    md_bitset_t s = make(&r, talloc);
    // Allocate after the result so a stale result would be overwritten
    uint8_t* junk = (uint8_t*)md_temp_alloc(temp, 1 << 16);
    memset(junk, 0xAB, 1 << 16);
    EXPECT_TRUE(matches(s, &r));
    md_temp_end(temp);
}

// Exhaustive over span edges: every beg/end in a window around word boundaries, for a RUN and for a BITS set
// with the same span (a hole punched in the middle). The random tests rarely land exactly on a boundary.
UTEST(bitset, span_edges) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const uint32_t off = 1u << 20;
    int failures = 0;
    for (uint32_t b = off; b < off + 192; ++b) {
        for (uint32_t e = b + 1; e <= b + 200; ++e) {
            md_bitset_t sets[2] = {md_bitset_range(b, e), {0}};
            const uint32_t hole = b + (e - b) / 2;
            const bool has_hole = e - b >= 3;
            if (has_hole) sets[1] = md_bitset_andnot(sets[0], md_bitset_range(hole, hole + 1), alloc);

            for (int v = 0; v < (has_hole ? 2 : 1); ++v) {
                const md_bitset_t s = sets[v];
                bool ok = md_bitset_validate(s, true) && (v == 0 ? md_bitset_is_run(s) : !md_bitset_is_run(s));
                ok = ok && s.beg == b && s.end == e;
                uint64_t n = 0;
                for (uint32_t base = (b & ~63u) - 64; ok && base < e + 128; base += 64) {
                    uint64_t m = 0;
                    for (uint32_t k = 0; k < 64; ++k) {
                        const uint32_t i = base + k;
                        const bool member = i >= b && i < e && !(v == 1 && i == hole);
                        m |= (uint64_t)member << k;
                        n += member;
                    }
                    ok = md_bitset_word(s, base) == m;
                }
                md_bitset_iter_t it = md_bitset_iter(s);
                uint64_t seen = 0;
                while (md_bitset_iter_next(&it)) ++seen;
                ok = ok && seen == n && md_bitset_count(s) == n;
                ok = ok && md_bitset_count_range(s, b, e) == n && md_bitset_count_range(s, 0, UINT32_MAX - 64) == n;
                if (!ok) ++failures;
            }
            if (has_hole) md_bitset_free(&sets[1], alloc);
        }
    }
    EXPECT_EQ(0, failures);
}

// A mask against a reference, under a random sequence of every mutating operation. The domain is narrower than
// what the operations reach, so the clamping to it is exercised too.
UTEST(bitset, mask) {
    count_alloc_t counter; md_allocator_i* alloc = count_alloc_init(&counter);
    enum { U = 4096 };
    static uint8_t ref[U];
    static ref_t rs;
    const uint32_t base = 1u << 16;   // keeps every index well away from 0
    for (int round = 0; round < 60; ++round) {
        memset(ref, 0, sizeof(ref));
        const uint32_t dlo = rnd_range(0, 200), dhi = U - rnd_range(0, 200);    // the domain, within [0, U)
        md_bitmask_t m = md_bitmask_create(base + dlo, base + dhi, alloc);
        for (int step = 0; step < 60; ++step) {
            const uint32_t x = rnd_range(0, U), y = rnd_range(x, U + 1 > x + 300 ? x + 300 : U + 1);
            const int op = rnd() % 14;
            rs.base = base + rnd_range(0, U - REF_N);
            gen(rs.ref, rnd() % SHAPE_COUNT);
            md_bitset_t s = (op >= 6 && op <= 10) ? make(&rs, alloc) : md_bitset_empty_set();
            #define IN_DOM(i) ((i) >= dlo && (i) < dhi)
            switch (op) {
            case 0: md_bitmask_set_bit(&m, base + x); if (IN_DOM(x)) ref[x] = 1; break;
            case 1: md_bitmask_clear_bit(&m, base + x); ref[x] = 0; break;
            case 2: md_bitmask_set_range(&m, base + x, base + y); for (uint32_t i = x; i < y; ++i) if (IN_DOM(i)) ref[i] = 1; break;
            case 3: md_bitmask_clear_range(&m, base + x, base + y); for (uint32_t i = x; i < y; ++i) ref[i] = 0; break;
            case 4: md_bitmask_invert_range(&m, base + x, base + y); for (uint32_t i = x; i < y; ++i) if (IN_DOM(i)) ref[i] ^= 1; break;
            case 5: if (step % 17 == 0) { md_bitmask_clear(&m); memset(ref, 0, sizeof(ref)); } break;
            case 6: md_bitmask_or(&m, s);     for (uint32_t i = 0; i < U; ++i) ref[i] |= IN_DOM(i) && ref_get(&rs, base + i); break;
            case 7: md_bitmask_and(&m, s);    for (uint32_t i = 0; i < U; ++i) ref[i] &= ref_get(&rs, base + i); break;
            case 8: md_bitmask_andnot(&m, s); for (uint32_t i = 0; i < U; ++i) ref[i] &= !ref_get(&rs, base + i); break;
            case 9: md_bitmask_xor(&m, s);    for (uint32_t i = 0; i < U; ++i) ref[i] ^= IN_DOM(i) && ref_get(&rs, base + i); break;
            case 10: md_bitmask_assign(&m, s); for (uint32_t i = 0; i < U; ++i) ref[i] = IN_DOM(i) && ref_get(&rs, base + i); break;
            // With a view of the mask itself as the operand
            case 11: md_bitmask_or(&m, md_bitmask_view(&m)); break;
            case 12: md_bitmask_assign(&m, md_bitmask_view(&m)); break;
            case 13: if (step % 13 == 0) { md_bitmask_andnot(&m, md_bitmask_view(&m)); memset(ref, 0, sizeof(ref)); } break;
            }
            #undef IN_DOM
            md_bitset_free(&s, alloc);
        }
        // Readers, the view and the extracted set agree with the reference
        uint64_t count = 0;
        bool ok = true;
        for (uint32_t i = 0; i < U; ++i) {
            count += ref[i];
            ok = ok && md_bitmask_test_bit(&m, base + i) == (bool)ref[i];
        }
        ok = ok && !md_bitmask_test_bit(&m, base - 1) && !md_bitmask_test_bit(&m, base + U);
        ok = ok && md_bitmask_count(&m) == count && md_bitmask_empty(&m) == (count == 0);
        const uint32_t qb = rnd_range(0, U), qe = rnd_range(qb, U + 1);
        uint64_t qc = 0; for (uint32_t i = qb; i < qe; ++i) qc += ref[i];
        ok = ok && md_bitmask_count_range(&m, base + qb, base + qe) == qc;
        ok = ok && md_bitmask_test_any_range(&m, base + qb, base + qe) == (qc > 0);
        ok = ok && md_bitmask_test_all_range(&m, base + qb, base + qe) == (qc == qe - qb);

        md_bitset_iter_t it = md_bitmask_iter(&m);
        uint64_t seen = 0;
        while (md_bitset_iter_next(&it)) { ok = ok && it.idx >= base && it.idx < base + U && ref[it.idx - base]; ++seen; }
        ok = ok && seen == count;

        const md_bitset_t v = md_bitmask_view(&m);
        ok = ok && md_bitset_validate(v, false) && md_bitset_count(v) == count;
        for (uint32_t i = 0; ok && i < U; ++i) ok = md_bitset_test(v, base + i) == (bool)ref[i];

        md_bitset_t f = md_bitmask_extract(&m, alloc);
        ok = ok && md_bitset_validate(f, true) && md_bitset_equal(f, v) && md_bitset_hash64(f, 3) == md_bitset_hash64(v, 3);
        EXPECT_TRUE(ok);
        md_bitset_free(&f, alloc);
        md_bitmask_free(&m);
    }
    EXPECT_EQ(0, counter.live_bytes); EXPECT_EQ(0, counter.live_allocs);
}

// Words are allocated as they are first written, and not at all for what is never written. Once reserved,
// writing within the reservation never moves them.
UTEST(bitset, mask_allocation) {
    count_alloc_t counter; md_allocator_i* alloc = count_alloc_init(&counter);
    md_bitmask_t m = md_bitmask_create(0, MD_BITSET_MAX_INDEX, alloc);
    EXPECT_EQ(0, counter.live_allocs);
    EXPECT_TRUE(md_bitmask_empty(&m));
    md_bitmask_set_bit(&m, 3000000000u);
    md_bitmask_set_bit(&m, 3000000100u);
    EXPECT_LE(counter.live_bytes, 64 * 8);   // a few words, not 3 billion bits
    EXPECT_EQ(2u, (unsigned)md_bitmask_count(&m));
    md_bitset_t v = md_bitmask_view(&m);
    EXPECT_EQ(3000000000u, v.beg);
    EXPECT_EQ(3000000101u, v.end);

    // Clearing keeps the words, and is O(1): nothing is a member any more
    md_bitmask_clear(&m);
    EXPECT_TRUE(md_bitmask_empty(&m));
    EXPECT_FALSE(md_bitmask_test_bit(&m, 3000000000u));
    EXPECT_TRUE(md_bitset_empty(md_bitmask_view(&m)));
    md_bitmask_free(&m);
    EXPECT_EQ(0, counter.live_bytes);

    md_bitmask_init(&m, 0, 100000, alloc);
    md_bitmask_reserve(&m, 0, 100000);
    const uint64_t* words = m.words;
    for (uint32_t i = 0; i < 100000; i += 7) md_bitmask_set_bit(&m, i);
    md_bitmask_set_range(&m, 0, 100000);
    md_bitmask_clear(&m);
    md_bitmask_set_bit(&m, 99999);
    EXPECT_TRUE(m.words == words);
    // Outside the domain: dropped
    md_bitmask_set_bit(&m, 100000);
    md_bitmask_set_range(&m, 99990, 200000);
    EXPECT_EQ(10u, (unsigned)md_bitmask_count(&m));

    // A narrower domain drops the members outside it
    md_bitmask_set_domain(&m, 99995, 100000);
    EXPECT_EQ(5u, (unsigned)md_bitmask_count(&m));
    EXPECT_FALSE(md_bitmask_test_bit(&m, 99994));
    md_bitmask_free(&m);
    EXPECT_EQ(0, counter.live_bytes); EXPECT_EQ(0, counter.live_allocs);
}

// A RUN and a BITS set with the same members are equal and hash the same: a view of a full mask is BITS
UTEST(bitset, equality_by_members) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_bitmask_t m = md_bitmask_create(0, 1000, alloc);
    md_bitmask_set_range(&m, 100, 300);
    const md_bitset_t v = md_bitmask_view(&m);
    EXPECT_FALSE(md_bitset_is_run(v));
    EXPECT_TRUE(md_bitset_validate(v, false));
    EXPECT_FALSE(md_bitset_validate(v, true));
    EXPECT_TRUE(md_bitset_equal(v, md_bitset_range(100, 300)));
    EXPECT_EQ(md_bitset_hash64(v, 11), md_bitset_hash64(md_bitset_range(100, 300), 11));
    md_bitset_t e = md_bitmask_extract(&m, alloc);
    EXPECT_TRUE(md_bitset_is_run(e));
    md_bitmask_clear_bit(&m, 200);
    EXPECT_FALSE(md_bitset_equal(md_bitmask_view(&m), md_bitset_range(100, 300)));
    EXPECT_NE(md_bitset_hash64(md_bitmask_view(&m), 11), md_bitset_hash64(md_bitset_range(100, 300), 11));
    // Zero words at the ends of what was written are not part of the view
    md_bitmask_clear_range(&m, 100, 299);
    const md_bitset_t one = md_bitmask_view(&m);
    EXPECT_EQ(299u, one.beg);
    EXPECT_EQ(300u, one.end);
    EXPECT_TRUE(md_bitset_validate(one, false));
    md_bitmask_free(&m);
}

UTEST(bitset, replace_assign) {
    count_alloc_t counter; md_allocator_i* alloc = count_alloc_init(&counter);
    static ref_t ra, rb;
    for (int round = 0; round < 200; ++round) {
        ra.base = 8192; rb.base = 8192 + rnd_range(0, 512);
        gen(ra.ref, rnd() % SHAPE_COUNT);
        gen(rb.ref, rnd() % SHAPE_COUNT);
        md_bitset_t b = make(&rb, alloc);
        md_bitset_t a = make(&ra, alloc);
        md_bitset_t expect = md_bitset_empty_set();
        switch (round % 3) {
        case 0: expect = md_bitset_or(a, b, alloc);  md_bitset_replace(&a, md_bitset_or(a, b, alloc), alloc); break;
        case 1: expect = md_bitset_copy(b, alloc);   md_bitset_assign(&a, b, alloc); md_bitset_assign(&a, a, alloc); break;
        case 2: expect = md_bitset_copy(a, alloc);   md_bitset_replace(&a, a, alloc); break;
        }
        EXPECT_TRUE(md_bitset_validate(a, true));
        EXPECT_TRUE(md_bitset_equal(a, expect));
        md_bitset_free(&a, alloc);
        md_bitset_free(&b, alloc);
        md_bitset_free(&expect, alloc);
    }
    EXPECT_EQ(0, counter.live_bytes); EXPECT_EQ(0, counter.live_allocs);
}

// The serialized form is the one the former md_bitfield wrote. These strings were produced by
// md_bitfield_serialize + base64 before it was removed: they must read back, and the same sets must write
// the same bytes (workspace files hold masks in this form).
static bool golden_round_trip(const char* b64, md_bitset_t expected) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const size_t len = strlen(b64);
    uint8_t raw[1024];
    const size_t raw_len = md_base64_decode(raw, b64, len);
    md_bitset_t s = {0};
    bool ok = raw_len > 0 && md_bitset_deserialize(&s, raw, raw_len, alloc);
    ok = ok && md_bitset_equal(s, expected);

    uint8_t out[1024];
    ok = ok && md_bitset_serialize_size_in_bytes(expected) <= sizeof(out);
    const size_t out_len = ok ? md_bitset_serialize(out, expected) : 0;
    char enc[2048];
    const size_t enc_len = out_len ? md_base64_encode(enc, out, out_len) : 0;
    ok = ok && enc_len == len && memcmp(enc, b64, len) == 0;
    md_bitset_free(&s, alloc);
    return ok;
}

UTEST(bitset, serialization_matches_former_bitfield) {
    md_allocator_i* alloc = md_get_heap_allocator();

    const uint32_t bits[] = { 0, 1, 5, 63, 64, 65, 1000, 4095 };
    md_bitset_t a = md_bitset_from_indices(bits, ARRAY_SIZE(bits), alloc);
    EXPECT_TRUE(golden_round_trip("KQMAAAABAAcAIwBgAACAQA/gaAAAAeA0cQQAAAAAgA==", a));

    EXPECT_TRUE(golden_round_trip("JAEAAAD/IAAAAOAtAAQAAAAAAA==", md_bitset_range(0, 32)));

    md_bitset_t c = md_bitset_or(md_bitset_range(600, 2000), md_bitset_range(5000, 5001), alloc);
    EXPECT_TRUE(golden_round_trip("KQQAAQACgAMACQDgAgAA/+BlAOACeeAjAA4BAAAAAAAAAAAAAAAAAAA=", c));

    md_bitset_free(&a, alloc);
    md_bitset_free(&c, alloc);
}

// Round trip over every shape and offset the format can address
UTEST(bitset, serialization_round_trip) {
    md_allocator_i* alloc = md_get_heap_allocator();
    static ref_t r;
    for (int round = 0; round < NUM_ROUNDS; ++round) {
        r.base = pick_base() & 0xFFFFF;  // the format addresses 2^24 bits
        gen(r.ref, round % SHAPE_COUNT);
        md_bitset_t s = make(&r, alloc);
        const size_t cap = md_bitset_serialize_size_in_bytes(s);
        if (md_bitset_empty(s)) { EXPECT_EQ((size_t)0, cap); continue; }
        uint8_t* bytes = (uint8_t*)malloc(cap);
        const size_t n = md_bitset_serialize(bytes, s);
        EXPECT_GT(n, (size_t)0);
        md_bitset_t d = {0};
        EXPECT_TRUE(md_bitset_deserialize(&d, bytes, n, alloc));
        EXPECT_TRUE(md_bitset_equal(d, s));
        free(bytes);
        md_bitset_free(&s, alloc);
        md_bitset_free(&d, alloc);
    }
}
