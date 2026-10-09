#pragma once

#include <core/md_intrinsics.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitfield_t;

// md_bitset_t: an immutable set of indices in [0, 2^32 - 64), held as a single span.
//
// The span [beg, end) is trimmed to the first and last member, and the set takes one of three forms:
//
//   EMPTY   beg == end (canonically {0, 0, NULL})
//   RUN     words == NULL, beg < end      every index in [beg, end) is a member; no payload
//   BITS    words != NULL                 bit k of words[] is index (beg & ~63) + k; words cover
//                                         [beg & ~63, align_up(end, 64)), bits outside [beg, end) are zero
//
// The form is canonical: a set whose members are all of [beg, end) is always a RUN, never BITS, and an
// empty set is always {0}. Two sets are equal exactly when beg, end and the words are equal.
//
// The value is 16 bytes and is copied freely. A BITS payload is one allocation of 1 + num_words u64:
// words[-1] holds the member count, so md_bitset_count is O(1) for every form.
//
// LIFETIME. A set does not own or know its allocator. Producers take the allocator of the result; readers
// take the set alone. Copying the struct aliases the payload, which is fine while the payload's allocator
// outlives the copy (e.g. a frame value referring to a set in the IR arena). A value that must outlive
// its payload's allocator is deep copied with md_bitset_copy. The payload is never written after the set
// is produced.
typedef struct md_bitset_t {
    uint32_t beg;
    uint32_t end;
    const uint64_t* words;
} md_bitset_t;

// The iterator caches one 64-index window, so a loop over the members costs one ctz per member and one
// window fetch per 64 indices, whatever the form.
typedef struct md_bitset_iter_t {
    md_bitset_t set;
    uint64_t mask;   // remaining members of the current window
    uint32_t base;   // first index of the current window, multiple of 64
    uint32_t end;    // iteration stops before this index
    uint32_t idx;    // current member
} md_bitset_iter_t;

// GPU / serialized form. Position independent: word_off indexes a packed u64 word buffer (u32 index
// 2*word_off when read as u32 pairs). word_off == MD_BITSET_RUN marks a RUN; count is the member count.
#define MD_BITSET_RUN 0xFFFFFFFFu
typedef struct md_bitset_gpu_t {
    uint32_t beg;
    uint32_t end;
    uint32_t word_off;
    uint32_t count;
} md_bitset_gpu_t;

#ifdef __cplusplus
extern "C" {
#endif

// --- Inline readers ---

static inline md_bitset_t md_bitset_empty_set(void) {
    md_bitset_t s = {0, 0, NULL};
    return s;
}

// The set of all indices in [beg, end). Needs no allocation.
static inline md_bitset_t md_bitset_range(uint32_t beg, uint32_t end) {
    md_bitset_t s = {0, 0, NULL};
    if (beg < end) { s.beg = beg; s.end = end; }
    return s;
}

static inline bool md_bitset_empty(md_bitset_t s) { return s.beg == s.end; }
static inline bool md_bitset_is_run(md_bitset_t s) { return s.beg < s.end && s.words == NULL; }

static inline uint64_t md_bitset_count(md_bitset_t s) {
    return s.words ? s.words[-1] : (uint64_t)(s.end - s.beg);
}

// Number of payload words (0 for EMPTY and RUN)
static inline size_t md_bitset_num_words(md_bitset_t s) {
    return s.words ? (size_t)((((uint64_t)s.end + 63) >> 6) - (s.beg >> 6)) : 0;
}

static inline bool md_bitset_test(md_bitset_t s, uint32_t i) {
    if (i - s.beg >= s.end - s.beg) return false;  // outside the span, EMPTY included (unsigned wrap)
    if (!s.words) return true;
    const uint32_t k = i - (s.beg & ~63u);
    return (s.words[k >> 6] >> (k & 63)) & 1;
}

// Membership of the 64 indices [base, base + 64) as a mask; base must be a multiple of 64.
// Uniform over all three forms, so a consumer can process a set a window at a time.
static inline uint64_t md_bitset_word(md_bitset_t s, uint32_t base) {
    if (base >= s.end || (uint64_t)base + 64 <= s.beg) return 0;
    if (s.words) return s.words[(base - (s.beg & ~63u)) >> 6];
    uint64_t m = ~0ull;
    if (s.beg > base) m &= ~0ull << (s.beg - base);
    if (s.end - base < 64) m &= ~(~0ull << (s.end - base));
    return m;
}

static inline md_bitset_iter_t md_bitset_iter_range(md_bitset_t s, uint32_t beg, uint32_t end) {
    md_bitset_iter_t it = {s, 0, 0, 0, 0};
    beg = beg > s.beg ? beg : s.beg;
    end = end < s.end ? end : s.end;
    if (beg < end) {
        it.base = beg & ~63u;
        it.end  = end;
        it.mask = md_bitset_word(s, it.base) & (~0ull << (beg - it.base));
    }
    return it;
}

static inline md_bitset_iter_t md_bitset_iter(md_bitset_t s) {
    return md_bitset_iter_range(s, s.beg, s.end);
}

// Advance to the next member; returns false when there are none left.
//   md_bitset_iter_t it = md_bitset_iter(s);
//   while (md_bitset_iter_next(&it)) { uint32_t i = it.idx; ... }
static inline bool md_bitset_iter_next(md_bitset_iter_t* it) {
    while (!it->mask) {
        if (it->end - it->base <= 64) return false;  // also true for an empty iterator (end == base == 0)
        it->base += 64;
        it->mask = md_bitset_word(it->set, it->base);
    }
    const uint32_t i = it->base + (uint32_t)ctz64(it->mask);
    if (i >= it->end) { it->mask = 0; return false; }
    it->mask &= it->mask - 1;
    it->idx = i;
    return true;
}

// Reference implementation of the shader side test against the packed form.
static inline bool md_bitset_gpu_test(md_bitset_gpu_t s, const uint64_t* words, uint32_t i) {
    if (i - s.beg >= s.end - s.beg) return false;
    if (s.word_off == MD_BITSET_RUN) return true;
    const uint32_t k = i - (s.beg & ~63u);
    return (words[s.word_off + (k >> 6)] >> (k & 63)) & 1;
}

// --- Producers: the result's payload (if any) is allocated with alloc ---

// Build from a dense mask: bit k of words[] is index base + k. base must be a multiple of 64.
md_bitset_t md_bitset_from_words(const uint64_t* words, size_t num_words, uint32_t base, struct md_allocator_i* alloc);

// Build from indices, in any order, duplicates allowed.
md_bitset_t md_bitset_from_indices(const uint32_t* indices, size_t num_indices, struct md_allocator_i* alloc);

md_bitset_t md_bitset_and   (md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);
md_bitset_t md_bitset_or    (md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);
md_bitset_t md_bitset_andnot(md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);
md_bitset_t md_bitset_xor   (md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);

// Complement within the domain [beg, end). A set carries no universe, so the caller supplies it.
md_bitset_t md_bitset_not   (md_bitset_t s, uint32_t beg, uint32_t end, struct md_allocator_i* alloc);

// Deep copy into alloc (a RUN or EMPTY copies without allocating).
md_bitset_t md_bitset_copy  (md_bitset_t s, struct md_allocator_i* alloc);

// Release the payload of a set produced with alloc. Only meaningful for allocators that free
// individually; sets in arenas go away with the arena. Resets *s to EMPTY.
void md_bitset_free(md_bitset_t* s, struct md_allocator_i* alloc);

// --- Builder ---
//
// The set itself is immutable, so anything that accumulates (a procedure setting the atoms it matches, a union
// over many sets, a visualization collecting atoms) goes through a builder and produces the set once, at the
// end. Folding with md_bitset_or instead allocates every intermediate: n operands cost O(n * span) time and,
// in an arena, O(n * span) memory.
//
// A builder is dense over a DOMAIN [beg, end) fixed at init, the universe or the span of a context. Members
// outside the domain are dropped, so a builder over a context's span restricts to that context for free.
// It tracks the words it has touched: finishing and resetting cost O(touched), not O(domain), which is what
// makes one builder reusable across many small contexts. Its words come from scratch, typically a temp arena
// distinct from the allocator of the result.
typedef struct md_bitset_builder_t {
    uint64_t* words;    // words[k] holds the indices [base + 64 k, base + 64 (k + 1))
    uint32_t  base;     // beg & ~63
    uint32_t  beg;      // domain
    uint32_t  end;
    uint32_t  num_words;
    uint32_t  lo;       // touched words [lo, hi), empty when lo >= hi
    uint32_t  hi;
} md_bitset_builder_t;

#define MD_BITSET_MAX_INDEX 0xFFFFFFC0u  // indices are < 2^32 - 64

void        md_bitset_builder_init   (md_bitset_builder_t* b, uint32_t beg, uint32_t end, struct md_allocator_i* scratch);
void        md_bitset_builder_free   (md_bitset_builder_t* b, struct md_allocator_i* scratch);  // only for scratch which frees individually
void        md_bitset_builder_reset  (md_bitset_builder_t* b);                                  // empty again, O(touched)

void        md_bitset_builder_set        (md_bitset_builder_t* b, uint32_t i);
void        md_bitset_builder_set_range  (md_bitset_builder_t* b, uint32_t beg, uint32_t end);
void        md_bitset_builder_set_indices(md_bitset_builder_t* b, const uint32_t* indices, size_t num_indices);
void        md_bitset_builder_or         (md_bitset_builder_t* b, md_bitset_t s);   // b |= s
void        md_bitset_builder_and        (md_bitset_builder_t* b, md_bitset_t s);   // b &= s
void        md_bitset_builder_andnot     (md_bitset_builder_t* b, md_bitset_t s);   // b &= ~s

bool        md_bitset_builder_test   (const md_bitset_builder_t* b, uint32_t i);
// The set built so far, in canonical form, its payload (if any) in alloc. The builder is left as it is.
md_bitset_t md_bitset_builder_finish (const md_bitset_builder_t* b, struct md_allocator_i* alloc);

// --- Readers (non-inline) ---

uint64_t md_bitset_count_range   (md_bitset_t s, uint32_t beg, uint32_t end);
bool     md_bitset_test_any_range(md_bitset_t s, uint32_t beg, uint32_t end);
bool     md_bitset_test_all_range(md_bitset_t s, uint32_t beg, uint32_t end);

bool     md_bitset_equal (md_bitset_t a, md_bitset_t b);
uint64_t md_bitset_hash64(md_bitset_t s, uint64_t seed);

// Writes up to cap member indices in ascending order, returns the number written.
size_t   md_bitset_extract_indices(uint32_t* out, size_t cap, md_bitset_t s);

// Checks the canonical-form invariants (for tests and debug asserts).
bool     md_bitset_validate(md_bitset_t s);

// --- md_bitfield interop (migration) ---

md_bitset_t md_bitset_from_bitfield(const struct md_bitfield_t* bf, struct md_allocator_i* alloc);
void        md_bitset_to_bitfield  (struct md_bitfield_t* dst, md_bitset_t s);

// --- GPU packing ---

// Number of u64 words needed to pack the sets.
size_t md_bitset_pack_num_words(const md_bitset_t* sets, size_t num_sets);

// Write num_sets descriptors to out_sets and their payloads, back to back, to out_words.
void   md_bitset_pack(md_bitset_gpu_t* out_sets, uint64_t* out_words, const md_bitset_t* sets, size_t num_sets);

#ifdef __cplusplus
}
#endif
