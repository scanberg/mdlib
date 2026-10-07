#pragma once

#include <core/md_intrinsics.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;

// Indices are below 2^32 - 64, so the base of a 64-index window plus 64 never wraps
#define MD_BITSET_MAX_INDEX 0xFFFFFFC0u

// Two types:
//
//   md_bitset_t   A read-only set of indices: the type sets are passed around and stored as. 16 bytes, copied
//                 freely. A range is a set without any memory (md_bitset_range).
//   md_bitmask_t  A dense, mutable set over a domain: the type a set is edited in place in (a selection, a hover
//                 highlight, the frames evaluated so far) or built in, index by index. md_bitmask_view hands out
//                 its members as an md_bitset_t without copying.

// ---------------------------------------------------------------------------------------------------------
// md_bitset_t
//
// Two forms:
//   RUN   words == NULL   every index in [beg, end) is a member; beg == end is the empty set
//   BITS  words != NULL   bit k of words[] is index (beg & ~63) + k. The words cover [beg & ~63, align_up(end, 64)),
//                         bits outside [beg, end) are zero.
// In either form, beg and end - 1 are members (unless the set is empty), so empty, first and last are O(1).
//
// A set does not own anything. It refers to words produced into an allocator by the operations below, which are
// never written again, or to the words of a mask (md_bitmask_view), valid until that mask is next modified,
// cleared or freed. Copying the struct copies the reference. Whatever keeps a set longer than what it refers to
// lives takes its own copy (md_bitset_copy) into an allocator that lives long enough.
//
// The sets the producers below return are also tidy: a set whose members are all of [beg, end) is a RUN, and the
// empty set is {0}. Readers do not rely on it: equality and hashing go by the members, not the form.
typedef struct md_bitset_t {
    uint32_t beg;
    uint32_t end;
    const uint64_t* words;
} md_bitset_t;

// One 64-index window at a time; a loop over the members costs one ctz per member and one fetch per window.
typedef struct md_bitset_iter_t {
    const uint64_t* words;  // NULL: every index in [beg, end) is a member
    uint32_t word_base;     // index of bit 0 of words[0]
    uint32_t beg;           // first index considered
    uint32_t end;           // iteration stops before this index
    uint32_t base;          // first index of the current window, multiple of 64
    uint64_t mask;          // remaining members of the current window
    uint32_t idx;           // current member
} md_bitset_iter_t;

// GPU form. word_off indexes a packed u64 word buffer (u32 index 2*word_off when read as u32 pairs).
// word_off == MD_BITSET_RUN marks a RUN; count is the member count.
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

// The members of the window [base, base + 64) of a set given by its parts (base a multiple of 64)
static inline uint64_t md_bitset__window(const uint64_t* words, uint32_t word_base, uint32_t beg, uint32_t end, uint32_t base) {
    if (base >= end || (uint64_t)base + 64 <= beg) return 0;
    if (words) return words[(base - word_base) >> 6];
    uint64_t m = ~0ull;
    if (beg > base) m &= ~0ull << (beg - base);
    if (end - base < 64) m &= ~(~0ull << (end - base));
    return m;
}

// --- Inline readers ---

static inline md_bitset_t md_bitset_empty_set(void) {
    md_bitset_t s = {0, 0, NULL};
    return s;
}

// The set of all indices in [beg, end). Needs no memory.
static inline md_bitset_t md_bitset_range(uint32_t beg, uint32_t end) {
    md_bitset_t s = {0, 0, NULL};
    if (beg < end) { s.beg = beg; s.end = end; }
    return s;
}

static inline bool md_bitset_empty(md_bitset_t s) { return s.beg >= s.end; }
static inline bool md_bitset_is_run(md_bitset_t s) { return s.beg < s.end && s.words == NULL; }

// First member and one past the last (both 0 for the empty set)
static inline uint32_t md_bitset_first(md_bitset_t s) { return s.beg < s.end ? s.beg : 0; }
static inline uint32_t md_bitset_last_end(md_bitset_t s) { return s.beg < s.end ? s.end : 0; }

// Number of words a BITS set refers to (0 for a RUN)
static inline size_t md_bitset_num_words(md_bitset_t s) {
    return (s.words && s.beg < s.end) ? (size_t)((((uint64_t)s.end + 63) >> 6) - (s.beg >> 6)) : 0;
}

static inline bool md_bitset_test(md_bitset_t s, uint32_t i) {
    if (i - s.beg >= s.end - s.beg || s.beg >= s.end) return false;  // outside the span (unsigned wrap)
    if (!s.words) return true;
    const uint32_t k = i - (s.beg & ~63u);
    return (s.words[k >> 6] >> (k & 63)) & 1;
}

// Members of the window [base, base + 64) as a mask; base must be a multiple of 64.
static inline uint64_t md_bitset_word(md_bitset_t s, uint32_t base) {
    if (s.beg >= s.end) return 0;
    return md_bitset__window(s.words, s.beg & ~63u, s.beg, s.end, base);
}

// Iterate the members within [beg, end)
static inline md_bitset_iter_t md_bitset_iter_range(md_bitset_t s, uint32_t beg, uint32_t end) {
    md_bitset_iter_t it = {0};
    if (s.beg >= s.end) return it;
    beg = beg > s.beg ? beg : s.beg;
    end = end < s.end ? end : s.end;
    if (beg < end) {
        it.words = s.words;
        it.word_base = s.beg & ~63u;
        it.beg  = beg;
        it.end  = end;
        it.base = beg & ~63u;
        it.mask = md_bitset__window(it.words, it.word_base, beg, end, it.base) & (~0ull << (beg - it.base));
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
        it->mask = md_bitset__window(it->words, it->word_base, it->beg, it->end, it->base);
    }
    const uint32_t i = it->base + (uint32_t)ctz64(it->mask);
    if (i >= it->end) { it->mask = 0; return false; }
    it->mask &= it->mask - 1;
    it->idx = i;
    return true;
}

// Reference implementation of the shader side test against the packed form.
static inline bool md_bitset_gpu_test(md_bitset_gpu_t s, const uint64_t* words, uint32_t i) {
    if (i - s.beg >= s.end - s.beg || s.beg >= s.end) return false;
    if (s.word_off == MD_BITSET_RUN) return true;
    const uint32_t k = i - (s.beg & ~63u);
    return (words[s.word_off + (k >> 6)] >> (k & 63)) & 1;
}

// --- Readers ---

uint64_t md_bitset_count           (md_bitset_t s);  // O(1) for a RUN, a popcount over the words otherwise
uint64_t md_bitset_count_range     (md_bitset_t s, uint32_t beg, uint32_t end);
bool     md_bitset_test_any_range  (md_bitset_t s, uint32_t beg, uint32_t end);
bool     md_bitset_test_all_range  (md_bitset_t s, uint32_t beg, uint32_t end);

// By the members: a RUN and a BITS set with the same members are equal and hash the same
bool     md_bitset_equal (md_bitset_t a, md_bitset_t b);
uint64_t md_bitset_hash64(md_bitset_t s, uint64_t seed);

// Writes up to cap member indices in ascending order, returns the number written.
size_t   md_bitset_extract_indices(uint32_t* out, size_t cap, md_bitset_t s);
// The same, for the members remaining in an iterator (e.g. md_bitset_iter_range), into int32 storage.
size_t   md_bitset_iter_extract_indices(int32_t* out, size_t cap, md_bitset_iter_t it);

// Checks the invariants above (for tests and debug asserts). tidy additionally requires the tidy form.
bool     md_bitset_validate(md_bitset_t s, bool tidy);

// --- Producers: the result is tidy, its words (if any) allocated with alloc ---

// From a dense mask: bit k of words[] is index base + k. base must be a multiple of 64.
md_bitset_t md_bitset_from_words(const uint64_t* words, size_t num_words, uint32_t base, struct md_allocator_i* alloc);

// From indices, in any order, duplicates allowed.
md_bitset_t md_bitset_from_indices(const uint32_t* indices, size_t num_indices, struct md_allocator_i* alloc);

md_bitset_t md_bitset_and   (md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);
md_bitset_t md_bitset_or    (md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);
md_bitset_t md_bitset_andnot(md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);
md_bitset_t md_bitset_xor   (md_bitset_t a, md_bitset_t b, struct md_allocator_i* alloc);

// Complement within [beg, end). A set carries no universe, so the caller supplies it.
md_bitset_t md_bitset_not   (md_bitset_t s, uint32_t beg, uint32_t end, struct md_allocator_i* alloc);

// An independent copy into alloc (a RUN copies without allocating). This is how a set is kept.
md_bitset_t md_bitset_copy  (md_bitset_t s, struct md_allocator_i* alloc);

// Release the words of a set a producer allocated with alloc (never a view of a mask). Only meaningful for
// allocators that free individually; sets in arenas go with the arena. Resets *s to empty.
void md_bitset_free(md_bitset_t* s, struct md_allocator_i* alloc);

// For a set kept in one allocator and replaced over time (a stored selection, a filter result):
// replace takes ownership of value, which must have been produced with alloc; assign stores a copy of src.
// Either releases the previous value, which therefore must have been produced with alloc.
void md_bitset_replace(md_bitset_t* dst, md_bitset_t value, struct md_allocator_i* alloc);
void md_bitset_assign (md_bitset_t* dst, md_bitset_t src,   struct md_allocator_i* alloc);

// --- GPU packing ---

// Number of u64 words needed to pack the sets.
size_t md_bitset_pack_num_words(const md_bitset_t* sets, size_t num_sets);

// Write num_sets descriptors to out_sets and their words, back to back, to out_words.
void   md_bitset_pack(md_bitset_gpu_t* out_sets, uint64_t* out_words, const md_bitset_t* sets, size_t num_sets);

// --- Serialization ---
// The byte format is the one the former md_bitfield used (512-bit blocks, FastLZ compressed), so masks stored
// before md_bitset replaced it (e.g. in workspace files) still read, and the same set writes the same bytes.

size_t md_bitset_serialize_size_in_bytes(md_bitset_t s);           // upper bound, 0 for an empty set
size_t md_bitset_serialize(void* dst, md_bitset_t s);              // returns bytes written, 0 on failure or empty
bool   md_bitset_deserialize(md_bitset_t* dst, const void* src, size_t num_bytes, struct md_allocator_i* alloc);

// ---------------------------------------------------------------------------------------------------------
// md_bitmask_t
//
// A dense, mutable set over a domain [beg, end), fixed at init: members outside it are dropped, so a mask over a
// context's span restricts to it. Edits are in place: setting or testing an index is O(1) and allocates nothing
// once the words are there.
//
// The mask holds the words it has written, contiguously, and allocates them as they are first written (or all of
// the domain up front, md_bitmask_reserve). Clearing is O(1), and counting, iterating, viewing and extracting cost
// O(words written), not O(domain). Growing to words not yet allocated moves the words: a view, or an iterator over
// one, is invalidated by any change to the mask.
//
// A mask is not synchronized. A reader on another thread takes a copy (md_bitmask_extract) under the lock the
// writer holds.
typedef struct md_bitmask_t {
    uint64_t* words;                // words[k] holds the indices [64 (g0 + k), 64 (g0 + k + 1))
    struct md_allocator_i* alloc;
    uint32_t beg;                   // domain
    uint32_t end;
    uint32_t g0;                    // word index of words[0]
    uint32_t cap;                   // words allocated
    uint32_t lo;                    // words [lo, hi) are written; any other allocated word is undefined, and zero
    uint32_t hi;                    // as far as the members go. lo == hi: no members
} md_bitmask_t;

void         md_bitmask_init  (md_bitmask_t* m, uint32_t beg, uint32_t end, struct md_allocator_i* alloc);
md_bitmask_t md_bitmask_create(uint32_t beg, uint32_t end, struct md_allocator_i* alloc);
void         md_bitmask_free  (md_bitmask_t* m);

// Allocates the words of [beg, end) (within the domain) up front, so writing within it never moves the words
void md_bitmask_reserve   (md_bitmask_t* m, uint32_t beg, uint32_t end);
// Changes the domain, dropping the members outside the new one
void md_bitmask_set_domain(md_bitmask_t* m, uint32_t beg, uint32_t end);

void md_bitmask_clear(md_bitmask_t* m);  // no members, O(1)

void md_bitmask_set_bit     (md_bitmask_t* m, uint32_t idx);
void md_bitmask_clear_bit   (md_bitmask_t* m, uint32_t idx);
void md_bitmask_set_range   (md_bitmask_t* m, uint32_t beg, uint32_t end);
void md_bitmask_clear_range (md_bitmask_t* m, uint32_t beg, uint32_t end);
void md_bitmask_invert_range(md_bitmask_t* m, uint32_t beg, uint32_t end);   // complement within [beg, end)
void md_bitmask_set_indices (md_bitmask_t* m, const uint32_t* indices, size_t num_indices);

// With a set as the operand. The set may be a view of the mask itself (m |= view(m) is fine).
void md_bitmask_assign(md_bitmask_t* m, md_bitset_t s);   // m = s (within the domain)
void md_bitmask_or    (md_bitmask_t* m, md_bitset_t s);   // m |= s
void md_bitmask_and   (md_bitmask_t* m, md_bitset_t s);   // m &= s
void md_bitmask_andnot(md_bitmask_t* m, md_bitset_t s);   // m &= ~s
void md_bitmask_xor   (md_bitmask_t* m, md_bitset_t s);   // m ^= s

static inline bool md_bitmask_test_bit(const md_bitmask_t* m, uint32_t idx) {
    const uint32_t g = idx >> 6;
    return g - m->lo < m->hi - m->lo && ((m->words[g - m->g0] >> (idx & 63)) & 1);
}

bool     md_bitmask_empty           (const md_bitmask_t* m);
uint64_t md_bitmask_count           (const md_bitmask_t* m);
uint64_t md_bitmask_count_range     (const md_bitmask_t* m, uint32_t beg, uint32_t end);
bool     md_bitmask_test_any_range  (const md_bitmask_t* m, uint32_t beg, uint32_t end);
bool     md_bitmask_test_all_range  (const md_bitmask_t* m, uint32_t beg, uint32_t end);

// The members as a set referring to the mask's words, valid until the mask is next modified, cleared or freed.
// Its span is trimmed to the first and last member; it is not tidy (a full span is not turned into a RUN).
md_bitset_t md_bitmask_view   (const md_bitmask_t* m);
// The members as an independent, tidy set allocated with alloc
md_bitset_t md_bitmask_extract(const md_bitmask_t* m, struct md_allocator_i* alloc);

static inline md_bitset_iter_t md_bitmask_iter(const md_bitmask_t* m) {
    return md_bitset_iter(md_bitmask_view(m));
}

#ifdef __cplusplus
}
#endif
