#include <core/md_pair_frames.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_os.h>
#include <core/md_sort.h>

#include <math.h>

// Rows, one per entry ever true, in blocks of BLOCK_ROWS that never move once made: an entry, then its frames' bits.
// The blocks are reached through a directory of DIR_SIZE pages of DIR_SIZE blocks each, made as they are needed, so
// nothing a reader reaches moves while a writer adds. Writers take the lock; readers do not, and read what is
// published: count after the rows it counts, a frame's bit in added after its rows.
#define BLOCK_ROWS 1024u
#define DIR_SIZE   1024u

typedef struct block_t {
    uint32_t entry[BLOCK_ROWS];
    uint64_t bits[];            // [BLOCK_ROWS][words]
} block_t;

struct md_pair_frames_t {
    struct md_allocator_i* alloc;
    uint32_t            rows, cols, num_frames, words;
    uint64_t*           added;      // [words]
    block_t**           dir[DIR_SIZE];
    volatile uint64_t   count;      // the rows made
    // Entry + 1 to its row, by open addressing; the writers' alone
    uint32_t*           map_key;
    uint32_t*           map_row;
    uint32_t            map_cap;    // a power of two, or 0
    md_mutex_t          lock;
};

static inline size_t block_bytes(uint32_t words) {
    return sizeof(block_t) + (size_t)BLOCK_ROWS * words * sizeof(uint64_t);
}

static inline block_t* block_of(const md_pair_frames_t* pf, uint64_t row) {
    const uint64_t b = row / BLOCK_ROWS;
    return pf->dir[b / DIR_SIZE][b % DIR_SIZE];
}

static inline volatile uint64_t* row_bits(const md_pair_frames_t* pf, uint64_t row) {
    return (volatile uint64_t*)block_of(pf, row)->bits + (size_t)(row % BLOCK_ROWS) * pf->words;
}

static inline uint32_t slot_of(uint32_t key, uint32_t cap) {
    return (uint32_t)(((uint64_t)key * 0x9E3779B97F4A7C15ull) >> 32) & (cap - 1);
}

static bool map_find(const md_pair_frames_t* pf, uint32_t entry, uint32_t* row) {
    if (!pf->map_cap) return false;
    const uint32_t key = entry + 1;
    for (uint32_t s = slot_of(key, pf->map_cap);; s = (s + 1) & (pf->map_cap - 1)) {
        if (pf->map_key[s] == key) { *row = pf->map_row[s]; return true; }
        if (pf->map_key[s] == 0) return false;
    }
}

static void map_put(uint32_t* keys, uint32_t* rows, uint32_t cap, uint32_t key, uint32_t row) {
    uint32_t s = slot_of(key, cap);
    while (keys[s]) s = (s + 1) & (cap - 1);
    keys[s] = key;
    rows[s] = row;
}

// Room in the map for one more, at most half full
static bool map_reserve(md_pair_frames_t* pf, uint64_t num) {
    if ((num + 1) * 2 <= pf->map_cap) return true;
    const uint32_t cap = pf->map_cap ? pf->map_cap * 2 : 1024;
    if (cap < pf->map_cap) return false;
    uint32_t* keys = md_alloc(pf->alloc, (size_t)cap * sizeof(uint32_t));
    uint32_t* rows = md_alloc(pf->alloc, (size_t)cap * sizeof(uint32_t));
    if (!keys || !rows) {
        if (keys) md_free(pf->alloc, keys, (size_t)cap * sizeof(uint32_t));
        if (rows) md_free(pf->alloc, rows, (size_t)cap * sizeof(uint32_t));
        return false;
    }
    md_memset(keys, 0, (size_t)cap * sizeof(uint32_t));
    for (uint32_t s = 0; s < pf->map_cap; ++s) if (pf->map_key[s]) map_put(keys, rows, cap, pf->map_key[s], pf->map_row[s]);
    if (pf->map_cap) {
        md_free(pf->alloc, pf->map_key, (size_t)pf->map_cap * sizeof(uint32_t));
        md_free(pf->alloc, pf->map_row, (size_t)pf->map_cap * sizeof(uint32_t));
    }
    pf->map_key = keys;
    pf->map_row = rows;
    pf->map_cap = cap;
    return true;
}

// Row `row`, the next, for entry: its block made when it is the first of one, zeroed; not counted yet
static bool make_row(md_pair_frames_t* pf, uint64_t row, uint32_t entry) {
    if (!map_reserve(pf, row)) return false;
    const uint64_t b = row / BLOCK_ROWS;
    if (b / DIR_SIZE >= DIR_SIZE) return false;
    block_t** page = pf->dir[b / DIR_SIZE];
    if (!page) {
        page = md_alloc(pf->alloc, DIR_SIZE * sizeof(block_t*));
        if (!page) return false;
        md_memset(page, 0, DIR_SIZE * sizeof(block_t*));
        pf->dir[b / DIR_SIZE] = page;
    }
    if (!page[b % DIR_SIZE]) {
        block_t* blk = md_alloc(pf->alloc, block_bytes(pf->words));
        if (!blk) return false;
        md_memset(blk, 0, block_bytes(pf->words));
        page[b % DIR_SIZE] = blk;
    }
    page[b % DIR_SIZE]->entry[row % BLOCK_ROWS] = entry;
    map_put(pf->map_key, pf->map_row, pf->map_cap, entry + 1, (uint32_t)row);
    return true;
}

md_pair_frames_t* md_pair_frames_create(uint32_t rows, uint32_t cols, uint32_t num_frames, struct md_allocator_i* alloc) {
    const uint64_t total = (uint64_t)rows * cols;
    if (!alloc || total == 0 || total >= 0xFFFFFFFFull || num_frames == 0) return NULL;
    md_pair_frames_t* pf = md_alloc(alloc, sizeof(md_pair_frames_t));
    if (!pf) return NULL;
    md_memset(pf, 0, sizeof(*pf));
    pf->alloc = alloc;
    pf->rows = rows;
    pf->cols = cols;
    pf->num_frames = num_frames;
    pf->words = (num_frames + 63) / 64;
    pf->added = md_alloc(alloc, (size_t)pf->words * sizeof(uint64_t));
    if (!pf->added || !md_mutex_init(&pf->lock)) {
        if (pf->added) md_free(alloc, pf->added, (size_t)pf->words * sizeof(uint64_t));
        md_free(alloc, pf, sizeof(md_pair_frames_t));
        return NULL;
    }
    md_memset(pf->added, 0, (size_t)pf->words * sizeof(uint64_t));
    return pf;
}

void md_pair_frames_free(md_pair_frames_t* pf) {
    if (!pf) return;
    struct md_allocator_i* alloc = pf->alloc;
    for (uint32_t p = 0; p < DIR_SIZE; ++p) {
        if (!pf->dir[p]) continue;
        for (uint32_t b = 0; b < DIR_SIZE; ++b) if (pf->dir[p][b]) md_free(alloc, pf->dir[p][b], block_bytes(pf->words));
        md_free(alloc, pf->dir[p], DIR_SIZE * sizeof(block_t*));
    }
    if (pf->map_cap) {
        md_free(alloc, pf->map_key, (size_t)pf->map_cap * sizeof(uint32_t));
        md_free(alloc, pf->map_row, (size_t)pf->map_cap * sizeof(uint32_t));
    }
    md_free(alloc, pf->added, (size_t)pf->words * sizeof(uint64_t));
    md_mutex_destroy(&pf->lock);
    md_free(alloc, pf, sizeof(md_pair_frames_t));
}

md_pair_frames_t* md_pair_frames_copy(const md_pair_frames_t* src, struct md_allocator_i* alloc) {
    if (!src) return NULL;
    md_pair_frames_t* pf = md_pair_frames_create(src->rows, src->cols, src->num_frames, alloc);
    if (!pf) return NULL;
    md_pair_frames_t* s = (md_pair_frames_t*)src;    // its lock: no frame is added while it is copied
    md_mutex_lock(&s->lock);
    const uint64_t n = md_atomic_load64(&s->count);
    bool ok = true;
    for (uint64_t r = 0; r < n && ok; ++r) {
        ok = make_row(pf, r, block_of(s, r)->entry[r % BLOCK_ROWS]);
        if (ok) md_memcpy((void*)row_bits(pf, r), (const void*)row_bits(s, r), (size_t)pf->words * sizeof(uint64_t));
    }
    md_memcpy(pf->added, s->added, (size_t)pf->words * sizeof(uint64_t));
    md_mutex_unlock(&s->lock);
    if (!ok) { md_pair_frames_free(pf); return NULL; }
    pf->count = n;
    return pf;
}

uint32_t md_pair_frames_rows(const md_pair_frames_t* pf) { return pf ? pf->rows : 0; }
uint32_t md_pair_frames_cols(const md_pair_frames_t* pf) { return pf ? pf->cols : 0; }
uint32_t md_pair_frames_num_frames(const md_pair_frames_t* pf) { return pf ? pf->num_frames : 0; }

bool md_pair_frames_add(md_pair_frames_t* pf, uint32_t frame, const uint32_t* entries, size_t count) {
    if (!pf || frame >= pf->num_frames || (count && !entries)) return false;
    const uint64_t total = (uint64_t)pf->rows * pf->cols;
    const uint32_t w = frame / 64;
    const uint64_t bit = 1ull << (frame % 64);
    md_mutex_lock(&pf->lock);
    volatile uint64_t* aw = (volatile uint64_t*)&pf->added[w];
    if (md_atomic_load64(aw) & bit) {
        md_mutex_unlock(&pf->lock);
        return false;
    }
    uint64_t n = md_atomic_load64(&pf->count);
    bool ok = true;
    for (size_t i = 0; i < count; ++i) {
        const uint32_t e = entries[i];
        if (e >= total) continue;
        uint32_t r;
        if (!map_find(pf, e, &r)) {
            if (!make_row(pf, n, e)) { ok = false; break; }
            r = (uint32_t)n++;
        }
        volatile uint64_t* rw = row_bits(pf, r) + w;
        md_atomic_store64(rw, md_atomic_load64(rw) | bit);
    }
    // The rows, then the frame: a reader that sees the frame sees all of it
    md_atomic_store64(&pf->count, n);
    if (ok) md_atomic_store64(aw, md_atomic_load64(aw) | bit);
    md_mutex_unlock(&pf->lock);
    return ok;
}

bool md_pair_frames_has(const md_pair_frames_t* pf, uint32_t frame) {
    if (!pf || frame >= pf->num_frames) return false;
    return (md_atomic_load64((const volatile uint64_t*)&pf->added[frame / 64]) >> (frame % 64)) & 1;
}

md_bitset_t md_pair_frames_added(const md_pair_frames_t* pf, struct md_allocator_i* alloc) {
    if (!pf || !alloc) return md_bitset_empty_set();
    uint64_t* words = md_alloc(alloc, (size_t)pf->words * sizeof(uint64_t));
    if (!words) return md_bitset_empty_set();
    for (uint32_t w = 0; w < pf->words; ++w) words[w] = md_atomic_load64((const volatile uint64_t*)&pf->added[w]);
    const md_bitset_t s = md_bitset_from_words(words, pf->words, 0, alloc);
    md_free(alloc, words, (size_t)pf->words * sizeof(uint64_t));
    return s;
}

size_t md_pair_frames_count(const md_pair_frames_t* pf) {
    return pf ? (size_t)md_atomic_load64(&pf->count) : 0;
}

bool md_pair_frames_frame(uint8_t* out, const md_pair_frames_t* pf, uint32_t frame) {
    if (!out || !md_pair_frames_has(pf, frame)) return false;
    md_memset(out, 0, (size_t)pf->rows * pf->cols);
    const uint64_t n = md_atomic_load64(&pf->count);
    const uint32_t w = frame / 64;
    const uint64_t bit = 1ull << (frame % 64);
    for (uint64_t r = 0; r < n; ++r) {
        if (md_atomic_load64(row_bits(pf, r) + w) & bit) out[block_of(pf, r)->entry[r % BLOCK_ROWS]] = 1;
    }
    return true;
}

// The frames counted, into sel[words]: those added among frames, as they are now. A row's bits are read through them,
// so a frame added meanwhile is in neither the counts nor what they are divided by. Returns how many, with the words
// holding them in [lo, hi).
static size_t select_frames(uint64_t* sel, uint32_t* lo, uint32_t* hi, const md_pair_frames_t* pf, md_bitset_t frames) {
    size_t num = 0;
    *lo = pf->words;
    *hi = 0;
    for (uint32_t w = 0; w < pf->words; ++w) {
        sel[w] = md_atomic_load64((const volatile uint64_t*)&pf->added[w]) & md_bitset_word(frames, w * 64);
        if (!sel[w]) continue;
        num += md_popcnt64(sel[w]);
        if (w < *lo) *lo = w;
        *hi = w + 1;
    }
    return num;
}

static inline uint32_t row_count(const md_pair_frames_t* pf, uint64_t row, const uint64_t* sel, uint32_t lo, uint32_t hi) {
    const volatile uint64_t* bits = row_bits(pf, row);
    uint64_t c = 0;
    for (uint32_t w = lo; w < hi; ++w) c += md_popcnt64(md_atomic_load64(bits + w) & sel[w]);
    return (uint32_t)c;
}

size_t md_pair_frames_occupancy(float* out, const md_pair_frames_t* pf, md_bitset_t frames) {
    if (!out || !pf) return 0;
    const size_t total = (size_t)pf->rows * pf->cols;
    uint64_t* sel = md_alloc(pf->alloc, (size_t)pf->words * sizeof(uint64_t));
    if (!sel) return 0;
    uint32_t lo, hi;
    const size_t num = select_frames(sel, &lo, &hi, pf, frames);
    if (num == 0) {
        for (size_t i = 0; i < total; ++i) out[i] = NAN;
        md_free(pf->alloc, sel, (size_t)pf->words * sizeof(uint64_t));
        return 0;
    }
    md_memset(out, 0, total * sizeof(float));
    const uint64_t n = md_atomic_load64(&pf->count);
    const float scale = 1.0f / (float)num;
    for (uint64_t r = 0; r < n; ++r) {
        out[block_of(pf, r)->entry[r % BLOCK_ROWS]] = (float)row_count(pf, r, sel, lo, hi) * scale;
    }
    md_free(pf->alloc, sel, (size_t)pf->words * sizeof(uint64_t));
    return num;
}

size_t md_pair_frames_index(uint32_t** entries, uint32_t** slots, const md_pair_frames_t* pf, struct md_allocator_i* alloc) {
    if (!entries || !slots || !alloc) return 0;
    md_array_shrink(*entries, 0);
    md_array_shrink(*slots, 0);
    if (!pf) return 0;
    const uint64_t n = md_atomic_load64(&pf->count);
    if (n == 0) return 0;
    // As entry << 32 | slot, sorted: the entries ascending, each once
    uint64_t* keys = md_alloc(alloc, (size_t)n * sizeof(uint64_t));
    if (!keys) return 0;
    for (uint64_t r = 0; r < n; ++r) keys[r] = (uint64_t)block_of(pf, r)->entry[r % BLOCK_ROWS] << 32 | r;
    md_sort_radix_u64(keys, (size_t)n);
    md_array_resize(*entries, (size_t)n, alloc);
    md_array_resize(*slots, (size_t)n, alloc);
    for (uint64_t i = 0; i < n; ++i) {
        (*entries)[i] = (uint32_t)(keys[i] >> 32);
        (*slots)[i] = (uint32_t)keys[i];
    }
    md_free(alloc, keys, (size_t)n * sizeof(uint64_t));
    return (size_t)n;
}

size_t md_pair_frames_slot_counts(uint32_t* counts, const md_pair_frames_t* pf, const uint32_t* slots, size_t num, md_bitset_t frames) {
    if (!pf || (num && (!counts || !slots))) return 0;
    md_temp_scope_t temp = md_temp_begin();
    uint64_t* sel = md_temp_alloc_array(temp, uint64_t, pf->words);
    uint32_t lo, hi;
    const size_t nf = select_frames(sel, &lo, &hi, pf, frames);
    const uint64_t n = md_atomic_load64(&pf->count);
    for (size_t k = 0; k < num; ++k) counts[k] = nf && slots[k] < n ? row_count(pf, slots[k], sel, lo, hi) : 0;
    md_temp_end(temp);
    return nf;
}

md_bitset_t md_pair_frames_slot_frames(const md_pair_frames_t* pf, uint32_t slot, struct md_allocator_i* alloc) {
    if (!pf || !alloc || slot >= md_atomic_load64(&pf->count)) return md_bitset_empty_set();
    uint64_t* words = md_alloc(alloc, (size_t)pf->words * sizeof(uint64_t));
    if (!words) return md_bitset_empty_set();
    const volatile uint64_t* bits = row_bits(pf, slot);
    for (uint32_t w = 0; w < pf->words; ++w) words[w] = md_atomic_load64(bits + w);
    const md_bitset_t s = md_bitset_from_words(words, pf->words, 0, alloc);
    md_free(alloc, words, (size_t)pf->words * sizeof(uint64_t));
    return s;
}

md_bitset_t md_pair_frames_entry(const md_pair_frames_t* pf, uint32_t i, uint32_t j, struct md_allocator_i* alloc) {
    if (!pf || !alloc || i >= pf->rows || j >= pf->cols) return md_bitset_empty_set();
    md_pair_frames_t* p = (md_pair_frames_t*)pf;      // the map is read under the lock
    md_mutex_lock(&p->lock);
    uint32_t r;
    const bool found = map_find(p, i * pf->cols + j, &r);
    md_bitset_t s = md_bitset_empty_set();
    if (found) {
        uint64_t* words = md_alloc(alloc, (size_t)pf->words * sizeof(uint64_t));
        if (words) {
            const volatile uint64_t* bits = row_bits(p, r);
            for (uint32_t w = 0; w < pf->words; ++w) words[w] = md_atomic_load64(bits + w);
            s = md_bitset_from_words(words, pf->words, 0, alloc);
            md_free(alloc, words, (size_t)pf->words * sizeof(uint64_t));
        }
    }
    md_mutex_unlock(&p->lock);
    return s;
}
