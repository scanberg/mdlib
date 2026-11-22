#include "utest.h"
#include "test_executor.h"

#include <core/md_pair_frames.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_common.h>

#include <math.h>
#include <stdlib.h>

// Against a plain matrix per frame: the frames added in a scrambled order, then each frame, the fractions over frames
// given as a run and as bits, and an entry's frames, as the matrices give them

typedef struct { md_allocator_i iface; int64_t live_bytes; int64_t live_allocs; } pf_count_alloc_t;
static void* pf_count_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file; (void)line;
    pf_count_alloc_t* c = (pf_count_alloc_t*)inst;
    if (ptr) { c->live_bytes -= (int64_t)old_size; c->live_allocs -= 1; }
    if (new_size == 0) { free(ptr); return NULL; }
    c->live_bytes += (int64_t)new_size; c->live_allocs += 1;
    return realloc(ptr, new_size);
}

enum { PF_ROWS = 7, PF_COLS = 5, PF_FRAMES = 1300 };

// A truth of entry e at frame f: some entries never, some always, the rest now and then
static bool pf_truth(uint32_t f, uint32_t e) {
    if (e % 11 == 3) return false;
    if (e == 17) return true;
    uint64_t h = ((uint64_t)f * 2654435761u) ^ ((uint64_t)e * 40503u);
    h ^= h >> 13; h *= 0x5bd1e995; h ^= h >> 15;
    return (h % 7) < (e % 4);
}

static void pf_add_frame(md_pair_frames_t* pf, uint32_t f) {
    uint32_t entries[PF_ROWS * PF_COLS * 2 + 1];
    size_t n = 0;
    for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) if (pf_truth(f, e)) { entries[n++] = e; if (e % 5 == 0) entries[n++] = e; }
    entries[n++] = PF_ROWS * PF_COLS + 4;   // out of range: left out
    md_pair_frames_add(pf, f, entries, n);
}

static void pf_check(int* utest_result, const md_pair_frames_t* pf, md_allocator_i* alloc) {
    // Each frame
    uint8_t m[PF_ROWS * PF_COLS];
    size_t ever = 0;
    for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) {
        bool any = false;
        for (uint32_t f = 0; f < PF_FRAMES && !any; ++f) any = pf_truth(f, e);
        ever += any;
    }
    EXPECT_EQ(ever, md_pair_frames_count(pf));
    for (uint32_t f = 0; f < PF_FRAMES; f += 37) {
        ASSERT_TRUE(md_pair_frames_frame(m, pf, f));
        for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) EXPECT_EQ(pf_truth(f, e) ? 1 : 0, m[e]);
    }
    // Fractions over a run of frames, and over every third frame
    float occ[PF_ROWS * PF_COLS];
    const md_bitset_t run = md_bitset_range(70, 977);
    ASSERT_EQ(907u, md_pair_frames_occupancy(occ, pf, run));
    for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) {
        uint32_t c = 0;
        for (uint32_t f = 70; f < 977; ++f) c += pf_truth(f, e);
        EXPECT_NEAR((float)c / 907.0f, occ[e], 1e-6f);
    }
    uint32_t idx[PF_FRAMES];
    size_t ni = 0;
    for (uint32_t f = 1; f < PF_FRAMES; f += 3) idx[ni++] = f;
    const md_bitset_t thirds = md_bitset_from_indices(idx, ni, alloc);
    ASSERT_EQ(ni, md_pair_frames_occupancy(occ, pf, thirds));
    for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) {
        uint32_t c = 0;
        for (size_t k = 0; k < ni; ++k) c += pf_truth(idx[k], e);
        EXPECT_NEAR((float)c / (float)ni, occ[e], 1e-6f);
    }
    // Past the frames: none, and NaN
    EXPECT_EQ(0u, md_pair_frames_occupancy(occ, pf, md_bitset_range(PF_FRAMES, PF_FRAMES + 100)));
    EXPECT_TRUE(isnan(occ[0]));
    // The same as counts of the entries true in some: over the run, and at a frame
    uint32_t* entries = NULL;
    uint32_t* counts = NULL;
    ASSERT_EQ(907u, md_pair_frames_counts(&entries, &counts, pf, run, alloc));
    size_t k = 0;
    for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) {
        uint32_t c = 0;
        for (uint32_t f = 70; f < 977; ++f) c += pf_truth(f, e);
        if (!c) continue;
        ASSERT_LT(k, md_array_size(entries));
        EXPECT_EQ(e, entries[k]);
        EXPECT_EQ(c, counts[k]);
        ++k;
    }
    EXPECT_EQ(k, md_array_size(entries));
    EXPECT_EQ(k, md_array_size(counts));
    for (uint32_t f = 0; f < PF_FRAMES; f += 101) {
        ASSERT_EQ(1u, md_pair_frames_counts(&entries, &counts, pf, md_bitset_range(f, f + 1), alloc));
        k = 0;
        for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) {
            if (!pf_truth(f, e)) continue;
            ASSERT_LT(k, md_array_size(entries));
            EXPECT_EQ(e, entries[k]);
            EXPECT_EQ(1u, counts[k]);
            ++k;
        }
        EXPECT_EQ(k, md_array_size(entries));
    }
    EXPECT_EQ(0u, md_pair_frames_counts(&entries, &counts, pf, md_bitset_range(PF_FRAMES, PF_FRAMES + 100), alloc));
    EXPECT_EQ(0u, md_array_size(entries));
    EXPECT_EQ(0u, md_array_size(counts));
    // An entry's frames
    for (uint32_t e = 0; e < PF_ROWS * PF_COLS; ++e) {
        const md_bitset_t s = md_pair_frames_entry(pf, e / PF_COLS, e % PF_COLS, alloc);
        uint32_t c = 0;
        for (uint32_t f = 0; f < PF_FRAMES; ++f) {
            c += pf_truth(f, e);
            if (f % 97 == 0) EXPECT_EQ(pf_truth(f, e), md_bitset_test(s, f));
        }
        EXPECT_EQ((uint64_t)c, md_bitset_count(s));
    }
    const md_bitset_t added = md_pair_frames_added(pf, alloc);
    EXPECT_EQ((uint64_t)PF_FRAMES, md_bitset_count(added));
}

UTEST(pair_frames, against_matrices) {
    pf_count_alloc_t ca = {.iface = {.inst = NULL, .realloc = pf_count_realloc}};
    ca.iface.inst = (struct md_allocator_o*)&ca;
    md_allocator_i* alloc = &ca.iface;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));

    EXPECT_TRUE(md_pair_frames_create(0, 5, 10, alloc) == NULL);
    EXPECT_TRUE(md_pair_frames_create(7, 5, 0, alloc) == NULL);
    EXPECT_TRUE(md_pair_frames_create(1u << 16, 1u << 16, 10, alloc) == NULL);

    md_pair_frames_t* pf = md_pair_frames_create(PF_ROWS, PF_COLS, PF_FRAMES, alloc);
    ASSERT_TRUE(pf != NULL);
    uint8_t m[PF_ROWS * PF_COLS];
    EXPECT_FALSE(md_pair_frames_frame(m, pf, 0));
    // A stride coprime with the frames: every frame once, scrambled
    for (uint32_t k = 0; k < PF_FRAMES; ++k) pf_add_frame(pf, (k * 389u) % PF_FRAMES);
    // Once each; past the end, nothing
    const uint32_t one = 0;
    EXPECT_FALSE(md_pair_frames_add(pf, 5, &one, 1));
    EXPECT_FALSE(md_pair_frames_add(pf, PF_FRAMES, &one, 1));
    pf_check(utest_result, pf, arena);

    // A copy holds the same
    md_pair_frames_t* cp = md_pair_frames_copy(pf, alloc);
    ASSERT_TRUE(cp != NULL);
    pf_check(utest_result, cp, arena);
    md_pair_frames_free(cp);
    md_pair_frames_free(pf);
    EXPECT_EQ(0, ca.live_bytes);
    EXPECT_EQ(0, ca.live_allocs);
    md_arena_allocator_destroy(arena);
}

// Many rows, past a page of blocks and the map's first sizes: a matrix of 3000 × 3000 with entries true here and there
UTEST(pair_frames, many_rows) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_pair_frames_t* pf = md_pair_frames_create(3000, 3000, 100, heap);
    ASSERT_TRUE(pf != NULL);
    enum { N = 20000 };
    uint32_t* entries = malloc(N * sizeof(uint32_t));
    for (uint32_t f = 0; f < 100; ++f) {
        // Distinct within a frame: k, or past every frame's k, scattered by a multiplier coprime with 9,000,000
        for (uint32_t k = 0; k < N; ++k) entries[k] = (uint32_t)(((uint64_t)(k % 3 ? k + (f + 1) * N : k) * 7919u) % 9000000u);
        ASSERT_TRUE(md_pair_frames_add(pf, f, entries, N));
    }
    float* occ = malloc(9000000u * sizeof(float));
    ASSERT_EQ(100u, md_pair_frames_occupancy(occ, pf, md_bitset_range(0, 100)));
    // k ≡ 0 mod 3: the same entry every frame
    for (uint32_t k = 0; k < N; k += 3 * 101) EXPECT_EQ(1.0f, occ[((uint64_t)k * 7919u) % 9000000u]);
    size_t ever = 0;
    double sum = 0;
    for (uint32_t i = 0; i < 9000000u; ++i) { ever += occ[i] > 0.0f; sum += occ[i]; }
    EXPECT_EQ(ever, md_pair_frames_count(pf));
    EXPECT_NEAR((double)N, sum, 1e-2);
    // As counts: every entry ever true, ascending, as the fractions have them
    uint32_t* entries_of = NULL;
    uint32_t* counts = NULL;
    ASSERT_EQ(100u, md_pair_frames_counts(&entries_of, &counts, pf, md_bitset_range(0, 100), heap));
    ASSERT_EQ(ever, md_array_size(entries_of));
    for (size_t i = 0; i < ever; ++i) {
        if (i) EXPECT_LT(entries_of[i - 1], entries_of[i]);
        EXPECT_NEAR(occ[entries_of[i]], (float)counts[i] / 100.0f, 1e-6f);
    }
    md_array_free(entries_of, heap);
    md_array_free(counts, heap);
    free(occ);
    free(entries);
    md_pair_frames_free(pf);
}

// Frames from several threads at once, while another reads: the same as from one
typedef struct { md_pair_frames_t* pf; } pf_job_t;
static void pf_job(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    pf_job_t* j = arg;
    for (uint32_t f = beg; f < end; ++f) {
        pf_add_frame(j->pf, f);
        float occ[PF_ROWS * PF_COLS];
        md_pair_frames_occupancy(occ, j->pf, md_bitset_range(0, PF_FRAMES));
        for (int e = 0; e < PF_ROWS * PF_COLS; ++e) if (!(occ[e] >= 0.0f && occ[e] <= 1.0f)) md_pair_frames_add(j->pf, 0, NULL, 0);
    }
}

UTEST(pair_frames, threads) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_allocator_i* arena = md_arena_allocator_create(heap, MD_MEGABYTES(1));
    md_pair_frames_t* pf = md_pair_frames_create(PF_ROWS, PF_COLS, PF_FRAMES, heap);
    ASSERT_TRUE(pf != NULL);
    pf_job_t job = {pf};
    ASSERT_TRUE(md_parallel_for(test_executor_of_threads(8), PF_FRAMES, 7, pf_job, &job));
    pf_check(utest_result, pf, arena);
    md_pair_frames_free(pf);
    md_arena_allocator_destroy(arena);
}
