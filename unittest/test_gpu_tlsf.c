#include "utest.h"

#include <core/md_allocator.h>
#include "../src/core/md_gpu_tlsf.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* CPU-only tests of the sub-allocator behind md_gpu_malloc. No device needed. */

#define GRAN 256u

static uint64_t tlsf_rng(uint64_t* s) {
    *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17;
    return *s;
}

typedef struct { md_tlsf_node_t* node; uint64_t req; } live_t;

static int live_cmp(const void* a, const void* b) {
    const md_tlsf_node_t* x = ((const live_t*)a)->node;
    const md_tlsf_node_t* y = ((const live_t*)b)->node;
    if (x->region != y->region) return (uintptr_t)x->region < (uintptr_t)y->region ? -1 : 1;
    return x->offset < y->offset ? -1 : (x->offset > y->offset ? 1 : 0);
}

UTEST(gpu_tlsf, exact_fit_region_is_found) {
    md_tlsf_t t;
    md_tlsf_init(&t, md_get_heap_allocator(), GRAN);
    int region;
    /* A size just above a class boundary: the rounded-up search skips its
       class, so only the exact-class fallback finds it. */
    const uint64_t size = (1u << 20) + 3 * GRAN;
    ASSERT_TRUE(md_tlsf_add_region(&t, &region, size) != NULL);
    md_tlsf_node_t* n = md_tlsf_alloc(&t, size);
    ASSERT_TRUE(n != NULL);
    EXPECT_EQ(0u, (unsigned)n->offset);
    EXPECT_EQ(size, n->size);
    EXPECT_TRUE(md_tlsf_alloc(&t, GRAN) == NULL);
    md_tlsf_node_t* f = md_tlsf_free(&t, n);
    EXPECT_TRUE(md_tlsf_region_empty(f));
    md_tlsf_remove_region(&t, f);
    EXPECT_EQ(0u, (unsigned)t.free_bytes);
    EXPECT_EQ(0u, (unsigned)t.fl_bitmap);
    md_tlsf_destroy(&t);
}

UTEST(gpu_tlsf, split_and_merge_restore_the_region) {
    md_tlsf_t t;
    md_tlsf_init(&t, md_get_heap_allocator(), GRAN);
    int region;
    ASSERT_TRUE(md_tlsf_add_region(&t, &region, 64 * GRAN) != NULL);

    md_tlsf_node_t* a = md_tlsf_alloc(&t, 1);          /* rounds to GRAN */
    md_tlsf_node_t* b = md_tlsf_alloc(&t, 10 * GRAN);
    md_tlsf_node_t* c = md_tlsf_alloc(&t, 3 * GRAN + 1);
    ASSERT_TRUE(a && b && c);
    EXPECT_EQ(GRAN, (unsigned)a->size);
    EXPECT_EQ(4 * GRAN, (unsigned)c->size);
    EXPECT_EQ(0u, (unsigned)(a->offset % GRAN));
    EXPECT_EQ(0u, (unsigned)(c->offset % GRAN));
    EXPECT_EQ((64 - 15) * GRAN, (unsigned)t.free_bytes);

    /* Free the middle first: it cannot merge, the others then absorb it. */
    EXPECT_FALSE(md_tlsf_region_empty(md_tlsf_free(&t, b)));
    EXPECT_FALSE(md_tlsf_region_empty(md_tlsf_free(&t, a)));
    md_tlsf_node_t* f = md_tlsf_free(&t, c);
    EXPECT_TRUE(md_tlsf_region_empty(f));
    EXPECT_EQ(64 * GRAN, (unsigned)f->size);
    md_tlsf_destroy(&t);
}

UTEST(gpu_tlsf, randomized_against_invariants) {
    md_tlsf_t t;
    md_tlsf_init(&t, md_get_heap_allocator(), GRAN);
    static int regions[3];
    const uint64_t region_size[3] = {1u << 20, 1u << 18, 3u << 19};
    uint64_t total = 0;
    for (int i = 0; i < 3; ++i) {
        ASSERT_TRUE(md_tlsf_add_region(&t, &regions[i], region_size[i]) != NULL);
        total += region_size[i];
    }

    enum { MAX_LIVE = 2048 };
    live_t* live = (live_t*)malloc(MAX_LIVE * sizeof(live_t));
    ASSERT_TRUE(live != NULL);
    int n_live = 0;
    uint64_t rng = 0x9E3779B97F4A7C15ull;
    uint64_t live_bytes = 0;
    int failures = 0;

    for (int iter = 0; iter < 20000; ++iter) {
        const bool do_alloc = n_live == 0 || (n_live < MAX_LIVE && (tlsf_rng(&rng) % 100) < 55);
        if (do_alloc) {
            uint64_t req = 1 + tlsf_rng(&rng) % (tlsf_rng(&rng) % 8 == 0 ? 200000 : 8000);
            md_tlsf_node_t* n = md_tlsf_alloc(&t, req);
            if (!n) { ++failures; continue; }
            EXPECT_FALSE(n->is_free);
            EXPECT_TRUE(n->size >= req);
            EXPECT_TRUE(n->size < ((req + GRAN - 1) & ~(uint64_t)(GRAN - 1)) + GRAN);
            EXPECT_EQ(0u, (unsigned)(n->offset % GRAN));
            live[n_live].node = n;
            live[n_live].req  = req;
            ++n_live;
            live_bytes += n->size;
        } else {
            int i = (int)(tlsf_rng(&rng) % (uint64_t)n_live);
            live_bytes -= live[i].node->size;
            md_tlsf_free(&t, live[i].node);
            live[i] = live[--n_live];
        }
        EXPECT_EQ(total, t.free_bytes + live_bytes);

        if (iter % 997 == 0 && n_live > 1) {
            qsort(live, (size_t)n_live, sizeof(live_t), live_cmp);
            for (int i = 1; i < n_live; ++i) {
                const md_tlsf_node_t* p = live[i - 1].node;
                const md_tlsf_node_t* q = live[i].node;
                if (p->region == q->region) EXPECT_TRUE(p->offset + p->size <= q->offset);
            }
        }
    }
    EXPECT_TRUE(failures < 20000);

    int empty = 0;
    for (int i = 0; i < n_live; ++i) {
        if (md_tlsf_region_empty(md_tlsf_free(&t, live[i].node))) ++empty;
    }
    if (n_live == 0) empty = 3;
    EXPECT_EQ(3, empty);
    EXPECT_EQ(total, t.free_bytes);
    free(live);
    md_tlsf_destroy(&t);
}
