#pragma once

// Executors for the tests (md_executor.h): mdlib owns no threads, its tests do.
//
//   test_executor_threads    num_threads real threads, the calling one of them, taking pieces in turn
//   test_executor_shuffled   one thread, the pieces in a scrambled order and handed thread_idx in turn: a procedure
//                            whose result depends on the order of its pieces, or on which piece had which thread_idx,
//                            shows it here deterministically

#include <core/md_executor.h>
#include <core/md_os.h>
#include <core/md_common.h>

#include <stdint.h>
#include <stdbool.h>

#define TEST_EXECUTOR_MAX_THREADS 64

typedef struct test_executor_t {
    md_executor_t      ex;
    volatile int32_t   cancel;      // set: no more pieces are handed out, and parallel_for returns false
    volatile uint64_t  pieces;      // handed out so far, over every loop
} test_executor_t;

typedef struct test_executor_loop_t {
    test_executor_t*  te;
    uint32_t          count, grain;
    md_range_fn       body;
    void*             arg;
    volatile uint64_t next;         // the next piece
    volatile uint64_t skipped;
} test_executor_loop_t;

typedef struct test_executor_thread_t {
    test_executor_loop_t* loop;
    uint32_t              idx;
} test_executor_thread_t;

static void test_executor_thread_run(void* data) {
    test_executor_thread_t* t = (test_executor_thread_t*)data;
    test_executor_loop_t* l = t->loop;
    const uint64_t num_pieces = (l->count + l->grain - 1) / l->grain;
    for (;;) {
        const uint64_t p = md_atomic_fetch_add64(&l->next, 1);
        if (p >= num_pieces) break;
        if (l->te->cancel) { md_atomic_fetch_add64(&l->skipped, 1); continue; }
        md_atomic_fetch_add64(&l->te->pieces, 1);
        const uint32_t beg = (uint32_t)(p * l->grain);
        const uint32_t end = (uint32_t)((p + 1) * l->grain < l->count ? (p + 1) * l->grain : l->count);
        l->body(l->arg, beg, end, t->idx);
    }
}

static bool test_executor_threads_for(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    test_executor_t* te = (test_executor_t*)user;
    test_executor_loop_t loop = {.te = te, .count = count, .grain = grain, .body = body, .arg = arg};
    test_executor_thread_t th[TEST_EXECUTOR_MAX_THREADS];
    md_thread_t* handles[TEST_EXECUTOR_MAX_THREADS] = {0};
    const uint32_t n = te->ex.num_threads;
    for (uint32_t i = 0; i < n; ++i) { th[i].loop = &loop; th[i].idx = i; }
    for (uint32_t i = 1; i < n; ++i) handles[i] = md_thread_create(test_executor_thread_run, &th[i]);
    test_executor_thread_run(&th[0]);
    for (uint32_t i = 1; i < n; ++i) {
        if (handles[i]) md_thread_join(handles[i]);
    }
    return md_atomic_load64(&loop.skipped) == 0;
}

static bool test_executor_shuffled_for(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    test_executor_t* te = (test_executor_t*)user;
    const uint64_t num_pieces = ((uint64_t)count + grain - 1) / grain;
    // Every piece once, in the order of a stride coprime with the number of pieces
    uint64_t stride = num_pieces / 2 + 1;
    while (num_pieces > 1) {
        uint64_t a = stride, b = num_pieces;
        while (b) { const uint64_t r = a % b; a = b; b = r; }
        if (a == 1) break;
        ++stride;
    }
    bool all = true;
    for (uint64_t k = 0; k < num_pieces; ++k) {
        if (te->cancel) { all = false; break; }
        const uint64_t p = (k * stride + 3) % num_pieces;
        const uint32_t beg = (uint32_t)(p * grain);
        const uint32_t end = (uint32_t)((p + 1) * grain < count ? (p + 1) * grain : count);
        md_atomic_fetch_add64(&te->pieces, 1);
        body(arg, beg, end, (uint32_t)(k % te->ex.num_threads));
    }
    return all;
}

static inline test_executor_t test_executor_threads(uint32_t num_threads) {
    test_executor_t te = {0};
    te.ex.parallel_for = test_executor_threads_for;
    te.ex.num_threads = num_threads < 1 ? 1 : (num_threads > TEST_EXECUTOR_MAX_THREADS ? TEST_EXECUTOR_MAX_THREADS : num_threads);
    return te;
}

static inline test_executor_t test_executor_shuffled(uint32_t num_threads) {
    test_executor_t te = {0};
    te.ex.parallel_for = test_executor_shuffled_for;
    te.ex.num_threads = num_threads < 1 ? 1 : num_threads;
    return te;
}

// The executor of te, with its user set: te must stay where it is while ex is used
static inline const md_executor_t* test_executor(test_executor_t* te) {
    te->ex.user = te;
    return &te->ex;
}

// An executor of n real threads for as long as the test runs: one per n in each translation unit
static inline const md_executor_t* test_executor_of_threads(uint32_t n) {
    static test_executor_t pool[TEST_EXECUTOR_MAX_THREADS + 1];
    n = n < 1 ? 1 : (n > TEST_EXECUTOR_MAX_THREADS ? TEST_EXECUTOR_MAX_THREADS : n);
    if (!pool[n].ex.parallel_for) pool[n] = test_executor_threads(n);
    return test_executor(&pool[n]);
}
