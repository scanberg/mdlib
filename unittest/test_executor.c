#include "utest.h"
#include "test_executor.h"

#include <string.h>

// md_parallel_for and the executor contract (md_executor.h), as the tests' own executors keep it

typedef struct visit_t {
    volatile uint64_t* hits;      // per index
    volatile uint64_t* by_thread; // per thread_idx
    uint32_t           num_threads;
    volatile uint64_t  bad_idx;
} visit_t;

static void visit_body(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    visit_t* v = (visit_t*)arg;
    if (thread_idx >= v->num_threads) { md_atomic_fetch_add64(&v->bad_idx, 1); return; }
    for (uint32_t i = beg; i < end; ++i) md_atomic_fetch_add64(&v->hits[i], 1);
    md_atomic_fetch_add64(&v->by_thread[thread_idx], end - beg);
}

static void check_visits(int* utest_result, const md_executor_t* ex, uint32_t count, uint32_t grain) {
    uint64_t hits[1000] = {0};
    uint64_t by_thread[TEST_EXECUTOR_MAX_THREADS] = {0};
    ASSERT_LE(count, 1000u);
    visit_t v = {.hits = hits, .by_thread = by_thread, .num_threads = md_executor_num_threads(ex)};
    EXPECT_TRUE(md_parallel_for(ex, count, grain, visit_body, &v));
    EXPECT_EQ(0u, (uint32_t)v.bad_idx);
    uint32_t wrong = 0;
    for (uint32_t i = 0; i < count; ++i) wrong += hits[i] != 1;
    EXPECT_EQ(0u, wrong);
    uint64_t sum = 0;
    for (uint32_t t = 0; t < v.num_threads; ++t) sum += by_thread[t];
    EXPECT_EQ((uint64_t)count, sum);
}

UTEST(executor, null_runs_on_the_calling_thread) {
    EXPECT_EQ(1u, md_executor_num_threads(NULL));
    check_visits(utest_result, NULL, 1000, 7);
    // nothing to do: done
    EXPECT_TRUE(md_parallel_for(NULL, 0, 1, visit_body, NULL));
}

UTEST(executor, every_index_once_on_threads) {
    test_executor_t te = test_executor_threads(4);
    const md_executor_t* ex = test_executor(&te);
    EXPECT_EQ(4u, md_executor_num_threads(ex));
    check_visits(utest_result, ex, 1000, 1);
    check_visits(utest_result, ex, 1000, 64);
    check_visits(utest_result, ex, 3, 8);
}

UTEST(executor, every_index_once_shuffled) {
    test_executor_t te = test_executor_shuffled(3);
    const md_executor_t* ex = test_executor(&te);
    check_visits(utest_result, ex, 1000, 1);
    check_visits(utest_result, ex, 997, 10);
    check_visits(utest_result, ex, 1, 1);
}

UTEST(executor, a_cancelled_loop_says_so) {
    test_executor_t te = test_executor_threads(4);
    const md_executor_t* ex = test_executor(&te);
    te.cancel = 1;
    uint64_t hits[100] = {0};
    uint64_t by_thread[TEST_EXECUTOR_MAX_THREADS] = {0};
    visit_t v = {.hits = hits, .by_thread = by_thread, .num_threads = 4};
    EXPECT_FALSE(md_parallel_for(ex, 100, 1, visit_body, &v));
    uint64_t sum = 0;
    for (int i = 0; i < 100; ++i) sum += hits[i];
    EXPECT_EQ(0u, (uint32_t)sum);
}
