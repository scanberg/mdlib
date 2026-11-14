#pragma once

// What the benchmarks with executables of their own share (md_load_bench, md_infer_bench): a thread pool for the
// executor, the page faults of the process, a hash and the order of timings. Each includes it once.

#include <core/md_executor.h>
#include <core/md_os.h>
#include <core/md_common.h>

#include <stdint.h>
#include <stdbool.h>

#if MD_PLATFORM_WINDOWS
#   ifndef WIN32_LEAN_AND_MEAN
#       define WIN32_LEAN_AND_MEAN
#   endif
#   include <windows.h>
#   include <psapi.h>
#else
#   include <sys/resource.h>
#endif

#define MAX_REPS 64
#define MAX_THREADS 64

// The benchmark's own threads (mdlib owns none): a pool which keeps them, as an application's task system does. Workers
// sleep on a semaphore of their own between loops; a loop wakes them, takes pieces on the calling thread as well, and
// waits for them. Threads made and joined for every loop cost Windows about half a millisecond a loop at 8 threads, and a
// load runs hundreds of loops: that would be measured, and is no load's.
typedef struct bench_pool_t {
    md_executor_t     ex;
    uint32_t          num_workers;          // threads besides the calling one
    md_thread_t*      thread[MAX_THREADS];
    md_semaphore_t    start[MAX_THREADS];
    md_semaphore_t    done;
    volatile int32_t  quit;
    // The loop
    uint32_t          count, grain;
    md_range_fn       body;
    void*             arg;
    volatile uint64_t next;
} bench_pool_t;

typedef struct bench_worker_t {
    bench_pool_t* pool;
    uint32_t      idx;
} bench_worker_t;

static bench_worker_t bench_workers[MAX_THREADS];

static void bench_pool_run(bench_pool_t* pool, uint32_t thread_idx) {
    const uint64_t num_pieces = ((uint64_t)pool->count + pool->grain - 1) / pool->grain;
    for (;;) {
        const uint64_t p = md_atomic_fetch_add64(&pool->next, 1);
        if (p >= num_pieces) break;
        const uint32_t beg = (uint32_t)(p * pool->grain);
        const uint32_t end = (uint32_t)MD_MIN((p + 1) * pool->grain, (uint64_t)pool->count);
        pool->body(pool->arg, beg, end, thread_idx);
    }
}

static void bench_worker_main(void* data) {
    bench_worker_t* w = (bench_worker_t*)data;
    bench_pool_t* pool = w->pool;
    for (;;) {
        md_semaphore_aquire(&pool->start[w->idx]);
        if (pool->quit) break;
        bench_pool_run(pool, w->idx + 1);
        md_semaphore_release(&pool->done);
    }
}

static bool bench_pool_for(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    bench_pool_t* pool = (bench_pool_t*)user;
    if (count == 0) return true;
    pool->count = count;
    pool->grain = grain ? grain : 1;
    pool->body  = body;
    pool->arg   = arg;
    pool->next  = 0;
    // As many workers as there are pieces beyond the calling thread's
    const uint64_t num_pieces = ((uint64_t)count + pool->grain - 1) / pool->grain;
    const uint32_t woken = (uint32_t)MD_MIN((uint64_t)pool->num_workers, num_pieces - 1);
    for (uint32_t i = 0; i < woken; ++i) md_semaphore_release(&pool->start[i]);
    bench_pool_run(pool, 0);
    for (uint32_t i = 0; i < woken; ++i) md_semaphore_aquire(&pool->done);
    return true;
}

static void bench_pool_init(bench_pool_t* pool, uint32_t num_threads) {
    md_memset(pool, 0, sizeof(*pool));
    pool->num_workers = num_threads > 1 ? MD_MIN(num_threads, MAX_THREADS) - 1 : 0;
    pool->ex = (md_executor_t){ .parallel_for = bench_pool_for, .user = pool, .num_threads = pool->num_workers + 1 };
    md_semaphore_init(&pool->done, 0);
    for (uint32_t i = 0; i < pool->num_workers; ++i) {
        md_semaphore_init(&pool->start[i], 0);
        bench_workers[i] = (bench_worker_t){ pool, i };
        pool->thread[i] = md_thread_create(bench_worker_main, &bench_workers[i]);
    }
}

static void bench_pool_free(bench_pool_t* pool) {
    pool->quit = 1;
    for (uint32_t i = 0; i < pool->num_workers; ++i) md_semaphore_release(&pool->start[i]);
    for (uint32_t i = 0; i < pool->num_workers; ++i) {
        md_thread_join(pool->thread[i]);
        md_semaphore_destroy(&pool->start[i]);
    }
    md_semaphore_destroy(&pool->done);
}

// The page faults of the process so far: soft and hard (Windows counts them as one)
static uint64_t page_faults(void) {
#if MD_PLATFORM_WINDOWS
    PROCESS_MEMORY_COUNTERS pmc = { .cb = sizeof(pmc) };
    return GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)) ? (uint64_t)pmc.PageFaultCount : 0;
#else
    struct rusage ru;
    return getrusage(RUSAGE_SELF, &ru) == 0 ? (uint64_t)(ru.ru_minflt + ru.ru_majflt) : 0;
#endif
}

static uint64_t fnv(uint64_t h, const void* p, size_t n) {
    const uint8_t* b = (const uint8_t*)p;
    for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; }
    return h;
}

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}
