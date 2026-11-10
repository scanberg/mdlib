#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

// md_executor_t: how a procedure of mdlib spreads its work over threads it does not own.
//
// mdlib owns no threads and never waits on any. Where a procedure's work is parallel inside it (phases with serial
// merges between them, where a split API would expose its internals), it takes an executor and hands it its loops:
// the caller decides how they run. NULL runs them on the calling thread, and is always correct.
//
// Where the work splits at the API instead, the caller parallelizes around the calls (per frame, per structure), or
// the procedure is split into prepare / run a range / finish, when the caller must decide when pieces run (over
// several frames, progressively). See docs/threading.md.
//
// The contract a procedure may rely on, and an executor must keep:
//
//   - body is called over [0, count) in half open pieces of at least grain (but the last), in any order, on any of
//     up to num_threads threads at once, and parallel_for returns once every call has returned.
//   - thread_idx is in [0, num_threads). No two pieces running at the same time share one, so a procedure can keep
//     scratch or accumulators per thread_idx (allocated for num_threads before the loop, merged after it).
//   - parallel_for returns false if it stopped handing out pieces (the caller's work was cancelled): some of
//     [0, count) was then never called, and the procedure must treat its results as incomplete.
//
// And the one a procedure must keep:
//
//   - A body does not issue parallel work of its own (no nested md_parallel_for) or wait for any. A thread which waits
//     may run other pieces meanwhile, of the same loop too, and they would share the waiting piece's thread_idx.
//     A procedure's phases run one after the other, each one level deep.
//   - The result does not depend on which piece ran on which thread_idx, where it is to be reproducible: outputs are
//     addressed by index, and per thread accumulators are merged in a fixed order or are exact (integers).

typedef void (*md_range_fn)(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx);

typedef struct md_executor_t {
    bool   (*parallel_for)(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg);
    void*    user;
    uint32_t num_threads;   // the bound of thread_idx: at least 1
} md_executor_t;

#ifdef __cplusplus
extern "C" {
#endif

// The number of thread_idx a loop on ex may use: 1 for NULL
static inline uint32_t md_executor_num_threads(const md_executor_t* ex) {
    return (ex && ex->parallel_for && ex->num_threads > 0) ? ex->num_threads : 1;
}

// Runs body over [0, count) on ex, or on the calling thread (thread_idx 0) for NULL. False if ex stopped handing out
// pieces (cancelled).
static inline bool md_parallel_for(const md_executor_t* ex, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    if (count == 0) return true;
    if (!ex || !ex->parallel_for) {
        body(arg, 0, count, 0);
        return true;
    }
    return ex->parallel_for(ex->user, count, grain > 0 ? grain : 1, body, arg);
}

#ifdef __cplusplus
}
#endif
