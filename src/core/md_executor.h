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
    // Optional: whether the work the executor runs for was asked to stop. md_run_lanes asks at its syncs and claims,
    // on any of its lanes' threads, so an executor which keeps this per thread carries it into its pieces.
    bool   (*cancelled)(void* user);
    // EXPERIMENT (lanes): a procedure with a wide version and a fork-join version runs the fork-join one
    bool     fork_join;
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

// ### Lanes: wide sections (EXPERIMENT)
//
// md_run_lanes(ex, max_lanes, fn, arg) runs fn(lane, arg) on up to max_lanes threads at once: one function, run by
// every lane of a group, with its locals on each lane's own stack. Work is split inside it (MD_LANE_FOR), serial parts
// are narrowed to a lane (md_lane_idx(lane) == 0), and the lanes meet at syncs between them. No phase is a function
// of its own, and nothing is handed to another thread: the control flow of the whole procedure reads top to bottom.
//
// Membership is elastic. The lanes are the pieces of one md_parallel_for on ex, and a piece joins the group when it
// starts, until the group's first sync closes it; a piece which starts after that returns at once. So a lane only
// ever waits for lanes which have started, and an executor which runs pieces one after the other (or none at all
// but on the calling thread) still runs the group, with one lane. Hence:
//
//   - md_lane_idx is known from the start; md_lane_count only after the first sync, as is md_lane_split.
//   - Lanes which join late join only if the first sync has not been reached: open a wide function with a loop
//     (MD_LANE_FOR), and do its setup before going wide, not on lane 0 behind a broadcast at the top.
//   - Every lane reaches the same syncs in the same order: a sync on one side of a branch which not every lane takes
//     waits forever. Every sync records where it was called from, and a group which meets at different places stops
//     with an assertion instead. A body of MD_LANE_FOR is left by its end only (no break or return): the loop's last
//     claim is its closing sync.
//   - A lane does not wait for anything but its group: no md_parallel_for, and no other wait which could run pieces.
//   - Results are addressed by piece, not by lane, where they are to be reproducible: how many lanes a group has
//     depends on how soon the executor's threads started.
//
// Cancellation is decided at syncs, for every lane at once: once the executor says its work was cancelled (or a claim
// saw it, and handed out no more pieces), the next sync returns false on every lane, and md_lane_cancelled is true
// from then on.
//
// What md_run_lanes asks of an executor beyond md_parallel_for's contract: a piece which has started runs on its
// thread until it returns, and the thread runs nothing else meanwhile (a lane waits at syncs of its own, never in the
// executor). Every executor which runs a piece as a plain call keeps it.
//
// EXPERIMENT: md_executor_t.fork_join selects the fork-join version of a procedure which has both (the spatial
// acceleration structure's cells, the tridiagonalization of md_linalg_sym_eigen_ex), to compare them on one executor.

typedef struct md_lane_t md_lane_t;
typedef struct md_lane_range_t { uint32_t beg, end; } md_lane_range_t;
typedef void (*md_lane_fn)(md_lane_t* lane, void* arg);

#define MD_LANES_MAX 256

// Runs fn on a group of up to max_lanes lanes of ex (one, on the calling thread, for NULL). False if cancelled.
bool md_run_lanes(const md_executor_t* ex, uint32_t max_lanes, md_lane_fn fn, void* arg);

uint32_t md_lane_idx(const md_lane_t* lane);       // [0, md_lane_count)
uint32_t md_lane_count(const md_lane_t* lane);     // after the first sync
bool     md_lane_cancelled(const md_lane_t* lane); // as decided at the last sync

// A barrier for the group. False when the work was cancelled (on every lane alike).
bool     md_lane_sync_at(md_lane_t* lane, int site);
// value of lane src, on every lane (two syncs). src 0 before the first sync, the only lane known to be in the group.
uint64_t md_lane_bcast_u64_at(md_lane_t* lane, uint64_t value, uint32_t src, int site);
// The next piece of a loop over [0, count) claimed by the lane: false, after the loop's closing sync, when none is left
bool     md_lane_claim_at(md_lane_t* lane, md_lane_range_t* out_range, uint32_t count, uint32_t grain, int site);
// The lane's share of [0, count), split evenly by md_lane_count: for uniform work only, after the first sync
md_lane_range_t md_lane_split(const md_lane_t* lane, uint32_t count);

#define md_lane_sync(lane)                  md_lane_sync_at((lane), __LINE__)
#define md_lane_bcast_u64(lane, value, src) md_lane_bcast_u64_at((lane), (value), (src), __LINE__)
#define md_lane_bcast_ptr(lane, ptr, src)   ((void*)(uintptr_t)md_lane_bcast_u64_at((lane), (uint64_t)(uintptr_t)(ptr), (src), __LINE__))

// for each piece r ([r.beg, r.end), at most grain long) of [0, count) claimed by this lane, then a sync
#define MD_LANE_FOR(lane, count, grain, r) \
    for (md_lane_range_t r = {0, 0}; md_lane_claim_at((lane), &r, (count), (grain), __LINE__); )

#ifdef __cplusplus
}
#endif
