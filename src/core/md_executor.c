#include <core/md_executor.h>
#include <core/md_common.h>
#include <core/md_os.h>

#include <stdio.h>

// Lanes (md_executor.h): a group of up to MD_LANES_MAX lanes, which are the pieces of one md_parallel_for.
//
// Membership: one 64 bit word holds CLOSED, the number of lanes which joined and the number which arrived at the first
// sync. A lane joins (a compare exchange) while CLOSED is not set; the first sync completes when every lane which joined
// has arrived, and the lane which completes it sets CLOSED in the same exchange. So a lane either joins before the group
// closes, and the first sync waits for it, or sees it closed and leaves: nobody waits for a lane which has not started.
//
// Syncs after the first count arrivals against the fixed count. A lane waits by spinning for a while, then blocks on a
// semaphore: it registers as a sleeper, looks at the generation again, and sleeps; the lane which completes a sync bumps
// the generation and then posts one token per registered sleeper. Both sides are sequentially consistent, so a sleeper
// either sees the new generation or is counted. A token nobody took (its sleeper saw the generation in time) wakes a
// later waiter, which looks and sleeps again. Syncs of odd and even generations have a semaphore each: the lane which
// completes a sync may reach the next one and sleep before the sleepers of the last one have woken, and on the same
// semaphore it could take one of their tokens. Two generations on cannot happen while any of them still sleeps.

#ifndef MD_LANES_CHECK
#define MD_LANES_CHECK 1        // every sync records its site, and a group which meets at different ones stops
#endif

#define LANE_SPIN_NS 50000.0    // how long a lane spins at a sync before it sleeps

#if MD_COMPILER_MSVC
#   include <intrin.h>
#endif

// ### Atomics on naturally aligned words shared by the lanes

#if MD_COMPILER_GCC || MD_COMPILER_CLANG
static inline uint32_t a_load32(const volatile uint32_t* p)              { return __atomic_load_n(p, __ATOMIC_ACQUIRE); }
static inline uint32_t a_load32_sc(const volatile uint32_t* p)           { return __atomic_load_n(p, __ATOMIC_SEQ_CST); }
static inline void     a_store32(volatile uint32_t* p, uint32_t v)       { __atomic_store_n(p, v, __ATOMIC_RELEASE); }
static inline uint32_t a_add32(volatile uint32_t* p, uint32_t v)         { return __atomic_fetch_add(p, v, __ATOMIC_SEQ_CST); }
static inline uint32_t a_xchg32(volatile uint32_t* p, uint32_t v)        { return __atomic_exchange_n(p, v, __ATOMIC_SEQ_CST); }
static inline uint64_t a_load64(const volatile uint64_t* p)              { return __atomic_load_n(p, __ATOMIC_SEQ_CST); }
static inline uint64_t a_add64(volatile uint64_t* p, uint64_t v)         { return __atomic_fetch_add(p, v, __ATOMIC_SEQ_CST); }
// The value before: the exchange happened if it is expected
static inline uint64_t a_cas64(volatile uint64_t* p, uint64_t expected, uint64_t desired) {
    __atomic_compare_exchange_n(p, &expected, desired, false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
    return expected;
}
static inline void cpu_relax(void) {
#   if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#   elif defined(__aarch64__) || defined(__arm__)
    __asm__ volatile("yield");
#   endif
}
#elif MD_COMPILER_MSVC
// The interlocked operations are full barriers. Loads and stores: plain on x64 (whose ordering is that of acquire and
// release), with a barrier on ARM64.
#   if defined(_M_ARM64)
#       define LANE_HW_BARRIER() __dmb(_ARM64_BARRIER_ISH)
#   else
#       define LANE_HW_BARRIER() ((void)0)
#   endif
static inline uint32_t a_load32(const volatile uint32_t* p) {
    const uint32_t v = (uint32_t)__iso_volatile_load32((const volatile __int32*)p);
    LANE_HW_BARRIER();
    _ReadWriteBarrier();
    return v;
}
static inline uint32_t a_load32_sc(const volatile uint32_t* p)           { return (uint32_t)_InterlockedOr((volatile long*)p, 0); }
static inline void a_store32(volatile uint32_t* p, uint32_t v) {
    _ReadWriteBarrier();
    LANE_HW_BARRIER();
    __iso_volatile_store32((volatile __int32*)p, (__int32)v);
}
static inline uint32_t a_add32(volatile uint32_t* p, uint32_t v)         { return (uint32_t)_InterlockedExchangeAdd((volatile long*)p, (long)v); }
static inline uint32_t a_xchg32(volatile uint32_t* p, uint32_t v)        { return (uint32_t)_InterlockedExchange((volatile long*)p, (long)v); }
static inline uint64_t a_load64(const volatile uint64_t* p)              { return (uint64_t)_InterlockedCompareExchange64((volatile __int64*)p, 0, 0); }
static inline uint64_t a_add64(volatile uint64_t* p, uint64_t v)         { return (uint64_t)_InterlockedExchangeAdd64((volatile __int64*)p, (__int64)v); }
static inline uint64_t a_cas64(volatile uint64_t* p, uint64_t expected, uint64_t desired) {
    return (uint64_t)_InterlockedCompareExchange64((volatile __int64*)p, (__int64)desired, (__int64)expected);
}
static inline void cpu_relax(void) {
#   if defined(_M_ARM64)
    __yield();
#   else
    _mm_pause();
#   endif
}
#endif

// ### The group

#define LANE_CLOSED      (1ull << 63)
#define LANE_JOINED(s)   ((uint32_t)((s) >> 32) & 0x7FFFFFFFu)
#define LANE_ARRIVED(s)  ((uint32_t)(s))

typedef struct lane_group_t {
    const md_executor_t* ex;
    md_lane_fn  fn;
    void*       arg;
    uint32_t    max_lanes;

    volatile uint64_t open;     // LANE_CLOSED | joined << 32 | arrived at the first sync
    uint32_t    count;          // the lanes of the group, fixed by the first sync
    volatile uint32_t arrived;  // at the sync under way (after the first)
    volatile uint32_t gen;      // syncs completed
    volatile uint32_t sleepers[2]; // lanes registered to sleep, by the parity of the generation they wait on
    volatile uint64_t next;     // the loop under way: its next index, reset at every sync
    volatile uint32_t skipped;  // a claim saw the work cancelled and handed out no more pieces
    uint32_t    cancelled;      // decided at the last sync: the same for every lane until the next one
    uint64_t    slot;           // the value of a broadcast
    md_tick_t   spin_ticks;
    md_semaphore_t sem[2];
#if MD_LANES_CHECK
    int         site[MD_LANES_MAX];
#endif
} lane_group_t;

struct md_lane_t {
    lane_group_t* g;
    uint32_t      idx;
    uint32_t      num_syncs;
};

static bool lane_join(lane_group_t* g, md_lane_t* out_lane) {
    uint64_t s = a_load64(&g->open);
    for (;;) {
        if (s & LANE_CLOSED) return false;
        const uint64_t prev = a_cas64(&g->open, s, s + (1ull << 32));
        if (prev == s) break;
        s = prev;
    }
    out_lane->g = g;
    out_lane->idx = LANE_JOINED(s);
    out_lane->num_syncs = 0;
    MD_ASSERT(out_lane->idx < g->max_lanes);
    return true;
}

static void lane_wait(lane_group_t* g, uint32_t gen) {
    if (a_load32(&g->gen) != gen) return;
    const md_tick_t t0 = md_tick_now();
    for (uint32_t i = 1;; ++i) {
        cpu_relax();
        if (a_load32(&g->gen) != gen) return;
        if ((i & 63) == 0 && md_tick_now() - t0 > g->spin_ticks) break;
    }
    const uint32_t p = gen & 1;
    for (;;) {
        a_add32(&g->sleepers[p], 1);
        if (a_load32_sc(&g->gen) != gen) return;
        md_semaphore_aquire(&g->sem[p]);
        if (a_load32(&g->gen) != gen) return;
    }
}

static void lane_release(lane_group_t* g, uint32_t gen) {
    a_xchg32(&g->gen, gen + 1);
    const uint32_t p = gen & 1;
    const uint32_t n = a_xchg32(&g->sleepers[p], 0);
    if (n) md_semaphore_release_n(&g->sem[p], n);
}

static bool lane_sync(md_lane_t* lane, int site) {
    lane_group_t* g = lane->g;
#if MD_LANES_CHECK
    g->site[lane->idx] = site;
#else
    (void)site;
#endif
    const uint32_t gen = a_load32(&g->gen);
    bool last;
    if (lane->num_syncs++ == 0) {
        // The first: completes when every lane which joined has arrived, and closes the group
        uint64_t s = a_load64(&g->open);
        for (;;) {
            uint64_t n = s + 1;
            last = LANE_ARRIVED(n) == LANE_JOINED(n);
            if (last) n |= LANE_CLOSED;
            const uint64_t prev = a_cas64(&g->open, s, n);
            if (prev == s) break;
            s = prev;
        }
        if (last) g->count = LANE_JOINED(s);
    } else {
        last = a_add32(&g->arrived, 1) + 1 == g->count;
    }

    if (last) {
#if MD_LANES_CHECK
        for (uint32_t i = 1; i < g->count; ++i) {
            if (g->site[i] != g->site[0]) {
                fprintf(stderr, "md_lane_sync: lane 0 met the group at line %d, lane %u at line %d (-1: the end of the wide function)\n",
                        g->site[0], i, g->site[i]);
                md_assert_impl(__FILE__, __LINE__, __func__, "every lane of a group reaches the same syncs");
                break;
            }
        }
#endif
        if (!g->cancelled) {
            const md_executor_t* ex = g->ex;
            if (a_load32(&g->skipped) || (ex && ex->cancelled && ex->cancelled(ex->user))) g->cancelled = 1;
        }
        g->next = 0;
        g->arrived = 0;
        lane_release(g, gen);
    } else {
        lane_wait(g, gen);
    }
    return !g->cancelled;
}

static void lane_run(md_lane_t* lane) {
    lane_group_t* g = lane->g;
    g->fn(lane, g->arg);
#if MD_LANES_CHECK
    // Every lane leaves the function at its end: a lane which skipped a sync the others wait at arrives here instead
    lane_sync(lane, -1);
#endif
}

static void lane_piece(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    lane_group_t* g = (lane_group_t*)arg;
    for (uint32_t i = beg; i < end; ++i) {
        md_lane_t lane;
        if (!lane_join(g, &lane)) return;   // closed, and stays so
        lane_run(&lane);
    }
}

bool md_run_lanes(const md_executor_t* ex, uint32_t max_lanes, md_lane_fn fn, void* arg) {
    MD_ASSERT(fn);
    max_lanes = MD_CLAMP(max_lanes, 1u, (uint32_t)MD_LANES_MAX);
    max_lanes = MD_MIN(max_lanes, md_executor_num_threads(ex));

    lane_group_t g;
    md_memset(&g, 0, sizeof(g));
    g.ex  = ex;
    g.fn  = fn;
    g.arg = arg;
    g.max_lanes = max_lanes;

    if (max_lanes == 1) {
        // One lane, on the calling thread: no sync ever waits
        md_lane_t lane;
        lane_join(&g, &lane);
        lane_run(&lane);
        return !g.cancelled;
    }

    g.spin_ticks = (md_tick_t)(LANE_SPIN_NS / md_tick_to_nanoseconds(1));
    md_semaphore_init(&g.sem[0], 0);
    md_semaphore_init(&g.sem[1], 0);
    // False from the executor says pieces were skipped: lanes which did not join, which the group did without
    md_parallel_for(ex, max_lanes, 1, lane_piece, &g);
    md_semaphore_destroy(&g.sem[0]);
    md_semaphore_destroy(&g.sem[1]);
    // No lane at all: the executor handed out no piece (cancelled before any started)
    return LANE_JOINED(g.open) > 0 && !g.cancelled;
}

uint32_t md_lane_idx(const md_lane_t* lane) {
    return lane->idx;
}

uint32_t md_lane_count(const md_lane_t* lane) {
    MD_ASSERT(lane->num_syncs > 0 && "md_lane_count: known from the group's first sync");
    return lane->g->count;
}

bool md_lane_cancelled(const md_lane_t* lane) {
    return lane->g->cancelled != 0;
}

bool md_lane_sync_at(md_lane_t* lane, int site) {
    return lane_sync(lane, site);
}

uint64_t md_lane_bcast_u64_at(md_lane_t* lane, uint64_t value, uint32_t src, int site) {
    MD_ASSERT((lane->num_syncs > 0 || src == 0) && "md_lane_bcast: from lane 0 before the first sync");
    lane_group_t* g = lane->g;
    if (lane->idx == src) g->slot = value;
    lane_sync(lane, site);
    value = g->slot;
    lane_sync(lane, site);      // nobody writes the slot again before every lane has read it
    return value;
}

bool md_lane_claim_at(md_lane_t* lane, md_lane_range_t* out_range, uint32_t count, uint32_t grain, int site) {
    lane_group_t* g = lane->g;
    if (grain == 0) grain = 1;
    if (!g->cancelled && !a_load32(&g->skipped)) {
        const md_executor_t* ex = g->ex;
        if (ex && ex->cancelled && ex->cancelled(ex->user)) {
            a_store32(&g->skipped, 1);
        } else {
            const uint64_t beg = a_add64(&g->next, grain);
            if (beg < count) {
                out_range->beg = (uint32_t)beg;
                out_range->end = (uint32_t)MD_MIN(beg + grain, (uint64_t)count);
                return true;
            }
        }
    }
    lane_sync(lane, site);
    return false;
}

md_lane_range_t md_lane_split(const md_lane_t* lane, uint32_t count) {
    const uint32_t n = md_lane_count(lane);
    const uint32_t q = count / n, r = count % n, i = lane->idx;
    const uint32_t beg = i * q + MD_MIN(i, r);
    return (md_lane_range_t){ beg, beg + q + (i < r ? 1u : 0u) };
}
