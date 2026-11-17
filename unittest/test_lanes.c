#include "utest.h"
#include "test_executor.h"

#include <string.h>

// Lanes (md_executor.h): md_run_lanes, its syncs, broadcasts and loops, on the tests' executors and on none

// Every index of a few loops once, a broadcast, and a sum of doubles whose rounding depends on the order of its terms,
// addressed by piece: the same whatever the lanes
typedef struct wide_t {
    uint32_t          n, grain;
    const double*     x;
    double*           part;         // [pieces]
    volatile uint64_t hits[4096];
    volatile uint64_t lanes;        // lanes which ran
    volatile uint64_t bad;          // a lane which saw something it should not have
    uint64_t          token;        // lane 0's broadcast
    uint32_t          max_lanes;
    double            total[MD_LANES_MAX];
} wide_t;

static void wide_fn(md_lane_t* lane, void* arg) {
    wide_t* w = (wide_t*)arg;
    md_atomic_fetch_add64(&w->lanes, 1);
    if (md_lane_idx(lane) >= w->max_lanes) md_atomic_fetch_add64(&w->bad, 1);

    MD_LANE_FOR(lane, w->n, w->grain, r) {
        double acc = 0;
        for (uint32_t i = r.beg; i < r.end; ++i) {
            acc += w->x[i];
            md_atomic_fetch_add64(&w->hits[i], 1);
        }
        w->part[r.beg / w->grain] = acc;
    }
    const uint32_t count = md_lane_count(lane);
    if (count == 0 || count > w->max_lanes) md_atomic_fetch_add64(&w->bad, 1);

    // Serial, on lane 0, then everyone sees it
    uint64_t t = 0;
    if (md_lane_idx(lane) == 0) t = 0x1234567890ABCDEFull + count;
    t = md_lane_bcast_u64(lane, t, 0);
    if (t != 0x1234567890ABCDEFull + count) md_atomic_fetch_add64(&w->bad, 1);

    // Duplicated on every lane: the sum of the pieces in order
    const uint32_t pieces = (w->n + w->grain - 1) / w->grain;
    double total = 0;
    for (uint32_t p = 0; p < pieces; ++p) total += w->part[p];
    w->total[md_lane_idx(lane)] = total;

    // A second loop over the same range, and the even split
    MD_LANE_FOR(lane, w->n, 1, r) {
        for (uint32_t i = r.beg; i < r.end; ++i) md_atomic_fetch_add64(&w->hits[i], 1);
    }
    const md_lane_range_t s = md_lane_split(lane, w->n);
    for (uint32_t i = s.beg; i < s.end; ++i) md_atomic_fetch_add64(&w->hits[i], 1);
    md_lane_sync(lane);
}

static void check_wide(int* utest_result, const md_executor_t* ex, uint32_t max_lanes, uint32_t n, uint32_t grain) {
    static double x[4096], part[4096];
    for (uint32_t i = 0; i < n; ++i) x[i] = 1.0 / (1.0 + i);
    double ref = 0;
    for (uint32_t p = 0; p * grain < n; ++p) {
        double acc = 0;
        for (uint32_t i = p * grain; i < MD_MIN(n, (p + 1) * grain); ++i) acc += x[i];
        ref += acc;
    }
    wide_t* w = (wide_t*)calloc(1, sizeof(wide_t));
    w->n = n; w->grain = grain; w->x = x; w->part = part;
    w->max_lanes = MD_MIN(max_lanes, md_executor_num_threads(ex));
    EXPECT_TRUE(md_run_lanes(ex, max_lanes, wide_fn, w));
    EXPECT_GE(w->lanes, 1u);
    EXPECT_LE(w->lanes, (uint64_t)w->max_lanes);
    EXPECT_EQ(0u, (uint32_t)w->bad);
    uint32_t wrong = 0;
    for (uint32_t i = 0; i < n; ++i) wrong += w->hits[i] != 3;
    EXPECT_EQ(0u, wrong);
    for (uint32_t l = 0; l < w->lanes; ++l) EXPECT_EQ(ref, w->total[l]);
    free(w);
}

UTEST(lanes, null_runs_one_lane) {
    check_wide(utest_result, NULL, 8, 1000, 7);
}

UTEST(lanes, threads) {
    for (uint32_t t = 1; t <= 8; ++t) {
        test_executor_t te = test_executor_threads(t);
        const md_executor_t* ex = test_executor(&te);
        check_wide(utest_result, ex, t, 4096, 64);
        check_wide(utest_result, ex, t, 3, 8);
        check_wide(utest_result, ex, 2, 1000, 1);
    }
}

UTEST(lanes, shuffled_is_one_lane) {
    test_executor_t te = test_executor_shuffled(4);
    const md_executor_t* ex = test_executor(&te);
    check_wide(utest_result, ex, 4, 1000, 10);
}

// Many syncs in a row, each a phase whose writes every lane checks after it: lanes which outrun a sync, or a sync
// which lets a lane through early, show as a wrong value. Some phases are long on one lane, so others sleep.
typedef struct phases_t {
    volatile uint64_t value[MD_LANES_MAX];
    volatile uint64_t bad;
    uint32_t          num_phases;
} phases_t;

static void phases_fn(md_lane_t* lane, void* arg) {
    phases_t* p = (phases_t*)arg;
    MD_LANE_FOR(lane, 1, 1, r) { (void)r; }     // the first sync: the group is closed from here
    const uint32_t me = md_lane_idx(lane), n = md_lane_count(lane);
    for (uint32_t k = 1; k <= p->num_phases; ++k) {
        p->value[me] = k;
        if (k % 97 == 0 && me == k % n) md_thread_sleep(1);
        md_lane_sync(lane);
        for (uint32_t i = 0; i < n; ++i) {
            if (p->value[i] != k) md_atomic_fetch_add64(&p->bad, 1);
        }
        md_lane_sync(lane);
    }
}

UTEST(lanes, syncs_hold) {
    for (uint32_t t = 2; t <= 8; t *= 2) {
        test_executor_t te = test_executor_threads(t);
        phases_t p = {.num_phases = 2000};
        EXPECT_TRUE(md_run_lanes(test_executor(&te), t, phases_fn, &p));
        EXPECT_EQ(0u, (uint32_t)p.bad);
    }
}

// Syncs at which lanes sleep, back to back: the lane which completes one reaches the next and sleeps there while the
// sleepers of the last one have not yet woken. A wakeup taken by the wrong sleeper hangs this.
static void sleepers_fn(md_lane_t* lane, void* arg) {
    phases_t* p = (phases_t*)arg;
    MD_LANE_FOR(lane, 1, 1, r) { (void)r; }
    const uint32_t me = md_lane_idx(lane), n = md_lane_count(lane);
    for (uint32_t k = 1; k <= p->num_phases; ++k) {
        p->value[me] = k;
        if (me == k % n) {
            // Busy past the others' spin, so that they sleep at this sync; then straight on to the next one
            const md_tick_t t0 = md_tick_now();
            while (md_tick_to_nanoseconds(md_tick_now() - t0) < 300000.0) {}
        }
        md_lane_sync(lane);
        for (uint32_t i = 0; i < n; ++i) {
            if (p->value[i] != k) md_atomic_fetch_add64(&p->bad, 1);
        }
    }
}

UTEST(lanes, sleepers_wake) {
    for (uint32_t t = 2; t <= 8; t *= 2) {
        test_executor_t te = test_executor_threads(t);
        phases_t p = {.num_phases = 1000};
        EXPECT_TRUE(md_run_lanes(test_executor(&te), t, sleepers_fn, &p));
        EXPECT_EQ(0u, (uint32_t)p.bad);
    }
    // More lanes than cores, syncing as fast as they can: lanes are descheduled at syncs and sleep there
    test_executor_t te = test_executor_threads(16);
    phases_t p = {.num_phases = 20000};
    EXPECT_TRUE(md_run_lanes(test_executor(&te), 16, phases_fn, &p));
    EXPECT_EQ(0u, (uint32_t)p.bad);
}

// Cancelled part way: every lane sees it at the same sync, and md_run_lanes says so
typedef struct cancel_t {
    test_executor_t*  te;
    volatile uint64_t said[2];      // syncs which returned true / false
    volatile uint64_t after;        // lanes which went on past the sync which said so
} cancel_t;

static void cancel_fn(md_lane_t* lane, void* arg) {
    cancel_t* c = (cancel_t*)arg;
    MD_LANE_FOR(lane, 100, 1, r) { (void)r; }
    for (int k = 0; k < 10; ++k) {
        if (k == 3 && md_lane_idx(lane) == 0) c->te->cancel = 1;
        const bool ok = md_lane_sync(lane);
        md_atomic_fetch_add64(&c->said[ok ? 0 : 1], 1);
        if (!ok) return;
    }
    md_atomic_fetch_add64(&c->after, 1);
}

UTEST(lanes, cancelled_at_a_sync) {
    test_executor_t te = test_executor_threads(4);
    cancel_t c = {.te = &te};
    EXPECT_FALSE(md_run_lanes(test_executor(&te), 4, cancel_fn, &c));
    EXPECT_EQ(0u, (uint32_t)c.after);
    EXPECT_GT(c.said[1], 0u);
    // Each lane said true at the same syncs (k = 0 .. 2), then false once
    EXPECT_EQ(c.said[0] % 3, 0u);
    EXPECT_EQ(c.said[0] / 3, c.said[1]);
}

UTEST(lanes, cancelled_before_it_starts) {
    test_executor_t te = test_executor_threads(4);
    te.cancel = 1;
    cancel_t c = {.te = &te};
    EXPECT_FALSE(md_run_lanes(test_executor(&te), 4, cancel_fn, &c));
}

