// md_lanes_bench: the procedures with a wide version (md_executor.h, lanes), wide against fork-join, over threads.
//
//   md_lanes_bench [--reps n] [--threads a,b,c] [--only sync|eigen|cells] [--eigen n,n] [--cells n,n]
//
//   --reps <n>       timed runs per case, thread count and version: the best and the median are reported (default 5)
//   --threads <list> thread counts (default 1, 2, 4, ... and every logical core)
//   --only <what>    one of the three parts below
//   --eigen <list>   matrix sizes of md_linalg_sym_eigen_ex (default 500,1000,2000)
//   --cells <list>   point counts of md_spatial_acc_init (default 100000,1000000,4000000)
//
// Both versions run on the same pool (bench_common.h: workers which sleep on a semaphore each between loops, woken by
// every loop), in turns, so that what differs between them is the version alone:
//
//   fork-join  each phase of the procedure a loop of its own (md_parallel_for), the serial parts between them on the
//              calling thread: a wake of the workers and a wait for them per phase
//   wide       the procedure one function run by every lane (md_run_lanes), the phases its loops, the serial parts on
//              lane 0, syncs between them: one wake of the workers per procedure
//
// sync   what a phase costs with nothing in it: an empty md_parallel_for, against an md_lane_sync in a wide function
// eigen  md_linalg_sym_eigen_ex of a random symmetric matrix. Only its tridiagonalization has a wide version (the steps
//        with more than one piece: n - 384 of them, two syncs each); the divide and conquer and the back transformation
//        are fork-join in both, so the difference between the two is the tridiagonalization's alone.
// cells  md_spatial_acc_init of uniform random points at the density of water (0.1 per A^3), cutoff 6 A, in a periodic
//        box, the cell extent given (as for the frames of a trajectory) so no estimate is made. Its wide version is the
//        sort of the points into cells: a loop of fractional coordinates, then a pass or three of a count, the places
//        (serial, lane 0) and a scatter. The tiers after it are serial in both.
//
// Each line is a thread count: best and median of each version in ms, and fork-join / wide of the medians (above 1:
// wide is faster). Every version and count must produce the same result to the bit: a hash of it, '!' where it
// differs from the first. The header names the benchmark's revision (bumped with every change of it) and the mdlib
// commit it was built from.

#include <core/md_linalg.h>
#include <core/md_spatial_acc.h>
#include <core/md_unitcell.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_executor.h>
#include <core/md_common.h>
#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake
#include "bench_common.h"                 // The thread pool

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// The benchmark's own revision: bumped with every change of it, so that a result names the build it came from
#define MD_LANES_BENCH_REV "R1"

#define MAX_COUNTS 32
#define MAX_SIZES  16

static int reps = 5;
static char note[64];               // printed after a line's result, by a case which has something to say
static int counts[MAX_COUNTS];
static int num_counts = 0;

static double median(double* v, int n) {
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    return (n & 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static double best_of(const double* v, int n) {
    double b = v[0];
    for (int i = 1; i < n; ++i) b = MD_MIN(b, v[i]);
    return b;
}

static void parse_list(const char* s, int* out, int* num, int cap) {
    *num = 0;
    for (const char* c = s; *c && *num < cap; ) {
        const int n = atoi(c);
        if (n > 0) out[(*num)++] = n;
        while (*c && *c != ',') ++c;
        if (*c == ',') ++c;
    }
}

// A run of one version of a case: the time in ms, and a hash of what it produced
typedef double (*case_fn)(const md_executor_t* ex, void* data, uint64_t* out_hash);

static void run_case(const char* title, case_fn fn, void* data) {
    printf("%s\n", title);
    printf("  threads   fork-join best / median     wide best / median   fj/wide   result\n");
    uint64_t first = 0;
    bool have_first = false;
    for (int c = 0; c < num_counts; ++c) {
        static bench_pool_t pool;
        bench_pool_init(&pool, (uint32_t)counts[c]);
        md_executor_t fj = pool.ex;
        fj.fork_join = true;
        md_executor_t wide = pool.ex;
        wide.fork_join = false;
        double t_fj[MAX_REPS], t_wide[MAX_REPS];
        uint64_t h_fj = 0, h_wide = 0;
        bool same = true;
        note[0] = '\0';
        // One untimed run of each first, then the timed ones in turns
        fn(&fj, data, &h_fj);
        fn(&wide, data, &h_wide);
        same &= h_fj == h_wide;
        for (int r = 0; r < reps; ++r) {
            uint64_t h;
            t_fj[r] = fn(&fj, data, &h);
            same &= h == h_fj;
            t_wide[r] = fn(&wide, data, &h);
            same &= h == h_fj;
        }
        if (!have_first) { first = h_fj; have_first = true; }
        same &= h_fj == first;
        const double bf = best_of(t_fj, reps), bw = best_of(t_wide, reps);
        const double mf = median(t_fj, reps), mw = median(t_wide, reps);
        printf("  %7d   %9.3f / %-9.3f   %9.3f / %-9.3f   %6.2f   %016llx%s%s\n", counts[c], bf, mf, bw, mw, mf / mw,
               (unsigned long long)h_fj, same ? "" : " !", note);
        bench_pool_free(&pool);
    }
    printf("\n");
}

// ### sync: an empty phase

#define SYNC_PHASES 10000

static void empty_body(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)arg; (void)beg; (void)end; (void)thread_idx;
}

typedef struct sync_data_t {
    double   us;
    uint32_t lanes;
} sync_data_t;

static void spin_us(double us) {
    const md_tick_t t0 = md_tick_now();
    while (md_tick_to_nanoseconds(md_tick_now() - t0) < us * 1000.0) {}
}

static void empty_wide(md_lane_t* lane, void* arg) {
    sync_data_t* d = (sync_data_t*)arg;
    // Work enough for every worker to wake and join before the first sync closes the group, then the syncs alone
    MD_LANE_FOR(lane, 64, 1, r) { spin_us(100); }
    const md_tick_t t0 = md_tick_now();
    for (int k = 0; k < SYNC_PHASES; ++k) md_lane_sync(lane);
    if (md_lane_idx(lane) == 0) {
        d->us = md_tick_to_milliseconds(md_tick_now() - t0) * 1000.0 / SYNC_PHASES;
        d->lanes = md_lane_count(lane);
    }
}

static double sync_case(const md_executor_t* ex, void* data, uint64_t* out_hash) {
    (void)data;
    *out_hash = 0;
    if (ex->fork_join) {
        const md_tick_t t0 = md_tick_now();
        for (int k = 0; k < SYNC_PHASES; ++k) md_parallel_for(ex, md_executor_num_threads(ex), 1, empty_body, NULL);
        return md_tick_to_milliseconds(md_tick_now() - t0) * 1000.0 / SYNC_PHASES;
    }
    sync_data_t d = {0};
    md_run_lanes(ex, md_executor_num_threads(ex), empty_wide, &d);
    snprintf(note, sizeof(note), "  (%u lanes)", d.lanes);
    return d.us;
}

// ### eigen

typedef struct eigen_data_t {
    size_t  n;
    double* M;
    double* V;
    double* w;
} eigen_data_t;

static double eigen_case(const md_executor_t* ex, void* data, uint64_t* out_hash) {
    eigen_data_t* e = (eigen_data_t*)data;
    memcpy(e->V, e->M, sizeof(double) * e->n * e->n);
    const md_tick_t t0 = md_tick_now();
    const bool ok = md_linalg_sym_eigen_ex(e->V, e->w, e->n, ex);
    const double ms = md_tick_to_milliseconds(md_tick_now() - t0);
    uint64_t h = fnv(1469598103934665603ull, e->w, sizeof(double) * e->n);
    h = fnv(h, e->V, sizeof(double) * e->n * e->n);
    *out_hash = ok ? h : 0;
    return ms;
}

// ### cells

typedef struct cells_data_t {
    size_t        n;
    float*        xyz;
    md_unitcell_t cell;
} cells_data_t;

static double cells_case(const md_executor_t* ex, void* data, uint64_t* out_hash) {
    cells_data_t* c = (cells_data_t*)data;
    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    const md_tick_t t0 = md_tick_now();
    md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .xyz = c->xyz, .count = c->n, .cutoff = 6.0, .cell_ext = 6.0,
                                                         .unitcell = &c->cell, .executor = ex });
    const double ms = md_tick_to_milliseconds(md_tick_now() - t0);
    uint64_t h = fnv(1469598103934665603ull, acc.elem_x, acc.num_elems * sizeof(float));
    h = fnv(h, acc.elem_y, acc.num_elems * sizeof(float));
    h = fnv(h, acc.elem_z, acc.num_elems * sizeof(float));
    h = fnv(h, acc.elem_idx, acc.num_elems * sizeof(uint32_t));
    h = fnv(h, acc.cell_off, (acc.num_cells + 1) * sizeof(uint32_t));
    *out_hash = h;
    md_spatial_acc_free(&acc);
    return ms;
}

int main(int argc, char** argv) {
    const char* only = NULL;
    int eigen_sizes[MAX_SIZES] = { 500, 1000, 2000 };
    int num_eigen = 3;
    int cell_sizes[MAX_SIZES] = { 100000, 1000000, 4000000 };
    int num_cells = 3;

    md_os_sys_info_t si = {0};
    md_os_sys_info_query(&si);
    const int cores = si.num_virtual_cores > 0 ? si.num_virtual_cores : 1;

    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--reps") && i + 1 < argc) {
            reps = atoi(argv[++i]);
        } else if (!strcmp(a, "--threads") && i + 1 < argc) {
            parse_list(argv[++i], counts, &num_counts, MAX_COUNTS);
            for (int c = 0; c < num_counts; ++c) counts[c] = MD_MIN(counts[c], MAX_THREADS);
        } else if (!strcmp(a, "--only") && i + 1 < argc) {
            only = argv[++i];
        } else if (!strcmp(a, "--eigen") && i + 1 < argc) {
            parse_list(argv[++i], eigen_sizes, &num_eigen, MAX_SIZES);
        } else if (!strcmp(a, "--cells") && i + 1 < argc) {
            parse_list(argv[++i], cell_sizes, &num_cells, MAX_SIZES);
        } else {
            printf("usage: %s [--reps n] [--threads a,b,c] [--only sync|eigen|cells] [--eigen n,n] [--cells n,n]\n", argv[0]);
            return a[1] == 'h' ? 0 : 1;
        }
    }
    reps = MD_CLAMP(reps, 1, MAX_REPS);
    if (num_counts == 0) {
        for (int n = 1; n < cores && num_counts < MAX_COUNTS - 1; n *= 2) counts[num_counts++] = n;
        counts[num_counts++] = MD_MIN(cores, MAX_THREADS);
    }

    printf("md_lanes_bench %s: wide against fork-join\n", MD_LANES_BENCH_REV);
    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("CPU %d logical cores | %d timed run%s per version and count\n\n", cores, reps, reps == 1 ? "" : "s");

    if (!only || !strcmp(only, "sync")) {
        run_case("sync: one empty phase, in microseconds (fork-join: md_parallel_for; wide: md_lane_sync)", sync_case, NULL);
    }

    if (!only || !strcmp(only, "eigen")) {
        for (int s = 0; s < num_eigen; ++s) {
            eigen_data_t e = { .n = (size_t)eigen_sizes[s] };
            e.M = malloc(sizeof(double) * e.n * e.n);
            e.V = malloc(sizeof(double) * e.n * e.n);
            e.w = malloc(sizeof(double) * e.n);
            unsigned st = 7u + (unsigned)e.n;
            for (size_t i = 0; i < e.n; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    st = st * 1664525u + 1013904223u;
                    const double v = (double)(st >> 8) / (double)(1u << 24) - 0.5;
                    e.M[i * e.n + j] = v;
                    e.M[j * e.n + i] = v;
                }
                e.M[i * e.n + i] += (double)i * 0.1;
            }
            char title[128];
            snprintf(title, sizeof(title), "eigen: md_linalg_sym_eigen_ex, n = %zu (ms)", e.n);
            run_case(title, eigen_case, &e);
            free(e.M); free(e.V); free(e.w);
        }
    }

    if (!only || !strcmp(only, "cells")) {
        for (int s = 0; s < num_cells; ++s) {
            cells_data_t c = { .n = (size_t)cell_sizes[s] };
            const double ext = cbrt((double)c.n / 0.1);
            c.cell = md_unitcell_from_extent(ext, ext, ext);
            c.xyz = malloc(sizeof(float) * 3 * c.n);
            unsigned st = 11u;
            for (size_t i = 0; i < 3 * c.n; ++i) {
                st = st * 1664525u + 1013904223u;
                c.xyz[i] = (float)((double)(st >> 8) / (double)(1u << 24) * ext);
            }
            char title[128];
            snprintf(title, sizeof(title), "cells: md_spatial_acc_init, %zu points (ms)", c.n);
            run_case(title, cells_case, &c);
            free(c.xyz);
        }
    }
    return 0;
}
