// md_infer_bench: how the topology of a system is inferred, and over threads.
//
//   md_infer_bench [--reps n] [--threads a,b,c] [--reuse-memory] [--bonds] [--stages] [file ...]
//
//   --reps <n>       timed runs per thread count, the best and the median are reported (default 5)
//   --threads <list> thread counts to run with (default 1, 2, 4, ... and every logical core)
//   --reuse-memory   every run into the same memory, reset between runs, so that only the first has pages to fault in:
//                    against a run without it, what the page faults of a run cost
//   --bonds          md_bond_infer alone, not md_system_infer
//   --stages         the time of each step of md_system_infer, the median over the timed runs, per thread count
//                    (mdi_infer_profile). An indented step is a part of the step after it (those of the bonds), or one
//                    of "branches", which run side by side: the sum of the steps is of those not indented.
//
// Without files it runs a default set from test_data (default_files).
//
// .gro, .pdb, .cif and .tpr files are loaded once, then each run infers into a copy of the system as loaded
// (md_system_infer; the copy is not timed), or its bonds alone (--bonds). Each line is a thread
// count: the best and median wall time of a run, the speedup over one thread and the page faults of the best run. What a
// run infers is summarized by a hash of the bonds, the chemistry, the molecules and their structures, the rings and the
// fragments (of the bonds alone with --bonds), which every count must agree on ('!' marks one that does not). The header
// names the benchmark's revision (bumped with every change of it) and the mdlib commit it was built from.

#include <md_system.h>
#include <md_infer.h>
#include <md_bond.h>
#include <md_ring.h>
#include <md_fragment.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_mmcif.h>
#include <md_tpr.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake
#include "bench_common.h"                 // The thread pool, the page faults
#include <md_infer_profile.h>             // The time of each step of md_system_infer (--stages)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// The benchmark's own revision: bumped with every change of it, so that a result names the build it came from
#define MD_INFER_BENCH_REV "R6"

#define MAX_COUNTS 32

// The steps of md_system_infer over the thread counts (--stages): the median time of each, per count
typedef struct stages_t {
    int         num_steps;
    const char* name[MDI_INFER_PROFILE_MAX_STEPS];
    double      ms[MDI_INFER_PROFILE_MAX_STEPS][MAX_REPS];  // Of each timed run of the current count
    double      median[MDI_INFER_PROFILE_MAX_STEPS][MAX_COUNTS];
    bool        seen[MDI_INFER_PROFILE_MAX_STEPS][MAX_COUNTS];
    bool        in_count[MDI_INFER_PROFILE_MAX_STEPS];          // Recorded by a run of the current count
} stages_t;

static int stage_index(stages_t* st, const char* name) {
    for (int i = 0; i < st->num_steps; ++i) if (!strcmp(st->name[i], name)) return i;
    if (st->num_steps == MDI_INFER_PROFILE_MAX_STEPS) return -1;
    st->name[st->num_steps] = name;
    return st->num_steps++;
}

// A count begins: none of its runs have any steps yet
static void stages_begin_count(stages_t* st) {
    md_memset(st->ms, 0, sizeof(st->ms));
    md_memset(st->in_count, 0, sizeof(st->in_count));
}

// A run's steps, those of the same name (written in two parts) summed
static void stages_add_run(stages_t* st, const mdi_infer_profile_t* p, int rep) {
    for (int i = 0; i < p->count; ++i) {
        const int s = stage_index(st, p->name[i]);
        if (s < 0) continue;
        st->ms[s][rep] += p->ms[i];
        st->in_count[s] = true;
    }
}

static void stages_end_count(stages_t* st, int count_idx, int reps) {
    for (int i = 0; i < st->num_steps; ++i) {
        if (!st->in_count[i]) continue;
        double v[MAX_REPS];
        md_memcpy(v, st->ms[i], (size_t)reps * sizeof(double));
        qsort(v, (size_t)reps, sizeof(double), cmp_double);
        st->median[i][count_idx] = v[reps / 2];
        st->seen[i][count_idx] = true;
    }
}

static void stages_print(const stages_t* st, const int counts[], int num_counts) {
    printf("  %-24s", "step (median ms)");
    for (int c = 0; c < num_counts; ++c) printf(" %8d", counts[c]);
    printf("\n");
    double total[MAX_COUNTS] = {0};
    for (int i = 0; i < st->num_steps; ++i) {
        printf("  %-24s", st->name[i]);
        for (int c = 0; c < num_counts; ++c) {
            if (st->seen[i][c]) printf(" %8.1f", st->median[i][c]);
            else                printf(" %8s", "-");
            if (st->name[i][0] != ' ') total[c] += st->median[i][c];
        }
        printf("\n");
    }
    printf("  %-24s", "sum of steps");
    for (int c = 0; c < num_counts; ++c) printf(" %8.1f", total[c]);
    printf("\n");
}

static uint64_t bonds_hash(uint64_t h, const md_bond_data_t* b) {
    h = fnv(h, &b->count, sizeof(b->count));
    if (b->count) {
        h = fnv(h, b->pairs, b->count * sizeof(md_atom_pair_t));
        h = fnv(h, b->flags, b->count * sizeof(md_bond_flags_t));
    }
    return h;
}

// What md_bond_infer gives: the same bytes as bonds_hash of the bonds of one range, so the hash is as it was
static uint64_t bond_list_hash(uint64_t h, const md_bond_list_t* b) {
    const size_t count = md_bond_list_count(b);
    h = fnv(h, &count, sizeof(count));
    if (count) {
        h = fnv(h, b->pairs, count * sizeof(md_atom_pair_t));
        h = fnv(h, b->flags, count * sizeof(md_bond_flags_t));
    }
    return h;
}

// What md_system_infer gave the system
static uint64_t inferred_hash(const md_system_t* sys) {
    uint64_t h = bonds_hash(1469598103934665603ull, &sys->bond);
    const size_t na = sys->atom.count;
    if (sys->atom.flags)          h = fnv(h, sys->atom.flags, na * sizeof(md_atom_flags_t));
    if (sys->atom.formal_charge)  h = fnv(h, sys->atom.formal_charge, na);
    if (sys->atom.hydrogen_count) h = fnv(h, sys->atom.hydrogen_count, na);
    const md_molecule_data_t* m = &sys->molecule;
    h = fnv(h, &m->count, sizeof(m->count));
    if (m->count) h = fnv(h, m->offset, (m->count + 1) * sizeof(uint32_t));
    if (m->count && m->atom_idx) h = fnv(h, m->atom_idx, na * sizeof(int32_t));
    h = fnv(h, &m->structure_count, sizeof(m->structure_count));
    if (m->structure_idx) h = fnv(h, m->structure_idx, m->count * sizeof(int32_t));
    const size_t nr = md_ring_data_count(&sys->ring);
    for (size_t r = 0; r < nr; ++r) h = fnv(h, md_ring_atoms(&sys->ring, r), md_ring_size(&sys->ring, r) * sizeof(int32_t));
    h = fnv(h, &sys->fragment.count, sizeof(sys->fragment.count));
    return h;
}

static bool load(md_system_t* sys, md_system_state_t* st, md_str_t path) {
    if (md_str_ends_with(path, MD_STR_LIT(".gro"))) return md_gro_system_load_from_file(sys, st, path, NULL) == MD_LOAD_OK;
    if (md_str_ends_with(path, MD_STR_LIT(".pdb"))) return md_pdb_system_load_from_file(sys, st, path, NULL) == MD_LOAD_OK;
    if (md_str_ends_with(path, MD_STR_LIT(".cif"))) return md_mmcif_system_load_from_file(sys, st, path, NULL) == MD_LOAD_OK;
    if (md_str_ends_with(path, MD_STR_LIT(".tpr"))) return md_tpr_system_load_from_file(sys, st, path, NULL) == MD_LOAD_OK;
    return false;
}

static void quiet(struct md_logger_o* inst, enum md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
}

// Run when no files are given: a spread of sizes, formats and chemistry
// (solvated systems, protein, DNA, lipids, coarse grained), relative to test_data
static const char* default_files[] = {
    "nvt_prod.gro",
    "centered.gro",
    "inside-md-pullout.gro",
    "nucl-dna.gro",
    "tryptophan-md.gro",
    "dppc64.pdb",
    "tubulin-A-B.pdb",
    "8g7u.cif",
    "1fez.cif",
    "tpr/peptide_tip3p.tpr",
    "martini/martini3_membrane.tpr",
};

#define NUM_DEFAULT_FILES (sizeof(default_files) / sizeof(default_files[0]))

static void usage(const char* exe) {
    printf("usage: %s [--reps n] [--threads a,b,c] [--reuse-memory] [--bonds] [--stages] [file ...]\n", exe);
    printf("without files it runs the default set in %s\n", MD_BENCHMARK_DATA_DIR);
}

int main(int argc, char** argv) {
    int reps = 5;
    bool reuse = false;
    bool bonds_only = false;
    bool stages = false;
    int counts[MAX_COUNTS];
    int num_counts = 0;
    const char* files[64];
    int num_files = 0;

    md_os_sys_info_t si = {0};
    md_os_sys_info_query(&si);
    const int cores = si.num_virtual_cores > 0 ? si.num_virtual_cores : 1;

    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--reps") && i + 1 < argc) {
            reps = atoi(argv[++i]);
        } else if (!strcmp(a, "--threads") && i + 1 < argc) {
            for (const char* c = argv[++i]; *c && num_counts < MAX_COUNTS; ) {
                const int n = atoi(c);
                if (n > 0) counts[num_counts++] = MD_MIN(n, MAX_THREADS);
                while (*c && *c != ',') ++c;
                if (*c == ',') ++c;
            }
        } else if (!strcmp(a, "--reuse-memory")) {
            reuse = true;
        } else if (!strcmp(a, "--bonds")) {
            bonds_only = true;
        } else if (!strcmp(a, "--stages")) {
            stages = true;
        } else if (a[0] == '-') {
            usage(argv[0]);
            return a[1] == 'h' ? 0 : 1;
        } else if (num_files < 64) {
            files[num_files++] = a;
        }
    }
    reps = MD_CLAMP(reps, 1, MAX_REPS);
    if (num_counts == 0) {
        for (int n = 1; n < cores && num_counts < MAX_COUNTS - 1; n *= 2) counts[num_counts++] = n;
        counts[num_counts++] = MD_MIN(cores, MAX_THREADS);
    }
    static char default_paths[NUM_DEFAULT_FILES][1024];
    if (num_files == 0) {
        for (size_t i = 0; i < NUM_DEFAULT_FILES; ++i) {
            snprintf(default_paths[i], sizeof(default_paths[i]), "%s/%s", MD_BENCHMARK_DATA_DIR, default_files[i]);
            files[num_files++] = default_paths[i];
        }
    }

    // Errors only: the notes of the steps on each of many runs are not the point
    static md_logger_i logger = { .log = quiet };
    md_log_unregister(default_logger);
    md_log_register(&logger);

    printf("md_infer_bench %s: %s over threads\n", MD_INFER_BENCH_REV, bonds_only ? "md_bond_infer" : "md_system_infer");
    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("CPU %d logical cores | %d timed run%s per count | %s\n\n", cores, reps, reps == 1 ? "" : "s",
           reuse ? "into the same memory (faulted in once)" : "into fresh memory each run");

    for (int f = 0; f < num_files; ++f) {
        // The system as loaded, which every run starts from
        md_allocator_i* base_alloc = md_vm_arena_create(MD_GIGABYTES(64));
        md_system_t base = { .alloc = base_alloc };
        md_system_state_t st = { .alloc = base_alloc };
        if (!load(&base, &st, md_str_from_cstr(files[f]))) {
            fprintf(stderr, "%s: cannot be loaded (.gro, .pdb, .cif or .tpr)\n", files[f]);
            md_vm_arena_destroy(base_alloc);
            continue;
        }

        printf("%s (%zu atoms)\n", files[f], base.atom.count);
        printf("  threads      best    median   speedup    faults   result\n");
        double base_ms = 0;
        uint64_t first_hash = 0;
        static stages_t stage;
        md_memset(&stage, 0, sizeof(stage));
        const bool profile = stages && !bonds_only;
        for (int c = 0; c < num_counts; ++c) {
            static bench_pool_t pool;
            bench_pool_init(&pool, (uint32_t)counts[c]);
            const md_executor_t* ex = counts[c] > 1 ? &pool.ex : NULL;
            double ms[MAX_REPS];
            uint64_t hash = 0;
            size_t num_bonds = 0;
            uint64_t best_faults = 0;
            double best_ms = 1e300;
            bool ok = true;
            if (profile) stages_begin_count(&stage);
            md_allocator_i* shared = reuse ? md_vm_arena_create(MD_GIGABYTES(64)) : NULL;
            // One untimed run first: the allocators and the caches as a run after another finds them
            for (int r = -1; r < reps && ok; ++r) {
                md_allocator_i* alloc = shared ? shared : md_vm_arena_create(MD_GIGABYTES(64));
                if (shared) md_vm_arena_reset(shared);
                md_system_t sys = { .alloc = alloc };
                md_bond_list_t bonds = {0};
                if (!bonds_only) ok = md_system_copy(&sys, &base);
                mdi_infer_profile_t prof = {0};
                if (profile) mdi_infer_profile = &prof;
                const uint64_t f0 = page_faults();
                const md_tick_t t0 = md_tick_now();
                if (ok) ok = bonds_only ? md_bond_infer(&bonds, &base, &st, alloc, ex) : md_system_infer(&sys, &st, ex);
                const md_tick_t t1 = md_tick_now();
                const uint64_t f1 = page_faults();
                mdi_infer_profile = NULL;
                if (ok && r >= 0 && profile) stages_add_run(&stage, &prof, r);
                if (ok && r >= 0) {
                    ms[r] = md_tick_to_milliseconds(t1 - t0);
                    if (ms[r] < best_ms) {
                        best_ms = ms[r];
                        best_faults = f1 - f0;
                    }
                }
                if (ok && r == -1) {
                    hash = bonds_only ? bond_list_hash(1469598103934665603ull, &bonds) : inferred_hash(&sys);
                    num_bonds = bonds_only ? md_bond_list_count(&bonds) : sys.bond.count;
                }
                if (!shared) md_vm_arena_destroy(alloc);
            }
            if (shared) md_vm_arena_destroy(shared);
            bench_pool_free(&pool);
            if (!ok) {
                printf("  %7d   failed\n", counts[c]);
                continue;
            }
            if (profile) stages_end_count(&stage, c, reps);
            qsort(ms, (size_t)reps, sizeof(double), cmp_double);
            const double best = ms[0], median = ms[reps / 2];
            if (c == 0) {
                base_ms = best;
                first_hash = hash;
            }
            printf("  %7d %9.1f %9.1f %8.2fx %9llu   %zu bonds %016llx%s\n", counts[c], best, median, base_ms / best,
                   (unsigned long long)best_faults, num_bonds, (unsigned long long)hash, hash == first_hash ? "" : " !");
        }
        if (profile) {
            printf("\n");
            stages_print(&stage, counts, num_counts);
        }
        printf("\n");
        md_vm_arena_destroy(base_alloc);
    }
    return 0;
}
