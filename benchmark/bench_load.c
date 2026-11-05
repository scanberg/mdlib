// md_load_bench: how a structure file loads over threads.
//
//   md_load_bench [--reps n] [--threads a,b,c] [--str] file ...
//
//   --reps <n>       timed loads per thread count, the best and the median are reported (default 5)
//   --threads <list> thread counts to load with (default 1, 2, 4, ... and every logical core)
//   --str            the file read into memory first and loaded from the string: the parse alone, not the reads
//
// .pdb files load with md_pdb_system_init_from_file, .gro with md_gro_system_init_from_file. Each line is a thread
// count: the best and median wall time of a load, the speedup over one thread and the bytes of the file per second. The
// system each count loads is summarized by a hash of its atoms, residues, bonds and types, which every count must agree
// on ('!' marks one that does not). The header names the mdlib commit the benchmark was built from.

#include <md_pdb.h>
#include <md_gro.h>
#include <md_system.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>

#include "../unittest/test_executor.h"   // the benchmark's own threads: mdlib owns none
#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_REPS 64
#define MAX_COUNTS 32

static uint64_t fnv(uint64_t h, const void* p, size_t n) {
    const uint8_t* b = (const uint8_t*)p;
    for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; }
    return h;
}

static uint64_t system_hash(const md_system_t* sys, const md_system_state_t* st) {
    uint64_t h = 1469598103934665603ull;
    const size_t na = sys->atom.count, nr = sys->residue.count;
    h = fnv(h, &na, sizeof(na));
    h = fnv(h, st->xyz, na * sizeof(vec3_t));
    h = fnv(h, sys->atom.type_idx, na * sizeof(md_particle_type_idx_t));
    for (size_t t = 0; t < sys->atom.type.count; ++t) {
        const str_t name = md_particle_type_name(&sys->atom.type, t);
        h = fnv(h, name.ptr, name.len);
        h = fnv(h, &sys->atom.type.z[t], sizeof(sys->atom.type.z[t]));
    }
    if (nr) h = fnv(h, sys->residue.atom_offset, (nr + 1) * sizeof(uint32_t));
    h = fnv(h, sys->residue.seq_id, nr * sizeof(md_sequence_id_t));
    h = fnv(h, sys->bond.pairs, sys->bond.count * sizeof(md_atom_pair_t));
    if (sys->instance.count) h = fnv(h, sys->instance.res_offset, (sys->instance.count + 1) * sizeof(uint32_t));
    return h;
}

static bool load(md_system_t* sys, md_system_state_t* st, str_t path, str_t text, bool gro, const md_executor_t* ex) {
    if (text.len) {
        return gro ? md_gro_system_init_from_str(sys, st, text, ex) : md_pdb_system_init_from_str(sys, st, text, MD_PDB_OPTION_NONE, ex);
    }
    return gro ? md_gro_system_init_from_file(sys, st, path, ex) : md_pdb_system_init_from_file(sys, st, path, MD_PDB_OPTION_NONE, ex);
}

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static void quiet(struct md_logger_o* inst, enum md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
}

int main(int argc, char** argv) {
    int reps = 5;
    bool from_str = false;
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
                if (n > 0) counts[num_counts++] = MIN(n, TEST_EXECUTOR_MAX_THREADS);
                while (*c && *c != ',') ++c;
                if (*c == ',') ++c;
            }
        } else if (!strcmp(a, "--str")) {
            from_str = true;
        } else if (a[0] == '-') {
            printf("usage: %s [--reps n] [--threads a,b,c] [--str] file ...\n", argv[0]);
            return a[1] == 'h' ? 0 : 1;
        } else if (num_files < 64) {
            files[num_files++] = a;
        }
    }
    reps = CLAMP(reps, 1, MAX_REPS);
    if (num_counts == 0) {
        for (int n = 1; n < cores && num_counts < MAX_COUNTS - 1; n *= 2) counts[num_counts++] = n;
        counts[num_counts++] = MIN(cores, TEST_EXECUTOR_MAX_THREADS);
    }
    if (num_files == 0) {
        printf("usage: %s [--reps n] [--threads a,b,c] [--str] file ...\n", argv[0]);
        return 1;
    }

    // Errors only: the loaders' notes on each of many loads are not the point
    static md_logger_i logger = { .log = quiet };
    md_log_unregister(default_logger);
    md_log_register(&logger);

    printf("md_load_bench: structure files loaded over threads\n");
    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("CPU %d logical cores | %d timed load%s per count | %s\n\n", cores, reps, reps == 1 ? "" : "s",
           from_str ? "from a string (the parse, not the reads)" : "from the file");

    for (int f = 0; f < num_files; ++f) {
        const str_t path = str_from_cstr(files[f]);
        const bool gro = str_ends_with(path, STR_LIT(".gro"));
        md_file_info_t info = {0};
        if (!md_file_info_extract_from_path(path, &info)) {
            fprintf(stderr, "%s: cannot be read\n", files[f]);
            continue;
        }
        str_t text = {0};
        if (from_str) text = load_textfile(path, md_get_heap_allocator());
        const double mb = (double)info.size / (1024.0 * 1024.0);

        printf("%s (%.1f MB)\n", files[f], mb);
        printf("  threads      best    median   speedup      MB/s   system\n");
        double base = 0;
        uint64_t first_hash = 0;
        for (int c = 0; c < num_counts; ++c) {
            test_executor_t te = test_executor_threads((uint32_t)counts[c]);
            const md_executor_t* ex = counts[c] > 1 ? test_executor(&te) : NULL;
            double ms[MAX_REPS];
            uint64_t hash = 0;
            size_t num_atoms = 0;
            bool ok = true;
            // One untimed load first: the page cache and the allocators as a load after another finds them
            for (int r = -1; r < reps && ok; ++r) {
                md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(64));
                md_system_t sys = { .alloc = alloc };
                md_system_state_t st = { .alloc = alloc };
                const md_tick_t t0 = md_tick_now();
                ok = load(&sys, &st, path, text, gro, ex);
                const md_tick_t t1 = md_tick_now();
                if (ok && r >= 0) ms[r] = md_tick_to_milliseconds(t1 - t0);
                if (ok && r == -1) {
                    hash = system_hash(&sys, &st);
                    num_atoms = sys.atom.count;
                }
                md_vm_arena_destroy(alloc);
            }
            if (!ok) {
                printf("  %7d   failed to load\n", counts[c]);
                continue;
            }
            qsort(ms, (size_t)reps, sizeof(double), cmp_double);
            const double best = ms[0], median = ms[reps / 2];
            if (c == 0) {
                base = best;
                first_hash = hash;
            }
            printf("  %7d %9.1f %9.1f %8.2fx %9.0f   %zu atoms %016llx%s\n", counts[c], best, median, base / best, mb / (best * 1e-3),
                   num_atoms, (unsigned long long)hash, hash == first_hash ? "" : " !");
        }
        printf("\n");
        if (text.ptr) str_free(text, md_get_heap_allocator());
    }
    return 0;
}
