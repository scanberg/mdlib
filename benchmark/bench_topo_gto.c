// md_topo_gto_bench
//
// Timings of the certified critical point search on a GTO density (md_topo_compute_extremum_graph_gto and
// its GPU sweep), per phase, on the QM files in test_data. Meant to be run as is, so that numbers from
// different machines and builds compare:
//
//   md_topo_gto_bench [options] [dataset ...]
//
//   --rho <v>       rho_min (default 1e-4, viamd's default)
//   --reps <n>      timed runs per configuration, the median is reported (default 3)
//   --threads <n>   CPU worker threads (default 0: every logical core)
//   --cpu / --gpu   only that path (default: both, GPU when there is a device)
//   --profile       also a GPU run that waits for every kernel to time it (per-kernel times, GEMM FLOP/s)
//   --gemm <v>      GPU GEMM tiling (default 0; the variants are listed by --gemm-sweep)
//   --gemm-sweep    per dataset, every GEMM tiling: median wall time, GEMM kernel time, and whether the
//                   result is identical to tiling 0's (it should be: the sums run in the same order)
//   --data <dir>    test_data directory (default: the source tree's)
//   --verbose       keep mdlib's info and debug log lines (default: errors only)
//
// A dataset is a path relative to the data directory or an existing file: *.molden, or *.h5 (VeloxChem,
// HDF5 builds). Without any, every built-in dataset that exists is run. The GPU path runs once untimed
// first (driver-side compilation, allocations). The GPU result is compared with the CPU one:
// same critical points (type, position within 1e-5 Bohr) and the same number of graph edges.

#include <md_system.h>
#include <md_gto.h>
#include <md_topo.h>
#include <md_molden.h>
#if defined(MD_HDF5)
#include <md_vlx.h>
#endif
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_log.h>
#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bench_rev.h"     // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#ifndef MD_BENCHMARK_DATA_DIR
#define MD_BENCHMARK_DATA_DIR "test_data"
#endif

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define MAX_REPS 32

// The .h5 inputs need a build with MD_ENABLE_HDF5=ON (which defines MD_HDF5); without it they are listed as skipped.
static const char* builtin_datasets[] = {
    "molden/h2o_ccpvdz.molden",
    "vlx/h2o.h5",
    "vlx/acro-xps.h5",
    "vlx/amide.h5",
    "vlx/mol.h5",
};

#if defined(MD_HDF5)
#define HAVE_H5 1
#else
#define HAVE_H5 0
#endif

typedef struct input_t {
    md_allocator_i*   alloc;
    md_system_t       sys;
    md_system_state_t state;
    md_gto_basis_t    basis;
    double*           density;
    float*            xyz;          // Bohr, per basis atom
    size_t            num_atoms;
    size_t            num_ao;
} input_t;

typedef struct run_t {
    double ms_total;
    md_topo_gto_info_t info;
    md_topo_extremum_graph_t graph;
} run_t;

static bool ends_with(const char* s, const char* suffix) {
    const size_t n = strlen(s), m = strlen(suffix);
    return n >= m && strcmp(s + n - m, suffix) == 0;
}

static bool file_exists(const char* path) {
    FILE* f = fopen(path, "rb");
    if (f) fclose(f);
    return f != NULL;
}

static bool input_load(input_t* in, const char* path) {
    memset(in, 0, sizeof(*in));
    in->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    in->sys = (md_system_t){ .alloc = in->alloc };
    in->state = (md_system_state_t){ .alloc = in->alloc };
    const str_t file = str_from_cstr(path);
    bool ok = false;
    if (ends_with(path, ".molden")) {
        ok = md_molden_system_init_from_file(&in->sys, &in->state, file);
    }
#if defined(MD_HDF5)
    else if (ends_with(path, ".h5")) {
        ok = md_vlx_system_init_from_file(&in->sys, &in->state, file);
    }
#endif
    else {
        printf("%s: unsupported format (.molden or .h5)\n", path);
        return false;
    }
    if (!ok) { fprintf(stderr, "%s: could not be read\n", path); return false; }

    if (!md_gto_basis_extract_attributes(&in->basis, &in->sys.attributes, in->alloc)) {
        fprintf(stderr, "%s: no GTO basis\n", path);
        return false;
    }
    in->num_ao = md_gto_basis_num_ao(&in->basis);
    in->num_atoms = md_gto_basis_num_atoms(&in->basis);

    // QTAIM is defined on the total density; a file with one spin channel publishes only alpha, which is then the total.
    const md_attribute_t* d = md_attributes_find(&in->sys.attributes, STR_LIT("orbital/total/density"));
    if (!d) d = md_attributes_find(&in->sys.attributes, STR_LIT("orbital/alpha/density"));
    if (!d || d->format.rank != 2 || d->format.shape[0] != in->num_ao || d->format.shape[1] != in->num_ao) {
        fprintf(stderr, "%s: no density matrix matching the basis (%zu AOs)\n", path, in->num_ao);
        return false;
    }
    in->density = (double*)md_alloc(in->alloc, sizeof(double) * in->num_ao * in->num_ao);
    if (md_attribute_extract_f64(in->density, in->num_ao * in->num_ao, d, md_attribute_slice_all(), md_unit_none()) != in->num_ao * in->num_ao) {
        fprintf(stderr, "%s: density matrix could not be read\n", path);
        return false;
    }

    const md_attribute_t* c = md_attributes_find(&in->sys.attributes, STR_LIT("qm/atom/coordinate"));   // Angstrom
    double* xyz = (double*)md_alloc(in->alloc, sizeof(double) * 3 * MAX(in->num_atoms, 1));
    if (!c || md_attribute_extract_f64(xyz, 3 * in->num_atoms, c, md_attribute_slice_all(), md_unit_none()) < 3 * in->num_atoms) {
        fprintf(stderr, "%s: no coordinates for the %zu basis atoms\n", path, in->num_atoms);
        return false;
    }
    in->xyz = (float*)md_alloc(in->alloc, sizeof(float) * 3 * MAX(in->num_atoms, 1));
    for (size_t i = 0; i < 3 * in->num_atoms; ++i) in->xyz[i] = (float)(xyz[i] * ANGSTROM_TO_BOHR);
    return true;
}

static void input_free(input_t* in) {
    if (in->alloc) md_arena_allocator_destroy(in->alloc);
    memset(in, 0, sizeof(*in));
}

// One run; the graph is kept for the comparisons.
static bool run_once(run_t* r, const md_topo_gto_desc_t* desc, void* stream) {
    memset(r, 0, sizeof(*r));
    r->graph.alloc = md_get_heap_allocator();
    const md_tick_t t0 = md_tick_now();
    bool ok;
#if MD_ENABLE_GPU
    if (stream) ok = md_topo_compute_extremum_graph_gto_gpu(&r->graph, &r->info, desc, (md_gpu_stream_t)stream);
    else
#endif
    ok = md_topo_compute_extremum_graph_gto(&r->graph, &r->info, desc);
    (void)stream;
    r->ms_total = md_tick_to_milliseconds(md_tick_now() - t0);
    return ok;
}

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

// Median over the runs of the double at byte offset 'off' in run_t ('info' members included).
static double median_at(const run_t* runs, int n, size_t off) {
    double v[MAX_REPS];
    for (int i = 0; i < n; ++i) v[i] = *(const double*)((const char*)&runs[i] + off);
    qsort(v, n, sizeof(double), cmp_double);
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
#define MED(field) median_at(runs, nrep, offsetof(run_t, field))

static bool same_bytes(const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b) {
    if (a->num_vertices != b->num_vertices || a->num_edges != b->num_edges) return false;
    if (a->num_vertices && (memcmp(a->vertices, b->vertices, sizeof(md_topo_vert_t) * a->num_vertices) ||
                            memcmp(a->types, b->types, sizeof(*a->types) * a->num_vertices))) return false;
    return !a->num_edges || !memcmp(a->edges, b->edges, sizeof(md_topo_edge_t) * a->num_edges);
}

// Same critical points (type, position within tol Bohr), matched one to one, and the same edge count.
static bool same_cps(const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b, double tol, double* out_max_d) {
    *out_max_d = 0.0;
    if (a->num_vertices != b->num_vertices || a->num_edges != b->num_edges) return false;
    bool* used = (bool*)calloc(MAX(b->num_vertices, 1), 1);
    bool ok = true;
    for (uint32_t i = 0; i < a->num_vertices && ok; ++i) {
        double best = HUGE_VAL; uint32_t bj = 0;
        for (uint32_t j = 0; j < b->num_vertices; ++j) {
            if (used[j] || a->types[i] != b->types[j]) continue;
            const double dx = a->vertices[i].x - b->vertices[j].x, dy = a->vertices[i].y - b->vertices[j].y, dz = a->vertices[i].z - b->vertices[j].z;
            const double dd = sqrt(dx * dx + dy * dy + dz * dz);
            if (dd < best) { best = dd; bj = j; }
        }
        if (best > tol) ok = false;
        else { used[bj] = true; if (best > *out_max_d) *out_max_d = best; }
    }
    free(used);
    return ok;
}

// Errors only: the readers and md_gpu log a lot at debug level, which would bury the tables.
static void quiet_log(struct md_logger_o* inst, md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
}
static md_logger_i quiet_logger = { NULL, quiet_log };

typedef struct result_t {
    const char* name;
    double cpu_ms, gpu_ms;
} result_t;

static void print_row(const char* mode, const run_t* runs, int nrep, bool gpu) {
    const md_topo_gto_info_t* i0 = &runs[0].info;
    uint32_t cnt[MD_TOPO_NUM_TYPES] = {0};
    md_topo_count_vertex_types(cnt, &runs[0].graph);
    const double sweep = MED(info.ms_sweep);
    const double wait = gpu ? MED(info.ms_sweep_gpu_wait) : 0.0, polish = gpu ? MED(info.ms_sweep_polish) : 0.0;
    const double cpu_levels = MED(info.ms_sweep_cpu);
    const double host = gpu ? sweep - wait - polish - cpu_levels : 0.0;
    printf("  %-4s %10.1f %8.1f %10.1f", mode, MED(ms_total), MED(info.ms_setup), sweep);
    if (gpu) printf(" %10.1f %8.1f %9.1f %8.1f", wait, polish, cpu_levels, host > 0 ? host : 0.0);
    else     printf(" %10s %8s %9.1f %8s", "-", "-", cpu_levels, "-");
    printf(" %8.1f %7.1f | %9llu %6u %5u | %u/%u/%u/%u %s\n", MED(info.ms_separatrices), MED(info.ms_clusters),
        (unsigned long long)i0->num_box_evals, gpu ? i0->num_escalated_boxes : 0u, gpu ? i0->num_gpu_dispatches : 0u,
        cnt[MD_TOPO_MAXIMUM], cnt[MD_TOPO_SPLIT_SADDLE], cnt[MD_TOPO_JOIN_SADDLE], cnt[MD_TOPO_MINIMUM],
        i0->complete ? "complete" : (i0->cancelled ? "cancelled" : "INCOMPLETE"));
}

int main(int argc, char** argv) {
    double rho_min = 1.0e-4;
    int reps = 3, threads = 0;
    bool want_cpu = true, want_gpu = true, profile = false, verbose = false, gemm_sweep = false;
    uint32_t gemm_variant = 0;
    const char* data_dir = MD_BENCHMARK_DATA_DIR;
    const char* names[64];
    int num_names = 0;

    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (!strcmp(a, "--rho") && i + 1 < argc) rho_min = atof(argv[++i]);
        else if (!strcmp(a, "--reps") && i + 1 < argc) reps = atoi(argv[++i]);
        else if (!strcmp(a, "--threads") && i + 1 < argc) threads = atoi(argv[++i]);
        else if (!strcmp(a, "--data") && i + 1 < argc) data_dir = argv[++i];
        else if (!strcmp(a, "--cpu")) { want_cpu = true; want_gpu = false; }
        else if (!strcmp(a, "--gpu")) { want_gpu = true; want_cpu = false; }
        else if (!strcmp(a, "--profile")) profile = true;
        else if (!strcmp(a, "--gemm") && i + 1 < argc) gemm_variant = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--gemm-sweep")) gemm_sweep = true;
        else if (!strcmp(a, "--verbose")) verbose = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            printf("usage: %s [--rho v] [--reps n] [--threads n] [--cpu|--gpu] [--profile] [--gemm v] [--gemm-sweep] [--data dir] [--verbose] [dataset ...]\n", argv[0]);
            return 0;
        }
        else if (a[0] == '-') { fprintf(stderr, "unknown option %s (see --help)\n", a); return 1; }
        else if (num_names < 64) names[num_names++] = a;
    }
    if (reps < 1) reps = 1;
    if (reps > MAX_REPS) reps = MAX_REPS;
    if (num_names == 0) {
        for (size_t i = 0; i < sizeof(builtin_datasets) / sizeof(builtin_datasets[0]); ++i) names[num_names++] = builtin_datasets[i];
    }

    if (!verbose && default_logger) {
        md_log_unregister(default_logger);
        md_log_register(&quiet_logger);
    }

    // --- environment
    md_os_sys_info_t si = {0};
    md_os_sys_info_query(&si);
#if defined(NDEBUG)
    const char* build = "release";
#else
    const char* build = "DEBUG (timings not representative)";
#endif
    void* stream = NULL;
    char gpu_name[300] = "none";
#if MD_ENABLE_GPU
    md_gpu_device_t dev = NULL;
    if (want_gpu) {
        dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .label = "md_topo_gto_bench" });
        if (dev) {
            md_gpu_device_info_t di = {0};
            md_gpu_device_info(dev, &di);
            snprintf(gpu_name, sizeof(gpu_name), "%s%s", di.name, di.is_discrete ? "" : " (integrated)");
            stream = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
            md_topo_gpu_initialize(dev);    // kernel creation is not part of any run
        } else {
            snprintf(gpu_name, sizeof(gpu_name), "none (%s)", md_gpu_last_error());
        }
    }
#else
    snprintf(gpu_name, sizeof(gpu_name), "none (MD_ENABLE_GPU=OFF)");
#endif
    if (want_gpu && !stream) want_gpu = false;

    printf("md_topo_gto_bench: certified critical points of the GTO density\n");
    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("build %s | CPU %d logical cores, %s | GPU %s\n", build, si.num_virtual_cores,
           threads > 0 ? "threads as given" : "all used", gpu_name);
    printf("rho_min %.3g | h_min 1e-4 | separatrices traced | %d timed run%s per configuration (GPU after one warm-up), median shown\n",
           rho_min, reps, reps == 1 ? "" : "s");
#if MD_ENABLE_GPU
    if (want_gpu) {
        if (gemm_variant >= md_topo_gto_gpu_gemm_variant_count()) {
            fprintf(stderr, "--gemm %u: there are %u tilings (0..%u)\n", gemm_variant, md_topo_gto_gpu_gemm_variant_count(),
                    md_topo_gto_gpu_gemm_variant_count() - 1);
            return 1;
        }
        printf("GEMM tiling %s\n", md_topo_gto_gpu_gemm_variant_name(gemm_variant));
    }
#endif
    printf("data %s\n", data_dir);
    printf("formats .molden%s\n\n", HAVE_H5 ? ", .h5 (VeloxChem)" : " only: this mdlib was configured without MD_ENABLE_HDF5, so the .h5 inputs are skipped");

    result_t results[64];
    int num_results = 0;

    for (int di = 0; di < num_names; ++di) {
        char path[1024];
        if (!HAVE_H5 && ends_with(names[di], ".h5")) { printf("%s: skipped, needs MD_ENABLE_HDF5=ON\n\n", names[di]); continue; }
        if (file_exists(names[di])) snprintf(path, sizeof(path), "%s", names[di]);
        else snprintf(path, sizeof(path), "%s/%s", data_dir, names[di]);
        if (!file_exists(path)) { printf("%s: not found at %s, skipped\n\n", names[di], path); continue; }

        input_t in;
        if (!input_load(&in, path)) { input_free(&in); printf("\n"); continue; }
        printf("%s: %zu atoms, %zu shells, %zu primitives, %zu Cartesian AOs\n", names[di], in.num_atoms,
               (size_t)in.basis.num_shells, (size_t)in.basis.num_primitives, in.num_ao);
        printf("  %-4s %10s %8s %10s %10s %8s %9s %8s %8s %7s | %9s %6s %5s | %s\n", "path", "total ms", "setup", "sweep",
               "gpu", "polish", "cpu lvls", "host", "separ.", "clust.", "cubes", "escal.", "disp.", "max/bcp/rcp/ccp");

        md_topo_gto_desc_t desc = {
            .basis = &in.basis,
            .atom_xyz = in.xyz,
            .density_matrix = in.density,
            .rho_min = rho_min,
            .h_min = 1.0e-4,
            .trace_separatrices = true,
            .num_threads = (uint32_t)threads,
            .gpu_gemm_variant = gemm_variant,
        };

        result_t res = { names[di], -1.0, -1.0 };
        run_t cpu_ref;  memset(&cpu_ref, 0, sizeof(cpu_ref));
        bool have_cpu = false;
        run_t* runs = (run_t*)calloc(reps, sizeof(run_t));

        for (int pass = 0; pass < 2; ++pass) {
            const bool gpu = pass == 1;
            if (gpu ? !want_gpu : !want_cpu) continue;
            void* s = gpu ? stream : NULL;
            if (gpu) {
                // first use of the kernels on this input: driver-side compilation, allocations
                run_t warm;
                run_once(&warm, &desc, s);
                md_topo_extremum_graph_free(&warm.graph);
            }
            bool ok = true, deterministic = true;
            for (int r = 0; r < reps; ++r) {
                ok &= run_once(&runs[r], &desc, s);
                if (r > 0 && !same_bytes(&runs[0].graph, &runs[r].graph)) deterministic = false;
            }
            const int nrep = reps;
            print_row(gpu ? "GPU" : "CPU", runs, nrep, gpu);
            if (!ok) printf("       the run reported failure\n");
            if (!deterministic) printf("       NOT DETERMINISTIC: the runs differ\n");
            if (gpu && !runs[0].info.used_gpu) printf("       the GPU was not used (fell back to the CPU)\n");
            {
                const md_topo_gto_info_t* i0 = &runs[0].info;
                const uint64_t ev = i0->num_box_evals, sk = i0->num_children_skipped;
                printf("       %llu children excluded by their parent's expansion (%.0f%% of %llu)",
                       (unsigned long long)sk, ev + sk ? 100.0 * (double)sk / (double)(ev + sk) : 0.0, (unsigned long long)(ev + sk));
                if (gpu && i0->num_gpu_batches) printf(", %llu GPU batches, %.0f local AOs on average", (unsigned long long)i0->num_gpu_batches,
                                                       (double)i0->num_gpu_rows / (double)i0->num_gpu_batches);
                printf("\n");
            }
            if (gpu) res.gpu_ms = MED(ms_total); else res.cpu_ms = MED(ms_total);
            if (!gpu) {
                cpu_ref = runs[0];
                runs[0].graph = (md_topo_extremum_graph_t){0};
                have_cpu = true;
            } else if (have_cpu) {
                double dmax = 0.0;
                const bool same = same_cps(&cpu_ref.graph, &runs[0].graph, 1.0e-5, &dmax);
                printf("       GPU vs CPU: %s (max position difference %.2g Bohr)\n", same ? "same critical points and edge count" : "MISMATCH", dmax);
            }
            for (int r = 0; r < reps; ++r) md_topo_extremum_graph_free(&runs[r].graph);
        }

#if MD_ENABLE_GPU
        if (profile && want_gpu) {
            desc.profile_gpu_kernels = true;
            run_t warm;
            run_once(&warm, &desc, stream);
            md_topo_extremum_graph_free(&warm.graph);
            for (int r = 0; r < reps; ++r) run_once(&runs[r], &desc, stream);
            const int nrep = reps;
            const double ao = MED(info.ms_gpu_ao), gemm = MED(info.ms_gpu_gemm), epi = MED(info.ms_gpu_epilogue), dec = MED(info.ms_gpu_decide);
            const double wait = MED(info.ms_sweep_gpu_wait);
            const double k = ao + gemm + epi + dec;
            const double gflop = runs[0].info.gpu_gemm_flop * 1.0e-9;
            printf("  GPU kernels, each waited for (sweep %.1f ms in this mode):\n", MED(info.ms_sweep));
            printf("       ao %.1f ms (%.0f%%) | gemm %.1f ms (%.0f%%, %.0f GFLOP/s) | epilogue %.1f ms (%.0f%%) | decide %.1f ms (%.0f%%) | uploads, readbacks, waits %.1f ms\n",
                   ao, 100.0 * ao / (k > 0 ? k : 1), gemm, 100.0 * gemm / (k > 0 ? k : 1), gemm > 0 ? gflop / (gemm * 1.0e-3) : 0.0,
                   epi, 100.0 * epi / (k > 0 ? k : 1), dec, 100.0 * dec / (k > 0 ? k : 1), wait - k > 0 ? wait - k : 0.0);
            for (int r = 0; r < reps; ++r) md_topo_extremum_graph_free(&runs[r].graph);
            desc.profile_gpu_kernels = false;
        }

        if (gemm_sweep && want_gpu) {
            // Per tiling: a warm-up (kernel creation, driver compilation), 'reps' normal runs for the wall
            // time, 'reps' profiled runs for the GEMM kernel time. GFLOP/s are counted at tiling 0's
            // operation count (the useful work; larger tiles pad more), so they compare across tilings.
            const uint32_t nv = md_topo_gto_gpu_gemm_variant_count();
            md_topo_extremum_graph_t ref = {0};
            uint64_t ref_cubes = 0;
            double ref_flop = 0.0, ref_total = 0.0, ref_gemm = 0.0;
            printf("  GEMM tilings: median total ms | GEMM kernel ms, each waited for | vs tiling 0\n");
            for (uint32_t v = 0; v < nv; ++v) {
                desc.gpu_gemm_variant = v;
                desc.profile_gpu_kernels = false;
                run_t warm;
                run_once(&warm, &desc, stream);
                md_topo_extremum_graph_free(&warm.graph);
                bool ok = true;
                for (int r = 0; r < reps; ++r) ok &= run_once(&runs[r], &desc, stream);
                const int nrep = reps;
                const double total = MED(ms_total);
                const uint64_t cubes = runs[0].info.num_box_evals;
                bool same = true;
                if (v == 0) {
                    ref = runs[0].graph;
                    runs[0].graph = (md_topo_extremum_graph_t){0};
                    ref_cubes = cubes;
                } else {
                    same = cubes == ref_cubes && same_bytes(&ref, &runs[0].graph);
                }
                for (int r = 0; r < reps; ++r) md_topo_extremum_graph_free(&runs[r].graph);

                desc.profile_gpu_kernels = true;
                for (int r = 0; r < reps; ++r) ok &= run_once(&runs[r], &desc, stream);
                const double gemm = MED(info.ms_gpu_gemm);
                if (v == 0) { ref_flop = runs[0].info.gpu_gemm_flop; ref_total = total; ref_gemm = gemm; }
                for (int r = 0; r < reps; ++r) md_topo_extremum_graph_free(&runs[r].graph);

                char cmp[96];
                if (v == 0) snprintf(cmp, sizeof(cmp), "reference");
                else snprintf(cmp, sizeof(cmp), "total %.2fx, gemm %.2fx, %s", total > 0 ? ref_total / total : 0.0,
                              gemm > 0 ? ref_gemm / gemm : 0.0, same ? "same result" : "RESULT DIFFERS");
                printf("    %8.1f ms | gemm %7.1f ms %6.0f GFLOP/s | %-38s | %s%s\n", total, gemm,
                       gemm > 0 ? ref_flop * 1.0e-9 / (gemm * 1.0e-3) : 0.0, cmp, md_topo_gto_gpu_gemm_variant_name(v),
                       ok ? "" : " (a run reported failure)");
                if (!same) printf("        %llu cubes (tiling 0: %llu)\n", (unsigned long long)cubes, (unsigned long long)ref_cubes);
            }
            md_topo_extremum_graph_free(&ref);
            desc.gpu_gemm_variant = gemm_variant;
            desc.profile_gpu_kernels = false;
        }
#endif
        if (have_cpu) md_topo_extremum_graph_free(&cpu_ref.graph);
        free(runs);
        input_free(&in);
        printf("\n");
        if (num_results < 64) results[num_results++] = res;
    }

    // --- summary
    printf("summary (total wall time per search, median)\n");
    printf("  %-28s %12s %12s %9s\n", "dataset", "CPU ms", "GPU ms", "CPU/GPU");
    for (int i = 0; i < num_results; ++i) {
        const result_t* r = &results[i];
        char c[32] = "-", g[32] = "-", x[32] = "-";
        if (r->cpu_ms >= 0) snprintf(c, sizeof(c), "%.1f", r->cpu_ms);
        if (r->gpu_ms >= 0) snprintf(g, sizeof(g), "%.1f", r->gpu_ms);
        if (r->cpu_ms > 0 && r->gpu_ms > 0) snprintf(x, sizeof(x), "%.2fx", r->cpu_ms / r->gpu_ms);
        printf("  %-28s %12s %12s %9s\n", r->name, c, g, x);
    }

#if MD_ENABLE_GPU
    if (dev) {
        md_topo_gpu_shutdown();
        md_gpu_device_destroy(dev);
    }
#endif
    return 0;
}
