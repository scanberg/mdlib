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
//   --no-warmup     no GPU warm-up run before the timed ones (for inputs too large to run twice; the first
//                   run then includes the driver's kernel compilation)
//   --dump <prefix> write each run's critical points to <prefix>.cpu.txt / <prefix>.gpu.txt ("type x y z rho")
//   --compare <f>   list the critical points that differ from a --dump file (another run, build or machine)
//   --profile       also a GPU run that waits for every kernel to time it (per-kernel times, GEMM FLOP/s),
//                   and the kernels' driver statistics (registers, shared memory) where the driver gives them
//   --device <sel>  GPU adapter: an index from --list-devices or a case-insensitive part of its name
//                   ("intel", "3080"); without it MD_GPU_DEVICE is honoured, then --prefer
//   --prefer <p>    high-performance (default: discrete first) or low-power (integrated first)
//   --list-devices  list the GPU adapters and exit
//   --data <dir>    test_data directory (default: the source tree's)
//   --verbose       keep mdlib's info and debug log lines (default: errors only)
//
// A dataset is a path relative to the data directory or an existing file: *.molden, or *.h5 (VeloxChem,
// HDF5 builds). Without any, every built-in dataset that exists is run. The GPU path runs once untimed
// first (driver-side compilation, allocations). The GPU result is compared with the CPU one:
// same critical points (type, position within 1e-5 Bohr), and the edges both have (with the separatrices
// traced in fp32 on the GPU a few can differ, where a path runs close to the boundary between two basins).
// Each run prints hashes of its critical points and edges, and each input a fingerprint of the geometry,
// basis and density matrix, so runs on other machines or builds compare too.

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
    str_t dpath = STR_LIT("orbital/total/density");
    const md_attribute_t* d = md_attributes_find(&in->sys.attributes, dpath);
    if (!d) { dpath = STR_LIT("orbital/alpha/density"); d = md_attributes_find(&in->sys.attributes, dpath); }
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

// FNV-1a of the bytes: compares results across runs and builds (a run's critical points and edges as printed).
static uint64_t fnv64(uint64_t h, const void* p, size_t n) {
    const uint8_t* b = (const uint8_t*)p;
    for (size_t i = 0; i < n; ++i) h = (h ^ b[i]) * 1099511628211ull;
    return h;
}

// Same critical points (type, position within tol); map[i]: a's vertex i in b.
static bool same_cps(const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b, double tol, double* out_max_d, uint32_t* map) {
    *out_max_d = 0.0;
    if (a->num_vertices != b->num_vertices) return false;
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
        else { used[bj] = true; map[i] = bj; if (best > *out_max_d) *out_max_d = best; }
    }
    free(used);
    return ok;
}

static const char* cp_type_name(int t) {
    return t == MD_TOPO_MAXIMUM ? "max" : t == MD_TOPO_SPLIT_SADDLE ? "bcp" : t == MD_TOPO_JOIN_SADDLE ? "rcp" : t == MD_TOPO_MINIMUM ? "ccp" : "?";
}

// Prints the critical points of a that have no partner in b (same type, within tol Bohr; one to one),
// at most 'limit' of them. Returns how many there are.
static uint32_t list_unmatched(const char* label, const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b, double tol, uint32_t limit) {
    bool* used = (bool*)calloc(MAX(b->num_vertices, 1), 1);
    uint32_t n = 0;
    for (uint32_t i = 0; i < a->num_vertices; ++i) {
        double best = HUGE_VAL; uint32_t bj = 0;
        for (uint32_t j = 0; j < b->num_vertices; ++j) {
            if (used[j] || a->types[i] != b->types[j]) continue;
            const double dx = a->vertices[i].x - b->vertices[j].x, dy = a->vertices[i].y - b->vertices[j].y, dz = a->vertices[i].z - b->vertices[j].z;
            const double dd = sqrt(dx * dx + dy * dy + dz * dz);
            if (dd < best) { best = dd; bj = j; }
        }
        if (best <= tol) { used[bj] = true; continue; }
        if (n++ < limit) printf("         %s only: %s at (%.6f, %.6f, %.6f) Bohr, rho %.6g, nearest of its type %.3g Bohr away\n", label,
                                cp_type_name(a->types[i]), a->vertices[i].x, a->vertices[i].y, a->vertices[i].z, a->vertices[i].value, best);
    }
    free(used);
    return n;
}

// --dump / --compare: one critical point per line, "type x y z rho".
static void dump_cps(const char* path, const md_topo_extremum_graph_t* g) {
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "--dump: cannot write %s\n", path); return; }
    for (uint32_t i = 0; i < g->num_vertices; ++i)
        fprintf(f, "%d %.9g %.9g %.9g %.9g\n", (int)g->types[i], g->vertices[i].x, g->vertices[i].y, g->vertices[i].z, g->vertices[i].value);
    fclose(f);
}

static bool load_cps(const char* path, md_topo_extremum_graph_t* g) {
    FILE* f = fopen(path, "r");
    if (!f) return false;
    uint32_t cap = 1024, n = 0;
    md_topo_vert_t* v = (md_topo_vert_t*)malloc(sizeof(md_topo_vert_t) * cap);
    md_topo_critical_point_type_t* t = (md_topo_critical_point_type_t*)malloc(sizeof(*t) * cap);
    int ty; float x, y, z, r;
    while (fscanf(f, "%d %f %f %f %f", &ty, &x, &y, &z, &r) == 5) {
        if (n == cap) { cap *= 2; v = (md_topo_vert_t*)realloc(v, sizeof(*v) * cap); t = (md_topo_critical_point_type_t*)realloc(t, sizeof(*t) * cap); }
        v[n] = (md_topo_vert_t){ x, y, z, r };
        t[n] = (md_topo_critical_point_type_t)ty;
        n++;
    }
    fclose(f);
    memset(g, 0, sizeof(*g));
    g->num_vertices = n;
    g->vertices = v;
    g->types = t;
    return true;
}

static int u64_cmp(const void* pa, const void* pb) {
    const uint64_t x = *(const uint64_t*)pa, y = *(const uint64_t*)pb;
    return x < y ? -1 : (x > y ? 1 : 0);
}

// Edges of a (its vertices mapped into b's numbering) that b has too.
static uint32_t common_edges(const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b, const uint32_t* map) {
    uint64_t* ka = (uint64_t*)malloc(sizeof(uint64_t) * MAX(a->num_edges, 1));
    uint64_t* kb = (uint64_t*)malloc(sizeof(uint64_t) * MAX(b->num_edges, 1));
    for (uint32_t i = 0; i < a->num_edges; ++i) ka[i] = ((uint64_t)map[a->edges[i].from] << 32) | map[a->edges[i].to];
    for (uint32_t i = 0; i < b->num_edges; ++i) kb[i] = ((uint64_t)b->edges[i].from << 32) | b->edges[i].to;
    qsort(ka, a->num_edges, sizeof(uint64_t), u64_cmp);
    qsort(kb, b->num_edges, sizeof(uint64_t), u64_cmp);
    uint32_t n = 0;
    for (uint32_t i = 0, j = 0; i < a->num_edges && j < b->num_edges; ) {
        if (ka[i] == kb[j]) { n++; i++; j++; }
        else if (ka[i] < kb[j]) i++;
        else j++;
    }
    free(ka);
    free(kb);
    return n;
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
    bool want_cpu = true, want_gpu = true, profile = false, verbose = false;
    bool kernel_stats_printed = false;
    const char* dump_path = NULL;      // --dump <prefix>: <prefix>.cpu.txt / .gpu.txt
    const char* compare_path = NULL;   // --compare <file>: a --dump file from another run
    bool warmup = true;
    const char* device_sel = NULL;
    bool list_devices = false;
#if MD_ENABLE_GPU
    md_gpu_device_preference_t preference = MD_GPU_DEVICE_PREFER_DEFAULT;
#endif
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
        else if (!strcmp(a, "--no-warmup")) warmup = false;
        else if (!strcmp(a, "--dump") && i + 1 < argc) dump_path = argv[++i];
        else if (!strcmp(a, "--compare") && i + 1 < argc) compare_path = argv[++i];
        else if (!strcmp(a, "--device") && i + 1 < argc) device_sel = argv[++i];
        else if (!strcmp(a, "--list-devices")) list_devices = true;
        else if (!strcmp(a, "--prefer") && i + 1 < argc) {
            const char* pv = argv[++i];
#if MD_ENABLE_GPU
            if (!strcmp(pv, "high-performance")) preference = MD_GPU_DEVICE_PREFER_HIGH_PERFORMANCE;
            else if (!strcmp(pv, "low-power")) preference = MD_GPU_DEVICE_PREFER_LOW_POWER;
            else { fprintf(stderr, "--prefer %s: high-performance or low-power\n", pv); return 1; }
#else
            (void)pv;
#endif
        }
        else if (!strcmp(a, "--verbose")) verbose = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            printf("usage: %s [--rho v] [--reps n] [--threads n] [--cpu|--gpu] [--no-warmup] [--dump prefix] [--compare file] [--profile] [--device sel] [--prefer high-performance|low-power] [--list-devices] [--data dir] [--verbose] [dataset ...]\n", argv[0]);
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
    if (list_devices) {
        static const char* type_name[] = { "other", "discrete", "integrated", "virtual", "cpu" };
        md_gpu_adapter_info_t adapters[16];
        const uint32_t n = md_gpu_enumerate_adapters(adapters, 16);
        if (n == 0) printf("no GPU adapters (%s)\n", md_gpu_last_error() ? md_gpu_last_error() : "none found");
        for (uint32_t i = 0; i < n && i < 16; ++i) {
            const md_gpu_adapter_info_t* ad = &adapters[i];
            printf("%u: %s (0x%04X:0x%04X, %s)%s%s\n", i, ad->name, ad->vendor_id, ad->device_id,
                   (unsigned)ad->type < 5 ? type_name[ad->type] : "?", ad->usable ? "" : ", not usable: ", ad->usable ? "" : ad->missing);
        }
        return 0;
    }
    md_gpu_device_t dev = NULL;
    if (want_gpu) {
        dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .adapter = device_sel, .preference = preference, .label = "md_topo_gto_bench" });
        if (!dev && device_sel) {
            fprintf(stderr, "--device %s: %s (see --list-devices)\n", device_sel, md_gpu_last_error());
            return 1;
        }
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
    (void)device_sel;
    if (list_devices) { printf("no GPU support in this build (MD_ENABLE_GPU=OFF)\n"); return 0; }
    snprintf(gpu_name, sizeof(gpu_name), "none (MD_ENABLE_GPU=OFF)");
#endif
    if (want_gpu && !stream) want_gpu = false;

    printf("md_topo_gto_bench: certified critical points of the GTO density\n");
    printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
    printf("build %s | CPU %d logical cores, %s | GPU %s\n", build, si.num_virtual_cores,
           threads > 0 ? "threads as given" : "all used", gpu_name);
    printf("rho_min %.3g | h_min 1e-4 | separatrices traced | %d timed run%s per configuration (GPU %s), median shown\n",
           rho_min, reps, reps == 1 ? "" : "s", warmup ? "after one warm-up" : "without a warm-up");
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
        {
            // the input as the search sees it (geometry, basis, D): equal fingerprints, equal inputs, whatever the file
            uint64_t h = fnv64(1469598103934665603ull, in.xyz, sizeof(float) * 3 * in.num_atoms);
            h = fnv64(h, in.basis.shells, sizeof(*in.basis.shells) * in.basis.num_shells);
            h = fnv64(h, in.basis.alpha, sizeof(*in.basis.alpha) * in.basis.num_primitives);
            h = fnv64(h, in.basis.coeff, sizeof(*in.basis.coeff) * in.basis.num_primitives);
            h = fnv64(h, in.density, sizeof(double) * in.num_ao * in.num_ao);
            printf("  input fingerprint %016llx (geometry, basis, density matrix)\n", (unsigned long long)h);
        }
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
        };

        result_t res = { names[di], -1.0, -1.0 };
        run_t cpu_ref;  memset(&cpu_ref, 0, sizeof(cpu_ref));
        bool have_cpu = false;
        run_t* runs = (run_t*)calloc(reps, sizeof(run_t));

        for (int pass = 0; pass < 2; ++pass) {
            const bool gpu = pass == 1;
            if (gpu ? !want_gpu : !want_cpu) continue;
            void* s = gpu ? stream : NULL;
            if (gpu && warmup) {
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
                printf(", local AOs <= %u", i0->local_ao_bound);
                if (gpu && i0->gpu_density_bytes) {
                    const double dense = 4.0 * (double)in.num_ao * (double)in.num_ao;
                    printf(", D on the GPU %s %.3f MB (dense %.3f)", (double)i0->gpu_density_bytes == dense ? "dense" : "as tiles",
                           i0->gpu_density_bytes / 1e6, dense / 1e6);
                }
                printf("\n");
                {
                    const md_topo_extremum_graph_t* g0 = &runs[0].graph;
                    uint64_t hv = fnv64(1469598103934665603ull, g0->vertices, sizeof(md_topo_vert_t) * g0->num_vertices);
                    hv = fnv64(hv, g0->types, sizeof(*g0->types) * g0->num_vertices);
                    const uint64_t he = fnv64(1469598103934665603ull, g0->edges, sizeof(md_topo_edge_t) * g0->num_edges);
                    printf("       result hashes (to compare runs and builds): critical points %016llx, edges %016llx\n",
                           (unsigned long long)hv, (unsigned long long)he);
                    if (dump_path) {
                        char fn[1100];
                        snprintf(fn, sizeof(fn), "%s.%s.txt", dump_path, gpu ? "gpu" : "cpu");
                        dump_cps(fn, g0);
                        printf("       critical points written to %s\n", fn);
                    }
                    if (compare_path) {
                        md_topo_extremum_graph_t ref;
                        if (!load_cps(compare_path, &ref)) {
                            printf("       --compare: cannot read %s\n", compare_path);
                        } else {
                            printf("       against %s (%u critical points; position tolerance 1e-4 Bohr, as written):\n", compare_path, ref.num_vertices);
                            const uint32_t a = list_unmatched("there", &ref, g0, 1.0e-4, 40);
                            const uint32_t b = list_unmatched("here", g0, &ref, 1.0e-4, 40);
                            printf("         %u only there, %u only here\n", a, b);
                            free(ref.vertices);
                            free(ref.types);
                        }
                    }
                }
                printf("       screening: local radius %.2f Bohr (far %.2f, rung %u), at an atom %.0f local AOs (far %.0f), Hessian tail sampled %.2g",
                       i0->screening_radius_local, i0->screening_radius_far, i0->screening_rung, i0->screening_aos_local, i0->screening_aos_far,
                       i0->screening_tail_sampled);
                if (gpu) printf(", charged to a GPU batch at most %.2g", i0->screening_tail_max);
                printf("\n");
                if (i0->num_separatrix_paths) {
                    printf("       separatrices: %u paths (%u reached no CP), %llu point evaluations (the longest path %u)",
                           i0->num_separatrix_paths, i0->num_separatrix_open, (unsigned long long)i0->num_separatrix_evals, i0->max_separatrix_evals);
                    if (i0->separatrices_on_gpu) {
                        printf(", on the GPU in %u dispatches (half of the paths ended after %.1f s, 90%% after %.1f s, 99%% after %.1f s)",
                               i0->num_separatrix_dispatches, i0->ms_separatrix_gpu_ended[0] * 1e-3, i0->ms_separatrix_gpu_ended[1] * 1e-3,
                               i0->ms_separatrix_gpu_ended[2] * 1e-3);
                        if (i0->num_separatrix_cpu_paths) printf(", %u paths on the CPU (%.0f local AOs per evaluation)", i0->num_separatrix_cpu_paths, i0->separatrix_mean_local_aos);
                        printf("\n");
                    } else {
                        printf(", %.0f local AOs per evaluation (block %.0f)\n", i0->separatrix_mean_local_aos, i0->separatrix_mean_block_aos);
                    }
                }
            }
            if (gpu) res.gpu_ms = MED(ms_total); else res.cpu_ms = MED(ms_total);
            if (!gpu) {
                cpu_ref = runs[0];
                runs[0].graph = (md_topo_extremum_graph_t){0};
                have_cpu = true;
            } else if (have_cpu) {
                double dmax = 0.0;
                uint32_t* map = (uint32_t*)calloc(MAX(cpu_ref.graph.num_vertices, 1), sizeof(uint32_t));
                const bool same = same_cps(&cpu_ref.graph, &runs[0].graph, 1.0e-5, &dmax, map);
                if (!same) {
                    printf("       GPU vs CPU: critical points MISMATCH (%u vs %u)\n", runs[0].graph.num_vertices, cpu_ref.graph.num_vertices);
                    list_unmatched("CPU", &cpu_ref.graph, &runs[0].graph, 1.0e-5, 40);
                    list_unmatched("GPU", &runs[0].graph, &cpu_ref.graph, 1.0e-5, 40);
                } else {
                    const uint32_t ne_c = cpu_ref.graph.num_edges, ne_g = runs[0].graph.num_edges;
                    const uint32_t nc = common_edges(&cpu_ref.graph, &runs[0].graph, map);
                    printf("       GPU vs CPU: same critical points (max position difference %.2g Bohr); ", dmax);
                    if (nc == ne_c && nc == ne_g) printf("same edges (%u)\n", nc);
                    else printf("edges differ: %u common, %u only on the CPU, %u only on the GPU%s\n", nc, ne_c - nc, ne_g - nc,
                                runs[0].info.separatrices_on_gpu ? " (separatrices traced in fp32 on the GPU)" : " (MISMATCH)");
                }
                free(map);
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
            printf("       host: forming and launching chunks %.1f ms | outcomes and child batches %.1f ms\n",
                   MED(info.ms_sweep_host_launch), MED(info.ms_sweep_host_outcomes));
            if (!kernel_stats_printed) {   // the same for every input: once
                kernel_stats_printed = true;
                bool any = false;
                for (uint32_t i = 0; i < md_topo_gto_gpu_kernel_count(); ++i) {
                    const char* name = NULL;
                    md_gpu_kernel_info_t ki;
                    if (!md_topo_gto_gpu_kernel_info(md_gpu_stream_device(stream), i, &name, &ki)) continue;
                    if (!ki.stats[0]) continue;
                    if (!any) printf("  GPU kernels as the driver compiled them (registers and shared memory decide the groups per SM):\n");
                    any = true;
                    printf("       %-11s %4u threads: %s\n", name, ki.group_size[0] * ki.group_size[1] * ki.group_size[2], ki.stats);
                }
                if (!any) printf("  GPU kernels as the driver compiled them: the driver reports no statistics\n");
            }
            for (int r = 0; r < reps; ++r) md_topo_extremum_graph_free(&runs[r].graph);
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
