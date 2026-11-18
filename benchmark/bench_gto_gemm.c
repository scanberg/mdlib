// md_gto_gemm_bench
//
// The evaluations on the shared GEMM core (gto_gemm_core.slang) against their reference kernels, on the GPU, on
// VeloxChem files:
//
//   density   MD_GTO_GPU_DENSITY_ALGO_REFERENCE, _GEMM (Phi, then a symmetric GEMM on the core)
//   ELF, LOL  MD_GTO_GPU_LOCALIZATION_ALGO_REFERENCE (a loop over the AO pairs per voxel), _GEMM (Phi with
//             gradients, a symmetric GEMM of four rows a voxel per spin on the core, then the finish)
//   orbitals  the sum of psi^2 of the --orbitals highest occupied: MD_GTO_GPU_ORBITAL_ALGO_REFERENCE, _SHELL (8 a
//             pass), _GEMM (Phi, then Phi x C on the core)
//   dmat      the density matrix's tiles built from the SCF density's orbitals (md_gto_gpu_dmat_build, on the core),
//             timed as the GPU time of the build's submission (the orbitals' upload and the kernel).
//
// MD_GTO_GEMM_PROFILE=1 (with --verbose) splits each GEMM run into its passes (count+plan, phi+dgather, gemm), each
// waited for: the totals are then not the normal ones.
//
// as viamd evaluates them: through a dmat (md_gto_gpu_dmat, margin 1 bohr), ELF and LOL with rho_min 0 and
// MD_QM_ELF_D_MIN, the density's basis at 1e-6 and ELF's at 1e-7 (viamd's ES_LOCALIZATION_CUTOFF). Each
// configuration runs once untimed (the kernels' creation, the driver's compilation), then --reps times; the time
// is the launch to the stream's completion (the counting pass, its read back, every piece). Each output is
// compared with the first algorithm's of its kind, so that a wrong result on some GPU shows next to its time.
//
//   md_gto_gemm_bench [options] [file ...]
//
//   --h <bohr>       grid spacing (default 0.25)
//   --margin <bohr>  the grid's reach past the atoms (default 4)
//   --plane          also a grid of one layer through the centre (the 2D kernels, as viamd's slices)
//   --spins 2        ELF and LOL from two spins (D/2 each: the same values, the open-shell path's cost)
//   --only <kind>    density, elf, lol, orbitals or dmat (repeatable; default all five)
//   --orbitals <n>   orbitals of the orbitals' run (default 32)
//   --no-reference   skip the reference algorithms (slow on large systems)
//   --reps <n>       timed runs per configuration, min and median shown (default 3)
//   --cutoff-density <c>, --cutoff-elf <c>
//   --device <sel>   GPU adapter: an index from --list-devices or part of its name; --list-devices
//   --data <dir>     test_data directory (default: the source tree's)
//   --verbose        keep mdlib's log (the GEMM plans: batches, pieces, multiply-adds)
//
// A file is a path relative to <data>/vlx or an existing file. Without any: mol.h5, then protein.h5 if present.

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <md_qm.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_grid.h>
#include <core/md_gpu.h>

#include "bench_rev.h"      // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Bumped with every change to this program's measurement or output
#define MD_GTO_GEMM_BENCH_REV "gemm-bench r7 (the kernels pruned: one per pass for both block shapes, ELF/LOL's 64-wide tile everywhere; r5-r6: --variants)"

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define MAX_REPS 32

typedef struct src_t {
    md_allocator_i* arena;
    md_gto_basis_t  basis;
    float*          xyz;        // bohr
    size_t          num_atoms, num_ao;
    double*         D;          // total density [num_ao][num_ao]
    double*         mo;         // the orbitals' run: num_mo rows of num_ao coefficients, the highest occupied
    size_t          num_mo;
    double*         occ_rows;   // the dmat's run: the SCF density's orbitals (md_qm_density_orbitals), heap
    double*         occ_w;
    size_t          num_occ;
    md_float3_t     lo, hi;     // the atoms' bounds
} src_t;

static void src_free(src_t* s) {
    if (s->occ_rows) md_free(md_get_heap_allocator(), s->occ_rows, sizeof(double) * s->num_occ * s->num_ao);
    if (s->occ_w) md_free(md_get_heap_allocator(), s->occ_w, sizeof(double) * s->num_occ);
    if (s->arena) md_arena_allocator_destroy(s->arena);
    memset(s, 0, sizeof(*s));
}

static bool src_load(src_t* s, const char* path, size_t want_mo) {
    memset(s, 0, sizeof(*s));
    s->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = { .alloc = s->arena };
    md_system_state_t state = { .alloc = s->arena };
    if (!md_vlx_system_init_from_file(&sys, &state, md_str_from_cstr(path), NULL)) return false;
    const md_attribute_t* c = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* d = md_attributes_find(&sys.attributes, MD_STR_LIT("orbital/total/density"));
    double scale = 1.0;
    if (!d) {
        d = md_attributes_find(&sys.attributes, MD_STR_LIT("orbital/alpha/density"));
        scale = 2.0;    // a closed shell's alpha density: half the total
    }
    if (!c || !d || !md_gto_basis_extract_attributes(&s->basis, &sys.attributes, s->arena)) {
        fprintf(stderr, "%s: no basis, geometry or density\n", path);
        return false;
    }
    s->num_atoms = md_attribute_value_count(&c->format);
    double* xyz = md_arena_allocator_push(s->arena, sizeof(double) * 3 * s->num_atoms);
    md_attribute_extract_f64(xyz, 3 * s->num_atoms, &sys.attributes, c, md_attribute_slice_all(), md_unit_none());
    s->xyz = md_arena_allocator_push(s->arena, sizeof(float) * 3 * s->num_atoms);
    s->lo = md_float3_set1(FLT_MAX);
    s->hi = md_float3_set1(-FLT_MAX);
    for (size_t i = 0; i < s->num_atoms; ++i) {
        const md_float3_t p = md_float3_set((float)(xyz[3 * i + 0] * ANGSTROM_TO_BOHR), (float)(xyz[3 * i + 1] * ANGSTROM_TO_BOHR), (float)(xyz[3 * i + 2] * ANGSTROM_TO_BOHR));
        memcpy(s->xyz + 3 * i, p.elem, sizeof(float) * 3);
        s->lo = md_float3_min(s->lo, p);
        s->hi = md_float3_max(s->hi, p);
    }
    s->num_ao = md_gto_basis_num_ao(&s->basis);
    s->D = md_arena_allocator_push(s->arena, sizeof(double) * s->num_ao * s->num_ao);
    md_attribute_extract_f64(s->D, s->num_ao * s->num_ao, &sys.attributes, d, md_attribute_slice_all(), md_unit_none());
    if (scale != 1.0) for (size_t i = 0; i < s->num_ao * s->num_ao; ++i) s->D[i] *= scale;

    // the dmat's run: the orbitals of the total density
    {
        size_t nao = 0;
        s->num_occ = md_qm_density_orbitals(&s->occ_rows, &s->occ_w, &nao, &sys, MD_STR_LIT("orbital/total/density"), md_get_heap_allocator());
        if (nao != s->num_ao) s->num_occ = 0;
    }

    // the orbitals' run: the want_mo highest occupied (alpha) orbitals
    const md_attribute_t* mo  = md_attributes_find(&sys.attributes, MD_STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* occ = md_attributes_find(&sys.attributes, MD_STR_LIT("orbital/alpha/occupation"));
    if (mo && occ && want_mo > 0) {
        const size_t n_occ = md_attribute_value_count(&occ->format);
        double* o = md_arena_allocator_push(s->arena, sizeof(double) * n_occ);
        md_attribute_extract_f64(o, n_occ, &sys.attributes, occ, md_attribute_slice_all(), md_unit_none());
        size_t homo = 0;
        for (size_t i = n_occ; i-- > 0;) if (o[i] > 0.0) { homo = i; break; }
        s->num_mo = want_mo < homo + 1 ? want_mo : homo + 1;
        s->mo = md_arena_allocator_push(s->arena, sizeof(double) * s->num_mo * s->num_ao);
        for (size_t k = 0; k < s->num_mo; ++k) {
            md_attribute_extract_f64(s->mo + k * s->num_ao, s->num_ao, &sys.attributes, mo, md_attribute_slice_1((uint32_t)(homo + 1 - s->num_mo + k)), md_unit_none());
        }
    }
    return true;
}

// What one evaluation needs on the device: a basis at a cutoff, the atoms, the dmats (one per spin), the volume
typedef struct gpu_t {
    md_gpu_stream_t    st;
    md_gto_gpu_basis_t basis;
    md_gpu_addr_t      atoms;
    md_gto_gpu_dmat_t  dmat[2];
    md_gpu_addr_t      coeff;       // the orbitals' run
    md_gpu_texture_t   tex;
    md_gpu_mem_t       rb;
    md_grid_t          grid;
    size_t             n;
} gpu_t;

static void gpu_free(gpu_t* g) {
    if (!g->st) return;
    for (int k = 0; k < 2; ++k) if (g->dmat[k]) md_gto_gpu_dmat_destroy(g->st, g->dmat[k]);
    if (g->coeff) md_gpu_free(g->st, g->coeff);
    if (g->tex) md_gpu_texture_destroy(g->tex);
    if (g->rb.gpu) md_gpu_free(g->st, g->rb.gpu);
    if (g->atoms) md_gpu_free(g->st, g->atoms);
    if (g->basis) md_gto_gpu_basis_destroy(g->st, g->basis);
    md_gpu_stream_sync(g->st);
    memset(g, 0, sizeof(*g));
}

static bool gpu_init(gpu_t* g, md_gpu_stream_t st, const src_t* s, const md_grid_t* grid, double cutoff, int spins, bool orbitals) {
    memset(g, 0, sizeof(*g));
    g->st = st;
    g->grid = *grid;
    g->n = (size_t)grid->dim[0] * grid->dim[1] * grid->dim[2];
    g->basis = md_gto_gpu_basis_create(st, &(md_gto_gpu_basis_desc_t){ .basis = &s->basis, .cutoff = cutoff });
    g->atoms = md_gpu_malloc(st, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(s->num_atoms)).gpu;
    if (!g->basis || !g->atoms) return false;
    float* p = (float*)md_gpu_upload_begin(st, g->atoms, md_gto_gpu_atom_buffer_size(s->num_atoms));
    md_gto_gpu_atom_pack(p, s->xyz, 0, s->num_atoms);
    md_gpu_upload_end(st);
    if (orbitals) {
        const size_t nc = md_gto_gpu_basis_num_cgtos(g->basis);
        g->coeff = md_gpu_malloc(st, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_mo(s->num_mo, nc)).gpu;
        if (!g->coeff || !s->num_mo) return false;
        const double** rows = (const double**)malloc(sizeof(double*) * s->num_mo);
        for (size_t k = 0; k < s->num_mo; ++k) rows[k] = s->mo + k * s->num_ao;
        float* c = (float*)md_gpu_upload_begin(st, g->coeff, md_gto_gpu_coeff_size_mo(s->num_mo, nc));
        md_gto_gpu_coeff_pack_mo(c, rows, NULL, s->num_mo, nc);
        md_gpu_upload_end(st);
        free(rows);
        spins = 0;      // no matrices
    }
    double* Dk = NULL;
    if (spins == 2) {
        Dk = (double*)malloc(sizeof(double) * s->num_ao * s->num_ao);
        for (size_t i = 0; i < s->num_ao * s->num_ao; ++i) Dk[i] = 0.5 * s->D[i];
    }
    for (int k = 0; k < spins; ++k) {
        g->dmat[k] = md_gto_gpu_dmat_create(st, &(md_gto_gpu_dmat_desc_t){ .basis = &s->basis, .atom_xyz = s->xyz, .cutoff = cutoff, .margin = 1.0 });
        if (!g->dmat[k] || !md_gto_gpu_dmat_set_matrix(st, g->dmat[k], Dk ? Dk : s->D)) { free(Dk); return false; }
    }
    free(Dk);
    g->tex = md_gpu_texture_create(st, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)grid->dim[0], .height = (uint32_t)grid->dim[1], .depth_or_layers = (uint32_t)grid->dim[2],
    });
    g->rb = md_gpu_malloc(st, MD_GPU_MEM_HOST_READ, sizeof(float) * g->n);
    md_gpu_stream_sync(st);
    return g->tex && g->rb.cpu;
}

typedef enum { KIND_DENSITY, KIND_ELF, KIND_LOL, KIND_ORBITALS, KIND_DMAT, KIND_COUNT } kind_t;
static const char* kind_name[KIND_COUNT] = { "density", "ELF", "LOL", "orbitals", "dmat" };

// One evaluation of kind with algorithm (density: md_gto_gpu_density_algo_t, orbitals: md_gto_gpu_orbital_algo_t, else
// md_gto_gpu_localization_algo_t)
static const src_t* run_src = NULL;     // the dmat's run: its orbitals
static double run_gpu_ms = -1.0;        // the dmat's run: its submission's GPU time, < 0 where the queue cannot time

static void run_once(gpu_t* g, kind_t kind, int algo, int spins, size_t num_mo) {
    if (kind == KIND_DMAT) {
        (void)algo;
        md_gpu_stream_sync(g->st);
        md_gto_gpu_dmat_build(g->st, g->dmat[0], run_src->occ_rows, run_src->occ_w, run_src->num_occ);
        const md_gpu_sync_t sync = md_gpu_stream_record(g->st);
        md_gpu_sync_wait(sync);
        double sec = 0.0;
        run_gpu_ms = -1.0;
        for (int i = 0; i < 50; ++i) {
            if (md_gpu_sync_gpu_seconds(sync, &sec)) { run_gpu_ms = sec * 1.0e3; break; }
            md_thread_sleep(1);
        }
    } else if (kind == KIND_ORBITALS) {
        md_gto_gpu_orbital_desc_t desc = {
            .basis = g->basis, .atom_xyz = g->atoms, .coeff = g->coeff, .out_tex = g->tex, .grid = &g->grid,
            .sample_offset = {0.5f, 0.5f, 0.5f}, .num_orbitals = num_mo, .eval_mode = MD_GTO_EVAL_MODE_PSI_SQUARED,
            .op = MD_GTO_OP_SET, .algorithm = (md_gto_gpu_orbital_algo_t)algo,
        };
        md_gto_gpu_orbital_launch(g->st, &desc);
    } else if (kind == KIND_DENSITY) {
        md_gto_gpu_density_desc_t desc = {
            .basis = g->basis, .atom_xyz = g->atoms, .dmat = g->dmat[0], .out_tex = g->tex, .grid = &g->grid,
            .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = (md_gto_gpu_density_algo_t)algo,
        };
        md_gto_gpu_density_launch(g->st, &desc);
    } else {
        md_gto_gpu_localization_desc_t desc = {
            .basis = g->basis, .atom_xyz = g->atoms, .num_spins = (uint32_t)spins, .dmat = { g->dmat[0], spins == 2 ? g->dmat[1] : NULL },
            .kind = kind == KIND_LOL ? MD_GTO_GPU_LOL : MD_GTO_GPU_ELF, .rho_min = 0.0, .d_min = MD_QM_ELF_D_MIN,
            .out_tex = g->tex, .grid = &g->grid, .sample_offset = {0.5f, 0.5f, 0.5f},
            .algorithm = (md_gto_gpu_localization_algo_t)algo,
        };
        md_gto_gpu_localization_launch(g->st, &desc);
    }
}

static bool file_exists(const char* path) {
    FILE* f = fopen(path, "rb");
    if (f) fclose(f);
    return f != NULL;
}

static int dcmp(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : x > y;
}

static void quiet_log(struct md_logger_o* inst, md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
}
static md_logger_i quiet_logger = { NULL, quiet_log };

// Every algorithm of each kind asked for, on grid; the first of a kind is the one the others are compared with
static void bench_grid(md_gpu_stream_t st, const src_t* s, const md_grid_t* grid, const char* label, const bool* kinds, bool reference,
                       int reps, int spins, double cutoff_density, double cutoff_elf) {
    const size_t n = (size_t)grid->dim[0] * grid->dim[1] * grid->dim[2];
    printf("  %s: %dx%dx%d = %zu voxels of %.3g bohr\n", label, grid->dim[0], grid->dim[1], grid->dim[2], n, grid->spacing.x);
    float* first = (float*)malloc(sizeof(float) * n);
    float* out = (float*)malloc(sizeof(float) * n);
    for (int kind = 0; kind < KIND_COUNT; ++kind) {
        if (!kinds[kind]) continue;
        const int sp = kind == KIND_DENSITY ? 1 : spins;
        gpu_t g;
        if (kind == KIND_ORBITALS && !s->num_mo) continue;
        if (kind == KIND_DMAT && !s->num_occ) continue;
        run_src = s;
        if (!gpu_init(&g, st, s, grid, (kind == KIND_ELF || kind == KIND_LOL) ? cutoff_elf : cutoff_density, sp, kind == KIND_ORBITALS)) {
            fprintf(stderr, "    %s: could not set up the evaluation\n", kind_name[kind]);
            gpu_free(&g);
            continue;
        }
        static const int density_algos[] = { MD_GTO_GPU_DENSITY_ALGO_REFERENCE, MD_GTO_GPU_DENSITY_ALGO_GEMM };
        static const char* density_names[] = { "reference", "gemm (core)" };
        static const int loc_algos[] = { MD_GTO_GPU_LOCALIZATION_ALGO_REFERENCE, MD_GTO_GPU_LOCALIZATION_ALGO_GEMM };
        static const char* loc_names[] = { "reference", "gemm (core)" };
        static const int mo_algos[] = { MD_GTO_GPU_ORBITAL_ALGO_REFERENCE, MD_GTO_GPU_ORBITAL_ALGO_SHELL, MD_GTO_GPU_ORBITAL_ALGO_GEMM };
        static const char* mo_names[] = { "reference", "shell", "gemm (core)" };
        static const int dmat_algos[] = { 0 };
        static const char* dmat_names[] = { "build (core)" };
        const int* algos = kind == KIND_DENSITY ? density_algos : kind == KIND_ORBITALS ? mo_algos : kind == KIND_DMAT ? dmat_algos : loc_algos;
        const char** names = kind == KIND_DENSITY ? density_names : kind == KIND_ORBITALS ? mo_names : kind == KIND_DMAT ? dmat_names : loc_names;
        const int num_algos = kind == KIND_ORBITALS ? 3 : kind == KIND_DMAT ? 1 : 2;
        if (kind == KIND_ORBITALS) printf("    (orbitals: the sum of psi^2 of the %zu highest occupied)\n", s->num_mo);
        if (kind == KIND_DMAT) {
            md_gto_gpu_dmat_info_t di;
            md_gto_gpu_dmat_info(&di, g.dmat[0]);
            printf("    (dmat: %zu orbitals, %zu of %zu tile pairs; times: the build's GPU time)\n",
                   s->num_occ, di.num_tiles, di.num_tiles_total);
        }
        double base_ms = 0.0;
        bool have_first = false;
        for (int ai = (reference || kind == KIND_DMAT) ? 0 : 1; ai < num_algos; ++ai) {
            run_once(&g, (kind_t)kind, algos[ai], sp, s->num_mo);      // untimed: kernels, the driver's compilation
            md_gpu_stream_sync(st);
            double ms[MAX_REPS];
            for (int r = 0; r < reps; ++r) {
                const md_tick_t t0 = md_tick_now();
                run_once(&g, (kind_t)kind, algos[ai], sp, s->num_mo);
                md_gpu_stream_sync(st);
                ms[r] = md_tick_to_milliseconds(md_tick_now() - t0);
                if (kind == KIND_DMAT && run_gpu_ms >= 0.0) ms[r] = run_gpu_ms;
            }
            qsort(ms, (size_t)reps, sizeof(double), dcmp);
            const double med = ms[reps / 2];
            md_gpu_copy_from_texture(st, g.rb.gpu, g.tex, NULL);
            md_gpu_stream_sync(st);
            memcpy(have_first ? out : first, g.rb.cpu, sizeof(float) * n);
            char cmp[160] = "";
            if (have_first) {
                double dmax = 0.0, fmax_ = 0.0;
                size_t bad = 0;
                for (size_t i = 0; i < n; ++i) {
                    if (!(out[i] == out[i])) { ++bad; continue; }
                    dmax = fmax(dmax, fabs((double)out[i] - (double)first[i]));
                    fmax_ = fmax(fmax_, fabs((double)first[i]));
                }
                if (kind == KIND_DENSITY || kind == KIND_ORBITALS || kind == KIND_DMAT) snprintf(cmp, sizeof(cmp), "max |diff| / max %.2e", fmax_ > 0 ? dmax / fmax_ : dmax);
                else snprintf(cmp, sizeof(cmp), "max |diff| %.2e", dmax);
                if (bad) snprintf(cmp + strlen(cmp), sizeof(cmp) - strlen(cmp), ", %zu NaN", bad);
            } else {
                have_first = true;
                base_ms = med;
                snprintf(cmp, sizeof(cmp), "(compared with)");
            }
            printf("    %-8s %-13s %10.2f ms min %10.2f ms median  %6.2fx  %.1f Mvox/s  %s\n", kind_name[kind], names[ai], ms[0], med,
                   base_ms > 0 ? base_ms / med : 1.0, (double)n / (med * 1.0e3), cmp);
        }
        gpu_free(&g);
    }
    free(first);
    free(out);
}

int main(int argc, char** argv) {
    double h = 0.25, margin = 4.0, cutoff_density = 1.0e-6, cutoff_elf = 1.0e-7;
    int reps = 3, spins = 1;
    size_t want_mo = 32;
    bool plane = false, reference = true, verbose = false, list_devices = false;
    bool kinds[KIND_COUNT] = { false, false, false };
    bool any_kind = false;
    const char* device_sel = NULL;
    const char* data_dir = MD_BENCHMARK_DATA_DIR;
    const char* names[64];
    int num_names = 0;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if      (!strcmp(a, "--h") && i + 1 < argc) h = atof(argv[++i]);
        else if (!strcmp(a, "--margin") && i + 1 < argc) margin = atof(argv[++i]);
        else if (!strcmp(a, "--reps") && i + 1 < argc) reps = atoi(argv[++i]);
        else if (!strcmp(a, "--spins") && i + 1 < argc) spins = atoi(argv[++i]) == 2 ? 2 : 1;
        else if (!strcmp(a, "--orbitals") && i + 1 < argc) want_mo = (size_t)atoi(argv[++i]);
        else if (!strcmp(a, "--cutoff-density") && i + 1 < argc) cutoff_density = atof(argv[++i]);
        else if (!strcmp(a, "--cutoff-elf") && i + 1 < argc) cutoff_elf = atof(argv[++i]);
        else if (!strcmp(a, "--plane")) plane = true;
        else if (!strcmp(a, "--no-reference")) reference = false;
        else if (!strcmp(a, "--only") && i + 1 < argc) {
            const char* k = argv[++i];
            if      (!strcmp(k, "density")) kinds[KIND_DENSITY] = true;
            else if (!strcmp(k, "elf"))     kinds[KIND_ELF] = true;
            else if (!strcmp(k, "lol"))     kinds[KIND_LOL] = true;
            else if (!strcmp(k, "orbitals")) kinds[KIND_ORBITALS] = true;
            else if (!strcmp(k, "dmat"))    kinds[KIND_DMAT] = true;
            else { fprintf(stderr, "--only %s: density, elf, lol, orbitals or dmat\n", k); return 1; }
            any_kind = true;
        }
        else if (!strcmp(a, "--device") && i + 1 < argc) device_sel = argv[++i];
        else if (!strcmp(a, "--list-devices")) list_devices = true;
        else if (!strcmp(a, "--data") && i + 1 < argc) data_dir = argv[++i];
        else if (!strcmp(a, "--verbose")) verbose = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            printf("usage: %s [--h bohr] [--margin bohr] [--plane] [--spins 1|2] [--only density|elf|lol|orbitals|dmat]... [--orbitals n] [--no-reference] [--reps n] "
                   "[--cutoff-density c] [--cutoff-elf c] [--device sel] [--list-devices] [--data dir] [--verbose] [file ...]\n", argv[0]);
            return 0;
        }
        else if (a[0] == '-') { fprintf(stderr, "unknown option %s (see --help)\n", a); return 1; }
        else if (num_names < 64) names[num_names++] = a;
    }
    if (!any_kind) for (int k = 0; k < KIND_COUNT; ++k) kinds[k] = true;
    if (reps < 1) reps = 1;
    if (reps > MAX_REPS) reps = MAX_REPS;
    if (!verbose && default_logger) {
        md_log_unregister(default_logger);
        md_log_register(&quiet_logger);
    }

    if (list_devices) {
        md_gpu_adapter_info_t adapters[16];
        const uint32_t na = md_gpu_enumerate_adapters(adapters, 16);
        for (uint32_t i = 0; i < na && i < 16; ++i) printf("%u: %s%s%s\n", i, adapters[i].name, adapters[i].usable ? "" : ", not usable: ", adapters[i].usable ? "" : adapters[i].missing);
        return 0;
    }
    md_gpu_device_t dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .adapter = device_sel, .label = "md_gto_gemm_bench" });
    if (!dev) {
        fprintf(stderr, "no GPU device: %s\n", md_gpu_last_error());
        return 1;
    }
    md_gpu_device_info_t di = {0};
    md_gpu_device_info(dev, &di);
    md_gpu_stream_t st = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_initialize(dev);
#if defined(NDEBUG)
    const char* build = "release";
#else
    const char* build = "DEBUG (timings not representative)";
#endif
    printf("md_gto_gemm_bench %s, mdlib %s\n", MD_GTO_GEMM_BENCH_REV, MD_BENCH_SOURCE_REV);
    printf("  %s\n  %s\n  %s\n  ELF/LOL reference %s\n", MD_GTO_GPU_GEMM_CORE_REVISION, MD_GTO_GPU_GEMM_PIECES_REVISION, MD_GTO_GPU_PLANE_REVISION, MD_GTO_GPU_LOCALIZATION_REVISION);
    printf("build %s | GPU %s%s (driver %s)%s\n", build, di.name, di.is_discrete ? "" : " (integrated)", di.driver,
           di.type == MD_GPU_DEVICE_TYPE_CPU ? " -- a CPU emulation: correctness only, these times mean nothing" : "");
    if (getenv("MD_GTO_GEMM_PROFILE")) printf("MD_GTO_GEMM_PROFILE is set: every GEMM pass is waited for, so the GEMM totals are slower than normal\n");
    printf("cutoffs: density and orbitals %.0e, ELF/LOL %.0e | ELF/LOL rho_min 0, d_min %.0e, %d spin%s | dmat, margin 1 bohr | %d timed runs after one untimed, min and median\n\n",
           cutoff_density, cutoff_elf, MD_QM_ELF_D_MIN, spins, spins == 2 ? "s (D/2 each)" : " (closed shell)", reps);

    static const char* defaults[] = { "mol.h5", "protein.h5" };
    const bool given = num_names > 0;
    if (!given) { names[0] = defaults[0]; names[1] = defaults[1]; num_names = 2; }
    for (int fi = 0; fi < num_names; ++fi) {
        char path[1024];
        snprintf(path, sizeof(path), "%s", names[fi]);
        if (!file_exists(path)) snprintf(path, sizeof(path), "%s/vlx/%s", data_dir, names[fi]);
        if (!file_exists(path)) {
            if (given) fprintf(stderr, "%s: not found\n", names[fi]);
            continue;
        }
        src_t s;
        if (!src_load(&s, path, want_mo)) { src_free(&s); continue; }
        printf("%s: %zu atoms, %zu AOs\n", names[fi], s.num_atoms, s.num_ao);
        const md_float3_t lo = md_float3_sub1(s.lo, (float)margin), hi = md_float3_add1(s.hi, (float)margin);
        const md_grid_t cube = {
            .orientation = md_float3x3_ident(),
            .origin  = lo,
            .spacing = md_float3_set1((float)h),
            .dim     = { (int)ceil((hi.x - lo.x) / h), (int)ceil((hi.y - lo.y) / h), (int)ceil((hi.z - lo.z) / h) },
        };
        bench_grid(st, &s, &cube, "volume", kinds, reference, reps, spins, cutoff_density, cutoff_elf);
        if (plane) {
            md_grid_t pl = cube;
            pl.origin.z = 0.5f * (lo.z + hi.z);
            pl.dim[2] = 1;
            bench_grid(st, &s, &pl, "plane through the centre", kinds, reference, reps, spins, cutoff_density, cutoff_elf);
        }
        printf("\n");
        src_free(&s);
    }
    md_gto_gpu_shutdown();
    md_gpu_device_destroy(dev);
    return 0;
}
