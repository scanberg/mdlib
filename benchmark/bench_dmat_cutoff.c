// md_dmat_cutoff_bench
//
// What the cutoff of the one tile representation of the density matrix (md_gto_dmat_t) costs and loses: for each
// cutoff, the tiles it keeps (margin 1 bohr, viamd's), their bytes in float and double and the CPU build's time, and
// rho, ELF and LOL at points evaluated from them (the AOs screened at the same cutoff) against a reference made the
// same way at a negligible cutoff (1e-12, margin 0). CPU only, in double: the loss is the cutoff's, not float's.
//
//   points    uniform in the atoms' box +- 2 bohr; the errors are taken where the reference rho exceeds 1e-2, 1e-3
//             and 1e-4 (ELF's surfaces lie at rho of 1e-3 and above; below 1e-4 ELF is mostly its d_min floor)
//   ELF, LOL  md_qm_elf / md_qm_lol per spin (an open shell's alpha and beta tiles), rho_min 0, d_min MD_QM_ELF_D_MIN
//             (viamd's): the largest and the 99th percentile |difference|
//   rho       the largest |difference| / rho
//
//   md_dmat_cutoff_bench [--cutoffs a,b,..] [--points n] [--threads n] [--data dir] [file ...]
//
//   --cutoffs  default 1e-6,5e-7,3e-7,2e-7,1e-7
//   --points   default 20000
//
// A file is a path relative to <data>/vlx or an existing file. Without any: protein.h5 if present, else mol.h5.

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <md_gto_dmat.h>
#include <md_qm.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include "bench_rev.h"      // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake
#include "../unittest/test_executor.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Bumped with every change to this program's measurement or output
#define MD_DMAT_CUTOFF_BENCH_REV "dmat-cutoff r1"

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define MAX_CUTOFFS 12
#define NUM_RHO 3

static const double rho_thresholds[NUM_RHO] = { 1.0e-2, 1.0e-3, 1.0e-4 };

static double ms_since(md_tick_t t0) { return md_tick_to_milliseconds(md_tick_now() - t0); }

typedef struct src_t {
    md_allocator_i* arena;
    md_gto_basis_t  basis;
    float*          xyz;            // bohr
    size_t          num_atoms, num_ao;
    int             spins;          // 1: a closed shell, the total density; 2: alpha and beta
    double*         rows[2];        // each spin's (or the total's) orbitals, num_mo x num_ao, the basis order
    double*         w[2];
    size_t          num_mo[2];
} src_t;

static void src_free(src_t* s) {
    for (int k = 0; k < 2; ++k) {
        if (s->rows[k]) md_free(md_get_heap_allocator(), s->rows[k], sizeof(double) * s->num_mo[k] * s->num_ao);
        if (s->w[k]) md_free(md_get_heap_allocator(), s->w[k], sizeof(double) * s->num_mo[k]);
    }
    if (s->arena) md_arena_allocator_destroy(s->arena);
    memset(s, 0, sizeof(*s));
}

static bool src_load(src_t* s, const char* path) {
    memset(s, 0, sizeof(*s));
    s->arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = { .alloc = s->arena };
    md_system_state_t state = { .alloc = s->arena };
    if (!md_vlx_system_init_from_file(&sys, &state, md_str_from_cstr(path), NULL)) return false;
    const md_attribute_t* c = md_attributes_find(&sys.attributes, MD_STR_LIT("qm/atom/coordinate"));
    if (!c || !md_gto_basis_extract_attributes(&s->basis, &sys.attributes, s->arena)) {
        fprintf(stderr, "%s: no basis or geometry\n", path);
        return false;
    }
    s->num_atoms = md_attribute_value_count(&c->format);
    double* xyz = md_arena_allocator_push(s->arena, sizeof(double) * 3 * s->num_atoms);
    md_attribute_extract_f64(xyz, 3 * s->num_atoms, &sys.attributes, c, md_attribute_slice_all(), md_unit_none());
    s->xyz = md_arena_allocator_push(s->arena, sizeof(float) * 3 * s->num_atoms);
    for (size_t i = 0; i < 3 * s->num_atoms; ++i) s->xyz[i] = (float)(xyz[i] * ANGSTROM_TO_BOHR);
    // a closed shell (the same orbitals and occupations in both spins) as its total density, as viamd does
    size_t na = 0, nb = 0;
    double *ra = NULL, *wa = NULL, *rb = NULL, *wb = NULL;
    const size_t ma = md_qm_density_orbitals(&ra, &wa, &na, &sys, MD_STR_LIT("orbital/alpha/density"), md_get_heap_allocator());
    const size_t mb = md_qm_density_orbitals(&rb, &wb, &nb, &sys, MD_STR_LIT("orbital/beta/density"), md_get_heap_allocator());
    if (ma == 0 || mb == 0 || na != md_gto_basis_num_ao(&s->basis) || nb != na) {
        fprintf(stderr, "%s: no orbitals of the ground state density\n", path);
        return false;
    }
    s->num_ao = na;
    const bool closed = ma == mb && memcmp(ra, rb, sizeof(double) * ma * na) == 0 && memcmp(wa, wb, sizeof(double) * ma) == 0;
    s->rows[0] = ra; s->w[0] = wa; s->num_mo[0] = ma;
    s->rows[1] = rb; s->w[1] = wb; s->num_mo[1] = mb;
    s->spins = closed ? 1 : 2;
    if (closed) for (size_t k = 0; k < ma; ++k) wa[k] *= 2.0;     // the total
    return true;
}

typedef struct field_t {
    double* rho[2];
    double* grad[2];
    double* tau[2];
} field_t;

static void field_alloc(field_t* f, size_t n) {
    for (int k = 0; k < 2; ++k) {
        f->rho[k] = malloc(sizeof(double) * n);
        f->grad[k] = malloc(sizeof(double) * 3 * n);
        f->tau[k] = malloc(sizeof(double) * n);
    }
}

static void field_free(field_t* f) {
    for (int k = 0; k < 2; ++k) { free(f->rho[k]); free(f->grad[k]); free(f->tau[k]); }
}

typedef struct run_t {
    uint64_t tiles;
    double   ms_init, ms_build, ms_eval;
} run_t;

// rho, its gradient and tau of each spin at the points, from tiles at the cutoff
static bool evaluate(field_t* f, run_t* r, const src_t* s, double cutoff, double margin, const float* pts, size_t n, const md_executor_t* ex) {
    memset(r, 0, sizeof(*r));
    for (int k = 0; k < s->spins; ++k) {
        md_gto_dmat_t dm;
        md_tick_t t0 = md_tick_now();
        const md_gto_dmat_desc_t desc = { .basis = &s->basis, .atom_xyz = s->xyz, .cutoff = cutoff, .margin = margin };
        if (!md_gto_dmat_init(&dm, &desc, md_get_heap_allocator())) return false;
        r->ms_init += ms_since(t0);
        t0 = md_tick_now();
        if (!md_gto_dmat_build(&dm, s->rows[k], s->num_ao, 1, s->w[k], s->num_mo[k], ex)) { md_gto_dmat_free(&dm); return false; }
        r->ms_build += ms_since(t0);
        t0 = md_tick_now();
        const size_t got = md_gto_dmat_evaluate(f->rho[k], f->grad[k], f->tau[k], &dm, &s->basis, s->xyz, 0, pts, n, 0, ex);
        r->ms_eval += ms_since(t0);
        r->tiles = dm.num_tiles;
        md_gto_dmat_free(&dm);
        if (got != n) return false;
    }
    return true;
}

static void elf_lol(double* elf, double* lol, double* rho, const field_t* f, int spins, size_t n) {
    for (size_t p = 0; p < n; ++p) {
        const double* g = f->grad[0] + 3 * p;
        const double g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        if (spins == 1) {
            rho[p] = f->rho[0][p];
            elf[p] = md_qm_elf_closed_shell(f->rho[0][p], g2, f->tau[0][p], 0.0, MD_QM_ELF_D_MIN);
            lol[p] = md_qm_lol_closed_shell(f->rho[0][p], f->tau[0][p], 0.0);
        } else {
            const double* gb = f->grad[1] + 3 * p;
            const double g2b = gb[0] * gb[0] + gb[1] * gb[1] + gb[2] * gb[2];
            rho[p] = f->rho[0][p] + f->rho[1][p];
            elf[p] = md_qm_elf(f->rho[0][p], g2, f->tau[0][p], f->rho[1][p], g2b, f->tau[1][p], 0.0, MD_QM_ELF_D_MIN);
            lol[p] = md_qm_lol(f->rho[0][p], f->tau[0][p], f->rho[1][p], f->tau[1][p], 0.0);
        }
    }
}

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : x > y;
}

// The largest and the 99th percentile of |a - ref| where rho_ref > thr; NAN where no point is
static void err_stats(double* out_max, double* out_p99, const double* a, const double* ref, const double* rho_ref, double thr, size_t n, double* scratch) {
    size_t m = 0;
    for (size_t p = 0; p < n; ++p) if (rho_ref[p] > thr) scratch[m++] = fabs(a[p] - ref[p]);
    if (m == 0) { *out_max = *out_p99 = NAN; return; }
    qsort(scratch, m, sizeof(double), cmp_double);
    *out_max = scratch[m - 1];
    *out_p99 = scratch[(size_t)floor(0.99 * (double)(m - 1))];
}

typedef struct opts_t {
    double      cutoffs[MAX_CUTOFFS];
    int         num_cutoffs;
    uint32_t    points, threads;
    const char* data;
} opts_t;

static void run(const char* path, const opts_t* o, const md_executor_t* ex) {
    src_t s;
    if (!src_load(&s, path)) { src_free(&s); return; }
    const size_t n = o->points;
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (size_t a = 0; a < s.num_atoms; ++a) for (int k = 0; k < 3; ++k) {
        lo[k] = fminf(lo[k], s.xyz[3 * a + k] - 2.0f);
        hi[k] = fmaxf(hi[k], s.xyz[3 * a + k] + 2.0f);
    }
    float* pts = malloc(sizeof(float) * 3 * n);
    uint64_t x = 0x9E3779B97F4A7C15ull;
    for (size_t p = 0; p < 3 * n; ++p) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        const int k = (int)(p % 3);
        pts[p] = lo[k] + (hi[k] - lo[k]) * (float)((double)(x >> 11) * (1.0 / 9007199254740992.0));
    }
    const size_t occ = s.num_mo[0];
    printf("%s: %zu atoms, %zu AOs, %s, %zu orbitals%s | %zu points in its box +- 2 bohr\n", path, s.num_atoms, s.num_ao,
           s.spins == 1 ? "closed shell" : "open shell", occ, s.spins == 1 ? " (the total)" : " (alpha)", n);

    field_t ref, f;
    field_alloc(&ref, n);
    field_alloc(&f, n);
    double* elf_r = malloc(sizeof(double) * n), *lol_r = malloc(sizeof(double) * n), *rho_r = malloc(sizeof(double) * n);
    double* elf = malloc(sizeof(double) * n), *lol = malloc(sizeof(double) * n), *rho = malloc(sizeof(double) * n);
    double* scratch = malloc(sizeof(double) * n);
    run_t rr;
    if (!evaluate(&ref, &rr, &s, 1.0e-12, 0.0, pts, n, ex)) { fprintf(stderr, "the reference failed\n"); goto done; }
    elf_lol(elf_r, lol_r, rho_r, &ref, s.spins, n);
    const double total = (double)((s.num_ao + 15) / 16) * (double)((s.num_ao + 15) / 16 + 1) / 2;
    size_t cnt[NUM_RHO] = {0};
    for (size_t p = 0; p < n; ++p) for (int t = 0; t < NUM_RHO; ++t) cnt[t] += rho_r[p] > rho_thresholds[t];
    printf("  reference 1e-12, margin 0: %llu of %.0f tile pairs (%.1f MB in double) | build %.0f ms, %zu points %.0f ms | points with rho > 1e-2, 1e-3, 1e-4: %zu, %zu, %zu\n",
           (unsigned long long)rr.tiles, total, rr.tiles * 2048.0 / (1024.0 * 1024.0), rr.ms_build, n, rr.ms_eval, cnt[0], cnt[1], cnt[2]);
    printf("  %-8s %9s %6s %8s %8s %8s %8s | %-21s %-21s %-21s | %-21s %-21s %-21s | %9s %9s %9s\n", "cutoff", "tiles", "frac", "MB f32",
           "MB f64", "init ms", "build ms", "ELF >1e-2 max / p99", "ELF >1e-3", "ELF >1e-4", "LOL >1e-2 max / p99", "LOL >1e-3", "LOL >1e-4",
           "rho >1e-2", "rho >1e-3", "rho >1e-4");
    for (int c = 0; c < o->num_cutoffs; ++c) {
        run_t r;
        if (!evaluate(&f, &r, &s, o->cutoffs[c], 1.0, pts, n, ex)) { fprintf(stderr, "cutoff %.0e failed\n", o->cutoffs[c]); continue; }
        elf_lol(elf, lol, rho, &f, s.spins, n);
        printf("  %-8.1e %9llu %6.3f %8.1f %8.1f %8.0f %8.0f |", o->cutoffs[c], (unsigned long long)r.tiles, r.tiles / total,
               r.tiles * 1024.0 / (1024.0 * 1024.0), r.tiles * 2048.0 / (1024.0 * 1024.0), r.ms_init, r.ms_build);
        for (int q = 0; q < 2; ++q) {
            for (int t = 0; t < NUM_RHO; ++t) {
                double mx, p99;
                err_stats(&mx, &p99, q == 0 ? elf : lol, q == 0 ? elf_r : lol_r, rho_r, rho_thresholds[t], n, scratch);
                printf(" %9.2e / %9.2e", mx, p99);
            }
            printf(" |");
        }
        for (int t = 0; t < NUM_RHO; ++t) {
            double m = 0.0;
            for (size_t p = 0; p < n; ++p) if (rho_r[p] > rho_thresholds[t]) m = fmax(m, fabs(rho[p] - rho_r[p]) / rho_r[p]);
            printf(" %9.2e", m);
        }
        printf("\n");
    }
done:
    printf("\n");
    free(pts); free(scratch);
    free(elf_r); free(lol_r); free(rho_r); free(elf); free(lol); free(rho);
    field_free(&ref);
    field_free(&f);
    src_free(&s);
}

static void quiet_log(struct md_logger_o* inst, md_log_type_t type, const char* msg) {
    (void)inst;
    if (type == MD_LOG_TYPE_ERROR) fprintf(stderr, "%s\n", msg);
}
static md_logger_i quiet_logger = { NULL, quiet_log };

static bool file_exists(const char* p) {
    FILE* f = fopen(p, "rb");
    if (f) fclose(f);
    return f != NULL;
}

int main(int argc, char** argv) {
    opts_t o = { .points = 20000, .data = MD_BENCHMARK_DATA_DIR };
    const double def_cutoffs[] = { 1e-6, 5e-7, 3e-7, 2e-7, 1e-7 };
    o.num_cutoffs = 5;
    memcpy(o.cutoffs, def_cutoffs, sizeof(def_cutoffs));
    const char* files[32];
    int num_files = 0;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        const bool more = i + 1 < argc;
        if (!strcmp(a, "--cutoffs") && more) {
            o.num_cutoffs = 0;
            char buf[512];
            snprintf(buf, sizeof(buf), "%s", argv[++i]);
            for (char* t = strtok(buf, ","); t && o.num_cutoffs < MAX_CUTOFFS; t = strtok(NULL, ",")) o.cutoffs[o.num_cutoffs++] = atof(t);
        }
        else if (!strcmp(a, "--points") && more)  o.points = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--threads") && more) o.threads = (uint32_t)atoi(argv[++i]);
        else if (!strcmp(a, "--data") && more)    o.data = argv[++i];
        else if (a[0] == '-') {
            printf("usage: %s [--cutoffs a,b,..] [--points n] [--threads n] [--data dir] [file ...]\n", argv[0]);
            return 1;
        }
        else if (num_files < 32) files[num_files++] = a;
    }
    if (default_logger) {
        md_log_unregister(default_logger);
        md_log_register(&quiet_logger);
    }
    md_os_sys_info_t si = {0};
    md_os_sys_info_query(&si);
    test_executor_t exec = test_executor_threads(o.threads > 0 ? o.threads : (si.num_virtual_cores > 0 ? si.num_virtual_cores : 1));
    const md_executor_t* ex = test_executor(&exec);
    printf("md_dmat_cutoff_bench %s, mdlib %s | %u threads\n\n", MD_DMAT_CUTOFF_BENCH_REV, MD_BENCH_SOURCE_REV, ex->num_threads);

    char path[1024];
    if (num_files == 0) {
        snprintf(path, sizeof(path), "%s/vlx/protein.h5", o.data);
        if (!file_exists(path)) snprintf(path, sizeof(path), "%s/vlx/mol.h5", o.data);
        run(path, &o, ex);
    }
    for (int i = 0; i < num_files; ++i) {
        if (file_exists(files[i])) snprintf(path, sizeof(path), "%s", files[i]);
        else snprintf(path, sizeof(path), "%s/vlx/%s", o.data, files[i]);
        run(path, &o, ex);
    }
    return 0;
}
