// EXPERIMENT: a density matrix in 16 x 16 tiles of the AO pairs that meet (whose shells both reach some
// point at the evaluation cutoff), built on the GPU from the occupied orbitals (shaders/gto_dmat.slang).
//
// The pattern: shells ordered by radius class (2 bohr), then by atom (the file's order, along a chain),
// an atom's shells together; their AOs in that order cut into tiles of 16; a pair of tiles is kept when
// any of its shell pairs meets. Pairs that never meet are never read by an evaluation with the same
// cutoff, so leaving them out is that cutoff and nothing more.
//
//   gto_dmat.tiles_match_orbitals   every tile of mol.h5's SCF density against the double product
//   gto_dmat.benchmark              GTO_DMAT_FILE=<file in test_data/vlx> (protein.h5): the host steps
//                                   per frame (pattern, packing, upload) and the kernel, sparse and
//                                   with every tile (dense); GTO_DMAT_EPS (1e-6), GTO_DMAT_MARGIN (0),
//                                   GTO_DMAT_REPS (5)

#include "utest.h"
#include "vlx_test_util.h"

#include <md_gto.h>
#include <md_qm.h>
#include <core/md_os.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#include "gto_dmat_shaders.inl"

#define DM_T 16

typedef struct dm_args_t {          // mirrors DmatArgs
    uint32_t      num_mo, num_items, items_x, _pad0;
    md_gpu_addr_t P, w, items, cols, tiles;
} dm_args_t;

typedef struct dm_src_t {
    vlx_test_t     t;
    md_gto_basis_t basis;
    double*        xyz;             // bohr, per atom
    size_t         num_atoms, num_ao, num_mo;
    double*        rows;            // num_mo x num_ao, the basis order
    double*        w;
} dm_src_t;

typedef struct dm_pattern_t {
    size_t    n, nt, ntile;         // AOs, tiles per edge, tile pairs kept
    int*      perm;                 // tiled AO -> basis AO (-1 past n)
    size_t*   row_ptr;              // per row tile: its first entry
    uint32_t* cols;                 // per entry: its column tile (>= the row tile)
    uint32_t* items;                // 4 per work item: row tile, first entry, entries, 0
    size_t    num_items;
} dm_pattern_t;

static bool dm_load(dm_src_t* s, const char* file) {
    memset(s, 0, sizeof(*s));
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    if (!vlx_test_load(&s->t, str_from_cstr(path), MEGABYTES(1024))) return false;
    if (!qm_test_basis(&s->basis, &s->t)) return false;
    const md_attribute_t* ca = qm_test_attr(&s->t, STR_LIT("qm/atom/coordinate"));
    if (!ca) return false;
    s->num_atoms = md_attribute_element_count(&ca->format) / 3;
    s->xyz = (double*)malloc(sizeof(double) * 3 * s->num_atoms);
    md_attribute_extract_f64(s->xyz, 3 * s->num_atoms, ca, md_attribute_slice_all(), md_unit_none());
    for (size_t i = 0; i < 3 * s->num_atoms; ++i) s->xyz[i] *= QM_TEST_ANGSTROM_TO_BOHR;
    s->num_mo = md_qm_density_orbitals(&s->rows, &s->w, &s->num_ao, &s->t.sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
    return s->num_mo > 0 && s->num_ao == md_gto_basis_num_ao(&s->basis);
}

static void dm_free(dm_src_t* s) {
    if (s->rows) md_free(md_get_heap_allocator(), s->rows, sizeof(double) * s->num_mo * s->num_ao);
    if (s->w) md_free(md_get_heap_allocator(), s->w, sizeof(double) * s->num_mo);
    free(s->xyz);
    qm_test_free(&s->t);
}

static const uint64_t* dm_sort_key;
static int dm_cmp(const void* a, const void* b) {
    const uint64_t x = dm_sort_key[*(const int*)a], y = dm_sort_key[*(const int*)b];
    return x < y ? -1 : x > y;
}

// dense: every tile pair (the same kernel over the whole triangle)
static void dm_pattern(dm_pattern_t* p, const dm_src_t* s, double eps, double margin, bool dense) {
    memset(p, 0, sizeof(*p));
    const md_gto_basis_t* b = &s->basis;
    const size_t ns = b->num_shells;
    int* nc = (int*)malloc(sizeof(int) * ns);
    size_t* ao0 = (size_t*)malloc(sizeof(size_t) * ns);
    double* rad = (double*)malloc(sizeof(double) * ns);
    size_t n = 0;
    for (size_t sh = 0; sh < ns; ++sh) {
        const md_gto_shell_t* g = &b->shells[sh];
        nc[sh] = (int)md_gto_num_cart_ao(g->l);
        ao0[sh] = n;
        n += (size_t)nc[sh];
        // every AO of the shell below eps beyond the radius (as the iso kernels' cutoff)
        double nmax = 0.0, amin = 1e300;
        for (int q = 0; q < nc[sh]; ++q) { int i, j, k; md_gto_cart_ijk(&i, &j, &k, g->l, (uint32_t)q); nmax = fmax(nmax, fabs(md_gto_cart_norm_factor(i, j, k))); }
        for (uint32_t q = 0; q < g->num_primitives; ++q) amin = fmin(amin, b->alpha[g->primitive_offset + q]);
        double rr = sqrt(g->l / (2.0 * amin));
        for (;; rr += 0.05) {
            double t = 0.0;
            for (uint32_t q = 0; q < g->num_primitives; ++q) t += fabs(b->coeff[g->primitive_offset + q]) * exp(-b->alpha[g->primitive_offset + q] * rr * rr);
            if (nmax * pow(rr, g->l) * t < eps || rr > 1.0e3) break;
        }
        rad[sh] = rr;
    }
    // the tiled order: radius class, then atom, then the shell
    uint64_t* key = (uint64_t*)malloc(sizeof(uint64_t) * ns);
    int* order = (int*)malloc(sizeof(int) * ns);
    for (size_t sh = 0; sh < ns; ++sh) {
        key[sh] = ((uint64_t)(rad[sh] / 2.0) << 40) | ((uint64_t)b->shells[sh].atom_idx << 20) | (uint64_t)sh;
        order[sh] = (int)sh;
    }
    dm_sort_key = key;
    qsort(order, ns, sizeof(int), dm_cmp);
    size_t* tao0 = (size_t*)malloc(sizeof(size_t) * ns);
    p->n = n;
    p->nt = (n + DM_T - 1) / DM_T;
    p->perm = (int*)malloc(sizeof(int) * p->nt * DM_T);
    for (size_t i = 0; i < p->nt * DM_T; ++i) p->perm[i] = -1;
    size_t acc = 0;
    for (size_t i = 0; i < ns; ++i) {
        const int sh = order[i];
        tao0[sh] = acc;
        for (int q = 0; q < nc[sh]; ++q) p->perm[acc + (size_t)q] = (int)(ao0[sh] + (size_t)q);
        acc += (size_t)nc[sh];
    }
    const size_t nt = p->nt;
    uint8_t* on = (uint8_t*)calloc(nt * nt, 1);
    if (dense) {
        for (size_t I = 0; I < nt; ++I) for (size_t J = I; J < nt; ++J) on[I * nt + J] = 1;
    } else {
        for (size_t sh = 0; sh < ns; ++sh) {
            const double* a = s->xyz + 3 * b->shells[sh].atom_idx;
            for (size_t u = sh; u < ns; ++u) {
                const double* c = s->xyz + 3 * b->shells[u].atom_idx;
                const double dx = a[0] - c[0], dy = a[1] - c[1], dz = a[2] - c[2], lim = rad[sh] + rad[u] + margin;
                if (dx * dx + dy * dy + dz * dz > lim * lim) continue;
                for (size_t i = tao0[sh] / DM_T; i <= (tao0[sh] + (size_t)nc[sh] - 1) / DM_T; ++i)
                    for (size_t j = tao0[u] / DM_T; j <= (tao0[u] + (size_t)nc[u] - 1) / DM_T; ++j)
                        on[(i < j ? i : j) * nt + (i < j ? j : i)] = 1;
            }
        }
    }
    p->row_ptr = (size_t*)calloc(nt + 1, sizeof(size_t));
    for (size_t I = 0; I < nt; ++I) {
        size_t c = 0;
        for (size_t J = I; J < nt; ++J) c += on[I * nt + J];
        p->row_ptr[I + 1] = p->row_ptr[I] + c;
        p->num_items += (c + 7) / 8;
    }
    p->ntile = p->row_ptr[nt];
    p->cols = (uint32_t*)malloc(sizeof(uint32_t) * (p->ntile + 1));
    p->items = (uint32_t*)malloc(sizeof(uint32_t) * 4 * (p->num_items + 1));
    size_t e = 0, it = 0;
    for (size_t I = 0; I < nt; ++I) {
        for (size_t J = I; J < nt; ++J) if (on[I * nt + J]) p->cols[e++] = (uint32_t)J;
        for (size_t e0 = p->row_ptr[I]; e0 < p->row_ptr[I + 1]; e0 += 8) {
            const size_t cnt = p->row_ptr[I + 1] - e0 < 8 ? p->row_ptr[I + 1] - e0 : 8;
            p->items[4 * it + 0] = (uint32_t)I;
            p->items[4 * it + 1] = (uint32_t)e0;
            p->items[4 * it + 2] = (uint32_t)cnt;
            p->items[4 * it + 3] = 0;
            ++it;
        }
    }
    free(on); free(tao0); free(order); free(key); free(rad); free(ao0); free(nc);
}

static void dm_pattern_free(dm_pattern_t* p) {
    free(p->perm); free(p->row_ptr); free(p->cols); free(p->items);
    memset(p, 0, sizeof(*p));
}

// the orbitals by tile, as float: P[(tile * num_mo + k) * 16 + m]
static void dm_pack(float* P, float* w, const dm_src_t* s, const dm_pattern_t* p) {
    const size_t r = s->num_mo, n = s->num_ao;
    for (size_t t = 0; t < p->nt; ++t) {
        for (size_t k = 0; k < r; ++k) {
            const double* row = s->rows + k * n;
            float* dst = P + (t * r + k) * DM_T;
            for (int m = 0; m < DM_T; ++m) {
                const int ao = p->perm[t * DM_T + (size_t)m];
                dst[m] = ao >= 0 ? (float)row[ao] : 0.0f;
            }
        }
    }
    for (size_t k = 0; k < r; ++k) w[k] = (float)s->w[k];
}

typedef struct dm_gpu_t {
    md_gpu_device_t dev;
    md_gpu_stream_t s;
    md_gpu_kernel_t k;
} dm_gpu_t;

static bool dm_gpu_init(dm_gpu_t* g) {
    memset(g, 0, sizeof(*g));
    md_gpu_device_desc_t dd = {0};
    dd.label = "gto_dmat";
    g->dev = md_gpu_device_create(&dd);
    if (!g->dev) return false;
    g->s = md_gpu_stream_default(g->dev, MD_GPU_STREAM_COMPUTE);
    const md_gpu_kernel_desc_t kd = md_shader_gto_dmat_dmat_tiles_kernel();
    g->k = md_gpu_kernel_create(g->dev, &kd);
    return g->k != NULL;
}

static void dm_gpu_free(dm_gpu_t* g) {
    if (g->k) md_gpu_kernel_destroy(g->k);
    if (g->dev) md_gpu_device_destroy(g->dev);
    memset(g, 0, sizeof(*g));
}

typedef struct dm_run_t {
    double t_pack, t_upload, t_kernel;   // seconds (the kernel: the median of reps)
    double max_err, max_ref;             // over the entries checked
    size_t checked;
    bool   ok;
} dm_run_t;

// Packs, uploads and builds the tiles reps times; checks every check_stride-th tile against the double product
static dm_run_t dm_run(dm_gpu_t* g, const dm_src_t* s, const dm_pattern_t* p, int reps, size_t check_stride) {
    dm_run_t out = {0};
    const size_t r = s->num_mo;
    const size_t bytes_P = sizeof(float) * p->nt * r * DM_T;
    float* P = (float*)malloc(bytes_P);
    float* w = (float*)malloc(sizeof(float) * r);
    md_tick_t t0 = md_tick_now();
    dm_pack(P, w, s, p);
    out.t_pack = md_tick_to_seconds(md_tick_now() - t0);

    md_gpu_addr_t dP = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, bytes_P).gpu;
    md_gpu_addr_t dw = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, sizeof(float) * r).gpu;
    md_gpu_addr_t di = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * 4 * p->num_items).gpu;
    md_gpu_addr_t dc = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * p->ntile).gpu;
    md_gpu_mem_t  dt = md_gpu_malloc(g->s, MD_GPU_MEM_HOST_READ, sizeof(float) * 256 * p->ntile);
    if (!dP || !dw || !di || !dc || !dt.gpu) { free(P); free(w); return out; }
    md_gpu_stream_sync(g->s);
    t0 = md_tick_now();
    md_gpu_upload(g->s, dP, P, bytes_P);
    md_gpu_upload(g->s, dw, w, sizeof(float) * r);
    md_gpu_stream_sync(g->s);
    out.t_upload = md_tick_to_seconds(md_tick_now() - t0);
    md_gpu_upload(g->s, di, p->items, sizeof(uint32_t) * 4 * p->num_items);
    md_gpu_upload(g->s, dc, p->cols, sizeof(uint32_t) * p->ntile);

    const uint32_t gx = p->num_items < 65535 ? (uint32_t)p->num_items : 65535u;
    const uint32_t gy = (uint32_t)((p->num_items + gx - 1) / gx);
    dm_args_t a = {(uint32_t)r, (uint32_t)p->num_items, gx, 0, dP, dw, di, dc, dt.gpu};
    double times[64];
    const int nr = reps < 1 ? 1 : (reps > 64 ? 64 : reps);
    out.ok = MD_GPU_LAUNCH(g->s, g->k, md_gpu_grid(gx, gy, 1), a);     // warm up (pipelines, caches)
    md_gpu_stream_sync(g->s);
    for (int i = 0; i < nr && out.ok; ++i) {
        t0 = md_tick_now();
        out.ok = MD_GPU_LAUNCH(g->s, g->k, md_gpu_grid(gx, gy, 1), a);
        md_gpu_stream_sync(g->s);
        times[i] = md_tick_to_seconds(md_tick_now() - t0);
    }
    for (int i = 1; i < nr; ++i) for (int j = i; j > 0 && times[j] < times[j - 1]; --j) { const double x = times[j]; times[j] = times[j - 1]; times[j - 1] = x; }
    out.t_kernel = times[nr / 2];

    // the tiles against the double product of the same orbitals
    const float* tiles = (const float*)dt.cpu;
    const size_t n = s->num_ao;
    for (size_t I = 0; I < p->nt && out.ok; ++I) {
        for (size_t e = p->row_ptr[I]; e < p->row_ptr[I + 1]; ++e) {
            if (e % check_stride) continue;
            const size_t J = p->cols[e];
            for (int a0 = 0; a0 < DM_T; ++a0) for (int c0 = 0; c0 < DM_T; ++c0) {
                const int i = p->perm[I * DM_T + (size_t)a0], j = p->perm[J * DM_T + (size_t)c0];
                if (i < 0 || j < 0) continue;
                double ref = 0.0;
                for (size_t k = 0; k < r; ++k) ref += s->w[k] * s->rows[k * n + (size_t)i] * s->rows[k * n + (size_t)j];
                out.max_err = fmax(out.max_err, fabs((double)tiles[e * 256 + (size_t)a0 * DM_T + (size_t)c0] - ref));
                out.max_ref = fmax(out.max_ref, fabs(ref));
                out.checked++;
            }
        }
    }
    md_gpu_free(g->s, dP); md_gpu_free(g->s, dw); md_gpu_free(g->s, di); md_gpu_free(g->s, dc); md_gpu_free(g->s, dt.gpu);
    md_gpu_stream_sync(g->s);
    free(P); free(w);
    return out;
}

UTEST(gto_dmat, tiles_match_orbitals) {
    dm_gpu_t g;
    if (!dm_gpu_init(&g)) { dm_gpu_free(&g); UTEST_SKIP("no GPU device"); }
    dm_src_t s;
    ASSERT_TRUE(dm_load(&s, "mol.h5"));
    for (int dense = 0; dense < 2; ++dense) {
        dm_pattern_t p;
        dm_pattern(&p, &s, 1.0e-6, 0.0, dense != 0);
        const dm_run_t run = dm_run(&g, &s, &p, 1, 1);
        printf("  mol.h5 %s: %zu AOs, %zu orbitals, %zu of %zu tile pairs; max |err| %.3g of max |D| %.3g\n", dense ? "dense " : "sparse",
               s.num_ao, s.num_mo, p.ntile, p.nt * (p.nt + 1) / 2, run.max_err, run.max_ref);
        EXPECT_TRUE(run.ok);
        EXPECT_TRUE(run.checked > 0);
        EXPECT_LE(run.max_err, 1.0e-5 * run.max_ref);
        dm_pattern_free(&p);
    }
    dm_free(&s);
    dm_gpu_free(&g);
}

UTEST(gto_dmat, benchmark) {
    const char* file = getenv("GTO_DMAT_FILE");
    if (!file || !file[0]) UTEST_SKIP("set GTO_DMAT_FILE=<file in test_data/vlx>");
    const double eps = getenv("GTO_DMAT_EPS") ? atof(getenv("GTO_DMAT_EPS")) : 1.0e-6;
    const double margin = getenv("GTO_DMAT_MARGIN") ? atof(getenv("GTO_DMAT_MARGIN")) : 0.0;
    const int reps = getenv("GTO_DMAT_REPS") ? atoi(getenv("GTO_DMAT_REPS")) : 5;
    dm_gpu_t g;
    if (!dm_gpu_init(&g)) { dm_gpu_free(&g); UTEST_SKIP("no GPU device"); }
    md_gpu_device_info_t info;
    md_gpu_device_info(g.dev, &info);
    dm_src_t s;
    md_tick_t t0 = md_tick_now();
    ASSERT_TRUE(dm_load(&s, file));
    printf("\n  device: %s | %s: %zu atoms, %zu AOs, %zu occupied orbitals (loaded in %.1f s) | eps %g, margin %g bohr\n",
           info.name, file, s.num_atoms, s.num_ao, s.num_mo, md_tick_to_seconds(md_tick_now() - t0), eps, margin);
    for (int dense = 0; dense < 2; ++dense) {
        dm_pattern_t p;
        t0 = md_tick_now();
        dm_pattern(&p, &s, eps, margin, dense != 0);
        const double t_pattern = md_tick_to_seconds(md_tick_now() - t0);
        // check about 2000 tiles whatever the size
        const size_t stride = p.ntile > 2000 ? p.ntile / 2000 : 1;
        const dm_run_t run = dm_run(&g, &s, &p, reps, stride);
        const double fma = (double)p.ntile * 256.0 * (double)s.num_mo;
        printf("  %s: %zu of %zu tile pairs (%.1f%%), %.1f MB of tiles | host: pattern %.3f s, packing %.3f s, upload %.1f ms (%.1f MB)"
               " | kernel %.2f ms (median of %d), %.3g multiply-adds, %.1f G/s | max |err| %.2g of max |D| %.3g (%zu entries)\n",
               dense ? "dense " : "sparse", p.ntile, p.nt * (p.nt + 1) / 2, 100.0 * (double)p.ntile / (double)(p.nt * (p.nt + 1) / 2),
               256.0 * 4.0 * (double)p.ntile / 1e6, t_pattern, run.t_pack, run.t_upload * 1e3, 4.0 * (double)(p.nt * DM_T * s.num_mo) / 1e6,
               run.t_kernel * 1e3, reps, fma, fma / run.t_kernel / 1e9, run.max_err / run.max_ref, run.max_ref, run.checked);
        EXPECT_TRUE(run.ok);
        EXPECT_LE(run.max_err, 1.0e-4 * run.max_ref);
        dm_pattern_free(&p);
    }
    dm_free(&s);
    dm_gpu_free(&g);
}

#endif // MD_ENABLE_GPU
