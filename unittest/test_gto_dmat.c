// md_gto_gpu_dmat: a density matrix in tiles of the AO pairs that meet, built on the GPU from orbitals or
// filled from a matrix, and evaluated by md_gto_gpu_density_launch.
//
//   gto_dmat.density_matches_packed   h2o and mol: the volume through a dmat (from the orbitals, and from
//                                     the matrix) against the packed matrix, every algorithm
//   gto_dmat.pattern_follows_motion   the pattern holds within margin / 2 and not beyond
//   gto_dmat.density_job_matches_launch  mol: md_gto_gpu_density_job_begin, polled, bitwise the launch (the dmat
//                                     and the packed matrix, every algorithm); a job ended before it is
//                                     issued writes nothing
//   gto_dmat.benchmark                GTO_DMAT_FILE=<file in test_data/vlx> (protein.h5): the pattern, the
//                                     build per frame (packing, upload, kernel), and a sub-volume through
//                                     the dmat against the packed matrix; GTO_DMAT_REPS (5),
//                                     GTO_DMAT_MARGIN (0), GTO_DMAT_GRID (32: voxels per edge, 0.4 bohr)

#include "utest.h"
#include "vlx_test_util.h"

#include <md_gto.h>
#include <md_qm.h>
#include <core/md_os.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>

typedef struct dm_src_t {
    vlx_test_t     t;
    md_gto_basis_t basis;
    float*         xyz;             // bohr, per atom
    size_t         num_atoms, num_ao, num_mo;
    double*        rows;            // the SCF density's orbitals: num_mo x num_ao, the basis order
    double*        w;
    double*        D;               // the same as a matrix (double), num_ao x num_ao
} dm_src_t;

static bool dm_load(dm_src_t* s, const char* file) {
    memset(s, 0, sizeof(*s));
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    if (!vlx_test_load(&s->t, str_from_cstr(path), MEGABYTES(1024))) return false;
    if (!qm_test_basis(&s->basis, &s->t)) return false;
    s->num_atoms = s->t.sys.atom.count;
    s->xyz = (float*)malloc(sizeof(float) * 3 * s->num_atoms);
    if (vlx_test_atom_xyz_bohr(s->xyz, 3 * s->num_atoms, &s->t) != s->num_atoms) return false;
    s->num_mo = md_qm_density_orbitals(&s->rows, &s->w, &s->num_ao, &s->t.sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
    if (s->num_mo == 0 || s->num_ao != md_gto_basis_num_ao(&s->basis)) return false;
    s->D = (double*)malloc(sizeof(double) * s->num_ao * s->num_ao);
    return s->D && md_gto_density_from_orbitals(s->D, s->rows, s->num_ao, 1, s->w, s->num_mo, s->num_ao);
}

static void dm_free(dm_src_t* s) {
    if (s->rows) md_free(md_get_heap_allocator(), s->rows, sizeof(double) * s->num_mo * s->num_ao);
    if (s->w) md_free(md_get_heap_allocator(), s->w, sizeof(double) * s->num_mo);
    free(s->D);
    free(s->xyz);
    qm_test_free(&s->t);
}

typedef struct dm_gpu_t {
    md_gpu_device_t    dev;
    md_gpu_stream_t    s;
    md_gto_gpu_basis_t basis;
    md_gpu_addr_t      atoms, coeff;     // the packed matrix
    md_gpu_texture_t   tex;
    md_gpu_mem_t       rb;
    md_grid_t          grid;
    size_t             num_vox;
} dm_gpu_t;

static bool dm_gpu_init(dm_gpu_t* g, const dm_src_t* s, const md_grid_t* grid) {
    memset(g, 0, sizeof(*g));
    g->dev = md_gpu_device_create(NULL);
    if (!g->dev) return false;
    g->s = md_gpu_stream_default(g->dev, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_initialize(g->dev);
    g->basis = md_gto_gpu_basis_create(g->s, &(md_gto_gpu_basis_desc_t){ .basis = &s->basis, .cutoff = 1.0e-6 });
    g->atoms = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(s->num_atoms)).gpu;
    g->coeff = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_density(s->num_ao)).gpu;
    if (!g->basis || !g->atoms || !g->coeff) return false;
    float* p = (float*)md_gpu_upload_begin(g->s, g->atoms, md_gto_gpu_atom_buffer_size(s->num_atoms));
    md_gto_gpu_atom_pack(p, s->xyz, 0, s->num_atoms);
    md_gpu_upload_end(g->s);
    p = (float*)md_gpu_upload_begin(g->s, g->coeff, md_gto_gpu_coeff_size_density(s->num_ao));
    md_gto_gpu_coeff_pack_density(p, s->D, s->num_ao);
    md_gpu_upload_end(g->s);
    g->grid = *grid;
    g->num_vox = (size_t)grid->dim[0] * grid->dim[1] * grid->dim[2];
    g->tex = md_gpu_texture_create(g->s, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)grid->dim[0], .height = (uint32_t)grid->dim[1], .depth_or_layers = (uint32_t)grid->dim[2],
    });
    g->rb = md_gpu_malloc(g->s, MD_GPU_MEM_HOST_READ, sizeof(float) * g->num_vox);
    return g->tex && g->rb.cpu;
}

static void dm_gpu_free(dm_gpu_t* g) {
    if (g->s) {
        md_gpu_free(g->s, g->rb.gpu);
        if (g->tex) md_gpu_texture_destroy(g->tex);
        md_gpu_free(g->s, g->coeff);
        md_gpu_free(g->s, g->atoms);
        if (g->basis) md_gto_gpu_basis_destroy(g->s, g->basis);
        md_gpu_stream_sync(g->s);
        md_gto_gpu_shutdown();
    }
    if (g->dev) md_gpu_device_destroy(g->dev);
    memset(g, 0, sizeof(*g));
}

// The volume into out (num_vox floats) through the packed matrix (dmat NULL) or the dmat
static double dm_eval(dm_gpu_t* g, md_gto_gpu_dmat_t dmat, md_gto_gpu_density_algo_t algo, float* out) {
    md_gto_gpu_density_desc_t desc = {
        .basis = g->basis, .atom_xyz = g->atoms, .coeff = dmat ? 0 : g->coeff, .dmat = dmat, .out_tex = g->tex, .grid = &g->grid,
        .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = algo,
    };
    md_gpu_stream_sync(g->s);
    const md_tick_t t0 = md_tick_now();
    md_gto_gpu_density_launch(g->s, &desc);
    md_gpu_stream_sync(g->s);
    const double t = md_tick_to_seconds(md_tick_now() - t0);
    md_gpu_copy_from_texture(g->s, g->rb.gpu, g->tex, NULL);
    md_gpu_stream_sync(g->s);
    memcpy(out, g->rb.cpu, sizeof(float) * g->num_vox);
    return t;
}

static double dm_rel_diff(const float* a, const float* ref, size_t n, double* out_ref_max) {
    double d = 0.0, m = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (!(a[i] == a[i]) || !(ref[i] == ref[i])) return DBL_MAX;
        d = fmax(d, fabs((double)a[i] - (double)ref[i]));
        m = fmax(m, fabs((double)ref[i]));
    }
    if (out_ref_max) *out_ref_max = m;
    return m > 0.0 ? d / m : d;
}

// A grid of dim^3 voxels of h bohr around the molecule's centre
static md_grid_t dm_grid(const dm_src_t* s, int dim, float h) {
    double c[3] = {0, 0, 0};
    for (size_t a = 0; a < s->num_atoms; ++a) for (int k = 0; k < 3; ++k) c[k] += s->xyz[3 * a + k];
    for (int k = 0; k < 3; ++k) c[k] /= (double)s->num_atoms;
    const float half = 0.5f * h * (float)dim;
    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin  = vec3_set((float)c[0] - half, (float)c[1] - half, (float)c[2] - half),
        .spacing = vec3_set(h, h, h),
        .dim     = { dim, dim, dim },
    };
    return grid;
}

UTEST(gto_dmat, density_matches_packed) {
    static const char* files[] = {"h2o.h5", "mol.h5"};
    for (int fi = 0; fi < 2; ++fi) {
        dm_src_t s;
        ASSERT_TRUE(dm_load(&s, files[fi]));
        const md_grid_t grid = dm_grid(&s, 37, 0.3f);
        dm_gpu_t g;
        if (!dm_gpu_init(&g, &s, &grid)) { dm_gpu_free(&g); dm_free(&s); UTEST_SKIP("no GPU device"); }
        md_gto_gpu_dmat_t dm = md_gto_gpu_dmat_create(g.s, &(md_gto_gpu_dmat_desc_t){
            .basis = &s.basis, .atom_xyz = s.xyz, .cutoff = 1.0e-6 });
        ASSERT_TRUE(dm != NULL);
        md_gto_gpu_dmat_info_t info;
        md_gto_gpu_dmat_info(&info, dm);
        float* ref = (float*)malloc(sizeof(float) * g.num_vox);
        float* out = (float*)malloc(sizeof(float) * g.num_vox);
        for (int algo = MD_GTO_GPU_DENSITY_ALGO_REFERENCE; algo < MD_GTO_GPU_DENSITY_ALGO_COUNT; ++algo) {
            dm_eval(&g, NULL, (md_gto_gpu_density_algo_t)algo, ref);
            double ref_max = 0.0;
            for (int from = 0; from < 2; ++from) {
                ASSERT_TRUE(from == 0 ? md_gto_gpu_dmat_build(g.s, dm, s.rows, s.w, s.num_mo) : md_gto_gpu_dmat_set_matrix(g.s, dm, s.D));
                dm_eval(&g, dm, (md_gto_gpu_density_algo_t)algo, out);
                const double rel = dm_rel_diff(out, ref, g.num_vox, &ref_max);
                printf("  %s, %zu AOs, %zu of %zu tile pairs, algorithm %d, from the %s: max |diff| / max |rho| = %.2e (max %.3g)\n",
                       files[fi], info.num_ao, info.num_tiles, info.num_tiles_total, algo, from == 0 ? "orbitals" : "matrix  ", rel, ref_max);
                EXPECT_LT(rel, 1.0e-5);
                EXPECT_GT(ref_max, 1.0);
            }
        }
        free(out); free(ref);
        md_gto_gpu_dmat_destroy(g.s, dm);
        dm_gpu_free(&g);
        dm_free(&s);
    }
}

UTEST(gto_dmat, density_job_matches_launch) {
    dm_src_t s;
    ASSERT_TRUE(dm_load(&s, "mol.h5"));
    const md_grid_t grid = dm_grid(&s, 37, 0.3f);
    dm_gpu_t g;
    if (!dm_gpu_init(&g, &s, &grid)) { dm_gpu_free(&g); dm_free(&s); UTEST_SKIP("no GPU device"); }
    md_gto_gpu_dmat_t dm = md_gto_gpu_dmat_create(g.s, &(md_gto_gpu_dmat_desc_t){ .basis = &s.basis, .atom_xyz = s.xyz, .cutoff = 1.0e-6 });
    ASSERT_TRUE(dm != NULL);
    ASSERT_TRUE(md_gto_gpu_dmat_build(g.s, dm, s.rows, s.w, s.num_mo));
    float* ref = (float*)malloc(sizeof(float) * g.num_vox);
    float* out = (float*)malloc(sizeof(float) * g.num_vox);
    for (int algo = MD_GTO_GPU_DENSITY_ALGO_REFERENCE; algo < MD_GTO_GPU_DENSITY_ALGO_COUNT; ++algo) {
        for (int with_dmat = 0; with_dmat < 2; ++with_dmat) {
            dm_eval(&g, with_dmat ? dm : NULL, (md_gto_gpu_density_algo_t)algo, ref);
            // the volume cleared (so that a job that wrote nothing would show), then the job, polled without
            // waiting until it is issued
            md_gpu_addr_t zero = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, sizeof(float) * g.num_vox).gpu;
            ASSERT_TRUE(zero != 0);
            md_gpu_memset(g.s, zero, 0, sizeof(float) * g.num_vox);
            ASSERT_TRUE(md_gpu_copy_to_texture(g.s, g.tex, NULL, zero));
            md_gpu_free(g.s, zero);
            md_gto_gpu_density_desc_t desc = {
                .basis = g.basis, .atom_xyz = g.atoms, .coeff = with_dmat ? 0 : g.coeff, .dmat = with_dmat ? dm : NULL, .out_tex = g.tex,
                .grid = &g.grid, .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = (md_gto_gpu_density_algo_t)algo,
            };
            md_gto_gpu_job_t* job = md_gto_gpu_density_job_begin(g.s, &desc);
            ASSERT_TRUE(job != NULL);
            int polls = 0;
            while (!md_gto_gpu_job_poll(job, false)) { ++polls; md_thread_sleep(1); }
            md_gto_gpu_job_end(job);
            md_gpu_copy_from_texture(g.s, g.rb.gpu, g.tex, NULL);
            md_gpu_stream_sync(g.s);
            memcpy(out, g.rb.cpu, sizeof(float) * g.num_vox);
            size_t differ = 0;
            for (size_t i = 0; i < g.num_vox; ++i) differ += memcmp(&out[i], &ref[i], sizeof(float)) != 0;
            printf("  algorithm %d, %s: issued after %d polls, %zu of %zu voxels differ from the launch\n", algo, with_dmat ? "dmat  " : "packed", polls, differ, g.num_vox);
            EXPECT_EQ(differ, (size_t)0);
        }
    }
    // ended before it is issued: the volume keeps what it held (the last evaluation above)
    {
        md_gto_gpu_density_desc_t desc = {
            .basis = g.basis, .atom_xyz = g.atoms, .dmat = dm, .out_tex = g.tex, .grid = &g.grid,
            .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_ADD, .algorithm = MD_GTO_GPU_DENSITY_ALGO_GEMM,
        };
        md_gto_gpu_job_t* job = md_gto_gpu_density_job_begin(g.s, &desc);
        ASSERT_TRUE(job != NULL);
        md_gto_gpu_job_end(job);
        md_gpu_copy_from_texture(g.s, g.rb.gpu, g.tex, NULL);
        md_gpu_stream_sync(g.s);
        size_t differ = 0;
        for (size_t i = 0; i < g.num_vox; ++i) differ += memcmp(&((const float*)g.rb.cpu)[i], &out[i], sizeof(float)) != 0;
        EXPECT_EQ(differ, (size_t)0);
    }
    free(out); free(ref);
    md_gto_gpu_dmat_destroy(g.s, dm);
    dm_gpu_free(&g);
    dm_free(&s);
}

UTEST(gto_dmat, pattern_follows_motion) {
    dm_src_t s;
    ASSERT_TRUE(dm_load(&s, "mol.h5"));
    md_gpu_device_t dev = md_gpu_device_create(NULL);
    if (!dev) { dm_free(&s); UTEST_SKIP("no GPU device"); }
    md_gpu_stream_t st = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_dmat_t dm = md_gto_gpu_dmat_create(st, &(md_gto_gpu_dmat_desc_t){
        .basis = &s.basis, .atom_xyz = s.xyz, .cutoff = 1.0e-6, .margin = 1.0 });
    ASSERT_TRUE(dm != NULL);
    float* moved = (float*)malloc(sizeof(float) * 3 * s.num_atoms);
    memcpy(moved, s.xyz, sizeof(float) * 3 * s.num_atoms);
    EXPECT_TRUE(md_gto_gpu_dmat_pattern_valid(dm, moved, 0));
    moved[3 * (s.num_atoms / 2) + 1] += 0.45f;         // within margin / 2
    EXPECT_TRUE(md_gto_gpu_dmat_pattern_valid(dm, moved, 0));
    moved[3 * (s.num_atoms / 2) + 1] += 0.1f;          // beyond
    EXPECT_FALSE(md_gto_gpu_dmat_pattern_valid(dm, moved, 0));
    free(moved);
    md_gto_gpu_dmat_destroy(st, dm);
    md_gpu_stream_sync(st);
    md_gpu_device_destroy(dev);
    dm_free(&s);
}

UTEST(gto_dmat, benchmark) {
    const char* file = getenv("GTO_DMAT_FILE");
    if (!file || !file[0]) UTEST_SKIP("set GTO_DMAT_FILE=<file in test_data/vlx>");
    const int reps = getenv("GTO_DMAT_REPS") ? atoi(getenv("GTO_DMAT_REPS")) : 5;
    const double margin = getenv("GTO_DMAT_MARGIN") ? atof(getenv("GTO_DMAT_MARGIN")) : 0.0;
    const int gdim = getenv("GTO_DMAT_GRID") ? atoi(getenv("GTO_DMAT_GRID")) : 32;
    dm_src_t s;
    md_tick_t t0 = md_tick_now();
    ASSERT_TRUE(dm_load(&s, file));
    const double t_load = md_tick_to_seconds(md_tick_now() - t0);
    const md_grid_t grid = dm_grid(&s, gdim, 0.4f);
    dm_gpu_t g;
    if (!dm_gpu_init(&g, &s, &grid)) { dm_gpu_free(&g); dm_free(&s); UTEST_SKIP("no GPU device"); }
    md_gpu_device_info_t dinfo;
    md_gpu_device_info(g.dev, &dinfo);
    printf("\n  device: %s | %s: %zu atoms, %zu AOs, %zu occupied orbitals (loaded, and D in double on the CPU, in %.1f s) | margin %g bohr\n",
           dinfo.name, file, s.num_atoms, s.num_ao, s.num_mo, t_load, margin);

    t0 = md_tick_now();
    md_gto_gpu_dmat_t dm = md_gto_gpu_dmat_create(g.s, &(md_gto_gpu_dmat_desc_t){
        .basis = &s.basis, .atom_xyz = s.xyz, .cutoff = 1.0e-6, .margin = margin });
    md_gpu_stream_sync(g.s);
    const double t_create = md_tick_to_seconds(md_tick_now() - t0);
    ASSERT_TRUE(dm != NULL);
    md_gto_gpu_dmat_info_t info;
    md_gto_gpu_dmat_info(&info, dm);
    printf("  md_gto_gpu_dmat %s\n", info.revision);
    printf("  pattern: %zu of %zu tile pairs (%.1f%%), %.1f MB on the device (the packed triangle: %.1f MB); made in %.3f s\n",
           info.num_tiles, info.num_tiles_total, 100.0 * (double)info.num_tiles / (double)info.num_tiles_total, (double)info.bytes / 1e6,
           (double)md_gto_gpu_coeff_size_density(s.num_ao) / 1e6, t_create);

    // the build per frame: host (packing and upload) and the kernel, separately
    double tb[64], tk[64];
    const int nr = reps < 1 ? 1 : (reps > 64 ? 64 : reps);
    ASSERT_TRUE(md_gto_gpu_dmat_build(g.s, dm, s.rows, s.w, s.num_mo));   // warm up
    md_gpu_stream_sync(g.s);
    for (int i = 0; i < nr; ++i) {
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_gpu_dmat_build(g.s, dm, s.rows, s.w, s.num_mo));
        const md_tick_t t1 = md_tick_now();
        md_gpu_stream_sync(g.s);
        tb[i] = md_tick_to_seconds(md_tick_now() - t0);
        tk[i] = md_tick_to_seconds(md_tick_now() - t1);
    }
    for (int i = 1; i < nr; ++i) for (int j = i; j > 0 && tb[j] < tb[j - 1]; --j) { double x = tb[j]; tb[j] = tb[j - 1]; tb[j - 1] = x; x = tk[j]; tk[j] = tk[j - 1]; tk[j - 1] = x; }
    const double fma = (double)info.num_tiles * 256.0 * (double)s.num_mo;
    printf("  build from the orbitals (median of %d): %.1f ms, of it %.1f ms after the host returned (the kernel and the upload in flight);"
           " %.3g multiply-adds, %.1f G/s overall; %.1f MB of orbitals\n",
           nr, tb[nr / 2] * 1e3, tk[nr / 2] * 1e3, fma, fma / tb[nr / 2] / 1e9, 4.0 * (double)s.num_mo * (double)(16 * ((s.num_ao + 15) / 16)) / 1e6);
    t0 = md_tick_now();
    ASSERT_TRUE(md_gto_gpu_dmat_set_matrix(g.s, dm, s.D));
    md_gpu_stream_sync(g.s);
    printf("  fill from a matrix on the host: %.1f ms\n", md_tick_to_seconds(md_tick_now() - t0) * 1e3);

    // a sub-volume around the centre, the packed matrix against the dmat (from the orbitals)
    ASSERT_TRUE(md_gto_gpu_dmat_build(g.s, dm, s.rows, s.w, s.num_mo));
    float* ref = (float*)malloc(sizeof(float) * g.num_vox);
    float* out = (float*)malloc(sizeof(float) * g.num_vox);
    double tp = dm_eval(&g, NULL, MD_GTO_GPU_DENSITY_ALGO_GEMM, ref);
    tp = fmin(tp, dm_eval(&g, NULL, MD_GTO_GPU_DENSITY_ALGO_GEMM, ref));
    double td = dm_eval(&g, dm, MD_GTO_GPU_DENSITY_ALGO_GEMM, out);
    td = fmin(td, dm_eval(&g, dm, MD_GTO_GPU_DENSITY_ALGO_GEMM, out));
    double ref_max = 0.0;
    const double rel = dm_rel_diff(out, ref, g.num_vox, &ref_max);
    printf("  %d^3 voxels of 0.4 bohr (GEMM): packed %.1f ms, dmat %.1f ms; max |diff| / max |rho| = %.2e (max %.3g)\n",
           gdim, tp * 1e3, td * 1e3, rel, ref_max);
    EXPECT_LT(rel, 1.0e-4);
    free(out); free(ref);
    md_gto_gpu_dmat_destroy(g.s, dm);
    dm_gpu_free(&g);
    dm_free(&s);
}

#endif // MD_ENABLE_GPU
