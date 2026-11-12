// md_gto_gpu_dmat: a density matrix in tiles of the AO pairs that meet, built on the GPU from orbitals or
// filled from a matrix, and evaluated by md_gto_gpu_density_launch.
//
//   gto_dmat.density_matches_packed   h2o and mol: the volume through a dmat (from the orbitals, and from
//                                     the matrix) against the packed matrix, every algorithm
//   gto_dmat.pattern_follows_motion   the pattern holds within margin / 2 and not beyond
//   gto_dmat.density_job_matches_launch  mol: md_gto_gpu_density_job_begin, polled, bitwise the launch (the dmat
//                                     and the packed matrix, every algorithm); a job ended before it is
//                                     issued writes nothing
//   gto_localization.matches_cpu      ELF and LOL by their grid kernel (md_gto_gpu_localization_launch) against
//                                     md_qm_elf_evaluate in double at the voxels: h2o and mol (closed shells; the
//                                     packed matrix and a dmat), the OH radical (UHF, open shell)
//   gto_gpu.plane_matches_volume      a grid of one layer (the 2D kernels, 32x16x1) against the first layer of the
//                                     same grid 8 deep (3D, 8x8x8): density, orbitals, ELF and LOL
//   gto_localization.matches_independent  ELF and LOL against PySCF's AOs with Multiwfn's formulas
//                                     (unittest/elf_ref_data, tools/elf_reference.py), out into the tails
//   gto_localization.pieces_match_one_launch  mol: ELF and LOL in pieces of a group each (piece_work), as a launch
//                                     and as a job polled over many polls, bitwise the launch in one piece; the
//                                     job's progress reaches 1
//   gto_localization.tails            out to 9 bohr beyond the atoms, no floor of rho: ELF with d_min (Multiwfn's
//                                     convention) and without it (within its noise, md_qm_elf_noise), LOL, against
//                                     the matrices in double
//   gto_dmat.density_pieces_match_launch  mol: the GEMM in pieces of a few blocks each (piece_work), as a launch
//                                     and as a job polled over many polls, bitwise the launch in one piece; the
//                                     job's progress reaches 1, and a job dropped midway lets the stream go on
//   gto_dmat.benchmark                GTO_DMAT_FILE=<file in test_data/vlx> (protein.h5): the pattern, the
//                                     build per frame (packing, upload, kernel), and a sub-volume through
//                                     the dmat against the packed matrix; GTO_DMAT_REPS (5),
//                                     GTO_DMAT_MARGIN (0), GTO_DMAT_GRID (32: voxels per edge, 0.4 bohr)

#include "utest.h"
#include "vlx_test_util.h"

#include <md_gto.h>
#include <md_qm.h>
#include <md_molden.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

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
    if (!vlx_test_load(&s->t, md_str_from_cstr(path), MD_MEGABYTES(1024))) return false;
    if (!qm_test_basis(&s->basis, &s->t)) return false;
    s->num_atoms = s->t.sys.atom.count;
    s->xyz = (float*)malloc(sizeof(float) * 3 * s->num_atoms);
    if (vlx_test_atom_xyz_bohr(s->xyz, 3 * s->num_atoms, &s->t) != s->num_atoms) return false;
    s->num_mo = md_qm_density_orbitals(&s->rows, &s->w, &s->num_ao, &s->t.sys, MD_STR_LIT("orbital/total/density"), md_get_heap_allocator());
    if (s->num_mo == 0 || s->num_ao != md_gto_basis_num_ao(&s->basis)) return false;
    s->D = (double*)malloc(sizeof(double) * s->num_ao * s->num_ao);
    return s->D && md_gto_density_from_orbitals(s->D, s->rows, s->num_ao, 1, s->w, s->num_mo, s->num_ao, NULL);
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
        .orientation = md_float3x3_ident(),
        .origin  = md_float3_set((float)c[0] - half, (float)c[1] - half, (float)c[2] - half),
        .spacing = md_float3_set(h, h, h),
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
            // Bitwise, every algorithm: the reference kernel lists a workgroup's screened CGTOs in the order of their
            // index (gto_screen.slang), no longer in the order its waves win an atomic (Metal: the last bits run to run)
            size_t differ = 0;
            for (size_t i = 0; i < g.num_vox; ++i) differ += memcmp(&out[i], &ref[i], sizeof(float)) != 0;
            double ref_max = 0.0;
            const double rel = dm_rel_diff(out, ref, g.num_vox, &ref_max);
            printf("  algorithm %d, %s: issued after %d polls, %zu of %zu voxels differ from the launch, max |diff| / max |rho| = %.2e\n",
                   algo, with_dmat ? "dmat  " : "packed", polls, differ, g.num_vox, rel);
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

UTEST(gto_dmat, density_pieces_match_launch) {
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
    // a piece of about two blocks of 64 AOs: hundreds of pieces, batches split by their blocks
    const uint64_t piece_work = (uint64_t)2 * 512 * 64 * 64;
    for (int with_dmat = 0; with_dmat < 2; ++with_dmat) {
        dm_eval(&g, with_dmat ? dm : NULL, MD_GTO_GPU_DENSITY_ALGO_GEMM, ref);
        md_gto_gpu_density_desc_t desc = {
            .basis = g.basis, .atom_xyz = g.atoms, .coeff = with_dmat ? 0 : g.coeff, .dmat = with_dmat ? dm : NULL, .out_tex = g.tex,
            .grid = &g.grid, .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = MD_GTO_GPU_DENSITY_ALGO_GEMM,
            .piece_work = piece_work,
        };
        for (int as_job = 0; as_job < 2; ++as_job) {
            md_gpu_addr_t zero = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, sizeof(float) * g.num_vox).gpu;
            ASSERT_TRUE(zero != 0);
            md_gpu_memset(g.s, zero, 0, sizeof(float) * g.num_vox);
            ASSERT_TRUE(md_gpu_copy_to_texture(g.s, g.tex, NULL, zero));
            md_gpu_free(g.s, zero);
            int polls = 0;
            float progress = 0.0f;
            if (as_job) {
                md_gto_gpu_job_t* job = md_gto_gpu_density_job_begin(g.s, &desc);
                ASSERT_TRUE(job != NULL);
                float last = 0.0f;
                while (!md_gto_gpu_job_poll(job, false)) {
                    const float f = md_gto_gpu_job_progress(job);
                    EXPECT_GE(f, last);     // never back
                    last = f;
                    ++polls;
                }
                md_gpu_stream_sync(g.s);
                progress = md_gto_gpu_job_progress(job);
                md_gto_gpu_job_end(job);
            } else {
                md_gto_gpu_density_launch(g.s, &desc);
            }
            md_gpu_copy_from_texture(g.s, g.rb.gpu, g.tex, NULL);
            md_gpu_stream_sync(g.s);
            memcpy(out, g.rb.cpu, sizeof(float) * g.num_vox);
            size_t differ = 0;
            for (size_t i = 0; i < g.num_vox; ++i) differ += memcmp(&out[i], &ref[i], sizeof(float)) != 0;
            printf("  %s, %s: %d polls, progress %.3f, %zu of %zu voxels differ from the launch in one piece\n",
                   with_dmat ? "dmat  " : "packed", as_job ? "job   " : "launch", polls, progress, differ, g.num_vox);
            EXPECT_EQ(differ, (size_t)0);
            if (as_job) {
                EXPECT_GT(polls, 0);    // more pieces than are kept in flight
                EXPECT_EQ(progress, 1.0f);
            }
        }
    }
    // dropped midway: what is in flight completes, nothing more is issued, and the stream goes on
    {
        md_gto_gpu_density_desc_t desc = {
            .basis = g.basis, .atom_xyz = g.atoms, .dmat = dm, .out_tex = g.tex, .grid = &g.grid,
            .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET, .algorithm = MD_GTO_GPU_DENSITY_ALGO_GEMM, .piece_work = piece_work,
        };
        md_gto_gpu_job_t* job = md_gto_gpu_density_job_begin(g.s, &desc);
        ASSERT_TRUE(job != NULL);
        md_gpu_stream_sync(g.s);
        EXPECT_FALSE(md_gto_gpu_job_poll(job, false));     // the counts planned, the first pieces issued
        EXPECT_LT(md_gto_gpu_job_progress(job), 1.0f);
        md_gto_gpu_job_end(job);
        dm_eval(&g, dm, MD_GTO_GPU_DENSITY_ALGO_GEMM, out);
        size_t differ = 0;
        for (size_t i = 0; i < g.num_vox; ++i) differ += memcmp(&out[i], &ref[i], sizeof(float)) != 0;
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


// ---------------------------------------------------------------------------
// ELF and LOL on a grid
// ---------------------------------------------------------------------------

typedef struct loc_src_t {
    vlx_test_t     t;
    md_gto_basis_t basis;
    float*         xyz;             // bohr, per atom
    size_t         num_atoms, num_ao;
    int            num_spins;       // 1: D[0] the total density of a closed shell; 2: alpha and beta
    double*        D[2];
} loc_src_t;

static bool loc_load(loc_src_t* s, const char* file) {
    memset(s, 0, sizeof(*s));
    char path[512];
    const bool molden = strstr(file, ".molden") != NULL;
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/%s/%s", molden ? "molden" : "vlx", file);
    if (molden) {
        qm_test_init(&s->t, MD_MEGABYTES(256));
        if (!md_molden_system_init_from_file(&s->t.sys, &s->t.state, md_str_from_cstr(path))) return false;
    } else if (!vlx_test_load(&s->t, md_str_from_cstr(path), MD_MEGABYTES(1024))) {
        return false;
    }
    if (!qm_test_basis(&s->basis, &s->t)) return false;
    s->num_atoms = s->t.sys.atom.count;
    s->num_ao = md_gto_basis_num_ao(&s->basis);
    s->xyz = (float*)malloc(sizeof(float) * 3 * s->num_atoms);
    if (vlx_test_atom_xyz_bohr(s->xyz, 3 * s->num_atoms, &s->t) != s->num_atoms) return false;
    // a closed shell as md_qm_elf_evaluate tells it: the same orbitals and occupations in both spins
    md_allocator_i* heap = md_get_heap_allocator();
    double *ar = NULL, *aw = NULL, *br = NULL, *bw = NULL;
    size_t ana = 0, bna = 0;
    const size_t an = md_qm_density_orbitals(&ar, &aw, &ana, &s->t.sys, MD_STR_LIT("orbital/alpha/density"), heap);
    const size_t bn = md_qm_density_orbitals(&br, &bw, &bna, &s->t.sys, MD_STR_LIT("orbital/beta/density"), heap);
    const bool closed = an == bn && ana == bna && memcmp(ar, br, sizeof(double) * an * ana) == 0 && memcmp(aw, bw, sizeof(double) * an) == 0;
    bool ok = an > 0 && bn > 0 && ana == s->num_ao;
    s->num_spins = closed ? 1 : 2;
    for (int k = 0; ok && k < s->num_spins; ++k) {
        s->D[k] = (double*)calloc(s->num_ao * s->num_ao, sizeof(double));
        const double* rows = k == 0 ? ar : br;
        const double* w = k == 0 ? aw : bw;
        ok = s->D[k] && md_gto_density_from_orbitals(s->D[k], rows, s->num_ao, 1, w, k == 0 ? an : bn, s->num_ao, NULL);
        if (ok && closed) for (size_t i = 0; i < s->num_ao * s->num_ao; ++i) s->D[k][i] *= 2.0;   // the total
    }
    if (ar) md_free(heap, ar, sizeof(double) * an * ana);
    if (aw) md_free(heap, aw, sizeof(double) * an);
    if (br) md_free(heap, br, sizeof(double) * bn * bna);
    if (bw) md_free(heap, bw, sizeof(double) * bn);
    return ok;
}

static void loc_free(loc_src_t* s) {
    free(s->D[0]);
    free(s->D[1]);
    free(s->xyz);
    qm_test_free(&s->t);
}

// ELF (kind ELF) or LOL into out at the grid's voxel centres, by the kernel, from the packed matrices or dmats
static double loc_d_min = 0.0;   // the d_min the kernel is launched with (md_qm_elf; 0: none)
static uint64_t loc_piece_work = 0;  // the work of a submission (0: the default)
static bool loc_as_job = false;      // as a job (md_gto_gpu_localization_job_begin), polled without waiting
static int loc_job_polls = 0;
static float loc_job_progress = 0.0f;

static bool loc_eval_gpu(float* out, const loc_src_t* s, const md_grid_t* grid, md_gto_gpu_localization_t kind, bool use_dmat, double rho_min, double* out_seconds) {
    md_gpu_device_t dev = md_gpu_device_create(NULL);
    if (!dev) return false;
    md_gpu_stream_t st = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_initialize(dev);
    const double cutoff = 1.0e-8;
    md_gto_gpu_basis_t basis = md_gto_gpu_basis_create(st, &(md_gto_gpu_basis_desc_t){ .basis = &s->basis, .cutoff = cutoff });
    md_gpu_addr_t atoms = md_gpu_malloc(st, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(s->num_atoms)).gpu;
    float* p = (float*)md_gpu_upload_begin(st, atoms, md_gto_gpu_atom_buffer_size(s->num_atoms));
    md_gto_gpu_atom_pack(p, s->xyz, 0, s->num_atoms);
    md_gpu_upload_end(st);
    md_gto_gpu_localization_desc_t desc = {
        .basis = basis, .atom_xyz = atoms, .num_spins = (uint32_t)s->num_spins, .kind = kind, .rho_min = rho_min,
        .d_min = loc_d_min, .grid = grid, .sample_offset = {0.5f, 0.5f, 0.5f}, .piece_work = loc_piece_work,
    };
    md_gpu_addr_t coeff[2] = {0, 0};
    for (int k = 0; k < s->num_spins; ++k) {
        if (use_dmat) {
            desc.dmat[k] = md_gto_gpu_dmat_create(st, &(md_gto_gpu_dmat_desc_t){ .basis = &s->basis, .atom_xyz = s->xyz, .cutoff = cutoff });
            md_gto_gpu_dmat_set_matrix(st, desc.dmat[k], s->D[k]);
        } else {
            coeff[k] = md_gpu_malloc(st, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_density(s->num_ao)).gpu;
            p = (float*)md_gpu_upload_begin(st, coeff[k], md_gto_gpu_coeff_size_density(s->num_ao));
            md_gto_gpu_coeff_pack_density(p, s->D[k], s->num_ao);
            md_gpu_upload_end(st);
            desc.coeff[k] = coeff[k];
        }
    }
    const size_t n = (size_t)grid->dim[0] * grid->dim[1] * grid->dim[2];
    desc.out_tex = md_gpu_texture_create(st, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)grid->dim[0], .height = (uint32_t)grid->dim[1], .depth_or_layers = (uint32_t)grid->dim[2],
    });
    md_gpu_mem_t rb = md_gpu_malloc(st, MD_GPU_MEM_HOST_READ, sizeof(float) * n);
    // timed: a second launch, the kernel made by the first
    if (out_seconds && basis && atoms && desc.out_tex) md_gto_gpu_localization_launch(st, &desc);
    md_gpu_stream_sync(st);
    const md_tick_t t0 = md_tick_now();
    bool ok = basis && atoms && desc.out_tex && rb.cpu;
    if (ok && loc_as_job) {
        md_gto_gpu_job_t* job = md_gto_gpu_localization_job_begin(st, &desc);
        ok = job != NULL;
        loc_job_polls = 0;
        while (ok && !md_gto_gpu_job_poll(job, false)) loc_job_polls++;
        md_gpu_stream_sync(st);
        loc_job_progress = md_gto_gpu_job_progress(job);
        md_gto_gpu_job_end(job);
    } else if (ok) {
        ok = md_gto_gpu_localization_launch(st, &desc);
    }
    md_gpu_stream_sync(st);
    if (out_seconds) *out_seconds = md_tick_to_seconds(md_tick_now() - t0);
    if (ok) {
        md_gpu_copy_from_texture(st, rb.gpu, desc.out_tex, NULL);
        md_gpu_stream_sync(st);
        memcpy(out, rb.cpu, sizeof(float) * n);
    }
    md_gpu_free(st, rb.gpu);
    if (desc.out_tex) md_gpu_texture_destroy(desc.out_tex);
    for (int k = 0; k < 2; ++k) {
        if (coeff[k]) md_gpu_free(st, coeff[k]);
        if (desc.dmat[k]) md_gto_gpu_dmat_destroy(st, desc.dmat[k]);
    }
    md_gpu_free(st, atoms);
    if (basis) md_gto_gpu_basis_destroy(st, basis);
    md_gpu_stream_sync(st);
    md_gto_gpu_shutdown();
    md_gpu_device_destroy(dev);
    return ok;
}

UTEST(gto_localization, pieces_match_one_launch) {
    loc_src_t s;
    if (!loc_load(&s, "mol.h5")) { loc_free(&s); ASSERT_TRUE(false); }
    double c[3] = {0, 0, 0};
    for (size_t a = 0; a < s.num_atoms; ++a) for (int k = 0; k < 3; ++k) c[k] += s.xyz[3 * a + k] / (double)s.num_atoms;
    const int dim = 28;
    const float h = 0.35f;
    const md_grid_t grid = {
        .orientation = md_float3x3_ident(),
        .origin  = md_float3_set((float)c[0] - 0.5f * h * dim, (float)c[1] - 0.5f * h * dim, (float)c[2] - 0.5f * h * dim),
        .spacing = md_float3_set(h, h, h),
        .dim     = { dim, dim, dim },
    };
    const size_t n = (size_t)dim * dim * dim;
    float* ref = (float*)malloc(sizeof(float) * n);
    float* out = (float*)malloc(sizeof(float) * n);
    for (int kind = 0; kind < 2; ++kind) {
        const md_gto_gpu_localization_t k = kind ? MD_GTO_GPU_LOL : MD_GTO_GPU_ELF;
        loc_piece_work = (uint64_t)1 << 62;     // one piece
        loc_as_job = false;
        if (!loc_eval_gpu(ref, &s, &grid, k, true, 0.0, NULL)) { free(ref); free(out); loc_free(&s); UTEST_SKIP("no GPU device"); }
        for (int as_job = 0; as_job < 2; ++as_job) {
            loc_piece_work = 1;                 // a group each: 4^3 pieces
            loc_as_job = as_job != 0;
            ASSERT_TRUE(loc_eval_gpu(out, &s, &grid, k, true, 0.0, NULL));
            size_t differ = 0;
            for (size_t i = 0; i < n; ++i) differ += memcmp(&out[i], &ref[i], sizeof(float)) != 0;
            printf("  %s %s: %zu of %zu voxels differ from one launch (" MD_GTO_GPU_LOCALIZATION_REVISION ")\n", kind ? "LOL" : "ELF",
                   as_job ? "job   " : "launch", differ, n);
            if (as_job) printf("    %d polls, progress %.3f\n", loc_job_polls, loc_job_progress);
            EXPECT_EQ(differ, (size_t)0);
            if (as_job) {
                EXPECT_GT(loc_job_polls, 0);
                EXPECT_EQ(loc_job_progress, 1.0f);
            }
        }
    }
    loc_piece_work = 0;
    loc_as_job = false;
    free(ref); free(out);
    loc_free(&s);
}

UTEST(gto_localization, matches_cpu) {
    static const char* files[] = {"h2o.h5", "mol.h5", "oh_631gs_uhf.molden"};
    const double rho_min = 1.0e-4;
    for (int fi = 0; fi < 3; ++fi) {
        loc_src_t s;
        if (!loc_load(&s, files[fi])) { loc_free(&s); ASSERT_TRUE(false); }
        // a cube of 28^3 voxels of 0.35 bohr about the molecule's centre
        loc_src_t* ps = &s;
        double c[3] = {0, 0, 0};
        for (size_t a = 0; a < ps->num_atoms; ++a) for (int k = 0; k < 3; ++k) c[k] += ps->xyz[3 * a + k] / (double)ps->num_atoms;
        const int dim = 28;
        const float h = 0.35f;
        const md_grid_t grid = {
            .orientation = md_float3x3_ident(),
            .origin  = md_float3_set((float)c[0] - 0.5f * h * dim, (float)c[1] - 0.5f * h * dim, (float)c[2] - 0.5f * h * dim),
            .spacing = md_float3_set(h, h, h),
            .dim     = { dim, dim, dim },
        };
        const size_t n = (size_t)dim * dim * dim;
        float* pts = (float*)malloc(sizeof(float) * 3 * n);
        for (int z = 0; z < dim; ++z) for (int y = 0; y < dim; ++y) for (int x = 0; x < dim; ++x) {
            const size_t i = ((size_t)z * dim + y) * dim + x;
            pts[3 * i + 0] = grid.origin.x + ((float)x + 0.5f) * h;
            pts[3 * i + 1] = grid.origin.y + ((float)y + 0.5f) * h;
            pts[3 * i + 2] = grid.origin.z + ((float)z + 0.5f) * h;
        }
        double* elf = (double*)malloc(sizeof(double) * n);
        double* lol = (double*)malloc(sizeof(double) * n);
        ASSERT_EQ(md_qm_elf_evaluate(elf, lol, pts, n, 0, &s.t.sys, s.xyz, 0, rho_min, 0.0), n);
        float* g = (float*)malloc(sizeof(float) * n);
        for (int kind = 0; kind < 2; ++kind) {
            for (int dmat = 0; dmat < 2; ++dmat) {
                double sec = 0.0;
                if (!loc_eval_gpu(g, &s, &grid, kind ? MD_GTO_GPU_LOL : MD_GTO_GPU_ELF, dmat != 0, rho_min, &sec)) {
                    free(g); free(elf); free(lol); free(pts); loc_free(&s);
                    UTEST_SKIP("no GPU device");
                }
                const double* ref = kind ? lol : elf;
                double worst = 0.0, sum = 0.0;
                size_t nz = 0, edge = 0;
                for (size_t i = 0; i < n; ++i) {
                    if (ref[i] == 0.0 && g[i] == 0.0) continue;
                    // where one side is at rho_min and the other not (the floor itself, in float), not compared
                    if ((ref[i] == 0.0) != (g[i] == 0.0f)) { ++edge; continue; }
                    const double d = fabs((double)g[i] - ref[i]);
                    worst = fmax(worst, d);
                    sum += d;
                    ++nz;
                }
                printf("  %-20s %s %-6s: %zu voxels, max |gpu - cpu| %.2e, mean %.2e (%zu at the rho_min floor), kernel %.0f ms (" MD_GTO_GPU_LOCALIZATION_REVISION ")\n", files[fi],
                       kind ? "LOL" : "ELF", dmat ? "dmat" : "packed", nz, worst, nz ? sum / nz : 0.0, edge, 1000.0 * sec);
                // measured on lavapipe: ELF within 7.8e-5 (mol), LOL within 2.4e-6, the same through a dmat
                EXPECT_GT(nz, n / 20);
                EXPECT_LT(worst, kind ? 2.0e-5 : 5.0e-4);
                EXPECT_LT(edge, n / 1000 + 1);
            }
        }
        free(g); free(elf); free(lol); free(pts);
        loc_free(&s);
    }
}
// A grid of one layer is evaluated by the 2D kernels (blocks of 32x16x1 voxels, 16x16x1 for ELF and LOL): what they
// give is the first layer of the same grid eight layers deep, by the 3D ones (8x8x8), within the difference of what the
// two screen out. A turned grid at a plane's spacing, as a plane's is; the density (packed and dmat, every algorithm),
// orbitals (1: shell, 4: shell of 8, 12: GEMM) and ELF, LOL. With the times of one layer by the 2D kernels and by the
// 3D ones (MD_GTO_GPU_NO_2D), each a second launch: on llvmpipe 6x faster for the density's GEMM path and the
// orbitals, 1.3x for ELF and LOL, whose cost is the pairs of AOs at each voxel; the reference density has no 2D kernel.
typedef struct plane_pair_t {
    md_grid_t        g3, g2;
    md_gpu_texture_t t3, t2;
    md_gpu_mem_t     rb3, rb2;
    size_t           n3, n2;
    float            *v3, *v2;
} plane_pair_t;

static md_grid_t plane_grid(const float* xyz, size_t num_atoms, int layers) {
    double c[3] = {0, 0, 0};
    for (size_t a = 0; a < num_atoms; ++a) for (int k = 0; k < 3; ++k) c[k] += xyz[3 * a + k] / (double)num_atoms;
    // turned by 0.5 rad about (1, 2, 3) / |(1, 2, 3)|
    const double ax[3] = {1.0 / sqrt(14.0), 2.0 / sqrt(14.0), 3.0 / sqrt(14.0)}, an = 0.5, co = cos(an), si = sin(an);
    md_float3x3_t R;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        const double k[3][3] = {{0, -ax[2], ax[1]}, {ax[2], 0, -ax[0]}, {-ax[1], ax[0], 0}};
        double kk = 0; for (int m = 0; m < 3; ++m) kk += k[i][m] * k[m][j];
        R.elem[j][i] = (float)((i == j ? 1.0 : 0.0) + si * k[i][j] + (1.0 - co) * kk);   // column j, row i
    }
    // a plane's spacing (0.02 to 0.05 Angstrom: 0.04 to 0.09 bohr) over a few bohr: GTO_PLANE_H, _W, _H to measure others
    const float h = getenv("GTO_PLANE_H") ? (float)atof(getenv("GTO_PLANE_H")) : 0.07f;
    const int w  = getenv("GTO_PLANE_W") ? atoi(getenv("GTO_PLANE_W")) : 96;
    const int hh = getenv("GTO_PLANE_HH") ? atoi(getenv("GTO_PLANE_HH")) : 80;
    // the corner: the centre less half the extent along the first two axes, the first layer through the centre
    md_grid_t g = { .orientation = R, .spacing = md_float3_set(h, h, h), .dim = { w, hh, layers } };
    md_float3_t o = md_float3_set((float)c[0], (float)c[1], (float)c[2]);
    for (int k = 0; k < 3; ++k) o.elem[k] -= 0.5f * h * (float)w * R.elem[0][k] + 0.5f * h * (float)hh * R.elem[1][k];
    g.origin = o;
    return g;
}

static bool plane_pair_init(plane_pair_t* p, md_gpu_stream_t st, const float* xyz, size_t num_atoms) {
    memset(p, 0, sizeof(*p));
    p->g3 = plane_grid(xyz, num_atoms, 8);
    p->g2 = plane_grid(xyz, num_atoms, 1);
    p->n3 = (size_t)p->g3.dim[0] * p->g3.dim[1] * p->g3.dim[2];
    p->n2 = (size_t)p->g2.dim[0] * p->g2.dim[1];
    for (int k = 0; k < 2; ++k) {
        const md_grid_t* g = k ? &p->g2 : &p->g3;
        md_gpu_texture_t t = md_gpu_texture_create(st, &(md_gpu_texture_desc_t){
            .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
            .width = (uint32_t)g->dim[0], .height = (uint32_t)g->dim[1], .depth_or_layers = (uint32_t)g->dim[2] });
        if (k) p->t2 = t; else p->t3 = t;
    }
    p->rb3 = md_gpu_malloc(st, MD_GPU_MEM_HOST_READ, sizeof(float) * p->n3);
    p->rb2 = md_gpu_malloc(st, MD_GPU_MEM_HOST_READ, sizeof(float) * p->n2);
    p->v3 = (float*)malloc(sizeof(float) * p->n3);
    p->v2 = (float*)malloc(sizeof(float) * p->n2);
    return p->t3 && p->t2 && p->rb3.cpu && p->rb2.cpu;
}

static void plane_pair_read(plane_pair_t* p, md_gpu_stream_t st) {
    md_gpu_copy_from_texture(st, p->rb3.gpu, p->t3, NULL);
    md_gpu_copy_from_texture(st, p->rb2.gpu, p->t2, NULL);
    md_gpu_stream_sync(st);
    memcpy(p->v3, p->rb3.cpu, sizeof(float) * p->n3);
    memcpy(p->v2, p->rb2.cpu, sizeof(float) * p->n2);
}

// max |one layer - first layer| / max |first layer|
static double plane_pair_diff(const plane_pair_t* p, double* out_max) {
    return dm_rel_diff(p->v2, p->v3, p->n2, out_max);
}

static void plane_pair_free(plane_pair_t* p, md_gpu_stream_t st) {
    if (p->t3) md_gpu_texture_destroy(p->t3);
    if (p->t2) md_gpu_texture_destroy(p->t2);
    md_gpu_free(st, p->rb3.gpu);
    md_gpu_free(st, p->rb2.gpu);
    free(p->v3); free(p->v2);
}

UTEST(gto_gpu, plane_matches_volume) {
    printf("  (the kernels' times, for measurements: " MD_GTO_GPU_PLANE_REVISION ")\n");
    static const char* files[] = {"h2o.h5", "mol.h5"};
    for (int fi = 0; fi < 2; ++fi) {
        dm_src_t s;
        ASSERT_TRUE(dm_load(&s, files[fi]));
        const md_grid_t any = dm_grid(&s, 8, 0.3f);
        dm_gpu_t g;
        if (!dm_gpu_init(&g, &s, &any)) { dm_gpu_free(&g); dm_free(&s); UTEST_SKIP("no GPU device"); }
        plane_pair_t pp;
        ASSERT_TRUE(plane_pair_init(&pp, g.s, s.xyz, s.num_atoms));
        md_gto_gpu_dmat_t dm = md_gto_gpu_dmat_create(g.s, &(md_gto_gpu_dmat_desc_t){ .basis = &s.basis, .atom_xyz = s.xyz, .cutoff = 1.0e-6 });
        ASSERT_TRUE(dm && md_gto_gpu_dmat_set_matrix(g.s, dm, s.D));
        double m = 0.0, rel = 0.0;
        // the density
        for (int algo = MD_GTO_GPU_DENSITY_ALGO_REFERENCE; algo < MD_GTO_GPU_DENSITY_ALGO_COUNT; ++algo) {
            for (int use_dm = 0; use_dm < 2; ++use_dm) {
                double ms[3];
                for (int k = 0; k < 3; ++k) {
                    if (k == 2) md_os_env_set("MD_GTO_GPU_NO_2D", "1");   // one layer again, by the 3D kernels: what it cost before
                    md_gto_gpu_density_desc_t d = {
                        .basis = g.basis, .atom_xyz = g.atoms, .coeff = use_dm ? 0 : g.coeff, .dmat = use_dm ? dm : NULL,
                        .out_tex = k ? pp.t2 : pp.t3, .grid = k ? &pp.g2 : &pp.g3, .sample_offset = {0.5f, 0.5f, 0.5f},
                        .op = MD_GTO_OP_SET, .algorithm = (md_gto_gpu_density_algo_t)algo,
                    };
                    if (k == 2) plane_pair_read(&pp, g.s);   // (what is compared: the 2D kernels')
                    md_gto_gpu_density_launch(g.s, &d);      // once to make its kernels, then timed
                    md_gpu_stream_sync(g.s);
                    const md_tick_t t0 = md_tick_now();
                    md_gto_gpu_density_launch(g.s, &d);
                    md_gpu_stream_sync(g.s);
                    ms[k] = md_tick_to_milliseconds(md_tick_now() - t0);
                    if (k == 2) md_os_env_set("MD_GTO_GPU_NO_2D", NULL);
                }
                rel = plane_pair_diff(&pp, &m);
                printf("  %-8s density, algorithm %d, %s: max |plane - layer| / max = %.2e (max %.3g); 1 layer %.1f ms (3D kernels %.1f ms)\n",
                       files[fi], algo, use_dm ? "dmat  " : "packed", rel, m, ms[1], ms[2]);
                EXPECT_LT(rel, 1.0e-5);
                EXPECT_GT(m, 1.0e-3);
            }
        }
        // orbitals: the first n occupied, as psi (1) or psi^2 (more)
        static const size_t counts[] = {1, 4, 12};
        for (int ci = 0; ci < 3; ++ci) {
            const size_t n = counts[ci] < s.num_mo ? counts[ci] : s.num_mo;
            const double** rows = (const double**)malloc(sizeof(double*) * n);
            for (size_t i = 0; i < n; ++i) rows[i] = s.rows + (s.num_mo - 1 - i) * s.num_ao;   // from the highest
            md_gpu_addr_t c = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_mo(n, s.num_ao)).gpu;
            float* dst = (float*)md_gpu_upload_begin(g.s, c, md_gto_gpu_coeff_size_mo(n, s.num_ao));
            md_gto_gpu_coeff_pack_mo(dst, rows, NULL, n, s.num_ao);
            md_gpu_upload_end(g.s);
            double ms[3];
            for (int k = 0; k < 3; ++k) {
                if (k == 2) { plane_pair_read(&pp, g.s); md_os_env_set("MD_GTO_GPU_NO_2D", "1"); }
                md_gto_gpu_orbital_desc_t d = {
                    .basis = g.basis, .atom_xyz = g.atoms, .coeff = c, .out_tex = k ? pp.t2 : pp.t3, .grid = k ? &pp.g2 : &pp.g3,
                    .sample_offset = {0.5f, 0.5f, 0.5f}, .num_orbitals = n,
                    .eval_mode = n == 1 ? MD_GTO_EVAL_MODE_PSI : MD_GTO_EVAL_MODE_PSI_SQUARED, .op = MD_GTO_OP_SET,
                };
                md_gto_gpu_orbital_launch(g.s, &d);     // once to make its kernels, then timed
                md_gpu_stream_sync(g.s);
                const md_tick_t t0 = md_tick_now();
                md_gto_gpu_orbital_launch(g.s, &d);
                md_gpu_stream_sync(g.s);
                ms[k] = md_tick_to_milliseconds(md_tick_now() - t0);
                if (k == 2) md_os_env_set("MD_GTO_GPU_NO_2D", NULL);
            }
            rel = plane_pair_diff(&pp, &m);
            printf("  %-8s %zu orbital%s: max |plane - layer| / max = %.2e (max %.3g); 1 layer %.1f ms (3D kernels %.1f ms)\n",
                   files[fi], n, n == 1 ? " " : "s", rel, m, ms[1], ms[2]);
            // the shell kernel screens by |C| |phi| per block: other blocks, other shells at the cutoff (the basis's, 1e-6):
            // a difference of that order, not relative to the values (2.2e-6 measured, mol, one orbital)
            EXPECT_LT(rel * m, 1.0e-5);
            EXPECT_GT(m, 1.0e-3);
            md_gpu_stream_sync(g.s);
            md_gpu_free(g.s, c);
            free(rows);
        }
        md_gto_gpu_dmat_destroy(g.s, dm);
        plane_pair_free(&pp, g.s);
        dm_gpu_free(&g);
        dm_free(&s);
    }
    // ELF and LOL (loc_eval_gpu makes a device of its own each time)
    loc_src_t ls;
    ASSERT_TRUE(loc_load(&ls, "mol.h5"));
    const md_grid_t g3 = plane_grid(ls.xyz, ls.num_atoms, 8), g2 = plane_grid(ls.xyz, ls.num_atoms, 1);
    const size_t n3 = (size_t)g3.dim[0] * g3.dim[1] * g3.dim[2], n2 = (size_t)g2.dim[0] * g2.dim[1];
    float* v3 = (float*)malloc(sizeof(float) * n3);
    float* v2 = (float*)malloc(sizeof(float) * n2);
    for (int kind = 0; kind < 2; ++kind) {
        const md_gto_gpu_localization_t kd = kind ? MD_GTO_GPU_LOL : MD_GTO_GPU_ELF;
        double ms[3] = {0, 0, 0};
        float* v2b = (float*)malloc(sizeof(float) * n2);
        md_os_env_set("MD_GTO_GPU_NO_2D", "1");
        const bool ok3 = loc_eval_gpu(v2b, &ls, &g2, kd, true, 1.0e-4, &ms[2]);
        md_os_env_set("MD_GTO_GPU_NO_2D", NULL);
        free(v2b);
        if (!ok3 || !loc_eval_gpu(v3, &ls, &g3, kd, true, 1.0e-4, &ms[0]) || !loc_eval_gpu(v2, &ls, &g2, kd, true, 1.0e-4, &ms[1])) {
            free(v3); free(v2); loc_free(&ls);
            UTEST_SKIP("no GPU device");
        }
        double d = 0.0;
        for (size_t i = 0; i < n2; ++i) d = fmax(d, fabs((double)v2[i] - (double)v3[i]));
        printf("  mol.h5   %s: max |plane - layer| = %.2e; 1 layer %.1f ms (3D kernels %.1f ms)\n", kind ? "LOL" : "ELF", d, 1000.0 * ms[1], 1000.0 * ms[2]);
        EXPECT_LT(d, kind ? 1.0e-5 : 1.0e-3);
    }
    free(v3); free(v2);
    loc_free(&ls);
}

// ELF and LOL against an independent implementation: unittest/elf_ref_data/<name>.ref.txt, made by
// tools/elf_reference.py - the AOs and their gradients by PySCF (double) from the same molden file, rho, grad rho and
// tau of each spin from them, and ELF, LOL by Multiwfn's formulas (function.f90, ELF_LOL: the Kohout-Savin spin form,
// Fc_pol 4.557799872; with and without its D + 1e-5, ELF_addminimal), along 5 lines from the molecule's centre out
// to 12 bohr. md_qm_elf_evaluate (double) and the grid kernel (float, a grid of one row along each line) against
// them, by decade of rho. Where rho is large both agree to rounding; in the tails ELF is a small difference of large
// terms (md_qm_elf_noise), and the basis is float in mdlib: there they part by about that noise.
typedef struct elf_ref_t { double x, y, z, rho, elf, elf_min, lol; } elf_ref_t;

static size_t elf_ref_read(elf_ref_t* out, size_t cap, const char* name) {
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_SOURCE_DIR "/elf_ref_data/%s.ref.txt", name);
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    char line[512];
    size_t n = 0;
    while (n < cap && fgets(line, sizeof(line), f)) {
        if (line[0] == '#') continue;
        elf_ref_t r;
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &r.x, &r.y, &r.z, &r.rho, &r.elf, &r.elf_min, &r.lol) == 7) out[n++] = r;
    }
    fclose(f);
    return n;
}

UTEST(gto_localization, matches_independent) {
    static const char* names[] = {"h2o_ccpvdz", "oh_631gs_uhf", "ammonia_sto3g", "h2o_ccpvqz_cart"};
    enum { LINES = 5, PER_LINE = 49, NB = 13 };
    for (int fi = 0; fi < 4; ++fi) {
        elf_ref_t ref[LINES * PER_LINE];
        ASSERT_EQ(elf_ref_read(ref, LINES * PER_LINE, names[fi]), (size_t)(LINES * PER_LINE));
        char file[128];
        snprintf(file, sizeof(file), "%s.molden", names[fi]);
        loc_src_t s;
        if (!loc_load(&s, file)) { loc_free(&s); ASSERT_TRUE(false); }
        const size_t n = LINES * PER_LINE;
        float* pts = (float*)malloc(sizeof(float) * 3 * n);
        for (size_t i = 0; i < n; ++i) { pts[3*i] = (float)ref[i].x; pts[3*i+1] = (float)ref[i].y; pts[3*i+2] = (float)ref[i].z; }
        // in double on the CPU, where it can be (not with g functions: their derivatives need h, md_qm_kinetic_orbitals)
        double* elf   = (double*)calloc(n, sizeof(double));
        double* elf_m = (double*)calloc(n, sizeof(double));
        double* lol   = (double*)calloc(n, sizeof(double));
        const bool cpu = md_qm_kinetic_available(&s.t.sys);
        if (cpu) {
            ASSERT_EQ(md_qm_elf_evaluate(elf, lol, pts, n, 0, &s.t.sys, s.xyz, 0, 0.0, 0.0), n);
            ASSERT_EQ(md_qm_elf_evaluate(elf_m, NULL, pts, n, 0, &s.t.sys, s.xyz, 0, 0.0, MD_QM_ELF_D_MIN), n);
        }
        // the kernel, a grid of one row along each line (voxel centres on the points)
        float* g_elf   = (float*)calloc(n, sizeof(float));
        float* g_elf_m = (float*)calloc(n, sizeof(float));
        float* g_lol   = (float*)calloc(n, sizeof(float));
        bool gpu = true;
        for (int l = 0; l < LINES && gpu; ++l) {
            const elf_ref_t* r = ref + l * PER_LINE;
            double d[3] = { r[1].x - r[0].x, r[1].y - r[0].y, r[1].z - r[0].z };
            const double h = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
            for (int k = 0; k < 3; ++k) d[k] /= h;
            // two more axes, orthonormal
            double a[3] = { fabs(d[0]) < 0.9 ? 1.0 : 0.0, fabs(d[0]) < 0.9 ? 0.0 : 1.0, 0.0 };
            double e1[3] = { d[1]*a[2] - d[2]*a[1], d[2]*a[0] - d[0]*a[2], d[0]*a[1] - d[1]*a[0] };
            const double l1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
            for (int k = 0; k < 3; ++k) e1[k] /= l1;
            const double e2[3] = { d[1]*e1[2] - d[2]*e1[1], d[2]*e1[0] - d[0]*e1[2], d[0]*e1[1] - d[1]*e1[0] };
            md_grid_t grid = { .spacing = md_float3_set((float)h, (float)h, (float)h), .dim = { PER_LINE, 1, 1 } };
            for (int k = 0; k < 3; ++k) {
                grid.orientation.elem[0][k] = (float)d[k];
                grid.orientation.elem[1][k] = (float)e1[k];
                grid.orientation.elem[2][k] = (float)e2[k];
            }
            const double c[3] = { r[0].x, r[0].y, r[0].z };
            for (int k = 0; k < 3; ++k) grid.origin.elem[k] = (float)(c[k] - 0.5 * h * (d[k] + e1[k] + e2[k]));
            loc_d_min = 0.0;
            gpu = loc_eval_gpu(g_elf + l * PER_LINE, &s, &grid, MD_GTO_GPU_ELF, true, 0.0, NULL) &&
                  loc_eval_gpu(g_lol + l * PER_LINE, &s, &grid, MD_GTO_GPU_LOL, true, 0.0, NULL);
            loc_d_min = MD_QM_ELF_D_MIN;
            gpu = gpu && loc_eval_gpu(g_elf_m + l * PER_LINE, &s, &grid, MD_GTO_GPU_ELF, true, 0.0, NULL);
            loc_d_min = 0.0;
        }
        double e_cpu[NB] = {0}, e_gpu[NB] = {0}, m_cpu[NB] = {0}, m_gpu[NB] = {0}, l_err[NB] = {0}, v[NB] = {0}, vm[NB] = {0};
        size_t cnt[NB] = {0};
        double worst_cpu = 0.0, worst_gpu = 0.0, worst_m_cpu = 0.0, worst_m_gpu = 0.0, worst_lol = 0.0;
        for (size_t i = 0; i < n; ++i) {
            int b = ref[i].rho > 0 ? (int)floor(1.0 - log10(ref[i].rho)) : NB - 1;
            b = b < 0 ? 0 : (b >= NB ? NB - 1 : b);
            const double dc = cpu ? fabs(elf[i] - ref[i].elf) : 0.0, dg = gpu ? fabs((double)g_elf[i] - ref[i].elf) : 0.0;
            const double mc = cpu ? fabs(elf_m[i] - ref[i].elf_min) : 0.0, mg = gpu ? fabs((double)g_elf_m[i] - ref[i].elf_min) : 0.0;
            const double dl = fmax(cpu ? fabs(lol[i] - ref[i].lol) : 0.0, gpu ? fabs((double)g_lol[i] - ref[i].lol) : 0.0);
            e_cpu[b] = fmax(e_cpu[b], dc); e_gpu[b] = fmax(e_gpu[b], dg);
            m_cpu[b] = fmax(m_cpu[b], mc); m_gpu[b] = fmax(m_gpu[b], mg);
            l_err[b] = fmax(l_err[b], dl);
            v[b] += ref[i].elf; vm[b] += ref[i].elf_min; cnt[b]++;
            if (ref[i].rho > 1.0e-4) { worst_cpu = fmax(worst_cpu, dc); worst_gpu = fmax(worst_gpu, dg); }
            worst_m_cpu = fmax(worst_m_cpu, mc); worst_m_gpu = fmax(worst_m_gpu, mg);
            worst_lol = fmax(worst_lol, dl);
        }
        printf("  %-16s against PySCF + Multiwfn's formulas%s (rho decade: points; ELF: mean, max |err| cpu, gpu; with d_min: mean, max |err| cpu, gpu; LOL max |err|)\n",
               names[fi], cpu ? "" : " (no CPU path: g functions)");
        for (int b = 0; b < NB; ++b) if (cnt[b]) {
            printf("    1e%+03d %4zu  %.3f %.1e %.1e   %.3f %.1e %.1e   %.1e\n", 1 - b, cnt[b], v[b] / cnt[b], e_cpu[b], e_gpu[b],
                   vm[b] / cnt[b], m_cpu[b], m_gpu[b], l_err[b]);
        }
        // the function as defined (d_min 0): where rho > 1e-4, the same to rounding (the basis is float in mdlib); with
        // d_min, everywhere
        EXPECT_LT(worst_cpu, 1.0e-4);
        EXPECT_LT(worst_m_cpu, 1.0e-4);
        if (gpu) {
            EXPECT_LT(worst_gpu, 1.0e-3);
            EXPECT_LT(worst_m_gpu, 1.0e-3);
        }
        EXPECT_LT(worst_lol, 1.0e-5);
        free(g_elf); free(g_elf_m); free(g_lol); free(elf); free(elf_m); free(lol); free(pts);
        loc_free(&s);
    }
}

// ELF and LOL out into the tails, where rho vanishes: the kernel with no floor of rho, against the matrices in double
// (tails_dmat_terms: every AO, no cutoff) at every 3rd voxel of a grid out to 9 bohr beyond the outermost atom. With
// d_min (MD_QM_ELF_D_MIN, Multiwfn's convention) ELF is the exact one to within 1e-3 everywhere; without, its error is
// within 30 times md_qm_elf_noise (plus 1e-4): the estimate holds. LOL is within 1e-5 with no floor at all. By decade
// of rho: what ELF is (mean) with d_min 0 and with it, the error without (max, mean, the largest ratio to the noise)
// and with it (max).

// rho, grad rho, tau of one spin's matrix D at a point, all in double: the AOs from the basis (its float exponents
// and coefficients, as the kernel has them) without a cutoff
static void tails_dmat_terms(double* rho, double g[3], double* tau, const loc_src_t* s, const double* D, const float p[3], double* ao, double* aog) {
    size_t n = 0;
    for (uint32_t si = 0; si < s->basis.num_shells; ++si) {
        const md_gto_shell_t sh = s->basis.shells[si];
        const double d[3] = { (double)p[0] - s->xyz[3 * sh.atom_idx], (double)p[1] - s->xyz[3 * sh.atom_idx + 1], (double)p[2] - s->xyz[3 * sh.atom_idx + 2] };
        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        double rad = 0, drad = 0;   // sum c e^(-a r^2), and sum -2 a c e^(-a r^2)
        for (uint32_t k = 0; k < sh.num_primitives; ++k) {
            const double a = s->basis.alpha[sh.primitive_offset + k], c = s->basis.coeff[sh.primitive_offset + k];
            const double e = c * exp(-a * r2);
            rad += e; drad += -2.0 * a * e;
        }
        const int l = (int)sh.l;
        for (int i = l; i >= 0; --i) for (int j = l - i; j >= 0; --j) {
            const int k = l - i - j;
            const int e3[3] = {i, j, k};
            const double f = md_gto_cart_norm_factor(i, j, k);
            double m = 1, dm[3];
            for (int q = 0; q < 3; ++q) m *= pow(d[q], e3[q]);
            for (int q = 0; q < 3; ++q) {
                double o = 1;
                for (int w = 0; w < 3; ++w) if (w != q) o *= pow(d[w], e3[w]);
                dm[q] = (e3[q] > 0 ? e3[q] * pow(d[q], e3[q] - 1) : 0.0) * o;
            }
            ao[n] = f * m * rad;
            for (int q = 0; q < 3; ++q) aog[3 * n + q] = f * (dm[q] * rad + m * drad * d[q]);
            ++n;
        }
    }
    double r = 0, gx = 0, gy = 0, gz = 0, t = 0;
    for (size_t i = 0; i < n; ++i) {
        double u = 0, v[3] = {0, 0, 0};
        for (size_t j = 0; j < n; ++j) {
            const double Dij = D[i * n + j];
            u += Dij * ao[j];
            v[0] += Dij * aog[3 * j]; v[1] += Dij * aog[3 * j + 1]; v[2] += Dij * aog[3 * j + 2];
        }
        r += ao[i] * u;
        gx += 2 * aog[3 * i] * u; gy += 2 * aog[3 * i + 1] * u; gz += 2 * aog[3 * i + 2] * u;
        t += 0.5 * (aog[3 * i] * v[0] + aog[3 * i + 1] * v[1] + aog[3 * i + 2] * v[2]);
    }
    *rho = r; g[0] = gx; g[1] = gy; g[2] = gz; *tau = t;
}

static void loc_ref_dmat(double* elf, double* elf_min, double* lol, double* rho, double* noise, const loc_src_t* s, const float* p, double* ao) {
    double r[2] = {0, 0}, g[2][3] = {{0}}, t[2] = {0, 0};
    for (int k = 0; k < s->num_spins; ++k) tails_dmat_terms(&r[k], g[k], &t[k], s, s->D[k], p, ao, ao + s->num_ao);
    if (s->num_spins == 1) {   // two equal halves
        r[0] *= 0.5; t[0] *= 0.5; for (int q = 0; q < 3; ++q) g[0][q] *= 0.5;
        r[1] = r[0]; t[1] = t[0]; for (int q = 0; q < 3; ++q) g[1][q] = g[0][q];
    }
    const double g2a = g[0][0] * g[0][0] + g[0][1] * g[0][1] + g[0][2] * g[0][2];
    const double g2b = g[1][0] * g[1][0] + g[1][1] * g[1][1] + g[1][2] * g[1][2];
    *elf = md_qm_elf(r[0], g2a, t[0], r[1], g2b, t[1], 0.0, 0.0);
    *elf_min = md_qm_elf(r[0], g2a, t[0], r[1], g2b, t[1], 0.0, MD_QM_ELF_D_MIN);
    *lol = md_qm_lol(r[0], t[0], r[1], t[1], 0.0);
    *rho = r[0] + r[1];
    *noise = md_qm_elf_noise(r[0], t[0], r[1], t[1]);
}

UTEST(gto_localization, tails) {
    static const char* files[] = {"h2o.h5", "mol.h5", "oh_631gs_uhf.molden"};
    for (int fi = 0; fi < 3; ++fi) {
        loc_src_t s;
        if (!loc_load(&s, files[fi])) { loc_free(&s); ASSERT_TRUE(false); }
        double c[3] = {0, 0, 0}, ext = 0;
        for (size_t a = 0; a < s.num_atoms; ++a) for (int k = 0; k < 3; ++k) c[k] += s.xyz[3 * a + k] / (double)s.num_atoms;
        for (size_t a = 0; a < s.num_atoms; ++a) for (int k = 0; k < 3; ++k) ext = fmax(ext, fabs(s.xyz[3 * a + k] - c[k]));
        const int dim = 32;
        const float h = (float)(2.0 * (ext + 9.0) / dim);
        const md_grid_t grid = {
            .orientation = md_float3x3_ident(),
            .origin  = md_float3_set((float)c[0] - 0.5f * h * dim, (float)c[1] - 0.5f * h * dim, (float)c[2] - 0.5f * h * dim),
            .spacing = md_float3_set(h, h, h),
            .dim     = { dim, dim, dim },
        };
        const size_t n = (size_t)dim * dim * dim;
        float* elf_raw = (float*)malloc(sizeof(float) * n);
        float* elf_min = (float*)malloc(sizeof(float) * n);
        float* lol = (float*)malloc(sizeof(float) * n);
        loc_d_min = 0.0;
        bool ok = loc_eval_gpu(elf_raw, &s, &grid, MD_GTO_GPU_ELF, true, 0.0, NULL) && loc_eval_gpu(lol, &s, &grid, MD_GTO_GPU_LOL, true, 0.0, NULL);
        loc_d_min = MD_QM_ELF_D_MIN;
        ok = ok && loc_eval_gpu(elf_min, &s, &grid, MD_GTO_GPU_ELF, true, 0.0, NULL);
        loc_d_min = 0.0;
        if (!ok) { free(elf_raw); free(elf_min); free(lol); loc_free(&s); UTEST_SKIP("no GPU device"); }

        double* ao = (double*)malloc(sizeof(double) * 4 * s.num_ao);
        enum { NB = 14 };   // decades 1e1 .. 1e-12 and below
        double val[NB] = {0}, val_min[NB] = {0}, raw[NB] = {0}, rawsum[NB] = {0}, ratio[NB] = {0}, err_min[NB] = {0};
        size_t cnt[NB] = {0};
        double worst_min = 0, worst_lol = 0, worst_ratio = 0;
        for (size_t i = 0; i < n; i += 3) {
            const size_t x = i % dim, y = (i / dim) % dim, z = i / ((size_t)dim * dim);
            const float p[3] = { grid.origin.x + ((float)x + 0.5f) * h, grid.origin.y + ((float)y + 0.5f) * h, grid.origin.z + ((float)z + 0.5f) * h };
            double e, em, l, rho, noise;
            loc_ref_dmat(&e, &em, &l, &rho, &noise, &s, p, ao);
            const double r = fabs((double)elf_raw[i] - e);
            const double dm = fabs((double)elf_min[i] - em);
            worst_min = fmax(worst_min, dm);
            worst_lol = fmax(worst_lol, fabs((double)lol[i] - l));
            int b = rho > 0 ? (int)floor(1.0 - log10(rho)) : NB - 1;
            b = b < 0 ? 0 : (b >= NB ? NB - 1 : b);
            // where the noise is still below 0.2 (beyond it the values are noise)
            const double q = noise < 0.2 ? fmax(r - 1.0e-4, 0.0) / noise : 0.0;
            worst_ratio = fmax(worst_ratio, q);
            val[b] += e; val_min[b] += em; raw[b] = fmax(raw[b], r); rawsum[b] += r; ratio[b] = fmax(ratio[b], q);
            err_min[b] = fmax(err_min[b], dm); cnt[b]++;
        }
        printf("  %-20s with d_min within %.1e; without, error / noise at most %.1f; LOL within %.1e (" MD_GTO_GPU_LOCALIZATION_REVISION ")\n",
               files[fi], worst_min, worst_ratio, worst_lol);
        printf("    rho     voxels  ELF    with d_min  max err   mean err  err/noise  max err with d_min\n");
        for (int b = 0; b < NB; ++b) if (cnt[b]) {
            printf("    1e%+03d %6zu  %.3f  %.3f       %.2e  %.2e  %6.1f     %.2e\n", 1 - b, cnt[b], val[b] / cnt[b], val_min[b] / cnt[b],
                   raw[b], rawsum[b] / cnt[b], ratio[b], err_min[b]);
        }
        EXPECT_LT(worst_ratio, 30.0);
        EXPECT_LT(worst_min, 1.0e-3);
        EXPECT_LT(worst_lol, 1.0e-5);
        free(ao); free(elf_raw); free(elf_min); free(lol);
        loc_free(&s);
    }
}
#endif // MD_ENABLE_GPU
