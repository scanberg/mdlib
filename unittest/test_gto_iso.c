// PROTOTYPE - gridless isosurface extraction for Gaussian orbitals.
//
// The question this file answers: can an orbital isosurface be built directly from the analytic
// function, never manifesting a volume, with a guarantee that no part of the surface is missed?
//
// The pipeline:
//   1. Bound psi over an axis aligned box (four schemes, compared below).
//   2. A lazily built octree: a cell survives only if its bound admits psi = +c or psi = -c.
//      Shells whose contribution over a cell is provably tiny are dropped for the whole subtree,
//      their bound carried down as slack, so deep cells evaluate few shells.
//   3. Surface Nets on the surviving leaves (all at one size h, so no cracks).
//   4. Newton projection onto the true surface with the analytic gradient.
//   5. Tangential relaxation + reprojection.
//   6. Optional adaptive refinement: split edges whose midpoint is far from the surface.
//
// The bound is the whole game. Too loose wastes refinement; too tight silently loses surface and
// looks exactly like a correct result. So the tests here are built around catching that:
//   - every bound is fuzzed for containment against point samples,
//   - the sparse Surface Nets mesh must be IDENTICAL to dense Surface Nets over the same lattice
//     (a kept cell set that misses one crossing cell shows up as a count mismatch).
//
// Bound schemes, for psi = sum_s P_s(d) R_s(|d|^2), d = r - centre_s:
//   PRIM  per primitive, exact: x^i e^{-a x^2} factorizes per axis. Loose under cancellation.
//   GROUP per (atom, l) group: sum_n range(x^i y^j z^k over box) * range(R_n over [rho_min, rho_max]).
//         R_n is 1D; its range is exact from precomputed critical intervals, which removes the
//         cancellation between the shells of one atom (radial nodes).
//   T1    psi(c) + sum_d [d psi/dx_d](B) * [-r_d, r_d]          (mean value form)
//   T2    psi(c) + grad psi(c).delta + 1/2 delta^T H(B) delta   (Taylor model, order 2)
//         first order part is exact, only the Hessian enclosure is loose and it scales with r^2.
//   BEST  PRIM n GROUP n T1 n T2
//
// The GPU version (bottom of this file, shaders/gto_iso.slang) runs the shape the CPU experiments
// pointed at on md_gpu: GROUP bound octree down to 4^3 cell bricks, dense corners per brick,
// Surface Nets, projection, relaxation, triangles. gto_iso_gpu.matches_cpu holds it to the double
// precision CPU path on the same lattice.
//
// Environment:
//   GTO_ISO_VALIDATE=1  create the md_gpu device with Vulkan validation
//   GTO_ISO_REPORT=1    also run the scheme comparison and the full mesh pipeline report (~1 min)
//   GTO_ISO_OUT=<dir>   write OBJ meshes and CSV tables there
//   GTO_ISO_FULL=1      report + bigger fuzz counts, finer meshes (h 0.08), all orbitals (slow)
//   GTO_ISO_H=<bohr>    override the pipeline cell size
//   GTO_ISO_ONLY=a,b    run the pipeline for these orbitals only (e.g. amide_homo)
//
// Densities (gto_iso_density.*): SCF, attachment, detachment, difference and transition densities of
// the test files, each as products of orbitals (MOs, NTOs) and as an AO matrix the field splits
// itself; checked against the orbital evaluator, fuzzed for containment (both density bounds, with
// the octree's drop slack), held to dense Surface Nets, and timed per split (GTO_ISO_REPORT=1).

#include "utest.h"
#include "vlx_test_util.h"

#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_os.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <md_gto_iso.c>   // the implementation, internals included

static const char* iso_bound_name[ISO_BOUND_COUNT] = {"prim", "group", "taylor1", "taylor2", "best"};

// ---------------------------------------------------------------------------------------------
// Test utilities
// ---------------------------------------------------------------------------------------------

typedef struct iso_rng_t { uint64_t s; } iso_rng_t;

static inline uint64_t iso_rng_u64(iso_rng_t* r) {
    // splitmix64
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
static inline double iso_rng_f64(iso_rng_t* r) { return (double)(iso_rng_u64(r) >> 11) * (1.0 / 9007199254740992.0); }
static inline double iso_rng_range(iso_rng_t* r, double a, double b) { return a + (b - a) * iso_rng_f64(r); }
static inline double iso_rng_normal(iso_rng_t* r) {
    double u1 = fmax(iso_rng_f64(r), 1e-300), u2 = iso_rng_f64(r);
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

static int iso_cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static double iso_percentile(double* v, int n, double q) {
    if (n == 0) return NAN;
    qsort(v, n, sizeof(double), iso_cmp_double);
    int i = (int)floor(q * (n - 1) + 0.5);
    return v[i];
}

// ---------------------------------------------------------------------------------------------
// Mesh metrics: topology, triangle quality, deviation from the true surface
// ---------------------------------------------------------------------------------------------

typedef struct iso_metrics_t {
    int    nv, ne, nt;
    int    euler;
    int    components;
    int    boundary_edges;
    int    nonmanifold_edges;
    double min_angle, p1_angle, mean_min_angle, frac_below_20;
    double dev_max, dev_rms;        // at triangle centroids and edge midpoints, units of h
    double vert_res_max;            // units of h
    double orient_ok;               // fraction of triangles whose winding agrees with the field normal
    double dev_where[3];            // where dev_max occurred
    double dev_grad;                // |grad f| there
} iso_metrics_t;

static int iso_uf_find(int* p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

static void iso_mesh_metrics(iso_metrics_t* out, const iso_mesh_t* m, const iso_field_t* f) {
    memset(out, 0, sizeof(*out));
    const int nv = iso_mesh_nv(m), nt = iso_mesh_nt(m);
    out->nv = nv; out->nt = nt;

    iso_map_t edges;
    iso_map_init(&edges, (uint32_t)(nt * 2));
    int* parent = (int*)malloc(sizeof(int) * (nv + 1));
    for (int i = 0; i < nv; ++i) parent[i] = i;
    double* mins = (double*)malloc(sizeof(double) * (nt + 1));
    double dev2 = 0.0;
    long ndev = 0;
    double sum_min = 0.0;
    int below = 0;

    for (int t = 0; t < nt; ++t) {
        const uint32_t* v = m->tri + 3 * t;
        for (int e = 0; e < 3; ++e) {
            uint32_t a = v[e], b = v[(e + 1) % 3];
            if (a > b) { uint32_t x = a; a = b; b = x; }
            const uint64_t key = ((uint64_t)a << 32) | b;
            uint32_t c = 0;
            iso_map_get(&edges, key, &c);
            iso_map_put(&edges, key, c + 1);
            const int ra = iso_uf_find(parent, (int)a), rb = iso_uf_find(parent, (int)b);
            if (ra != rb) parent[ra] = rb;
        }
        // angles
        double mn = 180.0;
        for (int k = 0; k < 3; ++k) {
            const double* p0 = m->pos + 3 * v[k];
            const double* p1 = m->pos + 3 * v[(k + 1) % 3];
            const double* p2 = m->pos + 3 * v[(k + 2) % 3];
            const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const double w[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const double cu = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
            const double cw = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            double ca = (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / fmax(cu * cw, 1e-300);
            ca = fmax(-1.0, fmin(1.0, ca));
            mn = fmin(mn, acos(ca) * 180.0 / 3.14159265358979);
        }
        {
            const double* p0 = m->pos + 3 * v[0]; const double* p1 = m->pos + 3 * v[1]; const double* p2 = m->pos + 3 * v[2];
            const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const double w[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const double cr[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
            double dotn = 0.0;
            if (m->nrm) for (int k = 0; k < 3; ++k) for (int a = 0; a < 3; ++a) dotn += cr[a] * m->nrm[3 * v[k] + a];
            out->orient_ok += dotn > 0.0;
        }
        mins[t] = mn;
        sum_min += mn;
        below += mn < 20.0;
        // deviation at centroid and edge midpoints (first order distance |f|/|grad f|)
        for (int k = 0; k < 4; ++k) {
            double p[3];
            for (int a = 0; a < 3; ++a) {
                p[a] = (k == 3) ? (m->pos[3 * v[0] + a] + m->pos[3 * v[1] + a] + m->pos[3 * v[2] + a]) / 3.0
                                : 0.5 * (m->pos[3 * v[k] + a] + m->pos[3 * v[(k + 1) % 3] + a]);
            }
            double val, g[3];
            iso_eval(f, NULL, 0, p, 1, &val, g, NULL);
            const double gl = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            double d = fabs(m->sheet * val - m->iso) / fmax(gl, 1e-300) / m->h;
            if (d > 0.25) {
                // the first order distance overestimates where |grad f| is small (near nodes and
                // extrema): walk to the surface in small steps and take the distance travelled
                double q[3] = {p[0], p[1], p[2]}, nn[3];
                const double res = iso_project_point(f, q, m->sheet, m->iso, 0.02 * m->h, 400, nn);
                if (res < 1.0e-3 * m->h) d = fmin(d, sqrt(iso_dist2(p, q)) / m->h);
            }
            if (d > out->dev_max) {
                out->dev_max = d;
                memcpy(out->dev_where, p, sizeof(p));
                out->dev_grad = gl;
            }
            dev2 += d * d;
            ndev++;
        }
    }
    for (uint32_t i = 0; i < edges.cap; ++i) {
        if (!edges.keys[i]) continue;
        out->ne++;
        if (edges.vals[i] == 1) out->boundary_edges++;
        if (edges.vals[i] > 2)  out->nonmanifold_edges++;
    }
    for (int i = 0; i < nv; ++i) if (iso_uf_find(parent, i) == i) out->components++;
    out->euler = out->nv - out->ne + out->nt;
    out->orient_ok = nt ? out->orient_ok / nt : 0.0;
    out->dev_rms = ndev ? sqrt(dev2 / ndev) : 0.0;
    out->min_angle = nt ? iso_percentile(mins, nt, 0.0) : 0.0;
    out->p1_angle = nt ? iso_percentile(mins, nt, 0.01) : 0.0;
    out->mean_min_angle = nt ? sum_min / nt : 0.0;
    out->frac_below_20 = nt ? (double)below / nt : 0.0;
    if (m->res) for (int i = 0; i < nv; ++i) out->vert_res_max = fmax(out->vert_res_max, m->res[i] / m->h);
    free(mins); free(parent);
    iso_map_free(&edges);
}

static void iso_mesh_write_obj(FILE* fp, const iso_mesh_t* m, int vertex_base, const char* group) {
    const int nv = iso_mesh_nv(m), nt = iso_mesh_nt(m);
    fprintf(fp, "g %s\n", group);
    for (int i = 0; i < nv; ++i) fprintf(fp, "v %.6f %.6f %.6f\n", m->pos[3 * i], m->pos[3 * i + 1], m->pos[3 * i + 2]);
    if (m->nrm) for (int i = 0; i < nv; ++i) fprintf(fp, "vn %.5f %.5f %.5f\n", m->nrm[3 * i], m->nrm[3 * i + 1], m->nrm[3 * i + 2]);
    for (int t = 0; t < nt; ++t) {
        const int a = vertex_base + (int)m->tri[3 * t] + 1, b = vertex_base + (int)m->tri[3 * t + 1] + 1, c = vertex_base + (int)m->tri[3 * t + 2] + 1;
        if (m->nrm) fprintf(fp, "f %d//%d %d//%d %d//%d\n", a, a, b, b, c, c);
        else        fprintf(fp, "f %d %d %d\n", a, b, c);
    }
}


// ---------------------------------------------------------------------------------------------
// Test orbitals
// ---------------------------------------------------------------------------------------------

enum { ISO_MO_HOMO = -1, ISO_MO_LUMO = -2, ISO_MO_RANDOM = -3 };

// Loads one orbital of a VeloxChem file as a field. 'which' is an MO index or one of ISO_MO_*.
// ISO_MO_RANDOM: Gaussian random AO coefficients scaled to the HOMO's norm - lots of cancellation.
static bool iso_load(iso_field_t* f, const char* file, int which, const char* name) {
    vlx_test_t t = {0};
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    if (!vlx_test_load(&t, str_from_cstr(path), MEGABYTES(64))) return false;

    bool ok = false;
    const size_t num_atoms = t.sys.atom.count;
    float* xyz = (float*)malloc(sizeof(float) * 3 * num_atoms);
    md_gto_basis_t basis = {0};
    double* c = NULL;
    if (vlx_test_atom_xyz_bohr(xyz, 3 * num_atoms, &t) != num_atoms) goto done;
    if (!qm_test_basis(&basis, &t)) goto done;

    const size_t num_ao = md_gto_basis_num_ao(&basis);
    c = (double*)calloc(num_ao, sizeof(double));
    const size_t lumo = vlx_test_lumo_idx(&t, STR_LIT("orbital/alpha/occupation"));
    size_t idx = (which >= 0) ? (size_t)which : (which == ISO_MO_LUMO ? lumo : lumo - 1);
    if (qm_test_row(c, num_ao, &t, STR_LIT("orbital/alpha/coefficient"), idx) != num_ao) goto done;

    if (which == ISO_MO_RANDOM) {
        double n0 = 0.0, n1 = 0.0;
        iso_rng_t rng = {12345};
        for (size_t i = 0; i < num_ao; ++i) n0 += c[i] * c[i];
        for (size_t i = 0; i < num_ao; ++i) { c[i] = iso_rng_normal(&rng); n1 += c[i] * c[i]; }
        for (size_t i = 0; i < num_ao; ++i) c[i] *= sqrt(n0 / n1);
    }
    ok = iso_field_init(f, &basis, xyz, num_atoms, c, name);
done:
    free(xyz);
    free(c);
    qm_test_free(&t);
    return ok;
}

typedef struct iso_case_t { const char* file; int which; const char* name; double iso; } iso_case_t;

static const iso_case_t iso_cases[] = {
    {"h2o.h5",   ISO_MO_LUMO,   "h2o_lumo",     0.05},
    {"h2o.h5",   ISO_MO_HOMO,   "h2o_homo",     0.05},
    {"amide.h5", ISO_MO_HOMO,   "amide_homo",   0.05},
    {"amide.h5", ISO_MO_LUMO,   "amide_lumo",   0.05},
    {"amide.h5", 0,             "amide_core",   0.05},
    {"amide.h5", ISO_MO_RANDOM, "amide_random", 0.05},
    {"mol.h5",   ISO_MO_HOMO,   "mol_homo",     0.05},
    {"mol.h5",   ISO_MO_LUMO,   "mol_lumo",     0.05},
};
#define ISO_NUM_CASES (int)(sizeof(iso_cases) / sizeof(iso_cases[0]))

static bool iso_full(void) { const char* e = getenv("GTO_ISO_FULL"); return e && e[0] == '1'; }
static const char* iso_out_dir(void) { return getenv("GTO_ISO_OUT"); }

static FILE* iso_open_out(const char* name, const char* mode) {
    const char* dir = iso_out_dir();
    if (!dir) return NULL;
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    return fopen(path, mode);
}

// ---------------------------------------------------------------------------------------------
// The evaluator must agree with mdlib's, and its derivatives with finite differences
// ---------------------------------------------------------------------------------------------

UTEST(gto_iso, point_eval_matches_mdlib) {
    static const int pick[3] = {0, 2, 6};
    for (int ci = 0; ci < 3; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        // reference: mdlib's own float evaluator over the expanded primitives
        vlx_test_t t = {0};
        char path[512];
        snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", tc->file);
        ASSERT_TRUE(vlx_test_load(&t, str_from_cstr(path), MEGABYTES(64)));
        const size_t num_atoms = t.sys.atom.count;
        float* xyz = (float*)malloc(sizeof(float) * 3 * num_atoms);
        ASSERT_EQ(vlx_test_atom_xyz_bohr(xyz, 3 * num_atoms, &t), num_atoms);
        md_gto_basis_t basis = {0};
        ASSERT_TRUE(qm_test_basis(&basis, &t));
        const size_t num_ao = md_gto_basis_num_ao(&basis);
        double* c = (double*)calloc(num_ao, sizeof(double));
        const size_t lumo = vlx_test_lumo_idx(&t, STR_LIT("orbital/alpha/occupation"));
        ASSERT_EQ(qm_test_row(c, num_ao, &t, STR_LIT("orbital/alpha/coefficient"), tc->which == ISO_MO_LUMO ? lumo : lumo - 1), num_ao);

        iso_field_t f;
        ASSERT_TRUE(iso_field_init(&f, &basis, xyz, num_atoms, c, tc->name));

        const size_t num_gto = md_gto_pgto_count(&basis);
        md_gto_t* gtos = (md_gto_t*)malloc(sizeof(md_gto_t) * num_gto);
        const size_t n_gto = md_gto_expand_with_ao_coeffs(gtos, &basis, xyz, 0, c, 0.0);

        enum { N = 2000 };
        float* pts = (float*)malloc(sizeof(float) * 3 * N);
        float* ref = (float*)malloc(sizeof(float) * N);
        iso_rng_t rng = {7};
        for (int i = 0; i < N; ++i) {
            for (int a = 0; a < 3; ++a) pts[3 * i + a] = (float)iso_rng_range(&rng, f.atom_min[a] - 3.0, f.atom_max[a] + 3.0);
        }
        md_gto_xyz_evaluate(ref, pts, N, 0, gtos, n_gto, MD_GTO_EVAL_MODE_PSI);

        double max_ref = 0.0, max_err = 0.0, max_gerr = 0.0, max_herr = 0.0, max_g = 0.0, max_h = 0.0;
        for (int i = 0; i < N; ++i) {
            const double p[3] = {pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]};
            double v, g[3], h[6];
            iso_eval(&f, NULL, 0, p, 2, &v, g, h);
            max_ref = fmax(max_ref, fabs(ref[i]));
            max_err = fmax(max_err, fabs(v - ref[i]));

            // central differences: gradient from values, Hessian from gradients
            const double e = 1.0e-4;
            for (int a = 0; a < 3; ++a) {
                double pp[3] = {p[0], p[1], p[2]}, pm[3] = {p[0], p[1], p[2]};
                pp[a] += e; pm[a] -= e;
                double vp, vm, gp[3], gm[3];
                iso_eval(&f, NULL, 0, pp, 1, &vp, gp, NULL);
                iso_eval(&f, NULL, 0, pm, 1, &vm, gm, NULL);
                max_g = fmax(max_g, fabs(g[a]));
                max_gerr = fmax(max_gerr, fabs((vp - vm) / (2 * e) - g[a]));
                static const int hi[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
                for (int b = 0; b < 3; ++b) {
                    max_h = fmax(max_h, fabs(h[hi[a][b]]));
                    max_herr = fmax(max_herr, fabs((gp[b] - gm[b]) / (2 * e) - h[hi[a][b]]));
                }
            }
        }
        printf("  %-14s |psi|max %.3e  value err %.2e  grad err %.2e (|g|max %.2e)  hess err %.2e (|H|max %.2e)\n",
            tc->name, max_ref, max_err, max_gerr, max_g, max_herr, max_h);
        EXPECT_LT(max_err, 1.0e-5 * max_ref + 1.0e-6);
        EXPECT_LT(max_gerr, 1.0e-5 * max_g);
        EXPECT_LT(max_herr, 1.0e-4 * max_h);

        free(pts); free(ref); free(gtos); free(c); free(xyz);
        iso_field_free(&f);
        qm_test_free(&t);
    }
}

// The radial ranges R, R', R'' over [a, b] must contain dense samples of the function.
UTEST(gto_iso, radial_ranges_contain_samples) {
    static const int pick[3] = {2, 4, 5};
    for (int ci = 0; ci < 3; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_field_t f;
        ASSERT_TRUE(iso_load(&f, tc->file, tc->which, tc->name));
        iso_rng_t rng = {99};
        int violations = 0, num_crit = 0, num_res = 0;
        double worst = 0.0;
        for (int si = 0; si < f.num_shells; ++si) {
            const iso_shell_t* s = &f.shells[si];
            for (int n = 0; n < s->ncart; ++n) {
                if (!s->live[n]) continue;
                num_res++;
                for (int m = 0; m < ISO_NUM_CRIT; ++m) num_crit += s->crit_num[n][m];
                for (int trial = 0; trial < 60; ++trial) {
                    // log-uniform so the tight core region is sampled as well as the tail
                    double a = exp(iso_rng_range(&rng, log(1e-6), log(60.0)));
                    double b = a * exp(iso_rng_range(&rng, 0.0, 3.0));
                    if (iso_rng_f64(&rng) < 0.1) a = 0.0;
                    for (int m = 0; m < ISO_NUM_CRIT; ++m) {
                        const iv_t r = iso_g_range(&f, s, n, m, a, b);
                        const double scale = fmax(iv_mag(r), 1e-300);
                        for (int k = 0; k <= 400; ++k) {
                            const double rho = a + (b - a) * k / 400.0;
                            const double g = iso_g(s, n, m, rho);
                            const double tol = 1e-12 * scale;
                            if (g < r.lo - tol || g > r.hi + tol) {
                                violations++;
                                worst = fmax(worst, fmax(r.lo - g, g - r.hi) / scale);
                            }
                        }
                    }
                }
            }
        }
        printf("  %-14s input shells %d -> (atom,l) groups %d, radial components %d, critical intervals %d, violations %d (worst rel %.2e)\n",
            tc->name, f.num_input_shells, f.num_shells, num_res, num_crit, violations, worst);
        EXPECT_EQ(violations, 0);
        iso_field_free(&f);
    }
}

// ---------------------------------------------------------------------------------------------
// Containment fuzzing and tightness curves
// ---------------------------------------------------------------------------------------------

#define ISO_NUM_SIZES 8
static const double iso_sizes[ISO_NUM_SIZES] = {4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125};

UTEST(gto_iso, bounds_contain_samples) {
    const int boxes_per_size = iso_full() ? 2000 : 300;
    const int pts_per_box    = iso_full() ? 256 : 48;
    const int num_cases      = iso_full() ? ISO_NUM_CASES : 6;

    FILE* csv = iso_open_out("gto_iso_tightness.csv", "w");
    if (csv) fprintf(csv, "orbital,size,scheme,category,median,p90,count\n");

    int total_violations[ISO_BOUND_COUNT] = {0};
    long total_checks = 0;

    // pooled ratio statistics over all orbitals, near-surface boxes
    double* pool[ISO_NUM_SIZES][ISO_BOUND_COUNT];
    int     pool_n[ISO_NUM_SIZES][ISO_BOUND_COUNT];
    for (int z = 0; z < ISO_NUM_SIZES; ++z) for (int k = 0; k < ISO_BOUND_COUNT; ++k) {
        pool[z][k] = (double*)malloc(sizeof(double) * boxes_per_size * ISO_NUM_CASES);
        pool_n[z][k] = 0;
    }

    for (int ci = 0; ci < num_cases; ++ci) {
        const iso_case_t* tc = &iso_cases[ci];
        iso_field_t f;
        if (!iso_load(&f, tc->file, tc->which, tc->name)) { ASSERT_TRUE(false); }
        int* all = (int*)malloc(sizeof(int) * f.num_shells);
        for (int i = 0; i < f.num_shells; ++i) all[i] = i;

        const int num_atoms_guess = f.num_shells;  // centres to sample near: shell centres
        iso_rng_t rng = {1000 + ci};
        double* ratios[ISO_BOUND_COUNT];
        for (int k = 0; k < ISO_BOUND_COUNT; ++k) ratios[k] = (double*)malloc(sizeof(double) * boxes_per_size);

        for (int z = 0; z < ISO_NUM_SIZES; ++z) {
            int nr = 0;
            for (int bi = 0; bi < boxes_per_size; ++bi) {
                double cen[3];
                if (bi & 1) {
                    const iso_shell_t* s = &f.shells[iso_rng_u64(&rng) % num_atoms_guess];
                    for (int a = 0; a < 3; ++a) cen[a] = s->c[a] + 0.8 * iso_rng_normal(&rng);
                } else {
                    for (int a = 0; a < 3; ++a) cen[a] = iso_rng_range(&rng, f.atom_min[a] - 4.0, f.atom_max[a] + 4.0);
                }
                iso_box_t b;
                for (int a = 0; a < 3; ++a) {
                    const double half = 0.5 * iso_sizes[z] * iso_rng_range(&rng, 0.5, 1.0);
                    b.lo[a] = cen[a] - half;
                    b.hi[a] = cen[a] + half;
                }
                iso_encl_t enc;
                iso_enclose(&enc, &f, all, f.num_shells, 0.0, &b, 0x1F);

                // observed range: corners, centre, random interior points
                double omin = DBL_MAX, omax = -DBL_MAX;
                for (int k = 0; k < 9 + pts_per_box; ++k) {
                    double p[3];
                    for (int a = 0; a < 3; ++a) {
                        if (k < 8)       p[a] = (k >> a) & 1 ? b.hi[a] : b.lo[a];
                        else if (k == 8) p[a] = 0.5 * (b.lo[a] + b.hi[a]);
                        else             p[a] = iso_rng_range(&rng, b.lo[a], b.hi[a]);
                    }
                    double v;
                    iso_eval(&f, all, f.num_shells, p, 0, &v, NULL, NULL);
                    omin = fmin(omin, v);
                    omax = fmax(omax, v);
                    for (int s = 0; s < ISO_BOUND_COUNT; ++s) {
                        const double tol = 1e-13 * (enc.mag + fabs(v));
                        if (v < enc.e[s].lo - tol || v > enc.e[s].hi + tol) {
                            if (total_violations[s] < 5) {
                                printf("  VIOLATION %s %s size %.4f v=%.9e enc=[%.9e, %.9e]\n", tc->name, iso_bound_name[s], iso_sizes[z], v, enc.e[s].lo, enc.e[s].hi);
                            }
                            total_violations[s]++;
                        }
                    }
                    total_checks++;
                }
                // tightness only where it matters: boxes the surface passes through
                const bool near = (omin <= tc->iso && omax >= tc->iso) || (omin <= -tc->iso && omax >= -tc->iso);
                if (near && omax - omin > 0.0) {
                    for (int s = 0; s < ISO_BOUND_COUNT; ++s) {
                        const double r = iv_width(enc.e[s]) / (omax - omin);
                        ratios[s][nr] = r;
                        pool[z][s][pool_n[z][s]++] = r;
                    }
                    nr++;
                }
            }
            if (csv) {
                for (int s = 0; s < ISO_BOUND_COUNT; ++s) {
                    double tmp[4096];
                    const int n = nr < 4096 ? nr : 4096;
                    memcpy(tmp, ratios[s], sizeof(double) * n);
                    const double med = iso_percentile(tmp, n, 0.5), p90 = iso_percentile(tmp, n, 0.9);
                    fprintf(csv, "%s,%g,%s,near_surface,%g,%g,%d\n", tc->name, iso_sizes[z], iso_bound_name[s], med, p90, n);
                }
            }
        }
        for (int k = 0; k < ISO_BOUND_COUNT; ++k) free(ratios[k]);
        free(all);
        iso_field_free(&f);
    }
    if (csv) fclose(csv);

    printf("\n  Overestimation ratio width(bound)/width(sampled), boxes crossing the surface, median [p90]\n");
    printf("  %8s", "size");
    for (int s = 0; s < ISO_BOUND_COUNT; ++s) printf(" %20s", iso_bound_name[s]);
    printf("\n");
    for (int z = 0; z < ISO_NUM_SIZES; ++z) {
        printf("  %8.4f", iso_sizes[z]);
        for (int s = 0; s < ISO_BOUND_COUNT; ++s) {
            const double med = iso_percentile(pool[z][s], pool_n[z][s], 0.5);
            const double p90 = iso_percentile(pool[z][s], pool_n[z][s], 0.9);
            char buf[64];
            snprintf(buf, sizeof(buf), "%.3g [%.3g]", med, p90);
            printf(" %20s", buf);
        }
        printf("   n=%d\n", pool_n[z][0]);
    }
    printf("\n  containment checks %ld:", total_checks);
    for (int s = 0; s < ISO_BOUND_COUNT; ++s) printf(" %s=%d", iso_bound_name[s], total_violations[s]);
    printf(" violations\n");

    for (int z = 0; z < ISO_NUM_SIZES; ++z) for (int k = 0; k < ISO_BOUND_COUNT; ++k) free(pool[z][k]);
    for (int s = 0; s < ISO_BOUND_COUNT; ++s) EXPECT_EQ(total_violations[s], 0);
}

// ---------------------------------------------------------------------------------------------
// Cells the bound kept but whose corners show no crossing: is anything hiding inside?
// ---------------------------------------------------------------------------------------------


// 0: proven empty, 1: sign change found (a feature smaller than the cell), 2: unresolved at max depth
static int iso_hidden_rec(const iso_field_t* f, const int* all, const iso_box_t* b, int sheet, double iso, bool ref_inside, int depth, int max_depth) {
    iso_encl_t enc;
    iso_enclose(&enc, f, all, f->num_shells, 0.0, b, 1u << ISO_BOUND_SHELL);
    if (!iv_contains(enc.e[ISO_BOUND_SHELL], sheet * iso)) return 0;
    iso_enclose(&enc, f, all, f->num_shells, 0.0, b, 1u << ISO_BOUND_T2);
    if (!iv_contains(enc.e[ISO_BOUND_T2], sheet * iso)) return 0;
    for (int k = 0; k < 27; ++k) {
        const int idx[3] = {k % 3, (k / 3) % 3, k / 9};
        double p[3];
        for (int a = 0; a < 3; ++a) p[a] = b->lo[a] + 0.5 * idx[a] * (b->hi[a] - b->lo[a]);
        if ((sheet * iso_value(f, p) - iso > 0) != ref_inside) return 1;
    }
    if (depth == max_depth) return 2;
    int worst = 0;
    for (int c = 0; c < 8; ++c) {
        iso_box_t cb;
        for (int a = 0; a < 3; ++a) {
            const double mid = 0.5 * (b->lo[a] + b->hi[a]);
            cb.lo[a] = ((c >> a) & 1) ? mid : b->lo[a];
            cb.hi[a] = ((c >> a) & 1) ? b->hi[a] : mid;
        }
        const int r = iso_hidden_rec(f, all, &cb, sheet, iso, ref_inside, depth + 1, max_depth);
        if (r == 1) return 1;
        if (r == 2) worst = 2;
    }
    return worst;
}

typedef struct iso_hidden_cell_t { double p[3]; int sheet; } iso_hidden_cell_t;

// Every cell the octree kept but where Surface Nets sees no crossing must be proven empty, or it is
// reported. Cheap first: a Lipschitz test from a gradient bound per brick (|f| at the nearest corner
// exceeds what the gradient allows over half a cell diagonal). Only the remainder gets the
// recursive bound + sampling check.
typedef struct iso_guarantee_t {
    long candidates, by_trilinear, by_lipschitz, by_recursion, hidden, unresolved;
    double seconds;
    iso_hidden_cell_t* hidden_cells;   // md_array
} iso_guarantee_t;

static void iso_guarantee(iso_guarantee_t* r, const iso_octree_t* o, int k, iso_lattice_t* lat, const iso_leaf_t* cells, int nc, int max_depth) {
    memset(r, 0, sizeof(*r));
    md_tick_t t0 = md_tick_now();
    const iso_field_t* f = o->f;
    int* all = (int*)malloc(sizeof(int) * f->num_shells);
    for (int i = 0; i < f->num_shells; ++i) all[i] = i;

    // gradient bound per brick (all groups: a dropped group's gradient is not covered by the slack)
    const int nl = (int)md_array_size(o->leaves);
    iso_map_t brick_idx;
    iso_map_init(&brick_idx, (uint32_t)nl);
    double* L = (double*)malloc(sizeof(double) * (nl + 1));
    double* E = (double*)malloc(sizeof(double) * (nl + 1));   // trilinear interpolation error bound per cell
    const double bh = o->h;   // brick edge
    for (int li = 0; li < nl; ++li) {
        const iso_leaf_t* l = &o->leaves[li];
        iso_box_t b;
        const int idx[3] = {l->x, l->y, l->z};
        for (int a = 0; a < 3; ++a) { b.lo[a] = o->origin[a] + idx[a] * bh; b.hi[a] = b.lo[a] + bh; }
        iso_encl_t enc;
        iso_enclose(&enc, f, all, f->num_shells, 0.0, &b, (1u << ISO_BOUND_T1) | (1u << ISO_BOUND_T2));
        L[li] = enc.grad_mag;
        E[li] = lat->h * lat->h / 8.0 * (enc.hess_diag[0] + enc.hess_diag[1] + enc.hess_diag[2]);
        iso_map_put(&brick_idx, iso_key3(l->x, l->y, l->z), (uint32_t)li);
    }

    const double half_diag = 0.5 * sqrt(3.0) * lat->h;
    for (int ci = 0; ci < nc; ++ci) {
        const iso_leaf_t* c = &cells[ci];
        for (int s = 0; s < 2; ++s) {
            const int sheet = s == 0 ? 1 : -1;
            if (!(c->sheets & (1 << s))) continue;
            int mask = 0;
            double fmin_abs = DBL_MAX;
            for (int q = 0; q < 8; ++q) {
                const double fv = sheet * iso_lattice_corner(lat, c->x + (q & 1), c->y + ((q >> 1) & 1), c->z + ((q >> 2) & 1)) - o->iso;
                if (fv > 0) mask |= 1 << q;
                fmin_abs = fmin(fmin_abs, fabs(fv));
            }
            if (mask != 0 && mask != 0xFF) continue;
            r->candidates++;
            uint32_t bi = 0;
            iso_map_get(&brick_idx, iso_key3(c->x >> k, c->y >> k, c->z >> k), &bi);
            // all corners on one side: the trilinear interpolant stays there, and f differs from it
            // by at most sum_a h^2/8 max|f_aa| (tensor product linear interpolation error)
            if (fmin_abs > E[bi]) { r->by_trilinear++; continue; }
            if (fmin_abs > L[bi] * half_diag) { r->by_lipschitz++; continue; }

            iso_box_t b;
            const int idx[3] = {c->x, c->y, c->z};
            for (int a = 0; a < 3; ++a) { b.lo[a] = lat->origin[a] + idx[a] * lat->h; b.hi[a] = b.lo[a] + lat->h; }
            const int res = iso_hidden_rec(f, all, &b, sheet, o->iso, mask == 0xFF, 0, max_depth);
            if (res == 0) r->by_recursion++;
            else if (res == 1) {
                r->hidden++;
                iso_hidden_cell_t hc = {{0.5 * (b.lo[0] + b.hi[0]), 0.5 * (b.lo[1] + b.hi[1]), 0.5 * (b.lo[2] + b.hi[2])}, sheet};
                md_array_push(r->hidden_cells, hc, md_get_heap_allocator());
            }
            else r->unresolved++;
        }
    }
    free(L);
    free(E);
    iso_map_free(&brick_idx);
    free(all);
    r->seconds = iso_seconds_since(t0);
}

// ---------------------------------------------------------------------------------------------
// Sparse == dense: the octree may only discard cells without a crossing
// ---------------------------------------------------------------------------------------------

static void iso_mesh_checksum(const iso_mesh_t* m, double out[2]) {
    out[0] = out[1] = 0.0;
    const int nv = iso_mesh_nv(m);
    for (int i = 0; i < nv; ++i) { out[0] += m->pos[3 * i] + 2 * m->pos[3 * i + 1] + 3 * m->pos[3 * i + 2]; }
    out[1] = (double)(md_array_size(m->quad) / 4);
}

UTEST(gto_iso, sparse_surface_nets_equals_dense) {
    static const int pick[3] = {0, 2, 5};
    const double h = iso_full() ? 0.15 : 0.25;
    for (int ci = 0; ci < 3; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_field_t f;
        ASSERT_TRUE(iso_load(&f, tc->file, tc->which, tc->name));

        iso_octree_t o;
        iso_octree_build(&o, &f, tc->iso, h, ISO_BOUND_BEST, 0x1F);
        const int nl = (int)md_array_size(o.leaves);

        // dense lattice over the same root
        const int N = 1 << o.L;
        iso_leaf_t* dense = (iso_leaf_t*)malloc(sizeof(iso_leaf_t) * (size_t)N * N * N);
        size_t nd = 0;
        for (int z = 0; z < N; ++z) for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) dense[nd++] = (iso_leaf_t){x, y, z, 3};

        iso_lattice_t ls = {.f = &f, .h = o.h}, ld = {.f = &f, .h = o.h};
        memcpy(ls.origin, o.origin, sizeof(ls.origin));
        memcpy(ld.origin, o.origin, sizeof(ld.origin));
        iso_map_init(&ls.corner_map, (uint32_t)nl * 4);
        iso_map_init(&ld.corner_map, (uint32_t)((N + 1) * (N + 1) * (N + 1)));

        // the production path: GROUP bound only, octree stopped at 4^3 bricks
        iso_octree_t ob;
        iso_octree_build_ex(&ob, &f, tc->iso, h * 4, ISO_BOUND_SHELL, 0, 0.0);
        int nb = 0;
        iso_leaf_t* bcells = iso_bricks_to_cells(&ob, 2, &nb);
        iso_lattice_t lb = {.f = &f, .h = h};
        memcpy(lb.origin, ob.origin, sizeof(lb.origin));
        iso_map_init(&lb.corner_map, (uint32_t)nb * 2);

        md_tick_t t0 = md_tick_now();
        iso_mesh_t sp[2], de[2], br[2];
        for (int s = 0; s < 2; ++s) iso_surface_nets(&sp[s], &ls, o.leaves, nl, s == 0 ? 1 : -1, tc->iso);
        const double ts = iso_seconds_since(t0) + o.seconds;
        t0 = md_tick_now();
        for (int s = 0; s < 2; ++s) iso_surface_nets(&br[s], &lb, bcells, nb, s == 0 ? 1 : -1, tc->iso);
        const double tb = iso_seconds_since(t0) + ob.seconds;
        t0 = md_tick_now();
        for (int s = 0; s < 2; ++s) iso_surface_nets(&de[s], &ld, dense, (int)nd, s == 0 ? 1 : -1, tc->iso);
        const double td = iso_seconds_since(t0);

        printf("  %-13s h=%.3f lattice %d^3: leaf octree %ld evals %.2fs | bricks %ld evals %.2fs | dense %ld evals %.2fs\n",
            tc->name, h, N, ls.evals, ts, lb.evals, tb, ld.evals, td);
        for (int s = 0; s < 2; ++s) {
            double cs[2], cd[2], cb[2];
            iso_mesh_checksum(&sp[s], cs);
            iso_mesh_checksum(&de[s], cd);
            iso_mesh_checksum(&br[s], cb);
            printf("     sheet %c: leaf %d/%d, bricks %d/%d, dense %d/%d verts/quads, open edges %d/%d/%d\n",
                s == 0 ? '+' : '-', iso_mesh_nv(&sp[s]), (int)cs[1], iso_mesh_nv(&br[s]), (int)cb[1], iso_mesh_nv(&de[s]), (int)cd[1], sp[s].open_edges, br[s].open_edges, de[s].open_edges);
            EXPECT_EQ(iso_mesh_nv(&sp[s]), iso_mesh_nv(&de[s]));
            EXPECT_EQ(iso_mesh_nv(&br[s]), iso_mesh_nv(&de[s]));
            EXPECT_EQ(cs[1], cd[1]);
            EXPECT_EQ(cb[1], cd[1]);
            EXPECT_NEAR(cs[0], cd[0], 1e-9 * (1.0 + fabs(cd[0])));
            EXPECT_NEAR(cb[0], cd[0], 1e-9 * (1.0 + fabs(cd[0])));
            EXPECT_EQ(sp[s].open_edges, 0);
            EXPECT_EQ(br[s].open_edges, 0);
            iso_mesh_free(&sp[s]);
            iso_mesh_free(&de[s]);
            iso_mesh_free(&br[s]);
        }
        free(bcells);
        iso_map_free(&lb.corner_map);
        md_array_free(lb.corner_val, md_get_heap_allocator());
        iso_octree_free(&ob);
        free(dense);
        iso_map_free(&ls.corner_map); iso_map_free(&ld.corner_map);
        md_array_free(ls.corner_val, md_get_heap_allocator()); md_array_free(ld.corner_val, md_get_heap_allocator());
        iso_octree_free(&o);
        iso_field_free(&f);
    }
}

// ---------------------------------------------------------------------------------------------
// Which bound should drive the octree: work, kept leaves vs leaves that really cross
// ---------------------------------------------------------------------------------------------

static bool iso_report(void) { const char* e = getenv("GTO_ISO_REPORT"); return iso_full() || (e && e[0] == '1'); }

UTEST(gto_iso, octree_scheme_comparison) {
    if (!iso_report()) UTEST_SKIP("set GTO_ISO_REPORT=1 (about 45 s)");
    static const int pick[4] = {2, 4, 5, 6};
    typedef struct { int scheme; double t2_max; int brick; const char* label; } variant_t;
    static const variant_t variants[] = {
        {ISO_BOUND_PRIM,     0,     0, "prim"},
        {ISO_BOUND_SHELL,    0,     0, "group"},
        {ISO_BOUND_T2,       0,     0, "taylor2"},
        {ISO_BOUND_BEST,     0,     0, "best"},
        {ISO_SCHEME_CASCADE, 1e9,   0, "g>t2"},
        {ISO_BOUND_SHELL,    0,     1, "group b2"},
        {ISO_BOUND_SHELL,    0,     2, "group b4"},
        {ISO_BOUND_SHELL,    0,     3, "group b8"},
        {ISO_SCHEME_CASCADE, 1e9,   2, "g>t2 b4"},
    };
    const int num_variants = (int)(sizeof(variants) / sizeof(variants[0]));
    const double h = iso_full() ? 0.08 : 0.15;
    printf("\n  h = %.3f; bK = octree stops at bricks of K^3 cells, filled densely\n", h);
    printf("  %-13s %-10s %10s %12s %10s %10s %8s %8s %9s %7s %7s\n", "orbital", "bound", "enclosures", "group work", "cells", "crossing", "waste", "tree s", "corners", "eval s", "total s");
    for (int ci = 0; ci < 4; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_field_t f;
        ASSERT_TRUE(iso_load(&f, tc->file, tc->which, tc->name));
        int crossing_ref = -1;
        for (int k = 0; k < num_variants; ++k) {
            const variant_t* v = &variants[k];
            const uint32_t mask = v->scheme == ISO_BOUND_BEST ? 0x1F : v->scheme == ISO_SCHEME_CASCADE ? 0 : (1u << v->scheme);
            iso_octree_t o;
            iso_octree_build_ex(&o, &f, tc->iso, h * (1 << v->brick), v->scheme, mask, v->t2_max);
            int nc = 0;
            iso_leaf_t* cells = iso_bricks_to_cells(&o, v->brick, &nc);

            iso_lattice_t lat = {.f = &f, .h = h};
            memcpy(lat.origin, o.origin, sizeof(lat.origin));
            iso_map_init(&lat.corner_map, (uint32_t)nc * 2);
            int crossing = 0;
            md_tick_t tc0 = md_tick_now();
            for (int li = 0; li < nc; ++li) {
                const iso_leaf_t* c = &cells[li];
                bool any = false;
                for (int s = 0; s < 2 && !any; ++s) {
                    const int sheet = s == 0 ? 1 : -1;
                    int mask8 = 0;
                    for (int q = 0; q < 8; ++q) if (sheet * iso_lattice_corner(&lat, c->x + (q & 1), c->y + ((q >> 1) & 1), c->z + ((q >> 2) & 1)) - tc->iso > 0) mask8 |= 1 << q;
                    any = mask8 != 0 && mask8 != 0xFF;
                }
                crossing += any;
            }
            const double t_corner = iso_seconds_since(tc0);
            if (crossing_ref < 0) crossing_ref = crossing;
            EXPECT_EQ(crossing, crossing_ref);   // every variant must keep every crossing cell
            printf("  %-13s %-10s %10ld %12ld %10d %10d %7.1f%% %8.2f %9ld %7.2f %7.2f\n", tc->name, v->label, o.enclosures, o.group_work, nc, crossing,
                100.0 * (nc - crossing) / fmax(1, nc), o.seconds, lat.evals, t_corner, o.seconds + t_corner);
            free(cells);
            iso_map_free(&lat.corner_map);
            md_array_free(lat.corner_val, md_get_heap_allocator());
            iso_octree_free(&o);
        }
        iso_field_free(&f);
    }
}

// ---------------------------------------------------------------------------------------------
// The whole pipeline, with metrics at every stage
// ---------------------------------------------------------------------------------------------

static double iso_nearest_centre(const iso_field_t* f, const double p[3]) {
    double best = DBL_MAX;
    for (int i = 0; i < f->num_shells; ++i) best = fmin(best, iso_dist2(p, f->shells[i].c));
    return sqrt(best);
}

static void iso_print_metrics(const char* stage, const iso_metrics_t* m) {
    printf("     %-8s V %6d T %6d chi %3d comp %2d bnd %d nm %d | min angle %5.1f p1 %5.1f mean %5.1f | dev max %.4f rms %.5f h | orient %.3f\n",
        stage, m->nv, m->nt, m->euler, m->components, m->boundary_edges, m->nonmanifold_edges,
        m->min_angle, m->p1_angle, m->mean_min_angle, m->dev_max, m->dev_rms, m->orient_ok);
}

UTEST(gto_iso, mesh_pipeline) {
    if (!iso_report() && !getenv("GTO_ISO_ONLY")) UTEST_SKIP("set GTO_ISO_REPORT=1 (about 20 s)");
    const double h = getenv("GTO_ISO_H") ? atof(getenv("GTO_ISO_H")) : iso_full() ? 0.08 : 0.15;
    const int brick = 2;   // octree stops at 4^3 cell bricks, GROUP bound only (see octree_scheme_comparison)
    const int num_cases = iso_full() ? ISO_NUM_CASES : 6;
    const char* only = getenv("GTO_ISO_ONLY");
    FILE* csv = iso_open_out("gto_iso_pipeline.csv", "w");
    if (csv) fprintf(csv, "orbital,sheet,stage,nv,nt,euler,components,min_angle,p1_angle,mean_min_angle,dev_max_h,dev_rms_h,seconds\n");

    for (int ci = 0; ci < (only ? ISO_NUM_CASES : num_cases); ++ci) {
        const iso_case_t* tc = &iso_cases[ci];
        if (only && !strstr(only, tc->name)) continue;
        iso_field_t f;
        ASSERT_TRUE(iso_load(&f, tc->file, tc->which, tc->name));

        iso_octree_t o;
        iso_octree_build_ex(&o, &f, tc->iso, h * (1 << brick), ISO_BOUND_SHELL, 0, 0.0);
        int nc = 0;
        iso_leaf_t* cells = iso_bricks_to_cells(&o, brick, &nc);
        const int Lfine = o.L + brick;
        const long dense_evals = (long)((1 << Lfine) + 1) * ((1 << Lfine) + 1) * ((1 << Lfine) + 1);

        iso_lattice_t lat = {.f = &f, .h = h};
        memcpy(lat.origin, o.origin, sizeof(lat.origin));
        iso_map_init(&lat.corner_map, (uint32_t)nc * 2);

        long tree_nodes = 0;
        for (int l = 0; l <= o.L; ++l) tree_nodes += o.visited[l];
        printf("\n  %s  iso %.3f  h %.3f  root %.1f bohr = %d^3 cells  groups %d\n", tc->name, tc->iso, h, o.h * (1 << o.L), 1 << Lfine, f.num_shells);
        printf("     octree: %ld nodes, %d bricks of %d^3 kept, %.3fs\n", tree_nodes, (int)md_array_size(o.leaves), 1 << brick, o.seconds);

        md_tick_t t0 = md_tick_now();
        iso_mesh_t m[2];
        for (int s = 0; s < 2; ++s) iso_surface_nets(&m[s], &lat, cells, nc, s == 0 ? 1 : -1, tc->iso);
        const double t_sn = iso_seconds_since(t0);
        printf("     corner evaluations %ld (dense lattice %ld, %.0fx fewer), surface nets incl. evaluation %.3fs\n", lat.evals, dense_evals, (double)dense_evals / fmax(1, lat.evals), t_sn);

        iso_guarantee_t g;
        iso_guarantee(&g, &o, brick, &lat, cells, nc, 4);
        printf("     guarantee: %ld kept cells without a crossing -> %ld trilinear, %ld Lipschitz, %ld recursion, %ld HIDDEN, %ld unresolved at h/16 (%.3fs)\n",
            g.candidates, g.by_trilinear, g.by_lipschitz, g.by_recursion, g.hidden, g.unresolved, g.seconds);
        for (int i = 0; i < (int)md_array_size(g.hidden_cells) && i < 4; ++i) {
            const iso_hidden_cell_t* hc = &g.hidden_cells[i];
            // depth of the feature: how far past the iso value it reaches, as a distance (excess / |grad f|)
            double peak = -DBL_MAX, pk[3] = {0};
            for (int z = 0; z <= 16; ++z) for (int y = 0; y <= 16; ++y) for (int x = 0; x <= 16; ++x) {
                const double p[3] = {hc->p[0] + (x / 16.0 - 0.5) * h, hc->p[1] + (y / 16.0 - 0.5) * h, hc->p[2] + (z / 16.0 - 0.5) * h};
                const double v = hc->sheet * iso_value(&f, p);
                if (v > peak) { peak = v; memcpy(pk, p, sizeof(p)); }
            }
            double v, gr[3];
            iso_eval(&f, NULL, 0, pk, 1, &v, gr, NULL);
            const double depth = (peak - tc->iso) / fmax(sqrt(gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2]), 1e-12);
            printf("        hidden %c at (%.2f %.2f %.2f), %.3f bohr from the nearest atom: peak %.4f, depth ~%.3f h\n",
                hc->sheet > 0 ? '+' : '-', hc->p[0], hc->p[1], hc->p[2], iso_nearest_centre(&f, hc->p), peak, depth / h);
        }
        md_array_free(g.hidden_cells, md_get_heap_allocator());

        FILE* obj = NULL;
        if (iso_out_dir()) {
            char name[128];
            snprintf(name, sizeof(name), "gto_iso_%s.obj", tc->name);
            obj = iso_open_out(name, "w");
        }
        int vbase = 0;
        for (int s = 0; s < 2; ++s) {
            iso_mesh_t* mm = &m[s];
            if (iso_mesh_nv(mm) == 0) continue;
            const char sc = s == 0 ? '+' : '-';
            printf("   sheet %c (open edges %d)\n", sc, mm->open_edges);
            EXPECT_EQ(mm->open_edges, 0);
            iso_metrics_t met;
#define ISO_STAGE(stage_name, secs) do { iso_mesh_metrics(&met, mm, &f); iso_print_metrics(stage_name, &met); \
            if (csv) fprintf(csv, "%s,%c,%s,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.6f,%.6f,%.4f\n", tc->name, sc, stage_name, met.nv, met.nt, met.euler, met.components, met.min_angle, met.p1_angle, met.mean_min_angle, met.dev_max, met.dev_rms, (double)(secs)); } while (0)

            iso_mesh_triangulate(mm);
            ISO_STAGE("raw", 0.0);

            t0 = md_tick_now();
            iso_mesh_project(mm, &f, 8);
            iso_mesh_triangulate(mm);
            const double t_proj = iso_seconds_since(t0);
            ISO_STAGE("project", t_proj);
            const iso_metrics_t m_proj = met;

            t0 = md_tick_now();
            iso_mesh_relax(mm, &f, 3, 0.5);
            const double t_relax = iso_seconds_since(t0);
            ISO_STAGE("relax", t_relax);
            EXPECT_EQ(met.boundary_edges, 0);
            EXPECT_EQ(met.euler, m_proj.euler);   // relaxation must not change topology

            t0 = md_tick_now();
            int splits = 0;
            for (int pass = 0; pass < 2; ++pass) splits += iso_mesh_refine(mm, &f, 0.02 * h, 20.0);
            const double t_ref = iso_seconds_since(t0);
            ISO_STAGE("refine", t_ref);
            EXPECT_EQ(met.boundary_edges, 0);
            EXPECT_EQ(met.euler, m_proj.euler);
            t0 = md_tick_now();
            iso_mesh_relax(mm, &f, 2, 0.5);
            ISO_STAGE("relax2", iso_seconds_since(t0));
            EXPECT_EQ(met.euler, m_proj.euler);
            EXPECT_GT(met.orient_ok, 0.99);
            printf("     times: project %.3fs  relax %.3fs  refine %.3fs (%d splits); worst deviation at %.3f bohr from an atom, |grad f| %.3g there\n",
                t_proj, t_relax, t_ref, splits, iso_nearest_centre(&f, met.dev_where), met.dev_grad);
#undef ISO_STAGE
            if (obj) {
                iso_mesh_write_obj(obj, mm, vbase, s == 0 ? "positive" : "negative");
                vbase += iso_mesh_nv(mm);
            }
        }
        if (obj) fclose(obj);
        for (int s = 0; s < 2; ++s) iso_mesh_free(&m[s]);
        free(cells);
        iso_map_free(&lat.corner_map);
        md_array_free(lat.corner_val, md_get_heap_allocator());
        iso_octree_free(&o);
        iso_field_free(&f);
    }
    if (csv) fclose(csv);
}

// ---------------------------------------------------------------------------------------------
// What the eye gets: dense raycast (viamd's exact mode) vs the gridless mesh, on the same lattice
// ---------------------------------------------------------------------------------------------
//
// Dense: the hit point is the exact crossing of the TRILINEAR interpolant of the lattice, and the
// normal is a central difference one lattice step apart on that interpolant (isosurface.frag,
// gradient_view). Hardware filter quantization is not modelled, so this is the raycaster's best case.
// Mesh: vertices projected onto the true surface, normals analytic at the vertices and
// interpolated across triangles, flat triangles in between.
// Both are measured against the analytic surface: position error = distance to the true surface
// (Newton projection), normal error = angle to the true normal there.

static double iso_lat_trilinear(iso_lattice_t* lat, const double p[3]) {
    int c[3];
    double u[3];
    for (int a = 0; a < 3; ++a) {
        const double t = (p[a] - lat->origin[a]) / lat->h;
        c[a] = (int)floor(t);
        u[a] = t - c[a];
    }
    double v = 0.0;
    for (int k = 0; k < 8; ++k) {
        const int b[3] = {k & 1, (k >> 1) & 1, (k >> 2) & 1};
        const double w = (b[0] ? u[0] : 1 - u[0]) * (b[1] ? u[1] : 1 - u[1]) * (b[2] ? u[2] : 1 - u[2]);
        v += w * iso_lattice_corner(lat, c[0] + b[0], c[1] + b[1], c[2] + b[2]);
    }
    return v;
}

typedef struct iso_err_t { double* pos; double* ang; int n; } iso_err_t;

static void iso_err_push(iso_err_t* e, double pos, double ang) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_array_push(e->pos, pos, heap);
    md_array_push(e->ang, ang, heap);
    e->n++;
}

static void iso_err_print(const char* label, iso_err_t* e, double h) {
    double s2 = 0.0, sa = 0.0;
    for (int i = 0; i < e->n; ++i) { s2 += e->pos[i] * e->pos[i]; sa += e->ang[i]; }
    const double rms = sqrt(s2 / fmax(1, e->n));
    const double p99 = iso_percentile(e->pos, e->n, 0.99), pmax = iso_percentile(e->pos, e->n, 1.0);
    const double a_mean = sa / fmax(1, e->n);
    const double a99 = iso_percentile(e->ang, e->n, 0.99), amax = iso_percentile(e->ang, e->n, 1.0);
    printf("     %-15s position rms %.5f p99 %.5f max %.5f bohr (%.4f / %.4f / %.4f h) | normal mean %5.2f p99 %5.2f max %6.2f deg\n",
        label, rms, p99, pmax, rms / h, p99 / h, pmax / h, a_mean, a99, amax);
    md_array_free(e->pos, md_get_heap_allocator());
    md_array_free(e->ang, md_get_heap_allocator());
    memset(e, 0, sizeof(*e));
}

static double iso_angle_deg(const double a[3], const double b[3]) {
    const double la = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]), lb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    double c = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / fmax(la * lb, 1e-300);
    c = fmax(-1.0, fmin(1.0, c));
    return acos(c) * 180.0 / 3.14159265358979;
}

static void iso_mesh_sample_errors(iso_err_t* e, const iso_mesh_t* m, const iso_field_t* f, iso_rng_t* rng, int per_tri) {
    const int nt = iso_mesh_nt(m);
    for (int t = 0; t < nt; ++t) {
        const uint32_t* v = m->tri + 3 * t;
        for (int k = 0; k < per_tri; ++k) {
            double b0 = iso_rng_f64(rng), b1 = iso_rng_f64(rng);
            if (b0 + b1 > 1.0) { b0 = 1.0 - b0; b1 = 1.0 - b1; }
            const double b2 = 1.0 - b0 - b1;
            double p[3], n[3], q[3], nt_true[3];
            for (int a = 0; a < 3; ++a) {
                p[a] = b0 * m->pos[3 * v[0] + a] + b1 * m->pos[3 * v[1] + a] + b2 * m->pos[3 * v[2] + a];
                n[a] = b0 * m->nrm[3 * v[0] + a] + b1 * m->nrm[3 * v[1] + a] + b2 * m->nrm[3 * v[2] + a];
                q[a] = p[a];
            }
            iso_project_point(f, q, m->sheet, m->iso, 0.5 * m->h, 12, nt_true);
            iso_err_push(e, sqrt(iso_dist2(p, q)), iso_angle_deg(n, nt_true));
        }
    }
}

UTEST(gto_iso, visual_dense_vs_mesh) {
    if (!iso_report()) UTEST_SKIP("set GTO_ISO_REPORT=1");
    static const int pick[3] = {2, 7, 5};   // amide_homo, mol_lumo, amide_random
    static const double hs[3] = {0.15, 0.1, 0.05};
    for (int ci = 0; ci < 3; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_field_t f;
        ASSERT_TRUE(iso_load(&f, tc->file, tc->which, tc->name));
        for (int hi = 0; hi < 3; ++hi) {
            const double h = hs[hi];
            iso_octree_t o;
            iso_octree_build_ex(&o, &f, tc->iso, 4.0 * h, ISO_BOUND_SHELL, 0, 0.0);
            int nc = 0;
            iso_leaf_t* cells = iso_bricks_to_cells(&o, 2, &nc);
            iso_lattice_t lat = {.f = &f, .h = h};
            memcpy(lat.origin, o.origin, sizeof(lat.origin));
            iso_map_init(&lat.corner_map, (uint32_t)nc * 2);
            printf("\n  %s  h %.3f bohr\n", tc->name, h);

            iso_rng_t rng = {4242};
            iso_err_t dense = {0}, mesh = {0}, mesh_ref = {0};
            int tris = 0, tris_ref = 0;
            for (int s = 0; s < 2; ++s) {
                const int sheet = s == 0 ? 1 : -1;
                // dense: points on the trilinear isosurface, 8 per crossing cell
                for (int li = 0; li < nc; ++li) {
                    const iso_leaf_t* c = &cells[li];
                    if (!(c->sheets & (1 << s))) continue;
                    int mask = 0;
                    for (int k = 0; k < 8; ++k) if (sheet * iso_lattice_corner(&lat, c->x + (k & 1), c->y + ((k >> 1) & 1), c->z + ((k >> 2) & 1)) - tc->iso > 0) mask |= 1 << k;
                    if (mask == 0 || mask == 0xFF) continue;
                    const double lo[3] = {lat.origin[0] + c->x * h, lat.origin[1] + c->y * h, lat.origin[2] + c->z * h};
                    for (int k = 0; k < 8; ++k) {
                        double A[3], B[3], fa = 0, fb = 0;
                        int tries = 0;
                        do {
                            for (int a = 0; a < 3; ++a) { A[a] = lo[a] + h * iso_rng_f64(&rng); B[a] = lo[a] + h * iso_rng_f64(&rng); }
                            fa = sheet * iso_lat_trilinear(&lat, A) - tc->iso;
                            fb = sheet * iso_lat_trilinear(&lat, B) - tc->iso;
                        } while ((fa > 0) == (fb > 0) && ++tries < 64);
                        if ((fa > 0) == (fb > 0)) continue;
                        for (int it = 0; it < 50; ++it) {
                            double M[3];
                            for (int a = 0; a < 3; ++a) M[a] = 0.5 * (A[a] + B[a]);
                            const double fm = sheet * iso_lat_trilinear(&lat, M) - tc->iso;
                            if ((fm > 0) == (fa > 0)) { memcpy(A, M, sizeof(M)); fa = fm; } else memcpy(B, M, sizeof(M));
                        }
                        double p[3], q[3], g[3], nt_true[3];
                        for (int a = 0; a < 3; ++a) p[a] = q[a] = 0.5 * (A[a] + B[a]);
                        for (int a = 0; a < 3; ++a) {   // the raycaster's normal: central difference one lattice step apart
                            double pp[3] = {p[0], p[1], p[2]}, pm[3] = {p[0], p[1], p[2]};
                            pp[a] += h; pm[a] -= h;
                            g[a] = -sheet * (iso_lat_trilinear(&lat, pp) - iso_lat_trilinear(&lat, pm));
                        }
                        iso_project_point(&f, q, sheet, tc->iso, 0.5 * h, 12, nt_true);
                        iso_err_push(&dense, sqrt(iso_dist2(p, q)), iso_angle_deg(g, nt_true));
                    }
                }
                // mesh: the GPU path's stages (project, 3 relaxations), then optionally refined
                iso_mesh_t m;
                iso_surface_nets(&m, &lat, cells, nc, sheet, tc->iso);
                if (iso_mesh_nv(&m) == 0) { iso_mesh_free(&m); continue; }
                iso_mesh_triangulate(&m);
                iso_mesh_project(&m, &f, 8);
                iso_mesh_triangulate(&m);
                iso_mesh_relax(&m, &f, 3, 0.5);
                tris += iso_mesh_nt(&m);
                iso_mesh_sample_errors(&mesh, &m, &f, &rng, 4);
                for (int pass = 0; pass < 2; ++pass) iso_mesh_refine(&m, &f, 0.02 * h, 20.0);
                iso_mesh_relax(&m, &f, 2, 0.5);
                tris_ref += iso_mesh_nt(&m);
                iso_mesh_sample_errors(&mesh_ref, &m, &f, &rng, 4);
                iso_mesh_free(&m);
            }
            iso_err_print("dense raycast", &dense, h);
            char l1[64], l2[64];
            snprintf(l1, sizeof(l1), "mesh (%dk tri)", tris / 1000);
            snprintf(l2, sizeof(l2), "refined (%dk)", tris_ref / 1000);
            iso_err_print(l1, &mesh, h);
            iso_err_print(l2, &mesh_ref, h);

            free(cells);
            iso_map_free(&lat.corner_map);
            md_array_free(lat.corner_val, md_get_heap_allocator());
            iso_octree_free(&o);
        }
        iso_field_free(&f);
    }
}

// =============================================================================================
// Public API (md_gto_iso.h)
// =============================================================================================

typedef struct iso_api_case_t {
    vlx_test_t          t;
    md_gto_basis_t      basis;
    float*              xyz;
    double*             coeff;
    md_gto_iso_field_t* field;
} iso_api_case_t;

static bool iso_api_load(iso_api_case_t* c, const iso_case_t* tc) {
    memset(c, 0, sizeof(*c));
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", tc->file);
    if (!vlx_test_load(&c->t, str_from_cstr(path), MEGABYTES(64))) return false;
    const size_t num_atoms = c->t.sys.atom.count;
    c->xyz = (float*)malloc(sizeof(float) * 3 * num_atoms);
    if (vlx_test_atom_xyz_bohr(c->xyz, 3 * num_atoms, &c->t) != num_atoms || !qm_test_basis(&c->basis, &c->t)) return false;
    const size_t num_ao = md_gto_basis_num_ao(&c->basis);
    c->coeff = (double*)calloc(num_ao, sizeof(double));
    const size_t lumo = vlx_test_lumo_idx(&c->t, STR_LIT("orbital/alpha/occupation"));
    const size_t idx = tc->which >= 0 ? (size_t)tc->which : (tc->which == ISO_MO_LUMO ? lumo : lumo - 1);
    if (qm_test_row(c->coeff, num_ao, &c->t, STR_LIT("orbital/alpha/coefficient"), idx) != num_ao) return false;
    if (tc->which == ISO_MO_RANDOM) {
        double n0 = 0.0, n1 = 0.0;
        iso_rng_t rng = {12345};
        for (size_t i = 0; i < num_ao; ++i) n0 += c->coeff[i] * c->coeff[i];
        for (size_t i = 0; i < num_ao; ++i) { c->coeff[i] = iso_rng_normal(&rng); n1 += c->coeff[i] * c->coeff[i]; }
        for (size_t i = 0; i < num_ao; ++i) c->coeff[i] *= sqrt(n0 / n1);
    }
    c->field = md_gto_iso_field_create(&c->basis, c->xyz, 0, c->coeff, md_get_heap_allocator());
    return c->field != NULL;
}

static void iso_api_free(iso_api_case_t* c) {
    md_gto_iso_field_destroy(c->field);
    free(c->xyz);
    free(c->coeff);
    qm_test_free(&c->t);
    memset(c, 0, sizeof(*c));
}

// Dense Surface Nets over the cache's whole lattice, double precision: what the cache must reproduce.
static void iso_dense_counts(const md_gto_iso_field_t* field, const md_gto_iso_cache_info_t* info, double iso, int out_nv[2], int out_nq[2]) {
    const int N = (int)info->cells_per_axis;
    iso_leaf_t* dense = (iso_leaf_t*)malloc(sizeof(iso_leaf_t) * (size_t)N * N * N);
    size_t nd = 0;
    for (int z = 0; z < N; ++z) for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x) dense[nd++] = (iso_leaf_t){x, y, z, 3};
    iso_lattice_t lat = {.f = &field->f, .h = info->cell_size};
    memcpy(lat.origin, info->origin, sizeof(lat.origin));
    iso_map_init(&lat.corner_map, (uint32_t)((N + 1) * (N + 1) * (N + 1)));
    for (int s = 0; s < 2; ++s) {
        iso_mesh_t m;
        iso_surface_nets(&m, &lat, dense, (int)nd, s == 0 ? 1 : -1, iso);
        out_nv[s] = iso_mesh_nv(&m);
        out_nq[s] = (int)(md_array_size(m.quad) / 4);
        iso_mesh_free(&m);
    }
    iso_map_free(&lat.corner_map);
    md_array_free(lat.corner_val, md_get_heap_allocator());
    free(dense);
}

// An md_gto_iso_mesh_t sheet as an iso_mesh_t, for the metrics.
static void iso_api_sheet(iso_mesh_t* m, const md_gto_iso_mesh_t* am, int s, double iso, double h) {
    md_allocator_i* heap = md_get_heap_allocator();
    memset(m, 0, sizeof(*m));
    m->sheet = s == 0 ? 1 : -1; m->iso = iso; m->h = h;
    const uint32_t vo = am->sheet_vertex_offset[s];
    for (uint32_t i = 0; i < am->sheet_vertex_count[s]; ++i) {
        const double p[3] = {am->position[3 * (vo + i)], am->position[3 * (vo + i) + 1], am->position[3 * (vo + i) + 2]};
        const double n[3] = {am->normal[3 * (vo + i)], am->normal[3 * (vo + i) + 1], am->normal[3 * (vo + i) + 2]};
        md_array_push_array(m->pos, p, 3, heap);
        md_array_push_array(m->nrm, n, 3, heap);
    }
    for (uint32_t i = 0; i < am->sheet_index_count[s]; ++i) {
        const uint32_t v = am->index[am->sheet_index_offset[s] + i] - vo;
        md_array_push(m->tri, v, heap);
    }
}

UTEST(gto_iso_api, cache_serves_its_interval_exactly) {
    static const int pick[3] = {0, 2, 5};
    for (int ci = 0; ci < 3; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_api_case_t c;
        ASSERT_TRUE(iso_api_load(&c, tc));

        // coarse, so the dense reference stays quick; exactness does not depend on the spacing
        md_gto_iso_cache_desc_t cd = {.cell_size = 0.35, .iso_min = 0.02, .iso_max = 0.1};
        md_tick_t t0 = md_tick_now();
        md_gto_iso_cache_t* cache = md_gto_iso_cache_create(c.field, &cd, md_get_heap_allocator());
        const double t_cache = iso_seconds_since(t0);
        ASSERT_TRUE(cache != NULL);
        md_gto_iso_cache_info_t info;
        md_gto_iso_cache_info(&info, cache);
        printf("  %-13s cache [%.2f, %.2f] at h %.2f: %u bricks of %u^3 cells, %.1f MB (a dense float volume: %.1f MB), %.3fs\n",
            tc->name, info.iso_min, info.iso_max, info.cell_size, info.num_bricks, info.cells_per_axis,
            info.memory_bytes / 1048576.0, pow(info.cells_per_axis + 1, 3) * 4 / 1048576.0, t_cache);

        EXPECT_FALSE(md_gto_iso_cache_covers(cache, 0.01));
        EXPECT_FALSE(md_gto_iso_cache_covers(cache, 0.2));
        md_gto_iso_mesh_t bad;
        EXPECT_FALSE(md_gto_iso_mesh_extract(&bad, cache, &(md_gto_iso_mesh_desc_t){.iso = 0.2}, md_get_heap_allocator()));

        static const double isos[3] = {0.02, 0.05, 0.1};
        for (int k = 0; k < 3; ++k) {
            md_gto_iso_mesh_t raw;
            t0 = md_tick_now();
            ASSERT_TRUE(md_gto_iso_mesh_extract(&raw, cache, &(md_gto_iso_mesh_desc_t){.iso = isos[k]}, md_get_heap_allocator()));
            const double t_raw = iso_seconds_since(t0);
            int dnv[2], dnq[2];
            iso_dense_counts(c.field, &info, isos[k], dnv, dnq);
            printf("     iso %.2f raw: %u + %u vertices (dense lattice %d + %d), %u + %u triangles (dense %d + %d), %.3fs\n", isos[k],
                raw.sheet_vertex_count[0], raw.sheet_vertex_count[1], dnv[0], dnv[1],
                raw.sheet_index_count[0] / 3, raw.sheet_index_count[1] / 3, 2 * dnq[0], 2 * dnq[1], t_raw);
            for (int s = 0; s < 2; ++s) {
                EXPECT_EQ((int)raw.sheet_vertex_count[s], dnv[s]);
                EXPECT_EQ((int)raw.sheet_index_count[s], 6 * dnq[s]);
            }
            md_gto_iso_mesh_free(&raw);
        }
        md_gto_iso_cache_destroy(cache);
        iso_api_free(&c);
    }
}

UTEST(gto_iso_api, display_meshes_are_closed_and_on_the_surface) {
    static const int pick[4] = {0, 2, 3, 5};
    const double h = 0.15, iso = 0.05;
    for (int ci = 0; ci < 4; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_api_case_t c;
        ASSERT_TRUE(iso_api_load(&c, tc));
        md_tick_t t0 = md_tick_now();
        md_gto_iso_cache_t* cache = md_gto_iso_cache_create(c.field, &(md_gto_iso_cache_desc_t){.cell_size = h, .iso_min = iso}, md_get_heap_allocator());
        const double t_cache = iso_seconds_since(t0);
        ASSERT_TRUE(cache != NULL);

        md_gto_iso_mesh_t raw, q;
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_iso_mesh_extract(&raw, cache, &(md_gto_iso_mesh_desc_t){.iso = iso}, md_get_heap_allocator()));
        const double t_raw = iso_seconds_since(t0);
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_iso_mesh_extract(&q, cache, &(md_gto_iso_mesh_desc_t){.iso = iso, .project_iterations = 8, .relax_iterations = 3}, md_get_heap_allocator()));
        const double t_q = iso_seconds_since(t0);
        printf("  %-13s h %.2f: cache %.3fs, raw mesh %.3fs, projected + relaxed %.3fs, %u vertices, %u triangles\n",
            tc->name, h, t_cache, t_raw, t_q, q.vertex_count, q.index_count / 3);
        EXPECT_EQ(raw.vertex_count, q.vertex_count);
        EXPECT_EQ(raw.index_count, q.index_count);

        for (int s = 0; s < 2; ++s) {
            if (!q.sheet_vertex_count[s]) continue;
            iso_mesh_t mr, mq;
            iso_api_sheet(&mr, &raw, s, iso, h);
            iso_api_sheet(&mq, &q, s, iso, h);
            iso_metrics_t er, eq;
            iso_mesh_metrics(&er, &mr, &c.field->f);
            iso_mesh_metrics(&eq, &mq, &c.field->f);
            printf("     sheet %c raw      ", s == 0 ? '+' : '-'); iso_print_metrics("", &er);
            printf("     sheet %c display  ", s == 0 ? '+' : '-'); iso_print_metrics("", &eq);
            EXPECT_EQ(eq.boundary_edges, 0);
            EXPECT_EQ(eq.euler, er.euler);
            EXPECT_GT(er.orient_ok, 0.99);
            EXPECT_GT(eq.orient_ok, 0.99);
            EXPECT_LT(eq.dev_rms, 0.05);
            iso_mesh_free(&mr);
            iso_mesh_free(&mq);
        }
        md_gto_iso_mesh_free(&raw);
        md_gto_iso_mesh_free(&q);
        md_gto_iso_cache_destroy(cache);
        iso_api_free(&c);
    }
}

// =============================================================================================
// Density fields: rho = sum_t w_t a_t b_t
// =============================================================================================
//
// The densities viamd shows, built the two ways the API takes them:
//   products  what the file already has: occupied MOs with their occupations (SCF), NTO pairs
//             (transition density), NTOs squared with their weights (attachment, detachment,
//             difference)
//   matrix    the AO density matrix md_vlx publishes (orbital/total/density,
//             vlx/rsp/transition_density/*), split by the field itself
// and checked against each other and against the orbital evaluator.

enum { ISO_D_SCF, ISO_D_ATTACH, ISO_D_DETACH, ISO_D_DIFF, ISO_D_TRANS, ISO_D_COUNT };
static const char* iso_d_name[ISO_D_COUNT] = {"scf", "attach", "detach", "diff", "trans"};

typedef struct iso_dens_src_t {
    vlx_test_t     t;
    md_gto_basis_t basis;
    float*         xyz;
    size_t         num_atoms, num_ao, num_mo;
    double*        C;          // [num_mo][num_ao]
    double*        occ;        // alpha + beta occupation per MO
    size_t         num_nto;    // significant NTO pairs of state 1
    double*        lambda;     // [num_nto]
    double*        hole;       // [num_nto][num_ao]
    double*        particle;   // [num_nto][num_ao]
} iso_dens_src_t;

static void iso_dens_src_free(iso_dens_src_t* s) {
    free(s->xyz); free(s->C); free(s->occ); free(s->lambda); free(s->hole); free(s->particle);
    qm_test_free(&s->t);
    memset(s, 0, sizeof(*s));
}

static bool iso_dens_src_load_ex(iso_dens_src_t* s, const char* file, bool need_nto) {
    memset(s, 0, sizeof(*s));
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    if (!vlx_test_load(&s->t, str_from_cstr(path), MEGABYTES(512))) return false;
    s->num_atoms = s->t.sys.atom.count;
    s->xyz = (float*)malloc(sizeof(float) * 3 * s->num_atoms);
    if (vlx_test_atom_xyz_bohr(s->xyz, 3 * s->num_atoms, &s->t) != s->num_atoms || !qm_test_basis(&s->basis, &s->t)) return false;
    s->num_ao = md_gto_basis_num_ao(&s->basis);
    const md_attribute_t* ca = qm_test_attr(&s->t, STR_LIT("orbital/alpha/coefficient"));
    if (!ca) return false;
    s->num_mo = ca->format.shape[0];
    s->C = (double*)malloc(sizeof(double) * s->num_mo * s->num_ao);
    for (size_t m = 0; m < s->num_mo; ++m) {
        if (qm_test_row(s->C + m * s->num_ao, s->num_ao, &s->t, STR_LIT("orbital/alpha/coefficient"), m) != s->num_ao) return false;
    }
    s->occ = (double*)calloc(s->num_mo, sizeof(double));
    double* ob = (double*)calloc(s->num_mo, sizeof(double));
    qm_test_series(s->occ, s->num_mo, &s->t, STR_LIT("orbital/alpha/occupation"));
    qm_test_series(ob, s->num_mo, &s->t, STR_LIT("orbital/beta/occupation"));
    for (size_t m = 0; m < s->num_mo; ++m) s->occ[m] += ob[m];
    free(ob);

    // NTOs of the first excited state, the pairs with a weight above 1e-6
    const md_attribute_t* la = qm_test_attr(&s->t, STR_LIT("vlx/rsp/nto/lambda"));
    const md_attribute_t* pa = qm_test_attr(&s->t, STR_LIT("vlx/rsp/nto/particle/coefficient"));
    const md_attribute_t* ha = qm_test_attr(&s->t, STR_LIT("vlx/rsp/nto/hole/coefficient"));
    if (!la || !pa || !ha) return !need_nto;
    const size_t lmax = la->format.shape[1];
    double* lam = (double*)calloc(lmax, sizeof(double));
    if (md_attribute_extract_f64(lam, lmax, la, md_attribute_slice_1(0), md_unit_none()) != lmax) return false;
    s->lambda = (double*)calloc(lmax, sizeof(double));
    s->hole = (double*)calloc(lmax * s->num_ao, sizeof(double));
    s->particle = (double*)calloc(lmax * s->num_ao, sizeof(double));
    for (size_t k = 0; k < lmax; ++k) {
        if (!(lam[k] > 1.0e-6)) continue;
        const size_t n = s->num_nto++;
        s->lambda[n] = lam[k];
        if (md_attribute_extract_f64(s->hole + n * s->num_ao, s->num_ao, ha, md_attribute_slice_2(0, (uint32_t)k), md_unit_none()) != s->num_ao) return false;
        if (md_attribute_extract_f64(s->particle + n * s->num_ao, s->num_ao, pa, md_attribute_slice_2(0, (uint32_t)k), md_unit_none()) != s->num_ao) return false;
    }
    free(lam);
    return s->num_nto > 0 || !need_nto;
}

static bool iso_dens_src_load(iso_dens_src_t* s, const char* file) { return iso_dens_src_load_ex(s, file, true); }

// The products of one density kind: count terms, rows left/right (right NULL: squares).
typedef struct iso_dens_prod_t {
    size_t  count;
    double* left;
    double* right;
    double* w;
} iso_dens_prod_t;

static void iso_dens_prod_free(iso_dens_prod_t* p) { free(p->left); free(p->right); free(p->w); memset(p, 0, sizeof(*p)); }

static void iso_dens_products(iso_dens_prod_t* p, const iso_dens_src_t* s, int kind) {
    memset(p, 0, sizeof(*p));
    const size_t nao = s->num_ao;
    const size_t cap = s->num_mo + 2 * s->num_nto + 1;
    p->left = (double*)calloc(cap * nao, sizeof(double));
    p->w = (double*)calloc(cap, sizeof(double));
    switch (kind) {
    case ISO_D_SCF: {
        // BENCHMARK: GTO_ISO_SCF_EXTRA=N adds the first N virtual orbitals at weight 0.01 (a set of
        // as many more functions, for the wide kernels on a small system)
        int extra = getenv("GTO_ISO_SCF_EXTRA") ? atoi(getenv("GTO_ISO_SCF_EXTRA")) : 0;
        for (size_t m = 0; m < s->num_mo; ++m) {
            double w = s->occ[m];
            if (w == 0.0) {
                if (extra <= 0) continue;
                extra--;
                w = 0.01;
            }
            memcpy(p->left + p->count * nao, s->C + m * nao, sizeof(double) * nao);
            p->w[p->count++] = w;
        }
        break;
    }
    case ISO_D_ATTACH:
    case ISO_D_DETACH:
    case ISO_D_DIFF:
        for (size_t k = 0; k < s->num_nto; ++k) {
            if (kind != ISO_D_DETACH) {
                memcpy(p->left + p->count * nao, s->particle + k * nao, sizeof(double) * nao);
                p->w[p->count++] = s->lambda[k];
            }
            if (kind != ISO_D_ATTACH) {
                memcpy(p->left + p->count * nao, s->hole + k * nao, sizeof(double) * nao);
                p->w[p->count++] = kind == ISO_D_DIFF ? -s->lambda[k] : s->lambda[k];
            }
        }
        break;
    case ISO_D_TRANS:
        // rho_0n = sum_k sqrt(lambda_k) h_k p_k (up to the spin factor, which only scales it)
        p->right = (double*)calloc(cap * nao, sizeof(double));
        for (size_t k = 0; k < s->num_nto; ++k) {
            memcpy(p->left + p->count * nao, s->hole + k * nao, sizeof(double) * nao);
            memcpy(p->right + p->count * nao, s->particle + k * nao, sizeof(double) * nao);
            p->w[p->count++] = sqrt(s->lambda[k]);
        }
        break;
    }
}

// The same density as an AO matrix: sum_k w_k a_k b_k^T (non-symmetric for a transition density).
static double* iso_dens_matrix(const iso_dens_prod_t* p, size_t nao) {
    double* D = (double*)calloc(nao * nao, sizeof(double));
    for (size_t k = 0; k < p->count; ++k) {
        const double* a = p->left + k * nao;
        const double* b = p->right ? p->right + k * nao : a;
        for (size_t i = 0; i < nao; ++i) {
            const double wa = p->w[k] * a[i];
            if (wa == 0.0) continue;
            for (size_t j = 0; j < nao; ++j) D[i * nao + j] += wa * b[j];
        }
    }
    return D;
}

// max |A - B| / max |A|
static double iso_rel_diff(const double* a, const double* b, size_t n) {
    double m = 0.0, d = 0.0;
    for (size_t i = 0; i < n; ++i) { m = fmax(m, fabs(a[i])); d = fmax(d, fabs(a[i] - b[i])); }
    return d / fmax(m, 1e-300);
}

static md_gto_iso_field_t* iso_dens_field(const iso_dens_src_t* s, const iso_dens_prod_t* p, bool matrix) {
    if (!matrix) return md_gto_iso_field_create_products(&s->basis, s->xyz, 0, p->left, p->right, p->w, p->count, md_get_heap_allocator());
    double* D = iso_dens_matrix(p, s->num_ao);
    md_gto_iso_field_t* f = md_gto_iso_field_create_density(&s->basis, s->xyz, 0, D, md_get_heap_allocator());
    free(D);
    return f;
}

// Random points around the molecule, a third of them near atoms where densities are steep.
static void iso_dens_points(double* pts, int n, const iso_dens_src_t* s, iso_rng_t* rng) {
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    for (size_t a = 0; a < s->num_atoms; ++a) for (int k = 0; k < 3; ++k) { lo[k] = fmin(lo[k], s->xyz[3 * a + k]); hi[k] = fmax(hi[k], s->xyz[3 * a + k]); }
    for (int i = 0; i < n; ++i) {
        if (i % 3 == 0) {
            const size_t a = iso_rng_u64(rng) % s->num_atoms;
            for (int k = 0; k < 3; ++k) pts[3 * i + k] = s->xyz[3 * a + k] + 0.5 * iso_rng_normal(rng);
        } else {
            for (int k = 0; k < 3; ++k) pts[3 * i + k] = iso_rng_range(rng, lo[k] - 4.0, hi[k] + 4.0);
        }
    }
}

static const char* iso_dens_files[] = {"h2o.h5", "amide.h5"};

UTEST(gto_iso_density, matches_products_of_orbitals) {
    for (int fi = 0; fi < 2; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load(&s, iso_dens_files[fi]));
        const size_t nao = s.num_ao;

        // the products agree with the matrices md_vlx publishes for the same densities
        {
            iso_dens_prod_t p;
            iso_dens_products(&p, &s, ISO_D_SCF);
            double* D = iso_dens_matrix(&p, nao);
            size_t dim = 0;
            const double* ref = qm_test_matrix(&s.t, STR_LIT("orbital/total/density"), &dim);
            ASSERT_TRUE(ref != NULL);
            ASSERT_EQ(dim, nao);
            const double e = iso_rel_diff(ref, D, nao * nao);
            printf("  %-9s %zu AOs, %zu occupied, %zu NTO pairs; SCF products vs orbital/total/density: %.1e\n",
                iso_dens_files[fi], nao, p.count, s.num_nto, e);
            EXPECT_LT(e, 1.0e-10);
            free(D);
            iso_dens_prod_free(&p);
        }
        static const char* vlx_path[3] = {"vlx/rsp/transition_density/attachment", "vlx/rsp/transition_density/detachment", "vlx/rsp/transition_density/difference"};
        for (int k = 0; k < 3; ++k) {
            iso_dens_prod_t p;
            iso_dens_products(&p, &s, ISO_D_ATTACH + k);
            double* D = iso_dens_matrix(&p, nao);
            const md_attribute_t* a = qm_test_attr(&s.t, str_from_cstr(vlx_path[k]));
            ASSERT_TRUE(a != NULL);
            double* ref = (double*)malloc(sizeof(double) * nao * nao);
            ASSERT_EQ(md_attribute_extract_f64(ref, nao * nao, a, md_attribute_slice_1(0), md_unit_none()), nao * nao);
            const double e = iso_rel_diff(ref, D, nao * nao);
            printf("            %-6s NTO products vs %s: %.1e\n", iso_d_name[ISO_D_ATTACH + k], vlx_path[k], e);
            // the NTOs kept are those with lambda > 1e-6, so this agrees to about that
            EXPECT_LT(e, 1.0e-4);
            free(ref); free(D);
            iso_dens_prod_free(&p);
        }

        // point values: products field against the orbital evaluator, matrix field against products
        enum { N = 1500 };
        double* pts = (double*)malloc(sizeof(double) * 3 * N);
        iso_rng_t rng = {99};
        iso_dens_points(pts, N, &s, &rng);
        for (int kind = 0; kind < ISO_D_COUNT; ++kind) {
            iso_dens_prod_t p;
            iso_dens_products(&p, &s, kind);
            md_gto_iso_field_t* fp = iso_dens_field(&s, &p, false);
            md_gto_iso_field_t* fm = iso_dens_field(&s, &p, true);
            ASSERT_TRUE(fp != NULL && fm != NULL);
            md_gto_iso_field_info_t ip, im;
            md_gto_iso_field_info(&ip, fp);
            md_gto_iso_field_info(&im, fm);

            // reference: each function as an orbital field (the grouped evaluator, a separate code path)
            const int nfn = (int)(p.right ? 2 * p.count : p.count);
            iso_field_t* orb = (iso_field_t*)calloc(nfn, sizeof(iso_field_t));
            for (size_t k = 0; k < p.count; ++k) {
                ASSERT_TRUE(iso_field_init(&orb[p.right ? 2 * k : k], &s.basis, s.xyz, s.num_atoms, p.left + k * nao, "a"));
                if (p.right) ASSERT_TRUE(iso_field_init(&orb[2 * k + 1], &s.basis, s.xyz, s.num_atoms, p.right + k * nao, "b"));
            }
            double emax_p = 0.0, emax_m = 0.0, egrad = 0.0, vmax = 0.0, gmax = 0.0;
            for (int i = 0; i < N; ++i) {
                const double* x = pts + 3 * i;
                double ref = 0.0;
                for (size_t k = 0; k < p.count; ++k) {
                    const double a = iso_value(&orb[p.right ? 2 * k : k], x);
                    const double b = p.right ? iso_value(&orb[2 * k + 1], x) : a;
                    ref += p.w[k] * a * b;
                }
                double vp, vm, gp[3], gm[3];
                iso_eval(&fp->f, NULL, 0, x, 1, &vp, gp, NULL);
                iso_eval(&fm->f, NULL, 0, x, 1, &vm, gm, NULL);
                vmax = fmax(vmax, fabs(ref));
                emax_p = fmax(emax_p, fabs(vp - ref));
                emax_m = fmax(emax_m, fabs(vm - vp));
                // gradient by central differences
                const double hh = 1.0e-5;
                for (int a = 0; a < 3; ++a) {
                    double xp[3] = {x[0], x[1], x[2]}, xm[3] = {x[0], x[1], x[2]};
                    xp[a] += hh; xm[a] -= hh;
                    const double fd = (iso_value(&fp->f, xp) - iso_value(&fp->f, xm)) / (2.0 * hh);
                    egrad = fmax(egrad, fabs(fd - gp[a]));
                    gmax = fmax(gmax, fabs(gp[a]));
                }
            }
            printf("    %-7s %3u functions (matrix: %3u, %s, residual %.1e): |products - orbitals| %.1e, |matrix - products| %.1e, |grad - fd| %.1e (of max %.1e, %.1e)\n",
                iso_d_name[kind], ip.num_functions, im.num_functions,
                im.split == MD_GTO_ISO_SPLIT_LDL ? "LDL" : im.split == MD_GTO_ISO_SPLIT_EIGEN ? "eigen" : "?", im.residual,
                emax_p, emax_m, egrad, vmax, gmax);
            EXPECT_LT(emax_p, 1.0e-12 * fmax(1.0, vmax) + 1.0e-12);
            EXPECT_LT(emax_m, 1.0e-10 * fmax(1.0, vmax));
            EXPECT_LT(egrad, 1.0e-5 * fmax(1.0, gmax));
            EXPECT_LT(im.residual, 1.0e-8);
            // a semidefinite matrix is split by pivoted Cholesky into exactly its rank
            if (kind == ISO_D_SCF || kind == ISO_D_ATTACH || kind == ISO_D_DETACH) {
                EXPECT_EQ(im.split, (uint32_t)MD_GTO_ISO_SPLIT_LDL);
                EXPECT_EQ(im.num_functions, (uint32_t)p.count);
            }
            for (int k = 0; k < nfn; ++k) iso_field_free(&orb[k]);
            free(orb);
            md_gto_iso_field_destroy(fp);
            md_gto_iso_field_destroy(fm);
            iso_dens_prod_free(&p);
        }
        free(pts);
        iso_dens_src_free(&s);
    }
}

// The bound the octree uses, fuzzed for containment: random boxes of every size, plus the octree's
// own state - groups dropped on a parent box, carried as per-function slack into the child.
UTEST(gto_iso_density, bounds_contain_samples) {
    iso_density_build_groups = true;   // the GROUP form is fuzzed too
    const int boxes_per_size = iso_full() ? 200 : 40;
    const int samples = 64;
    for (int fi = 0; fi < 2; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load(&s, iso_dens_files[fi]));
        for (int kind = 0; kind < ISO_D_COUNT; ++kind) {
            for (int matrix = 0; matrix < 2; ++matrix) {
                iso_dens_prod_t p;
                iso_dens_products(&p, &s, kind);
                md_gto_iso_field_t* field = iso_dens_field(&s, &p, matrix);
                ASSERT_TRUE(field != NULL);
                const iso_field_t* f = &field->f;
                const int nf = f->dens->num_fn, ng = f->num_shells;
                iv_t* A = (iv_t*)malloc(sizeof(iv_t) * (nf + 1));
                iv_t* gint = (iv_t*)malloc(sizeof(iv_t) * (ng + 1));
                double* S = (double*)malloc(sizeof(double) * (nf + 1));
                double* slack = (double*)calloc(nf + 1, sizeof(double));
                int* all = (int*)malloc(sizeof(int) * (ng + 1));
                int* act = (int*)malloc(sizeof(int) * (ng + 1));
                for (int i = 0; i < ng; ++i) all[i] = i;
                iso_rng_t rng = {1234 + kind};
                const int nsh = f->dens->num_ao_shells;
                int* sall = (int*)malloc(sizeof(int) * (nsh + 1));
                int* sact = (int*)malloc(sizeof(int) * (nsh + 1));
                for (int k = 0; k < nsh; ++k) sall[k] = k;
                double* aomag = (double*)malloc(sizeof(double) * (f->dens->num_ao + 1));
                double* acc = (double*)malloc(sizeof(double) * 5 * (nf + 1));
                double* aslack = (double*)calloc(nf + 1, sizeof(double));
                long checks = 0, fails = 0;
                double width_ratio = 0.0, width_ratio_ao = 0.0;
                int nratio = 0;
                for (int si = 0; si < ISO_NUM_SIZES; ++si) {
                    const double size = iso_sizes[si];
                    for (int bi = 0; bi < boxes_per_size; ++bi) {
                        double c[3];
                        iso_dens_points(c, 1, &s, &rng);
                        if (bi & 1) { const size_t a = iso_rng_u64(&rng) % s.num_atoms; for (int k = 0; k < 3; ++k) c[k] = s.xyz[3 * a + k] + iso_rng_range(&rng, -1.5, 1.5); }
                        iso_box_t parent, child;
                        for (int k = 0; k < 3; ++k) { parent.lo[k] = c[k] - size; parent.hi[k] = c[k] + size; }
                        for (int k = 0; k < 3; ++k) { child.lo[k] = parent.lo[k] + (bi & (2 << k) ? size : 0.0); child.hi[k] = child.lo[k] + size; }
                        // the parent's pass: drop groups as the octree would at iso 0.01
                        for (int j = 0; j < nf; ++j) slack[j] = 0.0;
                        iso_density_bound(f, all, ng, slack, &parent, A, gint, S);
                        int na = 0;
                        for (int i = 0; i < ng; ++i) {
                            const double m = iv_mag(gint[i]);
                            if (m * S[f->shells[i].fn] < 1.0e-6 * 0.01) slack[f->shells[i].fn] += m; else act[na++] = i;
                        }
                        double* zero = (double*)calloc(nf + 1, sizeof(double));
                        const iv_t full = iso_density_bound(f, all, ng, zero, &child, A, gint, S);
                        const iv_t carried = iso_density_bound(f, act, na, slack, &child, A, gint, S);
                        // the AO form, with its own drop pass (whole shells) on the parent
                        const iso_density_t* dd = f->dens;
                        for (int j = 0; j < nf; ++j) aslack[j] = 0.0;
                        iso_density_bound_ao(f, sall, nsh, zero, &parent, A, aomag, S, acc);
                        int nas = 0;
                        for (int k = 0; k < nsh; ++k) {
                            const iso_ao_shell_t* a = &dd->ao_shells[k];
                            double effect = 0.0;
                            for (int c = 0; c < a->ncart; ++c) {
                                double e = 0.0;
                                for (int j = 0; j < nf; ++j) e += dd->Cabs[(size_t)(a->ao + c) * nf + j] * S[j];
                                effect += e * aomag[a->ao + c];
                            }
                            if (effect < 1.0e-6 * 0.01) {
                                for (int c = 0; c < a->ncart; ++c) for (int j = 0; j < nf; ++j) aslack[j] += dd->Cabs[(size_t)(a->ao + c) * nf + j] * aomag[a->ao + c];
                            } else {
                                sact[nas++] = k;
                            }
                        }
                        const iv_t ao_full = iso_density_bound_ao(f, sall, nsh, zero, &child, A, aomag, S, acc);
                        const iv_t ao_carried = iso_density_bound_ao(f, sact, nas, aslack, &child, A, aomag, S, acc);
                        free(zero);
                        double vlo = DBL_MAX, vhi = -DBL_MAX;
                        for (int k = 0; k < samples; ++k) {
                            double x[3];
                            for (int a = 0; a < 3; ++a) {
                                const double t = (k < 8) ? ((k >> a) & 1) : (k == 8 ? 0.5 : iso_rng_f64(&rng));
                                x[a] = child.lo[a] + t * (child.hi[a] - child.lo[a]);
                            }
                            const double v = iso_value(f, x);
                            vlo = fmin(vlo, v); vhi = fmax(vhi, v);
                            checks += 4;
                            const int bad = !iv_contains(full, v) + !iv_contains(carried, v) + !iv_contains(ao_full, v) + !iv_contains(ao_carried, v);
                            fails += bad;
                            if (bad && fails < 4) {
                                printf("    FAIL %s %s size %g v %.6e group [%.6e %.6e] carried [%.6e %.6e] ao [%.6e %.6e] carried [%.6e %.6e]\n", iso_d_name[kind], matrix ? "matrix" : "products", size, v,
                                    full.lo, full.hi, carried.lo, carried.hi, ao_full.lo, ao_full.hi, ao_carried.lo, ao_carried.hi);
                            }
                        }
                        if (vhi - vlo > 1e-8) {
                            width_ratio += log(iv_width(full) / (vhi - vlo));
                            width_ratio_ao += log(iv_width(ao_full) / (vhi - vlo));
                            nratio++;
                        }
                    }
                }
                printf("  %-9s %-6s %-8s %3d functions %5d groups: %ld checks, %ld outside; bound width / sampled range, geo. mean: group %.2f, AO %.2f\n",
                    iso_dens_files[fi], iso_d_name[kind], matrix ? "matrix" : "products", nf, ng, checks, fails, exp(width_ratio / fmax(1, nratio)), exp(width_ratio_ao / fmax(1, nratio)));
                EXPECT_EQ(fails, 0);
                free(sall); free(sact); free(aomag); free(acc); free(aslack);
                free(A); free(gint); free(S); free(slack); free(all); free(act);
                md_gto_iso_field_destroy(field);
                iso_dens_prod_free(&p);
            }
        }
        iso_dens_src_free(&s);
    }
    iso_density_build_groups = false;
}


// The cache over a density reproduces dense Surface Nets on its lattice, both sheets, for every iso
// value of its interval - for products and for a split matrix alike.
UTEST(gto_iso_density, cache_serves_its_interval_exactly) {
    typedef struct { int file, kind, matrix; double iso_min, iso_max; } dcase_t;
    static const dcase_t cases[] = {
        {0, ISO_D_SCF,   0, 0.01,  0.1},
        {0, ISO_D_TRANS, 1, 0.002, 0.02},
        {1, ISO_D_SCF,   1, 0.01,  0.1},
        {1, ISO_D_DIFF,  0, 0.002, 0.02},
        {1, ISO_D_TRANS, 0, 0.002, 0.02},
    };
    for (int ci = 0; ci < (int)(sizeof(cases) / sizeof(cases[0])); ++ci) {
        const dcase_t* dc = &cases[ci];
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load(&s, iso_dens_files[dc->file]));
        iso_dens_prod_t p;
        iso_dens_products(&p, &s, dc->kind);
        md_gto_iso_field_t* field = iso_dens_field(&s, &p, dc->matrix);
        ASSERT_TRUE(field != NULL);

        md_gto_iso_cache_desc_t cd = {.cell_size = 0.3, .iso_min = dc->iso_min, .iso_max = dc->iso_max};
        md_gto_iso_cache_t* cache = md_gto_iso_cache_create(field, &cd, md_get_heap_allocator());
        ASSERT_TRUE(cache != NULL);
        md_gto_iso_cache_info_t info;
        md_gto_iso_cache_info(&info, cache);
        printf("  %-9s %-6s %-8s [%.3f, %.3f] h %.2f: %u bricks of %u^3 cells, octree %.3fs, corners %.3fs\n",
            iso_dens_files[dc->file], iso_d_name[dc->kind], dc->matrix ? "matrix" : "products", info.iso_min, info.iso_max, info.cell_size,
            info.num_bricks, info.cells_per_axis, info.seconds_octree, info.seconds_corners);

        const double isos[3] = {dc->iso_min, sqrt(dc->iso_min * dc->iso_max), dc->iso_max};
        for (int k = 0; k < 3; ++k) {
            md_gto_iso_mesh_t raw;
            ASSERT_TRUE(md_gto_iso_mesh_extract(&raw, cache, &(md_gto_iso_mesh_desc_t){.iso = isos[k]}, md_get_heap_allocator()));
            int dnv[2], dnq[2];
            iso_dense_counts(field, &info, isos[k], dnv, dnq);
            printf("     iso %.4f: %u + %u vertices (dense lattice %d + %d), %u + %u triangles (dense %d + %d)\n", isos[k],
                raw.sheet_vertex_count[0], raw.sheet_vertex_count[1], dnv[0], dnv[1],
                raw.sheet_index_count[0] / 3, raw.sheet_index_count[1] / 3, 2 * dnq[0], 2 * dnq[1]);
            for (int sh = 0; sh < 2; ++sh) {
                EXPECT_EQ((int)raw.sheet_vertex_count[sh], dnv[sh]);
                EXPECT_EQ((int)raw.sheet_index_count[sh], 6 * dnq[sh]);
            }
            md_gto_iso_mesh_free(&raw);
        }
        md_gto_iso_cache_destroy(cache);
        md_gto_iso_field_destroy(field);
        iso_dens_prod_free(&p);
        iso_dens_src_free(&s);
    }
}

UTEST(gto_iso_density, display_meshes_are_closed_and_on_the_surface) {
    typedef struct { int file, kind; double iso; } dcase_t;
    static const dcase_t cases[] = {
        {1, ISO_D_SCF,   0.002},
        {1, ISO_D_SCF,   0.05},
        {1, ISO_D_DIFF,  0.002},
        {1, ISO_D_TRANS, 0.002},
    };
    const double h = 0.15;
    for (int ci = 0; ci < (int)(sizeof(cases) / sizeof(cases[0])); ++ci) {
        const dcase_t* dc = &cases[ci];
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load(&s, iso_dens_files[dc->file]));
        iso_dens_prod_t p;
        iso_dens_products(&p, &s, dc->kind);
        md_gto_iso_field_t* field = iso_dens_field(&s, &p, false);
        ASSERT_TRUE(field != NULL);
        md_tick_t t0 = md_tick_now();
        md_gto_iso_cache_t* cache = md_gto_iso_cache_create(field, &(md_gto_iso_cache_desc_t){.cell_size = h, .iso_min = dc->iso}, md_get_heap_allocator());
        const double t_cache = iso_seconds_since(t0);
        ASSERT_TRUE(cache != NULL);
        md_gto_iso_cache_info_t info;
        md_gto_iso_cache_info(&info, cache);

        md_gto_iso_mesh_t raw, q;
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_iso_mesh_extract(&raw, cache, &(md_gto_iso_mesh_desc_t){.iso = dc->iso}, md_get_heap_allocator()));
        const double t_raw = iso_seconds_since(t0);
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_iso_mesh_extract(&q, cache, &(md_gto_iso_mesh_desc_t){.iso = dc->iso, .project_iterations = 8, .relax_iterations = 3}, md_get_heap_allocator()));
        const double t_q = iso_seconds_since(t0);
        printf("  %-9s %-6s iso %.3f h %.2f: cache %.3fs (octree %.3fs, %u bricks, corners %.3fs), raw mesh %.3fs, display %.3fs, %u vertices, %u triangles\n",
            iso_dens_files[dc->file], iso_d_name[dc->kind], dc->iso, h, t_cache, info.seconds_octree, info.num_bricks, info.seconds_corners,
            t_raw, t_q, q.vertex_count, q.index_count / 3);
        EXPECT_EQ(raw.vertex_count, q.vertex_count);
        EXPECT_EQ(raw.index_count, q.index_count);
        EXPECT_GT(q.vertex_count, 0u);
        for (int sh = 0; sh < 2; ++sh) {
            if (!q.sheet_vertex_count[sh]) continue;
            iso_mesh_t mr, mq;
            iso_api_sheet(&mr, &raw, sh, dc->iso, h);
            iso_api_sheet(&mq, &q, sh, dc->iso, h);
            iso_metrics_t er, eq;
            iso_mesh_metrics(&er, &mr, &field->f);
            iso_mesh_metrics(&eq, &mq, &field->f);
            printf("     sheet %c raw      ", sh == 0 ? '+' : '-'); iso_print_metrics("", &er);
            printf("     sheet %c display  ", sh == 0 ? '+' : '-'); iso_print_metrics("", &eq);
            EXPECT_EQ(eq.boundary_edges, 0);
            EXPECT_EQ(eq.euler, er.euler);
            EXPECT_GT(er.orient_ok, 0.99);
            EXPECT_GT(eq.orient_ok, 0.99);
            EXPECT_LT(eq.dev_rms, 0.05);
            iso_mesh_free(&mr);
            iso_mesh_free(&mq);
        }
        md_gto_iso_mesh_free(&raw);
        md_gto_iso_mesh_free(&q);
        md_gto_iso_cache_destroy(cache);
        md_gto_iso_field_destroy(field);
        iso_dens_prod_free(&p);
        iso_dens_src_free(&s);
    }
}


// desc.defer_corners: the caller evaluates brick ranges, here from threads, and gets the same cache.
typedef struct iso_range_job_t { md_gto_iso_cache_t* cache; uint32_t first, count; } iso_range_job_t;
static void iso_range_job(void* p) { iso_range_job_t* j = (iso_range_job_t*)p; md_gto_iso_cache_evaluate(j->cache, j->first, j->count); }

UTEST(gto_iso_density, deferred_corners_from_threads) {
    iso_dens_src_t s;
    ASSERT_TRUE(iso_dens_src_load(&s, "amide.h5"));
    iso_dens_prod_t p;
    iso_dens_products(&p, &s, ISO_D_SCF);
    md_gto_iso_field_t* field = iso_dens_field(&s, &p, true);
    ASSERT_TRUE(field != NULL);
    md_gto_iso_cache_desc_t cd = {.cell_size = 0.2, .iso_min = 0.002, .iso_max = 0.02};
    md_gto_iso_cache_t* ref = md_gto_iso_cache_create(field, &cd, md_get_heap_allocator());
    cd.defer_corners = true;
    md_gto_iso_cache_t* def = md_gto_iso_cache_create(field, &cd, md_get_heap_allocator());
    ASSERT_TRUE(ref && def);
    md_gto_iso_cache_info_t info;
    md_gto_iso_cache_info(&info, def);
    enum { T = 4 };
    iso_range_job_t jobs[T];
    md_thread_t* th[T];
    const uint32_t per = (info.num_bricks + T - 1) / T;
    md_tick_t t0 = md_tick_now();
    for (int k = 0; k < T; ++k) {
        jobs[k] = (iso_range_job_t){def, k * per, per};
        th[k] = md_thread_create(iso_range_job, &jobs[k]);
    }
    for (int k = 0; k < T; ++k) md_thread_join(th[k]);
    printf("  %u bricks over %d threads: %.3fs (one thread: %.3fs)\n", info.num_bricks, T, iso_seconds_since(t0), ref->seconds_corners);
    ASSERT_EQ(ref->store.num_bricks, def->store.num_bricks);
    long diff = 0;
    for (size_t i = 0; i < 125 * (size_t)ref->store.num_bricks; ++i) {
        const float a = ref->store.corner[i], b = def->store.corner[i];
        diff += !((a == b) || (a != a && b != b));
    }
    EXPECT_EQ(diff, 0);
    md_gto_iso_mesh_t ma, mb;
    ASSERT_TRUE(md_gto_iso_mesh_extract(&ma, ref, &(md_gto_iso_mesh_desc_t){.iso = 0.005}, md_get_heap_allocator()));
    ASSERT_TRUE(md_gto_iso_mesh_extract(&mb, def, &(md_gto_iso_mesh_desc_t){.iso = 0.005}, md_get_heap_allocator()));
    EXPECT_EQ(ma.vertex_count, mb.vertex_count);
    EXPECT_EQ(ma.index_count, mb.index_count);
    md_gto_iso_mesh_free(&ma);
    md_gto_iso_mesh_free(&mb);
    md_gto_iso_cache_destroy(ref);
    md_gto_iso_cache_destroy(def);
    md_gto_iso_field_destroy(field);
    iso_dens_prod_free(&p);
    iso_dens_src_free(&s);
}

// An orbital as a one-function linear "density" (GTO_ISO_REPORT=1): does the AO-form bound, which a
// unified GPU pipeline would use for both kinds, keep as few bricks as the orbital GROUP bound?
UTEST(gto_iso_density, orbital_as_linear_term) {
    if (!iso_report()) { printf("  (GTO_ISO_REPORT=1 to run)\n"); return; }
    const double h = 0.15;
    for (int ci = 0; ci < ISO_NUM_CASES; ++ci) {
        const iso_case_t* tc = &iso_cases[ci];
        iso_api_case_t c;
        ASSERT_TRUE(iso_api_load(&c, tc));
        const iso_term_t term = {1.0, 0, -1};
        md_gto_iso_field_t* lin = iso_density_create(&c.basis, c.xyz, 0, c.coeff, 1, &term, 1, 0.0, ISO_FACTOR_NONE, md_get_heap_allocator(), "linear");
        ASSERT_TRUE(lin != NULL);
        iso_octree_t og, oa;
        iso_octree_build_range(&og, &c.field->f, tc->iso, tc->iso, 4.0 * h, ISO_BOUND_SHELL, 0, 0.0);
        iso_octree_build_range(&oa, &lin->f, tc->iso, tc->iso, 4.0 * h, ISO_BOUND_SHELL, 0, 0.0);
        long vg = 0, va = 0;
        for (int l = 0; l <= og.L; ++l) vg += og.visited[l];
        for (int l = 0; l <= oa.L; ++l) va += oa.visited[l];
        // the AO form must keep every brick the dense lattice needs: compare with the GROUP set
        iso_map_t gm;
        iso_map_init(&gm, (uint32_t)md_array_size(oa.leaves) + 1);
        for (size_t i = 0; i < md_array_size(oa.leaves); ++i) iso_map_put(&gm, iso_key3(oa.leaves[i].x, oa.leaves[i].y, oa.leaves[i].z), 1);
        int missing = 0;
        for (size_t i = 0; i < md_array_size(og.leaves); ++i) { uint32_t v; missing += !iso_map_get(&gm, iso_key3(og.leaves[i].x, og.leaves[i].y, og.leaves[i].z), &v); }
        printf("  %-13s GROUP %5d bricks (%6ld nodes, %.3fs) | AO form %5d bricks (%6ld nodes, %.3fs) | %+5.1f%% bricks, %d GROUP bricks not in AO set\n",
            tc->name, (int)md_array_size(og.leaves), vg, og.seconds, (int)md_array_size(oa.leaves), va, oa.seconds,
            100.0 * ((double)md_array_size(oa.leaves) / md_array_size(og.leaves) - 1.0), missing);
        iso_map_free(&gm);
        iso_octree_free(&og);
        iso_octree_free(&oa);
        md_gto_iso_field_destroy(lin);
        iso_api_free(&c);
    }
}

// How the way a density is split changes the cost (GTO_ISO_REPORT=1): the same SCF density as its
// canonical MOs (the eigenvectors of D in the S metric), as pivoted Cholesky factors, and as the
// Euclidean eigenvectors of D; difference and transition densities as NTO products and as a matrix.
static md_gto_iso_field_t* iso_dens_field_eigen(const iso_dens_src_t* s, const iso_dens_prod_t* p) {
    const int n = (int)s->num_ao;
    double* D = iso_dens_matrix(p, s->num_ao);
    double* V = (double*)malloc(sizeof(double) * n * n);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) V[i * n + j] = 0.5 * (D[i * n + j] + D[j * n + i]);
    double* ev = (double*)malloc(sizeof(double) * n);
    iso_eigen_sym(V, ev, n);
    double emax = 0.0;
    for (int k = 0; k < n; ++k) emax = fmax(emax, fabs(ev[k]));
    double* rows = (double*)malloc(sizeof(double) * n * n);
    iso_term_t* terms = (iso_term_t*)malloc(sizeof(iso_term_t) * n);
    int r = 0;
    for (int k = 0; k < n; ++k) {
        if (fabs(ev[k]) <= 1.0e-12 * emax) continue;
        for (int i = 0; i < n; ++i) rows[r * n + i] = V[i * n + k];
        terms[r] = (iso_term_t){ev[k], r, r};
        r++;
    }
    md_gto_iso_field_t* f = iso_density_create(&s->basis, s->xyz, 0, rows, r, terms, r, 0.0, ISO_FACTOR_EIGEN, md_get_heap_allocator(), "eigen");
    free(terms); free(rows); free(ev); free(V); free(D);
    return f;
}

UTEST(gto_iso_density, split_comparison) {
    if (!iso_report()) { printf("  (GTO_ISO_REPORT=1 to run)\n"); return; }
    iso_density_build_groups = true;
    static const char* files[] = {"amide.h5", "mol.h5"};
    const double h = 0.15;
    for (int fi = 0; fi < 2; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load(&s, files[fi]));
        printf("  %s: %zu atoms, %zu AOs, %zu NTO pairs\n", files[fi], s.num_atoms, s.num_ao, s.num_nto);
        typedef struct { int kind; int split; double iso; } row_t;   // split: 0 products, 1 matrix, 2 eigen
        static const row_t rows[] = {
            {ISO_D_SCF, 0, 0.002}, {ISO_D_SCF, 1, 0.002}, {ISO_D_SCF, 2, 0.002},
            {ISO_D_SCF, 0, 0.05},  {ISO_D_SCF, 1, 0.05},
            {ISO_D_DIFF, 0, 0.002}, {ISO_D_DIFF, 1, 0.002},
            {ISO_D_TRANS, 0, 0.002}, {ISO_D_TRANS, 1, 0.002},
        };
        static const char* split_name[3] = {"products", "LDL", "eigen"};
        printf("    %-6s %-8s %6s | %4s %6s %6s | %-5s %7s %8s %7s %9s | %7s %7s\n", "", "split", "iso", "fn", "groups", "setup", "bound", "octree", "visited", "bricks", "act/node", "corners", "per cnr");
        for (int ri = 0; ri < (int)(sizeof(rows) / sizeof(rows[0])); ++ri) {
            const row_t* r = &rows[ri];
            iso_dens_prod_t p;
            iso_dens_products(&p, &s, r->kind);
            md_tick_t t0 = md_tick_now();
            md_gto_iso_field_t* field = r->split == 2 ? iso_dens_field_eigen(&s, &p) : iso_dens_field(&s, &p, r->split == 1);
            const double t_setup = iso_seconds_since(t0);
            ASSERT_TRUE(field != NULL);
            md_gto_iso_field_info_t fin;
            md_gto_iso_field_info(&fin, field);

            for (int mode = 0; mode < 2; ++mode) {
            iso_dens_bound_default = mode == 0 ? ISO_DENS_BOUND_AO : ISO_DENS_BOUND_GROUP;
            iso_octree_t o;
            iso_octree_build_range(&o, &field->f, r->iso, r->iso, 4.0 * h, ISO_BOUND_SHELL, 0, 0.0);
            long visited = 0;
            for (int l = 0; l <= o.L; ++l) visited += o.visited[l];
            const int nb = (int)md_array_size(o.leaves);

            md_gto_iso_cache_t* cache = md_gto_iso_cache_create(field, &(md_gto_iso_cache_desc_t){.cell_size = h, .iso_min = r->iso}, md_get_heap_allocator());
            md_gto_iso_cache_info_t ci;
            md_gto_iso_cache_info(&ci, cache);
            const double corners = (double)nb * 125.0;
            printf("    %-6s %-8s %6.3f | %4u %6u %5.2fs | %-5s %6.3fs %8ld %7d %9.1f | %6.3fs %5.2fus\n",
                iso_d_name[r->kind], split_name[r->split], r->iso, fin.num_functions, fin.num_groups, t_setup, mode == 0 ? "AO" : "group",
                o.seconds, visited, nb, (double)o.group_work / fmax(1, visited), ci.seconds_corners, 1e6 * ci.seconds_corners / fmax(1, corners));
            md_gto_iso_cache_destroy(cache);
            iso_octree_free(&o);
            }
            iso_dens_bound_default = ISO_DENS_BOUND_AO;
            md_gto_iso_field_destroy(field);
            iso_dens_prod_free(&p);
        }
        iso_dens_src_free(&s);
    }
    iso_density_build_groups = false;
}

// =============================================================================================
// GPU: the same pipeline on md_gpu (shaders/gto_iso.slang)
// =============================================================================================
//
// The production shape the CPU experiments pointed at: GROUP bound only, octree stopped at 4^3
// cell bricks, bricks filled densely. One launch per octree level, then per brick / cell /
// vertex / quad. Every launch after the first is sized by a device-side count (iso_grid writes
// the indirect grid); the host reads back once, the brick count after the octree, to size the
// mesh buffers.

// =============================================================================================
// Average local ionization energy (md_qm_alie_orbitals): the weights, and I(r) from two product fields
// =============================================================================================
#include <md_qm.h>

typedef struct iso_alie_t {
    md_gto_iso_field_t* rho;    // sum n |phi|^2
    md_gto_iso_field_t* num;    // sum n |e| |phi|^2
    double homo;                // |e| of the least bound occupied orbital, hartree
    size_t count;
} iso_alie_t;

static bool iso_alie_make(iso_alie_t* a, const iso_dens_src_t* s) {
    memset(a, 0, sizeof(*a));
    double *rows = NULL, *w = NULL, *we = NULL;
    size_t na = 0;
    a->count = md_qm_alie_orbitals(&rows, &w, &we, &na, &s->t.sys, md_get_heap_allocator());
    if (!a->count || na != s->num_ao) return false;
    a->homo = DBL_MAX;
    for (size_t k = 0; k < a->count; ++k) a->homo = fmin(a->homo, we[k] / w[k]);
    a->rho = md_gto_iso_field_create_products(&s->basis, s->xyz, 0, rows, NULL, w, a->count, md_get_heap_allocator());
    a->num = md_gto_iso_field_create_products(&s->basis, s->xyz, 0, rows, NULL, we, a->count, md_get_heap_allocator());
    md_free(md_get_heap_allocator(), rows, sizeof(double) * a->count * na);
    md_free(md_get_heap_allocator(), w, sizeof(double) * a->count);
    md_free(md_get_heap_allocator(), we, sizeof(double) * a->count);
    return a->rho && a->num;
}

static void iso_alie_free(iso_alie_t* a) {
    if (a->rho) md_gto_iso_field_destroy(a->rho);
    if (a->num) md_gto_iso_field_destroy(a->num);
    memset(a, 0, sizeof(*a));
}

UTEST(qm_alie, orbitals_and_limits) {
    static const char* files[] = {"h2o.h5", "amide.h5"};
    for (int fi = 0; fi < 2; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load_ex(&s, files[fi], false));
        ASSERT_TRUE(md_qm_alie_available(&s.t.sys));
        double *rows = NULL, *w = NULL, *we = NULL;
        size_t na = 0;
        const size_t count = md_qm_alie_orbitals(&rows, &w, &we, &na, &s.t.sys, md_get_heap_allocator());
        ASSERT_GT(count, 0u);
        // the same orbitals and weights as the total density's
        double *drows = NULL, *dw = NULL;
        size_t dna = 0;
        const size_t dcount = md_qm_density_orbitals(&drows, &dw, &dna, &s.t.sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
        ASSERT_EQ(count, dcount);
        EXPECT_EQ(0, memcmp(rows, drows, sizeof(double) * count * na));
        EXPECT_EQ(0, memcmp(w, dw, sizeof(double) * count));
        // closed shell: n = 2 and the energy weights 2 |e| of the occupied orbitals in order
        double* e = (double*)calloc(s.num_mo, sizeof(double));
        ASSERT_EQ(qm_test_series(e, s.num_mo, &s.t, STR_LIT("orbital/alpha/energy")), s.num_mo);
        size_t k = 0;
        for (size_t m = 0; m < s.num_mo && k < count; ++m) {
            if (s.occ[m] == 0.0) continue;
            EXPECT_NEAR(w[k], 2.0, 1e-12);
            EXPECT_NEAR(we[k], 2.0 * fabs(e[m]), 1e-12);
            ++k;
        }
        EXPECT_EQ(k, count);
        free(e);
        md_free(md_get_heap_allocator(), rows, sizeof(double) * count * na);
        md_free(md_get_heap_allocator(), w, sizeof(double) * count);
        md_free(md_get_heap_allocator(), we, sizeof(double) * count);
        md_free(md_get_heap_allocator(), drows, sizeof(double) * dcount * dna);
        md_free(md_get_heap_allocator(), dw, sizeof(double) * dcount);

        // I(r): between |e_HOMO| and the deepest |e| everywhere, and on the 0.001 isodensity surface in the
        // usual 0.3 - 0.8 hartree. Not |e_HOMO| far out, as it would be with exact orbitals: in a gaussian basis
        // the tails are the most diffuse primitives' (and the evaluation cutoff's), so far from the molecule
        // I(r) means nothing (water 8.7 bohr out, rho 7e-12: 0.49 against |e_HOMO| 0.29)
        iso_alie_t a;
        ASSERT_TRUE(iso_alie_make(&a, &s));
        enum { N = 2000 };
        iso_rng_t rng = {11};
        double* pd = (double*)malloc(sizeof(double) * 3 * N);
        float* pts = (float*)malloc(sizeof(float) * 3 * N);
        iso_dens_points(pd, N, &s, &rng);
        for (int i = 0; i < 3 * N; ++i) pts[i] = (float)pd[i];
        double* r = (double*)malloc(sizeof(double) * N);
        double* q = (double*)malloc(sizeof(double) * N);
        md_gto_iso_field_evaluate(a.rho, pts, N, 0, r, NULL);
        md_gto_iso_field_evaluate(a.num, pts, N, 0, q, NULL);
        double imin = DBL_MAX, imax = 0.0;
        for (int i = 0; i < N; ++i) {
            if (!(r[i] > 1e-12)) continue;
            const double I = q[i] / r[i];
            imin = fmin(imin, I);
            imax = fmax(imax, I);
        }
        EXPECT_GE(imin, a.homo * (1.0 - 1e-6));
        // on the 0.001 au isodensity surface
        md_gto_iso_cache_t* cache = md_gto_iso_cache_create(a.rho, &(md_gto_iso_cache_desc_t){.cell_size = 0.3, .iso_min = 0.001}, md_get_heap_allocator());
        ASSERT_TRUE(cache != NULL);
        md_gto_iso_mesh_t m;
        ASSERT_TRUE(md_gto_iso_mesh_extract(&m, cache, &(md_gto_iso_mesh_desc_t){.iso = 0.001, .project_iterations = 8, .relax_iterations = 2}, md_get_heap_allocator()));
        const size_t nv = m.sheet_vertex_count[0];
        double* rv = (double*)malloc(sizeof(double) * nv);
        double* qv = (double*)malloc(sizeof(double) * nv);
        const float* sp = m.position + 3 * (size_t)m.sheet_vertex_offset[0];
        md_gto_iso_field_evaluate(a.rho, sp, nv, 0, rv, NULL);
        md_gto_iso_field_evaluate(a.num, sp, nv, 0, qv, NULL);
        double smin = DBL_MAX, smax = 0.0, rdev = 0.0;
        for (size_t i = 0; i < nv; ++i) {
            const double I = qv[i] / rv[i];
            smin = fmin(smin, I);
            smax = fmax(smax, I);
            rdev = fmax(rdev, fabs(rv[i] - 0.001) / 0.001);
        }
        printf("  %-9s %zu occupied, |e_HOMO| %.4f Eh: I(r) at %d points %.4f .. %.4f; on rho = 0.001 (%zu vertices, rho within %.0e) %.4f .. %.4f Eh (%.2f .. %.2f eV)\n",
               files[fi], a.count, a.homo, N, imin, imax, nv, rdev, smin, smax, smin * 27.211386, smax * 27.211386);
        EXPECT_GE(smin, a.homo * (1.0 - 1e-6));
        EXPECT_GT(smin, 0.2);
        EXPECT_LT(smax, 1.0);
        free(rv); free(qv);
        md_gto_iso_mesh_free(&m);
        md_gto_iso_cache_destroy(cache);
        free(r); free(q); free(pd); free(pts);
        iso_alie_free(&a);
        iso_dens_src_free(&s);
    }
}

// Local electron attachment energy (md_qm_attachment_orbitals): the bound virtuals, each closed shell spatial
// orbital twice (Brinck et al. 2016, eq. 10 and Table S1), E(r) <= 0, nothing when no virtual is bound. Checked
// against the paper's Table S4 (CH3Cl -0.28, CH2FCl -0.54, CHF2Cl -1.14, CF3Cl -2.45 eV at the chlorine sigma
// hole, 0.004 au) with B3LYP/6-31+G(d,p) Molden files from PySCF: -0.293, -0.542, -1.142, -2.451 eV.
UTEST(qm_leae, orbitals_and_sign) {
    static const char* files[] = {"h2o.h5", "amide.h5", "mol.h5"};
    for (int fi = 0; fi < 3; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load_ex(&s, files[fi], false));
        double* e = (double*)calloc(s.num_mo, sizeof(double));
        ASSERT_EQ(qm_test_series(e, s.num_mo, &s.t, STR_LIT("orbital/alpha/energy")), s.num_mo);
        size_t bound = 0;
        for (size_t m = 0; m < s.num_mo; ++m) bound += s.occ[m] == 0.0 && e[m] < 0.0;
        EXPECT_EQ(md_qm_attachment_available(&s.t.sys), bound > 0);
        double *rows = NULL, *w = NULL;
        size_t na = 0;
        const size_t count = md_qm_attachment_orbitals(&rows, &w, &na, &s.t.sys, md_get_heap_allocator());
        EXPECT_EQ(count, bound);
        printf("  %-9s %zu bound virtual orbitals", files[fi], bound);
        if (count) {
            // closed shell: each bound virtual once, weight 2 e, in order, its own coefficients
            size_t k = 0;
            for (size_t m = 0; m < s.num_mo && k < count; ++m) {
                if (!(s.occ[m] == 0.0 && e[m] < 0.0)) continue;
                EXPECT_NEAR(w[k], 2.0 * e[m], 1e-12);
                EXPECT_EQ(0, memcmp(rows + k * na, s.C + m * s.num_ao, sizeof(double) * na));
                printf(" %.3f eV", e[m] * 27.211386);
                ++k;
            }
            // E(r) <= 0 at random points
            double *dr = NULL, *dw = NULL;
            size_t dna = 0;
            const size_t dn = md_qm_density_orbitals(&dr, &dw, &dna, &s.t.sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
            md_gto_iso_field_t* rho = md_gto_iso_field_create_products(&s.basis, s.xyz, 0, dr, NULL, dw, dn, md_get_heap_allocator());
            md_gto_iso_field_t* num = md_gto_iso_field_create_products(&s.basis, s.xyz, 0, rows, NULL, w, count, md_get_heap_allocator());
            enum { N = 1000 };
            iso_rng_t rng = {5};
            double* pd = (double*)malloc(sizeof(double) * 3 * N);
            float* pts = (float*)malloc(sizeof(float) * 3 * N);
            iso_dens_points(pd, N, &s, &rng);
            for (int i = 0; i < 3 * N; ++i) pts[i] = (float)pd[i];
            double r[N], q[N];
            md_gto_iso_field_evaluate(rho, pts, N, 0, r, NULL);
            md_gto_iso_field_evaluate(num, pts, N, 0, q, NULL);
            int positive = 0;
            for (int i = 0; i < N; ++i) positive += q[i] > 0.0;
            EXPECT_EQ(positive, 0);
            md_gto_iso_field_destroy(rho);
            md_gto_iso_field_destroy(num);
            md_free(md_get_heap_allocator(), dr, sizeof(double) * dn * dna);
            md_free(md_get_heap_allocator(), dw, sizeof(double) * dn);
            md_free(md_get_heap_allocator(), rows, sizeof(double) * count * na);
            md_free(md_get_heap_allocator(), w, sizeof(double) * count);
            free(pd); free(pts);
        }
        printf("\n");
        free(e);
        iso_dens_src_free(&s);
    }
}

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#include "gto_iso_shaders.inl"

typedef struct iso_gpu_group_t {
    md_gpu_float4 c;
    md_gpu_float4 cut;
    md_gpu_uint4  info;
    md_gpu_uint4  info2;
} iso_gpu_group_t;

typedef struct iso_gpu_args_t {        // mirrors IsoArgs
    md_gpu_float4 origin_h;
    uint32_t      num_groups, level, brick_level, cap;
    float         iso, rel_eps, step_max, lambda;
    uint32_t      hash_mask, newton_iters, grid_mul, _pad0;
    md_gpu_addr_t groups, alpha, weights, residues, crit;
    md_gpu_addr_t nodes_in, nodes_out, count_in, count_out, grid_out;
    md_gpu_addr_t hash_keys, hash_vals, corners, cell_vert;
    md_gpu_addr_t vpos, vpos_out, vnrm, vcell, quads, tris, counters;
} iso_gpu_args_t;

enum { ISO_K_GRID, ISO_K_LEVEL, ISO_K_HASH, ISO_K_CORNERS, ISO_K_CELLS, ISO_K_QUADS, ISO_K_PROJECT, ISO_K_RELAX, ISO_K_TRIS, ISO_K_COUNT };

typedef struct iso_gpu_t {
    md_gpu_device_t dev;
    md_gpu_stream_t s;
    md_gpu_kernel_t k[ISO_K_COUNT];
} iso_gpu_t;

static bool iso_gpu_init(iso_gpu_t* g) {
    memset(g, 0, sizeof(*g));
    // GTO_ISO_VALIDATE=1: Vulkan validation layers through md_gpu (errors are logged)
    const char* v = getenv("GTO_ISO_VALIDATE");
    md_gpu_device_desc_t dd = {0};
    dd.enable_validation = v && v[0] == '1';
    dd.label = "gto_iso";
    g->dev = md_gpu_device_create(&dd);
    if (!g->dev) return false;
    g->s = md_gpu_stream_default(g->dev, MD_GPU_STREAM_COMPUTE);
    md_gto_gpu_initialize(g->dev);     // the D form's matrices in tiles are md_gto's (md_gto_gpu_dmat)
    const md_gpu_kernel_desc_t d[ISO_K_COUNT] = {
        md_shader_gto_iso_iso_grid_kernel(),
        md_shader_gto_iso_iso_octree_level_kernel(),
        md_shader_gto_iso_iso_brick_hash_kernel(),
        md_shader_gto_iso_iso_corners_kernel(),
        md_shader_gto_iso_iso_cells_kernel(),
        md_shader_gto_iso_iso_quads_kernel(),
        md_shader_gto_iso_iso_project_kernel(),
        md_shader_gto_iso_iso_relax_kernel(),
        md_shader_gto_iso_iso_triangulate_kernel(),
    };
    for (int i = 0; i < ISO_K_COUNT; ++i) {
        g->k[i] = md_gpu_kernel_create(g->dev, &d[i]);
        if (!g->k[i]) return false;
    }
    return true;
}

static void iso_gpu_free(iso_gpu_t* g) {
    if (g->dev) md_gto_gpu_shutdown();
    if (g->dev) md_gpu_device_destroy(g->dev);
    memset(g, 0, sizeof(*g));
}

// All device buffers of one extraction are temp allocations in the extraction's scope.
static md_gpu_addr_t iso_gpu_alloc(iso_gpu_t* g, size_t bytes) {
    return md_gpu_temp_alloc(g->s, MD_GPU_MEM_DEVICE, bytes ? bytes : 16).gpu;
}

static md_gpu_addr_t iso_gpu_upload(iso_gpu_t* g, const void* src, size_t bytes) {
    md_gpu_addr_t a = iso_gpu_alloc(g, bytes);
    if (a && bytes) md_gpu_upload(g->s, a, src, bytes);
    return a;
}

static md_gpu_addr_t iso_gpu_zeros(iso_gpu_t* g, size_t bytes) {
    md_gpu_addr_t a = iso_gpu_alloc(g, bytes);
    if (a) md_gpu_memset(g->s, a, 0, bytes ? bytes : 16);
    return a;
}

// Blocking readback (prototype convenience): device -> host-read memory -> dst.
static void iso_gpu_read(iso_gpu_t* g, void* dst, md_gpu_addr_t src, size_t bytes) {
    if (!bytes) return;
    md_gpu_mem_t m = md_gpu_malloc(g->s, MD_GPU_MEM_HOST_READ, bytes);
    md_gpu_copy(g->s, m.gpu, src, bytes);
    md_gpu_stream_sync(g->s);
    memcpy(dst, m.cpu, bytes);
    md_gpu_free(g->s, m.gpu);
}

static inline float iso_f32_down(double x) { float f = (float)x; return (double)f > x ? nextafterf(f, -INFINITY) : f; }
static inline float iso_f32_up(double x)   { float f = (float)x; return (double)f < x ? nextafterf(f,  INFINITY) : f; }

// Radius beyond which |psi_g| < eps (tail bound, valid past every radial maximum).
static double iso_shell_radius(const iso_shell_t* s, double eps) {
    double r = iso_shell_tail_r0(s);
    while (iso_shell_tail(s, r) > eps) r += 0.02;
    return r;
}

typedef struct iso_gpu_field_t {
    md_gpu_addr_t groups, alpha, weights, residues, crit;
    uint32_t      num_groups;
} iso_gpu_field_t;

static void iso_gpu_field_upload(iso_gpu_t* g, iso_gpu_field_t* out, const iso_field_t* f, double iso) {
    iso_gpu_group_t* groups = (iso_gpu_group_t*)calloc(f->num_shells + 1, sizeof(iso_gpu_group_t));
    float* alpha = NULL; float* w = NULL; md_gpu_uint4* residues = NULL; md_gpu_float4* crit = NULL;
    md_allocator_i* heap = md_get_heap_allocator();
    // Total slack from groups dropped by distance stays below 1e-3 iso.
    const double eps_drop = 1.0e-3 * iso / fmax(1, f->num_shells);
    for (int gi = 0; gi < f->num_shells; ++gi) {
        const iso_shell_t* s = &f->shells[gi];
        const double r_eval = iso_shell_radius(s, 1.0e-9);
        const double r_drop = iso_shell_radius(s, eps_drop);
        iso_gpu_group_t* G = &groups[gi];
        G->c   = (md_gpu_float4){(float)s->c[0], (float)s->c[1], (float)s->c[2], (float)(r_eval * r_eval)};
        // the drop test uses the box's nearest rho, so round the radius down and the epsilon up
        G->cut = (md_gpu_float4){iso_f32_down(r_drop * r_drop * (1.0 - 1e-6)), iso_f32_up(eps_drop), 0, 0};
        G->info  = (md_gpu_uint4){(uint32_t)s->l, (uint32_t)s->ncart, (uint32_t)s->nprim, (uint32_t)md_array_size(alpha)};
        G->info2 = (md_gpu_uint4){(uint32_t)md_array_size(w), (uint32_t)md_array_size(residues), 0, 0};
        for (int q = 0; q < s->nprim; ++q) { float a = (float)s->alpha[q]; md_array_push(alpha, a, heap); }
        for (int n = 0; n < s->ncart; ++n) {
            for (int q = 0; q < s->nprim; ++q) { float x = (float)s->w[n][q]; md_array_push(w, x, heap); }
            md_gpu_uint4 cp = {(uint32_t)(s->ijk[n][0] | (s->ijk[n][1] << 4) | (s->ijk[n][2] << 8)), (uint32_t)md_array_size(crit), (uint32_t)s->crit_num[n][0], s->live[n] ? 1u : 0u};
            for (int k = 0; k < s->crit_num[n][0]; ++k) {
                const iso_crit_t* c = &f->crit[s->crit_off[n][0] + k];
                md_gpu_float4 cf = {iso_f32_down(c->lo), iso_f32_up(c->hi), iso_f32_down(c->g.lo - 1e-7 * fabs(c->g.lo)), iso_f32_up(c->g.hi + 1e-7 * fabs(c->g.hi))};
                md_array_push(crit, cf, heap);
            }
            md_array_push(residues, cp, heap);
        }
    }
    out->num_groups = (uint32_t)f->num_shells;
    out->groups  = iso_gpu_upload(g, groups, sizeof(iso_gpu_group_t) * f->num_shells);
    out->alpha   = iso_gpu_upload(g, alpha, sizeof(float) * md_array_size(alpha));
    out->weights = iso_gpu_upload(g, w, sizeof(float) * md_array_size(w));
    out->residues = iso_gpu_upload(g, residues, sizeof(md_gpu_uint4) * md_array_size(residues));
    out->crit    = iso_gpu_upload(g, crit, sizeof(md_gpu_float4) * md_array_size(crit));
    free(groups);
    md_array_free(alpha, heap); md_array_free(w, heap); md_array_free(residues, heap); md_array_free(crit, heap);
}

typedef struct iso_gpu_mesh_t {
    uint32_t  nb, nv, nq, open_edges, overflow;
    uint32_t* bricks;      // 4 per brick: x, y, z, sheets
    float*    corners;     // 125 per brick
    float*    raw_pos;     // 4 per vertex, before projection
    float*    pos;         // 4 per vertex
    float*    nrm;         // 4 per vertex
    uint32_t* quads;       // 4 per quad
    uint32_t* tris;        // 6 per quad
    double    origin[3];
    double    h;
    int       L_brick;     // octree levels to the bricks
    double    t_octree, t_mesh;
} iso_gpu_mesh_t;

static void iso_gpu_mesh_free(iso_gpu_mesh_t* m) {
    free(m->bricks); free(m->corners); free(m->raw_pos); free(m->pos); free(m->nrm); free(m->quads); free(m->tris);
    memset(m, 0, sizeof(*m));
}

#define ISO_GPU_NODE_CAP (1u << 20)

enum { ISO_GPU_STAGE_SETUP, ISO_GPU_STAGE_OCTREE, ISO_GPU_STAGE_CORNERS, ISO_GPU_STAGE_NETS, ISO_GPU_STAGE_PROJECT, ISO_GPU_STAGE_RELAX, ISO_GPU_STAGE_TRIS, ISO_GPU_STAGE_COUNT };
static const char* iso_gpu_stage_name[ISO_GPU_STAGE_COUNT] = {"setup", "octree", "corners", "s.nets", "project", "relax", "tris"};

typedef struct iso_gpu_opts_t {
    int     relax_iters;
    bool    no_readback;       // benchmark: read back counters only
    double* stage_seconds;     // [ISO_GPU_STAGE_COUNT]: sync after every stage and time it (adds syncs)
} iso_gpu_opts_t;

static void iso_gpu_stage(iso_gpu_t* g, const iso_gpu_opts_t* o, int stage, md_tick_t* t) {
    if (!o->stage_seconds) return;
    md_gpu_stream_sync(g->s);
    const md_tick_t now = md_tick_now();
    o->stage_seconds[stage] = md_tick_to_seconds(now - *t);
    *t = now;
}

static bool iso_gpu_extract_ex(iso_gpu_t* g, const iso_field_t* f, double iso, double h, const iso_gpu_opts_t* opt, iso_gpu_mesh_t* out) {
    memset(out, 0, sizeof(*out));
    md_gpu_temp_t scope = md_gpu_temp_begin(g->s);
    (void)scope;
    md_tick_t t0 = md_tick_now();
    md_tick_t ts = t0;
    const int relax_iters = opt->relax_iters;

    iso_gpu_field_t fld;
    iso_gpu_field_upload(g, &fld, f, iso);

    const int Lb = iso_root(f, iso, 4.0 * h, out->origin);
    out->h = h;
    out->L_brick = Lb;
    if (Lb >= 10) { md_gpu_temp_end(g->s, scope); return false; }   // brick coordinates are packed in 10 bits

    iso_gpu_args_t a = {0};
    a.origin_h   = (md_gpu_float4){(float)out->origin[0], (float)out->origin[1], (float)out->origin[2], (float)h};
    a.num_groups = fld.num_groups;
    a.brick_level = (uint32_t)Lb;
    a.iso        = (float)iso;
    a.rel_eps    = 1.0e-4f;
    a.step_max   = (float)(0.5 * h);
    a.lambda     = 0.5f;
    a.groups = fld.groups; a.alpha = fld.alpha; a.weights = fld.weights; a.residues = fld.residues; a.crit = fld.crit;
    iso_gpu_stage(g, opt, ISO_GPU_STAGE_SETUP, &ts);

    // ---- octree, one launch per level, sized by the previous level's device-side count
    md_gpu_addr_t nodes[2] = {iso_gpu_alloc(g, sizeof(md_gpu_uint4) * ISO_GPU_NODE_CAP), iso_gpu_alloc(g, sizeof(md_gpu_uint4) * ISO_GPU_NODE_CAP)};
    md_gpu_addr_t counts   = iso_gpu_zeros(g, sizeof(uint32_t) * 16);
    md_gpu_addr_t grid     = iso_gpu_zeros(g, sizeof(uint32_t) * 4);
    md_gpu_addr_t counters = iso_gpu_zeros(g, sizeof(uint32_t) * 4);
    a.counters = counters;
    a.grid_out = grid;
    const md_gpu_uint4 root = {0, 0, 0, 3};
    const uint32_t one = 1;
    md_gpu_upload(g->s, nodes[0], &root, sizeof(root));
    md_gpu_upload(g->s, counts, &one, sizeof(one));
    for (int l = 0; l < Lb; ++l) {
        a.level     = (uint32_t)l;
        a.nodes_in  = nodes[l & 1];
        a.nodes_out = nodes[(l + 1) & 1];
        a.count_in  = counts + sizeof(uint32_t) * l;
        a.count_out = counts + sizeof(uint32_t) * (l + 1);
        a.cap       = ISO_GPU_NODE_CAP;
        a.grid_mul  = 1;
        MD_GPU_LAUNCH(g->s, g->k[ISO_K_GRID], md_gpu_grid(1, 1, 1), a);
        md_gpu_launch_indirect(g->s, g->k[ISO_K_LEVEL], grid, &a, sizeof(a));
    }
    const md_gpu_addr_t bricks = nodes[Lb & 1];
    const md_gpu_addr_t brick_count = counts + sizeof(uint32_t) * Lb;

    // the one readback: how many bricks, to size the mesh buffers
    uint32_t nb = 0, ctr[4] = {0};
    iso_gpu_read(g, &nb, brick_count, sizeof(nb));
    iso_gpu_read(g, ctr, counters, sizeof(ctr));
    out->t_octree = iso_seconds_since(t0);
    iso_gpu_stage(g, opt, ISO_GPU_STAGE_OCTREE, &ts);
    if (ctr[3] || nb > ISO_GPU_NODE_CAP) { out->overflow = ctr[3] + 1; md_gpu_temp_end(g->s, scope); return false; }
    out->nb = nb;
    if (nb == 0) { md_gpu_temp_end(g->s, scope); return true; }

    // ---- bricks: hash, corners, Surface Nets
    md_tick_t t1 = md_tick_now();
    uint32_t hcap = 64;
    while (hcap < 2 * nb) hcap <<= 1;
    const uint32_t cap_v = nb * 128u, cap_q = nb * 64u * 2u * 3u;
    a.nodes_in  = bricks;
    a.count_in  = brick_count;
    a.hash_mask = hcap - 1;
    a.hash_keys = iso_gpu_zeros(g, sizeof(uint32_t) * hcap);
    a.hash_vals = iso_gpu_alloc(g, sizeof(uint32_t) * hcap);
    a.corners   = iso_gpu_alloc(g, sizeof(float) * 125 * nb);
    a.cell_vert = iso_gpu_alloc(g, sizeof(uint32_t) * 128 * nb);
    a.vpos      = iso_gpu_alloc(g, sizeof(md_gpu_float4) * cap_v);
    a.vpos_out  = iso_gpu_alloc(g, sizeof(md_gpu_float4) * cap_v);
    a.vnrm      = iso_gpu_alloc(g, sizeof(md_gpu_float4) * cap_v);
    a.vcell     = iso_gpu_alloc(g, sizeof(md_gpu_uint4) * cap_v);
    a.quads     = iso_gpu_alloc(g, sizeof(md_gpu_uint4) * cap_q);
    a.tris      = iso_gpu_alloc(g, sizeof(uint32_t) * 6 * cap_q);

    a.cap = hcap;
    MD_GPU_LAUNCH(g->s, g->k[ISO_K_HASH], md_gpu_grid_for(g->k[ISO_K_HASH], nb, 1, 1), a);
    MD_GPU_LAUNCH(g->s, g->k[ISO_K_CORNERS], md_gpu_grid_for(g->k[ISO_K_CORNERS], nb * 128, 1, 1), a);
    iso_gpu_stage(g, opt, ISO_GPU_STAGE_CORNERS, &ts);
    a.cap = cap_v;
    MD_GPU_LAUNCH(g->s, g->k[ISO_K_CELLS], md_gpu_grid_for(g->k[ISO_K_CELLS], nb * 64, 1, 1), a);
    a.cap = cap_q;
    MD_GPU_LAUNCH(g->s, g->k[ISO_K_QUADS], md_gpu_grid_for(g->k[ISO_K_QUADS], nb * 64, 1, 1), a);

    iso_gpu_stage(g, opt, ISO_GPU_STAGE_NETS, &ts);
    // raw vertex positions, kept for the comparison with the CPU lattice
    md_gpu_addr_t raw = 0;
    if (!opt->no_readback) {
        raw = iso_gpu_alloc(g, sizeof(md_gpu_float4) * cap_v);
        md_gpu_copy(g->s, raw, a.vpos, sizeof(md_gpu_float4) * cap_v);
    }

    // ---- per vertex: projection, relaxation; per quad: triangles. Sized by device counts.
    iso_gpu_args_t v = a;
    v.count_in = counters;                 // vertices
    v.cap = cap_v;
    v.grid_mul = 1;
    MD_GPU_LAUNCH(g->s, g->k[ISO_K_GRID], md_gpu_grid(1, 1, 1), v);
    v.newton_iters = 8;
    md_gpu_launch_indirect(g->s, g->k[ISO_K_PROJECT], grid, &v, sizeof(v));
    iso_gpu_stage(g, opt, ISO_GPU_STAGE_PROJECT, &ts);
    v.newton_iters = 3;
    for (int it = 0; it < relax_iters; ++it) {
        md_gpu_launch_indirect(g->s, g->k[ISO_K_RELAX], grid, &v, sizeof(v));
        md_gpu_copy(g->s, v.vpos, v.vpos_out, sizeof(md_gpu_float4) * cap_v);
        md_gpu_launch_indirect(g->s, g->k[ISO_K_PROJECT], grid, &v, sizeof(v));
    }
    iso_gpu_stage(g, opt, ISO_GPU_STAGE_RELAX, &ts);
    iso_gpu_args_t q = a;
    q.count_in = counters + sizeof(uint32_t);   // quads
    q.cap = cap_q;
    q.grid_mul = 1;
    MD_GPU_LAUNCH(g->s, g->k[ISO_K_GRID], md_gpu_grid(1, 1, 1), q);
    md_gpu_launch_indirect(g->s, g->k[ISO_K_TRIS], grid, &q, sizeof(q));
    iso_gpu_stage(g, opt, ISO_GPU_STAGE_TRIS, &ts);

    iso_gpu_read(g, ctr, counters, sizeof(ctr));
    out->t_mesh = iso_seconds_since(t1);
    out->nv = ctr[0] < cap_v ? ctr[0] : cap_v;
    out->nq = ctr[1] < cap_q ? ctr[1] : cap_q;
    out->open_edges = ctr[2];
    out->overflow = ctr[3];
    if (opt->no_readback) { md_gpu_temp_end(g->s, scope); return true; }

    out->bricks  = (uint32_t*)malloc(sizeof(md_gpu_uint4) * nb);
    out->corners = (float*)malloc(sizeof(float) * 125 * nb);
    out->raw_pos = (float*)malloc(sizeof(float) * 4 * (out->nv + 1));
    out->pos     = (float*)malloc(sizeof(float) * 4 * (out->nv + 1));
    out->nrm     = (float*)malloc(sizeof(float) * 4 * (out->nv + 1));
    out->quads   = (uint32_t*)malloc(sizeof(uint32_t) * 4 * (out->nq + 1));
    out->tris    = (uint32_t*)malloc(sizeof(uint32_t) * 6 * (out->nq + 1));
    iso_gpu_read(g, out->bricks, bricks, sizeof(md_gpu_uint4) * nb);
    iso_gpu_read(g, out->corners, a.corners, sizeof(float) * 125 * nb);
    iso_gpu_read(g, out->raw_pos, raw, sizeof(float) * 4 * out->nv);
    iso_gpu_read(g, out->pos, a.vpos, sizeof(float) * 4 * out->nv);
    iso_gpu_read(g, out->nrm, a.vnrm, sizeof(float) * 4 * out->nv);
    iso_gpu_read(g, out->quads, a.quads, sizeof(uint32_t) * 4 * out->nq);
    iso_gpu_read(g, out->tris, a.tris, sizeof(uint32_t) * 6 * out->nq);
    md_gpu_temp_end(g->s, scope);
    return true;
}

static bool iso_gpu_extract(iso_gpu_t* g, const iso_field_t* f, double iso, double h, int relax_iters, iso_gpu_mesh_t* out) {
    const iso_gpu_opts_t opt = {.relax_iters = relax_iters};
    return iso_gpu_extract_ex(g, f, iso, h, &opt, out);
}

// One sheet of the GPU mesh as an iso_mesh_t (double), so the CPU metrics apply unchanged.
static void iso_gpu_mesh_sheet(iso_mesh_t* m, const iso_gpu_mesh_t* gm, int sheet, double iso, bool raw) {
    md_allocator_i* heap = md_get_heap_allocator();
    memset(m, 0, sizeof(*m));
    m->sheet = sheet; m->iso = iso; m->h = gm->h;
    uint32_t* remap = (uint32_t*)malloc(sizeof(uint32_t) * (gm->nv + 1));
    const float* P = raw ? gm->raw_pos : gm->pos;
    for (uint32_t i = 0; i < gm->nv; ++i) {
        remap[i] = UINT32_MAX;
        if ((P[4 * i + 3] > 0) != (sheet > 0)) continue;
        remap[i] = (uint32_t)iso_mesh_nv(m);
        const double p[3] = {P[4 * i], P[4 * i + 1], P[4 * i + 2]};
        md_array_push_array(m->pos, p, 3, heap);
        if (!raw) {
            const double n[3] = {gm->nrm[4 * i], gm->nrm[4 * i + 1], gm->nrm[4 * i + 2]};
            md_array_push_array(m->nrm, n, 3, heap);
            const double r = gm->nrm[4 * i + 3];
            md_array_push(m->res, r, heap);
        }
    }
    for (uint32_t qi = 0; qi < gm->nq; ++qi) {
        const uint32_t* q = gm->quads + 4 * qi;
        if (remap[q[0]] == UINT32_MAX) continue;
        const uint32_t qq[4] = {remap[q[0]], remap[q[1]], remap[q[2]], remap[q[3]]};
        md_array_push_array(m->quad, qq, 4, heap);
        const uint32_t* t = gm->tris + 6 * qi;
        const uint32_t tt[6] = {remap[t[0]], remap[t[1]], remap[t[2]], remap[t[3]], remap[t[4]], remap[t[5]]};
        md_array_push_array(m->tri, tt, 6, heap);
    }
    free(remap);
}

static const char* iso_gpu_no_device(void) {
    const char* e = md_gpu_last_error();
    return (e && e[0]) ? e : "No GPU device available";
}

// The GPU must find every cell the double precision CPU path finds, produce the same Surface
// Nets mesh over the same lattice, and project it onto the surface just as well.
UTEST(gto_iso_gpu, matches_cpu) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gpu_device_info_t info;
    md_gpu_device_info(g.dev, &info);
    printf("  device: %s (%s)\n", info.name, info.type == MD_GPU_DEVICE_TYPE_CPU ? "software" : "hardware");

    static const int pick[4] = {0, 2, 3, 5};
    const double h = 0.15;
    for (int ci = 0; ci < 4; ++ci) {
        const iso_case_t* tc = &iso_cases[pick[ci]];
        iso_field_t f;
        ASSERT_TRUE(iso_load(&f, tc->file, tc->which, tc->name));

        iso_gpu_mesh_t gm;
        const bool ok = iso_gpu_extract(&g, &f, tc->iso, h, 3, &gm);
        EXPECT_TRUE(ok);
        EXPECT_EQ(gm.overflow, 0u);
        EXPECT_EQ(gm.open_edges, 0u);
        if (!ok) { iso_field_free(&f); continue; }

        // CPU reference on the identical lattice: same root, GROUP bound, 4^3 bricks, double precision
        md_tick_t t0 = md_tick_now();
        iso_octree_t o;
        iso_octree_build_ex(&o, &f, tc->iso, 4.0 * h, ISO_BOUND_SHELL, 0, 0.0);
        int nc = 0;
        iso_leaf_t* cells = iso_bricks_to_cells(&o, 2, &nc);
        iso_lattice_t lat = {.f = &f, .h = h};
        memcpy(lat.origin, o.origin, sizeof(lat.origin));
        iso_map_init(&lat.corner_map, (uint32_t)nc * 2);
        iso_mesh_t cm[2];
        for (int s = 0; s < 2; ++s) iso_surface_nets(&cm[s], &lat, cells, nc, s == 0 ? 1 : -1, tc->iso);
        const double t_cpu_sn = iso_seconds_since(t0);
        EXPECT_EQ(o.L, gm.L_brick);
        for (int a = 0; a < 3; ++a) EXPECT_NEAR(o.origin[a], gm.origin[a], 1e-12);

        // 1. every brick holding a CPU crossing must be a GPU brick (the float bound stayed conservative)
        iso_map_t gb;
        iso_map_init(&gb, gm.nb);
        for (uint32_t b = 0; b < gm.nb; ++b) iso_map_put(&gb, iso_key3(gm.bricks[4 * b], gm.bricks[4 * b + 1], gm.bricks[4 * b + 2]), b);
        int missing = 0;
        double max_corner_err = 0.0, max_corner_abs = 0.0;
        for (int li = 0; li < (int)md_array_size(o.leaves); ++li) {
            const iso_leaf_t* l = &o.leaves[li];
            uint32_t b;
            const bool have = iso_map_get(&gb, iso_key3(l->x, l->y, l->z), &b);
            // does the CPU see a crossing in this brick?
            bool crossing = false;
            for (int c = 0; c < 64 && !crossing; ++c) {
                const int cx = l->x * 4 + (c & 3), cy = l->y * 4 + ((c >> 2) & 3), cz = l->z * 4 + (c >> 4);
                for (int s = 0; s < 2 && !crossing; ++s) {
                    int mask = 0;
                    for (int k = 0; k < 8; ++k) if ((s == 0 ? 1 : -1) * iso_lattice_corner(&lat, cx + (k & 1), cy + ((k >> 1) & 1), cz + ((k >> 2) & 1)) - tc->iso > 0) mask |= 1 << k;
                    crossing = mask != 0 && mask != 0xFF;
                }
            }
            if (crossing && !have) missing++;
            if (have) {   // 2. corner values: float GPU vs double CPU
                for (int k = 0; k < 125; ++k) {
                    const double cv = iso_lattice_corner(&lat, l->x * 4 + k % 5, l->y * 4 + (k / 5) % 5, l->z * 4 + k / 25);
                    max_corner_err = fmax(max_corner_err, fabs(cv - gm.corners[125 * b + k]));
                    max_corner_abs = fmax(max_corner_abs, fabs(cv));
                }
            }
        }
        EXPECT_EQ(missing, 0);
        EXPECT_LT(max_corner_err, 1e-5 * max_corner_abs);   // float vs double

        printf("\n  %-13s bricks gpu %u / cpu %d, corner |err| max %.2e (max |psi| %.2f)\n", tc->name, gm.nb, (int)md_array_size(o.leaves), max_corner_err, max_corner_abs);
        printf("     gpu (%s): octree %.3fs, mesh + project + 3x relax %.3fs | cpu Surface Nets only %.3fs\n",
            info.type == MD_GPU_DEVICE_TYPE_CPU ? "llvmpipe" : "device", gm.t_octree, gm.t_mesh, t_cpu_sn);

        FILE* obj = NULL;
        if (iso_out_dir()) {
            char name[128];
            snprintf(name, sizeof(name), "gto_iso_gpu_%s.obj", tc->name);
            obj = iso_open_out(name, "w");
        }
        int vbase = 0;
        // 3. the same Surface Nets mesh, and after projection + relaxation a mesh as good as the CPU's
        for (int s = 0; s < 2; ++s) {
            const int sheet = s == 0 ? 1 : -1;
            iso_mesh_t gr, gp;
            iso_gpu_mesh_sheet(&gr, &gm, sheet, tc->iso, true);
            iso_gpu_mesh_sheet(&gp, &gm, sheet, tc->iso, false);
            double cs_g[2], cs_c[2];
            iso_mesh_checksum(&gr, cs_g);
            iso_mesh_checksum(&cm[s], cs_c);
            const int nv_c = iso_mesh_nv(&cm[s]), nq_c = (int)(md_array_size(cm[s].quad) / 4);
            printf("     sheet %c: vertices gpu %d / cpu %d, quads gpu %d / cpu %d, raw position checksum rel diff %.1e\n",
                sheet > 0 ? '+' : '-', iso_mesh_nv(&gr), nv_c, (int)(md_array_size(gr.quad) / 4), nq_c, fabs(cs_g[0] - cs_c[0]) / fmax(1.0, fabs(cs_c[0])));
            EXPECT_EQ(iso_mesh_nv(&gr), nv_c);
            EXPECT_EQ((int)(md_array_size(gr.quad) / 4), nq_c);
            EXPECT_LT(fabs(cs_g[0] - cs_c[0]) / fmax(1.0, fabs(cs_c[0])), 1e-5);

            iso_metrics_t mc, mg;
            iso_mesh_triangulate(&cm[s]);
            iso_mesh_project(&cm[s], &f, 8);
            iso_mesh_triangulate(&cm[s]);
            iso_mesh_relax(&cm[s], &f, 3, 0.5);
            iso_mesh_metrics(&mc, &cm[s], &f);
            iso_mesh_metrics(&mg, &gp, &f);
            iso_print_metrics("cpu", &mc);
            iso_print_metrics("gpu", &mg);
            EXPECT_EQ(mg.euler, mc.euler);
            EXPECT_EQ(mg.components, mc.components);
            EXPECT_EQ(mg.boundary_edges, 0);
            EXPECT_GT(mg.orient_ok, 0.99);
            EXPECT_LT(mg.dev_rms, 2.0 * mc.dev_rms + 0.005);
            if (obj) {
                iso_mesh_write_obj(obj, &gp, vbase, sheet > 0 ? "positive" : "negative");
                vbase += iso_mesh_nv(&gp);
            }
            iso_mesh_free(&gr);
            iso_mesh_free(&gp);
            iso_mesh_free(&cm[s]);
        }
        if (obj) fclose(obj);
        iso_map_free(&gb);
        free(cells);
        iso_map_free(&lat.corner_map);
        md_array_free(lat.corner_val, md_get_heap_allocator());
        iso_octree_free(&o);
        iso_gpu_mesh_free(&gm);
        iso_field_free(&f);
    }
    iso_gpu_free(&g);
}
// ---------------------------------------------------------------------------------------------
// Benchmark (GTO_ISO_BENCH=1): gridless extraction on md_gpu against mdlib's dense GPU orbital
// evaluation of the same lattice, which is what an isosurface costs today before any surface
// extraction even starts.
//
//   GTO_ISO_BENCH=1 ./md_unittest --filter='gto_iso_gpu.benchmark'
//   GTO_ISO_BENCH_H=0.2,0.1,0.05    cell sizes (bohr), default 0.15,0.1,0.075,0.05
//   GTO_ISO_BENCH_REPS=7            timed repetitions (median reported), default 5
//
// "total" is one extraction with only the two syncs the pipeline needs (brick count, final
// counters), host side field packing included. "dense" evaluates psi at spacing h over the tight
// padded box (not the power of two octree root) with md_gto_gpu_orbital_launch (default algorithm) and no surface extraction at all; dense/tot
// is how many times longer that takes. The per stage columns come from a separate run that syncs after every stage, so they
// include launch and sync overhead and sum to more than the total.
// ---------------------------------------------------------------------------------------------

static int iso_parse_list(const char* s, double* out, int cap) {
    int n = 0;
    while (s && *s && n < cap) {
        char* end;
        const double v = strtod(s, &end);
        if (end == s) break;
        out[n++] = v;
        s = (*end == ',') ? end + 1 : end;
    }
    return n;
}

static double iso_median(double* v, int n) { return iso_percentile(v, n, 0.5); }

typedef struct iso_bench_case_t { const char* file; int which; const char* name; } iso_bench_case_t;

UTEST(gto_iso_gpu, benchmark) {
    const char* on = getenv("GTO_ISO_BENCH");
    if (!on || on[0] != '1') UTEST_SKIP("set GTO_ISO_BENCH=1");
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gpu_device_info_t info;
    md_gpu_device_info(g.dev, &info);
    md_gto_gpu_initialize(g.dev);

    double hs[16];
    int nh = iso_parse_list(getenv("GTO_ISO_BENCH_H"), hs, 16);
    if (nh == 0) { hs[0] = 0.15; hs[1] = 0.1; hs[2] = 0.075; hs[3] = 0.05; nh = 4; }
    const int reps = getenv("GTO_ISO_BENCH_REPS") ? atoi(getenv("GTO_ISO_BENCH_REPS")) : 5;
    const double iso = 0.05;

    static const iso_bench_case_t cases[] = {
        {"amide.h5",          ISO_MO_HOMO, "amide_homo"},
        {"mol.h5",            ISO_MO_HOMO, "mol_homo"},
        {"mol.h5",            ISO_MO_LUMO, "mol_lumo"},
        {"acro-xps.h5",       ISO_MO_HOMO, "acro_homo"},
        {"tq.scf.results.h5", ISO_MO_HOMO, "tq_homo"},
    };

    printf("\n  device: %s, %s | iso %.3f | median of %d\n", info.name, info.driver, iso, reps);
    printf("  %-11s %6s %6s %6s %7s %6s | %8s", "orbital", "h", "root", "bricks", "corners", "verts", "total ms");
    for (int st = 0; st < ISO_GPU_STAGE_COUNT; ++st) printf(" %7s", iso_gpu_stage_name[st]);
    printf(" | %9s %8s %14s\n", "dense ms", "dense/tot", "dense points");

    for (int ci = 0; ci < (int)(sizeof(cases) / sizeof(cases[0])); ++ci) {
        const iso_bench_case_t* bc = &cases[ci];
        vlx_test_t t = {0};
        char path[512];
        snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", bc->file);
        if (!vlx_test_load(&t, str_from_cstr(path), MEGABYTES(256))) { printf("  %-11s (no %s, skipped)\n", bc->name, bc->file); continue; }
        const size_t num_atoms = t.sys.atom.count;
        float* xyz = (float*)malloc(sizeof(float) * 3 * num_atoms);
        md_gto_basis_t basis = {0};
        if (vlx_test_atom_xyz_bohr(xyz, 3 * num_atoms, &t) != num_atoms || !qm_test_basis(&basis, &t)) {
            printf("  %-11s (could not read, skipped)\n", bc->name);
            free(xyz); qm_test_free(&t); continue;
        }
        const size_t num_ao = md_gto_basis_num_ao(&basis);
        double* c = (double*)calloc(num_ao, sizeof(double));
        const size_t lumo = vlx_test_lumo_idx(&t, STR_LIT("orbital/alpha/occupation"));
        qm_test_row(c, num_ao, &t, STR_LIT("orbital/alpha/coefficient"), bc->which == ISO_MO_LUMO ? lumo : lumo - 1);
        iso_field_t f;
        if (!iso_field_init(&f, &basis, xyz, num_atoms, c, bc->name)) { free(xyz); free(c); qm_test_free(&t); continue; }
        printf("  %-11s %zu atoms, %u shells, %zu AOs, %d (atom,l) groups\n", bc->name, num_atoms, basis.num_shells, num_ao, f.num_shells);

        // dense baseline: mdlib's GPU orbital evaluation into a 3D texture
        md_gto_gpu_basis_t gb = md_gto_gpu_basis_create(g.s, &(md_gto_gpu_basis_desc_t){ .basis = &basis, .cutoff = 1.0e-6 });
        md_gpu_mem_t atoms = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(num_atoms));
        md_gpu_mem_t coeff = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, md_gto_gpu_coeff_size_mo(1, num_ao));
        {
            float* p = (float*)md_gpu_upload_begin(g.s, atoms.gpu, md_gto_gpu_atom_buffer_size(num_atoms));
            md_gto_gpu_atom_pack(p, xyz, 0, num_atoms);
            md_gpu_upload_end(g.s);
            const double* rows[1] = {c};
            p = (float*)md_gpu_upload_begin(g.s, coeff.gpu, md_gto_gpu_coeff_size_mo(1, num_ao));
            md_gto_gpu_coeff_pack_mo(p, rows, NULL, 1, num_ao);
            md_gpu_upload_end(g.s);
        }

        for (int hi = 0; hi < nh; ++hi) {
            const double h = hs[hi];
            double origin[3];
            const int Lb = iso_root(&f, iso, 4.0 * h, origin);
            const int N = 4 << Lb;   // cells per axis

            // gridless: warm up, then untimed totals and per stage timings
            iso_gpu_mesh_t gm;
            iso_gpu_opts_t o = {.relax_iters = 3, .no_readback = true};
            if (!iso_gpu_extract_ex(&g, &f, iso, h, &o, &gm)) { printf("  %-11s %6.3f (failed: overflow or lattice too large)\n", bc->name, h); continue; }
            double tot[64], stage[ISO_GPU_STAGE_COUNT][64];
            for (int r = 0; r < reps && r < 64; ++r) {
                md_tick_t t0 = md_tick_now();
                iso_gpu_extract_ex(&g, &f, iso, h, &o, &gm);
                tot[r] = iso_seconds_since(t0) * 1e3;
            }
            double st[ISO_GPU_STAGE_COUNT];
            iso_gpu_opts_t os = {.relax_iters = 3, .no_readback = true, .stage_seconds = st};
            for (int r = 0; r < reps && r < 64; ++r) {
                iso_gpu_extract_ex(&g, &f, iso, h, &os, &gm);
                for (int k = 0; k < ISO_GPU_STAGE_COUNT; ++k) stage[k][r] = st[k] * 1e3;
            }
            const int nr = reps < 64 ? reps : 64;

            // dense: the tight padded box (not the octree's power of two cube) at spacing h,
            // skipped beyond 2^27 points (512 MB of float)
            double dense_ms = NAN;
            double blo[3], bhi[3];
            iso_surface_box(&f, iso, blo, bhi);
            int dd[3];
            for (int a = 0; a < 3; ++a) dd[a] = (int)ceil((bhi[a] - blo[a]) / h) + 1;
            const double nvox = (double)dd[0] * dd[1] * dd[2];
            if (nvox <= (double)(1u << 27)) {
                md_grid_t grid = {
                    .orientation = mat3_ident(),
                    .origin  = vec3_set((float)blo[0], (float)blo[1], (float)blo[2]),
                    .spacing = vec3_set((float)h, (float)h, (float)h),
                    .dim     = {dd[0], dd[1], dd[2]},
                };
                md_gpu_texture_t tex = md_gpu_texture_create(g.s, &(md_gpu_texture_desc_t){
                    .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
                    .width = (uint32_t)dd[0], .height = (uint32_t)dd[1], .depth_or_layers = (uint32_t)dd[2],
                });
                if (tex) {
                    md_gto_gpu_orbital_desc_t od = {
                        .basis = gb, .atom_xyz = atoms.gpu, .coeff = coeff.gpu, .out_tex = tex, .grid = &grid,
                        .num_orbitals = 1, .eval_mode = MD_GTO_EVAL_MODE_PSI, .op = MD_GTO_OP_SET,
                    };
                    md_gto_gpu_orbital_launch(g.s, &od);
                    md_gpu_stream_sync(g.s);
                    double dt[64];
                    for (int r = 0; r < nr; ++r) {
                        md_tick_t t0 = md_tick_now();
                        md_gto_gpu_orbital_launch(g.s, &od);
                        md_gpu_stream_sync(g.s);
                        dt[r] = iso_seconds_since(t0) * 1e3;
                    }
                    dense_ms = iso_median(dt, nr);
                    md_gpu_texture_destroy(tex);
                }
            }

            const long corners = (long)gm.nb * 125;
            printf("  %-11s %6.3f %4d^3 %6u %7ld %6u | %8.2f", bc->name, h, N, gm.nb, corners, gm.nv, iso_median(tot, nr));
            for (int k = 0; k < ISO_GPU_STAGE_COUNT; ++k) printf(" %7.2f", iso_median(stage[k], nr));
            if (isnan(dense_ms)) printf(" | %9s %8s %14s\n", "-", "-", "-");
            else printf(" | %9.2f %7.1fx %4dx%4dx%4d\n", dense_ms, dense_ms / iso_median(tot, nr), dd[0], dd[1], dd[2]);
        }
        md_gpu_free(g.s, atoms.gpu);
        md_gpu_free(g.s, coeff.gpu);
        md_gto_gpu_basis_destroy(g.s, gb);
        iso_field_free(&f);
        free(xyz); free(c);
        qm_test_free(&t);
    }
    md_gto_gpu_shutdown();
    iso_gpu_free(&g);
}

// =============================================================================================
// GPU library path (md_gto_iso_gpu_*, src/shaders/iso/gto_iso_gpu.slang) against the CPU reference
// =============================================================================================

typedef struct iso_gl_case_t { int kind; const char* name; double iso; } iso_gl_case_t;   // kind: -1 HOMO orbital, else ISO_D_*

static md_gto_iso_field_t* iso_gl_field(const iso_dens_src_t* s, int kind, iso_dens_prod_t* p) {
    memset(p, 0, sizeof(*p));
    if (kind < 0) {
        const size_t lumo = vlx_test_lumo_idx((vlx_test_t*)&s->t, STR_LIT("orbital/alpha/occupation"));
        const size_t mo = kind == -1 ? lumo - 1 : lumo;
        return md_gto_iso_field_create(&s->basis, s->xyz, 0, s->C + mo * s->num_ao, md_get_heap_allocator());
    }
    iso_dens_products(p, s, kind);
    return iso_dens_field(s, p, false);
}

// Converts sheet s of an md_gto_iso_mesh_t for the metrics (as iso_api_sheet).
static double iso_gl_seconds(md_tick_t t0) { return iso_seconds_since(t0); }

UTEST(gto_iso_gpu_lib, matches_cpu_reference) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    iso_dens_src_t s;
    // TEST: GTO_ISO_REF_FILE=mol.h5 compares another system (the cases keep their amide labels)
    ASSERT_TRUE(iso_dens_src_load(&s, getenv("GTO_ISO_REF_FILE") ? getenv("GTO_ISO_REF_FILE") : "amide.h5"));
    static const iso_gl_case_t cases[] = {
        {-1,          "amide_homo",  0.05},
        {ISO_D_SCF,   "amide_scf",   0.02},
        {ISO_D_DIFF,  "amide_diff",  0.004},
        {ISO_D_TRANS, "amide_trans", 0.004},
    };
    const double h = 0.25;
    for (int ci = 0; ci < (int)(sizeof(cases) / sizeof(cases[0])); ++ci) {
        const iso_gl_case_t* gc = &cases[ci];
        iso_dens_prod_t p;
        md_gto_iso_field_t* field = iso_gl_field(&s, gc->kind, &p);
        ASSERT_TRUE(field != NULL);

        // CPU reference: the cache's bricks and double precision corners (stored as float)
        md_gto_iso_cache_t* cache = md_gto_iso_cache_create(field, &(md_gto_iso_cache_desc_t){.cell_size = h, .iso_min = gc->iso}, md_get_heap_allocator());
        ASSERT_TRUE(cache != NULL);
        md_gto_iso_mesh_t cm;
        ASSERT_TRUE(md_gto_iso_mesh_extract(&cm, cache, &(md_gto_iso_mesh_desc_t){.iso = gc->iso, .project_iterations = 8, .relax_iterations = 3}, md_get_heap_allocator()));

        // GPU
        const md_gto_iso_field_t* fl[1] = {field};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
        ASSERT_TRUE(set != NULL);
        md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 8, .relax_iterations = 3};
        gd.iso[0] = gc->iso;
        isog_debug_t dbg = {0};
        isog_debug = &dbg;
        md_gto_iso_gpu_mesh_t gm;
        md_tick_t t0 = md_tick_now();
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &gd, &gm));
        const double t_gpu = iso_gl_seconds(t0);
        isog_debug = NULL;

        // same lattice
        EXPECT_NEAR(dbg.origin[0], cache->origin[0], 1e-9);
        EXPECT_EQ(dbg.brick_level, cache->brick_levels);

        // rigor: every cell whose double corners cross lies in a GPU brick with that sheet
        iso_map_t gmap;
        iso_map_init(&gmap, dbg.num_bricks + 1);
        uint32_t gpu_bricks = 0;
        for (uint32_t b = 0; b < dbg.num_bricks; ++b) {
            if (dbg.bricks[4 * b + 3] == 0) continue;   // a rejected sibling that kept its slot (the point list's order)
            iso_map_put(&gmap, iso_key3((int)dbg.bricks[4 * b], (int)dbg.bricks[4 * b + 1], (int)dbg.bricks[4 * b + 2]), b);
            gpu_bricks++;
        }
        long crossing = 0, missed = 0, compared = 0, sign_flips = 0;
        double max_err = 0.0, max_val = 0.0;
        const iso_brick_store_t* st = &cache->store;
        iso_lattice_t lat = {.f = &field->f, .h = h, .store = st};
        memcpy(lat.origin, cache->origin, sizeof(lat.origin));
        for (uint32_t b = 0; b < st->num_bricks; ++b) {
            const int32_t* bc = st->coord + 4 * b;
            uint32_t gb;
            const bool have = iso_map_get(&gmap, iso_key3(bc[0], bc[1], bc[2]), &gb);
            for (int cell = 0; cell < 64; ++cell) {
                const int x = bc[0] * 4 + (cell & 3), y = bc[1] * 4 + ((cell >> 2) & 3), z = bc[2] * 4 + (cell >> 4);
                for (int sh = 0; sh < 2; ++sh) {
                    const double sg = sh == 0 ? 1.0 : -1.0;
                    int pos = 0;
                    for (int k = 0; k < 8; ++k) pos += sg * iso_lattice_corner(&lat, x + (k & 1), y + ((k >> 1) & 1), z + ((k >> 2) & 1)) - gc->iso > 0.0;
                    if (pos == 0 || pos == 8) continue;
                    crossing++;
                    if (!have || !((dbg.bricks[4 * gb + 3] >> sh) & 1u)) missed++;
                }
            }
            if (!have) continue;
            // corner values where both have them
            for (int k = 0; k < 125; ++k) {
                const float cv = st->corner[125 * (size_t)b + k];
                if (cv != cv) continue;   // not owned by this CPU brick
                const float gv = dbg.corners[125 * (size_t)dbg.slots[gb] + k];
                compared++;
                max_err = fmax(max_err, fabs((double)gv - cv));
                max_val = fmax(max_val, fabs(cv));
                for (int sh = 0; sh < 2; ++sh) {
                    const double sg = sh == 0 ? 1.0 : -1.0;
                    sign_flips += (sg * cv - gc->iso > 0.0) != (sg * gv - gc->iso > 0.0);
                }
            }
        }
        iso_map_free(&gmap);

        // the GPU mesh against the double precision field
        md_gto_iso_mesh_t hm;
        ASSERT_TRUE(md_gto_iso_gpu_mesh_read(g.s, &gm, 0, &hm, md_get_heap_allocator()));
        printf("  %-12s h %.2f iso %.3f: bricks GPU %u CPU %u | crossing cells %ld, missed by GPU %ld | corners %ld: max |gpu-cpu| %.1e (max |f| %.1e), sign flips at iso %ld | %.2fs\n",
            gc->name, h, gc->iso, gpu_bricks, st->num_bricks, crossing, missed, compared, max_err, max_val, sign_flips, t_gpu);
        EXPECT_EQ(missed, 0);
        EXPECT_LT(max_err, 1.0e-5 * fmax(max_val, 1.0));
        EXPECT_EQ(gm.open_edges, 0u);
        for (int sh = 0; sh < 2; ++sh) {
            printf("     sheet %c: vertices GPU %u CPU %u, triangles GPU %u CPU %u\n", sh == 0 ? '+' : '-',
                hm.sheet_vertex_count[sh], cm.sheet_vertex_count[sh], hm.sheet_index_count[sh] / 3, cm.sheet_index_count[sh] / 3);
            // counts can differ only through corners whose float and double values straddle the iso value
            const double dv = fabs((double)hm.sheet_vertex_count[sh] - cm.sheet_vertex_count[sh]);
            EXPECT_LE(dv, 4.0 * sign_flips + 0.0);
            if (!hm.sheet_vertex_count[sh]) continue;
            iso_mesh_t mg, mc;
            iso_api_sheet(&mg, &hm, sh, gc->iso, h);
            iso_api_sheet(&mc, &cm, sh, gc->iso, h);
            iso_metrics_t eg, ec;
            iso_mesh_metrics(&eg, &mg, &field->f);
            iso_mesh_metrics(&ec, &mc, &field->f);
            printf("     sheet %c GPU  ", sh == 0 ? '+' : '-'); iso_print_metrics("", &eg);
            printf("     sheet %c CPU  ", sh == 0 ? '+' : '-'); iso_print_metrics("", &ec);
            EXPECT_EQ(eg.boundary_edges, 0);
            EXPECT_EQ(eg.euler, ec.euler);
            EXPECT_GT(eg.orient_ok, 0.99);
            // the same pipeline in float: as close to the surface as the double one
            EXPECT_LT(eg.dev_rms, 1.2 * ec.dev_rms + 0.005);
            EXPECT_LT(eg.vert_res_max, 1.0e-3);
            iso_mesh_free(&mg);
            iso_mesh_free(&mc);
        }
        md_gto_iso_mesh_free(&hm);
        md_gto_iso_gpu_mesh_free(g.s, &gm);
        md_gto_iso_gpu_set_destroy(g.s, set);
        free(dbg.bricks); free(dbg.slots); free(dbg.corners);
        md_gto_iso_mesh_free(&cm);
        md_gto_iso_cache_destroy(cache);
        md_gto_iso_field_destroy(field);
        iso_dens_prod_free(&p);
    }
    iso_dens_src_free(&s);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// Several fields in one set give each field exactly the mesh it gets alone.
UTEST(gto_iso_gpu_lib, fields_batch_like_alone) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    iso_dens_src_t s;
    ASSERT_TRUE(iso_dens_src_load(&s, "amide.h5"));
    const size_t lumo = vlx_test_lumo_idx(&s.t, STR_LIT("orbital/alpha/occupation"));
    enum { NF = 5 };
    md_gto_iso_field_t* fields[NF];
    double iso[NF];
    for (int k = 0; k < 4; ++k) { fields[k] = md_gto_iso_field_create(&s.basis, s.xyz, 0, s.C + (lumo - 2 + k) * s.num_ao, md_get_heap_allocator()); iso[k] = 0.05; }
    iso_dens_prod_t p;
    iso_dens_products(&p, &s, ISO_D_TRANS);
    fields[4] = iso_dens_field(&s, &p, false);
    iso[4] = 0.004;
    for (int k = 0; k < NF; ++k) ASSERT_TRUE(fields[k] != NULL);
    const double h = 0.25;
    md_gto_iso_gpu_set_t* all = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, (const md_gto_iso_field_t* const*)fields, NF, md_get_heap_allocator());
    ASSERT_TRUE(all != NULL);
    md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 4};
    for (int k = 0; k < NF; ++k) gd.iso[k] = iso[k];
    md_gto_iso_gpu_mesh_t gm;
    ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, all, &gd, &gm));
    printf("  %d fields together: %u bricks (union), %u corner slots, %u vertices\n", NF, gm.num_bricks, gm.num_corner_slots, gm.vertex_count);
    uint32_t sum_slots = 0;
    for (int k = 0; k < NF; ++k) {
        const md_gto_iso_field_t* one[1] = {fields[k]};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, one, 1, md_get_heap_allocator());
        md_gto_iso_gpu_desc_t d1 = {.cell_size = h, .project_iterations = 4};
        d1.iso[0] = iso[k];
        md_gto_iso_gpu_mesh_t m1;
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &d1, &m1));
        sum_slots += m1.num_bricks;
        printf("    field %d alone: %u bricks, vertices %u + %u (together %u + %u)\n", k, m1.num_bricks,
            m1.vertex_count_of[0][0], m1.vertex_count_of[0][1], gm.vertex_count_of[k][0], gm.vertex_count_of[k][1]);
        // the union root can be larger than a field's own, which moves its lattice: compare on the same lattice
        const bool same = m1.cells_per_axis == gm.cells_per_axis && m1.origin[0] == gm.origin[0] && m1.origin[1] == gm.origin[1] && m1.origin[2] == gm.origin[2];
        if (same) {
            EXPECT_EQ(m1.vertex_count_of[0][0], gm.vertex_count_of[k][0]);
            EXPECT_EQ(m1.vertex_count_of[0][1], gm.vertex_count_of[k][1]);
            EXPECT_EQ(m1.index_count_of[0][0], gm.index_count_of[k][0]);
            EXPECT_EQ(m1.index_count_of[0][1], gm.index_count_of[k][1]);
        } else {
            printf("      (another lattice: the union's root is larger)\n");
        }
        md_gto_iso_gpu_mesh_free(g.s, &m1);
        md_gto_iso_gpu_set_destroy(g.s, set);
    }
    EXPECT_EQ(gm.num_corner_slots <= sum_slots, true);
    md_gto_iso_gpu_mesh_free(g.s, &gm);
    md_gto_iso_gpu_set_destroy(g.s, all);
    for (int k = 0; k < NF; ++k) md_gto_iso_field_destroy(fields[k]);
    iso_dens_prod_free(&p);
    iso_dens_src_free(&s);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}


// A mesh's triangles by their corners' positions, each from its least corner (the winding kept), sorted:
// the same for two meshes that differ only in the order of their vertices and triangles (the outputs'
// order within a (field, sheet) range comes from atomics)
typedef struct iso_tri9_t { float v[9]; } iso_tri9_t;
static int iso_tri9_cmp(const void* pa, const void* pb) {
    const float* x = ((const iso_tri9_t*)pa)->v;
    const float* y = ((const iso_tri9_t*)pb)->v;
    for (int i = 0; i < 9; ++i) if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
    return 0;
}
static iso_tri9_t* iso_canon_tris(const float* pos4, const uint32_t* idx, uint32_t ni) {
    const uint32_t nt = ni / 3;
    iso_tri9_t* t = (iso_tri9_t*)malloc(sizeof(iso_tri9_t) * (nt + 1));
    for (uint32_t i = 0; i < nt; ++i) {
        int r = 0;      // the least corner first
        for (int c = 1; c < 3; ++c) {
            const float* pc = pos4 + 4 * idx[3 * i + c];
            const float* pr = pos4 + 4 * idx[3 * i + r];
            if (pc[0] < pr[0] || (pc[0] == pr[0] && (pc[1] < pr[1] || (pc[1] == pr[1] && pc[2] < pr[2])))) r = c;
        }
        for (int c = 0; c < 3; ++c) {
            const float* pc = pos4 + 4 * idx[3 * i + (r + c) % 3];
            for (int k = 0; k < 3; ++k) t[i].v[3 * c + k] = pc[k];
        }
    }
    qsort(t, nt, sizeof(iso_tri9_t), iso_tri9_cmp);
    return t;
}

// The extraction as a job (md_gto_iso_gpu_job_*): polled without waiting, with other work on the stream
// between the polls (a temp scope, a fill), it gives the mesh md_gto_iso_gpu_extract gives, to the bit;
// its host copy (desc.read_back) is that mesh; a job ended before it is done keeps no device memory.
UTEST(gto_iso_gpu_lib, job_matches_extract) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    iso_dens_src_t s;
    ASSERT_TRUE(iso_dens_src_load(&s, "amide.h5"));
    static const iso_gl_case_t cases[] = {
        {-1,          "amide_homo",  0.05},
        {ISO_D_SCF,   "amide_scf",   0.02},     // flooded (marked below): a readback per round
        {ISO_D_DIFF,  "amide_diff",  0.004},
    };
    for (int ci = 0; ci < 3; ++ci) {
        iso_dens_prod_t p;
        md_gto_iso_field_t* field = iso_gl_field(&s, cases[ci].kind, &p);
        ASSERT_TRUE(field != NULL);
        md_gto_iso_field_set_maxima_at_atoms(field, cases[ci].kind == ISO_D_SCF);
        const md_gto_iso_field_t* fl[1] = {field};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
        ASSERT_TRUE(set != NULL);
        md_gto_iso_gpu_desc_t gd = {.cell_size = 0.2, .project_iterations = 8, .relax_iterations = 3};
        gd.iso[0] = cases[ci].iso;

        md_gto_iso_gpu_mesh_t ref;
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &gd, &ref));
        const uint32_t nv = ref.vertex_count, ni = ref.index_count;
        float* rp = (float*)malloc(sizeof(float) * 4 * (nv + 1));
        float* rn = (float*)malloc(sizeof(float) * 4 * (nv + 1));
        uint32_t* ri = (uint32_t*)malloc(sizeof(uint32_t) * (ni + 1));
        ASSERT_TRUE(md_gto_iso_gpu_mesh_read_all(g.s, &ref, rp, rn, ri));

        md_gpu_sync_t sync = md_gpu_stream_record(g.s);
        md_gpu_sync_wait(sync);
        md_gpu_memory_stats_t m0;
        md_gpu_memory_stats(g.dev, MD_GPU_MEM_DEVICE, &m0);

        md_gto_iso_gpu_desc_t jd = gd;
        jd.read_back = true;
        md_gto_iso_gpu_job_t* job = md_gto_iso_gpu_job_begin(g.s, set, &jd);
        ASSERT_TRUE(job != NULL);
        md_gpu_mem_t scratch = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, 4096);
        uint32_t polls = 0;
        float prog = 0.0f;
        md_gto_iso_gpu_job_status_t st;
        while ((st = md_gto_iso_gpu_job_poll(job, 0.0)) == MD_GTO_ISO_GPU_JOB_RUNNING) {
            polls++;
            const float pr = md_gto_iso_gpu_job_progress(job);
            EXPECT_TRUE(pr >= prog);    // never goes back
            prog = pr;
            // other work on the stream between the polls
            md_gpu_temp_t t = md_gpu_temp_begin(g.s);
            md_gpu_mem_t tm = md_gpu_temp_alloc(g.s, MD_GPU_MEM_DEVICE, 1 << 16);
            if (tm.gpu) md_gpu_memset(g.s, tm.gpu, 0x5A, 1 << 16);
            md_gpu_temp_end(g.s, t);
            md_gpu_memset(g.s, scratch.gpu, 0, 4096);
            md_gpu_stream_flush(g.s);
        }
        ASSERT_EQ((int)st, (int)MD_GTO_ISO_GPU_JOB_DONE);
        EXPECT_EQ(md_gto_iso_gpu_job_progress(job), 1.0f);
        md_gto_iso_gpu_mesh_t jm;
        ASSERT_TRUE(md_gto_iso_gpu_job_end(job, &jm));
        printf("  %-12s %6u vertices: %u polls\n", cases[ci].name, nv, polls);
        EXPECT_TRUE(polls > 0);
        ASSERT_EQ(jm.vertex_count, nv);
        ASSERT_EQ(jm.index_count, ni);
        EXPECT_EQ(memcmp(jm.vertex_offset, ref.vertex_offset, sizeof(ref.vertex_offset)), 0);
        EXPECT_EQ(memcmp(jm.index_count_of, ref.index_count_of, sizeof(ref.index_count_of)), 0);
        if (nv > 0) {
            ASSERT_TRUE(jm.host_position != NULL && jm.host_normal != NULL && jm.host_index != NULL);
            // as the device has it (the job's own copy), and the same triangles as the blocking extraction's
            float* dp = (float*)malloc(sizeof(float) * 4 * (nv + 1));
            uint32_t* di = (uint32_t*)malloc(sizeof(uint32_t) * (ni + 1));
            ASSERT_TRUE(md_gto_iso_gpu_mesh_read_all(g.s, &jm, dp, NULL, di));
            EXPECT_EQ(memcmp(jm.host_position, dp, sizeof(float) * 4 * nv), 0);
            EXPECT_EQ(memcmp(jm.host_index, di, sizeof(uint32_t) * ni), 0);
            iso_tri9_t* ta = iso_canon_tris(rp, ri, ni);
            iso_tri9_t* tb = iso_canon_tris(jm.host_position, jm.host_index, ni);
            uint32_t diff = 0;
            for (uint32_t i = 0; i < ni / 3; ++i) diff += iso_tri9_cmp(&ta[i], &tb[i]) != 0;
            EXPECT_EQ(diff, 0u);
            if (diff) printf("    %u of %u triangles differ\n", diff, ni / 3);
            free(ta); free(tb); free(dp); free(di);
        }
        md_gto_iso_gpu_mesh_free(g.s, &jm);

        // in chunks (a long job's launches, each chunk followed by a progress marker): the same triangles, the
        // progress in more steps
        {
            const int saved = isog_chunks;
            isog_chunks = 6;
            md_gto_iso_gpu_job_t* jc = md_gto_iso_gpu_job_begin(g.s, set, &jd);
            ASSERT_TRUE(jc != NULL);
            float seen[256];
            uint32_t nseen = 0;
            while (md_gto_iso_gpu_job_poll(jc, 0.0) == MD_GTO_ISO_GPU_JOB_RUNNING) {
                const float pr = md_gto_iso_gpu_job_progress(jc);
                if ((nseen == 0 || pr != seen[nseen - 1]) && nseen < 256) seen[nseen++] = pr;
            }
            isog_chunks = saved;
            md_gto_iso_gpu_mesh_t cm;
            ASSERT_TRUE(md_gto_iso_gpu_job_end(jc, &cm));
            printf("    in 6 chunks: %u progress values seen:", nseen);
            for (uint32_t i = 0; i < nseen; ++i) printf(" %.2f", seen[i]);
            printf("\n");
            EXPECT_TRUE(nseen >= 8);
            ASSERT_EQ(cm.vertex_count, nv);
            ASSERT_EQ(cm.index_count, ni);
            if (nv > 0) {
                iso_tri9_t* ta = iso_canon_tris(rp, ri, ni);
                iso_tri9_t* tb = iso_canon_tris(cm.host_position, cm.host_index, ni);
                uint32_t diff = 0;
                for (uint32_t i = 0; i < ni / 3; ++i) diff += iso_tri9_cmp(&ta[i], &tb[i]) != 0;
                EXPECT_EQ(diff, 0u);
                free(ta); free(tb);
            }
            md_gto_iso_gpu_mesh_free(g.s, &cm);
        }

        // ended early, at each of its first few waits: nothing kept
        for (int stop = 0; stop < 4; ++stop) {
            md_gto_iso_gpu_job_t* jb = md_gto_iso_gpu_job_begin(g.s, set, &jd);
            ASSERT_TRUE(jb != NULL);
            for (int k = 0; k < stop; ++k) {
                md_gpu_sync_wait(md_gpu_stream_record(g.s));
                md_gto_iso_gpu_job_poll(jb, 0.0);
            }
            md_gto_iso_gpu_mesh_t none;
            const bool handed = md_gto_iso_gpu_job_end(jb, &none);
            if (handed) md_gto_iso_gpu_mesh_free(g.s, &none);
            else EXPECT_EQ(none.position, (md_gpu_addr_t)0);
        }
        md_gpu_free(g.s, scratch.gpu);
        md_gpu_sync_wait(md_gpu_stream_record(g.s));
        md_gpu_device_poll(g.dev);
        md_gpu_memory_stats_t m1;
        md_gpu_memory_stats(g.dev, MD_GPU_MEM_DEVICE, &m1);
        EXPECT_EQ(m1.allocations, m0.allocations);
        EXPECT_EQ(m1.bytes_in_use, m0.bytes_in_use);

        free(rp); free(rn); free(ri);
        md_gto_iso_gpu_mesh_free(g.s, &ref);
        md_gto_iso_gpu_set_destroy(g.s, set);
        md_gto_iso_field_destroy(field);
        iso_dens_prod_free(&p);
    }
    iso_dens_src_free(&s);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// One field given twice at two iso values (a density at two thresholds): the copies share their functions
// and each gives what it gives in a set of its own on the same lattice
UTEST(gto_iso_gpu_lib, repeated_field_shares_functions) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    iso_dens_src_t s;
    ASSERT_TRUE(iso_dens_src_load(&s, "amide.h5"));
    iso_dens_prod_t p;
    iso_dens_products(&p, &s, ISO_D_TRANS);
    md_gto_iso_field_t* f = iso_dens_field(&s, &p, false);
    ASSERT_TRUE(f != NULL);
    const double h = 0.25;
    // the smaller iso value decides the root box of both, so the single sets below see the same lattice
    const double iso[2] = {0.002, 0.004};
    const md_gto_iso_field_t* two[2] = {f, f};
    md_gto_iso_gpu_set_t* both = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, two, 2, md_get_heap_allocator());
    ASSERT_TRUE(both != NULL);
    md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 4};
    gd.iso[0] = iso[0];
    gd.iso[1] = iso[1];
    md_gto_iso_gpu_mesh_t gm;
    ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, both, &gd, &gm));
    for (int k = 0; k < 2; ++k) {
        const md_gto_iso_field_t* one[2] = {f, f};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, one, 2, md_get_heap_allocator());
        md_gto_iso_gpu_desc_t d1 = gd;
        d1.iso[1 - k] = iso[0];     // keeps the root of the pair
        d1.iso[k] = iso[k];
        md_gto_iso_gpu_mesh_t m1;
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &d1, &m1));
        EXPECT_EQ(m1.cells_per_axis, gm.cells_per_axis);
        EXPECT_EQ(m1.vertex_count_of[k][0], gm.vertex_count_of[k][0]);
        EXPECT_EQ(m1.vertex_count_of[k][1], gm.vertex_count_of[k][1]);
        EXPECT_EQ(m1.index_count_of[k][0], gm.index_count_of[k][0]);
        EXPECT_EQ(m1.index_count_of[k][1], gm.index_count_of[k][1]);
        md_gto_iso_gpu_mesh_free(g.s, &m1);
        md_gto_iso_gpu_set_destroy(g.s, set);
    }
    // and against a set holding the field once
    {
        const md_gto_iso_field_t* one[1] = {f};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, one, 1, md_get_heap_allocator());
        md_gto_iso_gpu_desc_t d1 = {.cell_size = h, .project_iterations = 4};
        d1.iso[0] = iso[0];
        md_gto_iso_gpu_mesh_t m1;
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &d1, &m1));
        EXPECT_EQ(m1.vertex_count_of[0][0], gm.vertex_count_of[0][0]);
        EXPECT_EQ(m1.index_count_of[0][1], gm.index_count_of[0][1]);
        printf("  repeated field: %u + %u vertices at %g and %g\n", gm.vertex_count_of[0][0] + gm.vertex_count_of[0][1],
               gm.vertex_count_of[1][0] + gm.vertex_count_of[1][1], iso[0], iso[1]);
        md_gto_iso_gpu_mesh_free(g.s, &m1);
        md_gto_iso_gpu_set_destroy(g.s, set);
    }
    // the whole mesh in one read, against the per field read
    {
        float*    p4 = (float*)malloc(sizeof(float) * 4 * gm.vertex_count);
        float*    n4 = (float*)malloc(sizeof(float) * 4 * gm.vertex_count);
        uint32_t* ix = (uint32_t*)malloc(sizeof(uint32_t) * gm.index_count);
        ASSERT_TRUE(md_gto_iso_gpu_mesh_read_all(g.s, &gm, p4, n4, ix));
        md_gto_iso_mesh_t m1;
        ASSERT_TRUE(md_gto_iso_gpu_mesh_read(g.s, &gm, 1, &m1, md_get_heap_allocator()));
        const uint32_t v0 = gm.vertex_offset[1][0], i0 = gm.index_offset[1][0];
        int bad = 0;
        for (uint32_t i = 0; i < m1.vertex_count; ++i) {
            for (int k = 0; k < 3; ++k) bad += m1.position[3 * i + k] != p4[4 * (v0 + i) + k] || m1.normal[3 * i + k] != n4[4 * (v0 + i) + k];
        }
        for (uint32_t i = 0; i < m1.index_count; ++i) bad += m1.index[i] + v0 != ix[i0 + i];
        EXPECT_EQ(0, bad);
        md_gto_iso_mesh_free(&m1);
        free(p4); free(n4); free(ix);
    }
    md_gto_iso_gpu_mesh_free(g.s, &gm);
    md_gto_iso_gpu_set_destroy(g.s, both);
    md_gto_iso_field_destroy(f);
    iso_dens_prod_free(&p);
    iso_dens_src_free(&s);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// GTO_ISO_BENCH=1: the GPU library path on whatever device md_gpu picks, against the production dense
// path at the same spacing (md_gto_gpu_density_launch / md_gto_gpu_orbital_launch into a volume over
// the tight surface box). The dense numbers are a volume only - no surface yet - so they are a lower
// bound for any volume-based isosurface. GTO_ISO_BENCH_H=0.2,0.15,... picks the spacings,
// GTO_ISO_BENCH_REPS the repetitions (median), GTO_ISO_BENCH_CPU=1 adds the CPU path (single thread).
typedef struct iso_bench_field_t {
    const char* file;
    const char* name;
    int         kind;       // -1 HOMO orbital, ISO_D_*, or -16: 16 orbitals around the gap (a grid of views)
    double      iso;
} iso_bench_field_t;

static double iso_bench_dense(iso_gpu_t* g, const iso_dens_src_t* s, const md_gto_iso_field_t* const* fields, const double* iso, int nfields,
                              const double* const* orb_rows, const double* D, double h, int reps, int out_dim[3]) {
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    for (int k = 0; k < nfields; ++k) {
        double bl[3], bh[3];
        iso_surface_box(&fields[k]->f, iso[k], bl, bh);
        for (int a = 0; a < 3; ++a) { lo[a] = fmin(lo[a], bl[a]); hi[a] = fmax(hi[a], bh[a]); }
    }
    int dd[3];
    for (int a = 0; a < 3; ++a) dd[a] = (int)ceil((hi[a] - lo[a]) / h) + 1;
    memcpy(out_dim, dd, sizeof(dd));
    if ((double)dd[0] * dd[1] * dd[2] > (double)(1u << 27)) return NAN;
    const size_t nao = s->num_ao;
    md_gto_gpu_basis_t gb = md_gto_gpu_basis_create(g->s, &(md_gto_gpu_basis_desc_t){ .basis = &s->basis, .cutoff = 1.0e-6 });
    md_gpu_mem_t atoms = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, md_gto_gpu_atom_buffer_size(s->num_atoms));
    float* p = (float*)md_gpu_upload_begin(g->s, atoms.gpu, md_gto_gpu_atom_buffer_size(s->num_atoms));
    md_gto_gpu_atom_pack(p, s->xyz, 0, s->num_atoms);
    md_gpu_upload_end(g->s);
    const size_t cbytes = D ? md_gto_gpu_coeff_size_density(nao) : md_gto_gpu_coeff_size_mo((size_t)nfields, nao);
    md_gpu_mem_t coeff = md_gpu_malloc(g->s, MD_GPU_MEM_DEVICE, cbytes);
    p = (float*)md_gpu_upload_begin(g->s, coeff.gpu, cbytes);
    if (D) md_gto_gpu_coeff_pack_density(p, D, nao);
    else   md_gto_gpu_coeff_pack_mo(p, orb_rows, NULL, (size_t)nfields, nao);
    md_gpu_upload_end(g->s);
    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin  = vec3_set((float)lo[0], (float)lo[1], (float)lo[2]),
        .spacing = vec3_set((float)h, (float)h, (float)h),
        .dim     = {dd[0], dd[1], dd[2]},
    };
    const int ntex = D ? 1 : nfields;
    md_gpu_texture_t tex[16] = {0};
    bool ok = true;
    for (int k = 0; k < ntex; ++k) {
        tex[k] = md_gpu_texture_create(g->s, &(md_gpu_texture_desc_t){
            .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
            .width = (uint32_t)dd[0], .height = (uint32_t)dd[1], .depth_or_layers = (uint32_t)dd[2],
        });
        ok &= tex[k] != NULL;
    }
    double med = NAN;
    if (ok) {
        double dt[64];
        for (int r = 0; r <= reps && r < 64; ++r) {
            md_tick_t t0 = md_tick_now();
            for (int k = 0; k < ntex; ++k) {
                if (D) {
                    md_gto_gpu_density_desc_t dd_ = {.basis = gb, .atom_xyz = atoms.gpu, .coeff = coeff.gpu, .out_tex = tex[k], .grid = &grid, .op = MD_GTO_OP_SET};
                    md_gto_gpu_density_launch(g->s, &dd_);
                } else {
                    // one orbital per volume: its row of the coefficient buffer
                    md_gto_gpu_orbital_desc_t od = {.basis = gb, .atom_xyz = atoms.gpu, .coeff = coeff.gpu + sizeof(float) * nao * (size_t)k,
                        .out_tex = tex[k], .grid = &grid, .num_orbitals = 1, .eval_mode = MD_GTO_EVAL_MODE_PSI, .op = MD_GTO_OP_SET};
                    md_gto_gpu_orbital_launch(g->s, &od);
                }
            }
            md_gpu_stream_sync(g->s);
            if (r > 0) dt[r - 1] = iso_seconds_since(t0) * 1e3;   // the first is a warm-up
        }
        med = iso_median(dt, reps < 64 ? reps : 63);
    }
    for (int k = 0; k < ntex; ++k) if (tex[k]) md_gpu_texture_destroy(tex[k]);
    md_gpu_free(g->s, atoms.gpu);
    md_gpu_free(g->s, coeff.gpu);
    md_gto_gpu_basis_destroy(g->s, gb);
    return med;
}

UTEST(gto_iso_gpu_lib, benchmark) {
    const char* on = getenv("GTO_ISO_BENCH");
    if (!on || on[0] != '1') UTEST_SKIP("set GTO_ISO_BENCH=1");
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gpu_device_info_t info;
    md_gpu_device_info(g.dev, &info);
    md_gto_gpu_initialize(g.dev);
    md_gto_iso_gpu_initialize(g.dev);
    double hs[16];
    int nh = iso_parse_list(getenv("GTO_ISO_BENCH_H"), hs, 16);
    if (nh == 0) { hs[0] = 0.2; hs[1] = 0.15; hs[2] = 0.1; nh = 3; }
    const int reps = getenv("GTO_ISO_BENCH_REPS") ? atoi(getenv("GTO_ISO_BENCH_REPS")) : 5;
    const bool with_cpu = getenv("GTO_ISO_BENCH_CPU") && getenv("GTO_ISO_BENCH_CPU")[0] == '1';
    // BENCHMARK (temporary): one extra run per case with the step by step profile; off with GTO_ISO_BENCH_PROFILE=0
    const bool with_profile = !(getenv("GTO_ISO_BENCH_PROFILE") && getenv("GTO_ISO_BENCH_PROFILE")[0] == '0');
    // only the cases whose name contains this (e.g. "scf", "mol ")
    const char* only = getenv("GTO_ISO_BENCH_ONLY");
    static md_gto_iso_gpu_profile_t prof;
    static char prof_text[32768];

    static const iso_bench_field_t cases[] = {
        {"amide.h5", "amide homo",   -1,          0.05},
        {"amide.h5", "amide scf",    ISO_D_SCF,   0.002},
        {"amide.h5", "amide diff",   ISO_D_DIFF,  0.002},
        {"amide.h5", "amide trans",  ISO_D_TRANS, 0.002},
        {"mol.h5",   "mol homo",     -1,          0.05},
        {"mol.h5",   "mol scf",      ISO_D_SCF,   0.002},
        {"mol.h5",   "mol diff",     ISO_D_DIFF,  0.002},
        {"mol.h5",   "mol trans",    ISO_D_TRANS, 0.002},
        {"mol.h5",   "mol 16 MOs",   -16,         0.05},
        {"tq.scf.results.h5",   "tq homo",    -1,        0.05},
        {"tq.scf.results.h5",   "tq scf",     ISO_D_SCF, 0.002},
        {"tq.scf.results.h5",   "tq 16 MOs",  -16,       0.05},
    };
    printf("\n  device: %s, %s | display meshes (project 8, relax 3) | median of %d\n", info.name, info.driver, reps);
    printf("  %-12s %5s | %4s %6s %6s %7s | %8s %7s %7s %7s %7s %7s %7s | %9s %8s %13s", "field", "h", "fn", "bricks", "slots", "verts",
        "total ms", "octree", "corners", "nets", "project", "relax", "output", "dense ms", "dense/tot", "dense grid");
    if (with_cpu) printf(" | %8s", "CPU s");
    printf("\n");

    for (int ci = 0; ci < (int)(sizeof(cases) / sizeof(cases[0])); ++ci) {
        const iso_bench_field_t* bc = &cases[ci];
        if (only && !strstr(bc->name, only)) continue;
        iso_dens_src_t s;
        const bool need_nto = bc->kind >= 0 && bc->kind != ISO_D_SCF;
        if (!iso_dens_src_load_ex(&s, bc->file, need_nto)) { printf("  %-12s (could not load %s)\n", bc->name, bc->file); iso_dens_src_free(&s); continue; }
        if (ci == 0 || strcmp(cases[ci - 1].file, bc->file) != 0 || (only && !strstr(cases[ci - 1].name, only))) {
            printf("  -- %s: %zu atoms, %u shells, %zu AOs, %zu MOs, %zu NTO pairs\n", bc->file, s.num_atoms, s.basis.num_shells, s.num_ao, s.num_mo, s.num_nto);
        }
        const size_t lumo = vlx_test_lumo_idx(&s.t, STR_LIT("orbital/alpha/occupation"));
        md_gto_iso_field_t* fields[16] = {0};
        const double* rows[16] = {0};
        double iso[16];
        int nf = 0;
        iso_dens_prod_t p = {0};
        double* D = NULL;
        if (bc->kind == -1 || bc->kind == -16) {
            nf = bc->kind == -1 ? 1 : 16;
            for (int k = 0; k < nf; ++k) {
                const size_t mo = bc->kind == -1 ? lumo - 1 : lumo - 8 + (size_t)k;
                rows[k] = s.C + mo * s.num_ao;
                fields[k] = md_gto_iso_field_create(&s.basis, s.xyz, 0, rows[k], md_get_heap_allocator());
                iso[k] = bc->iso;
            }
        } else {
            nf = 1;
            iso_dens_products(&p, &s, bc->kind);
            fields[0] = iso_dens_field(&s, &p, false);
            iso[0] = bc->iso;
            D = iso_dens_matrix(&p, s.num_ao);
            for (size_t i = 0; i < s.num_ao; ++i) for (size_t j = i + 1; j < s.num_ao; ++j) {
                const double v = 0.5 * (D[i * s.num_ao + j] + D[j * s.num_ao + i]);
                D[i * s.num_ao + j] = D[j * s.num_ao + i] = v;
            }
        }
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, (const md_gto_iso_field_t* const*)fields, (uint32_t)nf, md_get_heap_allocator());
        if (!set) { printf("  %-12s (set creation failed)\n", bc->name); continue; }
        uint32_t nfn = 0;
        for (int k = 0; k < nf; ++k) { md_gto_iso_field_info_t fi; md_gto_iso_field_info(&fi, fields[k]); nfn += fi.num_functions; }

        for (int hi = 0; hi < nh; ++hi) {
            const double h = hs[hi];
            md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 8, .relax_iterations = 3};
            for (int k = 0; k < nf; ++k) gd.iso[k] = iso[k];
            md_gto_iso_gpu_mesh_t gm;
            if (!md_gto_iso_gpu_extract(g.s, set, &gd, &gm)) { printf("  %-12s %5.3f (failed)\n", bc->name, h); continue; }
            md_gto_iso_gpu_mesh_free(g.s, &gm);
            double tot[64], st[MD_GTO_ISO_GPU_STAGE_COUNT][64];
            const int nr = reps < 64 ? reps : 64;
            for (int r = 0; r < nr; ++r) {
                md_gto_iso_gpu_extract(g.s, set, &gd, &gm);
                tot[r] = gm.seconds_total * 1e3;
                md_gto_iso_gpu_mesh_free(g.s, &gm);
            }
            gd.time_stages = true;
            for (int r = 0; r < nr; ++r) {
                md_gto_iso_gpu_extract(g.s, set, &gd, &gm);
                for (int k = 0; k < MD_GTO_ISO_GPU_STAGE_COUNT; ++k) st[k][r] = gm.seconds[k] * 1e3;
                if (r + 1 < nr) md_gto_iso_gpu_mesh_free(g.s, &gm);
            }
            int dd[3] = {0};
            const double dense = iso_bench_dense(&g, &s, (const md_gto_iso_field_t* const*)fields, iso, nf, rows, D, h, nr, dd);
            const double t = iso_median(tot, nr);
            printf("  %-12s %5.3f | %4u %6u %6u %7u | %8.2f", bc->name, h, nfn, gm.num_bricks, gm.num_corner_slots, gm.vertex_count, t);
            for (int k = 0; k < MD_GTO_ISO_GPU_STAGE_COUNT; ++k) printf(" %7.2f", iso_median(st[k], nr));
            if (isnan(dense)) printf(" | %9s %8s %13s", "-", "-", "-");
            else printf(" | %9.2f %7.2fx %4dx%3dx%3d", dense, dense / t, dd[0], dd[1], dd[2]);
            md_gto_iso_gpu_mesh_free(g.s, &gm);
            if (with_profile) {
                md_gto_iso_gpu_desc_t pd = gd;
                pd.time_stages = false;
                md_gto_iso_gpu_profile_set(&prof);
                const bool pok = md_gto_iso_gpu_extract(g.s, set, &pd, &gm);
                md_gto_iso_gpu_profile_set(NULL);
                if (pok) md_gto_iso_gpu_mesh_free(g.s, &gm);
                md_gto_iso_gpu_profile_print(&prof, prof_text, sizeof(prof_text));
            }
            if (with_cpu && nf == 1) {
                md_tick_t t0 = md_tick_now();
                md_gto_iso_cache_t* cache = md_gto_iso_cache_create(fields[0], &(md_gto_iso_cache_desc_t){.cell_size = h, .iso_min = iso[0]}, md_get_heap_allocator());
                md_gto_iso_mesh_t cm;
                md_gto_iso_mesh_extract(&cm, cache, &(md_gto_iso_mesh_desc_t){.iso = iso[0], .project_iterations = 8, .relax_iterations = 3}, md_get_heap_allocator());
                printf(" | %8.2f", iso_seconds_since(t0));
                md_gto_iso_mesh_free(&cm);
                md_gto_iso_cache_destroy(cache);
            }
            printf("\n");
            if (with_profile) printf("%s", prof_text);
        }
        md_gto_iso_gpu_set_destroy(g.s, set);
        for (int k = 0; k < nf; ++k) md_gto_iso_field_destroy(fields[k]);
        free(D);
        iso_dens_prod_free(&p);
        iso_dens_src_free(&s);
    }
    md_gto_iso_gpu_shutdown();
    md_gto_gpu_shutdown();
    iso_gpu_free(&g);
}

// ---------------------------------------------------------------------------------------------
// EXPERIMENT (GTO_ISO_TM=1): a Taylor model of rho per brick, as TOPO_ANALYTIC builds per cube
// (md_topo.c, cpg_box_eval / cpg_shell_rem): each AO expanded to order K-1 at the brick centre with a
// rigorous separable remainder E_K over the brick, rho's polynomial sum D T T exact (D multiplies exact
// centre data only), the remainder 2 sum |(D T)_mu| E_mu + sum |D| E E. Per lattice point of the brick
// the enclosure P(d) +- R; a brick whose 125 points all lie on one side of the iso value has no
// surface cell. Counts, over the bricks the GPU gave corners, which the model would drop, by what
// their corners show (surface: must never be dropped), in double on the CPU.
// ---------------------------------------------------------------------------------------------

#define TM_LMAX  4
#define TM_SEPA  7                     // 1D derivative orders 0..6
#define TM_EXTRA 2                     // exact 1D terms beyond the remainder's order

typedef struct tm_ctx_t {
    double kap[TM_LMAX + 1][TM_SEPA][TM_LMAX + TM_SEPA];   // d^a [x^i E] = sum_m kap[i][a][m] alpha^((m-i+a)/2) x^m E
    int    nmi;                                             // multi-indices of order <= K-1
    int    mi[20][3];
    double mi_fact[20];
} tm_ctx_t;

static void tm_init(tm_ctx_t* t, int K) {
    memset(t, 0, sizeof(*t));
    for (int i = 0; i <= TM_LMAX; ++i) {
        t->kap[i][0][i] = 1.0;
        for (int a = 0; a + 1 < TM_SEPA; ++a) {
            for (int m = 0; m <= i + a + 1 && m < TM_LMAX + TM_SEPA; ++m) {
                double v = (m + 1 <= i + a && m + 1 < TM_LMAX + TM_SEPA) ? (m + 1) * t->kap[i][a][m + 1] : 0.0;
                if (m >= 1) v -= 2.0 * t->kap[i][a][m - 1];
                t->kap[i][a + 1][m] = v;
            }
        }
    }
    static const double fact[4] = {1, 1, 2, 6};
    int n = 0;
    for (int o = 0; o < K; ++o)
        for (int i = o; i >= 0; --i)
            for (int j = o - i; j >= 0; --j) {
                const int k = o - i - j;
                t->mi[n][0] = i; t->mi[n][1] = j; t->mi[n][2] = k;
                t->mi_fact[n] = fact[i] * fact[j] * fact[k];
                n++;
            }
    t->nmi = n;
}

// 1D: centre derivatives val[i][a] (a <= TM_SEPA - 2) of x^i exp(-alpha x^2) at xc, and sups sup[i][a]
// over [xc - h, xc + h] for a >= 3 (as cpg_sep_1d)
static void tm_sep_1d(const tm_ctx_t* t, int l, double alpha, double xc, double h, double val[TM_LMAX + 1][TM_SEPA], double sup[TM_LMAX + 1][TM_SEPA]) {
    double g[TM_LMAX + TM_SEPA][TM_SEPA];
    const int imax = l + TM_SEPA - 2;
    const double E = exp(-alpha * xc * xc);
    double xp = 1.0;
    for (int i = 0; i <= imax; ++i) { g[i][0] = xp * E; xp *= xc; }
    for (int a = 0; a < TM_SEPA - 2; ++a)
        for (int i = 0; i < imax - a; ++i) {
            double v = -2.0 * alpha * g[i + 1][a];
            if (i >= 1) v += i * g[i - 1][a];
            g[i][a + 1] = v;
        }
    for (int i = 0; i <= l; ++i) for (int a = 0; a < TM_SEPA - 1; ++a) val[i][a] = g[i][a];
    const double x0 = xc - h, x1 = xc + h, ax0 = fabs(x0), ax1 = fabs(x1);
    const double tlo = (x0 <= 0.0 && x1 >= 0.0) ? 0.0 : fmin(ax0, ax1);
    const double thi = fmax(ax0, ax1);
    const double elo = exp(-alpha * tlo * tlo), ehi = exp(-alpha * thi * thi);
    const int mmax = l + TM_SEPA - 1;
    double Sm[TM_LMAX + TM_SEPA + 1], plo = 1.0, phi = 1.0;
    for (int m = 0; m <= mmax; ++m) {
        const double t2 = m / (2.0 * alpha);
        if (t2 <= tlo * tlo)      Sm[m] = plo * elo;
        else if (t2 >= thi * thi) Sm[m] = phi * ehi;
        else                      Sm[m] = pow(t2, 0.5 * m) * exp(-0.5 * m);
        plo *= tlo;
        phi *= thi;
    }
    for (int i = 0; i <= l; ++i)
        for (int a = 3; a < TM_SEPA; ++a) {
            double acc = 0.0;
            for (int m = (i + a) & 1; m <= i + a; m += 2) {
                const double k = t->kap[i][a][m];
                if (k != 0.0) acc += fabs(k) * pow(alpha, (m - i + a) / 2) * Sm[m];
            }
            sup[i][a] = acc;
        }
}

// The shell's AOs at the centre c: Taylor coefficients V[ci][mi] (derivatives of order <= K - 1) and the
// remainder bound E[ci] of |phi - T_{K-1} phi| over the cube c +- h (separable, as cpg_shell_rem)
static void tm_shell(const tm_ctx_t* t, const iso_ao_shell_t* s, const double c[3], double h, int K, double* V, double* E) {
    const int P = K - 1 + TM_EXTRA;
    static const double inv_fact[TM_SEPA] = {1.0, 1.0, 1.0 / 2, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720};
    for (int ci = 0; ci < s->ncart; ++ci) { E[ci] = 0.0; for (int m = 0; m < t->nmi; ++m) V[ci * 20 + m] = 0.0; }
    double val[3][TM_LMAX + 1][TM_SEPA], sup[3][TM_LMAX + 1][TM_SEPA];
    for (int q = 0; q < s->nprim; ++q) {
        const double al = s->alpha[q], cf = s->coef[q];
        for (int k = 0; k < 3; ++k) tm_sep_1d(t, s->l, al, c[k] - s->c[k], h, val[k], sup[k]);
        for (int ci = 0; ci < s->ncart; ++ci) {
            const int ii[3] = {s->ijk[ci][0], s->ijk[ci][1], s->ijk[ci][2]};
            const double w = cf * s->nrm[ci];
            for (int m = 0; m < t->nmi; ++m)
                V[ci * 20 + m] += w * val[0][ii[0]][t->mi[m][0]] * val[1][ii[1]][t->mi[m][1]] * val[2][ii[2]][t->mi[m][2]];
            // u_k[m]: |degree-m term| of axis k, the remainder at degree P + 1; sum over m0 + m1 + m2 >= K
            double u[3][TM_SEPA + 1];
            for (int k = 0; k < 3; ++k) {
                for (int m = 0; m <= P; ++m) u[k][m] = fabs(val[k][ii[k]][m]) * pow(h, m) * inv_fact[m];
                u[k][P + 1] = sup[k][ii[k]][P + 1] * pow(h, P + 1) * inv_fact[P + 1];
            }
            double b = 0.0;
            for (int m0 = 0; m0 <= P + 1; ++m0) for (int m1 = 0; m1 <= P + 1; ++m1) for (int m2 = 0; m2 <= P + 1; ++m2)
                if (m0 + m1 + m2 >= K) b += u[0][m0] * u[1][m1] * u[2][m2];
            E[ci] += fabs(w) * b;
        }
    }
}

UTEST(gto_iso_gpu_lib, taylor_brick_experiment) {
    if (!getenv("GTO_ISO_TM")) UTEST_SKIP("set GTO_ISO_TM=1 (EXPERIMENT)");
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    const char* files[] = {"amide.h5", "mol.h5", "tq.scf.results.h5"};
    const double h = getenv("GTO_ISO_TM_H") ? atof(getenv("GTO_ISO_TM_H")) : 0.2;
    const double iso = getenv("GTO_ISO_TM_ISO") ? atof(getenv("GTO_ISO_TM_ISO")) : 0.0025;
    for (int fi = 0; fi < 3; ++fi) {
        iso_dens_src_t s;
        if (!iso_dens_src_load(&s, files[fi])) { printf("  %s: not loaded\n", files[fi]); continue; }
        iso_dens_prod_t p;
        md_gto_iso_field_t* field = iso_gl_field(&s, ISO_D_SCF, &p);
        ASSERT_TRUE(field != NULL);
        const iso_density_t* d = field->f.dens;
        const int nao = d->num_ao;
        // D without the AO normalizations (the AO values below carry them)
        double* D = (double*)calloc((size_t)nao * nao, sizeof(double));
        for (int tt = 0; tt < d->num_terms; ++tt) {
            const iso_term_t* tm = &d->terms[tt];
            if (tm->b < 0) continue;
            const double w = tm->a == tm->b ? tm->w : 0.5 * tm->w;
            for (int i = 0; i < nao; ++i) {
                const double ca = d->C[(size_t)i * d->num_fn + tm->a], cb = d->C[(size_t)i * d->num_fn + tm->b];
                for (int j = 0; j < nao; ++j) {
                    const double da = d->C[(size_t)j * d->num_fn + tm->a], db = d->C[(size_t)j * d->num_fn + tm->b];
                    D[(size_t)i * nao + j] += tm->a == tm->b ? w * ca * da : w * (ca * db + cb * da);
                }
            }
        }
        // the GPU's bricks and corners
        const md_gto_iso_field_t* fl[1] = {field};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
        ASSERT_TRUE(set != NULL);
        md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 8, .relax_iterations = 3};
        gd.iso[0] = iso;
        isog_debug_t dbg = {0};
        isog_debug = &dbg;
        md_gto_iso_gpu_mesh_t gm;
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &gd, &gm));
        isog_debug = NULL;
        md_gto_iso_gpu_mesh_free(g.s, &gm);

        for (int K = 3; K <= 4; ++K) {
            tm_ctx_t t;
            tm_init(&t, K);
            const int nm = t.nmi;
            double* V = (double*)malloc(sizeof(double) * (size_t)nao * 20);
            double* E = (double*)malloc(sizeof(double) * (size_t)nao);
            double* Dv = (double*)malloc(sizeof(double) * (size_t)nao * 20);
            long cnt[3][2] = {{0}};          // [surface, inside, outside][kept, dropped by the model (box W)]
            long cntp[3] = {0};              // dropped with the remainder per point
            long bad = 0, outside_enc = 0;
            double rel_w = 0.0; long nrel = 0;
            for (uint32_t b = 0; b < dbg.num_bricks; ++b) {
                const uint32_t m = dbg.bricks[4 * b + 3];
                if (!(m & 1u)) continue;     // the + sheet of field 0
                const float* cv = dbg.corners + 125 * (size_t)dbg.slots[b];
                bool above = false, below = false;
                for (int k = 0; k < 125; ++k) { if (cv[k] - iso > 0.0) above = true; else below = true; }
                // surface: a cell with a sign change (corners straddle: some cell crosses)
                const int cls = (above && below) ? 0 : above ? 1 : 2;
                const double c[3] = {dbg.origin[0] + (dbg.bricks[4 * b] * 4.0 + 2.0) * h, dbg.origin[1] + (dbg.bricks[4 * b + 1] * 4.0 + 2.0) * h, dbg.origin[2] + (dbg.bricks[4 * b + 2] * 4.0 + 2.0) * h};
                const double hb = 2.0 * h;
                for (int sh = 0; sh < d->num_ao_shells; ++sh) tm_shell(&t, &d->ao_shells[sh], c, hb, K, V + (size_t)d->ao_shells[sh].ao * 20, E + d->ao_shells[sh].ao);
                // Dv = D V; DOT = V^T D V
                double DOT[20][20] = {{0}};
                for (int i = 0; i < nao; ++i) {
                    double acc[20] = {0};
                    const double* row = D + (size_t)i * nao;
                    for (int j = 0; j < nao; ++j) { const double dd = row[j]; if (dd == 0.0) continue; for (int a = 0; a < nm; ++a) acc[a] += dd * V[(size_t)j * 20 + a]; }
                    for (int a = 0; a < nm; ++a) Dv[(size_t)i * 20 + a] = acc[a];
                    for (int a = 0; a < nm; ++a) for (int bb = 0; bb < nm; ++bb) DOT[a][bb] += V[(size_t)i * 20 + a] * acc[bb];
                }
                double REE = 0.0, RW = 0.0;  // sum |D| E E; box: 2 sum E_mu W_mu
                for (int i = 0; i < nao; ++i) {
                    const double* row = D + (size_t)i * nao;
                    double s1 = 0.0;
                    for (int j = 0; j < nao; ++j) s1 += fabs(row[j]) * E[j];
                    REE += E[i] * s1;
                    double W = 0.0;
                    for (int a = 0; a < nm; ++a) {
                        const int o = t.mi[a][0] + t.mi[a][1] + t.mi[a][2];
                        W += fabs(Dv[(size_t)i * 20 + a]) * pow(hb, o) / t.mi_fact[a];
                    }
                    RW += 2.0 * E[i] * W;
                }
                bool drop_box = true, drop_pt = true;
                int side_box = 0, side_pt = 0;
                for (int k = 0; k < 125; ++k) {
                    const double dd[3] = {((k % 5) - 2) * h, (((k / 5) % 5) - 2) * h, ((k / 25) - 2) * h};
                    double mo[20];
                    for (int a = 0; a < nm; ++a) mo[a] = pow(dd[0], t.mi[a][0]) * pow(dd[1], t.mi[a][1]) * pow(dd[2], t.mi[a][2]) / t.mi_fact[a];
                    double P = 0.0;
                    for (int a = 0; a < nm; ++a) { double r = 0.0; for (int bb = 0; bb < nm; ++bb) r += DOT[a][bb] * mo[bb]; P += mo[a] * r; }
                    double Rp = REE;
                    for (int i = 0; i < nao; ++i) {
                        if (E[i] == 0.0) continue;
                        double w = 0.0;
                        for (int a = 0; a < nm; ++a) w += Dv[(size_t)i * 20 + a] * mo[a];
                        Rp += 2.0 * E[i] * fabs(w);
                    }
                    const double Rb = RW + REE;
                    // rigour: the corner (float) inside the enclosure
                    if (cv[k] < P - Rp - 1e-5 * fabs(P) - 1e-9 || cv[k] > P + Rp + 1e-5 * fabs(P) + 1e-9) outside_enc++;
                    const int sb = P - Rb > iso ? 1 : P + Rb < iso ? -1 : 0;
                    const int sp = P - Rp > iso ? 1 : P + Rp < iso ? -1 : 0;
                    if (sb == 0 || (side_box != 0 && sb != side_box)) drop_box = false;
                    if (sp == 0 || (side_pt != 0 && sp != side_pt)) drop_pt = false;
                    side_box = sb; side_pt = sp;
                    if (k == 62) { rel_w += Rb / fmax(fabs(P), 1e-30); nrel++; }
                }
                cnt[cls][drop_box ? 1 : 0]++;
                cntp[cls] += drop_pt;
                // inside bricks: the remainder by the shells' smallest exponent, against rho at the centre
                if (cls == 1 && getenv("GTO_ISO_TM_V") && K == 4 && DOT[0][0] < 6.0 * iso) {
                    static double part[6];
                    static long nb = 0;
                    static double sumP = 0, sumR = 0, sumPmin = 0;
                    double Pc = DOT[0][0];
                    for (int sh = 0; sh < d->num_ao_shells; ++sh) {
                        const iso_ao_shell_t* a = &d->ao_shells[sh];
                        double amin = 1e30; for (int q = 0; q < a->nprim; ++q) amin = fmin(amin, a->alpha[q]);
                        const int bin = amin < 0.2 ? 0 : amin < 0.5 ? 1 : amin < 1.0 ? 2 : amin < 3.0 ? 3 : amin < 10.0 ? 4 : 5;
                        for (int ci = 0; ci < a->ncart; ++ci) {
                            const int i = a->ao + ci;
                            double W = 0.0;
                            for (int aa = 0; aa < nm; ++aa) { const int o = t.mi[aa][0] + t.mi[aa][1] + t.mi[aa][2]; W += fabs(Dv[(size_t)i * 20 + aa]) * pow(hb, o) / t.mi_fact[aa]; }
                            part[bin] += 2.0 * E[i] * W / iso;
                        }
                    }
                    nb++; sumP += Pc / iso; sumR += (RW + REE) / iso;
                    if (nb % 100 == 0) printf("    inside bricks %ld: mean rho_c/iso %.2f, R/iso %.2f; R/iso by smallest exponent <0.2 %.2f <0.5 %.2f <1 %.2f <3 %.2f <10 %.2f >=10 %.2f\n",
                        nb, sumP / nb, sumR / nb, part[0] / nb, part[1] / nb, part[2] / nb, part[3] / nb, part[4] / nb, part[5] / nb);
                    (void)sumPmin;
                }
                if (cls == 0 && (drop_box || drop_pt)) bad++;
            }
            printf("  %-20s h %.2f iso %.4f K %d: bricks surface %ld, inside %ld (dropped %ld box, %ld per point), outside %ld (dropped %ld, %ld) | surface dropped %ld (must be 0), corners outside the enclosure %ld | mean R/|P| at the centre %.3f\n",
                files[fi], h, iso, K, cnt[0][0] + cnt[0][1], cnt[1][0] + cnt[1][1], cnt[1][1], cntp[1], cnt[2][0] + cnt[2][1], cnt[2][1], cntp[2], bad, outside_enc, nrel ? rel_w / nrel : 0.0);
            free(V); free(E); free(Dv);
        }
        md_gto_iso_gpu_set_destroy(g.s, set);
        free(dbg.bricks); free(dbg.slots); free(dbg.corners);
        free(D);
        md_gto_iso_field_destroy(field);
        iso_dens_prod_free(&p);
        iso_dens_src_free(&s);
    }
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// ---------------------------------------------------------------------------------------------
// EXPERIMENT (GTO_ISO_TM_MO=1): the same Taylor model for an orbital, psi = sum c_mu phi_mu: the
// polynomial sum c T exact, the remainder sum |c| E (one number per brick). A brick whose 125 points
// all lie in one band (above +iso, between, below -iso) has no surface cell of either sheet. Over the
// bricks the GPU's octree kept (they passed its bound), how many the model would drop besides.
// GTO_ISO_TM_FILES=a.h5,b.h5 picks the systems, GTO_ISO_TM_ISO the iso value (0.05).
// ---------------------------------------------------------------------------------------------

UTEST(gto_iso_gpu_lib, taylor_brick_mo_experiment) {
    if (!getenv("GTO_ISO_TM_MO")) UTEST_SKIP("set GTO_ISO_TM_MO=1 (EXPERIMENT)");
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    char flist[512];
    snprintf(flist, sizeof(flist), "%s", getenv("GTO_ISO_TM_FILES") ? getenv("GTO_ISO_TM_FILES") : "amide.h5,mol.h5,tq.scf.results.h5");
    const double h = getenv("GTO_ISO_TM_H") ? atof(getenv("GTO_ISO_TM_H")) : 0.2;
    const double iso = getenv("GTO_ISO_TM_ISO") ? atof(getenv("GTO_ISO_TM_ISO")) : 0.05;
    for (char* fname = strtok(flist, ","); fname; fname = strtok(NULL, ",")) {
        iso_dens_src_t s;
        if (!iso_dens_src_load(&s, fname)) { printf("  %s: not loaded\n", fname); continue; }
        for (int which = -1; which >= -2; --which) {
            iso_dens_prod_t p;
            md_gto_iso_field_t* field = iso_gl_field(&s, which, &p);
            ASSERT_TRUE(field != NULL);
            const size_t lumo = vlx_test_lumo_idx(&s.t, STR_LIT("orbital/alpha/occupation"));
            const double* cmo = s.C + (which == -1 ? lumo - 1 : lumo) * s.num_ao;
            // the orbital as a one-term product field: its AO shells and coefficients in the density's
            // convention (the Taylor coefficients below carry the normalizations)
            const double one = 1.0;
            md_gto_iso_field_t* aux = md_gto_iso_field_create_products(&s.basis, s.xyz, 0, cmo, NULL, &one, 1, md_get_heap_allocator());
            ASSERT_TRUE(aux != NULL);
            const iso_density_t* d = aux->f.dens;
            const int nao = d->num_ao;
            const md_gto_iso_field_t* fl[1] = {field};
            md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
            ASSERT_TRUE(set != NULL);
            md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 8, .relax_iterations = 3};
            gd.iso[0] = iso;
            isog_debug_t dbg = {0};
            isog_debug = &dbg;
            md_gto_iso_gpu_mesh_t gm;
            ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &gd, &gm));
            isog_debug = NULL;
            md_gto_iso_gpu_mesh_free(g.s, &gm);

            const bool hybrid = getenv("GTO_ISO_TM_HYBRID") != NULL;
            for (int K = 2; K <= 4; ++K) {
                tm_ctx_t t, tk[5];
                tm_init(&t, K);
                for (int k = 1; k <= K; ++k) tm_init(&tk[k], k);
                const int nm = t.nmi;
                double* V = (double*)malloc(sizeof(double) * (size_t)nao * 20);
                double* E = (double*)malloc(sizeof(double) * (size_t)nao);
                double* Vk = (double*)malloc(sizeof(double) * (size_t)nao * 20);
                double* Ek = (double*)malloc(sizeof(double) * (size_t)nao);
                long n_surf = 0, n_empty = 0, drop = 0, bad = 0, outside_enc = 0, reach = 0;
                long keep_hist[6] = {0};     // undropped empty bricks by max |psi| / iso over their corners
                double* rr = (double*)malloc(sizeof(double) * (dbg.num_bricks + 1));
                long nrr = 0;
                double rsum = 0.0;
                for (uint32_t b = 0; b < dbg.num_bricks; ++b) {
                    const uint32_t m = dbg.bricks[4 * b + 3];
                    if (!(m & 3u)) continue;
                    const float* cv = dbg.corners + 125 * (size_t)dbg.slots[b];
                    int band0 = 2, surf = 0;
                    for (int k = 0; k < 125; ++k) {
                        const int bd = cv[k] > iso ? 1 : cv[k] < -iso ? -1 : 0;
                        if (band0 == 2) band0 = bd; else if (bd != band0) surf = 1;
                    }
                    const double c[3] = {dbg.origin[0] + (dbg.bricks[4 * b] * 4.0 + 2.0) * h, dbg.origin[1] + (dbg.bricks[4 * b + 1] * 4.0 + 2.0) * h, dbg.origin[2] + (dbg.bricks[4 * b + 2] * 4.0 + 2.0) * h};
                    const double hb = 2.0 * h;
                    double T[20] = {0}, R = 0.0;
                    for (int sh = 0; sh < d->num_ao_shells; ++sh) {
                        const iso_ao_shell_t* a = &d->ao_shells[sh];
                        double dx = 0.0;     // the box's distance from the shell
                        for (int k = 0; k < 3; ++k) { const double e = fabs(c[k] - a->c[k]) - hb; if (e > 0.0) dx += e * e; }
                        if (dx > a->r2_cut) continue;
                        tm_shell(&t, a, c, hb, K, V + (size_t)a->ao * 20, E + a->ao);
                        if (hybrid) {
                            // per AO the order whose remainder is smallest (lower orders: the higher
                            // coefficients dropped, their terms inside the remainder)
                            for (int k = 1; k < K; ++k) {
                                tm_shell(&tk[k], a, c, hb, k, Vk + (size_t)a->ao * 20, Ek + a->ao);
                                for (int ci = 0; ci < a->ncart; ++ci) {
                                    const int i = a->ao + ci;
                                    if (Ek[i] < E[i]) {
                                        E[i] = Ek[i];
                                        for (int q = 0; q < nm; ++q) V[(size_t)i * 20 + q] = q < tk[k].nmi ? Vk[(size_t)i * 20 + q] : 0.0;
                                    }
                                }
                            }
                        }
                        for (int ci = 0; ci < a->ncart; ++ci) {
                            const double cc = d->C[(size_t)(a->ao + ci) * d->num_fn];
                            if (cc == 0.0) continue;
                            reach++;
                            for (int q = 0; q < nm; ++q) T[q] += cc * V[(size_t)(a->ao + ci) * 20 + q];
                            R += fabs(cc) * E[a->ao + ci];
                        }
                    }
                    rsum += R / iso;
                    bool dr = true;
                    int band = 2;
                    for (int k = 0; k < 125; ++k) {
                        const double dd[3] = {((k % 5) - 2) * h, (((k / 5) % 5) - 2) * h, ((k / 25) - 2) * h};
                        double P = 0.0;
                        for (int q = 0; q < nm; ++q) P += T[q] * pow(dd[0], t.mi[q][0]) * pow(dd[1], t.mi[q][1]) * pow(dd[2], t.mi[q][2]) / t.mi_fact[q];
                        if (cv[k] < P - R - 1e-5 * fabs(P) - 1e-9 || cv[k] > P + R + 1e-5 * fabs(P) + 1e-9) outside_enc++;
                        const int bd = P - R > iso ? 1 : (P + R < iso && P - R > -iso) ? 0 : P + R < -iso ? -1 : 2;
                        if (bd == 2 || (band != 2 && bd != band)) dr = false;
                        band = bd;
                    }
                    if (surf) { n_surf++; if (dr) bad++; }
                    else {
                        n_empty++;
                        if (dr) drop++;
                        else {
                            double mx = 0.0;
                            for (int k = 0; k < 125; ++k) mx = fmax(mx, fabs(cv[k]));
                            const double r = mx / iso;
                            keep_hist[r < 0.25 ? 0 : r < 0.5 ? 1 : r < 0.75 ? 2 : r < 1.0 ? 3 : r < 2.0 ? 4 : 5]++;
                        }
                    }
                    if (!surf) rr[nrr++] = R / iso;
                }
                // median R/iso over the empty bricks
                for (long i = 1; i < nrr; ++i) { double v = rr[i]; long j = i; while (j > 0 && rr[j - 1] > v) { rr[j] = rr[j - 1]; --j; } rr[j] = v; }
                const long nb = n_surf + n_empty;
                printf("  %-20s %s h %.2f iso %.3f K %d%s: bricks %ld, surface %ld, empty %ld, dropped %ld (%.0f%% of empty) | surface dropped %ld (must be 0), outside the enclosure %ld | AOs/brick %.0f, R/iso median %.3f\n",
                    fname, which == -1 ? "HOMO" : "LUMO", h, iso, K, hybrid ? " hybrid" : "", nb, n_surf, n_empty, drop, n_empty ? 100.0 * drop / n_empty : 0.0, bad, outside_enc, nb ? (double)reach / nb : 0.0, nrr ? rr[nrr / 2] : 0.0);
                printf("      kept empty bricks by max|psi|/iso: <0.25 %ld, <0.5 %ld, <0.75 %ld, <1 %ld, <2 %ld, >=2 %ld\n", keep_hist[0], keep_hist[1], keep_hist[2], keep_hist[3], keep_hist[4], keep_hist[5]);
                (void)rsum;
                free(V); free(E); free(Vk); free(Ek); free(rr);
            }
            md_gto_iso_gpu_set_destroy(g.s, set);
            free(dbg.bricks); free(dbg.slots); free(dbg.corners);
            md_gto_iso_field_destroy(aux);
            md_gto_iso_field_destroy(field);
            iso_dens_prod_free(&p);
        }
        iso_dens_src_free(&s);
    }
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// ---------------------------------------------------------------------------------------------
// TEST: the flood against the exact corners, over iso values up to where the
// pieces are atoms' small blobs: the meshes' vertex and triangle counts. The default (flood 1: SCF
// densities, marked as having their maxima at the atoms) must match exactly; the flood forced on the
// attachment and detachment densities (flood 2), and without the atoms' bricks as seeds, for what they
// miss. Amide by default; GTO_ISO_FLOOD_SWEEP_FILES=h2o.h5,amide.h5,mol.h5 for more systems.
// ---------------------------------------------------------------------------------------------

UTEST(gto_iso_gpu_lib, flood_sweep) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    md_gto_iso_gpu_switches_t sw0;
    md_gto_iso_gpu_switches(&sw0);
    char flist[512];
    snprintf(flist, sizeof(flist), "%s", getenv("GTO_ISO_FLOOD_SWEEP_FILES") ? getenv("GTO_ISO_FLOOD_SWEEP_FILES") : "amide.h5");
    const double isos[] = {0.0025, 0.02, 0.1, 0.3, 1.0, 3.0};
    const int kinds[] = {ISO_D_SCF, ISO_D_ATTACH, ISO_D_DETACH};
    const char* kname[] = {"scf", "attach", "detach"};
    long worse_atoms = 0, worse_forced = 0, worse_plain = 0;
    for (char* fname = strtok(flist, ","); fname; fname = strtok(NULL, ",")) {
        iso_dens_src_t s;
        if (!iso_dens_src_load(&s, fname)) { printf("  %s: not loaded\n", fname); continue; }
        for (int ki = 0; ki < 3; ++ki) {
            if (kinds[ki] != ISO_D_SCF && s.num_nto == 0) continue;
            iso_dens_prod_t p;
            md_gto_iso_field_t* field = iso_gl_field(&s, kinds[ki], &p);
            ASSERT_TRUE(field != NULL);
            md_gto_iso_field_set_maxima_at_atoms(field, kinds[ki] == ISO_D_SCF);
            const md_gto_iso_field_t* fl[1] = {field};
            md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
            ASSERT_TRUE(set != NULL);
            for (int ii = 0; ii < (int)(sizeof(isos) / sizeof(isos[0])); ++ii) {
                uint32_t nv[4] = {0}, nt[4] = {0};
                for (int mode = 0; mode < 4; ++mode) {      // exact, the default, flood forced, forced without the atoms
                    md_gto_iso_gpu_switches_t sw = sw0;
                    sw.flood = mode == 0 ? 0 : mode == 1 ? 1 : 2;
                    md_gto_iso_gpu_switches_set(&sw);
                    isog_flood_atoms = mode == 3 ? 0 : 1;
                    md_gto_iso_gpu_desc_t gd = {.cell_size = 0.2, .project_iterations = 8, .relax_iterations = 3};
                    gd.iso[0] = isos[ii];
                    md_gto_iso_gpu_mesh_t gm;
                    ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &gd, &gm));
                    nv[mode] = gm.vertex_count;
                    nt[mode] = gm.index_count / 3;
                    md_gto_iso_gpu_mesh_free(g.s, &gm);
                }
                isog_flood_atoms = 1;
                const bool ok1 = nv[1] == nv[0] && nt[1] == nt[0], ok2 = nv[2] == nv[0] && nt[2] == nt[0], ok3 = nv[3] == nv[0] && nt[3] == nt[0];
                worse_atoms += !ok1;
                worse_forced += !ok2;
                worse_plain += !ok3;
                printf("  %-14s %-6s iso %-7g exact %6u v %6u t | default %6u v%s | forced %6u v%s | forced, no atom seeds %6u v%s\n", fname, kname[ki], isos[ii],
                       nv[0], nt[0], nv[1], ok1 ? "" : " DIFFERS", nv[2], ok2 ? "" : " MISSED", nv[3], ok3 ? "" : " MISSED");
            }
            md_gto_iso_gpu_set_destroy(g.s, set);
            md_gto_iso_field_destroy(field);
            iso_dens_prod_free(&p);
        }
        iso_dens_src_free(&s);
    }
    md_gto_iso_gpu_switches_set(&sw0);
    printf("  cases where the default differs from the exact mesh: %ld; the flood forced: %ld, and without the atom seeds: %ld\n", worse_atoms, worse_forced, worse_plain);
    EXPECT_EQ(worse_atoms, 0);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// ---------------------------------------------------------------------------------------------
// EXPERIMENT (GTO_ISO_FLOOD=1): corners from seeds and a flood instead of every candidate brick. The
// candidates (bricks the octree keeps) are sampled on a coarse lattice (stride 4: the brick corners;
// stride 2: 27 points a brick); a candidate whose samples straddle the iso value is a seed and is
// evaluated in full; then, round by round, a candidate across a face of an evaluated brick whose 25
// face points hold a sign change on a face edge is evaluated too, until none is added. Surface bricks
// never evaluated are missed (their cells would be missing from the mesh). Simulated on the GPU's own
// candidates and corner values, per sheet.
// ---------------------------------------------------------------------------------------------

typedef struct flood_stats_t {
    long cand, surface, seeds, evaluated, missed, missed_cells, rounds, coarse_points, eval_points, cand_points;
} flood_stats_t;

// The bricks' face neighbour index (or -1) through a hash of the coordinates
static void flood_sim(const isog_debug_t* dbg, double iso, int sheet, int stride, flood_stats_t* st) {
    memset(st, 0, sizeof(*st));
    const uint32_t nb = dbg->num_bricks;
    iso_map_t map;
    iso_map_init(&map, nb + 1);
    int32_t* idx = (int32_t*)malloc(sizeof(int32_t) * (nb + 1));
    uint32_t nc = 0;
    for (uint32_t b = 0; b < nb; ++b) {
        idx[b] = -1;
        if (!((dbg->bricks[4 * b + 3] >> sheet) & 1u)) continue;
        iso_map_put(&map, iso_key3((int)dbg->bricks[4 * b], (int)dbg->bricks[4 * b + 1], (int)dbg->bricks[4 * b + 2]), b);
        idx[b] = (int32_t)nc++;
    }
    const double sg = sheet == 0 ? 1.0 : -1.0;
#define FL_V(b, x, y, z) (sg * dbg->corners[125 * (size_t)dbg->slots[b] + (x) + 5 * (y) + 25 * (z)] - iso > 0.0)
    uint8_t* state = (uint8_t*)calloc(nb + 1, 1);     // 1 evaluated
    uint32_t* queue = (uint32_t*)malloc(sizeof(uint32_t) * (nb + 1));
    uint32_t qn = 0;
    for (uint32_t b = 0; b < nb; ++b) {
        if (idx[b] < 0) continue;
        st->cand++;
        // a surface brick: a cell with a sign change
        long cells = 0;
        for (int c = 0; c < 64; ++c) {
            const int x = c & 3, y = (c >> 2) & 3, z = c >> 4;
            int pos = 0;
            for (int k = 0; k < 8; ++k) pos += FL_V(b, x + (k & 1), y + ((k >> 1) & 1), z + ((k >> 2) & 1));
            cells += pos != 0 && pos != 8;
        }
        if (cells) st->surface++;
        // seed: the coarse samples straddle
        int pos = 0, n = 0;
        for (int z = 0; z <= 4; z += stride) for (int y = 0; y <= 4; y += stride) for (int x = 0; x <= 4; x += stride) { pos += FL_V(b, x, y, z); n++; }
        if (pos != 0 && pos != n) { state[b] = 1; queue[qn++] = b; st->seeds++; }
    }
    // flood by rounds
    uint32_t head = 0;
    while (head < qn) {
        const uint32_t end = qn;
        st->rounds++;
        for (; head < end; ++head) {
            const uint32_t b = queue[head];
            const int bx = (int)dbg->bricks[4 * b], by = (int)dbg->bricks[4 * b + 1], bz = (int)dbg->bricks[4 * b + 2];
            for (int f = 0; f < 6; ++f) {
                const int ax = f >> 1, side = f & 1;      // the face at coordinate 0 or 4 along axis ax
                const int u = (ax + 1) % 3, v = (ax + 2) % 3;
                bool change = false;
                for (int i = 0; i <= 4 && !change; ++i) for (int j = 0; j <= 4 && !change; ++j) {
                    int p[3]; p[ax] = side ? 4 : 0; p[u] = i; p[v] = j;
                    const bool s0 = FL_V(b, p[0], p[1], p[2]);
                    if (i < 4) { int q[3] = {p[0], p[1], p[2]}; q[u]++; if (FL_V(b, q[0], q[1], q[2]) != s0) change = true; }
                    if (j < 4) { int q[3] = {p[0], p[1], p[2]}; q[v]++; if (FL_V(b, q[0], q[1], q[2]) != s0) change = true; }
                }
                if (!change) continue;
                int nbc[3] = {bx, by, bz};
                nbc[ax] += side ? 1 : -1;
                uint32_t nbi;
                if (nbc[0] < 0 || nbc[1] < 0 || nbc[2] < 0) continue;
                if (!iso_map_get(&map, iso_key3(nbc[0], nbc[1], nbc[2]), &nbi) || idx[nbi] < 0 || state[nbi]) continue;
                state[nbi] = 1;
                queue[qn++] = nbi;
            }
        }
    }
    st->evaluated = qn;
    // lattice points: every candidate's (today) against the evaluated bricks' plus the coarse samples
    {
        iso_map_t pm;
        iso_map_init(&pm, 64 * (size_t)(nc + 1));
        long all = 0, ev = 0, co = 0;
        for (int pass = 0; pass < 3; ++pass) {
            iso_map_free(&pm);
            iso_map_init(&pm, 64 * (size_t)(nc + 1));
            long n = 0;
            for (uint32_t b = 0; b < nb; ++b) {
                if (idx[b] < 0 || (pass == 1 && !state[b])) continue;
                const int st_ = pass == 2 ? stride : 1;
                for (int z = 0; z <= 4; z += st_) for (int y = 0; y <= 4; y += st_) for (int x = 0; x <= 4; x += st_) {
                    const int gx = (int)dbg->bricks[4 * b] * 4 + x, gy = (int)dbg->bricks[4 * b + 1] * 4 + y, gz = (int)dbg->bricks[4 * b + 2] * 4 + z;
                    uint32_t v;
                    const uint64_t key = iso_key3(gx, gy, gz);
                    if (!iso_map_get(&pm, key, &v)) { iso_map_put(&pm, key, 1); n++; }
                }
            }
            if (pass == 0) all = n; else if (pass == 1) ev = n; else co = n;
        }
        iso_map_free(&pm);
        st->coarse_points = co;
        st->eval_points = ev;
        st->cand_points = all;
    }
    for (uint32_t b = 0; b < nb; ++b) {
        if (idx[b] < 0 || state[b]) continue;
        long cells = 0;
        for (int c = 0; c < 64; ++c) {
            const int x = c & 3, y = (c >> 2) & 3, z = c >> 4;
            int pos = 0;
            for (int k = 0; k < 8; ++k) pos += FL_V(b, x + (k & 1), y + ((k >> 1) & 1), z + ((k >> 2) & 1));
            cells += pos != 0 && pos != 8;
        }
        if (cells) { st->missed++; st->missed_cells += cells; }
    }
#undef FL_V
    free(state); free(queue); free(idx);
    iso_map_free(&map);
}

UTEST(gto_iso_gpu_lib, flood_experiment) {
    if (!getenv("GTO_ISO_FLOOD_SURVEY")) UTEST_SKIP("set GTO_ISO_FLOOD_SURVEY=1 (what a flood would miss, simulated on the host)");
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    typedef struct { const char* file; bool nto; } fl_file_t;
    static const fl_file_t files[] = {{"amide.h5", true}, {"mol.h5", true}, {"h2o.h5", true}, {"tq.scf.results.h5", false},
                                      {"acro-xps.h5", false}, {"meth.scf.results.h5", false}, {"myjob.scf.h5", false}};
    typedef struct { int kind; const char* name; double iso; } fl_case_t;
    static const fl_case_t cases[] = {
        {ISO_D_SCF, "density", 0.0025}, {ISO_D_SCF, "density", 0.02}, {ISO_D_SCF, "density", 0.1},
        {-1, "HOMO", 0.05}, {-1, "HOMO", 0.02}, {-2, "LUMO", 0.05}, {-2, "LUMO", 0.02},
        {ISO_D_DIFF, "difference", 0.004}, {ISO_D_DIFF, "difference", 0.001}, {ISO_D_TRANS, "transition", 0.004}, {ISO_D_TRANS, "transition", 0.001},
    };
    const char* hs = getenv("GTO_ISO_FLOOD_H");
    double hl[8]; int nh = iso_parse_list(hs ? hs : "0.2,0.12", hl, 8);
    long tot_missed[2] = {0}, tot_cases = 0;
    for (size_t fi = 0; fi < sizeof(files) / sizeof(files[0]); ++fi) {
        iso_dens_src_t s;
        if (!iso_dens_src_load_ex(&s, files[fi].file, files[fi].nto)) { printf("  %s: not loaded\n", files[fi].file); continue; }
        printf("  %s: %zu atoms, %zu AOs\n", files[fi].file, s.num_atoms, s.num_ao);
        for (size_t ci = 0; ci < sizeof(cases) / sizeof(cases[0]); ++ci) {
            const fl_case_t* fc = &cases[ci];
            if ((fc->kind == ISO_D_DIFF || fc->kind == ISO_D_TRANS) && !files[fi].nto) continue;
            iso_dens_prod_t p;
            md_gto_iso_field_t* field = NULL;
            if (fc->kind == -2) {
                memset(&p, 0, sizeof(p));
                const size_t lumo = vlx_test_lumo_idx((vlx_test_t*)&s.t, STR_LIT("orbital/alpha/occupation"));
                field = md_gto_iso_field_create(&s.basis, s.xyz, 0, s.C + lumo * s.num_ao, md_get_heap_allocator());
            } else {
                field = iso_gl_field(&s, fc->kind, &p);
            }
            if (!field) continue;
            const md_gto_iso_field_t* fl[1] = {field};
            md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
            for (int hi = 0; hi < nh && set; ++hi) {
                md_gto_iso_gpu_desc_t gd = {.cell_size = hl[hi], .project_iterations = 8, .relax_iterations = 3};
                gd.iso[0] = fc->iso;
                isog_debug_t dbg = {0};
                isog_debug = &dbg;
                md_gto_iso_gpu_mesh_t gm;
                const bool ok = md_gto_iso_gpu_extract(g.s, set, &gd, &gm);
                isog_debug = NULL;
                if (ok) md_gto_iso_gpu_mesh_free(g.s, &gm);
                for (int sheet = 0; sheet < 2 && ok; ++sheet) {
                    for (int stride = 4; stride >= 2; stride -= 2) {
                        flood_stats_t st;
                        flood_sim(&dbg, fc->iso, sheet, stride, &st);
                        if (st.surface == 0 && st.cand == 0) continue;
                        printf("    %-10s iso %c%-6g h %.2f stride %d: candidates %5ld, surface %5ld | seeds %5ld, evaluated %5ld (%4.1f%%) in %ld rounds | points %4.1f%% (coarse %4.1f%%) | MISSED %ld bricks, %ld cells\n",
                               fc->name, sheet ? '-' : '+', fc->iso, hl[hi], stride, st.cand, st.surface, st.seeds, st.evaluated,
                               st.cand ? 100.0 * st.evaluated / st.cand : 0.0, st.rounds,
                               st.cand_points ? 100.0 * (st.eval_points + st.coarse_points) / st.cand_points : 0.0,
                               st.cand_points ? 100.0 * st.coarse_points / st.cand_points : 0.0, st.missed, st.missed_cells);
                        tot_missed[stride == 2] += st.missed;
                    }
                }
                tot_cases++;
                free(dbg.bricks); free(dbg.slots); free(dbg.corners);
            }
            if (set) md_gto_iso_gpu_set_destroy(g.s, set);
            md_gto_iso_field_destroy(field);
            iso_dens_prod_free(&p);
        }
        iso_dens_src_free(&s);
    }
    printf("  missed surface bricks over all cases: stride 4 %ld, stride 2 %ld (%ld extractions)\n", tot_missed[0], tot_missed[1], tot_cases);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}


// A large system's total density, as viamd extracts its mesh: the density matrix from the published virtual
// attribute, the field from it, the GPU set and the extraction, each timed. GTO_ISO_LARGE=<file in test_data/vlx>
// (a protein: 1152 atoms, 9072 AOs); GTO_ISO_LARGE_H (default 0.4), GTO_ISO_LARGE_ISO (default 0.0025).
#include <md_qm.h>
// md_gto_iso_gpu_evaluate: the fields of a set at arbitrary points against md_gto_iso_field_evaluate in double,
// through every vertex evaluator (orbitals, the tiled and the wide GEMM, the D form), from the host and on the
// device at a mesh's own vertices
typedef struct iso_ev_err_t { double val, grad, fmax, gmax; } iso_ev_err_t;

static iso_ev_err_t iso_ev_compare(const md_gto_iso_field_t* f, const float* xyz, size_t stride, size_t n, const float* val, const float* grad, size_t gstride) {
    double* rv = (double*)malloc(sizeof(double) * n);
    double* rg = (double*)malloc(sizeof(double) * 3 * n);
    md_gto_iso_field_evaluate(f, xyz, n, stride, rv, rg);
    iso_ev_err_t e = {0};
    for (size_t i = 0; i < n; ++i) {
        e.fmax = fmax(e.fmax, fabs(rv[i]));
        e.val = fmax(e.val, val[i] == val[i] ? fabs((double)val[i] - rv[i]) : INFINITY);
        for (int c = 0; c < 3; ++c) {
            e.gmax = fmax(e.gmax, fabs(rg[3 * i + c]));
            if (grad) { const float gc = grad[gstride * i + c]; e.grad = fmax(e.grad, gc == gc ? fabs((double)gc - rg[3 * i + c]) : INFINITY); }
        }
    }
    free(rv); free(rg);
    return e;
}

UTEST(gto_iso_gpu_lib, evaluate_matches_cpu) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    printf("  revision %s\n", ISOG_REVISION);
    static const char* files[] = {"amide.h5", "mol.h5"};
    for (int fi = 0; fi < 2; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load_ex(&s, files[fi], fi == 0));
        const size_t lumo = vlx_test_lumo_idx(&s.t, STR_LIT("orbital/alpha/occupation"));
        // 0, 1: HOMO, LUMO; 2: SCF density (products); 3: attachment - detachment of state 1 (signed; with NTOs)
        enum { NF = 4 };
        md_gto_iso_field_t* fields[NF];
        iso_dens_prod_t p[NF];
        memset(p, 0, sizeof(p));
        fields[0] = md_gto_iso_field_create(&s.basis, s.xyz, 0, s.C + (lumo - 1) * s.num_ao, md_get_heap_allocator());
        fields[1] = md_gto_iso_field_create(&s.basis, s.xyz, 0, s.C + lumo * s.num_ao, md_get_heap_allocator());
        iso_dens_products(&p[2], &s, ISO_D_SCF);
        fields[2] = iso_dens_field(&s, &p[2], false);
        iso_dens_products(&p[3], &s, ISO_D_DIFF);
        fields[3] = p[3].count ? iso_dens_field(&s, &p[3], false) : NULL;
        for (int k = 0; k < 3; ++k) ASSERT_TRUE(fields[k] != NULL);

        // random points, a third near the atoms, in no particular order
        enum { N = 4000 };
        iso_rng_t rng = {0x5eed + (uint64_t)fi};
        double* pd = (double*)malloc(sizeof(double) * 3 * N);
        float* pts = (float*)malloc(sizeof(float) * 3 * N);
        iso_dens_points(pd, N, &s, &rng);
        for (int i = 0; i < 3 * N; ++i) pts[i] = (float)pd[i];
        float* val = (float*)malloc(sizeof(float) * (NF + 1) * N);        // NF fields and a repeat
        float* grad = (float*)malloc(sizeof(float) * 3 * (NF + 1) * N);

        typedef struct { const char* name; int fields[NF]; int nf; int dform; } ev_case_t;
        const ev_case_t cases[] = {
            {"orbitals",          {0, 1},       2, 0},
            {"scf density",       {2},          1, 0},
            {"mixed",             {0, 2, 3, 1}, 4, 0},
            {"densities, D form", {2, 3},       2, 1},
        };
        for (int ci = 0; ci < (int)(sizeof(cases) / sizeof(cases[0])); ++ci) {
            const ev_case_t* c = &cases[ci];
            const md_gto_iso_field_t* fl[NF];
            int nf = 0;
            for (int k = 0; k < c->nf; ++k) if (fields[c->fields[k]]) fl[nf++] = fields[c->fields[k]];
            const int saved_dform = isog_dform;
            isog_dform = c->dform;
            md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, (uint32_t)nf, md_get_heap_allocator());
            ASSERT_TRUE(set != NULL);
            isog_args_t a;
            memset(&a, 0, sizeof(a));
            isog_args_set(&a, set);
            const isog_vkernel_t vk = isog_vertex_kernel(set, a.dform != 0u, false);
            // every field of the set, the last one twice (repeats are allowed)
            uint32_t which[NF + 1];
            for (int k = 0; k < nf; ++k) which[k] = (uint32_t)k;
            which[nf] = (uint32_t)(nf - 1);
            md_tick_t t0 = md_tick_now();
            const bool ok = md_gto_iso_gpu_evaluate(g.s, set, which, (uint32_t)nf + 1, pts, 0, N, val, grad);
            const double secs = md_tick_to_seconds(md_tick_now() - t0);
            EXPECT_TRUE(ok);
            if (ok) {
                for (int k = 0; k <= nf; ++k) {
                    const md_gto_iso_field_t* f = fl[which[k]];
                    const iso_ev_err_t e = iso_ev_compare(f, pts, 0, N, val + (size_t)k * N, grad + 3 * (size_t)k * N, 3);
                    printf("  %-9s %-18s %-24s field %d: max |f| %.2e, |df| %.1e (%.1e rel), max |grad| %.2e, |dgrad| %.1e (%.1e rel) | %.3f s\n",
                           files[fi], c->name, isog_kname[vk.k] ? isog_kname[vk.k] : "?", k, e.fmax, e.val, e.val / e.fmax, e.gmax, e.grad, e.grad / e.gmax, secs);
                    // float, and the AOs below the evaluation cutoff (1e-6) left out: relative and absolute parts
                    EXPECT_LT(e.val, 1.0e-5 * e.fmax + 1.0e-5);
                    EXPECT_LT(e.grad, 1.0e-5 * e.gmax + 1.0e-5);
                }
                // the repeated field: the same values
                EXPECT_EQ(0, memcmp(val + (size_t)(nf - 1) * N, val + (size_t)nf * N, sizeof(float) * N));
            }
            md_gto_iso_gpu_set_destroy(g.s, set);
            isog_dform = saved_dform;
        }

        // on the device: the SCF density at the vertices of the HOMO's mesh (stride 16), values only, then
        // values and gradients
        {
            const md_gto_iso_field_t* fl[2] = {fields[0], fields[2]};
            md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 2, md_get_heap_allocator());
            ASSERT_TRUE(set != NULL);
            md_gto_iso_gpu_desc_t gd = {.cell_size = 0.25, .project_iterations = 4, .relax_iterations = 2};
            gd.iso[0] = 0.05;
            gd.iso[1] = 0.02;
            md_gto_iso_gpu_mesh_t gm;
            ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, set, &gd, &gm));
            const uint32_t v0 = gm.vertex_offset[0][0];
            const uint32_t nv = gm.vertex_count_of[0][0] + gm.vertex_count_of[0][1];   // the HOMO's two sheets
            ASSERT_GT(nv, 0u);
            md_gpu_mem_t d = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, (sizeof(float) + sizeof(md_gpu_float4)) * nv + 16);
            md_gpu_mem_t h = md_gpu_malloc(g.s, MD_GPU_MEM_HOST_READ, (sizeof(float) + sizeof(md_gpu_float4)) * nv + 16);
            ASSERT_TRUE(d.gpu && h.gpu);
            md_gto_iso_gpu_eval_desc_t ed = {0};
            ed.xyz = gm.position + (md_gpu_addr_t)v0 * sizeof(md_gpu_float4);
            ed.xyz_stride = sizeof(md_gpu_float4);
            ed.num_xyz = nv;
            ed.num_fields = 1;
            ed.fields[0] = 1;
            ed.out_value = d.gpu;
            ed.out_failed = d.gpu + (sizeof(float) + sizeof(md_gpu_float4)) * nv;
            ASSERT_TRUE(md_gto_iso_gpu_evaluate_launch(g.s, set, &ed));
            md_gpu_copy(g.s, h.gpu, d.gpu, sizeof(float) * nv);
            md_gpu_copy(g.s, h.gpu + (sizeof(float) + sizeof(md_gpu_float4)) * nv, ed.out_failed, sizeof(uint32_t));
            md_gpu_stream_sync(g.s);
            float* only = (float*)malloc(sizeof(float) * nv);
            memcpy(only, h.cpu, sizeof(float) * nv);
            EXPECT_EQ(0u, *(const uint32_t*)((const char*)h.cpu + (sizeof(float) + sizeof(md_gpu_float4)) * nv));
            ed.out_gradient = d.gpu + sizeof(float) * nv;
            ASSERT_TRUE(md_gto_iso_gpu_evaluate_launch(g.s, set, &ed));
            md_gpu_copy(g.s, h.gpu, d.gpu, (sizeof(float) + sizeof(md_gpu_float4)) * nv);
            md_gpu_stream_sync(g.s);
            const float* hv = (const float*)h.cpu;
            const float* hg = (const float*)((const char*)h.cpu + sizeof(float) * nv);
            EXPECT_EQ(0, memcmp(only, hv, sizeof(float) * nv));
            for (uint32_t i = 0; i < nv; ++i) EXPECT_EQ(hv[i], hg[4 * i + 3]);
            float* pos = (float*)malloc(sizeof(float) * 4 * gm.vertex_count);
            ASSERT_TRUE(md_gto_iso_gpu_mesh_read_all(g.s, &gm, pos, NULL, NULL));
            const iso_ev_err_t e = iso_ev_compare(fields[2], pos + 4 * (size_t)v0, 4 * sizeof(float), nv, hv, hg, 4);
            printf("  %-9s scf density at %u HOMO vertices (device): max |f| %.2e, |df| %.1e (%.1e rel), |dgrad| %.1e (%.1e rel)\n",
                   files[fi], nv, e.fmax, e.val, e.val / e.fmax, e.grad, e.grad / e.gmax);
            EXPECT_LT(e.val, 1.0e-5 * e.fmax + 1.0e-5);
            EXPECT_LT(e.grad, 1.0e-5 * e.gmax + 1.0e-5);
            free(pos); free(only);
            md_gpu_free(g.s, d.gpu);
            md_gpu_free(g.s, h.gpu);
            md_gto_iso_gpu_mesh_free(g.s, &gm);
            md_gto_iso_gpu_set_destroy(g.s, set);
        }

        free(val); free(grad); free(pts); free(pd);
        for (int k = 0; k < NF; ++k) { if (fields[k]) md_gto_iso_field_destroy(fields[k]); iso_dens_prod_free(&p[k]); }
        iso_dens_src_free(&s);
    }
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

// The average local ionization energy at the vertices of the 0.001 isodensity surface, on the device: the
// density and the energy weighted density as two fields of one set, evaluated where the extraction left the
// vertices, against the CPU
UTEST(gto_iso_gpu_lib, alie_at_mesh_vertices) {
    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    printf("  revision %s\n", ISOG_REVISION);
    static const char* files[] = {"h2o.h5", "amide.h5", "mol.h5"};
    for (int fi = 0; fi < 3; ++fi) {
        iso_dens_src_t s;
        ASSERT_TRUE(iso_dens_src_load_ex(&s, files[fi], false));
        iso_alie_t a;
        ASSERT_TRUE(iso_alie_make(&a, &s));
        md_gto_iso_field_set_maxima_at_atoms(a.rho, true);
        const md_gto_iso_field_t* fl[2] = {a.rho, a.num};
        md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 2, md_get_heap_allocator());
        ASSERT_TRUE(set != NULL);
        // the surface of field 0 only: one sheet of the density
        md_gto_iso_gpu_set_t* surf = md_gto_iso_gpu_set_create(g.s, &s.basis, s.xyz, 0, fl, 1, md_get_heap_allocator());
        ASSERT_TRUE(surf != NULL);
        md_gto_iso_gpu_desc_t gd = {.cell_size = 0.3, .project_iterations = 8, .relax_iterations = 2};
        gd.iso[0] = 0.001;
        gd.sheets[0] = 1;
        md_gto_iso_gpu_mesh_t gm;
        ASSERT_TRUE(md_gto_iso_gpu_extract(g.s, surf, &gd, &gm));
        const uint32_t nv = gm.vertex_count;
        ASSERT_GT(nv, 0u);
        md_gpu_mem_t d = md_gpu_malloc(g.s, MD_GPU_MEM_DEVICE, sizeof(float) * 2 * nv + 16);
        md_gpu_mem_t h = md_gpu_malloc(g.s, MD_GPU_MEM_HOST_READ, sizeof(float) * 2 * nv + 16);
        ASSERT_TRUE(d.gpu && h.gpu);
        md_gto_iso_gpu_eval_desc_t ed = {0};
        ed.xyz = gm.position;
        ed.xyz_stride = sizeof(md_gpu_float4);
        ed.num_xyz = nv;
        ed.num_fields = 2;
        ed.fields[0] = 0;
        ed.fields[1] = 1;
        ed.out_value = d.gpu;
        ed.out_failed = d.gpu + sizeof(float) * 2 * nv;
        md_tick_t t0 = md_tick_now();
        ASSERT_TRUE(md_gto_iso_gpu_evaluate_launch(g.s, set, &ed));
        md_gpu_copy(g.s, h.gpu, d.gpu, sizeof(float) * 2 * nv + sizeof(uint32_t));
        md_gpu_stream_sync(g.s);
        const double secs = md_tick_to_seconds(md_tick_now() - t0);
        const float* hv = (const float*)h.cpu;
        EXPECT_EQ(0u, *(const uint32_t*)(hv + 2 * nv));
        float* pos = (float*)malloc(sizeof(float) * 4 * nv);
        ASSERT_TRUE(md_gto_iso_gpu_mesh_read_all(g.s, &gm, pos, NULL, NULL));
        double* rv = (double*)malloc(sizeof(double) * nv);
        double* qv = (double*)malloc(sizeof(double) * nv);
        md_gto_iso_field_evaluate(a.rho, pos, nv, 4 * sizeof(float), rv, NULL);
        md_gto_iso_field_evaluate(a.num, pos, nv, 4 * sizeof(float), qv, NULL);
        double err = 0.0, smin = DBL_MAX, smax = 0.0;
        for (uint32_t i = 0; i < nv; ++i) {
            const double Ic = qv[i] / rv[i];
            const double Ig = (double)hv[nv + i] / (double)hv[i];
            err = fmax(err, Ig == Ig ? fabs(Ig - Ic) : INFINITY);
            smin = fmin(smin, Ic);
            smax = fmax(smax, Ic);
        }
        printf("  %-9s %zu occupied: I(r) at %u vertices of rho = 0.001: %.4f .. %.4f Eh, GPU against CPU %.1e Eh | %.3f s\n",
               files[fi], a.count, nv, smin, smax, err, secs);
        EXPECT_LT(err, 1.0e-4);
        free(pos); free(rv); free(qv);
        md_gpu_free(g.s, d.gpu);
        md_gpu_free(g.s, h.gpu);
        md_gto_iso_gpu_mesh_free(g.s, &gm);
        md_gto_iso_gpu_set_destroy(g.s, surf);
        md_gto_iso_gpu_set_destroy(g.s, set);
        iso_alie_free(&a);
        iso_dens_src_free(&s);
    }
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
}

UTEST(gto_iso_gpu_lib, large_density) {
    const char* file = getenv("GTO_ISO_LARGE");
    if (!file || !file[0]) UTEST_SKIP("set GTO_ISO_LARGE=<file in test_data/vlx>");
    const double h = getenv("GTO_ISO_LARGE_H") ? atof(getenv("GTO_ISO_LARGE_H")) : 0.4;
    const double iso = getenv("GTO_ISO_LARGE_ISO") ? atof(getenv("GTO_ISO_LARGE_ISO")) : 0.0025;
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    md_tick_t t0 = md_tick_now();
    vlx_test_t t;
    ASSERT_TRUE(vlx_test_load(&t, str_from_cstr(path), MEGABYTES(1024)));
    const size_t num_atoms = t.sys.atom.count;
    float* xyz = (float*)malloc(sizeof(float) * 3 * num_atoms);
    md_gto_basis_t basis;
    ASSERT_EQ(vlx_test_atom_xyz_bohr(xyz, 3 * num_atoms, &t), num_atoms);
    ASSERT_TRUE(qm_test_basis(&basis, &t));
    const size_t nao = md_gto_basis_num_ao(&basis);
    printf("  %s: %zu atoms, %zu AOs; loaded in %.2f s\n", file, num_atoms, nao, md_tick_to_seconds(md_tick_now() - t0));

    // what viamd does: the field from the occupied orbitals (md_qm_density_orbitals), or with
    // GTO_ISO_LARGE_MATRIX=1 from the matrix (factorized again)
    const bool from_matrix = getenv("GTO_ISO_LARGE_MATRIX") && atoi(getenv("GTO_ISO_LARGE_MATRIX")) != 0;
    t0 = md_tick_now();
    const md_attribute_t* da = qm_test_attr(&t, STR_LIT("orbital/total/density"));
    ASSERT_TRUE(da != NULL);
    double* D = (double*)malloc(sizeof(double) * nao * nao);
    ASSERT_EQ(md_qm_extract_symmetric_f64(D, nao * nao, da, md_attribute_slice_all()), nao);
    double tr = 0.0;
    for (size_t i = 0; i < nao; ++i) tr += D[i * nao + i];
    printf("  density matrix (orbital/total/density): %.2f s (trace %.6g)\n", md_tick_to_seconds(md_tick_now() - t0), tr);

    md_gto_iso_field_t* field = NULL;
    t0 = md_tick_now();
    if (from_matrix) {
        field = md_gto_iso_field_create_density(&basis, xyz, 3 * sizeof(float), D, md_get_heap_allocator());
    } else {
        double *rows = NULL, *w = NULL;
        size_t na = 0;
        const size_t count = md_qm_density_orbitals(&rows, &w, &na, &t.sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
        ASSERT_TRUE(count > 0);
        ASSERT_EQ(na, nao);
        printf("  orbitals: %zu, %.2f s\n", count, md_tick_to_seconds(md_tick_now() - t0));
        field = md_gto_iso_field_create_products(&basis, xyz, 3 * sizeof(float), rows, NULL, w, count, md_get_heap_allocator());
        md_free(md_get_heap_allocator(), rows, sizeof(double) * count * na);
        md_free(md_get_heap_allocator(), w, sizeof(double) * count);
    }
    ASSERT_TRUE(field != NULL);
    md_gto_iso_field_set_maxima_at_atoms(field, true);
    md_gto_iso_field_info_t fi;
    md_gto_iso_field_info(&fi, field);
    printf("  field (%s): %.2f s, %u functions, split %d, residual %g\n", from_matrix ? "matrix" : "orbitals", md_tick_to_seconds(md_tick_now() - t0), fi.num_functions, (int)fi.split, fi.residual);
    free(D);

    iso_gpu_t g;
    if (!iso_gpu_init(&g)) { iso_gpu_free(&g); UTEST_SKIP(iso_gpu_no_device()); }
    md_gto_iso_gpu_initialize(g.dev);
    t0 = md_tick_now();
    const md_gto_iso_field_t* fl[1] = {field};
    md_gto_iso_gpu_set_t* set = md_gto_iso_gpu_set_create(g.s, &basis, xyz, 3 * sizeof(float), fl, 1, md_get_heap_allocator());
    md_gpu_stream_sync(g.s);
    printf("  GPU set: %.2f s%s\n", md_tick_to_seconds(md_tick_now() - t0), set ? "" : " FAILED");
    md_gpu_memory_stats_t ms;
    md_gpu_memory_stats(g.dev, MD_GPU_MEM_DEVICE, &ms);
    printf("  device memory: %.1f MB in use, %.1f MB peak, %.1f MB temp\n", ms.bytes_in_use / 1e6, ms.bytes_peak_in_use / 1e6, ms.bytes_temp / 1e6);
    if (set) {
        md_gto_iso_gpu_desc_t gd = {.cell_size = h, .project_iterations = 8, .relax_iterations = 3};
        gd.iso[0] = iso;
        static md_gto_iso_gpu_profile_t prof;
        static char text[1 << 16];
        md_gto_iso_gpu_profile_set(getenv("GTO_ISO_LARGE_PROFILE") ? &prof : NULL);
        t0 = md_tick_now();
        md_gto_iso_gpu_mesh_t m;
        const bool ok = md_gto_iso_gpu_extract(g.s, set, &gd, &m);
        printf("  extract (h %.2f, iso %g): %.2f s, %s: %u vertices, %u triangles\n", h, iso, md_tick_to_seconds(md_tick_now() - t0),
               ok ? "ok" : "FAILED", ok ? m.vertex_count : 0, ok ? m.index_count / 3 : 0);
        md_gto_iso_gpu_profile_set(NULL);
        if (getenv("GTO_ISO_LARGE_PROFILE")) { md_gto_iso_gpu_profile_print(&prof, text, sizeof(text)); printf("%s", text); }
        md_gpu_memory_stats(g.dev, MD_GPU_MEM_DEVICE, &ms);
        printf("  device memory: %.1f MB in use, %.1f MB peak, %.1f MB temp\n", ms.bytes_in_use / 1e6, ms.bytes_peak_in_use / 1e6, ms.bytes_temp / 1e6);
        EXPECT_TRUE(ok);
        if (ok) md_gto_iso_gpu_mesh_free(g.s, &m);
        md_gto_iso_gpu_set_destroy(g.s, set);
    }
    md_gto_iso_field_destroy(field);
    md_gto_iso_gpu_shutdown();
    iso_gpu_free(&g);
    free(xyz);
    qm_test_free(&t);
}

#endif // MD_ENABLE_GPU
