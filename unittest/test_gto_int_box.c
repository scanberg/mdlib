#include "utest.h"
#include "vlx_test_util.h"

#include <md_gto_int.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

#include <math.h>
#include <stdio.h>

// Enclosures of the electrostatic potential over boxes (md_gto_int_potential_box): that they hold
// what the potential is inside them, and how much of a lattice they leave to evaluate for the
// isosurfaces V = +-iso, against the cells that surface actually crosses.

typedef struct {
    md_gto_basis_t basis;
    float*  atom_xyz;   // bohr
    double* atom_z;
    size_t  num_atoms;
    double* D;
    size_t  num_ao;
} box_test_qm_t;

static bool box_test_load(box_test_qm_t* m, md_gto_int_charges_t* q, vlx_test_t* t, const char* file) {
    if (!vlx_test_load(t, md_str_from_cstr(file), MD_MEGABYTES(64))) return false;
    md_memset(m, 0, sizeof(*m));
    if (!qm_test_basis(&m->basis, t)) return false;
    m->num_ao = md_gto_basis_num_ao(&m->basis);
    const size_t cap = 3 * 1024;
    m->atom_xyz = (float*)md_alloc(t->alloc, sizeof(float) * cap);
    m->num_atoms = vlx_test_atom_xyz_bohr(m->atom_xyz, cap, t);
    if (!m->num_atoms) return false;
    m->atom_z = (double*)md_alloc(t->alloc, sizeof(double) * m->num_atoms);
    if (qm_test_series(m->atom_z, m->num_atoms, t, MD_STR_LIT("qm/atom/nuclear_charge")) != m->num_atoms) return false;
    size_t dim = 0;
    m->D = qm_test_matrix(t, MD_STR_LIT("orbital/total/density"), &dim);
    if (!m->D) m->D = qm_test_matrix(t, MD_STR_LIT("orbital/alpha/density"), &dim);
    if (!m->D || dim != m->num_ao) return false;
    md_gto_int_charges_desc_t desc = {
        .basis = &m->basis, .atom_xyz = m->atom_xyz, .density_matrix = m->D, .density_scale = -1.0,
        .point_xyz = m->atom_xyz, .point_charge = m->atom_z, .num_points = m->num_atoms,
    };
    return md_gto_int_charges_init(q, &desc, t->alloc);
}

static double box_test_V(const md_gto_int_charges_t* q, double x, double y, double z) {
    const float c[3] = { (float)x, (float)y, (float)z };
    double V;
    md_gto_int_potential_xyz(&V, NULL, c, 1, 0, q);
    return V;
}

static double box_test_rand(uint64_t* s) {
    *s = *s * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*s >> 11) / (double)(1ULL << 53);
}

// Random boxes around the molecule, from inside it to far out, sizes from a lattice cell to a few
// bohr: V sampled inside each must lie in its enclosure
UTEST(gto_int_box, encloses_the_potential) {
    vlx_test_t t;
    box_test_qm_t m;
    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(box_test_load(&m, &q, &t, MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"));
    uint64_t seed = 7;
    int tested = 0, skipped = 0;
    double worst_excess = 0.0;
    for (int b = 0; b < 400; ++b) {
        const double size = 0.1 * pow(40.0, box_test_rand(&seed));
        double lo[3], hi[3];
        for (int a = 0; a < 3; ++a) {
            lo[a] = -8.0 + 16.0 * box_test_rand(&seed);
            hi[a] = lo[a] + size;
        }
        double vlo, vhi, wlo, whi;
        if (!md_gto_int_potential_box(&vlo, &vhi, lo, hi, NULL, NULL, &q)) { skipped++; continue; }
        ASSERT_TRUE(md_gto_int_potential_box3(&wlo, &whi, lo, hi, &q));
        tested++;
        double smin = INFINITY, smax = -INFINITY;
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) for (int k = 0; k < 6; ++k) {
            const double x = lo[0] + (hi[0] - lo[0]) * i / 5.0, y = lo[1] + (hi[1] - lo[1]) * j / 5.0, z = lo[2] + (hi[2] - lo[2]) * k / 5.0;
            const double v = box_test_V(&q, x, y, z);
            smin = fmin(smin, v);
            smax = fmax(smax, v);
        }
        // float positions in the sampler: allow what moving a point by a float ulp does
        const double tol = 1e-9 + 1e-6 * fmax(fabs(smin), fabs(smax));
        EXPECT_LE(vlo, smin + tol);
        EXPECT_GE(vhi, smax - tol);
        EXPECT_LE(wlo, smin + tol);
        EXPECT_GE(whi, smax - tol);
        worst_excess = fmax(worst_excess, ((vhi - vlo) - (smax - smin)) / fmax(1e-12, smax - smin));
    }
    printf("enclosures: %d boxes, %d holding a nucleus skipped; widest relative excess over the sampled range %.2f\n", tested, skipped, worst_excess);
    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
}

// ---------------------------------------------------------------------------
// Culling a lattice
// ---------------------------------------------------------------------------
typedef struct {
    const md_gto_int_charges_t* q;
    double iso;
    double h;           // lattice spacing
    double origin[3];
    int    n;           // cells per axis of the root (a power of two)
    int    leaf_cells;  // cells per axis of a leaf box
    int    hessian_cells;   // a box of at most this many cells takes its parent's Hessian
    int    order;           // 2: md_gto_int_potential_box, 3: md_gto_int_potential_box3
    long   hessians;
    uint8_t* kept;      // [n^3] cells inside a kept leaf
    long   tests;
    long   leaves_kept;
    double seconds_tests;
} box_test_cull_t;

static bool box_test_meets(double lo, double hi, double iso) {
    return (lo <= iso && hi >= iso) || (lo <= -iso && hi >= -iso);
}

static void box_test_visit(box_test_cull_t* c, int cells, int ix, int iy, int iz, const double* Hlo_in, const double* Hhi_in) {
    double lo[3], hi[3];
    const int i0[3] = { ix, iy, iz };
    for (int a = 0; a < 3; ++a) {
        lo[a] = c->origin[a] + i0[a] * c->h;
        hi[a] = lo[a] + cells * c->h;
    }
    double vlo, vhi;
    double Hlo[6], Hhi[6];
    const double* hl = Hlo_in;
    const double* hh = Hhi_in;
    c->tests++;
    const md_tick_t t0 = md_tick_now();
    bool bounded = true;
    if (c->order == 3) {
        bounded = md_gto_int_potential_box3(&vlo, &vhi, lo, hi, c->q);
        hl = hh = NULL;
    } else if (!hl || cells > c->hessian_cells) {
        // its own: at the top, at the size it is passed down from, or below a box that had none
        c->hessians++;
        bounded = md_gto_int_hessian_box(Hlo, Hhi, lo, hi, c->q);
        hl = bounded ? Hlo : NULL;
        hh = bounded ? Hhi : NULL;
    }
    if (bounded && c->order != 3) bounded = md_gto_int_potential_box(&vlo, &vhi, lo, hi, hl, hh, c->q);
    c->seconds_tests += md_tick_to_seconds(md_tick_now() - t0);
    if (bounded && !box_test_meets(vlo, vhi, c->iso)) return;
    if (cells <= c->leaf_cells) {
        c->leaves_kept++;
        for (int z = iz; z < iz + cells; ++z) for (int y = iy; y < iy + cells; ++y) for (int x = ix; x < ix + cells; ++x) {
            c->kept[((size_t)z * c->n + y) * c->n + x] = 1;
        }
        return;
    }
    const int half = cells / 2;
    // the children take this box's Hessian once they are small enough
    const bool pass = bounded && hl && cells <= 2 * c->hessian_cells;
    for (int k = 0; k < 8; ++k) {
        box_test_visit(c, half, ix + (k & 1) * half, iy + ((k >> 1) & 1) * half, iz + ((k >> 2) & 1) * half, pass ? hl : NULL, pass ? hh : NULL);
    }
}

static void box_test_cull_report(const char* file, double iso, double h, double pad) {
    vlx_test_t t;
    box_test_qm_t m;
    md_gto_int_charges_t q = {0};
    if (!box_test_load(&m, &q, &t, file)) { printf("could not load %s\n", file); return; }

    double amin[3] = { INFINITY, INFINITY, INFINITY }, amax[3] = { -INFINITY, -INFINITY, -INFINITY };
    for (size_t i = 0; i < m.num_atoms; ++i) for (int a = 0; a < 3; ++a) {
        amin[a] = fmin(amin[a], m.atom_xyz[3 * i + a]);
        amax[a] = fmax(amax[a], m.atom_xyz[3 * i + a]);
    }
    double ext = 0;
    for (int a = 0; a < 3; ++a) ext = fmax(ext, amax[a] - amin[a] + 2 * pad);
    int n = 1;
    while (n * h < ext) n *= 2;
    double origin[3];
    for (int a = 0; a < 3; ++a) origin[a] = 0.5 * (amin[a] + amax[a]) - 0.5 * n * h;

    // The reference: every lattice point, every cell the surface crosses
    const size_t np = (size_t)(n + 1) * (n + 1) * (n + 1);
    float* V = (float*)malloc(sizeof(float) * np);
    const md_tick_t t0 = md_tick_now();
    for (int z = 0; z <= n; ++z) for (int y = 0; y <= n; ++y) for (int x = 0; x <= n; ++x) {
        V[((size_t)z * (n + 1) + y) * (n + 1) + x] = (float)box_test_V(&q, origin[0] + x * h, origin[1] + y * h, origin[2] + z * h);
    }
    const double seconds_full = md_tick_to_seconds(md_tick_now() - t0);
    const size_t nc = (size_t)n * n * n;
    uint8_t* crossing = (uint8_t*)calloc(nc, 1);
    long num_crossing = 0;
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        uint8_t above_p = 0, below_p = 0, above_m = 0, below_m = 0;
        for (int k = 0; k < 8; ++k) {
            const float v = V[((size_t)(z + ((k >> 2) & 1)) * (n + 1) + (y + ((k >> 1) & 1))) * (n + 1) + (x + (k & 1))];
            if (v >= iso) above_p = 1; else below_p = 1;
            if (v <= -iso) below_m = 1; else above_m = 1;
        }
        if ((above_p && below_p) || (above_m && below_m)) {
            crossing[((size_t)z * n + y) * n + x] = 1;
            num_crossing++;
        }
    }

    printf("\n%s: %zu atoms, %u gaussians, iso +-%.3g, lattice %d^3 cells of %.2f bohr, %zu points (%.1f s, %.1f us a point)\n",
           file, m.num_atoms, q.num_gaussians, iso, n, h, np, seconds_full, 1e6 * seconds_full / np);
    printf("  crossing cells %ld (%.3f%% of the lattice)\n", num_crossing, 100.0 * num_crossing / nc);
    const int leaves[] = { 4, 2, 1 };
    // The best any culling can do at each leaf size: the leaves that hold a crossing cell
    for (int L = 4; L >= 1; L /= 2) {
        long leaves = 0;
        for (int z = 0; z < n; z += L) for (int y = 0; y < n; y += L) for (int x = 0; x < n; x += L) {
            bool any = false;
            for (int k = 0; k < L * L * L && !any; ++k) any = crossing[((size_t)(z + k / (L * L)) * n + (y + (k / L) % L)) * n + (x + k % L)];
            leaves += any;
        }
        printf("  ideal leaf %d^3: %ld leaves, %ld cells (%.2fx crossing), about %.2f%% of the lattice's corners\n", L, leaves, leaves * L * L * L,
               (double)leaves * L * L * L / MD_MAX(1, num_crossing), 100.0 * leaves * (L + 1) * (L + 1) * (L + 1) / np);
    }
    const int orders[] = { 2, 3 };
    for (int li = 0; li < 3; ++li) for (int oi = 0; oi < 2; ++oi) {
        box_test_cull_t c = { .q = &q, .iso = iso, .h = h, .n = n, .leaf_cells = leaves[li], .order = orders[oi] };
        for (int a = 0; a < 3; ++a) c.origin[a] = origin[a];
        c.kept = (uint8_t*)calloc(nc, 1);
        box_test_visit(&c, n, 0, 0, 0, NULL, NULL);
        long kept = 0, missed = 0;
        for (size_t i = 0; i < nc; ++i) {
            kept += c.kept[i];
            if (crossing[i] && !c.kept[i]) missed++;
        }
        // points to evaluate: the corners of the kept cells
        long corners = 0;
        for (int z = 0; z <= n; ++z) for (int y = 0; y <= n; ++y) for (int x = 0; x <= n; ++x) {
            bool any = false;
            for (int k = 0; k < 8 && !any; ++k) {
                const int cx = x - (k & 1), cy = y - ((k >> 1) & 1), cz = z - ((k >> 2) & 1);
                if (cx < 0 || cy < 0 || cz < 0 || cx >= n || cy >= n || cz >= n) continue;
                any = c.kept[((size_t)cz * n + cy) * n + cx];
            }
            corners += any;
        }
        const double us_point = 1e6 * seconds_full / np;
        printf("  leaf %d^3, order %d: %6ld tests (%5ld Hessians, %4.0f point evals' worth in all a test), kept %6ld cells (%5.2fx crossing), %7ld corners (%5.2f%%), missed %ld; work %6.2f%% of the full lattice\n",
               leaves[li], orders[oi], c.tests, c.hessians, 1e6 * c.seconds_tests / c.tests / us_point,
               kept, (double)kept / MD_MAX(1, num_crossing), corners, 100.0 * corners / np,
               missed, 100.0 * (corners * us_point + c.seconds_tests * 1e6) / (np * us_point));
        free(c.kept);
    }
    free(crossing);
    free(V);
    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
}

UTEST(gto_int_box, culling_h2o) {
    box_test_cull_report(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5", 0.05, 0.2, 8.0);
}

UTEST(gto_int_box, culling_amide) {
    box_test_cull_report(MD_UNITTEST_DATA_DIR "/vlx/amide.h5", 0.05, 0.3, 8.0);
}
