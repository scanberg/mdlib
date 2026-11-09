#include <md_mbis.h>

#include <md_gto.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_executor.h>
#include <core/md_log.h>
#include <core/md_simd.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// The defaults, measured against a converged grid (100 radial x Lebedev 59): the charges within ~2e-5 e for molecules
// of H to S, ~5e-5 e in a dense system (a box of 125 waters), where the angular order is what decides it (Lebedev 29
// gave 1.2e-4 there). The threshold leaves ~3e-6 e of iteration error, an order below that.
#define MBIS_DEFAULT_RADIAL     48
#define MBIS_DEFAULT_ANGULAR    35
#define MBIS_DEFAULT_THRESHOLD  1.0e-6
#define MBIS_DEFAULT_MAX_ITER   1000

// Bound of |phi| beyond which a shell is not evaluated at a point.
#define MBIS_AO_EPS             1.0e-10
// A point whose density is below this is dropped: an atom's grid reaches ~14 bohr (48 radial points), a sphere
// of ~1e4 bohr^3, so at most ~1e-10 e an atom goes.
#define MBIS_RHO_FLOOR          1.0e-14
// A pro-atom shell is left out at a point where it is below this fraction of the pro-molecule density there: every
// weight then moves by at most that, relative (measured: ~2.5e-6 e on the charges; 1e-7 gave 1.8e-5 e).
#define MBIS_TERM_EPS           1.0e-8
// Stratmann, Scuseria and Frisch, Chem. Phys. Lett. 257, 213 (1996): the cell function is exactly 1 or 0 outside
// |mu| < a, which is what makes the screening below exact. They use 0.64; 0.5 makes the cells reach less far, and
// fewer of them overlap a point in dense matter, which halves the cost of the weights at the same accuracy
// (measured against the converged grid, molecules and a box of waters alike; 0.4 began to cost accuracy).
#define MBIS_SSF_A              0.5

#define MBIS_MAX_CHUNKS         64

// A sphere's points are batched by the nearest of these directions: the 6 axes and the 8 cube corners
#define MBIS_NUM_PATCHES        14

// The iterations: Anderson acceleration over the last MBIS_ANDERSON_DEPTH images; first on the coarse grid (every
// MBIS_COARSE_STRIDE-th radial point) until the residual is below MBIS_COARSE_THRESHOLD, unless started from a
// previous answer. Measured: depth 3 to 8, a stride of 4 or a coarse threshold of 1e-4 all within 10% of these.
#define MBIS_ANDERSON_DEPTH     5
#define MBIS_COARSE_STRIDE      2
#define MBIS_COARSE_THRESHOLD   1.0e-5

// Radial points closer to their nucleus than this fraction of the distance to the nearest atom use an angular
// grid of at most MBIS_INNER_ANGULAR (MBIS_INNER_LONE_ATOM bohr stands in for that distance for a lone atom).
// Measured on H to S: the charges move by less than 1e-6 e against the full grid everywhere, at ~40% fewer points.
#define MBIS_INNER_FRACTION     0.25
#define MBIS_INNER_ANGULAR      11
#define MBIS_INNER_LONE_ATOM    4.0

// exp(x) for every x of an array, n a multiple of 4, x <= 0 (below -708 taken as -708, ~3e-308). x = k ln2 + r,
// |r| <= ln2/2, Taylor to r^12 (relative error ~3e-16, as libm), 2^k assembled in the exponent bits. Four at a
// time, with the same arithmetic in every lane, so a value does not depend on where in the array it sits.
// Explicit rather than left to the compiler: a libm call per element does not vectorize on every compiler mdlib
// is built with, and the exponentials are most of what the density and the iterations cost.
static void mbis_exp(double* y, const double* x, size_t n) {
    const md_256d lo      = md_mm256_set1_pd(-708.0);
    const md_256d log2e   = md_mm256_set1_pd(1.4426950408889634);
    const md_256d shifter = md_mm256_set1_pd(6755399441055744.0);   // 1.5 * 2^52: rounds to an integer kept in the low bits
    const md_256d ln2_hi  = md_mm256_set1_pd(6.93147180369123816490e-01);
    const md_256d ln2_lo  = md_mm256_set1_pd(1.90821492927058770002e-10);
    const md_256i bias    = md_mm256_set1_epi64(0x4338000000000000LL - 1023);
    static const double c[13] = {1.0 / 479001600.0, 1.0 / 39916800.0, 1.0 / 3628800.0, 1.0 / 362880.0, 1.0 / 40320.0,
                                 1.0 / 5040.0, 1.0 / 720.0, 1.0 / 120.0, 1.0 / 24.0, 1.0 / 6.0, 0.5, 1.0, 1.0};
    for (size_t i = 0; i < n; i += 4) {
        const md_256d v = md_mm256_max_pd(md_mm256_loadu_pd(x + i), lo);
        const md_256d t = md_mm256_fmadd_pd(v, log2e, shifter);
        const md_256d k = md_mm256_sub_pd(t, shifter);
        const md_256d r = md_mm256_fnmadd_pd(k, ln2_lo, md_mm256_fnmadd_pd(k, ln2_hi, v));
        md_256d p = md_mm256_set1_pd(c[0]);
        for (int j = 1; j < 13; ++j) p = md_mm256_fmadd_pd(p, r, md_mm256_set1_pd(c[j]));
        const md_256i e = md_mm256_slli_epi64(md_mm256_sub_epi64(md_mm256_castpd_si256(t), bias), 52);
        md_mm256_storeu_pd(y + i, md_mm256_mul_pd(p, md_mm256_castsi256_pd(e)));
    }
}

// ---------------------------------------------------------------------------
// The basis, prepared for point evaluation
// ---------------------------------------------------------------------------

typedef struct mbis_shell_t {
    double   c[3];          // centre, bohr
    double   r2_cut;        // |phi| < MBIS_AO_EPS beyond sqrt(r2_cut)
    uint32_t l;
    uint32_t prim_off;
    uint32_t num_prim;
    uint32_t ao;            // first AO
    uint32_t num_ao;        // (l+1)(l+2)/2
} mbis_shell_t;

typedef struct mbis_cart_t {
    uint8_t ijk[15][3];
    double  nrm[15];
} mbis_cart_t;

// An angular grid, its points grouped into patches of neighbouring directions: a radial sphere's batches
typedef struct mbis_sphere_t {
    double*  u;             // [n][3] unit vectors, patch by patch
    double*  w;             // [n], summing to 4 pi
    uint32_t n;
    uint32_t num_patches;
    uint32_t patch_off[MBIS_NUM_PATCHES + 1];
    uint32_t max_patch;     // points in the largest patch
} mbis_sphere_t;

// Bound of |phi| over a shell's AOs at distance r: f(i,j,k) <= 1 and |x^i y^j z^k| <= r^l.
static double mbis_shell_bound(const double* alpha, const double* coeff, uint32_t n, uint32_t l, double r) {
    double s = 0.0;
    const double rl = pow(r, (double)l);
    for (uint32_t p = 0; p < n; ++p) s += fabs(coeff[p]) * rl * exp(-alpha[p] * r * r);
    return s;
}

// The radius beyond which the bound stays below eps: the bound decreases monotonically past the most
// diffuse primitive's peak, so step out from there and bisect.
static double mbis_shell_cutoff(const double* alpha, const double* coeff, uint32_t n, uint32_t l, double eps) {
    double amin = HUGE_VAL;
    for (uint32_t p = 0; p < n; ++p) amin = MD_MIN(amin, alpha[p]);
    if (!(amin > 0.0) || n == 0) return 0.0;
    double lo = sqrt(MD_MAX((double)l, 0.5) / (2.0 * amin));
    double hi = lo + 1.0;
    while (mbis_shell_bound(alpha, coeff, n, l, hi) > eps) {
        lo = hi;
        hi *= 2.0;
        if (hi > 1.0e4) return hi;
    }
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (mbis_shell_bound(alpha, coeff, n, l, mid) > eps) lo = mid; else hi = mid;
    }
    return hi;
}

// ---------------------------------------------------------------------------
// Quadrature
// ---------------------------------------------------------------------------

#include "md_mbis_lebedev.inl"

// Mura and Knowles, J. Chem. Phys. 104, 9848 (1996): r = -a ln(1 - x^3), x = i / (n + 1). The weight
// carries the r^2 of the volume element. a = 7 for the s block metals, 5 otherwise.
static void mbis_radial_grid(double* r, double* w, uint32_t n, uint8_t z) {
    const bool sblock = z == 3 || z == 4 || z == 11 || z == 12 || z == 19 || z == 20 || z == 37 || z == 38 || z == 55 || z == 56 || z == 87 || z == 88;
    const double a = sblock ? 7.0 : 5.0;
    for (uint32_t i = 0; i < n; ++i) {
        const double x = (i + 1.0) / (n + 1.0);
        const double x3 = x * x * x;
        r[i] = -a * log(1.0 - x3);
        w[i] = 3.0 * a * x * x / (1.0 - x3) / (n + 1.0) * r[i] * r[i];
    }
}

// The Lebedev grid integrating the spherical harmonics up to degree 'order' exactly (the smallest tabulated one that
// does), expanded from its orbits and sorted into patches.
static bool mbis_sphere_init(mbis_sphere_t* sph, uint32_t order) {
    size_t g = 0;
    while (g + 1 < MD_COUNTOF(mbis_lebedev) && mbis_lebedev[g].order < order) ++g;
    const uint32_t n = mbis_lebedev[g].num_points;
    double* u = (double*)malloc(sizeof(double) * 3 * n);
    double* w = (double*)malloc(sizeof(double) * n);
    sph->u = (double*)malloc(sizeof(double) * 3 * n);
    sph->w = (double*)malloc(sizeof(double) * n);
    if (!u || !w || !sph->u || !sph->w) { free(u); free(w); return false; }

    // Every orbit: its canonical point under all permutations and sign changes, each distinct point once
    uint32_t m = 0;
    static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (uint32_t o = 0; o < mbis_lebedev[g].num_orbits; ++o) {
        const double* c = mbis_lebedev[g].orbit[o];
        const double len = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);   // the tables carry 12 digits
        const uint32_t first = m;
        for (int pi = 0; pi < 6; ++pi) {
            for (int sg = 0; sg < 8; ++sg) {
                double v[3];
                for (int k = 0; k < 3; ++k) v[k] = ((sg >> k) & 1 ? -1.0 : 1.0) * c[perm[pi][k]] / len;
                bool seen = false;
                for (uint32_t q = first; q < m && !seen; ++q) {
                    seen = fabs(u[3 * q] - v[0]) < 1e-9 && fabs(u[3 * q + 1] - v[1]) < 1e-9 && fabs(u[3 * q + 2] - v[2]) < 1e-9;
                }
                if (seen) continue;
                if (m == n) { free(u); free(w); return false; }
                u[3 * m] = v[0]; u[3 * m + 1] = v[1]; u[3 * m + 2] = v[2];
                w[m] = 4.0 * 3.14159265358979323846 * c[3];
                ++m;
            }
        }
    }
    if (m != n) { free(u); free(w); return false; }

    // Patches: the nearest of the axes and cube corners, points in patch order
    const double cr = 1.0 / sqrt(3.0);
    const double dir[MBIS_NUM_PATCHES][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
        {cr, cr, cr}, {cr, cr, -cr}, {cr, -cr, cr}, {cr, -cr, -cr}, {-cr, cr, cr}, {-cr, cr, -cr}, {-cr, -cr, cr}, {-cr, -cr, -cr}};
    uint8_t patch[6000];
    if (n > MD_COUNTOF(patch)) { free(u); free(w); return false; }
    uint32_t count[MBIS_NUM_PATCHES] = {0};
    for (uint32_t j = 0; j < n; ++j) {
        int best = 0;
        double bd = -2.0;
        for (int k = 0; k < MBIS_NUM_PATCHES; ++k) {
            const double d = u[3 * j] * dir[k][0] + u[3 * j + 1] * dir[k][1] + u[3 * j + 2] * dir[k][2];
            if (d > bd + 1e-12) { bd = d; best = k; }
        }
        patch[j] = (uint8_t)best;
        count[best]++;
    }
    sph->n = n;
    sph->num_patches = MBIS_NUM_PATCHES;
    sph->patch_off[0] = 0;
    sph->max_patch = 0;
    for (int k = 0; k < MBIS_NUM_PATCHES; ++k) {
        sph->patch_off[k + 1] = sph->patch_off[k] + count[k];
        sph->max_patch = MD_MAX(sph->max_patch, count[k]);
    }
    uint32_t fill[MBIS_NUM_PATCHES];
    for (int k = 0; k < MBIS_NUM_PATCHES; ++k) fill[k] = sph->patch_off[k];
    for (uint32_t j = 0; j < n; ++j) {
        const uint32_t d = fill[patch[j]]++;
        memcpy(sph->u + 3 * d, u + 3 * j, sizeof(double) * 3);
        sph->w[d] = w[j];
    }
    free(u);
    free(w);
    return true;
}

// ---------------------------------------------------------------------------
// Becke partition (Stratmann-Scuseria-Frisch)
// ---------------------------------------------------------------------------

static inline double mbis_ssf_s(double mu) {
    if (mu <= -MBIS_SSF_A) return 1.0;
    if (mu >= MBIS_SSF_A) return 0.0;
    const double x = mu / MBIS_SSF_A, x2 = x * x;
    const double g = x * (35.0 + x2 * (-35.0 + x2 * (21.0 - 5.0 * x2))) / 16.0;
    return 0.5 * (1.0 - g);
}

static inline double mbis_dist(const double* a, const double* b) {
    const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

typedef struct mbis_neighbor_t {
    uint32_t atom;
    double   d;
} mbis_neighbor_t;

static int mbis_neighbor_cmp(const void* a, const void* b) {
    const mbis_neighbor_t* x = (const mbis_neighbor_t*)a;
    const mbis_neighbor_t* y = (const mbis_neighbor_t*)b;
    if (x->d < y->d) return -1;
    if (x->d > y->d) return 1;
    return (x->atom > y->atom) - (x->atom < y->atom);
}

// Weight of atom B's cell at p, which lies r from B. 'nb' holds, per atom, every atom sorted by distance from it
// (itself first). With k = (1+a)/(1-a) the SSF cells are exactly 0 or 1 beyond simple bounds, which keep the work
// per point independent of the size of the system:
// - Stratmann et al.: B's cell is 1 and every other 0 within 0.5 (1 - a) of the distance to B's nearest neighbour.
// - The cell of D vanishes once r_D >= k r_min, r_min the distance to the nearest atom: the atoms that may own a
//   cell at p are within k r_min of it, so within r + k r_min of B, a bound that tightens as r_min is found.
//   B's own vanishes as soon as one atom is k times closer.
// - Most of those are zeroed by the nearest atom alone (mu >= a), which is tested first.
// - s(mu_DC) = 1 once r_C >= k r_D, so the product for D runs over the atoms within (k + 1) r_D of D.
// The distance from the current point to each atom, computed once per point however many cells ask for it: entries
// are valid while their stamp is the point's.
typedef struct mbis_dist_cache_t {
    double*   r;        // [num_atoms]
    uint32_t* stamp;    // [num_atoms]
    uint32_t  cur;
} mbis_dist_cache_t;

static inline double mbis_dist_cached(mbis_dist_cache_t* dc, const double* p, const dvec3_t* xyz, uint32_t atom) {
    if (dc->stamp[atom] != dc->cur) {
        dc->stamp[atom] = dc->cur;
        dc->r[atom] = mbis_dist(p, xyz[atom].elem);
    }
    return dc->r[atom];
}

static double mbis_becke_weight(const double* p, double r, uint32_t b, const mbis_neighbor_t* nb_all, size_t na, const dvec3_t* xyz,
                                const double* inv_dist, uint32_t* cand, double* cand_r, mbis_dist_cache_t* dc) {
    const double k = (1.0 + MBIS_SSF_A) / (1.0 - MBIS_SSF_A);
    const mbis_neighbor_t* nb = nb_all + (size_t)b * na;
    if (na <= 1 || r <= 0.5 * (1.0 - MBIS_SSF_A) * nb[1].d) return 1.0;
    if (++dc->cur == 0) {   // the stamp wrapped: every entry is stale
        memset(dc->stamp, 0, sizeof(uint32_t) * na);
        dc->cur = 1;
    }
    dc->stamp[b] = dc->cur;
    dc->r[b] = r;

    double rmin = r;
    uint32_t c0 = b;
    double reach = r + k * rmin;
    size_t n = 0;
    for (size_t i = 0; i < na && nb[i].d <= reach; ++i) {
        const double rc = mbis_dist_cached(dc, p, xyz, nb[i].atom);
        if (rc * k <= r) return 0.0;
        if (rc < rmin) {
            rmin = rc;
            c0 = nb[i].atom;
            reach = r + k * rmin;
        }
        cand[n] = nb[i].atom;
        cand_r[n] = rc;
        ++n;
    }
    if (n <= 1) return 1.0;

    double cell_b = 0.0, sum = 0.0;
    for (size_t di = 0; di < n; ++di) {
        const double rd = cand_r[di];
        if (rd >= k * rmin) continue;
        const uint32_t d = cand[di];
        double cell = 0.0;
        if (d == c0 || (rd - rmin) * inv_dist[(size_t)d * na + c0] < MBIS_SSF_A) {
            const mbis_neighbor_t* nbd = nb_all + (size_t)d * na;
            const double lim = (k + 1.0) * rd;
            cell = 1.0;
            for (size_t j = 1; j < na && nbd[j].d < lim && cell > 0.0; ++j) {
                const double rc = mbis_dist_cached(dc, p, xyz, nbd[j].atom);
                if (rc >= k * rd) continue;
                cell *= mbis_ssf_s((rd - rc) / nbd[j].d);
            }
        }
        if (di == 0) {
            if (!(cell > 0.0)) return 0.0;
            cell_b = cell;
        }
        sum += cell;
    }
    return sum > 0.0 ? cell_b / sum : 0.0;
}

// ---------------------------------------------------------------------------
// The run
// ---------------------------------------------------------------------------

// The grid points of one atom kept for the iterations: radial and angular index, and the mass, the
// quadrature weight times the Becke weight times the density.
typedef struct mbis_atom_points_t {
    uint16_t* rad;
    uint16_t* ang;
    double*   mass;
    uint32_t  count;
    double    electrons;
} mbis_atom_points_t;

typedef struct mbis_thread_t {
    double*   phib;         // [batch size][cap_ao]: the AOs at a batch's points (0 where a shell is out of reach)
    double*   drow;         // [cap_ao + 4]: a row of D over a batch's columns
    double*   zrow;         // [cap_ao + 4]: zeros, the AOs of a point that is not there
    double*   brho;         // [batch size]: the density at a batch's points
    uint32_t* act_p;        // [batch size]: the points where a row's AO is not 0
    uint8_t*  reach;        // [batch size][num_shells]: which shells are evaluated at which of them
    double*   bx;           // [batch size + 4] x 7: the batch's points, and per shell r^2, the radial part,
    double*   by;           //   a primitive's exponent and exponential
    double*   bz;
    double*   br2;
    double*   bR;
    double*   barg;
    double*   bexp;
    uint32_t* shells;       // [num_shells]: a batch's shells
    uint32_t* loc;          // [num_shells]: where their AOs start in a point's row of phib
    uint32_t* run_beg;      // [num_shells]: a batch's runs of consecutive AOs: where they start in its rows,
    uint32_t* run_col;      // [num_shells]: their first AO,
    uint32_t* run_len;      // [num_shells]: and their length
    uint32_t* surv_j;       // [batch size]: a batch's points inside the cell
    double*   surv_w;       // [batch size]: and their cell weights
    mbis_dist_cache_t dc;   // for the cell weights
    uint32_t* cand;         // [num_atoms]
    double*   cand_r;       // [num_atoms]
    uint32_t* sub;          // [num_atoms]: the atoms whose pro-atoms reach a radial sphere
    double*   sub_rc2;      // [num_atoms]: and the squared distance beyond which they do not count
    uint32_t* term_k;       // [num_shells]: a point's pro-atom terms: their shell,
    double*   term_v;       // [num_shells + 4]: value,
    double*   term_x;       // [num_shells + 4]: ln(pre) - S r, whose exponential the value is,
    double*   term_r;       // [num_shells]: and distance to their atom
    size_t    cap_ao;
} mbis_thread_t;

typedef struct mbis_ctx_t {
    const md_mbis_desc_t* desc;
    const dvec3_t*        xyz;
    size_t                num_atoms;
    size_t                num_ao;

    mbis_shell_t*         shells;
    uint32_t*             atom_shell_off; // [num_atoms + 1]: atom a's shells are atom_shell[atom_shell_off[a] .. a + 1]
    uint32_t*             atom_shell;     // [num_shells]
    double                shell_reach;    // the longest any shell reaches
    uint32_t              num_shells;
    double*               alpha;
    double*               prim_r2;        // a primitive is left out beyond sqrt(prim_r2): the tight ones of a contraction
                                          // vanish long before the shell does
    double*               coeff;
    mbis_cart_t           cart[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];

    uint32_t              num_rad;
    double*               rad_r;          // [num_atoms][num_rad]
    double*               rad_w;          // [num_atoms][num_rad]
    mbis_sphere_t         sphere[2];      // the full angular grid, and the coarse one for the points near a nucleus
    double*               inner_r;        // [num_atoms]: radial points closer than this use sphere[1]

    mbis_neighbor_t*      nb;             // [num_atoms][num_atoms], sorted by distance from the row's atom
    double*               inv_dist;       // [num_atoms][num_atoms], 1 / |R_a - R_b| (0 on the diagonal)
    mbis_atom_points_t*   pts;            // [num_atoms]
    mbis_thread_t*        thr;            // [num_threads]
    bool                  oom;
    volatile MD_U64       atoms_done;     // grid phase, for the progress estimate

    // Pro-atoms
    uint32_t*             shell_off;      // [num_atoms + 1]
    uint32_t              num_pro;
    double*               N;              // [num_pro]
    double*               S;              // [num_pro], 1 / sigma
    double*               pre;            // [num_pro], N S^3 / (8 pi)
    double*               logpre;         // [num_pro], its logarithm (-HUGE_VAL for an empty shell)
    double*               atom_logpre;    // [num_atoms], max ln(pre) of the atom's shells
    double*               atom_smin;      // [num_atoms], min S of the atom's shells
    uint32_t              num_chunks;
    double*               acc;            // [num_chunks][2 num_pro]
    bool                  coarse;         // iterating on the coarse grid
    double                logpre_max;     // the largest ln(pre) of any shell, this iteration
    double                reach_bound;    // 1 / the smallest S of any shell
} mbis_ctx_t;

// The angular grid of atom b's radial point at r: near the nucleus the density is close to spherical and
// the cell function is 1, so a coarse grid does there.
static inline const mbis_sphere_t* mbis_sphere_at(const mbis_ctx_t* ctx, uint32_t b, double r) {
    return r < ctx->inner_r[b] ? &ctx->sphere[1] : &ctx->sphere[0];
}

static bool mbis_thread_reserve(mbis_thread_t* t, size_t n, size_t batch_size) {
    if (n <= t->cap_ao) return true;
    size_t cap = MD_MAX(n, t->cap_ao * 2);
    const size_t stride = (cap + 3) & ~(size_t)3;
    double* phib = (double*)realloc(t->phib, sizeof(double) * (batch_size * stride + 4));
    if (phib) t->phib = phib;
    double* drow = (double*)realloc(t->drow, sizeof(double) * (cap + 4));
    if (drow) t->drow = drow;
    double* zrow = (double*)realloc(t->zrow, sizeof(double) * (stride + 4));
    if (zrow) t->zrow = zrow;
    if (!phib || !drow || !zrow) return false;
    memset(t->zrow, 0, sizeof(double) * (stride + 4));
    t->cap_ao = cap;
    return true;
}

// The density on atom b's grid and the points worth keeping. A batch is one patch of one radial sphere: points
// close together, so the shells reaching them are few even in a large system. Per batch the AOs are evaluated a
// shell at a time over all its points, and phi^T D phi is contracted reading D in place.
static void mbis_grid_task(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    mbis_ctx_t* ctx = (mbis_ctx_t*)arg;
    mbis_thread_t* t = &ctx->thr[thread_idx];
    const double* D = ctx->desc->density_matrix;
    const size_t N = ctx->num_ao;

    for (uint32_t b = beg; b < end; ++b) {
        if (ctx->desc->cancel && *ctx->desc->cancel) return;
        const double* cb = ctx->xyz[b].elem;
        const double* rr = ctx->rad_r + (size_t)b * ctx->num_rad;
        const double* rw = ctx->rad_w + (size_t)b * ctx->num_rad;
        mbis_atom_points_t* out = &ctx->pts[b];
        const size_t cap = (size_t)ctx->num_rad * ctx->sphere[0].n;
        out->rad  = (uint16_t*)malloc(sizeof(uint16_t) * cap);
        out->ang  = (uint16_t*)malloc(sizeof(uint16_t) * cap);
        out->mass = (double*)malloc(sizeof(double) * cap);
        out->count = 0;
        out->electrons = 0.0;
        if (!out->rad || !out->ang || !out->mass) { ctx->oom = true; return; }

        for (uint32_t ir = 0; ir < ctx->num_rad; ++ir) {
            const double r = rr[ir];
            const mbis_sphere_t* sph = mbis_sphere_at(ctx, b, r);

            for (uint32_t patch = 0; patch < sph->num_patches; ++patch) {
                const uint32_t j0 = sph->patch_off[patch], j1 = sph->patch_off[patch + 1];

                // The patch's points inside b's cell. The cell weight is exactly 0 over most of an atom's outer grid
                // in a dense system and costs a fraction of the density, so it goes first; what is left lies in
                // b's own cell, a small region whatever the radius.
                uint32_t nsurv = 0;
                double lo[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL}, hi[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
                for (uint32_t j = j0; j < j1; ++j) {
                    const double* u = sph->u + 3 * j;
                    const double p[3] = {cb[0] + r * u[0], cb[1] + r * u[1], cb[2] + r * u[2]};
                    const double wb = mbis_becke_weight(p, r, b, ctx->nb, ctx->num_atoms, ctx->xyz, ctx->inv_dist, t->cand, t->cand_r, &t->dc);
                    if (!(wb > 0.0)) continue;
                    t->surv_j[nsurv] = j;
                    t->surv_w[nsurv] = wb;
                    ++nsurv;
                    for (int ax = 0; ax < 3; ++ax) { lo[ax] = MD_MIN(lo[ax], p[ax]); hi[ax] = MD_MAX(hi[ax], p[ax]); }
                }
                if (nsurv == 0) continue;
                const double cc[3] = {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])};
                const double ext[3] = {0.5 * (hi[0] - lo[0]), 0.5 * (hi[1] - lo[1]), 0.5 * (hi[2] - lo[2])};
                const double rsurv = sqrt(ext[0] * ext[0] + ext[1] * ext[1] + ext[2] * ext[2]);

                // Shells that reach those points, and where their AOs sit in a point's row of phib. Their atoms are
                // found down b's neighbour list: an atom d from b is at least d - |cc - b| from the batch centre,
                // so past d > |cc - b| + rsurv + the longest reach of any shell, none is.
                uint32_t ns = 0;
                size_t nl = 0;
                {
                    const double dcb = mbis_dist(cc, cb);
                    const mbis_neighbor_t* nbr = ctx->nb + (size_t)b * ctx->num_atoms;
                    for (size_t i = 0; i < ctx->num_atoms && nbr[i].d <= dcb + rsurv + ctx->shell_reach; ++i) {
                        const uint32_t a = nbr[i].atom;
                        for (uint32_t q = ctx->atom_shell_off[a]; q < ctx->atom_shell_off[a + 1]; ++q) {
                            const uint32_t s = ctx->atom_shell[q];
                            const mbis_shell_t* sh = &ctx->shells[s];
                            const double rc = sqrt(sh->r2_cut) + rsurv;
                            const double d[3] = {sh->c[0] - cc[0], sh->c[1] - cc[1], sh->c[2] - cc[2]};
                            if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > rc * rc) continue;
                            t->shells[ns] = s;
                            t->loc[ns] = (uint32_t)nl;
                            ++ns;
                            nl += sh->num_ao;
                        }
                    }
                }
                if (ns == 0) continue;
                if (!mbis_thread_reserve(t, nl, ctx->sphere[0].max_patch)) { ctx->oom = true; return; }

                // The AOs at every one of those points, shell by shell: a primitive's exponential at all of them in
                // one vector pass (the points padded to a multiple of 4 with copies of the last). Point i's values
                // are phib[i * ls + block column], ls = nl rounded up to 4 (the padding 0); reach[i * ns + a] says
                // whether shell a is evaluated there.
                const size_t ls = (nl + 3) & ~(size_t)3;
                const uint32_t np = (nsurv + 3) & ~3u;
                for (uint32_t i = 0; i < np; ++i) {
                    const double* u = sph->u + 3 * t->surv_j[MD_MIN(i, nsurv - 1)];
                    t->bx[i] = cb[0] + r * u[0];
                    t->by[i] = cb[1] + r * u[1];
                    t->bz[i] = cb[2] + r * u[2];
                }
                for (uint32_t a = 0; a < ns; ++a) {
                    const mbis_shell_t* sh = &ctx->shells[t->shells[a]];
                    bool any = false;
                    for (uint32_t i = 0; i < np; ++i) {
                        const double dx = t->bx[i] - sh->c[0], dy = t->by[i] - sh->c[1], dz = t->bz[i] - sh->c[2];
                        t->br2[i] = dx * dx + dy * dy + dz * dz;
                    }
                    for (uint32_t i = 0; i < nsurv; ++i) {
                        const bool in = t->br2[i] <= sh->r2_cut;
                        t->reach[i * ns + a] = in;
                        any |= in;
                        if (!in) memset(t->phib + (size_t)i * ls + t->loc[a], 0, sizeof(double) * sh->num_ao);
                    }
                    if (!any) continue;
                    for (uint32_t i = 0; i < np; ++i) t->bR[i] = 0.0;
                    for (uint32_t q = sh->prim_off; q < sh->prim_off + sh->num_prim; ++q) {
                        bool near = false;
                        for (uint32_t i = 0; i < nsurv; ++i) near |= t->br2[i] < ctx->prim_r2[q];
                        if (!near) continue;   // a tight primitive, far from every point
                        for (uint32_t i = 0; i < np; ++i) t->barg[i] = -ctx->alpha[q] * t->br2[i];
                        mbis_exp(t->bexp, t->barg, np);
                        for (uint32_t i = 0; i < np; ++i) t->bR[i] += ctx->coeff[q] * t->bexp[i];
                    }
                    const mbis_cart_t* ctab = &ctx->cart[sh->l];
                    for (uint32_t i = 0; i < nsurv; ++i) {
                        if (!t->reach[i * ns + a]) continue;
                        const double dd[3] = {t->bx[i] - sh->c[0], t->by[i] - sh->c[1], t->bz[i] - sh->c[2]};
                        double px[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], py[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], pz[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
                        px[0] = py[0] = pz[0] = 1.0;
                        for (uint32_t e = 1; e <= sh->l; ++e) { px[e] = px[e - 1] * dd[0]; py[e] = py[e - 1] * dd[1]; pz[e] = pz[e - 1] * dd[2]; }
                        double* dst = t->phib + (size_t)i * ls + t->loc[a];
                        for (uint32_t n = 0; n < sh->num_ao; ++n) {
                            dst[n] = ctab->nrm[n] * px[ctab->ijk[n][0]] * py[ctab->ijk[n][1]] * pz[ctab->ijk[n][2]] * t->bR[i];
                        }
                    }
                }

                // rho_p = sum_k phi_pk (D_kk phi_pk + 2 sum_{l > k} D_kl phi_pl): the upper triangle. A row of D is
                // gathered once per batch over the block's columns (which come in runs of consecutive AOs: shells in
                // AO order, a shell continuing the last run when its AOs follow on) into one contiguous row, then
                // dotted with the points where the row's AO is not 0, four at a time: each element of the row is
                // loaded once for four points, and nothing is stored but the four sums.
                uint32_t nrun = 0;
                for (uint32_t a = 0; a < ns; ++a) {
                    const mbis_shell_t* sh = &ctx->shells[t->shells[a]];
                    if (nrun > 0 && t->run_col[nrun - 1] + t->run_len[nrun - 1] == sh->ao) {
                        t->run_len[nrun - 1] += sh->num_ao;
                    } else {
                        t->run_beg[nrun] = t->loc[a];
                        t->run_col[nrun] = sh->ao;
                        t->run_len[nrun] = sh->num_ao;
                        ++nrun;
                    }
                }
                for (uint32_t i = 0; i < nsurv; ++i) {
                    memset(t->phib + (size_t)i * ls + nl, 0, sizeof(double) * (ls - nl));
                    t->brho[i] = 0.0;
                }
                memset(t->phib + (size_t)nsurv * ls, 0, sizeof(double) * 4);   // a vector load past the last row
                double* dl = t->drow;
                for (uint32_t ri = 0; ri < nrun; ++ri) {
                    for (uint32_t o = 0; o < t->run_len[ri]; ++o) {
                        const uint32_t lk = t->run_beg[ri] + o;
                        const uint32_t gk = t->run_col[ri] + o;
                        uint32_t nact = 0;
                        for (uint32_t si = 0; si < nsurv; ++si) {
                            if (t->phib[(size_t)si * ls + lk] != 0.0) t->act_p[nact++] = si;
                        }
                        if (nact == 0) continue;
                        // Row gk of D over the block's columns from lk on, 0 past the last
                        const double* drow = D + (size_t)gk * N;
                        for (uint32_t q = o; q < t->run_len[ri]; ++q) dl[t->run_beg[ri] + q] = drow[t->run_col[ri] + q];
                        for (uint32_t rj = ri + 1; rj < nrun; ++rj) {
                            const double* src = drow + t->run_col[rj];
                            double* dst = dl + t->run_beg[rj];
                            for (uint32_t q = 0; q < t->run_len[rj]; ++q) dst[q] = src[q];
                        }
                        dl[nl] = dl[nl + 1] = dl[nl + 2] = dl[nl + 3] = 0.0;
                        // Four points at a time: sum_{l > k} D_kl phi_pl as one vector pass over the row
                        for (uint32_t g = 0; g < nact; g += 4) {
                            const double* ph[4];
                            for (int q = 0; q < 4; ++q) ph[q] = g + q < nact ? t->phib + (size_t)t->act_p[g + q] * ls : t->zrow;
                            md_256d acc0 = md_mm256_setzero_pd(), acc1 = md_mm256_setzero_pd(), acc2 = md_mm256_setzero_pd(), acc3 = md_mm256_setzero_pd();
                            for (size_t c = lk + 1; c < nl; c += 4) {
                                const md_256d d = md_mm256_loadu_pd(dl + c);
                                acc0 = md_mm256_fmadd_pd(md_mm256_loadu_pd(ph[0] + c), d, acc0);
                                acc1 = md_mm256_fmadd_pd(md_mm256_loadu_pd(ph[1] + c), d, acc1);
                                acc2 = md_mm256_fmadd_pd(md_mm256_loadu_pd(ph[2] + c), d, acc2);
                                acc3 = md_mm256_fmadd_pd(md_mm256_loadu_pd(ph[3] + c), d, acc3);
                            }
                            double sum[4][4];
                            md_mm256_storeu_pd(sum[0], acc0);
                            md_mm256_storeu_pd(sum[1], acc1);
                            md_mm256_storeu_pd(sum[2], acc2);
                            md_mm256_storeu_pd(sum[3], acc3);
                            for (uint32_t q = 0; q < 4 && g + q < nact; ++q) {
                                const double f = ph[q][lk];
                                const double v = (sum[q][0] + sum[q][1]) + (sum[q][2] + sum[q][3]);
                                t->brho[t->act_p[g + q]] += f * (dl[lk] * f + 2.0 * v);
                            }
                        }
                    }
                }

                for (uint32_t si = 0; si < nsurv; ++si) {
                    const uint32_t j = t->surv_j[si];
                    const double wb = t->surv_w[si];
                    const double rho = t->brho[si];
                    if (!(rho > MBIS_RHO_FLOOR)) continue;

                    const double mass = rw[ir] * sph->w[j] * wb * rho;
                    if (!(mass > 0.0)) continue;
                    out->rad[out->count]  = (uint16_t)ir;
                    out->ang[out->count]  = (uint16_t)j;
                    out->mass[out->count] = mass;
                    out->count++;
                    out->electrons += mass;
                }
            }
        }
        // Allocated for every point of the grid, holding those kept: give the rest back
        if (out->count > 0) {
            uint16_t* rad = (uint16_t*)realloc(out->rad, sizeof(uint16_t) * out->count);
            uint16_t* ang = (uint16_t*)realloc(out->ang, sizeof(uint16_t) * out->count);
            double* mass = (double*)realloc(out->mass, sizeof(double) * out->count);
            if (rad) out->rad = rad;
            if (ang) out->ang = ang;
            if (mass) out->mass = mass;
        }
        const MD_U64 done = md_atomic_fetch_add64(&ctx->atoms_done, 1) + 1;
        if (ctx->desc->progress && thread_idx == 0) *ctx->desc->progress = 0.5f * (float)done / (float)ctx->num_atoms;
    }
}

// The pro-molecule at every kept point of the atoms of chunk c, and each shell's share of the density.
static void mbis_iter_task(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    mbis_ctx_t* ctx = (mbis_ctx_t*)arg;
    mbis_thread_t* t = &ctx->thr[thread_idx];
    const size_t na = ctx->num_atoms;

    for (uint32_t c = beg; c < end; ++c) {
        double* accN = ctx->acc + (size_t)c * 2 * ctx->num_pro;
        double* accM = accN + ctx->num_pro;
        memset(accN, 0, sizeof(double) * 2 * ctx->num_pro);
        const uint32_t b0 = (uint32_t)((uint64_t)na * c / ctx->num_chunks);
        const uint32_t b1 = (uint32_t)((uint64_t)na * (c + 1) / ctx->num_chunks);

        for (uint32_t b = b0; b < b1; ++b) {
            const mbis_atom_points_t* pts = &ctx->pts[b];
            const double* cb = ctx->xyz[b].elem;
            const double* rr = ctx->rad_r + (size_t)b * ctx->num_rad;
            const mbis_neighbor_t* nb = ctx->nb + (size_t)b * na;
            uint32_t cur_rad = UINT32_MAX;
            size_t nsub = 0;
            double floor_b = 0.0;
            const mbis_sphere_t* sph = NULL;

            for (uint32_t ip = 0; ip < pts->count; ++ip) {
                const uint32_t ir = pts->rad[ip];
                // The coarse grid: the odd radial points, half of them and a valid radial rule of its own (the
                // Mura-Knowles nodes x = i / (n + 1) at every other i), with twice the weight
                if (ctx->coarse && (ir % MBIS_COARSE_STRIDE) != MBIS_COARSE_STRIDE - 1) continue;
                const double r = rr[ir];
                if (ir != cur_rad) {
                    sph = mbis_sphere_at(ctx, b, r);
                    // b's own pro-atom is the same all over the sphere of radius r, and bounds the
                    // pro-molecule from below there: atoms whose every shell stays below a
                    // MBIS_TERM_EPS fraction of it are left out of this radial shell.
                    cur_rad = ir;
                    double rho_b = 0.0;
                    for (uint32_t k = ctx->shell_off[b]; k < ctx->shell_off[b + 1]; ++k) rho_b += ctx->pre[k] * exp(-ctx->S[k] * r);
                    floor_b = MBIS_TERM_EPS * rho_b;
                    const double lf = floor_b > 0.0 ? log(floor_b) : -HUGE_VAL;
                    // Every shell of atom a stays below exp(logpre_a - smin_a d) at distance d: beyond
                    // rc_a = (logpre_a - ln floor) / smin_a none of them counts. Down b's neighbour list the
                    // walk stops past the longest reach of any atom.
                    nsub = 0;
                    for (size_t i = 0; i < na && nb[i].d - r <= ctx->reach_bound * (ctx->logpre_max - lf); ++i) {
                        const uint32_t a = nb[i].atom;
                        const double rc = (ctx->atom_logpre[a] - lf) / ctx->atom_smin[a];
                        if (!(rc >= fabs(nb[i].d - r))) continue;   // |d - r|: the closest a gets to the sphere
                        t->sub[nsub] = a;
                        t->sub_rc2[nsub] = rc * rc;
                        ++nsub;
                    }
                }

                const double* u = sph->u + 3 * (size_t)pts->ang[ip];
                const double p[3] = {cb[0] + r * u[0], cb[1] + r * u[1], cb[2] + r * u[2]};

                size_t nt = 0;
                double rho0 = 0.0;
                for (size_t i = 0; i < nsub; ++i) {
                    const uint32_t a = t->sub[i];
                    double ra = r;
                    if (a != b) {
                        const double* ca = ctx->xyz[a].elem;
                        const double d[3] = {p[0] - ca[0], p[1] - ca[1], p[2] - ca[2]};
                        const double d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        if (d2 > t->sub_rc2[i]) continue;
                        ra = sqrt(d2);
                    }
                    for (uint32_t k = ctx->shell_off[a]; k < ctx->shell_off[a + 1]; ++k) {
                        t->term_k[nt] = k;
                        t->term_r[nt] = ra;
                        t->term_x[nt] = ctx->logpre[k] - ctx->S[k] * ra;
                        ++nt;
                    }
                }
                // The terms' values in one vector pass, then those below the floor out
                for (size_t q = nt; q < ((nt + 3) & ~(size_t)3); ++q) t->term_x[q] = -708.0;
                mbis_exp(t->term_v, t->term_x, (nt + 3) & ~(size_t)3);
                for (size_t q = 0; q < nt; ++q) {
                    if (t->term_v[q] < floor_b) t->term_v[q] = 0.0;
                    rho0 += t->term_v[q];
                }
                if (!(rho0 > 0.0)) continue;
                const double scale = (ctx->coarse ? (double)MBIS_COARSE_STRIDE : 1.0) * pts->mass[ip] / rho0;
                for (size_t q = 0; q < nt; ++q) {
                    const double f = scale * t->term_v[q];
                    accN[t->term_k[q]] += f;
                    accM[t->term_k[q]] += f * t->term_r[q];
                }
            }
        }
    }
}

static uint32_t mbis_num_shells(int electrons) {
    static const int noble[] = {2, 10, 18, 36, 54, 86, 118};
    uint32_t n = 1;
    for (size_t i = 0; i < MD_COUNTOF(noble) && electrons > noble[i]; ++i) ++n;
    return n;
}

static void mbis_ctx_free(mbis_ctx_t* ctx, uint32_t num_threads) {
    free(ctx->shells); free(ctx->atom_shell_off); free(ctx->atom_shell); free(ctx->alpha); free(ctx->prim_r2); free(ctx->coeff);
    free(ctx->rad_r); free(ctx->rad_w); free(ctx->inner_r);
    for (int i = 0; i < 2; ++i) { free(ctx->sphere[i].u); free(ctx->sphere[i].w); }
    free(ctx->nb); free(ctx->inv_dist);
    if (ctx->pts) {
        for (size_t i = 0; i < ctx->num_atoms; ++i) { free(ctx->pts[i].rad); free(ctx->pts[i].ang); free(ctx->pts[i].mass); }
        free(ctx->pts);
    }
    if (ctx->thr) {
        for (uint32_t i = 0; i < num_threads; ++i) {
            mbis_thread_t* t = &ctx->thr[i];
            free(t->zrow); free(t->brho); free(t->act_p); free(t->drow); free(t->phib); free(t->reach); free(t->bx); free(t->by); free(t->bz); free(t->br2); free(t->bR); free(t->barg); free(t->bexp); free(t->shells); free(t->loc); free(t->run_beg); free(t->run_col); free(t->run_len); free(t->cand); free(t->dc.r); free(t->dc.stamp); free(t->cand_r);
            free(t->sub); free(t->sub_rc2); free(t->surv_j); free(t->surv_w); free(t->term_k); free(t->term_v); free(t->term_x); free(t->term_r);
        }
        free(ctx->thr);
    }
    free(ctx->shell_off); free(ctx->N); free(ctx->S); free(ctx->pre); free(ctx->logpre);
    free(ctx->atom_logpre); free(ctx->atom_smin); free(ctx->acc);
}

bool md_mbis_compute(md_mbis_result_t* out, const md_mbis_desc_t* desc, struct md_allocator_i* alloc) {
    if (!out || !desc || !alloc || !desc->basis || !desc->atom_coord || !desc->atomic_number || !desc->density_matrix || desc->num_atoms == 0) {
        MD_LOG_ERROR("md_mbis_compute: missing argument");
        return false;
    }
    memset(out, 0, sizeof(*out));
    const md_gto_basis_t* basis = desc->basis;
    if (md_gto_basis_num_atoms(basis) > desc->num_atoms) {
        MD_LOG_ERROR("md_mbis_compute: the basis indexes %zu atoms, %zu were given", md_gto_basis_num_atoms(basis), desc->num_atoms);
        return false;
    }
    const uint32_t num_threads = md_executor_num_threads(desc->executor);

    mbis_ctx_t ctx = {0};
    ctx.desc = desc;
    ctx.xyz = desc->atom_coord;
    ctx.num_atoms = desc->num_atoms;
    ctx.num_ao = md_gto_basis_num_ao(basis);
    ctx.num_shells = basis->num_shells;
    ctx.num_rad = desc->radial_points ? desc->radial_points : MBIS_DEFAULT_RADIAL;
    const uint32_t order = desc->angular_order ? desc->angular_order : MBIS_DEFAULT_ANGULAR;
    const double threshold = desc->threshold > 0.0 ? desc->threshold : MBIS_DEFAULT_THRESHOLD;
    const uint32_t max_iter = desc->max_iterations ? desc->max_iterations : MBIS_DEFAULT_MAX_ITER;
    if (ctx.num_rad > UINT16_MAX || order > 2 * 255) {
        MD_LOG_ERROR("md_mbis_compute: grid too large");
        return false;
    }
    const size_t na = ctx.num_atoms;
    bool ok = false;

    // The basis
    ctx.shells = (mbis_shell_t*)calloc(ctx.num_shells ? ctx.num_shells : 1, sizeof(mbis_shell_t));
    ctx.alpha  = (double*)calloc(basis->num_primitives ? basis->num_primitives : 1, sizeof(double));
    ctx.prim_r2 = (double*)calloc(basis->num_primitives ? basis->num_primitives : 1, sizeof(double));
    ctx.coeff  = (double*)calloc(basis->num_primitives ? basis->num_primitives : 1, sizeof(double));
    if (!ctx.shells || !ctx.alpha || !ctx.prim_r2 || !ctx.coeff) goto oom;
    for (uint32_t p = 0; p < basis->num_primitives; ++p) { ctx.alpha[p] = basis->alpha[p]; ctx.coeff[p] = basis->coeff[p]; }
    for (uint32_t l = 0; l <= MD_GTO_MAX_ANGULAR_MOMENTUM; ++l) {
        for (uint32_t n = 0; n < md_gto_num_cart_ao(l); ++n) {
            int i, j, k;
            md_gto_cart_ijk(&i, &j, &k, l, n);
            ctx.cart[l].ijk[n][0] = (uint8_t)i; ctx.cart[l].ijk[n][1] = (uint8_t)j; ctx.cart[l].ijk[n][2] = (uint8_t)k;
            ctx.cart[l].nrm[n] = md_gto_cart_norm_factor(i, j, k);
        }
    }
    {
        uint32_t ao = 0;
        for (uint32_t s = 0; s < ctx.num_shells; ++s) {
            const md_gto_shell_t* bs = &basis->shells[s];
            if (bs->l > MD_GTO_MAX_ANGULAR_MOMENTUM) {
                MD_LOG_ERROR("md_mbis_compute: angular momentum %u is not supported", bs->l);
                goto fail;
            }
            mbis_shell_t* sh = &ctx.shells[s];
            sh->l = bs->l;
            sh->prim_off = bs->primitive_offset;
            sh->num_prim = bs->num_primitives;
            sh->ao = ao;
            sh->num_ao = md_gto_num_cart_ao(bs->l);
            ao += sh->num_ao;
            memcpy(sh->c, ctx.xyz[bs->atom_idx].elem, sizeof(sh->c));
            const double rc = mbis_shell_cutoff(ctx.alpha + sh->prim_off, ctx.coeff + sh->prim_off, sh->num_prim, sh->l, MBIS_AO_EPS);
            sh->r2_cut = rc * rc;
            ctx.shell_reach = MD_MAX(ctx.shell_reach, rc);
            // Each primitive to a share of the same bound, so the shell's stays below it
            for (uint32_t q = sh->prim_off; q < sh->prim_off + sh->num_prim; ++q) {
                const double rq = mbis_shell_cutoff(ctx.alpha + q, ctx.coeff + q, 1, sh->l, MBIS_AO_EPS / sh->num_prim);
                ctx.prim_r2[q] = rq * rq;
            }
        }
    }

    // Every atom's shells, in shell order
    ctx.atom_shell_off = (uint32_t*)calloc(na + 1, sizeof(uint32_t));
    ctx.atom_shell = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(ctx.num_shells, 1u));
    if (!ctx.atom_shell_off || !ctx.atom_shell) goto oom;
    for (uint32_t s = 0; s < ctx.num_shells; ++s) ctx.atom_shell_off[basis->shells[s].atom_idx + 1]++;
    for (size_t a = 0; a < na; ++a) ctx.atom_shell_off[a + 1] += ctx.atom_shell_off[a];
    {
        uint32_t* fill = (uint32_t*)malloc(sizeof(uint32_t) * (na + 1));
        if (!fill) goto oom;
        memcpy(fill, ctx.atom_shell_off, sizeof(uint32_t) * (na + 1));
        for (uint32_t s = 0; s < ctx.num_shells; ++s) ctx.atom_shell[fill[basis->shells[s].atom_idx]++] = s;
        free(fill);
    }

    // The grids
    ctx.rad_r = (double*)malloc(sizeof(double) * na * ctx.num_rad);
    ctx.rad_w = (double*)malloc(sizeof(double) * na * ctx.num_rad);
    ctx.inner_r = (double*)malloc(sizeof(double) * na);
    if (!ctx.rad_r || !ctx.rad_w || !ctx.inner_r) goto oom;
    for (size_t a = 0; a < na; ++a) mbis_radial_grid(ctx.rad_r + a * ctx.num_rad, ctx.rad_w + a * ctx.num_rad, ctx.num_rad, desc->atomic_number[a]);
    if (!mbis_sphere_init(&ctx.sphere[0], order) || !mbis_sphere_init(&ctx.sphere[1], MD_MIN(order, MBIS_INNER_ANGULAR))) goto oom;

    // Every atom's neighbours, nearest first (itself first)
    ctx.nb = (mbis_neighbor_t*)malloc(sizeof(mbis_neighbor_t) * na * na);
    ctx.inv_dist = (double*)malloc(sizeof(double) * na * na);
    if (!ctx.nb || !ctx.inv_dist) goto oom;
    for (size_t a = 0; a < na; ++a) {
        for (size_t b = 0; b < na; ++b) {
            const double d = mbis_dist(ctx.xyz[a].elem, ctx.xyz[b].elem);
            ctx.inv_dist[a * na + b] = (a != b && d > 0.0) ? 1.0 / d : 0.0;
        }
    }
    for (size_t b = 0; b < na; ++b) {
        mbis_neighbor_t* row = ctx.nb + b * na;
        for (size_t a = 0; a < na; ++a) row[a] = (mbis_neighbor_t){ (uint32_t)a, a == b ? -1.0 : mbis_dist(ctx.xyz[a].elem, ctx.xyz[b].elem) };
        qsort(row, na, sizeof(mbis_neighbor_t), mbis_neighbor_cmp);
        row[0].d = 0.0;
        ctx.inner_r[b] = MBIS_INNER_FRACTION * (na > 1 ? row[1].d : MBIS_INNER_LONE_ATOM);
    }

    // Pro-atoms: shells and the initial guess of Verstraelen et al. (as HORTON): equal populations,
    // exponents spaced geometrically from 2 n to 2.
    ctx.shell_off = (uint32_t*)malloc(sizeof(uint32_t) * (na + 1));
    if (!ctx.shell_off) goto oom;
    ctx.shell_off[0] = 0;
    for (size_t a = 0; a < na; ++a) {
        const double q = desc->nuclear_charge ? desc->nuclear_charge[a] : desc->atomic_number[a];
        ctx.shell_off[a + 1] = ctx.shell_off[a] + mbis_num_shells((int)lround(q));
    }
    ctx.num_pro = ctx.shell_off[na];
    ctx.N = (double*)malloc(sizeof(double) * ctx.num_pro);
    ctx.S = (double*)malloc(sizeof(double) * ctx.num_pro);
    ctx.pre = (double*)malloc(sizeof(double) * ctx.num_pro);
    ctx.logpre = (double*)malloc(sizeof(double) * ctx.num_pro);
    ctx.atom_logpre = (double*)malloc(sizeof(double) * na);
    ctx.atom_smin = (double*)malloc(sizeof(double) * na);
    ctx.num_chunks = (uint32_t)MD_MIN(na, (size_t)MBIS_MAX_CHUNKS);
    ctx.acc = (double*)malloc(sizeof(double) * 2 * ctx.num_pro * ctx.num_chunks);
    if (!ctx.N || !ctx.S || !ctx.pre || !ctx.logpre || !ctx.atom_logpre || !ctx.atom_smin || !ctx.acc) goto oom;
    for (size_t a = 0; a < na; ++a) {
        const double q = desc->nuclear_charge ? desc->nuclear_charge[a] : desc->atomic_number[a];
        const double el = MD_MAX(q, 0.0);
        const uint32_t ns = ctx.shell_off[a + 1] - ctx.shell_off[a];
        const double s0 = 2.0 * MD_MAX(el, 1.0);
        const double step = ns > 1 ? pow(2.0 / s0, 1.0 / (ns - 1)) : 1.0;
        for (uint32_t i = 0; i < ns; ++i) {
            ctx.N[ctx.shell_off[a] + i] = el / ns;
            ctx.S[ctx.shell_off[a] + i] = s0 * pow(step, (double)i);
        }
    }
    if (desc->initial) {
        // A previous partition of the same atoms, the last frame's say: much closer to the answer than the guess
        const md_mbis_result_t* in = desc->initial;
        bool same = in->num_atoms == na && in->num_shells == ctx.num_pro && in->shell_offset && in->shell_population && in->shell_width;
        for (size_t a = 0; same && a <= na; ++a) same = in->shell_offset[a] == ctx.shell_off[a];
        if (!same) {
            MD_LOG_ERROR("md_mbis_compute: the initial pro-atoms are not of these atoms");
            goto fail;
        }
        for (uint32_t k = 0; k < ctx.num_pro; ++k) {
            if (in->shell_population[k] >= 0.0 && in->shell_width[k] > 0.0) {
                ctx.N[k] = in->shell_population[k];
                ctx.S[k] = 1.0 / in->shell_width[k];
            }
        }
    }

    // Per thread scratch
    ctx.thr = (mbis_thread_t*)calloc(num_threads, sizeof(mbis_thread_t));
    ctx.pts = (mbis_atom_points_t*)calloc(na, sizeof(mbis_atom_points_t));
    if (!ctx.thr || !ctx.pts) goto oom;
    for (uint32_t i = 0; i < num_threads; ++i) {
        mbis_thread_t* t = &ctx.thr[i];
        t->shells    = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(ctx.num_shells, 1u));
        t->loc       = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(ctx.num_shells, 1u));
        t->run_beg   = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(ctx.num_shells, 1u));
        t->run_col   = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(ctx.num_shells, 1u));
        t->run_len   = (uint32_t*)malloc(sizeof(uint32_t) * MD_MAX(ctx.num_shells, 1u));
        const size_t bs = ctx.sphere[0].max_patch;
        t->reach = (uint8_t*)malloc(bs * MD_MAX(ctx.num_shells, 1u));
        t->brho  = (double*)malloc(sizeof(double) * bs);
        t->act_p = (uint32_t*)malloc(sizeof(uint32_t) * bs);
        if (!t->brho || !t->act_p) goto oom;
        double** barr[7] = {&t->bx, &t->by, &t->bz, &t->br2, &t->bR, &t->barg, &t->bexp};
        for (int q = 0; q < 7; ++q) { *barr[q] = (double*)malloc(sizeof(double) * (bs + 4)); if (!*barr[q]) goto oom; }
        if (!t->reach) goto oom;
        t->dc.r     = (double*)malloc(sizeof(double) * na);
        t->dc.stamp = (uint32_t*)calloc(na, sizeof(uint32_t));
        if (!t->dc.r || !t->dc.stamp) goto oom;
        t->cand   = (uint32_t*)malloc(sizeof(uint32_t) * na);
        t->cand_r = (double*)malloc(sizeof(double) * na);
        t->sub    = (uint32_t*)malloc(sizeof(uint32_t) * na);
        t->surv_j = (uint32_t*)malloc(sizeof(uint32_t) * ctx.sphere[0].max_patch);
        t->surv_w = (double*)malloc(sizeof(double) * ctx.sphere[0].max_patch);
        t->sub_rc2 = (double*)malloc(sizeof(double) * na);
        t->term_k = (uint32_t*)malloc(sizeof(uint32_t) * ctx.num_pro);
        t->term_v = (double*)malloc(sizeof(double) * (ctx.num_pro + 4));
        t->term_x = (double*)malloc(sizeof(double) * (ctx.num_pro + 4));
        t->term_r = (double*)malloc(sizeof(double) * ctx.num_pro);
        if (!t->shells || !t->loc || !t->run_beg || !t->run_col || !t->run_len || !t->cand || !t->cand_r || !t->sub || !t->sub_rc2 || !t->surv_j || !t->surv_w || !t->term_k || !t->term_v || !t->term_x || !t->term_r) goto oom;
    }

    // The density on the grid
    if (!md_parallel_for(desc->executor, (uint32_t)na, 1, mbis_grid_task, &ctx)) goto fail;
    if (ctx.oom) goto oom;
    if (desc->cancel && *desc->cancel) goto fail;
    if (desc->progress) *desc->progress = 0.5f;
    double electrons = 0.0;
    size_t num_points = 0;
    for (size_t a = 0; a < na; ++a) { electrons += ctx.pts[a].electrons; num_points += ctx.pts[a].count; }

    // The fixed point x = G(x) over x = (N, ln S), with Anderson acceleration (Walker and Ni, SIAM J. Numer. Anal.
    // 49, 1715 (2011)): the next iterate is the combination of the last MBIS_ANDERSON_DEPTH images that minimizes
    // the linearized residual. MBIS' own iteration converges linearly, slowly where an atom's shells trade
    // electrons between them; this takes ~5 times fewer steps to the same fixed point. A combination with a
    // population below zero is not taken: the plain step is, and the history restarts.
    const size_t P = ctx.num_pro, X = 2 * P;
    double* aa = (double*)malloc(sizeof(double) * X * (4 + 2 * MBIS_ANDERSON_DEPTH));
    if (!aa) goto oom;
    double* x  = aa;
    double* g  = x + X;
    double* fp = g + X;            // the previous residual
    double* gp = fp + X;           // the previous image
    double* dF = gp + X;           // [depth][X] differences of residuals
    double* dG = dF + X * MBIS_ANDERSON_DEPTH;
    uint32_t hist = 0, head = 0;
    for (size_t k = 0; k < P; ++k) { x[k] = ctx.N[k]; x[P + k] = log(ctx.S[k]); }

    uint32_t iter = 0;
    double change = HUGE_VAL, change0 = 0.0;
    bool have_prev = false;
    // Coarse to fine: the coarse grid's fixed point is within ~1e-4 of the full one's, so iterating there first
    // leaves the full grid only the last few steps
    const double coarse_threshold = MD_MAX(MBIS_COARSE_THRESHOLD, threshold);
    ctx.coarse = MBIS_COARSE_THRESHOLD > 0.0 && !desc->initial;   // a previous answer is closer than the coarse grid's
    while (iter < max_iter) {
        if (desc->cancel && *desc->cancel) { free(aa); goto fail; }
        for (size_t k = 0; k < P; ++k) { ctx.N[k] = x[k]; ctx.S[k] = exp(x[P + k]); }
        double logpre_max = -HUGE_VAL, smin_min = HUGE_VAL;
        for (size_t a = 0; a < na; ++a) {
            double lp = -HUGE_VAL, smin = HUGE_VAL;
            for (uint32_t k = ctx.shell_off[a]; k < ctx.shell_off[a + 1]; ++k) {
                ctx.pre[k] = ctx.N[k] * ctx.S[k] * ctx.S[k] * ctx.S[k] / (8.0 * 3.14159265358979323846);
                ctx.logpre[k] = ctx.pre[k] > 0.0 ? log(ctx.pre[k]) : -HUGE_VAL;
                lp = MD_MAX(lp, ctx.logpre[k]);
                smin = MD_MIN(smin, ctx.S[k]);
            }
            ctx.atom_logpre[a] = lp;
            ctx.atom_smin[a] = smin;
            logpre_max = MD_MAX(logpre_max, lp);
            smin_min = MD_MIN(smin_min, smin);
        }
        // No atom's terms count farther than (logpre_max - ln floor) / smin_min from it: the walk down a neighbour
        // list per radial sphere stops there
        ctx.logpre_max = logpre_max;
        ctx.reach_bound = 1.0 / smin_min;
        if (!md_parallel_for(desc->executor, ctx.num_chunks, 1, mbis_iter_task, &ctx)) { free(aa); goto fail; }

        // The image G(x), and the residual's largest entry: a population change in e, and the relative change
        // of a width (to first order)
        change = 0.0;
        for (size_t k = 0; k < P; ++k) {
            double n = 0.0, m = 0.0;
            for (uint32_t c = 0; c < ctx.num_chunks; ++c) {
                n += ctx.acc[(size_t)c * X + k];
                m += ctx.acc[(size_t)c * X + P + k];
            }
            g[k] = n;
            g[P + k] = (n > 0.0 && m > 0.0) ? log(3.0 * n / m) : x[P + k];
            change = MD_MAX(change, MD_MAX(fabs(g[k] - x[k]), fabs(g[P + k] - x[P + k])));
        }
        ++iter;
        if (iter == 1) change0 = change;
        if (desc->progress && change0 > threshold && change > 0.0) {
            // the residual falls about geometrically: how far it has come on a log scale
            const float f = (float)(log(change0 / change) / log(change0 / threshold));
            const float p = 0.5f + 0.5f * MD_CLAMP(f, 0.0f, 1.0f);
            if (p > *desc->progress) *desc->progress = p;
        }
        if (ctx.coarse && change < coarse_threshold) {
            // On to the full grid, from the coarse grid's answer; the history is of another map
            ctx.coarse = false;
            memcpy(x, g, sizeof(double) * X);
            have_prev = false;
            hist = 0;
            head = 0;
            continue;
        }
        if (!ctx.coarse && change < threshold) {
            memcpy(x, g, sizeof(double) * X);
            break;
        }

        // The history of differences, oldest overwritten first
        if (have_prev) {
            double* df = dF + (size_t)head * X;
            double* dg = dG + (size_t)head * X;
            for (size_t i = 0; i < X; ++i) {
                df[i] = (g[i] - x[i]) - fp[i];
                dg[i] = g[i] - gp[i];
            }
            head = (head + 1) % MBIS_ANDERSON_DEPTH;
            hist = MD_MIN(hist + 1, (uint32_t)MBIS_ANDERSON_DEPTH);
        }
        for (size_t i = 0; i < X; ++i) { fp[i] = g[i] - x[i]; gp[i] = g[i]; }
        have_prev = true;

        // gamma = argmin |f - dF gamma|, by the (regularized) normal equations; then x = g - dG gamma
        bool accelerated = false;
        if (hist > 0) {
            double A[MBIS_ANDERSON_DEPTH][MBIS_ANDERSON_DEPTH + 1];
            const uint32_t m = hist;
            double tr = 0.0;
            for (uint32_t i = 0; i < m; ++i) {
                const double* di = dF + (size_t)i * X;
                for (uint32_t j = 0; j <= i; ++j) {
                    const double* dj = dF + (size_t)j * X;
                    double v = 0.0;
                    for (size_t q = 0; q < X; ++q) v += di[q] * dj[q];
                    A[i][j] = A[j][i] = v;
                }
                double v = 0.0;
                for (size_t q = 0; q < X; ++q) v += di[q] * fp[q];
                A[i][m] = v;
                tr += A[i][i];
            }
            for (uint32_t i = 0; i < m; ++i) A[i][i] += 1.0e-12 * tr + 1.0e-300;
            // Gaussian elimination with partial pivoting
            bool singular = false;
            for (uint32_t c = 0; c < m && !singular; ++c) {
                uint32_t piv = c;
                for (uint32_t r = c + 1; r < m; ++r) if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
                if (!(fabs(A[piv][c]) > 0.0)) { singular = true; break; }
                if (piv != c) for (uint32_t q = 0; q <= m; ++q) { const double t = A[c][q]; A[c][q] = A[piv][q]; A[piv][q] = t; }
                for (uint32_t r = c + 1; r < m; ++r) {
                    const double f = A[r][c] / A[c][c];
                    for (uint32_t q = c; q <= m; ++q) A[r][q] -= f * A[c][q];
                }
            }
            if (!singular) {
                double gamma[MBIS_ANDERSON_DEPTH];
                for (int r = (int)m - 1; r >= 0; --r) {
                    double v = A[r][m];
                    for (uint32_t q = (uint32_t)r + 1; q < m; ++q) v -= A[r][q] * gamma[q];
                    gamma[r] = v / A[r][r];
                }
                accelerated = true;
                for (size_t i = 0; i < X; ++i) {
                    double v = g[i];
                    for (uint32_t j = 0; j < m; ++j) v -= gamma[j] * dG[(size_t)j * X + i];
                    x[i] = v;
                    if (!isfinite(v) || (i < P && v < 0.0)) { accelerated = false; break; }
                }
            }
        }
        if (!accelerated) {
            memcpy(x, g, sizeof(double) * X);
            hist = 0;
            head = 0;
        }
    }
    for (size_t k = 0; k < P; ++k) { ctx.N[k] = x[k]; ctx.S[k] = exp(x[P + k]); }
    free(aa);

    // The result
    out->num_atoms = na;
    out->num_shells = ctx.num_pro;
    out->charge = (double*)md_alloc(alloc, sizeof(double) * na);
    out->population = (double*)md_alloc(alloc, sizeof(double) * na);
    out->shell_offset = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * (na + 1));
    out->shell_population = (double*)md_alloc(alloc, sizeof(double) * ctx.num_pro);
    out->shell_width = (double*)md_alloc(alloc, sizeof(double) * ctx.num_pro);
    if (!out->charge || !out->population || !out->shell_offset || !out->shell_population || !out->shell_width) {
        md_mbis_result_free(out, alloc);
        goto oom;
    }
    memcpy(out->shell_offset, ctx.shell_off, sizeof(uint32_t) * (na + 1));
    for (uint32_t k = 0; k < ctx.num_pro; ++k) {
        out->shell_population[k] = ctx.N[k];
        out->shell_width[k] = 1.0 / ctx.S[k];
    }
    for (size_t a = 0; a < na; ++a) {
        double pop = 0.0;
        for (uint32_t k = ctx.shell_off[a]; k < ctx.shell_off[a + 1]; ++k) pop += ctx.N[k];
        const double q = desc->nuclear_charge ? desc->nuclear_charge[a] : desc->atomic_number[a];
        out->population[a] = pop;
        out->charge[a] = q - pop;
    }
    out->electron_count = electrons;
    out->num_points = (uint32_t)num_points;
    out->iterations = iter;
    out->change = change;
    out->converged = change < threshold;
    if (desc->progress) *desc->progress = 1.0f;
    ok = true;
    goto done;

oom:
    MD_LOG_ERROR("md_mbis_compute: out of memory");
fail:
    memset(out, 0, sizeof(*out));
done:
    mbis_ctx_free(&ctx, num_threads);
    return ok;
}

void md_mbis_result_free(md_mbis_result_t* r, struct md_allocator_i* alloc) {
    if (!r || !alloc) return;
    const size_t na = r->num_atoms;
    if (r->charge)           md_free(alloc, r->charge, sizeof(double) * na);
    if (r->population)       md_free(alloc, r->population, sizeof(double) * na);
    if (r->shell_offset)     md_free(alloc, r->shell_offset, sizeof(uint32_t) * (na + 1));
    if (r->shell_population) md_free(alloc, r->shell_population, sizeof(double) * r->num_shells);
    if (r->shell_width)      md_free(alloc, r->shell_width, sizeof(double) * r->num_shells);
    memset(r, 0, sizeof(*r));
}
