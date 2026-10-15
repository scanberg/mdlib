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

// ---------------------------------------------------------------------------------------------
// Intervals (no outward rounding; a relative slack is added where enclosures are finalized)
// ---------------------------------------------------------------------------------------------

typedef struct iv_t { double lo, hi; } iv_t;

static inline iv_t   iv(double lo, double hi)   { return (iv_t){lo, hi}; }
static inline iv_t   iv_pt(double x)            { return (iv_t){x, x}; }
static inline iv_t   iv_add(iv_t a, iv_t b)     { return (iv_t){a.lo + b.lo, a.hi + b.hi}; }
static inline iv_t   iv_neg(iv_t a)             { return (iv_t){-a.hi, -a.lo}; }
static inline double iv_mag(iv_t a)             { return fmax(fabs(a.lo), fabs(a.hi)); }
static inline double iv_width(iv_t a)           { return a.hi - a.lo; }
static inline iv_t   iv_widen(iv_t a, double e) { return (iv_t){a.lo - e, a.hi + e}; }
static inline bool   iv_contains(iv_t a, double x) { return a.lo <= x && x <= a.hi; }

static inline iv_t iv_scale(iv_t a, double s) {
    return s >= 0 ? (iv_t){a.lo * s, a.hi * s} : (iv_t){a.hi * s, a.lo * s};
}

static inline iv_t iv_mul(iv_t a, iv_t b) {
    const double p0 = a.lo * b.lo, p1 = a.lo * b.hi, p2 = a.hi * b.lo, p3 = a.hi * b.hi;
    return (iv_t){fmin(fmin(p0, p1), fmin(p2, p3)), fmax(fmax(p0, p1), fmax(p2, p3))};
}

static inline iv_t iv_hull(iv_t a, double x)  { return (iv_t){fmin(a.lo, x), fmax(a.hi, x)}; }
static inline iv_t iv_hull2(iv_t a, iv_t b)   { return (iv_t){fmin(a.lo, b.lo), fmax(a.hi, b.hi)}; }
static inline iv_t iv_isect(iv_t a, iv_t b)   { return (iv_t){fmax(a.lo, b.lo), fmin(a.hi, b.hi)}; }

// Exact range of x^n over [a.lo, a.hi]
static inline iv_t iv_pow(iv_t a, int n) {
    if (n == 0) return iv_pt(1.0);
    if (n == 1) return a;
    const double pl = pow(a.lo, n), ph = pow(a.hi, n);
    if (n & 1) return (iv_t){pl, ph};
    if (a.lo >= 0) return (iv_t){pl, ph};
    if (a.hi <= 0) return (iv_t){ph, pl};
    return (iv_t){0.0, fmax(pl, ph)};
}

// ---------------------------------------------------------------------------------------------
// Small utilities: rng, u64 -> u32 hash map, dynamic arrays via md_array
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

typedef struct iso_map_t {
    uint64_t* keys;   // stored as key+1, 0 = empty
    uint32_t* vals;
    uint32_t  cap;
    uint32_t  count;
} iso_map_t;

static inline uint64_t iso_hash64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33;
    return x;
}

static void iso_map_init(iso_map_t* m, uint32_t expected) {
    uint32_t cap = 64;
    while (cap < expected * 2) cap <<= 1;
    m->keys  = (uint64_t*)calloc(cap, sizeof(uint64_t));
    m->vals  = (uint32_t*)calloc(cap, sizeof(uint32_t));
    m->cap   = cap;
    m->count = 0;
}

static void iso_map_free(iso_map_t* m) { free(m->keys); free(m->vals); memset(m, 0, sizeof(*m)); }

static void iso_map_put(iso_map_t* m, uint64_t key, uint32_t val);

static void iso_map_grow(iso_map_t* m) {
    iso_map_t n;
    iso_map_init(&n, m->cap);
    for (uint32_t i = 0; i < m->cap; ++i) {
        if (m->keys[i]) iso_map_put(&n, m->keys[i] - 1, m->vals[i]);
    }
    iso_map_free(m);
    *m = n;
}

static void iso_map_put(iso_map_t* m, uint64_t key, uint32_t val) {
    if ((m->count + 1) * 2 > m->cap) iso_map_grow(m);
    uint32_t mask = m->cap - 1;
    uint32_t i = (uint32_t)iso_hash64(key) & mask;
    while (m->keys[i] && m->keys[i] != key + 1) i = (i + 1) & mask;
    if (!m->keys[i]) { m->keys[i] = key + 1; m->count++; }
    m->vals[i] = val;
}

static bool iso_map_get(const iso_map_t* m, uint64_t key, uint32_t* out) {
    if (!m->cap) return false;
    uint32_t mask = m->cap - 1;
    uint32_t i = (uint32_t)iso_hash64(key) & mask;
    while (m->keys[i]) {
        if (m->keys[i] == key + 1) { *out = m->vals[i]; return true; }
        i = (i + 1) & mask;
    }
    return false;
}

static inline uint64_t iso_key3(int x, int y, int z) {
    return ((uint64_t)(uint32_t)x) | ((uint64_t)(uint32_t)y << 21) | ((uint64_t)(uint32_t)z << 42);
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

static inline double iso_seconds_since(md_tick_t t0) { return md_tick_to_seconds(md_tick_now() - t0); }

// ---------------------------------------------------------------------------------------------
// The field: one MO, grouped by (atom, l)
// ---------------------------------------------------------------------------------------------
//
// Within one atom, every shell of angular momentum l shares the same Cartesian monomials, so the
// group sums to
//
//     psi_g(d) = sum_n  m_n(d) * R_n(|d|^2),   m_n = x^i y^j z^k,   R_n(rho) = sum_q w_nq exp(-a_q rho)
//
// with one radial function per Cartesian component n over the union of the group's exponents.
// This is where an MO cancels: 1s and 2s with opposite coefficients give R_n a radial node, and
// R_n is 1D, so its range over a rho interval is exact from its critical points.
// (Contractions on their own rarely cancel: segmented sets like def2 keep them single signed.)

#define ISO_MAX_PRIM 32
#define ISO_MAX_CART 15
#define ISO_NUM_CRIT 3   // critical intervals precomputed for g_0 (R), g_1 (R'), g_2 (R'')

// One critical interval of g_m: an interval in rho = |d|^2 that may contain a zero of g_{m+1},
// with a rigorous enclosure of g_m over it. Between critical intervals g_m is monotone.
typedef struct iso_crit_t {
    double lo, hi;
    iv_t   g;
} iso_crit_t;

typedef struct iso_shell_t {     // an (atom, l) group
    double  c[3];
    int     l, ncart, nprim;
    double  alpha[ISO_MAX_PRIM];
    double  w[ISO_MAX_CART][ISO_MAX_PRIM];   // MO coeff * f(i,j,k) * contraction coeff
    int8_t  ijk[ISO_MAX_CART][3];
    bool    live[ISO_MAX_CART];              // component has any nonzero weight
    double  r2_cut;                          // |group| < ISO_EVAL_EPS beyond this squared radius
    int     crit_off[ISO_MAX_CART][ISO_NUM_CRIT];
    int     crit_num[ISO_MAX_CART][ISO_NUM_CRIT];
} iso_shell_t;

typedef struct iso_field_t {
    iso_shell_t* shells;     // md_array
    iso_crit_t*  crit;       // md_array, shared storage
    int          num_shells;
    int          num_input_shells;
    double       atom_min[3], atom_max[3];
    char         name[64];
} iso_field_t;

#define ISO_EVAL_EPS 1.0e-12

static inline double iso_gw(const iso_shell_t* s, int q, int n, int m) {
    const double a = s->alpha[q], w = s->w[n][q];
    return (m == 0) ? w : (m == 1) ? -a * w : (m == 2) ? a * a * w : -a * a * a * w;
}

// g_m(rho) = d^m/drho^m R_n(rho)
static inline double iso_g(const iso_shell_t* s, int n, int m, double rho) {
    double r = 0.0;
    for (int q = 0; q < s->nprim; ++q) r += iso_gw(s, q, n, m) * exp(-s->alpha[q] * rho);
    return r;
}

// Naive enclosure of g_m over [u, v]: every term is monotone in rho, so each is exact on its own.
static inline iv_t iso_g_naive(const iso_shell_t* s, int n, int m, double u, double v) {
    iv_t r = iv_pt(0.0);
    for (int q = 0; q < s->nprim; ++q) {
        const double w = iso_gw(s, q, n, m);
        const double eu = w * exp(-s->alpha[q] * u), ev = w * exp(-s->alpha[q] * v);
        r = iv_add(r, iv(fmin(eu, ev), fmax(eu, ev)));
    }
    return r;
}

static void iso_crit_search(iso_field_t* f, const iso_shell_t* s, int n, int m, int first, double u, double v, int depth) {
    const iv_t d = iso_g_naive(s, n, m + 1, u, v);
    if (d.lo > 0.0 || d.hi < 0.0) return;   // g_{m+1} has no zero, g_m is monotone here
    if (depth >= 70 || (v - u) <= 1.0e-13 * (1.0 + v)) {
        iso_crit_t c = {u, v, iso_g_naive(s, n, m, u, v)};
        if ((int)md_array_size(f->crit) > first && md_array_last(f->crit)->hi >= u) {
            iso_crit_t* last = md_array_last(f->crit);
            last->hi = v;
            last->g  = iv_hull2(last->g, c.g);
        } else {
            md_array_push(f->crit, c, md_get_heap_allocator());
        }
        return;
    }
    const double mid = 0.5 * (u + v);
    iso_crit_search(f, s, n, m, first, u, mid, depth + 1);
    iso_crit_search(f, s, n, m, first, mid, v, depth + 1);
}

// Beyond rho_end the smallest exponent term of g_{m+1} dominates the sum of all the others, so
// g_{m+1} cannot change sign there. The ratio others/dominant only decreases with rho.
static double iso_crit_rho_end(const iso_shell_t* s, int n, int m1) {
    int qmin = -1;
    for (int q = 0; q < s->nprim; ++q) if (s->w[n][q] != 0.0 && (qmin < 0 || s->alpha[q] < s->alpha[qmin])) qmin = q;
    if (qmin < 0) return 0.0;
    const double amin = s->alpha[qmin];
    const double wmin = fabs(iso_gw(s, qmin, n, m1));
    for (double rho = 1.0; rho < 1.0e7; rho *= 2.0) {
        double others = 0.0;
        for (int q = 0; q < s->nprim; ++q) {
            if (q == qmin) continue;
            others += fabs(iso_gw(s, q, n, m1)) * exp(-(s->alpha[q] - amin) * rho);
        }
        if (others < 0.5 * wmin) return rho;
    }
    return 1.0e7;
}

// Range of g_m over [u, v], exact up to the width of a critical interval.
static inline iv_t iso_g_range(const iso_field_t* f, const iso_shell_t* s, int n, int m, double u, double v) {
    iv_t r = iv_pt(iso_g(s, n, m, u));
    r = iv_hull(r, iso_g(s, n, m, v));
    const iso_crit_t* c = f->crit + s->crit_off[n][m];
    for (int i = 0; i < s->crit_num[n][m]; ++i) {
        if (c[i].hi >= u && c[i].lo <= v) r = iv_hull2(r, c[i].g);
    }
    return r;
}

// Tail bound for |group| at distance >= r (valid once r is past every radial maximum).
static double iso_shell_tail(const iso_shell_t* s, double r) {
    double sum = 0.0;
    for (int n = 0; n < s->ncart; ++n) {
        double rad = 0.0;
        for (int q = 0; q < s->nprim; ++q) rad += fabs(s->w[n][q]) * exp(-s->alpha[q] * r * r);
        sum += rad;
    }
    return pow(r, s->l) * sum;
}

static double iso_shell_tail_r0(const iso_shell_t* s) {
    double amin = DBL_MAX;
    for (int q = 0; q < s->nprim; ++q) amin = fmin(amin, s->alpha[q]);
    return sqrt(s->l / (2.0 * amin));
}

static void iso_field_free(iso_field_t* f) {
    md_array_free(f->shells, md_get_heap_allocator());
    md_array_free(f->crit, md_get_heap_allocator());
    memset(f, 0, sizeof(*f));
}

// Builds the field from a basis, positions (bohr) and one Cartesian AO coefficient vector.
static bool iso_field_init(iso_field_t* f, const md_gto_basis_t* basis, const float* atom_xyz, size_t num_atoms, const double* ao_coeffs, const char* name) {
    memset(f, 0, sizeof(*f));
    snprintf(f->name, sizeof(f->name), "%s", name);
    md_allocator_i* heap = md_get_heap_allocator();
    f->num_input_shells = (int)basis->num_shells;

    for (int k = 0; k < 3; ++k) { f->atom_min[k] = DBL_MAX; f->atom_max[k] = -DBL_MAX; }
    for (size_t a = 0; a < num_atoms; ++a) {
        for (int k = 0; k < 3; ++k) {
            f->atom_min[k] = fmin(f->atom_min[k], atom_xyz[3 * a + k]);
            f->atom_max[k] = fmax(f->atom_max[k], atom_xyz[3 * a + k]);
        }
    }

    // group index per (atom, l)
    const int L1 = MD_GTO_MAX_ANGULAR_MOMENTUM + 1;
    int* group_of = (int*)malloc(sizeof(int) * num_atoms * L1);
    for (size_t i = 0; i < num_atoms * L1; ++i) group_of[i] = -1;

    size_t ao_idx = 0;
    bool ok = true;
    for (uint32_t si = 0; si < basis->num_shells; ++si) {
        const md_gto_shell_t* sh = &basis->shells[si];
        const int ncart = (int)md_gto_num_cart_ao(sh->l);
        if (sh->l > MD_GTO_MAX_ANGULAR_MOMENTUM) { ok = false; break; }

        bool any = false;
        for (int n = 0; n < ncart; ++n) any |= ao_coeffs[ao_idx + n] != 0.0;
        if (!any) { ao_idx += ncart; continue; }

        int* gi = &group_of[sh->atom_idx * L1 + sh->l];
        if (*gi < 0) {
            iso_shell_t g = {0};
            g.l = (int)sh->l;
            g.ncart = ncart;
            for (int k = 0; k < 3; ++k) g.c[k] = atom_xyz[3 * sh->atom_idx + k];
            for (int n = 0; n < ncart; ++n) {
                int i, j, k;
                md_gto_cart_ijk(&i, &j, &k, sh->l, (uint32_t)n);
                g.ijk[n][0] = (int8_t)i; g.ijk[n][1] = (int8_t)j; g.ijk[n][2] = (int8_t)k;
            }
            *gi = (int)md_array_size(f->shells);
            md_array_push(f->shells, g, heap);
        }
        iso_shell_t* g = &f->shells[*gi];
        for (uint32_t p = 0; p < sh->num_primitives; ++p) {
            const double a = basis->alpha[sh->primitive_offset + p];
            const double c = basis->coeff[sh->primitive_offset + p];
            int q = 0;
            while (q < g->nprim && g->alpha[q] != a) ++q;
            if (q == g->nprim) {
                if (g->nprim == ISO_MAX_PRIM) { ok = false; break; }
                g->alpha[g->nprim++] = a;
            }
            for (int n = 0; n < ncart; ++n) {
                const int i = g->ijk[n][0], j = g->ijk[n][1], k = g->ijk[n][2];
                g->w[n][q] += ao_coeffs[ao_idx + n] * md_gto_cart_norm_factor(i, j, k) * c;
            }
        }
        ao_idx += ncart;
    }
    free(group_of);
    if (!ok) return false;
    f->num_shells = (int)md_array_size(f->shells);

    for (int si = 0; si < f->num_shells; ++si) {
        iso_shell_t* s = &f->shells[si];
        for (int n = 0; n < s->ncart; ++n) {
            for (int q = 0; q < s->nprim; ++q) s->live[n] |= s->w[n][q] != 0.0;
        }
        double r = iso_shell_tail_r0(s);
        while (iso_shell_tail(s, r) > ISO_EVAL_EPS) r += 0.05;
        s->r2_cut = r * r;

        for (int n = 0; n < s->ncart; ++n) {
            for (int m = 0; m < ISO_NUM_CRIT; ++m) {
                const int first = (int)md_array_size(f->crit);
                s->crit_off[n][m] = first;
                if (s->live[n]) {
                    const double rho_end = iso_crit_rho_end(s, n, m + 1);
                    if (rho_end > 0.0) iso_crit_search(f, s, n, m, first, 0.0, rho_end, 0);
                }
                s->crit_num[n][m] = (int)md_array_size(f->crit) - first;
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------------------------
// Point evaluation: value, gradient, Hessian (xx, yy, zz, xy, xz, yz)
// ---------------------------------------------------------------------------------------------

static inline void iso_powers(double* out, double x, int l) {
    out[0] = 1.0;
    for (int i = 1; i <= l; ++i) out[i] = out[i - 1] * x;
}

// order 0: value, 1: + gradient, 2: + Hessian. Accumulates into val/grad/hess.
static inline void iso_shell_eval(const iso_shell_t* s, const double p[3], int order, double* val, double grad[3], double hess[6]) {
    static const int ia[6] = {0, 1, 2, 0, 0, 1};
    static const int ib[6] = {0, 1, 2, 1, 2, 2};
    const double d[3] = {p[0] - s->c[0], p[1] - s->c[1], p[2] - s->c[2]};
    const double rho = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

    double e[ISO_MAX_PRIM];
    for (int q = 0; q < s->nprim; ++q) e[q] = exp(-s->alpha[q] * rho);

    double px[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], py[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], pz[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
    iso_powers(px, d[0], s->l);
    iso_powers(py, d[1], s->l);
    iso_powers(pz, d[2], s->l);

    for (int n = 0; n < s->ncart; ++n) {
        if (!s->live[n]) continue;
        double R = 0, R1 = 0, R2 = 0;
        for (int q = 0; q < s->nprim; ++q) {
            const double t = s->w[n][q] * e[q];
            const double a = s->alpha[q];
            R  += t;
            R1 -= a * t;
            R2 += a * a * t;
        }
        const int i = s->ijk[n][0], j = s->ijk[n][1], k = s->ijk[n][2];
        const double M = px[i] * py[j] * pz[k];
        *val += M * R;
        if (order >= 1) {
            const double Mg[3] = {i ? i * px[i - 1] * py[j] * pz[k] : 0.0,
                                  j ? j * px[i] * py[j - 1] * pz[k] : 0.0,
                                  k ? k * px[i] * py[j] * pz[k - 1] : 0.0};
            const double gR[3] = {2 * d[0] * R1, 2 * d[1] * R1, 2 * d[2] * R1};
            for (int a = 0; a < 3; ++a) grad[a] += Mg[a] * R + M * gR[a];
            if (order >= 2) {
                const double Mh[6] = {i > 1 ? i * (i - 1) * px[i - 2] * py[j] * pz[k] : 0.0,
                                      j > 1 ? j * (j - 1) * px[i] * py[j - 2] * pz[k] : 0.0,
                                      k > 1 ? k * (k - 1) * px[i] * py[j] * pz[k - 2] : 0.0,
                                      (i && j) ? i * j * px[i - 1] * py[j - 1] * pz[k] : 0.0,
                                      (i && k) ? i * k * px[i - 1] * py[j] * pz[k - 1] : 0.0,
                                      (j && k) ? j * k * px[i] * py[j - 1] * pz[k - 1] : 0.0};
                for (int h = 0; h < 6; ++h) {
                    const int a = ia[h], b = ib[h];
                    const double HR = (a == b ? 2 * R1 : 0.0) + 4 * R2 * d[a] * d[b];
                    hess[h] += Mh[h] * R + Mg[a] * gR[b] + gR[a] * Mg[b] + M * HR;
                }
            }
        }
    }
}

// Full field at a point. active == NULL means all groups with the point evaluation cutoff applied.
static inline void iso_eval(const iso_field_t* f, const int* active, int num_active, const double p[3], int order, double* val, double grad[3], double hess[6]) {
    *val = 0;
    double g_dummy[3], h_dummy[6];
    if (!grad) grad = g_dummy;
    if (!hess) hess = h_dummy;
    grad[0] = grad[1] = grad[2] = 0;
    for (int h = 0; h < 6; ++h) hess[h] = 0;
    if (active) {
        for (int i = 0; i < num_active; ++i) iso_shell_eval(&f->shells[active[i]], p, order, val, grad, hess);
    } else {
        for (int i = 0; i < f->num_shells; ++i) {
            const iso_shell_t* s = &f->shells[i];
            const double dx = p[0] - s->c[0], dy = p[1] - s->c[1], dz = p[2] - s->c[2];
            if (dx * dx + dy * dy + dz * dz > s->r2_cut) continue;
            iso_shell_eval(s, p, order, val, grad, hess);
        }
    }
}

static inline double iso_value(const iso_field_t* f, const double p[3]) {
    double v;
    iso_eval(f, NULL, 0, p, 0, &v, NULL, NULL);
    return v;
}

// ---------------------------------------------------------------------------------------------
// Bounds over an axis aligned box
// ---------------------------------------------------------------------------------------------

enum {
    ISO_BOUND_PRIM = 0,
    ISO_BOUND_SHELL,
    ISO_BOUND_T1,
    ISO_BOUND_T2,
    ISO_BOUND_BEST,
    ISO_BOUND_COUNT
};

static const char* iso_bound_name[ISO_BOUND_COUNT] = {"prim", "group", "taylor1", "taylor2", "best"};

typedef struct iso_encl_t {
    iv_t   e[ISO_BOUND_COUNT];
    double mag;     // sum over shells of max|psi_s| over the box (scale for rounding slack)
    double grad_mag;// bound on |grad psi| over the box (when T1 is computed)
    double hess_diag[3]; // bounds on |d2 psi / dx_a^2| over the box (when T2 is computed)
} iso_encl_t;

// Exact range of t^n exp(-a t^2) over [b.lo, b.hi]
static inline iv_t iso_prim1d_range(int n, double a, iv_t b) {
    double lo = DBL_MAX, hi = -DBL_MAX;
#define ISO_CAND(t) do { const double _t = (t); const double _v = pow(_t, n) * exp(-a * _t * _t); lo = fmin(lo, _v); hi = fmax(hi, _v); } while (0)
    ISO_CAND(b.lo);
    ISO_CAND(b.hi);
    if (b.lo < 0.0 && b.hi > 0.0) ISO_CAND(0.0);
    if (n > 0) {
        const double t = sqrt(n / (2.0 * a));
        if (iv_contains(b, t))  ISO_CAND(t);
        if (iv_contains(b, -t)) ISO_CAND(-t);
    }
#undef ISO_CAND
    return iv(lo, hi);
}

typedef struct iso_box_t {
    double lo[3], hi[3];
} iso_box_t;

typedef struct iso_shell_box_t {
    iv_t   D[3];
    double rho_lo, rho_hi;
    iv_t   pw[3][MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
    double eu[ISO_MAX_PRIM], ev[ISO_MAX_PRIM];   // exp(-a_q rho_lo), exp(-a_q rho_hi): shared by every component and order
} iso_shell_box_t;

static inline void iso_shell_box(iso_shell_box_t* sb, const iso_shell_t* s, const iso_box_t* b) {
    sb->rho_lo = 0.0;
    sb->rho_hi = 0.0;
    for (int a = 0; a < 3; ++a) {
        sb->D[a] = iv(b->lo[a] - s->c[a], b->hi[a] - s->c[a]);
        const double l2 = sb->D[a].lo * sb->D[a].lo, h2 = sb->D[a].hi * sb->D[a].hi;
        if (!iv_contains(sb->D[a], 0.0)) sb->rho_lo += fmin(l2, h2);
        sb->rho_hi += fmax(l2, h2);
        for (int n = 0; n <= s->l; ++n) sb->pw[a][n] = iv_pow(sb->D[a], n);
    }
    for (int q = 0; q < s->nprim; ++q) {
        sb->eu[q] = exp(-s->alpha[q] * sb->rho_lo);
        sb->ev[q] = exp(-s->alpha[q] * sb->rho_hi);
    }
}

// iso_g_range over the box's [rho_lo, rho_hi], from the shared exponentials.
static inline iv_t iso_g_range_sb(const iso_field_t* f, const iso_shell_t* s, const iso_shell_box_t* sb, int n, int m) {
    double gu = 0.0, gv = 0.0;
    for (int q = 0; q < s->nprim; ++q) {
        const double w = iso_gw(s, q, n, m);
        gu += w * sb->eu[q];
        gv += w * sb->ev[q];
    }
    iv_t r = iv(fmin(gu, gv), fmax(gu, gv));
    const iso_crit_t* c = f->crit + s->crit_off[n][m];
    for (int i = 0; i < s->crit_num[n][m]; ++i) {
        if (c[i].hi >= sb->rho_lo && c[i].lo <= sb->rho_hi) r = iv_hull2(r, c[i].g);
    }
    return r;
}

static inline iv_t iso_mono(const iso_shell_box_t* sb, int i, int j, int k) {
    if (i < 0 || j < 0 || k < 0) return iv_pt(0.0);
    return iv_mul(iv_mul(sb->pw[0][i], sb->pw[1][j]), sb->pw[2][k]);
}

// GROUP bound of one (atom, l) group over the box: sum_n range(m_n) * range(R_n).
static inline iv_t iso_group_interval(const iso_field_t* f, const iso_shell_t* s, const iso_box_t* b) {
    iso_shell_box_t sb;
    iso_shell_box(&sb, s, b);
    iv_t r = iv_pt(0.0);
    for (int n = 0; n < s->ncart; ++n) {
        if (!s->live[n]) continue;
        const iv_t R = iso_g_range_sb(f, s, &sb, n, 0);
        r = iv_add(r, iv_mul(iso_mono(&sb, s->ijk[n][0], s->ijk[n][1], s->ijk[n][2]), R));
    }
    return r;
}

// Computes the enclosures selected by 'mask' (bit per ISO_BOUND_*) over the active shells.
// 'slack' is the carried bound of shells dropped higher up the tree.
static void iso_enclose(iso_encl_t* out, const iso_field_t* f, const int* active, int num_active, double slack, const iso_box_t* b, uint32_t mask) {
    const bool want_prim  = mask & (1u << ISO_BOUND_PRIM);
    const bool want_shell = mask & ((1u << ISO_BOUND_SHELL) | (1u << ISO_BOUND_BEST));
    const bool want_t1    = mask & ((1u << ISO_BOUND_T1) | (1u << ISO_BOUND_BEST));
    const bool want_t2    = mask & ((1u << ISO_BOUND_T2) | (1u << ISO_BOUND_BEST));

    iv_t prim = iv_pt(0.0), shell = iv_pt(0.0);
    iv_t G[3] = {iv_pt(0), iv_pt(0), iv_pt(0)};
    iv_t H[6] = {iv_pt(0), iv_pt(0), iv_pt(0), iv_pt(0), iv_pt(0), iv_pt(0)};
    double mag = 0.0;

    static const int ia[6] = {0, 1, 2, 0, 0, 1};
    static const int ib[6] = {0, 1, 2, 1, 2, 2};

    for (int si = 0; si < num_active; ++si) {
        const iso_shell_t* s = &f->shells[active[si]];
        iso_shell_box_t sb;
        iso_shell_box(&sb, s, b);

        if (want_prim) {
            iv_t pr[3][MD_GTO_MAX_ANGULAR_MOMENTUM + 1][ISO_MAX_PRIM];
            for (int a = 0; a < 3; ++a) for (int e = 0; e <= s->l; ++e) for (int q = 0; q < s->nprim; ++q) pr[a][e][q] = iso_prim1d_range(e, s->alpha[q], sb.D[a]);
            for (int n = 0; n < s->ncart; ++n) {
                if (!s->live[n]) continue;
                const int i = s->ijk[n][0], j = s->ijk[n][1], k = s->ijk[n][2];
                for (int q = 0; q < s->nprim; ++q) {
                    if (s->w[n][q] == 0.0) continue;
                    const iv_t t = iv_mul(iv_mul(pr[0][i][q], pr[1][j][q]), pr[2][k][q]);
                    prim = iv_add(prim, iv_scale(t, s->w[n][q]));
                }
            }
        }

        iv_t DD[6];
        if (want_t2) {
            for (int h = 0; h < 6; ++h) DD[h] = (ia[h] == ib[h]) ? iv_pow(sb.D[ia[h]], 2) : iv_mul(sb.D[ia[h]], sb.D[ib[h]]);
        }

        for (int n = 0; n < s->ncart; ++n) {
            if (!s->live[n]) continue;
            const int i = s->ijk[n][0], j = s->ijk[n][1], k = s->ijk[n][2];
            const iv_t M = iso_mono(&sb, i, j, k);
            const iv_t R = iso_g_range_sb(f, s, &sb, n, 0);
            const iv_t MR = iv_mul(M, R);
            mag += iv_mag(MR);
            if (want_shell) shell = iv_add(shell, MR);
            if (!(want_t1 || want_t2)) continue;

            const iv_t Mg[3] = {iv_scale(iso_mono(&sb, i - 1, j, k), i), iv_scale(iso_mono(&sb, i, j - 1, k), j), iv_scale(iso_mono(&sb, i, j, k - 1), k)};
            const iv_t R1 = iso_g_range_sb(f, s, &sb, n, 1);
            iv_t gR[3];
            for (int a = 0; a < 3; ++a) gR[a] = iv_scale(iv_mul(sb.D[a], R1), 2.0);
            if (want_t1) {
                for (int a = 0; a < 3; ++a) G[a] = iv_add(G[a], iv_add(iv_mul(Mg[a], R), iv_mul(M, gR[a])));
            }
            if (want_t2) {
                const iv_t R2 = iso_g_range_sb(f, s, &sb, n, 2);
                const iv_t Mh[6] = {iv_scale(iso_mono(&sb, i - 2, j, k), i * (i - 1)), iv_scale(iso_mono(&sb, i, j - 2, k), j * (j - 1)),
                                    iv_scale(iso_mono(&sb, i, j, k - 2), k * (k - 1)), iv_scale(iso_mono(&sb, i - 1, j - 1, k), i * j),
                                    iv_scale(iso_mono(&sb, i - 1, j, k - 1), i * k), iv_scale(iso_mono(&sb, i, j - 1, k - 1), j * k)};
                for (int h = 0; h < 6; ++h) {
                    const int a = ia[h], c = ib[h];
                    iv_t HR = iv_scale(iv_mul(DD[h], R2), 4.0);
                    if (a == c) HR = iv_add(HR, iv_scale(R1, 2.0));
                    iv_t t = iv_mul(Mh[h], R);
                    t = iv_add(t, iv_mul(Mg[a], gR[c]));
                    t = iv_add(t, iv_mul(gR[a], Mg[c]));
                    t = iv_add(t, iv_mul(M, HR));
                    H[h] = iv_add(H[h], t);
                }
            }
        }
    }

    double r[3], cen[3];
    for (int a = 0; a < 3; ++a) { r[a] = 0.5 * (b->hi[a] - b->lo[a]); cen[a] = 0.5 * (b->hi[a] + b->lo[a]); }

    double v0 = 0.0, g0[3] = {0};
    if (want_t1 || want_t2) iso_eval(f, active, num_active, cen, want_t2 ? 1 : 0, &v0, g0, NULL);

    for (int k = 0; k < ISO_BOUND_COUNT; ++k) out->e[k] = iv(-INFINITY, INFINITY);
    out->mag = mag;

#define ISO_FINAL(x) iv_widen((x), slack + 1.0e-12 * (mag + iv_width(x) + fabs(v0)))
    if (want_prim)  out->e[ISO_BOUND_PRIM]  = ISO_FINAL(prim);
    if (want_shell) out->e[ISO_BOUND_SHELL] = ISO_FINAL(shell);
    if (want_t1) {
        double w = 0.0;
        double g2 = 0.0;
        for (int a = 0; a < 3; ++a) { w += iv_mag(G[a]) * r[a]; g2 += iv_mag(G[a]) * iv_mag(G[a]); }
        out->e[ISO_BOUND_T1] = ISO_FINAL(iv(v0 - w, v0 + w));
        out->grad_mag = sqrt(g2) * (1.0 + 1.0e-12);
    }
    if (want_t2) {
        double lin = 0.0;
        for (int a = 0; a < 3; ++a) lin += fabs(g0[a]) * r[a];
        iv_t q = iv(v0 - lin, v0 + lin);
        for (int a = 0; a < 3; ++a) q = iv_add(q, iv_mul(iv_scale(H[a], 0.5), iv(0.0, r[a] * r[a])));
        for (int h = 3; h < 6; ++h) {
            const double rr = r[ia[h]] * r[ib[h]];
            const double m = iv_mag(H[h]) * rr;
            q = iv_add(q, iv(-m, m));
        }
        out->e[ISO_BOUND_T2] = ISO_FINAL(q);
        for (int a = 0; a < 3; ++a) out->hess_diag[a] = iv_mag(H[a]) * (1.0 + 1.0e-12);
    }
#undef ISO_FINAL
    if (mask & (1u << ISO_BOUND_BEST)) {
        iv_t e = iv_isect(iv_isect(out->e[ISO_BOUND_SHELL], out->e[ISO_BOUND_T1]), out->e[ISO_BOUND_T2]);
        if (want_prim) e = iv_isect(e, out->e[ISO_BOUND_PRIM]);
        out->e[ISO_BOUND_BEST] = e;
    }
}


// ---------------------------------------------------------------------------------------------
// Octree: lazily built from the bounds, never materializing a volume
// ---------------------------------------------------------------------------------------------

#define ISO_MAX_LEVEL 16
#define ISO_SCHEME_CASCADE 100   // GROUP first (free with the drop pass), then T2 (small boxes) or PRIM (large)

typedef struct iso_leaf_t { int32_t x, y, z; uint8_t sheets; } iso_leaf_t;

typedef struct iso_octree_t {
    const iso_field_t* f;
    double   iso;
    double   h;               // leaf edge
    double   origin[3];
    int      L;               // leaf level; root edge = h * 2^L
    int      scheme;          // ISO_BOUND_* deciding pruning, or ISO_SCHEME_CASCADE
    double   t2_max_size;     // cascade: T2 for boxes up to this edge, PRIM above
    uint32_t mask;            // bounds computed
    double   drop_eps;        // groups bounded below this are dropped for the subtree
    int*     active[ISO_MAX_LEVEL + 1];
    iso_leaf_t* leaves;       // md_array
    // stats
    long     visited[ISO_MAX_LEVEL + 1];
    long     kept[ISO_MAX_LEVEL + 1];
    long     group_work;      // sum over enclosures of active groups
    long     enclosures;      // enclosures beyond the GROUP pass every node gets
    long     decided_by[ISO_BOUND_COUNT];
    double   seconds;
} iso_octree_t;

static void iso_octree_visit(iso_octree_t* o, int level, int ix, int iy, int iz, const int* act, int nact, double slack) {
    const double size = o->h * (double)(1 << (o->L - level));
    iso_box_t b;
    b.lo[0] = o->origin[0] + ix * size; b.hi[0] = b.lo[0] + size;
    b.lo[1] = o->origin[1] + iy * size; b.hi[1] = b.lo[1] + size;
    b.lo[2] = o->origin[2] + iz * size; b.hi[2] = b.lo[2] + size;
    o->visited[level]++;

    // One pass: drop groups that are provably tiny here (their bound becomes slack for the whole
    // subtree), and sum the GROUP bound of the rest, which that test computes anyway.
    int* next = o->active[level];
    int n = 0;
    iv_t grp = iv_pt(0.0);
    double grp_mag = 0.0;
    for (int i = 0; i < nact; ++i) {
        const iv_t gi = iso_group_interval(o->f, &o->f->shells[act[i]], &b);
        const double m = iv_mag(gi);
        if (m < o->drop_eps) {
            slack += m;
        } else {
            next[n++] = act[i];
            grp = iv_add(grp, gi);
            grp_mag += m;
        }
    }
    o->group_work += nact;
    if (slack >= o->iso) return;   // cannot happen with sane drop_eps; keeps the guarantee honest

    uint8_t sheets = 0;
#define ISO_SHEETS(e) (uint8_t)((iv_contains((e), o->iso) ? 1 : 0) | (iv_contains((e), -o->iso) ? 2 : 0))
    if (n > 0) {
        const iv_t eg = iv_widen(grp, slack + 1.0e-12 * (grp_mag + iv_width(grp)));
        if (o->scheme == ISO_BOUND_SHELL) {
            sheets = ISO_SHEETS(eg);
        } else if (o->scheme == ISO_SCHEME_CASCADE) {
            sheets = ISO_SHEETS(eg);
            if (!sheets) {
                o->decided_by[ISO_BOUND_SHELL]++;
            } else {
                const int sc = (size <= o->t2_max_size) ? ISO_BOUND_T2 : ISO_BOUND_PRIM;
                iso_encl_t enc;
                iso_enclose(&enc, o->f, next, n, slack, &b, 1u << sc);
                o->enclosures++;
                o->group_work += n;
                sheets &= ISO_SHEETS(enc.e[sc]);
                if (!sheets) o->decided_by[sc]++;
            }
        } else {
            iso_encl_t enc;
            iso_enclose(&enc, o->f, next, n, slack, &b, o->mask);
            o->enclosures++;
            o->group_work += n;
            sheets = ISO_SHEETS(enc.e[o->scheme]);
        }
    }
#undef ISO_SHEETS
    if (!sheets) return;
    o->kept[level]++;

    if (level == o->L) {
        iso_leaf_t leaf = {ix, iy, iz, sheets};
        md_array_push(o->leaves, leaf, md_get_heap_allocator());
        return;
    }
    for (int c = 0; c < 8; ++c) {
        iso_octree_visit(o, level + 1, 2 * ix + (c & 1), 2 * iy + ((c >> 1) & 1), 2 * iz + ((c >> 2) & 1), next, n, slack);
    }
}

// The atom box padded until the sum of all group tails is below iso/2: no surface exists outside.
static void iso_surface_box(const iso_field_t* f, double iso, double lo[3], double hi[3]) {
    double r0 = 0.0;
    for (int i = 0; i < f->num_shells; ++i) r0 = fmax(r0, iso_shell_tail_r0(&f->shells[i]));
    double pad = r0;
    for (;; pad += 0.05) {
        double t = 0.0;
        for (int i = 0; i < f->num_shells; ++i) t += iso_shell_tail(&f->shells[i], pad);
        if (t < 0.5 * iso) break;
    }
    for (int a = 0; a < 3; ++a) { lo[a] = f->atom_min[a] - pad; hi[a] = f->atom_max[a] + pad; }
}

// Root: the cube of edge leaf_edge * 2^L centred on that box. Returns L and the root's origin.
static int iso_root(const iso_field_t* f, double iso, double leaf_edge, double origin[3]) {
    double lo[3], hi[3], ext = 0.0;
    iso_surface_box(f, iso, lo, hi);
    for (int a = 0; a < 3; ++a) ext = fmax(ext, hi[a] - lo[a]);
    int L = 0;
    while (leaf_edge * (double)(1 << L) < ext) L++;
    const double edge = leaf_edge * (double)(1 << L);
    for (int a = 0; a < 3; ++a) origin[a] = 0.5 * (lo[a] + hi[a]) - 0.5 * edge;
    return L;
}

static void iso_octree_build_ex(iso_octree_t* o, const iso_field_t* f, double iso, double h, int scheme, uint32_t mask, double t2_max_size) {
    memset(o, 0, sizeof(*o));
    o->f = f;
    o->iso = iso;
    o->h = h;
    o->scheme = scheme;
    o->t2_max_size = t2_max_size;
    o->mask = scheme == ISO_SCHEME_CASCADE ? 0 : (mask | (1u << scheme));
    o->drop_eps = 1.0e-4 * iso;

    o->L = iso_root(f, iso, h, o->origin);

    for (int l = 0; l <= o->L; ++l) o->active[l] = (int*)malloc(sizeof(int) * (f->num_shells + 1));
    int* all = (int*)malloc(sizeof(int) * (f->num_shells + 1));
    for (int i = 0; i < f->num_shells; ++i) all[i] = i;

    md_tick_t t0 = md_tick_now();
    iso_octree_visit(o, 0, 0, 0, 0, all, f->num_shells, 0.0);
    o->seconds = iso_seconds_since(t0);
    free(all);
}

static void iso_octree_build(iso_octree_t* o, const iso_field_t* f, double iso, double h, int scheme, uint32_t mask) {
    iso_octree_build_ex(o, f, iso, h, scheme, mask, 2.0);
}

static void iso_octree_free(iso_octree_t* o) {
    for (int l = 0; l <= ISO_MAX_LEVEL; ++l) free(o->active[l]);
    md_array_free(o->leaves, md_get_heap_allocator());
    memset(o, 0, sizeof(*o));
}

// Leaves of an octree built with leaf edge h * 2^k, expanded into the (2^k)^3 cells of edge h
// they contain. Same root, so the lattice is identical to the one of a tree built at h.
static iso_leaf_t* iso_bricks_to_cells(const iso_octree_t* o, int k, int* out_n) {
    const int B = 1 << k;
    const int nl = (int)md_array_size(o->leaves);
    iso_leaf_t* cells = (iso_leaf_t*)malloc(sizeof(iso_leaf_t) * (size_t)nl * B * B * B + 1);
    int n = 0;
    for (int li = 0; li < nl; ++li) {
        const iso_leaf_t* l = &o->leaves[li];
        for (int z = 0; z < B; ++z) for (int y = 0; y < B; ++y) for (int x = 0; x < B; ++x) {
            cells[n++] = (iso_leaf_t){l->x * B + x, l->y * B + y, l->z * B + z, l->sheets};
        }
    }
    *out_n = n;
    return cells;
}

// ---------------------------------------------------------------------------------------------
// Surface Nets over a sparse set of lattice cells, one sheet at a time (psi = +iso, psi = -iso)
// ---------------------------------------------------------------------------------------------

typedef struct iso_lattice_t {
    const iso_field_t* f;
    double     origin[3];
    double     h;
    iso_map_t  corner_map;
    double*    corner_val;    // md_array
    long       evals;
} iso_lattice_t;

static double iso_lattice_corner(iso_lattice_t* lat, int x, int y, int z) {
    const uint64_t key = iso_key3(x, y, z);
    uint32_t idx;
    if (iso_map_get(&lat->corner_map, key, &idx)) return lat->corner_val[idx];
    const double p[3] = {lat->origin[0] + x * lat->h, lat->origin[1] + y * lat->h, lat->origin[2] + z * lat->h};
    const double v = iso_value(lat->f, p);
    lat->evals++;
    iso_map_put(&lat->corner_map, key, (uint32_t)md_array_size(lat->corner_val));
    md_array_push(lat->corner_val, v, md_get_heap_allocator());
    return v;
}

typedef struct iso_mesh_t {
    int       sheet;          // +1 or -1
    double    iso;
    double    h;
    double*   pos;            // md_array, 3 per vertex
    double*   nrm;            // md_array, 3 per vertex, outward unit normal
    double*   res;            // md_array, 1 per vertex: |f|/|grad f| after the last projection
    uint32_t* quad;           // md_array, 4 per quad
    uint32_t* tri;            // md_array, 3 per triangle
    int       open_edges;     // crossing edges whose four cells were not all present (must be 0)
    int       ambiguous;      // cells the bound kept whose corners show no crossing
} iso_mesh_t;

static inline int iso_mesh_nv(const iso_mesh_t* m) { return (int)(md_array_size(m->pos) / 3); }
static inline int iso_mesh_nt(const iso_mesh_t* m) { return (int)(md_array_size(m->tri) / 3); }

static void iso_mesh_free(iso_mesh_t* m) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_array_free(m->pos, heap); md_array_free(m->nrm, heap); md_array_free(m->res, heap);
    md_array_free(m->quad, heap); md_array_free(m->tri, heap);
    memset(m, 0, sizeof(*m));
}

static const int iso_cube_edges[12][2] = {
    {0, 1}, {2, 3}, {4, 5}, {6, 7},   // x
    {0, 2}, {1, 3}, {4, 6}, {5, 7},   // y
    {0, 4}, {1, 5}, {2, 6}, {3, 7},   // z
};

// cells: (x,y,z) triples. sheet_bits: per cell, which sheets the cell may contain (bit 0: +, bit 1: -).
static void iso_surface_nets(iso_mesh_t* out, iso_lattice_t* lat, const iso_leaf_t* cells, int num_cells, int sheet, double iso) {
    md_allocator_i* heap = md_get_heap_allocator();
    memset(out, 0, sizeof(*out));
    out->sheet = sheet;
    out->iso = iso;
    out->h = lat->h;
    const uint8_t bit = sheet > 0 ? 1 : 2;

    iso_map_t cell_map;
    iso_map_init(&cell_map, (uint32_t)num_cells);

    // vertices
    for (int ci = 0; ci < num_cells; ++ci) {
        const iso_leaf_t* c = &cells[ci];
        if (!(c->sheets & bit)) continue;
        double fv[8];
        int mask = 0;
        for (int k = 0; k < 8; ++k) {
            fv[k] = sheet * iso_lattice_corner(lat, c->x + (k & 1), c->y + ((k >> 1) & 1), c->z + ((k >> 2) & 1)) - iso;
            if (fv[k] > 0) mask |= 1 << k;
        }
        if (mask == 0 || mask == 0xFF) { out->ambiguous++; continue; }
        double acc[3] = {0}, cnt = 0;
        for (int e = 0; e < 12; ++e) {
            const int a = iso_cube_edges[e][0], b = iso_cube_edges[e][1];
            if (((mask >> a) & 1) == ((mask >> b) & 1)) continue;
            const double t = fv[a] / (fv[a] - fv[b]);
            for (int d = 0; d < 3; ++d) {
                const double pa = (a >> d) & 1, pb = (b >> d) & 1;
                acc[d] += pa + t * (pb - pa);
            }
            cnt += 1;
        }
        const double p[3] = {lat->origin[0] + (c->x + acc[0] / cnt) * lat->h,
                             lat->origin[1] + (c->y + acc[1] / cnt) * lat->h,
                             lat->origin[2] + (c->z + acc[2] / cnt) * lat->h};
        iso_map_put(&cell_map, iso_key3(c->x, c->y, c->z), (uint32_t)iso_mesh_nv(out));
        md_array_push_array(out->pos, p, 3, heap);
    }

    // quads: one per crossing lattice edge, from the four cells around it
    for (int ci = 0; ci < num_cells; ++ci) {
        const iso_leaf_t* c = &cells[ci];
        uint32_t v0;
        if (!iso_map_get(&cell_map, iso_key3(c->x, c->y, c->z), &v0)) continue;
        const double f0 = sheet * iso_lattice_corner(lat, c->x, c->y, c->z) - iso;
        for (int a = 0; a < 3; ++a) {
            const int o[3] = {a == 0, a == 1, a == 2};
            const double f1 = sheet * iso_lattice_corner(lat, c->x + o[0], c->y + o[1], c->z + o[2]) - iso;
            if ((f0 > 0) == (f1 > 0)) continue;
            const int u = (a + 1) % 3, w = (a + 2) % 3;
            int q[4][3];
            for (int k = 0; k < 4; ++k) { q[k][0] = c->x; q[k][1] = c->y; q[k][2] = c->z; }
            q[1][u] -= 1;
            q[2][u] -= 1; q[2][w] -= 1;
            q[3][w] -= 1;
            uint32_t vi[4];
            bool ok = true;
            for (int k = 0; k < 4 && ok; ++k) ok = iso_map_get(&cell_map, iso_key3(q[k][0], q[k][1], q[k][2]), &vi[k]);
            if (!ok) { out->open_edges++; continue; }
            if (f0 <= 0) { const uint32_t t = vi[1]; vi[1] = vi[3]; vi[3] = t; }
            md_array_push_array(out->quad, vi, 4, heap);
        }
    }
    iso_map_free(&cell_map);
}

// ---------------------------------------------------------------------------------------------
// Mesh refinement against the analytic field
// ---------------------------------------------------------------------------------------------

// Newton projection of p onto sheet*psi = iso. Step clamped to max_step. Returns |f|/|grad f|,
// writes the outward unit normal.
static double iso_project_point(const iso_field_t* f, double p[3], int sheet, double iso, double max_step, int max_iter, double nrm[3]) {
    double res = 0.0;
    for (int it = 0; it <= max_iter; ++it) {
        double v, g[3];
        iso_eval(f, NULL, 0, p, 1, &v, g, NULL);
        const double fv = sheet * v - iso;
        for (int a = 0; a < 3; ++a) g[a] *= sheet;
        const double g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        if (g2 <= 1e-300) { res = INFINITY; break; }
        const double gl = sqrt(g2);
        for (int a = 0; a < 3; ++a) nrm[a] = -g[a] / gl;   // f > 0 inside, so outward is -grad f
        res = fabs(fv) / gl;
        if (it == max_iter || res < 1.0e-9) break;
        double step = fv / g2;
        const double len = fabs(step) * gl;
        if (len > max_step) step *= max_step / len;
        for (int a = 0; a < 3; ++a) p[a] -= step * g[a];
    }
    return res;
}

static void iso_mesh_project(iso_mesh_t* m, const iso_field_t* f, int max_iter) {
    md_allocator_i* heap = md_get_heap_allocator();
    const int nv = iso_mesh_nv(m);
    md_array_resize(m->nrm, 3 * nv, heap);
    md_array_resize(m->res, nv, heap);
    for (int i = 0; i < nv; ++i) {
        m->res[i] = iso_project_point(f, m->pos + 3 * i, m->sheet, m->iso, 0.5 * m->h, max_iter, m->nrm + 3 * i);
    }
}

static inline double iso_dist2(const double* a, const double* b) {
    const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
}

// Quads to triangles along the shorter diagonal of the (projected) quad.
static void iso_mesh_triangulate(iso_mesh_t* m) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_array_shrink(m->tri, 0);
    const int nq = (int)(md_array_size(m->quad) / 4);
    for (int i = 0; i < nq; ++i) {
        const uint32_t* q = m->quad + 4 * i;
        const double d02 = iso_dist2(m->pos + 3 * q[0], m->pos + 3 * q[2]);
        const double d13 = iso_dist2(m->pos + 3 * q[1], m->pos + 3 * q[3]);
        if (d02 <= d13) {
            const uint32_t t[6] = {q[0], q[1], q[2], q[0], q[2], q[3]};
            md_array_push_array(m->tri, t, 6, heap);
        } else {
            const uint32_t t[6] = {q[0], q[1], q[3], q[1], q[2], q[3]};
            md_array_push_array(m->tri, t, 6, heap);
        }
    }
}

// Vertex one-rings as CSR from the triangle edges.
static void iso_mesh_adjacency(const iso_mesh_t* m, int** out_off, int** out_adj) {
    const int nv = iso_mesh_nv(m), nt = iso_mesh_nt(m);
    int* deg = (int*)calloc(nv + 1, sizeof(int));
    iso_map_t seen;
    iso_map_init(&seen, (uint32_t)(nt * 3));
    uint32_t* pairs = (uint32_t*)malloc(sizeof(uint32_t) * nt * 6);
    int np = 0;
    for (int t = 0; t < nt; ++t) {
        for (int e = 0; e < 3; ++e) {
            uint32_t a = m->tri[3 * t + e], b = m->tri[3 * t + (e + 1) % 3];
            if (a > b) { uint32_t x = a; a = b; b = x; }
            const uint64_t key = ((uint64_t)a << 32) | b;
            uint32_t dummy;
            if (iso_map_get(&seen, key, &dummy)) continue;
            iso_map_put(&seen, key, 1);
            pairs[2 * np] = a; pairs[2 * np + 1] = b; np++;
            deg[a]++; deg[b]++;
        }
    }
    int* off = (int*)malloc(sizeof(int) * (nv + 1));
    off[0] = 0;
    for (int i = 0; i < nv; ++i) off[i + 1] = off[i] + deg[i];
    int* adj = (int*)malloc(sizeof(int) * (off[nv] + 1));
    for (int i = 0; i < nv; ++i) deg[i] = off[i];
    for (int i = 0; i < np; ++i) {
        adj[deg[pairs[2 * i]]++] = (int)pairs[2 * i + 1];
        adj[deg[pairs[2 * i + 1]]++] = (int)pairs[2 * i];
    }
    free(pairs); free(deg);
    iso_map_free(&seen);
    *out_off = off;
    *out_adj = adj;
}

// Tangential Laplacian relaxation: move toward the one-ring centroid within the tangent plane,
// then reproject. Jacobi style so the result does not depend on vertex order.
static void iso_mesh_relax(iso_mesh_t* m, const iso_field_t* f, int iterations, double lambda) {
    int *off, *adj;
    iso_mesh_adjacency(m, &off, &adj);
    const int nv = iso_mesh_nv(m);
    double* np = (double*)malloc(sizeof(double) * 3 * nv);
    for (int it = 0; it < iterations; ++it) {
        for (int i = 0; i < nv; ++i) {
            const double* p = m->pos + 3 * i;
            const double* n = m->nrm + 3 * i;
            double c[3] = {0}, k = off[i + 1] - off[i];
            for (int j = off[i]; j < off[i + 1]; ++j) for (int a = 0; a < 3; ++a) c[a] += m->pos[3 * adj[j] + a];
            if (k == 0) { memcpy(np + 3 * i, p, sizeof(double) * 3); continue; }
            double d[3] = {c[0] / k - p[0], c[1] / k - p[1], c[2] / k - p[2]};
            const double dn = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
            for (int a = 0; a < 3; ++a) np[3 * i + a] = p[a] + lambda * (d[a] - dn * n[a]);
        }
        memcpy(m->pos, np, sizeof(double) * 3 * nv);
        iso_mesh_project(m, f, 3);
    }
    free(np); free(off); free(adj);
}

// Adaptive refinement: split every edge whose projected midpoint lies further than tol from the
// chord midpoint, or whose end normals differ by more than max_angle. Red-green: a triangle with
// 1, 2 or 3 split edges becomes 2, 3 or 4 triangles, so the result is conforming without
// propagation. New vertices are projected.
static int iso_mesh_refine(iso_mesh_t* m, const iso_field_t* f, double tol, double max_angle_deg) {
    md_allocator_i* heap = md_get_heap_allocator();
    const double cos_max = cos(max_angle_deg * 3.14159265358979 / 180.0);
    const int nt = iso_mesh_nt(m);
    iso_map_t mid;
    iso_map_init(&mid, (uint32_t)(nt * 2));
    iso_map_t tested;
    iso_map_init(&tested, (uint32_t)(nt * 2));
    int splits = 0;

    for (int t = 0; t < nt; ++t) {
        for (int e = 0; e < 3; ++e) {
            uint32_t a = m->tri[3 * t + e], b = m->tri[3 * t + (e + 1) % 3];
            if (a > b) { uint32_t x = a; a = b; b = x; }
            const uint64_t key = ((uint64_t)a << 32) | b;
            uint32_t dummy;
            if (iso_map_get(&tested, key, &dummy)) continue;
            iso_map_put(&tested, key, 1);
            double p[3], n[3];
            for (int k = 0; k < 3; ++k) p[k] = 0.5 * (m->pos[3 * a + k] + m->pos[3 * b + k]);
            const double chord[3] = {p[0], p[1], p[2]};
            const double r = iso_project_point(f, p, m->sheet, m->iso, 0.5 * m->h, 6, n);
            const double* na = m->nrm + 3 * a;
            const double* nb = m->nrm + 3 * b;
            const double cosab = na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2];
            if (sqrt(iso_dist2(p, chord)) > tol || cosab < cos_max) {
                const uint32_t vi = (uint32_t)iso_mesh_nv(m);
                md_array_push_array(m->pos, p, 3, heap);
                md_array_push_array(m->nrm, n, 3, heap);
                md_array_push(m->res, r, heap);
                iso_map_put(&mid, key, vi);
                splits++;
            }
        }
    }
    iso_map_free(&tested);
    if (!splits) { iso_map_free(&mid); return 0; }

    uint32_t* tri = NULL;
    for (int t = 0; t < nt; ++t) {
        const uint32_t v[3] = {m->tri[3 * t], m->tri[3 * t + 1], m->tri[3 * t + 2]};
        uint32_t mv[3];
        bool has[3];
        int cnt = 0;
        for (int e = 0; e < 3; ++e) {
            uint32_t a = v[e], b = v[(e + 1) % 3];
            if (a > b) { uint32_t x = a; a = b; b = x; }
            has[e] = iso_map_get(&mid, ((uint64_t)a << 32) | b, &mv[e]);
            cnt += has[e];
        }
#define ISO_TRI(a, b, c) do { const uint32_t _t[3] = {(a), (b), (c)}; md_array_push_array(tri, _t, 3, heap); } while (0)
        if (cnt == 0) {
            ISO_TRI(v[0], v[1], v[2]);
        } else if (cnt == 3) {
            ISO_TRI(v[0], mv[0], mv[2]);
            ISO_TRI(v[1], mv[1], mv[0]);
            ISO_TRI(v[2], mv[2], mv[1]);
            ISO_TRI(mv[0], mv[1], mv[2]);
        } else if (cnt == 1) {
            const int e = has[0] ? 0 : has[1] ? 1 : 2;
            ISO_TRI(v[e], mv[e], v[(e + 2) % 3]);
            ISO_TRI(mv[e], v[(e + 1) % 3], v[(e + 2) % 3]);
        } else {
            // two split edges: the unsplit one is e; corner opposite to it is v[(e+2)%3]
            const int e = !has[0] ? 0 : !has[1] ? 1 : 2;
            const uint32_t A = v[e], B = v[(e + 1) % 3], C = v[(e + 2) % 3];
            const uint32_t mBC = mv[(e + 1) % 3], mCA = mv[(e + 2) % 3];
            ISO_TRI(C, mCA, mBC);
            // quad A, B, mBC, mCA: shorter diagonal
            if (iso_dist2(m->pos + 3 * A, m->pos + 3 * mBC) <= iso_dist2(m->pos + 3 * B, m->pos + 3 * mCA)) {
                ISO_TRI(A, B, mBC);
                ISO_TRI(A, mBC, mCA);
            } else {
                ISO_TRI(A, B, mCA);
                ISO_TRI(B, mBC, mCA);
            }
        }
#undef ISO_TRI
    }
    md_array_free(m->tri, heap);
    m->tri = tri;
    iso_map_free(&mid);
    return splits;
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
            const double d = fabs(m->sheet * val - m->iso) / fmax(gl, 1e-300) / m->h;
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
        int violations = 0, num_crit = 0, num_comp = 0;
        double worst = 0.0;
        for (int si = 0; si < f.num_shells; ++si) {
            const iso_shell_t* s = &f.shells[si];
            for (int n = 0; n < s->ncart; ++n) {
                if (!s->live[n]) continue;
                num_comp++;
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
            tc->name, f.num_input_shells, f.num_shells, num_comp, num_crit, violations, worst);
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
// GPU: the same pipeline on md_gpu (shaders/gto_iso.slang)
// =============================================================================================
//
// The production shape the CPU experiments pointed at: GROUP bound only, octree stopped at 4^3
// cell bricks, bricks filled densely. One launch per octree level, then per brick / cell /
// vertex / quad. Every launch after the first is sized by a device-side count (iso_grid writes
// the indirect grid); the host reads back once, the brick count after the octree, to size the
// mesh buffers.

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
    md_gpu_addr_t groups, alpha, weights, comps, crit;
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
    md_gpu_addr_t groups, alpha, weights, comps, crit;
    uint32_t      num_groups;
} iso_gpu_field_t;

static void iso_gpu_field_upload(iso_gpu_t* g, iso_gpu_field_t* out, const iso_field_t* f, double iso) {
    iso_gpu_group_t* groups = (iso_gpu_group_t*)calloc(f->num_shells + 1, sizeof(iso_gpu_group_t));
    float* alpha = NULL; float* w = NULL; md_gpu_uint4* comps = NULL; md_gpu_float4* crit = NULL;
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
        G->info2 = (md_gpu_uint4){(uint32_t)md_array_size(w), (uint32_t)md_array_size(comps), 0, 0};
        for (int q = 0; q < s->nprim; ++q) { float a = (float)s->alpha[q]; md_array_push(alpha, a, heap); }
        for (int n = 0; n < s->ncart; ++n) {
            for (int q = 0; q < s->nprim; ++q) { float x = (float)s->w[n][q]; md_array_push(w, x, heap); }
            md_gpu_uint4 cp = {(uint32_t)(s->ijk[n][0] | (s->ijk[n][1] << 4) | (s->ijk[n][2] << 8)), (uint32_t)md_array_size(crit), (uint32_t)s->crit_num[n][0], s->live[n] ? 1u : 0u};
            for (int k = 0; k < s->crit_num[n][0]; ++k) {
                const iso_crit_t* c = &f->crit[s->crit_off[n][0] + k];
                md_gpu_float4 cf = {iso_f32_down(c->lo), iso_f32_up(c->hi), iso_f32_down(c->g.lo - 1e-7 * fabs(c->g.lo)), iso_f32_up(c->g.hi + 1e-7 * fabs(c->g.hi))};
                md_array_push(crit, cf, heap);
            }
            md_array_push(comps, cp, heap);
        }
    }
    out->num_groups = (uint32_t)f->num_shells;
    out->groups  = iso_gpu_upload(g, groups, sizeof(iso_gpu_group_t) * f->num_shells);
    out->alpha   = iso_gpu_upload(g, alpha, sizeof(float) * md_array_size(alpha));
    out->weights = iso_gpu_upload(g, w, sizeof(float) * md_array_size(w));
    out->comps   = iso_gpu_upload(g, comps, sizeof(md_gpu_uint4) * md_array_size(comps));
    out->crit    = iso_gpu_upload(g, crit, sizeof(md_gpu_float4) * md_array_size(crit));
    free(groups);
    md_array_free(alpha, heap); md_array_free(w, heap); md_array_free(comps, heap); md_array_free(crit, heap);
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
    a.groups = fld.groups; a.alpha = fld.alpha; a.weights = fld.weights; a.comps = fld.comps; a.crit = fld.crit;
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
#endif // MD_ENABLE_GPU
