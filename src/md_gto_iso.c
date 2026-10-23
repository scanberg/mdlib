// md_gto_iso.c - isosurface meshes of orbitals and densities, straight from the basis. See md_gto_iso.h.
//
// How it works, and why it can be trusted
// ---------------------------------------
// An orbital: psi = sum_g sum_n x^i y^j z^k R_gn(|d|^2), one group g per (atom, l) and one radial
// function R_gn per Cartesian component (all shells of an atom with the same l share the
// monomials). R_gn is 1D, so its range over an interval of rho = |d|^2 is exact from precomputed
// critical intervals between which it is monotone. That makes a cheap and tight bound of psi over
// any box (the GROUP bound): sum_n range(monomial over box) * range(R_n over the box's rho interval).
//
// A density: rho = sum_t w_t a_t b_t over orbital-like functions (see "Density fields"). Each
// function is bounded over a box from the AO ranges of the shells that reach it, and the terms as
// interval products (exact squares). Where the products come from (occupied MOs, NTO pairs, or a
// density matrix split here) only decides how many functions there are.
//
// An octree is built top down from the root cube (padded atom box, outside which the tails provably
// stay below iso / 2). A node is kept only when its bound admits +iso or -iso; contributions whose
// bound over a node is below a small epsilon are dropped for the subtree and carried as slack, which
// keeps the bound rigorous. The octree stops at bricks of 4^3 cells: past that point, evaluating
// every corner of a brick is cheaper than bounding smaller boxes.
//
// Surface Nets on the kept cells is the same mesh as Surface Nets on the full lattice (no crossing
// cell can be discarded), then every vertex is moved onto the true surface by Newton steps along
// the analytic gradient.
//
// Other orbital bounds (per primitive, Taylor order 1 and 2) are here too: they are tighter on small
// boxes but cost more than they save in the octree. unittest/test_gto_iso.c compares them, fuzzes
// every bound for containment, and holds the sparse mesh to the dense one, for orbitals and for
// densities.

#include <md_gto_iso.h>

#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
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
static inline double iso_ipow(double x, int n) {
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

static inline iv_t iv_pow(iv_t a, int n) {
    if (n == 0) return iv_pt(1.0);
    if (n == 1) return a;
    const double pl = iso_ipow(a.lo, n), ph = iso_ipow(a.hi, n);
    if (n & 1) return (iv_t){pl, ph};
    if (a.lo >= 0) return (iv_t){pl, ph};
    if (a.hi <= 0) return (iv_t){ph, pl};
    return (iv_t){0.0, fmax(pl, ph)};
}

// ---------------------------------------------------------------------------------------------
// Small utilities: rng, u64 -> u32 hash map, dynamic arrays via md_array
// ---------------------------------------------------------------------------------------------

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
    int     fn;                              // the function the group belongs to (0 for an orbital)
    int     site;                            // the (atom, l) it sits on: groups of one site share c, l, alpha
} iso_shell_t;

struct iso_density_t;

// An orbital: the groups of one function. A density: the groups of all its functions, ordered by
// site, plus the density description (dens) that says how they combine; see "Density fields" below.
typedef struct iso_field_t {
    iso_shell_t* shells;     // md_array
    iso_crit_t*  crit;       // md_array, shared storage
    int          num_shells;
    int          num_input_shells;
    double       atom_min[3], atom_max[3];
    char         name[64];
    struct iso_density_t* dens;   // NULL for an orbital; owned
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

static void iso_density_free(struct iso_density_t* d);

static void iso_field_free(iso_field_t* f) {
    md_array_free(f->shells, md_get_heap_allocator());
    md_array_free(f->crit, md_get_heap_allocator());
    if (f->dens) iso_density_free(f->dens);
    memset(f, 0, sizeof(*f));
}

// Critical intervals of g_0 .. g_(orders-1) for every live component of the group.
static void iso_shell_finish(iso_field_t* f, iso_shell_t* s, int orders) {
    for (int n = 0; n < s->ncart; ++n) {
        s->live[n] = false;
        for (int q = 0; q < s->nprim; ++q) s->live[n] |= s->w[n][q] != 0.0;
    }
    double r = iso_shell_tail_r0(s);
    while (iso_shell_tail(s, r) > ISO_EVAL_EPS) r += 0.05;
    s->r2_cut = r * r;
    for (int n = 0; n < s->ncart; ++n) {
        for (int m = 0; m < ISO_NUM_CRIT; ++m) {
            const int first = (int)md_array_size(f->crit);
            s->crit_off[n][m] = first;
            if (s->live[n] && m < orders) {
                const double rho_end = iso_crit_rho_end(s, n, m + 1);
                if (rho_end > 0.0) iso_crit_search(f, s, n, m, first, 0.0, rho_end, 0);
            }
            s->crit_num[n][m] = (int)md_array_size(f->crit) - first;
        }
    }
}

// Builds the field from a basis, positions (bohr) and one Cartesian AO coefficient vector.
static bool iso_field_init(iso_field_t* f, const md_gto_basis_t* basis, const float* atom_xyz, size_t num_atoms, const double* ao_coeffs, const char* name) {
    memset(f, 0, sizeof(*f));
    md_snprintf(f->name, sizeof(f->name), "%s", name);
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
            g.site = (int)md_array_size(f->shells);
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

    for (int si = 0; si < f->num_shells; ++si) iso_shell_finish(f, &f->shells[si], ISO_NUM_CRIT);
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

static void iso_density_eval(const struct iso_density_t* d, const double p[3], int order, double* val, double grad[3]);

// Full field at a point. active == NULL means all groups with the point evaluation cutoff applied.
// A density ignores 'active' and has no Hessian (order <= 1).
static inline void iso_eval(const iso_field_t* f, const int* active, int num_active, const double p[3], int order, double* val, double grad[3], double hess[6]) {
    *val = 0;
    double g_dummy[3], h_dummy[6];
    if (!grad) grad = g_dummy;
    if (!hess) hess = h_dummy;
    grad[0] = grad[1] = grad[2] = 0;
    for (int h = 0; h < 6; ++h) hess[h] = 0;
    if (f->dens) {
        iso_density_eval(f->dens, p, order, val, grad);
        return;
    }
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
static inline iv_t iso_group_interval_sb(const iso_field_t* f, const iso_shell_t* s, const iso_shell_box_t* sb) {
    iv_t r = iv_pt(0.0);
    for (int n = 0; n < s->ncart; ++n) {
        if (!s->live[n]) continue;
        const iv_t R = iso_g_range_sb(f, s, sb, n, 0);
        r = iv_add(r, iv_mul(iso_mono(sb, s->ijk[n][0], s->ijk[n][1], s->ijk[n][2]), R));
    }
    return r;
}

static inline iv_t iso_group_interval(const iso_field_t* f, const iso_shell_t* s, const iso_box_t* b) {
    iso_shell_box_t sb;
    iso_shell_box(&sb, s, b);
    return iso_group_interval_sb(f, s, &sb);
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
// Density fields: rho = sum_t w_t a_t b_t
// ---------------------------------------------------------------------------------------------
//
// Every one-particle density of a GTO basis is a quadratic form rho(r) = phi(r)^T D phi(r) in the
// AO values. Only the symmetric part of D contributes (phi^T D phi = phi^T D^T phi), and any
// splitting of it into a short sum of products of orbital-like functions,
//
//     rho = sum_t w_t a_t(r) b_t(r),   a_t = A_t . phi,  b_t = B_t . phi,
//
// turns the density into something the orbital machinery bounds: every function is a linear
// combination of AOs, grouped by (atom, l) like an orbital, its GROUP bound over a box is an interval,
// and a product of two intervals - or the exact square of one - bounds each term. Where the products
// come from does not matter to the field:
//
//   * an SCF density is sum_i n_i psi_i^2 over its occupied MOs: no decomposition at all,
//   * the transition density of an excitation is sum_k s_k h_k p_k over its NTO pairs, and its
//     attachment / detachment densities are sums of squares of the same NTOs,
//   * a matrix alone (relaxed or difference densities, densities read from a file) is factored
//     here: pivoted LDL^T (pivoted Cholesky for a semidefinite D; its factors are localized, which
//     keeps the per-function bounds tight and lets the octree drop them away from where they live),
//     falling back to an eigendecomposition when an indefinite D makes 1x1 pivots unstable. The
//     residual of the factorization is bounded once, sum_ij |R_ij| max|phi_i| max|phi_j|, and added
//     to every bound, so a factored matrix is enclosed as rigorously as the products are.
//
// Point values come from the AOs (once per point) and the functions' AO coefficients, so the cost
// per point is the AO evaluation plus num_ao x num_functions multiply-adds.

typedef struct iso_term_t { double w; int a, b; } iso_term_t;   // a == b: a square, b < 0: linear (w a)

typedef struct iso_ao_shell_t {     // one contracted shell of the basis, for point evaluation
    double c[3];
    int    l, ncart, nprim, ao;
    double alpha[ISO_MAX_PRIM];
    double coef[ISO_MAX_PRIM];
    double nrm[ISO_MAX_CART];
    int8_t ijk[ISO_MAX_CART][3];
    double r2_cut;                  // every AO of the shell below ISO_EVAL_EPS beyond this
} iso_ao_shell_t;

enum {
    ISO_FACTOR_NONE = 0,            // products given by the caller
    ISO_FACTOR_LDL,                 // pivoted LDL^T of a matrix
    ISO_FACTOR_EIGEN,               // eigendecomposition of a matrix
};

typedef struct iso_density_t {
    int             num_fn, num_terms, num_ao, num_ao_shells;
    iso_term_t*     terms;
    double*         C;              // [num_ao][num_fn]: AO coefficients of every function
    double*         Cabs;           // |C|
    iso_ao_shell_t* ao_shells;
    iso_field_t     radial;         // one group per basis shell: its contracted radial function (ncart 1)
    double          slack;          // added to every bound: the factorization residual, if any
    int             factor;         // ISO_FACTOR_*
} iso_density_t;

#define ISO_DENSITY_STACK_FN 256

static void iso_density_free(iso_density_t* d) {
    free(d->terms);
    free(d->C);
    free(d->Cabs);
    free(d->ao_shells);
    md_array_free(d->radial.shells, md_get_heap_allocator());
    md_array_free(d->radial.crit, md_get_heap_allocator());
    free(d);
}

static void iso_density_eval(const iso_density_t* d, const double p[3], int order, double* val, double grad[3]) {
    const int nf = d->num_fn;
    double stack[4 * ISO_DENSITY_STACK_FN];
    double* buf = nf <= ISO_DENSITY_STACK_FN ? stack : (double*)malloc(sizeof(double) * 4 * nf);
    double* a  = buf;               // function values
    double* ga = buf + nf;          // gradients, [3][nf]
    const int nc = order >= 1 ? 4 : 1;
    memset(buf, 0, sizeof(double) * nc * nf);

    for (int si = 0; si < d->num_ao_shells; ++si) {
        const iso_ao_shell_t* s = &d->ao_shells[si];
        const double dd[3] = {p[0] - s->c[0], p[1] - s->c[1], p[2] - s->c[2]};
        const double rho = dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2];
        if (rho > s->r2_cut) continue;
        double R = 0.0, R1 = 0.0;
        for (int q = 0; q < s->nprim; ++q) {
            const double t = s->coef[q] * exp(-s->alpha[q] * rho);
            R  += t;
            R1 -= s->alpha[q] * t;
        }
        double px[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], py[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], pz[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
        iso_powers(px, dd[0], s->l);
        iso_powers(py, dd[1], s->l);
        iso_powers(pz, dd[2], s->l);
        for (int n = 0; n < s->ncart; ++n) {
            const int i = s->ijk[n][0], j = s->ijk[n][1], k = s->ijk[n][2];
            const double M = px[i] * py[j] * pz[k];
            const double phi = s->nrm[n] * M * R;
            const double* c = d->C + (size_t)(s->ao + n) * nf;
            for (int fi = 0; fi < nf; ++fi) a[fi] += c[fi] * phi;
            if (order >= 1) {
                const double g[3] = {
                    s->nrm[n] * ((i ? i * px[i - 1] * py[j] * pz[k] : 0.0) * R + M * 2.0 * dd[0] * R1),
                    s->nrm[n] * ((j ? j * px[i] * py[j - 1] * pz[k] : 0.0) * R + M * 2.0 * dd[1] * R1),
                    s->nrm[n] * ((k ? k * px[i] * py[j] * pz[k - 1] : 0.0) * R + M * 2.0 * dd[2] * R1),
                };
                for (int ax = 0; ax < 3; ++ax) {
                    double* gx = ga + ax * nf;
                    for (int fi = 0; fi < nf; ++fi) gx[fi] += c[fi] * g[ax];
                }
            }
        }
    }

    double v = 0.0, g[3] = {0, 0, 0};
    for (int t = 0; t < d->num_terms; ++t) {
        const iso_term_t* tm = &d->terms[t];
        if (tm->b < 0) {
            v += tm->w * a[tm->a];
            if (order >= 1) for (int ax = 0; ax < 3; ++ax) g[ax] += tm->w * ga[ax * nf + tm->a];
            continue;
        }
        v += tm->w * a[tm->a] * a[tm->b];
        if (order >= 1) {
            for (int ax = 0; ax < 3; ++ax) g[ax] += tm->w * (a[tm->a] * ga[ax * nf + tm->b] + a[tm->b] * ga[ax * nf + tm->a]);
        }
    }
    *val = v;
    if (grad) { grad[0] = g[0]; grad[1] = g[1]; grad[2] = g[2]; }
    if (buf != stack) free(buf);
}

// Upper bound of |x^i y^j z^k exp(-a r^2)| over all space: separable, each factor peaks at
// |t| = sqrt(i / 2a).
static double iso_cart_prim_max(int i, int j, int k, double a) {
    const int e[3] = {i, j, k};
    double m = 1.0;
    for (int ax = 0; ax < 3; ++ax) {
        if (e[ax] == 0) continue;
        const double t2 = e[ax] / (2.0 * a);
        m *= pow(t2, 0.5 * e[ax]) * exp(-0.5 * e[ax]);
    }
    return m;
}

// The basis as contracted shells for point evaluation, and max |phi_i| over all space per AO.
static bool iso_density_ao_shells(iso_density_t* d, const md_gto_basis_t* basis, const float* atom_xyz, double* out_ao_max) {
    d->num_ao_shells = (int)basis->num_shells;
    d->ao_shells = (iso_ao_shell_t*)calloc(basis->num_shells + 1, sizeof(iso_ao_shell_t));
    int ao = 0;
    for (uint32_t si = 0; si < basis->num_shells; ++si) {
        const md_gto_shell_t* sh = &basis->shells[si];
        iso_ao_shell_t* s = &d->ao_shells[si];
        if (sh->l > MD_GTO_MAX_ANGULAR_MOMENTUM || sh->num_primitives > ISO_MAX_PRIM) return false;
        s->l = (int)sh->l;
        s->ncart = (int)md_gto_num_cart_ao(sh->l);
        s->nprim = (int)sh->num_primitives;
        s->ao = ao;
        for (int k = 0; k < 3; ++k) s->c[k] = atom_xyz[3 * sh->atom_idx + k];
        double amin = DBL_MAX;
        for (int q = 0; q < s->nprim; ++q) {
            s->alpha[q] = basis->alpha[sh->primitive_offset + q];
            s->coef[q]  = basis->coeff[sh->primitive_offset + q];
            amin = fmin(amin, s->alpha[q]);
        }
        double nmax = 0.0;
        for (int n = 0; n < s->ncart; ++n) {
            int i, j, k;
            md_gto_cart_ijk(&i, &j, &k, sh->l, (uint32_t)n);
            s->ijk[n][0] = (int8_t)i; s->ijk[n][1] = (int8_t)j; s->ijk[n][2] = (int8_t)k;
            s->nrm[n] = md_gto_cart_norm_factor(i, j, k);
            nmax = fmax(nmax, fabs(s->nrm[n]));
            if (out_ao_max) {
                double m = 0.0;
                for (int q = 0; q < s->nprim; ++q) m += fabs(s->coef[q]) * iso_cart_prim_max(i, j, k, s->alpha[q]);
                out_ao_max[ao + n] = fabs(s->nrm[n]) * m;
            }
        }
        // past the radial maximum, r^l sum |c_q| exp(-a_q r^2) only decreases
        double r = sqrt(s->l / (2.0 * amin));
        for (;; r += 0.05) {
            double t = 0.0;
            for (int q = 0; q < s->nprim; ++q) t += fabs(s->coef[q]) * exp(-s->alpha[q] * r * r);
            if (nmax * pow(r, s->l) * t < ISO_EVAL_EPS * 1.0e-2) break;
        }
        s->r2_cut = r * r;
        ao += s->ncart;
    }
    d->num_ao = ao;
    return true;
}

// The contracted radial function of every basis shell as a group (ncart 1, l kept for the monomial
// powers), with its critical intervals: the AO range bounds use them.
static void iso_radial_init(iso_field_t* radial, const iso_ao_shell_t* shells, int n) {
    md_allocator_i* heap = md_get_heap_allocator();
    for (int si = 0; si < n; ++si) {
        const iso_ao_shell_t* a = &shells[si];
        iso_shell_t g = {0};
        g.l = a->l;
        g.ncart = 1;
        g.nprim = a->nprim;
        g.site = si;
        memcpy(g.c, a->c, sizeof(g.c));
        for (int q = 0; q < a->nprim; ++q) { g.alpha[q] = a->alpha[q]; g.w[0][q] = a->coef[q]; }
        md_array_push(radial->shells, g, heap);
    }
    radial->num_shells = n;
    for (int si = 0; si < n; ++si) iso_shell_finish(radial, &radial->shells[si], 1);
}

// Whether a density field also gets its functions as (site, function) groups, for the GROUP form of
// the bound (ISO_DENS_BOUND_GROUP). The AO form, which the octree uses, does not need them, and they
// are num_functions times the size of an orbital's; the tests turn them on to compare the two.
static bool iso_density_build_groups = false;

// A density field from num_fn functions (rows of num_ao AO coefficients) and terms over them.
// The groups are ordered by site, so a box is set up once per site and shared by its functions.
static bool iso_field_init_density(iso_field_t* f, const md_gto_basis_t* basis, const float* atom_xyz, size_t num_atoms,
                                   const double* rows, int num_fn, const iso_term_t* terms, int num_terms, double slack, int factor, const char* name) {
    memset(f, 0, sizeof(*f));
    md_snprintf(f->name, sizeof(f->name), "%s", name);
    md_allocator_i* heap = md_get_heap_allocator();
    f->num_input_shells = (int)basis->num_shells;
    for (int k = 0; k < 3; ++k) { f->atom_min[k] = DBL_MAX; f->atom_max[k] = -DBL_MAX; }
    for (size_t a = 0; a < num_atoms; ++a) {
        for (int k = 0; k < 3; ++k) {
            f->atom_min[k] = fmin(f->atom_min[k], atom_xyz[3 * a + k]);
            f->atom_max[k] = fmax(f->atom_max[k], atom_xyz[3 * a + k]);
        }
    }

    iso_density_t* d = (iso_density_t*)calloc(1, sizeof(iso_density_t));
    f->dens = d;
    d->num_fn = num_fn;
    d->num_terms = num_terms;
    d->slack = slack;
    d->factor = factor;
    d->terms = (iso_term_t*)malloc(sizeof(iso_term_t) * (num_terms + 1));
    memcpy(d->terms, terms, sizeof(iso_term_t) * num_terms);
    if (!iso_density_ao_shells(d, basis, atom_xyz, NULL)) return false;
    const int nao = d->num_ao;
    d->C = (double*)malloc(sizeof(double) * (size_t)nao * num_fn + 1);
    d->Cabs = (double*)malloc(sizeof(double) * (size_t)nao * num_fn + 1);
    for (int fi = 0; fi < num_fn; ++fi) {
        for (int i = 0; i < nao; ++i) {
            d->C[(size_t)i * num_fn + fi] = rows[(size_t)fi * nao + i];
            d->Cabs[(size_t)i * num_fn + fi] = fabs(rows[(size_t)fi * nao + i]);
        }
    }
    iso_radial_init(&d->radial, d->ao_shells, d->num_ao_shells);

    // sites: one per (atom, l), the union of the exponents of its shells
    const int L1 = MD_GTO_MAX_ANGULAR_MOMENTUM + 1;
    int* site_of = (int*)malloc(sizeof(int) * num_atoms * L1);
    for (size_t i = 0; i < num_atoms * L1; ++i) site_of[i] = -1;
    int* shell_site = (int*)malloc(sizeof(int) * (basis->num_shells + 1));
    iso_shell_t* sites = NULL;   // md_array: geometry and exponents only
    bool ok = true;
    for (uint32_t si = 0; si < basis->num_shells && ok; ++si) {
        const md_gto_shell_t* sh = &basis->shells[si];
        int* gi = &site_of[sh->atom_idx * L1 + sh->l];
        if (*gi < 0) {
            iso_shell_t g = {0};
            g.l = (int)sh->l;
            g.ncart = (int)md_gto_num_cart_ao(sh->l);
            for (int k = 0; k < 3; ++k) g.c[k] = atom_xyz[3 * sh->atom_idx + k];
            for (int n = 0; n < g.ncart; ++n) {
                int i, j, k;
                md_gto_cart_ijk(&i, &j, &k, sh->l, (uint32_t)n);
                g.ijk[n][0] = (int8_t)i; g.ijk[n][1] = (int8_t)j; g.ijk[n][2] = (int8_t)k;
            }
            *gi = (int)md_array_size(sites);
            g.site = *gi;
            md_array_push(sites, g, heap);
        }
        shell_site[si] = *gi;
        iso_shell_t* g = &sites[*gi];
        for (uint32_t p = 0; p < sh->num_primitives; ++p) {
            const double a = basis->alpha[sh->primitive_offset + p];
            int q = 0;
            while (q < g->nprim && g->alpha[q] != a) ++q;
            if (q == g->nprim) {
                if (g->nprim == ISO_MAX_PRIM) { ok = false; break; }
                g->alpha[g->nprim++] = a;
            }
        }
    }

    // one group per (site, function) with any weight
    const int num_sites = iso_density_build_groups ? (int)md_array_size(sites) : 0;
    for (int st = 0; st < num_sites && ok; ++st) {
        for (int fi = 0; fi < num_fn; ++fi) {
            iso_shell_t g = sites[st];
            g.fn = fi;
            bool any = false;
            int ao = 0;
            for (uint32_t si = 0; si < basis->num_shells; ++si) {
                const md_gto_shell_t* sh = &basis->shells[si];
                const int ncart = (int)md_gto_num_cart_ao(sh->l);
                if (shell_site[si] == st) {
                    for (uint32_t p = 0; p < sh->num_primitives; ++p) {
                        const double a = basis->alpha[sh->primitive_offset + p];
                        const double c = basis->coeff[sh->primitive_offset + p];
                        int q = 0;
                        while (g.alpha[q] != a) ++q;
                        for (int n = 0; n < ncart; ++n) {
                            const double x = d->C[(size_t)(ao + n) * num_fn + fi];
                            if (x == 0.0) continue;
                            g.w[n][q] += x * md_gto_cart_norm_factor(g.ijk[n][0], g.ijk[n][1], g.ijk[n][2]) * c;
                            any = true;
                        }
                    }
                }
                ao += ncart;
            }
            if (!any) continue;
            md_array_push(f->shells, g, heap);
        }
    }
    md_array_free(sites, heap);
    free(shell_site);
    free(site_of);
    if (!ok) return false;
    f->num_shells = (int)md_array_size(f->shells);
    // only g_0 is bounded for a density (the GROUP bound)
    for (int si = 0; si < f->num_shells; ++si) iso_shell_finish(f, &f->shells[si], 1);
    return true;
}

// Bound of |rho| at distance >= r from every atom (r past every radial maximum): every AO is at
// most nrm r^l sum_q |c_q| exp(-a_q r^2) there, and each function at most sum_i |C_i| times that.
//
// Asked at a few dozen r for one surface box, over every AO and function each time that is the AOs x
// functions of |C| read a few dozen times: for a protein's density (2147 orbitals over 9636 AOs) a
// third of a second, before every extraction. The shells' radial functions are those of the basis'
// few elements, though: the sum over each one's AOs of |nrm| |C| is taken once (iso_tail_init), and an
// r costs the distinct radial functions times the functions.
typedef struct iso_tail_t {
    int     nf, num;            // functions, distinct radial functions
    int*    radial;             // [num]: a shell of each (an index into radial.shells)
    double* G;                  // [num][nf]: sum over the shells of the function of sum over their AOs of |nrm| |C|
    double* fn_tail;            // [nf] scratch
} iso_tail_t;

static bool iso_radial_same(const iso_shell_t* a, const iso_shell_t* b) {
    if (a->l != b->l || a->ncart != b->ncart || a->nprim != b->nprim) return false;
    for (int q = 0; q < a->nprim; ++q) {
        if (a->alpha[q] != b->alpha[q]) return false;
        for (int n = 0; n < a->ncart; ++n) if (a->w[n][q] != b->w[n][q]) return false;
    }
    return true;
}

static bool iso_tail_init(iso_tail_t* t, const iso_density_t* d) {
    memset(t, 0, sizeof(*t));
    const int nf = d->num_fn, ns = d->num_ao_shells;
    t->nf = nf;
    int* group = (int*)malloc(sizeof(int) * (ns + 1));
    t->radial = (int*)malloc(sizeof(int) * (ns + 1));
    if (!group || !t->radial) { free(group); free(t->radial); t->radial = NULL; return false; }
    for (int si = 0; si < ns; ++si) {
        const iso_shell_t* r = &d->radial.shells[si];
        int u = 0;
        while (u < t->num && !iso_radial_same(&d->radial.shells[t->radial[u]], r)) ++u;
        if (u == t->num) t->radial[t->num++] = si;
        group[si] = u;
    }
    t->G = (double*)calloc((size_t)t->num * nf + 1, sizeof(double));
    t->fn_tail = (double*)malloc(sizeof(double) * (nf + 1));
    if (!t->G || !t->fn_tail) { free(group); free(t->radial); free(t->G); free(t->fn_tail); memset(t, 0, sizeof(*t)); return false; }
    for (int si = 0; si < ns; ++si) {
        const iso_ao_shell_t* a = &d->ao_shells[si];
        double* g = t->G + (size_t)group[si] * nf;
        for (int n = 0; n < a->ncart; ++n) {
            const double* ca = d->Cabs + (size_t)(a->ao + n) * nf;
            const double an = fabs(a->nrm[n]);
            for (int fi = 0; fi < nf; ++fi) g[fi] += ca[fi] * an;
        }
    }
    free(group);
    return true;
}

static void iso_tail_free(iso_tail_t* t) {
    free(t->radial);
    free(t->G);
    free(t->fn_tail);
    memset(t, 0, sizeof(*t));
}

static double iso_density_tail(const iso_field_t* f, double r, iso_tail_t* tl) {
    const iso_density_t* d = f->dens;
    const int nf = tl->nf;
    double* fn_tail = tl->fn_tail;
    for (int fi = 0; fi < nf; ++fi) fn_tail[fi] = 0.0;
    for (int u = 0; u < tl->num; ++u) {
        const double tu = iso_shell_tail(&d->radial.shells[tl->radial[u]], r);
        const double* g = tl->G + (size_t)u * nf;
        for (int fi = 0; fi < nf; ++fi) fn_tail[fi] += g[fi] * tu;
    }
    double t = d->slack;
    for (int k = 0; k < d->num_terms; ++k) t += fabs(d->terms[k].w) * fn_tail[d->terms[k].a] * (d->terms[k].b < 0 ? 1.0 : fn_tail[d->terms[k].b]);
    return t;
}

// ---------------------------------------------------------------------------------------------
// Splitting a symmetric matrix into signed squares: S = sum_k l_k c_k c_k^T + R
// ---------------------------------------------------------------------------------------------

// Pivoted LDL^T with 1x1 pivots, the largest |residual diagonal| first, until every residual diagonal
// is at rounding level. For a semidefinite S this is pivoted Cholesky. F[i * n + k] = c_k[i].
// Returns the number of factors, or -1 when the factors grow (an indefinite S that 1x1 pivots
// cannot handle, such as the symmetric part of a transition density with a zero diagonal block).
static int iso_factor_ldl(const double* S, int n, double* F, double* lam) {
    const double u = 0.5 * DBL_EPSILON;
    double dmax = 0.0, smax = 0.0;
    for (int i = 0; i < n; ++i) dmax = fmax(dmax, fabs(S[(size_t)i * n + i]));
    for (size_t i = 0; i < (size_t)n * n; ++i) smax = fmax(smax, fabs(S[i]));
    if (!(smax > 0.0)) return 0;
    const double tol = 64.0 * n * u * fmax(dmax, smax);
    double* d = (double*)malloc(sizeof(double) * n);
    double* w = (double*)malloc(sizeof(double) * n);
    bool* used = (bool*)calloc(n, sizeof(bool));
    for (int i = 0; i < n; ++i) d[i] = S[(size_t)i * n + i];
    int r = 0;
    bool grew = false;
    for (;;) {
        int p = -1;
        double best = tol;
        for (int i = 0; i < n; ++i) if (!used[i] && fabs(d[i]) > best) { best = fabs(d[i]); p = i; }
        if (p < 0) break;
        const double* fp = F + (size_t)p * n;
        for (int k = 0; k < r; ++k) w[k] = lam[k] * fp[k];
        double dp = S[(size_t)p * n + p];
        for (int k = 0; k < r; ++k) dp -= w[k] * fp[k];
        used[p] = true;
        d[p] = 0.0;
        if (!(fabs(dp) > tol)) continue;
        for (int i = 0; i < n; ++i) {
            double* fi = F + (size_t)i * n;
            if (i == p) { fi[r] = 1.0; continue; }
            if (used[i]) { fi[r] = 0.0; continue; }
            double v = S[(size_t)i * n + p];
            for (int k = 0; k < r; ++k) v -= w[k] * fi[k];
            fi[r] = v / dp;
            if (fabs(fi[r]) > 1.0e3) grew = true;
            d[i] -= dp * fi[r] * fi[r];
        }
        lam[r++] = dp;
        if (grew) break;
    }
    free(used); free(w); free(d);
    return grew ? -1 : r;
}

// Symmetric eigendecomposition: Householder tridiagonalization and the implicit QL algorithm
// (tred2 / tql2 as in EISPACK and JAMA). In: V = S (row major). Out: V's columns the eigenvectors,
// val the eigenvalues.
static void iso_eigen_sym(double* V, double* val, int n) {
#define ISO_V(i, j) V[(size_t)(i) * n + (j)]
    double* d = val;
    double* e = (double*)calloc(n, sizeof(double));
    for (int j = 0; j < n; ++j) d[j] = ISO_V(n - 1, j);
    for (int i = n - 1; i > 0; --i) {
        double scale = 0.0, h = 0.0;
        for (int k = 0; k < i; ++k) scale += fabs(d[k]);
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; ++j) { d[j] = ISO_V(i - 1, j); ISO_V(i, j) = 0.0; ISO_V(j, i) = 0.0; }
        } else {
            for (int k = 0; k < i; ++k) { d[k] /= scale; h += d[k] * d[k]; }
            double f = d[i - 1];
            double g = sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; ++j) e[j] = 0.0;
            for (int j = 0; j < i; ++j) {
                f = d[j];
                ISO_V(j, i) = f;
                g = e[j] + ISO_V(j, j) * f;
                for (int k = j + 1; k <= i - 1; ++k) { g += ISO_V(k, j) * d[k]; e[k] += ISO_V(k, j) * f; }
                e[j] = g;
            }
            f = 0.0;
            for (int j = 0; j < i; ++j) { e[j] /= h; f += e[j] * d[j]; }
            const double hh = f / (h + h);
            for (int j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (int j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; ++k) ISO_V(k, j) -= (f * e[k] + g * d[k]);
                d[j] = ISO_V(i - 1, j);
                ISO_V(i, j) = 0.0;
            }
        }
        d[i] = h;
    }
    for (int i = 0; i < n - 1; ++i) {
        ISO_V(n - 1, i) = ISO_V(i, i);
        ISO_V(i, i) = 1.0;
        const double h = d[i + 1];
        if (h != 0.0) {
            for (int k = 0; k <= i; ++k) d[k] = ISO_V(k, i + 1) / h;
            for (int j = 0; j <= i; ++j) {
                double g = 0.0;
                for (int k = 0; k <= i; ++k) g += ISO_V(k, i + 1) * ISO_V(k, j);
                for (int k = 0; k <= i; ++k) ISO_V(k, j) -= g * d[k];
            }
        }
        for (int k = 0; k <= i; ++k) ISO_V(k, i + 1) = 0.0;
    }
    for (int j = 0; j < n; ++j) { d[j] = ISO_V(n - 1, j); ISO_V(n - 1, j) = 0.0; }
    ISO_V(n - 1, n - 1) = 1.0;
    e[0] = 0.0;

    for (int i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    double f = 0.0, tst1 = 0.0;
    const double eps = DBL_EPSILON;
    for (int l = 0; l < n; ++l) {
        tst1 = fmax(tst1, fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps * tst1) break;
            m++;
        }
        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 60) break;
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (int i = l + 2; i < n; ++i) d[i] -= h;
                f += h;
                p = d[m];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0.0, s2 = 0.0;
                for (int i = m - 1; i >= l; --i) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (int k = 0; k < n; ++k) {
                        h = ISO_V(k, i + 1);
                        ISO_V(k, i + 1) = s * ISO_V(k, i) + c * h;
                        ISO_V(k, i) = c * ISO_V(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (fabs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0.0;
    }
    free(e);
#undef ISO_V
}

// max |R_ij| and sum_ij |R_ij| max|phi_i| max|phi_j| of R = S - sum_k l_k c_k c_k^T.
static double iso_factor_residual(const double* S, int n, const double* F, const double* lam, int r, const double* ao_max, double* out_rmax) {
    double slack = 0.0, rmax = 0.0;
    double* w = (double*)malloc(sizeof(double) * (r + 1));
    for (int i = 0; i < n; ++i) {
        const double* fi = F + (size_t)i * n;
        for (int k = 0; k < r; ++k) w[k] = lam[k] * fi[k];
        for (int j = 0; j < n; ++j) {
            const double* fj = F + (size_t)j * n;
            double v = S[(size_t)i * n + j];
            for (int k = 0; k < r; ++k) v -= w[k] * fj[k];
            rmax = fmax(rmax, fabs(v));
            slack += fabs(v) * ao_max[i] * ao_max[j];
        }
    }
    free(w);
    *out_rmax = rmax;
    return slack;
}

// Splits D (any n x n; its symmetric part is what a density sees) into r signed squares.
// rows: r x n factor rows (allocated), weights: r. Returns r (0: D is zero), and the rigorous
// residual bound sum_ij |R_ij| max|phi_i| max|phi_j| through out_slack.
static int iso_factor_density(const double* D, int n, const double* ao_max, double** out_rows, double** out_w, double* out_slack, int* out_factor) {
    double* S = (double*)malloc(sizeof(double) * (size_t)n * n);
    double smax = 0.0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        S[(size_t)i * n + j] = 0.5 * (D[(size_t)i * n + j] + D[(size_t)j * n + i]);
        smax = fmax(smax, fabs(S[(size_t)i * n + j]));
    }
    const double tol = 64.0 * n * 0.5 * DBL_EPSILON * smax;
    double* F = (double*)calloc((size_t)n * n, sizeof(double));   // F[i * n + k] = c_k[i]
    double* lam = (double*)calloc(n, sizeof(double));
    int factor = ISO_FACTOR_LDL;
    int r = iso_factor_ldl(S, n, F, lam);
    double rmax = 0.0, slack = 0.0;
    if (r >= 0) slack = iso_factor_residual(S, n, F, lam, r, ao_max, &rmax);
    // 1x1 pivots fail on growth, and stop early on a diagonal that is zero while the rest is not
    // (the symmetric part of a transition density between orbitals of different symmetry)
    if (r < 0 || rmax > 16.0 * tol) {
        // eigenvectors: always stable, but delocalized
        factor = ISO_FACTOR_EIGEN;
        memcpy(F, S, sizeof(double) * (size_t)n * n);
        double* ev = (double*)malloc(sizeof(double) * n);
        iso_eigen_sym(F, ev, n);
        double emax = 0.0;
        for (int k = 0; k < n; ++k) emax = fmax(emax, fabs(ev[k]));
        const double etol = 64.0 * n * 0.5 * DBL_EPSILON * emax;
        r = 0;
        for (int k = 0; k < n; ++k) {   // compact the significant columns
            if (!(fabs(ev[k]) > etol)) continue;
            for (int i = 0; i < n; ++i) F[(size_t)i * n + r] = F[(size_t)i * n + k];
            lam[r++] = ev[k];
        }
        free(ev);
        slack = iso_factor_residual(S, n, F, lam, r, ao_max, &rmax);
    }
    *out_rows = (double*)malloc(sizeof(double) * (size_t)(r + 1) * n);
    *out_w = (double*)malloc(sizeof(double) * (r + 1));
    for (int k = 0; k < r; ++k) {
        for (int i = 0; i < n; ++i) (*out_rows)[(size_t)k * n + i] = F[(size_t)i * n + k];
        (*out_w)[k] = lam[k];
    }
    *out_slack = slack * (1.0 + 1.0e-6) + DBL_MIN;
    *out_factor = factor;
    free(lam); free(F); free(S);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Octree: lazily built from the bounds, never materializing a volume
// ---------------------------------------------------------------------------------------------

#define ISO_MAX_LEVEL 16

enum { ISO_DENS_BOUND_AO = 0, ISO_DENS_BOUND_GROUP = 1 };   // the density bound the octree uses
static int iso_dens_bound_default = ISO_DENS_BOUND_AO;
#define ISO_SCHEME_CASCADE 100   // GROUP first (free with the drop pass), then T2 (small boxes) or PRIM (large)

typedef struct iso_leaf_t { int32_t x, y, z; uint8_t sheets; } iso_leaf_t;

typedef struct iso_octree_t {
    const iso_field_t* f;
    double   iso;             // the smallest iso value served (drop threshold, guarantee check)
    double   iso_lo, iso_hi;  // keep cells whose bound meets [iso_lo, iso_hi] or [-iso_hi, -iso_lo]
    double   h;               // leaf edge
    double   origin[3];
    int      L;               // leaf level; root edge = h * 2^L
    int      scheme;          // ISO_BOUND_* deciding pruning, or ISO_SCHEME_CASCADE
    double   t2_max_size;     // cascade: T2 for boxes up to this edge, PRIM above
    uint32_t mask;            // bounds computed
    double   drop_eps;        // groups bounded below this are dropped for the subtree
    int*     active[ISO_MAX_LEVEL + 1];
    // density fields: per function slack per level, and scratch for one node
    double*  fslack[ISO_MAX_LEVEL + 1];
    iv_t*    fiv;             // function intervals
    double*  fsens;           // d rho / d function, bounded over the node
    iv_t*    gint;            // group intervals of the node's active groups
    double*  aomag;           // AO form: max |phi_i| over the node, per AO
    double*  facc;            // AO form: accumulators, 5 per function
    int      dens_bound;      // ISO_DENS_BOUND_*
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
#define ISO_SHEETS(e) (uint8_t)((((e).lo <= o->iso_hi && (e).hi >= o->iso_lo) ? 1 : 0) | (((e).lo <= -o->iso_lo && (e).hi >= -o->iso_hi) ? 2 : 0))
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

// The same for a density: GROUP bounds summed per function, each function's interval widened by
// the bounds of its groups dropped higher up, then the terms as interval products (exact squares).
// A group is dropped for the subtree when its bound times the sensitivity of rho to its function
// (sum over its terms of |w| times the partner's magnitude) is below drop_eps.
// Bound of rho over the box from the active groups, each function widened by its carried slack.
// Scratch out: A (function intervals, num_fn), gint (per active group), S (sensitivities, num_fn).
static iv_t iso_density_bound(const iso_field_t* f, const int* act, int nact, const double* fslack, const iso_box_t* b,
                              iv_t* A, iv_t* gint, double* S) {
    const iso_density_t* d = f->dens;
    const int nf = d->num_fn;
    for (int j = 0; j < nf; ++j) A[j] = iv(-fslack[j], fslack[j]);
    iso_shell_box_t sb;
    int site = -1;
    for (int i = 0; i < nact; ++i) {
        const iso_shell_t* g = &f->shells[act[i]];
        if (g->site != site) { iso_shell_box(&sb, g, b); site = g->site; }
        const iv_t gi = iso_group_interval_sb(f, g, &sb);
        gint[i] = gi;
        A[g->fn] = iv_add(A[g->fn], gi);
    }
    for (int j = 0; j < nf; ++j) S[j] = 0.0;
    iv_t rho = iv_pt(0.0);
    double rmag = 0.0;
    for (int t = 0; t < d->num_terms; ++t) {
        const iso_term_t* tm = &d->terms[t];
        const double aw = fabs(tm->w);
        iv_t q;
        if (tm->b < 0) {
            q = A[tm->a];
            S[tm->a] += aw;
        } else if (tm->a == tm->b) {
            q = iv_pow(A[tm->a], 2);
            S[tm->a] += 2.0 * aw * iv_mag(A[tm->a]);
        } else {
            q = iv_mul(A[tm->a], A[tm->b]);
            S[tm->a] += aw * iv_mag(A[tm->b]);
            S[tm->b] += aw * iv_mag(A[tm->a]);
        }
        rho = iv_add(rho, iv_scale(q, tm->w));
        rmag += aw * iv_mag(q);
    }
    // rounding: the group sums carry a relative error far below 1e-12 of their magnitudes
    return iv_widen(rho, d->slack + 1.0e-12 * (rmag + iv_width(rho)));
}

// The AO form of the same bound (TOPO_ANALYTIC's idea of letting the coefficients multiply exact
// centre data): every active basis shell gives its AOs' range over the box, Phi_i, from its
// contracted radial function and the monomials, and its values at the box centre, phi_i(c). The
// functions then follow as two matrix-vector products with C and |C|, intersected:
//     a in  C mid(Phi)  +- |C| rad(Phi)                       (interval extension)
//     a in  C phi(c)    +- |C| max(Phi.hi - phi(c), phi(c) - Phi.lo)   (centred: cancellations of
//                                                               the centre values survive exactly)
// It sees no cancellation between the shells of one atom inside the radius, which the GROUP form
// does, but costs one pass over the active AOs instead of one per function and site.
// Scratch out: A (num_fn), aomag (num_ao: max |Phi_i|, for the drop pass), S (num_fn),
// acc (5 num_fn).
static iv_t iso_density_bound_ao(const iso_field_t* f, const int* act, int nact, const double* fslack, const iso_box_t* b,
                                 iv_t* A, double* aomag, double* S, double* acc) {
    const iso_density_t* d = f->dens;
    const int nf = d->num_fn;
    double* m1 = acc;
    double* r1 = acc + nf;
    double* m2 = acc + 2 * nf;
    double* r2 = acc + 3 * nf;
    double* mg = acc + 4 * nf;
    memset(acc, 0, sizeof(double) * 5 * nf);
    const double cen[3] = {0.5 * (b->lo[0] + b->hi[0]), 0.5 * (b->lo[1] + b->hi[1]), 0.5 * (b->lo[2] + b->hi[2])};
    for (int k = 0; k < nact; ++k) {
        const iso_ao_shell_t* a = &d->ao_shells[act[k]];
        const iso_shell_t* g = &d->radial.shells[act[k]];
        iso_shell_box_t sb;
        iso_shell_box(&sb, g, b);
        const iv_t R = iso_g_range_sb(&d->radial, g, &sb, 0, 0);
        const double dc[3] = {cen[0] - a->c[0], cen[1] - a->c[1], cen[2] - a->c[2]};
        const double rc = dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2];
        double Rc = 0.0;
        for (int q = 0; q < a->nprim; ++q) Rc += a->coef[q] * exp(-a->alpha[q] * rc);
        double px[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], py[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], pz[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
        iso_powers(px, dc[0], a->l);
        iso_powers(py, dc[1], a->l);
        iso_powers(pz, dc[2], a->l);
        for (int n = 0; n < a->ncart; ++n) {
            const int i = a->ijk[n][0], j = a->ijk[n][1], kk = a->ijk[n][2];
            const iv_t P = iv_scale(iv_mul(iso_mono(&sb, i, j, kk), R), a->nrm[n]);
            const double pc = a->nrm[n] * px[i] * py[j] * pz[kk] * Rc;
            const double mid = 0.5 * (P.lo + P.hi), rad = 0.5 * (P.hi - P.lo);
            const double dev = fmax(fmax(P.hi - pc, pc - P.lo), 0.0);
            const double mag = iv_mag(P);
            aomag[a->ao + n] = mag;
            const double* c  = d->C    + (size_t)(a->ao + n) * nf;
            const double* ca = d->Cabs + (size_t)(a->ao + n) * nf;
            for (int fi = 0; fi < nf; ++fi) {
                m1[fi] += c[fi] * mid;
                r1[fi] += ca[fi] * rad;
                m2[fi] += c[fi] * pc;
                r2[fi] += ca[fi] * dev;
                mg[fi] += ca[fi] * mag;
            }
        }
    }
    for (int fi = 0; fi < nf; ++fi) {
        const double e = fslack[fi] + 1.0e-12 * mg[fi];
        A[fi] = iv_isect(iv(m1[fi] - r1[fi] - e, m1[fi] + r1[fi] + e), iv(m2[fi] - r2[fi] - e, m2[fi] + r2[fi] + e));
    }
    for (int j = 0; j < nf; ++j) S[j] = 0.0;
    iv_t rho = iv_pt(0.0);
    double rmag = 0.0;
    for (int t = 0; t < d->num_terms; ++t) {
        const iso_term_t* tm = &d->terms[t];
        const double aw = fabs(tm->w);
        iv_t q;
        if (tm->b < 0) {
            q = A[tm->a];
            S[tm->a] += aw;
        } else if (tm->a == tm->b) {
            q = iv_pow(A[tm->a], 2);
            S[tm->a] += 2.0 * aw * iv_mag(A[tm->a]);
        } else {
            q = iv_mul(A[tm->a], A[tm->b]);
            S[tm->a] += aw * iv_mag(A[tm->b]);
            S[tm->b] += aw * iv_mag(A[tm->a]);
        }
        rho = iv_add(rho, iv_scale(q, tm->w));
        rmag += aw * iv_mag(q);
    }
    return iv_widen(rho, d->slack + 1.0e-12 * (rmag + iv_width(rho)));
}

static void iso_octree_visit_density_ao(iso_octree_t* o, int level, int ix, int iy, int iz, const int* act, int nact, const double* fslack) {
    const iso_field_t* f = o->f;
    const iso_density_t* d = f->dens;
    const int nf = d->num_fn;
    const double size = o->h * (double)(1 << (o->L - level));
    iso_box_t b;
    b.lo[0] = o->origin[0] + ix * size; b.hi[0] = b.lo[0] + size;
    b.lo[1] = o->origin[1] + iy * size; b.hi[1] = b.lo[1] + size;
    b.lo[2] = o->origin[2] + iz * size; b.hi[2] = b.lo[2] + size;
    o->visited[level]++;

    const iv_t rho = iso_density_bound_ao(f, act, nact, fslack, &b, o->fiv, o->aomag, o->fsens, o->facc);
    const double* S = o->fsens;
    o->group_work += nact;

    const uint8_t sheets = (uint8_t)(((rho.lo <= o->iso_hi && rho.hi >= o->iso_lo) ? 1 : 0) | ((rho.lo <= -o->iso_lo && rho.hi >= -o->iso_hi) ? 2 : 0));
    if (!sheets) return;
    o->kept[level]++;
    if (level == o->L) {
        iso_leaf_t leaf = {ix, iy, iz, sheets};
        md_array_push(o->leaves, leaf, md_get_heap_allocator());
        return;
    }

    // drop shells whose AOs move rho by less than drop_eps here; their bound becomes per-function slack
    int* next = o->active[level];
    double* ns = o->fslack[level];
    memcpy(ns, fslack, sizeof(double) * nf);
    int n = 0;
    for (int k = 0; k < nact; ++k) {
        const iso_ao_shell_t* a = &d->ao_shells[act[k]];
        double effect = 0.0;
        for (int c = 0; c < a->ncart; ++c) {
            const double* ca = d->Cabs + (size_t)(a->ao + c) * nf;
            double e = 0.0;
            for (int fi = 0; fi < nf; ++fi) e += ca[fi] * S[fi];
            effect += e * o->aomag[a->ao + c];
        }
        if (effect < o->drop_eps) {
            for (int c = 0; c < a->ncart; ++c) {
                const double* ca = d->Cabs + (size_t)(a->ao + c) * nf;
                const double m = o->aomag[a->ao + c];
                for (int fi = 0; fi < nf; ++fi) ns[fi] += ca[fi] * m;
            }
        } else {
            next[n++] = act[k];
        }
    }
    for (int c = 0; c < 8; ++c) {
        iso_octree_visit_density_ao(o, level + 1, 2 * ix + (c & 1), 2 * iy + ((c >> 1) & 1), 2 * iz + ((c >> 2) & 1), next, n, ns);
    }
}

static void iso_octree_visit_density(iso_octree_t* o, int level, int ix, int iy, int iz, const int* act, int nact, const double* fslack) {
    const iso_field_t* f = o->f;
    const int nf = f->dens->num_fn;
    const double size = o->h * (double)(1 << (o->L - level));
    iso_box_t b;
    b.lo[0] = o->origin[0] + ix * size; b.hi[0] = b.lo[0] + size;
    b.lo[1] = o->origin[1] + iy * size; b.hi[1] = b.lo[1] + size;
    b.lo[2] = o->origin[2] + iz * size; b.hi[2] = b.lo[2] + size;
    o->visited[level]++;

    const iv_t rho = iso_density_bound(f, act, nact, fslack, &b, o->fiv, o->gint, o->fsens);
    const double* S = o->fsens;
    o->group_work += nact;

    const uint8_t sheets = (uint8_t)(((rho.lo <= o->iso_hi && rho.hi >= o->iso_lo) ? 1 : 0) | ((rho.lo <= -o->iso_lo && rho.hi >= -o->iso_hi) ? 2 : 0));
    if (!sheets) return;
    o->kept[level]++;
    if (level == o->L) {
        iso_leaf_t leaf = {ix, iy, iz, sheets};
        md_array_push(o->leaves, leaf, md_get_heap_allocator());
        return;
    }

    int* next = o->active[level];
    double* ns = o->fslack[level];
    memcpy(ns, fslack, sizeof(double) * nf);
    int n = 0;
    for (int i = 0; i < nact; ++i) {
        const iso_shell_t* g = &f->shells[act[i]];
        const double m = iv_mag(o->gint[i]);
        if (m * S[g->fn] < o->drop_eps) {
            ns[g->fn] += m;
        } else {
            next[n++] = act[i];
        }
    }
    for (int c = 0; c < 8; ++c) {
        iso_octree_visit_density(o, level + 1, 2 * ix + (c & 1), 2 * iy + ((c >> 1) & 1), 2 * iz + ((c >> 2) & 1), next, n, ns);
    }
}

// The sum of all group tails at distance r (tl: a density's, iso_tail_init)
static double iso_surface_tail(const iso_field_t* f, double r, iso_tail_t* tl) {
    if (f->dens) return iso_density_tail(f, r, tl);
    double t = 0.0;
    for (int i = 0; i < f->num_shells; ++i) t += iso_shell_tail(&f->shells[i], r);
    return t;
}

// The atom box padded until the sum of all group tails is below iso/2: no surface exists outside.
static void iso_surface_box(const iso_field_t* f, double iso, double lo[3], double hi[3]) {
    double r0 = 0.0;
    for (int i = 0; i < f->num_shells; ++i) r0 = fmax(r0, iso_shell_tail_r0(&f->shells[i]));
    if (f->dens) for (int i = 0; i < f->dens->radial.num_shells; ++i) r0 = fmax(r0, iso_shell_tail_r0(&f->dens->radial.shells[i]));
    // the first pad = r0 + 0.05 k (k >= 0) with the tails below iso/2. Beyond r0 every tail decreases
    // (r^l exp(-a r^2) does for r^2 >= l / 2a), so k is found by doubling and bisection.
    iso_tail_t tl = {0};
    if (f->dens && !iso_tail_init(&tl, f->dens)) {
        // no memory for it: padded by what no Gaussian tail of a basis reaches past (exp(-0.05 * 20^2) = 2e-9)
        MD_LOG_ERROR("md_gto_iso: out of memory for the surface box");
        for (int a = 0; a < 3; ++a) { lo[a] = f->atom_min[a] - (r0 + 20.0); hi[a] = f->atom_max[a] + (r0 + 20.0); }
        return;
    }
#define ISO_TAIL_AT(k_) iso_surface_tail(f, r0 + 0.05 * (double)(k_), &tl)
    int k_hi = 0;
    while (ISO_TAIL_AT(k_hi) >= 0.5 * iso && r0 + 0.05 * k_hi <= 1.0e3) k_hi = k_hi ? 2 * k_hi : 1;
    int k_lo = k_hi / 2;    // tail(k_lo) >= iso/2 unless k_hi == 0
    if (k_hi > 0) {
        while (k_hi - k_lo > 1) {
            const int k = (k_lo + k_hi) / 2;
            if (ISO_TAIL_AT(k) < 0.5 * iso) k_hi = k; else k_lo = k;
        }
    }
#undef ISO_TAIL_AT
    const double pad = r0 + 0.05 * (double)k_hi;
    iso_tail_free(&tl);
    for (int a = 0; a < 3; ++a) { lo[a] = f->atom_min[a] - pad; hi[a] = f->atom_max[a] + pad; }
}

// The lattice is given in floats: its origin and its cell size are values a float holds. The GPU path
// forms its lattice points as origin + index * h in float, and an origin or h it cannot hold would move
// (or scale) its whole lattice against the one reported and the CPU's: by up to half an ulp of
// |origin| - the largest coordinate there is, the root's corner lying farthest from the molecule
// (6e-7 bohr at 16 bohr, 3e-4 of a core density's 45 au a tenth of a bohr from an oxygen). Held in
// floats, the two lattices differ only by the rounding of each point.
static inline double iso_lattice_float(double x) { return (double)(float)x; }

// Root over the box [lo, hi]: the cube of edge leaf_edge * 2^L centred on it. Returns L and the
// root's origin. The CPU and the GPU paths both take their root from here, so their lattices agree.
static int iso_root_box(const double lo[3], const double hi[3], double leaf_edge, double origin[3]) {
    double ext = 0.0;
    for (int a = 0; a < 3; ++a) ext = fmax(ext, hi[a] - lo[a]);
    int L = 0;
    while (leaf_edge * (double)(1 << L) < ext) L++;
    const double edge = leaf_edge * (double)(1 << L);
    for (int a = 0; a < 3; ++a) origin[a] = iso_lattice_float(0.5 * (lo[a] + hi[a]) - 0.5 * edge);
    return L;
}

// Root: the cube of edge leaf_edge * 2^L centred on the field's surface box.
static int iso_root(const iso_field_t* f, double iso, double leaf_edge, double origin[3]) {
    double lo[3], hi[3];
    iso_surface_box(f, iso, lo, hi);
    return iso_root_box(lo, hi, leaf_edge, origin);
}

static void iso_octree_build_range(iso_octree_t* o, const iso_field_t* f, double iso_lo, double iso_hi, double h, int scheme, uint32_t mask, double t2_max_size) {
    const double iso = iso_lo;
    memset(o, 0, sizeof(*o));
    o->f = f;
    o->iso = iso;
    o->iso_lo = iso_lo;
    o->iso_hi = iso_hi;
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
    if (f->dens) {
        const iso_density_t* d = f->dens;
        const int nf = d->num_fn;
        o->drop_eps = 1.0e-6 * iso;
        o->dens_bound = iso_dens_bound_default;
        for (int l = 0; l <= o->L; ++l) o->fslack[l] = (double*)malloc(sizeof(double) * (nf + 1));
        o->fiv   = (iv_t*)malloc(sizeof(iv_t) * (nf + 1));
        o->fsens = (double*)malloc(sizeof(double) * (nf + 1));
        o->gint  = (iv_t*)malloc(sizeof(iv_t) * (f->num_shells + 1));
        o->aomag = (double*)malloc(sizeof(double) * (d->num_ao + 1));
        o->facc  = (double*)malloc(sizeof(double) * 5 * (nf + 1));
        double* zero = (double*)calloc(nf + 1, sizeof(double));
        if (o->dens_bound == ISO_DENS_BOUND_AO || f->num_shells == 0) {
            for (int l = 0; l <= o->L; ++l) o->active[l] = (int*)realloc(o->active[l], sizeof(int) * (d->num_ao_shells + f->num_shells + 1));
            int* shells = (int*)malloc(sizeof(int) * (d->num_ao_shells + 1));
            for (int i = 0; i < d->num_ao_shells; ++i) shells[i] = i;
            iso_octree_visit_density_ao(o, 0, 0, 0, 0, shells, d->num_ao_shells, zero);
            free(shells);
        } else {
            iso_octree_visit_density(o, 0, 0, 0, 0, all, f->num_shells, zero);
        }
        free(zero);
    } else {
        iso_octree_visit(o, 0, 0, 0, 0, all, f->num_shells, 0.0);
    }
    o->seconds = iso_seconds_since(t0);
    free(all);
}

static void iso_octree_build_ex(iso_octree_t* o, const iso_field_t* f, double iso, double h, int scheme, uint32_t mask, double t2_max_size) {
    iso_octree_build_range(o, f, iso, iso, h, scheme, mask, t2_max_size);
}

static void iso_octree_build(iso_octree_t* o, const iso_field_t* f, double iso, double h, int scheme, uint32_t mask) {
    iso_octree_build_ex(o, f, iso, h, scheme, mask, 2.0);
}

static void iso_octree_free(iso_octree_t* o) {
    for (int l = 0; l <= ISO_MAX_LEVEL; ++l) { free(o->active[l]); free(o->fslack[l]); }
    free(o->fiv); free(o->fsens); free(o->gint); free(o->aomag); free(o->facc);
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

// Corner values of a set of 4^3 cell bricks, 125 per brick (faces shared between bricks are stored
// twice). The compact form a cache keeps: about 0.5 KB per brick.
typedef struct iso_brick_store_t {
    uint32_t  num_bricks;
    int32_t*  coord;          // 4 per brick: x, y, z in brick units, sheet mask
    float*    corner;         // 125 per brick, x fastest
    iso_map_t map;            // brick key -> index
} iso_brick_store_t;

// The brick lattice corner (x, y, z) is read from: of the up to 8 bricks holding it (a corner on a
// brick face is shared), the owner - the brick holding it at local 0..3 - when present, else the
// first present neighbour in a fixed order. Returns the brick index and the corner's local index
// there, or false. Only that brick has to evaluate the corner.
static bool iso_brick_store_owner(const iso_brick_store_t* s, int x, int y, int z, uint32_t* out_brick, int* out_local) {
    const int g[3] = {x, y, z};
    int b0[3], l0[3];
    for (int a = 0; a < 3; ++a) { b0[a] = g[a] >> 2; l0[a] = g[a] & 3; }
    for (int k = 0; k < 8; ++k) {
        int b[3], l[3];
        bool ok = true;
        for (int a = 0; a < 3 && ok; ++a) {
            const int alt = (k >> a) & 1;
            if (alt && l0[a] != 0) ok = false;
            b[a] = b0[a] - alt;
            l[a] = alt ? 4 : l0[a];
            if (b[a] < 0) ok = false;
        }
        if (!ok) continue;
        uint32_t bi;
        if (iso_map_get(&s->map, iso_key3(b[0], b[1], b[2]), &bi)) {
            *out_brick = bi;
            *out_local = l[0] + 5 * l[1] + 25 * l[2];
            return true;
        }
    }
    return false;
}

static bool iso_brick_store_corner(const iso_brick_store_t* s, int x, int y, int z, float* out) {
    uint32_t bi;
    int li;
    if (!iso_brick_store_owner(s, x, y, z, &bi, &li)) return false;
    *out = s->corner[125 * (size_t)bi + li];
    return true;
}

typedef struct iso_lattice_t {
    const iso_field_t* f;
    double     origin[3];
    double     h;
    iso_map_t  corner_map;
    double*    corner_val;    // md_array
    long       evals;
    const iso_brick_store_t* store;   // when set: corners come from the store, nothing is inserted
} iso_lattice_t;

static double iso_lattice_corner(iso_lattice_t* lat, int x, int y, int z) {
    if (lat->store) {
        float v;
        if (iso_brick_store_corner(lat->store, x, y, z, &v)) return v;
        const double p[3] = {lat->origin[0] + x * lat->h, lat->origin[1] + y * lat->h, lat->origin[2] + z * lat->h};
        lat->evals++;
        return iso_value(lat->f, p);   // outside every brick: no surface there, but stay correct
    }
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
    md_array_resize(m->nrm, (size_t)(3 * nv), heap);
    md_array_resize(m->res, (size_t)nv, heap);
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
// Corner values of one brick: the 5^3 lattice points of a 4^3 cell brick, x fastest
// ---------------------------------------------------------------------------------------------
//
// Both evaluators use separable exponentials, exp(-a |d|^2) = exp(-a dx^2) exp(-a dy^2) exp(-a dz^2):
// 15 per primitive and brick instead of 125, and skip what cannot reach the brick.

static inline double iso_box_dist2(const double c[3], const double lo[3], double edge) {
    double d2 = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double t = c[a] < lo[a] ? lo[a] - c[a] : (c[a] > lo[a] + edge ? c[a] - lo[a] - edge : 0.0);
        d2 += t * t;
    }
    return d2;
}

static void iso_orbital_brick(const iso_field_t* f, const double o[3], double h, const bool* need, double* out) {
    memset(out, 0, sizeof(double) * 125);
    const int L1 = MD_GTO_MAX_ANGULAR_MOMENTUM + 1;
    double P[3][5][MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
    double E[ISO_MAX_PRIM][3][5];
    double Eyz[ISO_MAX_PRIM][25];
    for (int si = 0; si < f->num_shells; ++si) {
        const iso_shell_t* s = &f->shells[si];
        if (iso_box_dist2(s->c, o, 4.0 * h) > s->r2_cut) continue;
        for (int a = 0; a < 3; ++a) {
            for (int k = 0; k < 5; ++k) {
                const double x = o[a] + k * h - s->c[a];
                P[a][k][0] = 1.0;
                for (int e = 1; e <= s->l; ++e) P[a][k][e] = P[a][k][e - 1] * x;
                for (int q = 0; q < s->nprim; ++q) E[q][a][k] = exp(-s->alpha[q] * x * x);
            }
        }
        for (int q = 0; q < s->nprim; ++q) for (int z = 0; z < 5; ++z) for (int y = 0; y < 5; ++y) Eyz[q][5 * z + y] = E[q][1][y] * E[q][2][z];
        for (int n = 0; n < s->ncart; ++n) {
            if (!s->live[n]) continue;
            const int i = s->ijk[n][0], j = s->ijk[n][1], k = s->ijk[n][2];
            for (int yz = 0; yz < 25; ++yz) {
                const int y = yz % 5, z = yz / 5;
                const double myz = P[1][y][j] * P[2][z][k];
                for (int x = 0; x < 5; ++x) {
                    if (!need[x + 5 * yz]) continue;
                    double R = 0.0;
                    for (int q = 0; q < s->nprim; ++q) R += s->w[n][q] * E[q][0][x] * Eyz[q][yz];
                    out[x + 5 * yz] += P[0][x][i] * myz * R;
                }
            }
        }
    }
    (void)L1;
}

// AO values at the brick's points (only shells that reach it, only the points it is asked for),
// the functions as one small matrix product with C, then the terms.
static void iso_density_brick(const iso_density_t* d, const double o[3], double h, const bool* need, double* out) {
    const int nf = d->num_fn;
    int pts[125], np = 0;
    for (int p = 0; p < 125; ++p) { out[p] = 0.0; if (need[p]) pts[np++] = p; }
    if (!np) return;
    double* phi = (double*)malloc(sizeof(double) * np * (d->num_ao + 1));   // [ao][np]
    int*    col = (int*)malloc(sizeof(int) * (d->num_ao + 1));
    double* a   = (double*)calloc((size_t)np * nf + 1, sizeof(double));      // [np][nf]
    int nc = 0;
    double P[3][5][MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
    double E[ISO_MAX_PRIM][3][5];
    double R[125];
    for (int si = 0; si < d->num_ao_shells; ++si) {
        const iso_ao_shell_t* s = &d->ao_shells[si];
        if (iso_box_dist2(s->c, o, 4.0 * h) > s->r2_cut) continue;
        for (int ax = 0; ax < 3; ++ax) {
            for (int k = 0; k < 5; ++k) {
                const double x = o[ax] + k * h - s->c[ax];
                P[ax][k][0] = 1.0;
                for (int e = 1; e <= s->l; ++e) P[ax][k][e] = P[ax][k][e - 1] * x;
                for (int q = 0; q < s->nprim; ++q) E[q][ax][k] = exp(-s->alpha[q] * x * x);
            }
        }
        for (int t = 0; t < np; ++t) {
            const int p = pts[t], x = p % 5, y = (p / 5) % 5, z = p / 25;
            double r = 0.0;
            for (int q = 0; q < s->nprim; ++q) r += s->coef[q] * E[q][0][x] * E[q][1][y] * E[q][2][z];
            R[t] = r;
        }
        for (int n = 0; n < s->ncart; ++n) {
            const int i = s->ijk[n][0], j = s->ijk[n][1], k = s->ijk[n][2];
            double* ph = phi + (size_t)nc * np;
            for (int t = 0; t < np; ++t) {
                const int p = pts[t], x = p % 5, y = (p / 5) % 5, z = p / 25;
                ph[t] = s->nrm[n] * P[0][x][i] * P[1][y][j] * P[2][z][k] * R[t];
            }
            col[nc++] = s->ao + n;
        }
    }
    for (int c = 0; c < nc; ++c) {
        const double* cr = d->C + (size_t)col[c] * nf;
        const double* ph = phi + (size_t)c * np;
        for (int t = 0; t < np; ++t) {
            const double v = ph[t];
            double* ap = a + (size_t)t * nf;
            for (int fi = 0; fi < nf; ++fi) ap[fi] += v * cr[fi];
        }
    }
    for (int t = 0; t < np; ++t) {
        const double* ap = a + (size_t)t * nf;
        double v = 0.0;
        for (int k = 0; k < d->num_terms; ++k) v += d->terms[k].w * ap[d->terms[k].a] * (d->terms[k].b < 0 ? 1.0 : ap[d->terms[k].b]);
        out[pts[t]] = v;
    }
    free(a); free(col); free(phi);
}

static void iso_field_brick(const iso_field_t* f, const double o[3], double h, const bool* need, double* out) {
    if (f->dens) iso_density_brick(f->dens, o, h, need, out);
    else         iso_orbital_brick(f, o, h, need, out);
}

// =============================================================================================
// Public API
// =============================================================================================

struct md_gto_iso_field_t {
    iso_field_t            f;
    struct md_allocator_i* alloc;
    double*                coeff;      // an orbital's AO coefficients (the GPU path takes them as a function)
    size_t                 num_coeff;
    float*                 dmat;       // a density matrix's symmetric part [num_ao][num_ao] (the GPU path's phi^T D phi)
    size_t                 dmat_dim;
    bool                   maxima_at_atoms;   // md_gto_iso_field_set_maxima_at_atoms
};

static float* iso_pack_atoms(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride, size_t* out_num_atoms) {
    const size_t num_atoms = md_gto_basis_num_atoms(basis);
    const size_t stride = atom_xyz_stride ? atom_xyz_stride : sizeof(float) * 3;
    float* xyz = (float*)malloc(sizeof(float) * 3 * (num_atoms + 1));
    for (size_t i = 0; i < num_atoms; ++i) {
        const float* p = (const float*)((const uint8_t*)atom_xyz + i * stride);
        xyz[3 * i + 0] = p[0];
        xyz[3 * i + 1] = p[1];
        xyz[3 * i + 2] = p[2];
    }
    *out_num_atoms = num_atoms;
    return xyz;
}

static md_gto_iso_field_t* iso_field_wrap(iso_field_t* f, bool ok, struct md_allocator_i* alloc, const char* who) {
    if (!ok) {
        MD_LOG_ERROR("%s: unsupported basis (angular momentum or contraction length)", who);
        iso_field_free(f);
        return NULL;
    }
    md_gto_iso_field_t* field = (md_gto_iso_field_t*)md_alloc(alloc, sizeof(md_gto_iso_field_t));
    memset(field, 0, sizeof(*field));
    field->f = *f;
    field->alloc = alloc;
    return field;
}

md_gto_iso_field_t* md_gto_iso_field_create(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                            const double* ao_coeffs, struct md_allocator_i* alloc) {
    if (!basis || !atom_xyz || !ao_coeffs || !alloc) {
        MD_LOG_ERROR("md_gto_iso_field_create: missing argument");
        return NULL;
    }
    size_t num_atoms;
    float* xyz = iso_pack_atoms(basis, atom_xyz, atom_xyz_stride, &num_atoms);
    iso_field_t f;
    const bool ok = iso_field_init(&f, basis, xyz, num_atoms, ao_coeffs, "");
    free(xyz);
    md_gto_iso_field_t* field = iso_field_wrap(&f, ok, alloc, "md_gto_iso_field_create");
    if (field) {
        field->num_coeff = md_gto_basis_num_ao(basis);
        field->coeff = (double*)md_alloc(alloc, sizeof(double) * (field->num_coeff + 1));
        memcpy(field->coeff, ao_coeffs, sizeof(double) * field->num_coeff);
    }
    return field;
}

// rows: num_fn x num_ao functions; terms over them.
static md_gto_iso_field_t* iso_density_create(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                              const double* rows, int num_fn, const iso_term_t* terms, int num_terms,
                                              double slack, int factor, struct md_allocator_i* alloc, const char* who) {
    size_t num_atoms;
    float* xyz = iso_pack_atoms(basis, atom_xyz, atom_xyz_stride, &num_atoms);
    iso_field_t f;
    const bool ok = iso_field_init_density(&f, basis, xyz, num_atoms, rows, num_fn, terms, num_terms, slack, factor, "");
    free(xyz);
    return iso_field_wrap(&f, ok, alloc, who);
}

md_gto_iso_field_t* md_gto_iso_field_create_products(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                     const double* left, const double* right, const double* weight, size_t count,
                                                     struct md_allocator_i* alloc) {
    if (!basis || !atom_xyz || !left || !weight || !alloc || count == 0) {
        MD_LOG_ERROR("md_gto_iso_field_create_products: missing argument");
        return NULL;
    }
    const size_t nao = md_gto_basis_num_ao(basis);
    const int nf = (int)(right ? 2 * count : count);
    double* rows = (double*)malloc(sizeof(double) * nao * nf);
    iso_term_t* terms = (iso_term_t*)malloc(sizeof(iso_term_t) * count);
    int n = 0, nt = 0;
    for (size_t k = 0; k < count; ++k) {
        if (weight[k] == 0.0) continue;
        memcpy(rows + (size_t)n * nao, left + k * nao, sizeof(double) * nao);
        const int a = n++;
        int b = a;
        if (right) {
            memcpy(rows + (size_t)n * nao, right + k * nao, sizeof(double) * nao);
            b = n++;
        }
        terms[nt++] = (iso_term_t){weight[k], a, b};
    }
    md_gto_iso_field_t* field = NULL;
    if (nt == 0) {
        MD_LOG_ERROR("md_gto_iso_field_create_products: every weight is zero");
    } else {
        field = iso_density_create(basis, atom_xyz, atom_xyz_stride, rows, n, terms, nt, 0.0, ISO_FACTOR_NONE, alloc, "md_gto_iso_field_create_products");
    }
    free(terms);
    free(rows);
    return field;
}

md_gto_iso_field_t* md_gto_iso_field_create_density(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                    const double* density_matrix, struct md_allocator_i* alloc) {
    if (!basis || !atom_xyz || !density_matrix || !alloc) {
        MD_LOG_ERROR("md_gto_iso_field_create_density: missing argument");
        return NULL;
    }
    size_t num_atoms;
    float* xyz = iso_pack_atoms(basis, atom_xyz, atom_xyz_stride, &num_atoms);
    iso_density_t tmp = {0};
    const int nao = (int)md_gto_basis_num_ao(basis);
    double* ao_max = (double*)malloc(sizeof(double) * (nao + 1));
    const bool ok = iso_density_ao_shells(&tmp, basis, xyz, ao_max);
    free(tmp.ao_shells);
    free(xyz);
    if (!ok) {
        free(ao_max);
        MD_LOG_ERROR("md_gto_iso_field_create_density: unsupported basis (angular momentum or contraction length)");
        return NULL;
    }
    double *rows = NULL, *w = NULL, slack = 0.0;
    int factor = 0;
    const int r = iso_factor_density(density_matrix, nao, ao_max, &rows, &w, &slack, &factor);
    free(ao_max);
    md_gto_iso_field_t* field = NULL;
    if (r <= 0) {
        MD_LOG_ERROR("md_gto_iso_field_create_density: the density matrix is zero");
    } else {
        iso_term_t* terms = (iso_term_t*)malloc(sizeof(iso_term_t) * (size_t)r);
        for (int k = 0; k < r; ++k) terms[k] = (iso_term_t){w[k], k, k};
        field = iso_density_create(basis, atom_xyz, atom_xyz_stride, rows, r, terms, r, slack, factor, alloc, "md_gto_iso_field_create_density");
        free(terms);
        // the matrix itself as well (symmetric part, float): the GPU path evaluates phi^T D phi from it
        if (field) {
            field->dmat = (float*)md_alloc(alloc, sizeof(float) * ((size_t)nao * nao + 1));
            if (field->dmat) {
                field->dmat_dim = (size_t)nao;
                for (int i = 0; i < nao; ++i) for (int j = 0; j < nao; ++j)
                    field->dmat[(size_t)i * nao + j] = (float)(0.5 * (density_matrix[(size_t)i * nao + j] + density_matrix[(size_t)j * nao + i]));
            }
        }
    }
    free(rows);
    free(w);
    return field;
}

void md_gto_iso_field_destroy(md_gto_iso_field_t* field) {
    if (!field) return;
    if (field->coeff) md_free(field->alloc, field->coeff, sizeof(double) * (field->num_coeff + 1));
    if (field->dmat) md_free(field->alloc, field->dmat, sizeof(float) * (field->dmat_dim * field->dmat_dim + 1));
    iso_field_free(&field->f);
    md_free(field->alloc, field, sizeof(*field));
}

void md_gto_iso_field_set_maxima_at_atoms(md_gto_iso_field_t* field, bool maxima_at_atoms) {
    if (field) field->maxima_at_atoms = maxima_at_atoms;
}

void md_gto_iso_field_info(md_gto_iso_field_info_t* out, const md_gto_iso_field_t* field) {
    if (!out) return;
    memset(out, 0, sizeof(*out));
    if (!field) return;
    const iso_field_t* f = &field->f;
    out->num_groups = (uint32_t)f->num_shells;
    out->memory_bytes = sizeof(*field) + md_array_size(f->shells) * sizeof(iso_shell_t) + md_array_size(f->crit) * sizeof(iso_crit_t);
    if (f->dens) {
        const iso_density_t* d = f->dens;
        out->num_groups = (uint32_t)d->num_ao_shells;
        out->memory_bytes += md_array_size(d->radial.shells) * sizeof(iso_shell_t) + md_array_size(d->radial.crit) * sizeof(iso_crit_t);
        out->kind = MD_GTO_ISO_FIELD_DENSITY;
        out->split = d->factor == ISO_FACTOR_LDL ? MD_GTO_ISO_SPLIT_LDL : d->factor == ISO_FACTOR_EIGEN ? MD_GTO_ISO_SPLIT_EIGEN : MD_GTO_ISO_SPLIT_NONE;
        out->num_functions = (uint32_t)d->num_fn;
        out->num_terms = (uint32_t)d->num_terms;
        out->residual = d->slack;
        out->memory_bytes += sizeof(*d) + sizeof(iso_term_t) * d->num_terms + 2 * sizeof(double) * (size_t)d->num_ao * d->num_fn
                           + sizeof(iso_ao_shell_t) * d->num_ao_shells;
    } else {
        out->kind = MD_GTO_ISO_FIELD_ORBITAL;
        out->num_functions = 1;
    }
}

void md_gto_iso_field_evaluate(const md_gto_iso_field_t* field, const float* xyz, size_t count, size_t stride, double* out_psi, double* out_grad) {
    if (!field || !xyz) return;
    const size_t st = stride ? stride : sizeof(float) * 3;
    for (size_t i = 0; i < count; ++i) {
        const float* q = (const float*)((const uint8_t*)xyz + i * st);
        const double p[3] = {q[0], q[1], q[2]};
        double v, g[3];
        iso_eval(&field->f, NULL, 0, p, out_grad ? 1 : 0, &v, g, NULL);
        if (out_psi) out_psi[i] = v;
        if (out_grad) { out_grad[3 * i] = g[0]; out_grad[3 * i + 1] = g[1]; out_grad[3 * i + 2] = g[2]; }
    }
}

struct md_gto_iso_cache_t {
    const iso_field_t*     f;
    struct md_allocator_i* alloc;
    double                 h, iso_min, iso_max;
    double                 origin[3];
    int                    brick_levels;   // octree levels from the root to the bricks
    iso_brick_store_t      store;
    double                 seconds_octree, seconds_corners;
};

md_gto_iso_cache_t* md_gto_iso_cache_create(const md_gto_iso_field_t* field, const md_gto_iso_cache_desc_t* desc, struct md_allocator_i* alloc) {
    if (!field || !desc || !alloc) {
        MD_LOG_ERROR("md_gto_iso_cache_create: missing argument");
        return NULL;
    }
    const double h = iso_lattice_float(desc->cell_size);
    const double iso_min = desc->iso_min;
    const double iso_max = desc->iso_max > iso_min ? desc->iso_max : iso_min;
    if (!(h > 0.0) || !(iso_min > 0.0)) {
        MD_LOG_ERROR("md_gto_iso_cache_create: cell_size and iso_min must be positive");
        return NULL;
    }

    iso_octree_t o;
    iso_octree_build_range(&o, &field->f, iso_min, iso_max, 4.0 * h, ISO_BOUND_SHELL, 0, 0.0);
    if (o.L + 2 > 20) {   // lattice coordinates are packed in 21 bits
        iso_octree_free(&o);
        MD_LOG_ERROR("md_gto_iso_cache_create: lattice too large for cell size %g", h);
        return NULL;
    }

    md_gto_iso_cache_t* c = (md_gto_iso_cache_t*)md_alloc(alloc, sizeof(md_gto_iso_cache_t));
    memset(c, 0, sizeof(*c));
    c->f = &field->f;
    c->alloc = alloc;
    c->h = h;
    c->iso_min = iso_min;
    c->iso_max = iso_max;
    c->brick_levels = o.L;
    c->seconds_octree = o.seconds;
    memcpy(c->origin, o.origin, sizeof(c->origin));

    const uint32_t nb = (uint32_t)md_array_size(o.leaves);
    iso_brick_store_t* s = &c->store;
    s->num_bricks = nb;
    s->coord  = (int32_t*)md_alloc(alloc, sizeof(int32_t) * 4 * (nb + 1));
    s->corner = (float*)md_alloc(alloc, sizeof(float) * 125 * (nb + 1));
    iso_map_init(&s->map, nb + 1);
    for (uint32_t b = 0; b < nb; ++b) {
        const iso_leaf_t* l = &o.leaves[b];
        s->coord[4 * b + 0] = l->x;
        s->coord[4 * b + 1] = l->y;
        s->coord[4 * b + 2] = l->z;
        s->coord[4 * b + 3] = l->sheets;
        iso_map_put(&s->map, iso_key3(l->x, l->y, l->z), b);
    }
    iso_octree_free(&o);
    if (!desc->defer_corners) {
        md_tick_t t_corners = md_tick_now();
        md_gto_iso_cache_evaluate(c, 0, nb);
        c->seconds_corners = iso_seconds_since(t_corners);
    }
    return c;
}

void md_gto_iso_cache_evaluate(md_gto_iso_cache_t* c, uint32_t first_brick, uint32_t num_bricks) {
    if (!c) return;
    const iso_brick_store_t* s = &c->store;
    const uint32_t end = first_brick + num_bricks < s->num_bricks ? first_brick + num_bricks : s->num_bricks;
    double v[125];
    bool need[125];
    for (uint32_t b = first_brick; b < end; ++b) {
        const int32_t* bc = s->coord + 4 * b;
        // only the corners lookups read from this brick: its own 4^3, and the shared ones on its far
        // faces whose owners are absent (the near-face corners of a present neighbour are its own)
        for (int k = 0; k < 125; ++k) {
            const int lx = k % 5, ly = (k / 5) % 5, lz = k / 25;
            if (lx < 4 && ly < 4 && lz < 4) { need[k] = true; continue; }
            uint32_t ob;
            int ol;
            need[k] = iso_brick_store_owner(s, bc[0] * 4 + lx, bc[1] * 4 + ly, bc[2] * 4 + lz, &ob, &ol) && ob == b && ol == k;
        }
        const double o[3] = {c->origin[0] + bc[0] * 4.0 * c->h, c->origin[1] + bc[1] * 4.0 * c->h, c->origin[2] + bc[2] * 4.0 * c->h};
        iso_field_brick(c->f, o, c->h, need, v);
        float* dst = s->corner + 125 * (size_t)b;
        for (int k = 0; k < 125; ++k) dst[k] = need[k] ? (float)v[k] : NAN;
    }
}

void md_gto_iso_cache_destroy(md_gto_iso_cache_t* c) {
    if (!c) return;
    md_free(c->alloc, c->store.coord, sizeof(int32_t) * 4 * (c->store.num_bricks + 1));
    md_free(c->alloc, c->store.corner, sizeof(float) * 125 * (c->store.num_bricks + 1));
    iso_map_free(&c->store.map);
    md_free(c->alloc, c, sizeof(*c));
}

bool md_gto_iso_cache_covers(const md_gto_iso_cache_t* c, double iso) {
    return c && iso >= c->iso_min && iso <= c->iso_max;
}

void md_gto_iso_cache_info(md_gto_iso_cache_info_t* out, const md_gto_iso_cache_t* c) {
    if (!out) return;
    memset(out, 0, sizeof(*out));
    if (!c) return;
    out->cell_size = c->h;
    out->iso_min = c->iso_min;
    out->iso_max = c->iso_max;
    memcpy(out->origin, c->origin, sizeof(out->origin));
    out->num_bricks = c->store.num_bricks;
    out->cells_per_axis = 4u << c->brick_levels;
    out->memory_bytes = sizeof(*c) + (size_t)c->store.num_bricks * (sizeof(int32_t) * 4 + sizeof(float) * 125)
                      + (size_t)c->store.map.cap * (sizeof(uint64_t) + sizeof(uint32_t));
    out->seconds_octree = c->seconds_octree;
    out->seconds_corners = c->seconds_corners;
}

bool md_gto_iso_mesh_extract(md_gto_iso_mesh_t* out, const md_gto_iso_cache_t* c, const md_gto_iso_mesh_desc_t* desc, struct md_allocator_i* alloc) {
    if (!out || !c || !desc || !alloc) {
        MD_LOG_ERROR("md_gto_iso_mesh_extract: missing argument");
        return false;
    }
    memset(out, 0, sizeof(*out));
    out->alloc = alloc;
    const double iso = desc->iso;
    if (!md_gto_iso_cache_covers(c, iso)) {
        MD_LOG_ERROR("md_gto_iso_mesh_extract: iso %g outside the cache's [%g, %g]", iso, c->iso_min, c->iso_max);
        return false;
    }
    const uint32_t sheets = desc->sheets ? desc->sheets : MD_GTO_ISO_SHEET_BOTH;

    // the cells of the bricks; their sheet masks say which sheets a cell may hold
    const uint32_t nb = c->store.num_bricks;
    iso_leaf_t* cells = (iso_leaf_t*)malloc(sizeof(iso_leaf_t) * 64 * (size_t)(nb + 1));
    int nc = 0;
    for (uint32_t b = 0; b < nb; ++b) {
        const int32_t* bc = c->store.coord + 4 * b;
        for (int k = 0; k < 64; ++k) {
            cells[nc++] = (iso_leaf_t){bc[0] * 4 + (k & 3), bc[1] * 4 + ((k >> 2) & 3), bc[2] * 4 + (k >> 4), (uint8_t)bc[3]};
        }
    }

    iso_lattice_t lat = {.f = c->f, .h = c->h, .store = &c->store};
    memcpy(lat.origin, c->origin, sizeof(lat.origin));

    iso_mesh_t m[2];
    memset(m, 0, sizeof(m));
    for (int s = 0; s < 2; ++s) {
        if (!(sheets & (1u << s))) continue;
        iso_surface_nets(&m[s], &lat, cells, nc, s == 0 ? 1 : -1, iso);
        if (iso_mesh_nv(&m[s]) == 0) continue;
        iso_mesh_project(&m[s], c->f, (int)desc->project_iterations);   // 0: normals only
        iso_mesh_triangulate(&m[s]);
        if (desc->relax_iterations) iso_mesh_relax(&m[s], c->f, (int)desc->relax_iterations, 0.5);
        for (uint32_t p = 0; p < desc->refine_passes; ++p) {
            if (!iso_mesh_refine(&m[s], c->f, 0.02 * c->h, 20.0)) break;
        }
    }
    free(cells);

    uint32_t nv = 0, ni = 0;
    for (int s = 0; s < 2; ++s) {
        out->sheet_vertex_offset[s] = nv;
        out->sheet_vertex_count[s]  = (uint32_t)iso_mesh_nv(&m[s]);
        out->sheet_index_offset[s]  = ni;
        out->sheet_index_count[s]   = (uint32_t)md_array_size(m[s].tri);
        nv += out->sheet_vertex_count[s];
        ni += out->sheet_index_count[s];
    }
    out->vertex_count = nv;
    out->index_count  = ni;
    out->position = (float*)md_alloc(alloc, sizeof(float) * 3 * (nv + 1));
    out->normal   = (float*)md_alloc(alloc, sizeof(float) * 3 * (nv + 1));
    out->index    = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * (ni + 1));
    for (int s = 0; s < 2; ++s) {
        const uint32_t vo = out->sheet_vertex_offset[s], io = out->sheet_index_offset[s];
        for (uint32_t i = 0; i < out->sheet_vertex_count[s]; ++i) {
            for (int a = 0; a < 3; ++a) {
                out->position[3 * (vo + i) + a] = (float)m[s].pos[3 * i + a];
                out->normal[3 * (vo + i) + a]   = m[s].nrm ? (float)m[s].nrm[3 * i + a] : 0.0f;
            }
        }
        for (uint32_t i = 0; i < out->sheet_index_count[s]; ++i) out->index[io + i] = vo + m[s].tri[i];
        iso_mesh_free(&m[s]);
    }
    return true;
}

void md_gto_iso_mesh_free(md_gto_iso_mesh_t* mesh) {
    if (!mesh || !mesh->alloc) return;
    md_free(mesh->alloc, mesh->position, sizeof(float) * 3 * (mesh->vertex_count + 1));
    md_free(mesh->alloc, mesh->normal, sizeof(float) * 3 * (mesh->vertex_count + 1));
    md_free(mesh->alloc, mesh->index, sizeof(uint32_t) * (mesh->index_count + 1));
    memset(mesh, 0, sizeof(*mesh));
}

#if MD_ENABLE_GPU
#include "md_gto_iso_gpu.inl"
#endif
