#include <md_mbis.h>

#include <md_gto.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_executor.h>
#include <core/md_log.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MBIS_DEFAULT_RADIAL     48
#define MBIS_DEFAULT_ANGULAR    35
#define MBIS_DEFAULT_THRESHOLD  1.0e-8
#define MBIS_DEFAULT_MAX_ITER   1000

// Bound of |phi| beyond which a shell is not evaluated at a point.
#define MBIS_AO_EPS             1.0e-12
// A point whose density is below this is dropped: an atom's grid reaches ~14 bohr (48 radial points), a sphere
// of ~1e4 bohr^3, so at most ~1e-10 e an atom goes.
#define MBIS_RHO_FLOOR          1.0e-14
// A pro-atom shell is left out at a point where it is below this fraction of the pro-molecule density
// there: every weight then moves by at most that, relative.
#define MBIS_TERM_EPS           1.0e-12
// Stratmann, Scuseria and Frisch, Chem. Phys. Lett. 257, 213 (1996): the cell function is exactly 1 or
// 0 outside |mu| < a, which is what makes the screening below exact.
#define MBIS_SSF_A              0.64

#define MBIS_MAX_CHUNKS         64
#define MBIS_ANDERSON_DEPTH     5

// Radial points closer to their nucleus than this fraction of the distance to the nearest atom use an angular
// grid of at most MBIS_INNER_ANGULAR (MBIS_INNER_LONE_ATOM bohr stands in for that distance for a lone atom).
// Measured on H to S: the charges move by less than 1e-6 e against the full grid everywhere, at ~40% fewer points.
#define MBIS_INNER_FRACTION     0.25
#define MBIS_INNER_ANGULAR      11
#define MBIS_INNER_LONE_ATOM    4.0

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

typedef struct mbis_sphere_t {
    double*  u;             // [n][3] unit vectors, band by band (one polar angle each)
    double*  w;             // [n], summing to 4 pi
    uint32_t n;
    uint32_t num_bands;
    uint32_t band_size;
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

// Gauss-Legendre nodes and weights on [-1, 1].
static void mbis_gauss_legendre(double* x, double* w, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        double z = cos(3.14159265358979323846 * (i + 0.75) / (n + 0.5));
        double dp = 1.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = 0.0;
            for (uint32_t k = 1; k <= n; ++k) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p2) / k;
            }
            dp = n * (z * p0 - p1) / (z * z - 1.0);
            const double dz = p0 / dp;
            z -= dz;
            if (fabs(dz) < 1.0e-15) break;
        }
        x[i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
}

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

// Gauss-Legendre in cos(theta) times a uniform grid in phi: exact for the spherical harmonics up to degree
// 'order' (odd), with (order + 1)^2 / 2 points.
static bool mbis_sphere_init(mbis_sphere_t* sph, uint32_t order) {
    const uint32_t n_theta = (order + 2) / 2;
    const uint32_t n_phi = 2 * n_theta;
    double gx[256], gw[256];
    if (n_theta > MD_COUNTOF(gx)) return false;
    sph->n = n_theta * n_phi;
    sph->num_bands = n_theta;
    sph->band_size = n_phi;
    sph->u = (double*)malloc(sizeof(double) * 3 * sph->n);
    sph->w = (double*)malloc(sizeof(double) * sph->n);
    if (!sph->u || !sph->w) return false;
    mbis_gauss_legendre(gx, gw, n_theta);
    const double two_pi = 2.0 * 3.14159265358979323846;
    uint32_t j = 0;
    for (uint32_t it = 0; it < n_theta; ++it) {
        const double ct = gx[it], st = sqrt(MD_MAX(0.0, 1.0 - ct * ct));
        for (uint32_t ip = 0; ip < n_phi; ++ip, ++j) {
            const double ph = two_pi * (ip + 0.5) / n_phi;
            sph->u[3 * j + 0] = st * cos(ph);
            sph->u[3 * j + 1] = st * sin(ph);
            sph->u[3 * j + 2] = ct;
            sph->w[j] = gw[it] * two_pi / n_phi;
        }
    }
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

// Weight of atom B's cell at p, which lies r from B. 'nb' is every atom sorted by distance from B (B first).
// With k = (1+a)/(1-a): s(mu_DC) = 1 once r_C >= k r_D, and the cell of D vanishes once r_D >= k r_min, so
// only atoms within k^2 r_min of p enter - all of them within r + k^2 r_min of B, which bounds the walk down
// the list as r_min shrinks - and B's own cell is 0 as soon as one atom is k times closer than B.
static double mbis_becke_weight(const double* p, double r, const mbis_neighbor_t* nb, size_t num_nb, const dvec3_t* xyz,
                                uint32_t* cand, double* cand_r) {
    const double k = (1.0 + MBIS_SSF_A) / (1.0 - MBIS_SSF_A);
    double rmin = r;
    double reach = r + k * k * rmin;
    size_t n = 0;
    for (size_t i = 0; i < num_nb && nb[i].d <= reach; ++i) {
        const double rc = i == 0 ? r : mbis_dist(p, xyz[nb[i].atom].elem);
        if (rc * k <= r) return 0.0;
        if (rc < rmin) {
            rmin = rc;
            reach = r + k * k * rmin;
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
        double cell = 1.0;
        for (size_t ci = 0; ci < n && cell > 0.0; ++ci) {
            if (ci == di || cand_r[ci] >= k * rd) continue;
            const double R = mbis_dist(xyz[cand[di]].elem, xyz[cand[ci]].elem);
            cell *= mbis_ssf_s((rd - cand_r[ci]) / R);
        }
        if (di == 0) cell_b = cell;
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
    double*   dloc;         // [cap_ao][cap_ao]: the density matrix over a batch's shells
    double*   phi;          // [cap_ao]: the non-negligible AOs at a point
    double*   tmp;          // [cap_ao]
    uint32_t* shells;       // [num_shells]: a batch's shells
    uint32_t* loc;          // [num_shells]: their first row in dloc
    uint32_t* run_beg;      // [num_shells]
    uint32_t* run_col;      // [num_shells]
    uint32_t* run_len;      // [num_shells]
    uint32_t* cand;         // [num_atoms]
    double*   cand_r;       // [num_atoms]
    uint32_t* sub;          // [num_atoms]
    uint32_t* term_k;       // [num_shells]
    double*   term_v;       // [num_shells]
    double*   term_r;       // [num_shells]
    size_t    cap_ao;
} mbis_thread_t;

typedef struct mbis_ctx_t {
    const md_mbis_desc_t* desc;
    const dvec3_t*        xyz;
    size_t                num_atoms;
    size_t                num_ao;

    mbis_shell_t*         shells;
    uint32_t              num_shells;
    double*               alpha;
    double*               coeff;
    mbis_cart_t           cart[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];

    uint32_t              num_rad;
    double*               rad_r;          // [num_atoms][num_rad]
    double*               rad_w;          // [num_atoms][num_rad]
    mbis_sphere_t         sphere[2];      // the full angular grid, and the coarse one for the points near a nucleus
    double*               inner_r;        // [num_atoms]: radial points closer than this use sphere[1]

    mbis_neighbor_t*      nb;             // [num_atoms][num_atoms], sorted by distance from the row's atom
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
    double*               atom_logpre;    // [num_atoms], max ln(pre) of the atom's shells
    double*               atom_smin;      // [num_atoms], min S of the atom's shells
    uint32_t              num_chunks;
    double*               acc;            // [num_chunks][2 num_pro]
} mbis_ctx_t;

// The angular grid of atom b's radial point at r: near the nucleus the density is close to spherical and
// the cell function is 1, so a coarse grid does there.
static inline const mbis_sphere_t* mbis_sphere_at(const mbis_ctx_t* ctx, uint32_t b, double r) {
    return r < ctx->inner_r[b] ? &ctx->sphere[1] : &ctx->sphere[0];
}

static bool mbis_thread_reserve(mbis_thread_t* t, size_t n) {
    if (n <= t->cap_ao) return true;
    size_t cap = MD_MAX(n, t->cap_ao * 2);
    double* dloc = (double*)realloc(t->dloc, sizeof(double) * cap * cap);
    if (dloc) t->dloc = dloc;
    double* phi  = (double*)realloc(t->phi,  sizeof(double) * cap);
    if (phi) t->phi = phi;
    double* tmp = (double*)realloc(t->tmp, sizeof(double) * cap);
    if (tmp) t->tmp = tmp;
    if (!dloc || !phi || !tmp) return false;
    t->cap_ao = cap;
    return true;
}

// The density on atom b's grid and the points worth keeping. A batch is one polar band of one radial sphere:
// the points on a circle, which is small enough that the shells reaching it are few even in a large system.
// The density matrix is gathered once per batch over those shells; per point, only the shells that are not
// negligible there enter phi^T D phi, over the upper triangle.
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

            for (uint32_t band = 0; band < sph->num_bands; ++band) {
                const uint32_t j0 = band * sph->band_size;
                // The band's circle: centre and radius
                const double ct = sph->u[3 * j0 + 2];
                const double cc[3] = {cb[0], cb[1], cb[2] + r * ct};
                const double rcirc = r * sqrt(MD_MAX(0.0, 1.0 - ct * ct));

                // Shells that reach it, and their offsets in the gathered block
                uint32_t ns = 0;
                size_t nl = 0;
                for (uint32_t s = 0; s < ctx->num_shells; ++s) {
                    const mbis_shell_t* sh = &ctx->shells[s];
                    const double rc = sqrt(sh->r2_cut) + rcirc;
                    const double d[3] = {sh->c[0] - cc[0], sh->c[1] - cc[1], sh->c[2] - cc[2]};
                    if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > rc * rc) continue;
                    t->shells[ns] = s;
                    t->loc[ns] = (uint32_t)nl;
                    ++ns;
                    nl += sh->num_ao;
                }
                if (ns == 0) continue;
                if (!mbis_thread_reserve(t, nl)) { ctx->oom = true; return; }

                // The density matrix restricted to those shells
                for (uint32_t a = 0; a < ns; ++a) {
                    const mbis_shell_t* sa = &ctx->shells[t->shells[a]];
                    for (uint32_t i = 0; i < sa->num_ao; ++i) {
                        const double* drow = D + (size_t)(sa->ao + i) * N;
                        double* dst = t->dloc + (size_t)(t->loc[a] + i) * nl;
                        for (uint32_t c = 0; c < ns; ++c) {
                            const mbis_shell_t* sc = &ctx->shells[t->shells[c]];
                            memcpy(dst + t->loc[c], drow + sc->ao, sizeof(double) * sc->num_ao);
                        }
                    }
                }

                for (uint32_t j = j0; j < j0 + sph->band_size; ++j) {
                    const double* u = sph->u + 3 * j;
                    const double p[3] = {cb[0] + r * u[0], cb[1] + r * u[1], cb[2] + r * u[2]};

                    // The AOs of the shells not negligible here, compacted, and where they sit in the block: run i
                    // is phi[run_beg[i] + q] at block row and column run_col[i] + q, q < run_len[i]
                    size_t nnz = 0;
                    uint32_t nrun = 0;
                    for (uint32_t a = 0; a < ns; ++a) {
                        const mbis_shell_t* sh = &ctx->shells[t->shells[a]];
                        const double dd[3] = {p[0] - sh->c[0], p[1] - sh->c[1], p[2] - sh->c[2]};
                        const double r2 = dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2];
                        if (r2 > sh->r2_cut) continue;
                        double R = 0.0;
                        for (uint32_t q = 0; q < sh->num_prim; ++q) R += ctx->coeff[sh->prim_off + q] * exp(-ctx->alpha[sh->prim_off + q] * r2);
                        double px[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], py[MD_GTO_MAX_ANGULAR_MOMENTUM + 1], pz[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];
                        px[0] = py[0] = pz[0] = 1.0;
                        for (uint32_t e = 1; e <= sh->l; ++e) { px[e] = px[e - 1] * dd[0]; py[e] = py[e - 1] * dd[1]; pz[e] = pz[e - 1] * dd[2]; }
                        const mbis_cart_t* ctab = &ctx->cart[sh->l];
                        // Runs of consecutive block columns: a shell continues the last run when it is next in the block
                        if (nrun > 0 && t->run_col[nrun - 1] + t->run_len[nrun - 1] == t->loc[a]) {
                            t->run_len[nrun - 1] += sh->num_ao;
                        } else {
                            t->run_beg[nrun] = (uint32_t)nnz;
                            t->run_col[nrun] = t->loc[a];
                            t->run_len[nrun] = sh->num_ao;
                            ++nrun;
                        }
                        for (uint32_t n = 0; n < sh->num_ao; ++n) {
                            t->phi[nnz++] = ctab->nrm[n] * px[ctab->ijk[n][0]] * py[ctab->ijk[n][1]] * pz[ctab->ijk[n][2]] * R;
                        }
                    }
                    if (nnz == 0) continue;

                    // rho = sum_l phi_l tmp_l, tmp_l = D_ll phi_l + 2 sum_{k < l} D_kl phi_k: the upper triangle, as
                    // axpys along the runs (which vectorize, where a dot product's reduction would not)
                    double* tmp = t->tmp;
                    memset(tmp, 0, sizeof(double) * nnz);
                    for (uint32_t ri = 0; ri < nrun; ++ri) {
                        for (uint32_t o = 0; o < t->run_len[ri]; ++o) {
                            const size_t k = t->run_beg[ri] + o;
                            const double f = t->phi[k];
                            const double f2 = 2.0 * f;
                            const double* drow = t->dloc + (size_t)(t->run_col[ri] + o) * nl;
                            tmp[k] += f * drow[t->run_col[ri] + o];
                            {   // the rest of its own run
                                const uint32_t len = t->run_len[ri] - o - 1;
                                double* dst = tmp + k + 1;
                                const double* src = drow + t->run_col[ri] + o + 1;
                                for (uint32_t q = 0; q < len; ++q) dst[q] += f2 * src[q];
                            }
                            for (uint32_t rj = ri + 1; rj < nrun; ++rj) {
                                double* dst = tmp + t->run_beg[rj];
                                const double* src = drow + t->run_col[rj];
                                const uint32_t len = t->run_len[rj];
                                for (uint32_t q = 0; q < len; ++q) dst[q] += f2 * src[q];
                            }
                        }
                    }
                    double rho = 0.0;
                    for (size_t k = 0; k < nnz; ++k) rho += t->phi[k] * tmp[k];
                    if (!(rho > MBIS_RHO_FLOOR)) continue;

                    const double wb = mbis_becke_weight(p, r, ctx->nb + (size_t)b * ctx->num_atoms, ctx->num_atoms, ctx->xyz, t->cand, t->cand_r);
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
                    nsub = 0;
                    for (size_t i = 0; i < na; ++i) {
                        const uint32_t a = nb[i].atom;
                        const double lo = fabs(nb[i].d - r);   // the closest a gets to the sphere
                        if (ctx->atom_logpre[a] - ctx->atom_smin[a] * lo < lf) continue;
                        t->sub[nsub++] = a;
                    }
                }

                const double* u = sph->u + 3 * (size_t)pts->ang[ip];
                const double p[3] = {cb[0] + r * u[0], cb[1] + r * u[1], cb[2] + r * u[2]};
                size_t nt = 0;
                double rho0 = 0.0;
                for (size_t i = 0; i < nsub; ++i) {
                    const uint32_t a = t->sub[i];
                    const double ra = (a == b) ? r : mbis_dist(p, ctx->xyz[a].elem);
                    for (uint32_t k = ctx->shell_off[a]; k < ctx->shell_off[a + 1]; ++k) {
                        const double v = ctx->pre[k] * exp(-ctx->S[k] * ra);
                        if (v < floor_b) continue;
                        t->term_k[nt] = k;
                        t->term_v[nt] = v;
                        t->term_r[nt] = ra;
                        ++nt;
                        rho0 += v;
                    }
                }
                if (!(rho0 > 0.0)) continue;
                const double scale = pts->mass[ip] / rho0;
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
    free(ctx->shells); free(ctx->alpha); free(ctx->coeff);
    free(ctx->rad_r); free(ctx->rad_w); free(ctx->inner_r);
    for (int i = 0; i < 2; ++i) { free(ctx->sphere[i].u); free(ctx->sphere[i].w); }
    free(ctx->nb);
    if (ctx->pts) {
        for (size_t i = 0; i < ctx->num_atoms; ++i) { free(ctx->pts[i].rad); free(ctx->pts[i].ang); free(ctx->pts[i].mass); }
        free(ctx->pts);
    }
    if (ctx->thr) {
        for (uint32_t i = 0; i < num_threads; ++i) {
            mbis_thread_t* t = &ctx->thr[i];
            free(t->phi); free(t->dloc); free(t->tmp); free(t->shells); free(t->loc); free(t->run_beg); free(t->run_col); free(t->run_len); free(t->cand); free(t->cand_r);
            free(t->sub); free(t->term_k); free(t->term_v); free(t->term_r);
        }
        free(ctx->thr);
    }
    free(ctx->shell_off); free(ctx->N); free(ctx->S); free(ctx->pre);
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
    ctx.coeff  = (double*)calloc(basis->num_primitives ? basis->num_primitives : 1, sizeof(double));
    if (!ctx.shells || !ctx.alpha || !ctx.coeff) goto oom;
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
        }
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
    if (!ctx.nb) goto oom;
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
    ctx.atom_logpre = (double*)malloc(sizeof(double) * na);
    ctx.atom_smin = (double*)malloc(sizeof(double) * na);
    ctx.num_chunks = (uint32_t)MD_MIN(na, (size_t)MBIS_MAX_CHUNKS);
    ctx.acc = (double*)malloc(sizeof(double) * 2 * ctx.num_pro * ctx.num_chunks);
    if (!ctx.N || !ctx.S || !ctx.pre || !ctx.atom_logpre || !ctx.atom_smin || !ctx.acc) goto oom;
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
        t->cand   = (uint32_t*)malloc(sizeof(uint32_t) * na);
        t->cand_r = (double*)malloc(sizeof(double) * na);
        t->sub    = (uint32_t*)malloc(sizeof(uint32_t) * na);
        t->term_k = (uint32_t*)malloc(sizeof(uint32_t) * ctx.num_pro);
        t->term_v = (double*)malloc(sizeof(double) * ctx.num_pro);
        t->term_r = (double*)malloc(sizeof(double) * ctx.num_pro);
        if (!t->shells || !t->loc || !t->run_beg || !t->run_col || !t->run_len || !t->cand || !t->cand_r || !t->sub || !t->term_k || !t->term_v || !t->term_r) goto oom;
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
    while (iter < max_iter) {
        if (desc->cancel && *desc->cancel) { free(aa); goto fail; }
        for (size_t k = 0; k < P; ++k) { ctx.N[k] = x[k]; ctx.S[k] = exp(x[P + k]); }
        for (size_t a = 0; a < na; ++a) {
            double lp = -HUGE_VAL, smin = HUGE_VAL;
            for (uint32_t k = ctx.shell_off[a]; k < ctx.shell_off[a + 1]; ++k) {
                ctx.pre[k] = ctx.N[k] * ctx.S[k] * ctx.S[k] * ctx.S[k] / (8.0 * 3.14159265358979323846);
                if (ctx.pre[k] > 0.0) lp = MD_MAX(lp, log(ctx.pre[k]));
                smin = MD_MIN(smin, ctx.S[k]);
            }
            ctx.atom_logpre[a] = lp;
            ctx.atom_smin[a] = smin;
        }
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
        if (change < threshold) {
            memcpy(x, g, sizeof(double) * X);
            break;
        }

        // The history of differences, oldest overwritten first
        if (iter > 1) {
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
