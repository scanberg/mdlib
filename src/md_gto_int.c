#include <md_gto_int.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>

#include <math.h>
#include <string.h>

#define GTO_INT_PI 3.14159265358979323846

// Hermite order needed to evaluate the field as well: one above the highest gaussian order.
#define GTO_INT_MAX_R_ORDER (MD_GTO_INT_MAX_ORDER + 1)
#define GTO_INT_MAX_R       (((GTO_INT_MAX_R_ORDER + 1) * (GTO_INT_MAX_R_ORDER + 2) * (GTO_INT_MAX_R_ORDER + 3)) / 6)
#define GTO_INT_MAX_H       (((MD_GTO_INT_MAX_ORDER + 1) * (MD_GTO_INT_MAX_ORDER + 2) * (MD_GTO_INT_MAX_ORDER + 3)) / 6)
#define GTO_INT_MAX_CART    (((MD_GTO_MAX_ANGULAR_MOMENTUM + 1) * (MD_GTO_MAX_ANGULAR_MOMENTUM + 2)) / 2)

// ---------------------------------------------------------------------------
// Hermite index
// ---------------------------------------------------------------------------
// Position of (t,u,v) in the coefficient order documented in md_gto_int.h: degree n = t+u+v
// ascending, then t descending, then u descending. Closed form, so no table has to be kept in
// step with the generated GPU code (tools/gen_gto_int_hermite.py uses the same expression).
static inline uint32_t herm_idx(int t, int u, int v) {
    const int n = t + u + v;
    return (uint32_t)(n * (n + 1) * (n + 2) / 6 + (n - t) * (n - t + 1) / 2 + (n - t - u));
}

// max over Y of |R~_tuv(Y)|, the Hermite Coulomb integral at unit exponent, in herm_idx order, with
// a 1% margin. R_tuv(p, X) = p^((t+u+v)/2) R~_tuv(sqrt(p) X), so this bounds every gaussian's
// potential anywhere in space from its coefficients alone. Sampled on [0,7]^3 at 0.025 (the
// maxima all lie within |Y| < 3; |R~| is even or odd per axis).
static const double gto_int_rmax[GTO_INT_MAX_H] = {
    1.010000e+00, 3.830785e-01, 3.830785e-01, 3.830785e-01, 6.733333e-01, 2.194139e-01,
    2.194139e-01, 6.733333e-01, 2.194139e-01, 6.733333e-01, 9.417070e-01, 4.149862e-01,
    4.149862e-01, 4.149862e-01, 1.427321e-01, 4.149862e-01, 9.417070e-01, 4.149862e-01,
    4.149862e-01, 9.417070e-01, 2.424000e+00, 6.385318e-01, 6.385318e-01, 8.080000e-01,
    2.004960e-01, 8.080000e-01, 6.385318e-01, 2.004960e-01, 2.004960e-01, 6.385318e-01,
    2.424000e+00, 6.385318e-01, 8.080000e-01, 6.385318e-01, 2.424000e+00, 4.905143e+00,
    1.832291e+00, 1.832291e+00, 1.243066e+00, 6.024508e-01, 1.243066e+00, 1.243066e+00,
    3.807891e-01, 3.807891e-01, 1.243066e+00, 1.832291e+00, 6.024508e-01, 3.807891e-01,
    6.024508e-01, 1.832291e+00, 4.905143e+00, 1.832291e+00, 1.243066e+00, 1.243066e+00,
    1.832291e+00, 4.905143e+00, 1.540853e+01, 4.064329e+00, 4.064329e+00, 3.672727e+00,
    1.276992e+00, 3.672727e+00, 2.204017e+00, 9.411368e-01, 9.411368e-01, 2.204017e+00,
    3.672727e+00, 9.411368e-01, 6.019017e-01, 9.411368e-01, 3.672727e+00, 4.064329e+00,
    1.276992e+00, 9.411368e-01, 9.411368e-01, 1.276992e+00, 4.064329e+00, 1.540853e+01,
    4.064329e+00, 3.672727e+00, 2.204017e+00, 3.672727e+00, 4.064329e+00, 1.540853e+01,
    4.476440e+01, 1.459289e+01, 1.459289e+01, 1.017422e+01, 4.426212e+00, 1.017422e+01,
    8.170864e+00, 2.776339e+00, 2.776339e+00, 8.170864e+00, 8.170864e+00, 2.330181e+00,
    1.689962e+00, 2.330181e+00, 8.170864e+00, 1.017422e+01, 2.776339e+00, 1.689962e+00,
    1.689962e+00, 2.776339e+00, 1.017422e+01, 1.459289e+01, 4.426212e+00, 2.776339e+00,
    2.330181e+00, 2.776339e+00, 4.426212e+00, 1.459289e+01, 4.476440e+01, 1.459289e+01,
    1.017422e+01, 8.170864e+00, 8.170864e+00, 1.017422e+01, 1.459289e+01, 4.476440e+01,
    1.885333e+02, 4.937815e+01, 4.937815e+01, 4.380127e+01, 1.429393e+01, 4.380127e+01,
    2.693333e+01, 1.091117e+01, 1.091117e+01, 2.693333e+01, 2.693333e+01, 7.618018e+00,
    6.024051e+00, 7.618018e+00, 2.693333e+01, 4.380127e+01, 1.091117e+01, 6.024051e+00,
    6.024051e+00, 1.091117e+01, 4.380127e+01, 4.937815e+01, 1.429393e+01, 1.091117e+01,
    7.618018e+00, 1.091117e+01, 1.429393e+01, 4.937815e+01, 1.885333e+02, 4.937815e+01,
    4.380127e+01, 2.693333e+01, 2.693333e+01, 4.380127e+01, 4.937815e+01, 1.885333e+02,
    3.786977e+01, 1.401564e+01, 1.616000e+01, 1.401564e+01, 3.786977e+01, 1.885333e+02,
    1.401564e+01, 2.693333e+01, 3.786977e+01,
};

// ---------------------------------------------------------------------------
// Boys function
// ---------------------------------------------------------------------------
// F_n(T) = int_0^1 t^2n exp(-T t^2) dt for n = 0..nmax, to double precision.
// Below T = 30: the series for F_nmax, which has only positive terms, then downward recursion,
// which is stable. Above: F_0 from erf and upward recursion, stable there because
// (2n+1) / 2T < 1 for every n this library needs.
static void gto_int_boys(double* F, int nmax, double T) {
    if (T < 1.0e-15) {
        for (int n = 0; n <= nmax; ++n) F[n] = 1.0 / (2 * n + 1);
        return;
    }
    const double e = exp(-T);
    if (T >= 30.0) {
        const double inv2T = 0.5 / T;
        F[0] = 0.5 * sqrt(GTO_INT_PI / T) * erf(sqrt(T));
        for (int n = 0; n < nmax; ++n) F[n + 1] = ((2 * n + 1) * F[n] - e) * inv2T;
        return;
    }
    double term = 1.0 / (2 * nmax + 1);
    double sum  = term;
    for (int k = 1; k < 200; ++k) {
        term *= 2.0 * T / (2 * nmax + 2 * k + 1);
        sum  += term;
        if (term < 1.0e-17 * sum) break;
    }
    F[nmax] = e * sum;
    for (int n = nmax - 1; n >= 0; --n) F[n] = (2.0 * T * F[n + 1] + e) / (2 * n + 1);
}

// ---------------------------------------------------------------------------
// Hermite Coulomb integrals
// ---------------------------------------------------------------------------
// R[herm_idx(t,u,v)] = R^0_tuv(p, X) for all t+u+v <= L, X = P - C, by the McMurchie-Davidson
// recursion
//     R^n_000       = (-2p)^n F_n(p |X|^2)
//     R^n_{t+1,u,v} = t R^{n+1}_{t-1,u,v} + X_x R^{n+1}_{t,u,v}     (likewise for u, v)
// run from the top auxiliary index down, so every level only reads the one above it.
static void gto_int_hermite_R(double* R, int L, double p, double x, double y, double z) {
    double F[GTO_INT_MAX_R_ORDER + 1] = {0};
    gto_int_boys(F, L, p * (x * x + y * y + z * z));

    double buf[2][GTO_INT_MAX_R];
    double* cur = buf[0];
    double* nxt = buf[1];

    const double m2p = -2.0 * p;
    double s = 1.0;
    for (int n = 0; n < L; ++n) s *= m2p;
    nxt[0] = s * F[L];

    for (int lvl = L - 1; lvl >= 0; --lvl) {
        s = 1.0;
        for (int n = 0; n < lvl; ++n) s *= m2p;
        cur[0] = s * F[lvl];
        const int top = L - lvl;
        uint32_t i = 1;
        for (int deg = 1; deg <= top; ++deg) {
            for (int t = deg; t >= 0; --t) {
                for (int u = deg - t; u >= 0; --u, ++i) {
                    const int v = deg - t - u;
                    double r;
                    if (t > 0) {
                        r = x * nxt[herm_idx(t - 1, u, v)];
                        if (t > 1) r += (t - 1) * nxt[herm_idx(t - 2, u, v)];
                    } else if (u > 0) {
                        r = y * nxt[herm_idx(t, u - 1, v)];
                        if (u > 1) r += (u - 1) * nxt[herm_idx(t, u - 2, v)];
                    } else {
                        r = z * nxt[herm_idx(t, u, v - 1)];
                        if (v > 1) r += (v - 1) * nxt[herm_idx(t, u, v - 2)];
                    }
                    cur[i] = r;
                }
            }
        }
        double* tmp = cur; cur = nxt; nxt = tmp;
    }
    MEMCPY(R, nxt, sizeof(double) * md_gto_int_num_hermite((uint32_t)L));
}

// ---------------------------------------------------------------------------
// Hermite expansion of a primitive product
// ---------------------------------------------------------------------------
// One axis: x_A^i x_B^j exp(-a x_A^2) exp(-b x_B^2) = sum_t E[i][j][t] Lambda_t(x; p, P), by
//     E^{00}_0      = exp(-mu X_AB^2)
//     E^{i+1,j}_t   = E^{ij}_{t-1} / 2p + X_PA E^{ij}_t + (t+1) E^{ij}_{t+1}
//     E^{i,j+1}_t   = E^{ij}_{t-1} / 2p + X_PB E^{ij}_t + (t+1) E^{ij}_{t+1}
typedef double gto_int_E_t[MD_GTO_MAX_ANGULAR_MOMENTUM + 1][MD_GTO_MAX_ANGULAR_MOMENTUM + 1][MD_GTO_INT_MAX_ORDER + 1];

static void gto_int_hermite_E(gto_int_E_t E, int la, int lb, double p, double PA, double PB, double K) {
    MEMSET(E, 0, sizeof(gto_int_E_t));
    const double o2p = 0.5 / p;
    E[0][0][0] = K;
    for (int i = 0; i <= la; ++i) {
        for (int j = 0; j <= lb; ++j) {
            if (i == 0 && j == 0) continue;
            // Raise i from (i-1, j) while possible, otherwise j from (0, j-1).
            const int pi = (i > 0) ? i - 1 : i;
            const int pj = (i > 0) ? j     : j - 1;
            const double X = (i > 0) ? PA : PB;
            const int top = pi + pj;   // highest t of the parent
            for (int t = 0; t <= i + j; ++t) {
                double v = 0.0;
                if (t <= top)     v += X * E[pi][pj][t];
                if (t > 0)        v += o2p * E[pi][pj][t - 1];
                if (t + 1 <= top) v += (t + 1) * E[pi][pj][t + 1];
                E[i][j][t] = v;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Building the distribution
// ---------------------------------------------------------------------------

typedef struct MD_HASHMAP_T(uint32_t) gto_int_map_t;

// md_hashmap indexes by the low bits of the key, and neither key here is spread out there: a site
// key's low bits are an exponent's mantissa, often mostly zero, and a pair key's are a site index.
// So keys are mixed first. The mix is a bijection, so distinct keys stay distinct; the two values
// the map reserves are moved aside, which could collide with another key with probability 2^-63
// (the same trade md_contact.c makes).
static inline uint64_t gto_int_mix_key(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebull;
    x ^= x >> 31;
    if (x >= MD_HASH_TOMBSTONE) x ^= (1ull << 63);
    return x;
}

static inline uint32_t gto_int_float_bits(float f) {
    uint32_t u;
    MEMCPY(&u, &f, sizeof(u));
    return u;
}

static inline const float* gto_int_xyz(const float* base, size_t stride, size_t i) {
    return (const float*)((const char*)base + i * stride);
}

typedef struct {
    double   P[3];
    double   p;
    uint32_t L;
    uint32_t coeff_offset;   // into the unsorted scratch coefficients
} gto_int_scratch_t;

static bool gto_int_build_gaussians(md_gto_int_charges_t* out, const md_gto_int_charges_desc_t* desc, md_allocator_i* alloc) {
    const md_gto_basis_t* basis = desc->basis;
    const size_t stride = desc->atom_xyz_stride ? desc->atom_xyz_stride : sizeof(float) * 3;
    const double scale  = desc->density_scale != 0.0 ? desc->density_scale : -1.0;
    const double threshold = MAX(desc->threshold, 0.0);
    const double* D = desc->density_matrix;
    const size_t N = md_gto_basis_num_ao(basis);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    // Shell offsets into the AO axis, and every shell's primitives mapped to a SITE: a unique
    // (atom, exponent) pair. Two primitives at the same site are the same function up to a
    // coefficient, so a pair of sites identifies a gaussian exactly. Keys are exact: the atom
    // index and the bits of the float exponent the basis stores.
    uint32_t* ao_off   = md_temp_alloc_array(temp, uint32_t, basis->num_shells);
    uint32_t* site_off = md_temp_alloc_array(temp, uint32_t, basis->num_shells);
    if (!ao_off || !site_off) goto done;

    uint32_t num_shell_prims = 0;
    {
        uint32_t ao = 0;
        for (uint32_t s = 0; s < basis->num_shells; ++s) {
            const md_gto_shell_t* sh = &basis->shells[s];
            if (sh->l > MD_GTO_MAX_ANGULAR_MOMENTUM) {
                MD_LOG_ERROR("md_gto_int: shell %u has angular momentum %u, above the supported %d", s, sh->l, MD_GTO_MAX_ANGULAR_MOMENTUM);
                goto done;
            }
            if ((size_t)sh->primitive_offset + sh->num_primitives > basis->num_primitives) {
                MD_LOG_ERROR("md_gto_int: shell %u spans primitives beyond the basis", s);
                goto done;
            }
            ao_off[s]   = ao;
            site_off[s] = num_shell_prims;
            ao += md_gto_num_cart_ao(sh->l);
            num_shell_prims += sh->num_primitives;
        }
    }

    uint32_t* site = md_temp_alloc_array(temp, uint32_t, num_shell_prims);
    if (!site) goto done;
    {
        gto_int_map_t site_map = { .allocator = temp_alloc };
        uint32_t num_sites = 0;
        for (uint32_t s = 0; s < basis->num_shells; ++s) {
            const md_gto_shell_t* sh = &basis->shells[s];
            for (uint32_t k = 0; k < sh->num_primitives; ++k) {
                const uint64_t key = gto_int_mix_key(((uint64_t)sh->atom_idx << 32) | gto_int_float_bits(basis->alpha[sh->primitive_offset + k]));
                uint32_t* found = md_hashmap_get(&site_map, key);
                if (found) {
                    site[site_off[s] + k] = *found;
                } else {
                    site[site_off[s] + k] = num_sites;
                    md_hashmap_add(&site_map, key, num_sites);
                    num_sites++;
                }
            }
        }
    }

    // Largest |D_mu,nu + D_nu,mu| of every shell pair's block, the density part of the prescreen.
    // Pass 1 finds the gaussians and their order, pass 2 accumulates their coefficients. Both visit
    // the same primitive pairs in the same order, so the prescreen decides identically in both.
    gto_int_map_t pair_map = { .allocator = temp_alloc };
    md_array(gto_int_scratch_t) gauss = NULL;
    double* scratch_coeff = NULL;
    size_t num_scratch_coeff = 0;

    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1) {
            const size_t num = md_array_size(gauss);
            for (size_t g = 0; g < num; ++g) {
                gauss[g].coeff_offset = (uint32_t)num_scratch_coeff;
                num_scratch_coeff += md_gto_int_num_hermite(gauss[g].L);
            }
            scratch_coeff = md_temp_alloc_zero_array(temp, double, MAX(num_scratch_coeff, 1));
            if (!scratch_coeff) goto done;
        }

        for (uint32_t sa = 0; sa < basis->num_shells; ++sa) {
            const md_gto_shell_t* A = &basis->shells[sa];
            const float* Af = gto_int_xyz(desc->atom_xyz, stride, A->atom_idx);
            const double Ar[3] = { Af[0], Af[1], Af[2] };
            const int la = (int)A->l;
            const uint32_t na = md_gto_num_cart_ao(A->l);

            for (uint32_t sb = sa; sb < basis->num_shells; ++sb) {
                const md_gto_shell_t* B = &basis->shells[sb];
                const float* Bf = gto_int_xyz(desc->atom_xyz, stride, B->atom_idx);
                const double Br[3] = { Bf[0], Bf[1], Bf[2] };
                const int lb = (int)B->l;
                const uint32_t nb = md_gto_num_cart_ao(B->l);
                const bool same = (sa == sb);

                // The weight of every AO pair of the block. A diagonal block is summed over its
                // full square, an off-diagonal one over (mu,nu) with both orders folded in.
                double w[GTO_INT_MAX_CART][GTO_INT_MAX_CART];
                double wmax = 0.0;
                for (uint32_t ca = 0; ca < na; ++ca) {
                    const size_t mu = ao_off[sa] + ca;
                    for (uint32_t cb = 0; cb < nb; ++cb) {
                        const size_t nu = ao_off[sb] + cb;
                        const double d = same ? D[mu * N + nu] : (D[mu * N + nu] + D[nu * N + mu]);
                        w[ca][cb] = d;
                        wmax = MAX(wmax, fabs(d));
                    }
                }
                if (wmax == 0.0) continue;

                const double AB[3] = { Ar[0] - Br[0], Ar[1] - Br[1], Ar[2] - Br[2] };
                const double R2 = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];

                for (uint32_t pa = 0; pa < A->num_primitives; ++pa) {
                    const double a  = basis->alpha[A->primitive_offset + pa];
                    const double ca = basis->coeff[A->primitive_offset + pa];
                    for (uint32_t pb = 0; pb < B->num_primitives; ++pb) {
                        const double b  = basis->alpha[B->primitive_offset + pb];
                        const double cb = basis->coeff[B->primitive_offset + pb];
                        const double p  = a + b;
                        const double mu = a * b / p;
                        const double K  = exp(-mu * R2);
                        const double P[3] = { (a * Ar[0] + b * Br[0]) / p, (a * Ar[1] + b * Br[1]) / p, (a * Ar[2] + b * Br[2]) / p };

                        // Prescreen: a generous bound on this pair's potential anywhere. It only
                        // has to be safe; the real screening is on the merged gaussian below.
                        const double cc = fabs(scale * ca * cb) * wmax * K;
                        if (cc == 0.0) continue;
                        if (threshold > 0.0) {
                            const double dPA = sqrt((P[0] - Ar[0]) * (P[0] - Ar[0]) + (P[1] - Ar[1]) * (P[1] - Ar[1]) + (P[2] - Ar[2]) * (P[2] - Ar[2]));
                            const double dPB = sqrt((P[0] - Br[0]) * (P[0] - Br[0]) + (P[1] - Br[1]) * (P[1] - Br[1]) + (P[2] - Br[2]) * (P[2] - Br[2]));
                            const double sp  = 1.0 + sqrt(p);
                            double poly = (double)(na * nb);
                            for (int k = 0; k < la; ++k) poly *= (1.0 + dPA) * sp;
                            for (int k = 0; k < lb; ++k) poly *= (1.0 + dPB) * sp;
                            if (cc * poly * (2.0 * GTO_INT_PI / p) * gto_int_rmax[0] * 100.0 < 1.0e-6 * threshold) continue;
                        }

                        uint32_t s0 = site[site_off[sa] + pa];
                        uint32_t s1 = site[site_off[sb] + pb];
                        if (s0 > s1) { uint32_t t = s0; s0 = s1; s1 = t; }
                        const uint64_t key = gto_int_mix_key(((uint64_t)s0 << 32) | s1);

                        if (pass == 0) {
                            uint32_t* found = md_hashmap_get(&pair_map, key);
                            if (found) {
                                gauss[*found].L = MAX(gauss[*found].L, (uint32_t)(la + lb));
                            } else {
                                gto_int_scratch_t g = { .P = { P[0], P[1], P[2] }, .p = p, .L = (uint32_t)(la + lb) };
                                md_hashmap_add(&pair_map, key, (uint32_t)md_array_size(gauss));
                                md_array_push(gauss, g, temp_alloc);
                            }
                            continue;
                        }

                        uint32_t* found = md_hashmap_get(&pair_map, key);
                        ASSERT(found);
                        double* h = scratch_coeff + gauss[*found].coeff_offset;

                        gto_int_E_t Ex, Ey, Ez;
                        gto_int_hermite_E(Ex, la, lb, p, P[0] - Ar[0], P[0] - Br[0], exp(-mu * AB[0] * AB[0]));
                        gto_int_hermite_E(Ey, la, lb, p, P[1] - Ar[1], P[1] - Br[1], exp(-mu * AB[1] * AB[1]));
                        gto_int_hermite_E(Ez, la, lb, p, P[2] - Ar[2], P[2] - Br[2], exp(-mu * AB[2] * AB[2]));

                        const double pref = scale * ca * cb * 2.0 * GTO_INT_PI / p;
                        for (uint32_t ia_c = 0; ia_c < na; ++ia_c) {
                            int ia, ja, ka;
                            md_gto_cart_ijk(&ia, &ja, &ka, A->l, ia_c);
                            const double fa = md_gto_cart_norm_factor(ia, ja, ka);
                            for (uint32_t ib_c = 0; ib_c < nb; ++ib_c) {
                                const double wab = w[ia_c][ib_c];
                                if (wab == 0.0) continue;
                                int ib, jb, kb;
                                md_gto_cart_ijk(&ib, &jb, &kb, B->l, ib_c);
                                const double fb = md_gto_cart_norm_factor(ib, jb, kb);
                                const double wt = pref * wab * fa * fb;
                                for (int t = 0; t <= ia + ib; ++t) {
                                    const double wx = wt * Ex[ia][ib][t];
                                    for (int u = 0; u <= ja + jb; ++u) {
                                        const double wxy = wx * Ey[ja][jb][u];
                                        for (int v = 0; v <= ka + kb; ++v) {
                                            h[herm_idx(t, u, v)] += wxy * Ez[ka][kb][v];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Screen on the merged gaussians and sort by order (a counting sort, stable).
    {
        const size_t num = md_array_size(gauss);
        double*   bound = md_temp_alloc_array(temp, double, MAX(num, 1));
        uint32_t  count[MD_GTO_INT_MAX_ORDER + 1] = {0};
        size_t    kept_coeffs = 0;
        if (!bound) goto done;

        for (size_t g = 0; g < num; ++g) {
            const double* h = scratch_coeff + gauss[g].coeff_offset;
            const double sp = sqrt(gauss[g].p);
            double B = 0.0;
            uint32_t i = 0;
            double pw = 1.0;   // p^(deg/2)
            for (uint32_t deg = 0; deg <= gauss[g].L; ++deg, pw *= sp) {
                const uint32_t n_deg = (deg + 1) * (deg + 2) / 2;
                for (uint32_t k = 0; k < n_deg; ++k, ++i) B += fabs(h[i]) * pw * gto_int_rmax[i];
            }
            bound[g] = B;
            if (B > 0.0 && B >= threshold) {
                count[gauss[g].L]++;
                kept_coeffs += md_gto_int_num_hermite(gauss[g].L);
            }
        }

        uint32_t num_kept = 0;
        out->order_offset[0] = 0;
        for (uint32_t L = 0; L <= MD_GTO_INT_MAX_ORDER; ++L) {
            out->order_offset[L + 1] = out->order_offset[L] + count[L];
            if (count[L]) out->max_order = L;
        }
        num_kept = out->order_offset[MD_GTO_INT_MAX_ORDER + 1];

        out->num_gaussians = num_kept;
        out->num_coeffs    = kept_coeffs;
        if (num_kept) {
            out->center       = md_alloc(alloc, sizeof(double) * 3 * num_kept);
            out->exponent     = md_alloc(alloc, sizeof(double) * num_kept);
            out->bound        = md_alloc(alloc, sizeof(double) * num_kept);
            out->coeff_offset = md_alloc(alloc, sizeof(uint32_t) * num_kept);
            out->coeff        = md_alloc(alloc, sizeof(double) * kept_coeffs);
            if (!out->center || !out->exponent || !out->bound || !out->coeff_offset || !out->coeff) goto done;

            // Coefficients are laid out in the sorted order too, so that those of order L are one
            // contiguous run and a gaussian's offset is implied by its index (the GPU relies on it).
            uint32_t next[MD_GTO_INT_MAX_ORDER + 1];
            size_t   coeff_next[MD_GTO_INT_MAX_ORDER + 1];
            size_t   run = 0;
            for (uint32_t L = 0; L <= MD_GTO_INT_MAX_ORDER; ++L) {
                next[L] = out->order_offset[L];
                coeff_next[L] = run;
                run += (size_t)count[L] * md_gto_int_num_hermite(L);
            }
            for (size_t g = 0; g < num; ++g) {
                if (!(bound[g] > 0.0 && bound[g] >= threshold)) continue;
                const uint32_t L = gauss[g].L;
                const uint32_t dst = next[L]++;
                const uint32_t nh = md_gto_int_num_hermite(L);
                out->center[dst * 3 + 0] = gauss[g].P[0];
                out->center[dst * 3 + 1] = gauss[g].P[1];
                out->center[dst * 3 + 2] = gauss[g].P[2];
                out->exponent[dst]       = gauss[g].p;
                out->bound[dst]          = bound[g];
                out->coeff_offset[dst]   = (uint32_t)coeff_next[L];
                MEMCPY(out->coeff + coeff_next[L], scratch_coeff + gauss[g].coeff_offset, sizeof(double) * nh);
                coeff_next[L] += nh;
            }
        }
    }
    result = true;

done:
    md_temp_end(temp);
    return result;
}

bool md_gto_int_charges_init(md_gto_int_charges_t* out, const md_gto_int_charges_desc_t* desc, md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(alloc);
    if (!desc) {
        return false;
    }
    MEMSET(out, 0, sizeof(*out));

    const bool has_density = desc->basis && desc->density_matrix && desc->basis->num_shells > 0;
    if (has_density) {
        if (!desc->atom_xyz) {
            MD_LOG_ERROR("md_gto_int_charges_init: a density was given without atom positions");
            return false;
        }
        if (!gto_int_build_gaussians(out, desc, alloc)) {
            md_gto_int_charges_free(out, alloc);
            return false;
        }
    }

    if (desc->num_points > 0) {
        if (!desc->point_xyz || !desc->point_charge) {
            MD_LOG_ERROR("md_gto_int_charges_init: num_points given without positions or charges");
            md_gto_int_charges_free(out, alloc);
            return false;
        }
        const size_t stride = desc->point_xyz_stride ? desc->point_xyz_stride : sizeof(float) * 3;
        out->point_xyz    = md_alloc(alloc, sizeof(double) * 3 * desc->num_points);
        out->point_charge = md_alloc(alloc, sizeof(double) * desc->num_points);
        if (!out->point_xyz || !out->point_charge) {
            md_gto_int_charges_free(out, alloc);
            return false;
        }
        for (size_t i = 0; i < desc->num_points; ++i) {
            const float* x = gto_int_xyz(desc->point_xyz, stride, i);
            out->point_xyz[i * 3 + 0] = x[0];
            out->point_xyz[i * 3 + 1] = x[1];
            out->point_xyz[i * 3 + 2] = x[2];
            out->point_charge[i] = desc->point_charge[i];
        }
        out->num_points = (uint32_t)desc->num_points;
    }
    return true;
}

void md_gto_int_charges_free(md_gto_int_charges_t* q, md_allocator_i* alloc) {
    if (!q) return;
    ASSERT(alloc);
    if (q->center)       md_free(alloc, q->center,       sizeof(double) * 3 * q->num_gaussians);
    if (q->exponent)     md_free(alloc, q->exponent,     sizeof(double) * q->num_gaussians);
    if (q->bound)        md_free(alloc, q->bound,        sizeof(double) * q->num_gaussians);
    if (q->coeff_offset) md_free(alloc, q->coeff_offset, sizeof(uint32_t) * q->num_gaussians);
    if (q->coeff)        md_free(alloc, q->coeff,        sizeof(double) * q->num_coeffs);
    if (q->point_xyz)    md_free(alloc, q->point_xyz,    sizeof(double) * 3 * q->num_points);
    if (q->point_charge) md_free(alloc, q->point_charge, sizeof(double) * q->num_points);
    MEMSET(q, 0, sizeof(*q));
}

// Work of one gaussian at one voxel, in recursion steps: the recursion itself (C(L+4,4) updates),
// the dot product, and a flat allowance for the setup and the Boys function.
static uint64_t gto_int_cost(uint32_t L) {
    return (uint64_t)((L + 1) * (L + 2) * (L + 3) * (L + 4) / 24) + md_gto_int_num_hermite(L) + 16;
}

uint64_t md_gto_int_charges_work_per_point(const md_gto_int_charges_t* q) {
    if (!q) return 0;
    uint64_t w = q->num_points;
    for (uint32_t L = 0; L <= MD_GTO_INT_MAX_ORDER; ++L) {
        w += (uint64_t)(q->order_offset[L + 1] - q->order_offset[L]) * gto_int_cost(L);
    }
    return w;
}

// Order of gaussian g, walking forward from the previous answer: gaussians are visited in order.
static inline uint32_t gto_int_order_of(const md_gto_int_charges_t* q, uint32_t g, uint32_t* L_hint) {
    uint32_t L = *L_hint;
    while (g >= q->order_offset[L + 1]) ++L;
    *L_hint = L;
    return L;
}

// ---------------------------------------------------------------------------
// Moments
// ---------------------------------------------------------------------------
// d^t/dP^t of int x^a exp(-p (x-P)^2) dx, divided by sqrt(pi/p), for a <= 2. Only Hermite terms of
// degree <= a contribute to a moment of order a, which is why moments are cheap.
static inline double gto_int_moment_1d(int a, int t, double P, double p) {
    switch (a) {
    case 0: return t == 0 ? 1.0 : 0.0;
    case 1: return t == 0 ? P : (t == 1 ? 1.0 : 0.0);
    case 2: return t == 0 ? P * P + 0.5 / p : (t == 1 ? 2.0 * P : (t == 2 ? 2.0 : 0.0));
    default: return 0.0;
    }
}

md_gto_int_moments_t md_gto_int_charges_moments(const md_gto_int_charges_t* q, const double origin[3]) {
    md_gto_int_moments_t m = {0};
    if (!q) return m;
    const double O[3] = { origin ? origin[0] : 0.0, origin ? origin[1] : 0.0, origin ? origin[2] : 0.0 };

    // (ax, ay, az) of the ten moments, in the order charge, dipole, second
    static const int mom[10][3] = {
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
        {2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2},
    };
    double acc[10] = {0};

    uint32_t L_hint = 0;
    for (uint32_t g = 0; g < q->num_gaussians; ++g) {
        const double p = q->exponent[g];
        const double P[3] = { q->center[g * 3 + 0] - O[0], q->center[g * 3 + 1] - O[1], q->center[g * 3 + 2] - O[2] };
        // h' carries 2 pi / p; the density coefficient times (pi/p)^1.5 is h' * p/(2 pi) * (pi/p)^1.5
        const double norm = p / (2.0 * GTO_INT_PI) * pow(GTO_INT_PI / p, 1.5);
        const double* h = q->coeff + q->coeff_offset[g];
        const uint32_t L = gto_int_order_of(q, g, &L_hint);
        const uint32_t top = MIN(L, 2u);
        uint32_t i = 0;
        for (uint32_t deg = 0; deg <= top; ++deg) {
            for (int t = (int)deg; t >= 0; --t) {
                for (int u = (int)deg - t; u >= 0; --u, ++i) {
                    const int v = (int)deg - t - u;
                    const double c = h[i] * norm;
                    for (int k = 0; k < 10; ++k) {
                        acc[k] += c * gto_int_moment_1d(mom[k][0], t, P[0], p)
                                    * gto_int_moment_1d(mom[k][1], u, P[1], p)
                                    * gto_int_moment_1d(mom[k][2], v, P[2], p);
                    }
                }
            }
        }
    }
    for (uint32_t k = 0; k < q->num_points; ++k) {
        const double c = q->point_charge[k];
        const double x = q->point_xyz[k * 3 + 0] - O[0];
        const double y = q->point_xyz[k * 3 + 1] - O[1];
        const double z = q->point_xyz[k * 3 + 2] - O[2];
        const double v[10] = { 1, x, y, z, x * x, x * y, x * z, y * y, y * z, z * z };
        for (int i = 0; i < 10; ++i) acc[i] += c * v[i];
    }

    m.charge = acc[0];
    for (int i = 0; i < 3; ++i) m.dipole[i] = acc[1 + i];
    for (int i = 0; i < 6; ++i) m.second[i] = acc[4 + i];
    return m;
}

// ---------------------------------------------------------------------------
// CPU evaluation
// ---------------------------------------------------------------------------

static double gto_int_eval_point(double field[3], const double C[3], const md_gto_int_charges_t* q) {
    double V = 0.0;
    double E[3] = {0, 0, 0};
    double R[GTO_INT_MAX_R];

    uint32_t L = 0;
    for (uint32_t g = 0; g < q->num_gaussians; ++g) {
        gto_int_order_of(q, g, &L);
        const double* P = q->center + g * 3;
        const double* h = q->coeff + q->coeff_offset[g];
        const int Lr = (int)L + (field ? 1 : 0);
        gto_int_hermite_R(R, Lr, q->exponent[g], P[0] - C[0], P[1] - C[1], P[2] - C[2]);

        uint32_t i = 0;
        double s = 0.0, ex = 0.0, ey = 0.0, ez = 0.0;
        for (int deg = 0; deg <= (int)L; ++deg) {
            for (int t = deg; t >= 0; --t) {
                for (int u = deg - t; u >= 0; --u, ++i) {
                    const int v = deg - t - u;
                    s += h[i] * R[i];
                    if (field) {
                        // E = -grad_C V, and d/dC R_tuv(P - C) = -R_{t+1,u,v} etc.
                        ex += h[i] * R[herm_idx(t + 1, u, v)];
                        ey += h[i] * R[herm_idx(t, u + 1, v)];
                        ez += h[i] * R[herm_idx(t, u, v + 1)];
                    }
                }
            }
        }
        V += s;
        E[0] += ex; E[1] += ey; E[2] += ez;
    }

    for (uint32_t k = 0; k < q->num_points; ++k) {
        const double d[3] = { C[0] - q->point_xyz[k * 3 + 0], C[1] - q->point_xyz[k * 3 + 1], C[2] - q->point_xyz[k * 3 + 2] };
        const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if (r2 < 1.0e-20) continue;
        const double inv_r = 1.0 / sqrt(r2);
        const double qr = q->point_charge[k] * inv_r;
        V += qr;
        if (field) {
            const double qr3 = qr * inv_r * inv_r;
            E[0] += qr3 * d[0]; E[1] += qr3 * d[1]; E[2] += qr3 * d[2];
        }
    }

    if (field) { field[0] = E[0]; field[1] = E[1]; field[2] = E[2]; }
    return V;
}

void md_gto_int_potential_xyz(double* out_potential, double* out_field, const float* xyz, size_t num_xyz, size_t xyz_stride,
                              const md_gto_int_charges_t* charges) {
    ASSERT(out_potential);
    ASSERT(xyz);
    ASSERT(charges);
    const size_t stride = xyz_stride ? xyz_stride : sizeof(float) * 3;
    for (size_t i = 0; i < num_xyz; ++i) {
        const float* x = gto_int_xyz(xyz, stride, i);
        const double C[3] = { x[0], x[1], x[2] };
        out_potential[i] = gto_int_eval_point(out_field ? out_field + i * 3 : NULL, C, charges);
    }
}

void md_gto_int_potential_grid_sub(float* out_values, const md_grid_t* grid, const float sample_offset[3],
                                   const int idx_off[3], const int idx_len[3], const md_gto_int_charges_t* charges) {
    ASSERT(out_values);
    ASSERT(grid);
    ASSERT(charges);
    const double so[3] = { sample_offset ? sample_offset[0] : 0.0, sample_offset ? sample_offset[1] : 0.0, sample_offset ? sample_offset[2] : 0.0 };

    for (int k = idx_off[2]; k < idx_off[2] + idx_len[2]; ++k) {
        for (int j = idx_off[1]; j < idx_off[1] + idx_len[1]; ++j) {
            for (int i = idx_off[0]; i < idx_off[0] + idx_len[0]; ++i) {
                const double l[3] = { (i + so[0]) * grid->spacing.elem[0], (j + so[1]) * grid->spacing.elem[1], (k + so[2]) * grid->spacing.elem[2] };
                double C[3];
                for (int a = 0; a < 3; ++a) {
                    // orientation columns are the grid axes in world space
                    C[a] = grid->origin.elem[a] + grid->orientation.elem[0][a] * l[0] + grid->orientation.elem[1][a] * l[1] + grid->orientation.elem[2][a] * l[2];
                }
                const size_t idx = ((size_t)k * grid->dim[1] + j) * grid->dim[0] + i;
                out_values[idx] = (float)gto_int_eval_point(NULL, C, charges);
            }
        }
    }
}

void md_gto_int_potential_grid(float* out_values, const md_grid_t* grid, const float sample_offset[3], const md_gto_int_charges_t* charges) {
    const int off[3] = {0, 0, 0};
    md_gto_int_potential_grid_sub(out_values, grid, sample_offset, off, grid->dim, charges);
}

// ---------------------------------------------------------------------------
// GPU
// ---------------------------------------------------------------------------
#if MD_ENABLE_GPU

#include <eval_gto_int_potential_shaders.inl>

// Must match eval_gto_int_potential.slang.
#define GTO_INT_BOYS_DT    0.1
#define GTO_INT_BOYS_TMAX  30.0
#define GTO_INT_BOYS_NN    16
#define GTO_INT_BOYS_ROWS  ((int)(GTO_INT_BOYS_TMAX / GTO_INT_BOYS_DT) + 2)

// Largest slab of voxels evaluated at once; its fixed point accumulator is 8 bytes per voxel.
#define GTO_INT_MAX_SLAB_VOXELS (1u << 21)
#define GTO_INT_DEFAULT_WORK    ((uint64_t)1 << 31)

/* Mirrors RootArgs in eval_gto_int_potential.slang. */
typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    box_off;
    md_gpu_uint4    box_dim;
    uint32_t        first;
    uint32_t        count;
    uint32_t        coeff_first;
    uint32_t        num_points;
    uint32_t        operation;
    uint32_t        _pad0;
    uint32_t        _pad1;
    uint32_t        _pad2;
    md_gpu_addr_t   gaussians;
    md_gpu_addr_t   coeffs;
    md_gpu_addr_t   boys;
    md_gpu_addr_t   points;
    md_gpu_addr_t   accum;
    md_gpu_storage_tex_t out_tex;
} gto_int_potential_args_t;

static md_gpu_device_t gto_int_gpu_device = NULL;
static md_gpu_kernel_t gto_int_k_order[MD_GTO_INT_MAX_ORDER + 1] = {0};
static md_gpu_kernel_t gto_int_k_resolve = NULL;

void md_gto_int_gpu_initialize(md_gpu_device_t device) {
    gto_int_gpu_device = device;
}

void md_gto_int_gpu_shutdown(void) {
    for (int L = 0; L <= MD_GTO_INT_MAX_ORDER; ++L) {
        if (gto_int_k_order[L]) { md_gpu_kernel_destroy(gto_int_k_order[L]); gto_int_k_order[L] = NULL; }
    }
    if (gto_int_k_resolve) { md_gpu_kernel_destroy(gto_int_k_resolve); gto_int_k_resolve = NULL; }
    gto_int_gpu_device = NULL;
}

static md_gpu_kernel_t gto_int_kernel(md_gpu_kernel_t* slot, md_gpu_kernel_desc_t desc) {
    if (!*slot && gto_int_gpu_device) {
        *slot = md_gpu_kernel_create(gto_int_gpu_device, &desc);
        if (!*slot) MD_LOG_ERROR("md_gto_int: failed to create kernel '%s': %s", desc.label ? desc.label : desc.entry_point, md_gpu_last_error());
    }
    return *slot;
}

// Created on first use, so a session that never sees a g function never compiles the l8 kernel.
static md_gpu_kernel_t gto_int_order_kernel(uint32_t L) {
    switch (L) {
    case 0: return gto_int_kernel(&gto_int_k_order[0], md_shader_eval_gto_int_potential_l0_kernel());
    case 1: return gto_int_kernel(&gto_int_k_order[1], md_shader_eval_gto_int_potential_l1_kernel());
    case 2: return gto_int_kernel(&gto_int_k_order[2], md_shader_eval_gto_int_potential_l2_kernel());
    case 3: return gto_int_kernel(&gto_int_k_order[3], md_shader_eval_gto_int_potential_l3_kernel());
    case 4: return gto_int_kernel(&gto_int_k_order[4], md_shader_eval_gto_int_potential_l4_kernel());
    case 5: return gto_int_kernel(&gto_int_k_order[5], md_shader_eval_gto_int_potential_l5_kernel());
    case 6: return gto_int_kernel(&gto_int_k_order[6], md_shader_eval_gto_int_potential_l6_kernel());
    case 7: return gto_int_kernel(&gto_int_k_order[7], md_shader_eval_gto_int_potential_l7_kernel());
    case 8: return gto_int_kernel(&gto_int_k_order[8], md_shader_eval_gto_int_potential_l8_kernel());
    default: return NULL;
    }
}

typedef struct md_gto_int_gpu_charges {
    md_gpu_addr_t buffer;
    size_t        off_boys;
    size_t        off_gaussians;
    size_t        off_coeffs;
    size_t        off_points;
    size_t        size;
    uint32_t      num_gaussians;
    uint32_t      num_points;
    uint32_t      max_order;
    uint32_t      order_offset[MD_GTO_INT_MAX_ORDER + 2];
    uint32_t      coeff_first[MD_GTO_INT_MAX_ORDER + 1];   // coefficient index of each order's first gaussian
    double        ref[3];   // device positions are stored relative to this point, see below
} md_gto_int_gpu_charges;

md_gto_int_gpu_charges_t md_gto_int_gpu_charges_create(md_gpu_stream_t stream, const md_gto_int_charges_t* q) {
    if (!stream || !q) {
        MD_LOG_ERROR("md_gto_int_gpu_charges_create: invalid input");
        return NULL;
    }

    md_gto_int_gpu_charges* g = md_alloc(md_get_heap_allocator(), sizeof(md_gto_int_gpu_charges));
    if (!g) return NULL;
    MEMSET(g, 0, sizeof(*g));
    g->num_gaussians = q->num_gaussians;
    g->num_points    = q->num_points;
    g->max_order     = q->max_order;
    MEMCPY(g->order_offset, q->order_offset, sizeof(g->order_offset));
    {
        uint32_t run = 0;
        for (uint32_t L = 0; L <= MD_GTO_INT_MAX_ORDER; ++L) {
            g->coeff_first[L] = run;
            run += (q->order_offset[L + 1] - q->order_offset[L]) * md_gto_int_num_hermite(L);
        }
    }

    // Positions go to the device RELATIVE TO THE CENTRE of the distribution, and the launch moves
    // the grid into the same frame in double before rounding. Near a nucleus the potential changes
    // by ~100 au per bohr, so the float resolution of an absolute coordinate (1e-6 bohr at 10 bohr,
    // 1e-5 at 100 - QM files are not always centred) is what limits the result there, not the
    // integrals.
    {
        double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
        for (uint32_t i = 0; i < q->num_gaussians; ++i) for (int k = 0; k < 3; ++k) {
            lo[k] = MIN(lo[k], q->center[i * 3 + k]); hi[k] = MAX(hi[k], q->center[i * 3 + k]);
        }
        for (uint32_t i = 0; i < q->num_points; ++i) for (int k = 0; k < 3; ++k) {
            lo[k] = MIN(lo[k], q->point_xyz[i * 3 + k]); hi[k] = MAX(hi[k], q->point_xyz[i * 3 + k]);
        }
        for (int k = 0; k < 3; ++k) g->ref[k] = (q->num_gaussians + q->num_points) ? 0.5 * (lo[k] + hi[k]) : 0.0;
    }

    // One allocation: Boys table | gaussians (float4) | coefficients | point charges (float4)
    const size_t boys_bytes  = sizeof(float) * GTO_INT_BOYS_ROWS * GTO_INT_BOYS_NN;
    const size_t gauss_bytes = sizeof(float) * 4 * MAX(q->num_gaussians, 1);
    const size_t coeff_bytes = sizeof(float) * MAX(q->num_coeffs, 1);
    const size_t point_bytes = sizeof(float) * 4 * MAX(q->num_points, 1);
    g->off_boys      = 0;
    g->off_gaussians = ALIGN_TO(g->off_boys      + boys_bytes,  256);
    g->off_coeffs    = ALIGN_TO(g->off_gaussians + gauss_bytes, 256);
    g->off_points    = ALIGN_TO(g->off_coeffs    + coeff_bytes, 256);
    g->size          = ALIGN_TO(g->off_points    + point_bytes, 256);

    g->buffer = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, g->size).gpu;
    if (!g->buffer) {
        MD_LOG_ERROR("md_gto_int_gpu_charges_create: failed to allocate %zu bytes", g->size);
        md_free(md_get_heap_allocator(), g, sizeof(*g));
        return NULL;
    }

    char* dst = (char*)md_gpu_upload_begin(stream, g->buffer, g->size);
    if (!dst) {
        md_gpu_free(stream, g->buffer);
        md_free(md_get_heap_allocator(), g, sizeof(*g));
        return NULL;
    }
    MEMSET(dst, 0, g->size);

    // The table the kernel's Boys function interpolates in, from the double precision one here.
    {
        float* boys = (float*)(dst + g->off_boys);
        double F[GTO_INT_BOYS_NN];
        for (int k = 0; k < GTO_INT_BOYS_ROWS; ++k) {
            gto_int_boys(F, GTO_INT_BOYS_NN - 1, k * GTO_INT_BOYS_DT);
            for (int n = 0; n < GTO_INT_BOYS_NN; ++n) boys[k * GTO_INT_BOYS_NN + n] = (float)F[n];
        }
    }
    {
        float* gs = (float*)(dst + g->off_gaussians);
        for (uint32_t i = 0; i < q->num_gaussians; ++i) {
            gs[i * 4 + 0] = (float)(q->center[i * 3 + 0] - g->ref[0]);
            gs[i * 4 + 1] = (float)(q->center[i * 3 + 1] - g->ref[1]);
            gs[i * 4 + 2] = (float)(q->center[i * 3 + 2] - g->ref[2]);
            gs[i * 4 + 3] = (float)q->exponent[i];
        }
        // md_gto_int_charges_init lays the coefficients out contiguously in gaussian order, so the
        // device copy needs no offsets: gaussian i of order L starts at
        // coeff_first[L] + (i - order_offset[L]) * num_hermite(L).
        float* cs = (float*)(dst + g->off_coeffs);
        for (size_t i = 0; i < q->num_coeffs; ++i) cs[i] = (float)q->coeff[i];
        float* ps = (float*)(dst + g->off_points);
        for (uint32_t i = 0; i < q->num_points; ++i) {
            ps[i * 4 + 0] = (float)(q->point_xyz[i * 3 + 0] - g->ref[0]);
            ps[i * 4 + 1] = (float)(q->point_xyz[i * 3 + 1] - g->ref[1]);
            ps[i * 4 + 2] = (float)(q->point_xyz[i * 3 + 2] - g->ref[2]);
            ps[i * 4 + 3] = (float)q->point_charge[i];
        }
    }
    md_gpu_upload_end(stream);
    return g;
}

void md_gto_int_gpu_charges_destroy(md_gpu_stream_t stream, md_gto_int_gpu_charges_t g) {
    if (!g) return;
    if (stream && g->buffer) md_gpu_free(stream, g->buffer);
    md_free(md_get_heap_allocator(), g, sizeof(*g));
}

// Grid index (with the sample offset) to world MINUS ref, column major as the kernel's mul()
// expects (the layout md_gto.c uses for its kernels). The translation is formed in double, so the
// frame shift costs no precision.
static void gto_int_index_to_world(md_gpu_float4x4* out, const md_grid_t* grid, const float so[3], const double ref[3]) {
    float (*m)[4] = (float(*)[4])out->m;
    for (int c = 0; c < 3; ++c) {
        for (int r = 0; r < 3; ++r) m[c][r] = grid->orientation.elem[c][r] * grid->spacing.elem[c];
        m[c][3] = 0.0f;
    }
    for (int r = 0; r < 3; ++r) {
        double t = (double)grid->origin.elem[r] - ref[r];
        for (int c = 0; c < 3; ++c) t += (double)so[c] * (double)grid->orientation.elem[c][r] * (double)grid->spacing.elem[c];
        m[3][r] = (float)t;
    }
    m[3][3] = 1.0f;
}

void md_gto_int_gpu_potential_launch(md_gpu_stream_t stream, const md_gto_int_gpu_potential_desc_t* desc) {
    if (!stream || !desc || !desc->charges || !desc->out_tex || !desc->grid) {
        MD_LOG_ERROR("md_gto_int_gpu_potential_launch: invalid input");
        return;
    }
    const md_gto_int_gpu_charges* g = desc->charges;
    const md_grid_t* grid = desc->grid;
    if (grid->dim[0] <= 0 || grid->dim[1] <= 0 || grid->dim[2] <= 0) return;

    const md_gpu_storage_tex_t out = md_gpu_texture_storage(desc->out_tex, 0);
    if (!out.handle) {
        MD_LOG_ERROR("md_gto_int_gpu_potential_launch: out_tex needs MD_GPU_TEX_STORAGE usage");
        return;
    }
    md_gpu_kernel_t k_resolve = gto_int_kernel(&gto_int_k_resolve, md_shader_eval_gto_int_potential_resolve_kernel());
    if (!k_resolve) {
        MD_LOG_ERROR("md_gto_int_gpu_potential_launch: kernels unavailable (md_gto_int_gpu_initialize not called?)");
        return;
    }
    for (uint32_t L = 0; L <= g->max_order; ++L) {
        if (g->order_offset[L + 1] > g->order_offset[L] && !gto_int_order_kernel(L)) return;
    }

    const uint64_t budget = desc->work_per_dispatch ? desc->work_per_dispatch : GTO_INT_DEFAULT_WORK;
    const uint32_t dx = (uint32_t)grid->dim[0], dy = (uint32_t)grid->dim[1], dz = (uint32_t)grid->dim[2];
    const uint64_t plane = (uint64_t)dx * dy;

    // Slabs of whole z planes: bounded by the accumulator size, and by the budget so that a single
    // gaussian of the highest order present still fits one dispatch.
    uint64_t slab_vox = MIN((uint64_t)GTO_INT_MAX_SLAB_VOXELS, budget / gto_int_cost(g->max_order));
    uint32_t slab_z   = (uint32_t)CLAMP(slab_vox / plane, 1, dz);

    gto_int_potential_args_t a = {0};
    gto_int_index_to_world(&a.index_to_world, grid, desc->sample_offset, g->ref);
    a.num_points = g->num_points;
    a.operation  = (uint32_t)desc->op;
    a.gaussians  = g->buffer + g->off_gaussians;
    a.coeffs     = g->buffer + g->off_coeffs;
    a.boys       = g->buffer + g->off_boys;
    a.points     = g->buffer + g->off_points;
    a.out_tex    = out;

    // Every dispatch below reads what the previous one wrote; say so whatever mode the caller's
    // stream is in, and leave it as it was.
    const md_gpu_ordering_t prev_ordering = md_gpu_stream_ordering(stream);
    md_gpu_stream_set_ordering(stream, MD_GPU_ORDER_IMPLICIT);

    md_gpu_temp_t scope = md_gpu_temp_begin(stream);
    const size_t accum_bytes = sizeof(uint32_t) * 2 * plane * slab_z;
    a.accum = md_gpu_temp_alloc(stream, MD_GPU_MEM_DEVICE, accum_bytes).gpu;
    if (!a.accum) {
        MD_LOG_ERROR("md_gto_int_gpu_potential_launch: failed to allocate %zu bytes of accumulator", accum_bytes);
    } else {
        for (uint32_t z0 = 0; z0 < dz; z0 += slab_z) {
            const uint32_t nz = MIN(slab_z, dz - z0);
            const uint64_t nvox = plane * nz;
            a.box_off = (md_gpu_uint4){ 0, 0, z0, 0 };
            a.box_dim = (md_gpu_uint4){ dx, dy, nz, 0 };

            md_gpu_memset(stream, a.accum, 0, (size_t)(sizeof(uint32_t) * 2 * nvox));

            for (uint32_t L = 0; L <= g->max_order; ++L) {
                const uint32_t begin = g->order_offset[L], end = g->order_offset[L + 1];
                if (begin == end) continue;
                md_gpu_kernel_t k = gto_int_order_kernel(L);
                const uint64_t per = MAX(1, budget / (nvox * gto_int_cost(L)));
                const uint32_t nh  = md_gto_int_num_hermite(L);
                const md_gpu_grid_t grid_groups = md_gpu_grid_for(k, dx, dy, nz);
                for (uint32_t first = begin; first < end; first += (uint32_t)MIN(per, (uint64_t)(end - first))) {
                    a.first       = first;
                    a.count       = (uint32_t)MIN(per, (uint64_t)(end - first));
                    a.coeff_first = g->coeff_first[L] + (first - begin) * nh;
                    md_gpu_launch(stream, k, grid_groups, &a, sizeof(a));
                }
            }

            md_gpu_launch(stream, k_resolve, md_gpu_grid_for(k_resolve, dx, dy, nz), &a, sizeof(a));
        }
    }
    md_gpu_temp_end(stream, scope);
    md_gpu_stream_set_ordering(stream, prev_ordering);
}

#endif
