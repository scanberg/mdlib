#include <md_gto.h>
#include <md_util.h>

#include <core/md_platform.h>

#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_simd.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_intrinsics.h>
#include <md_system.h>

#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

typedef struct {
    float coeff;
    float alpha;
    uint32_t ijkl;
} PGTO;

static inline void world_to_model_matrix(float out_mat[4][4], const md_grid_t* grid) {
    // There is no scaling applied in this transformation, only rotation and translation
    out_mat[0][0] = grid->orientation.elem[0][0];
    out_mat[0][1] = grid->orientation.elem[1][0];
    out_mat[0][2] = grid->orientation.elem[2][0];
    out_mat[0][3] = 0.0f;
    out_mat[1][0] = grid->orientation.elem[0][1];
    out_mat[1][1] = grid->orientation.elem[1][1];
    out_mat[1][2] = grid->orientation.elem[2][1];
    out_mat[1][3] = 0.0f;
    out_mat[2][0] = grid->orientation.elem[0][2];
    out_mat[2][1] = grid->orientation.elem[1][2];
    out_mat[2][2] = grid->orientation.elem[2][2];
    out_mat[2][3] = 0.0f;
    out_mat[3][0] = -grid->orientation.elem[0][0] * grid->origin.elem[0] - grid->orientation.elem[0][1] * grid->origin.elem[1] - grid->orientation.elem[0][2] * grid->origin.elem[2],
    out_mat[3][1] = -grid->orientation.elem[1][0] * grid->origin.elem[0] - grid->orientation.elem[1][1] * grid->origin.elem[1] - grid->orientation.elem[1][2] * grid->origin.elem[2],
    out_mat[3][2] = -grid->orientation.elem[2][0] * grid->origin.elem[0] - grid->orientation.elem[2][1] * grid->origin.elem[1] - grid->orientation.elem[2][2] * grid->origin.elem[2],
    out_mat[3][3] = 1.0f;
}

static inline void index_to_world_matrix(float out_mat[4][4], const md_grid_t* grid, const float sample_offset[3]) {
    out_mat[0][0] = grid->orientation.elem[0][0] * grid->spacing.elem[0];
    out_mat[0][1] = grid->orientation.elem[0][1] * grid->spacing.elem[0];
    out_mat[0][2] = grid->orientation.elem[0][2] * grid->spacing.elem[0];
    out_mat[0][3] = 0.0f;
    out_mat[1][0] = grid->orientation.elem[1][0] * grid->spacing.elem[1];
    out_mat[1][1] = grid->orientation.elem[1][1] * grid->spacing.elem[1];
    out_mat[1][2] = grid->orientation.elem[1][2] * grid->spacing.elem[1];
    out_mat[1][3] = 0.0f;
    out_mat[2][0] = grid->orientation.elem[2][0] * grid->spacing.elem[2];
    out_mat[2][1] = grid->orientation.elem[2][1] * grid->spacing.elem[2];
    out_mat[2][2] = grid->orientation.elem[2][2] * grid->spacing.elem[2];
    out_mat[2][3] = 0.0f;
    out_mat[3][0] = grid->origin.elem[0];
    out_mat[3][1] = grid->origin.elem[1];
    out_mat[3][2] = grid->origin.elem[2];
    out_mat[3][3] = 1.0f;

    // Incorporate a half voxel offset to move to voxel centers
    out_mat[3][0] += sample_offset[0] * (out_mat[0][0] + out_mat[1][0] + out_mat[2][0]);
    out_mat[3][1] += sample_offset[1] * (out_mat[0][1] + out_mat[1][1] + out_mat[2][1]);
    out_mat[3][2] += sample_offset[2] * (out_mat[0][2] + out_mat[1][2] + out_mat[2][2]);
}

// ---------------------------------------------------------------------------
// Spherical→Cartesian expansion tables
// Converts md_gto_basis_t (radial shells, pure radial coefficients) into the
// Cartesian layout used by evaluators.
// Coupling coefficients and Cartesian (l,m,n) tables ported from VeloxChem.
// ---------------------------------------------------------------------------

typedef uint8_t gto_lmn_t[3];

static const gto_lmn_t gto_S_lmn[1]  = {{0,0,0}};
static const gto_lmn_t gto_P_lmn[3]  = {{1,0,0},{0,1,0},{0,0,1}};
static const gto_lmn_t gto_D_lmn[6]  = {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
static const gto_lmn_t gto_F_lmn[10] = {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}};
static const gto_lmn_t gto_G_lmn[15] = {{4,0,0},{3,1,0},{3,0,1},{2,2,0},{2,1,1},{2,0,2},{1,3,0},{1,2,1},{1,1,2},{1,0,3},{0,4,0},{0,3,1},{0,2,2},{0,1,3},{0,0,4}};

#define GTO_d3  3.464101615137754587
#define GTO_f5  1.581138830084189666
#define GTO_f15 7.745966692414833770
#define GTO_f3  1.224744871391589049
#define GTO_g35 (4.0 * 5.916079783099616042)
#define GTO_g17 (4.0 * 4.183300132670377739)
#define GTO_g5  (4.0 * 2.236067977499789696)
#define GTO_g2  (4.0 * 1.581138830084189666)

static const double  gto_S_factors[] = {1.0};
static const uint8_t gto_S_indices[] = {0};
static const uint8_t gto_S_num_fac[] = {1};

static const double  gto_P_factors[] = {1.0, 1.0, 1.0};
static const uint8_t gto_P_indices[] = {1, 2, 0};
static const uint8_t gto_P_offsets[] = {0, 1, 2};
static const uint8_t gto_P_num_fac[] = {1, 1, 1};

static const double  gto_D_factors[] = {GTO_d3, GTO_d3, -1.0, -1.0, 2.0, GTO_d3, 0.5*GTO_d3, -0.5*GTO_d3};
static const uint8_t gto_D_indices[] = {1, 4, 0, 3, 5, 2, 0, 3};
static const uint8_t gto_D_offsets[] = {0, 1, 2, 5, 6};
static const uint8_t gto_D_num_fac[] = {1, 1, 3, 1, 2};

static const double  gto_F_factors[] = {3.0*GTO_f5, -GTO_f5, GTO_f15, 4.0*GTO_f3, -GTO_f3, -GTO_f3, 2.0, -3.0, -3.0, 4.0*GTO_f3, -GTO_f3, -GTO_f3, 0.5*GTO_f15, -0.5*GTO_f15, GTO_f5, -3.0*GTO_f5};
static const uint8_t gto_F_indices[] = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
static const uint8_t gto_F_offsets[] = {0, 2, 3, 6, 9, 12, 14};
static const uint8_t gto_F_num_fac[] = {2, 1, 3, 3, 3, 2, 2};

static const double gto_G_factors[] = {
    GTO_g35, -GTO_g35, 3.0*GTO_g17, -GTO_g17, 6.0*GTO_g5, -GTO_g5, -GTO_g5, 4.0*GTO_g2, -3.0*GTO_g2, -3.0*GTO_g2,
    8.0, 3.0, 3.0, 6.0, -24.0, -24.0, 4.0*GTO_g2, -3.0*GTO_g2, -3.0*GTO_g2, 3.0*GTO_g5,
    -3.0*GTO_g5, -0.5*GTO_g5, 0.5*GTO_g5, GTO_g17, -3.0*GTO_g17, 0.25*GTO_g35, 0.25*GTO_g35, -1.50*GTO_g35};
static const uint8_t gto_G_indices[] = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2, 7, 0, 10, 3};
static const uint8_t gto_G_offsets[] = {0, 2, 4, 7, 10, 16, 19, 23, 25};
static const uint8_t gto_G_num_fac[] = {2, 2, 3, 3, 6, 3, 4, 2, 3};

#undef GTO_d3
#undef GTO_f5
#undef GTO_f15
#undef GTO_f3
#undef GTO_g35
#undef GTO_g17
#undef GTO_g5
#undef GTO_g2

static inline int gto_sph_num_factors(int l, int isph) {
    switch (l) {
    case 0: return gto_S_num_fac[isph];
    case 1: return gto_P_num_fac[isph];
    case 2: return gto_D_num_fac[isph];
    case 3: return gto_F_num_fac[isph];
    case 4: return gto_G_num_fac[isph];
    default: ASSERT(false); return 0;
    }
}

static inline const double* gto_sph_factors(int l, int isph) {
    switch (l) {
    case 0: return gto_S_factors;
    case 1: return gto_P_factors + gto_P_offsets[isph];
    case 2: return gto_D_factors + gto_D_offsets[isph];
    case 3: return gto_F_factors + gto_F_offsets[isph];
    case 4: return gto_G_factors + gto_G_offsets[isph];
    default: ASSERT(false); return NULL;
    }
}

static inline const uint8_t* gto_sph_indices(int l, int isph) {
    switch (l) {
    case 0: return gto_S_indices;
    case 1: return gto_P_indices + gto_P_offsets[isph];
    case 2: return gto_D_indices + gto_D_offsets[isph];
    case 3: return gto_F_indices + gto_F_offsets[isph];
    case 4: return gto_G_indices + gto_G_offsets[isph];
    default: ASSERT(false); return NULL;
    }
}

static inline const gto_lmn_t* gto_cart_lmn(int l) {
    switch (l) {
    case 0: return gto_S_lmn;
    case 1: return gto_P_lmn;
    case 2: return gto_D_lmn;
    case 3: return gto_F_lmn;
    case 4: return gto_G_lmn;
    default: ASSERT(false); return NULL;
    }
}

// ---------------------------------------------------------------------------
// Cartesian AO helpers  (see the AO CONVENTION block in md_gto.h)
// ---------------------------------------------------------------------------

// (2n-1)!! with (-1)!! == 1.
static inline double gto_double_factorial_odd(int n) {
    double r = 1.0;
    for (int k = 2 * n - 1; k > 1; k -= 2) r *= (double)k;
    return r;
}

double md_gto_cart_norm_factor(int i, int j, int k) {
    double d = gto_double_factorial_odd(i) * gto_double_factorial_odd(j) * gto_double_factorial_odd(k);
    return 1.0 / sqrt(d);
}

bool md_gto_cart_ijk(int* out_i, int* out_j, int* out_k, uint32_t l, uint32_t cart_idx) {
    if (l > MD_GTO_MAX_ANGULAR_MOMENTUM || cart_idx >= md_gto_num_cart_ao(l)) return false;
    const gto_lmn_t* lmn = gto_cart_lmn((int)l);
    if (out_i) *out_i = lmn[cart_idx][0];
    if (out_j) *out_j = lmn[cart_idx][1];
    if (out_k) *out_k = lmn[cart_idx][2];
    return true;
}

// Per-(l, cart_idx) normalization factor, matching the lmn tables above.
static inline double gto_cart_norm(int l, int ic) {
    const gto_lmn_t* lmn = gto_cart_lmn(l);
    return md_gto_cart_norm_factor(lmn[ic][0], lmn[ic][1], lmn[ic][2]);
}

size_t md_gto_basis_num_ao(const md_gto_basis_t* basis) {
    if (!basis) return 0;
    size_t n = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        n += md_gto_num_cart_ao(basis->shells[si].l);
    }
    return n;
}

size_t md_gto_basis_ao_to_atom(uint32_t* out_atom_idx, const md_gto_basis_t* basis) {
    if (!out_atom_idx || !basis) return 0;
    size_t k = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const uint32_t n = md_gto_num_cart_ao(basis->shells[si].l);
        for (uint32_t q = 0; q < n; ++q) {
            out_atom_idx[k++] = basis->shells[si].atom_idx;
        }
    }
    return k;
}

size_t md_gto_basis_num_sph_ao(const md_gto_basis_t* basis) {
    if (!basis) return 0;
    size_t n = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        n += md_gto_num_sph_ao(basis->shells[si].l);
    }
    return n;
}

// Scatter one shell's spherical coefficients onto its Cartesian AOs.
// c_cart[m] = ( sum_isph c_sph[isph] * coupling(isph -> m) ) / f(m)
// The division by f undoes the normalization that the evaluator re-applies, so
// the represented function is unchanged.
static void gto_shell_sph_to_cart(double* out_cart, const double* in_sph, int l) {
    const int ncart = (int)md_gto_num_cart_ao((uint32_t)l);
    const int nsph  = (int)md_gto_num_sph_ao((uint32_t)l);

    for (int m = 0; m < ncart; ++m) out_cart[m] = 0.0;

    for (int isph = 0; isph < nsph; isph++) {
        const double   c      = in_sph[isph];
        const int      ncomp  = gto_sph_num_factors(l, isph);
        const double*  fcarts = gto_sph_factors(l, isph);
        const uint8_t* sidx   = gto_sph_indices(l, isph);
        for (int ic = 0; ic < ncomp; ic++) {
            out_cart[sidx[ic]] += c * fcarts[ic];
        }
    }

    for (int m = 0; m < ncart; ++m) out_cart[m] /= gto_cart_norm(l, m);
}

size_t md_gto_sph_to_cart_vector(double* out_cart, const double* in_sph, const md_gto_basis_t* basis) {
    if (!out_cart || !in_sph || !basis) return 0;
    size_t si_off = 0, ci_off = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const int l = (int)basis->shells[si].l;
        if (l > MD_GTO_MAX_ANGULAR_MOMENTUM) return 0;
        gto_shell_sph_to_cart(out_cart + ci_off, in_sph + si_off, l);
        si_off += md_gto_num_sph_ao((uint32_t)l);
        ci_off += md_gto_num_cart_ao((uint32_t)l);
    }
    return ci_off;
}

#define GTO_MAX_SPH_PER_SHELL  (2 * MD_GTO_MAX_ANGULAR_MOMENTUM + 1)
#define GTO_MAX_CART_PER_SHELL (((MD_GTO_MAX_ANGULAR_MOMENTUM + 1) * (MD_GTO_MAX_ANGULAR_MOMENTUM + 2)) / 2)
#define GTO_MAX_FACTORS        8

// The nonzero entries of one column of the per-shell transform, i.e. everything
// spherical function 'isph' contributes to:
//     c_cart[out_idx[t]] += out_w[t] * c_sph[isph]
// Returns the number of entries. A d function touches at most 3 monomials out of 6,
// so iterating these instead of a dense block skips most of the arithmetic.
static int gto_sph_cart_entries(uint8_t* out_idx, double* out_w, int l, int isph) {
    const int      ncomp  = gto_sph_num_factors(l, isph);
    const double*  fcarts = gto_sph_factors(l, isph);
    const uint8_t* sidx   = gto_sph_indices(l, isph);

    for (int t = 0; t < ncomp; ++t) {
        out_idx[t] = sidx[t];
        out_w[t]   = fcarts[t] / gto_cart_norm(l, sidx[t]);
    }
    return ncomp;
}

// M_cart = T^T * M_sph * T.
//
// T is block diagonal (a shell's spherical functions only touch that same shell's
// Cartesian AOs) and each block is sparse (8 nonzeros out of 5x6 for d, 16 of 7x10
// for f). Both are exploited:
//
//   pass 1  A = M[rows of shell I] * T     -- sparse scatter, one shell-row block
//   pass 2  out[rows of shell I] = T_I^T * A  -- AXPY over full rows
//
// Cost is ~3*N^2 rather than the ~30*N^2 of dense small-block multiplies, and both
// inner loops are contiguous, which matters a lot in unoptimized builds. Scratch is
// bounded by (2l+1) rows, a few hundred KB, independent of N.
//
// Symmetry is NOT assumed, so this stays correct for a non-symmetric input.
size_t md_gto_sph_to_cart_matrix(double* out_cart, const double* in_sph, const md_gto_basis_t* basis) {
    if (!out_cart || !in_sph || !basis) return 0;

    const size_t n_sph  = md_gto_basis_num_sph_ao(basis);
    const size_t n_cart = md_gto_basis_num_ao(basis);
    if (n_sph == 0 || n_cart == 0) return 0;

    for (uint32_t si = 0; si < basis->num_shells; si++) {
        if (basis->shells[si].l > MD_GTO_MAX_ANGULAR_MOMENTUM) return 0;
    }

    md_temp_scope_t temp = md_temp_begin();

    size_t* sph_off  = (size_t*)md_temp_alloc(temp, sizeof(size_t) * basis->num_shells);
    size_t* cart_off = (size_t*)md_temp_alloc(temp, sizeof(size_t) * basis->num_shells);
    // At most (2l+1) rows of the intermediate are live at a time.
    double* A        = (double*)md_temp_alloc(temp, sizeof(double) * GTO_MAX_SPH_PER_SHELL * n_cart);

    if (!sph_off || !cart_off || !A) {
        MD_LOG_ERROR("Failed to allocate scratch for spherical to Cartesian matrix conversion");
        md_temp_end(temp);
        return 0;
    }

    {
        size_t so = 0, co = 0;
        for (uint32_t si = 0; si < basis->num_shells; si++) {
            sph_off[si]  = so;
            cart_off[si] = co;
            so += md_gto_num_sph_ao(basis->shells[si].l);
            co += md_gto_num_cart_ao(basis->shells[si].l);
        }
    }

    uint8_t ent_idx[GTO_MAX_FACTORS];
    double  ent_w  [GTO_MAX_FACTORS];

    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const int li     = (int)basis->shells[si].l;
        const int nsph_i = (int)md_gto_num_sph_ao((uint32_t)li);
        const int ncrt_i = (int)md_gto_num_cart_ao((uint32_t)li);

        // Pass 1: right-multiply this shell's rows of M by T.
        // gto_shell_sph_to_cart() writes (and zeroes) exactly one shell block, so the
        // shell loop below covers every column of A without a separate clear.
        for (int p = 0; p < nsph_i; ++p) {
            const double* m_row = in_sph + (sph_off[si] + (size_t)p) * n_sph;
            double*       a_row = A + (size_t)p * n_cart;
            for (uint32_t sj = 0; sj < basis->num_shells; sj++) {
                gto_shell_sph_to_cart(a_row + cart_off[sj], m_row + sph_off[sj], (int)basis->shells[sj].l);
            }
        }

        // Pass 2: left-multiply by T_I^T, accumulating into the output rows.
        for (int a = 0; a < ncrt_i; ++a) {
            MEMSET(out_cart + (cart_off[si] + (size_t)a) * n_cart, 0, sizeof(double) * n_cart);
        }

        for (int p = 0; p < nsph_i; ++p) {
            const double* a_row = A + (size_t)p * n_cart;
            const int     ncomp = gto_sph_cart_entries(ent_idx, ent_w, li, p);

            for (int t = 0; t < ncomp; ++t) {
                double*      o_row = out_cart + (cart_off[si] + (size_t)ent_idx[t]) * n_cart;
                const double w     = ent_w[t];
                for (size_t c = 0; c < n_cart; ++c) {
                    o_row[c] += w * a_row[c];
                }
            }
        }
    }

    md_temp_end(temp);
    return n_cart;
}

// ---------------------------------------------------------------------------------------------
// Orbitals, spherical -> Cartesian, in blocks of orbitals
// ---------------------------------------------------------------------------------------------
// It is a transposition as much as a conversion (an AO major file has to become rows of orbitals), so
// it is blocked on both sides. A block is GTO_S2C_NB orbitals; per block the shells are taken a tile
// at a time, a tile being the shells that fit in GTO_S2C_CW Cartesian AOs. Per shell its spherical
// coefficients of the block are read - a run of GTO_S2C_NB values each when the input is contiguous
// along the orbitals - transformed into its Cartesian ones and put into the tile; the full tile then
// goes out as one run per orbital, with non-temporal stores, since nothing reads the output soon. So
// both sides move in runs of a KB or two rather than a value per row, and the output is not read for
// ownership first. Every input value is read once and every output value written once. The threads
// take blocks from a counter.
//
// The 9072 x 9072 coefficients of the protein test case (s, p and d shells), on a 2 core VM: about
// 0.3 s. Writing each shell's values straight into the rows of 64 orbitals instead of through a tile
// took twice that; the separate permutation and transposition this replaced took 2.0 s, and the
// conversion after them 1.5 s.
#define GTO_S2C_NB          256
#define GTO_S2C_CW          128
#define GTO_S2C_MAX_THREADS 64

typedef struct gto_s2c_entries_t {
    int     num;
    uint8_t idx[GTO_MAX_FACTORS];
    double  w  [GTO_MAX_FACTORS];
} gto_s2c_entries_t;

typedef struct gto_s2c_ctx_t {
    double*                  out;
    const double*            in;
    size_t                   mo_stride, ao_stride;
    const int*               perm;
    size_t                   num_mo, n_cart;
    const md_gto_shell_t*    shells;
    uint32_t                 num_shells;
    const uint32_t*          sph_off;
    const uint32_t*          cart_off;  // [num_shells + 1]
    const gto_s2c_entries_t* ent;       // [l][isph], GTO_MAX_SPH_PER_SHELL per l
    volatile uint64_t        next;
    volatile uint64_t        failed;
} gto_s2c_ctx_t;

// n values from src to dst, which is written and not read: non-temporal where it is 16 byte aligned
static inline void gto_s2c_stream(double* dst, const double* src, size_t n) {
    size_t i = 0;
    for (; i < n && ((uintptr_t)(dst + i) & 15u); ++i) dst[i] = src[i];
    for (; i + 2 <= n; i += 2) md_mm_stream_pd(dst + i, md_mm_loadu_pd(src + i));
    for (; i < n; ++i) dst[i] = src[i];
}

static void gto_s2c_work(void* user) {
    gto_s2c_ctx_t* c = (gto_s2c_ctx_t*)user;
    const uint64_t num_blocks = (c->num_mo + GTO_S2C_NB - 1) / GTO_S2C_NB;
    double* acc  = (double*)malloc(sizeof(double) * GTO_MAX_CART_PER_SHELL * GTO_S2C_NB);
    double* x    = (double*)malloc(sizeof(double) * GTO_S2C_NB);
    double* tile = (double*)malloc(sizeof(double) * GTO_S2C_NB * GTO_S2C_CW);
    if (!acc || !x || !tile) {
        md_atomic_fetch_add64(&c->failed, 1);
    } else for (;;) {
        const uint64_t b = md_atomic_fetch_add64(&c->next, 1);
        if (b >= num_blocks) break;
        const size_t m0 = (size_t)b * GTO_S2C_NB;
        const size_t nb = MIN((size_t)GTO_S2C_NB, c->num_mo - m0);
        for (uint32_t s0 = 0; s0 < c->num_shells; ) {
            // the shells of this tile: at least one (no shell has more than GTO_S2C_CW Cartesian AOs)
            const uint32_t c0 = c->cart_off[s0];
            uint32_t s1 = s0 + 1;
            while (s1 < c->num_shells && c->cart_off[s1 + 1] - c0 <= GTO_S2C_CW) ++s1;
            for (uint32_t si = s0; si < s1; ++si) {
                const uint32_t l    = c->shells[si].l;
                const int      nsph = (int)md_gto_num_sph_ao(l);
                const int      ncrt = (int)md_gto_num_cart_ao(l);
                MEMSET(acc, 0, sizeof(double) * (size_t)ncrt * GTO_S2C_NB);
                for (int s = 0; s < nsph; ++s) {
                    const size_t  a   = (size_t)c->sph_off[si] + (size_t)s;
                    const double* src = c->in + (size_t)(c->perm ? c->perm[a] : (int)a) * c->ao_stride + m0 * c->mo_stride;
                    const double* xs  = src;
                    if (c->mo_stride != 1) {
                        for (size_t j = 0; j < nb; ++j) x[j] = src[j * c->mo_stride];
                        xs = x;
                    }
                    const gto_s2c_entries_t* e = c->ent + l * GTO_MAX_SPH_PER_SHELL + (uint32_t)s;
                    for (int t = 0; t < e->num; ++t) {
                        const double w   = e->w[t];
                        double*      row = acc + (size_t)e->idx[t] * GTO_S2C_NB;
                        for (size_t j = 0; j < nb; ++j) row[j] += w * xs[j];
                    }
                }
                double* tl = tile + (c->cart_off[si] - c0);
                for (size_t j = 0; j < nb; ++j) {
                    for (int q = 0; q < ncrt; ++q) tl[j * GTO_S2C_CW + (size_t)q] = acc[(size_t)q * GTO_S2C_NB + j];
                }
            }
            const size_t cw = c->cart_off[s1] - c0;
            for (size_t j = 0; j < nb; ++j) gto_s2c_stream(c->out + (m0 + j) * c->n_cart + c0, tile + j * GTO_S2C_CW, cw);
            s0 = s1;
        }
    }
    md_mm_sfence();
    free(tile);
    free(x);
    free(acc);
}

size_t md_gto_sph_to_cart_orbitals(double* out_cart, const double* in_sph, size_t mo_stride, size_t ao_stride, const int* ao_perm, size_t num_mo, const md_gto_basis_t* basis) {
    if (!out_cart || !in_sph || !basis || num_mo == 0) return 0;
    const size_t n_sph  = md_gto_basis_num_sph_ao(basis);
    const size_t n_cart = md_gto_basis_num_ao(basis);
    if (n_sph == 0 || n_cart == 0) return 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        if (basis->shells[si].l > MD_GTO_MAX_ANGULAR_MOMENTUM) return 0;
    }

    gto_s2c_entries_t ent[(MD_GTO_MAX_ANGULAR_MOMENTUM + 1) * GTO_MAX_SPH_PER_SHELL];
    MEMSET(ent, 0, sizeof(ent));
    for (int l = 0; l <= MD_GTO_MAX_ANGULAR_MOMENTUM; ++l) {
        for (int s = 0; s < (int)md_gto_num_sph_ao((uint32_t)l); ++s) {
            gto_s2c_entries_t* e = ent + l * GTO_MAX_SPH_PER_SHELL + s;
            e->num = gto_sph_cart_entries(e->idx, e->w, l, s);
        }
    }

    uint32_t* off = (uint32_t*)malloc(sizeof(uint32_t) * (2 * (size_t)basis->num_shells + 1));
    if (!off) {
        MD_LOG_ERROR("md_gto_sph_to_cart_orbitals: out of memory");
        return 0;
    }
    uint32_t* sph_off  = off;
    uint32_t* cart_off = off + basis->num_shells;
    {
        uint32_t so = 0, co = 0;
        for (uint32_t si = 0; si < basis->num_shells; si++) {
            sph_off[si]  = so;
            cart_off[si] = co;
            so += md_gto_num_sph_ao(basis->shells[si].l);
            co += md_gto_num_cart_ao(basis->shells[si].l);
        }
        cart_off[basis->num_shells] = co;
    }

    gto_s2c_ctx_t ctx = {
        .out = out_cart, .in = in_sph, .mo_stride = mo_stride, .ao_stride = ao_stride, .perm = ao_perm,
        .num_mo = num_mo, .n_cart = n_cart, .shells = basis->shells, .num_shells = basis->num_shells,
        .sph_off = sph_off, .cart_off = cart_off, .ent = ent,
    };
    // one thread per core, the calling thread one of them, for anything worth it
    const uint64_t num_blocks = (num_mo + GTO_S2C_NB - 1) / GTO_S2C_NB;
    int num_threads = 1;
    md_os_sys_info_t info;
    if ((double)num_mo * (double)n_cart > 4.0e6 && md_os_sys_info_query(&info)) {
        num_threads = info.num_physical_cores > 0 ? info.num_physical_cores : info.num_virtual_cores;
    }
    num_threads = (int)CLAMP((uint64_t)num_threads, (uint64_t)1, MIN((uint64_t)GTO_S2C_MAX_THREADS, num_blocks));
    md_thread_t* threads[GTO_S2C_MAX_THREADS] = {0};
    for (int i = 1; i < num_threads; ++i) threads[i] = md_thread_create(gto_s2c_work, &ctx);
    gto_s2c_work(&ctx);
    for (int i = 1; i < num_threads; ++i) if (threads[i]) md_thread_join(threads[i]);
    const bool ok = md_atomic_load64(&ctx.failed) == 0 || md_atomic_load64(&ctx.next) >= num_blocks;
    if (!ok) MD_LOG_ERROR("md_gto_sph_to_cart_orbitals: out of memory");

    free(off);
    return ok ? num_mo : 0;
}

// ---------------------------------------------------------------------------------------------
// The density matrix of weighted orbitals, D = C^T W C
// ---------------------------------------------------------------------------------------------
// D in tiles of GTO_DM_NB x GTO_DM_NB, the upper triangle of them computed and mirrored on writing,
// each tile over panels of GTO_DM_KC orbitals: per panel the tile's rows (times the weights) and
// columns are packed into micro panels, and GTO_DM_MR x GTO_DM_NR blocks are summed in registers.
// The threads take tiles from a shared counter.
#define GTO_DM_MR 6
#define GTO_DM_NR 8
#define GTO_DM_NB 96            // a multiple of MR and NR
#define GTO_DM_KC 256
#define GTO_DM_MAX_THREADS 64

typedef struct gto_dm_ctx_t {
    double*         D;
    const double*   C;
    size_t          mo_stride, ao_stride;
    const uint32_t* mo;         // the orbitals of nonzero weight
    const double*   w;          // and their weights
    size_t          num_mo, num_ao;
    const uint32_t* tiles;      // (I, J) per tile, I <= J
    uint64_t        num_tiles;
    volatile uint64_t next;
    volatile uint64_t failed;
} gto_dm_ctx_t;

// Rows r0 .. r0 + nr - 1 (AOs) of panel k0 .. k0 + kc - 1 (orbitals), R at a time: P[(b kc + k) R + m]
// is AO r0 + b R + m of orbital k (times its weight), zero past nr.
static void gto_dm_pack(double* P, const gto_dm_ctx_t* c, size_t k0, size_t kc, size_t r0, size_t nr, size_t R, bool weighted) {
    const size_t nb = (nr + R - 1) / R;
    for (size_t b = 0; b < nb; ++b) {
        double* p = P + b * kc * R;
        const size_t m1 = MIN(R, nr - b * R);
        for (size_t k = 0; k < kc; ++k) {
            const double  s   = weighted ? c->w[k0 + k] : 1.0;
            const double* src = c->C + (size_t)c->mo[k0 + k] * c->mo_stride + (r0 + b * R) * c->ao_stride;
            size_t m = 0;
            for (; m < m1; ++m) p[k * R + m] = s * src[m * c->ao_stride];
            for (; m < R; ++m)  p[k * R + m] = 0.0;
        }
    }
}

// T[m][q] += sum_k a[k][m] b[k][q] for a 6 x 8 block (rows ldt apart)
static void gto_dm_micro(size_t kc, const double* a, const double* b, double* t, size_t ldt) {
    md_256d c00 = md_mm256_setzero_pd(), c01 = md_mm256_setzero_pd();
    md_256d c10 = md_mm256_setzero_pd(), c11 = md_mm256_setzero_pd();
    md_256d c20 = md_mm256_setzero_pd(), c21 = md_mm256_setzero_pd();
    md_256d c30 = md_mm256_setzero_pd(), c31 = md_mm256_setzero_pd();
    md_256d c40 = md_mm256_setzero_pd(), c41 = md_mm256_setzero_pd();
    md_256d c50 = md_mm256_setzero_pd(), c51 = md_mm256_setzero_pd();
    for (size_t k = 0; k < kc; ++k) {
        const md_256d b0 = md_mm256_loadu_pd(b);
        const md_256d b1 = md_mm256_loadu_pd(b + 4);
        md_256d x;
        x = md_mm256_set1_pd(a[0]); c00 = md_mm256_fmadd_pd(x, b0, c00); c01 = md_mm256_fmadd_pd(x, b1, c01);
        x = md_mm256_set1_pd(a[1]); c10 = md_mm256_fmadd_pd(x, b0, c10); c11 = md_mm256_fmadd_pd(x, b1, c11);
        x = md_mm256_set1_pd(a[2]); c20 = md_mm256_fmadd_pd(x, b0, c20); c21 = md_mm256_fmadd_pd(x, b1, c21);
        x = md_mm256_set1_pd(a[3]); c30 = md_mm256_fmadd_pd(x, b0, c30); c31 = md_mm256_fmadd_pd(x, b1, c31);
        x = md_mm256_set1_pd(a[4]); c40 = md_mm256_fmadd_pd(x, b0, c40); c41 = md_mm256_fmadd_pd(x, b1, c41);
        x = md_mm256_set1_pd(a[5]); c50 = md_mm256_fmadd_pd(x, b0, c50); c51 = md_mm256_fmadd_pd(x, b1, c51);
        a += GTO_DM_MR;
        b += GTO_DM_NR;
    }
#define GTO_DM_ADD(row, lo, hi) \
    md_mm256_storeu_pd(t + (row) * ldt,     md_mm256_add_pd(md_mm256_loadu_pd(t + (row) * ldt),     lo)); \
    md_mm256_storeu_pd(t + (row) * ldt + 4, md_mm256_add_pd(md_mm256_loadu_pd(t + (row) * ldt + 4), hi))
    GTO_DM_ADD(0, c00, c01);
    GTO_DM_ADD(1, c10, c11);
    GTO_DM_ADD(2, c20, c21);
    GTO_DM_ADD(3, c30, c31);
    GTO_DM_ADD(4, c40, c41);
    GTO_DM_ADD(5, c50, c51);
#undef GTO_DM_ADD
}

static void gto_dm_work(void* user) {
    gto_dm_ctx_t* c = (gto_dm_ctx_t*)user;
    const size_t NB = GTO_DM_NB, KC = GTO_DM_KC, n = c->num_ao;
    double* A = (double*)malloc(sizeof(double) * KC * NB);
    double* B = (double*)malloc(sizeof(double) * KC * NB);
    double* T = (double*)malloc(sizeof(double) * NB * NB);
    if (!A || !B || !T) {
        md_atomic_fetch_add64(&c->failed, 1);
    } else for (;;) {
        const uint64_t t = md_atomic_fetch_add64(&c->next, 1);
        if (t >= c->num_tiles) break;
        const size_t i0 = (size_t)c->tiles[2 * t] * NB, j0 = (size_t)c->tiles[2 * t + 1] * NB;
        const size_t ni = MIN(NB, n - i0), nj = MIN(NB, n - j0);
        MEMSET(T, 0, sizeof(double) * NB * NB);
        for (size_t k0 = 0; k0 < c->num_mo; k0 += KC) {
            const size_t kc = MIN(KC, c->num_mo - k0);
            gto_dm_pack(A, c, k0, kc, i0, ni, GTO_DM_MR, true);
            gto_dm_pack(B, c, k0, kc, j0, nj, GTO_DM_NR, false);
            for (size_t ib = 0; ib * GTO_DM_MR < ni; ++ib) {
                for (size_t jb = 0; jb * GTO_DM_NR < nj; ++jb) {
                    gto_dm_micro(kc, A + ib * kc * GTO_DM_MR, B + jb * kc * GTO_DM_NR, T + ib * GTO_DM_MR * NB + jb * GTO_DM_NR, NB);
                }
            }
        }
        // the upper triangle, mirrored (a diagonal tile's lower half is not taken from T, so that D
        // is exactly symmetric)
        for (size_t a = 0; a < ni; ++a) {
            for (size_t b = (i0 == j0 ? a : 0); b < nj; ++b) {
                const double v = T[a * NB + b];
                c->D[(i0 + a) * n + j0 + b] = v;
                c->D[(j0 + b) * n + i0 + a] = v;
            }
        }
    }
    free(T);
    free(B);
    free(A);
}

bool md_gto_density_from_orbitals(double* D, const double* C, size_t mo_stride, size_t ao_stride, const double* weight, size_t num_mo, size_t num_ao) {
    if (!D || !C || !weight || num_ao == 0) return false;
    uint32_t* mo = (uint32_t*)malloc(sizeof(uint32_t) * (num_mo + 1));
    double*   w  = (double*)malloc(sizeof(double) * (num_mo + 1));
    const size_t nt = (num_ao + GTO_DM_NB - 1) / GTO_DM_NB;
    const uint64_t num_tiles = (uint64_t)nt * (nt + 1) / 2;
    uint32_t* tiles = (uint32_t*)malloc(sizeof(uint32_t) * 2 * num_tiles);
    if (!mo || !w || !tiles) {
        free(tiles); free(w); free(mo);
        MD_LOG_ERROR("md_gto_density_from_orbitals: out of memory");
        return false;
    }
    size_t nk = 0;
    for (size_t k = 0; k < num_mo; ++k) {
        if (weight[k] == 0.0) continue;
        mo[nk] = (uint32_t)k;
        w[nk++] = weight[k];
    }
    uint64_t t = 0;
    for (size_t J = 0; J < nt; ++J) {
        for (size_t I = 0; I <= J; ++I) {
            tiles[2 * t] = (uint32_t)I;
            tiles[2 * t + 1] = (uint32_t)J;
            ++t;
        }
    }

    bool ok = true;
    if (nk == 0) {
        MEMSET(D, 0, sizeof(double) * num_ao * num_ao);
    } else {
        gto_dm_ctx_t ctx = {
            .D = D, .C = C, .mo_stride = mo_stride, .ao_stride = ao_stride, .mo = mo, .w = w,
            .num_mo = nk, .num_ao = num_ao, .tiles = tiles, .num_tiles = num_tiles,
        };
        // one thread per core, the calling thread one of them; small products on the calling thread alone
        int num_threads = 1;
        md_os_sys_info_t info;
        if ((double)nk * (double)num_ao * (double)num_ao > 4.0e8 && md_os_sys_info_query(&info)) {
            num_threads = info.num_physical_cores > 0 ? info.num_physical_cores : info.num_virtual_cores;
        }
        num_threads = (int)CLAMP((uint64_t)num_threads, (uint64_t)1, MIN((uint64_t)GTO_DM_MAX_THREADS, num_tiles));
        md_thread_t* threads[GTO_DM_MAX_THREADS] = {0};
        for (int i = 1; i < num_threads; ++i) threads[i] = md_thread_create(gto_dm_work, &ctx);
        gto_dm_work(&ctx);
        for (int i = 1; i < num_threads; ++i) if (threads[i]) md_thread_join(threads[i]);
        ok = md_atomic_load64(&ctx.failed) == 0 || md_atomic_load64(&ctx.next) >= num_tiles;
        if (!ok) MD_LOG_ERROR("md_gto_density_from_orbitals: out of memory");
    }
    free(tiles);
    free(w);
    free(mo);
    return ok;
}

// Count the number of CGTOs and PGTOs that will result from expanding a basis.
// One CGTO per Cartesian AO, one PGTO per (Cartesian AO, primitive) pair.
static void gto_basis_count(uint32_t* out_num_cgtos, uint32_t* out_num_pgtos,
    const md_gto_basis_t* basis)
{
    uint32_t nc = 0, np = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        uint32_t ncart  = md_gto_num_cart_ao(basis->shells[si].l);
        uint32_t nprims = basis->shells[si].num_primitives;
        nc += ncart;
        np += ncart * nprims;
    }
    *out_num_cgtos = nc;
    *out_num_pgtos = np;
}

size_t md_gto_basis_num_atoms(const md_gto_basis_t* basis) {
    if (!basis) return 0;
    uint32_t max_atom_idx = 0;
    for (uint32_t i = 0; i < basis->num_shells; ++i) {
        max_atom_idx = MAX(max_atom_idx, basis->shells[i].atom_idx);
    }
    return basis->num_shells ? (size_t)(max_atom_idx + 1) : 0;
}

static uint32_t gto_basis_num_atoms(const md_gto_basis_t* basis) {
    return (uint32_t)md_gto_basis_num_atoms(basis);
}

// Expand basis metadata into CGTO/PGTO arrays that do not depend on atom coordinates.
// out_cgto_atom_idx: [num_cgtos], one parent atom index per CGTO.
// out_cgto_r:        [num_cgtos], max primitive radius for screening.
// out_cgto_off_len:  [num_cgtos * 2], offset/length into out_pgto.
// out_pgto:          [num_pgtos], radial+angular primitive data.
static void gto_expand_basis_gpu_meta(
    uint32_t* out_cgto_atom_idx, float* out_cgto_r, uint32_t* out_cgto_off_len, PGTO* out_pgto,
    const md_gto_basis_t* basis, double cutoff)
{
    uint32_t ci = 0, pi = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        int l      = (int)shell->l;
        int ncart  = (int)md_gto_num_cart_ao(shell->l);
        int nprims = (int)shell->num_primitives;
        uint32_t   prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);

        for (int ic = 0; ic < ncart; ic++) {
            const int    lx = lmn[ic][0];
            const int    ly = lmn[ic][1];
            const int    lz = lmn[ic][2];
            const double nrm = gto_cart_norm(l, ic);

            uint32_t pi_beg = pi;
            double max_r = 0.0;
            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                float coeff_val = (float)(coef1 * nrm);
                double radius = md_gto_compute_radius_of_influence(lx, ly, lz, (double)coeff_val, (double)alpha, cutoff);
                max_r = MAX(max_r, radius);

                PGTO pgto = {
                    .coeff = coeff_val,
                    .alpha = alpha,
                    .ijkl = md_gto_pack_ijkl(lx, ly, lz, l),
                };
                out_pgto[pi++] = pgto;
            }

            out_cgto_atom_idx[ci] = shell->atom_idx;
            out_cgto_r[ci] = (float)max_r;
            out_cgto_off_len[2 * ci + 0] = pi_beg;
            out_cgto_off_len[2 * ci + 1] = pi - pi_beg;
            ci++;
        }
    }
}

static uint32_t gto_basis_num_shell_prims(const md_gto_basis_t* basis) {
    uint32_t n = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        n += basis->shells[si].num_primitives;
    }
    return n;
}

// Shell-level device data for the shell and GEMM kernels.
// out_shell_atom_idx: [num_shells]
// out_shell_r:        [num_shells] max radius of influence over all Cartesian components and
//                     primitives of the shell (the same per-primitive radii as the per-AO data)
// out_shell_info:     [num_shells * 4] (prim_offset, num_prims, l, ao_offset)
// out_shell_prim:     [num_shell_prims * 2] (coeff, alpha * log2(e)); the per-component
//                     normalisation f(i,j,k) is applied in the shader.
static void gto_expand_basis_gpu_shells(uint32_t* out_shell_atom_idx, float* out_shell_r,
    uint32_t* out_shell_info, float* out_shell_prim, const md_gto_basis_t* basis, double cutoff)
{
    const double LOG2E = 1.4426950408889634;
    uint32_t ao_off = 0, prim_off = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        const int l      = (int)shell->l;
        const int ncart  = (int)md_gto_num_cart_ao(shell->l);
        const int nprims = (int)shell->num_primitives;
        const uint32_t prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);

        double max_r = 0.0;
        for (int ic = 0; ic < ncart; ic++) {
            const double nrm = gto_cart_norm(l, ic);
            for (int ip = 0; ip < nprims; ip++) {
                const float alpha = basis->alpha[prim_base + ip];
                const float coeff = (float)(basis->coeff[prim_base + ip] * nrm);
                const double r = md_gto_compute_radius_of_influence(lmn[ic][0], lmn[ic][1], lmn[ic][2], (double)coeff, (double)alpha, cutoff);
                max_r = MAX(max_r, r);
            }
        }
        for (int ip = 0; ip < nprims; ip++) {
            out_shell_prim[2 * (prim_off + ip) + 0] = basis->coeff[prim_base + ip];
            out_shell_prim[2 * (prim_off + ip) + 1] = (float)((double)basis->alpha[prim_base + ip] * LOG2E);
        }
        out_shell_atom_idx[si] = shell->atom_idx;
        out_shell_r[si] = (float)max_r;
        out_shell_info[4 * si + 0] = prim_off;
        out_shell_info[4 * si + 1] = (uint32_t)nprims;
        out_shell_info[4 * si + 2] = (uint32_t)l;
        out_shell_info[4 * si + 3] = ao_off;
        prim_off += (uint32_t)nprims;
        ao_off   += (uint32_t)ncart;
    }
}

// Expand md_gto_basis_t into the flat Cartesian SoA arrays expected by the GPU shader.
// Spherical-to-Cartesian coupling factors (fcarts) are baked into pgto_coeff here.
// pgto_radius and cgto_xyzr.w are set from md_gto_compute_radius_of_influence(cutoff).
// Pass cutoff <= 0.0 to disable culling (sets all radii to FLT_MAX).
// All output arrays must be pre-allocated to num_cgtos / num_pgtos entries respectively.
static void gto_expand_basis(
    float* out_cgto_xyz, float* out_cgto_r, uint32_t* out_cgto_off_len, PGTO* out_pgto,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride, double cutoff)
{
    const size_t stride = atom_xyz_stride == 0 ? sizeof(float) * 3 : atom_xyz_stride;
    uint32_t ci = 0, pi = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        int l      = (int)shell->l;
        int ncart  = (int)md_gto_num_cart_ao(shell->l);
        int nprims = (int)shell->num_primitives;
        uint32_t   prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);
        const float* ap = (const float*)((const uint8_t*)atom_xyz + shell->atom_idx * stride);
        float ax = ap[0];
        float ay = ap[1];
        float az = ap[2];

        for (int ic = 0; ic < ncart; ic++) {
            const int    lx  = lmn[ic][0];
            const int    ly  = lmn[ic][1];
            const int    lz  = lmn[ic][2];
            const double nrm = gto_cart_norm(l, ic);

            uint32_t pi_beg = pi;

            double max_r = 0.0;
            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                float coeff_val = (float)(coef1 * nrm);
                double radius = md_gto_compute_radius_of_influence(lx, ly, lz, (double)coeff_val, (double)alpha, cutoff);
                max_r = MAX(max_r, radius);

                PGTO pgto = {
                    .coeff = coeff_val,
                    .alpha = alpha,
                    .ijkl = md_gto_pack_ijkl(lx, ly, lz, l),
                };
                out_pgto[pi++] = pgto;
            }

            out_cgto_xyz[ci * 3 + 0] = ax;
            out_cgto_xyz[ci * 3 + 1] = ay;
            out_cgto_xyz[ci * 3 + 2] = az;
            out_cgto_r[ci] = (float)max_r;
            out_cgto_off_len[2 * ci + 0] = pi_beg;
            out_cgto_off_len[2 * ci + 1] = pi - pi_beg;
            ci++;
        }
    }
}

static size_t density_matrix_upper_tri_size(size_t n) {
    return n * (n + 1) / 2;
}

// Convert a full N×N row-major double density matrix to compact upper-triangular float.
static void density_matrix_upper_tri_extract_float(float* out, const double* dm, size_t n) {
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            out[k++] = (float)dm[i * n + j];
        }
    }
}

// Extract the upper-triangular sub-block of a full N_full×N_full double density matrix
// for a subset of row/column indices given by old_indices[0..num_kept-1].
// D is symmetric so dm[gi * N_full + gj] is valid for any gi, gj pair.
static void density_matrix_sub_block_upper_tri_float(float* out, const double* dm,
    size_t N_full, const uint32_t* old_indices, uint32_t num_kept) {
    size_t k = 0;
    for (uint32_t i = 0; i < num_kept; i++) {
        uint32_t gi = old_indices[i];
        for (uint32_t j = i; j < num_kept; j++) {
            uint32_t gj = old_indices[j];
            out[k++] = (float)dm[(size_t)gi * N_full + gj];
        }
    }
}

// ---------------------------------------------------------------------------
// Sparse CGTO-pair evaluation: CPU-side pair list construction
// ---------------------------------------------------------------------------
//
// The density can be written as:
//   ρ(r) = Σ_{μ,ν} D_{μν} φ_μ(r) φ_ν(r)
//        = Σ_{μ≤ν} D̃_{μν} φ_μ(r) φ_ν(r)    where D̃_{μν} = D_{μμ} (diagonal)
//                                                             2·D_{μν} (off-diagonal)
//
// Rather than storing the full upper-triangular D matrix and iterating all
// pairs per workgroup, we build a sparse pair list that:
//   1. Morton-sorts CGTOs by center so spatially nearby CGTOs are contiguous.
//   2. Groups them into fixed-size batches of GTO_SPARSE_BATCH_SIZE.
//   3. Retains only pairs (μ,ν) within the same batch where:
//        |D_{μν}| >= dm_threshold   AND   dist(R_μ,R_ν) < r_μ + r_ν
//
// The GPU evaluator holds GTO_SPARSE_BATCH_SIZE phi values in its register
// file, evaluates only the CGTOs active for its 8×8×8 spatial region, and
// plows through the pair list without touching shared memory for phi.
//
// Cross-batch pairs are dropped. For spatially-sorted, localised systems
// these are exactly the long-range pairs with small |D_{μν}|, so the error
// is bounded by dm_threshold. For fully delocalised systems (conjugated,
// charge-transfer) this path degrades gracefully to a larger pair count.

// Number of CGTOs per batch = number of φ registers per GPU thread.
// Adjusting this trades register pressure against batch count.
#define GTO_SPARSE_BATCH_SIZE 64

// A single retained CGTO pair within a batch.
// i, j are LOCAL indices within the batch (0..GTO_SPARSE_BATCH_SIZE-1), j >= i.
// D is pre-doubled for off-diagonal pairs to account for both D_{μν} and D_{νμ}.
typedef struct {
    uint16_t i;
    uint16_t j;
    float    D;
} gto_cgto_pair_t;

// Descriptor for one batch of up to GTO_SPARSE_BATCH_SIZE CGTOs.
typedef struct {
    uint32_t cgto_start;   // First CGTO index in the Morton-sorted ordering
    uint32_t cgto_count;   // Number of CGTOs in this batch (≤ GTO_SPARSE_BATCH_SIZE)
    uint32_t pair_offset;  // Offset into the flat gto_cgto_pair_t array
    uint32_t pair_count;   // Number of pairs in this batch
    float    bbox_xyz[3];  // Bounding sphere center (centroid of member CGTO centers)
    float    bbox_r;       // Bounding sphere radius enclosing all member CGTO spheres
} gto_batch_t;

typedef struct {
    md_allocator_i* alloc;
    size_t num_batches;
    gto_batch_t* batches;          // [num_batches]
    gto_cgto_pair_t* pairs;        // [Σ batch.pair_count] all pairs for all batches
    size_t num_pairs;              // total number of pairs across all batches
    uint32_t* cgto_order;          // [num_cgtos] Morton permutation: cgto_order[new_idx] = original_idx
} gto_sparse_pair_list_t;

// Build the sparse CGTO-pair batch list.
//
// cgto_xyzr    - expanded CGTO data in original order, xyzr = center xyz + cutoff radius
// num_cgtos    - number of CGTOs
// density_matrix - full N×N row-major double matrix (symmetric)
// dm_threshold - pairs with |D_{μν}| < dm_threshold are dropped; pass 0.0 to keep all
// cgto_order   - caller-allocated uint32_t[num_cgtos]; filled with Morton permutation
//                cgto_order[new_idx] = original_idx
//
// out_batches and out_pairs point into temp memory allocated inside this function.
// The caller must copy them (or use them) before calling md_temp_set_pos_back.
//
// Returns the number of batches.
static void gto_build_sparse_pairs(
    gto_sparse_pair_list_t* list,
    const vec4_t*     cgto_xyzr,      // [num_cgtos] original order
    size_t            num_cgtos,
    const double*     density_matrix, // [num_cgtos × num_cgtos] row-major
    double            dm_threshold)
{
    ASSERT(list);
    ASSERT(list->alloc);
    ASSERT(cgto_xyzr);
    ASSERT(density_matrix);

    if (num_cgtos == 0) {
        return;
    }

    // Reset data
    //md_array_resize(list->cgto_order, num_cgtos, list->alloc);
    //md_array_resize(list->batches, num_batches, list->alloc);
    md_array_shrink(list->pairs, 0);
    list->num_batches = 0;
    list->num_pairs = 0;

    // Step 1: Morton-sort CGTO centers.
    // cgto_order[new_idx] = original_idx after the sort.
    //md_util_sort_spatial_xyz(list->cgto_order, (const float*)cgto_xyzr, sizeof(vec4_t), num_cgtos);

    //MEMSET(list->batches, 0, sizeof(gto_batch_t) * num_batches);

    for (size_t i = 0; i < num_cgtos; i++) {
        const float* vi = (const float*)&cgto_xyzr[i];
        float ix = vi[0], iy = vi[1], iz = vi[2], ir = vi[3];

        for (size_t j = i; j < num_cgtos; j++) {
            const float* vj = (const float*)&cgto_xyzr[j];
            float jx = vj[0], jy = vj[1], jz = vj[2], jr = vj[3];

            // Sphere intersection: drop pairs whose CGTO support never overlaps.
            float dx    = ix - jx, dy = iy - jy, dz = iz - jz;
            float dist2 = dx*dx + dy*dy + dz*dz;
            float rsum  = ir + jr;
            if (dist2 > rsum * rsum) continue;

            // Density matrix threshold.
            // Use original row-major layout; matrix is symmetric so [i*N+j] == [j*N+i].
            double dval = density_matrix[(size_t)i * num_cgtos + j];
            if (dm_threshold > 0.0 && fabs(dval) < dm_threshold) continue;

            gto_cgto_pair_t pair = {
                .i = (uint16_t)i,
                .j = (uint16_t)j,
                // Off-diagonal pairs: pre-double to fold in the symmetric D_{νμ} term.
                .D = (i == j) ? (float)dval : (float)(2.0 * dval),
            };
            md_array_push(list->pairs, pair, list->alloc);
        }
    }

    list->num_pairs = md_array_size(list->pairs);
}

size_t md_gto_pgto_count(const md_gto_basis_t* basis) {
    uint32_t nc, np;
    gto_basis_count(&nc, &np, basis);
    return np;
}

size_t md_gto_expand_with_ao_coeffs(md_gto_t* out, const md_gto_basis_t* basis,
    const float* atom_xyz, size_t atom_xyz_stride, const double* ao_coeffs, double cutoff)
{
    ASSERT(out);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(ao_coeffs);

    const size_t stride = atom_xyz_stride == 0 ? sizeof(float) * 3 : atom_xyz_stride;

    size_t num_gtos = 0;
    uint32_t cgto_idx = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        int l      = (int)shell->l;
        int ncart  = (int)md_gto_num_cart_ao(shell->l);
        int nprims = (int)shell->num_primitives;
        uint32_t prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);
        const float* ap = (const float*)((const uint8_t*)atom_xyz + shell->atom_idx * stride);
        float ax = ap[0];
        float ay = ap[1];
        float az = ap[2];

        for (int ic = 0; ic < ncart; ic++) {
            const double ao_coeff = ao_coeffs[cgto_idx++];
            const double nrm      = gto_cart_norm(l, ic);

            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                out[num_gtos++] = (md_gto_t){
                    .x      = ax,
                    .y      = ay,
                    .z      = az,
                    .coeff  = (float)(coef1 * nrm * ao_coeff),
                    .alpha  = alpha,
                    .cutoff = FLT_MAX,
                    .i      = (uint8_t)lmn[ic][0],
                    .j      = (uint8_t)lmn[ic][1],
                    .k      = (uint8_t)lmn[ic][2],
                    .l      = (uint8_t)l,
                };
            }
        }
    }

    if (cutoff > 0.0) {
        num_gtos = md_gto_cutoff_compute_and_filter(out, num_gtos, cutoff);
    }
    return num_gtos;
}

#if !MD_PLATFORM_OSX

#include <core/md_gl_util.h>
#include <gto_shaders.inl>
#include <GL/gl3w.h>

// This should be kept in sync with the define present in segment_and_attribute_to_group.comp
#define QUANTIZATION_SCALE_FACTOR 1.0e6

static GLuint get_gto_program(void) {
    static GLuint program = 0;
    if (!program) {
        GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
        if (md_gl_shader_compile(shader, (str_t){(const char*)eval_gto_comp, eval_gto_comp_size}, 0, 0)) {
            GLuint prog = glCreateProgram();
            if (md_gl_program_attach_and_link(prog, &shader, 1)) {
                program = prog;
            }
        }
        glDeleteShader(shader);
    }
    return program;
}

static GLuint get_gto_density_program(void) {
    static GLuint program = 0;
    if (!program) {
        GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
        if (md_gl_shader_compile(shader, (str_t){(const char*)eval_gto_density_comp, eval_gto_density_comp_size}, 0, 0)) {
            GLuint prog = glCreateProgram();
            if (md_gl_program_attach_and_link(prog, &shader, 1)) {
                program = prog;
            }
        }
        glDeleteShader(shader);
    }
    return program;
}

static GLuint get_gto_density_grad_program(void) {
    static GLuint program = 0;
    if (!program) {
        GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
        if (md_gl_shader_compile(shader, (str_t) { (const char*)eval_gto_density_grad_comp, eval_gto_density_grad_comp_size }, 0, 0)) {
            GLuint prog = glCreateProgram();
            if (md_gl_program_attach_and_link(prog, &shader, 1)) {
                program = prog;
            }
        }
        glDeleteShader(shader);
    }
    return program;
}


static GLuint get_buffer(size_t size) {
    GLuint id = 0;
    glCreateBuffers(1, &id);
    glBindBuffer(GL_ARRAY_BUFFER, id);
    glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return id;
}

static void free_buffer(GLuint id) {
    if (glIsBuffer(id)) {
        glDeleteBuffers(1, &id);
    }
}

static void gto_grid_evaluate_mo_GPU(uint32_t vol_tex, const md_grid_t* grid, md_gto_t* gtos, uint32_t* orb_offsets, float* orb_scaling, size_t num_orbs, md_gto_eval_mode_t mode, md_gto_op_t op, GLuint program) {
    ASSERT(grid);
    ASSERT(gtos);
    ASSERT(orb_offsets);
    ASSERT(orb_scaling);

    if (num_orbs == 0) {
        return;
    }

    md_gl_debug_push("EVAL ORBS");

    if (!glIsTexture(vol_tex)) {
        MD_LOG_ERROR("Invalid volume texture handle");
        return;
    }

    GLenum format = 0;
    if (glGetTextureLevelParameteriv) {
        glGetTextureLevelParameteriv(vol_tex,   0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);
    } else {
        glBindTexture(GL_TEXTURE_3D, vol_tex);
        glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, (GLint*)&format);
        glBindTexture(GL_TEXTURE_3D, 0);
    }

    switch (format) {
    case GL_R16F:
    case GL_R32F:
        break;
    default:
        // Not good
        MD_LOG_ERROR("Unrecognized internal format of supplied volume texture");
        goto done;
    }

    size_t num_gtos = orb_offsets[num_orbs];

    GLintptr   ssbo_gto_offset = 0;
    GLsizeiptr ssbo_gto_size   = sizeof(md_gto_t) * num_gtos;

    GLintptr   ssbo_orb_offset = ALIGN_TO(ssbo_gto_offset + ssbo_gto_size, 256);
    GLsizeiptr ssbo_orb_size   = sizeof(uint32_t) * (num_orbs + 1);

    GLintptr   ssbo_scl_offset = ALIGN_TO(ssbo_orb_offset + ssbo_orb_size, 256);
    GLsizeiptr ssbo_scl_size   = sizeof(float) * (num_orbs);

    size_t total_size = ALIGN_TO(ssbo_scl_offset + ssbo_scl_size, 256);
    GLuint ssbo = get_buffer(total_size);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_gto_offset, ssbo_gto_size, gtos);

    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_orb_offset, ssbo_orb_size, orb_offsets);
    // Fill last portion of buffer with point indices
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_scl_offset, ssbo_scl_size, orb_scaling);

    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, ssbo, ssbo_gto_offset, ssbo_gto_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 1, ssbo, ssbo_orb_offset, ssbo_orb_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 2, ssbo, ssbo_scl_offset, ssbo_scl_size);

    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);

    glUseProgram(program);

    float world_to_model[4][4];
    float index_to_world[4][4];
    const float sample_offset[3] = { 0.5f, 0.5f, 0.5f }; // Sample at voxel centers

    world_to_model_matrix(world_to_model, grid);
    index_to_world_matrix(index_to_world, grid, sample_offset);

    glUniformMatrix4fv(0, 1, GL_FALSE, (const float*)world_to_model);
    glUniformMatrix4fv(1, 1, GL_FALSE, (const float*)index_to_world);
    glUniform3fv(2, 1, grid->spacing.elem);
    glUniform1ui(3, (GLuint)num_orbs);
    glUniform1i(4, (GLint)mode);
    glUniform1ui(5, (GLuint)op);

    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_WRITE, format);

    int num_groups[3] = {
        DIV_UP(grid->dim[0], 8),
        DIV_UP(grid->dim[1], 8),
        DIV_UP(grid->dim[2], 8),
    };

    glDispatchCompute(num_groups[0], num_groups[1], num_groups[2]);

    glUseProgram(0);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT | GL_PIXEL_BUFFER_BARRIER_BIT);

    free_buffer(ssbo);
done:
    md_gl_debug_pop();
}

void md_gto_grid_evaluate_matrix_GPU(uint32_t vol_tex, const md_grid_t* grid,
    uint32_t num_cgtos, const vec4_t* cgto_xyzr, const uint32_t* cgto_off_len,
    uint32_t num_pgtos, const PGTO* pgto,
    const float* upper_triangular_matrix, size_t upper_triangular_len,
    bool include_gradients, md_gto_op_t op) {
    ASSERT(grid);
    ASSERT(upper_triangular_matrix);

    md_gl_debug_push("EVAL DENSITY");

    if (!glIsTexture(vol_tex)) {
        MD_LOG_ERROR("Invalid volume texture handle");
        return;
    }

    GLint format = 0;
    glBindTexture(GL_TEXTURE_3D, vol_tex);
    glGetTexLevelParameteriv(GL_TEXTURE_3D, 0, GL_TEXTURE_INTERNAL_FORMAT, &format);
    glBindTexture(GL_TEXTURE_3D, 0);

    GLuint program = 0;
    if (include_gradients) {
        switch (format) {
        case GL_RGBA16F:
        case GL_RGBA32F:
            break;
        default:
            // Not good
            MD_LOG_ERROR("Unrecognized internal format of supplied volume texture");
            goto done;
        }
        program = get_gto_density_grad_program();
    }
    else {
        switch (format) {
        case GL_R16F:
        case GL_R32F:
            break;
        default:
            // Not good
            MD_LOG_ERROR("Unrecognized internal format of supplied volume texture");
            goto done;
        }
        program = get_gto_density_program();
    }

    if (!program) {
        MD_LOG_ERROR("Program not found?!");
        goto done;
    }

    size_t matrix_dim = num_cgtos;

    typedef struct {
        mat4_t world_to_model;
        mat4_t index_to_world;
        vec4_t step;
        uint32_t D_matrix_dim;
        uint32_t operation;
        uint32_t _pad[2];
    } uniform_block_t;

    uniform_block_t ub_data = {0};
    world_to_model_matrix(ub_data.world_to_model.elem, grid);
    const float sample_offset[3] = { 0.5f, 0.5f, 0.5f }; // Sample at voxel centers
    index_to_world_matrix(ub_data.index_to_world.elem, grid, sample_offset);
    ub_data.step = vec4_from_vec3(grid->spacing, 0);
    ub_data.D_matrix_dim = (uint32_t)matrix_dim;
    ub_data.operation = (uint32_t)op;

    GLintptr   ssbo_cgto_xyzr_base      = 0;
    GLsizeiptr ssbo_cgto_xyzr_size      = sizeof(vec4_t) * num_cgtos;

    GLintptr   ssbo_cgto_off_len_base   = ALIGN_TO(ssbo_cgto_xyzr_base + ssbo_cgto_xyzr_size, 256);
    GLsizeiptr ssbo_cgto_off_len_size   = sizeof(uint32_t) * num_cgtos * 2;

    GLintptr   ssbo_pgto_base           = ALIGN_TO(ssbo_cgto_off_len_base + ssbo_cgto_off_len_size, 256);
    GLsizeiptr ssbo_pgto_size           = sizeof(PGTO) * num_pgtos;

    GLintptr   ssbo_matrix_base         = ALIGN_TO(ssbo_pgto_base + ssbo_pgto_size, 256);
    GLsizeiptr ssbo_matrix_size         = sizeof(float) * upper_triangular_len;

    GLintptr   ubo_base                 = ALIGN_TO(ssbo_matrix_base + ssbo_matrix_size, 256);
    GLsizeiptr ubo_size                 = sizeof(uniform_block_t);

    size_t total_size = ALIGN_TO(ubo_base + ubo_size, 256);
    GLuint buf = get_buffer(total_size);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buf);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_cgto_xyzr_base,    ssbo_cgto_xyzr_size,    cgto_xyzr);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_cgto_off_len_base, ssbo_cgto_off_len_size, cgto_off_len);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_pgto_base,         ssbo_pgto_size,         pgto);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, ssbo_matrix_base,       ssbo_matrix_size,       upper_triangular_matrix);

    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, buf, ssbo_cgto_xyzr_base,    ssbo_cgto_xyzr_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 1, buf, ssbo_cgto_off_len_base, ssbo_cgto_off_len_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 2, buf, ssbo_pgto_base,         ssbo_pgto_size);
    glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 3, buf, ssbo_matrix_base,       ssbo_matrix_size);

    glBindBuffer(GL_UNIFORM_BUFFER, buf);
    glBufferSubData(GL_UNIFORM_BUFFER, ubo_base, ubo_size, &ub_data);
    glBindBufferRange(GL_UNIFORM_BUFFER, 0, buf, ubo_base, ubo_size);

    glUseProgram(program);
    GLuint block = glGetUniformBlockIndex(program, "UniformBlock");
    glUniformBlockBinding(program, block, 0);

    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_WRITE, format);

    GLuint query;
    glGenQueries(1, &query);

    // Start timing
    glBeginQuery(GL_TIME_ELAPSED, query);

    int num_groups[3] = {
        DIV_UP(grid->dim[0], 8),
        DIV_UP(grid->dim[1], 8),
        DIV_UP(grid->dim[2], 8),
    };
    glDispatchCompute(num_groups[0], num_groups[1], num_groups[2]);

    // End timing
    glEndQuery(GL_TIME_ELAPSED);

    // Retrieve the result (blocking until GPU finishes)
    GLuint64 elapsedTime = 0;
    glGetQueryObjectui64v(query, GL_QUERY_RESULT, &elapsedTime); // nanoseconds

	MD_LOG_DEBUG("GTO Density evaluation of [%i,%i,%i] GPU time: %.3f ms", grid->dim[0], grid->dim[1], grid->dim[2], elapsedTime / 1e6);

    glUseProgram(0);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT | GL_SHADER_IMAGE_ACCESS_BARRIER_BIT | GL_PIXEL_BUFFER_BARRIER_BIT);

done:
    md_gl_debug_pop();
}

void md_gto_grid_evaluate_mo_GL(uint32_t vol_tex, const md_grid_t* grid, const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride, const double* mo_coeffs, double cutoff, md_gto_eval_mode_t mode, md_gto_op_t op) {
    ASSERT(grid);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(mo_coeffs);

    GLuint program = get_gto_program();
    if (!program) return;

    md_temp_scope_t temp = md_temp_begin();
    size_t max_gtos  = md_gto_pgto_count(basis);
    md_gto_t* gtos   = (md_gto_t*)md_temp_alloc(temp, sizeof(md_gto_t) * max_gtos);
    size_t num_gtos  = md_gto_expand_with_ao_coeffs(gtos, basis, atom_xyz, atom_xyz_stride, mo_coeffs, cutoff);

    if (num_gtos > 0) {
        uint32_t orb_offsets[2] = { 0, (uint32_t)num_gtos };
        float    orb_scaling[1] = { 1.0f };
        gto_grid_evaluate_mo_GPU(vol_tex, grid, gtos, orb_offsets, orb_scaling, 1, mode, op, program);
    }

    md_temp_end(temp);
}

void md_gto_grid_evaluate_multi_mo_GL(uint32_t vol_tex, const md_grid_t* grid, const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride, const double* mo_coeffs[], const double mo_scl[], size_t num_mos, double cutoff, md_gto_eval_mode_t mode, md_gto_op_t op) {
    ASSERT(grid);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(mo_coeffs);

    if (num_mos == 0) return;

    GLuint program = get_gto_program();
    if (!program) return;

    md_temp_scope_t temp = md_temp_begin();
    size_t max_gtos  = md_gto_pgto_count(basis);

    // Flat GTO buffer for all MOs concatenated
    md_gto_t* gtos        = (md_gto_t*)md_temp_alloc(temp, sizeof(md_gto_t) * max_gtos * num_mos);
    uint32_t* orb_offsets = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * (num_mos + 1));
    float*    orb_scaling = (float*)   md_temp_alloc(temp, sizeof(float)    *  num_mos);

    size_t total_gtos = 0;
    size_t num_orbs   = 0;
    orb_offsets[0]    = 0;

    for (size_t i = 0; i < num_mos; i++) {
        if (!mo_coeffs[i]) continue;
        size_t n = md_gto_expand_with_ao_coeffs(gtos + total_gtos, basis, atom_xyz, atom_xyz_stride, mo_coeffs[i], cutoff);
        if (n == 0) continue;
        total_gtos += n;
        orb_scaling[num_orbs] = (float)(mo_scl ? mo_scl[i] : 1.0);
        orb_offsets[++num_orbs] = (uint32_t)total_gtos;
    }

    if (total_gtos > 0) {
        gto_grid_evaluate_mo_GPU(vol_tex, grid, gtos, orb_offsets, orb_scaling, num_orbs, mode, op, program);
    }

    md_temp_end(temp);
}

void md_gto_grid_evaluate_density_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
    const double* density_matrix, bool include_gradients, md_gto_op_t op)
{
    ASSERT(basis);
    ASSERT(density_matrix);

    uint32_t num_cgtos, num_pgtos;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);

    md_temp_scope_t temp = md_temp_begin();
    float* upper_tri = (float*)md_temp_alloc(temp, sizeof(float) * density_matrix_upper_tri_size(num_cgtos));
    density_matrix_upper_tri_extract_float(upper_tri, density_matrix, num_cgtos);
    md_gto_grid_evaluate_density_packed_GL(vol_tex, grid, basis, atom_xyz, atom_xyz_stride, upper_tri, include_gradients, op);
    md_temp_end(temp);
}

void md_gto_grid_evaluate_density_packed_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
    const float* upper_tri, bool include_gradients, md_gto_op_t op)
{
    ASSERT(grid);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(upper_tri);

    uint32_t num_cgtos, num_pgtos;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);

    md_temp_scope_t temp = md_temp_begin();
    float*    cgto_xyz     = (float*)   md_temp_alloc(temp, sizeof(float)    * 3 * num_cgtos);
    float*    cgto_r       = (float*)   md_temp_alloc(temp, sizeof(float)    * 1 * num_cgtos);
    uint32_t* cgto_off_len = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_cgtos * 2);
    PGTO*     pgto         = (PGTO*)    md_temp_alloc(temp, sizeof(PGTO)     * num_pgtos);
    size_t    tri_len      = density_matrix_upper_tri_size(num_cgtos);

    gto_expand_basis(cgto_xyz, cgto_r, cgto_off_len, pgto, basis, atom_xyz, atom_xyz_stride, 1.0e-6);


    // Recombine into float4 for the GL path, which keeps its own xyzr SSBO layout.
    vec4_t* cgto_xyzr = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * num_cgtos);
    for (uint32_t i = 0; i < num_cgtos; ++i) {
        cgto_xyzr[i] = (vec4_t){cgto_xyz[i*3+0], cgto_xyz[i*3+1], cgto_xyz[i*3+2], cgto_r[i]};
    }

    md_gto_grid_evaluate_matrix_GPU(vol_tex, grid, num_cgtos, cgto_xyzr, cgto_off_len, num_pgtos, pgto, upper_tri, tri_len, include_gradients, op);

    md_temp_end(temp);
}

#else


// GPU-accelerated versions of the above evaluation functions.  See md_gto.c for details on the expected data layout and GPU buffer formats.
void md_gto_grid_evaluate_mo_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
    const double* mo_coeffs, double cutoff, md_gto_eval_mode_t mode, md_gto_op_t op) {
    (void)vol_tex; (void)grid; (void)basis; (void)atom_xyz; (void)atom_xyz_stride; (void)mo_coeffs; (void)cutoff; (void)mode; (void)op;
}

// mo_scl is optional and if null is supplied, then it is assumed that all orbitals should be scaled by 1.0 (i.e. no relative scaling between orbitals).
void md_gto_grid_evaluate_multi_mo_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
    const double* mo_coeffs[], const double mo_scl[], size_t num_mos, double cutoff, md_gto_eval_mode_t mode, md_gto_op_t op) {
    (void)vol_tex; (void)grid; (void)basis; (void)atom_xyz; (void)atom_xyz_stride; (void)mo_coeffs; (void)mo_scl; (void)num_mos; (void)cutoff; (void)mode; (void)op;
}

void md_gto_grid_evaluate_density_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
    const double* density_matrix, bool include_gradients, md_gto_op_t op) {
    (void)vol_tex; (void)grid; (void)basis; (void)atom_xyz; (void)atom_xyz_stride; (void)density_matrix; (void)include_gradients; (void)op;
}

void md_gto_grid_evaluate_density_packed_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
    const float* upper_tri, bool include_gradients, md_gto_op_t op) {
    (void)vol_tex; (void)grid; (void)basis; (void)atom_xyz; (void)atom_xyz_stride; (void)upper_tri; (void)include_gradients; (void)op;
}

#endif

#if MD_ENABLE_GPU

#include <core/md_gpu.h>
#include <gto_gpu_shaders.inl>


static md_gpu_kernel_t gto_k_density = NULL;
static md_gpu_kernel_t gto_k_mo      = NULL;

// Kernels other than the two reference ones are created on first use, so that only
// the path actually taken pays for its (possibly runtime) shader compilation.
//
// The set was chosen from measurements on NVIDIA GTX 1060 and L4, Apple M4 Pro and
// Intel HD 530 (branch gto_kernel_eval, results/): density through the two-pass GEMM
// with 128 voxels per group and 32- or 64-wide AO tiles; orbitals through the shell
// kernel with 4 voxels per thread (up to 8 orbitals) or the GEMM path (more).
static md_gpu_device_t gto_gpu_device = NULL;
static uint32_t        gto_gpu_vendor = 0;
static md_gpu_kernel_t gto_k_density_gemm_count   = NULL;
static md_gpu_kernel_t gto_k_density_gemm_phi     = NULL;
static md_gpu_kernel_t gto_k_density_gemm_dgather = NULL;
static md_gpu_kernel_t gto_k_density_gemm2[2] = {0};   // [32, 64 AO tile], 128 voxels per group
static md_gpu_kernel_t gto_k_mo_shell[2]      = {0};   // [1, 8 orbitals per pass], 4 voxels per thread
static md_gpu_kernel_t gto_k_mo_gemm          = NULL;  // 32 orbitals per tile
static md_gpu_kernel_t gto_k_dmat_tiles       = NULL;  // md_gto_gpu_dmat_build
// Held while a kernel is made: the evaluations and the dmat may be driven from different threads, each on a
// stream of its own, and both make kernels on first use
static md_mutex_t      gto_kernel_mutex;
static bool            gto_kernel_mutex_ok = false;

static void ensure_kernel(md_gpu_device_t device, md_gpu_kernel_t* slot, md_gpu_kernel_desc_t desc) {
    if (gto_kernel_mutex_ok) md_mutex_lock(&gto_kernel_mutex);
    if (!*slot) {
        *slot = md_gpu_kernel_create(device, &desc);
        if (!*slot) MD_LOG_ERROR("md_gto: failed to create kernel '%s': %s", desc.label, md_gpu_last_error());
    }
    if (gto_kernel_mutex_ok) md_mutex_unlock(&gto_kernel_mutex);
}

static void destroy_kernel(md_gpu_kernel_t* slot) {
    if (*slot) { md_gpu_kernel_destroy(*slot); *slot = NULL; }
}

void md_gto_gpu_initialize(md_gpu_device_t device) {
    if (!device) return;
    gto_gpu_device = device;
    if (!gto_kernel_mutex_ok) gto_kernel_mutex_ok = md_mutex_init(&gto_kernel_mutex);
    md_gpu_device_info_t info;
    gto_gpu_vendor = md_gpu_device_info(device, &info) ? info.vendor_id : 0;
    ensure_kernel(device, &gto_k_density, md_shader_eval_gto_density_main_kernel());
    ensure_kernel(device, &gto_k_mo,      md_shader_eval_gto_mo_main_kernel());
}

void md_gto_gpu_shutdown(void) {
    destroy_kernel(&gto_k_density);
    destroy_kernel(&gto_k_mo);
    destroy_kernel(&gto_k_density_gemm_count);
    destroy_kernel(&gto_k_density_gemm_phi);
    destroy_kernel(&gto_k_density_gemm_dgather);
    for (int i = 0; i < 2; ++i) {
        destroy_kernel(&gto_k_density_gemm2[i]);
        destroy_kernel(&gto_k_mo_shell[i]);
    }
    destroy_kernel(&gto_k_mo_gemm);
    destroy_kernel(&gto_k_dmat_tiles);
    gto_gpu_device = NULL;
}

static bool gto_density_gemm_kernels(void) {
    if (!gto_gpu_device) return false;
    ensure_kernel(gto_gpu_device, &gto_k_density_gemm_count,   md_shader_eval_gto_density_gemm_count_main_kernel());
    ensure_kernel(gto_gpu_device, &gto_k_density_gemm_phi,     md_shader_eval_gto_density_gemm_phi_main_kernel());
    ensure_kernel(gto_gpu_device, &gto_k_density_gemm_dgather, md_shader_eval_gto_density_gemm_dgather_main_kernel());
    return gto_k_density_gemm_count && gto_k_density_gemm_phi && gto_k_density_gemm_dgather;
}

// The GEMM pass, 128 voxels per group, with a 32- or 64-wide AO tile.
static md_gpu_kernel_t gto_density_gemm2_kernel(uint32_t gm) {
    const int im = gm == 32 ? 0 : gm == 64 ? 1 : -1;
    if (im < 0 || !gto_gpu_device) return NULL;
    md_gpu_kernel_t* slot = &gto_k_density_gemm2[im];
    if (!*slot) {
        ensure_kernel(gto_gpu_device, slot, im == 0 ? md_shader_eval_gto_density_gemm2_128x32_main_kernel()
                                                    : md_shader_eval_gto_density_gemm2_128x64_main_kernel());
    }
    return *slot;
}

/* Argument structs, mirroring the kernels in src/shaders/gto/. Both lead with two
   float4x4 at offset 0, where SPIR-V and MSL agree.
   tools/check_gpu_arg_layout.py verifies these against the compiled shaders. */
typedef struct {
    md_gpu_float4x4 world_to_model;
    md_gpu_float4x4 index_to_world;
    md_gpu_float4   step;
    md_gpu_uint4    grid_dim;
    uint32_t        num_cgtos;
    uint32_t        operation;
    uint32_t        dm_tiled;
    uint32_t        dm_nt;
    md_gpu_addr_t   cgto_atom_idx;
    md_gpu_addr_t   cgto_r;
    md_gpu_addr_t   cgto_off_len;
    md_gpu_addr_t   pgto;
    md_gpu_addr_t   atom_xyz;
    md_gpu_addr_t   D_matrix;
    md_gpu_addr_t   dm_tao;
    md_gpu_addr_t   dm_map;
    md_gpu_storage_tex_t out_tex;
} gto_density_args_t;

/* Mirrors RootArgs in eval_gto_density_gemm.slang (shared by its three entry points). */
typedef struct {
    md_gpu_float4x4 world_to_model;
    md_gpu_float4x4 index_to_world;
    md_gpu_float4   step;
    md_gpu_uint4    grid_dim;
    md_gpu_uint4    block_dim;
    uint32_t        num_shells;
    uint32_t        num_aos;
    uint32_t        operation;
    uint32_t        num_blocks;
    uint32_t        d_mode;
    uint32_t        num_mos;
    uint32_t        mo_mode;
    uint32_t        has_weights;
    uint32_t        dm_tiled;
    uint32_t        dm_nt;
    uint32_t        _dm_pad0;
    uint32_t        _dm_pad1;
    md_gpu_addr_t   shell_atom_idx;
    md_gpu_addr_t   shell_r;
    md_gpu_addr_t   shell_info;
    md_gpu_addr_t   shell_prim;
    md_gpu_addr_t   atom_xyz;
    md_gpu_addr_t   D_matrix;
    md_gpu_addr_t   block_counts;
    md_gpu_addr_t   block_table;
    md_gpu_addr_t   ao_list;
    md_gpu_addr_t   phi;
    md_gpu_addr_t   D_blocks;
    md_gpu_addr_t   mo_weights;
    md_gpu_addr_t   dm_tao;
    md_gpu_addr_t   dm_map;
    md_gpu_storage_tex_t out_tex;
} gto_density_gemm_args_t;

/* Mirrors DmatArgs in gto_dmat.slang. */
typedef struct {
    uint32_t      num_mo;
    uint32_t      num_items;
    uint32_t      items_x;
    uint32_t      ntao;
    md_gpu_addr_t P;
    md_gpu_addr_t w;
    md_gpu_addr_t items;
    md_gpu_addr_t cols;
    md_gpu_addr_t tiles;
} gto_dmat_args_t;

typedef struct {
    md_gpu_float4x4 world_to_model;
    md_gpu_float4x4 index_to_world;
    md_gpu_float4   step;
    md_gpu_uint4    grid_dim;
    uint32_t        num_cgtos;
    uint32_t        num_rows;
    uint32_t        mode;
    uint32_t        operation;
    md_gpu_addr_t   cgto_atom_idx;
    md_gpu_addr_t   cgto_r;
    md_gpu_addr_t   cgto_off_len;
    md_gpu_addr_t   pgto;
    md_gpu_addr_t   atom_xyz;
    md_gpu_addr_t   coeffs;
    md_gpu_storage_tex_t out_tex;
} gto_mo_args_t;

typedef struct {
    uint32_t num_cgtos;
    uint32_t num_pgtos;
    uint32_t num_atoms;
    uint32_t num_shells;
    uint32_t num_shell_prims;
    uint32_t off_cgto_atom_idx; // uint  × num_cgtos
    uint32_t off_cgto_r;        // float × num_cgtos
    uint32_t off_cgto_off_len;  // uint2 × num_cgtos
    uint32_t off_pgto;          // PGTO  × num_pgtos
    // Shell-level data, used by the shell and GEMM kernels.
    uint32_t off_shell_atom_idx; // uint  × num_shells
    uint32_t off_shell_r;        // float × num_shells
    uint32_t off_shell_info;     // uint4 × num_shells: (prim_offset, num_prims, l, ao_offset)
    uint32_t off_shell_prim;     // float2 × num_shell_prims: (coeff, alpha * log2(e))
    uint64_t total_size;
} md_gto_basis_layout_t;

static md_gto_basis_layout_t gto_basis_layout_compute(uint32_t num_cgtos, uint32_t num_pgtos, uint32_t num_atoms,
                                                      uint32_t num_shells, uint32_t num_shell_prims) {
    md_gto_basis_layout_t L = {0};
    L.num_cgtos = num_cgtos;
    L.num_pgtos = num_pgtos;
    L.num_atoms = num_atoms;
    L.num_shells = num_shells;
    L.num_shell_prims = num_shell_prims;

    L.off_cgto_atom_idx = 0;
    uint32_t end_cgto_atom_idx = L.off_cgto_atom_idx + (uint32_t)(sizeof(uint32_t) * num_cgtos);

    L.off_cgto_r       = (uint32_t)ALIGN_TO(end_cgto_atom_idx, 256);
    uint32_t end_cgto_r       = L.off_cgto_r       + (uint32_t)(sizeof(float)    * 1 * num_cgtos);

    L.off_cgto_off_len = (uint32_t)ALIGN_TO(end_cgto_r,       256);
    uint32_t end_cgto_off_len = L.off_cgto_off_len + (uint32_t)(sizeof(uint32_t) * 2 * num_cgtos);

    L.off_pgto         = (uint32_t)ALIGN_TO(end_cgto_off_len, 256);
    uint32_t end_pgto  = L.off_pgto + (uint32_t)(sizeof(PGTO) * num_pgtos);

    L.off_shell_atom_idx = (uint32_t)ALIGN_TO(end_pgto, 256);
    uint32_t end_shell_atom_idx = L.off_shell_atom_idx + (uint32_t)(sizeof(uint32_t) * num_shells);

    L.off_shell_r      = (uint32_t)ALIGN_TO(end_shell_atom_idx, 256);
    uint32_t end_shell_r = L.off_shell_r + (uint32_t)(sizeof(float) * num_shells);

    L.off_shell_info   = (uint32_t)ALIGN_TO(end_shell_r, 256);
    uint32_t end_shell_info = L.off_shell_info + (uint32_t)(sizeof(uint32_t) * 4 * num_shells);

    L.off_shell_prim   = (uint32_t)ALIGN_TO(end_shell_info, 256);
    uint32_t end_shell_prim = L.off_shell_prim + (uint32_t)(sizeof(float) * 2 * num_shell_prims);

    L.total_size       = (uint64_t)ALIGN_TO(end_shell_prim, 256);
    return L;
}

// ---------------------------------------------------------------------------
// md_gto_gpu_basis_t  —  device-local basis buffer (atom_idx, r, off_len, pgto)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// md_gto_gpu_basis_t  -  device-local basis allocation
// ---------------------------------------------------------------------------

typedef struct md_gto_gpu_basis {
    md_gpu_addr_t         buffer;
    md_gto_basis_layout_t layout;
} md_gto_gpu_basis;

md_gto_gpu_basis_t md_gto_gpu_basis_create(md_gpu_stream_t stream, const md_gto_gpu_basis_desc_t* desc) {
    ASSERT(stream);
    ASSERT(desc && desc->basis);
    const md_gto_basis_t* basis = desc->basis;

    uint32_t num_cgtos = 0, num_pgtos = 0;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);
    uint32_t num_atoms = gto_basis_num_atoms(basis);

    md_gto_basis_layout_t layout = gto_basis_layout_compute(num_cgtos, num_pgtos, num_atoms, basis->num_shells, gto_basis_num_shell_prims(basis));

    md_gpu_addr_t buf = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, (size_t)layout.total_size).gpu;
    if (!buf) {
        MD_LOG_ERROR("md_gto_gpu_basis_create: failed to allocate %zu bytes: %s",
                     (size_t)layout.total_size, md_gpu_last_error());
        return NULL;
    }

    md_gto_gpu_basis* gb = (md_gto_gpu_basis*)calloc(1, sizeof(md_gto_gpu_basis));
    if (!gb) { md_gpu_free(stream, buf); return NULL; }
    gb->buffer = buf;
    gb->layout = layout;

    const md_gto_basis_layout_t* L = &gb->layout;
    md_temp_scope_t temp = md_temp_begin();

    uint32_t* cgto_atom_idx = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * 1 * L->num_cgtos);
    float*    cgto_r        = (float*)   md_temp_alloc(temp, sizeof(float)    * 1 * L->num_cgtos);
    uint32_t* cgto_off_len  = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * 2 * L->num_cgtos);
    PGTO*     pgto          = (PGTO*)    md_temp_alloc(temp, sizeof(PGTO)         * L->num_pgtos);

    gto_expand_basis_gpu_meta(cgto_atom_idx, cgto_r, cgto_off_len, pgto, basis, desc->cutoff);

    uint32_t* shell_atom_idx = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * MAX(1, L->num_shells));
    float*    shell_r        = (float*)   md_temp_alloc(temp, sizeof(float)    * MAX(1, L->num_shells));
    uint32_t* shell_info     = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * 4 * MAX(1, L->num_shells));
    float*    shell_prim     = (float*)   md_temp_alloc(temp, sizeof(float)    * 2 * MAX(1, L->num_shell_prims));
    gto_expand_basis_gpu_shells(shell_atom_idx, shell_r, shell_info, shell_prim, basis, desc->cutoff);

    const size_t sz_atom_idx = sizeof(uint32_t) * 1 * L->num_cgtos;
    const size_t sz_r        = sizeof(float)    * 1 * L->num_cgtos;
    const size_t sz_off_len  = sizeof(uint32_t) * 2 * L->num_cgtos;
    const size_t sz_pgto     = sizeof(PGTO)         * L->num_pgtos;

    /* One path for both UMA and discrete: md_gpu hands back either the
       destination itself or staging, whichever is safe. */
    uint8_t* dst = (uint8_t*)md_gpu_upload_begin(stream, buf, (size_t)L->total_size);
    bool success = false;
    if (dst) {
        MEMCPY(dst + L->off_cgto_atom_idx, cgto_atom_idx, sz_atom_idx);
        MEMCPY(dst + L->off_cgto_r,        cgto_r,        sz_r);
        MEMCPY(dst + L->off_cgto_off_len,  cgto_off_len,  sz_off_len);
        MEMCPY(dst + L->off_pgto,          pgto,          sz_pgto);
        MEMCPY(dst + L->off_shell_atom_idx, shell_atom_idx, sizeof(uint32_t) * 1 * L->num_shells);
        MEMCPY(dst + L->off_shell_r,        shell_r,        sizeof(float)    * 1 * L->num_shells);
        MEMCPY(dst + L->off_shell_info,     shell_info,     sizeof(uint32_t) * 4 * L->num_shells);
        MEMCPY(dst + L->off_shell_prim,     shell_prim,     sizeof(float)    * 2 * L->num_shell_prims);
        success = md_gpu_upload_end(stream);
    } else {
        MD_LOG_ERROR("md_gto_gpu_basis_create: upload failed: %s", md_gpu_last_error());
    }

    md_temp_end(temp);
    if (!success) {
        md_gpu_free(stream, buf);
        free(gb);
        return NULL;
    }
    return gb;
}

void md_gto_gpu_basis_destroy(md_gpu_stream_t stream, md_gto_gpu_basis_t gb) {
    if (!gb) return;
    /* Stream-ordered: safe with work that uses the basis still in flight on
       `stream`. */
    md_gpu_free(stream, gb->buffer);
    free(gb);
}

md_gpu_addr_t md_gto_gpu_basis_buffer(md_gto_gpu_basis_t gb) {
    return gb ? gb->buffer : 0;
}

size_t md_gto_gpu_basis_num_cgtos(md_gto_gpu_basis_t gb) {
    if (!gb) return 0;
    return gb->layout.num_cgtos;
}

size_t md_gto_gpu_basis_num_pgtos(md_gto_gpu_basis_t gb) {
    if (!gb) return 0;
    return gb->layout.num_pgtos;
}

size_t md_gto_gpu_basis_num_atoms(md_gto_gpu_basis_t gb) {
    if (!gb) return 0;
    return gb->layout.num_atoms;
}

// ---------------------------------------------------------------------------
// md_gto_gpu_dmat_t  -  a density matrix in 16 x 16 tiles of the AO pairs that meet
// ---------------------------------------------------------------------------

#define GTO_DMAT_T    16u
#define GTO_DMAT_NONE 0xFFFFFFFFu
// The build's revision, reported by md_gto_gpu_dmat_info (the benchmark prints it, so a measurement says
// which kernel it is of). Bump with every change to gto_dmat.slang or to how it is launched.
#define GTO_DMAT_REVISION "r2 (64 threads per 16 x 128 block, 4 x 8 per thread, KC 16, next chunk in registers)"

struct md_gto_gpu_dmat {
    uint32_t      num_ao, num_atoms, nt;    // nt: tiles per edge
    uint64_t      num_tiles, num_items;
    double        margin;
    float*        xyz;                      // the positions the pattern was made for, 3 per atom
    int32_t*      perm;                     // tiled AO -> AO, -1 past num_ao (16 nt)
    uint64_t*     row_ptr;                  // per row tile: its first entry (nt + 1)
    uint32_t*     cols;                     // per entry: its column tile (>= the row tile)
    md_gpu_addr_t tao, map, tiles, items, cols_buf;
    md_gpu_addr_t P, w;                     // the orbitals of md_gto_gpu_dmat_build, grown as needed
    size_t        P_cap, w_cap;             // floats
};

// The shells' sort keys themselves, which hold the shell in their low bits: no state beside the arrays, so that
// patterns may be made on several threads at once
static int gto_dmat_cmp(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return x < y ? -1 : x > y;
}

static void gto_dmat_free_host(md_gto_gpu_dmat_t dm) {
    free(dm->xyz); free(dm->perm); free(dm->row_ptr); free(dm->cols);
}

md_gto_gpu_dmat_t md_gto_gpu_dmat_create(md_gpu_stream_t stream, const md_gto_gpu_dmat_desc_t* desc) {
    if (!stream || !desc || !desc->basis || !desc->atom_xyz || !(desc->cutoff > 0.0) || !(desc->margin >= 0.0)) {
        MD_LOG_ERROR("md_gto_gpu_dmat_create: invalid arguments (a basis, positions, a cutoff > 0, a margin >= 0)");
        return NULL;
    }
    const md_gto_basis_t* basis = desc->basis;
    const uint32_t ns = basis->num_shells;
    const uint32_t num_ao = (uint32_t)md_gto_basis_num_ao(basis);
    const uint32_t num_atoms = gto_basis_num_atoms(basis);
    if (ns == 0 || num_ao == 0) {
        MD_LOG_ERROR("md_gto_gpu_dmat_create: empty basis");
        return NULL;
    }
    if (ns > 0xFFFFFu || num_atoms > 0xFFFFFu) {
        MD_LOG_ERROR("md_gto_gpu_dmat_create: more than 2^20 shells or atoms");
        return NULL;
    }
    md_gto_gpu_dmat_t dm = (md_gto_gpu_dmat_t)calloc(1, sizeof(struct md_gto_gpu_dmat));
    if (!dm) return NULL;
    dm->num_ao = num_ao;
    dm->num_atoms = num_atoms;
    dm->margin = desc->margin;
    dm->nt = (num_ao + GTO_DMAT_T - 1) / GTO_DMAT_T;
    const uint32_t nt = dm->nt;

    md_temp_scope_t temp = md_temp_begin();
    // the shells' screening radii, as md_gto_gpu_basis_create makes them for this cutoff
    uint32_t* shell_atom_idx = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * ns);
    float*    shell_r        = (float*)   md_temp_alloc(temp, sizeof(float) * ns);
    uint32_t* shell_info     = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * 4 * ns);
    float*    shell_prim     = (float*)   md_temp_alloc(temp, sizeof(float) * 2 * MAX(1, gto_basis_num_shell_prims(basis)));
    uint64_t* key            = (uint64_t*)md_temp_alloc(temp, sizeof(uint64_t) * ns);
    uint32_t* order          = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * ns);
    uint32_t* tao0           = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * ns);
    uint32_t* tao            = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_ao);
    uint8_t*  on             = (uint8_t*) md_temp_alloc(temp, (size_t)nt * nt);
    uint32_t* map            = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * (size_t)nt * nt);
    dm->xyz     = (float*)malloc(sizeof(float) * 3 * MAX(1, num_atoms));
    dm->perm    = (int32_t*)malloc(sizeof(int32_t) * (size_t)nt * GTO_DMAT_T);
    dm->row_ptr = (uint64_t*)calloc((size_t)nt + 1, sizeof(uint64_t));
    bool ok = shell_atom_idx && shell_r && shell_info && shell_prim && key && order && tao0 && tao && on && map && dm->xyz && dm->perm && dm->row_ptr;
    if (ok) {
        gto_expand_basis_gpu_shells(shell_atom_idx, shell_r, shell_info, shell_prim, basis, desc->cutoff);
        const size_t stride = desc->atom_xyz_stride ? desc->atom_xyz_stride : sizeof(float) * 3;
        for (uint32_t a = 0; a < num_atoms; ++a) {
            const float* p = (const float*)((const uint8_t*)desc->atom_xyz + a * stride);
            dm->xyz[3 * a + 0] = p[0]; dm->xyz[3 * a + 1] = p[1]; dm->xyz[3 * a + 2] = p[2];
        }
        // the tiled order: radius class (2 bohr), then atom, then the shell; a shell's AOs stay together
        for (uint32_t s = 0; s < ns; ++s) {
            key[s] = ((uint64_t)(shell_r[s] * 0.5f) << 40) | ((uint64_t)basis->shells[s].atom_idx << 20) | s;
            order[s] = s;
        }
        qsort(key, ns, sizeof(uint64_t), gto_dmat_cmp);
        for (uint32_t i = 0; i < ns; ++i) order[i] = (uint32_t)(key[i] & 0xFFFFFu);
        for (size_t i = 0; i < (size_t)nt * GTO_DMAT_T; ++i) dm->perm[i] = -1;
        uint32_t acc = 0;
        for (uint32_t i = 0; i < ns; ++i) {
            const uint32_t s = order[i];
            const uint32_t nc = (uint32_t)md_gto_num_cart_ao(basis->shells[s].l);
            tao0[s] = acc;
            for (uint32_t c = 0; c < nc; ++c) {
                dm->perm[acc + c] = (int32_t)(shell_info[4 * s + 3] + c);
                tao[shell_info[4 * s + 3] + c] = acc + c;
            }
            acc += nc;
        }
        // the tile pairs of the shell pairs that meet (centres no further apart than their radii, plus margin)
        MEMSET(on, 0, (size_t)nt * nt);
        for (uint32_t s = 0; s < ns; ++s) {
            if (!(shell_r[s] > 0.0f)) continue;
            const float* p = dm->xyz + 3 * shell_atom_idx[s];
            const uint32_t ns_c = (uint32_t)md_gto_num_cart_ao(basis->shells[s].l);
            for (uint32_t u = s; u < ns; ++u) {
                if (!(shell_r[u] > 0.0f)) continue;
                const float* q = dm->xyz + 3 * shell_atom_idx[u];
                const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
                const double lim = (double)shell_r[s] + (double)shell_r[u] + desc->margin;
                if (dx * dx + dy * dy + dz * dz > lim * lim) continue;
                const uint32_t nu_c = (uint32_t)md_gto_num_cart_ao(basis->shells[u].l);
                for (uint32_t i = tao0[s] / GTO_DMAT_T; i <= (tao0[s] + ns_c - 1) / GTO_DMAT_T; ++i) {
                    for (uint32_t j = tao0[u] / GTO_DMAT_T; j <= (tao0[u] + nu_c - 1) / GTO_DMAT_T; ++j) {
                        on[(size_t)MIN(i, j) * nt + MAX(i, j)] = 1;
                    }
                }
            }
        }
        // rows of tiles: their column tiles, work items of up to 8, and the map
        for (uint32_t I = 0; I < nt; ++I) {
            uint64_t c = 0;
            for (uint32_t J = I; J < nt; ++J) c += on[(size_t)I * nt + J];
            dm->row_ptr[I + 1] = dm->row_ptr[I] + c;
            dm->num_items += (c + 7) / 8;
        }
        dm->num_tiles = dm->row_ptr[nt];
        dm->cols = (uint32_t*)malloc(sizeof(uint32_t) * (size_t)MAX(1, dm->num_tiles));
        uint32_t* items = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * 4 * (size_t)MAX(1, dm->num_items));
        ok = dm->cols && items;
        if (ok) {
            uint64_t e = 0, it = 0;
            for (size_t i = 0; i < (size_t)nt * nt; ++i) map[i] = GTO_DMAT_NONE;
            for (uint32_t I = 0; I < nt; ++I) {
                for (uint32_t J = I; J < nt; ++J) {
                    if (!on[(size_t)I * nt + J]) continue;
                    map[(size_t)I * nt + J] = (uint32_t)e;
                    dm->cols[e++] = J;
                }
                for (uint64_t e0 = dm->row_ptr[I]; e0 < dm->row_ptr[I + 1]; e0 += 8) {
                    items[4 * it + 0] = I;
                    items[4 * it + 1] = (uint32_t)e0;
                    items[4 * it + 2] = (uint32_t)MIN(8, dm->row_ptr[I + 1] - e0);
                    items[4 * it + 3] = 0;
                    ++it;
                }
            }
            dm->tao      = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * num_ao).gpu;
            dm->map      = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * (size_t)nt * nt).gpu;
            dm->items    = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * 4 * (size_t)MAX(1, dm->num_items)).gpu;
            dm->cols_buf = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * (size_t)MAX(1, dm->num_tiles)).gpu;
            dm->tiles    = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * 256 * (size_t)MAX(1, dm->num_tiles)).gpu;
            ok = dm->tao && dm->map && dm->items && dm->cols_buf && dm->tiles;
            if (ok) {
                ok = md_gpu_upload(stream, dm->tao, tao, sizeof(uint32_t) * num_ao)
                  && md_gpu_upload(stream, dm->map, map, sizeof(uint32_t) * (size_t)nt * nt)
                  && md_gpu_upload(stream, dm->items, items, sizeof(uint32_t) * 4 * (size_t)MAX(1, dm->num_items))
                  && md_gpu_upload(stream, dm->cols_buf, dm->cols, sizeof(uint32_t) * (size_t)MAX(1, dm->num_tiles))
                  && md_gpu_memset(stream, dm->tiles, 0, sizeof(float) * 256 * (size_t)MAX(1, dm->num_tiles));
            }
        }
    }
    md_temp_end(temp);
    if (!ok) {
        MD_LOG_ERROR("md_gto_gpu_dmat_create: out of memory");
        md_gto_gpu_dmat_destroy(stream, dm);
        return NULL;
    }
    return dm;
}

void md_gto_gpu_dmat_destroy(md_gpu_stream_t stream, md_gto_gpu_dmat_t dm) {
    if (!dm) return;
    md_gpu_free(stream, dm->tao);
    md_gpu_free(stream, dm->map);
    md_gpu_free(stream, dm->items);
    md_gpu_free(stream, dm->cols_buf);
    md_gpu_free(stream, dm->tiles);
    md_gpu_free(stream, dm->P);
    md_gpu_free(stream, dm->w);
    gto_dmat_free_host(dm);
    free(dm);
}

bool md_gto_gpu_dmat_pattern_valid(md_gto_gpu_dmat_t dm, const float* atom_xyz, size_t atom_xyz_stride) {
    if (!dm || !atom_xyz) return false;
    const size_t stride = atom_xyz_stride ? atom_xyz_stride : sizeof(float) * 3;
    const double lim = 0.25 * dm->margin * dm->margin;
    for (uint32_t a = 0; a < dm->num_atoms; ++a) {
        const float* p = (const float*)((const uint8_t*)atom_xyz + a * stride);
        const double dx = p[0] - dm->xyz[3 * a + 0], dy = p[1] - dm->xyz[3 * a + 1], dz = p[2] - dm->xyz[3 * a + 2];
        if (dx * dx + dy * dy + dz * dz > lim) return false;
    }
    return true;
}

static bool gto_dmat_grow(md_gpu_stream_t stream, md_gpu_addr_t* buf, size_t* cap, size_t floats) {
    if (*buf && *cap >= floats) return true;
    md_gpu_free(stream, *buf);
    *buf = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * MAX(1, floats)).gpu;
    *cap = *buf ? floats : 0;
    return *buf != 0;
}

// The kept orbitals (of nonzero weight) packed to float in the tiled AO order: P[k][t]. coeff(k, ao): orbital k's
// coefficient of AO ao, at coeff[k * mo_stride + ao * ao_stride]. Packed by blocks of 64 orbitals and 64 tiled AOs,
// so that orbitals stored as columns ([num_ao][num_mo], ao_stride = num_mo) are read along their rows.
static void gto_dmat_pack(float* P, const md_gto_gpu_dmat_t dm, const double* coeff, size_t mo_stride, size_t ao_stride,
                          const double* weight, size_t num_mo, size_t ntao) {
    enum { B = 64 };
    uint32_t* kept = (uint32_t*)malloc(sizeof(uint32_t) * (num_mo + 1));
    size_t r = 0;
    for (size_t k = 0; k < num_mo; ++k) if (weight[k] != 0.0) kept[r++] = (uint32_t)k;
    for (size_t t0 = 0; t0 < ntao; t0 += B) {
        const size_t t1 = MIN(ntao, t0 + B);
        for (size_t k0 = 0; k0 < r; k0 += B) {
            const size_t k1 = MIN(r, k0 + B);
            for (size_t t = t0; t < t1; ++t) {
                const int32_t ao = dm->perm[t];
                if (ao < 0) {
                    for (size_t k = k0; k < k1; ++k) P[k * ntao + t] = 0.0f;
                    continue;
                }
                const double* src = coeff + (size_t)ao * ao_stride;
                for (size_t k = k0; k < k1; ++k) P[k * ntao + t] = (float)src[(size_t)kept[k] * mo_stride];
            }
        }
    }
    free(kept);
}

static bool gto_dmat_build_packed(md_gpu_stream_t stream, md_gto_gpu_dmat_t dm, const double* coeff, size_t mo_stride, size_t ao_stride,
                                  const double* weight, size_t num_mo) {
    if (!stream || !dm || !coeff || !weight) {
        MD_LOG_ERROR("md_gto_gpu_dmat_build: invalid arguments");
        return false;
    }
    if (!gto_gpu_device) {
        MD_LOG_ERROR("md_gto_gpu_dmat_build: md_gto_gpu_initialize has not been called");
        return false;
    }
    ensure_kernel(gto_gpu_device, &gto_k_dmat_tiles, md_shader_gto_dmat_dmat_tiles_kernel());
    if (!gto_k_dmat_tiles) return false;
    size_t r = 0;
    for (size_t k = 0; k < num_mo; ++k) r += weight[k] != 0.0;
    if (r == 0 || dm->num_tiles == 0) {
        return md_gpu_memset(stream, dm->tiles, 0, sizeof(float) * 256 * (size_t)MAX(1, dm->num_tiles));
    }
    const size_t ntao = (size_t)dm->nt * GTO_DMAT_T;
    if (!gto_dmat_grow(stream, &dm->P, &dm->P_cap, r * ntao) || !gto_dmat_grow(stream, &dm->w, &dm->w_cap, r)) {
        MD_LOG_ERROR("md_gto_gpu_dmat_build: out of device memory (%zu orbitals)", r);
        return false;
    }
    float* P = (float*)md_gpu_upload_begin(stream, dm->P, sizeof(float) * r * ntao);
    if (!P) return false;
    gto_dmat_pack(P, dm, coeff, mo_stride, ao_stride, weight, num_mo, ntao);
    if (!md_gpu_upload_end(stream)) return false;
    float* w = (float*)md_gpu_upload_begin(stream, dm->w, sizeof(float) * r);
    if (!w) return false;
    size_t k1 = 0;
    for (size_t k = 0; k < num_mo; ++k) if (weight[k] != 0.0) w[k1++] = (float)weight[k];
    if (!md_gpu_upload_end(stream)) return false;

    const uint32_t gx = (uint32_t)MIN(dm->num_items, 65535u);
    const uint32_t gy = (uint32_t)DIV_UP(dm->num_items, gx);
    gto_dmat_args_t a = {
        .num_mo = (uint32_t)r, .num_items = (uint32_t)dm->num_items, .items_x = gx, .ntao = (uint32_t)ntao,
        .P = dm->P, .w = dm->w, .items = dm->items, .cols = dm->cols_buf, .tiles = dm->tiles,
    };
    return md_gpu_launch(stream, gto_k_dmat_tiles, md_gpu_grid(gx, gy, 1), &a, sizeof(a));
}

bool md_gto_gpu_dmat_build_columns(md_gpu_stream_t stream, md_gto_gpu_dmat_t dm, const double* coeff, size_t ld, const double* weight, size_t num_mo) {
    if (!dm || ld < num_mo) {
        MD_LOG_ERROR("md_gto_gpu_dmat_build_columns: invalid arguments");
        return false;
    }
    return gto_dmat_build_packed(stream, dm, coeff, 1, ld, weight, num_mo);
}

void md_gto_gpu_dmat_view(md_gto_gpu_dmat_view_t* out, md_gto_gpu_dmat_t dm) {
    if (!out) return;
    MEMSET(out, 0, sizeof(*out));
    if (!dm) return;
    out->num_ao = dm->num_ao;
    out->nt = dm->nt;
    out->num_tiles = dm->num_tiles;
    out->tao = dm->tao;
    out->map = dm->map;
    out->tiles = dm->tiles;
}

bool md_gto_gpu_dmat_build(md_gpu_stream_t stream, md_gto_gpu_dmat_t dm, const double* coeff, const double* weight, size_t num_mo) {
    if (!dm) {
        MD_LOG_ERROR("md_gto_gpu_dmat_build: invalid arguments");
        return false;
    }
    return gto_dmat_build_packed(stream, dm, coeff, dm->num_ao, 1, weight, num_mo);
}

bool md_gto_gpu_dmat_set_matrix(md_gpu_stream_t stream, md_gto_gpu_dmat_t dm, const double* matrix) {
    if (!stream || !dm || !matrix) {
        MD_LOG_ERROR("md_gto_gpu_dmat_set_matrix: invalid arguments");
        return false;
    }
    if (dm->num_tiles == 0) return true;
    float* dst = (float*)md_gpu_upload_begin(stream, dm->tiles, sizeof(float) * 256 * (size_t)dm->num_tiles);
    if (!dst) return false;
    const size_t n = dm->num_ao;
    for (uint32_t I = 0; I < dm->nt; ++I) {
        for (uint64_t e = dm->row_ptr[I]; e < dm->row_ptr[I + 1]; ++e) {
            const uint32_t J = dm->cols[e];
            float* t = dst + e * 256;
            for (uint32_t a = 0; a < GTO_DMAT_T; ++a) {
                const int32_t i = dm->perm[I * GTO_DMAT_T + a];
                for (uint32_t c = 0; c < GTO_DMAT_T; ++c) {
                    const int32_t j = dm->perm[J * GTO_DMAT_T + c];
                    t[a * GTO_DMAT_T + c] = (i < 0 || j < 0) ? 0.0f : (float)(0.5 * (matrix[(size_t)i * n + (size_t)j] + matrix[(size_t)j * n + (size_t)i]));
                }
            }
        }
    }
    return md_gpu_upload_end(stream);
}

void md_gto_gpu_dmat_info(md_gto_gpu_dmat_info_t* out, md_gto_gpu_dmat_t dm) {
    if (!out) return;
    MEMSET(out, 0, sizeof(*out));
    if (!dm) return;
    out->num_ao = dm->num_ao;
    out->num_tiles = (size_t)dm->num_tiles;
    out->num_tiles_total = (size_t)dm->nt * (dm->nt + 1) / 2;
    out->revision = GTO_DMAT_REVISION;
    out->bytes = sizeof(float) * 256 * (size_t)dm->num_tiles + sizeof(uint32_t) * ((size_t)dm->nt * dm->nt + dm->num_ao + 4 * (size_t)dm->num_items + (size_t)dm->num_tiles);
}

size_t md_gto_gpu_atom_buffer_size(size_t num_atoms) {
    return sizeof(float) * 4 * num_atoms;
}

void md_gto_gpu_atom_pack(float* dst_atom_xyzw, const float* atom_xyz, size_t atom_xyz_stride, size_t num_atoms) {
    ASSERT(dst_atom_xyzw && atom_xyz);
    const size_t stride = atom_xyz_stride == 0 ? sizeof(float) * 3 : atom_xyz_stride;
    if (stride == 16) {
        MEMCPY(dst_atom_xyzw, atom_xyz, sizeof(float) * 4 * num_atoms);
    } else {
        for (size_t i = 0; i < num_atoms; ++i) {
            const float* ap = (const float*)((const uint8_t*)atom_xyz + i * stride);
            dst_atom_xyzw[i * 4 + 0] = ap[0];
            dst_atom_xyzw[i * 4 + 1] = ap[1];
            dst_atom_xyzw[i * 4 + 2] = ap[2];
            dst_atom_xyzw[i * 4 + 3] = 0.0f;
        }
    }
}

// ---------------------------------------------------------------------------
// Coefficient buffer helpers
// ---------------------------------------------------------------------------

size_t md_gto_gpu_coeff_size_density(size_t num_cgtos) {
    size_t tri_len = (num_cgtos * (num_cgtos + 1)) / 2;
    return sizeof(float) * tri_len;
}

size_t md_gto_gpu_coeff_size_mo(size_t num_mos, size_t num_cgtos) {
    ASSERT(num_mos > 0);
    ASSERT(num_cgtos > 0);
    return sizeof(float) * num_mos * num_cgtos;
}

void md_gto_gpu_coeff_pack_density(float* dst, const double* density_matrix, size_t num_cgtos) {
    ASSERT(dst && density_matrix && num_cgtos > 0);
    density_matrix_upper_tri_extract_float(dst, density_matrix, num_cgtos);
}

void md_gto_gpu_coeff_pack_mo(float* dst, const double* const* mo_coeffs, const double* mo_scales, size_t num_mos, size_t num_cgtos) {
    ASSERT(dst && mo_coeffs && num_mos > 0 && num_cgtos > 0);
    for (size_t m = 0; m < num_mos; ++m) {
        double scale = mo_scales ? mo_scales[m] : 1.0;
        float* row = dst + m * num_cgtos;
        for (size_t c = 0; c < num_cgtos; ++c)
            row[c] = (float)(mo_coeffs[m][c] * scale);
    }
}

/* md_gto_gpu_coeff_upload_* are gone: md_gpu_upload_begin/end packs straight
   into the destination when that is safe and into staging when it is not, so a
   dedicated helper bought nothing. */

static void gto_fill_common_args(md_gpu_float4x4* w2m, md_gpu_float4x4* i2w, md_gpu_float4* step,
                                 md_gpu_uint4* grid_dim, const md_grid_t* grid, const float sample_offset[3]) {
    world_to_model_matrix((float(*)[4])w2m->m, grid);
    index_to_world_matrix((float(*)[4])i2w->m, grid, sample_offset);
    step->x = grid->spacing.elem[0];
    step->y = grid->spacing.elem[1];
    step->z = grid->spacing.elem[2];
    step->w = 0.0f;
    grid_dim->x = (uint32_t)grid->dim[0];
    grid_dim->y = (uint32_t)grid->dim[1];
    grid_dim->z = (uint32_t)grid->dim[2];
    grid_dim->w = 0;
}

static void gto_density_launch_reference(md_gpu_stream_t stream, const md_gto_gpu_density_desc_t* desc, md_gpu_storage_tex_t out) {
    if (!gto_k_density) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: kernel not initialized");
        return;
    }
    const md_gto_basis_layout_t* L = &desc->basis->layout;
    const md_gpu_addr_t base = desc->basis->buffer;

    gto_density_args_t a = {0};
    gto_fill_common_args(&a.world_to_model, &a.index_to_world, &a.step, &a.grid_dim, desc->grid, desc->sample_offset);
    a.num_cgtos = L->num_cgtos;
    a.operation = (uint32_t)desc->op;

    /* Device memory is a pointer, so the basis sub-ranges are plain
       arithmetic rather than (buffer, offset, usage) triples. */
    a.cgto_atom_idx = base + L->off_cgto_atom_idx;
    a.cgto_r        = base + L->off_cgto_r;
    a.cgto_off_len  = base + L->off_cgto_off_len;
    a.pgto          = base + L->off_pgto;
    a.atom_xyz      = desc->atom_xyz;
    a.D_matrix      = desc->dmat ? desc->dmat->tiles : desc->coeff;
    if (desc->dmat) {
        a.dm_tiled = 1;
        a.dm_nt    = desc->dmat->nt;
        a.dm_tao   = desc->dmat->tao;
        a.dm_map   = desc->dmat->map;
    }
    a.out_tex       = out;

    const md_gpu_grid_t g = md_gpu_grid_for(gto_k_density, a.grid_dim.x, a.grid_dim.y, a.grid_dim.z);
    md_gpu_launch(stream, gto_k_density, g, &a, sizeof(a));
}

#define GTO_GEMM_BLOCK_VOXELS 512u
#define GTO_GEMM_GROUP_VOXELS 128u   // voxels per group of the GEMM-pass kernels
#define GTO_GEMM_DEFAULT_SCRATCH ((size_t)256 << 20)

typedef struct {
    md_gpu_kernel_t large;      // kernel for blocks with more than small_n AOs (and all blocks if small is NULL)
    md_gpu_kernel_t small;      // optional kernel for blocks with at most small_n AOs
    uint32_t large_gp;          // voxels per group of 'large'
    uint32_t small_gp;
    uint32_t small_n;
    char     label[64];
} gto_gemm_cfg_t;

/* Automatic density GEMM configuration per GPU vendor (AO tile width, and the block
   size up to which the 32-wide kernel is used instead). Large systems (C60 on a fine
   grid, C60 with d/f shells, C240), time relative to the best configuration:
     64-wide + 32 for n <= 32:  GTX 1060 1.09, L4 1.00, M4 Pro 1.17, HD 530 1.00
     32-wide:                   GTX 1060 1.00, L4 1.18, M4 Pro 1.00, HD 530 1.25
   AMD is unmeasured and gets the NVIDIA choice. */
static void gto_gemm_auto_config(uint32_t* gm, int32_t* small_n) {
    switch (gto_gpu_vendor) {
    case 0x106B:   /* Apple */
        *gm = 32; *small_n = -1; break;
    case 0x10DE:   /* NVIDIA */
    case 0x8086:   /* Intel */
    case 0x1002:   /* AMD */
    default:
        *gm = 64; *small_n = 32; break;
    }
}

static bool gto_gemm_config(gto_gemm_cfg_t* cfg, const md_gto_gpu_density_desc_t* desc) {
    MEMSET(cfg, 0, sizeof(*cfg));
    uint32_t gm = desc->gemm_tile;
    int32_t  sn = desc->gemm_small_block;
    uint32_t agm; int32_t asn;
    gto_gemm_auto_config(&agm, &asn);
    if (gm == 0) gm = agm;
    if (sn == 0) sn = asn;
    if (gm == 32) sn = -1;   // already the small tile

    cfg->large    = gto_density_gemm2_kernel(gm);
    cfg->large_gp = GTO_GEMM_GROUP_VOXELS;
    if (!cfg->large) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: no GEMM kernel for a %u-wide AO tile (32 or 64)", gm);
        return false;
    }
    if (sn > 0) {
        cfg->small    = gto_density_gemm2_kernel(32);
        cfg->small_gp = GTO_GEMM_GROUP_VOXELS;
        cfg->small_n  = (uint32_t)sn;
        md_snprintf(cfg->label, sizeof(cfg->label), "%ux%u, %ux32 for n<=%d", GTO_GEMM_GROUP_VOXELS, gm, GTO_GEMM_GROUP_VOXELS, sn);
    } else {
        md_snprintf(cfg->label, sizeof(cfg->label), "%ux%u", GTO_GEMM_GROUP_VOXELS, gm);
    }
    return true;
}

// Fills the parts of the GEMM-path argument block shared by density and orbitals.
static void gto_gemm_args_init(gto_density_gemm_args_t* a, const md_grid_t* grid, const float sample_offset[3],
                               md_gto_gpu_basis_t gb, md_gpu_addr_t atom_xyz, md_gto_op_t op, md_gpu_storage_tex_t out) {
    MEMSET(a, 0, sizeof(*a));
    const md_gto_basis_layout_t* L = &gb->layout;
    const md_gpu_addr_t base = gb->buffer;
    gto_fill_common_args(&a->world_to_model, &a->index_to_world, &a->step, &a->grid_dim, grid, sample_offset);
    a->block_dim.x    = DIV_UP(a->grid_dim.x, 8);
    a->block_dim.y    = DIV_UP(a->grid_dim.y, 8);
    a->block_dim.z    = DIV_UP(a->grid_dim.z, 8);
    a->num_shells     = L->num_shells;
    a->num_aos        = L->num_cgtos;
    a->operation      = (uint32_t)op;
    a->shell_atom_idx = base + L->off_shell_atom_idx;
    a->shell_r        = base + L->off_shell_r;
    a->shell_info     = base + L->off_shell_info;
    a->shell_prim     = base + L->off_shell_prim;
    a->atom_xyz       = atom_xyz;
    a->out_tex        = out;
}

// The two-pass path: count, plan batches (blocks the calling thread for the counts),
// then per batch the Phi pass, optionally the D_b gather, and the consumer kernel(s)
// in cfg. `a` must be initialised (gto_gemm_args_init) with the consumer's fields set.
/* The GEMM path in two halves: the counting pass, whose per-block counts come back to the host, and the
   batches planned from them. md_gto_gpu_density_launch waits for the counts between the two on the calling
   thread; a job (md_gto_gpu_density_job_begin) issues the second half when they are there. */
typedef struct gto_gemm_pending_t {
    md_gpu_stream_t         stream;
    gto_density_gemm_args_t a;
    gto_gemm_cfg_t          cfg;
    size_t                  scratch_bytes;
    bool                    pregather_D;
    uint32_t                num_blocks;
    md_gpu_addr_t           counts;
    md_gpu_mem_t            counts_rb;
    md_gpu_sync_t           sync;       /* the counts on the host */
    md_tick_t               t_count0;
    double                  waited_ms;  /* the calling thread's wait for them, the launch's */
} gto_gemm_pending_t;

/* Pass 0: per-block counts, copied to the host to plan the batches. False when there is nothing to do or it
   could not be issued (logged). */
static bool gto_gemm_count(gto_gemm_pending_t* p, md_gpu_stream_t stream, const gto_density_gemm_args_t* ap, size_t scratch_bytes, bool pregather_D, const gto_gemm_cfg_t* cfgp) {
    MEMSET(p, 0, sizeof(*p));
    if (!gto_density_gemm_kernels()) {
        MD_LOG_ERROR("md_gto: GEMM-path kernels not available (md_gto_gpu_initialize not called?)");
        return false;
    }
    p->stream = stream;
    p->a = *ap;
    p->cfg = *cfgp;
    p->scratch_bytes = scratch_bytes;
    p->pregather_D = pregather_D;
    p->a.d_mode = pregather_D ? 1 : 0;

    const uint32_t num_blocks = p->a.block_dim.x * p->a.block_dim.y * p->a.block_dim.z;
    if (num_blocks == 0) return false;
    p->num_blocks = num_blocks;

    const size_t counts_size = sizeof(uint32_t) * 2 * num_blocks;
    p->counts = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, counts_size).gpu;
    p->counts_rb = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, counts_size);
    if (!p->counts || !p->counts_rb.cpu) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: failed to allocate block counts");
        md_gpu_free(stream, p->counts);
        md_gpu_free(stream, p->counts_rb.gpu);
        return false;
    }
    p->t_count0 = md_tick_now();
    p->a.block_counts = p->counts;
    p->a.num_blocks   = num_blocks;
    md_gpu_launch(stream, gto_k_density_gemm_count, md_gpu_grid(p->a.block_dim.x, p->a.block_dim.y, p->a.block_dim.z), &p->a, sizeof(p->a));
    md_gpu_copy(stream, p->counts_rb.gpu, p->counts, counts_size);
    p->sync = md_gpu_stream_record(stream);
    return true;
}

/* What the counting pass left, dropped without the rest (stream-ordered) */
static void gto_gemm_drop(gto_gemm_pending_t* p) {
    md_gpu_free(p->stream, p->counts);
    md_gpu_free(p->stream, p->counts_rb.gpu);
    p->counts = 0;
    p->counts_rb = (md_gpu_mem_t){0};
}

/* The batches, planned from the counts (which must be on the host: p->sync complete) */
static void gto_gemm_issue(gto_gemm_pending_t* p) {
    md_gpu_stream_t stream = p->stream;
    gto_density_gemm_args_t a = p->a;
    const gto_gemm_cfg_t cfg = p->cfg;
    const bool pregather_D = p->pregather_D;
    const size_t scratch_bytes = p->scratch_bytes;
    const uint32_t num_blocks = p->num_blocks;
    const md_gpu_mem_t counts_rb = p->counts_rb;
    const md_tick_t t_count0 = md_tick_now();
    md_gpu_free(stream, p->counts);
    p->counts = 0;

    /* Plan: non-empty blocks in linear order, batched so that each batch's Phi (and
       D_b) fits in the scratch budget. Within a batch the blocks for the small-tile
       kernel come first, so each kernel gets one contiguous range of the table.
       Empty blocks go into one trailing GEMM-only launch that just applies the
       operation with rho = 0. */
    const uint32_t* cnt = (const uint32_t*)counts_rb.cpu;
    double t_count = 0.0;
    const size_t budget = scratch_bytes ? scratch_bytes : GTO_GEMM_DEFAULT_SCRATCH;
    typedef struct { uint32_t first, count, num_small; uint64_t n_ao, n_db; } batch_t;
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* table   = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * 8 * num_blocks);
    uint32_t* order   = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_blocks);
    batch_t*  batches = (batch_t*) md_temp_alloc(temp, sizeof(batch_t) * (num_blocks + 1));
    uint32_t  num_batches = 0;
    uint32_t  num_entries = 0;
    uint64_t  max_ao = 0, max_db = 0;
    {
        uint32_t num_nonempty = 0;
        for (uint32_t b = 0; b < num_blocks; ++b) {
            if (cnt[2 * b + 1] != 0) order[num_nonempty++] = b;
        }
        uint32_t i0 = 0;
        while (i0 < num_nonempty) {
            /* Extend the batch while it fits (a single oversized block gets its own). */
            uint64_t bytes = 0;
            uint32_t i1 = i0;
            while (i1 < num_nonempty) {
                const uint64_t na = cnt[2 * order[i1] + 1];
                const uint64_t bb = sizeof(float) * (na * GTO_GEMM_BLOCK_VOXELS + (pregather_D ? na * na : 0));
                if (i1 > i0 && bytes + bb > budget) break;
                bytes += bb;
                i1++;
            }
            batch_t bt = { .first = num_entries };
            for (int pass = 0; pass < 2; ++pass) {
                for (uint32_t i = i0; i < i1; ++i) {
                    const uint32_t b  = order[i];
                    const uint32_t na = cnt[2 * b + 1];
                    const bool is_small = cfg.small && na <= cfg.small_n;
                    if (is_small != (pass == 0)) continue;
                    uint32_t* e = table + 8 * num_entries;
                    e[0] = b % a.block_dim.x;
                    e[1] = (b / a.block_dim.x) % a.block_dim.y;
                    e[2] = b / (a.block_dim.x * a.block_dim.y);
                    e[3] = na;
                    e[4] = (uint32_t)(bt.n_ao * GTO_GEMM_BLOCK_VOXELS);   // phi offset (floats)
                    e[5] = (uint32_t)bt.n_ao;                             // ao_list offset
                    e[6] = (uint32_t)bt.n_db;                             // D_b offset (floats)
                    e[7] = 0;
                    bt.count     += 1;
                    bt.num_small += is_small ? 1 : 0;
                    bt.n_ao      += na;
                    bt.n_db      += pregather_D ? (uint64_t)na * na : 0;
                    num_entries  += 1;
                }
            }
            max_ao = MAX(max_ao, bt.n_ao);
            max_db = MAX(max_db, bt.n_db);
            batches[num_batches++] = bt;
            i0 = i1;
        }
        /* Empty blocks */
        batch_t empty = { .first = num_entries };
        for (uint32_t b = 0; b < num_blocks; ++b) {
            if (cnt[2 * b + 1] != 0) continue;
            uint32_t* e = table + 8 * num_entries;
            e[0] = b % a.block_dim.x;
            e[1] = (b / a.block_dim.x) % a.block_dim.y;
            e[2] = b / (a.block_dim.x * a.block_dim.y);
            e[3] = 0; e[4] = 0; e[5] = 0; e[6] = 0; e[7] = 0;
            empty.count += 1;
            num_entries += 1;
        }
        batches[num_batches] = empty;   // stored one past the real batches
    }
    md_gpu_free(stream, counts_rb.gpu);

    const size_t table_size = sizeof(uint32_t) * 8 * MAX(1, num_entries);
    md_gpu_addr_t table_buf = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, table_size).gpu;
    md_gpu_addr_t phi_buf   = max_ao ? md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * GTO_GEMM_BLOCK_VOXELS * max_ao).gpu : 0;
    md_gpu_addr_t ao_buf    = max_ao ? md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(uint32_t) * max_ao).gpu : 0;
    md_gpu_addr_t db_buf    = max_db ? md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * max_db).gpu : 0;
    if (!table_buf || (max_ao && (!phi_buf || !ao_buf)) || (max_db && !db_buf)) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: failed to allocate GEMM scratch (%.1f MB)",
                     (double)(sizeof(float) * (GTO_GEMM_BLOCK_VOXELS * max_ao + max_db)) / (1024.0 * 1024.0));
        goto done;
    }
    md_gpu_upload(stream, table_buf, table, table_size);
    t_count = p->waited_ms + md_tick_to_milliseconds(md_tick_now() - t_count0);

    a.ao_list    = ao_buf;
    a.phi        = phi_buf;
    a.D_blocks   = db_buf;

    /* Diagnostics: MD_GTO_GEMM_PROFILE=1 synchronises after every pass and logs
       where the time goes. Never set in normal use. */
    const bool profile = getenv("MD_GTO_GEMM_PROFILE") != NULL;
    double t_phi = 0.0, t_gemm = 0.0;
    if (profile) md_gpu_stream_sync(stream);

    for (uint32_t i = 0; i <= num_batches; ++i) {
        const batch_t* bt = &batches[i];
        if (bt->count == 0) continue;
        if (i < num_batches) {
            md_tick_t t0 = md_tick_now();
            a.num_blocks  = bt->count;
            a.block_table = table_buf + sizeof(uint32_t) * 8 * bt->first;
            md_gpu_launch(stream, gto_k_density_gemm_phi, md_gpu_grid(bt->count, 1, 1), &a, sizeof(a));
            if (pregather_D) {
                md_gpu_launch(stream, gto_k_density_gemm_dgather, md_gpu_grid(bt->count, 1, 1), &a, sizeof(a));
            }
            if (profile) { md_gpu_stream_sync(stream); t_phi += md_tick_to_milliseconds(md_tick_now() - t0); }
        }
        md_tick_t t0 = md_tick_now();
        if (bt->num_small > 0) {
            a.num_blocks  = bt->num_small;
            a.block_table = table_buf + sizeof(uint32_t) * 8 * bt->first;
            md_gpu_launch(stream, cfg.small, md_gpu_grid(bt->num_small * (GTO_GEMM_BLOCK_VOXELS / cfg.small_gp), 1, 1), &a, sizeof(a));
        }
        const uint32_t num_large = bt->count - bt->num_small;
        if (num_large > 0) {
            a.num_blocks  = num_large;
            a.block_table = table_buf + sizeof(uint32_t) * 8 * (bt->first + bt->num_small);
            md_gpu_launch(stream, cfg.large, md_gpu_grid(num_large * (GTO_GEMM_BLOCK_VOXELS / cfg.large_gp), 1, 1), &a, sizeof(a));
        }
        if (profile) { md_gpu_stream_sync(stream); t_gemm += md_tick_to_milliseconds(md_tick_now() - t0); }
    }
    if (profile) {
        uint64_t total_ao = 0, total_small = 0;
        for (uint32_t i = 0; i < num_batches; ++i) { total_ao += batches[i].n_ao; total_small += batches[i].num_small; }
        const uint32_t nonempty = num_entries - batches[num_batches].count;
        MD_LOG_INFO("GTO %s (GEMM) profile: count+plan %.2f ms, phi%s %.2f ms, gemm[%s] %.2f ms | %u batches, avg %.0f AOs per non-empty block, %.0f%% small",
                    a.num_mos ? "orbitals" : "density", t_count, pregather_D ? "+dgather" : "", t_phi, cfg.label, t_gemm, num_batches,
                    (double)total_ao / MAX(1, nonempty), 100.0 * (double)total_small / MAX(1, nonempty));
    }

    MD_LOG_DEBUG("GTO %s (GEMM): %u blocks, %u non-empty, %u batches, max %.1f MB Phi per batch; the calling thread spent %.1f ms on the counts (waiting for them, and what the stream held before them, unless a job) and the plan",
                 a.num_mos ? "orbitals" : "density", num_blocks, num_entries - batches[num_batches].count, num_batches,
                 (double)(sizeof(float) * GTO_GEMM_BLOCK_VOXELS * max_ao) / (1024.0 * 1024.0), t_count);
done:
    md_gpu_free(stream, db_buf);
    md_gpu_free(stream, ao_buf);
    md_gpu_free(stream, phi_buf);
    md_gpu_free(stream, table_buf);
    md_temp_end(temp);
    p->counts_rb = (md_gpu_mem_t){0};
}

static void gto_gemm_execute(md_gpu_stream_t stream, gto_density_gemm_args_t* ap, size_t scratch_bytes, bool pregather_D, const gto_gemm_cfg_t* cfgp) {
    gto_gemm_pending_t p;
    if (!gto_gemm_count(&p, stream, ap, scratch_bytes, pregather_D, cfgp)) return;
    /* This is the one place the GEMM path blocks the calling thread */
    md_gpu_sync_wait(p.sync);
    p.waited_ms = md_tick_to_milliseconds(md_tick_now() - p.t_count0);
    gto_gemm_issue(&p);
}

static bool gto_density_gemm_setup(gto_density_gemm_args_t* a, gto_gemm_cfg_t* cfg, const md_gto_gpu_density_desc_t* desc, md_gpu_storage_tex_t out) {
    if (!gto_density_gemm_kernels()) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: GEMM kernels not available (md_gto_gpu_initialize not called?)");
        return false;
    }
    if (!gto_gemm_config(cfg, desc)) return false;
    gto_gemm_args_init(a, desc->grid, desc->sample_offset, desc->basis, desc->atom_xyz, desc->op, out);
    a->D_matrix = desc->dmat ? desc->dmat->tiles : desc->coeff;
    if (desc->dmat) {
        a->dm_tiled = 1;
        a->dm_nt    = desc->dmat->nt;
        a->dm_tao   = desc->dmat->tao;
        a->dm_map   = desc->dmat->map;
    }
    return true;
}

static void gto_density_launch_gemm(md_gpu_stream_t stream, const md_gto_gpu_density_desc_t* desc, md_gpu_storage_tex_t out) {
    gto_gemm_cfg_t cfg;
    gto_density_gemm_args_t a;
    if (!gto_density_gemm_setup(&a, &cfg, desc, out)) return;
    /* D_b is gathered once per block before the GEMM pass (d_mode 1), so the GEMM
       streams it with coalesced loads instead of gathering from the packed triangle. */
    gto_gemm_execute(stream, &a, desc->scratch_bytes, true, &cfg);
}

/* The storage the evaluation writes and its algorithm, or false (logged) */
static bool gto_density_check(md_gpu_storage_tex_t* out, md_gto_gpu_density_algo_t* algo, md_gpu_stream_t stream, const md_gto_gpu_density_desc_t* desc) {
    if (!stream || !desc || !desc->basis || !desc->atom_xyz || !(desc->coeff || desc->dmat) || !desc->out_tex || !desc->grid) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: invalid input");
        return false;
    }
    if (desc->dmat && desc->dmat->num_ao != desc->basis->layout.num_cgtos) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: the dmat spans %u AOs, the basis %u", desc->dmat->num_ao, desc->basis->layout.num_cgtos);
        return false;
    }
    *out = md_gpu_texture_storage(desc->out_tex, 0);
    if (!out->handle) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: out_tex needs MD_GPU_TEX_STORAGE usage");
        return false;
    }
    *algo = desc->algorithm;
    if (*algo == MD_GTO_GPU_DENSITY_ALGO_DEFAULT) {
        *algo = MD_GTO_GPU_DENSITY_ALGO_GEMM;   // 1.8-2.6x faster than the reference on large systems
    }
    if (*algo != MD_GTO_GPU_DENSITY_ALGO_REFERENCE && *algo != MD_GTO_GPU_DENSITY_ALGO_GEMM) {
        MD_LOG_ERROR("md_gto_gpu_density_launch: unknown algorithm %d", (int)*algo);
        return false;
    }
    return true;
}

void md_gto_gpu_density_launch(md_gpu_stream_t stream, const md_gto_gpu_density_desc_t* desc) {
    md_gpu_storage_tex_t out;
    md_gto_gpu_density_algo_t algo;
    if (!gto_density_check(&out, &algo, stream, desc)) return;
    if (algo == MD_GTO_GPU_DENSITY_ALGO_REFERENCE) {
        gto_density_launch_reference(stream, desc, out);
    } else {
        gto_density_launch_gemm(stream, desc, out);
    }
}

struct md_gto_gpu_job_t {
    gto_gemm_pending_t p;
    bool issued;        // all of it on the stream
};

md_gto_gpu_job_t* md_gto_gpu_density_job_begin(md_gpu_stream_t stream, const md_gto_gpu_density_desc_t* desc) {
    md_gpu_storage_tex_t out;
    md_gto_gpu_density_algo_t algo;
    if (!gto_density_check(&out, &algo, stream, desc)) return NULL;
    md_gto_gpu_job_t* job = (md_gto_gpu_job_t*)calloc(1, sizeof(md_gto_gpu_job_t));
    if (!job) return NULL;
    if (algo == MD_GTO_GPU_DENSITY_ALGO_REFERENCE) {
        gto_density_launch_reference(stream, desc, out);
        job->issued = true;
        return job;
    }
    gto_gemm_cfg_t cfg;
    gto_density_gemm_args_t a;
    // nothing to evaluate, or what could not be issued (logged): done, as the launch would be
    job->issued = !gto_density_gemm_setup(&a, &cfg, desc, out) || !gto_gemm_count(&job->p, stream, &a, desc->scratch_bytes, true, &cfg);
    return job;
}

bool md_gto_gpu_job_poll(md_gto_gpu_job_t* job, bool wait) {
    if (!job) return true;
    if (job->issued) return true;
    if (!md_gpu_sync_is_complete(job->p.sync)) {
        if (!wait) return false;
        const md_tick_t t0 = md_tick_now();
        md_gpu_sync_wait(job->p.sync);
        job->p.waited_ms = md_tick_to_milliseconds(md_tick_now() - t0);
    }
    gto_gemm_issue(&job->p);
    job->issued = true;
    return true;
}

void md_gto_gpu_job_end(md_gto_gpu_job_t* job) {
    if (!job) return;
    if (!job->issued) gto_gemm_drop(&job->p);
    free(job);
}

/* Mirrors RootArgs in eval_gto_mo_shell.slang. */
typedef struct {
    md_gpu_float4x4 world_to_model;
    md_gpu_float4x4 index_to_world;
    md_gpu_float4   step;
    md_gpu_uint4    grid_dim;
    uint32_t        num_shells;
    uint32_t        num_aos;
    uint32_t        num_mos;
    uint32_t        mode;
    uint32_t        operation;
    uint32_t        has_weights;
    uint32_t        coeff_screening;
    uint32_t        _pad0;
    md_gpu_addr_t   shell_atom_idx;
    md_gpu_addr_t   shell_r;
    md_gpu_addr_t   shell_info;
    md_gpu_addr_t   shell_prim;
    md_gpu_addr_t   atom_xyz;
    md_gpu_addr_t   coeffs;
    md_gpu_addr_t   weights;
    md_gpu_storage_tex_t out_tex;
} gto_mo_shell_args_t;

#define GTO_MO_SHELL_MB 8   // orbitals per pass of the multi-orbital shell kernel

/* The shell kernel with 4 voxels per thread, for one orbital (mb = 1) or up to
   GTO_MO_SHELL_MB per pass. 4 voxels per thread was the best or within 10% of it on
   every measured GPU for single orbitals (1 and 2 voxels per thread: 1.1-2.2x slower). */
static md_gpu_kernel_t gto_mo_shell_kernel(uint32_t mb) {
    if (!gto_gpu_device) return NULL;
    const int im = mb == 1 ? 0 : 1;
    md_gpu_kernel_t* slot = &gto_k_mo_shell[im];
    if (!*slot) {
        ensure_kernel(gto_gpu_device, slot, im == 0 ? md_shader_eval_gto_mo_shell_v4_m1_main_kernel()
                                                    : md_shader_eval_gto_mo_shell_v4_m8_main_kernel());
    }
    return *slot;
}

/* The GEMM consumer for orbitals: 32 orbitals per tile (64 was slower everywhere). */
static md_gpu_kernel_t gto_mo_gemm_kernel(void) {
    if (!gto_gpu_device) return NULL;
    ensure_kernel(gto_gpu_device, &gto_k_mo_gemm, md_shader_eval_gto_mo_gemm_128x32_main_kernel());
    return gto_k_mo_gemm;
}

static void gto_orbital_launch_reference(md_gpu_stream_t stream, const md_gto_gpu_orbital_desc_t* desc, md_gpu_storage_tex_t out) {
    if (!gto_k_mo) {
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: kernel not initialized");
        return;
    }
    if (desc->weights) {
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: the reference kernel does not support weights");
        return;
    }
    const md_gto_basis_layout_t* L = &desc->basis->layout;
    const md_gpu_addr_t base = desc->basis->buffer;

    gto_mo_args_t a = {0};
    gto_fill_common_args(&a.world_to_model, &a.index_to_world, &a.step, &a.grid_dim, desc->grid, desc->sample_offset);
    a.num_cgtos = L->num_cgtos;
    a.num_rows  = (uint32_t)desc->num_orbitals;
    a.mode      = (uint32_t)desc->eval_mode;
    a.operation = (uint32_t)desc->op;

    a.cgto_atom_idx = base + L->off_cgto_atom_idx;
    a.cgto_r        = base + L->off_cgto_r;
    a.cgto_off_len  = base + L->off_cgto_off_len;
    a.pgto          = base + L->off_pgto;
    a.atom_xyz      = desc->atom_xyz;
    a.coeffs        = desc->coeff;
    a.out_tex       = out;

    const md_gpu_grid_t g = md_gpu_grid_for(gto_k_mo, a.grid_dim.x, a.grid_dim.y, a.grid_dim.z);
    md_gpu_launch(stream, gto_k_mo, g, &a, sizeof(a));
}

static void gto_orbital_launch_shell(md_gpu_stream_t stream, const md_gto_gpu_orbital_desc_t* desc, md_gpu_storage_tex_t out) {
    const uint32_t mb = desc->num_orbitals == 1 ? 1 : GTO_MO_SHELL_MB;
    md_gpu_kernel_t k = gto_mo_shell_kernel(mb);
    if (!k) {
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: shell kernel not available (md_gto_gpu_initialize not called?)");
        return;
    }
    const md_gto_basis_layout_t* L = &desc->basis->layout;
    const md_gpu_addr_t base = desc->basis->buffer;

    gto_mo_shell_args_t a = {0};
    gto_fill_common_args(&a.world_to_model, &a.index_to_world, &a.step, &a.grid_dim, desc->grid, desc->sample_offset);
    a.num_shells      = L->num_shells;
    a.num_aos         = L->num_cgtos;
    a.num_mos         = (uint32_t)desc->num_orbitals;
    a.mode            = (uint32_t)desc->eval_mode;
    a.operation       = (uint32_t)desc->op;
    a.has_weights     = desc->weights ? 1 : 0;
    a.coeff_screening = desc->exact_screening ? 0 : 1;
    a.shell_atom_idx  = base + L->off_shell_atom_idx;
    a.shell_r         = base + L->off_shell_r;
    a.shell_info      = base + L->off_shell_info;
    a.shell_prim      = base + L->off_shell_prim;
    a.atom_xyz        = desc->atom_xyz;
    a.coeffs          = desc->coeff;
    a.weights         = desc->weights;
    a.out_tex         = out;

    const md_gpu_grid_t g = md_gpu_grid(DIV_UP(a.grid_dim.x, 8), DIV_UP(a.grid_dim.y, 8), DIV_UP(a.grid_dim.z, 8));
    md_gpu_launch(stream, k, g, &a, sizeof(a));
}

static void gto_orbital_launch_gemm(md_gpu_stream_t stream, const md_gto_gpu_orbital_desc_t* desc, md_gpu_storage_tex_t out) {
    gto_gemm_cfg_t cfg = {0};
    cfg.large    = gto_mo_gemm_kernel();
    cfg.large_gp = GTO_GEMM_GROUP_VOXELS;
    md_snprintf(cfg.label, sizeof(cfg.label), "mo %ux32", GTO_GEMM_GROUP_VOXELS);
    if (!cfg.large) {
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: GEMM kernel not available (md_gto_gpu_initialize not called?)");
        return;
    }

    gto_density_gemm_args_t a;
    gto_gemm_args_init(&a, desc->grid, desc->sample_offset, desc->basis, desc->atom_xyz, desc->op, out);
    a.D_matrix    = desc->coeff;
    a.num_mos     = (uint32_t)desc->num_orbitals;
    a.mo_mode     = (uint32_t)desc->eval_mode;
    a.has_weights = desc->weights ? 1 : 0;
    a.mo_weights  = desc->weights;
    gto_gemm_execute(stream, &a, desc->scratch_bytes, false, &cfg);
}

void md_gto_gpu_orbital_launch(md_gpu_stream_t stream, const md_gto_gpu_orbital_desc_t* desc) {
    if (!stream || !desc || !desc->basis || !desc->atom_xyz || !desc->coeff || !desc->out_tex || !desc->grid) {
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: invalid input");
        return;
    }
    if (desc->num_orbitals == 0) return;
    const md_gpu_storage_tex_t out = md_gpu_texture_storage(desc->out_tex, 0);
    if (!out.handle) {
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: out_tex needs MD_GPU_TEX_STORAGE usage");
        return;
    }

    md_gto_gpu_orbital_algo_t algo = desc->algorithm;
    if (algo == MD_GTO_GPU_ORBITAL_ALGO_DEFAULT) {
        /* The shell kernel re-evaluates the AOs once per GTO_MO_SHELL_MB orbitals; past
           that, evaluating them once and multiplying by the coefficients wins on most
           GPUs (for 32 orbitals on GTX 1060, M4 Pro and HD 530; not on L4, where the
           shell kernel was still 1.6-2.2x faster). The crossover is not measured. */
        algo = desc->num_orbitals <= GTO_MO_SHELL_MB ? MD_GTO_GPU_ORBITAL_ALGO_SHELL : MD_GTO_GPU_ORBITAL_ALGO_GEMM;
    }
    switch (algo) {
    case MD_GTO_GPU_ORBITAL_ALGO_REFERENCE: gto_orbital_launch_reference(stream, desc, out); break;
    case MD_GTO_GPU_ORBITAL_ALGO_SHELL:     gto_orbital_launch_shell(stream, desc, out);     break;
    case MD_GTO_GPU_ORBITAL_ALGO_GEMM:      gto_orbital_launch_gemm(stream, desc, out);      break;
    default:
        MD_LOG_ERROR("md_gto_gpu_orbital_launch: unknown algorithm %d", (int)algo);
        break;
    }
}

#endif // MD_ENABLE_GPU

static inline float fast_powf(float base, int exp) {
    float val = 1.0f;
    switch(exp) {
    case 4: val *= base; FALLTHROUGH;
    case 3: val *= base; FALLTHROUGH;
    case 2: val *= base; FALLTHROUGH;
    case 1: val *= base; FALLTHROUGH;
    case 0: break;
    }
    return val;
}

static inline double fast_pow(double base, int exp){
    double val = 1.0;
    switch(exp) {
    case 4: val *= base; FALLTHROUGH;
    case 3: val *= base; FALLTHROUGH;
    case 2: val *= base; FALLTHROUGH;
    case 1: val *= base; FALLTHROUGH;
    case 0: break;
    }
    return val;
}


static inline md_128 md_mm_fast_pow1(md_128 base, int exp) {
    switch (exp) {
    case 1:
        return base;
    case 2:
        return md_mm_mul_ps(base, base);
    case 3:
        return md_mm_mul_ps(base, md_mm_mul_ps(base, base));
    case 4: {
        md_128 squared = md_mm_mul_ps(base, base);
        return md_mm_mul_ps(squared, squared);
    }
    case 0:
    default:
        return md_mm_set1_ps(1.0f);
    }
}

static inline md_256 md_mm256_fast_pow1(md_256 base, int exp) {
    switch (exp) {
    case 1:
        return base;
    case 2:
        return md_mm256_mul_ps(base, base);
    case 3:
        return md_mm256_mul_ps(base, md_mm256_mul_ps(base, base));
    case 4: {
        md_256 squared = md_mm256_mul_ps(base, base);
        return md_mm256_mul_ps(squared, squared);
    }
    case 0:
    default:
        return md_mm256_set1_ps(1.0f);
    }
}

static inline md_128 md_mm_fast_pow(md_128 base1, md_128i exp) {
    md_128 base2 = md_mm_mul_ps(base1, base1);
    md_128 base3 = md_mm_mul_ps(base2, base1);
    md_128 base4 = md_mm_mul_ps(base2, base2);

    md_128 mask1 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(1)));
    md_128 mask2 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(2)));
    md_128 mask3 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(3)));
    md_128 mask4 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(4)));

    md_128 res = md_mm_set1_ps(1.0f);
    res = md_mm_select_ps(res, base1, mask1);
    res = md_mm_select_ps(res, base2, mask2);
    res = md_mm_select_ps(res, base3, mask3);
    res = md_mm_select_ps(res, base4, mask4);
    return res;
}

static inline md_256 md_mm256_fast_pow(md_256 base1, md_256i exp) {
    md_256 base2 = md_mm256_mul_ps(base1, base1);
    md_256 base3 = md_mm256_mul_ps(base2, base1);
    md_256 base4 = md_mm256_mul_ps(base2, base2);

    md_256 mask1 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(1)));
    md_256 mask2 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(2)));
    md_256 mask3 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(3)));
    md_256 mask4 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(4)));

    md_256 res = md_mm256_set1_ps(1.0f);
    res = md_mm256_select_ps(res, base1, mask1);
    res = md_mm256_select_ps(res, base2, mask2);
    res = md_mm256_select_ps(res, base3, mask3);
    res = md_mm256_select_ps(res, base4, mask4);
    return res;
}

#ifdef __AVX512F__
static inline __m512 md_mm512_fast_pow(__m512 base1, __m512i exp) {
    __m512 base2 = _mm512_mul_ps(base1,  base1);
    __m512 base3 = _mm512_mul_ps(base2, base1);
    __m512 base4 = _mm512_mul_ps(base2, base2);

    __mmask16 mask1 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(1), _MM_CMPINT_EQ);
    __mmask16 mask2 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(2), _MM_CMPINT_EQ);
    __mmask16 mask3 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(3), _MM_CMPINT_EQ);
    __mmask16 mask4 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(4), _MM_CMPINT_EQ);

    __m512 res = _mm512_set1_ps(1.0f);
    res = _mm512_mask_blend_ps(mask1, res, base1);
    res = _mm512_mask_blend_ps(mask2, res, base2);
    res = _mm512_mask_blend_ps(mask3, res, base3);
    res = _mm512_mask_blend_ps(mask4, res, base4);
    return res;
}
#endif

static inline void evaluate_grid_ref(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
        const int z_stride = iz * grid_dim[0] * grid_dim[1];
        for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
            const int y_stride = iy * grid_dim[0];
            for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
                const int x_stride = ix;

                float x = grid_origin[0] + ix * grid_step_x[0] + iy * grid_step_y[0] + iz * grid_step_z[0];
                float y = grid_origin[1] + ix * grid_step_x[1] + iy * grid_step_y[1] + iz * grid_step_z[1];
                float z = grid_origin[2] + ix * grid_step_x[2] + iy * grid_step_y[2] + iz * grid_step_z[2];

                double psi = 0.0;
                for (size_t i = 0; i < num_gtos; ++i) {
                    float px	= gtos[i].x;
                    float py	= gtos[i].y;
                    float pz	= gtos[i].z;
                    float alpha	= gtos[i].alpha;
                    float coeff	= gtos[i].coeff;
                    int   pi	= gtos[i].i;
                    int   pj	= gtos[i].j;
                    int   pk	= gtos[i].k;

                    float dx = x - px;
                    float dy = y - py;
                    float dz = z - pz;
                    float d2 = dx * dx + dy * dy + dz * dz;
                    float fx = powf(dx, (float)pi);
                    float fy = powf(dy, (float)pj);
                    float fz = powf(dz, (float)pk);
                    float exp_term = (alpha == 0.0f) ? 1.0f : expf(-alpha * d2);
                    float powxyz = fx * fy * fz;
                    float prod = coeff * powxyz * exp_term;
                    psi += prod;
                }

                if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                    psi *= psi;
                }

                int index = x_stride + y_stride + z_stride;
                grid_data[index] += (float)psi;
            }
        }
    }
}

#if defined(__AVX512F__) && defined(__AVX512DQ__)

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_512(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_256i vix = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
    const md_256  vxh = md_mm256_fmadd_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step[0]), md_mm256_set1_ps(grid_origin[0]));
    const __m512  vx  = _mm512_insertf32x8(_mm512_castps256_ps512(vxh), vxh, 1);
    const int x_stride = grid_idx_min[0];

    __m512 vpsi[8][4] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const __m512  px = _mm512_set1_ps(gtos[i].x);
        const __m512  py = _mm512_set1_ps(gtos[i].y);
        const __m512  pz = _mm512_set1_ps(gtos[i].z);
        const __m512  pc = _mm512_set1_ps(gtos[i].coeff);
        const __m512  pa = _mm512_set1_ps(-gtos[i].alpha); // Negate alpha here
        const __m512i pi = _mm512_set1_epi32(gtos[i].i);
        const __m512i pj = _mm512_set1_epi32(gtos[i].j);
        const __m512i pk = _mm512_set1_epi32(gtos[i].k);

        for (int iz = 0; iz < 8; ++iz) {
            float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
            __m512 vz = _mm512_set1_ps(z);
            for (int iy = 0; iy < 4; ++iy) {
                float y[2] = {
                    grid_origin[1] + (grid_idx_min[1] + iy * 2 + 0) * grid_step[1],
                    grid_origin[1] + (grid_idx_min[1] + iy * 2 + 1) * grid_step[1],
                };
                __m512 vy = _mm512_insertf32x8(_mm512_set1_ps(y[0]), _mm256_set1_ps(y[1]), 1);

                __m512 dx = _mm512_sub_ps(vx, px);
                __m512 dy = _mm512_sub_ps(vy, py);
                __m512 dz = _mm512_sub_ps(vz, pz);
                __m512 d2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
                __m512 fx = md_mm512_fast_pow(dx, pi);
                __m512 fy = md_mm512_fast_pow(dy, pj);
                __m512 fz = md_mm512_fast_pow(dz, pk);
                __m512 ex = md_mm512_exp_ps(_mm512_mul_ps(pa, d2));
                __m512 prod = _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(pc, fx), _mm512_mul_ps(fy, fz)), ex);

                vpsi[iz][iy] = _mm512_add_ps(vpsi[iz][iy], prod);
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 4; ++iy) {
            int y_stride[2] = {
                (grid_idx_min[1] + iy * 2 + 0) * grid_dim[0],
                (grid_idx_min[1] + iy * 2 + 1) * grid_dim[0],
            };
            int index[2] = {
                x_stride + y_stride[0] + z_stride,
                x_stride + y_stride[1] + z_stride,
            };

            md_512 psi = vpsi[iz][iy];
            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi = _mm512_mul_ps(psi, psi);
            }

            md_256 tpsi[2] = {
                _mm512_castps512_ps256(psi),
                _mm512_extractf32x8_ps(psi, 1),
            };

            md_mm256_storeu_ps(grid_data + index[0], tpsi[0]);
            md_mm256_storeu_ps(grid_data + index[1], tpsi[1]);
        }
    }
}

#endif

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_256(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_256i vix  = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
    const int x_stride = grid_idx_min[0];
    const md_256   vx  = md_mm256_fmadd_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step[0]), md_mm256_set1_ps(grid_origin[0]));

    // Operate on local block to avoid cache-line contention across threads
    md_256 vpsi[8][8] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const md_256  px = md_mm256_set1_ps(gtos[i].x);
        const md_256  py = md_mm256_set1_ps(gtos[i].y);
        const md_256  pz = md_mm256_set1_ps(gtos[i].z);
        const md_256  pc = md_mm256_set1_ps(gtos[i].coeff);
        const md_256  pa = md_mm256_set1_ps(-gtos[i].alpha); // Negate alpha here
        const md_256i pi = md_mm256_set1_epi32(gtos[i].i);
        const md_256i pj = md_mm256_set1_epi32(gtos[i].j);
        const md_256i pk = md_mm256_set1_epi32(gtos[i].k);

        // The terms in x are the same for the whole block and the terms in z for a row, so they are formed outside
        // the loops that do not change them (the compiler may not reorder the products to do it itself)
        const md_256 dx  = md_mm256_sub_ps(vx, px);
        const md_256 dx2 = md_mm256_mul_ps(dx, dx);
        const md_256 cfx = md_mm256_mul_ps(pc, md_mm256_fast_pow(dx, pi));

        for (int iz = 0; iz < 8; ++iz) {
            float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
            const md_256 dz    = md_mm256_sub_ps(md_mm256_set1_ps(z), pz);
            const md_256 dxz2  = md_mm256_fmadd_ps(dz, dz, dx2);
            const md_256 cfxfz = md_mm256_mul_ps(cfx, md_mm256_fast_pow(dz, pk));

            for (int iy = 0; iy < 8; ++iy) {
                float y = grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1];
                md_256 dy = md_mm256_sub_ps(md_mm256_set1_ps(y), py);
                md_256 d2 = md_mm256_fmadd_ps(dy, dy, dxz2);
                md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));
                md_256 fy = md_mm256_fast_pow(dy, pj);

                vpsi[iz][iy] = md_mm256_fmadd_ps(md_mm256_mul_ps(cfxfz, fy), ex, vpsi[iz][iy]);
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index = x_stride + y_stride + z_stride;
            md_256 psi = vpsi[iz][iy];

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi = md_mm256_mul_ps(psi, psi);
            }

            md_mm256_storeu_ps(grid_data + index, md_mm256_add_ps(md_mm256_loadu_ps(grid_data + index), psi));
        }
    }
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_8x8x8_256(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_256i vix = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
    const int x_stride = grid_idx_min[0];

    const md_256 gsx[3] = {
        md_mm256_add_ps(md_mm256_set1_ps(grid_origin[0]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[0]))),
        md_mm256_add_ps(md_mm256_set1_ps(grid_origin[1]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[1]))),
        md_mm256_add_ps(md_mm256_set1_ps(grid_origin[2]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[2]))),
    };

    const md_256 gsz[3] = {
        md_mm256_set1_ps(grid_step_z[0]),
        md_mm256_set1_ps(grid_step_z[1]),
        md_mm256_set1_ps(grid_step_z[2]),
    };

    const md_256 gsy[3] = {
        md_mm256_set1_ps(grid_step_y[0]),
        md_mm256_set1_ps(grid_step_y[1]),
        md_mm256_set1_ps(grid_step_y[2]),
    };

    // Operate on local block to avoid cache-line contention across threads
    md_256 vpsi[8][8] = {0};

    for (size_t gto_idx = 0; gto_idx < num_gtos; ++gto_idx) {
        const float px = gtos[gto_idx].x;
        const float py = gtos[gto_idx].y;
        const float pz = gtos[gto_idx].z;
        const float pc = gtos[gto_idx].coeff;
        const float pa = -gtos[gto_idx].alpha; // Negate alpha here
        const int pi = gtos[gto_idx].i;
        const int pj = gtos[gto_idx].j;
        const int pk = gtos[gto_idx].k;
    
        for (int iz = 0; iz < 8; ++iz) {
            const md_256 tz = md_mm256_cvtepi32_ps(md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[2]), md_mm256_set1_epi32(iz)));

            const md_256 xz[3] = {
                md_mm256_fmadd_ps(tz, gsz[0], gsx[0]),
                md_mm256_fmadd_ps(tz, gsz[1], gsx[1]),
                md_mm256_fmadd_ps(tz, gsz[2], gsx[2]),
            };

            for (int iy = 0; iy < 8; ++iy) {
                const md_256 ty = md_mm256_cvtepi32_ps(md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[1]), md_mm256_set1_epi32(iy)));

                md_256 vx = md_mm256_fmadd_ps(ty, gsy[0], xz[0]);
                md_256 vy = md_mm256_fmadd_ps(ty, gsy[1], xz[1]);
                md_256 vz = md_mm256_fmadd_ps(ty, gsy[2], xz[2]);

                md_256 dx = md_mm256_sub_ps(vx, md_mm256_set1_ps(px));
                md_256 dy = md_mm256_sub_ps(vy, md_mm256_set1_ps(py));
                md_256 dz = md_mm256_sub_ps(vz, md_mm256_set1_ps(pz));
                md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
                md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(md_mm256_set1_ps(pa), d2));
                md_256 fx = md_mm256_fast_pow1(dx, pi);
                md_256 fy = md_mm256_fast_pow1(dy, pj);
                md_256 fz = md_mm256_fast_pow1(dz, pk);

                md_256 prod_a = md_mm256_mul_ps(md_mm256_set1_ps(pc), fx);
                md_256 prod_b = md_mm256_mul_ps(fy, fz);

                vpsi[iz][iy] = md_mm256_fmadd_ps(md_mm256_mul_ps(prod_a, prod_b), ex, vpsi[iz][iy]);
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index = x_stride + y_stride + z_stride;
            md_256 psi = vpsi[iz][iy];

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi = md_mm256_mul_ps(psi, psi);
            }

            md_mm256_storeu_ps(grid_data + index, md_mm256_add_ps(md_mm256_loadu_ps(grid_data + index), psi));
        }
    }
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_128(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const md_128i vix[2] = {
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 0), md_mm_set_epi32(3,2,1,0)),
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 4), md_mm_set_epi32(3,2,1,0)),
    };
    const md_128   vx[2] = {
        md_mm_fmadd_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step[0]), md_mm_set1_ps(grid_origin[0])),
        md_mm_fmadd_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step[0]), md_mm_set1_ps(grid_origin[0])),
    };
    const int x_stride[2] = {
        grid_idx_min[0] + 0,
        grid_idx_min[0] + 4,
    };

    // Operate on local block to avoid cache-line contention across threads
    md_128 vpsi[8][8][2] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const md_128  px = md_mm_set1_ps(gtos[i].x);
        const md_128  py = md_mm_set1_ps(gtos[i].y);
        const md_128  pz = md_mm_set1_ps(gtos[i].z);
        const md_128  pc = md_mm_set1_ps(gtos[i].coeff);
        const md_128  pa = md_mm_set1_ps(-gtos[i].alpha); // Negate alpha here
        const md_128i pi = md_mm_set1_epi32(gtos[i].i);
        const md_128i pj = md_mm_set1_epi32(gtos[i].j);
        const md_128i pk = md_mm_set1_epi32(gtos[i].k);

        // The terms in x are the same for the whole block and the terms in z for a row, so they are formed outside
        // the loops that do not change them (the compiler may not reorder the products to do it itself)
        const md_128 dx[2]  = { md_mm_sub_ps(vx[0], px), md_mm_sub_ps(vx[1], px) };
        const md_128 dx2[2] = { md_mm_mul_ps(dx[0], dx[0]), md_mm_mul_ps(dx[1], dx[1]) };
        const md_128 cfx[2] = { md_mm_mul_ps(pc, md_mm_fast_pow(dx[0], pi)), md_mm_mul_ps(pc, md_mm_fast_pow(dx[1], pi)) };

        for (int iz = 0; iz < 8; ++iz) {
            const md_128 dz = md_mm_sub_ps(md_mm_set1_ps(grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2]), pz);
            const md_128 fz = md_mm_fast_pow(dz, pk);
            const md_128 dxz2[2]  = { md_mm_fmadd_ps(dz, dz, dx2[0]), md_mm_fmadd_ps(dz, dz, dx2[1]) };
            const md_128 cfxfz[2] = { md_mm_mul_ps(cfx[0], fz), md_mm_mul_ps(cfx[1], fz) };

            for (int iy = 0; iy < 8; ++iy) {
                const md_128 dy = md_mm_sub_ps(md_mm_set1_ps(grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1]), py);
                const md_128 fy = md_mm_fast_pow(dy, pj);

                for (int ix = 0; ix < 2; ++ix) {
                    md_128 d2 = md_mm_fmadd_ps(dy, dy, dxz2[ix]);
                    md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));

                    vpsi[iz][iy][ix] = md_mm_fmadd_ps(md_mm_mul_ps(cfxfz[ix], fy), ex, vpsi[iz][iy][ix]);
                }
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index[2] = {
                x_stride[0] + y_stride + z_stride,
                x_stride[1] + y_stride + z_stride,
            };

            md_128 psi[2] = {
                vpsi[iz][iy][0],
                vpsi[iz][iy][1],
            };

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi[0] = md_mm_mul_ps(psi[0], psi[0]);
                psi[1] = md_mm_mul_ps(psi[1], psi[1]);
            }

            md_mm_storeu_ps(grid_data + index[0], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[0]), psi[0]));
            md_mm_storeu_ps(grid_data + index[1], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[1]), psi[1]));
        }
    }
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_8x8x8_128(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    const int x_stride[2] = {
        grid_idx_min[0] + 0,
        grid_idx_min[0] + 4,
    };

    const md_128i vix[2] = {
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0])    , md_mm_set_epi32(3,2,1,0)),
        md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 4), md_mm_set_epi32(3,2,1,0)),
    };

    const md_128 gsx[2][3] = {
        {
            md_mm_add_ps(md_mm_set1_ps(grid_origin[0]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[0]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[1]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[1]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[2]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[2]))),
        },
        {
            md_mm_add_ps(md_mm_set1_ps(grid_origin[0]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[0]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[1]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[1]))),
            md_mm_add_ps(md_mm_set1_ps(grid_origin[2]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[2]))),
        },
    };

    const md_128 gsz[3] = {
        md_mm_set1_ps(grid_step_z[0]),
        md_mm_set1_ps(grid_step_z[1]),
        md_mm_set1_ps(grid_step_z[2]),
    };

    const md_128 gsy[3] = {
        md_mm_set1_ps(grid_step_y[0]),
        md_mm_set1_ps(grid_step_y[1]),
        md_mm_set1_ps(grid_step_y[2]),
    };

    // Operate on local block to avoid cache-line contention across threads
    md_128 vpsi[8][8][2] = {0};

    for (size_t i = 0; i < num_gtos; ++i) {
        const md_128  px = md_mm_set1_ps(gtos[i].x);
        const md_128  py = md_mm_set1_ps(gtos[i].y);
        const md_128  pz = md_mm_set1_ps(gtos[i].z);
        const md_128  pc = md_mm_set1_ps(gtos[i].coeff);
        const md_128  pa = md_mm_set1_ps(-gtos[i].alpha); // Negate alpha here
        const md_128i pi = md_mm_set1_epi32(gtos[i].i);
        const md_128i pj = md_mm_set1_epi32(gtos[i].j);
        const md_128i pk = md_mm_set1_epi32(gtos[i].k);

        for (int iz = 0; iz < 8; ++iz) {
            const md_128 tz = md_mm_cvtepi32_ps(md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[2]), md_mm_set1_epi32(iz)));
            const md_128 xz[2][3] = {
                {
                    md_mm_fmadd_ps(tz, gsz[0], gsx[0][0]),
                    md_mm_fmadd_ps(tz, gsz[1], gsx[0][1]),
                    md_mm_fmadd_ps(tz, gsz[2], gsx[0][2]),
                },
                {
                    md_mm_fmadd_ps(tz, gsz[0], gsx[1][0]),
                    md_mm_fmadd_ps(tz, gsz[1], gsx[1][1]),
                    md_mm_fmadd_ps(tz, gsz[2], gsx[1][2]),
                },
            };
            for (int iy = 0; iy < 8; ++iy) {
                const md_128 ty = md_mm_cvtepi32_ps(md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[1]), md_mm_set1_epi32(iy)));

                const md_128 vx[2] = {
                    md_mm_fmadd_ps(ty, gsy[0], xz[0][0]),
                    md_mm_fmadd_ps(ty, gsy[0], xz[1][0]),
                };
                const md_128 vy[2] = {
                    md_mm_fmadd_ps(ty, gsy[1], xz[0][1]),
                    md_mm_fmadd_ps(ty, gsy[1], xz[1][1]),
                };
                const md_128 vz[2] = {
                    md_mm_fmadd_ps(ty, gsy[2], xz[0][2]),
                    md_mm_fmadd_ps(ty, gsy[2], xz[1][2]),
                };

                for (int ix = 0; ix < 2; ++ix) {
                    md_128 dx = md_mm_sub_ps(vx[ix], px);
                    md_128 dy = md_mm_sub_ps(vy[ix], py);
                    md_128 dz = md_mm_sub_ps(vz[ix], pz);
                    md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
                    md_128 fx = md_mm_fast_pow(dx, pi);
                    md_128 fy = md_mm_fast_pow(dy, pj);
                    md_128 fz = md_mm_fast_pow(dz, pk);
                    md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));
                    md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);

                    vpsi[iz][iy][ix] = md_mm_add_ps(vpsi[iz][iy][ix], prod);
                }
            }
        }
    }

    // Write result block to memory
    for (int iz = 0; iz < 8; ++iz) {
        int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
        for (int iy = 0; iy < 8; ++iy) {
            int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
            int index[2] = {
                x_stride[0] + y_stride + z_stride,
                x_stride[1] + y_stride + z_stride,
            };

            md_128 psi[2] = {
                vpsi[iz][iy][0],
                vpsi[iz][iy][1],
            };

            if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
                psi[0] = md_mm_mul_ps(psi[0], psi[0]);
                psi[1] = md_mm_mul_ps(psi[1], psi[1]);
            }

            md_mm_storeu_ps(grid_data + index[0], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[0]), psi[0]));
            md_mm_storeu_ps(grid_data + index[1], md_mm_add_ps(md_mm_loadu_ps(grid_data + index[1]), psi[1]));
        }
    }
}

void md_gto_grid_evaluate_sub(float* out_values, const md_grid_t* grid, const int grid_idx_off[3], const int grid_idx_len[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    ASSERT(grid);
    ASSERT(gtos);

    const int* grid_idx_min = grid_idx_off;
    const int  grid_idx_max[3] = {
        grid_idx_off[0] + grid_idx_len[0],
        grid_idx_off[1] + grid_idx_len[1],
        grid_idx_off[2] + grid_idx_len[2],
    };

    //printf("Number of pgtos in volume region: %zu\n", gto.count);

    bool ortho =
        (grid->orientation.elem[0][1] == 0 && grid->orientation.elem[0][2] == 0) &&
        (grid->orientation.elem[1][0] == 0 && grid->orientation.elem[1][2] == 0) &&
        (grid->orientation.elem[2][0] == 0 && grid->orientation.elem[2][1] == 0);

    vec3_t step_x = vec3_mul1(grid->orientation.col[0], grid->spacing.x);
    vec3_t step_y = vec3_mul1(grid->orientation.col[1], grid->spacing.y);
    vec3_t step_z = vec3_mul1(grid->orientation.col[2], grid->spacing.z);

    // There are specialized versions for evaluating 8x8x8 subgrids
    // 8x8x8 Is a good chunk size to operate on as it probably fits in L1 Cache together with the GTOs
    // Then we vectorize over the spatial domain rather than the GTOs to get better register occupation
    if (grid_idx_len[0] == 8 && grid_idx_len[1] == 8 && grid_idx_len[2] == 8) {

#if defined(__AVX512F__) && defined(__AVX512DQ__) || defined (__AVX2__) || defined(__aarch64__) || defined(_M_ARM64)
        // @TODO: Implement real AVX512 path
        if (ortho) {
            evaluate_grid_ortho_8x8x8_256(out_values, grid_idx_min, grid->dim, grid->origin.elem, grid->spacing.elem, gtos, num_gtos, mode);
        } else {
            evaluate_grid_8x8x8_256(out_values, grid_idx_min, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
        }
#elif defined(__SSE2__)
        if (ortho) {
            evaluate_grid_ortho_8x8x8_128(out_values, grid_idx_min, grid->dim, grid->origin.elem, grid->spacing.elem, gtos, num_gtos, mode);
        } else {
            evaluate_grid_8x8x8_128(out_values, grid_idx_min, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
        }
#else
        evaluate_grid_ref(out_values, grid_idx_min, grid_idx_max, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
#endif
    } else {
        // Slowpath
        evaluate_grid_ref(out_values, grid_idx_min, grid_idx_max, grid->dim, grid->origin.elem, step_x.elem, step_y.elem, step_z.elem, gtos, num_gtos, mode);
    }
}

void md_gto_grid_evaluate(float* out_values, const md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    ASSERT(grid);
    ASSERT(gtos);

    int idx_off[3] = {0};
    int idx_len[3] = {0};

    float scl[3] = {
        (grid->orientation.elem[0][0] + grid->orientation.elem[0][1] + grid->orientation.elem[0][2]) * grid->spacing.elem[0],
        (grid->orientation.elem[1][0] + grid->orientation.elem[1][1] + grid->orientation.elem[1][2]) * grid->spacing.elem[1],
        (grid->orientation.elem[2][0] + grid->orientation.elem[2][1] + grid->orientation.elem[2][2]) * grid->spacing.elem[2]
    };

    md_temp_scope_t temp = md_temp_begin();
    md_gto_t* sub_gtos = (md_gto_t*)md_temp_alloc(temp, sizeof(md_gto_t) * num_gtos);

    for (idx_off[2] = 0; idx_off[2] < grid->dim[2]; idx_off[2] += 8) {
        idx_len[2] = MIN(8, grid->dim[2] - idx_off[2]);
        for (idx_off[1] = 0; idx_off[1] < grid->dim[1]; idx_off[1] += 8) {
            idx_len[1] = MIN(8, grid->dim[1] - idx_off[1]);
            for (idx_off[0] = 0; idx_off[0] < grid->dim[0]; idx_off[0] += 8) {
                idx_len[0] = MIN(8, grid->dim[0] - idx_off[0]);

                float aabb_min[3] = {
                    grid->origin.elem[0] + idx_off[0] * scl[0],
                    grid->origin.elem[1] + idx_off[1] * scl[1],
                    grid->origin.elem[2] + idx_off[2] * scl[2],
                };
                float aabb_max[3] = {
                    grid->origin.elem[0] + (idx_off[0] + idx_len[0]) * scl[0],
                    grid->origin.elem[1] + (idx_off[1] + idx_len[1]) * scl[1],
                    grid->origin.elem[2] + (idx_off[2] + idx_len[2]) * scl[2],
                };

                size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);
                md_gto_grid_evaluate_sub(out_values, grid, idx_off, idx_len, sub_gtos, num_sub_gtos, mode);
            }
        }
    }

    md_temp_end(temp);
}

// Evaluate GTOs over a set of passed in packed XYZ coordinates with a bytestride
static void evaluate_gtos(float* out_psi, const float* in_xyz, size_t num_xyz, size_t xyz_stride, const md_gto_t* in_gto, size_t num_gtos, md_gto_eval_mode_t mode) {
    for (size_t j = 0; j < num_xyz; ++j) {
        const float* xyz = (const float*)((const char*)in_xyz + j * xyz_stride);
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        double psi = 0.0;
        for (size_t i = 0; i < num_gtos; ++i) {
            double cutoff	= in_gto[i].cutoff;
            double rx		= x - in_gto[i].x;
            double ry		= y - in_gto[i].y;
            double rz		= z - in_gto[i].z;
            double r2		= rx * rx + ry * ry + rz * rz;
            if (r2 > cutoff * cutoff) {
                continue;
            }

            double alpha	= in_gto[i].alpha;
            double coeff	= in_gto[i].coeff;
            int   pi		= in_gto[i].i;
            int   pj		= in_gto[i].j;
            int   pk		= in_gto[i].k;

            double fx = pow(rx, pi);
            double fy = pow(ry, pj);
            double fz = pow(rz, pk);
            double powxyz = fx * fy * fz;
            double exp_term = alpha == 0 ? 1.0 : exp(-alpha * r2);

            double prod = coeff * powxyz * exp_term;
            psi += prod;
        }

        if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
            psi = psi * psi;
        }

        out_psi[j] = (float)psi;
    }
}

void md_gto_xyz_evaluate(float* out_psi, const float* in_xyz, size_t num_xyz, size_t stride, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
    if (!out_psi) {
        MD_LOG_ERROR("out_psi array is NULL!");
        return;
    }
    if (!in_xyz) {
        MD_LOG_ERROR("in_xyz base pointer is NULL!");
        return;
    }
    if (!gtos) {
        MD_LOG_ERROR("gtos is NULL!");
        return;
    }
    if (stride != 0 && stride < sizeof(float) * 3) {
        MD_LOG_ERROR("Invalid xyz stride: expected value >= 12 Bytes, got %zu", stride);
        return;
    }

    stride = (stride == 0) ? sizeof(float) * 3 : stride;
    evaluate_gtos(out_psi, in_xyz, num_xyz, stride, gtos, num_gtos, mode);
}

static inline double eval_G(double d, double C, int l, double neg_alpha) {
    return C * fast_pow(d, l) * exp(neg_alpha * d * d);
}

static inline void eval_G_and_G_prime(double* out_G, double* out_G_prime, double d, double C, int l, double neg_alpha) {
    double exp_term = exp(neg_alpha * d * d);
    *out_G		 = C * fast_pow(d, l) * exp_term;
    *out_G_prime = C * fast_pow(d, l-1) * (l + 2 * neg_alpha * d * d) * exp_term;
}

#define PRINT_RESULT 0

static double compute_distance_cutoff(double cutoff_value, int i, int j, int k, int l, double coeff, double alpha) {
    double d = 0.0;

    const double neg_alpha = -alpha;

    // Bake into single constant C
    const double C = fabs(coeff * sqrt((fast_pow(i,i) * fast_pow(j,j) * fast_pow(k,k)) / fast_pow(l,l)));

    // Compute maxima
    const double d_maxima = sqrt(l / (2.0 * fabs(neg_alpha)));

    // Check the contribution at the maxima
    const double y_max = eval_G(d_maxima, C, l, neg_alpha);
    if (y_max < cutoff_value) {
        d = 0.0;
        goto done;
    }

    // If we have an S-type orbital (l == 0) the expression collapses into an expression we can invert and evaluate
    if (l == 0) {
        double y = cutoff_value;
        double a = fabs(coeff) / y;
        double la = log(a);
        d = sqrt(fabs(la) / fabs(neg_alpha));
        goto done;
    }

    // If we end up here we need to perform a numerical search for a d value where the value G(d) < cutoff_value
    // We do not want to overestimate d which will result in a too large radius of influence for the PGTO.
    // And will directly negatively impact performance when evaluating on the grid.
    // Therefore we want to find d where G(d) < cutoff_value but within a tolerance of cutoff_value

    // Search parameters
    const double d_min = d_maxima + 0.001;
    const double d_max = d_maxima + 100.0;
    const double y_tol = cutoff_value * 0.001;
    const double d_tol = 1.0e-9;

    // Initial guess
    // This is the analytical solution for d^2/dx^2 G(d) = 0
    // Which should give us a value which corresponds to the point where we have the maximum negative slope
    d = 0.5 * sqrt(sqrt(neg_alpha*neg_alpha * (8*l + 1)) / (neg_alpha*neg_alpha) + (2*l + 1) / fabs(neg_alpha));

    // Newton-Rhapson iterative search
    for (int iter = 0; iter < 100; ++iter) {
        double y, yp;
        eval_G_and_G_prime(&y, &yp, d, C, l, neg_alpha);

        // Shift function so it intersects the x axis at the point we seek (with a bias towards values less than cutoff value)
        y = y - cutoff_value + y_tol;

        //printf ("d: %.10f, y: %.10f, yp: %.10f\n", d, y, yp);

        if (y < 0 && fabs(y) < y_tol) {
            //printf ("y tolerance met after %i iterations\n", iter);
            break;
        }

        if (fabs(yp) < DBL_EPSILON) {
            //printf ("Denominator is too small!\n");
            break;
        }

        double dn = d - y / yp;
        dn = CLAMP(dn, d_min, d_max);

        if (fabs(dn - d) < d_tol) {
            //printf ("d tolerance met after %i iterations\n", iter);
            break;
        }

        d = dn;
    }

done:
#if PRINT_RESULT
    if (d > 0.0) {
        printf("Cutoff dist and value: %15.5f, %15.12f\n", d, eval_G(d, C, l, neg_alpha));
    }
#endif
    return d;
}

double md_gto_compute_radius_of_influence(int i, int j, int k, double coeff, double alpha, double cutoff) {
    int l = i + j + k;
    return compute_distance_cutoff(cutoff, i, j, k, l, coeff, alpha);
}

size_t md_gto_cutoff_compute_and_filter(md_gto_t* gtos, size_t count, double value) {
    if (value == 0) {
        for (size_t i = 0; i < count; ++i) {
            gtos[i].cutoff = FLT_MAX;
        }
    } else {
        for (size_t i = 0; i < count;) {
            gtos[i].cutoff = (float)compute_distance_cutoff(value, gtos[i].i, gtos[i].j, gtos[i].k, gtos[i].l, gtos[i].coeff, gtos[i].alpha);
            if (gtos[i].cutoff == 0.0f) {
                gtos[i] = gtos[--count];
            } else {
                ++i;
            }
        }
    }
    return count;
}

size_t md_gto_aabb_test(md_gto_t* out_gtos, const float aabb_min[3], const float aabb_max[3], const md_gto_t* in_gtos, size_t in_num_gtos) {
    // Extract a subset of gtos that overlap with the evaluated subportion of the grid
    // @TODO: This can be vectorized, Let us pray to the compiler gods for now
    size_t num_gtos = 0;
    for (size_t i = 0; i < in_num_gtos; ++i) {
        float x  = in_gtos[i].x;
        float y  = in_gtos[i].y;
        float z  = in_gtos[i].z;
        float cutoff = in_gtos[i].cutoff;

        float cx = CLAMP(x, aabb_min[0], aabb_max[0]);
        float cy = CLAMP(y, aabb_min[1], aabb_max[1]);
        float cz = CLAMP(z, aabb_min[2], aabb_max[2]);

        float dx = x - cx;
        float dy = y - cy;
        float dz = z - cz;

        float d2 = dx * dx + dy * dy + dz * dz;

        if (d2 > cutoff * cutoff) {
            continue;
        }
        out_gtos[num_gtos++] = in_gtos[i];
    }
    return num_gtos;
}

// ---------------------------------------------------------------------------
// Basis from the attribute table
// ---------------------------------------------------------------------------

// Finds one of our own columns and hands back its storage. Every check is explicit because the
// table is open - anyone may publish under any path - so a basis/ path with the wrong type or rank
// is a mistake to refuse rather than to reinterpret.
static const void* gto_attr_column(size_t* out_count, const md_attributes_t* attributes, str_t path, md_attribute_type_t type) {
    const md_attribute_t* attr = md_attributes_find(attributes, path);
    const void* data = md_attribute_view(attr, type, 1, 1);
    if (attr && !data) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' is not the plain column a basis expects", STR_ARG(path));
    }
    if (data && out_count) *out_count = attr->format.shape[0];
    return data;
}

bool md_gto_basis_extract_attributes(md_gto_basis_t* out, const md_attributes_t* attributes, md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(alloc);

    if (!attributes) {
        return false;
    }

    MEMSET(out, 0, sizeof(*out));

    size_t num_atom_idx = 0, num_offset = 0, num_count = 0, num_angl = 0;
    const uint32_t* atom_idx = (const uint32_t*)gto_attr_column(&num_atom_idx, attributes, STR_LIT("basis/shell/atom_index"),        MD_ATTRIBUTE_TYPE_U32);
    const uint32_t* offset   = (const uint32_t*)gto_attr_column(&num_offset,   attributes, STR_LIT("basis/shell/primitive_offset"),  MD_ATTRIBUTE_TYPE_U32);
    const uint32_t* count    = (const uint32_t*)gto_attr_column(&num_count,    attributes, STR_LIT("basis/shell/primitive_count"),   MD_ATTRIBUTE_TYPE_U32);
    const uint32_t* angl     = (const uint32_t*)gto_attr_column(&num_angl,     attributes, STR_LIT("basis/shell/angular_momentum"),  MD_ATTRIBUTE_TYPE_U32);

    size_t num_alpha = 0, num_coeff = 0;
    const float* alpha = (const float*)gto_attr_column(&num_alpha, attributes, STR_LIT("basis/primitive/exponent"),    MD_ATTRIBUTE_TYPE_F32);
    const float* coeff = (const float*)gto_attr_column(&num_coeff, attributes, STR_LIT("basis/primitive/coefficient"), MD_ATTRIBUTE_TYPE_F32);

    if (!atom_idx || !offset || !count || !angl || !alpha || !coeff) {
        return false;
    }

    // These four are COLUMNS of one table, so they must be equal length - that is a property of
    // this particular group and not a general rule about siblings. Anything else is a table written
    // in pieces by more than one producer, and guessing which piece is current is not something to
    // attempt.
    const size_t num_shells = num_atom_idx;
    if (num_offset != num_shells || num_count != num_shells || num_angl != num_shells) {
        MD_LOG_ERROR("basis/shell columns disagree on length");
        return false;
    }
    const size_t num_primitives = num_alpha;
    if (num_coeff != num_primitives) {
        MD_LOG_ERROR("basis/primitive columns disagree on length");
        return false;
    }
    if (num_shells == 0 || num_primitives == 0) {
        return false;
    }

    // The shell list is stored as four columns rather than as a record, so that every value in the
    // table is a self describing quantity and a new field is additive. The cost is this one
    // interleave, paid once when a consumer builds a basis and not per evaluation.
    md_array_resize(out->shells, num_shells, alloc);
    for (size_t i = 0; i < num_shells; ++i) {
        if ((size_t)offset[i] + (size_t)count[i] > num_primitives) {
            MD_LOG_ERROR("Shell %zu spans primitives [%u,%u) beyond the %zu published", i, offset[i], offset[i] + count[i], num_primitives);
            md_array_free(out->shells, alloc);
            MEMSET(out, 0, sizeof(*out));
            return false;
        }
        out->shells[i] = (md_gto_shell_t){
            .atom_idx         = atom_idx[i],
            .primitive_offset = offset[i],
            .num_primitives   = count[i],
            .l                = angl[i],
        };
    }

    md_array_resize(out->alpha, num_primitives, alloc);
    md_array_resize(out->coeff, num_primitives, alloc);
    MEMCPY(out->alpha, alpha, num_primitives * sizeof(float));
    MEMCPY(out->coeff, coeff, num_primitives * sizeof(float));

    out->num_shells     = (uint32_t)num_shells;
    out->num_primitives = (uint32_t)num_primitives;
    return true;
}

void md_gto_basis_free(md_gto_basis_t* basis, md_allocator_i* alloc) {
    if (!basis) return;
    ASSERT(alloc);
    md_array_free(basis->shells, alloc);
    md_array_free(basis->alpha,  alloc);
    md_array_free(basis->coeff,  alloc);
    MEMSET(basis, 0, sizeof(*basis));
}
