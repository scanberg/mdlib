#include <md_gto.h>
#include <md_util.h>

#include <core/md_platform.h>
#include <core/md_gpu_util.h>

#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <stddef.h>

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

static inline void index_to_world_matrix(float out_mat[4][4], const md_grid_t* grid) {
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
    out_mat[3][0] += 0.5f * (out_mat[0][0] + out_mat[1][0] + out_mat[2][0]);
    out_mat[3][1] += 0.5f * (out_mat[0][1] + out_mat[1][1] + out_mat[2][1]);
    out_mat[3][2] += 0.5f * (out_mat[0][2] + out_mat[1][2] + out_mat[2][2]);
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

// Count the number of CGTOs and PGTOs that will result from expanding a basis.
static void gto_basis_count(uint32_t* out_num_cgtos, uint32_t* out_num_pgtos,
    const md_gto_basis_t* basis)
{
    uint32_t nc = 0, np = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        int l      = (int)basis->shells[si].l;
        int nsph   = 2 * l + 1;
        int nprims = (int)basis->shells[si].num_primitives;
        nc += (uint32_t)nsph;
        for (int isph = 0; isph < nsph; isph++) {
            np += (uint32_t)(gto_sph_num_factors(l, isph) * nprims);
        }
    }
    *out_num_cgtos = nc;
    *out_num_pgtos = np;
}

static uint32_t gto_basis_num_atoms(const md_gto_basis_t* basis) {
    uint32_t max_atom_idx = 0;
    for (uint32_t i = 0; i < basis->num_shells; ++i) {
        max_atom_idx = MAX(max_atom_idx, basis->shells[i].atom_idx);
    }
    return basis->num_shells ? (max_atom_idx + 1) : 0;
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
        int nsph   = 2 * l + 1;
        int nprims = (int)shell->num_primitives;
        uint32_t   prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);

        for (int isph = 0; isph < nsph; isph++) {
            int            ncomp  = gto_sph_num_factors(l, isph);
            const double*  fcarts = gto_sph_factors(l, isph);
            const uint8_t* sidx   = gto_sph_indices(l, isph);

            int lx[8], ly[8], lz[8];
            for (int ic = 0; ic < ncomp; ic++) {
                lx[ic] = lmn[sidx[ic]][0];
                ly[ic] = lmn[sidx[ic]][1];
                lz[ic] = lmn[sidx[ic]][2];
            }

            uint32_t pi_beg = pi;
            double max_r = 0.0;
            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                for (int ic = 0; ic < ncomp; ic++) {
                    float coeff_val = (float)(coef1 * fcarts[ic]);
                    double radius = md_gto_compute_radius_of_influence(lx[ic], ly[ic], lz[ic], (double)coeff_val, (double)alpha, cutoff);
                    max_r = MAX(max_r, radius);

                    PGTO pgto = {
                        .coeff = coeff_val,
                        .alpha = alpha,
                        .ijkl = md_gto_pack_ijkl(lx[ic], ly[ic], lz[ic], l),
                    };
                    out_pgto[pi++] = pgto;
                }
            }

            out_cgto_atom_idx[ci] = shell->atom_idx;
            out_cgto_r[ci] = (float)max_r;
            out_cgto_off_len[2 * ci + 0] = pi_beg;
            out_cgto_off_len[2 * ci + 1] = pi - pi_beg;
            ci++;
        }
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
        int nsph   = 2 * l + 1;
        int nprims = (int)shell->num_primitives;
        uint32_t   prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);
        const float* ap = (const float*)((const uint8_t*)atom_xyz + shell->atom_idx * stride);
        float ax = ap[0];
        float ay = ap[1];
        float az = ap[2];

        for (int isph = 0; isph < nsph; isph++) {
            int           ncomp   = gto_sph_num_factors(l, isph);
            const double* fcarts  = gto_sph_factors(l, isph);
            const uint8_t* sidx   = gto_sph_indices(l, isph);

            int lx[8], ly[8], lz[8];
            for (int ic = 0; ic < ncomp; ic++) {
                lx[ic] = lmn[sidx[ic]][0];
                ly[ic] = lmn[sidx[ic]][1];
                lz[ic] = lmn[sidx[ic]][2];
            }

            uint32_t pi_beg = pi;
            
            double max_r = 0.0;
            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                for (int ic = 0; ic < ncomp; ic++) {
                    float coeff_val = (float)(coef1 * fcarts[ic]);
                    double radius = md_gto_compute_radius_of_influence(lx[ic], ly[ic], lz[ic], (double)coeff_val, (double)alpha, cutoff);
                    max_r = MAX(max_r, radius);

                    PGTO pgto = {
                        .coeff = coeff_val,
                        .alpha = alpha,
                        .ijkl = md_gto_pack_ijkl(lx[ic], ly[ic], lz[ic], l),
                    };
                    out_pgto[pi++] = pgto;
                }
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

size_t md_gto_expand_with_mo(md_gto_t* out, const md_gto_basis_t* basis,
    const float* atom_xyz, size_t atom_xyz_stride, const double* mo_coeffs, double cutoff)
{
    ASSERT(out);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(mo_coeffs);

    const size_t stride = atom_xyz_stride == 0 ? sizeof(float) * 3 : atom_xyz_stride;

    size_t num_gtos = 0;
    uint32_t cgto_idx = 0;
    for (uint32_t si = 0; si < basis->num_shells; si++) {
        const md_gto_shell_t* shell = &basis->shells[si];
        int l      = (int)shell->l;
        int nsph   = 2 * l + 1;
        int nprims = (int)shell->num_primitives;
        uint32_t prim_base = shell->primitive_offset;
        const gto_lmn_t* lmn = gto_cart_lmn(l);
        const float* ap = (const float*)((const uint8_t*)atom_xyz + shell->atom_idx * stride);
        float ax = ap[0];
        float ay = ap[1];
        float az = ap[2];

        for (int isph = 0; isph < nsph; isph++) {
            double mo_coeff = mo_coeffs[cgto_idx++];
            int           ncomp  = gto_sph_num_factors(l, isph);
            const double* fcarts = gto_sph_factors(l, isph);
            const uint8_t* sidx  = gto_sph_indices(l, isph);

            for (int ip = 0; ip < nprims; ip++) {
                float alpha = basis->alpha[prim_base + ip];
                float coef1 = basis->coeff[prim_base + ip];
                for (int ic = 0; ic < ncomp; ic++) {
                    out[num_gtos++] = (md_gto_t){
                        .x      = ax,
                        .y      = ay,
                        .z      = az,
                        .coeff  = (float)(coef1 * fcarts[ic] * mo_coeff),
                        .alpha  = alpha,
                        .cutoff = FLT_MAX,
                        .i      = lmn[sidx[ic]][0],
                        .j      = lmn[sidx[ic]][1],
                        .k      = lmn[sidx[ic]][2],
                        .l      = (uint8_t)l,
                    };
                }
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

    world_to_model_matrix(world_to_model, grid);
    index_to_world_matrix(index_to_world, grid);

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
    index_to_world_matrix(ub_data.index_to_world.elem, grid);
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

    md_temp_t temp = md_temp_begin();
    size_t max_gtos  = md_gto_pgto_count(basis);
    md_gto_t* gtos   = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * max_gtos);
    size_t num_gtos  = md_gto_expand_with_mo(gtos, basis, atom_xyz, atom_xyz_stride, mo_coeffs, cutoff);

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

    md_temp_t temp = md_temp_begin();
    size_t max_gtos  = md_gto_pgto_count(basis);

    // Flat GTO buffer for all MOs concatenated
    md_gto_t* gtos        = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * max_gtos * num_mos);
    uint32_t* orb_offsets = (uint32_t*)md_temp_push(sizeof(uint32_t) * (num_mos + 1));
    float*    orb_scaling = (float*)   md_temp_push(sizeof(float)    *  num_mos);

    size_t total_gtos = 0;
    size_t num_orbs   = 0;
    orb_offsets[0]    = 0;

    for (size_t i = 0; i < num_mos; i++) {
        if (!mo_coeffs[i]) continue;
        size_t n = md_gto_expand_with_mo(gtos + total_gtos, basis, atom_xyz, atom_xyz_stride, mo_coeffs[i], cutoff);
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
    ASSERT(grid);
    ASSERT(basis);
    ASSERT(atom_xyz);
    ASSERT(density_matrix);

    uint32_t num_cgtos, num_pgtos;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);

    md_temp_t temp = md_temp_begin();
    float*    cgto_xyz     = (float*)   md_temp_push(sizeof(float)    * 3 * num_cgtos);
    float*    cgto_r       = (float*)   md_temp_push(sizeof(float)    * 1 * num_cgtos);
    uint32_t* cgto_off_len = (uint32_t*)md_temp_push(sizeof(uint32_t) * num_cgtos * 2);
    PGTO*     pgto         = (PGTO*)    md_temp_push(sizeof(PGTO)     * num_pgtos);
    size_t    tri_len      = density_matrix_upper_tri_size(num_cgtos);
    float*    upper_tri    = (float*)   md_temp_push(sizeof(float)    * tri_len);

    gto_expand_basis(cgto_xyz, cgto_r, cgto_off_len, pgto, basis, atom_xyz, atom_xyz_stride, 1.0e-6);
    density_matrix_upper_tri_extract_float(upper_tri, density_matrix, num_cgtos);

    // Recombine into float4 for the GL path, which keeps its own xyzr SSBO layout.
    vec4_t* cgto_xyzr = (vec4_t*)md_temp_push(sizeof(vec4_t) * num_cgtos);
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

#endif

#if MD_ENABLE_GPU

#include <core/md_gpu.h>
#include <gto_gpu_shaders.inl>
#include <gto_gpu_shaders_reflection.inl>

static md_gpu_compute_pipeline_t gto_pip_density   = NULL;
static md_gpu_compute_pipeline_t gto_pip_mo        = NULL;
#if defined(MD_GPU_BACKEND_VULKAN) || defined(MD_GPU_BACKEND_METAL)
static md_gpu_compute_pipeline_t gto_pip_mo_root   = NULL;
#endif

static md_gpu_compute_pipeline_t ensure_compute_pipeline(md_gpu_device_t device, md_gpu_compute_pipeline_t* pipeline, const void* blob_start, size_t blob_size, const char* name, uint32_t wg_x, uint32_t wg_y, uint32_t wg_z) {
    if (*pipeline == NULL) {
        md_gpu_compute_pipeline_desc_t desc = {
            .shader_bytes     = blob_start,
            .shader_byte_size = blob_size,
            .threadgroup_size = { wg_x, wg_y, wg_z },
        };
        *pipeline = md_gpu_compute_pipeline_create(device, &desc);
        if (*pipeline == NULL) {
            MD_LOG_ERROR("Failed to create compute pipeline: %s", name);
        }
    }
    return *pipeline;
}

void md_gto_gpu_initialize(md_gpu_device_t device) {
    if (!gto_pip_density) {
        gto_pip_density = ensure_compute_pipeline(device, &gto_pip_density, gto_eval_gto_density_start, gto_eval_gto_density_size(), "GTO density", 8, 8, 8);
        if (!gto_pip_density) {
            MD_LOG_ERROR("Failed to create GTO density compute pipeline");
        }
    }

    if (!gto_pip_mo) {
        gto_pip_mo = ensure_compute_pipeline(device, &gto_pip_mo, gto_eval_gto_mo_start, gto_eval_gto_mo_size(), "GTO MO", 8, 8, 8);
        if (!gto_pip_mo) {
            MD_LOG_ERROR("Failed to create GTO MO compute pipeline");
        }
    }

#if defined(MD_GPU_BACKEND_VULKAN) || defined(MD_GPU_BACKEND_METAL)
    if (!gto_pip_mo_root) {
        gto_pip_mo_root = ensure_compute_pipeline(device, &gto_pip_mo_root, gto_eval_gto_mo_root_start, gto_eval_gto_mo_root_size(), "GTO MO root", 8, 8, 8);
        if (!gto_pip_mo_root) {
            MD_LOG_ERROR("Failed to create GTO MO root compute pipeline");
        }
    }
#endif
}

void md_gto_gpu_shutdown(void) {
    if (gto_pip_density) { md_gpu_compute_pipeline_destroy(gto_pip_density); gto_pip_density = NULL; }
    if (gto_pip_mo)      { md_gpu_compute_pipeline_destroy(gto_pip_mo);      gto_pip_mo      = NULL; }
#if defined(MD_GPU_BACKEND_VULKAN) || defined(MD_GPU_BACKEND_METAL)
    if (gto_pip_mo_root) { md_gpu_compute_pipeline_destroy(gto_pip_mo_root); gto_pip_mo_root = NULL; }
#endif
}

// Internal layout descriptor for the basis GPU buffer.
// Coefficient and atom coordinate data live in separate caller-owned buffers.
//
// Regions:
//   cgto_atom_idx: uint   × num_cgtos
//   cgto_r:        float  × num_cgtos
//   cgto_off_len:  uint2  × num_cgtos
//   pgto:          PGTO   × num_pgtos
typedef struct {
    uint32_t num_cgtos;
    uint32_t num_pgtos;
    uint32_t num_atoms;
    uint32_t off_cgto_atom_idx; // uint  × num_cgtos
    uint32_t off_cgto_r;        // float × num_cgtos
    uint32_t off_cgto_off_len;  // uint2 × num_cgtos
    uint32_t off_pgto;          // PGTO  × num_pgtos
    uint64_t total_size;
} md_gto_basis_layout_t;

static md_gto_basis_layout_t gto_basis_layout_compute(uint32_t num_cgtos, uint32_t num_pgtos, uint32_t num_atoms) {
    md_gto_basis_layout_t L = {0};
    L.num_cgtos = num_cgtos;
    L.num_pgtos = num_pgtos;
    L.num_atoms = num_atoms;

    L.off_cgto_atom_idx = 0;
    uint32_t end_cgto_atom_idx = L.off_cgto_atom_idx + (uint32_t)(sizeof(uint32_t) * num_cgtos);

    L.off_cgto_r       = (uint32_t)ALIGN_TO(end_cgto_atom_idx, 256);
    uint32_t end_cgto_r       = L.off_cgto_r       + (uint32_t)(sizeof(float)    * 1 * num_cgtos);

    L.off_cgto_off_len = (uint32_t)ALIGN_TO(end_cgto_r,       256);
    uint32_t end_cgto_off_len = L.off_cgto_off_len + (uint32_t)(sizeof(uint32_t) * 2 * num_cgtos);

    L.off_pgto         = (uint32_t)ALIGN_TO(end_cgto_off_len, 256);
    uint32_t end_pgto  = L.off_pgto + (uint32_t)(sizeof(PGTO) * num_pgtos);

    L.total_size       = (uint64_t)ALIGN_TO(end_pgto, 256);
    return L;
}

// ---------------------------------------------------------------------------
// md_gto_gpu_basis_t  —  device-local basis buffer (atom_idx, r, off_len, pgto)
// ---------------------------------------------------------------------------

typedef struct md_gto_gpu_basis {
    md_gpu_device_t       device;
    md_gpu_buffer_t       buffer;
    md_gto_basis_layout_t layout;
} md_gto_gpu_basis;

// ---------------------------------------------------------------------------
// Internal staging helpers (local bump, not exposed)
// ---------------------------------------------------------------------------

static bool gto_local_staging_upload_at(md_gpu_device_t device, md_gpu_buffer_t dst, size_t dst_offset,
                                        const void* src, size_t size) {
    md_gpu_bump_alloc_t bump = {0};
    bump.device = device;
    if (!md_gpu_bump_ensure(&bump, size)) return false;
    md_gpu_alloc_t a = md_gpu_bump_push(&bump, size);
    MEMCPY(a.cpu, src, size);

    md_gpu_pass_buffer_usage_t buffers[2] = {
        { .buffer = bump.buffer, .usage = MD_GPU_PASS_RESOURCE_TRANSFER_SRC },
        { .buffer = dst,         .usage = MD_GPU_PASS_RESOURCE_TRANSFER_DST },
    };
    md_gpu_pass_desc_t pass_desc = {
        .label = "GTO staging upload",
        .buffers = buffers,
        .buffer_count = 2,
    };
    md_gpu_pass_t pass = md_gpu_pass_begin(device, &pass_desc);
    if (!pass) {
        md_gpu_bump_free(&bump);
        return false;
    }
    md_gpu_pass_copy_buffer(pass, bump.buffer, dst, size, a.offset, dst_offset);
    md_gpu_pass_id_t pass_id = md_gpu_pass_end(pass);
    md_gpu_pass_wait(device, pass_id);
    md_gpu_bump_free(&bump);
    return pass_id.value != 0;
}

md_gto_gpu_basis_t md_gto_gpu_basis_create(md_gpu_device_t device, const md_gto_gpu_basis_desc_t* desc) {
    ASSERT(device);
    ASSERT(desc && desc->basis);
    const md_gto_basis_t* basis = desc->basis;

    uint32_t num_cgtos = 0, num_pgtos = 0;
    gto_basis_count(&num_cgtos, &num_pgtos, basis);
    uint32_t num_atoms = gto_basis_num_atoms(basis);

    md_gto_basis_layout_t layout = gto_basis_layout_compute(num_cgtos, num_pgtos, num_atoms);
    md_gpu_buffer_desc_t buf_desc = { .size = layout.total_size, .flags = MD_GPU_BUFFER_NONE };
    md_gpu_buffer_t buf = md_gpu_buffer_create(device, &buf_desc);
    if (!buf) {
        MD_LOG_ERROR("md_gto_gpu_basis_create: failed to allocate device-local buffer (%zu bytes)", (size_t)layout.total_size);
        return NULL;
    }
    md_gto_gpu_basis* gb = (md_gto_gpu_basis*)calloc(1, sizeof(md_gto_gpu_basis));
    if (!gb) { md_gpu_buffer_destroy(buf); return NULL; }
    gb->device = device;
    gb->buffer = buf;
    gb->layout = layout;

    // Populate all basis regions immediately.
    const md_gto_basis_layout_t* L = &gb->layout;
    bool uma = (md_gpu_buffer_flags(gb->buffer) & MD_GPU_BUFFER_CPU_VISIBLE) != 0;
    bool success = false;

    md_temp_t temp = md_temp_begin();

    uint32_t* cgto_atom_idx = (uint32_t*)md_temp_push(sizeof(uint32_t) * L->num_cgtos);
    float*    cgto_r       = (float*)   md_temp_push(sizeof(float)    * 1 * L->num_cgtos);
    uint32_t* cgto_off_len = (uint32_t*)md_temp_push(sizeof(uint32_t) * 2 * L->num_cgtos);
    PGTO*     pgto         = (PGTO*)    md_temp_push(sizeof(PGTO)         * L->num_pgtos);

    gto_expand_basis_gpu_meta(cgto_atom_idx, cgto_r, cgto_off_len, pgto, basis, desc->cutoff);

    size_t sz_atom_idx = L->off_cgto_r       - L->off_cgto_atom_idx;
    size_t sz_r       = L->off_cgto_off_len - L->off_cgto_r;
    size_t sz_off_len = L->off_pgto         - L->off_cgto_off_len;
    size_t sz_pgto    = (size_t)L->total_size - L->off_pgto;

    if (uma) {
        uint8_t* dst = (uint8_t*)md_gpu_buffer_cpu_ptr(gb->buffer);
        MEMCPY(dst + L->off_cgto_atom_idx, cgto_atom_idx, sz_atom_idx);
        MEMCPY(dst + L->off_cgto_r,       cgto_r,       sz_r);
        MEMCPY(dst + L->off_cgto_off_len, cgto_off_len, sz_off_len);
        MEMCPY(dst + L->off_pgto,         pgto,         sz_pgto);
    } else {
        // Pack into a single temp staging buffer and upload in one command.
        size_t total = (size_t)L->total_size;
        uint8_t* tmp = (uint8_t*)md_temp_push(total);
        MEMCPY(tmp + L->off_cgto_atom_idx, cgto_atom_idx, sz_atom_idx);
        MEMCPY(tmp + L->off_cgto_r,       cgto_r,       sz_r);
        MEMCPY(tmp + L->off_cgto_off_len, cgto_off_len, sz_off_len);
        MEMCPY(tmp + L->off_pgto,         pgto,         sz_pgto);
        if (!gto_local_staging_upload_at(gb->device, gb->buffer, 0, tmp, total)) {
            MD_LOG_ERROR("md_gto_gpu_basis_create: staging upload failed");
            goto done;
        }
    }

    success = true;

done:
    md_temp_end(temp);
    if (!success) {
        md_gpu_buffer_destroy(buf);
        free(gb);
        return NULL;
    }

    return gb;
}

void md_gto_gpu_basis_destroy(md_gto_gpu_basis_t gb) {
    if (!gb) return;
    md_gpu_buffer_destroy(gb->buffer);
    free(gb);
}

md_gpu_buffer_t md_gto_gpu_basis_buffer(md_gto_gpu_basis_t gb) {
    if (!gb) return NULL;
    return gb->buffer;
}

uint32_t md_gto_gpu_basis_num_cgtos(md_gto_gpu_basis_t gb) {
    if (!gb) return 0;
    return gb->layout.num_cgtos;
}

uint32_t md_gto_gpu_basis_num_pgtos(md_gto_gpu_basis_t gb) {
    if (!gb) return 0;
    return gb->layout.num_pgtos;
}

uint32_t md_gto_gpu_basis_num_atoms(md_gto_gpu_basis_t gb) {
    if (!gb) return 0;
    return gb->layout.num_atoms;
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

void md_gto_gpu_coeff_upload_density_pass(md_gpu_pass_t pass, md_gpu_buffer_t coeff_buf,
    md_gpu_buffer_t src_buf, size_t src_offset, size_t num_cgtos) {
    if (!pass || !coeff_buf || !src_buf) return;
    size_t sz = md_gto_gpu_coeff_size_density(num_cgtos);
    md_gpu_pass_copy_buffer(pass, src_buf, coeff_buf, sz, src_offset, 0);
    md_gpu_pass_barrier_buffer(pass, coeff_buf, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);
}

void md_gto_gpu_coeff_upload_mo_pass(md_gpu_pass_t pass, md_gpu_buffer_t coeff_buf,
    md_gpu_buffer_t src_buf, size_t src_offset, size_t num_mos, size_t num_cgtos) {
    if (!pass || !coeff_buf || !src_buf) return;
    size_t sz = md_gto_gpu_coeff_size_mo(num_mos, num_cgtos);
    md_gpu_pass_copy_buffer(pass, src_buf, coeff_buf, sz, src_offset, 0);
    md_gpu_pass_barrier_buffer(pass, coeff_buf, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);
}

// ---------------------------------------------------------------------------
// GPU dispatch
// ---------------------------------------------------------------------------

void md_gto_gpu_density_pass_record(md_gpu_pass_t pass,
    md_gto_gpu_basis_t gb, md_gpu_buffer_t atom_buf, md_gpu_buffer_t coeff_buf,
    md_gpu_image_t image, const md_grid_t* grid, md_gto_op_t op)
{
    if (!pass || !gb || !atom_buf || !coeff_buf || !image || !grid) {
        MD_LOG_ERROR("md_gto_gpu_density_pass_record: invalid input");
        return;
    }

    const md_gto_basis_layout_t* L = &gb->layout;

    md_gpu_compute_pipeline_t pipeline = gto_pip_density;
    if (!pipeline) {
        MD_LOG_ERROR("md_gto_gpu_density_pass_record: compute pipeline not initialized");
        return;
    }

    typedef struct {
        float    world_to_model[4][4];
        float    index_to_world[4][4];
        float    step[4];
        uint32_t num_cgtos;
        uint32_t op;           /* low bits = SET/ADD/SUB/MAX/MIN, high bits = transforms */
        uint32_t _pad[2];
    } ubo_t;

    ubo_t ubo = {0};
    world_to_model_matrix(ubo.world_to_model, grid);
    index_to_world_matrix(ubo.index_to_world, grid);
    ubo.step[0]      = grid->spacing.elem[0];
    ubo.step[1]      = grid->spacing.elem[1];
    ubo.step[2]      = grid->spacing.elem[2];
    ubo.step[3]      = 0.0f;
    ubo.num_cgtos    = L->num_cgtos;
    ubo.op           = (uint32_t)op;

    size_t sz_cgto_atom_idx = L->off_cgto_r       - L->off_cgto_atom_idx;
    size_t sz_cgto_r        = L->off_cgto_off_len - L->off_cgto_r;
    size_t sz_cgto_off_len  = L->off_pgto         - L->off_cgto_off_len;
    size_t sz_pgto          = (size_t)L->total_size - L->off_pgto;
    size_t sz_atoms         = md_gto_gpu_atom_buffer_size(L->num_atoms);
    size_t tri_len = (L->num_cgtos * (L->num_cgtos + 1)) / 2;
    size_t sz_coeffs = sizeof(float) * tri_len;

    uint32_t wg_size[3] = {
        DIV_UP(grid->dim[0], 8),
        DIV_UP(grid->dim[1], 8),
        DIV_UP(grid->dim[2], 8),
    };

    md_gpu_pass_push_debug_group(pass, "GTO Density Pass");
    md_gpu_pass_bind_compute_pipeline(pass, pipeline);
    md_gpu_pass_push_constants(pass, &ubo, sizeof(ubo));
    md_gpu_pass_bind_buffer_range(pass, 0, gb->buffer, L->off_cgto_atom_idx,  sz_cgto_atom_idx);
    md_gpu_pass_bind_buffer_range(pass, 1, gb->buffer, L->off_cgto_r,         sz_cgto_r);
    md_gpu_pass_bind_buffer_range(pass, 2, gb->buffer, L->off_cgto_off_len,   sz_cgto_off_len);
    md_gpu_pass_bind_buffer_range(pass, 3, gb->buffer, L->off_pgto,           sz_pgto);
    md_gpu_pass_bind_buffer_range(pass, 4, atom_buf,   0,                     sz_atoms);
    md_gpu_pass_bind_buffer_range(pass, 5, coeff_buf,  0,                     sz_coeffs);
    md_gpu_pass_bind_image(pass, 0, image);
    md_gpu_pass_dispatch(pass, wg_size[0], wg_size[1], wg_size[2]);
    md_gpu_pass_pop_debug_group(pass);
}

void md_gto_gpu_mo_pass_record(md_gpu_pass_t pass,
    md_gto_gpu_basis_t gb, md_gpu_buffer_t atom_buf, md_gpu_buffer_t coeff_buf, size_t num_mos,
    md_gpu_image_t out_image, const md_grid_t* grid, md_gto_eval_mode_t eval_mode, md_gto_op_t op)
{
    if (!pass || !gb || !atom_buf || !coeff_buf || !out_image || !grid || num_mos == 0) {
        MD_LOG_ERROR("md_gto_gpu_mo_pass_record: invalid input");
        return;
    }

    const md_gto_basis_layout_t* L = &gb->layout;
    md_gpu_compute_pipeline_t pipeline = gto_pip_mo;
    if (!pipeline) {
        MD_LOG_ERROR("md_gto_gpu_mo_pass_record: compute pipeline not initialized");
        return;
    }

    typedef struct {
        float    world_to_model[4][4];
        float    index_to_world[4][4];
        float    step[4];
        uint32_t num_cgtos;
        uint32_t num_rows;
        uint32_t mode;
        uint32_t op;
    } ubo_t;

    ubo_t ubo = {0};
    world_to_model_matrix(ubo.world_to_model, grid);
    index_to_world_matrix(ubo.index_to_world, grid);
    ubo.step[0]   = grid->spacing.elem[0];
    ubo.step[1]   = grid->spacing.elem[1];
    ubo.step[2]   = grid->spacing.elem[2];
    ubo.step[3]   = 0.0f;
    ubo.num_cgtos = L->num_cgtos;
    ubo.num_rows  = (uint32_t)num_mos;
    ubo.mode      = (eval_mode == MD_GTO_EVAL_MODE_PSI_SQUARED) ? 1u : 0u;
    ubo.op        = (uint32_t)op;

    size_t sz_cgto_atom_idx = L->off_cgto_r         - L->off_cgto_atom_idx;
    size_t sz_cgto_r        = L->off_cgto_off_len   - L->off_cgto_r;
    size_t sz_cgto_off_len  = L->off_pgto           - L->off_cgto_off_len;
    size_t sz_pgto          = (size_t)L->total_size - L->off_pgto;
    size_t sz_atoms         = md_gto_gpu_atom_buffer_size(L->num_atoms);
    size_t sz_coeffs        = sizeof(float) * L->num_cgtos * num_mos;

    uint32_t wg_size[3] = {
        DIV_UP(grid->dim[0], 8),
        DIV_UP(grid->dim[1], 8),
        DIV_UP(grid->dim[2], 8),
    };

    md_gpu_pass_push_debug_group(pass, "GTO MO Pass");
    md_gpu_pass_bind_compute_pipeline(pass, pipeline);
    md_gpu_pass_push_constants(pass, &ubo, sizeof(ubo));
    md_gpu_pass_bind_buffer_range(pass, 0, gb->buffer, L->off_cgto_atom_idx, sz_cgto_atom_idx);
    md_gpu_pass_bind_buffer_range(pass, 1, gb->buffer, L->off_cgto_r,       sz_cgto_r);
    md_gpu_pass_bind_buffer_range(pass, 2, gb->buffer, L->off_cgto_off_len, sz_cgto_off_len);
    md_gpu_pass_bind_buffer_range(pass, 3, gb->buffer, L->off_pgto,         sz_pgto);
    md_gpu_pass_bind_buffer_range(pass, 4, atom_buf,   0,                   sz_atoms);
    md_gpu_pass_bind_buffer_range(pass, 5, coeff_buf,  0,                   sz_coeffs);
    md_gpu_pass_bind_image(pass, 0, out_image);
    md_gpu_pass_dispatch(pass, wg_size[0], wg_size[1], wg_size[2]);
    md_gpu_pass_pop_debug_group(pass);
}

void md_gto_gpu_mo_root_pass_record(md_gpu_pass_t pass,
    md_gto_gpu_basis_t gb, md_gpu_buffer_t atom_buf, md_gpu_buffer_t coeff_buf, size_t num_mos,
    md_gpu_image_t out_image, const md_grid_t* grid, md_gto_eval_mode_t eval_mode, md_gto_op_t op)
{
    if (!pass || !gb || !atom_buf || !coeff_buf || !out_image || !grid || num_mos == 0) {
        MD_LOG_ERROR("md_gto_gpu_mo_root_pass_record: invalid input");
        return;
    }

#if defined(MD_GPU_BACKEND_VULKAN) || defined(MD_GPU_BACKEND_METAL)
    const md_gto_basis_layout_t* L = &gb->layout;
    md_gpu_compute_pipeline_t pipeline = gto_pip_mo_root;
    if (!pipeline) {
        MD_LOG_ERROR("md_gto_gpu_mo_root_pass_record: compute pipeline not initialized");
        return;
    }

    uint64_t basis_addr = md_gpu_buffer_address(gb->buffer);
    uint64_t atom_addr  = md_gpu_buffer_address(atom_buf);
    uint64_t coeff_addr = md_gpu_buffer_address(coeff_buf);
    if (!basis_addr || !atom_addr || !coeff_addr) {
        MD_LOG_ERROR("md_gto_gpu_mo_root_pass_record: buffer device address unavailable");
        return;
    }

    STATIC_ASSERT(sizeof(gto_eval_gto_mo_root_args_t) <= MD_GPU_MAX_PUSH_CONSTANTS, "GTO root pass arguments exceed push constant budget");
    STATIC_ASSERT(gto_eval_gto_mo_root_args_binding_index == 0, "GTO root argument binding index changed");
    STATIC_ASSERT(gto_eval_gto_mo_root_thread_group_size_x == 8, "GTO root pass thread group X changed");
    STATIC_ASSERT(gto_eval_gto_mo_root_thread_group_size_y == 8, "GTO root pass thread group Y changed");
    STATIC_ASSERT(gto_eval_gto_mo_root_thread_group_size_z == 8, "GTO root pass thread group Z changed");
#if defined(MD_GPU_BACKEND_VULKAN)
    STATIC_ASSERT(gto_eval_gto_mo_root_args_binding_kind == MD_GPU_SHADER_BINDING_KIND_PUSH_CONSTANT_BUFFER, "GTO root argument binding kind changed");
#elif defined(MD_GPU_BACKEND_METAL)
    STATIC_ASSERT(gto_eval_gto_mo_root_args_binding_kind == MD_GPU_SHADER_BINDING_KIND_CONSTANT_BUFFER, "GTO root argument binding kind changed");
#endif

    gto_eval_gto_mo_root_args_t args = {0};
    world_to_model_matrix(args.world_to_model, grid);
    index_to_world_matrix(args.index_to_world, grid);
    args.step[0]        = grid->spacing.elem[0];
    args.step[1]        = grid->spacing.elem[1];
    args.step[2]        = grid->spacing.elem[2];
    args.step[3]        = 0.0f;
    args.grid_dim[0]    = (uint32_t)grid->dim[0];
    args.grid_dim[1]    = (uint32_t)grid->dim[1];
    args.grid_dim[2]    = (uint32_t)grid->dim[2];
    args.grid_dim[3]    = 0;
    args.num_cgtos      = L->num_cgtos;
    args.num_rows       = (uint32_t)num_mos;
    args.mode           = (eval_mode == MD_GTO_EVAL_MODE_PSI_SQUARED) ? 1u : 0u;
    args.operation      = (uint32_t)op;
    args.cgto_atom_idx  = basis_addr + L->off_cgto_atom_idx;
    args.cgto_r         = basis_addr + L->off_cgto_r;
    args.cgto_off_len   = basis_addr + L->off_cgto_off_len;
    args.pgto           = basis_addr + L->off_pgto;
    args.atom_xyz       = atom_addr;
    args.coeffs         = coeff_addr;

    uint32_t wg_size[3] = {
        DIV_UP(grid->dim[0], gto_eval_gto_mo_root_thread_group_size_x),
        DIV_UP(grid->dim[1], gto_eval_gto_mo_root_thread_group_size_y),
        DIV_UP(grid->dim[2], gto_eval_gto_mo_root_thread_group_size_z),
    };

    md_gpu_pass_push_debug_group(pass, "GTO MO Root Pass");
    md_gpu_pass_bind_compute_pipeline(pass, pipeline);
    md_gpu_pass_push_constants(pass, &args, sizeof(args));
    md_gpu_pass_bind_image(pass, 0, out_image);
    md_gpu_pass_dispatch(pass, wg_size[0], wg_size[1], wg_size[2]);
    md_gpu_pass_pop_debug_group(pass);
#else
    MD_LOG_ERROR("md_gto_gpu_mo_root_pass_record: root pointer shader is only compiled for Vulkan and Metal backends");
    (void)pass; (void)gb; (void)atom_buf; (void)coeff_buf; (void)num_mos; (void)out_image; (void)grid; (void)eval_mode; (void)op;
#endif
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
    res = md_mm_blendv_ps(res, base1, mask1);
    res = md_mm_blendv_ps(res, base2, mask2);
    res = md_mm_blendv_ps(res, base3, mask3);
    res = md_mm_blendv_ps(res, base4, mask4);
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
    res = md_mm256_blendv_ps(res, base1, mask1);
    res = md_mm256_blendv_ps(res, base2, mask2);
    res = md_mm256_blendv_ps(res, base3, mask3);
    res = md_mm256_blendv_ps(res, base4, mask4);
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

        for (int iz = 0; iz < 8; ++iz) {
            float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
            md_256 vz = md_mm256_set1_ps(z);

            for (int iy = 0; iy < 8; ++iy) {
                float y = grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1];
                md_256 vy = md_mm256_set1_ps(y);

                md_256 dx = md_mm256_sub_ps(vx, px);
                md_256 dy = md_mm256_sub_ps(vy, py);
                md_256 dz = md_mm256_sub_ps(vz, pz);
                md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
                md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));
                md_256 fx = md_mm256_fast_pow(dx, pi);
                md_256 fy = md_mm256_fast_pow(dy, pj);
                md_256 fz = md_mm256_fast_pow(dz, pk);

                md_256 prod_a = md_mm256_mul_ps(pc, fx);
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

        for (int iz = 0; iz < 8; ++iz) {
            md_128 vz = md_mm_set1_ps(grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2]);
            for (int iy = 0; iy < 8; ++iy) {
                md_128 vy = md_mm_set1_ps(grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1]);

                for (int ix = 0; ix < 2; ++ix) {
                    md_128 dx = md_mm_sub_ps(vx[ix], px);
                    md_128 dy = md_mm_sub_ps(vy,	 py);
                    md_128 dz = md_mm_sub_ps(vz,	 pz);
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

    md_temp_t temp = md_temp_begin();
    md_gto_t* sub_gtos = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * num_gtos);

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
