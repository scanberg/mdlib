#include "md_spatial_acc.h"

#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <md_util.h>
#include <md_types.h>

#include <float.h>

static inline float float_from_bits(uint32_t bits) {
    float f;
    MEMCPY(&f, &bits, sizeof(f));
    return f;
}

typedef md_128i ivec4_t;

static inline ivec4_t ivec4_set(int x, int y, int z, int w) {
    return simde_mm_set_epi32(w, z, y, x);
}

static inline ivec4_t ivec4_load(const int* v) {
    return simde_mm_loadu_si128((const md_128i*)v);
}

static inline void ivec4_store(int* v, ivec4_t a) {
    simde_mm_storeu_si128((md_128i*)v, a);
}

static inline ivec4_t ivec4_set1(int v) {
    return simde_mm_set1_epi32(v);
}

static inline ivec4_t ivec4_min(ivec4_t a, ivec4_t b) {
    return simde_mm_min_epi32(a, b);
}

static inline ivec4_t ivec4_max(ivec4_t a, ivec4_t b) {
    return simde_mm_max_epi32(a, b);
}

static inline ivec4_t ivec4_clamp(ivec4_t v, ivec4_t min, ivec4_t max) {
    return ivec4_max(ivec4_min(v, max), min);
}

static inline ivec4_t ivec4_from_vec4(vec4_t v) {
    return simde_mm_cvtps_epi32(v.m128);
}

static inline vec4_t vec4_from_ivec4(ivec4_t v) {
    vec4_t r;
    simde_mm_storeu_ps(r.elem, simde_mm_cvtepi32_ps(v));
    return r;
}

static inline ivec4_t ivec4_add(ivec4_t a, ivec4_t b) {
    return simde_mm_add_epi32(a, b);
}

static inline ivec4_t ivec4_sub(ivec4_t a, ivec4_t b) {
    return simde_mm_sub_epi32(a, b);
}

static inline ivec4_t ivec4_cmpgt(ivec4_t a, ivec4_t b) {
    return simde_mm_cmpgt_epi32(a, b);
}

static inline ivec4_t ivec4_cmplt(ivec4_t a, ivec4_t b) {
    return simde_mm_cmplt_epi32(b, a);
}

static inline ivec4_t ivec4_and(ivec4_t a, ivec4_t b) {
    return simde_mm_and_si128(a, b);
}

static inline ivec4_t ivec4_andnot(ivec4_t a, ivec4_t b) {
    return simde_mm_andnot_si128(b, a);
}

static inline ivec4_t ivec4_or(ivec4_t a, ivec4_t b) {
    return simde_mm_or_si128(a, b);
}

typedef struct {
    float x, y, z;
    uint32_t idx;
} elem_t;

void md_spatial_acc_free(md_spatial_acc_t* acc) {
    ASSERT(acc);
    if (acc->alloc) {
        md_allocator_i* a = acc->alloc;
        if (acc->elem_x)   md_free(a, acc->elem_x,   acc->num_elems * sizeof(float));
        if (acc->elem_y)   md_free(a, acc->elem_y,   acc->num_elems * sizeof(float));
        if (acc->elem_z)   md_free(a, acc->elem_z,   acc->num_elems * sizeof(float));
        if (acc->elem_idx) md_free(a, acc->elem_idx, acc->num_elems * sizeof(uint32_t));
        if (acc->cell_off) md_free(a, acc->cell_off, (acc->num_cells + 1) * sizeof(uint32_t));
    }
    MEMSET(acc, 0, sizeof(md_spatial_acc_t));
}

void md_spatial_acc_init(md_spatial_acc_t* acc, const float* in_x, const float* in_y, const float* in_z, size_t count, double cell_ext, const md_unitcell_t* unitcell) {
    ASSERT(acc);
    if (!acc->alloc) {
        MD_LOG_ERROR("Must have allocator set within spatial acc");
    }

    double A[3][3]  = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double Ai[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    uint32_t flags = 0;

    if (unitcell) {
        md_unitcell_basis_extract(A, unitcell);
        md_unitcell_inv_basis_extract(Ai, unitcell);
        flags = md_unitcell_flags(unitcell);
    }

    vec4_t coord_offset = vec4_zero();

    if ((flags & MD_UNITCELL_PBC_ALL) != MD_UNITCELL_PBC_ALL) {
        // Unit cell either missing or not periodic along one or more axis
        vec4_t aabb_min = {0}, aabb_max = {0};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, in_x, in_y, in_z, NULL, NULL, count);

        // Construct A and Ai from aabb extent
        vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);

        if ((flags & MD_UNITCELL_PBC_X) == 0) {
            coord_offset.x = -aabb_min.x;
            if (aabb_ext.x > 0.0f) {
                A[0][0]  = aabb_ext.x;
                Ai[0][0] = 1.0 / aabb_ext.x;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Y) == 0) {
            coord_offset.y = -aabb_min.y;
            if (aabb_ext.y > 0.0f) {
                A[1][1]  = aabb_ext.y;
                Ai[1][1] = 1.0 / aabb_ext.y;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Z) == 0) {
            coord_offset.z = -aabb_min.z;
            if (aabb_ext.z > 0.0f) {
                A[2][2]  = aabb_ext.z;
                Ai[2][2] = 1.0 / aabb_ext.z;
            }
        }
    }

    // Precompute metric G = A^T A
    double G00 = A[0][0] * A[0][0];
    double G11 = A[1][0] * A[1][0] + A[1][1] * A[1][1];
    double G22 = A[2][0] * A[2][0] + A[2][1] * A[2][1] + A[2][2] * A[2][2];
    double H01 = 2 * A[0][0] * A[1][0];
    double H02 = 2 * A[0][0] * A[2][0];
    double H12 = 2 * (A[1][0] * A[2][0] + A[1][1] * A[2][1]);

    // Choose grid resolution. Heuristic:  cell_dim ≈ |A|/CELL_EXT.
    // Using ~cutoff-ish spacing gives good pruning.
    const double CELL_EXT = cell_ext > 0.0 ? cell_ext : 6.0;

    // Estimate cell_dim by measuring the extents of the box vectors (norms of columns of A)
    // This is only a heuristic for bin counts; the grid is still in fractional space.
    double ax = sqrt(A[0][0] * A[0][0] + A[1][0] * A[1][0] + A[2][0] * A[2][0]);
    double by = sqrt(A[0][1] * A[0][1] + A[1][1] * A[1][1] + A[2][1] * A[2][1]);
    double cz = sqrt(A[0][2] * A[0][2] + A[1][2] * A[1][2] + A[2][2] * A[2][2]);

    uint32_t cell_dim[3] = {
        CLAMP((uint32_t)(ax / CELL_EXT + 0.5), 1, 1024),
        CLAMP((uint32_t)(by / CELL_EXT + 0.5), 1, 1024),
        CLAMP((uint32_t)(cz / CELL_EXT + 0.5), 1, 1024),
    };

    const uint32_t c0  = cell_dim[0];
    const uint32_t c01 = cell_dim[0] * cell_dim[1];
    const size_t num_cells = (size_t)cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Temporary arrays
    size_t pos = md_temp_get_pos();
    uint32_t* local_idx = (uint32_t*)md_temp_push(count * sizeof(uint32_t));
    uint32_t* cell_idx  = (uint32_t*)md_temp_push(count * sizeof(uint32_t));
    elem_t* scratch_s   = (elem_t*)md_temp_push(count * sizeof(elem_t));  // unsorted fractional coords

    // Persistent arrays
    size_t alloc_len = ALIGN_TO(count, 16);

    float*      element_x = md_alloc(acc->alloc, alloc_len * sizeof(float));
    float*      element_y = md_alloc(acc->alloc, alloc_len * sizeof(float));
    float*      element_z = md_alloc(acc->alloc, alloc_len * sizeof(float));
    uint32_t*   element_i = md_alloc(acc->alloc, alloc_len * sizeof(uint32_t));
    uint32_t* cell_offset = md_alloc(acc->alloc, (num_cells + 1) * sizeof(uint32_t));

    const vec4_t fcell_dim  = vec4_set(cell_dim[0], cell_dim[1], cell_dim[2], 0);
    const vec4_t fcell_max  = vec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);
    const ivec4_t icell_min = ivec4_set1(0);
    const ivec4_t icell_max = ivec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);

    ivec4_t cell_min = ivec4_set1(INT32_MAX);
    ivec4_t cell_max = ivec4_set1(INT32_MIN);

    mat4_t Ai4 = {
        (float)Ai[0][0], (float)Ai[0][1], (float)Ai[0][2], 0,
        (float)Ai[1][0], (float)Ai[1][1], (float)Ai[1][2], 0,
        (float)Ai[2][0], (float)Ai[2][1], (float)Ai[2][2], 0,
        0, 0, 0, 0,
    };

    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const vec4_t pbc_mask = vec4_set((flags & MD_UNITCELL_PBC_X) ? val : 0, (flags & MD_UNITCELL_PBC_Y) ? val : 0, (flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);

    // 1) Convert to fractional, wrap periodic axes into [0,1), bin to cells
    for (size_t i = 0; i < count; ++i) {
		uint32_t idx = (uint32_t)i;

        vec4_t r = {in_x[idx], in_y[idx], in_z[idx], 0};
        r = vec4_add(r, coord_offset);

        // Fractional coordinates
        vec4_t c = mat4_mul_vec4(Ai4, r);
        c = vec4_blend(c, vec4_fract(c), pbc_mask);

        // Bin to cell indices
        ivec4_t ic = ivec4_from_vec4(vec4_floor(vec4_mul(c, fcell_dim)));
        ic = ivec4_clamp(ic, icell_min, icell_max);

        cell_min = ivec4_min(cell_min, ic);
        cell_max = ivec4_max(cell_max, ic);

        uint32_t cell_coord[4];
        md_mm_storeu_epi32(cell_coord, ic);
        uint64_t ci = (uint64_t)cell_coord[2] * c01 + (uint64_t)cell_coord[1] * c0 + (uint64_t)cell_coord[0];

        ASSERT(ci < num_cells);

        local_idx[i] = cell_offset[ci]++;  // count for now
        cell_idx[i]  = (uint32_t)ci;

        // stash fractional coordinates
        scratch_s[i] = (elem_t){c.x, c.y, c.z, idx};
    }

    // 2) Prefix sum cell offsets
    uint32_t sum = 0;
    for (size_t ci = 0; ci <= num_cells; ++ci) {
        uint32_t len = cell_offset[ci];
        cell_offset[ci] = sum;
        sum += len;
    }

    // 3) Scatter fractional coords into 'elements' in cell order
    for (size_t i = 0; i < count; ++i) {
        uint32_t dst = cell_offset[cell_idx[i]] + local_idx[i];
        element_x[dst] = scratch_s[i].x;
        element_y[dst] = scratch_s[i].y;
        element_z[dst] = scratch_s[i].z;
        element_i[dst] = scratch_s[i].idx;
    }

    acc->elem_x = element_x;
    acc->elem_y = element_y;
    acc->elem_z = element_z;
    acc->elem_idx = element_i;
    acc->num_elems = count;
    MEMCPY(acc->cell_min, &cell_min, sizeof(acc->cell_min));
    MEMCPY(acc->cell_max, &cell_max, sizeof(acc->cell_max));
    MEMCPY(acc->cell_dim,  cell_dim, sizeof(acc->cell_dim));
    acc->cell_off = cell_offset;
    acc->num_cells = num_cells;
        
    acc->G00 = (float)G00;
    acc->G11 = (float)G11;
    acc->G22 = (float)G22;
    acc->H01 = (float)H01;
    acc->H02 = (float)H02;
    acc->H12 = (float)H12;

    acc->flags = flags;

    md_temp_set_pos_back(pos);
}

static const int FWD_NBRS[13][4] = {
    { 1, 0, 0, 0},  {-1, 1, 0, 0}, {0, 1, 0, 0}, { 1, 1, 0, 0}, {-1, -1, 1, 0}, {0, -1, 1, 0}, {1, -1, 1, 0},
    {-1, 0, 1, 0},  { 0, 0, 1, 0}, {1, 0, 1, 0}, {-1, 1, 1, 0}, { 0,  1, 1, 0}, {1, 1, 1, 0},
};

static inline md_128 distance_squared_tric_sse(md_128 dx, md_128 dy, md_128 dz, md_128 G00, md_128 G11, md_128 G22, md_128 H01, md_128 H02, md_128 H12) {
    md_128 dx2 = md_mm_mul_ps(dx, dx);
    md_128 dy2 = md_mm_mul_ps(dy, dy);
    md_128 dz2 = md_mm_mul_ps(dz, dz);

    md_128 dxy = md_mm_mul_ps(dx, dy);
    md_128 dxz = md_mm_mul_ps(dx, dz);
    md_128 dyz = md_mm_mul_ps(dy, dz);

    md_128 acc   = md_mm_fmadd_ps(G00, dx2, md_mm_fmadd_ps(G11, dy2, md_mm_mul_ps(G22, dz2)));
    md_128 cross = md_mm_fmadd_ps(H01, dxy, md_mm_fmadd_ps(H02, dxz, md_mm_mul_ps(H12, dyz)));
    md_128 d2 = md_mm_add_ps(acc, cross);
    
    return md_mm_add_ps(md_mm_add_ps(dx2, dy2), dz2);
}

static inline md_256 distance_squared_tric_avx(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22, md_256 H01, md_256 H02, md_256 H12) {
    md_256 dx2 = md_mm256_mul_ps(dx, dx);
    md_256 dy2 = md_mm256_mul_ps(dy, dy);
    md_256 dz2 = md_mm256_mul_ps(dz, dz);

    md_256 dxy = md_mm256_mul_ps(dx, dy);
    md_256 dxz = md_mm256_mul_ps(dx, dz);
    md_256 dyz = md_mm256_mul_ps(dy, dz);

    md_256 acc   = md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
    md_256 cross = md_mm256_fmadd_ps(H01, dxy, md_mm256_fmadd_ps(H02, dxz, md_mm256_mul_ps(H12, dyz)));
    md_256 d2 = md_mm256_add_ps(acc, cross);

	return md_mm256_add_ps(md_mm256_add_ps(dx2, dy2), dz2);
}

static inline md_128 distance_squared_ortho_sse(md_128 dx, md_128 dy, md_128 dz, md_128 G00, md_128 G11, md_128 G22) {
    md_128 dx2 = md_mm_mul_ps(dx, dx);
    md_128 dy2 = md_mm_mul_ps(dy, dy);
    md_128 dz2 = md_mm_mul_ps(dz, dz);
    return md_mm_fmadd_ps(G00, dx2, md_mm_fmadd_ps(G11, dy2, md_mm_mul_ps(G22, dz2)));
}

static inline md_256 distance_squared_ortho_avx(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22) {
    md_256 dx2 = md_mm256_mul_ps(dx, dx);
    md_256 dy2 = md_mm256_mul_ps(dy, dy);
    md_256 dz2 = md_mm256_mul_ps(dz, dz);
    return md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
}

// Macros for cell offsets/lengths in the sorted array
#define CELL_INDEX(x, y, z) ((size_t)z * c01 + (size_t)y * c0 + (size_t)x)
#define CELL_OFFSET(ci) (cell_offset[(ci)])
#define CELL_LENGTH(ci) (cell_offset[(ci) + 1] - cell_offset[(ci)])

static void for_each_in_neighboring_cells_triclinic(const md_spatial_acc_t* acc, md_spatial_acc_callback callback, void* user_param) {
	ASSERT(acc);

    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);

    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t num_cells = acc->num_cells;
    const uint32_t cdim[3] = {acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2]};
    const uint32_t cmin[3] = {acc->cell_min[0], acc->cell_min[1], acc->cell_min[2]};
    const uint32_t cmax[3] = {acc->cell_max[0], acc->cell_max[1], acc->cell_max[2]};
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

	ivec4_t cdim_v = ivec4_set(cdim[0], cdim[1], cdim[2], 0);

    ivec4_t pmask_v = ivec4_add(ivec4_set1(1), ivec4_set(
        (acc->flags & MD_UNITCELL_PBC_X) ? 0xFFFFFFFFU : 0,
        (acc->flags & MD_UNITCELL_PBC_Y) ? 0xFFFFFFFFU : 0,
        (acc->flags & MD_UNITCELL_PBC_Z) ? 0xFFFFFFFFU : 0,
		0));

    const bool periodic_x = acc->flags & MD_UNITCELL_PBC_X;
    const bool periodic_y = acc->flags & MD_UNITCELL_PBC_Y;
    const bool periodic_z = acc->flags & MD_UNITCELL_PBC_Z;

	// Temporary storage for distances
    float ij_dist2[8];

    for (uint32_t cz = cmin[2]; cz <= cmax[2]; ++cz) {
        for (uint32_t cy = cmin[1]; cy <= cmax[1]; ++cy) {
            for (uint32_t cx = cmin[0]; cx <= cmax[0]; ++cx) {

                const size_t ci    = CELL_INDEX(cx, cy, cz);
                const size_t off_i = CELL_OFFSET(ci);
                const size_t len_i = CELL_LENGTH(ci);

                if (len_i == 0) continue;

                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                // Self cell: only j > i
                for (size_t i = 0; i + 1 < len_i; ++i) {
                    const md_256 x = md_mm256_set1_ps(elem_i_x[i]);
                    const md_256 y = md_mm256_set1_ps(elem_i_y[i]);
                    const md_256 z = md_mm256_set1_ps(elem_i_z[i]);
                    const uint32_t idx = elem_i_idx[i];

                    for (size_t j = i + 1; j < len_i; j += 8) {
                        const md_256 dx = md_mm256_sub_ps(x, md_mm256_loadu_ps(elem_i_x + j));
                        const md_256 dy = md_mm256_sub_ps(y, md_mm256_loadu_ps(elem_i_y + j));
                        const md_256 dz = md_mm256_sub_ps(z, md_mm256_loadu_ps(elem_i_z + j));
                        const md_256 d2 = distance_squared_tric_avx(dx, dy, dz, G00, G11, G22, H01, H02, H12);

						size_t len = MIN(len_i - j, 8);

                        // Invoke callback
						callback(idx, elem_i_idx + j, d2, user_param);
                    }
                }

				ivec4_t c_v = ivec4_set(cx, cy, cz, 0);

                // Forward neighbors
                for (uint32_t n = 0; n < ARRAY_SIZE(FWD_NBRS); ++n) {
					ivec4_t fwd_v = ivec4_load(FWD_NBRS[n]);
                    
					ivec4_t n_v = ivec4_add(c_v, fwd_v);
					ivec4_t wrap_upper = ivec4_cmpgt(n_v, ivec4_sub(cdim_v, ivec4_set1(1)));
					ivec4_t wrap_lower = ivec4_cmplt(n_v, ivec4_set1(0));
					
                    // Check if we should skip (if non periodic)
					//bool skip = md_mm_movemask_epi8(ivec4_andnot(ivec4_or(wrap_upper, wrap_lower), pmask_v));
					//if (skip) continue;

					// Wrap indices
					n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
					n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

					// Calculate shifts for periodic wrapping
                    vec4_t shift_v = vec4_set1(0);
                    shift_v = vec4_add(shift_v, vec4_and(vec4_from_ivec4(wrap_lower), vec4_set1(1)));
                    shift_v = vec4_sub(shift_v, vec4_and(vec4_from_ivec4(wrap_upper), vec4_set1(1)));

                    int n[4], shift[4];
					ivec4_store(n, n_v);
					vec4_store(shift, shift_v);

					const md_256 shift_x = md_mm256_set1_ps(shift[0]);
					const md_256 shift_y = md_mm256_set1_ps(shift[1]);
					const md_256 shift_z = md_mm256_set1_ps(shift[2]);

                    const uint32_t cj    = CELL_INDEX(n[0], n[1], n[2]);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);

                    if (len_j == 0) continue;

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (size_t i = 0; i < len_i; ++i) {
                        const md_256 x = md_mm256_add_ps(md_mm256_set1_ps(elem_i_x[i]), shift_x);
                        const md_256 y = md_mm256_add_ps(md_mm256_set1_ps(elem_i_y[i]), shift_y);
                        const md_256 z = md_mm256_add_ps(md_mm256_set1_ps(elem_i_z[i]), shift_z);
                        const uint32_t idx = elem_i_idx[i];

                        for (size_t j = 0; j < len_j; j += 8) {
                            const md_256 dx = md_mm256_sub_ps(x, md_mm256_loadu_ps(elem_j_x + j));
                            const md_256 dy = md_mm256_sub_ps(y, md_mm256_loadu_ps(elem_j_y + j));
                            const md_256 dz = md_mm256_sub_ps(z, md_mm256_loadu_ps(elem_j_z + j));
                            const md_256 d2 = distance_squared_tric_avx(dx, dy, dz, G00, G11, G22, H01, H02, H12);

							size_t len = MIN(len_j - j, 8);

                            // Invoke callback
                            callback(idx, elem_j_idx + j, d2, user_param);
                        }
                    }
                }
            }
        }
    }
}

static void for_each_in_neighboring_cells_ortho(const md_spatial_acc_t* acc, md_spatial_acc_callback callback, void* user_param) {
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t num_cells = acc->num_cells;
    const uint32_t cdim[3] = {acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2]};
    const uint32_t cmin[3] = {acc->cell_min[0], acc->cell_min[1], acc->cell_min[2]};
    const uint32_t cmax[3] = {acc->cell_max[0], acc->cell_max[1], acc->cell_max[2]};
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    ivec4_t cdim_v = ivec4_set(cdim[0], cdim[1], cdim[2], 0);

    ivec4_t pmask_v = ivec4_add(ivec4_set1(1), ivec4_set(
        (acc->flags & MD_UNITCELL_PBC_X) ? 0xFFFFFFFFU : 0,
        (acc->flags & MD_UNITCELL_PBC_Y) ? 0xFFFFFFFFU : 0,
        (acc->flags & MD_UNITCELL_PBC_Z) ? 0xFFFFFFFFU : 0,
        0));

    const bool periodic_x = acc->flags & MD_UNITCELL_PBC_X;
    const bool periodic_y = acc->flags & MD_UNITCELL_PBC_Y;
    const bool periodic_z = acc->flags & MD_UNITCELL_PBC_Z;

    // Temporary storage for distances
    float ij_dist2[8];

    const md_256i inc = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    const md_256  inf = md_mm256_set1_ps(FLT_MAX);

    for (uint32_t cz = cmin[2]; cz <= cmax[2]; ++cz) {
        for (uint32_t cy = cmin[1]; cy <= cmax[1]; ++cy) {
            for (uint32_t cx = cmin[0]; cx <= cmax[0]; ++cx) {

                const uint32_t ci    = CELL_INDEX(cx, cy, cz);
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);

                if (len_i == 0) continue;

                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                md_256i v_j_max = md_mm256_set1_epi32(len_i-1);

                // Self cell: only j > i
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    const md_256 x = md_mm256_set1_ps(elem_i_x[i]);
                    const md_256 y = md_mm256_set1_ps(elem_i_y[i]);
                    const md_256 z = md_mm256_set1_ps(elem_i_z[i]);
                    const uint32_t idx = elem_i_idx[i];

                    md_256i v_j = md_mm256_add_epi32(md_mm256_set1_epi32(i + 1), inc);
                    for (uint32_t j = i + 1; j < len_i; j += 8) {
                        const md_256i cmp_mask_i = md_mm256_cmpgt_epi32(v_j, v_j_max);   // invalid = 0xFFFFFFFF
                        const md_256 invalid_mask = md_mm256_castsi256_ps(cmp_mask_i);

                        const md_256 dx = md_mm256_sub_ps(x, md_mm256_loadu_ps(elem_i_x + j));
                        const md_256 dy = md_mm256_sub_ps(y, md_mm256_loadu_ps(elem_i_y + j));
                        const md_256 dz = md_mm256_sub_ps(z, md_mm256_loadu_ps(elem_i_z + j));
                        const md_256 d2 = distance_squared_ortho_avx(dx, dy, dz, G00, G11, G22);

                        const md_256 d2_masked = md_mm256_blendv_ps(d2, inf, invalid_mask);

                        // Invoke callback
                        callback(idx, elem_i_idx + j, d2_masked, user_param);

                        v_j = md_mm256_add_epi32(v_j, md_mm256_set1_epi32(8));
                    }
                }

#if 1
                ivec4_t c_v = ivec4_set(cx, cy, cz, 0);

                // Forward neighbors
                for (uint32_t n = 0; n < ARRAY_SIZE(FWD_NBRS); ++n) {
                    ivec4_t fwd_v = ivec4_load(FWD_NBRS[n]);

                    ivec4_t n_v = ivec4_add(c_v, fwd_v);
                    ivec4_t wrap_upper = ivec4_cmpgt(n_v, ivec4_sub(cdim_v, ivec4_set1(1)));
                    ivec4_t wrap_lower = ivec4_cmplt(n_v, ivec4_set1(0));

                    // Check if we should skip (if non periodic)
                    //bool skip = md_mm_movemask_epi8(ivec4_andnot(ivec4_or(wrap_upper, wrap_lower), pmask_v)) != 0;
                    //if (skip) continue;

                    // Wrap indices
                    n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
                    n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

                    int n[4];
                    ivec4_store(n, n_v);

                    const uint32_t cj    = CELL_INDEX(n[0], n[1], n[2]);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);

                    if (len_j == 0) continue;

                    // Calculate shifts for periodic wrapping
                    ivec4_t shift_i = ivec4_sub(ivec4_and(wrap_lower, ivec4_set1(1)),
                                                ivec4_and(wrap_upper, ivec4_set1(1))); //  +1 for lower wrap, -1 for upper wrap

                    int shift_i_arr[4];
                    ivec4_store(shift_i_arr, shift_i);

                    const md_256 shift_x = md_mm256_set1_ps((float)shift_i_arr[0]);
                    const md_256 shift_y = md_mm256_set1_ps((float)shift_i_arr[1]);
                    const md_256 shift_z = md_mm256_set1_ps((float)shift_i_arr[2]);

                    v_j_max = md_mm256_set1_epi32(len_j-1);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (size_t i = 0; i < len_i; ++i) {
                        const md_256 x = md_mm256_add_ps(md_mm256_set1_ps(elem_i_x[i]), shift_x);
                        const md_256 y = md_mm256_add_ps(md_mm256_set1_ps(elem_i_y[i]), shift_y);
                        const md_256 z = md_mm256_add_ps(md_mm256_set1_ps(elem_i_z[i]), shift_z);
                        const uint32_t idx = elem_i_idx[i];

                        md_256i v_j = inc;
                        for (size_t j = 0; j < len_j; j += 8) {
                            const md_256i cmp_mask_i = md_mm256_cmpgt_epi32(v_j, v_j_max);   // invalid = 0xFFFFFFFF
                            const md_256 invalid_mask = md_mm256_castsi256_ps(cmp_mask_i);

                            const md_256 dx = md_mm256_sub_ps(x, md_mm256_loadu_ps(elem_j_x + j));
                            const md_256 dy = md_mm256_sub_ps(y, md_mm256_loadu_ps(elem_j_y + j));
                            const md_256 dz = md_mm256_sub_ps(z, md_mm256_loadu_ps(elem_j_z + j));
                            const md_256 d2 = distance_squared_ortho_avx(dx, dy, dz, G00, G11, G22);

                            const md_256 d2_masked = md_mm256_blendv_ps(d2, inf, invalid_mask);

                            // Invoke callback
                            callback(idx, elem_j_idx + j, d2_masked, user_param);

                            v_j = md_mm256_add_epi32(v_j, md_mm256_set1_epi32(8));
                        }
                    }
                }
                #endif
            }
        }
    }
}

#undef CELL_INDEX
#undef CELL_OFFSET
#undef CELL_LENGTH

void md_spatial_acc_for_each_in_neighboring_cells(const md_spatial_acc_t* acc, md_spatial_acc_callback callback, void* user_param) {
    ASSERT(acc);
	ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_in_neighboring_cells_triclinic(acc, callback, user_param);
    } else {
		for_each_in_neighboring_cells_ortho(acc, callback, user_param);
	}
}
