#include "md_spatial_acc.h"

#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_intrinsics.h>
#include <md_util.h>
#include <md_types.h>

#include <float.h>
#include <math.h>

#define SPATIAL_ACC_BUFLEN 1024
#define SPATIAL_ACC_MAX_NEIGHBOR_CELLS 5
#define SPATIAL_ACC_MAX_CELLS_PER_DIM 1024

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
    r.m128 = md_mm_cvtepi32_ps(v);
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
    return simde_mm_cmplt_epi32(a, b);
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

static inline ivec4_t ivec4_xor(ivec4_t a, ivec4_t b) {
    return simde_mm_xor_si128(a, b);
}

static inline bool ivec4_any(ivec4_t v) {
#if defined(__aarch64__)
    // Use NEON reduction: true if any 32-bit lane is non-zero.
    simde__m128i_private vp = simde__m128i_to_private(v);
    return vmaxvq_u32(vp.neon_u32) != 0u;
#else
    // Efficient on x86: ptest(a,a) via SIMDe
    return !simde_mm_testz_si128(v, v);
#endif
}

typedef struct {
    float x, y, z;
    uint32_t idx;
} elem_t;

static inline vec4_t md_coord_stream_load_vec4(const md_coord_stream_t* s, size_t i) {
    const size_t src = s->idx ? (size_t)s->idx[i] : i;

    if (s->layout == MD_COORD_STREAM_LAYOUT_SOA) {
        return vec4_set(s->soa.x[src], s->soa.y[src], s->soa.z[src], 0.0f);
    } else {
        const char* p = (const char*)s->aos.base + src * s->aos.stride;
        vec4_t v = {0};
        MEMCPY(&v, p, sizeof(float) * 3);
        return v;
    }
}

static inline int32_t md_coord_stream_load_idx(const md_coord_stream_t* s, size_t i) {
    return s->idx ? s->idx[i] : (int32_t)i;
}

void md_spatial_acc_free(md_spatial_acc_t* acc) {
    ASSERT(acc);
    if (acc->alloc) {
        if (acc->elem_x)   md_array_free(acc->elem_x,   acc->alloc);
        if (acc->elem_y)   md_array_free(acc->elem_y,   acc->alloc);
        if (acc->elem_z)   md_array_free(acc->elem_z,   acc->alloc);
        if (acc->elem_idx) md_array_free(acc->elem_idx, acc->alloc);
        if (acc->cell_off) md_array_free(acc->cell_off, acc->alloc);
    }
    MEMSET(acc, 0, sizeof(md_spatial_acc_t));
}

static void md_spatial_acc_reset(md_spatial_acc_t* acc) {
    ASSERT(acc);
    md_array_shrink(acc->elem_x, 0);
	md_array_shrink(acc->elem_y, 0);
	md_array_shrink(acc->elem_z, 0);
	md_array_shrink(acc->elem_idx, 0);
	md_array_shrink(acc->cell_off, 0);

    acc->num_cells = 0;
    acc->num_elems = 0;

    acc->G00 = acc->G11 = acc->G22 = 0.0f;
    acc->H01 = acc->H02 = acc->H12 = 0.0f;

	MEMSET(acc->cell_mask, 0, sizeof(acc->cell_mask));
    MEMSET(acc->cell_dim,  0, sizeof(acc->cell_dim));

	acc->flags = 0;
	MEMSET(acc->origin, 0, sizeof(acc->origin));
}

void md_spatial_acc_init(md_spatial_acc_t* acc, const md_coord_stream_t* stream, double in_cell_ext, const md_unitcell_t* in_unitcell, md_spatial_acc_flags_t in_flags) {
    ASSERT(acc);
    ASSERT(stream);

    if (!acc->alloc) {
        MD_LOG_ERROR("Must have allocator set within spatial acc");
        return;
    }

    if (in_flags & MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX) {
        if (!stream->idx) {
            MD_LOG_ERROR("Flag MD_SPATIAL_ACC_FLAG_STORE_SUPPLIED_IDX is set, but coordinate stream index is not supplied");
            return;
        }
    }

	// Reset acc, but to not free memory
    md_spatial_acc_reset(acc);

    if (stream->count == 0) {
		// Not a real error, but nothing to build. Leave acc in a valid empty state.
        return;
    }

    if (in_cell_ext <= 0.0) {
        // Fallback to some default value
        in_cell_ext = 6.0;
    }

    // Choose grid resolution. Heuristic:  cell_dim ≈ |A|/CELL_EXT.
    // Using ~cutoff-ish spacing gives good pruning.
    // We protect against too small cell_ext to avoid excessive memory usage and slowdown from too many cells.
    const double CELL_EXT = MAX(in_cell_ext, 3.0);

    double A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    uint32_t flags = 0;

    if (in_unitcell) {
        md_unitcell_A_extract_double(A, in_unitcell);
        md_unitcell_I_extract_double(I, in_unitcell);
        flags = md_unitcell_flags(in_unitcell);
    }

    vec4_t origin = vec4_zero();

    if ((flags & MD_UNITCELL_PBC_ALL) != MD_UNITCELL_PBC_ALL) {
        ASSERT((flags & MD_UNITCELL_TRICLINIC) == 0);
        // Unit cell either missing or not periodic along one or more axis
        vec4_t aabb_min = {0}, aabb_max = {0};
        for (size_t i = 0; i < stream->count; i++) {
            vec4_t v = md_coord_stream_load_vec4(stream, i);
            aabb_min = vec4_min(aabb_min, v);
            aabb_max = vec4_max(aabb_max, v);
        }

        // Construct A and I from aabb extent
        vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);

        // Round up to nearest N * CELL_EXT
        aabb_ext = vec4_mul_f(vec4_ceil(vec4_div_f(aabb_ext, (float)CELL_EXT)), (float)CELL_EXT);

        // Set min as center - half extent
        vec4_t aabb_center = vec4_mul_f(vec4_add(aabb_min, aabb_max), 0.5f);
        aabb_min = vec4_sub(aabb_center, vec4_mul_f(aabb_ext, 0.5f));
        aabb_max = vec4_add(aabb_center, vec4_mul_f(aabb_ext, 0.5f));

        if ((flags & MD_UNITCELL_PBC_X) == 0) {
            origin.x = aabb_min.x;
            if (aabb_ext.x > 0.0f) {
                A[0][0] = aabb_ext.x;
                I[0][0] = 1.0 / aabb_ext.x;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Y) == 0) {
            origin.y = aabb_min.y;
            if (aabb_ext.y > 0.0f) {
                A[1][1] = aabb_ext.y;
                I[1][1] = 1.0 / aabb_ext.y;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Z) == 0) {
            origin.z = aabb_min.z;
            if (aabb_ext.z > 0.0f) {
                A[2][2] = aabb_ext.z;
                I[2][2] = 1.0 / aabb_ext.z;
            }
        }
    }

    // Basis vectors are the COLUMNS of A (column-major storage: A[col][row])
    double a[3] = { A[0][0], A[0][1], A[0][2] };
    double b[3] = { A[1][0], A[1][1], A[1][2] };
    double c[3] = { A[2][0], A[2][1], A[2][2] };

    // Metric G = A^T A
    double G00 = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); // dot(a, a)
    double G11 = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]); // dot(b, b)
    double G22 = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]); // dot(c, c)
    double G01 = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]); // dot(a, b)
    double G02 = (a[0] * c[0] + a[1] * c[1] + a[2] * c[2]); // dot(a, c)
    double G12 = (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]); // dot(b, c)

    double H01 = 0.0, H02 = 0.0, H12 = 0.0;

    double norm_a = sqrt(G00);
    double norm_b = sqrt(G11);
    double norm_c = sqrt(G22);

    float inv_cell_ext[3] = {
        (float)(norm_a > 0.0 ? 1.0 / norm_a : 0.0),
        (float)(norm_b > 0.0 ? 1.0 / norm_b : 0.0),
        (float)(norm_c > 0.0 ? 1.0 / norm_c : 0.0),
    };

    if (flags & MD_UNITCELL_TRICLINIC) {
        H01 = 2.0 * G01;
        H02 = 2.0 * G02;
        H12 = 2.0 * G12;

        double det = G00 * (G11 * G22 - G12 * G12) - G01 * (G01 * G22 - G12 * G02) + G02 * (G01 * G12 - G11 * G02);

        if (det < DBL_EPSILON) {
            MD_LOG_ERROR("Degenerate unit cell / A matrix provided to spatial acc");
            md_spatial_acc_reset(acc);
            return;
        }

        double GI00 = (G11 * G22 - G12 * G12) / det;
        double GI11 = (G00 * G22 - G02 * G02) / det;
        double GI22 = (G00 * G11 - G01 * G01) / det;

        inv_cell_ext[0] = (float)sqrt(GI00);
        inv_cell_ext[1] = (float)sqrt(GI11);
        inv_cell_ext[2] = (float)sqrt(GI22);
    }

    // Estimate cell_dim by measuring the extents of the box vectors (norms of columns of A)
    // This is only a heuristic for bin counts; the grid is still in fractional space.
    uint32_t cell_dim[3] = {
        CLAMP((uint32_t)(norm_a / CELL_EXT), 1, SPATIAL_ACC_MAX_CELLS_PER_DIM),
        CLAMP((uint32_t)(norm_b / CELL_EXT), 1, SPATIAL_ACC_MAX_CELLS_PER_DIM),
        CLAMP((uint32_t)(norm_c / CELL_EXT), 1, SPATIAL_ACC_MAX_CELLS_PER_DIM),
    };

    const uint32_t c0  = cell_dim[0];
    const uint32_t c01 = cell_dim[0] * cell_dim[1];
    const size_t num_cells = (size_t)cell_dim[0] * cell_dim[1] * cell_dim[2];

#if DEBUG
    MD_LOG_DEBUG("cell_dim: %i %i %i", cell_dim[0], cell_dim[1], cell_dim[2]);
#endif

    // Temporary arrays
    const size_t temp_arena_page_size = MEGABYTES(4);
    md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), temp_arena_page_size);
    uint32_t* local_idx = (uint32_t*)md_arena_allocator_push(temp_arena, stream->count * sizeof(uint32_t));
    uint32_t* cell_idx  = (uint32_t*)md_arena_allocator_push(temp_arena, stream->count * sizeof(uint32_t));
    elem_t* scratch_s   = (elem_t*)  md_arena_allocator_push(temp_arena, stream->count * sizeof(elem_t));  // unsorted fractional coords

    // Resize / allocate persistent arrays
    size_t alloc_len = ALIGN_TO(stream->count, 16);

    md_array_resize(acc->elem_x, alloc_len, acc->alloc);
    md_array_resize(acc->elem_y, alloc_len, acc->alloc);
    md_array_resize(acc->elem_z, alloc_len, acc->alloc);
    md_array_resize(acc->elem_idx, alloc_len, acc->alloc);
    md_array_resize(acc->cell_off, num_cells + 1, acc->alloc);
    MEMSET(acc->cell_off, 0, (num_cells + 1) * sizeof(uint32_t));

    const vec4_t  fcell_dim = vec4_set((float)cell_dim[0], (float)cell_dim[1], (float)cell_dim[2], 0);
    const ivec4_t icell_min = ivec4_set1(0);
    const ivec4_t icell_max = ivec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);

    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const vec4_t pbc_mask = vec4_set((flags & MD_UNITCELL_PBC_X) ? val : 0, (flags & MD_UNITCELL_PBC_Y) ? val : 0, (flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);

    vec4_t vI[3] = {
        vec4_set((float)I[0][0], (float)I[0][1], (float)I[0][2], 0),
        vec4_set((float)I[1][0], (float)I[1][1], (float)I[1][2], 0),
        vec4_set((float)I[2][0], (float)I[2][1], (float)I[2][2], 0),
    };

	uint64_t cell_mask[3][16] = { 0 };

    // 1) Convert to fractional, wrap periodic axes into [0,1), bin to cells
    for (size_t i = 0; i < stream->count; ++i) {
        uint32_t idx = md_coord_stream_load_idx(stream, i);
        vec4_t r     = md_coord_stream_load_vec4(stream, i);

        // Fractional coordinates
        vec4_t s = vec4_linear_combine_3(vec4_sub(r, origin), vI);
        
        s = vec4_blend(s, vec4_fract(s), pbc_mask);

        // Bin to cell indices
        ivec4_t ic = ivec4_from_vec4(vec4_floor(vec4_mul(s, fcell_dim)));
        ic = ivec4_clamp(ic, icell_min, icell_max);

        uint32_t cell_coord[4];
        md_mm_storeu_epi32(cell_coord, ic);
        uint64_t ci = (uint64_t)cell_coord[2] * c01 + (uint64_t)cell_coord[1] * c0 + (uint64_t)cell_coord[0];

		cell_mask[0][cell_coord[0] / 64] |= (1ULL << (cell_coord[0] & 63));
		cell_mask[1][cell_coord[1] / 64] |= (1ULL << (cell_coord[1] & 63));
		cell_mask[2][cell_coord[2] / 64] |= (1ULL << (cell_coord[2] & 63));

        ASSERT(ci < num_cells);

        local_idx[i] = acc->cell_off[ci]++;  // count for now
        cell_idx[i]  = (uint32_t)ci;

        uint32_t elem_idx = (in_flags & MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX) ? idx : (uint32_t)i;

        // stash fractional coordinates
        scratch_s[i] = (elem_t){s.x, s.y, s.z, elem_idx};
    }

    // 2) Prefix sum cell offsets
    uint32_t sum = 0;
    for (size_t ci = 0; ci <= num_cells; ++ci) {
        uint32_t len = acc->cell_off[ci];
        acc->cell_off[ci] = sum;
        sum += len;
    }
    ASSERT(sum == stream->count);

    // 3) Scatter fractional coords into 'elements' in cell order
    for (size_t i = 0; i < stream->count; ++i) {
        uint32_t dst = acc->cell_off[cell_idx[i]] + local_idx[i];
        ASSERT(dst < stream->count);
        acc->elem_x[dst] = scratch_s[i].x;
        acc->elem_y[dst] = scratch_s[i].y;
        acc->elem_z[dst] = scratch_s[i].z;
        acc->elem_idx[dst] = scratch_s[i].idx;
    }

    acc->num_elems = stream->count;

    MEMCPY(acc->cell_mask, cell_mask, sizeof(acc->cell_mask));
    MEMCPY(acc->cell_dim,  cell_dim,  sizeof(acc->cell_dim));
    MEMCPY(acc->inv_cell_ext, inv_cell_ext, sizeof(acc->inv_cell_ext));
    acc->num_cells = num_cells;

    acc->G00 = (float)G00;
    acc->G11 = (float)G11;
    acc->G22 = (float)G22;
    acc->H01 = (float)H01;
    acc->H02 = (float)H02;
    acc->H12 = (float)H12;

	acc->A[0][0] = (float)A[0][0];
	acc->A[0][1] = (float)A[0][1];
	acc->A[0][2] = (float)A[0][2];

	acc->A[1][0] = (float)A[1][0];
	acc->A[1][1] = (float)A[1][1];
	acc->A[1][2] = (float)A[1][2];

	acc->A[2][0] = (float)A[2][0];
	acc->A[2][1] = (float)A[2][1];
	acc->A[2][2] = (float)A[2][2];

    acc->I[0][0] = (float)I[0][0];
    acc->I[0][1] = (float)I[0][1];
    acc->I[0][2] = (float)I[0][2];

    acc->I[1][0] = (float)I[1][0];
    acc->I[1][1] = (float)I[1][1];
    acc->I[1][2] = (float)I[1][2];

    acc->I[2][0] = (float)I[2][0];
    acc->I[2][1] = (float)I[2][1];
    acc->I[2][2] = (float)I[2][2];

    // Persist origin offset used to construct fractional frame
    acc->origin[0] = origin.x;
    acc->origin[1] = origin.y;
    acc->origin[2] = origin.z;

    acc->flags = flags;

    md_arena_allocator_destroy(temp_arena);
}

// Generate forward neighbor offsets for a 3D grid cell
// - out: user-provided array of size at least (ncell*2+1)^3 - 1
// - ncell: number of cells in the neighborhood (1 → 1-cell, 2 → 2-cell, etc.)
// - Each element in out is int[4] representing offset {dx, dy, dz, 0}
// Returns the number of neighbors written
static inline size_t generate_forward_neighbors4(int out[][4], const int ncell[3]) {
    size_t count = 0;

    for (int dz = 0; dz <= ncell[2]; ++dz) {        // memory order: z-major
        for (int dy = -ncell[1]; dy <= ncell[1]; ++dy) {
            for (int dx = -ncell[0]; dx <= ncell[0]; ++dx) {
                // skip the origin
                if (dx == 0 && dy == 0 && dz == 0) continue;

                // Forward neighbor condition: i < j
                if (dz > 0 || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0)) {
                    out[count][0] = dx;
                    out[count][1] = dy;
                    out[count][2] = dz;
					out[count][3] = 0;
                    count++;
                }
            }
        }
    }
    return count;
}

// Generate full neighbor offsets for a 3D grid cell
// - out: user-provided array of size at least (ncell*2+1)^3
// - ncell: number of cells in the neighborhood (1 → 1-cell, 2 → 2-cell, etc.)
// - Each element in out is int[4] representing offset {dx, dy, dz, 0}
// Returns the number of neighbors written
static inline size_t generate_neighbors4(int out[][4], const int ncell[3]) {
    size_t count = 0;
    for (int dz = -ncell[2]; dz <= ncell[2]; ++dz) {
        for (int dy = -ncell[1]; dy <= ncell[1]; ++dy) {
            for (int dx = -ncell[0]; dx <= ncell[0]; ++dx) {
                out[count][0] = dx;
                out[count][1] = dy;
                out[count][2] = dz;
                out[count][3] = 0;
                ++count;
            }
        }
    }
    return count;
}

static inline md_128 distance_squared_tri_128(md_128 dx, md_128 dy, md_128 dz, md_128 G00, md_128 G11, md_128 G22, md_128 H01, md_128 H02, md_128 H12) {
    md_128 dx2 = md_mm_mul_ps(dx, dx);
    md_128 dy2 = md_mm_mul_ps(dy, dy);
    md_128 dz2 = md_mm_mul_ps(dz, dz);

    md_128 dxy = md_mm_mul_ps(dx, dy);
    md_128 dxz = md_mm_mul_ps(dx, dz);
    md_128 dyz = md_mm_mul_ps(dy, dz);

    md_128 acc   = md_mm_fmadd_ps(G00, dx2, md_mm_fmadd_ps(G11, dy2, md_mm_mul_ps(G22, dz2)));
    md_128 cross = md_mm_fmadd_ps(H01, dxy, md_mm_fmadd_ps(H02, dxz, md_mm_mul_ps(H12, dyz)));
    return md_mm_add_ps(acc, cross);
}

static inline md_256 distance_squared_tri_256(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22, md_256 H01, md_256 H02, md_256 H12) {
    md_256 dx2 = md_mm256_mul_ps(dx, dx);
    md_256 dy2 = md_mm256_mul_ps(dy, dy);
    md_256 dz2 = md_mm256_mul_ps(dz, dz);

    md_256 dxy = md_mm256_mul_ps(dx, dy);
    md_256 dxz = md_mm256_mul_ps(dx, dz);
    md_256 dyz = md_mm256_mul_ps(dy, dz);

    md_256 acc   = md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
    md_256 cross = md_mm256_fmadd_ps(H01, dxy, md_mm256_fmadd_ps(H02, dxz, md_mm256_mul_ps(H12, dyz)));
    return md_mm256_add_ps(acc, cross);
}

static inline md_128 distance_squared_ort_128(md_128 dx, md_128 dy, md_128 dz, md_128 G00, md_128 G11, md_128 G22) {
    md_128 dx2 = md_mm_mul_ps(dx, dx);
    md_128 dy2 = md_mm_mul_ps(dy, dy);
    md_128 dz2 = md_mm_mul_ps(dz, dz);
    return md_mm_fmadd_ps(G00, dx2, md_mm_fmadd_ps(G11, dy2, md_mm_mul_ps(G22, dz2)));
}

static inline md_256 distance_squared_ort_256(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22) {
    md_256 dx2 = md_mm256_mul_ps(dx, dx);
    md_256 dy2 = md_mm256_mul_ps(dy, dy);
    md_256 dz2 = md_mm256_mul_ps(dz, dz);
    return md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
}

static inline int wrap_coord(int x, int N) {
    x += (x <  0) ? N : 0;
    x -= (x >= N) ? N : 0;
    return x;
}

static inline int isign(int a) {
    return (a > 0) - (a < 0);
}

static float calc_r2(double cutoff) {
    float r2 = (float)(cutoff * cutoff);
    return nextafterf(r2, r2 + 1.0f); // Round up to ensure we don't miss neighbors due to floating point precision
}

static inline vec4_t vec4_cart_to_fract(vec4_t in_c, const md_spatial_acc_t* acc) {
    const vec4_t t = vec4_set(acc->origin[0], acc->origin[1], acc->origin[2], 0);
    const vec4_t I[3] = {
        vec4_set(acc->I[0][0], acc->I[0][1], acc->I[0][2], 0),
        vec4_set(acc->I[1][0], acc->I[1][1], acc->I[1][2], 0),
        vec4_set(acc->I[2][0], acc->I[2][1], acc->I[2][2], 0),
    };
    return vec4_linear_combine_3(vec4_sub(in_c, t), I);
}

static inline void cart_to_fract(double out_s[3], const double in_c[3], const md_spatial_acc_t* acc) {
    // I is indexed as I[col][row]
    // fract = I * (cart - origin)
    const double x = in_c[0] - acc->origin[0];
    const double y = in_c[1] - acc->origin[1];
    const double z = in_c[2] - acc->origin[2];
    out_s[0] = acc->I[0][0] * x + acc->I[1][0] * y + acc->I[2][0] * z;
    out_s[1] = acc->I[0][1] * x + acc->I[1][1] * y + acc->I[2][1] * z;
    out_s[2] = acc->I[0][2] * x + acc->I[1][2] * y + acc->I[2][2] * z;
}

static inline void fract_to_cart(double out_x[3], const double in_s[3], const md_spatial_acc_t* acc) {
    // A is indexed as A[col][row]
    // cart = A * fract + origin
    out_x[0] = acc->A[0][0] * in_s[0] + acc->A[1][0] * in_s[1] + acc->A[2][0] * in_s[2] + acc->origin[0];
    out_x[1] = acc->A[0][1] * in_s[0] + acc->A[1][1] * in_s[1] + acc->A[2][1] * in_s[2] + acc->origin[1];
    out_x[2] = acc->A[0][2] * in_s[0] + acc->A[1][2] * in_s[1] + acc->A[2][2] * in_s[2] + acc->origin[2];
}

static inline void fract_to_cart_float(float out_x[3], const float in_s[3], const md_spatial_acc_t* acc) {
    // A is indexed as A[col][row]
    // cart = A * fract + origin
    out_x[0] = acc->A[0][0] * in_s[0] + acc->A[1][0] * in_s[1] + acc->A[2][0] * in_s[2] + acc->origin[0];
    out_x[1] = acc->A[0][1] * in_s[0] + acc->A[1][1] * in_s[1] + acc->A[2][1] * in_s[2] + acc->origin[1];
    out_x[2] = acc->A[0][2] * in_s[0] + acc->A[1][2] * in_s[1] + acc->A[2][2] * in_s[2] + acc->origin[2];
}

static inline void fract_to_cart_ort_256(
    md_256* out_cx, md_256* out_cy, md_256* out_cz,
    const md_256 in_sx, const md_256 in_sy, const md_256 in_sz,
    const md_256 A00, const md_256 A11, const md_256 A22,
    const md_256 O0,  const md_256 O1,  const md_256 O2)
{
    *out_cx = md_mm256_fmadd_ps(in_sx, A00, O0);
    *out_cy = md_mm256_fmadd_ps(in_sy, A11, O1);
    *out_cz = md_mm256_fmadd_ps(in_sz, A22, O2);
}

static inline void fract_to_cart_tri_256(
    md_256* out_cx, md_256* out_cy, md_256* out_cz,
    const md_256 in_sx, const md_256 in_sy, const md_256 in_sz,
    const md_256 A00, const md_256 A10, const md_256 A11, const md_256 A20, const md_256 A21, const md_256 A22,
    const md_256 O0,  const md_256 O1,  const md_256 O2)
{
    *out_cx = md_mm256_fmadd_ps(in_sx, A00, md_mm256_fmadd_ps(in_sy, A10, md_mm256_fmadd_ps(in_sz, A20, O0)));
    *out_cy = md_mm256_fmadd_ps(in_sy, A11, md_mm256_fmadd_ps(in_sz, A21, O1));
    *out_cz = md_mm256_fmadd_ps(in_sz, A22, O2);
}

static inline void cart_to_fract_ort_256(
    md_256* out_sx, md_256* out_sy, md_256* out_sz,
    const md_256 in_cx, const md_256 in_cy, const md_256 in_cz,
    const md_256 I00, const md_256 I11, const md_256 I22,
    const md_256 O0,  const md_256 O1,  const md_256 O2)
{
    // fract = I * (cart - origin)
    *out_sx = md_mm256_mul_ps(md_mm256_sub_ps(in_cx, O0), I00);
    *out_sy = md_mm256_mul_ps(md_mm256_sub_ps(in_cy, O1), I11);
    *out_sz = md_mm256_mul_ps(md_mm256_sub_ps(in_cz, O2), I22);
}

static inline void cart_to_fract_tri_256(
    md_256* out_sx, md_256* out_sy, md_256* out_sz,
    const md_256 in_cx, const md_256 in_cy, const md_256 in_cz,
    const md_256 I00, const md_256 I01, const md_256 I02,
    const md_256 I10, const md_256 I11, const md_256 I12,
    const md_256 I20, const md_256 I21, const md_256 I22,
    const md_256 O0,  const md_256 O1,  const md_256 O2)
{
    // fract = I * (cart - origin)
    const md_256 cx = md_mm256_sub_ps(in_cx, O0);
    const md_256 cy = md_mm256_sub_ps(in_cy, O1);
    const md_256 cz = md_mm256_sub_ps(in_cz, O2);
    *out_sx = md_mm256_fmadd_ps(cx, I00, md_mm256_fmadd_ps(cy, I01, md_mm256_mul_ps(cz, I02)));
    *out_sy = md_mm256_fmadd_ps(cx, I10, md_mm256_fmadd_ps(cy, I11, md_mm256_mul_ps(cz, I12)));
    *out_sz = md_mm256_fmadd_ps(cx, I20, md_mm256_fmadd_ps(cy, I21, md_mm256_mul_ps(cz, I22)));
}

static inline void batch_fract_to_cart_ort_256(float* x, float* y, float* z, size_t count, const md_spatial_acc_t* acc) {
	const md_256 A00 = md_mm256_set1_ps(acc->A[0][0]);
	const md_256 A11 = md_mm256_set1_ps(acc->A[1][1]);
	const md_256 A22 = md_mm256_set1_ps(acc->A[2][2]);

	const md_256 O0  = md_mm256_set1_ps(acc->origin[0]);
	const md_256 O1  = md_mm256_set1_ps(acc->origin[1]);
	const md_256 O2  = md_mm256_set1_ps(acc->origin[2]);

    for (size_t i = 0; i < count; i += 8) {
		md_256 sx = md_mm256_loadu_ps(x + i);
        md_256 sy = md_mm256_loadu_ps(y + i);
        md_256 sz = md_mm256_loadu_ps(z + i);

		md_256 cx, cy, cz;
        fract_to_cart_ort_256(&cx, &cy, &cz, sx, sy, sz, A00, A11, A22, O0, O1, O2);

		md_mm256_storeu_ps(x + i, cx);
		md_mm256_storeu_ps(y + i, cy);
		md_mm256_storeu_ps(z + i, cz);
    }
}

static inline void batch_fract_to_cart_tri_256(float* x, float* y, float* z, size_t count, const md_spatial_acc_t* acc) {
    const md_256 A00 = md_mm256_set1_ps(acc->A[0][0]);
    const md_256 A10 = md_mm256_set1_ps(acc->A[1][0]);
    const md_256 A11 = md_mm256_set1_ps(acc->A[1][1]);
    const md_256 A20 = md_mm256_set1_ps(acc->A[2][0]);
    const md_256 A21 = md_mm256_set1_ps(acc->A[2][1]);
    const md_256 A22 = md_mm256_set1_ps(acc->A[2][2]);
    const md_256 O0  = md_mm256_set1_ps(acc->origin[0]);
    const md_256 O1  = md_mm256_set1_ps(acc->origin[1]);
    const md_256 O2  = md_mm256_set1_ps(acc->origin[2]);

    for (size_t i = 0; i < count; i += 8) {
        md_256 sx = md_mm256_loadu_ps(x + i);
        md_256 sy = md_mm256_loadu_ps(y + i);
        md_256 sz = md_mm256_loadu_ps(z + i);

        md_256 cx, cy, cz;
        fract_to_cart_tri_256(&cx, &cy, &cz, sx, sy, sz, A00, A10, A11, A20, A21, A22, O0, O1, O2);

        md_mm256_storeu_ps(x + i, cx);
        md_mm256_storeu_ps(y + i, cy);
        md_mm256_storeu_ps(z + i, cz);
    }
}

// Macros for cell offsets/lengths in the sorted array
#define CELL_INDEX(x, y, z) ((size_t)z * c01 + (size_t)y * c0 + (size_t)x)
#define CELL_OFFSET(ci) (acc->cell_off[(ci)])
#define CELL_LENGTH(ci) (acc->cell_off[(ci) + 1] - acc->cell_off[(ci)])

#define POSSIBLY_INVOKE_CALLBACK_PAIR(estimated_count) \
    if (count + (estimated_count) >= SPATIAL_ACC_BUFLEN) { \
        callback(buf_i, buf_j, buf_d2, count, user_param); \
        count = 0; \
    } \

#define FLUSH_TAIL_PAIR() \
    if (count) { \
        callback(buf_i, buf_j, buf_d2, count, user_param); \
        count = 0; \
    } \

// The default convention for the callback is that we supply the coordinates as cartesian, not fractional.
// Therefore we perform the conversion here
#define POSSIBLY_INVOKE_CALLBACK_POINT_ORT(estimated_count) \
    if (count + (estimated_count) >= SPATIAL_ACC_BUFLEN) { \
        batch_fract_to_cart_ort_256(buf_x, buf_y, buf_z, count, acc); \
        callback(buf_i, buf_x, buf_y, buf_z, count, user_param); \
        count = 0; \
    } \

#define POSSIBLY_INVOKE_CALLBACK_POINT_FRACT_TRI(estimated_count) \
    if (count + (estimated_count) >= SPATIAL_ACC_BUFLEN) { \
        batch_fract_to_cart_tri_256(buf_x, buf_y, buf_z, count, acc); \
        callback(buf_i, buf_x, buf_y, buf_z, count, user_param); \
        count = 0; \
    } \

#define POSSIBLY_INVOKE_CALLBACK_POINT_CART_TRI(estimated_count) \
    if (count + (estimated_count) >= SPATIAL_ACC_BUFLEN) { \
        callback(buf_i, buf_x, buf_y, buf_z, count, user_param); \
        count = 0; \
    } \

#define FLUSH_TAIL_POINT_ORT() \
    if (count) { \
        batch_fract_to_cart_ort_256(buf_x, buf_y, buf_z, count, acc); \
        callback(buf_i, buf_x, buf_y, buf_z, count, user_param); \
        count = 0; \
    } \

#define FLUSH_TAIL_POINT_FRACT_TRI() \
    if (count) { \
        batch_fract_to_cart_tri_256(buf_x, buf_y, buf_z, count, acc); \
        callback(buf_i, buf_x, buf_y, buf_z, count, user_param); \
        count = 0; \
    } \

#define FLUSH_TAIL_POINT_CART_TRI() \
    if (count) { \
        callback(buf_i, buf_x, buf_y, buf_z, count, user_param); \
        count = 0; \
    } \

#define TEST_CELL_MASK(mask, coord) ((mask[(coord) / 64] & (1ULL << ((coord) % 64))) != 0)

static void for_each_internal_pair_in_neighboring_cells_triclinic(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param) {
    // Constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);

    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);

    // Intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    uint32_t buf_j[SPATIAL_ACC_BUFLEN];
    float    buf_d2[SPATIAL_ACC_BUFLEN];
	size_t count = 0;

    // Generate forward neighbor offsets
    const int ncell[3] = {1, 1, 1};
    int fwd_nbrs[16][4];
    size_t num_fwd_nbrs = generate_forward_neighbors4(fwd_nbrs, ncell);

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    // If we have a triclinic cell, the cell must be periodic in all dimensions.
    ASSERT(acc->flags & MD_UNITCELL_TRICLINIC);
    ASSERT((acc->flags & MD_UNITCELL_PBC_ALL) == (MD_UNITCELL_PBC_ALL));

    // --- Main cell loops ---
    for (uint32_t cz = 0; cz < cdim[2]; ++cz) {
		if (!TEST_CELL_MASK(acc->cell_mask[2], cz)) continue;  // skip empty z-slices
        for (uint32_t cy = 0; cy < cdim[1]; ++cy) {
			if (!TEST_CELL_MASK(acc->cell_mask[1], cy)) continue;  // skip empty y-slices
            for (uint32_t cx = 0; cx < cdim[0]; ++cx) {
				if (!TEST_CELL_MASK(acc->cell_mask[0], cx)) continue;  // skip empty x-slices

                const uint32_t ci    = CELL_INDEX(cx, cy, cz);
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);
                if (len_i == 0) continue;

                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_i - (i + 1), 8));

                    const md_256 v_xi    = md_mm256_set1_ps(elem_i_x[i]);
                    const md_256 v_yi    = md_mm256_set1_ps(elem_i_y[i]);
                    const md_256 v_zi    = md_mm256_set1_ps(elem_i_z[i]);
					const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                    md_256i v_j = md_mm256_add_epi32(md_mm256_set1_epi32(i + 1), inc);
                    for (uint32_t j = i + 1; j < len_i; j += 8) {
                        const md_256 v_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_i));

                        const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_i_x + j));
                        const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_i_y + j));
                        const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_i_z + j));
                        const md_256 v_d2 = distance_squared_tri_256(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

						const md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_i_idx + j));
						const int mask = md_mm256_movemask_ps(v_mask);

                        ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

						md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
						md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
						md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                        count += popcnt32(mask);

                        v_j = md_mm256_add_epi32(v_j, add8);
                    }
                }

                const ivec4_t c_v = ivec4_set(cx, cy, cz, 0);

                // --- Forward neighbors ---
                for (uint32_t n = 0; n < num_fwd_nbrs; ++n) {
                    const ivec4_t fwd_v = ivec4_load(fwd_nbrs[n]);
                    ivec4_t n_v = ivec4_add(c_v, fwd_v);

                    const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
                    const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);

                    // Apply wrapping
                    n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
                    n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

                    int n_arr[4];
                    ivec4_store(n_arr, n_v);

                    const uint32_t cj    = CELL_INDEX(n_arr[0], n_arr[1], n_arr[2]);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);
                    if (len_j == 0) continue;

                    // Compute shift vector (+1 for lower wrap, -1 for upper)
                    const ivec4_t shift_i = ivec4_sub(
                        ivec4_and(wrap_lower, ivec4_set1(1)),
                        ivec4_and(wrap_upper, ivec4_set1(1)));
                    const vec4_t shift_f = vec4_from_ivec4(shift_i);

                    const md_256 shift_x = md_mm256_set1_ps(shift_f.x);
                    const md_256 shift_y = md_mm256_set1_ps(shift_f.y);
                    const md_256 shift_z = md_mm256_set1_ps(shift_f.z);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_j, 8));

                        const md_256 v_xi = md_mm256_add_ps(md_mm256_set1_ps(elem_i_x[i]), shift_x);
                        const md_256 v_yi = md_mm256_add_ps(md_mm256_set1_ps(elem_i_y[i]), shift_y);
                        const md_256 v_zi = md_mm256_add_ps(md_mm256_set1_ps(elem_i_z[i]), shift_z);
						const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                        md_256i v_j = inc;
                        for (uint32_t j = 0; j < len_j; j += 8) {
                            const md_256  v_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));
                            const md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                            const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                            const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                            const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                            const md_256 v_d2 = distance_squared_tri_256(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                            const int mask = md_mm256_movemask_ps(v_mask);

                            ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                            md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
                            md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
                            md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                            count += popcnt32(mask);

                            v_j = md_mm256_add_epi32(v_j, add8);
                        }
                    }
                }
            }
        }
    }

    FLUSH_TAIL_PAIR();
}

static void for_each_internal_pair_in_neighboring_cells_ortho(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param) {

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    uint32_t buf_j[SPATIAL_ACC_BUFLEN];
    float    buf_d2[SPATIAL_ACC_BUFLEN];

    // Generate forward neighbor offsets
    const int ncell[3] = {1, 1, 1};
    int fwd_nbrs[16][4];
    size_t num_fwd_nbrs = generate_forward_neighbors4(fwd_nbrs, ncell);

    size_t count = 0;

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    const ivec4_t pmask_v = ivec4_set(
        (acc->flags & MD_UNITCELL_PBC_X) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Y) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Z) ? 0xFFFFFFFF : 0,
        0);

    // --- Main cell loops ---
    for (uint32_t cz = 0; cz < cdim[2]; ++cz) {
		if (!TEST_CELL_MASK(acc->cell_mask[2], cz)) continue;  // skip empty z-slices
        for (uint32_t cy = 0; cy < cdim[1]; ++cy) {
			if (!TEST_CELL_MASK(acc->cell_mask[1], cy)) continue;  // skip empty y-slices
            for (uint32_t cx = 0; cx < cdim[0]; ++cx) {
				if (!TEST_CELL_MASK(acc->cell_mask[0], cx)) continue;  // skip empty x-slices

                const uint32_t ci    = CELL_INDEX(cx, cy, cz);
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);
                if (len_i == 0) continue;

                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_i - (i + 1), 8));

                    const md_256 v_xi    = md_mm256_set1_ps(elem_i_x[i]);
                    const md_256 v_yi    = md_mm256_set1_ps(elem_i_y[i]);
                    const md_256 v_zi    = md_mm256_set1_ps(elem_i_z[i]);
					const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                    md_256i v_j = md_mm256_add_epi32(md_mm256_set1_epi32(i + 1), inc);
                    for (uint32_t j = i + 1; j < len_i; j += 8) {
                        const md_256 v_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_i));

                        const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_i_x + j));
                        const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_i_y + j));
                        const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_i_z + j));
                        const md_256 v_d2 = distance_squared_ort_256(v_dx, v_dy, v_dz, G00, G11, G22);

						const md_256i v_idxj = md_mm256_loadu_epi32(elem_i_idx + j);
						const int mask = md_mm256_movemask_ps(v_mask);

                        ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

						md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
						md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
						md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                        count += popcnt32(mask);

                        v_j = md_mm256_add_epi32(v_j, add8);
                    }
                }

                const ivec4_t c_v = ivec4_set(cx, cy, cz, 0);

                // --- Forward neighbors ---
                for (uint32_t n = 0; n < num_fwd_nbrs; ++n) {
                    const ivec4_t fwd_v = ivec4_load(fwd_nbrs[n]);
                    ivec4_t n_v = ivec4_add(c_v, fwd_v);

                    const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
                    const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);
                    const ivec4_t wrap_any   = ivec4_or(wrap_upper, wrap_lower);

                    // Skip nonperiodic wraps
                    if (ivec4_any(ivec4_andnot(wrap_any, pmask_v))) continue;

                    // Apply wrapping
                    n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
                    n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

                    int n_arr[4];
                    ivec4_store(n_arr, n_v);

                    const uint32_t cj    = CELL_INDEX(n_arr[0], n_arr[1], n_arr[2]);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);
                    if (len_j == 0) continue;

                    // Compute shift vector (+1 for lower wrap, -1 for upper)
                    const ivec4_t shift_i = ivec4_sub(
                        ivec4_and(wrap_lower, ivec4_set1(1)),
                        ivec4_and(wrap_upper, ivec4_set1(1)));
                    const vec4_t shift_f = vec4_from_ivec4(shift_i);

                    const md_256 shift_x = md_mm256_set1_ps(shift_f.x);
                    const md_256 shift_y = md_mm256_set1_ps(shift_f.y);
                    const md_256 shift_z = md_mm256_set1_ps(shift_f.z);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_j, 8));

                        const md_256 v_xi = md_mm256_add_ps(md_mm256_set1_ps(elem_i_x[i]), shift_x);
                        const md_256 v_yi = md_mm256_add_ps(md_mm256_set1_ps(elem_i_y[i]), shift_y);
                        const md_256 v_zi = md_mm256_add_ps(md_mm256_set1_ps(elem_i_z[i]), shift_z);
						const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                        md_256i v_j = inc;
                        for (uint32_t j = 0; j < len_j; j += 8) {
                            const md_256  v_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));
                            const md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                            const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                            const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                            const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                            const md_256 v_d2 = distance_squared_ort_256(v_dx, v_dy, v_dz, G00, G11, G22);

                            const int mask = md_mm256_movemask_ps(v_mask);

                            ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                            md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
                            md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
                            md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                            count += popcnt32(mask);

                            v_j = md_mm256_add_epi32(v_j, add8);
                        }
                    }
                }
            }
        }
    }

    FLUSH_TAIL_PAIR();
}

static void for_each_internal_pair_within_cutoff_triclinic(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    const int ncell[3] = {
        (int)ceil(cutoff * (double)acc->inv_cell_ext[0] * acc->cell_dim[0]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[1] * acc->cell_dim[1]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[2] * acc->cell_dim[2]),
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_pair_within_cutoff_ortho: cutoff too large for cell size");
        return;
    }

    float r2 = calc_r2(cutoff);

    // Initialize constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);
    const md_256 v_r2 = md_mm256_set1_ps(r2);

    // Support up to N-cell neighbor searches, optimal is 1-cell
    int neighbors[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];
    size_t num_neighbors = generate_forward_neighbors4(neighbors, ncell);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    uint32_t buf_j[SPATIAL_ACC_BUFLEN];
    float    buf_d2[SPATIAL_ACC_BUFLEN];

    size_t   count = 0;

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    // If we have a triclinic cell, the cell must be periodic in all dimensions.
    ASSERT(acc->flags & MD_UNITCELL_TRICLINIC);
    ASSERT((acc->flags & MD_UNITCELL_PBC_ALL) == (MD_UNITCELL_PBC_ALL));

    // --- Main cell loops ---
    for (uint32_t cz = 0; cz < cdim[2]; ++cz) {
		if (!TEST_CELL_MASK(acc->cell_mask[2], cz)) continue;  // skip empty z-slices
        for (uint32_t cy = 0; cy < cdim[1]; ++cy) {
			if (!TEST_CELL_MASK(acc->cell_mask[1], cy)) continue;  // skip empty y-slices
            for (uint32_t cx = 0; cx < cdim[0]; ++cx) {
				if (!TEST_CELL_MASK(acc->cell_mask[0], cx)) continue;  // skip empty x-slices

                const uint32_t ci    = CELL_INDEX(cx, cy, cz);
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);
                if (len_i == 0) continue;

                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_i - (i + 1), 8));

                    const md_256 v_xi    = md_mm256_set1_ps(elem_i_x[i]);
                    const md_256 v_yi    = md_mm256_set1_ps(elem_i_y[i]);
                    const md_256 v_zi    = md_mm256_set1_ps(elem_i_z[i]);
					const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                    md_256i v_j = md_mm256_add_epi32(md_mm256_set1_epi32(i + 1), inc);
                    for (uint32_t j = i + 1; j < len_i; j += 8) {
                        const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_i));

                        const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_i_x + j));
                        const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_i_y + j));
                        const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_i_z + j));
                        md_256 v_d2 = distance_squared_tri_256(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                        const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

                        // Fill buffers with results
                        const int mask = md_mm256_movemask_ps(v_mask);
                        if (mask) {
                            const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
						    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_i_idx + j));

						    v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
						    v_d2   = md_mm256_permutevar8x32_ps(v_d2,      v_idx_mask);

                            ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

						    md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
						    md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
						    md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                            count += popcnt32(mask);
                        }

                        v_j = md_mm256_add_epi32(v_j, add8);
                    }
                }

                const ivec4_t c_v = ivec4_set(cx, cy, cz, 0);

                // --- Forward neighbors ---
                for (uint32_t n = 0; n < num_neighbors; ++n) {
                    const ivec4_t fwd_v = ivec4_load(neighbors[n]);
                    ivec4_t n_v = ivec4_add(c_v, fwd_v);

                    const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
                    const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);

                    // Apply wrapping
                    n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
                    n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

                    int n_arr[4];
                    ivec4_store(n_arr, n_v);

                    const uint32_t cj    = CELL_INDEX(n_arr[0], n_arr[1], n_arr[2]);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);
                    if (len_j == 0) continue;

                    // Compute shift vector (+1 for lower wrap, -1 for upper)
                    const ivec4_t shift_i = ivec4_sub(
                        ivec4_and(wrap_lower, ivec4_set1(1)),
                        ivec4_and(wrap_upper, ivec4_set1(1)));
                    const vec4_t shift_f = vec4_from_ivec4(shift_i);

                    const md_256 shift_x = md_mm256_set1_ps(shift_f.x);
                    const md_256 shift_y = md_mm256_set1_ps(shift_f.y);
                    const md_256 shift_z = md_mm256_set1_ps(shift_f.z);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_j, 8));

                        const md_256 v_xi    = md_mm256_add_ps(md_mm256_set1_ps(elem_i_x[i]), shift_x);
                        const md_256 v_yi    = md_mm256_add_ps(md_mm256_set1_ps(elem_i_y[i]), shift_y);
                        const md_256 v_zi    = md_mm256_add_ps(md_mm256_set1_ps(elem_i_z[i]), shift_z);
						const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                        md_256i v_j = inc;
                        for (uint32_t j = 0; j < len_j; j += 8) {
                            const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));

                            const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                            const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                            const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                            md_256 v_d2 = distance_squared_tri_256(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                            const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

                            // Fill buffers with results
                            const int mask = md_mm256_movemask_ps(v_mask);
                            if (mask) {
                                const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                                md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                                v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                                v_d2   = md_mm256_permutevar8x32_ps(v_d2,      v_idx_mask);

                                ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                                md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
                                md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
                                md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                                count += popcnt32(mask);
                            }

                            v_j = md_mm256_add_epi32(v_j, add8);
                        }
                    }
                }
            }
        }
    }

    FLUSH_TAIL_PAIR();
}

static void for_each_internal_pair_within_cutoff_ortho(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    const int ncell[3] = {
        (int)ceil(cutoff * (double)acc->inv_cell_ext[0] * acc->cell_dim[0]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[1] * acc->cell_dim[1]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[2] * acc->cell_dim[2]),
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_pair_within_cutoff_ortho: cutoff too large for cell size");
        return;
    }

    const float r2 = calc_r2(cutoff);

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
	const md_256 v_r2 = md_mm256_set1_ps(r2);

	// Support up to 2-cell neighbor searches, optimal is 1-cell
    int neighbors[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];

    size_t num_neighbors = generate_forward_neighbors4(neighbors, ncell);

	// Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
	uint32_t buf_j[SPATIAL_ACC_BUFLEN];
	float    buf_d2[SPATIAL_ACC_BUFLEN];

    size_t   count = 0;

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    const ivec4_t pmask_v = ivec4_set(
        (acc->flags & MD_UNITCELL_PBC_X) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Y) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Z) ? 0xFFFFFFFF : 0,
        0);

    // --- Main cell loops ---
    for (uint32_t cz = 0; cz < cdim[2]; ++cz) {
		if (!TEST_CELL_MASK(acc->cell_mask[2], cz)) continue;  // skip empty z-slices
        for (uint32_t cy = 0; cy < cdim[1]; ++cy) {
			if (!TEST_CELL_MASK(acc->cell_mask[1], cy)) continue;  // skip empty y-slices
            for (uint32_t cx = 0; cx < cdim[0]; ++cx) {
				if (!TEST_CELL_MASK(acc->cell_mask[0], cx)) continue;  // skip empty x-slices

                const uint32_t ci    = CELL_INDEX(cx, cy, cz);
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);
                if (len_i == 0) continue;

                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_i - (i + 1), 8));

                    const md_256 v_xi    = md_mm256_set1_ps(elem_i_x[i]);
                    const md_256 v_yi    = md_mm256_set1_ps(elem_i_y[i]);
                    const md_256 v_zi    = md_mm256_set1_ps(elem_i_z[i]);
					const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                    md_256i v_j = md_mm256_add_epi32(md_mm256_set1_epi32(i + 1), inc);
                    for (uint32_t j = i + 1; j < len_i; j += 8) {

                        const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_i));

                        const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_i_x + j));
                        const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_i_y + j));
                        const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_i_z + j));
                        md_256 v_d2 = distance_squared_ort_256(v_dx, v_dy, v_dz, G00, G11, G22);

                        const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

						// Fill buffers with results
						const int mask = md_mm256_movemask_ps(v_mask);
                        if (mask) {
                            const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
						    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_i_idx + j));

						    v_idxj = md_mm256_permutevar8x32_epi32(v_idxj,  v_idx_mask);
						    v_d2   = md_mm256_permutevar8x32_ps(v_d2,       v_idx_mask);

                            ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

						    md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
						    md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
						    md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                            count += popcnt32(mask);
                        }

                        v_j = md_mm256_add_epi32(v_j, add8);
                    }
                }

                const ivec4_t c_v = ivec4_set(cx, cy, cz, 0);

                // --- Forward neighbors ---
                for (uint32_t n = 0; n < num_neighbors; ++n) {
                    const ivec4_t fwd_v = ivec4_load(neighbors[n]);
                    ivec4_t n_v = ivec4_add(c_v, fwd_v);

                    const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
                    const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);
                    const ivec4_t wrap_any   = ivec4_or(wrap_upper, wrap_lower);

                    // Skip nonperiodic wraps
                    if (ivec4_any(ivec4_andnot(wrap_any, pmask_v))) continue;

                    // Apply wrapping
                    n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
                    n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

                    int n_arr[4];
                    ivec4_store(n_arr, n_v);

                    const uint32_t cj    = CELL_INDEX(n_arr[0], n_arr[1], n_arr[2]);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);
                    if (len_j == 0) continue;

                    // Compute shift vector (+1 for lower wrap, -1 for upper)
                    const ivec4_t shift_i = ivec4_sub(
                        ivec4_and(wrap_lower, ivec4_set1(1)),
                        ivec4_and(wrap_upper, ivec4_set1(1)));
                    const vec4_t shift_f = vec4_from_ivec4(shift_i);

                    const md_256 shift_x = md_mm256_set1_ps(shift_f.x);
                    const md_256 shift_y = md_mm256_set1_ps(shift_f.y);
                    const md_256 shift_z = md_mm256_set1_ps(shift_f.z);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_j, 8));

                        const md_256 v_xi    = md_mm256_add_ps(md_mm256_set1_ps(elem_i_x[i]), shift_x);
                        const md_256 v_yi    = md_mm256_add_ps(md_mm256_set1_ps(elem_i_y[i]), shift_y);
                        const md_256 v_zi    = md_mm256_add_ps(md_mm256_set1_ps(elem_i_z[i]), shift_z);
						const md_256i v_idxi = md_mm256_set1_epi32(elem_i_idx[i]);

                        md_256i v_j = inc;
                        for (uint32_t j = 0; j < len_j; j += 8) {
                            const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));

                            const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                            const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                            const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                            md_256 v_d2 = distance_squared_ort_256(v_dx, v_dy, v_dz, G00, G11, G22);

                            const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

                            // Fill buffers with results
                            const int mask = md_mm256_movemask_ps(v_mask);
                            if (mask) {
                                const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                                md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                                v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                                v_d2   = md_mm256_permutevar8x32_ps(v_d2,      v_idx_mask);

                                ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                                md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
                                md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
                                md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                                count += popcnt32(mask);
                            }

                            v_j = md_mm256_add_epi32(v_j, add8);
                        }
                    }
                }
            }
        }
    }

    FLUSH_TAIL_PAIR();
}

static void for_each_external_pair_within_cutoff_triclinic(const md_spatial_acc_t* acc, const md_coord_stream_t* ext_stream, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags) {
    const int ncell[3] = {
        (int)ceil(cutoff * (double)acc->inv_cell_ext[0] * acc->cell_dim[0]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[1] * acc->cell_dim[1]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[2] * acc->cell_dim[2]),
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_external_pair_within_cutoff_ortho: cutoff too large for cell size");
        return;
    }

    if (flags & MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX) {
        if (!ext_stream->idx) {
            MD_LOG_ERROR("for_each_external_pair_within_cutoff_ortho: MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX flag is set but ext_stream->idx is NULL");
            return;
        }
    }

    const float r2 = calc_r2(cutoff);

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);
    const md_256 v_r2 = md_mm256_set1_ps(r2);

    // For storing neighbor cell coordinates around a single external point
    int neighbors[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];

    size_t num_neighbors = generate_neighbors4(neighbors, ncell);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    uint32_t buf_j[SPATIAL_ACC_BUFLEN];
    float    buf_d2[SPATIAL_ACC_BUFLEN];

    size_t   count = 0;

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    // If we have a triclinic cell, the cell must be periodic in all dimensions.
    ASSERT(acc->flags & MD_UNITCELL_TRICLINIC);
    ASSERT((acc->flags & MD_UNITCELL_PBC_ALL) == (MD_UNITCELL_PBC_ALL));

	const vec4_t cell_dim = vec4_set((float)cdim[0], (float)cdim[1], (float)cdim[2], 0);

    for (size_t ei = 0; ei < ext_stream->count; ++ei) {
        uint32_t idx = md_coord_stream_load_idx(ext_stream, ei);
        vec4_t r     = md_coord_stream_load_vec4(ext_stream, ei);

        vec4_t f = vec4_cart_to_fract(r, acc);
		ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f, cell_dim)));

        uint32_t idx_i = (flags & MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX) ? idx : (uint32_t)ei;
        const md_256i v_idxi = md_mm256_set1_epi32(idx_i);

        for (uint32_t n = 0; n < num_neighbors; ++n) {
            const ivec4_t fwd_v = ivec4_load(neighbors[n]);
            ivec4_t n_v = ivec4_add(c_v, fwd_v);

            const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
            const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);

            // Apply wrapping
            n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
            n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

            int n_arr[4];
            ivec4_store(n_arr, n_v);

            const uint32_t cj    = CELL_INDEX((uint32_t)n_arr[0], (uint32_t)n_arr[1], (uint32_t)n_arr[2]);
            const uint32_t off_j = CELL_OFFSET(cj);
            const uint32_t len_j = CELL_LENGTH(cj);
            if (len_j == 0) continue;

            const md_256i v_len_j = md_mm256_set1_epi32(len_j);

            // Compute shift vector (+1 for lower wrap, -1 for upper)
            const ivec4_t shift_i = ivec4_sub(
                ivec4_and(wrap_lower, ivec4_set1(1)),
                ivec4_and(wrap_upper, ivec4_set1(1))
            );

            // Shift external point by periodic image offset
			const vec4_t f_shift = vec4_add(f, vec4_from_ivec4(shift_i));

            const md_256 v_xi = md_mm256_set1_ps(f_shift.x);
            const md_256 v_yi = md_mm256_set1_ps(f_shift.y);
            const md_256 v_zi = md_mm256_set1_ps(f_shift.z);

            const float* elem_j_x = element_x + off_j;
            const float* elem_j_y = element_y + off_j;
            const float* elem_j_z = element_z + off_j;
            const uint32_t* elem_j_idx = element_i + off_j;

            POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_j, 8));

            md_256i v_j = inc;
            for (uint32_t j = 0; j < len_j; j += 8) {
                const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));

                const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                md_256 v_d2 = distance_squared_tri_256(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

                // Fill buffers with results
                const int mask = md_mm256_movemask_ps(v_mask);
                if (mask) {
                    const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                    v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                    v_d2   = md_mm256_permutevar8x32_ps(v_d2,      v_idx_mask);

                    ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                    md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
                    md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
                    md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                    count += popcnt32(mask);
                }

                v_j = md_mm256_add_epi32(v_j, add8);
            }
        }
    }

    FLUSH_TAIL_PAIR();
}

static void for_each_external_pair_within_cutoff_ortho(const md_spatial_acc_t* acc, const md_coord_stream_t* ext_stream, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags) {
    const int ncell[3] = {
        (int)ceil(cutoff * (double)acc->inv_cell_ext[0] * acc->cell_dim[0]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[1] * acc->cell_dim[1]),
        (int)ceil(cutoff * (double)acc->inv_cell_ext[2] * acc->cell_dim[2]),
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_external_pair_within_cutoff_ortho: cutoff too large for cell size");
        return;
    }

    if (flags & MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX) {
        if (!ext_stream->idx) {
            MD_LOG_ERROR("for_each_external_pair_within_cutoff_ortho: MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX flag is set but ext_stream->idx is NULL");
            return;
        }
    }

    const float r2 = calc_r2(cutoff);

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 v_r2 = md_mm256_set1_ps(r2);

    // For storing neighbor cell coordinates around a single external point
    int neighbors[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];

    size_t num_neighbors = generate_neighbors4(neighbors, ncell);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    uint32_t buf_j[SPATIAL_ACC_BUFLEN];
    float    buf_d2[SPATIAL_ACC_BUFLEN];

    size_t   count = 0;

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = cdim[0];
    const uint32_t c01 = cdim[0] * cdim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    const ivec4_t pmask_v = ivec4_set(
        (acc->flags & MD_UNITCELL_PBC_X) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Y) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Z) ? 0xFFFFFFFF : 0,
        0);

    vec4_t fract_mask;
    MEMCPY(&fract_mask, &pmask_v, sizeof(ivec4_t));  // avoid strict aliasing issues
	const vec4_t fcell_dim = vec4_set((float)cdim[0], (float)cdim[1], (float)cdim[2], 0);

    for (size_t ei = 0; ei < ext_stream->count; ++ei) {
        uint32_t idx = md_coord_stream_load_idx(ext_stream, ei);
        vec4_t r     = md_coord_stream_load_vec4(ext_stream, ei);

		vec4_t f = vec4_cart_to_fract(r, acc);
		f = vec4_blend(f, vec4_fract(f), fract_mask);
		ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f, fcell_dim)));

        uint32_t idx_i = (flags & MD_SPATIAL_ACC_FLAG_USE_SUPPLIED_IDX) ? idx : (uint32_t)ei;
        const md_256i v_idxi = md_mm256_set1_epi32(idx_i);

        for (uint32_t n = 0; n < num_neighbors; ++n) {
            const ivec4_t fwd_v = ivec4_load(neighbors[n]);
            ivec4_t n_v = ivec4_add(c_v, fwd_v);

            const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
            const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);
            const ivec4_t wrap_any   = ivec4_or(wrap_upper, wrap_lower);

            // Skip nonperiodic wraps
            if (ivec4_any(ivec4_andnot(wrap_any, pmask_v))) continue;

            // Apply wrapping
            n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
            n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));

            int n_arr[4];
            ivec4_store(n_arr, n_v);

            const uint32_t cj    = CELL_INDEX((uint32_t)n_arr[0], (uint32_t)n_arr[1], (uint32_t)n_arr[2]);
            const uint32_t off_j = CELL_OFFSET(cj);
            const uint32_t len_j = CELL_LENGTH(cj);
            if (len_j == 0) continue;

            const md_256i v_len_j = md_mm256_set1_epi32(len_j);

            // Compute shift vector (+1 for lower wrap, -1 for upper)
            const ivec4_t shift_i = ivec4_sub(
                ivec4_and(wrap_lower, ivec4_set1(1)),
                ivec4_and(wrap_upper, ivec4_set1(1)));

            // Shift external point by periodic image offset
			const vec4_t f_shift = vec4_add(f, vec4_from_ivec4(shift_i));

            const md_256 v_xi = md_mm256_set1_ps(f_shift.x);
            const md_256 v_yi = md_mm256_set1_ps(f_shift.y);
            const md_256 v_zi = md_mm256_set1_ps(f_shift.z);

            const float* elem_j_x = element_x + off_j;
            const float* elem_j_y = element_y + off_j;
            const float* elem_j_z = element_z + off_j;
            const uint32_t* elem_j_idx = element_i + off_j;

            POSSIBLY_INVOKE_CALLBACK_PAIR(ALIGN_TO(len_j, 8));

            md_256i v_j = inc;
            for (uint32_t j = 0; j < len_j; j += 8) {
                const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));

                const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                md_256 v_d2 = distance_squared_ort_256(v_dx, v_dy, v_dz, G00, G11, G22);

                const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

                // Fill buffers with results
                const int mask = md_mm256_movemask_ps(v_mask);
                if (mask) {
                    const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                    v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                    v_d2   = md_mm256_permutevar8x32_ps(v_d2,      v_idx_mask);

                    ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                    md_mm256_storeu_epi32(buf_i  + count,  v_idxi);
                    md_mm256_storeu_epi32(buf_j  + count,  v_idxj);
                    md_mm256_storeu_ps   (buf_d2 + count,  v_d2);

                    count += popcnt32(mask);
                }

                v_j = md_mm256_add_epi32(v_j, add8);
            }
        }
    }

    FLUSH_TAIL_PAIR();
}

static inline void cell_range_from_aabb_center_radius(
    int out_cmin[3],                 // inclusive
    int out_cmax[3],                 // exclusive (loop while ic < cmax)
    double out_fcen[3],              // fractional center (wrapped to [0,1) on periodic axes)
    double out_frad[3],              // fractional half-extents
    const double center[3],          // cartesian center
    const double radius[3],          // cartesian half-extents
    const md_spatial_acc_t* acc)
{
    const int pbc[3] = {
        (acc->flags & MD_UNITCELL_PBC_X) != 0,
        (acc->flags & MD_UNITCELL_PBC_Y) != 0,
        (acc->flags & MD_UNITCELL_PBC_Z) != 0,
    };
    
    double sc[3];
    cart_to_fract(sc, center, acc);

    for (int a = 0; a < 3; ++a) {
        if (pbc[a]) sc[a] = fract(sc[a]);
    }

    double cc[3];
    fract_to_cart(cc, sc, acc);

    // 2) Convert 8 corners to fractional, unwrap them near the center image, take component-wise min/max.
    double fmin[3] = { +DBL_MAX, +DBL_MAX, +DBL_MAX };
    double fmax[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };

    for (int iz = 0; iz < 2; ++iz) {
        const double z = cc[2] + (iz ? +radius[2] : -radius[2]);
        for (int iy = 0; iy < 2; ++iy) {
            const double y = cc[1] + (iy ? +radius[1] : -radius[1]);
            for (int ix = 0; ix < 2; ++ix) {
                const double x = cc[0] + (ix ? +radius[0] : -radius[0]);
                const double c[3] = { x, y, z };

                double s[3];
                cart_to_fract(s, c, acc);

                for (int a = 0; a < 3; ++a) {
                    fmin[a] = MIN(fmin[a], s[a]);
                    fmax[a] = MAX(fmax[a], s[a]);
                }
            }
        }
    }

    out_fcen[0] = sc[0];
    out_fcen[1] = sc[1];
    out_fcen[2] = sc[2];

    out_frad[0] = 0.5 * (fmax[0] - fmin[0]);
    out_frad[1] = 0.5 * (fmax[1] - fmin[1]);
    out_frad[2] = 0.5 * (fmax[2] - fmin[2]);

    // 3) Convert fractional bounds to cell index bounds in the (fractional) grid.
    for (int a = 0; a < 3; ++a) {
        const int dim = (int)acc->cell_dim[a];

        int cmin = (int)floor(fmin[a] * (double)dim);
        int cmax = (int)ceil (fmax[a] * (double)dim);

        if (cmax <= cmin) cmax = cmin + 1;  // ensure at least one cell

        if (!pbc[a]) {
            cmin = CLAMP(cmin, 0, dim);
            cmax = CLAMP(cmax, 0, dim);
            if (cmax <= cmin) cmax = MIN(cmin + 1, dim);
        }

        out_cmin[a] = cmin;
        out_cmax[a] = cmax;
    }
}
static void for_each_point_in_aabb_ortho(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);
    ASSERT((acc->flags & MD_UNITCELL_TRICLINIC) == 0);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    float    buf_x[SPATIAL_ACC_BUFLEN];
    float    buf_y[SPATIAL_ACC_BUFLEN];
    float    buf_z[SPATIAL_ACC_BUFLEN];

    size_t   count = 0;

    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    ASSERT(cdim[0] > 0 && cdim[1] > 0 && cdim[2] > 0);

    const uint32_t c0  = acc->cell_dim[0];
    const uint32_t c01 = acc->cell_dim[0] * acc->cell_dim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    const int pbc[3] = {
        (acc->flags & MD_UNITCELL_PBC_X) != 0,
        (acc->flags & MD_UNITCELL_PBC_Y) != 0,
        (acc->flags & MD_UNITCELL_PBC_Z) != 0,
    };

    int cmin[3], cmax[3];
    double fcen[3], frad[3];
    cell_range_from_aabb_center_radius(cmin, cmax, fcen, frad, aabb_cen, aabb_rad, acc);

    // Clamp the fractional cell extent
    for (int i = 0; i < 3; ++i) {
        frad[i] = MIN(frad[i], 0.5);
    }

    const md_256 v_fminx = md_mm256_set1_ps((float)(fcen[0] - frad[0]));
    const md_256 v_fminy = md_mm256_set1_ps((float)(fcen[1] - frad[1]));
    const md_256 v_fminz = md_mm256_set1_ps((float)(fcen[2] - frad[2]));

    const md_256 v_fmaxx = md_mm256_set1_ps((float)(fcen[0] + frad[0]));
    const md_256 v_fmaxy = md_mm256_set1_ps((float)(fcen[1] + frad[1]));
    const md_256 v_fmaxz = md_mm256_set1_ps((float)(fcen[2] + frad[2]));

    for (int icz = cmin[2]; icz < cmax[2]; ++icz) {
        int cz = pbc[2] ? wrap_coord(icz, (int)cdim[2]) : icz;
        int sz = isign(icz - cz);
        const md_256 shift_z = md_mm256_set1_ps((float)sz);

        for (int icy = cmin[1]; icy < cmax[1]; ++icy) {
            int cy = pbc[1] ? wrap_coord(icy, (int)cdim[1]) : icy;
            int sy = isign(icy - cy);
            const md_256 shift_y = md_mm256_set1_ps((float)sy);

            for (int icx = cmin[0]; icx < cmax[0]; ++icx) {
                int cx = pbc[0] ? wrap_coord(icx, (int)cdim[0]) : icx;
                int sx = isign(icx - cx);
                const md_256 shift_x = md_mm256_set1_ps((float)sx);

                const uint32_t ci  = CELL_INDEX((uint32_t)cx, (uint32_t)cy, (uint32_t)cz);
                const uint32_t off = CELL_OFFSET(ci);
                const uint32_t len = CELL_LENGTH(ci);
                if (len == 0) continue;

                const md_256i v_len = md_mm256_set1_epi32(len);

                const float* elem_x = acc->elem_x + off;
                const float* elem_y = acc->elem_y + off;
                const float* elem_z = acc->elem_z + off;
                const uint32_t* elem_idx = acc->elem_idx + off;

                POSSIBLY_INVOKE_CALLBACK_POINT_ORT(ALIGN_TO(len, 8));

                md_256i v_j = inc;
                for (uint32_t j = 0; j < len; j += 8) {
                    const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len));

                    md_256 v_xj = md_mm256_loadu_ps(elem_x + j);
                    md_256 v_yj = md_mm256_loadu_ps(elem_y + j);
                    md_256 v_zj = md_mm256_loadu_ps(elem_z + j);

                    v_xj = md_mm256_add_ps(v_xj, shift_x);
                    v_yj = md_mm256_add_ps(v_yj, shift_y);
                    v_zj = md_mm256_add_ps(v_zj, shift_z);

                    const md_256 v_mask_min = md_mm256_and_ps(
                        md_mm256_cmpge_ps(v_xj, v_fminx),
                        md_mm256_and_ps(
                            md_mm256_cmpge_ps(v_yj, v_fminy),
                            md_mm256_cmpge_ps(v_zj, v_fminz)));

                    const md_256 v_mask_max = md_mm256_and_ps(
                        md_mm256_cmple_ps(v_xj, v_fmaxx),
                        md_mm256_and_ps(
                            md_mm256_cmple_ps(v_yj, v_fmaxy),
                            md_mm256_cmple_ps(v_zj, v_fmaxz)));

                    const md_256 v_mask = md_mm256_and_ps(md_mm256_and_ps(v_mask_min, v_mask_max), j_mask);

                    const int mask = md_mm256_movemask_ps(v_mask);
                    if (mask) {
                        const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                        md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_idx + j));

                        v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                        v_xj   = md_mm256_permutevar8x32_ps(v_xj, v_idx_mask);
                        v_yj   = md_mm256_permutevar8x32_ps(v_yj, v_idx_mask);
                        v_zj   = md_mm256_permutevar8x32_ps(v_zj, v_idx_mask);

                        ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                        md_mm256_storeu_epi32(buf_i + count, v_idxj);
                        md_mm256_storeu_ps   (buf_x + count, v_xj);
                        md_mm256_storeu_ps   (buf_y + count, v_yj);
                        md_mm256_storeu_ps   (buf_z + count, v_zj);

                        count += popcnt32(mask);
                    }

                    v_j = md_mm256_add_epi32(v_j, add8);
                }
            }
        }
    }

    FLUSH_TAIL_POINT_ORT();
}

static void for_each_point_in_aabb_triclinic(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);
    ASSERT(acc->flags & MD_UNITCELL_TRICLINIC);
    ASSERT((acc->flags & MD_UNITCELL_PBC_ALL) == MD_UNITCELL_PBC_ALL);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    float    buf_x[SPATIAL_ACC_BUFLEN];
    float    buf_y[SPATIAL_ACC_BUFLEN];
    float    buf_z[SPATIAL_ACC_BUFLEN];
    size_t   count = 0;

    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    ASSERT(cdim[0] > 0 && cdim[1] > 0 && cdim[2] > 0);

    const uint32_t c0  = acc->cell_dim[0];
    const uint32_t c01 = acc->cell_dim[0] * acc->cell_dim[1];

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    int cmin[3], cmax[3];
    double fcen[3], frad[3];
    cell_range_from_aabb_center_radius(cmin, cmax, fcen, frad, aabb_cen, aabb_rad, acc);

    // Exact cartesian bounds (in the wrapped "center image")
    double cc[3];
    fract_to_cart(cc, fcen, acc);

    const md_256 v_cminx = md_mm256_set1_ps((float)(cc[0] - aabb_rad[0]));
    const md_256 v_cminy = md_mm256_set1_ps((float)(cc[1] - aabb_rad[1]));
    const md_256 v_cminz = md_mm256_set1_ps((float)(cc[2] - aabb_rad[2]));

    const md_256 v_cmaxx = md_mm256_set1_ps((float)(cc[0] + aabb_rad[0]));
    const md_256 v_cmaxy = md_mm256_set1_ps((float)(cc[1] + aabb_rad[1]));
    const md_256 v_cmaxz = md_mm256_set1_ps((float)(cc[2] + aabb_rad[2]));

    // Precompute frac->cart constants for triclinic conversion
    const md_256 A00 = md_mm256_set1_ps(acc->A[0][0]);
    const md_256 A10 = md_mm256_set1_ps(acc->A[1][0]);
    const md_256 A11 = md_mm256_set1_ps(acc->A[1][1]);
    const md_256 A20 = md_mm256_set1_ps(acc->A[2][0]);
    const md_256 A21 = md_mm256_set1_ps(acc->A[2][1]);
    const md_256 A22 = md_mm256_set1_ps(acc->A[2][2]);

    const md_256 O0  = md_mm256_set1_ps(acc->origin[0]);
    const md_256 O1  = md_mm256_set1_ps(acc->origin[1]);
    const md_256 O2  = md_mm256_set1_ps(acc->origin[2]);

    for (int icz = cmin[2]; icz < cmax[2]; ++icz) {
        int cz = wrap_coord(icz, (int)cdim[2]);
        int sz = isign(icz - cz);
        const md_256 shift_z = md_mm256_set1_ps((float)sz);

        for (int icy = cmin[1]; icy < cmax[1]; ++icy) {
            int cy = wrap_coord(icy, (int)cdim[1]);
            int sy = isign(icy - cy);
            const md_256 shift_y = md_mm256_set1_ps((float)sy);

            for (int icx = cmin[0]; icx < cmax[0]; ++icx) {
                int cx = wrap_coord(icx, (int)cdim[0]);
                int sx = isign(icx - cx);
                const md_256 shift_x = md_mm256_set1_ps((float)sx);

                const uint32_t ci  = CELL_INDEX((uint32_t)cx, (uint32_t)cy, (uint32_t)cz);
                const uint32_t off = CELL_OFFSET(ci);
                const uint32_t len = CELL_LENGTH(ci);
                if (len == 0) continue;

                const md_256i v_len = md_mm256_set1_epi32(len);

                const float* elem_x = acc->elem_x + off;
                const float* elem_y = acc->elem_y + off;
                const float* elem_z = acc->elem_z + off;
                const uint32_t* elem_idx = acc->elem_idx + off;

                POSSIBLY_INVOKE_CALLBACK_POINT_CART_TRI(ALIGN_TO(len, 8));

                md_256i v_j = inc;
                for (uint32_t j = 0; j < len; j += 8) {
                    const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len));

                    // Fractional coords (stored) + periodic image shift
                    md_256 v_sx = md_mm256_loadu_ps(elem_x + j);
                    md_256 v_sy = md_mm256_loadu_ps(elem_y + j);
                    md_256 v_sz = md_mm256_loadu_ps(elem_z + j);

                    v_sx = md_mm256_add_ps(v_sx, shift_x);
                    v_sy = md_mm256_add_ps(v_sy, shift_y);
                    v_sz = md_mm256_add_ps(v_sz, shift_z);

                    md_256 v_cx, v_cy, v_cz;
                    fract_to_cart_tri_256(&v_cx, &v_cy, &v_cz, v_sx, v_sy, v_sz, A00, A10, A11, A20, A21, A22, O0, O1, O2);

                    const md_256 v_cmask_min = md_mm256_and_ps(
                        md_mm256_cmpge_ps(v_cx, v_cminx),
                        md_mm256_and_ps(
                            md_mm256_cmpge_ps(v_cy, v_cminy),
                            md_mm256_cmpge_ps(v_cz, v_cminz)));

                    const md_256 v_cmask_max = md_mm256_and_ps(
                        md_mm256_cmple_ps(v_cx, v_cmaxx),
                        md_mm256_and_ps(
                            md_mm256_cmple_ps(v_cy, v_cmaxy),
                            md_mm256_cmple_ps(v_cz, v_cmaxz)));

                    const md_256 v_mask = md_mm256_and_ps(md_mm256_and_ps(v_cmask_min, v_cmask_max), j_mask);

                    const int mask = md_mm256_movemask_ps(v_mask);
                    if (mask) {
                        const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                        md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_idx + j));

                        v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                        v_sx   = md_mm256_permutevar8x32_ps(v_sx, v_idx_mask);
                        v_sy   = md_mm256_permutevar8x32_ps(v_sy, v_idx_mask);
                        v_sz   = md_mm256_permutevar8x32_ps(v_sz, v_idx_mask);

                        ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                        md_mm256_storeu_epi32(buf_i + count, v_idxj);
                        md_mm256_storeu_ps   (buf_x + count, v_sx);
                        md_mm256_storeu_ps   (buf_y + count, v_sy);
                        md_mm256_storeu_ps   (buf_z + count, v_sz);

                        count += popcnt32(mask);
                    }

                    v_j = md_mm256_add_epi32(v_j, add8);
                }
            }
        }
    }

    FLUSH_TAIL_POINT_CART_TRI();
}

static void for_each_point_in_sphere_ortho(const md_spatial_acc_t* acc, const double center[3], double radius, md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);

    const int ncell[3] = {
        (int)ceil(radius * (double)acc->inv_cell_ext[0] * acc->cell_dim[0]),
        (int)ceil(radius * (double)acc->inv_cell_ext[1] * acc->cell_dim[1]),
        (int)ceil(radius * (double)acc->inv_cell_ext[2] * acc->cell_dim[2]),
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_point_in_sphere_ortho: radius too large for cell size");
        return;
    }

    // Build neighbor offsets
    int neighbors[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];
    size_t num_neighbors = generate_neighbors4(neighbors, ncell);

    // Buffers
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    float    buf_x[SPATIAL_ACC_BUFLEN];
    float    buf_y[SPATIAL_ACC_BUFLEN];
    float    buf_z[SPATIAL_ACC_BUFLEN];
    size_t   count = 0;

    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = acc->cell_dim[0];
    const uint32_t c01 = acc->cell_dim[0] * acc->cell_dim[1];

    const float r2 = calc_r2(radius);

    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 v_r2 = md_mm256_set1_ps(r2);

    const ivec4_t cdim_v  = ivec4_set(cdim[0], cdim[1], cdim[2], 0);
    const ivec4_t cdim_1v = ivec4_sub(cdim_v, ivec4_set(1, 1, 1, 0));
    const ivec4_t zero_v  = ivec4_set1(0);

    const ivec4_t pmask_v = ivec4_set(
        (acc->flags & MD_UNITCELL_PBC_X) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Y) ? 0xFFFFFFFF : 0,
        (acc->flags & MD_UNITCELL_PBC_Z) ? 0xFFFFFFFF : 0,
        0);

    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const vec4_t pbc_mask = vec4_set((acc->flags & MD_UNITCELL_PBC_X) ? val : 0, (acc->flags & MD_UNITCELL_PBC_Y) ? val : 0, (acc->flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);
    
    const vec4_t fcell_dim = vec4_set((float)cdim[0], (float)cdim[1], (float)cdim[2], 0);
    const vec4_t r4 = vec4_set((float)center[0], (float)center[1], (float)center[2], 0);
    vec4_t f4 = vec4_cart_to_fract(r4, acc);
    f4 = vec4_blend(vec4_fract(f4), f4, pbc_mask);
    const ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f4, fcell_dim)));

    const md_256 v_xi = md_mm256_set1_ps(f4.x);
    const md_256 v_yi = md_mm256_set1_ps(f4.y);
    const md_256 v_zi = md_mm256_set1_ps(f4.z);

    // Iterate neighbor cells
    for (size_t n_idx = 0; n_idx < num_neighbors; ++n_idx) {
        const ivec4_t nbr = ivec4_load(neighbors[n_idx]);
        ivec4_t n_v = ivec4_add(c_v, nbr);
        
        const ivec4_t wrap_upper = ivec4_cmpgt(n_v, cdim_1v);
        const ivec4_t wrap_lower = ivec4_cmplt(n_v, zero_v);
        const ivec4_t wrap_any   = ivec4_or(wrap_upper, wrap_lower);
        
        // Skip nonperiodic wraps
        if (ivec4_any(ivec4_andnot(wrap_any, pmask_v))) continue;
        
        // Apply wrapping
        n_v = ivec4_add(n_v, ivec4_and(wrap_lower, cdim_v));
        n_v = ivec4_sub(n_v, ivec4_and(wrap_upper, cdim_v));
        
        int c[4];
        ivec4_store(c, n_v);

        const uint32_t ci  = CELL_INDEX((uint32_t)c[0], (uint32_t)c[1], (uint32_t)c[2]);
        const uint32_t off = CELL_OFFSET(ci);
        const uint32_t len = CELL_LENGTH(ci);
        if (len == 0) continue;

        const float* elem_x = acc->elem_x + off;
        const float* elem_y = acc->elem_y + off;
        const float* elem_z = acc->elem_z + off;
        const uint32_t* elem_idx = acc->elem_idx + off;

        // Compute shift vector (+1 for lower wrap, -1 for upper)
        const ivec4_t shift_i = ivec4_sub(
            ivec4_and(wrap_upper, ivec4_set1(1)),
            ivec4_and(wrap_lower, ivec4_set1(1))
        );

        const vec4_t f_shift = vec4_from_ivec4(shift_i);

        const md_256 shift_x = md_mm256_set1_ps(f_shift.x);
        const md_256 shift_y = md_mm256_set1_ps(f_shift.y);
        const md_256 shift_z = md_mm256_set1_ps(f_shift.z);

        POSSIBLY_INVOKE_CALLBACK_POINT_ORT(ALIGN_TO(len, 8));

        const md_256i v_len = md_mm256_set1_epi32(len);
        md_256i v_j = md_mm256_set_epi32(7,6,5,4,3,2,1,0);
        for (uint32_t j = 0; j < len; j += 8) {
            const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len));

            md_256 v_xj = md_mm256_loadu_ps(elem_x + j);
            md_256 v_yj = md_mm256_loadu_ps(elem_y + j);
            md_256 v_zj = md_mm256_loadu_ps(elem_z + j);

            v_xj = md_mm256_add_ps(v_xj, shift_x);
            v_yj = md_mm256_add_ps(v_yj, shift_y);
            v_zj = md_mm256_add_ps(v_zj, shift_z);

            const md_256 v_dx = md_mm256_sub_ps(v_xi, v_xj);
            const md_256 v_dy = md_mm256_sub_ps(v_yi, v_yj);
            const md_256 v_dz = md_mm256_sub_ps(v_zi, v_zj);

            md_256 v_d2 = distance_squared_ort_256(v_dx, v_dy, v_dz, G00, G11, G22);

            const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

            const int mask = md_mm256_movemask_ps(v_mask);
            if (mask) {
                const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_idx + j));

                v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                v_xj   = md_mm256_permutevar8x32_ps(v_xj, v_idx_mask);
                v_yj   = md_mm256_permutevar8x32_ps(v_yj, v_idx_mask);
                v_zj   = md_mm256_permutevar8x32_ps(v_zj, v_idx_mask);

                ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                md_mm256_storeu_epi32(buf_i + count, v_idxj);
                md_mm256_storeu_ps   (buf_x + count, v_xj);
                md_mm256_storeu_ps   (buf_y + count, v_yj);
                md_mm256_storeu_ps   (buf_z + count, v_zj);

                count += popcnt32(mask);
            }

            v_j = md_mm256_add_epi32(v_j, md_mm256_set1_epi32(8));
        }
    }

    FLUSH_TAIL_POINT_ORT();
}

static void for_each_point_in_sphere_triclinic(const md_spatial_acc_t* acc, const double center[3], double radius, md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc->flags & MD_UNITCELL_TRICLINIC && acc->flags & MD_UNITCELL_PBC_ALL);

    const int ncell[3] = {
        (int)ceil(radius * (double)acc->inv_cell_ext[0] * acc->cell_dim[0]),
        (int)ceil(radius * (double)acc->inv_cell_ext[1] * acc->cell_dim[1]),
        (int)ceil(radius * (double)acc->inv_cell_ext[2] * acc->cell_dim[2]),
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_point_in_sphere_triclinic: radius too large for cell size");
        return;
    }

    // Build neighbor offsets once
    int neighbors[SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS * SPATIAL_ACC_MAX_NEIGHBOR_CELLS][4];
    size_t num_neighbors = generate_neighbors4(neighbors, ncell);

    // Allocate intermediate buffers for passing to callback
    uint32_t buf_i[SPATIAL_ACC_BUFLEN];
    float    buf_x[SPATIAL_ACC_BUFLEN];
    float    buf_y[SPATIAL_ACC_BUFLEN];
    float    buf_z[SPATIAL_ACC_BUFLEN];

    size_t   count = 0;

    // Setup constants
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t c0  = acc->cell_dim[0];
    const uint32_t c01 = acc->cell_dim[0] * acc->cell_dim[1];

    const float r2 = calc_r2(radius);

    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);
    const md_256 v_r2 = md_mm256_set1_ps(r2);

    const md_256i add8 = md_mm256_set1_epi32(8);
    const md_256i inc  = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    // Fractional center and center cell
 const vec4_t r4 = vec4_set((float)center[0], (float)center[1], (float)center[2], 0);
	const vec4_t d4 = vec4_set((float)cdim[0], (float)cdim[1], (float)cdim[2], 0);
	const vec4_t f4 = vec4_fract(vec4_cart_to_fract(r4, acc));
    const ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f4, d4)));

    const md_256 v_xi = md_mm256_set1_ps(f4.x);
    const md_256 v_yi = md_mm256_set1_ps(f4.y);
    const md_256 v_zi = md_mm256_set1_ps(f4.z);

    for (uint32_t n_idx = 0; n_idx < num_neighbors; ++n_idx) {
        const ivec4_t nbr = ivec4_load(neighbors[n_idx]);
        ivec4_t cell_v = ivec4_add(c_v, nbr);

        const ivec4_t wrap_upper = ivec4_cmpgt(cell_v, ivec4_sub(ivec4_set((int)cdim[0], (int)cdim[1], (int)cdim[2], 0), ivec4_set1(1)));
        const ivec4_t wrap_lower = ivec4_cmplt(cell_v, ivec4_set1(0));

        // Apply wrapping (triclinic -> periodic in all dims)
        cell_v = ivec4_add(cell_v, ivec4_and(wrap_lower, ivec4_set((int)cdim[0], (int)cdim[1], (int)cdim[2], 0)));
        cell_v = ivec4_sub(cell_v, ivec4_and(wrap_upper, ivec4_set((int)cdim[0], (int)cdim[1], (int)cdim[2], 0)));

        int c[4];
        ivec4_store(c, cell_v);

        const uint32_t ci  = CELL_INDEX((uint32_t)c[0], (uint32_t)c[1], (uint32_t)c[2]);
        const uint32_t off = CELL_OFFSET(ci);
        const uint32_t len = CELL_LENGTH(ci);
        if (len == 0) continue;

        const md_256i v_len = md_mm256_set1_epi32(len);

        const float* elem_x = acc->elem_x + off;
        const float* elem_y = acc->elem_y + off;
        const float* elem_z = acc->elem_z + off;
        const uint32_t* elem_idx = acc->elem_idx + off;

        // shifting due to wrap: compute shift_i (+1 / -1) like other triclinic functions
        ivec4_t shift_i = ivec4_sub(
            ivec4_and(wrap_upper, ivec4_set1(1)),
            ivec4_and(wrap_lower, ivec4_set1(1))
        );
        vec4_t  shift_f = vec4_from_ivec4(shift_i);

		const md_256 shift_x = md_mm256_set1_ps(shift_f.x);
		const md_256 shift_y = md_mm256_set1_ps(shift_f.y);
		const md_256 shift_z = md_mm256_set1_ps(shift_f.z);

        POSSIBLY_INVOKE_CALLBACK_POINT_FRACT_TRI(ALIGN_TO(len, 8));

        md_256i v_j = inc;
        for (uint32_t j = 0; j < len; j += 8) {
            const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len));

            md_256 v_xj = md_mm256_loadu_ps(elem_x + j);
            md_256 v_yj = md_mm256_loadu_ps(elem_y + j);
            md_256 v_zj = md_mm256_loadu_ps(elem_z + j);

            v_xj = md_mm256_add_ps(v_xj, shift_x);
            v_yj = md_mm256_add_ps(v_yj, shift_y);
            v_zj = md_mm256_add_ps(v_zj, shift_z);

            const md_256 v_dx = md_mm256_sub_ps(v_xi, v_xj);
            const md_256 v_dy = md_mm256_sub_ps(v_yi, v_yj);
            const md_256 v_dz = md_mm256_sub_ps(v_zi, v_zj);

            md_256 v_d2 = distance_squared_tri_256(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

            const md_256 v_mask = md_mm256_and_ps(md_mm256_cmple_ps(v_d2, v_r2), j_mask);

            // Fill buffers with results
            int mask = md_mm256_movemask_ps(v_mask);
            if (mask) {
                const md_256i v_idx_mask = md_mm256_compression_mask_8x32(mask);
                md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_idx + j));

                v_idxj = md_mm256_permutevar8x32_epi32(v_idxj, v_idx_mask);
                v_xj   = md_mm256_permutevar8x32_ps(v_xj, v_idx_mask);
                v_yj   = md_mm256_permutevar8x32_ps(v_yj, v_idx_mask);
                v_zj   = md_mm256_permutevar8x32_ps(v_zj, v_idx_mask);

                ASSERT(count + 8 <= SPATIAL_ACC_BUFLEN);

                md_mm256_storeu_epi32(buf_i + count, v_idxj);
                md_mm256_storeu_ps   (buf_x + count, v_xj);
                md_mm256_storeu_ps   (buf_y + count, v_yj);
                md_mm256_storeu_ps   (buf_z + count, v_zj);

                count += popcnt32(mask);
            }

            v_j = md_mm256_add_epi32(v_j, add8);
        }
    }

    FLUSH_TAIL_POINT_CART_TRI();
}

#undef SPATIAL_ACC_BUFLEN
#undef CELL_INDEX
#undef CELL_OFFSET
#undef CELL_LENGTH

void md_spatial_acc_for_each_internal_pair_in_neighboring_cells(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param) {
    ASSERT(acc);
	ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_internal_pair_in_neighboring_cells_triclinic(acc, callback, user_param);
    } else {
		for_each_internal_pair_in_neighboring_cells_ortho(acc, callback, user_param);
	}
}

void md_spatial_acc_for_each_internal_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_internal_pair_within_cutoff_triclinic(acc, cutoff, callback, user_param);
    } else {
        for_each_internal_pair_within_cutoff_ortho(acc, cutoff, callback, user_param);
	}
}

// Iterate over external points against points within the spatial acceleration structure for a supplied cutoff
// The external points are not part of the spatial acceleration structure and will be represented in the callback as the 'i' indices and the internal points are the 'j' indices
void md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(const md_spatial_acc_t* acc, const md_coord_stream_t* ext_stream, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags) {
    ASSERT(acc);
    ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_external_pair_within_cutoff_triclinic(acc, ext_stream, cutoff, callback, user_param, flags);
    } else {
        for_each_external_pair_within_cutoff_ortho(acc, ext_stream, cutoff, callback, user_param, flags);
    }
}

void md_spatial_acc_for_each_point_in_aabb(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);

    for (int i = 0; i < 3; ++i) {
        if (aabb_rad[i] < 0) {
            MD_LOG_ERROR("md_spatial_acc_for_each_point_in_aabb: negative radius not allowed");
            return;
        }
    }
	
    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_point_in_aabb_triclinic(acc, aabb_cen, aabb_rad, callback, user_param);
    }
    else {
        for_each_point_in_aabb_ortho(acc, aabb_cen, aabb_rad, callback, user_param);
    }
}

void md_spatial_acc_for_each_point_in_sphere(const md_spatial_acc_t* acc, const double center[3], double radius, md_spatial_acc_point_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);
    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_point_in_sphere_triclinic(acc, center, radius, callback, user_param);
    }
    else {
        for_each_point_in_sphere_ortho(acc, center, radius, callback, user_param);
    }
}

# if 0
// Iterate over external points against points within the spatial acceleration structure in neighboring cells (1-cell neighborhood)
bool md_spatial_acc_for_each_external_point_in_neighboring_cells(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, const int32_t* ext_idx, size_t ext_count, md_spatial_acc_pair_callback_t callback, void* user_param) {
    // Fallback: use cutoff equal to smallest cell extent to approximate 1-cell neighborhood
    const double min_ext = MIN(acc->cell_ext[0], MIN(acc->cell_ext[1], acc->cell_ext[2]));
    if (min_ext == 0.0) return false;
    return md_spatial_acc_for_each_external_point_within_cutoff(acc, ext_x, ext_y, ext_z, ext_idx, ext_count, min_ext, callback, user_param);
}
#endif
