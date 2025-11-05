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

static const uint64_t avx_compress_lut[256] = {
    0x0000000000000000ULL,   0x0000000000000000ULL,   0x0000000000000001ULL,   0x0000000000000100ULL,
    0x0000000000000002ULL,   0x0000000000000200ULL,   0x0000000000000201ULL,   0x0000000000020100ULL,
    0x0000000000000003ULL,   0x0000000000000300ULL,   0x0000000000000301ULL,   0x0000000000030100ULL,
    0x0000000000000302ULL,   0x0000000000030200ULL,   0x0000000000030201ULL,   0x0000000003020100ULL,
    0x0000000000000004ULL,   0x0000000000000400ULL,   0x0000000000000401ULL,   0x0000000000040100ULL,
    0x0000000000000402ULL,   0x0000000000040200ULL,   0x0000000000040201ULL,   0x0000000004020100ULL,
    0x0000000000000403ULL,   0x0000000000040300ULL,   0x0000000000040301ULL,   0x0000000004030100ULL,
    0x0000000000040302ULL,   0x0000000004030200ULL,   0x0000000004030201ULL,   0x0000000403020100ULL,
    0x0000000000000005ULL,   0x0000000000000500ULL,   0x0000000000000501ULL,   0x0000000000050100ULL,
    0x0000000000000502ULL,   0x0000000000050200ULL,   0x0000000000050201ULL,   0x0000000005020100ULL,
    0x0000000000000503ULL,   0x0000000000050300ULL,   0x0000000000050301ULL,   0x0000000005030100ULL,
    0x0000000000050302ULL,   0x0000000005030200ULL,   0x0000000005030201ULL,   0x0000000503020100ULL,
    0x0000000000000504ULL,   0x0000000000050400ULL,   0x0000000000050401ULL,   0x0000000005040100ULL,
    0x0000000000050402ULL,   0x0000000005040200ULL,   0x0000000005040201ULL,   0x0000000504020100ULL,
    0x0000000000050403ULL,   0x0000000005040300ULL,   0x0000000005040301ULL,   0x0000000504030100ULL,
    0x0000000005040302ULL,   0x0000000504030200ULL,   0x0000000504030201ULL,   0x0000050403020100ULL,
    0x0000000000000006ULL,   0x0000000000000600ULL,   0x0000000000000601ULL,   0x0000000000060100ULL,
    0x0000000000000602ULL,   0x0000000000060200ULL,   0x0000000000060201ULL,   0x0000000006020100ULL,
    0x0000000000000603ULL,   0x0000000000060300ULL,   0x0000000000060301ULL,   0x0000000006030100ULL,
    0x0000000000060302ULL,   0x0000000006030200ULL,   0x0000000006030201ULL,   0x0000000603020100ULL,
    0x0000000000000604ULL,   0x0000000000060400ULL,   0x0000000000060401ULL,   0x0000000006040100ULL,
    0x0000000000060402ULL,   0x0000000006040200ULL,   0x0000000006040201ULL,   0x0000000604020100ULL,
    0x0000000000060403ULL,   0x0000000006040300ULL,   0x0000000006040301ULL,   0x0000000604030100ULL,
    0x0000000006040302ULL,   0x0000000604030200ULL,   0x0000000604030201ULL,   0x0000060403020100ULL,
    0x0000000000000605ULL,   0x0000000000060500ULL,   0x0000000000060501ULL,   0x0000000006050100ULL,
    0x0000000000060502ULL,   0x0000000006050200ULL,   0x0000000006050201ULL,   0x0000000605020100ULL,
    0x0000000000060503ULL,   0x0000000006050300ULL,   0x0000000006050301ULL,   0x0000000605030100ULL,
    0x0000000006050302ULL,   0x0000000605030200ULL,   0x0000000605030201ULL,   0x0000060503020100ULL,
    0x0000000000060504ULL,   0x0000000006050400ULL,   0x0000000006050401ULL,   0x0000000605040100ULL,
    0x0000000006050402ULL,   0x0000000605040200ULL,   0x0000000605040201ULL,   0x0000060504020100ULL,
    0x0000000006050403ULL,   0x0000000605040300ULL,   0x0000000605040301ULL,   0x0000060504030100ULL,
    0x0000000605040302ULL,   0x0000060504030200ULL,   0x0000060504030201ULL,   0x0006050403020100ULL,
    0x0000000000000007ULL,   0x0000000000000700ULL,   0x0000000000000701ULL,   0x0000000000070100ULL,
    0x0000000000000702ULL,   0x0000000000070200ULL,   0x0000000000070201ULL,   0x0000000007020100ULL,
    0x0000000000000703ULL,   0x0000000000070300ULL,   0x0000000000070301ULL,   0x0000000007030100ULL,
    0x0000000000070302ULL,   0x0000000007030200ULL,   0x0000000007030201ULL,   0x0000000703020100ULL,
    0x0000000000000704ULL,   0x0000000000070400ULL,   0x0000000000070401ULL,   0x0000000007040100ULL,
    0x0000000000070402ULL,   0x0000000007040200ULL,   0x0000000007040201ULL,   0x0000000704020100ULL,
    0x0000000000070403ULL,   0x0000000007040300ULL,   0x0000000007040301ULL,   0x0000000704030100ULL,
    0x0000000007040302ULL,   0x0000000704030200ULL,   0x0000000704030201ULL,   0x0000070403020100ULL,
    0x0000000000000705ULL,   0x0000000000070500ULL,   0x0000000000070501ULL,   0x0000000007050100ULL,
    0x0000000000070502ULL,   0x0000000007050200ULL,   0x0000000007050201ULL,   0x0000000705020100ULL,
    0x0000000000070503ULL,   0x0000000007050300ULL,   0x0000000007050301ULL,   0x0000000705030100ULL,
    0x0000000007050302ULL,   0x0000000705030200ULL,   0x0000000705030201ULL,   0x0000070503020100ULL,
    0x0000000000070504ULL,   0x0000000007050400ULL,   0x0000000007050401ULL,   0x0000000705040100ULL,
    0x0000000007050402ULL,   0x0000000705040200ULL,   0x0000000705040201ULL,   0x0000070504020100ULL,
    0x0000000007050403ULL,   0x0000000705040300ULL,   0x0000000705040301ULL,   0x0000070504030100ULL,
    0x0000000705040302ULL,   0x0000070504030200ULL,   0x0000070504030201ULL,   0x0007050403020100ULL,
    0x0000000000000706ULL,   0x0000000000070600ULL,   0x0000000000070601ULL,   0x0000000007060100ULL,
    0x0000000000070602ULL,   0x0000000007060200ULL,   0x0000000007060201ULL,   0x0000000706020100ULL,
    0x0000000000070603ULL,   0x0000000007060300ULL,   0x0000000007060301ULL,   0x0000000706030100ULL,
    0x0000000007060302ULL,   0x0000000706030200ULL,   0x0000000706030201ULL,   0x0000070603020100ULL,
    0x0000000000070604ULL,   0x0000000007060400ULL,   0x0000000007060401ULL,   0x0000000706040100ULL,
    0x0000000007060402ULL,   0x0000000706040200ULL,   0x0000000706040201ULL,   0x0000070604020100ULL,
    0x0000000007060403ULL,   0x0000000706040300ULL,   0x0000000706040301ULL,   0x0000070604030100ULL,
    0x0000000706040302ULL,   0x0000070604030200ULL,   0x0000070604030201ULL,   0x0007060403020100ULL,
    0x0000000000070605ULL,   0x0000000007060500ULL,   0x0000000007060501ULL,   0x0000000706050100ULL,
    0x0000000007060502ULL,   0x0000000706050200ULL,   0x0000000706050201ULL,   0x0000070605020100ULL,
    0x0000000007060503ULL,   0x0000000706050300ULL,   0x0000000706050301ULL,   0x0000070605030100ULL,
    0x0000000706050302ULL,   0x0000070605030200ULL,   0x0000070605030201ULL,   0x0007060503020100ULL,
    0x0000000007060504ULL,   0x0000000706050400ULL,   0x0000000706050401ULL,   0x0000070605040100ULL,
    0x0000000706050402ULL,   0x0000070605040200ULL,   0x0000070605040201ULL,   0x0007060504020100ULL,
    0x0000000706050403ULL,   0x0000070605040300ULL,   0x0000070605040301ULL,   0x0007060504030100ULL,
    0x0000070605040302ULL,   0x0007060504030200ULL,   0x0007060504030201ULL,   0x0706050403020100ULL,
};


typedef struct {
    float x, y, z;
    uint32_t idx;
} elem_t;

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

    MEMSET(acc->cell_dim, 0, sizeof(acc->cell_dim));
    MEMSET(acc->cell_min, 0, sizeof(acc->cell_min));
    MEMSET(acc->cell_max, 0, sizeof(acc->cell_max));

	acc->flags = 0;
    acc->origin[0] = acc->origin[1] = acc->origin[2] = 0.0f;
}

void md_spatial_acc_init(md_spatial_acc_t* acc, const float* in_x, const float* in_y, const float* in_z, const int32_t* in_idx, size_t count, double in_cell_ext, const md_unitcell_t* unitcell) {
    ASSERT(acc);
    if (!acc->alloc) {
        MD_LOG_ERROR("Must have allocator set within spatial acc");
        return;
    }

	// Reset acc, but to not free memory
    md_spatial_acc_reset(acc);

    double A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    uint32_t flags = 0;

    if (unitcell) {
        md_unitcell_basis_extract(A, unitcell);
        md_unitcell_inv_basis_extract(I, unitcell);
        flags = md_unitcell_flags(unitcell);
    }

    vec4_t coord_offset = vec4_zero();

    if ((flags & MD_UNITCELL_PBC_ALL) != MD_UNITCELL_PBC_ALL) {
        ASSERT((flags & MD_UNITCELL_TRICLINIC) == 0);
        // Unit cell either missing or not periodic along one or more axis
        vec4_t aabb_min = {0}, aabb_max = {0};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, in_x, in_y, in_z, NULL, NULL, count);

        // Construct A and I from aabb extent
        vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);

        if ((flags & MD_UNITCELL_PBC_X) == 0) {
            coord_offset.x = -aabb_min.x;
            if (aabb_ext.x > 0.0f) {
                A[0][0] = aabb_ext.x;
                I[0][0] = 1.0 / aabb_ext.x;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Y) == 0) {
            coord_offset.y = -aabb_min.y;
            if (aabb_ext.y > 0.0f) {
                A[1][1] = aabb_ext.y;
                I[1][1] = 1.0 / aabb_ext.y;
            }
        }
        if ((flags & MD_UNITCELL_PBC_Z) == 0) {
            coord_offset.z = -aabb_min.z;
            if (aabb_ext.z > 0.0f) {
                A[2][2] = aabb_ext.z;
                I[2][2] = 1.0 / aabb_ext.z;
            }
        }
    }

    // Precompute metric G = A^T A
    double G00 =     (A[0][0] * A[0][0] + A[1][0] * A[1][0] + A[2][0] * A[2][0]); // dot(a,a)
    double G11 =     (A[0][1] * A[0][1] + A[1][1] * A[1][1] + A[2][1] * A[2][1]); // dot(b,b)
    double G22 =     (A[0][2] * A[0][2] + A[1][2] * A[1][2] + A[2][2] * A[2][2]); // dot(c,c)
    double H01 = 2 * (A[0][0] * A[0][1] + A[1][0] * A[1][1] + A[2][0] * A[2][1]); // 2 dot(a,b)
    double H02 = 2 * (A[0][0] * A[0][2] + A[1][0] * A[1][2] + A[2][0] * A[2][2]); // 2 dot(a, c)
    double H12 = 2 * (A[0][1] * A[0][2] + A[1][1] * A[1][2] + A[2][1] * A[2][2]); // 2 dot(b, c)

    // Choose grid resolution. Heuristic:  cell_dim ≈ |A|/CELL_EXT.
    // Using ~cutoff-ish spacing gives good pruning.
    const double CELL_EXT = in_cell_ext > 0.0 ? in_cell_ext : 6.0;

    // Estimate cell_dim by measuring the extents of the box vectors (norms of columns of A)
    // This is only a heuristic for bin counts; the grid is still in fractional space.
    double ax = sqrt(A[0][0] * A[0][0] + A[1][0] * A[1][0] + A[2][0] * A[2][0]);
    double by = sqrt(A[0][1] * A[0][1] + A[1][1] * A[1][1] + A[2][1] * A[2][1]);
    double cz = sqrt(A[0][2] * A[0][2] + A[1][2] * A[1][2] + A[2][2] * A[2][2]);

    uint32_t cell_dim[3] = {
        CLAMP((uint32_t)(ax / CELL_EXT), 1, 1024),
        CLAMP((uint32_t)(by / CELL_EXT), 1, 1024),
        CLAMP((uint32_t)(cz / CELL_EXT), 1, 1024),
    };

	// Calculate the effective cell extents in real space based on cell_dim
    double cell_ext[3] = {
        ax / (double)cell_dim[0],
        by / (double)cell_dim[1],
        cz / (double)cell_dim[2],
	};

    const uint32_t c0  = cell_dim[0];
    const uint32_t c01 = cell_dim[0] * cell_dim[1];
    const size_t num_cells = (size_t)cell_dim[0] * cell_dim[1] * cell_dim[2];

    MD_LOG_DEBUG("cell_dim: %i %i %i", cell_dim[0], cell_dim[1], cell_dim[2]);

    // Temporary arrays
    const size_t temp_arena_page_size = MEGABYTES(4);
    md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), temp_arena_page_size);
    uint32_t* local_idx = (uint32_t*)md_arena_allocator_push(temp_arena, count * sizeof(uint32_t));
    uint32_t* cell_idx  = (uint32_t*)md_arena_allocator_push(temp_arena, count * sizeof(uint32_t));
    elem_t* scratch_s   = (elem_t*)  md_arena_allocator_push(temp_arena, count * sizeof(elem_t));  // unsorted fractional coords

    // Resize / allocate persistent arrays
    size_t alloc_len = ALIGN_TO(count, 16);

    md_array_resize(acc->elem_x, alloc_len, acc->alloc);
    md_array_resize(acc->elem_y, alloc_len, acc->alloc);
    md_array_resize(acc->elem_z, alloc_len, acc->alloc);
    md_array_resize(acc->elem_idx, alloc_len, acc->alloc);
    md_array_resize(acc->cell_off, num_cells + 1, acc->alloc);
    MEMSET(acc->cell_off, 0, (num_cells + 1) * sizeof(uint32_t));

    const vec4_t  fcell_dim = vec4_set((float)cell_dim[0], (float)cell_dim[1], (float)cell_dim[2], 0);
    const ivec4_t icell_min = ivec4_set1(0);
    const ivec4_t icell_max = ivec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);

    ivec4_t cell_min = ivec4_set1(INT32_MAX);
    ivec4_t cell_max = ivec4_set1(INT32_MIN);

    mat4_t I4 = {
        (float)I[0][0], (float)I[0][1], (float)I[0][2], 0,
        (float)I[1][0], (float)I[1][1], (float)I[1][2], 0,
        (float)I[2][0], (float)I[2][1], (float)I[2][2], 0,
        0, 0, 0, 0,
    };

    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const vec4_t pbc_mask = vec4_set((flags & MD_UNITCELL_PBC_X) ? val : 0, (flags & MD_UNITCELL_PBC_Y) ? val : 0, (flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);

    // 1) Convert to fractional, wrap periodic axes into [0,1), bin to cells
    for (size_t i = 0; i < count; ++i) {
		uint32_t idx = (in_idx) ? (uint32_t)in_idx[i] : (uint32_t)i;

        vec4_t r = {in_x[idx], in_y[idx], in_z[idx], 0};
        r = vec4_add(r, coord_offset);

        // Fractional coordinates
        vec4_t c = mat4_mul_vec4(I4, r);
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

        local_idx[i] = acc->cell_off[ci]++;  // count for now
        cell_idx[i]  = (uint32_t)ci;

        // stash fractional coordinates
        scratch_s[i] = (elem_t){c.x, c.y, c.z, idx};
    }

    // 2) Prefix sum cell offsets
    uint32_t sum = 0;
    for (size_t ci = 0; ci <= num_cells; ++ci) {
        uint32_t len = acc->cell_off[ci];
        acc->cell_off[ci] = sum;
        sum += len;
    }
    ASSERT(sum == count);

    // 3) Scatter fractional coords into 'elements' in cell order
    for (size_t i = 0; i < count; ++i) {
        uint32_t dst = acc->cell_off[cell_idx[i]] + local_idx[i];
        ASSERT(dst < count);
        acc->elem_x[dst] = scratch_s[i].x;
        acc->elem_y[dst] = scratch_s[i].y;
        acc->elem_z[dst] = scratch_s[i].z;
        acc->elem_idx[dst] = scratch_s[i].idx;
    }

    acc->num_elems = count;
    int tmp_min[4], tmp_max[4];
    ivec4_store(tmp_min, cell_min);
    ivec4_store(tmp_max, cell_max);

    MEMCPY(acc->cell_min, tmp_min,  sizeof(acc->cell_min));
    MEMCPY(acc->cell_max, tmp_max,  sizeof(acc->cell_max));
    MEMCPY(acc->cell_dim, cell_dim, sizeof(acc->cell_dim));
	acc->cell_ext[0] = (float)cell_ext[0];
	acc->cell_ext[1] = (float)cell_ext[1];
	acc->cell_ext[2] = (float)cell_ext[2];
    acc->num_cells = num_cells;
        
    acc->G00 = (float)G00;
    acc->G11 = (float)G11;
    acc->G22 = (float)G22;
    acc->H01 = (float)H01;
    acc->H02 = (float)H02;
    acc->H12 = (float)H12;

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
    acc->origin[0] = coord_offset.x;
    acc->origin[1] = coord_offset.y;
    acc->origin[2] = coord_offset.z;

    acc->flags = flags;

    md_arena_allocator_destroy(temp_arena);
}

// The 13 'forward' neighbor offsets in 3D
static const int FWD_NBRS[13][4] = {
    { 1, 0, 0, 0},  {-1, 1, 0, 0}, {0, 1, 0, 0}, { 1, 1, 0, 0}, {-1, -1, 1, 0}, {0, -1, 1, 0}, {1, -1, 1, 0},
    {-1, 0, 1, 0},  { 0, 0, 1, 0}, {1, 0, 1, 0}, {-1, 1, 1, 0}, { 0,  1, 1, 0}, {1, 1, 1, 0},
};

// Generate forward neighbor offsets for a 3D grid cell
// - out: user-provided array of size at least (ncell*2+1)^3 - 1
// - ncell: number of cells in the neighborhood (1 → 1-cell, 2 → 2-cell, etc.)
// - Each element in out is int[3] representing offset {dx, dy, dz}
// Returns the number of neighbors written
static inline size_t generate_forward_neighbors3(int out[][3], const int ncell[3]) {
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
                    count++;
                }
            }
        }
    }
    return count;
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

static inline md_128 distance_squared_tric_sse(md_128 dx, md_128 dy, md_128 dz, md_128 G00, md_128 G11, md_128 G22, md_128 H01, md_128 H02, md_128 H12) {
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

static inline md_256 distance_squared_tric_avx(md_256 dx, md_256 dy, md_256 dz, md_256 G00, md_256 G11, md_256 G22, md_256 H01, md_256 H02, md_256 H12) {
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

#define POSSIBLY_INVOKE_CALLBACK(estimated_count) \
    if (count + (estimated_count) >= SPATIAL_ACC_BUFLEN) { \
        callback(buf_i, buf_j, buf_d2, count, user_param); \
        count = 0; \
    } \

#define FLUSH_TAIL() \
    if (count) { \
        callback(buf_i, buf_j, buf_d2, count, user_param); \
        count = 0; \
    } \

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

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t cmin[3] = { acc->cell_min[0], acc->cell_min[1], acc->cell_min[2] };
    const uint32_t cmax[3] = { acc->cell_max[0], acc->cell_max[1], acc->cell_max[2] };
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

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_i - (i + 1), 8));

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
                        const md_256 v_d2 = distance_squared_tric_avx(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

						const md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_i_idx + j));
						int mask = md_mm256_movemask_ps(v_mask);

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
                for (uint32_t n = 0; n < ARRAY_SIZE(FWD_NBRS); ++n) {
                    const ivec4_t fwd_v = ivec4_load(FWD_NBRS[n]);
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

                    const float shift_xf = (float)simde_mm_extract_epi32(shift_i, 0);
                    const float shift_yf = (float)simde_mm_extract_epi32(shift_i, 1);
                    const float shift_zf = (float)simde_mm_extract_epi32(shift_i, 2);

                    const md_256 shift_x = md_mm256_set1_ps(shift_xf);
                    const md_256 shift_y = md_mm256_set1_ps(shift_yf);
                    const md_256 shift_z = md_mm256_set1_ps(shift_zf);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_j, 8));

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
                            const md_256 v_d2 = distance_squared_tric_avx(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                            int mask = md_mm256_movemask_ps(v_mask);

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

    FLUSH_TAIL();
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

    size_t count = 0;

    const float*    element_x = acc->elem_x;
    const float*    element_y = acc->elem_y;
    const float*    element_z = acc->elem_z;
    const uint32_t* element_i = acc->elem_idx;
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t cmin[3] = { acc->cell_min[0], acc->cell_min[1], acc->cell_min[2] };
    const uint32_t cmax[3] = { acc->cell_max[0], acc->cell_max[1], acc->cell_max[2] };
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

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_i - (i + 1), 8));

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
                        const md_256 v_d2 = distance_squared_ortho_avx(v_dx, v_dy, v_dz, G00, G11, G22);

						const md_256i v_idxj = md_mm256_loadu_epi32(elem_i_idx + j);
						int mask = md_mm256_movemask_ps(v_mask);

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
                for (uint32_t n = 0; n < ARRAY_SIZE(FWD_NBRS); ++n) {
                    const ivec4_t fwd_v = ivec4_load(FWD_NBRS[n]);
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

                    const float shift_xf = (float)simde_mm_extract_epi32(shift_i, 0);
                    const float shift_yf = (float)simde_mm_extract_epi32(shift_i, 1);
                    const float shift_zf = (float)simde_mm_extract_epi32(shift_i, 2);

                    const md_256 shift_x = md_mm256_set1_ps(shift_xf);
                    const md_256 shift_y = md_mm256_set1_ps(shift_yf);
                    const md_256 shift_z = md_mm256_set1_ps(shift_zf);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_j, 8));

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
                            const md_256 v_d2 = distance_squared_ortho_avx(v_dx, v_dy, v_dz, G00, G11, G22);

                            int mask = md_mm256_movemask_ps(v_mask);

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

    FLUSH_TAIL();
}

static bool for_each_internal_pair_within_cutoff_triclinic(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    int ncell[3] = {
        acc->cell_ext[0] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[0]) : 0,
        acc->cell_ext[1] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[1]) : 0,
        acc->cell_ext[2] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[2]) : 0,
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_pair_within_cutoff_ortho: cutoff too large for cell size");
        return false;
    }

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);
    const md_256 v_r2 = md_mm256_set1_ps((float)(cutoff * cutoff));

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
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t cmin[3] = { acc->cell_min[0], acc->cell_min[1], acc->cell_min[2] };
    const uint32_t cmax[3] = { acc->cell_max[0], acc->cell_max[1], acc->cell_max[2] };
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

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_i - (i + 1), 8));

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
                        md_256 v_d2 = distance_squared_tric_avx(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                        const md_256 v_mask = md_mm256_and_ps(md_mm256_cmplt_ps(v_d2, v_r2), j_mask);

                        // Fill buffers with results
                        int mask = md_mm256_movemask_ps(v_mask);

                        if (mask) {
						    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_i_idx + j));

                            uint64_t key = avx_compress_lut[mask];
                            md_256i v_perm_mask = simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));

						    v_idxj = simde_mm256_permutevar8x32_epi32(v_idxj, v_perm_mask);
						    v_d2   = simde_mm256_permutevar8x32_ps(v_d2,      v_perm_mask);

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

                    const float shift_xf = (float)simde_mm_extract_epi32(shift_i, 0);
                    const float shift_yf = (float)simde_mm_extract_epi32(shift_i, 1);
                    const float shift_zf = (float)simde_mm_extract_epi32(shift_i, 2);

                    const md_256 shift_x = md_mm256_set1_ps(shift_xf);
                    const md_256 shift_y = md_mm256_set1_ps(shift_yf);
                    const md_256 shift_z = md_mm256_set1_ps(shift_zf);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_j, 8));

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
                            md_256 v_d2 = distance_squared_tric_avx(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                            const md_256 v_mask = md_mm256_and_ps(md_mm256_cmplt_ps(v_d2, v_r2), j_mask);

                            // Fill buffers with results
                            int mask = md_mm256_movemask_ps(v_mask);

                            if (mask) {
                                md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                                uint64_t key = avx_compress_lut[mask];
                                md_256i v_perm_mask = simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));

                                v_idxj = simde_mm256_permutevar8x32_epi32(v_idxj, v_perm_mask);
                                v_d2   = simde_mm256_permutevar8x32_ps(v_d2,      v_perm_mask);

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

    FLUSH_TAIL();

    return true;
}

static bool for_each_internal_pair_within_cutoff_ortho(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    int ncell[3] = {
        acc->cell_ext[0] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[0]) : 0,
        acc->cell_ext[1] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[1]) : 0,
        acc->cell_ext[2] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[2]) : 0,
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_pair_within_cutoff_ortho: cutoff too large for cell size");
        return false;
    }

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
	const md_256 v_r2 = md_mm256_set1_ps((float)(cutoff * cutoff));

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
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t cmin[3] = { acc->cell_min[0], acc->cell_min[1], acc->cell_min[2] };
    const uint32_t cmax[3] = { acc->cell_max[0], acc->cell_max[1], acc->cell_max[2] };
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

                const md_256i v_len_i = md_mm256_set1_epi32(len_i);

                // --- Self cell: only j > i ---
                for (uint32_t i = 0; i < len_i - 1; ++i) {
                    POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_i - (i + 1), 8));

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
                        md_256 v_d2 = distance_squared_ortho_avx(v_dx, v_dy, v_dz, G00, G11, G22);

                        const md_256 v_mask = md_mm256_and_ps(md_mm256_cmplt_ps(v_d2, v_r2), j_mask);

						// Fill buffers with results
						int mask = md_mm256_movemask_ps(v_mask);

                        if (mask) {
						    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_i_idx + j));

                            uint64_t key = avx_compress_lut[mask];
                            const md_256i v_perm_mask = simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));

						    v_idxj = simde_mm256_permutevar8x32_epi32(v_idxj,  v_perm_mask);
						    v_d2   = simde_mm256_permutevar8x32_ps(v_d2,       v_perm_mask);

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

                    const float shift_xf = (float)simde_mm_extract_epi32(shift_i, 0);
                    const float shift_yf = (float)simde_mm_extract_epi32(shift_i, 1);
                    const float shift_zf = (float)simde_mm_extract_epi32(shift_i, 2);

                    const md_256 shift_x = md_mm256_set1_ps(shift_xf);
                    const md_256 shift_y = md_mm256_set1_ps(shift_yf);
                    const md_256 shift_z = md_mm256_set1_ps(shift_zf);

                    const md_256i v_len_j = md_mm256_set1_epi32(len_j);

                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    for (uint32_t i = 0; i < len_i; ++i) {
                        POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_j, 8));

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
                            md_256 v_d2 = distance_squared_ortho_avx(v_dx, v_dy, v_dz, G00, G11, G22);

                            const md_256 v_mask = md_mm256_and_ps(md_mm256_cmplt_ps(v_d2, v_r2), j_mask);

                            // Fill buffers with results
                            int mask = md_mm256_movemask_ps(v_mask);

                            if (mask) {
                                md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                                uint64_t key = avx_compress_lut[mask];
                                const md_256i v_perm_mask = simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));

                                v_idxj = simde_mm256_permutevar8x32_epi32(v_idxj, v_perm_mask);
                                v_d2   = simde_mm256_permutevar8x32_ps(v_d2,      v_perm_mask);

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

    FLUSH_TAIL();

    return true;
}

static bool for_each_external_point_within_cutoff_triclinic(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, const int32_t* ext_idx, size_t ext_count, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    const int ncell[3] = {
        acc->cell_ext[0] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[0]) : 0,
        acc->cell_ext[1] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[1]) : 0,
        acc->cell_ext[2] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[2]) : 0,
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_external_pair_within_cutoff_ortho: cutoff too large for cell size");
        return false;
    }

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 H01 = md_mm256_set1_ps(acc->H01);
    const md_256 H02 = md_mm256_set1_ps(acc->H02);
    const md_256 H12 = md_mm256_set1_ps(acc->H12);
    const md_256 v_r2 = md_mm256_set1_ps((float)(cutoff * cutoff));

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
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t cmin[3] = { acc->cell_min[0], acc->cell_min[1], acc->cell_min[2] };
    const uint32_t cmax[3] = { acc->cell_max[0], acc->cell_max[1], acc->cell_max[2] };
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

    const mat4_t I4 = {
        (float)acc->I[0][0], (float)acc->I[0][1], (float)acc->I[0][2], 0,
        (float)acc->I[1][0], (float)acc->I[1][1], (float)acc->I[1][2], 0,
        (float)acc->I[2][0], (float)acc->I[2][1], (float)acc->I[2][2], 0,
        0, 0, 0, 0,
    };


    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const vec4_t pbc_mask = vec4_set((acc->flags & MD_UNITCELL_PBC_X) ? val : 0, (acc->flags & MD_UNITCELL_PBC_Y) ? val : 0, (acc->flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);
	const vec4_t coord_offset = { acc->origin[0], acc->origin[1], acc->origin[2], 0 };
	const vec4_t fcell_dim = vec4_set((float)cdim[0], (float)cdim[1], (float)cdim[2], 0);

    for (size_t ei = 0; ei < ext_count; ++ei) {
        uint32_t idx = ext_idx ? (uint32_t)ext_idx[ei] : (uint32_t)ei;

        vec4_t r = {ext_x[idx], ext_y[idx], ext_z[idx], 0};
        r = vec4_add(r, coord_offset);
		vec4_t f = mat4_mul_vec4(I4, r);
		f = vec4_blend(f, vec4_fract(f), pbc_mask);
		ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f, fcell_dim)));

        const md_256i v_idxi = md_mm256_set1_epi32((uint32_t)idx);

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

            POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_j, 8));

            md_256i v_j = inc;
            for (uint32_t j = 0; j < len_j; j += 8) {
                const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));

                const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                md_256 v_d2 = distance_squared_tric_avx(v_dx, v_dy, v_dz, G00, G11, G22, H01, H02, H12);

                const md_256 v_mask = md_mm256_and_ps(md_mm256_cmplt_ps(v_d2, v_r2), j_mask);

                // Fill buffers with results
                int mask = md_mm256_movemask_ps(v_mask);

                if (mask) {
                    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                    uint64_t key = avx_compress_lut[mask];
                    const md_256i v_perm_mask = simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));

                    v_idxj = simde_mm256_permutevar8x32_epi32(v_idxj, v_perm_mask);
                    v_d2   = simde_mm256_permutevar8x32_ps(v_d2,      v_perm_mask);

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

    FLUSH_TAIL();

    return true;
}

static bool for_each_external_point_within_cutoff_ortho(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, const int32_t* ext_idx, size_t ext_count, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    const int ncell[3] = {
        acc->cell_ext[0] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[0]) : 0,
        acc->cell_ext[1] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[1]) : 0,
        acc->cell_ext[2] > 0 ? (int)ceil(cutoff / (double)acc->cell_ext[2]) : 0,
    };

    if (2 * ncell[0] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[1] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS || 2 * ncell[2] + 1 > SPATIAL_ACC_MAX_NEIGHBOR_CELLS) {
        MD_LOG_ERROR("for_each_external_pair_within_cutoff_ortho: cutoff too large for cell size");
        return false;
    }

    // Precompute constants
    const md_256 G00 = md_mm256_set1_ps(acc->G00);
    const md_256 G11 = md_mm256_set1_ps(acc->G11);
    const md_256 G22 = md_mm256_set1_ps(acc->G22);
    const md_256 v_r2 = md_mm256_set1_ps((float)(cutoff * cutoff));

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
    const uint32_t* cell_offset = acc->cell_off;
    const uint32_t cdim[3] = { acc->cell_dim[0], acc->cell_dim[1], acc->cell_dim[2] };
    const uint32_t cmin[3] = { acc->cell_min[0], acc->cell_min[1], acc->cell_min[2] };
    const uint32_t cmax[3] = { acc->cell_max[0], acc->cell_max[1], acc->cell_max[2] };
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

    const mat4_t I4 = {
        (float)acc->I[0][0], (float)acc->I[0][1], (float)acc->I[0][2], 0,
        (float)acc->I[1][0], (float)acc->I[1][1], (float)acc->I[1][2], 0,
        (float)acc->I[2][0], (float)acc->I[2][1], (float)acc->I[2][2], 0,
        0, 0, 0, 0,
    };


    float val;
    MEMSET(&val, 0xFF, sizeof(val));
    const vec4_t pbc_mask = vec4_set((acc->flags & MD_UNITCELL_PBC_X) ? val : 0, (acc->flags & MD_UNITCELL_PBC_Y) ? val : 0, (acc->flags & MD_UNITCELL_PBC_Z) ? val : 0, 0);
	const vec4_t coord_offset = { acc->origin[0], acc->origin[1], acc->origin[2], 0 };
	const vec4_t fcell_dim = vec4_set((float)cdim[0], (float)cdim[1], (float)cdim[2], 0);

    for (size_t ei = 0; ei < ext_count; ++ei) {
		uint32_t idx = ext_idx ? (uint32_t)ext_idx[ei] : (uint32_t)ei;

        vec4_t r = {ext_x[idx], ext_y[idx], ext_z[idx], 0};
        r = vec4_add(r, coord_offset);
		vec4_t f = mat4_mul_vec4(I4, r);
		f = vec4_blend(f, vec4_fract(f), pbc_mask);
		ivec4_t c_v = ivec4_from_vec4(vec4_floor(vec4_mul(f, fcell_dim)));

        const md_256i v_idxi = md_mm256_set1_epi32((uint32_t)idx);

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

            POSSIBLY_INVOKE_CALLBACK(ALIGN_TO(len_j, 8));

            md_256i v_j = inc;
            for (uint32_t j = 0; j < len_j; j += 8) {
                const md_256 j_mask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(v_j, v_len_j));

                const md_256 v_dx = md_mm256_sub_ps(v_xi, md_mm256_loadu_ps(elem_j_x + j));
                const md_256 v_dy = md_mm256_sub_ps(v_yi, md_mm256_loadu_ps(elem_j_y + j));
                const md_256 v_dz = md_mm256_sub_ps(v_zi, md_mm256_loadu_ps(elem_j_z + j));
                md_256 v_d2 = distance_squared_ortho_avx(v_dx, v_dy, v_dz, G00, G11, G22);

                const md_256 v_mask = md_mm256_and_ps(md_mm256_cmplt_ps(v_d2, v_r2), j_mask);

                // Fill buffers with results
                int mask = md_mm256_movemask_ps(v_mask);

                if (mask) {
                    md_256i v_idxj = md_mm256_loadu_si256((const md_256i*)(elem_j_idx + j));

                    uint64_t key = avx_compress_lut[mask];
                    const md_256i v_perm_mask = simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));

                    v_idxj = simde_mm256_permutevar8x32_epi32(v_idxj, v_perm_mask);
                    v_d2   = simde_mm256_permutevar8x32_ps(v_d2,      v_perm_mask);

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

    FLUSH_TAIL();

    return true;
}

#undef SPATIAL_ACC_BUFLEN
#undef CELL_INDEX
#undef CELL_OFFSET
#undef CELL_LENGTH

void md_spatial_acc_for_each_pair_in_neighboring_cells(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param) {
    ASSERT(acc);
	ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        for_each_internal_pair_in_neighboring_cells_triclinic(acc, callback, user_param);
    } else {
		for_each_internal_pair_in_neighboring_cells_ortho(acc, callback, user_param);
	}
}

bool md_spatial_acc_for_each_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    ASSERT(acc);
    ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        return for_each_internal_pair_within_cutoff_triclinic(acc, cutoff, callback, user_param);
    } else {
        return for_each_internal_pair_within_cutoff_ortho(acc, cutoff, callback, user_param);
	}
}

// Iterate over external points against points within the spatial acceleration structure for a supplied cutoff
// The external points are not part of the spatial acceleration structure and will be represented in the callback as the 'i' indices and the internal points are the 'j' indices
bool md_spatial_acc_for_each_external_point_within_cutoff(const md_spatial_acc_t* acc, const float* ext_x, const float* ext_y, const float* ext_z, const int32_t* ext_idx, size_t ext_count, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param) {
    ASSERT(acc);
	ASSERT(ext_x);
    ASSERT(ext_y);
    ASSERT(ext_z);
    ASSERT(callback);

    if (acc->flags & MD_UNITCELL_TRICLINIC) {
        return for_each_external_point_within_cutoff_triclinic(acc, ext_x, ext_y, ext_z, ext_idx, ext_count, cutoff, callback, user_param);
    } else {
        return for_each_external_point_within_cutoff_ortho(acc, ext_x, ext_y, ext_z, ext_idx, ext_count, cutoff, callback, user_param);
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
