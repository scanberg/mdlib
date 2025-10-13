#include "md_spatial_acc.h"

#include <core/md_log.h>
#include <md_util.h>
#include <md_types.h>

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

static inline ivec4_t ivec4_add(ivec4_t a, ivec4_t b) {
    return simde_mm_add_epi32(a, b);
}

static inline ivec4_t ivec4_sub(ivec4_t a, ivec4_t b) {
    return simde_mm_sub_epi32(a, b);
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

    float* element_x = md_alloc(alloc, alloc_len * sizeof(float));
    float* element_y = md_alloc(alloc, alloc_len * sizeof(float));
    float* element_z = md_alloc(alloc, alloc_len * sizeof(float));
    uint32_t* element_i = md_alloc(alloc, alloc_len * sizeof(uint32_t));

    uint32_t* cell_offset = (uint32_t*)md_alloc(alloc, (num_cells + 1) * sizeof(uint32_t));

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
        vec4_t r = {in_x[i], in_y[i], in_z[i], 0};
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
    uint32_t running = 0;
    for (size_t ci = 0; ci <= num_cells; ++ci) {
        uint32_t len = cell_offset[ci];
        cell_offset[ci] = running;
        running += len;
    }

    // 3) Scatter fractional coords into 'elements' in cell order
    for (size_t i = 0; i < count; ++i) {
        uint32_t dst = cell_offset[cell_idx[i]] + local_idx[i];
        element_x[dst] = scratch_s[i].x;
        element_y[dst] = scratch_s[i].y;
        element_z[dst] = scratch_s[i].z;
        element_i[dst] = scratch_s[i].idx;
    }

    int cmin[4];
    int cmax[4];
    ivec4_store(cmin, cell_min);
    ivec4_store(cmax, cell_max);

    acc->elem_x = element_x;
    acc->elem_y = element_y;
    acc->elem_z = element_z;
    acc->elem_idx = element_i;
    acc->num_elems = count;
    MEMCPY(acc->cell_min, cmin, sizeof(acc->cell_min));
    MEMCPY(acc->cell_max, cmax, sizeof(acc->cell_max));
    MEMCPY(acc->cell_dim, cell_dim, sizeof(acc->cell_dim));
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