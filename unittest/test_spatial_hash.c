#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_spatial_hash.h>
#include <core/md_str.h>
#include <core/md_intrinsics.h>
#include <core/md_os.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_molecule.h>
#include <md_util.h>

static bool iter_fn(const md_spatial_hash_elem_t* elem, void* user_param) {
    (void)elem;
    uint32_t* count = (uint32_t*)user_param;
    *count += 1;
    return true;
}

typedef struct {
    uint32_t exclude_idx;
    uint32_t count;
} iter_excl_data_t;

static bool iter_excl_fn(const md_spatial_hash_elem_t* elem, void* user_param) {
    (void)elem;
    iter_excl_data_t* data = user_param;
    if (elem->idx != data->exclude_idx) {
        data->count += 1;
    }
    return true;
}

static bool iter_batch_fn(const md_spatial_hash_elem_t* elem, int mask, void* user_param) {
    (void)elem;
    uint32_t* count = (uint32_t*)user_param;
    *count += popcnt32(mask);
    return true;
}

typedef struct {
    uint32_t i, j;
    double d2;
} dist_pair_t;

typedef struct {
    dist_pair_t* pairs;
    uint32_t count;
} iter_batch_pair_data_t;

static bool iter_batch_pairs_fn(const md_spatial_hash_elem_t* elem, md_256 d2, int mask, size_t i, void* user_param) {
    iter_batch_pair_data_t* data = user_param;
    
    float d2_raw[8];
    md_mm256_storeu_ps(d2_raw, d2);

    while (mask) {
        const int idx = ctz32(mask);
        mask = mask & ~(1 << idx);
        uint32_t j = elem[idx].idx;
        data->pairs[data->count].i = MIN((uint32_t)i,j);
        data->pairs[data->count].j = MAX((uint32_t)i,j);
        data->pairs[data->count].d2 = d2_raw[idx];
        data->count += 1;
    }

    return true;
}

UTEST(spatial_hash, small_periodic) {
    float x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    float y[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    md_unit_cell_t unit_cell = md_util_unit_cell_from_extent(10, 0, 0);
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(x, y, z, NULL, 10, &unit_cell, md_get_heap_allocator());
    ASSERT_TRUE(spatial_hash);
    
    uint32_t count = 0;
    md_spatial_hash_query(spatial_hash, (vec3_t){5,0,0}, 1.5f, iter_fn, &count);
    EXPECT_EQ(3, count);

    count = 0;
    md_spatial_hash_query(spatial_hash, (vec3_t){8.5f, 0, 0}, 3, iter_fn, &count);
    EXPECT_EQ(6, count);

    md_spatial_hash_free(spatial_hash);
}

UTEST(spatial_hash, big) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    const str_t pdb_file = STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &mol.unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);

    vec3_t pos = {24, 48, 24};
    float  rad = 10.0f;
    uint32_t count = 0;
    md_spatial_hash_query(spatial_hash, pos, rad, iter_fn, &count);

    EXPECT_NE(0, count);

    md_arena_allocator_destroy(alloc);
}

typedef struct elem_t {
    float x, y, z;
    uint32_t idx;
} elem_t;

typedef struct data_t {
    md_256 rx, ry, rz, r2;
    const elem_t* elem;
    const uint32_t* cell_offset;
} data_t;

static inline size_t test_elem_scalar(float x, float y, float z, float r2, const elem_t* elem, uint32_t len) {
    size_t result = 0;
    for (uint32_t i = 0; i < len; ++i) {
        float dx = elem[i].x - x;
        float dy = elem[i].y - y;
        float dz = elem[i].z - z;
        float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < r2) {
            result++;
        }
    }
    return result;
}

static inline size_t test_elem(float x, float y, float z, float r2, const elem_t* elem, uint32_t len) {
    size_t result = 0;

    md_256 rx = md_mm256_set1_ps(x);
    md_256 ry = md_mm256_set1_ps(y);
    md_256 rz = md_mm256_set1_ps(z);
    md_256 r2_vec = md_mm256_set1_ps(r2);

    for (uint32_t i = 0; i < len; i += 8) {
        md_256 vx, vy, vz;
        md_mm256_unpack_xyz_ps(&vx, &vy, &vz, (const float*)(elem + i), sizeof(elem_t));
        md_256 dx = md_mm256_sub_ps(vx, rx);
        md_256 dy = md_mm256_sub_ps(vy, ry);
        md_256 dz = md_mm256_sub_ps(vz, rz);
        md_256 d2 = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(dx, dx), md_mm256_mul_ps(dy, dy)), md_mm256_mul_ps(dz, dz));
        md_256 vmask = md_mm256_cmplt_ps(d2, r2_vec);

        const uint32_t remainder = MIN(len - i, 8);
        const uint32_t lane_mask = (1U << remainder) - 1U;
        const uint32_t mask = md_mm256_movemask_ps(vmask) & lane_mask;

        result += popcnt32(mask);
    }

    return result;
}

typedef void (*n2_callback_fn)(const float* j_x, const float* j_y, const float* j_z, const uint32_t* j_idx, md_256 d2, int mask, uint32_t i,
                               void* user_param);

static inline void test_elem_frac_ortho_simd2(md_256 x, md_256 y, md_256 z, md_256 r2,
                                              const float* elem_x, const float* elem_y, const float* elem_z, const uint32_t* elem_idx, int len,
                                              md_256 G00, md_256 G11, md_256 G22, md_256 H01, md_256 H02, md_256 H12,
                                              uint32_t i, dist_pair_t* out_pairs, size_t* count) {
    while (len > 0) {
        //md_256 vx, vy, vz;
        //md_mm256_unpack_xyz_ps(&vx, &vy, &vz, (const float*)elem, sizeof(elem_t));
        md_256 vx = md_mm256_loadu_ps(elem_x);
        md_256 vy = md_mm256_loadu_ps(elem_y);
        md_256 vz = md_mm256_loadu_ps(elem_z);

        md_256 dx = md_mm256_sub_ps(vx, x);
        md_256 dy = md_mm256_sub_ps(vy, y);
        md_256 dz = md_mm256_sub_ps(vz, z);

        // Orthonormal / diagonal G: d2 = dx^2 + dy^2 + dz^2 (or scaled diag)
        // If truly identity, you can drop the coeff multiplies entirely.
        md_256 dx2 = md_mm256_mul_ps(dx, dx);
        md_256 dy2 = md_mm256_mul_ps(dy, dy);
        md_256 dz2 = md_mm256_mul_ps(dz, dz);
        md_256 d2 = md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));

        md_256 vmask = md_mm256_cmplt_ps(d2, r2);

        const int remainder = MIN(len, 8);
        const int lane_mask = (1U << remainder) - 1U;
        const int mask = md_mm256_movemask_ps(vmask) & lane_mask;

        if (mask) {
            float d2_raw[8];
            md_mm256_storeu_ps(d2_raw, d2);
            int bits = mask;
            while (bits) {
                const int idx = ctz32(bits);
                bits = bits & ~(1 << idx);
                uint32_t j = elem_idx[idx];
                out_pairs[*count].i = MIN(i, j);
                out_pairs[*count].j = MAX(i, j);
                out_pairs[*count].d2 = d2_raw[idx];
                *count += 1;
            }
        }

        len -= 8;
        elem_x += 8;
        elem_y += 8;
        elem_z += 8;
        elem_idx += 8;
    };
}

static inline void test_elem_frac_ortho_simd(md_256 x, md_256 y, md_256 z, md_256 r2, const elem_t* elem, int len, md_256 G00, md_256 G11, md_256 G22,
                                             md_256 H01, md_256 H02, md_256 H12, uint32_t i, dist_pair_t* pairs, size_t* count) {
    while (len > 0) {
        md_256 vx, vy, vz;
        md_mm256_unpack_xyz_ps(&vx, &vy, &vz, (const float*)elem, sizeof(elem_t));

        md_256 dx = md_mm256_sub_ps(vx, x);
        md_256 dy = md_mm256_sub_ps(vy, y);
        md_256 dz = md_mm256_sub_ps(vz, z);

        // Orthonormal / diagonal G: d2 = dx^2 + dy^2 + dz^2 (or scaled diag)
        // If truly identity, you can drop the coeff multiplies entirely.
        md_256 dx2 = md_mm256_mul_ps(dx, dx);
        md_256 dy2 = md_mm256_mul_ps(dy, dy);
        md_256 dz2 = md_mm256_mul_ps(dz, dz);
        md_256 d2 = md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));

        md_256 vmask = md_mm256_cmplt_ps(d2, r2);

        const int remainder = MIN(len, 8);
        const int lane_mask = (1U << remainder) - 1U;
        const int mask = md_mm256_movemask_ps(vmask) & lane_mask;

        size_t cnt = popcnt32(mask);
        if (pairs && mask) {
            int bits = mask;
            float d2_raw[8];
            md_mm256_storeu_ps(d2_raw, d2);

            while (bits) {
                const int idx = ctz32(bits);
                bits = bits & ~(1 << idx);
                uint32_t j = elem[idx].idx;
                pairs[*count].i = MIN(i, j);
                pairs[*count].j = MAX(i, j);
                pairs[*count].d2 = d2_raw[idx];
                *count += 1;
            }
        } else {
            *count += cnt;
        }

        len -= 8;
        elem += 8;
    };
}

static inline void test_elem_frac_simd(md_256 x, md_256 y, md_256 z, md_256 r2, const elem_t* elem, int len,
    md_256 G00, md_256 G11, md_256 G22,
    md_256 H01, md_256 H02, md_256 H12,
    uint32_t i, dist_pair_t* pairs, size_t* count)
{
    while (len > 0) {
        md_256 vx, vy, vz;
        md_mm256_unpack_xyz_ps(&vx, &vy, &vz, (const float*)elem, sizeof(elem_t));

        md_256 dx = md_mm256_sub_ps(vx, x);
        md_256 dy = md_mm256_sub_ps(vy, y);
        md_256 dz = md_mm256_sub_ps(vz, z);

        // Full symmetric quadratic form:
        // d2 = G00*dx^2 + G11*dy^2 + G22*dz^2 + H01*(dx*dy) + H02*(dx*dz) + H12*(dy*dz)
        md_256 dx2 = md_mm256_mul_ps(dx, dx);
        md_256 dy2 = md_mm256_mul_ps(dy, dy);
        md_256 dz2 = md_mm256_mul_ps(dz, dz);

        md_256 dxy = md_mm256_mul_ps(dx, dy);
        md_256 dxz = md_mm256_mul_ps(dx, dz);
        md_256 dyz = md_mm256_mul_ps(dy, dz);

        md_256 acc   = md_mm256_fmadd_ps(G00, dx2, md_mm256_fmadd_ps(G11, dy2, md_mm256_mul_ps(G22, dz2)));
        md_256 cross = md_mm256_fmadd_ps(H01, dxy, md_mm256_fmadd_ps(H02, dxz, md_mm256_mul_ps(H12, dyz)));
        md_256 d2 = md_mm256_add_ps(acc, cross);

        md_256 vmask = md_mm256_cmplt_ps(d2, r2);

        const int remainder = MIN(len, 8);
        const int lane_mask = (1U << remainder) - 1U;
        const int mask = md_mm256_movemask_ps(vmask) & lane_mask;

        size_t cnt = popcnt32(mask);
        if (pairs && mask) {
            int bits = mask;
            float d2_raw[8];
            md_mm256_storeu_ps(d2_raw, d2);

            while (bits) {
                const int idx = ctz32(bits);
                bits = bits & ~(1 << idx);
                uint32_t j = elem[idx].idx;
                pairs[*count].i  = MIN(i, j);
                pairs[*count].j  = MAX(i, j);
                pairs[*count].d2 = d2_raw[idx];
                *count += 1;
            }
        } else {
            *count += cnt;
        }

        len  -= 8;
        elem += 8;
    };
}

// -----------------------------------------------------------------------------
// Neighbor topology (forward neighbors for self+13 scheme) and wrapping helper
// -----------------------------------------------------------------------------

static const int8_t FWD_NBRS[13][3] = {
    {1, 0, 0},  {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},  {-1, -1, 1}, {0, -1, 1}, {1, -1, 1},
    {-1, 0, 1}, {0, 0, 1},  {1, 0, 1}, {-1, 1, 1}, {0, 1, 1},   {1, 1, 1},
};

// Wrap one axis with optional PBC, returning wrapped index and image shift k∈{-1,0,+1}.
// If axis is non-periodic and n is OOB, returns -1 and k=0.
static inline int wrap_and_shift(int n, int dim, int periodic, int* k) {
    *k = 0;
    if (!periodic) {
        if (n < 0 || n >= dim) return -1;
        return n;
    }
    if (n < 0) {
        n += dim;
        *k = -1;
    } else if (n >= dim) {
        n -= dim;
        *k = +1;
    }
    return n;
}

// -----------------------------------------------------------------------------
// Main: count pairs within cutoff using cell list in fractional space
// -----------------------------------------------------------------------------

typedef md_128i ivec4_t;

static inline ivec4_t ivec4_set(int x, int y, int z, int w) {
    return simde_mm_set_epi32(w, z, y, x);
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

static size_t do_pairwise_periodic_triclinic(const float* in_x, const float* in_y, const float* in_z, size_t num_points, double cutoff,
                                             const md_unit_cell_t* unit_cell, dist_pair_t* out_pairs) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));  // 1 GiB arena

    mat3_t A = {0};
    mat3_t Ai = {0};
    uint32_t flags = 0;

    if (unit_cell) {
        A  = unit_cell->basis;
        Ai = unit_cell->inv_basis;
        flags = unit_cell->flags;
    }

    vec4_t coord_offset = vec4_zero();

    const int PBC_ALL = (MD_UNIT_CELL_FLAG_PBC_X | MD_UNIT_CELL_FLAG_PBC_Y | MD_UNIT_CELL_FLAG_PBC_Z);
    if (flags & PBC_ALL != PBC_ALL) {
        if (unit_cell) {
            // We cannot have a triclinic system that is partially periodic
            ASSERT(flags & MD_UNIT_CELL_FLAG_ORTHO);
        }
        // Unit cell either missing or not periodic along one or more axis
        vec4_t aabb_min = {0}, aabb_max = {0};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, in_x, in_y, in_z, NULL, NULL, num_points);

        // Construct A and Ai from aabb extent
        vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);

        if (!(flags & MD_UNIT_CELL_FLAG_PBC_X)) {
            coord_offset.x = -aabb_min.x;
            if (aabb_ext.x > 0.0f) {
                A.elem[0][0]  = aabb_ext.x;
                Ai.elem[0][0] = 1.0f / aabb_ext.x;
            }
        }
        if (!(flags & MD_UNIT_CELL_FLAG_PBC_Y)) {
            coord_offset.y = -aabb_min.y;
            if (aabb_ext.y > 0.0f) {
                A.elem[1][1]  = aabb_ext.y;
                Ai.elem[1][1] = 1.0f / aabb_ext.y;
            }
        }
        if (!(flags & MD_UNIT_CELL_FLAG_PBC_Z)) {
            coord_offset.z = -aabb_min.z;
            if (aabb_ext.z > 0.0f) {
                A.elem[2][2]  = aabb_ext.z;
                Ai.elem[2][2] = 1.0f / aabb_ext.z;
            }
        }
    }

    // Precompute metric G = A^T A
    const mat3_t G = mat3_mul(mat3_transpose(A), A);
    const mat4x3_t G43 = mat4x3_from_mat3(G);

    // Choose grid resolution. Heuristic:  cell_dim ≈ |A|/CELL_EXT.
    // Using ~cutoff-ish spacing gives good pruning. Tweak CELL_EXT if you like.
    const double CELL_EXT = cutoff > 0.0 ? cutoff : 6.0;

    // Estimate cell_dim by measuring the extents of the box vectors (norms of columns of A)
    // This is only a heuristic for bin counts; the grid is still in fractional space.
    double ax = sqrt((double)A.elem[0][0] * (double)A.elem[0][0] + (double)A.elem[1][0] * (double)A.elem[1][0] + (double)A.elem[2][0] * (double)A.elem[2][0]);
    double by = sqrt((double)A.elem[0][1] * (double)A.elem[0][1] + (double)A.elem[1][1] * (double)A.elem[1][1] + (double)A.elem[2][1] * (double)A.elem[2][1]);
    double cz = sqrt((double)A.elem[0][2] * (double)A.elem[0][2] + (double)A.elem[1][2] * (double)A.elem[1][2] + (double)A.elem[2][2] * (double)A.elem[2][2]);

    uint32_t cell_dim[3] = {
        CLAMP((uint32_t)(ax / CELL_EXT + 0.5), 1, 1024),
        CLAMP((uint32_t)(by / CELL_EXT + 0.5), 1, 1024),
        CLAMP((uint32_t)(cz / CELL_EXT + 0.5), 1, 1024),
    };

    const uint32_t c0  = cell_dim[0];
    const uint32_t c01 = cell_dim[0] * cell_dim[1];
    const size_t num_cells = (size_t)cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Temporary arrays
    uint32_t* cell_offset = (uint32_t*)md_vm_arena_push_zero(arena, (num_cells + 1) * sizeof(uint32_t));
    uint32_t* local_idx = (uint32_t*)md_vm_arena_push(arena, num_points * sizeof(uint32_t));
    uint32_t* cell_idx = (uint32_t*)md_vm_arena_push(arena, num_points * sizeof(uint32_t));
    elem_t* scratch_s = (elem_t*)md_vm_arena_push(arena, num_points * sizeof(elem_t));  // holds fractional coords before scatter

    md_timestamp_t t0 = md_time_current();

    const vec4_t  fcell_dim = vec4_set(cell_dim[0], cell_dim[1], cell_dim[2], 0);
    const vec4_t  fcell_max = vec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);
    const ivec4_t icell_min = ivec4_set1(0);
    const ivec4_t icell_max = ivec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);

    const ivec4_t c = ivec4_set(1, c0, c01, 0);

    ivec4_t cell_min = ivec4_set1(INT32_MAX);
    ivec4_t cell_max = ivec4_set1(INT32_MIN);

    mat4x3_t A43  = mat4x3_from_mat3(A);
    mat4x3_t Ai43 = mat4x3_from_mat3(Ai);

    // 1) Convert to fractional, wrap periodic axes into [0,1), bin to cells
    for (size_t i = 0; i < num_points; ++i) {
        vec4_t r = {in_x[i], in_y[i], in_z[i], 0};
        r = vec4_add(r, coord_offset);

        vec4_t s = mat4x3_mul_vec4(Ai43, r);
        s = vec4_fract(s);

        vec4_t f = vec4_mul(s, fcell_dim);
        f = vec4_floor(f);
        ivec4_t ic = ivec4_from_vec4(f);

        ic = ivec4_clamp(ic, icell_min, icell_max);

        cell_min = ivec4_min(cell_min, ic);
        cell_max = ivec4_max(cell_max, ic);

        uint32_t cell_coord[4];
        md_mm_storeu_epi32(cell_coord, ic);
        uint64_t ci = (uint64_t)cell_coord[2] * c01 + (uint64_t)cell_coord[1] * c0 + (uint64_t)cell_coord[0];

        ASSERT(ci < num_cells);

        local_idx[i] = cell_offset[ci]++;  // count for now
        cell_idx[i] = (uint32_t)ci;

        // stash fractional coordinates
        MEMCPY(&scratch_s[i], &s, sizeof(elem_t));
        scratch_s[i].idx = (uint32_t)i;
    }

    md_timestamp_t t1 = md_time_current();
    printf("1: %.2f\n", md_time_as_milliseconds(t1 - t0));

    // 2) Prefix sum cell offsets
    uint32_t running = 0;
    for (size_t ci = 0; ci <= num_cells; ++ci) {
        uint32_t len = cell_offset[ci];
        cell_offset[ci] = running;
        running += len;
    }

    size_t alloc_len = running + 8;
    float*    element_x = md_vm_arena_push(arena, alloc_len * sizeof(float));
    float*    element_y = md_vm_arena_push(arena, alloc_len * sizeof(float));
    float*    element_z = md_vm_arena_push(arena, alloc_len * sizeof(float));
    uint32_t* element_i = md_vm_arena_push(arena, alloc_len * sizeof(uint32_t));
    elem_t* elements = md_vm_arena_push(arena, alloc_len * sizeof(elem_t));  // sorted fractional coords

    md_timestamp_t t2 = md_time_current();
    printf("2: %.2f\n", md_time_as_milliseconds(t2 - t1));

    // 3) Scatter fractional coords into 'elements' in cell order
    for (size_t i = 0; i < num_points; ++i) {
        uint32_t dst = cell_offset[cell_idx[i]] + local_idx[i];
        //elements[dst] = scratch_s[i];
        element_x[dst] = scratch_s[i].x;
        element_y[dst] = scratch_s[i].y;
        element_z[dst] = scratch_s[i].z;
        element_i[dst] = scratch_s[i].idx;
    }

    md_timestamp_t t3 = md_time_current();
    printf("3: %.2f\n", md_time_as_milliseconds(t3 - t2));

    // 5) Pair counting
    const float r2_cut = cutoff * cutoff;
    size_t count = 0;

// Macros for cell offsets/lengths in the sorted array
#define CELL_INDEX(x, y, z) ((uint32_t)z * c01 + (uint32_t)y * c0 + (uint32_t)x)
#define CELL_OFFSET(ci) (cell_offset[(ci)])
#define CELL_LENGTH(ci) (cell_offset[(ci) + 1] - cell_offset[(ci)])

    int cmin[4], cmax[4];
    MEMCPY(cmin, &cell_min, sizeof(cmin));
    MEMCPY(cmax, &cell_max, sizeof(cmax));

    const md_256 G00 = md_mm256_set1_ps(G.elem[0][0]);
    const md_256 G11 = md_mm256_set1_ps(G.elem[1][1]);
    const md_256 G22 = md_mm256_set1_ps(G.elem[2][2]);

    const md_256 H01 = md_mm256_set1_ps(2.0 * G.elem[0][1]);
    const md_256 H02 = md_mm256_set1_ps(2.0 * G.elem[0][2]);
    const md_256 H12 = md_mm256_set1_ps(2.0 * G.elem[1][2]);

    const md_256 r2  = md_mm256_set1_ps(r2_cut);

    for (uint32_t cz = cmin[2]; cz <= cmax[2]; ++cz) {
        for (uint32_t cy = cmin[1]; cy <= cmax[1]; ++cy) {
            for (uint32_t cx = cmin[0]; cx <= cmax[0]; ++cx) {
                const uint32_t ci    = CELL_INDEX(cx, cy, cz);
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);
                //const elem_t* elem_i = elements + off_i;
                const float* elem_i_x      = element_x + off_i;
                const float* elem_i_y      = element_y + off_i;
                const float* elem_i_z      = element_z + off_i;
                const uint32_t* elem_i_idx = element_i + off_i;

                if (len_i == 0) continue;

                // Self cell: only j > i
                for (uint32_t a = 0; a < len_i; ++a) {
                    const md_256 x = md_mm256_set1_ps(elem_i_x[a]);
                    const md_256 y = md_mm256_set1_ps(elem_i_y[a]);
                    const md_256 z = md_mm256_set1_ps(elem_i_z[a]);
                    uint32_t i = elem_i_idx[a];

                    uint32_t off = a + 1;
                    uint32_t len = len_i - off;
                    
                    //test_elem_frac_ortho_simd(x, y, z, r2, elem_i + (a + 1), len_i - (a + 1), G00, G11, G22, H01, H02, H12, i, pairs, &result);
                    test_elem_frac_ortho_simd2(x, y, z, r2, elem_i_x + off, elem_i_y + off, elem_i_z + off, elem_i_idx + off, len, G00, G11, G22, H01,
                                               H02, H12, i, out_pairs, &count);
                }

                // Forward neighbors
                for (int n = 0; n < 13; ++n) {
                    int ndx = FWD_NBRS[n][0];
                    int ndy = FWD_NBRS[n][1];
                    int ndz = FWD_NBRS[n][2];

                    int nx = (int)cx + ndx;
                    int ny = (int)cy + ndy;
                    int nz = (int)cz + ndz;

                    int pnx = (nx + cell_dim[0]) % cell_dim[0];
                    int pny = (ny + cell_dim[1]) % cell_dim[1];
                    int pnz = (nz + cell_dim[2]) % cell_dim[2];

                    float px = 0.0f;
                    if (pnx < nx) px = -1.0f;
                    if (pnx > nx) px =  1.0f;

                    float py = 0.0f;
                    if (pny < ny) py = -1.0f;
                    if (pny > ny) py =  1.0f;

                    float pz = 0.0f;
                    if (pnz < nz) pz = -1.0f;
                    if (pnz > nz) pz =  1.0f;

                    const uint32_t cj    = CELL_INDEX(pnx, pny, pnz);
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);
                    const elem_t* elem_j = elements + off_j;
                    const float* elem_j_x = element_x + off_j;
                    const float* elem_j_y = element_y + off_j;
                    const float* elem_j_z = element_z + off_j;
                    const uint32_t* elem_j_idx = element_i + off_j;

                    if (len_j == 0) continue;

                    for (uint32_t a = 0; a < len_i; ++a) {
                        const md_256 x = md_mm256_set1_ps(elem_i_x[a] + px);
                        const md_256 y = md_mm256_set1_ps(elem_i_y[a] + py);
                        const md_256 z = md_mm256_set1_ps(elem_i_z[a] + pz);
                        const uint32_t i = elem_i_idx[a];

                        //test_elem_frac_ortho_simd(x, y, z, r2, elem_j, len_j, G00, G11, G22, H01, H02, H12, i, pairs, &result);
                        test_elem_frac_ortho_simd2(x, y, z, r2, elem_j_x, elem_j_y, elem_j_z, elem_j_idx, len_j, G00, G11, G22, H01, H02, H12, i, out_pairs, &count);
                    }
                }
            }
        }
    }

    md_timestamp_t t4 = md_time_current();
    printf("4: %.2f\n", md_time_as_milliseconds(t4 - t3));

#undef CELL_INDEX
#undef CELL_OFFSET
#undef CELL_LENGTH

    md_vm_arena_destroy(arena);
    return count;
}

static size_t do_brute_force(const float* in_x, const float* in_y, const float* in_z, size_t num_points, float cutoff, const md_unit_cell_t* unit_cell, dist_pair_t* pairs) {
    size_t count = 0;
    const float r2 = cutoff * cutoff;

    if (unit_cell && (unit_cell->flags & MD_UNIT_CELL_FLAG_PBC_ALL) == MD_UNIT_CELL_FLAG_PBC_ALL) {
        vec4_t ext = md_unit_cell_box_ext(unit_cell);
        vec4_t inv_ext = vec4_set(ext.x > 0.0f ? 1.0f / ext.x : 0.0f, ext.y > 0.0f ? 1.0f / ext.y : 0.0f, ext.z > 0.0f ? 1.0f / ext.z : 0.0f, 0.0f);

        for (size_t i = 0; i < num_points - 1; ++i) {
            float x = in_x[i];
            float y = in_y[i];
            float z = in_z[i];
            for (size_t j = i + 1; j < num_points; ++j) {
                float dx = in_x[j] - x;
                float dy = in_y[j] - y;
                float dz = in_z[j] - z;

                vec4_t d = vec4_min_image(vec4_set(dx, dy, dz, 0.0f), ext, inv_ext);
                float d2 = fmaf(d.x, d.x, fmaf(d.y, d.y, d.z * d.z));
                if (d2 <= r2) {
                    if (pairs) {
                        pairs[count].i = (uint32_t)i;
                        pairs[count].j = (uint32_t)j;
                        pairs[count].d2 = d2;
                    }
                    count += 1;
                }
            }
        }
    } else {
        for (size_t i = 0; i < num_points - 1; ++i) {
            float x = in_x[i];
            float y = in_y[i];
            float z = in_z[i];
            for (size_t j = i + 1; j < num_points; ++j) {
                float dx = in_x[j] - x;
                float dy = in_y[j] - y;
                float dz = in_z[j] - z;

                float d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= r2) {
                    if (pairs) {
                        pairs[count].i = (uint32_t)i;
                        pairs[count].j = (uint32_t)j;
                        pairs[count].d2 = d2;
                    }
                    count += 1;
                }
            }
        }
    }

    return count;
}

static inline double rnd_rng(double min, double max) {
    double r = ((double)rand() / (double)RAND_MAX);
    return r * (max - min) + min;
}

static int compare_dist_pair(void const* a, void const* b) {
    const dist_pair_t* pa = a;
    const dist_pair_t* pb = b;
    if (pa->i != pb->i) return pa->i - pb->i;
    if (pa->j != pb->j) return pa->j - pb->j;
    return 0;
}

static void n2_callback(const float* j_x, const float* j_y, const float* j_z, const uint32_t* j_idx, md_256 d2, int mask, uint32_t i, void* user_param) {
    size_t* count = (size_t*)user_param;
    *count += popcnt32(mask);
};

void insertion_sort_uint32(uint32_t* array, size_t n) {
    for (int64_t i = 1; i < (int64_t)n; i++) {
        int64_t j = i;
        while (j > 0 && array[j] < array[j - 1]) {
            uint32_t tmp = array[j];
            array[j] = array[j - 1];
            array[j - 1] = tmp;
            j--;
        }
    }
}

UTEST(spatial_hash, n2) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro_data, alloc));

    vec3_t* xyz = md_alloc(alloc, sizeof(vec3_t) * mol.atom.count);
    for (size_t i = 0; i < mol.atom.count; ++i) {
        xyz[i] = vec3_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]);
    }

    const vec4_t mask = md_unit_cell_pbc_mask(&mol.unit_cell);
    mat4x3_t I = mat4x3_from_mat3(mol.unit_cell.inv_basis);
    mat4x3_t M = mat4x3_from_mat3(mol.unit_cell.basis);

    {

#define test_count 2048
        float x[test_count];
        float y[test_count];
        float z[test_count];

        const double ext = 103.0;

        srand(0);
        for (size_t i = 0; i < test_count; ++i) {
            x[i] = rnd_rng(0.0, ext);
            y[i] = rnd_rng(0.0, ext);
            z[i] = rnd_rng(0.0, ext);
        }

        dist_pair_t bf_pairs[4096];
        dist_pair_t sh_pairs[4096];

        md_unit_cell_t test_cell = md_util_unit_cell_from_extent(ext, ext, ext);
        md_spatial_hash_t* sh = md_spatial_hash_create_soa(x, y, z, NULL, test_count, &test_cell, alloc);



        for (double rad = 4.0; rad <= 5.5; rad += 0.5) {
            size_t bf_count = do_brute_force(x, y, z, test_count, rad, &test_cell, bf_pairs);

            size_t sh_count = do_pairwise_periodic_triclinic(x, y, z, test_count, rad, &test_cell, sh_pairs);
            
            if (bf_count != sh_count) {
                qsort(bf_pairs, bf_count, sizeof(dist_pair_t), compare_dist_pair);
                qsort(sh_pairs, sh_count, sizeof(dist_pair_t), compare_dist_pair);
                
                printf("D2: %f\n", rad * rad);
                printf("Box ext: %f %f %f\n", ext, ext, ext);
                size_t i = 0, j = 0;
                for (; i < bf_count && j < sh_count;) {
                    int c = compare_dist_pair(&bf_pairs[i], &sh_pairs[j]);
                    if (c == 0) {
                        i++; j++;
                    } else if (c < 0) {
                        float x_i = x[bf_pairs[i].i];
                        float y_i = y[bf_pairs[i].i];
                        float z_i = z[bf_pairs[i].i];
                        float x_j = x[bf_pairs[i].j];
                        float y_j = y[bf_pairs[i].j];
                        float z_j = z[bf_pairs[i].j];

                        int cx_i = (int)floorf(x_i / rad);
                        int cy_i = (int)floorf(y_i / rad);
                        int cz_i = (int)floorf(z_i / rad);
                        int cx_j = (int)floorf(x_j / rad);
                        int cy_j = (int)floorf(y_j / rad);
                        int cz_j = (int)floorf(z_j / rad);

                        printf("BF only: %u %u %.5f\n", bf_pairs[i].i, bf_pairs[i].j, bf_pairs[i].d2);
                        printf("i coord: %f %f %f [%i,%i,%i]\n", x_i, y_i, z_i, cx_i, cy_i, cz_i);
                        printf("j coord: %f %f %f [%i,%i,%i]\n", x_j, y_j, z_j, cx_j, cy_j, cz_j);
                        i++;
                    } else {
                        printf("SH only: %u %u %.5f\n", sh_pairs[j].i, sh_pairs[j].j, sh_pairs[j].d2);
                        printf("i coord: %f %f %f\n", x[sh_pairs[j].i], y[sh_pairs[j].i], z[sh_pairs[j].i]);
                        printf("j coord: %f %f %f\n", x[sh_pairs[j].j], y[sh_pairs[j].j], z[sh_pairs[j].j]);
                        j++;
                    }
                }
                printf("WIERD\n");
            }
        }
        md_spatial_hash_free(sh);
#undef test_count
    }

#if 1
    md_unit_cell_t unit_cell = mol.unit_cell;
    // Clear pbc flags

    size_t expected_count = 3711880;
    if (false) {
        unit_cell.flags &= ~(MD_UNIT_CELL_FLAG_PBC_X | MD_UNIT_CELL_FLAG_PBC_Y | MD_UNIT_CELL_FLAG_PBC_Z);
        expected_count = 3701958;
    }

    // Custom implementation of pairwise periodic N^2
    dist_pair_t* pairs = md_alloc(alloc, sizeof(dist_pair_t) * 10000000);
    md_timestamp_t start = md_time_current();
    //size_t count = do_pairwise_periodic(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell);
    size_t custom_count = do_pairwise_periodic_triclinic(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell, pairs);
    md_timestamp_t end = md_time_current();
    printf("Custom: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, custom_count);
    if (custom_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, custom_count);
    }

    start = md_time_current();

    int* x_coords = md_alloc(alloc, sizeof(int) * mol.atom.count);
    for (size_t i = 0; i < mol.atom.count; ++i) {
        x_coords[i] = (uint32_t)(mol.atom.x[i] * 1000.0f);
    }
    insertion_sort_uint32(x_coords, mol.atom.count);

    end = md_time_current();
    printf("Radix sort: %f ms\n", md_time_as_milliseconds(end - start));
    
    // Current implementation of spatial hash for N^2
    
    start = md_time_current();
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &unit_cell, alloc);
    uint32_t count = 0;
    md_spatial_hash_query_multi_batch(spatial_hash, xyz, mol.atom.count, 5.0f, iter_batch_fn, &count);
    end = md_time_current();
    printf("Spatial hash: %f ms\n", md_time_as_milliseconds(end - start));
    size_t sh_count = count;
    ASSERT_TRUE(spatial_hash);
    EXPECT_EQ(expected_count, sh_count);
    if (sh_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, sh_count);
    }

    // Brute force
    start = md_time_current();
    size_t bf_count = do_brute_force(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell, NULL);
    end = md_time_current();
    printf("Brute force: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, bf_count);
    if (bf_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, bf_count);
    }
#endif

    md_arena_allocator_destroy(alloc);
}

struct spatial_hash {
    md_allocator_i* arena;
    md_molecule_t mol;
    vec3_t pbc_ext;
};

UTEST_F_SETUP(spatial_hash) {
    utest_fixture->arena = md_vm_arena_create(GIGABYTES(4));

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), utest_fixture->arena));
    ASSERT_TRUE(md_gro_molecule_init(&utest_fixture->mol, &gro_data, utest_fixture->arena));
    utest_fixture->pbc_ext = mat3_mul_vec3(utest_fixture->mol.unit_cell.basis, vec3_set1(1.f));
}

UTEST_F_TEARDOWN(spatial_hash) {
    md_vm_arena_destroy(utest_fixture->arena);
}

static inline float rnd() {
    return rand() / (float)RAND_MAX;
}

UTEST_F(spatial_hash, test_correctness_centered) {
    md_molecule_t* mol = &utest_fixture->mol;
    md_allocator_i* alloc = utest_fixture->arena;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, NULL, mol->atom.count, NULL, alloc);
    ASSERT_TRUE(spatial_hash);
    
    srand(31);
    
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            vec4_t c = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};

            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_ala) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &mol.unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(mol.unit_cell.basis, vec3_set1(1)), 0);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20.f;
        const float rad2 = radius * radius;

        int ref_count = 0;
        for (size_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_ala_vec3) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, alloc));

    vec3_t* xyz = md_array_create(vec3_t, mol.atom.count, alloc);
    for (int i = 0; i < mol.atom.count; ++i) {
        xyz[i] = vec3_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]);
    }

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_vec3(xyz, NULL, md_array_size(xyz), &mol.unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(mol.unit_cell.basis, vec3_set1(1)), 0);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_water) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro"), NULL, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, NULL, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(mol.unit_cell.basis, vec3_set1(1)), 0);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_periodic_centered) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    md_molecule_t* mol = &utest_fixture->mol;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, NULL, mol->atom.count, &mol->unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);
    
    srand(31);

    const int num_iter = 100;
    const vec4_t period = vec4_from_vec3(pbc_ext, 0);
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float radius = rnd() * 20.0f;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < mol->atom.count; ++i) {
            vec4_t c = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};

            if (vec4_periodic_distance_squared(vec4_from_vec3(pos, 0), c, period) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_periodic_water) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro"), NULL, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &mol.unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(mol.unit_cell.basis, vec3_set1(1)), 0);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_periodic_distance_squared(vec4_from_vec3(pos, 0), c, pbc_ext) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);

        if (count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, count);
        }
    }

    md_vm_arena_temp_end(temp);
}