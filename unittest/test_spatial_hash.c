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

// -----------------------------------------------------------------------------
// MIC + metric distance in fractional space
// -----------------------------------------------------------------------------

static inline float mic1(float ds, int periodic) {
    // Minimum image per-axis in fractional coordinates
    if (!periodic) return ds;  // no wrap on non-PBC axis
    ds -= roundf(ds);          // shift into [-0.5, 0.5)
    return ds;
}

static inline float dist2_metric_frac(float dsx, float dsy, float dsz, const float G[3][3]) {
    // d^2 = (Δs)^T G (Δs)
    float gx = G[0][0] * dsx + G[0][1] * dsy + G[0][2] * dsz;
    float gy = G[1][0] * dsx + G[1][1] * dsy + G[1][2] * dsz;
    float gz = G[2][0] * dsx + G[2][1] * dsy + G[2][2] * dsz;
    return dsx * gx + dsy * gy + dsz * gz;
}

static inline float d2_metric(vec4_t ds, mat4x3_t G) {
    vec4_t v = mat4x3_mul_vec4(G, ds);
    return vec4_dot(ds, v);
}

// Scalar kernels (drop in your AVX where these are called if you want)

static inline size_t test_elem_frac(float sx, float sy, float sz, const elem_t* elem, uint32_t len, const mat4x3_t* G, const int periodic[3],
                                    float r2_cut) {
    size_t cnt = 0;
    vec4_t s = vec4_set(sx, sy, sz, 0);
    for (uint32_t i = 0; i < len; ++i) {
        vec4_t t;
        MEMCPY(&t, elem + i, sizeof(vec4_t));
        t = vec4_blend_mask(t, vec4_zero(), MD_SIMD_BLEND_MASK(0, 0, 0, 1));
        vec4_t ds = vec4_sub(t, s);
        float d2 = d2_metric(ds, *G);
        cnt += (d2 <= r2_cut);
    }
    return cnt;
}

static inline size_t test_elem_frac_shift(float sx, float sy, float sz, int kx, int ky, int kz,  // image shift from neighbor wrap
                                          const elem_t* elem, uint32_t len, const mat3_t* G, const int periodic[3], float r2_cut) {
    size_t cnt = 0;
    for (uint32_t j = 0; j < len; ++j) {
        float dsx = (elem[j].x - sx) + (float)kx;
        float dsy = (elem[j].y - sy) + (float)ky;
        float dsz = (elem[j].z - sz) + (float)kz;
        // MIC after applying the integer image shift
        dsx = mic1(dsx, periodic[0]);
        dsy = mic1(dsy, periodic[1]);
        dsz = mic1(dsz, periodic[2]);
        float d2 = dist2_metric_frac(dsx, dsy, dsz, G);
        cnt += (d2 <= r2_cut);
    }
    return cnt;
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

static inline ivec4_t ivec4_clamp(ivec4_t v, ivec4_t min, ivec4_t max) {
    return simde_mm_max_epi32(simde_mm_min_epi32(v, max), min);
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

static inline ivec4_t ivec4_mul(ivec4_t a, ivec4_t b) {
    return simde_mm_mul_epi32(a, b);
}

static inline int ivec4_dot(ivec4_t a, ivec4_t b) {
    ivec4_t prod = ivec4_mul(a, b);
    // Shuffle and add to accumulate
    ivec4_t shuf = simde_mm_shuffle_epi32(prod, SIMDE_MM_SHUFFLE(2,3,0,1));
    ivec4_t sums = ivec4_add(prod, shuf);
    shuf = simde_mm_shuffle_epi32(sums, SIMDE_MM_SHUFFLE(1,0,3,2));
    sums = ivec4_add(sums, shuf);

    // Extract the result (all lanes now hold the dot product)
    return simde_mm_cvtsi128_si32(sums);
}

static size_t do_pairwise_periodic_triclinic(const float* in_x, const float* in_y, const float* in_z, size_t num_points, float cutoff,const md_unit_cell_t* unit_cell) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));  // 1 GiB arena

    mat4x3_t A = {0};
    mat4x3_t Ai = {0};
    uint32_t flags = 0;

    if (unit_cell) {
        A  = mat4x3_from_mat3(unit_cell->basis);
        Ai = mat4x3_from_mat3(unit_cell->inv_basis);
        flags = unit_cell->flags;
    }

    const int PBC_ALL = (MD_UNIT_CELL_FLAG_PBC_X | MD_UNIT_CELL_FLAG_PBC_Y | MD_UNIT_CELL_FLAG_PBC_Z);
    if (flags & PBC_ALL != PBC_ALL) {
        // Unit cell either missing or not periodic along one or more axis
        vec4_t aabb_min = {0}, aabb_max = {0};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, in_x, in_y, in_z, NULL, NULL, num_points);

        // Construct A and Ai from aabb extent (use 0
    }

    // Precompute metric G = A^T A
    const mat4x3_t G = mat4x3_from_mat3(mat3_mul(mat3_transpose(unit_cell->basis), unit_cell->basis));

    // Choose grid resolution. Heuristic:  cell_dim ≈ |A|/CELL_EXT.
    // Using ~cutoff-ish spacing gives good pruning. Tweak CELL_EXT if you like.
    const float CELL_EXT = cutoff > 0.f ? cutoff : 1.0f;

    // Estimate cell_dim by measuring the extents of the box vectors (norms of columns of A)
    // This is only a heuristic for bin counts; the grid is still in fractional space.
    float ax = sqrtf(A.elem[0][0] * A.elem[0][0] + A.elem[1][0] * A.elem[1][0] + A.elem[2][0] * A.elem[2][0]);
    float by = sqrtf(A.elem[0][1] * A.elem[0][1] + A.elem[1][1] * A.elem[1][1] + A.elem[2][1] * A.elem[2][1]);
    float cz = sqrtf(A.elem[0][2] * A.elem[0][2] + A.elem[1][2] * A.elem[1][2] + A.elem[2][2] * A.elem[2][2]);

    uint32_t cell_dim[3] = {
        CLAMP((uint32_t)(ax / CELL_EXT + 0.5f), 1, 1024),
        CLAMP((uint32_t)(by / CELL_EXT + 0.5f), 1, 1024),
        CLAMP((uint32_t)(cz / CELL_EXT + 0.5f), 1, 1024),
    };

    const uint32_t c0 = cell_dim[0];
    const uint32_t c01 = cell_dim[0] * cell_dim[1];
    const size_t num_cells = (size_t)cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Temporary arrays
    uint32_t* cell_offset = (uint32_t*)md_vm_arena_push_zero(arena, (num_cells + 1) * sizeof(uint32_t));
    uint32_t* local_idx = (uint32_t*)md_vm_arena_push(arena, num_points * sizeof(uint32_t));
    uint32_t* cell_idx = (uint32_t*)md_vm_arena_push(arena, num_points * sizeof(uint32_t));
    elem_t* scratch_s = (elem_t*)md_vm_arena_push(arena, num_points * sizeof(elem_t));  // holds fractional coords before scatter
    elem_t* elements = (elem_t*)md_vm_arena_push(arena, num_points * sizeof(elem_t));   // sorted fractional coords

    md_timestamp_t t0 = md_time_current();

    const vec4_t  fcell_dim = vec4_set(cell_dim[0], cell_dim[1], cell_dim[2], 0);
    const ivec4_t icell_min = ivec4_set1(0);
    const ivec4_t icell_max = ivec4_set(cell_dim[0] - 1, cell_dim[1] - 1, cell_dim[2] - 1, 0);

    const ivec4_t c = ivec4_set(1, c0, c01, 0);

    // 1) Convert to fractional, wrap periodic axes into [0,1), bin to cells
    for (size_t i = 0; i < num_points; ++i) {
        vec4_t r = {in_x[i], in_y[i], in_z[i], 0};
        vec4_t s = mat4x3_mul_vec4(Ai, r);
        s = vec4_fract(s);

        vec4_t f = vec4_floor(vec4_mul(s, fcell_dim));
        ivec4_t ic = ivec4_from_vec4(f);

        ic = ivec4_clamp(ic, icell_min, icell_max);

        uint32_t ci = ivec4_dot(ic, c);
        ASSERT(ci < num_cells);

        local_idx[i] = cell_offset[ci]++;  // count for now
        cell_idx[i] = ci;

        // stash fractional coordinates
        MEMCPY(&scratch_s[i], &s, sizeof(elem_t));
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

    md_timestamp_t t2 = md_time_current();
    printf("2: %.2f\n", md_time_as_milliseconds(t2 - t1));

    // 3) Scatter fractional coords into 'elements' in cell order
    for (size_t i = 0; i < num_points; ++i) {
        uint32_t dst = cell_offset[cell_idx[i]] + local_idx[i];
        elements[dst] = scratch_s[i];
    }

    md_timestamp_t t3 = md_time_current();
    printf("3: %.2f\n", md_time_as_milliseconds(t3 - t2));

    /*
    // 4) Precompute wrapped neighbor indices (+ image shifts) for every cell
    typedef struct {
        uint32_t idx;
        int8_t kx, ky, kz;
    } nb_t;
    nb_t* nb13 = (nb_t*)md_vm_arena_push(arena, num_cells * 13 * sizeof(nb_t));

    for (uint32_t cz = 0; cz < cell_dim[2]; ++cz) {
        for (uint32_t cy = 0; cy < cell_dim[1]; ++cy) {
            for (uint32_t cx = 0; cx < cell_dim[0]; ++cx) {
                const uint32_t ci = (uint32_t)cz * c01 + (uint32_t)cy * c0 + (uint32_t)cx;
                nb_t* dst = nb13 + ci * 13;
                for (int n = 0; n < 13; ++n) {
                    int nx = (int)cx + FWD_NBRS[n][0];
                    int ny = (int)cy + FWD_NBRS[n][1];
                    int nz = (int)cz + FWD_NBRS[n][2];
                    int kx = 0, ky = 0, kz = 0;

                    nx = wrap_and_shift(nx, (int)cell_dim[0], periodic[0], &kx);
                    ny = wrap_and_shift(ny, (int)cell_dim[1], periodic[1], &ky);
                    nz = wrap_and_shift(nz, (int)cell_dim[2], periodic[2], &kz);

                    nb_t e;
                    if (nx < 0 || ny < 0 || nz < 0) {
                        // Non-periodic axis went OOB → mark as invalid
                        e.idx = (uint32_t)~0u;
                        e.kx = e.ky = e.kz = 0;
                    } else {
                        e.idx = (uint32_t)nz * c01 + (uint32_t)ny * c0 + (uint32_t)nx;
                        e.kx = (int8_t)kx;
                        e.ky = (int8_t)ky;
                        e.kz = (int8_t)kz;
                    }
                    dst[n] = e;
                }
            }
        }
    }
    */

    md_timestamp_t t4 = md_time_current();
    printf("4: %.2f\n", md_time_as_milliseconds(t4 - t3));

    // 5) Pair counting
    const float r2_cut = cutoff * cutoff;
    size_t result = 0;

// Macros for cell offsets/lengths in the sorted array
#define CELL_OFFSET(ci) (cell_offset[(ci)])
#define CELL_LENGTH(ci) (cell_offset[(ci) + 1] - cell_offset[(ci)])

    for (uint32_t cz = 0; cz < cell_dim[2]; ++cz) {
        for (uint32_t cy = 0; cy < cell_dim[1]; ++cy) {
            for (uint32_t cx = 0; cx < cell_dim[0]; ++cx) {
                const uint32_t ci = (uint32_t)cz * c01 + (uint32_t)cy * c0 + (uint32_t)cx;
                const uint32_t off_i = CELL_OFFSET(ci);
                const uint32_t len_i = CELL_LENGTH(ci);
                const elem_t* cell_i = elements + off_i;

                // Self cell: only j > i
                for (uint32_t a = 0; a < len_i; ++a) {
                    float sx = cell_i[a].x, sy = cell_i[a].y, sz = cell_i[a].z;
                    result += test_elem_frac(sx, sy, sz, cell_i + (a + 1), len_i - (a + 1), &G, periodic, r2_cut);
                }

                /*
                // 13 forward neighbors (already wrapped)
                const nb_t* nb = nb13 + ci * 13;
                for (int n = 0; n < 13; ++n) {
                    if (nb[n].idx == (uint32_t)~0u) continue;  // invalid (non-PBC OOB)
                    const uint32_t cj = nb[n].idx;
                    const uint32_t off_j = CELL_OFFSET(cj);
                    const uint32_t len_j = CELL_LENGTH(cj);
                    const elem_t* cell_j = elements + off_j;
                    const int kx = nb[n].kx, ky = nb[n].ky, kz = nb[n].kz;

                    for (uint32_t a = 0; a < len_i; ++a) {
                        float sx = cell_i[a].x, sy = cell_i[a].y, sz = cell_i[a].z;
                        result += test_elem_frac_shift(sx, sy, sz, kx, ky, kz, cell_j, len_j, &G, periodic, r2_cut);
                    }
                }
                */
            }
        }
    }

    md_timestamp_t t5 = md_time_current();
    printf("5: %.2f\n", md_time_as_milliseconds(t5 - t4));

#undef CELL_OFFSET
#undef CELL_LENGTH

    md_vm_arena_destroy(arena);
    return result;
}

static size_t do_brute_force(const float* in_x, const float* in_y, const float* in_z, size_t num_points, float cutoff, const md_unit_cell_t* unit_cell) {
    size_t count = 0;
    const float r2 = cutoff * cutoff;

    if (unit_cell) {
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
                if ((d.x * d.x + d.y * d.y + d.z * d.z) < r2) {
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

                if ((dx * dx + dy * dy + dz * dz) < r2) {
                    count += 1;
                }
            }
        }
    }

    return count;
}

#define CELL_EXT (6.0f)

static size_t do_pairwise_periodic(const float* in_x, const float* in_y, const float* in_z, size_t num_points, float cutoff, const md_unit_cell_t* unit_cell) {
    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));

    #define PBC_ALL (MD_UNIT_CELL_FLAG_PBC_X | MD_UNIT_CELL_FLAG_PBC_Y | MD_UNIT_CELL_FLAG_PBC_Z)

    mat4x3_t Ai = mat4x3_from_mat3(unit_cell->inv_basis);
    mat4x3_t A  = mat4x3_from_mat3(unit_cell->basis);
    uint32_t flags = unit_cell->flags;
    vec3_t ext = md_unit_cell_extent(unit_cell);

    if ((flags & PBC_ALL) != PBC_ALL) {
        // We need to calculate the AABB of the points
        float aabb_min[3] = {0};
        float aabb_max[3] = {0};
        md_util_aabb_compute(aabb_min, aabb_max, in_x, in_y, in_z, NULL, NULL, num_points);

        if (flags & MD_UNIT_CELL_FLAG_PBC_X) {
            ext.x = aabb_max[0];
            Ai.elem[0][0] = ext.x > 0.0f ? 1.0f / ext.x : 0.0f;
             A.elem[0][0] = ext.x;
        }
        if (flags & MD_UNIT_CELL_FLAG_PBC_Y) {
            ext.y = aabb_max[1];
            Ai.elem[1][1] = ext.y > 0.0f ? 1.0f / ext.y : 0.0f;
             A.elem[1][1] = ext.y;
        }
        if (flags & MD_UNIT_CELL_FLAG_PBC_Z) {
            ext.z = aabb_max[2];
            Ai.elem[2][2] = ext.z > 0.0f ? 1.0f / ext.z : 0.0f;
             A.elem[2][2] = ext.z;
        }
    }

    printf("unit cell ext: %f %f %f \n", ext.x, ext.y, ext.z);

    md_timestamp_t start = md_time_current();

    uint32_t cell_dim[3] = {
        CLAMP((uint32_t)(ext.x / CELL_EXT + 0.5f), 1, 1024),
        CLAMP((uint32_t)(ext.y / CELL_EXT + 0.5f), 1, 1024),
        CLAMP((uint32_t)(ext.z / CELL_EXT + 0.5f), 1, 1024),
    };

    vec4_t cell_ext = {
        ext.x / cell_dim[0],
        ext.y / cell_dim[1],
        ext.z / cell_dim[2],
        1.0f,
    };

    printf("cell ext: %f %f %f \n", cell_ext.x, cell_ext.y, cell_ext.z);

    uint32_t c0  = cell_dim[0];
    uint32_t c01 = cell_dim[0] * cell_dim[1];

    size_t num_cells = cell_dim[0] * cell_dim[1] * cell_dim[2];
    size_t cell_offset_count = num_cells + 1;

    printf("Cell dim: %i %i %i \n", cell_dim[0], cell_dim[1], cell_dim[2]);

    uint32_t* cell_offset = md_vm_arena_push_zero(temp_arena, cell_offset_count * sizeof(uint32_t));
    uint32_t* local_idx   = md_vm_arena_push(temp_arena, num_points * sizeof(uint32_t));
    uint32_t* cell_idx    = md_vm_arena_push(temp_arena, num_points * sizeof(uint32_t));
    elem_t* elements      = md_vm_arena_push(temp_arena, num_points * sizeof(elem_t));

    if ((flags & PBC_ALL) == PBC_ALL) {
        // Handle full PBC case
        for (size_t i = 0; i < num_points; ++i) {
            // Convert to fractional: s = Ai * r
            vec4_t r = {in_x[i], in_y[i], in_z[i], 0.0f};
            vec4_t s = mat4x3_mul_vec4(Ai, r);
            s = vec4_fract(vec4_fract(s));

            // Map to cell indices in [0, cell_dim[d]-1]
            uint32_t cx = (uint32_t)CLAMP((int)floorf(s.x * cell_dim[0]), 0, (int)cell_dim[0] - 1);
            uint32_t cy = (uint32_t)CLAMP((int)floorf(s.y * cell_dim[1]), 0, (int)cell_dim[1] - 1);
            uint32_t cz = (uint32_t)CLAMP((int)floorf(s.z * cell_dim[2]), 0, (int)cell_dim[2] - 1);

            uint32_t ci = cz * c01 + cy * c0 + cx;
            ASSERT(ci < num_cells);

            // bump counts
            uint32_t li  = cell_offset[ci]++;
            local_idx[i] = li;
            cell_idx[i]  = ci;
        }
    }

    uint32_t max_cell_count = 0;

    // Prefix sum the cell offsets
    uint32_t offset = 0;
    for (size_t i = 0; i < cell_offset_count; ++i) {
        uint32_t length = cell_offset[i];
        max_cell_count = MAX(max_cell_count, length);
        cell_offset[i] = offset;
        offset += length;
    }

    printf("Max cell count: %i\n", max_cell_count);

    // Calculate final destination index and write data
    for (size_t i = 0; i < num_points; ++i) {
        uint32_t ci = cell_idx[i];
        uint32_t dst_idx = cell_offset[ci] + local_idx[i];
        elements[dst_idx] = (elem_t) {
            in_x[i],
            in_y[i],
            in_z[i],
            (uint32_t)i,
        };
    }

    md_timestamp_t end = md_time_current();

    printf("Cell creation took: %f ms\n", md_time_as_milliseconds(end - start));

    data_t data = { 0 };
    const float r2 = cutoff * cutoff;
    data.r2 = md_mm256_set1_ps(cutoff * cutoff);
    data.elem = elements;
    data.cell_offset = cell_offset;

    size_t result = 0;

        // Convert 3D -> 1D cell index
    #define CELL_INDEX(x, y, z) ((z) * c01 + (y) * c0 + (x))
    #define CELL_OFFSET(ci) (cell_offset[ci])
    #define CELL_LENGTH(ci) (cell_offset[ci + 1] - cell_offset[ci])
    #define TEST_ELEM(x, y, z, r2, elem, len) \
        test_elem(x, y, z, r2, elem, len)
    #define TEST_CELL(x, y, z, r2, ci) \
        TEST_ELEM(x, y, z, r2, elements + CELL_OFFSET(ci), CELL_LENGTH(ci))

    // Forward neighbor offsets in 3D grid (no PBC)
    static const int8_t FWD_NBRS[13][3] = {
        {1, 0, 0},
        {-1, 1, 0},
        {0, 1, 0},
        {1, 1, 0},
        {-1, -1, 1},
        {0, -1, 1},
        {1, -1, 1},
        {-1, 0, 1},
        {0, 0, 1},
        {1, 0, 1},
        {-1, 1, 1},
        {0, 1, 1},
        {1, 1, 1},
    };

    #if 0
    for (uint32_t ci = 0; ci < num_cells - 1; ++ci) {
        uint32_t offset = cell_offset[ci];
        uint32_t length = cell_offset[ci + 1] - offset;
        const elem_t* elem = elements + offset;
        for (size_t i = 0; i < length; ++i) {
            float x = elem[i].x;
            float y = elem[i].y;
            float z = elem[i].z;
            result += TEST_ELEM(x, y, z, r2, elem + (i + 1), length - (i + 1));

            for (uint32_t cj = ci + 1; cj < num_cells; ++cj) {
                result += TEST_CELL(x, y, z, r2, cj);
            }
        }
    }

    #else

    for (uint32_t cz = 0; cz < cell_dim[2]; ++cz) {
        for (uint32_t cy = 0; cy < cell_dim[1]; ++cy) {
            for (uint32_t cx = 0; cx < cell_dim[0]; ++cx) {
                uint32_t ci = CELL_INDEX(cx, cy, cz);
                uint32_t offset_i = cell_offset[ci];
                uint32_t length_i = cell_offset[ci + 1] - offset_i;
                const elem_t* cell_i = elements + offset_i;

                // Self-cell: only j > i
                for (uint32_t a = 0; a < length_i; ++a) {
                    float x = cell_i[a].x;
                    float y = cell_i[a].y;
                    float z = cell_i[a].z;
                    result += TEST_ELEM(x, y, z, r2, cell_i + a + 1, length_i - (a + 1));
                }

                // Forward neighbors
                for (int n = 0; n < 13; ++n) {
                    int nx = (int)cx + FWD_NBRS[n][0];
                    int ny = (int)cy + FWD_NBRS[n][1];
                    int nz = (int)cz + FWD_NBRS[n][2];

                    if (nx < 0 || ny < 0 || nz < 0 || nx >= (int)cell_dim[0] || ny >= (int)cell_dim[1] || nz >= (int)cell_dim[2]) {
                        continue;  // Skip out-of-bounds
                    }

                    uint32_t cj = CELL_INDEX(nx, ny, nz);
                    uint32_t offset_j = cell_offset[cj];
                    uint32_t length_j = cell_offset[cj + 1] - offset_j;
                    const elem_t* cell_j = elements + offset_j;

                    for (uint32_t a = 0; a < length_i; ++a) {
                        float x = cell_i[a].x;
                        float y = cell_i[a].y;
                        float z = cell_i[a].z;
                        result += TEST_ELEM(x, y, z, r2, cell_j, length_j);
                    }
                }
            }
        }
    }
    #endif

    #undef CELL_INDEX
    #undef CELL_OFFSET
    #undef CELL_LENGTH
    #undef TEST_ELEM
    #undef TEST_CELL

    md_vm_arena_destroy(temp_arena);

    return result;
}

UTEST(spatial_hash, n2) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro_data, alloc));

    const vec4_t mask = md_unit_cell_pbc_mask(&mol.unit_cell);
    mat4x3_t I = mat4x3_from_mat3(mol.unit_cell.inv_basis);
    mat4x3_t M = mat4x3_from_mat3(mol.unit_cell.basis);

    for (size_t i = 0; i < mol.atom.count; ++i) {
        vec4_t original = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
        vec4_t coord = mat4x3_mul_vec4(I, original);
        coord = vec4_fract(coord);
        coord = mat4x3_mul_vec4(M, coord);
        coord = vec4_blend(original, coord, mask);

        mol.atom.x[i] = coord.x;
        mol.atom.y[i] = coord.y;
        mol.atom.z[i] = coord.z;
    }

    md_unit_cell_t unit_cell = mol.unit_cell;
    // Clear pbc flags

    size_t expected_count = 3701955;
    if (true) {
        unit_cell.flags &= ~(MD_UNIT_CELL_FLAG_PBC_X | MD_UNIT_CELL_FLAG_PBC_Y | MD_UNIT_CELL_FLAG_PBC_Z);
        expected_count = 3701955;
    }

    // Custom implementation of pairwise periodic N^2
    md_timestamp_t start = md_time_current();
    //size_t count = do_pairwise_periodic(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell);
    size_t count = do_pairwise_periodic_triclinic(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell);
    md_timestamp_t end = md_time_current();
    printf("Custom: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, count);
    if (count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, count);
    }
    
    // Current implementation of spatial hash for N^2
    start = md_time_current();
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &unit_cell, alloc);
    iter_excl_data_t data = {0};
    for (uint32_t i = 0; i < mol.atom.count; ++i) {
        vec3_t pos = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]};
        data.exclude_idx = i;
        md_spatial_hash_query(spatial_hash, pos, 5.0f, iter_excl_fn, &data);
    }
    end = md_time_current();
    printf("Spatial hash: %f ms\n", md_time_as_milliseconds(end - start));
    size_t sh_count = data.count;
    ASSERT_TRUE(spatial_hash);
    EXPECT_EQ(expected_count, count);
    if (count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, count);
    }

    // Brute force
    start = md_time_current();
    size_t bf_count = do_brute_force(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell);
    end = md_time_current();
    printf("Brute force: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, bf_count);
    if (count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, count);
    }

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