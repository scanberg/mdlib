#include "md_spatial_hash.h"

#include "md_allocator.h"
#include "md_log.h"
#include "md_array.h"
#include "md_simd.h"
#include "md_util.h"

#include <float.h>
#include <string.h>
#include <stdint.h>

#define MD_SPATIAL_HASH_MAGIC 0xF67cbab6
#define CELL_EXT (6.0f)
#define INV_CELL_EXT (1.0f / CELL_EXT)
#define SQRT_CELL_EXT (8.4852813742385702f)

static inline vec3_t unpack_unorm_vec3(uint32_t p) {
    float x = ((p >>  0) & 0x3FF) / 1023.0f;
    float y = ((p >> 10) & 0x3FF) / 1023.0f;
    float z = ((p >> 20) & 0x3FF) / 1023.0f;
    return (vec3_t) {x,y,z};
}

static inline uint32_t pack_unorm_vec3(vec3_t v) {
    uint32_t x = (uint32_t)(v.x * 1023.0f + 0.5f);
    uint32_t y = (uint32_t)(v.y * 1023.0f + 0.5f);
    uint32_t z = (uint32_t)(v.z * 1023.0f + 0.5f);
    return (z << 20) | (y << 10) | x;
}

bool init(md_spatial_hash_t* hash, const float* in_x, const float* in_y, const float* in_z, size_t count, size_t stride, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(alloc);

    stride = MAX(stride, sizeof(float));

    md_allocator_i* temp_alloc = default_temp_allocator;
    uint32_t* cell_coord        = md_alloc(temp_alloc, sizeof(uint32_t) * count);
    uint32_t* packed_cell_coord = md_alloc(temp_alloc, sizeof(uint32_t) * count);
    uint32_t* local_idx         = md_alloc(temp_alloc, sizeof(uint32_t) * count);

    const vec4_t ext = vec4_from_vec3(pbc_ext, 0);
    const vec4_t ref = vec4_mul_f(ext, 0.5f);
    vec4_t aabb_min = {0};
    vec4_t aabb_max = {0};
    md_util_compute_aabb_periodic_xyz((vec3_t*)(&aabb_min), (vec3_t*)(&aabb_max), in_x, in_y, in_z, count, stride, pbc_ext);

    const int32_t cell_min[3] = {
        (int32_t)floorf(aabb_min.x * INV_CELL_EXT),
        (int32_t)floorf(aabb_min.y * INV_CELL_EXT),
        (int32_t)floorf(aabb_min.z * INV_CELL_EXT),
    };

    const int32_t cell_max[3] = {
        (int32_t)floorf(aabb_max.x * INV_CELL_EXT) + 1,
        (int32_t)floorf(aabb_max.y * INV_CELL_EXT) + 1,
        (int32_t)floorf(aabb_max.z * INV_CELL_EXT) + 1,
    };

    const int32_t cell_dim[3] = {
        MAX(1, cell_max[0] - cell_min[0]),
        MAX(1, cell_max[1] - cell_min[1]),
        MAX(1, cell_max[2] - cell_min[2]),
    };

    const uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    {
        size_t i = 0;
        const size_t simd_count = (count / md_simd_widthf) * md_simd_widthf;
        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef x = md_simd_loadf((const float*)((const char*)in_x + i * stride));
            md_simd_typef y = md_simd_loadf((const float*)((const char*)in_y + i * stride));
            md_simd_typef z = md_simd_loadf((const float*)((const char*)in_z + i * stride));

            x = md_simd_mulf(md_simd_deperiodizef(x, md_simd_set1f(ref.x), md_simd_set1f(ext.x)), md_simd_set1f(INV_CELL_EXT));
            y = md_simd_mulf(md_simd_deperiodizef(y, md_simd_set1f(ref.y), md_simd_set1f(ext.y)), md_simd_set1f(INV_CELL_EXT));
            z = md_simd_mulf(md_simd_deperiodizef(z, md_simd_set1f(ref.z), md_simd_set1f(ext.z)), md_simd_set1f(INV_CELL_EXT));

            md_simd_typef fx = md_simd_floorf(x);
            md_simd_typef fy = md_simd_floorf(y);
            md_simd_typef fz = md_simd_floorf(z);

            md_simd_typei ccx = md_simd_subi(md_simd_convert_f_to_i(fx), md_simd_set1i(cell_min[0]));
            md_simd_typei ccy = md_simd_subi(md_simd_convert_f_to_i(fy), md_simd_set1i(cell_min[1]));
            md_simd_typei ccz = md_simd_subi(md_simd_convert_f_to_i(fz), md_simd_set1i(cell_min[2]));

            md_simd_typei pcx = md_simd_convert_f_to_i(md_simd_fmaddf(md_simd_subf(x, fx), md_simd_set1f(1023.0f), md_simd_set1f(0.5f)));
            md_simd_typei pcy = md_simd_convert_f_to_i(md_simd_fmaddf(md_simd_subf(y, fy), md_simd_set1f(1023.0f), md_simd_set1f(0.5f)));
            md_simd_typei pcz = md_simd_convert_f_to_i(md_simd_fmaddf(md_simd_subf(z, fz), md_simd_set1f(1023.0f), md_simd_set1f(0.5f)));

            md_simd_typei cc = md_simd_ori(md_simd_shift_lefti(ccz, 20), md_simd_ori(md_simd_shift_lefti(ccy, 10), ccx));
            md_simd_typei pc = md_simd_ori(md_simd_shift_lefti(pcz, 20), md_simd_ori(md_simd_shift_lefti(pcy, 10), pcx));

            md_simd_storei((int*)cell_coord + i, cc);
            md_simd_storei((int*)packed_cell_coord + i, pc);
        }

        // Handle remainder
        for (; i < count; ++i) {
            vec4_t c = {
                *(const float*)((const char*)in_x + i * stride),
                *(const float*)((const char*)in_y + i * stride),
                *(const float*)((const char*)in_z + i * stride),
                0
            };
            c = vec4_mul_f(vec4_deperiodize(c, ref, ext), INV_CELL_EXT);

            vec4_t f = vec4_floor(c);
            vec4_t p = vec4_fmadd(vec4_sub(c, f), vec4_set1(1023.0f), vec4_set1(0.5f));
            uint32_t cc =
                (((int32_t)f.z - cell_min[2]) << 20) |
                (((int32_t)f.y - cell_min[1]) << 10) |
                 ((int32_t)f.x - cell_min[0]);

            cell_coord[i]        = cc;
            packed_cell_coord[i] = ((uint32_t)p.z << 20) | ((uint32_t)p.y << 10) | (uint32_t)p.x;
        }
    }

    // Allocate needed data
    uint32_t* cells         = md_alloc(alloc, sizeof(uint32_t) * cell_count);
    uint32_t* coords        = md_alloc(alloc, sizeof(uint32_t) * count);
    uint32_t* original_idx  = md_alloc(alloc, sizeof(uint32_t) * count);

    MEMSET(cells, 0, sizeof(uint32_t) * cell_count);

    const uint32_t cell_dim_01 = cell_dim[0] * cell_dim[1];

    for (size_t i = 0; i < count; ++i) {
        uint32_t cc = cell_coord[i];
        uint32_t cz = ((cc >> 20) & 0x3FF);
        uint32_t cy = ((cc >> 10) & 0x3FF);
        uint32_t cx = ((cc >>  0) & 0x3FF);
        uint32_t idx = cz * cell_dim_01 + cy * cell_dim[0] + cx;
        ASSERT(idx < cell_count);

        local_idx[i] = cells[idx]++;
        ASSERT(cells[idx] < 1024 && "Too many entities per cell");
    }

    uint32_t offset = cells[0];
    for (uint32_t i = 1; i < cell_count; ++i) {
        uint32_t length = cells[i];
        cells[i] = (offset << 10) | length;
        offset += length;
    }

#if DEBUG
    uint32_t cell = cells[cell_count-1];
    ASSERT((cell >> 10) + (cell & 1023) == count);
#endif

    // Write data to correct location
    for (size_t i = 0; i < count; ++i) {
        uint32_t cc = cell_coord[i];
        uint32_t cz = ((cc >> 20) & 0x3FF);
        uint32_t cy = ((cc >> 10) & 0x3FF);
        uint32_t cx = ((cc >>  0) & 0x3FF);
        uint32_t cell_idx = cz * cell_dim_01 + cy * cell_dim[0] + cx;
        uint32_t dst_idx = (cells[cell_idx] >> 10) + local_idx[i];
        coords[dst_idx]  = packed_cell_coord[i];
        original_idx[dst_idx] = (uint32_t)i;
    }

    md_free(temp_alloc, local_idx,          sizeof(uint32_t) * count);
    md_free(temp_alloc, packed_cell_coord,  sizeof(uint32_t) * count);
    md_free(temp_alloc, cell_coord,         sizeof(uint32_t) * count);

    hash->cells    = cells;
    hash->coords   = coords;
    hash->indices  = original_idx;
    hash->alloc    = alloc;
    MEMCPY(hash->cell_min, cell_min, sizeof(cell_min));
    hash->coord_count = (uint32_t)count;
    MEMCPY(hash->cell_dim, cell_dim, sizeof(cell_dim));
    hash->magic = MD_SPATIAL_HASH_MAGIC;
    hash->pbc_ext = vec4_from_vec3(pbc_ext, 0);

    return true;
}

bool md_spatial_hash_init(md_spatial_hash_t* hash, const vec3_t* pos, int64_t count, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(pos);
    ASSERT(alloc);

    if (count < 0) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid count");
        return false;
    }

    return init(hash, &pos->x, &pos->y, &pos->z, (size_t)count, sizeof(vec3_t), pbc_ext, alloc);
}

bool md_spatial_hash_init_soa(md_spatial_hash_t* hash, const float* x, const float* y, const float* z, int64_t count, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(alloc);

    if (count < 0) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid count");
        return false;
    }

    return init(hash, x, y, z, (size_t)count, 0, pbc_ext, alloc);
}

bool md_spatial_hash_free(md_spatial_hash_t* hash) {
    ASSERT(hash);
    ASSERT(hash->magic == MD_SPATIAL_HASH_MAGIC);
    ASSERT(hash->alloc);
    uint32_t cell_count = hash->cell_dim[0] * hash->cell_dim[1] * hash->cell_dim[2];
    uint32_t coord_count = hash->coord_count;
    md_free(hash->alloc, hash->cells,   sizeof(uint32_t) * cell_count);
    md_free(hash->alloc, hash->coords,  sizeof(uint32_t) * coord_count);
    md_free(hash->alloc, hash->indices, sizeof(uint32_t) * coord_count);
    return true;
}

// Non periodic version (simpler)
static inline bool query_pos_rad(const md_spatial_hash_t* hash, vec3_t position, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    vec4_t pos = vec4_from_vec3(position, 0);
    float rad = radius;
    float rad2 = radius * radius;
    //float cell_rad2 = (radius + SQRT_CELL_EXT) * (radius + SQRT_CELL_EXT);

    const int32_t* cell_min = hash->cell_min;
    const int32_t* cell_dim = hash->cell_dim;

    int32_t cell_beg[3] = {
        CLAMP((int32_t)(floorf((pos.x - rad) * INV_CELL_EXT)) - cell_min[0], 0, cell_dim[0] - 1),
        CLAMP((int32_t)(floorf((pos.y - rad) * INV_CELL_EXT)) - cell_min[1], 0, cell_dim[1] - 1),
        CLAMP((int32_t)(floorf((pos.z - rad) * INV_CELL_EXT)) - cell_min[2], 0, cell_dim[2] - 1)
    };

    int32_t cell_end[3] = {
        CLAMP((int32_t)(ceilf((pos.x + rad) * INV_CELL_EXT)) - cell_min[0], 0, cell_dim[0]),
        CLAMP((int32_t)(ceilf((pos.y + rad) * INV_CELL_EXT)) - cell_min[1], 0, cell_dim[1]),
        CLAMP((int32_t)(ceilf((pos.z + rad) * INV_CELL_EXT)) - cell_min[2], 0, cell_dim[2])
    };

    int32_t cc[3];
    vec4_t cc_min = {0};
    for (cc[2] = cell_beg[2]; cc[2] < cell_end[2]; ++cc[2]) {
        cc_min.z = (cell_min[2] + cc[2]) * CELL_EXT;
        for (cc[1] = cell_beg[1]; cc[1] < cell_end[1]; ++cc[1]) {
            cc_min.y = (cell_min[1] + cc[1]) * CELL_EXT;
            for (cc[0] = cell_beg[0]; cc[0] < cell_end[0]; ++cc[0]) {
                cc_min.x = (cell_min[0] + cc[0]) * CELL_EXT;
                //vec3_t cc_mid = vec3_add_f(cc_min, 0.5f * CELL_EXT);
                // Check the distance from cell centrum to the point, for early discard of entire cell
                //if (vec3_distance_squared(cc_mid, xyz) < cell_rad2) {
                    uint32_t cell_idx = cc[2] * cell_dim[1] * cell_dim[0] + cc[1] * cell_dim[0] + cc[0];
                    uint32_t cell_data = hash->cells[cell_idx];
                    uint32_t cell_offset = cell_data >> 10;
                    uint32_t cell_length = cell_data & 1023;
                    for (uint32_t i = cell_offset; i < cell_offset + cell_length; ++i) {
                        vec4_t p = vec4_add(cc_min, vec4_mul_f(vec4_from_vec3(unpack_unorm_vec3(hash->coords[i]), 0), CELL_EXT));
                        const float d2 = vec4_distance_squared(p, pos);
                        if (d2 < rad2) {
                            if (!iter(hash->indices[i], vec3_from_vec4(p), user_param)) {
                                return false;
                            }
                        }
                    }
                //}
            }
        }
    }
    return true;
}

static inline void inc_cell(int* cell_idx, int cell_beg, int cell_dim, int cell_pbc) {
    *cell_idx += 1;
    int cc = cell_beg + *cell_idx;
    if (cc >= cell_dim) {
        *cell_idx += cell_pbc - cell_dim;
    }
}

// Periodic version (more advanced)
// The periodic domain is aligned to the cells at the origin, but not at the other end since we use a fixed cell-size
// Furthermore, the cell domain usually only spans a subrange of the periodic domain
// This makes the query non-trivial as in the non periodic case.

static inline bool query_pos_rad_periodic(const md_spatial_hash_t* hash, vec3_t position, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    const vec4_t ref = vec4_mul_f(hash->pbc_ext, 0.5f);
    const vec4_t pbc_ext = hash->pbc_ext;
    vec4_t pos = vec4_from_vec3(position, 0);
    float rad = radius;
    float rad2 = radius * radius;
    //float cell_rad2 = (radius + SQRT_CELL_EXT) * (radius + SQRT_CELL_EXT);

    const int32_t cell_min[3] = { hash->cell_min[0], hash->cell_min[1], hash->cell_min[2] };
    const int32_t cell_dim[3] = { hash->cell_dim[0], hash->cell_dim[1], hash->cell_dim[2] };

    int32_t cell_pbc[3] = {
        (int32_t)(pbc_ext.x * INV_CELL_EXT) + 1,
        (int32_t)(pbc_ext.y * INV_CELL_EXT) + 1,
        (int32_t)(pbc_ext.z * INV_CELL_EXT) + 1,
    };

    // Deperiodize the extent of the search (pos +/- rad) with respect to the periodic domain
    vec4_t cell_pos_rad_min = vec4_floor(vec4_mul_f(vec4_deperiodize(vec4_sub_f(pos, rad), ref, pbc_ext), INV_CELL_EXT));

    // cell_pos_min is within the correct period, now we can skip forward to the first cell which actually contains data
    // We use the convention of [0,dim[ rather than [min,min+dim[,
    int32_t cell_beg[3] = {
        MAX((int32_t)(cell_pos_rad_min.x) - cell_min[0], 0),
        MAX((int32_t)(cell_pos_rad_min.y) - cell_min[1], 0),
        MAX((int32_t)(cell_pos_rad_min.z) - cell_min[2], 0),
    };

    // Extent in cells to search
    int32_t cell_ext[3] = {
        MIN((int32_t)(ceilf((pos.x + rad) * INV_CELL_EXT) - floorf((pos.x - rad) * INV_CELL_EXT)) + 1, cell_pbc[0] - 1) + 1,
        MIN((int32_t)(ceilf((pos.y + rad) * INV_CELL_EXT) - floorf((pos.y - rad) * INV_CELL_EXT)) + 1, cell_pbc[1] - 1) + 1,
        MIN((int32_t)(ceilf((pos.z + rad) * INV_CELL_EXT) - floorf((pos.z - rad) * INV_CELL_EXT)) + 1, cell_pbc[2] - 1) + 1,
    };

    // If cell beg is outside of occupied cell domain [cell_min, cell_dim[ we skip forward to next period and accomodate for that jump
    for (int i = 0; i < 3; ++i) {
        if (cell_beg[i] >= cell_dim[i]) {
            // How far ahead we need skip in order to end up within the occupied cell domain of the next period
            int delta = cell_pbc[i] - cell_beg[i];
            cell_beg[i] = (cell_beg[i] + delta) % cell_pbc[i];
            cell_ext[i] -= delta;
            if (cell_ext[i] < 0) return true;
        }
    }

    int32_t ci[3] = {0};
    int32_t cc[3] = {0};
    vec4_t cc_min = {0};
    for (ci[2] = 0; ci[2] < cell_ext[2]; inc_cell(&ci[2], cell_beg[2], cell_dim[2], cell_pbc[2])) {
        cc[2] = (cell_beg[2] + ci[2]) % cell_pbc[2];
        ASSERT(cc[2] < cell_dim[2]);
        cc_min.z = (cell_min[2] + cc[2]) * CELL_EXT;
        for (ci[1] = 0; ci[1] < cell_ext[1]; inc_cell(&ci[1], cell_beg[1], cell_dim[1], cell_pbc[1])) {
            cc[1] = (cell_beg[1] + ci[1]) % cell_pbc[1];
            ASSERT(cc[1] < cell_dim[1]);
            cc_min.y = (cell_min[1] + cc[1]) * CELL_EXT;
            for (ci[0] = 0; ci[0] < cell_ext[0]; inc_cell(&ci[0], cell_beg[0], cell_dim[0], cell_pbc[0])) {
                cc[0] = (cell_beg[0] + ci[0]) % cell_pbc[0];
                ASSERT(cc[0] < cell_dim[0]);
                cc_min.x = (cell_min[0] + cc[0]) * CELL_EXT;
                //vec3_t cc_mid = vec3_add_f(cc_min, 0.5f * CELL_EXT);
                // Check the distance from cell centrum to the point, for early discard of entire cell
                //if (vec3_distance_squared(cc_mid, xyz) < cell_rad2) {
                uint32_t cell_idx = cc[2] * cell_dim[1] * cell_dim[0] + cc[1] * cell_dim[0] + cc[0];
                uint32_t cell_data = hash->cells[cell_idx];
                uint32_t cell_offset = cell_data >> 10;
                uint32_t cell_length = cell_data & 1023;
                for (uint32_t i = cell_offset; i < cell_offset + cell_length; ++i) {
                    vec4_t p = vec4_add(cc_min, vec4_mul_f(vec4_from_vec3(unpack_unorm_vec3(hash->coords[i]), 0), CELL_EXT));
                    const float d2 = vec4_periodic_distance_squared(p, pos, pbc_ext);
                    if (d2 < rad2) {
                        if (!iter(hash->indices[i], vec3_from_vec4(p), user_param)) {
                            return false;
                        }
                    }
                }
                //}
            }
        }
    }
    return true;
}

bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(spatial_hash);
    ASSERT(iter);

    // This is for non periodic lookup
    if (vec4_equal(spatial_hash->pbc_ext, vec4_zero())) {
        return query_pos_rad(spatial_hash, pos, radius, iter, user_param);
    } else {
        return query_pos_rad_periodic(spatial_hash, pos, radius, iter, user_param);
    }
}
