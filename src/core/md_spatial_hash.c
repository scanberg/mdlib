#include "md_spatial_hash.h"

#include "md_allocator.h"
#include "md_log.h"
#include "md_array.inl"
#include "md_simd.h"

#include <float.h>
#include <string.h>

#define MD_SPATIAL_HASH_MAGIC 0xF67cbab6

typedef struct md_spatial_hash_cell_t {
    uint32_t offset;
    uint32_t length;
} md_spatial_hash_cell_t;

/*
// Store as relative cell-coordinates quantized to 16 bits: Über precision
typedef union md_spatial_hash_coord_t {
    struct {
        uint16_t x;
        uint16_t y;
        uint16_t z;
        uint16_t w;
    };
    uint16_t elem[4];
} md_spatial_hash_coord_t;

static inline vec4_t decompress_coord(const uint16_t compressed_coord[4], vec4_t cell_min, float cell_ext) {
    const float scl = cell_ext / (1 << 16);
    vec4_t coord = {.mm128 = md_simd_set_f128(compressed_coord[0], compressed_coord[1], compressed_coord[2], compressed_coord[3])};
    return vec4_add(cell_min, vec4_mul_f(coord, scl));
}

static inline md_spatial_hash_coord_t compress_coord(vec4_t coord, vec4_t cell_min, float cell_ext) {
    const float scl = (1 << 16) / cell_ext;
    coord = vec4_sub(coord, cell_min);
    coord = vec4_mul_f(coord, scl);
    coord = vec4_add_f(coord, 0.5f);
    return (md_spatial_hash_coord_t){(uint16_t)coord.x, (uint16_t)coord.y, (uint16_t)coord.z, (uint16_t)coord.w};
}
*/

static inline float get_coord(const float* in_coord, int64_t stride, int64_t idx) {
    return *(const float*)((const char*)in_coord + stride * idx);
}

static inline vec4_t get_vec4(const float* in_coords[3], int64_t stride, int64_t idx) {
    return (vec4_t) {
        *(const float*)((const char*)in_coords[0] + stride * idx),
        in_coords[1] ? *(const float*)((const char*)in_coords[1] + stride * idx) : 0,
        in_coords[2] ? *(const float*)((const char*)in_coords[2] + stride * idx) : 0,
        0,
    };
}

bool md_spatial_hash_init(md_spatial_hash_t* hash, const md_spatial_hash_args_t* args) {
    ASSERT(hash);
    ASSERT(args);
    ASSERT(args->alloc);

    md_allocator_i* alloc = args->alloc;
    md_allocator_i* temp_alloc = args->temp_alloc ? args->temp_alloc : default_temp_allocator;

    const float* in_coords[3] = {
        args->coords.x,
        args->coords.y,
        args->coords.z,
    };
    int64_t stride = args->coords.stride;
    int64_t count = args->coords.count;
    int dim = 0;
    float cell_ext = args->cell_ext;

    if (in_coords[0] != NULL) {
        dim = 1;
        if (in_coords[1] != NULL) {
            dim = 2;
            if (in_coords[1] != NULL) {
                dim = 3;
            }
        }
    }

    if (count == 0) {
        // This is valid, but no data is given
        memset(hash, 0, sizeof(md_spatial_hash_t));
        hash->magic = MD_SPATIAL_HASH_MAGIC;
        return true;
    }

    if (count < 0) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid count");
        return false;
    }

    if (dim == 0) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid input coordinate dimensions");
        return false;
    }

    if (stride != 0) {
        // Check the stride
        if (stride < 0 || stride % 4 != 0) {
            md_printf(MD_LOG_TYPE_ERROR, "Invalid stride");
            return false;
        }
    } else {
        // if stride is zero, we set it to assume packed floats, i.e. 4 bytes
        stride = 4;   
    }

    if (cell_ext < 0) {
        md_printf(MD_LOG_TYPE_ERROR, "Cell Ext parameter (%g) must be set to zero or positive.", cell_ext);
        return false;
    }

    if (hash->coords != 0 || hash->cells != 0) {
        md_print(MD_LOG_TYPE_ERROR, "Spatial Hash structure is not cleared, potentially leaking memory here.");
    }

    vec4_t aabb_min = vec4_from_float( FLT_MAX);
    vec4_t aabb_max = vec4_from_float(-FLT_MAX);

    for (int64_t i = 0; i < count; ++i) {
        vec4_t c = get_vec4(in_coords, stride, i);
        aabb_min = vec4_min(aabb_min, c);
        aabb_max = vec4_max(aabb_max, c);
    }

    vec4_t aabb_mid = vec4_lerp(aabb_min, aabb_max, 0.5f);
    vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);

    if (cell_ext == 0) {
        // Try to pick a good cell_ext
        // Assume that we roughly have a uniform distribution of points in the space
        // We would ideally have less than ~32 points per cell, this is a crude guestimate
        // There is a balance between the overhead of traversing cells compared to the cost of actually checking the points.
        // So this is something which should be measured.

        // Compute the approximate area for a target number of points
        float spatial_unit = 1;
        for (int i = 0; i < dim; ++i) {
            spatial_unit *= (aabb_ext.elem[i] > 0 ? aabb_ext.elem[i] : 1);
        }    

        const float target = 32;
        const float spatial_unit_per_point = spatial_unit / (float)count;
        const float approx_spatial_units = target * spatial_unit_per_point;

        switch (dim) {
        case 1:
            cell_ext = MAX(approx_spatial_units, aabb_ext.x / 1024.0f);
            break;
        case 2:
            cell_ext = MAX(sqrtf(approx_spatial_units), MAX(aabb_ext.x, aabb_ext.y) / 1024.0f);
            break;
        case 3:
            cell_ext = MAX(cbrtf(approx_spatial_units), MAX(aabb_ext.x, MAX(aabb_ext.y, aabb_ext.z)) / 1024.0f);
            break;
        default:
            ASSERT(false);
        }
    }

    uint32_t cell_dim[3] = {
        (uint32_t)(aabb_ext.x / cell_ext) + 1,
        (uint32_t)(aabb_ext.y / cell_ext) + 1,
        (uint32_t)(aabb_ext.z / cell_ext) + 1,
    };

    ASSERT(cell_dim[0] > 0);
    ASSERT(cell_dim[1] > 0);
    ASSERT(cell_dim[2] > 0);

    vec4_t   full_ext = {cell_dim[0] * cell_ext, cell_dim[1] * cell_ext, cell_dim[2] * cell_ext, 0};
    vec4_t   cell_min = vec4_sub(aabb_mid, vec4_mul_f(full_ext, 0.5f));
    uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    ASSERT(full_ext.x >= aabb_ext.x && full_ext.y >= aabb_ext.y && full_ext.z >= aabb_ext.z);
    ASSERT(cell_min.x <= aabb_min.x && cell_min.y <= aabb_min.y && cell_min.z <= aabb_min.z);

    // Allocate needed data
    md_spatial_hash_cell_t* cells = 0;
    float* coords = 0;
    uint32_t* original_idx = 0;
    md_array_resize(cells, cell_count, alloc);
    md_array_resize(coords, count * dim, alloc);
    md_array_resize(original_idx, count, alloc);

    memset(cells, 0, cell_count * sizeof(md_spatial_hash_cell_t));

    uint32_t* cell_idx  = 0;
    uint32_t* local_idx = 0;
    md_array_resize(cell_idx,  count, temp_alloc);
    md_array_resize(local_idx, count, temp_alloc);

    const float inv_cell_ext = 1.0f / cell_ext;

    switch (dim) {
    case 1:
        for (uint32_t i = 0; i < (uint32_t)count; ++i) {
            uint32_t idx = (uint32_t)((get_coord(in_coords[0], stride, i) - cell_min.x) * inv_cell_ext);
            cell_idx[i]  = idx;
            local_idx[i] = cells[idx].length++;
            ASSERT(cells[idx].length != 0 && "hash cell has wrapped, too many entities per cell");
        }
        break;
    case 2:
        for (uint32_t i = 0; i < (uint32_t)count; ++i) {
            uint32_t cell_coord[2] = {
                (uint32_t)((get_coord(in_coords[0], stride, i) - cell_min.x) * inv_cell_ext),
                (uint32_t)((get_coord(in_coords[1], stride, i) - cell_min.y) * inv_cell_ext),
            };

            uint32_t idx = cell_coord[1] * cell_dim[0] + cell_coord[0];
            cell_idx[i]  = idx;
            local_idx[i] = cells[idx].length++;
            ASSERT(cells[idx].length != 0 && "hash cell has wrapped, too many entities per cell");
        }
        break;
    case 3:
        for (uint32_t i = 0; i < (uint32_t)count; ++i) {
            uint32_t cell_coord[3] = {
                (uint32_t)((get_coord(in_coords[0], stride, i) - cell_min.x) * inv_cell_ext),
                (uint32_t)((get_coord(in_coords[1], stride, i) - cell_min.y) * inv_cell_ext),
                (uint32_t)((get_coord(in_coords[2], stride, i) - cell_min.z) * inv_cell_ext),
            };

            uint32_t idx = cell_coord[2] * cell_dim[0] * cell_dim[1] + cell_coord[1] * cell_dim[0] + cell_coord[0];
            cell_idx[i]  = idx;
            local_idx[i] = cells[idx].length++;
            ASSERT(cells[idx].length != 0 && "hash cell has wrapped, too many entities per cell");
        }
        break;
    default:
        ASSERT(false);
    }

    for (uint32_t i = 1; i < cell_count; ++i) {
        cells[i].offset = cells[i-1].offset + cells[i-1].length;
    }

#if DEBUG
    md_spatial_hash_cell_t* cell = md_array_last(cells);
    ASSERT(cell);
    ASSERT(cell->offset + cell->length == count);
#endif

    for (int64_t i = 0; i < count; ++i) {
        uint32_t dst_idx = cells[cell_idx[i]].offset + local_idx[i];
        ASSERT(dst_idx < count);

        for (int j = 0; j < dim; ++j) {
            coords[dst_idx * dim + j] = get_coord(in_coords[j], stride, i);
        }
        original_idx[dst_idx] = (uint32_t)i;
    }

    md_array_free(local_idx, temp_alloc);
    md_array_free(cell_idx,  temp_alloc);

    hash->cells    = cells;
    hash->coords   = coords;
    hash->indices  = original_idx;
    hash->alloc    = alloc;
    hash->cell_min = vec3_from_vec4(cell_min);
    hash->cell_ext = cell_ext;
    memcpy(hash->cell_dim, cell_dim, sizeof(cell_dim));
    hash->magic = MD_SPATIAL_HASH_MAGIC;

    return true;
}

bool md_spatial_hash_free(md_spatial_hash_t* hash) {
    ASSERT(hash);
    ASSERT(hash->magic == MD_SPATIAL_HASH_MAGIC);
    md_array_free(hash->cells, hash->alloc);
    md_array_free(hash->coords, hash->alloc);
    md_array_free(hash->indices, hash->alloc);
    return true;
}

bool md_spatial_hash_query(const md_spatial_hash_t* hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    int dim = (int)(md_array_size(hash->coords) / md_array_size(hash->indices));
    ASSERT(dim > 0);

    vec4_t xyz = vec4_from_vec3(pos, 0);
    float rad = radius;
    float rad2 = radius * radius;

    vec4_t cell_min = vec4_from_vec3(hash->cell_min, 0);
    float  cell_ext = hash->cell_ext;
    int32_t cell_dim[3] = {hash->cell_dim[0], hash->cell_dim[1], hash->cell_dim[2]};

    vec4_t min_xyz = vec4_div_f(vec4_sub(vec4_sub_f(xyz, rad), cell_min), cell_ext);
    uint32_t min_coord[3] = {
        CLAMP((int32_t)(min_xyz.x), 0, (int32_t)cell_dim[0] - 1),
        CLAMP((int32_t)(min_xyz.y), 0, (int32_t)cell_dim[1] - 1),
        CLAMP((int32_t)(min_xyz.z), 0, (int32_t)cell_dim[2] - 1)
    };

    vec4_t max_xyz = vec4_div_f(vec4_sub(vec4_add_f(xyz, rad), cell_min), cell_ext);
    uint32_t max_coord[3] = {
        CLAMP((int32_t)(max_xyz.x), 0, (int32_t)cell_dim[0] - 1),
        CLAMP((int32_t)(max_xyz.y), 0, (int32_t)cell_dim[1] - 1),
        CLAMP((int32_t)(max_xyz.z), 0, (int32_t)cell_dim[2] - 1)
    };

    uint32_t cell_coord[3];
    for (cell_coord[2] = min_coord[2]; cell_coord[2] <= max_coord[2]; ++cell_coord[2]) {
        for (cell_coord[1] = min_coord[1]; cell_coord[1] <= max_coord[1]; ++cell_coord[1]) {
            for (cell_coord[0] = min_coord[0]; cell_coord[0] <= max_coord[0]; ++cell_coord[0]) {
                uint32_t cell_idx = cell_coord[2] * cell_dim[1] * cell_dim[0] + cell_coord[1] * cell_dim[0] + cell_coord[0];
                md_spatial_hash_cell_t cell = hash->cells[cell_idx];
                //vec4_t local_cell_min = vec4_add(cell_min, (vec4_t){cell_coord[0] * cell_ext, cell_coord[1] * cell_ext, cell_coord[2] * cell_ext, 0});
                for (uint32_t i = cell.offset; i < cell.offset + cell.length; ++i) {
                    //vec4_t p = decompress_coord(hash->coords + i * dim, local_cell_min, cell_ext);
                    vec4_t p = {0};
                    memcpy(p.elem, hash->coords + i * dim, sizeof(float) * dim);
                    if (vec4_distance_squared(p, xyz) < rad2) {
                        if (!iter(hash->indices[i], vec3_from_vec4(p), user_param)) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

static inline vec4_t vec4_deperiodize(vec4_t pos, vec4_t ref, vec4_t pbc_ext) {
    vec4_t d = vec4_sub(ref, pos);
    vec4_t m = (vec4_t){ .mm128 = md_simd_cmp_gt_f128(md_simd_abs_f128(d.mm128), md_simd_mul_f128(pbc_ext.mm128, md_simd_set1_f128(0.5f))) };
    vec4_t t = (vec4_t){ .mm128 = md_simd_blend_f128(md_simd_copysign_f128(pbc_ext.mm128, d.mm128), md_simd_zero_f128(), m.mm128) };
    return vec4_add(pos, t);
}

bool md_spatial_hash_query_periodic(const md_spatial_hash_t* hash, vec3_t in_pos, float in_rad, vec3_t in_pbc_min, vec3_t in_pbc_max, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    int dim = (int)(md_array_size(hash->coords) / md_array_size(hash->indices));
    ASSERT(dim > 0);

    vec4_t xyz = vec4_from_vec3(in_pos, 0);
    float rad = in_rad;
    float rad2 = rad * rad;

    //vec4_t pbc_ext = vec4_sub(vec4_from_vec3(in_pbc_max, 0), vec4_from_vec3(in_pbc_min, 0));

    vec4_t cell_min = vec4_from_vec3(hash->cell_min, 0);
    float  cell_ext = hash->cell_ext;
    int32_t cell_dim[3] = {hash->cell_dim[0], hash->cell_dim[1], hash->cell_dim[2]};

    vec4_t min_xyz = vec4_div_f(vec4_sub(vec4_sub_f(xyz, rad), cell_min), cell_ext);
    uint32_t min_coord[3] = {
        CLAMP((int32_t)(min_xyz.x), 0, (int32_t)cell_dim[0] - 1),
        CLAMP((int32_t)(min_xyz.y), 0, (int32_t)cell_dim[1] - 1),
        CLAMP((int32_t)(min_xyz.z), 0, (int32_t)cell_dim[2] - 1)
    };

    vec4_t max_xyz = vec4_div_f(vec4_sub(vec4_add_f(xyz, rad), cell_min), cell_ext);
    uint32_t max_coord[3] = {
        CLAMP((int32_t)(max_xyz.x), 0, (int32_t)cell_dim[0] - 1),
        CLAMP((int32_t)(max_xyz.y), 0, (int32_t)cell_dim[1] - 1),
        CLAMP((int32_t)(max_xyz.z), 0, (int32_t)cell_dim[2] - 1)
    };

    uint32_t cell_coord[3];
    for (cell_coord[2] = min_coord[2]; cell_coord[2] <= max_coord[2]; ++cell_coord[2]) {
        for (cell_coord[1] = min_coord[1]; cell_coord[1] <= max_coord[1]; ++cell_coord[1]) {
            for (cell_coord[0] = min_coord[0]; cell_coord[0] <= max_coord[0]; ++cell_coord[0]) {
                uint32_t cell_idx = cell_coord[2] * cell_dim[1] * cell_dim[0] + cell_coord[1] * cell_dim[0] + cell_coord[0];
                md_spatial_hash_cell_t cell = hash->cells[cell_idx];
                //vec4_t local_cell_min = vec4_add(cell_min, (vec4_t){cell_coord[0] * cell_ext, cell_coord[1] * cell_ext, cell_coord[2] * cell_ext, 0});
                for (uint32_t i = cell.offset; i < cell.offset + cell.length; ++i) { 
                    //vec4_t p = vec4_deperiodize(decompress_coord(hash->coords + i * dim, local_cell_min, cell_ext), xyz, pbc_ext);
                    vec4_t p = { 0 };
                    memcpy(p.elem, hash->coords + i * dim, sizeof(float)* dim);
                    if (vec4_distance_squared(p, xyz) < rad2) {
                        if (!iter(hash->indices[i], vec3_from_vec4(p), user_param)) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}