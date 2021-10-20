#include "md_spatial_hash.h"

#include "md_allocator.h"
#include "md_log.h"
#include "md_array.inl"

#include <float.h>
#include <string.h>

#define MD_SPATIAL_HASH_MAGIC 0xF67cbab672578178
#define MD_SPATIAL_HASH_DEFAULT_CELL_EXT 6.0
#define MD_SPATIAL_HASH_MIN_CELL_EXT 3.0

typedef struct md_spatial_hash_cell_t {
    uint32_t offset : 22; // ~4.2 M
    uint32_t length : 10; // 1024
} md_spatial_hash_cell_t;

// Store as relative cell-coordinates quantized to 10 bits
typedef struct md_spatial_hash_coord_t {
    uint32_t x : 10;    // 1024
    uint32_t y : 10;    // 1024
    uint32_t z : 10;    // 1024
    uint32_t unused : 2;
} md_spatial_hash_coord_t;

static inline vec3_t decompress_coord(md_spatial_hash_coord_t coord, vec3_t cell_min, float cell_ext) {
    const float scl = cell_ext / (1 << 10);
    return vec3_add(cell_min, vec3_mul_f((vec3_t){(float)coord.x, (float)coord.y, (float)coord.z}, scl));
}

static inline md_spatial_hash_coord_t compress_coord(vec3_t xyz, vec3_t cell_min, float cell_ext) {
    const float scl = (1 << 10) / cell_ext;
    xyz = vec3_sub(xyz, cell_min);
    xyz = vec3_mul_f(xyz, scl);
    xyz = vec3_add_f(xyz, 0.5f);
    md_spatial_hash_coord_t coord = {(uint32_t)xyz.x, (uint32_t)xyz.y, (uint32_t)xyz.z, 0};
    return coord;
}

bool md_spatial_hash_init(md_spatial_hash_t* hash, const float* in_x, const float* in_y, const float* in_z, int64_t count, float cell_ext, struct md_allocator_i* alloc) {
    ASSERT(hash);
    ASSERT(alloc);

    if (count == 0) return false;
    if (count < 0) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid count");
        return false;
    }

    if (cell_ext < 0) {
        md_printf(MD_LOG_TYPE_ERROR, "Cell Extent parameter (%g) must be positive.", cell_ext);
        return false;
    }

    if (cell_ext == 0) cell_ext = MD_SPATIAL_HASH_DEFAULT_CELL_EXT;
    cell_ext = MAX(cell_ext, MD_SPATIAL_HASH_MIN_CELL_EXT);

    if (hash->coords != 0 || hash->cells != 0) {
        md_print(MD_LOG_TYPE_ERROR, "Spatial Hash structure is not set to zero, leaking potential memory here.");
    }

    vec3_t aabb_min = (vec3_t) { FLT_MAX,  FLT_MAX,  FLT_MAX};
    vec3_t aabb_max = (vec3_t) {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    for (int64_t i = 0; i < count; ++i) {
        aabb_min.x = MIN(aabb_min.x, in_x[i]);
        aabb_min.y = MIN(aabb_min.y, in_y[i]);
        aabb_min.z = MIN(aabb_min.z, in_z[i]);

        aabb_max.x = MAX(aabb_max.x, in_x[i]);
        aabb_max.y = MAX(aabb_max.y, in_y[i]);
        aabb_max.z = MAX(aabb_max.z, in_z[i]);
    }

    // Try to pick a good cell_ext
    vec3_t aabb_mid = vec3_lerp(aabb_min, aabb_max, 0.5f);
    vec3_t aabb_ext = vec3_sub(aabb_max, aabb_min);

    uint32_t cell_dim[3] = {
        (uint32_t)(aabb_ext.x / cell_ext) + 1,
        (uint32_t)(aabb_ext.y / cell_ext) + 1,
        (uint32_t)(aabb_ext.z / cell_ext) + 1
    };

    ASSERT(cell_dim[0] > 0);
    ASSERT(cell_dim[1] > 0);
    ASSERT(cell_dim[2] > 0);

    vec3_t  full_ext = {cell_dim[0] * cell_ext, cell_dim[1] * cell_ext, cell_dim[2] * cell_ext};
    vec3_t  cell_min = vec3_sub(aabb_mid, vec3_mul_f(full_ext, 0.5f));
    uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Allocate needed data
    md_spatial_hash_cell_t*  cells = 0;
    md_spatial_hash_coord_t* coords = 0;
    uint32_t* original_idx = 0;
    md_array_resize(cells, cell_count, alloc);
    md_array_resize(coords, count, alloc);
    md_array_resize(original_idx, count, alloc);

    memset(cells, 0, cell_count * sizeof(md_spatial_hash_cell_t));

    uint32_t* cell_idx  = 0;
    uint32_t* local_idx = 0;
    md_array_resize(cell_idx,  count, default_temp_allocator);
    md_array_resize(local_idx, count, default_temp_allocator);

    const float inv_cell_ext = 1.0f / cell_ext;

    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        uint32_t cell_coord[3] = {
            (uint32_t)((in_x[i] - cell_min.x) * inv_cell_ext),
            (uint32_t)((in_y[i] - cell_min.y) * inv_cell_ext),
            (uint32_t)((in_z[i] - cell_min.z) * inv_cell_ext)
        };

        uint32_t idx = cell_coord[2] * cell_dim[0] * cell_dim[1] + cell_coord[1] * cell_dim[0] + cell_coord[0];
        cell_idx[i]  = idx;
        local_idx[i] = cells[idx].length++;
        ASSERT(cells[idx].length != 0 && "hash cell has wrapped, too many entities per cell");
    }

    for (uint32_t i = 1; i < cell_count; ++i) {
        cells[i].offset = cells[i-1].offset + cells[i-1].length;
    }

    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        uint32_t dst_idx = cells[cell_idx[i]].offset + local_idx[i];
        uint32_t cell_coord[3] = {
            (uint32_t)((in_x[i] - cell_min.x) * inv_cell_ext),
            (uint32_t)((in_y[i] - cell_min.y) * inv_cell_ext),
            (uint32_t)((in_z[i] - cell_min.z) * inv_cell_ext)
        };
        vec3_t xyz = {in_x[i], in_y[i], in_z[i]};
        vec3_t local_cell_min = vec3_add(cell_min, (vec3_t){cell_coord[0] * cell_ext, cell_coord[1] * cell_ext, cell_coord[2] * cell_ext});

        ASSERT(dst_idx < count);
        coords[dst_idx] = compress_coord(xyz, local_cell_min, cell_ext);
        original_idx[dst_idx] = i;
    }

    md_array_free(local_idx, default_temp_allocator);
    md_array_free(cell_idx,  default_temp_allocator);

    md_spatial_hash_cell_t* cell = md_array_last(cells);
    ASSERT(cell);
    ASSERT(cell->offset + cell->length == count);

    hash->cells  = cells;
    hash->coords = coords;
    hash->index  = original_idx;
    hash->cell_min = cell_min;
    hash->cell_ext = cell_ext;
    hash->cell_dim[0] = cell_dim[0];
    hash->cell_dim[1] = cell_dim[1];
    hash->cell_dim[2] = cell_dim[2];
    hash->magic = MD_SPATIAL_HASH_MAGIC;

    return true;
}

// Allot of this is redundant, clean this up at some point
bool md_spatial_hash_init_vec3(md_spatial_hash_t* hash, const vec3_t* xyz, int64_t count, float cell_ext, struct md_allocator_i* alloc) {
    ASSERT(hash);
    ASSERT(alloc);

    if (count == 0) return false;
    if (count < 0) {
        md_print(MD_LOG_TYPE_ERROR, "Invalid count");
        return false;
    }

    if (cell_ext < 0) {
        md_printf(MD_LOG_TYPE_ERROR, "Cell Extent parameter (%g) must be positive.", cell_ext);
        return false;
    }

    if (cell_ext == 0) cell_ext = MD_SPATIAL_HASH_DEFAULT_CELL_EXT;

    if (hash->coords != 0 || hash->cells != 0) {
        md_print(MD_LOG_TYPE_ERROR, "Spatial Hash structure is not set to zero, leaking potential memory here.");
    }

    vec3_t aabb_min = (vec3_t) { FLT_MAX,  FLT_MAX,  FLT_MAX};
    vec3_t aabb_max = (vec3_t) {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    for (int64_t i = 0; i < count; ++i) {
        aabb_min.x = MIN(aabb_min.x, xyz[i].x);
        aabb_min.y = MIN(aabb_min.y, xyz[i].y);
        aabb_min.z = MIN(aabb_min.z, xyz[i].z);

        aabb_max.x = MAX(aabb_max.x, xyz[i].x);
        aabb_max.y = MAX(aabb_max.y, xyz[i].y);
        aabb_max.z = MAX(aabb_max.z, xyz[i].z);
    }

    // Try to pick a good cell_ext
    vec3_t aabb_mid = vec3_lerp(aabb_min, aabb_max, 0.5f);
    vec3_t aabb_ext = vec3_sub(aabb_max, aabb_min);

    uint32_t cell_dim[3] = {
        (uint32_t)(aabb_ext.x / cell_ext) + 1,
        (uint32_t)(aabb_ext.y / cell_ext) + 1,
        (uint32_t)(aabb_ext.z / cell_ext) + 1
    };

    ASSERT(cell_dim[0] > 0);
    ASSERT(cell_dim[1] > 0);
    ASSERT(cell_dim[2] > 0);

    vec3_t  full_ext = {cell_dim[0] * cell_ext, cell_dim[1] * cell_ext, cell_dim[2] * cell_ext};
    vec3_t  cell_min = vec3_sub(aabb_mid, vec3_mul_f(full_ext, 0.5f));
    uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Allocate needed data
    md_spatial_hash_cell_t*  cells = 0;
    md_spatial_hash_coord_t* coords = 0;
    uint32_t* original_idx = 0;
    md_array_resize(cells, cell_count, alloc);
    md_array_resize(coords, count, alloc);
    md_array_resize(original_idx, count, alloc);

    memset(cells, 0, cell_count * sizeof(md_spatial_hash_cell_t));

    uint32_t* cell_idx  = 0;
    uint32_t* local_idx = 0;
    md_array_resize(cell_idx,  count, default_temp_allocator);
    md_array_resize(local_idx, count, default_temp_allocator);

    const float inv_cell_ext = 1.0f / cell_ext;

    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        uint32_t cell_coord[3] = {
            (uint32_t)((xyz[i].x - cell_min.x) * inv_cell_ext),
            (uint32_t)((xyz[i].y - cell_min.y) * inv_cell_ext),
            (uint32_t)((xyz[i].z - cell_min.z) * inv_cell_ext)
        };

        uint32_t idx = cell_coord[2] * cell_dim[0] * cell_dim[1] + cell_coord[1] * cell_dim[0] + cell_coord[0];
        cell_idx[i]  = idx;
        local_idx[i] = cells[idx].length++;
        ASSERT(cells[idx].length != 0 && "hash cell has wrapped, too many entities per cell");
    }

    for (uint32_t i = 1; i < cell_count; ++i) {
        cells[i].offset = cells[i-1].offset + cells[i-1].length;
    }

    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        uint32_t dst_idx = cells[cell_idx[i]].offset + local_idx[i];
        uint32_t cell_coord[3] = {
            (uint32_t)((xyz[i].x - cell_min.x) * inv_cell_ext),
            (uint32_t)((xyz[i].y - cell_min.y) * inv_cell_ext),
            (uint32_t)((xyz[i].z - cell_min.z) * inv_cell_ext)
        };
        vec3_t local_cell_min = vec3_add(cell_min, (vec3_t){cell_coord[0] * cell_ext, cell_coord[1] * cell_ext, cell_coord[2] * cell_ext});

        ASSERT(dst_idx < count);
        coords[dst_idx] = compress_coord(xyz[i], local_cell_min, cell_ext);
        original_idx[dst_idx] = i;
    }

    md_array_free(local_idx, default_temp_allocator);
    md_array_free(cell_idx,  default_temp_allocator);

    md_spatial_hash_cell_t* cell = md_array_last(cells);
    ASSERT(cell);
    ASSERT(cell->offset + cell->length == count);

    hash->cells  = cells;
    hash->coords = coords;
    hash->index  = original_idx;
    hash->cell_min = cell_min;
    hash->cell_ext = cell_ext;
    hash->cell_dim[0] = cell_dim[0];
    hash->cell_dim[1] = cell_dim[1];
    hash->cell_dim[2] = cell_dim[2];
    hash->magic = MD_SPATIAL_HASH_MAGIC;

    return true;
}

bool md_spatial_hash_free(md_spatial_hash_t* hash, struct md_allocator_i* alloc) {
    ASSERT(hash);
    ASSERT(alloc);
    ASSERT(hash->magic == MD_SPATIAL_HASH_MAGIC);
    md_array_free(hash->cells, alloc);
    md_array_free(hash->coords, alloc);
    return true;
}

bool md_spatial_hash_query(const md_spatial_hash_t* hash, vec4_t pos_rad, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    vec3_t xyz = vec3_from_vec4(pos_rad);
    float rad = pos_rad.w;
    float rad2 = pos_rad.w * pos_rad.w;

    vec3_t min_xyz = vec3_div_f(vec3_sub(vec3_sub_f(xyz, rad), hash->cell_min), hash->cell_ext);
    uint32_t min_coord[3] = {
        CLAMP((int32_t)(min_xyz.x), 0, (int32_t)hash->cell_dim[0] - 1),
        CLAMP((int32_t)(min_xyz.y), 0, (int32_t)hash->cell_dim[1] - 1),
        CLAMP((int32_t)(min_xyz.z), 0, (int32_t)hash->cell_dim[2] - 1)
    };

    vec3_t max_xyz = vec3_div_f(vec3_sub(vec3_add_f(xyz, rad), hash->cell_min), hash->cell_ext);
    uint32_t max_coord[3] = {
        CLAMP((int32_t)(max_xyz.x), 0, (int32_t)hash->cell_dim[0] - 1),
        CLAMP((int32_t)(max_xyz.y), 0, (int32_t)hash->cell_dim[1] - 1),
        CLAMP((int32_t)(max_xyz.z), 0, (int32_t)hash->cell_dim[2] - 1)
    };

    uint32_t cell_coord[3];
    for (cell_coord[2] = min_coord[2]; cell_coord[2] <= max_coord[2]; ++cell_coord[2]) {
        for (cell_coord[1] = min_coord[1]; cell_coord[1] <= max_coord[1]; ++cell_coord[1]) {
            for (cell_coord[0] = min_coord[0]; cell_coord[0] <= max_coord[0]; ++cell_coord[0]) {
                uint32_t cell_idx = cell_coord[2] * hash->cell_dim[1] * hash->cell_dim[0] + cell_coord[1] * hash->cell_dim[0] + cell_coord[0];
                md_spatial_hash_cell_t cell = hash->cells[cell_idx];
                vec3_t local_cell_min = vec3_add(hash->cell_min, (vec3_t){cell_coord[0] * hash->cell_ext, cell_coord[1] * hash->cell_ext, cell_coord[2] * hash->cell_ext});
                for (uint32_t i = cell.offset; i < cell.offset + cell.length; ++i) {
                    vec3_t p = decompress_coord(hash->coords[i], local_cell_min, hash->cell_ext);
                    if (vec3_distance_squared(p, xyz) < rad2) {
                        if (!iter(hash->index[i], p, user_param)) {
                            return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

bool md_spatial_hash_query_n(const md_spatial_hash_t* hash, vec4_t* pos_rad, int64_t count, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(pos_rad);
    ASSERT(iter);

    bool result = true;

    // At some point, it will be faster to just brute force check all particles individually
    for (int64_t i = 0; i < count; ++i) {
        result = md_spatial_hash_query(hash, pos_rad[i], iter, user_param);
        if (!result) break;
    }

    // The strategy here is to use two bitfields with the length of the total number of cells
    // Then we initialize the first by setting all bits for all cells which have a length > 0
    // Then we iterate over pos_rad and set all cell bits in the second field which correspond to the cells it 'touches'
    // Then we perform an AND operation on the two bitfields
    // As a last step we iterate over the bitset and only check the marked cells for potential...
    // 
    // This will only work if the 

    return result;
}