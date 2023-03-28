#include "md_spatial_hash.h"

#include "md_allocator.h"
#include "md_arena_allocator.h"
#include "md_log.h"
#include "md_array.h"
#include "md_simd.h"
#include "md_util.h"
#include "md_bitfield.h"

#include <float.h>
#include <string.h>
#include <stdint.h>

#define MD_SPATIAL_HASH_MAGIC 0xF67cbab6
#define CELL_EXT (6.0f)
#define INV_CELL_EXT (1.0f / CELL_EXT)
#define SQRT_CELL_EXT (8.4852813742385702f)
#define LENGTH_BITS 10
#define MAX_CELL_DIM 100000

#define SCL (CELL_EXT / 1023.f)

#define UNPACK_CELL_COORD(p) vec4_set( ((p >>  0) & 0x3FF) * SCL, ((p >> 10) & 0x3FF) * SCL, ((p >> 20) & 0x3FF) * SCL, 0)

static void compute_aabb(vec4_t* out_aabb_min, vec4_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, size_t count, size_t stride, vec4_t pbc_ext) {
    ASSERT(md_simd_f32_width == md_simd_i32_width);

    md_simd_f32_t vx_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vy_min = md_simd_set1_f32(+FLT_MAX);
    md_simd_f32_t vz_min = md_simd_set1_f32(+FLT_MAX);

    md_simd_f32_t vx_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vy_max = md_simd_set1_f32(-FLT_MAX);
    md_simd_f32_t vz_max = md_simd_set1_f32(-FLT_MAX);

    const vec4_t ext = pbc_ext;
    const vec4_t ref = vec4_mul_f(ext, 0.5f);
    const size_t simd_count = ROUND_DOWN(count, md_simd_f32_width);

    size_t i = 0;
    for (; i < simd_count; i += md_simd_f32_width) {
        md_simd_f32_t x = md_simd_load_f32((const float*)((const char*)in_x + i * stride));
        md_simd_f32_t y = md_simd_load_f32((const float*)((const char*)in_y + i * stride));
        md_simd_f32_t z = md_simd_load_f32((const float*)((const char*)in_z + i * stride));

        x = md_simd_deperiodize(x, md_simd_set1_f32(ref.x), md_simd_set1_f32(ext.x));
        y = md_simd_deperiodize(y, md_simd_set1_f32(ref.y), md_simd_set1_f32(ext.y));
        z = md_simd_deperiodize(z, md_simd_set1_f32(ref.z), md_simd_set1_f32(ext.z));

        vx_min = md_simd_min(vx_min, x);
        vy_min = md_simd_min(vy_min, y);
        vz_min = md_simd_min(vz_min, z);

        vx_max = md_simd_max(vx_max, x);
        vy_max = md_simd_max(vy_max, y);
        vz_max = md_simd_max(vz_max, z);
    }

    vec4_t aabb_min = {
        md_simd_hmin(vx_min),
        md_simd_hmin(vy_min),
        md_simd_hmin(vz_min),
        0
    };

    vec4_t aabb_max = {
        md_simd_hmax(vx_max),
        md_simd_hmax(vy_max),
        md_simd_hmax(vz_max),
        0
    };

    // Handle remainder
    for (; i < count; ++i) {
        vec4_t c = {
            *(const float*)((const char*)in_x + i * stride),
            *(const float*)((const char*)in_y + i * stride),
            *(const float*)((const char*)in_z + i * stride),
            0
        };
        c = vec4_deperiodize(c, ref, ext);

        aabb_min = vec4_min(aabb_min, c);
        aabb_max = vec4_max(aabb_max, c);
    }

    *out_aabb_min = aabb_min;
    *out_aabb_max = aabb_max;
}

static void compute_aabb_indexed(vec4_t* out_aabb_min, vec4_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, const int32_t* indices, size_t index_count, size_t stride, vec4_t pbc_ext) {
    const vec4_t ext = pbc_ext;
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min = vec4_set1(FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);

    // Handle remainder
    for (size_t i = 0; i < index_count; ++i) {
        const int32_t idx = indices[i];
        vec4_t c = {
            *(const float*)((const char*)in_x + idx * stride),
            *(const float*)((const char*)in_y + idx * stride),
            *(const float*)((const char*)in_z + idx * stride),
            0
        };
        c = vec4_deperiodize(c, ref, ext);

        aabb_min = vec4_min(aabb_min, c);
        aabb_max = vec4_max(aabb_max, c);
    }

    *out_aabb_min = aabb_min;
    *out_aabb_max = aabb_max;
}

static void compute_cell_coords(int32_t cell_min[3], int32_t cell_dim[3], uint32_t* cell_index, uint32_t* fract_coord, const float* in_x, const float* in_y, const float* in_z, size_t count, size_t stride, vec3_t pbc_ext) {
    ASSERT(md_simd_f32_width == md_simd_i32_width);

    const vec4_t ext = vec4_from_vec3(pbc_ext, 0);
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min, aabb_max;
    compute_aabb(&aabb_min, &aabb_max, in_x, in_y, in_z, count, stride, ext);

    // Compute the cell indices
    cell_min[0] = (int32_t)floorf(aabb_min.x * INV_CELL_EXT);
    cell_min[1] = (int32_t)floorf(aabb_min.y * INV_CELL_EXT);
    cell_min[2] = (int32_t)floorf(aabb_min.z * INV_CELL_EXT);

    const int32_t cell_max[3] = {
        (int32_t)floorf(aabb_max.x * INV_CELL_EXT) + 1,
        (int32_t)floorf(aabb_max.y * INV_CELL_EXT) + 1,
        (int32_t)floorf(aabb_max.z * INV_CELL_EXT) + 1,
    };

    cell_dim[0] = MAX(1, cell_max[0] - cell_min[0]);
    cell_dim[1] = MAX(1, cell_max[1] - cell_min[1]);
    cell_dim[2] = MAX(1, cell_max[2] - cell_min[2]);

    const size_t simd_count = ROUND_DOWN(count, md_simd_f32_width);
    const int32_t cell_dim_01 = cell_dim[0] * cell_dim[1];
    
    size_t i = 0;
    for (; i < simd_count; i += md_simd_f32_width) {
        md_simd_f32_t x = md_simd_load_f32((const float*)((const char*)in_x + i * stride));
        md_simd_f32_t y = md_simd_load_f32((const float*)((const char*)in_y + i * stride));
        md_simd_f32_t z = md_simd_load_f32((const float*)((const char*)in_z + i * stride));

        x = md_simd_deperiodize(x, md_simd_set1_f32(ref.x), md_simd_set1_f32(ext.x));
        y = md_simd_deperiodize(y, md_simd_set1_f32(ref.y), md_simd_set1_f32(ext.y));
        z = md_simd_deperiodize(z, md_simd_set1_f32(ref.z), md_simd_set1_f32(ext.z));

        x = md_simd_mul(x, md_simd_set1_f32(INV_CELL_EXT));
        y = md_simd_mul(y, md_simd_set1_f32(INV_CELL_EXT));
        z = md_simd_mul(z, md_simd_set1_f32(INV_CELL_EXT));

        md_simd_f32_t fx = md_simd_floor(x);
        md_simd_f32_t fy = md_simd_floor(y);
        md_simd_f32_t fz = md_simd_floor(z);

        md_simd_i32_t ccx = md_simd_sub(md_simd_convert_f32(fx), md_simd_set1_i32(cell_min[0]));
        md_simd_i32_t ccy = md_simd_sub(md_simd_convert_f32(fy), md_simd_set1_i32(cell_min[1]));
        md_simd_i32_t ccz = md_simd_sub(md_simd_convert_f32(fz), md_simd_set1_i32(cell_min[2]));

        md_simd_i32_t pcx = md_simd_convert_f32(md_simd_fmadd(md_simd_sub(x, fx), md_simd_set1_f32(1023.0f), md_simd_set1_f32(0.5f)));
        md_simd_i32_t pcy = md_simd_convert_f32(md_simd_fmadd(md_simd_sub(y, fy), md_simd_set1_f32(1023.0f), md_simd_set1_f32(0.5f)));
        md_simd_i32_t pcz = md_simd_convert_f32(md_simd_fmadd(md_simd_sub(z, fz), md_simd_set1_f32(1023.0f), md_simd_set1_f32(0.5f)));

        //md_simd_i32_t cc = md_simd_or(md_simd_shift_left(ccz, 20), md_simd_or(md_simd_shift_left(ccy, 10), ccx));
        
        // This is a bit of a hack to get around the multiplication for integers in SSE/AVX
        // The multiplication operations (of integers) of element sizes most often store the result in 2x element size
        // This means we loose our lanes if we perform the operation with epi32 -> result in epi64
        // Unless we resort to the mullo instruction which have half the throughput and twice the latency which is not ideal
        // Therefore we resport to 16-bit operations which yield 32-bit results.
        // 16-bits should be sufficient for the sizes which we operate on
#ifdef __AVX2__
        ccz = (md_simd_i32_t){_mm256_madd_epi16(ccz.m256i, _mm256_set1_epi32(cell_dim_01))};
        ccy = (md_simd_i32_t){_mm256_madd_epi16(ccy.m256i, _mm256_set1_epi32(cell_dim[0]))};
#elif defined(__AVX__)
        // Have to emulate the wide instruction for AVX
        ccz = (md_simd_i32_t){MD_SIMD_DOUBLE_PUMP2_SI128(ccz.m256i, _mm256_set1_epi32(cell_dim_01), _mm_madd_epi16)};
        ccy = (md_simd_i32_t){MD_SIMD_DOUBLE_PUMP2_SI128(ccy.m256i, _mm256_set1_epi32(cell_dim[0]), _mm_madd_epi16)};
#elif defined(__x86_64__)
        // SSE2
        ccz = (md_simd_i32_t){_mm_madd_epi16(ccz.m128i, _mm_set1_epi32(cell_dim_01))};
        ccy = (md_simd_i32_t){_mm_madd_epi16(ccy.m128i, _mm_set1_epi32(cell_dim[0]))};
#else
#error "Unsupported architecture :<"
#endif

        md_simd_i32_t ci = md_simd_add(ccz, md_simd_add(ccy, ccx));
        md_simd_i32_t pc = md_simd_or(md_simd_shift_left(pcz, 20), md_simd_or(md_simd_shift_left(pcy, 10), pcx));

        md_simd_store((int*)cell_index  + i, ci);
        md_simd_store((int*)fract_coord + i, pc);
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
        uint32_t ci = ((int32_t)f.z - cell_min[2]) * cell_dim_01 + ((int32_t)f.y - cell_min[1]) * cell_dim[0] + ((int32_t)f.x - cell_min[0]);

        cell_index[i] = ci;
        fract_coord[i] = ((uint32_t)p.z << 20) | ((uint32_t)p.y << 10) | (uint32_t)p.x;
    }
}

static void compute_cell_coords_indexed(int32_t cell_min[3], int32_t cell_dim[3], uint32_t* cell_index, uint32_t* fract_coord, const float* in_x, const float* in_y, const float* in_z, const int* indices, size_t count, size_t stride, vec3_t pbc_ext) {
    ASSERT(md_simd_f32_width == md_simd_i32_width);

    const vec4_t ext = vec4_from_vec3(pbc_ext, 0);
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min, aabb_max;
    compute_aabb_indexed(&aabb_min, &aabb_max, in_x, in_y, in_z, indices, count, stride, ext);

    // Compute the cell indices
    cell_min[0] = (int32_t)floorf(aabb_min.x * INV_CELL_EXT);
    cell_min[1] = (int32_t)floorf(aabb_min.y * INV_CELL_EXT);
    cell_min[2] = (int32_t)floorf(aabb_min.z * INV_CELL_EXT);

    const int32_t cell_max[3] = {
        (int32_t)floorf(aabb_max.x * INV_CELL_EXT) + 1,
        (int32_t)floorf(aabb_max.y * INV_CELL_EXT) + 1,
        (int32_t)floorf(aabb_max.z * INV_CELL_EXT) + 1,
    };

    cell_dim[0] = MAX(1, cell_max[0] - cell_min[0]);
    cell_dim[1] = MAX(1, cell_max[1] - cell_min[1]);
    cell_dim[2] = MAX(1, cell_max[2] - cell_min[2]);

    const int32_t cell_dim_01 = cell_dim[0] * cell_dim[1];

    size_t dst_idx = 0;
    for (size_t i = 0; i < count; ++i) {
        int src_idx = indices[i];
        vec4_t c = {
            *(const float*)((const char*)in_x + src_idx * stride),
            *(const float*)((const char*)in_y + src_idx * stride),
            *(const float*)((const char*)in_z + src_idx * stride),
            0
        };
        c = vec4_mul_f(vec4_deperiodize(c, ref, ext), INV_CELL_EXT);

        const vec4_t f = vec4_floor(c);
        const vec4_t p = vec4_fmadd(vec4_sub(c, f), vec4_set1(1023.0f), vec4_set1(0.5f));
        const uint32_t ci = ((int32_t)f.z - cell_min[2]) * cell_dim_01 + ((int32_t)f.y - cell_min[1]) * cell_dim[0] + ((int32_t)f.x - cell_min[0]);

        cell_index [dst_idx] = ci;
        fract_coord[dst_idx] = ((uint32_t)p.z << 20) | ((uint32_t)p.y << 10) | (uint32_t)p.x;

        dst_idx += 1;
    }

    ASSERT(dst_idx == count);
}

static bool init(md_spatial_hash_t* hash, const float* in_x, const float* in_y, const float* in_z, const int* indices, size_t count, size_t stride, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash);
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(alloc);

    bool result = false;

    stride = MAX(stride, sizeof(float));

    md_allocator_i* temp_alloc = default_allocator;
    size_t bytes = sizeof(uint32_t) * count * 3;
    void* mem = md_alloc(temp_alloc, bytes);

    uint32_t* fract_cell_coord  = (uint32_t*)mem;
    uint32_t* cell_index        = (uint32_t*)mem + count;
    uint32_t* local_idx         = (uint32_t*)mem + count * 2;

    int32_t cell_min[3];
    int32_t cell_dim[3];

    if (indices) {
        compute_cell_coords_indexed(cell_min, cell_dim, cell_index, fract_cell_coord, in_x, in_y, in_z, indices, count, stride, pbc_ext);
    } else {
        compute_cell_coords(cell_min, cell_dim, cell_index, fract_cell_coord, in_x, in_y, in_z, count, stride, pbc_ext);
    }

    if (cell_dim[0] > MAX_CELL_DIM ||
        cell_dim[1] > MAX_CELL_DIM ||
        cell_dim[2] > MAX_CELL_DIM)
    {
        MD_LOG_ERROR("Spatial hash cell dimension is too large: {%i %i %i}", cell_dim[0], cell_dim[1], cell_dim[2]);
        goto done;
    }

    const uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Allocate needed data
    uint32_t* cells         = md_alloc(alloc, sizeof(uint32_t) * cell_count);
    uint32_t* coords        = md_alloc(alloc, sizeof(uint32_t) * count);
    uint32_t* original_idx  = md_alloc(alloc, sizeof(uint32_t) * count);

    MEMSET(cells, 0, sizeof(uint32_t) * cell_count);

    for (size_t i = 0; i < count; ++i) {
        uint32_t idx = cell_index[i];
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
    if (indices) {
        for (size_t i = 0; i < count; ++i) {
            uint32_t dst_idx        = (cells[cell_index[i]] >> 10) + local_idx[i];
            coords[dst_idx]         = fract_cell_coord[i];
            original_idx[dst_idx]   = (uint32_t)indices[i];
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            uint32_t dst_idx        = (cells[cell_index[i]] >> 10) + local_idx[i];
            coords[dst_idx]         = fract_cell_coord[i];
            original_idx[dst_idx]   = (uint32_t)i;
        }
    }

#if 0
    for (int32_t z = 0; z < cell_dim[2]; z++) {
        for (int32_t y = 0; y < cell_dim[1]; y++) {
            for (int32_t x = 0; x < cell_dim[0]; x++) {
                const int32_t ci = z * cell_dim[1] * cell_dim[0] + y * cell_dim[0] + x;
                const uint32_t cell_data = cells[ci];
                const uint32_t cell_offset = cell_data >> 10;
                const uint32_t cell_length = cell_data & 1023;
                vec4_t cc_min = vec4_set((cell_min[0] + x) * CELL_EXT, (cell_min[1] + y) * CELL_EXT, (cell_min[2] + z) * CELL_EXT, 0);
                for (uint32_t i = cell_offset; i < cell_offset + cell_length; ++i) {
                    const vec4_t p = vec4_add(cc_min, UNPACK_CELL_COORD(coords[i]));

                    int32_t idx = original_idx[i];
                    const vec4_t ref_pos = vec4_set(
                            *(const float*)((const char*)in_x + idx * stride),
                            *(const float*)((const char*)in_y + idx * stride),
                            *(const float*)((const char*)in_z + idx * stride),
                            0);
                
                    const float d2 = vec4_periodic_distance_squared(p, ref_pos, vec4_from_vec3(pbc_ext, 0.0f));
                    if (d2 > 0.1f) {
                        ASSERT(false);
                    }
                }
            }
        }
    }
#endif

    hash->cells    = cells;
    hash->coords   = coords;
    hash->indices  = original_idx;
    hash->alloc    = alloc;
    MEMCPY(hash->cell_min, cell_min, sizeof(cell_min));
    hash->coord_count = (uint32_t)count;
    MEMCPY(hash->cell_dim, cell_dim, sizeof(cell_dim));
    hash->magic = MD_SPATIAL_HASH_MAGIC;
    hash->pbc_ext = vec4_from_vec3(pbc_ext, 0);

    result = true;
done:
    md_free(temp_alloc, mem, bytes);
    return result;
}

bool md_spatial_hash_init(md_spatial_hash_t* hash, const vec3_t* coord, int64_t count, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(coord);
    ASSERT(alloc);

    if (count < 0) {
        MD_LOG_ERROR("Invalid count");
        return false;
    }

    return init(hash, &coord->x, &coord->y, &coord->z, NULL, (size_t)count, sizeof(vec3_t), pbc_ext, alloc);
}

bool md_spatial_hash_init_indexed(md_spatial_hash_t* hash, const vec3_t* coord, const int32_t* indices, int64_t index_count, vec3_t pbc_ext, struct md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(coord);
    ASSERT(alloc);

    if (index_count < 0) {
        MD_LOG_ERROR("Invalid count");
        return false;
    }

    return init(hash, &coord->x, &coord->y, &coord->z, indices, (size_t)index_count, sizeof(vec3_t), pbc_ext, alloc);
}

bool md_spatial_hash_init_soa(md_spatial_hash_t* hash, const float* x, const float* y, const float* z, int64_t count, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(alloc);

    if (count < 0) {
        MD_LOG_ERROR("Invalid count");
        return false;
    }

    return init(hash, x, y, z, NULL, (size_t)count, sizeof(float), pbc_ext, alloc);
}

bool md_spatial_hash_init_indexed_soa(md_spatial_hash_t* hash, const float* x, const float* y, const float* z, const int32_t* indices, int64_t index_count, vec3_t pbc_ext, md_allocator_i* alloc) {
    ASSERT(hash); 
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);
    ASSERT(alloc);

    if (index_count < 0) {
        MD_LOG_ERROR("Invalid count");
        return false;
    }

    return init(hash, x, y, z, indices, (size_t)index_count, sizeof(float), pbc_ext, alloc);
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

    int32_t cc[3] = {0};
    vec4_t cc_min = {0};
    for (cc[2] = cell_beg[2]; cc[2] < cell_end[2]; ++cc[2]) {
        cc_min.z = (cell_min[2] + cc[2]) * CELL_EXT;
        for (cc[1] = cell_beg[1]; cc[1] < cell_end[1]; ++cc[1]) {
            cc_min.y = (cell_min[1] + cc[1]) * CELL_EXT;
            for (cc[0] = cell_beg[0]; cc[0] < cell_end[0]; ++cc[0]) {
                cc_min.x = (cell_min[0] + cc[0]) * CELL_EXT;
                const uint32_t cell_idx = cc[2] * cell_dim[1] * cell_dim[0] + cc[1] * cell_dim[0] + cc[0];
                const uint32_t cell_data = hash->cells[cell_idx];
                const uint32_t cell_offset = cell_data >> 10;
                const uint32_t cell_length = cell_data & 1023;
                for (uint32_t i = cell_offset; i < cell_offset + cell_length; ++i) {
                    const vec4_t p = vec4_add(cc_min, UNPACK_CELL_COORD(hash->coords[i]));
                    const float d2 = vec4_distance_squared(p, pos);
                    if (d2 < rad2) {
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
                const uint32_t cell_idx = cc[2] * cell_dim[1] * cell_dim[0] + cc[1] * cell_dim[0] + cc[0];
                const uint32_t cell_data = hash->cells[cell_idx];
                const uint32_t cell_offset = cell_data >> 10;
                const uint32_t cell_length = cell_data & 1023;
                for (uint32_t i = cell_offset; i < cell_offset + cell_length; ++i) {
                    const vec4_t p = vec4_add(cc_min, UNPACK_CELL_COORD(hash->coords[i]));
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

bool validate_spatial_hash(const md_spatial_hash_t* spatial_hash) {
    if (!spatial_hash) {
        MD_LOG_ERROR("spatial_hash is null");
        return false;
    }
    
    if (spatial_hash->magic != MD_SPATIAL_HASH_MAGIC) {
        MD_LOG_ERROR("spatial_hash has invalid magic");
        return false;
    }

    return true;
}

bool md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(spatial_hash);
    ASSERT(iter);

    if (!validate_spatial_hash(spatial_hash)) {
        return false;
    }

    // This is for non periodic lookup
    if (vec4_equal(spatial_hash->pbc_ext, vec4_zero())) {
        return query_pos_rad(spatial_hash, pos, radius, iter, user_param);
    } else {
        return query_pos_rad_periodic(spatial_hash, pos, radius, iter, user_param);
    }
}

typedef struct {
    md_array(uint32_t) arr;
    md_allocator_i* alloc;
} idx_data_t;

bool idx_fn(uint32_t idx, vec3_t pos, void* user_param) {
    (void)pos;
    md_bitfield_t* bf = (md_bitfield_t*)user_param;
    md_bitfield_set_bit(bf, idx);
    return true;
}

int64_t md_spatial_hash_query_idx(int32_t* buf, int64_t cap, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius) {
    md_bitfield_t bf = md_bitfield_create(default_temp_allocator);
    md_bitfield_reserve_range(&bf, 0, spatial_hash->coord_count);

    md_spatial_hash_query(spatial_hash, pos, radius, idx_fn, &bf);
#if DEBUG
    if ((int64_t)md_bitfield_popcount(&bf) > cap) {
        MD_LOG_DEBUG("Size exceeds the capacity, some elements will not be written");
    }
#endif

    return md_bitfield_extract_indices(buf, cap, &bf);
}

void md_spatial_hash_query_bits(md_bitfield_t* bf, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius) {
    md_spatial_hash_query(spatial_hash, pos, radius, idx_fn, bf);
}
