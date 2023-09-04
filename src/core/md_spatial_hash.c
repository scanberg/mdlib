#include <core/md_spatial_hash.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_simd.h>
#include <core/md_bitfield.h>

#include <md_util.h>

#include <float.h>
#include <string.h>
#include <stdint.h>

#define MD_SPATIAL_HASH_MAGIC 0xF67cbab6
#define CELL_EXT (6.0f)
#define INV_CELL_EXT (1.0f / CELL_EXT)
#define SQRT_CELL_EXT (8.4852813742385702f)
#define LENGTH_BITS 10
#define LENGTH_MASK ((1 << LENGTH_BITS) - 1)
#define OFFSET_BITS (32 - LENGTH_BITS)
#define LENGTH_CAP (1 << LENGTH_BITS)
#define MAX_CELL_DIM 1023

#define SCL (CELL_EXT / 1023.f)

#define UNPACK_CELL_COORD(p) vec4_set( ((p >>  0) & 0x3FF) * SCL, ((p >> 10) & 0x3FF) * SCL, ((p >> 20) & 0x3FF) * SCL, 0)

typedef md_spatial_hash_elem_t elem_t;

typedef uint32_t cell_t;

typedef struct md_spatial_hash_t {
    vec4_t pbc_ext;
    int32_t cell_min[3];
    uint32_t elem_count;
    int32_t cell_dim[3];
    uint32_t magic;
    void* data;
    md_allocator_i* alloc;
} md_spatial_hash_t;

/*
struct md_spatial_acc_t {
    vec4_t cell_ext;
    md_unit_cell_t unit_cell;
    md_array(elem_t) elems;
	md_array(cell_t) cells;
    md_allocator_i* alloc;
    int32_t cell_min[3];
    int32_t cell_dim[3];
};
*/

static void compute_aabb_vec3(vec4_t* out_aabb_min, vec4_t* out_aabb_max, const vec3_t* in_xyz, const int32_t* indices, int64_t count, vec4_t pbc_ext) {
    const vec4_t ext = pbc_ext;
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    if (count == 0) {
        *out_aabb_min = vec4_zero();
        *out_aabb_max = vec4_zero();
        return;
    }

    vec4_t aabb_min = vec4_set1(+FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);
    
    for (int64_t i = 0; i < count; ++i) {
        int32_t idx = indices ? indices[i] : (int32_t)i;
        vec4_t c = vec4_from_vec3(in_xyz[idx], 0);
        c = vec4_deperiodize(c, ref, ext);
        aabb_min = vec4_min(aabb_min, c);
        aabb_max = vec4_max(aabb_max, c);
    }

    *out_aabb_min = aabb_min;
    *out_aabb_max = aabb_max;
}

static void compute_aabb_soa(vec4_t* out_aabb_min, vec4_t* out_aabb_max, const float* in_x, const float* in_y, const float* in_z, const int32_t* indices, int64_t count, vec4_t pbc_ext) {
    const vec4_t ext = pbc_ext;
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    if (count == 0) {
        *out_aabb_min = vec4_zero();
        *out_aabb_max = vec4_zero();
        return;
    }

    vec4_t aabb_min = vec4_set1(FLT_MAX);
    vec4_t aabb_max = vec4_set1(-FLT_MAX);

    for (int64_t i = 0; i < count; ++i) {
        const int32_t idx = indices ? indices[i] : (int32_t)i;
        vec4_t c = { in_x[idx], in_y[idx], in_z[idx], 0 };
        c = vec4_deperiodize(c, ref, ext);
        aabb_min = vec4_min(aabb_min, c);
        aabb_max = vec4_max(aabb_max, c);
    }

    *out_aabb_min = aabb_min;
    *out_aabb_max = aabb_max;
}

/*
md_spatial_acc_t* md_spatial_acc_create_vec3(const vec3_t* in_xyz, const int32_t* indices, int64_t count, const md_unit_cell_t* unit_cell, md_allocator_i* alloc) {
    ASSERT(in_xyz);
    ASSERT(alloc);
    ASSERT(count >= 0);

    if (count == 0) {
        MD_LOG_INFO("Count is zero, no spatial acceleration structure was created");
        return NULL;
    }

    uint32_t* cell_index  = md_alloc(md_temp_allocator, sizeof(uint32_t) * count);
    uint16_t* local_index = md_alloc(md_temp_allocator, sizeof(uint16_t) * count);

    ASSERT(cell_index);
    ASSERT(local_index);

    md_spatial_acc_t acc = {0};
    acc.alloc = alloc;

    vec4_t pbc_ext = {0};
    if (unit_cell) {
        acc.unit_cell = *unit_cell;
        pbc_ext = vec4_from_vec3(mat3_mul_vec3(unit_cell->basis, vec3_set1(1.0f)), 0.0f);
    }

    const vec4_t ext = pbc_ext;
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min, aabb_max;
    compute_aabb_vec3(&aabb_min, &aabb_max, in_xyz, indices, count, ext);
    const vec4_t pad = vec4_set(0.1f, 0.1f, 0.1f, 0.0f);
    aabb_min = vec4_sub(aabb_min, pad);
    aabb_max = vec4_add(aabb_max, pad);

    // We want to have cells ~CELL_EXT 
    vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);
    vec4_t cell_est = vec4_ceil(vec4_div_f(aabb_ext, CELL_EXT));
    acc.cell_ext = vec4_div(aabb_ext, cell_est);

    acc.cell_min[0] = (int)floorf(aabb_min.x / acc.cell_ext.x);
    acc.cell_min[1] = (int)floorf(aabb_min.y / acc.cell_ext.y);
    acc.cell_min[2] = (int)floorf(aabb_min.z / acc.cell_ext.z);

    int cell_max[3];
    cell_max[0] = (int)ceilf(aabb_max.x / acc.cell_ext.x) + 1;
    cell_max[1] = (int)ceilf(aabb_max.y / acc.cell_ext.y) + 1;
    cell_max[2] = (int)ceilf(aabb_max.z / acc.cell_ext.z) + 1;

    acc.cell_dim[0] = (int)MAX(1, cell_max[0] - acc.cell_min[0]);
    acc.cell_dim[1] = (int)MAX(1, cell_max[1] - acc.cell_min[1]);
    acc.cell_dim[2] = (int)MAX(1, cell_max[2] - acc.cell_min[2]);

    const int cd_0  = acc.cell_dim[0];
    const int cd_01 = acc.cell_dim[0] * acc.cell_dim[1];
    for (int64_t i = 0; i < count; ++i) {
        const int32_t idx = indices ? indices[i] : (int32_t)i;
        vec4_t coord = vec4_from_vec3(in_xyz[idx], 0);
        coord = vec4_deperiodize(coord, ref, ext);

        vec4_t cell = vec4_div(coord, acc.cell_ext);
        vec4_t cc = vec4_floor(cell);

        uint32_t cx = ((int)cc.x - acc.cell_min[0]);
        uint32_t cy = ((int)cc.y - acc.cell_min[1]);
        uint32_t cz = ((int)cc.z - acc.cell_min[2]);
        cell_index[i] = (uint32_t)(cx + cy * cd_0 + cz * cd_01);
    }

    const int64_t cell_count = acc.cell_dim[0] * acc.cell_dim[1] * acc.cell_dim[2];
    acc.cells = md_array_create(cell_t, cell_count, alloc);
    MEMSET(acc.cells, 0, sizeof(uint32_t) * cell_count);

    for (int64_t i = 0; i < count; ++i) {
        uint32_t idx = cell_index[i];
        ASSERT(idx < cell_count);
        local_index[i] = (uint16_t)acc.cells[idx]++;
        ASSERT(acc.cells[idx] < LENGTH_CAP && "Too many entities per cell");
    }

    // Prefix sum the offset
    uint32_t offset = acc.cells[0];
    for (int64_t i = 1; i < cell_count; ++i) {
        const uint32_t length = acc.cells[i];
        acc.cells[i] = (offset << 10) | length;
        offset += length;
    }

#if DEBUG
    uint32_t cell = acc.cells[cell_count - 1];
    ASSERT((cell >> 10) + (cell & 1023) == count);
#endif

    acc.elems = md_array_create(elem_t, count, alloc);

    for (int64_t i = 0; i < count; ++i) {
        const int64_t ci = cell_index[i];
        const int64_t dst_idx = (acc.cells[ci] >> 10) + local_index[i];
        const int64_t src_idx = indices ? indices[i] : i;
        vec4_t coord = vec4_deperiodize(vec4_from_vec3(in_xyz[src_idx], 0), ref, ext);
        acc.elems[dst_idx] = (elem_t){ coord.x, coord.y, coord.z, (uint32_t)src_idx };
        //acc.elems[dst_idx] = (elem_t){ (uint32_t)src_idx };
    }

    md_free(md_temp_allocator, local_index, sizeof(uint16_t) * count);
    md_free(md_temp_allocator, cell_index,  sizeof(uint32_t) * count);

    md_spatial_acc_t* ptr = md_alloc(alloc, sizeof(md_spatial_acc_t));
    MEMCPY(ptr, &acc, sizeof(md_spatial_acc_t));

    return ptr;
}

md_spatial_acc_t* md_spatial_acc_create_soa(const float* in_x, const float* in_y, const float* in_z, const int32_t* indices, int64_t count, const md_unit_cell_t* unit_cell, md_allocator_i* alloc) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(alloc);
    ASSERT(count >= 0);

    if (count == 0) {
        MD_LOG_INFO("Count is zero, no spatial acceleration structure was created");
        return NULL;
    }

    uint32_t* cell_index  = md_alloc(md_temp_allocator, sizeof(uint32_t) * count);
    uint16_t* local_index = md_alloc(md_temp_allocator, sizeof(uint16_t) * count);

    ASSERT(cell_index);
    ASSERT(local_index);

    md_spatial_acc_t acc = {0};
    acc.alloc = alloc;

    vec4_t pbc_ext = {0};
    if (unit_cell) {
        acc.unit_cell = *unit_cell;
        pbc_ext = vec4_from_vec3(mat3_mul_vec3(unit_cell->basis, vec3_set1(1.0f)), 0.0f);
    }

    const vec4_t ext = pbc_ext;
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min, aabb_max;
    compute_aabb_soa(&aabb_min, &aabb_max, in_x, in_y, in_z, indices, count, ext);
    //const vec4_t pad = vec4_set(0.1f, 0.1f, 0.1f, 0.0f);
    //aabb_min = vec4_sub(aabb_min, pad);
    //aabb_max = vec4_add(aabb_max, pad);

    // We want to have cells ~CELL_EXT 
    vec4_t aabb_ext = vec4_sub(aabb_max, aabb_min);
    vec4_t cell_est = vec4_ceil(vec4_div_f(aabb_ext, CELL_EXT));
    acc.cell_ext = vec4_blend(vec4_div(aabb_ext, cell_est), vec4_set1(CELL_EXT), vec4_cmp_eq(cell_est, vec4_zero()));

    vec4_t c_min = vec4_floor(vec4_div(aabb_min, acc.cell_ext));
    vec4_t c_max = vec4_ceil(vec4_div(aabb_max,  acc.cell_ext));

    acc.cell_min[0] = (int)c_min.x;
    acc.cell_min[1] = (int)c_min.y;
    acc.cell_min[2] = (int)c_min.z;

    int cell_max[3];
    cell_max[0] = (int)c_max.x + 1;
    cell_max[1] = (int)c_max.y + 1;
    cell_max[2] = (int)c_max.z + 1;

    acc.cell_dim[0] = (int)MAX(1, cell_max[0] - acc.cell_min[0]);
    acc.cell_dim[1] = (int)MAX(1, cell_max[1] - acc.cell_min[1]);
    acc.cell_dim[2] = (int)MAX(1, cell_max[2] - acc.cell_min[2]);

    const int64_t cell_count = acc.cell_dim[0] * acc.cell_dim[1] * acc.cell_dim[2];
    acc.cells = md_array_create(cell_t, cell_count, alloc);
    MEMSET(acc.cells, 0, sizeof(uint32_t) * cell_count);

    const uint32_t cd_0  = acc.cell_dim[0];
    const uint32_t cd_01 = acc.cell_dim[0] * acc.cell_dim[1];
    for (int64_t i = 0; i < count; ++i) {
        const int32_t idx = indices ? indices[i] : (int32_t)i;
        vec4_t coord = vec4_set(in_x[idx], in_y[idx], in_z[idx], 0);
        coord = vec4_deperiodize(coord, ref, ext);

        vec4_t cell = vec4_div(coord, acc.cell_ext);
        vec4_t cc = vec4_floor(cell);

        uint32_t cx = ((int)cc.x - acc.cell_min[0]);
        uint32_t cy = ((int)cc.y - acc.cell_min[1]);
        uint32_t cz = ((int)cc.z - acc.cell_min[2]);
        uint32_t ci = (uint32_t)(cx + cy * cd_0 + cz * cd_01);
        uint16_t li = (uint16_t)acc.cells[ci]++;
        
        cell_index[i]  = ci;
        local_index[i] = li;
    }

    // Prefix sum the offset
    uint32_t offset = acc.cells[0];
    for (int64_t i = 1; i < cell_count; ++i) {
        const uint32_t length = acc.cells[i];
        acc.cells[i] = (offset << 10) | length;
        offset += length;
    }

#if DEBUG
    uint32_t cell = acc.cells[cell_count - 1];
    ASSERT((cell >> LENGTH_BITS) + (cell & LENGTH_MASK) == count);
#endif

    acc.elems = md_array_create(elem_t, count, alloc);
    
    for (int64_t i = 0; i < count; ++i) {
        const int64_t ci = cell_index[i];
        const int64_t dst_idx = (acc.cells[ci] >> 10) + local_index[i];
        const int64_t src_idx = indices ? indices[i] : i;
        vec4_t coord = vec4_deperiodize(vec4_set(in_x[src_idx], in_y[src_idx], in_z[src_idx], 0), ref, ext);
        acc.elems[dst_idx] = (elem_t){ coord.x, coord.y, coord.z, (uint32_t)src_idx };
        //acc.elems[dst_idx] = (elem_t){ (uint32_t)src_idx };
    }

    md_free(md_temp_allocator, local_index, sizeof(uint16_t) * count);
    md_free(md_temp_allocator, cell_index,  sizeof(uint32_t) * count);

    md_spatial_acc_t* ptr = md_alloc(alloc, sizeof(md_spatial_acc_t));
    MEMCPY(ptr, &acc, sizeof(md_spatial_acc_t));

    return ptr;
}

void md_spatial_acc_free(md_spatial_acc_t* acc) {
    if (acc) {
        ASSERT(acc->alloc);
        md_allocator_i* alloc = acc->alloc;
		md_array_free(acc->cells, alloc);
        md_array_free(acc->elems, alloc);
        MEMSET(acc, 0, sizeof(md_spatial_acc_t));
		md_free(alloc, acc, sizeof(md_spatial_acc_t));
    }
}
*/

/*
md_spatial_acc_iter_t md_spatial_acc_query(const md_spatial_acc_t* acc, vec3_t pos, float radius) {
    ASSERT(acc);
    md_spatial_acc_iter_t it = {0};
    it.acc = acc;

    const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(acc->unit_cell.basis, vec3_set1(1)), 0.0f);
    const vec4_t pbc_ref = vec4_mul_f(pbc_ext, 0.5f);

    const vec4_t cell_min = vec4_div(vec4_deperiodize(vec4_sub_f(vec4_from_vec3(pos, 0), radius), pbc_ref, pbc_ext), acc->cell_ext);
    const vec4_t cell_max = vec4_div(vec4_deperiodize(vec4_add_f(vec4_from_vec3(pos, 0), radius), pbc_ref, pbc_ext), acc->cell_ext);
    
    it.cc_beg[0] = MAX((int32_t)cell_min.x, acc->cell_min[0]);
    it.cc_beg[1] = MAX((int32_t)cell_min.y, acc->cell_min[1]);
    it.cc_beg[2] = MAX((int32_t)cell_min.z, acc->cell_min[2]);

    int32_t cc_end[3] = {
        (int32_t)cell_max.x + 1,
        (int32_t)cell_max.y + 1,
        (int32_t)cell_max.z + 1,
    };
    
    // Don't clamp here, because we want to handle periodic conditions as well.
    it.cc_ext[0] = cc_end[0] - it.cc_beg[0];
    it.cc_ext[1] = cc_end[1] - it.cc_beg[1];
    it.cc_ext[2] = cc_end[2] - it.cc_beg[2];

    // Set start coordinate 
    it.cc[0] = it.cc_beg[0] - 1;  // -1 because we want to start with the first cell after we perform next on the iterator
    it.cc[1] = it.cc_beg[1];
    it.cc[2] = it.cc_beg[2];
    
    it.cell_idx = 0;
    it.cell_idx_end = 0;
    it.elem_idx = 0;

    it.query_pos = pos;
    it.query_rad2 = radius * radius;
    
    return it;
}

static inline bool next_cell(md_spatial_acc_iter_t* iter) {
    iter->cc[0] += 1;

    if (iter->cc[0] < iter->cc_ext[0]) {
        return true;
    }
    
    iter->cc[0] = iter->cc_beg[0];
    iter->cc[1] += 1;

    if (iter->cc[1] < iter->cc_ext[1]) {
        return true;
    }

    iter->cc[1] = iter->cc_beg[1];
    iter->cc[2] += 1;
    
    return (iter->cc[2] != iter->cc_ext[2]);
}

static inline elem_t* get_next_potential_elem(md_spatial_acc_iter_t* iter) {
    if (iter->cell_idx + 1 < iter->cell_idx_end) {
        iter->cell_idx += 1;
    } else if (next_cell(iter)) {
            // Goto next cell
            const uint32_t cx = (iter->cc_beg[0] + iter->cc[0]) - iter->acc->cell_min[0];
            const uint32_t cy = (iter->cc_beg[1] + iter->cc[1]) - iter->acc->cell_min[1];
            const uint32_t cz = (iter->cc_beg[2] + iter->cc[2]) - iter->acc->cell_min[2];
            const uint32_t ci = cz * (iter->acc->cell_dim[0] * iter->acc->cell_dim[1]) + cy * iter->acc->cell_dim[0] + cx;
            const uint32_t cell_value = iter->acc->cells[ci];
            iter->cell_idx_end = (cell_value >> LENGTH_BITS) + (cell_value & LENGTH_MASK);
            iter->cell_idx     = (cell_value >> LENGTH_BITS);
    } else {
        return NULL;
    }
    
    return &iter->acc->elems[iter->cell_idx];
}

bool md_spatial_acc_iter_next(md_spatial_acc_iter_t* iter) {
    ASSERT(iter);
    ASSERT(iter->acc);
    
    elem_t* elem;
    while ((elem = get_next_potential_elem(iter))) {
        vec4_t p1 = vec4_from_vec3(elem->coord, 0);
        vec4_t p2 = vec4_from_vec3(iter->query_pos, 0);
        const float d2 = vec4_distance_squared(p1, p2);
        if (d2 < iter->query_rad2) {
            iter->elem_idx = elem->idx;
            iter->elem_pos = elem->coord;
            return true;
        }
    }
    
    return false;
}
*/

md_spatial_hash_t* md_spatial_hash_create_vec3(const vec3_t* in_xyz, const int32_t* in_indices, int64_t count, const md_unit_cell_t* unit_cell, md_allocator_i* alloc) {
    ASSERT(in_xyz);
    ASSERT(alloc);

    if (count < 0) {
        MD_LOG_ERROR("Invalid count");
        return NULL;
    }

    md_spatial_hash_t* hash = NULL;

    md_allocator_i* temp_alloc = md_heap_allocator;
    size_t temp_bytes = sizeof(uint32_t) * count * 2;
    void* temp_mem = md_alloc(temp_alloc, temp_bytes);

    uint32_t* cell_index = (uint32_t*)temp_mem + count * 0;
    uint32_t* local_idx  = (uint32_t*)temp_mem + count * 1;

    int32_t cell_min[3];
    int32_t cell_dim[3];

    ASSERT(md_simd_f32_width == md_simd_i32_width);
    const vec4_t ext = unit_cell ? vec4_from_vec3(mat3_mul_vec3(unit_cell->basis, vec3_set1(1.0f)), 0) : vec4_zero();
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min, aabb_max;
    compute_aabb_vec3(&aabb_min, &aabb_max, in_xyz, in_indices, count, ext);

    const vec4_t c_min = vec4_floor(vec4_div_f(aabb_min, CELL_EXT));
    const vec4_t c_max = vec4_floor(vec4_div_f(aabb_max, CELL_EXT));

    // Compute the cell indices
    cell_min[0] = (int32_t)c_min.x;
    cell_min[1] = (int32_t)c_min.y;
    cell_min[2] = (int32_t)c_min.z;

    const int32_t cell_max[3] = {
        (int32_t)c_max.x + 1,
        (int32_t)c_max.y + 1,
        (int32_t)c_max.z + 1,
    };

    cell_dim[0] = MAX(1, cell_max[0] - cell_min[0]);
    cell_dim[1] = MAX(1, cell_max[1] - cell_min[1]);
    cell_dim[2] = MAX(1, cell_max[2] - cell_min[2]);

    if (cell_dim[0] > MAX_CELL_DIM ||
        cell_dim[1] > MAX_CELL_DIM ||
        cell_dim[2] > MAX_CELL_DIM)
    {
        MD_LOG_ERROR("Spatial hash cell dimension is too large: {%i %i %i}", cell_dim[0], cell_dim[1], cell_dim[2]);
        goto done;
    }

    const uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Allocate needed data
    const int64_t tot_bytes = sizeof(md_spatial_hash_t) + sizeof(elem_t) * ROUND_UP(count, md_simd_f32_width) + sizeof(cell_t) * cell_count;
    void* mem = md_alloc(alloc, tot_bytes);
    void* data = (char*)mem + sizeof(md_spatial_hash_t);

    hash = mem;
    elem_t* elem_data = data;
    cell_t* cell_data = (cell_t*)((char*)data + sizeof(elem_t) * ROUND_UP(count, md_simd_f32_width));
    MEMSET(cell_data, 0, sizeof(cell_t) * cell_count);

    const int32_t cell_dim_01 = cell_dim[0] * cell_dim[1];

    // Handle remainder
    for (int64_t i = 0; i < count; ++i) {
        const int64_t idx = in_indices ? in_indices[i] : i;
        const vec4_t coord = vec4_from_vec3(in_xyz[idx], 0);
        vec4_t cell = vec4_mul_f(vec4_deperiodize(coord, ref, ext), INV_CELL_EXT);
        vec4_t whole = vec4_floor(cell);
        
        int32_t cz = MAX(0, (int32_t)whole.z - cell_min[2]);
        int32_t cy = MAX(0, (int32_t)whole.y - cell_min[1]);
        int32_t cx = MAX(0, (int32_t)whole.x - cell_min[0]);
        uint32_t ci = cz * cell_dim_01 + cy * cell_dim[0] + cx;

        cell_index[i] = ci;
        ASSERT(ci < cell_count);

        local_idx[i]  = cell_data[ci]++;
        ASSERT(cell_data[ci] < LENGTH_CAP && "Too many entities per cell");
    }

    uint32_t offset = cell_data[0];
    for (uint32_t i = 1; i < cell_count; ++i) {
        const uint32_t length = cell_data[i];
        if (length) {
            cell_data[i] = (offset << 10) | length;
        }
        offset += length;
    }

    for (int64_t i = 0; i < count; ++i) {
        const int64_t cell_idx = cell_index[i];
        const int64_t src_idx  = in_indices ? in_indices[i] : i;
        const int64_t dst_idx  = (cell_data[cell_idx] >> LENGTH_BITS) + local_idx[i];
        const vec4_t coord     = vec4_deperiodize(vec4_from_vec3(in_xyz[src_idx], 0), ref, ext);
        elem_data[dst_idx]     = (elem_t){coord.x, coord.y, coord.z, (uint32_t)src_idx};
    }

    hash->pbc_ext = ext;
    MEMCPY(hash->cell_min, cell_min, sizeof(cell_min));
    hash->elem_count = (uint32_t)count;
    MEMCPY(hash->cell_dim, cell_dim, sizeof(cell_dim));
    hash->magic = MD_SPATIAL_HASH_MAGIC;
    hash->data  = data;
    hash->alloc = alloc;

done:
    md_free(temp_alloc, temp_mem, temp_bytes);
    return hash;
}

md_spatial_hash_t* md_spatial_hash_create_soa(const float* in_x, const float* in_y, const float* in_z, const int32_t* in_indices, int64_t count, const md_unit_cell_t* unit_cell, md_allocator_i* alloc) {
    ASSERT(in_x);
    ASSERT(in_y);
    ASSERT(in_z);
    ASSERT(alloc);

    if (count < 0) {
        MD_LOG_ERROR("Invalid count");
        return NULL;
    }

    md_spatial_hash_t* hash = NULL;

    md_allocator_i* temp_alloc = md_heap_allocator;
    size_t temp_bytes = sizeof(uint32_t) * count * 2;
    void* temp_mem = md_alloc(temp_alloc, temp_bytes);

    uint32_t* cell_index = (uint32_t*)temp_mem + count * 0;
    uint32_t* local_idx  = (uint32_t*)temp_mem + count * 1;

    int32_t cell_min[3];
    int32_t cell_dim[3];

    ASSERT(md_simd_f32_width == md_simd_i32_width);

    const vec4_t ext = unit_cell ? vec4_from_vec3(mat3_mul_vec3(unit_cell->basis, vec3_set1(1.0f)), 0) : vec4_zero();
    const vec4_t ref = vec4_mul_f(ext, 0.5f);

    vec4_t aabb_min, aabb_max;
    compute_aabb_soa(&aabb_min, &aabb_max, in_x, in_y, in_z, in_indices, count, ext);

    vec4_t c_min = vec4_floor(vec4_div_f(aabb_min, CELL_EXT));
    vec4_t c_max = vec4_floor(vec4_div_f(aabb_max, CELL_EXT));

    // Compute the cell indices
    cell_min[0] = (int32_t)c_min.x;
    cell_min[1] = (int32_t)c_min.y;
    cell_min[2] = (int32_t)c_min.z;

    const int32_t cell_max[3] = {
        (int32_t)c_max.x + 1,
        (int32_t)c_max.y + 1,
        (int32_t)c_max.z + 1,
    };

    cell_dim[0] = MAX(1, cell_max[0] - cell_min[0]);
    cell_dim[1] = MAX(1, cell_max[1] - cell_min[1]);
    cell_dim[2] = MAX(1, cell_max[2] - cell_min[2]);

    if (cell_dim[0] > MAX_CELL_DIM ||
        cell_dim[1] > MAX_CELL_DIM ||
        cell_dim[2] > MAX_CELL_DIM)
    {
        MD_LOG_ERROR("Spatial hash cell dimension is too large: {%i %i %i}", cell_dim[0], cell_dim[1], cell_dim[2]);
        goto done;
    }

    const uint32_t cell_count = cell_dim[0] * cell_dim[1] * cell_dim[2];

    // Allocate needed data
    const int64_t tot_bytes = sizeof(md_spatial_hash_t) + sizeof(elem_t) * ROUND_UP(count, md_simd_f32_width) + sizeof(cell_t) * cell_count;
    void* mem = md_alloc(alloc, tot_bytes);
    void* data = (char*)mem + sizeof(md_spatial_hash_t);

    hash = mem;
    elem_t* elem_data = data;
    cell_t* cell_data = (cell_t*)((char*)data + sizeof(elem_t) * ROUND_UP(count, md_simd_f32_width));
    MEMSET(cell_data, 0, sizeof(cell_t) * cell_count);

    const int32_t cell_dim_01 = cell_dim[0] * cell_dim[1];

    /*
    size_t i = 0;
    const size_t simd_count = ROUND_DOWN(count, md_simd_f32_width);
    for (; i < simd_count; i += md_simd_f32_width) {
    md_simd_f32_t x = md_simd_load_f32((const float*)((const char*)in_x + i * stride));
    md_simd_f32_t y = md_simd_load_f32((const float*)((const char*)in_y + i * stride));
    md_simd_f32_t z = md_simd_load_f32((const float*)((const char*)in_z + i * stride));

    x = md_simd_deperiodize(x, md_simd_set1_f32(ref.x), md_simd_set1_f32(ext.x));
    y = md_simd_deperiodize(y, md_simd_set1_f32(ref.y), md_simd_set1_f32(ext.y));
    z = md_simd_deperiodize(z, md_simd_set1_f32(ref.z), md_simd_set1_f32(ext.z));

    // Cell relative coordinates
    md_simd_f32_t cx = md_simd_mul(x, md_simd_set1_f32(INV_CELL_EXT));
    md_simd_f32_t cy = md_simd_mul(y, md_simd_set1_f32(INV_CELL_EXT));
    md_simd_f32_t cz = md_simd_mul(z, md_simd_set1_f32(INV_CELL_EXT));

    md_simd_f32_t flx = md_simd_floor(cx);
    md_simd_f32_t fly = md_simd_floor(cy);
    md_simd_f32_t flz = md_simd_floor(cz);

    md_simd_i32_t ccx = md_simd_sub(md_simd_convert_f32(flx), md_simd_set1_i32(cell_min[0]));
    md_simd_i32_t ccy = md_simd_sub(md_simd_convert_f32(fly), md_simd_set1_i32(cell_min[1]));
    md_simd_i32_t ccz = md_simd_sub(md_simd_convert_f32(flz), md_simd_set1_i32(cell_min[2]));

    md_simd_f32_t frx = md_simd_sub(cx, flx);
    md_simd_f32_t fry = md_simd_sub(cy, fly);
    md_simd_f32_t frz = md_simd_sub(cz, flz);

    md_simd_i32_t pcx = md_simd_convert_f32(md_simd_fmadd(frx, md_simd_set1_f32(1023.0f), md_simd_set1_f32(0.5f)));
    md_simd_i32_t pcy = md_simd_convert_f32(md_simd_fmadd(fry, md_simd_set1_f32(1023.0f), md_simd_set1_f32(0.5f)));
    md_simd_i32_t pcz = md_simd_convert_f32(md_simd_fmadd(frz, md_simd_set1_f32(1023.0f), md_simd_set1_f32(0.5f)));

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

    md_simd_f32_t xx = _mm256_cvtepi32_ps(pcx.m256i);
    md_simd_f32_t yy = _mm256_cvtepi32_ps(pcy.m256i);
    md_simd_f32_t zz = _mm256_cvtepi32_ps(pcz.m256i);

    md_simd_f32_t qx = md_simd_add(md_simd_mul(flx, md_simd_set1_f32(CELL_EXT)), md_simd_mul(xx, md_simd_set1_f32(SCL)));
    md_simd_f32_t qy = md_simd_add(md_simd_mul(fly, md_simd_set1_f32(CELL_EXT)), md_simd_mul(yy, md_simd_set1_f32(SCL)));
    md_simd_f32_t qz = md_simd_add(md_simd_mul(flz, md_simd_set1_f32(CELL_EXT)), md_simd_mul(zz, md_simd_set1_f32(SCL)));

    md_simd_f32_t dx = md_simd_abs(md_simd_sub(x, qx));
    md_simd_f32_t dy = md_simd_abs(md_simd_sub(y, qy));
    md_simd_f32_t dz = md_simd_abs(md_simd_sub(z, qz));

    md_simd_store((int*)cell_index  + i, ci);
    md_simd_store((int*)fract_coord + i, pc);
    }
    */

    // Handle remainder
    for (int64_t i = 0; i < count; ++i) {
        const int64_t idx = in_indices ? in_indices[i] : i;
        const vec4_t coord = vec4_set(in_x[idx], in_y[idx], in_z[idx], 0);
        vec4_t cell = vec4_mul_f(vec4_deperiodize(coord, ref, ext), INV_CELL_EXT);

        vec4_t whole = vec4_floor(cell);
        //vec4_t fract = vec4_sub(cell, whole);
        //vec4_t packed = vec4_fmadd(fract, vec4_set1(1023.0f), vec4_set1(0.5f));
        int32_t cz = MAX(0, (int32_t)whole.z - cell_min[2]);
        int32_t cy = MAX(0, (int32_t)whole.y - cell_min[1]);
        int32_t cx = MAX(0, (int32_t)whole.x - cell_min[0]);
        uint32_t ci = cz * cell_dim_01 + cy * cell_dim[0] + cx;

        cell_index[i] = ci;
        ASSERT(ci < cell_count);

        local_idx[i]  = cell_data[ci]++;
        ASSERT(cell_data[ci] < LENGTH_CAP && "Too many entities per cell");
        
        //fract_coord[i] = ((uint32_t)packed.z << 20) | ((uint32_t)packed.y << 10) | (uint32_t)packed.x;
    }

    uint32_t offset = cell_data[0];
    for (uint32_t i = 1; i < cell_count; ++i) {
        const uint32_t length = cell_data[i];
        if (length) {
            cell_data[i] = (offset << 10) | length;
        }
        offset += length;
    }

    for (int64_t i = 0; i < count; ++i) {
        const int64_t cell_idx  = cell_index[i];
        const int64_t src_idx   = in_indices ? in_indices[i] : i;
        const int64_t dst_idx   = (cell_data[cell_idx] >> LENGTH_BITS) + local_idx[i];
        const vec4_t coord      = vec4_deperiodize(vec4_set(in_x[src_idx], in_y[src_idx], in_z[src_idx], 0), ref, ext);
        elem_data[dst_idx]      = (elem_t){coord.x, coord.y, coord.z, (uint32_t)src_idx};
    }
    
    hash->pbc_ext = ext;
    MEMCPY(hash->cell_min, cell_min, sizeof(cell_min));
    hash->elem_count = (uint32_t)count;
    MEMCPY(hash->cell_dim, cell_dim, sizeof(cell_dim));
    hash->magic = MD_SPATIAL_HASH_MAGIC;
    hash->data  = data;
    hash->alloc = alloc;

done:
    md_free(temp_alloc, temp_mem, temp_bytes);
    return hash;
}

void md_spatial_hash_free(md_spatial_hash_t* hash) {
    ASSERT(hash);
    ASSERT(hash->magic == MD_SPATIAL_HASH_MAGIC);
    ASSERT(hash->alloc);
    md_allocator_i* alloc = hash->alloc;
    const int64_t size = sizeof(md_spatial_hash_t) + sizeof(cell_t) * hash->cell_dim[0] * hash->cell_dim[1] * hash->cell_dim[2] + sizeof(elem_t) * hash->elem_count;
    // Only zero the hash, not the actual array fields
    MEMSET(hash, 0, sizeof(md_spatial_hash_t));
    md_free(alloc, hash, size);
}

// Non periodic version (simpler)
static inline void query_pos_rad(const md_spatial_hash_t* hash, vec3_t position, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    vec4_t pos = vec4_from_vec3(position, 0);
    float rad = radius;
    float rad2 = radius * radius;

    const int32_t* cell_min = hash->cell_min;
    const int32_t* cell_dim = hash->cell_dim;
    const elem_t* elems = hash->data;
    const cell_t* cells = (const cell_t*)((const char*)hash->data + sizeof(elem_t) * ROUND_UP(hash->elem_count, md_simd_f32_width));

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
    for (cc[2] = cell_beg[2]; cc[2] < cell_end[2]; ++cc[2]) {
        for (cc[1] = cell_beg[1]; cc[1] < cell_end[1]; ++cc[1]) {
            for (cc[0] = cell_beg[0]; cc[0] < cell_end[0]; ++cc[0]) {
                const uint32_t cell_idx = cc[2] * cell_dim[1] * cell_dim[0] + cc[1] * cell_dim[0] + cc[0];
                const uint32_t cell_data = cells[cell_idx];
                if (cell_data == 0) {
                    continue;
                }
                const uint32_t cell_offset = cell_data >> LENGTH_BITS;
                const uint32_t cell_length = cell_data  & LENGTH_MASK;
                for (uint32_t i = cell_offset; i < cell_offset + cell_length; ++i) {
                    const elem_t* elem = &elems[i];
                    const vec4_t p = vec4_from_vec3(elem->xyz, 0);
                    const float d2 = vec4_distance_squared(p, pos);
                    if (d2 < rad2) {
                        if (!iter(elem, user_param)) {
                            return;
                        }
                    }
                }
            }
        }
    }
}

static inline void query_pos_rad_batch(const md_spatial_hash_t* hash, vec3_t position, float radius, md_spatial_hash_batch_iter_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    const vec4_t pos = vec4_from_vec3(position, 0);
    float rad2 = radius * radius;

    const int32_t* cell_min = hash->cell_min;
    const int32_t* cell_dim = hash->cell_dim;
    const elem_t* elems = hash->data;
    const cell_t* cells = (const cell_t*)((const char*)hash->data + sizeof(elem_t) * ROUND_UP(hash->elem_count, md_simd_f32_width));

    const int32_t cc_beg_x = (int32_t)(floorf((pos.x - radius) * INV_CELL_EXT)) - cell_min[0];
    const int32_t cc_beg_y = (int32_t)(floorf((pos.y - radius) * INV_CELL_EXT)) - cell_min[1];
    const int32_t cc_beg_z = (int32_t)(floorf((pos.z - radius) * INV_CELL_EXT)) - cell_min[2];

    const int32_t cc_end_x = (int32_t)(ceilf((pos.x + radius) * INV_CELL_EXT)) - cell_min[0];
    const int32_t cc_end_y = (int32_t)(ceilf((pos.y + radius) * INV_CELL_EXT)) - cell_min[1];
    const int32_t cc_end_z = (int32_t)(ceilf((pos.z + radius) * INV_CELL_EXT)) - cell_min[2];

    const int32_t cc_beg[3] = {
        CLAMP(cc_beg_x, 0, cell_dim[0] - 1),
        CLAMP(cc_beg_y, 0, cell_dim[1] - 1),
        CLAMP(cc_beg_z, 0, cell_dim[2] - 1)
    };

    const int32_t cc_end[3] = {
        CLAMP(cc_end_x, 0, cell_dim[0]),
        CLAMP(cc_end_y, 0, cell_dim[1]),
        CLAMP(cc_end_z, 0, cell_dim[2])
    };

    const int32_t cd_0  = cell_dim[0];
    const int32_t cd_01 = cell_dim[0] * cell_dim[1];

    int32_t cc[3];
    for (cc[2] = cc_beg[2]; cc[2] < cc_end[2]; ++cc[2]) {
        const int32_t cz = cc[2] * cd_01;
        for (cc[1] = cc_beg[1]; cc[1] < cc_end[1]; ++cc[1]) {
            const int32_t czy = cz + cc[1] * cd_0;
            for (cc[0] = cc_beg[0]; cc[0] < cc_end[0]; ++cc[0]) {
                const uint32_t cell_idx  = czy + cc[0];
                const uint32_t cell_data = cells[cell_idx];

                if (cell_data == 0) {
                    continue;
                }
                const uint32_t cell_offset = cell_data >> LENGTH_BITS;
                const uint32_t cell_length = cell_data  & LENGTH_MASK;

                int len = cell_length;
                const elem_t* elem = elems + cell_offset;

                while (len > 0) {
                    md_simd_f32_t vx,vy,vz;
                    md_simd_unpack_xyz(&vx, &vy, &vz, (const float*)elem, sizeof(elem_t));
                    md_simd_f32_t dx = md_simd_sub(vx, md_simd_set1_f32(pos.x));
                    md_simd_f32_t dy = md_simd_sub(vy, md_simd_set1_f32(pos.y));
                    md_simd_f32_t dz = md_simd_sub(vz, md_simd_set1_f32(pos.z));
                    md_simd_f32_t d2 = md_simd_add(md_simd_add(md_simd_mul(dx, dx), md_simd_mul(dy, dy)), md_simd_mul(dz, dz));
                    md_simd_f32_t vmask = md_simd_cmp_lt(d2, md_simd_set1_f32(rad2));

                    const int step = MIN(len, md_simd_f32_width);
                    const int lane_mask = (1 << step) - 1;
                    const int mask = md_simd_movemask(vmask) & lane_mask;

                    if (!iter(elem, mask, user_param)) {
                        return;
                    }

                    len  -= md_simd_f32_width;
                    elem += md_simd_f32_width;
                };
            }
        }
    }
}

static inline void inc_cell(int* ci, int beg, int dim, int pbc) {
    *ci += 1;
    int cc = beg + *ci;
    if (cc >= dim) {
        *ci += pbc - dim;
    }
}

static void inc_cell2(int* cc, int max, int jmp) {
    *cc += 1;
    if (*cc == max) {
        *cc += jmp;
    }
}

// Periodic version (more complex)
// The periodic domain is aligned to the cells at the origin, but not at the other end since we use a fixed cell-size
// Furthermore, the cell domain usually only spans a subrange of the periodic domain
// This makes the query non-trivial as in the non periodic case.

static inline void query_pos_rad_periodic(const md_spatial_hash_t* hash, vec3_t position, float rad, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    const vec4_t ref = vec4_mul_f(hash->pbc_ext, 0.5f);
    const vec4_t pbc_ext = hash->pbc_ext;
    const vec4_t pos = vec4_from_vec3(position, 0);
    float rad2 = rad * rad;

    const int32_t* cell_min = hash->cell_min;
    const int32_t* cell_dim = hash->cell_dim;
    const elem_t* elems = hash->data;
    const cell_t* cells = (const cell_t*)((const char*)hash->data + sizeof(elem_t) * ROUND_UP(hash->elem_count, md_simd_f32_width));

    const int cell_max[3] = {
        cell_min[0] + cell_dim[0],
        cell_min[1] + cell_dim[1],
        cell_min[2] + cell_dim[2],
    };

    const int cell_pbc[3] = {
        MAX((int)ceilf(pbc_ext.x * INV_CELL_EXT), 1),
        MAX((int)ceilf(pbc_ext.y * INV_CELL_EXT), 1),
        MAX((int)ceilf(pbc_ext.z * INV_CELL_EXT), 1),
    };

    // Deperiodize the extent of the search (pos +/- rad) with respect to the periodic domain
    const vec4_t pos_min = vec4_deperiodize(vec4_sub_f(pos, rad), ref, pbc_ext);
    const vec4_t pos_max = vec4_add_f(pos_min, 2.0f * rad);

    const vec4_t cell_pos_min = vec4_floor(vec4_mul_f(pos_min, INV_CELL_EXT));
    const vec4_t cell_pos_max = vec4_ceil (vec4_mul_f(pos_max, INV_CELL_EXT));

    // cell_pos_min is within the correct period, now we can skip forward to the first cell which actually contains data
    // We use the convention of [0,dim[ rather than [min,min+dim[,
    int cell_beg[3] = {
        MAX((int)cell_pos_min.x, cell_min[0]),
        MAX((int)cell_pos_min.y, cell_min[1]),
        MAX((int)cell_pos_min.z, cell_min[2]),
    };

    // Extent in cells to search
    int cell_ext[3] = {
        MAX(MIN((int)cell_pos_max.x - cell_beg[0], cell_pbc[0] - 1) + 1, 1),
        MAX(MIN((int)cell_pos_max.y - cell_beg[1], cell_pbc[1] - 1) + 1, 1),
        MAX(MIN((int)cell_pos_max.z - cell_beg[2], cell_pbc[2] - 1) + 1, 1),
    };

    // How far ahead we need skip in order to end up within the occupied cell domain of the next period
    const int cell_jmp[3] = {
        cell_pbc[0] - cell_dim[0],
        cell_pbc[1] - cell_dim[1],
        cell_pbc[2] - cell_dim[2],
    };
    
    // If cell beg is outside of occupied cell domain [cell_min, cell_dim[ we skip forward to next period and accomodate for that jump
    if (cell_beg[0] >= cell_max[0]) {
        const int jmp = cell_pbc[0] - cell_beg[0] + cell_min[0];
        cell_beg[0] = cell_min[0];
        cell_ext[0] -= jmp;
        if (cell_ext[0] <= 0) return;
    }
    if (cell_beg[1] >= cell_max[1]) {
        const int jmp = cell_pbc[1] - cell_beg[1] + cell_min[1];
        cell_beg[1] = cell_min[1];
        cell_ext[1] -= jmp;
        if (cell_ext[1] <= 0) return;
    }
    if (cell_beg[2] >= cell_max[2]) {
        const int jmp = cell_pbc[2] - cell_beg[2] + cell_min[2];
        cell_beg[2] = cell_min[2];
        cell_ext[2] -= jmp;
        if (cell_ext[2] <= 0) return;
    }

    int cell_end[3] = {
        cell_beg[0] + cell_ext[0],
        cell_beg[1] + cell_ext[1],
        cell_beg[2] + cell_ext[2],
    };
    
    const int cd_0  = cell_dim[0];
    const int cd_01 = cell_dim[0] * cell_dim[1];

    for (int cz = cell_beg[2]; cz < cell_end[2]; inc_cell2(&cz, cell_max[2], cell_jmp[2])) {
        const int cdz = (cz % cell_pbc[2]);
        const int ciz = cdz - cell_min[2];
        ASSERT(0 <= ciz && ciz < cell_dim[2]);
        const int idx_z = ciz * cd_01;
        for (int cy = cell_beg[1]; cy < cell_end[1]; inc_cell2(&cy, cell_max[1], cell_jmp[1])) {
            const int cdy = (cy % cell_pbc[1]);
            const int ciy = cdy - cell_min[1];
            ASSERT(0 <= ciy && ciy < cell_dim[1]);
            const int idx_yz = idx_z + ciy * cd_0;
            for (int cx = cell_beg[0]; cx < cell_end[0]; inc_cell2(&cx, cell_max[0], cell_jmp[0])) {
                const int cdx = (cx % cell_pbc[0]);
                const int cix = cdx - cell_min[0];
                ASSERT(0 <= cix && cix < cell_dim[0]);

                const uint32_t cell_idx = idx_yz + cix;
                const uint32_t cell_data = cells[cell_idx];
                if (cell_data == 0) {
                    continue;
                }

                const int cell_offset = cell_data >> LENGTH_BITS;
                const int cell_length = cell_data & LENGTH_MASK;

                const elem_t* cell_elem = elems + cell_offset;               
                
                for (int i = 0; i < cell_length; ++i) {
                    const elem_t* elem = cell_elem + i;
                    const vec4_t p = vec4_from_vec3(elem->xyz, 0);
                    const float d2 = vec4_periodic_distance_squared(p, pos, pbc_ext);
                    if (d2 < rad2) {
                        if (!iter(elem, user_param)) {
                            return;
                        }
                    }
                }
            }
        }
    }
}

static inline int mini(int x, int y) {
    return x < y ? x : y;
}

static inline void query_pos_rad_periodic_batch(const md_spatial_hash_t* hash, vec3_t position, float rad, md_spatial_hash_batch_iter_fn iter, void* user_param) {
    ASSERT(hash);
    ASSERT(iter);

    const vec4_t ref = vec4_mul_f(hash->pbc_ext, 0.5f);
    const vec4_t pbc_ext = hash->pbc_ext;
    const vec4_t pos = vec4_from_vec3(position, 0);
    float rad2 = rad * rad;

    const int32_t* cell_min = hash->cell_min;
    const int32_t* cell_dim = hash->cell_dim;
    const elem_t* elems = hash->data;
    const cell_t* cells = (const cell_t*)((const char*)hash->data + sizeof(elem_t) * ROUND_UP(hash->elem_count, md_simd_f32_width));

    const int32_t cell_max[3] = {
        cell_min[0] + cell_dim[0],
        cell_min[1] + cell_dim[1],
        cell_min[2] + cell_dim[2],
    };

    const int32_t cell_pbc[3] = {
        MAX((int32_t)ceilf(pbc_ext.x * INV_CELL_EXT), 1),
        MAX((int32_t)ceilf(pbc_ext.y * INV_CELL_EXT), 1),
        MAX((int32_t)ceilf(pbc_ext.z * INV_CELL_EXT), 1),
    };

    // Deperiodize the extent of the search (pos +/- rad) with respect to the periodic domain
    const vec4_t pos_min = vec4_deperiodize(vec4_sub_f(pos, rad), ref, pbc_ext);
    const vec4_t pos_max = vec4_add_f(pos_min, 2.0f * rad);

    const vec4_t cell_pos_min = vec4_floor(vec4_mul_f(pos_min, INV_CELL_EXT));
    const vec4_t cell_pos_max = vec4_ceil (vec4_mul_f(pos_max, INV_CELL_EXT));

    // cell_pos_min is within the correct period, now we can skip forward to the first cell which actually contains data (cell_min)
    int32_t cell_beg[3] = {
        MAX((int32_t)cell_pos_min.x, cell_min[0]),
        MAX((int32_t)cell_pos_min.y, cell_min[1]),
        MAX((int32_t)cell_pos_min.z, cell_min[2]),
    };

    // Extent in cells to search
    int32_t cell_ext[3] = {
        MAX(MIN((int32_t)cell_pos_max.x - cell_beg[0], cell_pbc[0] - 1) + 1, 1),
        MAX(MIN((int32_t)cell_pos_max.y - cell_beg[1], cell_pbc[1] - 1) + 1, 1),
        MAX(MIN((int32_t)cell_pos_max.z - cell_beg[2], cell_pbc[2] - 1) + 1, 1),
    };

    // How far ahead we need skip in order to end up within the occupied cell domain of the next period
    const int32_t cell_jmp[3] = {
        cell_pbc[0] - cell_dim[0],
        cell_pbc[1] - cell_dim[1],
        cell_pbc[2] - cell_dim[2],
    };

    // If cell beg is outside of occupied cell domain [cell_min, cell_dim[ we skip forward to next period and accomodate for that jump
    if (cell_beg[0] >= cell_max[0]) {
        const int32_t jmp = cell_pbc[0] - cell_beg[0] + cell_min[0];
        cell_beg[0] = cell_min[0];
        cell_ext[0] -= jmp;
        if (cell_ext[0] <= 0) return;
    }
    if (cell_beg[1] >= cell_max[1]) {
        const int32_t jmp = cell_pbc[1] - cell_beg[1] + cell_min[1];
        cell_beg[1] = cell_min[1];
        cell_ext[1] -= jmp;
        if (cell_ext[1] <= 0) return;
    }
    if (cell_beg[2] >= cell_max[2]) {
        const int32_t jmp = cell_pbc[2] - cell_beg[2] + cell_min[2];
        cell_beg[2] = cell_min[2];
        cell_ext[2] -= jmp;
        if (cell_ext[2] <= 0) return;
    }

    const int32_t cell_end[3] = {
        cell_beg[0] + cell_ext[0],
        cell_beg[1] + cell_ext[1],
        cell_beg[2] + cell_ext[2],
    };

    const int32_t cd_0  = cell_dim[0];
    const int32_t cd_01 = cell_dim[0] * cell_dim[1];

    const md_simd_f32_t rx = md_simd_set1_f32(pos.x);
    const md_simd_f32_t ry = md_simd_set1_f32(pos.y);
    const md_simd_f32_t rz = md_simd_set1_f32(pos.z);

    const md_simd_f32_t px = md_simd_set1_f32(pbc_ext.x);
    const md_simd_f32_t py = md_simd_set1_f32(pbc_ext.y);
    const md_simd_f32_t pz = md_simd_set1_f32(pbc_ext.z);

    const md_simd_f32_t rpx = pbc_ext.x ? md_simd_div(md_simd_set1_f32(1.0f), px) : md_simd_zero_f32();
    const md_simd_f32_t rpy = pbc_ext.y ? md_simd_div(md_simd_set1_f32(1.0f), py) : md_simd_zero_f32();
    const md_simd_f32_t rpz = pbc_ext.z ? md_simd_div(md_simd_set1_f32(1.0f), pz) : md_simd_zero_f32();

    for (int32_t cz = cell_beg[2]; cz < cell_end[2]; inc_cell2(&cz, cell_max[2], cell_jmp[2])) {
        const int32_t ciz = (cz % cell_pbc[2]) - cell_min[2];
        ASSERT(0 <= ciz && ciz < cell_dim[2]);
        const int32_t idx_z = ciz * cd_01;
        for (int32_t cy = cell_beg[1]; cy < cell_end[1]; inc_cell2(&cy, cell_max[1], cell_jmp[1])) {
            const int32_t ciy = (cy % cell_pbc[1]) - cell_min[1];
            ASSERT(0 <= ciy && ciy < cell_dim[1]);
            const int32_t idx_yz = idx_z + ciy * cd_0;
            for (int32_t cx = cell_beg[0]; cx < cell_end[0]; inc_cell2(&cx, cell_max[0], cell_jmp[0])) {
                const int32_t cix = (cx % cell_pbc[0]) - cell_min[0];
                ASSERT(0 <= cix && cix < cell_dim[0]);

                const uint32_t cell_idx = idx_yz + cix;
                const uint32_t cell_data = cells[cell_idx];
                if (cell_data == 0) {
                    continue;
                }

                const int32_t cell_offset = cell_data >> LENGTH_BITS;
                const int32_t cell_length = cell_data  & LENGTH_MASK;

                int32_t len = cell_length;

                const elem_t* elem = elems + cell_offset;

                while (len > 0) {
                    md_simd_f32_t vx,vy,vz;
                    md_simd_unpack_xyz(&vx, &vy, &vz, (const float*)elem, sizeof(elem_t));

                    md_simd_f32_t dx = md_simd_minimum_image(md_simd_sub(vx, rx), px, rpx);
                    md_simd_f32_t dy = md_simd_minimum_image(md_simd_sub(vy, ry), py, rpy);
                    md_simd_f32_t dz = md_simd_minimum_image(md_simd_sub(vz, rz), pz, rpz);
                    md_simd_f32_t d2 = md_simd_add(md_simd_add(md_simd_mul(dx, dx), md_simd_mul(dy, dy)), md_simd_mul(dz, dz));
                    md_simd_f32_t vmask = md_simd_cmp_lt(d2, md_simd_set1_f32(rad2));

                    const int32_t step = mini(len, md_simd_f32_width);
                    const int32_t lane_mask = (1 << step) - 1;
                    const int32_t mask = md_simd_movemask(vmask) & lane_mask;

                    if (!iter(elem, mask, user_param)) {
                        return;
                    }
                    len  -= md_simd_f32_width;
                    elem += md_simd_f32_width;
                };
            }
        }
    }
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

void md_spatial_hash_query(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_iterator_fn iter, void* user_param) {
    ASSERT(spatial_hash);
    ASSERT(iter);

    if (!validate_spatial_hash(spatial_hash)) {
        return;
    }

    // This is for non periodic lookup
    if (vec4_equal(spatial_hash->pbc_ext, vec4_zero())) {
        query_pos_rad(spatial_hash, pos, radius, iter, user_param);
    } else {
        query_pos_rad_periodic(spatial_hash, pos, radius, iter, user_param);
    }
}

void md_spatial_hash_query_batch(const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius, md_spatial_hash_batch_iter_fn iter, void* user_param) {
    ASSERT(spatial_hash);
    ASSERT(iter);

    if (!validate_spatial_hash(spatial_hash)) {
        return;
    }

    // This is for non periodic lookup
    if (vec4_equal(spatial_hash->pbc_ext, vec4_zero())) {
        query_pos_rad_batch(spatial_hash, pos, radius, iter, user_param);
    } else {
        query_pos_rad_periodic_batch(spatial_hash, pos, radius, iter, user_param);
    }
}

bool idx_fn(const md_spatial_hash_elem_t* elem_arr, int mask, void* user_param) {
    md_bitfield_t* bf = (md_bitfield_t*)user_param;
    while (mask) {
        const int idx = ctz32(mask);
        md_bitfield_set_bit(bf, elem_arr[idx].idx);
        mask = mask & ~(1 << idx);
    }
    return true;
}

int64_t md_spatial_hash_query_idx(int32_t* buf, int64_t cap, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius) {
    md_bitfield_t bf = md_bitfield_create(md_temp_allocator);
    md_bitfield_reserve_range(&bf, 0, spatial_hash->elem_count);

    md_spatial_hash_query_batch(spatial_hash, pos, radius, idx_fn, &bf);
#if DEBUG
    if ((int64_t)md_bitfield_popcount(&bf) > cap) {
        MD_LOG_DEBUG("Size exceeds the capacity, some elements will not be written");
    }
#endif

    return md_bitfield_extract_indices(buf, cap, &bf);
}

void md_spatial_hash_query_bits(md_bitfield_t* bf, const md_spatial_hash_t* spatial_hash, vec3_t pos, float radius) {
    md_spatial_hash_query_batch(spatial_hash, pos, radius, idx_fn, bf);
}
