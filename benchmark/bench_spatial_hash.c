#include "ubench.h"

#include <md_gro.h>
#include <md_molecule.h>
#include <core/md_intrinsics.h>
#include <core/md_spatial_hash.h>
#include <core/md_arena_allocator.h>

struct spatial_hash {
    md_allocator_i* arena;
    md_system_t mol;
    md_spatial_hash_t* pbc_sh;
    vec3_t pbc_ext;

    vec3_t* coords;
    
    md_spatial_hash_t* reg_sh;
};

UBENCH_F_SETUP(spatial_hash) {
    ubench_fixture->arena = md_vm_arena_create(GIGABYTES(4));
    md_allocator_i* alloc = ubench_fixture->arena;

    md_gro_data_t gro_data = {0};
    md_gro_data_parse_file(&gro_data, STR_LIT(MD_BENCHMARK_DATA_DIR "/centered.gro"), alloc);
    md_gro_molecule_init(&ubench_fixture->mol, &gro_data, alloc);
    ubench_fixture->pbc_sh = md_spatial_hash_create_soa(ubench_fixture->mol.atom.x, ubench_fixture->mol.atom.y, ubench_fixture->mol.atom.z, NULL, ubench_fixture->mol.atom.count, &ubench_fixture->mol.unit_cell, alloc);
    ubench_fixture->reg_sh = md_spatial_hash_create_soa(ubench_fixture->mol.atom.x, ubench_fixture->mol.atom.y, ubench_fixture->mol.atom.z, NULL, ubench_fixture->mol.atom.count, NULL, alloc);

    md_system_t* mol = &ubench_fixture->mol;

    ubench_fixture->coords = md_vm_arena_push(alloc, mol->atom.count * sizeof(vec3_t));
    for (size_t i = 0; i < mol->atom.count; ++i) {
        ubench_fixture->coords[i] = vec3_set(mol->atom.x[i], mol->atom.y[i], mol->atom.z[i]);
    }
}

UBENCH_F_TEARDOWN(spatial_hash) {
    md_vm_arena_destroy(ubench_fixture->arena);
}

UBENCH_EX_F(spatial_hash, init_non_periodic) {
    md_system_t* mol = &ubench_fixture->mol;
    md_allocator_i* alloc = ubench_fixture->arena;

    UBENCH_DO_BENCHMARK() {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(alloc);
        UBENCH_DO_NOTHING(md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, NULL, mol->atom.count, NULL, alloc));
        md_vm_arena_temp_end(temp);
    }
}

UBENCH_EX_F(spatial_hash, init_periodic) {
    md_system_t* mol = &ubench_fixture->mol;
    md_allocator_i* alloc = ubench_fixture->arena;

    UBENCH_DO_BENCHMARK() {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(alloc);
        UBENCH_DO_NOTHING(md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, NULL, mol->atom.count, &mol->unit_cell, alloc));
        md_vm_arena_temp_end(temp);
    }
}

static bool func(const md_spatial_hash_elem_t* elem, void* user_data) {
    (void)elem;
    uint32_t* count = user_data;
    (*count)++;
    return true;
}

#define COUNT 1000
#define RAD 10.0f

UBENCH_EX_F(spatial_hash, query_non_periodic) {
    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        md_spatial_hash_query_multi(ubench_fixture->reg_sh, ubench_fixture->coords, COUNT, RAD, func, &count);
        UBENCH_DO_NOTHING(&count);
    }
}

UBENCH_EX_F(spatial_hash, query_periodic) {
    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        md_spatial_hash_query_multi(ubench_fixture->pbc_sh, ubench_fixture->coords, COUNT, RAD, func, &count);
        UBENCH_DO_NOTHING(&count);
    }
}

static inline bool batch_func(const md_spatial_hash_elem_t* elem, int mask, void* user_data) {
    (void)elem;
    uint32_t* count = user_data;
    *count += popcnt32(mask);
    return true;
}

UBENCH_EX_F(spatial_hash, query_non_periodic_batch) {
    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        md_spatial_hash_query_multi_batch(ubench_fixture->reg_sh, ubench_fixture->coords, COUNT, RAD, batch_func, &count);
        UBENCH_DO_NOTHING(&count);
    }
}

UBENCH_EX_F(spatial_hash, query_periodic_batch) {
    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        md_spatial_hash_query_multi_batch(ubench_fixture->pbc_sh, ubench_fixture->coords, COUNT, RAD, batch_func, &count);
        UBENCH_DO_NOTHING(&count);
    }
}

UBENCH_EX_F(spatial_hash, query_bruteforce) {
    const md_system_t* mol = &ubench_fixture->mol;

    uint32_t count = 0;
    const float RAD2 = RAD * RAD;
    UBENCH_DO_BENCHMARK() {
        const vec4_t pbc_ext = vec4_from_vec3(mat3_diag(mol->unit_cell.basis), 0);
        for (size_t i = 0; i < COUNT; ++i) {
            const vec4_t p = vec4_set(mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0);
            for (size_t j = 0; j < mol->atom.count; ++j) {
                const vec4_t c = vec4_set(mol->atom.x[j], mol->atom.y[j], mol->atom.z[j], 0);
                if (vec4_periodic_distance_squared(p, c, pbc_ext) < RAD2) {
                    count += 1;
                }
            }
        }
        UBENCH_DO_NOTHING(&count);
    }
}

/*
static inline void do_it(float* out_xyzw, const float* in_x, const float* in_y, const float* in_z, const float in_w, size_t count) {
    __m128 wwww = _mm_set1_ps(in_w);
    for (size_t i = 0; i < count; i += 4) {
        __m128 xxxx = _mm_load_ps(in_x + i);
        __m128 yyyy = _mm_load_ps(in_y + i);
        __m128 zzzz = _mm_load_ps(in_z + i);

        __m128 _Tmp0 = _mm_shuffle_ps((xxxx), (yyyy), 0x44);
        __m128 _Tmp2 = _mm_shuffle_ps((xxxx), (yyyy), 0xEE);
        __m128 _Tmp1 = _mm_shuffle_ps((zzzz), (wwww), 0x44);
        __m128 _Tmp3 = _mm_shuffle_ps((zzzz), (wwww), 0xEE);

        __m128 xyzw0 = _mm_shuffle_ps(_Tmp0, _Tmp1, 0x88);
        __m128 xyzw1 = _mm_shuffle_ps(_Tmp0, _Tmp1, 0xDD);
        __m128 xyzw2 = _mm_shuffle_ps(_Tmp2, _Tmp3, 0x88);
        __m128 xyzw3 = _mm_shuffle_ps(_Tmp2, _Tmp3, 0xDD);

        _mm_store_ps(out_xyzw + i * 4 + 0, xyzw0);
        _mm_store_ps(out_xyzw + i * 4 + 1, xyzw1);
        _mm_store_ps(out_xyzw + i * 4 + 2, xyzw2);
        _mm_store_ps(out_xyzw + i * 4 + 3, xyzw3);
    }
}

static inline void do_it2(vec4_t* out_xyzw, const float* in_x, const float* in_y, const float* in_z, const float in_w, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        out_xyzw[i] = vec4_set(in_x[i], in_y[i], in_z[i], in_w);
    }
}

UBENCH_EX_F(spatial_hash, testink_1) {
    md_system_t* mol = &ubench_fixture->mol;
    md_allocator_i* alloc = &ubench_fixture->alloc;
    vec3_t pbc_ext = ubench_fixture->pbc_ext;

    md_spatial_hash_t spatial_hash = {0};

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&ubench_fixture->arena);
    vec4_t* data = md_vm_arena_push(&ubench_fixture->arena, sizeof(vec4_t) * mol->atom.count);
    UBENCH_DO_BENCHMARK() {
        do_it(data, mol->atom.x, mol->atom.y, mol->atom.z, 1.0f, mol->atom.count);
    }

    md_vm_arena_temp_end(temp);
}

UBENCH_EX_F(spatial_hash, testink_2) {
    md_system_t* mol = &ubench_fixture->mol;
    md_allocator_i* alloc = &ubench_fixture->alloc;
    vec3_t pbc_ext = ubench_fixture->pbc_ext;

    md_spatial_hash_t spatial_hash = {0};

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&ubench_fixture->arena);
    vec4_t* data = md_vm_arena_push(&ubench_fixture->arena, sizeof(vec4_t) * mol->atom.count);
    UBENCH_DO_BENCHMARK() {
        do_it2(data, mol->atom.x, mol->atom.y, mol->atom.z, 1.0f, mol->atom.count);
    }

    md_vm_arena_temp_end(temp);
}
*/
