#include "ubench.h"

#include <md_gro.h>
#include <md_molecule.h>
#include <core/md_spatial_hash.h>
#include <core/md_arena_allocator.h>

struct spatial_hash {
    md_vm_arena_t arena;
    md_allocator_i alloc;
    md_molecule_t mol;
    md_spatial_hash_t* sh;
    vec3_t pbc_ext;
};

UBENCH_F_SETUP(spatial_hash) {
    md_vm_arena_init(&ubench_fixture->arena, GIGABYTES(4));

    ubench_fixture->alloc = md_vm_arena_create_interface(&ubench_fixture->arena);

    md_gro_data_t gro_data = {0};
    md_gro_data_parse_file(&gro_data, STR(MD_BENCHMARK_DATA_DIR "/centered.gro"), &ubench_fixture->alloc);
    md_gro_molecule_init(&ubench_fixture->mol, &gro_data, &ubench_fixture->alloc);
    ubench_fixture->sh = md_spatial_hash_create_soa(ubench_fixture->mol.atom.x, ubench_fixture->mol.atom.y, ubench_fixture->mol.atom.z, NULL, ubench_fixture->mol.atom.count, &ubench_fixture->mol.unit_cell, &ubench_fixture->alloc);
}

UBENCH_F_TEARDOWN(spatial_hash) {
    md_vm_arena_free(&ubench_fixture->arena);
}

UBENCH_EX_F(spatial_hash, init) {
    md_molecule_t* mol = &ubench_fixture->mol;
    md_allocator_i* alloc = &ubench_fixture->alloc;
    vec3_t pbc_ext = ubench_fixture->pbc_ext;


    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&ubench_fixture->arena);
    UBENCH_DO_BENCHMARK() {
        md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, NULL, mol->atom.count, &mol->unit_cell, alloc);
    }

    md_vm_arena_temp_end(temp);
}

static bool func(uint32_t idx, vec3_t pos, void* param) {
    (void)idx;
    (void)pos;
    uint32_t* count = param;
    (*count)++;
    return true;
}

UBENCH_EX_F(spatial_hash, query) {
    const vec3_t pos[10] = {
        {24, 8,  4},
        {14, 1,  1},
        {44, 57, 30},
        {24, 48, 5},
        {24, 32, 8},
        {14, 10, 34},
        {04, 19, 44},
        {24, 0,  14},
        {34, 100,24},
        {23, 29, 1}
    };
    const float rad = 10.0f;
    
    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        for (int i = 0; i < 10000; ++i) {
            md_spatial_hash_query(ubench_fixture->sh, pos[i % 10], rad, func, &count);
        }
    }
}

UBENCH_EX_F(spatial_hash, bruteforce_query) {
    const vec3_t pos[10] = {
        {24, 8,  4},
        {14, 1,  1},
        {44, 57, 30},
        {24, 48, 5},
        {24, 32, 8},
        {14, 10, 34},
        {04, 19, 44},
        {24, 0,  14},
        {34, 100,24},
        {23, 29, 1}
    };
    const float rad = 10.0f;

    const md_molecule_t* mol = &ubench_fixture->mol;

    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        for (int iter = 0; iter < 10000; ++iter) {
            const vec4_t p = vec4_from_vec3(pos[iter % 10], 0);
            for (int64_t i = 0; i < mol->atom.count; ++i) {
                const vec4_t c = vec4_set(mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0);
                if (vec4_distance_squared(p, c) < rad * rad) {
                    count += 1;
                    UBENCH_DO_NOTHING(&count);
                }
            }
        }
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
    md_molecule_t* mol = &ubench_fixture->mol;
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
    md_molecule_t* mol = &ubench_fixture->mol;
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