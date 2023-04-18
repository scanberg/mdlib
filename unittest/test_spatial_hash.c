#include "utest.h"

#include "core/md_allocator.h"
#include "core/md_arena_allocator.h"
#include "core/md_bitfield.h"
#include "core/md_spatial_hash.h"
#include "core/md_str.h"
#include "core/md_os.h"
#include <md_pdb.h>
#include <md_gro.h>
#include <md_molecule.h>
#include <md_util.h>

static bool func(const md_spatial_hash_elem_t* elem, void* param) {
    (void)elem;
    uint32_t* count = param;
    (*count)++;
    return true;
}

UTEST(spatial_hash, small_periodic) {
    float x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    float y[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    md_unit_cell_t unit_cell = md_util_unit_cell_ortho(10, 0, 0);
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(x, y, z, NULL, 10, &unit_cell, default_allocator);
    ASSERT_TRUE(spatial_hash);
    
    uint32_t count = 0;
    md_spatial_hash_query(spatial_hash, (vec3_t){5,0,0}, 1.5f, func, &count);
    EXPECT_EQ(3, count);

    count = 0;
    md_spatial_hash_query(spatial_hash, (vec3_t){8.5f, 0, 0}, 3, func, &count);
    EXPECT_EQ(6, count);

    md_spatial_hash_free(spatial_hash);
}

UTEST(spatial_hash, big) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(64));
    const str_t pdb_file = STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &mol.unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);

    vec3_t pos = {24, 48, 24};
    float  rad = 10.0f;
    uint32_t count = 0;
    md_spatial_hash_query(spatial_hash, pos, rad, func, &count);

    EXPECT_NE(0, count);

    md_arena_allocator_destroy(alloc);
}

struct spatial_hash {
    md_vm_arena_t arena;
    md_allocator_i alloc;
    md_molecule_t mol;
    vec3_t pbc_ext;
};

UTEST_F_SETUP(spatial_hash) {
    md_vm_arena_init(&utest_fixture->arena, GIGABYTES(4));

    utest_fixture->alloc = md_vm_arena_create_interface(&utest_fixture->arena);

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR(MD_UNITTEST_DATA_DIR "/centered.gro"), &utest_fixture->alloc));
    ASSERT_TRUE(md_gro_molecule_init(&utest_fixture->mol, &gro_data, &utest_fixture->alloc));
    utest_fixture->pbc_ext = mat3_mul_vec3(utest_fixture->mol.unit_cell.basis, vec3_set1(1.f));
}

UTEST_F_TEARDOWN(spatial_hash) {
    md_vm_arena_free(&utest_fixture->arena);
}

static inline float rnd() {
    return rand() / (float)RAND_MAX;
}

static bool iter_fn(const md_spatial_hash_elem_t* elem, void* user_param) {
    uint32_t* count = (uint32_t*)user_param;
    *count += 1;
    return true;
}

static bool iter_batch_fn(const md_spatial_hash_elem_t* elem, int mask, void* user_param) {
    (void)elem;
    uint32_t* count = (uint32_t*)user_param;
    *count += popcnt32(mask);
    return true;
}

UTEST_F(spatial_hash, test_correctness_centered) {
    md_molecule_t* mol = &utest_fixture->mol;
    md_allocator_i* alloc = &utest_fixture->alloc;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol->atom.x, mol->atom.y, mol->atom.z, NULL, mol->atom.count, NULL, alloc);
    ASSERT_TRUE(spatial_hash);
    
    srand(31);
    
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float rad = rnd() * 20;

        int ref_count = 0;
        const float rad2 = rad * rad;
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            vec4_t c = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};

            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        
        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, rad, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, rad, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);
        if (ref_count != batch_count) {
            int new_count = 0;
            md_spatial_hash_query_batch(spatial_hash, pos, rad, iter_batch_fn, &new_count);
        }
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_ala) {
    md_allocator_i* alloc = &utest_fixture->alloc;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), alloc));

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
        for (int64_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_ala_vec3) {
    md_allocator_i* alloc = &utest_fixture->alloc;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), alloc));

    vec3_t* xyz = md_array_create(vec3_t, mol.atom.count, alloc);
    for (int i = 0; i < mol.atom.count; ++i) {
        xyz[i] = vec3_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]);
    }

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_vec3(xyz, NULL, md_array_size(xyz), &mol.unit_cell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = vec4_from_vec3(mat3_mul_vec3(mol.unit_cell.basis, vec3_set1(1)), 0);

    //srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (int64_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_water) {
    md_allocator_i* alloc = &utest_fixture->alloc;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/water.gro"), alloc));

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
        for (int64_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_periodic_centered) {
    md_allocator_i* alloc = &utest_fixture->alloc;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);

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
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            vec4_t c = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};

            if (vec4_periodic_distance_squared(vec4_from_vec3(pos, 0), c, period) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);
    }

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_periodic_water) {
    md_allocator_i* alloc = &utest_fixture->alloc;
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);

    md_molecule_t mol;
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/water.gro"), alloc));

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
        for (int64_t i = 0; i < mol.atom.count; ++i) {
            const vec4_t c = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0};
            if (vec4_periodic_distance_squared(vec4_from_vec3(pos, 0), c, pbc_ext) < rad2) {
                ref_count += 1;
            }
        }

        int count = 0;
        md_spatial_hash_query(spatial_hash, pos, radius, iter_fn, &count);
        EXPECT_EQ(ref_count, count);

        int batch_count = 0;
        md_spatial_hash_query_batch(spatial_hash, pos, radius, iter_batch_fn, &batch_count);
        EXPECT_EQ(ref_count, batch_count);
    }

    md_vm_arena_temp_end(temp);
}