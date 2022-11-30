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

static bool func(uint32_t idx, vec3_t pos, void* param) {
    (void)idx;
    (void)pos;
    uint32_t* count = param;
    (*count)++;
    return true;
}

UTEST(spatial_hash, small_periodic) {
    float x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    float y[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    vec3_t pbc_ext = {10,0,0};
    md_spatial_hash_t spatial_hash = {0};
    ASSERT_TRUE(md_spatial_hash_init_soa(&spatial_hash, x, y, z, 10, pbc_ext, default_allocator));

    vec3_t pos = {5, 0, 0};
    float  rad = 1.5f;
    uint32_t count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, pos, rad, func, &count));
    EXPECT_EQ(3, count);

    count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, (vec3_t){8.5f, 0, 0}, 3, func, &count));
    EXPECT_EQ(6, count);

    md_spatial_hash_free(&spatial_hash);
}

UTEST(spatial_hash, big) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(64));
    const str_t pdb_file = MAKE_STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    vec3_t pbc_ext = {0,0,0};
    md_spatial_hash_t spatial_hash = {0};
    ASSERT_TRUE(md_spatial_hash_init_soa(&spatial_hash, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, pbc_ext, alloc));

    vec3_t pos = {24, 48, 24};
    float  rad = 10.0f;
    uint32_t count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, pos, rad, func, &count));

    EXPECT_LT(0, count);

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
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, MAKE_STR(MD_UNITTEST_DATA_DIR "/centered.gro"), &utest_fixture->alloc));
    ASSERT_TRUE(md_gro_molecule_init(&utest_fixture->mol, &gro_data, &utest_fixture->alloc));
    utest_fixture->pbc_ext = mat3_mul_vec3(utest_fixture->mol.coord_frame, vec3_set1(1.f));
}

UTEST_F_TEARDOWN(spatial_hash) {
    md_vm_arena_free(&utest_fixture->arena);
}

UTEST_F(spatial_hash, perf_test_init) {
    md_molecule_t* mol = &utest_fixture->mol;
    md_allocator_i* alloc = &utest_fixture->alloc;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t spatial_hash = {0};

    const uint32_t num_iter = 100;
    timestamp_t t0 = md_os_time_current();
    for (uint32_t i = 0; i < num_iter; ++i) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);
        md_spatial_hash_init_soa(&spatial_hash, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.count, pbc_ext, alloc);
        md_vm_arena_temp_end(temp);
    }
    timestamp_t t1 = md_os_time_current();
    double t = md_os_time_as_milliseconds(t1 - t0);

    printf("Avg. time taken to generate spatial hash structure of %i atoms: %.3fms\n", (int)mol->atom.count, t / (double)num_iter);
}

static inline float rnd() {
    return rand() / (float)RAND_MAX;
}

typedef struct param_t {
    uint32_t *count;
    md_bitfield_t* ref_bf;
} param_t;

static bool iter_fn(uint32_t idx, vec3_t coord, void* user_param) {
    param_t* param = user_param; 
    uint32_t* count = param->count;
    *count += 1;
    return true;
}

UTEST_F(spatial_hash, test_correctness_non_periodic) {
    md_molecule_t* mol = &utest_fixture->mol;
    md_allocator_i* alloc = &utest_fixture->alloc;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t spatial_hash = {0};
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);
    md_spatial_hash_init_soa(&spatial_hash, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.count, (vec3_t){0,0,0}, alloc);
    
    srand(31);

    /*
    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);
    md_bitfield_reserve_range(&bf, 0, mol->atom.count);
    */

    timestamp_t t = 0;
    const uint32_t num_iter = 100;
    for (uint32_t iter = 0; iter < 100; ++iter) {
        //md_bitfield_clear(&bf);
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float radius = rnd() * 30;

        uint32_t ref_count = 0;
        const float rad2 = radius * radius;
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            vec4_t c = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};

            if (vec4_distance_squared(vec4_from_vec3(pos, 0), c) < rad2) {
                ref_count += 1;
                //md_bitfield_set_bit(&bf, i);
            }
        }

        uint32_t count = 0;
        param_t param = {
            .count = &count,
            //.ref_bf = &bf,
        };
        timestamp_t t0 = md_os_time_current();
        md_spatial_hash_query(&spatial_hash, pos, radius, iter_fn, &param);
        timestamp_t t1 = md_os_time_current();
        t += (t1 - t0);

        // We have quantization artifacts since the coordinates are compressed into 10-bits per dimension relative to the cell.
        // This means we will have some straddling cases that are wither just within or outside of the search radius, thus we cannot match the reference exactly.
        int delta = (int)ref_count - (int)count;
        EXPECT_LT(ABS(delta), 2);
    }

    printf("Avg. time taken per query: %.4fms\n", md_os_time_as_milliseconds(t) / (double)num_iter);

    md_vm_arena_temp_end(temp);
}

UTEST_F(spatial_hash, test_correctness_periodic) {
    md_molecule_t* mol = &utest_fixture->mol;
    md_allocator_i* alloc = &utest_fixture->alloc;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t spatial_hash = {0};
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&utest_fixture->arena);
    md_spatial_hash_init_soa(&spatial_hash, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.count, pbc_ext, alloc);

    srand(31);

    /*
    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);
    md_bitfield_reserve_range(&bf, 0, mol->atom.count);
    */

    timestamp_t t = 0;
    const uint32_t num_iter = 100;
    const vec4_t period = vec4_from_vec3(pbc_ext, 0);
    for (uint32_t iter = 0; iter < num_iter; ++iter) {
        //md_bitfield_clear(&bf);
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float radius = rnd() * 30;

        uint32_t ref_count = 0;
        const float rad2 = radius * radius;
        for (int64_t i = 0; i < mol->atom.count; ++i) {
            vec4_t c = {mol->atom.x[i], mol->atom.y[i], mol->atom.z[i], 0};

            if (vec4_periodic_distance_squared(vec4_from_vec3(pos, 0), c, period) < rad2) {
                ref_count += 1;
                //md_bitfield_set_bit(&bf, i);
            }
        }

        uint32_t count = 0;
        param_t param = {
            .count = &count,
            //.ref_bf = &bf,
        };
        timestamp_t t0 = md_os_time_current();
        md_spatial_hash_query(&spatial_hash, pos, radius, iter_fn, &param);
        timestamp_t t1 = md_os_time_current();
        t += (t1 - t0);

        // We have quantization artifacts since the coordinates are compressed into 10-bits per dimension relative to the cell.
        // This means we will have some straddling cases that are wither just within or outside of the search radius, thus we cannot match the reference exactly.
        int delta = (int)ref_count - (int)count;
        EXPECT_LT(ABS(delta), 2);
    }

    printf("Avg. time taken per query: %.4fms\n", md_os_time_as_milliseconds(t) / (double)num_iter);

    md_vm_arena_temp_end(temp);
}