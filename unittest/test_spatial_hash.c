#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_spatial_hash.h>
#include <core/md_str.h>
#include <core/md_intrinsics.h>
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

    md_unit_cell_t unit_cell = md_util_unit_cell_from_extent(10, 0, 0);
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(x, y, z, NULL, 10, &unit_cell, md_get_heap_allocator());
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
    md_spatial_hash_query(spatial_hash, pos, rad, func, &count);

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

static inline size_t test_elem(md_256 rx, md_256 ry, md_256 rz, md_256 r2, const elem_t* elem, uint32_t len) {
    size_t result = 0;

    for (uint32_t i = 0; i < len; i += 8) {
        md_256 vx, vy, vz;
        md_mm256_unpack_xyz_ps(&vx, &vy, &vz, (const float*)(elem + i), sizeof(elem_t));
        md_256 dx = md_mm256_sub_ps(vx, rx);
        md_256 dy = md_mm256_sub_ps(vy, ry);
        md_256 dz = md_mm256_sub_ps(vz, rz);
        md_256 d2 = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(dx, dx), md_mm256_mul_ps(dy, dy)), md_mm256_mul_ps(dz, dz));
        md_256 vmask = md_mm256_cmplt_ps(d2, r2);

        const uint32_t remainder = MIN(len - i, 8);
        const uint32_t lane_mask = (1U << remainder) - 1U;
        const uint32_t mask = md_mm256_movemask_ps(vmask) & lane_mask;

        result += popcnt32(mask);
    }

    return result;
}

static inline size_t test_cell(const data_t* data, uint32_t cell_idx) {
    uint32_t offset = data->cell_offset[cell_idx];
    uint32_t length = data->cell_offset[cell_idx + 1] - offset;

    return test_elem(data->rx, data->ry, data->rz, data->r2, data->elem + offset, length);
}

#define CELL_EXT (6.0f)

static size_t do_pairwise_periodic(const float* in_x, const float* in_y, const float* in_z, size_t num_points, float cutoff, const md_unit_cell_t* unit_cell) {
    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));

    vec3_t ext = md_unit_cell_extent(unit_cell);

    printf("unit cell ext: %f %f %f \n", ext.x, ext.y, ext.z);

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

    vec4_t inv_cell_ext = vec4_div(vec4_set1(1.0f), cell_ext);

    uint32_t c0  = cell_dim[0];
    uint32_t c01 = cell_dim[0] * cell_dim[1];

    size_t num_cells = cell_dim[0] * cell_dim[1] * cell_dim[2];
    size_t cell_offset_count = num_cells + 1;

    printf("Cell dim: %i %i %i \n", cell_dim[0], cell_dim[1], cell_dim[2]);

    uint32_t* cell_offset = md_vm_arena_push_zero(temp_arena, cell_offset_count * sizeof(uint32_t));
    uint32_t* local_idx   = md_vm_arena_push(temp_arena, num_points * sizeof(uint32_t));
    uint32_t* cell_idx    = md_vm_arena_push(temp_arena, num_points * sizeof(uint32_t));
    elem_t* elements      = md_vm_arena_push(temp_arena, num_points * sizeof(elem_t));

    const vec4_t mask = md_unit_cell_pbc_mask(unit_cell);
    mat4x3_t I = mat4x3_from_mat3(unit_cell->inv_basis);
    mat4x3_t M = mat4x3_from_mat3(unit_cell->basis);

    // Calculate cell and local cell indices
    for (size_t i = 0; i < num_points; ++i) {
        vec4_t coord = {in_x[i], in_y[i], in_z[i], 0};
        coord = vec4_mul(coord, inv_cell_ext);

        uint32_t cx = (uint32_t)coord.x;
        uint32_t cy = (uint32_t)coord.y;
        uint32_t cz = (uint32_t)coord.z;
        uint32_t ci = cz * c01 + cy * c0 + cx;
        if (ci >= num_cells) {
            printf("Invalid index with cell coords: %i %i %i\n", cx, cy, cz);
        }
        ASSERT(ci < num_cells);
        uint32_t li = cell_offset[ci]++;
        local_idx[i] = li;
        cell_idx[i]  = ci;
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
            i,
        };
    }

    data_t data = { 0 };
    data.r2 = md_mm256_set1_ps(cutoff * cutoff);
    data.elem = elements;
    data.cell_offset = cell_offset;

    size_t result = 0;

    for (uint32_t ci = 0; ci < num_cells - 1; ++ci) {
        uint32_t offset = cell_offset[ci];
        uint32_t length = cell_offset[ci + 1] - offset;
        const elem_t* elem = elements + offset;
        for (size_t i = 0; i < length; ++i) {
            data.rx = md_mm256_set1_ps(elem[i].x);
            data.ry = md_mm256_set1_ps(elem[i].y);
            data.rz = md_mm256_set1_ps(elem[i].z);
            result += test_elem(data.rx, data.ry, data.rz, data.r2, elem + i, length - i);

            for (uint32_t cj = ci + 1; cj < num_cells; ++cj) {
                result += test_cell(&data, cj);
            }
        }
    }

    md_vm_arena_destroy(temp_arena);

    return result;
}

UTEST(spatial_hash, n2_custom) {
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

    size_t count = do_pairwise_periodic(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &mol.unit_cell);

    EXPECT_EQ(30, count);

    md_arena_allocator_destroy(alloc);
}

UTEST(spatial_hash, n2_brute_force) {
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

    const float r2 = 5.0f * 5.0f;

    size_t count = 0;
    for (size_t i = 0; i < mol.atom.count - 1; ++i) {
        vec4_t xi = {mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], 0.0f};
        for (size_t j = i + 1; j < mol.atom.count; ++j) {
            vec4_t xj = {mol.atom.x[j], mol.atom.y[j], mol.atom.z[j], 0.0f};
            if (vec4_distance_squared(xi, xj) < r2) {
                count += 1;
            }
        }
    }

    EXPECT_EQ(30, count);

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