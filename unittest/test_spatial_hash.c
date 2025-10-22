#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_spatial_hash.h>
#include <core/md_spatial_acc.h>
#include <core/md_str.h>
#include <core/md_intrinsics.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_system.h>
#include <md_util.h>

typedef struct {
    uint32_t i, j;
    double d2;
} dist_pair_t;

typedef struct {
    dist_pair_t* pairs;
    uint32_t count;
} spatial_acc_data_t;

static void spatial_acc_neighbor_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    uint32_t* count = (uint32_t*)user_param;

    const md_256i v_len = md_mm256_set1_epi32((int)num_pairs);
    const md_256i v_8 = md_mm256_set1_epi32(8);
    const md_256 v_r2 = md_mm256_set1_ps(25.0f);  // 5.0^2
    md_256i v_k = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);

    for (size_t k = 0; k < num_pairs; k += 8) {
        md_256i k_mask = md_mm256_cmplt_epi32(v_k, v_len);
        md_256 v_d2 = md_mm256_loadu_ps(ij_dist2 + k);
        md_256 d2_mask = md_mm256_cmplt_ps(v_d2, v_r2);
        md_256 v_mask = md_mm256_and_ps(md_mm256_castsi256_ps(k_mask), d2_mask);

        *count += popcnt32(md_mm256_movemask_ps(v_mask));
        v_k = md_mm256_add_epi32(v_k, v_8);
    }
}

static void spatial_acc_cutoff_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    uint32_t* count = (uint32_t*)user_param;
    *count += (uint32_t)num_pairs;
}

static bool iter_fn(const md_spatial_hash_elem_t* elem, void* user_param) {
    (void)elem;
    uint32_t* count = (uint32_t*)user_param;
    *count += 1;
    return true;
}

typedef struct {
    uint32_t exclude_idx;
    uint32_t count;
} iter_excl_data_t;

static bool iter_excl_fn(const md_spatial_hash_elem_t* elem, void* user_param) {
    (void)elem;
    iter_excl_data_t* data = user_param;
    if (elem->idx != data->exclude_idx) {
        data->count += 1;
    }
    return true;
}

static bool iter_batch_fn(const md_spatial_hash_elem_t* elem, int mask, void* user_param) {
    (void)elem;
    uint32_t* count = (uint32_t*)user_param;
    *count += popcnt32(mask);
    return true;
}

typedef struct {
    dist_pair_t* pairs;
    uint32_t count;
} iter_batch_pair_data_t;

static bool iter_batch_pairs_fn(const md_spatial_hash_elem_t* elem, md_256 d2, int mask, size_t i, void* user_param) {
    iter_batch_pair_data_t* data = user_param;
    
    float d2_raw[8];
    md_mm256_storeu_ps(d2_raw, d2);

    while (mask) {
        const int idx = ctz32(mask);
        mask = mask & ~(1 << idx);
        uint32_t j = elem[idx].idx;
        data->pairs[data->count].i = MIN((uint32_t)i,j);
        data->pairs[data->count].j = MAX((uint32_t)i,j);
        data->pairs[data->count].d2 = d2_raw[idx];
        data->count += 1;
    }

    return true;
}

UTEST(spatial_hash, small_periodic) {
    float x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    float y[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    md_unitcell_t unit_cell = md_unitcell_from_extent(10, 0, 0);
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(x, y, z, NULL, 10, &unit_cell, md_get_heap_allocator());
    ASSERT_TRUE(spatial_hash);
    
    uint32_t count = 0;
    md_spatial_hash_query(spatial_hash, (vec3_t){5,0,0}, 1.5f, iter_fn, &count);
    EXPECT_EQ(3, count);

    count = 0;
    md_spatial_hash_query(spatial_hash, (vec3_t){8.5f, 0, 0}, 3, iter_fn, &count);
    EXPECT_EQ(6, count);

    md_spatial_hash_free(spatial_hash);
}

UTEST(spatial_hash, big) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    const str_t pdb_file = STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_system_t mol = {0};
    ASSERT_TRUE(md_pdb_system_init(&mol, &pdb_data, MD_PDB_OPTION_NONE, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(mol.atom.x, mol.atom.y, mol.atom.z, NULL, mol.atom.count, &mol.unitcell, alloc);
    ASSERT_TRUE(spatial_hash);

    vec3_t pos = {24, 48, 24};
    float  rad = 10.0f;
    uint32_t count = 0;
    md_spatial_hash_query(spatial_hash, pos, rad, iter_fn, &count);

    EXPECT_NE(0, count);

    md_arena_allocator_destroy(alloc);
}

static size_t do_brute_force(const float* in_x, const float* in_y, const float* in_z, size_t num_points, float cutoff, const md_unitcell_t* cell, dist_pair_t* pairs) {
    size_t count = 0;
    const float r2 = cutoff * cutoff;

    uint32_t flags = md_unitcell_flags(cell);

    if ((flags & MD_UNITCELL_PBC_ALL) == MD_UNITCELL_PBC_ALL) {
        vec4_t ext = md_unitcell_diag_vec4(cell);
        vec4_t inv_ext = vec4_set(ext.x > 0.0f ? 1.0f / ext.x : 0.0f, ext.y > 0.0f ? 1.0f / ext.y : 0.0f, ext.z > 0.0f ? 1.0f / ext.z : 0.0f, 0.0f);

        for (size_t i = 0; i < num_points - 1; ++i) {
            float x = in_x[i];
            float y = in_y[i];
            float z = in_z[i];
            for (size_t j = i + 1; j < num_points; ++j) {
                float dx = in_x[j] - x;
                float dy = in_y[j] - y;
                float dz = in_z[j] - z;

                vec4_t d = vec4_min_image(vec4_set(dx, dy, dz, 0.0f), ext, inv_ext);
                float d2 = fmaf(d.x, d.x, fmaf(d.y, d.y, d.z * d.z));
                if (d2 <= r2) {
                    if (pairs) {
                        pairs[count].i = (uint32_t)i;
                        pairs[count].j = (uint32_t)j;
                        pairs[count].d2 = d2;
                    }
                    count += 1;
                }
            }
        }
    } else {
        for (size_t i = 0; i < num_points - 1; ++i) {
            float x = in_x[i];
            float y = in_y[i];
            float z = in_z[i];
            for (size_t j = i + 1; j < num_points; ++j) {
                float dx = in_x[j] - x;
                float dy = in_y[j] - y;
                float dz = in_z[j] - z;

                float d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= r2) {
                    if (pairs) {
                        pairs[count].i = (uint32_t)i;
                        pairs[count].j = (uint32_t)j;
                        pairs[count].d2 = d2;
                    }
                    count += 1;
                }
            }
        }
    }

    return count;
}

static inline double rnd_rng(double min, double max) {
    double r = ((double)rand() / (double)RAND_MAX);
    return r * (max - min) + min;
}

static int compare_dist_pair(void const* a, void const* b) {
    const dist_pair_t* pa = a;
    const dist_pair_t* pb = b;
    if (pa->i != pb->i) return pa->i - pb->i;
    if (pa->j != pb->j) return pa->j - pb->j;
    return 0;
}

static void n2_callback(const float* j_x, const float* j_y, const float* j_z, const uint32_t* j_idx, md_256 d2, int mask, uint32_t i, void* user_param) {
    size_t* count = (size_t*)user_param;
    *count += popcnt32(mask);
};

UTEST(spatial_hash, n2) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), alloc));

    md_system_t sys = {0};
    ASSERT_TRUE(md_gro_system_init(&sys, &gro_data, alloc));

    vec3_t* xyz = md_alloc(alloc, sizeof(vec3_t) * sys.atom.count);
    for (size_t i = 0; i < sys.atom.count; ++i) {
        xyz[i] = vec3_set(sys.atom.x[i], sys.atom.y[i], sys.atom.z[i]);
    }

    const vec4_t mask = md_unitcell_pbc_mask_vec4(&sys.unitcell);
    mat4x3_t I = mat4x3_from_mat3(md_unitcell_inv_basis_mat3(&sys.unitcell));
    mat4x3_t M = mat4x3_from_mat3(md_unitcell_basis_mat3(&sys.unitcell));

    {

#define test_count 2048
        float x[test_count];
        float y[test_count];
        float z[test_count];

        const double ext[3]  = {103.0, 103.0, 103.0};
	    const double min_rad = 3.0;
		const double max_rad = 6.0;

        srand(0);
        for (size_t i = 0; i < test_count; ++i) {
            x[i] = rnd_rng(0.0, ext[0]);
            y[i] = rnd_rng(0.0, ext[1]);
            z[i] = rnd_rng(0.0, ext[2]);
        }

        dist_pair_t bf_pairs[4096];
        dist_pair_t sa_pairs[4096];

        md_unitcell_t test_cell = md_unitcell_from_extent(ext[0], ext[1], ext[2]);

        md_spatial_acc_t sa = { .alloc = alloc };
		md_spatial_acc_init(&sa, x, y, z, test_count, max_rad, &test_cell);

        for (double rad = 3.0; rad <= 6.0; rad += 0.5) {
            size_t bf_count = do_brute_force(x, y, z, test_count, rad, &test_cell, bf_pairs);

            uint32_t sa_count = 0;
            md_spatial_acc_for_each_pair_within_cutoff(&sa, rad, spatial_acc_cutoff_callback, &sa_count);

            EXPECT_EQ(bf_count, sa_count);
            if (bf_count != sa_count) {
                printf("wierd expected: %zu, but got: %u\n", bf_count, sa_count);
            }
#if 0
            if (bf_count != sh_count) {
                qsort(bf_pairs, bf_count, sizeof(dist_pair_t), compare_dist_pair);
                qsort(sh_pairs, sh_count, sizeof(dist_pair_t), compare_dist_pair);
                
                printf("D2: %f\n", rad * rad);
                printf("Box ext: %f %f %f\n", ext, ext, ext);
                size_t i = 0, j = 0;
                for (; i < bf_count && j < sh_count;) {
                    int c = compare_dist_pair(&bf_pairs[i], &sh_pairs[j]);
                    if (c == 0) {
                        i++; j++;
                    } else if (c < 0) {
                        float x_i = x[bf_pairs[i].i];
                        float y_i = y[bf_pairs[i].i];
                        float z_i = z[bf_pairs[i].i];
                        float x_j = x[bf_pairs[i].j];
                        float y_j = y[bf_pairs[i].j];
                        float z_j = z[bf_pairs[i].j];

                        int cx_i = (int)floorf(x_i / rad);
                        int cy_i = (int)floorf(y_i / rad);
                        int cz_i = (int)floorf(z_i / rad);
                        int cx_j = (int)floorf(x_j / rad);
                        int cy_j = (int)floorf(y_j / rad);
                        int cz_j = (int)floorf(z_j / rad);

                        printf("BF only: %u %u %.5f\n", bf_pairs[i].i, bf_pairs[i].j, bf_pairs[i].d2);
                        printf("i coord: %f %f %f [%i,%i,%i]\n", x_i, y_i, z_i, cx_i, cy_i, cz_i);
                        printf("j coord: %f %f %f [%i,%i,%i]\n", x_j, y_j, z_j, cx_j, cy_j, cz_j);
                        i++;
                    } else {
                        printf("SH only: %u %u %.5f\n", sh_pairs[j].i, sh_pairs[j].j, sh_pairs[j].d2);
                        printf("i coord: %f %f %f\n", x[sh_pairs[j].i], y[sh_pairs[j].i], z[sh_pairs[j].i]);
                        printf("j coord: %f %f %f\n", x[sh_pairs[j].j], y[sh_pairs[j].j], z[sh_pairs[j].j]);
                        j++;
                    }
                }
                printf("WIERD\n");
            }
#endif
        }
#undef test_count
    }

#if 1
    md_unitcell_t cell = sys.unitcell;
    // Clear pbc flags

    size_t expected_count = 3711880;
    if (false) {
        cell.flags &= ~(MD_UNITCELL_PBC_X | MD_UNITCELL_PBC_Y | MD_UNITCELL_PBC_Z);
        expected_count = 3701958;
    }

    md_timestamp_t start, end;

    // Custom implementation of pairwise periodic N^2
#if 0
    dist_pair_t* pairs = md_alloc(alloc, sizeof(dist_pair_t) * 10000000);
    md_timestamp_t start = md_time_current();
    //size_t count = do_pairwise_periodic(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 5.0f, &unit_cell);
    size_t custom_count = do_pairwise_periodic_triclinic(sys.atom.x, sys.atom.y, sys.atom.z, sys.atom.count, 5.0f, &cell, pairs);
    md_timestamp_t end = md_time_current();
    printf("Custom: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, custom_count);
    if (custom_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, custom_count);
    }
    #endif
    
    uint32_t count = 0;

    // Current implementation of spatial hash for N^2
    start = md_time_current();
    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count, &cell, alloc);
    md_spatial_hash_query_multi_batch(spatial_hash, xyz, sys.atom.count, 5.0f, iter_batch_fn, &count);
    end = md_time_current();
    printf("Spatial hash: %f ms\n", md_time_as_milliseconds(end - start));
    size_t sh_count = count;
    ASSERT_TRUE(spatial_hash);
    EXPECT_EQ(expected_count, sh_count);
    if (sh_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, sh_count);
    }

    count = 0;

    // Spatial acc implementation
    start = md_time_current();
    md_spatial_acc_t acc = {.alloc = alloc};
    md_spatial_acc_init(&acc, sys.atom.x, sys.atom.y, sys.atom.z, sys.atom.count, 5.0, &cell);
    //md_spatial_acc_for_each_pair_within_cutoff(&acc, 5.0f, spatial_acc_pair_callback, &count);
    md_spatial_acc_for_each_pair_in_neighboring_cells(&acc, spatial_acc_neighbor_callback, &count);
    end = md_time_current();
    size_t sa_count = count;
    //end = md_time_current();
    printf("Spatial acc cell neighborhood: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, sa_count);
    if (sa_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, sa_count);
    }

    // Brute force
    start = md_time_current();
    size_t bf_count = do_brute_force(sys.atom.x, sys.atom.y, sys.atom.z, sys.atom.count, 5.0f, &cell, NULL);
    end = md_time_current();
    printf("Brute force: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, bf_count);
    if (bf_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, bf_count);
    }
#endif

    md_arena_allocator_destroy(alloc);
}

struct spatial_hash {
    md_allocator_i* arena;
    md_system_t sys;
    vec3_t pbc_ext;
};

UTEST_F_SETUP(spatial_hash) {
    utest_fixture->arena = md_vm_arena_create(GIGABYTES(4));

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), utest_fixture->arena));
    ASSERT_TRUE(md_gro_system_init(&utest_fixture->sys, &gro_data, utest_fixture->arena));
    utest_fixture->pbc_ext = md_unitcell_diag_vec3(&utest_fixture->sys.unitcell);
}

UTEST_F_TEARDOWN(spatial_hash) {
    md_vm_arena_destroy(utest_fixture->arena);
}

static inline float rnd() {
    return rand() / (float)RAND_MAX;
}

UTEST_F(spatial_hash, test_correctness_centered) {
    md_system_t* sys = &utest_fixture->sys;
    md_allocator_i* alloc = utest_fixture->arena;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(sys->atom.x, sys->atom.y, sys->atom.z, NULL, sys->atom.count, NULL, alloc);
    ASSERT_TRUE(spatial_hash);
    
    srand(31);
    
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(utest_fixture->arena);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (int64_t i = 0; i < sys->atom.count; ++i) {
            vec4_t c = {sys->atom.x[i], sys->atom.y[i], sys->atom.z[i], 0};

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

    md_system_t sys;
    ASSERT_TRUE(md_pdb_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count, &sys.unitcell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = md_unitcell_diag_vec4(&sys.unitcell);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20.f;
        const float rad2 = radius * radius;

        int ref_count = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const vec4_t c = {sys.atom.x[i], sys.atom.y[i], sys.atom.z[i], 0};
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

    md_system_t sys;
    ASSERT_TRUE(md_pdb_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, alloc));

    vec3_t* xyz = md_array_create(vec3_t, sys.atom.count, alloc);
    for (int i = 0; i < sys.atom.count; ++i) {
        xyz[i] = vec3_set(sys.atom.x[i], sys.atom.y[i], sys.atom.z[i]);
    }

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_vec3(xyz, NULL, md_array_size(xyz), &sys.unitcell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = md_unitcell_diag_vec4(&sys.unitcell);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const vec4_t c = {sys.atom.x[i], sys.atom.y[i], sys.atom.z[i], 0};
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

    md_system_t sys;
    ASSERT_TRUE(md_gro_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro"), NULL, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count, NULL, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = md_unitcell_diag_vec4(&sys.unitcell);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const vec4_t c = {sys.atom.x[i], sys.atom.y[i], sys.atom.z[i], 0};
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

    md_system_t* sys = &utest_fixture->sys;
    vec3_t pbc_ext = utest_fixture->pbc_ext;

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(sys->atom.x, sys->atom.y, sys->atom.z, NULL, sys->atom.count, &sys->unitcell, alloc);
    ASSERT_TRUE(spatial_hash);
    
    srand(31);

    const int num_iter = 100;
    const vec4_t period = vec4_from_vec3(pbc_ext, 0);
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), pbc_ext);
        float radius = rnd() * 20.0f;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < sys->atom.count; ++i) {
            vec4_t c = {sys->atom.x[i], sys->atom.y[i], sys->atom.z[i], 0};

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

    md_system_t sys;
    ASSERT_TRUE(md_gro_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro"), NULL, alloc));

    md_spatial_hash_t* spatial_hash = md_spatial_hash_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count, &sys.unitcell, alloc);
    ASSERT_TRUE(spatial_hash);

    const vec4_t pbc_ext = md_unitcell_diag_vec4(&sys.unitcell);

    srand(31);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        vec3_t pos = vec3_mul(vec3_set(rnd(), rnd(), rnd()), vec3_from_vec4(pbc_ext));
        float radius = rnd() * 20;

        int ref_count = 0;
        const float rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const vec4_t c = {sys.atom.x[i], sys.atom.y[i], sys.atom.z[i], 0};
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