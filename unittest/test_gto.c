#include "utest.h"

#include <core/md_arena_allocator.h>
#include <md_vlx.h>
#include <md_cube.h>

#include <md_gto.h>

#include <float.h>

// Conversion from Ångström to Bohr
#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903

static void init(md_grid_t* grid, md_gto_t** gtos, size_t* num_gtos, int vol_dim, str_t filename, md_allocator_i* arena) {
    md_vlx_t* vlx = md_vlx_create(arena);
    md_vlx_parse_file(vlx, filename);

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);

	const dvec3_t* coords = md_vlx_atom_coordinates(vlx);
    for (size_t i = 0; i < md_vlx_number_of_atoms(vlx); ++i) {
        vec3_t c = {(float)coords[i].x, (float)coords[i].y, (float)coords[i].z};
        min_box = vec3_min(min_box, c);
        max_box = vec3_max(max_box, c);
    }

    min_box = vec3_sub_f(min_box, 2.0f);
    max_box = vec3_add_f(max_box, 2.0f);

    min_box = vec3_mul_f(min_box, ANGSTROM_TO_BOHR);
    max_box = vec3_mul_f(max_box, ANGSTROM_TO_BOHR);

    size_t bytes = sizeof(float) * vol_dim * vol_dim * vol_dim;
    float* vol_data = md_arena_allocator_push(arena, bytes);
    MEMSET(vol_data, 0, bytes);

    vec3_t step = vec3_div_f(vec3_sub(max_box, min_box), (float)vol_dim);

    *grid = (md_grid_t) {
        .data = vol_data,
        .dim = {vol_dim, vol_dim, vol_dim},
        .origin = {min_box.x, min_box.y, min_box.z},
        .step_x = {step.x, 0, 0},
        .step_y = {0, step.y, 0},
        .step_z = {0, 0, step.z},
    };

    *num_gtos = md_vlx_mo_gto_count(vlx);
    *gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * *num_gtos);
    md_vlx_mo_gto_extract(*gtos, vlx, 120, MD_VLX_MO_TYPE_ALPHA);
    md_gto_cutoff_compute(*gtos, *num_gtos, 1.0e-6);
}

UTEST(gto, evaluate_grid) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    
    md_grid_t grid;
    md_gto_t* gtos;
    size_t num_gtos;

    int vol_dim = 64;
    init(&grid, &gtos, &num_gtos, vol_dim, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out"), arena);

    md_grid_t ref_grid = grid;
    ref_grid.data = md_alloc(arena, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    MEMSET(grid.data, 0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);
    MEMSET(ref_grid.data, 0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    md_gto_t* sub_gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);

    // Number of NxNxN blocks in each dimension
    int num_blk[3] = {
        grid.dim[0] / BLK_DIM,
        grid.dim[1] / BLK_DIM,
        grid.dim[2] / BLK_DIM,
    };

    int num_tot_blk = num_blk[0] * num_blk[1] * num_blk[2];

    float step_ortho[3] = {
        grid.step_x[0],
        grid.step_y[1],
        grid.step_z[2],
    };

    int bi[3] = {0};

    for (; bi[2] < num_blk[2]; ++bi[2]) {
        for (; bi[1] < num_blk[1]; ++bi[1]) {
            for (; bi[0] < num_blk[0]; ++bi[0]) {
                int off_idx[3] = { bi[0] * BLK_DIM, bi[1] * BLK_DIM, bi[2] * BLK_DIM };
                int len_idx[3] = { BLK_DIM, BLK_DIM, BLK_DIM };

                int beg_idx[3] = {off_idx[0], off_idx[1], off_idx[2]};
                int end_idx[3] = {off_idx[0] + len_idx[0], off_idx[1] + len_idx[1], off_idx[2] + len_idx[2]};

                float aabb_min[3] = {
                    grid.origin[0] + beg_idx[0] * grid.step_x[0],
                    grid.origin[1] + beg_idx[1] * grid.step_y[1],
                    grid.origin[2] + beg_idx[2] * grid.step_z[2],
                };
                float aabb_max[3] = {
                    grid.origin[0] + end_idx[0] * grid.step_x[0],
                    grid.origin[1] + end_idx[1] * grid.step_y[1],
                    grid.origin[2] + end_idx[2] * grid.step_z[2],
                };

                size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);

                md_gto_grid_evaluate_sub(&grid, beg_idx, end_idx, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
            }
        }
    }

    md_arena_allocator_destroy(arena);
}

// compares vlx file of specific orbital and compares it with a cube file, returns the abs max difference between the two
static double compare_vlx_and_cube(const md_vlx_t* vlx, size_t mo_idx, double cutoff_value, const md_cube_t* cube, md_allocator_i* arena) {

    md_grid_t grid = {
        .data = md_arena_allocator_push(arena, sizeof(float) * cube->data.num_x * cube->data.num_y * cube->data.num_z),
        .dim = {cube->data.num_x, cube->data.num_y, cube->data.num_z},
        .origin = {cube->origin[0], cube->origin[1], cube->origin[2]},
        .step_x = {cube->xaxis[0], cube->xaxis[1], cube->xaxis[2]},
        .step_y = {cube->yaxis[0], cube->yaxis[1], cube->yaxis[2]},
        .step_z = {cube->zaxis[0], cube->zaxis[1], cube->zaxis[2]},
    };

    size_t num_gtos = md_vlx_mo_gto_count(vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mo_gto_extract(gtos, vlx, mo_idx, MD_VLX_MO_TYPE_ALPHA);
    md_gto_cutoff_compute(gtos, num_gtos, 0);

    size_t count = grid.dim[0] * grid.dim[1] * grid.dim[2];

    float* psi  = md_arena_allocator_push(arena, sizeof(float)  * count);
    vec3_t* xyz = md_arena_allocator_push(arena, sizeof(vec3_t) * count);
    for (int z = 0; z < grid.dim[2]; ++z) {
        for (int y = 0; y < grid.dim[1]; ++y) {
            for (int x = 0; x < grid.dim[0]; ++x) {
                int i = z * grid.dim[1] * grid.dim[0] + y * grid.dim[0] + x;
                xyz[i] = vec3_set(
                    grid.origin[0] + grid.step_x[0] * x,
                    grid.origin[1] + grid.step_y[1] * y,
                    grid.origin[2] + grid.step_z[2] * z);
            }
        }
    }

    md_gto_xyz_evaluate(psi, (float*)xyz, count, sizeof(vec3_t), gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    MEMSET(grid.data, 0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    md_gto_grid_evaluate(&grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    double xyz_sum  = 0.0;
    double grid_sum = 0.0;
    double cube_sum = 0.0;

    double max_delta  = 0.0;
    int    max_idx[3] = {0};

    for (int iz = 0; iz < grid.dim[2]; ++iz) {
        for (int iy = 0; iy < grid.dim[1]; ++iy) {
            for (int ix = 0; ix < grid.dim[0]; ++ix) {
                int grid_idx = iz * grid.dim[0] * grid.dim[1] + iy * grid.dim[0] + ix;
                int cube_idx = ix * grid.dim[1] * grid.dim[2] + iy * grid.dim[2] + iz;
                double g_val = grid.data[grid_idx];
                double c_val = cube->data.val[cube_idx];
                double p_val = psi[grid_idx];
                grid_sum += g_val;
                cube_sum += c_val;
                xyz_sum  += p_val;

                double delta = fabs(g_val - c_val);
                if (delta > max_delta) {
                    max_delta = delta;
                    max_idx[0] = ix;
                    max_idx[1] = iy;
                    max_idx[2] = iz;
                }
            }
        }
    }

    printf("Max delta: %g at [%i,%i,%i]\n", max_delta, max_idx[0], max_idx[1], max_idx[2]);
    printf("GRID SUM: %.5f\n", grid_sum);
    printf("CUBE SUM: %.5f\n", cube_sum);
    printf("XYZ  SUM: %.5f\n", xyz_sum);

    return max_delta;

    md_arena_allocator_destroy(arena);
}

UTEST(gto, h2o) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube_lumo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o_lumo.cube"), arena));

    md_vlx_t* vlx = md_vlx_create(arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, 0.0, &cube_lumo, arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);

    md_arena_allocator_destroy(arena);
}

#if 0
UTEST(gto, amide) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube_lumo = {0};
    md_cube_t cube_homo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide_lumo.cube"), arena));
    ASSERT_TRUE(md_cube_file_load(&cube_homo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide_homo.cube"), arena));

    md_vlx_t* vlx = md_vlx_create(arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx);
    size_t homo_idx = md_vlx_scf_homo_idx(vlx);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, 0.0, &cube_lumo, arena);
    double max_delta_homo = compare_vlx_and_cube(vlx, homo_idx, 0.0, &cube_homo, arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);
    EXPECT_LT(max_delta_homo, 1.0E-4);

    md_arena_allocator_destroy(arena);
}

UTEST(gto, ne) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube_lumo = {0};
    md_cube_t cube_homo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/ne_lumo.cube"), arena));
    ASSERT_TRUE(md_cube_file_load(&cube_homo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/ne_homo.cube"), arena));

    md_vlx_t* vlx = md_vlx_create(arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/ne.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx);
    size_t homo_idx = md_vlx_scf_homo_idx(vlx);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, 0.0, &cube_lumo, arena);
    double max_delta_homo = compare_vlx_and_cube(vlx, homo_idx, 0.0, &cube_homo, arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);
    EXPECT_LT(max_delta_homo, 1.0E-4);

    md_arena_allocator_destroy(arena);
}

UTEST(gto, myjob) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube_lumo = {0};
    md_cube_t cube_homo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/myjob_lumo.cube"), arena));
    ASSERT_TRUE(md_cube_file_load(&cube_homo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/myjob_homo.cube"), arena));

    md_vlx_t* vlx = md_vlx_create(arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/myjob.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx);
    size_t homo_idx = md_vlx_scf_homo_idx(vlx);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, 0.0, &cube_lumo, arena);
    double max_delta_homo = compare_vlx_and_cube(vlx, homo_idx, 0.0, &cube_homo, arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);
    EXPECT_LT(max_delta_homo, 1.0E-4);

    md_arena_allocator_destroy(arena);
}
#endif
