#include "utest.h"

#include <core/md_arena_allocator.h>
#include <md_vlx.h>
#include <md_gto.c>

// Conversion from Ångström to Bohr
const float factor = 1.0 / 0.529177210903;
#define BLK_DIM 8

static void init(md_grid_t* grid, md_gto_t** gtos, size_t* num_gtos, int vol_dim, md_allocator_i* arena) {
    md_vlx_data_t vlx = {0};
    md_vlx_data_parse_file(&vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out"), arena);

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);

    for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
        vec3_t c = {(float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]};
        min_box = vec3_min(min_box, c);
        max_box = vec3_max(max_box, c);
    }

    min_box = vec3_sub_f(min_box, 2.0f);
    max_box = vec3_add_f(max_box, 2.0f);

    min_box = vec3_mul_f(min_box, factor);
    max_box = vec3_mul_f(max_box, factor);

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

    *num_gtos = md_vlx_mol_pgto_count(&vlx);
    *gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * *num_gtos);
    md_vlx_mol_pgto_extract(*gtos, &vlx, 120);
    md_gto_cutoff_compute(*gtos, *num_gtos, 1.0e-6);
}

UTEST(gto, evaluate_grid) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    
    md_grid_t grid;
    md_gto_t* gtos;
    size_t num_gtos;

    int vol_dim = 256;
    init(&grid, &gtos, &num_gtos, vol_dim, arena);

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

    for (int i = 0; i < num_tot_blk; ++i) {
        // Determine block index from linear input index i
        int blk_x = i % num_blk[0];
        int blk_y = (i / num_blk[0]) % num_blk[1];
        int blk_z = i / (num_blk[0] * num_blk[1]);

        int off_idx[3] = { blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM };
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

        evaluate_grid_ref(ref_grid.data, beg_idx, end_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);

        const float eps = 1.0e-5f;

#ifdef __SSE2__
        evaluate_grid_ortho_8x8x8_128(grid.data, beg_idx, grid.dim, grid.origin, step_ortho, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }

        evaluate_grid_8x8x8_128(grid.data, beg_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }
#endif

#ifdef __AVX2__
        evaluate_grid_ortho_8x8x8_256(grid.data, beg_idx, grid.dim, grid.origin, step_ortho, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }

        evaluate_grid_8x8x8_256(grid.data, beg_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }
#endif
    }

    md_arena_allocator_destroy(arena);
}