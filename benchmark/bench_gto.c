#include "ubench.h"

#include <md_vlx.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>

#include <md_gto.h>
//#include <md_gto.c>

#include <float.h>

#define BLK_DIM 8

UBENCH_EX(gto, evaluate_grid) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    
    md_vlx_data_t vlx = {0};
    md_vlx_data_parse_file(&vlx, STR_LIT(MD_BENCHMARK_DATA_DIR "/vlx/mol.out"), arena);

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);

    for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
        vec3_t c = {(float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]};
        min_box = vec3_min(min_box, c);
        max_box = vec3_max(max_box, c);
    }

    min_box = vec3_sub_f(min_box, 2.0f);
    max_box = vec3_add_f(max_box, 2.0f);

    // Conversion from Ångström to Bohr
    const float factor = 1.0 / 0.529177210903;

    min_box = vec3_mul_f(min_box, factor);
    max_box = vec3_mul_f(max_box, factor);

    int vol_dim = 256;
    size_t bytes = sizeof(float) * vol_dim * vol_dim * vol_dim;
    float* vol_data = md_arena_allocator_push(arena, bytes);
    MEMSET(vol_data, 0, bytes);

    vec3_t step = vec3_div_f(vec3_sub(max_box, min_box), (float)vol_dim);

    md_grid_t grid = {
        .data = vol_data,
        .dim = {vol_dim, vol_dim, vol_dim},
        .origin = {min_box.x, min_box.y, min_box.z},
        .step_x = {step.x, 0, 0},
        .step_y = {0, step.y, 0},
        .step_z = {0, 0, step.z},
    };

    size_t num_gtos = md_vlx_mo_gto_count(&vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mo_gto_extract(gtos, &vlx, 120);
    md_gto_cutoff_compute(gtos, num_gtos, 1.0e-6);

    md_gto_t* sub_gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);

    int num_blk = vol_dim / BLK_DIM;
    int num_tot_blk = num_blk * num_blk * num_blk;

    UBENCH_DO_BENCHMARK() {
        for (int i = 0; i < num_tot_blk; ++i) {
            // Determine block (x/y/z) index from linear index (i)
            int blk_z = i & (num_blk - 1);
            int blk_y = (i / num_blk) & (num_blk - 1);
            int blk_x = i / (num_blk * num_blk);

            const int off_idx[3] = {blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM};
            const int len_idx[3] = {BLK_DIM, BLK_DIM, BLK_DIM};

            float aabb_min[3] = {
                grid.origin[0] + off_idx[0] * grid.step_x[0],
                grid.origin[1] + off_idx[1] * grid.step_y[1],
                grid.origin[2] + off_idx[2] * grid.step_z[2],
            };
            float aabb_max[3] = {
                grid.origin[0] + (off_idx[0] + len_idx[0]) * grid.step_x[0],
                grid.origin[1] + (off_idx[1] + len_idx[1]) * grid.step_y[1],
                grid.origin[2] + (off_idx[2] + len_idx[2]) * grid.step_z[2],
            };

            size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);
            
            md_gto_grid_evaluate_sub(&grid, off_idx, len_idx, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
            //evaluate_grid_ortho_8x8x8_256(grid.data, off_idx, len_idx, grid.origin, step.elem, sub_pgtos, num_sub_pgtos, MD_GTO_EVAL_MODE_PSI);
            //evaluate_grid_8x8x8_256(grid.data, off_idx, len_idx, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        }
    }

    UBENCH_DO_NOTHING(&grid);

    md_arena_allocator_destroy(arena);
}
