#include <core/md_grid.h>
#include <core/md_common.h>
#include <core/md_vec_math.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_grid_extract_points(float* out_xyz, const md_grid_t* grid) {
    MD_ASSERT(out_xyz);
    MD_ASSERT(grid);

    md_f32x4x4_t index_to_world = md_grid_index_to_world(grid);

    for (int z = 0; z < grid->dim[2]; ++z) {
        for (int y = 0; y < grid->dim[1]; ++y) {
            for (int x = 0; x < grid->dim[0]; ++x) {
                int i = (z * grid->dim[1] + y) * grid->dim[0] + x;
                md_f32x4_t v = md_f32x4_set((float)x, (float)y, (float)z, 1.0f);
                md_f32x4_t w = md_f32x4x4_mul_f32x4(index_to_world, v);
                out_xyz[i*3 + 0] = w.x;
                out_xyz[i*3 + 1] = w.y;
                out_xyz[i*3 + 2] = w.z;
            }
        }
    }
}

md_f32x4x4_t md_grid_model_to_world(const md_grid_t* grid) {
    MD_ASSERT(grid);
    md_f32x4x4_t R = md_f32x4x4_from_f32x3x3(grid->orientation);
    md_f32x4x4_t T = md_f32x4x4_translate_f32x3(grid->origin);
    return md_f32x4x4_mul(T, R);
}

md_f32x4x4_t md_grid_world_to_model(const md_grid_t* grid) {
    MD_ASSERT(grid);
    // These are the inverse matrices
    md_f32x4x4_t R = md_f32x4x4_from_f32x3x3(md_f32x3x3_transpose(grid->orientation));
    md_f32x4x4_t T = md_f32x4x4_translate_f32x3(md_f32x3_mul1(grid->origin, -1.0f));
    return md_f32x4x4_mul(R, T);
}

md_f32x4x4_t md_grid_index_to_world(const md_grid_t* grid) {
    MD_ASSERT(grid);

    md_f32x4x4_t S = md_f32x4x4_scale_f32x3(grid->spacing);
    md_f32x4x4_t R = md_f32x4x4_from_f32x3x3(grid->orientation);
    md_f32x4x4_t T = md_f32x4x4_translate_f32x3(grid->origin);

    return md_f32x4x4_mul(T, md_f32x4x4_mul(R, S));
}

md_f32x4x4_t md_grid_world_to_index(const md_grid_t* grid) {
    MD_ASSERT(grid);

    // These are the inverse matrices
    md_f32x4x4_t S = md_f32x4x4_scale_f32x3(md_f32x3_div(md_f32x3_set1(1.0f), grid->spacing));
    md_f32x4x4_t R = md_f32x4x4_from_f32x3x3(md_f32x3x3_transpose(grid->orientation));
    md_f32x4x4_t T = md_f32x4x4_translate_f32x3(md_f32x3_mul1(grid->origin, -1.0f));

    return md_f32x4x4_mul(S, md_f32x4x4_mul(R, T));
}

#ifdef __cplusplus
}
#endif
