#include <core/md_grid.h>
#include <core/md_common.h>
#include <core/md_vec_math.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_grid_extract_points(float* out_xyz, const md_grid_t* grid) {
    MD_ASSERT(out_xyz);
    MD_ASSERT(grid);

    md_float4x4_t index_to_world = md_grid_index_to_world(grid);

    for (int z = 0; z < grid->dim[2]; ++z) {
        for (int y = 0; y < grid->dim[1]; ++y) {
            for (int x = 0; x < grid->dim[0]; ++x) {
                int i = (z * grid->dim[1] + y) * grid->dim[0] + x;
                md_float4_t v = md_float4_set((float)x, (float)y, (float)z, 1.0f);
                md_float4_t w = md_float4x4_mul_float4(index_to_world, v);
                out_xyz[i*3 + 0] = w.x;
                out_xyz[i*3 + 1] = w.y;
                out_xyz[i*3 + 2] = w.z;
            }
        }
    }
}

md_float4x4_t md_grid_model_to_world(const md_grid_t* grid) {
    MD_ASSERT(grid);
    md_float4x4_t R = md_float4x4_from_float3x3(grid->orientation);
    md_float4x4_t T = md_float4x4_translate_float3(grid->origin);
    return md_float4x4_mul(T, R);
}

md_float4x4_t md_grid_world_to_model(const md_grid_t* grid) {
    MD_ASSERT(grid);
    // These are the inverse matrices
    md_float4x4_t R = md_float4x4_from_float3x3(md_float3x3_transpose(grid->orientation));
    md_float4x4_t T = md_float4x4_translate_float3(md_float3_mul1(grid->origin, -1.0f));
    return md_float4x4_mul(R, T);
}

md_float4x4_t md_grid_index_to_world(const md_grid_t* grid) {
    MD_ASSERT(grid);

    md_float4x4_t S = md_float4x4_scale_float3(grid->spacing);
    md_float4x4_t R = md_float4x4_from_float3x3(grid->orientation);
    md_float4x4_t T = md_float4x4_translate_float3(grid->origin);

    return md_float4x4_mul(T, md_float4x4_mul(R, S));
}

md_float4x4_t md_grid_world_to_index(const md_grid_t* grid) {
    MD_ASSERT(grid);

    // These are the inverse matrices
    md_float4x4_t S = md_float4x4_scale_float3(md_float3_div(md_float3_set1(1.0f), grid->spacing));
    md_float4x4_t R = md_float4x4_from_float3x3(md_float3x3_transpose(grid->orientation));
    md_float4x4_t T = md_float4x4_translate_float3(md_float3_mul1(grid->origin, -1.0f));

    return md_float4x4_mul(S, md_float4x4_mul(R, T));
}

#ifdef __cplusplus
}
#endif
