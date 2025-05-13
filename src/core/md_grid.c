#include <core/md_grid.h>
#include <core/md_common.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_grid_extract_points(float* out_xyz, const md_grid_t* grid) {
    ASSERT(out_xyz);
    ASSERT(grid);

    mat4_t index_to_world = md_grid_index_to_world(grid);

    for (int z = 0; z < grid->dim[2]; ++z) {
        for (int y = 0; y < grid->dim[1]; ++y) {
            for (int x = 0; x < grid->dim[0]; ++x) {
                int i = (z * grid->dim[1] + y) * grid->dim[0] + x;
                vec4_t v = vec4_set((float)x, (float)y, (float)z, 1.0f);
                vec4_t w = mat4_mul_vec4(index_to_world, v);
                out_xyz[i*3 + 0] = w.x;
                out_xyz[i*3 + 1] = w.y;
                out_xyz[i*3 + 2] = w.z;
            }
        }
    }
}

mat4_t md_grid_model_to_world(const md_grid_t* grid) {
    ASSERT(grid);
    mat4_t R = mat4_from_mat3(grid->orientation);
    mat4_t T = mat4_translate_vec3(grid->origin);
    return mat4_mul(T, R);
}

mat4_t md_grid_world_to_model(const md_grid_t* grid) {
    ASSERT(grid);
    // These are the inverse matrices
    mat4_t R = mat4_from_mat3(mat3_transpose(grid->orientation));
    mat4_t T = mat4_translate_vec3(vec3_mul_f(grid->origin, -1.0f));
    return mat4_mul(R, T);
}

mat4_t md_grid_index_to_world(const md_grid_t* grid) {
    ASSERT(grid);

    mat4_t S = mat4_scale_vec3(grid->spacing);
    mat4_t R = mat4_from_mat3(grid->orientation);
    mat4_t T = mat4_translate_vec3(grid->origin);

    return mat4_mul(T, mat4_mul(R, S));
}

mat4_t md_grid_world_to_index(const md_grid_t* grid) {
    ASSERT(grid);

    // These are the inverse matrices
    mat4_t S = mat4_scale_vec3(vec3_div(vec3_set1(1.0f), grid->spacing));
    mat4_t R = mat4_from_mat3(mat3_transpose(grid->orientation));
    mat4_t T = mat4_translate_vec3(vec3_mul_f(grid->origin, -1.0f));

    return mat4_mul(S, mat4_mul(R, T));
}

#ifdef __cplusplus
extern "C" {
#endif