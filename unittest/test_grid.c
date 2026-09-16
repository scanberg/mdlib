#include "utest.h"

#include <core/md_grid.h>
#include <core/md_vec_math.h>
#include <core/md_allocator.h>

#include <math.h>

/* md_grid carries four transforms that are two inverse pairs, and the inverses are written out by
 * hand rather than computed - world_to_model transposes the orientation and negates the origin
 * instead of inverting a matrix. That is correct only while the orientation stays a rotation, and
 * it is the kind of shortcut that silently stops holding. So every test here uses a grid that is
 * rotated, off-origin and anisotropically spaced: with an identity orientation a transpose and an
 * inverse are the same thing, and the tests would pass without checking anything. */
static md_grid_t test_grid(void) {
    /* 30 degrees about an axis that is not one of the basis vectors, so no row or column of the
     * orientation is trivial. */
    const vec3_t axis = vec3_normalize(vec3_set(1.0f, 2.0f, -0.5f));

    md_grid_t grid = {
        .orientation = mat3_angle_axis(0.5235987756f, axis),
        .origin      = {{ -3.25f, 7.5f, 0.125f }},
        .spacing     = {{ 0.5f, 0.25f, 1.5f }},
        .dim         = { 4, 3, 2 },
    };
    return grid;
}

static float max_abs_diff(mat4_t A, mat4_t B) {
    float m = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const float d = fabsf(A.elem[i][j] - B.elem[i][j]);
            if (d > m) m = d;
        }
    }
    return m;
}

UTEST(grid, num_points_is_the_product_of_the_dimensions) {
    md_grid_t grid = test_grid();
    EXPECT_EQ((size_t)(4 * 3 * 2), md_grid_num_points(&grid));

    /* A null grid answers 0 rather than dereferencing, which is what lets a caller size a buffer
     * before it has decided on a grid. */
    EXPECT_EQ((size_t)0, md_grid_num_points(NULL));

    /* An empty dimension collapses the whole thing, it does not wrap or saturate. */
    grid.dim[1] = 0;
    EXPECT_EQ((size_t)0, md_grid_num_points(&grid));
}

UTEST(grid, extent_and_center_follow_spacing_and_dimensions) {
    const md_grid_t grid = test_grid();

    const vec3_t extent = md_grid_extent(&grid);
    EXPECT_NEAR(4 * 0.5f,  extent.x, 1.0e-6f);
    EXPECT_NEAR(3 * 0.25f, extent.y, 1.0e-6f);
    EXPECT_NEAR(2 * 1.5f,  extent.z, 1.0e-6f);

    /* The centre is the origin plus half the extent, in world units and without the orientation
     * applied - worth pinning, because reading it as the centre of the rotated box is a natural
     * mistake and the two only agree for an axis-aligned grid. */
    const vec3_t center = md_grid_center(&grid);
    EXPECT_NEAR(-3.25f  + 1.0f,   center.x, 1.0e-6f);
    EXPECT_NEAR( 7.5f   + 0.375f, center.y, 1.0e-6f);
    EXPECT_NEAR( 0.125f + 1.5f,   center.z, 1.0e-6f);

    const vec3_t origin = md_grid_origin(&grid);
    EXPECT_NEAR(grid.origin.x, origin.x, 0.0f);
    EXPECT_NEAR(grid.origin.y, origin.y, 0.0f);
    EXPECT_NEAR(grid.origin.z, origin.z, 0.0f);
}

UTEST(grid, model_and_world_transforms_are_inverses) {
    const md_grid_t grid = test_grid();

    const mat4_t to_world = md_grid_model_to_world(&grid);
    const mat4_t to_model = md_grid_world_to_model(&grid);

    EXPECT_LT(max_abs_diff(mat4_mul(to_world, to_model), mat4_ident()), 1.0e-5f);
    EXPECT_LT(max_abs_diff(mat4_mul(to_model, to_world), mat4_ident()), 1.0e-5f);

    /* The hand-written inverse must also agree with a computed one. */
    EXPECT_LT(max_abs_diff(to_model, mat4_inverse(to_world)), 1.0e-5f);

    /* Model space carries no scaling, so a unit step in model space is a unit step in world space
     * however the grid is spaced. */
    const vec4_t a = mat4_mul_vec4(to_world, vec4_set(0, 0, 0, 1));
    const vec4_t b = mat4_mul_vec4(to_world, vec4_set(1, 0, 0, 1));
    EXPECT_NEAR(1.0f, vec3_length(vec3_sub(vec3_from_vec4(b), vec3_from_vec4(a))), 1.0e-5f);
}

UTEST(grid, index_and_world_transforms_are_inverses) {
    const md_grid_t grid = test_grid();

    const mat4_t to_world = md_grid_index_to_world(&grid);
    const mat4_t to_index = md_grid_world_to_index(&grid);

    EXPECT_LT(max_abs_diff(mat4_mul(to_world, to_index), mat4_ident()), 1.0e-5f);
    EXPECT_LT(max_abs_diff(mat4_mul(to_index, to_world), mat4_ident()), 1.0e-5f);
    EXPECT_LT(max_abs_diff(to_index, mat4_inverse(to_world)), 1.0e-5f);

    /* Index space does carry the spacing, and per axis: one voxel along x is spacing.x of world
     * distance, not the length of the spacing vector. Anisotropic spacing is what separates a
     * correct implementation from one that scales uniformly. */
    const vec4_t o  = mat4_mul_vec4(to_world, vec4_set(0, 0, 0, 1));
    const vec4_t dx = mat4_mul_vec4(to_world, vec4_set(1, 0, 0, 1));
    const vec4_t dy = mat4_mul_vec4(to_world, vec4_set(0, 1, 0, 1));
    const vec4_t dz = mat4_mul_vec4(to_world, vec4_set(0, 0, 1, 1));
    EXPECT_NEAR(0.50f, vec3_length(vec3_sub(vec3_from_vec4(dx), vec3_from_vec4(o))), 1.0e-5f);
    EXPECT_NEAR(0.25f, vec3_length(vec3_sub(vec3_from_vec4(dy), vec3_from_vec4(o))), 1.0e-5f);
    EXPECT_NEAR(1.50f, vec3_length(vec3_sub(vec3_from_vec4(dz), vec3_from_vec4(o))), 1.0e-5f);

    /* Index (0,0,0) is the origin. */
    EXPECT_NEAR(grid.origin.x, o.x, 1.0e-5f);
    EXPECT_NEAR(grid.origin.y, o.y, 1.0e-5f);
    EXPECT_NEAR(grid.origin.z, o.z, 1.0e-5f);
}

UTEST(grid, extracted_points_match_the_index_transform) {
    const md_grid_t grid = test_grid();
    const size_t    num  = md_grid_num_points(&grid);
    ASSERT_EQ((size_t)24, num);

    float xyz[24 * 3];
    md_grid_extract_points(xyz, &grid);

    /* The extraction is just index_to_world applied over the lattice, and x varies fastest. The
     * ordering matters to every consumer that indexes the volume alongside the points, so it is
     * asserted rather than assumed. */
    const mat4_t to_world = md_grid_index_to_world(&grid);
    for (int z = 0; z < grid.dim[2]; ++z) {
        for (int y = 0; y < grid.dim[1]; ++y) {
            for (int x = 0; x < grid.dim[0]; ++x) {
                const int i = (z * grid.dim[1] + y) * grid.dim[0] + x;
                const vec4_t w = mat4_mul_vec4(to_world, vec4_set((float)x, (float)y, (float)z, 1.0f));
                EXPECT_NEAR(w.x, xyz[i * 3 + 0], 1.0e-5f);
                EXPECT_NEAR(w.y, xyz[i * 3 + 1], 1.0e-5f);
                EXPECT_NEAR(w.z, xyz[i * 3 + 2], 1.0e-5f);
            }
        }
    }

    /* And the round trip closes: every extracted point maps back to the integer index it came from. */
    const mat4_t to_index = md_grid_world_to_index(&grid);
    for (size_t i = 0; i < num; ++i) {
        const vec4_t idx = mat4_mul_vec4(to_index, vec4_set(xyz[i*3+0], xyz[i*3+1], xyz[i*3+2], 1.0f));
        EXPECT_NEAR(roundf(idx.x), idx.x, 1.0e-3f);
        EXPECT_NEAR(roundf(idx.y), idx.y, 1.0e-3f);
        EXPECT_NEAR(roundf(idx.z), idx.z, 1.0e-3f);
    }
}
