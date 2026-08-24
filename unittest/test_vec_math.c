#include "utest.h"

#include <core/md_vec_math.h>

UTEST(vec_math, mat3) {
	mat3_t M = {
		1,2,3,
		4,3,2,
		3,4,1,
	};

	vec3_t v[2] = {
		{1,1,1},
		{2,2,2},
	};

	float x[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2};
	float y[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2};
	float z[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2};

	float x_[17];
	float y_[17];
	float z_[17];

	mat3_batch_transform(x_, y_, z_, x, y, z, 17, M);
	mat3_batch_transform_inplace(x, y, z, 17, M);

	v[0] = mat3_mul_vec3(M, v[0]);
	v[1] = mat3_mul_vec3(M, v[1]);

	EXPECT_EQ(8, x[0]);
	EXPECT_EQ(9, y[0]);
	EXPECT_EQ(6, z[0]);

	EXPECT_EQ(8, x_[0]);
	EXPECT_EQ(9, y_[0]);
	EXPECT_EQ(6, z_[0]);

	EXPECT_EQ(8, v[0].x);
	EXPECT_EQ(9, v[0].y);
	EXPECT_EQ(6, v[0].z);

	EXPECT_EQ(16, x[16]);
	EXPECT_EQ(18, y[16]);
	EXPECT_EQ(12, z[16]);

	EXPECT_EQ(16, x_[16]);
	EXPECT_EQ(18, y_[16]);
	EXPECT_EQ(12, z_[16]);

	EXPECT_EQ(16, v[1].x);
	EXPECT_EQ(18, v[1].y);
	EXPECT_EQ(12, v[1].z);
}

UTEST(vec_math, mat4) {
	mat4_t M = {
		1,2,3,4,
		4,3,2,1,
		3,4,1,2,
		4,3,2,1,
	};

	vec4_t v[2] = {
		{1,1,1,1},
		{2,2,2,1},
	};

	float x[2] = {1,2};
	float y[2] = {1,2};
	float z[2] = {1,2};

	float x_[2];
	float y_[2];
	float z_[2];

	mat4_batch_transform(x_, y_, z_, x, y, z, 1.0f, 2, M);
	mat4_batch_transform_inplace(x, y, z, 1.0f, 2, M);

	v[0] = mat4_mul_vec4(M, v[0]);
	v[1] = mat4_mul_vec4(M, v[1]);

	EXPECT_EQ(12, x[0]);
	EXPECT_EQ(12, y[0]);
	EXPECT_EQ(8,  z[0]);

	EXPECT_EQ(12, x_[0]);
	EXPECT_EQ(12, y_[0]);
	EXPECT_EQ(8,  z_[0]);

	EXPECT_EQ(12, v[0].x);
	EXPECT_EQ(12, v[0].y);
	EXPECT_EQ(8,  v[0].z);

	EXPECT_EQ(20, x[1]);
	EXPECT_EQ(21, y[1]);
	EXPECT_EQ(14, z[1]);

	EXPECT_EQ(20, x_[1]);
	EXPECT_EQ(21, y_[1]);
	EXPECT_EQ(14, z_[1]);

	EXPECT_EQ(20, v[1].x);
	EXPECT_EQ(21, v[1].y);
	EXPECT_EQ(14, v[1].z);
}

// mat4x3_t is THREE columns of four rows. The product combines those three columns using v.x, v.y
// and v.z only: v.w takes no part, and there is no fourth column to reach for. The scalar fallback
// used to read M.elem[3][..], which is past the end of the struct, so this pins whichever path the
// build selects to the same reference computation.
UTEST(vec_math, mat4x3_mul_vec4) {
	mat4x3_t M;
	M.col[0] = vec4_set( 1.0f,  2.0f,  3.0f,  4.0f);
	M.col[1] = vec4_set( 5.0f,  6.0f,  7.0f,  8.0f);
	M.col[2] = vec4_set( 9.0f, 10.0f, 11.0f, 12.0f);   // a non zero w row, so a path that drops it is caught

	const vec4_t v = vec4_set(0.5f, -1.5f, 2.0f, 1000.0f);   // v.w is large on purpose: it must not leak in
	const vec4_t r = mat4x3_mul_vec4(M, v);

	for (int row = 0; row < 4; ++row) {
		const float expect = M.elem[0][row] * v.x + M.elem[1][row] * v.y + M.elem[2][row] * v.z;
		EXPECT_NEAR(r.elem[row], expect, 1.0e-5f);
	}

	// changing v.w must change nothing at all
	const vec4_t r2 = mat4x3_mul_vec4(M, vec4_set(v.x, v.y, v.z, -7.0f));
	for (int row = 0; row < 4; ++row) {
		EXPECT_NEAR(r2.elem[row], r.elem[row], 1.0e-6f);
	}

	// built from a mat3, it has to agree with mat3_mul_vec3 on xyz and leave w at zero
	mat3_t A = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	const mat4x3_t B = mat4x3_from_mat3(A);
	const vec3_t   e = mat3_mul_vec3(A, vec3_set(v.x, v.y, v.z));
	const vec4_t   g = mat4x3_mul_vec4(B, v);
	EXPECT_NEAR(g.x, e.x, 1.0e-5f);
	EXPECT_NEAR(g.y, e.y, 1.0e-5f);
	EXPECT_NEAR(g.z, e.z, 1.0e-5f);
	EXPECT_NEAR(g.w, 0.0f, 1.0e-6f);
}
