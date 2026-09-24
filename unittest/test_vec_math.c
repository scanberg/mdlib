#include "utest.h"

#include <core/md_vec_math.h>

// The batch forms against the single vector products, on 17 atoms so both the eight wide path and
// the scalar tail are covered, out of place and in place.
UTEST(vec_math, mat3) {
	mat3_t M = {
		1,2,3,
		4,3,2,
		3,4,1,
	};

	vec3_t xyz[17], out[17];
	for (int i = 0; i < 17; ++i) xyz[i] = vec3_set(1.0f + i, 2.0f - 0.5f * i, 0.25f * i);

	mat3_batch_transform(out, xyz, 17, M);
	for (int i = 0; i < 17; ++i) {
		const vec3_t r = mat3_mul_vec3(M, xyz[i]);
		EXPECT_NEAR(r.x, out[i].x, 1.0e-4f);
		EXPECT_NEAR(r.y, out[i].y, 1.0e-4f);
		EXPECT_NEAR(r.z, out[i].z, 1.0e-4f);
	}
	mat3_batch_transform_inplace(xyz, 17, M);
	for (int i = 0; i < 17; ++i) {
		EXPECT_EQ(out[i].x, xyz[i].x);
		EXPECT_EQ(out[i].y, xyz[i].y);
		EXPECT_EQ(out[i].z, xyz[i].z);
	}

	vec3_t one = {1,1,1};
	one = mat3_mul_vec3(M, one);
	EXPECT_EQ(8, one.x);
	EXPECT_EQ(9, one.y);
	EXPECT_EQ(6, one.z);
}

UTEST(vec_math, mat4) {
	mat4_t M = {
		1,2,3,4,
		4,3,2,1,
		3,4,1,2,
		4,3,2,1,
	};

	vec3_t xyz[17], out[17];
	for (int i = 0; i < 17; ++i) xyz[i] = vec3_set(1.0f + i, 2.0f - 0.5f * i, 0.25f * i);

	mat4_batch_transform(out, xyz, 1.0f, 17, M);
	for (int i = 0; i < 17; ++i) {
		const vec4_t r = mat4_mul_vec4(M, vec4_from_vec3(xyz[i], 1.0f));
		EXPECT_NEAR(r.x, out[i].x, 1.0e-4f);
		EXPECT_NEAR(r.y, out[i].y, 1.0e-4f);
		EXPECT_NEAR(r.z, out[i].z, 1.0e-4f);
	}
	mat4_batch_transform_inplace(xyz, 1.0f, 17, M);
	for (int i = 0; i < 17; ++i) {
		EXPECT_EQ(out[i].x, xyz[i].x);
		EXPECT_EQ(out[i].y, xyz[i].y);
		EXPECT_EQ(out[i].z, xyz[i].z);
	}
}

UTEST(vec_math, translate) {
	vec3_t xyz[19], out[19];
	for (int i = 0; i < 19; ++i) xyz[i] = vec3_set((float)i, 100.0f + i, -3.0f * i);
	const vec3_t t = {0.5f, -2.0f, 7.0f};
	vec3_batch_translate(out, xyz, 19, t);
	for (int i = 0; i < 19; ++i) {
		EXPECT_EQ(xyz[i].x + t.x, out[i].x);
		EXPECT_EQ(xyz[i].y + t.y, out[i].y);
		EXPECT_EQ(xyz[i].z + t.z, out[i].z);
	}
	vec3_batch_translate_inplace(xyz, 19, t);
	for (int i = 0; i < 19; ++i) {
		EXPECT_EQ(out[i].x, xyz[i].x);
		EXPECT_EQ(out[i].y, xyz[i].y);
		EXPECT_EQ(out[i].z, xyz[i].z);
	}
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
