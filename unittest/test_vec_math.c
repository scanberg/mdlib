#include "utest.h"

#include <core/md_vec_math.h>

// The batch forms against the single vector products, on 17 atoms so both the eight wide path and
// the scalar tail are covered, out of place and in place.
UTEST(vec_math, mat3) {
	md_f32x3x3_t M = {
		1,2,3,
		4,3,2,
		3,4,1,
	};

	md_f32x3_t xyz[17], out[17];
	for (int i = 0; i < 17; ++i) xyz[i] = md_f32x3_set(1.0f + i, 2.0f - 0.5f * i, 0.25f * i);

	md_f32x3x3_batch_transform(out, xyz, 17, M);
	for (int i = 0; i < 17; ++i) {
		const md_f32x3_t r = md_f32x3x3_mul_f32x3(M, xyz[i]);
		EXPECT_NEAR(r.x, out[i].x, 1.0e-4f);
		EXPECT_NEAR(r.y, out[i].y, 1.0e-4f);
		EXPECT_NEAR(r.z, out[i].z, 1.0e-4f);
	}
	md_f32x3x3_batch_transform_inplace(xyz, 17, M);
	for (int i = 0; i < 17; ++i) {
		EXPECT_EQ(out[i].x, xyz[i].x);
		EXPECT_EQ(out[i].y, xyz[i].y);
		EXPECT_EQ(out[i].z, xyz[i].z);
	}

	md_f32x3_t one = {1,1,1};
	one = md_f32x3x3_mul_f32x3(M, one);
	EXPECT_EQ(8, one.x);
	EXPECT_EQ(9, one.y);
	EXPECT_EQ(6, one.z);
}

UTEST(vec_math, mat4) {
	md_f32x4x4_t M = {
		1,2,3,4,
		4,3,2,1,
		3,4,1,2,
		4,3,2,1,
	};

	md_f32x3_t xyz[17], out[17];
	for (int i = 0; i < 17; ++i) xyz[i] = md_f32x3_set(1.0f + i, 2.0f - 0.5f * i, 0.25f * i);

	md_f32x4x4_batch_transform(out, xyz, 1.0f, 17, M);
	for (int i = 0; i < 17; ++i) {
		const md_f32x4_t r = md_f32x4x4_mul_f32x4(M, md_f32x4_from_f32x3(xyz[i], 1.0f));
		EXPECT_NEAR(r.x, out[i].x, 1.0e-4f);
		EXPECT_NEAR(r.y, out[i].y, 1.0e-4f);
		EXPECT_NEAR(r.z, out[i].z, 1.0e-4f);
	}
	md_f32x4x4_batch_transform_inplace(xyz, 1.0f, 17, M);
	for (int i = 0; i < 17; ++i) {
		EXPECT_EQ(out[i].x, xyz[i].x);
		EXPECT_EQ(out[i].y, xyz[i].y);
		EXPECT_EQ(out[i].z, xyz[i].z);
	}
}

UTEST(vec_math, translate) {
	md_f32x3_t xyz[19], out[19];
	for (int i = 0; i < 19; ++i) xyz[i] = md_f32x3_set((float)i, 100.0f + i, -3.0f * i);
	const md_f32x3_t t = {0.5f, -2.0f, 7.0f};
	md_f32x3_batch_translate(out, xyz, 19, t);
	for (int i = 0; i < 19; ++i) {
		EXPECT_EQ(xyz[i].x + t.x, out[i].x);
		EXPECT_EQ(xyz[i].y + t.y, out[i].y);
		EXPECT_EQ(xyz[i].z + t.z, out[i].z);
	}
	md_f32x3_batch_translate_inplace(xyz, 19, t);
	for (int i = 0; i < 19; ++i) {
		EXPECT_EQ(out[i].x, xyz[i].x);
		EXPECT_EQ(out[i].y, xyz[i].y);
		EXPECT_EQ(out[i].z, xyz[i].z);
	}
}

// md_f32x4x3_t is THREE columns of four rows. The product combines those three columns using v.x, v.y
// and v.z only: v.w takes no part, and there is no fourth column to reach for. The scalar fallback
// used to read M.elem[3][..], which is past the end of the struct, so this pins whichever path the
// build selects to the same reference computation.
UTEST(vec_math, md_f32x4x3_mul_f32x4) {
	md_f32x4x3_t M;
	M.col[0] = md_f32x4_set( 1.0f,  2.0f,  3.0f,  4.0f);
	M.col[1] = md_f32x4_set( 5.0f,  6.0f,  7.0f,  8.0f);
	M.col[2] = md_f32x4_set( 9.0f, 10.0f, 11.0f, 12.0f);   // a non zero w row, so a path that drops it is caught

	const md_f32x4_t v = md_f32x4_set(0.5f, -1.5f, 2.0f, 1000.0f);   // v.w is large on purpose: it must not leak in
	const md_f32x4_t r = md_f32x4x3_mul_f32x4(M, v);

	for (int row = 0; row < 4; ++row) {
		const float expect = M.elem[0][row] * v.x + M.elem[1][row] * v.y + M.elem[2][row] * v.z;
		EXPECT_NEAR(r.elem[row], expect, 1.0e-5f);
	}

	// changing v.w must change nothing at all
	const md_f32x4_t r2 = md_f32x4x3_mul_f32x4(M, md_f32x4_set(v.x, v.y, v.z, -7.0f));
	for (int row = 0; row < 4; ++row) {
		EXPECT_NEAR(r2.elem[row], r.elem[row], 1.0e-6f);
	}

	// built from a mat3, it has to agree with md_f32x3x3_mul_f32x3 on xyz and leave w at zero
	md_f32x3x3_t A = {
		1,2,3,
		4,5,6,
		7,8,9,
	};
	const md_f32x4x3_t B = md_f32x4x3_from_f32x3x3(A);
	const md_f32x3_t   e = md_f32x3x3_mul_f32x3(A, md_f32x3_set(v.x, v.y, v.z));
	const md_f32x4_t   g = md_f32x4x3_mul_f32x4(B, v);
	EXPECT_NEAR(g.x, e.x, 1.0e-5f);
	EXPECT_NEAR(g.y, e.y, 1.0e-5f);
	EXPECT_NEAR(g.z, e.z, 1.0e-5f);
	EXPECT_NEAR(g.w, 0.0f, 1.0e-6f);
}

// md_f32x4_cross is the cross product of the xyz, as md_f32x3_cross gives it, whatever the w of its operands, and its own w is
// exactly 0: a w that is not representable as a product (1/3) would leave the rounding of a.w * b.w there were it not
// cleared, and md_f32x4_dot and md_f32x4_length take w into account.
UTEST(vec_math, md_f32x4_cross) {
	const md_f32x3_t a3[4] = { {1, 0, 0}, {0.3f, -1.7f, 2.5f}, {-4.0f, 0.125f, 9.0f}, {1.0e-3f, 2.0e3f, -7.0f} };
	const md_f32x3_t b3[4] = { {0, 1, 0}, {2.2f, 0.4f, -1.1f}, {3.5f, -6.0f, 0.75f}, {5.0f, -1.0e-2f, 3.0e2f} };
	const float  w[3]  = { 0.0f, 1.0f, 1.0f / 3.0f };
	for (int i = 0; i < 4; ++i) {
		const md_f32x3_t e = md_f32x3_cross(a3[i], b3[i]);
		const float tol = 1.0e-6f * (1.0f + md_f32x3_length(a3[i]) * md_f32x3_length(b3[i]));
		for (int k = 0; k < 3; ++k) {
			const md_f32x4_t c = md_f32x4_cross(md_f32x4_from_f32x3(a3[i], w[k]), md_f32x4_from_f32x3(b3[i], -w[k] * 7.0f));
			EXPECT_NEAR(c.x, e.x, tol);
			EXPECT_NEAR(c.y, e.y, tol);
			EXPECT_NEAR(c.z, e.z, tol);
			EXPECT_EQ(0.0f, c.w);
		}
	}
	// Right handed: x cross y is z
	const md_f32x4_t z = md_f32x4_cross(md_f32x4_set(1, 0, 0, 0), md_f32x4_set(0, 1, 0, 0));
	EXPECT_EQ(0.0f, z.x);
	EXPECT_EQ(0.0f, z.y);
	EXPECT_EQ(1.0f, z.z);
}
