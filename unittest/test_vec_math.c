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
