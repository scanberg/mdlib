#include "utest.h"

#include <md_cube.h>
#include <core/md_allocator.h>

UTEST(cube, read_benzene) {
	md_cube_t cube = {0};
	bool result = md_cube_file_load(&cube, STR(MD_UNITTEST_DATA_DIR "/benzene-pot.cube"), default_allocator);
	ASSERT_TRUE(result);

	EXPECT_NEAR(-9.053142, cube.origin[0], 1e-5);
	EXPECT_NEAR(-9.077350, cube.origin[1], 1e-5);
	EXPECT_NEAR(-6.512871, cube.origin[2], 1e-5);

	EXPECT_EQ(12, cube.atom.count);
	EXPECT_EQ(89, cube.data.num_x);
	EXPECT_EQ(89, cube.data.num_y);
	EXPECT_EQ(64, cube.data.num_z);
	EXPECT_EQ(1,  cube.data.num_m);

	EXPECT_EQ(6, cube.atom.number[0]);
	EXPECT_NEAR( 6.000000, cube.atom.charge[0],		1e-5);
	EXPECT_NEAR(-1.902215, cube.atom.coord[0][0],	1e-5);
	EXPECT_NEAR( 0.000005, cube.atom.coord[0][2],	1e-5);
	
	EXPECT_EQ(1, cube.atom.number[11]);
	EXPECT_NEAR( 1.000000, cube.atom.charge[11],	1e-5);
	EXPECT_NEAR(-4.514684, cube.atom.coord[11][0],	1e-5);
	EXPECT_NEAR(-0.000106, cube.atom.coord[11][2],	1e-5);
	
	EXPECT_NEAR(3.37638E-04, cube.data.val[0],		1e-5);
	EXPECT_NEAR(3.51617E-04, cube.data.val[89*89*64*1 - 1], 1e-5);

	EXPECT_STREQ(" Title Card Required potential", str_ptr(cube.title));
	EXPECT_STREQ(" Electrostatic potential from Total SCF Density", str_ptr(cube.comment));
	
	md_cube_free(&cube, default_allocator);
}

UTEST(cube, read_A2B2) {
	md_cube_t cube = {0};
	bool result = md_cube_file_load(&cube, STR(MD_UNITTEST_DATA_DIR "/A2B2_State2-GS.cube"), default_allocator);
	EXPECT_TRUE(result);

	EXPECT_NEAR(-16.311785, cube.origin[0], 1e-5);
	EXPECT_NEAR(-21.859915, cube.origin[1], 1e-5);
	EXPECT_NEAR(-21.637702, cube.origin[2], 1e-5);

	EXPECT_EQ(144, cube.atom.count);
	EXPECT_EQ(66, cube.data.num_x);
	EXPECT_EQ(88, cube.data.num_y);
	EXPECT_EQ(88, cube.data.num_z);
	EXPECT_EQ(1,  cube.data.num_m);

	EXPECT_EQ(8, cube.atom.number[0]);
	EXPECT_NEAR( 8.000000, cube.atom.charge[0],		1e-5);
	EXPECT_NEAR(-7.172985, cube.atom.coord[0][0],	1e-5);
	EXPECT_NEAR(11.547489, cube.atom.coord[0][2],	1e-5);
	
	EXPECT_EQ(1, cube.atom.number[143]);
	EXPECT_NEAR( 1.000000, cube.atom.charge[143],	1e-5);
	EXPECT_NEAR(-4.068669, cube.atom.coord[143][0],	1e-5);
	EXPECT_NEAR(-2.197243, cube.atom.coord[143][2],	1e-5);

	EXPECT_EQ(1, cube.data.num_m);
	EXPECT_EQ(1, cube.data.id[0]);

	EXPECT_NEAR(-2.02559E-11, cube.data.val[0],		1e-5);
	EXPECT_NEAR(5.31713E-10, cube.data.val[66*88*88*1 - 1], 1e-5);

	EXPECT_STREQ("     comment MO=411 ||   comment MO=412 ||   comment MO=409 ||   comment MO=410",  str_ptr(cube.title));
	EXPECT_STREQ("     MO coefficients scaled by    0.515986 +   MO coefficien +   MO coefficients", str_ptr(cube.comment));

	md_cube_free(&cube, default_allocator);
}

static bool cmp_cube_v3(md_cube_v3 a, md_cube_v3 b) {
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

UTEST(cube, serialize_deserialize) {
	md_cube_t cube_a = {0};
	str_t path = STR(MD_UNITTEST_DATA_DIR "/A2B2_State2-GS.cube");
	ASSERT_TRUE(md_cube_file_load(&cube_a, path, default_allocator));
	
	str_t str = md_cube_serialize(&cube_a, default_allocator);
    ASSERT_FALSE(str_empty(str));

	md_cube_t cube_b = {0};
	ASSERT_TRUE(md_cube_deserialize(&cube_b, str, default_allocator));

	EXPECT_STREQ(cube_a.title.ptr, cube_b.title.ptr);
	EXPECT_STREQ(cube_a.comment.ptr, cube_b.comment.ptr);

	EXPECT_TRUE(cmp_cube_v3(cube_a.origin, cube_b.origin));
	EXPECT_TRUE(cmp_cube_v3(cube_a.xaxis, cube_b.xaxis));
	EXPECT_TRUE(cmp_cube_v3(cube_a.yaxis, cube_b.yaxis));
	EXPECT_TRUE(cmp_cube_v3(cube_a.zaxis, cube_b.zaxis));

	EXPECT_EQ(cube_a.data.num_x, cube_b.data.num_x);
	EXPECT_EQ(cube_a.data.num_y, cube_b.data.num_y);
	EXPECT_EQ(cube_a.data.num_z, cube_b.data.num_z);
	EXPECT_EQ(cube_a.data.num_m, cube_b.data.num_m);

	if (cube_a.data.id) {
		EXPECT_NE(cube_b.data.id, NULL);
		int size_a = cube_a.data.num_x * cube_a.data.num_y * cube_a.data.num_z * cube_a.data.num_m;
		int size_b = cube_b.data.num_x * cube_b.data.num_y * cube_b.data.num_z * cube_b.data.num_m;
		ASSERT_EQ(size_a, size_b);
		for (int i = 0; i < size_a; ++i) {
			EXPECT_NEAR(cube_a.data.val[i], cube_b.data.val[i], 1.0e-9);
		}
	}

	ASSERT_EQ(cube_a.atom.count, cube_b.atom.count);
	for (int i = 0; i < cube_a.atom.count; ++i) {
        EXPECT_EQ(cube_a.atom.number[i], cube_b.atom.number[i]);
        EXPECT_NEAR(cube_a.atom.charge[i], cube_b.atom.charge[i], 1.0e-9);
        EXPECT_TRUE(cmp_cube_v3(cube_a.atom.coord[i], cube_b.atom.coord[i]));
	}

	str_free(str, default_allocator);
	md_cube_free(&cube_a, default_allocator);
	md_cube_free(&cube_b, default_allocator);
}
