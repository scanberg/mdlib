#include "utest.h"

#include "qm_test_util.h"

#include <md_cube.h>
#include <md_io.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_unit.h>
#include <core/md_file.h>

#include <stdio.h>

UTEST(cube, read_benzene) {
	md_cube_t cube = {0};
	bool result = md_cube_file_load(&cube, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/benzene-pot.cube"), md_get_heap_allocator());
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

	EXPECT_STREQ(" Title Card Required potential", md_str_ptr(cube.title));
	EXPECT_STREQ(" Electrostatic potential from Total SCF Density", md_str_ptr(cube.comment));
	
	md_cube_free(&cube, md_get_heap_allocator());
}

UTEST(cube, read_A2B2) {
	md_cube_t cube = {0};
	bool result = md_cube_file_load(&cube, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/A2B2_State2-GS.cube"), md_get_heap_allocator());
	ASSERT_TRUE(result);

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

	ASSERT_TRUE(cube.data.val);
	EXPECT_NEAR(-2.02559E-11, cube.data.val[0],		1e-5);
	EXPECT_NEAR(5.31713E-10, cube.data.val[66*88*88*1 - 1], 1e-5);

	EXPECT_STREQ("     comment MO=411 ||   comment MO=412 ||   comment MO=409 ||   comment MO=410",  md_str_ptr(cube.title));
	EXPECT_STREQ("     MO coefficients scaled by    0.515986 +   MO coefficien +   MO coefficients", md_str_ptr(cube.comment));

	md_cube_free(&cube, md_get_heap_allocator());
}

UTEST(cube, nonexistent_file) {
	md_cube_t cube = {0};
	bool result = md_cube_file_load(&cube, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.cube"), md_get_heap_allocator());
	EXPECT_FALSE(result);
	
	// Should be safe to free even when loading failed
	md_cube_free(&cube, md_get_heap_allocator());
}

UTEST(cube, empty_path) {
	md_cube_t cube = {0};
	md_str_t empty_path = {0,0};
	bool result = md_cube_file_load(&cube, empty_path, md_get_heap_allocator());
	EXPECT_FALSE(result);
	
	md_cube_free(&cube, md_get_heap_allocator());
}

// A negative number of points means the header's lengths are in Angstrom (Gaussian's convention, cubegen writes it);
// they are kept in bohr. The count was taken as given and the cube refused as one of no points.
UTEST(cube, angstrom_header) {
    const md_str_t text = MD_STR_LIT(
        "title\n"
        "comment\n"
        "    1    1.000000    0.000000    0.000000\n"
        "   -2    0.529177    0.000000    0.000000\n"
        "   -1    0.000000    0.529177    0.000000\n"
        "   -1    0.000000    0.000000    0.529177\n"
        "    8    8.000000    0.529177    0.000000    0.000000\n"
        "  1.0E-01  2.0E-01\n");
    md_cube_t cube = {0};
    ASSERT_TRUE(md_cube_deserialize(&cube, text, md_get_heap_allocator()));
    EXPECT_EQ(2, cube.data.num_x);
    EXPECT_EQ(1, cube.data.num_y);
    EXPECT_NEAR(1.0 / 0.52917721092, cube.origin[0], 1.0e-4);
    EXPECT_NEAR(1.0, cube.xaxis[0], 1.0e-5);
    EXPECT_NEAR(1.0, cube.atom.coord[0][0], 1.0e-5);
    EXPECT_NEAR(0.2f, cube.data.val[1], 1.0e-7);
    md_cube_free(&cube, md_get_heap_allocator());
}

// A grid without atoms, and a grid of two values per point (NVAL 2): the first was refused for having no atoms.
UTEST(cube, no_atoms_and_values_per_point) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const md_str_t text = MD_STR_LIT(
        "no atoms\n"
        "two values per point\n"
        "    0    0.000000    0.000000    0.000000    2\n"
        "    2    0.500000    0.000000    0.000000\n"
        "    1    0.000000    0.500000    0.000000\n"
        "    2    0.000000    0.000000    0.500000\n"
        "  1.0E+00  2.0E+00  3.0E+00  4.0E+00  5.0E+00  6.0E+00\n"
        "  7.0E+00  8.0E+00\n");
    md_cube_t a = {0};
    ASSERT_TRUE(md_cube_deserialize(&a, text, alloc));
    EXPECT_EQ(0, a.atom.count);
    EXPECT_EQ(2, a.data.num_m);
    EXPECT_EQ(8.0f, a.data.val[7]);
    md_cube_free(&a, alloc);
}

// A cube as a system: its atoms in Å, and its grid published as a script volume is, its origin the corner of the first
// voxel
UTEST(cube, system_init) {
    qm_test_t t;
    qm_test_init(&t, MD_MEGABYTES(64));
    const md_str_t path = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/benzene-pot.cube");
    EXPECT_EQ(MD_IO_FORMAT_CUBE, md_io_format_from_file(path));
    ASSERT_TRUE(md_cube_system_init_from_file(&t.sys, &t.state, path, NULL));
    ASSERT_EQ(12u, t.sys.atom.count);
    EXPECT_EQ(6, md_atom_atomic_number(&t.sys.atom, 0));
    EXPECT_NEAR(-1.902215 * 0.52917721092, t.state.xyz[0].x, 1e-5);

    const md_attribute_t* v = qm_test_attr(&t, MD_STR_LIT("cube/benzene-pot"));
    ASSERT_TRUE(v);
    EXPECT_EQ(MD_ATTRIBUTE_TYPE_F32, v->format.type);
    EXPECT_EQ(3u, v->format.rank);
    EXPECT_EQ(89u, v->format.shape[0]);
    EXPECT_EQ(64u, v->format.shape[2]);
    EXPECT_NEAR(3.37638E-04, ((const float*)v->data)[0], 1e-9);

    double o[3], b[9], xyz[36], idx[12];
    const md_attribute_t* ao = qm_test_attr(&t, MD_STR_LIT("cube/benzene-pot/origin"));
    const md_attribute_t* ab = qm_test_attr(&t, MD_STR_LIT("cube/benzene-pot/basis"));
    ASSERT_TRUE(ao && ab);
    ASSERT_EQ(3u, md_attribute_extract_f64(o, 3, &t.sys.attributes, ao, md_attribute_slice_all(), md_unit_angstrom()));
    ASSERT_EQ(9u, md_attribute_extract_f64(b, 9, &t.sys.attributes, ab, md_attribute_slice_all(), md_unit_angstrom()));
    // The first point less half a step along each axis
    EXPECT_NEAR(-9.053142 * 0.52917721092 - 0.5 * (b[0] + b[3] + b[6]), o[0], 1e-5);
    EXPECT_NEAR(0.0, b[1], 1e-9);

    const md_attribute_t* ac = qm_test_attr(&t, MD_STR_LIT("cube/benzene-pot/atom/coordinate"));
    const md_attribute_t* ai = qm_test_attr(&t, MD_STR_LIT("cube/benzene-pot/atom/system_index"));
    ASSERT_TRUE(ac && ai);
    ASSERT_EQ(36u, md_attribute_extract_f64(xyz, 36, &t.sys.attributes, ac, md_attribute_slice_all(), md_unit_angstrom()));
    EXPECT_NEAR(t.state.xyz[11].x, xyz[33], 1e-5);
    ASSERT_EQ(12u, md_attribute_extract_f64(idx, 12, &t.sys.attributes, ai, md_attribute_slice_all(), md_unit_none()));
    EXPECT_EQ(11.0, idx[11]);
    qm_test_free(&t);
}

// A cube of the loaded system's atoms supplements it; another's does not. A cube of several values a point publishes
// them outermost.
UTEST(cube, supplement) {
    qm_test_t t;
    qm_test_init(&t, MD_MEGABYTES(64));
    const md_str_t benzene = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/benzene-pot.cube");
    const md_str_t a2b2    = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/A2B2_State2-GS.cube");
    ASSERT_TRUE(md_cube_system_init_from_file(&t.sys, &t.state, benzene, NULL));
    EXPECT_TRUE(md_cube_system_is_file_supplemental(&t.sys, benzene));
    EXPECT_FALSE(md_cube_system_is_file_supplemental(&t.sys, a2b2));
    EXPECT_FALSE(md_cube_system_supplement_from_file(&t.sys, a2b2, NULL));
    EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("cube/A2B2_State2-GS")));
    ASSERT_TRUE(md_cube_system_supplement_from_file(&t.sys, benzene, NULL));
    EXPECT_TRUE(qm_test_has(&t, MD_STR_LIT("cube/benzene-pot")));
    qm_test_free(&t);

    qm_test_init(&t, MD_MEGABYTES(4));
    const md_str_t text = MD_STR_LIT(
        "two values per point\n"
        "\n"
        "   -1    0.000000    0.000000    0.000000\n"
        "    2    0.500000    0.000000    0.000000\n"
        "    1    0.000000    0.500000    0.000000\n"
        "    2    0.000000    0.000000    0.500000\n"
        "    8    8.000000    0.000000    0.000000    0.000000\n"
        "    2   11   12\n"
        "  1.0E+00  2.0E+00  3.0E+00  4.0E+00  5.0E+00  6.0E+00\n"
        "  7.0E+00  8.0E+00\n");
    const char* tmp = "cube_test_two.cube";
    md_file_t f = {0};
    ASSERT_TRUE(md_file_open(&f, md_str_from_cstr(tmp), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    md_file_write(f, text.ptr, text.len);
    md_file_close(&f);
    ASSERT_TRUE(md_cube_system_init_from_file(&t.sys, &t.state, md_str_from_cstr(tmp), NULL));
    const md_attribute_t* v = qm_test_attr(&t, MD_STR_LIT("cube/cube_test_two"));
    ASSERT_TRUE(v);
    ASSERT_EQ(4u, v->format.rank);
    EXPECT_EQ(2u, v->format.shape[0]);
    const float* d = (const float*)v->data;
    // [x][y][z][m] = 1..8: of the first value 1 3 5 7, of the second 2 4 6 8
    EXPECT_EQ(3.0f, d[1]);
    EXPECT_EQ(2.0f, d[4]);
    double id[2];
    const md_attribute_t* ai = qm_test_attr(&t, MD_STR_LIT("cube/cube_test_two/id"));
    ASSERT_TRUE(ai);
    ASSERT_EQ(2u, md_attribute_extract_f64(id, 2, &t.sys.attributes, ai, md_attribute_slice_all(), md_unit_none()));
    EXPECT_EQ(12.0, id[1]);
    remove(tmp);
    qm_test_free(&t);
}
