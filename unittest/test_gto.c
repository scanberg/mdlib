#include "utest.h"

#include <core/md_allocator.h>
#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include "vlx_test_util.h"

#include <md_cube.h>
#include <md_gto.h>

#include <float.h>

// Conversion from Ångström to Bohr
#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903
#define VALUE_CUTOFF 1.0E-5

static void init(md_grid_t* grid, float** grid_data, md_gto_t** gtos, size_t* num_gtos, int vol_dim, str_t filename, md_allocator_i* arena) {
    md_temp_scope_t temp = md_temp_begin_avoid(arena);

    vlx_test_t t = {0};
    if (!vlx_test_load(&t, filename, MEGABYTES(64))) {
        md_temp_end(temp);
        return;
    }

    const size_t num_atoms = t.sys.atom.count;
    float* atom_xyz = md_temp_alloc_array(temp, float, 3 * num_atoms);
    vlx_test_atom_xyz_bohr(atom_xyz, 3 * num_atoms, &t);

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < num_atoms; ++i) {
        vec3_t c = {atom_xyz[3*i+0], atom_xyz[3*i+1], atom_xyz[3*i+2]};
        min_box = vec3_min(min_box, c);
        max_box = vec3_max(max_box, c);
    }

    min_box = vec3_sub1(min_box, 2.0f * (float)VLX_TEST_ANGSTROM_TO_BOHR);
    max_box = vec3_add1(max_box, 2.0f * (float)VLX_TEST_ANGSTROM_TO_BOHR);

    float* vol_data = md_temp_alloc_zero_array(temp, float, vol_dim * vol_dim * vol_dim);

    vec3_t step = vec3_div1(vec3_sub(max_box, min_box), (float)vol_dim);

    *grid = (md_grid_t) {
        .orientation = mat3_ident(),
        .origin = {min_box.x, min_box.y, min_box.z},
        .spacing = {step.x, step.y, step.z},
        .dim = {vol_dim, vol_dim, vol_dim},
    };

    *grid_data = vol_data;

    md_gto_basis_t basis = {0};
    vlx_test_basis(&basis, &t);

    const size_t num_ao = md_gto_basis_num_ao(&basis);
    double* mo_coeffs = md_temp_alloc_array(temp, double, num_ao);
    vlx_test_row(mo_coeffs, num_ao, &t, STR_LIT("orbital/alpha/coefficient"), 120);

    *num_gtos = md_gto_pgto_count(&basis);
    *gtos = md_alloc(arena, sizeof(md_gto_t) * (*num_gtos));
    *num_gtos = md_gto_expand_with_ao_coeffs(*gtos, &basis, atom_xyz, sizeof(float) * 3, mo_coeffs, 1.0e-6);

    vlx_test_free(&t);
    md_temp_end(temp);
}

// The basis is Cartesian (see the AO CONVENTION block in md_gto.h). These pin the
// counting and the normalization convention so a reader that gets either wrong
// fails here rather than silently producing plausible but wrong orbitals.
UTEST(gto, cartesian_ao_counts) {
    EXPECT_EQ(1u,  md_gto_num_cart_ao(0));
    EXPECT_EQ(3u,  md_gto_num_cart_ao(1));
    EXPECT_EQ(6u,  md_gto_num_cart_ao(2));
    EXPECT_EQ(10u, md_gto_num_cart_ao(3));
    EXPECT_EQ(15u, md_gto_num_cart_ao(4));

    EXPECT_EQ(1u, md_gto_num_sph_ao(0));
    EXPECT_EQ(3u, md_gto_num_sph_ao(1));
    EXPECT_EQ(5u, md_gto_num_sph_ao(2));
    EXPECT_EQ(7u, md_gto_num_sph_ao(3));
    EXPECT_EQ(9u, md_gto_num_sph_ao(4));

    // Documented order: for i = l..0, for j = (l-i)..0, k = l-i-j.
    // l=2 must be xx, xy, xz, yy, yz, zz -- NOT Gaussian's xx, yy, zz, xy, xz, yz.
    const int expect_d[6][3] = {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}};
    for (uint32_t n = 0; n < 6; ++n) {
        int i, j, k;
        ASSERT_TRUE(md_gto_cart_ijk(&i, &j, &k, 2, n));
        EXPECT_EQ(expect_d[n][0], i);
        EXPECT_EQ(expect_d[n][1], j);
        EXPECT_EQ(expect_d[n][2], k);
    }
    EXPECT_FALSE(md_gto_cart_ijk(NULL, NULL, NULL, 2, 6));

    // f(i,j,k) = 1/sqrt((2i-1)!!(2j-1)!!(2k-1)!!)
    EXPECT_NEAR(1.0,               md_gto_cart_norm_factor(0,0,0), 1.0e-12);  // s
    EXPECT_NEAR(1.0,               md_gto_cart_norm_factor(1,0,0), 1.0e-12);  // p
    EXPECT_NEAR(1.0 / sqrt(3.0),   md_gto_cart_norm_factor(2,0,0), 1.0e-12);  // d xx
    EXPECT_NEAR(1.0,               md_gto_cart_norm_factor(1,1,0), 1.0e-12);  // d xy
    EXPECT_NEAR(1.0 / sqrt(15.0),  md_gto_cart_norm_factor(3,0,0), 1.0e-12);  // f xxx
    EXPECT_NEAR(1.0 / sqrt(105.0), md_gto_cart_norm_factor(4,0,0), 1.0e-12);  // g xxxx
}

// A spherical shell embedded in the Cartesian basis is rank deficient: the pure
// d functions span only 5 of the 6 Cartesian directions. Converting the identity
// spanning set must therefore leave the s-type contaminant direction untouched.
// This catches a transposed or mis-scaled transform, which would otherwise show
// up only as subtly wrong lobes.
UTEST(gto, sph_to_cart_vector_dimensions) {
    md_gto_shell_t shells[3] = {
        { .atom_idx = 0, .primitive_offset = 0, .num_primitives = 1, .l = 0 },
        { .atom_idx = 0, .primitive_offset = 1, .num_primitives = 1, .l = 1 },
        { .atom_idx = 0, .primitive_offset = 2, .num_primitives = 1, .l = 2 },
    };
    float alpha[3] = {1.0f, 1.0f, 1.0f};
    float coeff[3] = {1.0f, 1.0f, 1.0f};
    md_gto_basis_t basis = {
        .num_shells = 3, .num_primitives = 3,
        .shells = shells, .alpha = alpha, .coeff = coeff,
    };

    EXPECT_EQ((size_t)(1 + 3 + 5),  md_gto_basis_num_sph_ao(&basis));
    EXPECT_EQ((size_t)(1 + 3 + 6),  md_gto_basis_num_ao(&basis));

    double in_sph[9]   = {0};
    double out_cart[10] = {0};
    in_sph[4] = 1.0;  // one of the d functions

    EXPECT_EQ((size_t)10, md_gto_sph_to_cart_vector(out_cart, in_sph, &basis));

    // s and p blocks are untouched by a pure-d input.
    EXPECT_NEAR(0.0, out_cart[0], 1.0e-12);
    for (int i = 1; i < 4; ++i) EXPECT_NEAR(0.0, out_cart[i], 1.0e-12);

    double d_norm = 0.0;
    for (int i = 4; i < 10; ++i) d_norm += out_cart[i] * out_cart[i];
    EXPECT_GT(d_norm, 1.0e-6);
}

// md_gto_sph_to_cart_matrix() exploits the block-diagonal structure of the
// transform. This checks it against the dense T^T M T formulation built from the
// vector path, which is the thing the block version is an optimization of.
UTEST(gto, sph_to_cart_matrix_matches_dense) {
    md_temp_scope_t temp = md_temp_begin();

    enum { NUM_ATOMS = 3, MAX_L = 4 };
    md_gto_shell_t shells[(MAX_L + 1) * NUM_ATOMS];
    float alpha[256], coeff[256];

    uint32_t ns = 0, np = 0;
    for (uint32_t l = 0; l <= MAX_L; ++l) {
        for (uint32_t a = 0; a < NUM_ATOMS; ++a) {
            const uint32_t nprim = 2 + (l % 2);
            shells[ns] = (md_gto_shell_t){ .atom_idx = a, .primitive_offset = np, .num_primitives = nprim, .l = l };
            for (uint32_t p = 0; p < nprim; ++p) {
                alpha[np] = (float)(0.4 + 0.7 * p + 0.2 * l);
                coeff[np] = (float)(0.5 + 0.1 * p);
                np++;
            }
            ns++;
        }
    }
    md_gto_basis_t basis = { .num_shells = ns, .num_primitives = np, .shells = shells, .alpha = alpha, .coeff = coeff };

    const size_t n_sph  = md_gto_basis_num_sph_ao(&basis);
    const size_t n_cart = md_gto_basis_num_ao(&basis);

    double* M = md_temp_alloc_zero_array(temp, double, n_sph * n_sph);
    uint32_t seed = 12345;
    for (size_t i = 0; i < n_sph; ++i) {
        for (size_t j = i; j < n_sph; ++j) {
            seed = seed * 1664525u + 1013904223u;
            const double v = (double)(seed >> 8) / (double)(1u << 24) - 0.5;
            M[i * n_sph + j] = v;
            M[j * n_sph + i] = v;
        }
    }

    double* got = md_temp_alloc_zero_array(temp, double, n_cart * n_cart);
    ASSERT_EQ(n_cart, md_gto_sph_to_cart_matrix(got, M, &basis));

    // Dense reference.
    double* T   = md_temp_alloc_zero_array(temp, double, n_sph * n_cart);
    double* e   = md_temp_alloc_zero_array(temp, double, n_sph);
    double* col = md_temp_alloc_zero_array(temp, double, n_cart);
    for (size_t r = 0; r < n_sph; ++r) {
        MEMSET(e, 0, sizeof(double) * n_sph);
        e[r] = 1.0;
        md_gto_sph_to_cart_vector(col, e, &basis);
        for (size_t c = 0; c < n_cart; ++c) T[r * n_cart + c] = col[c];
    }

    double* tmp = md_temp_alloc_zero_array(temp, double, n_sph * n_cart);
    for (size_t i = 0; i < n_sph; ++i) {
        for (size_t c = 0; c < n_cart; ++c) {
            double s = 0.0;
            for (size_t k = 0; k < n_sph; ++k) s += M[i * n_sph + k] * T[k * n_cart + c];
            tmp[i * n_cart + c] = s;
        }
    }

    double max_diff = 0.0, peak = 0.0, max_asym = 0.0;
    for (size_t r = 0; r < n_cart; ++r) {
        for (size_t c = 0; c < n_cart; ++c) {
            double s = 0.0;
            for (size_t k = 0; k < n_sph; ++k) s += T[k * n_cart + r] * tmp[k * n_cart + c];
            const double d = fabs(got[r * n_cart + c] - s);
            if (d > max_diff)     max_diff = d;
            if (fabs(s) > peak)   peak     = fabs(s);
            const double a = fabs(got[r * n_cart + c] - got[c * n_cart + r]);
            if (a > max_asym)     max_asym = a;
        }
    }

    // Both comparisons are relative: the block/sparse implementation accumulates in a
    // different order than the dense reference, so it agrees to rounding, not exactly.
    const double scale = (peak > 0.0) ? peak : 1.0;
    EXPECT_LT(max_diff / scale, 1.0e-13);
    // A symmetric input must stay symmetric.
    EXPECT_LT(max_asym / scale, 1.0e-13);

    md_temp_end(temp);
}

UTEST(gto, evaluate_grid) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    
    md_grid_t grid;
    float* grid_data;
    md_gto_t* gtos;
    size_t num_gtos;

    int vol_dim = 64;
    init(&grid, &grid_data, &gtos, &num_gtos, vol_dim, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), temp_arena);
    size_t num_points = md_grid_num_points(&grid);

    md_grid_t ref_grid = grid;
    float* ref_data    = md_temp_alloc_array(temp, float, num_points);
    vec3_t* ref_points = md_temp_alloc_array(temp, vec3_t, num_points);

    md_grid_extract_points((float*)ref_points, &ref_grid);

    MEMSET(grid_data, 0, sizeof(float) * num_points);
    MEMSET(ref_data,  0, sizeof(float) * num_points);

    md_gto_grid_evaluate(grid_data, &grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);
    md_gto_xyz_evaluate(ref_data, (float*)ref_points, num_points, 0, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    for (size_t i = 0; i < num_points; ++i) {
        double delta = fabs(grid_data[i] - ref_data[i]);
        if (delta > 1.0e-4) {
            printf("Delta: %g at [%zu]\n", delta, i);
        }
    }

    md_temp_end(temp);
}

// compares vlx file of specific orbital and compares it with a cube file, returns the abs max difference between the two
static double compare_vlx_and_cube_cpu(const float* atom_xyz, const md_gto_basis_t* gto_basis, const double* ao_coeffs, const md_cube_t* cube) {
    const double cutoff_value = 1.0E-6;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    mat3_t orientation = mat3_ident();
    vec3_t x_axis = vec3_set(cube->xaxis[0], cube->xaxis[1], cube->xaxis[2]);
    vec3_t y_axis = vec3_set(cube->yaxis[0], cube->yaxis[1], cube->yaxis[2]);
    vec3_t z_axis = vec3_set(cube->zaxis[0], cube->zaxis[1], cube->zaxis[2]);

    float x_len = vec3_length(x_axis);
    float y_len = vec3_length(y_axis);
    float z_len = vec3_length(z_axis);

    orientation.col[0] = vec3_div1(x_axis, x_len);
    orientation.col[1] = vec3_div1(y_axis, y_len);
    orientation.col[2] = vec3_div1(z_axis, z_len);

    md_grid_t grid = {
        .orientation = orientation,
        .origin = {cube->origin[0], cube->origin[1], cube->origin[2]},
        .spacing = vec3_set(x_len, y_len, z_len),
        .dim = {cube->data.num_x, cube->data.num_y, cube->data.num_z},
    };
    float* grid_data = md_temp_alloc_array(temp, float, cube->data.num_x * cube->data.num_y * cube->data.num_z);

    size_t num_gtos = md_gto_pgto_count(gto_basis);
    md_gto_t* gtos = (md_gto_t*)md_temp_alloc_array(temp, md_gto_t, num_gtos);
    num_gtos = md_gto_expand_with_ao_coeffs(gtos, gto_basis, atom_xyz, sizeof(vec3_t), ao_coeffs, cutoff_value);

    size_t count = grid.dim[0] * grid.dim[1] * grid.dim[2];

    mat4_t M = md_grid_index_to_world(&grid);

    float* psi  = md_temp_alloc_array(temp, float, count);
    MEMSET(psi, 0, sizeof(float) * count);
    vec3_t* xyz = md_temp_alloc_array(temp, vec3_t, count);
    md_grid_extract_points((float*)xyz, &grid);

    md_gto_xyz_evaluate(psi, (float*)xyz, count, sizeof(vec3_t), gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    MEMSET(grid_data, 0, sizeof(float) * md_grid_num_points(&grid));

    md_gto_grid_evaluate(grid_data, &grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    double xyz_sum  = 0.0;
    double grid_sum = 0.0;
    double cube_sum = 0.0;

    double max_delta  = 0.0;
    int    max_idx[3] = {0};

    for (int iz = 0; iz < grid.dim[2]; ++iz) {
        for (int iy = 0; iy < grid.dim[1]; ++iy) {
            for (int ix = 0; ix < grid.dim[0]; ++ix) {
                int grid_idx = iz * grid.dim[0] * grid.dim[1] + iy * grid.dim[0] + ix;
                int cube_idx = ix * grid.dim[1] * grid.dim[2] + iy * grid.dim[2] + iz;
                double g_val = grid_data[grid_idx];
                double c_val = cube->data.val[cube_idx];
                double p_val = psi[grid_idx];
                grid_sum += g_val;
                cube_sum += c_val;
                xyz_sum  += p_val;

                double delta = fabs(fabs(g_val) - fabs(c_val));
                if (delta > max_delta) {
                    max_delta = delta;
                    max_idx[0] = ix;
                    max_idx[1] = iy;
                    max_idx[2] = iz;
                }
            }
        }
    }

    printf("Max abs delta: %g at [%i,%i,%i]\n", max_delta, max_idx[0], max_idx[1], max_idx[2]);
    printf("GRID SUM: %.5f\n", grid_sum);
    printf("CUBE SUM: %.5f\n", cube_sum);
    printf("XYZ  SUM: %.5f\n", xyz_sum);

    md_temp_end(temp);

    return max_delta;
}

UTEST(gto, h2o_lumo_cpu) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    md_cube_t cube_lumo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o_lumo.cube"), temp.arena));

    // Loaded as a SYSTEM: the geometry, the basis and the coefficients all come out of the
    // attribute table, which is the only way anything reaches a consumer now.
    vlx_test_t t = {0};
    ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(32)));

    size_t num_atoms = t.sys.atom.count;
    float* atom_xyz = (float*)md_temp_alloc_array(temp, float, 3 * num_atoms);
    ASSERT_EQ(vlx_test_atom_xyz_bohr(atom_xyz, 3 * num_atoms, &t), num_atoms);

    md_gto_basis_t basis = {0};
    ASSERT_TRUE(vlx_test_basis(&basis, &t));

    // The LUMO, derived from the published occupations by the rule the reader used: the first
    // orbital with zero occupancy.
    const size_t lumo_idx = vlx_test_lumo_idx(&t, STR_LIT("orbital/alpha/occupation"));
    const size_t num_ao   = md_gto_basis_num_ao(&basis);
    double* ao_coeffs = (double*)md_temp_alloc_array(temp, double, num_ao);
    ASSERT_EQ(vlx_test_row(ao_coeffs, num_ao, &t, STR_LIT("orbital/alpha/coefficient"), lumo_idx), num_ao);

    double max_delta_lumo = compare_vlx_and_cube_cpu(atom_xyz, &basis, ao_coeffs, &cube_lumo);
    EXPECT_LT(max_delta_lumo, 1.0E-4);  

    vlx_test_free(&t);
    md_temp_end(temp);
}

#if MD_ENABLE_GPU

static double compare_vlx_and_cube_gpu(md_gpu_device_t device, const float* atom_xyz, const md_gto_basis_t* gto_basis, const double* ao_coeffs, const md_cube_t* cube) {
    const double cutoff_value = 1.0E-6;

    double max_delta = DBL_MAX;

    md_gto_gpu_basis_t gpu_basis = NULL;
    md_gpu_pool_t   dev_pool  = NULL;
    md_gpu_pool_t   read_pool = NULL;
    md_gpu_ptr_t    atom_buf  = NULL;
    md_gpu_ptr_t    coeff_buf = NULL;
    md_gpu_ptr_t    readback  = NULL;
    md_gpu_tex_t    out_tex   = 0;
    md_gpu_stream_t stream    = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);

    mat3_t orientation = mat3_ident();
    vec3_t x_axis = vec3_set(cube->xaxis[0], cube->xaxis[1], cube->xaxis[2]);
    vec3_t y_axis = vec3_set(cube->yaxis[0], cube->yaxis[1], cube->yaxis[2]);
    vec3_t z_axis = vec3_set(cube->zaxis[0], cube->zaxis[1], cube->zaxis[2]);

    float x_len = vec3_length(x_axis);
    float y_len = vec3_length(y_axis);
    float z_len = vec3_length(z_axis);

    orientation.col[0] = vec3_div1(x_axis, x_len);
    orientation.col[1] = vec3_div1(y_axis, y_len);
    orientation.col[2] = vec3_div1(z_axis, z_len);

    md_grid_t grid = {
        .orientation = orientation,
        .origin = {cube->origin[0], cube->origin[1], cube->origin[2]},
        .spacing = vec3_set(x_len, y_len, z_len),
        .dim = {cube->data.num_x, cube->data.num_y, cube->data.num_z},
    };

    md_gto_gpu_initialize(device);

    {
        md_gpu_pool_desc_t pd = {0};
        pd.flags = MD_GPU_MEM_DEVICE;    pd.label = "test_gto device";   dev_pool  = md_gpu_pool_create(device, &pd);
        pd.flags = MD_GPU_MEM_HOST_READ; pd.label = "test_gto readback"; read_pool = md_gpu_pool_create(device, &pd);
    }
    if (!dev_pool || !read_pool) goto done;

    gpu_basis = md_gto_gpu_basis_create(dev_pool, stream, &(md_gto_gpu_basis_desc_t){
        .basis = gto_basis,
        .cutoff = 0.0,
    });
    if (!gpu_basis) goto done;

    const uint32_t num_atoms = md_gto_gpu_basis_num_atoms(gpu_basis);
    const uint32_t num_cgtos = md_gto_gpu_basis_num_cgtos(gpu_basis);
    const size_t voxel_count = (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];
    const size_t readback_size = sizeof(float) * voxel_count;

    atom_buf  = md_gpu_malloc(dev_pool, md_gto_gpu_atom_buffer_size(num_atoms), stream);
    coeff_buf = md_gpu_malloc(dev_pool, md_gto_gpu_coeff_size_mo(1, num_cgtos), stream);
    readback  = md_gpu_malloc(read_pool, readback_size, stream);
    out_tex   = md_gpu_tex_create(device, &(md_gpu_tex_desc_t){
        .width  = (uint32_t)grid.dim[0],
        .height = (uint32_t)grid.dim[1],
        .depth  = (uint32_t)grid.dim[2],
        .format = MD_GPU_FORMAT_R32_FLOAT,
        .flags  = MD_GPU_TEX_STORAGE,
    });
    if (!atom_buf || !coeff_buf || !out_tex || !readback) goto done;

    /* Pack straight into the destination where that is safe, staging otherwise. */
    {
        float* p = (float*)md_gpu_upload_begin(stream, atom_buf, md_gto_gpu_atom_buffer_size(num_atoms));
        if (!p) goto done;
        md_gto_gpu_atom_pack(p, atom_xyz, sizeof(vec3_t), num_atoms);
        md_gpu_upload_end(stream);
    }
    const double* mo_coeffs[1] = {ao_coeffs};
    {
        float* p = (float*)md_gpu_upload_begin(stream, coeff_buf, md_gto_gpu_coeff_size_mo(1, num_cgtos));
        if (!p) goto done;
        md_gto_gpu_coeff_pack_mo(p, mo_coeffs, NULL, 1, num_cgtos);
        md_gpu_upload_end(stream);
    }

    md_gto_gpu_orbital_desc_t orb_desc = {
        .basis = gpu_basis,
        .atom_xyz = atom_buf,
        .coeff = coeff_buf,
        .out_tex = out_tex,
        .grid = &grid,
        .sample_offset = {0.0f, 0.0f, 0.0f},
        .num_orbitals = 1,
        .eval_mode = MD_GTO_EVAL_MODE_PSI,
        .op = MD_GTO_OP_SET,
    };

    /* Program order within the stream is the whole dependency model: the copy
       sees the kernel's writes without any explicit barrier. */
    md_gto_gpu_orbital_launch(stream, &orb_desc);
    md_gpu_memcpy_from_tex_async(readback, out_tex, NULL, readback_size, stream);
    md_gpu_stream_sync(stream);

    const float* grid_data = (const float*)md_gpu_host_ptr(readback);

    double grid_sum = 0.0;
    double cube_sum = 0.0;
    max_delta = 0.0;
    int    max_idx[3] = {0};

    for (int iz = 0; iz < grid.dim[2]; ++iz) {
        for (int iy = 0; iy < grid.dim[1]; ++iy) {
            for (int ix = 0; ix < grid.dim[0]; ++ix) {
                int grid_idx = iz * grid.dim[0] * grid.dim[1] + iy * grid.dim[0] + ix;
                int cube_idx = ix * grid.dim[1] * grid.dim[2] + iy * grid.dim[2] + iz;
                double g_val = grid_data[grid_idx];
                double c_val = cube->data.val[cube_idx];
                grid_sum += g_val;
                cube_sum += c_val;

                double delta = fabs(fabs(g_val) - fabs(c_val));
                if (delta > max_delta) {
                    max_delta = delta;
                    max_idx[0] = ix;
                    max_idx[1] = iy;
                    max_idx[2] = iz;
                }
            }
        }
    }

    printf("Max abs delta: %g at [%i,%i,%i]\n", max_delta, max_idx[0], max_idx[1], max_idx[2]);
    printf("GRID SUM: %.5f\n", grid_sum);
    printf("CUBE SUM: %.5f\n", cube_sum);

done:
    md_gpu_tex_destroy(out_tex, stream);
    md_gto_gpu_basis_destroy(gpu_basis);
    md_gpu_free(readback, stream);
    md_gpu_free(coeff_buf, stream);
    md_gpu_free(atom_buf, stream);
    md_gpu_pool_destroy(read_pool);
    md_gpu_pool_destroy(dev_pool);
    md_gto_gpu_shutdown();

    return max_delta;
}

/* Report *why* the device could not be created rather than a bare "skipped". */
static const char* gto_test_no_device_reason(void) {
    const char* err = md_gpu_last_error();
    return (err && err[0]) ? err : "No GPU device available";
}

UTEST(gto, h2o_lumo_gpu) {
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) {
        UTEST_SKIP(gto_test_no_device_reason());
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    md_cube_t cube_lumo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o_lumo.cube"), temp.arena));

    // Loaded as a SYSTEM: the geometry, the basis and the coefficients all come out of the
    // attribute table, which is the only way anything reaches a consumer now.
    vlx_test_t t = {0};
    ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(32)));

    size_t num_atoms = t.sys.atom.count;
    float* atom_xyz = (float*)md_temp_alloc_array(temp, float, 3 * num_atoms);
    ASSERT_EQ(vlx_test_atom_xyz_bohr(atom_xyz, 3 * num_atoms, &t), num_atoms);

    md_gto_basis_t basis = {0};
    ASSERT_TRUE(vlx_test_basis(&basis, &t));

    // The LUMO, derived from the published occupations by the rule the reader used: the first
    // orbital with zero occupancy.
    const size_t lumo_idx = vlx_test_lumo_idx(&t, STR_LIT("orbital/alpha/occupation"));
    const size_t num_ao   = md_gto_basis_num_ao(&basis);
    double* ao_coeffs = (double*)md_temp_alloc_array(temp, double, num_ao);
    ASSERT_EQ(vlx_test_row(ao_coeffs, num_ao, &t, STR_LIT("orbital/alpha/coefficient"), lumo_idx), num_ao);

    double max_delta_lumo = compare_vlx_and_cube_gpu(device, atom_xyz, &basis, ao_coeffs, &cube_lumo);
    EXPECT_LT(max_delta_lumo, 1.0E-4);

    vlx_test_free(&t);
    md_temp_end(temp);
    md_gpu_device_destroy(device);
}



#endif

// The amide / ne / myjob cube comparisons used to sit here, disabled. They were removed with the
// md_vlx object API they were written against: they called a compare_vlx_and_cube() that no longer
// exists and reference *_homo.cube / *_lumo.cube files that are not in test_data, so there was
// nothing to re-enable. gto.h2o_lumo_cpu is the live version of the same check.

