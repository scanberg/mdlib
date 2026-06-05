#include "utest.h"

#include <core/md_allocator.h>
#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include <md_vlx.h>
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
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    md_vlx_t* vlx = md_vlx_create(temp_arena);
    bool read = md_vlx_parse_file(vlx, filename);

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);

	const dvec3_t* coords = md_vlx_atom_coordinates(vlx);
    for (size_t i = 0; i < md_vlx_number_of_atoms(vlx); ++i) {
        vec3_t c = {(float)coords[i].x, (float)coords[i].y, (float)coords[i].z};
        min_box = vec3_min(min_box, c);
        max_box = vec3_max(max_box, c);
    }

    min_box = vec3_sub1(min_box, 2.0f);
    max_box = vec3_add1(max_box, 2.0f);

    min_box = vec3_mul1(min_box, ANGSTROM_TO_BOHR);
    max_box = vec3_mul1(max_box, ANGSTROM_TO_BOHR);

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
    md_vlx_gto_basis_extract(&basis, vlx, temp_arena);

    size_t num_atoms = md_vlx_number_of_atoms(vlx);
    float* atom_xyz = md_temp_alloc_array(temp, float, 3 * num_atoms);
    for (size_t i = 0; i < num_atoms; i++) {
        atom_xyz[3*i+0] = (float)(coords[i].x * ANGSTROM_TO_BOHR);
        atom_xyz[3*i+1] = (float)(coords[i].y * ANGSTROM_TO_BOHR);
        atom_xyz[3*i+2] = (float)(coords[i].z * ANGSTROM_TO_BOHR);
    }

    size_t num_ao = md_vlx_scf_number_of_atomic_orbitals(vlx);
    double* mo_coeffs = md_temp_alloc_array(temp, double, num_ao);
    md_vlx_scf_mo_coefficients_extract(mo_coeffs, vlx, 120, MD_VLX_SPIN_ALPHA);

    *num_gtos = md_gto_pgto_count(&basis);
    *gtos = md_alloc(arena, sizeof(md_gto_t) * (*num_gtos));
    *num_gtos = md_gto_expand_with_ao_coeffs(*gtos, &basis, atom_xyz, sizeof(vec3_t), mo_coeffs, 1.0e-6);

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

    md_vlx_t* vlx = md_vlx_create(temp.arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5")));

    size_t num_atoms = md_vlx_number_of_atoms(vlx);
    const dvec3_t* vlx_coords = md_vlx_atom_coordinates(vlx);
    float* atom_xyz = (float*)md_temp_alloc_array(temp, float, 3 * num_atoms);
    for (size_t i = 0; i < num_atoms; i++) {
        atom_xyz[3*i+0] = (float)(vlx_coords[i].x * ANGSTROM_TO_BOHR);
        atom_xyz[3*i+1] = (float)(vlx_coords[i].y * ANGSTROM_TO_BOHR);
        atom_xyz[3*i+2] = (float)(vlx_coords[i].z * ANGSTROM_TO_BOHR);
    }

    md_gto_basis_t basis = {0};
    md_vlx_gto_basis_extract(&basis, vlx, temp_arena);
    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_ALPHA);
    const double* ao_coeffs = md_vlx_scf_mo_coefficients(vlx, lumo_idx, MD_VLX_SPIN_ALPHA);

    double max_delta_lumo = compare_vlx_and_cube_cpu(atom_xyz, &basis, ao_coeffs, &cube_lumo);
    EXPECT_LT(max_delta_lumo, 1.0E-4);  

    md_temp_end(temp);
}

#if MD_ENABLE_GPU

static double compare_vlx_and_cube_gpu(md_gpu_device_t device, const float* atom_xyz, const md_gto_basis_t* gto_basis, const double* ao_coeffs, const md_cube_t* cube) {
    const double cutoff_value = 1.0E-6;

    double max_delta = DBL_MAX;

    md_gto_gpu_basis_t gpu_basis = NULL;
    md_gpu_buffer_t atom_buf = NULL;
    md_gpu_buffer_t coeff_buf = NULL;
    md_gpu_image_t out_image = NULL;
    md_gpu_buffer_t readback_buf = NULL;

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

    gpu_basis = md_gto_gpu_basis_create(device, &(md_gto_gpu_basis_desc_t){
        .basis = gto_basis,
        .cutoff = 0.0,
    });
    if (!gpu_basis) goto done;

    const uint32_t num_atoms = md_gto_gpu_basis_num_atoms(gpu_basis);
    const uint32_t num_cgtos = md_gto_gpu_basis_num_cgtos(gpu_basis);
    const size_t voxel_count = (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];
    const size_t readback_size = sizeof(float) * voxel_count;

    atom_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = md_gto_gpu_atom_buffer_size(num_atoms),
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    coeff_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = md_gto_gpu_coeff_size_mo(1, num_cgtos),
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    out_image = md_gpu_image_create(device, &(md_gpu_image_desc_t){
        .width  = (uint32_t)grid.dim[0],
        .height = (uint32_t)grid.dim[1],
        .depth  = (uint32_t)grid.dim[2],
        .format = MD_GPU_IMAGE_FORMAT_R32_FLOAT,
        .flags  = MD_GPU_IMAGE_STORAGE,
    });
    readback_buf = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){
        .size  = readback_size,
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    if (!atom_buf || !coeff_buf || !out_image || !readback_buf) goto done;

    md_gto_gpu_atom_pack((float*)md_gpu_buffer_cpu_ptr(atom_buf), atom_xyz, sizeof(vec3_t), num_atoms);
    const double* mo_coeffs[1] = {ao_coeffs};
    md_gto_gpu_coeff_pack_mo((float*)md_gpu_buffer_cpu_ptr(coeff_buf), mo_coeffs, NULL, 1, num_cgtos);

    md_gpu_queue_t queue = md_gpu_queue_compute(device);
    if (!queue) goto done;

    md_gpu_cmd_t cmd = md_gpu_cmd_begin(queue, "gto cube gpu eval");
    if (!cmd) goto done;

    md_gto_gpu_orbital_desc_t orb_desc = {
        .basis = gpu_basis,
        .atom_xyz = atom_buf,
        .coeff = coeff_buf,
        .out_image = out_image,
        .grid = &grid,
        .sample_offset = {0.0f, 0.0f, 0.0f},
        .num_orbitals = 1,
        .eval_mode = MD_GTO_EVAL_MODE_PSI,
        .op = MD_GTO_OP_SET,
    };

    md_gto_gpu_orbital_record(cmd, &orb_desc);
    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_TRANSFER);
    md_gpu_cmd_copy_image_region_to_buffer(cmd, out_image, (md_gpu_image_region_t){0}, readback_buf, 0);
    md_gpu_cmd_end(cmd);

    md_gpu_event_t event = md_gpu_queue_submit_one(queue, cmd);
    md_gpu_event_wait(event);

    const float* grid_data = (const float*)md_gpu_buffer_cpu_ptr(readback_buf);

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
    md_gpu_buffer_destroy(readback_buf);
    md_gpu_image_destroy(out_image);
    md_gpu_buffer_destroy(coeff_buf);
    md_gpu_buffer_destroy(atom_buf);
    md_gto_gpu_basis_destroy(gpu_basis);
    md_gto_gpu_shutdown();

    return max_delta;
}

UTEST(gto, h2o_lumo_gpu) {
    md_gpu_device_t device = md_gpu_device_create();
    if (!device) {
        UTEST_SKIP("No GPU device available");
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    md_cube_t cube_lumo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o_lumo.cube"), temp.arena));

    md_vlx_t* vlx = md_vlx_create(temp.arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5")));

    size_t num_atoms = md_vlx_number_of_atoms(vlx);
    const dvec3_t* vlx_coords = md_vlx_atom_coordinates(vlx);
    float* atom_xyz = (float*)md_temp_alloc_array(temp, float, 3 * num_atoms);
    for (size_t i = 0; i < num_atoms; i++) {
        atom_xyz[3*i+0] = (float)(vlx_coords[i].x * ANGSTROM_TO_BOHR);
        atom_xyz[3*i+1] = (float)(vlx_coords[i].y * ANGSTROM_TO_BOHR);
        atom_xyz[3*i+2] = (float)(vlx_coords[i].z * ANGSTROM_TO_BOHR);
    }

    md_gto_basis_t basis = {0};
    md_vlx_gto_basis_extract(&basis, vlx, temp_arena);
    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_ALPHA);
    const double* ao_coeffs = md_vlx_scf_mo_coefficients(vlx, lumo_idx, MD_VLX_SPIN_ALPHA);

    double max_delta_lumo = compare_vlx_and_cube_gpu(device, atom_xyz, &basis, ao_coeffs, &cube_lumo);
    EXPECT_LT(max_delta_lumo, 1.0E-4);

    md_temp_end(temp);
    md_gpu_device_destroy(device);
}



#endif

#if 0
UTEST(gto, amide) {
    md_temp_scope_t temp = md_temp_begin();

    md_cube_t cube_lumo = {0};
    md_cube_t cube_homo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide_lumo.cube"), temp.arena));
    ASSERT_TRUE(md_cube_file_load(&cube_homo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide_homo.cube"), temp.arena));

    md_vlx_t* vlx = md_vlx_create(temp.arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_ALPHA);
    size_t homo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, VALUE_CUTOFF, &cube_lumo, temp.arena);
    double max_delta_homo = compare_vlx_and_cube(vlx, homo_idx, VALUE_CUTOFF, &cube_homo, temp.arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);
    EXPECT_LT(max_delta_homo, 1.0E-4);

    md_temp_end(temp);
}

UTEST(gto, ne) {
    md_temp_scope_t temp = md_temp_begin();

    md_cube_t cube_lumo = {0};
    md_cube_t cube_homo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/ne_lumo.cube"), temp.arena));
    ASSERT_TRUE(md_cube_file_load(&cube_homo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/ne_homo.cube"), temp.arena));

    md_vlx_t* vlx = md_vlx_create(temp.arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/ne.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_ALPHA);
    size_t homo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, VALUE_CUTOFF, &cube_lumo, temp.arena);
    double max_delta_homo = compare_vlx_and_cube(vlx, homo_idx, VALUE_CUTOFF, &cube_homo, temp.arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);
    EXPECT_LT(max_delta_homo, 1.0E-4);

    md_temp_end(temp);
}

UTEST(gto, myjob) {
    md_temp_scope_t temp = md_temp_begin();

    md_cube_t cube_lumo = {0};
    md_cube_t cube_homo = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_lumo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/myjob_lumo.cube"), temp.arena));
    ASSERT_TRUE(md_cube_file_load(&cube_homo, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/myjob_homo.cube"), temp.arena));

    md_vlx_t* vlx = md_vlx_create(temp.arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/myjob.out")));

    size_t lumo_idx = md_vlx_scf_lumo_idx(vlx, MD_VLX_SPIN_ALPHA);
    size_t homo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);

    double max_delta_lumo = compare_vlx_and_cube(vlx, lumo_idx, VALUE_CUTOFF, &cube_lumo, temp.arena);
    double max_delta_homo = compare_vlx_and_cube(vlx, homo_idx, VALUE_CUTOFF, &cube_homo, temp.arena);

    EXPECT_LT(max_delta_lumo, 1.0E-4);
    EXPECT_LT(max_delta_homo, 1.0E-4);

    md_temp_end(temp);
}
#endif
