#include "utest.h"

#include <md_vlx.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

#include <md_vlx.c>

static const int ref_iter[] = {
    1,2,3,4,5,6,7,8,9,10,11,12
};
static const double ref_ener_tot[] = {
    -444.502770501781,
    -444.502826985951,
    -444.518269328072,
    -444.518408360201,
    -444.518497179843,
    -444.518500090662,
    -444.518500772981,
    -444.518500782306,
    -444.518500783079,
    -444.518500783141,
    -444.518500783156,
    -444.518500783158,
};
static const double ref_ener_change[] = {
    0.0000000000,
    -0.0000564842,
    -0.0154423421,
    -0.0001390321,
    -0.0000888196,
    -0.0000029108,
    -0.0000006823,
    -0.0000000093,
    -0.0000000008,
    -0.0000000001,
    -0.0000000000,
    -0.0000000000,
};
static const double ref_grad_norm[] = {
    0.33661898,
    0.33937059,
    0.04132164,
    0.02546891,
    0.00508750,
    0.00230342,
    0.00028359,
    0.00008491,
    0.00002482,
    0.00000991,
    0.00000406,
    0.00000094,
};
static const double ref_max_grad[] = {
    0.01168705,
    0.01156775,
    0.00163688,
    0.00092076,
    0.00023894,
    0.00009228,
    0.00000913,
    0.00000327,
    0.00000086,
    0.00000045,
    0.00000018,
    0.00000003,
};
static const double ref_density_change[] = {
    0.00000000,
    0.31704185,
    0.17166350,
    0.02943676,
    0.01204202,
    0.00304038,
    0.00107142,
    0.00015964,
    0.00003900,
    0.00001305,
    0.00000550,
    0.00000235,
};

UTEST(vlx, vlx_parse) {
	md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());
    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out"));
    ASSERT_TRUE(result);

#if 1
    EXPECT_EQ(0,  md_vlx_molecular_charge(vlx));
    EXPECT_EQ(26, md_vlx_number_of_atoms(vlx));
    EXPECT_EQ(1,  md_vlx_spin_multiplicity(vlx));
    EXPECT_EQ(41, md_vlx_number_of_alpha_electrons(vlx));
    EXPECT_EQ(41, md_vlx_number_of_beta_electrons(vlx));

	const dvec3_t* coords = md_vlx_atom_coordinates(vlx);
    ASSERT_TRUE(coords != NULL);

    EXPECT_NEAR(-3.259400000000, coords[0].x, 1.0e-5);
    EXPECT_NEAR( 0.145200000000, coords[0].y, 1.0e-5);
    EXPECT_NEAR(-0.048400000000, coords[0].z, 1.0e-5);

    EXPECT_TRUE(str_eq(md_vlx_basis_set_ident(vlx), STR_LIT("DEF2-SVP")));

    int num_iter = md_vlx_scf_history_size(vlx);
	const double* energy = md_vlx_scf_history_energy(vlx);
    const double* energy_diff = md_vlx_scf_history_energy_diff(vlx);
	const double* grad_norm = md_vlx_scf_history_gradient_norm(vlx);
	const double* max_grad = md_vlx_scf_history_max_gradient(vlx);
	const double* density_diff = md_vlx_scf_history_density_diff(vlx);

    ASSERT_EQ(12, num_iter);
    for (size_t i = 0; i < num_iter; ++i) {
        EXPECT_NEAR(ref_ener_tot[i], energy[i], 1.0e-5);
        EXPECT_NEAR(ref_ener_change[i], energy_diff[i], 1.0e-5);
        EXPECT_NEAR(ref_grad_norm[i], grad_norm[i], 1.0e-5);
        EXPECT_NEAR(ref_max_grad[i], max_grad[i], 1.0e-5);
        EXPECT_NEAR(ref_density_change[i], density_diff[i], 1.0e-5);
    }

    // @TODO: Test RSP

    md_vlx_destroy(vlx);
#endif
}

UTEST(vlx, correctness) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_vlx_t* vlx = md_vlx_create(arena);
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.out")));

    const int vol_dim = 80;

    size_t bytes = sizeof(float) * vol_dim * vol_dim * vol_dim;
    float* vol_data = md_arena_allocator_push(arena, bytes);
    MEMSET(vol_data, 0, bytes);

    float* ref_data = md_arena_allocator_push(arena, bytes);
    MEMSET(ref_data, 0, bytes);

    md_grid_t grid = (md_grid_t) {
        .data = vol_data,
        .dim = {vol_dim, vol_dim, vol_dim},
        .origin = {-5.744767, -5.744767, -5.522177},
        .step_x = {0.143619, 0, 0},
        .step_y = {0, 0.143619, 0},
        .step_z = {0, 0, 0.143619},
    };

    size_t mo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_MO_TYPE_ALPHA);

    size_t num_gtos = md_vlx_mo_gto_count(vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mo_gto_extract(gtos, vlx, mo_idx, MD_VLX_MO_TYPE_ALPHA);

    size_t cap_phi = num_gtos;
    double* phi = (double*)md_arena_allocator_push(arena, sizeof(double) * cap_phi);

    size_t num_mo_coeffs = number_of_mo_coefficients(&vlx->scf.alpha);
    double* mo_coeffs = md_arena_allocator_push(arena, sizeof(double) * num_mo_coeffs);
    extract_mo_coefficients(mo_coeffs, &vlx->scf.alpha, mo_idx);

    int beg_idx[3] = {0, 0, 0};
    int end_idx[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

    md_gto_grid_evaluate(&grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    for (int iz = 0; iz < grid.dim[2]; ++iz) {
        double z = grid.origin[2] + grid.step_z[2] * iz;
        z *= BOHR_TO_ANGSTROM;
        for (int iy = 0; iy < grid.dim[1]; ++iy) {
            double y = grid.origin[1] + grid.step_y[1] * iy;
            y *= BOHR_TO_ANGSTROM;
            for (int ix = 0; ix < grid.dim[0]; ++ix) {
                double x = grid.origin[0] + grid.step_x[0] * ix;
                x *= BOHR_TO_ANGSTROM;

                size_t num_phi = compPhiAtomicOrbitals(phi, cap_phi, vlx->atom_coordinates, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set, x, y, z);

                ASSERT(num_phi == num_mo_coeffs);
                double psi = 0.0;
                for (size_t i = 0; i < num_mo_coeffs; ++i) {
                    psi += mo_coeffs[i] * phi[i]; 
                }

                int idx = iz * grid.dim[1] * grid.dim[0] + iy * grid.dim[0] + ix;
                float grid_psi = grid.data[idx];
                ref_data[idx] = (float)psi;
            }
        }
    }

    md_arena_allocator_destroy(arena);
}

UTEST(vlx, minimal_example) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.out");

    md_vlx_t* vlx = md_vlx_create(arena);
    if (!md_vlx_parse_file(vlx, path)) {
        MD_LOG_ERROR("Could not parse VLX file '"STR_FMT"'", STR_ARG(path));
        return;
    }

    // The volume dimensions which we aim to sample molecular orbital over
    const int vol_dim = 80;

    // The molecular orbital index we aim to sample
    size_t mo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_MO_TYPE_ALPHA);

    // Extract GTOs
    size_t num_gtos = md_vlx_mo_gto_count(vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mo_gto_extract(gtos, vlx, mo_idx, MD_VLX_MO_TYPE_ALPHA);

    // Compute radial cutoffs for the GTOs
    md_gto_cutoff_compute(gtos, num_gtos, 1.0e-6);

    // Calculate bounding box (AABB)
    vec3_t min_aabb = vec3_set1( FLT_MAX);
    vec3_t max_aabb = vec3_set1(-FLT_MAX);

    for (size_t i = 0; i < num_gtos; ++i) {
        vec3_t coord = vec3_set(gtos[i].x, gtos[i].y, gtos[i].z);
        min_aabb = vec3_min(min_aabb, coord);
        max_aabb = vec3_max(max_aabb, coord);
    }

    // Add some padding
    const float pad = 6.0f;
    min_aabb = vec3_sub_f(min_aabb, pad);
    max_aabb = vec3_add_f(max_aabb, pad);

    printf("min_box: %g %g %g \n", min_aabb.x, min_aabb.y, min_aabb.z);
    printf("max_box: %g %g %g \n", max_aabb.x, max_aabb.y, max_aabb.z);

    vec3_t ext_aabb = vec3_sub(max_aabb, min_aabb);
    vec3_t step_size = vec3_div_f(ext_aabb, (float)vol_dim);

    // Shift origin by half a voxel such that the samples are constructed from the center of each voxel
    vec3_t origin = vec3_add(min_aabb, vec3_mul_f(step_size, 0.5f));

    // Allocate data for storing the result
    size_t bytes = sizeof(float) * vol_dim * vol_dim * vol_dim;
    float* vol_data = md_arena_allocator_push(arena, bytes);
    MEMSET(vol_data, 0, bytes);

    // Setup the grid structure that control how we aim to sample over space
    md_grid_t grid = (md_grid_t) {
        .data = vol_data,
        .dim = {vol_dim, vol_dim, vol_dim},
        .origin = {origin.x, origin.y, origin.z},
        .step_x = {step_size.x, 0, 0},
        .step_y = {0, step_size.y, 0},
        .step_z = {0, 0, step_size.z},
    };

    // Evaluate the GTOs over the supplied grid
    md_gto_grid_evaluate(&grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    // Define a subrange over the center of the volume
    int range_min = vol_dim / 2 - vol_dim / 8;
    int range_max = vol_dim / 2 + vol_dim / 8;

    // Print out the values over a subrange of the volume
    for (int iz = range_min; iz < range_max; ++iz) {
        for (int iy = range_min; iy < range_max; ++iy) {
            for (int ix = range_min; ix < range_max; ++ix) {
                int idx = iz * grid.dim[1] * grid.dim[0] + iy * grid.dim[0] + ix;
                //printf("%g ", grid.data[idx]);
            }
        }
    }

    md_arena_allocator_destroy(arena);
}

#if 0
UTEST(vlx, scf_results_h2o) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_vlx_scf_results_t scf = {0};
    bool result = md_vlx_read_scf_results(&scf, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o-rcs.scf.results.h5"), arena);
    EXPECT_TRUE(result);

    ASSERT_EQ(3, scf.atom.number_of_atoms);
    EXPECT_EQ(8, scf.atom.nuclear_charges[0]);
    EXPECT_EQ(1, scf.atom.nuclear_charges[1]);
    EXPECT_EQ(1, scf.atom.nuclear_charges[2]);

    EXPECT_EQ(0.0, scf.atom.coordinates[0].x);
    EXPECT_EQ(0.0, scf.atom.coordinates[0].y);
    EXPECT_EQ(0.0, scf.atom.coordinates[0].z);

    EXPECT_EQ(0.0, scf.atom.coordinates[1].x);
    EXPECT_EQ(0.0, scf.atom.coordinates[1].y);
    EXPECT_EQ(3.3925119279731675, scf.atom.coordinates[1].z);

    EXPECT_EQ(0.0, scf.molecular_charge);

    bool empty = str_empty(scf.potfile_text);

    md_arena_allocator_destroy(arena);
}

UTEST(vlx, scf_results_meth) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_vlx_scf_results_t scf = {0};
    bool result = md_vlx_read_scf_results(&scf, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/meth.scf.results.h5"), arena);
    EXPECT_TRUE(result);

    md_arena_allocator_destroy(arena);
}
#endif

UTEST(vlx, mol_out) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out"));
    EXPECT_TRUE(result);

    md_vlx_destroy(vlx);
}

UTEST(vlx, scf_results) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/tq.scf.results.h5"));
    EXPECT_TRUE(result);

    md_vlx_destroy(vlx);
}

UTEST(vlx, h5_pure) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-rsp.h5"));
    EXPECT_TRUE(result);

    md_vlx_destroy(vlx);
}