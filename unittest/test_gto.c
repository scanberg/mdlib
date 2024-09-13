#include "utest.h"

#include <core/md_arena_allocator.h>
#include <md_vlx.h>
#include <md_cube.h>

#include <md_gto.h>

// Conversion from Ångström to Bohr
const float factor = 1.0 / 0.529177210903;
#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.529177210903

static void init(md_grid_t* grid, md_gto_t** gtos, size_t* num_gtos, int vol_dim, str_t filename, md_allocator_i* arena) {
    md_vlx_data_t vlx = {0};
    md_vlx_data_parse_file(&vlx, filename, arena);

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);

    for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
        vec3_t c = {(float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]};
        min_box = vec3_min(min_box, c);
        max_box = vec3_max(max_box, c);
    }

    min_box = vec3_sub_f(min_box, 2.0f);
    max_box = vec3_add_f(max_box, 2.0f);

    min_box = vec3_mul_f(min_box, factor);
    max_box = vec3_mul_f(max_box, factor);

    size_t bytes = sizeof(float) * vol_dim * vol_dim * vol_dim;
    float* vol_data = md_arena_allocator_push(arena, bytes);
    MEMSET(vol_data, 0, bytes);

    vec3_t step = vec3_div_f(vec3_sub(max_box, min_box), (float)vol_dim);

    *grid = (md_grid_t) {
        .data = vol_data,
        .dim = {vol_dim, vol_dim, vol_dim},
        .origin = {min_box.x, min_box.y, min_box.z},
        .step_x = {step.x, 0, 0},
        .step_y = {0, step.y, 0},
        .step_z = {0, 0, step.z},
    };

    *num_gtos = md_vlx_mol_gto_count(&vlx);
    *gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * *num_gtos);
    md_vlx_mol_gto_extract(*gtos, &vlx, 120);
    md_gto_cutoff_compute(*gtos, *num_gtos, 1.0e-6);
}

UTEST(gto, evaluate_grid) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    
    md_grid_t grid;
    md_gto_t* gtos;
    size_t num_gtos;

    int vol_dim = 64;
    init(&grid, &gtos, &num_gtos, vol_dim, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out"), arena);

    md_grid_t ref_grid = grid;
    ref_grid.data = md_alloc(arena, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    MEMSET(grid.data, 0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);
    MEMSET(ref_grid.data, 0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    md_gto_t* sub_gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);

    // Number of NxNxN blocks in each dimension
    int num_blk[3] = {
        grid.dim[0] / BLK_DIM,
        grid.dim[1] / BLK_DIM,
        grid.dim[2] / BLK_DIM,
    };

    int num_tot_blk = num_blk[0] * num_blk[1] * num_blk[2];

    float step_ortho[3] = {
        grid.step_x[0],
        grid.step_y[1],
        grid.step_z[2],
    };

    for (int i = 0; i < num_tot_blk; ++i) {
        // Determine block index from linear input index i
        int blk_x = i % num_blk[0];
        int blk_y = (i / num_blk[0]) % num_blk[1];
        int blk_z = i / (num_blk[0] * num_blk[1]);

        int off_idx[3] = { blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM };
        int len_idx[3] = { BLK_DIM, BLK_DIM, BLK_DIM };

        int beg_idx[3] = {off_idx[0], off_idx[1], off_idx[2]};
        int end_idx[3] = {off_idx[0] + len_idx[0], off_idx[1] + len_idx[1], off_idx[2] + len_idx[2]};

        float aabb_min[3] = {
            grid.origin[0] + beg_idx[0] * grid.step_x[0],
            grid.origin[1] + beg_idx[1] * grid.step_y[1],
            grid.origin[2] + beg_idx[2] * grid.step_z[2],
        };
        float aabb_max[3] = {
            grid.origin[0] + end_idx[0] * grid.step_x[0],
            grid.origin[1] + end_idx[1] * grid.step_y[1],
            grid.origin[2] + end_idx[2] * grid.step_z[2],
        };

        size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);

        md_gto_grid_evaluate_sub(&grid, beg_idx, end_idx, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

        //evaluate_grid_ref(ref_grid.data, beg_idx, end_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);

        const float eps = 1.0e-5f;

#if 0
#ifdef __SSE2__
        evaluate_grid_ortho_8x8x8_128(grid.data, beg_idx, grid.dim, grid.origin, step_ortho, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }

        evaluate_grid_8x8x8_128(grid.data, beg_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }
#endif

#ifdef __AVX2__
        evaluate_grid_ortho_8x8x8_256(grid.data, beg_idx, grid.dim, grid.origin, step_ortho, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }

        evaluate_grid_8x8x8_256(grid.data, beg_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        for (int z = beg_idx[2]; z < end_idx[2]; ++z) {
            for (int y = beg_idx[1]; y < end_idx[1]; ++y) {
                for (int x = beg_idx[0]; x < end_idx[0]; ++x) {
                    int idx = z * grid.dim[0] * grid.dim[1] + y * grid.dim[0] + x;
                    ASSERT_NEAR(ref_grid.data[idx], grid.data[idx], eps);
                }
            }
        }
#endif
#endif
    }

    md_arena_allocator_destroy(arena);
}

UTEST(gto, amide) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    int vol_dim = 80;
    md_vlx_data_t vlx = {0};
    md_vlx_data_parse_file(&vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.out"), arena);

    md_grid_t grid = {
        .data = md_arena_allocator_push(arena, sizeof(float) * vol_dim * vol_dim * vol_dim),
        .dim = {vol_dim, vol_dim, vol_dim},
        .origin = {-7.952157f, -8.874910f, -7.201935f},
        .step_x = {0.218627f, 0, 0},
        .step_y = {0, 0.215167f, 0},
        .step_z = {0, 0, 0.178142f},
    };

    size_t mo_idx = vlx.scf.lumo_idx;

    size_t num_gtos = md_vlx_nto_gto_count(&vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mol_gto_extract(gtos, &vlx, mo_idx);
    md_gto_cutoff_compute(gtos, num_gtos, 0);

    MEMSET(grid.data,     0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    md_gto_grid_evaluate(&grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    double sum = 0.0;

    size_t count = grid.dim[0] * grid.dim[1] * grid.dim[2];
    for (size_t i = 0; i < count; ++i) {
        sum += grid.data[i] * grid.data[i];
    }

    EXPECT_NEAR(119.308, sum, 1.0e-3);

    md_arena_allocator_destroy(arena);
}

UTEST(gto, h2o) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube = {0};
    ASSERT_TRUE(md_cube_file_load(&cube, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o_lumo.cube"), arena));

    md_vlx_data_t vlx = {0};
    ASSERT_TRUE(md_vlx_data_parse_file(&vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.out"), arena));

    const float scl = 1.0f;

    size_t mo_idx = vlx.scf.lumo_idx;

    /*
    size_t num_rows = vlx.scf.alpha.orbitals.dim[0];
    size_t num_cols = vlx.scf.alpha.orbitals.dim[1];

    for (size_t i = 0; i < num_rows; ++i) {
        printf("mo[%zu]: %g\n", i, vlx.scf.alpha.orbitals.data[i * num_cols + mo_idx]);
    }
    */

    md_grid_t grid = {
        .data = md_arena_allocator_push(arena, sizeof(float) * cube.data.num_x * cube.data.num_y * cube.data.num_z),
        .dim = {cube.data.num_x, cube.data.num_y, cube.data.num_z},
        .origin = {cube.origin[0] * scl, cube.origin[1] * scl, cube.origin[2] * scl},
        .step_x = {cube.xaxis[0] * scl, cube.xaxis[1] * scl, cube.xaxis[2] * scl},
        .step_y = {cube.yaxis[0] * scl, cube.yaxis[1] * scl, cube.yaxis[2] * scl},
        .step_z = {cube.zaxis[0] * scl, cube.zaxis[1] * scl, cube.zaxis[2] * scl},
    };

    size_t num_gtos = md_vlx_mol_gto_count(&vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mol_gto_extract(gtos, &vlx, mo_idx);
    md_gto_cutoff_compute(gtos, num_gtos, 0);

    size_t count = grid.dim[0] * grid.dim[1] * grid.dim[2];

    float* psi = md_arena_allocator_push(arena, sizeof(float) * count);
    vec3_t* xyz = md_arena_allocator_push(arena, sizeof(vec3_t) * count);
    for (int z = 0; z < grid.dim[2]; ++z) {
        for (int y = 0; y < grid.dim[1]; ++y) {
            for (int x = 0; x < grid.dim[0]; ++x) {
                int i = z * grid.dim[1] * grid.dim[0] + y * grid.dim[0] + x;
                xyz[i] = vec3_set(
                    grid.origin[0] + grid.step_x[0] * x,
                    grid.origin[1] + grid.step_y[1] * y,
                    grid.origin[2] + grid.step_z[2] * z);
            }
        }
    }

    md_gto_xyz_evaluate(psi, (float*)xyz, count, sizeof(vec3_t), gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    MEMSET(grid.data, 0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    int beg_idx[3] = {0, 0, 0};
    int end_idx[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

    md_gto_grid_evaluate(&grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    double xyz_sum  = 0.0;
    double grid_sum = 0.0;
    double cube_sum = 0.0;

    double max_delta = 0.0;
    int    max_idx = 0;

    for (int iz = 0; iz < grid.dim[2]; ++iz) {
        for (int iy = 0; iy < grid.dim[1]; ++iy) {
            for (int ix = 0; ix < grid.dim[0]; ++ix) {
                int grid_idx = iz * grid.dim[0] * grid.dim[1] + iy * grid.dim[0] + ix;
                int cube_idx = ix * grid.dim[1] * grid.dim[2] + iy * grid.dim[2] + iz;
                double g_val = grid.data[grid_idx] * grid.data[grid_idx];
                double c_val = cube.data.val[cube_idx] * cube.data.val[cube_idx];
                double p_val = psi[grid_idx] * psi[grid_idx];
                grid_sum += g_val;
                cube_sum += c_val;
                xyz_sum  += p_val;

                double delta = fabs(g_val - c_val);
                if (delta > max_delta) {
                    max_delta = delta;
                    max_idx = grid_idx;
                }
            }
        }
    }


    int ix =  max_idx %  grid.dim[0];
    int iy = (max_idx /  grid.dim[0]) % grid.dim[1];
    int iz =  max_idx / (grid.dim[0]  * grid.dim[1]);

    printf("Max i: %i\n", max_idx);
    printf("Max delta: %g at [%i,%i,%i]\n", max_delta, ix, iy, iz);

    EXPECT_LT(max_delta, 4.0e-5);

    printf("H2O GRID SUM: %.3f\n", grid_sum);
    printf("H2O CUBE SUM: %.3f\n", cube_sum);
    printf("H2O XYZ  SUM: %.3f\n", xyz_sum);

    ASSERT_EQ(grid.dim[0], cube.data.num_x);
    ASSERT_EQ(grid.dim[1], cube.data.num_y);
    ASSERT_EQ(grid.dim[2], cube.data.num_z);
    ASSERT_EQ(1, cube.data.num_m);

    EXPECT_NEAR(grid_sum, cube_sum, 1.0e-2);
    EXPECT_NEAR(xyz_sum,  cube_sum, 1.0e-2);

    md_arena_allocator_destroy(arena);
}

// Single contracted basis function
typedef struct basis_set_func_t {
    uint8_t  type; // Azimuthal Quantum Number
    uint8_t  param_count;
    uint16_t param_offset;
} basis_set_func_t;

typedef struct basis_set_basis_t {
    uint8_t  max_type;
    uint8_t  basis_func_count;
    uint16_t basis_func_offset;
} basis_set_basis_t;

typedef struct basis_set_t {
    str_t label;
    struct {
        size_t count;
        double* exponents;
        double* normalization_coefficients;
    } param;

    struct {
        size_t count;
        basis_set_func_t* data;
    } basis_func;

    // The atom basis entries are implicitly stored in the order of atomic numbers
    // 0 would be a NULL entry, 1 corresponds to Hydrogen, 2 corresponds to Helium etc.
    struct {
        size_t count;
        basis_set_basis_t* data;
    } atom_basis;
} basis_set_t;

typedef struct basis_func_range_t {
    int beg;
    int end;
} basis_func_range_t;

static basis_set_basis_t* basis_set_get_atom_basis(const basis_set_t* basis_set, int atomic_number) {
    if (atomic_number < basis_set->atom_basis.count) {
        return basis_set->atom_basis.data + atomic_number;
    }
    return NULL;
}

static basis_func_range_t basis_get_atomic_angl_basis_func_range(const basis_set_t* basis_set, int atomic_number, int angl) {
    basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
    ASSERT(atom_basis);

    basis_func_range_t range = {0};

    int beg = atom_basis->basis_func_offset;
    int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
    for (int i = beg; i < end; ++i) {
        int type = basis_set->basis_func.data[i].type;
        if (type == angl) {
            range.beg = (range.end == 0) ? i : range.beg;
            range.end = i + 1;
        }
    }

    return range;
}

typedef struct basis_func_t {
    int type;
    int count;
    double* exponents;
    double* normalization_coefficients;
} basis_func_t;

static basis_func_t get_basis_func(const basis_set_t* basis_set, int basis_func_idx) {
    basis_set_func_t func = basis_set->basis_func.data[basis_func_idx];
    return (basis_func_t) {
        .type = func.type,
            .count = func.param_count,
            .exponents = basis_set->param.exponents + func.param_offset,
            .normalization_coefficients = basis_set->param.normalization_coefficients + func.param_offset,
    };
}

#if 0
UTEST(gto, co_lumo) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube = {0};
    ASSERT_TRUE(md_cube_file_load(&cube, STR_LIT("e:/data/veloxchem/system/co/lumo.cube"), arena));

    md_vlx_data_t vlx = {0};
    ASSERT_TRUE(md_vlx_data_parse_file(&vlx, STR_LIT("e:/data/veloxchem/system/co/co.out"), arena));

    size_t mo_idx = vlx.scf.lumo_idx;

    size_t num_rows = vlx.scf.alpha.orbitals.dim[0];
    size_t num_cols = vlx.scf.alpha.orbitals.dim[1];

    double* mo = md_arena_allocator_push(arena, sizeof(double) * num_rows);

    for (size_t i = 0; i < num_rows; ++i) {
        printf("mo[%zu]: %g\n", i, vlx.scf.alpha.orbitals.data[i * num_cols + mo_idx]);
    }

    int elem[2] = {6, 8};
    for (int i = 0; i < ARRAY_SIZE(elem); ++i) {
        printf ("elem_id: %i\n", elem[i]);
        for (int angl = 0; angl < 5; ++angl) {
            basis_func_range_t range = basis_get_atomic_angl_basis_func_range(vlx.basis.basis_set, elem[i], angl);
            if (range.end != range.beg) {
                printf ("  angl: %i\n", angl);
                for (int funcidx = range.beg; funcidx < range.end; funcidx++) {
                    basis_func_t basis_func = get_basis_func(vlx.basis.basis_set, funcidx);

                    // process primitives
                    int nprims = basis_func.count;
                    const double* exp = basis_func.exponents;
                    const double* fac = basis_func.normalization_coefficients;
                    for (int j = 0; j < nprims; ++j) {
                        printf("    exp: %g, fac: %g\n", exp[j], fac[j]);
                    }
                }
            }
        }
    }

    const float scl = 1.0f;

    md_grid_t grid = {
        .data = md_arena_allocator_push(arena, sizeof(float) * cube.data.num_x * cube.data.num_y * cube.data.num_z),
        .dim = {cube.data.num_x, cube.data.num_y, cube.data.num_z},
        .origin = {cube.origin[0] * scl, cube.origin[1] * scl, cube.origin[2] * scl},
        .step_x = {cube.xaxis[0] * scl, cube.xaxis[1] * scl, cube.xaxis[2] * scl},
        .step_y = {cube.yaxis[0] * scl, cube.yaxis[1] * scl, cube.yaxis[2] * scl},
        .step_z = {cube.zaxis[0] * scl, cube.zaxis[1] * scl, cube.zaxis[2] * scl},
    };

    size_t num_gtos = md_vlx_mol_gto_count(&vlx);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_mol_gto_extract(gtos, &vlx, mo_idx);

    MEMSET(grid.data,     0, sizeof(float) * grid.dim[0] * grid.dim[1] * grid.dim[2]);

    int beg_idx[3] = {0, 0, 0};
    int end_idx[3] = {grid.dim[0], grid.dim[1], grid.dim[2]};

    evaluate_grid_ref(grid.data, beg_idx, end_idx, grid.dim, grid.origin, grid.step_x, grid.step_y, grid.step_z, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

    double grid_sum = 0.0;
    double cube_sum = 0.0;

    size_t count = grid.dim[0] * grid.dim[1] * grid.dim[2];
    for (size_t i = 0; i < count; ++i) {
        grid_sum += grid.data[i] * grid.data[i];
        cube_sum += cube.data.val[i] * cube.data.val[i];
    }

    printf("NE GRID SUM: %.3f\n", grid_sum);
    printf("NE CUBE SUM: %.3f\n", cube_sum);

    ASSERT_EQ(grid.dim[0], cube.data.num_x);
    ASSERT_EQ(grid.dim[1], cube.data.num_y);
    ASSERT_EQ(grid.dim[2], cube.data.num_z);
    ASSERT_EQ(1, cube.data.num_m);

    for (int x = 0; x < grid.dim[0]; x++) {
        for (int y = 0; y < grid.dim[1]; y++) {
            for (int z = 0; z < grid.dim[2]; z++) {
                int cube_idx = x * cube.data.num_y * cube.data.num_z + y * cube.data.num_z + z;
                int grid_idx = z * grid.dim[1] * grid.dim[0] + y * grid.dim[0] + x;
                double cube_val = cube.data.val[cube_idx];
                double grid_val = grid.data[grid_idx];

                double delta = fabs(cube_val - grid_val);
                if (delta > 1.0e-5) {
                    while(0) {};            
                }
            }
        }
    }

    md_arena_allocator_destroy(arena);
}

UTEST(gto, amide_nto) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_cube_t cube_h = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_h,   STR_LIT("e:/data/veloxchem/system/tq-mod/tq_S2_NTO_H1.cube"), arena));

    md_cube_t cube_p = {0};
    ASSERT_TRUE(md_cube_file_load(&cube_p,   STR_LIT("e:/data/veloxchem/system/tq-mod/tq_S2_NTO_P1.cube"), arena));

    md_vlx_data_t vlx = {0};
    ASSERT_TRUE(md_vlx_data_parse_file(&vlx, STR_LIT("e:/data/veloxchem/system/tq-mod/tq.out"), arena));

    size_t num_rows = vlx.scf.alpha.orbitals.dim[0];
    size_t num_cols = vlx.scf.alpha.orbitals.dim[1];

    md_grid_t grid_h = {
        .data = md_arena_allocator_push(arena, sizeof(float) * cube_h.data.num_x * cube_h.data.num_y * cube_h.data.num_z),
        .dim = {cube_h.data.num_x, cube_h.data.num_y, cube_h.data.num_z},
        .origin = {cube_h.origin[0], cube_h.origin[1], cube_h.origin[2]},
        .step_x = {cube_h.xaxis[0] , cube_h.xaxis[1] , cube_h.xaxis[2] },
        .step_y = {cube_h.yaxis[0] , cube_h.yaxis[1] , cube_h.yaxis[2] },
        .step_z = {cube_h.zaxis[0] , cube_h.zaxis[1] , cube_h.zaxis[2] },
    };

    md_grid_t grid_p = {
        .data = md_arena_allocator_push(arena, sizeof(float) * cube_p.data.num_x * cube_p.data.num_y * cube_p.data.num_z),
        .dim    = {cube_p.data.num_x, cube_p.data.num_y, cube_p.data.num_z},
        .origin = {cube_p.origin[0] , cube_p.origin[1] , cube_p.origin[2]},
        .step_x = {cube_p.xaxis[0]  , cube_p.xaxis[1]  , cube_p.xaxis[2] },
        .step_y = {cube_p.yaxis[0]  , cube_p.yaxis[1]  , cube_p.yaxis[2] },
        .step_z = {cube_p.zaxis[0]  , cube_p.zaxis[1]  , cube_p.zaxis[2] },
    };

    size_t num_gtos = md_vlx_mol_gto_count(&vlx);
    md_gto_t* gtos_h = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_nto_gto_extract(gtos_h, &vlx, 0, 0, MD_VLX_NTO_TYPE_HOLE);

    md_gto_t* gtos_p = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    md_vlx_nto_gto_extract(gtos_p, &vlx, 0, 0, MD_VLX_NTO_TYPE_PARTICLE);

    MEMSET(grid_h.data,     0, sizeof(float) * grid_h.dim[0] * grid_h.dim[1] * grid_h.dim[2]);
    MEMSET(grid_p.data,     0, sizeof(float) * grid_p.dim[0] * grid_p.dim[1] * grid_p.dim[2]);

    int beg_idx[3] = {0, 0, 0};

    evaluate_grid_ref(grid_h.data, beg_idx, grid_h.dim, grid_h.dim, grid_h.origin, grid_h.step_x, grid_h.step_y, grid_h.step_z, gtos_h, num_gtos, MD_GTO_EVAL_MODE_PSI);
    evaluate_grid_ref(grid_p.data, beg_idx, grid_p.dim, grid_p.dim, grid_p.origin, grid_p.step_x, grid_p.step_y, grid_p.step_z, gtos_p, num_gtos, MD_GTO_EVAL_MODE_PSI);

    size_t count = grid_h.dim[0] * grid_h.dim[1] * grid_h.dim[2];

    double grid_sum_h = 0.0;
    double grid_sum_p = 0.0;
    double cube_sum_h = 0.0;
    double cube_sum_p = 0.0;

    for (size_t i = 0; i < count; ++i) {
        grid_sum_p += grid_p.data[i] * grid_p.data[i];
        grid_sum_h += grid_h.data[i] * grid_h.data[i];
        cube_sum_p += cube_p.data.val[i] * cube_p.data.val[i];
        cube_sum_h += cube_h.data.val[i] * cube_h.data.val[i];
    }

    printf("GRID SUM H: %g\n", grid_sum_h);
    printf("CUBE SUM H: %g\n", cube_sum_h);
    printf("GRID SUM P: %g\n", grid_sum_p);
    printf("CUBE SUM P: %g\n", cube_sum_p);

    ASSERT_EQ(grid_h.dim[0], cube_h.data.num_x);
    ASSERT_EQ(grid_h.dim[1], cube_h.data.num_y);
    ASSERT_EQ(grid_h.dim[2], cube_h.data.num_z);
    ASSERT_EQ(1, cube_h.data.num_m);

    md_arena_allocator_destroy(arena);
}
#endif