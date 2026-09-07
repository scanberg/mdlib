#include "ubench.h"

#include <md_vlx.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>

#include <md_gto.h>
//#include <md_gto.c>

#include <float.h>

#define BLK_DIM 8
#define ANGSTROM_TO_BOHR 1.0 / 0.529177210903

UBENCH_EX(gto, evaluate_grid) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    
    // md_vlx has no reader object: a file is loaded into a system and everything it carried - the
    // geometry, the basis, the MO coefficients - is read back out of that system's attribute table.
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    if (!md_vlx_system_init_from_file(&sys, &state, STR_LIT(MD_BENCHMARK_DATA_DIR "/vlx/mol.out"))) {
        MD_LOG_ERROR("Could not load benchmark vlx file");
        return;
    }

    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);

    const md_attribute_t* coord_attr = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/coordinate"));
    const size_t num_atoms = coord_attr ? md_attribute_value_count(&coord_attr->format) : 0;
    double* coords = (double*)md_arena_allocator_push(arena, sizeof(double) * 3 * num_atoms);
    md_attribute_extract_f64(coords, 3 * num_atoms, coord_attr, md_unit_none());

    vec3_t* atom_xyz = (vec3_t*)md_arena_allocator_push(arena, sizeof(vec3_t) * num_atoms);
    for (size_t i = 0; i < num_atoms; ++i) {
        atom_xyz[i] = vec3_set((float)(coords[3*i+0] * ANGSTROM_TO_BOHR), (float)(coords[3*i+1] * ANGSTROM_TO_BOHR), (float)(coords[3*i+2] * ANGSTROM_TO_BOHR));
        min_box = vec3_min(min_box, atom_xyz[i]);
        max_box = vec3_max(max_box, atom_xyz[i]);
    }

    min_box = vec3_sub1(min_box, 2.0f);
    max_box = vec3_add1(max_box, 2.0f);

    // Conversion from Ångström to Bohr
    const float factor = ANGSTROM_TO_BOHR;

    min_box = vec3_mul1(min_box, factor);
    max_box = vec3_mul1(max_box, factor);

    int vol_dim = 256;
    size_t bytes = sizeof(float) * vol_dim * vol_dim * vol_dim;
    float* vol_data = md_arena_allocator_push(arena, bytes);
    MEMSET(vol_data, 0, bytes);

    vec3_t step = vec3_div1(vec3_sub(max_box, min_box), (float)vol_dim);

    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin = min_box,
        .spacing = step,
        .dim = {vol_dim, vol_dim, vol_dim},
    };

    md_gto_basis_t basis = {0};
    md_gto_basis_extract_attributes(&basis, &sys.attributes, arena);

    size_t num_gtos = md_gto_pgto_count(&basis);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);

    size_t num_aos = md_gto_basis_num_ao(&basis);
    double* mo_coeffs = (double*)md_arena_allocator_push(arena, sizeof(double) * num_aos);
    const md_attribute_t* coeff_attr = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_slice_t mo_slice = md_attribute_slice_1(120);
    md_attribute_extract_slice_f64(mo_coeffs, num_aos, coeff_attr, &mo_slice, md_unit_none());

    md_gto_expand_with_ao_coeffs(gtos, &basis, (const float*)atom_xyz, sizeof(vec3_t), mo_coeffs, 1.0e-6);

    md_gto_t* sub_gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);

    int num_blk = vol_dim / BLK_DIM;
    int num_tot_blk = num_blk * num_blk * num_blk;

    UBENCH_DO_BENCHMARK() {
        for (int i = 0; i < num_tot_blk; ++i) {
            // Determine block (x/y/z) index from linear index (i)
            int blk_z = i & (num_blk - 1);
            int blk_y = (i / num_blk) & (num_blk - 1);
            int blk_x = i / (num_blk * num_blk);

            const int off_idx[3] = {blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM};
            const int len_idx[3] = {BLK_DIM, BLK_DIM, BLK_DIM};

            float aabb_min[3] = {
                grid.origin.x + off_idx[0] * grid.spacing.x,
                grid.origin.y + off_idx[1] * grid.spacing.y,
                grid.origin.z + off_idx[2] * grid.spacing.z,
            };
            float aabb_max[3] = {
                grid.origin.x + (off_idx[0] + len_idx[0]) * grid.spacing.x,
                grid.origin.y + (off_idx[1] + len_idx[1]) * grid.spacing.y,
                grid.origin.z + (off_idx[2] + len_idx[2]) * grid.spacing.z,
            };

            size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);
            
            md_gto_grid_evaluate_sub(vol_data, &grid, off_idx, len_idx, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
            //evaluate_grid_ortho_8x8x8_256(grid.data, off_idx, len_idx, grid.origin, step.elem, sub_pgtos, num_sub_pgtos, MD_GTO_EVAL_MODE_PSI);
            //evaluate_grid_8x8x8_256(grid.data, off_idx, len_idx, grid.origin, grid.step_x, grid.step_y, grid.step_z, sub_gtos, num_sub_gtos, MD_GTO_EVAL_MODE_PSI);
        }
    }

    UBENCH_DO_NOTHING(&grid);

    md_arena_allocator_destroy(arena);
}
