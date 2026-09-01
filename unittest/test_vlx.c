#include "utest.h"

#include <md_vlx.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

#include <md_vlx.c>

static const double ref_ener_tot = -444.518500783179;

UTEST(vlx, parse) {
	md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());
    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"));
    ASSERT_TRUE(result);

    EXPECT_EQ(0,  md_vlx_molecular_charge(vlx));
    EXPECT_EQ(26, md_vlx_number_of_atoms(vlx));
    EXPECT_EQ(1,  md_vlx_spin_multiplicity(vlx));
    EXPECT_EQ(41, md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_ALPHA));
    EXPECT_EQ(41, md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_BETA));

	const dvec3_t* coords = md_vlx_atom_coordinates(vlx);
    ASSERT_TRUE(coords != NULL);

    EXPECT_NEAR(-3.259400000000, coords[0].x, 1.0e-5);
    EXPECT_NEAR( 0.145200000000, coords[0].y, 1.0e-5);
    EXPECT_NEAR(-0.048400000000, coords[0].z, 1.0e-5);

    EXPECT_TRUE(str_eq(md_vlx_basis_set_ident(vlx), STR_LIT("DEF2-SVP")));

    int num_iter = md_vlx_scf_history_size(vlx);
	const double* energy = md_vlx_scf_history_energy(vlx);

    ASSERT_TRUE(num_iter > 0);
    EXPECT_NEAR(ref_ener_tot, energy[num_iter - 1], 1.0e-5);

    // @TODO: Test RSP

    md_vlx_destroy(vlx);
}

UTEST(vlx, minimal_example) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5");

    md_vlx_t* vlx = md_vlx_create(arena);
    if (!md_vlx_parse_file(vlx, path)) {
        MD_LOG_ERROR("Could not parse VLX file '"STR_FMT"'", STR_ARG(path));
        return;
    }

    const dvec3_t* vlx_coords = md_vlx_atom_coordinates(vlx);
    vec3_t* atom_xyz = (vec3_t*)md_arena_allocator_push(arena, sizeof(vec3_t) * md_vlx_number_of_atoms(vlx));

    for (size_t i = 0; i < md_vlx_number_of_atoms(vlx); i++) {
        atom_xyz[i] = vec3_set(
            (float)(vlx_coords[i].x * ANGSTROM_TO_BOHR),
            (float)(vlx_coords[i].y * ANGSTROM_TO_BOHR),
            (float)(vlx_coords[i].z * ANGSTROM_TO_BOHR)
        );
    }

    // The volume dimensions which we aim to sample molecular orbital over
    const int vol_dim = 80;

    // The molecular orbital index we aim to sample
    size_t mo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);

    // Extract GTOs
    md_gto_basis_t basis = {0};
    md_vlx_gto_basis_extract(&basis, vlx, arena);

    size_t num_gtos = md_gto_pgto_count(&basis);
    md_gto_t* gtos = (md_gto_t*)md_arena_allocator_push(arena, sizeof(md_gto_t) * num_gtos);
    
    size_t num_aos = md_vlx_scf_number_of_atomic_orbitals(vlx);
    double* mo_coeffs = (double*)md_arena_allocator_push(arena, sizeof(double) * num_aos);
    md_vlx_scf_mo_coefficients_extract(mo_coeffs, vlx, mo_idx, MD_VLX_SPIN_ALPHA);

    md_gto_expand_with_ao_coeffs(gtos, &basis, (const float*)atom_xyz, sizeof(vec3_t), mo_coeffs, 1.0e-6);

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
    min_aabb = vec3_sub1(min_aabb, pad);
    max_aabb = vec3_add1(max_aabb, pad);

    printf("min_box: %g %g %g \n", min_aabb.x, min_aabb.y, min_aabb.z);
    printf("max_box: %g %g %g \n", max_aabb.x, max_aabb.y, max_aabb.z);

    vec3_t ext_aabb = vec3_sub(max_aabb, min_aabb);
    vec3_t step_size = vec3_div1(ext_aabb, (float)vol_dim);

    // Shift origin by half a voxel such that the samples are constructed from the center of each voxel
    vec3_t origin = vec3_add(min_aabb, vec3_mul1(step_size, 0.5f));

    // Allocate data for storing the result
    float* vol_data = md_arena_allocator_push(arena, sizeof(float) * vol_dim * vol_dim * vol_dim);
    MEMSET(vol_data, 0, sizeof(float) * vol_dim * vol_dim * vol_dim);

    // Setup the grid structure that control how we aim to sample over space
    md_grid_t grid = (md_grid_t) {
        .orientation = mat3_ident(),
        .origin = origin,
        .spacing = step_size,
        .dim = {vol_dim, vol_dim, vol_dim},
    };

    // Evaluate the GTOs over the supplied grid
    md_gto_grid_evaluate(vol_data, &grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);

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

UTEST(vlx, mol) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"));
    EXPECT_TRUE(result);

    md_vlx_destroy(vlx);
}

UTEST(vlx, scf_results) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/tq.scf.results.h5"));
    EXPECT_TRUE(result);

    md_vlx_destroy(vlx);
}

UTEST(vlx, acro_rsp) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-rsp.h5"));
    EXPECT_TRUE(result);

    md_vlx_destroy(vlx);
}

#if 0
UTEST(vlx, acro_xps) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    // Acrolein (C C C O H H H H) with an '/xps' group alongside the existing '/rsp' group.
    // The two carbons at index 0 and 1 share a delocalized core hole at the same energy.
    bool result = md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-xps.h5"));
    EXPECT_TRUE(result);

    EXPECT_TRUE(md_vlx_has_xps(vlx));
    EXPECT_EQ(4, md_vlx_xps_count(vlx));
    EXPECT_EQ(2, md_vlx_xps_group_count(vlx));

    // Entries in the file are unordered; the parser sorts by (element, ionization_energy).
    const md_vlx_xps_entry_t* e = md_vlx_xps_entries(vlx);
    ASSERT_TRUE(e != NULL);
    for (size_t i = 1; i < md_vlx_xps_count(vlx); ++i) {
        EXPECT_TRUE(e[i - 1].element < e[i].element ||
                   (e[i - 1].element == e[i].element && e[i - 1].ionization_energy <= e[i].ionization_energy));
    }

    // Carbon: 3 entries, contiguous, lowest energy first, the delocalized pair leading.
    const md_vlx_xps_group_t* c = md_vlx_xps_group_by_element(vlx, 6);
    ASSERT_TRUE(c != NULL);
    EXPECT_EQ(3, c->count);
    EXPECT_EQ(e, c->entries);   // first group is a view onto the start of the flat array
    EXPECT_TRUE(c->entries[0].is_delocalized);
    EXPECT_TRUE(c->entries[1].is_delocalized);
    EXPECT_FALSE(c->entries[2].is_delocalized);
    EXPECT_NEAR(291.02, c->entries[0].ionization_energy, 1.0e-9);
    EXPECT_NEAR(293.51, c->entries[2].ionization_energy, 1.0e-9);
    EXPECT_NEAR(0.9644, c->entries[2].contribution, 1.0e-9);
    EXPECT_EQ(2, c->entries[2].atom_index);
    EXPECT_EQ(6, c->entries[2].mo_index);

    // Oxygen: single entry, and its view must start where carbon's ends.
    const md_vlx_xps_group_t* o = md_vlx_xps_group_by_element(vlx, 8);
    ASSERT_TRUE(o != NULL);
    EXPECT_EQ(1, o->count);
    EXPECT_EQ(c->entries + c->count, o->entries);
    EXPECT_EQ(3, o->entries[0].atom_index);
    EXPECT_NEAR(538.74, o->entries[0].ionization_energy, 1.0e-9);

    // Element not present in the calculation.
    EXPECT_TRUE(md_vlx_xps_group_by_element(vlx, 1) == NULL);
    EXPECT_TRUE(md_vlx_xps_group_by_index(vlx, 2) == NULL);

    md_vlx_destroy(vlx);
}
#endif

UTEST(vlx, acro_rsp_has_no_xps) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

    EXPECT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-rsp.h5")));
    EXPECT_FALSE(md_vlx_has_xps(vlx));
    EXPECT_EQ(0, md_vlx_xps_count(vlx));
    EXPECT_EQ(0, md_vlx_xps_group_count(vlx));
    EXPECT_TRUE(md_vlx_xps_entries(vlx) == NULL);

    md_vlx_destroy(vlx);
}

// 'orbital/{alpha,beta}/density' is never stored: it is reconstructed on demand from the
// coefficient and occupation attributes (vlx_build_occupation_density_matrix). This is the check
// that the reconstruction is actually the same matrix the file carried, not just plausible.
UTEST(vlx, orbital_density_matches_stored) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5")));

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_vlx_system_init_from_data(&sys, &state, vlx));
    md_vlx_publish_attributes(&sys, vlx);

    const size_t num_ao = md_vlx_scf_number_of_atomic_orbitals(vlx);
    ASSERT_TRUE(num_ao > 0);

    double* reconstructed = (double*)md_alloc(alloc, sizeof(double) * num_ao * num_ao);

    const md_attribute_t* alpha_density = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/density"));
    ASSERT_TRUE(alpha_density != NULL);
    EXPECT_EQ(alpha_density->storage, MD_ATTRIBUTE_STORAGE_VIRTUAL);
    EXPECT_EQ(alpha_density->format.shape[0], (uint32_t)num_ao);
    EXPECT_EQ(alpha_density->format.shape[1], (uint32_t)num_ao);
    ASSERT_EQ(md_attribute_extract_f64(reconstructed, num_ao * num_ao, alpha_density, md_unit_none()), num_ao * num_ao);

    const double* stored_alpha = md_vlx_scf_density_matrix_data(vlx, MD_VLX_SPIN_ALPHA);
    ASSERT_TRUE(stored_alpha != NULL);
    ASSERT_EQ(md_vlx_scf_density_matrix_size(vlx), num_ao);

    double max_diff = 0.0;
    for (size_t i = 0; i < num_ao * num_ao; ++i) {
        max_diff = MAX(max_diff, fabs(reconstructed[i] - stored_alpha[i]));
    }
    EXPECT_NEAR(0.0, max_diff, 1.0e-8);

    // h2o.h5 is a restricted calculation, so beta is alpha: it is published as an ALIAS rather
    // than as a second reconstruction, which is what makes reading either name give the same
    // matrix for no extra work. Absent, it would be a hole in the table for every consumer that
    // asks for a spin by name.
    const md_attribute_t* beta_density = md_attributes_find(&sys.attributes, STR_LIT("orbital/beta/density"));
    ASSERT_TRUE(beta_density != NULL);
    EXPECT_EQ(beta_density->storage, MD_ATTRIBUTE_STORAGE_ALIAS);
    EXPECT_TRUE(md_attribute_same_data(beta_density, alpha_density));

    // And it reads through alpha's provider. An alias of a VIRTUAL attribute used to return nothing
    // at all here, silently, which is what made 'orbital/total/density' produce zero values.
    double* via_beta = (double*)md_alloc(alloc, sizeof(double) * num_ao * num_ao);
    ASSERT_EQ(md_attribute_extract_f64(via_beta, num_ao * num_ao, beta_density, md_unit_none()), num_ao * num_ao);
    for (size_t i = 0; i < num_ao * num_ao; ++i) {
        ASSERT_NEAR(reconstructed[i], via_beta[i], 1.0e-12);
    }

    md_vlx_destroy(vlx);
    md_arena_allocator_destroy(alloc);
}

// The transition densities: attachment, detachment and their difference. There is nothing to
// compare them against - the file carries no reference matrices - so what is checked here is the
// MECHANICS and the invariants that hold whatever the numbers are.
//
// They are the one place in the table where a slice is load bearing rather than a convenience. The
// attribute is {S,A,A} and every plane costs a full reconstruction from the response eigenvectors,
// so asking for one state has to mean reconstructing one state. Both paths are exercised below.
static void check_transition_density_attributes(int* utest_result, str_t file) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());
    ASSERT_TRUE(md_vlx_parse_file(vlx, file));

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_vlx_system_init_from_data(&sys, &state, vlx));
    md_vlx_publish_attributes(&sys, vlx);

    const size_t num_states = md_vlx_rsp_number_of_excited_states(vlx);
    ASSERT_TRUE(num_states > 0);

    // None of these files names num_core/num_valence/num_virtual, so the split the solution vectors
    // are indexed by is derived from the SCF occupations and checked against the vector's own
    // length. Without it every reconstruction below declines - which is exactly what it used to do.
    ASSERT_TRUE(vlx->rsp.num_virtual > 0);
    const size_t nocc = vlx->rsp.num_core > 0 ? vlx->rsp.num_core : vlx->rsp.num_valence;
    const size_t nvir = vlx->rsp.num_virtual;
    ASSERT_TRUE(nocc > 0);

    const size_t vec_len = vlx->rsp.solution_matrix.size[1];
    EXPECT_TRUE(vec_len == nocc * nvir || vec_len == 2 * nocc * nvir);

    str_t paths[3] = {
        STR_LIT("vlx/rsp/transition_density/attachment"),
        STR_LIT("vlx/rsp/transition_density/detachment"),
        STR_LIT("vlx/rsp/transition_density/difference"),
    };

    const md_attribute_t* attr[3] = {0};
    for (int i = 0; i < 3; ++i) {
        attr[i] = md_attributes_find(&sys.attributes, paths[i]);
        ASSERT_TRUE(attr[i] != NULL);

        // {S,A,A}: state outermost, then the AO x AO matrix. Square is load bearing downstream -
        // the GL and GPU density paths pack only the upper triangle.
        EXPECT_EQ(attr[i]->storage, MD_ATTRIBUTE_STORAGE_VIRTUAL);
        EXPECT_EQ(attr[i]->format.type, MD_ATTRIBUTE_TYPE_F64);
        ASSERT_EQ(attr[i]->format.rank, 3u);
        EXPECT_EQ(attr[i]->format.shape[0], (uint32_t)num_states);
        EXPECT_EQ(attr[i]->format.shape[1], attr[i]->format.shape[2]);

        // A virtual attribute hands out no resident storage to write into.
        EXPECT_TRUE(md_attributes_data(&sys.attributes, attr[i]->id, MD_ATTRIBUTE_TYPE_F64) == NULL);
    }

    const size_t num_ao = attr[0]->format.shape[1];
    ASSERT_TRUE(num_ao > 0);
    EXPECT_EQ(num_ao, md_vlx_scf_number_of_atomic_orbitals(vlx));

    const size_t plane = num_ao * num_ao;
    double* mat[3] = {0};

    // One state at a time, which is the shape a representation asks in.
    for (int i = 0; i < 3; ++i) {
        const md_attribute_slice_t slice = md_attribute_slice_1(0);

        md_attribute_format_t sliced = {0};
        ASSERT_TRUE(md_attribute_slice_format(&sliced, attr[i], &slice));
        EXPECT_EQ(sliced.rank, 2u);
        EXPECT_EQ(sliced.shape[0], (uint32_t)num_ao);
        EXPECT_EQ(sliced.shape[1], (uint32_t)num_ao);
        ASSERT_EQ(md_attribute_slice_count(attr[i], &slice), plane);

        mat[i] = (double*)md_alloc(alloc, sizeof(double) * plane);
        ASSERT_EQ(md_attribute_extract_slice_f64(mat[i], plane, attr[i], &slice, md_unit_none()), plane);
    }

    double max_asym = 0.0;
    double max_value = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t r = 0; r < num_ao; ++r) {
            for (size_t c = 0; c < num_ao; ++c) {
                const double v = mat[i][r * num_ao + c];
                ASSERT_TRUE(v == v);                    // not NaN
                ASSERT_TRUE(fabs(v) < DBL_MAX);         // not inf
                max_asym  = MAX(max_asym,  fabs(v - mat[i][c * num_ao + r]));
                max_value = MAX(max_value, fabs(v));
            }
        }
    }

    // Symmetric by construction (vlx_symmetrize_square), and the GL/GPU density paths read only the
    // upper triangle - so an asymmetric matrix would be silently half consumed.
    EXPECT_NEAR(0.0, max_asym, 1.0e-12);

    // Reconstructed from a real excitation, so not the zero matrix a declining provider would have
    // been indistinguishable from before the extract started reporting a short write.
    EXPECT_TRUE(max_value > 1.0e-8);

    // The difference IS attachment minus detachment. Cheap to state and it pins the one relation
    // between the three that no reference values are needed to know.
    double max_diff = 0.0;
    for (size_t i = 0; i < plane; ++i) {
        max_diff = MAX(max_diff, fabs(mat[2][i] - (mat[0][i] - mat[1][i])));
    }
    EXPECT_NEAR(0.0, max_diff, 1.0e-12);

    // The whole attribute, every state at once. Same values in state 0's plane as the slice gave,
    // which is the property that makes a slice an optimisation rather than a second code path.
    for (int i = 0; i < 3; ++i) {
        const size_t total = num_states * plane;
        ASSERT_EQ(md_attribute_element_count(&attr[i]->format), total);

        double* all = (double*)md_alloc(alloc, sizeof(double) * total);
        ASSERT_EQ(md_attribute_extract_f64(all, total, attr[i], md_unit_none()), total);

        for (size_t v = 0; v < plane; ++v) {
            ASSERT_NEAR(mat[i][v], all[v], 1.0e-12);
        }
    }

    // A state past the end selects nothing rather than reading past the array.
    const md_attribute_slice_t past_end = md_attribute_slice_1((uint32_t)num_states);
    EXPECT_EQ(md_attribute_slice_count(attr[0], &past_end), 0u);
    EXPECT_EQ(md_attribute_extract_slice_f64(mat[0], plane, attr[0], &past_end, md_unit_none()), 0u);

    md_vlx_destroy(vlx);
    md_arena_allocator_destroy(alloc);
}

// Two files on purpose: h2o is one excited state in a small basis, amide is 110 atomic orbitals -
// the size the "provider wrote 0 of 12100" failure was reported at, and the one where a per state
// reconstruction is expensive enough that the slice has to mean what it says.
UTEST(vlx, transition_density_attributes_h2o) {
    check_transition_density_attributes(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"));
}

UTEST(vlx, transition_density_attributes_amide) {
    check_transition_density_attributes(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.h5"));
}

// The same reconstruction reached the other way, through the reader's own entry point. If these two
// ever disagree the attribute is not a view of the file's data but a second implementation of it.
static void check_transition_density_matches_direct(int* utest_result, str_t file) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());
    ASSERT_TRUE(md_vlx_parse_file(vlx, file));

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_vlx_system_init_from_data(&sys, &state, vlx));
    md_vlx_publish_attributes(&sys, vlx);

    const size_t dim = md_vlx_rsp_transition_density_matrix_size(vlx, 0);
    ASSERT_TRUE(dim > 0);
    const size_t plane = dim * dim;

    const md_vlx_transition_type_t types[3] = {
        MD_VLX_TRANSITION_ATTACHMENT,
        MD_VLX_TRANSITION_DETACHMENT,
        MD_VLX_TRANSITION_DIFFERENCE,
    };
    str_t paths[3] = {
        STR_LIT("vlx/rsp/transition_density/attachment"),
        STR_LIT("vlx/rsp/transition_density/detachment"),
        STR_LIT("vlx/rsp/transition_density/difference"),
    };

    double* direct = (double*)md_alloc(alloc, sizeof(double) * plane);
    double* viaattr = (double*)md_alloc(alloc, sizeof(double) * plane);

    for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(md_vlx_rsp_transition_density_matrix_extract(direct, vlx, 0, types[i]), dim);

        const md_attribute_t* a = md_attributes_find(&sys.attributes, paths[i]);
        ASSERT_TRUE(a != NULL);
        ASSERT_EQ(a->format.shape[1], (uint32_t)dim);

        const md_attribute_slice_t slice = md_attribute_slice_1(0);
        ASSERT_EQ(md_attribute_extract_slice_f64(viaattr, plane, a, &slice, md_unit_none()), plane);

        double max_diff = 0.0;
        for (size_t v = 0; v < plane; ++v) {
            max_diff = MAX(max_diff, fabs(direct[v] - viaattr[v]));
        }
        EXPECT_NEAR(0.0, max_diff, 1.0e-12);
    }

    md_vlx_destroy(vlx);
    md_arena_allocator_destroy(alloc);
}

UTEST(vlx, transition_density_attribute_matches_direct_extract) {
    check_transition_density_matches_direct(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"));
}

UTEST(vlx, transition_density_attribute_matches_direct_extract_amide) {
    check_transition_density_matches_direct(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.h5"));
}

// The combined spin densities, which are derivations over derivations: alpha and beta are each
// computed on demand, and these read both. That is legal because the graph stays acyclic, and it is
// the case that broke - beta is an ALIAS of alpha in a restricted calculation, and reading an alias
// of a computed attribute used to return nothing at all, silently.
//
// A restricted calculation gives the invariants for free and needs no reference values: the total
// is exactly twice alpha, and the difference is exactly zero.
UTEST(vlx, combined_spin_densities_h2o) {
    md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());
    ASSERT_TRUE(md_vlx_parse_file(vlx, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5")));
    ASSERT_EQ(md_vlx_scf_type(vlx), MD_VLX_SCF_RESTRICTED);

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_vlx_system_init_from_data(&sys, &state, vlx));
    md_vlx_publish_attributes(&sys, vlx);

    const size_t num_ao = md_vlx_scf_number_of_atomic_orbitals(vlx);
    ASSERT_TRUE(num_ao > 0);
    const size_t plane = num_ao * num_ao;

    str_t paths[4] = {
        STR_LIT("orbital/alpha/density"),
        STR_LIT("orbital/beta/density"),
        STR_LIT("orbital/total/density"),
        STR_LIT("orbital/difference/density"),
    };

    double* mat[4] = {0};
    for (int i = 0; i < 4; ++i) {
        const md_attribute_t* a = md_attributes_find(&sys.attributes, paths[i]);
        ASSERT_TRUE(a != NULL);
        ASSERT_EQ(a->format.rank, 2u);
        EXPECT_EQ(a->format.shape[0], (uint32_t)num_ao);
        EXPECT_EQ(a->format.shape[1], (uint32_t)num_ao);

        mat[i] = (double*)md_alloc(alloc, sizeof(double) * plane);
        ASSERT_EQ(md_attribute_extract_f64(mat[i], plane, a, md_unit_none()), plane);
    }

    double max_total = 0.0, max_diff = 0.0, max_beta = 0.0, magnitude = 0.0;
    for (size_t i = 0; i < plane; ++i) {
        max_beta  = MAX(max_beta,  fabs(mat[1][i] - mat[0][i]));            // beta IS alpha
        max_total = MAX(max_total, fabs(mat[2][i] - 2.0 * mat[0][i]));      // total = alpha + beta
        max_diff  = MAX(max_diff,  fabs(mat[3][i]));                        // difference = 0
        magnitude = MAX(magnitude, fabs(mat[0][i]));
    }
    EXPECT_NEAR(0.0, max_beta,  1.0e-12);
    EXPECT_NEAR(0.0, max_total, 1.0e-12);
    EXPECT_NEAR(0.0, max_diff,  1.0e-12);

    // Guards the invariants above against the case where every matrix is zero, which would satisfy
    // all three and is exactly what a declining provider used to produce.
    EXPECT_TRUE(magnitude > 1.0e-8);

    md_vlx_destroy(vlx);
    md_arena_allocator_destroy(alloc);
}
