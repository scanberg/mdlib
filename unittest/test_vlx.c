#include "utest.h"

#include "vlx_test_util.h"

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <float.h>

// Everything a VeloxChem file carries reaches a consumer as an md_system_t and its ATTRIBUTE TABLE:
// md_vlx.h has three entry points and no reader object to ask questions of. So these tests read what
// a real consumer reads, by the same paths, and a value that cannot be reached this way is a value
// no consumer can reach either - which is the property they exist to hold.

static const double ref_ener_tot = -444.518500783179;

UTEST(vlx, parse) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), MEGABYTES(64)));

	// The atoms land in the system itself, not in the table.
	EXPECT_EQ(26u, t.sys.atom.count);
	EXPECT_EQ(26u, t.state.num_atoms);

	EXPECT_EQ(0.0, vlx_test_scalar(&t, STR_LIT("vlx/molecular_charge"), -1.0));
	EXPECT_EQ(1.0, vlx_test_scalar(&t, STR_LIT("vlx/spin_multiplicity"), -1.0));
	EXPECT_EQ(41.0, vlx_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), -1.0));
	EXPECT_EQ(41.0, vlx_test_scalar(&t, STR_LIT("vlx/electron_count/beta"), -1.0));

	EXPECT_TRUE(str_eq(vlx_test_string(&t, STR_LIT("vlx/basis_set")), STR_LIT("DEF2-SVP")));

	// The QM geometry, in Angstrom, as the calculation was run at.
	const md_attribute_t* coord = vlx_test_attr(&t, STR_LIT("qm/atom/coordinate"));
	ASSERT_TRUE(coord != NULL);
	ASSERT_EQ(md_attribute_components(&coord->format), 3u);
	ASSERT_EQ(md_attribute_value_count(&coord->format), 26u);

	double xyz[3 * 26] = {0};
	ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_unit_none()), ARRAY_SIZE(xyz));
	EXPECT_NEAR(-3.259400000000, xyz[0], 1.0e-5);
	EXPECT_NEAR( 0.145200000000, xyz[1], 1.0e-5);
	EXPECT_NEAR(-0.048400000000, xyz[2], 1.0e-5);

	// The system's own state is that geometry too, narrowed to float.
	EXPECT_NEAR(xyz[0], (double)t.state.x[0], 1.0e-5);
	EXPECT_NEAR(xyz[1], (double)t.state.y[0], 1.0e-5);
	EXPECT_NEAR(xyz[2], (double)t.state.z[0], 1.0e-5);

	const size_t num_iter = vlx_test_count(&t, STR_LIT("vlx/scf/history/energy"));
	ASSERT_TRUE(num_iter > 0);

	double* energy = (double*)md_alloc(t.alloc, sizeof(double) * num_iter);
	ASSERT_EQ(vlx_test_series(energy, num_iter, &t, STR_LIT("vlx/scf/history/energy")), num_iter);
	EXPECT_NEAR(ref_ener_tot, energy[num_iter - 1], 1.0e-5);

	vlx_test_free(&t);
}

// The facts about a calculation that are TEXT and nothing else, plus the counts that used to be
// reachable only by asking the reader. If one of these stops being published, nothing downstream can
// tell a restricted run from an unrestricted one except by guessing from which arrays are shared.
UTEST(vlx, run_description_is_published) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(16)));

	EXPECT_TRUE(str_eq(vlx_test_string(&t, STR_LIT("vlx/scf/type")), STR_LIT("restricted")));
	EXPECT_TRUE(vlx_test_scalar(&t, STR_LIT("vlx/spin_multiplicity"), 0.0) == 1.0);
	EXPECT_TRUE(vlx_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), 0.0) > 0.0);
	EXPECT_TRUE(vlx_test_scalar(&t, STR_LIT("vlx/nuclear_repulsion_energy"), 0.0) > 0.0);

	// A response calculation, so the type is there and names which one.
	EXPECT_TRUE(str_eq(vlx_test_string(&t, STR_LIT("vlx/rsp/type")), STR_LIT("linear")));

	// No geometry optimisation in this file, so nothing under vlx/opt. An absent path is how a
	// consumer learns a block is missing - not a zero it would have to interpret.
	EXPECT_FALSE(vlx_test_has(&t, STR_LIT("vlx/opt/type")));
	EXPECT_FALSE(vlx_test_has(&t, STR_LIT("vlx/opt/state_index")));
	EXPECT_FALSE(vlx_test_has(&t, STR_LIT("vlx/opt/irc_ts_index")));

	vlx_test_free(&t);
}

// The Z/Y sign convention, pinned.
//
// mdlib DERIVES the NTOs from the response solution vector, as the SVD of T = Z - Y, where Z and
// Y name the two halves of that vector as VeloxChem writes them. VeloxChem also publishes its own
// NTO eigenvalues in 'rsp/nto_lambdas', computed independently by its own code. If the two agree,
// mdlib is combining the halves the way VeloxChem does; if a future change flips a sign, these
// numbers separate immediately and loudly.
//
// The reference values below were read straight out of the test files' nto_lambdas datasets.
// Verified 2026-09-02 that Z+Y and Z alone do NOT reproduce them (h2o: 1.0419 and 1.0006 for the
// leading value against the stored 0.9606), so this is a real discriminator and not a test that
// would pass under any convention.
UTEST(vlx, nto_lambdas_match_the_file) {
	struct {
		const char* path;
		size_t count;
		double lambda[4];
	} cases[] = {
		{ MD_UNITTEST_DATA_DIR "/vlx/h2o.h5",   4, { 0.960565601258, 0.000257984966, 0.000199122711, 0.000157261013 } },
		{ MD_UNITTEST_DATA_DIR "/vlx/amide.h5", 4, { 0.951067058354, 0.001038273356, 0.000565841191, 0.000215796789 } },
	};

	for (size_t c = 0; c < ARRAY_SIZE(cases); ++c) {
		vlx_test_t t = {0};
		ASSERT_TRUE(vlx_test_load(&t, str_from_cstr(cases[c].path), MEGABYTES(64)));

		// {S,Lmax}: one row per excited state, padded to the widest row with zeros.
		const md_attribute_t* a = vlx_test_attr(&t, STR_LIT("vlx/rsp/nto/lambda"));
		ASSERT_TRUE(a != NULL);
		ASSERT_EQ(a->format.rank, 2u);

		const size_t num_lambdas = a->format.shape[1];
		ASSERT_TRUE(num_lambdas >= cases[c].count);

		double* lambdas = (double*)md_alloc(t.alloc, sizeof(double) * num_lambdas);
		ASSERT_EQ(vlx_test_row(lambdas, num_lambdas, &t, STR_LIT("vlx/rsp/nto/lambda"), 0), num_lambdas);

		for (size_t i = 0; i < cases[c].count; ++i) {
			EXPECT_NEAR(cases[c].lambda[i], lambdas[i], 1.0e-6);
		}

		// Lambdas come out largest first, and they do not sum to one: the transition density is
		// not normalized the way the solution vector is. Anything downstream that treats them as
		// shares has to renormalize, which is what the charge transfer diagram does. The zero pad
		// on a short row sorts with them rather than against them, which is the whole reason zero
		// is the honest pad for a weight.
		for (size_t i = 1; i < num_lambdas; ++i) {
			EXPECT_TRUE(lambdas[i] <= lambdas[i - 1] + 1.0e-12);
		}

		vlx_test_free(&t);
	}
}

// What sampling a molecular orbital over a grid looks like end to end, with no reader in the loop:
// load the file into a system, rebuild the basis from the format neutral basis/ attributes, take one
// MO's coefficient row, evaluate. This doubles as the worked example for the whole interface.
UTEST(vlx, minimal_example) {
	vlx_test_t t = {0};
	str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5");
	ASSERT_TRUE(vlx_test_load(&t, path, MEGABYTES(32)));

	md_allocator_i* arena = t.alloc;

	const size_t num_atoms = t.sys.atom.count;
	float* atom_xyz = (float*)md_alloc(arena, sizeof(float) * 3 * num_atoms);
	ASSERT_EQ(vlx_test_atom_xyz_bohr(atom_xyz, 3 * num_atoms, &t), num_atoms);

	// The volume dimensions which we aim to sample molecular orbital over
	const int vol_dim = 80;

	// The molecular orbital index we aim to sample - the HOMO, one below the first empty orbital.
	const size_t lumo_idx = vlx_test_lumo_idx(&t, STR_LIT("orbital/alpha/occupation"));
	ASSERT_TRUE(lumo_idx > 0);
	const size_t mo_idx = lumo_idx - 1;

	// Extract the GTO basis from the attributes the loader published
	md_gto_basis_t basis = {0};
	ASSERT_TRUE(vlx_test_basis(&basis, &t));

	size_t num_gtos = md_gto_pgto_count(&basis);
	md_gto_t* gtos = (md_gto_t*)md_alloc(arena, sizeof(md_gto_t) * num_gtos);

	const size_t num_aos = md_gto_basis_num_ao(&basis);
	double* mo_coeffs = (double*)md_alloc(arena, sizeof(double) * num_aos);
	ASSERT_EQ(vlx_test_row(mo_coeffs, num_aos, &t, STR_LIT("orbital/alpha/coefficient"), mo_idx), num_aos);

	md_gto_expand_with_ao_coeffs(gtos, &basis, atom_xyz, sizeof(float) * 3, mo_coeffs, 1.0e-6);

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

	vec3_t ext_aabb = vec3_sub(max_aabb, min_aabb);
	vec3_t step_size = vec3_div1(ext_aabb, (float)vol_dim);

	// Shift origin by half a voxel such that the samples are constructed from the center of each voxel
	vec3_t origin = vec3_add(min_aabb, vec3_mul1(step_size, 0.5f));

	// Allocate data for storing the result
	float* vol_data = (float*)md_alloc(arena, sizeof(float) * vol_dim * vol_dim * vol_dim);
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

	// An occupied orbital is not the zero function; this is what catches a coefficient row that
	// came back empty rather than wrong.
	double max_value = 0.0;
	for (int i = 0; i < vol_dim * vol_dim * vol_dim; ++i) {
		max_value = MAX(max_value, fabs((double)vol_data[i]));
	}
	EXPECT_TRUE(max_value > 1.0e-6);

	vlx_test_free(&t);
}

UTEST(vlx, mol) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), MEGABYTES(64)));
	vlx_test_free(&t);
}

UTEST(vlx, scf_results) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/tq.scf.results.h5"), MEGABYTES(64)));
	vlx_test_free(&t);
}

UTEST(vlx, acro_rsp) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-xps.h5"), MEGABYTES(64)));
	vlx_test_free(&t);
}

// XPS is a delta-SCF property, not a response property, so it may coexist with any response type or
// with none. A file without it publishes nothing under vlx/xps - and an absent path, rather than an
// empty array, is what a consumer tests.
//
// Reads h2o.h5, which is a response calculation with no core-hole states, because acro-rsp.h5 is
// not in test_data - this test and vlx.acro_rsp both used to name it and both used to fail on that.
UTEST(vlx, file_without_xps_publishes_none) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(64)));

	EXPECT_FALSE(vlx_test_has(&t, STR_LIT("vlx/xps/ionization_energy")));
	EXPECT_FALSE(vlx_test_has(&t, STR_LIT("vlx/xps/element")));

	md_attribute_id_t ids[8];
	EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t.sys.attributes, STR_LIT("vlx/xps")), 0u);

	vlx_test_free(&t);
}

// XPS entries, published as one attribute per field of the record over a shared {C} index space.
// The per element grouping is NOT published: entries are laid out as contiguous runs of equal
// element, so a consumer scans vlx/xps/element for the run it wants. That derivation is what this
// pins, because it is the thing that breaks silently if the sort ever changes.
UTEST(vlx, acro_xps) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-xps.h5"), MEGABYTES(64)));

	const size_t count = vlx_test_count(&t, STR_LIT("vlx/xps/ionization_energy"));
	ASSERT_TRUE(count > 0);

	double* energy       = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* element      = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* contribution = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* atom_index   = (double*)md_alloc(t.alloc, sizeof(double) * count);

	ASSERT_EQ(vlx_test_series(energy,       count, &t, STR_LIT("vlx/xps/ionization_energy")), count);
	ASSERT_EQ(vlx_test_series(element,      count, &t, STR_LIT("vlx/xps/element")),           count);
	ASSERT_EQ(vlx_test_series(contribution, count, &t, STR_LIT("vlx/xps/contribution")),      count);
	ASSERT_EQ(vlx_test_series(atom_index,   count, &t, STR_LIT("vlx/xps/atom_index")),        count);

	// Sorted by (element, ionization energy), which is what makes one element's states a contiguous
	// run a consumer can find by scanning.
	for (size_t i = 1; i < count; ++i) {
		EXPECT_TRUE(element[i - 1] < element[i] ||
					(element[i - 1] == element[i] && energy[i - 1] <= energy[i] + 1.0e-12));
	}

	for (size_t i = 0; i < count; ++i) {
		EXPECT_TRUE(energy[i] > 0.0);
		EXPECT_TRUE(contribution[i] >= 0.0 && contribution[i] <= 1.0 + 1.0e-12);
		EXPECT_TRUE(atom_index[i] >= 0.0 && atom_index[i] < (double)t.sys.atom.count);
	}

	vlx_test_free(&t);
}

// 'orbital/alpha/density' is never stored: it is reconstructed on demand from the coefficient and
// occupation attributes published beside it. This is the check that the reconstruction is the matrix
// those two describe - recomputed here from the same published attributes, independently of the
// provider - and not merely a plausible one.
//
// NOT checked by tr(D S) against the electron count, which is the obvious test and a wrong one here:
// 'basis/overlap' is the Cartesian embedding of a spherical overlap (25 Cartesian AOs for 24
// spherical ones in this file), so it is rank deficient and its diagonal is not unity - tr(S) is 91,
// not 25. tr(D S) comes out 5.36 against 5 alpha electrons for both the reconstructed density AND
// the one the file stores, so the discrepancy is in the overlap conversion and nothing a density
// test can pin. See the AO CONVENTION block in md_gto.h.
UTEST(vlx, orbital_density_matches_coefficients_and_occupations) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(16)));

	const md_attribute_t* coeff = vlx_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
	const md_attribute_t* occup = vlx_test_attr(&t, STR_LIT("orbital/alpha/occupation"));
	const md_attribute_t* dens  = vlx_test_attr(&t, STR_LIT("orbital/alpha/density"));
	ASSERT_TRUE(coeff != NULL && occup != NULL && dens != NULL);

	// {M,A} coefficients, {M} occupations, {A,A} density - the three have to agree on M and A or the
	// reconstruction is indexing a matrix it does not belong to.
	ASSERT_EQ(coeff->format.rank, 2u);
	const size_t num_mo = coeff->format.shape[0];
	const size_t num_ao = coeff->format.shape[1];
	ASSERT_EQ(md_attribute_element_count(&occup->format), num_mo);
	ASSERT_EQ(dens->format.shape[0], (uint32_t)num_ao);
	ASSERT_EQ(dens->format.shape[1], (uint32_t)num_ao);

	EXPECT_EQ(dens->storage, MD_ATTRIBUTE_STORAGE_VIRTUAL);
	EXPECT_TRUE(md_attributes_data(&t.sys.attributes, dens->id, MD_ATTRIBUTE_TYPE_F64) == NULL);

	const size_t plane = num_ao * num_ao;
	double* C = (double*)md_alloc(t.alloc, sizeof(double) * num_mo * num_ao);
	double* occ = (double*)md_alloc(t.alloc, sizeof(double) * num_mo);
	double* D = (double*)md_alloc(t.alloc, sizeof(double) * plane);
	double* ref = (double*)md_alloc(t.alloc, sizeof(double) * plane);

	ASSERT_EQ(md_attribute_extract_f64(C, num_mo * num_ao, coeff, md_unit_none()), num_mo * num_ao);
	ASSERT_EQ(md_attribute_extract_f64(occ, num_mo, occup, md_unit_none()), num_mo);
	ASSERT_EQ(md_attribute_extract_f64(D, plane, dens, md_unit_none()), plane);

	// D = sum over molecular orbitals of occ_i * c_i c_i^T, which is the definition and not a
	// restatement of the provider: it reads the same two attributes any consumer would.
	MEMSET(ref, 0, sizeof(double) * plane);
	double occ_sum = 0.0;
	for (size_t mo = 0; mo < num_mo; ++mo) {
		occ_sum += occ[mo];
		if (occ[mo] == 0.0) continue;
		const double* c = C + mo * num_ao;
		for (size_t i = 0; i < num_ao; ++i) {
			for (size_t j = 0; j < num_ao; ++j) {
				ref[i * num_ao + j] += occ[mo] * c[i] * c[j];
			}
		}
	}

	// The occupations account for every alpha electron the file says the calculation had.
	EXPECT_NEAR(vlx_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), -1.0), occ_sum, 1.0e-9);

	double max_diff = 0.0, max_asym = 0.0, magnitude = 0.0;
	for (size_t i = 0; i < num_ao; ++i) {
		for (size_t j = 0; j < num_ao; ++j) {
			max_diff  = MAX(max_diff,  fabs(D[i * num_ao + j] - ref[i * num_ao + j]));
			max_asym  = MAX(max_asym,  fabs(D[i * num_ao + j] - D[j * num_ao + i]));
			magnitude = MAX(magnitude, fabs(D[i * num_ao + j]));
		}
	}
	EXPECT_NEAR(0.0, max_diff, 1.0e-12);

	// The density evaluation path packs only the upper triangle, so an asymmetric matrix would be
	// silently half consumed.
	EXPECT_NEAR(0.0, max_asym, 1.0e-12);

	// Guards both against the all zeros case, which every difference above would also satisfy.
	EXPECT_TRUE(magnitude > 1.0e-8);

	vlx_test_free(&t);
}

// The combined spin densities, which are derivations over derivations: alpha and beta are each
// computed on demand, and these read both. That is legal because the graph stays acyclic, and it is
// the case that broke - beta is an ALIAS of alpha in a restricted calculation, and reading an alias
// of a computed attribute used to return nothing at all, silently.
//
// A restricted calculation gives the invariants for free and needs no reference values: the total
// is exactly twice alpha, and the difference is exactly zero.
UTEST(vlx, combined_spin_densities_h2o) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(16)));

	EXPECT_TRUE(str_eq(vlx_test_string(&t, STR_LIT("vlx/scf/type")), STR_LIT("restricted")));

	const md_attribute_t* alpha = vlx_test_attr(&t, STR_LIT("orbital/alpha/density"));
	ASSERT_TRUE(alpha != NULL);
	ASSERT_EQ(alpha->format.rank, 2u);
	const size_t num_ao = alpha->format.shape[0];
	const size_t plane  = num_ao * num_ao;

	str_t paths[4] = {
		STR_LIT("orbital/alpha/density"),
		STR_LIT("orbital/beta/density"),
		STR_LIT("orbital/total/density"),
		STR_LIT("orbital/difference/density"),
	};

	double* mat[4] = {0};
	for (int i = 0; i < 4; ++i) {
		const md_attribute_t* a = vlx_test_attr(&t, paths[i]);
		ASSERT_TRUE(a != NULL);
		ASSERT_EQ(a->format.rank, 2u);
		EXPECT_EQ(a->format.shape[0], (uint32_t)num_ao);
		EXPECT_EQ(a->format.shape[1], (uint32_t)num_ao);

		mat[i] = (double*)md_alloc(t.alloc, sizeof(double) * plane);
		ASSERT_EQ(md_attribute_extract_f64(mat[i], plane, a, md_unit_none()), plane);
	}

	// Restricted: beta is a second NAME for alpha, not a second reconstruction.
	const md_attribute_t* beta = vlx_test_attr(&t, STR_LIT("orbital/beta/density"));
	EXPECT_EQ(beta->storage, MD_ATTRIBUTE_STORAGE_ALIAS);
	EXPECT_TRUE(md_attribute_same_data(beta, alpha));

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

	vlx_test_free(&t);
}

// The transition densities: attachment, detachment and their difference. There is nothing to
// compare them against - the file carries no reference matrices - so what is checked here is the
// MECHANICS and the invariants that hold whatever the numbers are.
//
// They are the one place in the table where a slice is load bearing rather than a convenience. The
// attribute is {S,A,A} and every plane costs a full reconstruction from the response eigenvectors,
// so asking for one state has to mean reconstructing one state. Both paths are exercised below.
static void check_transition_density_attributes(int* utest_result, str_t file) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, file, MEGABYTES(64)));

	const size_t num_states = vlx_test_count(&t, STR_LIT("vlx/rsp/oscillator_strength"));
	ASSERT_TRUE(num_states > 0);

	str_t paths[3] = {
		STR_LIT("vlx/rsp/transition_density/attachment"),
		STR_LIT("vlx/rsp/transition_density/detachment"),
		STR_LIT("vlx/rsp/transition_density/difference"),
	};

	const md_attribute_t* attr[3] = {0};
	for (int i = 0; i < 3; ++i) {
		attr[i] = vlx_test_attr(&t, paths[i]);
		ASSERT_TRUE(attr[i] != NULL);

		// {S,A,A}: state outermost, then the AO x AO matrix. Square is load bearing downstream -
		// the GL and GPU density paths pack only the upper triangle.
		EXPECT_EQ(attr[i]->storage, MD_ATTRIBUTE_STORAGE_VIRTUAL);
		EXPECT_EQ(attr[i]->format.type, MD_ATTRIBUTE_TYPE_F64);
		ASSERT_EQ(attr[i]->format.rank, 3u);
		EXPECT_EQ(attr[i]->format.shape[0], (uint32_t)num_states);
		EXPECT_EQ(attr[i]->format.shape[1], attr[i]->format.shape[2]);

		// A virtual attribute hands out no resident storage to write into.
		EXPECT_TRUE(md_attributes_data(&t.sys.attributes, attr[i]->id, MD_ATTRIBUTE_TYPE_F64) == NULL);
	}

	const size_t num_ao = attr[0]->format.shape[1];
	ASSERT_TRUE(num_ao > 0);

	// The AO axis is the one the coefficients and the overlap live on; if these disagreed the
	// matrices would be reconstructed against a basis they do not belong to.
	const md_attribute_t* coeff = vlx_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
	ASSERT_TRUE(coeff != NULL);
	EXPECT_EQ(coeff->format.shape[1], (uint32_t)num_ao);

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

		mat[i] = (double*)md_alloc(t.alloc, sizeof(double) * plane);
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

	// Symmetric by construction, and the GL/GPU density paths read only the upper triangle - so an
	// asymmetric matrix would be silently half consumed.
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

		double* all = (double*)md_alloc(t.alloc, sizeof(double) * total);
		ASSERT_EQ(md_attribute_extract_f64(all, total, attr[i], md_unit_none()), total);

		for (size_t v = 0; v < plane; ++v) {
			ASSERT_NEAR(mat[i][v], all[v], 1.0e-12);
		}
	}

	// A state past the end selects nothing rather than reading past the array.
	const md_attribute_slice_t past_end = md_attribute_slice_1((uint32_t)num_states);
	EXPECT_EQ(md_attribute_slice_count(attr[0], &past_end), 0u);
	EXPECT_EQ(md_attribute_extract_slice_f64(mat[0], plane, attr[0], &past_end, md_unit_none()), 0u);

	vlx_test_free(&t);
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

// The NTO coefficient vectors are resident and the transition densities are computed, but both are
// built from the same solution vectors - so the leading particle NTO of a state has to be a vector
// in the same AO space, of the same length, as everything else the table publishes over AOs. This is
// what catches the two coming apart on the lambda padding axis.
UTEST(vlx, nto_coefficients_share_the_ao_axis) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(32)));

	const md_attribute_t* lambda   = vlx_test_attr(&t, STR_LIT("vlx/rsp/nto/lambda"));
	const md_attribute_t* particle = vlx_test_attr(&t, STR_LIT("vlx/rsp/nto/particle/coefficient"));
	const md_attribute_t* hole     = vlx_test_attr(&t, STR_LIT("vlx/rsp/nto/hole/coefficient"));
	const md_attribute_t* coeff    = vlx_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
	ASSERT_TRUE(lambda != NULL && particle != NULL && hole != NULL && coeff != NULL);

	// {S,Lmax,A} against the weights' {S,Lmax} and the MO coefficients' {M,A}.
	ASSERT_EQ(particle->format.rank, 3u);
	EXPECT_EQ(particle->format.shape[0], lambda->format.shape[0]);
	EXPECT_EQ(particle->format.shape[1], lambda->format.shape[1]);
	EXPECT_EQ(particle->format.shape[2], coeff->format.shape[1]);

	EXPECT_EQ(hole->format.shape[0], particle->format.shape[0]);
	EXPECT_EQ(hole->format.shape[1], particle->format.shape[1]);
	EXPECT_EQ(hole->format.shape[2], particle->format.shape[2]);

	// The leading pair of the first state carries the weight the file names, so its vectors are not
	// the zero pad a short row is filled with.
	const size_t num_ao = particle->format.shape[2];
	double* vec = (double*)md_alloc(t.alloc, sizeof(double) * num_ao);

	const md_attribute_slice_t leading = md_attribute_slice_2(0, 0);
	ASSERT_EQ(md_attribute_extract_slice_f64(vec, num_ao, particle, &leading, md_unit_none()), num_ao);

	double norm = 0.0;
	for (size_t i = 0; i < num_ao; ++i) {
		norm += vec[i] * vec[i];
	}
	EXPECT_TRUE(norm > 1.0e-8);

	vlx_test_free(&t);
}
