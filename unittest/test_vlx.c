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

	EXPECT_EQ(0.0, qm_test_scalar(&t, STR_LIT("vlx/molecular_charge"), -1.0));
	EXPECT_EQ(1.0, qm_test_scalar(&t, STR_LIT("vlx/spin_multiplicity"), -1.0));
	EXPECT_EQ(41.0, qm_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), -1.0));
	EXPECT_EQ(41.0, qm_test_scalar(&t, STR_LIT("vlx/electron_count/beta"), -1.0));

	EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("vlx/basis_set")), STR_LIT("DEF2-SVP")));

	// The QM geometry, in Angstrom, as the calculation was run at.
	const md_attribute_t* coord = qm_test_attr(&t, STR_LIT("qm/atom/coordinate"));
	ASSERT_TRUE(coord != NULL);
	ASSERT_EQ(coord->format.components, 3u);
	ASSERT_EQ(md_attribute_value_count(&coord->format), 26u);

	double xyz[3 * 26] = {0};
	ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_attribute_slice_all(), md_unit_none()), ARRAY_SIZE(xyz));
	EXPECT_NEAR(-3.259400000000, xyz[0], 1.0e-5);
	EXPECT_NEAR( 0.145200000000, xyz[1], 1.0e-5);
	EXPECT_NEAR(-0.048400000000, xyz[2], 1.0e-5);

	// The system's own state is that geometry too, narrowed to float.
	EXPECT_NEAR(xyz[0], (double)t.state.xyz[0].x, 1.0e-5);
	EXPECT_NEAR(xyz[1], (double)t.state.xyz[0].y, 1.0e-5);
	EXPECT_NEAR(xyz[2], (double)t.state.xyz[0].z, 1.0e-5);

	const size_t num_iter = qm_test_count(&t, STR_LIT("vlx/scf/history/energy"));
	ASSERT_TRUE(num_iter > 0);

	double* energy = (double*)md_alloc(t.alloc, sizeof(double) * num_iter);
	ASSERT_EQ(qm_test_series(energy, num_iter, &t, STR_LIT("vlx/scf/history/energy")), num_iter);
	EXPECT_NEAR(ref_ener_tot, energy[num_iter - 1], 1.0e-5);

	qm_test_free(&t);
}

// The facts about a calculation that are TEXT and nothing else, plus the counts that used to be
// reachable only by asking the reader. If one of these stops being published, nothing downstream can
// tell a restricted run from an unrestricted one except by guessing from which arrays are shared.
UTEST(vlx, run_description_is_published) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(16)));

	EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("vlx/scf/type")), STR_LIT("restricted")));
	EXPECT_TRUE(qm_test_scalar(&t, STR_LIT("vlx/spin_multiplicity"), 0.0) == 1.0);
	EXPECT_TRUE(qm_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), 0.0) > 0.0);
	EXPECT_TRUE(qm_test_scalar(&t, STR_LIT("vlx/nuclear_repulsion_energy"), 0.0) > 0.0);

	// A response calculation, so the type is there and names which one.
	EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("vlx/rsp/type")), STR_LIT("linear")));

	// No geometry optimisation in this file, so nothing under vlx/opt. An absent path is how a
	// consumer learns a block is missing - not a zero it would have to interpret.
	EXPECT_FALSE(qm_test_has(&t, STR_LIT("vlx/opt/type")));
	EXPECT_FALSE(qm_test_has(&t, STR_LIT("vlx/opt/state_index")));
	EXPECT_FALSE(qm_test_has(&t, STR_LIT("vlx/opt/irc_ts_index")));

	qm_test_free(&t);
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
		const md_attribute_t* a = qm_test_attr(&t, STR_LIT("vlx/rsp/nto/lambda"));
		ASSERT_TRUE(a != NULL);
		ASSERT_EQ(a->format.rank, 2u);

		const size_t num_lambdas = a->format.shape[1];
		ASSERT_TRUE(num_lambdas >= cases[c].count);

		double* lambdas = (double*)md_alloc(t.alloc, sizeof(double) * num_lambdas);
		ASSERT_EQ(qm_test_row(lambdas, num_lambdas, &t, STR_LIT("vlx/rsp/nto/lambda"), 0), num_lambdas);

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

		qm_test_free(&t);
	}
}

// What sampling a molecular orbital over a grid looks like end to end, with no reader in the loop:
// load the file into a system, rebuild the basis from the format neutral basis/ attributes, take one
// MO's coefficient row, evaluate. This doubles as the worked example for the whole interface.
UTEST(vlx, minimal_example) {
	vlx_test_t t = {0};
	str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5");
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
	ASSERT_TRUE(qm_test_basis(&basis, &t));

	size_t num_gtos = md_gto_pgto_count(&basis);
	md_gto_t* gtos = (md_gto_t*)md_alloc(arena, sizeof(md_gto_t) * num_gtos);

	const size_t num_aos = md_gto_basis_num_ao(&basis);
	double* mo_coeffs = (double*)md_alloc(arena, sizeof(double) * num_aos);
	ASSERT_EQ(qm_test_row(mo_coeffs, num_aos, &t, STR_LIT("orbital/alpha/coefficient"), mo_idx), num_aos);

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

	qm_test_free(&t);
}

UTEST(vlx, mol) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), MEGABYTES(64)));
	qm_test_free(&t);
}

UTEST(vlx, scf_results) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/tq.scf.results.h5"), MEGABYTES(64)));
	qm_test_free(&t);
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

	EXPECT_FALSE(qm_test_has(&t, STR_LIT("vlx/xps/ionization_energy")));
	EXPECT_FALSE(qm_test_has(&t, STR_LIT("vlx/xps/element")));

	md_attribute_id_t ids[8];
	EXPECT_EQ(md_attributes_query(ids, ARRAY_SIZE(ids), &t.sys.attributes, STR_LIT("vlx/xps")), 0u);

	qm_test_free(&t);
}

// XPS entries, published as one attribute per field of the record over a shared {C} index space.
// The per element grouping is NOT published: entries are laid out as contiguous runs of equal
// element, so a consumer scans vlx/xps/element for the run it wants. That derivation is what this
// pins, because it is the thing that breaks silently if the sort ever changes.
UTEST(vlx, acro_xps) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-xps.h5"), MEGABYTES(64)));

	const size_t count = qm_test_count(&t, STR_LIT("vlx/xps/ionization_energy"));
	ASSERT_TRUE(count > 0);

	double* energy       = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* element      = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* contribution = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* atom_index   = (double*)md_alloc(t.alloc, sizeof(double) * count);

	ASSERT_EQ(qm_test_series(energy,       count, &t, STR_LIT("vlx/xps/ionization_energy")), count);
	ASSERT_EQ(qm_test_series(element,      count, &t, STR_LIT("vlx/xps/element")),           count);
	ASSERT_EQ(qm_test_series(contribution, count, &t, STR_LIT("vlx/xps/contribution")),      count);
	ASSERT_EQ(qm_test_series(atom_index,   count, &t, STR_LIT("vlx/xps/atom_index")),        count);

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

	qm_test_free(&t);
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

	const md_attribute_t* coeff = qm_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
	const md_attribute_t* occup = qm_test_attr(&t, STR_LIT("orbital/alpha/occupation"));
	const md_attribute_t* dens  = qm_test_attr(&t, STR_LIT("orbital/alpha/density"));
	ASSERT_TRUE(coeff != NULL && occup != NULL && dens != NULL);

	// {M,A} coefficients, {M} occupations, {A,A} density - the three have to agree on M and A or the
	// reconstruction is indexing a matrix it does not belong to.
	ASSERT_EQ(coeff->format.rank, 2u);
	const size_t num_mo = coeff->format.shape[0];
	const size_t num_ao = coeff->format.shape[1];
	ASSERT_EQ(md_attribute_element_count(&occup->format), num_mo);
	ASSERT_EQ(dens->format.shape[0], (uint32_t)num_ao);
	ASSERT_EQ(dens->format.shape[1], (uint32_t)num_ao);

	EXPECT_TRUE(md_attribute_is_virtual(dens));
	EXPECT_TRUE(md_attributes_data(&t.sys.attributes, dens->id, MD_ATTRIBUTE_TYPE_F64) == NULL);

	const size_t plane = num_ao * num_ao;
	double* C = (double*)md_alloc(t.alloc, sizeof(double) * num_mo * num_ao);
	double* occ = (double*)md_alloc(t.alloc, sizeof(double) * num_mo);
	double* D = (double*)md_alloc(t.alloc, sizeof(double) * plane);
	double* ref = (double*)md_alloc(t.alloc, sizeof(double) * plane);

	ASSERT_EQ(md_attribute_extract_f64(C, num_mo * num_ao, coeff, md_attribute_slice_all(), md_unit_none()), num_mo * num_ao);
	ASSERT_EQ(md_attribute_extract_f64(occ, num_mo, occup, md_attribute_slice_all(), md_unit_none()), num_mo);
	ASSERT_EQ(md_attribute_extract_f64(D, plane, dens, md_attribute_slice_all(), md_unit_none()), plane);

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
	EXPECT_NEAR(qm_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), -1.0), occ_sum, 1.0e-9);

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

	qm_test_free(&t);
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

	EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("vlx/scf/type")), STR_LIT("restricted")));

	const md_attribute_t* alpha = qm_test_attr(&t, STR_LIT("orbital/alpha/density"));
	ASSERT_TRUE(alpha != NULL);
	ASSERT_EQ(alpha->format.rank, 2u);
	const size_t num_ao = alpha->format.shape[0];
	const size_t plane  = num_ao * num_ao;

	str_t paths[4] = {
		STR_INIT("orbital/alpha/density"),
		STR_INIT("orbital/beta/density"),
		STR_INIT("orbital/total/density"),
		STR_INIT("orbital/difference/density"),
	};

	double* mat[4] = {0};
	for (int i = 0; i < 4; ++i) {
		const md_attribute_t* a = qm_test_attr(&t, paths[i]);
		ASSERT_TRUE(a != NULL);
		ASSERT_EQ(a->format.rank, 2u);
		EXPECT_EQ(a->format.shape[0], (uint32_t)num_ao);
		EXPECT_EQ(a->format.shape[1], (uint32_t)num_ao);

		mat[i] = (double*)md_alloc(t.alloc, sizeof(double) * plane);
		ASSERT_EQ(md_attribute_extract_f64(mat[i], plane, a, md_attribute_slice_all(), md_unit_none()), plane);
	}

	// Restricted: beta is a second NAME for alpha, not a second reconstruction.
	const md_attribute_t* beta = qm_test_attr(&t, STR_LIT("orbital/beta/density"));
	EXPECT_TRUE(md_attribute_is_alias(beta));
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

	qm_test_free(&t);
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

	const size_t num_states = qm_test_count(&t, STR_LIT("vlx/rsp/oscillator_strength"));
	ASSERT_TRUE(num_states > 0);

	str_t paths[3] = {
		STR_INIT("vlx/rsp/transition_density/attachment"),
		STR_INIT("vlx/rsp/transition_density/detachment"),
		STR_INIT("vlx/rsp/transition_density/difference"),
	};

	const md_attribute_t* attr[3] = {0};
	for (int i = 0; i < 3; ++i) {
		attr[i] = qm_test_attr(&t, paths[i]);
		ASSERT_TRUE(attr[i] != NULL);

		// {S,A,A}: state outermost, then the AO x AO matrix. Square is load bearing downstream -
		// the GL and GPU density paths pack only the upper triangle.
		EXPECT_TRUE(md_attribute_is_virtual(attr[i]));
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
	const md_attribute_t* coeff = qm_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
	ASSERT_TRUE(coeff != NULL);
	EXPECT_EQ(coeff->format.shape[1], (uint32_t)num_ao);

	const size_t plane = num_ao * num_ao;
	double* mat[3] = {0};

	// One state at a time, which is the shape a representation asks in.
	for (int i = 0; i < 3; ++i) {
		const md_attribute_slice_t slice = md_attribute_slice_1(0);

		md_attribute_format_t sliced = {0};
		ASSERT_TRUE(md_attribute_slice_format(&sliced, attr[i], slice));
		EXPECT_EQ(sliced.rank, 2u);
		EXPECT_EQ(sliced.shape[0], (uint32_t)num_ao);
		EXPECT_EQ(sliced.shape[1], (uint32_t)num_ao);
		ASSERT_EQ(md_attribute_slice_count(attr[i], slice), plane);

		mat[i] = (double*)md_alloc(t.alloc, sizeof(double) * plane);
		ASSERT_EQ(md_attribute_extract_f64(mat[i], plane, attr[i], slice, md_unit_none()), plane);
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
		ASSERT_EQ(md_attribute_extract_f64(all, total, attr[i], md_attribute_slice_all(), md_unit_none()), total);

		for (size_t v = 0; v < plane; ++v) {
			ASSERT_NEAR(mat[i][v], all[v], 1.0e-12);
		}
	}

	// A state past the end selects nothing rather than reading past the array.
	const md_attribute_slice_t past_end = md_attribute_slice_1((uint32_t)num_states);
	EXPECT_EQ(md_attribute_slice_count(attr[0], past_end), 0u);
	EXPECT_EQ(md_attribute_extract_f64(mat[0], plane, attr[0], past_end, md_unit_none()), 0u);

	qm_test_free(&t);
}

// Two files on purpose: h2o is one excited state in a small basis, amide is 110 atomic orbitals -
// the size the "provider wrote 0 of 12100" failure was reported at, and the one where a per state
// reconstruction is expensive enough that the slice has to mean what it says.
// basis/overlap and the coefficients beside it have to describe ONE basis, and these are the two
// numbers that say whether they do.
//
// A REGRESSION TEST WITH A DATE ON IT. Until 2026-09-15 the published overlap was the file's own
// spherical S pushed through md_gto_sph_to_cart_matrix, which computes T^T S T. That is the right
// transform for a DENSITY, which is built from coefficients, and the wrong one for an OVERLAP,
// which is built from the basis functions themselves. Nothing asserted either number, so nothing
// noticed: on this file the orthonormality below was 4.1e+02 and the electron count 10.72 for a ten
// electron molecule. Both look like data. It is integrated from the published basis now - see
// md_qm_publish_overlap, which also says why no other conversion would have worked.
UTEST(vlx, overlap_describes_the_same_basis_as_the_coefficients) {
	const char* files[] = {
		MD_UNITTEST_DATA_DIR "/vlx/h2o.h5",
		MD_UNITTEST_DATA_DIR "/vlx/mol.h5",
	};
	for (size_t i = 0; i < ARRAY_SIZE(files); ++i) {
		vlx_test_t t = {0};
		ASSERT_TRUE(vlx_test_load(&t, str_from_cstr(files[i]), MEGABYTES(256)));

		// The molecular orbitals are orthonormal against it, or the two disagree about the basis.
		EXPECT_LT(qm_test_orthonormality(&t, STR_LIT("orbital/alpha/coefficient")), 1.0e-5);

		// ...and tr(D S) is the electron count the file itself states. This is the assertion the
		// occupation convention shows up in, which orthonormality is completely blind to.
		const double alpha = qm_test_scalar(&t, STR_LIT("vlx/electron_count/alpha"), -1.0);
		const double beta  = qm_test_scalar(&t, STR_LIT("vlx/electron_count/beta"),  -1.0);
		ASSERT_TRUE(alpha > 0.0);
		ASSERT_TRUE(beta  > 0.0);

		EXPECT_NEAR(alpha,        qm_test_electron_count(&t, STR_LIT("orbital/alpha/density")), 1.0e-4);
		EXPECT_NEAR(alpha + beta, qm_test_electron_count(&t, STR_LIT("orbital/total/density")), 1.0e-4);
		EXPECT_NEAR(0.0,          qm_test_electron_count(&t, STR_LIT("orbital/difference/density")), 1.0e-9);

		qm_test_free(&t);
	}
}

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

	const md_attribute_t* lambda   = qm_test_attr(&t, STR_LIT("vlx/rsp/nto/lambda"));
	const md_attribute_t* particle = qm_test_attr(&t, STR_LIT("vlx/rsp/nto/particle/coefficient"));
	const md_attribute_t* hole     = qm_test_attr(&t, STR_LIT("vlx/rsp/nto/hole/coefficient"));
	const md_attribute_t* coeff    = qm_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
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
	ASSERT_EQ(md_attribute_extract_f64(vec, num_ao, particle, leading, md_unit_none()), num_ao);

	double norm = 0.0;
	for (size_t i = 0; i < num_ao; ++i) {
		norm += vec[i] * vec[i];
	}
	EXPECT_TRUE(norm > 1.0e-8);

	qm_test_free(&t);
}

// ---------------------------------------------------------------------------
// POLARIZABLE EMBEDDING
//
// No .h5 of an embedding run is checked in, so these make one: a copy of h2o.h5 whose potfile_text
// holds a potential the way VeloxChem stores it - PyFraME JSON as write_pe_jsonfile writes it
// (json.dump with indent 4), read back through '\n'.join(f.readlines()), which doubles every
// newline, into a fixed length string dataset of shape (1,), as np.bytes_([text]) makes it. The copy
// is written beside h2o.h5, where its basis set is, and removed again as soon as it has been read -
// before any assertion that could end the test early.
// ---------------------------------------------------------------------------

#include <hdf5.h>
#include <stdio.h>	// remove
#include <core/md_os.h>
#include <core/md_str_builder.h>
#include <md_pot.h>
#include <md_util.h>
#include "system_invariants.h"

#define VLX_PE_DIR MD_UNITTEST_DATA_DIR "/vlx/"

static bool vlx_test_write_file(str_t path, const void* data, size_t size) {
	md_file_t out = {0};
	if (!md_file_open(&out, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
	const bool ok = md_file_write(out, data, size) == size;
	md_file_close(&out);
	return ok;
}

// A copy of h2o.h5 with 'text' as its potfile_text
static bool vlx_test_write_pe_h5(str_t dst, str_t text) {
	md_allocator_i* heap = md_get_heap_allocator();
	md_file_t in = {0};
	if (!md_file_open(&in, STR_LIT(VLX_PE_DIR "h2o.h5"), MD_FILE_READ)) return false;
	const size_t size = (size_t)md_file_size(in);
	void* bytes = md_alloc(heap, size);
	const bool read = md_file_read(in, bytes, size) == size;
	md_file_close(&in);
	const bool written = read && vlx_test_write_file(dst, bytes, size);
	md_free(heap, bytes, size);
	if (!written) return false;

	char path[1024];
	str_copy_to_char_buf(path, sizeof(path), dst);
	hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0) return false;

	const size_t len = text.len ? text.len : 1;
	char* buf = md_alloc(heap, len);
	MEMSET(buf, 0, len);
	if (text.len) MEMCPY(buf, text.ptr, text.len);

	const bool deleted = H5Lexists(file, "potfile_text", H5P_DEFAULT) <= 0 || H5Ldelete(file, "potfile_text", H5P_DEFAULT) >= 0;
	hid_t type = H5Tcopy(H5T_C_S1);
	H5Tset_size(type, len);
	H5Tset_strpad(type, H5T_STR_NULLPAD);
	const hsize_t dims[1] = { 1 };
	hid_t space = H5Screate_simple(1, dims, NULL);
	hid_t dset  = deleted ? H5Dcreate2(file, "potfile_text", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
	const bool ok = dset >= 0 && H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) >= 0;
	if (dset >= 0) H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(type);
	H5Fclose(file);
	md_free(heap, buf, len);
	return ok;
}

// What f.readlines() joined by '\n' makes of a text: every newline doubled
static str_t vlx_test_as_read_back(str_t json, md_allocator_i* alloc) {
	md_strb_t sb = md_strb_create(alloc);
	for (size_t i = 0; i < json.len; ++i) {
		md_strb_push_char(&sb, json.ptr[i]);
		if (json.ptr[i] == '\n') md_strb_push_char(&sb, '\n');
	}
	return str_copy(md_strb_to_str(sb), alloc);
}

// The k-th row of a fragment type in a per type section of a .pot, or NULL
static const md_pot_charge_t* vlx_test_pot_charge(const md_pot_t* pot, str_t type, size_t k) {
	for (size_t i = 0; i < pot->num_charges; ++i) {
		if (str_eq(pot->charges[i].fragment_name, type) && k-- == 0) return &pot->charges[i];
	}
	return NULL;
}

static const md_pot_polarizability_t* vlx_test_pot_polarizability(const md_pot_t* pot, str_t type, size_t k) {
	for (size_t i = 0; i < pot->num_polarizabilities; ++i) {
		if (str_eq(pot->polarizabilities[i].fragment_name, type) && k-- == 0) return &pot->polarizabilities[i];
	}
	return NULL;
}

// A .pot as write_pe_jsonfile turns it into JSON: fragments in residue number order, numbered from 1
// and named for their residue; each site with the charge and polarizability rows of its type by its
// position in the fragment, zero where the type has none; coordinates in bohr; the polarizability
// behind the (0,0) and (0,1) blocks of order [1,1]; the exclusions the fragment's own sites.
static str_t vlx_test_pe_json(const md_pot_t* pot, md_allocator_i* alloc) {
	const double bohr_in_angstroms = 0.529177210903;	// VeloxChem's
	const double to_angstrom = md_pot_unit_to_angstrom(pot->unit);

	int max_id = 0;
	for (size_t i = 0; i < pot->num_sites; ++i) max_id = MAX(max_id, pot->sites[i].fragment_id);

	md_strb_t sb = md_strb_create(alloc);
	md_strb_push_cstr(&sb,
		"{\n"
		"    \"quantum_subsystems\": [\n"
		"        {\n"
		"            \"nuclei\": []\n"
		"        }\n"
		"    ],\n"
		"    \"classical_subsystems\": [\n"
		"        {\n"
		"            \"classical_fragments\": [");
	size_t num_frag = 0, num_atoms = 0;
	for (int id = 0; id <= max_id; ++id) {
		size_t n = 0;
		str_t type = {0};
		for (size_t i = 0; i < pot->num_sites; ++i) {
			if (pot->sites[i].fragment_id != id) continue;
			if (n++ == 0) type = pot->sites[i].fragment_name;
		}
		if (n == 0) continue;
		num_frag += 1;
		md_strb_fmt(&sb, "%s\n                {\n                    \"index\": %zu,\n                    \"name\": \"" STR_FMT "\",\n                    \"atoms\": [",
			num_frag > 1 ? "," : "", num_frag, STR_ARG(type));
		const size_t first = num_atoms + 1;
		size_t k = 0;
		for (size_t i = 0; i < pot->num_sites; ++i) {
			const md_pot_site_t* s = &pot->sites[i];
			if (s->fragment_id != id) continue;
			const md_pot_charge_t* q = vlx_test_pot_charge(pot, type, k);
			const md_pot_polarizability_t* p = vlx_test_pot_polarizability(pot, type, k);
			const double* a = p ? p->alpha : (const double[6]){0};
			md_strb_fmt(&sb, "%s\n                        {\n"
				"                            \"index\": %zu,\n"
				"                            \"element\": \"%s\",\n"
				"                            \"coordinate\": [\n"
				"                                %.17g,\n"
				"                                %.17g,\n"
				"                                %.17g\n"
				"                            ],\n"
				"                            \"multipoles\": {\n"
				"                                \"elements\": [\n"
				"                                    %.17g\n"
				"                                ]\n"
				"                            },\n"
				"                            \"exclusions\": [",
				k ? "," : "", num_atoms + 1, s->element,
				s->coord[0] * to_angstrom / bohr_in_angstroms, s->coord[1] * to_angstrom / bohr_in_angstroms, s->coord[2] * to_angstrom / bohr_in_angstroms,
				q ? q->charge : 0.0);
			for (size_t e = 0; e < n; ++e) {
				md_strb_fmt(&sb, "%s\n                                %zu", e ? "," : "", first + e);
			}
			md_strb_fmt(&sb, "\n                            ],\n"
				"                            \"polarizabilities\": {\n"
				"                                \"elements\": [\n"
				"                                    0.0,\n                                    0.0,\n                                    0.0,\n                                    0.0,\n"
				"                                    %.17g,\n                                    %.17g,\n                                    %.17g,\n"
				"                                    %.17g,\n                                    %.17g,\n                                    %.17g\n"
				"                                ],\n"
				"                                \"order\": [\n                                    1,\n                                    1\n                                ]\n"
				"                            }\n"
				"                        }",
				a[0], a[1], a[2], a[3], a[4], a[5]);
			num_atoms += 1;
			k += 1;
		}
		md_strb_push_cstr(&sb, "\n                    ]\n                }");
	}
	md_strb_push_cstr(&sb, "\n            ]\n        }\n    ]\n}");
	return vlx_test_as_read_back(md_strb_to_str(sb), alloc);
}

// Read on the bits: the library is built with fast math, where NAN does not compare as itself
static bool vlx_test_absent(double v) {
	uint64_t u;
	memcpy(&u, &v, sizeof(u));
	return (u & 0x7fffffffffffffffull) > 0x7ff0000000000000ull;
}

// The reference potential as VeloxChem embeds it: 50 waters by residue number, 39 of them
// polarizable (HOH_pe) and 11 not (HOH_npe), interleaved.
UTEST(vlx, pe_environment_becomes_atoms_of_the_system) {
	md_allocator_i* heap = md_get_heap_allocator();
	md_pot_t pot = {0};
	ASSERT_TRUE(md_pot_parse_file(&pot, STR_LIT(MD_UNITTEST_DATA_DIR "/pot/water_pe_npe.pot"), heap));
	const str_t json = vlx_test_pe_json(&pot, heap);
	md_pot_free(&pot, heap);

	const str_t h5 = STR_LIT(VLX_PE_DIR "unittest_pe_reference.h5");
	const bool written = vlx_test_write_pe_h5(h5, json);
	str_free(json, heap);
	ASSERT_TRUE(written);
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MEGABYTES(64));
	remove(h5.ptr);
	ASSERT_TRUE(loaded);

	const size_t num_qm = 3, num_mm = 150, num_atoms = num_qm + num_mm;
	ASSERT_EQ(num_atoms, t.sys.atom.count);
	ASSERT_EQ(num_atoms, t.state.num_atoms);

	// The QM atom domain is untouched, and its atoms are still the system's first
	EXPECT_EQ(num_qm, qm_test_count(&t, STR_LIT("qm/atom/atomic_number")));
	EXPECT_FALSE(qm_test_has(&t, STR_LIT("qm/atom/system_index")));
	EXPECT_EQ(8, md_atom_atomic_number(&t.sys.atom, 0));
	EXPECT_EQ(1, md_atom_atomic_number(&t.sys.atom, 1));

	// The sites follow, named for their elements - the JSON keeps no atom names - and back in
	// Angstrom from the bohr it stores
	EXPECT_EQ(8, md_atom_atomic_number(&t.sys.atom, num_qm));
	EXPECT_EQ(1, md_atom_atomic_number(&t.sys.atom, num_qm + 2));
	EXPECT_TRUE(str_eq(md_atom_name(&t.sys.atom, num_qm),     STR_LIT("O")));
	EXPECT_TRUE(str_eq(md_atom_name(&t.sys.atom, num_qm + 2), STR_LIT("H")));
	EXPECT_NEAR(-5.672, t.state.xyz[num_qm].x, 1e-5);
	EXPECT_NEAR( 2.390, t.state.xyz[num_qm].y, 1e-5);
	EXPECT_NEAR(-4.911, t.state.xyz[num_qm].z, 1e-5);
	// Residue 50 is last in the JSON, though not in the .pot
	EXPECT_NEAR(-8.589, t.state.xyz[num_atoms - 1].x, 1e-5);
	EXPECT_NEAR( 1.786, t.state.xyz[num_atoms - 1].z, 1e-5);

	// One component per fragment, named for its residue and numbered by its fragment number, after
	// one for the QM region - components cover every atom or none
	ASSERT_EQ(1u + 50u, t.sys.component.count);
	EXPECT_VALID_SYSTEM(&t.sys);
	EXPECT_TRUE(str_eq(md_component_name(&t.sys.component, 0), STR_LIT("QM")));
	EXPECT_EQ(0, md_system_component_find_by_atom_idx(&t.sys, 0));
	EXPECT_EQ(1, md_system_component_find_by_atom_idx(&t.sys, num_qm));
	const md_urange_t first = md_system_component_atom_range(&t.sys, 1);
	EXPECT_EQ(num_qm, first.beg);
	EXPECT_EQ(num_qm + 3, first.end);
	for (size_t c = 1; c < t.sys.component.count; ++c) {
		EXPECT_TRUE(str_eq(md_component_name(&t.sys.component, c), STR_LIT("HOH")));
		EXPECT_EQ((int)c, md_component_seq_id(&t.sys.component, c));
		EXPECT_EQ(MD_COMPONENT_KIND_WATER, md_system_component_kind(&t.sys, c));
	}

	// The parameters, per atom over the system, and nothing at all for the QM atoms. Residue 8 is
	// the first HOH_npe: its charges are its own, and it has no polarizabilities
	double q[153], a[153];
	ASSERT_EQ(num_atoms, qm_test_series(q, num_atoms, &t, STR_LIT("atom/charge")));
	ASSERT_EQ(num_atoms, qm_test_series(a, num_atoms, &t, STR_LIT("atom/polarizability")));
	for (size_t i = 0; i < num_qm; ++i) {
		EXPECT_TRUE(vlx_test_absent(q[i]));
		EXPECT_TRUE(vlx_test_absent(a[i]));
	}
	EXPECT_NEAR(-0.67444, q[num_qm + 0], 1e-12);
	EXPECT_NEAR( 0.33722, q[num_qm + 1], 1e-12);
	EXPECT_NEAR( 0.33722, q[num_qm + 2], 1e-12);
	EXPECT_NEAR(5.73935, a[num_qm + 0], 1e-12);
	EXPECT_NEAR(2.30839, a[num_qm + 1], 1e-12);
	const md_urange_t npe = md_system_component_atom_range(&t.sys, 8);
	ASSERT_EQ(3, npe.end - npe.beg);
	EXPECT_NEAR(-0.834, q[npe.beg], 1e-12);
	EXPECT_NEAR( 0.417, q[npe.beg + 2], 1e-12);
	EXPECT_EQ(0.0, a[npe.beg]);
	EXPECT_EQ(0.0, a[npe.beg + 2]);

	size_t num_polarizable = 0;
	double sum = 0.0;
	for (size_t i = num_qm; i < num_atoms; ++i) {
		sum += q[i];
		num_polarizable += a[i] > 0.0;
	}
	EXPECT_NEAR(0.0, sum, 1e-6);			// Neutral waters
	EXPECT_EQ(39u * 3u, num_polarizable);

	const md_attribute_t* charge = qm_test_attr(&t, STR_LIT("atom/charge"));
	ASSERT_TRUE(charge != NULL);
	EXPECT_TRUE(md_unit_equal(charge->unit, md_unit_elementary_charge()));

	// The file's own per atom columns run over the system's atoms too, with no value for the sites
	double z[153];
	ASSERT_EQ(num_atoms, qm_test_series(z, num_atoms, &t, STR_LIT("atom/nuclear_charges")));
	EXPECT_EQ(8.0, z[0]);
	EXPECT_EQ(1.0, z[2]);
	EXPECT_TRUE(vlx_test_absent(z[num_qm]));
	EXPECT_TRUE(vlx_test_absent(z[num_atoms - 1]));

	// And what an application infers from that once it is loaded: every water an instance of its own,
	// the QM region one more
	ASSERT_TRUE(md_util_system_infer(&t.sys, &t.state, MD_UTIL_INFER_ALL));
	EXPECT_VALID_SYSTEM(&t.sys);
	EXPECT_EQ(1u + 50u, t.sys.instance.count);
	EXPECT_EQ(1u + 50u, t.sys.structure.count);

	qm_test_free(&t);
}

// Written by hand to reach what the reference does not: a site without an element, an anisotropic
// tensor, the (1,1) block on its own, higher multipoles behind a charge, a site with neither, a
// fragment without sites, a name without the tag, fragment numbers out of order, two classical
// subsystems, and Windows line ends.
static const char vlx_test_custom_json[] =
	"{\r\n"
	"  \"quantum_subsystems\": [{\"nuclei\": []}],\r\n"
	"  \"classical_subsystems\": [\r\n"
	"    {\"classical_fragments\": [\r\n"
	"      {\"index\": 7, \"name\": \"NA_npe\", \"atoms\": [\r\n"
	"        {\"index\": 1, \"element\": \"Na\", \"coordinate\": [0.0, 0.0, 20.0], \"multipoles\": {\"elements\": [1.0]}, \"exclusions\": [1]}\r\n"
	"      ]},\r\n"
	"      {\"index\": 5, \"name\": \"HOH_pe\", \"atoms\": [\r\n"
	"        {\"index\": 2, \"element\": \"O\", \"coordinate\": [0.0, 0.0, 10.0], \"multipoles\": {\"elements\": [-0.8]},\r\n"
	"         \"polarizabilities\": {\"elements\": [0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 6.0, 0.0, 6.0], \"order\": [1, 1]}},\r\n"
	"        {\"index\": 3, \"element\": \"H\", \"coordinate\": [1.0, 0.0, 10.0], \"multipoles\": {\"elements\": [0.4, 0.1, 0.2, 0.3]},\r\n"
	"         \"polarizabilities\": {\"elements\": [2.0, 0.0, 0.0, 2.0, 0.0, 2.0]}},\r\n"
	"        {\"index\": 4, \"element\": \"H\", \"coordinate\": [-1.0, 0.0, 10.0], \"multipoles\": {\"elements\": [0.5]},\r\n"
	"         \"polarizabilities\": {\"elements\": [0, 0, 0, 0, 2, 0, 0, 2, 0, 2], \"order\": [1, 1]}},\r\n"
	"        {\"index\": 5, \"element\": \"X\", \"coordinate\": [0.5, 0.0, 10.0], \"multipoles\": {\"elements\": [-0.5]},\r\n"
	"         \"polarizabilities\": {\"elements\": [0.0, 0.0, 0.0, 0.0, 1.0, 0.3, 0.2, 2.0, 0.1, 6.0], \"order\": [1, 1]}}\r\n"
	"      ]},\r\n"
	"      {\"index\": 6, \"name\": \"EMPTY\", \"atoms\": []}\r\n"
	"    ]},\r\n"
	"    {\"classical_fragments\": [\r\n"
	"      {\"index\": 1, \"name\": \"CL\", \"atoms\": [\r\n"
	"        {\"index\": 6, \"element\": \"Cl\", \"coordinate\": [0.0, 5.0, 0.0]}\r\n"
	"      ]}\r\n"
	"    ]}\r\n"
	"  ]\r\n"
	"}\r\n";

UTEST(vlx, pe_environment_in_full) {
	const str_t h5 = STR_LIT(VLX_PE_DIR "unittest_pe_custom.h5");
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, (str_t){ vlx_test_custom_json, sizeof(vlx_test_custom_json) - 1 }));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MEGABYTES(64));
	remove(h5.ptr);
	ASSERT_TRUE(loaded);

	ASSERT_EQ(3u + 6u, t.sys.atom.count);

	const double bohr = 0.5291772109029999;
	EXPECT_NEAR(20.0 * bohr, t.state.xyz[3].z, 1e-5);
	EXPECT_NEAR( 1.0 * bohr, t.state.xyz[5].x, 1e-5);
	EXPECT_NEAR( 5.0 * bohr, t.state.xyz[8].y, 1e-5);

	EXPECT_EQ(11, md_atom_atomic_number(&t.sys.atom, 3));
	EXPECT_EQ(17, md_atom_atomic_number(&t.sys.atom, 8));
	// The expansion point is a virtual site, not an atom of an unknown element
	EXPECT_EQ(0, md_atom_atomic_number(&t.sys.atom, 7));
	EXPECT_EQ(MD_PARTICLE_VIRTUAL_SITE, md_atom_particle_kind(&t.sys.atom, 7));
	EXPECT_EQ(MD_PARTICLE_ATOM, md_atom_particle_kind(&t.sys.atom, 4));

	// The fragment without sites is no component; the numbers are the file's
	ASSERT_EQ(4u, t.sys.component.count);
	EXPECT_VALID_SYSTEM(&t.sys);
	EXPECT_TRUE(str_eq(md_component_name(&t.sys.component, 1), STR_LIT("NA")));
	EXPECT_TRUE(str_eq(md_component_name(&t.sys.component, 2), STR_LIT("HOH")));
	EXPECT_TRUE(str_eq(md_component_name(&t.sys.component, 3), STR_LIT("CL")));
	EXPECT_EQ(7, md_component_seq_id(&t.sys.component, 1));
	EXPECT_EQ(5, md_component_seq_id(&t.sys.component, 2));
	EXPECT_EQ(1, md_component_seq_id(&t.sys.component, 3));
	EXPECT_EQ(MD_COMPONENT_KIND_ION, md_system_component_kind(&t.sys, 1));
	const md_urange_t hoh = md_system_component_atom_range(&t.sys, 2);
	EXPECT_EQ(4, hoh.beg);
	EXPECT_EQ(8, hoh.end);

	double q[9], a[9];
	ASSERT_EQ(9u, qm_test_series(q, 9, &t, STR_LIT("atom/charge")));
	ASSERT_EQ(9u, qm_test_series(a, 9, &t, STR_LIT("atom/polarizability")));
	EXPECT_NEAR( 1.0, q[3], 1e-12);
	EXPECT_NEAR(-0.8, q[4], 1e-12);
	EXPECT_NEAR( 0.4, q[5], 1e-12);	// The charge, the dipole after it plays no part
	EXPECT_NEAR( 0.5, q[6], 1e-12);
	EXPECT_NEAR(-0.5, q[7], 1e-12);
	EXPECT_EQ(0.0, q[8]);				// No multipoles: no charge
	EXPECT_EQ(0.0, a[3]);				// No polarizabilities: embedded non-polarizably
	EXPECT_NEAR(6.0, a[4], 1e-12);
	EXPECT_NEAR(2.0, a[5], 1e-12);		// The (1,1) block on its own
	EXPECT_NEAR(2.0, a[6], 1e-12);
	EXPECT_NEAR(3.0, a[7], 1e-12);		// (1 + 2 + 6) / 3, the off diagonal plays no part
	EXPECT_EQ(0.0, a[8]);

	qm_test_free(&t);
}

// What write_pe_jsonfile could not have written is not the potential the calculation ran with, so
// none of it is used, and a potential in CPPE's format is not read. Neither is a failed load.
UTEST(vlx, pe_environment_left_out_when_it_cannot_be_the_one) {
	static const char* texts[] = {
		// Not JSON, or cut short
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"HOH_pe\", \"atoms\": [",
		"{\"classical_subsystems\": [}",
		// No environment
		"{}",
		"{\"quantum_subsystems\": [{\"nuclei\": []}]}",
		"{\"classical_subsystems\": {}}",
		"{\"classical_subsystems\": [{\"classical_fragments\": []}]}",
		// Sites that are not as written
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"A\", \"atoms\": ["
			"{\"element\": \"O\", \"coordinate\": [0.0, 0.0]}]}]}]}",
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"A\", \"atoms\": ["
			"{\"coordinate\": [0.0, 0.0, 0.0]}]}]}]}",
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"A\", \"atoms\": ["
			"{\"element\": \"O\", \"coordinate\": [0.0, \"0.0\", 0.0]}]}]}]}",
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"A\", \"atoms\": ["
			"{\"element\": \"O\", \"coordinate\": [0.0, 0.0, 0.0], \"multipoles\": {\"elements\": [\"q\"]}}]}]}]}",
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"A\", \"atoms\": ["
			"{\"element\": \"O\", \"coordinate\": [0.0, 0.0, 0.0], \"polarizabilities\": {\"elements\": [1, 2, 3, 4, 5, 6, 7, 8, 9]}}]}]}]}",
		"{\"classical_subsystems\": [{\"classical_fragments\": [{\"index\": 1, \"name\": \"A\", \"atoms\": ["
			"{\"element\": \"O\", \"coordinate\": [0.0, 0.0, 0.0], \"polarizabilities\": {\"elements\": [0, 0, 0, 0, 1, 0, 0, 1, 0, 1], \"order\": [2, 2]}}]}]}]}",
		// CPPE's own format, which a run embedded through CPPE stored
		"@COORDINATES\n3\nAA\nO -5.672 2.390 -4.911\nH -6.099 2.789 -4.153\nH -6.322 2.449 -5.611\n@MULTIPOLES\nORDER 0\n3\n1 -0.67444\n2 0.33722\n3 0.33722\n",
	};
	for (size_t i = 0; i < ARRAY_SIZE(texts); ++i) {
		const str_t h5 = STR_LIT(VLX_PE_DIR "unittest_pe_bad.h5");
		ASSERT_TRUE(vlx_test_write_pe_h5(h5, str_from_cstr(texts[i])));
		vlx_test_t t = {0};
		const bool loaded = vlx_test_load(&t, h5, MEGABYTES(64));
		remove(h5.ptr);

		EXPECT_TRUE_MSG(loaded, texts[i]);
		EXPECT_EQ_MSG(3u, t.sys.atom.count, texts[i]);
		EXPECT_EQ_MSG(0u, t.sys.component.count, texts[i]);
		EXPECT_FALSE_MSG(qm_test_has(&t, STR_LIT("atom/charge")), texts[i]);
		EXPECT_EQ_MSG(3u, qm_test_count(&t, STR_LIT("atom/nuclear_charges")), texts[i]);
		qm_test_free(&t);
	}
}

// A supplemental load leaves the atoms of the system it supplements alone, so it adds no sites
UTEST(vlx, pe_environment_not_added_by_a_supplemental_load) {
	const str_t h5 = STR_LIT(VLX_PE_DIR "unittest_pe_supplement.h5");
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, (str_t){ vlx_test_custom_json, sizeof(vlx_test_custom_json) - 1 }));

	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, STR_LIT(VLX_PE_DIR "h2o.h5"), MEGABYTES(64));
	const bool supplemented = loaded && md_vlx_system_supplement_from_file(&t.sys, h5);
	remove(h5.ptr);

	EXPECT_TRUE(supplemented);
	EXPECT_EQ(3u, t.sys.atom.count);
	EXPECT_EQ(3u, t.state.num_atoms);
	EXPECT_FALSE(qm_test_has(&t, STR_LIT("atom/charge")));

	qm_test_free(&t);
}
