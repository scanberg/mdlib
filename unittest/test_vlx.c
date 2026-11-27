#include "utest.h"

#include "vlx_test_util.h"

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_grid.h>
#include <md_gto_iso.h>
#include <md_qm.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <math.h>
#include <float.h>

// Everything a VeloxChem file carries reaches a consumer as an md_system_t and its ATTRIBUTE TABLE:
// md_vlx.h has three entry points and no reader object to ask questions of. So these tests read what
// a real consumer reads, by the same paths, and a value that cannot be reached this way is a value
// no consumer can reach either - which is the property they exist to hold.

static const double ref_ener_tot = -444.518500783179;

UTEST(vlx, parse) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), MD_MEGABYTES(64)));

	// The atoms land in the system itself, not in the table.
	EXPECT_EQ(26u, t.sys.atom.count);
	EXPECT_EQ(26u, t.state.num_atoms);

	EXPECT_EQ(0.0, qm_test_scalar(&t, MD_STR_LIT("vlx/molecular_charge"), -1.0));
	EXPECT_EQ(1.0, qm_test_scalar(&t, MD_STR_LIT("vlx/spin_multiplicity"), -1.0));
	EXPECT_EQ(41.0, qm_test_scalar(&t, MD_STR_LIT("vlx/electron_count/alpha"), -1.0));
	EXPECT_EQ(41.0, qm_test_scalar(&t, MD_STR_LIT("vlx/electron_count/beta"), -1.0));

	EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("vlx/basis_set")), MD_STR_LIT("DEF2-SVP")));

	// The QM geometry, in Angstrom, as the calculation was run at.
	const md_attribute_t* coord = qm_test_attr(&t, MD_STR_LIT("qm/atom/coordinate"));
	ASSERT_TRUE(coord != NULL);
	ASSERT_EQ(coord->format.components, 3u);
	ASSERT_EQ(md_attribute_value_count(&coord->format), 26u);

	double xyz[3 * 26] = {0};
	ASSERT_EQ(md_attribute_extract_f64(xyz, MD_COUNTOF(xyz), &t.sys.attributes, coord, md_attribute_slice_all(), md_unit_none()), MD_COUNTOF(xyz));
	EXPECT_NEAR(-3.259400000000, xyz[0], 1.0e-5);
	EXPECT_NEAR( 0.145200000000, xyz[1], 1.0e-5);
	EXPECT_NEAR(-0.048400000000, xyz[2], 1.0e-5);

	// The system's own state is that geometry too, narrowed to float.
	EXPECT_NEAR(xyz[0], (double)t.state.xyz[0].x, 1.0e-5);
	EXPECT_NEAR(xyz[1], (double)t.state.xyz[0].y, 1.0e-5);
	EXPECT_NEAR(xyz[2], (double)t.state.xyz[0].z, 1.0e-5);

	const size_t num_iter = qm_test_count(&t, MD_STR_LIT("vlx/scf/history/energy"));
	ASSERT_TRUE(num_iter > 0);

	double* energy = (double*)md_alloc(t.alloc, sizeof(double) * num_iter);
	ASSERT_EQ(qm_test_series(energy, num_iter, &t, MD_STR_LIT("vlx/scf/history/energy")), num_iter);
	EXPECT_NEAR(ref_ener_tot, energy[num_iter - 1], 1.0e-5);

	qm_test_free(&t);
}

// The facts about a calculation that are TEXT and nothing else, plus the counts that used to be
// reachable only by asking the reader. If one of these stops being published, nothing downstream can
// tell a restricted run from an unrestricted one except by guessing from which arrays are shared.
UTEST(vlx, run_description_is_published) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MD_MEGABYTES(16)));

	EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("vlx/scf/type")), MD_STR_LIT("restricted")));
	EXPECT_TRUE(qm_test_scalar(&t, MD_STR_LIT("vlx/spin_multiplicity"), 0.0) == 1.0);
	EXPECT_TRUE(qm_test_scalar(&t, MD_STR_LIT("vlx/electron_count/alpha"), 0.0) > 0.0);
	EXPECT_TRUE(qm_test_scalar(&t, MD_STR_LIT("vlx/nuclear_repulsion_energy"), 0.0) > 0.0);

	// A response calculation, so the type is there and names which one.
	EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("vlx/rsp/type")), MD_STR_LIT("linear")));

	// No geometry optimisation in this file, so nothing under vlx/opt. An absent path is how a
	// consumer learns a block is missing - not a zero it would have to interpret.
	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("vlx/opt/type")));
	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("vlx/opt/state_index")));
	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("vlx/opt/irc_ts_index")));

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

	for (size_t c = 0; c < MD_COUNTOF(cases); ++c) {
		vlx_test_t t = {0};
		ASSERT_TRUE(vlx_test_load(&t, md_str_from_cstr(cases[c].path), MD_MEGABYTES(64)));

		// {S,Lmax}: one row per excited state, padded to the widest row with zeros.
		const md_attribute_t* a = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/lambda"));
		ASSERT_TRUE(a != NULL);
		ASSERT_EQ(a->format.rank, 2u);

		const size_t num_lambdas = a->format.shape[1];
		ASSERT_TRUE(num_lambdas >= cases[c].count);

		double* lambdas = (double*)md_alloc(t.alloc, sizeof(double) * num_lambdas);
		ASSERT_EQ(qm_test_row(lambdas, num_lambdas, &t, MD_STR_LIT("vlx/rsp/nto/lambda"), 0), num_lambdas);

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
	md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5");
	ASSERT_TRUE(vlx_test_load(&t, path, MD_MEGABYTES(32)));

	md_allocator_i* arena = t.alloc;

	const size_t num_atoms = t.sys.atom.count;
	float* atom_xyz = (float*)md_alloc(arena, sizeof(float) * 3 * num_atoms);
	ASSERT_EQ(vlx_test_atom_xyz_bohr(atom_xyz, 3 * num_atoms, &t), num_atoms);

	// The volume dimensions which we aim to sample molecular orbital over
	const int vol_dim = 80;

	// The molecular orbital index we aim to sample - the HOMO, one below the first empty orbital.
	const size_t lumo_idx = vlx_test_lumo_idx(&t, MD_STR_LIT("orbital/alpha/occupation"));
	ASSERT_TRUE(lumo_idx > 0);
	const size_t mo_idx = lumo_idx - 1;

	// Extract the GTO basis from the attributes the loader published
	md_gto_basis_t basis = {0};
	ASSERT_TRUE(qm_test_basis(&basis, &t));

	const size_t num_aos = md_gto_basis_num_ao(&basis);
	double* mo_coeffs = (double*)md_alloc(arena, sizeof(double) * num_aos);
	ASSERT_EQ(qm_test_row(mo_coeffs, num_aos, &t, MD_STR_LIT("orbital/alpha/coefficient"), mo_idx), num_aos);

	// The orbital as a field over the basis (md_gto_iso.h): what is evaluated at any point
	md_gto_iso_field_t* field = md_gto_iso_field_create(&basis, atom_xyz, sizeof(float) * 3, mo_coeffs, arena);
	ASSERT_TRUE(field != NULL);

	// Calculate bounding box (AABB) of the atoms
	md_float3_t min_aabb = md_float3_set1( FLT_MAX);
	md_float3_t max_aabb = md_float3_set1(-FLT_MAX);

	for (size_t i = 0; i < num_atoms; ++i) {
		md_float3_t coord = md_float3_set(atom_xyz[3 * i + 0], atom_xyz[3 * i + 1], atom_xyz[3 * i + 2]);
		min_aabb = md_float3_min(min_aabb, coord);
		max_aabb = md_float3_max(max_aabb, coord);
	}

	// Add some padding
	const float pad = 6.0f;
	min_aabb = md_float3_sub1(min_aabb, pad);
	max_aabb = md_float3_add1(max_aabb, pad);

	md_float3_t ext_aabb = md_float3_sub(max_aabb, min_aabb);
	md_float3_t step_size = md_float3_div1(ext_aabb, (float)vol_dim);

	// Shift origin by half a voxel such that the samples are constructed from the center of each voxel
	md_float3_t origin = md_float3_add(min_aabb, md_float3_mul1(step_size, 0.5f));

	// Setup the grid structure that control how we aim to sample over space
	md_grid_t grid = (md_grid_t) {
		.orientation = md_float3x3_ident(),
		.origin = origin,
		.spacing = step_size,
		.dim = {vol_dim, vol_dim, vol_dim},
	};

	// Evaluate the orbital at the points of the grid
	const size_t num_points = md_grid_num_points(&grid);
	md_float3_t* points = (md_float3_t*)md_alloc(arena, sizeof(md_float3_t) * num_points);
	double* vol_data = (double*)md_alloc(arena, sizeof(double) * num_points);
	md_grid_extract_points((float*)points, &grid);
	md_gto_iso_field_evaluate(field, (const float*)points, num_points, sizeof(md_float3_t), vol_data, NULL);
	md_gto_iso_field_destroy(field);

	// An occupied orbital is not the zero function; this is what catches a coefficient row that
	// came back empty rather than wrong.
	double max_value = 0.0;
	for (size_t i = 0; i < num_points; ++i) {
		max_value = MD_MAX(max_value, fabs(vol_data[i]));
	}
	EXPECT_TRUE(max_value > 1.0e-6);

	qm_test_free(&t);
}

UTEST(vlx, mol) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), MD_MEGABYTES(64)));
	qm_test_free(&t);
}

UTEST(vlx, scf_results) {
	vlx_test_t t = {0};
	EXPECT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/tq.scf.results.h5"), MD_MEGABYTES(64)));
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
	EXPECT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MD_MEGABYTES(64)));

	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("vlx/xps/ionization_energy")));
	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("vlx/xps/element")));

	md_attribute_id_t ids[8];
	EXPECT_EQ(md_attributes_query(ids, MD_COUNTOF(ids), &t.sys.attributes, MD_STR_LIT("vlx/xps")), 0u);

	qm_test_free(&t);
}

// XPS entries, published as one attribute per field of the record over a shared {C} index space.
// The per element grouping is NOT published: entries are laid out as contiguous runs of equal
// element, so a consumer scans vlx/xps/element for the run it wants. That derivation is what this
// pins, because it is the thing that breaks silently if the sort ever changes.
UTEST(vlx, acro_xps) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/acro-xps.h5"), MD_MEGABYTES(64)));

	const size_t count = qm_test_count(&t, MD_STR_LIT("vlx/xps/ionization_energy"));
	ASSERT_TRUE(count > 0);

	double* energy       = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* element      = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* contribution = (double*)md_alloc(t.alloc, sizeof(double) * count);
	double* atom_index   = (double*)md_alloc(t.alloc, sizeof(double) * count);

	ASSERT_EQ(qm_test_series(energy,       count, &t, MD_STR_LIT("vlx/xps/ionization_energy")), count);
	ASSERT_EQ(qm_test_series(element,      count, &t, MD_STR_LIT("vlx/xps/element")),           count);
	ASSERT_EQ(qm_test_series(contribution, count, &t, MD_STR_LIT("vlx/xps/contribution")),      count);
	ASSERT_EQ(qm_test_series(atom_index,   count, &t, MD_STR_LIT("vlx/xps/atom_index")),        count);

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

// One dataset of an .h5 file straight through HDF5, as double (heap; dims[0..rank) and *out_rank set): what a
// test compares the reader against, read without it
static double* vlx_test_h5_dataset(size_t dims[2], int* out_rank, const char* path, const char* name) {
	double* out = NULL;
	hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file < 0) return NULL;
	hid_t dset  = H5Lexists(file, name, H5P_DEFAULT) > 0 ? H5Dopen2(file, name, H5P_DEFAULT) : -1;
	hid_t space = dset >= 0 ? H5Dget_space(dset) : -1;
	hsize_t d[2] = {1, 1};
	const int rank = space >= 0 ? H5Sget_simple_extent_dims(space, d, NULL) : -1;
	if (rank == 1 || rank == 2) {
		out = (double*)malloc(sizeof(double) * d[0] * d[1]);
		if (out && H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out) < 0) {
			free(out);
			out = NULL;
		}
		dims[0] = (size_t)d[0];
		dims[1] = rank == 2 ? (size_t)d[1] : 1;
		*out_rank = rank;
	}
	if (space >= 0) H5Sclose(space);
	if (dset >= 0) H5Dclose(dset);
	H5Fclose(file);
	return out;
}

// VeloxChem's own SCF density, D_alpha, which no reader publishes: a density is made from the orbitals and their
// occupations (md_qm_density_factors). This is the check that nothing is lost by that, and the one place the file's
// matrix is read at all: in the file's own AO order and spherical basis, D_alpha is sum_k n_k C_k C_k^T over the
// orbitals and occupations the file stores beside it (C_alpha [AO][MO], as VeloxChem writes it). That those orbitals
// reach md_gto's Cartesian basis intact is what the orthonormality and electron count tests check.
UTEST(vlx, file_density_is_the_orbitals) {
	static const char* files[] = { MD_UNITTEST_DATA_DIR "/vlx/h2o.h5", MD_UNITTEST_DATA_DIR "/vlx/amide.h5" };
	for (int fi = 0; fi < 2; ++fi) {
		size_t dd[2] = {0}, cd[2] = {0}, od[2] = {0};
		int dr = 0, cr = 0, orank = 0;
		double* D   = vlx_test_h5_dataset(dd, &dr, files[fi], "scf/D_alpha");
		double* C   = vlx_test_h5_dataset(cd, &cr, files[fi], "scf/C_alpha");
		double* occ = vlx_test_h5_dataset(od, &orank, files[fi], "scf/occ_alpha");
		ASSERT_TRUE(D && C && occ);
		ASSERT_EQ(2, dr);
		ASSERT_EQ(2, cr);
		const size_t na = dd[0], nm = cd[1];
		ASSERT_EQ(na, dd[1]);
		ASSERT_EQ(na, cd[0]);
		ASSERT_EQ(nm, od[0]);

		double max_diff = 0.0, magnitude = 0.0;
		for (size_t i = 0; i < na; ++i) {
			for (size_t j = 0; j < na; ++j) {
				double v = 0.0;
				for (size_t k = 0; k < nm; ++k) v += occ[k] * C[i * nm + k] * C[j * nm + k];
				max_diff  = MD_MAX(max_diff, fabs(D[i * na + j] - v));
				magnitude = MD_MAX(magnitude, fabs(D[i * na + j]));
			}
		}
		printf("  %s: %zu AOs, %zu orbitals; D_alpha against sum_k n_k C_k C_k^T: max |diff| %.1e (max |D| %.3g)\n",
			files[fi], na, nm, max_diff, magnitude);
		EXPECT_LT(max_diff, 1.0e-12 * MD_MAX(magnitude, 1.0));
		// guards against the all zeros case, which the difference would also satisfy
		EXPECT_TRUE(magnitude > 1.0e-8);
		free(D);
		free(C);
		free(occ);
	}
}

// The SCF densities are names, not attributes: orbital/{alpha,beta,total,difference}/density, which
// md_qm_density_factors makes from the orbitals and occupations. A restricted calculation gives the invariants for
// free and needs no reference values: beta is alpha (beta's orbitals a second name for alpha's), the total is exactly
// twice alpha, and the difference is zero - no factors at all.
UTEST(vlx, scf_density_names_h2o) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MD_MEGABYTES(16)));
	EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("vlx/scf/type")), MD_STR_LIT("restricted")));

	static const md_str_t paths[4] = {
		MD_STR_INIT("orbital/alpha/density"),
		MD_STR_INIT("orbital/beta/density"),
		MD_STR_INIT("orbital/total/density"),
		MD_STR_INIT("orbital/difference/density"),
	};
	for (int i = 0; i < 4; ++i) {
		EXPECT_TRUE(md_qm_density_factors_available(&t.sys, paths[i]));
		EXPECT_FALSE(qm_test_has(&t, paths[i]));	// no AO matrix is published
	}

	// alpha: its occupied orbitals, weighted by their occupations, which add up to the alpha electrons
	double *rows[4] = {0}, *w[4] = {0};
	size_t count[4] = {0}, na[4] = {0};
	for (int i = 0; i < 4; ++i) {
		count[i] = md_qm_density_factors(&rows[i], &w[i], &na[i], &t.sys, paths[i], NULL, t.alloc);
	}
	ASSERT_TRUE(count[0] > 0);
	double occ_sum = 0.0;
	for (size_t k = 0; k < count[0]; ++k) occ_sum += w[0][k];
	EXPECT_NEAR(qm_test_scalar(&t, MD_STR_LIT("vlx/electron_count/alpha"), -1.0), occ_sum, 1.0e-9);

	// beta is alpha, and the total the same orbitals at twice the weight
	ASSERT_EQ(count[0], count[1]);
	ASSERT_EQ(count[0], count[2]);
	for (size_t k = 0; k < count[0]; ++k) {
		EXPECT_EQ(w[0][k], w[1][k]);
		EXPECT_EQ(2.0 * w[0][k], w[2][k]);
		for (size_t a = 0; a < na[0]; ++a) {
			EXPECT_EQ(rows[0][k * na[0] + a], rows[1][k * na[0] + a]);
			EXPECT_EQ(rows[0][k * na[0] + a], rows[2][k * na[0] + a]);
		}
	}
	// the difference: nothing
	EXPECT_EQ(0u, count[3]);

	qm_test_free(&t);
}

// An excited state's attachment, detachment and difference densities, published as their factors
// (vlx/rsp/transition_density/<name>/factor {S,R,A}, virtual, and /weight {S,R}, resident), never as matrices: what
// md_qm_density_factors hands out sums to the densities formed here from the amplitudes the long way (C_v^T t^T t C_v,
// C_o^T t t^T C_o, a dense reference only this test makes). A slice is load bearing: one state's factors are computed
// for one state. Both paths are exercised below, and the physics each must keep: both densities hold |t|^2 electrons
// (tr(D S), the orbitals orthonormal in S), and the difference none.
static double* td_reference(const vlx_test_t* t, uint32_t state, int which, size_t* out_num_ao) {
	const size_t nocc = (size_t)(qm_test_scalar(t, MD_STR_LIT("vlx/rsp/num_core"), 0.0) > 0.0 ? qm_test_scalar(t, MD_STR_LIT("vlx/rsp/num_core"), 0.0)
	                                                                                          : qm_test_scalar(t, MD_STR_LIT("vlx/rsp/num_valence"), 0.0));
	const size_t nvir = (size_t)qm_test_scalar(t, MD_STR_LIT("vlx/rsp/num_virtual"), 0.0);
	const md_attribute_t* sol = qm_test_attr(t, MD_STR_LIT("vlx/rsp/solution_matrix"));
	const md_attribute_t* ca = qm_test_attr(t, MD_STR_LIT("orbital/alpha/coefficient"));
	if (!sol || !ca || !nocc || !nvir) return NULL;
	const size_t vec_len = sol->format.shape[1], nmo = ca->format.shape[0], A = ca->format.shape[1];
	double* row = (double*)malloc(sizeof(double) * vec_len);
	double* C = (double*)malloc(sizeof(double) * nmo * A);
	double* amp = (double*)malloc(sizeof(double) * nocc * nvir);
	double* P = (double*)calloc(nmo * nmo, sizeof(double));     // the density in the MO basis
	double* D = (double*)calloc(A * A, sizeof(double));
	md_attribute_extract_f64(row, vec_len, &t->sys.attributes, sol, md_attribute_slice_1(state), md_unit_none());
	md_attribute_extract_f64(C, nmo * A, &t->sys.attributes, ca, md_attribute_slice_all(), md_unit_none());
	for (size_t k = 0; k < nocc * nvir; ++k) amp[k] = vec_len == 2 * nocc * nvir ? row[k] - row[nocc * nvir + k] : row[k];
	const double s_att = which == 1 ? 0.0 : 1.0, s_det = which == 0 ? 0.0 : which == 1 ? 1.0 : -1.0;     // difference: attachment - detachment
	for (size_t a = 0; a < nvir; ++a) for (size_t b = 0; b < nvir; ++b) {
		double v = 0.0;
		for (size_t i = 0; i < nocc; ++i) v += amp[i * nvir + a] * amp[i * nvir + b];
		P[(nocc + a) * nmo + nocc + b] += s_att * v;
	}
	for (size_t i = 0; i < nocc; ++i) for (size_t j = 0; j < nocc; ++j) {
		double v = 0.0;
		for (size_t a = 0; a < nvir; ++a) v += amp[i * nvir + a] * amp[j * nvir + a];
		P[i * nmo + j] += s_det * v;
	}
	for (size_t p = 0; p < nmo; ++p) for (size_t q = 0; q < nmo; ++q) {
		const double v = P[p * nmo + q];
		if (v == 0.0) continue;
		for (size_t mu = 0; mu < A; ++mu) for (size_t nu = 0; nu < A; ++nu) D[mu * A + nu] += v * C[p * A + mu] * C[q * A + nu];
	}
	free(row); free(C); free(amp); free(P);
	*out_num_ao = A;
	return D;
}

static void check_transition_density_attributes(int* utest_result, md_str_t file) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, file, MD_MEGABYTES(64)));

	const size_t num_states = qm_test_count(&t, MD_STR_LIT("vlx/rsp/oscillator_strength"));
	ASSERT_TRUE(num_states > 0);
	const md_attribute_t* coeff = qm_test_attr(&t, MD_STR_LIT("orbital/alpha/coefficient"));
	ASSERT_TRUE(coeff != NULL);
	const size_t A = coeff->format.shape[1];
	size_t sn = 0;
	double* S = qm_test_matrix(&t, MD_STR_LIT("basis/overlap"), &sn);
	ASSERT_TRUE(S != NULL);
	ASSERT_EQ(A, sn);

	static const char* names[3] = { "attachment", "detachment", "difference" };
	for (int k = 0; k < 3; ++k) {
		char fp[96], wp[96], dp[96];
		snprintf(dp, sizeof(dp), "vlx/rsp/transition_density/%s", names[k]);
		snprintf(fp, sizeof(fp), "%s/factor", dp);
		snprintf(wp, sizeof(wp), "%s/weight", dp);
		const md_attribute_t* fa = qm_test_attr(&t, md_str_from_cstr(fp));
		const md_attribute_t* wa = qm_test_attr(&t, md_str_from_cstr(wp));
		ASSERT_TRUE(fa != NULL && wa != NULL);
		EXPECT_TRUE(qm_test_attr(&t, md_str_from_cstr(dp)) == NULL);     // no matrix
		EXPECT_TRUE(md_qm_density_factors_available(&t.sys, md_str_from_cstr(dp)));
		// {S,R,A} virtual, its weights {S,R} resident; the AO axis the coefficients'
		EXPECT_TRUE(md_attribute_is_virtual(fa));
		EXPECT_FALSE(md_attribute_is_virtual(wa));
		ASSERT_EQ(fa->format.rank, 3u);
		ASSERT_EQ(wa->format.rank, 2u);
		EXPECT_EQ(fa->format.shape[0], (uint32_t)num_states);
		EXPECT_EQ(fa->format.shape[1], wa->format.shape[1]);
		EXPECT_EQ(fa->format.shape[2], (uint32_t)A);
		const size_t R = fa->format.shape[1];

		for (uint32_t st = 0; st < (uint32_t)num_states && st < 3; ++st) {
			const md_attribute_slice_t slice = md_attribute_slice_1(st);
			double *rows = NULL, *w = NULL;
			size_t na = 0;
			const size_t n = md_qm_density_factors(&rows, &w, &na, &t.sys, md_str_from_cstr(dp), &slice, t.alloc);
			ASSERT_GT(n, 0u);
			ASSERT_LE(n, R);
			ASSERT_EQ(A, na);
			size_t rn = 0;
			double* ref = td_reference(&t, st, k, &rn);
			ASSERT_TRUE(ref != NULL);
			ASSERT_EQ(A, rn);
			double emax = 0.0, rmax = 0.0, nel = 0.0, nel_ref = 0.0;
			for (size_t mu = 0; mu < A; ++mu) for (size_t nu = 0; nu < A; ++nu) {
				double v = 0.0;
				for (size_t r = 0; r < n; ++r) v += w[r] * rows[r * A + mu] * rows[r * A + nu];
				emax = fmax(emax, fabs(v - ref[mu * A + nu]));
				rmax = fmax(rmax, fabs(ref[mu * A + nu]));
				nel += v * S[mu * A + nu];
				nel_ref += ref[mu * A + nu] * S[mu * A + nu];
			}
			printf("  " MD_STR_FMT " state %u %-10s %zu of %zu factors: %.1e of max %.3g; tr(D S) %.3e (reference %.3e)\n", MD_STR_ARG(file), st, names[k],
				   n, R, emax / rmax, rmax, nel, nel_ref);
			EXPECT_TRUE(rmax > 1.0e-8);
			EXPECT_LT(emax, 1.0e-12 * rmax);
			EXPECT_NEAR(nel_ref, nel, 1.0e-10);
			// attachment and detachment hold as many electrons, |t|^2: to within the orbitals' orthonormality in the
			// overlap integrated from the basis (overlap_describes_the_same_basis_as_the_coefficients: 1e-5)
			if (k == 2) { EXPECT_NEAR(0.0, nel, 1.0e-5); }
			free(ref);
			md_free(t.alloc, rows, sizeof(double) * n * A);
			md_free(t.alloc, w, sizeof(double) * n);
		}

		// the whole attribute, every state at once: state 0's plane as its slice gave it
		const size_t plane = R * A, total = num_states * plane;
		double* one = (double*)md_alloc(t.alloc, sizeof(double) * plane);
		double* all = (double*)md_alloc(t.alloc, sizeof(double) * total);
		ASSERT_EQ(md_attribute_extract_f64(one, plane, &t.sys.attributes, fa, md_attribute_slice_1(0), md_unit_none()), plane);
		ASSERT_EQ(md_attribute_extract_f64(all, total, &t.sys.attributes, fa, md_attribute_slice_all(), md_unit_none()), total);
		EXPECT_EQ(0, memcmp(one, all, sizeof(double) * plane));
		// a state past the end selects nothing rather than reading past the array
		const md_attribute_slice_t past_end = md_attribute_slice_1((uint32_t)num_states);
		EXPECT_EQ(md_attribute_slice_count(fa, past_end), 0u);
		EXPECT_EQ(md_attribute_extract_f64(one, plane, &t.sys.attributes, fa, past_end, md_unit_none()), 0u);
	}
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
	for (size_t i = 0; i < MD_COUNTOF(files); ++i) {
		vlx_test_t t = {0};
		ASSERT_TRUE(vlx_test_load(&t, md_str_from_cstr(files[i]), MD_MEGABYTES(256)));

		// The molecular orbitals are orthonormal against it, or the two disagree about the basis.
		EXPECT_LT(qm_test_orthonormality(&t, MD_STR_LIT("orbital/alpha/coefficient")), 1.0e-5);

		// ...and tr(D S) is the electron count the file itself states. This is the assertion the
		// occupation convention shows up in, which orthonormality is completely blind to.
		const double alpha = qm_test_scalar(&t, MD_STR_LIT("vlx/electron_count/alpha"), -1.0);
		const double beta  = qm_test_scalar(&t, MD_STR_LIT("vlx/electron_count/beta"),  -1.0);
		ASSERT_TRUE(alpha > 0.0);
		ASSERT_TRUE(beta  > 0.0);

		EXPECT_NEAR(alpha,        qm_test_electron_count(&t, MD_STR_LIT("orbital/alpha/density")), 1.0e-4);
		EXPECT_NEAR(alpha + beta, qm_test_electron_count(&t, MD_STR_LIT("orbital/total/density")), 1.0e-4);
		EXPECT_NEAR(0.0,          qm_test_electron_count(&t, MD_STR_LIT("orbital/difference/density")), 1.0e-9);

		qm_test_free(&t);
	}
}

UTEST(vlx, transition_density_attributes_h2o) {
	check_transition_density_attributes(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"));
}

UTEST(vlx, transition_density_attributes_amide) {
	check_transition_density_attributes(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.h5"));
}

// The NTO coefficient vectors are resident and the transition densities are computed, but both are
// built from the same solution vectors - so the leading particle NTO of a state has to be a vector
// in the same AO space, of the same length, as everything else the table publishes over AOs. This is
// what catches the two coming apart on the lambda padding axis.
UTEST(vlx, nto_coefficients_share_the_ao_axis) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MD_MEGABYTES(32)));

	const md_attribute_t* lambda   = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/lambda"));
	const md_attribute_t* particle = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/particle/coefficient"));
	const md_attribute_t* hole     = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/hole/coefficient"));
	const md_attribute_t* coeff    = qm_test_attr(&t, MD_STR_LIT("orbital/alpha/coefficient"));
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
	ASSERT_EQ(md_attribute_extract_f64(vec, num_ao, &t.sys.attributes, particle, leading, md_unit_none()), num_ao);

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
// No .h5 of an embedding run is checked in, so these make one: h2o.h5 with a potential named in its
// SCF settings exactly as VeloxChem writes it (scf/potfile, a scalar UTF-8 string holding the path the
// run was given). It is written beside h2o.h5, where its basis set is, and removed again as soon as
// it has been read - before any assertion that could end the test early.
// ---------------------------------------------------------------------------

#include <stdio.h>	// remove
#include <core/md_file.h>
#include <md_infer.h>
#include <md_filter.h>
#include <md_script.h>
#include <core/md_bitset.h>
#include <string.h>
#include "system_invariants.h"

#define VLX_PE_DIR MD_UNITTEST_DATA_DIR "/vlx/"

static bool vlx_test_write_file(md_str_t path, const void* data, size_t size) {
	md_file_t out = {0};
	if (!md_file_open(&out, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
	const bool ok = md_file_write(out, data, size) == size;
	md_file_close(&out);
	return ok;
}

// A copy of h2o.h5 whose SCF settings name 'potfile'
static bool vlx_test_write_pe_h5(md_str_t dst, const char* potfile) {
	md_file_t in = {0};
	if (!md_file_open(&in, MD_STR_LIT(VLX_PE_DIR "h2o.h5"), MD_FILE_READ)) return false;
	const size_t size = (size_t)md_file_size(in);
	void* bytes = md_alloc(md_get_heap_allocator(), size);
	const bool read = md_file_read(in, bytes, size) == size;
	md_file_close(&in);
	const bool written = read && vlx_test_write_file(dst, bytes, size);
	md_free(md_get_heap_allocator(), bytes, size);
	if (!written) return false;

	char path[1024];
	md_str_copy_to_char_buf(path, sizeof(path), dst);
	hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0) return false;
	hid_t scf   = H5Gopen(file, "scf", H5P_DEFAULT);
	hid_t type  = H5Tcopy(H5T_C_S1);
	H5Tset_size(type, H5T_VARIABLE);
	H5Tset_cset(type, H5T_CSET_UTF8);
	hid_t space = H5Screate(H5S_SCALAR);
	hid_t dset  = scf >= 0 ? H5Dcreate2(scf, "potfile", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
	const bool ok = dset >= 0 && H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &potfile) >= 0;
	if (dset >= 0) H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(type);
	if (scf >= 0) H5Gclose(scf);
	H5Fclose(file);
	return ok;
}

// Read on the bits: the library is built with fast math, where NAN does not compare as itself
static bool vlx_test_absent(double v) {
	uint64_t u;
	memcpy(&u, &v, sizeof(u));
	return (u & 0x7fffffffffffffffull) > 0x7ff0000000000000ull;
}

// The reference potential beside h2o.h5's own folder, named relative to it the way a run directory
// copied as a whole would name it: 50 waters, the 39 polarizable ones (117 sites) listed first and the
// 11 non-polarizable ones after them, though numbered among them (8, 12, 16, ..., 49). As atoms they
// follow their numbers: water n is the sites num_qm + 3 (n - 1) and on.
UTEST(vlx, pe_environment_becomes_atoms_of_the_system) {
	const md_str_t h5 = MD_STR_LIT(VLX_PE_DIR "unittest_pe_relative.h5");
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, "../pot/water_pe_npe.pot"));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MD_MEGABYTES(64));
	remove(h5.ptr);
	ASSERT_TRUE(loaded);

	const size_t num_qm = 3, num_mm = 150, num_atoms = num_qm + num_mm;
	ASSERT_EQ(num_atoms, t.sys.atom.count);
	ASSERT_EQ(num_atoms, t.state.num_atoms);

	// The QM atom domain is untouched, and its atoms are still the system's first
	EXPECT_EQ(num_qm, qm_test_count(&t, MD_STR_LIT("qm/atom/atomic_number")));
	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("qm/atom/system_index")));
	EXPECT_EQ(8, md_atom_atomic_number(&t.sys.atom, 0));
	EXPECT_EQ(1, md_atom_atomic_number(&t.sys.atom, 1));

	// The sites follow, with their elements, names and coordinates (already Angstrom in this file)
	EXPECT_EQ(8, md_atom_atomic_number(&t.sys.atom, num_qm));
	EXPECT_TRUE(md_str_eq(md_atom_name(&t.sys.atom, num_qm),     MD_STR_LIT("OW")));
	EXPECT_TRUE(md_str_eq(md_atom_name(&t.sys.atom, num_qm + 2), MD_STR_LIT("H2")));
	EXPECT_NEAR(-5.672, t.state.xyz[num_qm].x, 1e-5);
	EXPECT_NEAR( 2.390, t.state.xyz[num_qm].y, 1e-5);
	EXPECT_NEAR(-4.911, t.state.xyz[num_qm].z, 1e-5);
	EXPECT_NEAR(-11.932, t.state.xyz[num_qm + 21].x, 1e-5);	// Water 8, the file's first non-polarizable one
	EXPECT_NEAR(-5.908, t.state.xyz[num_qm + 146].x, 1e-5);		// Water 49, the file's last site
	EXPECT_NEAR(-8.589, t.state.xyz[num_atoms - 1].x, 1e-5);	// Water 50, polarizable

	// One component per fragment, named for its residue and numbered by its fragment number, after
	// one for the QM region - components cover every atom or none.
	ASSERT_EQ(1u + 50u, t.sys.residue.count);
	EXPECT_VALID_SYSTEM(&t.sys);
	// What the fragments are is inferred from their atoms and bonds, as for a structure file
	ASSERT_TRUE(md_system_infer(&t.sys, &t.state, NULL));
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, 0), MD_STR_LIT("QM")));
	EXPECT_EQ(0, md_system_residue_find_by_atom_idx(&t.sys, 0));
	EXPECT_EQ(1, md_system_residue_find_by_atom_idx(&t.sys, num_qm));
	md_urange_t first = md_system_residue_atom_range(&t.sys, 1);
	EXPECT_EQ(num_qm, first.beg);
	EXPECT_EQ(num_qm + 3, first.end);
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, 1), MD_STR_LIT("HOH")));
	EXPECT_EQ(1, md_system_residue_seq_id(&t.sys, 1));
	const md_residue_idx_t npe = md_system_residue_find_by_atom_idx(&t.sys, num_qm + 21);
	ASSERT_TRUE(npe >= 0);
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, npe), MD_STR_LIT("HOH")));
	EXPECT_EQ(8, md_system_residue_seq_id(&t.sys, npe));
	for (size_t c = 1; c < t.sys.residue.count; ++c) {
		EXPECT_EQ((int)c, md_system_residue_seq_id(&t.sys, c));
		EXPECT_EQ(MD_RESIDUE_KIND_WATER, md_system_residue_kind(&t.sys, c));
	}

	// The parameters, per atom over the system: the type's rows by position in the fragment, and
	// nothing at all for the QM atoms
	double q[153], a[153];
	ASSERT_EQ(num_atoms, qm_test_series(q, num_atoms, &t, MD_STR_LIT("atom/charge")));
	ASSERT_EQ(num_atoms, qm_test_series(a, num_atoms, &t, MD_STR_LIT("atom/polarizability")));
	for (size_t i = 0; i < num_qm; ++i) {
		EXPECT_TRUE(vlx_test_absent(q[i]));
		EXPECT_TRUE(vlx_test_absent(a[i]));
	}
	EXPECT_NEAR(-0.67444, q[num_qm + 0], 1e-8);
	EXPECT_NEAR( 0.33722, q[num_qm + 1], 1e-8);
	EXPECT_NEAR( 0.33722, q[num_qm + 2], 1e-8);
	EXPECT_NEAR(-0.83400, q[num_qm + 21], 1e-8);
	EXPECT_NEAR( 0.41700, q[num_qm + 146], 1e-8);
	EXPECT_NEAR( 0.33722, q[num_atoms - 1], 1e-8);
	EXPECT_NEAR(5.73935, a[num_qm + 0], 1e-8);
	EXPECT_NEAR(2.30839, a[num_qm + 1], 1e-8);
	EXPECT_EQ(0.0, a[num_qm + 21]);	// HOH_npe has no polarizabilities: embedded non-polarizably
	EXPECT_EQ(0.0, a[num_qm + 146]);
	EXPECT_NEAR(2.30839, a[num_atoms - 1], 1e-8);

	double sum = 0.0;
	for (size_t i = num_qm; i < num_atoms; ++i) sum += q[i];
	EXPECT_NEAR(0.0, sum, 1e-6);	// Neutral waters

	const md_attribute_t* charge = qm_test_attr(&t, MD_STR_LIT("atom/charge"));
	ASSERT_TRUE(charge != NULL);
	EXPECT_TRUE(md_unit_equal(charge->unit, md_unit_elementary_charge()));

	// The file's own per atom columns run over the system's atoms too, with no value for the sites
	double z[153];
	ASSERT_EQ(num_atoms, qm_test_series(z, num_atoms, &t, MD_STR_LIT("atom/nuclear_charges")));
	EXPECT_EQ(8.0, z[0]);
	EXPECT_EQ(1.0, z[2]);
	EXPECT_TRUE(vlx_test_absent(z[num_qm]));
	EXPECT_TRUE(vlx_test_absent(z[num_atoms - 1]));

	// And what an application infers from that once it is loaded: every water a part of its own, the QM region one
	// more. The file gives no units: no instances.
	ASSERT_TRUE(md_system_infer(&t.sys, &t.state, NULL));
	EXPECT_VALID_SYSTEM(&t.sys);
	EXPECT_EQ(0u, t.sys.instance.count);
	EXPECT_EQ(1u + 50u, t.sys.part.count);
	EXPECT_EQ(1u + 50u, t.sys.molecule.count);

	qm_test_free(&t);
}

// A potential written by hand to reach what the reference file does not: atomic units, a site
// without an element, an anisotropic tensor, and a fragment type the file lists out of order. Named
// by an absolute path from the machine the run was on, with the file itself copied beside the .h5,
// which is VeloxChem's own fallback.
static const char vlx_test_custom_pot[] =
	"@environment\n"
	"units: au\n"
	"xyz:\n"
	"Na   0.0  0.0 20.0  NA_npe 4 NA\n"
	"O    0.0  0.0 10.0  HOH_pe 5 OW\n"
	"H    1.0  0.0 10.0  HOH_pe 5 HW1\n"
	"H   -1.0  0.0 10.0  HOH_pe 5 HW2\n"
	"X    0.5  0.0 10.0  HOH_pe 5 X1\n"
	"@end\n"
	"@charges\n"
	"O   -0.8  HOH_pe\n"
	"H    0.4  HOH_pe\n"
	"Na   1.0  NA_npe\n"
	"H    0.5  HOH_pe\n"
	"X   -0.5  HOH_pe\n"
	"@end\n"
	"@polarizabilities\n"
	"O   6.0 0.0 0.0 6.0 0.0 6.0  HOH_pe\n"
	"H   2.0 0.0 0.0 2.0 0.0 2.0  HOH_pe\n"
	"H   2.0 0.0 0.0 2.0 0.0 2.0  HOH_pe\n"
	"X   1.0 0.3 0.2 2.0 0.1 6.0  HOH_pe\n"
	"@end\n";

UTEST(vlx, pe_environment_by_file_name_in_atomic_units) {
	const md_str_t h5  = MD_STR_LIT(VLX_PE_DIR "unittest_pe_custom.h5");
	const md_str_t pot = MD_STR_LIT(VLX_PE_DIR "unittest_pe_custom.pot");
	ASSERT_TRUE(vlx_test_write_file(pot, vlx_test_custom_pot, sizeof(vlx_test_custom_pot) - 1));
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, "/cluster/scratch/run42/unittest_pe_custom.pot"));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MD_MEGABYTES(64));
	remove(h5.ptr);
	remove(pot.ptr);
	ASSERT_TRUE(loaded);

	ASSERT_EQ(3u + 5u, t.sys.atom.count);

	const double bohr = 0.5291772109029999;
	EXPECT_NEAR(20.0 * bohr, t.state.xyz[3].z, 1e-5);
	EXPECT_NEAR( 1.0 * bohr, t.state.xyz[5].x, 1e-5);

	// The expansion point is a virtual site, not an atom of an unknown element
	EXPECT_EQ(0, md_atom_atomic_number(&t.sys.atom, 7));
	EXPECT_EQ(MD_PARTICLE_VIRTUAL_SITE, md_atom_particle_kind(&t.sys.atom, 7));
	EXPECT_EQ(MD_PARTICLE_ATOM, md_atom_particle_kind(&t.sys.atom, 4));

	ASSERT_EQ(3u, t.sys.residue.count);
	EXPECT_VALID_SYSTEM(&t.sys);
	ASSERT_TRUE(md_system_infer(&t.sys, &t.state, NULL));
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, 1), MD_STR_LIT("NA")));
	EXPECT_EQ(4, md_system_residue_seq_id(&t.sys, 1));
	EXPECT_EQ(MD_RESIDUE_KIND_ION, md_system_residue_kind(&t.sys, 1));
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, 2), MD_STR_LIT("HOH")));

	double q[8], a[8];
	ASSERT_EQ(8u, qm_test_series(q, 8, &t, MD_STR_LIT("atom/charge")));
	ASSERT_EQ(8u, qm_test_series(a, 8, &t, MD_STR_LIT("atom/polarizability")));
	EXPECT_NEAR( 1.0, q[3], 1e-12);
	EXPECT_NEAR(-0.8, q[4], 1e-12);
	EXPECT_NEAR( 0.4, q[5], 1e-12);
	EXPECT_NEAR( 0.5, q[6], 1e-12);	// The HOH_pe rows interleaved with NA_npe's still go in order
	EXPECT_NEAR(-0.5, q[7], 1e-12);
	EXPECT_EQ(0.0, a[3]);
	EXPECT_NEAR(6.0, a[4], 1e-12);
	EXPECT_NEAR(3.0, a[7], 1e-12);	// (1 + 2 + 6) / 3, the off diagonal plays no part

	qm_test_free(&t);
}

// The QM region and its environment as script selections: 'qm' the calculation's own atoms, 'environment'
// the sites of the embedding, both one selection per component, and either a compile error where its
// region is not there
UTEST(vlx, pe_regions_are_selectable) {
	const md_str_t h5  = MD_STR_LIT(VLX_PE_DIR "unittest_pe_regions.h5");
	const md_str_t pot = MD_STR_LIT(VLX_PE_DIR "unittest_pe_regions.pot");
	ASSERT_TRUE(vlx_test_write_file(pot, vlx_test_custom_pot, sizeof(vlx_test_custom_pot) - 1));
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, "unittest_pe_regions.pot"));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MD_MEGABYTES(64));
	remove(h5.ptr);
	remove(pot.ptr);
	ASSERT_TRUE(loaded);
	ASSERT_EQ(3u + 5u, t.sys.atom.count);

	for (size_t i = 0; i < t.sys.atom.count; ++i) {
		EXPECT_EQ(i < 3, md_atom_is_qm(&t.sys.atom, i));
	}

	const md_filter_context_t fctx = {&t.sys, &t.state, NULL};
	md_bitset_t bf = {0};

	ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("qm"), &fctx, t.alloc).success);
	EXPECT_EQ(3u, md_bitset_count(bf));
	EXPECT_EQ(3u, md_bitset_count_range(bf, 0, 3));

	ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("environment"), &fctx, t.alloc).success);
	EXPECT_EQ(5u, md_bitset_count(bf));
	EXPECT_EQ(5u, md_bitset_count_range(bf, 3, 8));

	ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("not qm"), &fctx, t.alloc).success);
	EXPECT_EQ(5u, md_bitset_count(bf));

	// One selection per fragment: NA, and HOH with its expansion point
	md_array(md_bitset_t) arr = 0;
	ASSERT_TRUE(md_filter_sets(&arr, MD_STR_LIT("environment"), &fctx, t.alloc).success);
	ASSERT_EQ(2u, md_array_size(arr));
	EXPECT_EQ(1u, md_bitset_count(arr[0]));
	EXPECT_EQ(4u, md_bitset_count(arr[1]));

	arr = 0;
	ASSERT_TRUE(md_filter_sets(&arr, MD_STR_LIT("qm"), &fctx, t.alloc).success);
	ASSERT_EQ(1u, md_array_size(arr));
	EXPECT_EQ(3u, md_bitset_count(arr[0]));

	// They compose like any other selection, and take a context like the residue selectors do
	EXPECT_TRUE(md_filter(&bf, MD_STR_LIT("within(100, qm) and environment"), &fctx, t.alloc).success);
	EXPECT_EQ(5u, md_bitset_count(bf));
	EXPECT_TRUE(md_filter(&bf, MD_STR_LIT("environment in resname('NA')"), &fctx, t.alloc).success);
	EXPECT_EQ(1u, md_bitset_count(bf));
	EXPECT_TRUE(md_bitset_test(bf, 3));

	qm_test_free(&t);
}

// The reference potential, 50 waters around a QM water: the examples of the script reference compile
// against it, and the environment is one selection per water
UTEST(vlx, pe_regions_in_a_script) {
	const md_str_t h5 = MD_STR_LIT(VLX_PE_DIR "unittest_pe_script.h5");
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, "../pot/water_pe_npe.pot"));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MD_MEGABYTES(64));
	remove(h5.ptr);
	ASSERT_TRUE(loaded);

	md_script_ir_t* ir = md_script_ir_create(t.alloc);
	const md_str_t src = MD_STR_LIT(
		"d = distance_min(qm(), water());\n"
		"near_qm = within(5, qm());\n"
		"first_shell = within(3.5, qm()) and environment();\n"
		"n_fragments = count(environment(), \"residue\");\n");
	EXPECT_TRUE(md_script_ir_compile_from_source(ir, src, &t.sys, NULL));
	EXPECT_TRUE(md_script_ir_valid(ir));
	for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
		printf("%.*s\n", (int)md_script_ir_errors(ir)[i].text.len, md_script_ir_errors(ir)[i].text.ptr);
	}
	md_script_ir_free(ir);

	const md_filter_context_t fctx = {&t.sys, &t.state, NULL};
	md_array(md_bitset_t) arr = 0;
	ASSERT_TRUE(md_filter_sets(&arr, MD_STR_LIT("environment"), &fctx, t.alloc).success);
	EXPECT_EQ(50u, md_array_size(arr));

	qm_test_free(&t);
}

// The fragments become components in the order of their numbers, whatever order the file lists them
// in. VeloxChem writes the polarizable fragments first and the rest after, so a peptide that the
// polarizable region cuts through arrives scattered: here residue 11 (polarizable), a water, then
// residues 10 and 12. In file order that is no molecule and no backbone; in number order it is a
// tripeptide, one instance with one backbone over all three residues. The parameters go with their
// sites, so residue 11 keeps the polarizable type's charges in the middle of the chain.
UTEST(vlx, pe_fragments_follow_their_numbers) {
	static const char pot_text[] =
		"@environment\n"
		"units: angstrom\n"
		"xyz:\n"
		"N    21.463    0.376   -2.396  GLY_pe 11 N\n"
		"C    21.899    0.981   -3.649  GLY_pe 11 CA\n"
		"C    21.768    2.500   -3.602  GLY_pe 11 C\n"
		"O    22.693    3.219   -3.981  GLY_pe 11 O\n"
		"O    30.000    0.000    0.000  HOH_pe 997 OW\n"
		"H    30.957    0.000    0.000  HOH_pe 997 HW1\n"
		"H    29.760    0.927    0.000  HOH_pe 997 HW2\n"
		"N    20.000    0.000    0.000  GLY_npe 10 N\n"
		"C    21.458    0.000    0.000  GLY_npe 10 CA\n"
		"C    22.009    0.711   -1.231  GLY_npe 10 C\n"
		"O    22.910    1.543   -1.121  GLY_npe 10 O\n"
		"N    20.618    2.976   -3.137  GLY_npe 12 N\n"
		"C    20.364    4.408   -3.041  GLY_npe 12 CA\n"
		"C    21.421    5.099   -2.187  GLY_npe 12 C\n"
		"O    21.958    6.137   -2.575  GLY_npe 12 O\n"
		"@end\n"
		"@charges\n"
		"N   -0.41  GLY_pe\n"
		"C    0.02  GLY_pe\n"
		"C    0.53  GLY_pe\n"
		"O   -0.49  GLY_pe\n"
		"O   -0.80  HOH_pe\n"
		"H    0.40  HOH_pe\n"
		"H    0.40  HOH_pe\n"
		"N   -0.31  GLY_npe\n"
		"C    0.12  GLY_npe\n"
		"C    0.43  GLY_npe\n"
		"O   -0.39  GLY_npe\n"
		"@end\n"
		"@polarizabilities\n"
		"N   7.0 0.0 0.0 7.0 0.0 7.0  GLY_pe\n"
		"C   8.0 0.0 0.0 8.0 0.0 8.0  GLY_pe\n"
		"C   9.0 0.0 0.0 9.0 0.0 9.0  GLY_pe\n"
		"O   6.5 0.0 0.0 6.5 0.0 6.5  GLY_pe\n"
		"O   6.0 0.0 0.0 6.0 0.0 6.0  HOH_pe\n"
		"H   2.0 0.0 0.0 2.0 0.0 2.0  HOH_pe\n"
		"H   2.0 0.0 0.0 2.0 0.0 2.0  HOH_pe\n"
		"@end\n";

	const md_str_t h5  = MD_STR_LIT(VLX_PE_DIR "unittest_pe_order.h5");
	const md_str_t pot = MD_STR_LIT(VLX_PE_DIR "unittest_pe_order.pot");
	ASSERT_TRUE(vlx_test_write_file(pot, pot_text, sizeof(pot_text) - 1));
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, "unittest_pe_order.pot"));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MD_MEGABYTES(64));
	remove(h5.ptr);
	remove(pot.ptr);
	ASSERT_TRUE(loaded);
	ASSERT_EQ(3u + 15u, t.sys.atom.count);

	// QM, then residues 10, 11, 12 and the water
	ASSERT_EQ(5u, t.sys.residue.count);
	EXPECT_VALID_SYSTEM(&t.sys);
	const int seq[5] = {0, 10, 11, 12, 997};
	const uint32_t beg[5] = {0, 3, 7, 11, 15};
	for (size_t c = 0; c < 5; ++c) {
		EXPECT_EQ(seq[c], md_system_residue_seq_id(&t.sys, c));
		EXPECT_EQ(beg[c], md_system_residue_atom_range(&t.sys, c).beg);
	}
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, 2), MD_STR_LIT("GLY")));
	EXPECT_TRUE(md_str_eq(md_system_residue_name(&t.sys, 4), MD_STR_LIT("HOH")));

	// A fragment's sites keep the order the file gives them in
	EXPECT_TRUE(md_str_eq(md_atom_name(&t.sys.atom, 3), MD_STR_LIT("N")));
	EXPECT_TRUE(md_str_eq(md_atom_name(&t.sys.atom, 6), MD_STR_LIT("O")));
	EXPECT_NEAR(20.000, t.state.xyz[3].x, 1e-5);
	EXPECT_NEAR(21.463, t.state.xyz[7].x, 1e-5);
	EXPECT_NEAR(20.618, t.state.xyz[11].x, 1e-5);
	EXPECT_NEAR(30.000, t.state.xyz[15].x, 1e-5);

	double q[18], a[18];
	ASSERT_EQ(18u, qm_test_series(q, 18, &t, MD_STR_LIT("atom/charge")));
	ASSERT_EQ(18u, qm_test_series(a, 18, &t, MD_STR_LIT("atom/polarizability")));
	EXPECT_NEAR(-0.31, q[3],  1e-12);	// Residue 10, non-polarizable
	EXPECT_NEAR( 0.43, q[5],  1e-12);
	EXPECT_NEAR(-0.41, q[7],  1e-12);	// Residue 11, polarizable
	EXPECT_NEAR( 0.53, q[9],  1e-12);
	EXPECT_NEAR(-0.49, q[10], 1e-12);
	EXPECT_NEAR(-0.31, q[11], 1e-12);	// Residue 12
	EXPECT_NEAR(-0.80, q[15], 1e-12);	// The water
	EXPECT_NEAR( 0.40, q[17], 1e-12);
	EXPECT_EQ(0.0, a[3]);
	EXPECT_NEAR(7.0, a[7],  1e-12);
	EXPECT_NEAR(9.0, a[9],  1e-12);
	EXPECT_EQ(0.0, a[14]);
	EXPECT_NEAR(6.0, a[15], 1e-12);

	// What viamd infers from it: the three residues one peptide, with one backbone through them all
	ASSERT_TRUE(md_system_infer(&t.sys, &t.state, NULL));
	EXPECT_VALID_SYSTEM(&t.sys);
	const md_part_idx_t part = md_system_part_find_by_residue_idx(&t.sys, 1);
	ASSERT_TRUE(part >= 0);
	EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(&t.sys, part));
	EXPECT_EQ(3u, md_system_part_residue_range(&t.sys, part).end - md_system_part_residue_range(&t.sys, part).beg);
	ASSERT_EQ(1u, t.sys.protein_backbone.range.count);
	ASSERT_EQ(3u, t.sys.protein_backbone.segment.count);
	for (size_t s = 0; s < 3; ++s) {
		EXPECT_EQ((int)(1 + s), t.sys.protein_backbone.segment.res_idx[s]);
	}

	qm_test_free(&t);
}

// A calculation without an embedding is QM throughout: 'qm' is all of it and 'environment' is no selection at all
UTEST(vlx, qm_without_environment) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(VLX_PE_DIR "h2o.h5"), MD_MEGABYTES(64)));

	const md_filter_context_t fctx = {&t.sys, &t.state, NULL};
	md_bitset_t bf = {0};
	ASSERT_TRUE(md_filter(&bf, MD_STR_LIT("qm"), &fctx, t.alloc).success);
	EXPECT_EQ(3u, md_bitset_count(bf));

	const md_filter_result_t res = md_filter(&bf, MD_STR_LIT("environment"), &fctx, t.alloc);
	EXPECT_FALSE(res.success);
	EXPECT_TRUE(strstr(res.error, "QM throughout") != NULL);

	qm_test_free(&t);
}

// What VeloxChem would have refused is not the potential the calculation ran with, so none of it is
// used. A potential that cannot be found leaves the environment out. Neither is a failed load.
UTEST(vlx, pe_environment_left_out_when_it_cannot_be_the_one) {
	static const char bad_pot[] =
		"@environment\n"
		"xyz:\n"
		"O  0.0 0.0 0.0  HOH_pe 1 OW\n"
		"H  1.0 0.0 0.0  HOH_pe 1 H1\n"
		"H -1.0 0.0 0.0  HOH_pe 1 H2\n"
		"@end\n"
		"@charges\n"
		"O -0.8 HOH_pe\n"
		"H  0.4 HOH_pe\n"
		"@end\n";

	const md_str_t pot     = MD_STR_LIT(VLX_PE_DIR "unittest_pe_bad.pot");
	const md_str_t h5_bad  = MD_STR_LIT(VLX_PE_DIR "unittest_pe_bad.h5");
	const md_str_t h5_none = MD_STR_LIT(VLX_PE_DIR "unittest_pe_missing.h5");
	ASSERT_TRUE(vlx_test_write_file(pot, bad_pot, sizeof(bad_pot) - 1));
	ASSERT_TRUE(vlx_test_write_pe_h5(h5_bad,  "unittest_pe_bad.pot"));
	ASSERT_TRUE(vlx_test_write_pe_h5(h5_none, "unittest_pe_not_there.pot"));

	vlx_test_t bad = {0}, none = {0};
	const bool loaded_bad  = vlx_test_load(&bad,  h5_bad,  MD_MEGABYTES(64));
	const bool loaded_none = vlx_test_load(&none, h5_none, MD_MEGABYTES(64));
	remove(h5_bad.ptr);
	remove(h5_none.ptr);
	remove(pot.ptr);

	EXPECT_TRUE(loaded_bad);
	EXPECT_TRUE(loaded_none);
	EXPECT_EQ(3u, bad.sys.atom.count);
	EXPECT_EQ(3u, none.sys.atom.count);
	EXPECT_EQ(0u, bad.sys.residue.count);
	EXPECT_FALSE(qm_test_has(&bad,  MD_STR_LIT("atom/charge")));
	EXPECT_FALSE(qm_test_has(&none, MD_STR_LIT("atom/charge")));
	EXPECT_EQ(3u, qm_test_count(&none, MD_STR_LIT("atom/nuclear_charges")));

	qm_test_free(&bad);
	qm_test_free(&none);
}

// A supplemental load leaves the atoms of the system it supplements alone, so it adds no sites
UTEST(vlx, pe_environment_not_added_by_a_supplemental_load) {
	const md_str_t h5 = MD_STR_LIT(VLX_PE_DIR "unittest_pe_supplement.h5");
	ASSERT_TRUE(vlx_test_write_pe_h5(h5, "../pot/water_pe_npe.pot"));

	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, MD_STR_LIT(VLX_PE_DIR "h2o.h5"), MD_MEGABYTES(64));
	const bool supplemented = loaded && md_vlx_system_supplement_from_file(&t.sys, h5, NULL);
	remove(h5.ptr);

	EXPECT_TRUE(supplemented);
	EXPECT_EQ(3u, t.sys.atom.count);
	EXPECT_EQ(3u, t.state.num_atoms);
	EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("atom/charge")));
	// Without a map the supplement's atoms are the system's own, and they stay its QM region
	for (size_t i = 0; i < t.sys.atom.count; ++i) {
		EXPECT_TRUE(md_atom_is_qm(&t.sys.atom, i));
	}

	qm_test_free(&t);
}

// ---------------------------------------------------------------------------
// INDUCED DIPOLES
//
// A fika embedding run of VeloxChem keeps the induced dipoles it converged to,
// scf/embedding_induced_dipoles [P][3] in e bohr: one per POLARIZABLE site of the environment, in
// the order of the .pot's rows - not the order the sites become atoms in, which follows the fragment
// numbers. Over the system's atoms a non-polarizable site has a zero dipole.
// ---------------------------------------------------------------------------

// A copy of h2o.h5 naming 'potfile', with scf/embedding_induced_dipoles [rows][cols] beside it
static bool vlx_test_write_induced_dipoles_h5(md_str_t dst, const char* potfile, const double* mu, size_t rows, size_t cols) {
	if (!vlx_test_write_pe_h5(dst, potfile)) return false;

	char path[1024];
	md_str_copy_to_char_buf(path, sizeof(path), dst);
	hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0) return false;
	hid_t scf = H5Gopen(file, "scf", H5P_DEFAULT);
	const hsize_t dims[2] = { rows, cols };
	hid_t space = H5Screate_simple(2, dims, NULL);
	hid_t dset  = scf >= 0 ? H5Dcreate2(scf, "embedding_induced_dipoles", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
	const bool ok = dset >= 0 && H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mu) >= 0;
	if (dset >= 0) H5Dclose(dset);
	H5Sclose(space);
	if (scf >= 0) H5Gclose(scf);
	H5Fclose(file);
	return ok;
}

// Three waters, listed out of their numbers the way VeloxChem writes them (polarizable first): 3, 2,
// then 1. The atoms follow the numbers, 1, 2, 3; the dipoles follow the file, 3, 2. Water 1 is _npe
// and so not polarizable, though its type has polarizabilities, and the second hydrogen of a _pe
// water has none: the polarizable sites are the oxygen and the first hydrogen of waters 3 and 2.
static const char vlx_test_induced_pot[] =
	"@environment\n"
	"units: angstrom\n"
	"xyz:\n"
	"O   3.000  0.000  0.000  HOH_pe 3 OW\n"
	"H   3.957  0.000  0.000  HOH_pe 3 H1\n"
	"H   2.760  0.927  0.000  HOH_pe 3 H2\n"
	"O   6.000  0.000  0.000  HOH_pe 2 OW\n"
	"H   6.957  0.000  0.000  HOH_pe 2 H1\n"
	"H   5.760  0.927  0.000  HOH_pe 2 H2\n"
	"O   9.000  0.000  0.000  HOH_npe 1 OW\n"
	"H   9.957  0.000  0.000  HOH_npe 1 H1\n"
	"H   8.760  0.927  0.000  HOH_npe 1 H2\n"
	"@end\n"
	"@charges\n"
	"O  -0.80  HOH_pe\n"
	"H   0.40  HOH_pe\n"
	"H   0.40  HOH_pe\n"
	"O  -0.83  HOH_npe\n"
	"H   0.415 HOH_npe\n"
	"H   0.415 HOH_npe\n"
	"@end\n"
	"@polarizabilities\n"
	"O   6.0 0.0 0.0 6.0 0.0 6.0  HOH_pe\n"
	"H   2.0 0.0 0.0 2.0 0.0 2.0  HOH_pe\n"
	"H   0.0 0.0 0.0 0.0 0.0 0.0  HOH_pe\n"
	"O   5.0 0.0 0.0 5.0 0.0 5.0  HOH_npe\n"
	"H   1.0 0.0 0.0 1.0 0.0 1.0  HOH_npe\n"
	"H   1.0 0.0 0.0 1.0 0.0 1.0  HOH_npe\n"
	"@end\n";

// The dipoles of the polarizable sites in file order, each naming its site: (fragment number, index in
// the fragment, 0.5)
static const double vlx_test_induced_mu[4][3] = {
	{3, 0, 0.5},   {3, 1, 0.5},
	{2, 0, 0.5},   {2, 1, 0.5},
};

UTEST(vlx, pe_induced_dipoles_land_on_their_sites) {
	const md_str_t h5  = MD_STR_LIT(VLX_PE_DIR "unittest_pe_induced.h5");
	const md_str_t pot = MD_STR_LIT(VLX_PE_DIR "unittest_pe_induced.pot");
	ASSERT_TRUE(vlx_test_write_file(pot, vlx_test_induced_pot, sizeof(vlx_test_induced_pot) - 1));
	ASSERT_TRUE(vlx_test_write_induced_dipoles_h5(h5, "unittest_pe_induced.pot", &vlx_test_induced_mu[0][0], 4, 3));
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, h5, MD_MEGABYTES(64));
	remove(h5.ptr);
	remove(pot.ptr);
	ASSERT_TRUE(loaded);

	const size_t num_qm = 3, num_atoms = num_qm + 9;
	ASSERT_EQ(num_atoms, t.sys.atom.count);

	// As the file holds them, over the polarizable sites
	const md_attribute_t* raw = qm_test_attr(&t, MD_STR_LIT("vlx/scf/embedding/induced_dipole"));
	ASSERT_TRUE(raw != NULL);
	EXPECT_EQ(3u, raw->format.components);
	EXPECT_EQ(1u, raw->format.rank);
	EXPECT_EQ(4u, md_attribute_value_count(&raw->format));
	EXPECT_TRUE(md_unit_equal(raw->unit, md_unit_elementary_charge_bohr()));

	// Over the system's atoms: each atom's dipole is the one its own site was given
	const md_attribute_t* attr = qm_test_attr(&t, MD_STR_LIT("atom/induced_dipole"));
	ASSERT_TRUE(attr != NULL);
	EXPECT_EQ(3u, attr->format.components);
	EXPECT_EQ(1u, attr->format.rank);
	ASSERT_EQ(num_atoms, md_attribute_value_count(&attr->format));
	EXPECT_TRUE(md_unit_equal(attr->unit, md_unit_elementary_charge_bohr()));

	double mu[3 * 12];
	ASSERT_EQ(3 * num_atoms, qm_test_series(mu, MD_COUNTOF(mu), &t, MD_STR_LIT("atom/induced_dipole")));
	for (size_t i = 0; i < 3 * num_qm; ++i) {
		EXPECT_TRUE(vlx_test_absent(mu[i]));	// The QM atoms are no sites
	}
	for (size_t i = num_qm; i < num_atoms; ++i) {
		const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&t.sys, i);
		ASSERT_TRUE(res >= 0);
		const int    frag  = md_system_residue_seq_id(&t.sys, res);
		const size_t local = i - md_system_residue_atom_range(&t.sys, res).beg;
		if (frag == 1 || local == 2) {
			EXPECT_EQ(0.0, mu[3 * i + 0]);
			EXPECT_EQ(0.0, mu[3 * i + 1]);
			EXPECT_EQ(0.0, mu[3 * i + 2]);
		} else {
			EXPECT_EQ((double)frag,  mu[3 * i + 0]);
			EXPECT_EQ((double)local, mu[3 * i + 1]);
			EXPECT_EQ(0.5,           mu[3 * i + 2]);
		}
	}

	qm_test_free(&t);
}

// The run's induced dipoles count the polarizable sites of the potential it ran with. A potential
// beside the .h5 with another number of them is not that one - one per SITE is not that either - and its environment is left out - while the dipoles,
// which are the file's own, are still published over its sites. A dataset that is not a 3-vector per
// site is not read at all. Neither is a failed load.
UTEST(vlx, pe_induced_dipoles_must_match_the_potential) {
	const md_str_t pot      = MD_STR_LIT(VLX_PE_DIR "unittest_pe_induced_count.pot");
	const md_str_t h5_count = MD_STR_LIT(VLX_PE_DIR "unittest_pe_induced_count.h5");
	const md_str_t h5_shape = MD_STR_LIT(VLX_PE_DIR "unittest_pe_induced_shape.h5");
	ASSERT_TRUE(vlx_test_write_file(pot, vlx_test_induced_pot, sizeof(vlx_test_induced_pot) - 1));
	static const double per_site[9][3] = {0};	// One per site, the non-polarizable ones included
	ASSERT_TRUE(vlx_test_write_induced_dipoles_h5(h5_count, "unittest_pe_induced_count.pot", &per_site[0][0], 9, 3));
	ASSERT_TRUE(vlx_test_write_induced_dipoles_h5(h5_shape, "unittest_pe_induced_count.pot", &vlx_test_induced_mu[0][0], 6, 2));	// 12 values, as [4][3] holds

	vlx_test_t count = {0}, shape = {0};
	const bool loaded_count = vlx_test_load(&count, h5_count, MD_MEGABYTES(64));
	const bool loaded_shape = vlx_test_load(&shape, h5_shape, MD_MEGABYTES(64));
	remove(h5_count.ptr);
	remove(h5_shape.ptr);
	remove(pot.ptr);

	EXPECT_TRUE(loaded_count);
	EXPECT_EQ(3u, count.sys.atom.count);
	EXPECT_EQ(0u, count.sys.residue.count);
	EXPECT_FALSE(qm_test_has(&count, MD_STR_LIT("atom/charge")));
	EXPECT_FALSE(qm_test_has(&count, MD_STR_LIT("atom/induced_dipole")));
	const md_attribute_t* raw = qm_test_attr(&count, MD_STR_LIT("vlx/scf/embedding/induced_dipole"));
	ASSERT_TRUE(raw != NULL);
	EXPECT_EQ(9u, md_attribute_value_count(&raw->format));

	// Without a usable count there is nothing to check the potential against: it is added on its own merits
	EXPECT_TRUE(loaded_shape);
	EXPECT_EQ(3u + 9u, shape.sys.atom.count);
	EXPECT_TRUE(qm_test_has(&shape, MD_STR_LIT("atom/charge")));
	EXPECT_FALSE(qm_test_has(&shape, MD_STR_LIT("atom/induced_dipole")));
	EXPECT_FALSE(qm_test_has(&shape, MD_STR_LIT("vlx/scf/embedding/induced_dipole")));

	qm_test_free(&count);
	qm_test_free(&shape);
}

// A run of the experimental embedding branch (CAM-B3LYP/def2-SVP, 10 QM atoms, FIKA embedding) with
// its potential where the run named it, pot_frames/pe_frame_000000.pot beside the .h5: 769 sites, of
// which the 48 of fragments 6 (LYS), 8 (ILE) and 9 (GLY) are polarizable - one per induced dipole.
UTEST(vlx, pe_induced_dipoles_from_file) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/pe_induced_dipoles.h5"), MD_MEGABYTES(64)));

	const size_t num_qm = 10, num_atoms = num_qm + 769;
	ASSERT_EQ(num_atoms, t.sys.atom.count);

	// As written
	const md_attribute_t* raw = qm_test_attr(&t, MD_STR_LIT("vlx/scf/embedding/induced_dipole"));
	ASSERT_TRUE(raw != NULL);
	EXPECT_EQ(3u, raw->format.components);
	ASSERT_EQ(48u, md_attribute_value_count(&raw->format));
	EXPECT_TRUE(md_unit_equal(raw->unit, md_unit_elementary_charge_bohr()));

	double mu[3 * 48];
	ASSERT_EQ(MD_COUNTOF(mu), qm_test_series(mu, MD_COUNTOF(mu), &t, MD_STR_LIT("vlx/scf/embedding/induced_dipole")));
	EXPECT_NEAR( 0.0558321,  mu[0], 1e-7);
	EXPECT_NEAR( 0.00781969, mu[1], 1e-8);
	EXPECT_NEAR( 0.05448396, mu[2], 1e-8);
	EXPECT_NEAR( 0.00200192, mu[3 * 47 + 0], 1e-8);
	EXPECT_NEAR(-0.1647142,  mu[3 * 47 + 1], 1e-7);
	EXPECT_NEAR( 0.02673146, mu[3 * 47 + 2], 1e-8);

	// Over the system: the dipoles in order over the atoms of fragments 6, 8 and 9, zero on every
	// other site, nothing on the QM atoms
	double* atom_mu = (double*)md_alloc(t.alloc, sizeof(double) * 3 * num_atoms);
	double* polar   = (double*)md_alloc(t.alloc, sizeof(double) * num_atoms);
	ASSERT_EQ(3 * num_atoms, qm_test_series(atom_mu, 3 * num_atoms, &t, MD_STR_LIT("atom/induced_dipole")));
	ASSERT_EQ(num_atoms, qm_test_series(polar, num_atoms, &t, MD_STR_LIT("atom/polarizability")));
	size_t next = 0;
	for (size_t i = 0; i < num_atoms; ++i) {
		if (i < num_qm) {
			EXPECT_TRUE(vlx_test_absent(atom_mu[3 * i]));
			continue;
		}
		const int frag = md_system_residue_seq_id(&t.sys, md_system_residue_find_by_atom_idx(&t.sys, i));
		const bool polarizable = frag == 6 || frag == 8 || frag == 9;
		EXPECT_EQ(polarizable, polar[i] > 0.0);
		for (int c = 0; c < 3; ++c) {
			EXPECT_EQ(polarizable ? mu[3 * next + c] : 0.0, atom_mu[3 * i + c]);
		}
		next += polarizable;
	}
	EXPECT_EQ(48u, next);

	// LYS 6 opens the polarizable dipoles with its backbone N
	const md_residue_idx_t lys = md_system_residue_find_by_atom_idx(&t.sys, num_qm + 109);
	EXPECT_EQ(6, md_system_residue_seq_id(&t.sys, lys));
	EXPECT_TRUE(md_str_eq(md_atom_name(&t.sys.atom, num_qm + 109), MD_STR_LIT("N")));
	EXPECT_EQ(mu[0], atom_mu[3 * (num_qm + 109)]);

	qm_test_free(&t);
}

// ---------------------------------------------------------------------------
// RESPONSE VECTORS AND NATURAL TRANSITION ORBITALS, ON DEMAND
//
// Nothing of the response vectors is read at load: vlx/rsp/solution_matrix reads the states asked for
// from the file, and the NTOs are computed a state at a time from it. acro-xps.h5 keeps one dataset per
// state (rsp/S1, ...), pe_induced_dipoles.h5 the rows of rsp/full_solutions_matrix; a copy of the latter
// with the matrix chunked and deflated is not one run of bytes and is read through HDF5.
// ---------------------------------------------------------------------------

// One state's response vector straight from the file, through HDF5 (heap, *out_len values)
static double* vlx_test_h5_solution(size_t* out_len, const char* path, uint32_t state) {
	double* out = NULL;
	hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file < 0) return NULL;
	char name[32];
	snprintf(name, sizeof(name), "rsp/S%u", state + 1);
	const bool per_state = H5Lexists(file, name, H5P_DEFAULT) > 0;
	hid_t dset  = H5Dopen2(file, per_state ? name : "rsp/full_solutions_matrix", H5P_DEFAULT);
	hid_t space = dset >= 0 ? H5Dget_space(dset) : -1;
	hsize_t dims[2] = {0};
	const int rank = space >= 0 ? H5Sget_simple_extent_dims(space, dims, NULL) : -1;
	const hsize_t len = per_state ? dims[0] : dims[1];
	if (rank == (per_state ? 1 : 2) && len > 0 && (per_state || dims[0] > state)) {
		out = (double*)malloc(sizeof(double) * len);
		hid_t mem = H5Screate_simple(1, &len, NULL);
		if (!per_state) {
			const hsize_t start[2] = { state, 0 }, count[2] = { 1, len };
			H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL, count, NULL);
		}
		if (H5Dread(dset, H5T_NATIVE_DOUBLE, mem, per_state ? H5S_ALL : space, H5P_DEFAULT, out) < 0) {
			free(out);
			out = NULL;
		}
		H5Sclose(mem);
		*out_len = (size_t)len;
	}
	if (space >= 0) H5Sclose(space);
	if (dset >= 0) H5Dclose(dset);
	H5Fclose(file);
	return out;
}

// A copy of 'src' with the 2D dataset 'name' rewritten chunked by row, and deflated where HDF5 has it: not one run of
// bytes, so read through HDF5
static bool vlx_test_write_chunked_h5(md_str_t src, md_str_t dst, const char* name) {
	md_file_t in = {0};
	if (!md_file_open(&in, src, MD_FILE_READ)) return false;
	const size_t size = (size_t)md_file_size(in);
	void* bytes = md_alloc(md_get_heap_allocator(), size);
	const bool read = md_file_read(in, bytes, size) == size;
	md_file_close(&in);
	const bool written = read && vlx_test_write_file(dst, bytes, size);
	md_free(md_get_heap_allocator(), bytes, size);
	if (!written) return false;

	char path[1024];
	md_str_copy_to_char_buf(path, sizeof(path), dst);
	hid_t file = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0) return false;
	bool ok = false;
	hid_t dset = H5Dopen2(file, name, H5P_DEFAULT);
	if (dset >= 0) {
		hid_t space = H5Dget_space(dset);
		hsize_t dims[2] = {0};
		if (H5Sget_simple_extent_dims(space, dims, NULL) == 2) {
			double* data = (double*)malloc(sizeof(double) * dims[0] * dims[1]);
			ok = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) >= 0;
			H5Dclose(dset);
			dset = -1;
			ok = ok && H5Ldelete(file, name, H5P_DEFAULT) >= 0;
			hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
			const hsize_t chunk[2] = { 1, dims[1] };
			H5Pset_chunk(dcpl, 2, chunk);
			if (H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0) H5Pset_deflate(dcpl, 6);
			hid_t out = ok ? H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT) : -1;
			ok = out >= 0 && H5Dwrite(out, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) >= 0;
			if (out >= 0) H5Dclose(out);
			H5Pclose(dcpl);
			free(data);
		}
		H5Sclose(space);
		if (dset >= 0) H5Dclose(dset);
	}
	H5Fclose(file);
	return ok;
}

// Every state's vector, by state, by element and whole, is the file's, and the attribute holds none of it
static void vlx_test_check_solutions(int32_t* utest_result, const vlx_test_t* t, const char* file) {
	const md_attribute_t* sol = qm_test_attr(t, MD_STR_LIT("vlx/rsp/solution_matrix"));
	ASSERT_TRUE(sol != NULL);
	ASSERT_TRUE(sol->data == NULL && sol->virt.provider != NULL);
	const uint32_t num_states = sol->format.shape[0];
	const size_t   len        = sol->format.shape[1];
	ASSERT_TRUE(num_states > 0 && len > 0);

	double* row   = (double*)malloc(sizeof(double) * len);
	double* whole = (double*)malloc(sizeof(double) * len * num_states);
	ASSERT_EQ(md_attribute_extract_f64(whole, len * num_states, &t->sys.attributes, sol, md_attribute_slice_all(), md_unit_none()), len * num_states);
	for (uint32_t s = 0; s < num_states; ++s) {
		size_t ref_len = 0;
		double* ref = vlx_test_h5_solution(&ref_len, file, s);
		ASSERT_TRUE(ref != NULL);
		ASSERT_EQ(ref_len, len);
		ASSERT_EQ(md_attribute_extract_f64(row, len, &t->sys.attributes, sol, md_attribute_slice_1(s), md_unit_none()), len);
		EXPECT_EQ(0, memcmp(row, ref, sizeof(double) * len));
		EXPECT_EQ(0, memcmp(whole + s * len, ref, sizeof(double) * len));
		double elem = 0.0;
		ASSERT_EQ(md_attribute_extract_f64(&elem, 1, &t->sys.attributes, sol, md_attribute_slice_2(s, (uint32_t)(len / 2)), md_unit_none()), 1u);
		EXPECT_EQ(ref[len / 2], elem);
		free(ref);
	}
	free(row);
	free(whole);
}

UTEST(vlx, response_vectors_are_read_on_demand) {
	const char* files[] = { MD_UNITTEST_DATA_DIR "/vlx/acro-xps.h5", MD_UNITTEST_DATA_DIR "/vlx/pe_induced_dipoles.h5" };
	for (size_t f = 0; f < MD_COUNTOF(files); ++f) {
		vlx_test_t t = {0};
		ASSERT_TRUE(vlx_test_load(&t, md_str_from_cstr(files[f]), MD_MEGABYTES(64)));
		vlx_test_check_solutions(utest_result, &t, files[f]);
		qm_test_free(&t);
	}
}

UTEST(vlx, response_vectors_chunked_are_read_through_hdf5) {
	const md_str_t plain_path = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/pe_induced_dipoles.h5");
	const md_str_t chunk_path = MD_STR_LIT(VLX_PE_DIR "unittest_rsp_chunked.h5");
	ASSERT_TRUE(vlx_test_write_chunked_h5(plain_path, chunk_path, "rsp/full_solutions_matrix"));

	vlx_test_t plain = {0}, chunked = {0};
	const bool loaded = vlx_test_load(&plain, plain_path, MD_MEGABYTES(64)) && vlx_test_load(&chunked, chunk_path, MD_MEGABYTES(64));
	if (loaded) {
		vlx_test_check_solutions(utest_result, &chunked, chunk_path.ptr);
	}

	// the same NTO weights, every state
	const md_attribute_t* la = loaded ? qm_test_attr(&plain, MD_STR_LIT("vlx/rsp/nto/lambda")) : NULL;
	const md_attribute_t* lb = loaded ? qm_test_attr(&chunked, MD_STR_LIT("vlx/rsp/nto/lambda")) : NULL;
	size_t n = la ? md_attribute_element_count(&la->format) : 0;
	double* a = (double*)malloc(sizeof(double) * (n + 1));
	double* b = (double*)malloc(sizeof(double) * (n + 1));
	const size_t got_a = la ? md_attribute_extract_f64(a, n, &plain.sys.attributes, la, md_attribute_slice_all(), md_unit_none()) : 0;
	const size_t got_b = lb ? md_attribute_extract_f64(b, n, &chunked.sys.attributes, lb, md_attribute_slice_all(), md_unit_none()) : 0;
	remove(chunk_path.ptr);
	ASSERT_TRUE(loaded);
	ASSERT_TRUE(n > 0);
	ASSERT_EQ(got_a, n);
	ASSERT_EQ(got_b, n);
	EXPECT_EQ(0, memcmp(a, b, sizeof(double) * n));
	free(a);
	free(b);
	qm_test_free(&plain);
	qm_test_free(&chunked);
}

// ---------------------------------------------------------------------------
// ORBITALS, ON DEMAND
//
// orbital/{alpha,beta}/coefficient are read from C_alpha / C_beta when asked for and converted as they are read; none
// is held. Every way of asking - all of them, one, a run (md_attribute_slice_rows), one element - gives the same
// values, and a file whose C_alpha is chunked (read through HDF5) or that was rewritten since it was loaded gives
// them too.
// ---------------------------------------------------------------------------

// The orbitals as read every way against 'ref', the whole matrix (num_mo x num_ao), when given; the whole read
// otherwise. Returns the whole matrix (heap), for a caller to compare files with.
static double* vlx_test_check_orbitals(int32_t* utest_result, const vlx_test_t* t, const double* ref) {
	const md_attribute_t* c = qm_test_attr(t, MD_STR_LIT("orbital/alpha/coefficient"));
	if (!c || c->data != NULL || c->format.rank != 2) {
		EXPECT_TRUE(c != NULL && c->data == NULL && c->format.rank == 2);
		return NULL;
	}
	const uint32_t nm = c->format.shape[0], na = c->format.shape[1];
	const md_attributes_t* table = &t->sys.attributes;
	double* whole = (double*)malloc(sizeof(double) * nm * na);
	EXPECT_EQ(md_attribute_extract_f64(whole, (size_t)nm * na, table, c, md_attribute_slice_all(), md_unit_none()), (size_t)nm * na);
	if (ref) {
		EXPECT_EQ(0, memcmp(whole, ref, sizeof(double) * nm * na));
	}
	double* row = (double*)malloc(sizeof(double) * nm * na);
	for (uint32_t k = 0; k < nm; k += 3) {
		EXPECT_EQ(md_attribute_extract_f64(row, na, table, c, md_attribute_slice_1(k), md_unit_none()), (size_t)na);
		EXPECT_EQ(0, memcmp(row, whole + (size_t)k * na, sizeof(double) * na));
		double e = 0.0;
		EXPECT_EQ(md_attribute_extract_f64(&e, 1, table, c, md_attribute_slice_2(k, na / 2), md_unit_none()), 1u);
		EXPECT_EQ(whole[(size_t)k * na + na / 2], e);
	}
	// runs: the first five, one from the middle to the end, all of them
	const uint32_t runs[3][2] = { {0, MD_MIN(5u, nm)}, {nm / 2, nm - nm / 2}, {0, nm} };
	for (int r = 0; r < 3; ++r) {
		const md_attribute_slice_t sl = md_attribute_slice_rows(runs[r][0], runs[r][1]);
		md_attribute_format_t f = {0};
		EXPECT_TRUE(md_attribute_slice_format(&f, c, sl));
		EXPECT_EQ(2u, f.rank);
		EXPECT_EQ(runs[r][1], f.shape[0]);
		EXPECT_EQ(na, f.shape[1]);
		const size_t len = (size_t)runs[r][1] * na;
		EXPECT_EQ(md_attribute_extract_f64(row, len, table, c, sl, md_unit_none()), len);
		EXPECT_EQ(0, memcmp(row, whole + (size_t)runs[r][0] * na, sizeof(double) * len));
	}
	// past the end: refused
	EXPECT_EQ(0u, md_attribute_extract_f64(row, (size_t)2 * na, table, c, md_attribute_slice_rows(nm - 1, 2), md_unit_none()));
	free(row);
	return whole;
}

UTEST(vlx, orbitals_are_read_on_demand) {
	static const char* files[] = { MD_UNITTEST_DATA_DIR "/vlx/h2o.h5", MD_UNITTEST_DATA_DIR "/vlx/amide.h5", MD_UNITTEST_DATA_DIR "/vlx/mol.h5" };
	for (int f = 0; f < 3; ++f) {
		vlx_test_t t = {0};
		ASSERT_TRUE(vlx_test_load(&t, md_str_from_cstr(files[f]), MD_MEGABYTES(64)));
		double* whole = vlx_test_check_orbitals(utest_result, &t, NULL);
		EXPECT_TRUE(whole != NULL);
		// a restricted calculation's beta orbitals are alpha's, by a second name
		const md_attribute_t* b = qm_test_attr(&t, MD_STR_LIT("orbital/beta/coefficient"));
		EXPECT_TRUE(b != NULL && md_attribute_same_data(b, qm_test_attr(&t, MD_STR_LIT("orbital/alpha/coefficient"))));
		free(whole);
		qm_test_free(&t);
	}
}

UTEST(vlx, orbitals_chunked_are_read_through_hdf5) {
	const md_str_t plain_path = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5");
	const md_str_t chunk_path = MD_STR_LIT(VLX_PE_DIR "unittest_orbitals_chunked.h5");
	ASSERT_TRUE(vlx_test_write_chunked_h5(plain_path, chunk_path, "scf/C_alpha"));
	vlx_test_t plain = {0}, chunked = {0};
	const bool loaded = vlx_test_load(&plain, plain_path, MD_MEGABYTES(16)) && vlx_test_load(&chunked, chunk_path, MD_MEGABYTES(16));
	double* ref = loaded ? vlx_test_check_orbitals(utest_result, &plain, NULL) : NULL;
	double* got = ref ? vlx_test_check_orbitals(utest_result, &chunked, ref) : NULL;
	remove(chunk_path.ptr);
	EXPECT_TRUE(loaded && ref && got);
	free(ref);
	free(got);
	qm_test_free(&plain);
	qm_test_free(&chunked);
}

// The file rewritten in place after it was loaded, holding the same calculation laid out anew (C_alpha chunked): the
// provider takes the dataset's place again and reads the same orbitals
UTEST(vlx, orbitals_follow_a_rewritten_file) {
	const md_str_t plain_path = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5");
	const md_str_t copy_path  = MD_STR_LIT(VLX_PE_DIR "unittest_orbitals_rewritten.h5");
	ASSERT_TRUE(vlx_test_write_chunked_h5(plain_path, copy_path, "scf/S"));   // a copy, its orbitals as they were
	vlx_test_t t = {0};
	const bool loaded = vlx_test_load(&t, copy_path, MD_MEGABYTES(16));
	double* before = loaded ? vlx_test_check_orbitals(utest_result, &t, NULL) : NULL;
	const bool rewritten = loaded && vlx_test_write_chunked_h5(plain_path, copy_path, "scf/C_alpha");
	double* after = rewritten && before ? vlx_test_check_orbitals(utest_result, &t, before) : NULL;
	remove(copy_path.ptr);
	EXPECT_TRUE(loaded && rewritten && before && after);
	free(before);
	free(after);
	qm_test_free(&t);
}

// The NTOs are the singular value decomposition of the transition X = Z - Y, so with every pair kept
// the weights sum to |X|^2, and the hole and particle orbitals are each orthonormal in the overlap
// metric, as far as the MO coefficients are (C S C^T = 1, and the singular vectors are orthonormal).
// The first state read again after all the others is the same to the bit.
UTEST(vlx, nto_are_the_transition_decomposed) {
	vlx_test_t t = {0};
	ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/pe_induced_dipoles.h5"), MD_MEGABYTES(64)));
	const md_attribute_t* la  = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/lambda"));
	const md_attribute_t* ha  = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/hole/coefficient"));
	const md_attribute_t* pa  = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/nto/particle/coefficient"));
	const md_attribute_t* sol = qm_test_attr(&t, MD_STR_LIT("vlx/rsp/solution_matrix"));
	ASSERT_TRUE(la && ha && pa && sol);
	ASSERT_TRUE(la->data == NULL && ha->data == NULL && pa->data == NULL);

	const size_t nocc = (size_t)qm_test_scalar(&t, MD_STR_LIT("vlx/rsp/num_valence"), 0.0);
	const size_t nvir = (size_t)qm_test_scalar(&t, MD_STR_LIT("vlx/rsp/num_virtual"), 0.0);
	const uint32_t num_states = la->format.shape[0];
	const size_t lmax = la->format.shape[1], num_ao = ha->format.shape[2], len = sol->format.shape[1];
	ASSERT_TRUE(num_states > 1);
	ASSERT_EQ(lmax, nocc < nvir ? nocc : nvir); // every pair kept
	ASSERT_TRUE(len == 2 * nocc * nvir);

	size_t dim = 0;
	const double* S = qm_test_matrix(&t, MD_STR_LIT("basis/overlap"), &dim);
	ASSERT_TRUE(S != NULL);
	ASSERT_EQ(dim, num_ao);
	// the bound the orbitals they are built from meet: this file's are orthonormal to some 4e-7
	const double mo_error = qm_test_orthonormality(&t, MD_STR_LIT("orbital/alpha/coefficient"));
	ASSERT_LT(mo_error, 1.0e-5);

	double* lam   = (double*)malloc(sizeof(double) * lmax);
	double* x     = (double*)malloc(sizeof(double) * len);
	double* orb   = (double*)malloc(sizeof(double) * lmax * num_ao);
	double* first = (double*)malloc(sizeof(double) * lmax * num_ao);
	double* So    = (double*)malloc(sizeof(double) * num_ao);
	for (uint32_t s = 0; s < num_states; ++s) {
		ASSERT_EQ(md_attribute_extract_f64(lam, lmax, &t.sys.attributes, la, md_attribute_slice_1(s), md_unit_none()), lmax);
		ASSERT_EQ(md_attribute_extract_f64(x, len, &t.sys.attributes, sol, md_attribute_slice_1(s), md_unit_none()), len);
		double norm2 = 0.0, sum = 0.0;
		for (size_t k = 0; k < nocc * nvir; ++k) norm2 += (x[k] - x[nocc * nvir + k]) * (x[k] - x[nocc * nvir + k]);
		for (size_t k = 0; k < lmax; ++k) {
			sum += lam[k];
			if (k > 0) EXPECT_TRUE(lam[k] <= lam[k - 1]);
		}
		EXPECT_NEAR(norm2, sum, 1.0e-12 * norm2);

		const md_attribute_t* orbs[2] = { ha, pa };
		for (int o = 0; o < 2; ++o) {
			ASSERT_EQ(md_attribute_extract_f64(orb, lmax * num_ao, &t.sys.attributes, orbs[o], md_attribute_slice_1(s), md_unit_none()), lmax * num_ao);
			if (s == 0 && o == 0) memcpy(first, orb, sizeof(double) * lmax * num_ao);
			double worst = 0.0;
			for (size_t k = 0; k < lmax; ++k) {
				if (!(lam[k] > 0.0)) continue;
				for (size_t a = 0; a < num_ao; ++a) {
					double v = 0.0;
					for (size_t b = 0; b < num_ao; ++b) v += S[a * num_ao + b] * orb[k * num_ao + b];
					So[a] = v;
				}
				for (size_t l = 0; l < lmax; ++l) {
					if (!(lam[l] > 0.0)) continue;
					double v = 0.0;
					for (size_t a = 0; a < num_ao; ++a) v += orb[l * num_ao + a] * So[a];
					const double err = fabs(v - (k == l ? 1.0 : 0.0));
					worst = err > worst ? err : worst;
				}
			}
			EXPECT_LT(worst, 10.0 * mo_error + 1.0e-12);
		}

		// one pair's orbital alone is that row of the state's
		double one[4096];
		ASSERT_TRUE(num_ao <= MD_COUNTOF(one));
		ASSERT_EQ(md_attribute_extract_f64(one, num_ao, &t.sys.attributes, pa, md_attribute_slice_2(s, 0), md_unit_none()), num_ao);
		EXPECT_EQ(0, memcmp(one, orb, sizeof(double) * num_ao));
	}

	// the first state again, after the others
	ASSERT_EQ(md_attribute_extract_f64(orb, lmax * num_ao, &t.sys.attributes, ha, md_attribute_slice_1(0), md_unit_none()), lmax * num_ao);
	EXPECT_EQ(0, memcmp(orb, first, sizeof(double) * lmax * num_ao));

	free(lam);
	free(x);
	free(orb);
	free(first);
	free(So);
	qm_test_free(&t);
}
