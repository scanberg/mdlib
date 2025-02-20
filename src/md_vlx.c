#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>

#include <md_util.h>
#include <md_molecule.h>
#include <md_gto.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <float.h>

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999

#define HARTREE_TO_EV 27.2114079527

typedef enum {
	VLX_FLAG_CORE = 1,
	VLX_FLAG_SCF  = 2,
	VLX_FLAG_RSP  = 4,
	VLX_FLAG_ALL  = -1,
} vlx_flags_t;

/*

COMMENTS (Robin):

This file is meant to cover the VeloxChem file format.
Alot of the functionality for constructing and extracting basis functions and gaussian type orbitals is extracted from the VeloxChem source code.
https://github.com/VeloxChem/VeloxChem which is released under the LGPL-3.0 license.

*/

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
	str_t identifier;
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
	// 0 is a NULL entry, 1 = Hydrogen, 2 = Helium etc.
	struct {
		size_t count;
		basis_set_basis_t* data;
	} atom_basis;
} basis_set_t;

// New format

typedef struct md_vlx_1d_data_t {
	size_t  size;
	double* data;
} md_vlx_1d_data_t;

typedef struct md_vlx_2d_data_t {
	size_t  size[2];
	double* data;
} md_vlx_2d_data_t;

typedef struct md_vlx_orbital_t {
	md_vlx_2d_data_t coefficients;
	md_vlx_1d_data_t energy;
	md_vlx_1d_data_t occupancy;
	// In VeloxChem there are two additional fields present in the case of molecular orbitals which are D, and F
	// D may correspond to density
	// F likely corresponds to the Fock Matrix
} md_vlx_orbital_t;

typedef struct md_vlx_scf_history_t {
	size_t  number_of_iterations;
	double* density_diff;
	double* energy_diff;
	double* energy;
	double* gradient_norm;
	double* max_gradient;
} md_vlx_scf_history_t;

// Self Consistent Field
typedef struct md_vlx_scf_t {
	md_vlx_scf_type_t type;
	size_t homo_idx;
	size_t lumo_idx;

	double energy;
	dvec3_t ground_state_dipole_moment;

	md_vlx_orbital_t alpha;
	md_vlx_orbital_t beta;

	md_vlx_2d_data_t S;
	md_vlx_scf_history_t history;
} md_vlx_scf_t;

typedef struct md_vlx_rsp_t {
	size_t   number_of_excited_states;
	dvec3_t* electric_transition_dipoles;
	dvec3_t* magnetic_transition_dipoles;
	dvec3_t* velocity_transition_dipoles;
	double* rotatory_strengths;		// unit = 10^-40 cgs
	double* oscillator_strengths;
	double* absorption_ev;
	md_vlx_orbital_t* nto;
} md_vlx_rsp_t;

typedef struct md_vlx_t {
	str_t  basis_set_ident;
	str_t  dft_func_label;
	str_t  potfile_text;

	size_t number_of_atoms;
	size_t number_of_alpha_electrons;
	size_t number_of_beta_electrons;

	double molecular_charge;
	double nuclear_repulsion;
	size_t spin_multiplicity;

	// Arrays (length = number_of_atoms)
	dvec3_t* atom_coordinates;
	md_element_t* atomic_numbers;

	// Data blocks
	md_vlx_scf_t scf;
	md_vlx_rsp_t rsp;
	basis_set_t basis_set;

	struct md_allocator_i* arena;
} md_vlx_t;

static int char_to_angular_momentum_type(int c) {
	switch (c) {
	case 'S': return 0;
	case 'P': return 1;
	case 'D': return 2;
	case 'F': return 3;
	case 'G': return 4;
	default: return -1;
	}
}

static inline basis_set_basis_t* basis_set_get_atom_basis(const basis_set_t* basis_set, int atomic_number) {
	if (atomic_number < basis_set->atom_basis.count) {
		return basis_set->atom_basis.data + atomic_number;
	}
	return NULL;
}

static inline int compute_max_angular_momentum(const basis_set_t* basis_set, const md_element_t* atomic_numbers, size_t count) {
	int max_angl = 0;
	for (size_t i = 0; i < count; ++i) {
		const basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_numbers[i]);
		max_angl = MAX(max_angl, (int)atom_basis->max_type);
	}
	return max_angl;
}

#define d3  3.464101615137754587
#define f5  1.581138830084189666
#define f15 7.745966692414833770
#define f3  1.224744871391589049
#define g35 4.0 * 5.916079783099616042
#define g17 4.0 * 4.183300132670377739
#define g5  4.0 * 2.236067977499789696
#define g2  4.0 * 1.581138830084189666

static const double		S_factors[] = {1.0};
static const uint8_t    S_indices[] = {0};
static const uint8_t	S_num_fac[] = {1};

static const double		P_factors[] = {1.0, 1.0, 1.0};
static const uint8_t	P_indices[] = {1, 2, 0};
static const uint8_t	P_offsets[] = {0, 1, 2};
static const uint8_t	P_num_fac[] = {1, 1, 1};

static const double		D_factors[] = {d3, d3, -1.0, -1.0, 2.0, d3, 0.5 * d3, -0.5 * d3};
static const uint8_t    D_indices[] = {1, 4, 0, 3, 5, 2, 0, 3};
static const uint8_t    D_offsets[] = {0, 1, 2, 5, 6};
static const uint8_t	D_num_fac[] = {1, 1, 3, 1, 2};

static const double		F_factors[] = {3.0 * f5, -f5, f15, 4.0 * f3, -f3, -f3, 2.0, -3.0, -3.0, 4.0 * f3, -f3, -f3, 0.5 * f15, -0.5 * f15, f5, -3.0 * f5};
static const uint8_t    F_indices[] = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
static const uint8_t	F_offsets[] = {0, 2, 3, 6, 9, 12, 14};
static const uint8_t	F_num_fac[] = {2, 1, 3, 3, 3, 2, 2};

static const double		G_factors[] = {
	g35, -g35, 3.0 * g17, -g17, 6.0 * g5, -g5, -g5, 4.0 * g2, -3.0 * g2, -3.0 * g2,
	8.0, 3.0, 3.0, 6.0, -24.0, -24.0, 4.0 * g2, -3.0 * g2, -3.0 * g2, 3.0 * g5,
	-3.0 * g5, -0.5 * g5, 0.5 * g5,  g17,  -3.0 * g17, 0.25 * g35, 0.25 * g35, -1.50 * g35};
static const uint8_t    G_indices[] = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2, 7, 0, 10, 3};
static const uint8_t	G_offsets[] = {0, 2, 4, 7, 10, 16, 19, 23, 25};
static const uint8_t	G_num_fac[] = {2, 2, 3, 3, 6, 3, 4, 2, 3};

#undef d3
#undef f5
#undef f15
#undef f3
#undef g35
#undef g17
#undef g5 
#undef g2 

static inline int spherical_momentum_num_components(int angl) {
	return angl * 2 + 1;
}

static inline int spherical_momentum_num_factors(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph < ARRAY_SIZE(S_num_fac));
		return S_num_fac[isph];
	case 1:
		ASSERT(isph < ARRAY_SIZE(P_num_fac));
		return P_num_fac[isph];
	case 2:
		ASSERT(isph < ARRAY_SIZE(D_num_fac));
		return D_num_fac[isph];
	case 3:
		ASSERT(isph < ARRAY_SIZE(F_num_fac));
		return F_num_fac[isph];
	case 4:
		ASSERT(isph < ARRAY_SIZE(G_num_fac));
		return G_num_fac[isph];
	default:
		ASSERT(false);
		return 0;
	}
}

static inline const double* spherical_momentum_factors(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph == 0);
		return S_factors;
	case 1:
		ASSERT(isph < ARRAY_SIZE(P_offsets));
		return P_factors + P_offsets[isph];
	case 2:
		ASSERT(isph < ARRAY_SIZE(D_offsets));
		return D_factors + D_offsets[isph];
	case 3:
		ASSERT(isph < ARRAY_SIZE(F_offsets));
		return F_factors + F_offsets[isph];
	case 4:
		ASSERT(isph < ARRAY_SIZE(G_offsets));
		return G_factors + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

static inline const uint8_t* spherical_momentum_indices(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph == 0);
		return S_indices;
	case 1:
		ASSERT(isph < ARRAY_SIZE(P_offsets));
		return P_indices + P_offsets[isph];
	case 2:
		ASSERT(isph < ARRAY_SIZE(D_offsets));
		return D_indices + D_offsets[isph];
	case 3:
		ASSERT(isph < ARRAY_SIZE(F_offsets));
		return F_indices + F_offsets[isph];
	case 4:
		ASSERT(isph < ARRAY_SIZE(G_offsets));
		return G_indices + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

typedef uint8_t lmn_t[3];

// S: 0
static const lmn_t S_lmn[1] = {{0,0,0}};
// P: x y z
static const lmn_t P_lmn[3] = {{1,0,0}, {0,1,0}, {0,0,1}};
// D: xx xy xz yy yz zz
static const lmn_t D_lmn[6] = {{2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2}};
// F: xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
static const lmn_t F_lmn[10] = {{3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1}, {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}};
// G: xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz
static const lmn_t G_lmn[15] = {{4,0,0}, {3,1,0}, {3,0,1}, {2,2,0}, {2,1,1}, {2,0,2}, {1,3,0}, {1,2,1}, {1,1,2}, {1,0,3}, {0,4,0}, {0,3,1}, {0,2,2}, {0,1,3}, {0,0,4}};

static inline const lmn_t* cartesian_angular_momentum(int angl) {
	switch (angl) {
	case 0: return S_lmn;
	case 1: return P_lmn;
	case 2: return D_lmn;
	case 3: return F_lmn;
	case 4: return G_lmn;
	default: ASSERT(false); return NULL;
	}
}

typedef struct basis_func_range_t {
	int beg;
	int end;
} basis_func_range_t;

static basis_func_range_t basis_get_atomic_angl_basis_func_range(const basis_set_t* basis_set, int atomic_number, int angl) {
	basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
	basis_func_range_t range = {0};

	if (atom_basis) {
		int beg = atom_basis->basis_func_offset;
		int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
		for (int i = beg; i < end; ++i) {
			int type = basis_set->basis_func.data[i].type;
			if (type == angl) {
				range.beg = (range.end == 0) ? i : range.beg;
				range.end = i + 1;
			}
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

static inline basis_func_t get_basis_func(const basis_set_t* basis_set, int basis_func_idx) {
	basis_set_func_t func = basis_set->basis_func.data[basis_func_idx];
	return (basis_func_t) {
		.type = func.type,
		.count = func.param_count,
		.exponents = basis_set->param.exponents + func.param_offset,
		.normalization_coefficients = basis_set->param.normalization_coefficients + func.param_offset,
	};
}

static size_t compPhiAtomicOrbitals(double* out_phi, size_t phi_cap,
	const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t num_atoms,
	const basis_set_t* basis_set,
	double xp,
	double yp,
	double zp)
{
	int natoms = (int)num_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, num_atoms);

	size_t count = 0;

	// azimuthal quantum number: s,p,d,f,...
	for (int aoidx = 0, angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		const lmn_t* lmn = cartesian_angular_momentum(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// prepare Cartesian components (Maximum number of components should be 6 here for the currently supported basis set)
			double lx[8];
			double ly[8];
			double lz[8];
			int			      ncomp = spherical_momentum_num_factors(angl, isph);
			const double*	fcarts  = spherical_momentum_factors(angl, isph);
			const uint8_t*	indices = spherical_momentum_indices(angl, isph);

			for (int icomp = 0; icomp < ncomp; icomp++) {
				int cartind = indices[icomp];

				lx[icomp] = lmn[cartind][0];
				ly[icomp] = lmn[cartind][1];
				lz[icomp] = lmn[cartind][2];
			}

			// go through atoms

			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				// process coordinates
				// Conversion from Ångström to Bohr
				double rx = (xp - atom_coordinates[atomidx].x) * ANGSTROM_TO_BOHR;
				double ry = (yp - atom_coordinates[atomidx].y) * ANGSTROM_TO_BOHR;
				double rz = (zp - atom_coordinates[atomidx].z) * ANGSTROM_TO_BOHR;
				double r2 = rx*rx + ry*ry + rz*rz;

				// process atomic orbitals
				int idelem = atomic_numbers[atomidx];

				basis_func_range_t range = basis_get_atomic_angl_basis_func_range(basis_set, idelem, angl);
				for (int funcidx = range.beg; funcidx < range.end; funcidx++, aoidx++) {
					double phiao = 0.0;

					basis_func_t bf = get_basis_func(basis_set, funcidx);

					// process primitives
					for (int iprim = 0; iprim < bf.count; iprim++) {
						double expon = exp(-bf.exponents[iprim] * r2);
						double coef1 = bf.normalization_coefficients[iprim];

						// transform from Cartesian to spherical harmonics
						for (int icomp = 0; icomp < ncomp; icomp++) {
							double coef2 = coef1 * fcarts[icomp];
							double powxyz = pow(rx, lx[icomp]) * pow(ry, ly[icomp]) * pow(rz, lz[icomp]);
							phiao += coef2 * powxyz * expon;
						}
					}

					out_phi[count++] = phiao;
					if (count == phi_cap) {
						return count;
					}
				}
			}
		}
	}

	return count;
}

static size_t vlx_pgto_count(const md_vlx_t* vlx) {
	int natoms = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	size_t count = 0;

	// azimuthal quantum number: s,p,d,f,...
	for (int aoidx = 0, angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			int	ncomp = spherical_momentum_num_factors(angl, isph);
			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				int idelem = vlx->atomic_numbers[atomidx];

				// process atomic orbitals
				basis_func_range_t range = basis_get_atomic_angl_basis_func_range(&vlx->basis_set, idelem, angl);
				for (int funcidx = range.beg; funcidx < range.end; funcidx++, aoidx++) {
					// process primitives
					basis_func_t basis_func = get_basis_func(&vlx->basis_set, funcidx);
					count += basis_func.count * ncomp;
				}
			}
		}
	}

	return count;
}

static size_t extract_pgto_data(md_gto_t* pgtos, const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set, const double* mo_coeffs) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	size_t count = 0;
	size_t mo_coeff_idx = 0;

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		const lmn_t* lmn = cartesian_angular_momentum(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// prepare Cartesian components (Maximum number of components should be 6 here for the currently supported basis sets)
			int lx[8];
			int ly[8];
			int lz[8];
			int			      ncomp	= spherical_momentum_num_factors(angl, isph);
			const double*	fcarts  = spherical_momentum_factors(angl, isph);
			const uint8_t*	indices = spherical_momentum_indices(angl, isph);

			for (int icomp = 0; icomp < ncomp; icomp++) {
				int cartind = indices[icomp];

				lx[icomp] = lmn[cartind][0];
				ly[icomp] = lmn[cartind][1];
				lz[icomp] = lmn[cartind][2];
			}

			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				// process coordinates
				// Conversion from Ångström to Bohr
				float x = (float)(atom_coordinates[atomidx].x * ANGSTROM_TO_BOHR);
				float y = (float)(atom_coordinates[atomidx].y * ANGSTROM_TO_BOHR);
				float z = (float)(atom_coordinates[atomidx].z * ANGSTROM_TO_BOHR);

				int idelem = atomic_numbers[atomidx];

				// process atomic orbitals
				basis_func_range_t range = basis_get_atomic_angl_basis_func_range(basis_set, idelem, angl);
				for (int funcidx = range.beg; funcidx < range.end; funcidx++) {
					const double mo_coeff = mo_coeffs[mo_coeff_idx++];

					// process primitives
					basis_func_t basis_func = get_basis_func(basis_set, funcidx);
					ASSERT(basis_func.type == angl);
					const int        nprims = basis_func.count;
					const double* exponents = basis_func.exponents;
					const double* normcoefs = basis_func.normalization_coefficients;

					for (int iprim = 0; iprim < nprims; iprim++) {
						double alpha = exponents[iprim];
						double coef1 = normcoefs[iprim];

						// transform from Cartesian to spherical harmonics
						for (int icomp = 0; icomp < ncomp; icomp++) {
							double fcart = fcarts[icomp];
							double coeff = coef1 * fcart * mo_coeff;

							pgtos[count].x		= x;
							pgtos[count].y		= y;
							pgtos[count].z		= z;
							pgtos[count].coeff  = (float)coeff;
							pgtos[count].alpha  = (float)alpha; 
							pgtos[count].cutoff = FLT_MAX;
							pgtos[count].i		= (uint8_t)lx[icomp];
							pgtos[count].j		= (uint8_t)ly[icomp];
							pgtos[count].k		= (uint8_t)lz[icomp];
							pgtos[count].l		= (uint8_t)angl;

							count += 1;
						}
					}
				}
			}
		}
	}

	return count;
}

static double compute_overlap(basis_func_t func, size_t i, size_t j) {
	const double fab  = 1.0 / (func.exponents[i] + func.exponents[j]);
	const double fab2 = fab * fab;
	const double ovl = func.normalization_coefficients[i] * func.normalization_coefficients[j] * pow(PI * fab, 1.5);

	switch (func.type) {
	case 0: return ovl;
	case 1: return 0.5 * fab * ovl;
	case 2: return 3.0 * fab2 * ovl;
	case 3: return 7.5 * fab2 * fab * ovl;
	case 4: return 420.0 * fab2 * fab2 * ovl;
	case 5: return 1890.0 * fab2 * fab2 * fab * ovl;
	case 6: return 41580.0 * fab2 * fab2 * fab2 * ovl;
	default:
		ASSERT(false);
		return 0;
	}
}

static void rescale_basis_func(basis_func_t func) {
	const double fpi = 2.0 / PI;

	for (size_t i = 0; i < func.count; i++) {
		func.normalization_coefficients[i] *= pow(func.exponents[i] * fpi, 0.75);
	}

	if (func.type < 0 || 6 < func.type) {
		MD_LOG_DEBUG("Invalid basis function type supplied in rescaling");
		return;
	}

	static const double f_table[] = {
		0,
		2.0,
		1.15470053837925152902, // 2.0 / sqrt(3.0)
		1.03279555898864450271, // 4.0 / sqrt(15.0)
		0.19518001458970663587, // 2.0 / sqrt(105.0)
		0.13012000972647109058, // 4.0 / sqrt(945.0)
		0.03923265908909997910, // 4.0 / sqrt(10395.0)
	};

	double f = f_table[func.type];
	double e = (double)func.type * 0.5;

	for (size_t i = 0; i < func.count; i++) {
		func.normalization_coefficients[i] *= pow(f * func.exponents[i], e);
	}
}

static void normalize_basis_set(basis_set_t* basis_set) {
	for (size_t func_idx = 0; func_idx < basis_set->basis_func.count; ++func_idx) {
		basis_func_t func = get_basis_func(basis_set, (int)func_idx);
		// uncontracted basis, set expansion coeficient to 1.0
		if (func.count == 1) func.normalization_coefficients[0] = 1.0;

		// normalize primitive GBFs
		rescale_basis_func(func);

		// compute overlap
		double ovl = 0.0;
		for (size_t i = 0; i < func.count; i++) {
			ovl += compute_overlap(func, i, i);
			for (size_t j = i + 1; j < func.count; j++) {
				ovl += 2.0 * compute_overlap(func, i, j);
			}
		}

		// renormalize primitive BFs
		ovl = 1.0 / sqrt(ovl);
		for (size_t i = 0; i < func.count; i++) {
			func.normalization_coefficients[i] *= ovl;
		}
	}
}

static bool parse_basis_set(basis_set_t* basis_set, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	MEMSET(basis_set, 0, sizeof(basis_set_t));

	str_t line;
	str_t tok[4];

	// Insert null_basis element for index 0
	const basis_set_basis_t null_basis = {0};

	basis_set_basis_t* curr_atom_basis = NULL;
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (!num_tok) continue;

		if (num_tok == 2 && str_eq(tok[0], STR_LIT("@BASIS_SET"))) {
			basis_set->identifier = str_copy(tok[1], alloc);
		}
		else if (num_tok == 2 && str_eq(tok[0], STR_LIT("@ATOMBASIS"))) {
			int atomic_number = md_util_element_lookup_ignore_case(tok[1]);
			if (atomic_number == 0) {
				MD_LOG_ERROR("Unrecognized element '" STR_FMT "' in basis set", STR_ARG(tok[1]));
				return false;
			}
			basis_set_basis_t atom_basis = {
				.max_type = 0,
				.basis_func_count = 0,
				.basis_func_offset = (uint16_t)basis_set->basis_func.count,
			};

			// Grow the array and fill in slots with null_basis
			while (md_array_size(basis_set->atom_basis.data) < atomic_number) {
				md_array_push(basis_set->atom_basis.data, null_basis, alloc);
			}

			curr_atom_basis = md_array_push(basis_set->atom_basis.data, atom_basis, alloc);
			basis_set->atom_basis.count = md_array_size(basis_set->atom_basis.data);
		}
		else if (num_tok == 1 && str_eq(tok[0], STR_LIT("@END"))) {
			curr_atom_basis = NULL;
		}
		else if (num_tok == 3) {
			int type = char_to_angular_momentum_type(tok[0].ptr[0]);
			if (type == -1) {
				MD_LOG_ERROR("Unrecognized angular momentum type '" STR_FMT "' in basis set", STR_ARG(tok[0]));
				return false;
			}
			int count = (int)parse_int(tok[1]);
			if (count == 0 || count > 255) {
				MD_LOG_ERROR("Invalid number of coefficients in atom basis in basis set");
				return false;
			}

			if (curr_atom_basis == NULL) {
				MD_LOG_ERROR("No atom basis has been defined for supplied coefficients");
				return false;
			}

			// We have a new basis function for the current atom basis
			basis_set_func_t basis_func = {
				.type = (uint8_t)type,
				.param_count = (uint8_t)count,
				.param_offset = (uint16_t)basis_set->param.count,
			};
			md_array_push(basis_set->basis_func.data, basis_func, alloc);
			basis_set->basis_func.count += 1;

			for (int i = 0; i < count; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader)) {
					MD_LOG_ERROR("Failed to parse coefficients in atom basis function");
					return false;
				}
				num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
				if (num_tok != 2) {
					MD_LOG_ERROR("Unexpected number of parameters in atom basis function coefficients");
					return false;
				}
				char buf[64];
				str_copy_to_char_buf(buf, sizeof(buf), tok[0]);

				double exponent = parse_float(tok[0]);
				double coeff    = parse_float(tok[1]);
				
				md_array_push(basis_set->param.exponents, exponent, alloc);
				md_array_push(basis_set->param.normalization_coefficients, coeff, alloc);
				basis_set->param.count += 1;
			}

			curr_atom_basis->basis_func_count += 1;
			curr_atom_basis->max_type = MAX(curr_atom_basis->max_type, (uint8_t)type);
		}
	}

	return true;
}

static bool parse_vlx_geom(md_vlx_t* vlx, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	bool result = false;

	md_buffered_reader_skip_line(reader); // ====...
	md_buffered_reader_skip_line(reader); // *empty*
	md_buffered_reader_skip_line(reader); // Atom  Coordinate X  Coordinate Y  Coordinate Z
	md_buffered_reader_skip_line(reader); // *empty* 

	typedef struct {
		md_label_t sym;
		double x, y, z;
	} field_t;

	size_t temp_pos = md_temp_get_pos();
	md_allocator_i* temp_alloc = md_get_temp_allocator();
	md_array(field_t) fields = 0;

	str_t line;
	str_t tok[8];
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok == 4) {
			field_t f = {0};
			f.sym = make_label(tok[0]);
			f.x   = parse_float(tok[1]);
			f.y   = parse_float(tok[2]);
			f.z   = parse_float(tok[3]);

			md_array_push(fields, f, temp_alloc);
		} else if (num_tok == 0) {
			// Assume valid end here upon empty line
			break;
		} else {
			MD_LOG_ERROR("Unexpected number of tokens in geometry section, expected 4, got (%zu)");
			goto done;
		}
	}

	size_t count = md_array_size(fields);
	if (count == 0) {
		MD_LOG_ERROR("No atomic coordinates found");
		goto done;
	}

	// If we end up here, we may read the following lines in order
	// These are however not present in geometry optimizations
	
	if (md_buffered_reader_peek_line(&line, reader) && str_begins_with(str_trim(line), STR_LIT("Molecular charge"))) {
		// Molecular charge            : int                                                                 
		// Spin multiplicity           : int                                                               
		// Number of atoms             : int                                                               
		// Number of alpha electrons   : int                                                               
		// Number of beta  electrons   : int   
		str_t field_ident[] = {
			STR_LIT("Molecular charge            :"),                                                              
			STR_LIT("Spin multiplicity           :"),                                                            
			STR_LIT("Number of atoms             :"),                                                            
			STR_LIT("Number of alpha electrons   :"),                                                            
			STR_LIT("Number of beta  electrons   :"),
		};
		int64_t field_vals[ARRAY_SIZE(field_ident)];

		size_t loc;
		for (size_t i = 0; i < ARRAY_SIZE(field_ident); ++i) {
			str_t ident = field_ident[i];
			if (!md_buffered_reader_extract_line(&line, reader) || !str_find_str(&loc, line, ident)) {
				MD_LOG_ERROR("Failed to parse line: '"STR_FMT"'", STR_ARG(field_ident[i]));
				goto done;
			}
			str_t value_str = str_trim(str_substr(line, loc + str_len(ident), SIZE_MAX));
			field_vals[i] = parse_int(value_str);
		}

		if ((size_t)field_vals[2] != count) {
			MD_LOG_ERROR("Incorrect number of atoms parsed, expected %zu entries, parsed %zu.", (size_t)field_vals[2], count);
			goto done;
		}

		// Copy data
		vlx->molecular_charge			= (int)field_vals[0];
		vlx->spin_multiplicity			= (int)field_vals[1];
		vlx->number_of_atoms			= (size_t)field_vals[2];
		vlx->number_of_alpha_electrons	= (size_t)field_vals[3];
		vlx->number_of_beta_electrons	= (size_t)field_vals[4];
	}

	md_array_grow(vlx->atomic_numbers, count, alloc);
	md_array_grow(vlx->atom_coordinates, count, alloc);

	for (size_t i = 0; i < count; ++i) {
		md_element_t elem = md_util_element_lookup_ignore_case(LBL_TO_STR(fields[i].sym));
		if (elem == 0) {
			MD_LOG_ERROR("Unrecognized element '%s' in geometry", fields[i].sym);
			goto done;
		}
		vlx->atomic_numbers[i] = md_util_element_lookup_ignore_case(LBL_TO_STR(fields[i].sym));
		vlx->atom_coordinates[i] = (dvec3_t){ fields[i].x, fields[i].y, fields[i].z };
	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);
	return result;
}

static bool parse_vlx_basis(md_vlx_t* vlx, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	md_buffered_reader_skip_line(reader); // ====...
	md_buffered_reader_skip_line(reader); // *empty*

	str_t line;
	str_t tok[8];
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok == 2 && str_eq(tok[0], STR_LIT("Basis:"))) {
			vlx->basis_set_ident = str_copy(tok[1], alloc);
			return true;
		} else if (num_tok == 1 && str_begins_with(tok[0], STR_LIT("====="))) {
			// Parsed into next section >.<
			return false;
		}
	}
	return false;
}

static bool parse_vlx_scf(md_vlx_t* vlx, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[8];
	str_t line;

	md_buffered_reader_skip_line(reader); // =====

	int mask = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);

		if (str_begins_with(line, STR_LIT("Wave Function Model"))) {
			size_t loc = 0;
			if (str_find_char(&loc, line, ':')) {
				str_t scf_type = str_trim_beg(str_substr(line, loc+1, SIZE_MAX));
				if (str_begins_with(scf_type, STR_LIT("Spin-Restricted Open-Shell"))) {
					vlx->scf.type = MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL;
				} else if (str_begins_with(scf_type, STR_LIT("Spin-Restricted"))) {
					vlx->scf.type = MD_VLX_SCF_TYPE_RESTRICTED;
				} else if (str_begins_with(scf_type, STR_LIT("Spin-Unrestricted"))) {
					vlx->scf.type = MD_VLX_SCF_TYPE_UNRESTRICTED;
				} else {
					vlx->scf.type = MD_VLX_SCF_TYPE_UNKNOWN;
					MD_LOG_ERROR("Unexpected Wave Function Model: '"STR_FMT"'", STR_ARG(scf_type));
					return false;
				}
				mask |= 1;
			}
		} else if (str_begins_with(line, STR_LIT("Iter. |")) && str_ends_with(line, STR_LIT("| Density Change"))) {
			md_buffered_reader_skip_line(reader); // -----
			// Parse table, start tokenization
			while (md_buffered_reader_extract_line(&line, reader)) {
				size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
				if (num_tok == 6) {
					double energy_tot		= parse_float(tok[1]);
					double energy_change	= parse_float(tok[2]);
					double gradient_norm	= parse_float(tok[3]);
					double max_gradient		= parse_float(tok[4]);
					double density_change	= parse_float(tok[5]);

					md_array_push(vlx->scf.history.energy,		  energy_tot, alloc);
					md_array_push(vlx->scf.history.energy_diff,	  energy_change, alloc);
					md_array_push(vlx->scf.history.gradient_norm, gradient_norm, alloc);
					md_array_push(vlx->scf.history.max_gradient,  max_gradient, alloc);
					md_array_push(vlx->scf.history.density_diff,  density_change, alloc);
					vlx->scf.history.number_of_iterations += 1;

					mask |= 2;
				} else if (num_tok == 0) {
					// Assume valid end here upon empty line
					break;
				} else {
					MD_LOG_ERROR("Unexpected number of tokens in scf energy iteration section, expected 6, got (%zu)");
					return false;
				}
			}
		} else if (str_eq(line, STR_LIT("Ground State Dipole Moment"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_extract_line(&line, reader);
			line = str_trim_beg(line);

			if (str_begins_with(line, STR_LIT("*** Warning:"))) {
				md_buffered_reader_skip_line(reader);
				md_buffered_reader_skip_line(reader);
				md_buffered_reader_skip_line(reader);
			}

			dvec3_t vec = { 0 };
			for (int i = 0; i < 3; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 6) {
					MD_LOG_ERROR("Failed to parse SCF Ground State Dipole Moment, incomplete fields!");
					return false;
				}
				vec.elem[i] = parse_float(tok[2]);
			}
			vlx->scf.ground_state_dipole_moment.x = vec.x;
			vlx->scf.ground_state_dipole_moment.y = vec.y;
			vlx->scf.ground_state_dipole_moment.z = vec.z;
			mask |= 4;
		} else if (str_begins_with(line, STR_LIT("====="))) {
			// We've read too far and into the next section
			MD_LOG_ERROR("Failed to parse SCF section, some fields are missing");
			return false;
		}
		if (mask == 7) {
			return true;
		}
	}

	return false;
}

static bool parse_vlx_rsp_dipole_moments(dvec3_t* moments, size_t num_excited_states, md_buffered_reader_t* reader) {
	str_t tok[8];
	str_t line;

	for (size_t i = 0; i < num_excited_states; ++i) {
		if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 6) {
			MD_LOG_ERROR("Unexpected number of tokens in rsp dipole moment table");
			return false;
		}
		
		str_t first = str_join(tok[0], tok[1]);
		if (!str_eq(first, STR_LIT("Excited State"))) {
			MD_LOG_ERROR("Unexpected entry in rsp dipole moment table, expected 'Excited State', got '"STR_FMT"'", first);
			return false;
		}

		str_t ident = tok[2];
		if (ident.len > 0 && ident.ptr[ident.len-1] == ':') {
			ident.len -= 1;
		}

		moments[i].x = parse_float(tok[3]);
		moments[i].y = parse_float(tok[4]);
		moments[i].z = parse_float(tok[5]);
	}

	return true;
}

static bool parse_vlx_rsp(md_vlx_t* vlx, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[16];
	str_t line;

	md_buffered_reader_skip_line(reader); // =====

	int mask = 0;
	while ( md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);
		if (str_empty(line)) continue;

		if (str_begins_with(line, STR_LIT("Number of States")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			vlx->rsp.number_of_excited_states = parse_int(tok[4]);
			mask |= 1;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole length, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(vlx->rsp.electric_transition_dipoles, vlx->rsp.number_of_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(vlx->rsp.electric_transition_dipoles, vlx->rsp.number_of_excited_states, reader)) {
				return false;
			}
			mask |= 2;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole velocity, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(vlx->rsp.velocity_transition_dipoles, vlx->rsp.number_of_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(vlx->rsp.velocity_transition_dipoles, vlx->rsp.number_of_excited_states, reader)) {
				return false;
			}
			mask |= 4;
		} else if (str_eq(line, STR_LIT("Magnetic Transition Dipole Moments (a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(vlx->rsp.magnetic_transition_dipoles, vlx->rsp.number_of_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(vlx->rsp.magnetic_transition_dipoles, vlx->rsp.number_of_excited_states, reader)) {
				return false;
			}
			mask |= 8;
		}
		else if (str_eq(line, STR_LIT("One-Photon Absorption"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_array_resize(vlx->rsp.oscillator_strengths, vlx->rsp.number_of_excited_states, alloc);
			md_array_resize(vlx->rsp.absorption_ev, vlx->rsp.number_of_excited_states, alloc);
			for (size_t i = 0; i < vlx->rsp.number_of_excited_states; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 9) {
					MD_LOG_ERROR("Unexpected number of tokens in entry when parsing One-Photon Absorption");
					return false;
				}
				vlx->rsp.absorption_ev[i] = parse_float(tok[5]);
				vlx->rsp.oscillator_strengths[i] = parse_float(tok[8]);
			}
			mask |= 16;
		} else if (str_eq(line, STR_LIT("Electronic Circular Dichroism"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_array_resize(vlx->rsp.rotatory_strengths, vlx->rsp.number_of_excited_states, alloc);
			for (size_t i = 0; i < vlx->rsp.number_of_excited_states; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 9) {
					MD_LOG_ERROR("Unexpected number of tokens in entry when parsing Electronic Circular Dichroism");
					return false;
				}
				vlx->rsp.rotatory_strengths[i] = parse_float(tok[6]);
			}
			mask |= 32;
		} else if (str_begins_with(line, STR_LIT("===="))) {
			// Parsed into next section
			return false;
		}
		if (mask == 63) {
			return true;
		}
	}
	return false;
}

static bool vlx_parse_out(md_vlx_t* vlx, md_buffered_reader_t* reader) {
	str_t line;
	while (md_buffered_reader_extract_line(&line, reader)) {
		str_t str = str_trim(line);
		if (str_eq(str, STR_LIT("Molecular Geometry (Angstroms)"))) {
			if (!parse_vlx_geom(vlx, reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse geometry");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Molecular Basis (Atomic Basis)"))) {
			if (!parse_vlx_basis(vlx, reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse basis");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Self Consistent Field Driver Setup"))) {
			if (!parse_vlx_scf(vlx, reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse SCF section");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Linear Response EigenSolver Setup"))) {
			if (!parse_vlx_rsp(vlx, reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse RSP section");
				return false;
			}
		}
	}

	return true;
}

static bool h5_read_scalar(void* buf, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
	if (exists == 0) {
		return false;
	}

	// Open the dataset containing the double value
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool result = false;

	// Read the dataset into the 'value' variable
	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if (status != 0) {
		MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
		goto done;
	}

	result = true;
done:
	H5Dclose(dataset_id);

	return result;
}

static bool h5_read_str(str_t* str, hid_t file_id, const char* field_name, md_allocator_i* alloc) {
	bool result = false;

	htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
	if (exists == 0) {
		return false;
	}

	// Open the dataset
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	// Get the datatype and space
	hid_t datatype_id = H5Dget_type(dataset_id);  // Get datatype
	hid_t space_id = H5Dget_space(dataset_id);    // Get dataspace

	// Determine size of string (assume variable-length string)
	size_t size = H5Tget_size(datatype_id);
	str_t data = str_alloc(size, alloc);
	herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)data.ptr);
	if (status != 0) {
		MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
		str_free(data, alloc);
		goto done;
	}
	*str = data;

	result = true;
done:

	// Close HDF5 resources
	H5Tclose(datatype_id);
	H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static size_t h5_read_cstr(char* out_str, size_t str_cap, hid_t file_id, const char* field_name) {
	htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
	if (exists == 0) {
		return 0;
	}

	// Open the dataset
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return 0;
	}

	size_t result = 0;

	// Get the datatype and space
	hid_t datatype_id = H5Dget_type(dataset_id);  // Get datatype

	// Determine size of string (assume variable-length string)
	size_t size = H5Tget_size(datatype_id);
	if (size < str_cap) {
		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_str);
		if (status != 0) {
			MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
			goto done;
		}
		result = size;
	}
done:

	// Close HDF5 resources
	H5Tclose(datatype_id);
	H5Dclose(dataset_id);

	return result;
}

static int h5_read_dataset_dims(size_t* dims, int max_dims, hid_t file_id, const char* field_name) {
	htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
	if (exists == 0) {
		return false;
	}

	// Open the dataset
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return -1;
	}

	// Get the datatype and space
	hid_t space_id = H5Dget_space(dataset_id);    // Get dataspace

	// Determine size of string (assume variable-length string)
	int ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim < 0) {
		MD_LOG_ERROR("Failed to get number of dimensions for H5 dataset: '%s'", field_name);
		goto done;
	}

	if (ndim > max_dims) {
		MD_LOG_ERROR("Too many dimensions in data");
		goto done;
	}

	ndim = H5Sget_simple_extent_dims(space_id, dims, 0);

done:
	H5Sclose(space_id);
	H5Dclose(dataset_id);

	return ndim;
}

static bool h5_read_dataset_data(void* out_data, const size_t dims[], int num_dim, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	ASSERT(out_data);
	ASSERT(dims);

	htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
	if (exists == 0) {
		return false;
	}

	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool result = false;
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 space");
		goto done;
	}

	int ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim != num_dim) {
		MD_LOG_ERROR("Unexpected number of dimensions when reading dataset, got %i, expected %i", ndim, num_dim);
		goto done;
	}

	if (ndim > 8) {
		MD_LOG_ERROR("Too many dimensions in data");
		goto done;
	}

	size_t temp_dims[8];

	ndim = H5Sget_simple_extent_dims(space_id, temp_dims, 0);
	ASSERT(ndim == num_dim);

	for (int i = 0; i < ndim; ++i) {
		if (temp_dims[i] != dims[i]) {
			MD_LOG_ERROR("Incompatible dims");
			goto done;
		}
	}

	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_data);

	if (status != 0) {
		MD_LOG_ERROR("An error occured when reading H5 data");
		goto done;
	}

	result = true;
done:
	H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

// Data extraction procedures

static bool h5_read_scf_data(md_vlx_t* vlx, hid_t handle) {
	char scf_type[64] = {0};
	if (!h5_read_cstr(scf_type, sizeof(scf_type), handle, "scf_type")) {
		return false;
	}

	if (str_eq_cstr(STR_LIT("restricted"), scf_type)) {
		vlx->scf.type = MD_VLX_SCF_TYPE_RESTRICTED;
	} else if (str_eq_cstr(STR_LIT("restricted_openshell"), scf_type)) {
		vlx->scf.type = MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL;
	} else if (str_eq_cstr(STR_LIT("unrestricted"), scf_type)) {
		vlx->scf.type = MD_VLX_SCF_TYPE_UNRESTRICTED;
	} else {
		vlx->scf.type = MD_VLX_SCF_TYPE_UNKNOWN;
		MD_LOG_ERROR("Unrecognized scf type present in h5 scf section: '%s'", scf_type);
		return false;
	}

	if (!h5_read_scalar(&vlx->scf.energy, handle, H5T_NATIVE_DOUBLE, "scf_energy")) {
		return false;
	}

	size_t dim[2];
	h5_read_dataset_dims(dim, 2, handle, "C_alpha");

	md_array_resize(vlx->scf.alpha.coefficients.data, dim[0] * dim[1], vlx->arena);
	MEMCPY(vlx->scf.alpha.coefficients.size, dim, sizeof(dim));

	md_array_resize(vlx->scf.alpha.energy.data, dim[1], vlx->arena);
	vlx->scf.alpha.energy.size = dim[1];

	md_array_resize(vlx->scf.alpha.occupancy.data, dim[1], vlx->arena);
	vlx->scf.alpha.occupancy.size = dim[1];

	// Extract alpha data
	if (!h5_read_dataset_data(vlx->scf.alpha.coefficients.data, vlx->scf.alpha.coefficients.size, 2, handle, H5T_NATIVE_DOUBLE, "C_alpha")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.alpha.energy.data, &vlx->scf.alpha.energy.size, 1, handle, H5T_NATIVE_DOUBLE, "E_alpha")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.alpha.occupancy.data, &vlx->scf.alpha.occupancy.size, 1, handle, H5T_NATIVE_DOUBLE, "occ_alpha")) {
		return false;
	}

	if (vlx->scf.type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
		md_array_resize(vlx->scf.beta.coefficients.data, dim[0] * dim[1], vlx->arena);
		MEMCPY(vlx->scf.beta.coefficients.size, dim, sizeof(dim));

		md_array_resize(vlx->scf.beta.energy.data, dim[1], vlx->arena);
		vlx->scf.beta.energy.size = dim[1];

		md_array_resize(vlx->scf.beta.occupancy.data, dim[1], vlx->arena);
		vlx->scf.beta.occupancy.size = dim[1];

		// Extract beta data
		if (!h5_read_dataset_data(vlx->scf.beta.coefficients.data, vlx->scf.beta.coefficients.size, 2, handle, H5T_NATIVE_DOUBLE, "C_beta")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.beta.energy.data, &vlx->scf.beta.energy.size, 1, handle, H5T_NATIVE_DOUBLE, "E_beta")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, &vlx->scf.beta.occupancy.size, 1, handle, H5T_NATIVE_DOUBLE, "occ_beta")) {
			return false;
		}
	} else {
		// Shallow copy fields from Alpha
		MEMCPY(&vlx->scf.beta, &vlx->scf.alpha, sizeof(md_vlx_orbital_t));

		if (vlx->scf.type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
			vlx->scf.beta.occupancy.data = 0;
			md_array_resize(vlx->scf.beta.occupancy.data, vlx->scf.beta.occupancy.size, vlx->arena);
			if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, &vlx->scf.beta.occupancy.size, 1, handle, H5T_NATIVE_DOUBLE, "occ_beta")) {
				return false;
			}
		}
	}

	md_array_resize(vlx->scf.S.data, dim[0] * dim[1], vlx->arena);
	MEMCPY(vlx->scf.S.size, dim, sizeof(dim));

	if (!h5_read_dataset_data(vlx->scf.S.data, vlx->scf.S.size, 2, handle, H5T_NATIVE_DOUBLE, "S")) {
		return false;
	}

	// The ground state dipole moment is not present in all versions
	size_t dipole_dim[1] = { 3 };
	if (!h5_read_dataset_data(&vlx->scf.ground_state_dipole_moment, dipole_dim, 1, handle, H5T_NATIVE_DOUBLE, "dipole_moment")) {
		//return false;
	}

	size_t scf_hist_iter = 0;
	if (!h5_read_dataset_dims(&scf_hist_iter, 1, handle, "scf_history_diff_density")) {
		return false;
	}

	if (scf_hist_iter > 0) {
		vlx->scf.history.number_of_iterations = scf_hist_iter;
		md_array_resize(vlx->scf.history.density_diff, scf_hist_iter, vlx->arena);
		md_array_resize(vlx->scf.history.energy, scf_hist_iter, vlx->arena);
		md_array_resize(vlx->scf.history.energy_diff, scf_hist_iter, vlx->arena);
		md_array_resize(vlx->scf.history.gradient_norm, scf_hist_iter, vlx->arena);
		md_array_resize(vlx->scf.history.max_gradient, scf_hist_iter, vlx->arena);
	}

	if (!h5_read_dataset_data(vlx->scf.history.density_diff, &scf_hist_iter, 1, handle, H5T_NATIVE_DOUBLE, "scf_history_diff_density")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.energy_diff, &scf_hist_iter, 1, handle, H5T_NATIVE_DOUBLE, "scf_history_diff_energy")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.energy, &scf_hist_iter, 1, handle, H5T_NATIVE_DOUBLE, "scf_history_energy")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.gradient_norm, &scf_hist_iter, 1, handle, H5T_NATIVE_DOUBLE, "scf_history_gradient_norm")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.max_gradient, &scf_hist_iter, 1, handle, H5T_NATIVE_DOUBLE, "scf_history_max_gradient")) {
		return false;
	}

	return true;
}

static bool h5_read_rsp_data(md_vlx_t* vlx, hid_t handle) {
	if (!h5_read_scalar(&vlx->rsp.number_of_excited_states, handle, H5T_NATIVE_HSIZE, "number_of_states")) {
		return false;
	}

	if (vlx->rsp.number_of_excited_states == 0) {
		// Wierd
		MD_LOG_ERROR("No excited states listed in response section of VeloxChem h5 file");
		return false;
	}

	// Allocate data
	md_array_resize(vlx->rsp.nto, vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.nto, 0, vlx->rsp.number_of_excited_states * sizeof(md_vlx_orbital_t));

	md_array_resize(vlx->rsp.electric_transition_dipoles, vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.electric_transition_dipoles, 0, vlx->rsp.number_of_excited_states * sizeof(dvec3_t));

	md_array_resize(vlx->rsp.magnetic_transition_dipoles, vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.magnetic_transition_dipoles, 0, vlx->rsp.number_of_excited_states * sizeof(dvec3_t));

	md_array_resize(vlx->rsp.velocity_transition_dipoles, vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.velocity_transition_dipoles, 0, vlx->rsp.number_of_excited_states * sizeof(dvec3_t));

	md_array_resize(vlx->rsp.absorption_ev,	vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.absorption_ev, 0, vlx->rsp.number_of_excited_states * sizeof(double));

	md_array_resize(vlx->rsp.oscillator_strengths, vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.oscillator_strengths, 0, vlx->rsp.number_of_excited_states * sizeof(double));

	md_array_resize(vlx->rsp.rotatory_strengths, vlx->rsp.number_of_excited_states, vlx->arena);
	MEMSET(vlx->rsp.rotatory_strengths, 0, vlx->rsp.number_of_excited_states * sizeof(double));

	// NTO data
	char buf[32];
	for (size_t i = 0; i < vlx->rsp.number_of_excited_states; ++i) {
		int idx = (int)(i + 1);

		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_orbitals", idx);

		uint64_t dim[2];
		if (!h5_read_dataset_dims(dim, 2, handle, buf)) {
			return false;
		}
		if (dim[0] == 0 || dim[1] == 0) {
			MD_LOG_ERROR("Invalid dimensions in NTO orbitals");
			return false;
		}

		md_array_resize(vlx->rsp.nto[i].coefficients.data, dim[0] * dim[1], vlx->arena);
		MEMSET(vlx->rsp.nto[i].coefficients.data, 0, md_array_bytes(vlx->rsp.nto[i].coefficients.data));
		MEMCPY(vlx->rsp.nto[i].coefficients.size, dim, sizeof(dim));

		md_array_resize(vlx->rsp.nto[i].energy.data, dim[1], vlx->arena);
		MEMSET(vlx->rsp.nto[i].energy.data, 0, md_array_bytes(vlx->rsp.nto[i].energy.data));
		vlx->rsp.nto[i].energy.size = dim[1];

		md_array_resize(vlx->rsp.nto[i].occupancy.data, dim[1], vlx->arena);
		MEMSET(vlx->rsp.nto[i].occupancy.data, 0, md_array_bytes(vlx->rsp.nto[i].occupancy.data));
		vlx->rsp.nto[i].occupancy.size = dim[1];

		if (!h5_read_dataset_data(vlx->rsp.nto[i].coefficients.data, vlx->rsp.nto[i].coefficients.size, 2, handle, H5T_NATIVE_DOUBLE, buf)) {
			return false;
		}

		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_occupations", idx);
		if (!h5_read_dataset_data(vlx->rsp.nto[i].occupancy.data, &vlx->rsp.nto[i].occupancy.size, 1, handle, H5T_NATIVE_DOUBLE, buf)) {
			return false;
		}

		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_energies", idx);
		if (!h5_read_dataset_data(vlx->rsp.nto[i].energy.data, &vlx->rsp.nto[i].energy.size, 1, handle, H5T_NATIVE_DOUBLE, buf)) {
			return false;
		}


	}

	// Dipoles
	size_t dipole_dim[2] = { vlx->rsp.number_of_excited_states, 3 };
	if (!h5_read_dataset_data(vlx->rsp.electric_transition_dipoles, dipole_dim, 2, handle, H5T_NATIVE_DOUBLE, "electric_transition_dipoles")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->rsp.magnetic_transition_dipoles, dipole_dim, 2, handle, H5T_NATIVE_DOUBLE, "magnetic_transition_dipoles")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->rsp.velocity_transition_dipoles, dipole_dim, 2, handle, H5T_NATIVE_DOUBLE, "velocity_transition_dipoles")) {
		return false;
	}

	// Abs, rot and osc
	size_t dim[1] = { vlx->rsp.number_of_excited_states };
	if (!h5_read_dataset_data(vlx->rsp.absorption_ev, dim, 1, handle, H5T_NATIVE_DOUBLE, "eigenvalues")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->rsp.oscillator_strengths, dim, 1, handle, H5T_NATIVE_DOUBLE, "oscillator_strengths")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->rsp.rotatory_strengths, dim, 1, handle, H5T_NATIVE_DOUBLE, "rotatory_strengths")) {
		return false;
	}

	// Convert Atomic units (Hartree) to eV
	if (vlx->rsp.absorption_ev) {
		for (size_t i = 0; i < vlx->rsp.number_of_excited_states; ++i) {
			vlx->rsp.absorption_ev[i] *= HARTREE_TO_EV;
		}
	}

	return true;
}

static bool h5_read_core_data(md_vlx_t* vlx, hid_t handle) {
	ASSERT(vlx);

	if (!h5_read_str(&vlx->basis_set_ident, handle, "basis_set", vlx->arena)) {
		return false;
	}

	if (!h5_read_str(&vlx->dft_func_label, handle, "dft_func_label", vlx->arena)) {
		return false;
	}

	if (!h5_read_scalar(&vlx->molecular_charge, handle, H5T_NATIVE_DOUBLE, "molecular_charge")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->nuclear_repulsion, handle, H5T_NATIVE_DOUBLE, "nuclear_repulsion")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_alpha_electrons, handle, H5T_NATIVE_INT64, "number_of_alpha_electrons")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_atoms, handle, H5T_NATIVE_INT64, "number_of_atoms")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_beta_electrons, handle, H5T_NATIVE_INT64, "number_of_beta_electrons")) {
		return false;
	}

	if (!h5_read_str(&vlx->potfile_text, handle, "potfile_text", vlx->arena)) {
		return false;
	}

	if (!h5_read_scalar(&vlx->spin_multiplicity, handle, H5T_NATIVE_INT64, "spin_multiplicity")) {
		return false;
	}

	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("Number of atoms is zero");
		return false;
	}

	md_array_resize(vlx->atom_coordinates, vlx->number_of_atoms, vlx->arena);
	MEMSET(vlx->atom_coordinates, 0, md_array_bytes(vlx->atom_coordinates));
	size_t coord_dim[2] = { vlx->number_of_atoms, 3 };
	if (!h5_read_dataset_data(vlx->atom_coordinates, coord_dim, 2, handle, H5T_NATIVE_DOUBLE, "atom_coordinates")) {
		return false;
	}

	md_array_resize(vlx->atomic_numbers, vlx->number_of_atoms, vlx->arena);
	MEMSET(vlx->atomic_numbers, 0, md_array_bytes(vlx->atomic_numbers));
	size_t atomic_numbers_dim[1] = { vlx->number_of_atoms };
	if (!h5_read_dataset_data(vlx->atomic_numbers, atomic_numbers_dim, 1, handle, H5T_NATIVE_UINT8, "nuclear_charges")) {
		return false;
	}

	// Convert Atomic units to Ångström
	if (vlx->atom_coordinates) {
		for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
			vlx->atom_coordinates[i].x *= BOHR_TO_ANGSTROM;
			vlx->atom_coordinates[i].y *= BOHR_TO_ANGSTROM;
			vlx->atom_coordinates[i].z *= BOHR_TO_ANGSTROM;
		}
	}

	return true;
}

static bool vlx_read_scf_results(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	ASSERT(vlx);

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	bool result = false;

	if (flags & VLX_FLAG_CORE) {
		if (!h5_read_core_data(vlx, file_id)) {
			goto done;
		}
	}
	
	if (flags & VLX_FLAG_SCF) {
		if (!h5_read_scf_data(vlx, file_id)) {
			goto done;
		}
	}

	result = true;
done:
	H5Fclose(file_id);

	return result;
}

// This is the newest version of the file format where everything is contained within a single h5 file
static bool vlx_read_h5_file(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	ASSERT(vlx);

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	bool result = false;

	if (flags & VLX_FLAG_CORE) {
		if (!h5_read_core_data(vlx, file_id)) {
			goto done;
		}
	}

	// SCF
	if (flags & VLX_FLAG_SCF) {
		if (H5Lexists(file_id, "scf", H5P_DEFAULT) > 0) {
			hid_t scf_id = H5Gopen(file_id, "scf", H5P_DEFAULT);
			if (scf_id != H5I_INVALID_HID) {
				result = h5_read_scf_data(vlx, scf_id);
				H5Gclose(scf_id);
				if (!result) goto done;
			}
		}
	}

	// RSP
	if (flags & VLX_FLAG_RSP) {
		if (H5Lexists(file_id, "rsp", H5P_DEFAULT) > 0) {
			hid_t rsp_id = H5Gopen(file_id, "rsp", H5P_DEFAULT);
			if (rsp_id != H5I_INVALID_HID) {
				result = h5_read_rsp_data(vlx, rsp_id);
				H5Gclose(rsp_id);
				if (!result) goto done;
			}
		}
	}

	result = true;
done:
	H5Fclose(file_id);

	return result;
}

static bool vlx_parse_out_file(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	size_t temp_pos = md_temp_get_pos();
	size_t cap = KILOBYTES(16);
	char* buf = md_temp_push(cap);

	bool result = false;

	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
	str_t line;
	while (md_buffered_reader_extract_line(&line, &reader)) {
		str_t str = str_trim(line);
		if (flags & VLX_FLAG_CORE) {
			if (str_eq(str, STR_LIT("Molecular Geometry (Angstroms)"))) {
				if (!parse_vlx_geom(vlx, &reader, vlx->arena)) {
					MD_LOG_ERROR("Failed to parse geometry");
					goto done;
				}
			} else if (str_eq(str, STR_LIT("Molecular Basis (Atomic Basis)"))) {
				if (!parse_vlx_basis(vlx, &reader, vlx->arena)) {
					MD_LOG_ERROR("Failed to parse basis");
					goto done;
				}
			}
		}
		if (flags & VLX_FLAG_SCF) {
			if(str_eq(str, STR_LIT("Self Consistent Field Driver Setup"))) {
				if (!parse_vlx_scf(vlx, &reader, vlx->arena)) {
					MD_LOG_ERROR("Failed to parse SCF section");
					goto done;
				}
			}
		}
		if (flags & VLX_FLAG_RSP) {
			if (str_eq(str, STR_LIT("Linear Response EigenSolver Setup"))) {
				if (!parse_vlx_rsp(vlx, &reader, vlx->arena)) {
					MD_LOG_ERROR("Failed to parse RSP section");
					goto done;
				}
			}
		}
	}
	md_file_close(file);

	str_t base_file = {0};
	extract_file_path_without_ext(&base_file, filename);

	md_strb_t sb = md_strb_create(md_get_temp_allocator());

	if (flags & VLX_FLAG_SCF) {
		// Attempt to read scf data
		md_strb_reset(&sb);
		md_strb_fmt(&sb, STR_FMT ".scf.h5", STR_ARG(base_file));
		if (md_path_is_valid(md_strb_to_str(sb))) {
			// Open an existing file
			hid_t file_id = H5Fopen(md_strb_to_cstr(sb), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id == H5I_INVALID_HID) {
				MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(md_strb_to_str(sb)));
				goto done;
			}
			size_t dim[2];
			if (!h5_read_dataset_dims(dim, 2, file_id, "alpha_orbitals")) {
				goto done;
			}

			md_array_resize(vlx->scf.alpha.coefficients.data, dim[0] * dim[1], vlx->arena);
			MEMCPY(vlx->scf.alpha.coefficients.size, dim, sizeof(dim));

			md_array_resize(vlx->scf.alpha.energy.data, dim[1], vlx->arena);
			vlx->scf.alpha.energy.size = dim[1];

			md_array_resize(vlx->scf.alpha.occupancy.data, dim[1], vlx->arena);
			vlx->scf.alpha.occupancy.size = dim[1];

			if (!h5_read_dataset_data(vlx->scf.alpha.coefficients.data, vlx->scf.alpha.coefficients.size, 2, file_id, H5T_NATIVE_DOUBLE, "alpha_orbitals")) {
				goto done;
			}
			if (!h5_read_dataset_data(vlx->scf.alpha.energy.data, &vlx->scf.alpha.energy.size, 1, file_id, H5T_NATIVE_DOUBLE, "alpha_energies")) {
				goto done;
			}
			if (!h5_read_dataset_data(vlx->scf.alpha.occupancy.data, &vlx->scf.alpha.occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "alpha_occupations")) {
				goto done;
			}

			if (vlx->scf.type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
				md_array_resize(vlx->scf.beta.coefficients.data, dim[0] * dim[1], vlx->arena);
				MEMCPY(vlx->scf.beta.coefficients.size, dim, sizeof(dim));

				md_array_resize(vlx->scf.beta.energy.data, dim[1], vlx->arena);
				vlx->scf.beta.energy.size = dim[1];

				md_array_resize(vlx->scf.beta.occupancy.data, dim[1], vlx->arena);
				vlx->scf.beta.occupancy.size = dim[1];

				// Extract beta data
				if (!h5_read_dataset_data(vlx->scf.beta.coefficients.data, vlx->scf.beta.coefficients.size, 2, file_id, H5T_NATIVE_DOUBLE, "beta_orbitals")) {
					goto done;
				}
				if (!h5_read_dataset_data(vlx->scf.beta.energy.data, &vlx->scf.beta.energy.size, 1, file_id, H5T_NATIVE_DOUBLE, "beta_energies")) {
					goto done;
				}
				if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, &vlx->scf.beta.occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "beta_occupations")) {
					goto done;
				}
			} else {
				// Shallow copy fields from Alpha
				MEMCPY(&vlx->scf.beta, &vlx->scf.alpha, sizeof(md_vlx_orbital_t));

				if (vlx->scf.type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
					vlx->scf.beta.occupancy.data = 0;
					md_array_resize(vlx->scf.beta.occupancy.data, vlx->scf.beta.occupancy.size, vlx->arena);
					if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, &vlx->scf.beta.occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "beta_occupations")) {
						goto done;
					}
				}
			}
		}
	}

	if (flags & VLX_FLAG_RSP) {
		for (int i = 0; i < (int)vlx->rsp.number_of_excited_states; ++i) {
			md_strb_reset(&sb);
			md_strb_fmt(&sb, STR_FMT "_S%i_NTO.h5", STR_ARG(base_file), i + 1);
			if (md_path_is_valid(md_strb_to_str(sb))) {
				// Open an existing file
				hid_t file_id = H5Fopen(md_strb_to_cstr(sb), H5F_ACC_RDONLY, H5P_DEFAULT);
				if (file_id == H5I_INVALID_HID) {
					MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(md_strb_to_str(sb)));
					goto done;
				}

				uint64_t dim[2];
				if (!h5_read_dataset_dims(dim, 2, file_id, "alpha_orbitals")) {
					goto done;
				}
				if (dim[0] == 0 || dim[1] == 0) {
					MD_LOG_ERROR("Invalid dimensions in NTO orbitals");
					goto done;
				}

				md_vlx_orbital_t* nto = md_array_push(vlx->rsp.nto, (md_vlx_orbital_t){0}, vlx->arena);

				md_array_resize(nto->coefficients.data, dim[0] * dim[1], vlx->arena);
				MEMCPY(nto->coefficients.size, dim, sizeof(dim));

				md_array_resize(nto->occupancy.data, dim[1], vlx->arena);
				nto->occupancy.size = dim[1];

				if (!h5_read_dataset_data(nto->coefficients.data, nto->coefficients.size, 2, file_id, H5T_NATIVE_DOUBLE, "alpha_orbitals")) {
					goto done;
				}
				if (!h5_read_dataset_data(nto->occupancy.data, &nto->occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "alpha_occupations")) {
					goto done;
				}
			} else {
				MD_LOG_INFO("The veloxchem object specified %zu excited states, but the matching NTOs could not be found.", vlx->rsp.number_of_excited_states);
				if (vlx->rsp.nto) {
					md_array_free(vlx->rsp.nto, vlx->arena);
					vlx->rsp.nto = NULL;
					break;
				}
			}
		}
	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);
	return result;
}

// Internal version to control what portions to load
static bool vlx_parse_file(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	size_t temp_pos = md_temp_get_pos();

	bool result = false;

	if (str_ends_with(filename, STR_LIT(".out"))) {
		if (!vlx_parse_out_file(vlx, filename, flags)) {
			goto done;
		}
	} else if (str_ends_with(filename, STR_LIT(".scf.results.h5"))) {
		if (!vlx_read_scf_results(vlx, filename, flags)) {
			goto done;
		}
	} else if (str_ends_with(filename, STR_LIT(".h5"))) {
		if (!vlx_read_h5_file(vlx, filename, flags)) {
			goto done;
		}
	} else {
		MD_LOG_DEBUG("Unsupported file format");
		goto done;
	}

	if (!str_empty(vlx->basis_set_ident)) {
		size_t cap = KILOBYTES(16);
		char*  buf = md_temp_push(cap);
		md_strb_t sb = md_strb_create(md_get_temp_allocator());

		md_strb_fmt(&sb, "%s/" STR_FMT, MD_VLX_BASIS_FOLDER, STR_ARG(vlx->basis_set_ident));
		md_file_o* basis_file = md_file_open(md_strb_to_str(sb), MD_FILE_READ | MD_FILE_BINARY);
		if (basis_file) {
			md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, basis_file);
			bool parse_result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
			md_file_close(basis_file);
			if (!parse_result) {
				MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
				goto done;
			}
			normalize_basis_set(&vlx->basis_set);
		} else {
			// Attempt to read basis set file from same folder as file
			str_t folder = { 0 };
			if (!extract_folder_path(&folder, filename)) {
				MD_LOG_ERROR("An error occured when extracting the path to supplied file");
				return false;
			}
			md_strb_reset(&sb);
			md_strb_push_str(&sb, folder);
			md_strb_push_str(&sb, vlx->basis_set_ident);
			basis_file = md_file_open(md_strb_to_str(sb), MD_FILE_READ | MD_FILE_BINARY);
			if (basis_file) {
				md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, basis_file);
				bool parse_result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
				md_file_close(basis_file);
				if (!parse_result) {
					MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
					goto done;
				}
				normalize_basis_set(&vlx->basis_set);
			}
			else {
				MD_LOG_ERROR("Could not find basis file corresponding to identifier: '"STR_FMT"'", STR_ARG(vlx->basis_set_ident));
				goto done;
			}
		}
	}

	// Identify homo and lumo
	if (vlx->scf.alpha.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.alpha.occupancy.size; ++i) {
			if (vlx->scf.alpha.occupancy.data[i] == 0.0) {
				vlx->scf.homo_idx = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.lumo_idx = i;
				break;
			}
		}
	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);

	return result;
}

// Extract Natural Transition Orbitals PGTOs
size_t md_vlx_nto_gto_count(const md_vlx_t* vlx) {
	return vlx_pgto_count(vlx);
}

static inline void extract_row(double* dst, const md_vlx_2d_data_t* data, size_t row_idx) {
	ASSERT(dst);
	ASSERT(data);
	ASSERT(row_idx < data->size[0]);

	size_t num_cols = data->size[1];
	for (size_t i = 0; i < num_cols; ++i) {
		dst[i] = data->data[row_idx * num_cols + i];
	}
}

static inline void extract_col(double* dst, const md_vlx_2d_data_t* data, size_t col_idx) {
	ASSERT(dst);
	ASSERT(data);
	ASSERT(col_idx < data->size[1]);

	for (size_t i = 0; i < data->size[0]; ++i) {
		dst[i] = data->data[i * data->size[1] + col_idx];
	}
}

static inline size_t number_of_molecular_orbitals(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline size_t number_of_mo_coefficients(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline void extract_mo_coefficients(double* out_coeff, const md_vlx_orbital_t* orb, size_t mo_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);

	ASSERT(mo_idx < number_of_mo_coefficients(orb));

	extract_col(out_coeff, &orb->coefficients, mo_idx);
}

static inline size_t number_of_atomic_orbitals(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[0];
}

static inline size_t number_of_ao_coefficients(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[0];
}

static inline void extract_ao_coefficients(double* out_coeff, const md_vlx_orbital_t* orb, size_t ao_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);
	ASSERT(ao_idx < orb->coefficients.size[0]);

	ASSERT(ao_idx < number_of_ao_coefficients(orb));

	extract_row(out_coeff, &orb->coefficients, ao_idx);
}

bool md_vlx_nto_gto_extract(md_gto_t* pgtos, const md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type) {
	ASSERT(pgtos);
	ASSERT(vlx);

	if (nto_idx >= vlx->rsp.number_of_excited_states) {
		MD_LOG_ERROR("Invalid nto index!");
		return false;
	}

	if (!vlx->rsp.nto) {
		MD_LOG_ERROR("Veloxchem data is missing NTO data");
		return false;
	}

	// The lambda values are stored symmetrically around homo/lumo
	size_t max_lambda_idx = vlx->scf.homo_idx;

	if (max_lambda_idx == 0) {
		MD_LOG_ERROR("Internal error: Incorrect max lambda");
		return false;
	}

	if (lambda_idx >= max_lambda_idx) {
		MD_LOG_ERROR("Invalid lambda index!");
		return false;
	}

	int64_t mo_idx = 0;
	if (type == MD_VLX_NTO_TYPE_PARTICLE) {
		mo_idx = (int64_t)vlx->scf.lumo_idx + (int64_t)lambda_idx;
	} else if (type == MD_VLX_NTO_TYPE_HOLE) {
		mo_idx = (int64_t)vlx->scf.homo_idx - (int64_t)lambda_idx;
	} else {
		MD_LOG_ERROR("Invalid NTO type!");
		return false;
	}

	const md_vlx_orbital_t* orb = &vlx->rsp.nto[nto_idx];

	size_t temp_pos = md_temp_get_pos();
	size_t num_mo_coeffs = number_of_mo_coefficients(orb);
	double* mo_coeffs = md_temp_push(sizeof(double) * num_mo_coeffs);

	extract_mo_coefficients(mo_coeffs, orb, mo_idx);
	extract_pgto_data(pgtos, vlx->atom_coordinates, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set, mo_coeffs);

	md_temp_set_pos_back(temp_pos);
	return true;
}

size_t md_vlx_mo_gto_count(const md_vlx_t* vlx) {
	ASSERT(vlx);
	// @NOTE: This needs to be modified in the case of Unrestricted Open Shell type.
	// In such case, we need to expand the pgtos with the contribution of Beta electrons
	return vlx_pgto_count(vlx);
}

bool md_vlx_mo_gto_extract(md_gto_t* gtos, const md_vlx_t* vlx, size_t mo_idx, md_vlx_mo_type_t type) {
	ASSERT(gtos);
	ASSERT(vlx);

	const md_vlx_orbital_t* orb = 0;
	if (type == MD_VLX_MO_TYPE_ALPHA) {
		orb = &vlx->scf.alpha;
	} else if (type == MD_VLX_MO_TYPE_BETA) {
		orb = &vlx->scf.beta;
	} else {
		MD_LOG_ERROR("Invalid MO type!");
		return false;
	}

	if (mo_idx >= number_of_molecular_orbitals(orb)) {
		MD_LOG_ERROR("Invalid mo index!");
		return false;
	}

	size_t temp_pos = md_temp_get_pos();
	size_t num_mo_coeffs = number_of_mo_coefficients(orb);
	double* mo_coeffs = md_temp_push(sizeof(double) * num_mo_coeffs);

	extract_mo_coefficients(mo_coeffs, orb, mo_idx);
	extract_pgto_data(gtos, vlx->atom_coordinates, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set, mo_coeffs);

	md_temp_set_pos_back(temp_pos);
	return true;
}

size_t md_vlx_rsp_number_of_excited_states(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.number_of_excited_states;
	}
	return 0;
}

const dvec3_t* md_vlx_rsp_electric_transition_dipole_moments(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.electric_transition_dipoles;
	}
	return NULL;
}

const dvec3_t* md_vlx_rsp_magnetic_transition_dipole_moments(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.magnetic_transition_dipoles;
	}
	return NULL;
}

const dvec3_t* md_vlx_rsp_velocity_transition_dipole_moments(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.velocity_transition_dipoles;
	}
	return NULL;
}

const double* md_vlx_rsp_rotatory_strengths(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.rotatory_strengths;
	}
	return NULL;
}

const double* md_vlx_rsp_oscillator_strengths(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.oscillator_strengths;
	}
	return NULL;
}

const double* md_vlx_rsp_absorption_ev(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.absorption_ev;
	}
	return NULL;
}

const double* md_vlx_rsp_nto_occupancy(const md_vlx_t* vlx, size_t nto_idx) {
	if (vlx) {
		if (vlx->rsp.nto && nto_idx < vlx->rsp.number_of_excited_states) {
			return vlx->rsp.nto[nto_idx].occupancy.data;
		}
	}
	return NULL;
}

const double* md_vlx_rsp_nto_energy(const md_vlx_t* vlx, size_t nto_idx) {
	if (vlx) {
		if (vlx->rsp.nto && nto_idx < vlx->rsp.number_of_excited_states) {
			return vlx->rsp.nto[nto_idx].energy.data;
		}
	}
	return NULL;
}

md_vlx_t* md_vlx_create(md_allocator_i* backing) {
	ASSERT(backing);
	md_allocator_i* arena = md_arena_allocator_create(backing, MEGABYTES(1));
	ASSERT(arena);
	md_vlx_t* vlx = md_alloc(arena, sizeof(md_vlx_t));
	if (!vlx) {
		MD_LOG_ERROR("Failed to allocate memory for veloxchem object");
		return vlx;
	}
	MEMSET(vlx, 0, sizeof(md_vlx_t));
	vlx->arena = arena;
	return vlx;
}

void md_vlx_reset(md_vlx_t* vlx) {
	if (vlx) {
		ASSERT(vlx->arena);
		md_allocator_i* arena = vlx->arena;
		md_arena_allocator_reset(arena);
		MEMSET(vlx, 0, sizeof(md_vlx_t));
		vlx->arena = arena;
	}
}

void md_vlx_destroy(md_vlx_t* vlx) {
	if (vlx) {
		md_arena_allocator_destroy(vlx->arena);
	} else {
		MD_LOG_DEBUG("Attempt to destroy NULL vlx object");
	}
}

bool md_vlx_molecule_init(md_molecule_t* mol, const md_vlx_t* vlx, md_allocator_i* alloc) {
	ASSERT(mol);
	ASSERT(vlx);

	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("The veloxchem object contains no atoms");
		return false;
	}

	mol->atom.count = vlx->number_of_atoms;
	md_array_resize(mol->atom.element, mol->atom.count, alloc);
	md_array_resize(mol->atom.type, mol->atom.count, alloc);
	md_array_resize(mol->atom.x, mol->atom.count, alloc);
	md_array_resize(mol->atom.y, mol->atom.count, alloc);
	md_array_resize(mol->atom.z, mol->atom.count, alloc);

	for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
		mol->atom.element[i] = vlx->atomic_numbers[i];
		mol->atom.type[i] = make_label(md_util_element_symbol(vlx->atomic_numbers[i]));
		mol->atom.x[i] = (float)vlx->atom_coordinates[i].x;
		mol->atom.y[i] = (float)vlx->atom_coordinates[i].y;
		mol->atom.z[i] = (float)vlx->atom_coordinates[i].z;
	}

	return true;
}

static bool vlx_mol_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
	(void)mol;
	(void)str;
	(void)arg;
	(void)alloc;
	MD_LOG_ERROR("This is not implemented yeti");
	return false;
}

static bool vlx_mol_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
	(void)arg;
	md_vlx_t* vlx = md_vlx_create(md_get_heap_allocator());

	bool success = false;
	if (vlx_parse_file(vlx, filename, VLX_FLAG_CORE)) {
		success = md_vlx_molecule_init(mol, vlx, alloc);
	}

	md_vlx_destroy(vlx);
	return success;
}

static md_molecule_loader_i vlx_loader = {
	vlx_mol_init_from_str,
	vlx_mol_init_from_file
};

// Externally visible procedures

md_molecule_loader_i* md_vlx_molecule_api(void) {
	return &vlx_loader;
}

size_t md_vlx_number_of_atoms(const md_vlx_t* vlx) {
	if (vlx) return vlx->number_of_atoms;
	return 0;
}

size_t md_vlx_number_of_alpha_electrons(const md_vlx_t* vlx) {
	if (vlx) return vlx->number_of_alpha_electrons;
	return 0;
}

size_t md_vlx_number_of_beta_electrons(const md_vlx_t* vlx) {
	if (vlx) return vlx->number_of_beta_electrons;
	return 0;
}

double md_vlx_molecular_charge(const md_vlx_t* vlx) {
	if (vlx) return vlx->molecular_charge;
	return 0;
}

double md_vlx_nuclear_repulsion(const md_vlx_t* vlx) {
	if (vlx) return vlx->nuclear_repulsion;
	return 0;
}

size_t md_vlx_spin_multiplicity(const md_vlx_t* vlx) {
	if (vlx) return vlx->spin_multiplicity;
	return 0;
}

str_t md_vlx_basis_set_ident(const md_vlx_t* vlx) {
	if (vlx) return vlx->basis_set_ident;
	return (str_t){0};
}

str_t md_vlx_dft_func_label(const md_vlx_t* vlx) {
	if (vlx) return vlx->dft_func_label;
	return (str_t){0};
}

str_t md_vlx_potfile(const md_vlx_t* vlx) {
	if (vlx) return vlx->potfile_text;
	return (str_t){0};
}

const dvec3_t* md_vlx_atom_coordinates(const md_vlx_t* vlx) {
	if (vlx) return vlx->atom_coordinates;
	return NULL;
}

const uint8_t* md_vlx_atomic_numbers(const md_vlx_t* vlx) {
	if (vlx) return vlx->atomic_numbers;
	return NULL;
}

md_vlx_scf_type_t md_vlx_scf_type(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.type;
	return MD_VLX_SCF_TYPE_UNKNOWN;
}

dvec3_t md_vlx_scf_ground_state_dipole_moment(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.ground_state_dipole_moment;
	return (dvec3_t){0};
}

size_t md_vlx_scf_homo_idx(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.homo_idx;
	}
	return 0;
}

size_t md_vlx_scf_lumo_idx(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.lumo_idx;
	}
	return 0;
}

size_t md_vlx_scf_number_of_molecular_orbitals(const md_vlx_t* vlx) {
	if (vlx) {
		return number_of_molecular_orbitals(&vlx->scf.alpha);
	}
	return 0;
}

const double* md_vlx_scf_mo_occupancy(const md_vlx_t* vlx, md_vlx_mo_type_t type) {
	if (vlx) {
		if (type == MD_VLX_MO_TYPE_ALPHA) {
			return vlx->scf.alpha.occupancy.data;
		} 
		else if (type == MD_VLX_MO_TYPE_BETA) {
			return vlx->scf.beta.occupancy.data;
		}
	}
	return NULL;
}

const double* md_vlx_scf_mo_energy(const md_vlx_t* vlx, md_vlx_mo_type_t type) {
	if (vlx) {
		if (type == MD_VLX_MO_TYPE_ALPHA) {
			return vlx->scf.alpha.energy.data;
		}
		else if (type == MD_VLX_MO_TYPE_BETA) {
			return vlx->scf.beta.energy.data;
		}
	}
	return NULL;
}

// SCF History
size_t md_vlx_scf_history_size(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.number_of_iterations;
	return 0;
}

const double* md_vlx_scf_history_energy(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.energy;
	return NULL;
}

const double* md_vlx_scf_history_energy_diff(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.energy_diff;
	return NULL;
}

const double* md_vlx_scf_history_density_diff(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.density_diff;
	return NULL;
}

const double* md_vlx_scf_history_gradient_norm(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.gradient_norm;
	return NULL;
}

const double* md_vlx_scf_history_max_gradient(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.max_gradient;
	return NULL;
}

bool md_vlx_parse_file(md_vlx_t* vlx, str_t filename) {
	return vlx_parse_file(vlx, filename, VLX_FLAG_ALL);
}
