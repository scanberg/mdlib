#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <core/md_simd.h>
#include <md_util.h>
#include <md_molecule.h>
#include <md_gto.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <float.h>

enum {
	VLX_FLAG_GEOM  = 1,
	VLX_FLAG_BASIS = 2,
	VLX_FLAG_SCF   = 4,
	VLX_FLAG_RSP   = 8,
	VLX_FLAG_ALL   = -1,
};

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
	// 0 would be a NULL entry, 1 corresponds to Hydrogen, 2 corresponds to Helium etc.
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
	md_vlx_2d_data_t C; // Coefficients
	md_vlx_2d_data_t D; // ?
	md_vlx_1d_data_t E; // Energy
	md_vlx_2d_data_t F; // Fock-matrix
	md_vlx_1d_data_t occupancy;
} md_vlx_orbital_t;

typedef struct md_vlx_scf_history_t {
	size_t  number_of_steps;
	double* density_diff;
	double* energy_diff;
	double* energy;
	double* gradient_norm;
	double* max_gradient;
} md_vlx_scf_history_t;

// Self Consistent Field
typedef struct md_vlx_scf2_t {
	str_t type;
	uint32_t homo_idx;
	uint32_t lumo_idx;

	dvec3_t ground_state_dipole_moment;

	md_vlx_orbital_t alpha;
	md_vlx_orbital_t beta;

	md_vlx_2d_data_t S;
	md_vlx_scf_history_t history;
} md_vlx_scf2_t;

typedef struct md_vlx_rsp2_t {
	size_t number_of_exited_states;

	dvec3_t* electronic_dipole_moment;
	dvec3_t* magnetic_dipole_moment;

	md_vlx_2d_data_t* nto_orb;
	md_vlx_1d_data_t* nto_occ;
} md_vlx_rsp2_t;

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
	md_vlx_scf2_t scf;
	md_vlx_rsp2_t rsp;
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

typedef struct vlx_molecule_t {
	size_t num_atoms;
	md_element_t* atomic_number;
	double* coord_x;
	double* coord_y;
	double* coord_z;
} vlx_molecule_t;

static basis_set_basis_t* basis_set_get_atom_basis(const basis_set_t* basis_set, int atomic_number) {
	if (atomic_number < basis_set->atom_basis.count) {
		return basis_set->atom_basis.data + atomic_number;
	}
	return NULL;
}

static int compute_max_angular_momentum(const basis_set_t* basis_set, const md_element_t* atomic_numbers, size_t count) {
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

static size_t compPhiAtomicOrbitals(double* out_phi, size_t phi_cap, const vlx_molecule_t* molecule,
	const basis_set_t* basis_set,
	double xp,
	double yp,
	double zp)
{
	int natoms = (int)molecule->num_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, molecule->atomic_number, molecule->num_atoms);

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
				const double factor = 1.0 / 0.529177210903;
				double rx = (xp - molecule->coord_x[atomidx]) * factor;
				double ry = (yp - molecule->coord_y[atomidx]) * factor;
				double rz = (zp - molecule->coord_z[atomidx]) * factor;
				double r2 = rx * rx + ry * ry + rz * rz;

				// process atomic orbitals
				int idelem = molecule->atomic_number[atomidx];

				basis_func_range_t range = basis_get_atomic_angl_basis_func_range(basis_set, idelem, angl);
				for (int funcidx = range.beg; funcidx < range.end; funcidx++, aoidx++) {
					double phiao = 0.0;

					basis_func_t basis_func = get_basis_func(basis_set, funcidx);

					// process primitives
					int nprims = basis_func.count;
					const double* exponents = basis_func.exponents;
					const double* normcoefs = basis_func.normalization_coefficients;

					for (int iprim = 0; iprim < nprims; iprim++) {
						double expon = exp(-exponents[iprim] * r2);
						double coef1 = normcoefs[iprim];

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

static size_t vlx_pgto_count(const md_vlx_data_t* vlx_data) {
	int natoms = (int)vlx_data->geom.num_atoms;
	int max_angl = compute_max_angular_momentum(vlx_data->basis.basis_set, vlx_data->geom.atomic_number, vlx_data->geom.num_atoms);

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
				int idelem = vlx_data->geom.atomic_number[atomidx];

				// process atomic orbitals
				basis_func_range_t range = basis_get_atomic_angl_basis_func_range(vlx_data->basis.basis_set, idelem, angl);
				for (int funcidx = range.beg; funcidx < range.end; funcidx++, aoidx++) {
					// process primitives
					basis_func_t basis_func = get_basis_func(vlx_data->basis.basis_set, funcidx);
					count += basis_func.count * ncomp;
				}
			}
		}
	}

	return count;
}

static size_t extract_pgto_data(md_gto_t* pgtos, const vlx_molecule_t* molecule, const basis_set_t* basis_set, const double* mo_coeffs) {
	int natoms = (int)molecule->num_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, molecule->atomic_number, molecule->num_atoms);

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
				const double factor = 1.0 / 0.529177210903;
				float x = (float)(molecule->coord_x[atomidx] * factor);
				float y = (float)(molecule->coord_y[atomidx] * factor);
				float z = (float)(molecule->coord_z[atomidx] * factor);

				int idelem = molecule->atomic_number[atomidx];

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

static bool parse_vlx_geom(md_vlx_geom_t* geom, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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
			MD_LOG_ERROR("Incorrect number of atoms parsed, expected (%zu) entries, parsed (%zu).", geom->num_atoms, count);
			goto done;
		}

		// Copy data
		geom->molecular_charge		= (int)field_vals[0];
		geom->spin_multiplicity		= (int)field_vals[1];
		geom->num_atoms				= (size_t)field_vals[2];
		geom->num_alpha_electrons	= (size_t)field_vals[3];
		geom->num_beta_electrons	= (size_t)field_vals[4];
	}

	md_array_grow(geom->atomic_number, count, alloc);
	md_array_grow(geom->atom_symbol, count, alloc);
	md_array_grow(geom->coord_x, count, alloc);
	md_array_grow(geom->coord_y, count, alloc);
	md_array_grow(geom->coord_z, count, alloc);

	for (size_t i = 0; i < count; ++i) {
		geom->atomic_number[i] = md_util_element_lookup_ignore_case(LBL_TO_STR(fields[i].sym));
		MEMCPY(&geom->atom_symbol[i], &fields[i].sym, sizeof(md_label_t));
		geom->coord_x[i] = fields[i].x;
		geom->coord_y[i] = fields[i].y;
		geom->coord_z[i] = fields[i].z;
	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);
	return result;
}

static bool parse_vlx_basis(md_vlx_basis_t* basis, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	md_buffered_reader_skip_line(reader); // ====...
	md_buffered_reader_skip_line(reader); // *empty*

	str_t line;
	str_t tok[8];
	int mask = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok == 2 && str_eq(tok[0], STR_LIT("Basis:"))) {
			basis->ident = str_copy(tok[1], alloc);
			mask |= 1;
		} else if (num_tok == 5 && str_eq(tok[3], STR_LIT(":"))) {
			str_t first = str_join(tok[0], tok[2]);
			if (str_eq(first, STR_LIT("Contracted Basis Functions"))) {
				basis->num_contracted_basis_functions = parse_int(tok[4]);
				mask |= 2;
			} else if (str_eq(first, STR_LIT("Primitive Basis Functions"))) {
				basis->num_primitive_basis_functions = parse_int(tok[4]);
				mask |= 4;
				break;
			}
		} else if (num_tok == 1 && str_begins_with(tok[0], STR_LIT("====="))) {
			// Parsed into next section >.<
			return false;
		}
	}

	return mask == 7;
}

static bool parse_vlx_scf(md_vlx_scf_t* scf, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[8];
	str_t line;

	md_buffered_reader_skip_line(reader); // =====

	int mask = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);
		if (str_begins_with(line, STR_LIT("Iter. |")) && str_ends_with(line, STR_LIT("| Density Change"))) {
			md_buffered_reader_skip_line(reader); // -----
			// Parse table, start tokenization
			while (md_buffered_reader_extract_line(&line, reader)) {
				size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
				if (num_tok == 6) {
					int    iteration		= (int)parse_int(tok[0]);
					double energy_tot		= parse_float(tok[1]);
					double energy_change	= parse_float(tok[2]);
					double gradient_norm	= parse_float(tok[3]);
					double max_gradient		= parse_float(tok[4]);
					double density_change	= parse_float(tok[5]);

					md_array_push(scf->iter.iteration, iteration, alloc);
					md_array_push(scf->iter.energy_total,  energy_tot, alloc);
					md_array_push(scf->iter.energy_change, energy_change, alloc);
					md_array_push(scf->iter.gradient_norm, gradient_norm, alloc);
					md_array_push(scf->iter.max_gradient, max_gradient, alloc);
					md_array_push(scf->iter.density_change, density_change, alloc);
					scf->iter.count += 1;

					mask |= 1;
				} else if (num_tok == 0) {
					// Assume valid end here upon empty line
					break;
				} else {
					MD_LOG_ERROR("Unexpected number of tokens in scf energy iteration section, expected 6, got (%zu)");
					return false;
				}
			}
		} else if (str_begins_with(line, STR_LIT("Total Energy")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			scf->total_energy = parse_float(tok[3]);
			mask |= 2;
		} else if (str_begins_with(line, STR_LIT("Electronic Energy")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			scf->electronic_energy = parse_float(tok[3]);
			mask |= 4;
		} else if (str_begins_with(line, STR_LIT("Nuclear Repulsion Energy")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 6) {
			scf->nuclear_repulsion_energy = parse_float(tok[4]);
			mask |= 8;
		} else if (str_begins_with(line, STR_LIT("Gradient Norm")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			scf->gradient_norm = parse_float(tok[3]);
			mask |= 16;
		} else if (str_eq(line, STR_LIT("Ground State Dipole Moment"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // *empty*
			double vec[3];
			for (int i = 0; i < 3; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 6) {
					MD_LOG_ERROR("Failed to parse SCF Ground State Dipole Moment, incomplete fields!");
					return false;
				}
				vec[i] = parse_float(tok[2]);
			}
			scf->ground_state_dipole_moment.ident = STR_LIT("Ground State");
			scf->ground_state_dipole_moment.x = vec[0];
			scf->ground_state_dipole_moment.y = vec[1];
			scf->ground_state_dipole_moment.z = vec[2];
			mask |= 32;
		} else if (str_begins_with(line, STR_LIT("====="))) {
			// We've read too far and into the next section
			MD_LOG_ERROR("Failed to parse SCF section, some fields are missing");
			return false;
		}
		if (mask == 63) {
			return true;
		}
	}

	return false;
}

static bool parse_vlx_rsp_dipole_moments(md_vlx_dipole_moment_t* moments, size_t num_excited_states, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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

		moments[i].ident = str_copy(ident, alloc);
		moments[i].x = parse_float(tok[3]);
		moments[i].y = parse_float(tok[4]);
		moments[i].z = parse_float(tok[5]);
	}

	return true;
}

static bool parse_vlx_rsp(md_vlx_rsp_t* rsp, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[10];
	str_t line;

	md_buffered_reader_skip_line(reader); // =====

	int mask = 0;
	while ( md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);
		if (str_empty(line)) continue;

		if (str_begins_with(line, STR_LIT("Number of States")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			rsp->num_excited_states = parse_int(tok[4]);
			mask |= 1;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole length, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->electronic_transition_length, rsp->num_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(rsp->electronic_transition_length, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 2;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole velocity, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->electronic_transition_velocity, rsp->num_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(rsp->electronic_transition_velocity, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 4;
		} else if (str_eq(line, STR_LIT("Magnetic Transition Dipole Moments (a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->magnetic_transition, rsp->num_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(rsp->magnetic_transition, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 8;
		} else if (str_eq(line, STR_LIT("One-Photon Absorption"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_array_resize(rsp->absorption_ev, rsp->num_excited_states, alloc);
			md_array_resize(rsp->absorption_osc_str, rsp->num_excited_states, alloc);
			for (size_t i = 0; i < rsp->num_excited_states; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || !extract_tokens(tok, ARRAY_SIZE(tok), &line) == 9) {
					MD_LOG_ERROR("Unexpected number of tokens in entry when parsing One-Photon Absorption");
					return false;
				}
				rsp->absorption_ev[i] = parse_float(tok[5]);
				rsp->absorption_osc_str[i] = parse_float(tok[8]);
			}
			mask |= 16;
		} else if (str_eq(line, STR_LIT("Electronic Circular Dichroism"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_array_resize(rsp->electronic_circular_dichroism_cgs, rsp->num_excited_states, alloc);
			for (size_t i = 0; i < rsp->num_excited_states; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || !extract_tokens(tok, ARRAY_SIZE(tok), &line) == 9) {
					MD_LOG_ERROR("Unexpected number of tokens in entry when parsing Electronic Circular Dichroism");
					return false;
				}
				rsp->electronic_circular_dichroism_cgs[i] = parse_float(tok[6]);
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

bool vlx_parse_out(md_vlx_data_t* vlx, md_buffered_reader_t* reader) {
	str_t line;
	while (md_buffered_reader_extract_line(&line, reader)) {
		str_t str = str_trim(line);
		if (str_eq(str, STR_LIT("Molecular Geometry (Angstroms)"))) {
			if (!parse_vlx_geom(&vlx->geom, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse geometry");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Molecular Basis (Atomic Basis)"))) {
			if (!parse_vlx_basis(&vlx->basis, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse basis");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Self Consistent Field Driver Setup"))) {
			if (!parse_vlx_scf(&vlx->scf, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse SCF section");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Linear Response EigenSolver Setup"))) {
			if (!parse_vlx_rsp(&vlx->rsp, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse RSP section");
				return false;
			}
		}
	}

	return true;
}

bool vlx_load_orbital_h5_data(md_vlx_orbitals_t* orb, str_t filename, str_t ident, md_allocator_i* alloc) {
	ASSERT(orb);

	hid_t  file_id = 0, dataset_id = 0, space_id = 0; /* identifiers */
	herr_t status = 0;

	bool result = false;

	/* Open an existing file. */
	file_id = H5Fopen(str_beg(filename), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id > 0) {

		char lbl[128];
		snprintf(lbl, sizeof(lbl), "/" STR_FMT "_orbitals", STR_ARG(ident));

		dataset_id = H5Dopen(file_id, lbl, H5P_DEFAULT);
		if (dataset_id > 0) {
			space_id = H5Dget_space(dataset_id);
			if (space_id > 0) {
				int ndim = H5Sget_simple_extent_ndims(space_id);
				if (ndim != 2) {
					MD_LOG_ERROR("Unexpected number of dimensions when reading h5 orbitals");
					goto done_orb;
				}

				hsize_t dims[2];
				ndim = H5Sget_simple_extent_dims(space_id, dims, 0);
				if (ndim != 2) {
					MD_LOG_ERROR("Unexpected number of dimensions when reading h5 orbitals");
					goto done_orb;
				}		

				double* data = md_alloc(alloc, sizeof(double) * dims[0] * dims[1]);
				if (!data) {
					MD_LOG_ERROR("An error occured when allocating data for h5 orbital");
					goto done_orb;
				}

				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

				if (status != 0) {
					MD_LOG_ERROR("An error occured when reading h5 orbital");
					md_free(alloc, data, sizeof(double) * dims[0] * dims[1]);
					goto done_orb;
				}

				orb->orbitals.dim[0] = dims[0];
				orb->orbitals.dim[1] = dims[1];
				orb->orbitals.data   = data;

			done_orb:
				H5Sclose(space_id);
			}
			H5Dclose(dataset_id);
		}

		snprintf(lbl, sizeof(lbl), "/" STR_FMT "_energies", STR_ARG(ident));

		dataset_id = H5Dopen(file_id, lbl, H5P_DEFAULT);
		if (dataset_id > 0) {
			space_id = H5Dget_space(dataset_id);
			if (space_id > 0) {
				int ndim = H5Sget_simple_extent_ndims(space_id);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in h5 energies");
					goto done_ener;
				}

				hsize_t dim;
				ndim = H5Sget_simple_extent_dims(space_id, &dim, 0);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in h5 energies");
					goto done_ener;
				}		

				double* data = md_alloc(alloc, sizeof(double) * dim);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

				if (status != 0) {
					MD_LOG_ERROR("Could not read h5 energies");
					goto done_ener;
				}

				orb->energies.count = dim;
				orb->energies.data = data;

				done_ener:
				H5Sclose(space_id);
			}
			H5Dclose(dataset_id);
		}

		snprintf(lbl, sizeof(lbl), "/" STR_FMT "_occupations", STR_ARG(ident));

		dataset_id = H5Dopen(file_id, lbl, H5P_DEFAULT);
		if (dataset_id > 0) {
			space_id = H5Dget_space(dataset_id);
			if (space_id > 0) {
				int ndim = H5Sget_simple_extent_ndims(space_id);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in h5 occupations");
					goto done_occ;
				}

				hsize_t dim;
				ndim = H5Sget_simple_extent_dims(space_id, &dim, 0);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in h5 occupations");
					goto done_occ;
				}

				double* data = md_alloc(alloc, sizeof(double) * dim);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

				if (status != 0) {
					MD_LOG_ERROR("Could not read h5 occupations");
					goto done_occ;
				}

				orb->occupations.count = dim;
				orb->occupations.data = data;
			done_occ:
				H5Sclose(space_id);
			}
		}
		status = H5Fclose(file_id);
	}

	result = orb->orbitals.dim[0] != 0 && orb->energies.count != 0 && orb->occupations.count != 0;

	return result;
}

bool md_vlx_read_scalar(void* buf, hid_t file_id, hid_t mem_type_id, const char* field_name) {
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

bool md_vlx_read_str_field(str_t* str, hid_t file_id, const char* field_name, md_allocator_i* alloc) {
	bool result = false;

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

bool md_vlx_read_multi_dim_field_fixed(void* out_data, const size_t dims[], int num_dim, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	ASSERT(out_data);
	ASSERT(dims);

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

bool md_vlx_read_multi_dim_field(void** out_data, size_t* out_dims, int num_dim, hid_t file_id, hid_t mem_type_id, const char* field_name, md_allocator_i* alloc) {
	ASSERT(out_data);
	ASSERT(out_dims);
	
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

	ndim = H5Sget_simple_extent_dims(space_id, out_dims, 0);
	ASSERT(ndim == num_dim);

	size_t dim_size = out_dims[0];
	for (int i = 1; i < ndim; ++i) {
		dim_size *= out_dims[i];
	}

	double* data = md_alloc(alloc, sizeof(double) * dim_size);
	if (!data) {
		MD_LOG_ERROR("An error occured when allocating data for H5 dataset");
		MEMSET(out_dims, 0, sizeof(size_t) * num_dim);
		goto done;
	}

	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	if (status != 0) {
		MD_LOG_ERROR("An error occured when reading H5 orbital");
		md_free(alloc, data, sizeof(double) * dim_size);
		MEMSET(out_dims, 0, sizeof(size_t) * num_dim);
		goto done;
	}

	*out_data = data;

	result = true;
done:
	H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

/*
bool md_vlx_read_scf_results(md_vlx_scf_results_t* scf, str_t filename, md_allocator_i* alloc) {
	ASSERT(scf);

	hid_t  file_id = 0;

	bool result = false;
	char buf[2048];

	// Ensure a zero terminated string for interfacing to HDF5
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	// Open an existing file
	file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	if (!md_vlx_read_str_field(&scf->type, file_id, "scf_type", alloc)) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->energy, file_id, H5T_NATIVE_DOUBLE, "scf_energy")) {
		return false;
	}

	if (!md_vlx_read_str_field(&scf->basis_set, file_id, "basis_set", alloc)) {
		return false;
	}

	if (!md_vlx_read_str_field(&scf->dft_func_label, file_id, "dft_func_label", alloc)) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->nuclear_repulsion, file_id, H5T_NATIVE_DOUBLE, "nuclear_repulsion")) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->number_of_alpha_electrons, file_id, H5T_NATIVE_INT64, "number_of_alpha_electrons")) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->atom.number_of_atoms, file_id, H5T_NATIVE_INT64, "number_of_atoms")) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->number_of_beta_electrons, file_id, H5T_NATIVE_INT64, "number_of_beta_electrons")) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->molecular_charge, file_id, H5T_NATIVE_DOUBLE, "molecular_charge")) {
		return false;
	}

	if (!md_vlx_read_str_field(&scf->potfile_text, file_id, "potfile_text", alloc)) {
		return false;
	}

	if (!md_vlx_read_scalar(&scf->spin_multiplicity, file_id, H5T_NATIVE_INT64, "spin_multiplicity")) {
		return false;
	}

	size_t coord_dim[2];
	if (!md_vlx_read_multi_dim_field(&scf->atom.coordinates, coord_dim, 2, file_id, H5T_NATIVE_DOUBLE, "atom_coordinates", alloc)) {
		return false;
	}

	size_t nuclear_charge_dim;
	if (!md_vlx_read_multi_dim_field(&scf->atom.nuclear_charges, &nuclear_charge_dim, 1, file_id, H5T_NATIVE_INT32, "nuclear_charges", alloc)) {
		return false;
	}

	// Alpha data
	if (!md_vlx_read_multi_dim_field(&scf->alpha.C.data, scf->alpha.C.size, 2, file_id, H5T_NATIVE_DOUBLE, "C_alpha", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->alpha.D.data, scf->alpha.D.size, 2, file_id, H5T_NATIVE_DOUBLE, "D_alpha", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->alpha.E.data, &scf->alpha.E.size, 1, file_id, H5T_NATIVE_DOUBLE, "E_alpha", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->alpha.F.data, scf->alpha.F.size, 2, file_id, H5T_NATIVE_DOUBLE, "F_alpha", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->alpha.occupancy.data, &scf->alpha.occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "occ_alpha", alloc)) {
		return false;
	}

	// Beta
	if (!md_vlx_read_multi_dim_field(&scf->beta.C.data, scf->beta.C.size, 2, file_id, H5T_NATIVE_DOUBLE, "C_beta", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->beta.D.data, scf->beta.D.size, 2, file_id, H5T_NATIVE_DOUBLE, "D_beta", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->beta.E.data, &scf->beta.E.size, 1, file_id, H5T_NATIVE_DOUBLE, "E_beta", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->beta.F.data, scf->beta.F.size, 2, file_id, H5T_NATIVE_DOUBLE, "F_beta", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->beta.occupancy.data, &scf->beta.occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "occ_beta", alloc)) {
		return false;
	}

	if (!md_vlx_read_multi_dim_field(&scf->S.data, scf->S.size, 2, file_id, H5T_NATIVE_DOUBLE, "S", alloc)) {
		return false;
	}

	size_t dipole_dim[1] = {3};
	if (!md_vlx_read_multi_dim_field_fixed(&scf->dipole_moment, dipole_dim, 1, file_id, H5T_NATIVE_DOUBLE, "dipole_moment")) {
		return false;
	}

	size_t scf_dim;
	if (!md_vlx_read_multi_dim_field(&scf->history.diff_density,  &scf_dim, 1, file_id, H5T_NATIVE_DOUBLE, "scf_history_diff_density", alloc)) {
		return false;
	}
	if (!md_vlx_read_multi_dim_field(&scf->history.diff_energy,   &scf_dim, 1, file_id, H5T_NATIVE_DOUBLE, "scf_history_diff_energy", alloc)) {
		return false;
	}
	if (!md_vlx_read_multi_dim_field(&scf->history.energy,	      &scf_dim, 1, file_id, H5T_NATIVE_DOUBLE, "scf_history_energy", alloc)) {
		return false;
	}
	if (!md_vlx_read_multi_dim_field(&scf->history.gradient_norm, &scf_dim, 1, file_id, H5T_NATIVE_DOUBLE, "scf_history_gradient_norm", alloc)) {
		return false;
	}
	if (!md_vlx_read_multi_dim_field(&scf->history.max_gradient, &scf_dim, 1, file_id, H5T_NATIVE_DOUBLE, "scf_history_max_gradient", alloc)) {
		return false;
	}
	
	herr_t status = H5Fclose(file_id);

	return true;
}
*/

void md_vlx_data_free(md_vlx_data_t* data) {
	ASSERT(data);
	if (data->alloc) {
		md_arena_allocator_destroy(data->alloc);
	}
	MEMSET(data, 0, sizeof(md_vlx_data_t));
}

// Extract Natural Transition Orbitals PGTOs
size_t md_vlx_nto_gto_count(const md_vlx_data_t* vlx) {
	return vlx_pgto_count(vlx);
}

bool md_vlx_nto_gto_extract(md_gto_t* pgtos, const md_vlx_data_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type) {
	ASSERT(pgtos);
	ASSERT(vlx);

	if (nto_idx >= vlx->rsp.num_excited_states) {
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

	size_t num_rows = vlx->rsp.nto[nto_idx].orbitals.dim[0];
	size_t num_cols = vlx->rsp.nto[nto_idx].orbitals.dim[1];

	size_t temp_pos = md_temp_get_pos();
	size_t num_mo_coeffs = num_cols;
	double* mo_coeffs = md_temp_push(sizeof(double) * num_mo_coeffs);

	for (size_t i = 0; i < num_mo_coeffs; ++i) {
		mo_coeffs[i] = vlx->rsp.nto[nto_idx].orbitals.data[i * num_cols + mo_idx];
	}

	vlx_molecule_t mol = {
		.num_atoms = vlx->geom.num_atoms,
		.atomic_number = vlx->geom.atomic_number,
		.coord_x = vlx->geom.coord_x,
		.coord_y = vlx->geom.coord_y,
		.coord_z = vlx->geom.coord_z,
	};

	extract_pgto_data(pgtos, &mol, vlx->basis.basis_set, mo_coeffs);

	md_temp_set_pos_back(temp_pos);
	return true;
}

size_t md_vlx_mo_gto_count(const md_vlx_data_t* vlx) {
	ASSERT(vlx);
	// @NOTE: This needs to be modified in the case of Unrestricted Open Shell type.
	// In such case, we need to expand the pgtos with the contribution of Beta electrons
	return vlx_pgto_count(vlx);
}

bool md_vlx_mo_gto_extract(md_gto_t* pgtos, const md_vlx_data_t* vlx, size_t mo_idx) {
	ASSERT(pgtos);
	ASSERT(vlx);

	size_t num_rows = vlx->scf.alpha.orbitals.dim[0];
	size_t num_cols = vlx->scf.alpha.orbitals.dim[1];

	if (mo_idx >= num_rows) {
		MD_LOG_ERROR("Invalid mo index!");
		return false;
	}

	size_t temp_pos = md_temp_get_pos();
	size_t num_mo_coeffs = vlx->scf.alpha.orbitals.dim[1];
	double* mo_coeffs = md_temp_push(sizeof(double) * num_mo_coeffs);

	for (size_t i = 0; i < num_mo_coeffs; ++i) {
		mo_coeffs[i] = vlx->scf.alpha.orbitals.data[i * num_cols + mo_idx];
	}

	vlx_molecule_t mol = {
		.num_atoms = vlx->geom.num_atoms,
		.atomic_number = vlx->geom.atomic_number,
		.coord_x = vlx->geom.coord_x,
		.coord_y = vlx->geom.coord_y,
		.coord_z = vlx->geom.coord_z,
	};

	extract_pgto_data(pgtos, &mol, vlx->basis.basis_set, mo_coeffs);

	md_temp_set_pos_back(temp_pos);
	return true;
}

static bool vlx_data_parse_str(md_vlx_data_t* vlx, str_t str, md_allocator_i* alloc) {
	MEMSET(vlx, 0, sizeof(md_vlx_data_t));
	vlx->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));
	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	return vlx_parse_out(vlx, &reader);
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

void md_vlx_destroy(md_vlx_t* vlx) {
	if (vlx) {
		md_arena_allocator_destroy(vlx->arena);
	} else {
		MD_LOG_DEBUG("Attempt to destroy NULL vlx object");
	}
}

bool md_vlx_data_parse_str(md_vlx_data_t* vlx, str_t str, md_allocator_i* alloc) {
	return vlx_data_parse_str(vlx, str, alloc);
}

static bool vlx_data_parse_file(md_vlx_data_t* vlx, str_t filename, md_allocator_i* alloc) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	MEMSET(vlx, 0, sizeof(md_vlx_data_t));
	vlx->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));

	size_t cap = KILOBYTES(16);
	char* buf = md_temp_push(cap);
	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

	bool result = vlx_parse_out(vlx, &reader);
	md_file_close(file);

	if (!result) return false;

	if (!str_empty(vlx->basis.ident)) {
		md_strb_t sb = md_strb_create(md_get_temp_allocator());
		md_strb_fmt(&sb, "%s/" STR_FMT, MD_VLX_BASIS_FOLDER, STR_ARG(vlx->basis.ident));
		str_t basis_path = md_strb_to_str(sb);
		md_file_o* basis_file = md_file_open(basis_path, MD_FILE_READ | MD_FILE_BINARY);
		if (basis_file) {
			md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, file);
			vlx->basis.basis_set = md_alloc(vlx->alloc, sizeof(basis_set_t));
			if (!parse_basis_set(vlx->basis.basis_set, &basis_reader, vlx->alloc)) {
				MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
				return false;
			}

			normalize_basis_set(vlx->basis.basis_set);

			md_file_close(basis_file);
		}
	}

	str_t base_file = {0};
	if (!extract_file_path_without_ext(&base_file, filename)) {
		MD_LOG_ERROR("Failed to extract base file path");
		return false;
	}
	md_strb_t sb = md_strb_create(md_get_temp_allocator());

	if (vlx->scf.nuclear_repulsion_energy > 0) {
		md_strb_fmt(&sb, STR_FMT ".scf.h5", STR_ARG(base_file));
		str_t scf_path = md_strb_to_str(sb);
		if (md_path_is_valid(scf_path)) {
			if (!vlx_load_orbital_h5_data(&vlx->scf.alpha, scf_path, STR_LIT("alpha"), vlx->alloc)) {
				MD_LOG_ERROR("Failed to load orbital h5 parameters from file '" STR_FMT "'", STR_ARG(scf_path));
				return false;
			}
			for (size_t i = 0; i < vlx->scf.alpha.occupations.count; ++i) {
				if (vlx->scf.alpha.occupations.data[i] == 0.0) {
					vlx->scf.homo_idx = (size_t)MAX(0, (int64_t)i - 1);
					vlx->scf.lumo_idx = i;
					break;
				}
			}
		}
	}

	if (vlx->rsp.num_excited_states > 0) {
		for (int i = 1; i <= (int)vlx->rsp.num_excited_states; ++i) {
			md_strb_reset(&sb);
			md_strb_fmt(&sb, STR_FMT "_S%i_NTO.h5", STR_ARG(base_file), i);
			str_t nto_path = md_strb_to_str(sb);
			if (md_path_is_valid(nto_path)) {
				md_vlx_orbitals_t nto = {0};
				if (!vlx_load_orbital_h5_data(&nto, nto_path, STR_LIT("alpha"), vlx->alloc)) {
					MD_LOG_ERROR("Failed to load NTO h5 parameters from file '" STR_FMT "'", STR_ARG(nto_path));
					return false;
				}
				md_array_push(vlx->rsp.nto, nto, vlx->alloc);
			}
		}
	}

	return true;
}

bool md_vlx_data_parse_file(md_vlx_data_t* vlx, str_t filename, md_allocator_i* alloc) {
	return vlx_data_parse_file(vlx, filename, alloc);
}

bool md_vlx_molecule_init(md_molecule_t* mol, const md_vlx_data_t* vlx, md_allocator_i* alloc) {
	ASSERT(mol);
	ASSERT(vlx);

	if (vlx->geom.num_atoms == 0) {
		MD_LOG_ERROR("The veloxchem data object contains no atoms");
		return false;
	}

	mol->atom.count = vlx->geom.num_atoms;
	md_array_resize(mol->atom.type, mol->atom.count, alloc);
	md_array_resize(mol->atom.x, mol->atom.count, alloc);
	md_array_resize(mol->atom.y, mol->atom.count, alloc);
	md_array_resize(mol->atom.z, mol->atom.count, alloc);

	for (size_t i = 0; i < vlx->geom.num_atoms; ++i) {
		mol->atom.type[i] = vlx->geom.atom_symbol[i];
		mol->atom.x[i] = (float)vlx->geom.coord_x[i];
		mol->atom.y[i] = (float)vlx->geom.coord_y[i];
		mol->atom.z[i] = (float)vlx->geom.coord_z[i];
	}

	return true;
}

static bool vlx_mol_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
	(void)arg;
	md_vlx_data_t vlx = {0};
	bool success = false;
	if (vlx_data_parse_str(&vlx, str, md_get_heap_allocator())) {
		success = md_vlx_molecule_init(mol, &vlx, alloc);
	}
	md_vlx_data_free(&vlx);

	return success;
}

static bool vlx_mol_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
	(void)arg;
	md_vlx_data_t vlx = {0};
	bool success = false;
	if (vlx_data_parse_file(&vlx, filename, md_get_heap_allocator())) {
		success = md_vlx_molecule_init(mol, &vlx, alloc);
	}
	md_vlx_data_free(&vlx);

	return success;
}

static md_molecule_loader_i vlx_loader = {
	vlx_mol_init_from_str,
	vlx_mol_init_from_file
};

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

size_t md_vlx_spin_multilicity(const md_vlx_t* vlx) {
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

str_t md_vlx_scf_type(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.type;
	return (str_t){0};
}

dvec3_t md_vlx_scf_ground_state_dipole_moment(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.ground_state_dipole_moment;
	return (dvec3_t){0};
}
// SCF History
size_t md_vlx_scf_history_size(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.number_of_steps;
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

bool md_vlx_parse_scf_results_h5(md_vlx_t* vlx, str_t filename) {
	return false;
}

bool md_vlx_parse_out_file(md_vlx_t* vlx, str_t filename) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	MEMSET(vlx, 0, sizeof(md_vlx_data_t));

	md_vlx_geom_t geom = {0};
	md_vlx_basis_t basis = {0};
	md_vlx_scf_t scf = {0};
	md_vlx_rsp_t rsp = {0};

	size_t cap = KILOBYTES(16);
	char* buf = md_temp_push(cap);
	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
	str_t line;
	while (md_buffered_reader_extract_line(&line, &reader)) {
		str_t str = str_trim(line);
		if (str_eq(str, STR_LIT("Molecular Geometry (Angstroms)"))) {
			if (!parse_vlx_geom(&geom, &reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse geometry");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Molecular Basis (Atomic Basis)"))) {
			if (!parse_vlx_basis(&basis, &reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse basis");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Self Consistent Field Driver Setup"))) {
			if (!parse_vlx_scf(&scf, &reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse SCF section");
				return false;
			}
		} else if (str_eq(str, STR_LIT("Linear Response EigenSolver Setup"))) {
			if (!parse_vlx_rsp(&rsp, &reader, vlx->arena)) {
				MD_LOG_ERROR("Failed to parse RSP section");
				return false;
			}
		}
	}
	md_file_close(file);

	md_strb_t sb = md_strb_create(md_get_temp_allocator());

	if (!str_empty(vlx->basis_set_ident)) {
		md_strb_fmt(&sb, "%s/" STR_FMT, MD_VLX_BASIS_FOLDER, STR_ARG(vlx->basis_set_ident));
		md_file_o* basis_file = md_file_open(md_strb_to_str(sb), MD_FILE_READ | MD_FILE_BINARY);
		if (basis_file) {
			md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, file);
			bool result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
			md_file_close(basis_file);
			if (!result) {
				MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
				return false;
			}
			normalize_basis_set(&vlx->basis_set);
		} else {
			// Attempt to read basis set file from same folder as file
			str_t folder = {0};
			if (!extract_folder_path(&folder, filename)) {
				MD_LOG_ERROR("An error occured when extracting the path to supplied file");
				return false;
			}
			md_strb_reset(&sb);
			md_strb_push_str(&sb, folder);
			md_strb_push_str(&sb, vlx->basis_set_ident);
			basis_file = md_file_open(md_strb_to_str(sb), MD_FILE_READ | MD_FILE_BINARY);
			if (basis_file) {
				md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, file);
				bool result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
				if (!result) {
					MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
					return false;
				}
				normalize_basis_set(&vlx->basis_set);
			} else {
				MD_LOG_ERROR("Could not find basis file corresponding to identifier: '"STR_FMT"'", STR_ARG(vlx->basis_set_ident));
			}
		}
	}

	str_t base_file = {0};
	if (!extract_file_path_without_ext(&base_file, filename)) {
		MD_LOG_ERROR("Failed to extract base file path");
		return false;
	}

	{
		// Attempt to read scf data
		md_strb_reset(&sb);
		md_strb_fmt(&sb, STR_FMT ".scf.h5", STR_ARG(base_file));
		str_t scf_path = md_strb_to_str(sb);
		if (md_path_is_valid(scf_path)) {
			// Open an existing file
			hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id == H5I_INVALID_HID) {
				MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
				return false;
			}
			if (!md_vlx_read_multi_dim_field(&vlx->scf.alpha.C.data, vlx->scf.alpha.C.size, 2, file_id, H5T_NATIVE_DOUBLE, "alpha_orbitals", vlx->arena)) {
				return false;
			}
			if (!md_vlx_read_multi_dim_field(&vlx->scf.alpha.E.data, &vlx->scf.alpha.E.size, 1, file_id, H5T_NATIVE_DOUBLE, "alpha_energies", vlx->arena)) {
				return false;
			}
			if (!md_vlx_read_multi_dim_field(&vlx->scf.alpha.occupancy.data, &vlx->scf.alpha.occupancy.size, 1, file_id, H5T_NATIVE_DOUBLE, "alpha_occupation", vlx->arena)) {
				return false;
			}
			// Identify homo and lumo
			for (size_t i = 0; i < vlx->scf.alpha.occupancy.size; ++i) {
				if (vlx->scf.alpha.occupancy.data[i] == 0.0) {
					vlx->scf.homo_idx = (uint32_t)MAX(0, (int32_t)i - 1);
					vlx->scf.lumo_idx = (uint32_t)i;
					break;
				}
			}
		}
	}

	if (vlx->rsp.number_of_exited_states > 0) {
		md_array_resize(vlx->rsp.nto_orb, vlx->rsp.number_of_exited_states, vlx->arena);
		md_array_resize(vlx->rsp.nto_occ, vlx->rsp.number_of_exited_states, vlx->arena);
		MEMSET(vlx->rsp.nto_orb, 0, sizeof(md_vlx_2d_data_t) * vlx->rsp.number_of_exited_states);
		MEMSET(vlx->rsp.nto_occ, 0, sizeof(md_vlx_1d_data_t) * vlx->rsp.number_of_exited_states);

		for (int i = 1; i <= (int)vlx->rsp.number_of_exited_states; ++i) {
			md_strb_reset(&sb);
			md_strb_fmt(&sb, STR_FMT "_S%i_NTO.h5", STR_ARG(base_file), i);
			str_t nto_path = md_strb_to_str(sb);
			if (md_path_is_valid(nto_path)) {
				// Open an existing file
				hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
				if (file_id == H5I_INVALID_HID) {
					MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
					return false;
				}
				md_vlx_2d_data_t nto_orb = {0};
				md_vlx_1d_data_t nto_occ = {0};

				int idx = i - 1;
				if (!md_vlx_read_multi_dim_field(&(vlx->rsp.nto_orb[idx].data),  vlx->rsp.nto_orb[idx].size, 2, file_id, H5T_NATIVE_DOUBLE, "alpha_orbitals", vlx->arena)) {
					return false;
				}
				if (!md_vlx_read_multi_dim_field(&(vlx->rsp.nto_occ[idx].data), &vlx->rsp.nto_occ[idx].size, 1, file_id, H5T_NATIVE_DOUBLE, "alpha_occupation", vlx->arena)) {
					return false;
				}
			}
		}
	}

	return true;
}

bool md_vlx_parse_file(md_vlx_t* vlx, str_t filename) {
	if (str_ends_with(filename, STR_LIT(".out"))) {
		
	} else if (str_ends_with(filename, STR_LIT(".scf.results.h5"))) {
	
	} else if (str_ends_with(filename, STR_LIT(".h5"))) {
	
	} else {
		MD_LOG_DEBUG("Unsupported file format");
		return false;
	}

	return true;
}
