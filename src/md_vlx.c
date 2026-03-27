#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>

#include <md_util.h>
#include <md_system.h>
#include <md_trajectory.h>
#include <md_gto.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <float.h>

#define VLX_MAGIC 0x87b716a78cfb2813

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999

#define HARTREE_TO_EV 27.2114079527

typedef enum {
	VLX_FLAG_CORE = 1,
	VLX_FLAG_SCF  = 2,
	VLX_FLAG_RSP  = 4,
	VLX_FLAG_VIB  = 8,
	VLX_FLAG_OPT  = 16,
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

#define MD_VLX_MAX_RANK 5

typedef struct md_ndarray_f64_t {
	uint32_t rank;
	uint32_t size[MD_VLX_MAX_RANK];
	double*  data;
} md_ndarray_f64_t;

static inline size_t md_ndarray_f64_size(const md_ndarray_f64_t* arr) {
	size_t total = 1;
	for (uint32_t i = 0; i < arr->rank; ++i) {
		total *= arr->size[i];
	}
	return total;
}

typedef struct md_2darray_f64_t {
	size_t  size[2];
	double* data;
} md_2darray_f64_t;

typedef struct md_ndarray_i32_t {
	uint32_t rank;
	uint32_t size[MD_VLX_MAX_RANK];
	int32_t* data;
} md_ndarray_i32_t;

static inline size_t md_ndarray_i32_size(const md_ndarray_i32_t* arr) {
	size_t total = 1;
	for (uint32_t i = 0; i < arr->rank; ++i) {
		total *= arr->size[i];
	}
	return total;
}

typedef struct md_vlx_orbital_t {
    md_ndarray_f64_t coefficients;	// [num_frames][num_orbitals][num_ao]
    md_ndarray_f64_t density;		// [num_frames][num_ao][num_ao]
    md_ndarray_f64_t energy;		// [num_frames][num_orbitals]
    md_ndarray_f64_t occupancy;		// [num_frames][num_orbitals]

    int* homo_idx;					// [num_frames]
	int* lumo_idx;					// [num_frames]
} md_vlx_orbital_t;

typedef struct md_vlx_scf_history_t {
	md_ndarray_i32_t number_of_iterations;
	md_ndarray_f64_t density_diff;
	md_ndarray_f64_t energy_diff;
	md_ndarray_f64_t energy;
	md_ndarray_f64_t gradient_norm;
	md_ndarray_f64_t max_gradient;
} md_vlx_scf_history_t;

// Self Consistent Field
typedef struct md_vlx_scf_t {
	md_vlx_scf_type_t type;

	md_ndarray_f64_t energy; // [num_frames]
	md_ndarray_f64_t ground_state_dipole_moment; // [num_frames(optional)][3]

	md_vlx_orbital_t alpha;
	md_vlx_orbital_t beta;

	md_ndarray_f64_t resp_charges;
	md_ndarray_f64_t S;

	md_vlx_scf_history_t history;
} md_vlx_scf_t;

typedef struct md_vlx_cpp_t {
	size_t number_of_frequencies;
	double* frequencies;
	double* sigmas;
	double* delta_epsilon;
} md_vlx_cpp_t;

typedef struct md_vlx_rsp_t {
	size_t   number_of_excited_states;
	dvec3_t* electric_transition_dipoles;
	dvec3_t* magnetic_transition_dipoles;
	dvec3_t* velocity_transition_dipoles;
	double* rotatory_strengths;		// unit = 10^-40 cgs
	double* oscillator_strengths;
	double* absorption_ev;
	md_vlx_orbital_t* nto;
	md_vlx_cpp_t cpp;
} md_vlx_rsp_t;

typedef struct md_vlx_vib_t {
	size_t number_of_normal_modes;
	double* force_constants;
	double* ir_intensities;
	double* frequencies;
	double* reduced_masses;
	dvec3_t** normal_modes;
} md_vlx_vib_t;

typedef struct md_vlx_opt_t {
	size_t number_of_steps;
	double* nuclear_repulsion_energies;
	double* energies;
	dvec3_t* coordinates;
} md_vlx_opt_t;

typedef struct md_vlx_t {
	basis_set_t basis_set;

	str_t  basis_set_ident;
	str_t  dft_func_label;
	str_t  potfile_text;

	size_t number_of_atoms;
	size_t number_of_alpha_electrons;
	size_t number_of_beta_electrons;

    size_t number_of_frames; // For trajectory data, 1 for static calculations, >1 for trajectories
    int*   frame_ids;        // optional

	double molecular_charge;
	double nuclear_repulsion_energy;
	int    spin_multiplicity;

    // [num_frames][num_atoms][3]
	md_ndarray_f64_t atom_coordinates;

	// Data blocks
	md_vlx_scf_t scf;
	md_vlx_rsp_t rsp;
	md_vlx_vib_t vib;
	md_vlx_opt_t opt;

	// Atomic orbital data represented as GTOs
	md_gto_data_t gto_data;

	md_element_t* atomic_numbers;
	int* ao_to_atom_idx;    // Maps atomic orbitals to atom indices

	struct md_allocator_i* arena;
} md_vlx_t;

static void orbital_identify_homo_lumo(int* out_homo_idx, int* out_lumo_idx, const double* occ_data, size_t occ_size) {
	ASSERT(out_homo_idx);
	ASSERT(out_lumo_idx);

	for (size_t i = 0; i < occ_size; ++i) {
		if (occ_data[i] > 0.5) {
			*out_homo_idx = (int)i;
		} else {
			*out_lumo_idx = (int)i;
			break;
		}
	}
}

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
	if (atomic_number < (int)basis_set->atom_basis.count) {
		return basis_set->atom_basis.data + atomic_number;
	}
	return NULL;
}

static inline int compute_max_angular_momentum(const basis_set_t* basis_set, const md_element_t* atomic_numbers, size_t count) {
	ASSERT(basis_set);
	ASSERT(atomic_numbers);
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

static size_t basis_set_extract_atomic_basis_func_angl(basis_func_t* out_funcs, size_t cap_funcs, const basis_set_t* basis_set, int atomic_number, int angl) {
    size_t count = 0;

    basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
    if (atom_basis) {
        int beg = atom_basis->basis_func_offset;
        int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
        for (int i = beg; i < end; ++i) {
            if (count == cap_funcs) break;
            if (basis_set->basis_func.data[i].type == angl) {
                out_funcs[count++] = get_basis_func(basis_set, i);
            }
        }
    }

    return count;
}

// This is a ported reference implementation from VeloxChem found in VisualizationDriver.cpp
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

	basis_func_t basis_funcs[128];

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

				size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++, aoidx++) {
					double phiao = 0.0;

					basis_func_t bf = basis_funcs[funcidx];

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

	basis_func_t basis_funcs[128];

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			int	ncomp = spherical_momentum_num_factors(angl, isph);
			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				int idelem = vlx->atomic_numbers[atomidx];

				// process atomic orbitals
                size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), &vlx->basis_set, idelem, angl);
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
					// process primitives
					count += basis_funcs[funcidx].count * ncomp;
				}
			}
		}
	}

	return count;
}

static size_t extract_ao_to_atom_idx(int* out_ao_to_atom, const md_atomic_number_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	size_t count = 0;

    basis_func_t basis_funcs[128];

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// int	ncomp = spherical_momentum_num_factors(angl, isph);

			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				int idelem = atomic_numbers[atomidx];
				size_t num_ao = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);

				for (size_t iao = 0; iao < num_ao; iao++) {
					out_ao_to_atom[count] = atomidx;
					count += 1;
				}
			}
		}
	}
	return count;
}

static size_t extract_pgto_data(md_gto_t* out_gtos, int* out_atom_idx, const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set, const double* mo_coeffs) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	size_t count = 0;
	size_t mo_coeff_idx = 0;

    basis_func_t basis_funcs[128];

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

				size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);

				// process atomic orbitals
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
					const double mo_coeff = mo_coeffs ? mo_coeffs[mo_coeff_idx++] : 1.0;

					// process primitives
					basis_func_t basis_func = basis_funcs[funcidx];
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

							out_gtos[count].x		= x;
							out_gtos[count].y		= y;
							out_gtos[count].z		= z;
							out_gtos[count].coeff	= (float)coeff;
							out_gtos[count].alpha	= (float)alpha; 
							out_gtos[count].cutoff	= FLT_MAX;
							out_gtos[count].i		= (uint8_t)lx[icomp];
							out_gtos[count].j		= (uint8_t)ly[icomp];
							out_gtos[count].k		= (uint8_t)lz[icomp];
							out_gtos[count].l		= (uint8_t)angl;

							if (out_atom_idx) {
								out_atom_idx[count] = atomidx;
							}

							count += 1;
						}
					}
				}
			}
		}
	}
	return count;
}

static void extract_gto_data(struct md_gto_data_t* out_data, const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set, md_allocator_i* alloc) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	// uint32_t coeff_idx = 0; // same as the cgto_idx

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

				basis_func_t basis_funcs[128];
				size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);

				// process atomic orbitals
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
                    vec4_t   cgto_xyzr = {x, y, z, FLT_MAX};
                    uint32_t cgto_offset = (uint32_t)out_data->num_pgtos;

					// process primitives
					basis_func_t basis_func = basis_funcs[funcidx];
					ASSERT(basis_func.type == angl);
					const int        nprims = basis_func.count;
					const double* exponents = basis_func.exponents;
					const double* normcoefs = basis_func.normalization_coefficients;

					size_t new_pgto_count = out_data->num_pgtos + (size_t)nprims * (size_t)ncomp;
					md_array_ensure(out_data->pgto_alpha,  new_pgto_count, alloc);
					md_array_ensure(out_data->pgto_coeff,  new_pgto_count, alloc);
					md_array_ensure(out_data->pgto_radius, new_pgto_count, alloc);
					md_array_ensure(out_data->pgto_ijkl,   new_pgto_count, alloc);

					for (int iprim = 0; iprim < nprims; iprim++) {
						double alpha = exponents[iprim];
						double coef1 = normcoefs[iprim];

						// transform from Cartesian to spherical harmonics
						for (int icomp = 0; icomp < ncomp; icomp++) {
							uint32_t packed_ijkl = md_gto_pack_ijkl(lx[icomp], ly[icomp], lz[icomp], angl);
							md_array_push_no_grow(out_data->pgto_alpha, (float)alpha);
							md_array_push_no_grow(out_data->pgto_coeff, (float)(coef1 * fcarts[icomp]));
							md_array_push_no_grow(out_data->pgto_radius, FLT_MAX);
							md_array_push_no_grow(out_data->pgto_ijkl,  packed_ijkl);
							out_data->num_pgtos += 1;
						}
					}
					md_array_push(out_data->cgto_xyzr, cgto_xyzr, alloc);
					md_array_push(out_data->cgto_offset, cgto_offset, alloc);
					out_data->num_cgtos += 1;
				}
			}
		}
	}

	if (out_data->num_cgtos > 0) {
		md_array_push(out_data->cgto_offset, (uint32_t)out_data->num_pgtos, alloc);
	}
}

static inline double compute_overlap(basis_func_t func, int i, int j) {
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

	for (int i = 0; i < func.count; i++) {
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

	for (int i = 0; i < func.count; i++) {
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
		for (int i = 0; i < func.count; i++) {
			ovl += compute_overlap(func, i, i);
			for (int j = i + 1; j < func.count; j++) {
				ovl += 2.0 * compute_overlap(func, i, j);
			}
		}

		// renormalize primitive BFs
		ovl = 1.0 / sqrt(ovl);
		for (int i = 0; i < func.count; i++) {
			func.normalization_coefficients[i] *= ovl;
		}
	}
}

static bool parse_basis_set(basis_set_t* basis_set, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	MEMSET(basis_set, 0, sizeof(basis_set_t));

	str_t line;
	str_t tok[4];
	size_t line_count = 0;

	// Insert null_basis element for index 0
	const basis_set_basis_t null_basis = {0};

	basis_set_basis_t* curr_atom_basis = NULL;
	while (md_buffered_reader_extract_line(&line, reader)) {
		line_count += 1;
		str_t line_original = line;
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (!num_tok) continue;

		if ((num_tok == 2 || num_tok == 3) && str_eq(tok[0], STR_LIT("@BASIS_SET"))) {
			// In the case of renamed identifiers using the alias table, there is an additional token which provides the original identifier (prepended with an !)
			str_t ident = {0};
			if (num_tok == 2)
				ident = tok[1];
			else {
				ident = str_substr(tok[2], 1, SIZE_MAX);
			}
			MD_LOG_DEBUG("Parsing Basis Set with identifier: '" STR_FMT "'", STR_ARG(ident));
			basis_set->identifier = str_copy(ident, alloc);
		}
		else if (num_tok == 2 && str_eq(tok[0], STR_LIT("@ATOMBASIS"))) {
			int atomic_number = md_atomic_number_from_symbol(tok[1], true);
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
			while ((int)md_array_size(basis_set->atom_basis.data) < atomic_number) {
				md_array_push(basis_set->atom_basis.data, null_basis, alloc);
			}

			md_array_push(basis_set->atom_basis.data, atom_basis, alloc);
			curr_atom_basis = md_array_last(basis_set->atom_basis.data);

			basis_set->atom_basis.count = md_array_size(basis_set->atom_basis.data);
		}
		else if (num_tok == 1 && str_eq(tok[0], STR_LIT("@END"))) {
			curr_atom_basis = NULL;
		}
		else if (num_tok == 3) {
			int type = char_to_angular_momentum_type(tok[0].ptr[0]);
			if (type == -1) {
				MD_LOG_ERROR("Unrecognized angular momentum type '" STR_FMT "' in basis set", STR_ARG(tok[0]));
				MD_LOG_ERROR("This occured on line %zu: '" STR_FMT "'", line_count, STR_ARG(line_original));

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

static bool h5_read_scalar(void* buf, hid_t file_id, hid_t mem_type_id, const char* field_name) {
    bool result = false;

    htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
    if (exists == 0) {
        return false;
    }

    hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) {
        MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
        return false;
    }

    hid_t space_id = H5Dget_space(dataset_id);
    int ndims = H5Sget_simple_extent_ndims(space_id);

    herr_t status = -1;
    if (ndims == 2) {
        // 2D array: select [0,0]
        hsize_t offset[2] = {0, 0};
        hsize_t count[2]  = {1, 1};
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        hid_t memspace_id = H5Screate_simple(2, count, NULL);
        status = H5Dread(dataset_id, mem_type_id, memspace_id, space_id, H5P_DEFAULT, buf);
        H5Sclose(memspace_id);
    } else if (ndims == 1) {
        hsize_t offset[1] = {0};
        hsize_t count[1]  = {1};
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        hid_t memspace_id = H5Screate_simple(1, count, NULL);
        status = H5Dread(dataset_id, mem_type_id, memspace_id, space_id, H5P_DEFAULT, buf);
        H5Sclose(memspace_id);
    } else if (ndims == 0) {
        status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    } else {
        MD_LOG_ERROR("Unsupported scalar dataset rank: %d", ndims);
        goto done;
    }

    if (status != 0) {
        MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
        goto done;
    }

    result = true;
done:
    H5Sclose(space_id);
    H5Dclose(dataset_id);
    return result;
}

static bool h5_read_str(str_t* str, hid_t file_id, const char* field_name, md_allocator_i* alloc) {
    bool result = false;

    htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
    if (exists == 0) {
        return false;
    }

    hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
    if (dataset_id == H5I_INVALID_HID) {
        MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
        return false;
    }

    hid_t datatype_id = H5Dget_type(dataset_id);
    hid_t space_id = H5Dget_space(dataset_id);

    int ndims = H5Sget_simple_extent_ndims(space_id);
    size_t size = H5Tget_size(datatype_id);
    str_t data = str_alloc(size, alloc);

    herr_t status = -1;
    if (ndims == 2) {
        // 2D array: select [0,0]
        hsize_t offset[2] = {0, 0};
        hsize_t count[2]  = {1, 1};
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        hid_t memspace_id = H5Screate_simple(2, count, NULL);
        status = H5Dread(dataset_id, datatype_id, memspace_id, space_id, H5P_DEFAULT, (char*)data.ptr);
        H5Sclose(memspace_id);
    } else if (ndims == 1) {
        hsize_t offset[1] = {0};
        hsize_t count[1]  = {1};
        H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        hid_t memspace_id = H5Screate_simple(1, count, NULL);
        status = H5Dread(dataset_id, datatype_id, memspace_id, space_id, H5P_DEFAULT, (char*)data.ptr);
        H5Sclose(memspace_id);
    } else if (ndims == 0) {
        status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)data.ptr);
    } else {
        MD_LOG_ERROR("Unsupported string dataset rank: %d", ndims);
        str_free(data, alloc);
        goto done;
    }

    if (status != 0) {
        MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
        str_free(data, alloc);
        goto done;
    }
    *str = data;
    result = true;
done:
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
	hid_t space_id = H5Dget_space(dataset_id);

	// Get the rank of the field
    int rank = H5Sget_simple_extent_ndims(H5Dget_space(dataset_id));

	// If the rank is greater than 1, we just read the first element [0,0,...] as a string
	if (rank > 1) {
		hsize_t offset[8] = {0};
		hsize_t count[8] = {1};
		for (int i = 0; i < rank; ++i) count[i] = 1;
		H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		hid_t memspace_id = H5Screate_simple(rank, count, NULL);
		size_t size = H5Tget_size(datatype_id);
		if (size < str_cap) {
			herr_t status = H5Dread(dataset_id, datatype_id, memspace_id, space_id, H5P_DEFAULT, out_str);
			H5Sclose(memspace_id);
			if (status != 0) {
				MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
				goto done;
			}
			result = size;
			goto done;
		}
		H5Sclose(memspace_id);
		goto done;
	}

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
	H5Sclose(space_id);
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

	ndim = H5Sget_simple_extent_dims(space_id, (hsize_t*)dims, 0);

done:
	H5Sclose(space_id);
	H5Dclose(dataset_id);

	return ndim;
}

static bool h5_check_dataset_exists(hid_t file_id, const char* field_name) {
	htri_t exists = H5Lexists(file_id, field_name, H5P_DEFAULT);
	return exists > 0;
}

static bool h5_read_dataset_data(void* out_data, size_t num_samples, hid_t file_id, hid_t mem_type_id, const char* field_name) {
    ASSERT(out_data);

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

    hsize_t total_points = H5Sget_simple_extent_npoints(space_id);

    if (num_samples > (size_t)total_points) {
        MD_LOG_ERROR("Requested more samples than available in dataset");
        goto done;
    }

    herr_t status = -1;
    if (num_samples == (size_t)total_points) {
        // Read the whole dataset as usual
        status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_data);
    } else {
		// Create temporary buffer to read all into
        size_t temp_bytes = total_points * H5Tget_size(mem_type_id);
        void* temp_buffer = malloc(temp_bytes);
        status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_buffer);
        if (status == 0) {
            MEMCPY(out_data, temp_buffer, num_samples * H5Tget_size(mem_type_id));
        }
        free(temp_buffer);
    }

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

static bool h5_extract_ndarray_f64(md_ndarray_f64_t* out_array, hid_t file_id, const char* field_name, md_allocator_i* arena) {
	size_t dims[8];
	int ndim = h5_read_dataset_dims(dims, 8, file_id, field_name);
	if (ndim < 0) {
		return false;
	}
	ASSERT(ndim < MD_VLX_MAX_RANK);

	size_t total_size = 1;
	for (int i = 0; i < ndim; i++) {
		total_size *= dims[i];
        out_array->size[i] = (uint32_t)dims[i];
	}
	out_array->rank = ndim;

	md_array_resize(out_array->data, total_size, arena);
	return h5_read_dataset_data(out_array->data, total_size, file_id, H5T_NATIVE_DOUBLE, field_name);
}

static bool h5_extract_ndarray_i32(md_ndarray_i32_t* out_array, hid_t file_id, const char* field_name, md_allocator_i* arena) {
	size_t dims[8];
	int ndim = h5_read_dataset_dims(dims, 8, file_id, field_name);
	if (ndim < 0) {
		return false;
	}
	ASSERT(ndim < MD_VLX_MAX_RANK);

	size_t total_size = 1;
	for (int i = 0; i < ndim; i++) {
		total_size *= dims[i];
        out_array->size[i] = (uint32_t)dims[i];
	}
	out_array->rank = ndim;

	md_array_resize(out_array->data, total_size, arena);
	return h5_read_dataset_data(out_array->data, total_size, file_id, H5T_NATIVE_INT32, field_name);
}

// Data extraction procedures

static bool h5_read_scf_data(md_vlx_scf_t* scf, hid_t handle, md_allocator_i* arena) {
	char scf_type[64] = {0};
	if (!h5_read_cstr(scf_type, sizeof(scf_type), handle, "scf_type")) {
		return false;
	}

	if (str_eq_cstr(STR_LIT("restricted"), scf_type)) {
		scf->type = MD_VLX_SCF_TYPE_RESTRICTED;
	} else if (str_eq_cstr(STR_LIT("restricted_openshell"), scf_type)) {
		scf->type = MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL;
	} else if (str_eq_cstr(STR_LIT("unrestricted"), scf_type)) {
		scf->type = MD_VLX_SCF_TYPE_UNRESTRICTED;
	} else {
		scf->type = MD_VLX_SCF_TYPE_UNKNOWN;
		MD_LOG_ERROR("Unrecognized scf type present in h5 scf section: '%s'", scf_type);
		return false;
	}

	if (!h5_extract_ndarray_f64(&scf->energy, handle, "scf_energy", arena)) {
		return false;
	}

    size_t num_frames = md_ndarray_f64_size(&scf->energy);

	/*
	size_t dim[2] = {0};
	h5_read_dataset_dims(dim, 2, handle, "C_alpha");

	// Density dimensions (May differ from dim is always square)
	size_t den_dim[2] = {0};
    h5_read_dataset_dims(den_dim, 2, handle, "D_alpha");

	md_array_resize(scf->alpha.coefficients.data, dim[0] * dim[1], arena);
	MEMCPY(scf->alpha.coefficients.size, dim, sizeof(dim));

	md_array_resize(scf->alpha.energy.data, dim[1], arena);
	scf->alpha.energy.size = dim[1];
	md_array_resize(scf->alpha.occupancy.data, dim[1], arena);
	scf->alpha.occupancy.size = dim[1];

	md_array_resize(scf->alpha.density.data, den_dim[0] * den_dim[1], arena);
    MEMCPY(scf->alpha.density.size, den_dim, sizeof(den_dim));
	*/

	// Extract alpha data
	if (!h5_extract_ndarray_f64(&scf->alpha.coefficients, handle, "C_alpha", arena)) {
		return false;
	}
	if (!h5_extract_ndarray_f64(&scf->alpha.energy, handle, "E_alpha", arena)) {
		return false;
	}
	if (!h5_extract_ndarray_f64(&scf->alpha.occupancy, handle, "occ_alpha", arena)) {
		return false;
	}
    if (!h5_extract_ndarray_f64(&scf->alpha.density, handle, "D_alpha", arena)) {
        return false;
    }

	md_array_resize(scf->alpha.homo_idx, num_frames, arena);
	md_array_resize(scf->alpha.lumo_idx, num_frames, arena);

	// Identify homo and lumo
	for (size_t i = 0; i < num_frames; ++i) {
		// Create slices for occupancy data
		size_t num_occ = scf->alpha.occupancy.size[scf->alpha.occupancy.rank - 1];
		const double* occ_data = scf->alpha.occupancy.data + i * num_occ;
		orbital_identify_homo_lumo(scf->alpha.homo_idx + i, scf->alpha.lumo_idx + i, occ_data, num_occ);
	}

	if (scf->type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
		/*
		md_array_resize(scf->beta.coefficients.data, dim[0] * dim[1], arena);
		MEMCPY(scf->beta.coefficients.size, dim, sizeof(dim));

		md_array_resize(scf->beta.energy.data, dim[1], arena);
		scf->beta.energy.size = dim[1];

		md_array_resize(scf->beta.occupancy.data, dim[1], arena);
		scf->beta.occupancy.size = dim[1];

		md_array_resize(scf->beta.density.data, den_dim[0] * den_dim[0], arena);
        MEMCPY(scf->beta.density.size, den_dim, sizeof(den_dim));
		*/

		// Extract beta data
		if (!h5_extract_ndarray_f64(&scf->beta.coefficients, handle, "C_beta", arena)) {
			return false;
		}
		if (!h5_extract_ndarray_f64(&scf->beta.energy, handle, "E_beta", arena)) {
			return false;
		}
		if (!h5_extract_ndarray_f64(&scf->beta.occupancy, handle, "occ_beta", arena)) {
			return false;
		}
        if (!h5_extract_ndarray_f64(&scf->beta.density, handle, "D_beta", arena)) {
            return false;
        }
	} else {
		// Shallow copy fields from Alpha
		MEMCPY(&scf->beta, &scf->alpha, sizeof(md_vlx_orbital_t));
		if (scf->type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
			scf->beta.occupancy.data = 0;
            MEMSET(&scf->beta.occupancy, 0, sizeof(scf->beta.occupancy));
			if (!h5_extract_ndarray_f64(&scf->beta.occupancy, handle, "occ_beta", arena)) {
				return false;
			}

			scf->beta.homo_idx = 0;
			scf->beta.lumo_idx = 0;
			md_array_resize(scf->beta.homo_idx, num_frames, arena);
			md_array_resize(scf->beta.lumo_idx, num_frames, arena);

			// Identify homo and lumo
			for (size_t i = 0; i < num_frames; ++i) {
				// Create slices for occupancy data
				size_t num_occ = scf->beta.occupancy.size[scf->alpha.occupancy.rank - 1];
				const double* occ_data = scf->beta.occupancy.data + i * num_occ;
				orbital_identify_homo_lumo(scf->beta.homo_idx + i, scf->beta.lumo_idx + i, occ_data, num_occ);
			}
		}
	}



	// S matrix is overlap (notice dimension is the same as D)
	//md_array_resize(scf->S.data, den_dim[0] * den_dim[1], arena);
    //MEMCPY(scf->S.size, den_dim, sizeof(den_dim));

	if (!h5_extract_ndarray_f64(&scf->S, handle, "S", arena)) {
		return false;
	}

	// The ground state dipole moment is not present in all versions
	if (!h5_extract_ndarray_f64(&scf->ground_state_dipole_moment, handle, "dipole_moment", arena)) {
		//return false;
	}

    // Extract SCF history data
	if (!h5_extract_ndarray_f64(&scf->history.density_diff, handle, "scf_history_diff_density", arena)) {
		return false;
	}
	if (!h5_extract_ndarray_f64(&scf->history.energy_diff, handle, "scf_history_diff_energy", arena)) {
		return false;
	}
	if (!h5_extract_ndarray_f64(&scf->history.energy, handle, "scf_history_energy", arena)) {
		return false;
	}
	if (!h5_extract_ndarray_f64(&scf->history.gradient_norm, handle, "scf_history_gradient_norm", arena)) {
		return false;
	}
	if (!h5_extract_ndarray_f64(&scf->history.max_gradient, handle, "scf_history_max_gradient", arena)) {
		return false;
	}

	if (!h5_extract_ndarray_f64(&scf->resp_charges, handle, "charges_resp", arena)) {

	}
	/*
	
	
	{
		size_t charge_resp_dim;
		if (h5_read_dataset_dims(&charge_resp_dim, 1, handle, "charges_resp")) {
			md_array_resize(scf->resp_charges, charge_resp_dim, arena);
			if (!h5_read_dataset_data(scf->resp_charges, md_array_size(scf->resp_charges), handle, H5T_NATIVE_DOUBLE, "charges_resp")) {
				MD_LOG_ERROR("Could not read charges_resp");
				return false;
			}
		}
	}
	*/

	return true;
}

static bool h5_read_nto_data(md_vlx_rsp_t* rsp, hid_t handle, md_allocator_i* arena) {
	ASSERT(rsp);
	char buf[64];
	size_t dim[2];

	// Test to read first NTO entry to get dimensions, NTO is optional and may not exist
	if (!h5_read_dataset_dims(dim, 2, handle, "NTO_S1_alpha_orbitals")) {
		// No NTO data present, do not fail here
		return true;
	}

	md_array_resize(rsp->nto, rsp->number_of_excited_states, arena);
	MEMSET(rsp->nto, 0, rsp->number_of_excited_states * sizeof(md_vlx_orbital_t));

	for (size_t i = 0; i < rsp->number_of_excited_states; ++i) {
		int idx = (int)(i + 1);

		/*
		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_orbitals", idx);
		if (!h5_read_dataset_dims(dim, 2, handle, buf)) {
			return false;
		}
		if (dim[0] == 0 || dim[1] == 0) {
			MD_LOG_ERROR("Invalid dimensions in NTO orbitals");
			return false;
		}
		*/

		/*
		md_array_resize(rsp->nto[i].coefficients.data, dim[0] * dim[1], arena);
		MEMSET(rsp->nto[i].coefficients.data, 0, md_array_bytes(rsp->nto[i].coefficients.data));
		MEMCPY(rsp->nto[i].coefficients.size, dim, sizeof(dim));

		md_array_resize(rsp->nto[i].energy.data, dim[1], arena);
		MEMSET(rsp->nto[i].energy.data, 0, md_array_bytes(rsp->nto[i].energy.data));
		rsp->nto[i].energy.size = dim[1];

		md_array_resize(rsp->nto[i].occupancy.data, dim[1], arena);
		MEMSET(rsp->nto[i].occupancy.data, 0, md_array_bytes(rsp->nto[i].occupancy.data));
		rsp->nto[i].occupancy.size = dim[1];

		if (!h5_read_dataset_data(rsp->nto[i].coefficients.data, md_array_size(rsp->nto[i].coefficients.data), handle, H5T_NATIVE_DOUBLE, buf)) {
			return false;
		}

		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_occupations", idx);
		if (!h5_read_dataset_data(rsp->nto[i].occupancy.data, md_array_size(rsp->nto[i].occupancy.data), handle, H5T_NATIVE_DOUBLE, buf)) {
			return false;
		}

		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_energies", idx);
		if (!h5_read_dataset_data(rsp->nto[i].energy.data, md_array_size(rsp->nto[i].energy.data), handle, H5T_NATIVE_DOUBLE, buf)) {
			return false;
		}

		orbital_identify_homo_lumo(&rsp->nto[i]);
		// The NTO data is a bit shoehorned into the orbital representation and the homo lumo indices do not have a clear equivalence.
		// The occupancy information is instead the eigenvalue weights (lambdas) and are centered around the set values.
		rsp->nto[i].homo_idx = rsp->nto[i].homo_idx / 2;
		rsp->nto[i].lumo_idx = rsp->nto[i].homo_idx + 1;
		*/
	}

	return true;
}

static bool h5_read_rsp_data(md_vlx_rsp_t* rsp, hid_t handle, md_allocator_i* arena) {
	h5_read_scalar(&rsp->number_of_excited_states, handle, H5T_NATIVE_HSIZE, "number_of_states");
	if (rsp->number_of_excited_states > 0) {
		// Standard Linear Response data, allocate and read

		// Allocate data
		md_array_resize(rsp->electric_transition_dipoles, rsp->number_of_excited_states, arena);
		MEMSET(rsp->electric_transition_dipoles, 0, rsp->number_of_excited_states * sizeof(dvec3_t));

		md_array_resize(rsp->magnetic_transition_dipoles, rsp->number_of_excited_states, arena);
		MEMSET(rsp->magnetic_transition_dipoles, 0, rsp->number_of_excited_states * sizeof(dvec3_t));

		md_array_resize(rsp->velocity_transition_dipoles, rsp->number_of_excited_states, arena);
		MEMSET(rsp->velocity_transition_dipoles, 0, rsp->number_of_excited_states * sizeof(dvec3_t));

		md_array_resize(rsp->absorption_ev,	rsp->number_of_excited_states, arena);
		MEMSET(rsp->absorption_ev, 0, rsp->number_of_excited_states * sizeof(double));

		md_array_resize(rsp->oscillator_strengths, rsp->number_of_excited_states, arena);
		MEMSET(rsp->oscillator_strengths, 0, rsp->number_of_excited_states * sizeof(double));

		md_array_resize(rsp->rotatory_strengths, rsp->number_of_excited_states, arena);
		MEMSET(rsp->rotatory_strengths, 0, rsp->number_of_excited_states * sizeof(double));

		// NTO data
		if (h5_check_dataset_exists(handle, "NTO_S1_alpha_orbitals")) {
			if (!h5_read_nto_data(rsp, handle, arena)) {
				return false;
			}
		} else {
			// Check for 'nto' folder inside
            hid_t nto_group = H5Gopen(handle, "nto", H5P_DEFAULT);
			if (nto_group >= 0) {
				bool result = h5_read_nto_data(rsp, nto_group, arena);
				H5Gclose(nto_group);
                if (!result) {
                    return false;
				}
            }
		}

		// Dipoles
		size_t num_dipole_points = rsp->number_of_excited_states * 3;
		if (!h5_read_dataset_data(rsp->electric_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "electric_transition_dipoles")) {
			return false;
		}
		if (!h5_read_dataset_data(rsp->magnetic_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "magnetic_transition_dipoles")) {
			return false;
		}
		if (!h5_read_dataset_data(rsp->velocity_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "velocity_transition_dipoles")) {
			return false;
		}

		// Abs, rot and osc
		if (!h5_read_dataset_data(rsp->absorption_ev, md_array_size(rsp->absorption_ev), handle, H5T_NATIVE_DOUBLE, "eigenvalues")) {
			return false;
		}
		if (!h5_read_dataset_data(rsp->oscillator_strengths, md_array_size(rsp->oscillator_strengths), handle, H5T_NATIVE_DOUBLE, "oscillator_strengths")) {
			return false;
		}
		if (!h5_read_dataset_data(rsp->rotatory_strengths, md_array_size(rsp->rotatory_strengths), handle, H5T_NATIVE_DOUBLE, "rotatory_strengths")) {
			return false;
		}

		// Convert Atomic units (Hartree) to eV
		if (rsp->absorption_ev) {
			for (size_t i = 0; i < rsp->number_of_excited_states; ++i) {
				rsp->absorption_ev[i] *= HARTREE_TO_EV;
			}
		}
	}

	// CPP data is optional, only read if present
	size_t dim;
	if (h5_read_dataset_dims(&dim, 1, handle, "frequencies")) {
		rsp->cpp.number_of_frequencies = dim;
		md_array_resize(rsp->cpp.frequencies, dim, arena);
		if (!h5_read_dataset_data(rsp->cpp.frequencies, md_array_size(rsp->cpp.frequencies), handle, H5T_NATIVE_DOUBLE, "frequencies")) {
			// Frequencies has to be present if cpp section is present, fail if missing
			return false;
		}
		
		if (h5_check_dataset_exists(handle, "sigma")) {
			md_array_resize(rsp->cpp.sigmas, dim, arena);
			if (!h5_read_dataset_data(rsp->cpp.sigmas, md_array_size(rsp->cpp.sigmas), handle, H5T_NATIVE_DOUBLE, "sigma")) {
				return false;
			}
		}

		if (h5_check_dataset_exists(handle, "delta-epsilon")) {
			md_array_resize(rsp->cpp.delta_epsilon, dim, arena);
			if (!h5_read_dataset_data(rsp->cpp.delta_epsilon, md_array_size(rsp->cpp.delta_epsilon), handle, H5T_NATIVE_DOUBLE, "delta-epsilon")) {
				return false;
			}
		}
	}

	return true;
}

static bool h5_read_vib_data(md_vlx_vib_t* vib, hid_t handle, size_t number_of_atoms, md_allocator_i* arena) {
	size_t number_of_modes = 0;

	// Attempt to read number_of_modes (Available in new format)
	if (!h5_read_scalar(&number_of_modes, handle, H5T_NATIVE_HSIZE, "number_of_modes")) {
		// Fallback (Old format, read force_constant dims to get number of modes)
		size_t dim[2];
		int num_dim = h5_read_dataset_dims(dim, 2, handle, "force_constants");
		if (num_dim <= 0) {
			return false;
		}
		// This is a fix because the input data in one version is supplied as a 2D object
		number_of_modes = (num_dim == 1) ? dim[0] : dim[1];
	}

	if (number_of_modes == 0) {
		return false;
	}

	vib->number_of_normal_modes = number_of_modes;

	md_array_resize(vib->force_constants, number_of_modes, arena);
	if (!h5_read_dataset_data(vib->force_constants, md_array_size(vib->force_constants), handle, H5T_NATIVE_DOUBLE, "force_constants")) {
		return false;
	}

	md_array_resize(vib->ir_intensities, number_of_modes, arena);
	if (!h5_read_dataset_data(vib->ir_intensities, md_array_size(vib->ir_intensities), handle, H5T_NATIVE_DOUBLE, "ir_intensities")) {
		return false;
	}

	md_array_resize(vib->frequencies, number_of_modes, arena);
	if (!h5_read_dataset_data(vib->frequencies, md_array_size(vib->frequencies), handle, H5T_NATIVE_DOUBLE, "vib_frequencies")) {
		return false;
	}

	md_array_resize(vib->reduced_masses, number_of_modes, arena);
	if (!h5_read_dataset_data(vib->reduced_masses, md_array_size(vib->reduced_masses), handle, H5T_NATIVE_DOUBLE, "reduced_masses")) {
		return false;
	}

	// Check if "normal_modes" is a group or dataset
	hid_t obj_info = H5Oopen(handle, "normal_modes", H5P_DEFAULT);
	if (obj_info == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open 'normal_modes' object");
		return false;
	}

	// Read normal modes (group or dataset)
	hid_t obj_type = H5Iget_type(obj_info);
	if (obj_type == H5I_GROUP) {
        hid_t normal_modes_id = H5Gopen(handle, "normal_modes", H5P_DEFAULT);
        if (normal_modes_id != H5I_INVALID_HID) {
			if (number_of_atoms == 0) {
				MD_LOG_ERROR("Missing number of atoms, is required for normal modes");
				return false;
			}

            char lbl[32];
            for (size_t i = 0; i < vib->number_of_normal_modes; ++i) {
                snprintf(lbl, sizeof(lbl), "%zu", i + 1);

                dvec3_t* data = md_array_create(dvec3_t, number_of_atoms, arena);
                MEMSET(data, 0, sizeof(dvec3_t) * number_of_atoms);
                if (!h5_read_dataset_data(data, 3 * number_of_atoms, normal_modes_id, H5T_NATIVE_DOUBLE, lbl)) {
                    MD_LOG_ERROR("Failed to extract dataset in '%s' normal mode", lbl);
                    md_array_free(data, arena);
                    return false;
                }

                // Success, append ata
                md_array_push(vib->normal_modes, data, arena);
            }
        }
	} else if (obj_type == H5I_DATASET) {
		// Handle dataset case
		// Iterate over outer dimension in dataset, which should be [number_of_normal_modes][number_of_atoms][3]

		size_t data_dim[3];
		int num_dim = h5_read_dataset_dims(data_dim, 3, handle, "normal_modes");

		// Assert expected dimensions
		if (num_dim != 3 || data_dim[0] != vib->number_of_normal_modes || data_dim[1] != number_of_atoms || data_dim[2] != 3) {
			MD_LOG_ERROR("Unexpected dimensions in normal_modes dataset");
			H5Oclose(obj_info);
			return false;
		}

		size_t num_points = data_dim[0] * data_dim[1] * data_dim[2];
		double* raw_data = md_array_create(double, num_points, arena);
		if (!h5_read_dataset_data(raw_data, num_points, handle, H5T_NATIVE_DOUBLE, "normal_modes")) {
			MD_LOG_ERROR("Failed to read normal_modes dataset");
			md_array_free(raw_data, arena);
			H5Oclose(obj_info);
			return false;
		}

		// Set the pointers to each normal mode (within raw_data)
		dvec3_t* base_ptr = (dvec3_t*)raw_data;
		for (size_t i = 0; i < vib->number_of_normal_modes; ++i) {
			dvec3_t* mode_data = base_ptr + (i * number_of_atoms);
			md_array_push(vib->normal_modes, mode_data, arena);
		}
	} else {
		MD_LOG_ERROR("Unrecognized object type for 'normal_modes'");
		H5Oclose(obj_info);
		return false;
	}

	return true;
}

static bool h5_read_opt_data(md_vlx_opt_t* opt, hid_t handle, size_t number_of_atoms, md_allocator_i* arena) {
	// @TODO(This will likely be exposed as its own variable in the future, for now we extract the length from one of the fields)
	size_t dim[3];
	int num_dim;
	
	num_dim = h5_read_dataset_dims(dim, 3, handle, "nuclear_repulsion_energies");
	if (num_dim <= 0) {
		return false;
	}

	// This is a fix because the input data in one version is supplied as a 2D object
	size_t len = dim[0];
	opt->number_of_steps = len;

	md_array_resize(opt->nuclear_repulsion_energies, len, arena);
	if (!h5_read_dataset_data(opt->nuclear_repulsion_energies, md_array_size(opt->nuclear_repulsion_energies), handle, H5T_NATIVE_DOUBLE, "nuclear_repulsion_energies")) {
		return false;
	}

	const char* energy_ident = "";
	const char* coord_ident  = "";

	if (h5_check_dataset_exists(handle, "opt_energies")) {
		energy_ident = "opt_energies";
		coord_ident  = "opt_coordinates_au";
	}
	else if (h5_check_dataset_exists(handle, "scan_energies")) {
		energy_ident = "scan_energies";
		coord_ident  = "scan_coordinates_au";
	}

	if (strlen(energy_ident) == 0 || strlen(coord_ident) == 0) {
		MD_LOG_ERROR("No optimization or scan data found in HDF5 file");
		return false;
	}

	md_array_resize(opt->energies, len, arena);
	if (!h5_read_dataset_data(opt->energies, md_array_size(opt->energies), handle, H5T_NATIVE_DOUBLE, energy_ident)) {
		return false;
	}

	num_dim = h5_read_dataset_dims(dim, 3, handle, coord_ident);
	if (dim[0] != len || dim[1] != number_of_atoms || dim[2] != 3) {
		MD_LOG_ERROR("Inconsistent or invalid opt_coordinates dimensions");
		return false;
	}

	md_array_resize(opt->coordinates, dim[0] * dim[1], arena);
	if (!h5_read_dataset_data(opt->coordinates, md_array_size(opt->coordinates) * 3, handle, H5T_NATIVE_DOUBLE, coord_ident)) {
		return false;
	}

	if (opt->coordinates) {
		for (size_t i = 0; i < dim[0] * dim[1]; ++i) {
			opt->coordinates[i] = dvec3_mul_f(opt->coordinates[i], BOHR_TO_ANGSTROM);
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

	if (!h5_read_scalar(&vlx->nuclear_repulsion_energy, handle, H5T_NATIVE_DOUBLE, "nuclear_repulsion")) {
		return false;
	}

	size_t dims[8] = {0};
	if (!h5_read_dataset_dims(dims, ARRAY_SIZE(dims), handle, "nuclear_repulsion")) {
		return false;
	}

	vlx->number_of_frames = dims[0] > 0 ? dims[0] : 1;

	if (!h5_read_scalar(&vlx->number_of_alpha_electrons, handle, H5T_NATIVE_UINT64, "number_of_alpha_electrons")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_atoms, handle, H5T_NATIVE_UINT64, "number_of_atoms")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_beta_electrons, handle, H5T_NATIVE_UINT64, "number_of_beta_electrons")) {
		return false;
	}

	if (!h5_read_str(&vlx->potfile_text, handle, "potfile_text", vlx->arena)) {
		return false;
	}

	if (!h5_read_scalar(&vlx->spin_multiplicity, handle, H5T_NATIVE_INT32, "spin_multiplicity")) {
		return false;
	}

	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("Number of atoms is zero");
		return false;
	}

    if (!h5_extract_ndarray_f64(&vlx->atom_coordinates, handle, "atom_coordinates", vlx->arena)) {
		return false;
	}
    // We expect coordinates to be rank 2 or 3, and have inner most dimension of 3, and outer dimension matching number of atoms
	if ((vlx->atom_coordinates.rank != 2 && vlx->atom_coordinates.rank != 3) ||
		vlx->atom_coordinates.size[vlx->atom_coordinates.rank-1] != 3 ||
		vlx->atom_coordinates.size[vlx->atom_coordinates.rank-2] != vlx->number_of_atoms) {
		MD_LOG_ERROR("Inconsistent or invalid atom_coordinates dimensions");
		return false;
	}

	if (vlx->number_of_frames > 1) {
		if (vlx->atom_coordinates.rank != 3 ||
			vlx->atom_coordinates.size[0] != vlx->number_of_frames)
		{
            MD_LOG_ERROR("Inconsistent or invalid atom_coordinates dimensions for multiple frames");
			return false;
		}
	}

	// Convert Atomic units to Ångström
	if (vlx->atom_coordinates.data) {
		size_t tot_len = md_ndarray_f64_size(&vlx->atom_coordinates);
		for (size_t i = 0; i < tot_len; ++i) {
            vlx->atom_coordinates.data[i] *= BOHR_TO_ANGSTROM;
		}
	}

	md_array_resize(vlx->atomic_numbers, vlx->number_of_atoms, vlx->arena);
	MEMSET(vlx->atomic_numbers, 0, md_array_bytes(vlx->atomic_numbers));
	if (!h5_read_dataset_data(vlx->atomic_numbers, md_array_size(vlx->atomic_numbers), handle, H5T_NATIVE_UINT8, "nuclear_charges")) {
		return false;
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
		if (!h5_read_scf_data(&vlx->scf, file_id, vlx->arena)) {
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
				result = h5_read_scf_data(&vlx->scf, scf_id, vlx->arena);
				H5Gclose(scf_id);
				if (!result) goto done;
			}
		}
	}

	// VIB
    if (flags & VLX_FLAG_VIB) {
        if (H5Lexists(file_id, "vib", H5P_DEFAULT) > 0) {
            hid_t vib_id = H5Gopen(file_id, "vib", H5P_DEFAULT);
            if (vib_id != H5I_INVALID_HID) {
                result = h5_read_vib_data(&vlx->vib, vib_id, vlx->number_of_atoms, vlx->arena);
                H5Gclose(vib_id);
                if (!result) goto done;
            }
        }
    }

	// OPT
	if (flags & VLX_FLAG_OPT) {
		if (H5Lexists(file_id, "opt", H5P_DEFAULT) > 0) {
			hid_t opt_id = H5Gopen(file_id, "opt", H5P_DEFAULT);
			if (opt_id != H5I_INVALID_HID) {
				result = h5_read_opt_data(&vlx->opt, opt_id, vlx->number_of_atoms, vlx->arena);
				H5Gclose(opt_id);
				if (!result) goto done;
			}
		}
	}

	// RSP
	if (flags & VLX_FLAG_RSP) {
		if (H5Lexists(file_id, "rsp", H5P_DEFAULT) > 0) {
			hid_t rsp_id = H5Gopen(file_id, "rsp", H5P_DEFAULT);
			if (rsp_id != H5I_INVALID_HID) {
				result = h5_read_rsp_data(&vlx->rsp, rsp_id, vlx->arena);
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

#define BAKE_STR(str) {str "", sizeof(str) - 1}

static inline str_t resolve_basis_set_ident(str_t input) {
	struct map_t {
        str_t in;
        str_t out;
	};

	static const struct map_t alias_table[] = {
        {BAKE_STR("6-31G*"),			BAKE_STR("6-31G_D_")},
		{BAKE_STR("6-31G**"),			BAKE_STR("6-31G_D,P_")},
		{BAKE_STR("6-31+G*"),			BAKE_STR("6-31+G_D_")},
        {BAKE_STR("6-31+G**"),			BAKE_STR("6-31+G_D,P_")},
		{BAKE_STR("6-31++G*"),			BAKE_STR("6-31++G_D_")},
		{BAKE_STR("6-31++G**"),			BAKE_STR("6-31++G_D,P_")},
        {BAKE_STR("6-311G*"),			BAKE_STR("6-311G_D_")},
		{BAKE_STR("6-311G**"),			BAKE_STR("6-311G_D,P_")},
		{BAKE_STR("6-311+G*"),			BAKE_STR("6-311+G_D_")},
        {BAKE_STR("6-311+G**"),			BAKE_STR("6-311+G_D,P_")},
		{BAKE_STR("6-311++G*"),			BAKE_STR("6-311++G_D_")},
		{BAKE_STR("6-311++G**"),		BAKE_STR("6-311++G_D,P_")},
        {BAKE_STR("6-31G(2DF,P)"),		BAKE_STR("6-31G_2DF,P_")},
		{BAKE_STR("6-31G(3DF,3PD)"),	BAKE_STR("6-31G_3DF,3PD_")},
		{BAKE_STR("6-311G(2DF,2PD)"),	BAKE_STR("6-311G_2DF,2PD_")},
        {BAKE_STR("6-311+G(2D,P)"),		BAKE_STR("6-311+G_2D,P_")},
		{BAKE_STR("6-311++G(2D,2P)"),	BAKE_STR("6-311++G_2D,2P_")},
		{BAKE_STR("6-311++G(3DF,3PD)"),	BAKE_STR("6-311++G_3DF,3PD_")},
        {BAKE_STR("DEF2-SV(P)"),		BAKE_STR("DEF2-SV_P_")},
    };

	for (size_t i = 0; i < ARRAY_SIZE(alias_table); ++i) {
        if (str_eq(input, alias_table[i].in)) {
            return alias_table[i].out;
		}
	}

	return input;
}

#undef BAKE_STR

// Internal version to control what portions to load
static bool vlx_parse_file(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	size_t temp_pos = md_temp_get_pos();

	bool result = false;

	if (str_ends_with(filename, STR_LIT(".scf.results.h5"))) {
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

		str_t ident = resolve_basis_set_ident(vlx->basis_set_ident);
		MD_LOG_DEBUG("Basis set ident: '" STR_FMT "'", STR_ARG(ident));

		char exe_buf[1024];
		str_t exe_path = {exe_buf, md_path_write_exe(exe_buf, sizeof(exe_buf))};

		str_t exe_dir = {0};
		if (!extract_folder_path(&exe_dir, exe_path)) {
			MD_LOG_ERROR("Failed to extract executable directory");
		}

		md_strb_fmt(&sb, STR_FMT "%s/" STR_FMT, STR_ARG(exe_dir), MD_VLX_BASIS_FOLDER, STR_ARG(ident));
		str_t basis_filepath = md_strb_to_str(sb);
		md_file_t basis_file = {0};
		if (md_file_open(&basis_file, basis_filepath, MD_FILE_READ)) {
			MD_LOG_DEBUG("Attempting to parse VLX basis set from file: '" STR_FMT "'", STR_ARG(basis_filepath));
			md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, basis_file);
			bool parse_result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
			md_file_close(&basis_file);
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
			md_strb_push_str(&sb, ident);
			basis_filepath = md_strb_to_str(sb);
			if (md_file_open(&basis_file, basis_filepath, MD_FILE_READ)) {
				MD_LOG_DEBUG("Attempting to parse VLX basis set from file: '" STR_FMT "'", STR_ARG(basis_filepath));
				md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, basis_file);
				bool parse_result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
				md_file_close(&basis_file);
				if (!parse_result) {
					MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
					goto done;
				}
				normalize_basis_set(&vlx->basis_set);
			}
			else {
                MD_LOG_ERROR("Could not find basis file corresponding to identifier: '" STR_FMT "'", STR_ARG(ident));
				goto done;
			}
		}
	}

	if (vlx->number_of_atoms > 0 && vlx->scf.type != MD_VLX_SCF_TYPE_UNKNOWN) {
		// Extract ao_to_atom_idx map
		size_t N = md_vlx_scf_number_of_atomic_orbitals(vlx);
		if (N > 0) {
			md_array_resize(vlx->ao_to_atom_idx, N, vlx->arena);
			extract_ao_to_atom_idx(vlx->ao_to_atom_idx, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set);
		}

		const dvec3_t* coords = (const dvec3_t*)vlx->atom_coordinates.data;
		extract_gto_data(&vlx->gto_data, coords, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set, vlx->arena);
	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);

	return result;
}

// Extract Natural Transition Orbitals PGTOs
size_t md_vlx_nto_gto_count(const md_vlx_t* vlx) {
	return vlx->gto_data.num_pgtos;
}

static inline void extract_row(double* dst, const md_2darray_f64_t* data, size_t row_idx) {
	ASSERT(dst);
	ASSERT(data);
	ASSERT(row_idx < data->size[0]);

	size_t num_cols = data->size[1];
	for (size_t i = 0; i < num_cols; ++i) {
		dst[i] = data->data[row_idx * num_cols + i];
	}
}

static inline void extract_col(double* dst, const md_2darray_f64_t* data, size_t col_idx) {
	ASSERT(dst);
	ASSERT(data);
	ASSERT(col_idx < data->size[1]);

	for (size_t i = 0; i < data->size[0]; ++i) {
		dst[i] = data->data[i * data->size[1] + col_idx];
	}
}

static inline size_t orb_number_of_snapshots(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
    if (orb->coefficients.rank == 2) {
		return 1; // No snapshot dimension, treat as single snapshot
	} else if (orb->coefficients.rank == 3) {
		return orb->coefficients.size[0];
	}
	// Error
	MD_LOG_ERROR("malformed orb object");
	return 0;
}

static inline size_t number_of_molecular_orbitals(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[orb->coefficients.rank-2];
}

static inline size_t number_of_mo_coefficients(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[orb->coefficients.rank-1];
}

static inline md_2darray_f64_t create_2d_slice(md_ndarray_f64_t ndarray, size_t slice_idx) {
	ASSERT(ndarray.rank >= 2);
	md_2darray_f64_t arr = {
		.data = ndarray.data + slice_idx * ndarray.size[ndarray.rank - 2] * ndarray.size[ndarray.rank - 1],
		.size = {ndarray.size[ndarray.rank - 1], ndarray.size[ndarray.rank - 2]}
	};
	return arr;
}

static inline void extract_mo_coefficients(double* out_coeff, const md_vlx_orbital_t* orb, size_t snapshot_idx, size_t mo_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);
	ASSERT(mo_idx < number_of_molecular_orbitals(orb));
    ASSERT(snapshot_idx < orb_number_of_snapshots(orb));

    md_2darray_f64_t arr = create_2d_slice(orb->coefficients, snapshot_idx);
	extract_col(out_coeff, &arr, mo_idx);
}

static inline size_t number_of_atomic_orbitals(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[orb->coefficients.rank-1];
}

static inline size_t number_of_ao_coefficients(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[orb->coefficients.rank-2];
}

static inline void extract_ao_coefficients(double* out_coeff, const md_vlx_orbital_t* orb, size_t ao_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);
	ASSERT(ao_idx < number_of_atomic_orbitals(orb));

    md_2darray_f64_t arr = create_2d_slice(orb->coefficients, 0);
	extract_row(out_coeff, &arr, ao_idx);
}	

static size_t extract_gtos(md_gto_t* out_gtos, const md_gto_data_t* ao_data, const double* mo_coeffs, double value_cutoff) {
	size_t count = 0;
	for (size_t i = 0; i < ao_data->num_cgtos; ++i) {
		for (size_t j = ao_data->cgto_offset[i]; j < ao_data->cgto_offset[i+1]; ++j) {
			int pi,pj,pk,pl;
			md_gto_unpack_ijkl(ao_data->pgto_ijkl[j], &pi, &pj, &pk, &pl);
			double radius = md_gto_compute_radius_of_influence(pi, pj, pk, ao_data->pgto_coeff[j], ao_data->pgto_alpha[j], value_cutoff);
			if (radius == 0.0) {
				continue; // Skip GTOs with zero radius of influence given this cutoff_value
			}

			out_gtos[count].x = ao_data->cgto_xyzr[i].x;
			out_gtos[count].y = ao_data->cgto_xyzr[i].y;
			out_gtos[count].z = ao_data->cgto_xyzr[i].z;
			out_gtos[count].coeff = (float)(mo_coeffs[i] * ao_data->pgto_coeff[j]);
			out_gtos[count].alpha = ao_data->pgto_alpha[j];
			out_gtos[count].cutoff = (float)radius;
			out_gtos[count].i = (uint8_t)pi;
			out_gtos[count].j = (uint8_t)pj;
			out_gtos[count].k = (uint8_t)pk;
			out_gtos[count].l = (uint8_t)pl;
			out_gtos[count]._pad = (uint32_t)i; // Store mo coeff idx here
			count += 1;
		}
	}
	return count;
}

bool md_vlx_nto_gto_extract(md_gto_t* out_gtos, const md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type) {
	ASSERT(out_gtos);
	ASSERT(vlx);

	if (nto_idx >= vlx->rsp.number_of_excited_states) {
		MD_LOG_ERROR("Invalid nto index!");
		return false;
	}

	if (!vlx->rsp.nto) {
		MD_LOG_ERROR("Veloxchem data is missing NTO data");
		return false;
	}

	const md_vlx_orbital_t* orb = &vlx->rsp.nto[nto_idx];

	// The lambda values are stored symmetrically around homo/lumo
	size_t max_lambda_idx = orb->homo_idx[0];

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
		mo_idx = (int64_t)orb->lumo_idx[0] + (int64_t)lambda_idx;
	} else if (type == MD_VLX_NTO_TYPE_HOLE) {
		mo_idx = (int64_t)orb->homo_idx[0] - (int64_t)lambda_idx;
	} else {
		MD_LOG_ERROR("Invalid NTO type!");
		return false;
	}

	size_t temp_pos = md_temp_get_pos();
	size_t num_mo_coeffs = number_of_mo_coefficients(orb);
	double* mo_coeffs = md_temp_push(sizeof(double) * num_mo_coeffs);
	const dvec3_t* atom_xyz = (const dvec3_t*)vlx->atom_coordinates.data;

	extract_mo_coefficients(mo_coeffs, orb, 0, mo_idx);
	extract_gtos(out_gtos, &vlx->gto_data, mo_coeffs, 0.0);


	md_temp_set_pos_back(temp_pos);
	return true;
}

size_t md_vlx_mo_gto_count(const md_vlx_t* vlx) {
	ASSERT(vlx);
	return vlx->gto_data.num_pgtos;
}

// Attempts to compute fitting volume dimensions given an input extent and a suggested number of samples per length unit
static inline void compute_dim(int out_dim[3], vec3_t in_ext, float samples_per_unit_length) {
	out_dim[0] = CLAMP(ALIGN_TO((int)(in_ext.x * samples_per_unit_length), 8), 8, 512);
	out_dim[1] = CLAMP(ALIGN_TO((int)(in_ext.y * samples_per_unit_length), 8), 8, 512);
	out_dim[2] = CLAMP(ALIGN_TO((int)(in_ext.z * samples_per_unit_length), 8), 8, 512);
}

static inline vec3_t closest_point_in_aabb(vec3_t p, vec3_t aabb_min, vec3_t aabb_max) {
    return vec3_clamp(p, aabb_min, aabb_max);
}

static inline float monomial_bound(vec3_t center, int powers[3], vec3_t aabb_min, vec3_t aabb_max) {
    vec3_t bounds = {0};
    bounds.x = MAX(fabsf(powf(aabb_min.x - center.x, (float)powers[0])), fabsf(powf(aabb_max.x - center.x, (float)powers[0])));
    bounds.y = MAX(fabsf(powf(aabb_min.y - center.y, (float)powers[1])), fabsf(powf(aabb_max.y - center.y, (float)powers[1])));
    bounds.z = MAX(fabsf(powf(aabb_min.z - center.z, (float)powers[2])), fabsf(powf(aabb_max.z - center.z, (float)powers[2])));
    return bounds.x * bounds.y * bounds.z;
}

static inline void gaussian_product_center(vec3_t* out_P, float* out_gamma, float alpha1, vec3_t A, float alpha2, vec3_t B) {
    *out_gamma = alpha1 + alpha2;
    *out_P = vec3_div_f(vec3_add(vec3_mul_f(A, alpha1), vec3_mul_f(B, alpha2)), *out_gamma);
}

# if 0
static double estimate_contracted_pair_bound(vec3_t A, const struct pgto_t* ao1, size_t num_ao1, vec3_t B, const struct pgto_t* ao2, size_t num_ao2, double D_mu_nu, vec3_t aabb_min, vec3_t aabb_max) {
    double max_bound = 0.0;
    for (size_t i = 0; i < num_ao1; ++i) {
    	float alpha1 = ao1[i].alpha;
    	float coeff1 = ao1[i].coeff;

        for (size_t j = 0; j < num_ao2; ++j) {
	    	float alpha2 = ao2[j].alpha;
	    	float coeff2 = ao2[j].coeff;

            vec3_t P;
            float gamma;
            gaussian_product_center(&P, &gamma, alpha1, A, alpha2, B);

            //printf("P: %f %f %f, gamma: %f\n", P.x, P.y, P.z, gamma);

            vec3_t r_min = closest_point_in_aabb(P, aabb_min, aabb_max);
            float dist2 = vec3_distance_squared(r_min, P);
            double exp_bound = exp(-gamma * dist2);

            // Polynomial Bound
            int powers[3] = {
            	ao1[i].i + ao2[j].i,
            	ao1[i].j + ao2[j].j,
            	ao1[i].k + ao2[j].k,
            };
            double poly_bound = monomial_bound(P, powers, aabb_min, aabb_max);

            // Primitive pair contribution
            double prim_bound = fabs(coeff1 * coeff2) * poly_bound * exp_bound;
            max_bound += prim_bound;
        }
    }

    return fabs(D_mu_nu) * max_bound;
}
#endif

size_t md_vlx_mo_gto_extract(md_gto_t gtos[], const md_vlx_t* vlx, size_t snapshot_idx, size_t mo_idx, md_vlx_mo_type_t type, double value_cutoff) {
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

#if 0
	{
		vec3_t aabb_min = vec3_set1( FLT_MAX);
		vec3_t aabb_max = vec3_set1(-FLT_MAX);
		double cutoff = 1.0e-6;

		for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
			vec3_t coord = { vlx->atom_coordinates[i].x, vlx->atom_coordinates[i].y, vlx->atom_coordinates[i].z };
			aabb_min = vec3_min(aabb_min, coord);
			aabb_max = vec3_max(aabb_max, coord);
		}

		aabb_min = vec3_mul_f(aabb_min, ANGSTROM_TO_BOHR);
		aabb_max = vec3_mul_f(aabb_max, ANGSTROM_TO_BOHR);
		const float pad = 5.0f;
		aabb_min = vec3_sub_f(aabb_min, pad);
		aabb_max = vec3_sub_f(aabb_max, pad);
		vec3_t ext = vec3_sub(aabb_max, aabb_min);

		const float samples_per_unit_length = 4 * BOHR_TO_ANGSTROM;

		int dim[3];
		compute_dim(dim, ext, samples_per_unit_length);

		vec3_t voxel_ext = vec3_div(ext, vec3_set(dim[0], dim[1], dim[2]));

		const double* den_matrix = vlx->scf.alpha.density.data;
		size_t den_dim = vlx->scf.alpha.density.size[0];

		printf("Den_dim: %zu\n", den_dim);
		printf("dimensions of volume: %i %i %i\n", dim[0], dim[1], dim[2]);

		md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

		size_t non_zero_blocks = 0;
		size_t min_CGTOs_per_region = 1024;
		size_t max_CGTOs_per_region = 0;

		size_t sparse_matrix_bytes = 0;

		if (den_matrix && den_dim > 0) {
			for (size_t i = 0; i < ao_data.num_cgtos; ++i) {
				printf("CGTO %zu\n", i);
				for (size_t j = ao_data.cgtos[i].pgto_offset; j < ao_data.cgtos[i].pgto_offset + ao_data.cgtos[i].pgto_count; ++j) {
					struct pgto_t pgto = ao_data.pgtos[j];
					printf("\t %i [%i %i %i]   %12.6f %12.6f\n", pgto.l, pgto.i, pgto.j, pgto.k, pgto.coeff, pgto.alpha);
				}
			}

			md_array(uint32_t) candidates = 0;
			md_array(float)  max_cgto_radius = 0;

			md_array_resize(max_cgto_radius, ao_data.num_cgtos, arena);

			for (size_t i = 0; i < ao_data.num_cgtos; ++i) {
				float max_radius = 0.0f;
				for (size_t j = ao_data.cgtos[i].pgto_offset; j < ao_data.cgtos[i].pgto_offset + ao_data.cgtos[i].pgto_count; ++j) {
					ao_data.pgtos[j].radius = md_gto_compute_radius_of_influence(ao_data.pgtos[j].i, ao_data.pgtos[j].j, ao_data.pgtos[j].k, ao_data.pgtos[j].coeff, ao_data.pgtos[j].alpha, cutoff);
					max_radius = MAX(max_radius, ao_data.pgtos[j].radius);
				}
				max_cgto_radius[i] = max_radius;
			}

			// Iterate over all 8x8x8 subregions
			for (int z = 0; z < dim[2]; z += 8) {
				for (int y = 0; y < dim[1]; y += 8) {
					for (int x = 0; x < dim[0]; x += 8) {
						vec3_t local_min = vec3_add(aabb_min, vec3_mul(vec3_set(x + 0, y + 0, z + 0), voxel_ext));
						vec3_t local_max = vec3_add(aabb_min, vec3_mul(vec3_set(x + 8, y + 8, z + 8), voxel_ext));
						size_t num_cgtos = 0;

						for (size_t i = 0; i < ao_data.num_cgtos; ++i) {
							vec3_t coord   = vec3_load(&ao_data.cgtos[i].x);
							vec3_t clamped = vec3_clamp(coord, local_min, local_max);
							float r2 = max_cgto_radius[i] * max_cgto_radius[i];
							if (vec3_distance_squared(coord, clamped) < r2) {
								num_cgtos += 1;
							}
						}

						min_CGTOs_per_region = MIN(min_CGTOs_per_region, num_cgtos);
						max_CGTOs_per_region = MAX(max_CGTOs_per_region, num_cgtos);
					}
				}
			}
		}

		size_t dense_matrix_bytes = den_dim * den_dim * sizeof(float);

		size_t total = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2] / 512;
		printf("Number of 8x8x8 blocks in total: %zu\n", total);
		printf("Number of populated 8x8x8 blocks: %zu \n", non_zero_blocks);
		printf("Min CGTOs per region: %zu \n", min_CGTOs_per_region);
		printf("Max CGTOs per region: %zu \n", max_CGTOs_per_region);

		printf("Dense matrix bytes: %zu \n", dense_matrix_bytes);
		printf("Spares matrix bytes: %zu \n", sparse_matrix_bytes);

		md_vm_arena_destroy(arena);
	}
#endif

	size_t temp_pos = md_temp_get_pos();
	size_t num_mo_coeffs = number_of_mo_coefficients(orb);
	double* mo_coeffs = md_temp_push(sizeof(double) * num_mo_coeffs);
	extract_mo_coefficients(mo_coeffs, orb, snapshot_idx, mo_idx);

	size_t count = extract_gtos(gtos, &vlx->gto_data, mo_coeffs, value_cutoff);

	md_temp_set_pos_back(temp_pos);
	return count;
}

const double* md_vlx_scf_resp_charges(const md_vlx_t* vlx) {
	if (vlx) {
		if (vlx->scf.resp_charges.data) {
			return vlx->scf.resp_charges.data;
		}
	}
	return NULL;
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

size_t md_vlx_rsp_cpp_number_of_frequencies(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.cpp.number_of_frequencies;
	}
	return 0;
}

const double* md_vlx_rsp_cpp_frequencies(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.cpp.frequencies;
	}
	return NULL;
}

const double* md_vlx_rsp_cpp_sigma(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.cpp.sigmas;
	}
	return NULL;
}

const double* md_vlx_rsp_cpp_delta_epsilon(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.cpp.delta_epsilon;
	}
	return NULL;
}

bool md_vlx_rsp_has_nto(const md_vlx_t* vlx) {
	return vlx && vlx->rsp.nto;
}

const double* md_vlx_rsp_nto_occupancy(const md_vlx_t* vlx, size_t nto_idx) {
	if (vlx) {
		if (vlx->rsp.nto && nto_idx < vlx->rsp.number_of_excited_states) {
			return vlx->rsp.nto[nto_idx].occupancy.data;
		}
	}
	return NULL;
}

const double* md_vlx_rsp_nto_lambdas(const md_vlx_t* vlx, size_t nto_idx) {
	if (vlx) {
		if (vlx->rsp.nto && nto_idx < vlx->rsp.number_of_excited_states) {
			return vlx->rsp.nto[nto_idx].occupancy.data + vlx->rsp.nto[nto_idx].lumo_idx[0];
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

size_t md_vlx_rsp_nto_number_of_ao_coefficients(const struct md_vlx_t* vlx, size_t nto_idx) {
	if (vlx) {
		if (vlx->rsp.nto && nto_idx < vlx->rsp.number_of_excited_states) {
			return vlx->rsp.nto[nto_idx].coefficients.size[0];
		}
	}
	return 0;
}

size_t md_vlx_rsp_nto_extract_coefficients(double* out_ao, const md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type) {
	if (vlx) {
		if (vlx->rsp.nto && nto_idx < vlx->rsp.number_of_excited_states) {
			const md_vlx_orbital_t* orb = &vlx->rsp.nto[nto_idx];
			// We return the data pointer to the LUMO orbital index

			int64_t mo_idx = 0;
			if (type == MD_VLX_NTO_TYPE_PARTICLE) {
				mo_idx = (int64_t)orb->lumo_idx[0] + (int64_t)lambda_idx;
			} else if (type == MD_VLX_NTO_TYPE_HOLE) {
				mo_idx = (int64_t)orb->homo_idx[0] - (int64_t)lambda_idx;
			} else {
				MD_LOG_ERROR("Invalid NTO type!");
				return 0;
			}

			if (mo_idx < 0 || mo_idx >= (int64_t)md_vlx_scf_number_of_molecular_orbitals(vlx)) {
				MD_LOG_ERROR("lambda_idx out of bounds");
				return 0;
			}

			if (out_ao) {
				extract_mo_coefficients(out_ao, orb, 0, mo_idx);
			}
			return number_of_mo_coefficients(orb);
		}
	}
	return 0;
}

size_t md_vlx_vib_number_of_normal_modes(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->vib.number_of_normal_modes;
	}
	return 0;
}

const double* md_vlx_vib_ir_intensities(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->vib.ir_intensities;
	}
	return NULL;
}

const double* md_vlx_vib_frequencies(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->vib.frequencies;
    }
    return NULL;
}

const double* md_vlx_vib_reduced_masses(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->vib.reduced_masses;
    }
    return NULL;
}

const double* md_vlx_vib_force_constants(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->vib.force_constants;
    }
    return NULL;
}

const dvec3_t* md_vlx_vib_normal_mode(const struct md_vlx_t* vlx, size_t idx) {
	if (vlx) {
		if (vlx->vib.normal_modes && idx < vlx->vib.number_of_normal_modes) {
			return vlx->vib.normal_modes[idx];
		}
	}
	return NULL;
}

// OPT
size_t md_vlx_opt_number_of_steps(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->opt.number_of_steps;
	}
	return 0;
}

// Returns atom coordinates for a given optimization step
const dvec3_t* md_vlx_opt_coordinates(const struct md_vlx_t* vlx, size_t opt_idx) {
	if (vlx) {
		if (vlx->opt.coordinates && 0 <= opt_idx && opt_idx < vlx->opt.number_of_steps) {
			const size_t stride = vlx->number_of_atoms;
			return vlx->opt.coordinates + stride * opt_idx;
		}
	}
	return NULL;
}

const double* md_vlx_opt_energies(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->opt.energies;
	}
	return NULL;
}

const double* md_vlx_opt_nuclear_repulsion_energies(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->opt.nuclear_repulsion_energies;
	}
	return NULL;
}

/*
size_t md_vlx_vib_num_raman_activity(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->vib.num_raman;
	}
    return 0;
}

const double* md_vlx_vib_raman_activity(const md_vlx_t* vlx, size_t idx) {
    if (vlx) {
        return vlx->vib.raman_activity;
	}
    return NULL;
}
*/

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

typedef struct vlx_traj_t {
    uint64_t magic; // For validation

	size_t num_frames;
	size_t num_atoms;

    float* atom_x; // [num_frames][num_atoms]
    float* atom_y; // [num_frames][num_atoms]
    float* atom_z; // [num_frames][num_atoms]

    double* frame_times; // [num_frames]

	md_allocator_i* alloc;
} vlx_traj_t;

static bool vlx_traj_get_header(const struct md_trajectory_o* inst, md_trajectory_header_t* out_header) {
	ASSERT(inst);
    ASSERT(out_header);
    vlx_traj_t* traj = (vlx_traj_t*)inst;
    if (traj && traj->magic == VLX_MAGIC) {
		out_header->num_frames  = traj->num_frames;
		out_header->num_atoms   = traj->num_atoms;
		out_header->frame_times = traj->frame_times;
		out_header->time_unit   = md_unit_none();
		return true;
	}

    return false;
}

static bool vlx_traj_reader_load_frame(struct md_trajectory_reader_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* out_header, float* out_x, float* out_y, float* out_z) {
    vlx_traj_t* traj = (vlx_traj_t*)inst;
	if (traj && traj->magic == VLX_MAGIC) {
        if (frame_idx < 0 || (size_t)frame_idx >= traj->num_frames) {
			MD_LOG_ERROR("Frame index out of bounds");
			return false;
		}
        // Extract data for the requested frame
		if (out_x) {
            MEMCPY(out_x, traj->atom_x + frame_idx * traj->num_atoms, sizeof(float) * traj->num_atoms);
		}
        if (out_y) {
			MEMCPY(out_y, traj->atom_y + frame_idx * traj->num_atoms, sizeof(float) * traj->num_atoms);
		}
        if (out_z) {
			MEMCPY(out_z, traj->atom_z + frame_idx * traj->num_atoms, sizeof(float) * traj->num_atoms);
		}
		if (out_header) {
            out_header->index = frame_idx;
            out_header->num_atoms = traj->num_atoms;
            out_header->timestamp = traj->frame_times[frame_idx];
			out_header->unitcell = md_unitcell_none();
		}
		return true;
	}
	return false;
}

static void vlx_traj_reader_free(md_trajectory_reader_i* reader) {
	// No dynamic resources to free in this example, but if there were, we'd clean them up here.
    MEMSET(reader, 0, sizeof(md_trajectory_reader_i));
}


static bool vlx_traj_reader_init(md_trajectory_reader_i* reader, struct md_trajectory_o* inst) {
    vlx_traj_t* traj = (vlx_traj_t*)inst;
    if (traj && traj->magic == VLX_MAGIC) {
		reader->inst = (struct md_trajectory_reader_o*)traj;
		reader->load_frame = vlx_traj_reader_load_frame;
        reader->free = vlx_traj_reader_free;
		return true;
	}
	return false;
}

static void vlx_traj_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    vlx_traj_t* vlx_traj = (vlx_traj_t*)traj->inst;
    if (vlx_traj->magic != VLX_MAGIC) {
        MD_LOG_ERROR("VLX: Cannot free trajectory, is not a valid VLX trajectory.");
        ASSERT(false);
        return;
    }
    MEMSET(traj, 0, sizeof(md_trajectory_i));
    md_arena_allocator_destroy(vlx_traj->alloc);
}

md_trajectory_i* md_vlx_create_trajectory(const md_vlx_t* vlx, md_allocator_i* ext_alloc, uint32_t flags) {
	ASSERT(ext_alloc);
	if (!vlx) {
		MD_LOG_ERROR("vlx pointer is NULL");
		return NULL;
	}
    if (vlx->number_of_atoms == 0 || vlx->number_of_frames == 0) {
        MD_LOG_INFO("vlx object contains no trajectory data");
		return NULL;
	}

    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

	md_trajectory_i* traj = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(vlx_traj_t));
	if (!traj) {
		MD_LOG_ERROR("Failed to allocate memory for trajectory interface");
		return NULL;
	}

	MEMSET(traj, 0, sizeof(md_trajectory_i) + sizeof(vlx_traj_t));
    vlx_traj_t* vlx_traj = (vlx_traj_t*)(traj + 1);

	// Extract coordinate information from the vlx object and populate vlx_traj
    vlx_traj->magic = VLX_MAGIC;
    vlx_traj->num_frames = vlx->number_of_frames;
    vlx_traj->num_atoms = vlx->number_of_atoms;
    vlx_traj->alloc = alloc;

    size_t num_points = vlx_traj->num_frames * vlx_traj->num_atoms;

	md_array_resize(vlx_traj->atom_x, num_points, alloc);
    md_array_resize(vlx_traj->atom_y, num_points, alloc);
    md_array_resize(vlx_traj->atom_z, num_points, alloc);
    md_array_resize(vlx_traj->frame_times, vlx_traj->num_frames, alloc);

	ASSERT(3 * num_points == md_ndarray_f64_size(&vlx->atom_coordinates));
    for (size_t i = 0; i < num_points; ++i) {
		vlx_traj->atom_x[i] = (float)vlx->atom_coordinates.data[i * 3 + 0];
		vlx_traj->atom_y[i] = (float)vlx->atom_coordinates.data[i * 3 + 1];
		vlx_traj->atom_z[i] = (float)vlx->atom_coordinates.data[i * 3 + 2];
	}

	for (size_t i = 0; i < vlx_traj->num_frames; ++i) {
		vlx_traj->frame_times[i] = (double)i; // Placeholder: using frame index as time, replace with actual time if available
    }
	
    traj->inst = (struct md_trajectory_o*)vlx_traj;
    traj->free = vlx_traj_free;
	traj->get_header  = vlx_traj_get_header;
    traj->init_reader = vlx_traj_reader_init;

	return traj;
}

bool md_vlx_system_init_from_data(md_system_t* sys, const md_vlx_t* vlx) {
	ASSERT(sys);
	ASSERT(vlx);
	
	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("The veloxchem object contains no atoms");
		return false;
	}

	if (!sys->alloc) {
		MD_LOG_ERROR("System allocator is not set");
		return false;
    }

	md_system_reset(sys);

	size_t capacity = ROUND_UP(vlx->number_of_atoms, 16);

	md_array_resize(sys->atom.x,		capacity, sys->alloc);
	md_array_resize(sys->atom.y,		capacity, sys->alloc);
	md_array_resize(sys->atom.z,		capacity, sys->alloc);
    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags,    capacity, sys->alloc);

	MEMSET(sys->atom.x,			0, md_array_bytes(sys->atom.x));
	MEMSET(sys->atom.y,			0, md_array_bytes(sys->atom.y));
	MEMSET(sys->atom.z,			0, md_array_bytes(sys->atom.z));
	MEMSET(sys->atom.type_idx,  0, md_array_bytes(sys->atom.type_idx));
    MEMSET(sys->atom.flags,		0, md_array_bytes(sys->atom.flags));

    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    const dvec3_t* coords = (const dvec3_t*)vlx->atom_coordinates.data;
	for (size_t i = 0; i < vlx->number_of_atoms; ++i) {	
		sys->atom.x[i] = (float)coords[i].x;
		sys->atom.y[i] = (float)coords[i].y;
		sys->atom.z[i] = (float)coords[i].z;
		
		md_atomic_number_t z = vlx->atomic_numbers[i];
		str_t sym  = md_atomic_number_symbol(z);
        float mass = md_atomic_number_mass(z);
		float radius = md_atomic_number_vdw_radius(z);
		uint32_t color = md_atomic_number_cpk_color(z);

		md_atom_type_idx_t type_idx = md_atom_type_find_or_add(&sys->atom.type, sym, z, mass, radius, color, 0, sys->alloc);
		sys->atom.type_idx[i] = type_idx;
	}

	sys->atom.count = vlx->number_of_atoms;

    if (vlx->number_of_frames > 1) {
		sys->trajectory = md_vlx_create_trajectory(vlx, sys->alloc, 0);
	}

	return true;
}

bool md_vlx_system_init_from_file(md_system_t* sys, str_t filename) {
	ASSERT(sys);

    md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
	md_vlx_t* vlx = md_vlx_create(temp_arena);

	bool success = vlx_parse_file(vlx, filename, VLX_FLAG_CORE) && md_vlx_system_init_from_data(sys, vlx);

	md_arena_allocator_destroy(temp_arena);
	return success;
}

// Externally visible procedures

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

double md_vlx_nuclear_repulsion_energy(const md_vlx_t* vlx) {
	if (vlx) return vlx->nuclear_repulsion_energy;
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
	if (vlx) return (const dvec3_t*)vlx->atom_coordinates.data;
	return NULL;
}

const uint8_t* md_vlx_atomic_numbers(const md_vlx_t* vlx) {
	if (vlx) return vlx->atomic_numbers;
	return NULL;
}

const int* md_vlx_ao_to_atom_idx(const md_vlx_t* vlx) {
	if (vlx) return vlx->ao_to_atom_idx;
	return NULL;
}

md_vlx_scf_type_t md_vlx_scf_type(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.type;
	return MD_VLX_SCF_TYPE_UNKNOWN;
}

dvec3_t md_vlx_scf_ground_state_dipole_moment(const md_vlx_t* vlx) {
	if (vlx) {
        if (vlx->scf.ground_state_dipole_moment.data) {
			return *((dvec3_t*)(vlx->scf.ground_state_dipole_moment.data));
		}
	}
	return (dvec3_t){0};
}

size_t md_vlx_scf_homo_idx(const md_vlx_t* vlx, md_vlx_mo_type_t type) {
	if (vlx) {
		if (type == MD_VLX_MO_TYPE_ALPHA) {
			if (vlx->scf.alpha.homo_idx) {
				return vlx->scf.alpha.homo_idx[0];
			}
		} else if (type == MD_VLX_MO_TYPE_BETA) {
			if (vlx->scf.beta.homo_idx) {
				return vlx->scf.beta.homo_idx[0];
			}
		}
	}
	return 0;
}

size_t md_vlx_scf_lumo_idx(const md_vlx_t* vlx, md_vlx_mo_type_t type) {
	if (vlx) {
		if (type == MD_VLX_MO_TYPE_ALPHA) {
			if (vlx->scf.alpha.lumo_idx) {
				return vlx->scf.alpha.lumo_idx[0];
			}
		} else if (type == MD_VLX_MO_TYPE_BETA) {
			if (vlx->scf.beta.lumo_idx) {
				return vlx->scf.beta.lumo_idx[0];
			}
		}
	}
	return 0;
}

size_t md_vlx_scf_number_of_atomic_orbitals(const md_vlx_t* vlx) {
	if (vlx) {
		return number_of_atomic_orbitals(&vlx->scf.alpha);
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

bool md_vlx_scf_extract_gto_data(md_gto_data_t* out_gto_data, const md_vlx_t* vlx, double cutoff_value, md_allocator_i* alloc) {
	if (vlx) {
        ASSERT(out_gto_data);
		for (size_t i = 0; i < vlx->gto_data.num_cgtos; ++i) {
            double cgto_radius = 0.0;
			uint32_t cgto_offset = (uint32_t)out_gto_data->num_pgtos;
			for (size_t j = vlx->gto_data.cgto_offset[i]; j < vlx->gto_data.cgto_offset[i + 1]; ++j) {
				int pi, pj, pk, pl;
				md_gto_unpack_ijkl(vlx->gto_data.pgto_ijkl[j], &pi, &pj, &pk, &pl);
				double radius = md_gto_compute_radius_of_influence(pi, pj, pk, vlx->gto_data.pgto_coeff[j], vlx->gto_data.pgto_alpha[j], cutoff_value);
				if (radius == 0.0) {
					continue;
				}

                float pgto_coeff = (float)vlx->gto_data.pgto_coeff[j];
				float pgto_alpha = (float)vlx->gto_data.pgto_alpha[j];
				float pgto_radius = (float)radius;
				uint32_t pgto_ijkl = vlx->gto_data.pgto_ijkl[j];

				md_array_push(out_gto_data->pgto_coeff, pgto_coeff, alloc);
				md_array_push(out_gto_data->pgto_alpha, pgto_alpha, alloc);
				md_array_push(out_gto_data->pgto_radius, pgto_radius, alloc);
				md_array_push(out_gto_data->pgto_ijkl, pgto_ijkl, alloc);
				out_gto_data->num_pgtos += 1;

				cgto_radius = MAX(cgto_radius, radius);
            }
            // Push cgtos regardless of radius (otherwise indexing will be off)
            vec4_t cgto_xyzr = {vlx->gto_data.cgto_xyzr[i].x, vlx->gto_data.cgto_xyzr[i].y, vlx->gto_data.cgto_xyzr[i].z, (float)cgto_radius};
            md_array_push(out_gto_data->cgto_xyzr,   cgto_xyzr,   alloc);
            md_array_push(out_gto_data->cgto_offset, cgto_offset, alloc);
			
			out_gto_data->num_cgtos += 1;
        }
		// push final offset
		md_array_push(out_gto_data->cgto_offset, (uint32_t)out_gto_data->num_pgtos, alloc);

		return true;
	}
	return false;
}

static inline size_t get_matrix_index(size_t i, size_t j, size_t N) {
	size_t row = (i < j) ? i : j;
	size_t col = (i < j) ? j : i;
	size_t row_offset = row * (2 * N - row + 1) / 2;
	return row_offset + (col - row);
}

// The overlap matrix is a square, symmetric matrix [N][N], this returns the length N
size_t  md_vlx_scf_overlap_matrix_size(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.S.size[0];
	}
	return 0;
}

const double* md_vlx_scf_overlap_matrix_data(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.S.data;
	}
	return NULL;
}

// Returns the size of the density matrix in number of elements.
// This only contains the upper triangular part of the matrix, due to symmetry.
size_t md_vlx_scf_upper_triangular_density_matrix_size(const md_vlx_t* vlx) {
	if (vlx) {
        size_t N = vlx->scf.alpha.density.size[vlx->scf.alpha.density.rank - 1];
		return (N > 0) ? (N * (N + 1)) / 2 : 0;
	}
	return 0;
}

// Populates the provided array with the density matrix data.
bool md_vlx_scf_extract_upper_triangular_density_matrix_data(float* out_values, const md_vlx_t* vlx, md_vlx_mo_type_t type) {
	if (vlx) {
		size_t dim = 0;
        const double* density_data = NULL;
        if (type == MD_VLX_MO_TYPE_ALPHA) {
			dim = vlx->scf.alpha.density.size[vlx->scf.alpha.density.rank - 1];
			density_data = vlx->scf.alpha.density.data;
		} else if (type == MD_VLX_MO_TYPE_BETA) {
            dim = vlx->scf.beta.density.size[vlx->scf.beta.density.rank - 1];
			density_data = vlx->scf.beta.density.data;
		} else {
			MD_LOG_ERROR("Invalid MO type for density matrix extraction!");
			return false;
        }

		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = i; j < dim; ++j) {
				size_t idx = get_matrix_index(i, j, dim);
                out_values[idx] = (float)density_data[i * dim + j];  // Convert to float
            }
		}
		return true;
	}
	return false;
}

// Get the regular density matrix size N in (N x N)
size_t md_vlx_scf_density_matrix_size(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.alpha.density.size[vlx->scf.alpha.density.rank - 1];
	}
	return 0;
}

// Extracts the full density matrix into a square matrix representation
bool md_vlx_scf_extract_density_matrix_data(float* out_values, const struct md_vlx_t* vlx, md_vlx_mo_type_t type) {
	if (vlx) {
		size_t dim = vlx->scf.alpha.density.size[vlx->scf.alpha.density.rank - 1];
		const double* density_data = NULL;
		if (type == MD_VLX_MO_TYPE_ALPHA) {
			density_data = vlx->scf.alpha.density.data;
		} else if (type == MD_VLX_MO_TYPE_BETA) {
			density_data = vlx->scf.beta.density.data;
		} else {
			MD_LOG_ERROR("Invalid MO type for density matrix extraction!");
			return false;
		}
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				out_values[i * dim + j] = (float)density_data[i * dim + j];  // Convert to float
			}
		}
		return true;
	}
	return false;
}

// SCF History
size_t md_vlx_scf_history_size(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.number_of_iterations.data[vlx->scf.history.number_of_iterations.rank - 1];
	return 0;
}

const double* md_vlx_scf_history_energy(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.energy.data;
	return NULL;
}

const double* md_vlx_scf_history_energy_diff(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.energy_diff.data;
	return NULL;
}

const double* md_vlx_scf_history_density_diff(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.density_diff.data;
	return NULL;
}

const double* md_vlx_scf_history_gradient_norm(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.gradient_norm.data;
	return NULL;
}

const double* md_vlx_scf_history_max_gradient(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.history.max_gradient.data;
	return NULL;
}

bool md_vlx_parse_file(md_vlx_t* vlx, str_t filename) {
	return vlx_parse_file(vlx, filename, VLX_FLAG_ALL);
}
