#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <core/md_hash.h>

#include <md_system.h>
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
    md_vlx_2d_data_t density;
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
	size_t homo_idx[2];
	size_t lumo_idx[2];

	double energy;
	dvec3_t ground_state_dipole_moment;

	md_vlx_orbital_t alpha;
	md_vlx_orbital_t beta;

	double* resp_charges;

	md_vlx_2d_data_t S;
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

	double molecular_charge;
	double nuclear_repulsion_energy;
	size_t spin_multiplicity;

	// Arrays (length = number_of_atoms)
	dvec3_t* atom_coordinates;

	md_vlx_atomic_property_t* atomic_properties; // Optional data, may be NULL, length is number of atomic properties

	// Data blocks
	md_vlx_scf_t scf;
	md_vlx_rsp_t rsp;
	md_vlx_vib_t vib;
	md_vlx_opt_t opt;


	md_element_t* atomic_numbers;
	int* ao_to_atom_idx;    // Maps atomic orbitals to atom indices (shell order)
	int* local_to_global_atom_idx; // Maps local atom indices to global system indices for subsystems. NULL if not a subsystem.
	// ao_remap[shell_ao_idx] = vlx_ao_idx
	// Maps from shell order (angl→atom→func→isph) to VeloxChem matrix row order (angl→isph→atom→func).
	// Built once after the basis set is parsed; used to permute C, D, S matrices into shell order.
	int* ao_remap;

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

// Build a permutation table that maps from shell order (angl→atom→func→isph)
// to VeloxChem matrix row order (angl→isph→atom→func).
// ao_remap[shell_ao_idx] = vlx_ao_idx.
// Returns the total number of AOs (length of the table), or 0 on failure.
static size_t build_ao_remap(int** out_remap, const md_vlx_t* vlx, md_allocator_i* alloc) {
	ASSERT(out_remap);
	ASSERT(vlx);
	ASSERT(alloc);

	int natoms   = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	// First, count total AOs
	size_t num_aos = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			basis_func_t basis_funcs[128];
			size_t num_funcs = basis_set_extract_atomic_basis_func_angl(
				basis_funcs, ARRAY_SIZE(basis_funcs), &vlx->basis_set, idelem, angl);
			num_aos += (size_t)nsph * num_funcs;
		}
	}

	if (num_aos == 0) return 0;

	int* remap = (int*)md_alloc(alloc, sizeof(int) * num_aos);
	MEMSET(remap, 0, sizeof(int) * num_aos);

	// VeloxChem matrix AO index: angl → isph → atom → func
	// Shell AO index:            angl → atom → func → isph
	//
	// We walk shell order (outer loop) and record where each entry maps in VLX order.
	// vlx_ao_start[angl][isph][atom] is needed. We pre-compute the vlx base offset per
	// (angl, isph, atom) by walking the vlx ordering once.

	// Compute vlx_base[angl][isph][atomidx] = starting VLX AO index for that group.
	// Flat: vlx_base[angl * (max_angl+1) * natoms + isph * natoms + atomidx]
	// But isph is up to 2*max_angl+1. Use a temp VLA-style allocation.
	int max_nsph = 2 * max_angl + 1;
	size_t temp_pos = md_temp_get_pos();
	int* vlx_base = (int*)md_temp_push(sizeof(int) * (max_angl + 1) * max_nsph * natoms);
	int* vlx_num  = (int*)md_temp_push(sizeof(int) * (max_angl + 1) * natoms);
	MEMSET(vlx_base, 0, sizeof(int) * (max_angl + 1) * max_nsph * natoms);
	MEMSET(vlx_num,  0, sizeof(int) * (max_angl + 1) * natoms);

	// Count funcs per (angl, atom) so we know stride for vlx_base
	for (int angl = 0; angl <= max_angl; angl++) {
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			basis_func_t bf[128];
			size_t nf = basis_set_extract_atomic_basis_func_angl(bf, ARRAY_SIZE(bf), &vlx->basis_set, idelem, angl);
			vlx_num[angl * natoms + atomidx] = (int)nf;
		}
	}

	// Walk VLX order to fill vlx_base
	int vlx_idx = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int isph = 0; isph < nsph; isph++) {
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				vlx_base[angl * max_nsph * natoms + isph * natoms + atomidx] = vlx_idx;
				vlx_idx += vlx_num[angl * natoms + atomidx];
			}
		}
	}

	// Now walk shell order and fill remap
	int shell_idx = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int nfuncs = vlx_num[angl * natoms + atomidx];
			for (int funcidx = 0; funcidx < nfuncs; funcidx++) {
				for (int isph = 0; isph < nsph; isph++, shell_idx++) {
					int base = vlx_base[angl * max_nsph * natoms + isph * natoms + atomidx];
					remap[shell_idx] = base + funcidx;
				}
			}
		}
	}

	md_temp_set_pos_back(temp_pos);
	*out_remap = remap;
	return num_aos;
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
		vlx->molecular_charge			= (double)field_vals[0];
		vlx->spin_multiplicity			= (size_t)field_vals[1];
		vlx->number_of_atoms			= (size_t)field_vals[2];
		vlx->number_of_alpha_electrons	= (size_t)field_vals[3];
		vlx->number_of_beta_electrons	= (size_t)field_vals[4];
	}

	md_array_grow(vlx->atomic_numbers, count, alloc);
	md_array_grow(vlx->atom_coordinates, count, alloc);

	for (size_t i = 0; i < count; ++i) {
		md_atomic_number_t anum = md_atomic_number_from_symbol(LBL_TO_STR(fields[i].sym), false);
		if (anum == 0) {
			MD_LOG_ERROR("Unrecognized element '%s' in geometry", fields[i].sym);
			goto done;
		}
		vlx->atomic_numbers[i] = anum;
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

	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	hid_t datatype_id = H5Dget_type(dataset_id);
	hid_t space_id = H5Dget_space(dataset_id);

	if (datatype_id == H5I_INVALID_HID || space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 datatype/space for dataset: '%s'", field_name);
		goto done;
	}

	if (H5Tget_class(datatype_id) != H5T_STRING) {
		MD_LOG_ERROR("H5 dataset is not a string: '%s'", field_name);
		goto done;
	}

	if (H5Tis_variable_str(datatype_id)) {
		char* tmp = NULL;
		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
		if (status != 0) {
			MD_LOG_ERROR("Failed to read variable-length string for H5 dataset: '%s'", field_name);
			goto done;
		}

		size_t len = tmp ? strlen(tmp) : 0;
		str_t data = str_alloc(len, alloc);
		if (len > 0) {
			MEMCPY((char*)data.ptr, tmp, len);
		}
		*str = data;

		if (tmp) {
			H5free_memory(tmp);
		}
	} else {
		const size_t raw_len = H5Tget_size(datatype_id);
		char* tmp = md_alloc(md_get_heap_allocator(), raw_len + 1);
		if (!tmp) {
			MD_LOG_ERROR("Failed to allocate temporary buffer for H5 dataset: '%s'", field_name);
			goto done;
		}
		MEMSET(tmp, 0, raw_len + 1);

		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		if (status != 0) {
			MD_LOG_ERROR("Failed to read fixed-length string for H5 dataset: '%s'", field_name);
			md_free(md_get_heap_allocator(), tmp, raw_len + 1);
			goto done;
		}

		size_t len = strnlen(tmp, raw_len);
		str_t data = str_alloc(len, alloc);
		if (len > 0) {
			MEMCPY((char*)data.ptr, tmp, len);
		}
		*str = data;

		md_free(md_get_heap_allocator(), tmp, raw_len + 1);
	}

	result = true;

done:
	if (datatype_id != H5I_INVALID_HID) H5Tclose(datatype_id);
	if (space_id != H5I_INVALID_HID) H5Sclose(space_id);
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

	hsize_t num_points = H5Sget_simple_extent_npoints(space_id);

	if (num_points != num_samples) {
		MD_LOG_ERROR("Unexpected number of points when reading dataset, got %i, expected %i", num_points, num_samples);
		goto done;
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

// Scan through an h5 group (no recursion) and check for datasets with the attribute 'atomic_property', if found, this will contain the label of the property.
// If read, we will attempt to extract this as an array of doubles (if it has the length of number of atoms).
static bool h5_read_atomic_properties(md_vlx_t* vlx, hid_t group_handle) {
	H5G_info_t info = { 0 };
	if (H5Gget_info(group_handle, &info) < 0) {
		MD_LOG_ERROR("Failed to get group info when reading atomic properties");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		ssize_t size = H5Gget_objname_by_idx(group_handle, i, name_buf, sizeof(name_buf));
		if (size < 0) {
			continue;
		}
		H5G_obj_t type = H5Gget_objtype_by_idx(group_handle, i);
		MD_LOG_DEBUG("obj_name: '%s', obj_type: %i", name_buf, type);

		// Ensure that the type is a dataset, if not we skip
		if (type != H5G_DATASET) {
			continue;
		}
		hid_t dataset_id = H5Dopen(group_handle, name_buf, H5P_DEFAULT);
		if (dataset_id == H5I_INVALID_HID) {
			continue;
		}

		if (!H5Aexists(dataset_id, "atomic_property")) {
			continue;
		}

		hid_t attr_id = H5Aopen(dataset_id, "atomic_property", H5P_DEFAULT);
		if (attr_id == H5I_INVALID_HID) {
			H5Dclose(dataset_id);
			continue;
		}
		char property_label[256] = { 0 };
		hid_t attr_type = H5Aget_type(attr_id);

		// Check if attribute type is variable length string or fixed length string, and read accordingly
		if (H5Tis_variable_str(attr_type)) {
			char* var_str;
			H5Aread(attr_id, attr_type, &var_str);
			strncpy(property_label, var_str, sizeof(property_label));
			H5free_memory(var_str);
		} else if (H5Tget_class(attr_type) == H5T_STRING) {
			H5Aread(attr_id, attr_type, property_label);
		} else {
			MD_LOG_ERROR("Unexpected attribute type for 'atomic_property' attribute in dataset '%s'", name_buf);
			goto done;
		}

		// If we end up here we did not fail, but we may not have read a property label either, so we check if we got something valid
		if (property_label[0] == '\0') {
			// Not an error, just use the field name as the property label
			strncpy(property_label, name_buf, sizeof(property_label));
		}
		property_label[sizeof(property_label) - 1] = '\0';

		// We have a property label, we attempt to read the dataset as an array of doubles with the length of number of atoms
		hid_t space_id = H5Dget_space(dataset_id);
		if (space_id == H5I_INVALID_HID) {
			MD_LOG_ERROR("Failed to get dataspace for dataset '%s'", name_buf);
			goto done;
		}

		int num_dims = H5Sget_simple_extent_ndims(space_id);
		if (num_dims < 0) {
			MD_LOG_ERROR("Failed to get number of dimensions for dataset '%s'", name_buf);
			H5Sclose(space_id);
			goto done;
		}
		if (num_dims > 2) {
			MD_LOG_ERROR("Too many dimensions for atomic property dataset '%s', expected at most 3, got %i", name_buf, num_dims);
			H5Sclose(space_id);
			goto done;
		}

		size_t dims[2] = { 0 };
		H5Sget_simple_extent_dims(space_id, (hsize_t*)dims, 0);

		// We expect the inner most dimension to be the number of atoms. Otherwise we skip.
		if (dims[num_dims - 1] != vlx->number_of_atoms) {
			MD_LOG_ERROR("Unexpected size of innermost dimension for atomic property dataset '%s', expected %zu, got %zu", name_buf, vlx->number_of_atoms, dims[num_dims - 1]);
			H5Sclose(space_id);
			continue;
		}

		size_t num_points = H5Sget_simple_extent_npoints(space_id);
		
		H5Sclose(space_id);

		// Construct a unique uint64_t key for this property.
		uint64_t key = md_hash64(name_buf, sizeof(name_buf), 0);
		
		md_vlx_atomic_property_t property = {
			 .label = str_copy_cstr(property_label, vlx->arena),
			 .key = key,
			 .dim[0] = vlx->number_of_atoms,
			 .dim[1] = num_dims > 1 ? dims[0] : 1,
			 .data = NULL,
		};

		md_array_resize(property.data, num_points, vlx->arena);
		herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, property.data);
		if (status != 0) {
			MD_LOG_ERROR("Failed to read data for atomic property dataset '%s'", name_buf);
			goto done;
		}

		md_array_push(vlx->atomic_properties, property, vlx->arena);
	done:
		H5Tclose(attr_type);
		H5Aclose(attr_id);
		H5Dclose(dataset_id);
	}
	return true;
}

// ---------------------------------------------------------------------------
// AO permutation helpers
// Reorder AO-indexed matrices from VeloxChem order into shell order in-place.
// ---------------------------------------------------------------------------

// Permute rows of an AO×MO matrix (num_ao rows, num_mo columns, row-major).
// remap[shell_ao] = vlx_ao  =>  dst_row[shell_ao] = src_row[vlx_ao]
// Permute rows of a [num_ao x num_mo] matrix according to remap and transpose to [num_mo x num_ao].
// On entry  mat is [num_ao][num_mo] in VeloxChem AO order.
// On return mat is [num_mo][num_ao] in shell order — each MO is a contiguous row.
static void ao_permute(double* mat, size_t num_ao, size_t num_mo, const int* remap) {
	size_t temp_pos = md_temp_get_pos();
	double* tmp = (double*)md_temp_push(sizeof(double) * num_ao * num_mo);
	MEMCPY(tmp, mat, sizeof(double) * num_ao * num_mo);
	for (size_t mo = 0; mo < num_mo; mo++) {
		for (size_t ao = 0; ao < num_ao; ao++) {
			mat[mo * num_ao + ao] = tmp[(size_t)remap[ao] * num_mo + mo];
		}
	}
	md_temp_set_pos_back(temp_pos);
}

// Permute both rows and columns of a square AO×AO matrix (num_ao × num_ao, row-major).
static void ao_permute_square(double* mat, size_t num_ao, const int* remap) {
	size_t temp_pos = md_temp_get_pos();
	double* tmp = (double*)md_temp_push(sizeof(double) * num_ao * num_ao);
	MEMCPY(tmp, mat, sizeof(double) * num_ao * num_ao);
	for (size_t i = 0; i < num_ao; i++) {
		size_t si = (size_t)remap[i];
		for (size_t j = 0; j < num_ao; j++) {
			size_t sj = (size_t)remap[j];
			mat[i * num_ao + j] = tmp[si * num_ao + sj];
		}
	}
	md_temp_set_pos_back(temp_pos);
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

	// Density dimensions (May differ from dim is always square)
	size_t den_dim[2];
    h5_read_dataset_dims(den_dim, 2, handle, "D_alpha");

	md_array_resize(vlx->scf.alpha.coefficients.data, dim[0] * dim[1], vlx->arena);
	MEMCPY(vlx->scf.alpha.coefficients.size, dim, sizeof(dim));

	md_array_resize(vlx->scf.alpha.energy.data, dim[1], vlx->arena);
	vlx->scf.alpha.energy.size = dim[1];

	md_array_resize(vlx->scf.alpha.occupancy.data, dim[1], vlx->arena);
	vlx->scf.alpha.occupancy.size = dim[1];

	md_array_resize(vlx->scf.alpha.density.data, den_dim[0] * den_dim[1], vlx->arena);
    MEMCPY(vlx->scf.alpha.density.size, den_dim, sizeof(den_dim));

	// Extract alpha data
	if (!h5_read_dataset_data(vlx->scf.alpha.coefficients.data, md_array_size(vlx->scf.alpha.coefficients.data), handle, H5T_NATIVE_DOUBLE, "C_alpha")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.alpha.energy.data, md_array_size(vlx->scf.alpha.energy.data), handle, H5T_NATIVE_DOUBLE, "E_alpha")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.alpha.occupancy.data, md_array_size(vlx->scf.alpha.occupancy.data), handle, H5T_NATIVE_DOUBLE, "occ_alpha")) {
		return false;
	}
    if (!h5_read_dataset_data(vlx->scf.alpha.density.data, md_array_size(vlx->scf.alpha.density.data), handle, H5T_NATIVE_DOUBLE, "D_alpha")) {
        return false;
    }

	if (vlx->scf.type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
		md_array_resize(vlx->scf.beta.coefficients.data, dim[0] * dim[1], vlx->arena);
		MEMCPY(vlx->scf.beta.coefficients.size, dim, sizeof(dim));

		md_array_resize(vlx->scf.beta.energy.data, dim[1], vlx->arena);
		vlx->scf.beta.energy.size = dim[1];

		md_array_resize(vlx->scf.beta.occupancy.data, dim[1], vlx->arena);
		vlx->scf.beta.occupancy.size = dim[1];

		md_array_resize(vlx->scf.beta.density.data, den_dim[0] * den_dim[0], vlx->arena);
        MEMCPY(vlx->scf.beta.density.size, den_dim, sizeof(den_dim));

		// Extract beta data
		if (!h5_read_dataset_data(vlx->scf.beta.coefficients.data, md_array_size(vlx->scf.beta.coefficients.data), handle, H5T_NATIVE_DOUBLE, "C_beta")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.beta.energy.data, md_array_size(vlx->scf.beta.energy.data), handle, H5T_NATIVE_DOUBLE, "E_beta")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, md_array_size(vlx->scf.beta.occupancy.data), handle, H5T_NATIVE_DOUBLE, "occ_beta")) {
			return false;
		}
        if (!h5_read_dataset_data(vlx->scf.beta.density.data, md_array_size(vlx->scf.beta.density.data), handle, H5T_NATIVE_DOUBLE, "D_beta")) {
            return false;
        }
	} else {
		// Shallow copy fields from Alpha
		MEMCPY(&vlx->scf.beta, &vlx->scf.alpha, sizeof(md_vlx_orbital_t));
		if (vlx->scf.type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
			vlx->scf.beta.occupancy.data = 0;
			md_array_resize(vlx->scf.beta.occupancy.data, vlx->scf.beta.occupancy.size, vlx->arena);
			if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, md_array_size(vlx->scf.beta.occupancy.data), handle, H5T_NATIVE_DOUBLE, "occ_beta")) {
				return false;
			}
		}
	}

	// S matrix is overlap (notice dimension is the same as D)
	md_array_resize(vlx->scf.S.data, den_dim[0] * den_dim[1], vlx->arena);
    MEMCPY(vlx->scf.S.size, den_dim, sizeof(den_dim));

	if (!h5_read_dataset_data(vlx->scf.S.data, md_array_size(vlx->scf.S.data), handle, H5T_NATIVE_DOUBLE, "S")) {
		return false;
	}

	// The ground state dipole moment is not present in all versions
	if (!h5_read_dataset_data(&vlx->scf.ground_state_dipole_moment, 3, handle, H5T_NATIVE_DOUBLE, "dipole_moment")) {
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

	if (!h5_read_dataset_data(vlx->scf.history.density_diff, scf_hist_iter, handle, H5T_NATIVE_DOUBLE, "scf_history_diff_density")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.energy_diff, scf_hist_iter, handle, H5T_NATIVE_DOUBLE, "scf_history_diff_energy")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.energy, scf_hist_iter, handle, H5T_NATIVE_DOUBLE, "scf_history_energy")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.gradient_norm, scf_hist_iter, handle, H5T_NATIVE_DOUBLE, "scf_history_gradient_norm")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.history.max_gradient, scf_hist_iter, handle, H5T_NATIVE_DOUBLE, "scf_history_max_gradient")) {
		return false;
	}

	{
		size_t charge_resp_dim;
		if (h5_read_dataset_dims(&charge_resp_dim, 1, handle, "charges_resp")) {
			md_array_resize(vlx->scf.resp_charges, charge_resp_dim, vlx->arena);
			if (!h5_read_dataset_data(vlx->scf.resp_charges, md_array_size(vlx->scf.resp_charges), handle, H5T_NATIVE_DOUBLE, "charges_resp")) {
				MD_LOG_ERROR("Could not read charges_resp");
				return false;
			}
		}
	}

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

		snprintf(buf, sizeof(buf), "NTO_S%i_alpha_orbitals", idx);
		if (!h5_read_dataset_dims(dim, 2, handle, buf)) {
			return false;
		}
		if (dim[0] == 0 || dim[1] == 0) {
			MD_LOG_ERROR("Invalid dimensions in NTO orbitals");
			return false;
		}

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
	}

	return true;
}

static bool h5_read_rsp_data(md_vlx_t* vlx, hid_t handle) {
	h5_read_scalar(&vlx->rsp.number_of_excited_states, handle, H5T_NATIVE_HSIZE, "number_of_states");
	if (vlx->rsp.number_of_excited_states > 0) {
		// Standard Linear Response data, allocate and read

		// Allocate data
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
		if (h5_check_dataset_exists(handle, "NTO_S1_alpha_orbitals")) {
			if (!h5_read_nto_data(&vlx->rsp, handle, vlx->arena)) {
				return false;
			}
		} else {
			// Check for 'nto' folder inside
			if (H5Lexists(handle, "nto", H5P_DEFAULT)) {
				hid_t nto_group = H5Gopen(handle, "nto", H5P_DEFAULT);
				if (nto_group >= 0) {
					bool result = h5_read_nto_data(&vlx->rsp, nto_group, vlx->arena);
					H5Gclose(nto_group);
					if (!result) {
						return false;
					}
				}
			}
		}

		// Dipoles
		size_t num_dipole_points = vlx->rsp.number_of_excited_states * 3;
		if (!h5_read_dataset_data(vlx->rsp.electric_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "electric_transition_dipoles")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->rsp.magnetic_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "magnetic_transition_dipoles")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->rsp.velocity_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "velocity_transition_dipoles")) {
			return false;
		}

		// Abs, rot and osc
		if (!h5_read_dataset_data(vlx->rsp.absorption_ev, md_array_size(vlx->rsp.absorption_ev), handle, H5T_NATIVE_DOUBLE, "eigenvalues")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->rsp.oscillator_strengths, md_array_size(vlx->rsp.oscillator_strengths), handle, H5T_NATIVE_DOUBLE, "oscillator_strengths")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->rsp.rotatory_strengths, md_array_size(vlx->rsp.rotatory_strengths), handle, H5T_NATIVE_DOUBLE, "rotatory_strengths")) {
			return false;
		}

		// Convert Atomic units (Hartree) to eV
		if (vlx->rsp.absorption_ev) {
			for (size_t i = 0; i < vlx->rsp.number_of_excited_states; ++i) {
				vlx->rsp.absorption_ev[i] *= HARTREE_TO_EV;
			}
		}
	}

	// CPP data is optional, only read if present
	size_t dim;
	if (h5_read_dataset_dims(&dim, 1, handle, "frequencies")) {
		vlx->rsp.cpp.number_of_frequencies = dim;
		md_array_resize(vlx->rsp.cpp.frequencies, dim, vlx->arena);
		if (!h5_read_dataset_data(vlx->rsp.cpp.frequencies, md_array_size(vlx->rsp.cpp.frequencies), handle, H5T_NATIVE_DOUBLE, "frequencies")) {
			// Frequencies has to be present if cpp section is present, fail if missing
			return false;
		}
		
		if (h5_check_dataset_exists(handle, "sigma")) {
			md_array_resize(vlx->rsp.cpp.sigmas, dim, vlx->arena);
			if (!h5_read_dataset_data(vlx->rsp.cpp.sigmas, md_array_size(vlx->rsp.cpp.sigmas), handle, H5T_NATIVE_DOUBLE, "sigma")) {
				return false;
			}
		}

		if (h5_check_dataset_exists(handle, "delta-epsilon")) {
			md_array_resize(vlx->rsp.cpp.delta_epsilon, dim, vlx->arena);
			if (!h5_read_dataset_data(vlx->rsp.cpp.delta_epsilon, md_array_size(vlx->rsp.cpp.delta_epsilon), handle, H5T_NATIVE_DOUBLE, "delta-epsilon")) {
				return false;
			}
		}
	}

	return true;
}

static bool h5_read_vib_data(md_vlx_t* vlx, hid_t handle) {
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

	vlx->vib.number_of_normal_modes = number_of_modes;

	md_array_resize(vlx->vib.force_constants, number_of_modes, vlx->arena);
	if (!h5_read_dataset_data(vlx->vib.force_constants, md_array_size(vlx->vib.force_constants), handle, H5T_NATIVE_DOUBLE, "force_constants")) {
		return false;
	}

	md_array_resize(vlx->vib.ir_intensities, number_of_modes, vlx->arena);
	if (!h5_read_dataset_data(vlx->vib.ir_intensities, md_array_size(vlx->vib.ir_intensities), handle, H5T_NATIVE_DOUBLE, "ir_intensities")) {
		return false;
	}

	md_array_resize(vlx->vib.frequencies, number_of_modes, vlx->arena);
	if (!h5_read_dataset_data(vlx->vib.frequencies, md_array_size(vlx->vib.frequencies), handle, H5T_NATIVE_DOUBLE, "vib_frequencies")) {
		return false;
	}

	md_array_resize(vlx->vib.reduced_masses, number_of_modes, vlx->arena);
	if (!h5_read_dataset_data(vlx->vib.reduced_masses, md_array_size(vlx->vib.reduced_masses), handle, H5T_NATIVE_DOUBLE, "reduced_masses")) {
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
			if (vlx->number_of_atoms == 0) {
				MD_LOG_ERROR("Missing number of atoms, is required for normal modes");
				return false;
			}

            char lbl[32];
            for (size_t i = 0; i < vlx->vib.number_of_normal_modes; ++i) {
                snprintf(lbl, sizeof(lbl), "%zu", i + 1);

                dvec3_t* data = md_array_create(dvec3_t, vlx->number_of_atoms, vlx->arena);
                MEMSET(data, 0, sizeof(dvec3_t) * vlx->number_of_atoms);

                if (!h5_read_dataset_data(data, 3 * vlx->number_of_atoms, normal_modes_id, H5T_NATIVE_DOUBLE, lbl)) {
                    MD_LOG_ERROR("Failed to extract dataset in '%s' normal mode", lbl);
                    md_array_free(data, vlx->arena);
                    return false;
                }

                // Success, append ata
                md_array_push(vlx->vib.normal_modes, data, vlx->arena);
            }
        }
	} else if (obj_type == H5I_DATASET) {
		// Handle dataset case
		// Iterate over outer dimension in dataset, which should be [number_of_normal_modes][number_of_atoms][3]

		size_t data_dim[3];
		int num_dim = h5_read_dataset_dims(data_dim, 3, handle, "normal_modes");

		// Assert expected dimensions
		if (num_dim != 3 || data_dim[0] != vlx->vib.number_of_normal_modes || data_dim[1] != vlx->number_of_atoms || data_dim[2] != 3) {
			MD_LOG_ERROR("Unexpected dimensions in normal_modes dataset");
			H5Oclose(obj_info);
			return false;
		}

		size_t num_points = data_dim[0] * data_dim[1] * data_dim[2];
		double* raw_data = md_array_create(double, num_points, vlx->arena);
		if (!h5_read_dataset_data(raw_data, num_points, handle, H5T_NATIVE_DOUBLE, "normal_modes")) {
			MD_LOG_ERROR("Failed to read normal_modes dataset");
			md_array_free(raw_data, vlx->arena);
			H5Oclose(obj_info);
			return false;
		}

		// Set the pointers to each normal mode (within raw_data)
		dvec3_t* base_ptr = (dvec3_t*)raw_data;
		for (size_t i = 0; i < vlx->vib.number_of_normal_modes; ++i) {
			dvec3_t* mode_data = base_ptr + (i * vlx->number_of_atoms);
			md_array_push(vlx->vib.normal_modes, mode_data, vlx->arena);
		}
	} else {
		MD_LOG_ERROR("Unrecognized object type for 'normal_modes'");
		H5Oclose(obj_info);
		return false;
	}

	return true;
}

static bool h5_read_opt_data(md_vlx_t* vlx, hid_t handle) {
	// @TODO(This will likely be exposed as its own variable in the future, for now we extract the length from one of the fields)
	size_t dim[3];
	int num_dim;
	
	num_dim = h5_read_dataset_dims(dim, 3, handle, "nuclear_repulsion_energies");
	if (num_dim <= 0) {
		return false;
	}

	// This is a fix because the input data in one version is supplied as a 2D object
	size_t len = dim[0];
	vlx->opt.number_of_steps = len;

	md_array_resize(vlx->opt.nuclear_repulsion_energies, len, vlx->arena);
	if (!h5_read_dataset_data(vlx->opt.nuclear_repulsion_energies, md_array_size(vlx->opt.nuclear_repulsion_energies), handle, H5T_NATIVE_DOUBLE, "nuclear_repulsion_energies")) {
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

	md_array_resize(vlx->opt.energies, len, vlx->arena);
	if (!h5_read_dataset_data(vlx->opt.energies, md_array_size(vlx->opt.energies), handle, H5T_NATIVE_DOUBLE, energy_ident)) {
		return false;
	}

	num_dim = h5_read_dataset_dims(dim, 3, handle, coord_ident);
	if (dim[0] != len || dim[1] != vlx->number_of_atoms || dim[2] != 3) {
		MD_LOG_ERROR("Inconsistent or invalid opt_coordinates dimensions");
		return false;
	}

	md_array_resize(vlx->opt.coordinates, dim[0] * dim[1], vlx->arena);
	if (!h5_read_dataset_data(vlx->opt.coordinates, md_array_size(vlx->opt.coordinates) * 3, handle, H5T_NATIVE_DOUBLE, coord_ident)) {
		return false;
	}

	if (vlx->opt.coordinates) {
		for (size_t i = 0; i < dim[0] * dim[1]; ++i) {
			vlx->opt.coordinates[i] = dvec3_mul_f(vlx->opt.coordinates[i], BOHR_TO_ANGSTROM);
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
	if (!h5_read_dataset_data(vlx->atom_coordinates, md_array_size(vlx->atom_coordinates) * 3, handle, H5T_NATIVE_DOUBLE, "atom_coordinates")) {
		return false;
	}

	md_array_resize(vlx->atomic_numbers, vlx->number_of_atoms, vlx->arena);
	MEMSET(vlx->atomic_numbers, 0, md_array_bytes(vlx->atomic_numbers));
	if (!h5_read_dataset_data(vlx->atomic_numbers, md_array_size(vlx->atomic_numbers), handle, H5T_NATIVE_UINT8, "nuclear_charges")) {
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

	if (h5_check_dataset_exists(handle, "qm_atom_indices")) {
        md_array_resize(vlx->local_to_global_atom_idx, vlx->number_of_atoms, vlx->arena);
		if (!h5_read_dataset_data(vlx->local_to_global_atom_idx, vlx->number_of_atoms, handle, H5T_NATIVE_INT32, "qm_atom_indices")) {
			return false;
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
		h5_read_atomic_properties(vlx, file_id);
	}

	// SCF
	if (flags & VLX_FLAG_SCF) {
		if (H5Lexists(file_id, "scf", H5P_DEFAULT) > 0) {
			hid_t scf_id = H5Gopen(file_id, "scf", H5P_DEFAULT);
			if (scf_id != H5I_INVALID_HID) {
				result = h5_read_scf_data(vlx, scf_id);
				h5_read_atomic_properties(vlx, scf_id);
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
                result = h5_read_vib_data(vlx, vib_id);
				h5_read_atomic_properties(vlx, vib_id);
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
				result = h5_read_opt_data(vlx, opt_id);
				h5_read_atomic_properties(vlx, opt_id);
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
				result = h5_read_rsp_data(vlx, rsp_id);
				h5_read_atomic_properties(vlx, rsp_id);
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
	md_file_t file = {0};
	if (!md_file_open(&file, filename, MD_FILE_READ)) {
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
	md_file_close(&file);

	if (vlx->number_of_atoms == 0 || str_empty(vlx->basis_set_ident) || vlx->atomic_numbers == 0) {
		MD_LOG_ERROR("Failed to parse required information in out file");
		goto done;
	}

	str_t base_file = {0};
	extract_file_path_without_ext(&base_file, filename);

	md_strb_t sb = md_strb_create(md_get_temp_allocator());

	if (flags & VLX_FLAG_SCF) {
		// Attempt to read scf data
		md_strb_reset(&sb);
		md_strb_fmt(&sb, STR_FMT ".scf.results.h5", STR_ARG(base_file));
		if (md_path_is_valid(md_strb_to_str(sb))) {
			// Open an existing file
			hid_t file_id = H5Fopen(md_strb_to_cstr(sb), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (file_id == H5I_INVALID_HID) {
				MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(md_strb_to_str(sb)));
				goto done;
			}
			if (!h5_read_scf_data(vlx, file_id)) {
				goto done;
			}
		} else {
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

				if (!h5_read_dataset_data(vlx->scf.alpha.coefficients.data, md_array_size(vlx->scf.alpha.coefficients.data), file_id, H5T_NATIVE_DOUBLE, "alpha_orbitals")) {
					goto done;
				}
				if (!h5_read_dataset_data(vlx->scf.alpha.energy.data, md_array_size(vlx->scf.alpha.energy.data), file_id, H5T_NATIVE_DOUBLE, "alpha_energies")) {
					goto done;
				}
				if (!h5_read_dataset_data(vlx->scf.alpha.occupancy.data, md_array_size(vlx->scf.alpha.occupancy.data), file_id, H5T_NATIVE_DOUBLE, "alpha_occupations")) {
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
					if (!h5_read_dataset_data(vlx->scf.beta.coefficients.data, md_array_size(vlx->scf.beta.coefficients.data), file_id, H5T_NATIVE_DOUBLE, "beta_orbitals")) {
						goto done;
					}
					if (!h5_read_dataset_data(vlx->scf.beta.energy.data, md_array_size(vlx->scf.beta.energy.data), file_id, H5T_NATIVE_DOUBLE, "beta_energies")) {
						goto done;
					}
					if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, md_array_size(vlx->scf.beta.occupancy.data), file_id, H5T_NATIVE_DOUBLE, "beta_occupations")) {
						goto done;
					}
				} else {
					// Shallow copy fields from Alpha
					MEMCPY(&vlx->scf.beta, &vlx->scf.alpha, sizeof(md_vlx_orbital_t));

					if (vlx->scf.type == MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL) {
						vlx->scf.beta.occupancy.data = 0;
						md_array_resize(vlx->scf.beta.occupancy.data, vlx->scf.beta.occupancy.size, vlx->arena);
						if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, md_array_size(vlx->scf.beta.occupancy.data), file_id, H5T_NATIVE_DOUBLE, "beta_occupations")) {
							goto done;
						}
					}
				}
			}
		}
	}
	if (flags & VLX_FLAG_RSP) {
		// Test if we can read the first NTO file
		md_strb_reset(&sb);
		md_strb_fmt(&sb, STR_FMT "_S1_NTO.h5", STR_ARG(base_file));
		if (H5Fopen(md_strb_to_cstr(sb), H5F_ACC_RDONLY, H5P_DEFAULT) != H5I_INVALID_HID) {
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

					size_t dim[2];
					if (!h5_read_dataset_dims(dim, 2, file_id, "alpha_orbitals")) {
						goto done;
					}
					if (dim[0] == 0 || dim[1] == 0) {
						MD_LOG_ERROR("Invalid dimensions in NTO orbitals");
						goto done;
					}

					md_array_push(vlx->rsp.nto, (md_vlx_orbital_t) { 0 }, vlx->arena);
					md_vlx_orbital_t* nto = md_array_last(vlx->rsp.nto);

					md_array_resize(nto->coefficients.data, dim[0] * dim[1], vlx->arena);
					MEMCPY(nto->coefficients.size, dim, sizeof(dim));

					md_array_resize(nto->occupancy.data, dim[1], vlx->arena);
					nto->occupancy.size = dim[1];

					if (!h5_read_dataset_data(nto->coefficients.data, md_array_size(nto->coefficients.data), file_id, H5T_NATIVE_DOUBLE, "alpha_orbitals")) {
						goto done;
					}
					if (!h5_read_dataset_data(nto->occupancy.data, md_array_size(nto->occupancy.data), file_id, H5T_NATIVE_DOUBLE, "alpha_occupations")) {
						goto done;
					}

					// Permute+transpose into [num_mo][num_ao] shell order
					if (vlx->ao_remap) {
						size_t nao = nto->coefficients.size[0];
						size_t nmo = nto->coefficients.size[1];
						ao_permute(nto->coefficients.data, nao, nmo, vlx->ao_remap);
						nto->coefficients.size[0] = nmo;
						nto->coefficients.size[1] = nao;
					}
				}
				else {
					MD_LOG_INFO("The veloxchem object specified %zu excited states, but the matching NTOs could not be found.", vlx->rsp.number_of_excited_states);
					if (vlx->rsp.nto) {
						md_array_free(vlx->rsp.nto, vlx->arena);
						vlx->rsp.nto = NULL;
						break;
					}
				}
			}
		}
	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);
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

	// Build the AO remap table and apply it to all loaded matrices.
	// This must happen after the basis set has been successfully resolved,
	// since build_ao_remap() requires basis topology to be valid.
	if (vlx->ao_remap == NULL && vlx->basis_set.atom_basis.count > 0) {
		size_t num_ao = build_ao_remap(&vlx->ao_remap, vlx, vlx->arena);
		if (num_ao > 0 && vlx->ao_remap) {
			// SCF: permute+transpose C [num_ao x num_mo] -> [num_mo x num_ao] in shell order
			if (vlx->scf.alpha.coefficients.data && num_ao == vlx->scf.alpha.coefficients.size[0]) {
				size_t num_mo_a = vlx->scf.alpha.coefficients.size[1];
				ao_permute(vlx->scf.alpha.coefficients.data, num_ao, num_mo_a, vlx->ao_remap);
				vlx->scf.alpha.coefficients.size[0] = num_mo_a;
				vlx->scf.alpha.coefficients.size[1] = num_ao;
			}
			if (vlx->scf.alpha.density.data && num_ao == vlx->scf.alpha.density.size[0]) {
				ao_permute_square(vlx->scf.alpha.density.data, num_ao, vlx->ao_remap);
			}
			if (vlx->scf.type == MD_VLX_SCF_TYPE_UNRESTRICTED) {
				if (vlx->scf.beta.coefficients.data && num_ao == vlx->scf.beta.coefficients.size[0]) {
					size_t num_mo_b = vlx->scf.beta.coefficients.size[1];
					ao_permute(vlx->scf.beta.coefficients.data, num_ao, num_mo_b, vlx->ao_remap);
					vlx->scf.beta.coefficients.size[0] = num_mo_b;
					vlx->scf.beta.coefficients.size[1] = num_ao;
				}
				if (vlx->scf.beta.density.data && num_ao == vlx->scf.beta.density.size[0]) {
					ao_permute_square(vlx->scf.beta.density.data, num_ao, vlx->ao_remap);
				}
			}
			if (vlx->scf.S.data && num_ao == vlx->scf.S.size[0]) {
				ao_permute_square(vlx->scf.S.data, num_ao, vlx->ao_remap);
			}
			// NTO: permute+transpose each [num_ao x num_mo] -> [num_mo x num_ao] in shell order
			if (vlx->rsp.nto) {
				for (size_t ni = 0; ni < vlx->rsp.number_of_excited_states; ni++) {
					md_vlx_orbital_t* nto = &vlx->rsp.nto[ni];
					if (!nto->coefficients.data) continue;
					if (num_ao != nto->coefficients.size[0]) continue;
					size_t nmo = nto->coefficients.size[1];
					ao_permute(nto->coefficients.data, num_ao, nmo, vlx->ao_remap);
					nto->coefficients.size[0] = nmo;
					nto->coefficients.size[1] = num_ao;
				}
			}
		}
	}

	// Identify homo and lumo
	if (vlx->scf.alpha.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.alpha.occupancy.size; ++i) {
			if (vlx->scf.alpha.occupancy.data[i] == 0.0) {
				vlx->scf.homo_idx[0] = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.lumo_idx[0] = i;
				break;
			}
		}
	}

	if (vlx->scf.beta.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.beta.occupancy.size; ++i) {
			if (vlx->scf.beta.occupancy.data[i] == 0.0) {
				vlx->scf.homo_idx[1] = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.lumo_idx[1] = i;
				break;
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

	}

	result = true;
done:
	md_temp_set_pos_back(temp_pos);

	return result;
}

static inline void extract_row(double* dst, const md_vlx_2d_data_t* data, size_t row_idx) {
        ASSERT(dst);
        ASSERT(data);
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

static inline size_t number_of_atomic_orbitals(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[0];
}

static inline size_t number_of_ao_coefficients(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline void extract_ao_coefficients(double* out_coeff, const md_vlx_orbital_t* orb, size_t ao_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);
	ASSERT(ao_idx < number_of_atomic_orbitals(orb));

	extract_row(out_coeff, &orb->coefficients, ao_idx);
}

const double* md_vlx_scf_resp_charges(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.resp_charges;
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
			return vlx->rsp.nto[nto_idx].occupancy.data + vlx->scf.lumo_idx[0];
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

size_t md_vlx_vib_number_of_normal_modes(const md_vlx_t* vlx) {
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

size_t md_vlx_atomic_property_count(const md_vlx_t* vlx) {
	if (vlx) {
		return md_array_size(vlx->atomic_properties);
	}
	return 0;
}

const md_vlx_atomic_property_t* md_vlx_atomic_property_by_index(const md_vlx_t* vlx, size_t idx) {
	if (vlx) {
		if (idx < md_array_size(vlx->atomic_properties)) {
			return &vlx->atomic_properties[idx];
		}
	}
	return NULL;
}

const md_vlx_atomic_property_t* md_vlx_atomic_property_by_key(const md_vlx_t* vlx, uint64_t key) {
	if (vlx) {
		for (size_t i = 0; i < md_array_size(vlx->atomic_properties); ++i) {
			if (vlx->atomic_properties[i].key == key) {
				return &vlx->atomic_properties[i];
			}
		}
	}
	return NULL;
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

	for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
		sys->atom.x[i] = (float)vlx->atom_coordinates[i].x;
		sys->atom.y[i] = (float)vlx->atom_coordinates[i].y;
		sys->atom.z[i] = (float)vlx->atom_coordinates[i].z;
		
		md_atomic_number_t z = vlx->atomic_numbers[i];
		str_t sym  = md_atomic_number_symbol(z);
        float mass = md_atomic_number_mass(z);
		float radius = md_atomic_number_vdw_radius(z);
		uint32_t color = md_atomic_number_cpk_color(z);

		md_atom_type_idx_t type_idx = md_atom_type_find_or_add(&sys->atom.type, sym, z, mass, radius, color, 0, sys->alloc);
		sys->atom.type_idx[i] = type_idx;
	}

	sys->atom.count = vlx->number_of_atoms;

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

// Attempt to open the file and check if it can supplement the existing system with QM data
bool md_vlx_system_is_file_supplemental(const md_system_t* sys, str_t filename) {
	ASSERT(sys);

	// Simple check here, we just check for the existence of a couple of fields in h5 file.

	str_t ext = {0};
	if (extract_ext(&ext, filename)) {
        if (!str_eq_ignore_case(ext, STR_LIT("h5")) && !str_eq_ignore_case(ext, STR_LIT("hdf5"))) {
			// Unsupported file extension
			return false;
		}
	}

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
	
    if (h5_check_dataset_exists(file_id, "qm_atom_indices") &&
		h5_check_dataset_exists(file_id, "nuclear_charges"))
	{
        md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
		// Derive atomic numbers from nuclear charges and verify with qm_atom_indices that they are consistent with the loaded system.
        // This is a heuristic check, but it should be sufficient to determine if the file contains data that can supplement the existing system.
		
		size_t vlx_num_atoms = 0;
		if (!h5_read_scalar(&vlx_num_atoms, file_id, H5T_NATIVE_UINT64, "number_of_atoms")) {
			return false;
		}

		md_array(int) qm_atom_indices = md_array_create(int, vlx_num_atoms, temp_arena);
		md_array(int) nuclear_charges = md_array_create(int, vlx_num_atoms, temp_arena);
		
        if (h5_read_dataset_data(qm_atom_indices, vlx_num_atoms, file_id, H5T_NATIVE_INT32, "qm_atom_indices") &&
			h5_read_dataset_data(nuclear_charges, vlx_num_atoms, file_id, H5T_NATIVE_INT32, "nuclear_charges"))
		{
			bool match = true;
			for (size_t i = 0; i < vlx_num_atoms; ++i) {
                int idx = qm_atom_indices[i];
				int z = md_atom_atomic_number(&sys->atom, idx);
                if (z != nuclear_charges[i]) {
					match = false;
					break;
				}
			}
			result = match;
		
		}
		md_arena_allocator_destroy(temp_arena);
	}

done:
	H5Fclose(file_id);

	return result;
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
	if (vlx) return vlx->atom_coordinates;
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

const int* md_vlx_local_to_global_atom_idx(const md_vlx_t* vlx) {
	if (vlx) return vlx->local_to_global_atom_idx;
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

size_t md_vlx_scf_homo_idx(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.homo_idx[0];
		} else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.homo_idx[1];
		}
	}
	return 0;
}

size_t md_vlx_scf_lumo_idx(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.lumo_idx[0];
		} else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.lumo_idx[1];
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

const double* md_vlx_scf_mo_occupancy(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.occupancy.data;
		} 
		else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.beta.occupancy.data;
		}
	}
	return NULL;
}

const double* md_vlx_scf_mo_energy(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.energy.data;
		}
		else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.beta.energy.data;
		}
	}
	return NULL;
}

bool md_vlx_gto_basis_extract(md_gto_basis_t* out, const md_vlx_t* vlx, md_allocator_i* alloc) {
	if (!vlx || !out) return false;
	ASSERT(alloc);

	MEMSET(out, 0, sizeof(*out));

	int natoms   = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	// Emit one shell per contracted function, ordered angl -> atom -> func.
	// MO coefficient vectors stored in vlx are already permuted to this order
	// by build_ao_remap() applied at parse time.
	for (int angl = 0; angl <= max_angl; angl++) {
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			basis_func_t basis_funcs[128];
			size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(
				basis_funcs, ARRAY_SIZE(basis_funcs), &vlx->basis_set, idelem, angl);

			for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
				const basis_func_t* bf = &basis_funcs[funcidx];
				md_gto_shell_t shell = {
					.atom_idx         = (uint32_t)atomidx,
					.primitive_offset = out->num_primitives,
					.l                = (uint16_t)angl,
					.num_primitives   = (uint16_t)bf->count,
				};
				md_array_push(out->shells, shell, alloc);
				out->num_shells++;

				for (int ip = 0; ip < bf->count; ip++) {
					md_array_push(out->alpha, (float)bf->exponents[ip], alloc);
					md_array_push(out->coeff, (float)bf->normalization_coefficients[ip], alloc);
					out->num_primitives++;
				}
			}
		}
	}
	return out->num_shells > 0;
}

// Returns a direct pointer to the AO coefficient vector for MO mo_idx.
// The matrix is stored [num_mo][num_ao] after permutation and transpose at load time,
// so each MO's coefficients are contiguous and in shell order.
const double* md_vlx_scf_mo_coefficients(const md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t type) {
	if (!vlx) return NULL;
	const md_vlx_orbital_t* orb = (type == MD_VLX_SPIN_ALPHA) ? &vlx->scf.alpha :
								  (type == MD_VLX_SPIN_BETA)  ? &vlx->scf.beta  : NULL;
	if (!orb || !orb->coefficients.data) return NULL;
	size_t num_mo = orb->coefficients.size[0];
	size_t num_ao = orb->coefficients.size[1];
	if (mo_idx >= num_mo || num_ao == 0) return NULL;
	return orb->coefficients.data + mo_idx * num_ao;
}

// Returns a direct pointer to the AO coefficient vector for NTO lambda_idx / type.
// Layout is [num_mo][num_ao] after permutation at load time.
const double* md_vlx_rsp_nto_coefficients(const md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type) {
	if (!vlx || !vlx->rsp.nto) return NULL;
	if (nto_idx >= vlx->rsp.number_of_excited_states) return NULL;
	const md_vlx_orbital_t* nto = &vlx->rsp.nto[nto_idx];
	if (!nto->coefficients.data) return NULL;
	size_t num_mo = nto->coefficients.size[0];
	size_t num_ao = nto->coefficients.size[1];
	// NTO particle orbitals are stored in the first half, hole in the second half
	size_t offset = (type == MD_VLX_NTO_TYPE_HOLE) ? (num_mo / 2) : 0;
	size_t idx = offset + lambda_idx;
	if (idx >= num_mo || num_ao == 0) return NULL;
	return nto->coefficients.data + idx * num_ao;
}

// Deprecated extraction wrappers kept for backward compatibility.
size_t md_vlx_scf_mo_coefficients_extract(double* out, const md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t type) {
	const double* src = md_vlx_scf_mo_coefficients(vlx, mo_idx, type);
	if (!src) return 0;
	size_t num_ao = vlx->scf.alpha.coefficients.size[1];
	if (type == MD_VLX_SPIN_BETA) num_ao = vlx->scf.beta.coefficients.size[1];
	if (out) MEMCPY(out, src, sizeof(double) * num_ao);
	return num_ao;
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

// Get the regular density matrix size N in (N x N)
size_t md_vlx_scf_density_matrix_size(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.alpha.density.size[0];
	}
	return 0;
}

// Extracts the full density matrix into a square matrix representation
bool md_vlx_scf_extract_density_matrix_data_float(float* out_values, const struct md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		size_t dim = vlx->scf.alpha.density.size[0];
		const double* density_data = NULL;
		if (type == MD_VLX_SPIN_ALPHA) {
			density_data = vlx->scf.alpha.density.data;
		} else if (type == MD_VLX_SPIN_BETA) {
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

const double* md_vlx_scf_density_matrix_data(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.density.data;
		} else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.beta.density.data;
		} else {
			MD_LOG_ERROR("Invalid MO type for density matrix extraction!");
			return NULL;
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
