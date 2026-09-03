#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <core/md_hash.h>

#include <md_system.h>
#include <md_gto.h>
#include <md_util.h>

#include <hdf5.h>

// Internal only: what an atomic property looks like between reading it out of the h5 file and
// publishing it as an attribute on the system. Nothing outside this file sees it, and once the
// remaining vlx payload moves into the attribute table it goes away entirely.
typedef struct md_vlx_atomic_property_t {
	str_t   label;      // Display text, as authored in the file
	str_t   name;       // Dataset name in the h5 file. The identity, and what the attribute path is built from
	int     num_dims;   // Dimensions exactly as the h5 dataspace gave them: row major, atom axis innermost
	size_t  dim[2];
	double* data;       // num_dims dims worth of values, in that same order
} md_vlx_atomic_property_t;

// Internal only, and for the same reason as md_vlx_atomic_property_t above: what a density
// property looks like between reading it out of the h5 file and publishing it as an attribute.
typedef struct md_vlx_density_property_t {
	str_t    label;     // Display text, as authored in the file
	str_t    name;      // Dataset name in the h5 file. The identity, and what the attribute path is built from
	uint64_t key;
	size_t 	 dim[2];
	double*  data;      // dim[0] * dim[1] values, row major
} md_vlx_density_property_t;

#include <float.h>
#include <math.h>
#include <stdlib.h>	// qsort

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999

#define VLX_NTO_POWER_ITERATIONS 256
#define VLX_NTO_EIGENVALUE_EPSILON 1.0e-14
#define VLX_NTO_CONVERGENCE_EPSILON 1.0e-10

typedef enum {
	VLX_FLAG_CORE = 1,
	VLX_FLAG_SCF  = 2,
	VLX_FLAG_RSP  = 4,
	VLX_FLAG_VIB  = 8,
	VLX_FLAG_OPT  = 16,
	VLX_FLAG_XPS  = 32,
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
	size_t homo_idx;
	size_t lumo_idx;
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

	double energy;
	dvec3_t ground_state_dipole_moment;

	md_vlx_orbital_t alpha;
	md_vlx_orbital_t beta;

	double* resp_charges;

	md_vlx_2d_data_t S;
	md_vlx_scf_history_t history;
} md_vlx_scf_t;

// Internal, parse time only. Offsets rather than pointers, because md_array_push reallocs the entry
// array and would dangle any pointer stored into a group before the array is final.
typedef struct vlx_xps_group_internal_t {
	md_element_t element;
	uint32_t     offset;	// Offset into md_vlx_xps_t::entries
	uint32_t     count;
} vlx_xps_group_internal_t;

typedef struct md_vlx_xps_t {
	md_vlx_xps_entry_t*       entries;			// Flat, sorted by (element, ionization_energy)
	vlx_xps_group_internal_t* groups_internal;	// Built during parse
	md_vlx_xps_group_t*       groups;			// Public view, materialized by vlx_xps_finalize()
} md_vlx_xps_t;

typedef struct md_vlx_rsp_t {
	md_vlx_rsp_type_t type;

	double c6;

	size_t   number_of_frequencies;
	size_t   num_core;
	size_t   num_valence;
	size_t   num_virtual;

	dvec3_t* electric_transition_dipoles;
	dvec3_t* magnetic_transition_dipoles;
	dvec3_t* velocity_transition_dipoles;

	double* frequencies;			// unit = eV
	double* rotatory_strengths;		// unit = 10^-40 cgs
	double* oscillator_strengths;	// Linear and TPA, unitless (peaks to be broadened, absorption)

	double* sigmas;					// CPP only, unit = eV (broadened absorption)
	double* optical_rotations;		// CPP only, unit = deg dm^-1 (broadened optical rotation)
	double* delta_epsilons;			// CPP only, unit = 10^-3 (broadened circular dichroism)

	double* tpa_strengths_linear;	// TPA only
	double* tpa_strengths_circular;	// TPA only
	double* cross_sections;			// TPA only

	// RIXS
	struct {
		// Number of incomming photons (P) == length of the 'photon_energies' dataset
		size_t num_incomming_photons;

		// Number of final (valence) states (F) == first dimension of the 2D datasets below.
		// NOTE: this is generally NOT equal to the number of core-excited states (C), which is
		// what 'number_of_frequencies' holds for RIXS. The 2D arrays are [F][P] row-major.
		size_t num_final_states;

		double* cross_sections;			// [F][P] row-major, a.u.
		double* photon_energies;		// [P], a.u.
		double* elastic_cross_sections;	// [P], a.u.
		double* emission_energies;		// [F][P] row-major, a.u.
		double* energy_losses;			// [F][P] row-major, a.u.
		double* scattering_amplitude_re;
		double* scattering_amplitude_im;

		// These are regular arrays with length num_frequencies (C)
		double* core_eigenvalues;
		double* core_osc_strengths;

		double  gamma_fwhm_ev;
	} rixs;

	// Linear only, Should have dimensions [number_of_frequencies][num_occ * num_vir * 2]
	md_vlx_2d_data_t solution_matrix;
} md_vlx_rsp_t;

typedef struct md_vlx_vib_t {
	size_t number_of_normal_modes;
	double* force_constants;
	double* ir_intensities;
	double* frequencies;
	double* reduced_masses;
	dvec3_t** normal_modes;

	double* tpa_trans_linear;
	double* tpa_trans_circular;

	double* tpa_reduced_gamma_re;
    double* tpa_reduced_gamma_im;
	double* tpa_reduced_cross_section;

    size_t num_external_frequencies;
    double* external_frequencies;

    // Raman activities is multidimensional, with dimensions [number_of_external_frequencies][number_of_normal_modes]
	double* raman_activities;
} md_vlx_vib_t;

typedef struct md_vlx_opt_t {
    md_vlx_opt_type_t type;
	size_t state_index;
	size_t ts_index;
	size_t number_of_steps;
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
	md_vlx_density_property_t* density_properties; // Optional data, may be NULL, length is number of density properties

	// Data blocks
	md_vlx_scf_t scf;
	md_vlx_rsp_t rsp;
	md_vlx_vib_t vib;
	md_vlx_opt_t opt;
	md_vlx_xps_t xps;

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
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(S_num_fac));
		return S_num_fac[isph];
	case 1:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(P_num_fac));
		return P_num_fac[isph];
	case 2:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(D_num_fac));
		return D_num_fac[isph];
	case 3:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(F_num_fac));
		return F_num_fac[isph];
	case 4:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(G_num_fac));
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
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(P_offsets));
		return P_factors + P_offsets[isph];
	case 2:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(D_offsets));
		return D_factors + D_offsets[isph];
	case 3:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(F_offsets));
		return F_factors + F_offsets[isph];
	case 4:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(G_offsets));
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
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(P_offsets));
		return P_indices + P_offsets[isph];
	case 2:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(D_offsets));
		return D_indices + D_offsets[isph];
	case 3:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(F_offsets));
		return F_indices + F_offsets[isph];
	case 4:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(G_offsets));
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

static size_t basis_set_count_atomic_basis_func(const basis_set_t* basis_set, int atomic_number, int angl) {
    size_t count = 0;
    basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
    if (atom_basis) {
        int beg = atom_basis->basis_func_offset;
        int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
        for (int i = beg; i < end; ++i) {
            if (basis_set->basis_func.data[i].type == angl) {
                count++;
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
					if (out_ao_to_atom) {
						out_ao_to_atom[count] = atomidx;
					}
					count += 1;
				}
			}
		}
	}
	return count;
}

static size_t compute_basis_num_atomic_orbitals(const md_vlx_t* vlx) {
	ASSERT(vlx);
	if (!vlx->basis_set.atom_basis.count || !vlx->atomic_numbers || vlx->number_of_atoms == 0) {
		return 0;
	}
	return extract_ao_to_atom_idx(NULL, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set);
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
static bool build_ao_remap(int* out_remap, size_t capacity, const md_vlx_t* vlx) {
	ASSERT(out_remap);
	ASSERT(vlx);

	int natoms   = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	// First, count total AOs
	size_t num_aos = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			size_t num_funcs = basis_set_count_atomic_basis_func(&vlx->basis_set, idelem, angl);
			num_aos += (size_t)nsph * num_funcs;
		}
	}

	if (num_aos != capacity) {
		MD_LOG_ERROR("Capacity of remap table did not match the number of AOs found in the basis set");
		return false;
	}

	MEMSET(out_remap, 0, sizeof(int) * num_aos);

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
	md_temp_scope_t temp = md_temp_begin();
	int* vlx_base = (int*)md_temp_alloc(temp, sizeof(int) * (max_angl + 1) * max_nsph * natoms);
	int* vlx_num  = (int*)md_temp_alloc(temp, sizeof(int) * (max_angl + 1) * natoms);
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
					out_remap[shell_idx] = base + funcidx;
				}
			}
		}
	}

	md_temp_end(temp);
	return true;
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

// HDF5 prints a full diagnostic stack to stderr for every failed call. This reader
// probes for optional data constantly and reports its own failures through MD_LOG_*
// with the offending field name, so the automatic handler is pure noise. Silence it
// for the duration of a parse and restore whatever the host application had set --
// mdlib must not leave global HDF5 state modified.
typedef struct h5_error_scope_t {
	H5E_auto2_t func;
	void*       client_data;
} h5_error_scope_t;

static h5_error_scope_t h5_error_scope_begin(void) {
	h5_error_scope_t scope = {0};
	H5Eget_auto2(H5E_DEFAULT, &scope.func, &scope.client_data);
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	return scope;
}

static void h5_error_scope_end(h5_error_scope_t scope) {
	H5Eset_auto2(H5E_DEFAULT, scope.func, scope.client_data);
}

// H5Lexists() only tolerates a missing *final* path component. If an intermediate
// group is absent it fails outright and, with the default error handler installed,
// dumps a full diagnostic stack to stderr. Nearly every probe in this file is for
// optional data (a file without TPA has no "tpa_strengths" group at all), so walk
// the path one component at a time and stop at the first miss.
//
// Returns false both for "not present" and for a genuine error; callers here only
// care whether the data is usable.
static bool h5_link_exists(hid_t loc_id, const char* path) {
	if (!path || !*path) return false;

	const size_t len = strlen(path);

	char buf[512];
	if (len >= sizeof(buf)) {
		MD_LOG_ERROR("HDF5 path exceeds %zu characters: '%s'", sizeof(buf) - 1, path);
		return false;
	}
	MEMCPY(buf, path, len + 1);

	for (size_t i = 0; i < len; ++i) {
		if (buf[i] != '/') continue;
		if (i == 0 || buf[i - 1] == '/') continue;  // leading or repeated separator

		buf[i] = '\0';
		const htri_t exists = H5Lexists(loc_id, buf, H5P_DEFAULT);
		buf[i] = '/';

		if (exists <= 0) return false;
	}

	return H5Lexists(loc_id, buf, H5P_DEFAULT) > 0;
}

static H5I_type_t h5_get_object_type(hid_t loc_id, const char* name) {
    if (!h5_link_exists(loc_id, name)) {
        return H5I_UNINIT;
    }

    hid_t obj_id = H5Oopen(loc_id, name, H5P_DEFAULT);
    if (obj_id < 0) {
        return H5I_BADID;
    }

    H5I_type_t type = H5Iget_type(obj_id);
    H5Oclose(obj_id);
    return type;
}

// Number of elements in a dataset's dataspace. Returns false on failure.
//
// Every read below uses H5S_ALL for both the memory and file selection, which makes
// H5Dread write *the whole dataset* into the caller's buffer. Any read into a fixed
// size destination must therefore check this first, or a file with an unexpected
// shape silently overruns the destination.
static bool h5_dataset_num_elements(hsize_t* out_count, hid_t dataset_id, const char* field_name) {
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 dataspace for dataset: '%s'", field_name);
		return false;
	}

	const hssize_t npoints = H5Sget_simple_extent_npoints(space_id);
	H5Sclose(space_id);

	if (npoints < 0) {
		MD_LOG_ERROR("Failed to query element count for H5 dataset: '%s'", field_name);
		return false;
	}

	*out_count = (hsize_t)npoints;
	return true;
}

// Reads a dataset expected to hold exactly one element. Returns false if it holds
// anything else, rather than overrunning 'buf'.
static bool h5_read_scalar(void* buf, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	if (!h5_link_exists(file_id, field_name)) {
		return false;
	}

	// Open the dataset containing the double value
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool result = false;

	hsize_t num_elem = 0;
	if (!h5_dataset_num_elements(&num_elem, dataset_id, field_name)) {
		goto done;
	}
	if (num_elem != 1) {
		MD_LOG_ERROR("Expected a scalar H5 dataset for '%s', got %llu elements", field_name, (unsigned long long)num_elem);
		goto done;
	}

	// Read the dataset into the 'value' variable
	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if (status < 0) {
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

	if (!h5_link_exists(file_id, field_name)) {
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

	{
		// Both branches below read with H5S_ALL into storage for a single string.
		hsize_t num_elem = 0;
		if (!h5_dataset_num_elements(&num_elem, dataset_id, field_name)) {
			goto done;
		}
		if (num_elem != 1) {
			MD_LOG_ERROR("Expected a single string in H5 dataset '%s', got %llu", field_name, (unsigned long long)num_elem);
			goto done;
		}
	}

	if (H5Tis_variable_str(datatype_id)) {
		char* tmp = NULL;
		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
		if (status < 0) {
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
		bool fixed_result = false;
		const size_t raw_len = H5Tget_size(datatype_id);
		md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
		char* tmp = md_temp_alloc(temp_scope, raw_len + 1);
		if (!tmp) {
			MD_LOG_ERROR("Failed to allocate temporary buffer for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}
		MEMSET(tmp, 0, raw_len + 1);

		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read fixed-length string for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}

		size_t len = strnlen(tmp, raw_len);
		str_t data = str_alloc(len, alloc);
		if (len > 0) {
			MEMCPY((char*)data.ptr, tmp, len);
		}
		*str = data;
		fixed_result = true;

	fixed_done:
		md_temp_end(temp_scope);
		if (!fixed_result) goto done;
	}

	result = true;

done:
	if (datatype_id != H5I_INVALID_HID) H5Tclose(datatype_id);
	if (space_id != H5I_INVALID_HID) H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static size_t h5_read_cstr(char* out_str, size_t str_cap, hid_t file_id, const char* field_name) {
    ASSERT(out_str);
	ASSERT(str_cap > 0);
	out_str[0] = '\0';

	if (!h5_link_exists(file_id, field_name)) {
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
	if (datatype_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 datatype for dataset: '%s'", field_name);
		goto done;
	}

	if (H5Tget_class(datatype_id) != H5T_STRING) {
		MD_LOG_ERROR("H5 dataset is not a string: '%s'", field_name);
		goto done;
	}

	{
		// Both branches below read with H5S_ALL into storage for a single string.
		hsize_t num_elem = 0;
		if (!h5_dataset_num_elements(&num_elem, dataset_id, field_name)) {
			goto done;
		}
		if (num_elem != 1) {
			MD_LOG_ERROR("Expected a single string in H5 dataset '%s', got %llu", field_name, (unsigned long long)num_elem);
			goto done;
		}
	}

	if (H5Tis_variable_str(datatype_id)) {
		// Variable-length string
		char* tmp = NULL;
		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read variable-length string for H5 dataset: '%s'", field_name);
			goto done;
		}

		if (tmp) {
			size_t len = strlen(tmp);
			result = MIN(len, str_cap - 1);
			MEMCPY(out_str, tmp, result);
			out_str[result] = '\0';
			H5free_memory(tmp);
		}
	} else {
      // Fixed-length strings may not be null-terminated, so read into a temporary buffer first.
		bool fixed_result = false;
		size_t size = H5Tget_size(datatype_id);
		md_temp_scope_t temp = md_temp_begin();

		char* tmp = md_temp_alloc(temp, size + 1);
		if (!tmp) {
			MD_LOG_ERROR("Failed to allocate temporary buffer for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}
		MEMSET(tmp, 0, size + 1);

		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read fixed-length string for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}

		size_t len = strnlen(tmp, size);
		result = MIN(len, str_cap - 1);
		MEMCPY(out_str, tmp, result);
		out_str[result] = '\0';
		fixed_result = true;

	fixed_done:
		md_temp_end(temp);
		if (!fixed_result) goto done;
	}

done:
	// Close HDF5 resources
  if (datatype_id != H5I_INVALID_HID) H5Tclose(datatype_id);
	H5Dclose(dataset_id);

	return result;
}

#define H5_MAX_RANK 32

// Writes the extent of each dimension into dims[0 .. rank-1] and returns the rank.
//
// Returns 0 for "not present" and for every failure. Never returns a negative value:
// most call sites test this truthily, and a -1 would read as success.
static int h5_read_dataset_dims(size_t* dims, int max_dims, hid_t file_id, const char* field_name) {
	ASSERT(dims);

	if (max_dims <= 0 || max_dims > H5_MAX_RANK) {
		MD_LOG_ERROR("Invalid max_dims (%i) requested for H5 dataset: '%s'", max_dims, field_name);
		return 0;
	}

	if (!h5_link_exists(file_id, field_name)) {
		return 0;
	}

	// Open the dataset
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return 0;
	}

	int   result   = 0;
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 dataspace for dataset: '%s'", field_name);
		goto done;
	}

	int ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim < 0) {
		MD_LOG_ERROR("Failed to get number of dimensions for H5 dataset: '%s'", field_name);
		goto done;
	}

	if (ndim > max_dims) {
		MD_LOG_ERROR("H5 dataset '%s' has rank %i, caller supplied room for %i", field_name, ndim, max_dims);
		goto done;
	}

	// hsize_t is always 64-bit while size_t is not, so read into a correctly typed
	// buffer and widen/narrow explicitly rather than aliasing the caller's array.
	hsize_t extent[H5_MAX_RANK] = {0};
	ndim = H5Sget_simple_extent_dims(space_id, extent, NULL);
	if (ndim < 0) {
		MD_LOG_ERROR("Failed to get dimensions for H5 dataset: '%s'", field_name);
		goto done;
	}

	for (int i = 0; i < ndim; ++i) {
		if (extent[i] > (hsize_t)SIZE_MAX) {
			MD_LOG_ERROR("H5 dataset '%s' dimension %i does not fit in size_t", field_name, i);
			goto done;
		}
		dims[i] = (size_t)extent[i];
	}

	result = ndim;

done:
	if (space_id != H5I_INVALID_HID) H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static bool h5_check_dataset_exists(hid_t file_id, const char* field_name) {
	return h5_link_exists(file_id, field_name);
}

// Reads a string attribute into a caller-supplied buffer, always null terminated
// and always truncated to fit. Returns false if the attribute is absent, is not a
// string, or cannot be read; every handle opened here is released on every path.
//
// The size check on the fixed-length branch matters: the stored size is whatever
// the writer chose, so reading it straight into the destination overflows for any
// attribute longer than the caller's buffer.
static bool h5_read_string_attribute(char* out_buf, size_t out_cap, hid_t obj_id, const char* attr_name) {
	ASSERT(out_buf);
	ASSERT(out_cap > 0);
	out_buf[0] = '\0';

	// htri_t: negative on error, and negative is truthy.
	if (H5Aexists(obj_id, attr_name) <= 0) {
		return false;
	}

	hid_t attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT);
	if (attr_id == H5I_INVALID_HID) {
		return false;
	}

	bool  result    = false;
	hid_t attr_type = H5Aget_type(attr_id);
	if (attr_type == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query type of attribute '%s'", attr_name);
		goto done;
	}

	if (H5Tget_class(attr_type) != H5T_STRING) {
		MD_LOG_ERROR("Attribute '%s' is not a string", attr_name);
		goto done;
	}

	if (H5Tis_variable_str(attr_type)) {
		char* var_str = NULL;
		if (H5Aread(attr_id, attr_type, &var_str) < 0) {
			MD_LOG_ERROR("Failed to read variable-length string attribute '%s'", attr_name);
			goto done;
		}
		// On success HDF5 may still hand back NULL for an empty string.
		if (var_str) {
			const size_t len = MIN(strlen(var_str), out_cap - 1);
			MEMCPY(out_buf, var_str, len);
			out_buf[len] = '\0';
			H5free_memory(var_str);
		}
		result = true;
	} else {
		const size_t size = H5Tget_size(attr_type);
		if (size == 0) {
			MD_LOG_ERROR("Attribute '%s' has zero size", attr_name);
			goto done;
		}

		// Fixed-length strings are not necessarily null terminated.
		md_temp_scope_t temp = md_temp_begin();
		char* tmp = md_temp_alloc(temp, size + 1);
		if (tmp) {
			MEMSET(tmp, 0, size + 1);
			if (H5Aread(attr_id, attr_type, tmp) >= 0) {
				const size_t len = MIN(strnlen(tmp, size), out_cap - 1);
				MEMCPY(out_buf, tmp, len);
				out_buf[len] = '\0';
				result = true;
			} else {
				MD_LOG_ERROR("Failed to read fixed-length string attribute '%s'", attr_name);
			}
		}
		md_temp_end(temp);
	}

done:
	if (attr_type != H5I_INVALID_HID) H5Tclose(attr_type);
	H5Aclose(attr_id);
	return result;
}

// Checks a square AO-basis matrix for symmetry and forces it if it is close but not
// exact. Returns true if it was already symmetric to tolerance.
//
// Consumers of AO density matrices in this codebase read only the upper triangle, so
// an asymmetric matrix is not merely inaccurate -- half of it is discarded without a
// trace. Anything beyond rounding is reported with the offending magnitude so it is
// visible rather than absorbed.
static bool vlx_report_and_enforce_symmetry(double* mat, size_t dim, const char* label) {
	ASSERT(mat);

	double max_asym = 0.0;
	double max_abs  = 0.0;
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = i + 1; j < dim; ++j) {
			const double a = mat[i * dim + j];
			const double b = mat[j * dim + i];
			max_asym = MAX(max_asym, fabs(a - b));
			max_abs  = MAX(max_abs, MAX(fabs(a), fabs(b)));
		}
	}

	// Scale-relative, so this does not fire on accumulated rounding in a large matrix.
	const double tolerance = 1.0e-10 * (max_abs > 0.0 ? max_abs : 1.0);
	if (max_asym <= tolerance) {
		return true;
	}

	MD_LOG_INFO("Density property '%s' is not symmetric (max deviation %g, largest element %g). "
				"Symmetrizing: the density evaluation path only reads the upper triangle.",
				label, max_asym, max_abs);

	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = i + 1; j < dim; ++j) {
			const double value = 0.5 * (mat[i * dim + j] + mat[j * dim + i]);
			mat[i * dim + j] = value;
			mat[j * dim + i] = value;
		}
	}
	return false;
}

typedef bool (*h5_group_visit_cb_t)(md_vlx_t* vlx, hid_t group_handle, const char* group_path, void* user_data);

// Guards against pathological or cyclic (soft/external link) hierarchies. VeloxChem
// files nest a handful of levels; anything deeper is not something we should follow.
#define H5_MAX_GROUP_DEPTH 32

static bool h5_visit_groups_recursive_impl(md_vlx_t* vlx, hid_t group_handle, const char* group_path, h5_group_visit_cb_t callback, void* user_data, int depth) {
	ASSERT(vlx);
	ASSERT(group_path);
	ASSERT(callback);

	if (!callback(vlx, group_handle, group_path, user_data)) {
		return false;
	}

	if (depth >= H5_MAX_GROUP_DEPTH) {
		MD_LOG_ERROR("HDF5 group nesting exceeds %i levels at '%s', not descending further", H5_MAX_GROUP_DEPTH, group_path);
		return true;
	}

	H5G_info_t info = { 0 };
	if (H5Gget_info(group_handle, &info) < 0) {
		MD_LOG_ERROR("Failed to get group info when traversing HDF5 groups");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		// H5Gget_objname_by_idx / H5Gget_objtype_by_idx are the deprecated 1.6 API and
		// are compiled out when HDF5 is built without the compatibility layer.
		// H5Lget_name_by_idx returns the length excluding the terminator, or negative.
		const ssize_t name_len = H5Lget_name_by_idx(group_handle, ".", H5_INDEX_NAME, H5_ITER_INC, i, name_buf, sizeof(name_buf), H5P_DEFAULT);
		if (name_len < 0) {
			continue;
		}
		if ((size_t)name_len >= sizeof(name_buf)) {
			// Truncated: opening the shortened name would resolve to the wrong object.
			MD_LOG_ERROR("Skipping HDF5 link under '%s' whose name exceeds %zu characters", group_path, sizeof(name_buf) - 1);
			continue;
		}

		if (h5_get_object_type(group_handle, name_buf) != H5I_GROUP) {
			continue;
		}

		hid_t child_group = H5Gopen(group_handle, name_buf, H5P_DEFAULT);
		if (child_group == H5I_INVALID_HID) {
			continue;
		}

		char child_path[512];
		if (strcmp(group_path, "/") == 0) {
			snprintf(child_path, sizeof(child_path), "/%s", name_buf);
		} else {
			snprintf(child_path, sizeof(child_path), "%s/%s", group_path, name_buf);
		}

		bool result = h5_visit_groups_recursive_impl(vlx, child_group, child_path, callback, user_data, depth + 1);
		H5Gclose(child_group);
		if (!result) {
			return false;
		}
	}

	return true;
}

static bool h5_visit_groups_recursive(md_vlx_t* vlx, hid_t group_handle, const char* group_path, h5_group_visit_cb_t callback, void* user_data) {
	return h5_visit_groups_recursive_impl(vlx, group_handle, group_path, callback, user_data, 0);
}

static bool h5_read_dataset_data(void* out_data, size_t num_samples, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	ASSERT(out_data);

	if (!h5_link_exists(file_id, field_name)) {
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

	const hssize_t num_points = H5Sget_simple_extent_npoints(space_id);
	if (num_points < 0) {
		MD_LOG_ERROR("Failed to query element count for H5 dataset: '%s'", field_name);
		goto done;
	}

	if ((hsize_t)num_points != (hsize_t)num_samples) {
		MD_LOG_ERROR("Unexpected number of points reading H5 dataset '%s', got %llu, expected %zu",
			field_name, (unsigned long long)num_points, num_samples);
		goto done;
	}

	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_data);

	if (status < 0) {
		MD_LOG_ERROR("Failed to read H5 dataset: '%s'", field_name);
		goto done;
	}

	result = true;
done:
	H5Sclose(space_id);
	H5Dclose(dataset_id);

    return result;
}

// HDF5 has no native complex type. h5py stores numpy complex arrays as a compound type with two
// floating point members, named 'r' and 'i' by default (configurable via h5py's 'complex_names'),
// laid out as interleaved (re, im) pairs on disk.
//
// This reads such a dataset and splits it into two separate, tightly packed arrays of doubles.
// HDF5 does the de-interleaving for us: a memory datatype that declares only one of the two members
// makes H5Dread gather just that component, so no interleaved temporary buffer is needed.
//
// out_real and out_imag must each hold num_samples doubles. Either may be NULL to skip that
// component. Plain (non-compound) real datasets are also accepted, in which case out_imag is zeroed.
static bool h5_read_complex_dataset_split(double* out_real, double* out_imag, size_t num_samples, hid_t file_id, const char* field_name) {
	if (!h5_link_exists(file_id, field_name)) {
		return false;
	}

	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool  result       = false;
	hid_t space_id     = H5I_INVALID_HID;
	hid_t file_type_id = H5I_INVALID_HID;
	hid_t real_type_id = H5I_INVALID_HID;
	hid_t imag_type_id = H5I_INVALID_HID;
	char* member_name[2] = { NULL, NULL };

	space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 space for dataset: '%s'", field_name);
		goto done;
	}

	hsize_t num_points = H5Sget_simple_extent_npoints(space_id);
	if (num_points != (hsize_t)num_samples) {
		MD_LOG_ERROR("Unexpected number of points when reading dataset '%s', got %i, expected %i", field_name, (int)num_points, (int)num_samples);
		goto done;
	}

	file_type_id = H5Dget_type(dataset_id);
	if (file_type_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to get H5 datatype for dataset: '%s'", field_name);
		goto done;
	}

	if (H5Tget_class(file_type_id) != H5T_COMPOUND) {
		// Not stored as complex, treat the data as purely real.
		if (out_real) {
			if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_real) < 0) {
				MD_LOG_ERROR("An error occured when reading H5 dataset: '%s'", field_name);
				goto done;
			}
		}
		if (out_imag) {
			MEMSET(out_imag, 0, num_samples * sizeof(double));
		}
		result = true;
		goto done;
	}

	if (H5Tget_nmembers(file_type_id) != 2) {
		MD_LOG_ERROR("Expected 2 members in the compound H5 datatype of dataset '%s'", field_name);
		goto done;
	}

	// The member names depend on how the file was written, so query them instead of assuming 'r'/'i'.
	member_name[0] = H5Tget_member_name(file_type_id, 0);
	member_name[1] = H5Tget_member_name(file_type_id, 1);
	if (!member_name[0] || !member_name[1]) {
		MD_LOG_ERROR("Failed to get the compound member names of H5 dataset: '%s'", field_name);
		goto done;
	}

	// Members are conventionally ordered (real, imaginary), but do not rely on it.
	int real_idx = (member_name[0][0] == 'i' || member_name[0][0] == 'I') ? 1 : 0;
	int imag_idx = 1 - real_idx;

	if (out_real) {
		real_type_id = H5Tcreate(H5T_COMPOUND, sizeof(double));
		if (real_type_id == H5I_INVALID_HID || H5Tinsert(real_type_id, member_name[real_idx], 0, H5T_NATIVE_DOUBLE) < 0) {
			MD_LOG_ERROR("Failed to create a memory datatype for the real part of H5 dataset: '%s'", field_name);
			goto done;
		}
		if (H5Dread(dataset_id, real_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_real) < 0) {
			MD_LOG_ERROR("An error occured when reading the real part of H5 dataset: '%s'", field_name);
			goto done;
		}
	}

	if (out_imag) {
		imag_type_id = H5Tcreate(H5T_COMPOUND, sizeof(double));
		if (imag_type_id == H5I_INVALID_HID || H5Tinsert(imag_type_id, member_name[imag_idx], 0, H5T_NATIVE_DOUBLE) < 0) {
			MD_LOG_ERROR("Failed to create a memory datatype for the imaginary part of H5 dataset: '%s'", field_name);
			goto done;
		}
		if (H5Dread(dataset_id, imag_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_imag) < 0) {
			MD_LOG_ERROR("An error occured when reading the imaginary part of H5 dataset: '%s'", field_name);
			goto done;
		}
	}

	result = true;
done:
	if (member_name[0]) H5free_memory(member_name[0]);
	if (member_name[1]) H5free_memory(member_name[1]);
	if (real_type_id != H5I_INVALID_HID) H5Tclose(real_type_id);
	if (imag_type_id != H5I_INVALID_HID) H5Tclose(imag_type_id);
	if (file_type_id != H5I_INVALID_HID) H5Tclose(file_type_id);
	if (space_id     != H5I_INVALID_HID) H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static bool h5_read_atomic_properties_in_group(md_vlx_t* vlx, hid_t group_handle, const char* group_path, void* user_data) {
	(void)group_path;
	(void)user_data;

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

		// Ensure that the type is a dataset, if not we skip
		if (type != H5G_DATASET) {
			continue;
		}
		hid_t dataset_id = H5Dopen(group_handle, name_buf, H5P_DEFAULT);
		if (dataset_id == H5I_INVALID_HID) {
			continue;
		}

		// The presence of this attribute is what marks a dataset as an atomic property.
		// A dataset that lacks it, or whose attribute is unreadable, is simply skipped
		// rather than aborting the whole traversal.
		char property_label[256] = { 0 };
		if (!h5_read_string_attribute(property_label, sizeof(property_label), dataset_id, "atomic_property")) {
			H5Dclose(dataset_id);
			continue;
		}

		if (property_label[0] == '\0') {
			// Not an error, just use the field name as the property label
			snprintf(property_label, sizeof(property_label), "%s", name_buf);
		}

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
		if (num_dims < 1 || num_dims > 2) {
			MD_LOG_ERROR("Unsupported rank for atomic property dataset '%s', expected 1 or 2, got %i", name_buf, num_dims);
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

		// The dimensions are kept exactly as the dataspace reported them. Re-spelling them here -
		// atoms first, variants second - is what previously left the struct disagreeing with the
		// buffer it describes, and every reader of it having to know that.
		md_vlx_atomic_property_t property = {
			 .label = str_copy_cstr(property_label, vlx->arena),
			 .name = str_copy_cstr(name_buf, vlx->arena),
			 .num_dims = num_dims,
			 .data = NULL,
		};
		for (int d = 0; d < num_dims; ++d) {
			property.dim[d] = dims[d];
		}

		md_array_resize(property.data, num_points, vlx->arena);
		herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, property.data);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read data for atomic property dataset '%s'", name_buf);
			goto done;
		}

		md_array_push(vlx->atomic_properties, property, vlx->arena);
	done:
		// The attribute handles are owned and released by h5_read_string_attribute().
		H5Dclose(dataset_id);
	}
	return true;
}

static bool h5_read_density_properties_in_group(md_vlx_t* vlx, hid_t group_handle, const char* group_path, void* user_data) {
	(void)group_path;
	(void)user_data;

	H5G_info_t info = { 0 };
	if (H5Gget_info(group_handle, &info) < 0) {
		MD_LOG_ERROR("Failed to get group info when reading density properties");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		ssize_t size = H5Gget_objname_by_idx(group_handle, i, name_buf, sizeof(name_buf));
		if (size < 0) {
			continue;
		}
		H5G_obj_t type = H5Gget_objtype_by_idx(group_handle, i);

		// Ensure that the type is a dataset, if not we skip
		if (type != H5G_DATASET) {
			continue;
		}
		hid_t dataset_id = H5Dopen(group_handle, name_buf, H5P_DEFAULT);
		if (dataset_id == H5I_INVALID_HID) {
			continue;
		}

		// As above: the attribute is the marker, so skip rather than abort.
		char property_label[256] = { 0 };
		if (!h5_read_string_attribute(property_label, sizeof(property_label), dataset_id, "density_property")) {
			H5Dclose(dataset_id);
			continue;
		}

		if (property_label[0] == '\0') {
			snprintf(property_label, sizeof(property_label), "%s", name_buf);
		}

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

		if (num_dims != 2) {
			MD_LOG_ERROR("Too many dimensions for density property dataset '%s', expected 2, got %i", name_buf, num_dims);
			H5Sclose(space_id);
			goto done;
		}

		size_t dims[2] = { 0 };
		H5Sget_simple_extent_dims(space_id, (hsize_t*)dims, 0);

		size_t num_aos = md_vlx_scf_number_of_atomic_orbitals(vlx);
		if (dims[0] != num_aos || dims[1] != num_aos) {
			MD_LOG_ERROR("Unexpected dimensions for density property dataset '%s', expected [%zu x %zu], got [%zu x %zu]", name_buf, num_aos, num_aos, dims[0], dims[1]);
			H5Sclose(space_id);
			continue;
		}

		size_t num_points = H5Sget_simple_extent_npoints(space_id);
		
		H5Sclose(space_id);

		// Construct a unique uint64_t key for this property.
		uint64_t key = md_hash64(name_buf, sizeof(name_buf), 0);
		
		md_vlx_density_property_t property = {
			 .label = str_copy_cstr(property_label, vlx->arena),
			 .name = str_copy_cstr(name_buf, vlx->arena),
			 .key = key,
			 .dim[0] = dims[0],
			 .dim[1] = dims[1],
			 .data = NULL,
		};

		md_array_resize(property.data, num_points, vlx->arena);
		herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, property.data);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read data for density property dataset '%s'", name_buf);
			goto done;
		}

		// Symmetry is a load-bearing assumption downstream: the GL/GPU density path
		// packs only the upper triangle (density_matrix_upper_tri_extract_float in
		// md_gto.c) and the lower half is never read. The SCF and transition densities
		// are symmetric by construction -- the latter is explicitly symmetrized in
		// vlx_rsp_extract_transition_density -- but density properties are a generic
		// pass-through from the file, so nothing has checked them until here.
		// Report and enforce rather than letting half the matrix be silently dropped.
		if (dims[0] == dims[1] && dims[0] > 1) {
			vlx_report_and_enforce_symmetry(property.data, dims[0], property_label);
		}

		md_array_push(vlx->density_properties, property, vlx->arena);
		MD_LOG_DEBUG("Read density property '%s' with dimensions [%zu x %zu]", property_label, dims[0], dims[1]);
	done:
		// The attribute handles are owned and released by h5_read_string_attribute().
		H5Dclose(dataset_id);
	}
	return true;
}

static bool h5_read_atomic_properties(md_vlx_t* vlx, hid_t group_handle) {
	return h5_visit_groups_recursive(vlx, group_handle, "/", h5_read_atomic_properties_in_group, NULL);
}

static bool h5_read_density_properties(md_vlx_t* vlx, hid_t group_handle) {
	return h5_visit_groups_recursive(vlx, group_handle, "/", h5_read_density_properties_in_group, NULL);
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
	md_temp_scope_t temp = md_temp_begin();
	double* tmp = md_temp_alloc_array(temp, double, num_ao * num_mo);
    if (!tmp) {
        MD_LOG_ERROR("Failed to allocate temporary buffer for AO permutation");
		goto done;
    }
	MEMCPY(tmp, mat, sizeof(double) * num_ao * num_mo);
	for (size_t mo = 0; mo < num_mo; mo++) {
		for (size_t ao = 0; ao < num_ao; ao++) {
			mat[mo * num_ao + ao] = tmp[(size_t)remap[ao] * num_mo + mo];
		}
	}
done:
	md_temp_end(temp);
}


// Permute both rows and columns of a square AO×AO matrix (num_ao × num_ao, row-major).
static void ao_permute_square(double* mat, size_t num_ao, const int* remap) {
	md_temp_scope_t temp = md_temp_begin();
	double* tmp = md_temp_alloc_array(temp, double, num_ao * num_ao);
    if (!tmp) {
        MD_LOG_ERROR("Failed to allocate temporary buffer for AO permutation");
		goto done;
    }
	MEMCPY(tmp, mat, sizeof(double) * num_ao * num_ao);
	for (size_t i = 0; i < num_ao; i++) {
		size_t si = (size_t)remap[i];
		for (size_t j = 0; j < num_ao; j++) {
			size_t sj = (size_t)remap[j];
			mat[i * num_ao + j] = tmp[si * num_ao + sj];
		}
	}
done:
	md_temp_end(temp);
}

// Permute columns of a [num_mo x num_ao] matrix according to remap.
static void ao_permute_cols(double* mat, size_t num_mo, size_t num_ao, const int* remap) {
	md_temp_scope_t temp = md_temp_begin();
	double* tmp = md_temp_alloc_array(temp, double, num_mo * num_ao);
	if (!tmp) {
		MD_LOG_ERROR("Failed to allocate temporary buffer for AO permutation");
		goto done;
	}
	MEMCPY(tmp, mat, sizeof(double) * num_mo * num_ao);
	for (size_t mo = 0; mo < num_mo; ++mo) {
		for (size_t ao = 0; ao < num_ao; ++ao) {
			mat[mo * num_ao + ao] = tmp[mo * num_ao + (size_t)remap[ao]];
		}
	}
done:
	md_temp_end(temp);
}

// ---------------------------------------------------------------------------
// Pure/spherical -> Cartesian AO conversion
// ---------------------------------------------------------------------------
// VeloxChem stores AO data in the pure (2l+1) basis. md_gto_basis_t consumes the
// Cartesian ((l+1)(l+2)/2) basis -- see the AO CONVENTION block in md_gto.h.
// Everything AO-indexed is converted once, here, at the end of parsing.
//
// ORDER MATTERS. This must run *after* ao_permute*() has put the matrices in
// shell order, and after every piece of format-internal math that touches the
// AO basis. The Cartesian embedding is rank deficient (the s-type contaminant of
// a d shell is unoccupied), so a converted overlap matrix is singular: anything
// that inverts or Lowdin-orthogonalizes S must have run already.

// [num_mo][n_sph] -> [num_mo][n_cart], reallocated from 'arena'.
static bool vlx_cart_convert_coeff(md_vlx_2d_data_t* mat, const md_gto_basis_t* basis,
	size_t n_sph, size_t n_cart, md_allocator_i* arena, const char* label)
{
	if (!mat->data) return true;

	if (mat->size[1] != n_sph) {
		MD_LOG_ERROR("%s: expected %zu spherical AO columns, got %zu", label, n_sph, mat->size[1]);
		return false;
	}

	const size_t num_mo = mat->size[0];
	double* dst = (double*)md_alloc(arena, sizeof(double) * num_mo * n_cart);
	if (!dst) {
		MD_LOG_ERROR("%s: failed to allocate Cartesian coefficient matrix", label);
		return false;
	}

	for (size_t mo = 0; mo < num_mo; ++mo) {
		if (md_gto_sph_to_cart_vector(dst + mo * n_cart, mat->data + mo * n_sph, basis) != n_cart) {
			MD_LOG_ERROR("%s: spherical to Cartesian conversion failed for MO %zu", label, mo);
			return false;
		}
	}

	mat->data    = dst;
	mat->size[1] = n_cart;
	return true;
}

// [n_sph][n_sph] -> [n_cart][n_cart], reallocated from 'arena'.
static bool vlx_cart_convert_square(md_vlx_2d_data_t* mat, const md_gto_basis_t* basis,
	size_t n_sph, size_t n_cart, md_allocator_i* arena, const char* label)
{
	if (!mat->data) return true;

	if (mat->size[0] != n_sph || mat->size[1] != n_sph) {
		MD_LOG_ERROR("%s: expected [%zu x %zu] spherical AO matrix, got [%zu x %zu]",
			label, n_sph, n_sph, mat->size[0], mat->size[1]);
		return false;
	}

	double* dst = (double*)md_alloc(arena, sizeof(double) * n_cart * n_cart);
	if (!dst) {
		MD_LOG_ERROR("%s: failed to allocate Cartesian matrix", label);
		return false;
	}

	if (md_gto_sph_to_cart_matrix(dst, mat->data, basis) != n_cart) {
		MD_LOG_ERROR("%s: spherical to Cartesian matrix conversion failed", label);
		return false;
	}

	mat->data    = dst;
	mat->size[0] = n_cart;
	mat->size[1] = n_cart;
	return true;
}

static bool vlx_convert_ao_data_to_cartesian(md_vlx_t* vlx) {
	md_temp_scope_t temp = md_temp_begin();
	md_allocator_i* temp_alloc = md_temp_allocator(temp);
	bool result = false;

	md_gto_basis_t basis = {0};
	if (!md_vlx_gto_basis_extract(&basis, vlx, temp_alloc)) {
		MD_LOG_ERROR("Failed to extract GTO basis for Cartesian AO conversion");
		goto done;
	}

	const size_t n_sph  = md_gto_basis_num_sph_ao(&basis);
	const size_t n_cart = md_gto_basis_num_ao(&basis);
	if (n_sph == 0 || n_cart == 0) {
		MD_LOG_ERROR("Cartesian AO conversion: empty basis (n_sph=%zu n_cart=%zu)", n_sph, n_cart);
		goto done;
	}

	// In the restricted case the beta orbital shares alpha's buffers (see the
	// struct memcpy in the parse path). Detect that so we convert once and
	// re-alias rather than converting the same memory twice.
	const bool beta_aliases_alpha =
		(vlx->scf.beta.coefficients.data != NULL &&
		 vlx->scf.beta.coefficients.data == vlx->scf.alpha.coefficients.data);

	if (!vlx_cart_convert_coeff(&vlx->scf.alpha.coefficients, &basis, n_sph, n_cart, vlx->arena, "Alpha orbital coefficients")) goto done;
	if (!vlx_cart_convert_square(&vlx->scf.alpha.density,     &basis, n_sph, n_cart, vlx->arena, "Alpha density")) goto done;

	if (beta_aliases_alpha) {
		MEMCPY(&vlx->scf.beta.coefficients, &vlx->scf.alpha.coefficients, sizeof(md_vlx_2d_data_t));
		MEMCPY(&vlx->scf.beta.density,      &vlx->scf.alpha.density,      sizeof(md_vlx_2d_data_t));
	} else {
		if (!vlx_cart_convert_coeff(&vlx->scf.beta.coefficients, &basis, n_sph, n_cart, vlx->arena, "Beta orbital coefficients")) goto done;
		if (!vlx_cart_convert_square(&vlx->scf.beta.density,     &basis, n_sph, n_cart, vlx->arena, "Beta density")) goto done;
	}

	if (!vlx_cart_convert_square(&vlx->scf.S, &basis, n_sph, n_cart, vlx->arena, "SCF overlap")) goto done;

	// Density properties are AO-basis [N][N] matrices read straight from the file.
	for (size_t i = 0; i < md_array_size(vlx->density_properties); ++i) {
		md_vlx_density_property_t* prop = &vlx->density_properties[i];
		if (!prop->data) continue;

		md_vlx_2d_data_t view = { .size = { prop->dim[0], prop->dim[1] }, .data = prop->data };
		if (!vlx_cart_convert_square(&view, &basis, n_sph, n_cart, vlx->arena, "Density property")) goto done;

		prop->data   = view.data;
		prop->dim[0] = view.size[0];
		prop->dim[1] = view.size[1];
	}

	// Derive the AO -> atom map from the shell list, so it cannot drift out of
	// step with the AO ordering the evaluator walks.
	md_array_resize(vlx->ao_to_atom_idx, n_cart, vlx->arena);
	{
		uint32_t* tmp_map = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * n_cart);
		if (!tmp_map || md_gto_basis_ao_to_atom(tmp_map, &basis) != n_cart) {
			MD_LOG_ERROR("Failed to derive the AO to atom map");
			goto done;
		}
		for (size_t i = 0; i < n_cart; ++i) {
			vlx->ao_to_atom_idx[i] = (int)tmp_map[i];
		}
	}

	result = true;
done:
	md_temp_end(temp);
	return result;
}

static bool validate_square_matrix_dims(const md_vlx_2d_data_t* data, const char* label) {
	ASSERT(data);
	ASSERT(label);

	if (data->size[0] == 0 || data->size[1] == 0) {
		MD_LOG_ERROR("%s matrix has invalid dimensions [%zu x %zu]", label, data->size[0], data->size[1]);
		return false;
	}
	if (data->size[0] != data->size[1]) {
		MD_LOG_ERROR("%s matrix must be square, got [%zu x %zu]", label, data->size[0], data->size[1]);
		return false;
	}
	return true;
}

static bool infer_num_mo_from_coeff_dims(size_t* num_mo, const size_t coeff_dim[2], size_t num_ao, const char* label) {
	ASSERT(num_mo);
	ASSERT(coeff_dim);
	ASSERT(label);

	if (coeff_dim[0] == 0 || coeff_dim[1] == 0) {
		MD_LOG_ERROR("%s coefficient matrix has invalid dimensions [%zu x %zu]", label, coeff_dim[0], coeff_dim[1]);
		return false;
	}

	if (coeff_dim[0] == num_ao) {
		*num_mo = coeff_dim[1];
		return true;
	}
	if (coeff_dim[1] == num_ao) {
		*num_mo = coeff_dim[0];
		return true;
	}

	MD_LOG_ERROR("%s coefficient matrix [%zu x %zu] does not contain AO dimension %zu", label, coeff_dim[0], coeff_dim[1], num_ao);
	return false;
}

static bool validate_orbital_canonical_layout(const md_vlx_orbital_t* orb, size_t num_ao, const char* label) {
	ASSERT(orb);
	ASSERT(label);

	if (orb->coefficients.data) {
		const size_t num_mo = orb->coefficients.size[0];
		if (num_mo == 0 || orb->coefficients.size[1] != num_ao) {
			MD_LOG_ERROR("%s coefficient matrix is not in canonical [MO x AO] layout, got [%zu x %zu], expected [num_mo x %zu]", label, orb->coefficients.size[0], orb->coefficients.size[1], num_ao);
			return false;
		}
		if (orb->energy.data && orb->energy.size != num_mo) {
			MD_LOG_ERROR("%s energy vector length mismatch, expected %zu, got %zu", label, num_mo, orb->energy.size);
			return false;
		}
		if (orb->occupancy.data && orb->occupancy.size != num_mo) {
			MD_LOG_ERROR("%s occupancy vector length mismatch, expected %zu, got %zu", label, num_mo, orb->occupancy.size);
			return false;
		}
	}

	if (orb->density.data) {
		if (!validate_square_matrix_dims(&orb->density, label)) {
			return false;
		}
		if (orb->density.size[0] != num_ao) {
			MD_LOG_ERROR("%s density matrix dimension mismatch, expected %zu, got %zu", label, num_ao, orb->density.size[0]);
			return false;
		}
	}

	return true;
}

static bool normalize_orbital_coefficients(md_vlx_orbital_t* orb, size_t num_ao, const int* remap, const char* label) {
	ASSERT(orb);
	ASSERT(label);

	if (!orb->coefficients.data) {
		return true;
	}

	const size_t rows = orb->coefficients.size[0];
	const size_t cols = orb->coefficients.size[1];
	if (rows == 0 || cols == 0) {
		MD_LOG_ERROR("%s coefficient matrix has invalid dimensions [%zu x %zu]", label, rows, cols);
		return false;
	}

	if (rows == num_ao) {
		const size_t num_mo = cols;
		if (remap) {
			ao_permute(orb->coefficients.data, num_ao, num_mo, remap);
		}
		orb->coefficients.size[0] = num_mo;
		orb->coefficients.size[1] = num_ao;
		return true;
	}

	if (cols == num_ao) {
		const size_t num_mo = rows;
		if (remap) {
			ao_permute_cols(orb->coefficients.data, num_mo, num_ao, remap);
		}
		return true;
	}

	MD_LOG_ERROR("%s coefficient matrix [%zu x %zu] does not contain AO dimension %zu", label, rows, cols, num_ao);
	return false;
}

static bool validate_scf_canonical_layout(const md_vlx_t* vlx) {
	ASSERT(vlx);

	if (vlx->scf.S.data) {
		if (!validate_square_matrix_dims(&vlx->scf.S, "SCF overlap")) {
			return false;
		}
	}

	if (vlx->scf.alpha.density.data) {
		if (!validate_square_matrix_dims(&vlx->scf.alpha.density, "Alpha density")) {
			return false;
		}
		if (!validate_orbital_canonical_layout(&vlx->scf.alpha, vlx->scf.alpha.density.size[0], "Alpha orbital")) {
			return false;
		}
	}

	if (vlx->scf.beta.density.data) {
		if (!validate_square_matrix_dims(&vlx->scf.beta.density, "Beta density")) {
			return false;
		}
		if (!validate_orbital_canonical_layout(&vlx->scf.beta, vlx->scf.beta.density.size[0], "Beta orbital")) {
			return false;
		}
	}

	if (vlx->scf.alpha.density.data && vlx->scf.S.data && vlx->scf.alpha.density.size[0] != vlx->scf.S.size[0]) {
		MD_LOG_ERROR("SCF overlap/AO dimension mismatch, alpha density is %zu and overlap is %zu", vlx->scf.alpha.density.size[0], vlx->scf.S.size[0]);
		return false;
	}
	if (vlx->scf.alpha.density.data && vlx->scf.beta.density.data && vlx->scf.alpha.density.size[0] != vlx->scf.beta.density.size[0]) {
		MD_LOG_ERROR("SCF alpha/beta AO dimension mismatch, alpha density is %zu and beta density is %zu", vlx->scf.alpha.density.size[0], vlx->scf.beta.density.size[0]);
		return false;
	}

	return true;
}


// Data extraction procedures
static bool h5_read_scf_data(md_vlx_t* vlx, hid_t handle) {
	char scf_type[64] = {0};
	if (!h5_read_cstr(scf_type, sizeof(scf_type), handle, "scf_type")) {
		return false;
	}

	if (str_eq_cstr(STR_LIT("restricted"), scf_type)) {
		vlx->scf.type = MD_VLX_SCF_RESTRICTED;
	} else if (str_eq_cstr(STR_LIT("restricted_openshell"), scf_type)) {
		vlx->scf.type = MD_VLX_SCF_RESTRICTED_OPENSHELL;
	} else if (str_eq_cstr(STR_LIT("unrestricted"), scf_type)) {
		vlx->scf.type = MD_VLX_SCF_UNRESTRICTED;
	} else {
		vlx->scf.type = MD_VLX_SCF_UNKNOWN;
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
	if (!validate_square_matrix_dims(&(md_vlx_2d_data_t){ .size = {den_dim[0], den_dim[1]}, .data = NULL }, "Alpha density")) {
		return false;
	}

	const size_t num_ao = den_dim[0];
	size_t num_mo = 0;
	if (!infer_num_mo_from_coeff_dims(&num_mo, dim, num_ao, "Alpha coefficient")) {
		return false;
	}

	md_array_resize(vlx->scf.alpha.coefficients.data, dim[0] * dim[1], vlx->arena);
	MEMCPY(vlx->scf.alpha.coefficients.size, dim, sizeof(dim));

	md_array_resize(vlx->scf.alpha.energy.data, num_mo, vlx->arena);
	vlx->scf.alpha.energy.size = num_mo;

	md_array_resize(vlx->scf.alpha.occupancy.data, num_mo, vlx->arena);
	vlx->scf.alpha.occupancy.size = num_mo;

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

	if (vlx->scf.type == MD_VLX_SCF_UNRESTRICTED) {
		size_t beta_dim[2];
		h5_read_dataset_dims(beta_dim, 2, handle, "C_beta");
		size_t beta_num_mo = 0;
		if (!infer_num_mo_from_coeff_dims(&beta_num_mo, beta_dim, num_ao, "Beta coefficient")) {
			return false;
		}

		size_t beta_den_dim[2];
		h5_read_dataset_dims(beta_den_dim, 2, handle, "D_beta");
		if (!validate_square_matrix_dims(&(md_vlx_2d_data_t){ .size = {beta_den_dim[0], beta_den_dim[1]}, .data = NULL }, "Beta density")) {
			return false;
		}
		if (beta_den_dim[0] != num_ao) {
			MD_LOG_ERROR("Alpha/Beta AO dimension mismatch, alpha density is %zu and beta density is %zu", num_ao, beta_den_dim[0]);
			return false;
		}

		md_array_resize(vlx->scf.beta.coefficients.data, beta_dim[0] * beta_dim[1], vlx->arena);
		MEMCPY(vlx->scf.beta.coefficients.size, beta_dim, sizeof(beta_dim));

		md_array_resize(vlx->scf.beta.energy.data, beta_num_mo, vlx->arena);
		vlx->scf.beta.energy.size = beta_num_mo;

		md_array_resize(vlx->scf.beta.occupancy.data, beta_num_mo, vlx->arena);
		vlx->scf.beta.occupancy.size = beta_num_mo;

		md_array_resize(vlx->scf.beta.density.data, beta_den_dim[0] * beta_den_dim[1], vlx->arena);
			MEMCPY(vlx->scf.beta.density.size, beta_den_dim, sizeof(beta_den_dim));

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
		if (vlx->scf.type == MD_VLX_SCF_RESTRICTED_OPENSHELL) {
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

	// NOTE: H5Lexists returns htri_t -- negative on error, which is truthy. Test
	// explicitly, or a failed lookup is taken as "present" and the read proceeds.
	if (h5_link_exists(handle, "scf_history")) {
		// Extract new history format. This contains individual groups for each iteration labeled '0' ... 'N'.
		// Each group contains scalar datasets for the values of that iteration.
		// The individual datasets are named:
		// - 'diff_density'
		// - 'diff_energy'
		// - 'energy'
		// - 'gradient_norm'
		// - 'max_gradient'

		hid_t scf_history = H5Gopen(handle, "scf_history", H5P_DEFAULT);
		if (scf_history < 0) {
			return false;
		}

		// Determine number of iterations by counting number of groups. Assume that all groups are iterations.
		hsize_t h5_num_links = 0;
		if (H5Gget_num_objs(scf_history, &h5_num_links) < 0) {
			return false;
		}
		size_t num_links = (size_t)h5_num_links;

		md_array_resize(vlx->scf.history.density_diff,	num_links, vlx->arena);
		md_array_resize(vlx->scf.history.energy_diff,	num_links, vlx->arena);
		md_array_resize(vlx->scf.history.energy,		num_links, vlx->arena);
		md_array_resize(vlx->scf.history.gradient_norm, num_links, vlx->arena);
		md_array_resize(vlx->scf.history.max_gradient,  num_links, vlx->arena);

		size_t num_iter = 0;
		for (size_t i = 0; i < num_links; ++i) {
			// We check and read them in iteration order
			char name_buf[64];
			snprintf(name_buf, sizeof(name_buf), "%zu", i);
			if (H5Lexists(scf_history, name_buf, H5P_DEFAULT) > 0) {
				hid_t iter_group = H5Gopen(scf_history, name_buf, H5P_DEFAULT);
				if (iter_group < 0) {
					return false;
				}
				if (!h5_read_dataset_data(&vlx->scf.history.density_diff[i], 1, iter_group, H5T_NATIVE_DOUBLE, "diff_density")) {
					return false;
				}
				if (!h5_read_dataset_data(&vlx->scf.history.energy_diff[i], 1, iter_group, H5T_NATIVE_DOUBLE, "diff_energy")) {
					return false;
				}
				if (!h5_read_dataset_data(&vlx->scf.history.energy[i], 1, iter_group, H5T_NATIVE_DOUBLE, "energy")) {
					return false;
				}
				if (!h5_read_dataset_data(&vlx->scf.history.gradient_norm[i], 1, iter_group, H5T_NATIVE_DOUBLE, "gradient_norm")) {
					return false;
				}
				if (!h5_read_dataset_data(&vlx->scf.history.max_gradient[i], 1, iter_group, H5T_NATIVE_DOUBLE, "max_gradient")) {
					return false;
				}
				H5Gclose(iter_group);
				num_iter += 1;
			}
		}
		vlx->scf.history.number_of_iterations = num_iter;
	} else if (h5_link_exists(handle, "scf_history_energy")) {
		size_t scf_hist_len = 0;
		if (!h5_read_dataset_dims(&scf_hist_len, 1, handle, "scf_history_energy")) {
			return false;
		}

		if (scf_hist_len > 0) {
			vlx->scf.history.number_of_iterations = scf_hist_len;
			md_array_resize(vlx->scf.history.density_diff, scf_hist_len, vlx->arena);
			md_array_resize(vlx->scf.history.energy, scf_hist_len, vlx->arena);
			md_array_resize(vlx->scf.history.energy_diff, scf_hist_len, vlx->arena);
			md_array_resize(vlx->scf.history.gradient_norm, scf_hist_len, vlx->arena);
			md_array_resize(vlx->scf.history.max_gradient, scf_hist_len, vlx->arena);
		}

		if (!h5_read_dataset_data(vlx->scf.history.density_diff, scf_hist_len, handle, H5T_NATIVE_DOUBLE, "scf_history_diff_density")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.history.energy_diff, scf_hist_len, handle, H5T_NATIVE_DOUBLE, "scf_history_diff_energy")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.history.energy, scf_hist_len, handle, H5T_NATIVE_DOUBLE, "scf_history_energy")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.history.gradient_norm, scf_hist_len, handle, H5T_NATIVE_DOUBLE, "scf_history_gradient_norm")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.history.max_gradient, scf_hist_len, handle, H5T_NATIVE_DOUBLE, "scf_history_max_gradient")) {
			return false;
		}
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

static bool h5_read_optional_1d_data(md_vlx_1d_data_t* out_data, hid_t handle, const char* field_name, md_allocator_i* arena) {
	ASSERT(out_data);
	ASSERT(field_name);
	ASSERT(arena);

	if (!h5_check_dataset_exists(handle, field_name)) {
		return true;
	}

	size_t dim[2] = {0};
	int num_dims = h5_read_dataset_dims(dim, 2, handle, field_name);
	if (num_dims <= 0) {
		MD_LOG_ERROR("Invalid dimensions in H5 vector dataset '%s'", field_name);
		return false;
	}

	size_t sample_count = 1;
	for (int dim_idx = 0; dim_idx < num_dims; ++dim_idx) {
		sample_count *= dim[dim_idx];
	}

	if (sample_count == 0) {
		MD_LOG_ERROR("Empty H5 vector dataset '%s'", field_name);
		return false;
	}

	md_array_resize(out_data->data, sample_count, arena);
	MEMSET(out_data->data, 0, md_array_bytes(out_data->data));
	out_data->size = sample_count;

	if (!h5_read_dataset_data(out_data->data, sample_count, handle, H5T_NATIVE_DOUBLE, field_name)) {
		return false;
	}

	return true;
}

// The occupied/virtual split the response solution vectors are indexed by.
//
// VeloxChem writes num_core/num_valence/num_virtual only for some calculations - none of the files
// in test_data carries them - and without that split a solution vector is an undifferentiated run
// of amplitudes that nothing can be reconstructed from. So when the file is silent, derive it: an
// ordinary valence excitation spans every occupied orbital and every virtual one, which the SCF
// occupations already say.
//
// Derived, never assumed: the split has to reproduce the solution vector's own length (amp_count,
// or twice it when the vector carries both X and Y), and it is adopted only when it does. That
// check is what makes this safe for a core excitation, where the file DOES name num_core and the
// derived valence split would be wrong - there, the stored values are kept and this does nothing.
static void vlx_rsp_infer_occupied_virtual_split(md_vlx_t* vlx) {
	ASSERT(vlx);

	if (vlx->rsp.num_core > 0 || vlx->rsp.num_valence > 0) {
		return;     // the file said so
	}
	const size_t vec_len = vlx->rsp.solution_matrix.size[1];
	if (!vlx->rsp.solution_matrix.data || vec_len == 0) {
		return;
	}

	const double* occ    = vlx->scf.alpha.occupancy.data;
	const size_t  num_mo = vlx->scf.alpha.occupancy.size;
	if (!occ || num_mo == 0) {
		MD_LOG_ERROR("Response data has no occupied/virtual split and no SCF occupations to derive one from");
		return;
	}

	size_t nocc = 0;
	for (size_t i = 0; i < num_mo; ++i) {
		if (occ[i] > 0.0) nocc += 1;
	}
	const size_t nvir = num_mo - nocc;

	if (nocc == 0 || nvir == 0) {
		MD_LOG_ERROR("Cannot derive an occupied/virtual split: %zu of %zu orbitals are occupied", nocc, num_mo);
		return;
	}

	const size_t amp_count = nocc * nvir;
	if (vec_len != amp_count && vec_len != 2 * amp_count) {
		MD_LOG_ERROR("Derived %zu occupied x %zu virtual orbitals, which does not explain a solution vector of %zu values (expected %zu or %zu)",
			nocc, nvir, vec_len, amp_count, 2 * amp_count);
		return;
	}

	vlx->rsp.num_valence = nocc;
	vlx->rsp.num_virtual = nvir;
	MD_LOG_DEBUG("Derived the response occupied/virtual split from the SCF occupations: %zu x %zu", nocc, nvir);
}

static bool h5_read_rsp_data(md_vlx_t* vlx, hid_t handle) {

	h5_read_scalar(&vlx->rsp.number_of_frequencies, handle, H5T_NATIVE_HSIZE, "number_of_states");
	if (vlx->rsp.number_of_frequencies > 0) {
		// Standard Linear Response data, allocate and read
		vlx->rsp.type = MD_VLX_RSP_LINEAR;

		if (h5_check_dataset_exists(handle, "num_core")) {
			h5_read_scalar(&vlx->rsp.num_core, handle, H5T_NATIVE_INT64, "num_core");
		}

		if (h5_check_dataset_exists(handle, "num_valence")) {
			h5_read_scalar(&vlx->rsp.num_valence, handle, H5T_NATIVE_INT64, "num_valence");
		}

		if (h5_check_dataset_exists(handle, "num_virtual")) {
			h5_read_scalar(&vlx->rsp.num_virtual, handle, H5T_NATIVE_INT64, "num_virtual");
		}

        if (h5_check_dataset_exists(handle, "eigenvalues")) {
			md_array_resize(vlx->rsp.frequencies,	vlx->rsp.number_of_frequencies, vlx->arena);
			MEMSET(vlx->rsp.frequencies, 0, vlx->rsp.number_of_frequencies * sizeof(double));
			if (!h5_read_dataset_data(vlx->rsp.frequencies, md_array_size(vlx->rsp.frequencies), handle, H5T_NATIVE_DOUBLE, "eigenvalues")) {
				return false;
			}
        }

		// Response eigenvectors used to derive NTOs and attachment/detachment matrices.
		if (h5_check_dataset_exists(handle, "S1")) {
            size_t dims[2] = { vlx->rsp.number_of_frequencies, 0 };
            h5_read_dataset_dims(&dims[1], 1, handle, "S1");

            size_t len = dims[0] * dims[1];
			if (len > 0) {
				md_array_resize(vlx->rsp.solution_matrix.data, len, vlx->arena);
                vlx->rsp.solution_matrix.size[0] = dims[0];
                vlx->rsp.solution_matrix.size[1] = dims[1];
				char field_name[16];
				for (size_t state_idx = 0; state_idx < vlx->rsp.number_of_frequencies; ++state_idx) {
					snprintf(field_name, sizeof(field_name), "S%zu", state_idx + 1);
					double* dst = vlx->rsp.solution_matrix.data + state_idx * vlx->rsp.solution_matrix.size[1];
					if (!h5_read_dataset_data(dst, vlx->rsp.solution_matrix.size[1], handle, H5T_NATIVE_DOUBLE, field_name)) {
						return false;
					}
				}
			}
		} else if (h5_check_dataset_exists(handle, "full_solutions_matrix")) {
            size_t dims[2] = { 0 };
            h5_read_dataset_dims(dims, 2, handle, "full_solutions_matrix");
            size_t len = dims[0] * dims[1];
			if (len > 0) {
				md_array_resize(vlx->rsp.solution_matrix.data, len, vlx->arena);
                vlx->rsp.solution_matrix.size[0] = dims[0];
                vlx->rsp.solution_matrix.size[1] = dims[1];
                if (!h5_read_dataset_data(vlx->rsp.solution_matrix.data, len, handle, H5T_NATIVE_DOUBLE, "full_solutions_matrix")) {
					return false;
				}
			}
		}
	}

	// After the solution matrix is in hand, so its length is available to check a derived split
	// against, and after the SCF block was read - which vlx_read_h5_file guarantees.
	vlx_rsp_infer_occupied_virtual_split(vlx);

	if (vlx->rsp.type == MD_VLX_RSP_UNKNOWN) {
		// No standard response data, check for other types of response data by looking for type field
		if (h5_check_dataset_exists(handle, "rsp_type")) {
			char type_buf[32] = { 0 };
			h5_read_cstr(type_buf, sizeof(type_buf), handle, "rsp_type");
			if (strncmp(type_buf, "cpp", sizeof(type_buf)) == 0) {
				vlx->rsp.type = MD_VLX_RSP_CPP;
			} else if (strncmp(type_buf, "c6", sizeof(type_buf)) == 0) {
				vlx->rsp.type = MD_VLX_RSP_C6;
			} else if (strncmp(type_buf, "tpa_transition", sizeof(type_buf)) == 0) {
				vlx->rsp.type = MD_VLX_RSP_TPA_TRANSITION;
			} else if (strncmp(type_buf, "tpa", sizeof(type_buf)) == 0) {
				vlx->rsp.type = MD_VLX_RSP_TPA;
			} else if (strncmp(type_buf, "rixs", sizeof(type_buf)) == 0) {
                vlx->rsp.type = MD_VLX_RSP_RIXS;
			}
		}

		if (vlx->rsp.type == MD_VLX_RSP_C6) {
			// Try to read the c6 field
			if (h5_check_dataset_exists(handle, "c6")) {
                if (!h5_read_dataset_data(&vlx->rsp.c6, 1, handle, H5T_NATIVE_DOUBLE, "c6")) {
                    MD_LOG_ERROR("Could not read c6 dataset");
                    return false;
                }
			}
		} else if (vlx->rsp.type == MD_VLX_RSP_CPP || vlx->rsp.type == MD_VLX_RSP_TPA) {
			size_t dim;
			if (h5_read_dataset_dims(&dim, 1, handle, "frequencies")) {
				vlx->rsp.number_of_frequencies = dim;
				md_array_resize(vlx->rsp.frequencies, dim, vlx->arena);
				if (!h5_read_dataset_data(vlx->rsp.frequencies, md_array_size(vlx->rsp.frequencies), handle, H5T_NATIVE_DOUBLE, "frequencies")) {
					return false;
				}
			}
		} else if (vlx->rsp.type == MD_VLX_RSP_TPA_TRANSITION) {
			size_t dim;
			if (h5_read_dataset_dims(&dim, 1, handle, "photon_energies")) {
				vlx->rsp.number_of_frequencies = dim;
				md_array_resize(vlx->rsp.frequencies, dim, vlx->arena);
				if (!h5_read_dataset_data(vlx->rsp.frequencies, md_array_size(vlx->rsp.frequencies), handle, H5T_NATIVE_DOUBLE, "photon_energies")) {
					return false;
				}
			}
		} else if (vlx->rsp.type == MD_VLX_RSP_RIXS) {
			size_t dim;
			if (h5_read_dataset_dims(&dim, 1, handle, "core_eigenvalues")) {
				vlx->rsp.number_of_frequencies = dim;
				md_array_resize(vlx->rsp.frequencies, dim, vlx->arena);
				if (!h5_read_dataset_data(vlx->rsp.frequencies, md_array_size(vlx->rsp.frequencies), handle, H5T_NATIVE_DOUBLE, "core_eigenvalues")) {
					return false;
				}
			}

			if (h5_read_dataset_dims(&dim, 1, handle, "photon_energies")) {
				vlx->rsp.rixs.num_incomming_photons = dim;
				md_array_resize(vlx->rsp.rixs.photon_energies, dim, vlx->arena);
				if (!h5_read_dataset_data(vlx->rsp.rixs.photon_energies, md_array_size(vlx->rsp.rixs.photon_energies), handle, H5T_NATIVE_DOUBLE, "photon_energies")) {
					return false;
				}
			}
		}
	}

	if (vlx->rsp.number_of_frequencies > 0) {
		// Dipoles
		size_t num_dipole_points = vlx->rsp.number_of_frequencies * 3;
		if (h5_check_dataset_exists(handle, "electric_transition_dipoles")) {
			md_array_resize(vlx->rsp.electric_transition_dipoles, vlx->rsp.number_of_frequencies, vlx->arena);
			MEMSET(vlx->rsp.electric_transition_dipoles, 0, md_array_bytes(vlx->rsp.electric_transition_dipoles));
			if (!h5_read_dataset_data(vlx->rsp.electric_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "electric_transition_dipoles")) {
				md_array_free(vlx->rsp.electric_transition_dipoles, vlx->arena);
				vlx->rsp.electric_transition_dipoles = NULL;
			}
		}
		
		if (h5_check_dataset_exists(handle, "magnetic_transition_dipoles")) {
			md_array_resize(vlx->rsp.magnetic_transition_dipoles, vlx->rsp.number_of_frequencies, vlx->arena);
			MEMSET(vlx->rsp.magnetic_transition_dipoles, 0, md_array_bytes(vlx->rsp.magnetic_transition_dipoles));
			if (!h5_read_dataset_data(vlx->rsp.magnetic_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "magnetic_transition_dipoles")) {
				md_array_free(vlx->rsp.magnetic_transition_dipoles, vlx->arena);
				vlx->rsp.magnetic_transition_dipoles = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "velocity_transition_dipoles")) {
			md_array_resize(vlx->rsp.velocity_transition_dipoles, vlx->rsp.number_of_frequencies, vlx->arena);
			MEMSET(vlx->rsp.velocity_transition_dipoles, 0, md_array_bytes(vlx->rsp.velocity_transition_dipoles));
			if (!h5_read_dataset_data(vlx->rsp.velocity_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "velocity_transition_dipoles")) {
				md_array_free(vlx->rsp.velocity_transition_dipoles, vlx->arena);
				vlx->rsp.velocity_transition_dipoles = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "tpa_strengths/circular")) {
            md_array_resize(vlx->rsp.tpa_strengths_circular, vlx->rsp.number_of_frequencies, vlx->arena);
            MEMSET(vlx->rsp.tpa_strengths_circular, 0, md_array_bytes(vlx->rsp.tpa_strengths_circular));
            if (!h5_read_dataset_data(vlx->rsp.tpa_strengths_circular, vlx->rsp.number_of_frequencies, handle, H5T_NATIVE_DOUBLE, "tpa_strengths/circular")) {
                md_array_free(vlx->rsp.tpa_strengths_circular, vlx->arena);
                vlx->rsp.tpa_strengths_circular = NULL;
            }
		}

        if (h5_check_dataset_exists(handle, "tpa_strengths/linear")) {
            md_array_resize(vlx->rsp.tpa_strengths_linear, vlx->rsp.number_of_frequencies, vlx->arena);
            MEMSET(vlx->rsp.tpa_strengths_linear, 0, md_array_bytes(vlx->rsp.tpa_strengths_linear));
            if (!h5_read_dataset_data(vlx->rsp.tpa_strengths_linear, vlx->rsp.number_of_frequencies, handle, H5T_NATIVE_DOUBLE, "tpa_strengths/linear")) {
                md_array_free(vlx->rsp.tpa_strengths_linear, vlx->arena);
                vlx->rsp.tpa_strengths_linear = NULL;
            }
        }

		if (vlx->rsp.type == MD_VLX_RSP_RIXS) {
			// RIXS involves three independent dimensions:
			//   C = number of core-excited (intermediate) states  -> rsp.number_of_frequencies
			//   F = number of final (valence-excited) states       -> rixs.num_final_states
			//   P = number of incoming photon energies            -> rixs.num_incomming_photons
			// The core states are summed over coherently inside the scattering amplitude and never
			// appear as an output dimension, so F is completely unrelated to C. The 2D datasets are
			// stored row-major as [F][P]. Derive F from a representative dataset rather than assuming.
			static const char* rixs_2d_fields[] = { "cross_sections", "energy_losses", "emission_energies" };
			for (size_t i = 0; i < ARRAY_SIZE(rixs_2d_fields); ++i) {
				size_t dim[2] = { 0 };
				if (h5_read_dataset_dims(dim, (int)ARRAY_SIZE(dim), handle, rixs_2d_fields[i]) == 2 && dim[0] > 0 && dim[1] > 0) {
					vlx->rsp.rixs.num_final_states = dim[0];

					if (vlx->rsp.rixs.num_incomming_photons == 0) {
						// 'photon_energies' was missing or unreadable, recover P from the column count.
						vlx->rsp.rixs.num_incomming_photons = dim[1];
					} else if (dim[1] != vlx->rsp.rixs.num_incomming_photons) {
						MD_LOG_ERROR("RIXS: dataset '%s' has %i columns, expected %i incoming photon energies",
							rixs_2d_fields[i], (int)dim[1], (int)vlx->rsp.rixs.num_incomming_photons);
						vlx->rsp.rixs.num_final_states = 0;
					}
					break;
				}
			}

			// Element count shared by all the [F][P] datasets below.
			const size_t num_2d_elem = vlx->rsp.rixs.num_final_states * vlx->rsp.rixs.num_incomming_photons;

			if (num_2d_elem > 0 && h5_check_dataset_exists(handle, "cross_sections")) {
				md_array_resize(vlx->rsp.rixs.cross_sections, num_2d_elem, vlx->arena);
				MEMSET(vlx->rsp.rixs.cross_sections, 0, md_array_bytes(vlx->rsp.rixs.cross_sections));
				if (!h5_read_dataset_data(vlx->rsp.rixs.cross_sections, md_array_size(vlx->rsp.rixs.cross_sections), handle, H5T_NATIVE_DOUBLE, "cross_sections")) {
					md_array_free(vlx->rsp.rixs.cross_sections, vlx->arena);
					vlx->rsp.rixs.cross_sections = NULL;
				}
			}

			if (h5_check_dataset_exists(handle, "core_osc_strengths")) {
				md_array_resize(vlx->rsp.rixs.core_osc_strengths, vlx->rsp.number_of_frequencies, vlx->arena);
				MEMSET(vlx->rsp.rixs.core_osc_strengths, 0, md_array_bytes(vlx->rsp.rixs.core_osc_strengths));
				if (!h5_read_dataset_data(vlx->rsp.rixs.core_osc_strengths, md_array_size(vlx->rsp.rixs.core_osc_strengths), handle, H5T_NATIVE_DOUBLE, "core_osc_strengths")) {
					md_array_free(vlx->rsp.rixs.core_osc_strengths, vlx->arena);
					vlx->rsp.rixs.core_osc_strengths = NULL;
				}
			}

			if (h5_check_dataset_exists(handle, "elastic_cross_sections")) {
				md_array_resize(vlx->rsp.rixs.elastic_cross_sections, vlx->rsp.rixs.num_incomming_photons, vlx->arena);
				MEMSET(vlx->rsp.rixs.elastic_cross_sections, 0, md_array_bytes(vlx->rsp.rixs.elastic_cross_sections));
				if (!h5_read_dataset_data(vlx->rsp.rixs.elastic_cross_sections, md_array_size(vlx->rsp.rixs.elastic_cross_sections), handle, H5T_NATIVE_DOUBLE, "elastic_cross_sections")) {
					md_array_free(vlx->rsp.rixs.elastic_cross_sections, vlx->arena);
					vlx->rsp.rixs.elastic_cross_sections = NULL;
				}
			}

			if (num_2d_elem > 0 && h5_check_dataset_exists(handle, "emission_energies")) {
				md_array_resize(vlx->rsp.rixs.emission_energies, num_2d_elem, vlx->arena);
				MEMSET(vlx->rsp.rixs.emission_energies, 0, md_array_bytes(vlx->rsp.rixs.emission_energies));
				if (!h5_read_dataset_data(vlx->rsp.rixs.emission_energies, md_array_size(vlx->rsp.rixs.emission_energies), handle, H5T_NATIVE_DOUBLE, "emission_energies")) {
					md_array_free(vlx->rsp.rixs.emission_energies, vlx->arena);
					vlx->rsp.rixs.emission_energies = NULL;
				}
			}

			if (num_2d_elem > 0 && h5_check_dataset_exists(handle, "energy_losses")) {
				md_array_resize(vlx->rsp.rixs.energy_losses, num_2d_elem, vlx->arena);
				MEMSET(vlx->rsp.rixs.energy_losses, 0, md_array_bytes(vlx->rsp.rixs.energy_losses));
				if (!h5_read_dataset_data(vlx->rsp.rixs.energy_losses, md_array_size(vlx->rsp.rixs.energy_losses), handle, H5T_NATIVE_DOUBLE, "energy_losses")) {
					md_array_free(vlx->rsp.rixs.energy_losses, vlx->arena);
					vlx->rsp.rixs.energy_losses = NULL;
				}
			}

			if (h5_check_dataset_exists(handle, "photon_energies")) {
				md_array_resize(vlx->rsp.rixs.photon_energies, vlx->rsp.rixs.num_incomming_photons, vlx->arena);
				MEMSET(vlx->rsp.rixs.photon_energies, 0, md_array_bytes(vlx->rsp.rixs.photon_energies));
				if (!h5_read_dataset_data(vlx->rsp.rixs.photon_energies, md_array_size(vlx->rsp.rixs.photon_energies), handle, H5T_NATIVE_DOUBLE, "photon_energies")) {
					md_array_free(vlx->rsp.rixs.photon_energies, vlx->arena);
					vlx->rsp.rixs.photon_energies = NULL;
				}
			}

			if (h5_check_dataset_exists(handle, "gamma_fwhm_ev")) {
				h5_read_scalar(&vlx->rsp.rixs.gamma_fwhm_ev, handle, H5T_NATIVE_DOUBLE, "gamma_fwhm_ev");
			}

			if (h5_check_dataset_exists(handle, "scattering_amplitudes")) {
				// The scattering amplitudes are complex and have the shape [F][P][3][3], where F is
				// the number of final states, P the number of incoming photon energies and the
				// trailing 3x3 is the Cartesian scattering amplitude tensor.
				// Derive the element count from the dataset itself rather than assuming a rank.
				size_t dim[4] = {0};
				int ndim = h5_read_dataset_dims(dim, (int)ARRAY_SIZE(dim), handle, "scattering_amplitudes");

				// h5_read_dataset_dims reports the rank of the dataset, which may exceed the number of
				// entries it actually wrote, so clamp before iterating.
				if (ndim > (int)ARRAY_SIZE(dim)) {
					MD_LOG_ERROR("Unexpected rank of H5 dataset 'scattering_amplitudes'");
					ndim = 0;
				}

				size_t num_elem = (ndim > 0) ? 1 : 0;
				for (int i = 0; i < ndim; ++i) {
					num_elem *= dim[i];
				}

				// Sanity check against the dimensions derived above: [F][P][3][3] == F * P * 9.
				if (num_elem > 0 && num_2d_elem > 0 && num_elem != num_2d_elem * 9) {
					MD_LOG_ERROR("RIXS: 'scattering_amplitudes' holds %i elements, expected %i (%i final states x %i photon energies x 3 x 3)",
						(int)num_elem, (int)(num_2d_elem * 9), (int)vlx->rsp.rixs.num_final_states, (int)vlx->rsp.rixs.num_incomming_photons);
				}

				if (num_elem > 0) {
					md_array_resize(vlx->rsp.rixs.scattering_amplitude_re, num_elem, vlx->arena);
					md_array_resize(vlx->rsp.rixs.scattering_amplitude_im, num_elem, vlx->arena);

					MEMSET(vlx->rsp.rixs.scattering_amplitude_re, 0, md_array_bytes(vlx->rsp.rixs.scattering_amplitude_re));
					MEMSET(vlx->rsp.rixs.scattering_amplitude_im, 0, md_array_bytes(vlx->rsp.rixs.scattering_amplitude_im));

					// Splits the interleaved complex data on disk into two separate arrays.
					if (!h5_read_complex_dataset_split(vlx->rsp.rixs.scattering_amplitude_re, vlx->rsp.rixs.scattering_amplitude_im, num_elem, handle, "scattering_amplitudes")) {
						md_array_free(vlx->rsp.rixs.scattering_amplitude_re, vlx->arena);
						md_array_free(vlx->rsp.rixs.scattering_amplitude_im, vlx->arena);
						vlx->rsp.rixs.scattering_amplitude_re = NULL;
						vlx->rsp.rixs.scattering_amplitude_im = NULL;
					}
				}
			}
		} else {
			if (h5_check_dataset_exists(handle, "cross_sections")) {
				if (vlx->rsp.type != MD_VLX_RSP_RIXS) {
					md_array_resize(vlx->rsp.cross_sections, vlx->rsp.number_of_frequencies, vlx->arena);
					MEMSET(vlx->rsp.cross_sections, 0, md_array_bytes(vlx->rsp.cross_sections));
					if (!h5_read_dataset_data(vlx->rsp.cross_sections, vlx->rsp.number_of_frequencies, handle, H5T_NATIVE_DOUBLE, "cross_sections")) {
						md_array_free(vlx->rsp.cross_sections, vlx->arena);
						vlx->rsp.cross_sections = NULL;
					}
				}
			}
		}

		if (h5_check_dataset_exists(handle, "sigma")) {
			md_array_resize(vlx->rsp.sigmas, vlx->rsp.number_of_frequencies, vlx->arena);
            MEMSET(vlx->rsp.sigmas, 0, md_array_bytes(vlx->rsp.sigmas));
			if (!h5_read_dataset_data(vlx->rsp.sigmas, md_array_size(vlx->rsp.sigmas), handle, H5T_NATIVE_DOUBLE, "sigma")) {
                md_array_free(vlx->rsp.sigmas, vlx->arena);
                vlx->rsp.sigmas = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "optical-rotation")) {
			md_array_resize(vlx->rsp.optical_rotations, vlx->rsp.number_of_frequencies, vlx->arena);
            MEMSET(vlx->rsp.optical_rotations, 0, md_array_bytes(vlx->rsp.optical_rotations));
			if (!h5_read_dataset_data(vlx->rsp.optical_rotations, md_array_size(vlx->rsp.optical_rotations), handle, H5T_NATIVE_DOUBLE, "optical-rotation")) {
                md_array_free(vlx->rsp.optical_rotations, vlx->arena);
                vlx->rsp.optical_rotations = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "delta-epsilon")) {
			md_array_resize(vlx->rsp.delta_epsilons, vlx->rsp.number_of_frequencies, vlx->arena);
            MEMSET(vlx->rsp.delta_epsilons, 0, md_array_bytes(vlx->rsp.delta_epsilons));
			if (!h5_read_dataset_data(vlx->rsp.delta_epsilons, md_array_size(vlx->rsp.delta_epsilons), handle, H5T_NATIVE_DOUBLE, "delta-epsilon")) {
                md_array_free(vlx->rsp.delta_epsilons, vlx->arena);
                vlx->rsp.delta_epsilons = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "oscillator_strengths")) {
			md_array_resize(vlx->rsp.oscillator_strengths, vlx->rsp.number_of_frequencies, vlx->arena);
			MEMSET(vlx->rsp.oscillator_strengths, 0, md_array_bytes(vlx->rsp.oscillator_strengths));
			if (!h5_read_dataset_data(vlx->rsp.oscillator_strengths, md_array_size(vlx->rsp.oscillator_strengths), handle, H5T_NATIVE_DOUBLE, "oscillator_strengths")) {
                md_array_free(vlx->rsp.oscillator_strengths, vlx->arena);
                vlx->rsp.oscillator_strengths = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "rotatory_strengths")) {
			md_array_resize(vlx->rsp.rotatory_strengths, vlx->rsp.number_of_frequencies, vlx->arena);
			MEMSET(vlx->rsp.rotatory_strengths, 0, md_array_bytes(vlx->rsp.rotatory_strengths));
			if (!h5_read_dataset_data(vlx->rsp.rotatory_strengths, md_array_size(vlx->rsp.rotatory_strengths), handle, H5T_NATIVE_DOUBLE, "rotatory_strengths")) {
                md_array_free(vlx->rsp.rotatory_strengths, vlx->arena);
                vlx->rsp.rotatory_strengths = NULL;
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

	if (h5_check_dataset_exists(handle, "force_constants")) {
		md_array_resize(vlx->vib.force_constants, number_of_modes, vlx->arena);
		if (!h5_read_dataset_data(vlx->vib.force_constants, md_array_size(vlx->vib.force_constants), handle, H5T_NATIVE_DOUBLE, "force_constants")) {
			return false;
		}
	}

	if (h5_check_dataset_exists(handle, "ir_intensities")) {
		md_array_resize(vlx->vib.ir_intensities, number_of_modes, vlx->arena);
		if (!h5_read_dataset_data(vlx->vib.ir_intensities, md_array_size(vlx->vib.ir_intensities), handle, H5T_NATIVE_DOUBLE, "ir_intensities")) {
			return false;
		}
	}

	if (h5_check_dataset_exists(handle, "vib_frequencies")) {
		md_array_resize(vlx->vib.frequencies, number_of_modes, vlx->arena);
		if (!h5_read_dataset_data(vlx->vib.frequencies, md_array_size(vlx->vib.frequencies), handle, H5T_NATIVE_DOUBLE, "vib_frequencies")) {
			return false;
		}
	}

	if (h5_check_dataset_exists(handle, "reduced_masses")) {
		md_array_resize(vlx->vib.reduced_masses, number_of_modes, vlx->arena);
		if (!h5_read_dataset_data(vlx->vib.reduced_masses, md_array_size(vlx->vib.reduced_masses), handle, H5T_NATIVE_DOUBLE, "reduced_masses")) {
			return false;
		}
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

	size_t number_of_external_frequencies = 0;
	if (h5_read_scalar(&number_of_external_frequencies, handle, H5T_NATIVE_HSIZE, "number_of_external_frequencies")) {
        vlx->vib.num_external_frequencies = number_of_external_frequencies;

        if (h5_check_dataset_exists(handle, "external_frequencies")) {
            md_array_resize(vlx->vib.external_frequencies, number_of_external_frequencies, vlx->arena);
            if (!h5_read_dataset_data(vlx->vib.external_frequencies, md_array_size(vlx->vib.external_frequencies), handle, H5T_NATIVE_DOUBLE, "external_frequencies")) {
                return false;
            }
        }

		if (h5_check_dataset_exists(handle, "raman_activities")) {
            md_array_resize(vlx->vib.raman_activities, number_of_external_frequencies * number_of_modes, vlx->arena);
            if (!h5_read_dataset_data(vlx->vib.raman_activities, md_array_size(vlx->vib.raman_activities), handle, H5T_NATIVE_DOUBLE, "raman_activities")) {
                return false;
            }
		}
	}

	return true;
}

// This procedure is an abstraction to help read illformed groups containing a collection of individual datasets containing one value each.
// The expected names of the datasets are '0', '1', ... up to the expected count. The dataset prefix is used for error messages only, it does not have to be present in the actual dataset names.
static bool h5_extract_group_as_array_double(md_array(double)* out_data, hid_t group_handle, md_allocator_i* arena) {
	ASSERT(out_data);
	ASSERT(arena);

	hsize_t num_links = 0;
	if (H5Gget_num_objs(group_handle, &num_links) < 0) {
		MD_LOG_ERROR("Failed to get number of links in group");
		return false;
	}

	size_t expected_count = (size_t)num_links;
	md_array_ensure(*out_data, expected_count, arena);

	for (size_t i = 0; i < expected_count; ++i) {
		char dataset_name[64];
		snprintf(dataset_name, sizeof(dataset_name), "%zu", i);
		if (!h5_check_dataset_exists(group_handle, dataset_name)) {
			MD_LOG_ERROR("Expected dataset '%s' in group not found", dataset_name);
			return false;
		}

		double value;
		if (!h5_read_dataset_data(&value, 1, group_handle, H5T_NATIVE_DOUBLE, dataset_name)) {
			MD_LOG_ERROR("Failed to read dataset '%s' in group", dataset_name);
			return false;
		}
		md_array_push_no_grow(*out_data, value);
	}
	return true;
}

static bool h5_extract_as_array_double(md_array(double)* out_data, hid_t handle, const char* name, md_allocator_i* arena) {
	ASSERT(out_data);
	ASSERT(name);
	ASSERT(arena);

	bool result = false;
	
	H5I_type_t obj_type = h5_get_object_type(handle, name);
	if (obj_type == H5I_GROUP) {
		hid_t obj_handle = H5Oopen(handle, name, H5P_DEFAULT);
		if (obj_handle < 0) {
			MD_LOG_ERROR("Failed to open object '%s'", name);
			return false;
		}
		result = h5_extract_group_as_array_double(out_data, obj_handle, arena);
		H5Oclose(obj_handle);
	} else if (obj_type == H5I_DATASET) {
		size_t dim;
		int num_dim = h5_read_dataset_dims(&dim, 1, handle, name);
		if (num_dim <= 0) {
			MD_LOG_ERROR("Invalid dimensions in dataset '%s'", name);
			goto done;
		}
		md_array_resize(*out_data, dim, arena);
		MEMSET(*out_data, 0, md_array_bytes(*out_data));
		result = h5_read_dataset_data(*out_data, dim, handle, H5T_NATIVE_DOUBLE, name);
	} else {
		MD_LOG_ERROR("Unrecognized object type for '%s'", name);
	}

done:
	return result;
}

static bool h5_read_opt_data(md_vlx_t* vlx, hid_t handle) {
	const md_vlx_opt_type_t opt_types[] = { MD_VLX_OPT_GEOMETRY, MD_VLX_OPT_CONSTRAINED, MD_VLX_OPT_IRC };
	const char* valid_prefixes[] = { "opt", "scan", "irc" };
	const char* energy_ident = NULL;
	const char* coord_ident = NULL;
	H5I_type_t  coord_type = H5I_BADID;
    md_vlx_opt_type_t opt_type = MD_VLX_OPT_UNKNOWN;

	for (size_t i = 0; i < ARRAY_SIZE(valid_prefixes); ++i) {
		H5I_type_t type = -1;
		char energy_name[32];
		char coord_name[32];

		snprintf(energy_name, sizeof(energy_name), "%s_energies", valid_prefixes[i]);
		snprintf(coord_name, sizeof(coord_name), "%s_coordinates_au", valid_prefixes[i]);

		type = h5_get_object_type(handle, coord_name);
		if (type == H5I_DATASET || type == H5I_GROUP) {
			coord_ident = coord_name;
			coord_type = type;
		}

		type = h5_get_object_type(handle, energy_name);
		if (type == H5I_DATASET || type == H5I_GROUP) {
			energy_ident = energy_name;
			opt_type = opt_types[i];
			break;
		}
	}

	if (energy_ident) {
		// Extract energies
		if (!h5_extract_as_array_double(&vlx->opt.energies, handle, energy_ident, vlx->arena)) {
			return false;
		}

		size_t len = md_array_size(vlx->opt.energies);
		vlx->opt.number_of_steps = len;

		if (h5_check_dataset_exists(handle, "state_index")) {
			if (!h5_read_scalar(&vlx->opt.state_index, handle, H5T_NATIVE_INT64, "state_index")) {
                return false;
            }
		}

        if (h5_check_dataset_exists(handle, "ts_index")) {
            if (!h5_read_scalar(&vlx->opt.ts_index, handle, H5T_NATIVE_INT64, "ts_index")) {
                return false;
            }
        }

		if (coord_ident && coord_type == H5I_DATASET) {
			// Extract coordinates
			size_t dim[4];
			int num_dim = h5_read_dataset_dims(dim, ARRAY_SIZE(dim), handle, coord_ident);
			if (num_dim <= 0) {
				MD_LOG_ERROR("Invalid dimensions in '%s'", coord_ident);
				return false;
			}

			if (dim[1] != vlx->number_of_atoms || dim[2] != 3) {
				MD_LOG_ERROR("Unexpected dimensions in '%s'", coord_ident);
				return false;
			}
		
			if (dim[0] != len) {
				MD_LOG_ERROR("Energy/coordinate step count mismatch between '%s' and '%s'", energy_ident, coord_ident);
				return false;
			}

			// Read coordinate data
			md_array_resize(vlx->opt.coordinates, dim[0] * dim[1], vlx->arena);
			if (!h5_read_dataset_data(vlx->opt.coordinates, md_array_size(vlx->opt.coordinates) * 3, handle, H5T_NATIVE_DOUBLE, coord_ident)) {
				return false;
			}

			if (vlx->opt.coordinates) {
				for (size_t i = 0; i < dim[0] * dim[1]; ++i) {
					vlx->opt.coordinates[i] = dvec3_mul1(vlx->opt.coordinates[i], BOHR_TO_ANGSTROM);
				}
			}
		}

        vlx->opt.type = opt_type;

		return true;
	}

	return false;
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

static const double* vlx_rsp_get_solution_vector(const md_vlx_t* vlx, size_t state_idx, size_t* out_nocc, size_t* out_nvir, size_t* out_amp_count, bool* out_has_y) {
	ASSERT(vlx);

	if (!vlx->rsp.solution_matrix.data || state_idx >= vlx->rsp.solution_matrix.size[0]) {
		return NULL;
	}

	const size_t nocc = vlx->rsp.num_core > 0 ? vlx->rsp.num_core : vlx->rsp.num_valence;
	if (nocc == 0) {
		return NULL;
	}

	const size_t nvir = vlx->rsp.num_virtual;
	if (nvir == 0) {
		return NULL;
	}

	size_t vec_size = vlx->rsp.solution_matrix.size[1];

	const size_t amp_count = nocc * nvir;
	bool has_y = false;
	if (vec_size == amp_count) {
		has_y = false;
	} else if (vec_size == 2 * amp_count) {
		has_y = true;
	} else {
		MD_LOG_ERROR("Unexpected response eigenvector length for state %zu: got %zu, expected %zu or %zu", state_idx + 1, vec_size, amp_count, 2 * amp_count);
		return NULL;
	}

	if (out_nocc) *out_nocc = nocc;
	if (out_nvir) *out_nvir = nvir;
	if (out_amp_count) *out_amp_count = amp_count;
	if (out_has_y) *out_has_y = has_y;

    return vlx->rsp.solution_matrix.data + state_idx * vec_size;
}

static void vlx_symmetrize_square(double* mat, size_t dim) {
	ASSERT(mat);
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = i + 1; j < dim; ++j) {
			const double value = 0.5 * (mat[i * dim + j] + mat[j * dim + i]);
			mat[i * dim + j] = value;
			mat[j * dim + i] = value;
		}
	}
}

static void vlx_transform_mo_density_to_ao(double* out_ao, const double* mo_density, const double* coeff, size_t mo_offset, size_t mo_count, size_t num_ao) {
	ASSERT(out_ao);
	ASSERT(mo_density);
	ASSERT(coeff);

	md_temp_scope_t temp = md_temp_begin();
	double* work = md_temp_alloc_array(temp, double, mo_count * num_ao);

	for (size_t i = 0; i < mo_count; ++i) {
		for (size_t ao = 0; ao < num_ao; ++ao) {
			double sum = 0.0;
			for (size_t j = 0; j < mo_count; ++j) {
				sum += mo_density[i * mo_count + j] * coeff[(mo_offset + j) * num_ao + ao];
			}
			work[i * num_ao + ao] = sum;
		}
	}

	for (size_t ao_i = 0; ao_i < num_ao; ++ao_i) {
		for (size_t ao_j = 0; ao_j < num_ao; ++ao_j) {
			double sum = 0.0;
			for (size_t i = 0; i < mo_count; ++i) {
				sum += coeff[(mo_offset + i) * num_ao + ao_i] * work[i * num_ao + ao_j];
			}
			out_ao[ao_i * num_ao + ao_j] = sum;
		}
	}

	md_temp_end(temp);
}

// The whole of the transition density reconstruction, in terms of a plain solution vector and a
// plain AO coefficient matrix rather than a md_vlx_t. This is what lets the attribute provider
// below reconstruct the same matrix with no vlx object in reach: everything it needs is one row of
// a response solution and the (already resident) alpha MO coefficients.
static bool vlx_build_transition_density_matrix(double* out_matrix, const double* solution_vector, size_t vec_len,
	size_t nocc, size_t nvir, const double* coeff, size_t num_ao, md_vlx_transition_type_t type)
{
	ASSERT(out_matrix);
	ASSERT(solution_vector);
	ASSERT(coeff);

	if (nocc == 0 || nvir == 0 || num_ao == 0) {
		return false;
	}

	const size_t amp_count = nocc * nvir;
	bool has_y;
	if (vec_len == amp_count) {
		has_y = false;
	} else if (vec_len == 2 * amp_count) {
		has_y = true;
	} else {
		MD_LOG_ERROR("Response solution vector holds %zu values, expected %zu or %zu (%zu occupied x %zu virtual)", vec_len, amp_count, 2 * amp_count, nocc, nvir);
		return false;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* detach_mo = md_temp_alloc_array(temp, double, nocc * nocc);
	double* attach_mo = md_temp_alloc_array(temp, double, nvir * nvir);
	double* detach_ao = NULL;
	MEMSET(detach_mo, 0, sizeof(double) * nocc * nocc);
	MEMSET(attach_mo, 0, sizeof(double) * nvir * nvir);

	for (size_t i = 0; i < nocc; ++i) {
		for (size_t j = i; j < nocc; ++j) {
			double value = 0.0;
			for (size_t a = 0; a < nvir; ++a) {
				const size_t ia = i * nvir + a;
				const size_t ja = j * nvir + a;
                const double z_i = solution_vector[ia];
                const double z_j = solution_vector[ja];
                const double y_i = has_y ? solution_vector[amp_count + ia] : 0.0;
                const double y_j = has_y ? solution_vector[amp_count + ja] : 0.0;
				const double t_i = z_i - y_i;
				const double t_j = z_j - y_j;
				value += t_i * t_j;
			}
			detach_mo[i * nocc + j] = value;
			detach_mo[j * nocc + i] = value;
		}
	}

	for (size_t a = 0; a < nvir; ++a) {
		for (size_t b = a; b < nvir; ++b) {
			double value = 0.0;
			for (size_t i = 0; i < nocc; ++i) {
				const size_t ia = i * nvir + a;
				const size_t ib = i * nvir + b;
                const double z_i = solution_vector[ia];
                const double z_j = solution_vector[ib];
                const double y_i = has_y ? solution_vector[amp_count + ia] : 0.0;
                const double y_j = has_y ? solution_vector[amp_count + ib] : 0.0;
				const double t_a = z_i - y_i;
				const double t_b = z_j - y_j;
				value += t_a * t_b;
			}
			attach_mo[a * nvir + b] = value;
			attach_mo[b * nvir + a] = value;
		}
	}

	if (type == MD_VLX_TRANSITION_DETACHMENT) {
		vlx_transform_mo_density_to_ao(out_matrix, detach_mo, coeff, 0, nocc, num_ao);
	} else {
		vlx_transform_mo_density_to_ao(out_matrix, attach_mo, coeff, nocc, nvir, num_ao);
		if (type == MD_VLX_TRANSITION_DIFFERENCE) {
			detach_ao = md_temp_alloc_array(temp, double, num_ao * num_ao);
			vlx_transform_mo_density_to_ao(detach_ao, detach_mo, coeff, 0, nocc, num_ao);
			for (size_t i = 0; i < num_ao * num_ao; ++i) {
				out_matrix[i] -= detach_ao[i];
			}
		}
	}

	vlx_symmetrize_square(out_matrix, num_ao);
	md_temp_end(temp);
	return true;
}

static bool vlx_rsp_extract_transition_density_matrix(double* out_matrix, const md_vlx_t* vlx, size_t state_idx, md_vlx_transition_type_t type) {
	ASSERT(out_matrix);
	ASSERT(vlx);

	size_t nocc = 0;
	size_t nvir = 0;
	size_t amp_count = 0;
	bool has_y = false;

	const double* solution_vector = vlx_rsp_get_solution_vector(vlx, state_idx, &nocc, &nvir, &amp_count, &has_y);
	if (!solution_vector) {
		return false;
	}

	const md_vlx_2d_data_t* coeff = &vlx->scf.alpha.coefficients;
	const size_t vec_len = has_y ? 2 * amp_count : amp_count;

	return vlx_build_transition_density_matrix(out_matrix, solution_vector, vec_len, nocc, nvir, coeff->data, coeff->size[1], type);
}

static double vlx_dot(const double* a, const double* b, size_t count) {
	double sum = 0.0;
	for (size_t i = 0; i < count; ++i) {
		sum += a[i] * b[i];
	}
	return sum;
}

static double vlx_normalize(double* vec, size_t count) {
	double norm = sqrt(vlx_dot(vec, vec, count));
	if (norm > 0.0) {
		const double inv_norm = 1.0 / norm;
		for (size_t i = 0; i < count; ++i) {
			vec[i] *= inv_norm;
		}
	}
	return norm;
}

static void vlx_orthogonalize(double* vec, const double* basis, size_t basis_count, size_t dim) {
	for (size_t basis_idx = 0; basis_idx < basis_count; ++basis_idx) {
		const double* b = basis + basis_idx * dim;
		const double projection = vlx_dot(vec, b, dim);
		for (size_t i = 0; i < dim; ++i) {
			vec[i] -= projection * b[i];
		}
	}
}

static bool vlx_symmetric_top_eigenpairs(double* out_values, double* out_vectors, const double* matrix, size_t dim, size_t pair_count) {
	ASSERT(out_values);
	ASSERT(out_vectors);
	ASSERT(matrix);

	if (dim == 0 || pair_count == 0) {
		return false;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* work = md_temp_alloc_array(temp, double, dim * dim);
	double* vec  = md_temp_alloc_array(temp, double, dim);
	double* next = md_temp_alloc_array(temp, double, dim);
	MEMCPY(work, matrix, sizeof(double) * dim * dim);
	MEMSET(out_values, 0, sizeof(double) * pair_count);
	MEMSET(out_vectors, 0, sizeof(double) * pair_count * dim);

	for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
		for (size_t i = 0; i < dim; ++i) {
			vec[i] = 1.0 + 0.013 * (double)(((i + 1) * (pair_idx + 3)) % 17);
		}
		vlx_orthogonalize(vec, out_vectors, pair_idx, dim);
		if (vlx_normalize(vec, dim) <= DBL_EPSILON) {
			vec[pair_idx % dim] = 1.0;
			vlx_orthogonalize(vec, out_vectors, pair_idx, dim);
			vlx_normalize(vec, dim);
		}

		for (size_t iter = 0; iter < VLX_NTO_POWER_ITERATIONS; ++iter) {
			for (size_t i = 0; i < dim; ++i) {
				double sum = 0.0;
				for (size_t j = 0; j < dim; ++j) {
					sum += work[i * dim + j] * vec[j];
				}
				next[i] = sum;
			}

			vlx_orthogonalize(next, out_vectors, pair_idx, dim);
			if (vlx_normalize(next, dim) <= DBL_EPSILON) {
				break;
			}

			double diff = 0.0;
			double neg_diff = 0.0;
			for (size_t i = 0; i < dim; ++i) {
				const double d = next[i] - vec[i];
				const double nd = next[i] + vec[i];
				diff += d * d;
				neg_diff += nd * nd;
				vec[i] = next[i];
			}
			if (sqrt(MIN(diff, neg_diff)) < VLX_NTO_CONVERGENCE_EPSILON) {
				break;
			}
		}

		double eigenvalue = 0.0;
		for (size_t i = 0; i < dim; ++i) {
			double row_sum = 0.0;
			for (size_t j = 0; j < dim; ++j) {
				row_sum += work[i * dim + j] * vec[j];
			}
			eigenvalue += vec[i] * row_sum;
		}

		if (eigenvalue < VLX_NTO_EIGENVALUE_EPSILON) {
			break;
		}

		out_values[pair_idx] = eigenvalue;
		MEMCPY(out_vectors + pair_idx * dim, vec, sizeof(double) * dim);

		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				work[i * dim + j] -= eigenvalue * vec[i] * vec[j];
			}
		}
	}

	md_temp_end(temp);
	return true;
}

static size_t vlx_rsp_extract_nto_from_solution(double* out_coefficients, double* out_lambdas, const md_vlx_t* vlx, size_t state_idx, md_vlx_nto_type_t type, size_t lambda_count) {
	ASSERT(vlx);

	if (lambda_count == 0 || (!out_coefficients && !out_lambdas)) {
		return 0;
	}

	size_t nocc = 0;
	size_t nvir = 0;
	size_t amp_count = 0;
	bool has_y = false;
    const double* solution_vector = vlx_rsp_get_solution_vector(vlx, state_idx, &nocc, &nvir, &amp_count, &has_y);
    if (!solution_vector) {
        return 0;
    }

	const size_t pair_count = MIN(MIN(nocc, nvir), lambda_count);
	if (pair_count == 0) {
		return 0;
	}

	const md_vlx_2d_data_t* scf_coeff = &vlx->scf.alpha.coefficients;
	const size_t num_ao = scf_coeff->size[1];
	const double* coeff = scf_coeff->data;

	if (out_coefficients) {
		MEMSET(out_coefficients, 0, sizeof(double) * lambda_count * num_ao);
	}
	if (out_lambdas) {
		MEMSET(out_lambdas, 0, sizeof(double) * lambda_count);
	}

	const bool use_left_gram = nocc <= nvir;
	const size_t small_dim = use_left_gram ? nocc : nvir;
	const size_t large_dim = use_left_gram ? nvir : nocc;

	size_t written_count = 0;
	md_temp_scope_t temp = md_temp_begin();
	double* transition = md_temp_alloc_array(temp, double, nocc * nvir);
	double* gram = md_temp_alloc_array(temp, double, small_dim * small_dim);
	double* eigenvalues = md_temp_alloc_array(temp, double, pair_count);
	double* small_vectors = md_temp_alloc_array(temp, double, pair_count * small_dim);
	double* large_vectors = md_temp_alloc_array(temp, double, pair_count * large_dim);
	MEMSET(gram, 0, sizeof(double) * small_dim * small_dim);
	MEMSET(large_vectors, 0, sizeof(double) * pair_count * large_dim);

	for (size_t i = 0; i < nocc; ++i) {
		for (size_t a = 0; a < nvir; ++a) {
			const size_t idx = i * nvir + a;
            const double z = solution_vector[idx];
            const double y = has_y ? solution_vector[amp_count + idx] : 0.0;
			transition[idx] = z - y;
		}
	}

	if (use_left_gram) {
		for (size_t i = 0; i < nocc; ++i) {
			for (size_t j = i; j < nocc; ++j) {
				double value = 0.0;
				for (size_t a = 0; a < nvir; ++a) {
					value += transition[i * nvir + a] * transition[j * nvir + a];
				}
				gram[i * nocc + j] = value;
				gram[j * nocc + i] = value;
			}
		}
	} else {
		for (size_t a = 0; a < nvir; ++a) {
			for (size_t b = a; b < nvir; ++b) {
				double value = 0.0;
				for (size_t i = 0; i < nocc; ++i) {
					value += transition[i * nvir + a] * transition[i * nvir + b];
				}
				gram[a * nvir + b] = value;
				gram[b * nvir + a] = value;
			}
		}
	}

	if (!vlx_symmetric_top_eigenpairs(eigenvalues, small_vectors, gram, small_dim, pair_count)) {
		goto done;
	}

	for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
		const double lambda = eigenvalues[pair_idx];
		if (lambda < VLX_NTO_EIGENVALUE_EPSILON) {
			break;
		}

		const double sigma = sqrt(lambda);
		if (sigma <= DBL_EPSILON) {
			break;
		}

		if (out_lambdas) {
			out_lambdas[pair_idx] = lambda;
		}

		if (out_coefficients) {
			const double* small = small_vectors + pair_idx * small_dim;
			double* large = large_vectors + pair_idx * large_dim;

			if (use_left_gram) {
				for (size_t a = 0; a < nvir; ++a) {
					double value = 0.0;
					for (size_t i = 0; i < nocc; ++i) {
						value += transition[i * nvir + a] * small[i];
					}
					large[a] = value / sigma;
				}
			} else {
				for (size_t i = 0; i < nocc; ++i) {
					double value = 0.0;
					for (size_t a = 0; a < nvir; ++a) {
						value += transition[i * nvir + a] * small[a];
					}
					large[i] = value / sigma;
				}
			}
			vlx_normalize(large, large_dim);

			const double* u = use_left_gram ? small : large;
			const double* v = use_left_gram ? large : small;
			double* out_coeff = out_coefficients + pair_idx * num_ao;

			if (type == MD_VLX_NTO_PARTICLE) {
				for (size_t ao = 0; ao < num_ao; ++ao) {
					double value = 0.0;
					for (size_t a = 0; a < nvir; ++a) {
						value += coeff[(nocc + a) * num_ao + ao] * v[a];
					}
					out_coeff[ao] = value;
				}
			} else if (type == MD_VLX_NTO_HOLE) {
				for (size_t ao = 0; ao < num_ao; ++ao) {
					double value = 0.0;
					for (size_t i = 0; i < nocc; ++i) {
						value += coeff[i * num_ao + ao] * u[i];
					}
					out_coeff[ao] = value;
				}
			} else {
				break;
			}
		}

		written_count++;
	}

done:
	md_temp_end(temp);
	return written_count;
}

static size_t vlx_rsp_extract_nto(double* out_coefficients, double* out_lambdas, const md_vlx_t* vlx, size_t state_idx, md_vlx_nto_type_t type, size_t lambda_count) {
	if (!vlx || state_idx >= vlx->rsp.number_of_frequencies || lambda_count == 0) {
		return 0;
	}

	if (vlx->rsp.solution_matrix.data && state_idx < vlx->rsp.solution_matrix.size[0]) {
		return vlx_rsp_extract_nto_from_solution(out_coefficients, out_lambdas, vlx, state_idx, type, lambda_count);
	}

	return 0;
}
static bool vlx_read_scf_results(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	ASSERT(vlx);

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	h5_error_scope_t h5_err_scope = h5_error_scope_begin();

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		h5_error_scope_end(h5_err_scope);
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
	h5_error_scope_end(h5_err_scope);

	return result;
}

// This is the newest version of the file format where everything is contained within a single h5 file
static bool h5_read_xps_data(md_vlx_t* vlx, hid_t handle);

static bool vlx_read_h5_file(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	ASSERT(vlx);

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	h5_error_scope_t h5_err_scope = h5_error_scope_begin();

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		h5_error_scope_end(h5_err_scope);
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

	// VIB
    if (flags & VLX_FLAG_VIB) {
        if (H5Lexists(file_id, "vib", H5P_DEFAULT) > 0) {
            hid_t vib_id = H5Gopen(file_id, "vib", H5P_DEFAULT);
            if (vib_id != H5I_INVALID_HID) {
                result = h5_read_vib_data(vlx, vib_id);
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
				H5Gclose(rsp_id);
				if (!result) goto done;
			}
		}
	}

	// XPS. Optional top level group, independent of the response block: XPS is delta-SCF, so it may
	// appear alongside any md_vlx_rsp_type_t or with no response data at all.
	if (flags & VLX_FLAG_XPS) {
		if (H5Lexists(file_id, "xps", H5P_DEFAULT) > 0) {
			hid_t xps_id = H5Gopen(file_id, "xps", H5P_DEFAULT);
			if (xps_id != H5I_INVALID_HID) {
				result = h5_read_xps_data(vlx, xps_id);
				H5Gclose(xps_id);
				if (!result) goto done;
			}
		}
	}

	if (flags & VLX_FLAG_CORE) {
		if (!h5_read_atomic_properties(vlx, file_id)) {
			goto done;
		}
	}

	if (flags & VLX_FLAG_SCF) {
		if (!h5_read_density_properties(vlx, file_id)) {
			goto done;
		}
	}

	result = true;
done:
	H5Fclose(file_id);
	h5_error_scope_end(h5_err_scope);

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

// XPS
//
// Layout of the optional '/xps' group:
//
//   /xps/<element symbol>/<n>/atom_index            scalar, integer
//                            /contribution          scalar, double
//                            /ionization_energy_ev  scalar, double, unit eV
//                            /is_delocalized        scalar, bool (h5py enum over int8)
//                            /mo_index              scalar, integer
//
// <n> is a flat '0'..'N-1' enumeration with no meaning beyond ordering, so it is discarded --
// entries are re-sorted by (element, ionization_energy) in vlx_xps_finalize().
//
// This pushes into vlx->xps.entries but does NOT finalize. vlx_parse_file() does that once, after
// every push, because finalizing takes pointers into an array that md_array_push may still realloc.
static bool h5_read_xps_data(md_vlx_t* vlx, hid_t handle) {
	ASSERT(vlx);

	H5G_info_t info = { 0 };
	if (H5Gget_info(handle, &info) < 0) {
		MD_LOG_ERROR("XPS: failed to get group info");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		if (H5Gget_objname_by_idx(handle, i, name_buf, sizeof(name_buf)) < 0) {
			continue;
		}
		if (H5Gget_objtype_by_idx(handle, i) != H5G_GROUP) {
			continue;
		}

		const md_element_t element = md_util_element_lookup(str_from_cstr(name_buf), true);
		if (element == 0) {
			MD_LOG_INFO("XPS: skipping group '%s', not a recognized element symbol", name_buf);
			continue;
		}

		hid_t elem_group = H5Gopen(handle, name_buf, H5P_DEFAULT);
		if (elem_group == H5I_INVALID_HID) {
			MD_LOG_ERROR("XPS: failed to open element group '%s'", name_buf);
			continue;
		}

		hsize_t num_links = 0;
		if (H5Gget_num_objs(elem_group, &num_links) < 0) {
			MD_LOG_ERROR("XPS: failed to count entries in element group '%s'", name_buf);
			H5Gclose(elem_group);
			continue;
		}

		// Indexed by name rather than by link index, so entries are visited in numeric order
		// regardless of how HDF5 chose to order the links.
		for (hsize_t j = 0; j < num_links; ++j) {
			char idx_buf[32];
			snprintf(idx_buf, sizeof(idx_buf), "%i", (int)j);
			if (H5Lexists(elem_group, idx_buf, H5P_DEFAULT) <= 0) {
				continue;
			}

			hid_t entry_group = H5Gopen(elem_group, idx_buf, H5P_DEFAULT);
			if (entry_group == H5I_INVALID_HID) {
				continue;
			}

			md_vlx_xps_entry_t entry = {
				.atom_index = -1,
				.mo_index   = -1,
				.element    = element,
			};

			// The only field without a sensible default: an entry with no energy has nothing to plot.
			if (!h5_read_dataset_data(&entry.ionization_energy, 1, entry_group, H5T_NATIVE_DOUBLE, "ionization_energy_ev")) {
				MD_LOG_ERROR("XPS: '%s/%s' has no 'ionization_energy_ev', skipping entry", name_buf, idx_buf);
				H5Gclose(entry_group);
				continue;
			}

			// Remaining fields are optional and keep their defaults if absent.
			h5_read_dataset_data(&entry.contribution, 1, entry_group, H5T_NATIVE_DOUBLE, "contribution");
			h5_read_dataset_data(&entry.atom_index,   1, entry_group, H5T_NATIVE_INT32,  "atom_index");
			h5_read_dataset_data(&entry.mo_index,     1, entry_group, H5T_NATIVE_INT32,  "mo_index");

			// h5py writes Python bools as an HDF5 enum with an int8 base; H5Dread converts
			// enum -> integer for us, so reading it as int8 works for both that and a plain int.
			int8_t is_delocalized = 0;
			if (h5_read_dataset_data(&is_delocalized, 1, entry_group, H5T_NATIVE_INT8, "is_delocalized")) {
				entry.is_delocalized = (is_delocalized != 0);
			}

			if (entry.atom_index >= 0 && (size_t)entry.atom_index >= vlx->number_of_atoms) {
				MD_LOG_INFO("XPS: '%s/%s' has out of range atom_index %i (%zu atoms), clearing it",
					name_buf, idx_buf, (int)entry.atom_index, vlx->number_of_atoms);
				entry.atom_index = -1;
			}

			md_array_push(vlx->xps.entries, entry, vlx->arena);
			H5Gclose(entry_group);
		}

		H5Gclose(elem_group);
	}

	return true;
}

static void vlx_xps_finalize(md_vlx_t* vlx);

// Internal version to control what portions to load
static bool vlx_parse_file(md_vlx_t* vlx, str_t filename, vlx_flags_t flags) {
	md_temp_scope_t temp = md_temp_begin();
	md_allocator_i* temp_alloc = md_temp_allocator(temp);

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

	// Must run after the last md_array_push into vlx->xps.entries, since it converts group offsets
	// into pointers into that array. No-op when the file had no '/xps' group.
	vlx_xps_finalize(vlx);

	if (!str_empty(vlx->basis_set_ident)) {
		size_t cap = KILOBYTES(16);
		char*  buf = md_temp_alloc_array(temp, char, cap);
		md_strb_t sb = md_strb_create(temp_alloc);

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
				goto done;
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
	if ((flags & VLX_FLAG_SCF) && vlx->basis_set.atom_basis.count > 0) {
		size_t num_ao = 0;
		if (vlx->scf.alpha.density.data) {
			num_ao = vlx->scf.alpha.density.size[0];
		} else if (vlx->scf.S.data) {
			num_ao = vlx->scf.S.size[0];
		} else if (vlx->scf.alpha.coefficients.data) {
			num_ao = compute_basis_num_atomic_orbitals(vlx);
			if (num_ao == 0) {
				MD_LOG_ERROR("Unable to infer AO dimension for SCF coefficient normalization");
				goto done;
			}
		}
        md_array_resize(vlx->ao_remap, num_ao, vlx->arena);
		if (!build_ao_remap(vlx->ao_remap, num_ao, vlx)) {
			MD_LOG_ERROR("Failed to build AO remap table");
			goto done;
		}

		if (num_ao > 0 && vlx->ao_remap) {
			// Normalize SCF coefficients to canonical [num_mo x num_ao] in shell order.
			if (!normalize_orbital_coefficients(&vlx->scf.alpha, num_ao, vlx->ao_remap, "Alpha orbital")) {
				goto done;
			}
			if (vlx->scf.alpha.density.data && num_ao == vlx->scf.alpha.density.size[0]) {
				ao_permute_square(vlx->scf.alpha.density.data, num_ao, vlx->ao_remap);
			}
			if (vlx->scf.type == MD_VLX_SCF_UNRESTRICTED) {
				if (!normalize_orbital_coefficients(&vlx->scf.beta, num_ao, vlx->ao_remap, "Beta orbital")) {
					goto done;
				}
				if (vlx->scf.beta.density.data && num_ao == vlx->scf.beta.density.size[0]) {
					ao_permute_square(vlx->scf.beta.density.data, num_ao, vlx->ao_remap);
				}
			}
			else {
				// memcpy again from alpha into beta as dims may have changed.
				MEMCPY(&vlx->scf.beta.coefficients, &vlx->scf.alpha.coefficients, sizeof(md_vlx_2d_data_t));
				MEMCPY(&vlx->scf.beta.density, &vlx->scf.alpha.density, sizeof(md_vlx_2d_data_t));
			}
			if (vlx->scf.S.data && num_ao == vlx->scf.S.size[0]) {
				ao_permute_square(vlx->scf.S.data, num_ao, vlx->ao_remap);
			}

			// Everything is now in shell order and no further AO-basis math is
			// performed, so this is the point to leave VeloxChem's pure/spherical
			// basis for the Cartesian one that md_gto_basis_t requires.
			if (!vlx_convert_ao_data_to_cartesian(vlx)) {
				MD_LOG_ERROR("Failed to convert AO data to the Cartesian basis");
				goto done;
			}
		}
	}

	if ((flags & VLX_FLAG_SCF) && !validate_scf_canonical_layout(vlx)) {
		goto done;
	}

	// Identify homo and lumo
	if (vlx->scf.alpha.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.alpha.occupancy.size; ++i) {
			if (vlx->scf.alpha.occupancy.data[i] == 0.0) {
				vlx->scf.alpha.homo_idx = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.alpha.lumo_idx = i;
				break;
			}
		}
	}

	if (vlx->scf.beta.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.beta.occupancy.size; ++i) {
			if (vlx->scf.beta.occupancy.data[i] == 0.0) {
				vlx->scf.beta.homo_idx = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.beta.lumo_idx = i;
				break;
			}
		}
	}

	// NOTE: ao_to_atom_idx is built inside vlx_convert_ao_data_to_cartesian(), derived
	// from the shell list so it matches the Cartesian AO ordering. The old
	// extract_ao_to_atom_idx() path produced spherical indices and is no longer used
	// here (it is still used to infer the pre-conversion AO count).

	result = true;
done:
	md_temp_end(temp);

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
	return orb->coefficients.size[0];
}

static inline size_t number_of_atomic_orbitals(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline size_t number_of_ao_coefficients(const md_vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline void extract_ao_coefficients(double* out_coeff, const md_vlx_orbital_t* orb, size_t ao_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);
	ASSERT(ao_idx < number_of_atomic_orbitals(orb));

	extract_col(out_coeff, &orb->coefficients, ao_idx);
}

const double* md_vlx_scf_resp_charges(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.resp_charges;
	}
	return NULL;
}

size_t md_vlx_rsp_number_of_excited_states(const md_vlx_t* vlx) {
	if (vlx) {
		if (vlx->rsp.type == MD_VLX_RSP_LINEAR) {
			return vlx->rsp.number_of_frequencies;
		}
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

md_vlx_rsp_type_t md_vlx_rsp_type(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.type;
	}
	return MD_VLX_RSP_UNKNOWN;
}

double md_vlx_c6_value(const md_vlx_t* vlx) {
	if (vlx && vlx->rsp.type == MD_VLX_RSP_C6) {
		return vlx->rsp.c6;
	}
	return 0.0;
}

size_t md_vlx_rsp_number_of_frequencies(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.number_of_frequencies;
	}
	return 0;
}

const double* md_vlx_rsp_frequencies(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.frequencies;
	}
	return NULL;
}

const double* md_vlx_rsp_sigma(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.sigmas;
	}
	return NULL;
}

const double* md_vlx_rsp_delta_epsilons(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->rsp.delta_epsilons;
    }
    return NULL;
}

const double* md_vlx_rsp_optical_rotations(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.optical_rotations;
	}
	return NULL;
}

const double* md_vlx_rsp_tpa_trans_linear(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.tpa_strengths_linear;
	}
	return NULL;
}

const double* md_vlx_rsp_tpa_trans_circular(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.tpa_strengths_circular;
	}
	return NULL;
}

const double* md_vlx_rsp_tpa_cross_sections(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.cross_sections;
	}
	return NULL;
}

const double* md_vlx_rsp_tpa_gamma_re(const md_vlx_t* vlx) {
	return NULL;
}

const double* md_vlx_rsp_tpa_gamma_im(const md_vlx_t* vlx) {
	return NULL;
}

// RIXS
// @TODO: The RIXS datasets 'photon_energies', 'elastic_cross_sections', 'emission_energies',
// 'energy_losses' and the 'gamma_fwhm_ev' attribute are not yet read in vlx_parse_rsp().
// 'num_incomming_photons' and 'num_final_states' also need to be assigned from the dataset dims.
// The accessors below are thin and will start returning valid data as soon as that is in place.

static inline bool vlx_is_rixs(const md_vlx_t* vlx) {
	return vlx && vlx->rsp.type == MD_VLX_RSP_RIXS;
}

size_t md_vlx_rsp_rixs_number_of_photon_energies(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.num_incomming_photons : 0;
}

size_t md_vlx_rsp_rixs_number_of_final_states(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.num_final_states : 0;
}

size_t md_vlx_rsp_rixs_number_of_core_states(const md_vlx_t* vlx) {
	// Core eigenvalues are currently stored in rsp.frequencies for RIXS
	return vlx_is_rixs(vlx) ? vlx->rsp.number_of_frequencies : 0;
}

const double* md_vlx_rsp_rixs_photon_energies(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.photon_energies : NULL;
}

const double* md_vlx_rsp_rixs_elastic_cross_sections(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.elastic_cross_sections : NULL;
}

const double* md_vlx_rsp_rixs_cross_sections(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.cross_sections : NULL;
}

const double* md_vlx_rsp_rixs_energy_losses(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.energy_losses : NULL;
}

const double* md_vlx_rsp_rixs_emission_energies(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.emission_energies : NULL;
}

const double* md_vlx_rsp_rixs_core_eigenvalues(const md_vlx_t* vlx) {
	if (!vlx_is_rixs(vlx)) return NULL;
	// Prefer the dedicated array if it has been populated, otherwise fall back to rsp.frequencies,
	// which is where 'core_eigenvalues' is currently parsed into.
	return vlx->rsp.rixs.core_eigenvalues ? vlx->rsp.rixs.core_eigenvalues : vlx->rsp.frequencies;
}

const double* md_vlx_rsp_rixs_core_osc_strengths(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.core_osc_strengths : NULL;
}

double md_vlx_rsp_rixs_gamma_fwhm_ev(const md_vlx_t* vlx) {
	return vlx_is_rixs(vlx) ? vlx->rsp.rixs.gamma_fwhm_ev : 0.0;
}

bool md_vlx_rsp_has_nto(const md_vlx_t* vlx) {
	if (!vlx) return false;
	if (vlx->rsp.solution_matrix.data && vlx->rsp.solution_matrix.size[0] == vlx->rsp.number_of_frequencies) {
		return true;
	}
	return false;
}

size_t md_vlx_rsp_nto_lambdas_extract(double* out_lambdas, const md_vlx_t* vlx, size_t state_idx, size_t lambda_count) {
	return vlx_rsp_extract_nto(NULL, out_lambdas, vlx, state_idx, MD_VLX_NTO_PARTICLE, lambda_count);
}

size_t md_vlx_rsp_nto_coefficients_extract(double* out_coefficients, double* out_lambdas, const md_vlx_t* vlx, size_t state_idx, md_vlx_nto_type_t type, size_t lambda_count) {
	return vlx_rsp_extract_nto(out_coefficients, out_lambdas, vlx, state_idx, type, lambda_count);
}

size_t md_vlx_rsp_transition_density_matrix_size(const md_vlx_t* vlx, size_t state_idx) {
	if (!vlx || state_idx >= vlx->rsp.number_of_frequencies) {
		return 0;
	}

	size_t nocc = 0;
	size_t nvir = 0;
	if (!vlx_rsp_get_solution_vector(vlx, state_idx, &nocc, &nvir, NULL, NULL)) {
		return 0;
	}
	(void)nocc;
	(void)nvir;

	return vlx->scf.alpha.coefficients.size[1];
}

size_t md_vlx_rsp_transition_density_matrix_extract(double* out_matrix, const md_vlx_t* vlx, size_t state_idx, md_vlx_transition_type_t type) {
	if (!out_matrix || !vlx || state_idx >= vlx->rsp.number_of_frequencies) {
		return 0;
	}

	const size_t dim = md_vlx_rsp_transition_density_matrix_size(vlx, state_idx);
	if (dim == 0) {
		return 0;
	}

	if (!vlx_rsp_extract_transition_density_matrix(out_matrix, vlx, state_idx, type)) {
		return 0;
	}

	return dim;
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

md_vlx_opt_type_t md_vlx_opt_type(const md_vlx_t* vlx) {
	return vlx->opt.type;
}

size_t md_vlx_opt_state_index(const md_vlx_t* vlx) {
	return vlx->opt.state_index;
}

size_t md_vlx_opt_irc_ts_index(const md_vlx_t* vlx) {
	return vlx->opt.ts_index;
}

size_t md_vlx_opt_number_of_steps(const struct md_vlx_t* vlx) {
	if (vlx) {
		return vlx->opt.number_of_steps;
	}
	return 0;
}

// Returns atom coordinates for a given optimization step
const dvec3_t* md_vlx_opt_coordinates(const struct md_vlx_t* vlx, size_t opt_idx) {
	if (vlx) {
		if (vlx->opt.coordinates && opt_idx < vlx->opt.number_of_steps) {
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

size_t md_vlx_vib_number_of_external_frequencies(const md_vlx_t* vlx) {
    if (vlx) {
        return vlx->vib.num_external_frequencies;
	}
    return 0;
}

const double* md_vlx_vib_external_frequencies(const md_vlx_t* vlx) {
	if (vlx) {
		return vlx->vib.external_frequencies;
	}
	return NULL;
}

const double* md_vlx_vib_raman_activities(const md_vlx_t* vlx, size_t idx) {
    if (vlx) {
        if (vlx->vib.raman_activities && idx < vlx->vib.num_external_frequencies) {
            return vlx->vib.raman_activities + idx * vlx->vib.number_of_normal_modes;
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

// XPS

static int vlx_xps_compare_entry(const void* a, const void* b) {
	const md_vlx_xps_entry_t* ea = (const md_vlx_xps_entry_t*)a;
	const md_vlx_xps_entry_t* eb = (const md_vlx_xps_entry_t*)b;
	if (ea->element != eb->element) {
		return (ea->element < eb->element) ? -1 : 1;
	}
	if (ea->ionization_energy != eb->ionization_energy) {
		return (ea->ionization_energy < eb->ionization_energy) ? -1 : 1;
	}
	return 0;
}

// Sorts entries into contiguous per element runs and materializes the public group views.
// MUST be the last operation that touches vlx->xps.entries: md_array_push reallocs, so every pointer
// written here is invalidated by any subsequent push.
static void vlx_xps_finalize(md_vlx_t* vlx) {
	const size_t num_entries = md_array_size(vlx->xps.entries);

	md_array_shrink(vlx->xps.groups_internal, 0);
	md_array_shrink(vlx->xps.groups, 0);

	if (num_entries == 0) {
		return;
	}

	qsort(vlx->xps.entries, num_entries, sizeof(md_vlx_xps_entry_t), vlx_xps_compare_entry);

	// Phase 1: derive the group runs, offsets only.
	for (size_t i = 0; i < num_entries; ++i) {
		vlx_xps_group_internal_t* last = md_array_last(vlx->xps.groups_internal);
		if (last && last->element == vlx->xps.entries[i].element) {
			last->count += 1;
		} else {
			vlx_xps_group_internal_t grp = {
				.element = vlx->xps.entries[i].element,
				.offset  = (uint32_t)i,
				.count   = 1,
			};
			md_array_push(vlx->xps.groups_internal, grp, vlx->arena);
		}
	}

	// Phase 2: offsets -> pointers. Safe only because entries is no longer being grown.
	const size_t num_groups = md_array_size(vlx->xps.groups_internal);
	md_array_resize(vlx->xps.groups, num_groups, vlx->arena);
	for (size_t i = 0; i < num_groups; ++i) {
		const vlx_xps_group_internal_t* src = &vlx->xps.groups_internal[i];
		vlx->xps.groups[i].element = src->element;
		vlx->xps.groups[i].count   = src->count;
		vlx->xps.groups[i].entries = vlx->xps.entries + src->offset;
	}
}

bool md_vlx_has_xps(const md_vlx_t* vlx) {
	return vlx && md_array_size(vlx->xps.entries) > 0;
}

size_t md_vlx_xps_count(const md_vlx_t* vlx) {
	return vlx ? md_array_size(vlx->xps.entries) : 0;
}

const md_vlx_xps_entry_t* md_vlx_xps_entries(const md_vlx_t* vlx) {
	return vlx ? vlx->xps.entries : NULL;
}

size_t md_vlx_xps_group_count(const md_vlx_t* vlx) {
	return vlx ? md_array_size(vlx->xps.groups) : 0;
}

const md_vlx_xps_group_t* md_vlx_xps_group_by_index(const md_vlx_t* vlx, size_t idx) {
	if (vlx && idx < md_array_size(vlx->xps.groups)) {
		return &vlx->xps.groups[idx];
	}
	return NULL;
}

const md_vlx_xps_group_t* md_vlx_xps_group_by_element(const md_vlx_t* vlx, md_element_t element) {
	if (vlx) {
		// Group count is the number of distinct elements in the calculation, typically 1-5.
		for (size_t i = 0; i < md_array_size(vlx->xps.groups); ++i) {
			if (vlx->xps.groups[i].element == element) {
				return &vlx->xps.groups[i];
			}
		}
	}
	return NULL;
}

// Builds an attribute path from a fixed group prefix and a name taken from the file. A '/' inside
// the name would silently introduce a group level in a namespace where the separator is the only
// structure there is, so it is folded to '_'. Returns an empty str_t when the name does not fit,
// which the callers treat as "skip this one" rather than as a reason to stop publishing.
static str_t vlx_attribute_path(char* buf, size_t cap, str_t group, str_t name) {
	int len = snprintf(buf, cap, STR_FMT "/" STR_FMT, STR_ARG(group), STR_ARG(name));
	if (len <= 0 || (size_t)len >= cap) {
		MD_LOG_ERROR("Attribute path '" STR_FMT "/" STR_FMT "' does not fit in %zu characters", STR_ARG(group), STR_ARG(name), cap - 1);
		return (str_t){0};
	}
	for (int c = (int)group.len + 1; c < len; ++c) {
		if (buf[c] == '/') buf[c] = '_';
	}
	return str_from_cstrn(buf, (size_t)len);
}

// Publishes one attribute under a path this publisher owns, replacing whatever was there - see
// md_attributes_replace on why that is what a producer wants.
static md_attribute_id_t vlx_publish(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_format_t format, const void* data, size_t byte_size) {
	return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
		.path      = path,
		.format    = format,
		.unit      = unit,
		.label     = label,
		.data      = data,
		.byte_size = byte_size,
	});
}

// The same, for an attribute computed through a provider instead of one copied in. The provider's
// user_data is 'sys' itself (see the transition density providers below), a borrowed pointer that
// needs no bookkeeping and outlives 'vlx'.
static md_attribute_id_t vlx_publish_virtual(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_format_t format, const md_attribute_virtual_t* virt) {
	return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
		.path   = path,
		.format = format,
		.unit   = unit,
		.label  = label,
		.virt   = virt,
	});
}

// rank 1 {N}, one scalar per element.
static md_attribute_id_t vlx_publish_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t count) {
	if (!values || count == 0) {
		return MD_ATTRIBUTE_INVALID;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count },
	};
	return vlx_publish(sys, path, label, unit, format, values, count * sizeof(double));
}

// Gives an already published attribute a second, format neutral name. Both names then read one
// datum - no copy, and a consumer of either is unaffected when the other appears or goes.
static md_attribute_id_t vlx_alias(md_system_t* sys, md_attribute_id_t target, str_t path) {
	if (target != MD_ATTRIBUTE_INVALID) {
		const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
		if (existing) {
			md_attributes_remove(&sys->attributes, existing->id);
		}
		return md_attributes_alias(&sys->attributes, target, path, (str_t){0}, (str_t){0});
	}
	return MD_ATTRIBUTE_INVALID;
}

// Publishes beta's copy of a per orbital series, or a SECOND NAME for alpha's when the two share
// storage. md_vlx shallow copies beta from alpha for anything but an unrestricted calculation, so
// comparing the POINTERS is what tells the cases apart - and it is the only test that gets the
// restricted open shell case right, where the orbitals are shared but the occupations are read
// separately. Switching on md_vlx_scf_type() instead would alias an occupation array that differs.
static md_attribute_id_t vlx_publish_or_alias(md_system_t* sys, md_attribute_id_t alpha_id, str_t path, str_t label, md_unit_t unit,
                                              const double* alpha_values, const double* beta_values, size_t count) {
	if (beta_values && beta_values == alpha_values) {
		return vlx_alias(sys, alpha_id, path);
	}
	return vlx_publish_series(sys, path, label, unit, beta_values, count);
}

// rank 0, a single scalar. The value is copied, so a local is fine.
static void vlx_publish_scalar(md_system_t* sys, str_t path, str_t label, md_unit_t unit, double value) {
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 0,
	};
	vlx_publish(sys, path, label, unit, format, &value, sizeof(double));
}

// rank 1 {N} of 3 component values. dvec3_t is three contiguous doubles, so the source array is
// already the interleaved layout an attribute stores and this is a straight copy.
static void vlx_publish_vec3_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const dvec3_t* values, size_t count) {
	if (!values || count == 0) {
		return;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 1, .shape = { (uint32_t)count },
	};
	vlx_publish(sys, path, label, unit, format, values, count * 3 * sizeof(double));
}

// A single string is rank 1 {1}, by the same rule that makes a single 3-vector rank 2 {1,3}. The
// descriptor carries the TEXT and the table stores a handle - see the STRINGS note in md_system.h.
static void vlx_publish_str(md_system_t* sys, str_t path, str_t label, str_t value) {
	if (str_empty(value)) {
		return;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { 1 },
	};
	vlx_publish(sys, path, label, md_unit_none(), format, &value, sizeof(str_t));
}

// rank 2 {A,B}, one scalar per (a,b), row major with b fastest - which is the layout md_vlx.h
// documents for the 2D response quantities, so no rearrangement happens here.
static void vlx_publish_matrix(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t rows, size_t cols) {
	if (!values || rows == 0 || cols == 0) {
		return;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)rows, (uint32_t)cols },
	};
	vlx_publish(sys, path, label, unit, format, values, rows * cols * sizeof(double));
}

// Publishes one COLUMN of an array of structs. The source is strided and an attribute is
// contiguous, so the values are gathered into the table's own storage through md_attributes_data
// rather than through a temporary which is then copied again.
//
// A struct of mixed types is not one attribute - a value has ONE type - so a record with six fields
// becomes six sibling paths over the same index space. That is what the ATTRIBUTES note means by
// independent quantities being sibling paths: the transposition from the file's row layout is the
// whole cost, and it happens once, here.
static void vlx_publish_column(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_type_t type, const void* base, size_t stride, size_t count) {
	if (!base || count == 0) {
		return;
	}

	md_attribute_format_t format = {
		.type = type, .components = 1, .rank = 1, .shape = { (uint32_t)count },
	};
	md_attribute_id_t id = vlx_publish(sys, path, label, unit, format, NULL, 0);
	if (id == MD_ATTRIBUTE_INVALID) {
		return;
	}

	uint8_t* dst = (uint8_t*)md_attributes_data(&sys->attributes, id, type);
	if (!dst) {
		md_attributes_remove(&sys->attributes, id);
		return;
	}

	const size_t   elem_size = md_attribute_type_size(type);
	const uint8_t* src       = (const uint8_t*)base;
	for (size_t i = 0; i < count; ++i) {
		MEMCPY(dst + i * elem_size, src + i * stride, elem_size);
	}
}

typedef const dvec3_t* (*vlx_vec3_row_fn)(const md_vlx_t* vlx, size_t row);

// rank 2 {A,B} of 3 component values, assembled from an accessor which hands back one row at a
// time. Created empty and filled through md_attributes_data, so each row is written straight into
// the table's storage instead of into a temporary which is then copied again.
static void vlx_publish_vec3_rows(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const md_vlx_t* vlx, size_t num_rows, size_t row_len, vlx_vec3_row_fn row_fn) {
	if (num_rows == 0 || row_len == 0 || !row_fn) {
		return;
	}

	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 2, .shape = { (uint32_t)num_rows, (uint32_t)row_len },
	};
	md_attribute_id_t id = vlx_publish(sys, path, label, unit, format, NULL, 0);
	if (id == MD_ATTRIBUTE_INVALID) {
		return;
	}

	double* dst = (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64);
	if (!dst) {
		md_attributes_remove(&sys->attributes, id);
		return;
	}

	for (size_t i = 0; i < num_rows; ++i) {
		const dvec3_t* row = row_fn(vlx, i);
		if (!row) {
			// A missing row would leave the attribute half written and zero filled, which reads as
			// data rather than as an absence. Better that the path is simply not there.
			MD_LOG_ERROR("Missing row %zu while publishing '" STR_FMT "'", i, STR_ARG(path));
			md_attributes_remove(&sys->attributes, id);
			return;
		}
		MEMCPY(dst + i * row_len * 3, row, row_len * 3 * sizeof(double));
	}
}

// The anchor of a dipole group: rank 0, one 3 component value, constant over whatever index space
// the group's vector has. Angstrom because it is a point in system space, unlike the moment.
//
// It is NOT replicated to match the vector's shape. Group members share an index space, not a
// shape; storing the same three numbers once per excited state would be N copies with nothing
// keeping them equal, to save a consumer one line.
static void vlx_publish_origin(md_system_t* sys, str_t path, dvec3_t origin) {
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 0,
	};
	vlx_publish(sys, path, (str_t){0}, md_unit_angstrom(), format, &origin, 3 * sizeof(double));
}

// The centre of charge, in Angstrom: where a dipole moment is drawn from. Nuclear charge weighted
// positions less the electronic contribution, per electron. Returns false when the file does not
// carry what it takes to compute one, in which case no dipole group is published at all - half a
// group is not a dipole anyone can draw.
static bool vlx_centre_of_charge(dvec3_t* out_angstrom, const md_vlx_t* vlx) {
	const size_t   num_atoms     = md_vlx_number_of_atoms(vlx);
	const dvec3_t* atom_coord    = md_vlx_atom_coordinates(vlx);
	const uint8_t* atomic_number = md_vlx_atomic_numbers(vlx);

	if (num_atoms == 0 || !atom_coord || !atomic_number) {
		return false;
	}

	const size_t num_electrons = md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_ALPHA) + md_vlx_number_of_electrons(vlx, MD_VLX_SPIN_BETA);
	if (num_electrons == 0) {
		return false;
	}

	// Coordinates are Angstrom while the moment is atomic units, so the nuclear term is taken to
	// bohr first and the result taken back at the end.
	double nx = 0.0, ny = 0.0, nz = 0.0;
	for (size_t i = 0; i < num_atoms; ++i) {
		const double z = (double)atomic_number[i];
		nx += atom_coord[i].x * ANGSTROM_TO_BOHR * z;
		ny += atom_coord[i].y * ANGSTROM_TO_BOHR * z;
		nz += atom_coord[i].z * ANGSTROM_TO_BOHR * z;
	}

	const dvec3_t moment = md_vlx_scf_ground_state_dipole_moment(vlx);
	const double  inv_ne = 1.0 / (double)num_electrons;

	out_angstrom->x = (nx - moment.x) * inv_ne * BOHR_TO_ANGSTROM;
	out_angstrom->y = (ny - moment.y) * inv_ne * BOHR_TO_ANGSTROM;
	out_angstrom->z = (nz - moment.z) * inv_ne * BOHR_TO_ANGSTROM;
	return true;
}

// Reconstructs one excited state's AO-basis transition density purely from attributes already on
// 'sys' - the response solution vectors and the alpha MO coefficients - so this runs with no vlx
// object in reach, and keeps working after one is torn down. 'type' is baked in by the three thin
// wrappers below, one per sibling path, the same shape as the MO coefficient split above.
//
// A whole (unsliced) request reconstructs every state one after another into 'dst'; expensive, but
// no more so than the vlx-based accessor doing the same loop, and slicing by state is how a caller
// avoids paying for states it does not need.
static size_t vlx_transition_density_provide(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_vlx_transition_type_t type) {
	md_system_t* sys = (md_system_t*)user_data;

	const md_attribute_t* sol   = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/solution_matrix"));
	const md_attribute_t* coeff = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/coefficient"));
	const md_attribute_t* core  = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/num_core"));
	const md_attribute_t* val   = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/num_valence"));
	const md_attribute_t* vir   = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/num_virtual"));
	if (!sol || !coeff || !core || !val || !vir) {
		MD_LOG_ERROR("'" STR_FMT "' is missing the response data it reconstructs from", STR_ARG(attr->path));
		return 0;
	}

	double num_core = 0.0, num_valence = 0.0, num_virtual = 0.0;
	md_attribute_extract_f64(&num_core,    1, core, md_unit_none());
	md_attribute_extract_f64(&num_valence, 1, val,  md_unit_none());
	md_attribute_extract_f64(&num_virtual, 1, vir,  md_unit_none());

	const size_t nocc       = num_core > 0.0 ? (size_t)num_core : (size_t)num_valence;
	const size_t nvir       = (size_t)num_virtual;
	const size_t num_states = sol->format.shape[0];
	const size_t vec_len    = sol->format.shape[1];
	const size_t num_ao     = coeff->format.shape[1];
	const size_t num_mo     = coeff->format.shape[0];

	// A provider that returns 0 is reported by the extract as "wrote 0 of N", which says nothing
	// about WHICH input was missing. Each of these is a distinct, actionable state, so each says so.
	if (nocc == 0 || nvir == 0 || num_ao == 0 || num_states == 0) {
		MD_LOG_ERROR("'" STR_FMT "': %zu occupied, %zu virtual, %zu atomic orbitals, %zu states - none of these may be zero",
			STR_ARG(attr->path), nocc, nvir, num_ao, num_states);
		return 0;
	}

	// The destination is sized by the CALLER's slice, and this writes num_ao^2 per state. Nothing
	// upstream guarantees the two agree - the shape was published from
	// md_vlx_scf_number_of_atomic_orbitals() while num_ao here comes off the coefficient matrix -
	// so disagreeing is an overrun, not a wrong picture. Check before writing a single value.
	const size_t states_written = (slice && slice->num_idx > 0) ? 1 : num_states;
	const size_t needed = states_written * num_ao * num_ao;
	if (needed != cap) {
		MD_LOG_ERROR("'" STR_FMT "': asked for %zu values, would write %zu (%zu state(s) x %zu x %zu atomic orbitals)",
			STR_ARG(attr->path), cap, needed, states_written, num_ao, num_ao);
		return 0;
	}

	md_temp_scope_t temp = md_temp_begin();
	size_t result = 0;

	double* coeff_data = md_temp_alloc_array(temp, double, num_mo * num_ao);
	double* row        = md_temp_alloc_array(temp, double, vec_len);
	if (!coeff_data || !row) {
		MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), num_mo * num_ao + vec_len);
		goto done;
	}

	// Every step below gets its own report. An '&&' chain here costs nothing to write and tells a
	// reader of the log only that the provider produced nothing, which is exactly the position this
	// was debugged from - and there are five separate attributes it can be let down by.
	if (md_attribute_extract_f64(coeff_data, num_mo * num_ao, coeff, md_unit_none()) != num_mo * num_ao) {
		MD_LOG_ERROR("'" STR_FMT "': could not read the %zu x %zu coefficients from '" STR_FMT "'",
			STR_ARG(attr->path), num_mo, num_ao, STR_ARG(coeff->path));
		goto done;
	}

	for (size_t s = 0; s < states_written; ++s) {
		const uint32_t state_idx = (slice && slice->num_idx > 0) ? slice->idx[0] : (uint32_t)s;
		const md_attribute_slice_t row_slice = md_attribute_slice_1(state_idx);

		if (md_attribute_extract_slice_f64(row, vec_len, sol, &row_slice, md_unit_none()) != vec_len) {
			MD_LOG_ERROR("'" STR_FMT "': could not read state %u's %zu element solution vector from '" STR_FMT "' (%u x %u)",
				STR_ARG(attr->path), state_idx, vec_len, STR_ARG(sol->path), sol->format.shape[0], sol->format.shape[1]);
			goto done;
		}
		if (!vlx_build_transition_density_matrix((double*)dst + s * num_ao * num_ao, row, vec_len, nocc, nvir, coeff_data, num_ao, type)) {
			MD_LOG_ERROR("'" STR_FMT "': could not build state %u's density from %zu occupied x %zu virtual over %zu atomic orbitals",
				STR_ARG(attr->path), state_idx, nocc, nvir, num_ao);
			goto done;
		}
	}
	result = cap;

done:
	md_temp_end(temp);
	return result;
}

static size_t vlx_transition_density_attachment_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	return vlx_transition_density_provide(dst, cap, attr, slice, user_data, MD_VLX_TRANSITION_ATTACHMENT);
}

static size_t vlx_transition_density_detachment_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	return vlx_transition_density_provide(dst, cap, attr, slice, user_data, MD_VLX_TRANSITION_DETACHMENT);
}

static size_t vlx_transition_density_difference_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	return vlx_transition_density_provide(dst, cap, attr, slice, user_data, MD_VLX_TRANSITION_DIFFERENCE);
}

// D[ao_i][ao_j] = sum_mo occ[mo] * C[mo][ao_i] * C[mo][ao_j] - the definition of the one particle
// density matrix in terms of the (possibly fractionally occupied) orbitals that produced it, so
// this reconstructs the ground state density EXACTLY from data already resident as attributes,
// with no separate density read needed. A zero occupied MO contributes nothing, so it is skipped
// rather than paid for.
static bool vlx_build_occupation_density_matrix(double* out_matrix, const double* coeff, const double* occ, size_t num_mo, size_t num_ao) {
	ASSERT(out_matrix);
	ASSERT(coeff);
	ASSERT(occ);

	if (num_mo == 0 || num_ao == 0) {
		return false;
	}

	MEMSET(out_matrix, 0, sizeof(double) * num_ao * num_ao);
	for (size_t mo = 0; mo < num_mo; ++mo) {
		const double w = occ[mo];
		if (w == 0.0) continue;
		const double* c = coeff + mo * num_ao;
		for (size_t ao_i = 0; ao_i < num_ao; ++ao_i) {
			const double wci = w * c[ao_i];
			if (wci == 0.0) continue;
			for (size_t ao_j = 0; ao_j < num_ao; ++ao_j) {
				out_matrix[ao_i * num_ao + ao_j] += wci * c[ao_j];
			}
		}
	}
	return true;
}

// Reconstructs one spin's ground state AO density from the MO coefficients and occupations this
// system already carries as attributes - the same 'no vlx object in reach' shape as the transition
// density providers above, sharing the same rationale: recomputing on demand costs one pass over
// the occupied orbitals instead of holding a second, redundant [A][A] copy alongside them.
static size_t vlx_scf_density_provide(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, md_vlx_spin_t spin) {
	md_system_t* sys = (md_system_t*)user_data;

	str_t coeff_path = spin == MD_VLX_SPIN_ALPHA ? STR_LIT("orbital/alpha/coefficient")        : STR_LIT("orbital/beta/coefficient");
	str_t occ_path   = spin == MD_VLX_SPIN_ALPHA ? STR_LIT("vlx/scf/orbital/alpha/occupation") : STR_LIT("vlx/scf/orbital/beta/occupation");

	const md_attribute_t* coeff = md_attributes_find(&sys->attributes, coeff_path);
	const md_attribute_t* occ   = md_attributes_find(&sys->attributes, occ_path);
	if (!coeff || !occ) {
		MD_LOG_ERROR("'" STR_FMT "' is missing the orbital data it reconstructs from", STR_ARG(attr->path));
		return 0;
	}

	const size_t num_mo = coeff->format.shape[0];
	const size_t num_ao = coeff->format.shape[1];
	if (occ->format.shape[0] != num_mo) {
		MD_LOG_ERROR("'" STR_FMT "': occupation holds %u values, coefficients hold %zu orbitals", STR_ARG(attr->path), occ->format.shape[0], num_mo);
		return 0;
	}

	if (num_ao * num_ao != cap) {
		MD_LOG_ERROR("'" STR_FMT "' was asked for %zu values and its coefficients span %zu atomic orbitals", STR_ARG(attr->path), cap, num_ao);
		return 0;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* coeff_data = md_temp_alloc_array(temp, double, num_mo * num_ao);
	double* occ_data   = md_temp_alloc_array(temp, double, num_mo);

	bool ok = coeff_data && occ_data;
	if (!ok) {
		MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), num_mo * num_ao + num_mo);
	}
	ok = ok && md_attribute_extract_f64(coeff_data, num_mo * num_ao, coeff, md_unit_none()) == num_mo * num_ao
	        && md_attribute_extract_f64(occ_data, num_mo, occ, md_unit_none()) == num_mo
	        && vlx_build_occupation_density_matrix((double*)dst, coeff_data, occ_data, num_mo, num_ao);

	md_temp_end(temp);
	return ok ? cap : 0;
}

static size_t vlx_scf_alpha_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice; // one whole {A,A} matrix, not indexed by anything a slice could fix
	return vlx_scf_density_provide(dst, cap, attr, user_data, MD_VLX_SPIN_ALPHA);
}

static size_t vlx_scf_beta_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice;
	return vlx_scf_density_provide(dst, cap, attr, user_data, MD_VLX_SPIN_BETA);
}

// alpha +/- beta. Both halves are themselves virtual, so this is a derivation over derivations -
// legal because the graph stays acyclic, and the case the acyclicity rule in md_system.h exists
// for. A restricted calculation gets this for free: beta is an ALIAS of alpha there, so the total
// comes out as twice alpha and the difference as zero without a special case anywhere.
//
// Each read rebuilds both matrices. That is the caller's to cache, and a representation already
// keys one on the volume hash it computes from its own settings.
static size_t vlx_scf_density_combine(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, double beta_scale) {
	md_system_t* sys = (md_system_t*)user_data;

	const md_attribute_t* alpha = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/density"));
	const md_attribute_t* beta  = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/density"));
	if (!alpha || !beta) {
		MD_LOG_ERROR("'" STR_FMT "' is missing a spin density it combines", STR_ARG(attr->path));
		return 0;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* beta_data = md_temp_alloc_array(temp, double, cap);

	bool ok = beta_data != NULL;
	if (!ok) {
		MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), cap);
	}
	if (ok && md_attribute_extract_f64((double*)dst, cap, alpha, md_unit_none()) != cap) {
		MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(alpha->path));
		ok = false;
	}
	if (ok && md_attribute_extract_f64(beta_data, cap, beta, md_unit_none()) != cap) {
		MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(beta->path));
		ok = false;
	}

	if (ok) {
		double* out = (double*)dst;
		for (size_t i = 0; i < cap; ++i) {
			out[i] += beta_scale * beta_data[i];
		}
	}

	md_temp_end(temp);
	return ok ? cap : 0;
}

static size_t vlx_scf_total_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice; // one whole {A,A} matrix, not indexed by anything a slice could fix
	return vlx_scf_density_combine(dst, cap, attr, user_data, 1.0);
}

static size_t vlx_scf_difference_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice;
	return vlx_scf_density_combine(dst, cap, attr, user_data, -1.0);
}

void md_vlx_publish_attributes(md_system_t* sys, const md_vlx_t* vlx) {
	ASSERT(sys);

	if (!vlx) {
		return;
	}
	if (!sys->attributes.alloc) {
		MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
		return;
	}

	// Not in md_unit.h because nothing outside quantum chemistry asks for them, and a unit which is
	// only ever constructed in one place is better constructed there than named globally.
	const md_unit_t hartree    = md_unit_hartree();
	const md_unit_t e_bohr     = md_unit_elementary_charge_bohr();
	const md_unit_t angstrom   = md_unit_angstrom();
	const md_unit_t wavenumber = md_unit_pow(md_unit_scl(md_unit_meter(), 1.0e-2), -1);      // cm^-1
	const md_unit_t km_per_mol = md_unit_div(md_unit_scl(md_unit_meter(), 1.0e3), md_unit_mole());
	const md_unit_t amu        = md_unit_scl(md_unit_kilogram(), 1.66053906660e-27);

	// A label is only carried where the leaf cannot spell it: a consumer prettifies the last path
	// segment when there is none, so "gradient_norm" needs no help and "ir_intensity" does.

	// ---- Molecule level scalars. rank 0 is a single value, not an array of one. ----
	vlx_publish_scalar(sys, STR_LIT("vlx/molecular_charge"),          STR_LIT("Molecular Charge"),			md_unit_none(), md_vlx_molecular_charge(vlx));
	vlx_publish_scalar(sys, STR_LIT("vlx/nuclear_repulsion_energy"),  STR_LIT("Nuclear Repulsion Energy"),	hartree,        md_vlx_nuclear_repulsion_energy(vlx));
	// The two facts about a calculation that are TEXT and nothing else - no consumer can derive them
	// from the columns, the way it can derive the SCF type from whether the spin channels share data.
	vlx_publish_str(sys, STR_LIT("vlx/basis_set"),      STR_LIT("Basis Set"),      md_vlx_basis_set_ident(vlx));
	vlx_publish_str(sys, STR_LIT("vlx/dft_functional"), STR_LIT("DFT Functional"), md_vlx_dft_func_label(vlx));

	if (md_vlx_rsp_type(vlx) == MD_VLX_RSP_C6) {
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/c6"), STR_LIT("C6 Coefficient"), md_unit_none(), md_vlx_c6_value(vlx));
	}

	// ---- The QM ATOM DOMAIN. ----
	//
	// The atoms this calculation covered, in ITS order and at ITS geometry. That is not the
	// system's atom set: a calculation can cover part of a loaded system - a chromophore inside a
	// protein - so the two spaces differ in length AND in order, and qm/atom/system_index is the
	// only bridge between them.
	//
	// Its own prefix, and NOT atom/*, which is the SYSTEM's atom domain: a consumer walking atom/
	// has no idea quantum chemistry exists and would index a QM column by system atom - silently,
	// and wrongly, on exactly the subset case this domain exists for. atom_property_query already
	// filters by extent and component count, but that filter passes a scalar QM column whenever the
	// two atom counts happen to agree, which is not a distinction worth resting on.
	//
	// Not under basis/ either, though the shell list indexes this space: the basis is one thing
	// defined OVER these atoms, and so are the normal modes below. A nuclear coordinate is not a
	// property of a basis set, and a file can carry a geometry without carrying a basis at all.
	// basis/shell/atom_index is an index INTO qm/atom, which is the relationship stated plainly.
	{
		const size_t num_qm_atoms = md_vlx_number_of_atoms(vlx);
		const md_element_t* atomic_number = md_vlx_atomic_numbers(vlx);

		if (num_qm_atoms > 0 && atomic_number) {
			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_U8, .components = 1, .rank = 1, .shape = { (uint32_t)num_qm_atoms },
			};
			vlx_publish(sys, STR_LIT("qm/atom/atomic_number"), STR_LIT("Atomic Number"), md_unit_none(),
						format, atomic_number, num_qm_atoms * sizeof(md_element_t));
		}

		// Angstrom, matching the system's own coordinates rather than the bohr the evaluator works
		// in: a consumer comparing this geometry against md_system_state_t should not have to convert
		// first. This is the geometry the CALCULATION was run at, which is not necessarily where the
		// system's atoms are now - a trajectory frame or an optimisation step moves them.
		vlx_publish_vec3_series(sys, STR_LIT("qm/atom/coordinate"), STR_LIT("Coordinate"), angstrom,
								md_vlx_atom_coordinates(vlx), num_qm_atoms);
	}

	// ---- SCF: one value per iteration of the convergence history. ----
	const size_t num_iter = md_vlx_scf_history_size(vlx);
	vlx_publish_series(sys, STR_LIT("vlx/scf/history/energy"),        STR_LIT("Energy"),				hartree,		md_vlx_scf_history_energy(vlx),			num_iter);
	vlx_publish_series(sys, STR_LIT("vlx/scf/history/energy_diff"),   STR_LIT("Energy Difference"),		hartree,        md_vlx_scf_history_energy_diff(vlx),	num_iter);
	vlx_publish_series(sys, STR_LIT("vlx/scf/history/density_diff"),  STR_LIT("Density Difference"),	md_unit_none(), md_vlx_scf_history_density_diff(vlx),	num_iter);
	vlx_publish_series(sys, STR_LIT("vlx/scf/history/gradient_norm"), STR_LIT("Gradient Norm"),			md_unit_none(), md_vlx_scf_history_gradient_norm(vlx),	num_iter);
	vlx_publish_series(sys, STR_LIT("vlx/scf/history/max_gradient"),  STR_LIT("Max Gradient"),			md_unit_none(), md_vlx_scf_history_max_gradient(vlx),	num_iter);

	// ---- SCF: one value per molecular orbital, per spin. Sibling paths rather than one {2,M}
	// attribute, because a beta set is either present or absent and never an index to loop over.
	const size_t num_mos = md_vlx_scf_number_of_molecular_orbitals(vlx);
	const double* alpha_energy_data = md_vlx_scf_mo_energy(vlx, MD_VLX_SPIN_ALPHA);
	const double* beta_energy_data  = md_vlx_scf_mo_energy(vlx, MD_VLX_SPIN_BETA);
	const double* alpha_occ_data    = md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_ALPHA);
	const double* beta_occ_data     = md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_BETA);

	const md_attribute_id_t alpha_energy = vlx_publish_series(sys, STR_LIT("vlx/scf/orbital/alpha/energy"),     STR_LIT("Energy"),		hartree,        alpha_energy_data,	num_mos);
	const md_attribute_id_t alpha_occ    = vlx_publish_series(sys, STR_LIT("vlx/scf/orbital/alpha/occupation"), STR_LIT("Occupation"), md_unit_none(), alpha_occ_data,		num_mos);

	// A restricted calculation has one set of orbitals; beta is a second name for it, not a second
	// copy. The open shell case shares the energies and not the occupations, and falls out of the
	// same pointer test without being special cased.
	const md_attribute_id_t beta_energy  = vlx_publish_or_alias(sys, alpha_energy, STR_LIT("vlx/scf/orbital/beta/energy"),      STR_LIT("Energy"),		hartree,        alpha_energy_data, beta_energy_data, num_mos);
	const md_attribute_id_t beta_occ     = vlx_publish_or_alias(sys, alpha_occ,    STR_LIT("vlx/scf/orbital/beta/occupation"),  STR_LIT("Occupation"), md_unit_none(), alpha_occ_data,    beta_occ_data,    num_mos);

	// An orbital's coefficients and density are already published under the neutral orbital/ tree,
	// because md_gto_basis_t is mdlib's representation and not this program's. Its energy and
	// occupation are the same objects described the same way by every QM code, so leaving them
	// only under vlx/ splits one set of orbitals across two namespaces.
	//
	// Aliasing rather than moving them: a consumer already reading the vlx/ path keeps working,
	// and one written against the neutral name works too. That is the whole reason a second name
	// is cheaper than a rename.
	vlx_alias(sys, alpha_energy, STR_LIT("orbital/alpha/energy"));
	vlx_alias(sys, alpha_occ,    STR_LIT("orbital/alpha/occupation"));
	vlx_alias(sys, beta_energy,  STR_LIT("orbital/beta/energy"));
	vlx_alias(sys, beta_occ,     STR_LIT("orbital/beta/occupation"));

	// ---- The GTO basis and the MO coefficients: everything an evaluator needs to turn an orbital
	// into samples on a grid, with no reader in the loop.
	//
	// These land on FORMAT NEUTRAL paths, unlike the vlx/ tree above, because they are not this
	// program's output: md_gto_basis_t is mdlib's own normalised Cartesian representation, and
	// md_vlx_gto_basis_extract is the conversion into it. Another QM reader fills the same paths
	// with the same meaning, which is the test a path has to pass to lose its prefix.
	//
	// The shell list is published as four columns rather than as md_gto_shell_t records. A record
	// in the table would be a struct layout contract, additive changes to it would silently
	// invalidate anything stored, and nothing here would be self describing. The price is one
	// interleave in md_gto_basis_extract_attributes, paid when a consumer builds a basis rather
	// than per evaluation.
	{
		md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
		md_gto_basis_t basis = {0};

		if (md_vlx_gto_basis_extract(&basis, vlx, md_temp_allocator(temp))) {
			const size_t num_shells     = basis.num_shells;
			const size_t num_primitives = basis.num_primitives;
			const size_t shell_stride   = sizeof(md_gto_shell_t);

			vlx_publish_column(sys, STR_LIT("basis/shell/atom_index"),       STR_LIT("Atom Index"),       md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->atom_idx,         shell_stride, num_shells);
			vlx_publish_column(sys, STR_LIT("basis/shell/primitive_offset"), STR_LIT("Primitive Offset"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->primitive_offset, shell_stride, num_shells);
			vlx_publish_column(sys, STR_LIT("basis/shell/primitive_count"),  STR_LIT("Primitive Count"),  md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->num_primitives,   shell_stride, num_shells);
			vlx_publish_column(sys, STR_LIT("basis/shell/angular_momentum"), STR_LIT("Angular Momentum"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->l,                shell_stride, num_shells);

			// Exponents are bohr^-2 and the contraction coefficients carry the shell's radial
			// normalisation; the per monomial factor is applied at evaluation. See the AO
			// CONVENTION block in md_gto.h - these values only mean anything against it.
			const md_unit_t inv_bohr_sq = md_unit_pow(md_unit_bohr_radius(), -2);
			vlx_publish_column(sys, STR_LIT("basis/primitive/exponent"),    STR_LIT("Exponent"),	inv_bohr_sq,    MD_ATTRIBUTE_TYPE_F32, basis.alpha, sizeof(float), num_primitives);
			vlx_publish_column(sys, STR_LIT("basis/primitive/coefficient"), STR_LIT("Coefficient"), md_unit_none(), MD_ATTRIBUTE_TYPE_F32, basis.coeff, sizeof(float), num_primitives);

			// The MO matrix is stored [M][A] contiguously and is already Cartesian - the spherical
			// to Cartesian conversion happens once at parse time - so this is a straight copy into
			// a rank 2 attribute, and a consumer takes one orbital with a slice extract.
			//
			// f64 and not f32: md_gto takes AO coefficients as double to keep the QM code's
			// precision at the boundary, and there is no point publishing them already narrowed.
			const size_t num_ao = md_vlx_scf_number_of_atomic_orbitals(vlx);
			if (num_ao > 0 && num_mos > 0) {
				if (num_ao != md_gto_basis_num_ao(&basis)) {
					MD_LOG_ERROR("MO coefficients span %zu AOs but the basis has %zu; not publishing them", num_ao, md_gto_basis_num_ao(&basis));
				} else {
					md_attribute_format_t format = {
						.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
						.shape = { (uint32_t)num_mos, (uint32_t)num_ao },
					};
					const size_t byte_size = num_mos * num_ao * sizeof(double);

					// Row zero is the base of the whole matrix.
					const double* alpha_coeff = md_vlx_scf_mo_coefficients(vlx, 0, MD_VLX_SPIN_ALPHA);
					const double* beta_coeff  = md_vlx_scf_mo_coefficients(vlx, 0, MD_VLX_SPIN_BETA);

					md_attribute_id_t alpha_coeff_id = MD_ATTRIBUTE_INVALID;
					if (alpha_coeff) {
						alpha_coeff_id = vlx_publish(sys, STR_LIT("orbital/alpha/coefficient"),	STR_LIT("Alpha Coefficient"),	md_unit_none(), format, alpha_coeff, byte_size);
					}
					// Beta previously went unpublished whenever it shared alpha's buffer, so a
					// restricted calculation simply had no orbital/beta/coefficient for anyone to
					// read. A second name costs nothing and makes it present and correct.
					if (beta_coeff) {
						if (beta_coeff == alpha_coeff) {
							vlx_alias(sys, alpha_coeff_id, STR_LIT("orbital/beta/coefficient"));
						} else {
							vlx_publish(sys, STR_LIT("orbital/beta/coefficient"),	STR_LIT("Beta Coefficient"),	md_unit_none(), format, beta_coeff, byte_size);
						}
					}

					// The AO overlap S, {A,A} and symmetric. It belongs to the BASIS and not to a spin
					// channel: both channels and every density in this table are expressed against this
					// one metric, so basis/ and not orbital/.
					//
					// It is in the same AO order and convention as the coefficients above - converted by
					// vlx_cart_convert_square and permuted by ao_permute_square, exactly as they are - so
					// the two can be used together without further ceremony. Read the AO CONVENTION block
					// in md_gto.h first: this is the CARTESIAN overlap, and the Cartesian embedding of a
					// spherical basis is rank deficient, so for any file that stored spherical data this
					// matrix is SINGULAR. That is fine for what it is wanted for - Mulliken partitioning
					// and tr(DS) - and fatal for anything that inverts or factorises it.
					const double* overlap = md_vlx_scf_overlap_matrix_data(vlx);
					if (overlap && md_vlx_scf_overlap_matrix_size(vlx) == num_ao) {
						md_attribute_format_t overlap_format = {
							.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
							.shape = { (uint32_t)num_ao, (uint32_t)num_ao },
						};
						vlx_publish(sys, STR_LIT("basis/overlap"), STR_LIT("AO Overlap"), md_unit_none(),
									overlap_format, overlap, num_ao * num_ao * sizeof(double));
					}

					// Ground state densities: computed on demand from the coefficients and occupations
					// just published above rather than kept as a second resident [A][A] copy. See
					// vlx_build_occupation_density_matrix() for why this is exact, not approximate.
					md_attribute_format_t density_format = {
						.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
						.shape = { (uint32_t)num_ao, (uint32_t)num_ao },
					};
					md_attribute_id_t alpha_density_id = MD_ATTRIBUTE_INVALID;
					if (alpha_coeff) {
						md_attribute_virtual_t alpha_density_virt = { .provider = vlx_scf_alpha_density_provider, .user_data = sys };
						alpha_density_id = vlx_publish_virtual(sys, STR_LIT("orbital/alpha/density"), STR_LIT("Alpha Density"), md_unit_none(), density_format, &alpha_density_virt);
					}
					if (beta_coeff) {
						// A density is built from coefficients AND occupations, so it is the same
						// density only when BOTH are shared. Restricted open shell shares the
						// orbitals and not the occupations, so beta gets its own provider there -
						// which now works, because the coefficients it reads exist as an alias.
						const bool same_density = (beta_coeff == alpha_coeff) &&
							(md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_BETA) == md_vlx_scf_mo_occupancy(vlx, MD_VLX_SPIN_ALPHA));

						if (same_density) {
							vlx_alias(sys, alpha_density_id, STR_LIT("orbital/beta/density"));
						} else {
							md_attribute_virtual_t beta_density_virt = { .provider = vlx_scf_beta_density_provider, .user_data = sys };
							vlx_publish_virtual(sys, STR_LIT("orbital/beta/density"), STR_LIT("Beta Density"), md_unit_none(), density_format, &beta_density_virt);
						}

						// The two the UI offers besides the spins. Published rather than combined at
						// the point of use, so that "which densities does this system have" is one
						// question asked of the table and not four cases in a consumer.
						md_attribute_virtual_t total_density_virt = { .provider = vlx_scf_total_density_provider,      .user_data = sys };
						md_attribute_virtual_t diff_density_virt  = { .provider = vlx_scf_difference_density_provider, .user_data = sys };
						vlx_publish_virtual(sys, STR_LIT("orbital/total/density"),      STR_LIT("Total Density"),       md_unit_none(), density_format, &total_density_virt);
						vlx_publish_virtual(sys, STR_LIT("orbital/difference/density"), STR_LIT("Spin Difference Density"), md_unit_none(), density_format, &diff_density_virt);
					}
				}
			}
		}

		md_temp_end(temp);
	}

	// ---- Density properties: AO basis {A,A} matrices carried through from the file as they were
	// found. Unlike the SCF densities above these are not derived from anything else the table
	// holds - there is nothing to compute them from - so they are resident, one path per property.
	//
	// The path is built from the DATASET NAME, not from the label and not from the index. An index
	// silently re-points at a different property whenever the set changes across a reload, and a
	// label is display text that two datasets are free to share.
	for (size_t i = 0; i < md_array_size(vlx->density_properties); ++i) {
		const md_vlx_density_property_t* prop = &vlx->density_properties[i];
		if (!prop || !prop->data || prop->dim[0] == 0 || prop->dim[1] == 0) {
			continue;
		}

		str_t name = str_empty(prop->name) ? prop->label : prop->name;
		if (str_empty(name)) {
			continue;
		}

		char path_buf[256];
		str_t path = vlx_attribute_path(path_buf, sizeof(path_buf), STR_LIT("vlx/density_property"), name);
		if (str_empty(path)) {
			continue;
		}

		// The label is what the file called it for a human, the path is its identity. When they are
		// the same word there is nothing for the label to add, and an absent one is a valid state.
		str_t label = str_eq(prop->label, name) ? (str_t){0} : prop->label;

		md_attribute_format_t format = {
			.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
			.shape = { (uint32_t)prop->dim[0], (uint32_t)prop->dim[1] },
		};
		vlx_publish(sys, path, label, md_unit_none(), format, prop->data, prop->dim[0] * prop->dim[1] * sizeof(double));
	}

	// ---- RSP: one value per excited state. ----
	const size_t num_states = md_vlx_rsp_number_of_excited_states(vlx);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/oscillator_strength"), STR_LIT("Oscillator Strength"), md_unit_none(), md_vlx_rsp_oscillator_strengths(vlx), num_states);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/rotatory_strength"),   STR_LIT("Rotatory Strength"),   md_unit_none(), md_vlx_rsp_rotatory_strengths(vlx),   num_states);

	// ---- RSP: the raw solution vectors and the occupied/virtual split they are indexed by. These
	// are not meant for direct consumption - a consumer wants the reconstructed density, not the
	// eigenvector it came from - they exist purely so vlx_transition_density_provide() below can
	// rebuild a transition density with no 'vlx' object in reach, which is the whole point of
	// storing them as attributes rather than reaching back into 'vlx' from the provider.
	if (vlx->rsp.solution_matrix.data && num_states > 0) {
		vlx_publish_matrix(sys, STR_LIT("vlx/rsp/solution_matrix"), (str_t){0}, md_unit_none(),
			vlx->rsp.solution_matrix.data, vlx->rsp.solution_matrix.size[0], vlx->rsp.solution_matrix.size[1]);
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/num_core"),    (str_t){0}, md_unit_none(), (double)vlx->rsp.num_core);
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/num_valence"), (str_t){0}, md_unit_none(), (double)vlx->rsp.num_valence);
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/num_virtual"), (str_t){0}, md_unit_none(), (double)vlx->rsp.num_virtual);

		// One virtual attribute per sibling density, each computed on demand from the solution
		// vectors and MO coefficients above. rank {S,A,A}: slice by state for one density, or take
		// the whole thing and pay for reconstructing every state - see the caveat on
		// vlx_transition_density_provide() about that cost.
		const size_t num_ao = md_vlx_scf_number_of_atomic_orbitals(vlx);
		if (num_ao > 0) {
			md_attribute_format_t density_format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 3,
				.shape = { (uint32_t)num_states, (uint32_t)num_ao, (uint32_t)num_ao },
			};
			md_attribute_virtual_t attach_virt = { .provider = vlx_transition_density_attachment_provider, .user_data = sys };
			md_attribute_virtual_t detach_virt = { .provider = vlx_transition_density_detachment_provider, .user_data = sys };
			md_attribute_virtual_t diff_virt   = { .provider = vlx_transition_density_difference_provider, .user_data = sys };
			vlx_publish_virtual(sys, STR_LIT("vlx/rsp/transition_density/attachment"), STR_LIT("Attachment Density"), md_unit_none(), density_format, &attach_virt);
			vlx_publish_virtual(sys, STR_LIT("vlx/rsp/transition_density/detachment"), STR_LIT("Detachment Density"), md_unit_none(), density_format, &detach_virt);
			vlx_publish_virtual(sys, STR_LIT("vlx/rsp/transition_density/difference"), STR_LIT("Difference Density"), md_unit_none(), density_format, &diff_virt);
		}
	}

	// ---- RSP: the frequency axis, and the quantities sampled over it. What the axis MEANS depends
	// on md_vlx_rsp_type: excitation energies for a linear response, a sampled grid for a complex
	// polarisation propagator run. Same numbers either way, which is why it is one path.
	const size_t num_freqs = md_vlx_rsp_number_of_frequencies(vlx);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/frequency"),            STR_LIT("Response Frequency"),       hartree,        md_vlx_rsp_frequencies(vlx),       num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/cpp/sigma"),            STR_LIT("Absorption Cross Section"), md_unit_none(), md_vlx_rsp_sigma(vlx),             num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/cpp/delta_epsilon"),    STR_LIT(u8"Δε"),					  md_unit_none(), md_vlx_rsp_delta_epsilons(vlx),    num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/cpp/optical_rotation"), STR_LIT("Optical Rotation"),         md_unit_none(), md_vlx_rsp_optical_rotations(vlx), num_freqs);

	// ---- RSP two photon absorption. ----
	vlx_publish_series(sys, STR_LIT("vlx/rsp/tpa/cross_section"),   STR_LIT("Cross Section"),           md_unit_none(), md_vlx_rsp_tpa_cross_sections(vlx),  num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/tpa/gamma_re"),        STR_LIT(u8"γ (Re)"),				md_unit_none(), md_vlx_rsp_tpa_gamma_re(vlx),        num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/tpa/gamma_im"),        STR_LIT(u8"γ (Im)"),				md_unit_none(), md_vlx_rsp_tpa_gamma_im(vlx),        num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/tpa/linear"),          STR_LIT("Linear Polarisation"),     md_unit_none(), md_vlx_rsp_tpa_trans_linear(vlx),    num_freqs);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/tpa/circular"),        STR_LIT("Circular Polarisation"),   md_unit_none(), md_vlx_rsp_tpa_trans_circular(vlx),  num_freqs);

	// ---- RSP RIXS. The 2D quantities are {F,P}: final state outermost, photon energy innermost,
	// exactly as md_vlx.h documents the buffers, so the shape is a transcription and not a choice.
	const size_t num_photon = md_vlx_rsp_rixs_number_of_photon_energies(vlx);
	const size_t num_final  = md_vlx_rsp_rixs_number_of_final_states(vlx);
	const size_t num_core   = md_vlx_rsp_rixs_number_of_core_states(vlx);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/rixs/photon_energy"),           STR_LIT("Photon Energy"),				hartree,		md_vlx_rsp_rixs_photon_energies(vlx),			num_photon);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/rixs/elastic_cross_section"),   STR_LIT("Elastic Cross Section"),		md_unit_none(),	md_vlx_rsp_rixs_elastic_cross_sections(vlx),	num_photon);
	vlx_publish_matrix(sys, STR_LIT("vlx/rsp/rixs/cross_section"),           STR_LIT("Cross Section"),				md_unit_none(), md_vlx_rsp_rixs_cross_sections(vlx),			num_final, num_photon);
	vlx_publish_matrix(sys, STR_LIT("vlx/rsp/rixs/energy_loss"),             STR_LIT("Energy Loss"),				hartree,        md_vlx_rsp_rixs_energy_losses(vlx),				num_final, num_photon);
	vlx_publish_matrix(sys, STR_LIT("vlx/rsp/rixs/emission_energy"),         STR_LIT("Emission Energy"),			hartree,        md_vlx_rsp_rixs_emission_energies(vlx),			num_final, num_photon);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/rixs/core_energy"),             STR_LIT("Core Energy"),				hartree,        md_vlx_rsp_rixs_core_eigenvalues(vlx),			num_core);
	vlx_publish_series(sys, STR_LIT("vlx/rsp/rixs/core_oscillator_strength"),STR_LIT("Core Oscillator Strength"),	md_unit_none(), md_vlx_rsp_rixs_core_osc_strengths(vlx),		num_core);
	if (num_core > 0) {
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/rixs/gamma_fwhm"), STR_LIT("Core-hole Lifetime Broadening"), md_unit_electronvolt(), md_vlx_rsp_rixs_gamma_fwhm_ev(vlx));
	}

	// ---- XPS: one entry per computed core-hole state. The file's record has six fields of four
	// different types, which is six sibling paths over one {C} index space rather than one
	// attribute - a value has one type. Note what this buys: a consumer plotting ionization energy
	// against contribution now hands ImPlot two CONTIGUOUS arrays instead of one base pointer and a
	// struct stride.
	//
	// The per element grouping is not published. It is derived: entries are laid out as contiguous
	// runs of equal element, so a consumer wanting one element's states scans vlx/xps/element for
	// its run. Publishing the runs as well would be two representations of one fact, with nothing
	// keeping them in agreement.
	if (md_vlx_has_xps(vlx)) {
		const size_t              num_xps = md_vlx_xps_count(vlx);
		const md_vlx_xps_entry_t* xps     = md_vlx_xps_entries(vlx);

		if (xps && num_xps > 0) {
			// The bool field is copied a byte at a time as U8; a wider bool would take the wrong
			// byte on a big endian target, so it is worth failing the build rather than the render.
			STATIC_ASSERT(sizeof(bool) == 1, "XPS is_delocalized is published as a single byte");

			const size_t stride = sizeof(md_vlx_xps_entry_t);
			vlx_publish_column(sys, STR_LIT("vlx/xps/ionization_energy"), STR_LIT("Ionization Energy"),		md_unit_electronvolt(), MD_ATTRIBUTE_TYPE_F64, &xps->ionization_energy, stride, num_xps);
			vlx_publish_column(sys, STR_LIT("vlx/xps/contribution"),      STR_LIT("Core MO Contribution"),  md_unit_none(),         MD_ATTRIBUTE_TYPE_F64, &xps->contribution,      stride, num_xps);
			vlx_publish_column(sys, STR_LIT("vlx/xps/atom_index"),        STR_LIT("Atom Index"),            md_unit_none(),         MD_ATTRIBUTE_TYPE_I32, &xps->atom_index,        stride, num_xps);
			vlx_publish_column(sys, STR_LIT("vlx/xps/mo_index"),          STR_LIT("MO Index"),              md_unit_none(),         MD_ATTRIBUTE_TYPE_I32, &xps->mo_index,          stride, num_xps);
			vlx_publish_column(sys, STR_LIT("vlx/xps/element"),           STR_LIT("Atomic Number"),         md_unit_none(),         MD_ATTRIBUTE_TYPE_U8,  &xps->element,           stride, num_xps);
			vlx_publish_column(sys, STR_LIT("vlx/xps/is_delocalized"),    STR_LIT("Is Delocalized"),		md_unit_none(),			MD_ATTRIBUTE_TYPE_U8,  &xps->is_delocalized,	stride, num_xps);
		}
	}

	// ---- NTO lambdas: the weights of the natural transition orbital pairs, one row per excited
	// state. The rows are RAGGED - a state has as many pairs as it has - so they are padded to the
	// widest row and the shape is {S,Lmax}. Zero is the honest pad here rather than a sentinel: a
	// lambda IS a weight, an absent pair carries none, and the consumer which already stops at a
	// 1e-3 cutoff stops in exactly the same place. Anything ragged enough that zero would be a real
	// value does not belong in a rectangular attribute at all.
	if (md_vlx_rsp_has_nto(vlx) && num_states > 0) {
		double row[MD_VLX_NTO_MAX_LAMBDAS];
		size_t max_lambdas = 0;
		for (size_t s = 0; s < num_states; ++s) {
			const size_t n = md_vlx_rsp_nto_lambdas_extract(row, vlx, s, ARRAY_SIZE(row));
			max_lambdas = MAX(max_lambdas, n);
		}

		if (max_lambdas > 0) {
			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
				.shape = { (uint32_t)num_states, (uint32_t)max_lambdas },
			};
			md_attribute_id_t id = vlx_publish(sys, STR_LIT("vlx/rsp/nto/lambda"), STR_LIT("NTO Weight"), md_unit_none(), format, NULL, 0);
			double* dst = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;

			if (dst) {
				// Created zeroed, so a short row needs nothing written past its own length.
				for (size_t s = 0; s < num_states; ++s) {
					const size_t n = md_vlx_rsp_nto_lambdas_extract(row, vlx, s, ARRAY_SIZE(row));
					MEMCPY(dst + s * max_lambdas, row, MIN(n, max_lambdas) * sizeof(double));
				}
			} else if (id != MD_ATTRIBUTE_INVALID) {
				md_attributes_remove(&sys->attributes, id);
			}
		}

		// The NTO coefficient vectors themselves, {S,Lmax,A}, padded on the lambda axis exactly as
		// the weights above so the two index the same space and a consumer moves one slider.
		//
		// RESIDENT rather than computed on demand, unlike the transition densities. The reason is
		// size, not principle: S x Lmax x A doubles is under a megabyte for a realistic file, while
		// a transition density is A^2 PER STATE. And it is what keeps this self contained - the NTO
		// math reads the vlx object, so a provider would have to close over the reader, which is
		// the one thing the port exists to avoid. Paying it once at load buys that outright.
		const size_t num_ao_nto = md_vlx_scf_number_of_atomic_orbitals(vlx);
		if (max_lambdas > 0 && num_ao_nto > 0) {
			const md_vlx_nto_type_t types[2] = { MD_VLX_NTO_PARTICLE, MD_VLX_NTO_HOLE };
			str_t paths[2] = { STR_LIT("vlx/rsp/nto/particle/coefficient"), STR_LIT("vlx/rsp/nto/hole/coefficient") };
			str_t labels[2] = { STR_LIT("Particle"), STR_LIT("Hole") };

			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 3,
				.shape = { (uint32_t)num_states, (uint32_t)max_lambdas, (uint32_t)num_ao_nto },
			};

			for (int t = 0; t < 2; ++t) {
				md_attribute_id_t id = vlx_publish(sys, paths[t], labels[t], md_unit_none(), format, NULL, 0);
				double* dst = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
				if (!dst) {
					if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, id);
					continue;
				}

				// One state at a time, straight into its own plane of the table's storage. The
				// extract writes [lambda_count][num_ao] which is exactly the plane's layout, and a
				// state with fewer pairs leaves the rest of its plane at the zero it was created
				// with.
				const size_t plane = max_lambdas * num_ao_nto;
				for (size_t s = 0; s < num_states; ++s) {
					md_vlx_rsp_nto_coefficients_extract(dst + s * plane, NULL, vlx, s, types[t], max_lambdas);
				}
			}
		}
	}

	// ---- VIB: one value per normal mode. ----
	const size_t num_modes = md_vlx_vib_number_of_normal_modes(vlx);
	vlx_publish_series(sys, STR_LIT("vlx/vib/frequency"),			STR_LIT("Frequency"),			wavenumber,		md_vlx_vib_frequencies(vlx),			num_modes);
	vlx_publish_series(sys, STR_LIT("vlx/vib/ir_intensity"),		STR_LIT("IR Intensity"),		km_per_mol,     md_vlx_vib_ir_intensities(vlx),			num_modes);
	vlx_publish_series(sys, STR_LIT("vlx/vib/reduced_mass"),		STR_LIT("Reduced Mass"),		amu,            md_vlx_vib_reduced_masses(vlx),			num_modes);
	vlx_publish_series(sys, STR_LIT("vlx/vib/force_constant"),		STR_LIT("Force Constant"),		md_unit_none(), md_vlx_vib_force_constants(vlx),		num_modes);
	vlx_publish_series(sys, STR_LIT("vlx/vib/external_frequency"),	STR_LIT("External Frequency"),	hartree,		md_vlx_vib_external_frequencies(vlx),	md_vlx_vib_number_of_external_frequencies(vlx));

	// The displacements are per atom, so the atom axis is the last index axis and the mode axis
	// leads: one mode's displacements are contiguous. This is the {M,N} case the ATTRIBUTES note in
	// md_system.h uses as its example, and it is why an atom axis is not always shape[0].
	vlx_publish_vec3_rows(sys, STR_LIT("qm/atom/normal_mode"), STR_LIT("Normal Mode"), md_unit_none(), vlx, num_modes, md_vlx_number_of_atoms(vlx), md_vlx_vib_normal_mode);

	// ---- OPT: one value per optimisation step, and the geometry at each one. ----
	const size_t num_steps = md_vlx_opt_number_of_steps(vlx);
	vlx_publish_series(sys,		STR_LIT("vlx/opt/energy"),		STR_LIT("Energy"),		hartree, md_vlx_opt_energies(vlx), num_steps);
	vlx_publish_vec3_rows(sys,	STR_LIT("vlx/opt/coordinate"),	STR_LIT("Coordinate"),	angstrom, vlx, num_steps, md_vlx_number_of_atoms(vlx), md_vlx_opt_coordinates);

	// ---- Dipole groups. Every vector gets an origin beside it in the SAME shape, because a dipole
	// has no index space of its own to anchor it - see the ANCHORING note in md_system.h. Publishing
	// a vector without its origin would leave a group nothing can draw, so all of this is skipped
	// together when the centre of charge cannot be computed.
	dvec3_t origin = {0, 0, 0};
	if (vlx_centre_of_charge(&origin, vlx)) {
		const dvec3_t ground_state = md_vlx_scf_ground_state_dipole_moment(vlx);
		vlx_publish_vec3_series(sys, STR_LIT("dipole/ground_state/vector"), STR_LIT("Ground State"), e_bohr, &ground_state, 1);
		vlx_publish_origin(sys, STR_LIT("dipole/ground_state/origin"), origin);

		// One per excited state. The magnetic and velocity forms are different quantities in
		// different units, which is exactly why unit sits on the attribute and not on the group.
		if (num_states > 0) {
			const dvec3_t* electric = md_vlx_rsp_electric_transition_dipole_moments(vlx);
			const dvec3_t* magnetic = md_vlx_rsp_magnetic_transition_dipole_moments(vlx);
			const dvec3_t* velocity = md_vlx_rsp_velocity_transition_dipole_moments(vlx);

			if (electric) {
				vlx_publish_vec3_series(sys, STR_LIT("dipole/electric_transition/vector"), STR_LIT("Electric Transition"), e_bohr, electric, num_states);
				vlx_publish_origin(sys, STR_LIT("dipole/electric_transition/origin"), origin);
			}
			if (magnetic) {
				vlx_publish_vec3_series(sys, STR_LIT("dipole/magnetic_transition/vector"), STR_LIT("Magnetic Transition"), md_unit_none(), magnetic, num_states);
				vlx_publish_origin(sys, STR_LIT("dipole/magnetic_transition/origin"), origin);
			}
			if (velocity) {
				vlx_publish_vec3_series(sys, STR_LIT("dipole/velocity_transition/vector"), STR_LIT("Velocity Transition"), md_unit_none(), velocity, num_states);
				vlx_publish_origin(sys, STR_LIT("dipole/velocity_transition/origin"), origin);
			}
		}
	}
}

bool md_vlx_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_vlx_t* vlx) {
	ASSERT(sys);
	ASSERT(state);
	ASSERT(vlx);
	
	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("The veloxchem object contains no atoms");
		return false;
	}

	if (!sys->alloc) {
		MD_LOG_ERROR("System allocator not set");
		return false;
    }

	if (!state || !state->alloc) {
		MD_LOG_ERROR("State allocator not set");
		return false;
	}

	md_system_reset(sys);
	md_system_state_init(state, vlx->number_of_atoms);

	size_t capacity = ROUND_UP(vlx->number_of_atoms, 16);

    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags,    capacity, sys->alloc);

	MEMSET(sys->atom.type_idx,  0, md_array_bytes(sys->atom.type_idx));
    MEMSET(sys->atom.flags,		0, md_array_bytes(sys->atom.flags));

    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

	for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
		state->x[i] = (float)vlx->atom_coordinates[i].x;
		state->y[i] = (float)vlx->atom_coordinates[i].y;
		state->z[i] = (float)vlx->atom_coordinates[i].z;
		
		md_atomic_number_t z = vlx->atomic_numbers[i];
		str_t sym  = md_atomic_number_symbol(z);
        float mass = md_atomic_number_mass(z);
		float radius = md_atomic_number_vdw_radius(z);
		uint32_t color = md_atomic_number_cpk_color(z);

		md_atom_type_idx_t type_idx = md_atom_type_find_or_add(&sys->atom.type, sym, z, mass, radius, color, 0, sys->alloc);
		sys->atom.type_idx[i] = type_idx;
	}

	sys->atom.count = vlx->number_of_atoms;
    state->num_atoms = sys->atom.count;

	// Publish the atomic properties into the system's attribute table at LOAD TIME. A consumer then
	// reads them straight out of sys->attributes and nothing has to keep the vlx object alive, or
	// ask it questions, to colour by a per atom quantity.
	//
	// The h5 dataspace and the attribute format are the same numbers in the same order:
	// h5_read_atomic_properties_in_group only accepts a dataset whose INNERMOST dimension is the
	// atom count, which is the attribute convention that the atom axis is the last index axis. So a
	// plain per atom property is rank 1 {N}, one with variants (excited states, spins, whatever the
	// file meant) is rank 2 {S,N}, and the values are scalars in both cases.
	for (size_t i = 0; i < md_array_size(vlx->atomic_properties); ++i) {
		const md_vlx_atomic_property_t* prop = &vlx->atomic_properties[i];

		str_t name = str_empty(prop->name) ? prop->label : prop->name;
		if (!prop->data || prop->num_dims < 1 || str_empty(name)) {
			continue;
		}

		md_attribute_format_t format = {
			.type       = MD_ATTRIBUTE_TYPE_F64,
			.components = 1,
			.rank       = (uint32_t)prop->num_dims,
		};
		size_t num_values = 1;
		for (int d = 0; d < prop->num_dims; ++d) {
			format.shape[d] = (uint32_t)prop->dim[d];
			num_values *= prop->dim[d];
		}

		char path_buf[256];
		str_t path = vlx_attribute_path(path_buf, sizeof(path_buf), STR_LIT("atom"), name);
		if (str_empty(path)) {
			continue;
		}

		// The label is what the file called it for a human, the path is its identity. When they are
		// the same word there is nothing for the label to add, and an absent one is a valid state.
		str_t label = str_eq(prop->label, name) ? (str_t){0} : prop->label;

		md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
			.path      = path,
			.format    = format,
			.unit      = md_unit_none(),
			.label     = label,
			.data      = prop->data,
			.byte_size = num_values * sizeof(double),
		});
	}

	return true;
}

// Which system atom each QM atom is, or nothing at all when the two spaces coincide.
//
// The file cannot decide this for itself: the same h5 carries a local-to-global map whether it is
// opened standalone - where the system IS the QM atoms and the map must NOT be applied - or against
// a larger system, where it must. What resolves it is WHICH ENTRY POINT WAS CALLED, which is why
// this takes the answer as an argument rather than trying to work it out.
//
// Publishing nothing is not the same as leaving it alone: a stale map from a previous load would
// send every evaluation to the wrong atoms, so the standalone case actively removes it.
static void vlx_publish_atom_system_index(md_system_t* sys, const md_vlx_t* vlx, bool supplemental) {
	ASSERT(sys);

	const str_t path = STR_LIT("qm/atom/system_index");
	// Declared before the if: this is C, not the C++ side of the tree.
	const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
	if (existing) {
		md_attributes_remove(&sys->attributes, existing->id);
	}
	if (!supplemental) {
		return;
	}

	const size_t num_qm_atoms = md_vlx_number_of_atoms(vlx);
	const int* local_to_global = md_vlx_local_to_global_atom_idx(vlx);
	if (num_qm_atoms == 0 || !local_to_global) {
		return;
	}

	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_U32, .components = 1, .rank = 1, .shape = { (uint32_t)num_qm_atoms },
	};
	md_attribute_id_t id = md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
		.path   = path,
			.format = format,
			.unit   = md_unit_none(),
			.label  = STR_LIT("System Atom Index"),
	});

	uint32_t* dst = (uint32_t*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_U32);
	if (!dst) {
		if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, id);
		return;
	}
	for (size_t i = 0; i < num_qm_atoms; ++i) {
		dst[i] = (uint32_t)local_to_global[i];
	}
}

bool md_vlx_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
	ASSERT(sys);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
	md_vlx_t* vlx = md_vlx_create(temp_arena);

	// VLX_FLAG_ALL rather than CORE, and the publish right here: the attribute table is populated
	// on the LOAD path, so a system carries its data the moment it is loaded and a consumer finds it
	// by asking the system. It used to be published by the veloxchem UI component, which meant the
	// data existed only because that component was compiled in and had parsed the same file a second
	// time - and meant no other reader could ever put comparable data in the same table.
	bool success = vlx_parse_file(vlx, filename, VLX_FLAG_ALL) && md_vlx_system_init_from_data(sys, state, vlx);
	if (success) {
		md_vlx_publish_attributes(sys, vlx);
		// Standalone: the system IS the QM atoms, so the map is cleared rather than written.
		vlx_publish_atom_system_index(sys, vlx, false);
	}

	md_temp_end(temp_scope);
	return success;
}

bool md_vlx_system_supplement_from_file(md_system_t* sys, str_t filename) {
	ASSERT(sys);

	if (!sys->attributes.alloc) {
		MD_LOG_ERROR("Cannot supplement a system which has not been initialised");
		return false;
	}

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
	md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
	md_vlx_t* vlx = md_vlx_create(temp_arena);

	// The atoms and the state are deliberately left alone - they belong to whatever loaded the
	// system first, and this file only adds to its table. Whether the file actually belongs to this
	// system is md_vlx_system_is_file_supplemental's question and the caller has already asked it.
	bool success = vlx_parse_file(vlx, filename, VLX_FLAG_ALL);
	if (success) {
		md_vlx_publish_attributes(sys, vlx);
		vlx_publish_atom_system_index(sys, vlx, true);
	}

	md_temp_end(temp_scope);
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

	h5_error_scope_t h5_err_scope = h5_error_scope_begin();

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		h5_error_scope_end(h5_err_scope);
		return false;
	}

	bool result = false;
	
    if (h5_check_dataset_exists(file_id, "qm_atom_indices") &&
		h5_check_dataset_exists(file_id, "nuclear_charges"))
	{
		md_temp_scope_t temp_scope = md_temp_begin();
		md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
		// Derive atomic numbers from nuclear charges and verify with qm_atom_indices that they are consistent with the loaded system.
        // This is a heuristic check, but it should be sufficient to determine if the file contains data that can supplement the existing system.
		
		size_t vlx_num_atoms = 0;
		if (!h5_read_scalar(&vlx_num_atoms, file_id, H5T_NATIVE_UINT64, "number_of_atoms")) {
			goto temp_done;
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
	temp_done:
		md_temp_end(temp_scope);
	}

	H5Fclose(file_id);
	h5_error_scope_end(h5_err_scope);

	return result;
}

// Externally visible procedures

size_t md_vlx_number_of_atoms(const md_vlx_t* vlx) {
	if (vlx) return vlx->number_of_atoms;
	return 0;
}

size_t md_vlx_number_of_electrons(const md_vlx_t* vlx, md_vlx_spin_t spin) {
	if (vlx) {
		if (spin == MD_VLX_SPIN_ALPHA) {
			return vlx->number_of_alpha_electrons;
		} else if (spin == MD_VLX_SPIN_BETA) {
			return vlx->number_of_beta_electrons;
		}
	}
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
	return MD_VLX_SCF_UNKNOWN;
}

dvec3_t md_vlx_scf_ground_state_dipole_moment(const md_vlx_t* vlx) {
	if (vlx) return vlx->scf.ground_state_dipole_moment;
	return (dvec3_t){0};
}

size_t md_vlx_scf_homo_idx(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.homo_idx;
		} else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.beta.homo_idx;
		}
	}
	return 0;
}

size_t md_vlx_scf_lumo_idx(const md_vlx_t* vlx, md_vlx_spin_t type) {
	if (vlx) {
		if (type == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.lumo_idx;
		} else if (type == MD_VLX_SPIN_BETA) {
			return vlx->scf.beta.lumo_idx;
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

size_t md_vlx_scf_number_of_occupied(const md_vlx_t* vlx, md_vlx_spin_t spin) {
	if (vlx) {
		if (spin == MD_VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.lumo_idx;
		} else if (spin == MD_VLX_SPIN_BETA) {
			return vlx->scf.beta.lumo_idx;
		}
	}
	return 0;
}

size_t md_vlx_scf_number_of_virtual(const md_vlx_t* vlx, md_vlx_spin_t spin) {
	if (vlx) {
		return md_vlx_scf_number_of_molecular_orbitals(vlx) - md_vlx_scf_number_of_occupied(vlx, spin);
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
					.num_primitives   = (uint32_t)bf->count,
					.l                = (uint32_t)angl,
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
const double* md_vlx_scf_mo_coefficients(const md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t spin) {
	if (!vlx) return NULL;
	const md_vlx_orbital_t* orb = (spin == MD_VLX_SPIN_ALPHA) ? &vlx->scf.alpha :
								  (spin == MD_VLX_SPIN_BETA)  ? &vlx->scf.beta  : NULL;
	if (!orb || !orb->coefficients.data) return NULL;
	size_t num_mo = orb->coefficients.size[0];
	size_t num_ao = orb->coefficients.size[1];
	if (mo_idx >= num_mo || num_ao == 0) return NULL;
	return orb->coefficients.data + mo_idx * num_ao;
}

// Deprecated extraction wrappers kept for backward compatibility.
size_t md_vlx_scf_mo_coefficients_extract(double* out, const md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t spin) {
	const double* src = md_vlx_scf_mo_coefficients(vlx, mo_idx, spin);
	if (!src) return 0;
	size_t num_ao = number_of_ao_coefficients(&vlx->scf.alpha); 
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
