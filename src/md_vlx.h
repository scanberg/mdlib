#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <md_types.h>
#include <md_gto.h>
#include <core/md_str.h>

struct md_allocator_i;
struct md_system_t;
struct md_system_loader_i;

typedef struct md_vlx_t md_vlx_t;

typedef struct md_vlx_atomic_property_t {
	str_t    label;
	uint64_t key;
	size_t   dim[2];	// dim[0] = number of atoms, dim[1] = 1 for scalar properties, N for multi-dimensional properties.
	double*  data;
} md_vlx_atomic_property_t;

typedef struct md_vlx_density_property_t {
	str_t    label;
	uint64_t key;
	size_t 	 dim[2];
	double*  data;
} md_vlx_density_property_t;

// XPS: one entry per computed core-hole state.
// VeloxChem performs the calculation per element, removing one core electron from each atom of that
// element in turn and taking the total energy difference against the ground state (delta-SCF).
// Field order is deliberate: 'ionization_energy' and 'contribution' are both double and adjacent so
// they can be handed straight to ImPlot as x/y with stride = sizeof(md_vlx_xps_entry_t).
typedef struct md_vlx_xps_entry_t {
	double       ionization_energy;	// unit: eV
	double       contribution;		// Fraction of the core MO localized on 'atom_index'
	int32_t      atom_index;		// Index into the molecular structure
	int32_t      mo_index;			// Index into the MO (orbital) arrays
	md_element_t element;			// Atomic number, matches the owning group
	bool         is_delocalized;	// Core hole is spread over several symmetry equivalent atoms
} md_vlx_xps_entry_t;

// A group is the set of entries belonging to a single element. 'entries' points into the flat array
// returned by md_vlx_xps_entries() -- it is a view, not a separate allocation, and shares its lifetime.
typedef struct md_vlx_xps_group_t {
	md_element_t              element;
	size_t                    count;
	const md_vlx_xps_entry_t* entries;	// Contiguous run, sorted by ionization_energy
} md_vlx_xps_group_t;

typedef enum {
	MD_VLX_SPIN_ALPHA = 0,
	MD_VLX_SPIN_BETA  = 1,
} md_vlx_spin_t;

typedef enum {
	MD_VLX_NTO_PARTICLE = 0,
	MD_VLX_NTO_HOLE = 1,
} md_vlx_nto_type_t;

typedef enum {
	MD_VLX_TRANSITION_ATTACHMENT = 0,
	MD_VLX_TRANSITION_DETACHMENT = 1,
	MD_VLX_TRANSITION_DIFFERENCE = 2,
} md_vlx_transition_type_t;

typedef enum {
	MD_VLX_SCF_UNKNOWN = 0,
	MD_VLX_SCF_RESTRICTED,
	MD_VLX_SCF_RESTRICTED_OPENSHELL,
	MD_VLX_SCF_UNRESTRICTED,
} md_vlx_scf_type_t;

typedef enum {
	MD_VLX_RSP_UNKNOWN = 0,
    MD_VLX_RSP_LINEAR,			// Linear response (frequency-dependent polarizabilities)
    MD_VLX_RSP_CPP,				// Complex Polarization Propagator
    MD_VLX_RSP_C6,				// Homomolecular C_6 value (in a.u.)
    MD_VLX_RSP_TPA,		        // Two-photon Absorption (TPA) cross-sections
    MD_VLX_RSP_TPA_TRANSITION,  // Two-Photon Absorption transition properties
    MD_VLX_RSP_RIXS,			// Resonant Inelastic X-ray Scattering (RIXS) cross-sections
} md_vlx_rsp_type_t;

// PES (Potential Energy Surface) operations (Geometry Optimizations)
typedef enum {
    MD_VLX_OPT_UNKNOWN = 0,
    MD_VLX_OPT_GEOMETRY,			// Local minimum optimization (Ground state or excited state)
    MD_VLX_OPT_CONSTRAINED,			// Optimization with geometric constraints
    MD_VLX_OPT_TS,					// First-order saddle point optimization (Transition State)
    MD_VLX_OPT_IRC,					// Reaction path optimization (Intrinsic Reaction Coordinate)
	MD_VLX_OPT_COUNT,
} md_vlx_opt_type_t;

#ifdef __cplusplus
extern "C" {
#endif

// Basic operations

struct md_vlx_t* md_vlx_create(struct md_allocator_i* backing);
void md_vlx_reset(struct md_vlx_t* vlx);
void md_vlx_destroy(struct md_vlx_t* vlx);

// Will parse data from .h5 files or .out which are produced by veloxchem
bool md_vlx_parse_file(struct md_vlx_t* vlx, str_t filename);

size_t md_vlx_number_of_atoms(const struct md_vlx_t* vlx);
size_t md_vlx_number_of_electrons(const struct md_vlx_t* vlx, md_vlx_spin_t spin);

double md_vlx_molecular_charge(const struct md_vlx_t* vlx);
double md_vlx_nuclear_repulsion_energy(const struct md_vlx_t* vlx);
size_t md_vlx_spin_multiplicity(const struct md_vlx_t* vlx);

str_t  md_vlx_basis_set_ident(const struct md_vlx_t* vlx);
str_t  md_vlx_dft_func_label(const struct md_vlx_t* vlx);
str_t  md_vlx_potfile(const struct md_vlx_t* vlx);

const dvec3_t* md_vlx_atom_coordinates(const struct md_vlx_t* vlx);
const uint8_t* md_vlx_atomic_numbers(const struct md_vlx_t* vlx);

// Maps AO index to atom index, length is number of atomic orbitals
const int* md_vlx_ao_to_atom_idx(const struct md_vlx_t* vlx);

// Maps local atom indices to full system indices.
// If the system is *NOT* a subsystem, NULL will be returned, and the full system indices are the same as the atom indices.
// Otherwise it will have the length of number of atoms, and for each atom index it will give the corresponding full system index that maps to it.
const int* md_vlx_local_to_global_atom_idx(const struct md_vlx_t* vlx);

// SCF
md_vlx_scf_type_t md_vlx_scf_type(const struct md_vlx_t* vlx);

dvec3_t md_vlx_scf_ground_state_dipole_moment(const struct md_vlx_t* vlx);
size_t  md_vlx_scf_homo_idx(const struct md_vlx_t* vlx, md_vlx_spin_t spin);
size_t  md_vlx_scf_lumo_idx(const struct md_vlx_t* vlx, md_vlx_spin_t spin);

size_t  md_vlx_scf_number_of_atomic_orbitals   (const struct md_vlx_t* vlx);
size_t  md_vlx_scf_number_of_molecular_orbitals(const struct md_vlx_t* vlx);

size_t md_vlx_scf_number_of_occupied(const struct md_vlx_t* vlx, md_vlx_spin_t spin);
size_t md_vlx_scf_number_of_virtual (const struct md_vlx_t* vlx, md_vlx_spin_t spin);

const double* md_vlx_scf_mo_occupancy(const struct md_vlx_t* vlx, md_vlx_spin_t spin);
const double* md_vlx_scf_mo_energy(const struct md_vlx_t* vlx, md_vlx_spin_t spin);

// ---------------------------------------------------------------------------
// Basis-centric API  (primary interface — use these in new code)
// ---------------------------------------------------------------------------

// Build a static md_gto_basis_t from vlx data. Call once at load time.
// Shells are emitted in order: angular momentum → atom → contracted function.
// All AO-indexed matrices stored in vlx are already reordered to match this
// shell order at parse time, so MO vectors and density matrices can be used
// directly without any permutation.
bool md_vlx_gto_basis_extract(md_gto_basis_t* out, const struct md_vlx_t* vlx, struct md_allocator_i* alloc);

// Direct pointer to the AO coefficient vector for MO mo_idx (shell order).
// The matrix is stored [num_mo][num_ao] internally; this returns a pointer to
// row mo_idx, which is a contiguous array of num_ao doubles.
// Returns NULL on failure (out-of-range mo_idx, missing data, etc.).
const double* md_vlx_scf_mo_coefficients(const struct md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t spin);

// Extract NTO AO coefficient vectors for the first lambda_count lambdas of one state/type.
// out_coefficients must hold lambda_count * md_vlx_scf_number_of_atomic_orbitals() doubles
// and is written as [lambda_count][num_ao] in shell order. out_lambdas is optional and,
// when non-NULL, must hold lambda_count doubles. Returns the number of lambdas written.
size_t md_vlx_rsp_nto_coefficients_extract(double* out_coefficients, double* out_lambdas, const struct md_vlx_t* vlx, size_t state_idx, md_vlx_nto_type_t type, size_t lambda_count);

// The overlap matrix (S) is a square, symmetric matrix [N][N], this returns the length N
size_t  md_vlx_scf_overlap_matrix_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_overlap_matrix_data(const struct md_vlx_t* vlx);

// The density matrix is a square, symmetric matrix [N][N], this returns the length N
size_t md_vlx_scf_density_matrix_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_density_matrix_data(const struct md_vlx_t* vlx, md_vlx_spin_t spin);

// ---------------------------------------------------------------------------
// Deprecated extraction API — prefer the direct pointer functions above
// ---------------------------------------------------------------------------

// Copy one MO coefficient vector into caller-supplied storage.
// out must hold at least md_vlx_scf_number_of_atomic_orbitals() doubles.
// Returns the number of AOs on success, 0 on failure.
size_t md_vlx_scf_mo_coefficients_extract(double* out, const struct md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t spin);

// SCF History
size_t		  md_vlx_scf_history_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_energy(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_energy_diff(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_density_diff(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_gradient_norm(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_max_gradient(const struct md_vlx_t* vlx);

// RSP (Standard Linear Response, TDDFT)
size_t md_vlx_rsp_number_of_excited_states(const struct md_vlx_t* vlx);

// RSP arrays with length of num excited states
const dvec3_t* md_vlx_rsp_electric_transition_dipole_moments(const struct md_vlx_t* vlx);
const dvec3_t* md_vlx_rsp_magnetic_transition_dipole_moments(const struct md_vlx_t* vlx);
const dvec3_t* md_vlx_rsp_velocity_transition_dipole_moments(const struct md_vlx_t* vlx);

const double*  md_vlx_rsp_rotatory_strengths(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_oscillator_strengths(const struct md_vlx_t* vlx);

// NEW UNIFIED RSP API
md_vlx_rsp_type_t md_vlx_rsp_type(const struct md_vlx_t* vlx);

// Homomolecular C_6 value (in a.u.), only valid for md_vlx_rsp_type() == MD_VLX_RSP_C6
double md_vlx_c6_value(const struct md_vlx_t* vlx);

size_t md_vlx_rsp_number_of_frequencies(const struct md_vlx_t* vlx);

// Frequencies are in a.u. (Hartree), and length is given by md_vlx_rsp_number_of_frequencies()
const double* md_vlx_rsp_frequencies(const struct md_vlx_t* vlx);

// CPP Absorption output in a.u.
const double* md_vlx_rsp_sigma(const struct md_vlx_t* vlx);

// CPP ECD output in delta epsilon [L mol^-1 cm^-1]
const double* md_vlx_rsp_delta_epsilons(const struct md_vlx_t* vlx);

const double* md_vlx_rsp_optical_rotations(const struct md_vlx_t* vlx);

// RSP TPA Transition (Valid output for md_vlx_rsp_type() == MD_VLX_RSP_TPA_TRANSITION)
const double* md_vlx_rsp_tpa_trans_linear(const struct md_vlx_t* vlx);
const double* md_vlx_rsp_tpa_trans_circular(const struct md_vlx_t* vlx);

// RSP TPA (Valid output only for md_vlx_rsp_type() == MD_VLX_RSP_TPA)
const double* md_vlx_rsp_tpa_cross_sections(const struct md_vlx_t* vlx);
const double* md_vlx_rsp_tpa_gamma_re(const struct md_vlx_t* vlx);
const double* md_vlx_rsp_tpa_gamma_im(const struct md_vlx_t* vlx);

// RSP RIXS (Valid output only for md_vlx_rsp_type() == MD_VLX_RSP_RIXS)
//
// Three independent dimensions are involved:
//   P = number of incoming photon energies      (md_vlx_rsp_rixs_number_of_photon_energies)
//   F = number of final (valence) states        (md_vlx_rsp_rixs_number_of_final_states)
//   C = number of core-excited (intermediate) states (md_vlx_rsp_rixs_number_of_core_states)
//
// The 2D quantities are stored ROW-MAJOR as [F][P], i.e. element (f,p) is at [f * P + p].
// This matches the numpy layout produced by RixsDriver (shape (num_final_states, num_photon_energies)).
size_t md_vlx_rsp_rixs_number_of_photon_energies(const struct md_vlx_t* vlx);
size_t md_vlx_rsp_rixs_number_of_final_states(const struct md_vlx_t* vlx);
size_t md_vlx_rsp_rixs_number_of_core_states(const struct md_vlx_t* vlx);

// Length P, unit: a.u. (Hartree). Incoming photon energies (omega).
const double* md_vlx_rsp_rixs_photon_energies(const struct md_vlx_t* vlx);

// Length P, unit: a.u. Elastic (Rayleigh) scattering cross-section per photon energy.
const double* md_vlx_rsp_rixs_elastic_cross_sections(const struct md_vlx_t* vlx);

// Row-major [F][P], unit: a.u. Inelastic scattering cross-sections.
const double* md_vlx_rsp_rixs_cross_sections(const struct md_vlx_t* vlx);

// Row-major [F][P], unit: a.u. Energy loss (== valence eigenvalue, independent of photon energy).
const double* md_vlx_rsp_rixs_energy_losses(const struct md_vlx_t* vlx);

// Row-major [F][P], unit: a.u. Emission energy (omega - valence eigenvalue).
const double* md_vlx_rsp_rixs_emission_energies(const struct md_vlx_t* vlx);

// Length C, unit: a.u. (Hartree). Core-excitation energies, used for the XAS side panel.
const double* md_vlx_rsp_rixs_core_eigenvalues(const struct md_vlx_t* vlx);

// Length C, unitless. Oscillator strengths of the core excitations.
const double* md_vlx_rsp_rixs_core_osc_strengths(const struct md_vlx_t* vlx);

// Core-hole lifetime broadening, unit: eV (FWHM). Used to broaden the XAS panel.
double md_vlx_rsp_rixs_gamma_fwhm_ev(const struct md_vlx_t* vlx);

// RSP NTO api

bool md_vlx_rsp_has_nto(const struct md_vlx_t* vlx);
size_t md_vlx_rsp_nto_lambdas_extract(double* out_lambdas, const md_vlx_t* vlx, size_t state_idx, size_t lambda_count);

// Extract transition-density AO matrix [N][N] in shell order from response eigenvectors.
// out_matrix must hold md_vlx_rsp_transition_density_matrix_size(vlx, state_idx)^2 doubles.
// Attachment = T^T T, Detachment = T T^T, Difference = Attachment - Detachment, where T = Z - Y.
size_t md_vlx_rsp_transition_density_matrix_size(const struct md_vlx_t* vlx, size_t state_idx);
size_t md_vlx_rsp_transition_density_matrix_extract(double* out_matrix, const struct md_vlx_t* vlx, size_t state_idx, md_vlx_transition_type_t type);

// VIB
// Many of the fields within the VIB portion has a length given by the degrees of freedom (D)
size_t md_vlx_vib_number_of_normal_modes(const struct md_vlx_t* vlx);

// Returns arrays of length D (Number of normal modes)
const double* md_vlx_vib_ir_intensities(const struct md_vlx_t* vlx); 	// Unit: km/mol
const double* md_vlx_vib_frequencies(const struct md_vlx_t* vlx);		// Unit: cm^-1
const double* md_vlx_vib_reduced_masses(const struct md_vlx_t* vlx);	// Unit: amu
const double* md_vlx_vib_force_constants(const struct md_vlx_t* vlx);	// Unit: mdyne/Å (millidyne / Ångström)

// OPT (Potential Energy Surface) Optimizations
// @TODO: Add constraints (when they appear in h5)

// Returns the type of optimization performed (if any) in the VLX file.
md_vlx_opt_type_t md_vlx_opt_type(const struct md_vlx_t* vlx);

// Returns the electronic state index of the optimization (if any) in the VLX file.
// 0 == ground state, 1 == first excited state, etc.
// This is only used if md_vlx_opt_type == MD_VLX_OPT_GEOMETRY to indicate which state was optimized.
size_t md_vlx_opt_state_index(const struct md_vlx_t* vlx); 

// Returns the index of the optimization that is used as the reference point to perform optimization forwards and backwards along the gradient.
// This is only used if md_vlx_opt_type == MD_VLX_OPT_IRC.
size_t md_vlx_opt_irc_ts_index(const struct md_vlx_t* vlx);

// Returns the number of optimization steps (length of the arrays below)
size_t md_vlx_opt_number_of_steps(const struct md_vlx_t* vlx);

// Returns atom coordinates for a given optimization step
const dvec3_t* md_vlx_opt_coordinates(const struct md_vlx_t* vlx, size_t opt_idx);
const double*  md_vlx_opt_energies(const struct md_vlx_t* vlx);

// The normal mode for each atom
// Returns array of length N (Number of atoms)
const dvec3_t* md_vlx_vib_normal_mode(const struct md_vlx_t* vlx, size_t idx);

size_t md_vlx_vib_number_of_external_frequencies(const struct md_vlx_t* vlx);

const double* md_vlx_vib_external_frequencies(const struct md_vlx_t* vlx);	// Unit: a.u. (Hartree)

// Returns array of length D (Number of normal modes)
const double* md_vlx_vib_raman_activities(const struct md_vlx_t* vlx, size_t idx);

// XPS
// NOTE: XPS is a delta-SCF property, not a response property, so it is a top level data block and may
// coexist with any md_vlx_rsp_type_t (or with none at all).

bool md_vlx_has_xps(const struct md_vlx_t* vlx);

// Flat view over every entry across all elements, sorted by (element, ionization_energy).
// Use this for the summed spectrum and for atom_index hit tests.
size_t                    md_vlx_xps_count(const struct md_vlx_t* vlx);
const md_vlx_xps_entry_t* md_vlx_xps_entries(const struct md_vlx_t* vlx);

// Per element views into the array above. Returns NULL if the index/element is not present.
size_t                    md_vlx_xps_group_count(const struct md_vlx_t* vlx);
const md_vlx_xps_group_t* md_vlx_xps_group_by_index(const struct md_vlx_t* vlx, size_t idx);
const md_vlx_xps_group_t* md_vlx_xps_group_by_element(const struct md_vlx_t* vlx, md_element_t element);

// Atomic properties
size_t md_vlx_atomic_property_count(const struct md_vlx_t* vlx);
const md_vlx_atomic_property_t* md_vlx_atomic_property_by_index(const md_vlx_t* vlx, size_t idx);
const md_vlx_atomic_property_t* md_vlx_atomic_property_by_key(const md_vlx_t* vlx, uint64_t key);

// Density properties (Various propperties that are defined using AO densities)
size_t md_vlx_density_property_count(const struct md_vlx_t* vlx);
const md_vlx_density_property_t* md_vlx_density_property_by_index(const md_vlx_t* vlx, size_t idx);
const md_vlx_density_property_t* md_vlx_density_property_by_key(const md_vlx_t* vlx, uint64_t key);

// SYSTEM
bool md_vlx_system_init_from_data(struct md_system_t* sys, const md_vlx_t* vlx);
bool md_vlx_system_init_from_file(struct md_system_t* sys, str_t filename);
bool md_vlx_system_is_file_supplemental(const struct md_system_t* sys, str_t filename);

#ifdef __cplusplus
}
#endif
