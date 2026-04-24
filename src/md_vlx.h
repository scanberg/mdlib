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

typedef enum {
	MD_VLX_SPIN_ALPHA = 0,
	MD_VLX_SPIN_BETA  = 1,
} md_vlx_spin_t;

typedef enum {
	MD_VLX_NTO_TYPE_PARTICLE = 0,
	MD_VLX_NTO_TYPE_HOLE = 1,
} md_vlx_nto_type_t;

typedef enum {
	MD_VLX_SCF_TYPE_UNKNOWN = 0,
	MD_VLX_SCF_TYPE_RESTRICTED,
	MD_VLX_SCF_TYPE_RESTRICTED_OPENSHELL,
	MD_VLX_SCF_TYPE_UNRESTRICTED,
} md_vlx_scf_type_t;

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
size_t md_vlx_number_of_alpha_electrons(const struct md_vlx_t* vlx);
size_t md_vlx_number_of_beta_electrons(const struct md_vlx_t* vlx);

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
size_t  md_vlx_scf_homo_idx(const struct md_vlx_t* vlx, md_vlx_spin_t type);
size_t  md_vlx_scf_lumo_idx(const struct md_vlx_t* vlx, md_vlx_spin_t type);

size_t  md_vlx_scf_number_of_atomic_orbitals   (const struct md_vlx_t* vlx);
size_t  md_vlx_scf_number_of_molecular_orbitals(const struct md_vlx_t* vlx);

const double* md_vlx_scf_mo_occupancy(const struct md_vlx_t* vlx, md_vlx_spin_t type);
const double* md_vlx_scf_mo_energy(const struct md_vlx_t* vlx, md_vlx_spin_t type);

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
const double* md_vlx_scf_mo_coefficients(const struct md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t type);

// Direct pointer to the AO coefficient vector for NTO lambda_idx / type (shell order).
// Returns NULL on failure.
const double* md_vlx_rsp_nto_coefficients(const struct md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type);

// The overlap matrix (S) is a square, symmetric matrix [N][N], this returns the length N
size_t  md_vlx_scf_overlap_matrix_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_overlap_matrix_data(const struct md_vlx_t* vlx);

// The density matrix is a square, symmetric matrix [N][N], this returns the length N
size_t md_vlx_scf_density_matrix_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_density_matrix_data(const struct md_vlx_t* vlx, md_vlx_spin_t type);

// ---------------------------------------------------------------------------
// Deprecated extraction API — prefer the direct pointer functions above
// ---------------------------------------------------------------------------

// Copy one MO coefficient vector into caller-supplied storage.
// out must hold at least md_vlx_scf_number_of_atomic_orbitals() doubles.
// Returns the number of AOs on success, 0 on failure.
size_t md_vlx_scf_mo_coefficients_extract(double* out, const struct md_vlx_t* vlx, size_t mo_idx, md_vlx_spin_t type);

// SCF History
size_t		  md_vlx_scf_history_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_energy(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_energy_diff(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_density_diff(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_gradient_norm(const struct md_vlx_t* vlx);
const double* md_vlx_scf_history_max_gradient(const struct md_vlx_t* vlx);

// SCF (Optional data)
// 
// If present, has the length of num of atoms
const double* md_vlx_scf_resp_charges(const struct md_vlx_t* vlx);

// RSP (Standard Linear Response, TD-DFT)
size_t md_vlx_rsp_number_of_excited_states(const struct md_vlx_t* vlx);

// RSP arrays with length of num excited states
const dvec3_t* md_vlx_rsp_electric_transition_dipole_moments(const struct md_vlx_t* vlx);
const dvec3_t* md_vlx_rsp_magnetic_transition_dipole_moments(const struct md_vlx_t* vlx);
const dvec3_t* md_vlx_rsp_velocity_transition_dipole_moments(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_rotatory_strengths(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_oscillator_strengths(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_absorption_ev(const struct md_vlx_t* vlx);

// RSP CPP (Complex Polarization Propagator)
size_t md_vlx_rsp_cpp_number_of_frequencies(const struct md_vlx_t* vlx);

// unit = eV
const double* md_vlx_rsp_cpp_frequencies(const struct md_vlx_t* vlx);

// Encodes unpolarized? response, i.e. the average of left and right polarized light
const double* md_vlx_rsp_cpp_sigma(const struct md_vlx_t* vlx);

// Encodes difference between left vs. right polarized light
const double* md_vlx_rsp_cpp_delta_epsilon(const struct md_vlx_t* vlx);

bool md_vlx_rsp_has_nto(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_nto_occupancy(const struct md_vlx_t* vlx, size_t nto_idx);
const double*  md_vlx_rsp_nto_lambdas(const md_vlx_t* vlx, size_t nto_idx);
const double*  md_vlx_rsp_nto_energy(const struct md_vlx_t* vlx, size_t nto_idx);

// VIB
// Many of the fields within the VIB portion has a length given by the degrees of freedom (D)
size_t md_vlx_vib_number_of_normal_modes(const struct md_vlx_t* vlx);

// Returns arrays of length D
const double* md_vlx_vib_ir_intensities(const struct md_vlx_t* vlx); 	// Unit: km/mol
const double* md_vlx_vib_frequencies(const struct md_vlx_t* vlx);		// Unit: cm^-1
const double* md_vlx_vib_reduced_masses(const struct md_vlx_t* vlx);	// Unit: amu
const double* md_vlx_vib_force_constants(const struct md_vlx_t* vlx);	// Unit: mdyne/Å (millidyne / Ångström)

// OPT
size_t md_vlx_opt_number_of_steps(const struct md_vlx_t* vlx);

// Returns atom coordinates for a given optimization step
const dvec3_t* md_vlx_opt_coordinates(const struct md_vlx_t* vlx, size_t opt_idx);
const double* md_vlx_opt_energies(const struct md_vlx_t* vlx);
const double* md_vlx_opt_nuclear_repulsion_energies(const struct md_vlx_t* vlx);

// The normal mode for each atom
// Returns array of length N (Number of atoms)
const dvec3_t* md_vlx_vib_normal_mode(const struct md_vlx_t* vlx, size_t idx);

size_t md_vlx_vib_num_raman_activity(const struct md_vlx_t* vlx);

// Returns array of length D
const double*  md_vlx_vib_raman_activity(const struct md_vlx_t* vlx, size_t idx);

// Atomic properties
size_t md_vlx_atomic_property_count(const struct md_vlx_t* vlx);
const md_vlx_atomic_property_t* md_vlx_atomic_property_by_index(const md_vlx_t* vlx, size_t idx);
const md_vlx_atomic_property_t* md_vlx_atomic_property_by_key(const md_vlx_t* vlx, uint64_t key);

// SYSTEM
bool md_vlx_system_init_from_data(struct md_system_t* sys, const md_vlx_t* vlx);
bool md_vlx_system_init_from_file(struct md_system_t* sys, str_t filename);
bool md_vlx_system_is_file_supplemental(const struct md_system_t* sys, str_t filename);

#ifdef __cplusplus
}
#endif
