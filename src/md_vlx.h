#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <md_types.h>
#include <md_gto.h>
#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;
struct md_system_loader_i;

typedef enum {
	MD_VLX_MO_TYPE_ALPHA = 0,
	MD_VLX_MO_TYPE_BETA = 1,
} md_vlx_mo_type_t;

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

typedef struct md_vlx_t md_vlx_t;

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

// SCF
md_vlx_scf_type_t md_vlx_scf_type(const struct md_vlx_t* vlx);
dvec3_t md_vlx_scf_ground_state_dipole_moment(const struct md_vlx_t* vlx);
size_t  md_vlx_scf_homo_idx(const struct md_vlx_t* vlx, md_vlx_mo_type_t type);
size_t  md_vlx_scf_lumo_idx(const struct md_vlx_t* vlx, md_vlx_mo_type_t type);

size_t  md_vlx_scf_number_of_atomic_orbitals   (const struct md_vlx_t* vlx);
size_t  md_vlx_scf_number_of_molecular_orbitals(const struct md_vlx_t* vlx);

const double* md_vlx_scf_mo_occupancy(const struct md_vlx_t* vlx, md_vlx_mo_type_t type);
const double* md_vlx_scf_mo_energy(const struct md_vlx_t* vlx, md_vlx_mo_type_t type);

// Atomic orbital GTO data
bool md_vlx_scf_extract_gto_data(md_gto_data_t* out_gto_data, const struct md_vlx_t* vlx, double cutoff_value, struct md_allocator_i* alloc);

// The overlap matrix (S) is a square, symmetric matrix [N][N], this returns the length N
size_t  md_vlx_scf_overlap_matrix_size(const struct md_vlx_t* vlx);
const double* md_vlx_scf_overlap_matrix_data(const struct md_vlx_t* vlx);

// Get the required element size of a compact upper triangular matrix representation
size_t md_vlx_scf_upper_triangular_density_matrix_size(const struct md_vlx_t* vlx);

// Extracts the elements of the matrix into a compact upper triangular matrix
bool md_vlx_scf_extract_upper_triangular_density_matrix_data(float* out_values, const struct md_vlx_t* vlx, md_vlx_mo_type_t type);

// Get the regular density matrix size N in (N x N)
size_t md_vlx_scf_density_matrix_size(const struct md_vlx_t* vlx);

// Extracts the full density matrix into a square matrix representation
bool md_vlx_scf_extract_density_matrix_data(float* out_values, const struct md_vlx_t* vlx, md_vlx_mo_type_t type);

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

// RSP
size_t md_vlx_rsp_number_of_excited_states(const struct md_vlx_t* vlx);

// RSP arrays with length of num excited states
const dvec3_t* md_vlx_rsp_electric_transition_dipole_moments(const struct md_vlx_t* vlx);
const dvec3_t* md_vlx_rsp_magnetic_transition_dipole_moments(const struct md_vlx_t* vlx);
const dvec3_t* md_vlx_rsp_velocity_transition_dipole_moments(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_rotatory_strengths(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_oscillator_strengths(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_absorption_ev(const struct md_vlx_t* vlx);

bool md_vlx_rsp_has_nto(const struct md_vlx_t* vlx);
const double*  md_vlx_rsp_nto_occupancy(const struct md_vlx_t* vlx, size_t nto_idx);
const double*  md_vlx_rsp_nto_lambdas(const md_vlx_t* vlx, size_t nto_idx);
const double*  md_vlx_rsp_nto_energy(const struct md_vlx_t* vlx, size_t nto_idx);

// Extracts a submatrix of the density matrix for a given NTO and atom
// out_matrix: Preallocated matrix of size (matrix_size x matrix_size)
// matrix_size: The size of the square matrix to extract
// vlx: A valid VeloxChem object
// nto_idx: The index of the NTO to extract from
// atom_idx: The atom index to extract the submatrix for
// type: The NTO type (Particle / Hole)
void md_vlx_rsp_nto_extract_density_submatrix_for_atom(float* out_matrix, size_t matrix_size, const struct md_vlx_t* vlx, size_t nto_idx, size_t atom_idx, md_vlx_nto_type_t type);

// returns the AO coefficients for the given NTO, lambda idx and type
size_t md_vlx_rsp_nto_extract_coefficients(double* out_ao, const struct md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type);

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

// Extract Natural Transition Orbitals GTOs
// nto_idx: The index of the excited state (0-based indexing)
// lambda_idx: The lambda component to extract (0-based indexing), 0 corresponds to the most significant index
size_t md_vlx_nto_gto_count(const md_vlx_t* vlx);
bool   md_vlx_nto_gto_extract(md_gto_t* gtos, const md_vlx_t* vlx, size_t nto_idx, size_t lambda_idx, md_vlx_nto_type_t type);

// Extract Molecular Orbital (MO) PGTOs

// This provides an upper limit to the number of gtos supplied
// Use this to reserve the data for the gtos
size_t md_vlx_mo_gto_count(const md_vlx_t* vlx);

// This extract the gtos into the supplied array
// Returns the number of gtos written
// gtos: an array to hold the extracted gtos
// vlx: a pointer to a valid VeloxChem object
// mo_idx: Molecular Orbital Index
// type: Molecular orbital type: Alpha / Beta
// value_cutoff: A cutoff value which will be used to calculate an effective radius of influence for the gtos (0 == no cutoff (Infinite radius of influence))
size_t md_vlx_mo_gto_extract(md_gto_t gtos[], const md_vlx_t* vlx, size_t mo_idx, md_vlx_mo_type_t type, double value_cutoff);

// MOLECULE
bool md_vlx_system_init(struct md_system_t* sys, const md_vlx_t* vlx, struct md_allocator_i* alloc);

struct md_system_loader_i* md_vlx_system_loader(void);


#ifdef __cplusplus
}
#endif
