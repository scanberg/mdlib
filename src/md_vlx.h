#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_types.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_molecule_t;
struct md_molecule_loader_i;
struct basis_set_t;
struct vec3_t;

typedef struct md_vlx_geom_t {
	size_t num_atoms;
	md_label_t* atom_symbol;
	md_element_t* atomic_number;
	double* coord_x;
	double* coord_y;
	double* coord_z;

	size_t num_alpha_electrons;
	size_t num_beta_electrons;
	int molecular_charge;
	int spin_multiplicity;
} md_vlx_geom_t;

typedef struct md_vlx_dipole_moment_t {
	str_t ident;
	double x;
	double y;
	double z;
} md_vlx_dipole_moment_t;

typedef struct md_vlx_basis_t {
	str_t ident;
	size_t num_contracted_basis_functions;
	size_t num_primitive_basis_functions;
	struct basis_set_t* basis_set;
} md_vlx_basis_t;

// Self Consistent Field
typedef struct md_vlx_scf_t {
	struct {
		size_t count;
		int* iteration;
		double* energy_total;
		double* energy_change;
		double* gradient_norm;
		double* max_gradient;
		double* density_change;
	} iter;

	// Converged parameters
	double total_energy;
	double electronic_energy;
	double nuclear_repulsion_energy;
	double gradient_norm;

	struct {
		struct {
			size_t count;
			double* data;
		} energies;

		struct {
			size_t count;
			double* data;
		} occupations;

		struct {
			size_t dim[2];
			double* data;
		} orbitals;
	} alpha;

	md_vlx_dipole_moment_t ground_state_dipole_moment;
} md_vlx_scf_t;

typedef struct md_vlx_rsp_t {
	size_t num_excited_states;
	str_t* ident;
	md_vlx_dipole_moment_t* electronic_transition_length;
	md_vlx_dipole_moment_t* electronic_transition_velocity;
	md_vlx_dipole_moment_t* magnetic_transition;
	double* absorption_ev;
	double* absorption_osc_str;
	double* electronic_circular_dichroism_cgs;
} md_vlx_rsp_t;

typedef struct md_vlx_data_t {
	md_vlx_geom_t  geom;
	md_vlx_basis_t basis;
	md_vlx_scf_t   scf;
	md_vlx_rsp_t   rsp;
	md_allocator_i* alloc;
} md_vlx_data_t;


typedef struct md_vlx_grid_t {
	float* data;
	int   dim[3];
	float origin[3];
	float stepsize[3];
} md_vlx_grid_t;

// RAW FUNCTIONS
bool md_vlx_data_parse_str(md_vlx_data_t* data, str_t string, struct md_allocator_i* alloc);
bool md_vlx_data_parse_file(md_vlx_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_vlx_data_free(md_vlx_data_t* data);

bool md_vlx_grid_evaluate(md_vlx_grid_t* grid, const md_vlx_geom_t* geom, const md_vlx_basis_t* basis, const double* mo_coeffs, size_t num_mo_coeffs);
bool md_vlx_grid_evaluate_sub(md_vlx_grid_t* grid, const int grid_idx_min[3], const int grid_idx_max[3], const md_vlx_geom_t* geom, const md_vlx_basis_t* basis, const double* mo_coeffs, size_t num_mo_coeffs);
bool md_vlx_get_mo(double* out_values, const vec3_t* coords, size_t num_coords, const md_vlx_geom_t* geom, const md_vlx_basis_t* basis, const double* mo_coeffs, size_t num_mo_coeffs);

// MOLECULE
bool md_vlx_molecule_init(struct md_molecule_t* mol, const md_vlx_data_t* data, struct md_allocator_i* alloc);

struct md_molecule_loader_i* md_vlx_molecule_api();

#ifdef __cplusplus
}
#endif