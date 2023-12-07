#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <md_util.h>
#include <md_trajectory.h>

struct md_molecule_t;

//All lammps units should have a 1:1 mapping to the md_molecule according to https://docs.lammps.org/2001/units.html

// We follow the convention of MDAnalysis and encode the atom fields within a string
// https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/topology/LAMMPSParser.py
// The only fields extracted are "id, resid, type, charge, x, y, z"
// Where "id, type, x, y, z" are required fields

#define MD_LAMMPS_ATOM_FORMAT_FULL   "id resid type charge x y z"
#define MD_LAMMPS_ATOM_FORMAT_ATOMIC "id type x y z"

#ifdef __cplusplus
extern "C" {
#endif

//Contains data about a single atom
typedef struct md_lammps_atom_t {
	int32_t id;
	int32_t resid;
	int32_t type;
	float charge;
	float x;
	float y;
	float z;
	float mass;
} md_lammps_atom_t;

typedef struct md_lammps_bond_t {
	int32_t id;
	int32_t type;
	int32_t atom_id[2];
} md_lammps_bond_t;

typedef struct md_lammps_angle_t {
	int32_t id;
	int32_t type;
	int32_t atom_id[3];
} md_lammps_angle_t;

typedef struct md_lammps_dihedral_t {
	int32_t id;
	int32_t type;
	int32_t atom_id[4];
} md_lammps_dihedral_t;

typedef struct md_lammps_cell_t {
	float xlo, xhi;
	float ylo, yhi;
	float zlo, zhi;
	float xy, xz, yz;
} md_lammps_cell_t;

typedef struct md_lammps_data_t {
	char title[256];

	int32_t num_atoms;
	int32_t num_bonds;
	int32_t num_angles;
	int32_t num_dihedrals;
	int32_t num_impropers;

	int32_t num_atom_types;
	int32_t num_bond_types;
	int32_t num_angle_types;
	int32_t num_dihedral_types;
	int32_t num_improper_types;

	md_lammps_cell_t cell;

	md_lammps_atom_t* atoms;
	md_lammps_bond_t* bonds;
	md_lammps_angle_t* angles;
	md_lammps_dihedral_t* dihedrals;
	md_lammps_dihedral_t* impropers;
} md_lammps_data_t;

bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, const char* atom_format, struct md_allocator_i* alloc);
bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, const char* atom_format, struct md_allocator_i* alloc);
void md_lammps_data_free(md_lammps_data_t* data, struct md_allocator_i* alloc);

// Molecule
bool md_lammps_molecule_init(struct md_molecule_t* mol, const md_lammps_data_t* lammps_data, struct md_allocator_i* alloc);

//struct md_molecule_loader_i* md_lammps_molecule_api();

//Trajectory
struct md_trajectory_i* md_lammps_trajectory_create(str_t filename, struct md_allocator_i* alloc);
void md_lammps_trajectory_free(struct md_trajectory_i* traj);

struct md_trajectory_loader_i* md_lammps_trajectory_loader();

#ifdef __cplusplus
}
#endif

