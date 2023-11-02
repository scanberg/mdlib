#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

//All lammps units should have a 1:1 mapping to the md_molecule according to https://docs.lammps.org/2001/units.html

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_molecule_t;
struct md_molecule_loader_i;

//Contains potential data format identifiers
typedef enum lammps_atom_data_format {
	ATOM_IDX,
	MOL_IDX,
	ATOM_TYPE,
	PARTIAL_CHARGE,
	ATOM_X,
	ATOM_Y,
	ATOM_Z,
	NX,
	NY,
	NZ,
	res_name,
	res_id
} lammps_atom_data_format;

static const lammps_atom_data_format DATA_FORMAT_FULL[] = {
	ATOM_IDX, MOL_IDX, ATOM_TYPE, PARTIAL_CHARGE, ATOM_X, ATOM_Y, ATOM_Z, NX, NY, NZ
};

//Contains ptr to array with a data_format definition
typedef struct data_format_t {
	const lammps_atom_data_format* ptr;
	int64_t len;
} data_format_t;

//Keyword for style that the data can have
typedef enum atom_style {
	full
}atom_style;

//Used to create the data_format with just a style keyword
bool get_data_format(data_format_t* format, atom_style style);

//

//Contains data about a single atom
typedef struct md_lammps_atom_t {
	int32_t atom_idx;
	int32_t mol_idx;
	int32_t atom_type;
	float partial_charge;
	float x;
	float y;
	float z;
	int32_t nx;
	int32_t ny;
	int32_t nz;
	//char res_name[8];
	//int32_t res_id; use mol_idx
} md_lammps_atom_t;

typedef struct md_lammps_atom_mass_t {
	int32_t atom_type;
	float mass;
} md_lammps_atom_mass_t;

float get_mass(md_lammps_atom_mass_t* masses, int32_t type, int32_t num_types);

typedef struct md_lammps_atom_bond_t {
	int64_t bond_idx;
	int32_t bond_type;
	int64_t first_atom_idx;
	int64_t second_atom_idx;
} md_lammps_atom_bond_t;

typedef struct md_lammps_cell {
	float x, y, z;
	float xy, xz, yz;
} md_lammps_cell;

typedef struct md_lammps_data_t {
	char title[256];

	md_lammps_cell cell_def;

	int64_t num_atoms;
	int32_t num_atom_types;
	int64_t num_bonds;
	int32_t num_bond_types;
	int64_t num_angles;
	int32_t num_angle_types;
	int64_t num_dihedrals;
	int32_t num_dihedral_types;

	md_lammps_atom_mass_t* atom_type_mass;
	md_lammps_atom_bond_t* bonds;

	md_lammps_atom_t* atom_data;
} md_lammps_data_t;

// Parse a text-blob as LAMMPS
bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, struct md_allocator_i* alloc, data_format_t* data_format);
bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, struct md_allocator_i* alloc, data_format_t* data_format);
void md_lammps_data_free(md_lammps_data_t* data, struct md_allocator_i* alloc);


// Molecule
bool md_lammps_molecule_init(struct md_molecule_t* mol, const md_lammps_data_t* lammps_data, struct md_allocator_i* alloc);

//struct md_molecule_loader_i* md_lammps_molecule_api();


#ifdef __cplusplus
}
#endif

