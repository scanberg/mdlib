#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_molecule_t;
struct md_molecule_loader_i;

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
	NZ
} lammps_atom_data_format;

typedef struct data_format_t {
	const lammps_atom_data_format* ptr;
	int64_t len;
} data_format_t;

typedef enum atom_style {
	full
}atom_style;

data_format_t* get_data_format(atom_style style);

typedef struct md_lammps_atom_t {
	int32_t atom_idx;
	int32_t mol_idx;
	int32_t atom_type;
	float partial_charge;
	float x;
	float y;
	float z;
} md_lammps_atom_t;

typedef struct md_lammps_atom_mass_t {
	int32_t atom_idx;
	float mass;
} md_lammps_atom_mass_t;

typedef struct md_lammps_data_t {
	char title[256];
	float cell_ext[3];

	md_lammps_atom_mass_t* atom_type_mass;

	// Field data
	int64_t num_atoms;
	md_lammps_atom_t* atom_data;
} md_lammps_data_t;

// Parse a text-blob as LAMMPS
bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, struct md_allocator_i* alloc, data_format_t* data_format);
bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, struct md_allocator_i* alloc, data_format_t* data_format);
void md_lammps_data_free(md_lammps_data_t* data, struct md_allocator_i* alloc);

/*
// Molecule
bool md_lammps_molecule_init(struct md_molecule_t* mol, const md_lammps_data_t* lammps_data, struct md_allocator_i* alloc);

struct md_molecule_loader_i* md_lammps_molecule_api();
*/

#ifdef __cplusplus
}
#endif

