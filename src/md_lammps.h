#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_trajectory_i;
struct md_system_t;

//All lammps units should have a 1:1 mapping to the md_molecule according to https://docs.lammps.org/2001/units.html

// We follow the convention of MDAnalysis and encode the atom fields within a string
// https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/topology/LAMMPSParser.py
// The only fields extracted are "id, resid, type, q, x, y, z"
// Where "id, type, x, y, z" are required fields

typedef enum {
	MD_LAMMPS_ATOM_FORMAT_UNKNOWN,
	MD_LAMMPS_ATOM_FORMAT_ANGLE,
	MD_LAMMPS_ATOM_FORMAT_ATOMIC,
	MD_LAMMPS_ATOM_FORMAT_BODY,
	MD_LAMMPS_ATOM_FORMAT_BOND,
	MD_LAMMPS_ATOM_FORMAT_CHARGE,
	MD_LAMMPS_ATOM_FORMAT_DIPOLE,
	MD_LAMMPS_ATOM_FORMAT_DPD,
	MD_LAMMPS_ATOM_FORMAT_EDPD,
	MD_LAMMPS_ATOM_FORMAT_MDPD,
	MD_LAMMPS_ATOM_FORMAT_ELECTRON,
	MD_LAMMPS_ATOM_FORMAT_ELLIPSOID,
	MD_LAMMPS_ATOM_FORMAT_FULL,
	MD_LAMMPS_ATOM_FORMAT_LINE,
	MD_LAMMPS_ATOM_FORMAT_MESO,
	MD_LAMMPS_ATOM_FORMAT_MOLECULAR,
	MD_LAMMPS_ATOM_FORMAT_PERI,
	MD_LAMMPS_ATOM_FORMAT_SMD,
	MD_LAMMPS_ATOM_FORMAT_SPHERE,
	MD_LAMMPS_ATOM_FORMAT_TEMPLATE,
	MD_LAMMPS_ATOM_FORMAT_TRI,
	MD_LAMMPS_ATOM_FORMAT_WAVEPACKET,
	MD_LAMMPS_ATOM_FORMAT_COUNT
} md_lammps_atom_format_t;

// Representing atom type data
typedef struct md_lammps_atom_type_t {
	int32_t id;
	float mass;
	float radius; // Approximated from LJ or equivalent (if present)
} md_lammps_atom_type_t;

// Contains data about a single atom
typedef struct md_lammps_atom_t {
	int32_t id;
	int32_t resid;
	int32_t type;
	float charge;
	float x;
	float y;
	float z;
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

	size_t num_atoms;
	size_t num_bonds;
	size_t num_angles;
	size_t num_dihedrals;
	size_t num_impropers;

	size_t num_atom_types;
	size_t num_bond_types;
	size_t num_angle_types;
	size_t num_dihedral_types;
	size_t num_improper_types;

	md_lammps_cell_t cell;

	md_lammps_atom_type_t* atom_types;
	md_lammps_atom_t* atoms;
	md_lammps_bond_t* bonds;
	md_lammps_angle_t* angles;
	md_lammps_dihedral_t* dihedrals;
	md_lammps_dihedral_t* impropers;
} md_lammps_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// Get name or string representation of a atom format
const char** md_lammps_atom_format_names(void);
const char** md_lammps_atom_format_strings(void);

// Try to find a predefined atom format from input
md_lammps_atom_format_t md_lammps_atom_format_from_str(str_t str);
md_lammps_atom_format_t md_lammps_atom_format_from_file(str_t filename);

// validate a arbitrary atom format encoded in a string (atom_format)
// Writes a human readable error message to out_buf if the format is invalid.
// the format must contain the fields: 'id type x y z'
bool md_lammps_validate_atom_format(char* out_buf, size_t buf_cap, const char* atom_format);

// DATA PARSING
// atom_format is optional and will be used if supplied
bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, const char* atom_format, struct md_allocator_i* alloc);
bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, const char* atom_format, struct md_allocator_i* alloc);
void md_lammps_data_free(md_lammps_data_t* data, struct md_allocator_i* alloc);

// SYSTEM
bool md_lammps_system_init_from_data(struct md_system_t* sys, const md_lammps_data_t* lammps_data);
bool md_lammps_system_init_from_file(struct md_system_t* sys, str_t filename, const char* atom_format);
bool md_lammps_system_init_from_str (struct md_system_t* sys, str_t str,	  const char* atom_format);

// TRAJECTORY
bool md_lammps_trajectory_attach_from_file(struct md_system_t* sys, str_t filename, uint32_t flags);
struct md_trajectory_i* md_lammps_trajectory_create(str_t filename, struct md_allocator_i* alloc, uint32_t flags);

#ifdef __cplusplus
}
#endif

