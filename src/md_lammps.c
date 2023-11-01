#include <md_lammps.h>

#include <md_util.h>

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_array.h>
#include <core/md_parse.h>

#include <string.h>

bool get_data_format(data_format_t* format, atom_style style) {
	switch (style)
	{
	default:
		MD_LOG_ERROR("Invalid atom_style");
		return false;
	case full:
		format->ptr = DATA_FORMAT_FULL;
		format->len = ARRAY_SIZE(DATA_FORMAT_FULL);
		return true;
	}
}


float get_mass(md_lammps_atom_mass_t* m, int32_t type, int32_t num_types)
{
	for (int32_t i = 0; i < num_types; i++) {
		if (m[i].atom_type == type) {
			return m[i].mass;
		}
	}
	MD_LOG_ERROR("atom type could not be found");
	return -1.0f;
}



static bool md_lammps_data_parse(md_lammps_data_t* data, md_buffered_reader_t* reader, struct md_allocator_i* alloc, data_format_t* data_format) {
	ASSERT(data);
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tokens[16];

	//lammps_atom_data_structure test_data_structure[] = { ATOM_IDX, MOL_IDX, ATOM_TYPE, PARTIAL_CHARGE, ATOM_COORD, ATOM_COORD, ATOM_COORD };

	//Read the title of the file
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse lammps title");
			return false;
		}
		str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

		//Skip empty line
		md_buffered_reader_skip_line(reader);
	}

	//Read the number of atoms
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num atoms");
			return false;
		}
		data->num_atoms = parse_int(str_trim(line));
		if (!data->num_atoms) {
			MD_LOG_ERROR("Failed to parse LAMMPS num atoms");
			return false;
		}

		md_array_resize(data->atom_data, data->num_atoms, alloc);
	}

	//Read num atom types
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num atom types");
			return false;
		}
		data->num_atom_types = (int32_t)parse_int(str_trim(line));
		if (!data->num_atom_types) {
			MD_LOG_ERROR("Failed to parse LAMMPS num atom types");
			return false;
		}
		md_array_resize(data->atom_type_mass, data->num_atom_types, alloc);
	}

	//Read num bonds
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num bonds");
			return false;
		}
		data->num_bonds = parse_int(str_trim(line));
		if (!data->num_bonds) {
			MD_LOG_ERROR("Failed to parse LAMMPS num bonds");
			return false;
		}

		md_array_resize(data->bonds, data->num_bonds, alloc);
	}

	//Read num bond types
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num bond types");
			return false;
		}
		data->num_bond_types = (int32_t)parse_int(str_trim(line));
		if (!data->num_bond_types) {
			MD_LOG_ERROR("Failed to parse LAMMPS num bond types");
			return false;
		}
	}

	//Read num angles
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num angles");
			return false;
		}
		data->num_angles = parse_int(str_trim(line));
		if (!data->num_angles) {
			MD_LOG_ERROR("Failed to parse LAMMPS num angles");
			return false;
		}
	}

	//Read num angle types
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num angle types");
			return false;
		}
		data->num_angle_types = (int32_t)parse_int(str_trim(line));
		if (!data->num_angle_types) {
			MD_LOG_ERROR("Failed to parse LAMMPS num angle types");
			return false;
		}
	}

	//Read num dihedrals
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num dihedrals");
			return false;
		}
		data->num_dihedrals = parse_int(str_trim(line));
		if (!data->num_dihedrals) {
			MD_LOG_ERROR("Failed to parse LAMMPS num dihedrals");
			return false;
		}
	}

	//Read num dihedral types
	{
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to parse LAMMPS num dihedral types");
			return false;
		}
		data->num_dihedral_types = (int32_t)parse_int(str_trim(line));
		if (!data->num_dihedral_types) {
			MD_LOG_ERROR("Failed to parse LAMMPS num dihedral types");
			return false;
		}
	}

	//Skip empty line
	md_buffered_reader_skip_line(reader);

	//Parse cell extent
	{
		//Cubic part
		for (int32_t i = 0; i < 3; i++) {
			if (!md_buffered_reader_extract_line(&line, reader)) {
				MD_LOG_ERROR("Could not read cell extent line");
				return false;
			}
			if (extract_tokens(tokens, 2, &line) != 2) {
				MD_LOG_ERROR("Wrong amount of extent tokens tokens");
				return false;
			}
			data->cell_ext[i] = (float)parse_float(tokens[1]) - (float)parse_float(tokens[0]);
		}
		
		//Triclinic part
		md_buffered_reader_extract_line(&line, reader);
		if (extract_tokens(tokens, 3, &line) == 3){
			data->xy = (float)parse_float(tokens[0]);
			data->xz = (float)parse_float(tokens[1]);
			data->yz = (float)parse_float(tokens[2]);
		}
		else {
			//It was cubic
			data->xy = 0;
			data->xz = 0;
			data->yz = 0;
		}
		
	}

	//Mass parsing
	{
		//Jump ahead to the Masses definition
		{
			do {

				if (!md_buffered_reader_extract_line(&line, reader)) {
					MD_LOG_ERROR("Could not find Masses line");
					return false;
				}
			} while (!str_equal_cstr_n(line, "Masses", 6));

			//Skip empty line
			md_buffered_reader_skip_line(reader);
		}

		//Start reading the Masses data
		for (int32_t i = 0; i < data->num_atom_types; i++) {
			if (!md_buffered_reader_extract_line(&line, reader)) {
				MD_LOG_ERROR("Failed to extract mass line");
				return false;
			}
			if (extract_tokens(tokens, 2, &line) != 2) {
				MD_LOG_ERROR("Wrong amount of mass tokens");
				return false;
			}

			md_lammps_atom_mass_t* mass = &data->atom_type_mass[i];

			mass->atom_type = (int32_t)parse_int(tokens[0]);
			mass->mass = (float)parse_float(tokens[1]);
		}
	}

	//Atom parsing
	{
		//Jump ahead to the Atoms definition
		{
			do {

				if (!md_buffered_reader_extract_line(&line, reader)) {
					MD_LOG_ERROR("Could not find Atoms line");
					return false;
				}
			} while (!str_equal_cstr_n(line, "Atoms", 5));

			//Skip empty line
			md_buffered_reader_skip_line(reader);
		}

		//Start reading the atom data
		for (int64_t i = 0; i < data->num_atoms; ++i) {
			if (!md_buffered_reader_extract_line(&line, reader)) {
				MD_LOG_ERROR("Failed to extract atom line");
				return false;
			}
			const int64_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
			if (num_tokens < data_format->len) {
				MD_LOG_ERROR("Failed to parse atom coordinates, expected at least %i tokens, got %i", (int)data_format->len, (int)num_tokens);
				return false;
			}

			md_lammps_atom_t* atom = &data->atom_data[i];

			/*
			ATOM_IDX,
			MOL_IDX,
			ATOM_TYPE,
			PARTIAL_CHARGE,
			ATOM_COORD
			*/

			for (int64_t d = 0; d < data_format->len; d++) {
				switch (data_format->ptr[d])
				{
				default:
					MD_LOG_ERROR("Failed to read atom data format");
					return false;
				case ATOM_IDX:
					atom->atom_idx = (int32_t)parse_int(tokens[d]);
					break;
				case MOL_IDX:
					atom->mol_idx = (int32_t)parse_int(tokens[d]);
					break;
				case ATOM_TYPE:
					atom->atom_type = (int32_t)parse_int(tokens[d]);
					break;
				case PARTIAL_CHARGE:
					atom->partial_charge = (float)parse_float(tokens[d]);
					break;
				case ATOM_X:
					atom->x = (float)parse_float(tokens[d]);
					break;
				case ATOM_Y:
					atom->y = (float)parse_float(tokens[d]);
					break;
				case ATOM_Z:
					atom->z = (float)parse_float(tokens[d]);
					break;
				case NX:
					atom->nx = (int32_t)parse_int(tokens[d]);
					break;
				case NY:
					atom->ny = (int32_t)parse_int(tokens[d]);
					break;
				case NZ:
					atom->nz = (int32_t)parse_int(tokens[d]);
					break;
				}
			}

			//str_copy_to_char_buf(atom->atom_type, sizeof(atom->atom_type), str_trim(str_substr(line, 5, 5)));
			//str_copy_to_char_buf(atom->atom_name, sizeof(atom->atom_name), str_trim(str_substr(line, 10, 5)));
		}
	}

	//Bonds parsing
	{
		//Jump ahead to the Bonds definition
		{
			do {

				if (!md_buffered_reader_extract_line(&line, reader)) {
					MD_LOG_ERROR("Could not find Bonds line");
					return false;
				}
			} while (!str_equal_cstr_n(line, "Bonds", 5));
			//Skip empty line
			md_buffered_reader_skip_line(reader);
		}

		//Start reading the bonds
		for (int32_t i = 0; i < data->num_bonds; i++) {
			if (!md_buffered_reader_extract_line(&line, reader)) {
				MD_LOG_ERROR("Failed to extract bond line");
				return false;
			}
			if (extract_tokens(tokens, 4, &line) != 4) {
				MD_LOG_ERROR("Wrong amount of bond tokens");
				return false;
			}

			md_lammps_atom_bond_t* bond = &data->bonds[i];

			bond->bond_idx = (int64_t)parse_int(tokens[0]);
			bond->bond_type = (int32_t)parse_int(tokens[1]);
			bond->first_atom_idx = (int64_t)parse_int(tokens[2]);
			bond->second_atom_idx = (int64_t)parse_int(tokens[3]);
		}
	}

	/*

	if (!md_buffered_reader_extract_line(&line, reader)) {
		MD_LOG_ERROR("Failed to extract unitcell line");
		return false;
	}

	const int64_t num_tokens = extract_float_tokens(tokens, ARRAY_SIZE(tokens), line);

	if (num_tokens != 3) {
		MD_LOG_ERROR("Failed to parse cell extent, expected 3 tokens, got %i", (int)num_tokens);
		return false;
	}

	data->cell_ext[0] = (float)parse_float(tokens[0]);
	data->cell_ext[1] = (float)parse_float(tokens[1]);
	data->cell_ext[2] = (float)parse_float(tokens[2]);

	*/

	return true;
}

bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, struct md_allocator_i* alloc, data_format_t* data_format) {
	ASSERT(data);
	ASSERT(alloc);

	md_buffered_reader_t line_reader = md_buffered_reader_from_str(str);
	return md_lammps_data_parse(data, &line_reader, alloc, data_format);
}

bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, struct md_allocator_i* alloc, data_format_t* data_format) {
	bool result = false;
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (file) {
		const int64_t cap = MEGABYTES(1);
		char* buf = md_alloc(md_heap_allocator, cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = md_lammps_data_parse(data, &line_reader, alloc, data_format);

		md_free(md_heap_allocator, buf, cap);
		md_file_close(file);
	}
	else {
		MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
	}
	return result;
}

void md_lammps_data_free(md_lammps_data_t* data, struct md_allocator_i* alloc) {
	ASSERT(data);
	if (data->atom_data) md_array_free(data->atom_data, alloc);
	MEMSET(data, 0, sizeof(md_lammps_data_t));
}

bool md_lammps_molecule_init(md_molecule_t* mol, const md_lammps_data_t* data, md_allocator_i* alloc)
{
	ASSERT(mol);
	ASSERT(data);
	ASSERT(alloc);

	MEMSET(mol, 0, sizeof(md_molecule_t));

	const int64_t num_atoms = data->num_atoms;

	md_array_ensure(mol->atom.x, num_atoms, alloc);
	md_array_ensure(mol->atom.y, num_atoms, alloc);
	md_array_ensure(mol->atom.z, num_atoms, alloc);
	md_array_ensure(mol->atom.residue_idx, num_atoms, alloc);

	int32_t cur_res_id = -1;
	for (int64_t i = 0; i < num_atoms; ++i) {
		const float x = data->atom_data[i].x * 10.0f; // convert from nm to Ångström
		const float y = data->atom_data[i].y * 10.0f; // convert from nm to Ångström
		const float z = data->atom_data[i].z * 10.0f; // convert from nm to Ångström

		int32_t res_id = data->atom_data[i].mol_idx;
		if (res_id != cur_res_id) {
			cur_res_id = res_id;
			md_residue_id_t id = res_id;
			md_range_t atom_range = { (uint32_t)mol->atom.count, (uint32_t)mol->atom.count };

			mol->residue.count += 1;
			md_array_push(mol->residue.id, id, alloc);
			md_array_push(mol->residue.atom_range, atom_range, alloc);
		}

		if (mol->residue.atom_range) md_array_last(mol->residue.atom_range)->end += 1;

		mol->atom.count += 1;

		//Set coordinates
		md_array_push(mol->atom.x, x, alloc);
		md_array_push(mol->atom.y, y, alloc);
		md_array_push(mol->atom.z, z, alloc);

		if (mol->residue.count) md_array_push(mol->atom.residue_idx, (md_residue_idx_t)(mol->residue.count - 1), alloc);

		//Set mass
		const float mass = get_mass(data->atom_type_mass, data->atom_data[i].atom_type, data->num_atom_types);
		md_array_push(mol->atom.mass, mass, alloc);
	}

	//Set elements
	md_array_resize(mol->atom.element, num_atoms, alloc);
	if (!md_util_element_from_mass(mol->atom.element, mol->atom.mass, num_atoms)) MD_LOG_ERROR("One or more masses are missing matching element");

	//Create unit cell
	mol->unit_cell = md_util_unit_cell_from_triclinic(data->cell_ext[0] * 10.0, data->cell_ext[1] * 10.0, data->cell_ext[2] * 10.0, data->xy * 10.0, data->xz * 10.0, data->yz * 10.0);
	return true;
}

bool md_lammps_init_from_str(md_molecule_t* mol, str_t str, md_allocator_i* alloc, data_format_t* format) {
	md_lammps_data_t data = { 0 };
	bool success = false;
	if (md_lammps_data_parse_str(&data, str, md_heap_allocator, format)) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_heap_allocator);

	return success;
}

bool md_lammps_init_from_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc, data_format_t* format) {
	md_lammps_data_t data = { 0 };
	bool success = false;
	if (md_lammps_data_parse_file(&data, filename, md_heap_allocator, format)) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_heap_allocator);

	return success;
}

/*
static md_molecule_loader_i lammps_api = {
	lammps_init_from_str,
	lammps_init_from_file,
};


md_molecule_loader_i* md_lammps_molecule_api() {
	return &lammps_api;
}
*/