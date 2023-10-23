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

static bool md_lammps_data_parse(md_lammps_data_t* data, md_buffered_reader_t* reader, struct md_allocator_i* alloc, data_format_t* data_format) {
	ASSERT(data);
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tokens[16];

	//lammps_atom_data_structure test_data_structure[] = { ATOM_IDX, MOL_IDX, ATOM_TYPE, PARTIAL_CHARGE, ATOM_COORD, ATOM_COORD, ATOM_COORD };

	//Read the title of the file
	if (!md_buffered_reader_extract_line(&line, reader)) {
		MD_LOG_ERROR("Failed to parse lammps title");
		return false;
	}
	str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

	//Skip empty line
	md_buffered_reader_skip_line(reader);

	//Read the number of atoms
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

	//Jump ahead to the Atoms definition
	int32_t atoms_line_counter = 0;
	do {
		md_buffered_reader_extract_line(&line, reader);
	} while (!str_equal_cstr_n(line, "Atoms", 5) && atoms_line_counter++ < 1000);

	//Skip empty line
	md_buffered_reader_skip_line(reader);

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
				atom->atom_idx = (int32_t)parse_int_wide(tokens[d].ptr, tokens[d].len);
				break;
			case MOL_IDX:
				atom->mol_idx = (int32_t)parse_int_wide(tokens[d].ptr, tokens[d].len);
				break;
			case ATOM_TYPE:
				atom->atom_type = (int32_t)parse_int_wide(tokens[d].ptr, tokens[d].len);
				break;
			case PARTIAL_CHARGE:
				atom->partial_charge = (float)parse_float_wide(tokens[d].ptr, tokens[d].len);
				break;
			case ATOM_X:
				atom->x = (float)parse_float_wide(tokens[d].ptr, tokens[d].len);
				break;
			case ATOM_Y:
				atom->y = (float)parse_float_wide(tokens[d].ptr, tokens[d].len);
				break;
			case ATOM_Z:
				atom->z = (float)parse_float_wide(tokens[d].ptr, tokens[d].len);
				break;
			case NX:
				atom->nx = (int32_t)parse_int_wide(tokens[d].ptr, tokens[d].len);
				break;
			case NY:
				atom->ny = (int32_t)parse_int_wide(tokens[d].ptr, tokens[d].len);
				break;
			case NZ:
				atom->nz = (int32_t)parse_int_wide(tokens[d].ptr, tokens[d].len);
				break;
			}
		}

		//str_copy_to_char_buf(atom->atom_type, sizeof(atom->atom_type), str_trim(str_substr(line, 5, 5)));
		//str_copy_to_char_buf(atom->atom_name, sizeof(atom->atom_name), str_trim(str_substr(line, 10, 5)));
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