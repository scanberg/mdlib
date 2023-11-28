#include <md_lammps.h>

//#include <md_util.h>

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
//#include <core/md_os.h>
#include <core/md_array.h>
#include <core/md_parse.h>
//#include <md_trajectory.h>

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


static float get_mass(md_lammps_atom_mass_t* m, int32_t type, int32_t num_types)
{
	for (int32_t i = 0; i < num_types; i++) {
		if (m[i].atom_type == type) {
			return m[i].mass;
		}
	}
	MD_LOG_ERROR("atom type could not be found");
	return 0.0f;
}



static bool md_lammps_data_parse(md_lammps_data_t* data, md_buffered_reader_t* reader, struct md_allocator_i* alloc, data_format_t* data_format) {
	ASSERT(data);
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tokens[16];


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

	//Parse cell definition
	{
		//Cubic part
		float cell_xyz[3] = { 0 };
		for (int32_t i = 0; i < 3; i++) {
			if (!md_buffered_reader_extract_line(&line, reader)) {
				MD_LOG_ERROR("Could not read cell extent line");
				return false;
			}
			if (extract_tokens(tokens, 2, &line) != 2) {
				MD_LOG_ERROR("Wrong amount of extent tokens tokens");
				return false;
			}
			cell_xyz[i] = (float)parse_float(tokens[1]) - (float)parse_float(tokens[0]);
		}
		data->cell_def.x = cell_xyz[0];
		data->cell_def.y = cell_xyz[1];
		data->cell_def.z = cell_xyz[2];

		//Triclinic part
		md_buffered_reader_extract_line(&line, reader);
		if (extract_tokens(tokens, 3, &line) == 3) {
			data->cell_def.xy = (float)parse_float(tokens[0]);
			data->cell_def.xz = (float)parse_float(tokens[1]);
			data->cell_def.yz = (float)parse_float(tokens[2]);
		}
		else {
			//It was cubic
			data->cell_def.xy = 0;
			data->cell_def.xz = 0;
			data->cell_def.yz = 0;
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
	md_array_ensure(mol->atom.flags, num_atoms, alloc);
	md_array_ensure(mol->atom.resid, num_atoms, alloc);
	md_array_ensure(mol->atom.mass, num_atoms, alloc);

	int32_t prev_res_id = -1;
	for (int64_t i = 0; i < num_atoms; ++i) {
		const float x = data->atom_data[i].x;
		const float y = data->atom_data[i].y;
		const float z = data->atom_data[i].z;
		const md_residue_id_t res_id = data->atom_data[i].mol_idx;
		const float mass = get_mass(data->atom_type_mass, data->atom_data[i].atom_type, data->num_atom_types);
		int flags = 0;

		if (prev_res_id != res_id) {
			flags |= MD_FLAG_RES_BEG;

			if (prev_res_id != -1) {
				*md_array_last(mol->atom.flags) |= MD_FLAG_RES_END;
			}
			prev_res_id = res_id;
		}

		md_array_push(mol->atom.x, x, alloc);
		md_array_push(mol->atom.y, y, alloc);
		md_array_push(mol->atom.z, z, alloc);
		md_array_push(mol->atom.flags, flags, alloc);
		md_array_push(mol->atom.resid, res_id, alloc);
		md_array_push(mol->atom.mass, mass, alloc);
	}

	mol->atom.count = num_atoms;

	//Set elements
	md_array_resize(mol->atom.element, num_atoms, alloc);
	if (!md_util_element_from_mass(mol->atom.element, mol->atom.mass, num_atoms)) {
		MD_LOG_ERROR("One or more masses are missing matching element");
	}

	//Create unit cell
	mol->unit_cell = md_util_unit_cell_from_triclinic(data->cell_def.x, data->cell_def.y, data->cell_def.z, data->cell_def.xy, data->cell_def.xz, data->cell_def.yz);

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


//Cant use interface as the format parameter is needed as well
/*
static md_molecule_loader_i lammps_api = {
	lammps_init_from_str,
	lammps_init_from_file,
};


md_molecule_loader_i* md_lammps_molecule_api() {
	return &lammps_api;
}
*/

//Reads data that is useful later when we want to parse a frame from the trajectory
static bool md_lammps_trajectory_parse(md_lammps_trajectory_t* traj, md_buffered_reader_t* reader, struct md_allocator_i* alloc) {
	ASSERT(traj);
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tokens[6];
	int64_t row = 0; //We start on row 0. When you read + handle a row, or skip a row, row should be increased.
	int64_t num_frames = 0;
	int64_t num_atoms = 0;
	md_unit_cell_t unit_cell = { 0 };

	int64_t array_max_size = 10000;
	//Create an array with the frame times. We make for 10 000 frames as we dont know the max number before parsing
	md_array(double) frame_times = md_array_create(double, array_max_size, alloc);

	int64_t* frame_offsets = 0;

	

	//Setup frames and first frame

	md_buffered_reader_skip_line(reader);
	row++;
	md_buffered_reader_extract_line(&line, reader); //Read from row 1
	frame_times[0] = (double)parse_float(line);
	num_frames++;
	row++;


	//Parse num of atoms
	{
		md_buffered_reader_skip_line(reader);
		row++;
		md_buffered_reader_extract_line(&line, reader); //Read num of atoms
		num_atoms = parse_int(line);
		
		row++;
	}

	//Parse unit_cell definition
	{
		float cell_extent[3] = 0;
		float cell_tri[3] = 0; //Triclinic

		md_buffered_reader_extract_line(&line, reader); //Read BOX BOUNDS
		if (str_equal_cstr(line, "ITEM: BOX BOUNDS pp pp pp")) {
			//Cubic
			row++;
			for (int8_t i = 0; i < 3; i++) {
				md_buffered_reader_extract_line(&line, reader);
				int64_t num_tokens = extract_tokens(tokens, 4, &line);
				cell_extent[i] = parse_float(tokens[1]) - parse_float(tokens[0]);
				row++;
			}
			unit_cell = md_util_unit_cell_from_extent(cell_extent[0], cell_extent[1], cell_extent[2]);
		}
		else if (str_equal_cstr(line, "ITEM: BOX BOUNDS xy xz yz pp pp pp")) {
			//Triclinic
			row++;
			for (int8_t i = 0; i < 3; i++) {
				md_buffered_reader_extract_line(&line, reader);
				int64_t num_tokens = extract_tokens(tokens, 4, &line);
				cell_extent[i] = parse_float(tokens[1]) - parse_float(tokens[0]);
				cell_tri[i] = parse_float(tokens[2]);
				row++;
			}
			unit_cell = md_util_unit_cell_from_triclinic(cell_extent[0], cell_extent[1], cell_extent[2], cell_tri[0], cell_tri[1], cell_tri[2]);
		}
		else {
			MD_LOG_ERROR("Could not correctly parse BOX BOUND");
			return false;
		}
	}


	while (md_buffered_reader_extract_line(&line, reader)) {

		if (str_equal_cstr_n(line, "ITEM: TIMESTEP", 14)) {
			row++;
			md_buffered_reader_extract_line(&line, reader);
			num_frames++;
			if (num_frames > array_max_size) {
				MD_LOG_ERROR("num of frames > %i, max array size to low", array_max_size);
				return false;
			}
			frame_times[num_frames - 1] = (double)parse_float(line);
			row++;
		} 
		else if (str_equal_cstr_n(line, "ITEM: ATOMS", 11)) {
			row++; //We increment before, as the offset should indicate the first row with atom information
			md_array_push(frame_offsets, row, alloc);

			//TODO: Implement this instead
			// This is a bit nasty, we want to get the correct offset to the beginning of the current line.
			// Therefore we need to do some pointer arithmetic because just using the length of the line may not get us
			// all the way back in case there were skipped \r characters.
			const int64_t offset = md_buffered_reader_tellg(reader) - (reader->str.ptr - line.ptr);
		}
		else {
			//Nothing of interest, just increment the row
			row++;
		}
	}

	md_array_resize(frame_times, num_frames, alloc);
	//Lammps trajectory uses femtoseconds
	//traj->header.time_unit.mult = 1.0e-15;
	//traj->header.num_frames = num_frames;
	//traj->header.frame_times = frame_times;
	//traj->header.num_atoms = num_atoms;

	traj->unit_cell = unit_cell;
	traj->frame_offsets = frame_offsets;

	traj->magic = MD_LAMMPS_TRAJ_MAGIC;
	//traj->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	//pdb->filesize = filesize;
	//pdb->frame_offsets = offsets;
	traj->allocator = alloc;
	traj->mutex = md_mutex_create();
	traj->header = (md_trajectory_header_t){
		.num_frames = num_frames,
		.num_atoms = num_atoms,
		.max_frame_data_size = max_frame_size,
		.time_unit.mult = 1.0e-15,
		.frame_times = frame_times,
	};
	//pdb->unit_cell = cell;

	return true;
}

bool md_lammps_trajectory_parse_file(md_lammps_trajectory_t* traj, str_t filename, struct md_allocator_i* alloc) {
	bool result = false;
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (file) {
		const int64_t cap = MEGABYTES(1);
		char* buf = md_alloc(md_heap_allocator, cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = md_lammps_trajectory_parse(traj, &line_reader, alloc);

		md_free(md_heap_allocator, buf, cap);
		md_file_close(file);
	}
	else {
		MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
	}
	return result;
}

//Loads the frame data from a given frame index, reading only the relevant part file for that frame.
bool lammps_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
	ASSERT(inst);

	md_lammps_trajectory_t* lammps_traj = (md_lammps_trajectory_t*)inst;
	if (lammps_traj->magic != MD_LAMMPS_TRAJ_MAGIC) {
		MD_LOG_ERROR("Error when decoding frame coord, xtc magic did not match");
		return false;
	}

	// Should this be exposed?
	md_allocator_i* alloc = md_temp_allocator;

	bool result = true;
	const int64_t frame_size = lammps_fetch_frame_data(inst, frame_idx, NULL);
	if (frame_size > 0) { //This seems to check first that there actually is data to be read, before writing to frame_data
		// This is a borderline case if one should use the md_temp_allocator as the raw frame size could potentially be several megabytes...
		void* frame_data = md_alloc(alloc, frame_size);
		const int64_t read_size = pdb_fetch_frame_data(inst, frame_idx, frame_data);
		if (read_size != frame_size) {
			MD_LOG_ERROR("Failed to read the expected size");
			md_free(alloc, frame_data, frame_size);
			return false;
		}

		//Decode is what actually writes the data to x, y and z which are float pointers
		result = lammps_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
		md_free(alloc, frame_data, frame_size);
	}

	return result;
}

static bool try_read_cache(str_t cache_file, md_array(int64_t)* offsets, int64_t* num_atoms, md_unit_cell_t* cell, md_allocator_i* alloc) {
	md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
	bool result = false;
	if (file) {
		int64_t num_offsets = 0;
		uint32_t magic = 0;
		uint32_t version = 0;

		if (md_file_read(file, &magic, sizeof(magic)) != sizeof(magic) || magic != MD_LAMMPS_CACHE_MAGIC) {
			MD_LOG_ERROR("Failed to read offset cache, magic was incorrect or corrupt");
			goto done;
		}

		if (md_file_read(file, &version, sizeof(version)) != sizeof(version) || version != MD_LAMMPS_CACHE_VERSION) {
			MD_LOG_ERROR("Failed to read offset cache, version was incorrect");
			goto done;
		}

		if (md_file_read(file, num_atoms, sizeof(num_atoms)) != sizeof(num_atoms) || num_atoms == 0) {
			MD_LOG_ERROR("Failed to read offset cache, number of atoms was zero or corrupt");
			goto done;
		}

		if (md_file_read(file, cell, sizeof(md_unit_cell_t)) != sizeof(md_unit_cell_t)) {
			MD_LOG_ERROR("Failed to read offset cache, cell was corrupt");
			goto done;
		}

		if (md_file_read(file, &num_offsets, sizeof(num_offsets)) != sizeof(num_offsets) || num_offsets == 0) {
			MD_LOG_ERROR("Failed to read offset cache, number of frames was zero or corrupted");
			goto done;
		}

		md_array_resize(*offsets, num_offsets, alloc);

		const int64_t offset_bytes = md_array_bytes(*offsets);
		if (md_file_read(file, *offsets, offset_bytes) != offset_bytes) {
			MD_LOG_ERROR("Failed to read offset cache, offsets are incomplete");
			md_array_free(*offsets, alloc);
			goto done;
		}

		result = true;
	done:
		md_file_close(file);
	}
	return result;
}

static bool write_cache(str_t cache_file, const md_array(int64_t) offsets, int64_t num_offsets, int64_t num_atoms, const md_unit_cell_t* cell) {
	bool result = false;
	md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
	if (file) {
		const uint32_t magic = MD_LAMMPS_CACHE_MAGIC;
		const uint32_t version = MD_LAMMPS_CACHE_VERSION;
		const int64_t offset_bytes = num_offsets * sizeof(int64_t);

		if (md_file_write(file, &magic, sizeof(uint32_t)) != sizeof(uint32_t) ||
			md_file_write(file, &version, sizeof(uint32_t)) != sizeof(uint32_t) ||
			md_file_write(file, &num_atoms, sizeof(int64_t)) != sizeof(int64_t) ||
			md_file_write(file, cell, sizeof(md_unit_cell_t)) != sizeof(md_unit_cell_t) ||
			md_file_write(file, &num_offsets, sizeof(int64_t)) != sizeof(int64_t) ||
			md_file_write(file, offsets, offset_bytes) != offset_bytes)
		{
			MD_LOG_ERROR("Failed to write lammps cache");
			goto done;
		}

		result = true;

	done:
		md_file_close(file);
	}
	else {
		MD_LOG_ERROR("Failed to write offset cache, could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
	}

	return result;
}

void lammps_trajectory_free(struct md_trajectory_o* inst) {
	ASSERT(inst);
	md_lammps_trajectory_t* lammps_traj = (md_lammps_trajectory_t*)inst;
	if (lammps_traj->file) md_file_close(lammps_traj->file);
	if (lammps_traj->frame_offsets) md_array_free(lammps_traj->frame_offsets, lammps_traj->allocator);
	if (lammps_traj->header.frame_times) md_array_free(lammps_traj->header.frame_times, lammps_traj->allocator);
	md_mutex_destroy(&lammps_traj->mutex);
}

md_trajectory_i* md_lammps_trajectory_create(str_t filename, struct md_allocator_i* alloc) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file for LAMMPS trajectory");
		return false;
	}

	int64_t filesize = md_file_size(file);
	md_file_close(file);

	char buf[1024] = "";
	int len = snprintf(buf, sizeof(buf), "%.*s.cache", (int)filename.len, filename.ptr);
	str_t cache_file = { buf, len };

	md_unit_cell_t cell = { 0 };
	int64_t num_atoms = 0;
	int64_t* offsets = 0;

	//At this point, things probably changes as we are reading from a specific lammps trajectory data file

	if (!try_read_cache(cache_file, &offsets, &num_atoms, &cell, alloc)) { //If the cache file does not already exist
		md_lammps_trajectory_t data = { 0 };
		if (!md_lammps_trajectory_parse_file(&data, filename, md_heap_allocator)) {
			//We could not parse the data file
			return false;
		}

		for (int64_t i = 0; i < data.num_models; ++i) {
			md_array_push(offsets, data.models[i].byte_offset, alloc);
		}
		md_array_push(offsets, filesize, alloc);

		if (data.num_cryst1 > 0) {
			if (data.num_cryst1 > 1) {
				md_log(MD_LOG_TYPE_INFO, "The PDB file contains multiple CRYST1 entries, will pick the first one for determining the simulation box");
			}
			// If it is in fact a box, that will be handled as well
			cell = md_util_unit_cell_from_extent_and_angles(data.cryst1[0].a, data.cryst1[0].b, data.cryst1[0].c, data.cryst1[0].alpha, data.cryst1[0].beta, data.cryst1[0].gamma);
		}

		write_cache(cache_file, offsets, md_array_size(offsets), num_atoms, &cell);
		md_pdb_data_free(&data, md_heap_allocator);

	}

	const int64_t num_frames = md_array_size(offsets) - 1;

	int64_t max_frame_size = 0;
	for (int64_t i = 0; i < num_frames; ++i) {
		const int64_t beg = offsets[i + 0];
		const int64_t end = offsets[i + 1];
		const int64_t frame_size = end - beg;
		max_frame_size = MAX(max_frame_size, frame_size);
	}

	md_array(double) frame_times = md_array_create(double, num_frames, alloc);
	for (int64_t i = 0; i < num_frames; ++i) {
		frame_times[i] = (double)i;
	}

	void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(pdb_trajectory_t));
	ASSERT(mem);
	MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(pdb_trajectory_t));

	md_trajectory_i* traj = mem;
	pdb_trajectory_t* pdb = (pdb_trajectory_t*)(traj + 1);

	pdb->magic = MD_PDB_TRAJ_MAGIC;
	pdb->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	pdb->filesize = filesize;
	pdb->frame_offsets = offsets;
	pdb->allocator = alloc;
	pdb->mutex = md_mutex_create();
	pdb->header = (md_trajectory_header_t){
		.num_frames = md_array_size(offsets) - 1,
		.num_atoms = num_atoms,
		.max_frame_data_size = max_frame_size,
		.time_unit = {0},
		.frame_times = frame_times,
	};
	pdb->unit_cell = cell;

	traj->inst = (struct md_trajectory_o*)pdb;
	traj->get_header = pdb_get_header;
	traj->load_frame = pdb_load_frame;
	traj->fetch_frame_data = pdb_fetch_frame_data;
	traj->decode_frame_data = pdb_decode_frame_data;

	return traj;
}