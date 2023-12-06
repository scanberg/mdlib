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

static const char* field_names[] = {
	"id",
	"resid",
	"type",
	"charge",
	"x",
	"y",
	"z"
};

static int required_fields[] = {0,2,4,5,6};

typedef struct mass_entry_t {
	int32_t type;
	float   mass;
} mass_entry_t;

static int interpret_format(int mappings[ARRAY_SIZE(field_names)], const char* atom_format) {
	MEMSET(mappings, -1, sizeof(int) * ARRAY_SIZE(field_names));

	str_t str = str_trim(str_from_cstr(atom_format));
	str_t tokens[16];
	const int64_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &str);
	for (int64_t i = 0; i < num_tokens; ++i) {
		const str_t token = tokens[i];
		for (int64_t j = 0; j < ARRAY_SIZE(field_names); ++j) {
			if (str_equal(token, str_from_cstr(field_names[j]))) {
				if (mappings[j] != -1) {
					MD_LOG_ERROR("Lammps atom format: Contains duplicate field name '%.*s'", token.len, token.ptr);
					return 0;
				}
				mappings[j] = (int)i;
				break;
			}
		}
	}

	// Check required mappings
	for (int i = 0; i < ARRAY_SIZE(required_fields); ++i) {
		int idx = required_fields[i];
		if (mappings[idx] == -1) {
			MD_LOG_ERROR("Lammps atom format: Missing required identifier '%s'", field_names[idx]);
			return 0;
		}
	}

	return (int)num_tokens;
}

static bool parse_atoms(md_lammps_atom_t out_atoms[], int32_t num_atoms, md_buffered_reader_t* reader, const int mappings[]) {
	ASSERT(mappings[0] != -1);
	ASSERT(mappings[2] != -1);
	ASSERT(mappings[4] != -1);
	ASSERT(mappings[5] != -1);
	ASSERT(mappings[6] != -1);

	str_t tok[16];
	str_t line;
	int32_t read_atoms = 0;
	while (read_atoms < num_atoms && md_buffered_reader_extract_line(&line, reader)) {
		int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 5) {
			MD_LOG_ERROR("Failed to parse atom line, expected at least 5 tokens, got %i", (int)num_tok);
			return false;
		}

		md_lammps_atom_t* atom = &out_atoms[read_atoms++];
		atom->id = (int32_t)parse_int(tok[mappings[0]]);
		atom->resid = mappings[1] != -1 ? (int32_t)parse_int(tok[mappings[1]]) : -1;
		atom->type = (int32_t)parse_int(tok[mappings[2]]);
		atom->charge = mappings[3] != -1 ? (float)parse_float(tok[mappings[3]]) : 0.0f;
		atom->x = (float)parse_float(tok[mappings[4]]);
		atom->y = (float)parse_float(tok[mappings[5]]);
		atom->z = (float)parse_float(tok[mappings[6]]);
	}

	return true;
}

static bool parse_bonds(md_lammps_bond_t out_bonds[], int32_t bond_cap, md_buffered_reader_t* reader) {
	str_t tok[4];
	str_t line;
	int32_t num_bonds = 0;
	while (num_bonds < bond_cap && md_buffered_reader_extract_line(&line, reader)) {
		int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 4) {
			MD_LOG_ERROR("Failed to parse bond line, expected 4 tokens, got %i", (int)num_tok);
			return false;
		}

		md_lammps_bond_t* bond = &out_bonds[num_bonds++];
		bond->id = (int32_t)parse_int(tok[0]);
		bond->type = (int32_t)parse_int(tok[1]);
		bond->atom_id[0] = (int32_t)parse_int(tok[2]);
		bond->atom_id[1] = (int32_t)parse_int(tok[3]);
	}

	return true;
}

static bool parse_angles(md_lammps_angle_t out_angles[], int32_t angle_cap, md_buffered_reader_t* reader) {
	str_t tok[8];
	str_t line;
	int32_t num_angles = 0;
	while (num_angles < angle_cap && md_buffered_reader_extract_line(&line, reader)) {
		int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 5) {
			MD_LOG_ERROR("Failed to parse angle line, expected 5 tokens, got %i", (int)num_tok);
			return false;
		}
	
		md_lammps_angle_t* angle = &out_angles[num_angles++];
		angle->id = (int32_t)parse_int(tok[0]);
		angle->type = (int32_t)parse_int(tok[1]);
		angle->atom_id[0] = (int32_t)parse_int(tok[2]);
		angle->atom_id[1] = (int32_t)parse_int(tok[3]);
		angle->atom_id[2] = (int32_t)parse_int(tok[4]);
	}

	return true;
}

static bool parse_dihedrals(md_lammps_dihedral_t out_dihedrals[], int32_t dihedral_cap, md_buffered_reader_t* reader) {
	str_t tok[8];
	str_t line;
	int32_t num_dihedrals = 0;
	while (num_dihedrals < dihedral_cap && md_buffered_reader_extract_line(&line, reader)) {
		int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 6) {
			MD_LOG_ERROR("Failed to parse dihedral line, expected 6 tokens, got %i", (int)num_tok);
			return false;
		}
	
		md_lammps_dihedral_t* dihedral = &out_dihedrals[num_dihedrals++];
		dihedral->id = (int32_t)parse_int(tok[0]);
		dihedral->type = (int32_t)parse_int(tok[1]);
		dihedral->atom_id[0] = (int32_t)parse_int(tok[2]);
		dihedral->atom_id[1] = (int32_t)parse_int(tok[3]);
		dihedral->atom_id[2] = (int32_t)parse_int(tok[4]);
		dihedral->atom_id[3] = (int32_t)parse_int(tok[5]);
	}

	return true;
}

static bool parse_masses(md_array(float)* mass_type_table, int32_t num_atom_types, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[4];
	str_t line;
	int32_t num_mass = 0;
	while (num_mass < num_atom_types && md_buffered_reader_extract_line(&line, reader)) {
		int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 2) {
			MD_LOG_ERROR("Failed to parse mass line, expected 2 tokens, got %i", (int)num_tok);
			return false;
		}
		int   type = (int)parse_int(tok[0]);
		float mass = (float)parse_float(tok[1]);
		if (type >= md_array_size(*mass_type_table)) {
			md_array_resize(*mass_type_table, type, alloc);
		}
		(*mass_type_table)[type] = mass;
		num_mass += 1;
	}
	return true;
}

static bool md_lammps_data_parse(md_lammps_data_t* data, md_buffered_reader_t* reader, const char* format, struct md_allocator_i* alloc) {
	ASSERT(data);
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tok[16];

	int mappings[ARRAY_SIZE(field_names)];
	int num_fields = interpret_format(mappings, format);
	if (num_fields == 0) {
		return false;
	}

	// Read the title of the file
	if (!md_buffered_reader_extract_line(&line, reader)) {
		MD_LOG_ERROR("Failed to parse lammps title");
		return false;
	}

	MEMSET(data, 0, sizeof(md_lammps_data_t));

	md_array(float) mass_table = 0;

	str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

	// Parse headers and sections
	while (md_buffered_reader_extract_line(&line, reader)) {
		const int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok > 0 && str_equal(tok[0], STR("Atoms"))) {
			if (!data->num_atoms) {
				MD_LOG_ERROR("Encountered Atom entries, but number of atoms were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->atoms, data->num_atoms, alloc);
			if (!parse_atoms(data->atoms, data->num_atoms, reader, mappings)) {
				return false;
			}
			if (mass_table) {
				for (int64_t i = 0; i < data->num_atoms; ++i) {
					int32_t type = data->atoms[i].type;
					data->atoms[i].mass = type < md_array_size(mass_table) ? mass_table[data->atoms[i].type] : 0.0f;
				}
			}
		} else if (num_tok > 0 && str_equal(tok[0], STR("Bonds"))) {
			if (!data->num_bonds) {
				MD_LOG_ERROR("Encountered Bond entries, but number of bonds were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->bonds, data->num_bonds, alloc);
			if (!parse_bonds(data->bonds, data->num_bonds, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_equal(tok[0], STR("Angles"))) {
			if (!data->num_angles) {
				MD_LOG_ERROR("Encountered Angle entries, but number of angles were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_equal(tok[0], STR("Dihedrals"))) {
			if (!data->num_dihedrals) {
				MD_LOG_ERROR("Encountered Dihedral entries, but number of dihedrals were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_equal(tok[0], STR("Masses"))) {
			if (!data->num_atom_types) {
				MD_LOG_ERROR("Encountered Mass entries, but number of atom types were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			mass_table = md_array_create(float, data->num_atom_types + 1, md_temp_allocator);
			MEMSET(mass_table, 0, md_array_bytes(mass_table));
			if (!parse_masses(&mass_table, data->num_atom_types, reader, md_temp_allocator)) {
				return false;
			}
		} else if (num_tok == 2 && is_int(tok[0])) {
			if (str_equal(tok[1], STR("atoms"))) {
				data->num_atoms = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("bonds"))) {
				data->num_bonds = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("angles"))) {
				data->num_angles = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("dihedrals"))) {
				data->num_dihedrals = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("impropers"))) {
				data->num_impropers = (int32_t)parse_int(tok[0]);
			}
		} else if (num_tok == 3 && str_equal(tok[2], STR("types")) && is_int(tok[0])) {
			if (str_equal(tok[1], STR("atom"))) {
				data->num_atom_types = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("bond"))) {
				data->num_bond_types = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("angle"))) {
				data->num_angle_types = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("dihedral"))) {
				data->num_dihedral_types = (int32_t)parse_int(tok[0]);
			} else if (str_equal(tok[1], STR("improper"))) {
				data->num_improper_types = (int32_t)parse_int(tok[0]);
			} 
		} else if (num_tok == 4 && str_equal(tok[2], STR("xlo")) && str_equal(tok[3], STR("xhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.xlo = (float)parse_float(tok[0]);
			data->cell.xhi = (float)parse_float(tok[1]);
		} else if (num_tok == 4 && str_equal(tok[2], STR("ylo")) && str_equal(tok[3], STR("yhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.ylo = (float)parse_float(tok[0]);
			data->cell.yhi = (float)parse_float(tok[1]);
		} else if (num_tok == 4 && str_equal(tok[2], STR("zlo")) && str_equal(tok[3], STR("zhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.zlo = (float)parse_float(tok[0]);
			data->cell.zhi = (float)parse_float(tok[1]);
		} else if (num_tok == 6 && str_equal(tok[3], STR("xy")) && str_equal(tok[4], STR("xz")) && str_equal(tok[5], STR("yz"))) {
			if (!is_float(tok[0]) || !is_float(tok[1]) || !is_float(tok[2])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.xy = (float)parse_float(tok[0]);
			data->cell.xz = (float)parse_float(tok[1]);
			data->cell.yz = (float)parse_float(tok[2]);
		}
	}

	return true;
}

bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, const char* format, struct md_allocator_i* alloc) {
	ASSERT(data);
	ASSERT(alloc);

	md_buffered_reader_t line_reader = md_buffered_reader_from_str(str);
	return md_lammps_data_parse(data, &line_reader, format, alloc);
}

bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, const char* format, struct md_allocator_i* alloc) {
	bool result = false;
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (file) {
		const int64_t cap = MEGABYTES(1);
		char* buf = md_alloc(md_heap_allocator, cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = md_lammps_data_parse(data, &line_reader, format, alloc);

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
	if (data->atoms) md_array_free(data->atoms, alloc);
	if (data->bonds) md_array_free(data->bonds, alloc);
	if (data->angles) md_array_free(data->angles, alloc);
	if (data->dihedrals) md_array_free(data->dihedrals, alloc);
	if (data->impropers) md_array_free(data->impropers, alloc);
	MEMSET(data, 0, sizeof(md_lammps_data_t));
}

bool md_lammps_molecule_init(md_molecule_t* mol, const md_lammps_data_t* data, md_allocator_i* alloc)
{
	ASSERT(mol);
	ASSERT(data);
	ASSERT(alloc);

	MEMSET(mol, 0, sizeof(md_molecule_t));

	const int64_t num_atoms = data->num_atoms;

	md_array_resize(mol->atom.x,	 num_atoms, alloc);
	md_array_resize(mol->atom.y,	 num_atoms, alloc);
	md_array_resize(mol->atom.z,	 num_atoms, alloc);
	md_array_resize(mol->atom.flags, num_atoms, alloc);
	md_array_resize(mol->atom.resid, num_atoms, alloc);
	md_array_resize(mol->atom.mass,  num_atoms, alloc);

	int32_t prev_res_id = -1;
	for (int64_t i = 0; i < num_atoms; ++i) {
		mol->atom.x[i] = data->atoms[i].x;
		mol->atom.y[i] = data->atoms[i].y;
		mol->atom.z[i] = data->atoms[i].z;
		mol->atom.resid[i] = data->atoms[i].resid;
		mol->atom.mass[i] = data->atoms[i].mass;
		mol->atom.flags[i] = 0;

		if (prev_res_id != mol->atom.resid[i]) {
			mol->atom.flags[i] |= MD_FLAG_RES_BEG;

			if (prev_res_id != -1) {
				*md_array_last(mol->atom.flags) |= MD_FLAG_RES_END;
			}
			prev_res_id = mol->atom.resid[i];
		}
	}

	mol->atom.count = num_atoms;

	//Set elements
	md_array_resize(mol->atom.element, num_atoms, alloc);
	if (!md_util_element_from_mass(mol->atom.element, mol->atom.mass, num_atoms)) {
		MD_LOG_ERROR("One or more masses are missing matching element");
	}

	//Create unit cell
	float M[3][3];
	M[0][0] = data->cell.xhi - data->cell.xlo;
	M[1][1] = data->cell.yhi - data->cell.ylo;
	M[2][2] = data->cell.zhi - data->cell.zlo;
	M[1][0] = data->cell.xy;
	M[2][0] = data->cell.xz;
	M[2][1] = data->cell.yz;
	mol->unit_cell = md_util_unit_cell_from_matrix(M);
	
	return true;
}

bool md_lammps_init_from_str(md_molecule_t* mol, str_t str, md_allocator_i* alloc, const char* format) {
	md_lammps_data_t data = { 0 };
	bool success = false;
	if (md_lammps_data_parse_str(&data, str, format, md_heap_allocator)) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_heap_allocator);

	return success;
}

bool md_lammps_init_from_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc, const char* format) {
	md_lammps_data_t data = { 0 };
	bool success = false;
	if (md_lammps_data_parse_file(&data, filename, format, md_heap_allocator)) {
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

typedef struct lammps_trajectory_t {
	int64_t* frame_offsets;
	int64_t num_frame_offsets;

	uint64_t magic;
	md_file_o* file;
	uint64_t filesize;
	md_unit_cell_t unit_cell;
	md_trajectory_header_t header;
	md_allocator_i* allocator;
	md_mutex_t mutex;
} lammps_trajectory_t;

//Reads data that is useful later when we want to parse a frame from the trajectory
bool lammps_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
	lammps_trajectory_t* dataPtr = (lammps_trajectory_t*)inst;
	ASSERT(dataPtr);
	ASSERT(dataPtr->magic == MD_LAMMPS_TRAJ_MAGIC);
	ASSERT(header);

	*header = dataPtr->header;
	return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
// Returns size in bytes of frame, frame_data_ptr is optional and is the destination to write the frame data to.
int64_t lammps_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
	lammps_trajectory_t* traj_data = (lammps_trajectory_t*)inst;
	ASSERT(traj_data);
	ASSERT(traj_data->magic == MD_LAMMPS_TRAJ_MAGIC);

	if (!traj_data->filesize) {
		MD_LOG_ERROR("File size is zero");
		return 0;
	}

	if (!traj_data->file) {
		MD_LOG_ERROR("File handle is NULL");
		return 0;
	}

	if (!traj_data->frame_offsets) {
		MD_LOG_ERROR("Frame offsets is empty");
		return 0;
	}

	if (!(0 <= frame_idx && frame_idx < (int64_t)md_array_size(traj_data->frame_offsets) - 1)) {
		MD_LOG_ERROR("Frame index is out of range");
		return 0;
	}

	const int64_t beg = traj_data->frame_offsets[frame_idx + 0];
	const int64_t end = traj_data->frame_offsets[frame_idx + 1];
	const int64_t frame_size = end - beg;

	if (frame_data_ptr) {
		// Store the index to the frame since this is generally not found within the actual frame data
		int64_t* ptr = (int64_t*)frame_data_ptr;
		ptr[0] = frame_idx;
		ASSERT(traj_data->file);
		md_mutex_lock(&traj_data->mutex);
		md_file_seek(traj_data->file, beg, MD_FILE_BEG);
		const int64_t bytes_read = md_file_read(traj_data->file, &ptr[1], frame_size);
		md_mutex_unlock(&traj_data->mutex);
		ASSERT(frame_size == bytes_read);
	}

	return frame_size;
}

bool lammps_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
	ASSERT(inst);
	ASSERT(frame_data_ptr);
	ASSERT(frame_data_size);

	str_t tokens[32];
	int64_t num_tokens = 0;
	int64_t num_token_types = 0;

	const int64_t* ptr = frame_data_ptr;
	int32_t frame_idx = (int32_t)ptr[0];


	lammps_trajectory_t* traj_data = (lammps_trajectory_t*)inst;
	if (traj_data->magic != MD_LAMMPS_TRAJ_MAGIC) {
		MD_LOG_ERROR("Error when decoding frame header, lammps magic did not match");
		return false;
	}
	if (frame_idx < 0 || frame_idx >= traj_data->header.num_frames) {
		MD_LOG_ERROR("Error when decoding frame data, corrupt frame index");
		return false;
	}

	str_t str = { .ptr = (const char*)(ptr + 1), .len = frame_data_size };
	str_t line;

	//Parse the ITEM: ATOMS to see what type if coordinates we are working with
	str_t coord_type = { 0 };
	int32_t coord_start = -1;
	if (str_extract_line(&line, &str)) {
		if (str_equal_cstr_n(line, "ITEM: ATOMS", 11)) {
			num_tokens = extract_tokens(tokens, 32, &line);
			num_token_types = num_tokens - 2;
			if (num_token_types < 3) {
				MD_LOG_ERROR("Fewer than 3 token types in ITEM: ATOMS line");
				return false;
			}
			for (int32_t i = 0; i < num_tokens; i++) {
				//Use the first complete coord combo. x,y,z,xs,ys picks x,y,z and x,y,xs,ys,zs picks xs,ys,zs
				if (
					str_equal_cstr(tokens[i], "x") && str_equal_cstr(tokens[i + 1], "y") && str_equal_cstr(tokens[i + 2], "z") ||
					str_equal_cstr(tokens[i], "xs") && str_equal_cstr(tokens[i + 1], "ys") && str_equal_cstr(tokens[i + 2], "zs") ||
					str_equal_cstr(tokens[i], "xu") && str_equal_cstr(tokens[i + 1], "yu") && str_equal_cstr(tokens[i + 2], "zu") ||
					str_equal_cstr(tokens[i], "xs") && str_equal_cstr(tokens[i + 1], "ys") && str_equal_cstr(tokens[i + 2], "zs")) {
					coord_type = tokens[i];
					coord_start = i - 2; //Skip ITEM: ATOMS tokens
					break;
				}
			}
			if (coord_start == -1) {
				MD_LOG_ERROR("Could not parse coord definition in ITEM: ATOMS line");
				return false;
			}
			//step = (int32_t)parse_int(str_substr(line, 10, 4));
		}
	}

	//Extract coordinates and scale them if the type requires it. We leave unwrapped coordinates as is
	int64_t i = 0;
	while (str_extract_line(&line, &str) && i < traj_data->header.num_atoms) {
		num_tokens = extract_tokens(tokens, 32, &line);
		if (coord_start + 2 + 1 > num_tokens) {
			MD_LOG_ERROR("Not enough tokens in atom line");
			return false;
		}
		if (str_equal_cstr(coord_type, "x") || str_equal_cstr(coord_type, "xu")) { //Unscaled coordinates
			x[i] = (float)parse_float(tokens[coord_start]);
			y[i] = (float)parse_float(tokens[coord_start + 1]);
			z[i] = (float)parse_float(tokens[coord_start + 2]);
		}
		else if (str_equal_cstr(coord_type, "xs") || str_equal_cstr(coord_type, "xsu")) { //Scaled coordinates
			mat3_t unit_cell = traj_data->unit_cell.basis;
			vec3_t scale_vec = {
				(float)parse_float(tokens[coord_start]),
				(float)parse_float(tokens[coord_start + 1]),
				(float)parse_float(tokens[coord_start + 2]),
			};
			vec3_t result_vec = mat3_mul_vec3(unit_cell, scale_vec);
			x[i] = result_vec.x;
			y[i] = result_vec.y;
			z[i] = result_vec.z;
		}
		i++;
	}

	if (header) {
		header->num_atoms = traj_data->header.num_atoms;
		header->index = frame_idx;
		header->timestamp = (double)frame_idx;
		header->unit_cell = traj_data->unit_cell;
	}

	return true;
}

static bool md_lammps_trajectory_parse(int64_t* num_atoms, md_unit_cell_t* unit_cell, int64_t* num_frame_offsets, md_array(int64_t)* frame_offsets, md_array(int64_t)* frame_times, md_buffered_reader_t* reader, struct md_allocator_i* alloc) {
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tokens[6];
	int64_t num_frames = 0;


	//Setup frames and first frame
	md_buffered_reader_skip_line(reader);
	md_buffered_reader_extract_line(&line, reader);
	md_array_push(*frame_times, parse_int(line), alloc);
	num_frames++;


	//Parse num of atoms
	{
		md_buffered_reader_skip_line(reader);
		md_buffered_reader_extract_line(&line, reader); //Read num of atoms
		*num_atoms = parse_int(line);
	}

	//Parse unit_cell definition
	{
		float cell_extent[3];
		float cell_tri[3]; //Triclinic

		md_buffered_reader_extract_line(&line, reader); //Read BOX BOUNDS
		if (str_equal_cstr(line, "ITEM: BOX BOUNDS pp pp pp")) {
			//Cubic
			for (int8_t i = 0; i < 3; i++) {
				md_buffered_reader_extract_line(&line, reader);
				extract_tokens(tokens, 4, &line);
				cell_extent[i] = (float)(parse_float(tokens[1]) - parse_float(tokens[0]));
			}
			*unit_cell = md_util_unit_cell_from_extent(cell_extent[0], cell_extent[1], cell_extent[2]);
		}
		else if (str_equal_cstr(line, "ITEM: BOX BOUNDS xy xz yz pp pp pp")) {
			//Triclinic
			for (int8_t i = 0; i < 3; i++) {
				md_buffered_reader_extract_line(&line, reader);
				extract_tokens(tokens, 4, &line);
				cell_extent[i] = (float)(parse_float(tokens[1]) - parse_float(tokens[0]));
				cell_tri[i] = (float)parse_float(tokens[2]);
			}
			*unit_cell = md_util_unit_cell_from_triclinic(cell_extent[0], cell_extent[1], cell_extent[2], cell_tri[0], cell_tri[1], cell_tri[2]);
		}
		else {
			MD_LOG_ERROR("Could not correctly parse BOX BOUND");
			return false;
		}
	}

	while (md_buffered_reader_extract_line(&line, reader)) {

		if (str_equal_cstr_n(line, "ITEM: TIMESTEP", 14)) {
			md_buffered_reader_extract_line(&line, reader);
			num_frames++;
			md_array_push(*frame_times, parse_int(line), alloc);
		}
		else if (str_equal_cstr_n(line, "ITEM: ATOMS", 11)) {

			// This is a bit nasty, we want to get the correct offset to the beginning of the current line.
			// Therefore we need to do some pointer arithmetic because just using the length of the line may not get us
			// all the way back in case there were skipped \r characters.
			const int64_t offset = md_buffered_reader_tellg(reader) - (reader->str.ptr - line.ptr);
			md_array_push(*frame_offsets, offset, alloc);
		}
	}

	//We add the end of the file to frame_offsets, so frame_offset size = num_frames + 1

	const int64_t end_of_file = md_buffered_reader_tellg(reader);
	md_array_push(*frame_offsets, end_of_file, alloc);
	*num_frame_offsets = (int64_t)md_array_size(*frame_offsets);

	return true;
}

bool md_lammps_trajectory_parse_file(int64_t* num_atoms, md_unit_cell_t* unit_cell, int64_t* num_frame_offsets, md_array(int64_t)* frame_offsets, md_array(int64_t)* frame_times, str_t filename, struct md_allocator_i* alloc) {
	bool result = false;
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (file) {
		const int64_t cap = MEGABYTES(1);
		char* buf = md_alloc(md_heap_allocator, cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = md_lammps_trajectory_parse(num_atoms, unit_cell, num_frame_offsets, frame_offsets, frame_times, &line_reader, alloc);

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

	lammps_trajectory_t* lammps_traj = (lammps_trajectory_t*)inst;
	if (lammps_traj->magic != MD_LAMMPS_TRAJ_MAGIC) {
		MD_LOG_ERROR("Error when decoding frame coord, xtc magic did not match");
		return false;
	}

	// Should this be exposed?
	md_allocator_i* alloc = md_temp_allocator;

	bool result = true;
	const int64_t frame_size = lammps_fetch_frame_data(inst, frame_idx, NULL);
	if (frame_size > 0) { //This check first that there actually is data to be read, before writing to frame_data
		// This is a borderline case if one should use the md_temp_allocator as the raw frame size could potentially be several megabytes...
		void* frame_data = md_alloc(alloc, frame_size);
		const int64_t read_size = lammps_fetch_frame_data(inst, frame_idx, frame_data);
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

static bool try_read_cache(str_t cache_file, int64_t* num_atoms, md_unit_cell_t* unit_cell, int64_t* num_frame_offsets, md_array(int64_t)* frame_offsets, md_array(int64_t)* frame_times, int64_t filesize, md_allocator_i* alloc) {
	md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
	bool result = false;
	if (file) {
		uint32_t version = 0;
		uint32_t magic = 0;
		int64_t loc_num_atoms = 0;
		md_unit_cell_t loc_unit_cell = { 0 };
		int64_t num_offsets = 0;
		md_array(int64_t) loc_frame_offsets = 0;
		md_array(int64_t) loc_frame_times = 0;


		int64_t loc_filesize = 0;


		if (md_file_read(file, &magic, sizeof(magic)) != sizeof(magic) || magic != MD_LAMMPS_CACHE_MAGIC) {
			MD_LOG_ERROR("Failed to read offset cache, magic was incorrect or corrupt");
			goto done;
		}

		if (md_file_read(file, &version, sizeof(version)) != sizeof(version) || version != MD_LAMMPS_CACHE_VERSION) {
			MD_LOG_ERROR("Failed to read offset cache, version was incorrect");
			goto done;
		}

		if (md_file_read(file, &loc_num_atoms, sizeof(loc_num_atoms)) != sizeof(loc_num_atoms) || loc_num_atoms == 0) {
			MD_LOG_ERROR("Failed to read offset cache, number of atoms was zero or corrupt");
			goto done;
		}

		if (md_file_read(file, &loc_unit_cell, sizeof(md_unit_cell_t)) != sizeof(md_unit_cell_t)) {
			MD_LOG_ERROR("Failed to read offset cache, cell was corrupt");
			goto done;
		}

		if (md_file_read(file, &num_offsets, sizeof(num_offsets)) != sizeof(num_offsets) || num_offsets == 0) {
			MD_LOG_ERROR("Failed to read offset cache, number of frames was zero or corrupted");
			goto done;
		}

		md_array_resize(loc_frame_offsets, num_offsets, alloc);
		const int64_t offset_bytes = md_array_bytes(loc_frame_offsets);
		if (md_file_read(file, loc_frame_offsets, offset_bytes) != offset_bytes) {
			MD_LOG_ERROR("Failed to read offset cache, offsets are incomplete");
			md_array_free(loc_frame_offsets, alloc);
			goto done;
		}

		int64_t num_frames = num_offsets - 1;
		md_array_resize(loc_frame_times, num_frames, alloc);
		const int64_t frame_times_bytes = md_array_bytes(loc_frame_times);
		if (md_file_read(file, loc_frame_times, frame_times_bytes) != frame_times_bytes) {
			MD_LOG_ERROR("Failed to read frame times cache, frame times are incomplete");
			md_array_free(loc_frame_times, alloc);
			goto done;
		}

		if (md_file_read(file, &loc_filesize, sizeof(loc_filesize)) != sizeof(loc_filesize) || loc_filesize != filesize) {
			MD_LOG_ERROR("Failed to read offset cache, filesize is not the same");
			goto done;
		}

		//Only overwrite if all tests passed
		*num_atoms = loc_num_atoms;
		*unit_cell = loc_unit_cell;
		*num_frame_offsets = num_offsets;
		*frame_offsets = loc_frame_offsets;
		*frame_times = loc_frame_times;

		result = true;
	done:
		md_file_close(file);
	}
	return result;
}

static bool write_cache(str_t cache_file, int64_t* num_atoms, md_unit_cell_t* unit_cell, md_array(int64_t)* frame_offsets, md_array(int64_t)* frame_times, int64_t* filesize) {
	bool result = false;
	md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
	if (file) {
		const uint32_t magic = MD_LAMMPS_CACHE_MAGIC;
		const uint32_t version = MD_LAMMPS_CACHE_VERSION;

		const int64_t num_frame_offsets = (int64_t)md_array_size(*frame_offsets);
		const int64_t offset_bytes = md_array_bytes(*frame_offsets);

		const int64_t frame_times_bytes = md_array_bytes(*frame_times);


		if (md_file_write(file, &magic, sizeof(uint32_t)) != sizeof(uint32_t) ||
			md_file_write(file, &version, sizeof(uint32_t)) != sizeof(uint32_t) ||
			md_file_write(file, num_atoms, sizeof(int64_t)) != sizeof(int64_t) ||
			md_file_write(file, unit_cell, sizeof(md_unit_cell_t)) != sizeof(md_unit_cell_t) ||
			md_file_write(file, &num_frame_offsets, sizeof(int64_t)) != sizeof(int64_t) ||
			md_file_write(file, *frame_offsets, offset_bytes) != offset_bytes ||
			md_file_write(file, *frame_times, frame_times_bytes) != frame_times_bytes ||
			md_file_write(file, filesize, sizeof(int64_t)) != sizeof(int64_t))
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

void lammps_trajectory_data_free(struct lammps_trajectory_t* lammps_traj) {
	ASSERT(lammps_traj);
	if (lammps_traj->file) md_file_close(lammps_traj->file);
	if (lammps_traj->frame_offsets) md_array_free(lammps_traj->frame_offsets, lammps_traj->allocator);
	if (lammps_traj->header.frame_times) md_array_free(lammps_traj->header.frame_times, lammps_traj->allocator);
	md_mutex_destroy(&lammps_traj->mutex);
}

void lammps_trajectory_free(struct md_trajectory_o* inst) {
	ASSERT(inst);
	lammps_trajectory_t* lammps_traj = (lammps_trajectory_t*)inst;
	lammps_trajectory_data_free(lammps_traj);
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

	int64_t num_atoms = 0;
	md_unit_cell_t cell = { 0 };
	int64_t num_frame_offsets = 0;
	md_array(int64_t) offsets = 0;
	md_array(int64_t) frame_times = 0;

	if (!try_read_cache(cache_file, &num_atoms, &cell, &num_frame_offsets, &offsets, &frame_times, filesize, alloc)) { //If the cache file does not exist, we create one
		if (!md_lammps_trajectory_parse_file(&num_atoms, &cell, &num_frame_offsets, &offsets, &frame_times, filename, md_heap_allocator)) {
			//We could not parse the data file
			md_array_free(offsets, alloc);
			md_array_free(frame_times, alloc);
			return false;
		}

		if (!write_cache(cache_file, &num_atoms, &cell, &offsets, &frame_times, &filesize)) {
			MD_LOG_ERROR("Could not write cache");
			md_array_free(offsets, alloc);
			md_array_free(frame_times, alloc);
			return false;
		}
	}

	int64_t num_frames = (int64_t)md_array_size(offsets) - 1;

	int64_t max_frame_size = 0;
	//Calculate the max frame size
	for (int64_t i = 0; i < num_frames; i++) {
		const int64_t beg = offsets[i + 0];
		const int64_t end = offsets[i + 1];
		const int64_t frame_size = end - beg;
		max_frame_size = MAX(max_frame_size, frame_size);
	}

	md_array(double) double_frame_times = md_array_create(double, num_frames, alloc);
	for (int64_t i = 0; i < num_frames; i++) {
		double_frame_times[i] = (double)frame_times[i];
	}

	void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(lammps_trajectory_t));
	ASSERT(mem);
	MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(lammps_trajectory_t));

	md_trajectory_i* traj = mem;
	lammps_trajectory_t* traj_data = (lammps_trajectory_t*)(traj + 1);

	//Can I set dataPtr to traj_data here?
	traj_data->magic = MD_LAMMPS_TRAJ_MAGIC;
	traj_data->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	traj_data->filesize = filesize;

	traj_data->frame_offsets = offsets;

	traj_data->allocator = alloc;
	traj_data->mutex = md_mutex_create();

	traj_data->header = (md_trajectory_header_t){
		.num_frames = num_frames,
		.num_atoms = num_atoms,
		.max_frame_data_size = max_frame_size,
		.time_unit = md_unit_femtosecond(),
		.frame_times = double_frame_times,
	};
	traj_data->unit_cell = cell;

	traj->inst = (struct md_trajectory_o*)traj_data;
	traj->get_header = lammps_get_header;
	traj->load_frame = lammps_load_frame;
	traj->fetch_frame_data = lammps_fetch_frame_data;
	traj->decode_frame_data = lammps_decode_frame_data;

	return traj;
}

void md_lammps_trajectory_free(md_trajectory_i* traj) {
	ASSERT(traj);
	ASSERT(traj->inst);
	lammps_trajectory_t* lammps_data = (lammps_trajectory_t*)traj->inst;
	if (lammps_data->magic != MD_LAMMPS_TRAJ_MAGIC) {
		MD_LOG_ERROR("Trajectory is not a valid LAMMPS trajectory.");
		ASSERT(false);
		return;
	}

	md_allocator_i* alloc = lammps_data->allocator;
	lammps_trajectory_free(traj->inst);
	md_free(alloc, traj, sizeof(md_trajectory_i) + sizeof(lammps_trajectory_t));
}

static md_trajectory_loader_i lammps_traj_loader = {
	md_lammps_trajectory_create,
	md_lammps_trajectory_free,
};

md_trajectory_loader_i* md_lammps_trajectory_loader() {
	return &lammps_traj_loader;
}