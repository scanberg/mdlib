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

static size_t interpret_format(int mappings[ARRAY_SIZE(field_names)], const char* atom_format) {
	MEMSET(mappings, -1, sizeof(int) * ARRAY_SIZE(field_names));

	str_t str = str_trim(str_from_cstr(atom_format));
	str_t tokens[16];
	const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &str);
	for (size_t i = 0; i < num_tokens; ++i) {
		const str_t token = tokens[i];
		for (size_t j = 0; j < ARRAY_SIZE(field_names); ++j) {
			if (str_eq(token, str_from_cstr(field_names[j]))) {
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
	for (size_t i = 0; i < ARRAY_SIZE(required_fields); ++i) {
		int idx = required_fields[i];
		if (mappings[idx] == -1) {
			MD_LOG_ERROR("Lammps atom format: Missing required identifier '%s'", field_names[idx]);
			return 0;
		}
	}

	return num_tokens;
}

static bool parse_atoms(md_lammps_atom_t out_atoms[], size_t num_atoms, md_buffered_reader_t* reader, const int mappings[]) {
	ASSERT(mappings[0] != -1);
	ASSERT(mappings[2] != -1);
	ASSERT(mappings[4] != -1);
	ASSERT(mappings[5] != -1);
	ASSERT(mappings[6] != -1);

	str_t tok[16];
	str_t line;
	size_t read_atoms = 0;
	while (read_atoms < num_atoms && md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
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

static bool parse_bonds(md_lammps_bond_t out_bonds[], size_t bond_cap, md_buffered_reader_t* reader) {
	str_t tok[4];
	str_t line;
	size_t num_bonds = 0;
	while (num_bonds < bond_cap && md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
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

static bool parse_angles(md_lammps_angle_t out_angles[], size_t angle_cap, md_buffered_reader_t* reader) {
	str_t tok[8];
	str_t line;
	size_t num_angles = 0;
	while (num_angles < angle_cap && md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
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

static bool parse_dihedrals(md_lammps_dihedral_t out_dihedrals[], size_t dihedral_cap, md_buffered_reader_t* reader) {
	str_t tok[8];
	str_t line;
	size_t num_dihedrals = 0;
	while (num_dihedrals < dihedral_cap && md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
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

static bool parse_masses(md_array(float)* mass_type_table, size_t num_atom_types, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[4];
	str_t line;
	size_t num_mass = 0;
	while (num_mass < num_atom_types && md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 2) {
			MD_LOG_ERROR("Failed to parse mass line, expected 2 tokens, got %i", (int)num_tok);
			return false;
		}
		int   type = (int)parse_int(tok[0]);
		float mass = (float)parse_float(tok[1]);
		if (type >= (int)md_array_size(*mass_type_table)) {
			md_array_resize(*mass_type_table, (size_t)type, alloc);
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
	size_t num_fields = interpret_format(mappings, format);
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
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok > 0 && str_eq(tok[0], STR("Atoms"))) {
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
				for (size_t i = 0; i < data->num_atoms; ++i) {
					int32_t type = data->atoms[i].type;
					data->atoms[i].mass = type < (int)md_array_size(mass_table) ? mass_table[data->atoms[i].type] : 0.0f;
				}
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR("Bonds"))) {
			if (!data->num_bonds) {
				MD_LOG_ERROR("Encountered Bond entries, but number of bonds were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->bonds, data->num_bonds, alloc);
			if (!parse_bonds(data->bonds, data->num_bonds, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR("Angles"))) {
			if (!data->num_angles) {
				MD_LOG_ERROR("Encountered Angle entries, but number of angles were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR("Dihedrals"))) {
			if (!data->num_dihedrals) {
				MD_LOG_ERROR("Encountered Dihedral entries, but number of dihedrals were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR("Masses"))) {
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
			if (str_eq(tok[1], STR("atoms"))) {
				data->num_atoms = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("bonds"))) {
				data->num_bonds = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("angles"))) {
				data->num_angles = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("dihedrals"))) {
				data->num_dihedrals = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("impropers"))) {
				data->num_impropers = (int32_t)parse_int(tok[0]);
			}
		} else if (num_tok == 3 && str_eq(tok[2], STR("types")) && is_int(tok[0])) {
			if (str_eq(tok[1], STR("atom"))) {
				data->num_atom_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("bond"))) {
				data->num_bond_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("angle"))) {
				data->num_angle_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("dihedral"))) {
				data->num_dihedral_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR("improper"))) {
				data->num_improper_types = (int32_t)parse_int(tok[0]);
			} 
		} else if (num_tok == 4 && str_eq(tok[2], STR("xlo")) && str_eq(tok[3], STR("xhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.xlo = (float)parse_float(tok[0]);
			data->cell.xhi = (float)parse_float(tok[1]);
		} else if (num_tok == 4 && str_eq(tok[2], STR("ylo")) && str_eq(tok[3], STR("yhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.ylo = (float)parse_float(tok[0]);
			data->cell.yhi = (float)parse_float(tok[1]);
		} else if (num_tok == 4 && str_eq(tok[2], STR("zlo")) && str_eq(tok[3], STR("zhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.zlo = (float)parse_float(tok[0]);
			data->cell.zhi = (float)parse_float(tok[1]);
		} else if (num_tok == 6 && str_eq(tok[3], STR("xy")) && str_eq(tok[4], STR("xz")) && str_eq(tok[5], STR("yz"))) {
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
		const size_t cap = MEGABYTES(1);
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

	const size_t num_atoms = data->num_atoms;

	md_array_resize(mol->atom.x,	 num_atoms, alloc);
	md_array_resize(mol->atom.y,	 num_atoms, alloc);
	md_array_resize(mol->atom.z,	 num_atoms, alloc);
	md_array_resize(mol->atom.flags, num_atoms, alloc);
	md_array_resize(mol->atom.resid, num_atoms, alloc);
	md_array_resize(mol->atom.mass,  num_atoms, alloc);

	int32_t prev_res_id = -1;
	for (size_t i = 0; i < num_atoms; ++i) {
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