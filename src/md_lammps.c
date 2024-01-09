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

#define MD_LAMMPS_MOLECULE_LOADER_ARG_TYPE 0x341293abc8273650

static const char* atom_format_name[MD_LAMMPS_ATOM_FORMAT_COUNT] = {
	"unknown",
	"angle",
	"atomic",
	"body",
	"bond",
	"charge",
	"dipole",
	"dpd",
	"edpd",
	"mdpd",
	"electron",
	"ellipsoid",
	"full",
	"line",
	"meso",
	"molecular",
	"peri",
	"smd",
	"sphere",
	"template",
	"tri",
	"wavepacket",
};

static const char* atom_format_string[MD_LAMMPS_ATOM_FORMAT_COUNT] = {
	"", // unknown
	"id resid type x y z", // angle
	"id type x y z", // atomic
	"id type bodyflag mass x y z", // body
	"id resid type x y z", // bond
	"id type q x y z", // charge
	"id type q x y z", // dipole
	"id type theta x y z", // dpd
	"id type temp cv x y z", // edpd
	"id type rho x y z", // mdpd
	"id type q spin radius x y z", // electron
	"id type ellispoidflag density x y z", // ellipsoid
	"id resid type q x y z", // full
	"id resid type lineflag density x y z", // line
	"id type rho e cv x y z", // meso
	"id resid type x y z", // molecular
	"id type volume density x y z", // peri
	"id type molecule volume mass kernel-radius contact-radius x0 y0 z0 x y z", // smd
	"id type diameter density x y z", // sphere
	"id resid template-index template-atom type x y z", // template
	"id resid type triangleflag density x y z", // tri
	"id type q spin eradius etag cs_re cs_im x y z", // wavepacket
};

enum {
	TYPE_UNKNOWN,
	TYPE_INT,
	TYPE_FLOAT,
};

enum {
	ATOM_FIELD_UNKNOWN,
	ATOM_FIELD_ID,
	ATOM_FIELD_TYPE,
	ATOM_FIELD_X,
	ATOM_FIELD_Y,
	ATOM_FIELD_Z,
	ATOM_FIELD_RESID,
	ATOM_FIELD_Q,
	ATOM_FIELD_MASS,
	ATOM_FIELD_COUNT
};

static const char* atom_field_name[ATOM_FIELD_COUNT] = {
	"unknown",
	"id",
	"type",
	"x",
	"y",
	"z",
	"resid",
	"q",
	"mass",
};

static const int atom_field_type[ATOM_FIELD_COUNT] = {
	TYPE_UNKNOWN,
	TYPE_INT,
	TYPE_INT,
	TYPE_FLOAT,
	TYPE_FLOAT,
	TYPE_FLOAT,
	TYPE_INT,
	TYPE_FLOAT,
	TYPE_FLOAT,
};

static int atom_field_required[] = {1,2,3,4,5};

static size_t interpret_format(int mappings[ATOM_FIELD_COUNT], const char* atom_format, char* err_buf, size_t err_cap) {
	MEMSET(mappings, -1, sizeof(int) * ATOM_FIELD_COUNT);

	str_t str = str_trim(str_from_cstr(atom_format));
	if (str_empty(str)) {
		if (err_buf) {
			snprintf(err_buf, err_cap, "Lammps atom format: Empty format string");
		} else {
			MD_LOG_ERROR("Lammps atom format: Empty format string");
		}
		return 0;
	}

	str_t tokens[16];
	const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &str);
	for (size_t i = 0; i < num_tokens; ++i) {
		const str_t token = tokens[i];
		for (size_t j = 1; j < ATOM_FIELD_COUNT; ++j) {
			if (str_eq(token, str_from_cstr(atom_field_name[j]))) {
				if (mappings[j] != -1) {
					if (err_buf) {
						snprintf(err_buf, err_cap, "Lammps atom format: Contains duplicate field name '"STR_FMT"'", STR_ARG(token));
					} else {
						MD_LOG_ERROR("Lammps atom format: Contains duplicate field name '"STR_FMT"'", STR_ARG(token));
					}
					return 0;
				}
				mappings[j] = (int)i;
				break;
			}
		}
	}

	// Check required mappings
	for (size_t i = 0; i < ARRAY_SIZE(atom_field_required); ++i) {
		int idx = atom_field_required[i];
		if (mappings[idx] == -1) {
			if (err_buf) {
				snprintf(err_buf, err_cap, "Lammps atom format: Missing required identifier '%s'", atom_field_name[idx]);
			} else {
				MD_LOG_ERROR("Lammps atom format: Missing required identifier '%s'", atom_field_name[idx]);
			}
			return 0;
		}
	}

	return num_tokens;
}

// Only detects Atomic or Full which are predefined strings that we can return
static md_lammps_atom_format_t detect_atom_format(md_buffered_reader_t* reader) {
	str_t tok[16];
	str_t atom_lines[8];
	str_t line;
	str_t hint = {0};
	size_t atom_line_count = 0;
	size_t line_count = 0;
	md_lammps_atom_format_t format = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;

	while (md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok > 0) {
			if (str_eq(tok[0], STR_LIT("Atoms"))) {
				if (num_tok == 3 && str_eq(tok[1], STR_LIT("#"))) {
					hint = tok[2];
				}

				// Read empty line
				md_buffered_reader_extract_line(&line, reader);
				if (!str_empty(line)) {
					// Expected empty line here
					return format;
				}

				// We found Atoms entry, now read some lines
				for (size_t i = 0; i < ARRAY_SIZE(atom_lines); ++i) {
					md_buffered_reader_extract_line(&line, reader);
					if (str_empty(line)) {
						break;
					}
					atom_lines[atom_line_count++] = line;
				}
				break;
			}
		}

		if (line_count++ > 256) {
			MD_LOG_ERROR("Failed to detect atom format, could not find 'Atoms' entry");
			return format;
		}
	}

	if (atom_line_count > 0) {
		if (!str_empty(hint)) {
			for (size_t i = 0; i < MD_LAMMPS_ATOM_FORMAT_COUNT; ++i) {
				if (str_eq_cstr(hint, atom_format_name[i])) {
					format = i;
					break;
				}
			}
			if (!format) {
				MD_LOG_ERROR("Could not detect atom format, using hint '"STR_FMT"'", STR_ARG(hint));
			}
		} else {
			size_t num_tokens = extract_tokens(tok, ARRAY_SIZE(tok), &atom_lines[0]);
			if (num_tokens == 5) {
				format = MD_LAMMPS_ATOM_FORMAT_ATOMIC;
			} else if (num_tokens == 8 && is_int(tok[5]) && is_int(tok[6]) && is_int(tok[7]))  {
				format = MD_LAMMPS_ATOM_FORMAT_ATOMIC;
			}
		}

		// Verify format against the lines we read
		if (format) {
			int mappings[ATOM_FIELD_COUNT];
			size_t num_fields = interpret_format(mappings, atom_format_string[format], 0, 0);
			if (num_fields) {
				for (size_t i = 0; i < atom_line_count; ++i) {
					str_t atom_line = atom_lines[i];
					const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &atom_line);
					if (num_tok < num_fields) {
						MD_LOG_ERROR("Failed to detect atom format, expected %zu tokens, got %zu", num_fields, num_tok);
						format = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
						break;
					}

					for (size_t j = 0; j < ATOM_FIELD_COUNT; ++j) {
						if (mappings[j] != -1) {
							const int idx = mappings[j];
							str_t str = tok[idx];

							switch (atom_field_type[j]) {
							case TYPE_INT: {
								if (!is_int(str)) {
									MD_LOG_ERROR("Failed to detect atom format, expected int at line %zu, token %zu", i, j);
									format = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
									break;
								}
							} break;
							case TYPE_FLOAT: {
								if (!is_float(str)) {
									MD_LOG_ERROR("Failed to detect atom format, expected float at line %zu, token %zu", i, j);
									format = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
									break;
								}
							} break;
							default:
								format = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
								break;
							}
						} 
					}
				}
				return format;
			}
		} 
	}

	return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
}

static bool parse_atoms(md_lammps_atom_t out_atoms[], size_t num_atoms, md_buffered_reader_t* reader, const int mappings[]) {
	ASSERT(mappings[ATOM_FIELD_ID]	 != -1);
	ASSERT(mappings[ATOM_FIELD_TYPE] != -1);
	ASSERT(mappings[ATOM_FIELD_X]	 != -1);
	ASSERT(mappings[ATOM_FIELD_Y]	 != -1);
	ASSERT(mappings[ATOM_FIELD_Z]	 != -1);

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
		atom->id = (int32_t)parse_int(tok[mappings[ATOM_FIELD_ID]]);
		atom->resid = mappings[ATOM_FIELD_RESID] != -1 ? (int32_t)parse_int(tok[mappings[ATOM_FIELD_RESID]]) : -1;
		atom->type = (int32_t)parse_int(tok[mappings[ATOM_FIELD_TYPE]]);
		atom->charge = mappings[ATOM_FIELD_Q] != -1 ? (float)parse_float(tok[mappings[ATOM_FIELD_Q]]) : 0.0f;
		atom->mass = mappings[ATOM_FIELD_MASS] != -1 ? (float)parse_float(tok[mappings[ATOM_FIELD_MASS]]) : 0.0f;
		atom->x = (float)parse_float(tok[mappings[ATOM_FIELD_X]]);
		atom->y = (float)parse_float(tok[mappings[ATOM_FIELD_Y]]);
		atom->z = (float)parse_float(tok[mappings[ATOM_FIELD_Z]]);
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

size_t md_lammps_atom_format_count() {
	return MD_LAMMPS_ATOM_FORMAT_COUNT;
}

const char** md_lammps_atom_format_names() {
	return atom_format_name;
}

const char** md_lammps_atom_format_strings() {
	return atom_format_string;
}

md_lammps_atom_format_t md_lammps_atom_format_from_str(str_t str) {
	md_buffered_reader_t line_reader = md_buffered_reader_from_str(str);
	return detect_atom_format(&line_reader);
}

md_lammps_atom_format_t md_lammps_atom_format_from_file(str_t str) {
	md_file_o* file = md_file_open(str, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Could not open file '%.*s'", str.len, str.ptr);
		return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
	}

	const size_t cap = KILOBYTES(4);
	char* buf = md_temp_push(cap);

	md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
	md_lammps_atom_format_t format = detect_atom_format(&line_reader);

	md_temp_pop(cap);
	md_file_close(file);

	return format;
}

static bool md_lammps_data_parse(md_lammps_data_t* data, md_buffered_reader_t* reader, const char* format, struct md_allocator_i* alloc) {
	ASSERT(data);
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tok[16];

	int mappings[ARRAY_SIZE(atom_field_name)] = {0};
	size_t num_fields = interpret_format(mappings, format, 0, 0);
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
		if (num_tok > 0 && str_eq(tok[0], STR_LIT("Atoms"))) {
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
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Bonds"))) {
			if (!data->num_bonds) {
				MD_LOG_ERROR("Encountered Bond entries, but number of bonds were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->bonds, data->num_bonds, alloc);
			if (!parse_bonds(data->bonds, data->num_bonds, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Angles"))) {
			if (!data->num_angles) {
				MD_LOG_ERROR("Encountered Angle entries, but number of angles were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Dihedrals"))) {
			if (!data->num_dihedrals) {
				MD_LOG_ERROR("Encountered Dihedral entries, but number of dihedrals were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Masses"))) {
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
			if (str_eq(tok[1], STR_LIT("atoms"))) {
				data->num_atoms = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("bonds"))) {
				data->num_bonds = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("angles"))) {
				data->num_angles = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("dihedrals"))) {
				data->num_dihedrals = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("impropers"))) {
				data->num_impropers = (int32_t)parse_int(tok[0]);
			}
		} else if (num_tok == 3 && str_eq(tok[2], STR_LIT("types")) && is_int(tok[0])) {
			if (str_eq(tok[1], STR_LIT("atom"))) {
				data->num_atom_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("bond"))) {
				data->num_bond_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("angle"))) {
				data->num_angle_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("dihedral"))) {
				data->num_dihedral_types = (int32_t)parse_int(tok[0]);
			} else if (str_eq(tok[1], STR_LIT("improper"))) {
				data->num_improper_types = (int32_t)parse_int(tok[0]);
			} 
		} else if (num_tok == 4 && str_eq(tok[2], STR_LIT("xlo")) && str_eq(tok[3], STR_LIT("xhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.xlo = (float)parse_float(tok[0]);
			data->cell.xhi = (float)parse_float(tok[1]);
		} else if (num_tok == 4 && str_eq(tok[2], STR_LIT("ylo")) && str_eq(tok[3], STR_LIT("yhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.ylo = (float)parse_float(tok[0]);
			data->cell.yhi = (float)parse_float(tok[1]);
		} else if (num_tok == 4 && str_eq(tok[2], STR_LIT("zlo")) && str_eq(tok[3], STR_LIT("zhi"))) {
			if (!is_float(tok[0]) || !is_float(tok[1])) {
				MD_LOG_ERROR("Failed to parse cell definition");
				return false;
			}
			data->cell.zlo = (float)parse_float(tok[0]);
			data->cell.zhi = (float)parse_float(tok[1]);
		} else if (num_tok == 6 && str_eq(tok[3], STR_LIT("xy")) && str_eq(tok[4], STR_LIT("xz")) && str_eq(tok[5], STR_LIT("yz"))) {
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

bool md_lammps_validate_atom_format(const char* format, char* err_buf, size_t err_cap) {
	int mappings[ARRAY_SIZE(atom_field_name)] = {0};
	return interpret_format(mappings, format, err_buf, err_cap) > 0;
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

bool md_lammps_molecule_init(md_molecule_t* mol, const md_lammps_data_t* data, md_allocator_i* alloc) {
	ASSERT(mol);
	ASSERT(data);
	ASSERT(alloc);

	MEMSET(mol, 0, sizeof(md_molecule_t));
	const size_t num_atoms = data->num_atoms;

	md_array_resize(mol->atom.x,	 num_atoms, alloc);
	md_array_resize(mol->atom.y,	 num_atoms, alloc);
	md_array_resize(mol->atom.z,	 num_atoms, alloc);
	md_array_resize(mol->atom.mass,  num_atoms, alloc);

	bool has_resid = false;
	if (num_atoms > 0 && data->atoms[0].resid != -1) {
		md_array_resize(mol->atom.resid, num_atoms, alloc);
		has_resid = true;
	}

	for (size_t i = 0; i < num_atoms; ++i) {
		mol->atom.x[i] = data->atoms[i].x;
		mol->atom.y[i] = data->atoms[i].y;
		mol->atom.z[i] = data->atoms[i].z;
		mol->atom.mass[i] = data->atoms[i].mass;
		if (has_resid) {
			mol->atom.resid[i] = data->atoms[i].resid;
		}
	}

	mol->atom.count = num_atoms;

	//Set elements
	md_array_resize(mol->atom.element, num_atoms, alloc);
	if (!md_util_element_from_mass(mol->atom.element, mol->atom.mass, num_atoms)) {
		MD_LOG_ERROR("One or more masses are missing matching element");
	}

	//Create unit cell
	float M[3][3] = {0};
	M[0][0] = data->cell.xhi - data->cell.xlo;
	M[1][1] = data->cell.yhi - data->cell.ylo;
	M[2][2] = data->cell.zhi - data->cell.zlo;
	M[1][0] = data->cell.xy;
	M[2][0] = data->cell.xz;
	M[2][1] = data->cell.yz;
	mol->unit_cell = md_util_unit_cell_from_matrix(M);
	
	return true;
}

static bool lammps_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
	if (!arg) {
		MD_LOG_ERROR("Missing required argument for lammps molecule loader");
		return false;
	}

	const md_lammps_molecule_loader_arg_t* lammps_arg = (const md_lammps_molecule_loader_arg_t*)arg;
	if (lammps_arg->type != MD_LAMMPS_MOLECULE_LOADER_ARG_TYPE) {
		MD_LOG_ERROR("Invalid argument type for lammps molecule loader");
		return false;
	}

	const char* format = lammps_arg->atom_format_str;

	md_lammps_data_t data = { 0 };
	bool success = false;
	if (md_lammps_data_parse_str(&data, str, format, md_heap_allocator)) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_heap_allocator);

	return success;
}

static bool lammps_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
	if (!arg) {
		MD_LOG_ERROR("Missing required argument for lammps molecule loader");
		return false;
	}

	const md_lammps_molecule_loader_arg_t* lammps_arg = (const md_lammps_molecule_loader_arg_t*)arg;
	if (lammps_arg->type != MD_LAMMPS_MOLECULE_LOADER_ARG_TYPE) {
		MD_LOG_ERROR("Invalid argument type for lammps molecule loader");
		return false;
	}

	const char* format = lammps_arg->atom_format_str;

	md_lammps_data_t data = { 0 };
	bool success = false;
	if (md_lammps_data_parse_file(&data, filename, format, md_heap_allocator)) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_heap_allocator);

	return success;
}

md_lammps_molecule_loader_arg_t md_lammps_molecule_loader_arg(const char* atom_format_str) {
	md_lammps_molecule_loader_arg_t arg = { 0 };
	arg.type = MD_LAMMPS_MOLECULE_LOADER_ARG_TYPE;
	arg.atom_format_str = atom_format_str;
	return arg;
}

//Cant use interface as the format parameter is needed as well
static md_molecule_loader_i lammps_api = {
	lammps_init_from_str,
	lammps_init_from_file,
};

md_molecule_loader_i* md_lammps_molecule_api() {
	return &lammps_api;
}
