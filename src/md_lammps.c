#include <md_lammps.h>

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_parse.h>

#include <string.h>

#define MD_LAMMPS_TRAJ_MAGIC 0x2312ad7b78a9bc20
#define MD_LAMMPS_CACHE_MAGIC 0x89172bab
#define MD_LAMMPS_CACHE_VERSION 15
#define MD_LAMMPS_MOLECULE_LOADER_ARG_TYPE 0x341293abc8273650

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

// There are more fields such as velocities and forces, but for now we only care about coordinate data
// @TODO: Support more fields and auxiliary data
enum {
	COORD_FIELD_X,
	COORD_FIELD_Y,
	COORD_FIELD_Z,
	COORD_FIELD_IX,
	COORD_FIELD_IY,
	COORD_FIELD_IZ,
	COORD_FIELD_COUNT,
};

enum {
	COORD_FLAG_NONE      = 0,
	COORD_FLAG_CARTESIAN = 1,
	COORD_FLAG_SCALED    = 2,
	COORD_FLAG_UNWRAP    = 4,
};

typedef struct coord_mappings_t {
	int8_t  id_idx;
	int8_t  coord_idx[3];
	int8_t  image_idx[3];
	int8_t  num_coord_tokens;	// This encodes the expected number of tokens in a coordinate entry
	int8_t  flags;
	int8_t  _pad[7];
} coord_mappings_t;

typedef struct lammps_trajectory_t {
	uint64_t magic;
	int64_t* frame_offsets;

	md_file_o* file;
	md_trajectory_header_t header;
	coord_mappings_t coord_mappings;

	md_allocator_i* allocator;
	md_mutex_t mutex;
} lammps_trajectory_t;

typedef struct lammps_cache_t {
	md_trajectory_cache_header_t header;
	int64_t* frame_offsets;
	int64_t* frame_times;
	coord_mappings_t coord_mappings;
} lammps_cache_t;

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

static const int atom_field_required[] = {1,2,3,4,5};

#if 0
static const char* coord_field_name[] = {
	"unknown",
	"x",
	"y",
	"z",
	"xs",
	"ys",
	"zs",
	"xu",
	"yu",
	"zu",
	"xsu",
	"ysu",
	"zsu",
	"ix",
	"iy",
	"iz",
};
#endif

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
	md_lammps_atom_format_t format = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;

	while (md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok > 0) {
			if (str_eq(tok[0], STR_LIT("Atoms"))) {
				if (num_tok >= 3 && str_eq(tok[1], STR_LIT("#"))) {
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

		if (md_buffered_reader_tellg(reader) > MEGABYTES(1)) {
			// If we cannot find the Atom entry within the first megabyte, we give up
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
			if (num_tokens == 5 || (num_tokens == 8 && is_int(tok[5]) && is_int(tok[6]) && is_int(tok[7]))) {
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

		md_lammps_atom_t* atom = out_atoms + read_atoms;
		atom->id = (int32_t)parse_int(tok[mappings[ATOM_FIELD_ID]]);
		atom->resid = mappings[ATOM_FIELD_RESID] != -1 ? (int32_t)parse_int(tok[mappings[ATOM_FIELD_RESID]]) : -1;
		atom->type = (int32_t)parse_int(tok[mappings[ATOM_FIELD_TYPE]]);
		atom->charge = mappings[ATOM_FIELD_Q] != -1 ? (float)parse_float(tok[mappings[ATOM_FIELD_Q]]) : 0.0f;
		atom->mass = mappings[ATOM_FIELD_MASS] != -1 ? (float)parse_float(tok[mappings[ATOM_FIELD_MASS]]) : 0.0f;
		atom->x = (float)parse_float(tok[mappings[ATOM_FIELD_X]]);
		atom->y = (float)parse_float(tok[mappings[ATOM_FIELD_Y]]);
		atom->z = (float)parse_float(tok[mappings[ATOM_FIELD_Z]]);

		read_atoms += 1;
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

static size_t parse_masses(float* mass_type_table, size_t mass_type_capacity, size_t expected_count, md_buffered_reader_t* reader) {
	str_t tok[4];
	str_t line;
	size_t extracted_count = 0;
	for (size_t i = 0; i < expected_count; ++i) {
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to extract mass line");
			return 0;
		}
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 2) {
			MD_LOG_ERROR("Failed to parse mass line, expected 2 tokens, got %zu", num_tok);
			return 0;
		}
		int   type = (int)parse_int(tok[0]);
		float mass = (float)parse_float(tok[1]);
		if (type < 0 || (size_t)type >= mass_type_capacity) {
			MD_LOG_ERROR("Invalid atom type index in Masses: %d (capacity %zu)", type, mass_type_capacity);
			return 0;
		}
		mass_type_table[type] = mass;
		extracted_count += 1;
	}

	return extracted_count;
}

size_t md_lammps_atom_format_count(void) {
	return MD_LAMMPS_ATOM_FORMAT_COUNT;
}

const char** md_lammps_atom_format_names(void) {
	return atom_format_name;
}

const char** md_lammps_atom_format_strings(void) {
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

	const size_t cap = MEGABYTES(1);
	char* buf = md_temp_push(cap);

	md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
	md_lammps_atom_format_t format = detect_atom_format(&line_reader);

	md_temp_pop(cap);
	md_file_close(file);

	return format;
}

static int compare_atom(const void* a, const void* b) {
	const md_lammps_atom_t* atom_a = (const md_lammps_atom_t*)a;
	const md_lammps_atom_t* atom_b = (const md_lammps_atom_t*)b;
	return atom_a->id - atom_b->id;
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

	double atom_type_mass_table[512] = {0};

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
			// Sort atoms by id
			qsort(data->atoms, data->num_atoms, sizeof(md_lammps_atom_t), compare_atom);

			for (size_t i = 0; i < data->num_atoms; ++i) {
				int32_t type = data->atoms[i].type;
				data->atoms[i].mass = type < (int)ARRAY_SIZE(atom_type_mass_table) ? atom_type_mass_table[type] : 0.0f;
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
			md_array_resize(data->dihedrals, data->num_dihedrals, alloc);
			if (!parse_dihedrals(data->dihedrals, data->num_dihedrals, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Masses"))) {
			if (!data->num_atom_types) {
				MD_LOG_ERROR("Encountered Mass entries, but number of atom types were not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			if (parse_masses(atom_type_mass_table, ARRAY_SIZE(atom_type_mass_table), data->num_atom_types, reader) != data->num_atom_types) {
				MD_LOG_ERROR("Number of masses in table did not match the number of atom types");
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
		char* buf = md_alloc(md_get_heap_allocator(), cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = md_lammps_data_parse(data, &line_reader, format, alloc);

		md_free(md_get_heap_allocator(), buf, cap);
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
	const size_t capacity = ROUND_UP(data->num_atoms, 16);

	md_array_resize(mol->atom.type,  capacity, alloc);
	md_array_resize(mol->atom.x,	 capacity, alloc);
	md_array_resize(mol->atom.y,	 capacity, alloc);
	md_array_resize(mol->atom.z,	 capacity, alloc);
	md_array_resize(mol->atom.mass,  capacity, alloc);

	bool has_resid = false;
	if (data->num_atoms > 0 && data->atoms[0].resid != -1) {
		md_array_resize(mol->atom.resid, capacity, alloc);
		has_resid = true;
	}

	for (size_t i = 0; i < data->num_atoms; ++i) {
		mol->atom.type[i].len = (uint8_t)snprintf(mol->atom.type[i].buf, sizeof(mol->atom.type[i].buf), "%i", data->atoms[i].type);
		mol->atom.x[i] = data->atoms[i].x - data->cell.xlo;
		mol->atom.y[i] = data->atoms[i].y - data->cell.ylo;
		mol->atom.z[i] = data->atoms[i].z - data->cell.zlo;
		mol->atom.mass[i] = data->atoms[i].mass;
		if (has_resid) {
			mol->atom.resid[i] = data->atoms[i].resid;
		}
	}

	//Set elements using conservative mass→element mapping
	md_array_resize(mol->atom.element, capacity, alloc);
	if (!md_util_lammps_element_from_mass(mol->atom.element, mol->atom.mass, data->num_atoms)) {
		// CG/reduced-units detected or mapping failed, leave elements as 0
		MD_LOG_DEBUG("LAMMPS data appears to be coarse-grained or reduced-units, elements left unassigned");
	}

	mol->atom.count = data->num_atoms;

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
	if (md_lammps_data_parse_str(&data, str, format, md_get_heap_allocator())) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_get_heap_allocator());

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
	if (md_lammps_data_parse_file(&data, filename, format, md_get_heap_allocator())) {
		success = md_lammps_molecule_init(mol, &data, alloc);
	}
	md_lammps_data_free(&data, md_get_heap_allocator());

	return success;
}

md_lammps_molecule_loader_arg_t md_lammps_molecule_loader_arg(const char* atom_format_str) {
	md_lammps_molecule_loader_arg_t arg = { 0 };
	arg.type = MD_LAMMPS_MOLECULE_LOADER_ARG_TYPE;
	arg.atom_format_str = atom_format_str;
	return arg;
}

static md_molecule_loader_i lammps_api = {
	lammps_init_from_str,
	lammps_init_from_file,
};

md_molecule_loader_i* md_lammps_molecule_api(void) {
	return &lammps_api;
}

// TRAJECTORY OPERATIONS

//Reads data that is useful later when we want to parse a frame from the trajectory
bool lammps_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
	lammps_trajectory_t* traj = (lammps_trajectory_t*)inst;
	ASSERT(traj);
	ASSERT(traj->magic == MD_LAMMPS_TRAJ_MAGIC);
	ASSERT(header);

	*header = traj->header;
	return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
// Returns size in bytes of frame, frame_data_ptr is optional and is the destination to write the frame data to.
size_t lammps_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
	lammps_trajectory_t* traj_data = (lammps_trajectory_t*)inst;
	ASSERT(traj_data);
	ASSERT(traj_data->magic == MD_LAMMPS_TRAJ_MAGIC);

	if (!traj_data->file) {
		MD_LOG_ERROR("File handle is NULL");
		return 0;
	}

	if (!traj_data->frame_offsets) {
		MD_LOG_ERROR("Frame offsets is empty");
		return 0;
	}
	
	if (frame_idx < 0 || (int64_t)traj_data->header.num_frames <= frame_idx) {
		MD_LOG_ERROR("Frame index is out of range");
		return 0;
	}

	const int64_t beg = traj_data->frame_offsets[frame_idx + 0];
	const int64_t end = traj_data->frame_offsets[frame_idx + 1];
	const size_t frame_size = (size_t)MAX(0, end - beg);
	const size_t total_size = sizeof(int64_t) + frame_size;

	if (frame_data_ptr) {
		// Store the index to the frame since this is generally not found within the actual frame data
		int64_t* ptr = (int64_t*)frame_data_ptr;
		ptr[0] = frame_idx;

		ASSERT(traj_data->file);
		md_mutex_lock(&traj_data->mutex);
		md_file_seek(traj_data->file, beg, MD_FILE_BEG);
		const size_t bytes_read = md_file_read(traj_data->file, &ptr[1], frame_size);
		(void)bytes_read;
		md_mutex_unlock(&traj_data->mutex);
		ASSERT(frame_size == bytes_read);
	}

	return total_size;
}

typedef struct {
	double xlo, xhi, xy;
	double ylo, yhi, xz;
	double zlo, zhi, yz;
} box_bounds_t;

typedef struct {
	int64_t timestep;
	size_t num_atoms;
	box_bounds_t box_bounds;
} header_t;

static bool parse_box_bounds(box_bounds_t* box_bounds, md_buffered_reader_t* reader) {
	ASSERT(reader);
	str_t line;

	if (!md_buffered_reader_extract_line(&line, reader)) {
		MD_LOG_ERROR("Could not extract box bounds");
		return false;
	}

	if (!str_eq_cstr_n(line, "ITEM: BOX BOUNDS", 16)) {
		MD_LOG_ERROR("Unexpected beginning of line: '" STR_FMT "' expected ITEM: BOX BOUNDS", STR_ARG(line));
		return false;
	}

	line = str_trim_beg(str_substr(line, 16, SIZE_MAX));

	// Should either match "pp pp pp" or "xy xz yz pp pp pp"
	if (!str_eq_cstr_n(line, "pp pp pp", 8) && !str_eq_cstr_n(line, "xy xz yz pp pp pp", 17)) {
		MD_LOG_ERROR("Unrecognized format in ITEM: BOX BOUNDS: '" STR_FMT "'", STR_ARG(line));
		return false;
	}

	str_t tokens[4];
	double values[3][3] = {0};
	for (size_t i = 0; i < 3; ++i) {
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to extract box bounds");
			return false;
		}
		size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
		if (num_tokens < 2) {
			MD_LOG_ERROR("Failed to extract box bounds");
			return false;
		}
		for (size_t j = 0; j < num_tokens; ++j) {
			values[i][j] = parse_float(tokens[j]);
		}
	}

	if (box_bounds) {
		MEMCPY(box_bounds, values, sizeof(box_bounds_t));
	}

	return true;
}

static bool parse_header(header_t* header, md_buffered_reader_t* reader) {
	ASSERT(reader);
	str_t line;

	if (!md_buffered_reader_peek_line(&line, reader) ||
		!str_eq_cstr_n(line, "ITEM: TIMESTEP", 14))
	{
		return false;
	}

	md_buffered_reader_skip_line(reader);
	if (md_buffered_reader_extract_line(&line, reader)) {
		if (header) {
			header->timestep = parse_int(line);
		} else if (!is_int(line)) {
			MD_LOG_ERROR("Could not extract timestep");
			return false;
		}
	} else {
		MD_LOG_ERROR("Could not extract timestep");
		return false;
	}

	if (md_buffered_reader_peek_line(&line, reader) &&
		str_eq_cstr_n(line, "ITEM: NUMBER OF ATOMS", 21) &&
		md_buffered_reader_skip_line(reader) &&
		md_buffered_reader_extract_line(&line, reader))
	{
		if (header) {
			header->num_atoms = parse_int(line);
		} else if (!is_int(line)) {
			MD_LOG_ERROR("Could not extract number of atoms");
			return false;
		}
	} else {
		MD_LOG_ERROR("Could not extract number of atoms");
		return false;
	}

	if (md_buffered_reader_peek_line(&line, reader) &&
		str_eq_cstr_n(line, "ITEM: BOX BOUNDS", 16))
	{
		if (!parse_box_bounds(header ? &header->box_bounds : NULL, reader)) {
			return false;
		}
	} else {
		MD_LOG_ERROR("Could not extract box bounds");
		return false;
	}

	return true;
}

// num_atom_tokens returns the number of expected tokens for the atom coordinates
static bool parse_coord_mappings(coord_mappings_t* mappings, str_t str) {
	ASSERT(mappings);

	str_t tokens[32];
	size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &str);
	ASSERT(str_eq(tokens[0], STR_LIT("ITEM:")));
	ASSERT(str_eq(tokens[1], STR_LIT("ATOMS")));

	if (num_tokens < 2) {
		return false;
	}

	const str_t* labels = tokens + 2;
	size_t num_labels = num_tokens - 2;

	MEMSET(mappings, -1, sizeof(coord_mappings_t));

	for (size_t i = 0; i < num_labels; ++i) {
		// @NOTE: We only try to match the first occurence of the coordinate mapping
		// This is because there is an implicit order to the mappings that we prefer to use
		if (str_eq(labels[i], STR_LIT("id"))) {
			mappings->id_idx = (int8_t)i;
		}
		if (mappings->flags == -1) {
			if (str_eq(labels[i], STR_LIT("x")) && str_eq(labels[i+1], STR_LIT("y")) && str_eq(labels[i+2], STR_LIT("z"))) {
				mappings->coord_idx[0] = (int8_t)i;
				mappings->coord_idx[1] = (int8_t)i + 1;
				mappings->coord_idx[2] = (int8_t)i + 2;
				mappings->flags = COORD_FLAG_CARTESIAN;
				i += 2;
			}
			else if (str_eq(labels[i], STR_LIT("xs")) && str_eq(labels[i+1], STR_LIT("ys")) && str_eq(labels[i+2], STR_LIT("zs"))) {
				mappings->coord_idx[0] = (int8_t)i;
				mappings->coord_idx[1] = (int8_t)i + 1;
				mappings->coord_idx[2] = (int8_t)i + 2;
				mappings->flags = COORD_FLAG_SCALED;
				i += 2;
			}
			else if (str_eq(labels[i], STR_LIT("xu")) && str_eq(labels[i+1], STR_LIT("yu")) && str_eq(labels[i+2], STR_LIT("zu"))) {
				mappings->coord_idx[0] = (int8_t)i;
				mappings->coord_idx[1] = (int8_t)i + 1;
				mappings->coord_idx[2] = (int8_t)i + 2;
				mappings->flags = COORD_FLAG_CARTESIAN;
				i += 2;
			}
			else if (str_eq(labels[i], STR_LIT("xsu")) && str_eq(labels[i+1], STR_LIT("ysu")) && str_eq(labels[i+2], STR_LIT("zsu"))) {
				mappings->coord_idx[0] = (int8_t)i;
				mappings->coord_idx[1] = (int8_t)i + 1;
				mappings->coord_idx[2] = (int8_t)i + 2;
				mappings->flags = COORD_FLAG_SCALED;
				i += 2;
			}
		} else {
			if (str_eq(labels[i], STR_LIT("ix")) && str_eq(labels[i+1], STR_LIT("iy")) && str_eq(labels[i+2], STR_LIT("iz"))) {
				mappings->image_idx[0] = (int8_t)i;
				mappings->image_idx[1] = (int8_t)i + 1;
				mappings->image_idx[2] = (int8_t)i + 2;
				mappings->flags = COORD_FLAG_UNWRAP;
				i += 2;
			}
		}
	}

	ASSERT(num_labels < 127);
	mappings->num_coord_tokens = (int8_t)num_labels;

	if (mappings->id_idx != -1 && mappings->flags != -1) {
		if (mappings->flags & COORD_FLAG_UNWRAP) {
			return mappings->image_idx[0] != -1;
		}
		return true;
	}

	return false;
}

typedef struct {
	int32_t id;
	float x, y, z;
} id_xyz_t;

int compare_id_xyz(const void* a, const void* b) {
	const id_xyz_t* id_xyz_a = (const id_xyz_t*)a;
	const id_xyz_t* id_xyz_b = (const id_xyz_t*)b;
	return id_xyz_a->id - id_xyz_b->id;
}

bool lammps_decode_frame_data(struct md_trajectory_o* inst, const void* data_ptr, size_t data_size, md_trajectory_frame_header_t* out_frame_header, float* out_x, float* out_y, float* out_z) {
	ASSERT(inst);
	ASSERT(data_ptr);
	ASSERT(data_size);

	str_t tokens[32];
	int64_t timestep = 0;
	int64_t frame_idx = ((int64_t*)data_ptr)[0];
	md_unit_cell_t unit_cell = {0};

	bool output_header = out_frame_header != NULL;
	bool output_coords = out_x != NULL && out_y != NULL && out_z != NULL;

	lammps_trajectory_t* traj_data = (lammps_trajectory_t*)inst;
	if (traj_data->magic != MD_LAMMPS_TRAJ_MAGIC) {
		MD_LOG_ERROR("Error when decoding frame header, lammps magic did not match");
		return false;
	}
	if (frame_idx < 0 || frame_idx >= (int64_t)traj_data->header.num_frames) {
		MD_LOG_ERROR("Error when decoding frame data, corrupt frame index");
		return false;
	}

	str_t str = { .ptr = (const char*)(data_ptr) + sizeof(int64_t), .len = data_size - sizeof(int64_t) };
	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	str_t line;

	header_t header;
	if (!parse_header(&header, &reader)) {
		MD_LOG_ERROR("Could not parse header");
		return false;
	}

	// https://docs.lammps.org/Howto_triclinic.html
	double xlo = header.box_bounds.xlo - MIN(0.0, MIN(header.box_bounds.xy, MIN(header.box_bounds.xz, header.box_bounds.yz)));
	double xhi = header.box_bounds.xhi - MAX(0.0, MAX(header.box_bounds.xy, MAX(header.box_bounds.xz, header.box_bounds.yz)));
	double xz  = header.box_bounds.xz;
	double ylo = header.box_bounds.ylo - MIN(0.0, header.box_bounds.yz);
	double yhi = header.box_bounds.yhi - MAX(0.0, header.box_bounds.yz);
	double xy  = header.box_bounds.xy;
	double zlo = header.box_bounds.zlo;
	double zhi = header.box_bounds.zhi;
	double yz  = header.box_bounds.yz;

	double xlen = xhi - xlo;
	double ylen = yhi - ylo;
	double zlen = zhi - zlo;

	unit_cell = md_util_unit_cell_from_triclinic(xlen, ylen, zlen, xy, xz, yz);

	// transform matrix to apply
	mat4_t M = mat4_translate(-(float)xlo, -(float)ylo, -(float)zlo);
	if (traj_data->coord_mappings.flags & COORD_FLAG_SCALED) {
		// Scaling
		M = mat4_from_mat3(unit_cell.basis);
	}

	if (output_coords) {
		if (!md_buffered_reader_extract_line(&line, &reader) || !str_eq_cstr_n(line, "ITEM: ATOMS", 11)) {
			MD_LOG_ERROR("Expected ITEM: ATOMS after header");
			return false;
		}
		size_t line_count = 0;
		size_t expected_num_tokens = traj_data->coord_mappings.num_coord_tokens;

		// We need to store the coordinates in a temporary buffer since we need to sort them by id
		id_xyz_t* id_xyz = md_alloc(md_get_heap_allocator(), sizeof(id_xyz_t) * header.num_atoms);

		while (md_buffered_reader_extract_line(&line, &reader) && line_count < header.num_atoms) {
			size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
			if (num_tokens != expected_num_tokens) {
				MD_LOG_ERROR("Unexpected number of tokens in ITEM: ATOMS line, got %zu, expected %zu", num_tokens, expected_num_tokens);
				return false;
			}

			int32_t id = (int32_t)parse_int(tokens[traj_data->coord_mappings.id_idx]);
			vec4_t coord = {
				(float)parse_float(tokens[traj_data->coord_mappings.coord_idx[0]]),
				(float)parse_float(tokens[traj_data->coord_mappings.coord_idx[1]]),
				(float)parse_float(tokens[traj_data->coord_mappings.coord_idx[2]]),
				1.0f
			};
			coord = mat4_mul_vec4(M, coord);

			if (traj_data->coord_mappings.flags & COORD_FLAG_UNWRAP) {
				int64_t ix = parse_int(tokens[traj_data->coord_mappings.image_idx[0]]);
				int64_t iy = parse_int(tokens[traj_data->coord_mappings.image_idx[1]]);
				int64_t iz = parse_int(tokens[traj_data->coord_mappings.image_idx[2]]);
				vec4_t trans = {
					(float)(ix * xlen),
					(float)(iy * ylen),
					(float)(iz * zlen),
					1.0f
				};
				coord = vec4_add(coord, trans);
			}

			id_xyz[line_count] = (id_xyz_t){ id, coord.x, coord.y, coord.z };
			line_count += 1;
		}

		// Sort atoms by id
		qsort(id_xyz, line_count, sizeof(id_xyz_t), compare_id_xyz);

		for (size_t i = 0; i < line_count; ++i) {
			out_x[i] = id_xyz[i].x;
			out_y[i] = id_xyz[i].y;
			out_z[i] = id_xyz[i].z;
		}
	}

	if (output_header) {
		out_frame_header->num_atoms = header.num_atoms;
		out_frame_header->index = frame_idx;
		out_frame_header->timestamp = (double)timestep;
		out_frame_header->unit_cell = unit_cell;
	}

	return true;
}

// Parse and validate the trajectory data and record offsets into the file for each frame
static bool lammps_trajectory_parse(lammps_cache_t* cache, md_buffered_reader_t* reader, struct md_allocator_i* alloc) {
	ASSERT(reader);
	ASSERT(alloc);
	str_t line;
	str_t tokens[32];
	size_t num_atoms = 0;
	size_t num_frames = 0;
	coord_mappings_t mappings = {0};

	while (md_buffered_reader_peek_line(&line, reader)) {
		const int64_t offset = md_buffered_reader_tellg(reader);
		header_t header;
		if (!parse_header(&header, reader)) {
			break;
		}

		if (num_atoms == 0) {
			num_atoms = header.num_atoms;
		} else if (num_atoms != header.num_atoms) {
			MD_LOG_ERROR("Number of atoms does not match between frames");
			return false;
		}

		if (!md_buffered_reader_peek_line(&line, reader) ||
			!str_eq_cstr_n(line, "ITEM: ATOMS", 11)) {
			MD_LOG_ERROR("Expected ITEM: ATOMS after header");
			return false;
		}

		if (num_frames == 0) {
			// Parse coord mappings
			if (!parse_coord_mappings(&mappings, line)) {
				MD_LOG_ERROR("Could not parse coord mappings");
				return false;
			}
		}

		// Skip ITEM: ATOMS line
		md_buffered_reader_skip_line(reader);

		// In theory, we could just skip the lines, but we want to validate something in the atom coordinate section
		for (size_t i = 0; i < header.num_atoms; ++i) {
			if (!md_buffered_reader_extract_line(&line, reader)) {
				MD_LOG_ERROR("Could not extract atom line");
				return false;
			}
			size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
			// Expect these to match for all frames
			if (num_tokens != (size_t)mappings.num_coord_tokens) {
				MD_LOG_ERROR("Number of tokens in ITEM: ATOMS line does not match between frames");
				return false;
			}
		}

		md_array_push(cache->frame_offsets, offset, alloc);
		md_array_push(cache->frame_times, header.timestep, alloc);
		num_frames += 1;
	}

	if (num_frames == 0) {
		return false;
	}

	cache->header.num_frames = num_frames;
	cache->header.num_atoms = num_atoms;
	cache->coord_mappings = mappings;
	//We add the end of the file to frame_offsets, so frame_offset size = num_frames + 1

	const int64_t end_of_file = md_buffered_reader_tellg(reader);
	md_array_push(cache->frame_offsets, end_of_file, alloc);

	return true;
}

static bool lammps_trajectory_parse_file(lammps_cache_t* cache, str_t filename, struct md_allocator_i* alloc) {
	bool result = false;
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (file) {
		const int64_t cap = MEGABYTES(1);
		char* buf = md_alloc(md_get_heap_allocator(), cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = lammps_trajectory_parse(cache, &line_reader, alloc);

		md_free(md_get_heap_allocator(), buf, cap);
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

	bool result = true;
	const size_t frame_size = lammps_fetch_frame_data(inst, frame_idx, NULL);
	if (frame_size > 0) { //This check first that there actually is data to be read, before writing to frame_data
		md_allocator_i* alloc = md_get_heap_allocator();

		void* frame_data = md_alloc(alloc, frame_size);
		ASSERT(frame_data);

		const size_t read_size = lammps_fetch_frame_data(inst, frame_idx, frame_data);
		if (read_size != frame_size) {
			MD_LOG_ERROR("Failed to read the expected size");
			goto done;
		}

		//Decode is what actually writes the data to x, y and z which are float pointers
		result = lammps_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
done:
		md_free(alloc, frame_data, frame_size);
	}

	return result;
}

static bool try_read_cache(str_t cache_file, lammps_cache_t* cache, size_t traj_num_bytes, md_allocator_i* alloc) {
	ASSERT(cache);
	ASSERT(alloc);

	md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
	bool result = false;
	if (file) {
		
		if (md_file_read(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
			MD_LOG_ERROR("LAMMPS trajectory cache: failed to read header");
			goto done;
		}

		if (cache->header.magic != MD_LAMMPS_CACHE_MAGIC) {
			MD_LOG_ERROR("LAMMPS trajectory cache: magic was incorrect or corrupt");
			goto done;
		}
		if (cache->header.version != MD_LAMMPS_CACHE_VERSION) {
			MD_LOG_INFO("LAMMPS trajectory cache: version mismatch, expected %i, got %i", MD_LAMMPS_CACHE_VERSION, (int)cache->header.version);
		}
		if (cache->header.num_bytes != traj_num_bytes) {
			MD_LOG_INFO("LAMMPS trajectory cache: trajectory size mismatch, expected %i, got %i", (int)traj_num_bytes, (int)cache->header.num_bytes);
		}
		if (cache->header.num_atoms == 0) {
			MD_LOG_ERROR("LAMMPS trajectory cache: num atoms was zero");
			goto done;
		}
		if (cache->header.num_frames == 0) {
			MD_LOG_ERROR("LAMMPS trajectory cache: num frames was zero");
			goto done;
		}

		size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
		md_array_resize(cache->frame_offsets, cache->header.num_frames + 1, alloc);
		if (md_file_read(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
			MD_LOG_ERROR("LAMMPS trajectory cache: Failed to read offset data");
			md_free(alloc, cache->frame_offsets, offset_bytes);
			goto done;
		}
		size_t num_frame_offsets = md_array_size(cache->frame_offsets);
		if (num_frame_offsets != cache->header.num_frames + 1) {
			MD_LOG_ERROR("LAMMPS trajectory: Read frame offset array size is not correct");
		}

		size_t times_bytes = (cache->header.num_frames) * sizeof(int64_t);
		md_array_resize(cache->frame_times, cache->header.num_frames, alloc);
		if (md_file_read(file, cache->frame_times, times_bytes) != times_bytes) {
			MD_LOG_ERROR("LAMMPS trajectory cache: Failed to read frame times data");
			md_free(alloc, cache->frame_times, times_bytes);
			goto done;
		}
		size_t num_times = md_array_size(cache->frame_times);
		if (num_times != cache->header.num_frames) {
			MD_LOG_ERROR("LAMMPS trajectory: Read frame times array size is not correct");
		}

		if (md_file_read(file, &cache->coord_mappings, sizeof(cache->coord_mappings)) != sizeof(cache->coord_mappings) || cache->coord_mappings.flags == COORD_FLAG_NONE) {
			MD_LOG_ERROR("Failed to read coord type cache, not valid");
			goto done;
		}

		// Test position in file, we expect to be at the end of the file
		if (md_file_tell(file) != (int64_t)md_file_size(file)) {
			MD_LOG_ERROR("PDB trajectory cache: file position was not at the end of the file");
			md_free(alloc, cache->frame_offsets, offset_bytes);
			md_free(alloc, cache->frame_times, times_bytes);
			goto done;
		}

		result = true;
	done:
		md_file_close(file);
	}
	return result;
}

static bool write_cache(lammps_cache_t* cache, str_t cache_file) {
	bool result = false;

	md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("LAMMPS trajectory cache: could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
		return false;
	}

	if (md_file_write(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
		MD_LOG_ERROR("LAMMPS trajectory cache: failed to write header");
		goto done;
	}

	size_t num_frame_offsets = md_array_size(cache->frame_offsets);
	if (num_frame_offsets != cache->header.num_frames + 1) {
		MD_LOG_ERROR("Read frame offset array size is not correct");
	}
	size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
	if (md_file_write(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
		MD_LOG_ERROR("LAMMPS trajectory cache: failed to write frame offsets");
		goto done;
	}

	size_t times_bytes = (cache->header.num_frames) * sizeof(int64_t);
	if (md_file_write(file, cache->frame_times, times_bytes) != times_bytes) {
		MD_LOG_ERROR("LAMMPS trajectory cache: failed to write frame times");
		goto done;
	}

	if (md_file_write(file, &cache->coord_mappings, sizeof(cache->coord_mappings)) != sizeof(cache->coord_mappings)) {
		MD_LOG_ERROR("LAMMPS trajectory cache: failed to write coord type");
		goto done;
	}

	result = true;

done:
	md_file_close(file);
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

md_trajectory_i* md_lammps_trajectory_create(str_t filename, struct md_allocator_i* ext_alloc, uint32_t flags) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file for LAMMPS trajectory");
		return false;
	}

	int64_t filesize = md_file_size(file);
	md_file_close(file);

	md_strb_t sb = md_strb_create(md_get_temp_allocator());
	md_strb_fmt(&sb, STR_FMT ".cache", STR_ARG(filename));
	str_t cache_file = md_strb_to_str(sb);

	lammps_cache_t cache = {0};
	md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

	if (!try_read_cache(cache_file, &cache, filesize, alloc)) {
		//If the cache file does not exist, we create one
		if (!lammps_trajectory_parse_file(&cache, filename, alloc)) {
			MD_LOG_ERROR("LAMMPS trajectory could not be read from file");
			return false;
		}

		cache.header.magic     = MD_LAMMPS_CACHE_MAGIC;
		cache.header.version   = MD_LAMMPS_CACHE_VERSION;
		cache.header.num_bytes = filesize;

		if (!(flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
			// If we fail to write the cache, that's ok, we can inform about it, but do not halt
			if (write_cache(&cache, cache_file)) {
				MD_LOG_INFO("LAMMPS: Successfully created cache file for '" STR_FMT "'", STR_ARG(cache_file));
			}
		}
	}

	size_t max_frame_size = 0;
	//Calculate the max frame size
	for (size_t i = 0; i < cache.header.num_frames; i++) {
		const int64_t beg = cache.frame_offsets[i + 0];
		const int64_t end = cache.frame_offsets[i + 1];
		const size_t frame_size = MAX(0, end - beg);
		max_frame_size = MAX(max_frame_size, frame_size);
	}

	md_array(double) frame_times = md_array_create(double, cache.header.num_frames, alloc);
	for (size_t i = 0; i < cache.header.num_frames; i++) {
		frame_times[i] = (double)cache.frame_times[i];
	}

	void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(lammps_trajectory_t));
	ASSERT(mem);
	MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(lammps_trajectory_t));

	md_trajectory_i* traj = mem;
	lammps_trajectory_t* traj_data = (lammps_trajectory_t*)(traj + 1);

	traj_data->magic = MD_LAMMPS_TRAJ_MAGIC;
	traj_data->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	traj_data->frame_offsets = cache.frame_offsets;
	traj_data->allocator = alloc;
	traj_data->mutex = md_mutex_create();

	traj_data->header = (md_trajectory_header_t){
		.num_frames = cache.header.num_frames,
		.num_atoms = cache.header.num_atoms,
		.max_frame_data_size = max_frame_size,
		.time_unit = md_unit_femtosecond(),
		.frame_times = frame_times,
	};

	traj_data->coord_mappings = cache.coord_mappings;

	traj->inst = (struct md_trajectory_o*)traj_data;
	traj->get_header = lammps_get_header;
	traj->load_frame = lammps_load_frame;
	//traj->fetch_frame_data = lammps_fetch_frame_data;
	//traj->decode_frame_data = lammps_decode_frame_data;

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

md_trajectory_loader_i* md_lammps_trajectory_loader(void) {
	return &lammps_traj_loader;
}
