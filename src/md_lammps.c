#include "md_lammps.h"

#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_parse.h>

#define MD_LAMMPS_TRAJ_MAGIC 0x2312ad7b78a9bc20
#define MD_LAMMPS_TRAJ_READER_MAGIC 0x2312ad7b78a9bc21
#define MD_LAMMPS_CACHE_MAGIC 0x89172bab
#define MD_LAMMPS_CACHE_VERSION 18
#define MD_LAMMPS_SYSTEM_LOADER_ARG_TYPE 0x341293abc8273650


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

typedef struct lammps_cache_t {
	md_run_cache_header_t header;
	int64_t* frame_offsets;
	// The LAMMPS TIMESTEP of each frame
	int64_t* frame_steps;
	coord_mappings_t coord_mappings;
	// The box of each frame, 9 floats, row i box vector i (after the coord mappings on disk)
	float* frame_cells;
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

static size_t parse_masses(md_lammps_atom_type_t types[], size_t num_atom_types, md_buffered_reader_t* reader) {
	str_t tok[4];
	str_t line;
	size_t count = 0;
	for (size_t i = 0; i < num_atom_types; ++i) {
		if (!md_buffered_reader_extract_line(&line, reader)) {
			MD_LOG_ERROR("Failed to extract mass line");
			goto done;
		}
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok < 2) {
			MD_LOG_ERROR("Failed to parse mass line, expected 2 tokens, got %zu", num_tok);
			goto done;
		}
		int   type = (int)parse_int(tok[0]);
		float mass = (float)parse_float(tok[1]);
		if (type < 1) {
			MD_LOG_ERROR("Invalid atom type index in Masses: %d", type);
			goto done;
		}

		if (types[i].id == 0) {
			types[i].id = type;
		} else {
			if (types[i].id != type) {
				MD_LOG_ERROR("Wonky type in Masses: %d, expected %d", type, types[i].id);
				goto done;
			}
		}

		types[i].mass = mass;
		md_atomic_number_t z = md_atomic_number_infer_from_mass(mass);
		if (z != MD_Z_X) {
			types[i].radius = md_atomic_number_vdw_radius(z);
		}
		count += 1;
	}

done:
	return count;
}

// --- LJ radius ---
static double lj_radius(double sigma) {
	return 0.5 * pow(2.0, 1.0/6.0) * sigma; // ≈ 0.561231 * sigma
}

// --- Morse radius ---
static double morse_radius(double r0) {
	return 0.5 * r0;
}

// --- Buckingham radius (numerical solve for rmin) ---
static double buckingham_radius(double A, double rho, double C) {
	// Solve f(r) = -A/rho * exp(-r/rho) + 6C/r^7 = 0
	// Bracket between [rlo, rhi] and use bisection
	double rlo = 1e-6;
	double rhi = 10.0 * rho;
	for (int iter = 0; iter < 200; iter++) {
		double mid = 0.5 * (rlo + rhi);
		double fmid = -(A/rho) * exp(-mid/rho) + 6.0*C/pow(mid,7);
		double flo  = -(A/rho) * exp(-rlo/rho) + 6.0*C/pow(rlo,7);
		if (fmid == 0.0) { rlo = rhi = mid; break; }
		if (fmid * flo < 0.0) {
			rhi = mid;
		} else {
			rlo = mid;
		}
	}
	double rmin = 0.5 * (rlo + rhi);
	return 0.5 * rmin;
}

static double estimate_radius(str_t style, const double* params, int nparams) {
	if (str_eq_cstr_n(style, "lj", 2) && nparams >= 2) {
		double sigma = params[1];
		return lj_radius(sigma);
	}
	if (str_eq_cstr_n(style, "morse", 5) && nparams >= 3) {
		double r0 = params[2];
		return morse_radius(r0);
	}
	if (str_eq_cstr_n(style, "buck", 4) && nparams >= 3) {
		double A   = params[0];
		double rho = params[1];
		double C   = params[2];
		return buckingham_radius(A, rho, C);
	}
	return 0.0; // unknown / unsupported
}

static bool parse_pair_coeffs(md_lammps_atom_type_t types[], size_t num_atom_types, md_buffered_reader_t* reader, str_t hint) {
	str_t tok[8];
	double params[8];

	if (!str_empty(hint)) {
		str_t style = {0};
		bool hybrid = false;
		if (str_eq_cstr_n(hint, "hybrid ", 7)) {
			hybrid = true;
		} else {
			style = hint;
		}

		str_t line;
		for (size_t i = 0; i < num_atom_types; ++i) {
			if (md_buffered_reader_extract_line(&line, reader)) {
				size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);

				int type = (int)parse_int(tok[0]);
				size_t tok_idx = 1;
				if (hybrid) {
					// In hybrid schemes, the style is given as the second parameter
					style = tok[1];
					tok_idx = 2;
				}

				// Extract params
				int num_param = 0;
				for (size_t j = tok_idx; j < num_tok; ++j) {
					params[num_param++] = parse_float(tok[j]);
				}

				if (types[i].id == 0) {
					types[i].id = type;
				} else {
					if (types[i].id != type) {
						MD_LOG_ERROR("Unmatching type in Pair Coeffs: %d, expected %d", type, types[i].id);
						return false;
					}
				}

				if (types[i].radius == 0) {
					types[i].radius = (float)estimate_radius(style, params, num_param);
				}
			}
		}
	}
	return true;
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
	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	return detect_atom_format(&reader);
}

md_lammps_atom_format_t md_lammps_atom_format_from_file(str_t str) {
	md_file_t file = {0};
	if (!md_file_open(&file, str, MD_FILE_READ)) {
		MD_LOG_ERROR("Could not open file '%.*s'", str.len, str.ptr);
		return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
	}

	md_temp_scope_t temp = md_temp_begin();
	const size_t cap = MEGABYTES(1);
	char* buf = md_temp_alloc(temp, cap);

	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
	md_lammps_atom_format_t format = detect_atom_format(&reader);

	md_temp_end(temp);
	md_file_close(&file);

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

	str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

	// Parse headers and sections
	while (md_buffered_reader_extract_line(&line, reader)) {
		const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok > 0 && str_eq(tok[0], STR_LIT("Atoms"))) {
			if (!data->num_atoms) {
				MD_LOG_ERROR("Encountered Atom entries, but number of atoms was not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->atoms, data->num_atoms, alloc);
			if (!parse_atoms(data->atoms, data->num_atoms, reader, mappings)) {
				return false;
			}
			// Sort atoms by id
			qsort(data->atoms, data->num_atoms, sizeof(md_lammps_atom_t), compare_atom);
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Bonds"))) {
			if (!data->num_bonds) {
				MD_LOG_ERROR("Encountered Bond entries, but number of bonds was not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->bonds, data->num_bonds, alloc);
			if (!parse_bonds(data->bonds, data->num_bonds, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Angles"))) {
			if (!data->num_angles) {
				MD_LOG_ERROR("Encountered Angle entries, but number of angles was not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!parse_angles(data->angles, data->num_angles, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Dihedrals"))) {
			if (!data->num_dihedrals) {
				MD_LOG_ERROR("Encountered Dihedral entries, but number of dihedrals was not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->dihedrals, data->num_dihedrals, alloc);
			if (!parse_dihedrals(data->dihedrals, data->num_dihedrals, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Impropers"))) {
			if (!data->num_impropers) {
				MD_LOG_ERROR("Encountered Impropers entries, but number of impropers was not set or zero");
				return false;
			}
			md_buffered_reader_skip_line(reader);
			md_array_resize(data->impropers, data->num_impropers, alloc);
			if (!parse_dihedrals(data->impropers, data->num_impropers, reader)) {
				return false;
			}
		} else if (num_tok > 0 && str_eq(tok[0], STR_LIT("Masses"))) {
			if (!data->num_atom_types) {
				MD_LOG_ERROR("Encountered Mass entries, but number of atom types was not set or zero");
				return false;
			}

			if (!data->atom_types) {
				md_array_resize(data->atom_types, data->num_atom_types, alloc);
				MEMSET(data->atom_types, 0, md_array_bytes(data->atom_types));
			}
			
			md_buffered_reader_skip_line(reader);
			if (parse_masses(data->atom_types, data->num_atom_types, reader) != data->num_atom_types) {
				MD_LOG_ERROR("Number of entries in Masses section did not match the number of atom types");
				return false;
			}
		} else if (num_tok > 1 && str_eq(tok[0], STR_LIT("Pair")) && str_eq(tok[1], STR_LIT("Coeffs"))) {
			if (!data->num_atom_types) {
				MD_LOG_ERROR("Encountered Pair Coeff entries, but number of atom types was not set or zero");
				return false;
			}

			if (!data->atom_types) {
				md_array_resize(data->atom_types, data->num_atom_types, alloc);
				MEMSET(data->atom_types, 0, md_array_bytes(data->atom_types));
			}

			char buf[128];
			str_t hint = {0};
			if (num_tok >= 4 && str_eq(tok[2], STR_LIT("#"))) {
				// Create amalgamation of tokens
				const char* beg = str_beg(tok[3]);
				const char* end = str_end(tok[num_tok-1]);
				str_t str = {beg, end-beg};
				str_copy_to_char_buf(buf, ARRAY_SIZE(buf), str);
				hint = str_from_cstr(buf);
			}

			md_buffered_reader_skip_line(reader);
			if (!parse_pair_coeffs(data->atom_types, data->num_atom_types, reader, hint)) {
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

bool md_lammps_validate_atom_format(char* err_buf, size_t err_cap, const char* format) {
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
	md_file_t file = {0};
	if (md_file_open(&file, filename, MD_FILE_READ)) {
		const size_t cap = MEGABYTES(1);
		md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
		char* buf = md_temp_alloc(temp_scope, cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = md_lammps_data_parse(data, &line_reader, format, alloc);

		md_temp_end(temp_scope);
		md_file_close(&file);
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
	if (data->atom_types) md_array_free(data->atom_types, alloc);
	MEMSET(data, 0, sizeof(md_lammps_data_t));
}

bool md_lammps_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_lammps_data_t* data) {
	ASSERT(sys);
	ASSERT(state);
	ASSERT(data);

	if (!sys->alloc) {
		MD_LOG_ERROR("System allocator not set");
		return false;
	}

	if (!state || !state->alloc) {
		MD_LOG_ERROR("State allocator not set");
		return false;
	}

	md_system_reset(sys);
	md_system_state_init(state, 0);

	const size_t capacity = ROUND_UP(data->num_atoms, 16);

	md_array_ensure(state->xyz,		capacity, state->alloc);
	md_array_ensure(sys->atom.type_idx, capacity, sys->alloc);
	md_array_ensure(sys->atom.flags,    capacity, sys->alloc);

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);

	bool has_resid = (data->num_atoms > 0 && data->atoms[0].resid != -1);

	// Reset atom data and initialize to default value
	sys->atom.type.count = 0;
	md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0, 0, 0, 0, sys->alloc);

	float* type_masses			= md_temp_alloc_array(temp_scope, float, data->num_atom_types);
	md_atomic_number_t* type_z	= md_temp_alloc_array(temp_scope, md_atomic_number_t, data->num_atom_types);
	int* type_map				= md_temp_alloc_array(temp_scope, int, data->num_atom_types);
	for (size_t i = 0; i < data->num_atom_types; ++i) {
		type_masses[i] = data->atom_types[i].mass;
    }
	
	size_t num_successfully_mapped_types = md_util_element_from_mass(type_z, type_masses, data->num_atom_types);

	// @TODO: Do something with the variable above to propagate information of possible Coarse-grained system
	(void)num_successfully_mapped_types;
	
	for (size_t i = 0; i < data->num_atom_types; ++i) {
		char buf[8];
		int len = snprintf(buf, sizeof(buf), "type_%i", data->atom_types[i].id);
		
		str_t type_id = {buf, len};
		float mass   = data->atom_types[i].mass;
		float radius = data->atom_types[i].radius;
		uint32_t color = type_z[i] != MD_Z_X ? md_atomic_number_cpk_color(type_z[i]) : 0xFF808080; // Gray for unknown elements
		
		type_map[i] = md_atom_type_find_or_add(&sys->atom.type, type_id, type_z[i], mass, radius, color, 0, sys->alloc);
	}

	int prev_resid = -1;
	for (size_t i = 0; i < data->num_atoms; ++i) {
		md_atom_type_idx_t type_idx = 0;
		for (size_t j = 0; j < data->num_atom_types; ++j) {
			if (data->atoms[i].type == data->atom_types[j].id) {
				type_idx = (md_atom_type_idx_t)type_map[j];
			}
		}

		int resid = data->atoms[i].resid;
		if (has_resid) {
			if (resid != prev_resid) {
				char buf[8];
				int len = snprintf(buf, sizeof(buf), "comp_%i", resid);
				str_t name = {buf, len};
				md_array_push(sys->component.atom_offset, (uint32_t)sys->atom.count, sys->alloc);
				md_array_push(sys->component.seq_id, resid, sys->alloc);
				md_array_push(sys->component.name, make_label(name), sys->alloc);
				md_array_push(sys->component.flags, 0, sys->alloc);
				sys->component.count += 1;
			}
		}

		md_array_push_no_grow(sys->atom.type_idx, type_idx);
		md_array_push_no_grow(state->xyz, vec3_set(data->atoms[i].x - data->cell.xlo, data->atoms[i].y - data->cell.ylo, data->atoms[i].z - data->cell.zlo));
		md_array_push_no_grow(sys->atom.flags, 0);
		sys->atom.count +=1;

		prev_resid = resid;
	}

	if (has_resid) {
		md_array_push(sys->component.atom_offset, (uint32_t)sys->atom.count, sys->alloc); // Final sentinel
		// No point in trying to infer residue flags as it uses atom names / labels and residue names as hints
	}

	// The per atom charge the 'charge' and 'full' atom styles carry. It was parsed and then dropped;
	// nothing in mdlib reads a partial charge, so it belongs in the attribute table rather than in
	// md_atom_data_t. A style without a q column leaves every value at 0, which is uniform and
	// therefore publishes nothing - so this needs no guard on the atom style.
	{
		md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
		float* charge = (float*)md_temp_alloc(temp, sizeof(float) * MAX(data->num_atoms, (size_t)1));
		if (charge) {
			for (size_t i = 0; i < data->num_atoms; ++i) {
				charge[i] = data->atoms[i].charge;
			}
			md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/charge"), md_unit_none(), 1, charge, data->num_atoms);
		}
		md_temp_end(temp);
	}

	// Create unit cell
	double x = data->cell.xhi - data->cell.xlo;
	double y = data->cell.yhi - data->cell.ylo;
	double z = data->cell.zhi - data->cell.zlo;
	double xy = data->cell.xy;
	double xz = data->cell.xz;
	double yz = data->cell.yz;
    state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);

	// Create bonds
	if (data->num_bonds > 0) {
		md_array_ensure(sys->bond.pairs, data->num_bonds, sys->alloc);
		md_array_ensure(sys->bond.flags, data->num_bonds, sys->alloc);

		for (size_t i = 0; i < data->num_bonds; ++i) {
			int32_t atom_id0 = data->bonds[i].atom_id[0];
			int32_t atom_id1 = data->bonds[i].atom_id[1];
			// LAMMPS atom ids are 1-based
			ASSERT(atom_id0 >= 1 && atom_id0 <= (int32_t)sys->atom.count);
			ASSERT(atom_id1 >= 1 && atom_id1 <= (int32_t)sys->atom.count);
			md_atom_pair_t pair = {
				(uint32_t)(MIN(atom_id0, atom_id1) - 1),
				(uint32_t)(MAX(atom_id0, atom_id1) - 1),
			};
			md_bond_flags_t flag = 0;
            md_array_push_no_grow(sys->bond.pairs, pair);
            md_array_push_no_grow(sys->bond.flags, flag);

			sys->bond.count += 1;
		}
    }

	md_temp_end(temp_scope);

    ASSERT(md_array_size(state->xyz) == sys->atom.count);
    state->num_atoms = sys->atom.count;

	return true;
}

bool md_lammps_system_init_from_str(md_system_t* sys, md_system_state_t* state, str_t str, const char* atom_format) {
	if (!atom_format) {
        md_lammps_atom_format_t format = md_lammps_atom_format_from_str(str);
		if (format == MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
			MD_LOG_ERROR("Could not detect atom format from supplied string");
			return false;
        }
        atom_format = atom_format_string[format];
    }

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

	md_lammps_data_t data = { 0 };
	bool success = md_lammps_data_parse_str(&data, str, atom_format, temp_alloc) && md_lammps_system_init_from_data(sys, state, &data);

	md_lammps_data_free(&data, temp_alloc);
    md_temp_end(temp_scope);
	return success;
}

bool md_lammps_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename, const char* atom_format) {
	if (!atom_format) {
        md_lammps_atom_format_t format = md_lammps_atom_format_from_file(filename);
		if (format == MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
			MD_LOG_ERROR("Could not detect atom format from file '" STR_FMT "'", STR_ARG(filename));
			return false;
        }
        atom_format = atom_format_string[format];
    }

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

	md_lammps_data_t data = { 0 };
	bool success = md_lammps_data_parse_file(&data, filename, atom_format, temp_alloc) && md_lammps_system_init_from_data(sys, state, &data);

	md_lammps_data_free(&data, temp_alloc);
    md_temp_end(temp_scope);
	return success;
}

// TRAJECTORY OPERATIONS

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

// The cell a frame's box bounds describe, and the corner it starts at.
// https://docs.lammps.org/Howto_triclinic.html
static md_unitcell_t lammps_cell_from_bounds(const box_bounds_t* bb, double lo[3]) {
	const double xlo = bb->xlo - MIN(0.0, MIN(bb->xy, MIN(bb->xz, bb->yz)));
	const double xhi = bb->xhi - MAX(0.0, MAX(bb->xy, MAX(bb->xz, bb->yz)));
	const double ylo = bb->ylo - MIN(0.0, bb->yz);
	const double yhi = bb->yhi - MAX(0.0, bb->yz);
	const double zlo = bb->zlo;
	const double zhi = bb->zhi;
	if (lo) {
		lo[0] = xlo;
		lo[1] = ylo;
		lo[2] = zlo;
	}
	return md_unitcell_from_basis_parameters(xhi - xlo, yhi - ylo, zhi - zlo, bb->xy, bb->xz, bb->yz);
}

// One frame's text, from its ITEM: TIMESTEP on: the cell, and the coordinates in id order written
// to out_x[i * stride] and so on - planar with stride 1, packed with the three pointers one apart and
// stride 3.
static bool lammps_decode_frame_text(const coord_mappings_t* mappings, str_t str, size_t* out_num_atoms, md_unitcell_t* out_cell, float* out_x, float* out_y, float* out_z, size_t stride) {
	bool result = false;
	str_t tokens[32];
	const bool output_coords = out_x != NULL && out_y != NULL && out_z != NULL;

	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	str_t line;

	header_t header;
	if (!parse_header(&header, &reader)) {
		MD_LOG_ERROR("Could not parse header");
		return false;
	}

	double lo[3];
	const md_unitcell_t cell = lammps_cell_from_bounds(&header.box_bounds, lo);
	const double xlen = cell.x;
	const double ylen = cell.y;
	const double zlen = cell.z;

	// transform matrix to apply
	mat4_t M = mat4_translate(-(float)lo[0], -(float)lo[1], -(float)lo[2]);
	if (mappings->flags & COORD_FLAG_SCALED) {
		// Scaling
		mat3_t A;
		md_unitcell_A_extract_float(A.elem, &cell);
		M = mat4_from_mat3(A);
	}

	if (output_coords) {
		if (!md_buffered_reader_extract_line(&line, &reader) || !str_eq_cstr_n(line, "ITEM: ATOMS", 11)) {
			MD_LOG_ERROR("Expected ITEM: ATOMS after header");
			goto done;
		}
		size_t line_count = 0;
		size_t expected_num_tokens = mappings->num_coord_tokens;
		bool coords_result = false;

		// We need to store the coordinates in a temporary buffer since we need to sort them by id
		md_temp_scope_t temp_scope = md_temp_begin();
		id_xyz_t* id_xyz = md_temp_alloc_array(temp_scope, id_xyz_t, header.num_atoms);

		while (md_buffered_reader_extract_line(&line, &reader) && line_count < header.num_atoms) {
			size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
			if (num_tokens != expected_num_tokens) {
				MD_LOG_ERROR("Unexpected number of tokens in ITEM: ATOMS line, got %zu, expected %zu", num_tokens, expected_num_tokens);
				goto coords_done;
			}

			int32_t id = (int32_t)parse_int(tokens[mappings->id_idx]);
			vec4_t coord = {
				(float)parse_float(tokens[mappings->coord_idx[0]]),
				(float)parse_float(tokens[mappings->coord_idx[1]]),
				(float)parse_float(tokens[mappings->coord_idx[2]]),
				1.0f
			};
			coord = mat4_mul_vec4(M, coord);

			if (mappings->flags & COORD_FLAG_UNWRAP) {
				int64_t ix = parse_int(tokens[mappings->image_idx[0]]);
				int64_t iy = parse_int(tokens[mappings->image_idx[1]]);
				int64_t iz = parse_int(tokens[mappings->image_idx[2]]);
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
			out_x[i * stride] = id_xyz[i].x;
			out_y[i * stride] = id_xyz[i].y;
			out_z[i * stride] = id_xyz[i].z;
		}
		coords_result = (line_count == header.num_atoms);
		if (!coords_result) {
			MD_LOG_ERROR("The frame holds %zu atom lines of the %zu it declares", line_count, header.num_atoms);
		}

	coords_done:
		md_temp_end(temp_scope);
		if (!coords_result) goto done;
	}

	if (out_num_atoms) {
		*out_num_atoms = header.num_atoms;
	}

	if (out_cell) {
		*out_cell = cell;
	}

	result = true;
done:
	return result;
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
		md_array_push(cache->frame_steps, header.timestep, alloc);
		{
			const md_unitcell_t cell = lammps_cell_from_bounds(&header.box_bounds, NULL);
			float A[3][3];
			md_unitcell_A_extract_float(A, &cell);
			md_array_push_array(cache->frame_cells, &A[0][0], 9, alloc);
		}
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

static bool lammps_trajectory_parse_file(lammps_cache_t* cache, str_t filename, struct md_allocator_i* alloc, struct md_allocator_i* avoid) {
	bool result = false;
	md_file_t file = {0};
	if (md_file_open(&file, filename, MD_FILE_READ)) {
		const int64_t cap = MEGABYTES(1);
		// Clear of what alloc was made from, not of alloc itself: an arena on a temp allocator would
		// otherwise be rewound with the scratch.
		md_temp_scope_t temp_scope = md_temp_begin_avoid(avoid);
		char* buf = md_temp_alloc(temp_scope, cap);

		md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
		result = lammps_trajectory_parse(cache, &line_reader, alloc);

		md_temp_end(temp_scope);
		md_file_close(&file);
	}
	else {
		MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
	}
	return result;
}


// The cache beside path when it was made from the file as it is now: header, num_frames + 1 offsets,
// num_frames steps, the coord mappings, 9 cell floats per frame.
static bool try_read_cache(lammps_cache_t* cache, str_t path, md_allocator_i* alloc) {
	ASSERT(cache);
	ASSERT(alloc);
	md_file_t file = {0};
	if (!md_run_cache_open(&file, &cache->header, path, MD_LAMMPS_CACHE_MAGIC, MD_LAMMPS_CACHE_VERSION)) {
		return false;
	}
	const size_t n = cache->header.num_frames;
	md_array_resize(cache->frame_offsets, n + 1, alloc);
	md_array_resize(cache->frame_steps,   n,     alloc);
	md_array_resize(cache->frame_cells,   n * 9, alloc);
	const bool ok =
		md_file_read(file, cache->frame_offsets, (n + 1) * sizeof(int64_t)) == (n + 1) * sizeof(int64_t) &&
		md_file_read(file, cache->frame_steps, n * sizeof(int64_t)) == n * sizeof(int64_t) &&
		md_file_read(file, &cache->coord_mappings, sizeof(cache->coord_mappings)) == sizeof(cache->coord_mappings) &&
		cache->coord_mappings.flags != COORD_FLAG_NONE &&
		md_file_read(file, cache->frame_cells, n * 9 * sizeof(float)) == n * 9 * sizeof(float) &&
		md_file_tell(file) == (int64_t)md_file_size(file);
	if (!ok) {
		MD_LOG_ERROR("The LAMMPS cache beside '" STR_FMT "' is incomplete", STR_ARG(path));
		md_array_free(cache->frame_offsets, alloc);
		md_array_free(cache->frame_steps,   alloc);
		md_array_free(cache->frame_cells,   alloc);
		cache->frame_offsets = NULL;
		cache->frame_steps   = NULL;
		cache->frame_cells   = NULL;
	}
	md_file_close(&file);
	return ok;
}

static bool write_cache(const lammps_cache_t* cache, str_t path, const md_file_info_t* scanned) {
	const size_t n = cache->header.num_frames;
	if (md_array_size(cache->frame_offsets) != n + 1 || md_array_size(cache->frame_steps) != n || md_array_size(cache->frame_cells) != n * 9) {
		MD_LOG_ERROR("The LAMMPS index of '" STR_FMT "' is inconsistent; no cache is written", STR_ARG(path));
		return false;
	}
	md_file_t file = {0};
	if (!md_run_cache_create(&file, path, scanned, MD_LAMMPS_CACHE_MAGIC, MD_LAMMPS_CACHE_VERSION, cache->header.num_atoms, n)) {
		return false;
	}
	const bool ok =
		md_file_write(file, cache->frame_offsets, (n + 1) * sizeof(int64_t)) == (n + 1) * sizeof(int64_t) &&
		md_file_write(file, cache->frame_steps, n * sizeof(int64_t)) == n * sizeof(int64_t) &&
		md_file_write(file, &cache->coord_mappings, sizeof(cache->coord_mappings)) == sizeof(cache->coord_mappings) &&
		md_file_write(file, cache->frame_cells, n * 9 * sizeof(float)) == n * 9 * sizeof(float);
	if (!ok) {
		MD_LOG_ERROR("Failed to write the LAMMPS cache beside '" STR_FMT "'", STR_ARG(path));
	}
	md_file_close(&file);
	return ok;
}

// Where each frame is, its step and its box, and how atom lines are laid out: from the cache beside
// the file when that is current, from a scan of the file otherwise (writing the cache unless told
// not to). From alloc; the scan's scratch keeps clear of avoid, which is what alloc was made from.
static bool lammps_index_load(lammps_cache_t* cache, str_t filename, md_run_flags_t flags, md_allocator_i* alloc, md_allocator_i* avoid) {
	MEMSET(cache, 0, sizeof(*cache));

	if (try_read_cache(cache, filename, alloc)) {
		return true;
	}
	MEMSET(cache, 0, sizeof(*cache));

	// The file as it is before the scan is what the cache is stamped with
	md_file_info_t scanned = {0};
	if (!md_file_info_extract_from_path(filename, &scanned)) {
		MD_LOG_ERROR("Failed to open file for LAMMPS trajectory");
		return false;
	}
	if (!lammps_trajectory_parse_file(cache, filename, alloc, avoid)) {
		MD_LOG_ERROR("LAMMPS trajectory could not be read from file");
		return false;
	}

	if (!(flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
		// A cache that cannot be written only costs the next load a scan
		write_cache(cache, filename, &scanned);
	}
	return true;
}

// ### RUN ###

enum {
	LAMMPS_LAYOUT_ID, LAMMPS_LAYOUT_X, LAMMPS_LAYOUT_Y, LAMMPS_LAYOUT_Z,
	LAMMPS_LAYOUT_IX, LAMMPS_LAYOUT_IY, LAMMPS_LAYOUT_IZ,
	LAMMPS_LAYOUT_NUM_TOKENS, LAMMPS_LAYOUT_FLAGS, LAMMPS_LAYOUT_COUNT
};

// <run>/atom/position: one read of the frame's text, the atoms put in id order into the caller's
// buffer. Atom lines are in no particular order in a dump, so one atom costs the whole frame.
static size_t lammps_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
	const md_system_t* sys = (const md_system_t*)user_data;
	ASSERT(sys);
	if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

	md_run_source_t src;
	if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;

	char buf[512];
	const md_attribute_t* layout_attr = md_attributes_find(&sys->attributes, md_run_path(buf, sizeof(buf), src.run, STR_LIT("source/layout")));
	if (!layout_attr || layout_attr->format.type != MD_ATTRIBUTE_TYPE_I32 || !layout_attr->data || md_attribute_value_count(&layout_attr->format) != LAMMPS_LAYOUT_COUNT) {
		MD_LOG_ERROR("LAMMPS: the run '" STR_FMT "' has lost its layout", STR_ARG(src.run));
		return 0;
	}
	const int32_t* layout = (const int32_t*)layout_attr->data;
	coord_mappings_t mappings = {0};
	mappings.id_idx           = (int8_t)layout[LAMMPS_LAYOUT_ID];
	mappings.coord_idx[0]     = (int8_t)layout[LAMMPS_LAYOUT_X];
	mappings.coord_idx[1]     = (int8_t)layout[LAMMPS_LAYOUT_Y];
	mappings.coord_idx[2]     = (int8_t)layout[LAMMPS_LAYOUT_Z];
	mappings.image_idx[0]     = (int8_t)layout[LAMMPS_LAYOUT_IX];
	mappings.image_idx[1]     = (int8_t)layout[LAMMPS_LAYOUT_IY];
	mappings.image_idx[2]     = (int8_t)layout[LAMMPS_LAYOUT_IZ];
	mappings.num_coord_tokens = (int8_t)layout[LAMMPS_LAYOUT_NUM_TOKENS];
	mappings.flags            = (int8_t)layout[LAMMPS_LAYOUT_FLAGS];

	const uint32_t frame = slice->idx[0];
	const size_t N = attr->format.shape[1];
	size_t first = 0, count = N;
	if (slice->num_idx == 2) {
		if (slice->idx[1] >= N) return 0;
		first = slice->idx[1];
		count = 1;
	}
	if (cap != count * 3) return 0;

	const size_t frame_size = (size_t)src.size[frame];
	md_temp_scope_t temp = md_temp_begin();
	size_t written = 0;
	char*  text = md_temp_alloc(temp, MAX(frame_size, 1));
	float* xyz  = (count == N) ? (float*)dst : md_temp_alloc(temp, N * 3 * sizeof(float));
	if (text && xyz && md_attribute_io_read_at(io, src.path, src.offset[frame], text, frame_size) == frame_size) {
		size_t num_atoms = 0;
		if (lammps_decode_frame_text(&mappings, (str_t){ text, frame_size }, &num_atoms, NULL, xyz + 0, xyz + 1, xyz + 2, 3) && num_atoms == N) {
			if (count != N) {
				MEMCPY(dst, xyz + first * 3, 3 * sizeof(float));
			}
			written = cap;
		}
	} else {
		MD_LOG_ERROR("LAMMPS: Failed to read frame %u from '" STR_FMT "'", frame, STR_ARG(src.path));
	}
	md_temp_end(temp);
	return written;
}

bool md_lammps_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
	ASSERT(sys);
	char path_buf[4096];
	const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
	const str_t path = { path_buf, path_len };

	md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
	bool result = false;

	lammps_cache_t index;
	if (path_len == 0 || !lammps_index_load(&index, path, flags, arena, md_get_heap_allocator())) {
		goto done;
	}

	const size_t F = index.header.num_frames;
	double*  times = md_alloc(arena, F * sizeof(double));
	int64_t* sizes = md_alloc(arena, F * sizeof(int64_t));
	for (size_t i = 0; i < F; ++i) {
		times[i] = (double)i;
		sizes[i] = index.frame_offsets[i + 1] - index.frame_offsets[i];
	}

	const md_attribute_virtual_t virt = { .provider = lammps_position_provider, .user_data = sys };
	const md_run_desc_t desc = {
		.num_frames    = F,
		.num_atoms     = index.header.num_atoms,
		.time          = times,
		.time_unit     = md_unit_none(),    // a dump has TIMESTEP but never dt: ordinals
		.step          = index.frame_steps,
		.unitcell      = index.frame_cells,
		.source_path   = path,
		.source_offset = index.frame_offsets,
		.source_size   = sizes,
		.position_virt = &virt,
	};
	if (!md_run_publish(sys, run, &desc)) {
		goto done;
	}

	const coord_mappings_t* m = &index.coord_mappings;
	const int32_t layout[LAMMPS_LAYOUT_COUNT] = {
		[LAMMPS_LAYOUT_ID] = m->id_idx,
		[LAMMPS_LAYOUT_X]  = m->coord_idx[0], [LAMMPS_LAYOUT_Y]  = m->coord_idx[1], [LAMMPS_LAYOUT_Z]  = m->coord_idx[2],
		[LAMMPS_LAYOUT_IX] = m->image_idx[0], [LAMMPS_LAYOUT_IY] = m->image_idx[1], [LAMMPS_LAYOUT_IZ] = m->image_idx[2],
		[LAMMPS_LAYOUT_NUM_TOKENS] = m->num_coord_tokens,
		[LAMMPS_LAYOUT_FLAGS]      = m->flags,
	};
	char buf[512];
	if (!md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
		.path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/layout")),
		.format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { LAMMPS_LAYOUT_COUNT } },
		.unit = md_unit_none(),
		.description = STR_INIT("Columns of an atom line: id, x, y, z, ix, iy, iz; token count; cartesian, scaled, unwrapped"),
		.data = layout, .byte_size = sizeof(layout)})) {
		md_attributes_remove_prefix(&sys->attributes, run);
		goto done;
	}
	result = true;

done:
	md_arena_allocator_destroy(arena);
	return result;
}

