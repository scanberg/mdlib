#include "md_lammps.h"

#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_hash.h>

#define MD_LAMMPS_TRAJ_MAGIC 0x2312ad7b78a9bc20
#define MD_LAMMPS_TRAJ_READER_MAGIC 0x2312ad7b78a9bc21
#define MD_LAMMPS_CACHE_MAGIC 0x89172bab
#define MD_LAMMPS_CACHE_VERSION 19
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
	COORD_FLAG_UNWRAP    = 4,	// Wrapped coordinates (x, xs) with image flags (ix iy iz): unwrapped by them
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
	"amoeba",
	"bpm/sphere",
	"dielectric",
	"hybrid",
	"rheo",
	"rheo/thermal",
	"sph",
	"spin",
	"tdpd",
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
	"id resid type q x y z", // amoeba
	"id resid type diameter density x y z", // bpm/sphere
	"id resid type q x y z", // dielectric (mux muy muz area ed em epsilon curvature follow)
	"id type x y z", // hybrid (the columns of the sub-styles follow)
	"id type status rho x y z", // rheo
	"id type status rho energy x y z", // rheo/thermal
	"id type rho esph cv x y z", // sph
	"id type x y z", // spin (spx spy spz sp follow)
	"id type x y z", // tdpd (cc1 ... follow)
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
	MD_MEMSET(mappings, -1, sizeof(int) * ATOM_FIELD_COUNT);

	str_t str = str_trim(str_from_cstr(atom_format));
	if (str_empty(str)) {
		if (err_buf) {
			md_snprintf(err_buf, err_cap, "Lammps atom format: Empty format string");
		} else {
			MD_LOG_ERROR("Lammps atom format: Empty format string");
		}
		return 0;
	}

	str_t tokens[16];
	const size_t num_tokens = extract_tokens(tokens, MD_COUNTOF(tokens), &str);
	for (size_t i = 0; i < num_tokens; ++i) {
		const str_t token = tokens[i];
		for (size_t j = 1; j < ATOM_FIELD_COUNT; ++j) {
			if (str_eq(token, str_from_cstr(atom_field_name[j]))) {
				if (mappings[j] != -1) {
					if (err_buf) {
						md_snprintf(err_buf, err_cap, "Lammps atom format: Contains duplicate field name '"STR_FMT"'", STR_ARG(token));
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
	for (size_t i = 0; i < MD_COUNTOF(atom_field_required); ++i) {
		int idx = atom_field_required[i];
		if (mappings[idx] == -1) {
			if (err_buf) {
				md_snprintf(err_buf, err_cap, "Lammps atom format: Missing required identifier '%s'", atom_field_name[idx]);
			} else {
				MD_LOG_ERROR("Lammps atom format: Missing required identifier '%s'", atom_field_name[idx]);
			}
			return 0;
		}
	}

	return num_tokens;
}

// ### DATA FILE ###
//
// A data file is a title line, header lines ("3 atoms", "2 atom types", "0.0 10.0 xlo xhi", ...), then sections, each a
// keyword line ("Atoms", "Masses", "Pair Coeffs", ...), a blank line and a line for each of its entries. A '#' at the
// start of a line or after a blank begins a comment, which a section keyword line may carry to say how its entries
// are laid out ("Atoms # full", "Pair Coeffs # lj/cut"). An atom type is a number, or a label of the Atom Type Labels
// section where it has one. https://docs.lammps.org/read_data.html

enum {
	SEC_NONE,
	SEC_ATOMS,
	SEC_MASSES,
	SEC_PAIR_COEFFS,
	SEC_BONDS,
	SEC_ANGLES,
	SEC_DIHEDRALS,
	SEC_IMPROPERS,
	SEC_ATOM_TYPE_LABELS,
	SEC_OTHER,		// Of no interest here: skipped
};

static const struct {
	const char* keyword;
	int section;
} lammps_sections[] = {
	{ "Atoms",						SEC_ATOMS },
	{ "Masses",						SEC_MASSES },
	{ "Pair Coeffs",				SEC_PAIR_COEFFS },
	{ "Bonds",						SEC_BONDS },
	{ "Angles",						SEC_ANGLES },
	{ "Dihedrals",					SEC_DIHEDRALS },
	{ "Impropers",					SEC_IMPROPERS },
	{ "Atom Type Labels",			SEC_ATOM_TYPE_LABELS },
	{ "Velocities",					SEC_OTHER },
	{ "Ellipsoids",					SEC_OTHER },
	{ "Lines",						SEC_OTHER },
	{ "Triangles",					SEC_OTHER },
	{ "Bodies",						SEC_OTHER },
	{ "PairIJ Coeffs",				SEC_OTHER },
	{ "Bond Coeffs",				SEC_OTHER },
	{ "Angle Coeffs",				SEC_OTHER },
	{ "Dihedral Coeffs",			SEC_OTHER },
	{ "Improper Coeffs",			SEC_OTHER },
	{ "BondBond Coeffs",			SEC_OTHER },
	{ "BondAngle Coeffs",			SEC_OTHER },
	{ "MiddleBondTorsion Coeffs",	SEC_OTHER },
	{ "EndBondTorsion Coeffs",		SEC_OTHER },
	{ "AngleTorsion Coeffs",		SEC_OTHER },
	{ "AngleAngleTorsion Coeffs",	SEC_OTHER },
	{ "BondBond13 Coeffs",			SEC_OTHER },
	{ "AngleAngle Coeffs",			SEC_OTHER },
	{ "Bond Type Labels",			SEC_OTHER },
	{ "Angle Type Labels",			SEC_OTHER },
	{ "Dihedral Type Labels",		SEC_OTHER },
	{ "Improper Type Labels",		SEC_OTHER },
	{ "CMAP",						SEC_OTHER },
};

// The content of a line of a data file: up to its comment ('#' at the start or after a blank), trimmed. Its comment,
// after the '#' and trimmed, in *comment where given.
static str_t lammps_content(str_t line, str_t* comment) {
	const char* beg = line.ptr;
	const char* end = line.ptr + line.len;
	const char* c = beg;
	while (c < end) {
		c = (const char*)memchr(c, '#', (size_t)(end - c));
		if (!c) { c = end; break; }
		if (c == beg || is_whitespace(c[-1])) break;
		++c;
	}
	if (comment) {
		*comment = c < end ? str_trim((str_t){ c + 1, (size_t)(end - c - 1) }) : (str_t){0};
	}
	return str_trim((str_t){ beg, (size_t)(c - beg) });
}

// The section a keyword line begins (its content, the words of the keyword a blank apart), SEC_NONE for any other line
static int lammps_section_kind(str_t content) {
	if (content.len == 0 || !is_alpha(content.ptr[0])) return SEC_NONE;
	for (size_t i = 0; i < MD_COUNTOF(lammps_sections); ++i) {
		const char* k = lammps_sections[i].keyword;
		const char* c = content.ptr;
		const char* e = content.ptr + content.len;
		while (*k && c < e) {
			if (*k == ' ') {
				if (!is_whitespace(*c)) break;
				while (c < e && is_whitespace(*c)) ++c;
				++k;
			} else if (*k == *c) {
				++k;
				++c;
			} else {
				break;
			}
		}
		if (*k == '\0' && c == e) return lammps_sections[i].section;
	}
	return SEC_NONE;
}

// The next entry of a section: the content of its next line which has any, false at the end of the file or at the
// keyword line of the next section (which is left to be read).
static bool lammps_section_line(md_reader_i* reader, str_t* content) {
	str_t line;
	const char* next;
	while (md_reader_find_line(reader, &line, &next)) {
		const str_t c = lammps_content(line, NULL);
		if (c.len == 0) {
			reader->cur = next;
			continue;
		}
		if (is_alpha(c.ptr[0]) && lammps_section_kind(c) != SEC_NONE) {
			return false;
		}
		reader->cur = next;
		*content = c;
		return true;
	}
	return false;
}

// The atom type of a token: its number, or the type its label is (of labels, a map of the hash of the label to the
// type), 0 when it is neither
static int32_t lammps_atom_type(str_t tok, const md_hashmap32_t* labels) {
	int64_t v = 0;
	if (md_parse_i64(&v, tok) == tok.len) return (int32_t)v;
	if (labels && labels->num_buckets) {
		const uint32_t* found = md_hashmap_get(labels, md_hash64(tok.ptr, tok.len, 0));
		if (found) return (int32_t)*found;
	}
	return 0;
}

// The highest column of the fields the mappings map, plus one: the tokens of an atom line which are read
static size_t lammps_atom_columns(const int mappings[ATOM_FIELD_COUNT]) {
	int max = -1;
	for (size_t i = 1; i < ATOM_FIELD_COUNT; ++i) max = MD_MAX(max, mappings[i]);
	return (size_t)(max + 1);
}

// Whether an atom line (its content) is of the atom format of mappings: the number columns are numbers, the type a
// number within [1, num_types] (when num_types is known) or a label. Of exactly the columns of the format, or of
// three image flags more, when exact; at least as many otherwise.
static bool lammps_atom_line_fits(str_t content, const int mappings[ATOM_FIELD_COUNT], size_t num_fields, size_t num_types, bool exact) {
	str_t tok[24];
	const size_t n = extract_tokens(tok, MD_COUNTOF(tok), &content);
	if (n < num_fields) return false;
	if (exact) {
		if (n != num_fields && n != num_fields + 3) return false;
		if (n == num_fields + 3 && !(is_int(tok[n - 3]) && is_int(tok[n - 2]) && is_int(tok[n - 1]))) return false;
	}
	for (size_t j = 1; j < ATOM_FIELD_COUNT; ++j) {
		if (mappings[j] < 0) continue;
		const str_t t = tok[mappings[j]];
		if (j == ATOM_FIELD_TYPE) {
			if (is_int(t)) {
				const int64_t type = parse_int(t);
				if (type < 1 || (num_types && (size_t)type > num_types)) return false;
			} else if (!is_alpha(t.ptr[0])) {
				return false;
			}
		} else if (atom_field_type[j] == TYPE_INT) {
			if (!is_int(t)) return false;
		} else if (atom_field_type[j] == TYPE_FLOAT) {
			if (!is_float(t)) return false;
		}
	}
	return true;
}

// The atom format of a data file: the one its Atoms keyword line names ("Atoms # full", "Atoms # hybrid sphere
// dipole"), else the first of full, molecular, charge and atomic which the first lines of the section fit exactly
// (with or without image flags).
static md_lammps_atom_format_t detect_atom_format(md_reader_i* reader) {
	str_t line;
	str_t atom_lines[8];
	char atom_line_buf[8][512];
	size_t atom_line_count = 0;
	char hint_buf[64] = {0};
	str_t hint = {0};
	size_t num_types = 0;
	bool header = true;

	while (md_reader_extract_line(reader, &line)) {
		str_t comment;
		const str_t content = lammps_content(line, &comment);
		const int sec = lammps_section_kind(content);
		if (sec == SEC_NONE) {
			if (header) {
				str_t tok[4];
				str_t c = content;
				if (extract_tokens(tok, MD_COUNTOF(tok), &c) == 3 && str_eq(tok[1], STR_LIT("atom")) && str_eq(tok[2], STR_LIT("types"))) {
					num_types = (size_t)MD_MAX(0, parse_int(tok[0]));
				}
			}
		} else if (sec == SEC_ATOMS) {
			str_t word;
			if (extract_token(&word, &comment)) {
				const size_t len = str_copy_to_char_buf(hint_buf, sizeof(hint_buf), word);
				hint = (str_t){ hint_buf, len };
			}
			// Copies: a file reader may move what it has read as it reads more
			str_t atom_line;
			while (atom_line_count < MD_COUNTOF(atom_lines) && lammps_section_line(reader, &atom_line)) {
				const size_t len = str_copy_to_char_buf(atom_line_buf[atom_line_count], sizeof(atom_line_buf[0]), atom_line);
				atom_lines[atom_line_count] = (str_t){ atom_line_buf[atom_line_count], len };
				atom_line_count += 1;
			}
			break;
		} else {
			header = false;
		}

		if (md_reader_tell(reader) > MD_MEGABYTES(64)) {
			MD_LOG_ERROR("Failed to detect atom format, could not find 'Atoms' entry");
			return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
		}
	}

	if (atom_line_count == 0) {
		MD_LOG_ERROR("Failed to detect atom format, the file has no Atoms section");
		return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
	}

	md_lammps_atom_format_t candidates[4] = { MD_LAMMPS_ATOM_FORMAT_FULL, MD_LAMMPS_ATOM_FORMAT_MOLECULAR, MD_LAMMPS_ATOM_FORMAT_CHARGE, MD_LAMMPS_ATOM_FORMAT_ATOMIC };
	size_t num_candidates = MD_COUNTOF(candidates);
	bool exact = true;
	if (!str_empty(hint)) {
		md_lammps_atom_format_t named = MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
		for (size_t i = 1; i < MD_LAMMPS_ATOM_FORMAT_COUNT; ++i) {
			if (str_eq_cstr(hint, atom_format_name[i])) {
				named = (md_lammps_atom_format_t)i;
				break;
			}
		}
		if (named == MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
			MD_LOG_ERROR("LAMMPS: Unknown atom style '" STR_FMT "' in the Atoms section", STR_ARG(hint));
			return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
		}
		candidates[0] = named;
		num_candidates = 1;
		exact = false;
	}

	for (size_t k = 0; k < num_candidates; ++k) {
		int mappings[ATOM_FIELD_COUNT];
		const size_t num_fields = interpret_format(mappings, atom_format_string[candidates[k]], 0, 0);
		bool fits = num_fields > 0;
		for (size_t i = 0; i < atom_line_count && fits; ++i) {
			fits = lammps_atom_line_fits(atom_lines[i], mappings, num_fields, num_types, exact);
		}
		if (fits) {
			return candidates[k];
		}
	}

	if (!str_empty(hint)) {
		MD_LOG_ERROR("LAMMPS: The lines of the Atoms section are not of the atom style '" STR_FMT "' it names", STR_ARG(hint));
	} else {
		MD_LOG_ERROR("LAMMPS: Could not tell the atom style of the Atoms section (name it as 'Atoms # <style>')");
	}
	return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
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

// The radii of the types of the Pair Coeffs section, from the pair style its keyword line names ("Pair Coeffs # lj/cut",
// "Pair Coeffs # hybrid lj/cut morse", which names it on each line): where the type has none from its mass
static bool parse_pair_coeffs(md_lammps_atom_type_t types[], size_t num_types, md_reader_i* reader, str_t hint, const md_hashmap32_t* labels) {
	str_t tok[8];
	double params[8];

	str_t style = hint;
	const bool hybrid = str_eq_cstr_n(hint, "hybrid", 6) && (hint.len == 6 || is_whitespace(hint.ptr[6]));

	str_t line;
	for (size_t i = 0; i < num_types; ++i) {
		if (!lammps_section_line(reader, &line)) {
			MD_LOG_ERROR("LAMMPS: The Pair Coeffs section ends after %zu of %zu atom types", i, num_types);
			return false;
		}
		if (str_empty(hint)) continue;

		const size_t num_tok = extract_tokens(tok, MD_COUNTOF(tok), &line);
		const int32_t type = lammps_atom_type(tok[0], labels);
		if (type < 1 || (size_t)type > num_types) {
			MD_LOG_ERROR("LAMMPS: Invalid atom type '" STR_FMT "' in Pair Coeffs", STR_ARG(tok[0]));
			return false;
		}
		size_t tok_idx = 1;
		if (hybrid) {
			// In hybrid schemes, the style is given as the second parameter
			if (num_tok < 2) continue;
			style = tok[1];
			tok_idx = 2;
		}

		int num_param = 0;
		for (size_t j = tok_idx; j < num_tok; ++j) {
			params[num_param++] = parse_float(tok[j]);
		}

		md_lammps_atom_type_t* t = &types[type - 1];
		if (t->radius == 0) {
			t->radius = (float)estimate_radius(style, params, num_param);
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
	md_reader_i reader = md_reader_from_str(str);
	return detect_atom_format(&reader);
}

md_lammps_atom_format_t md_lammps_atom_format_from_file(str_t str) {
	md_file_t file = {0};
	if (!md_file_open(&file, str, MD_FILE_READ)) {
		MD_LOG_ERROR("Could not open file '%.*s'", str.len, str.ptr);
		return MD_LAMMPS_ATOM_FORMAT_UNKNOWN;
	}

	md_temp_scope_t temp = md_temp_begin();
	const size_t cap = MD_MEGABYTES(1);

	md_file_reader_t reader;
	md_file_reader_init(&reader, file, cap, md_temp_allocator(temp));
	md_lammps_atom_format_t format = detect_atom_format(&reader.r);

	md_file_reader_close(&reader);
	md_temp_end(temp);
	md_file_close(&file);

	return format;
}

static int compare_atom(const void* a, const void* b) {
	const md_lammps_atom_t* atom_a = (const md_lammps_atom_t*)a;
	const md_lammps_atom_t* atom_b = (const md_lammps_atom_t*)b;
	return (atom_a->id > atom_b->id) - (atom_a->id < atom_b->id);
}

// A header line (its content), into data. False for a malformed one of the box.
static bool lammps_header_line(md_lammps_data_t* data, str_t content) {
	str_t tok[8];
	const size_t num_tok = extract_tokens(tok, MD_COUNTOF(tok), &content);
	if (num_tok == 2 && is_int(tok[0])) {
		const size_t n = (size_t)MD_MAX(0, parse_int(tok[0]));
		if      (str_eq(tok[1], STR_LIT("atoms")))		data->num_atoms = n;
		else if (str_eq(tok[1], STR_LIT("bonds")))		data->num_bonds = n;
		else if (str_eq(tok[1], STR_LIT("angles")))		data->num_angles = n;
		else if (str_eq(tok[1], STR_LIT("dihedrals")))	data->num_dihedrals = n;
		else if (str_eq(tok[1], STR_LIT("impropers")))	data->num_impropers = n;
	} else if (num_tok == 3 && str_eq(tok[2], STR_LIT("types")) && is_int(tok[0])) {
		const size_t n = (size_t)MD_MAX(0, parse_int(tok[0]));
		if      (str_eq(tok[1], STR_LIT("atom")))		data->num_particle_types = n;
		else if (str_eq(tok[1], STR_LIT("bond")))		data->num_bond_types = n;
		else if (str_eq(tok[1], STR_LIT("angle")))		data->num_angle_types = n;
		else if (str_eq(tok[1], STR_LIT("dihedral")))	data->num_dihedral_types = n;
		else if (str_eq(tok[1], STR_LIT("improper")))	data->num_improper_types = n;
	} else if (num_tok == 4 && (str_eq(tok[2], STR_LIT("xlo")) || str_eq(tok[2], STR_LIT("ylo")) || str_eq(tok[2], STR_LIT("zlo")))) {
		const char axis = tok[2].ptr[0];
		if (tok[3].len != 3 || tok[3].ptr[0] != axis || !str_eq_cstr_n(str_substr(tok[3], 1, 2), "hi", 2)) return true;
		if (!is_float(tok[0]) || !is_float(tok[1])) {
			MD_LOG_ERROR("Failed to parse cell definition");
			return false;
		}
		const float lo = (float)parse_float(tok[0]);
		const float hi = (float)parse_float(tok[1]);
		if (axis == 'x') { data->cell.xlo = lo; data->cell.xhi = hi; }
		if (axis == 'y') { data->cell.ylo = lo; data->cell.yhi = hi; }
		if (axis == 'z') { data->cell.zlo = lo; data->cell.zhi = hi; }
	} else if (num_tok == 6 && str_eq(tok[3], STR_LIT("xy")) && str_eq(tok[4], STR_LIT("xz")) && str_eq(tok[5], STR_LIT("yz"))) {
		if (!is_float(tok[0]) || !is_float(tok[1]) || !is_float(tok[2])) {
			MD_LOG_ERROR("Failed to parse cell definition");
			return false;
		}
		data->cell.xy = (float)parse_float(tok[0]);
		data->cell.xz = (float)parse_float(tok[1]);
		data->cell.yz = (float)parse_float(tok[2]);
	} else if ((num_tok == 4 && (str_eq(tok[3], STR_LIT("avec")) || str_eq(tok[3], STR_LIT("bvec")) || str_eq(tok[3], STR_LIT("cvec")))) ||
			   (num_tok == 5 && str_eq(tok[3], STR_LIT("abc")) && str_eq(tok[4], STR_LIT("origin")))) {
		// A general triclinic box
		if (!is_float(tok[0]) || !is_float(tok[1]) || !is_float(tok[2])) {
			MD_LOG_ERROR("Failed to parse cell definition");
			return false;
		}
		float* v = str_eq(tok[3], STR_LIT("avec")) ? data->cell.avec : str_eq(tok[3], STR_LIT("bvec")) ? data->cell.bvec : str_eq(tok[3], STR_LIT("cvec")) ? data->cell.cvec : data->cell.origin;
		v[0] = (float)parse_float(tok[0]);
		v[1] = (float)parse_float(tok[1]);
		v[2] = (float)parse_float(tok[2]);
		data->cell.general = true;
	}
	return true;
}

// The types, of the header's number of them: type i + 1 at i, where the sections of types put what they hold
static void lammps_types_init(md_lammps_data_t* data, md_allocator_i* alloc) {
	if (data->particle_types || data->num_particle_types == 0) return;
	md_array_resize(data->particle_types, data->num_particle_types, alloc);
	MD_MEMSET(data->particle_types, 0, md_array_bytes(data->particle_types));
	for (size_t i = 0; i < data->num_particle_types; ++i) {
		data->particle_types[i].id = (int32_t)(i + 1);
	}
}

// The ids of the atoms of the lines of a section of n lines, and its type: tokens 0 and 1, then the atoms
static bool lammps_parse_topology(int32_t* out, size_t stride, size_t num_atom_ids, size_t n, md_reader_i* reader, const char* section) {
	str_t tok[8];
	str_t line;
	for (size_t i = 0; i < n; ++i) {
		if (!lammps_section_line(reader, &line)) {
			MD_LOG_ERROR("LAMMPS: The %s section ends after %zu of its %zu entries", section, i, n);
			return false;
		}
		const size_t num_tok = extract_tokens(tok, 2 + num_atom_ids, &line);
		if (num_tok < 2 + num_atom_ids) {
			MD_LOG_ERROR("LAMMPS: Failed to parse a line of the %s section, expected %zu tokens, got %zu", section, 2 + num_atom_ids, num_tok);
			return false;
		}
		int32_t* o = out + i * stride;
		for (size_t j = 0; j < 2 + num_atom_ids; ++j) {
			o[j] = (int32_t)parse_int(tok[j]);
		}
	}
	return true;
}

static bool md_lammps_data_parse(md_lammps_data_t* data, md_reader_i* reader, const char* format, struct md_allocator_i* alloc) {
	MD_ASSERT(data);
	MD_ASSERT(reader);
	MD_ASSERT(alloc);
	str_t line;

	int mappings[ATOM_FIELD_COUNT];
	const size_t num_fields = interpret_format(mappings, format, 0, 0);
	if (num_fields == 0) {
		return false;
	}
	const size_t num_columns = lammps_atom_columns(mappings);

	// Read the title of the file
	if (!md_reader_extract_line(reader, &line)) {
		MD_LOG_ERROR("Failed to parse lammps title");
		return false;
	}

	MD_MEMSET(data, 0, sizeof(md_lammps_data_t));

	str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

	md_temp_scope_t temp = md_temp_begin_avoid(alloc);
	md_hashmap32_t labels = { .allocator = md_temp_allocator(temp) };
	bool result = false;
	bool header = true;
	bool has_atoms = false;

	const char* next;
	while (md_reader_find_line(reader, &line, &next)) {
		reader->cur = next;
		str_t comment;
		const str_t content = lammps_content(line, &comment);
		if (content.len == 0) continue;

		const int sec = lammps_section_kind(content);
		if (sec == SEC_NONE) {
			// A header line, or a line of a section which is not read
			if (header && !lammps_header_line(data, content)) goto done;
			continue;
		}
		header = false;
		lammps_types_init(data, alloc);

		switch (sec) {
		case SEC_ATOM_TYPE_LABELS:
			for (size_t i = 0; i < data->num_particle_types; ++i) {
				str_t tok[2];
				if (!lammps_section_line(reader, &line) || extract_tokens(tok, 2, &line) != 2) {
					MD_LOG_ERROR("LAMMPS: The Atom Type Labels section ends after %zu of %zu atom types", i, data->num_particle_types);
					goto done;
				}
				const int64_t type = parse_int(tok[0]);
				if (!is_int(tok[0]) || type < 1 || (size_t)type > data->num_particle_types) {
					MD_LOG_ERROR("LAMMPS: Invalid atom type '" STR_FMT "' in Atom Type Labels", STR_ARG(tok[0]));
					goto done;
				}
				str_copy_to_char_buf(data->particle_types[type - 1].label, sizeof(data->particle_types[type - 1].label), tok[1]);
				md_hashmap_add(&labels, md_hash64(tok[1].ptr, tok[1].len, 0), (uint32_t)type);
			}
			break;
		case SEC_MASSES:
			if (!data->num_particle_types) {
				MD_LOG_ERROR("Encountered Mass entries, but number of atom types was not set or zero");
				goto done;
			}
			for (size_t i = 0; i < data->num_particle_types; ++i) {
				str_t tok[2];
				if (!lammps_section_line(reader, &line) || extract_tokens(tok, 2, &line) != 2) {
					MD_LOG_ERROR("Number of entries in Masses section did not match the number of atom types");
					goto done;
				}
				const int32_t type = lammps_atom_type(tok[0], &labels);
				if (type < 1 || (size_t)type > data->num_particle_types) {
					MD_LOG_ERROR("Invalid atom type index in Masses: '" STR_FMT "'", STR_ARG(tok[0]));
					goto done;
				}
				md_lammps_atom_type_t* t = &data->particle_types[type - 1];
				t->mass = (float)parse_float(tok[1]);
				const md_atomic_number_t z = md_atomic_number_infer_from_mass(t->mass);
				if (z != MD_Z_X) {
					t->radius = md_atomic_number_vdw_radius(z);
				}
			}
			break;
		case SEC_PAIR_COEFFS:
			if (!data->num_particle_types) {
				MD_LOG_ERROR("Encountered Pair Coeff entries, but number of atom types was not set or zero");
				goto done;
			}
			if (!parse_pair_coeffs(data->particle_types, data->num_particle_types, reader, comment, &labels)) {
				goto done;
			}
			break;
		case SEC_ATOMS: {
			if (!data->num_atoms) {
				MD_LOG_ERROR("Encountered Atom entries, but number of atoms was not set or zero");
				goto done;
			}
			md_array_resize(data->atoms, data->num_atoms, alloc);
			bool sorted = true;
			for (size_t i = 0; i < data->num_atoms; ++i) {
				str_t tok[24];
				if (!lammps_section_line(reader, &line)) {
					MD_LOG_ERROR("LAMMPS: The Atoms section ends after %zu of %zu atoms", i, data->num_atoms);
					goto done;
				}
				const str_t atom_line = line;
				if (extract_tokens(tok, num_columns, &line) < num_columns) {
					MD_LOG_ERROR("LAMMPS: Failed to parse atom line '" STR_FMT "', expected at least %zu columns (%s)", STR_ARG(atom_line), num_columns, format);
					goto done;
				}
				md_lammps_atom_t* atom = &data->atoms[i];
				atom->id     = (int32_t)parse_int(tok[mappings[ATOM_FIELD_ID]]);
				atom->resid  = mappings[ATOM_FIELD_RESID] != -1 ? (int32_t)parse_int(tok[mappings[ATOM_FIELD_RESID]]) : -1;
				atom->type   = lammps_atom_type(tok[mappings[ATOM_FIELD_TYPE]], &labels);
				atom->charge = mappings[ATOM_FIELD_Q] != -1 ? (float)parse_float(tok[mappings[ATOM_FIELD_Q]]) : 0.0f;
				atom->x      = (float)parse_float(tok[mappings[ATOM_FIELD_X]]);
				atom->y      = (float)parse_float(tok[mappings[ATOM_FIELD_Y]]);
				atom->z      = (float)parse_float(tok[mappings[ATOM_FIELD_Z]]);
				if (atom->type < 1 || (data->num_particle_types && (size_t)atom->type > data->num_particle_types)) {
					MD_LOG_ERROR("LAMMPS: Invalid atom type '" STR_FMT "' of atom %d", STR_ARG(tok[mappings[ATOM_FIELD_TYPE]]), atom->id);
					goto done;
				}
				sorted &= i == 0 || data->atoms[i - 1].id < atom->id;
			}
			if (!sorted) {
				qsort(data->atoms, data->num_atoms, sizeof(md_lammps_atom_t), compare_atom);
			}
			for (size_t i = 1; i < data->num_atoms; ++i) {
				if (data->atoms[i].id == data->atoms[i - 1].id) {
					MD_LOG_ERROR("LAMMPS: The atom id %d is given twice", data->atoms[i].id);
					goto done;
				}
			}
			has_atoms = true;
		} break;
		case SEC_BONDS:
			md_array_resize(data->bonds, data->num_bonds, alloc);
			if (!lammps_parse_topology((int32_t*)data->bonds, 4, 2, data->num_bonds, reader, "Bonds")) goto done;
			break;
		case SEC_ANGLES:
			md_array_resize(data->angles, data->num_angles, alloc);
			if (!lammps_parse_topology((int32_t*)data->angles, 5, 3, data->num_angles, reader, "Angles")) goto done;
			break;
		case SEC_DIHEDRALS:
			md_array_resize(data->dihedrals, data->num_dihedrals, alloc);
			if (!lammps_parse_topology((int32_t*)data->dihedrals, 6, 4, data->num_dihedrals, reader, "Dihedrals")) goto done;
			break;
		case SEC_IMPROPERS:
			md_array_resize(data->impropers, data->num_impropers, alloc);
			if (!lammps_parse_topology((int32_t*)data->impropers, 6, 4, data->num_impropers, reader, "Impropers")) goto done;
			break;
		default:
			// Its lines are skipped as those of no section: none of them is a keyword line
			break;
		}
	}

	if (!has_atoms) {
		MD_LOG_ERROR("LAMMPS: The data file has no Atoms section");
		goto done;
	}
	result = true;
done:
	md_temp_end(temp);
	return result;
}

bool md_lammps_validate_atom_format(char* err_buf, size_t err_cap, const char* format) {
	int mappings[MD_COUNTOF(atom_field_name)] = {0};
	return interpret_format(mappings, format, err_buf, err_cap) > 0;
}

bool md_lammps_data_parse_str(md_lammps_data_t* data, str_t str, const char* format, struct md_allocator_i* alloc) {
	MD_ASSERT(data);
	MD_ASSERT(alloc);

	md_reader_i line_reader = md_reader_from_str(str);
	return md_lammps_data_parse(data, &line_reader, format, alloc);
}

bool md_lammps_data_parse_file(md_lammps_data_t* data, str_t filename, const char* format, struct md_allocator_i* alloc) {
	bool result = false;
	md_file_t file = {0};
	if (md_file_open(&file, filename, MD_FILE_READ)) {
		const size_t cap = MD_MEGABYTES(1);
		md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);

		md_file_reader_t line_reader;
		md_file_reader_init(&line_reader, file, cap, md_temp_allocator(temp_scope));
		result = md_lammps_data_parse(data, &line_reader.r, format, alloc);

		md_file_reader_close(&line_reader);
		md_temp_end(temp_scope);
		md_file_close(&file);
	}
	else {
		MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
	}
	return result;
}

void md_lammps_data_free(md_lammps_data_t* data, struct md_allocator_i* alloc) {
	MD_ASSERT(data);
	if (data->atoms) md_array_free(data->atoms, alloc);
	if (data->bonds) md_array_free(data->bonds, alloc);
	if (data->angles) md_array_free(data->angles, alloc);
	if (data->dihedrals) md_array_free(data->dihedrals, alloc);
	if (data->impropers) md_array_free(data->impropers, alloc);
	if (data->particle_types) md_array_free(data->particle_types, alloc);
	MD_MEMSET(data, 0, sizeof(md_lammps_data_t));
}


// Chemical formula in Hill order (C, H, then the other elements alphabetically), empty when an atom has no element
static str_t lammps_formula(char* buf, size_t cap, const md_system_t* sys, md_urange_t range) {
	uint32_t count[256] = {0};
	for (uint32_t i = range.beg; i < range.end; ++i) {
		const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
		if (z == 0) return (str_t){0};
		count[z] += 1;
	}
	int len = 0;
	const md_atomic_number_t first[2] = { MD_Z_C, MD_Z_H };
	for (int k = 0; k < 2 && count[MD_Z_C] > 0; ++k) {
		const md_atomic_number_t z = first[k];
		if (count[z] == 0) continue;
		const str_t sym = md_atomic_number_symbol(z);
		len += md_snprintf(buf + len, cap - len, STR_FMT, STR_ARG(sym));
		if (count[z] > 1 && len < (int)cap) len += md_snprintf(buf + len, cap - len, "%u", count[z]);
		if (len >= (int)cap) return (str_t){0};
		count[z] = 0;
	}
	for (;;) {
		// The next element alphabetically
		int best = -1;
		for (int z = 1; z < 256; ++z) {
			if (count[z] == 0) continue;
			if (best == -1 || str_cmp_lex(md_atomic_number_symbol((md_atomic_number_t)z), md_atomic_number_symbol((md_atomic_number_t)best)) < 0) best = z;
		}
		if (best == -1) break;
		const str_t sym = md_atomic_number_symbol((md_atomic_number_t)best);
		len += md_snprintf(buf + len, cap - len, STR_FMT, STR_ARG(sym));
		if (count[best] > 1 && len < (int)cap) len += md_snprintf(buf + len, cap - len, "%u", count[best]);
		if (len >= (int)cap) return (str_t){0};
		count[best] = 0;
	}
	return (str_t){ buf, (size_t)len };
}

// The molecules of a LAMMPS data file are its residues. Each is an instance, of a molecule type per distinct sequence of
// particle types, which the file implies without naming: the type is named by its formula. Each instance is one
// molecule whatever its bonds, as the file has it. What its parts are is inferred (md_util_system_infer_parts).
static void lammps_molecules_as_instances(md_system_t* sys) {
	md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
	md_hashmap32_t type_map = { .allocator = md_temp_allocator(temp) };

	for (size_t c = 0; c < sys->residue.count; ++c) {
		const md_urange_t range = md_system_residue_atom_range(sys, c);
		const uint64_t key = md_hash64(sys->atom.type_idx + range.beg, (range.end - range.beg) * sizeof(md_particle_type_idx_t), range.end - range.beg);

		md_molecule_type_idx_t type_idx = -1;
		const uint32_t* found = md_hashmap_get(&type_map, key);
		if (found) {
			type_idx = (md_molecule_type_idx_t)*found;
		} else {
			char buf[64];
			str_t name = lammps_formula(buf, sizeof(buf), sys, range);
			if (str_empty(name)) {
				const int len = md_snprintf(buf, sizeof(buf), "molecule type %zu", sys->molecule.type.count + 1);
				name = (str_t){ buf, (size_t)len };
			}
			type_idx = md_molecule_type_data_push(&sys->molecule.type, name, MD_MOLECULE_TYPE_FLAG_NONE, sys->alloc);
			md_hashmap_add(&type_map, key, (uint32_t)type_idx);
		}

		md_array_push(sys->instance.id, (md_label_t){0}, sys->alloc);
		md_array_push(sys->instance.auth_id, (md_label_t){0}, sys->alloc);
		md_array_push(sys->instance.res_offset, (uint32_t)c, sys->alloc);
		md_array_push(sys->instance.part_type_idx, -1, sys->alloc);
		md_array_push(sys->instance.mol_type_idx, type_idx, sys->alloc);
		sys->instance.count += 1;
	}
	if (sys->instance.count) {
		md_array_push(sys->instance.res_offset, (uint32_t)sys->residue.count, sys->alloc);
	}
	md_temp_end(temp);
}

bool md_lammps_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_lammps_data_t* data) {
	MD_ASSERT(sys);
	MD_ASSERT(state);
	MD_ASSERT(data);

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

	const size_t capacity = MD_ROUND_UP(data->num_atoms, 16);

	md_array_ensure(state->xyz,		capacity, state->alloc);
	md_array_ensure(sys->atom.type_idx, capacity, sys->alloc);
	md_array_ensure(sys->atom.flags,    capacity, sys->alloc);

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);

	bool has_resid = (data->num_atoms > 0 && data->atoms[0].resid != -1);

	// Reset atom data and initialize to default value
	sys->atom.type.count = 0;
	md_particle_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0, 0, 0, 0, sys->alloc);

	// The types: type i + 1 at i, or none where the file has no Masses section (nor labels) to make them of
	const size_t num_types = data->particle_types ? data->num_particle_types : 0;
	float* type_masses			= md_temp_alloc_array(temp_scope, float, MD_MAX(num_types, 1));
	md_atomic_number_t* type_z	= md_temp_alloc_array(temp_scope, md_atomic_number_t, MD_MAX(num_types, 1));
	md_particle_type_idx_t* type_map = md_temp_alloc_array(temp_scope, md_particle_type_idx_t, MD_MAX(num_types, 1));
	for (size_t i = 0; i < num_types; ++i) {
		type_masses[i] = data->particle_types[i].mass;
	}
	
	// The element of a type is told by its mass, else (no mass given) by its label
	md_util_element_from_mass(type_z, type_masses, num_types);
	for (size_t i = 0; i < num_types; ++i) {
		const str_t label = str_from_cstr(data->particle_types[i].label);
		if (type_z[i] == MD_Z_X && type_masses[i] == 0.0f && !str_empty(label)) {
			type_z[i] = md_atomic_number_infer_from_label(label, (str_t){0}, 0);
		}
	}

	for (size_t i = 0; i < num_types; ++i) {
		// Named by its label, else by its number; each its own type, though two names may read alike as cut
		char buf[32];
		const str_t label = str_from_cstr(data->particle_types[i].label);
		const int len = str_empty(label) ? md_snprintf(buf, sizeof(buf), "type_%i", data->particle_types[i].id) : md_snprintf(buf, sizeof(buf), STR_FMT, STR_ARG(label));
		str_t type_id = {buf, (size_t)len};
		float mass   = data->particle_types[i].mass > 0 ? data->particle_types[i].mass : md_atomic_number_mass(type_z[i]);
		float radius = data->particle_types[i].radius > 0 ? data->particle_types[i].radius : md_atomic_number_vdw_radius(type_z[i]);
		uint32_t color = type_z[i] != MD_Z_X ? md_atomic_number_cpk_color(type_z[i]) : 0xFF808080; // Gray for unknown elements
		
		type_map[i] = md_particle_type_add(&sys->atom.type, type_id, type_id, type_z[i], mass, radius, color, 0, sys->alloc);
	}

	// A general triclinic box turns the coordinates into the orientation md_unitcell_t holds
	double R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
	double origin[3] = { data->cell.xlo, data->cell.ylo, data->cell.zlo };
	bool rotate = false;
	if (data->cell.general) {
		const double basis[3][3] = {
			{ data->cell.avec[0], data->cell.avec[1], data->cell.avec[2] },
			{ data->cell.bvec[0], data->cell.bvec[1], data->cell.bvec[2] },
			{ data->cell.cvec[0], data->cell.cvec[1], data->cell.cvec[2] },
		};
		rotate = md_unitcell_from_basis_vectors(&state->unitcell, R, basis);
		origin[0] = data->cell.origin[0];
		origin[1] = data->cell.origin[1];
		origin[2] = data->cell.origin[2];
	} else {
		state->unitcell = md_unitcell_from_basis_parameters(data->cell.xhi - data->cell.xlo, data->cell.yhi - data->cell.ylo, data->cell.zhi - data->cell.zlo, data->cell.xy, data->cell.xz, data->cell.yz);
	}

	int prev_resid = -1;
	size_t num_bad_types = 0;
	for (size_t i = 0; i < data->num_atoms; ++i) {
		const int32_t type = data->atoms[i].type;
		md_particle_type_idx_t type_idx = 0;
		if (1 <= type && (size_t)type <= num_types) {
			type_idx = type_map[type - 1];
		} else if (num_types) {
			num_bad_types += 1;
		}

		int resid = data->atoms[i].resid;
		if (has_resid) {
			if (resid != prev_resid || i == 0) {
				md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, sys->alloc);
				md_array_push(sys->residue.seq_id, resid, sys->alloc);
				// No name in the file: the generic one (a number would not fit a residue name past 99, and a formula
				// may read as a residue of its own, "C" a cytidine)
				md_array_push(sys->residue.name, make_label(STR_LIT("MOL")), sys->alloc);
				md_array_push(sys->residue.flags, 0, sys->alloc);
				sys->residue.count += 1;
			}
		}

		const double p[3] = { data->atoms[i].x - origin[0], data->atoms[i].y - origin[1], data->atoms[i].z - origin[2] };
		const vec3_t xyz = rotate ?
			vec3_set((float)(R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2]), (float)(R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2]), (float)(R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2])) :
			vec3_set((float)p[0], (float)p[1], (float)p[2]);

		md_array_push_no_grow(sys->atom.type_idx, type_idx);
		md_array_push_no_grow(state->xyz, xyz);
		md_array_push_no_grow(sys->atom.flags, 0);
		sys->atom.count +=1;

		prev_resid = resid;
	}
	if (num_bad_types) {
		MD_LOG_ERROR("LAMMPS: %zu atoms are of no atom type of the file", num_bad_types);
	}

	if (has_resid) {
		md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, sys->alloc); // Final sentinel
		lammps_molecules_as_instances(sys);
	}

	// The per atom charge the 'charge' and 'full' atom styles carry. It was parsed and then dropped;
	// nothing in mdlib reads a partial charge, so it belongs in the attribute table rather than in
	// md_atom_data_t. A style without a q column leaves every value at 0, which is uniform and
	// therefore publishes nothing - so this needs no guard on the atom style.
	{
		md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
		float* charge = (float*)md_temp_alloc(temp, sizeof(float) * MD_MAX(data->num_atoms, (size_t)1));
		if (charge) {
			for (size_t i = 0; i < data->num_atoms; ++i) {
				charge[i] = data->atoms[i].charge;
			}
			md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/charge"), md_unit_none(), 1, charge, data->num_atoms);
		}
		md_temp_end(temp);
	}

	// Create bonds: atom ids are the atoms' in order (sorted), 1 to N but where atoms have been taken away
	if (data->num_bonds > 0 && data->bonds) {
		md_array_ensure(sys->bond.pairs, data->num_bonds, sys->alloc);
		md_array_ensure(sys->bond.flags, data->num_bonds, sys->alloc);

		const size_t N = data->num_atoms;
		const bool dense = N > 0 && data->atoms[0].id == 1 && data->atoms[N - 1].id == (int32_t)N;
		size_t num_bad = 0;
		for (size_t i = 0; i < data->num_bonds; ++i) {
			uint32_t idx[2];
			for (int k = 0; k < 2; ++k) {
				const int32_t id = data->bonds[i].atom_id[k];
				idx[k] = UINT32_MAX;
				if (dense) {
					if (1 <= id && (size_t)id <= N) idx[k] = (uint32_t)(id - 1);
				} else {
					size_t lo = 0, hi = N;
					while (lo < hi) {
						const size_t mid = (lo + hi) / 2;
						if (data->atoms[mid].id < id) lo = mid + 1; else hi = mid;
					}
					if (lo < N && data->atoms[lo].id == id) idx[k] = (uint32_t)lo;
				}
			}
			if (idx[0] == UINT32_MAX || idx[1] == UINT32_MAX || idx[0] == idx[1]) {
				num_bad += 1;
				continue;
			}
			md_atom_pair_t pair = {{ MD_MIN(idx[0], idx[1]), MD_MAX(idx[0], idx[1]) }};
			md_bond_flags_t flag = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
			md_array_push_no_grow(sys->bond.pairs, pair);
			md_array_push_no_grow(sys->bond.flags, flag);
			sys->bond.count += 1;
		}
		if (num_bad) {
			MD_LOG_ERROR("LAMMPS: %zu bonds are of atoms the file does not have; they are left out", num_bad);
		}
		md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
	}

	// The atoms have no names to classify the molecules by, but water and monatomic ions are told by their atoms
	md_util_system_infer_residue_flags(sys);

	md_temp_end(temp_scope);

	MD_ASSERT(md_array_size(state->xyz) == sys->atom.count);
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
//
// A dump file (dump atom or custom): frames of "ITEM: TIMESTEP", "ITEM: NUMBER OF ATOMS", "ITEM: BOX BOUNDS" and
// "ITEM: ATOMS <columns>", each followed by its lines, where "ITEM: UNITS" (the first frame) and "ITEM: TIME" (with
// dump_modify time yes) may come first. https://docs.lammps.org/dump.html

typedef struct {
	double xlo, xhi, xy;
	double ylo, yhi, xz;
	double zlo, zhi, yz;
	// A general triclinic box ("abc origin"): its vectors and origin instead
	bool general;
	double abc[3][3];
	double origin[3];
} box_bounds_t;

typedef struct {
	int64_t timestep;
	size_t num_atoms;
	box_bounds_t box_bounds;
} header_t;

// Each of the three boundary flags of a dimension: p, f, s or m, twice ("pp", "fm", ...)
static bool is_boundary(str_t tok) {
	if (tok.len != 2) return false;
	for (int i = 0; i < 2; ++i) {
		const char c = tok.ptr[i];
		if (c != 'p' && c != 'f' && c != 's' && c != 'm') return false;
	}
	return true;
}

static bool parse_box_bounds(box_bounds_t* box_bounds, md_reader_i* reader) {
	MD_ASSERT(reader);
	str_t line;

	if (!md_reader_extract_line(reader, &line)) {
		MD_LOG_ERROR("Could not extract box bounds");
		return false;
	}

	if (!str_eq_cstr_n(line, "ITEM: BOX BOUNDS", 16)) {
		MD_LOG_ERROR("Unexpected beginning of line: '" STR_FMT "' expected ITEM: BOX BOUNDS", STR_ARG(line));
		return false;
	}

	// "pp pp pp", "xy xz yz pp pp pp" or "abc origin pp pp pp", of any boundaries
	str_t rest = str_substr(line, 16, SIZE_MAX);
	str_t tok[8];
	const size_t num_tok = extract_tokens(tok, MD_COUNTOF(tok), &rest);
	size_t skip = 0;
	bool tilt = false, general = false;
	if (num_tok >= 3 && str_eq(tok[0], STR_LIT("xy")) && str_eq(tok[1], STR_LIT("xz")) && str_eq(tok[2], STR_LIT("yz"))) {
		tilt = true;
		skip = 3;
	} else if (num_tok >= 2 && str_eq(tok[0], STR_LIT("abc")) && str_eq(tok[1], STR_LIT("origin"))) {
		general = true;
		skip = 2;
	}
	if (num_tok != skip + 3 || !is_boundary(tok[skip]) || !is_boundary(tok[skip + 1]) || !is_boundary(tok[skip + 2])) {
		MD_LOG_ERROR("Unrecognized format in ITEM: BOX BOUNDS: '" STR_FMT "'", STR_ARG(line));
		return false;
	}

	const size_t num_values = general ? 4 : tilt ? 3 : 2;
	double values[3][4] = {0};
	for (size_t i = 0; i < 3; ++i) {
		if (!md_reader_extract_line(reader, &line)) {
			MD_LOG_ERROR("Failed to extract box bounds");
			return false;
		}
		const size_t n = extract_tokens(tok, num_values, &line);
		if (n < num_values) {
			MD_LOG_ERROR("Failed to extract box bounds");
			return false;
		}
		for (size_t j = 0; j < num_values; ++j) {
			values[i][j] = parse_float(tok[j]);
		}
	}

	if (box_bounds) {
		MD_MEMSET(box_bounds, 0, sizeof(*box_bounds));
		if (general) {
			box_bounds->general = true;
			for (size_t i = 0; i < 3; ++i) {
				box_bounds->abc[i][0] = values[i][0];
				box_bounds->abc[i][1] = values[i][1];
				box_bounds->abc[i][2] = values[i][2];
				box_bounds->origin[i] = values[i][3];
			}
		} else {
			box_bounds->xlo = values[0][0]; box_bounds->xhi = values[0][1]; box_bounds->xy = values[0][2];
			box_bounds->ylo = values[1][0]; box_bounds->yhi = values[1][1]; box_bounds->xz = values[1][2];
			box_bounds->zlo = values[2][0]; box_bounds->zhi = values[2][1]; box_bounds->yz = values[2][2];
		}
	}

	return true;
}

// The header of a frame, from its first line on. UNITS and TIME before the TIMESTEP are skipped.
static bool parse_header(header_t* header, md_reader_i* reader) {
	MD_ASSERT(reader);
	str_t line;

	for (;;) {
		if (!md_reader_peek_line(reader, &line)) return false;
		if (str_eq_cstr_n(line, "ITEM: UNITS", 11) || (str_eq_cstr_n(line, "ITEM: TIME", 10) && !str_eq_cstr_n(line, "ITEM: TIMESTEP", 14))) {
			md_reader_skip_line(reader);
			md_reader_skip_line(reader);
			continue;
		}
		break;
	}

	if (!str_eq_cstr_n(line, "ITEM: TIMESTEP", 14)) {
		return false;
	}

	md_reader_skip_line(reader);
	if (md_reader_extract_line(reader, &line)) {
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

	if (md_reader_peek_line(reader, &line) &&
		str_eq_cstr_n(line, "ITEM: NUMBER OF ATOMS", 21) &&
		md_reader_skip_line(reader) &&
		md_reader_extract_line(reader, &line))
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

	if (md_reader_peek_line(reader, &line) &&
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

// The columns of the ITEM: ATOMS line str: of the id, of the coordinates (x y z, xs ys zs, xu yu zu or xsu ysu zsu, the
// first kind the line has) and of the image flags, which unwrap wrapped coordinates
static bool parse_coord_mappings(coord_mappings_t* mappings, str_t str) {
	MD_ASSERT(mappings);

	str_t tokens[130];
	const size_t num_tokens = extract_tokens(tokens, MD_COUNTOF(tokens), &str);
	if (num_tokens < 2 || !str_eq(tokens[0], STR_LIT("ITEM:")) || !str_eq(tokens[1], STR_LIT("ATOMS"))) {
		return false;
	}

	const str_t* labels = tokens + 2;
	const size_t num_labels = num_tokens - 2;
	if (num_labels > 127) {
		MD_LOG_ERROR("LAMMPS: More columns than 127 in ITEM: ATOMS");
		return false;
	}

	MD_MEMSET(mappings, -1, sizeof(coord_mappings_t));
	mappings->flags = COORD_FLAG_NONE;

	static const struct { const char* name[3]; int flags; bool wrapped; } kinds[] = {
		{ { "x",   "y",   "z"   }, COORD_FLAG_CARTESIAN, true  },
		{ { "xs",  "ys",  "zs"  }, COORD_FLAG_SCALED,    true  },
		{ { "xu",  "yu",  "zu"  }, COORD_FLAG_CARTESIAN, false },
		{ { "xsu", "ysu", "zsu" }, COORD_FLAG_SCALED,    false },
	};
	int first_col = INT32_MAX;
	bool wrapped = false;
	int image[3] = { -1, -1, -1 };
	for (size_t k = 0; k < MD_COUNTOF(kinds); ++k) {
		int col[3] = { -1, -1, -1 };
		for (size_t i = 0; i < num_labels; ++i) {
			for (int d = 0; d < 3; ++d) {
				if (col[d] == -1 && str_eq_cstr(labels[i], kinds[k].name[d])) col[d] = (int)i;
			}
		}
		if (col[0] >= 0 && col[1] >= 0 && col[2] >= 0 && col[0] < first_col) {
			first_col = col[0];
			mappings->coord_idx[0] = (int8_t)col[0];
			mappings->coord_idx[1] = (int8_t)col[1];
			mappings->coord_idx[2] = (int8_t)col[2];
			mappings->flags = (int8_t)kinds[k].flags;
			wrapped = kinds[k].wrapped;
		}
	}
	for (size_t i = 0; i < num_labels; ++i) {
		if (str_eq(labels[i], STR_LIT("id")) && mappings->id_idx == -1) mappings->id_idx = (int8_t)i;
		if (str_eq(labels[i], STR_LIT("ix")) && image[0] == -1) image[0] = (int)i;
		if (str_eq(labels[i], STR_LIT("iy")) && image[1] == -1) image[1] = (int)i;
		if (str_eq(labels[i], STR_LIT("iz")) && image[2] == -1) image[2] = (int)i;
	}
	if (wrapped && image[0] >= 0 && image[1] >= 0 && image[2] >= 0) {
		mappings->image_idx[0] = (int8_t)image[0];
		mappings->image_idx[1] = (int8_t)image[1];
		mappings->image_idx[2] = (int8_t)image[2];
		mappings->flags |= COORD_FLAG_UNWRAP;
	}

	mappings->num_coord_tokens = (int8_t)num_labels;

	return mappings->id_idx != -1 && mappings->flags != COORD_FLAG_NONE;
}

// The columns of an atom line that are read: up to the last of those the mappings name, plus one
static size_t coord_mappings_columns(const coord_mappings_t* m) {
	int max = m->id_idx;
	for (int d = 0; d < 3; ++d) {
		max = MD_MAX(max, m->coord_idx[d]);
		if (m->flags & COORD_FLAG_UNWRAP) max = MD_MAX(max, m->image_idx[d]);
	}
	return (size_t)(max + 1);
}

typedef struct {
	int32_t id;
	float x, y, z;
} id_xyz_t;

static int compare_id_xyz(const void* a, const void* b) {
	const id_xyz_t* id_xyz_a = (const id_xyz_t*)a;
	const id_xyz_t* id_xyz_b = (const id_xyz_t*)b;
	return (id_xyz_a->id > id_xyz_b->id) - (id_xyz_a->id < id_xyz_b->id);
}

// The cell a frame's box bounds describe, the corner it starts at, and the rotation of the coordinates into the
// orientation md_unitcell_t holds (the identity but for a general triclinic box). R may be NULL.
// https://docs.lammps.org/Howto_triclinic.html: the bounds of a triclinic box are those of the box with its tilts,
// xlo_bound = xlo + MD_MIN(0, xy, xz, xy + xz), xhi_bound = xhi + MD_MAX(0, xy, xz, xy + xz), ylo_bound = ylo + MD_MIN(0, yz),
// yhi_bound = yhi + MD_MAX(0, yz).
static md_unitcell_t lammps_cell_from_bounds(const box_bounds_t* bb, double lo[3], double R[3][3]) {
	if (bb->general) {
		md_unitcell_t cell;
		md_unitcell_from_basis_vectors(&cell, R, bb->abc);
		if (lo) {
			lo[0] = bb->origin[0];
			lo[1] = bb->origin[1];
			lo[2] = bb->origin[2];
		}
		return cell;
	}
	const double xlo = bb->xlo - MD_MIN(0.0, MD_MIN(bb->xy, MD_MIN(bb->xz, bb->xy + bb->xz)));
	const double xhi = bb->xhi - MD_MAX(0.0, MD_MAX(bb->xy, MD_MAX(bb->xz, bb->xy + bb->xz)));
	const double ylo = bb->ylo - MD_MIN(0.0, bb->yz);
	const double yhi = bb->yhi - MD_MAX(0.0, bb->yz);
	const double zlo = bb->zlo;
	const double zhi = bb->zhi;
	if (lo) {
		lo[0] = xlo;
		lo[1] = ylo;
		lo[2] = zlo;
	}
	if (R) {
		MD_MEMSET(R, 0, sizeof(double) * 9);
		R[0][0] = R[1][1] = R[2][2] = 1.0;
	}
	return md_unitcell_from_basis_parameters(xhi - xlo, yhi - ylo, zhi - zlo, bb->xy, bb->xz, bb->yz);
}

// One frame's text, from its first ITEM on: the cell, and the coordinates in id order written to out_x[i * stride] and
// so on - planar with stride 1, packed with the three pointers one apart and stride 3. Relative to the corner of the box
// (as those of the data file are), turned as its cell is (a general triclinic box), and unwrapped by the image flags
// where the coordinates are wrapped and the frame has them.
static bool lammps_decode_frame_text(const coord_mappings_t* mappings, str_t str, size_t* out_num_atoms, md_unitcell_t* out_cell, float* out_x, float* out_y, float* out_z, size_t stride) {
	bool result = false;
	str_t tokens[128];
	const bool output_coords = out_x != NULL && out_y != NULL && out_z != NULL;

	md_reader_i reader = md_reader_from_str(str);
	str_t line;

	header_t header;
	if (!parse_header(&header, &reader)) {
		MD_LOG_ERROR("Could not parse header");
		return false;
	}

	double lo[3];
	double R[3][3];
	const md_unitcell_t cell = lammps_cell_from_bounds(&header.box_bounds, lo, R);
	double A[3][3];
	md_unitcell_A_extract_double(A, &cell);

	if (output_coords) {
		if (!md_reader_extract_line(&reader, &line) || !str_eq_cstr_n(line, "ITEM: ATOMS", 11)) {
			MD_LOG_ERROR("Expected ITEM: ATOMS after header");
			goto done;
		}
		const size_t num_columns = coord_mappings_columns(mappings);
		if (num_columns > MD_COUNTOF(tokens)) goto done;
		size_t line_count = 0;
		bool coords_result = false;

		// We need to store the coordinates in a temporary buffer since we need to sort them by id
		md_temp_scope_t temp_scope = md_temp_begin();
		id_xyz_t* id_xyz = md_temp_alloc_array(temp_scope, id_xyz_t, MD_MAX(header.num_atoms, 1));
		bool sorted = true;

		while (line_count < header.num_atoms && md_reader_extract_line(&reader, &line)) {
			const size_t num_tokens = extract_tokens(tokens, num_columns, &line);
			if (num_tokens < num_columns) {
				MD_LOG_ERROR("Unexpected number of tokens in ITEM: ATOMS line, got %zu, expected at least %zu", num_tokens, num_columns);
				goto coords_done;
			}

			const int32_t id = (int32_t)parse_int(tokens[mappings->id_idx]);
			const double c[3] = {
				parse_float(tokens[mappings->coord_idx[0]]),
				parse_float(tokens[mappings->coord_idx[1]]),
				parse_float(tokens[mappings->coord_idx[2]]),
			};
			double p[3];
			if (mappings->flags & COORD_FLAG_SCALED) {
				// Fractions of the cell vectors (the columns of A)
				for (int k = 0; k < 3; ++k) p[k] = A[0][k] * c[0] + A[1][k] * c[1] + A[2][k] * c[2];
			} else {
				const double d[3] = { c[0] - lo[0], c[1] - lo[1], c[2] - lo[2] };
				for (int k = 0; k < 3; ++k) p[k] = R[k][0] * d[0] + R[k][1] * d[1] + R[k][2] * d[2];
			}

			if (mappings->flags & COORD_FLAG_UNWRAP) {
				const double img[3] = {
					(double)parse_int(tokens[mappings->image_idx[0]]),
					(double)parse_int(tokens[mappings->image_idx[1]]),
					(double)parse_int(tokens[mappings->image_idx[2]]),
				};
				for (int k = 0; k < 3; ++k) p[k] += A[0][k] * img[0] + A[1][k] * img[1] + A[2][k] * img[2];
			}

			id_xyz[line_count] = (id_xyz_t){ id, (float)p[0], (float)p[1], (float)p[2] };
			sorted &= line_count == 0 || id_xyz[line_count - 1].id < id;
			line_count += 1;
		}

		// Sort atoms by id
		if (!sorted) {
			qsort(id_xyz, line_count, sizeof(id_xyz_t), compare_id_xyz);
		}

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

// The frames of a dump: where each begins, its step and cell, and the columns of its atom lines (those of the first
// frame, which all must have). The atom lines are skipped, not read: they are read when a frame is.
static bool lammps_trajectory_parse(lammps_cache_t* cache, md_reader_i* reader, struct md_allocator_i* alloc) {
	MD_ASSERT(reader);
	MD_ASSERT(alloc);
	str_t line;
	size_t num_atoms = 0;
	size_t num_frames = 0;
	coord_mappings_t mappings = {0};
	int64_t frames_end = 0;

	for (;;) {
		const int64_t offset = md_reader_tell(reader);
		header_t header;
		if (!parse_header(&header, reader)) {
			break;
		}

		if (num_frames == 0) {
			num_atoms = header.num_atoms;
		} else if (num_atoms != header.num_atoms) {
			MD_LOG_ERROR("Number of atoms does not match between frames");
			return false;
		}

		if (!md_reader_peek_line(reader, &line) ||
			!str_eq_cstr_n(line, "ITEM: ATOMS", 11)) {
			MD_LOG_ERROR("Expected ITEM: ATOMS after header");
			return false;
		}

		if (num_frames == 0) {
			if (!parse_coord_mappings(&mappings, line)) {
				MD_LOG_ERROR("Could not parse coord mappings");
				return false;
			}
		}

		// Skip ITEM: ATOMS line, and the atom lines
		md_reader_skip_line(reader);
		if (md_reader_skip_lines(reader, header.num_atoms) != header.num_atoms) {
			MD_LOG_ERROR("The last frame of the dump is cut short: it is left out");
			break;
		}

		md_array_push(cache->frame_offsets, offset, alloc);
		md_array_push(cache->frame_steps, header.timestep, alloc);
		{
			const md_unitcell_t cell = lammps_cell_from_bounds(&header.box_bounds, NULL, NULL);
			float A[3][3];
			md_unitcell_A_extract_float(A, &cell);
			md_array_push_array(cache->frame_cells, &A[0][0], 9, alloc);
		}
		num_frames += 1;
		frames_end = md_reader_tell(reader);
	}

	if (num_frames == 0) {
		return false;
	}

	cache->header.num_frames = num_frames;
	cache->header.num_atoms = num_atoms;
	cache->coord_mappings = mappings;
	//We add the end of the last frame to frame_offsets, so frame_offset size = num_frames + 1
	md_array_push(cache->frame_offsets, frames_end, alloc);

	return true;
}

static bool lammps_trajectory_parse_file(lammps_cache_t* cache, str_t filename, struct md_allocator_i* alloc, struct md_allocator_i* avoid) {
	bool result = false;
	md_file_t file = {0};
	if (md_file_open(&file, filename, MD_FILE_READ)) {
		const int64_t cap = MD_MEGABYTES(1);
		// Clear of what alloc was made from, not of alloc itself: an arena on a temp allocator would
		// otherwise be rewound with the scratch.
		md_temp_scope_t temp_scope = md_temp_begin_avoid(avoid);

		md_file_reader_t line_reader;
		md_file_reader_init(&line_reader, file, cap, md_temp_allocator(temp_scope));
		result = lammps_trajectory_parse(cache, &line_reader.r, alloc);

		md_file_reader_close(&line_reader);
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
	MD_ASSERT(cache);
	MD_ASSERT(alloc);
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
	MD_MEMSET(cache, 0, sizeof(*cache));

	if (try_read_cache(cache, filename, alloc)) {
		return true;
	}
	MD_MEMSET(cache, 0, sizeof(*cache));

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
	MD_ASSERT(sys);
	if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

	md_run_source_t src;
	if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;

	const md_attribute_t* layout_attr = md_attributes_find_in(&sys->attributes, src.run, STR_LIT("source/layout"));
	const int32_t* layout = (const int32_t*)md_attribute_view(layout_attr, MD_ATTRIBUTE_TYPE_I32, 1, 1);
	if (!layout || layout_attr->format.shape[0] != LAMMPS_LAYOUT_COUNT) {
		MD_LOG_ERROR("LAMMPS: the run '" STR_FMT "' has lost its layout", STR_ARG(src.run));
		return 0;
	}
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
	char*  text = md_temp_alloc(temp, MD_MAX(frame_size, 1));
	float* xyz  = (count == N) ? (float*)dst : md_temp_alloc(temp, N * 3 * sizeof(float));
	if (text && xyz && md_attribute_io_read_at(io, src.path, src.offset[frame], text, frame_size) == frame_size) {
		size_t num_atoms = 0;
		if (lammps_decode_frame_text(&mappings, (str_t){ text, frame_size }, &num_atoms, NULL, xyz + 0, xyz + 1, xyz + 2, 3) && num_atoms == N) {
			if (count != N) {
				MD_MEMCPY(dst, xyz + first * 3, 3 * sizeof(float));
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
	MD_ASSERT(sys);
	char path_buf[4096];
	const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
	const str_t path = { path_buf, path_len };

	md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
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

