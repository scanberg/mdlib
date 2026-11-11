#include <md_xyz.h>

#include <md_system.h>
#include <md_util.h>
#include <md_unitcell.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// THE DIALECTS
//
// Each model (frame) is a line with the number of atoms, then a comment line, then a line per atom: its element
// (symbol or atomic number), x, y and z, after which a line may hold more columns, which are not read.
//
// Extended XYZ (https://github.com/libAtoms/extxyz): the comment line holds key=value pairs. Properties says what
// the columns of an atom line are, name:type:count for each, in any order ("species:S:1:pos:R:3:forces:R:3"): the
// element is species (or Z, an atomic number), the coordinates pos. Lattice is the cell, three vectors in any
// orientation. A comment with a Lattice but no Properties is the default columns, species and pos.
//
// Tinker (.xyz, .arc): the first line holds the number of atoms and the comment, and there is no comment line but
// for a line of the cell (a, b, c, alpha, beta, gamma) where the box is periodic. An atom line is its index, its name
// (an atom name of the force field: "CT", "HW", "C1"), x, y, z, its atom type and the indices of the atoms bonded to
// it. An ARC file is Tinker frames one after another.
//
// A plain XYZ comment line of six numbers is taken as a cell line as well (ARC).

#define MD_XYZ_CACHE_MAGIC      0x8265485749172bab
#define MD_XYZ_CACHE_VERSION    6
#define MD_XYZ_MOL_MAGIC        0x285ada29078a9bc8

enum {
    XYZ_TINKER          = 1,
    XYZ_ARC             = 2,    // A cell line: Tinker's periodic box, or six numbers as the comment of a plain file
    XYZ_EXTENDED        = 4,
    XYZ_STORE_COMMENT   = 8,
    XYZ_ELEMENT_Z       = 16,   // Extended: the element column is Z, atomic numbers
    XYZ_NO_ELEMENT      = 32,   // Extended: no column says what element an atom is
};

// The columns of an atom line, which the flags carry: of the element and of x (y and z follow it)
#define XYZ_COLUMNS(elem, pos)  (((uint32_t)(elem) << 8) | ((uint32_t)(pos) << 16))
#define XYZ_ELEM_COL(flags)     (((flags) >> 8) & 0xFF)
#define XYZ_POS_COL(flags)      (((flags) >> 16) & 0xFF)
#define XYZ_MAX_COL             250

// The connectivity an atom line of Tinker holds at most (Tinker's own limit is 8)
#define XYZ_MAX_CONN ((int)(sizeof(((md_xyz_coordinate_t*)0)->connectivity) / sizeof(int)))

static inline bool is_unsigned_int(md_str_t str) {
    str = md_str_trim(str);
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    while (c < end) {
        if (!md_ascii_is_digit(*c)) return false;
        ++c;
    }
    return true;
}

// The name of an atom of Tinker: a letter first ("C", "CT", "HW1", "N3+")
static inline bool is_atom_name(md_str_t str) {
    return str.len > 0 && md_ascii_is_alpha(str.ptr[0]);
}

// The value of key=value in a comment line of extended XYZ, the key in any case and at the start or after a blank.
// The value as written: in its quotes, braces or brackets where it has them, else up to the next blank.
static bool xyz_find_key(md_str_t* out_value, md_str_t line, md_str_t key) {
    const char* beg = line.ptr;
    const char* end = line.ptr + line.len;
    for (const char* c = beg; c + key.len < end; ++c) {
        if (c > beg && !md_ascii_is_whitespace(c[-1])) continue;
        if (c[key.len] != '=' || !md_str_eq_ignore_case((md_str_t){c, key.len}, key)) continue;
        const char* v = c + key.len + 1;
        const char* e = v;
        if (v < end && (*v == '"' || *v == '{' || *v == '[')) {
            const char open = *v;
            const char close = open == '"' ? '"' : open == '{' ? '}' : ']';
            int depth = 0;
            for (e = v + 1; e < end; ++e) {
                if (*e == close && depth == 0) { ++e; break; }
                if (open != '"' && *e == open) ++depth;
                if (open != '"' && *e == close) --depth;
            }
        } else {
            while (e < end && !md_ascii_is_whitespace(*e)) ++e;
        }
        *out_value = (md_str_t){ v, (size_t)(e - v) };
        return true;
    }
    return false;
}

// The columns of the atom lines of extended XYZ, from the Properties of its comment line, into flags
static bool xyz_extended_columns(uint32_t* flags, md_str_t comment) {
    md_str_t props = MD_STR_LIT("species:S:1:pos:R:3");
    md_str_t value = props;
    if (xyz_find_key(&value, comment, MD_STR_LIT("Properties"))) {
        if (value.len >= 2 && value.ptr[0] == '"' && value.ptr[value.len - 1] == '"') {
            value = md_str_substr(value, 1, value.len - 2);
        }
        props = value;
    }

    int col = 0, species = -1, z = -1, pos = -1;
    md_str_t field[3];
    for (;;) {
        size_t n = 0;
        while (n < 3 && md_parse_next_token_delim(&field[n], &props, ':')) ++n;
        if (n == 0) break;
        if (n != 3 || !is_unsigned_int(field[2]) || field[1].len != 1) {
            MD_LOG_ERROR("XYZ: Invalid Properties in the extended XYZ comment line: '" MD_STR_FMT "'", MD_STR_ARG(value));
            return false;
        }
        const int count = (int)parse_int(field[2]);
        const char type = field[1].ptr[0];
        if (count == 1 && (type == 'S' || type == 's') && md_str_eq_ignore_case(field[0], MD_STR_LIT("species")) && species < 0) {
            species = col;
        } else if (count == 1 && (type == 'I' || type == 'i') && md_str_eq(field[0], MD_STR_LIT("Z")) && z < 0) {
            z = col;
        } else if (count == 3 && (type == 'R' || type == 'r') && (md_str_eq_ignore_case(field[0], MD_STR_LIT("pos")) || md_str_eq_ignore_case(field[0], MD_STR_LIT("positions"))) && pos < 0) {
            pos = col;
        }
        col += count;
        if (col > XYZ_MAX_COL) break;
    }

    if (pos < 0 || pos > XYZ_MAX_COL) {
        MD_LOG_ERROR("XYZ: The Properties of the extended XYZ comment line have no positions (pos:R:3)");
        return false;
    }
    int elem = species;
    if (elem < 0 && z >= 0) {
        elem = z;
        *flags |= XYZ_ELEMENT_Z;
    }
    if (elem < 0 || elem > XYZ_MAX_COL) {
        elem = 0;
        *flags |= XYZ_NO_ELEMENT;
    }
    *flags |= XYZ_COLUMNS(elem, pos);
    return true;
}

// A number with a Fortran exponent, "1.5D-03", as some older programs write them
static bool xyz_parse_fortran_f64(double* out, md_str_t tok) {
    char buf[64];
    if (tok.len >= sizeof(buf)) return false;
    md_memcpy(buf, tok.ptr, tok.len);
    for (size_t i = 0; i < tok.len; ++i) {
        if (buf[i] == 'D' || buf[i] == 'd') buf[i] = 'e';
    }
    return md_parse_f64(out, (md_str_t){ buf, tok.len }) == tok.len;
}

// A coordinate: a number, of a Fortran exponent as well
static bool xyz_is_float(md_str_t tok) {
    double v;
    return is_float(tok) || xyz_parse_fortran_f64(&v, tok);
}

static inline bool extract_flags(uint32_t* flags, md_reader_i* reader) {
    MD_ASSERT(flags);    

    // Extract first three lines
    md_str_t lines[3];
    bool extract_lines = (
        md_reader_extract_line(reader, &lines[0]) &&
        md_reader_extract_line(reader, &lines[1]) && 
        md_reader_extract_line(reader, &lines[2]));

    // Reset reader back to its original position
    md_reader_rewind(reader);
    
    if (!extract_lines) {
        MD_LOG_ERROR("Invalid format for XYZ: Failed to extract first three lines to determine format");
        return false;
    }

    {
        // Test if first line has an unsigned integer as first token (Should be universally applicable to all XYZ formats
        md_str_t token;
        md_str_t line = lines[0];
        if (!md_parse_next_token(&token, &line) || !is_unsigned_int(token)) {
            MD_LOG_ERROR("Invalid format for XYZ: Missing coordinate count");
            return false;
        }
    }

    // Extended XYZ: the comment line says what the columns are
    {
        md_str_t value;
        if (xyz_find_key(&value, lines[1], MD_STR_LIT("Properties")) || xyz_find_key(&value, lines[1], MD_STR_LIT("Lattice"))) {
            *flags |= XYZ_EXTENDED;
            if (!xyz_extended_columns(flags, lines[1])) {
                return false;
            }
            md_str_t tok;
            md_str_t line = lines[2];
            size_t n = 0;
            while (md_parse_next_token(&tok, &line)) ++n;
            if (n < XYZ_POS_COL(*flags) + 3) {
                MD_LOG_ERROR("Invalid format for XYZ: The atom lines have fewer columns than the Properties of the extended XYZ comment line say");
                return false;
            }
            return true;
        }
    }

    {
        // Test if we have an ARC trajectory
        // Second line should then contain exactly 6 floats
        md_str_t token[8];
        md_str_t line = lines[1];
        const int64_t num_tokens = md_parse_next_tokens(token, MD_COUNTOF(token), &line);

        if (num_tokens == 6 &&
            is_float(token[0]) && is_float(token[1]) && is_float(token[2]) &&
            is_float(token[3]) && is_float(token[4]) && is_float(token[5]))
        {
            *flags |= XYZ_ARC;
        }
    }

    {
        // Determine coordinate structure
        // Traditional XYZ only holds 4 fields (atomic number or element symbol and coordinates)
        // Tinker holds additional fields (atom index, atom name, atom type and connectivity): at least 6 fields
        md_str_t token[8];
        md_str_t line = lines[2];
        int64_t num_tokens = md_parse_next_tokens(token, MD_COUNTOF(token), &line);

        if (num_tokens >= 6 &&
            is_unsigned_int(token[0]) &&
            (is_atom_name(token[1]) || is_unsigned_int(token[1])) &&
            xyz_is_float(token[2]) && xyz_is_float(token[3]) && xyz_is_float(token[4]) &&
            is_unsigned_int(token[5]))
        {
            *flags |= XYZ_TINKER | XYZ_COLUMNS(1, 2);
            return true;
        } else if (num_tokens >= 4 &&
            (is_atom_name(token[0]) || is_unsigned_int(token[0])) &&
            xyz_is_float(token[1]) && xyz_is_float(token[2]) && xyz_is_float(token[3]))
        {
            // Ordinary XYZ, of which further columns are not read
            *flags |= XYZ_COLUMNS(0, 1);
            return true;
        }
    }

    MD_LOG_ERROR("Unrecognized XYZ format");
    return false;
}

// @NOTE(Robin):
// Following the information available at https://github.com/libAtoms/extxyz
// Which I assume is the correct specification for the file format?
// It seems that Lattices can be encoded in a number of ways:
// New style: string which is enclosed in quotes containing 9 whitespace separated floats
// Old style: [] enclosed list of 3 vectors, each vector is a list of 3 floats: [[1,2,3],[4,5,6],[7,8,9]]
// Old style: {} enclosed list 9 elements separated by whitespace: {1 2 3 4 5 6 7 8 9}
// The cell as written: row i is vector i, in whatever orientation.

static inline bool extract_extxyz_cell(double cell[3][3], md_str_t line) {
    md_str_t value;
    if (!xyz_find_key(&value, line, MD_STR_LIT("Lattice"))) {
        // Lattice information not found, since Lattice is optional, this is not an error
        return true;
    }

    if (value.len < 2) {
        MD_LOG_ERROR("XYZ: Missing lattice information");
        return false;
    }

    md_str_t tok[9];
    size_t num_tok = 0;
    const char open = value.ptr[0];
    const char close = value.ptr[value.len - 1];
    md_str_t inner = md_str_substr(value, 1, value.len - 2);

    if ((open == '"' && close == '"') || (open == '{' && close == '}')) {
        num_tok = md_parse_next_tokens(tok, MD_COUNTOF(tok), &inner);
        if (md_parse_next_token(&(md_str_t){0}, &inner)) {
            MD_LOG_ERROR("XYZ: Invalid number of tokens in Lattice encoding");
            return false;
        }
    } else if (open == '[' && close == ']') {
        // This is the cheeky one: [[1,2,3],[4,5,6],[7,8,9]]
        for (int i = 0; i < 3; ++i) {
            const char* b = inner.ptr;
            const char* e = inner.ptr + inner.len;
            while (b < e && *b != '[') ++b;
            const char* c = b;
            while (c < e && *c != ']') ++c;
            if (c >= e) {
                MD_LOG_ERROR("XYZ: Failed to extract Lattice vector");
                return false;
            }
            md_str_t vec = { b + 1, (size_t)(c - b - 1) };
            inner = (md_str_t){ c + 1, (size_t)(e - c - 1) };
            size_t num_sub_toks = md_parse_next_tokens_delim(tok + i*3, 3, &vec, ',');
            if (num_sub_toks != 3 || !md_str_empty(md_str_trim(vec))) {
                MD_LOG_ERROR("XYZ: Failed to extract Lattice vector");
                return false;
            }
            tok[i*3+0] = md_str_trim(tok[i*3+0]);
            tok[i*3+1] = md_str_trim(tok[i*3+1]);
            tok[i*3+2] = md_str_trim(tok[i*3+2]);
            num_tok += num_sub_toks;
        }
    } else {
        MD_LOG_ERROR("XYZ: Unrecognized Lattice encoding");
        return false;
    }

    double val[9] = {0};
    for (size_t i = 0; i < num_tok; ++i) {
        if (md_parse_f64(&val[i], tok[i]) != tok[i].len) {
            MD_LOG_ERROR("XYZ: Invalid number in Lattice: '" MD_STR_FMT "'", MD_STR_ARG(tok[i]));
            return false;
        }
    }
    md_memset(cell, 0, sizeof(double) * 9);
    if (num_tok == 3) {
        // The diagonal
        for (size_t i = 0; i < 3; ++i) cell[i][i] = val[i];
    } else if (num_tok == 9) {
        for (size_t i = 0; i < 9; ++i) cell[i/3][i%3] = val[i];
    } else {
        MD_LOG_ERROR("XYZ: Invalid number of tokens in Lattice encoding");
        return false;
    }

    return true;
}

// The element of an atom of the symbol or name sym. A plain XYZ file holds element symbols, in any case and some
// with a number after them ("C", "CL", "C12", "Cl1"); Tinker holds atom names of its force field ("CA" is a carbon)
static md_atomic_number_t xyz_element(md_str_t sym, bool tinker) {
    md_atomic_number_t z = md_util_element_lookup(sym, false);
    if (z) return z;
    if (!tinker) {
        z = md_util_element_lookup(sym, true);
        if (z) return z;
    }
    return md_atomic_number_infer_from_label(sym, (md_str_t){0}, tinker ? 3 : 0);
}

// The element column of a line: a symbol or name, else an atomic number
static void xyz_set_element(md_xyz_coordinate_t* coord, md_str_t tok, uint32_t flags) {
    if ((flags & XYZ_ELEMENT_Z) || !md_ascii_is_alpha(tok.ptr[0])) {
        const int64_t z = parse_int(tok);
        coord->atomic_number = (0 < z && z < 119) ? (int)z : 0;
        const md_str_t sym = coord->atomic_number ? md_util_element_symbol((md_element_t)coord->atomic_number) : tok;
        md_str_copy_to_char_buf(coord->element_symbol, sizeof(coord->element_symbol), sym);
    } else {
        md_str_copy_to_char_buf(coord->element_symbol, sizeof(coord->element_symbol), tok);
        coord->atomic_number = xyz_element(tok, flags & XYZ_TINKER);
    }
}

static inline bool extract_coord(md_xyz_coordinate_t* coord, md_str_t line, uint32_t flags) {
    MD_ASSERT(coord);

    const md_str_t original = line;

    if (flags & XYZ_TINKER) {
        // Index, name, x, y, z, then the atom type and the bonded atoms
        md_str_t tok[6 + 9];
        const size_t num_tokens = md_parse_next_tokens(tok, MD_COUNTOF(tok), &line);
        if (num_tokens < 5) {
            MD_LOG_ERROR("Invalid number of tokens in Tinker XYZ coordinate when parsing line: '" MD_STR_FMT "', expected >= 5, got %i", MD_STR_ARG(original), (int)num_tokens);
            return false;
        }
        coord->atom_index = (int)parse_int(tok[0]);
        xyz_set_element(coord, tok[1], flags);
        coord->x = (float)parse_float(tok[2]);
        coord->y = (float)parse_float(tok[3]);
        coord->z = (float)parse_float(tok[4]);
        if (num_tokens > 5) {
            coord->particle_type = (int)parse_int(tok[5]);
        }
        for (size_t i = 6; i < num_tokens; ++i) {
            if (!is_unsigned_int(tok[i])) {
                MD_LOG_ERROR("Invalid connectivity information in XYZ file when parsing line: '" MD_STR_FMT "', token: '" MD_STR_FMT "'", MD_STR_ARG(original), MD_STR_ARG(tok[i]));
                return false;
            }
            if (i - 6 < (size_t)XYZ_MAX_CONN) {
                coord->connectivity[i - 6] = (int)parse_int(tok[i]);
            }
        }
        if (md_parse_next_token(&(md_str_t){0}, &line)) {
            MD_LOG_ERROR("More bonded atoms than %d in Tinker XYZ line: '" MD_STR_FMT "'", XYZ_MAX_CONN, MD_STR_ARG(original));
            return false;
        }
        return true;
    }

    // The columns which are read, of those there are
    const size_t elem_col = (flags & XYZ_NO_ELEMENT) ? SIZE_MAX : XYZ_ELEM_COL(flags);
    const size_t pos_col  = XYZ_POS_COL(flags);
    const size_t need = MD_MAX(pos_col + 3, elem_col != SIZE_MAX ? elem_col + 1 : 0);
    md_str_t elem = {0};
    md_str_t pos[3];
    md_str_t tok;
    for (size_t col = 0; col < need; ++col) {
        if (!md_parse_next_token(&tok, &line)) {
            MD_LOG_ERROR("Invalid number of tokens in XYZ coordinate when parsing line: '" MD_STR_FMT "', expected >= %zu, got %zu", MD_STR_ARG(original), need, col);
            return false;
        }
        if (col == elem_col) elem = tok;
        if (col - pos_col < 3) pos[col - pos_col] = tok;
    }
    if (elem.len) {
        xyz_set_element(coord, elem, flags);
    }
    double v[3];
    for (int k = 0; k < 3; ++k) {
        if (md_parse_f64(&v[k], pos[k]) != pos[k].len && !xyz_parse_fortran_f64(&v[k], pos[k])) {
            MD_LOG_ERROR("Invalid coordinate '" MD_STR_FMT "' in XYZ line: '" MD_STR_FMT "'", MD_STR_ARG(pos[k]), MD_STR_ARG(original));
            return false;
        }
    }
    coord->x = (float)v[0];
    coord->y = (float)v[1];
    coord->z = (float)v[2];
    return true;
}

// The header of a model: the number of atoms, and the comment and cell. Its cell as the file gives it (row i vector i),
// zero without one.
static inline bool xyz_parse_model_header(md_xyz_model_t* model, md_reader_i* reader, uint32_t flags, size_t* coord_count, md_allocator_i* alloc) {
    MD_ASSERT(model);
    MD_ASSERT(reader);
    MD_ASSERT(coord_count);

    md_str_t line = {0};
    md_str_t tokens[8];

    if (!md_reader_extract_line(reader, &line)) {
        return false;
    }

    line = md_str_trim(line);
    if (md_str_empty(line)) {
        return false;
    }

    // The number of atoms leads the first line
    md_str_t rest = line;
    md_str_t count_tok;
    if (!md_parse_next_token(&count_tok, &rest) || !is_unsigned_int(count_tok)) {
        MD_LOG_ERROR("Failed to extract coordinate count in XYZ header: '" MD_STR_FMT "'", MD_STR_ARG(line));
        return false;
    }
    const int64_t count = parse_int(count_tok);
    if (count <= 0) {
        MD_LOG_ERROR("Failed to extract coordinate count in XYZ header");
        return false;
    }

    // Tinker has the comment after the number of atoms
    if ((flags & XYZ_TINKER) && (flags & XYZ_STORE_COMMENT)) {
        const md_str_t comment = md_str_trim(rest);
        if (!md_str_empty(comment)) {
            MD_ASSERT(alloc);
            model->comment = md_alloc_str(alloc, comment);
        }
    }

    // A comment line of its own, or Tinker's line of the cell
    if (!(flags & XYZ_TINKER) || (flags & XYZ_ARC)) {
        if (!md_reader_extract_line(reader, &line)) {
            MD_LOG_ERROR("XYZ: The model ends after its number of atoms");
            return false;
        }

        if (!(flags & XYZ_TINKER) && (flags & XYZ_STORE_COMMENT) && !md_str_empty(line)) {
            MD_ASSERT(alloc);
            model->comment = md_alloc_str(alloc, line);
        }
    }

    if (flags & XYZ_ARC) {
        const size_t num_tokens = md_parse_next_tokens(tokens, MD_COUNTOF(tokens), &line);
        
        if (num_tokens != 6) {
            MD_LOG_ERROR("Unexpected number of tokens (%zu) when parsing XYZ Arc Cell data", num_tokens);
            return false;
        }

        const double a = parse_float(tokens[0]);
        const double b = parse_float(tokens[1]);
        const double c = parse_float(tokens[2]);
        const double alpha = parse_float(tokens[3]);
        const double beta  = parse_float(tokens[4]);
        const double gamma = parse_float(tokens[5]);

        md_unitcell_t cell = md_unitcell_from_extent_and_angles(a, b, c, alpha, beta, gamma);
        md_unitcell_A_extract_double(model->cell, &cell);
    } else if (flags & XYZ_EXTENDED) {
        if (!extract_extxyz_cell(model->cell, line)) {
            MD_LOG_ERROR("Failed to extract cell data from line");
            return false;
        }
    }

    *coord_count = (size_t)count;

    return true;
}

// The cell of a model as md_unitcell_t holds it, and the rotation of its coordinates into it (the identity where the
// cell is in that orientation already, as one of Tinker always is): see md_unitcell_from_basis_vectors.
static bool xyz_model_cell(md_unitcell_t* out_cell, double R[3][3], const double cell[3][3]) {
    const double basis[3][3] = {
        { cell[0][0], cell[0][1], cell[0][2] },
        { cell[1][0], cell[1][1], cell[1][2] },
        { cell[2][0], cell[2][1], cell[2][2] },
    };
    return md_unitcell_from_basis_vectors(out_cell, R, basis);
}

static inline vec3_t xyz_rotate(double R[3][3], float x, float y, float z) {
    return vec3_set(
        (float)(R[0][0] * x + R[0][1] * y + R[0][2] * z),
        (float)(R[1][0] * x + R[1][1] * y + R[1][2] * z),
        (float)(R[2][0] * x + R[2][1] * y + R[2][2] * z));
}

static bool xyz_parse(md_xyz_data_t* data, md_reader_i* reader, md_allocator_i* alloc, bool stop_after_first_model) {
    MD_ASSERT(data);
    MD_ASSERT(reader);
    MD_ASSERT(alloc);
    
    uint32_t flags = 0;
    if (!extract_flags(&flags, reader)) {
        MD_LOG_ERROR("Parse XYZ: Invalid format");
        return false;
    }

    size_t expected_count = 0;
    size_t byte_offset = 0;

    for (;;) {
        md_xyz_model_t mdl = {0};
        if (!xyz_parse_model_header(&mdl, reader, flags | XYZ_STORE_COMMENT, &expected_count, alloc)) {
            break;
        }
        md_array_ensure(data->coordinates, md_array_size(data->coordinates) + expected_count, alloc);
        
        mdl.byte_offset = byte_offset;
        mdl.beg_coord_index = (uint32_t)md_array_size(data->coordinates);
        for (size_t i = 0; i < expected_count; ++i) {
            md_str_t line;
            md_xyz_coordinate_t coord = {0};
            if (md_reader_extract_line(reader, &line) &&
                extract_coord(&coord, line, flags))
            {
                md_array_push(data->coordinates, coord, alloc);
            } else {
                MD_LOG_ERROR("Parse XYZ, Failed to parse coordinate");
                return false;
            }
        }

        mdl.end_coord_index = (uint32_t)md_array_size(data->coordinates);
        md_array_push(data->models, mdl, alloc);

        if (stop_after_first_model) {
            break;
        }
        byte_offset = (size_t)md_reader_tell(reader);
    }

    data->num_coordinates = md_array_size(data->coordinates);
    data->num_models = md_array_size(data->models);

    if (data->num_models == 0) {
        MD_LOG_ERROR("Parse XYZ: No model in the file");
        return false;
    }

    return true;
}

// PUBLIC PROCEDURES

// The structure of XYZ files is as follows:
// Number of Coordinates (integer)
// Comment (string)
// Atom Coordinates: Symbol/Atomic Number (2 char or int), x (float), y (float), z (float)

// The aim is to identify the single line with number of coordinates and treat this as
// the header for each 'model' which could constitute a frame within an animation
// (given that the number of coordinates are the same)

bool md_xyz_data_parse_str(md_xyz_data_t* data, md_str_t str, struct md_allocator_i* alloc) {
    md_reader_i reader = md_reader_from_str(str);
    return xyz_parse(data, &reader, alloc, false);
}

bool md_xyz_data_parse_file(md_xyz_data_t* data, md_str_t filename, struct md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Parse XYZ: Failed to open file '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    size_t buf_cap = MD_MEGABYTES(1);
    
    md_file_reader_t reader;
    md_file_reader_init(&reader, file, buf_cap, md_temp_allocator(temp));
    bool result = xyz_parse(data, &reader.r, alloc, false);
    
    md_file_reader_close(&reader);
    md_temp_end(temp);
    md_file_close(&file);

    return result;
}

void md_xyz_data_free(md_xyz_data_t* data, struct md_allocator_i* alloc) {
    MD_ASSERT(data);
    if (data->coordinates) md_array_free(data->coordinates, alloc);
    if (data->models) {
        for (size_t i = 0; i < data->num_models; ++i) {
            if (!md_str_empty(data->models[i].comment)) {
                md_free_str(alloc, data->models[i].comment);
            }
        }
        md_array_free(data->models, alloc);
    }
    md_memset(data, 0, sizeof(md_xyz_data_t));
}

// The bonds of Tinker's connectivity, of the atoms [beg, end) of data: each pair once (it is given by both its atoms),
// by the indices the file gives its atoms (1 to N in order, but which are looked up all the same)
static void xyz_tinker_bonds(md_system_t* sys, const md_xyz_data_t* data, size_t beg, size_t end, md_allocator_i* temp) {
    const size_t N = end - beg;
    // The atom of each index of the file, where the indices are 1 to N in order, which is all Tinker writes
    bool sequential = true;
    for (size_t i = 0; i < N; ++i) {
        if (data->coordinates[beg + i].atom_index != (int)(i + 1)) { sequential = false; break; }
    }
    md_hashmap32_t index_map = { .allocator = temp };
    if (!sequential) {
        for (size_t i = 0; i < N; ++i) {
            md_hashmap_add(&index_map, (uint64_t)(uint32_t)data->coordinates[beg + i].atom_index, (uint32_t)i);
        }
    }

    md_array(uint64_t) keys = 0;
    for (size_t i = 0; i < N; ++i) {
        const md_xyz_coordinate_t* c = &data->coordinates[beg + i];
        for (int k = 0; k < XYZ_MAX_CONN && c->connectivity[k] > 0; ++k) {
            uint32_t j = UINT32_MAX;
            if (sequential) {
                if ((size_t)c->connectivity[k] <= N) j = (uint32_t)(c->connectivity[k] - 1);
            } else {
                const uint32_t* found = md_hashmap_get(&index_map, (uint64_t)(uint32_t)c->connectivity[k]);
                if (found) j = *found;
            }
            if (j == UINT32_MAX || j == i) continue;
            md_array_push(keys, ((uint64_t)MD_MIN(i, j) << 32) | MD_MAX(i, j), temp);
        }
    }
    const size_t num_keys = md_array_size(keys);
    if (num_keys == 0) return;

    uint32_t* order = md_alloc(temp, num_keys * sizeof(uint32_t));
    md_util_argsort_radix_u64(order, keys, num_keys);
    const md_bond_flags_t file = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_FILE);
    md_array_ensure(sys->bond.pairs, num_keys, sys->alloc);
    md_array_ensure(sys->bond.flags, num_keys, sys->alloc);
    for (size_t i = 0; i < num_keys; ++i) {
        const uint64_t key = keys[order[i]];
        if (i > 0 && key == keys[order[i - 1]]) continue;
        const md_atom_pair_t pair = {{ (md_atom_idx_t)(key >> 32), (md_atom_idx_t)(key & 0xFFFFFFFF) }};
        md_array_push(sys->bond.pairs, pair, sys->alloc);
        md_array_push(sys->bond.flags, file, sys->alloc);
        sys->bond.count += 1;
    }
    md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
}

// What makes a system of the atoms of a model, one atom at a time
typedef struct xyz_builder_t {
    md_system_t* sys;
    md_system_state_t* state;
    md_hashmap32_t type_map;
    // The types of the last symbols seen, by a hash of the symbol and element: most lines repeat one of a few
    uint64_t cache_key[16];
    md_particle_type_idx_t cache_type[16];
    bool rotate;
    double R[3][3];
} xyz_builder_t;

static void xyz_builder_begin(xyz_builder_t* b, md_system_t* sys, md_system_state_t* state, size_t num_atoms, const double cell[3][3], md_allocator_i* temp) {
    md_memset(b, 0, sizeof(*b));
    b->sys = sys;
    b->state = state;
    b->type_map.allocator = temp;

    md_system_reset(sys);
    md_system_state_init(state, 0);

    const size_t reserve_size = MD_ALIGN_TO(num_atoms, 16);
    md_array_ensure(state->xyz, reserve_size, state->alloc);
    md_array_ensure(sys->atom.type_idx, reserve_size, sys->alloc);
    md_array_ensure(sys->atom.flags, reserve_size, sys->alloc);

    // Setup particle types including unknown type
    md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    // A cell in another orientation than md_unitcell_t holds turns the coordinates along with it
    b->rotate = xyz_model_cell(&state->unitcell, b->R, cell);
}

static void xyz_builder_add(xyz_builder_t* b, const md_xyz_coordinate_t* c) {
    md_system_t* sys = b->sys;
    uint32_t sym_word;
    md_memcpy(&sym_word, c->element_symbol, sizeof(sym_word));
    const uint64_t key = ((uint64_t)sym_word << 32) | (uint32_t)c->atomic_number;
    const size_t slot = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 60);
    md_particle_type_idx_t type;
    if (b->cache_key[slot] == key && b->cache_type[slot]) {
        type = b->cache_type[slot];
    } else {
        // The symbol as written, not the whole buffer: the rest of it is zeros, which a name would carry along
        size_t symbol_len = 0;
        while (symbol_len < sizeof(c->element_symbol) && c->element_symbol[symbol_len]) ++symbol_len;
        const md_str_t symbol = { c->element_symbol, symbol_len };
        type = md_util_particle_type_find_or_add(sys, &b->type_map, symbol, (md_atomic_number_t)c->atomic_number);
        b->cache_key[slot] = key;
        b->cache_type[slot] = type;
    }

    sys->atom.count += 1;
    md_array_push_no_grow(b->state->xyz, b->rotate ? xyz_rotate(b->R, c->x, c->y, c->z) : vec3_set(c->x, c->y, c->z));
    md_array_push_no_grow(sys->atom.flags, 0);
    md_array_push_no_grow(sys->atom.type_idx, type);
}

static void xyz_builder_end(xyz_builder_t* b) {
    MD_ASSERT(md_array_size(b->state->xyz) == b->sys->atom.count);
    b->state->num_atoms = b->sys->atom.count;
}

bool md_xyz_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_xyz_data_t* data, md_xyz_options_t options) {
    MD_ASSERT(sys);
    MD_ASSERT(state);
    MD_ASSERT(data);
    (void)options;

    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }

    if (!state || !state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    if (data->num_models == 0 || !data->models) {
        MD_LOG_ERROR("XYZ: No model to make a system of");
        return false;
    }

    // If we have more than one model, it is a trajectory: the structure is its first model
    const size_t beg_coord_index = data->models[0].beg_coord_index;
    const size_t end_coord_index = data->models[0].end_coord_index;

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    xyz_builder_t b;
    xyz_builder_begin(&b, sys, state, end_coord_index - beg_coord_index, data->models[0].cell, md_temp_allocator(temp));

    bool tinker_bonds = false;
    for (size_t i = beg_coord_index; i < end_coord_index; ++i) {
        xyz_builder_add(&b, &data->coordinates[i]);
        tinker_bonds |= data->coordinates[i].connectivity[0] > 0;
    }

    // Tinker's connectivity is the bonds
    if (tinker_bonds) {
        xyz_tinker_bonds(sys, data, beg_coord_index, end_coord_index, md_temp_allocator(temp));
    }

    xyz_builder_end(&b);
    md_temp_end(temp);
    return true;
}

// The system of the first model the reader reads. Its atoms go into the system as they are read, but for Tinker's,
// whose bonds are made of all of them: those are read as data first.
static bool xyz_system_load(md_system_t* sys, md_system_state_t* state, md_reader_i* reader, md_xyz_options_t options) {
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state || !state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    uint32_t flags = 0;
    if (!extract_flags(&flags, reader)) {
        MD_LOG_ERROR("Parse XYZ: Invalid format");
        goto done;
    }

    if (flags & XYZ_TINKER) {
        md_xyz_data_t data = {0};
        result = xyz_parse(&data, reader, temp_alloc, true) && md_xyz_system_init_from_data(sys, state, &data, options);
        goto done;
    }

    md_xyz_model_t model = {0};
    size_t count = 0;
    if (!xyz_parse_model_header(&model, reader, flags, &count, NULL)) {
        MD_LOG_ERROR("Parse XYZ: No model in the file");
        goto done;
    }

    xyz_builder_t b;
    xyz_builder_begin(&b, sys, state, count, model.cell, temp_alloc);
    for (size_t i = 0; i < count; ++i) {
        md_str_t line;
        md_xyz_coordinate_t coord;
        md_memset(coord.element_symbol, 0, sizeof(coord.element_symbol));
        coord.atomic_number = 0;
        if (!md_reader_extract_line(reader, &line) || !extract_coord(&coord, line, flags)) {
            MD_LOG_ERROR("Parse XYZ, Failed to parse coordinate");
            md_system_reset(sys);
            md_system_state_init(state, 0);
            goto done;
        }
        xyz_builder_add(&b, &coord);
    }
    xyz_builder_end(&b);
    result = true;

done:
    md_temp_end(temp);
    return result;
}

bool md_xyz_system_init_from_str(md_system_t* sys, md_system_state_t* state, md_str_t str, md_xyz_options_t options) {
    MD_ASSERT(sys);
    md_reader_i reader = md_reader_from_str(str);
    // The structure is the first model: of a trajectory, the rest are its frames (md_xyz_system_publish_run)
    return xyz_system_load(sys, state, &reader, options);
}

bool md_xyz_system_init_from_file(md_system_t* sys, md_system_state_t* state, md_str_t filename, md_xyz_options_t options) {
    MD_ASSERT(sys);
    
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Failed to open file '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_file_reader_t reader;
    md_file_reader_init(&reader, file, MD_MEGABYTES(1), md_temp_allocator(temp_scope));
    
    // The structure is the first model: of a trajectory, the rest are its frames (md_xyz_system_publish_run)
    const bool result = xyz_system_load(sys, state, &reader.r, options);

    md_file_reader_close(&reader);
    md_temp_end(temp_scope);
    md_file_close(&file);

    return result;
}

typedef struct xyz_cache_t {
    md_run_cache_header_t header;
    int64_t* offsets;   // num_frames + 1
    float*   cells;     // 9 per frame, Angstrom, row i box vector i, turned as md_unitcell_t holds it; zero without a cell
} xyz_cache_t;

// The cache beside path when it was made from the file as it is now: header, num_frames + 1 offsets,
// 9 cell floats per frame.
static bool try_read_cache(xyz_cache_t* cache, md_str_t path, md_allocator_i* alloc) {
    MD_ASSERT(cache);
    MD_ASSERT(alloc);
    md_file_t file = {0};
    if (!md_run_cache_open(&file, &cache->header, path, MD_XYZ_CACHE_MAGIC, MD_XYZ_CACHE_VERSION)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    const size_t cell_bytes   = cache->header.num_frames * 9 * sizeof(float);
    cache->offsets = md_alloc(alloc, offset_bytes);
    cache->cells   = md_alloc(alloc, cell_bytes);
    const bool ok =
        md_file_read(file, cache->offsets, offset_bytes) == offset_bytes &&
        md_file_read(file, cache->cells,   cell_bytes)   == cell_bytes &&
        md_file_tell(file) == (int64_t)md_file_size(file);
    if (!ok) {
        MD_LOG_ERROR("The XYZ cache beside '" MD_STR_FMT "' is incomplete", MD_STR_ARG(path));
        md_free(alloc, cache->offsets, offset_bytes);
        md_free(alloc, cache->cells,   cell_bytes);
        cache->offsets = NULL;
        cache->cells   = NULL;
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const xyz_cache_t* cache, md_str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    if (!md_run_cache_create(&file, path, scanned, MD_XYZ_CACHE_MAGIC, MD_XYZ_CACHE_VERSION, cache->header.num_atoms, cache->header.num_frames)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    const size_t cell_bytes   = cache->header.num_frames * 9 * sizeof(float);
    const bool ok =
        md_file_write(file, cache->offsets, offset_bytes) == offset_bytes &&
        md_file_write(file, cache->cells,   cell_bytes)   == cell_bytes;
    if (!ok) {
        MD_LOG_ERROR("Failed to write the XYZ cache beside '" MD_STR_FMT "'", MD_STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// The file's layout, where each model starts and the cell each states: from the cache beside the file
// when that is current, from a scan of the file otherwise (writing the cache unless told not to). The scan reads the
// header of each model and skips its atom lines.
// From alloc; the scan's scratch keeps clear of avoid, which is what alloc was made from.
static bool xyz_index_load(xyz_cache_t* cache, uint32_t* out_xyz_flags, md_str_t filename, md_run_flags_t run_flags, md_allocator_i* alloc, md_allocator_i* avoid) {
    md_memset(cache, 0, sizeof(*cache));

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Failed to open file for XYZ trajectory");
        return false;
    }

    md_temp_scope_t temp_scope = md_temp_begin_avoid(avoid);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    bool result = false;

    md_file_reader_t reader;
    md_file_reader_init(&reader, file, MD_MEGABYTES(1), temp_alloc);

    uint32_t xyz_flags = 0;
    if (!extract_flags(&xyz_flags, &reader.r)) {
        MD_LOG_ERROR("Failed to determine format for XYZ trajectory");
        goto done;
    }
    // Comments are never kept for a frame
    xyz_flags &= ~XYZ_STORE_COMMENT;
    *out_xyz_flags = xyz_flags;

    // The file as it is before the scan is what a new cache is stamped with
    md_file_info_t scanned = {0};
    md_file_info_extract(file, &scanned);
    const int64_t filesize = (int64_t)scanned.size;

    if (try_read_cache(cache, filename, alloc)) {
        result = true;
        goto done;
    }
    md_memset(cache, 0, sizeof(*cache));

    md_array(int64_t) offsets = 0;
    md_array(float)   cells   = 0;
    size_t num_atoms = 0;
    int64_t frames_end = 0;
    for (;;) {
        const int64_t offset = md_reader_tell(&reader.r);
        md_xyz_model_t model = {0};
        size_t count = 0;
        if (!xyz_parse_model_header(&model, &reader.r, xyz_flags, &count, NULL)) {
            break;
        }
        if (offsets && count != num_atoms) {
            MD_LOG_ERROR("The XYZ file models are not of equal length and cannot be read as a trajectory");
            goto done;
        }
        num_atoms = count;
        if (md_reader_skip_lines(&reader.r, count) != count) {
            MD_LOG_ERROR("The last model of the XYZ file is cut short: it is left out of the trajectory");
            break;
        }
        md_unitcell_t cell;
        xyz_model_cell(&cell, NULL, model.cell);
        float A[3][3] = {0};
        md_unitcell_A_extract_float(A, &cell);
        md_array_push(offsets, offset, temp_alloc);
        md_array_push_array(cells, &A[0][0], 9, temp_alloc);
        frames_end = md_reader_tell(&reader.r);
    }
    const size_t num_frames = md_array_size(offsets);
    if (num_frames <= 1) {
        md_log(MD_LOG_TYPE_INFO, "The XYZ file did not contain multiple entries and cannot be read as a trajectory");
        goto done;
    }

    cache->header.num_atoms = num_atoms;
    cache->header.num_frames = num_frames;
    cache->offsets = md_alloc(alloc, sizeof(int64_t) * (num_frames + 1));
    cache->cells   = md_alloc(alloc, sizeof(float) * 9 * num_frames);
    md_memcpy(cache->offsets, offsets, sizeof(int64_t) * num_frames);
    md_memcpy(cache->cells, cells, sizeof(float) * 9 * num_frames);
    // The last frame ends with its last atom line
    cache->offsets[num_frames] = MD_MIN(filesize, frames_end);

    if (!(run_flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        // A cache that cannot be written only costs the next load a parse
        write_cache(cache, filename, &scanned);
    }
    result = true;

done:
    md_file_reader_close(&reader);
    md_temp_end(temp_scope);
    md_file_close(&file);
    return result;
}

// ### RUN ###

// The coordinates of atoms [first, first + count) of one frame's text, packed xyz, turned with the frame's cell where
// md_unitcell_t holds it in another orientation. Returns how many were written.
static size_t xyz_parse_frame_coords(float* xyz, md_str_t text, uint32_t flags, size_t first, size_t count) {
    md_reader_i reader = md_reader_from_str(text);
    md_xyz_model_t model = {0};
    size_t coord_count = 0;
    if (!xyz_parse_model_header(&model, &reader, flags, &coord_count, NULL)) {
        return 0;
    }
    md_unitcell_t cell;
    double R[3][3];
    const bool rotate = xyz_model_cell(&cell, R, model.cell);

    // Only the columns of the coordinates are read
    const uint32_t col_flags = (flags & ~(uint32_t)XYZ_TINKER) | XYZ_NO_ELEMENT;

    size_t written = 0;
    if (first >= coord_count) return 0;
    if (md_reader_skip_lines(&reader, first) != first) return 0;
    md_str_t line;
    for (size_t i = first; written < count && i < coord_count && md_reader_extract_line(&reader, &line); ++i) {
        md_xyz_coordinate_t coord = {0};
        if (!extract_coord(&coord, line, col_flags)) {
            MD_LOG_ERROR("Error when decoding coordinate");
            return written;
        }
        const vec3_t p = rotate ? xyz_rotate(R, coord.x, coord.y, coord.z) : vec3_set(coord.x, coord.y, coord.z);
        xyz[written * 3 + 0] = p.x;
        xyz[written * 3 + 1] = p.y;
        xyz[written * 3 + 2] = p.z;
        written += 1;
    }
    return written;
}

static size_t xyz_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    MD_ASSERT(sys);
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, MD_STR_LIT("atom/position"))) return 0;

    const int32_t* layout = (const int32_t*)md_attribute_view(md_attributes_find_in(&sys->attributes, src.run, MD_STR_LIT("source/layout")), MD_ATTRIBUTE_TYPE_I32, 1, 0);
    if (!layout) {
        MD_LOG_ERROR("XYZ: the run '" MD_STR_FMT "' has lost its layout", MD_STR_ARG(src.run));
        return 0;
    }
    const uint32_t flags = (uint32_t)layout[0];

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
    char* text = md_temp_alloc(temp, MD_MAX(frame_size, 1));
    if (text && md_attribute_io_read_at(io, src.path, src.offset[frame], text, frame_size) == frame_size) {
        if (xyz_parse_frame_coords((float*)dst, (md_str_t){ text, frame_size }, flags, first, count) == count) {
            written = cap;
        } else {
            MD_LOG_ERROR("XYZ: frame %u of '" MD_STR_FMT "' has fewer than %zu atoms", frame, MD_STR_ARG(src.path), first + count);
        }
    } else {
        MD_LOG_ERROR("XYZ: Failed to read frame %u from '" MD_STR_FMT "'", frame, MD_STR_ARG(src.path));
    }
    md_temp_end(temp);
    return written;
}

bool md_xyz_system_publish_run(md_system_t* sys, md_str_t filename, md_str_t run, uint32_t flags) {
    MD_ASSERT(sys);
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const md_str_t path = { path_buf, path_len };

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    bool result = false;

    xyz_cache_t index;
    uint32_t xyz_flags = 0;
    if (path_len == 0 || !xyz_index_load(&index, &xyz_flags, path, flags, arena, md_get_heap_allocator())) {
        goto done;
    }

    const size_t F = index.header.num_frames;
    double*  times = md_alloc(arena, F * sizeof(double));
    int64_t* sizes = md_alloc(arena, F * sizeof(int64_t));
    for (size_t i = 0; i < F; ++i) {
        times[i] = (double)i;
        sizes[i] = index.offsets[i + 1] - index.offsets[i];
    }

    const md_attribute_virtual_t virt = { .provider = xyz_position_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = index.header.num_atoms,
        .time          = times,
        .time_unit     = md_unit_none(),    // nothing in the file says when: ordinals
        .unitcell      = index.cells,
        .source_path   = path,
        .source_offset = index.offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    if (!md_run_publish(sys, run, &desc)) {
        goto done;
    }

    char buf[512];
    const int32_t layout = (int32_t)xyz_flags;
    if (!md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/layout")),
        .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 0 },
        .unit = md_unit_none(), .description = MD_STR_INIT("Which XYZ dialect the frames are written in"),
        .data = &layout, .byte_size = sizeof(layout)})) {
        md_attributes_remove_prefix(&sys->attributes, run);
        goto done;
    }
    result = true;

done:
    md_arena_allocator_destroy(arena);
    return result;
}

#ifdef __cplusplus
}
#endif
