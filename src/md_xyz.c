#include <md_xyz.h>

#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_XYZ_CACHE_MAGIC      0x8265485749172bab
#define MD_XYZ_CACHE_VERSION    5
#define MD_XYZ_MOL_MAGIC        0x285ada29078a9bc8

enum {
    XYZ_TINKER          = 1,
    XYZ_ARC             = 2,
    XYZ_EXTENDED		= 4,
    XYZ_STORE_COMMENT   = 8,
};

// The opaque blob
typedef struct xyz_molecule_t {
    uint64_t magic;
    md_allocator_i* allocator;
} xyz_molecule_t;

// We massage the beg and end indices here to correspond the xyz specification
// This makes our life easier when specifying all the different ranges
static inline int32_t extract_int(str_t line, size_t beg, size_t end) {
    if (line.len < end) return 0;
    return (int32_t)parse_int(str_trim(str_substr(line, beg, end-beg)));
}

static inline float extract_float(str_t line, size_t beg, size_t end) {
    if (line.len < end) return 0.0f;
    return (float)parse_float(str_trim(str_substr(line, beg, end-beg)));
}

static inline bool is_unsigned_int(str_t str) {
    str = str_trim(str);
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    while (c < end) {
        if (!is_digit(*c)) return false;
        ++c;
    }
    return true;
}

static inline bool is_string(str_t str) {
    str = str_trim(str);

    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    while (c < end) {
        if (!is_alpha(*c)) return false;
        ++c;
    }
    return true;
}

static inline bool extract_flags(uint32_t* flags, md_buffered_reader_t* reader) {
    ASSERT(flags);    

    // Extract first three lines
    str_t lines[3];
    bool extract_lines = (
        md_buffered_reader_extract_line(&lines[0], reader) &&
        md_buffered_reader_extract_line(&lines[1], reader) && 
        md_buffered_reader_extract_line(&lines[2], reader));

    // Reset reader back to its original position
    md_buffered_reader_reset(reader);
    
    if (!extract_lines) {
        MD_LOG_ERROR("Invalid format for XYZ: Failed to extract first three lines to determine format");
        return false;
    }

    {
        // Test if first line has an unsigned integer as first token (Should be universally applicable to all XYZ formats
        str_t token;
        str_t line = lines[0];
        if (!extract_token(&token, &line) || !is_unsigned_int(token)) {
            MD_LOG_ERROR("Invalid format for XYZ: Missing coordinate count");
            return false;
        }
    }

    // Test for extended XYZ format
    if (str_find_str(NULL, lines[1], STR_LIT("Properties="))) {
        *flags |= XYZ_EXTENDED;
    } else {    
        // Test if we have an ARC trajectory
        // Second line should then contain exactly 6 floats
        str_t token[8];
        str_t line = lines[1];
        const int64_t num_tokens = extract_tokens(token, ARRAY_SIZE(token), &line);

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
        // Tinker holds additional fields (atom index, atom type and connectivity) at least 6 fields
        str_t token[8];
        str_t line = lines[2];
        int64_t num_tokens = extract_tokens(token, ARRAY_SIZE(token), &line);

        if (num_tokens >= 6 &&
            is_unsigned_int(token[0]) &&
            (is_string(token[1]) || is_unsigned_int(token[1])) &&
            is_float(token[2]) && is_float(token[3]) && is_float(token[4]) &&
            is_unsigned_int(token[5]))
        {
            *flags |= XYZ_TINKER;
            return true;
        } else if (num_tokens >= 4 &&
            (is_string(token[0]) || is_unsigned_int(token[0])) &&
            is_float(token[1]) && is_float(token[2]) && is_float(token[3]))
        {
            // Ordinary XYZ (with potential bond info)
            return true;
        }
    }

    MD_LOG_ERROR("Unrecognized XYZ format");
    return false;
}

static inline str_t extract_quoted_substr(str_t in_str) {
    const char* beg = in_str.ptr;
    const char* end = in_str.ptr + in_str.len;

    str_t result = {0};

    if (beg >= end) return result;

    while (beg < end && *beg != '\"') {
        ++beg;
    }

    if (beg >= end) return result;

    ++beg;
    const char* c = beg;
    while (c < end && *c != '\"') {
        ++c;
    }

    if (c < end) {
        result.ptr = beg;
        result.len = (int64_t)(c - beg);
    }

    return result;
}

static inline str_t extract_balanced_substr(str_t in_str, char beg_char, char end_char) {
    const char* beg = in_str.ptr;
    const char* end = in_str.ptr + in_str.len;

    ASSERT(beg_char != end_char);

    str_t result = {0};

    if (beg >= end) return result;

    while (beg < end && *beg != beg_char) {
        ++beg;
    }

    if (beg >= end) return result;

    ++beg;
    const char* c = beg;
    int depth = 1;
    while (c < end && depth > 0) {
        if (*c == beg_char) {
            ++depth;
        }
        else if (*c == end_char) {
            --depth;
        }
        ++c;
    }

    if (depth == 0) {
        result.ptr = beg;
        result.len = (int64_t)(c - beg - 1);
    }

    return result;
}

// @NOTE(Robin):
// Following the information available at https://github.com/libAtoms/extxyz
// Which I assume is the correct specification for the file format?
// It seems that Lattices can be encoded in a number of ways:
// New style: string which is enclosed in quotes containing 9 whitespace separated floats
// Old style: [] enclosed list of 3 vectors, each vector is a list of 3 floats: [[1,2,3],[4,5,6],[7,8,9]]
// Old style: {} enclosed list 9 elements separated by whitespace: {1 2 3 4 5 6 7 8 9}

static inline bool extract_extxyz_cell(float cell[3][3], str_t line) {
    const str_t pattern = STR_LIT("Lattice=");
    size_t loc;
    if (!str_find_str(&loc, line, pattern)) {
        // Lattice information not found, since Lattice is optional, this is not an error
        return true;
    }

    line = str_substr(line, loc + pattern.len, SIZE_MAX);
    if (line.len == 0) {
        MD_LOG_ERROR("Missing lattice information");
        return false;
    }

    str_t tok[9];
    size_t num_tok = 0;

    if (line.ptr[0] == '\"') {
        line = extract_quoted_substr(line);
        if (str_empty(line)) {
            MD_LOG_ERROR("XYZ: Failed to extract quoted string");
            return false;
        }
        num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
    } else if (line.ptr[0] == '{') {
        line = extract_balanced_substr(line, '{', '}');
        if (str_empty(line)) {
            MD_LOG_ERROR("XYZ: Failed to extract curly-brace string");
            return false;
        }
        num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
    } else if (line.ptr[0] == '[') {
        // This is the cheeky one
        line = extract_balanced_substr(line, '[', ']');
        if (str_empty(line)) {
            MD_LOG_ERROR("XYZ: Failed to extract bracketed string");
            return false;
        }
        for (int i = 0; i < 3; ++i) {
            str_t vec = extract_balanced_substr(line, '[', ']');
            if (str_empty(vec)) {
                MD_LOG_ERROR("XYZ: Failed to extract bracketed string");
                return false;
            }
            line = str_substr(line, vec.len + 2, SIZE_MAX);
            size_t num_sub_toks = extract_tokens_delim(tok + i*3, ARRAY_SIZE(tok) - i*3, &vec, ',');
            if (num_sub_toks != 3) {
                MD_LOG_ERROR("XYZ: Failed to extract Lattice vector");
                return false;
            }
            tok[i*3+0] = str_trim(tok[i*3+0]);
            tok[i*3+1] = str_trim(tok[i*3+1]);
            tok[i*3+2] = str_trim(tok[i*3+2]);
            num_tok += num_sub_toks;
        }
    } else {
        MD_LOG_ERROR("XYZ: Unrecognized Lattice encoding");
        return false;
    }

    if (num_tok == 3) {
        // Assume these encode the diagonal
        for (size_t i = 0; i < 3; ++i) {
            if (is_float(tok[i])) {
                cell[i][i] = (float)parse_float(tok[i]);
            } else {
                return false;
            }
        }
    } else if (num_tok == 9) {
        for (size_t i = 0; i < num_tok; ++i) {
            if (is_float(tok[i])) {
                cell[i/3][i%3] = (float)parse_float(tok[i]);
            } else {
                return false;
            }
        }
    } else {
        MD_LOG_ERROR("XYZ: Invalid number of tokens in Lattice encoding");
        return false;
    }

    return true;
}

static inline bool extract_coord(md_xyz_coordinate_t* coord, str_t line, uint32_t flags) {
    ASSERT(coord);

    str_t original = line;

    str_t tokens[16];
    const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);

    if (num_tokens < 4) {
        MD_LOG_ERROR("Invalid number of tokens in XYZ coordinate when parsing line: '%.*s', expected >= 4, got %i", (int)original.len, original.ptr, (int)num_tokens);
        return false;
    }

    int tok_idx = 0;
    if ((flags & XYZ_TINKER) && num_tokens > 4) {
        coord->atom_index = (int)parse_int(tokens[tok_idx++]);
    }

    ASSERT(!str_empty(tokens[tok_idx]));
    if (is_alpha(tokens[tok_idx].ptr[0])) {
        size_t len = str_copy_to_char_buf(coord->element_symbol, sizeof(coord->element_symbol), tokens[tok_idx++]);
        str_t symbol = {coord->element_symbol, len};
        coord->atomic_number = md_util_element_lookup(symbol, false);
    } else {
        coord->atomic_number = (int)parse_int(tokens[tok_idx++]);
    }

    coord->x = (float)parse_float(tokens[tok_idx++]);
    coord->y = (float)parse_float(tokens[tok_idx++]);
    coord->z = (float)parse_float(tokens[tok_idx++]);
    
    if (!(flags & XYZ_EXTENDED) && tok_idx < (int)num_tokens) {
        coord->atom_type = (int)parse_int(tokens[tok_idx++]);

        // Connectivity follows atom_type
        for (int i = tok_idx; i < (int)num_tokens; ++i) {
            if (!is_digit(tokens[i].ptr[0])) {
                MD_LOG_ERROR("Invalid connectivity information in XYZ file when parsing line: '%.*s', token: '%.*s', number of tokens: %i", (int)original.len, original.ptr, (int)tokens[i].len, tokens[i].ptr, (int)num_tokens);
                return false;
            }
            coord->connectivity[i-tok_idx] = (int)parse_int(tokens[i]);
        }
    }

    return true;
}

static inline bool xyz_parse_model_header(md_xyz_model_t* model, md_buffered_reader_t* reader, uint32_t flags, size_t* coord_count, md_allocator_i* alloc) {
    ASSERT(model);
    ASSERT(reader);
    ASSERT(coord_count);

    str_t line = {0};
    str_t tokens[8];
    int64_t count = 0;

    if (model->byte_offset == 261828) {
        while(0);
    }

    if (!md_buffered_reader_extract_line(&line, reader)) {
        return false;
    }

    line = str_trim(line);
    if (str_empty(line)) {
        return false;
    }

    // Parse data from first line, we only need the first two tokens even though more may exist
    const size_t num_tok = extract_tokens(tokens, 2, &line);
    if (num_tok) {
        count = parse_int(tokens[0]);
    }

    if (count <= 0) {
        MD_LOG_ERROR("Failed to extract coordinate count in XYZ header");  
        return false;
    }


    if ((flags & XYZ_TINKER)) {
        // Comment encoded after the first token
        str_t comment = str_substr(line, (ptrdiff_t)(tokens[1].ptr - line.ptr), SIZE_MAX);

        if ((flags & XYZ_STORE_COMMENT) && !str_empty(comment)) {
            ASSERT(alloc);
            model->comment = str_copy(comment, alloc);
        }
    }

    // Arc files are also flagged as tinker, but only pure tinker files uses 1 line for its header.
    bool pure_tinker = (flags & XYZ_TINKER) && !(flags & XYZ_ARC);
    if (!pure_tinker) {
        // Tinker is the only format that only uses 1 line for the header
        // The others uses 2 lines
        if (!md_buffered_reader_extract_line(&line, reader)) {
            MD_LOG_ERROR("Failed to extract extra line");
        }

        str_t comment = line;
        if ((flags & XYZ_STORE_COMMENT) && !str_empty(comment)) {
            ASSERT(alloc);
            model->comment = str_copy(comment, alloc);
        }
    }

    if (flags & XYZ_ARC) {
        const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
        
        if (num_tokens != 6) {
            MD_LOG_ERROR("Unexpected number of tokens (%zu) when parsing XYZ Arc Cell data", num_tokens);
            return false;
        }

        double extent[3];
        double angle[3];

        extent[0] = (float)parse_float(tokens[0]);
        extent[1] = (float)parse_float(tokens[1]);
        extent[2] = (float)parse_float(tokens[2]);

        angle[0] = (float)parse_float(tokens[3]);
        angle[1] = (float)parse_float(tokens[4]);
        angle[2] = (float)parse_float(tokens[5]);

        md_unitcell_t cell = md_unitcell_from_extent_and_angles(extent[0], extent[1], extent[2], angle[0], angle[1], angle[2]);
        mat3_t A = { 0 };
        md_unitcell_A_extract_float(A.elem, &cell);
        MEMCPY(model->cell, A.elem, sizeof(model->cell));
    } else if (flags & XYZ_EXTENDED) {
        // Extract cell data from line
        if (!extract_extxyz_cell(model->cell, line)) {
            MD_LOG_ERROR("Failed to extract cell data from line");
            return false;
        }
    }

    *coord_count = count;

    return true;
}

/*
static inline int32_t xyz_parse_model_coordinates(md_xyz_data_t* data, str_t* str, int32_t count, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(str);
    ASSERT(count);
    ASSERT(alloc);

    int32_t i = 0;
    str_t line;
    for (; i < count; ++i) {
        if (!str_extract_line(&line, str))
            break;
        md_xyz_coordinate_t coord;
        if (extract_coord(&coord, line)) {
            md_array_push(data->coordinates, coord, alloc);
        } else {
            MD_LOG_ERROR("Failed to parse model coordinate");
            break;
        }
    }

    return i;
}
*/

bool xyz_parse(md_xyz_data_t* data, md_buffered_reader_t* reader, md_allocator_i* alloc, bool stop_after_first_model) {
    ASSERT(data);
    ASSERT(reader);
    ASSERT(alloc);
    
    uint32_t flags = 0;
    if (!extract_flags(&flags, reader)) {
        MD_LOG_ERROR("Parse XYZ: Invalid format");
        return false;
    }

    size_t expected_count = 0;
    md_xyz_model_t mdl = {0};
    size_t byte_offset = 0;

    while (xyz_parse_model_header(&mdl, reader, flags | XYZ_STORE_COMMENT, &expected_count, alloc)) {
        md_array_ensure(data->coordinates, md_array_size(data->coordinates) + expected_count, alloc);
        
        mdl.byte_offset = byte_offset;
        mdl.beg_coord_index = (uint32_t)md_array_size(data->coordinates);
        for (size_t i = 0; i < expected_count; ++i) {
            str_t line;
            md_xyz_coordinate_t coord = {0};
            if (md_buffered_reader_extract_line(&line, reader) &&
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
        byte_offset = (size_t)md_buffered_reader_tellg(reader);
    }

    data->num_coordinates = md_array_size(data->coordinates);
    data->num_models = md_array_size(data->models);

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

bool md_xyz_data_parse_str(md_xyz_data_t* data, str_t str, struct md_allocator_i* alloc) {
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);
    return xyz_parse(data, &reader, alloc, false);
}

bool md_xyz_data_parse_file(md_xyz_data_t* data, str_t filename, struct md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Parse XYZ: Failed to open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    size_t buf_cap = MEGABYTES(1);
    char* buf = md_temp_alloc(temp, buf_cap);
    
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, buf_cap, file);
    bool result = xyz_parse(data, &reader, alloc, false);
    
    md_temp_end(temp);
    md_file_close(&file);

    return result;
}

void md_xyz_data_free(md_xyz_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    if (data->coordinates) md_array_free(data->coordinates, alloc);
    if (data->models) {
        for (size_t i = 0; i < data->num_models; ++i) {
            if (!str_empty(data->models[i].comment)) {
                str_free(data->models[i].comment, alloc);
            }
        }
        md_array_free(data->models, alloc);
    }
    MEMSET(data, 0, sizeof(md_xyz_data_t));
}

bool md_xyz_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_xyz_data_t* data, md_xyz_options_t options) {
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

    size_t beg_coord_index = 0;
    size_t end_coord_index = data->num_coordinates;

    // if we have more than one model, interperet it as a trajectory and only load the first model
    if (data->num_models > 0) {
        // Limit the scope of atom coordinate entries if we have a trajectory (only consider first model)
        beg_coord_index = data->models[0].beg_coord_index;
        end_coord_index = data->models[0].end_coord_index;
    }

    const size_t num_atoms = end_coord_index - beg_coord_index;
    const size_t reserve_size = ALIGN_TO(num_atoms, 16);

    md_array_ensure(state->xyz, reserve_size, state->alloc);
    md_array_ensure(sys->atom.type_idx, reserve_size, sys->alloc);

    // Setup atom types including unknown type
    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    for (size_t i = beg_coord_index; i < end_coord_index; ++i) {
        float x = data->coordinates[i].x;
        float y = data->coordinates[i].y;
        float z = data->coordinates[i].z;
        str_t atom_symbol = {data->coordinates[i].element_symbol, sizeof(data->coordinates[i].element_symbol)};
        md_atomic_number_t atomic_number = (md_atomic_number_t)data->coordinates[i].atomic_number;
        float mass = md_atomic_number_mass(atomic_number);
        float radius = md_atomic_number_vdw_radius(atomic_number);
        uint32_t color = md_atomic_number_cpk_color(atomic_number);
        md_atom_type_idx_t atom_type_idx = md_atom_type_find_or_add(&sys->atom.type, atom_symbol, atomic_number, mass, radius, color, 0, sys->alloc);

        sys->atom.count += 1;
        md_array_push(state->xyz, vec3_set(x, y, z), state->alloc);
        md_array_push(sys->atom.flags, 0, sys->alloc);
        md_array_push(sys->atom.type_idx, atom_type_idx, sys->alloc);
    }

    state->unitcell = md_unitcell_from_matrix_float(MD_AS_CONST_MAT3(data->models[0].cell));

    ASSERT(md_array_size(state->xyz) == sys->atom.count);
    state->num_atoms = sys->atom.count;

    return true;
}

bool md_xyz_system_init_from_str(md_system_t* sys, md_system_state_t* state, str_t str, md_xyz_options_t options) {
    ASSERT(sys);
    
    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);

    md_xyz_data_t data = {0};
    bool result = xyz_parse(&data, &reader, temp_arena, false) && md_xyz_system_init_from_data(sys, state, &data, options);

    md_temp_end(temp_scope);
    return result;
}

bool md_xyz_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename, md_xyz_options_t options) {
    ASSERT(sys);
    
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Failed to open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
    size_t buf_cap = MEGABYTES(1);
    char* buf = md_temp_alloc(temp_scope, buf_cap);

    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, buf_cap, file);
    
    md_xyz_data_t data = {0};
    bool result = xyz_parse(&data, &reader, temp_arena, false) && md_xyz_system_init_from_data(sys, state, &data, options);

    // Several frames are a trajectory as well: md_xyz_system_publish_run makes them a run.

    md_temp_end(temp_scope);
    md_file_close(&file);

    return result;
}

typedef struct xyz_cache_t {
    md_run_cache_header_t header;
    int64_t* offsets;   // num_frames + 1
    float*   cells;     // 9 per frame, Angstrom, row i box vector i; zero without a cell
} xyz_cache_t;

// The cache beside path when it was made from the file as it is now: header, num_frames + 1 offsets,
// 9 cell floats per frame.
static bool try_read_cache(xyz_cache_t* cache, str_t path, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);
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
        MD_LOG_ERROR("The XYZ cache beside '" STR_FMT "' is incomplete", STR_ARG(path));
        md_free(alloc, cache->offsets, offset_bytes);
        md_free(alloc, cache->cells,   cell_bytes);
        cache->offsets = NULL;
        cache->cells   = NULL;
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const xyz_cache_t* cache, str_t path, const md_file_info_t* scanned) {
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
        MD_LOG_ERROR("Failed to write the XYZ cache beside '" STR_FMT "'", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// The file's layout, where each model starts and the cell each states: from the cache beside the file
// when that is current, from a parse of the file otherwise (writing the cache unless told not to).
// From alloc; the parse's scratch keeps clear of avoid, which is what alloc was made from.
static bool xyz_index_load(xyz_cache_t* cache, uint32_t* out_xyz_flags, str_t filename, md_run_flags_t run_flags, md_allocator_i* alloc, md_allocator_i* avoid) {
    MEMSET(cache, 0, sizeof(*cache));

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Failed to open file for XYZ trajectory");
        return false;
    }

    uint32_t xyz_flags = 0;
    {
        char buf[1024];
        size_t len = md_file_read(file, buf, sizeof(buf));
        md_buffered_reader_t reader = md_buffered_reader_from_str((str_t){buf, len});
        if (!extract_flags(&xyz_flags, &reader)) {
            MD_LOG_ERROR("Failed to determine format for XYZ trajectory");
            md_file_close(&file);
            return false;
        }
    }
    // Comments are never kept for a frame
    xyz_flags &= ~XYZ_STORE_COMMENT;
    *out_xyz_flags = xyz_flags;

    // The file as it is before the parse is what a new cache is stamped with
    md_file_info_t scanned = {0};
    md_file_info_extract(file, &scanned);
    const int64_t filesize = (int64_t)scanned.size;
    md_file_close(&file);

    if (try_read_cache(cache, filename, alloc)) {
        return true;
    }
    MEMSET(cache, 0, sizeof(*cache));

    md_temp_scope_t temp_scope = md_temp_begin_avoid(avoid);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    bool result = false;

    md_xyz_data_t data = {0};
    if (!md_xyz_data_parse_file(&data, filename, temp_alloc)) {
        goto done;
    }
    if (data.num_models <= 1) {
        md_log(MD_LOG_TYPE_INFO, "The XYZ file did not contain multiple entries and cannot be read as a trajectory");
        goto done;
    }

    // Validate the models, pick the atom count in the first model and ensure that all other models have the same number of atoms
    const int64_t num_atoms = data.models[0].end_coord_index - data.models[0].beg_coord_index;
    for (size_t i = 1; i < data.num_models; ++i) {
        const int64_t length = data.models[i].end_coord_index - data.models[i].beg_coord_index;
        if (length != num_atoms) {
            MD_LOG_ERROR("The XYZ file models are not of equal length and cannot be read as a trajectory");
            goto done;
        }
    }

    cache->header.num_atoms = num_atoms;
    cache->header.num_frames = data.num_models;
    cache->offsets = md_alloc(alloc, sizeof(int64_t) * (cache->header.num_frames + 1));
    cache->cells   = md_alloc(alloc, sizeof(float) * 9 * cache->header.num_frames);

    for (size_t i = 0; i < data.num_models; ++i) {
        cache->offsets[i] = data.models[i].byte_offset;
        MEMCPY(cache->cells + i * 9, data.models[i].cell, 9 * sizeof(float));
    }
    cache->offsets[data.num_models] = filesize;

    if (!(run_flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        // A cache that cannot be written only costs the next load a parse
        write_cache(cache, filename, &scanned);
    }
    result = true;

done:
    md_temp_end(temp_scope);
    return result;
}

// ### RUN ###

// The coordinates of atoms [first, first + count) of one frame's text, packed xyz. Returns how many
// were written.
static size_t xyz_parse_frame_coords(float* xyz, str_t text, uint32_t flags, size_t first, size_t count) {
    md_buffered_reader_t reader = md_buffered_reader_from_str(text);
    md_xyz_model_t model = {0};
    size_t coord_count = 0;
    if (!xyz_parse_model_header(&model, &reader, flags, &coord_count, NULL)) {
        return 0;
    }

    size_t i = 0;
    size_t written = 0;
    str_t line;
    str_t tokens[8];
    while (written < count && i < coord_count && md_buffered_reader_extract_line(&line, &reader)) {
        if (line.len < 6) continue;
        const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);
        if (num_tokens < 4) {
            MD_LOG_ERROR("Error when decoding coordinate");
            return written;
        }
        // Tinker lines lead with an index and a name, plain ones with the element alone
        const size_t offset = num_tokens > 4 ? 2 : 1;
        if (i >= first) {
            xyz[written * 3 + 0] = (float)parse_float(tokens[offset + 0]);
            xyz[written * 3 + 1] = (float)parse_float(tokens[offset + 1]);
            xyz[written * 3 + 2] = (float)parse_float(tokens[offset + 2]);
            written += 1;
        }
        i += 1;
    }
    return written;
}

static size_t xyz_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    ASSERT(sys);
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;

    char buf[512];
    const md_attribute_t* layout = md_attributes_find(&sys->attributes, md_run_path(buf, sizeof(buf), src.run, STR_LIT("source/layout")));
    if (!layout || layout->format.type != MD_ATTRIBUTE_TYPE_I32 || !layout->data) {
        MD_LOG_ERROR("XYZ: the run '" STR_FMT "' has lost its layout", STR_ARG(src.run));
        return 0;
    }
    const uint32_t flags = (uint32_t)((const int32_t*)layout->data)[0];

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
    char* text = md_temp_alloc(temp, MAX(frame_size, 1));
    if (text && md_attribute_io_read_at(io, src.path, src.offset[frame], text, frame_size) == frame_size) {
        if (xyz_parse_frame_coords((float*)dst, (str_t){ text, frame_size }, flags, first, count) == count) {
            written = cap;
        } else {
            MD_LOG_ERROR("XYZ: frame %u of '" STR_FMT "' has fewer than %zu atoms", frame, STR_ARG(src.path), first + count);
        }
    } else {
        MD_LOG_ERROR("XYZ: Failed to read frame %u from '" STR_FMT "'", frame, STR_ARG(src.path));
    }
    md_temp_end(temp);
    return written;
}

bool md_xyz_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const str_t path = { path_buf, path_len };

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
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
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/layout")),
        .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 0 },
        .unit = md_unit_none(), .description = STR_LIT("Which XYZ dialect the frames are written in"),
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
