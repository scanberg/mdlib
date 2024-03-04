#include <md_xyz.h>

#include <md_molecule.h>
#include <md_trajectory.h>
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
#define MD_XYZ_CACHE_VERSION    2
#define MD_XYZ_MOL_MAGIC        0x285ada29078a9bc8
#define MD_XYZ_TRAJ_MAGIC       0x2312ad7b78a9bc78

enum {
    XYZ_TINKER          = 1,
    XYZ_ARC             = 2,
    XYZ_EXTENDED		= 4,
};

// The opaque blob
typedef struct xyz_molecule_t {
    uint64_t magic;
    md_allocator_i* allocator;
} xyz_molecule_t;

typedef struct xyz_trajectory_t {
    uint64_t magic;
    md_file_o* file;
    int64_t* frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;
    uint32_t flags;
} xyz_trajectory_t;

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
    int64_t num_tok = 0;

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
            int64_t num_sub_toks = extract_tokens_delim(tok + i*3, ARRAY_SIZE(tok) - i*3, &vec, ',');
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
        for (int64_t i = 0; i < 3; ++i) {
            if (is_float(tok[i])) {
				cell[i][i] = (float)parse_float(tok[i]);
			} else {
				return false;
			}
		}
    } else if (num_tok == 9) {
        for (int64_t i = 0; i < num_tok; ++i) {
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
    const int64_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);

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
        str_copy_to_char_buf(coord->element_symbol, sizeof(coord->element_symbol), tokens[tok_idx++]);
    } else {
        coord->atomic_number = (int)parse_int(tokens[tok_idx++]);
    }

    coord->x = (float)parse_float(tokens[tok_idx++]);
    coord->y = (float)parse_float(tokens[tok_idx++]);
    coord->z = (float)parse_float(tokens[tok_idx++]);
    
    if (!(flags & XYZ_EXTENDED) && tok_idx < num_tokens) {
        coord->atom_type = (int)parse_int(tokens[tok_idx++]);

        // Connectivity follows atom_type
        for (int i = tok_idx; i < num_tokens; ++i) {
            if (!is_digit(tokens[i].ptr[0])) {
                MD_LOG_ERROR("Invalid connectivity information in XYZ file when parsing line: '%.*s', token: '%.*s', number of tokens: %i", (int)original.len, original.ptr, (int)tokens[i].len, tokens[i].ptr, (int)num_tokens);
                return false;
            }
            coord->connectivity[i-tok_idx] = (int)parse_int(tokens[i]);
        }
    }

    return true;
}

static inline bool xyz_parse_model_header(md_xyz_model_t* model, md_buffered_reader_t* reader, uint32_t flags, size_t* coord_count) {
    ASSERT(model);
    ASSERT(reader);
    ASSERT(coord_count);

    str_t line[2] = {0};
    str_t tokens[8];
    str_t comment = {0};
    int32_t count = 0;

    if (!md_buffered_reader_extract_line(&line[0], reader)) {
        return false;
    }

    str_t trimmed = str_trim(line[0]);
    if (trimmed.len == 0) {
        return false;
    }

    // Check for exact match here
    // Arc files are also flagged as tinker, but only pure tinker files uses 1 line for its header.
    if (!(flags == XYZ_TINKER)) {
        // Tinker is the only format that only uses 1 line for the header
        // The others uses 2 lines
        if (!md_buffered_reader_extract_line(&line[1], reader)) {
            MD_LOG_ERROR("Failed to extract extra line");
        }
    }

    // Parse data from first line, we only need the first two tokens even though more may exist
    const size_t num_tok = extract_tokens(tokens, 2, &line[0]);
    if (num_tok) {
        count = (int32_t)parse_int(tokens[0]);
    }

    if (count <= 0) {
        MD_LOG_ERROR("Failed to extract coordinate count in XYZ header");  
        return false;
    }

    if (flags & XYZ_TINKER) {
        // Comment encoded after the first token
        comment = str_substr(line[0], (int64_t)(tokens[1].ptr - line[0].ptr), SIZE_MAX);
    } else {
        // Second line is the comment
        comment = line[1];
    }

    if (comment.len > 0) {

        str_copy_to_char_buf(model->comment, sizeof(model->comment), comment);
    }

    if (flags & XYZ_ARC) {
        const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line[1]);
        
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

        md_unit_cell_t cell = md_util_unit_cell_from_extent_and_angles(extent[0], extent[1], extent[2], angle[0], angle[1], angle[2]);
        MEMCPY(model->cell, cell.basis.elem, sizeof(model->cell));
    } else if (flags & XYZ_EXTENDED) {
		// Extract cell data from comment
		if (!extract_extxyz_cell(model->cell, comment)) {
			MD_LOG_ERROR("Failed to extract cell data from comment");
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

bool xyz_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    ASSERT(xyz);
    ASSERT(xyz->magic == MD_XYZ_TRAJ_MAGIC);
    ASSERT(header);

    *header = xyz->header;
    return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
// Returns size in bytes of frame, frame_data_ptr is optional and is the destination to write the frame data to.
size_t xyz_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    ASSERT(xyz);
    ASSERT(xyz->magic == MD_XYZ_TRAJ_MAGIC);

    if (!xyz->file) {
        MD_LOG_ERROR("File handle is NULL");
        return 0;
    }

    if (!xyz->frame_offsets) {
        MD_LOG_ERROR("Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || (int64_t)xyz->header.num_frames <= frame_idx) {
        MD_LOG_ERROR("Frame index is out of range");
        return 0;
    }

    const int64_t beg = xyz->frame_offsets[frame_idx + 0];
    const int64_t end = xyz->frame_offsets[frame_idx + 1];
    const size_t frame_size = (size_t)(end - beg);
    const size_t total_size = sizeof(int64_t) + frame_size;

    if (frame_data_ptr) {
        // Store the index to the frame since this is generally not found within the actual frame data
        int64_t* ptr = (int64_t*)frame_data_ptr;
        ptr[0] = frame_idx;

        ASSERT(xyz->file);
        md_mutex_lock(&xyz->mutex);
        md_file_seek(xyz->file, beg, MD_FILE_BEG);
        const size_t bytes_read = md_file_read(xyz->file, &ptr[1], frame_size);
        (void)bytes_read;
        md_mutex_unlock(&xyz->mutex);
        ASSERT(frame_size == bytes_read);
    }

    return total_size;
}

bool xyz_decode_frame_data(struct md_trajectory_o* inst, const void* data_ptr, size_t data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    if (!data_ptr) {
        MD_LOG_ERROR("Data pointer is NULL");
		return false;
    }

    if (data_size == 0) {
        MD_LOG_ERROR("Data size is zero");
        return false;
    }

    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    if (xyz->magic != MD_XYZ_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame header, xyz magic did not match");
        return false;
    }

    const int64_t step = *((int64_t*)data_ptr);
    if (step < 0 || step >= (int64_t)xyz->header.num_frames) {
        MD_LOG_ERROR("Error when decoding frame data, corrupt frame index");
        return false;
    }

    str_t str = { .ptr = (const char*)(data_ptr) + sizeof(int64_t), .len = data_size - sizeof(int64_t) };
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);

    md_xyz_model_t model = {0};
    size_t coord_count = 0;
    if (!xyz_parse_model_header(&model, &reader, xyz->flags, &coord_count)) {
        MD_LOG_ERROR("Error when decoding header");
        return false;
    }

    size_t i = 0;
    str_t line;
    str_t tokens[8];
    while (md_buffered_reader_extract_line(&line, &reader) && i < xyz->header.num_atoms) {
        if (line.len < 6) continue;

        const size_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);

        if (num_tokens < 4) {
            MD_LOG_ERROR("Error when decoding coordinate");
            return false;
        }
        
        int offset = 1;
        if (num_tokens > 4) {
            offset += 1;
        }
        
        if (x) x[i] = (float)parse_float(tokens[offset + 0]);
        if (y) y[i] = (float)parse_float(tokens[offset + 1]);
        if (z) z[i] = (float)parse_float(tokens[offset + 2]);

        i += 1;
    }

    if (header) {
        header->num_atoms = i;
        header->index = step;
        header->timestamp = (double)(step); // This information is missing from xyz trajectories
        //header->unit_cell = md_util_unit_cell_from_extent_and_angles(model.cell_extent[0], model.cell_extent[1], model.cell_extent[2], model.cell_angle[0], model.cell_angle[1], model.cell_angle[2]);
        header->unit_cell = md_util_unit_cell_from_matrix(model.cell);
    }

    return true;
}

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
    int64_t byte_offset = 0;

    while (xyz_parse_model_header(&mdl, reader, flags, &expected_count)) {
        md_array_ensure(data->coordinates, md_array_size(data->coordinates) + expected_count, alloc);
        
        mdl.byte_offset = byte_offset;
        mdl.beg_coord_index = (int32_t)md_array_size(data->coordinates);
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
        mdl.end_coord_index = (int32_t)md_array_size(data->coordinates);
        md_array_push(data->models, mdl, alloc);

        if (stop_after_first_model) {
            break;
        }
        byte_offset = md_buffered_reader_tellg(reader);
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
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int64_t cap = MEGABYTES(1);
        char* buf = md_alloc(md_get_heap_allocator(), cap);
        ASSERT(buf);
        
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        result = xyz_parse(data, &reader, alloc, false);
        
        md_free(md_get_heap_allocator(), buf, cap);
        md_file_close(file);
    } else {
        MD_LOG_ERROR("Parse XYZ: Failed to open file '%.*s'", (int)filename.len, filename.ptr);
    }
    return result;
}

void md_xyz_data_free(md_xyz_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    if (data->coordinates) md_array_free(data->coordinates, alloc);
    if (data->models) md_array_free(data->models, alloc);
    MEMSET(data, 0, sizeof(md_xyz_data_t));
}

bool md_xyz_molecule_init(md_molecule_t* mol, const md_xyz_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(data);
    ASSERT(alloc);

    int64_t beg_coord_index = 0;
    int64_t end_coord_index = data->num_coordinates;

    // if we have more than one model, interperet it as a trajectory and only load the first model
    if (data->num_models > 0) {
        // Limit the scope of atom coordinate entries if we have a trajectory (only consider first model)
        beg_coord_index = data->models[0].beg_coord_index;
        end_coord_index = data->models[0].end_coord_index;
    }

    MEMSET(mol, 0, sizeof(md_molecule_t));

    const size_t num_atoms = end_coord_index - beg_coord_index;

    md_array_ensure(mol->atom.x, num_atoms, alloc);
    md_array_ensure(mol->atom.y, num_atoms, alloc);
    md_array_ensure(mol->atom.z, num_atoms, alloc);
    md_array_ensure(mol->atom.element, num_atoms, alloc);

    for (int64_t i = beg_coord_index; i < end_coord_index; ++i) {
        float x = data->coordinates[i].x;
        float y = data->coordinates[i].y;
        float z = data->coordinates[i].z;
        str_t atom_type = {data->coordinates[i].element_symbol, sizeof(data->coordinates[i].element_symbol)};
        md_element_t element = (md_element_t)data->coordinates[i].atomic_number;

        mol->atom.count += 1;
        md_array_push(mol->atom.x, x, alloc);
        md_array_push(mol->atom.y, y, alloc);
        md_array_push(mol->atom.z, z, alloc);
        md_array_push(mol->atom.element, element, alloc);
        md_array_push(mol->atom.type, make_label(atom_type), alloc);
        md_array_push(mol->atom.flags, 0, alloc);
    }

    mol->unit_cell = md_util_unit_cell_from_matrix(data->models[0].cell);

    return true;
}

static bool xyz_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
    (void)arg;

    md_xyz_data_t data = {0};
    
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);
    bool result = xyz_parse(&data, &reader, md_get_heap_allocator(), true);
    result = result && md_xyz_molecule_init(mol, &data, alloc);
    md_xyz_data_free(&data, md_get_heap_allocator());
    
    return result;
}

static bool xyz_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_xyz_data_t data = {0};
    
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int64_t cap = MEGABYTES(1);
        char* buf = md_alloc(md_get_heap_allocator(), cap);
        ASSERT(buf);

        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        result = xyz_parse(&data, &reader, md_get_heap_allocator(), true);
        result = result && md_xyz_molecule_init(mol, &data, alloc);

        md_xyz_data_free(&data, md_get_heap_allocator());
        md_free(md_get_heap_allocator(), buf, cap);
        md_file_close(file);
    } else {
        MD_LOG_ERROR("Init XYZ from file: Failed to open file '%.*s'", (int)filename.len, filename.ptr);
    }
    return result;
}

static md_molecule_loader_i xyz_molecule_api = {
    xyz_init_from_str,
    xyz_init_from_file,
};

md_molecule_loader_i* md_xyz_molecule_api() {
    return &xyz_molecule_api;
}

bool xyz_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    if (xyz->magic != MD_XYZ_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    bool result = true;
    const size_t frame_size = xyz_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        md_allocator_i* alloc = frame_size > md_temp_allocator_max_allocation_size() ? md_get_heap_allocator() : md_get_temp_allocator();
        void* frame_data = md_alloc(alloc, frame_size);
        const size_t read_size = xyz_fetch_frame_data(inst, frame_idx, frame_data);
        if (read_size != frame_size) {
            MD_LOG_ERROR("Failed to read the expected size");
            md_free(alloc, frame_data, frame_size);
            return false;
        }

        result = xyz_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
        md_free(alloc, frame_data, frame_size);
    }

    return result;
}

typedef struct xyz_cache_t {
    md_trajectory_cache_header_t header;
    int64_t* offsets;
} xyz_cache_t;

static bool try_read_cache(xyz_cache_t* cache, str_t cache_file, size_t traj_num_bytes, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);

    bool result = false;
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        if (md_file_read(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
            MD_LOG_ERROR("XYZ trajectory cache: failed to read header");
            goto done;
        }

        if (cache->header.magic != MD_XYZ_CACHE_MAGIC) {
        	MD_LOG_ERROR("XYZ trajectory cache: magic was incorrect or corrupt");
        	goto done;
        }
        if (cache->header.version != MD_XYZ_CACHE_VERSION) {
            MD_LOG_INFO("XYZ trajectory cache: version mismatch, expected %i, got %i", MD_XYZ_CACHE_VERSION, (int)cache->header.version);
        }
        if (cache->header.num_bytes != traj_num_bytes) {
			MD_LOG_INFO("XYZ trajectory cache: trajectory size mismatch, expected %zu, got %zu", traj_num_bytes, cache->header.num_bytes);
		}
        if (cache->header.num_atoms == 0) {
			MD_LOG_ERROR("XYZ trajectory cache: num atoms was zero");
			goto done;
		}
        if (cache->header.num_frames == 0) {
            MD_LOG_ERROR("XYZ trajectory cache: num frames was zero");
            goto done;
        }

        const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
        cache->offsets = md_alloc(alloc, offset_bytes);
        if (md_file_read(file, cache->offsets, offset_bytes) != offset_bytes) {
            MD_LOG_ERROR("Failed to read offset cache, offsets are incomplete");
            md_free(alloc, cache->offsets, offset_bytes);
            goto done;
        }

        result = true;
    done:
        md_file_close(file);
    }
    return result;
}

static bool write_cache(const xyz_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_INFO("XYZ trajectory cache: could not open file '"STR_FMT"'", STR_ARG(cache_file));
        return false;
    }

    if (md_file_write(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
		MD_LOG_ERROR("XYZ trajectory cache: failed to write header");
		goto done;
	}

    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    if (md_file_write(file, cache->offsets, offset_bytes) != offset_bytes) {
        MD_LOG_ERROR("Failed to write offset cache, offsets");
        goto done;
    }

    result = true;

done:
    md_file_close(file);
    return result;
}

md_trajectory_i* md_xyz_trajectory_create(str_t filename, md_allocator_i* ext_alloc, uint32_t traj_flags) {
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (!file) {
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
            return false;
        }
    }

    int64_t filesize = md_file_size(file);
    md_file_close(file);

    char buf[1024] = "";
    int len = snprintf(buf, sizeof(buf), "%.*s.cache", (int)filename.len, filename.ptr);
    str_t cache_file = {buf, (size_t)len};

    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    xyz_cache_t cache = {0};
    if (!try_read_cache(&cache, cache_file, filesize, alloc)) {
        md_allocator_i* temp_alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

        bool result = false;
        md_xyz_data_t data = {0};
        if (!md_xyz_data_parse_file(&data, filename, temp_alloc)) {
            goto cleanup;
        }

        if (data.num_models <= 1) {
            md_log(MD_LOG_TYPE_INFO, "The XYZ file did not contain multiple entries and cannot be read as a trajectory");
            md_xyz_data_free(&data, temp_alloc);
            goto cleanup;
        }

        // Validate the models, pick the atom count in the first model and ensure that all other models have the same number of atoms
        const int64_t num_atoms = data.models[0].end_coord_index - data.models[0].beg_coord_index;
        for (size_t i = 1; i < data.num_models; ++i) {
            const int64_t length = data.models[i].end_coord_index - data.models[i].beg_coord_index;
            if (length != num_atoms) {
                MD_LOG_ERROR("The XYZ file models are not of equal length and cannot be read as a trajectory");
                goto cleanup;
            }
        }

        cache.header.magic = MD_XYZ_CACHE_MAGIC;
        cache.header.version = MD_XYZ_CACHE_VERSION;
        cache.header.num_bytes = filesize;
        cache.header.num_atoms = num_atoms;
        cache.header.num_frames = data.num_models;
        cache.offsets = md_alloc(alloc, sizeof(int64_t) * (cache.header.num_frames + 1));

        for (size_t i = 0; i < data.num_models; ++i) {
            cache.offsets[i] = data.models[i].byte_offset;
        }
        cache.offsets[data.num_models] = filesize;

        if (!(traj_flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
            // If we fail to write the cache, that's ok, we can inform about it, but do not halt
            if (write_cache(&cache, cache_file)) {
                MD_LOG_INFO("XYZ: Successfully created cache file for '" STR_FMT "'", STR_ARG(cache_file));
            }
        }

        result = true;
    cleanup:
        md_arena_allocator_destroy(temp_alloc);
        if (!result) {
            return false;
        }
    }

    size_t max_frame_size = 0;
    for (size_t i = 0; i < cache.header.num_frames; ++i) {
        const size_t frame_size = (size_t)MAX(0, cache.offsets[i + 1] - cache.offsets[i]);
        max_frame_size = MAX(max_frame_size, frame_size);
    }

    void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));
    ASSERT(mem);
    MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));

    md_trajectory_i* traj = mem;
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)(traj + 1);

    md_array(double) frame_times = md_array_create(double, cache.header.num_frames, alloc);
    for (size_t i = 0; i < cache.header.num_frames; ++i) {
        frame_times[i] = (double)i;
    }

    xyz->magic = MD_XYZ_TRAJ_MAGIC;
    xyz->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    xyz->frame_offsets = cache.offsets;
    xyz->allocator = alloc;
    xyz->mutex = md_mutex_create();
    xyz->header = (md_trajectory_header_t) {
        .num_frames = cache.header.num_frames,
        .num_atoms = cache.header.num_atoms,
        .max_frame_data_size = max_frame_size,
        .time_unit = {0},
        .frame_times = frame_times,
    };
    xyz->flags = xyz_flags;

    traj->inst = (struct md_trajectory_o*)xyz;
    traj->get_header = xyz_get_header;
    traj->load_frame = xyz_load_frame;
    //traj->fetch_frame_data = xyz_fetch_frame_data;
    //traj->decode_frame_data = xyz_decode_frame_data;

    return traj;
}

void md_xyz_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)traj->inst;
    if (xyz->magic != MD_XYZ_TRAJ_MAGIC) {
        MD_LOG_ERROR("Trajectory is not a valid xyz trajectory.");
        ASSERT(false);
        return;
    }
    
    if (xyz->file) md_file_close(xyz->file);
    md_mutex_destroy(&xyz->mutex);
    md_arena_allocator_destroy(xyz->allocator);
}

static md_trajectory_loader_i xyz_traj_loader = {
    md_xyz_trajectory_create,
    md_xyz_trajectory_free,
};

md_trajectory_loader_i* md_xyz_trajectory_loader() {
    return &xyz_traj_loader;
}

#ifdef __cplusplus
}
#endif
