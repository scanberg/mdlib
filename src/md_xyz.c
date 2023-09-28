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

#define MD_XYZ_CACHE_MAGIC 0x8265485749172bab
#define MD_XYZ_MOL_MAGIC   0x285ada29078a9bc8
#define MD_XYZ_TRAJ_MAGIC  0x2312ad7b78a9bc78

enum {
    XYZ_TINKER          = 1,
    XYZ_ARC             = 2,
};

// The opaque blob
typedef struct xyz_molecule_t {
    uint64_t magic;
    md_allocator_i* allocator;
} xyz_molecule_t;

typedef struct xyz_trajectory_t {
    uint64_t magic;
    md_file_o* file;
    uint64_t filesize;
    int64_t* frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;
    uint32_t flags;
} xyz_trajectory_t;

// We massage the beg and end indices here to correspond the xyz specification
// This makes our life easier when specifying all the different ranges
static inline int32_t extract_int(str_t line, int64_t beg, int64_t end) {
    if (line.len < end) return 0;
    return (int32_t)parse_int(str_trim(str_substr(line, beg, end-beg)));
}

static inline float extract_float(str_t line, int64_t beg, int64_t end) {
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
        // Test if first line has an unsigned integer as first token (Should be universally applicable to all XYZ + ARC formats
        str_t token;
        str_t line = lines[0];
        if (!extract_token(&token, &line) || !is_unsigned_int(token)) {
            MD_LOG_ERROR("Invalid format for XYZ: Missing coordinate count");
            return false;
        }
    }

    {
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
        const int64_t num_tokens = extract_tokens(token, ARRAY_SIZE(token), &line);

        if (num_tokens >= 6){
            if (is_unsigned_int(token[0]) &&
                (is_string(token[1]) || is_unsigned_int(token[1])) &&
                is_float(token[2]) && is_float(token[3]) && is_float(token[4]) &&
                is_unsigned_int(token[5]))
            {
                *flags |= XYZ_TINKER;
            } else {
                MD_LOG_ERROR("Invalid format for XYZ: Unknown variation of Tinker format");
                return false;
            }
        }
        else if (num_tokens == 4) {
            if (!(is_string(token[0]) || is_unsigned_int(token[0])) ||
                !(is_float(token[1]) && is_float(token[2]) && is_float(token[3])))
            {
                MD_LOG_ERROR("Invalid format for XYZ: Unknown variation of format");
                return false;
            }
        } else {
            MD_LOG_ERROR("Invalid format for XYZ");
            return false;
        }
    }

    return true;
}

static inline bool extract_coord(md_xyz_coordinate_t* coord, str_t line) {
    ASSERT(coord);

    str_t original = line;

    str_t tokens[16];
    const int64_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);

    if (num_tokens < 4) {
        MD_LOG_ERROR("Invalid number of tokens in XYZ coordinate when parsing line: '%.*s', expected >= 4, got %i", (int)original.len, original.ptr, (int)num_tokens);
        return false;
    }

    int tok_idx = 0;
    if (num_tokens > 4) {
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
    
    if (tok_idx < num_tokens) {
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

static inline bool xyz_parse_model_header(md_xyz_model_t* model, md_buffered_reader_t* reader, uint32_t flags, int32_t* coord_count) {
    ASSERT(model);
    ASSERT(reader);
    ASSERT(coord_count);

    str_t line[2];
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
    const int64_t num_tok = extract_tokens(tokens, 2, &line[0]);
    if (num_tok) {
        count = (int32_t)parse_int(tokens[0]);
    }

    if (count <= 0) {
        MD_LOG_ERROR("Failed to extract coordinate count in XYZ header");  
        return false;
    }

    if (flags & XYZ_TINKER) {
        // Comment encoded after the first token
        comment = str_substr(line[0], (int64_t)(tokens[1].ptr - line[0].ptr), -1);
    } else {
        // Second line is the comment
        comment = line[1];
    }

    if (comment.len > 0) {
        str_copy_to_char_buf(model->comment, sizeof(model->comment), comment);
    }

    if (flags & XYZ_ARC) {
        const int64_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line[1]);
        
        if (num_tokens != 6) {
            MD_LOG_ERROR("Unexpected number of tokens (%i) when parsing XYZ Arc Cell data");
            return false;
        }
        model->cell_extent[0] = (float)parse_float(tokens[0]);
        model->cell_extent[1] = (float)parse_float(tokens[1]);
        model->cell_extent[2] = (float)parse_float(tokens[2]);

        model->cell_angle[0] = (float)parse_float(tokens[3]);
        model->cell_angle[1] = (float)parse_float(tokens[4]);
        model->cell_angle[2] = (float)parse_float(tokens[5]);
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
int64_t xyz_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    ASSERT(xyz);
    ASSERT(xyz->magic == MD_XYZ_TRAJ_MAGIC);

    if (!xyz->filesize) {
        MD_LOG_ERROR("File size is zero");
        return 0;
    }

    if (!xyz->file) {
        MD_LOG_ERROR("File handle is NULL");
        return 0;
    }

    if (!xyz->frame_offsets) {
        MD_LOG_ERROR("Frame offsets is empty");
        return 0;
    }

    if (!(0 <= frame_idx && frame_idx < (int64_t)md_array_size(xyz->frame_offsets) - 1)) {
        MD_LOG_ERROR("Frame index is out of range");
        return 0;
    }

    const int64_t beg = xyz->frame_offsets[frame_idx + 0];
    const int64_t end = xyz->frame_offsets[frame_idx + 1];
    const int64_t frame_size = end - beg;
    const int64_t total_size = sizeof(int64_t) + frame_size;

    if (frame_data_ptr) {
        // Store the index to the frame since this is generally not found within the actual frame data
        int64_t* ptr = (int64_t*)frame_data_ptr;
        ptr[0] = frame_idx;

        ASSERT(xyz->file);
        md_mutex_lock(&xyz->mutex);
        md_file_seek(xyz->file, beg, MD_FILE_BEG);
        const int64_t bytes_read = md_file_read(xyz->file, &ptr[1], frame_size);
        md_mutex_unlock(&xyz->mutex);
        ASSERT(frame_size == bytes_read);
    }

    return total_size;
}

bool xyz_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);
    ASSERT(frame_data_ptr);
    ASSERT(frame_data_size);

    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    if (xyz->magic != MD_XYZ_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame header, xyz magic did not match");
        return false;
    }

    const int64_t* ptr = frame_data_ptr;
    int32_t step = (int32_t)ptr[0];
    if (step < 0 || step >= xyz->header.num_frames) {
        MD_LOG_ERROR("Error when decoding frame data, corrupt frame index");
        return false;
    }

    str_t str = { .ptr = (const char*)(ptr+1), .len = frame_data_size };
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);

    md_xyz_model_t model = {0};
    int32_t coord_count = 0;
    if (!xyz_parse_model_header(&model, &reader, xyz->flags, &coord_count)) {
        MD_LOG_ERROR("Error when decoding header");
        return false;
    }

    int64_t i = 0;
    str_t line;
    str_t tokens[8];
    while (md_buffered_reader_extract_line(&line, &reader) && i < xyz->header.num_atoms) {
        if (line.len < 6) continue;

        const int64_t num_tokens = extract_tokens(tokens, ARRAY_SIZE(tokens), &line);

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
        header->unit_cell = md_util_unit_cell_from_extent_and_angles(model.cell_extent[0], model.cell_extent[1], model.cell_extent[2], model.cell_angle[0], model.cell_angle[1], model.cell_angle[2]);
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

    int32_t expected_count = 0;
    md_xyz_model_t mdl = {0};
    int64_t byte_offset = 0;

    while (xyz_parse_model_header(&mdl, reader, flags, &expected_count)) {
        md_array_ensure(data->coordinates, md_array_size(data->coordinates) + expected_count, alloc);
        
        mdl.byte_offset = byte_offset;
        mdl.beg_coord_index = (int32_t)md_array_size(data->coordinates);
        for (int32_t i = 0; i < expected_count; ++i) {
            str_t line;
            md_xyz_coordinate_t coord = {0};
            if (md_buffered_reader_extract_line(&line, reader) &&
                extract_coord(&coord, line))
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
        char* buf = md_alloc(md_heap_allocator, cap);
        ASSERT(buf);
        
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        result = xyz_parse(data, &reader, alloc, false);
        
        md_free(md_heap_allocator, buf, cap);
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

    const int64_t num_atoms = end_coord_index - beg_coord_index;

    md_array_ensure(mol->atom.x, num_atoms, alloc);
    md_array_ensure(mol->atom.y, num_atoms, alloc);
    md_array_ensure(mol->atom.z, num_atoms, alloc);
    md_array_ensure(mol->atom.element, num_atoms, alloc);

    for (int64_t i = beg_coord_index; i < end_coord_index; ++i) {
        float x = data->coordinates[i].x;
        float y = data->coordinates[i].y;
        float z = data->coordinates[i].z;
        str_t atom_name = {data->coordinates[i].element_symbol, sizeof(data->coordinates[i].element_symbol)};
        md_element_t element = (md_element_t)data->coordinates[i].atomic_number;

        mol->atom.count += 1;
        md_array_push(mol->atom.x, x, alloc);
        md_array_push(mol->atom.y, y, alloc);
        md_array_push(mol->atom.z, z, alloc);
        md_array_push(mol->atom.element, element, alloc);
        md_array_push(mol->atom.name, make_label(atom_name), alloc);
        md_array_push(mol->atom.flags, 0, alloc);
    }

    return true;
}

static bool xyz_init_from_str(md_molecule_t* mol, str_t str, md_allocator_i* alloc) {
    md_xyz_data_t data = {0};
    
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);
    bool result = xyz_parse(&data, &reader, md_heap_allocator, true);
    result = result && md_xyz_molecule_init(mol, &data, alloc);
    md_xyz_data_free(&data, md_heap_allocator);
    
    return result;
}

static bool xyz_init_from_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc) {
    md_xyz_data_t data = {0};
    
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int64_t cap = MEGABYTES(1);
        char* buf = md_alloc(md_heap_allocator, cap);
        ASSERT(buf);

        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        result = xyz_parse(&data, &reader, md_heap_allocator, true);
        result = result && md_xyz_molecule_init(mol, &data, alloc);

        md_xyz_data_free(&data, md_heap_allocator);
        md_free(md_heap_allocator, buf, cap);
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

    // Should this be exposed?
    md_allocator_i* alloc = md_temp_allocator;

    bool result = true;
    const int64_t frame_size = xyz_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        // This is a borderline case if one should use the md_temp_allocator as the raw frame size could potentially be several megabytes...
        void* frame_data = md_alloc(alloc, frame_size);
        const int64_t read_size = xyz_fetch_frame_data(inst, frame_idx, frame_data);
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

static bool try_read_cache(str_t cache_file, int64_t** offsets, int64_t* num_atoms, md_allocator_i* alloc) {
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int64_t read_bytes = 0;
        int64_t num_offsets = 0;
        uint64_t magic = 0;

        bool result = true;

        read_bytes = md_file_read(file, &magic, sizeof(uint64_t));
        if (read_bytes != sizeof(uint64_t) || magic != MD_XYZ_CACHE_MAGIC) {
            MD_LOG_ERROR("Failed to read offset cache, magic was incorrect or corrupt");
            result = false;
            goto done;
        }

        read_bytes = md_file_read(file, num_atoms, sizeof(int64_t));
        if (read_bytes != sizeof(int64_t) || num_atoms == 0) {
            MD_LOG_ERROR("Failed to read offset cache, number of atoms was zero or corrupt");
            result = false;
            goto done;
        }

        read_bytes = md_file_read(file, &num_offsets, sizeof(int64_t));
        if (read_bytes != sizeof(int64_t) || num_offsets == 0) {
            MD_LOG_ERROR("Failed to read offset cache, number of frames was zero or corrupted");
            result = false;
            goto done;
        }

        int64_t* tmp_offsets = 0;
        md_array_resize(tmp_offsets, num_offsets, alloc);

        read_bytes = md_file_read(file, tmp_offsets, num_offsets * sizeof(int64_t));
        if (read_bytes != (int64_t)(num_offsets * sizeof(int64_t))) {
            MD_LOG_ERROR("Failed to read offset cache, offsets are incomplete");
            md_array_free(tmp_offsets, alloc);
            result = false;
            goto done;
        }

        *offsets = tmp_offsets;
    done:
        md_file_close(file);
        return result;
    }

    return false;
}

static bool write_cache(str_t cache_file, int64_t* offsets, int64_t num_atoms) {
    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (file) {
        int64_t num_offsets = md_array_size(offsets);
        int64_t written_bytes = 0;
        const uint64_t magic = MD_XYZ_CACHE_MAGIC;

        written_bytes = md_file_write(file, &magic, sizeof(uint64_t));
        ASSERT(written_bytes == sizeof(uint64_t));

        written_bytes = md_file_write(file, &num_atoms, sizeof(int64_t));
        ASSERT(written_bytes == sizeof(int64_t));

        written_bytes = md_file_write(file, &num_offsets, sizeof(int64_t));
        ASSERT(written_bytes == sizeof(int64_t));

        written_bytes = md_file_write(file, offsets, num_offsets * sizeof(int64_t));
        ASSERT(written_bytes == (int64_t)(num_offsets * sizeof(int64_t)));

        md_file_close(file);
        return true;
    }

    MD_LOG_ERROR("Failed to write offset cache, could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
    return false;
}

void xyz_trajectory_free(struct md_trajectory_o* inst) {
    ASSERT(inst);
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    if (xyz->file) md_file_close(xyz->file);
    if (xyz->frame_offsets) md_array_free(xyz->frame_offsets, xyz->allocator);
    md_mutex_destroy(&xyz->mutex);
}

md_trajectory_i* md_xyz_trajectory_create(str_t filename, struct md_allocator_i* alloc) {
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_ERROR("Failed to open file for XYZ trajectory");
        return false;
    }

    uint32_t flags = 0;
    {
        char buf[1024];
        int64_t len = md_file_read(file, buf, sizeof(buf));
        md_buffered_reader_t reader = md_buffered_reader_from_str((str_t){buf, len});
        if (!extract_flags(&flags, &reader)) {
            MD_LOG_ERROR("Failed to determine format for XYZ trajectory");
            return false;
        }
    }

    int64_t filesize = md_file_size(file);
    md_file_close(file);

    char buf[1024] = "";
    int len = snprintf(buf, sizeof(buf), "%.*s.cache", (int)filename.len, filename.ptr);
    str_t cache_file = {buf, len};

    int64_t num_atoms = 0;
    int64_t* offsets = 0;

    if (!try_read_cache(cache_file, &offsets, &num_atoms, alloc)) {
        md_xyz_data_t data = {0};
        if (!md_xyz_data_parse_file(&data, filename, md_heap_allocator)) {
            return false;
        }

        if (data.num_models <= 1) {
            md_log(MD_LOG_TYPE_INFO, "The xyz file did not contain multiple entries and cannot be read as a trajectory");
            md_xyz_data_free(&data, md_heap_allocator);
            return false;
        }

        {
            // Validate the models
            const int64_t ref_length = data.models[0].end_coord_index - data.models[0].beg_coord_index;
            for (int64_t i = 1; i < data.num_models; ++i) {
                const int64_t length = data.models[i].end_coord_index - data.models[i].beg_coord_index;
                if (length != ref_length) {
                    MD_LOG_ERROR("The xyz file models are not of equal length and cannot be read as a trajectory");
                    md_xyz_data_free(&data, md_heap_allocator);
                    return false;
                }
            }
            num_atoms = ref_length;
        }

        for (int64_t i = 0; i < data.num_models; ++i) {
            md_array_push(offsets, data.models[i].byte_offset, alloc);
        }
        md_array_push(offsets, filesize, alloc);

        write_cache(cache_file, offsets, num_atoms);
    }

    int64_t max_frame_size = 0;
    for (int64_t i = 0; i < md_array_size(offsets) - 1; ++i) {
        const int64_t beg = offsets[i + 0];
        const int64_t end = offsets[i + 1];
        const int64_t frame_size = end - beg;
        max_frame_size = MAX(max_frame_size, frame_size);
    }

    void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));
    ASSERT(mem);
    MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));

    md_trajectory_i* traj = mem;
    xyz_trajectory_t* xyz = (xyz_trajectory_t*)(traj + 1);

    xyz->magic = MD_XYZ_TRAJ_MAGIC;
    xyz->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    xyz->filesize = filesize;
    xyz->frame_offsets = offsets;
    xyz->allocator = alloc;
    xyz->mutex = md_mutex_create();
    xyz->header = (md_trajectory_header_t) {
        .num_frames = md_array_size(offsets) - 1,
        .num_atoms = num_atoms,
        .max_frame_data_size = max_frame_size,
        .time_unit = {0},
    };
    xyz->flags = flags;

    traj->inst = (struct md_trajectory_o*)xyz;
    traj->get_header = xyz_get_header;
    traj->load_frame = xyz_load_frame;
    traj->fetch_frame_data = xyz_fetch_frame_data;
    traj->decode_frame_data = xyz_decode_frame_data;

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
    
    md_allocator_i* alloc = xyz->allocator;
    xyz_trajectory_free(traj->inst);
    md_free(alloc, traj, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));
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
