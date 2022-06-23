#include "md_xyz.h"

#include <core/md_common.h>
#include <core/md_array.inl>
#include <core/md_str.h>
#include <core/md_file.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_sync.h>
#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_util.h>

#include <string.h> // memcpy
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_XYZ_MOL_MAGIC  0x285ada29078a9bc8
#define MD_XYZ_TRAJ_MAGIC 0x2312ad7b78a9bc78

enum {
    XYZ_ELEMENT_SYMBOL  = 1,
    XYZ_ATOMIC_NUMBER   = 2,
    XYZ_TINKER          = 4,
    XYZ_TINKER_ARC      = 8,
};

typedef struct xyz_span_t {
    int beg;
    int end;
} xyz_span_t;

typedef struct xyz_format_t {
    uint32_t flags;
    xyz_span_t element;
    xyz_span_t x;
    xyz_span_t y;
    xyz_span_t z;
    xyz_span_t atom_index;
    xyz_span_t atom_type;

    // Header info
    xyz_span_t coord_count;
    xyz_span_t cell_ext[3];
    xyz_span_t cell_angle[3];
} xyz_format_t;

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
    xyz_format_t format;
} xyz_trajectory_t;

// We massage the beg and end indices here to correspond the xyz specification
// This makes our life easier when specifying all the different ranges
static inline int32_t extract_int(str_t line, int64_t beg, int64_t end) {
    if (line.len < end) return 0;
    return (int32_t)parse_int(trim_whitespace(substr(line, beg, end-beg)));
}

static inline float extract_float(str_t line, int64_t beg, int64_t end) {
    if (line.len < end) return 0.0f;
    return (float)parse_float(trim_whitespace(substr(line, beg, end-beg)));
}

static inline bool is_float(str_t str) {
    str = trim_whitespace(str);
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    if (*c == '-') ++c;
    while (c < end && is_digit(*c)) ++c;
    if (*c != '.') return false;
    ++c;
    while (c < end && is_digit(*c)) ++c;
    return c == end;
}

static inline bool is_unsigned_int(str_t str) {
    str = trim_whitespace(str);
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
    str = trim_whitespace(str);

    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    if (c >= end) return false;
    while (c < end) {
        if (!is_alpha(*c)) return false;
        ++c;
    }
    return true;
}

static inline bool extract_format(xyz_format_t* format, str_t str) {
    ASSERT(format);
    if (str_empty(str)) return false;

    // Extract first three lines
    str_t lines[3];
    extract_line(&lines[0], &str);
    extract_line(&lines[1], &str);
    extract_line(&lines[2], &str);

#define BEG(x) (int32_t)(x.ptr - base)
#define END(x) (int32_t)(x.ptr - base + x.len)

    {
        // Test if first line has an unsigned integer as first token (Should be universally applicable to all XYZ + ARC formats
        str_t token;
        str_t line = lines[0];
        if (!extract_next_token(&token, &line) || !is_unsigned_int(token)) {
            md_print(MD_LOG_TYPE_ERROR, "Invalid format for XYZ: Missing coordinate count");
            return false;
        }

        const char* base = lines[0].ptr;
        format->coord_count.beg = 0;
        format->coord_count.end = END(token);
    }

    {
        // Test if we have an ARC trajectory
        // Second line should then contain exactly 6 floats
        int token_count = 0;
        str_t token[16];
        str_t line = lines[1];
        while (token_count < ARRAY_SIZE(token) && extract_next_token(&token[token_count], &line)) {
            if (!is_float(token[token_count])) {
                break;
            }
            token_count += 1;
        }
        if (token_count == 6) {
            const char* base = lines[1].ptr;
            format->flags |= XYZ_TINKER_ARC;
            format->cell_ext[0].beg = END(token[0]) - 11;
            format->cell_ext[0].end = END(token[0]);
            format->cell_ext[1].beg = END(token[1]) - 11;
            format->cell_ext[1].end = END(token[1]);
            format->cell_ext[2].beg = END(token[2]) - 11;
            format->cell_ext[2].end = END(token[2]);

            format->cell_angle[0].beg = END(token[3]) - 11;
            format->cell_angle[0].end = END(token[3]);
            format->cell_angle[1].beg = END(token[4]) - 11;
            format->cell_angle[1].end = END(token[4]);
            format->cell_angle[2].beg = END(token[5]) - 11;
            format->cell_angle[2].end = END(token[5]);
        }
    }

    // Determine coordinate structure
    // Traditional XYZ only holds 4 fields (atomic number or element symbol and coordinates)
    // Tinker holds additional fields (atom index, atom type and connectivity) at least 6 fields
    int token_count = 0;
    str_t token[16];
    str_t line = lines[2];
    while (token_count < ARRAY_SIZE(token) && extract_next_token(&token[token_count], &line)) {
        token_count += 1;
    }
    if (token_count >= 6){
        if (is_unsigned_int(token[0]) &&
            (is_string(token[1]) || is_unsigned_int(token[1])) &&
            is_float(token[2]) && is_float(token[3]) && is_float(token[4]) &&
            is_unsigned_int(token[5])) {
            format->flags |= XYZ_TINKER;

            const char* base = lines[2].ptr;

            // Determine tinker spans
            format->atom_index.beg = 0;
            format->atom_index.end = BEG(token[1]);

            format->element.beg = END(token[0]);
            format->element.end = END(token[2]) - 11;

            if (is_string(token[1])) {
                format->flags |= XYZ_ELEMENT_SYMBOL;
            } else if (is_unsigned_int(token[1])) {
                format->flags |= XYZ_ATOMIC_NUMBER;
            } else {
                md_print(MD_LOG_TYPE_ERROR, "Invalid format for element");
                return false;
            }

            format->x.beg = END(token[2]) - 11;
            format->x.end = END(token[2]);

            format->y.beg = END(token[3]) - 11;
            format->y.end = END(token[3]);

            format->z.beg = END(token[4]) - 11;
            format->z.end = END(token[4]);

            format->atom_type.beg = END(token[5]) - 6;
            format->atom_type.end = END(token[5]);

        } else {
            md_print(MD_LOG_TYPE_ERROR, "Invalid format for XYZ: Unknown variation of Tinker format?");
            return false;
        }
    }
    else if (token_count == 4) {
        if ((is_string(token[0]) || is_unsigned_int(token[0])) &&
            is_float(token[1]) && is_float(token[2]) && is_float(token[3])) {
            
            const char* base = lines[2].ptr;

            // Determine standard xyz spans
            format->element.beg = 0;
            format->element.end = END(token[1]) - 11;

            if (is_string(token[0])) {
                format->flags |= XYZ_ELEMENT_SYMBOL;
            } else if (is_unsigned_int(token[0])) {
                format->flags |= XYZ_ATOMIC_NUMBER;
            } else {
                md_print(MD_LOG_TYPE_ERROR, "Invalid format for element");
                return false;
            }

            format->x.beg = END(token[1]) - 11;
            format->x.end = END(token[1]);

            format->y.beg = END(token[2]) - 11;
            format->y.end = END(token[2]);

            format->z.beg = END(token[3]) - 11;
            format->z.end = END(token[3]);
        } else {
            md_print(MD_LOG_TYPE_ERROR, "Invalid format for XYZ: Unknown variation of Tinker format?");
            return false;
        }
    } else {
        md_print(MD_LOG_TYPE_ERROR, "Invalid format for XYZ");
        return false;
    }

#undef BEG
#undef END

    return true;
}

static inline md_xyz_coordinate_t extract_coord(str_t line, const xyz_format_t* format) {
    ASSERT(format);
    md_xyz_coordinate_t coord = { 0 };
    
    if (format->flags & XYZ_ELEMENT_SYMBOL) {
        str_t element = trim_whitespace(substr(line, format->element.beg, format->element.end - format->element.beg));
        strncpy(coord.element_symbol, element.ptr, MIN(element.len, (int)sizeof(coord.element_symbol)));
    } else if (format->flags & XYZ_ATOMIC_NUMBER) {
        coord.atomic_number = extract_int(line, format->element.beg, format->element.end);
    } else {
        md_print(MD_LOG_TYPE_ERROR, "Undetermined format for element type in XYZ coordinate");
        goto done;
    }

    coord.x = extract_float(line, format->x.beg, format->x.end);
    coord.y = extract_float(line, format->y.beg, format->y.end);
    coord.z = extract_float(line, format->z.beg, format->z.end);

    if (format->flags & XYZ_TINKER) {
        coord.atom_type = extract_int(line, format->atom_type.beg, format->atom_type.end);
        // Connectivity follows atom_type
        str_t token;
        line = substr(line, format->atom_type.end, -1);
        int connection_count = 0;
        while (extract_next_token(&token, &line)) {
            token = trim_whitespace(token);
            if (token.len > 0) {
                int connection = (int)parse_int(token);
                if (connection) {
                    coord.connectivity[connection_count++] = connection;
                } else {
                    md_print(MD_LOG_TYPE_ERROR, "Invalid connectivity information in XYZ file");
                    goto done;
                }
            }
        }
    }

    done:
    return coord;
}

static inline bool xyz_parse_model_header(md_xyz_model_t* model, str_t* str, const xyz_format_t* format, int32_t* coord_count) {
    ASSERT(model);
    ASSERT(str);
    ASSERT(format);
    ASSERT(coord_count);

    if (str->len == 0) return false;

    str_t line[2];

    if (!extract_line(&line[0], str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to extract header line");
        return false;
    }

    int32_t count = extract_int(line[0], format->coord_count.beg, format->coord_count.end);
    if (count <= 0) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to parse coordinate count in header");
        return false;
    }

    if (!(format->flags & XYZ_TINKER) || (format->flags & XYZ_TINKER_ARC)) {
        if (!extract_line(&line[1], str)) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to extract extra line");  
            return false;
        }
    }

    str_t comment = {0};
    if (format->flags & XYZ_TINKER) {
        // Comment is in the first line following the coordinate count
        comment = substr(line[0], format->coord_count.end, -1);
    } else {
        // Comment is the second line
        comment = line[1];
    }

    if (comment.len > 0) {
        strncpy(model->comment, comment.ptr, MIN(sizeof(model->comment)-1, (size_t)comment.len));
    }

    if (format->flags & XYZ_TINKER_ARC) {
        model->cell_extent[0] = extract_float(line[1], format->cell_ext[0].beg, format->cell_ext[0].end);
        model->cell_extent[1] = extract_float(line[1], format->cell_ext[1].beg, format->cell_ext[1].end);
        model->cell_extent[2] = extract_float(line[1], format->cell_ext[2].beg, format->cell_ext[2].end);

        model->cell_angle[0] = extract_float(line[1], format->cell_angle[0].beg, format->cell_angle[0].end);
        model->cell_angle[1] = extract_float(line[1], format->cell_angle[1].beg, format->cell_angle[1].end);
        model->cell_angle[2] = extract_float(line[1], format->cell_angle[2].beg, format->cell_angle[2].end);
    }

    *coord_count = count;

    return true;
}

static inline int32_t xyz_parse_model_coordinates(md_xyz_data_t* data, str_t* str, const xyz_format_t* format, int32_t count, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(str);
    ASSERT(format);
    ASSERT(count);
    ASSERT(alloc);

    int32_t i = 0;
    str_t line;
    for (; i < count; ++i) {
        if (!extract_line(&line, str))
            break;
        md_array_push(data->coordinates, extract_coord(line, format), alloc);
    }

    return i;
}

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
        md_print(MD_LOG_TYPE_ERROR, "File size is zero");
        return 0;
    }

    if (!xyz->file) {
        md_print(MD_LOG_TYPE_ERROR, "File handle is NULL");
        return 0;
    }

    if (!xyz->frame_offsets) {
        md_print(MD_LOG_TYPE_ERROR, "Frame offsets is empty");
        return 0;
    }

    if (!(0 <= frame_idx && frame_idx < (int64_t)md_array_size(xyz->frame_offsets) - 1)) {
        md_print(MD_LOG_TYPE_ERROR, "Frame index is out of range");
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
        md_print(MD_LOG_TYPE_ERROR, "Error when decoding frame header, xyz magic did not match");
        return false;
    }

    const int64_t* ptr = frame_data_ptr;
    int32_t step = (int32_t)ptr[0];
    if (step < 0 || step >= xyz->header.num_frames) {
        md_print(MD_LOG_TYPE_ERROR, "Error when decoding frame data, corrupt frame index");
        return false;
    }

    str_t str = { .ptr = (const char*)(ptr+1), .len = frame_data_size };

    md_xyz_model_t model = {0};
    int32_t coord_count = 0;
    if (!xyz_parse_model_header(&model, &str, &xyz->format, &coord_count)) {
        md_print(MD_LOG_TYPE_ERROR, "Error when decoding header");
        return false;
    }


    int64_t i = 0;
    str_t line;
    while (extract_line(&line, &str) && i < xyz->header.num_atoms) {
        if (line.len < 6) continue;

        float coord_x = extract_float(line, xyz->format.x.beg, xyz->format.x.end);
        float coord_y = extract_float(line, xyz->format.y.beg, xyz->format.y.end);
        float coord_z = extract_float(line, xyz->format.z.beg, xyz->format.z.end);
        if (x) x[i] = coord_x;
        if (y) y[i] = coord_y;
        if (z) z[i] = coord_z;

        i += 1;
    }

    if (header) {
        header->num_atoms = i;
        header->index = step;
        header->timestamp = (double)(step); // This information is missing from xyz trajectories
        mat3_t box = {0};
        if (model.cell_extent[0] != 0) {
            if (model.cell_angle[0] != 90 && model.cell_angle[1] != 90 && model.cell_angle[2] != 90) {
                box = md_util_compute_unit_cell_basis(model.cell_extent[0], model.cell_extent[1], model.cell_extent[2], model.cell_angle[0], model.cell_angle[1], model.cell_angle[2]);
            } else {
                box.elem[0][0] = model.cell_extent[0];
                box.elem[1][1] = model.cell_extent[1];
                box.elem[2][2] = model.cell_extent[2];
            }
        }
        memcpy(header->box, &box, sizeof(header->box));
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

bool md_xyz_data_parse_str(md_xyz_data_t* data, str_t str, struct md_allocator_i* alloc) {
    xyz_format_t format = {0};
    if (extract_format(&format, str)) {
        int32_t expected_count = 0;
        md_xyz_model_t model = {0};
        const char* base_offset = str.ptr;
        int64_t byte_offset = 0;
        while (xyz_parse_model_header(&model, &str, &format, &expected_count)) {
            model.beg_coord_index = (int32_t)md_array_size(data->coordinates);
            model.end_coord_index = (int32_t)md_array_size(data->coordinates);
            model.byte_offset = byte_offset;

            md_xyz_model_t* last_model = md_array_last(data->models);
            if (last_model) {
                last_model->end_coord_index = (int32_t)md_array_size(data->coordinates);
            }
            int32_t real_count = xyz_parse_model_coordinates(data, &str, &format, expected_count, alloc);
            if (real_count != expected_count) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to read expected coordinate count in XYZ");
                return false;
            }
            byte_offset = str.ptr - base_offset;
        }
    }

    md_xyz_model_t* last_model = md_array_last(data->models);
    if (last_model) {
        last_model->end_coord_index = (int32_t)md_array_size(data->coordinates);
    }
    data->num_coordinates = md_array_size(data->coordinates);
    data->num_models = md_array_size(data->models);

    return true;
}

bool md_xyz_data_parse_file(md_xyz_data_t* data, str_t filename, struct md_allocator_i* alloc) {
    buffered_reader_t reader = {0};
    bool result = false;

    if (buffered_reader_init(&reader, KILOBYTES(32), filename, alloc)) {
        str_t chunk;
        int64_t reader_offset = buffered_reader_get_offset(&reader);
        if (buffered_reader_get_chunk(&chunk, &reader)) {
            const char* chunk_base = chunk.ptr;
            xyz_format_t format = {0};
            if (extract_format(&format, chunk)) {
                int32_t expected_count = 0;

                while (true) {
                    if (expected_count == 0) {
                        const int64_t chunk_offset = chunk.ptr - chunk_base;
                        md_xyz_model_t model = {0};
                        if (!xyz_parse_model_header(&model, &chunk, &format, &expected_count)) {
                            md_print(MD_LOG_TYPE_ERROR, "Failed to read xyz header");
                            goto done;
                        }
                        model.beg_coord_index = (int32_t)md_array_size(data->coordinates);
                        model.end_coord_index = (int32_t)md_array_size(data->coordinates);
                        model.byte_offset = reader_offset + chunk_offset;

                        md_xyz_model_t* last_model = md_array_last(data->models);
                        if (last_model) {
                            last_model->end_coord_index = (int32_t)md_array_size(data->coordinates);
                        }

                        md_array_push(data->models, model, alloc);
                    }
                    int32_t read_count = xyz_parse_model_coordinates(data, &chunk, &format, expected_count, alloc);
                    expected_count -= read_count;

                    if (chunk.len == 0) {
                        reader_offset = buffered_reader_get_offset(&reader);
                        if (!buffered_reader_get_chunk(&chunk, &reader)) {
                            break;
                        }
                        chunk_base = chunk.ptr;
                    }
                }
            }
        }

        result = true;
        done:
        buffered_reader_free(&reader, alloc);
    }

    md_xyz_model_t* last_model = md_array_last(data->models);
    if (last_model) {
        last_model->end_coord_index = (int32_t)md_array_size(data->coordinates);
    }
    data->num_coordinates = md_array_size(data->coordinates);
    data->num_models = md_array_size(data->models);

    return result;
}

void md_xyz_data_free(md_xyz_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    if (data->coordinates) md_array_free(data->coordinates, alloc);
    if (data->models) md_array_free(data->models, alloc);
    memset(data, 0, sizeof(md_xyz_data_t));
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

    if (mol->inst) {
        md_print(MD_LOG_TYPE_DEBUG, "molecule inst object is not zero, potentially leaking memory when clearing");
    }

    memset(mol, 0, sizeof(md_molecule_t));

    mol->inst = md_alloc(alloc, sizeof(xyz_molecule_t));
    xyz_molecule_t* inst = (xyz_molecule_t*)mol->inst;
    memset(inst, 0, sizeof(xyz_molecule_t));
    
    inst->magic = MD_XYZ_MOL_MAGIC;
    inst->allocator = md_arena_allocator_create(alloc, KILOBYTES(64));

    const int64_t num_atoms = end_coord_index - beg_coord_index;

    // Change allocator to arena so we can very easily free it later.
    alloc = inst->allocator;

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

    md_util_postprocess_molecule(mol, alloc);

    return true;
}

bool md_xyz_molecule_free(md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(mol->inst);
    ASSERT(alloc);

    xyz_molecule_t* inst = (xyz_molecule_t*)mol->inst;
    if (inst->magic != MD_XYZ_MOL_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "XYZ magic did not match!");
        return false;
    }

    md_arena_allocator_destroy(inst->allocator);
    md_free(alloc, inst, sizeof(xyz_molecule_t));
    memset(mol, 0, sizeof(md_molecule_t));

    return true;
}

static bool xyz_init_from_str(md_molecule_t* mol, str_t str, md_allocator_i* alloc) {
    md_xyz_data_t data = {0};

    md_xyz_data_parse_str(&data, str, default_allocator);
    bool success = md_xyz_molecule_init(mol, &data, alloc);
    md_xyz_data_free(&data, default_allocator);
    
    return success;
}

static bool xyz_init_from_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc) {
    buffered_reader_t reader = {0};
    md_allocator_i* temp_alloc = default_allocator;
    bool result = false;

    if (buffered_reader_init(&reader, KILOBYTES(64), filename, temp_alloc)) {
        md_xyz_data_t data = {0};
        str_t chunk;
        if (buffered_reader_get_chunk(&chunk, &reader)) {
            xyz_format_t format = {0};
            if (extract_format(&format, chunk)) {
                int32_t expected_count = 0;
                md_xyz_model_t model = {0};
                if (!xyz_parse_model_header(&model, &chunk, &format, &expected_count)) {
                    md_print(MD_LOG_TYPE_ERROR, "Failed to read xyz header");
                    goto done;
                }

                model.beg_coord_index = 0;
                model.end_coord_index = 0;
                md_array_push(data.models, model, temp_alloc);

                while (true) {
                    int32_t read_count = xyz_parse_model_coordinates(&data, &chunk, &format, expected_count, temp_alloc);
                    expected_count -= read_count;
                    if (expected_count == 0) break;
                    
                    if (!buffered_reader_get_chunk(&chunk, &reader)) {
                        goto done;
                    }
                }
            }
        }

        md_xyz_model_t* last_model = md_array_last(data.models);
        if (last_model) {
            last_model->end_coord_index = (int32_t)md_array_size(data.coordinates);
        }
        data.num_coordinates = md_array_size(data.coordinates);
        data.num_models = md_array_size(data.models);
        result = md_xyz_molecule_init(mol, &data, alloc);
        done:
        buffered_reader_free(&reader, temp_alloc);
    }
    return result;
}

static md_molecule_api xyz_molecule_api = {
    xyz_init_from_str,
    xyz_init_from_file,
    md_xyz_molecule_free
};

md_molecule_api* md_xyz_molecule_api() {
    return &xyz_molecule_api;
}

bool xyz_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    xyz_trajectory_t* xyz = (xyz_trajectory_t*)inst;
    if (xyz->magic != MD_XYZ_TRAJ_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    // Should this be exposed?
    md_allocator_i* alloc = default_temp_allocator;

    bool result = true;
    const int64_t frame_size = xyz_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        // This is a borderline case if one should use the default_temp_allocator as the raw frame size could potentially be several megabytes...
        void* frame_data = md_alloc(alloc, frame_size);
        const int64_t read_size = xyz_fetch_frame_data(inst, frame_idx, frame_data);
        if (read_size != frame_size) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read the expected size");
            md_free(alloc, frame_data, frame_size);
            return false;
        }

        result = xyz_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
        md_free(alloc, frame_data, frame_size);
    }

    return result;
}

#define CACHE_MAGIC 0x8265485749172bab
static bool try_read_cache(str_t cache_file, int64_t** offsets, int64_t* num_atoms, md_allocator_i* alloc) {
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int64_t read_bytes = 0;
        int64_t num_offsets = 0;
        uint64_t magic = 0;

        bool result = true;

        read_bytes = md_file_read(file, &magic, sizeof(uint64_t));
        if (read_bytes != sizeof(uint64_t) || magic != CACHE_MAGIC) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, magic was incorrect or corrupt");
            result = false;
            goto done;
        }

        read_bytes = md_file_read(file, num_atoms, sizeof(int64_t));
        if (read_bytes != sizeof(int64_t) || num_atoms == 0) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, number of atoms was zero or corrupt");
            result = false;
            goto done;
        }

        read_bytes = md_file_read(file, &num_offsets, sizeof(int64_t));
        if (read_bytes != sizeof(int64_t) || num_offsets == 0) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, number of frames was zero or corrupted");
            result = false;
            goto done;
        }

        int64_t* tmp_offsets = 0;
        md_array_resize(tmp_offsets, num_offsets, alloc);

        read_bytes = md_file_read(file, tmp_offsets, num_offsets * sizeof(int64_t));
        if (read_bytes != (int64_t)(num_offsets * sizeof(int64_t))) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, offsets are incomplete");
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
        const uint64_t magic = CACHE_MAGIC;

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

    md_printf(MD_LOG_TYPE_ERROR, "Failed to write offset cache, could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
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
        md_print(MD_LOG_TYPE_ERROR, "Failed to open file for XYZ trajectory");
        return false;
    }

    xyz_format_t format = {0};
    {
        char buf[1024];
        int64_t len = md_file_read(file, buf, sizeof(buf));
        if (!extract_format(&format, (str_t){buf, len})) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to interpret format for XYZ trajectory");
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
        if (!md_xyz_data_parse_file(&data, filename, default_allocator)) {
            return false;
        }

        if (data.num_models <= 1) {
            md_print(MD_LOG_TYPE_INFO, "The xyz file did not contain multiple entries and cannot be read as a trajectory");
            md_xyz_data_free(&data, default_allocator);
            return false;
        }

        {
            // Validate the models
            const int64_t ref_length = data.models[0].end_coord_index - data.models[0].beg_coord_index;
            for (int64_t i = 1; i < data.num_models; ++i) {
                const int64_t length = data.models[i].end_coord_index - data.models[i].beg_coord_index;
                if (length != ref_length) {
                    md_print(MD_LOG_TYPE_ERROR, "The xyz file models are not of equal length and cannot be read as a trajectory");
                    md_xyz_data_free(&data, default_allocator);
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
    memset(mem, 0, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));

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
    xyz->format = format;

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
        md_printf(MD_LOG_TYPE_ERROR, "Trajectory is not a valid xyz trajectory.");
        ASSERT(false);
        return;
    }
    
    md_allocator_i* alloc = xyz->allocator;
    xyz_trajectory_free(traj->inst);
    md_free(alloc, traj, sizeof(md_trajectory_i) + sizeof(xyz_trajectory_t));
}

static md_trajectory_api xyz_traj_api = {
    md_xyz_trajectory_create,
    md_xyz_trajectory_free,
};

md_trajectory_api* md_xyz_trajectory_api() {
    return &xyz_traj_api;
}

#ifdef __cplusplus
}
#endif
