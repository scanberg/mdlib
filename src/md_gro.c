#include "md_gro.h"

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_array.inl>
#include <md_util.h>

#include <string.h>

#define MD_GRO_MOL_MAGIC  0xbabcebf677bfaf67

typedef struct label {
    char str[8];
} label_t;

typedef struct gro_molecule {
    uint64_t magic;
    label_t* labels;
    struct md_allocator_i* allocator;
} gro_molecule_t;

static inline bool ranges_overlap(md_range_t a, md_range_t b) {
    return (a.beg < b.end && b.beg < a.end);
}

static inline int64_t compute_position_field_width(str_t line) {
    // First float starts at offset 20, count
    if (line.len > 20) {
        const char* c = line.ptr + 20;
        const char* end = line.ptr + line.len;
        while (c != end && *c != '\n' && *c == ' ') ++c;
        while (c != end && *c != '\n' && *c != ' ') ++c;
        if (c != end) {
            return (int64_t)(c - (line.ptr + 20));
        }
    }
    return 0;
}

static inline bool parse_header(str_t* str, md_gro_data_t* data) {
    str_t line = {0};
    if (!extract_line(&line, str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read title");
        return false;
    }

    line = trim_whitespace(line);
    const int64_t len = MIN((int64_t)ARRAY_SIZE(data->title) - 1, line.len);
    strncpy(data->title, line.ptr, len);
    data->title[len] = '\0';

    if (!extract_line(&line, str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read number of atoms");
        return false;
    }

    data->num_atoms = parse_int(trim_whitespace(line));
    const int64_t min_bounds = 1;
    const int64_t max_bounds = 10000000;
    if (data->num_atoms < min_bounds || max_bounds < data->num_atoms) {
        md_printf(MD_LOG_TYPE_ERROR, "Number of atoms (%lli) outside of plausible range [%lli, %lli]", data->num_atoms, min_bounds, max_bounds);
        return false;
    }

    return true;
}

static inline str_t parse_atom_data(str_t str, md_gro_data_t* data, int64_t pos_field_width, int64_t* read_count, int64_t read_target, md_allocator_i* alloc) {
    str_t line = {0};
    while ((*read_count < read_target) && extract_line(&line, &str)) {
        int64_t res_id  = parse_int(trim_whitespace(substr(line, 0, 5)));
        str_t res_name  = trim_whitespace(substr(line, 5, 5));
        str_t atom_name = trim_whitespace(substr(line, 10, 5));
        double x = parse_float(trim_whitespace(substr(line, 20 + 0 * pos_field_width, pos_field_width))) * 10.0; // nm -> �
        double y = parse_float(trim_whitespace(substr(line, 20 + 1 * pos_field_width, pos_field_width))) * 10.0;
        double z = parse_float(trim_whitespace(substr(line, 20 + 2 * pos_field_width, pos_field_width))) * 10.0;

        md_gro_atom_t atom = {0};
        atom.res_id = (int32_t)res_id;
        strncpy(atom.res_name, res_name.ptr, MIN(res_name.len, (int64_t)ARRAY_SIZE(atom.res_name) - 1));
        strncpy(atom.atom_name, atom_name.ptr, MIN(atom_name.len, (int64_t)ARRAY_SIZE(atom.atom_name) - 1));
        atom.x = (float)x;
        atom.y = (float)y;
        atom.z = (float)z;
        md_array_push(data->atom_data, atom, alloc);
        *read_count += 1;
    }

    return str;
}

static inline bool parse_unitcell(str_t str, md_gro_data_t* data) {
    ASSERT(data);

    str_t line = {0,0};
    if (!extract_line(&line, &str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to extract line for unit cell.");
        return false;
    }

    data->cell_ext[0] = (float)parse_float(trim_whitespace(substr(line, 0, 10)));
    data->cell_ext[1] = (float)parse_float(trim_whitespace(substr(line, 10, 10)));
    data->cell_ext[2] = (float)parse_float(trim_whitespace(substr(line, 20, 10)));

    return true;
}

bool md_gro_data_parse_str(md_gro_data_t* data, str_t str, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    
    if (!parse_header(&str, data)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to parse gro header!");
        return false;
    }

    const int64_t pos_field_width = compute_position_field_width(str);

    int64_t num_atoms_read = 0;
    str = parse_atom_data(str, data, pos_field_width, &num_atoms_read, data->num_atoms, alloc);

    if (num_atoms_read != data->num_atoms) {
        md_printf(MD_LOG_TYPE_ERROR, "There was a descrepancy between the number of atoms read (%lli) compared to the number of atoms specified in the header (%lli)", num_atoms_read, data->num_atoms);
        return false;
    }

    if (!parse_unitcell(str, data)) {
        return false;
    }

    return true;
}

bool md_gro_data_parse_file(md_gro_data_t* data, str_t filename, struct md_allocator_i* alloc) {
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        bool success = false;
        const int64_t buf_size = KILOBYTES(64ULL);
        char* buf = md_alloc(default_temp_allocator, buf_size);
        
        int64_t pos_field_width = 0;
        int64_t num_atoms_read = 0;

        // Read header explicitly
        {
            const int64_t header_size = KILOBYTES(1);
            str_t str = {.ptr = buf, .len = md_file_read(file, buf, header_size) };

            if (!parse_header(&str, data)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to parse gro header!");
                goto done;
            }

            pos_field_width = compute_position_field_width(str);
            if (!pos_field_width) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to determine position field width");
                goto done;
            }

            // Once we have read the header, we set the file position to just after the header
            const int64_t seek_target = str.ptr - buf;
            if (!md_file_seek(file, seek_target, MD_FILE_BEG)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to seek in file");
                goto done;
            }
        }

        int64_t bytes_read = 0;
        int64_t buf_offset = 0;
        str_t str = {.ptr = buf};

        // Read fields in chunks
        while ((bytes_read = md_file_read(file, buf + buf_offset, buf_size - buf_offset)) > 0) {
            str.len = bytes_read + buf_offset;

            // Have we read all bytes in file?
            if (bytes_read < buf_size - buf_offset) {
                str = parse_atom_data(str, data, pos_field_width, &num_atoms_read, data->num_atoms, alloc);
                break;
            } else {
                // We want to make sure we send complete lines for parsing so we locate the last '\n'
                const int64_t last_new_line = rfind_char(str, '\n');
                if (last_new_line == -1) {
                    md_print(MD_LOG_TYPE_ERROR, "Unexpected structure within gro file, missing new lines?");
                    goto done;
                }
                ASSERT(str.ptr[last_new_line] == '\n');
                str.len = last_new_line + 1;
                parse_atom_data(str, data, pos_field_width, &num_atoms_read, data->num_atoms, alloc);
            }

            // Copy remainder to beginning of buffer
            buf_offset = buf_size - str.len;
            memcpy(buf, str.ptr + str.len, buf_offset);
        }

        if (num_atoms_read != data->num_atoms) {
            md_printf(MD_LOG_TYPE_ERROR, "There was a descrepancy between the number of atoms read (%lli) compared to the number of atoms specified in the header (%lli)", num_atoms_read, data->num_atoms);
            goto done;
        }

        if (!parse_unitcell(str, data)) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to parse unitcell in gro file");
            goto done;
        }

        success = true;
done:
        md_file_close(file);
        return success;
    }
    md_printf(MD_LOG_TYPE_ERROR, "Could not open file '%.*s'", filename.len, filename.ptr);
    return false;
}

void md_gro_data_free(md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    if (data->atom_data) md_array_free(data->atom_data, alloc);
    memset(data, 0, sizeof(md_gro_data_t));
}

/*
static inline str_t add_label(str_t str, gro_molecule_t* inst) {
    ASSERT(inst);
    for (int64_t i = 0; i < md_array_size(inst->labels); ++i) {
        str_t label = { .ptr = inst->labels[i].str, .len = strnlen(inst->labels[i].str, ARRAY_SIZE(inst->labels[i].str)) };
        if (compare_str(str, label)) {
            return label;
        }
    }
    label_t label = {0};
    memcpy(label.str, str.ptr, MIN(str.len, (int64_t)ARRAY_SIZE(label.str) - 1));
    return (str_t){md_array_push(inst->labels, label, inst->allocator)->str, str.len};
}
*/

bool md_gro_molecule_init(struct md_molecule_t* mol, const md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(data);
    ASSERT(alloc);

    if (mol->inst) {
        md_print(MD_LOG_TYPE_DEBUG, "molecule inst object is not zero, potentially leaking memory when clearing");
    }

    memset(mol, 0, sizeof(md_molecule_t));

    mol->inst = md_alloc(alloc, sizeof(gro_molecule_t));
    gro_molecule_t* inst = (gro_molecule_t*)mol->inst;
    memset(inst, 0, sizeof(gro_molecule_t));

    inst->magic = MD_GRO_MOL_MAGIC;
    inst->allocator = md_arena_allocator_create(alloc, KILOBYTES(16));
    inst->labels = 0;

    const int64_t num_atoms = data->num_atoms;

    alloc = inst->allocator;

    md_array_ensure(mol->atom.x, num_atoms, alloc);
    md_array_ensure(mol->atom.y, num_atoms, alloc);
    md_array_ensure(mol->atom.z, num_atoms, alloc);
    md_array_ensure(mol->atom.name, num_atoms, alloc);
    md_array_ensure(mol->atom.residue_idx, num_atoms, alloc);

    int32_t cur_res_id = -1;
    for (int64_t i = 0; i < num_atoms; ++i) {
        const float x = data->atom_data[i].x;
        const float y = data->atom_data[i].y;
        const float z = data->atom_data[i].z;
        str_t atom_name = (str_t){data->atom_data[i].atom_name, strnlen(data->atom_data[i].atom_name, ARRAY_SIZE(data->atom_data[i].atom_name))};
        str_t res_name = (str_t){data->atom_data[i].res_name, strnlen(data->atom_data[i].res_name, ARRAY_SIZE(data->atom_data[i].res_name))};

        int32_t res_id = data->atom_data[i].res_id;
        if (res_id != cur_res_id) {
            cur_res_id = res_id;
            md_residue_id_t id = res_id;
            md_range_t atom_range = {(uint32_t)mol->atom.count, (uint32_t)mol->atom.count};

            mol->residue.count += 1;
            md_array_push(mol->residue.name, make_label(res_name), alloc);
            md_array_push(mol->residue.id, id, alloc);
            md_array_push(mol->residue.atom_range, atom_range, alloc);
        }

        if (mol->residue.atom_range) md_array_last(mol->residue.atom_range)->end += 1;

        mol->atom.count += 1;
        md_array_push(mol->atom.x, x, alloc);
        md_array_push(mol->atom.y, y, alloc);
        md_array_push(mol->atom.z, z, alloc);
        md_array_push(mol->atom.name, make_label(atom_name), alloc);

        if (mol->residue.count) md_array_push(mol->atom.residue_idx, (md_residue_idx_t)(mol->residue.count - 1), alloc);
    }

    mol->coord_frame = (mat3_t) {
        data->cell_ext[0] * 10.0f, 0, 0,
        0, data->cell_ext[1] * 10.0f, 0,
        0, 0, data->cell_ext[2] * 10.0f,
    };

    md_util_postprocess_molecule(mol, alloc);

    return true;
}

bool md_gro_molecule_free(struct md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(mol->inst);
    ASSERT(alloc);

    gro_molecule_t* inst = (gro_molecule_t*)mol->inst;
    if (inst->magic != MD_GRO_MOL_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "GRO magic did not match!");
        return false;
    }

    md_arena_allocator_destroy(inst->allocator);
    md_free(alloc, inst, sizeof(gro_molecule_t));
    memset(mol, 0, sizeof(md_molecule_t));

    return true;
}

static bool gro_init_from_str(md_molecule_t* mol, str_t str, md_allocator_i* alloc) {
    md_gro_data_t data = {0};
    bool success = false;
    if (md_gro_data_parse_str(&data, str, default_allocator)) {
        success = md_gro_molecule_init(mol, &data, alloc);
    }
    md_gro_data_free(&data, default_allocator);

    return success;
}

static bool gro_init_from_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc) {
    md_gro_data_t data = {0};
    bool success = false;
    if (md_gro_data_parse_file(&data, filename, default_allocator)) {
        success = md_gro_molecule_init(mol, &data, alloc);
    }
    md_gro_data_free(&data, default_allocator);

    return success;
}

static md_molecule_api gro_api = {
    gro_init_from_str,
    gro_init_from_file,
    md_gro_molecule_free,
};

md_molecule_api* md_gro_molecule_api() {
    return &gro_api;
}