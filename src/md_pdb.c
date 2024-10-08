#include <md_pdb.h>

#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_PDB_TRAJ_MAGIC 0x2312ad7b78a9bc78
#define MD_PDB_CACHE_MAGIC 0x89172bab124545
#define MD_PDB_CACHE_VERSION 15
#define MD_PDB_PARSE_BIOMT 1

// The opaque blob
typedef struct pdb_trajectory_t {
    uint64_t magic;
    md_file_o* file;
    int64_t* frame_offsets;
    md_unit_cell_t unit_cell;                // For pdb trajectories we have a static cell
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;
} pdb_trajectory_t;

static inline void copy_str_field(char* dst, size_t dst_size, str_t line, size_t beg, size_t end) {
    str_copy_to_char_buf(dst, dst_size, str_trim(str_substr(line, beg - 1, end-beg + 1)));
}

// We massage the beg and end indices here to correspond the pdb specification
// This makes our life easier when specifying all the different ranges
static inline int32_t extract_int(str_t line, size_t beg, size_t end) {
    if (line.len < end) return 0;
    return (int32_t)parse_int(str_substr(line, beg - 1, end-beg + 1));
}

static inline float extract_float(str_t line, size_t beg, size_t end) {
    if (line.len < end) return 0.0f;
    str_t str = str_trim(str_substr(line, beg - 1, end-beg + 1));
    return (float)parse_float_wide(str.ptr, str.len);
}

static inline char extract_char(str_t line, size_t idx) {
    if (line.len < idx) return 0;
    return line.ptr[idx - 1];
}

static inline md_pdb_coordinate_t extract_coord(str_t line) {
    char chain_id = extract_char(line, 22);
    if (!is_alpha(chain_id)) {
		chain_id = ' ';
	}
    md_pdb_coordinate_t coord = {
        .atom_serial = extract_int(line, 7, 11),
        .atom_name = {0},
        .alt_loc = extract_char(line, 17),
        .res_name = {0},
        .chain_id = chain_id,
        .res_seq = extract_int(line, 23, 26),
        .icode = extract_char(line, 27),
        .x = extract_float(line, 31, 38),
        .y = extract_float(line, 39, 46),
        .z = extract_float(line, 47, 54),
        .occupancy = extract_float(line, 55, 60),
        .temp_factor = extract_float(line, 61, 66),
        .element = {0},
        .charge = {0},
    };
    copy_str_field(coord.atom_name, sizeof(coord.atom_name),    line, 13, 16);
    copy_str_field(coord.res_name,  sizeof(coord.res_name),     line, 18, 21);
    copy_str_field(coord.element,   sizeof(coord.element),      line, 77, 78);
    copy_str_field(coord.charge,    sizeof(coord.charge),       line, 79, 80);
    return coord;
}

static inline md_pdb_helix_t extract_helix(str_t line) {
    md_pdb_helix_t helix = {
        .serial_number = extract_int(line, 8, 10),
        .id = {0},
        .init_res_name = {0},
        .init_chain_id = extract_char(line, 20),
        .init_res_seq = extract_int(line, 22, 25),
        .init_res_i_code = extract_char(line, 26),
        .end_res_name = {0},
        .end_chain_id = extract_char(line, 32),
        .end_res_seq = extract_int(line, 34, 37),
        .end_res_i_code = extract_char(line, 38),
        .helix_class = extract_int(line, 39, 40),
        .comment = {0},
        .length = extract_int(line, 72, 76),
    };
    copy_str_field(helix.id,            sizeof(helix.id),               line, 12, 14);
    copy_str_field(helix.init_res_name, sizeof(helix.init_res_name),    line, 16, 18);
    copy_str_field(helix.end_res_name,  sizeof(helix.end_res_name),     line, 28, 30);
    copy_str_field(helix.comment,       sizeof(helix.comment),          line, 41, 70);
    return helix;
}

static inline md_pdb_sheet_t extract_sheet(str_t line) {
    md_pdb_sheet_t sheet = {
        .strand = extract_int(line, 8, 10),
        .id = {0},
        .num_strands = extract_int(line, 15, 16),
        .init_res_name = {0},
        .init_chain_id = extract_char(line, 20),
        .init_res_seq = extract_int(line, 22, 25),
        .init_res_i_code = extract_char(line, 26),
        .end_res_name = {0},
        .end_chain_id = extract_char(line, 32),
        .end_res_seq = extract_int(line, 34, 37),
        .end_res_i_code = extract_char(line, 38),
        .sense = extract_int(line, 39, 40),
        .cur_atom_name = {0},
        .cur_res_name = {0},
        .cur_chain_id = extract_char(line, 50),
        .cur_res_seq = extract_int(line, 51, 54),
        .cur_res_i_code = extract_char(line, 55),
        .prev_atom_name = {0},
        .prev_res_name = {0},
        .prev_chain_id = extract_char(line, 65),
        .prev_res_seq = extract_int(line, 66, 69),
        .prev_res_i_code = extract_char(line, 70),
    };
    copy_str_field(sheet.id,            sizeof(sheet.id),               line, 12, 14);
    copy_str_field(sheet.init_res_name, sizeof(sheet.init_res_name),    line, 18, 20);
    copy_str_field(sheet.end_res_name,  sizeof(sheet.end_res_name),     line, 29, 31);
    copy_str_field(sheet.cur_atom_name, sizeof(sheet.cur_atom_name),    line, 42, 45);
    copy_str_field(sheet.cur_res_name,  sizeof(sheet.cur_res_name),     line, 46, 48);
    copy_str_field(sheet.prev_atom_name,sizeof(sheet.prev_atom_name),   line, 57, 60);
    copy_str_field(sheet.prev_res_name, sizeof(sheet.prev_res_name),    line, 61, 63);
    return sheet;
}

static inline md_pdb_connect_t extract_connect(str_t line) {
    line = str_trim(line);
    md_pdb_connect_t c = {
        .atom_serial = {
            [0] = extract_int(line, 7, 11),
            [1] = extract_int(line, 12, 16),
            [2] = (line.len > 20) ? extract_int(line, 17, 21) : 0,
            [3] = (line.len > 25) ? extract_int(line, 22, 26) : 0,
            [4] = (line.len > 30) ? extract_int(line, 27, 31) : 0,
        },
    };
    return c;
}

static inline struct md_pdb_cryst1_t extract_cryst1(str_t line) {
    line = str_trim(line);
    md_pdb_cryst1_t c = {
        .a = extract_float(line, 7,  15),
        .b = extract_float(line, 16, 24),
        .c = extract_float(line, 25, 33),
        .alpha = extract_float(line, 34, 40),
        .beta  = extract_float(line, 41, 47),
        .gamma = extract_float(line, 48, 54),
        .space_group = {0},
        .z = extract_int(line, 67, 70),
    };
    copy_str_field(c.space_group, ARRAY_SIZE(c.space_group), line, 56, 66);
    return c;
}

static inline void append_connect(md_pdb_connect_t* connect, const md_pdb_connect_t* other) {
    MEMCPY(&connect->atom_serial[5], &other->atom_serial[1], 3 * sizeof(int32_t));
}

bool pdb_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    pdb_trajectory_t* pdb = (pdb_trajectory_t*)inst;
    ASSERT(pdb);
    ASSERT(pdb->magic == MD_PDB_TRAJ_MAGIC);
    ASSERT(header);

    *header = pdb->header;
    return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
// Returns size in bytes of frame, frame_data_ptr is optional and is the destination to write the frame data to.
size_t pdb_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
    pdb_trajectory_t* pdb = (pdb_trajectory_t*)inst;
    ASSERT(pdb);
    ASSERT(pdb->magic == MD_PDB_TRAJ_MAGIC);

    if (!pdb->file) {
        MD_LOG_ERROR("File handle is NULL");
        return 0;
    }

    if (!pdb->frame_offsets) {
        MD_LOG_ERROR("Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || (int64_t)pdb->header.num_frames <= frame_idx) {
        MD_LOG_ERROR("Frame index is out of range");
        return 0;
    }

    const int64_t beg = pdb->frame_offsets[frame_idx + 0];
    const int64_t end = pdb->frame_offsets[frame_idx + 1];
    const size_t frame_size = (size_t)MAX(0, end - beg);

    if (frame_data_ptr) {
        ASSERT(pdb->file);
        md_mutex_lock(&pdb->mutex);
        md_file_seek(pdb->file, beg, MD_FILE_BEG);
        const size_t bytes_read = md_file_read(pdb->file, frame_data_ptr, frame_size);
        (void)bytes_read;
        md_mutex_unlock(&pdb->mutex);
        ASSERT(frame_size == bytes_read);
    }

    return frame_size;
}

bool pdb_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, size_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);
    ASSERT(frame_data_ptr);
    ASSERT(frame_data_size);

    pdb_trajectory_t* pdb = (pdb_trajectory_t*)inst;
    if (pdb->magic != MD_PDB_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame header, pdb magic did not match");
        return false;
    }

    str_t str = { .ptr = (char*)frame_data_ptr, .len = frame_data_size };
    str_t line;

    int32_t step = 0;
    if (str_extract_line(&line, &str)) {
        if (str_eq_cstr_n(line, "MODEL", 5)) {
            step = (int32_t)parse_int(str_substr(line, 10, 4));
        }
    }

    size_t i = 0;
    while (str_extract_line(&line, &str) && i < pdb->header.num_atoms) {
        if (line.len < 6) continue;
        if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
            if (x) { x[i] = extract_float(line, 31, 38); }
            if (y) { y[i] = extract_float(line, 39, 46); }
            if (z) { z[i] = extract_float(line, 47, 54); }

            i += 1;
        }
    }

    if (header) {
        header->num_atoms = i;
        header->index = step;
        header->timestamp = (double)(step-1); // This information is missing from PDB trajectories
        header->unit_cell = pdb->unit_cell;
    }

    return true;
}

// PUBLIC PROCEDURES

bool pdb_parse(md_pdb_data_t* data, md_buffered_reader_t* reader, struct md_allocator_i* alloc, bool stop_after_first_model) {
    ASSERT(data);
    ASSERT(reader);
    ASSERT(alloc);

    str_t line;
    while (md_buffered_reader_extract_line(&line, reader)) {
        if (line.len < 6) continue;
        if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
            md_pdb_coordinate_t coord = extract_coord(line);
            coord.flags |= line.ptr[0] == 'H' ? MD_PDB_COORD_FLAG_HETATM : 0;
            md_array_push(data->atom_coordinates, coord, alloc);
        }
        else if (str_eq_cstr_n(line, "TER", 3)) {
			md_pdb_coordinate_t* last = md_array_last(data->atom_coordinates);
            if (last) {
				last->flags |= MD_PDB_COORD_FLAG_TERMINATOR;
            }
        }
        else if (str_eq_cstr_n(line, "HELIX", 5)) {
            md_pdb_helix_t helix = extract_helix(line);
            md_array_push(data->helices, helix, alloc);
        }
        else if (str_eq_cstr_n(line, "SHEET", 5)) {
            md_pdb_sheet_t sheet = extract_sheet(line);
            md_array_push(data->sheets, sheet, alloc);
        }
        else if (str_eq_cstr_n(line, "CONECT", 6)) {
            md_pdb_connect_t connect = extract_connect(line);
            md_pdb_connect_t* last = md_array_last(data->connections);
            if (last && last->atom_serial[0] == connect.atom_serial[0]) {
                append_connect(last, &connect);
            }
            else {
                md_array_push(data->connections, connect, alloc);
            }
        }
        else if (str_eq_cstr_n(line, "MODEL", 5)) {
            if (stop_after_first_model && md_array_size(data->models) > 0) {
                break;
            }
            const int32_t num_coords = (int32_t)md_array_size(data->atom_coordinates);
            // This is a bit nasty, we want to get the correct offset to the beginning of the current line.
            // Therefore we need to do some pointer arithmetic because just using the length of the line may not get us
            // all the way back in case there were skipped \r characters.
            const int64_t offset = (size_t)MAX(0, md_buffered_reader_tellg(reader) - (reader->str.ptr - line.ptr));
            md_pdb_model_t model = {
                .serial = (int32_t)parse_int(str_substr(line, 6, SIZE_MAX)),
                .beg_atom_serial = num_coords > 0 ? data->atom_coordinates[num_coords-1].atom_serial : 1,
                .beg_atom_index = num_coords,
                .byte_offset = offset,
            };
            md_array_push(data->models, model, alloc);
        }
        else if (str_eq_cstr_n(line, "ENDMDL", 6)) {
            const int32_t num_coords = (int32_t)md_array_size(data->atom_coordinates);
            md_pdb_model_t* model = md_array_last(data->models);
            if (model) {
                model->end_atom_serial = num_coords > 0 ? data->atom_coordinates[num_coords-1].atom_serial : 1;
                model->end_atom_index = num_coords;
            }
        }
        else if (str_eq_cstr_n(line, "CRYST1", 6)) {
            md_pdb_cryst1_t cryst1 = extract_cryst1(line);
            md_array_push(data->cryst1, cryst1, alloc);
        }
        else if (str_eq_cstr_n(line, "REMARK 350", 10)) {
            if (str_eq_cstr(str_substr(line, 11, 12), "BIOMOLECULE:")) {
                md_pdb_assembly_t assembly = {0};
                assembly.id = (int32_t)parse_int(str_substr(line, 23, SIZE_MAX));
                assembly.transform_offset = (int32_t)md_array_size(data->transforms);
                md_array_push(data->assemblies, assembly, alloc);
            }
            else if (str_eq_cstr(str_substr(line, 11, 30), "APPLY THE FOLLOWING TO CHAINS:") ||
                     str_eq_cstr(str_substr(line, 11, 30), "                   AND CHAINS:")) {
                md_pdb_assembly_t* assembly = md_array_last(data->assemblies);
                if (assembly) {
                    str_t chains = str_trim(str_substr(line, 41, 30));
                    str_t chain;

                    // Find the next index of insertion for chain
                    size_t next_idx = 0;
                    while (next_idx < ARRAY_SIZE(assembly->apply_to_chains) && assembly->apply_to_chains[next_idx] != 0) {
                        next_idx += 1;
                    }
                    while (extract_token_delim(&chain, &chains, ',') && next_idx < ARRAY_SIZE(assembly->apply_to_chains)) {
                        assembly->apply_to_chains[next_idx] = str_trim(chain).ptr[0];
                        next_idx += 1;
                    }
                } else {
                    md_log(MD_LOG_TYPE_DEBUG, "Error while parsing PDB assembly, assembly is missing");
                }
            }
            else if (str_eq_cstr(str_substr(line, 13, 5), "BIOMT")) {
                md_pdb_assembly_t* assembly = md_array_last(data->assemblies);
                if (assembly) {
                    const int32_t row_idx = char_to_digit(line.ptr[18]) - 1;
                    ASSERT(-1 < row_idx && row_idx < 4);

                    if (row_idx == 0) {
                        mat4_t transform = mat4_ident();
                        md_array_push(data->transforms, transform, alloc);
                        assembly->transform_count += 1;
                    }
                    mat4_t* transform = md_array_last(data->transforms);
                    if (transform) {
                        transform->elem[0][row_idx] = (float)parse_float(str_trim(str_substr(line, 23, 10)));
                        transform->elem[1][row_idx] = (float)parse_float(str_trim(str_substr(line, 33, 10)));
                        transform->elem[2][row_idx] = (float)parse_float(str_trim(str_substr(line, 43, 10)));
                        transform->elem[3][row_idx] = (float)parse_float(str_trim(str_substr(line, 53, SIZE_MAX)));
                    } else {
                        md_log(MD_LOG_TYPE_DEBUG, "Error while parsing PDB assembly, transform is missing");
                    }
                } else {
                    md_log(MD_LOG_TYPE_DEBUG, "Error while parsing PDB assembly, assembly is missing");
                }
            }
        }
    }

    data->num_models            = md_array_size(data->models);
    data->num_atom_coordinates  = md_array_size(data->atom_coordinates);
    data->num_helices           = md_array_size(data->helices);
    data->num_sheets            = md_array_size(data->sheets);
    data->num_connections       = md_array_size(data->connections);
    data->num_cryst1            = md_array_size(data->cryst1);
    data->num_assemblies        = md_array_size(data->assemblies);
    data->num_transforms        = md_array_size(data->transforms);

    return true;
}

bool md_pdb_data_parse_str(md_pdb_data_t* data, str_t str, md_allocator_i* alloc) {
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);
    return pdb_parse(data, &reader, alloc, false);
}

bool md_pdb_data_parse_file(md_pdb_data_t* data, str_t filename, md_allocator_i* alloc) {
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        const size_t cap = MEGABYTES(1);
        char* buf = md_alloc(md_get_heap_allocator(), cap);
        
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        result = pdb_parse(data, &reader, alloc, false);
        
        md_free(md_get_heap_allocator(), buf, cap);
        md_file_close(file);
    } else {
        MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
    }
    return result;
}

void md_pdb_data_free(md_pdb_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    if (data->cryst1)               md_array_free(data->cryst1, alloc);
    if (data->models)               md_array_free(data->models, alloc);
    if (data->atom_coordinates)     md_array_free(data->atom_coordinates, alloc);
    if (data->helices)              md_array_free(data->helices, alloc);
    if (data->sheets)               md_array_free(data->sheets, alloc);
    if (data->connections)          md_array_free(data->connections, alloc);
    if (data->assemblies)           md_array_free(data->assemblies, alloc);
}

md_chain_idx_t find_chain_idx(char chain_id, md_molecule_t* mol) {
    for (md_chain_idx_t i = 0; i < (md_chain_idx_t)mol->chain.count; ++i) {
        if (mol->chain.id[i].buf[0] == chain_id) return i;
    }
    return -1;
}

bool md_pdb_molecule_init(md_molecule_t* mol, const md_pdb_data_t* data, md_pdb_options_t options, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(data);
    ASSERT(alloc);

    bool result = false;

    md_allocator_i* temp_alloc = md_vm_arena_create(GIGABYTES(4));

    size_t beg_atom_index = 0;
    size_t end_atom_index = data->num_atom_coordinates;

    // if we have more than one model, interperet it as a trajectory and only load the first model
    if (data->num_models > 0 && !(options & MD_PDB_OPTION_CONCAT_MODELS)) {
        // Limit the scope of atom coordinate entries if we have a trajectory (only consider first model)
        beg_atom_index = data->models[0].beg_atom_index;
        end_atom_index = data->models[0].end_atom_index;
    }

    MEMSET(mol, 0, sizeof(md_molecule_t));

    size_t capacity = ROUND_UP(end_atom_index - beg_atom_index, 16);

    md_array(md_label_t) chain_ids = 0;
    md_array(md_range_t) chain_ranges = 0;

    md_array_ensure(mol->atom.x,       capacity, alloc);
    md_array_ensure(mol->atom.y,       capacity, alloc);
    md_array_ensure(mol->atom.z,       capacity, alloc);
    md_array_ensure(mol->atom.element, capacity, alloc);
    md_array_ensure(mol->atom.type,    capacity, alloc);
    md_array_ensure(mol->atom.resid,   capacity, alloc);
    md_array_ensure(mol->atom.resname, capacity, alloc);
    md_array_ensure(mol->atom.flags,   capacity, alloc);

    int32_t prev_res_id = -1;
    char prev_chain_id = -1;

    for (size_t i = beg_atom_index; i < end_atom_index; ++i) {
        float x = data->atom_coordinates[i].x;
        float y = data->atom_coordinates[i].y;
        float z = data->atom_coordinates[i].z;
        
        str_t atom_type = str_from_cstrn(data->atom_coordinates[i].atom_name, sizeof(data->atom_coordinates[i].atom_name));
        str_t res_name  = str_from_cstrn(data->atom_coordinates[i].res_name,  sizeof(data->atom_coordinates[i].res_name));
        md_residue_id_t res_id = data->atom_coordinates[i].res_seq;
        char chain_id = data->atom_coordinates[i].chain_id;
        md_flags_t flags = 0;
        md_element_t element = 0;

        if (data->atom_coordinates[i].element[0]) {
            // If the element is available, we directly try to use that to lookup the element
            str_t element_str = str_from_cstrn(data->atom_coordinates[i].element, sizeof(data->atom_coordinates[i].element));
            element = md_util_element_lookup_ignore_case(element_str);
            if (element == 0) {
                MD_LOG_INFO("Failed to lookup element with string '" STR_FMT "'", STR_ARG(element_str));
            }
        }

        if (data->atom_coordinates[i].flags & MD_PDB_COORD_FLAG_TERMINATOR) {
            flags |= MD_FLAG_CHAIN_END;
        }

        if (chain_id != ' ') {
            // Populate chain_ids with empty ' '
            while (md_array_size(chain_ids) < mol->atom.count) {
                md_array_push(chain_ids, (md_label_t){0}, temp_alloc);
            }
            str_t chain_str = { &data->atom_coordinates[i].chain_id, 1 };
            md_array_push(chain_ids, make_label(chain_str), temp_alloc);

            if (chain_id != prev_chain_id) {
                if (prev_chain_id != -1) {
                    *md_array_last(mol->atom.flags) |= MD_FLAG_CHAIN_END;
                }
                flags |= MD_FLAG_CHAIN_BEG;
                md_array_push(chain_ranges, ((md_range_t){(int)i, (int)i}), temp_alloc);
                prev_chain_id = chain_id;
            }

            ASSERT(chain_ranges);
            md_array_last(chain_ranges)->end += 1;
        }

        if (res_id != prev_res_id) {
            if (prev_res_id != -1) {
				*md_array_last(mol->atom.flags) |= MD_FLAG_RES_END;
			}

            flags |= MD_FLAG_RES_BEG;
            prev_res_id = res_id;
        }

        mol->atom.count += 1;
        md_array_push_no_grow(mol->atom.x, x);
        md_array_push_no_grow(mol->atom.y, y);
        md_array_push_no_grow(mol->atom.z, z);
        md_array_push_no_grow(mol->atom.element, element);
        md_array_push_no_grow(mol->atom.type, make_label(atom_type));
        md_array_push_no_grow(mol->atom.flags, flags);
        md_array_push_no_grow(mol->atom.resname, make_label(res_name));
        md_array_push_no_grow(mol->atom.resid, res_id);
    }

    if (chain_ids) {
        while (md_array_size(chain_ids) < mol->atom.count) {
            md_array_push(chain_ids, (md_label_t){0}, temp_alloc);
        }
        // Copy chain ids if they are filled (meaning they had some form of entry)
        md_array_push_array(mol->chain.id, chain_ids, md_array_size(chain_ids), alloc);
    }

    if (mol->atom.count > 0) {
        md_flags_t flags = MD_FLAG_RES_END;
        if (mol->chain.id && !str_eq(LBL_TO_STR(mol->chain.id[mol->atom.count - 1]), STR_LIT(" "))) {
            flags |= MD_FLAG_CHAIN_END;
        }
        mol->atom.flags[mol->atom.count - 1] |= flags;
    }

    if (data->num_cryst1 > 0) {
        // Use first crystal
        const md_pdb_cryst1_t* cryst = &data->cryst1[0];
        if (cryst->a == 1 && cryst->b == 1 && cryst->c == 1 && cryst->alpha == 90 && cryst->beta == 90 && cryst->gamma == 90) {
            // This is a special case:
            // This is the identity matrix, and in such case, we assume there is no unit cell (no periodic boundary conditions)
        } else {
            mol->unit_cell = md_util_unit_cell_from_extent_and_angles(cryst->a, cryst->b, cryst->c, cryst->alpha, cryst->beta, cryst->gamma);
        }
    };

    // Create instances from assemblies
    for (size_t aidx = 0; aidx < data->num_assemblies; ++aidx) {
        const md_pdb_assembly_t* assembly = &data->assemblies[aidx];
        md_label_t instance_label = {0};
        instance_label.len = (uint8_t)snprintf(instance_label.buf, sizeof(instance_label.buf), "ASM_%i", assembly->id);

        for (uint32_t tidx = assembly->transform_offset; tidx < assembly->transform_offset + assembly->transform_count; ++tidx) {
            if (mat4_equal(data->transforms[tidx], mat4_ident())) {
                // Do not add the identity transform as an instance, that is implicit in the structure
                // Only add the additional instances
                continue;
            }
            md_range_t instance_range = {0,0};
            for (size_t cidx = 0; cidx < ARRAY_SIZE(data->assemblies[aidx].apply_to_chains); ++cidx) {
                // A transform can be applied to multiple chains,
                // We extract the consecutive ranges of the chains and possibly split this into multiple instances with the same label and transform

                char chain_id = data->assemblies[aidx].apply_to_chains[cidx];
                if (chain_id == 0) {
                    break;
                }

                int chain_idx = -1;
                for (size_t i = 0; i < md_array_size(chain_ranges); ++i) {
                    if (mol->chain.id[chain_ranges[i].beg].buf[0] == chain_id) {
						chain_idx = (int)i;
						break;
					}
                }
                if (chain_idx != -1) {
                    if (instance_range.beg == 0 && instance_range.end == 0) {
                        instance_range = chain_ranges[chain_idx];
                    } else {
                        // Append if possible
                        if (instance_range.end == chain_ranges[chain_idx].beg) {
                            instance_range.end = chain_ranges[chain_idx].end;
                        } else {
                            // Discontinous range, we need to commit and reset the range
                            md_array_push(mol->instance.atom_range, instance_range, alloc);
                            md_array_push(mol->instance.label, instance_label, alloc);
                            md_array_push(mol->instance.transform, data->transforms[tidx], alloc);
                            instance_range = chain_ranges[chain_idx];
                        }
                    }
                } else {
                    // Fatal error
                    MD_LOG_ERROR("Malformed PDB dataset, an assembly refers to a non existing chain id");
                    goto done;
                }
            }

            md_array_push(mol->instance.atom_range, instance_range, alloc);
            md_array_push(mol->instance.label, instance_label, alloc);
            md_array_push(mol->instance.transform, data->transforms[tidx], alloc);
        }
    }

    // Do not use the connectivity information from the PDB, it is often not complete and only provides the connectivity for HETATM
    // Instead, we leave the connectivity information empty and the postprocessing figure out connectivity.
#if 0
    for (int64_t j = 0; j < data->num_connections; ++j) {
        const md_pdb_connect_t* c = &data->connections[j];
        const int32_t a = c->atom_serial[0] - 1;
        if (a < 0 || (int32_t)mol->atom.count <= a) {
            goto error;
        }

		for (int32_t i = 1; i < (int32_t)ARRAY_SIZE(c->atom_serial); ++i) {
            const int32_t b = c->atom_serial[i] - 1;
            if (b == -1) break;
            if ((int32_t)mol->atom.count <= b) {
                goto error;
            }
            md_bond_t bond = {a, b};
            md_array_push(mol->bonds, bond, alloc);
        }
        continue;
    error:
        MD_LOG_ERROR("Malformed PDB dataset, connection refers to a non existing atom");
        return false;
    }
#endif

    mol->instance.count = md_array_size(mol->instance.transform);
    ASSERT(md_array_size(mol->instance.label) == mol->instance.count);
    ASSERT(md_array_size(mol->instance.atom_range) == mol->instance.count);

    result = true;
done:
    md_vm_arena_destroy(temp_alloc);
    return result;
}

static bool pdb_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_pdb_data_t data = {0};

    md_pdb_data_parse_str(&data, str, md_get_heap_allocator());
    bool success = md_pdb_molecule_init(mol, &data, MD_PDB_OPTION_NONE, alloc);
    md_pdb_data_free(&data, md_get_heap_allocator());
    
    return success;
}

static bool pdb_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        const size_t cap = MEGABYTES(1);
        char *buf = md_alloc(md_get_heap_allocator(), cap);
        ASSERT(buf);
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        
        md_pdb_data_t data = {0};
        bool parse   = pdb_parse(&data, &reader, md_get_heap_allocator(), true);
        bool success = parse && md_pdb_molecule_init(mol, &data, MD_PDB_OPTION_NONE, alloc);
        
        md_pdb_data_free(&data, md_get_heap_allocator());
        md_free(md_get_heap_allocator(), buf, cap);
        
        return success;
    }
    MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
    return false;
}

static md_molecule_loader_i pdb_molecule_api = {
    pdb_init_from_str,
    pdb_init_from_file,
};

md_molecule_loader_i* md_pdb_molecule_api(void) {
    return &pdb_molecule_api;
}

bool pdb_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    pdb_trajectory_t* pdb = (pdb_trajectory_t*)inst;
    if (pdb->magic != MD_PDB_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    bool result = false;
    const size_t frame_size = pdb_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        md_allocator_i* alloc = md_get_heap_allocator();
        void* frame_data = md_alloc(alloc, frame_size);
        ASSERT(frame_data);

        const size_t read_size = pdb_fetch_frame_data(inst, frame_idx, frame_data);
        if (read_size != frame_size) {
            MD_LOG_ERROR("Failed to read the expected size");
            goto done;
        }

        result = pdb_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
    done:
        md_free(alloc, frame_data, frame_size);
    }

    return result;
}

typedef struct pdb_cache_t {
    md_trajectory_cache_header_t header;
    md_unit_cell_t cell;
    int64_t* frame_offsets;
} pdb_cache_t;

static bool try_read_cache(pdb_cache_t* cache, str_t cache_file, size_t traj_num_bytes, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);

    bool result = false;
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        if (md_file_read(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
            MD_LOG_ERROR("PDB trajectory cache: failed to read header");
            goto done;
        }

        if (cache->header.magic != MD_PDB_CACHE_MAGIC) {
            MD_LOG_ERROR("PDB trajectory cache: magic was incorrect or corrupt");
            goto done;
        }
        if (cache->header.version != MD_PDB_CACHE_VERSION) {
            MD_LOG_INFO("PDB trajectory cache: version mismatch, expected %i, got %i", MD_PDB_CACHE_VERSION, (int)cache->header.version);
        }
        if (cache->header.num_bytes != traj_num_bytes) {
            MD_LOG_INFO("PDB trajectory cache: trajectory size mismatch, expected %zu, got %zu", traj_num_bytes, cache->header.num_bytes);
        }
        if (cache->header.num_atoms == 0) {
            MD_LOG_ERROR("PDB trajectory cache: num atoms was zero");
            goto done;
        }
        if (cache->header.num_frames == 0) {
            MD_LOG_ERROR("PDB trajectory cache: num frames was zero");
            goto done;
        }

        if (md_file_read(file, &cache->cell, sizeof(cache->cell)) != sizeof(cache->cell)) {
			MD_LOG_ERROR("PDB trajectory cache: failed to read unit cell");
			goto done;
		}

        const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
        cache->frame_offsets = md_alloc(alloc, offset_bytes);
        if (md_file_read(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
            MD_LOG_ERROR("PDB trajectory cache: Failed to read offset data");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            goto done;
        }

        // Test position in file, we expect to be at the end of the file
        if (md_file_tell(file) != (int64_t)md_file_size(file)) {
            MD_LOG_ERROR("PDB trajectory cache: file position was not at the end of the file");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            goto done;
        }

        result = true;
    done:
        md_file_close(file);
    }
    return result;
}

static bool write_cache(const pdb_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_INFO("PDB trajectory cache: could not open file '"STR_FMT"'", STR_ARG(cache_file));
        return false;
    }

    if (md_file_write(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
        MD_LOG_ERROR("PDB trajectory cache: failed to write header");
        goto done;
    }

    if (md_file_write(file, &cache->cell, sizeof(cache->cell)) != sizeof(cache->cell)) {
    	MD_LOG_ERROR("PDB trajectory cache: failed to write unit cell");
    	goto done;
    }

    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    if (md_file_write(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
        MD_LOG_ERROR("PDB trajectory cache: failed to write frame offsets");
        goto done;
    }

    result = true;

done:
    md_file_close(file);
    return result;
}

md_trajectory_i* md_pdb_trajectory_create(str_t filename, struct md_allocator_i* ext_alloc, uint32_t flags) {
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_ERROR("Failed to open file for PDB trajectory");
        return false;
    }

    const size_t filesize = md_file_size(file);
    md_file_close(file);

    md_strb_t sb = md_strb_create(md_get_temp_allocator());
    md_strb_fmt(&sb, STR_FMT ".cache", STR_ARG(filename));
    str_t cache_file = md_strb_to_str(sb);

    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    pdb_cache_t cache = {0};
    if (!try_read_cache(&cache, cache_file, filesize, alloc)) {
        md_allocator_i* temp_alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

        bool result = false;
        md_pdb_data_t data = {0};
        if (!md_pdb_data_parse_file(&data, filename, temp_alloc)) {
            goto cleanup;
        }

        if (data.num_models <= 1) {
            MD_LOG_INFO("The PDB file does not contain multiple model entries and cannot be read as a trajectory");
            goto cleanup;
        }

        // Validate the models
        const size_t num_atoms = (size_t)MAX(0, data.models[0].end_atom_index - data.models[0].beg_atom_index);
        if (!num_atoms) {
			MD_LOG_ERROR("The PDB file models are empty and cannot be read as a trajectory");
			goto cleanup;
		}
        for (size_t i = 1; i < data.num_models; ++i) {
            const size_t length = (size_t)MAX(0, data.models[i].end_atom_index - data.models[i].beg_atom_index);
            if (length && length != num_atoms) {
                MD_LOG_ERROR("The PDB file models are empty or not of equal length and cannot be read as a trajectory");
                goto cleanup;
            }
        }
        
        cache.header.magic = MD_PDB_CACHE_MAGIC;
        cache.header.version = MD_PDB_CACHE_VERSION;
        cache.header.num_bytes = filesize;
        cache.header.num_atoms = num_atoms;
        cache.header.num_frames = data.num_models;

        cache.frame_offsets = md_alloc(alloc, (cache.header.num_frames + 1) * sizeof(int64_t));
        for (size_t i = 0; i < cache.header.num_frames; ++i) {
            cache.frame_offsets[i] = data.models[i].byte_offset;
        }
        cache.frame_offsets[cache.header.num_frames] = filesize;

        if (data.num_cryst1 > 0) {
            if (data.num_cryst1 > 1) {
                md_log(MD_LOG_TYPE_INFO, "The PDB file contains multiple CRYST1 entries, will pick the first one for determining the simulation box");
            }
            // If it is in fact a box, that will be handled as well
            cache.cell = md_util_unit_cell_from_extent_and_angles(data.cryst1[0].a, data.cryst1[0].b, data.cryst1[0].c, data.cryst1[0].alpha, data.cryst1[0].beta, data.cryst1[0].gamma);
        }

        if (!(flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
            // If we fail to write the cache, that's ok, we can inform about it, but do not halt
            if (write_cache(&cache, cache_file)) {
                MD_LOG_INFO("PDB: Successfully created cache file for '" STR_FMT "'", STR_ARG(cache_file));
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
        const size_t frame_size = (size_t)MAX(0, cache.frame_offsets[i + 1] - cache.frame_offsets[i]);
        max_frame_size = MAX(max_frame_size, frame_size);
    }

    md_array(double) frame_times = md_array_create(double, cache.header.num_frames, alloc);
    for (size_t i = 0; i < cache.header.num_frames; ++i) {
        frame_times[i] = (double)i;
    }

    void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(pdb_trajectory_t));
    ASSERT(mem);
    MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(pdb_trajectory_t));

    md_trajectory_i* traj = mem;
    pdb_trajectory_t* pdb = (pdb_trajectory_t*)(traj + 1);

    pdb->magic = MD_PDB_TRAJ_MAGIC;
    pdb->file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    pdb->frame_offsets = cache.frame_offsets;
    pdb->allocator = alloc;
    pdb->mutex = md_mutex_create();
    pdb->header = (md_trajectory_header_t) {
        .num_frames = cache.header.num_frames,
        .num_atoms = cache.header.num_atoms,
        .max_frame_data_size = max_frame_size,
        .time_unit = {0},
        .frame_times = frame_times,
    };
    pdb->unit_cell = cache.cell;

    traj->inst = (struct md_trajectory_o*)pdb;
    traj->get_header = pdb_get_header;
    traj->load_frame = pdb_load_frame;
    //traj->fetch_frame_data = pdb_fetch_frame_data;
    //traj->decode_frame_data = pdb_decode_frame;

    return traj;
}

void md_pdb_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    pdb_trajectory_t* pdb = (pdb_trajectory_t*)traj->inst;
    if (pdb->magic != MD_PDB_TRAJ_MAGIC) {
        MD_LOG_ERROR("Trajectory is not a valid PDB trajectory.");
        ASSERT(false);
        return;
    }
    
    if (pdb->file) md_file_close(pdb->file);
    md_mutex_destroy(&pdb->mutex);
    md_arena_allocator_destroy(pdb->allocator);
}

static md_trajectory_loader_i pdb_traj_loader = {
    md_pdb_trajectory_create,
    md_pdb_trajectory_free,
};

md_trajectory_loader_i* md_pdb_trajectory_loader(void) {
    return &pdb_traj_loader;
}

#ifdef __cplusplus
}
#endif
