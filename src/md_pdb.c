#include <md_pdb.h>

#include <md_system.h>
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
#include <core/md_hash.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

static md_trajectory_i* md_pdb_trajectory_create(str_t filename, struct md_allocator_i* ext_alloc, md_trajectory_flags_t flags);

#define MD_PDB_TRAJ_MAGIC 0x2312ad7b78a9bc78
#define MD_PDB_TRAJ_READER_MAGIC 0x2312ad7b78a9bc79
#define MD_PDB_CACHE_MAGIC 0x89172bab124545
#define MD_PDB_CACHE_VERSION 15
#define MD_PDB_PARSE_BIOMT 1

// The opaque blob
typedef struct pdb_trajectory_t {
    uint64_t magic;
	str_t filepath;
    int64_t* frame_offsets;
    md_unitcell_t unitcell;                // For pdb trajectories we have a static cell
    md_trajectory_header_t header;
    md_allocator_i* allocator;
} pdb_trajectory_t;

typedef struct pdb_reader_t {
    uint64_t magic;
    md_file_t file;
    const pdb_trajectory_t* traj;
    md_array(uint8_t) frame_data;
    md_allocator_i* arena;
} pdb_reader_t;

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

// Parse atom serial or residue sequence number, which can be "*****" for missing values
static inline int32_t parse_id(str_t line, size_t beg, size_t end) {
	size_t len = end - beg + 1;
    str_t str = str_trim(str_substr(line, beg - 1, len));
    if (str_eq_n(str, STR_LIT("*****"), len)) {
        return INT_MAX;
    }
    bool all_digits = true;
    for (size_t i = 0; i < str.len; ++i) {
        if (!is_digit(str.ptr[i])) {
            all_digits = false;
        }
    }
    if (all_digits)
	    return (int32_t)parse_int(str);
    else
		return (int32_t)parse_hex(str);
}

static inline md_pdb_coordinate_t extract_coord(str_t line) {
    char chain_id = extract_char(line, 22);
    if (!is_alpha(chain_id)) {
		chain_id = ' ';
	}
    md_pdb_coordinate_t coord = {
        .atom_serial = parse_id(line, 7, 11),
        .atom_name = {0},
        .alt_loc = extract_char(line, 17),
        .res_name = {0},
        .chain_id = chain_id,
        .res_seq = parse_id(line, 23, 26),
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
static size_t pdb_fetch_frame_data(const pdb_trajectory_t* pdb, md_file_t file, int64_t frame_idx, void* frame_data_ptr) {
    ASSERT(pdb);
    ASSERT(pdb->magic == MD_PDB_TRAJ_MAGIC);

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
    size_t frame_size = (size_t)MAX(0, end - beg);

    if (frame_data_ptr) {
        if (!md_file_valid(file)) {
            MD_LOG_ERROR("Failed to open file '" STR_FMT "'", STR_ARG(pdb->filepath));
            return 0;
	    }

        size_t bytes_read = md_file_read_at(file, beg, frame_data_ptr, frame_size);
		frame_size = (bytes_read == frame_size) ? frame_size : 0;
    }

    return frame_size;
}

static bool pdb_decode_frame_data(const pdb_trajectory_t* pdb, const void* frame_data, size_t frame_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(frame_data);
    ASSERT(frame_size);

    if (pdb->magic != MD_PDB_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame header, pdb magic did not match");
        return false;
    }

    str_t str = { .ptr = (char*)frame_data, .len = frame_size };
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
        header->unitcell = pdb->unitcell;
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
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Could not open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    const size_t cap = MEGABYTES(1);
    char* buf = md_alloc(md_get_heap_allocator(), cap);
        
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
    bool result = pdb_parse(data, &reader, alloc, false);
        
    md_free(md_get_heap_allocator(), buf, cap);
    md_file_close(&file);

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

bool md_pdb_system_init_from_data(md_system_t* sys, const md_pdb_data_t* data, md_pdb_options_t options) {
    ASSERT(sys);
    ASSERT(data);

    if (!sys->alloc) {
       MD_LOG_ERROR("System allocator not set");
       return false;
    }

    md_system_reset(sys);

    md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    size_t beg_atom_index = 0;
    size_t end_atom_index = data->num_atom_coordinates;

    // if we have more than one model, interperet it as a trajectory and only load the first model
    if (data->num_models > 0 && !(options & MD_PDB_OPTION_CONCAT_MODELS)) {
        // Limit the scope of atom coordinate entries if we have a trajectory (only consider first model)
        beg_atom_index = data->models[0].beg_atom_index;
        end_atom_index = data->models[0].end_atom_index;
    }

    bool result = false;

    size_t capacity = ROUND_UP(end_atom_index - beg_atom_index, 16);

    md_array(str_t) atom_name = 0;
    md_array_ensure(atom_name, capacity, temp_arena);
    // Keep track of the asymmetric unit id for each component, this serves as a basis for determining entities and instances
    md_array(str_t) comp_auth_asym_ids = 0;

    md_array_ensure(sys->atom.x,        capacity, sys->alloc);
    md_array_ensure(sys->atom.y,        capacity, sys->alloc);
    md_array_ensure(sys->atom.z,        capacity, sys->alloc);
    md_array_ensure(sys->atom.type_idx, capacity, sys->alloc);
    md_array_ensure(sys->atom.flags,    capacity, sys->alloc);

	// Add default unknown atom type at index 0
    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

	uint64_t prev_comp_key = 0;
    size_t num_unassigned_atom_types = 0;

    bool terminator = false;
    for (size_t i = beg_atom_index; i < end_atom_index; ++i) {
		char alt_loc = data->atom_coordinates[i].alt_loc;
		if (alt_loc != ' ' && alt_loc != 'A') continue; // Ignore alt loc entries unless 'A' == first alt loc

        float x = data->atom_coordinates[i].x;
        float y = data->atom_coordinates[i].y;
        float z = data->atom_coordinates[i].z;
        
        str_t atom_id  = str_from_cstrn(data->atom_coordinates[i].atom_name, sizeof(data->atom_coordinates[i].atom_name));
        str_t res_name = str_from_cstrn(data->atom_coordinates[i].res_name,  sizeof(data->atom_coordinates[i].res_name));
        md_sequence_id_t seq_id = data->atom_coordinates[i].res_seq;
        if (seq_id == INT_MAX) {
			// Missing residue sequence id, assign one based on previous residue
			seq_id = sys->component.seq_id ? sys->component.seq_id[sys->component.count - 1] + 1 : 1;
		}
        char chain_id = data->atom_coordinates[i].chain_id;
        md_flags_t flags = (data->atom_coordinates[i].flags & MD_PDB_COORD_FLAG_HETATM) ? MD_FLAG_HETERO : 0;
        md_atomic_number_t atomic_number = 0;
        md_atom_type_idx_t atom_type_idx = 0;

        if (data->atom_coordinates[i].element[0]) {
            // If the element is available, use that to lookup the element
            str_t sym = str_from_cstrn(data->atom_coordinates[i].element, sizeof(data->atom_coordinates[i].element));
            atomic_number = md_atomic_number_from_symbol(sym, true);
            float mass   = md_atomic_number_mass(atomic_number);
            float radius = md_atomic_number_vdw_radius(atomic_number);
            uint32_t color = md_atomic_number_cpk_color(atomic_number);
            atom_type_idx = md_atom_type_find_or_add(&sys->atom.type, atom_id, atomic_number, mass, radius, color, 0, sys->alloc);
        } else {
            num_unassigned_atom_types += 1;
        }

		uint64_t comp_key = md_hash64_str(res_name, md_hash64(&seq_id, sizeof(seq_id), chain_id));

        if (comp_key != prev_comp_key || terminator) {
            // New residue
            // Propagate HETERO flag to residue
            md_flags_t comp_flags = flags;

            sys->component.count += 1;
            md_array_push(sys->component.atom_offset, (uint32_t)sys->atom.count, sys->alloc);
            md_array_push(sys->component.name,   make_label(res_name), sys->alloc);
            md_array_push(sys->component.seq_id, seq_id, sys->alloc);
            md_array_push(sys->component.flags,  comp_flags, sys->alloc);

            str_t asym_id = str_trim(str_from_cstrn(&data->atom_coordinates[i].chain_id, 1));
            if (terminator) {
                // Trigger a new component for the next residue even if the chain id is the same, to correctly capture chain breaks
                asym_id.len = 2;
            }
            md_array_push(comp_auth_asym_ids, asym_id, temp_arena);
        }

        sys->atom.count += 1;
        md_array_push_no_grow(atom_name, atom_id);
        md_array_push_no_grow(sys->atom.x, x);
        md_array_push_no_grow(sys->atom.y, y);
        md_array_push_no_grow(sys->atom.z, z);
        md_array_push_no_grow(sys->atom.flags, flags);
        md_array_push_no_grow(sys->atom.type_idx, atom_type_idx);

		prev_comp_key = comp_key;

        terminator = (data->atom_coordinates[i].flags & MD_PDB_COORD_FLAG_TERMINATOR) != 0;
    }
    md_array_push(sys->component.atom_offset, (uint32_t)sys->atom.count, sys->alloc);  // Final sentinel

    if (data->num_cryst1 > 0) {
        // Use first crystal
        const md_pdb_cryst1_t* cryst = &data->cryst1[0];
        if (cryst->a == 1 && cryst->b == 1 && cryst->c == 1 && cryst->alpha == 90 && cryst->beta == 90 && cryst->gamma == 90) {
            // This is a special case:
            // This is the identity matrix, and in such case, we assume there is no unit cell (no periodic boundary conditions)
        } else {
            sys->unitcell = md_unitcell_from_extent_and_angles(cryst->a, cryst->b, cryst->c, cryst->alpha, cryst->beta, cryst->gamma);
        }
    };

    if (num_unassigned_atom_types > 0) {
        md_util_system_infer_atom_types(sys, atom_name, sys->alloc);
    }
	md_util_system_infer_covalent_bonds(sys, sys->alloc);
    md_util_system_infer_comp_flags(sys);
    md_util_system_infer_entity_and_instance(sys, comp_auth_asym_ids, sys->alloc);

    /*
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
            md_urange_t instance_range = {0,0};
            for (size_t cidx = 0; cidx < ARRAY_SIZE(data->assemblies[aidx].apply_to_chains); ++cidx) {
                // A transform can be applied to multiple chains,
                // We extract the consecutive ranges of the chains and possibly split this into multiple instances with the same label and transform

                char chain_id = data->assemblies[aidx].apply_to_chains[cidx];
                if (chain_id == 0) {
                    break;
                }

				// Find the chain index
                int chain_idx = -1;
                for (size_t i = 0; i < num_chains; ++i) {
                    int chain_id_i = chain_ids[i];
                    if (chain_id_i == chain_id) {
						chain_idx = (int)i;
						break;
					}
                }
                if (chain_idx != -1) {
                    if (instance_range.beg == 0 && instance_range.end == 0) {
                        instance_range = chain_atom_ranges[chain_idx];
                    } else {
                        // Append if possible
                        if (instance_range.end == chain_atom_ranges[chain_idx].beg) {
                            instance_range.end = chain_atom_ranges[chain_idx].end;
                        } else {
                            // Discontinous range, we need to commit and reset the range
                            md_array_push(mol->assembly.atom_range, instance_range, mol_alloc);
                            md_array_push(mol->assembly.label,      instance_label, mol_alloc);
                            md_array_push(mol->assembly.transform,  data->transforms[tidx], mol_alloc);
                            instance_range = chain_atom_ranges[chain_idx];
                        }
                    }
                } else {
                    // Fatal error
                    MD_LOG_ERROR("Malformed PDB dataset, an assembly refers to a non existing chain id");
                    goto done;
                }
            }

            md_array_push(mol->assembly.atom_range, instance_range, mol_alloc);
            md_array_push(mol->assembly.label, instance_label, mol_alloc);
            md_array_push(mol->assembly.transform, data->transforms[tidx], mol_alloc);
        }
    }

    mol->assembly.count = md_array_size(mol->assembly.transform);
    ASSERT(md_array_size(mol->assembly.label) == mol->assembly.count);
    ASSERT(md_array_size(mol->assembly.atom_range) == mol->assembly.count); 
    */

    result = true;
    md_arena_allocator_destroy(temp_arena);
    return result;
}

bool md_pdb_system_init_from_str(md_system_t* sys, str_t str, md_pdb_options_t options) {
    md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_pdb_data_t data = {0};
    md_pdb_data_parse_str(&data, str, temp_arena);
    bool success = md_pdb_system_init_from_data(sys, &data, options);
    
    md_arena_allocator_destroy(temp_arena);
    return success;
}

bool md_pdb_system_init_from_file(md_system_t* sys, str_t filename, md_pdb_options_t options) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Could not open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    md_allocator_i* temp_arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    const size_t cap = MEGABYTES(1);
    char *buf = md_arena_allocator_push(temp_arena, cap);
    ASSERT(buf);
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
    
    md_pdb_data_t data = {0};
    bool success = pdb_parse(&data, &reader, temp_arena, false) && md_pdb_system_init_from_data(sys, &data, options);

    // If the file contained multiple models, interpret as a trajectory and attach one.
    if (success && data.num_models > 1) {
        md_trajectory_flags_t traj_flags = MD_TRAJECTORY_FLAG_NONE;
        if (options & MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE) {
            traj_flags |= MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE;
        }
        md_trajectory_i* traj = md_pdb_trajectory_create(filename, sys->alloc, traj_flags);
        if (traj) {
            md_system_attach_trajectory(sys, traj);
        }
    }

    md_arena_allocator_destroy(temp_arena);

    return success;
}

bool pdb_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    (void)inst;
    (void)frame_idx;
    (void)header;
    (void)x;
    (void)y;
    (void)z;
    return false;
}

static bool pdb_reader_load_frame(struct md_trajectory_reader_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    pdb_reader_t* reader = (pdb_reader_t*)inst;
    ASSERT(reader->magic == MD_PDB_TRAJ_READER_MAGIC);

    const pdb_trajectory_t* pdb = reader->traj;
    if (pdb->magic != MD_PDB_TRAJ_MAGIC) {
        MD_LOG_ERROR("Error when decoding frame coord, pdb magic did not match");
        return false;
    }

    bool result = false;
    const size_t frame_size = pdb_fetch_frame_data(pdb, reader->file, frame_idx, NULL);
    if (frame_size > 0) {
        md_array_ensure(reader->frame_data, frame_size, reader->arena);
        const size_t read_size = pdb_fetch_frame_data(pdb, reader->file, frame_idx, reader->frame_data);
        if (read_size != frame_size) {
            MD_LOG_ERROR("Failed to read the expected size");
            return false;
        }

        result = pdb_decode_frame_data(pdb, reader->frame_data, frame_size, header, x, y, z);
    }

    return result;
}

static void pdb_trajectory_reader_free(struct md_trajectory_reader_i* reader) {
    if (!reader) {
        return;
    }

    pdb_reader_t* inst = (pdb_reader_t*)reader->inst;
    if (inst) {
        ASSERT(inst->magic == MD_PDB_TRAJ_READER_MAGIC);
        if (md_file_valid(inst->file)) {
            md_file_close(&inst->file);
        }
        md_arena_allocator_destroy(inst->arena);
    }

    MEMSET(reader, 0, sizeof(*reader));
}

static bool pdb_trajectory_reader_init(md_trajectory_reader_i* reader, struct md_trajectory_o* traj_inst) {
    ASSERT(reader);
    ASSERT(traj_inst);

    pdb_trajectory_t* pdb = (pdb_trajectory_t*)traj_inst;
    ASSERT(pdb->magic == MD_PDB_TRAJ_MAGIC);

    md_file_t file = {0};
    if (!md_file_open(&file, pdb->filepath, MD_FILE_READ)) {
        MD_LOG_ERROR("Failed to open file '" STR_FMT "'", STR_ARG(pdb->filepath));
        return false;
    }

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    pdb_reader_t* inst = md_alloc(arena, sizeof(pdb_reader_t));
    MEMSET(inst, 0, sizeof(pdb_reader_t));
    inst->magic = MD_PDB_TRAJ_READER_MAGIC;
    inst->file = file;
    inst->traj = pdb;
    inst->arena = arena;

    MEMSET(reader, 0, sizeof(*reader));
    reader->inst = (struct md_trajectory_reader_o*)inst;
    reader->free = pdb_trajectory_reader_free;
    reader->load_frame = pdb_reader_load_frame;
    return true;
}

typedef struct pdb_cache_t {
    md_trajectory_cache_header_t header;
    md_unitcell_t cell;
    int64_t* frame_offsets;
} pdb_cache_t;

static bool try_read_cache(pdb_cache_t* cache, str_t cache_file, size_t traj_num_bytes, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);

    bool result = false;
    md_file_t file = {0};
    if (md_file_open(&file, cache_file, MD_FILE_READ)) {
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
        md_file_close(&file);
    }
    return result;
}

static bool write_cache(const pdb_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_t file = {0};
    if (!md_file_open(&file, cache_file, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
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
    md_file_close(&file);
    return result;
}

static void md_pdb_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    pdb_trajectory_t* pdb = (pdb_trajectory_t*)traj->inst;
    if (pdb->magic != MD_PDB_TRAJ_MAGIC) {
        MD_LOG_ERROR("Trajectory is not a valid PDB trajectory.");
        ASSERT(false);
        return;
    }
    
    md_arena_allocator_destroy(pdb->allocator);
}

static md_trajectory_i* md_pdb_trajectory_create(str_t filename, struct md_allocator_i* ext_alloc, md_trajectory_flags_t flags) {
	md_file_info_t file_info = { 0 };
    if (!md_file_info_extract_from_path(filename, &file_info)) {
        MD_LOG_ERROR("Failed to extract file info from path '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    const size_t filesize = file_info.size;

	char cache_path_buf[4096];
	snprintf(cache_path_buf, sizeof(cache_path_buf), STR_FMT ".cache", STR_ARG(filename));
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
            cache.cell = md_unitcell_from_extent_and_angles(data.cryst1[0].a, data.cryst1[0].b, data.cryst1[0].c, data.cryst1[0].alpha, data.cryst1[0].beta, data.cryst1[0].gamma);
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
	pdb->filepath = str_copy(filename, alloc);
    pdb->frame_offsets = cache.frame_offsets;
    pdb->allocator = alloc;
    pdb->header = (md_trajectory_header_t) {
        .num_frames = cache.header.num_frames,
        .num_atoms = cache.header.num_atoms,
        .time_unit = {0},
        .frame_times = frame_times,
    };
    pdb->unitcell = cache.cell;

    traj->inst = (struct md_trajectory_o*)pdb;
    traj->free = md_pdb_trajectory_free;
    traj->get_header = pdb_get_header;
    traj->init_reader = pdb_trajectory_reader_init;

    return traj;
}

#ifdef __cplusplus
}
#endif
