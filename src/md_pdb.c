#include <md_pdb.h>

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
#include <core/md_vec_math.h>
#include <core/md_hash.h>

#include <hybrid_36_c.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_PDB_CACHE_MAGIC 0x89172bab124545
#define MD_PDB_CACHE_VERSION 17
#define MD_PDB_PARSE_BIOMT 1

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
    return (float)parse_float(str_substr(line, beg - 1, end-beg + 1));
}

// A column the format writes with a fixed number of decimals ("%8.3f" for coordinates, "%6.2f" for
// occupancy and temperature factor): read without a scan when it is laid out that way, and parsed like
// any other column when it is not
static inline float extract_fixed_float(str_t line, size_t beg, size_t end, size_t decimals) {
    if (line.len < end) return 0.0f;
    const str_t field = str_substr(line, beg - 1, end-beg + 1);
    float value;
    return md_parse_fixed_f32(&value, field, decimals) ? value : (float)parse_float(field);
}

static inline char extract_char(str_t line, size_t idx) {
    if (line.len < idx) return 0;
    return line.ptr[idx - 1];
}

// An atom serial or a residue sequence number: a decimal number, past the width of its field hybrid-36 (A0000, then
// a0000: a letter first). 0 where the field holds neither, counted in num_bad: "*****" past the width, which many
// writers put there, or the hexadecimal VMD writes past 99999 ("186a0" for 100000), which is not hybrid-36 (a number
// which begins with a digit is decimal) and cannot be told from decimal where it has no letter (18700, 100096). Such a
// file is read by the order of its atoms; only CONECT records refer to atoms by their serials.
static inline int32_t parse_id(str_t line, size_t beg, size_t end, uint32_t* num_bad) {
	uint32_t len = (uint32_t)(end - beg) + 1;
    int32_t  res = 0;
    if (hy36decode(len, line.ptr + beg - 1, len, &res)) {
        if (num_bad) *num_bad += 1;
        return 0;
    }
    return res;
}

static inline md_pdb_coordinate_t extract_coord(str_t line, uint32_t* num_bad_ids) {
    char chain_id = extract_char(line, 22);
    if (!is_alpha(chain_id)) {
		chain_id = ' ';
	}
    md_pdb_coordinate_t coord = {
        .atom_serial = parse_id(line, 7, 11, num_bad_ids),
        .atom_name = {0},
        .alt_loc = extract_char(line, 17),
        .res_name = {0},
        .chain_id = chain_id,
        .res_seq = parse_id(line, 23, 26, num_bad_ids),
        .icode = extract_char(line, 27),
        .x = extract_fixed_float(line, 31, 38, 3),
        .y = extract_fixed_float(line, 39, 46, 3),
        .z = extract_fixed_float(line, 47, 54, 3),
        .occupancy = extract_fixed_float(line, 55, 60, 2),
        .temp_factor = extract_fixed_float(line, 61, 66, 2),
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

// A serial of a CONECT record, 0 where there is none. Serials are hybrid-36 as those of the atoms are, past 99999 atoms.
static inline int32_t extract_connect_serial(str_t line, size_t beg, size_t end) {
    if (line.len < end) return 0;
    if (str_empty(str_trim(str_substr(line, beg - 1, end - beg + 1)))) return 0;
    return parse_id(line, beg, end, NULL);
}

static inline md_pdb_connect_t extract_connect(str_t line) {
    line = str_trim(line);
    md_pdb_connect_t c = {
        .atom_serial = {
            [0] = extract_connect_serial(line, 7, 11),
            [1] = extract_connect_serial(line, 12, 16),
            [2] = extract_connect_serial(line, 17, 21),
            [3] = extract_connect_serial(line, 22, 26),
            [4] = extract_connect_serial(line, 27, 31),
        },
    };
    return c;
}

// The record of the cell: CRYST1, and CRYST with a blank where the 1 is, as some tools write it (the Conformational
// Analysis Tools of CHARMM-GUI's trajectories, with all angles 0)
static inline bool is_cryst1_record(str_t line) {
    return str_eq_cstr_n(line, "CRYST", 5) && (line.len == 5 || line.ptr[5] == '1' || line.ptr[5] == ' ');
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
    // No angles at all is a rectangular box whose angles were left out, not a degenerate cell
    if (c.alpha == 0 && c.beta == 0 && c.gamma == 0) {
        c.alpha = c.beta = c.gamma = 90;
    }
    return c;
}

// A further CONECT record of the same atom: its bonded atoms go into the free slots, as many as there are
static inline void append_connect(md_pdb_connect_t* connect, const md_pdb_connect_t* other) {
    size_t n = 1;
    while (n < ARRAY_SIZE(connect->atom_serial) && connect->atom_serial[n]) ++n;
    for (size_t k = 1; k < ARRAY_SIZE(other->atom_serial) && other->atom_serial[k] && n < ARRAY_SIZE(connect->atom_serial); ++k) {
        connect->atom_serial[n++] = other->atom_serial[k];
    }
}

// PUBLIC PROCEDURES

bool pdb_parse(md_pdb_data_t* data, md_buffered_reader_t* reader, struct md_allocator_i* alloc, bool stop_after_first_model) {
    ASSERT(data);
    ASSERT(reader);
    ASSERT(alloc);

    str_t line;
    int32_t biomolecule = 0;    // REMARK 350: the one whose groups of chains follow
    uint32_t num_bad_ids = 0;   // Serials and residue numbers which are not numbers of the format (parse_id)
    while (md_buffered_reader_extract_line(&line, reader)) {
        if (line.len < 6) continue;
        if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
            md_pdb_coordinate_t coord = extract_coord(line, &num_bad_ids);
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
        else if (is_cryst1_record(line)) {
            md_pdb_cryst1_t cryst1 = extract_cryst1(line);
            md_array_push(data->cryst1, cryst1, alloc);
        }
        else if (str_eq_cstr_n(line, "REMARK 350", 10)) {
            // A biomolecule is one or more groups of chains, each followed by the BIOMT transforms which copy them
            const str_t rest = str_substr(line, 11, SIZE_MAX);
            if (str_eq_cstr_n(rest, "BIOMOLECULE:", 12)) {
                biomolecule = (int32_t)parse_int(str_trim(str_substr(rest, 12, SIZE_MAX)));
            }
            else if (str_eq_cstr_n(rest, "APPLY THE FOLLOWING TO CHAINS:", 30) || str_eq_cstr_n(str_trim(rest), "AND CHAINS:", 11)) {
                // A group begins with the chains to apply to, which AND CHAINS lines continue
                if (str_eq_cstr_n(rest, "APPLY", 5)) {
                    const md_pdb_assembly_t group = {
                        .id = biomolecule,
                        .chain_offset = (uint32_t)md_array_size(data->assembly_chains),
                        .transform_offset = (uint32_t)md_array_size(data->transforms),
                    };
                    md_array_push(data->assemblies, group, alloc);
                }
                md_pdb_assembly_t* group = md_array_last(data->assemblies);
                if (group) {
                    size_t colon = 0;
                    str_find_char(&colon, rest, ':');
                    str_t chains = str_substr(rest, colon + 1, SIZE_MAX);
                    str_t chain;
                    while (extract_token_delim(&chain, &chains, ',')) {
                        chain = str_trim(chain);
                        if (str_empty(chain)) continue;
                        md_array_push(data->assembly_chains, chain.ptr[0], alloc);
                        group->chain_count += 1;
                    }
                } else {
                    MD_LOG_DEBUG("REMARK 350: chains of no biomolecule");
                }
            }
            else if (str_eq_cstr_n(str_trim(rest), "BIOMT", 5) && line.len > 18) {
                // REMARK 350   BIOMTr   n  m[r][1] m[r][2] m[r][3]  t[r], the rows of transform n of the group
                md_pdb_assembly_t* group = md_array_last(data->assemblies);
                const int32_t row = char_to_digit(line.ptr[18]) - 1;
                str_t fields = str_substr(line, 19, SIZE_MAX);
                str_t tok[5];
                const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &fields);
                if (!group || row < 0 || row > 2 || num_tok < 5) {
                    MD_LOG_DEBUG("REMARK 350: a malformed BIOMT record");
                } else {
                    if (row == 0) {
                        md_array_push(data->transforms, mat4_ident(), alloc);
                        group->transform_count += 1;
                    }
                    mat4_t* transform = md_array_last(data->transforms);
                    for (int c = 0; c < 3; ++c) {
                        transform->elem[c][row] = (float)parse_float(tok[1 + c]);
                    }
                    transform->elem[3][row] = (float)parse_float(tok[4]);
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
    data->num_assembly_chains   = md_array_size(data->assembly_chains);
    data->num_transforms        = md_array_size(data->transforms);
    if (num_bad_ids > 0) {
        MD_LOG_DEBUG("%u atom serials or residue numbers are not numbers of the PDB format ('*****' past the width of the field, "
            "or hexadecimal as VMD writes past 99999): the atoms are taken in their order", num_bad_ids);
    }

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
    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    char* buf = md_temp_alloc(temp_scope, cap);
        
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
    bool result = pdb_parse(data, &reader, alloc, false);
        
    md_temp_end(temp_scope);
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
    if (data->assembly_chains)      md_array_free(data->assembly_chains, alloc);
    if (data->transforms)           md_array_free(data->transforms, alloc);
}

// The assemblies of REMARK 350: the chains of a group of a biomolecule are the instances of their ids (a chain, and
// the molecules it holds: its ligands, its waters). A transform in two groups (the identity of each) is one operator,
// so that the copy it makes is one, named by its number in the first group which has it.
static void pdb_assemblies(md_system_t* sys, const md_pdb_data_t* data, md_allocator_i* temp_alloc) {
    if (data->num_assemblies == 0 || data->num_transforms == 0) return;

    md_array(mat4_t)     transform = 0;
    md_array(md_label_t) oper_id = 0;
    md_array(md_assembly_gen_t) gen = 0;

    for (size_t a = 0; a < data->num_assemblies; ++a) {
        const md_pdb_assembly_t* group = data->assemblies + a;
        if (group->chain_offset + group->chain_count > data->num_assembly_chains || group->transform_offset + group->transform_count > data->num_transforms) {
            MD_LOG_ERROR("REMARK 350: a group of chains out of range");
            continue;
        }

        md_array(md_instance_idx_t) inst = 0;
        for (uint32_t c = 0; c < group->chain_count; ++c) {
            const char chain = data->assembly_chains[group->chain_offset + c];
            bool found = false;
            for (size_t i = 0; i < sys->instance.count; ++i) {
                const str_t auth = md_instance_auth_id(&sys->instance, i);
                if (auth.len == 1 && auth.ptr[0] == chain) {
                    md_array_push(inst, (md_instance_idx_t)i, temp_alloc);
                    found = true;
                }
            }
            if (!found) {
                MD_LOG_DEBUG("REMARK 350: biomolecule %d names chain '%c', which has no atoms", group->id, chain);
            }
        }

        md_array(uint32_t) oper = 0;
        for (uint32_t t = 0; t < group->transform_count; ++t) {
            const mat4_t M = data->transforms[group->transform_offset + t];
            size_t k = 0;
            for (; k < md_array_size(transform); ++k) {
                if (MEMCMP(&transform[k], &M, sizeof(mat4_t)) == 0) break;
            }
            if (k == md_array_size(transform)) {
                char buf[16];
                const int len = snprintf(buf, sizeof(buf), "%u", t + 1);
                md_array_push(transform, M, temp_alloc);
                md_array_push(oper_id, make_label((str_t){buf, (size_t)MAX(len, 0)}), temp_alloc);
            }
            md_array_push(oper, (uint32_t)k, temp_alloc);
        }

        char buf[16];
        const int len = snprintf(buf, sizeof(buf), "%d", group->id);
        const md_assembly_gen_t row = {
            .assembly_id = make_label((str_t){buf, (size_t)MAX(len, 0)}),
            .inst        = inst,
            .inst_count  = md_array_size(inst),
            .oper        = oper,
            .oper_count  = md_array_size(oper),
        };
        md_array_push(gen, row, temp_alloc);
    }

    md_system_assemblies_init(sys, gen, md_array_size(gen), transform, oper_id, md_array_size(transform));
}

bool md_pdb_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_pdb_data_t* data, md_pdb_options_t options) {
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

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

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

    // Accumulated over the atoms actually kept, not over the coordinate records: altLoc
    // entries are skipped below, so a column indexed by record would be misaligned from the
    // first alternate onwards.
    md_array(float) atom_occupancy = 0;
    md_array(float) atom_b_factor  = 0;
    md_array_ensure(atom_occupancy, capacity, temp_arena);
    md_array_ensure(atom_b_factor,  capacity, temp_arena);
    // Keep track of the asymmetric unit id for each residue, this serves as a basis for determining entities and instances
    md_array(str_t) res_auth_asym_ids = 0;

    md_array_ensure(state->xyz,           capacity, state->alloc);
    md_array_ensure(sys->atom.type_idx, capacity, sys->alloc);
    md_array_ensure(sys->atom.flags,    capacity, sys->alloc);

	// Add default unknown particle type at index 0
    md_particle_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
    md_hashmap32_t type_map = { .allocator = temp_arena };

	uint64_t prev_res_key = 0;
    size_t num_unassigned_particle_types = 0;

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
			seq_id = sys->residue.seq_id ? sys->residue.seq_id[sys->residue.count - 1] + 1 : 1;
		}
        char chain_id = data->atom_coordinates[i].chain_id;
        md_atomic_number_t atomic_number = 0;
        md_particle_type_idx_t particle_type_idx = 0;

        if (data->atom_coordinates[i].element[0]) {
            // If the element is available, use that to lookup the element
            str_t sym = str_from_cstrn(data->atom_coordinates[i].element, sizeof(data->atom_coordinates[i].element));
            atomic_number = md_atomic_number_from_symbol(sym, true);
            particle_type_idx = md_util_particle_type_find_or_add(sys, &type_map, atom_id, atomic_number);
        } else {
            num_unassigned_particle_types += 1;
        }

        // A residue is its name, its sequence number and its insertion code in its chain: 100, 100A and 100B are three
        // residues of the same number (the insertions of Kabat numbered antibodies), also when they share a name
        const char icode = data->atom_coordinates[i].icode;
        uint64_t res_key = md_hash64_str(res_name, md_hash64(&icode, sizeof(icode), md_hash64(&seq_id, sizeof(seq_id), chain_id)));

        if (res_key != prev_res_key || terminator) {
            // New residue. A TER record ends the residue before it, whatever the residue after it is called. Whether it
            // also ends a chain is left to the chain ids and the bonds: TER is written inconsistently, and a chain
            // continues over gaps anyway (see md_util_system_infer_entity_and_instance).
            sys->residue.count += 1;
            md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, sys->alloc);
            md_array_push(sys->residue.name,     make_label(res_name), sys->alloc);
            md_array_push(sys->residue.seq_id, seq_id, sys->alloc);
            md_array_push(sys->residue.flags,    MD_RESIDUE_FLAG_NONE, sys->alloc);

            str_t asym_id = str_trim(str_from_cstrn(&data->atom_coordinates[i].chain_id, 1));
            md_array_push(res_auth_asym_ids, asym_id, temp_arena);
        }

        sys->atom.count += 1;

        md_array_push_no_grow(atom_name, atom_id);
        md_array_push_no_grow(state->xyz, vec3_set(x, y, z));
        md_array_push_no_grow(sys->atom.flags, MD_ATOM_FLAG_NONE);
        md_array_push_no_grow(sys->atom.type_idx, particle_type_idx);
        md_array_push_no_grow(atom_occupancy, data->atom_coordinates[i].occupancy);
        md_array_push_no_grow(atom_b_factor,  data->atom_coordinates[i].temp_factor);

		prev_res_key = res_key;

        terminator = (data->atom_coordinates[i].flags & MD_PDB_COORD_FLAG_TERMINATOR) != 0;
    }
    md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, sys->alloc);    // Final sentinel

    ASSERT(md_array_size(state->xyz) == sys->atom.count);
    state->num_atoms = sys->atom.count;

    // Occupancy is a dimensionless fraction; the b factor is a mean square displacement. The paths
    // and the units are md_system's to define, not this loader's - md_mmcif publishes the same two
    // out of the same file family and they have to agree.
    md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/occupancy"), md_unit_none(),                     1, atom_occupancy, sys->atom.count);
    md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/b_factor"),  md_unit_pow(md_unit_angstrom(), 2), 1, atom_b_factor,  sys->atom.count);

    if (data->num_cryst1 > 0) {
        // Use first crystal
        const md_pdb_cryst1_t* cryst = &data->cryst1[0];
        if (cryst->a == 1 && cryst->b == 1 && cryst->c == 1 && cryst->alpha == 90 && cryst->beta == 90 && cryst->gamma == 90) {
            // This is a special case:
            // This is the identity matrix, and in such case, we assume there is no unit cell (no periodic boundary conditions)
        } else {
            state->unitcell = md_unitcell_from_extent_and_angles(cryst->a, cryst->b, cryst->c, cryst->alpha, cryst->beta, cryst->gamma);
        }
    };
    
    if (num_unassigned_particle_types > 0) {
        md_util_system_infer_particle_types(sys, atom_name);
    }

    // The bonds of the CONECT records, of origin MD_BOND_ORIGIN_FILE: inference keeps them and adds the bonds they leave
    // out (CONECT records usually cover the ligands only, in a coarse grained structure from martinize2 all of it).
    // Their serials name atoms of the first model; alternate locations other than the first are not loaded.
    if (data->num_connections > 0) {
        md_hashmap32_t serial_to_atom = { .allocator = temp_arena };
        md_hashmap_reserve(&serial_to_atom, sys->atom.count);
        size_t atom_idx = 0;
        for (size_t i = beg_atom_index; i < end_atom_index; ++i) {
            const char alt_loc = data->atom_coordinates[i].alt_loc;
            if (alt_loc != ' ' && alt_loc != 'A') continue;
            const uint64_t key = (uint64_t)(uint32_t)data->atom_coordinates[i].atom_serial;
            // 0: a serial which is not a number of the format (parse_id), which no record can refer to
            if (key && !md_hashmap_get(&serial_to_atom, key)) md_hashmap_add(&serial_to_atom, key, (uint32_t)atom_idx);
            atom_idx += 1;
        }
        ASSERT(atom_idx == sys->atom.count);

        // Each bond is in the records of both its atoms, and some writers repeat it for its bond order
        md_hashmap32_t seen = { .allocator = temp_arena };
        md_hashmap_reserve(&seen, data->num_connections * 4);
        const md_bond_flags_t file = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_FILE);
        for (size_t i = 0; i < data->num_connections; ++i) {
            const int32_t* serial = data->connections[i].atom_serial;
            if (!serial[0]) continue;
            const uint32_t* a = md_hashmap_get(&serial_to_atom, (uint64_t)(uint32_t)serial[0]);
            if (!a) continue;
            for (size_t k = 1; k < ARRAY_SIZE(data->connections[i].atom_serial) && serial[k]; ++k) {
                const uint32_t* b = md_hashmap_get(&serial_to_atom, (uint64_t)(uint32_t)serial[k]);
                if (!b || *a == *b) continue;
                const uint32_t lo = MIN(*a, *b), hi = MAX(*a, *b);
                const uint64_t pair_key = ((uint64_t)lo << 32) | hi;
                if (md_hashmap_get(&seen, pair_key)) continue;
                md_hashmap_add(&seen, pair_key, 1);
                md_atom_pair_t pair = {{ (md_atom_idx_t)lo, (md_atom_idx_t)hi }};
                md_array_push(sys->bond.pairs, pair, sys->alloc);
                md_array_push(sys->bond.flags, file, sys->alloc);
                sys->bond.count += 1;
            }
        }
    }

    // NOTE: infer_res_flags below consumes bond connectivity, so this cannot simply move
    // out into md_util_system_infer without moving res_flags with it.
    md_util_system_infer_covalent_bonds(sys, state);
    md_util_system_infer_residue_flags(sys);
    md_util_system_infer_entity_and_instance(sys, res_auth_asym_ids);

    // The assemblies name the chains by the ids they have in the file, which are those of the instances
    pdb_assemblies(sys, data, temp_arena);

    result = true;
    md_temp_end(temp_scope);
    return result;
}

bool md_pdb_system_init_from_str(md_system_t* sys, md_system_state_t* state, str_t str, md_pdb_options_t options) {
    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    md_pdb_data_t data = {0};
    md_pdb_data_parse_str(&data, str, temp_arena);
    bool success = md_pdb_system_init_from_data(sys, state, &data, options);
    
    md_temp_end(temp_scope);
    return success;
}

bool md_pdb_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename, md_pdb_options_t options) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Could not open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
    const size_t cap = MEGABYTES(1);
    char *buf = md_temp_alloc(temp_scope, cap);
    ASSERT(buf);
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
    
    md_pdb_data_t data = {0};
    bool success = pdb_parse(&data, &reader, temp_arena, false) && md_pdb_system_init_from_data(sys, state, &data, options);

    // Several models are a trajectory as well: md_pdb_system_publish_run makes them a run.

    md_temp_end(temp_scope);

    return success;
}

typedef struct pdb_cache_t {
    md_run_cache_header_t header;
    md_unitcell_t cell;
    int64_t* frame_offsets;
} pdb_cache_t;

// The cache beside path when it was made from the file as it is now: header, cell, num_frames + 1
// model offsets.
static bool try_read_cache(pdb_cache_t* cache, str_t path, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);
    md_file_t file = {0};
    if (!md_run_cache_open(&file, &cache->header, path, MD_PDB_CACHE_MAGIC, MD_PDB_CACHE_VERSION)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    cache->frame_offsets = md_alloc(alloc, offset_bytes);
    const bool ok =
        md_file_read(file, &cache->cell, sizeof(cache->cell)) == sizeof(cache->cell) &&
        md_file_read(file, cache->frame_offsets, offset_bytes) == offset_bytes &&
        md_file_tell(file) == (int64_t)md_file_size(file);
    if (!ok) {
        MD_LOG_ERROR("The PDB cache beside '" STR_FMT "' is incomplete", STR_ARG(path));
        md_free(alloc, cache->frame_offsets, offset_bytes);
        cache->frame_offsets = NULL;
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const pdb_cache_t* cache, str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    if (!md_run_cache_create(&file, path, scanned, MD_PDB_CACHE_MAGIC, MD_PDB_CACHE_VERSION, cache->header.num_atoms, cache->header.num_frames)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    const bool ok =
        md_file_write(file, &cache->cell, sizeof(cache->cell)) == sizeof(cache->cell) &&
        md_file_write(file, cache->frame_offsets, offset_bytes) == offset_bytes;
    if (!ok) {
        MD_LOG_ERROR("Failed to write the PDB cache beside '" STR_FMT "'", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// Where each model is and the file's cell, from the cache beside the file when that is current and
// from a parse of the file otherwise (writing the cache unless told not to). From alloc; the parse's
// scratch keeps clear of avoid, which is what alloc was made from.
static bool pdb_index_load(pdb_cache_t* cache, str_t filename, md_run_flags_t flags, md_allocator_i* alloc, md_allocator_i* avoid) {
    MEMSET(cache, 0, sizeof(*cache));

    if (try_read_cache(cache, filename, alloc)) {
        return true;
    }
    MEMSET(cache, 0, sizeof(*cache));

    // The file as it is before the parse is what the cache is stamped with
    md_file_info_t scanned = { 0 };
    if (!md_file_info_extract_from_path(filename, &scanned)) {
        MD_LOG_ERROR("Failed to extract file info from path '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    const size_t filesize = scanned.size;

    md_temp_scope_t temp_scope = md_temp_begin_avoid(avoid);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    bool success = false;

    md_pdb_data_t data = {0};
    if (!md_pdb_data_parse_file(&data, filename, temp_alloc)) {
        goto done;
    }

    if (data.num_models <= 1) {
        MD_LOG_INFO("The PDB file does not contain multiple model entries and cannot be read as a trajectory");
        goto done;
    }

    // Validate the models
    const size_t num_atoms = (size_t)MAX(0, data.models[0].end_atom_index - data.models[0].beg_atom_index);
    if (!num_atoms) {
        MD_LOG_ERROR("The PDB file models are empty and cannot be read as a trajectory");
        goto done;
    }
    for (size_t i = 1; i < data.num_models; ++i) {
        const size_t length = (size_t)MAX(0, data.models[i].end_atom_index - data.models[i].beg_atom_index);
        if (length && length != num_atoms) {
            MD_LOG_ERROR("The PDB file models are empty or not of equal length and cannot be read as a trajectory");
            goto done;
        }
    }

    cache->header.num_atoms = num_atoms;
    cache->header.num_frames = data.num_models;

    cache->frame_offsets = md_alloc(alloc, (cache->header.num_frames + 1) * sizeof(int64_t));
    for (size_t i = 0; i < cache->header.num_frames; ++i) {
        cache->frame_offsets[i] = data.models[i].byte_offset;
    }
    cache->frame_offsets[cache->header.num_frames] = filesize;

    if (data.num_cryst1 > 0) {
        if (data.num_cryst1 > 1) {
            md_log(MD_LOG_TYPE_INFO, "The PDB file contains multiple CRYST1 entries, will pick the first one for determining the simulation box");
        }
        // If it is in fact a box, that will be handled as well
        cache->cell = md_unitcell_from_extent_and_angles(data.cryst1[0].a, data.cryst1[0].b, data.cryst1[0].c, data.cryst1[0].alpha, data.cryst1[0].beta, data.cryst1[0].gamma);
    }

    if (!(flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        // A cache that cannot be written only costs the next load a parse
        write_cache(cache, filename, &scanned);
    }
    success = true;

done:
    md_temp_end(temp_scope);
    return success;
}

// ### RUN ###

// The coordinates of atoms [first, first + count) of one model, packed xyz, straight off its ATOM and
// HETATM records. Returns how many were written.
static size_t pdb_parse_model_coords(float* xyz, str_t str, size_t first, size_t count) {
    str_t line;
    // The model's leading MODEL record
    str_extract_line(&line, &str);

    size_t i = 0;
    size_t written = 0;
    while (written < count && str_extract_line(&line, &str)) {
        if (line.len < 6) continue;
        if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
            if (i >= first) {
                xyz[written * 3 + 0] = extract_float(line, 31, 38);
                xyz[written * 3 + 1] = extract_float(line, 39, 46);
                xyz[written * 3 + 2] = extract_float(line, 47, 54);
                written += 1;
            }
            i += 1;
        }
    }
    return written;
}

// <run>/atom/position: one read of the model's text, parsed into the caller's buffer. One atom asks
// for the text up to that atom only.
static size_t pdb_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    ASSERT(sys);
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;

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
        if (pdb_parse_model_coords((float*)dst, (str_t){ text, frame_size }, first, count) == count) {
            written = cap;
        } else {
            MD_LOG_ERROR("PDB: model %u of '" STR_FMT "' has fewer than %zu atoms", frame, STR_ARG(src.path), first + count);
        }
    } else {
        MD_LOG_ERROR("PDB: Failed to read model %u from '" STR_FMT "'", frame, STR_ARG(src.path));
    }
    md_temp_end(temp);
    return written;
}

bool md_pdb_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const str_t path = { path_buf, path_len };

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;

    pdb_cache_t index;
    if (path_len == 0 || !pdb_index_load(&index, path, flags, arena, md_get_heap_allocator())) {
        goto done;
    }

    const size_t F = index.header.num_frames;
    const size_t N = index.header.num_atoms;
    double*  times = md_alloc(arena, F * sizeof(double));
    int64_t* sizes = md_alloc(arena, F * sizeof(int64_t));
    float*   boxes = md_alloc(arena, F * 9 * sizeof(float));

    // One CRYST1 for the file, so one cell for every model; zero without one
    float A[3][3] = {0};
    if (index.cell.flags) {
        md_unitcell_A_extract_float(A, &index.cell);
    }
    for (size_t i = 0; i < F; ++i) {
        times[i] = (double)i;
        sizes[i] = index.frame_offsets[i + 1] - index.frame_offsets[i];
        MEMCPY(boxes + i * 9, A, sizeof(A));
    }

    const md_attribute_virtual_t virt = { .provider = pdb_position_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = times,
        .time_unit     = md_unit_none(),    // a model has no time: ordinals
        .unitcell      = boxes,
        .source_path   = path,
        .source_offset = index.frame_offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    result = md_run_publish(sys, run, &desc);

done:
    md_arena_allocator_destroy(arena);
    return result;
}

bool md_pdb_system_write_state_to_file(md_file_t file, const struct md_system_t* sys, const md_system_state_t* state, const int32_t* atom_indices, size_t num_atoms, int model_num) {
    if (model_num > 0) {
        md_file_printf(file, "MODEL     %4d\n", model_num);
    }

	md_instance_idx_t prev_inst_idx = INT_MIN;
    for (size_t i = 0; i < num_atoms; ++i) {
        int idx = atom_indices ? atom_indices[i] : (int)i;

        char atom_id[6] = "";
        {
            const char* errmsg = hy36encode(5, idx + 1, atom_id);
            if (errmsg) {
                MD_LOG_DEBUG("Failed to encode atom index '%d', error: %s. This atom will be written with a default atom id of 1", idx + 1, errmsg);
            }
        }
        
        // Get element symbol
		md_atomic_number_t atomic_nr = md_atom_atomic_number(&sys->atom, idx);

		char element[3] = "";
        str_copy_to_char_buf(element, sizeof(element), md_atomic_number_symbol(atomic_nr));
        
        // Get atom name (use element symbol if not available)
		char name[5] = "";
		str_copy_to_char_buf(name, sizeof(name), md_atom_name(&sys->atom, idx));
        
        // Get residue name
		md_residue_idx_t res_idx = md_residue_find_by_atom_idx(&sys->residue, idx);
		char resname[5] = "";
		str_copy_to_char_buf(resname, sizeof(resname), md_residue_name(&sys->residue, res_idx));

		md_sequence_id_t res_seq = md_residue_seq_id(&sys->residue, res_idx);
        
        // Get chain ID from instance
        char chain_id[4] = " ";
        md_instance_idx_t inst_idx = md_system_instance_find_by_atom_idx(sys, idx);
        str_copy_to_char_buf(chain_id, sizeof(chain_id), md_instance_auth_id(&sys->instance, inst_idx));

        char res_seq_str[5] = "";
        {
            const char* errmsg = hy36encode(4, res_seq, res_seq_str);
            if (errmsg) {
                MD_LOG_DEBUG("Failed to encode residue sequence id '%d' for atom index %d, error: %s.", res_seq, idx, errmsg);
            }
        }

        // PDB format:
        // ATOM serial name altLoc resName chainID resSeq iCode x y z occupancy tempFactor element charge
        md_file_printf(file,
            "%-6s%5s %-4s%1c%3s %1c%4s%1c   "
            "%8.3f%8.3f%8.3f"
            "%6.2f%6.2f          "
            "%2s%2s\n",
            "ATOM",
            atom_id,            // serial number
            name,               // atom name
            ' ',                // altLoc
            resname,            // residue name
            chain_id[0],        // chain ID
            res_seq_str,        // residue sequence number
            ' ',                // iCode
            state->xyz[idx].x,      // x coordinate
            state->xyz[idx].y,      // y coordinate
            state->xyz[idx].z,      // z coordinate
            1.0,                // occupancy
            0.0,                // bfactor
            element,            // element symbol
            ""                  // charge
        );

		if (prev_inst_idx != INT_MIN && inst_idx != prev_inst_idx) {
			md_file_printf(file, "TER\n");
		}
        prev_inst_idx = inst_idx;
    }
    
    if (model_num > 0) {
        md_file_printf(file, "ENDMDL\n");
    } else {
        md_file_printf(file, "END\n");
    }

    return true;
}

#ifdef __cplusplus
}
#endif
