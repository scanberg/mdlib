#include <md_gro.h>

#include <md_util.h>

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_array.h>
#include <core/md_parse.h>
#include <core/md_hash.h>
#include <core/md_executor.h>

#include <core/md_reader.h>

#include <string.h>

#define MD_GRO_MOL_MAGIC  0xbabcebf677bfaf67
#define NM_TO_ANGSTROM 10.0f

typedef struct gro_molecule {
    uint64_t magic;
    struct md_allocator_i* allocator;
} gro_molecule_t;

// This is a specialized version of tokenization that only extracts float tokens.
// It is designed to avoid the issue where floats are not separated by whitespace.
static size_t extract_float_tokens(str_t* tok_arr, size_t tok_cap, str_t str) {
    size_t n = 0;
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;
    
    while (n < tok_cap && c < end) {
        // Find beg of float token
        while (c < end && !(is_digit(*c) || *c == '-'))
            ++c;

        if (c == end) {
            break;
        }
        const char* beg = c;
        ++c; // Skip the first found character (may be -, then we don't have to check against that again).
        
        // Find end of float token
        while (c < end && (is_digit(*c) || *c == '.'))
            ++c;
        if (beg != c) {
            tok_arr[n++] = (str_t) {beg, c-beg};
        }
    }
    
    return n;
}

// The columns of the coordinates, as GROMACS writes them: three fields of 'width' characters from
// column 20, with 'decimals' decimals - "%8.3f" - and the velocities in three more with one decimal
// more - "%8.4f" - both wider by as many decimals as were asked for beyond that. As GROMACS reads them,
// the width is the distance between the first two decimal points of an atom line, and the decimals
// five fewer: the point is always the fifth character of a field. False when the line does not show
// such columns.
static bool gro_coordinate_columns(size_t* width, size_t* decimals, str_t line) {
    const str_t coords = str_substr(line, 20, SIZE_MAX);
    size_t dot[3];
    size_t from = 0;
    for (int i = 0; i < 3; ++i) {
        size_t loc;
        if (!str_find_char(&loc, str_substr(coords, from, SIZE_MAX), '.')) return false;
        dot[i] = from + loc;
        from = dot[i] + 1;
    }
    const size_t w = dot[1] - dot[0];
    if (dot[0] != 4 || w < 6 || dot[2] - dot[1] != w) return false;
    *width = w;
    *decimals = w - 5;
    return true;
}

// The coordinates of an atom line, in nm, and its velocities in nm/ps (NAN when the line has none). A line that keeps
// to the columns (when there are, as the first atom line shows) is read field by field without a scan; one that does
// not is read by its tokens. False when it was read by position, as it had too few tokens.
static bool gro_parse_coordinates(float xyz[3], float vel[3], str_t line, bool columns, size_t width, size_t decimals) {
    vel[0] = vel[1] = vel[2] = NAN;
    if (columns && line.len >= 20 + 3 * width) {
        const char* f = line.ptr + 20;
        bool read = md_parse_fixed_f32(&xyz[0], (str_t){ f,             width }, decimals) &&
                    md_parse_fixed_f32(&xyz[1], (str_t){ f + width,     width }, decimals) &&
                    md_parse_fixed_f32(&xyz[2], (str_t){ f + 2 * width, width }, decimals);
        // Then nothing, or the three velocities
        if (read && !str_empty(str_trim(str_substr(line, 20 + 3 * width, SIZE_MAX)))) {
            read = line.len >= 20 + 6 * width &&
                   md_parse_fixed_f32(&vel[0], (str_t){ f + 3 * width, width }, decimals + 1) &&
                   md_parse_fixed_f32(&vel[1], (str_t){ f + 4 * width, width }, decimals + 1) &&
                   md_parse_fixed_f32(&vel[2], (str_t){ f + 5 * width, width }, decimals + 1);
            if (!read) {
                vel[0] = vel[1] = vel[2] = NAN;
            }
        }
        if (read) return true;
    }

    str_t tokens[16];
    const int64_t num_tokens = extract_float_tokens(tokens, ARRAY_SIZE(tokens), str_substr(line, 20, SIZE_MAX));
    if (num_tokens < 3) {
        // Fallback to fixed width format
        xyz[0] = (float)parse_float(str_substr(line, 20, 8));
        xyz[1] = (float)parse_float(str_substr(line, 28, 8));
        xyz[2] = (float)parse_float(str_substr(line, 36, 8));
        return false;
    }
    xyz[0] = (float)parse_float(tokens[0]);
    xyz[1] = (float)parse_float(tokens[1]);
    xyz[2] = (float)parse_float(tokens[2]);

    // Velocities are all three or none: a line carrying only part of them is malformed, and
    // guessing which residue was meant is worse than reading none.
    if (num_tokens >= 6) {
        vel[0] = (float)parse_float(tokens[3]);
        vel[1] = (float)parse_float(tokens[4]);
        vel[2] = (float)parse_float(tokens[5]);
    }
    return true;
}

static bool md_gro_data_parse(md_gro_data_t* data, md_reader_i* reader, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(reader);
    ASSERT(alloc);
    str_t line;
    str_t tokens[16];

    if (!md_reader_extract_line(reader, &line)) {
        MD_LOG_ERROR("Failed to parse gro title");
        return false;
    }
    str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

    if (!md_reader_extract_line(reader, &line)) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        return false;
    }
    
    data->num_atoms = parse_int(str_trim(line));
    if (!data->num_atoms) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        return false;
    }

    md_array_resize(data->atom_data, data->num_atoms, alloc);

	bool warn_about_fixed_width = false;

    // The columns as the first atom line shows them. A line that keeps to them is read field by
    // field without a scan; one that does not is read by its tokens.
    size_t width = 0, decimals = 0;
    bool columns = false;

    for (size_t i = 0; i < data->num_atoms; ++i) {
        if (!md_reader_extract_line(reader, &line)) {
            MD_LOG_ERROR("Failed to extract atom line");
            return false;
        }
        md_gro_atom_t* atom = &data->atom_data[i];

        if (i == 0) {
            columns = gro_coordinate_columns(&width, &decimals, line);
        }

        float xyz[3], vel[3];
        if (!gro_parse_coordinates(xyz, vel, line, columns, width, decimals)) {
            warn_about_fixed_width = true;
        }
        atom->x  = xyz[0]; atom->y  = xyz[1]; atom->z  = xyz[2];
        atom->vx = vel[0]; atom->vy = vel[1]; atom->vz = vel[2];

        atom->res_id = (int32_t)parse_int(str_trim(str_substr(line, 0, 5)));
        
        str_copy_to_char_buf(atom->res_name,  sizeof(atom->res_name),  str_trim(str_substr(line,  5, 5)));
        str_copy_to_char_buf(atom->atom_name, sizeof(atom->atom_name), str_trim(str_substr(line, 10, 5)));
    }

    if (warn_about_fixed_width) {
        MD_LOG_DEBUG("Parsed GRO file using fixed-width format. This may indicate an issue with the file formatting.");
	}

    if (!md_reader_extract_line(reader, &line)) {
        MD_LOG_ERROR("Failed to extract unitcell line");
        return false;
    }

    const int64_t num_tokens = extract_float_tokens(tokens, ARRAY_SIZE(tokens), line);

    if (num_tokens != 3 && num_tokens != 9) {
        MD_LOG_ERROR("Failed to parse box, expected 3 or 9 tokens");
        return false;
    }
    
    data->box[0][0] = (float)parse_float(tokens[0]);
    data->box[1][1] = (float)parse_float(tokens[1]);
    data->box[2][2] = (float)parse_float(tokens[2]);

    if (num_tokens > 3) {
        data->box[0][1] = (float)parse_float(tokens[3]);
        data->box[0][2] = (float)parse_float(tokens[4]);
        data->box[1][0] = (float)parse_float(tokens[5]);
        data->box[1][2] = (float)parse_float(tokens[6]);
        data->box[2][0] = (float)parse_float(tokens[7]);
        data->box[2][1] = (float)parse_float(tokens[8]);
    }
    
    return true;
}

bool md_gro_data_parse_str(md_gro_data_t* data, str_t str, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);

    md_reader_i line_reader = md_reader_from_str(str);
    return md_gro_data_parse(data, &line_reader, alloc);
}

bool md_gro_data_parse_file(md_gro_data_t* data, str_t filename, struct md_allocator_i* alloc) {
    bool result = false;
    md_file_t file = {0};
    if (md_file_open(&file, filename, MD_FILE_READ)) {
        const int64_t cap = MEGABYTES(1);
        md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
        
        md_file_reader_t line_reader;
        md_file_reader_init(&line_reader, file, cap, md_temp_allocator(temp_scope));
        result = md_gro_data_parse(data, &line_reader.r, alloc);
        
        md_file_reader_close(&line_reader);
        md_temp_end(temp_scope);
        md_file_close(&file);
    } else {
        MD_LOG_ERROR("Could not open file '%.*s'", filename.len, filename.ptr);
    }
    return result;
}

void md_gro_data_free(md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    if (data->atom_data) md_array_free(data->atom_data, alloc);
    MEMSET(data, 0, sizeof(md_gro_data_t));
}

// ---------------------------------------------------------------------------------------------------------------------
// THE SYSTEM
//
// What a load gives is the core: the atoms, where they are, their types, the residues and the cell. Nothing is derived
// from the coordinates here, neither the bonds nor what the residues are by them: that is md_util_system_infer's, which
// a caller runs when (and where, on what thread) it wants.

static bool gro_system_begin(md_system_t* sys, md_system_state_t* state, size_t num_atoms) {
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state || !state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }
    md_system_reset(sys);
    if (!md_system_state_init(state, num_atoms)) return false;

    const size_t capacity = ROUND_UP(num_atoms, 16);
    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags, capacity, sys->alloc);
    MEMSET(sys->atom.type_idx, 0, capacity * sizeof(md_particle_type_idx_t));
    MEMSET(sys->atom.flags, 0, capacity * sizeof(md_atom_flags_t));
    sys->atom.count = num_atoms;

    sys->atom.type.count = 0;
    md_particle_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
    return true;
}

static inline void gro_push_residue(md_system_t* sys, uint32_t atom_beg, md_label_t name, md_sequence_id_t seq_id) {
    sys->residue.count += 1;
    md_array_push(sys->residue.atom_offset, atom_beg, sys->alloc);
    md_array_push(sys->residue.name, name, sys->alloc);
    md_array_push(sys->residue.seq_id, seq_id, sys->alloc);
    md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, sys->alloc);
}

// The velocities, in the file's own nm/ps rather than converted the way the coordinates are, because a velocity has a
// unit and the attribute carries it: a consumer wanting Angstrom per picosecond asks the extract for it and gets a
// checked conversion, which is strictly better than a silent one here that no one downstream can see happened.
static void gro_publish_velocities(md_system_t* sys, const float* velocity, size_t num_atoms) {
    const md_unit_t nm_per_ps = md_unit_div(md_unit_nanometer(), md_unit_picosecond());
    md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/velocity"), nm_per_ps, 3, velocity, num_atoms);
}

static bool gro_parse_box(md_unitcell_t* cell, str_t line) {
    str_t tokens[16];
    const int64_t num_tokens = extract_float_tokens(tokens, ARRAY_SIZE(tokens), line);
    if (num_tokens != 3 && num_tokens != 9) {
        MD_LOG_ERROR("Failed to parse box, expected 3 or 9 tokens");
        return false;
    }
    float box[3][3] = {0};
    box[0][0] = (float)parse_float(tokens[0]);
    box[1][1] = (float)parse_float(tokens[1]);
    box[2][2] = (float)parse_float(tokens[2]);
    if (num_tokens > 3) {
        box[0][1] = (float)parse_float(tokens[3]);
        box[0][2] = (float)parse_float(tokens[4]);
        box[1][0] = (float)parse_float(tokens[5]);
        box[1][2] = (float)parse_float(tokens[6]);
        box[2][0] = (float)parse_float(tokens[7]);
        box[2][1] = (float)parse_float(tokens[8]);
    }
    // nm to Ångström
    for (int i = 0; i < 3; ++i) {
        box[i][0] *= NM_TO_ANGSTROM;
        box[i][1] *= NM_TO_ANGSTROM;
        box[i][2] *= NM_TO_ANGSTROM;
    }
    *cell = md_unitcell_from_matrix_float(MD_AS_CONST_MAT3(box));
    return true;
}

bool md_gro_system_init_from_data(struct md_system_t* sys, md_system_state_t* state, const md_gro_data_t* data) {
    ASSERT(sys);
    ASSERT(data);
    if (!gro_system_begin(sys, state, data->num_atoms)) return false;

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_label_t* atom_names = md_temp_alloc_array(temp, md_label_t, data->num_atoms);

	uint64_t prev_res_key = 0;
    for (size_t i = 0; i < data->num_atoms; ++i) {
        const md_gro_atom_t* atom = &data->atom_data[i];
        str_t res_name = str_from_cstrn(atom->res_name, sizeof(atom->res_name));
		uint64_t res_key = md_hash64_str(res_name, (uint64_t)atom->res_id);
        if (i == 0 || res_key != prev_res_key) {
            gro_push_residue(sys, (uint32_t)i, make_label(res_name), atom->res_id);
		}
		prev_res_key = res_key;
        atom_names[i] = make_label(str_from_cstrn(atom->atom_name, sizeof(atom->atom_name)));
        state->xyz[i] = vec3_set(atom->x * NM_TO_ANGSTROM, atom->y * NM_TO_ANGSTROM, atom->z * NM_TO_ANGSTROM);
    }
	md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, sys->alloc);    // Final sentinel

    float* velocity = md_temp_alloc_array(temp, float, data->num_atoms * 3);
    if (velocity) {
        for (size_t i = 0; i < data->num_atoms; ++i) {
            velocity[i * 3 + 0] = data->atom_data[i].vx;
            velocity[i * 3 + 1] = data->atom_data[i].vy;
            velocity[i * 3 + 2] = data->atom_data[i].vz;
        }
        gro_publish_velocities(sys, velocity, data->num_atoms);
    }

    float box[3][3] = {0};
    for (int i = 0; i < 3; ++i) {
        box[i][0] = data->box[i][0] * NM_TO_ANGSTROM;
        box[i][1] = data->box[i][1] * NM_TO_ANGSTROM;
        box[i][2] = data->box[i][2] * NM_TO_ANGSTROM;
    }
    state->unitcell = md_unitcell_from_matrix_float(MD_AS_CONST_MAT3(box));

    md_util_system_infer_particle_types_from_labels(sys, atom_names);

    md_temp_end(temp);
    state->num_atoms = sys->atom.count;
    return true;
}

// ---------------------------------------------------------------------------------------------------------------------
// Straight from the text: the lines of the atoms are read a block at a time, each block in pieces at line boundaries,
// on as many threads as the executor has. A piece counts its lines first, which places its atoms; then it reads them
// into the system's own arrays where they go. The residues are where the residue fields of a line differ from those of
// the line before: each piece keeps the ones it finds, and they are appended in order once the pieces are done.

#define GRO_RES_FIELDS 10   // The residue number and name: columns 0 to 10

// NAN, by its bits: mdlib builds with fast math
static inline bool gro_absent(float v) {
    uint32_t u;
    MEMCPY(&u, &v, sizeof(u));
    return (u & 0x7FFFFFFFu) > 0x7F800000u;
}

typedef struct gro_residue_start_t {
    uint32_t         atom;
    md_sequence_id_t seq_id;
    md_label_t       name;
} gro_residue_start_t;

typedef struct gro_piece_t {
    md_text_piece_t range;      // of the block
    size_t   num_lines;
    size_t   first_line;        // of the block
    size_t   num_residues;      // found, at res_starts + first_line
    size_t   end_of_atoms;      // the byte after the last atom line, if the atoms end in this piece
    bool     atoms_end;
    bool     loose;             // a line was read by position
} gro_piece_t;

typedef struct gro_block_t {
    const char*           text;
    size_t                len;
    gro_piece_t*          pieces;
    size_t                atom_beg;         // of the block's first line
    size_t                num_atoms;        // of the file: the lines after them are not atoms
    str_t                 prev_res;         // the residue fields of the line before the block, empty if none
    gro_residue_start_t*  res_starts;       // [lines of the block]
    bool                  columns;
    size_t                width, decimals;
    vec3_t*               xyz;
    float*                velocity;         // NULL when the file has none
    md_label_t*           atom_names;
} gro_block_t;

static void gro_count_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    gro_block_t* b = (gro_block_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        gro_piece_t* piece = &b->pieces[p];
        size_t n = md_text_count_newlines(b->text + piece->range.beg, b->text + piece->range.end);
        // The last line of the text may not end in one
        if (piece->range.end == b->len && piece->range.end > piece->range.beg && b->text[piece->range.end - 1] != '\n') ++n;
        piece->num_lines = n;
    }
}

static inline bool gro_same_residue(str_t a, str_t b) {
    const size_t na = MIN(a.len, GRO_RES_FIELDS);
    const size_t nb = MIN(b.len, GRO_RES_FIELDS);
    return na == nb && MEMCMP(a.ptr, b.ptr, na) == 0;
}

static void gro_parse_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    gro_block_t* b = (gro_block_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        gro_piece_t* piece = &b->pieces[p];
        piece->num_residues = 0;
        piece->atoms_end = false;
        piece->loose = false;

        // The residue fields of the line before: in the piece before, or the block before
        str_t prev = piece->range.beg > 0 ? md_text_line_before(b->text, piece->range.beg) : b->prev_res;

        gro_residue_start_t* res = b->res_starts + piece->first_line;
        md_text_lines_t lines = md_text_lines(b->text + piece->range.beg, b->text + piece->range.end);
        size_t atom = b->atom_beg + piece->first_line;
        str_t line;
        for (;;) {
            if (atom >= b->num_atoms) {
                if (lines.beg < lines.end) {
                    piece->atoms_end = true;
                    piece->end_of_atoms = (size_t)(lines.beg - b->text);
                }
                break;
            }
            if (!md_text_lines_next(&lines, &line)) break;

            if (atom == 0 || !gro_same_residue(prev, line)) {
                res->atom   = (uint32_t)atom;
                res->seq_id = (md_sequence_id_t)parse_int(str_trim(str_substr(line, 0, 5)));
                res->name   = make_label(str_trim(str_substr(line, 5, 5)));
                ++res;
            }
            prev = line;

            float xyz[3], vel[3];
            if (!gro_parse_coordinates(xyz, vel, line, b->columns, b->width, b->decimals)) piece->loose = true;
            b->xyz[atom] = vec3_set(xyz[0] * NM_TO_ANGSTROM, xyz[1] * NM_TO_ANGSTROM, xyz[2] * NM_TO_ANGSTROM);
            if (b->velocity) {
                b->velocity[atom * 3 + 0] = vel[0];
                b->velocity[atom * 3 + 1] = vel[1];
                b->velocity[atom * 3 + 2] = vel[2];
            }
            b->atom_names[atom] = make_label(str_trim(str_substr(line, 10, 5)));
            ++atom;
        }
        if (!piece->atoms_end && atom >= b->num_atoms) {
            piece->atoms_end = true;
            piece->end_of_atoms = piece->range.end;
        }
        piece->num_residues = (size_t)(res - (b->res_starts + piece->first_line));
    }
}

static bool gro_load(md_system_t* sys, md_system_state_t* state, md_reader_i* src, const md_executor_t* ex) {
    ASSERT(sys);
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    // What the load holds besides the system, in a temp scope rewound at the end: the executor's tasks allocate
    // nothing (gro_count_lines, gro_parse_lines), so their threads' scopes and this one do not meet
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* scratch = md_temp_allocator(temp);
    bool result = false;
    size_t num_atoms = 0;
    gro_block_t b = {0};
    md_array(gro_piece_t) pieces = 0;
    md_array(md_text_piece_t) ranges = 0;
    md_array(gro_residue_start_t) res_starts = 0;  // of a block, at the place of the line of each piece
    md_array(gro_residue_start_t) residues = 0;    // of the load
    char prev_res[GRO_RES_FIELDS];
    bool loose = false;
    size_t atom = 0;
    str_t line;

    if (!md_reader_extract_line(src, &line)) {
        MD_LOG_ERROR("Failed to parse gro title");
        goto done;
    }
    if (!md_reader_extract_line(src, &line)) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        goto done;
    }
    const int64_t num_atoms_in_file = parse_int(str_trim(line));
    if (num_atoms_in_file <= 0 || num_atoms_in_file >= UINT32_MAX) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        goto done;
    }
    num_atoms = (size_t)num_atoms_in_file;

    if (!gro_system_begin(sys, state, num_atoms)) goto done;
    // Room for what the passes hold, so that nothing grows (and copies) in the arena: the residues of the load (no
    // more than the atoms), the residues of a pass (no more than its lines, of 20 characters at least)
    md_array_ensure(residues, num_atoms, scratch);
    md_array_ensure(res_starts, MD_TEXT_PASS_SIZE / 20 + 16, scratch);

    // The columns, and whether there are velocities, as the first atom line shows them. GROMACS writes the velocities
    // of all atoms or of none: a file whose first atom has none is read without.
    b.num_atoms = num_atoms;
    b.xyz = state->xyz;
    {
        str_t first = {0};
        md_reader_peek_line(src, &first);
        b.columns = gro_coordinate_columns(&b.width, &b.decimals, first);
        float xyz[3], vel[3];
        gro_parse_coordinates(xyz, vel, first, b.columns, b.width, b.decimals);
        if (!gro_absent(vel[0])) {
            b.velocity = md_alloc(scratch, num_atoms * 3 * sizeof(float));
        }
    }
    b.atom_names = md_alloc(scratch, num_atoms * sizeof(md_label_t));

    const uint32_t num_threads = md_executor_num_threads(ex);

    while (atom < num_atoms) {
        str_t block = md_reader_block(src);
        if (block.len == 0) {
            MD_LOG_ERROR("The file ends after %zu of its %zu atoms", atom, num_atoms);
            goto done;
        }
        // A pass over at most MD_TEXT_PASS_SIZE: the lines counted, then read, of text the count left in the caches
        block.len = md_text_pass_len(block.ptr, block.len);

        const size_t num_pieces = md_text_num_pieces(block.len, num_threads);
        md_array_resize(pieces, num_pieces, scratch);
        md_array_resize(ranges, num_pieces, scratch);
        md_text_split(ranges, num_pieces, block.ptr, block.len);
        for (size_t p = 0; p < num_pieces; ++p) pieces[p] = (gro_piece_t){ .range = ranges[p] };

        b.text = block.ptr;
        b.len = block.len;
        b.pieces = pieces;
        b.atom_beg = atom;

        if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, gro_count_lines, &b)) {
            goto done;  // Cancelled
        }
        size_t num_lines = 0;
        for (size_t p = 0; p < num_pieces; ++p) {
            pieces[p].first_line = num_lines;
            num_lines += pieces[p].num_lines;
        }
        md_array_resize(res_starts, num_lines, scratch);
        b.res_starts = res_starts;

        if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, gro_parse_lines, &b)) {
            goto done;  // Cancelled
        }

        // The residues in order, and where the atoms end
        size_t consumed = block.len;
        for (size_t p = 0; p < num_pieces; ++p) {
            const gro_piece_t* piece = &pieces[p];
            md_array_push_array(residues, res_starts + piece->first_line, piece->num_residues, scratch);
            loose |= piece->loose;
            if (piece->atoms_end) {
                consumed = piece->end_of_atoms;
                break;
            }
        }
        atom += MIN(num_lines, num_atoms - atom);

        // The residue fields of the last atom line, for the next block
        if (atom < num_atoms) {
            const str_t last = md_text_line_before(block.ptr, consumed);
            const size_t n = MIN(last.len, sizeof(prev_res));
            MEMCPY(prev_res, last.ptr, n);
            b.prev_res = (str_t){ prev_res, n };
        }
        md_reader_advance(src, consumed);
    }
    // The residues, at once: arrays which grow in the system's allocator (an arena) leave their copies behind
    {
        const size_t num_res = md_array_size(residues);
        md_array_resize(sys->residue.atom_offset, num_res + 1, sys->alloc);
        md_array_resize(sys->residue.name,   num_res, sys->alloc);
        md_array_resize(sys->residue.seq_id, num_res, sys->alloc);
        md_array_resize(sys->residue.flags,  num_res, sys->alloc);
        for (size_t r = 0; r < num_res; ++r) {
            sys->residue.atom_offset[r] = residues[r].atom;
            sys->residue.name[r]   = residues[r].name;
            sys->residue.seq_id[r] = residues[r].seq_id;
            sys->residue.flags[r]  = MD_RESIDUE_FLAG_NONE;
        }
        sys->residue.atom_offset[num_res] = (uint32_t)num_atoms;    // Final sentinel
        sys->residue.count = num_res;
    }

    if (loose) {
        MD_LOG_DEBUG("Parsed GRO file using fixed-width format. This may indicate an issue with the file formatting.");
    }

    if (!md_reader_extract_line(src, &line)) {
        MD_LOG_ERROR("Failed to extract unitcell line");
        goto done;
    }
    if (!gro_parse_box(&state->unitcell, line)) goto done;

    if (b.velocity) {
        gro_publish_velocities(sys, b.velocity, num_atoms);
    }
    md_util_system_infer_particle_types_from_labels(sys, b.atom_names);
    state->num_atoms = num_atoms;
    result = true;

done:
    // The scratch goes with the scope
    md_temp_end(temp);
    if (!result) md_system_reset(sys);
    return result;
}

bool md_gro_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename, const md_executor_t* ex) {
    // The reader's buffer in a temp scope, below the one of the load (gro_load)
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_file_reader_t fr;
    if (!md_file_reader_open(&fr, filename, MD_TEXT_BLOCK_SIZE, md_temp_allocator(temp))) {
        MD_LOG_ERROR("Could not open file '%.*s'", (int)filename.len, filename.ptr);
        md_temp_end(temp);
        return false;
    }
    const bool result = gro_load(sys, state, &fr.r, ex);
    md_file_reader_close(&fr);
    md_temp_end(temp);
    return result;
}

bool md_gro_system_init_from_str(md_system_t* sys, md_system_state_t* state, str_t str, const md_executor_t* ex) {
    md_reader_i src = md_reader_from_str(str);
    return gro_load(sys, state, &src, ex);
}
