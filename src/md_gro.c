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
#include <core/md_parse_simd.h>
#include <core/md_hash.h>
#include <core/md_executor.h>

#include <core/md_reader.h>

#include <string.h>

#define NM_TO_ANGSTROM 10.0f

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
    const int64_t num_tokens = extract_float_tokens(tokens, MD_COUNTOF(tokens), str_substr(line, 20, SIZE_MAX));
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

// The coordinates and velocities of an atom line in the columns GROMACS writes by default, "%8.3f" and "%8.4f": three
// fields at a time (md_parse_fixed_f32x3_w8), what gro_parse_coordinates reads by those columns. False where it would
// read the line otherwise, or would find something there that is not three velocities: the caller reads the line with it.
static MD_FORCE_INLINE bool gro_parse_coordinates_w8(float xyz[3], float vel[3], str_t line) {
    if (line.len < 20 + 24 || !md_parse_fixed_f32x3_w8(xyz, line.ptr + 20, 3)) return false;
    if (line.len == 20 + 24) {
        vel[0] = vel[1] = vel[2] = NAN;
        return true;
    }
    return line.len >= 20 + 48 && md_parse_fixed_f32x3_w8(vel, line.ptr + 44, 4);
}

// ---------------------------------------------------------------------------------------------------------------------
// THE SYSTEM
//
// What a load gives is the core: the atoms, where they are, their types, the residues and the cell. Nothing is derived
// from the coordinates here, neither the bonds nor what the residues are by them: that is md_util_system_infer's, which
// a caller runs when (and where, on what thread) it wants.
//
// Straight from the text, a block at a time (read over the threads of the executor where the source reads at offsets,
// a file), each block in passes and each pass in pieces at line boundaries, a few for each thread (md_reader.h: blocks
// and pieces). Every line is an atom up to the number of them the file gives, so the line of an atom is its index. Three
// loops over the pieces of a pass, in parallel:
//
//   count  the lines of a piece
//   parse  each atom where it goes (its line), and the residues which begin in the piece, at its first line in the
//          residues of the pass: a residue begins where the residue fields of a line differ from those of the line before
//   place  the residues of each piece where they go in those of the load
//
// Between them the calling thread works with the pieces, never with their atoms or residues.

#define GRO_RES_FIELDS 10   // The residue number and name: columns 0 to 10

// NAN, by its bits: mdlib builds with fast math
static inline bool gro_absent(float v) {
    uint32_t u;
    MD_MEMCPY(&u, &v, sizeof(u));
    return (u & 0x7FFFFFFFu) > 0x7F800000u;
}

// The residue fields of a line (its first GRO_RES_FIELDS characters, fewer of a shorter line) as words, and how many
// there are: two lines are of the same residue where these are the same
typedef struct gro_res_fields_t {
    uint64_t lo;
    uint16_t hi;
    uint16_t len;
} gro_res_fields_t;

static MD_FORCE_INLINE gro_res_fields_t gro_res_fields(str_t line) {
    gro_res_fields_t f = { 0, 0, (uint16_t)MD_MIN(line.len, GRO_RES_FIELDS) };
    if (line.len >= GRO_RES_FIELDS) {
        MD_MEMCPY(&f.lo, line.ptr, 8);
        MD_MEMCPY(&f.hi, line.ptr + 8, 2);
    } else if (line.len) {
        char c[GRO_RES_FIELDS] = {0};
        MD_MEMCPY(c, line.ptr, line.len);
        MD_MEMCPY(&f.lo, c, 8);
        MD_MEMCPY(&f.hi, c + 8, 2);
    }
    return f;
}

static MD_FORCE_INLINE bool gro_res_fields_eq(gro_res_fields_t a, gro_res_fields_t b) {
    return a.lo == b.lo && a.hi == b.hi && a.len == b.len;
}

// The field of 5 characters at offset as a label, without the whitespace around it: make_label(str_trim()) of it. As one
// word where the line goes on for 8 characters from there and the field has no control character (whose blanks are
// spaces alone, which a byte above ' ' bounds: no scan); else as make_label(str_trim()).
static MD_FORCE_INLINE md_label_t gro_label5(str_t line, size_t offset) {
    if (line.len >= offset + 8) {
        uint64_t w;
        MD_MEMCPY(&w, line.ptr + offset, 8);
        const uint64_t field = 0x8080808080ull;                             // the high bits of its 5 bytes
        const uint64_t x = w & 0x7F7F7F7F7F7F7F7Full;
        // 0x80 in each byte at or above ' ' (or with its high bit set): none below it, no control character
        if ((((x + 0x6060606060606060ull) | w) & field) == field) {
            // 0x80 in each byte above ' '
            const uint64_t mark = ((x + 0x5F5F5F5F5F5F5F5Full) | w) & field;
            md_label_t lbl = {0};
            if (mark) {
                const uint32_t lo = md_ctz64(mark) >> 3;
                const uint32_t hi = (63 - md_clz64(mark)) >> 3;
                const uint64_t v = (w >> (8 * lo)) & (~UINT64_C(0) >> (64 - 8 * (hi - lo + 1)));
                MD_MEMCPY(lbl.buf, &v, sizeof(v));
            }
            return lbl;
        }
    }
    return make_label(str_trim(str_substr(line, offset, 5)));
}

// The residue number of a line (columns 0 to 4): parse_int(str_trim()) of it. Blanks and then digits up to its end
// without a scan of what is around; else as parse_int(str_trim()).
static MD_FORCE_INLINE int64_t gro_res_number(str_t line) {
    if (line.len >= 5) {
        const char* c = line.ptr;
        uint32_t i = 0;
        while (i < 5 && c[i] == ' ') ++i;
        int64_t value = 0;
        uint32_t j = i;
        while (j < 5 && is_digit(c[j])) value = value * 10 + (c[j++] - '0');
        if (j == 5) return value;
    }
    return parse_int(str_trim(str_substr(line, 0, 5)));
}

// The velocities, in the file's own nm/ps rather than converted the way the coordinates are, because a velocity has a
// unit and the attribute carries it: a consumer wanting Angstrom per picosecond asks the extract for it and gets a
// checked conversion, which is strictly better than a silent one here that no one downstream can see happened.
static inline md_unit_t gro_velocity_unit(void) {
    return md_unit_div(md_unit_nanometer(), md_unit_picosecond());
}

static bool gro_parse_box(md_unitcell_t* cell, str_t line) {
    str_t tokens[16];
    const int64_t num_tokens = extract_float_tokens(tokens, MD_COUNTOF(tokens), line);
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

typedef struct gro_residue_start_t {
    uint32_t         atom;
    md_sequence_id_t seq_id;
    md_label_t       name;
} gro_residue_start_t;

typedef struct gro_piece_t {
    md_text_piece_t range;      // of the pass
    size_t   num_lines;
    size_t   first_line;        // of the pass
    size_t   num_atoms;         // its lines which are atoms: those before the end of the atoms
    size_t   num_res;           // the residues which begin in it, at pass_res + first_line
    size_t   first_res;         // of the load
    size_t   end_of_atoms;      // the byte after its last atom line, where the atoms end in it
    bool     loose;             // a line was read by position
} gro_piece_t;

// What a load has: the pass it reads, and where what it reads goes
typedef struct gro_load_t {
    const char*           text;             // the pass
    size_t                len;
    gro_piece_t*          pieces;
    size_t                atom_beg;         // the atom of the pass's first line
    size_t                num_atoms;        // of the file: the lines after them are not atoms
    str_t                 prev;             // the residue fields of the atom line before the pass, empty if none
    bool                  columns_w8;       // the columns GROMACS writes by default: three fields of 8 at a time
    gro_residue_start_t*  pass_res;         // [lines of the pass]
    gro_residue_start_t*  res;              // of the load
    bool                  columns;
    size_t                width, decimals;
    vec3_t*               xyz;
    float*                velocity;         // NULL when the file has none
    md_label_t*           atom_names;
    md_particle_type_idx_t* type_idx;       // none yet: typed by the names once all are read
    md_atom_flags_t*      atom_flags;
} gro_load_t;

static void gro_count_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    gro_load_t* L = (gro_load_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        gro_piece_t* piece = &L->pieces[p];
        size_t n = md_text_count_newlines(L->text + piece->range.beg, L->text + piece->range.end);
        // The last line of the text may not end in one
        if (piece->range.end == L->len && piece->range.end > piece->range.beg && L->text[piece->range.end - 1] != '\n') ++n;
        piece->num_lines = n;
    }
}

static void gro_parse_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    gro_load_t* L = (gro_load_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        gro_piece_t* piece = &L->pieces[p];
        piece->num_res = 0;
        piece->loose = false;
        piece->end_of_atoms = piece->range.end;

        // The residue fields of the line before: in the piece before, or the pass before
        gro_res_fields_t prev = gro_res_fields(piece->range.beg > 0 ? md_text_line_before(L->text, piece->range.beg) : L->prev);

        gro_residue_start_t* res = L->pass_res + piece->first_line;
        md_text_lines_t lines = md_text_lines(L->text + piece->range.beg, L->text + piece->range.end);
        size_t atom = L->atom_beg + piece->first_line;
        const size_t atom_end = atom + piece->num_atoms;
        str_t line;
        while (atom < atom_end && md_text_lines_next(&lines, &line)) {
            const gro_res_fields_t fields = gro_res_fields(line);
            if (atom == 0 || !gro_res_fields_eq(prev, fields)) {
                res->atom   = (uint32_t)atom;
                res->seq_id = (md_sequence_id_t)gro_res_number(line);
                res->name   = gro_label5(line, 5);
                ++res;
            }
            prev = fields;

            float xyz[3], vel[3];
            if (!(L->columns_w8 && gro_parse_coordinates_w8(xyz, vel, line)) &&
                !gro_parse_coordinates(xyz, vel, line, L->columns, L->width, L->decimals)) piece->loose = true;
            L->xyz[atom] = vec3_set(xyz[0] * NM_TO_ANGSTROM, xyz[1] * NM_TO_ANGSTROM, xyz[2] * NM_TO_ANGSTROM);
            if (L->velocity) {
                L->velocity[atom * 3 + 0] = vel[0];
                L->velocity[atom * 3 + 1] = vel[1];
                L->velocity[atom * 3 + 2] = vel[2];
            }
            L->atom_names[atom] = gro_label5(line, 10);
            L->type_idx[atom] = 0;
            L->atom_flags[atom] = 0;
            ++atom;
        }
        if (piece->num_atoms < piece->num_lines) {
            piece->end_of_atoms = (size_t)(lines.beg - L->text);
        }
        piece->num_res = (size_t)(res - (L->pass_res + piece->first_line));
    }
}

static void gro_place_residues(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    gro_load_t* L = (gro_load_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        const gro_piece_t* piece = &L->pieces[p];
        MD_MEMCPY(L->res + piece->first_res, L->pass_res + piece->first_line, piece->num_res * sizeof(gro_residue_start_t));
    }
}

// The residues of the system from where they begin, a chunk of GRO_RES_CHUNK at a time
#define GRO_RES_CHUNK 16384

typedef struct gro_residues_t {
    md_residue_data_t* residue;
    const gro_residue_start_t* starts;
    size_t num_res;
} gro_residues_t;

static void gro_build_residues(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const gro_residues_t* rb = (const gro_residues_t*)arg;
    const size_t r_end = MD_MIN((size_t)end * GRO_RES_CHUNK, rb->num_res);
    for (size_t r = (size_t)beg * GRO_RES_CHUNK; r < r_end; ++r) {
        rb->residue->atom_offset[r] = rb->starts[r].atom;
        rb->residue->name[r]   = rb->starts[r].name;
        rb->residue->seq_id[r] = rb->starts[r].seq_id;
        rb->residue->flags[r]  = MD_RESIDUE_FLAG_NONE;
    }
}

// Room for count elements in an array of the load, which grows by doubling past what was reserved
#define GRO_ROOM(arr, count, alloc) do { \
        if ((count) > md_array_capacity(arr)) md_array_ensure((arr), MD_MAX((size_t)(count), 2 * md_array_capacity(arr)), (alloc)); \
        md_array_resize((arr), (count), (alloc)); \
    } while (0)

// A pass of len bytes at text, of the atoms from L->atom_beg: false if the executor was cancelled. Where the atoms end in
// it, consumed is where; else all of it.
static bool gro_read_pass(gro_load_t* L, md_array(gro_piece_t)* pieces, md_array(md_text_piece_t)* ranges, md_array(gro_residue_start_t)* pass_res,
                          size_t* num_res, bool* loose, size_t* consumed, const char* text, size_t len, md_allocator_i* temp, const md_executor_t* ex) {
    const size_t num_pieces = md_text_num_pieces(len, md_executor_num_threads(ex));
    md_array_resize(*pieces, num_pieces, temp);
    md_array_resize(*ranges, num_pieces, temp);
    md_text_split(*ranges, num_pieces, text, len);
    for (size_t p = 0; p < num_pieces; ++p) (*pieces)[p] = (gro_piece_t){ .range = (*ranges)[p] };
    L->text = text;
    L->len = len;
    L->pieces = *pieces;

    if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, gro_count_lines, L)) return false;
    // The lines of each piece, and which of them are atoms
    size_t num_lines = 0;
    for (size_t p = 0; p < num_pieces; ++p) {
        gro_piece_t* piece = &(*pieces)[p];
        piece->first_line = num_lines;
        const size_t first_atom = L->atom_beg + num_lines;
        piece->num_atoms = first_atom < L->num_atoms ? MD_MIN(piece->num_lines, L->num_atoms - first_atom) : 0;
        num_lines += piece->num_lines;
    }
    GRO_ROOM(*pass_res, num_lines, temp);
    L->pass_res = *pass_res;

    if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, gro_parse_lines, L)) return false;
    // Where the residues of each piece go, and where the atoms end
    *consumed = len;
    bool ended = false;
    size_t pass_atoms = 0;
    for (size_t p = 0; p < num_pieces; ++p) {
        gro_piece_t* piece = &(*pieces)[p];
        piece->first_res = *num_res;
        *num_res += piece->num_res;
        *loose |= piece->loose;
        pass_atoms += piece->num_atoms;
        if (!ended && piece->num_atoms < piece->num_lines) {
            *consumed = piece->end_of_atoms;
            ended = true;
        }
    }
    if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, gro_place_residues, L)) return false;
    L->atom_beg += pass_atoms;
    return true;
}

static bool gro_load(md_system_t* sys, md_system_state_t* state, md_reader_i* src, const md_executor_t* ex) {
    MD_ASSERT(sys);
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

    // What the load holds besides the system, in a temp scope rewound at the end: the executor's tasks allocate
    // nothing, so their threads' scopes and this one do not meet
    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp = md_temp_allocator(temp_scope);
    const size_t pass_size = md_text_pass_size(md_executor_num_threads(ex));
    bool result = false;
    bool loose = false;
    size_t num_res = 0;
    gro_load_t L = {0};
    md_array(gro_piece_t) pieces = 0;
    md_array(md_text_piece_t) ranges = 0;
    md_array(gro_residue_start_t) pass_res = 0;
    md_array(gro_residue_start_t) res = 0;
    char prev_res[GRO_RES_FIELDS];
    md_attribute_id_t velocity_id = MD_ATTRIBUTE_INVALID;
    str_t line;

    // The title and the number of atoms, of the first block
    str_t block = md_reader_block_parallel(src, ex);
    const char* c = block.ptr;
    const char* e = block.ptr + block.len;
    if (c == e) {
        MD_LOG_ERROR("Failed to parse gro title");
        goto done;
    }
    md_text_next_line(&c, e);
    if (c == e) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        goto done;
    }
    const int64_t num_atoms_in_file = parse_int(str_trim(md_text_next_line(&c, e)));
    if (num_atoms_in_file <= 0 || num_atoms_in_file >= UINT32_MAX) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        goto done;
    }
    const size_t num_atoms = (size_t)num_atoms_in_file;
    size_t taken = (size_t)(c - block.ptr);

    if (!md_system_state_init(state, num_atoms)) goto done;
    {
        // The atoms' own are written as they are read (gro_parse_lines), the padding past them here
        const size_t capacity = MD_ROUND_UP(num_atoms, 16);
        md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
        md_array_resize(sys->atom.flags, capacity, sys->alloc);
        MD_MEMSET(sys->atom.type_idx + num_atoms, 0, (capacity - num_atoms) * sizeof(md_particle_type_idx_t));
        MD_MEMSET(sys->atom.flags + num_atoms, 0, (capacity - num_atoms) * sizeof(md_atom_flags_t));
        sys->atom.count = num_atoms;
        md_particle_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
    }
    // Room for the residues of the load (no more than the atoms) and of a pass (no more than its lines, of 20 characters
    // at least), so that nothing grows (and copies) in the arena
    md_array_ensure(res, num_atoms, temp);
    md_array_ensure(pass_res, pass_size / 20 + 16, temp);

    // The columns, and whether there are velocities, as the first atom line shows them. GROMACS writes the velocities
    // of all atoms or of none: a file whose first atom has none is read without.
    L.num_atoms = num_atoms;
    L.xyz = state->xyz;
    L.res = res;
    L.atom_names = md_alloc(temp, num_atoms * sizeof(md_label_t));
    L.type_idx = sys->atom.type_idx;
    L.atom_flags = sys->atom.flags;
    {
        const char* fc = c;
        const str_t first = md_text_next_line(&fc, e);
        L.columns = gro_coordinate_columns(&L.width, &L.decimals, first);
        L.columns_w8 = L.columns && L.width == 8;
        float xyz[3], vel[3];
        gro_parse_coordinates(xyz, vel, first, L.columns, L.width, L.decimals);
        if (!gro_absent(vel[0])) {
            // Read into the attribute itself, committed once all are read
            L.velocity = md_attributes_reserve_atom_column(&sys->attributes, STR_LIT("atom/velocity"), gro_velocity_unit(), 3, num_atoms, &velocity_id);
        }
    }

    for (;;) {
        // The block a pass at a time: the lines counted, then read, of text the count left in the caches
        while (L.atom_beg < num_atoms && taken < block.len) {
            const char* text = block.ptr + taken;
            const size_t len = md_text_pass_len_of(text, block.len - taken, pass_size);
            size_t consumed = len;
            if (!gro_read_pass(&L, &pieces, &ranges, &pass_res, &num_res, &loose, &consumed, text, len, temp, ex)) goto done;  // Cancelled
            taken += consumed;
            // The residue fields of the last atom line, for the next pass
            if (L.atom_beg < num_atoms) {
                const str_t last = md_text_line_before(block.ptr, taken);
                const size_t n = MD_MIN(last.len, sizeof(prev_res));
                MD_MEMCPY(prev_res, last.ptr, n);
                L.prev = (str_t){ prev_res, n };
            }
        }
        md_reader_advance(src, taken);
        if (L.atom_beg >= num_atoms) break;
        block = md_reader_block_parallel(src, ex);
        taken = 0;
        if (block.len == 0) {
            MD_LOG_ERROR("The file ends after %zu of its %zu atoms", L.atom_beg, num_atoms);
            goto done;
        }
    }

    if (loose) {
        MD_LOG_DEBUG("Parsed GRO file using fixed-width format. This may indicate an issue with the file formatting.");
    }
    if (!md_reader_extract_line(src, &line)) {
        MD_LOG_ERROR("Failed to extract unitcell line");
        goto done;
    }
    if (!gro_parse_box(&state->unitcell, line)) goto done;

    // The residues
    md_array_resize(sys->residue.atom_offset, num_res + 1, sys->alloc);
    md_array_resize(sys->residue.name,   num_res, sys->alloc);
    md_array_resize(sys->residue.seq_id, num_res, sys->alloc);
    md_array_resize(sys->residue.flags,  num_res, sys->alloc);
    sys->residue.count = num_res;
    {
        gro_residues_t rb = { &sys->residue, res, num_res };
        if (!md_parallel_for(ex, (uint32_t)MD_DIV_UP(num_res, GRO_RES_CHUNK), 1, gro_build_residues, &rb)) goto done;    // Cancelled
    }
    sys->residue.atom_offset[num_res] = (uint32_t)num_atoms;    // Final sentinel

    if (L.velocity) {
        md_attributes_commit_atom_column(&sys->attributes, velocity_id);
    }
    md_util_system_infer_particle_types_from_labels(sys, L.atom_names, ex);
    state->num_atoms = num_atoms;
    result = true;

done:
    // The scratch goes with the scope
    md_temp_end(temp_scope);
    if (!result) {
        md_system_reset(sys);
        md_system_state_init(state, 0);
    }
    return result;
}

#undef GRO_ROOM

bool md_gro_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename, const md_executor_t* ex) {
    // The reader's buffer in a temp scope, below the one of the load (gro_load): a few passes of it
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_file_reader_t fr;
    const size_t cap = MD_MAX((size_t)MD_TEXT_BLOCK_SIZE, 2 * md_text_pass_size(md_executor_num_threads(ex)));
    if (!md_file_reader_open(&fr, filename, cap, md_temp_allocator(temp))) {
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
