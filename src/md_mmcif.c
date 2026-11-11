#include "md_mmcif.h"

#include <core/md_str.h>
#include <core/md_file.h>
#include <core/md_allocator.h>
#include <core/md_parse.h>
#include <core/md_parse_simd.h>
#include <core/md_reader.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>
#include <core/md_executor.h>
#include <core/md_vec_math.h>
#include <md_util.h>
#include <md_system.h>
#include <md_element.h>

#include <math.h>
#include <inttypes.h>


// Enumerate fields in _atom_site
enum {
    ATOM_SITE_GROUP_PDB,
    ATOM_SITE_ID,
    ATOM_SITE_TYPE_SYMBOL,
    ATOM_SITE_LABEL_ATOM_ID,
    ATOM_SITE_LABEL_ALT_ID,
    ATOM_SITE_LABEL_COMP_ID,
    ATOM_SITE_LABEL_ASYM_ID,
    ATOM_SITE_LABEL_ENTITY_ID,
    ATOM_SITE_LABEL_SEQ_ID,
    ATOM_SITE_PDBX_PDB_INS_CODE,
    ATOM_SITE_CARTN_X,
    ATOM_SITE_CARTN_Y,
    ATOM_SITE_CARTN_Z,
    ATOM_SITE_OCCUPANCY,
    ATOM_SITE_B_ISO_OR_EQUIV,
    ATOM_SITE_PDBX_FORMAL_CHARGE,
    ATOM_SITE_AUTH_SEQ_ID,
    ATOM_SITE_AUTH_COMP_ID,
    ATOM_SITE_AUTH_ASYM_ID,
    ATOM_SITE_AUTH_ATOM_ID,
    ATOM_SITE_PDBX_PDB_MODEL_NUM,
    ATOM_SITE_COUNT
};

static const md_str_t atom_site_labels[] = {
    MD_STR_INIT("group_PDB"),
    MD_STR_INIT("id"),
    MD_STR_INIT("type_symbol"),
    MD_STR_INIT("label_atom_id"),
    MD_STR_INIT("label_alt_id"),
    MD_STR_INIT("label_comp_id"),
    MD_STR_INIT("label_asym_id"),
    MD_STR_INIT("label_entity_id"),
    MD_STR_INIT("label_seq_id"),
    MD_STR_INIT("pdbx_PDB_ins_code"),
    MD_STR_INIT("Cartn_x"),
    MD_STR_INIT("Cartn_y"),
    MD_STR_INIT("Cartn_z"),
    MD_STR_INIT("occupancy"),
    MD_STR_INIT("B_iso_or_equiv"),
    MD_STR_INIT("pdbx_formal_charge"),
    MD_STR_INIT("auth_seq_id"),
    MD_STR_INIT("auth_comp_id"),
    MD_STR_INIT("auth_asym_id"),
    MD_STR_INIT("auth_atom_id"),
    MD_STR_INIT("pdbx_PDB_model_num"),
};

static const int required_atom_site_fields[] = {
    ATOM_SITE_ID,
    ATOM_SITE_TYPE_SYMBOL,
    ATOM_SITE_LABEL_ATOM_ID,
    ATOM_SITE_LABEL_ALT_ID,
    ATOM_SITE_LABEL_COMP_ID,
    ATOM_SITE_LABEL_ASYM_ID,
    ATOM_SITE_LABEL_ENTITY_ID,
    ATOM_SITE_LABEL_SEQ_ID,
    ATOM_SITE_CARTN_X,
    ATOM_SITE_CARTN_Y,
    ATOM_SITE_CARTN_Z,
};

enum {
    ATOM_SITE_GROUP_PDB_ATOM,
    ATOM_SITE_GROUP_PDB_HETATM,
};

typedef uint8_t atom_site_group_pdb_t;

typedef struct {

    char label_atom_id[8];
    char label_comp_id[8];
    char label_entity_id[8];
//    char pdbx_PDB_ins_code;

    float x, y, z;
    int label_seq_id;
    int auth_seq_id;

    // The optional per atom columns. Kept as floats even for the integral formal charge, because
    // that is what they are published as - one scalar per atom - and converting once here beats
    // three columns that each need a different read. NAN marks a column the file does not define,
    // so a partially populated one is still publishable and an absent one is not silently a zero.
    float occupancy;
    float b_iso_or_equiv;
    float formal_charge;

    // As long as a label holds (md_label_t): large entries have asym ids of four characters, which the assemblies name
    char label_asym_id[8];
    char auth_asym_id[8];
    char type_symbol[3];
    char label_alt_id;
} mmcif_atom_site_entry_t;

typedef struct {
    double a, b, c;           // Lengths of cell edges
    double alpha, beta, gamma; // Angles between cell edges
} mmcif_cell_params_t;

enum {
    ENTITY_TYPE_UNKNOWN,
    ENTITY_TYPE_BRANCHED,
    ENTITY_TYPE_MACROLIDE,
    ENTITY_TYPE_NON_POLYMER,
    ENTITY_TYPE_POLYMER,
    ENTITY_TYPE_WATER,
};

static const md_str_t entity_type_str[] = {
    MD_STR_INIT(""),
    MD_STR_INIT("branched"),
    MD_STR_INIT("macrolide"),
    MD_STR_INIT("non-polymer"),
    MD_STR_INIT("polymer"),
    MD_STR_INIT("water"),
};

typedef int mmcif_entity_type_t;

static inline mmcif_entity_type_t mmcif_entity_type_from_str(md_str_t type) {
    for (size_t i = 0; i < MD_COUNTOF(entity_type_str); ++i) {
        if (md_str_eq(type, entity_type_str[i])) {
            return (mmcif_entity_type_t)i;
        }
    }
    return ENTITY_TYPE_UNKNOWN;
}

enum {
    ENTITY_POLY_TYPE_UNKNOWN,
    ENTITY_POLY_TYPE_CYCLIC_PSEUDO_PEPTIDE, // Cyclic protein
    ENTITY_POLY_TYPE_OTHER,
    ENTITY_POLY_TYPE_PEPTIDE_NUCLEIC_ACID,  // Protein like
    ENTITY_POLY_TYPE_POLYDEOXYRIBONUCLEOTIDE,  // DNA
    ENTITY_POLY_TYPE_POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID,  // DNA/RNA hybrid
    ENTITY_POLY_TYPE_POLYPEPTIDE_D,                                      // D-chiral protein
    ENTITY_POLY_TYPE_POLYPEPTIDE_L,                                      // L-chiral protein
    ENTITY_POLY_TYPE_POLYRIBONUCLEOTIDE         // RNA
};

static const md_str_t entity_poly_type_str[] = {
    MD_STR_INIT(""),
    MD_STR_INIT("cyclic-pseudo-peptide"),
    MD_STR_INIT("other"),
    MD_STR_INIT("peptide nucleic acid"),
    MD_STR_INIT("polydeoxyribonucleotide"),
    MD_STR_INIT("polydeoxyribonucleotide/polyribonucleotide hybrid"),
    MD_STR_INIT("polypeptide(D)"),
    MD_STR_INIT("polypeptide(L)"),
    MD_STR_INIT("polyribonucleotide"),
};

typedef int mmcif_entity_poly_type_t;

static inline mmcif_entity_poly_type_t mmcif_entity_poly_type_from_str(md_str_t type) {
    for (size_t i = 0; i < MD_COUNTOF(entity_poly_type_str); ++i) {
        if (md_str_eq(type, entity_poly_type_str[i])) {
            return (mmcif_entity_poly_type_t)i;
        }
    }
    return ENTITY_POLY_TYPE_UNKNOWN;
}

typedef struct {
    md_str_t id;
    mmcif_entity_type_t type;
    mmcif_entity_poly_type_t poly_type;
    md_str_t description;
} mmcif_entity_t;

static md_part_kind_t mmcif_entity_kind(const mmcif_entity_t* entity) {
    switch (entity->type) {
    case ENTITY_TYPE_POLYMER:
        switch (entity->poly_type) {
        case ENTITY_POLY_TYPE_POLYPEPTIDE_L:
        case ENTITY_POLY_TYPE_POLYPEPTIDE_D:
        case ENTITY_POLY_TYPE_CYCLIC_PSEUDO_PEPTIDE:
            return MD_PART_KIND_PEPTIDE;
        case ENTITY_POLY_TYPE_POLYDEOXYRIBONUCLEOTIDE:
            return MD_PART_KIND_DNA;
        case ENTITY_POLY_TYPE_POLYRIBONUCLEOTIDE:
            return MD_PART_KIND_RNA;
        case ENTITY_POLY_TYPE_POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID:
            return MD_PART_KIND_NUCLEIC;
        case ENTITY_POLY_TYPE_PEPTIDE_NUCLEIC_ACID:
        case ENTITY_POLY_TYPE_OTHER:
            return MD_PART_KIND_POLYMER;
        default:
            // Classified from its residues once they are (md_part_infer)
            return MD_PART_KIND_UNKNOWN;
        }
    case ENTITY_TYPE_BRANCHED:
        return MD_PART_KIND_BRANCHED;
    case ENTITY_TYPE_WATER:
        return MD_PART_KIND_WATER;
    case ENTITY_TYPE_NON_POLYMER:
    case ENTITY_TYPE_MACROLIDE:
        return MD_PART_KIND_NON_POLYMER;
    default:
        return MD_PART_KIND_UNKNOWN;
    }
}

static inline mmcif_entity_t* mmcif_entity_find(md_array(mmcif_entity_t) entities, md_str_t id) {
    for (mmcif_entity_t* it = md_array_beg(entities); it != md_array_end(entities); ++it) {
        if (md_str_eq(it->id, id)) return it;
    }
    return NULL;
}

static inline mmcif_entity_t* mmcif_entity_find_or_create(md_array(mmcif_entity_t)* entities, md_str_t id, md_allocator_i* alloc) {
    mmcif_entity_t* it = mmcif_entity_find(*entities, id);
    if (it) return it;

    mmcif_entity_t ent = {
        .id = md_alloc_str(alloc, id),
    };
    md_array_push(*entities, ent, alloc);
    return md_array_last(*entities);
}

// An optional _atom_site column: absent from this file, or present and holding one of mmCIF's two
// null tokens ('.' unknown, '?' not applicable). NAN for all three, which the publisher below reads
// as "this atom has no value" without confusing it with a real zero.
static inline bool float_bits_nan(float v) {
    uint32_t u;
    md_memcpy(&u, &v, sizeof(u));
    return (u & 0x7F800000u) == 0x7F800000u && (u & 0x007FFFFFu) != 0;
}

static inline float mmcif_optional_float(const md_str_t tok[], const int table[], int field) {
    if (table[field] == -1) {
        return NAN;
    }
    md_str_t t = tok[table[field]];
    if (t.len == 0 || (t.len == 1 && (t.ptr[0] == '.' || t.ptr[0] == '?'))) {
        return NAN;
    }
    return (float)parse_float(t);
}

static inline int parse_int_with_default(md_str_t tok, int def) {
    if (tok.len == 0 || tok.ptr[0] == '.' || tok.ptr[0] == '?') {
        return def;
    }
    return (int)parse_int(tok);
}

typedef struct {
    md_reader_i* reader;
    md_str_writer_t sb;
    md_str_t line;
    md_str_t peek;
    bool has_peek;
    bool line_start;
    bool in_loop;
    bool plain;         // The token fetched last is as the line has it (not a quoted string or a text field)
    bool peek_plain;    // and the one peeked
    bool pending_start; // line is a line given back whole (mmcif_parse_atom_site): its first token is at a line start
} mmcif_parse_state_t;

// Returns true if token was extracted,
// False if EOF / stream
static bool mmcif_fetch_token(md_str_t* out_tok, mmcif_parse_state_t* s) {
    md_str_t tok = {0};
    s->line_start = s->pending_start;
    s->pending_start = false;
    s->plain = true;
    while (true) {
        if (md_str_empty(s->line)) {
            if (!md_reader_extract_line(s->reader, &s->line)) {
                return false;
            }
            s->line_start = true;
        }
        if (!md_parse_next_token(&tok, &s->line) || md_str_len(tok) == 0 || (tok.len > 0 && tok.ptr[0] == '#')) {
            s->line = (md_str_t){0};
            continue;
        }

        char first_char = tok.ptr[0];
        if (first_char == '\'' || first_char == '\"') {
            s->plain = false;
            char delim = tok.ptr[0];
            if (tok.len > 1 && tok.ptr[tok.len-1] == delim) {
                tok.ptr +=1;
                tok.len -=2;
            } else {
                size_t loc = 0;
                if ((loc = md_str_find_char(s->line, delim)) == MD_STR_NPOS) {
                    MD_LOG_ERROR("Unmatched string");
                    return false;
                }
                tok.ptr  += 1;
                tok.len  += loc;
                s->line.ptr += (loc+1);
                s->line.len -= (loc+1);
            }
        } else if (first_char == ';') {
            s->plain = false;
            md_str_writer_reset(&s->sb);
            md_str_writer_push_str(&s->sb, md_str_substr(tok, 1, SIZE_MAX));
            md_str_writer_push_str(&s->sb, s->line);

            bool term = false;
            while (md_reader_extract_line(s->reader, &s->line)) {
                if (s->line.len > 0 && s->line.ptr[0] == ';') {
                    term = true;
                    s->line = (md_str_t) {0};
                    break;
                }
                md_str_writer_push_char(&s->sb, '\n');
                md_str_writer_push_str(&s->sb, s->line);
            }
            if (!term) {
                MD_LOG_ERROR("Unterminated multiline string");
                return false;
            }
            tok = md_str_writer_to_str(&s->sb);
        }
        break;
    }

    *out_tok = tok;
    return true;
}

static inline bool mmcif_next_token(md_str_t* tok, mmcif_parse_state_t* state) {
    MD_ASSERT(tok);
    MD_ASSERT(state);
    if (state->has_peek) {
        *tok = state->peek;
        state->has_peek = false;
        return true;
    }
    return mmcif_fetch_token(tok, state);
}

static inline bool mmcif_peek_token(md_str_t* tok, mmcif_parse_state_t* state) {
    MD_ASSERT(tok);
    MD_ASSERT(state);
    if (state->has_peek) {
        *tok = state->peek;
        return true;
    }
    if (!mmcif_fetch_token(&state->peek, state)) {
        return false;
    }
    state->has_peek = true;
    state->peek_plain = state->plain;
    *tok = state->peek;
    return true;
}

// This is a bit wierd, since were mixing lines and token semantics.
// But will return the current line within the state that is about to be tokenized
static inline md_str_t mmcif_cur_line(mmcif_parse_state_t* s) {
    MD_ASSERT(s);     
    return s->line;
}

// Advances to the next line, which is not part of a multiline string token
static inline bool mmcif_advance_to_next_line(mmcif_parse_state_t* s) {
    MD_ASSERT(s);
    md_str_t tok = {0};
    while (mmcif_peek_token(&tok, s)) {
        if (s->line_start) return true;
        mmcif_next_token(&tok, s);
    }
    return false; // EOF
}

// Advances to a new line which starts with a valid control or item token
static inline bool mmcif_advance_to_next_control_token(mmcif_parse_state_t* s) {
    MD_ASSERT(s);

    md_str_t tok = {0};
    while (mmcif_peek_token(&tok, s)) {
        if (tok.len > 0 && s->line_start) {
            const char c0 = tok.ptr[0];

            // Category and control keywords:
            if (c0 == '_' || 
                md_str_eq(tok, MD_STR_LIT("loop_")) ||
                md_str_begins_with(tok, MD_STR_LIT("data_")) ||
                md_str_begins_with(tok, MD_STR_LIT("save_"))) 
            {
                return true;
            }
        }

        mmcif_next_token(&tok, s);
    }

    return false; // EOF
}



typedef struct {
    size_t num_fields;
    size_t num_rows;
    md_str_t* headers;
    md_str_t* values;
} mmcif_section_t;

// Help accessors for data in sections
static inline md_str_t mmcif_section_header(const mmcif_section_t* sec, size_t col) {
    if (!sec || col >= sec->num_fields) return (md_str_t){0};
    return sec->headers[col];
}

static inline md_str_t mmcif_section_value(const mmcif_section_t* sec, size_t row, size_t col) {
    if (!sec || col >= sec->num_fields || row >= sec->num_rows) return (md_str_t){0};
    return sec->values[row * sec->num_fields + col];
}

static inline bool mmcif_is_item(md_str_t t) { return t.len && t.ptr[0] == '_'; }
static inline bool mmcif_is_control(md_str_t t) {
    return (t.len >= 5 &&
        (md_str_eq(t,MD_STR_LIT("loop_")) ||
         md_str_eq(t,MD_STR_LIT("data_")) ||
         md_str_eq(t,MD_STR_LIT("save_")) ||
         md_str_eq(t,MD_STR_LIT("stop_"))));
}

static md_str_t mmcif_category(md_str_t item) {
    if (!mmcif_is_item(item)) return (md_str_t){0};
    md_str_t tail = md_str_substr(item,1,SIZE_MAX);
    size_t dot;
    if ((dot = md_str_find_char(tail, '.')) != MD_STR_NPOS) return (md_str_t){tail.ptr, dot};
    return (md_str_t){0};
}

// This is an internal procedure, we do not care about cleaning anything up
static bool mmcif_parse_section(mmcif_section_t* sec, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    MD_ASSERT(sec);
    MD_ASSERT(state);
    MD_ASSERT(alloc);
    md_memset(sec, 0, sizeof(*sec));

    if (state->in_loop) {
        // --------------- Collect headers ---------------
        md_str_t t;
        while (mmcif_peek_token(&t, state) && mmcif_is_item(t)) {
            if (!mmcif_next_token(&t, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure while reading loop headers");
                return false;
            }
            md_array_push(sec->headers, md_alloc_str(alloc, t), alloc);
        }
        size_t ncols = md_array_size(sec->headers);
        if (ncols == 0) {
            MD_LOG_ERROR("Loop section without headers");
            return false;
        }

        // --------------- Collect rows ---------------
        size_t values_read = 0;
        while (true) {
            md_str_t look;
            if (!mmcif_peek_token(&look, state)) {
                // EOF -> ok end
                break;
            }
            if (mmcif_is_control(look) || mmcif_is_item(look)) {
                // Reached structural boundary
                break;
            }
            // Read one row
            size_t got = 0;
            for (; got < ncols; ++got) {
                md_str_t v;
                if (!mmcif_next_token(&v, state)) {
                    MD_LOG_ERROR("Unexpected EOF in middle of loop row");
                    return false;
                }
                if (mmcif_is_control(v) || mmcif_is_item(v)) {
                    // Structural token unexpectedly inside a row -> error
                    MD_LOG_ERROR("Premature loop termination inside row (got %zu/%zu values)", got, ncols);
                    return false;
                }
                md_array_push(sec->values, md_alloc_str(alloc, v), alloc);
            }
            if (got != ncols) {
                // Should not happen due to logic above, but guard anyway
                MD_LOG_ERROR("Incomplete row (got %zu of %zu columns)", got, ncols);
                return false;
            }
            values_read += got;
        }

        if (values_read % ncols != 0) {
            MD_LOG_ERROR("Internal inconsistency: values_read %% ncols != 0 (%zu %% %zu)", values_read, ncols);
            return false;
        }

        sec->num_fields = ncols;
        sec->num_rows   = (ncols) ? values_read / ncols : 0;

        if (sec->num_rows == 0) {
            MD_LOG_ERROR("Loop section has headers but no data rows");
            return false;
        }
    } else {
        // --------------- Unlooped key/value section ---------------
        md_str_t first;
        if (!mmcif_peek_token(&first, state) || !mmcif_is_item(first)) {
            // Nothing to parse (not an error? choose: treat as failure to signal absence)
            return false;
        }
        md_str_t cat = mmcif_category(first);
        if (md_str_empty(cat)) {
            MD_LOG_ERROR("Malformed item (missing category segment): " MD_STR_FMT, MD_STR_ARG(first));
            return false;
        }

        while (true) {
            md_str_t key;
            if (!mmcif_peek_token(&key, state)) break;
            if (mmcif_is_control(key) || !mmcif_is_item(key)) break;

            // Category change ends this section
            md_str_t this_cat = mmcif_category(key);
            if (!md_str_empty(this_cat) && !md_str_eq(this_cat, cat)) break;

            // Consume key
            if (!mmcif_next_token(&key, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure consuming item key");
                return false;
            }
            md_array_push(sec->headers, md_alloc_str(alloc, key), alloc);

            // Value (optional)
            md_str_t val;
            if (!mmcif_peek_token(&val, state) || mmcif_is_control(val) || mmcif_is_item(val)) {
                // Missing value -> store empty
                return false;
            }
            if (!mmcif_next_token(&val, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure consuming item value");
                return false;
            }
            md_array_push(sec->values, md_alloc_str(alloc, val), alloc);
        }

        sec->num_fields = md_array_size(sec->headers);
        sec->num_rows   = (sec->num_fields > 0) ? 1 : 0;

        if (sec->num_fields == 0) {
            MD_LOG_ERROR("Unlooped section contained no items");
            return false;
        }
        if (md_array_size(sec->values) != sec->num_fields) {
            MD_LOG_ERROR("Key/value count mismatch (%zu headers vs %zu values)",
                (size_t)sec->num_fields, (size_t)md_array_size(sec->values));
            return false;
        }
    }

    return true;
}

static inline int mmcif_section_find_col_idx(const mmcif_section_t* sec, md_str_t field_label) {
    for (size_t i = 0; i < sec->num_fields; ++i) {
        if (md_str_eq(field_label, sec->headers[i])) {
            return (int)i;
        }
    }
    return -1;
}

static bool mmcif_parse_entity(md_array(mmcif_entity_t)* entities, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    MD_ASSERT(entities);
    MD_ASSERT(state);
    MD_ASSERT(alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    bool success = false;
    
    mmcif_section_t sec = {0};
    if (!mmcif_parse_section(&sec, state, temp_alloc)) {
        MD_LOG_ERROR("Failed to parse _entity section");
        goto done;
    }

    int id_col   = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_entity.id"));
    int type_col = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_entity.type"));
    int desc_col = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_entity.pdbx_description"));

    if (id_col == -1) {
        MD_LOG_ERROR("Missing entity field 'id'");
        goto done;
    }
    if (type_col == -1) {
        MD_LOG_ERROR("Missing entity field 'type'");
        goto done;
    }
    if (desc_col == -1) {
        MD_LOG_ERROR("Missing entity field 'pdbx_description'");
        goto done;
    }

    for (size_t i = 0; i < sec.num_rows; ++i) {
        md_str_t id = mmcif_section_value(&sec, i, id_col);
        if (!md_str_empty(id)) {
            mmcif_entity_t* ent = mmcif_entity_find_or_create(entities, id, alloc);
            MD_ASSERT(ent);
            ent->type = mmcif_entity_type_from_str(mmcif_section_value(&sec, i, type_col));
            ent->description = md_alloc_str(alloc, mmcif_section_value(&sec, i, desc_col));
        }
    }

    success = true;
done:
    md_temp_end(temp);
    return success;
}

// We only populate the entity_poly type within the given entities
static bool mmcif_parse_entity_poly(md_array(mmcif_entity_t)* entities, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    MD_ASSERT(entities);
    MD_ASSERT(state);
    MD_ASSERT(alloc);

    bool success = false;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    mmcif_section_t sec = {0};
    if (!mmcif_parse_section(&sec, state, temp_alloc)) {
        MD_LOG_ERROR("Failed to parse _entity_poly section");
        goto done;
    }

    int id_col   = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_entity_poly.entity_id"));
    int type_col = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_entity_poly.type"));

    if (id_col == -1) {
        MD_LOG_ERROR("Missing entity field 'entity_id'");
        goto done;
    }
    if (type_col == -1) {
        MD_LOG_ERROR("Missing entity field 'type'");
        goto done;
    }

    for (size_t i = 0; i < sec.num_rows; ++i) {
        md_str_t id = mmcif_section_value(&sec, i, id_col);
        if (!md_str_empty(id)) {
            mmcif_entity_t* ent = mmcif_entity_find_or_create(entities, id, alloc);
            MD_ASSERT(ent);
            ent->poly_type = mmcif_entity_poly_type_from_str(mmcif_section_value(&sec, i, type_col));
        }
    }

    success = true;
done:
    md_temp_end(temp);
    return success;
}

// ASSEMBLIES. A row of _pdbx_struct_assembly_gen: the asyms (label_asym_id) its operators copy
typedef struct {
    md_str_t assembly_id;
    md_str_t oper_expression;
    md_str_t asym_id_list;
} mmcif_assembly_gen_t;

// An operator of _pdbx_struct_oper_list
typedef struct {
    md_str_t  id;
    md_f32x4x4_t transform;
} mmcif_oper_t;

static bool mmcif_parse_assembly_gen(md_array(mmcif_assembly_gen_t)* rows, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool success = false;

    mmcif_section_t sec = {0};
    if (!mmcif_parse_section(&sec, state, temp_alloc)) {
        goto done;
    }
    const int id_col   = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_pdbx_struct_assembly_gen.assembly_id"));
    const int expr_col = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_pdbx_struct_assembly_gen.oper_expression"));
    const int asym_col = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_pdbx_struct_assembly_gen.asym_id_list"));
    if (id_col == -1 || expr_col == -1 || asym_col == -1) {
        MD_LOG_ERROR("_pdbx_struct_assembly_gen lacks assembly_id, oper_expression or asym_id_list");
        goto done;
    }
    for (size_t i = 0; i < sec.num_rows; ++i) {
        const mmcif_assembly_gen_t row = {
            .assembly_id     = md_alloc_str(alloc, mmcif_section_value(&sec, i, id_col)),
            .oper_expression = md_alloc_str(alloc, mmcif_section_value(&sec, i, expr_col)),
            .asym_id_list    = md_alloc_str(alloc, mmcif_section_value(&sec, i, asym_col)),
        };
        md_array_push(*rows, row, alloc);
    }
    success = true;
done:
    md_temp_end(temp);
    return success;
}

static bool mmcif_parse_oper_list(md_array(mmcif_oper_t)* opers, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool success = false;

    mmcif_section_t sec = {0};
    if (!mmcif_parse_section(&sec, state, temp_alloc)) {
        goto done;
    }
    const int id_col = mmcif_section_find_col_idx(&sec, MD_STR_LIT("_pdbx_struct_oper_list.id"));
    int m_col[3][3];
    int v_col[3];
    for (int r = 0; r < 3; ++r) {
        char name[64];
        for (int c = 0; c < 3; ++c) {
            const int len = snprintf(name, sizeof(name), "_pdbx_struct_oper_list.matrix[%d][%d]", r + 1, c + 1);
            m_col[r][c] = mmcif_section_find_col_idx(&sec, (md_str_t){name, (size_t)len});
        }
        const int len = snprintf(name, sizeof(name), "_pdbx_struct_oper_list.vector[%d]", r + 1);
        v_col[r] = mmcif_section_find_col_idx(&sec, (md_str_t){name, (size_t)len});
    }
    if (id_col == -1) {
        MD_LOG_ERROR("_pdbx_struct_oper_list lacks id");
        goto done;
    }
    for (int r = 0; r < 3; ++r) {
        if (v_col[r] == -1 || m_col[r][0] == -1 || m_col[r][1] == -1 || m_col[r][2] == -1) {
            MD_LOG_ERROR("_pdbx_struct_oper_list lacks the elements of its matrices or vectors");
            goto done;
        }
    }
    for (size_t i = 0; i < sec.num_rows; ++i) {
        // Row r column c of the matrix, column major as md_f32x4x4_t: elem[c][r]
        mmcif_oper_t op = {.id = md_alloc_str(alloc, mmcif_section_value(&sec, i, id_col)), .transform = md_f32x4x4_ident()};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                op.transform.elem[c][r] = (float)parse_float(mmcif_section_value(&sec, i, m_col[r][c]));
            }
            op.transform.elem[3][r] = (float)parse_float(mmcif_section_value(&sec, i, v_col[r]));
        }
        md_array_push(*opers, op, alloc);
    }
    success = true;
done:
    md_temp_end(temp);
    return success;
}

// The operators of an oper_expression, indices of opers: "1", "1,2,5", "1-4", "(1-4)" or a product of such groups
// "(X0)(1-60)", "(1-60)(61-88)": an operator for each combination of one of each group, the transform of the product
// in the order written (that of the last group applied first, as Mol* reads it), appended to opers. ids are numbers
// in a range ("1-4"), any name otherwise. False for an expression which names an operator not in opers.
static bool mmcif_oper_expression(md_array(uint32_t)* out, md_str_t expr, md_array(mmcif_oper_t)* opers, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* ta = md_temp_allocator(temp);
    bool result = false;

    // The groups: the parenthesized parts, or the whole expression
    md_array(md_str_t) groups = 0;
    expr = md_str_trim(expr);
    if (md_str_contains_char(expr, '(')) {
        while (!md_str_empty(expr)) {
            size_t beg, end;
            if (((beg = md_str_find_char(expr, '(')) == MD_STR_NPOS) || ((end = md_str_find_char(expr, ')')) == MD_STR_NPOS) || end < beg) {
                MD_LOG_ERROR("Unbalanced operator expression: '" MD_STR_FMT "'", MD_STR_ARG(expr));
                goto done;
            }
            md_array_push(groups, md_str_substr(expr, beg + 1, end - beg - 1), ta);
            expr = md_str_trim(md_str_substr(expr, end + 1, SIZE_MAX));
        }
    } else {
        md_array_push(groups, expr, ta);
    }

    // The operators of each group
    const size_t group_count = md_array_size(groups);
    md_array(uint32_t)* group_ops = md_temp_alloc_zero(temp, MD_MAX(group_count, 1) * sizeof(md_array(uint32_t)));
    const size_t base_count = md_array_size(*opers);
    for (size_t g = 0; g < group_count; ++g) {
        md_str_t list = groups[g];
        md_str_t item;
        while (md_parse_next_token_delim(&item, &list, ',')) {
            item = md_str_trim(item);
            if (md_str_empty(item)) continue;
            md_str_t lo_str, hi_str;
            int64_t lo = 0, hi = -1;
            if (md_str_cut(&lo_str, &hi_str, item, md_str_find_char(item, '-')) && lo_str.len > 0) {
                lo = parse_int(md_str_trim(lo_str));
                hi = parse_int(md_str_trim(hi_str));
            }
            if (hi >= lo) {
                for (int64_t v = lo; v <= hi; ++v) {
                    size_t k = 0;
                    for (; k < base_count; ++k) {
                        const md_str_t id = (*opers)[k].id;
                        if (!md_str_empty(id) && md_ascii_is_digit(id.ptr[0]) && parse_int(id) == v) break;
                    }
                    if (k == base_count) {
                        MD_LOG_ERROR("Operator %" PRId64 " of the expression is not in _pdbx_struct_oper_list", v);
                        goto done;
                    }
                    md_array_push(group_ops[g], (uint32_t)k, ta);
                }
            } else {
                size_t k = 0;
                for (; k < base_count; ++k) {
                    if (md_str_eq((*opers)[k].id, item)) break;
                }
                if (k == base_count) {
                    MD_LOG_ERROR("Operator '" MD_STR_FMT "' of the expression is not in _pdbx_struct_oper_list", MD_STR_ARG(item));
                    goto done;
                }
                md_array_push(group_ops[g], (uint32_t)k, ta);
            }
        }
        if (md_array_size(group_ops[g]) == 0) {
            MD_LOG_ERROR("An empty group in an operator expression");
            goto done;
        }
    }

    if (group_count == 1) {
        md_array_push_array(*out, group_ops[0], md_array_size(group_ops[0]), alloc);
    } else {
        // Every combination, the first group varying slowest
        size_t total = 1;
        for (size_t g = 0; g < group_count; ++g) {
            total *= md_array_size(group_ops[g]);
            if (total > 1000000) {
                MD_LOG_ERROR("An operator expression of more than a million operators");
                goto done;
            }
        }
        size_t* idx = md_temp_alloc_zero(temp, group_count * sizeof(size_t));
        for (size_t n = 0; n < total; ++n) {
            mmcif_oper_t op = {.transform = md_f32x4x4_ident()};
            char name[64];
            int len = 0;
            for (size_t g = 0; g < group_count; ++g) {
                const mmcif_oper_t* f = *opers + group_ops[g][idx[g]];
                op.transform = md_f32x4x4_mul(op.transform, f->transform);
                len += snprintf(name + len, sizeof(name) - len, "%s" MD_STR_FMT, g ? "x" : "", MD_STR_ARG(f->id));
                len = MD_MIN(len, (int)sizeof(name) - 1);
            }
            op.id = md_alloc_str(alloc, (md_str_t){name, (size_t)len});
            md_array_push(*out, (uint32_t)md_array_size(*opers), alloc);
            md_array_push(*opers, op, alloc);
            for (size_t g = group_count; g-- > 0; ) {
                if (++idx[g] < md_array_size(group_ops[g])) break;
                idx[g] = 0;
            }
        }
    }
    result = true;
done:
    md_temp_end(temp);
    return result;
}

// The assemblies of the system from the rows of _pdbx_struct_assembly_gen and the operators of _pdbx_struct_oper_list.
// An asym names every instance of its id: a polymer chain, or each of the molecules of a non-polymer asym (its waters).
static void mmcif_assemblies(md_system_t* sys, const md_array(mmcif_assembly_gen_t) rows, md_array(mmcif_oper_t) opers, md_allocator_i* temp_alloc) {
    const size_t row_count = md_array_size(rows);
    if (row_count == 0 || md_array_size(opers) == 0) return;

    md_array(md_assembly_gen_t) gen = 0;
    // The operators of rows with the same expression, which the products of one keep the same
    md_array(md_str_t) expr_seen = 0;
    md_array(md_array(uint32_t)) expr_ops = 0;

    for (size_t r = 0; r < row_count; ++r) {
        md_array(md_instance_idx_t) inst = 0;
        md_str_t list = rows[r].asym_id_list;
        md_str_t asym;
        while (md_parse_next_token_delim(&asym, &list, ',')) {
            asym = md_str_trim(asym);
            if (md_str_empty(asym)) continue;
            bool found = false;
            for (size_t i = 0; i < sys->instance.count; ++i) {
                if (md_str_eq(md_system_instance_id(sys, i), asym)) {
                    md_array_push(inst, (md_instance_idx_t)i, temp_alloc);
                    found = true;
                }
            }
            if (!found) {
                // An asym of no atoms the file kept (one of alternate locations only), or of those left out
                MD_LOG_DEBUG("Assembly '" MD_STR_FMT "' names asym '" MD_STR_FMT "', which has no atoms", MD_STR_ARG(rows[r].assembly_id), MD_STR_ARG(asym));
            }
        }

        const md_str_t expr = md_str_trim(rows[r].oper_expression);
        md_array(uint32_t) ops = 0;
        for (size_t e = 0; e < md_array_size(expr_seen); ++e) {
            if (md_str_eq(expr_seen[e], expr)) {
                ops = expr_ops[e];
                break;
            }
        }
        if (!ops) {
            if (!mmcif_oper_expression(&ops, expr, &opers, temp_alloc)) {
                MD_LOG_ERROR("Assembly '" MD_STR_FMT "': could not read the operator expression '" MD_STR_FMT "'", MD_STR_ARG(rows[r].assembly_id), MD_STR_ARG(expr));
                continue;
            }
            md_array_push(expr_seen, expr, temp_alloc);
            md_array_push(expr_ops, ops, temp_alloc);
        }

        const md_assembly_gen_t g = {
            .assembly_id = md_label_make(rows[r].assembly_id),
            .inst        = inst,
            .inst_count  = md_array_size(inst),
            .oper        = ops,
            .oper_count  = md_array_size(ops),
        };
        md_array_push(gen, g, temp_alloc);
    }

    const size_t oper_count = md_array_size(opers);
    md_f32x4x4_t*     transform = md_alloc(temp_alloc, oper_count * sizeof(md_f32x4x4_t));
    md_label_t* id        = md_alloc(temp_alloc, oper_count * sizeof(md_label_t));
    for (size_t i = 0; i < oper_count; ++i) {
        transform[i] = opers[i].transform;
        id[i] = md_label_make(opers[i].id);
    }
    md_system_assemblies_init(sys, gen, md_array_size(gen), transform, id, oper_count);
}

// The tokens of a row which is a line of its own (from beg, its first token, to end), as mmcif_fetch_token gives them:
// whitespace delimited, a quoted one without its quotes. False where the line is not num tokens, or one of them is not
// read that way (a comment, a text field, a quoted string with whitespace in it): the tokenizer reads the row instead.
//
// The whitespace of the line a chunk of 64 bytes at a time (md_parse_whitespace_mask64; no chunk from end on), and of each
// chunk where the tokens begin (a byte which is not whitespace after one which is) and where they end (whitespace after
// a byte which is not), all of them at once: the n:th token is from the n:th beginning to the n:th end. Not one token
// after another, each from where the last ended, which is a chain of as many steps as the row has tokens.
// The whitespace of the chunk at c, the bytes from end on counted as such
MD_TEXT_NO_SANITIZE static inline uint64_t mmcif_ws_chunk(const char* c, const char* end) {
    uint64_t ws = md_parse_whitespace_mask64(c);
    if (end - c < 64) ws |= ~UINT64_C(0) << (end - c);
    return ws;
}

// What mmcif_fetch_token reads otherwise when a token begins with it: 1 a comment or a text field, 2 a quoted string
static const uint8_t mmcif_first_char[256] = { ['#'] = 1, [';'] = 1, ['\''] = 2, ['"'] = 2 };

// tok has room for MMCIF_SPLIT_SLACK tokens past num, which the loops below write ahead of the count
#define MMCIF_SPLIT_SLACK 4

MD_TEXT_NO_SANITIZE static bool mmcif_split_row(md_str_t* tok, size_t num, const char* beg, const char* end) {
    size_t n = 0;
    uint32_t first = 0;         // Of the tokens' first characters, what mmcif_first_char says of them
    size_t n_end = 0;
    const char* chunk = (const char*)((uintptr_t)beg & ~(uintptr_t)63);
    uint64_t before = (UINT64_C(1) << (beg - chunk)) - 1;      // The bytes before beg count as whitespace
    uint64_t carry = 1;                                         // Whether the byte before the chunk is whitespace
    for (; chunk < end; chunk += 64) {
        const uint64_t ws = mmcif_ws_chunk(chunk, end) | before;
        before = 0;
        const uint64_t prev = (ws << 1) | carry;                // Of each byte, whether the one before it is whitespace
        carry = ws >> 63;
        uint64_t starts = ~ws & prev;
        uint64_t ends = ws & ~prev;
        const size_t count = (size_t)md_popcnt64(starts);
        if (n + count > num) return false;                      // More tokens than columns
        // Four at a time, a branch for four (a bit past the last is no bit: md_ctz64(0) is past the chunk, written to the
        // room past the count, and a later chunk's token there writes over it)
        for (size_t k = 0; k < count; k += 4) {
            for (int j = 0; j < 4; ++j) {
                const char* t = chunk + (starts ? md_ctz64(starts) : 0);
                tok[n + k + j].ptr = t;
                first |= starts ? mmcif_first_char[(uint8_t)*t] : 0;
                starts &= starts - 1;
            }
        }
        n += count;
        // Of tokens which began here or before
        const size_t count_end = (size_t)md_popcnt64(ends);
        for (size_t k = 0; k < count_end; k += 4) {
            for (int j = 0; j < 4; ++j) {
                const size_t i = n_end + k + j;
                if (i < n) tok[i].len = (size_t)(chunk + (ends ? md_ctz64(ends) : 0) - tok[i].ptr);
                ends &= ends - 1;
            }
        }
        n_end += count_end;
    }
    // The last one ends with the line, where the line ends with a chunk
    if (n_end < n) {
        tok[n_end].len = (size_t)(end - tok[n_end].ptr);
    }
    if (n != num) return false;
    if (first & 1) return false;
    // As mmcif_fetch_token: a quoted token without its quotes, where it ends in its quote
    for (size_t i = 0; first && i < num; ++i) {
        const char f = tok[i].ptr[0];
        if (f == '\'' || f == '"') {
            if (tok[i].len < 2 || tok[i].ptr[tok[i].len - 1] != f) return false;
            tok[i].ptr += 1;
            tok[i].len -= 2;
        }
    }
    return true;
}

// The columns of _atom_site: where each field is among the tokens of a row (-1 for none), and how many there are
typedef struct mmcif_atom_site_cols_t {
    int    table[ATOM_SITE_COUNT];
    size_t num_cols;
    bool   have_auth_seq_id;
    bool   have_auth_asym_id;
} mmcif_atom_site_cols_t;

// The entry of a row of _atom_site, from its tokens
static inline void mmcif_atom_site_entry(mmcif_atom_site_entry_t* out, const md_str_t tok[], const mmcif_atom_site_cols_t* c) {
    const int* table = c->table;
    //md_str_t group_PDB       = tok[table[ATOM_SITE_GROUP_PDB]];
    md_str_t type_symbol     = tok[table[ATOM_SITE_TYPE_SYMBOL]];
    md_str_t label_atom_id   = tok[table[ATOM_SITE_LABEL_ATOM_ID]];
    md_str_t label_alt_id    = tok[table[ATOM_SITE_LABEL_ALT_ID]];
    md_str_t label_asym_id   = tok[table[ATOM_SITE_LABEL_ASYM_ID]];
    md_str_t label_comp_id   = tok[table[ATOM_SITE_LABEL_COMP_ID]];
    md_str_t label_entity_id = tok[table[ATOM_SITE_LABEL_ENTITY_ID]];
    md_str_t label_seq_id    = tok[table[ATOM_SITE_LABEL_SEQ_ID]];
    md_str_t auth_seq_id     = c->have_auth_seq_id ? tok[table[ATOM_SITE_AUTH_SEQ_ID]] : (md_str_t){0};
    md_str_t auth_asym_id    = c->have_auth_asym_id ? tok[table[ATOM_SITE_AUTH_ASYM_ID]] : (md_str_t){0};
    md_str_t cartn_x         = tok[table[ATOM_SITE_CARTN_X]];
    md_str_t cartn_y         = tok[table[ATOM_SITE_CARTN_Y]];
    md_str_t cartn_z         = tok[table[ATOM_SITE_CARTN_Z]];

    const int default_value = INT32_MIN;

    mmcif_atom_site_entry_t entry = {
        .x = (float)parse_float(cartn_x),
        .y = (float)parse_float(cartn_y),
        .z = (float)parse_float(cartn_z),
        .label_seq_id = parse_int_with_default(label_seq_id, default_value),
        .auth_seq_id = parse_int_with_default(auth_seq_id, default_value),
        .label_alt_id   = label_alt_id.len > 0 ? label_alt_id.ptr[0] : '.',
        .occupancy      = mmcif_optional_float(tok, table, ATOM_SITE_OCCUPANCY),
        .b_iso_or_equiv = mmcif_optional_float(tok, table, ATOM_SITE_B_ISO_OR_EQUIV),
        .formal_charge  = mmcif_optional_float(tok, table, ATOM_SITE_PDBX_FORMAL_CHARGE),
    };

    md_str_copy_to_char_buf(entry.label_atom_id,   sizeof(entry.label_atom_id),    label_atom_id);
    md_str_copy_to_char_buf(entry.label_asym_id,   sizeof(entry.label_asym_id),    label_asym_id);
    md_str_copy_to_char_buf(entry.label_comp_id,   sizeof(entry.label_comp_id),    label_comp_id);
    md_str_copy_to_char_buf(entry.label_entity_id, sizeof(entry.label_entity_id),  label_entity_id);
    md_str_copy_to_char_buf(entry.type_symbol,     sizeof(entry.type_symbol),      type_symbol);
    md_str_copy_to_char_buf(entry.auth_asym_id,    sizeof(entry.auth_asym_id),     auth_asym_id);

    *out = entry;
}

// A line which is a row of its own, as the rows loop reads one: its tokens (mmcif_split_row), where it is neither
// quoted at its first token nor an item or a control word (what comes after the rows)
static inline bool mmcif_plain_row(md_str_t* tok, size_t num_cols, md_str_t line) {
    const char* b = line.ptr;
    const char* e = line.ptr + line.len;
    while (b < e && md_ascii_is_whitespace(*b)) ++b;
    return b < e && *b != '\'' && *b != '"' && mmcif_split_row(tok, num_cols, b, e) && !mmcif_is_item(tok[0]) && !mmcif_is_control(tok[0]);
}

// ### THE ROWS OVER THE THREADS ###
//
// The rows of _atom_site are most of a large file, and independent of one another: they are read a block of whole lines
// at a time from the reader (md_reader_block_parallel), each block in passes, each pass in pieces at line boundaries (a
// few for each thread, md_text_split). Two loops over the pieces of a pass (md_parallel_for), each a fork and a join:
//
//   count  the lines of each piece; their sum before a piece is where its rows go among the entries
//   parse  the rows of each piece into their places, up to its first line which is not a row of its own (a comment, a
//          row over lines, a quoted first token, what comes after the rows): where the piece stops
//
// Then, in the order of the pieces, the rows up to the first stop are the entries the pass adds; the reader is taken to
// that line, which the tokenizer reads (as it reads any), and the blocks go on after it. The entries are placed by their
// lines, so they are the same whichever thread read which piece. A pass of one piece is not counted: its rows go after
// those before it.

#define MMCIF_BLOCK_MAX_COLS 128     // A row of more columns is read by the tokenizer

typedef struct mmcif_rows_piece_t {
    md_text_piece_t range;      // of the pass
    size_t num_lines;           // the count
    size_t first;               // where its rows go, among the entries of the pass
    size_t num_rows;            // its lines which are rows, up to the first which is not
    size_t stop;                // that line (of the pass), range.end where all are rows
} mmcif_rows_piece_t;

typedef struct mmcif_rows_t {
    const char* text;           // the pass
    size_t len;
    mmcif_rows_piece_t* pieces;
    mmcif_atom_site_entry_t* entries;   // of the pass
    const mmcif_atom_site_cols_t* cols;
} mmcif_rows_t;

static void mmcif_rows_count(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const mmcif_rows_t* R = (const mmcif_rows_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        mmcif_rows_piece_t* piece = &R->pieces[p];
        size_t n = md_text_count_newlines(R->text + piece->range.beg, R->text + piece->range.end);
        // The last line of the text may not end in one
        if (piece->range.end == R->len && piece->range.end > piece->range.beg && R->text[R->len - 1] != '\n') ++n;
        piece->num_lines = n;
    }
}

static void mmcif_rows_parse(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const mmcif_rows_t* R = (const mmcif_rows_t*)arg;
    md_str_t tok[MMCIF_BLOCK_MAX_COLS + MMCIF_SPLIT_SLACK];
    for (uint32_t p = beg; p < end; ++p) {
        mmcif_rows_piece_t* piece = &R->pieces[p];
        mmcif_atom_site_entry_t* out = R->entries + piece->first;
        size_t n = 0;
        piece->stop = piece->range.end;
        md_text_lines_t lines = md_text_lines(R->text + piece->range.beg, R->text + piece->range.end);
        md_str_t line;
        while (md_text_lines_next(&lines, &line)) {
            if (!mmcif_plain_row(tok, R->cols->num_cols, line)) {
                piece->stop = (size_t)(line.ptr - R->text);
                break;
            }
            mmcif_atom_site_entry(&out[n++], tok, R->cols);
        }
        piece->num_rows = n;
    }
}

typedef enum {
    MMCIF_ROWS_END,         // the source ended
    MMCIF_ROWS_STOPPED,     // at a line which is not a row of its own: the reader is at it
    MMCIF_ROWS_CANCELLED,
} mmcif_rows_result_t;

// The rows which follow in the reader, as many as are rows of their own, added to the entries
static mmcif_rows_result_t mmcif_rows_from_blocks(md_array(mmcif_atom_site_entry_t)* entries, const mmcif_atom_site_cols_t* cols, md_reader_i* reader,
                                                  md_allocator_i* alloc, const md_executor_t* ex) {
    const uint32_t num_threads = md_executor_num_threads(ex);
    const size_t pass_size = md_text_pass_size(num_threads);
    const size_t max_pieces = md_text_num_pieces(pass_size, num_threads);
    md_text_piece_t*    ranges = md_alloc(alloc, max_pieces * sizeof(md_text_piece_t));
    mmcif_rows_piece_t* pieces = md_alloc(alloc, max_pieces * sizeof(mmcif_rows_piece_t));
    // Rows of fewer bytes than a character and a blank for each column there are not
    const size_t min_row = 2 * MD_MAX(cols->num_cols, 1);

    for (;;) {
        const md_str_t block = md_reader_block_parallel(reader, ex);
        if (block.len == 0) return MMCIF_ROWS_END;
        size_t taken = 0;
        while (taken < block.len) {
            const char* text = block.ptr + taken;
            const size_t len = md_text_pass_len_of(text, block.len - taken, pass_size);
            const size_t num_pieces = MD_MIN(md_text_num_pieces(len, num_threads), max_pieces);
            md_text_split(ranges, num_pieces, text, len);
            mmcif_rows_t R = { .text = text, .len = len, .pieces = pieces, .cols = cols };
            for (size_t p = 0; p < num_pieces; ++p) pieces[p] = (mmcif_rows_piece_t){ .range = ranges[p] };

            const size_t base = md_array_size(*entries);
            if (num_pieces == 1) {
                md_array_ensure(*entries, base + len / min_row + 1, alloc);
            } else {
                if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, mmcif_rows_count, &R)) return MMCIF_ROWS_CANCELLED;
                size_t num_lines = 0;
                for (size_t p = 0; p < num_pieces; ++p) {
                    pieces[p].first = num_lines;
                    num_lines += pieces[p].num_lines;
                }
                md_array_ensure(*entries, base + num_lines, alloc);
            }
            R.entries = *entries + base;
            if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, mmcif_rows_parse, &R)) return MMCIF_ROWS_CANCELLED;

            // The rows up to the first line which is not one
            size_t num_rows = 0;
            size_t stop = len;
            for (size_t p = 0; p < num_pieces; ++p) {
                num_rows += pieces[p].num_rows;
                if (pieces[p].stop < pieces[p].range.end) {
                    stop = pieces[p].stop;
                    break;
                }
            }
            md_array_resize(*entries, base + num_rows, alloc);
            taken += stop;
            if (stop < len) {
                md_reader_advance(reader, taken);
                return MMCIF_ROWS_STOPPED;
            }
        }
        md_reader_advance(reader, taken);
    }
}

static bool mmcif_parse_atom_site(md_array(mmcif_atom_site_entry_t)* atom_entries, mmcif_parse_state_t* state, md_allocator_i* alloc, const md_executor_t* ex) {
    MD_ASSERT(atom_entries);
    MD_ASSERT(state);
    MD_ASSERT(alloc);
    MD_ASSERT(state->in_loop == true);

    mmcif_atom_site_cols_t cols = {0};
    int* table = cols.table;
    md_memset(table, -1, sizeof(cols.table));
    size_t num_cols = 0;

    // Read category headers
    md_str_t item;
    while (mmcif_peek_token(&item, state) && mmcif_is_item(item) && md_str_begins_with(item, MD_STR_LIT("_atom_site."))) {
        if (!mmcif_next_token(&item, state)) {
            MD_LOG_ERROR("Unexpected tokenizer failure while reading _atom_site headers");
            return false;
        }
        md_str_t field = md_str_substr(item, 11, SIZE_MAX);
        for (size_t i = 0; i < MD_COUNTOF(atom_site_labels); ++i) {
            if (table[i] != -1) {
                continue;
            }
            if (md_str_eq(field, atom_site_labels[i])) {
                table[i] = (int)num_cols;
                break;
            }
        }
        num_cols += 1;
    }

    for (size_t i = 0; i < MD_COUNTOF(required_atom_site_fields); ++i) {
        if (table[required_atom_site_fields[i]] == -1) {
            MD_LOG_ERROR("Missing required column in _atom_site: "MD_STR_FMT, MD_STR_ARG(atom_site_labels[required_atom_site_fields[i]]));
            return false;
        }
    }

    cols.num_cols = num_cols;
    cols.have_auth_seq_id  = table[ATOM_SITE_AUTH_SEQ_ID]  != -1;
    cols.have_auth_asym_id = table[ATOM_SITE_AUTH_ASYM_ID] != -1;

    // One token per column of a row. A file carries as many _atom_site columns as its producer chose to write, so
    // the row is sized by them rather than by a fixed bound.
    md_str_t* tok = md_alloc(alloc, sizeof(md_str_t) * (num_cols + MMCIF_SPLIT_SLACK));
    md_memset(tok, 0, sizeof(md_str_t) * (num_cols + MMCIF_SPLIT_SLACK));

    // Where the tokenizer has nothing of a line left, the rows which follow are read in blocks (mmcif_rows_from_blocks),
    // up to a line which is not a row of its own: that line is the tokenizer's, given back to it whole, which reads it
    // as it reads any (a comment, a row over lines, what comes after the rows). Rows of more columns than the blocks
    // read are read a line at a time.
    const bool blocks = num_cols <= MMCIF_BLOCK_MAX_COLS;
    bool give_back = false;
    for (;;) {
        bool split = false;
        if (!state->has_peek && md_str_empty(state->line)) {
            if (blocks && !give_back) {
                const mmcif_rows_result_t res = mmcif_rows_from_blocks(atom_entries, &cols, state->reader, alloc, ex);
                if (res == MMCIF_ROWS_CANCELLED) return false;
                if (res == MMCIF_ROWS_END) break;
                give_back = true;
                continue;
            }
            give_back = false;
            md_str_t line;
            if (!md_reader_extract_line(state->reader, &line)) break;      // The end of the file
            split = mmcif_plain_row(tok, num_cols, line);
            if (!split) {
                state->line = line;
                state->pending_start = true;
            }
        }

        if (!split) {
            if (!mmcif_advance_to_next_line(state)) break;
            md_str_t peek;
            if (!mmcif_peek_token(&peek, state)) {
                MD_LOG_ERROR("Unexpected EOF while reading _atom_site entries");
            }
            if (md_str_empty(peek) || mmcif_is_item(peek) || mmcif_is_control(peek)) {
                // End of _atom_site section
                break;
            }
            // From the token peeked (as the line has it) to the end of the line, which is what the tokenizer has left
            if (state->has_peek && state->peek_plain && state->line_start && state->line.ptr) {
                split = mmcif_split_row(tok, num_cols, peek.ptr, state->line.ptr + state->line.len);
                if (split) {
                    state->has_peek = false;
                    state->line = (md_str_t){0};
                }
            }
        }
        for (size_t i = 0; !split && i < num_cols; ++i) {
            if (!mmcif_next_token(&tok[i], state)) {
                MD_LOG_ERROR("Unexpected EOF while reading _atom_site entry");
                return false;
            }
        }
         
        mmcif_atom_site_entry_t entry;
        mmcif_atom_site_entry(&entry, tok, &cols);
        md_array_push(*atom_entries, entry, alloc);
    }

    return true;
}

enum {
    CELL_FIELD_A = 1,
    CELL_FIELD_B = 2,
    CELL_FIELD_C = 4,
    CELL_FIELD_ALPHA = 8,
    CELL_FIELD_BETA = 16,
    CELL_FIELD_GAMMA = 32,
};

static bool mmcif_parse_cell(mmcif_cell_params_t* cell, mmcif_parse_state_t* state) {
    MD_ASSERT(cell);
    MD_ASSERT(state);
    MD_ASSERT(state->in_loop == false);

    uint32_t field_flags = 0;

    while (mmcif_advance_to_next_control_token(state)) {
        md_str_t tok;
        if (!mmcif_peek_token(&tok, state)) {
            MD_LOG_ERROR("Unexpected EOF");
            return false;
        }

        if (md_str_begins_with(tok, MD_STR_LIT("_cell."))) {
            if (!mmcif_next_token(&tok, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure consuming _cell token");
                return false;
            }

            md_str_t val_tok;
            if (!mmcif_next_token(&val_tok, state)) {
                MD_LOG_ERROR("Unexpected EOF consuming _cell value");
                return false;
            }

            double val = parse_float(val_tok);

            if (md_str_eq_cstr(tok, "_cell.angle_alpha")) {
                cell->alpha = val;
                field_flags |= CELL_FIELD_ALPHA;
            } else if (md_str_eq_cstr(tok, "_cell.angle_beta")) {
                cell->beta = val;
                field_flags |= CELL_FIELD_BETA;
            } else if (md_str_eq_cstr(tok, "_cell.angle_gamma")) {
                cell->gamma = val;
                field_flags |= CELL_FIELD_GAMMA;
            } else if (md_str_eq_cstr(tok, "_cell.length_a")) {
                cell->a = val;
                field_flags |= CELL_FIELD_A;
            } else if (md_str_eq_cstr(tok, "_cell.length_b")) {
                cell->b = val;
                field_flags |= CELL_FIELD_B;
            } else if (md_str_eq_cstr(tok, "_cell.length_c")) {
                cell->c = val;
                field_flags |= CELL_FIELD_C;
            }
        } else {
            break;
        }
    }

    return (field_flags == (CELL_FIELD_A | CELL_FIELD_B | CELL_FIELD_C | CELL_FIELD_ALPHA | CELL_FIELD_BETA | CELL_FIELD_GAMMA));
}

static bool mmcif_parse(md_system_t* sys, md_system_state_t* out_state, md_reader_i* reader, md_allocator_i* alloc, const md_executor_t* ex) {
    bool result = false;
    bool atom_site_parsed = false;
    bool cell_parsed = false;
    bool entity_parsed = false;
    bool entity_poly_parsed = false;
    bool assembly_gen_parsed = false;
    bool oper_list_parsed = false;

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    md_array(mmcif_atom_site_entry_t) atom_entries = 0;
    md_array(mmcif_entity_t) entities = 0;
    md_array(mmcif_assembly_gen_t) assembly_gen = 0;
    md_array(mmcif_oper_t) opers = 0;

    md_array_ensure(atom_entries, 1024, temp_arena);
    md_array_ensure(entities, 4, temp_arena);

    mmcif_cell_params_t cell = {0};
    mmcif_parse_state_t state = {
        .reader = reader,
        .sb = md_str_writer_create(temp_arena),
    };

    while (mmcif_advance_to_next_control_token(&state)) {
        md_str_t tok;
        if (!mmcif_peek_token(&tok, &state)) {
            MD_LOG_ERROR("Unexpected EOF");
            goto done;
        }

        if (md_str_eq(tok, MD_STR_LIT("loop_"))) {
            state.in_loop = true;
            mmcif_next_token(&tok, &state);
            continue;
        }

        if (md_str_begins_with(tok, MD_STR_LIT("_atom_site."))) {
            if (!mmcif_parse_atom_site(&atom_entries, &state, temp_arena, ex)) {
                MD_LOG_ERROR("Failed to parse _atom_site section");
                goto done;
            }
            atom_site_parsed = true;
        } else if (md_str_begins_with(tok, MD_STR_LIT("_cell."))) {
            if (!mmcif_parse_cell(&cell, &state)) {
                MD_LOG_ERROR("Failed to parse _cell section");
                goto done;
            }
            cell_parsed = true;
        } else if (md_str_begins_with(tok, MD_STR_LIT("_entity."))) {
            if (!mmcif_parse_entity(&entities, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _entity section");
                goto done;
            }
            entity_parsed = true;
        } else if (md_str_begins_with(tok, MD_STR_LIT("_entity_poly."))) {
            if (!mmcif_parse_entity_poly(&entities, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _entity_poly section");
                goto done;
            }
            entity_poly_parsed = true;
        } else if (md_str_begins_with(tok, MD_STR_LIT("_pdbx_struct_assembly_gen."))) {
            // The assemblies are read as far as they can be: without them the structure is still the structure
            if (!mmcif_parse_assembly_gen(&assembly_gen, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _pdbx_struct_assembly_gen section");
            }
            assembly_gen_parsed = true;
        } else if (md_str_begins_with(tok, MD_STR_LIT("_pdbx_struct_oper_list."))) {
            if (!mmcif_parse_oper_list(&opers, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _pdbx_struct_oper_list section");
            }
            oper_list_parsed = true;
        } else {
            mmcif_next_token(&tok, &state);
        }

        if (atom_site_parsed && cell_parsed && entity_parsed && entity_poly_parsed && assembly_gen_parsed && oper_list_parsed) {
            break;
        }

        state.in_loop = false;
    }

    if (entity_parsed) {
        // The entities are the part types the file gives, added first so that the instances can name them
        size_t num_entities = md_array_size(entities);
        for (size_t i = 0; i < num_entities; ++i) {
            md_part_type_data_push(&sys->part.type, entities[i].id, md_part_type_flags_set_kind(MD_PART_TYPE_FLAG_NONE, mmcif_entity_kind(&entities[i])), entities[i].description, alloc);
        }
    }
    
    // Populate molecule from parsed data
    if (atom_site_parsed) {
        size_t num_atoms = md_array_size(atom_entries);
        size_t reserve_size = MD_ALIGN_TO(num_atoms, 16);

        md_array_ensure(out_state->xyz, reserve_size, out_state->alloc);
        md_array_ensure(sys->atom.type_idx, reserve_size, alloc);
        md_array_ensure(sys->atom.flags, reserve_size, alloc);

        md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, alloc); // Ensure that index 0 is always unknown
        md_hashmap32_t type_map = { .allocator = temp_arena };

        // Gathered alongside the coordinates rather than from atom_entries afterwards, because the
        // loop SKIPS alternate locations - so entry i and atom i are not the same index, and a
        // column built from the entries would be shifted against the atoms it describes.
        float* occupancy      = (float*)md_alloc(temp_arena, sizeof(float) * reserve_size);
        float* b_iso          = (float*)md_alloc(temp_arena, sizeof(float) * reserve_size);
        float* formal_charge  = (float*)md_alloc(temp_arena, sizeof(float) * reserve_size);

        // Of the atom before: the fields of its residue and its instance as words (the entries' fields are zero past
        // their characters, so the same characters are the same words), and its entity. A residue begins where the
        // residue's differ, an instance where the asym's does.
        bool     have_prev = false;
        uint64_t prev_asym = 0, prev_comp = 0, prev_entity = 0;
        int      prev_seq = 0;
        md_part_type_idx_t entity_idx = -1;

        // The type of an atom by its name and element symbol, of which a structure has few (the atoms of each residue
        // of a kind have the same): looked up once for each, by their words, before the element is found by its
        // symbol and the type by its name
        typedef struct {
            uint64_t name;
            uint32_t symbol;
            md_particle_type_idx_t type;
        } mmcif_type_seen_t;
        md_array(mmcif_type_seen_t) seen = 0;
        md_hashmap32_t seen_map = { .allocator = temp_arena };


        for (size_t i = 0; i < num_atoms; ++i) {
            // Ignore alt loc entries unless 'A'
            if (atom_entries[i].label_alt_id != '.' &&
                atom_entries[i].label_alt_id != 'A') continue;

            const mmcif_atom_site_entry_t* entry = &atom_entries[i];
            uint64_t asym, comp, entity, name;
            uint32_t symbol_word = 0;
            MD_STATIC_ASSERT(sizeof(entry->label_asym_id) == 8 && sizeof(entry->label_comp_id) == 8 && sizeof(entry->label_entity_id) == 8 &&
                          sizeof(entry->label_atom_id) == 8 && sizeof(entry->type_symbol) == 3, "Fields of 8 characters (a word), the symbol of 3");
            md_memcpy(&asym,   entry->label_asym_id,   8);
            md_memcpy(&comp,   entry->label_comp_id,   8);
            md_memcpy(&entity, entry->label_entity_id, 8);
            md_memcpy(&name,   entry->label_atom_id,   8);
            md_memcpy(&symbol_word, entry->type_symbol, 3);
            const int seq_id = entry->label_seq_id != INT32_MIN ? entry->label_seq_id : entry->auth_seq_id;

            if (!have_prev || entity != prev_entity) {
                entity_idx = md_system_part_type_find_by_id(sys, md_str_from_cstrn(entry->label_entity_id, sizeof(entry->label_entity_id)));
            }

            md_particle_type_idx_t particle_type_idx = 0;
            {
                const uint64_t key = (name ^ ((uint64_t)symbol_word << 40 | symbol_word)) * 0x9E3779B97F4A7C15ull;
                const uint32_t* s_idx = md_hashmap_get(&seen_map, key);
                if (s_idx && seen[*s_idx].name == name && seen[*s_idx].symbol == symbol_word) {
                    particle_type_idx = seen[*s_idx].type;
                } else {
                    md_str_t symbol = md_str_from_cstrn(entry->type_symbol, sizeof(entry->type_symbol));
                    md_str_t atom_id = md_str_from_cstrn(entry->label_atom_id, sizeof(entry->label_atom_id));
                    md_atomic_number_t atomic_number = md_atomic_number_from_symbol(symbol, true);
                    particle_type_idx = md_util_particle_type_find_or_add(sys, &type_map, atom_id, atomic_number);
                    if (!s_idx) {
                        md_hashmap_add(&seen_map, key, (uint32_t)md_array_size(seen));
                        md_array_push(seen, ((mmcif_type_seen_t){ name, symbol_word, particle_type_idx }), temp_arena);
                    }
                }
            }

            const bool new_inst = !have_prev || asym != prev_asym;
            if (new_inst || comp != prev_comp || seq_id != prev_seq) {
                md_str_t comp_id = md_str_from_cstrn(entry->label_comp_id, MD_COUNTOF(entry->label_comp_id));
                md_str_t label_asym_id = md_str_from_cstrn(entry->label_asym_id, sizeof(entry->label_asym_id));
                md_str_t auth_asym_id = md_str_from_cstrn(entry->auth_asym_id, sizeof(entry->auth_asym_id));
                // The chain id is the author's (md_chain_data_t). A file which does not give it calls its chains by their
                // labels only, and those are then its chain ids.
                if (md_str_empty(auth_asym_id) || md_str_eq(auth_asym_id, MD_STR_LIT("?")) || md_str_eq(auth_asym_id, MD_STR_LIT("."))) {
                    auth_asym_id = label_asym_id;
                }

                // --- New residue boundary ---
                // Start residue (store atom offset before adding any residue atoms)
                md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, alloc);
                md_array_push(sys->residue.name, md_label_make(comp_id), alloc);
                md_array_push(sys->residue.seq_id, seq_id, alloc);
                md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, alloc);

                // Instance handling. An instance is what the file gives: an asym, a polymer chain or a branched
                // glycan. The asym of a non-polymer entity may hold many molecules (the waters of a chain): each is an
                // instance of its own, sharing the asym id.
                const md_part_kind_t entity_kind = entity_idx != -1 ? md_system_part_type_kind(sys, entity_idx) : MD_PART_KIND_UNKNOWN;
                const bool one_per_residue = entity_kind == MD_PART_KIND_WATER || entity_kind == MD_PART_KIND_NON_POLYMER;
                if (sys->instance.count == 0 || new_inst || one_per_residue) {
                    md_array_push(sys->instance.id, md_label_make(label_asym_id), alloc);
                    md_array_push(sys->instance.auth_id, md_label_make(auth_asym_id), alloc);
                    md_array_push(sys->instance.res_offset, (uint32_t)sys->residue.count, alloc);
                    md_array_push(sys->instance.part_type_idx, entity_idx, alloc);
                    md_array_push(sys->instance.mol_type_idx, -1, alloc);
                    sys->instance.count += 1;
                }
                sys->residue.count += 1;
            }
 
            have_prev = true;
            prev_asym = asym;
            prev_comp = comp;
            prev_entity = entity;
            prev_seq = seq_id;
 
            md_array_push_no_grow(out_state->xyz, md_f32x3_set(atom_entries[i].x, atom_entries[i].y, atom_entries[i].z));
            md_array_push_no_grow(sys->atom.type_idx, particle_type_idx);
            md_array_push_no_grow(sys->atom.flags, MD_ATOM_FLAG_NONE);
            if (occupancy && b_iso && formal_charge) {
                occupancy[sys->atom.count]     = atom_entries[i].occupancy;
                b_iso[sys->atom.count]         = atom_entries[i].b_iso_or_equiv;
                formal_charge[sys->atom.count] = atom_entries[i].formal_charge;
            }
            sys->atom.count += 1;
        }

        MD_ASSERT(md_array_size(out_state->xyz) == sys->atom.count);
        out_state->num_atoms = sys->atom.count;

        // The same three paths and units md_pdb publishes, because they are the same quantities out
        // of the same file family and a consumer must not have to ask which reader produced a
        // structure. A column the file omits is all NAN, which is uniform, so it publishes nothing.
        if (occupancy && b_iso && formal_charge) {
            md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/occupancy"),     md_unit_none(),                     1, occupancy,     sys->atom.count);
            md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/b_factor"),      md_unit_pow(md_unit_angstrom(), 2), 1, b_iso,         sys->atom.count);
            md_attributes_publish_atom_column(&sys->attributes, MD_STR_LIT("atom/formal_charge"), md_unit_none(),                     1, formal_charge, sys->atom.count);

            // The formal charges given by the file are the system's (md_atom_data_t) until perception replaces them,
            // which keeps them by reading the column above. Built with fast math, so NAN (no value) is told by its
            // bits and not by q != q.
            bool any_charge = false;
            for (size_t i = 0; i < sys->atom.count && !any_charge; ++i) any_charge = !float_bits_nan(formal_charge[i]);
            if (any_charge) {
                md_array_resize(sys->atom.formal_charge, sys->atom.count, alloc);
                for (size_t i = 0; i < sys->atom.count; ++i) {
                    const float q = formal_charge[i];
                    sys->atom.formal_charge[i] = float_bits_nan(q) ? 0 : (int8_t)MD_CLAMP(lroundf(q), -8, 8);
                }
            }
        }

        if (sys->residue.atom_offset) {
            md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, alloc);    // Final sentinel
        }
        if (sys->instance.res_offset) {
            md_array_push(sys->instance.res_offset, (uint32_t)sys->residue.count, alloc);     // Final sentinel
        }
    }

    out_state->unitcell = md_unitcell_none();
    if (cell_parsed) {
        out_state->unitcell = md_unitcell_from_extent_and_angles(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
    }

    // The bonds, what the residues are by them and the chains are md_util_system_infer's (MD_UTIL_INFER_BOND_BIT,
    // _RESIDUE_BIT, _PART_BIT): the loader gives what the file says.

    // The assemblies name the instances by their asym ids
    mmcif_assemblies(sys, assembly_gen, opers, temp_arena);

    result = sys->atom.count > 0;

done:
    md_temp_end(temp_scope);
    return result;
}

bool md_mmcif_system_init_from_str(md_system_t* sys, md_system_state_t* state, md_str_t str, const md_executor_t* ex) {
    MD_ASSERT(sys);
    MD_ASSERT(state);

    if (md_str_empty(str)) {
        MD_LOG_ERROR("Input string is empty");
        return false;
    }

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

    md_reader_i reader = md_reader_from_str(str);
    bool result = mmcif_parse(sys, state, &reader, sys->alloc, ex);
    state->num_atoms = sys->atom.count;
    return result;
}

bool md_mmcif_system_init_from_file(md_system_t* sys, md_system_state_t* state, md_str_t filename, const md_executor_t* ex) {
    MD_ASSERT(sys);
    MD_ASSERT(state);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Could not open file '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }

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

    // The reader's buffer in a temp scope below the one of the load (mmcif_parse): blocks of the rows of _atom_site, a
    // few passes of them (mmcif_rows_from_blocks)
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    const size_t cap = MD_MAX((size_t)MD_TEXT_BLOCK_SIZE, 2 * md_text_pass_size(md_executor_num_threads(ex)));

    md_file_reader_t reader;
    md_file_reader_init(&reader, file, cap, md_temp_allocator(temp));
    bool result = mmcif_parse(sys, state, &reader.r, sys->alloc, ex);
    state->num_atoms = sys->atom.count;

    md_file_reader_close(&reader);
    md_temp_end(temp);
    md_file_close(&file);

    return result;
}
