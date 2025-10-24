#include "md_mmcif.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <md_util.h>
#include <md_system.h>

#define BAKE(str) {str, sizeof(str)-1}

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

static const str_t atom_site_labels[] = {
    BAKE("group_PDB"),
    BAKE("id"),
    BAKE("type_symbol"),
    BAKE("label_atom_id"),
    BAKE("label_alt_id"),
    BAKE("label_comp_id"),
    BAKE("label_asym_id"),
    BAKE("label_entity_id"),
    BAKE("label_seq_id"),
    BAKE("pdbx_PDB_ins_code"),
    BAKE("Cartn_x"),
    BAKE("Cartn_y"),
    BAKE("Cartn_z"),
    BAKE("occupancy"),
    BAKE("B_iso_or_equiv"),
    BAKE("pdbx_formal_charge"),
    BAKE("auth_seq_id"),
    BAKE("auth_comp_id"),
    BAKE("auth_asym_id"),
    BAKE("auth_atom_id"),
    BAKE("pdbx_PDB_model_num"),
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
//    int8_t pdbx_formal_charge;

    float x, y, z;
    int label_seq_id;
    int auth_seq_id;
//    float occupancy;
//    float b_iso_or_equiv;
    char label_asym_id[4];
    char auth_asym_id[4];
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

static const str_t entity_type_str[] = {
    BAKE(""),
    BAKE("branched"),
    BAKE("macrolide"),
    BAKE("non-polymer"),
    BAKE("polymer"),
    BAKE("water"),
};

typedef int mmcif_entity_type_t;

static inline mmcif_entity_type_t mmcif_entity_type_from_str(str_t type) {
    for (size_t i = 0; i < ARRAY_SIZE(entity_type_str); ++i) {
        if (str_eq(type, entity_type_str[i])) {
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

static const str_t entity_poly_type_str[] = {
    BAKE(""),
    BAKE("cyclic-pseudo-peptide"),
    BAKE("other"),
    BAKE("peptide nucleic acid"),
    BAKE("polydeoxyribonucleotide"),
    BAKE("polydeoxyribonucleotide/polyribonucleotide hybrid"),
    BAKE("polypeptide(D)"),
    BAKE("polypeptide(L)"),
    BAKE("polyribonucleotide"),
};

typedef int mmcif_entity_poly_type_t;

static inline mmcif_entity_poly_type_t mmcif_entity_poly_type_from_str(str_t type) {
    for (size_t i = 0; i < ARRAY_SIZE(entity_poly_type_str); ++i) {
        if (str_eq(type, entity_poly_type_str[i])) {
            return (mmcif_entity_poly_type_t)i;
        }
    }
    return ENTITY_POLY_TYPE_UNKNOWN;
}

typedef struct {
    str_t id;
    mmcif_entity_type_t type;
    mmcif_entity_poly_type_t poly_type;
    str_t description;
} mmcif_entity_t;

static inline mmcif_entity_t* mmcif_entity_find(md_array(mmcif_entity_t) entities, str_t id) {
    for (mmcif_entity_t* it = md_array_beg(entities); it != md_array_end(entities); ++it) {
        if (str_eq(it->id, id)) return it;
    }
    return NULL;
}

static inline mmcif_entity_t* mmcif_entity_find_or_create(md_array(mmcif_entity_t)* entities, str_t id, md_allocator_i* alloc) {
    mmcif_entity_t* it = mmcif_entity_find(*entities, id);
    if (it) return it;

    mmcif_entity_t ent = {
        .id = str_copy(id, alloc),
    };
    md_array_push(*entities, ent, alloc);
    return md_array_last(*entities);
}

static inline int parse_int_with_default(str_t tok, int def) {
    if (tok.len == 0 || tok.ptr[0] == '.' || tok.ptr[0] == '?') {
        return def;
    }
    return (int)parse_int(tok);
}

typedef struct {
    md_buffered_reader_t* reader;
    md_strb_t sb;
    str_t line;
    str_t peek;
    bool has_peek;
    bool line_start;
    bool in_loop;
} mmcif_parse_state_t;

// Returns true if token was extracted,
// False if EOF / stream
static bool mmcif_fetch_token(str_t* out_tok, mmcif_parse_state_t* s) {
    str_t tok = {0};
    s->line_start = false;
    while (true) {
        if (str_empty(s->line)) {
            if (!md_buffered_reader_extract_line(&s->line, s->reader)) {
                return false;
            }
            s->line_start = true;
        }
        if (!extract_token(&tok, &s->line) || str_len(tok) == 0 || (tok.len > 0 && tok.ptr[0] == '#')) {
            s->line = (str_t){0};
            continue;
        }

        char first_char = tok.ptr[0];
        if (first_char == '\'' || first_char == '\"') {
            char delim = tok.ptr[0];
            if (tok.len > 1 && tok.ptr[tok.len-1] == delim) {
                tok.ptr +=1;
                tok.len -=2;
            } else {
                size_t loc = 0;
                if (!str_find_char(&loc, s->line, delim)) {
                    MD_LOG_ERROR("Unmatched string");
                    return false;
                }
                tok.ptr  += 1;
                tok.len  += loc;
                s->line.ptr += (loc+1);
                s->line.len -= (loc+1);
            }
        } else if (first_char == ';') {
            md_strb_reset(&s->sb);
            md_strb_push_str(&s->sb, str_substr(tok, 1, SIZE_MAX));
            md_strb_push_str(&s->sb, s->line);

            bool term = false;
            while (md_buffered_reader_extract_line(&s->line, s->reader)) {
                if (s->line.len > 0 && s->line.ptr[0] == ';') {
                    term = true;
                    s->line = (str_t) {0};
                    break;
                }
                md_strb_push_char(&s->sb, '\n');
                md_strb_push_str (&s->sb, s->line);
            }
            if (!term) {
                MD_LOG_ERROR("Unterminated multiline string");
                return false;
            }
            tok = md_strb_to_str(s->sb);
        }
        break;
    }

    *out_tok = tok;
    return true;
}

static inline bool mmcif_next_token(str_t* tok, mmcif_parse_state_t* state) {
    ASSERT(tok);
    ASSERT(state);
    if (state->has_peek) {
        *tok = state->peek;
        state->has_peek = false;
        return true;
    }
    return mmcif_fetch_token(tok, state);
}

static inline bool mmcif_peek_token(str_t* tok, mmcif_parse_state_t* state) {
    ASSERT(tok);
    ASSERT(state);
    if (state->has_peek) {
        *tok = state->peek;
        return true;
    }
    if (!mmcif_fetch_token(&state->peek, state)) {
        return false;
    }
    state->has_peek = true;
    *tok = state->peek;
    return true;
}

// This is a bit wierd, since were mixing lines and token semantics.
// But will return the current line within the state that is about to be tokenized
static inline str_t mmcif_cur_line(mmcif_parse_state_t* s) {
    ASSERT(s);     
    return s->line;
}

// Advances to the next line, which is not part of a multiline string token
static inline bool mmcif_advance_to_next_line(mmcif_parse_state_t* s) {
    ASSERT(s);
    str_t tok = {0};
    while (mmcif_peek_token(&tok, s)) {
        if (s->line_start) return true;
        mmcif_next_token(&tok, s);
    }
    return false; // EOF
}

// Advances to a new line which starts with a valid control or item token
static inline bool mmcif_advance_to_next_control_token(mmcif_parse_state_t* s) {
    ASSERT(s);

    str_t tok = {0};
    while (mmcif_peek_token(&tok, s)) {
        if (tok.len > 0 && s->line_start) {
            const char c0 = tok.ptr[0];

            // Category and control keywords:
            if (c0 == '_' || 
                str_eq(tok, STR_LIT("loop_")) ||
                str_eq_cstr_n(tok, "data_", 5) ||
                str_eq_cstr_n(tok, "save_", 5)) 
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
    str_t* headers;
    str_t* values;
} mmcif_section_t;

// Help accessors for data in sections
static inline str_t mmcif_section_header(const mmcif_section_t* sec, size_t col) {
    if (!sec || col >= sec->num_fields) return (str_t){0};
    return sec->headers[col];
}

static inline str_t mmcif_section_value(const mmcif_section_t* sec, size_t row, size_t col) {
    if (!sec || col >= sec->num_fields || row >= sec->num_rows) return (str_t){0};
    return sec->values[row * sec->num_fields + col];
}

static inline bool mmcif_is_item(str_t t) { return t.len && t.ptr[0] == '_'; }
static inline bool mmcif_is_control(str_t t) {
    return (t.len >= 5 &&
        (str_eq(t,STR_LIT("loop_")) ||
         str_eq(t,STR_LIT("data_")) ||
         str_eq(t,STR_LIT("save_")) ||
         str_eq(t,STR_LIT("stop_"))));
}

static str_t mmcif_category(str_t item) {
    if (!mmcif_is_item(item)) return (str_t){0};
    str_t tail = str_substr(item,1,SIZE_MAX);
    size_t dot;
    if (str_find_char(&dot, tail, '.')) return (str_t){tail.ptr, dot};
    return (str_t){0};
}

// This is an internal procedure, we do not care about cleaning anything up
static bool mmcif_parse_section(mmcif_section_t* sec, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    ASSERT(sec);
    ASSERT(state);
    ASSERT(alloc);
    MEMSET(sec, 0, sizeof(*sec));

    if (state->in_loop) {
        // --------------- Collect headers ---------------
        str_t t;
        while (mmcif_peek_token(&t, state) && mmcif_is_item(t)) {
            if (!mmcif_next_token(&t, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure while reading loop headers");
                return false;
            }
            md_array_push(sec->headers, str_copy(t, alloc), alloc);
        }
        size_t ncols = md_array_size(sec->headers);
        if (ncols == 0) {
            MD_LOG_ERROR("Loop section without headers");
            return false;
        }

        // --------------- Collect rows ---------------
        size_t values_read = 0;
        while (true) {
            str_t look;
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
                str_t v;
                if (!mmcif_next_token(&v, state)) {
                    MD_LOG_ERROR("Unexpected EOF in middle of loop row");
                    return false;
                }
                if (mmcif_is_control(v) || mmcif_is_item(v)) {
                    // Structural token unexpectedly inside a row -> error
                    MD_LOG_ERROR("Premature loop termination inside row (got %zu/%zu values)", got, ncols);
                    return false;
                }
                md_array_push(sec->values, str_copy(v, alloc), alloc);
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
        str_t first;
        if (!mmcif_peek_token(&first, state) || !mmcif_is_item(first)) {
            // Nothing to parse (not an error? choose: treat as failure to signal absence)
            return false;
        }
        str_t cat = mmcif_category(first);
        if (str_empty(cat)) {
            MD_LOG_ERROR("Malformed item (missing category segment): " STR_FMT, STR_ARG(first));
            return false;
        }

        while (true) {
            str_t key;
            if (!mmcif_peek_token(&key, state)) break;
            if (mmcif_is_control(key) || !mmcif_is_item(key)) break;

            // Category change ends this section
            str_t this_cat = mmcif_category(key);
            if (!str_empty(this_cat) && !str_eq(this_cat, cat)) break;

            // Consume key
            if (!mmcif_next_token(&key, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure consuming item key");
                return false;
            }
            md_array_push(sec->headers, str_copy(key, alloc), alloc);

            // Value (optional)
            str_t val;
            if (!mmcif_peek_token(&val, state) || mmcif_is_control(val) || mmcif_is_item(val)) {
                // Missing value -> store empty
                return false;
            }
            if (!mmcif_next_token(&val, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure consuming item value");
                return false;
            }
            md_array_push(sec->values, str_copy(val, alloc), alloc);
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

static inline int mmcif_section_find_col_idx(const mmcif_section_t* sec, str_t field_label) {
    for (size_t i = 0; i < sec->num_fields; ++i) {
        if (str_eq(field_label, sec->headers[i])) {
            return (int)i;
        }
    }
    return -1;
}

static bool mmcif_parse_entity(md_array(mmcif_entity_t)* entities, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    ASSERT(entities);
    ASSERT(state);
    ASSERT(alloc);

    md_allocator_i* temp_alloc = md_get_temp_allocator();
    size_t temp_pos = md_temp_get_pos();

    bool success = false;
    
    mmcif_section_t sec = {0};
    if (!mmcif_parse_section(&sec, state, temp_alloc)) {
        MD_LOG_ERROR("Failed to parse _entity section");
        goto done;
    }

    int id_col   = mmcif_section_find_col_idx(&sec, STR_LIT("_entity.id"));
    int type_col = mmcif_section_find_col_idx(&sec, STR_LIT("_entity.type"));
    int desc_col = mmcif_section_find_col_idx(&sec, STR_LIT("_entity.pdbx_description"));

    if (id_col == -1) {
        MD_LOG_ERROR("Missing entity field 'id'");
        return false;
    }
    if (type_col == -1) {
        MD_LOG_ERROR("Missing entity field 'type'");
        return false;
    }
    if (desc_col == -1) {
        MD_LOG_ERROR("Missing entity field 'pdbx_description'");
        return false;
    }

    for (size_t i = 0; i < sec.num_rows; ++i) {
        str_t id = mmcif_section_value(&sec, i, id_col);
        if (!str_empty(id)) {
            mmcif_entity_t* ent = mmcif_entity_find_or_create(entities, id, alloc);
            ASSERT(ent);
            ent->type = mmcif_entity_type_from_str(mmcif_section_value(&sec, i, type_col));
            ent->description = str_copy(mmcif_section_value(&sec, i, desc_col), alloc);
        }
    }

    success = true;
done:
    md_temp_set_pos_back(temp_pos);
    return success;
}

// We only populate the entity_poly type within the given entities
static bool mmcif_parse_entity_poly(md_array(mmcif_entity_t)* entities, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    ASSERT(entities);
    ASSERT(state);
    ASSERT(alloc);

    bool success = false;

    md_allocator_i* temp_alloc = md_get_temp_allocator();
    size_t temp_pos = md_temp_get_pos();

    mmcif_section_t sec = {0};
    if (!mmcif_parse_section(&sec, state, temp_alloc)) {
        MD_LOG_ERROR("Failed to parse _entity_poly section");
        goto done;
    }

    int id_col   = mmcif_section_find_col_idx(&sec, STR_LIT("_entity_poly.entity_id"));
    int type_col = mmcif_section_find_col_idx(&sec, STR_LIT("_entity_poly.type"));

    if (id_col == -1) {
        MD_LOG_ERROR("Missing entity field 'entity_id'");
        return false;
    }
    if (type_col == -1) {
        MD_LOG_ERROR("Missing entity field 'type'");
        return false;
    }

    for (size_t i = 0; i < sec.num_rows; ++i) {
        str_t id = mmcif_section_value(&sec, i, id_col);
        if (!str_empty(id)) {
            mmcif_entity_t* ent = mmcif_entity_find_or_create(entities, id, alloc);
            ASSERT(ent);
            ent->poly_type = mmcif_entity_poly_type_from_str(mmcif_section_value(&sec, i, type_col));
        }
    }

    success = true;
done:
    md_temp_set_pos_back(temp_pos);
    return success;
}

static bool mmcif_parse_atom_site(md_array(mmcif_atom_site_entry_t)* atom_entries, mmcif_parse_state_t* state, md_allocator_i* alloc) {
    ASSERT(atom_entries);
    ASSERT(state);
    ASSERT(state->in_loop == true);

    int table[ATOM_SITE_COUNT];
    MEMSET(table, -1, sizeof(table));
    size_t num_cols = 0;

    bool success = false;

    // Read category headers
    str_t item;
    while (mmcif_peek_token(&item, state) && mmcif_is_item(item) && str_eq_cstr_n(item, "_atom_site.", 11)) {
        if (!mmcif_next_token(&item, state)) {
            MD_LOG_ERROR("Unexpected tokenizer failure while reading _atom_site headers");
            return false;
        }
        str_t field = str_substr(item, 11, SIZE_MAX);
        for (size_t i = 0; i < ARRAY_SIZE(atom_site_labels); ++i) {
            if (table[i] != -1) {
                continue;
            }
            if (str_eq(field, atom_site_labels[i])) {
                table[i] = (int)num_cols;
                break;
            }
        }
        num_cols += 1;
    }

    for (size_t i = 0; i < ARRAY_SIZE(required_atom_site_fields); ++i) {
        if (table[required_atom_site_fields[i]] == -1) {
            MD_LOG_ERROR("Missing required column in _atom_site: "STR_FMT, STR_ARG(atom_site_labels[required_atom_site_fields[i]]));
            return false;
        }
    }

    bool have_auth_seq_id  = table[ATOM_SITE_AUTH_SEQ_ID]  != -1;
    bool have_auth_asym_id = table[ATOM_SITE_AUTH_ASYM_ID] != -1;

    str_t tok[64] = {0};
    while (mmcif_advance_to_next_line(state)) {
        str_t peek;
        if (!mmcif_peek_token(&peek, state)) {
            MD_LOG_ERROR("Unexpected EOF while reading _atom_site entries");
        }
        if (str_empty(peek) || mmcif_is_item(peek) || mmcif_is_control(peek)) {
            // End of _atom_site section
            break;
        }

        for (size_t i = 0; i < num_cols; ++i) {
            if (!mmcif_next_token(&tok[i], state)) {
                MD_LOG_ERROR("Unexpected EOF while reading _atom_site entry");
                return false;
            }
        }
         
        str_t group_PDB       = tok[table[ATOM_SITE_GROUP_PDB]];
        str_t type_symbol     = tok[table[ATOM_SITE_TYPE_SYMBOL]];
        str_t label_atom_id   = tok[table[ATOM_SITE_LABEL_ATOM_ID]];
        str_t label_alt_id    = tok[table[ATOM_SITE_LABEL_ALT_ID]];
        str_t label_asym_id   = tok[table[ATOM_SITE_LABEL_ASYM_ID]];
        str_t label_comp_id   = tok[table[ATOM_SITE_LABEL_COMP_ID]];
        str_t label_entity_id = tok[table[ATOM_SITE_LABEL_ENTITY_ID]];
        str_t label_seq_id    = tok[table[ATOM_SITE_LABEL_SEQ_ID]];
        str_t auth_seq_id     = have_auth_seq_id ? tok[table[ATOM_SITE_AUTH_SEQ_ID]] : (str_t){0};
        str_t auth_asym_id    = have_auth_asym_id ? tok[table[ATOM_SITE_AUTH_ASYM_ID]] : (str_t){0};
        str_t cartn_x         = tok[table[ATOM_SITE_CARTN_X]];
        str_t cartn_y         = tok[table[ATOM_SITE_CARTN_Y]];
        str_t cartn_z         = tok[table[ATOM_SITE_CARTN_Z]];

        const int default_value = INT32_MIN;

        mmcif_atom_site_entry_t entry = {
            .x = (float)parse_float(cartn_x),
            .y = (float)parse_float(cartn_y),
            .z = (float)parse_float(cartn_z),
            .label_seq_id = parse_int_with_default(label_seq_id, default_value),
            .auth_seq_id = parse_int_with_default(auth_seq_id, default_value),
            .label_alt_id   = label_alt_id.len > 0 ? label_alt_id.ptr[0] : '.',
        };

        str_copy_to_char_buf(entry.label_atom_id,   sizeof(entry.label_atom_id),    label_atom_id);
        str_copy_to_char_buf(entry.label_asym_id,   sizeof(entry.label_asym_id),    label_asym_id);
        str_copy_to_char_buf(entry.label_comp_id,   sizeof(entry.label_comp_id),    label_comp_id);
        str_copy_to_char_buf(entry.label_entity_id, sizeof(entry.label_entity_id),  label_entity_id);
        str_copy_to_char_buf(entry.type_symbol,     sizeof(entry.type_symbol),      type_symbol);
        str_copy_to_char_buf(entry.auth_asym_id,    sizeof(entry.auth_asym_id),     auth_asym_id);

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
    ASSERT(cell);
    ASSERT(state);
    ASSERT(state->in_loop == false);

    uint32_t field_flags = 0;

    while (mmcif_advance_to_next_control_token(state)) {
        str_t tok;
        if (!mmcif_peek_token(&tok, state)) {
            MD_LOG_ERROR("Unexpected EOF");
            return false;
        }

        if (str_eq_cstr_n(tok, "_cell.", 6)) {
            if (!mmcif_next_token(&tok, state)) {
                MD_LOG_ERROR("Unexpected tokenizer failure consuming _cell token");
                return false;
            }

            str_t val_tok;
            if (!mmcif_next_token(&val_tok, state)) {
                MD_LOG_ERROR("Unexpected EOF consuming _cell value");
                return false;
            }

            double val = parse_float(val_tok);

            if (str_eq_cstr(tok, "_cell.angle_alpha")) {
                cell->alpha = val;
                field_flags |= CELL_FIELD_ALPHA;
            } else if (str_eq_cstr(tok, "_cell.angle_beta")) {
                cell->beta = val;
                field_flags |= CELL_FIELD_BETA;
            } else if (str_eq_cstr(tok, "_cell.angle_gamma")) {
                cell->gamma = val;
                field_flags |= CELL_FIELD_GAMMA;
            } else if (str_eq_cstr(tok, "_cell.length_a")) {
                cell->a = val;
                field_flags |= CELL_FIELD_A;
            } else if (str_eq_cstr(tok, "_cell.length_b")) {
                cell->b = val;
                field_flags |= CELL_FIELD_B;
            } else if (str_eq_cstr(tok, "_cell.length_c")) {
                cell->c = val;
                field_flags |= CELL_FIELD_C;
            }
        } else {
            break;
        }
    }

    return (field_flags == (CELL_FIELD_A | CELL_FIELD_B | CELL_FIELD_C | CELL_FIELD_ALPHA | CELL_FIELD_BETA | CELL_FIELD_GAMMA));
}

static bool mmcif_parse(md_system_t* sys, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    bool atom_site_parsed = false;
    bool cell_parsed = false;
    bool entity_parsed = false;
    bool entity_poly_parsed = false;

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    md_array(mmcif_atom_site_entry_t) atom_entries = 0;
    md_array(mmcif_entity_t) entities = 0;

    md_array_ensure(atom_entries, 1024, temp_arena);
    md_array_ensure(entities, 4, temp_arena);

    MEMSET(sys, 0, sizeof(md_system_t));

    mmcif_cell_params_t cell = {0};
    mmcif_parse_state_t state = {
        .reader = reader,
        .sb = md_strb_create(temp_arena),
    };

    while (mmcif_advance_to_next_control_token(&state)) {
        str_t tok;
        if (!mmcif_peek_token(&tok, &state)) {
            MD_LOG_ERROR("Unexpected EOF");
            return false;
        }

        if (str_eq(tok, STR_LIT("loop_"))) {
            state.in_loop = true;
            mmcif_next_token(&tok, &state);
            continue;
        }

        if (str_eq_cstr_n(tok, "_atom_site.", 11)) {
            if (!mmcif_parse_atom_site(&atom_entries, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _atom_site section");
                return false;
            }
            atom_site_parsed = true;
        } else if (str_eq_cstr_n(tok, "_cell.", 6)) {
            if (!mmcif_parse_cell(&cell, &state)) {
                MD_LOG_ERROR("Failed to parse _cell section");
                return false;
            }
            cell_parsed = true;
        } else if (str_eq_cstr_n(tok, "_entity.", 8)) {
            if (!mmcif_parse_entity(&entities, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _entity section");
                return false;
            }
            entity_parsed = true;
        } else if (str_eq_cstr_n(tok, "_entity_poly.", 13)) {
            if (!mmcif_parse_entity_poly(&entities, &state, temp_arena)) {
                MD_LOG_ERROR("Failed to parse _entity_poly section");
                return false;
            }
            entity_poly_parsed = true;
        } else {
            mmcif_next_token(&tok, &state);
        }

        if (atom_site_parsed && cell_parsed && entity_parsed && entity_poly_parsed) {
            break;
        }

        state.in_loop = false;
    }

    if (entity_parsed) {
        // Add entities first so that we can reference them when adding components
        size_t num_entities = md_array_size(entities);
        md_array_ensure(sys->entity.id, num_entities, alloc);
        md_array_ensure(sys->entity.flags, num_entities, alloc);
        md_array_ensure(sys->entity.description, num_entities, alloc);
        sys->entity.count = num_entities;

        for (size_t i = 0; i < num_entities; ++i) {
            // Determine entity flags
            md_flags_t flags = 0;

            str_t desc = entities[i].description;
            bool ion = str_find_str(NULL, desc, STR_LIT("ion")) ||
                       str_find_str(NULL, desc, STR_LIT("ION")) ||
                       str_find_str(NULL, desc, STR_LIT("Ion"));

            if (entities[i].type == ENTITY_TYPE_POLYMER) {
                flags |= MD_FLAG_POLYMER;
                if (entities[i].poly_type == ENTITY_POLY_TYPE_POLYPEPTIDE_L) {
                    flags |= MD_FLAG_AMINO_ACID | MD_FLAG_ISOMER_L;
                } else if (entities[i].poly_type == ENTITY_POLY_TYPE_POLYPEPTIDE_D) {
                    flags |= MD_FLAG_AMINO_ACID | MD_FLAG_ISOMER_D;
                } else if (entities[i].poly_type == ENTITY_POLY_TYPE_POLYRIBONUCLEOTIDE ||
                    entities[i].poly_type == ENTITY_POLY_TYPE_POLYDEOXYRIBONUCLEOTIDE ||
                    entities[i].poly_type == ENTITY_POLY_TYPE_POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID) {
                    flags |= MD_FLAG_NUCLEOTIDE;
                }
            } else if (entities[i].type == ENTITY_TYPE_WATER) {
                flags |= MD_FLAG_WATER;
            } else if (ion) {
                flags |= MD_FLAG_ION;
            } else {
                flags |= MD_FLAG_HETERO;
            }

            md_array_push_no_grow(sys->entity.id, make_label(entities[i].id));
            md_array_push_no_grow(sys->entity.flags, flags);
            md_array_push_no_grow(sys->entity.description, str_copy(entities[i].description, alloc));
        }
    }

    // Populate molecule from parsed data
    if (atom_site_parsed) {
        size_t num_atoms = md_array_size(atom_entries);
        size_t reserve_size = ALIGN_TO(num_atoms, 16);

        md_array_ensure(sys->atom.x, reserve_size, alloc);
        md_array_ensure(sys->atom.y, reserve_size, alloc);
        md_array_ensure(sys->atom.z, reserve_size, alloc);
        md_array_ensure(sys->atom.type_idx, reserve_size, alloc);
        md_array_ensure(sys->atom.flags, reserve_size, alloc);

        md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, alloc);  // Ensure that index 0 is always unknown

        uint64_t prev_comp_key = 0; // Key of active componenent
        uint64_t prev_inst_key = 0; // Key of active instance

        for (size_t i = 0; i < num_atoms; ++i) {
            // Ignore alt loc entries unless 'A'
            if (atom_entries[i].label_alt_id != '.' &&
                atom_entries[i].label_alt_id != 'A') continue;

            str_t symbol = str_from_cstrn(atom_entries[i].type_symbol, sizeof(atom_entries[i].type_symbol));
            str_t atom_id = str_from_cstrn(atom_entries[i].label_atom_id, sizeof(atom_entries[i].label_atom_id));
            str_t entity_id = str_from_cstrn(atom_entries[i].label_entity_id, sizeof(atom_entries[i].label_entity_id));
            str_t auth_asym_id = str_from_cstrn(atom_entries[i].auth_asym_id, sizeof(atom_entries[i].auth_asym_id));
            str_t label_asym_id = str_from_cstrn(atom_entries[i].label_asym_id, sizeof(atom_entries[i].label_asym_id));

            md_entity_idx_t entity_idx = md_entity_find_by_id(&sys->entity, entity_id);

		    uint64_t inst_key = md_hash64_str(label_asym_id, 0);

            str_t comp_id = str_from_cstrn(atom_entries[i].label_comp_id, ARRAY_SIZE(atom_entries[i].label_comp_id));
            int seq_id = atom_entries[i].label_seq_id != INT32_MIN ? atom_entries[i].label_seq_id : atom_entries[i].auth_seq_id;
            uint64_t comp_key = md_hash64_str(comp_id, (uint64_t)seq_id ^ inst_key);

            md_atomic_number_t atomic_number = md_atomic_number_from_symbol(symbol, true);
            float mass = md_atomic_number_mass(atomic_number);
            float radius = md_atomic_number_vdw_radius(atomic_number);
            md_atom_type_idx_t atom_type_idx = md_atom_type_find_or_add(&sys->atom.type, atom_id, atomic_number, mass, radius, alloc);

            md_flags_t flags = md_entity_flags(&sys->entity, entity_idx);

            if (comp_key != prev_comp_key) {
                // --- New residue boundary ---
                static const uint32_t comp_flag_filter = MD_FLAG_HETERO | MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE | MD_FLAG_WATER | MD_FLAG_ISOMER_L | MD_FLAG_ISOMER_D;
                md_flags_t comp_flags = flags & comp_flag_filter;

                // Start residue (store atom offset before adding any residue atoms)
                md_array_push(sys->comp.atom_offset, (uint32_t)sys->atom.count, alloc);
                md_array_push(sys->comp.name,   make_label(comp_id), alloc);
                md_array_push(sys->comp.seq_id, seq_id, alloc);
                md_array_push(sys->comp.flags,  comp_flags, alloc);

                // Instance handling
                if (inst_key != prev_inst_key) {
                    // Start a new instance
                    md_label_t inst_id      = make_label(label_asym_id);
                    md_label_t inst_auth_id = make_label(auth_asym_id);

                    md_array_push(sys->inst.id, inst_id, alloc);
                    md_array_push(sys->inst.auth_id, inst_auth_id, alloc);
                    md_array_push(sys->inst.comp_offset, (uint32_t)sys->comp.count, alloc);
                    md_array_push(sys->inst.entity_idx, entity_idx, alloc);
                    sys->inst.count += 1;
                }
                sys->comp.count += 1;
            }
 
            prev_inst_key = inst_key;
            prev_comp_key = comp_key;
 
            sys->atom.count += 1;
            md_array_push_no_grow(sys->atom.x, atom_entries[i].x);
            md_array_push_no_grow(sys->atom.y, atom_entries[i].y);
            md_array_push_no_grow(sys->atom.z, atom_entries[i].z);
            md_array_push_no_grow(sys->atom.type_idx, atom_type_idx);
            md_array_push_no_grow(sys->atom.flags, flags);
        }

        if (sys->comp.atom_offset) {
            md_array_push(sys->comp.atom_offset, (uint32_t)sys->atom.count, alloc);  // Final sentinel
        }
        if (sys->inst.comp_offset) {
            md_array_push(sys->inst.comp_offset, (uint32_t)sys->comp.count, alloc);  // Final sentinel
        }
    }

    if (cell_parsed) {
        sys->unitcell = md_unitcell_from_extent_and_angles(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
    } else {
        sys->unitcell = md_unitcell_none();
    }

    return sys->atom.count > 0;
}

static bool mmcif_init_from_str(md_system_t* sys, str_t str, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);
    MEMSET(sys, 0, sizeof(md_system_t));
    return mmcif_parse(sys, &reader, alloc);
}

static bool mmcif_init_from_file(md_system_t* sys, str_t filename, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);

    if (file) {
        const size_t pos = md_temp_get_pos();
        const size_t cap = MEGABYTES(1);
        void* buf = md_temp_push(cap);

        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        MEMSET(sys, 0, sizeof(md_system_t));
        result = mmcif_parse(sys, &reader, alloc);

        md_temp_set_pos_back(pos);
        md_file_close(file);
    }
    return result;
}

static md_system_loader_i api = {
    .init_from_str = mmcif_init_from_str,
    .init_from_file = mmcif_init_from_file,
};

struct md_system_loader_i* md_mmcif_system_loader(void) {
    return &api;
}
