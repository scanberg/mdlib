#include "md_mmcif.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <md_util.h>
#include <md_molecule.h>

#define BAKE(str) {str, sizeof(str)-1}

// Enumerate fields in _atom_site
enum {
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
    int id;
    int label_entity_id;
    int label_seq_id;
//    int pdbx_PDB_model_num;

    char label_atom_id[5];
    char label_comp_id[5];
    char type_symbol[4];
    char label_alt_id;
    char label_asym_id;
    atom_site_group_pdb_t group_PDB;
//    char pdbx_PDB_ins_code;
//    int8_t pdbx_formal_charge;

    float x, y, z;
//    float occupancy;
//    float b_iso_or_equiv;
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

typedef uint8_t mmcif_entity_type_t;

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
    BAKE("polydeoxyribonucleotide"),
    BAKE("polydeoxyribonucleotide/polyribonucleotide hybrid"),
    BAKE("polypeptide(D)"),
    BAKE("polypeptide(L)"),
    BAKE("polyribonucleotide"),
};

typedef uint8_t mmcif_entity_poly_type_t;

typedef struct {
    int id;
    mmcif_entity_type_t type;
    mmcif_entity_poly_type_t poly_type;
    // perhaps add description field later
} mmcif_entity_t;

typedef struct {
    int id;
    int entity_id;
} mmcif_struct_asym_t;

static inline int parse_int_with_default(str_t tok, int def) {
    if (tok.len == 0 || tok.ptr[0] == '.' || tok.ptr[0] == '?') {
        return def;
    }
    return (int)parse_int(tok);
}

static inline str_t append_token(str_t cur, str_t new_tok) {
    const char* beg = str_beg(cur);
    const char* end = str_end(new_tok);
    return (str_t){beg, end - beg};
}

// Extracts tokens from entries that potentially contain multi-line strings
// Since it potentially contains multi-line strings, we need to allocate new memory for
// the tokens as the buffered reader may flush itself while extracting new lines
// Returns the number of extracted tokens
static size_t mmcif_extract_entry_tokens(str_t out_tokens[], size_t num_entry_fields, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    str_t line;
    size_t count = 0;

    md_strb_t sb = {.alloc = alloc};

    while (count < num_entry_fields && md_buffered_reader_extract_line(&line, reader)) {
        // End of entries
        if (line.len > 0 && line.ptr[0] == '#') {
            break;
        }

        // Multiline comment
        if (line.len > 0 && line.ptr[0] == ';') {
            md_strb_reset(&sb);
            md_strb_push_str(&sb, str_substr(line, 1, SIZE_MAX));

            while (md_buffered_reader_extract_line(&line, reader)) {
                if (line.len > 0 && line.ptr[0] == ';') break;
                md_strb_push_str(&sb, line);
            }
            out_tokens[count++] = str_copy(md_strb_to_str(sb), alloc);
        }

        // Standard tokenization of the line
        str_t tok;
        while (count < num_entry_fields && extract_token(&tok, &line)) {
            if (tok.ptr[0] == '\'' || tok.ptr[0] == '\"') {
                size_t loc;
                if (!str_find_char(&loc, line, tok.ptr[0])) {
                    MD_LOG_ERROR("Unpaired quotes in entry!");
                    goto done;
                }
                tok.len  += loc - 1;
                line.ptr += loc;
                line.len -= loc;
            }
            out_tokens[count++] = str_copy(tok, alloc);
        }
    }

done:
    return count;
}

enum {
    ENTITY_ID,
    ENTITY_TYPE,
    ENTITY_SRC_METHOD,
    ENTITY_COUNT,
};

static const str_t entity_labels[] = {
    BAKE("id"),
    BAKE("type"),
    BAKE("src_method"),
};

static bool mmcif_parse_entity(md_array(mmcif_entity_t)* entities, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    ASSERT(entities);
    ASSERT(reader);
    str_t line;
    str_t tok[16];
    int table[ENTITY_COUNT];
    MEMSET(table, -1, sizeof(table));
    int num_cols = 0;

    bool success = false;

    while (md_buffered_reader_peek_line(&line, reader)) {
        line = str_trim(line);
        if (str_eq_cstr_n(line, "_entity.", 8)) {
            str_t field = str_substr(line, 8, SIZE_MAX);
            for (int i = 0; i < ENTITY_COUNT; ++i) {
                if (table[i] != -1) {
                    continue;
                }
                if (str_eq(field, entity_labels[i])) {
                    table[i] = num_cols;
                    break;
                }
            }
            num_cols += 1;
            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

    md_allocator_i* temp_alloc = md_get_temp_allocator();
    size_t temp_pos = md_temp_get_pos();

    // Ensure that we have all the required fields
    for (size_t i = 0; i < ARRAY_SIZE(table); ++i) {
        if (table[i] == -1) {
            MD_LOG_ERROR("Missing required field in _entity: "STR_FMT, STR_ARG(entity_labels[i]));
            goto done;
        }
    }

    while (md_buffered_reader_peek_line(&line, reader)) {
        if (str_len(line) > 0 && line.ptr[0] == '#') {
            break;
        }

        size_t num_tokens = mmcif_extract_entry_tokens(tok, num_cols, reader, temp_alloc);
        if (num_tokens != num_cols) {
            MD_LOG_ERROR("Too few tokens in entry when parsing entity section");
            goto done;
        }

        int id = parse_int_with_default(tok[table[ENTITY_ID]], -1);
        int type = ENTITY_TYPE_UNKNOWN;
        for (int i = 0; i < (int)ARRAY_SIZE(entity_type_str); ++i) {
            if (str_eq(tok[table[ENTITY_TYPE]], entity_type_str[i])) {
                type = i;
            }
        }

        mmcif_entity_t entity = {
            .id = id,
            .type = type,
        };

        md_array_push(*entities, entity, alloc);
    }

    success = true;
done:
    md_temp_set_pos_back(temp_pos);
    return success;
}

enum {
    ENTITY_POLY_ENTITY_ID = 0,
    ENTITY_POLY_TYPE = 1,
    ENTITY_POLY_COUNT,
};

static const str_t entity_poly_labels[] = {
    BAKE("entity_id"),
    BAKE("type"),
};

// We only populate the entity_poly type within the given entities
static size_t mmcif_parse_entity_poly(mmcif_entity_t entities[], size_t num_entities, md_buffered_reader_t* reader) {
    ASSERT(reader);
    str_t line;
    str_t tok[16];

    int table[ENTITY_POLY_COUNT];
    MEMSET(table, -1, sizeof(table));
    int num_cols = 0;

    bool success = false;

    while (md_buffered_reader_peek_line(&line, reader)) {
        line = str_trim(line);
        if (str_eq_cstr_n(line, "_entity_poly.", 13)) {
            str_t field = str_substr(line, 13, SIZE_MAX);
            for (int i = 0; i < ENTITY_POLY_COUNT; ++i) {
                if (table[i] != -1) {
                    continue;
                }
                if (str_eq(field, entity_labels[i])) {
                    table[i] = num_cols;
                    break;
                }
            }
            num_cols += 1;
            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

    md_allocator_i* temp_alloc = md_get_temp_allocator();
    size_t temp_pos = md_temp_get_pos();

    // Ensure that we have all the required fields
    size_t max_tok = 0;
    for (size_t i = 0; i < ARRAY_SIZE(table); ++i) {
        if (table[i] == -1) {
            MD_LOG_ERROR("Missing required field in _entity_poly: "STR_FMT, STR_ARG(entity_labels[i]));
            goto done;
        }
        max_tok = MAX(max_tok, (size_t)table[i] + 1);
    }

    while (md_buffered_reader_peek_line(&line, reader)) {
        if (str_len(line) > 0 && line.ptr[0] == '#') {
            break;
        }

        size_t num_tokens = mmcif_extract_entry_tokens(tok, num_cols, reader, temp_alloc);
        if (num_tokens != num_cols) {
            MD_LOG_ERROR("Too few tokens in entry when parsing entity section");
            goto done;
        }

        int id = parse_int_with_default(tok[table[ENTITY_POLY_ENTITY_ID]], -1);
        int poly_type = ENTITY_POLY_TYPE_UNKNOWN;
        for (int i = 0; i < (int)ARRAY_SIZE(entity_poly_type_str); ++i) {
            if (str_eq(tok[table[ENTITY_POLY_TYPE]], entity_poly_type_str[i])) {
                poly_type = i;
                break;
            }
        }

        if (id == -1) {
            MD_LOG_ERROR("Invalid entity_id in _entity_poly entry");
            goto done;
        }

        if (poly_type == ENTITY_POLY_TYPE_UNKNOWN) {
            MD_LOG_ERROR("Invalid type in _entity_poly entry");
        }

        // Attempt to assign type
        for (size_t i = 0; i < num_entities; ++i) {
            if (entities[i].id = id) {
                entities[i].poly_type = (mmcif_entity_poly_type_t)poly_type;
            }
        }
    }

    success = true;
done:
    md_temp_set_pos_back(temp_pos);
    return success;
}

static bool mmcif_parse_atom_site(md_array(mmcif_atom_site_entry_t)* atom_entries, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    ASSERT(atom_entries);
    ASSERT(reader);

    str_t line;
    str_t tok[32];

    int table[ATOM_SITE_COUNT];
    MEMSET(table, -1, sizeof(table));
    int num_cols = 0;

    bool success = false;

    while (md_buffered_reader_peek_line(&line, reader)) {
        line = str_trim(line);
        if (str_eq_cstr_n(line, "_atom_site.", 11)) {
            str_t field = str_substr(line, 11, SIZE_MAX);
            for (int i = 0; i < (int)ARRAY_SIZE(atom_site_labels); ++i) {
                if (table[i] != -1) {
                    continue;
                }
                if (str_eq(field, atom_site_labels[i])) {
                    table[i] = num_cols;
                    break;
                }
            }
            num_cols += 1;
            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

    // Ensure that we have all the required fields
    size_t max_tok = table[ATOM_SITE_AUTH_SEQ_ID] > 0 ? (size_t)table[ATOM_SITE_AUTH_SEQ_ID] + 1 : 0;
    for (size_t i = 0; i < ARRAY_SIZE(required_atom_site_fields); ++i) {
        if (table[required_atom_site_fields[i]] == -1) {
            MD_LOG_ERROR("Missing required column in _atom_site: "STR_FMT, STR_ARG(atom_site_labels[required_atom_site_fields[i]]));
            goto done;
        }
        max_tok = MAX(max_tok, (size_t)table[required_atom_site_fields[i]] + 1);
    }

    while (md_buffered_reader_peek_line(&line, reader)) {
        if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
            const size_t num_tokens = extract_tokens(tok, max_tok, &line);
            if (num_tokens < max_tok) {
                MD_LOG_ERROR("Too few tokens in line: "STR_FMT, STR_ARG(line));
                goto done;
            }
            
            str_t alt_id    = tok[table[ATOM_SITE_LABEL_ALT_ID]];
            str_t asym_id   = tok[table[ATOM_SITE_LABEL_ASYM_ID]];

            const int default_value = -INT32_MAX;

            mmcif_atom_site_entry_t entry = {
                .id              = parse_int_with_default(tok[table[ATOM_SITE_LABEL_ATOM_ID]],   default_value),
                .label_entity_id = parse_int_with_default(tok[table[ATOM_SITE_LABEL_ENTITY_ID]], default_value),
                .label_seq_id    = parse_int_with_default(tok[table[ATOM_SITE_LABEL_SEQ_ID]],    default_value),

                .label_alt_id    =  alt_id.len == 1 ?  alt_id.ptr[0] : ' ',
                .label_asym_id   = asym_id.len == 1 ? asym_id.ptr[0] : ' ',
                .group_PDB       = line.ptr[0] == 'H' ? ATOM_SITE_GROUP_PDB_HETATM : ATOM_SITE_GROUP_PDB_ATOM,
                
                .x = (float)parse_float(tok[table[ATOM_SITE_CARTN_X]]),
                .y = (float)parse_float(tok[table[ATOM_SITE_CARTN_Y]]),
                .z = (float)parse_float(tok[table[ATOM_SITE_CARTN_Z]]),
            };

            str_copy_to_char_buf(entry.type_symbol,   sizeof(entry.type_symbol),   tok[table[ATOM_SITE_TYPE_SYMBOL]]);
            str_copy_to_char_buf(entry.label_atom_id, sizeof(entry.label_atom_id), tok[table[ATOM_SITE_LABEL_ATOM_ID]]);
            str_copy_to_char_buf(entry.label_comp_id, sizeof(entry.label_comp_id), tok[table[ATOM_SITE_LABEL_COMP_ID]]);

            md_array_push(*atom_entries, entry, alloc);

            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

    success = true;
done:
    return success;
}

<<<<<<< HEAD
static bool mmcif_parse_cell(mmcif_cell_params_t* cell, md_buffered_reader_t* reader) {
=======
static bool mmcif_parse_cell(md_unitcell_t* cell, md_buffered_reader_t* reader) {
>>>>>>> 10c9697 (unit_cell -> unitcell and moved unitcell functionality into types.h from util)
    ASSERT(cell);
    ASSERT(reader);
    str_t line;

    while (md_buffered_reader_peek_line(&line, reader)) {
        line = str_trim(line);
        if (str_eq_cstr_n(line, "_cell.", 6)) {
            str_t tok[2];
            const size_t num_tokens = extract_tokens(tok, ARRAY_SIZE(tok), &line);
            if (num_tokens == 0) {
                break;
            }
            if (str_eq(tok[0], STR_LIT("#"))) {
                break;
            }

            if (str_eq_cstr(tok[0], "_cell.angle_alpha")) {
                cell->alpha = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.angle_beta")) {
                cell->beta = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.angle_gamma")) {
                cell->gamma = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.length_a")) {
                cell->a = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.length_b")) {
                cell->b = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.length_c")) {
                cell->c = parse_float(tok[1]);
            }

            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

<<<<<<< HEAD
=======
    if (cell_flags == 63) {
        if (param[0] == 90.0 && param[1] == 90.0 && param[2] == 90.0 && param[3] == 1.0 && param[4] == 1.0 && param[5] == 1.0) {
            // This is the identity matrix, and in such case, we assume there is no unit cell (no periodic boundary conditions)
        } else {
            *cell = md_unitcell_from_extent_and_angles(param[3], param[4], param[5], param[0], param[1], param[2]);
        }
        return true;
    }
>>>>>>> 10c9697 (unit_cell -> unitcell and moved unitcell functionality into types.h from util)
    return false;
}

static bool mmcif_parse(md_molecule_t* mol, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    bool atom_site_parsed = false;
    bool cell_parsed = false;
    bool entity_parsed = false;
    bool entity_poly_parsed = false;

    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(4));

    md_array(mmcif_atom_site_entry_t) atom_entries = 0;
    md_array(mmcif_entity_t) entities = 0;

    md_array_ensure(atom_entries, 1024, temp_arena);
    md_array_ensure(entities, 4, temp_arena);

    mmcif_cell_params_t cell = {0};

    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        if (line.len > 0) {
            line = str_trim(line);
            if (str_eq_cstr_n(line, "_atom_site.", 11)) {
                if (!mmcif_parse_atom_site(&atom_entries, reader, temp_arena)) {
                    MD_LOG_ERROR("Failed to parse _atom_site section");
                    return false;
                }
                atom_site_parsed = true;
            } else if (str_eq_cstr_n(line, "_cell.", 6)) {
<<<<<<< HEAD
                if (!mmcif_parse_cell(&cell, reader)) {
                    MD_LOG_ERROR("Failed to parse _cell section");
=======
                if (!mmcif_parse_cell(&mol->unitcell, reader)) {
                    MD_LOG_ERROR("Failed to parse _cell");
>>>>>>> 10c9697 (unit_cell -> unitcell and moved unitcell functionality into types.h from util)
                    return false;
                }
                cell_parsed = true;
            } else if (str_eq_cstr_n(line, "_entity.", 8)) {
                if (!mmcif_parse_entity(&entities, reader, temp_arena)) {
                    MD_LOG_ERROR("Failed to parse _entity section");
                    return false;
                }
                entity_parsed = true;
            } else if (str_eq_cstr_n(line, "_entity_poly.", 13)) {
                if (!mmcif_parse_entity_poly(entities, md_array_size(entities), reader)) {
                    MD_LOG_ERROR("Failed to parse _entity_poly section");
                    return false;
                }
                entity_poly_parsed = true;
            }

            if (atom_site_parsed && cell_parsed && entity_parsed && entity_poly_parsed) {
                break;
            }
        }
        md_buffered_reader_skip_line(reader);
    }

    // Populate molecule from parsed data
    if (atom_site_parsed) {
        size_t num_atoms = md_array_size(atom_entries);
        size_t reserve_size = ALIGN_TO(num_atoms, 16);

        md_array_ensure(mol->atom.x, reserve_size, alloc);
        md_array_ensure(mol->atom.y, reserve_size, alloc);
        md_array_ensure(mol->atom.z, reserve_size, alloc);
        md_array_ensure(mol->atom.type_idx, reserve_size, alloc);
        md_array_ensure(mol->atom.flags, reserve_size, alloc);

        md_array(str_t) atom_resname = 0;
        md_array(int)   atom_resid = 0;
        md_array_ensure(atom_resname, reserve_size, temp_arena);
        md_array_ensure(atom_resid,   reserve_size, temp_arena);

        mol->atom.type_data.count = 0;
        md_atom_type_find_or_add(&mol->atom.type_data, STR_LIT("Unknown"), 0, 0.0f, 0.0f, alloc);  // Ensure that index 0 is always unknown

        for (size_t i = 0; i < num_atoms; ++i) {
            if (atom_entries[i].label_alt_id != ' ') continue;

            str_t symbol = str_from_cstrn(atom_entries[i].type_symbol, ARRAY_SIZE(atom_entries[i].type_symbol));
            str_t atom_id = str_from_cstrn(atom_entries[i].label_atom_id, ARRAY_SIZE(atom_entries[i].label_atom_id));
            md_atomic_number_t atomic_number = md_atomic_number_from_symbol(symbol);
            float mass = md_atomic_mass(atomic_number);
            float radius = md_vdw_radius(atomic_number);
            md_atom_type_idx_t atom_type_idx = md_atom_type_find_or_add(&mol->atom.type_data, atom_id, atomic_number, mass, radius, alloc);

            md_flags_t flags = 0;
            flags |= atom_entries[i].group_PDB == ATOM_SITE_GROUP_PDB_HETATM ? MD_FLAG_HETATM : 0;
            

            mol->atom.count += 1;
            md_array_push_no_grow(mol->atom.x, atom_entries[i].x);
            md_array_push_no_grow(mol->atom.y, atom_entries[i].y);
            md_array_push_no_grow(mol->atom.z, atom_entries[i].z);
            md_array_push_no_grow(mol->atom.type_idx, atom_type_idx);
            md_array_push_no_grow(mol->atom.flags, 0);

            md_array_push_no_grow(atom_resname, str_from_cstrn(atom_entries[i].label_comp_id, ARRAY_SIZE(atom_entries[i].label_comp_id)));
            md_array_push_no_grow(atom_resid, atom_entries[i].label_seq_id);
        }
        md_util_init_residue_data(&mol->residue, mol->atom.flags, atom_resid, atom_resname, mol->atom.count, alloc);
    }

    if (entity_parsed) {
        // Add protein chains
        // Set flags for 
    }

    if (cell_parsed) {
        mol->unit_cell = md_util_unit_cell_from_extent_and_angles(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
    }

    return mol->atom.count > 0;
}

static bool mmcif_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);
    MEMSET(mol, 0, sizeof(md_molecule_t));
    return mmcif_parse(mol, &reader, alloc);
}

static bool mmcif_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);

    if (file) {
        const size_t pos = md_temp_get_pos();
        const size_t cap = MEGABYTES(1);
        void* buf = md_temp_push(cap);

        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        MEMSET(mol, 0, sizeof(md_molecule_t));
        result = mmcif_parse(mol, &reader, alloc);

        md_temp_set_pos_back(pos);
        md_file_close(file);
    }
    return result;
}

static md_molecule_loader_i api = {
    .init_from_str = mmcif_init_from_str,
    .init_from_file = mmcif_init_from_file,
};

struct md_molecule_loader_i* md_mmcif_molecule_api(void) {
    return &api;
}
