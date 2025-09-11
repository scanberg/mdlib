#include "md_mmcif.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_arena_allocator.h>
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

static const int atom_field_required[] = {
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

static bool mmcif_parse_atom_site(md_atom_data_t* atom, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    ASSERT(atom);
    ASSERT(reader);

    const size_t temp_pos = md_temp_get_pos();

    str_t line;
    str_t tok[32];

    int table[ATOM_SITE_COUNT];
    MEMSET(table, -1, sizeof(table));
    int num_cols = 0;
    int num_atoms = 0;

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

    // Assert that we have all the required fields
    size_t max_tok = table[ATOM_SITE_AUTH_SEQ_ID] > 0 ? (size_t)table[ATOM_SITE_AUTH_SEQ_ID] + 1 : 0;
    for (size_t i = 0; i < ARRAY_SIZE(atom_field_required); ++i) {
        if (table[atom_field_required[i]] == -1) {
            MD_LOG_ERROR("Missing required column in _atom_site: "STR_FMT, STR_ARG(atom_site_labels[atom_field_required[i]]));
            return false;
        }
        max_tok = MAX(max_tok, (size_t)table[atom_field_required[i]] + 1);
    }

    // Construct hashmap to accelerate common symbols for lookup
    const uint32_t num_elem = 119;
    md_hashmap32_t elem_map = { .allocator = md_get_temp_allocator() };
    md_hashmap_reserve(&elem_map, 256);
    const str_t* symbols = md_util_element_symbols();
    for (uint32_t i = 0; i < num_elem; ++i) {
        str_t sym = symbols[i];
        if (str_len(symbols[i]) > 1) {
            char lbl[4] = {0};
            str_copy_to_char_buf(lbl, sizeof(lbl), sym);
            convert_to_upper(lbl, sizeof(lbl));
            md_hashmap_add(&elem_map, md_hash64(sym.ptr, sym.len, 0), i);
            md_hashmap_add(&elem_map, md_hash64(lbl, 2, 0), i);

        } else {
            md_hashmap_add(&elem_map, md_hash64(sym.ptr, sym.len, 0), i);
        }
    }

    int32_t prev_entity_id = -1;
    while (md_buffered_reader_peek_line(&line, reader)) {
        if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
            const size_t num_tokens = extract_tokens(tok, max_tok, &line);
            if (num_tokens < max_tok) {
                MD_LOG_ERROR("Too few tokens in line: "STR_FMT, STR_ARG(line));
                goto done;
            }

            md_flags_t flags = 0;
            int32_t entity_id = -1;
            
            str_t sym = tok[table[ATOM_SITE_TYPE_SYMBOL]];
            md_element_t elem = 0;
            
            // Prefer _atom_site.type_symbol if not '.'
            if (sym.len > 0 && sym.ptr[0] != '.') {
                uint32_t* entry = md_hashmap_get(&elem_map, md_hash64(sym.ptr, sym.len, 0));
                elem = entry ? (md_element_t)*entry : 0;
            }

            md_label_t type = make_label(tok[table[ATOM_SITE_LABEL_ATOM_ID]]);

            str_t altloc = tok[table[ATOM_SITE_LABEL_ALT_ID]];
            if (altloc.len > 0 && altloc.ptr[0] != '.') {
                // Skip Altloc if its not A
                if (altloc.ptr[0] != 'A') {
                    goto next;
                }
            }

            md_label_t resname = make_label(tok[table[ATOM_SITE_LABEL_COMP_ID]]);
            md_label_t chain_id = make_label(tok[table[ATOM_SITE_LABEL_ASYM_ID]]);

            str_t str = tok[table[ATOM_SITE_LABEL_ENTITY_ID]];
            if (str.len > 0 && str.ptr[0] != '.') {
                entity_id = (int32_t)parse_int(str);
            }
            
            md_residue_id_t res_id = -INT32_MAX;
            str_t seq_id = tok[table[ATOM_SITE_LABEL_SEQ_ID]];
            if (seq_id.len > 0 && seq_id.ptr[0] != '.') {
                res_id = (int32_t)parse_int(seq_id);
            } else if (table[ATOM_SITE_AUTH_SEQ_ID] != -1) {
                seq_id = tok[table[ATOM_SITE_AUTH_SEQ_ID]];
                if (seq_id.len > 0 && seq_id.ptr[0] != '.') {
                    res_id = (int32_t)parse_int(seq_id);
                }
            }
            
            float x = (float)parse_float(tok[table[ATOM_SITE_CARTN_X]]);
            float y = (float)parse_float(tok[table[ATOM_SITE_CARTN_Y]]);
            float z = (float)parse_float(tok[table[ATOM_SITE_CARTN_Z]]);

            if (tok[0].ptr[0] == 'H') {
                flags |= MD_FLAG_HETATM;
            }
            if (entity_id != prev_entity_id) {
                flags |= MD_FLAG_CHAIN_BEG;
                if (atom->flags) {
                    *md_array_last(atom->flags) |= MD_FLAG_CHAIN_END;
                }
                prev_entity_id = entity_id;
            }

            md_array_push(atom->element, elem,     alloc);
            md_array_push(atom->type,	 type,     alloc);
            md_array_push(atom->x,		 x,	       alloc);
            md_array_push(atom->y,		 y,		   alloc);
            md_array_push(atom->z,		 z,        alloc);
            md_array_push(atom->flags,	 flags,	   alloc);
            md_array_push(atom->resid,	 res_id,   alloc);
            md_array_push(atom->resname, resname,  alloc);
            md_array_push(atom->chainid, chain_id, alloc);

            num_atoms += 1;
            next:
            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

    if (num_atoms > 0) {
        size_t capacity = ROUND_UP(num_atoms, 16);
        md_array_ensure(atom->element,  capacity, alloc);
        md_array_ensure(atom->type,     capacity, alloc);
        md_array_ensure(atom->x,        capacity, alloc);
        md_array_ensure(atom->y,        capacity, alloc);
        md_array_ensure(atom->z,        capacity, alloc);
        md_array_ensure(atom->flags,    capacity, alloc);
        md_array_ensure(atom->resid,    capacity, alloc);
        md_array_ensure(atom->resname,  capacity, alloc);
        md_array_ensure(atom->chainid,  capacity, alloc);
    }

    // Fill in missing elements using hash-backed inference
    for (size_t i = 0; i < (size_t)num_atoms; ++i) {
        if (atom->element[i] == 0) {
            // Use hash-backed inference for missing elements
            str_t atom_label = LBL_TO_STR(atom->type[i]);
            
            // Try direct element lookup first
            md_element_t elem = md_util_element_lookup_ignore_case(atom_label);
            
            // If that fails, use the hash-backed inference from md_util_element_guess
            if (elem == 0) {
                // Create a temporary molecule structure for the single atom to use element_guess
                md_molecule_t temp_mol = {0};
                temp_mol.atom.count = 1;
                temp_mol.atom.type = &atom->type[i];
                temp_mol.atom.resname = &atom->resname[i];
                temp_mol.atom.flags = atom->flags ? &atom->flags[i] : NULL;
                
                md_element_t temp_element = 0;
                if (md_util_element_guess(&temp_element, 1, &temp_mol)) {
                    elem = temp_element;
                }
            }
            
            atom->element[i] = elem;
        }
    }

    atom->count = num_atoms;
done:
    md_temp_set_pos_back(temp_pos);
    return atom->count > 0;
}

static bool mmcif_parse_cell(md_unit_cell_t* cell, md_buffered_reader_t* reader) {
    ASSERT(cell);
    ASSERT(reader);
    str_t line;
    uint32_t cell_flags = 0;
    double param[6] = {0};

    while (md_buffered_reader_peek_line(&line, reader)) {
        line = str_trim(line);
        if (str_eq_cstr_n(line, "_cell.", 6)) {
            str_t tok[2];
            const int64_t num_tokens = extract_tokens(tok, ARRAY_SIZE(tok), &line);
            if (num_tokens == 0) {
                break;
            }

            if (str_eq_cstr(tok[0], "_cell.angle_alpha")) {
                cell_flags |= 1;
                param[0] = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.angle_beta")) {
                cell_flags |= 2;
                param[1] = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.angle_gamma")) {
                cell_flags |= 4;
                param[2] = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.length_a")) {
                cell_flags |= 8;
                param[3] = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.length_b")) {
                cell_flags |= 16;
                param[4] = parse_float(tok[1]);
            } else if (str_eq_cstr(tok[0], "_cell.length_c")) {
                cell_flags |= 32;
                param[5] = parse_float(tok[1]);
            }

            md_buffered_reader_skip_line(reader);
        } else {
            break;
        }
    }

    if (cell_flags == 63) {
        if (param[0] == 90.0 && param[1] == 90.0 && param[2] == 90.0 && param[3] == 1.0 && param[4] == 1.0 && param[5] == 1.0) {
            // This is the identity matrix, and in such case, we assume there is no unit cell (no periodic boundary conditions)
        } else {
            *cell = md_util_unit_cell_from_extent_and_angles(param[3],param[4],param[5],param[0],param[1],param[2]);
        }
        return true;
    }
    return false;
}

static bool mmcif_parse(md_molecule_t* mol, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    bool atom_site_found = false;

    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        if (line.len) {
            line = str_trim(line);

            if (str_eq_cstr_n(line, "_atom_site.", 11)) {
                if (!mmcif_parse_atom_site(&mol->atom, reader, alloc)) {
                    MD_LOG_ERROR("Failed to parse _atom_site");
                    return false;
                }
                atom_site_found = true;
            } else if (str_eq_cstr_n(line, "_cell.", 6)) {
                if (!mmcif_parse_cell(&mol->unit_cell, reader)) {
                    MD_LOG_ERROR("Failed to parse _cell");
                    return false;
                }
            }
        }
        md_buffered_reader_skip_line(reader);
    }

    return atom_site_found;
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
