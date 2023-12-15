#include "md_mmcif.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_log.h>
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

static const int required_fields[] = {
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
	for (size_t i = 0; i < ARRAY_SIZE(required_fields); ++i) {
		if (table[required_fields[i]] == -1) {
			MD_LOG_ERROR("Missing required column in _atom_site: "STR_FMT, STR_ARG(atom_site_labels[required_fields[i]]));
			return false;
		}
	}

	int32_t prev_entity_id = -1;
	while (md_buffered_reader_peek_line(&line, reader)) {
		if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
			const int64_t num_tokens = extract_tokens(tok, ARRAY_SIZE(tok), &line);
			if (num_tokens < num_cols) {
				MD_LOG_ERROR("Too few tokens in line: "STR_FMT, STR_ARG(line));
				return false;
			}

			md_flags_t flags = 0;
			int32_t entity_id = -1;
			
			md_element_t element = md_util_element_lookup(tok[table[ATOM_SITE_TYPE_SYMBOL]]);
			md_array_push(atom->element, element, alloc);

			md_label_t type = make_label(tok[table[ATOM_SITE_LABEL_ATOM_ID]]);
			md_array_push(atom->type, type, alloc);
			
			md_label_t resname = make_label(tok[table[ATOM_SITE_LABEL_COMP_ID]]);
			md_array_push(atom->resname, resname, alloc);
			
			md_label_t chain_id = make_label(tok[table[ATOM_SITE_LABEL_ASYM_ID]]);
			md_array_push(atom->chainid, chain_id, alloc);

			str_t str = tok[table[ATOM_SITE_LABEL_ENTITY_ID]];
			if (str.len > 0 && str.ptr[0] != '.') {
				entity_id = (int32_t)parse_int(str);
			}
			
			{
				md_residue_id_t res_id = -INT32_MAX;
				str_t id = tok[table[ATOM_SITE_LABEL_SEQ_ID]];
				if (id.len > 0 && id.ptr[0] != '.') {
					res_id = (int32_t)parse_int(id);
				} else if (table[ATOM_SITE_AUTH_SEQ_ID] != -1) {
					id = tok[table[ATOM_SITE_AUTH_SEQ_ID]];
					if (id.len > 0 && id.ptr[0] != '.') {
						res_id = (int32_t)parse_int(id);
					}
				}
				md_array_push(atom->resid, res_id, alloc);
			}
			
			float x = (float)parse_float(tok[table[ATOM_SITE_CARTN_X]]);
			md_array_push(atom->x, x, alloc);
			
			float y = (float)parse_float(tok[table[ATOM_SITE_CARTN_Y]]);
			md_array_push(atom->y, y, alloc);
			
			float z = (float)parse_float(tok[table[ATOM_SITE_CARTN_Z]]);
			md_array_push(atom->z, z, alloc);
			

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

			md_array_push(atom->flags, flags, alloc);

			num_atoms += 1;
			md_buffered_reader_skip_line(reader);
		} else {
			break;
		}
	}

	atom->count = num_atoms;

	return num_atoms > 0;
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
	bool atom_site = false;

	str_t line;

	int64_t line_count = 0;
	while (md_buffered_reader_peek_line(&line, reader)) {
		if (line.len) {
			line = str_trim(line);

			if (str_eq_cstr_n(line, "_atom_site.", 11)) {
				if (!mmcif_parse_atom_site(&mol->atom, reader, alloc)) {
					MD_LOG_ERROR("Failed to parse _atom_site");
					return false;
				}
				atom_site = true;
			} else if (str_eq_cstr_n(line, "_cell.", 6)) {
				if (!mmcif_parse_cell(&mol->unit_cell, reader)) {
					MD_LOG_ERROR("Failed to parse _cell");
					return false;
				}
			}
		}

		line_count += 1;
		md_buffered_reader_skip_line(reader);
	}

	return atom_site;
}

static bool mmcif_init_from_str(md_molecule_t* mol, str_t str, md_allocator_i* alloc) {
	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	MEMSET(mol, 0, sizeof(md_molecule_t));
	return mmcif_parse(mol, &reader, alloc);
}

static bool mmcif_init_from_file(md_molecule_t* mol, str_t filename, md_allocator_i* alloc) {
	bool result = false;
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);

	if (file) {
		const int64_t cap = MEGABYTES(1);
		void* buf = md_alloc(md_heap_allocator, cap);

		md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
		MEMSET(mol, 0, sizeof(md_molecule_t));
		result = mmcif_parse(mol, &reader, alloc);

		md_free(md_heap_allocator, buf, cap);
		md_file_close(file);
	}
	return result;
}

static md_molecule_loader_i api = {
	.init_from_str = mmcif_init_from_str,
	.init_from_file = mmcif_init_from_file,
};

struct md_molecule_loader_i* md_mmcif_molecule_api() {
	return &api;
}