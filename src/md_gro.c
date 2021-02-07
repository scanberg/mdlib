#include "md_gro.h"

#include <core/common.h>
#include <core/str_util.h>
#include <md_allocator.h>
#include <md_log.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ATOM_GROW_CAPACITY 32
#define RES_GROW_CAPACITY  8
internal inline int32_t compute_position_field_width(str_t line) {
    // First float starts at offset 20, count
    if (line.len > 20) {
        const char* c = line.str + 20;
        const char* end = line.str + line.len;
        while (c != end && *c != '\n' && *c == ' ') ++c;
        while (c != end && *c != '\n' && *c != ' ') ++c;
        if (c != end) {
            return (int32_t)(c - (line.str + 20));
        }
    }
    return 0;
}

// GRO files are have a very simple structure and the atom entry is very similar to XYZ or pdb, but with fewer fields.
// First line is a header
// Second line is the number of atoms present
// Following N lines contain the atom data, where N is the number of atoms
// Last line holds the simulation box dimensions
/*
md_error md_parse_gro(const md_parse_gro_desc* desc) {
    ASSERT(desc);

    md_molecule* mol = desc->molecule;
    struct md_log_i* err_log = NULL;
    struct md_allocator_i* alloc = NULL;

    if (desc->optional.context) {
        err_log = &desc->optional.context->error_log;
    }

    if (!desc->str) {
        log_error(err_log, MD_LOG_SEVERITY_ERROR, "Input string was NULL");
        return MD_PARSE_ERROR_NULL_PARAMETER;
    }
    if (!desc->molecule) {
        log_error(err_log, MD_LOG_SEVERITY_ERROR, "Molecule was NULL");
        return MD_PARSE_ERROR_NULL_PARAMETER;
    }
    memset(mol, 0, sizeof(*mol));

    // Use str as our string stream and modify it as we parse
    str_t str = {desc->str, desc->str + desc->str_len};
    str_t line;
    str_t header;
    str_t atom_count_str;
    int pos_width = 0;

    if (!extract_line(&header, &str)) {           // header
        log_error(err_log, MD_LOG_SEVERITY_ERROR, "Could not read header.");
        return MD_PARSE_ERROR_MALFORMED_INPUT;
    }

    if (!extract_line(&atom_count_str, &str)) {    // num atoms
        log_error(err_log, MD_LOG_SEVERITY_ERROR, "Could not read number of atoms.");
        return MD_PARSE_ERROR_MALFORMED_INPUT;
    }

    uint32_t atom_count = parse_int(atom_count_str);
    if (atom_count == 0) {
        log_error(err_log, MD_LOG_SEVERITY_ERROR, "Could not parse number of atoms");
        return false;
    }

    if (peek_line(&line, &str)) {
        pos_width = compute_position_field_width(line);
    }

    if (pos_width == 0) {
        log_error(err_log, MD_LOG_SEVERITY_ERROR, "Could not identify internal line format of gro file");
        return MD_PARSE_ERROR_MALFORMED_INPUT;
    }

    const md_atom_field atom_fields = MD_ATOM_FIELD_POS | MD_ATOM_FIELD_NAME | MD_ATOM_FIELD_RESIDUE_IDX;
    realloc_atom_fields(mol, atom_fields, atom_count, alloc);
    
    uint32_t res_capacity = RES_GROW_CAPACITY;
    const md_residue_field residue_fields = MD_RESIDUE_FIELD_ATOM_RANGE | MD_RESIDUE_FIELD_NAME | MD_RESIDUE_FIELD_ID;
    realloc_residue_fields(mol, residue_fields, res_capacity, alloc);
    
    uint32_t res_count = 0;
    uint32_t cur_res_id = -1;
    for (int i = 0; i < atom_count; ++i) {
        if (!extract_line(&line, &str)) {
            log_error(err_log, MD_LOG_SEVERITY_ERROR, "Failed to extract line, too few lines for the number of atoms specified");
            return MD_PARSE_ERROR_MALFORMED_INPUT;
        }

        int res_id = parse_int(str_sub(line, 0, 5));
        str_t name = trim_whitespace(str_sub(line, 10, 5));

        float x, y, z;
        mol->atom.x[i] = parse_float(trim_whitespace(str_sub(line, 20 + 0 * pos_width, pos_width))) * 10.0f; // nm -> Å
        mol->atom.y[i] = parse_float(trim_whitespace(str_sub(line, 20 + 1 * pos_width, pos_width))) * 10.0f;
        mol->atom.z[i] = parse_float(trim_whitespace(str_sub(line, 20 + 2 * pos_width, pos_width))) * 10.0f;
        mol->atom.residue_idx[i] = res_count;

        mol->residue.atom_range[res_count].end += 1;
        
        if (res_id != cur_res_id) { 
            cur_res_id = res_id;
            str_t resname = trim_whitespace(str_sub(line, 5, 5));
            mol->residue.atom_range[res_count].beg = i;
            mol->residue.atom_range[res_count].end = i;
            ++res_count;
            if (res_count == res_capacity) {
                res_capacity += RES_GROW_CAPACITY;
                realloc_residue_fields(mol, residue_fields, res_capacity, alloc);
            }
        }
    }

    // Simulation box extent
    float bx, by, bz;
    bx = parse_float(trim_whitespace(str_sub(line, 0, 8)))  * 10.0f;
    by = parse_float(trim_whitespace(str_sub(line, 10, 8))) * 10.0f;
    bz = parse_float(trim_whitespace(str_sub(line, 20, 8))) * 10.0f;

    mol->atom.count = atom_count;
    mol->residue.count = res_count;

    return MD_PARSE_SUCCESS;
}
*/


/*
md_error md_parse_pdb(const md_parse_pdb_desc* desc) {

    return MD_PARSE_SUCCESS;
}
*/

#ifdef __cplusplus
}
#endif