#include "md_gro.h"

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_array.inl>
#include <md_util.h>

#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline bool ranges_overlap(md_range a, md_range b) {
    return (a.beg < b.end && b.beg < a.end);
}

static inline int64_t compute_position_field_width(str_t line) {
    // First float starts at offset 20, count
    if (line.len > 20) {
        const char* c = line.ptr + 20;
        const char* end = line.ptr + line.len;
        while (c != end && *c != '\n' && *c == ' ') ++c;
        while (c != end && *c != '\n' && *c != ' ') ++c;
        if (c != end) {
            return (int64_t)(c - (line.ptr + 20));
        }
    }
    return 0;
}

static inline bool parse_header(str_t* str, struct md_gro_data_t* data) {
    str_t line = {0};
    if (!extract_line(&line, str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read title");
        return false;
    }

    line = trim_whitespace(line);
    const int64_t len = MIN(ARRAY_SIZE(data->title) - 1, line.len);
    strncpy(data->title, line.ptr, len);
    data->title[len] = '\0';

    if (!extract_line(&line, str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read number of atoms");
        return false;
    }

    data->num_atoms = parse_int(trim_whitespace(line));
    const int64_t min_bounds = 1;
    const int64_t max_bounds = 10000000;
    if (data->num_atoms < min_bounds || max_bounds < data->num_atoms) {
        md_printf(MD_LOG_TYPE_ERROR, "Number of atoms (%lli) outside of plausible range [%lli, %lli]", data->num_atoms, min_bounds, max_bounds);
        return false;
    }

    return true;
}

static inline void parse_atom_data(str_t* str, struct md_gro_data_t* data, int64_t pos_field_width, int64_t* read_count, int64_t read_target, md_allocator_i* alloc) {
    str_t line = {0};
    while (*read_count < read_target && extract_line(&line, str)) {
        int64_t res_id  = parse_int(trim_whitespace(substr(*str, 0, 5)));
        str_t res_name  = trim_whitespace(substr(*str, 5, 5));
        str_t atom_name = trim_whitespace(substr(*str, 10, 5));
        double x = parse_float(substr(*str, 20 + 0 * pos_field_width, pos_field_width)) * 10.0; // nm -> Å
        double y = parse_float(substr(*str, 20 + 1 * pos_field_width, pos_field_width)) * 10.0;
        double z = parse_float(substr(*str, 20 + 2 * pos_field_width, pos_field_width)) * 10.0;

        md_gro_atom_t atom = {0};
        atom.res_id = (int32_t)res_id;
        strncpy(atom.res_name, res_name.ptr, MIN(res_name.len, ARRAY_SIZE(atom.res_name) - 1));
        strncpy(atom.atom_name, atom_name.ptr, MIN(atom_name.len, ARRAY_SIZE(atom.atom_name) - 1));
        atom.x = (float)x;
        atom.y = (float)y;
        atom.z = (float)z;
        md_array_push(data->atom_data, atom, alloc);
        ++read_count;
    }
}

static inline bool parse_unitcell(str_t* str, struct md_gro_data_t* data) {
    ASSERT(str);
    ASSERT(data);

    str_t line = {0};
    if (!extract_line(&line, str)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to extract line for unit cell.");
        return false;
    }

    data->unit_cell[0] = (float)parse_float(trim_whitespace(substr(*str, 0, 10)));
    data->unit_cell[1] = (float)parse_float(trim_whitespace(substr(*str, 10, 10)));
    data->unit_cell[2] = (float)parse_float(trim_whitespace(substr(*str, 20, 10)));

    return true;
}

bool md_gro_data_parse_str(str_t str, struct md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    
    if (!parse_header(&str, data)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to parse gro header!");
        return false;
    }

    const int64_t pos_field_width = compute_position_field_width(str);

    int64_t num_atoms_read = 0;
    parse_atom_data(&str, data, pos_field_width, &num_atoms_read, data->num_atoms, alloc);

    if (num_atoms_read != data->num_atoms) {
        md_printf(MD_LOG_TYPE_ERROR, "There was a descrepancy between the number of atoms read (%lli) compared to the number of atoms specified in the header (%lli)", num_atoms_read, data->num_atoms);
        return false;
    }

    if (!parse_unitcell(&str, data)) {
        return false;
    }

    return true;
}

bool md_gro_data_parse_file(str_t filename, struct md_gro_data_t* data, struct md_allocator_i* alloc) {
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        const uint64_t buf_size = KILOBYTES(64ULL);
        char* buf_ptr = md_alloc(default_temp_allocator, buf_size);
        
        int64_t pos_field_width = 0;
        int64_t num_atoms_read = 0;

        uint64_t file_offset = 0;
        uint64_t bytes_read = 0;
        uint64_t read_offset = 0;
        uint64_t read_size = buf_size;

        // Read header explicitly
        {
            const int64_t header_read_size = KILOBYTES(1);
            bytes_read = md_file_read(file, buf_ptr + read_offset, header_read_size);
            str_t str = {.ptr = buf_ptr, .len = bytes_read};

            if (!parse_header(&str, data)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to parse gro header!");
                md_file_close(file);
                return false;
            }

            pos_field_width = compute_position_field_width(str);
            if (!pos_field_width) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to determine position field width");
                md_file_close(file);
                return false;
            }

            // Once we have read the header, we set the file position to just after the header
            const int64_t seek_target = str.ptr - buf_ptr;
            if (!md_file_seek(file, seek_target, MD_FILE_BEG)) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to seek in file");
                md_file_close(file);
                return false;
            }
        }

        // Read fields in chunks
        while ((bytes_read = md_file_read(file, buf_ptr + read_offset, read_size)) == read_size) {
            str_t str = {.ptr = buf_ptr, .len = bytes_read};

            // We want to make sure we send complete lines for parsing so we locate the last '\n'
            const int64_t last_new_line = rfind_char(str, '\n');
            ASSERT(str.ptr[last_new_line] == '\n');
            if (last_new_line == -1) {
                md_print(MD_LOG_TYPE_ERROR, "Unexpected structure within gro file, missing new lines?");
                md_file_close(file);
                return false;
            }
            str.len = last_new_line + 1;

            parse_atom_data(&str, data, pos_field_width, &num_atoms_read, data->num_atoms, alloc);

            // Copy remainder to beginning of buffer
            read_offset = buf_size - str.len;
            memcpy(buf_ptr, str.ptr + str.len, read_offset);
            read_size = buf_size - read_offset;
            file_offset += read_size;
        }

        // Handle remainder
        str_t str = {.ptr = buf_ptr, .len = bytes_read + read_offset}; // add read_offset here since we also have some portion which was copied
        parse_atom_data(&str, data, pos_field_width, &num_atoms_read, data->num_atoms, alloc);

        if (num_atoms_read != data->num_atoms) {
            md_printf(MD_LOG_TYPE_ERROR, "There was a descrepancy between the number of atoms read (%lli) compared to the number of atoms specified in the header (%lli)", num_atoms_read, data->num_atoms);
            md_file_close(file);
            return false;
        }

        if (!parse_unitcell(&str, data)) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to parse unitcell in gro file");
            md_file_close(file);
            return false;
        }

        md_file_close(file);
        return true;
    }
    md_printf(MD_LOG_TYPE_ERROR, "Could not open file '%.*s'", filename.len, filename.ptr);
    return false;
}

void md_gro_data_free(struct md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(data);
    if (data->atom_data) md_array_free(data->atom_data, alloc);
    memset(data, 0, sizeof(struct md_gro_data_t));
}

static inline str_t add_label(str_t str, md_gro_molecule_t* mol) {
    for (int64_t i = 0; i < md_array_size(mol->storage.labels); ++i) {
        if (strncmp(str.ptr, mol->storage.labels[i].str, str.len) == 0) {
            return (str_t){mol->storage.labels[i].str, str.len};
        }
    }
    md_gro_label_t label = {0};
    memcpy(label.str, str.ptr, MIN(str.len, ARRAY_SIZE(label.str) - 1));
    md_array_push(mol->storage.labels, label, mol->storage.arena);
    return (str_t){label.str, str.len};
}

bool md_gro_molecule_init(struct md_gro_molecule_t* gro, const struct md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(gro);
    ASSERT(data);
    ASSERT(alloc);

    memset(gro, 0, sizeof(md_gro_molecule_t));

    const int64_t num_atoms = data->num_atoms;

    md_array_ensure(gro->mol.atom.x, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.y, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.z, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.element, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.name, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.radius, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.mass, num_atoms, alloc);
    md_array_ensure(gro->mol.atom.flags, num_atoms, alloc);

    gro->storage.arena = md_arena_allocator_create(alloc, 4096);

    int32_t cur_res_id = -1;
    char cur_chain_id = -1;

    for (int64_t i = 0; i < num_atoms; ++i) {
        const float x = data->atom_data[i].x;
        const float y = data->atom_data[i].y;
        const float z = data->atom_data[i].z;
        str_t atom_name = (str_t){data->atom_data[i].atom_name, strnlen(data->atom_data[i].atom_name, ARRAY_SIZE(data->atom_data[i].atom_name))};
        atom_name = add_label(atom_name, gro);

        str_t res_name = (str_t){data->atom_data[i].res_name, strnlen(data->atom_data[i].res_name, ARRAY_SIZE(data->atom_data[i].res_name))};

        md_element element = md_util_decode_element(atom_name, res_name);
        if (element == 0) {
            md_printf(MD_LOG_TYPE_INFO, "Failed to decode element with atom name '%s' and residue name '%s'", data->atom_data[i].atom_name, data->atom_data[i].res_name);
        }

        const float radius = md_util_vdw_radius(element);
        const float mass = md_util_atomic_mass(element);

        int32_t res_id = data->atom_data[i].res_id;
        if (res_id != cur_res_id) {
            cur_res_id = res_id;

            res_name = add_label(res_name, gro);
            md_residue_id id = res_id;
            md_range atom_range = {gro->mol.atom.count, gro->mol.atom.count};

            gro->mol.residue.count += 1;
            md_array_push(gro->mol.residue.name, res_name.ptr, alloc);
            md_array_push(gro->mol.residue.id, id, alloc);
            md_array_push(gro->mol.residue.atom_range, atom_range, alloc);
        }

        if (gro->mol.residue.atom_range) md_array_last(gro->mol.residue.atom_range)->end += 1;

        gro->mol.atom.count += 1;
        md_array_push(gro->mol.atom.x, x, alloc);
        md_array_push(gro->mol.atom.y, y, alloc);
        md_array_push(gro->mol.atom.z, z, alloc);
        md_array_push(gro->mol.atom.element, element, alloc);
        md_array_push(gro->mol.atom.name, atom_name.ptr, alloc);
        md_array_push(gro->mol.atom.radius, radius, alloc);
        md_array_push(gro->mol.atom.mass, mass, alloc);
        md_array_push(gro->mol.atom.flags, 0, alloc);

        if (gro->mol.residue.count) md_array_push(gro->mol.atom.residue_idx, (md_residue_idx)(gro->mol.residue.count - 1), alloc);
    }

    {
        // Compute backbone
        md_array_ensure(gro->mol.residue.backbone_atoms, gro->mol.residue.count, alloc);
        md_util_backbone_args_t args = {
            .atom = {
                .count = gro->mol.atom.count,
                .name = gro->mol.atom.name,
        },
        .residue = {
                .count = gro->mol.residue.count,
                .atom_range = gro->mol.residue.atom_range,
        },
        };
        md_util_extract_backbone(gro->mol.residue.backbone_atoms, &args);
    }

    {
        // Compute secondary structures
        if (gro->mol.chain.count) {
            md_array_ensure(gro->mol.residue.secondary_structure, gro->mol.residue.count, alloc);
            md_util_secondary_structure_args_t args = {
                .atom = {
                    .count = gro->mol.atom.count,
                    .x = gro->mol.atom.x,
                    .y = gro->mol.atom.y,
                    .z = gro->mol.atom.z,
            },
            .residue = {
                    .count = gro->mol.residue.count,
                    .backbone_atoms = gro->mol.residue.backbone_atoms,
            },
            .chain = {
                    .count = gro->mol.chain.count,
                    .residue_range = gro->mol.chain.residue_range
            }
            };
            md_util_compute_secondary_structure(gro->mol.residue.secondary_structure, &args);
        }
    }

    {
        md_array_ensure(gro->mol.residue.internal_covalent_bond_range, gro->mol.residue.count, alloc);
        md_array_ensure(gro->mol.residue.complete_covalent_bond_range, gro->mol.residue.count, alloc);

        // Use heuristical method of finding covalent bonds
        md_util_covalent_bond_args_t args = {
            .atom = {
                .count = gro->mol.atom.count,
                .x = gro->mol.atom.x,
                .y = gro->mol.atom.y,
                .z = gro->mol.atom.z,
                .element = gro->mol.atom.element
        },
            .residue = {
                .count = gro->mol.residue.count,
                .atom_range = gro->mol.residue.atom_range,
                .internal_bond_range = gro->mol.residue.internal_covalent_bond_range,
                .complete_bond_range = gro->mol.residue.complete_covalent_bond_range,
        },
        };
        gro->mol.covalent_bond.bond = md_util_extract_covalent_bonds(&args, alloc);
        gro->mol.covalent_bond.count = md_array_size(gro->mol.covalent_bond.bond);
    }

    {
        // Compute artificial chains and label them A - Z
        if (gro->mol.residue.complete_covalent_bond_range) {
            const char* chain_labels[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
            int64_t chain_idx = 0;
            int64_t beg_res_idx = 0;
            for (int64_t res_idx = 0; res_idx < gro->mol.residue.count - 1; ++res_idx) {
                if (!ranges_overlap(gro->mol.residue.complete_covalent_bond_range[res_idx],
                                    gro->mol.residue.complete_covalent_bond_range[res_idx+1])) {
                    const int64_t end_res_idx = res_idx;
                    const str_t chain_id = {.ptr = chain_labels[chain_idx], .len = 1};
                    const md_range res_range  = {beg_res_idx, end_res_idx};
                    const md_range atom_range = {gro->mol.residue.atom_range[beg_res_idx].beg, gro->mol.residue.atom_range[end_res_idx].end};

                    md_array_push(gro->mol.chain.id, add_label(chain_id, gro).ptr, alloc);
                    md_array_push(gro->mol.chain.residue_range, res_range, alloc);
                    md_array_push(gro->mol.chain.atom_range, atom_range, alloc);
                    
                    chain_idx = (chain_idx + 1) % ARRAY_SIZE(chain_labels);
                }
            }

            const int64_t end_res_idx = gro->mol.residue.count;
            const str_t chain_id = {.ptr = chain_labels[chain_idx], .len = 1};
            const md_range res_range  = {beg_res_idx, end_res_idx};
            const md_range atom_range = {gro->mol.residue.atom_range[beg_res_idx].beg, gro->mol.residue.atom_range[end_res_idx].end};

            md_array_push(gro->mol.chain.id, add_label(chain_id, gro).ptr, alloc);
            md_array_push(gro->mol.chain.residue_range, res_range, alloc);
            md_array_push(gro->mol.chain.atom_range, atom_range, alloc);

            gro->mol.chain.count = chain_idx;

            md_array_ensure(gro->mol.atom.chain_idx, gro->mol.chain.count, alloc);
            for (int64_t ci = 0; ci < gro->mol.chain.count; ++ci) {
                for (int64_t ai = (int64_t)gro->mol.chain.atom_range[ci].beg; ai < (int64_t)gro->mol.chain.atom_range[ci].end; ++ai) {
                    gro->mol.atom.chain_idx[ai] = ci;
                }
            }
        }
    }

    return true;
}

bool md_gro_molecule_free(struct md_gro_molecule_t* mol, struct md_allocator_i* alloc) {

}

#ifdef __cplusplus
}
#endif