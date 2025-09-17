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

#include <string.h>

#define MD_GRO_MOL_MAGIC  0xbabcebf677bfaf67

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

static bool md_gro_data_parse(md_gro_data_t* data, md_buffered_reader_t* reader, struct md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(reader);
    ASSERT(alloc);
    str_t line;
    str_t tokens[16];

    if (!md_buffered_reader_extract_line(&line, reader)) {
        MD_LOG_ERROR("Failed to parse gro title");
        return false;
    }
    str_copy_to_char_buf(data->title, sizeof(data->title), str_trim(line));

    if (!md_buffered_reader_extract_line(&line, reader)) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        return false;
    }
    
    data->num_atoms = parse_int(str_trim(line));
    if (!data->num_atoms) {
        MD_LOG_ERROR("Failed to parse gro num atoms");
        return false;
    }

    md_array_resize(data->atom_data, data->num_atoms, alloc);

    for (size_t i = 0; i < data->num_atoms; ++i) {
        if (!md_buffered_reader_extract_line(&line, reader)) {
            MD_LOG_ERROR("Failed to extract atom line");
            return false;
        }
        const int64_t num_tokens = extract_float_tokens(tokens, ARRAY_SIZE(tokens), str_substr(line, 20, SIZE_MAX));
        if (num_tokens < 3) {
            MD_LOG_ERROR("Failed to parse atom coordinates, expected at least 3 tokens, got %i", (int)num_tokens);
            return false;
        }
        
        md_gro_atom_t* atom = &data->atom_data[i];
        
        atom->res_id = (int32_t)parse_int(str_trim(str_substr(line, 0, 5)));
        atom->x = (float)parse_float_wide(tokens[0].ptr, tokens[0].len);
        atom->y = (float)parse_float_wide(tokens[1].ptr, tokens[1].len);
        atom->z = (float)parse_float_wide(tokens[2].ptr, tokens[2].len);
        
        str_copy_to_char_buf(atom->res_name,  sizeof(atom->res_name),  str_trim(str_substr(line,  5, 5)));
        str_copy_to_char_buf(atom->atom_name, sizeof(atom->atom_name), str_trim(str_substr(line, 10, 5)));
    }

    if (!md_buffered_reader_extract_line(&line, reader)) {
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

    md_buffered_reader_t line_reader = md_buffered_reader_from_str(str);
    return md_gro_data_parse(data, &line_reader, alloc);
}

bool md_gro_data_parse_file(md_gro_data_t* data, str_t filename, struct md_allocator_i* alloc) {
    bool result = false;
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        const int64_t cap = MEGABYTES(1);
        char* buf = md_alloc(md_get_heap_allocator(), cap);
        
        md_buffered_reader_t line_reader = md_buffered_reader_from_file(buf, cap, file);
        result = md_gro_data_parse(data, &line_reader, alloc);
        
        md_free(md_get_heap_allocator(), buf, cap);
        md_file_close(file);
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

bool md_gro_molecule_init(struct md_molecule_t* mol, const md_gro_data_t* data, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(data);
    ASSERT(alloc);

    MEMSET(mol, 0, sizeof(md_molecule_t));
    md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));

    const size_t capacity = ROUND_UP(data->num_atoms, 16);

    md_array_ensure(mol->atom.x, capacity, alloc);
    md_array_ensure(mol->atom.y, capacity, alloc);
    md_array_ensure(mol->atom.z, capacity, alloc);
    md_array_ensure(mol->atom.type_idx, capacity, alloc);
    md_array_ensure(mol->atom.flags, capacity, alloc);

    md_residue_id_t*    residue_ids   = md_vm_arena_push_array(temp_arena, md_residue_id_t, capacity);
    str_t*              residue_names = md_vm_arena_push_array(temp_arena, str_t, capacity);

    mol->atom.type.count = 0;
    md_atom_type_find_or_add(&mol->atom.type, STR_LIT("Unknown"), 0, 0.0f, 0.0f, alloc);

    for (size_t i = 0; i < data->num_atoms; ++i) {
        const float x = data->atom_data[i].x * 10.0f; // convert from nm to Ångström
        const float y = data->atom_data[i].y * 10.0f; // convert from nm to Ångström
        const float z = data->atom_data[i].z * 10.0f; // convert from nm to Ångström
        str_t atom_name = str_from_cstrn(data->atom_data[i].atom_name, sizeof(data->atom_data[i].atom_name));
        str_t res_name  = str_from_cstrn(data->atom_data[i].res_name,  sizeof(data->atom_data[i].res_name));
        md_residue_id_t res_id = data->atom_data[i].res_id;
        md_atomic_number_t atomic_number = md_atom_infer_atomic_number(atom_name, res_name);

        float mass = md_util_element_atomic_mass(atomic_number);
        float radius = md_util_element_vdw_radius(atomic_number);

        md_atom_type_idx_t type_idx = md_atom_type_find_or_add(&mol->atom.type, atom_name, atomic_number, mass, radius, alloc);

        mol->atom.count += 1;
        md_array_push_no_grow(mol->atom.x, x);
        md_array_push_no_grow(mol->atom.y, y);
        md_array_push_no_grow(mol->atom.z, z);
        md_array_push_no_grow(mol->atom.type_idx, type_idx);
        md_array_push_no_grow(mol->atom.flags, 0);

        residue_ids[i] = res_id;
        residue_names[i] = res_name;
    }

    md_util_residue_infer(&mol->residue, mol->atom.flags, residue_ids, residue_names, mol->atom.count, alloc);
    md_util_residue_infer_flags(&mol->residue, mol->atom.flags, &mol->atom);

    float box[3][3];
    MEMCPY(&box, data->box, sizeof(mat3_t));
    // convert from nm to Ångström
    for (int i = 0; i < 3; ++i) {
        box[i][0] *= 10.0f;
        box[i][1] *= 10.0f;
        box[i][2] *= 10.0f;
    }

    mol->unit_cell = md_util_unit_cell_from_matrix(box);

    md_vm_arena_destroy(temp_arena);

    return true;
}

static bool gro_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_gro_data_t data = {0};
    bool success = false;
    if (md_gro_data_parse_str(&data, str, md_get_heap_allocator())) {
        success = md_gro_molecule_init(mol, &data, alloc);
    }
    md_gro_data_free(&data, md_get_heap_allocator());

    return success;
}

static bool gro_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
    (void)arg;
    md_gro_data_t data = {0};
    bool success = false;
    if (md_gro_data_parse_file(&data, filename, md_get_heap_allocator())) {
        success = md_gro_molecule_init(mol, &data, alloc);
    }
    md_gro_data_free(&data, md_get_heap_allocator());

    return success;
}

static md_molecule_loader_i gro_api = {
    gro_init_from_str,
    gro_init_from_file,
};

md_molecule_loader_i* md_gro_molecule_api(void) {
    return &gro_api;
}
