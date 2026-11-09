#include <md_ndx.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_str.h>

static bool ndx_parse(md_ndx_data_t* data, md_reader_i* r, md_allocator_i* alloc) {
    md_memset(data, 0, sizeof(md_ndx_data_t));
    md_array(md_ndx_group_t) groups = 0;
    md_array(uint32_t) atoms = 0;
    size_t line_num = 0;
    bool ok = true;

    md_str_t line;
    while (ok && md_reader_extract_line(r, &line)) {
        ++line_num;
        size_t semi;
        if (md_str_find_char(&semi, line, ';')) line = md_str_substr(line, 0, semi);
        line = md_str_trim(line);
        if (md_str_empty(line)) continue;

        if (line.ptr[0] == '[') {
            size_t close;
            if (!md_str_find_char(&close, line, ']')) {
                MD_LOG_ERROR("NDX: Line %zu: a group header without ']'", line_num);
                ok = false;
                break;
            }
            // The atoms of the group before go with it
            if (md_array_size(groups)) {
                md_ndx_group_t* g = md_array_last(groups);
                g->count = md_array_size(atoms);
                if (g->count) md_array_push_array(g->atom_idx, atoms, g->count, alloc);
                md_array_shrink(atoms, 0);
            }
            const md_str_t name = md_str_trim(md_str_substr(line, 1, close - 1));
            md_ndx_group_t g = { .name = md_str_empty(name) ? (md_str_t){0} : md_str_copy(name, alloc) };
            md_array_push(groups, g, alloc);
            continue;
        }

        if (!md_array_size(groups)) {
            MD_LOG_ERROR("NDX: Line %zu: atom numbers before the first group header", line_num);
            ok = false;
            break;
        }
        md_str_t tok;
        while (extract_token(&tok, &line)) {
            int64_t v = 0;
            if (md_parse_i64(&v, tok) != tok.len || v < 1 || v > (int64_t)MD_BITSET_MAX_INDEX) {
                MD_LOG_ERROR("NDX: Line %zu: '" MD_STR_FMT "' is not an atom number", line_num, MD_STR_ARG(tok));
                ok = false;
                break;
            }
            const uint32_t idx = (uint32_t)(v - 1);
            md_array_push(atoms, idx, md_get_heap_allocator());
            if (idx + 1 > data->max_atom) data->max_atom = idx + 1;
        }
    }

    if (ok && md_array_size(groups)) {
        md_ndx_group_t* g = md_array_last(groups);
        g->count = md_array_size(atoms);
        if (g->count) md_array_push_array(g->atom_idx, atoms, g->count, alloc);
    }
    md_array_free(atoms, md_get_heap_allocator());

    data->groups = groups;
    data->num_groups = md_array_size(groups);
    if (!ok) {
        md_ndx_data_free(data, alloc);
        return false;
    }
    if (!data->num_groups) {
        MD_LOG_ERROR("NDX: No groups");
        return false;
    }
    return true;
}

bool md_ndx_data_parse_str(md_ndx_data_t* data, md_str_t str, md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    md_reader_i r = md_reader_from_str(str);
    return ndx_parse(data, &r, alloc);
}

bool md_ndx_data_parse_file(md_ndx_data_t* data, md_str_t filename, md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_file_reader_t fr;
    bool result = false;
    if (md_file_reader_open(&fr, filename, MD_KILOBYTES(64), md_temp_allocator(temp))) {
        result = ndx_parse(data, &fr.r, alloc);
        md_file_reader_close(&fr);
    } else {
        MD_LOG_ERROR("NDX: Could not open '" MD_STR_FMT "'", MD_STR_ARG(filename));
    }
    md_temp_end(temp);
    return result;
}

void md_ndx_data_free(md_ndx_data_t* data, md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    for (size_t i = 0; i < md_array_size(data->groups); ++i) {
        if (!md_str_empty(data->groups[i].name)) md_str_free(data->groups[i].name, alloc);
        md_array_free(data->groups[i].atom_idx, alloc);
    }
    md_array_free(data->groups, alloc);
    md_memset(data, 0, sizeof(md_ndx_data_t));
}

bool md_ndx_system_supplement(md_system_t* sys, const md_ndx_data_t* data) {
    MD_ASSERT(sys);
    MD_ASSERT(data);
    if (data->max_atom > sys->atom.count) {
        MD_LOG_ERROR("NDX: The index file names atom %u, the system has %zu atoms: it is of another system", data->max_atom, sys->atom.count);
        return false;
    }
    md_allocator_i* heap = md_get_heap_allocator();
    for (size_t i = 0; i < data->num_groups; ++i) {
        const md_ndx_group_t* g = &data->groups[i];
        md_bitset_t atoms = md_bitset_from_indices(g->atom_idx, g->count, heap);
        md_system_group_add(sys, g->name, atoms, 0);
        md_bitset_free(&atoms, heap);
    }
    MD_LOG_INFO("NDX: %zu groups", data->num_groups);
    return true;
}

bool md_ndx_system_supplement_from_file(md_system_t* sys, md_str_t filename) {
    MD_ASSERT(sys);
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_ndx_data_t data = {0};
    bool result = md_ndx_data_parse_file(&data, filename, arena) && md_ndx_system_supplement(sys, &data);
    md_arena_allocator_destroy(arena);
    return result;
}
