#include <md_csv.h>
#include <md_system.h>
#include <core/md_unit.h>
#include <stdio.h>

#include "core/md_allocator.h"
#include "core/md_os.h"
#include "core/md_log.h"
#include "core/md_str_builder.h"
#include "core/md_array.h"
#include <core/md_parse.h>

static bool parse(md_csv_t* csv, md_buffered_reader_t* reader, struct md_allocator_i* alloc) {
    ASSERT(csv);
    ASSERT(reader);
    ASSERT(alloc);
    
    str_t line;
    if (!md_buffered_reader_extract_line(&line, reader)) {
        MD_LOG_ERROR("CSV: Failed to read any line");
        return false;
    }

    if (csv->field_values || csv->field_names) {
        MD_LOG_DEBUG("CSV: potential memory leak, csv structure was not empty");
    }

    MEMSET(csv, 0, sizeof(md_csv_t));

    // Read first row explicitly to determine the expected number of columns
    // Also potentially extract field names from first row
    str_t tok;
    str_t first_line = line;
    bool has_field_names = true;
    size_t num_fields = 0;
    while (extract_token_delim(&tok, &first_line, ',')) {
        tok = str_trim(tok);
        if (is_float(tok)) {
            has_field_names = false;
        }
        num_fields += 1;
    }

    if (num_fields == 0) {
        MD_LOG_ERROR("CSV: No fields found");
        return false;
    }
    
    if (has_field_names) {
        MD_LOG_INFO("CSV: First row contains non-numeric values, assuming field names");
        // Read first line as field names
        for (size_t i = 0; i < num_fields; ++i) {
            extract_token_delim(&tok, &line, ',');
            str_t name = str_copy(str_trim(tok), alloc);
            md_array_push(csv->field_names, name, alloc);
        }
    }

    for (size_t i = 0; i < num_fields; ++i) {
        md_array_push(csv->field_values, NULL, alloc);
    }

    while (md_buffered_reader_extract_line(&line, reader)) {
        size_t i;
        for (i = 0; i < num_fields; ++i) {
            if (!extract_token_delim(&tok, &line, ',')) {
                break;
            }
            tok = str_trim(tok);
            if (is_float(tok)) {
                md_array_push(csv->field_values[i], (float)parse_float(tok), alloc);
            } else {
                MD_LOG_ERROR("CSV: Unable to parse float from token: '"STR_FMT"'", STR_ARG(tok));
                return false;
            }
        }
        if (i < num_fields || !str_empty(line)) {
            MD_LOG_ERROR("CSV: Number of columns in row %zu does not match the first row", num_fields);
            return false;
        }
    }

    csv->num_fields = num_fields;
    csv->num_values = md_array_size(csv->field_values[0]);

    return true;
}

bool md_csv_parse_str (md_csv_t* csv, str_t in_str, struct md_allocator_i* alloc) {
	md_buffered_reader_t reader = md_buffered_reader_from_str(in_str);
    return parse(csv, &reader, alloc);
}

bool md_csv_parse_file(md_csv_t* csv, str_t in_path, struct md_allocator_i* alloc) {
    md_file_t file = {0};
    if (md_file_open(&file, in_path, MD_FILE_READ)) {
        size_t cap = MEGABYTES(1);
        md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
        char* buf = md_temp_alloc(temp_scope, cap);
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        bool result = parse(csv, &reader, alloc);
        md_temp_end(temp_scope);
        md_file_close(&file);
        return result;
    } else {
        MD_LOG_ERROR("CSV: Failed to open file '"STR_FMT"'", STR_ARG(in_path));
    }
    return false;
}

static void write(md_strb_t* sb, const float* const field_values[], const str_t field_names[], size_t num_fields, size_t num_values) {
    if (field_names) {
        for (size_t i = 0; i < num_fields; ++i) {
            md_strb_fmt(sb, STR_FMT, STR_ARG(field_names[i]));
            const char c = i < num_fields - 1 ? ',' : '\n';
            md_strb_push_char(sb, c);
        }
    }
    for (size_t row = 0; row < num_values; ++row) {
        for (size_t col = 0; col < num_fields; ++col) {
            md_strb_fmt(sb, "%f", field_values[col][row]);
            const char c = col < num_fields - 1 ? ',' : '\n';
            md_strb_push_char(sb, c);
        }
    }
}

str_t md_csv_write_to_str (const float* const field_values[], const str_t field_names[], size_t num_fields, size_t num_values, struct md_allocator_i* alloc) {
    str_t result = {0};
    if (field_values && num_fields > 0 && num_values > 0) {
        md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
        md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
        md_strb_t sb = md_strb_create(temp_alloc);
        write(&sb, field_values, field_names, num_fields, num_values);
        result = str_copy(md_strb_to_str(sb), alloc);
        md_temp_end(temp_scope);
    }
    return result;
}

bool md_csv_write_to_file(const float* const field_values[], const str_t field_names[], size_t num_fields, size_t num_values, str_t path) {
    if (field_values && num_fields > 0 && num_values > 0) {
        md_file_t file = {0};
        if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            md_temp_scope_t temp_scope = md_temp_begin();
            md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
            md_strb_t sb = md_strb_create(temp_alloc);
            write(&sb, field_values, field_names, num_fields, num_values);
            str_t str = md_strb_to_str(sb);
            const size_t written_bytes = md_file_write(file, str.ptr, str.len);
            md_temp_end(temp_scope);
            md_file_close(&file);
            
            if (written_bytes == str.len) {
                return true;
            } else {
                MD_LOG_ERROR("CSV: Unexpected error, some bytes were not written");
            }
        } else {
            MD_LOG_ERROR("CSV: File could not be opened for writing: '%.*s'", (int)path.len, path.ptr);
        }
    }
    return false;
}

void md_csv_free(md_csv_t* csv, struct md_allocator_i* alloc) {
    ASSERT(csv);
    ASSERT(alloc);

    if (csv->field_names) {
        for (size_t i = 0; i < md_array_size(csv->field_names); ++i) {
            str_free(csv->field_names[i], alloc);
        }
    }
    
    for (size_t i = 0; i < md_array_size(csv->field_values); ++i) {
        md_array_free(csv->field_values[i], alloc);
    }
    MEMSET(csv, 0, sizeof(md_csv_t));
}

// ### RUN ###

// The unit in parentheses in a label - "Time (ps)", "Energy (kJ/mol)" - or none
static md_unit_t csv_label_unit(str_t label) {
    size_t beg, end;
    md_unit_t unit = md_unit_none();
    if (str_find_char(&beg, label, '(') && str_find_char(&end, label, ')') && end > beg) {
        md_unit_t parsed;
        if (md_unit_parse(&parsed, str_substr(label, beg + 1, end - beg - 1))) {
            unit = parsed;
        }
    }
    return unit;
}

// "<kind>/<file stem>", the stem folded to lower case letters, digits and '_'
static str_t csv_group(char* buf, size_t cap, const char* kind, str_t filename) {
    str_t file = filename;
    extract_file(&file, filename);
    size_t dot;
    if (str_rfind_char(&dot, file, '.') && dot > 0) {
        file = str_substr(file, 0, dot);
    }
    size_t len = (size_t)snprintf(buf, cap, "%s/", kind);
    const size_t base = len;
    bool sep = false;
    for (size_t i = 0; i < file.len && len + 2 < cap; ++i) {
        char c = file.ptr[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9'))) { sep = (len > base); continue; }
        if (sep) { buf[len++] = '_'; sep = false; }
        buf[len++] = c;
    }
    if (len == base) len += (size_t)snprintf(buf + len, cap - len, "data");
    buf[len] = '\0';
    return (str_t){ buf, len };
}

bool md_csv_system_supplement_from_file(struct md_system_t* sys, str_t filename, str_t run) {
    ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool result = false;

    md_csv_t csv = {0};
    if (!md_csv_parse_file(&csv, filename, alloc) || csv.num_fields == 0 || csv.num_values == 0) {
        MD_LOG_ERROR("CSV: failed to read '" STR_FMT "'", STR_ARG(filename));
        goto done;
    }

    {
        // The first column is time when its name says so; each column's unit is in its name.
        const bool has_time = csv.field_names && str_eq_cstr_n_ignore_case(str_trim(csv.field_names[0]), "time", 4);
        const size_t first = has_time ? 1 : 0;
        const size_t C = csv.num_fields - first;
        if (C == 0) {
            MD_LOG_ERROR("CSV: '" STR_FMT "' holds a time column and nothing else", STR_ARG(filename));
            goto done;
        }
        str_t*        names   = md_temp_alloc(temp, C * sizeof(str_t));
        md_unit_t*    units   = md_temp_alloc(temp, C * sizeof(md_unit_t));
        const float** columns = md_temp_alloc(temp, C * sizeof(float*));
        for (size_t k = 0; k < C; ++k) {
            names[k]   = csv.field_names ? csv.field_names[first + k] : (str_t){0};
            units[k]   = csv.field_names ? csv_label_unit(csv.field_names[first + k]) : md_unit_none();
            columns[k] = csv.field_values[first + k];
        }
        double* time = NULL;
        if (has_time) {
            time = md_temp_alloc(temp, csv.num_values * sizeof(double));
            for (size_t i = 0; i < csv.num_values; ++i) time[i] = csv.field_values[0][i];
        }

        char group_buf[256];
        char path_buf[4096];
        const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
        const md_run_series_desc_t desc = {
            .group       = csv_group(group_buf, sizeof(group_buf), "csv", filename),
            .num_rows    = csv.num_values,
            .num_columns = C,
            .time        = time,
            .time_unit   = has_time ? csv_label_unit(csv.field_names[0]) : md_unit_none(),
            .names       = names,
            .units       = units,
            .columns     = columns,
            .source_path = (str_t){ path_buf, path_len },
        };
        result = md_run_publish_series(sys, run, &desc);
    }

done:
    md_temp_end(temp);
    return result;
}
