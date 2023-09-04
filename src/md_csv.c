#include <md_csv.h>

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
    bool has_field_names = false;
    int64_t num_fields = 0;
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
        for (int64_t i = 0; i < num_fields; ++i) {
            extract_token_delim(&tok, &line, ',');
            tok = str_trim(tok);
            md_array_push(csv->field_names, tok, alloc);
        }
    }

    for (int64_t i = 0; i < num_fields; ++i) {
        md_array_push(csv->field_values, NULL, alloc);
    }

    while (md_buffered_reader_extract_line(&line, reader)) {
        int64_t i;
        for (i = 0; i < num_fields; ++i) {
            if (extract_token_delim(&tok, &line, ',')) {
                break;
            }
            tok = str_trim(tok);
            if (is_float(tok)) {
                md_array_push(csv->field_values[i], (float)parse_float(tok), alloc);
            } else {
                MD_LOG_ERROR("CSV: Unable to parse float from token: '%.*s'", (int)tok.len, tok.ptr);
                return false;
            }
        }
        if (i < num_fields || !str_empty(line)) {
            MD_LOG_ERROR("CSV: Number of columns in row %i does not match the first row", (int)num_fields);
            return false;
        }
    }

    return true;
}

bool md_csv_parse_str (md_csv_t* csv, str_t in_str, struct md_allocator_i* alloc) {
	md_buffered_reader_t reader = md_buffered_reader_from_str(in_str);
    return parse(csv, &reader, alloc);
}

bool md_csv_parse_file(md_csv_t* csv, str_t in_path, struct md_allocator_i* alloc) {
    md_file_o* file = md_file_open(in_path, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        size_t cap = MEGABYTES(1);
        char* buf = md_alloc(md_heap_allocator, cap);
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        bool result = parse(csv, &reader, alloc);
        md_free(md_heap_allocator, buf, cap);
        return result;
    } else {
        MD_LOG_ERROR("CSV: Failed to open file '%.*s'", (int)in_path.len, in_path.ptr);
    }
    return false;
}

static void write(md_strb_t* sb, const float* field_values[], const str_t field_names[], int64_t num_fields, int64_t num_values) {
    if (field_names) {
        for (int64_t i = 0; i < num_fields; ++i) {
            md_strb_fmt(sb, "%.*s", (int)field_names[i].len, field_names[i].ptr);
            const char c = i < num_fields - 1 ? ',' : '\n';
            md_strb_push_char(sb, c);
        }
    }
    for (int row = 0; row < num_values; ++row) {
        for (int64_t col = 0; col < num_fields; ++col) {
            md_strb_fmt(sb, "%f", field_values[col][row]);
            const char c = col < num_fields - 1 ? ',' : '\n';
            md_strb_push_char(sb, c);
        }
    }
}

str_t md_csv_write_to_str (const float* field_values[], const str_t field_names[], int64_t num_fields, int64_t num_values, struct md_allocator_i* alloc) {
    str_t result = {0};
    if (field_values && num_fields > 0 && num_values > 0) {
        md_strb_t sb = md_strb_create(md_heap_allocator);
        write(&sb, field_values, field_names, num_fields, num_values);
        result = str_copy(md_strb_to_str(&sb), alloc);
        md_strb_free(&sb);
    }
    return result;
}

bool md_csv_write_to_file(const float* field_values[], const str_t field_names[], int64_t num_fields, int64_t num_values, str_t path) {
    if (field_values && num_fields > 0 && num_values > 0) {
        md_file_o* file = md_file_open(path, MD_FILE_WRITE | MD_FILE_BINARY);
        if (file) {
            md_strb_t sb = md_strb_create(md_heap_allocator);
            write(&sb, field_values, field_names, num_fields, num_values);
            str_t str = md_strb_to_str(&sb);
            const int64_t written_bytes = md_file_write(file, str.ptr, str.len);
            md_strb_free(&sb);
            md_file_close(file);
            
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
        for (int64_t i = 0; i < md_array_size(csv->field_names); ++i) {
            str_free(csv->field_names[i], alloc);
        }
    }
    
    for (int64_t i = 0; i < md_array_size(csv->field_values); ++i) {
        md_array_free(csv->field_values[i], alloc);
    }
    MEMSET(csv, 0, sizeof(md_csv_t));
}