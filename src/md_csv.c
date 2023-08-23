#include <md_xvg.h>

#include "core/md_allocator.h"
#include "core/md_os.h"
#include "core/md_log.h"
#include "core/md_str_builder.h"
#include "core/md_array.h"
#include <core/md_parse.h>

static md_array(float) parse(int* out_num_cols, int* out_num_rows, md_buffered_reader_t* reader, struct md_allocator_i* alloc) {
    str_t line;
    if (!md_buffered_reader_extract_line(&line, reader)) {
        MD_LOG_ERROR("CSV: Failed to read any line");
        return NULL;
    }
    int num_cols = 0;
    int num_rows = 0;
    md_array(float) data = 0;

    // Count number of columns
    {
        str_t tok;
        while (extract_token_delim(&tok, &line, ',')) {
            ++num_cols;
            if (is_float(tok)) {
                md_array_push(data, (float)parse_float(tok), alloc);
            } else {
                MD_LOG_ERROR("CSV: Unexpected token in data: '%.*s'", (int)tok.len, tok.ptr);
                goto done;
            }
        }
    }
    num_rows = 1;

    while (md_buffered_reader_extract_line(&line, reader)) {
        str_t tok;
        int i;
        for (i = 0; i < num_cols; ++i) {
            if (extract_token_delim(&tok, &line, ',')) {
                break;
            }
            if (is_float(tok)) {
                md_array_push(data, (float)parse_float(tok), alloc);
            } else {
                MD_LOG_ERROR("CSV: Unexpected token in data: '%.*s'", (int)tok.len, tok.ptr);
                goto done;
            }
        }
        if (i < num_cols || !str_empty(line)) {
            MD_LOG_ERROR("CSV: Number of columns in row %i does not match the first row", num_rows);
            goto done;
        }
        num_rows++;
    }

    *out_num_cols = num_cols;
    *out_num_rows = num_rows;
    done:
    return data;
}

md_array(float) md_csv_parse_str(int* out_num_cols, int* out_num_rows, str_t in_str, struct md_allocator_i* alloc) {
	md_buffered_reader_t reader = md_buffered_reader_from_str(in_str);
    return parse(out_num_cols, out_num_rows, &reader, alloc);
}

md_array(float) md_csv_parse_file(int* out_num_cols, int* out_num_rows, str_t in_path, struct md_allocator_i* alloc) {
    md_file_o* file = md_file_open(in_path, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        size_t cap = MEGABYTES(1);
        char* buf = md_alloc(default_allocator, cap);
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);
        md_array(float) result = parse(out_num_cols, out_num_rows, &reader, alloc);
        md_free(default_allocator, buf, cap);
        return result;
    }
    return NULL;
}

static void write(md_strb_t* sb, const float* data, int num_cols, int num_rows) {
    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {
            md_strb_fmt(sb, "%f", data[row * num_cols + col]);
            if (col < num_cols - 1) {
                md_strb_push_char(sb, ',');
            }
        }
        md_strb_push_char(sb, '\n');
    }
}

str_t md_csv_write_to_str (const float* data, int num_cols, int num_rows, struct md_allocator_i* alloc) {
    str_t result = {0};
    if (data && num_cols > 0 && num_rows > 0) {
        md_strb_t sb = md_strb_create(default_allocator);
        write(&sb, data, num_cols, num_rows);
        result = str_copy(md_strb_to_str(&sb), alloc);
        md_strb_free(&sb);
    }
    return result;
}

bool md_csv_write_to_file(const float* data, int num_cols, int num_rows, str_t path) {
    if (data && num_cols > 0 && num_rows > 0) {
        md_file_o* file = md_file_open(path, MD_FILE_WRITE | MD_FILE_BINARY);
        if (file) {
            md_strb_t sb = md_strb_create(default_allocator);
            write(&sb, data, num_cols, num_rows);
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
