#include <md_csv.h>
#include <md_io.h>
#include <md_system.h>
#include <md_run.h>
#include <core/md_unit.h>
#include <stdio.h>
#include <math.h>

#include "core/md_allocator.h"
#include "core/md_file.h"
#include "core/md_path.h"
#include "core/md_format.h"
#include "core/md_log.h"
#include "core/md_writer.h"
#include "core/md_array.h"
#include <core/md_parse.h>
#include <core/md_reader.h>

// A field without the blanks and the quotes around it
static md_str_t csv_field(md_str_t tok) {
    tok = md_str_trim(tok);
    if (tok.len >= 2 && ((tok.ptr[0] == '"' && tok.ptr[tok.len - 1] == '"') || (tok.ptr[0] == '\'' && tok.ptr[tok.len - 1] == '\''))) {
        tok = md_str_trim(md_str_substr(tok, 1, tok.len - 2));
    }
    return tok;
}

// The next field of a row, up to the delimiter or the end, which may be empty ("1,,3"; "1," ends in an empty one):
// false past the last one. *more is true at the start of a row.
static bool csv_next_field(md_str_t* out, md_str_t* row, char delim, bool* more) {
    if (!*more) return false;
    size_t i = 0;
    while (i < row->len && row->ptr[i] != delim) ++i;
    *out = (md_str_t){ row->ptr, i };
    if (i < row->len) {
        *row = md_str_substr(*row, i + 1, SIZE_MAX);
    } else {
        *row = (md_str_t){ row->ptr + row->len, 0 };
        *more = false;
    }
    return true;
}

// The next line which holds anything: blank lines and comments ('#') are passed over
static bool csv_next_line(md_reader_i* reader, md_str_t* line) {
    while (md_reader_extract_line(reader, line)) {
        const md_str_t t = md_str_trim(*line);
        if (t.len && t.ptr[0] != '#') return true;
    }
    return false;
}

// The rows of the file: a first row of names (any field of it not a number) or of values, which says how many
// columns there are, then rows of values. The columns are separated by commas, or semicolons or tabs when the first
// row has no comma (as spreadsheets with a decimal comma write them). An empty field is a value that is missing:
// NaN. Blank lines and comments ('#') are passed over.
static bool parse(md_csv_t* csv, md_reader_i* reader, struct md_allocator_i* alloc) {
    MD_ASSERT(csv);
    MD_ASSERT(reader);
    MD_ASSERT(alloc);
    
    md_str_t line;
    if (!csv_next_line(reader, &line)) {
        MD_LOG_ERROR("CSV: Failed to read any line");
        return false;
    }

    if (csv->field_values || csv->field_names) {
        MD_LOG_DEBUG("CSV: potential memory leak, csv structure was not empty");
    }

    md_memset(csv, 0, sizeof(md_csv_t));

    const char delim = md_str_contains_char(line, ',') ? ',' : md_str_contains_char(line, ';') ? ';' : md_str_contains_char(line, '\t') ? '\t' : ',';

    // Read first row explicitly to determine the expected number of columns
    // Also potentially extract field names from first row
    md_str_t tok;
    md_str_t first_line = md_str_trim(line);
    bool has_field_names = false;
    size_t num_fields = 0;
    bool more = true;
    while (csv_next_field(&tok, &first_line, delim, &more)) {
        tok = csv_field(tok);
        if (!md_parse_is_float(tok)) {
            has_field_names = true;
        }
        num_fields += 1;
    }

    if (num_fields == 0) {
        MD_LOG_ERROR("CSV: No fields found");
        return false;
    }

    for (size_t i = 0; i < num_fields; ++i) {
        md_array_push(csv->field_values, NULL, alloc);
    }

    bool first = true;
    do {
        if (first && has_field_names) {
            MD_LOG_INFO("CSV: First row contains non-numeric values, assuming field names");
            md_str_t row = md_str_trim(line);
            more = true;
            for (size_t i = 0; i < num_fields; ++i) {
                csv_next_field(&tok, &row, delim, &more);
                md_array_push(csv->field_names, md_alloc_str(alloc, csv_field(tok)), alloc);
            }
            first = false;
            continue;
        }
        first = false;
        md_str_t row = md_str_trim(line);
        more = true;
        size_t i;
        for (i = 0; i < num_fields; ++i) {
            if (!csv_next_field(&tok, &row, delim, &more)) {
                break;
            }
            tok = csv_field(tok);
            if (md_str_empty(tok)) {
                md_array_push(csv->field_values[i], NAN, alloc);
            } else if (md_parse_is_float(tok)) {
                md_array_push(csv->field_values[i], (float)md_parse_float(tok), alloc);
            } else {
                MD_LOG_ERROR("CSV: Unable to parse float from token: '"MD_STR_FMT"'", MD_STR_ARG(tok));
                return false;
            }
        }
        if (i < num_fields || more) {
            MD_LOG_ERROR("CSV: A row has %s columns than the first row's %zu", i < num_fields ? "fewer" : "more", num_fields);
            return false;
        }
    } while (csv_next_line(reader, &line));

    csv->num_fields = num_fields;
    csv->num_values = md_array_size(csv->field_values[0]);

    return true;
}

bool md_csv_parse_str (md_csv_t* csv, md_str_t in_str, struct md_allocator_i* alloc) {
	md_reader_i reader = md_reader_from_str(in_str);
    return parse(csv, &reader, alloc);
}

static bool md_csv_parse_file(md_csv_t* csv, md_str_t in_path, struct md_allocator_i* alloc) {
    md_file_t file = {0};
    if (md_file_open(&file, in_path, MD_FILE_READ)) {
        size_t cap = MD_MEGABYTES(1);
        md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
        md_file_reader_t reader;
        md_file_reader_init(&reader, file, cap, md_temp_allocator(temp_scope));
        bool result = parse(csv, &reader.r, alloc);
        md_file_reader_close(&reader);
        md_temp_end(temp_scope);
        md_file_close(&file);
        return result;
    } else {
        MD_LOG_ERROR("CSV: Failed to open file '"MD_STR_FMT"'", MD_STR_ARG(in_path));
    }
    return false;
}

void md_csv_write(md_writer_i* w, const float* const field_values[], const md_str_t field_names[], size_t num_fields, size_t num_values) {
    MD_ASSERT(w);
    if (!field_values || num_fields == 0) {
        return;
    }
    if (field_names) {
        for (size_t col = 0; col < num_fields; ++col) {
            md_write_str(w, field_names[col]);
            md_write_char(w, col < num_fields - 1 ? ',' : '\n');
        }
    }
    for (size_t row = 0; row < num_values; ++row) {
        for (size_t col = 0; col < num_fields; ++col) {
            md_write_f64(w, field_values[col][row], 6);     // as "%f" writes it
            md_write_char(w, col < num_fields - 1 ? ',' : '\n');
        }
    }
}

md_str_t md_csv_write_to_str (const float* const field_values[], const md_str_t field_names[], size_t num_fields, size_t num_values, struct md_allocator_i* alloc) {
    md_str_t result = {0};
    if (field_values && num_fields > 0 && num_values > 0) {
        md_str_writer_t sw = md_str_writer_create(alloc);
        md_csv_write(&sw.w, field_values, field_names, num_fields, num_values);
        result = md_str_writer_release(&sw);
    }
    return result;
}

bool md_csv_write_to_file(const float* const field_values[], const md_str_t field_names[], size_t num_fields, size_t num_values, md_str_t path) {
    if (field_values && num_fields > 0 && num_values > 0) {
        md_file_t file = {0};
        if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            // Streamed through a staging buffer, rather than built whole in memory first
            md_file_writer_t fw = md_file_writer_create(file, NULL, 0);
            md_csv_write(&fw.w, field_values, field_names, num_fields, num_values);
            const bool written = md_file_writer_close(&fw);
            md_file_close(&file);

            if (written) {
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
    MD_ASSERT(csv);
    MD_ASSERT(alloc);

    if (csv->field_names) {
        for (size_t i = 0; i < md_array_size(csv->field_names); ++i) {
            md_free_str(alloc, csv->field_names[i]);
        }
    }
    
    for (size_t i = 0; i < md_array_size(csv->field_values); ++i) {
        md_array_free(csv->field_values[i], alloc);
    }
    md_memset(csv, 0, sizeof(md_csv_t));
}

// ### RUN ###

// The unit in parentheses in a label - "Time (ps)", "Energy (kJ/mol)" - or none
static md_unit_t csv_label_unit(md_str_t label) {
    size_t beg, end;
    md_unit_t unit = md_unit_none();
    if (((beg = md_str_find_char(label, '(')) != MD_STR_NPOS) && ((end = md_str_find_char(label, ')')) != MD_STR_NPOS) && end > beg) {
        md_unit_t parsed;
        if (md_unit_parse(&parsed, md_str_substr(label, beg + 1, end - beg - 1))) {
            unit = parsed;
        }
    }
    return unit;
}

// "<kind>/<file stem>", the stem folded to lower case letters, digits and '_'
static md_str_t csv_group(char* buf, size_t cap, const char* kind, md_str_t filename) {
    md_str_t file = md_path_stem(filename);
    size_t len = (size_t)md_snprintf(buf, cap, "%s/", kind);
    const size_t base = len;
    bool sep = false;
    for (size_t i = 0; i < file.len && len + 2 < cap; ++i) {
        char c = file.ptr[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9'))) { sep = (len > base); continue; }
        if (sep) { buf[len++] = '_'; sep = false; }
        buf[len++] = c;
    }
    if (len == base) len += (size_t)md_snprintf(buf + len, cap - len, "data");
    buf[len] = '\0';
    return (md_str_t){ buf, len };
}

bool md_csv_system_supplement_from_file(struct md_system_t* sys, md_str_t filename, const md_load_desc_t* desc) {
    const md_str_t run = desc ? desc->run : (md_str_t){0};
    MD_ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool result = false;

    md_csv_t csv = {0};
    if (!md_csv_parse_file(&csv, filename, alloc) || csv.num_fields == 0 || csv.num_values == 0) {
        MD_LOG_ERROR("CSV: failed to read '" MD_STR_FMT "'", MD_STR_ARG(filename));
        goto done;
    }

    {
        // The first column is time when its name says so; each column's unit is in its name.
        const bool has_time = csv.field_names && md_str_begins_with_ignore_case(md_str_trim(csv.field_names[0]), MD_STR_LIT("time"));
        const size_t first = has_time ? 1 : 0;
        const size_t C = csv.num_fields - first;
        if (C == 0) {
            MD_LOG_ERROR("CSV: '" MD_STR_FMT "' holds a time column and nothing else", MD_STR_ARG(filename));
            goto done;
        }
        md_str_t*        names   = md_temp_alloc(temp, C * sizeof(md_str_t));
        md_unit_t*    units   = md_temp_alloc(temp, C * sizeof(md_unit_t));
        const float** columns = md_temp_alloc(temp, C * sizeof(float*));
        for (size_t k = 0; k < C; ++k) {
            names[k]   = csv.field_names ? csv.field_names[first + k] : (md_str_t){0};
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
            .source_path = (md_str_t){ path_buf, path_len },
        };
        result = md_run_publish_series(sys, run, &desc);
    }

done:
    md_temp_end(temp);
    return result;
}
