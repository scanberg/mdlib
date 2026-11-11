#include <md_xvg.h>
#include <md_system.h>
#include <md_run.h>
#include <core/md_unit.h>
#include <stdio.h>

#include "core/md_allocator.h"
#include "core/md_file.h"
#include "core/md_path.h"
#include "core/md_format.h"
#include "core/md_log.h"
#include "core/md_writer.h"
#include "core/md_array.h"
#include <core/md_parse.h>
#include <core/md_reader.h>

#include <time.h>

md_str_t md_xvg_format_header(md_str_t title, md_str_t xaxis_label, md_str_t yaxis_label, size_t num_legends, const md_str_t* legends, struct md_allocator_i* str_alloc) {
	MD_ASSERT(str_alloc);
	
	md_temp_scope_t temp_scope = md_temp_begin_avoid(str_alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
	md_str_writer_t sb = md_str_writer_create(temp_alloc);

	time_t t;
	struct tm* info;
	time(&t);
	info = localtime(&t);
	
	md_str_writer_fmt(&sb, "# This file was created %s", asctime(info));
	md_str_writer_fmt(&sb, "# Created by:\n");
	md_str_writer_fmt(&sb, "# VIAMD \n");

	if (!md_str_empty(title)) {
		md_str_writer_fmt(&sb, "@    title \"%.*s\"\n", (int)title.len, title.ptr);
	} else {
		md_str_writer_fmt(&sb, "@    title \"VIAMD export\"\n");
	}
	md_str_writer_fmt(&sb, "@    xaxis  label \"%.*s\"\n", (int)xaxis_label.len, xaxis_label.ptr);
	md_str_writer_fmt(&sb, "@    yaxis  label \"%.*s\"\n", (int)yaxis_label.len, yaxis_label.ptr);
	md_str_writer_fmt(&sb, "@ TYPE xy\n");
	md_str_writer_fmt(&sb, "@ view 0.15, 0.15, 0.75, 0.85\n");
	md_str_writer_fmt(&sb, "@ legend on\n");
	md_str_writer_fmt(&sb, "@ legend box on\n");
	md_str_writer_fmt(&sb, "@ legend loctype view\n");
	md_str_writer_fmt(&sb, "@ legend 0.78, 0.8\n");
	md_str_writer_fmt(&sb, "@ legend length 2\n");
	
	for (size_t i = 0; i < num_legends; ++i) {
		md_str_writer_fmt(&sb, "@ s%d legend \"%.*s\"\n", (int)i, (int)legends[i].len, legends[i].ptr);
	}

	md_str_t result = md_alloc_str(str_alloc, md_str_writer_to_str(&sb));

	md_temp_end(temp_scope);
	return result;
}

md_str_t md_xvg_format(md_str_t header, size_t num_fields, size_t num_values, const float* const field_values[], struct md_allocator_i* str_alloc) {
	MD_ASSERT(str_alloc);
	md_str_t str = {0};

	if (num_fields == 0) {
        MD_LOG_ERROR("XVG: num_fields is zero");
		return str;
	}

	if (num_values == 0) {
        MD_LOG_ERROR("XVG: num_values is zero");
		return str;
    }

	if (!field_values) {
        MD_LOG_ERROR("XVG: field_values is null even though num_fields > 0");
		return str;
    }

	md_temp_scope_t temp_scope = md_temp_begin_avoid(str_alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
	md_str_writer_t sb = {0};
	md_str_writer_init(&sb, temp_alloc);
	md_str_writer_push_str(&sb, header);

	for (size_t i = 0; i < num_values; ++i) {
		for (size_t j = 0; j < num_fields; ++j) {
			md_str_writer_fmt(&sb, "%12.6f", field_values[j][i]);
			const char c = (j < num_fields - 1) ? ' ' : '\n';
			md_str_writer_push_char(&sb, c);
		}
	}

	str = md_alloc_str(str_alloc, md_str_writer_to_str(&sb));
	md_temp_end(temp_scope);
	return str;
}

md_str_t md_xvg_to_str(const md_xvg_t* xvg, struct md_allocator_i* alloc) {
	MD_ASSERT(alloc);
	md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
	md_str_t header = md_xvg_format_header(xvg->header_info.title, xvg->header_info.xaxis_label, xvg->header_info.yaxis_label, md_array_size(xvg->header_info.legends), xvg->header_info.legends, temp_alloc);
	md_str_t result = md_xvg_format(header, md_array_size(xvg->fields), xvg->fields ? md_array_size(xvg->fields[0]) : 0, (const float* const*)xvg->fields, alloc);
	md_temp_end(temp_scope);
	return result;
}

md_str_t remove_quotes(md_str_t str) {
	const char* c = str.ptr;
	const char* beg = str.ptr;
	const char* end = str.ptr + str.len;
	
	while (c != end) {
		if (*c == '"') {
			beg = ++c;
			break;
		}
		++c;
	}
	while (c != end) {
		if (*c == '"') {
			end = c;
			break;
		}
		++c;
	}
	return (md_str_t) {beg, end-beg};
}

md_str_t concat_tokens(md_str_t first, md_str_t last) {
	md_str_t str = {first.ptr, (last.ptr + last.len) - first.ptr};
	return str;
}

bool parse_header_info(md_xvg_header_info_t* header_info, md_reader_i* reader, md_allocator_i* alloc) {
	md_str_t line;

	bool result = false;
	md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
	md_str_writer_t sb = md_str_writer_create(temp_alloc);

	// Comments and meta lines, and the blank lines between them
	while (md_reader_peek_line(reader, &line)) {
		const md_str_t t = md_str_trim(line);
		if (t.len && t.ptr[0] != '#' && t.ptr[0] != '@') break;
		md_reader_skip_line(reader);
		if (!t.len) continue;
		md_str_writer_push_str(&sb, t);
		md_str_writer_push_char(&sb, '\n');
	}

	if (md_str_writer_len(&sb) == 0) {
		MD_LOG_ERROR("XVG: Expected header section");
		goto done;
	}

	// Extract string from reader
	md_str_t header = md_alloc_str(alloc, md_str_writer_to_str(&sb));
	header_info->header = header;

	// read comment
	while (md_parse_peek_line(&line, &header)) {
		if (line.ptr[0] != '#') break;
		md_parse_skip_line(&header);
	}
	
	header_info->comment = (md_str_t){header_info->header.ptr, line.ptr - header_info->header.ptr};
	
	md_str_t tok[16];
	
	// Read meta
	while (md_parse_next_line(&line, &header)) {
		int64_t num_tok = md_parse_next_tokens(tok, MD_COUNTOF(tok), &line);
		if (num_tok >= 3 && md_str_eq(tok[0], MD_STR_LIT("@"))) {
			if		  (md_str_eq(tok[1], MD_STR_LIT("title"))) {
				header_info->title = remove_quotes(concat_tokens(tok[2], tok[num_tok - 1]));
			} else if (md_str_eq(tok[1], MD_STR_LIT("xaxis"))) {
				header_info->xaxis_label = remove_quotes(concat_tokens(tok[3], tok[num_tok - 1]));
			} else if (md_str_eq(tok[1], MD_STR_LIT("yaxis"))) {
				header_info->yaxis_label = remove_quotes(concat_tokens(tok[3], tok[num_tok - 1]));
			} else if (md_str_eq(tok[2], MD_STR_LIT("legend"))) {
				md_array_push(header_info->legends, remove_quotes(concat_tokens(tok[3], tok[num_tok - 1])), alloc);
			}
		}
	}

	header_info->num_legends = md_array_size(header_info->legends);
	result = true;
done:
	md_temp_end(temp_scope);
	return result;
}

bool parse(md_xvg_t* xvg, md_reader_i* reader, md_allocator_i* alloc) {
	if (!parse_header_info(&xvg->header_info, reader, alloc)) {
		return false;
	}

	// The values of the first set: lines of numbers, the first of which says how many columns there are. Blank
	// lines, comments ('#', also after the values of a line) and meta lines ('@', which GROMACS may write between
	// sets) are passed over; a '&' ends the set, and what follows it, another set, is not read.
	md_array(md_array(float)) fields = 0;
	size_t num_fields = 0;
	md_str_t line;
	while (md_reader_extract_line(reader, &line)) {
		line = md_str_substr(line, 0, md_str_find_char(line, '#'));
		line = md_str_trim(line);
		if (!line.len || line.ptr[0] == '@') continue;
		if (line.ptr[0] == '&') break;

		md_str_t tok = {0};
		if (!fields) {
			// The first line of values makes the fields: every following one is expected to have as many
			while (md_parse_next_token(&tok, &line)) {
				if (!is_float(tok)) {
					MD_LOG_ERROR("XVG: Invalid float among field values: '"MD_STR_FMT"'", MD_STR_ARG(tok));
					return false;
				}
				md_array(float) values = 0;
				md_array_push(values, (float)parse_float(tok), alloc);
				md_array_push(fields, values, alloc);
			}
			num_fields = md_array_size(fields);
			continue;
		}
		for (size_t i = 0; i < num_fields; ++i) {
			if (!md_parse_next_token(&tok, &line)) {
				MD_LOG_ERROR("XVG: To few entries in field line, expected %i, got %i", (int)num_fields, (int)i);
				return false;
			}

			if (!is_float(tok)) {
				MD_LOG_ERROR("XVG: Invalid float: '"MD_STR_FMT"'", MD_STR_ARG(tok));
				return false;
			}

			md_array_push(fields[i], (float)parse_float(tok), alloc);
		}
	}

	if (fields == 0) {
		MD_LOG_ERROR("XVG: Missing field values");
		return false;
	}
	
	xvg->num_fields = num_fields;
	xvg->num_values = fields ? md_array_size(fields[0]) : 0;
	xvg->fields = fields;
	
	return true;
}

bool md_xvg_parse_str(md_xvg_t* xvg, md_str_t str, md_allocator_i* alloc) {
	MD_ASSERT(alloc);

	md_reader_i reader = md_reader_from_str(str);
	return parse(xvg, &reader, alloc);
}

bool md_xvg_parse_file(md_xvg_t* xvg, md_str_t path, md_allocator_i* alloc) {
	MD_ASSERT(alloc);

	md_file_t file = {0};

	if (!md_file_open(&file, path, MD_FILE_READ)) {
		MD_LOG_ERROR("XVG: Failed to deserialize file, file could not be opened '"MD_STR_FMT"'", MD_STR_ARG(path));
		return false;
	}

	const size_t cap = MD_MEGABYTES(1);
	md_temp_scope_t temp = md_temp_begin_avoid(alloc);
	md_file_reader_t reader;
	md_file_reader_init(&reader, file, cap, md_temp_allocator(temp));

	bool result = parse(xvg, &reader.r, alloc);
	md_file_reader_close(&reader);
	md_temp_end(temp);
	md_file_close(&file);
	return result;
}

void md_xvg_free(md_xvg_t* xvg, md_allocator_i* alloc) {
	MD_ASSERT(xvg);
	MD_ASSERT(alloc);
	if (!md_str_empty(xvg->header_info.header)) {
		md_free_str(alloc, xvg->header_info.header);
	}
	md_array_free(xvg->header_info.legends, alloc);
	for (size_t i = 0; i < md_array_size(xvg->fields); ++i) {
		md_array_free(xvg->fields[i], alloc);
	}
	md_array_free(xvg->fields, alloc);
}

// ### RUN ###

// The unit in parentheses in a label - "Time (ps)", "Energy (kJ/mol)" - or none
static md_unit_t xvg_label_unit(md_str_t label) {
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
static md_str_t xvg_group(char* buf, size_t cap, const char* kind, md_str_t filename) {
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

bool md_xvg_system_supplement_from_file(struct md_system_t* sys, md_str_t filename, md_str_t run) {
    MD_ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool result = false;

    md_xvg_t xvg = {0};
    if (!md_xvg_parse_file(&xvg, filename, alloc) || xvg.num_fields == 0 || xvg.num_values == 0) {
        MD_LOG_ERROR("XVG: failed to read '" MD_STR_FMT "'", MD_STR_ARG(filename));
        goto done;
    }

    {
        // The first column is time when the x axis says so; every other column shares the y axis unit.
        const bool has_time = md_str_begins_with_ignore_case(md_str_trim(xvg.header_info.xaxis_label), MD_STR_LIT("time"));
        const size_t first = has_time ? 1 : 0;
        const size_t C = xvg.num_fields - first;
        if (C == 0) {
            MD_LOG_ERROR("XVG: '" MD_STR_FMT "' holds a time column and nothing else", MD_STR_ARG(filename));
            goto done;
        }
        const md_unit_t y_unit = xvg_label_unit(xvg.header_info.yaxis_label);

        md_str_t*       names   = md_temp_alloc(temp, C * sizeof(md_str_t));
        md_unit_t*   units   = md_temp_alloc(temp, C * sizeof(md_unit_t));
        const float** columns = md_temp_alloc(temp, C * sizeof(float*));
        for (size_t k = 0; k < C; ++k) {
            const size_t field = first + k;
            columns[k] = xvg.fields[field];
            if (!has_time && field == 0) {
                names[k] = xvg.header_info.xaxis_label;
                units[k] = xvg_label_unit(xvg.header_info.xaxis_label);
                continue;
            }
            const size_t legend = field - 1;
            names[k] = (legend < md_array_size(xvg.header_info.legends)) ? xvg.header_info.legends[legend] :
                       (C == 1 ? xvg.header_info.yaxis_label : (md_str_t){0});
            units[k] = y_unit;
        }

        double* time = NULL;
        if (has_time) {
            time = md_temp_alloc(temp, xvg.num_values * sizeof(double));
            for (size_t i = 0; i < xvg.num_values; ++i) time[i] = xvg.fields[0][i];
        }

        char group_buf[256];
        char path_buf[4096];
        const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
        const md_run_series_desc_t desc = {
            .group       = xvg_group(group_buf, sizeof(group_buf), "xvg", filename),
            .num_rows    = xvg.num_values,
            .num_columns = C,
            .time        = time,
            .time_unit   = has_time ? xvg_label_unit(xvg.header_info.xaxis_label) : md_unit_none(),
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
