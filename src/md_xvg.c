#include <md_xvg.h>

#include "core/md_allocator.h"
#include "core/md_os.h"
#include "core/md_log.h"
#include "core/md_str_builder.h"
#include "core/md_array.h"
#include <core/md_parse.h>

#include <time.h>

static bool open_file(md_file_o* file, str_t path) {
	ASSERT(file);

	if (str_empty(path)) {
		MD_LOG_ERROR("XVG: Path was empty");
		return false;
	}

	file = md_file_open(path, MD_FILE_READ);
	if (!file) {
		MD_LOG_ERROR("XVG: Could not open file: '%.*s'", path.len, path.ptr);
		return false;
	}

	return true;
}

str_t md_xvg_format_header(str_t title, str_t xaxis_label, str_t yaxis_label, size_t num_legends, const str_t* legends, struct md_allocator_i* str_alloc) {
	ASSERT(str_alloc);
	
	md_strb_t sb = md_strb_create(md_heap_allocator);

	time_t t;
	struct tm* info;
	time(&t);
	info = localtime(&t);
	
	md_strb_fmt(&sb, "# This file was created %s", asctime(info));
	md_strb_fmt(&sb, "# Created by:\n");
	md_strb_fmt(&sb, "# VIAMD \n");

	if (!str_empty(title)) {
		md_strb_fmt(&sb, "@    title \"%.*s\"\n", (int)title.len, title.ptr);
	} else {
		md_strb_fmt(&sb, "@    title \"VIAMD export\"\n");
	}
	md_strb_fmt(&sb, "@    xaxis  label \"%.*s\"\n", (int)xaxis_label.len, xaxis_label.ptr);
	md_strb_fmt(&sb, "@    yaxis  label \"%.*s\"\n", (int)yaxis_label.len, yaxis_label.ptr);
	md_strb_fmt(&sb, "@ TYPE xy\n");
	md_strb_fmt(&sb, "@ view 0.15, 0.15, 0.75, 0.85\n");
	md_strb_fmt(&sb, "@ legend on\n");
	md_strb_fmt(&sb, "@ legend box on\n");
	md_strb_fmt(&sb, "@ legend loctype view\n");
	md_strb_fmt(&sb, "@ legend 0.78, 0.8\n");
	md_strb_fmt(&sb, "@ legend length 2\n");
	
	for (size_t i = 0; i < num_legends; ++i) {
		md_strb_fmt(&sb, "@ s%d legend \"%.*s\"\n", (int)i, (int)legends[i].len, legends[i].ptr);
	}

	str_t result = str_copy(md_strb_to_str(sb), str_alloc);

	md_strb_free(&sb);
	return result;
}

str_t md_xvg_format(str_t header, size_t num_fields, size_t num_values, const float* const field_values[], struct md_allocator_i* str_alloc) {
	ASSERT(str_alloc);
	str_t str = {0};

	if (num_fields == 0) {
        MD_LOG_ERROR("XVG: num_fields is zero");
        goto done;
	}

	if (num_values == 0) {
        MD_LOG_ERROR("XVG: num_values is zero");
        goto done;
    }

	if (!field_values) {
        MD_LOG_ERROR("XVG: field_values is null even though num_fields > 0");
        goto done;
    }

	md_strb_t sb = {0};
	md_strb_init(&sb, md_heap_allocator);
	md_strb_push_str(&sb, header);

	for (size_t i = 0; i < num_values; ++i) {
		for (size_t j = 0; j < num_fields; ++j) {
			md_strb_fmt(&sb, "%12.6f", field_values[j][i]);
			const char c = (j < num_fields - 1) ? ' ' : '\n';
			md_strb_push_char(&sb, c);
		}
	}

	str = str_copy(md_strb_to_str(sb), str_alloc);
	md_strb_free(&sb);
done:
	return str;
}

str_t md_xvg_to_str(const md_xvg_t* xvg, struct md_allocator_i* alloc) {
	ASSERT(alloc);
	str_t header = md_xvg_format_header(xvg->header_info.title, xvg->header_info.xaxis_label, xvg->header_info.yaxis_label, md_array_size(xvg->header_info.legends), xvg->header_info.legends, md_heap_allocator);
	str_t result = md_xvg_format(header, md_array_size(xvg->fields), xvg->fields ? md_array_size(xvg->fields[0]) : 0, (const float* const*)xvg->fields, alloc);
	str_free(header, md_heap_allocator);
	return result;
}

str_t remove_quotes(str_t str) {
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
	return (str_t) {beg, end-beg};
}

str_t concat_tokens(str_t first, str_t last) {
	str_t str = {first.ptr, (last.ptr + last.len) - first.ptr};
	return str;
}

bool parse_header_info(md_xvg_header_info_t* header_info, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t line;

	bool result = false;
	md_strb_t sb = md_strb_create(md_heap_allocator);

	while (md_buffered_reader_peek_line(&line, reader)) {
		if (line.ptr[0] != '#' && line.ptr[0] != '@') break;
		md_buffered_reader_skip_line(reader);
		md_strb_push_str(&sb, line);
		md_strb_push_char(&sb, '\n');
	}

	if (md_strb_len(sb) == 0) {
		MD_LOG_ERROR("XVG: Expected header section");
		goto done;
	}

	// Extract string from reader
	str_t header = str_copy(md_strb_to_str(sb), alloc);
	header_info->header = header;

	// read comment
	while (str_peek_line(&line, &header)) {
		if (line.ptr[0] != '#') break;
		str_skip_line(&header);
	}
	
	header_info->comment = (str_t){header_info->header.ptr, line.ptr - header_info->header.ptr};
	
	str_t tok[16];
	
	// Read meta
	while (str_extract_line(&line, &header)) {
		int64_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok && str_eq(tok[0], STR_LIT("@"))) {
			if		  (str_eq(tok[1], STR_LIT("title"))) {
				header_info->title = remove_quotes(concat_tokens(tok[2], tok[num_tok - 1]));
			} else if (str_eq(tok[1], STR_LIT("xaxis"))) {
				header_info->xaxis_label = remove_quotes(concat_tokens(tok[3], tok[num_tok - 1]));
			} else if (str_eq(tok[1], STR_LIT("yaxis"))) {
				header_info->yaxis_label = remove_quotes(concat_tokens(tok[3], tok[num_tok - 1]));
			} else if (str_eq(tok[2], STR_LIT("legend"))) {
				md_array_push(header_info->legends, remove_quotes(concat_tokens(tok[3], tok[num_tok - 1])), alloc);
			}
		}
	}

	header_info->num_legends = md_array_size(header_info->legends);
	result = true;
done:
	md_strb_free(&sb);
	return result;
}

bool parse(md_xvg_t* xvg, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	if (!parse_header_info(&xvg->header_info, reader, alloc)) {
		return false;
	}

	str_t line;
	if (!md_buffered_reader_extract_line(&line, reader)) {
		MD_LOG_ERROR("XVG: Expected field values");
		return false;
	}

	md_array(md_array(float)) fields = 0;

	// Go over the first line and make this our reference field line.
	// Every following line is expected to contain as many values as this line.
	str_t tok = {0};
	while (extract_token(&tok, &line)) {
		if (!is_float(tok)) {
			MD_LOG_ERROR("XVG: Invalid float among field values: '"STR_FMT"'", STR_ARG(tok));
			return false;
		}
		md_array(float) values = 0;
		md_array_push(values, (float)parse_float(tok), alloc);
		md_array_push(fields, values, alloc);
	}

	const size_t num_fields = md_array_size(fields);
	
	while (md_buffered_reader_extract_line(&line, reader)) {
		for (size_t i = 0; i < num_fields; ++i) {
			if (!extract_token(&tok, &line)) {
				MD_LOG_ERROR("XVG: To few entries in field line, expected %i, got %i", (int)num_fields, i);
				return false;
			}

			if (!is_float(tok)) {
				MD_LOG_ERROR("XVG: Invalid float: '"STR_FMT"'", STR_ARG(tok));
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

bool md_xvg_parse_str(md_xvg_t* xvg, str_t str, md_allocator_i* alloc) {
	ASSERT(alloc);

	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	return parse(xvg, &reader, alloc);
}

bool md_xvg_parse_file(md_xvg_t* xvg, str_t path, md_allocator_i* alloc) {
	ASSERT(alloc);

	md_file_o* file = md_file_open(path, MD_FILE_READ);

	if (!file) {
		MD_LOG_ERROR("XVG: Failed to deserialize file, file could not be opened '"STR_FMT"'", STR_ARG(path));
		return false;
	}

	const size_t cap = MEGABYTES(1);
	char* buf = md_alloc(md_heap_allocator, cap);
	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

	bool result = parse(xvg, &reader, alloc);
	md_file_close(file);
	return result;
}

void md_xvg_free(md_xvg_t* xvg, md_allocator_i* alloc) {
	ASSERT(xvg);
	ASSERT(alloc);
	if (!str_empty(xvg->header_info.header)) {
		str_free(xvg->header_info.header, alloc);
	}
	md_array_free(xvg->header_info.legends, alloc);
	md_array_free(xvg->fields, alloc);
}
