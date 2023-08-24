#include <md_xvg.h>

#include "core/md_allocator.h"
#include "core/md_os.h"
#include "core/md_log.h"
#include "core/md_str_builder.h"
#include "core/md_array.h"
#include <core/md_parse.h>

#include <time.h>

float* md_xvg_row(md_xvg_t* xvg, int row_index) {
	ASSERT(xvg);
	ASSERT(xvg->values);
	return xvg->values + xvg->num_cols * row_index;
}

float md_xvg_value(md_xvg_t* xvg, int col_index, int row_index) {
	ASSERT(xvg);
	ASSERT(xvg->values);
	return xvg->values[xvg->num_cols * row_index + col_index];
}

bool md_xvg_valid(const md_xvg_t* xvg) {
	if (!xvg) {
		MD_LOG_ERROR("XVG: object is NULL");
		return false;
	}

	if (xvg->num_rows < 0) {
        MD_LOG_ERROR("XVG: num_cols is negative");
        return false;
	}

	if (xvg->num_cols < 0) {
		MD_LOG_ERROR("XVG: num_rows is negative");
		return false;
	}

	if (xvg->num_rows > 0 && xvg->num_cols > 0 && !xvg->values) {
        MD_LOG_ERROR("XVG: values are null even though num_fields > 0");
        return false;
    }

	return true;
}

static bool open_file(md_file_o* file, str_t path) {
	ASSERT(file);

	if (str_empty(path)) {
		MD_LOG_ERROR("XVG: Path was empty");
		return false;
	}

	file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("XVG: Could not open file: '%.*s'", path.len, path.ptr);
		return false;
	}

	return true;
}

str_t md_xvg_serialize_to_str(const md_xvg_t* xvg, struct md_allocator_i* alloc) {
    ASSERT(alloc);
	str_t str = {0};

	if (!md_xvg_valid(xvg)) {
		MD_LOG_ERROR("XVG: xvg-object is invalid");
		goto done;
	}

	md_strb_t sb = {0};
	md_strb_init(&sb, default_allocator);

	time_t t;
	struct tm* info;
	time(&t);
	info = localtime(&t);
	
	md_strb_fmt(&sb, "# This file was created %s\n", asctime(info));
	md_strb_fmt(&sb, "# Created by:\n");
	md_strb_fmt(&sb, "# VIAMD \n");

	if (str_empty(xvg->meta_data)) {
		md_strb_fmt(&sb, "@    title \"VIAMD export\"\n");
		md_strb_fmt(&sb, "@    xaxis  label \"Time\"\n");
		md_strb_fmt(&sb, "@    yaxis  label \"Values\"\n");
		md_strb_fmt(&sb, "@ TYPE xy\n");
		md_strb_fmt(&sb, "@ view 0.15, 0.15, 0.75, 0.85\n");
		md_strb_fmt(&sb, "@ legend on\n");
		md_strb_fmt(&sb, "@ legend box on\n");
		md_strb_fmt(&sb, "@ legend loctype view\n");
		md_strb_fmt(&sb, "@ legend 0.78, 0.8\n");
		md_strb_fmt(&sb, "@ legend length 2\n");
	} else {
		md_strb_push_str(&sb, xvg->meta_data);
        if (str_end(xvg->meta_data)[-1] != '\n')
			md_strb_push_char(&sb, '\n');
	}

	for (int i=0; i < xvg->num_rows; ++i) {
		for (int j=0; j < xvg->num_cols; ++j) {
			md_strb_fmt(&sb, "%12.6f", xvg->values + j * xvg->num_rows + i);
            const char c = (j < xvg->num_cols - 1) ? ' ' : '\n';
			md_strb_push_char(&sb, c);
		}
	}
	
	str = str_copy(md_strb_to_str(&sb), alloc);
	md_strb_free(&sb);
done:
	return str;
}

bool deserialize(md_xvg_t* xvg, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	md_allocator_i* temp_alloc = default_allocator;
	md_strb_t sb = md_strb_create(temp_alloc);

	str_t line;
	while (md_buffered_reader_extract_line(&line, reader)) {
		if (line.len > 0 && (line.ptr[0] == '#' || line.ptr[0] == '@')) {
			if (line.ptr[0] == '@') {
				md_strb_push_str(&sb, line);
				md_strb_push_char(&sb, '\n');
			}
		} else {
			break;
		}
	}
	if (line.len == 0) {
		MD_LOG_ERROR("XVG: Expected field values");
		return false;
	}

	md_array(float) values = 0;
	int num_cols = 0;

	// Go over the first line and make this our reference field line.
	// Every following line is expected to contain as many values as this line.
	str_t tok;
	while (extract_token(&tok, &line)) {
		if (!is_float(tok)) {
			MD_LOG_ERROR("XVG: Invalid float among field values: '%*s'", (int)tok.len, tok.ptr);
			goto done;
		}
		num_cols += 1;
		float val = (float)parse_float(tok);
		md_array_push(values, val, temp_alloc);
	}

	int num_rows = 1;

	while (md_buffered_reader_extract_line(&line, reader)) {
		for (int i = 0; i < num_cols; ++i) {
			if (!extract_token(&tok, &line)) {
				MD_LOG_ERROR("XVG: To few entries in field line, expected %i, got %i: '%*s'", num_cols, i);
				goto done;
			}

			if (!is_float(tok)) {
				MD_LOG_ERROR("XVG: Invalid float: '%*s'", (int)tok.len, tok.ptr);
				goto done;
			}

			float val = (float)parse_float(tok);
			md_array_push(values, val, temp_alloc);
		}
		num_rows += 1;
	}

	if (values == 0) {
		MD_LOG_ERROR("XVG: Missing field values", (int)tok.len, tok.ptr);
		goto done;
	}

	// We done, and everything went ok
	md_xvg_init(xvg, num_cols, num_rows, values, md_strb_to_str(&sb), alloc);
done:
	md_array_free(values, temp_alloc);
	md_strb_free(&sb);

	return true;
}

bool md_xvg_deserialize_from_str(md_xvg_t* xvg, str_t str, md_allocator_i* alloc) {
	ASSERT(alloc);

	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	return deserialize(xvg, &reader, alloc);
}

bool md_xvg_init_from_file(md_xvg_t* xvg, str_t path, md_allocator_i* alloc) {
	ASSERT(alloc);

	md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

	if (!file) {
		MD_LOG_ERROR("XVG: Failed to deserialize file, file could not be opened '%.*s'", (int)path.len, path.ptr);
		return false;
	}

	int64_t cap = MEGABYTES(1);
	char* buf = md_alloc(default_allocator, cap);
	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

	bool result = deserialize(xvg, &reader, alloc);
	md_file_close(file);
	return result;
}

bool md_xvg_serialize_to_file(const md_xvg_t* xvg, str_t path) {
	if (!md_xvg_valid(xvg)) {
		MD_LOG_ERROR("XVG: Not valid");
		return false;
	}

	bool success = false;
	md_file_o* file = NULL;
	if (open_file(file, path)) {
		str_t str = md_xvg_serialize_to_str(xvg, default_allocator);
		if (!str_empty(str)) {
			success = md_file_write(file, str_ptr(str), str_len(str)) == str_len(str);
			str_free(str, default_allocator);
		}
		md_file_close(file);
	}

	return success;
}

bool md_xvg_init(md_xvg_t* xvg, int num_cols, int num_rows, const float* values, str_t meta_data, struct md_allocator_i* alloc) {
	ASSERT(xvg);
	ASSERT(alloc);
	
	if (num_cols <= 0) {
		MD_LOG_ERROR("XVG: Invalid number of columns");
		return false;
	}
	if (num_rows <= 0) {
		MD_LOG_ERROR("XVG: Invalid number of rows");
		return false;
	}

	if (xvg->values != NULL) {
		MD_LOG_DEBUG("XVG: values is not empty before initialization, potential memory leak");
	}

	const int64_t num_bytes =num_cols * num_rows * sizeof(float);
	xvg->values = (float*)md_alloc(alloc, num_bytes);
	ASSERT(xvg->values);

	xvg->num_cols = num_cols;
	xvg->num_rows = num_rows;

	if (!str_empty(meta_data)) {
		xvg->meta_data = str_copy(meta_data, alloc);
	} else {
		xvg->meta_data = (str_t){0};
	}

	if (values) {
		MEMCPY(xvg->values, values, num_bytes);
	} else {
		MEMSET(xvg->values, 0, num_bytes);
	}

	return true;
}

void md_xvg_free(md_xvg_t* xvg, md_allocator_i* alloc) {
	if (md_xvg_valid(xvg)) {
		if (!str_empty(xvg->meta_data)) {
			str_free(xvg->meta_data, alloc);
		}
		const int64_t num_bytes = xvg->num_cols * xvg->num_rows * sizeof(float);
		md_free(alloc, xvg->values, num_bytes);
	}
}
