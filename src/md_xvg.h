#pragma once

#include <core/md_str.h>
#include <core/md_array.h>

struct md_allocator_i;

// Utils for xmgrace multicolum files (XVG)

typedef struct md_xvg_header_info_t {
	str_t header;	// Complete header
	str_t comment;
	str_t meta;
	str_t title;
	str_t xaxis_label;
	str_t yaxis_label;

	int64_t num_legends;
	md_array(str_t) legends;
} md_xvg_header_info_t;

typedef struct md_xvg_t {
	md_xvg_header_info_t header_info;
	int64_t num_fields;
	int64_t num_values;
	md_array(md_array(float)) fields;
} md_xvg_t;

#ifdef __cplusplus
extern "C" {
#endif

bool md_xvg_parse_str (md_xvg_t* xvg, str_t str, struct md_allocator_i* alloc);
bool md_xvg_parse_file(md_xvg_t* xvg, str_t path_to_file, struct md_allocator_i* alloc);

void md_xvg_free(md_xvg_t* xvg, struct md_allocator_i* alloc);

// Format a header string for xvg files. cap signifies the maximum number of characters that can be written to buf.
// Returns the number of characters written to buf.
str_t	md_xvg_format_header(str_t title, str_t xaxis_label, str_t yaxis_label, int64_t num_legends, const str_t* legends, struct md_allocator_i* str_alloc);
str_t	md_xvg_format		(str_t header, int64_t num_fields, int64_t num_values, const float* const field_values[], struct md_allocator_i* str_alloc);

str_t	md_xvg_to_str(const md_xvg_t* xvg, struct md_allocator_i* str_alloc);

#ifdef __cplusplus
}
#endif
