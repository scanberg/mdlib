#pragma once

#include <core/md_str.h>

struct md_allocator_i;

// Utils for xmgrace multicolum files (XVG)
// 
// Num rows dictates how many rows there are.
// Num cols dictates how many columns there are.
// The values are stored in a contigous array in row major format.
// This means that each row is stored linearly in memory followed by next row.
// meta_data is the data encoded in the @ section, it is unfiltered
typedef struct md_xvg_t {
	int num_rows;
	int num_cols;
	float* values;
	str_t meta_data;
} md_xvg_t;

#ifdef __cplusplus
extern "C" {
#endif

float* md_xvg_row(md_xvg_t* xvg, int row_index);
float  md_xvg_value(md_xvg_t* xvg, int col_index, int row_index);

str_t md_xvg_serialize_to_str(const md_xvg_t* xvg, struct md_allocator_i* str_alloc);
bool  md_xvg_deserialize_from_str(md_xvg_t* xvg, str_t data, struct md_allocator_i* alloc);
	
bool  md_xvg_serialize_to_file(const md_xvg_t* xvg, str_t path_to_file);
bool  md_xvg_init_from_file(md_xvg_t* xvg, str_t path_to_file, struct md_allocator_i* alloc);

bool  md_xvg_valid(const md_xvg_t* xvg);

// Initialize the memory of an xvg struct allocating the sufficient space required.
// values are optional and if it is supplied it is interpreted as row major and expected to have the
// same length as num_cols * num_rows. If values are NULL, then the memory will be allocated and set to zero.
bool  md_xvg_init(md_xvg_t* xvg, int num_cols, int num_rows, const float* values, str_t meta_data, struct md_allocator_i* alloc);
void  md_xvg_free(md_xvg_t* xvg, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
