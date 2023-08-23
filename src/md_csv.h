#pragma once

#include <core/md_str.h>
#include <core/md_array.h>

struct md_allocator_i;

// Utils for comma separated value files (CSV)

#ifdef __cplusplus
extern "C" {
#endif

md_array(float) md_csv_parse_str(int* out_num_cols, int* out_num_rows, str_t in_str, struct md_allocator_i* alloc);
md_array(float) md_csv_parse_file(int* out_num_cols, int* out_num_rows, str_t in_path, struct md_allocator_i* alloc);

str_t md_csv_write_to_str (const float* data, int num_cols, int num_rows, struct md_allocator_i* str_alloc);
bool  md_csv_write_to_file(const float* data, int num_cols, int num_rows, str_t path);


#ifdef __cplusplus
}
#endif
