#pragma once

#include <core/md_str.h>

struct md_allocator_i;

// Utils for comma separated value files (CSV)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_csv_t {
    int64_t num_fields;
    int64_t num_values;
    
    str_t*  field_names;    // optional, if not NULL, then should have length num_fields
    float** field_values;   // length num_fields
} md_csv_t;

// The result is an array of fields each containing the values
bool md_csv_parse_str (md_csv_t* csv, str_t in_str,  struct md_allocator_i* alloc);
bool md_csv_parse_file(md_csv_t* csv, str_t in_path, struct md_allocator_i* alloc);

str_t md_csv_write_to_str (const float* fields[], const str_t field_names[], int64_t num_fields, int64_t num_values, struct md_allocator_i* str_alloc);
bool  md_csv_write_to_file(const float* fields[], const str_t field_names[], int64_t num_fields, int64_t num_values, str_t path);

// Free csv structures created from the csv_parse functions
void md_csv_free(md_csv_t* csv, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
