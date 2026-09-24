#pragma once

#include <core/md_str.h>

struct md_allocator_i;
struct md_system_t;

// Utils for comma separated value files (CSV)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_csv_t {
    size_t num_fields;
    size_t num_values;
    
    str_t*  field_names;    // optional, if not NULL, then should have length num_fields
    float** field_values;   // length num_fields
} md_csv_t;

// The result is an array of fields each containing the values
bool md_csv_parse_str (md_csv_t* csv, str_t in_str,  struct md_allocator_i* alloc);
bool md_csv_parse_file(md_csv_t* csv, str_t in_path, struct md_allocator_i* alloc);

str_t md_csv_write_to_str (const float* const field_values[], const str_t field_names[], size_t num_fields, size_t num_values, struct md_allocator_i* str_alloc);
bool  md_csv_write_to_file(const float* const field_values[], const str_t field_names[], size_t num_fields, size_t num_values, str_t path_to_file);

// Free csv structures created from the csv_parse functions
void md_csv_free(md_csv_t* csv, struct md_allocator_i* alloc);

// The file's columns as a series along the run "run/<name>" (see md_run_publish_series in
// md_system.h), in the group "csv/<file stem>": with a first column named time, its values are the
// group's own time axis and every run frame must find its time there; otherwise the rows are the
// run's frames. A column's unit is the one in parentheses in its name, "Distance (nm)". Read in the
// script with attr("csv/<stem>/<column>").
bool md_csv_system_supplement_from_file(struct md_system_t* sys, str_t filename, str_t run);

#ifdef __cplusplus
}
#endif
