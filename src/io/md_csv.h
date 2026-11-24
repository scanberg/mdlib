#pragma once

#include <core/md_str.h>
#include <md_io.h>

struct md_allocator_i;
struct md_system_t;
struct md_system_state_t;
struct md_load_desc_t;
struct md_writer_i;

// Utils for comma separated value files (CSV)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_csv_t {
    size_t num_fields;
    size_t num_values;
    
    md_str_t*  field_names;    // optional, if not NULL, then should have length num_fields
    float** field_values;   // length num_fields
} md_csv_t;

// The result is an array of fields each containing the values
bool md_csv_parse_str (md_csv_t* csv, md_str_t in_str,  struct md_allocator_i* alloc);

// Writes the fields as columns, a row per value, with field_names (if not NULL) as the first row.
// Values are written as "%f" writes them.
void  md_csv_write(struct md_writer_i* w, const float* const field_values[], const md_str_t field_names[], size_t num_fields, size_t num_values);

// The same, as a string allocated from str_alloc (see md_alloc_str), and to a file which is created or truncated
md_str_t md_csv_write_to_str (const float* const field_values[], const md_str_t field_names[], size_t num_fields, size_t num_values, struct md_allocator_i* str_alloc);
bool  md_csv_write_to_file(const float* const field_values[], const md_str_t field_names[], size_t num_fields, size_t num_values, md_str_t path_to_file);

// Free csv structures created from the csv_parse functions
void md_csv_free(md_csv_t* csv, struct md_allocator_i* alloc);

// Its columns as a series along the run desc->run, "run/<name>" (see md_run_publish_series in
// md_system.h), in the group "csv/<file stem>": with a first column named time, its values are the
// group's own time axis and every run frame must find its time there; otherwise the rows are the
// run's frames. A column's unit is the one in parentheses in its name, "Distance (nm)". Read in the
// script with attr("csv/<stem>/<column>").
// Data of a system's atoms (md_io.h, LOAD): onto a system with the run, MD_LOAD_MISMATCH when the file is not along it.
md_load_result_t md_csv_system_load_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t filename, const struct md_load_desc_t* desc);

// Whether the file is along the run desc->run of the system, changing nothing (md_io.h, CHECK)
md_load_result_t md_csv_system_check_file(const struct md_system_t* sys, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif
