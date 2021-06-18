#ifndef _MD_SCRIPT_H_
#define _MD_SCRIPT_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_bitfield.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule;
struct md_trajectory;
struct md_allocator;

typedef enum md_script_unit {
    MD_SCRIPT_UNIT_NONE,
    MD_SCRIPT_UNIT_ANGSTROM,
    MD_SCRIPT_UNIT_DEGREES,
} md_script_unit_t;

typedef enum md_script_property_type {
    MD_SCRIPT_TYPE_NONE,
    MD_SCRIPT_TYPE_TEMPORAL,
    MD_SCRIPT_TYPE_DISTRIBUTION,
    MD_SCRIPT_TYPE_VOLUME
} md_script_property_type_t;

struct md_script_ir_o;
struct md_script_eval_o;

typedef struct md_script_error {
    // Data for indicating where the error occured within the script string
    int32_t line;      // Line number
    int32_t offset;    // String pointer
    int32_t length;    // Length in characters

    // Human readable error message
    str_t error;
} md_script_error_t;

// Opaque immediate representation (compilation result)
typedef struct md_script_ir {
    struct md_script_ir_o* o; // opaque internal data

    int64_t num_errors;
    const md_script_error_t* errors;
} md_script_ir_t;

typedef struct md_script_scalar {
    int64_t num_values;
    float*  values;     // if the underlying type is an aggregation, this will be the mean
    float*  variance;   // optional, only computed if the values are computed as an aggregate
    
    float   min_value;
    float   max_value;
} md_script_scalar_t;

typedef struct md_script_distribution {
    int64_t num_values;
    float* values;
    float* variance;
    
    float min_value;
    float max_value;

    float min_range;
    float max_range;

} md_script_distribution_t;

typedef struct md_script_volume {
    // Extent of the volume
    float ext[3];

    int64_t num_values; // should be the same as dim[0] * dim[1] * dim[2]
    float* values;
    float* variance;    // Optional, only if aggregated, num_values in length
    
    float min_value;
    float max_value;

} md_script_volume_t;

typedef struct md_script_property {
    str_t ident;
    md_script_property_type_t type;
    md_script_unit_t unit;
    int32_t dim[4];         // This gives the source data dimension

    // These are optional, and one or more may be allocated depending on the type of the property.
    md_script_scalar_t*         temporal;
    md_script_distribution_t*   distribution;
    md_script_volume_t*         volume;
} md_script_property_t;

// Opaque container for the evaluation result
typedef struct md_script_eval_result {
    struct md_script_eval_o* o; // opaque internal data

    int64_t num_properties;
    const md_script_property_t* properties;
} md_script_eval_result_t;

// ### COMPILE ###
typedef struct md_script_ir_compile_args {
    str_t src;
    const struct md_molecule* mol;
    struct md_allocator* alloc;
} md_script_ir_compile_args_t;

bool md_script_ir_compile(struct md_script_ir* ir, struct md_script_ir_compile_args args);
bool md_script_ir_free(struct md_script_ir* ir);

// ### EVALUATE ###
typedef struct md_script_eval_args {
    const struct md_script_ir* ir;
    const struct md_molecule* mol;
    const struct md_trajectory* traj;

    // Optional, set to zero if all frames should be evaluated
    md_bitfield_t* frame_mask;

    struct md_allocator* alloc;
} md_script_eval_args_t;

bool md_script_eval(struct md_script_eval_result* result, struct md_script_eval_args args);
bool md_script_eval_free(struct md_script_eval_result* result);

#ifdef __cplusplus
}
#endif

#endif