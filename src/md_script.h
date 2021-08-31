#ifndef _MD_SCRIPT_H_
#define _MD_SCRIPT_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_bitfield.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule_t;
struct md_trajectory_i;
struct md_allocator_i;
struct mat4_t;

typedef enum md_script_unit_t {
    MD_SCRIPT_UNIT_NONE,
    MD_SCRIPT_UNIT_ANGSTROM,
    MD_SCRIPT_UNIT_DEGREES,
} md_script_unit_t;

typedef enum md_script_property_type_t {
    MD_SCRIPT_PROPERTY_TYPE_NONE,
    MD_SCRIPT_PROPERTY_TYPE_TEMPORAL,
    MD_SCRIPT_PROPERTY_TYPE_DISTRIBUTION,
    MD_SCRIPT_PROPERTY_TYPE_VOLUME
} md_script_property_type_t;

struct md_script_ir_o;
struct md_script_eval_o;

typedef struct md_script_error_t {
    // Data for indicating where the error occured within the script string
    int32_t line;      // Line number
    int32_t col_beg;    // String pointer
    int32_t length;    // Length in characters

    // Human readable error message
    str_t error;
} md_script_error_t;

// Opaque immediate representation (compilation result)
typedef struct md_script_ir_t {
    struct md_script_ir_o* o; // opaque internal data

    int64_t num_errors;
    const md_script_error_t* errors;
} md_script_ir_t;

typedef struct md_script_property_data_t {
    int64_t num_values;
    float*  values;     // if the underlying type is an aggregation, this will be the mean
    float*  variance;   // optional, only computed if the values are computed as an aggregate
    
    float   min_value;
    float   max_value;

    float   min_range[3];
    float   max_range[3];
} md_script_property_data_t;

// If we have a volume computed from a spatial distribution function, this will hold the metadata
// Which is required to render the reference structures for the volume.
typedef struct md_script_sdf_metadata_t {
    struct {
        int64_t count;
        md_exp_bitfield_t* atoms;
        struct mat4_t* model_to_volume_matrices;
    } reference_structures;
    md_exp_bitfield_t* target_atoms;
} md_script_sdf_metadata_t;

typedef struct md_script_property_t {
    str_t ident;
    md_script_property_type_t type;
    md_script_unit_t unit;
    int32_t dim[4];         // This gives the source data dimension
    md_script_property_data_t data;
    md_script_sdf_metadata_t* sdf_meta; // This is only provided if the property involves a sdf computation
} md_script_property_t;

// Opaque container for the evaluation result
typedef struct md_script_eval_result_t {
    struct md_script_eval_o* o; // opaque internal data

    int64_t num_properties;
    md_script_property_t* properties;
} md_script_eval_result_t;

// ### COMPILE ###
typedef struct md_script_ir_compile_args_t {
    str_t src;
    const struct md_molecule_t* mol;
    struct md_allocator_i* alloc;
} md_script_ir_compile_args_t;

bool md_script_ir_compile(md_script_ir_t* ir, md_script_ir_compile_args_t args);
bool md_script_ir_free(md_script_ir_t* ir);

// ### EVALUATE ###
typedef struct md_script_eval_args_t {
    const struct md_script_ir_t* ir;
    const struct md_molecule_t* mol;
    const struct md_trajectory_i* traj;

    // Optional, set this mask to mask out which frames should be evaluated.
    md_bitfield_t* frame_mask;

    struct md_allocator_i* alloc;
} md_script_eval_args_t;

bool md_script_eval(md_script_eval_result_t* result, md_script_eval_args_t args);
bool md_script_eval_free(md_script_eval_result_t* result);

// ### VISUALIZE ###

typedef struct md_script_token_t {
    const void* o;

    int32_t line;
    int32_t col_beg;
    int32_t col_end;
    int32_t depth;

    str_t text;
} md_script_token_t;

struct md_script_tokens_o;
typedef struct md_script_tokens_t {
    struct md_script_tokens_o* o;

    int64_t num_tokens;
    md_script_token_t* tokens;
} md_script_tokens_t;

struct md_script_visualization_o;
typedef struct md_script_visualization_t {
    struct md_script_visualization_o* o;

    struct {
        int64_t count;
        float* pos;         // Position xyz
    } vertex;

    struct {
        int64_t count;      // number of points
        uint16_t* idx;      // vertex indices;
    } point;

    struct {
        int64_t count;      // number of lines
        uint16_t* idx;      // p0 p1 vertex indices
    } line;

    struct {
        int64_t count;      // number of triangles
        uint16_t* idx;      // p0 p1 p2 vertex indices
    } triangle;

    struct {
        int64_t count;
        float* pos_rad; // xyzr xyzr
    } sphere;

    md_exp_bitfield_t atom_mask;
} md_script_visualization_t;

bool md_script_tokens_init(md_script_tokens_t* tokens, const struct md_script_ir_t* ir, struct md_allocator_i* alloc);
bool md_script_tokens_free(md_script_tokens_t* tokens);

typedef struct md_script_visualization_args_t {
    const struct md_script_token_t* token;
    const struct md_script_ir_t* ir;
    const struct md_molecule_t* mol;
    struct md_allocator_i* alloc;
} md_script_visualization_args_t;

bool md_script_visualization_init(md_script_visualization_t* vis, struct md_script_visualization_args_t args);
bool md_script_visualization_free(md_script_visualization_t* vis);

#ifdef __cplusplus
}
#endif

#endif