#ifndef _MD_SCRIPT_H_
#define _MD_SCRIPT_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <core/md_vec_math.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule_t;
struct md_trajectory_i;
struct md_allocator_i;
struct md_frame_cache_t;

// This is represents something which can be visualized
struct md_script_vis_token_t;

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

typedef struct md_script_token_t {
    const struct md_script_vis_token_t* vis_token;

    int32_t line;
    int32_t col_beg;
    int32_t col_end;
    int32_t depth;

    str_t text;
} md_script_token_t;

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
    uint64_t fingerprint; // Unique ID to compare against to see if your representation is up to date.

    int64_t num_errors;
    const md_script_error_t* errors;

    int64_t num_tokens;
    const md_script_token_t* tokens;
} md_script_ir_t;

typedef struct md_script_property_data_aggregate_t {
    int64_t num_values;
    float*  mean;       
    float*  variance;
} md_script_property_data_aggregate_t;

typedef struct md_script_property_data_t {
    int32_t dim[4];     // Dimension of values, they are exposed packed in a linear array

    int64_t num_values; // Raw 1D length of values
    float*  values;     // Raw linear access to values, check dim for the dimensions of the data

    md_script_property_data_aggregate_t* aggregate; // optional, only computed if the values are computed as an aggregate
    
    float   min_value;
    float   max_value;

    float   min_range[4];   // min range in each dimension
    float   max_range[4];   // max range in each dimension
} md_script_property_data_t;

typedef struct md_script_property_t {
    str_t ident;
    md_script_property_type_t type;
    md_script_unit_t unit;
    md_script_property_data_t data;

    const struct md_script_vis_token_t* vis_token; // For visualization of the property
} md_script_property_t;

// Opaque container for the evaluation result
typedef struct md_script_eval_t {
    struct md_script_eval_o* o; // opaque internal data
    uint64_t fingerprint; // Unique ID to compare against to see if your representation is up to date.

    int64_t num_properties;
    md_script_property_t* properties;
} md_script_eval_t;

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

    // Optional, if this is set, this will be used to load trajectory frames.
    struct md_frame_cache_t* frame_cache;

    // Optional, set this mask to mask out which frames that should be evaluated.
    md_exp_bitfield_t* filter_mask;
} md_script_eval_args_t;

// Allocate the data for properties within the evaluation
// We need to pass the number of frames we want the data to hold
// Should be performed as soon as the IR has changed.
bool md_script_eval_alloc(md_script_eval_t* eval, int64_t num_frames, const md_script_ir_t* ir, struct md_allocator_i* alloc);

// Compute properties
// Must be performed after the properties has been allocated.
bool md_script_eval_compute(md_script_eval_t* eval, md_script_eval_args_t args);

bool md_script_eval_free(md_script_eval_t* eval);

// These is meant to be used if the evaulation runs in its own thread (which is recommended)
int md_script_eval_completed_frame_count(const md_script_eval_t* eval);
void md_script_eval_interrupt(const md_script_eval_t* eval);

// ### VISUALIZE ###

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
        float* pos_rad;     // xyzr
    } sphere;

    // This is a bit of a shoe-horn case where we want to visualize the superimposed structures and the atoms involved
    // in computing an SDF, therefore this requires transformation matrices as well as the involved structures
    struct {
        int64_t count;
        mat4_t* matrices;
        md_exp_bitfield_t* structures;
        float extent;
    } sdf;

    md_exp_bitfield_t atom_mask;
} md_script_visualization_t;

typedef struct md_script_visualization_args_t {
    const struct md_script_vis_token_t* token;
    const struct md_script_ir_t* ir;
    const struct md_molecule_t* mol;
    const struct md_trajectory_i* traj;

    // Optional, if this is set, this will be used to load trajectory frames
    struct md_frame_cache_t* frame_cache;

    struct md_allocator_i* alloc;
} md_script_visualization_args_t;

bool md_script_visualization_init(md_script_visualization_t* vis, struct md_script_visualization_args_t args);
bool md_script_visualization_free(md_script_visualization_t* vis);

#ifdef __cplusplus
}
#endif

#endif