#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>

struct md_bitfield_t;
struct md_molecule_t;
struct md_trajectory_i;
struct md_allocator_i;

// This is represents an opaque token of something which can be visualized
struct md_script_vis_token_t;

typedef enum md_script_property_flags_t {
    MD_SCRIPT_PROPERTY_FLAG_TEMPORAL        = 0x0001,
    MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION    = 0x0002,
    MD_SCRIPT_PROPERTY_FLAG_VOLUME          = 0x0004,
    MD_SCRIPT_PROPERTY_FLAG_PERIODIC        = 0x0008,
    MD_SCRIPT_PROPERTY_FLAG_SDF             = 0x0010,
} md_script_property_flags_t;

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
    int32_t col_beg;   // Column offset
    int32_t length;    // Length in characters

    // Human readable error message
    str_t error;
} md_script_error_t;

// Opaque Immediate Representation (compilation result)
typedef struct md_script_ir_t md_script_ir_t;

// Opaque object for evaluation result
typedef struct md_script_eval_t md_script_eval_t;

typedef struct md_script_aggregate_t {
    int64_t num_values;
    float*  population_mean;
    float*  population_var;
    float*  population_min;
    float*  population_max;
} md_script_aggregate_t;

/*
typedef struct md_script_distribution_t {
    uint32_t num_bins;
} md_script_distribution_t;
*/

typedef struct md_script_property_data_t {
    int32_t dim[4];     // Dimension of values, they are exposed packed in a linear array

    int64_t num_values; // Raw 1D length of values
    float*  values;     // Raw linear access to values, check dim for the dimensions of the data

    md_script_aggregate_t* aggregate; // optional, only computed if the values are computed as an aggregate

    float   min_value;      // min total value
    float   max_value;      // max total value

    float   min_range[4];   // min range in each dimension
    float   max_range[4];   // max range in each dimension

    char    unit[32];

    uint64_t fingerprint; // Unique ID to compare against to see if your version is up to date.
} md_script_property_data_t;

typedef struct md_script_property_t {
    str_t ident;
    
    md_script_property_flags_t flags;
    md_script_property_data_t  data;

    const struct md_script_vis_token_t* vis_token; // For visualization of the property
} md_script_property_t;

struct md_script_visualization_o;
typedef struct md_script_visualization_t {
    struct md_script_visualization_o* o;

    // Geometry
    struct {
        int64_t count;
        vec3_t* pos;        // Position xyz
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
        vec4_t* pos_rad;     // xyzr
    } sphere;

    // This is a bit of a shoe-horn case where we want to visualize the superimposed structures and the atoms involved
    // in computing an SDF, therefore this requires transformation matrices as well as the involved structures
    struct {
        int64_t count;
        mat4_t* matrices;
        struct md_bitfield_t* structures;
        float extent;
    } sdf;

    // Atoms
    struct {
        struct md_bitfield_t* atom_masks;
    } structures;

    struct md_bitfield_t* atom_mask;
    
} md_script_visualization_t;

enum {
    MD_SCRIPT_VISUALIZE_DEFAULT     = 0, // Default is to visualize everything, equivalent to all flags
    MD_SCRIPT_VISUALIZE_GEOMETRY    = 1,
    MD_SCRIPT_VISUALIZE_ATOMS       = 2,
    MD_SCRIPT_VISUALIZE_SDF         = 4,
};

typedef uint32_t md_script_visualization_flags_t;

typedef struct md_script_visualization_args_t {
    const struct md_script_vis_token_t* token;
    const struct md_script_ir_t* ir;
    const struct md_molecule_t* mol;
    const struct md_trajectory_i* traj;
    struct md_allocator_i* alloc;
    md_script_visualization_flags_t flags;
} md_script_visualization_args_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### IR ###
// src      : source code to compile
// mol      : molecule
// alloc    : allocator
// ctx_ir   : provide a context of identifiers and expressions [optional]
md_script_ir_t* md_script_ir_create(str_t src, const struct md_molecule_t* mol, struct md_allocator_i* alloc, const md_script_ir_t* ctx_ir);
void md_script_ir_free(md_script_ir_t* ir);

int64_t md_script_ir_num_errors(const md_script_ir_t* ir);
const md_script_error_t* md_script_ir_errors(const md_script_ir_t* ir);

int64_t md_script_ir_num_tokens(const md_script_ir_t* ir);
const md_script_token_t* md_script_ir_tokens(const md_script_ir_t* ir);

bool md_script_ir_valid(const md_script_ir_t* ir);

// ### EVALUATE ###

// Allocate and initialize the data for properties within the evaluation
// We need to pass the number of frames we want the data to hold
// Should be performed as soon as the IR has changed.
md_script_eval_t* md_script_eval_create(int64_t num_frames, const md_script_ir_t* ir, struct md_allocator_i* alloc);

void md_script_eval_free(md_script_eval_t* eval);

// Compute properties
// Must be performed after the eval_init
// eval     : evaluation object to hold result
// ir       : holds an IR of the script to be evaluated
// mol      : molecule
// traj     : trajectory
// frame_mask [OPTIONAL] : A mask which holds the frames which should be evaluated. Supply this if only a subset should be evaluated.
bool md_script_eval_compute(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_molecule_t* mol, const struct md_trajectory_i* traj, const struct md_bitfield_t* frame_mask);

int64_t md_script_eval_num_properties(const md_script_eval_t* eval);
const md_script_property_t* md_script_eval_properties(const md_script_eval_t* eval);

// Compile and evaluate a single property from a string
//bool md_script_compile_and_eval_property(md_script_property_t* prop, str_t expr, const struct md_molecule_t* mol, struct md_allocator_i* alloc, const md_script_ir_t* ctx_ir, char* err_str, int64_t err_cap);


uint32_t md_script_eval_num_completed_frames(const md_script_eval_t* eval);
float    md_script_eval_completed_fraction(const md_script_eval_t* eval);
void     md_script_eval_interrupt(md_script_eval_t* eval);

// ### VISUALIZE ###
bool md_script_visualization_init(md_script_visualization_t* vis, struct md_script_visualization_args_t args);
bool md_script_visualization_free(md_script_visualization_t* vis);

#ifdef __cplusplus
}
#endif
