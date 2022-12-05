#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_unit.h>

struct md_bitfield_t;
struct md_molecule_t;
struct md_trajectory_i;
struct md_allocator_i;

/*

This is currently a giant TURD waiting to be polished into something more user friendly.
The API is a mess which has not been cleaned up for quite some time and more or less reflects all the quirks and 'Hacks'
in order to make it work within VIAMD
@TODO: Clean this shit up!

*/


// This is represents an opaque marker of something which can be visualized
//struct md_script_vis_token_t;
typedef struct md_script_vis_payload_t md_script_vis_payload_t;

typedef enum md_script_property_flags_t {
    MD_SCRIPT_PROPERTY_FLAG_TEMPORAL        = 0x0001,
    MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION    = 0x0002,
    MD_SCRIPT_PROPERTY_FLAG_VOLUME          = 0x0004,
    MD_SCRIPT_PROPERTY_FLAG_PERIODIC        = 0x0008,
    MD_SCRIPT_PROPERTY_FLAG_SDF             = 0x0010,
} md_script_property_flags_t;


// This represents the byte offset into the script source
typedef struct md_script_range_marker_t {
    int beg;
    int end;

} md_script_range_marker_t;

typedef struct md_script_vis_token_t {
    md_script_range_marker_t range;
    int depth;
    str_t text;
    const md_script_vis_payload_t* payload;
} md_script_vis_token_t;

typedef struct md_script_error_t {
    md_script_range_marker_t range;
    str_t text;
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

typedef struct md_script_property_data_t {
    int32_t dim[4];     // Dimension of values, they are exposed packed in a linear array

    int64_t num_values; // Raw 1D length of values
    float*  values;     // Raw linear access to values, check dim for the dimensions of the data

    md_script_aggregate_t* aggregate; // optional, only computed if the values are computed as an aggregate

    float   min_value;      // min total value
    float   max_value;      // max total value

    float   min_range[4];   // min range in each dimension
    float   max_range[4];   // max range in each dimension

    md_unit_t unit;

    uint64_t fingerprint; // Unique ID to compare against to see if your version is up to date.
} md_script_property_data_t;

typedef struct md_script_property_t {
    str_t ident;
    
    md_script_property_flags_t flags;
    md_script_property_data_t  data;

    const struct md_script_vis_payload_t* vis_payload; // For visualization of the property
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
        uint32_t* idx;      // vertex indices;
    } point;

    struct {
        int64_t count;      // number of lines
        uint32_t* idx;      // p0 p1 vertex indices
    } line;

    struct {
        int64_t count;      // number of triangles
        uint32_t* idx;      // p0 p1 p2 vertex indices
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
    const struct md_script_vis_payload_t* payload;
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

md_script_ir_t* md_script_ir_create(struct md_allocator_i* alloc);
void md_script_ir_free(md_script_ir_t* ir);

void md_script_ir_clear(md_script_ir_t* ir);

typedef struct md_script_bitfield_identifier_t {
    str_t identifier_name;
    const struct md_bitfield_t* bitfield;
} md_script_bitfield_identifier_t;
bool md_script_ir_add_bitfield_identifiers(md_script_ir_t* ir, const md_script_bitfield_identifier_t* bitfield_identifiers, int64_t count);

bool md_script_ir_compile_source(md_script_ir_t* ir, str_t src, const struct md_molecule_t* mol, const md_script_ir_t* ctx_ir);

int64_t md_script_ir_num_errors(const md_script_ir_t* ir);
const md_script_error_t* md_script_ir_errors(const md_script_ir_t* ir);

int64_t md_script_ir_num_vis_tokens(const md_script_ir_t* ir);
const md_script_vis_token_t* md_script_ir_vis_tokens(const md_script_ir_t* ir);

bool md_script_ir_valid(const md_script_ir_t* ir);

uint64_t md_script_ir_fingerprint(const md_script_ir_t* ir);

// ### EVALUATE ###
// This API is a dumpster-fire currently and should REALLY be simplified.

// Allocate and initialize the data for properties within the evaluation
// We need to pass the number of frames we want the data to hold
// Should be performed as soon as the IR has changed.
md_script_eval_t* md_script_eval_create(int64_t num_frames, const md_script_ir_t* ir, struct md_allocator_i* alloc);

void md_script_eval_free(md_script_eval_t* eval);

// Clear before evaluating (computing data) frames
void md_script_eval_clear(md_script_eval_t* eval);

// Compute properties
// Must be performed after the eval_init
// eval             : evaluation object to hold result
// ir               : holds an IR of the script to be evaluated
// mol              : molecule
// traj             : trajectory
// frame_(beg/end)  : range of frames [beg,end[ to evaluate 
bool md_script_eval_frame_range(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_molecule_t* mol, const struct md_trajectory_i* traj, uint32_t frame_beg, uint32_t frame_end);

// This is perhaps the granularity to operate on in a threaded context, in order to be able to interrupt evaluations without too much wait time.
//bool md_script_eval_frame(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_molecule_t* mol, const struct md_trajectory_i* traj, uint32_t frame_idx);

int64_t md_script_eval_num_properties(const md_script_eval_t* eval);
const md_script_property_t* md_script_eval_properties(const md_script_eval_t* eval);

// Compile and evaluate a single property from a string
//bool md_script_compile_and_eval_property(md_script_property_t* prop, str_t expr, const struct md_molecule_t* mol, struct md_allocator_i* alloc, const md_script_ir_t* ctx_ir, char* err_str, int64_t err_cap);

//uint32_t md_script_eval_num_frames_completed(const md_script_eval_t* eval);
//uint32_t md_script_eval_num_frames_total(const md_script_eval_t* eval);
const struct md_bitfield_t* md_script_eval_completed_frames(const md_script_eval_t* eval);
void     md_script_eval_interrupt(md_script_eval_t* eval);

// ### VISUALIZE ###
bool md_script_visualization_init(md_script_visualization_t* vis, struct md_script_visualization_args_t args);
bool md_script_visualization_free(md_script_visualization_t* vis);

// ### MISC ###

bool md_script_identifier_name_valid(str_t ident);

#ifdef __cplusplus
}
#endif
