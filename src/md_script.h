#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_unit.h>
#include <core/md_bitfield.h>
#include <core/md_array.h>

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

// This is represents an opaque marker of a portion of the AST-tree that can be used to invoke visualizations
typedef struct md_script_vis_payload_o md_script_vis_payload_o;

typedef enum md_script_property_flags_t {
    MD_SCRIPT_PROPERTY_FLAG_NONE            = 0,
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
    const md_script_vis_payload_o* payload;
} md_script_vis_token_t;

typedef struct md_log_token_t {
    md_script_range_marker_t range;
    str_t text;
    md_bitfield_t* context;
} md_log_token_t;

// Opaque Immediate Representation (compilation result)
typedef struct md_script_ir_t md_script_ir_t;

// Opaque object for evaluation result
typedef struct md_script_eval_t md_script_eval_t;

typedef struct md_script_aggregate_t {
    size_t  num_values;
    float*  population_mean;
    float*  population_var;
    vec2_t* population_ext; // min and max
    //float*  population_min;
    //float*  population_max;
} md_script_aggregate_t;

typedef struct md_script_property_data_t {
    int32_t dim[4];     // Dimension of values, they are exposed packed in a linear array

    size_t num_values; // Raw 1D length of values
    float*  values;     // Raw linear access to values, check dim for the dimensions of the data
    float*  weights;    // Optional, only applicable to distributions

    md_script_aggregate_t* aggregate; // optional, only computed if the values are computed as an aggregate

    float   min_value;      // min total value
    float   max_value;      // max total value

    float   min_range[2];   // min range in each dimension, [0] for temporal, [0]+[1] for distribution
    float   max_range[2];   // max range in each dimension, [0] for temporal, [0]+[1] for distribution

    md_unit_t unit[2];      // [0] is time/frame for temporal, x-value for distributions and 3D.
                            // [1] is the unit of the data value

    uint64_t fingerprint;   // Essentially a checksum of the data to compare against
} md_script_property_data_t;
 
typedef struct md_script_property_t {
    str_t ident;
    
    md_script_property_flags_t flags;
    md_script_property_data_t  data;

    const md_script_vis_payload_o* vis_payload; // For visualization of the property
} md_script_property_t;

typedef struct md_script_vis_vertex_t {
    vec3_t pos;
	uint32_t color;
} md_script_vis_vertex_t;

typedef struct md_script_vis_sphere_t {
    vec3_t pos;
    float radius;
    uint32_t color;
} md_script_vis_sphere_t;

typedef struct md_script_vis_t {
    uint64_t magic;
    struct md_allocator_i* alloc;

    md_array(md_script_vis_vertex_t) points;
    md_array(md_script_vis_vertex_t) lines;
    md_array(md_script_vis_vertex_t) triangles;
    md_array(md_script_vis_sphere_t) spheres;

    // This is a bit of a shoe-horn case where we want to visualize the superimposed structures and the atoms involved
    // in computing an SDF, therefore this requires transformation matrices as well as the involved structures
    struct {
        md_array(mat4_t) matrices;
        md_array(struct md_bitfield_t) structures;
        float extent;
    } sdf;

    md_array(struct md_bitfield_t) structures;

    // This is an all encompassing atom mask, when there may be no real 'substructures'
    struct md_bitfield_t atom_mask;
} md_script_vis_t;

enum {
    MD_SCRIPT_VISUALIZE_DEFAULT     = 0, // Default is to visualize everything, equivalent to all flags
    MD_SCRIPT_VISUALIZE_GEOMETRY    = 1,
    MD_SCRIPT_VISUALIZE_ATOMS       = 2,
    MD_SCRIPT_VISUALIZE_SDF         = 4,
};

typedef uint32_t md_script_vis_flags_t;

typedef struct md_script_vis_ctx_t {
    const struct md_script_ir_t* ir;
    const struct md_molecule_t* mol;
    const struct md_trajectory_i* traj;
} md_script_vis_ctx_t;

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

bool md_script_ir_add_identifier_bitfield(md_script_ir_t* ir, str_t ident, const struct md_bitfield_t* bf);

// Returns true if any expression in the script contains a reference to to an identifier by supplied name
bool md_script_ir_contains_identifier_reference(const md_script_ir_t* ir, str_t name);

bool md_script_ir_compile_from_source(md_script_ir_t* ir, str_t src, const struct md_molecule_t* mol, const struct md_trajectory_i* traj, const md_script_ir_t* ctx_ir);

size_t md_script_ir_num_errors(const md_script_ir_t* ir);
const md_log_token_t* md_script_ir_errors(const md_script_ir_t* ir);

size_t md_script_ir_num_warnings(const md_script_ir_t* ir);
const md_log_token_t* md_script_ir_warnings(const md_script_ir_t* ir);

size_t md_script_ir_num_vis_tokens(const md_script_ir_t* ir);
const md_script_vis_token_t* md_script_ir_vis_tokens(const md_script_ir_t* ir);

bool md_script_ir_valid(const md_script_ir_t* ir);

uint64_t md_script_ir_fingerprint(const md_script_ir_t* ir);

// Get identifiers within a script
size_t md_script_ir_num_identifiers(const md_script_ir_t* ir);
const str_t* md_script_ir_identifiers(const md_script_ir_t* ir);

// ### PROPERTIES ###
size_t       md_script_ir_property_count(const md_script_ir_t* ir);
const str_t* md_script_ir_property_names(const md_script_ir_t* ir);

// Extract properties which has the supplied flags set.
//size_t md_script_ir_property_id_filter_on_flags(str_t* out_names, size_t out_cap, const md_script_ir_t* ir, md_script_property_flags_t flags);

// If the name does not exist, it will return 0
md_script_property_flags_t      md_script_ir_property_flags(const md_script_ir_t* ir, str_t name);
const md_script_vis_payload_o*  md_script_ir_property_vis_payload(const md_script_ir_t* ir, str_t name);

// Returns the identifier of the payload
// May return an empty string if not valid
str_t md_script_payload_ident(const md_script_vis_payload_o* payload);

// Returns the major dimension (dim[0]) of the payload
// returns zero if the payload is not valid
size_t md_script_payload_dim(const md_script_vis_payload_o* payload);

// ### EVALUATE ###
// This API is a dumpster-fire currently and should REALLY be simplified.

// Allocate and initialize the data for properties within the evaluation
// We need to pass the number of frames we want the data to hold
// Should be performed as soon as the IR has changed.
md_script_eval_t* md_script_eval_create(size_t num_frames, const md_script_ir_t* ir, struct md_allocator_i* alloc);

uint64_t md_script_eval_ir_fingerprint(const md_script_eval_t* eval);

void md_script_eval_free(md_script_eval_t* eval);

// Clear before evaluating (computing data) frames
void md_script_eval_clear_data(md_script_eval_t* eval);

// Compute properties
// Must be performed after the eval_init
// eval             : evaluation object to hold result
// ir               : holds an IR of the script to be evaluated
// mol              : molecule
// traj             : trajectory
// frame_(beg/end)  : range of frames [beg,end[ to evaluate 
bool md_script_eval_frame_range(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_molecule_t* mol, const struct md_trajectory_i* traj, uint32_t frame_beg, uint32_t frame_end);

// Extract property data within in an evaluation
size_t md_script_eval_property_count(const md_script_eval_t* eval);
const md_script_property_data_t* md_script_eval_property_data(const md_script_eval_t* eval, str_t name);

// Get the frames encoded as a bitfield of the completed frames
size_t md_script_eval_frame_count(const md_script_eval_t* eval);
const struct md_bitfield_t* md_script_eval_frame_mask(const md_script_eval_t* eval);

// Interrupt the current evaluation
void md_script_eval_interrupt(md_script_eval_t* eval);

// ### VISUALIZE ###
void md_script_vis_init(md_script_vis_t* vis, struct md_allocator_i* alloc);
bool md_script_vis_free(md_script_vis_t* vis);

bool md_script_vis_clear(md_script_vis_t* vis);
bool md_script_vis_eval_payload(md_script_vis_t* vis, const md_script_vis_payload_o* payload, int subidx, const md_script_vis_ctx_t* ctx, md_script_vis_flags_t flags);

// ### MISC ###

bool md_script_identifier_name_valid(str_t ident);

#ifdef __cplusplus
}
#endif
