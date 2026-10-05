#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_unit.h>
#include <core/md_bitfield.h>
#include <core/md_array.h>

struct md_bitfield_t;
struct md_attributes_t;
struct md_system_t;
struct md_allocator_i;

/*

This is currently a giant TURD waiting to be polished into something more user friendly.
The API is a mess which has not been cleaned up for quite some time and more or less reflects all the quirks and 'Hacks'
in order to make it work within VIAMD
@TODO: Clean this shit up!

*/

// A part of a compiled script which can be visualized: a property, or a subexpression the editor marks.
// It names the part within the one compilation which made it, and is resolved against an IR when it is used.
// Against any other IR, also one compiled later from the same source, or once its own has been recompiled,
// it resolves to nothing rather than to whatever took its place. It holds no pointer: it can be kept for as
// long as one likes, copied and passed between threads. A zeroed one refers to nothing.
typedef struct md_script_vis_ref_t {
    uint64_t ir_id;     // the compilation which made it, 0 for none
    uint32_t node_idx;
} md_script_vis_ref_t;

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
    md_script_vis_ref_t ref;
} md_script_vis_token_t;

typedef struct md_log_token_t {
    md_script_range_marker_t range;
    str_t text;
    md_bitfield_t* context;
} md_log_token_t;

// Opaque Immediate Representation (compilation result)
// Once compiled, an IR is read-only: no evaluation or visualization writes to it, so any number of them can
// share one IR, on any threads. It has to outlive them, and must not be compiled into, cleared or freed while
// they run: that is up to its owner.
typedef struct md_script_ir_t md_script_ir_t;

// Opaque object for evaluation result
typedef struct md_script_eval_t md_script_eval_t;

typedef struct md_script_vis_vertex_t {
    vec3_t pos;
	uint32_t color;
    vec3_t normal;
    uint32_t picking_idx;
} md_script_vis_vertex_t;

typedef struct md_script_vis_sphere_t {
    vec3_t pos;
    float radius;
    uint32_t color;
} md_script_vis_sphere_t;

// A label placed in the scene. 'str' is always filled in and can be drawn as it stands.
// When the label is a QUANTITY, 'value' and 'unit' carry it in numeric form as well, so a viewer
// which shows values in units of its own choosing can format it rather than take the default
// spelling in 'str'. A label which is not a quantity leaves 'unit' as none.
typedef struct md_script_vis_text_t {
    vec3_t pos;
    str_t  str;
    double value;
    md_unit_t unit;
} md_script_vis_text_t;

typedef struct md_script_vis_t {
    uint64_t magic;
    struct md_allocator_i* alloc;

    md_array(md_script_vis_vertex_t) points;
    md_array(md_script_vis_vertex_t) lines;
    md_array(md_script_vis_vertex_t) triangles;
    md_array(md_script_vis_sphere_t) spheres;
    md_array(md_script_vis_text_t)   text;

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
    MD_SCRIPT_VISUALIZE_TEXT        = 8,
};

typedef uint32_t md_script_vis_flags_t;

typedef struct md_script_vis_ctx_t {
    const struct md_script_ir_t* ir;
    const struct md_system_t* sys;
    const struct md_system_state_t* state;
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

bool md_script_ir_compile_from_source(md_script_ir_t* ir, str_t src, const struct md_system_t* sys, const md_script_ir_t* ctx_ir);

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

// If the name does not exist, it will return 0
md_script_property_flags_t md_script_ir_property_flags(const md_script_ir_t* ir, str_t name);

// ### VISUALIZATION REFERENCES ###
// A property is found by name, which is what carries over from one compilation to the next: the reference
// made from it belongs to this compilation of ir. It refers to nothing if ir has no such property.
md_script_vis_ref_t md_script_ir_property_vis_ref(const md_script_ir_t* ir, str_t name);

// Whether ref was made by this compilation of ir, and so resolves in it
bool md_script_vis_ref_valid(const md_script_ir_t* ir, md_script_vis_ref_t ref);

static inline bool md_script_vis_ref_empty(md_script_vis_ref_t ref) {
    return ref.ir_id == 0;
}

static inline bool md_script_vis_ref_equal(md_script_vis_ref_t a, md_script_vis_ref_t b) {
    return a.ir_id == b.ir_id && a.node_idx == b.node_idx;
}

// The identifier of what ref refers to in ir
// An empty string if it has none, or ref does not resolve in ir
str_t md_script_vis_ref_ident(const md_script_ir_t* ir, md_script_vis_ref_t ref);

// The major dimension (dim[0]) of what ref refers to in ir
// This can be -1 if the dimension is not determined (i.e. dynamic length)
// Zero if ref does not resolve in ir
int md_script_vis_ref_dim(const md_script_ir_t* ir, md_script_vis_ref_t ref);

// ### EVALUATE ###
// This API is a dumpster-fire currently and should REALLY be simplified.

// One evaluation holds several FRAME SETS: each is an evaluation over its own frames, and the frames are
// evaluated once, over the union of them all.
//   Temporal properties are per frame, so they are shared: evaluated once for every frame of any set, kept
//   once (in the table of set 0), and taken over a set's frames by masking (md_script_eval_frame_set_completed).
//   Distributions and volumes are aggregated over frames, so each set has its own (in its own table).
// A frame is evaluated while some set has it and has not added it yet. Interrupted, nothing is lost: the
// frames not yet added remain to be evaluated.
typedef struct md_script_eval_desc_t {
    size_t num_frames;                          // frames the evaluation can hold
    size_t num_frame_sets;                      // 0 is one set of every frame
    const struct md_bitfield_t* const* frame_sets;   // num_frame_sets of them, NULL for every frame
} md_script_eval_desc_t;

// Allocate and initialize the data for properties within the evaluation
// Should be performed as soon as the IR has changed.
md_script_eval_t* md_script_eval_create_desc(const md_script_ir_t* ir, const md_script_eval_desc_t* desc, struct md_allocator_i* alloc);

// One frame set of every frame
md_script_eval_t* md_script_eval_create(size_t num_frames, const md_script_ir_t* ir, struct md_allocator_i* alloc);

// ### FRAME SETS ###
// Not while frames are being evaluated: interrupt and wait for the evaluation first.
size_t md_script_eval_frame_set_count(const md_script_eval_t* eval);

// Replaces the frames of set idx (NULL: every frame) and clears what it has aggregated. Temporal values
// already evaluated are kept: only the frames the set has not added remain to be evaluated.
bool md_script_eval_set_frame_set(md_script_eval_t* eval, size_t idx, const struct md_bitfield_t* frames);

const struct md_bitfield_t* md_script_eval_frame_set_frames(const md_script_eval_t* eval, size_t idx);

// The frames set idx has added so far. The mask to take temporal properties over it with: each was
// evaluated (md_script_eval_frame_mask) when it was added.
const struct md_bitfield_t* md_script_eval_frame_set_completed(const md_script_eval_t* eval, size_t idx);

// The table of set idx: for set 0 the same as md_script_eval_attributes, temporal properties included.
// Any other set publishes only its distributions and volumes, under the same paths.
const struct md_attributes_t* md_script_eval_frame_set_attributes(const md_script_eval_t* eval, size_t idx);

// The frames which remain to be evaluated, written to out. Returns their number.
size_t md_script_eval_pending_frames(md_script_eval_t* eval, struct md_bitfield_t* out);

uint64_t md_script_eval_ir_fingerprint(const md_script_eval_t* eval);

void md_script_eval_free(md_script_eval_t* eval);

// Clear before evaluating (computing data) frames
void md_script_eval_clear_data(md_script_eval_t* eval);

// Compute properties
// eval             : evaluation object to hold result
// ir               : holds an IR of the script to be evaluated
// sys              : the system, whose attribute table holds the run
// run              : "run/<name>", the run whose frames are evaluated (see RUNS in md_system.h)
// frames           : the frames to evaluate. Those with nothing left to do (see FRAME SETS) are skipped.
// Safe to call from several threads at once, also with overlapping frames (a frame is added to a set once):
// each call keeps its own extraction context, and with it its own open files.
bool md_script_eval_frames(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_system_t* sys, str_t run, const uint32_t* frames, size_t num_frames);

// The same, for the frames [beg,end[
bool md_script_eval_frame_range(md_script_eval_t* eval, const struct md_script_ir_t* ir, const struct md_system_t* sys, str_t run, uint32_t frame_beg, uint32_t frame_end);

// Number of properties held by the evaluation
size_t md_script_eval_property_count(const md_script_eval_t* eval);

// The evaluated properties, as ATTRIBUTES published under 'script/'. This is the only form the
// results take: the evaluation writes straight into this table and keeps no other copy.
//
//   script/<ident>            the values. Temporal properties carry MD_ATTRIBUTE_FLAG_TEMPORAL and
//                             are shaped {num_frames, population}; a distribution is {num_bins};
//                             a volume is {x, y, z}. The unit is the unit of the values.
//   script/<ident>/range      rank 0, 2 components (min, max). The range the script declared for
//                             the property, or the observed one when it declared none. Temporal:
//                             the range of the values. Distribution: the range of the bin axis.
//   script/<ident>/mean       temporal with a population only: the per frame summary over it
//   script/<ident>/variance
//   script/<ident>/extent     2 components, (min, max)
//   script/<ident>/weight     a distribution's per bin weight
//   script/<ident>/bin        a distribution's bin coordinates (virtual), in the bin axis unit
//   time                      the frame axis every temporal property is temporal along (see FRAME
//                             AXES in md_system.h). Frame ordinals without a unit: the evaluation
//                             knows how many frames it covers, not when they were taken.
//
// The table and the attributes in it are created with the evaluation and never added to or
// removed, so attribute pointers stay valid for the lifetime of the evaluation. The data is written
// in place while frames are evaluated; versions bump when md_script_eval_frame_range completes, so
// a consumer caching something derived from a property compares md_attributes_version.
const struct md_attributes_t* md_script_eval_attributes(const md_script_eval_t* eval);

// The frames whose temporal values are evaluated, of any set
size_t md_script_eval_frame_count(const md_script_eval_t* eval);
const struct md_bitfield_t* md_script_eval_frame_mask(const md_script_eval_t* eval);

// Interrupt the current evaluation. It stays interrupted until reset: reset before evaluating again.
void md_script_eval_interrupt(md_script_eval_t* eval);
void md_script_eval_reset_interrupt(md_script_eval_t* eval);

// ### VISUALIZE ###
void md_script_vis_init(md_script_vis_t* vis, struct md_allocator_i* alloc);
bool md_script_vis_free(md_script_vis_t* vis);

bool md_script_vis_clear(md_script_vis_t* vis);
// Visualizes what ref refers to, resolved against ctx->ir. False if it does not resolve there.
// subidx selects one element of an array (-1 for all of them)
bool md_script_vis_eval_ref(md_script_vis_t* vis, md_script_vis_ref_t ref, int subidx, const md_script_vis_ctx_t* ctx, md_script_vis_flags_t flags);

bool md_script_vis_eval_string(md_script_vis_t* vis, str_t str, const md_script_vis_ctx_t* ctx, md_script_vis_flags_t flags);

// ### MISC ###

bool md_script_identifier_name_valid(str_t ident);

// Reserved words of the language (and, or, xor, not, in, of, out)
size_t       md_script_num_keywords(void);
const str_t* md_script_keywords(void);

// Names the language defines: every built-in procedure (once, however many overloads it has) and the predefined
// constants. Meant for syntax highlighting and completion.
// Writes up to cap names to out and returns how many there are, so md_script_builtin_identifiers(NULL, 0) gives the
// size of the buffer to pass.
size_t md_script_builtin_identifiers(str_t* out, size_t cap);

// ### COMPLETION ###

typedef enum md_script_completion_kind_t {
    MD_SCRIPT_COMPLETION_KEYWORD = 1,
    MD_SCRIPT_COMPLETION_PROCEDURE,     // A built-in procedure
    MD_SCRIPT_COMPLETION_CONSTANT,      // A predefined constant, such as PI
    MD_SCRIPT_COMPLETION_VARIABLE,      // A name the script assigns to
    MD_SCRIPT_COMPLETION_PARAMETER,     // A parameter of the procedure being called, to give the argument by name
    MD_SCRIPT_COMPLETION_VALUE,         // A value the argument takes: a name found in the system (residue names,
                                        // elements, chains, ...) or one of a fixed set, such as the units of count()
} md_script_completion_kind_t;

typedef struct md_script_completion_t {
    str_t label;                        // What to show and to match what has been typed against: ALA, cutoff
    str_t text;                         // What replaces the range: "ALA" (quoted if the cursor is not in a string),
                                        // cutoff=
    md_script_completion_kind_t kind;
} md_script_completion_t;

typedef struct md_script_completions_t {
    md_script_range_marker_t range;     // The bytes of the source that a completion replaces: the identifier, or the
                                        // contents of the string, that the cursor is in. Empty at the cursor if none.
    str_t prefix;                       // The part of the range before the cursor: what has been typed so far
    size_t count;
    md_script_completion_t* items;
} md_script_completions_t;

// What can be written at the byte offset cursor of the script src, for completion in an editor.
// The source may be incomplete, as it is while being typed: this follows the tokens up to the cursor and does not
// compile anything. The items are all that is valid at the cursor, neither filtered by what has been typed nor
// ranked, which is left to the editor. Nothing is offered in comments and numbers, and inside a string only when the
// argument it is given to takes known values. Values from the system are only offered if sys is not NULL.
// Everything returned (items and strings) is allocated from alloc, which is meant to be a temporary arena.
md_script_completions_t md_script_complete(str_t src, int cursor, const struct md_system_t* sys, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
