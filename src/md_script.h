#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_unit.h>
#include <core/md_bitset.h>
#include <core/md_array.h>

struct md_attribute_t;
struct md_attributes_t;
struct md_system_t;
struct md_allocator_i;
struct md_writer_i;

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
    // A temporal property of periodic values, such as angles about an axis: a value is the same as itself plus a whole
    // period. script/<ident>/period and /seam say which period, and where a plot of the values is best cut (see
    // md_script_eval_attributes), and the summary over its population is a circular one.
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
    const md_bitset_t* context;
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
        md_array(md_bitset_t) structures;   // payloads in alloc
        float extent;
    } sdf;

    md_array(md_bitset_t) structures;       // payloads in alloc

    // The atoms of each element of the bitfields visualized, in order: of residue(:) one set per residue, of a
    // single bitfield one set (payloads in alloc). For telling the elements apart, as groups of the atom mask.
    md_array(md_bitset_t) groups;

    // This is an all encompassing atom mask, when there may be no real 'substructures' (payload in alloc)
    md_bitset_t atom_mask;
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

// Adds a constant identifier holding a set of atoms. The IR keeps its own copy.
bool md_script_ir_add_identifier_bitset(md_script_ir_t* ir, str_t ident, md_bitset_t set);

// Returns true if any expression in the script contains a reference to to an identifier by supplied name
bool md_script_ir_contains_identifier_reference(const md_script_ir_t* ir, str_t name);

bool md_script_ir_compile_from_source(md_script_ir_t* ir, str_t src, const struct md_system_t* sys, const md_script_ir_t* ctx_ir);

size_t md_script_ir_num_errors(const md_script_ir_t* ir);
const md_log_token_t* md_script_ir_errors(const md_script_ir_t* ir);

size_t md_script_ir_num_warnings(const md_script_ir_t* ir);
const md_log_token_t* md_script_ir_warnings(const md_script_ir_t* ir);

// Writes the errors and warnings of the last compilation to w, in the format of GCC and clang: in the order they
// appear in the source, each with its location (line:column, 1-based, '<script>' in place of a file name) and message,
// followed by the line(s) of the source it refers to with its range marked beneath. Ends with clang's summary:
//
//   <script>:2:5: error: Undefined identifier 'y'
//       2 | b = y + 2;
//         |     ^
//   1 error generated.
//
// Writes nothing if there is nothing to report. To see why a script did not compile, in a test or a tool:
//   md_script_ir_write_diagnostics(md_writer_stdout(), ir);
void md_script_ir_write_diagnostics(struct md_writer_i* w, const md_script_ir_t* ir);

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

// What a property is computed from: its expression, with every variable it uses standing for what that variable
// is computed from in turn, the constants and selections it reads by their contents, the attributes by their
// version, and the topology of the system it was compiled against. Not its name, nor anything else in the script.
// Two compilations agree on it exactly when the property's values agree, which is what lets an evaluation keep the
// values of a property across an edit of the script that does not touch it (see md_script_eval_create_from).
// Compared within one process only. 0 if there is no such property.
uint64_t md_script_ir_property_fingerprint(const md_script_ir_t* ir, str_t name);

// Where the property is named in the source it was compiled from: the bytes of its name where it is assigned. Empty
// if there is no such property.
md_script_range_marker_t md_script_ir_property_range(const md_script_ir_t* ir, str_t name);

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

// An evaluation holds every property over the frames it has evaluated. A temporal property holds a value per frame;
// a distribution or a volume is aggregated over the frames, and what is published (md_script_eval_attributes) is
// the aggregate over every frame evaluated so far. Over other frames an aggregate is asked for through a SLOT (see
// SLOTS below), and the frames are evaluated once for all of them. Interrupted, nothing is lost: what remains to be
// done stays to be done.
//
// A distribution is also kept per frame while that fits in frame_budget: an aggregate over any frames is then summed
// from what is kept, and nothing is evaluated for it. The distributions are given room in script order; one that
// does not fit is aggregated as frames are evaluated, like a volume.
#define MD_SCRIPT_EVAL_FRAME_BUDGET_DEFAULT ((size_t)256 << 20)

typedef struct md_script_eval_desc_t {
    size_t num_frames;      // frames the evaluation can hold
    size_t frame_budget;    // bytes for distributions kept per frame. 0 for MD_SCRIPT_EVAL_FRAME_BUDGET_DEFAULT,
                            // and any amount less than one distribution of num_frames frames keeps none.
} md_script_eval_desc_t;

// Allocate and initialize the data for properties within the evaluation
// Should be performed as soon as the IR has changed.
md_script_eval_t* md_script_eval_create_desc(const md_script_ir_t* ir, const md_script_eval_desc_t* desc, struct md_allocator_i* alloc);

// With the default budget (md_script_eval_desc_t)
md_script_eval_t* md_script_eval_create(size_t num_frames, const md_script_ir_t* ir, struct md_allocator_i* alloc);

// As md_script_eval_create_desc, taking over from prev every property ir computes from the same thing (see
// md_script_ir_property_fingerprint), whatever it is called now: its values, and the frames they are of, which
// are not evaluated again. And its slots, with what they hold. A slot of a property which changed is kept as asked
// for, empty, as long as there is still an aggregate of its name. Nothing over another number of frames, which is
// another run. A distribution taken over is kept the way it was kept, per frame or not, whatever the budget. prev must
// be over the same system and run and must not be evaluating; it is not changed, and may be freed after. The versions of what is published, and of the slots, carry on from prev's, so one of the new
// evaluation never equals one of prev.
md_script_eval_t* md_script_eval_create_from(const md_script_eval_t* prev, const md_script_ir_t* ir, const md_script_eval_desc_t* desc, struct md_allocator_i* alloc);

// The frames a property holds: a temporal property the ones it has evaluated, a distribution or a volume the ones its
// published aggregate holds. Each property its own: one taken over by md_script_eval_create_from holds its frames
// while a new one holds none yet. A copy taken under the evaluation's lock, allocated with alloc, so it can be read
// while the evaluation runs.
md_bitset_t md_script_eval_property_frames(const md_script_eval_t* eval, str_t ident, struct md_allocator_i* alloc);

// The frames which remain to be evaluated, for the properties and for the slots, allocated with alloc
md_bitset_t md_script_eval_pending_frames(md_script_eval_t* eval, struct md_allocator_i* alloc);

// What the property was computed from, as it was compiled when the evaluation was made (see
// md_script_ir_property_fingerprint). 0 if there is no such property.
uint64_t md_script_eval_property_fingerprint(const md_script_eval_t* eval, str_t ident);

// ### WHAT AN EDIT CHANGES ###
// How the properties of a compilation stand against an evaluation made from an earlier one: what
// md_script_eval_create_from(eval, ir) would take over, and what would then remain to be evaluated. Nothing is
// evaluated or copied to find out, so it can be asked after every compilation, and while eval is being evaluated.
typedef enum md_script_eval_status_t {
    MD_SCRIPT_EVAL_STATUS_NEW = 0,      // nothing to take over: the evaluation has no property computed alike, nor of its name
    MD_SCRIPT_EVAL_STATUS_CHANGED,      // the evaluation has a property of its name computed from something else, which is let go
    MD_SCRIPT_EVAL_STATUS_HELD,         // taken over from a property computed alike, with the frames it holds
} md_script_eval_status_t;

typedef struct md_script_eval_plan_t {
    md_script_eval_status_t status;
    str_t    from;          // HELD: the property it is taken over from, which is named otherwise when it was renamed
    uint32_t num_frames;    // the frames it would hold
    uint32_t num_pending;   // the frames it would have evaluated, for itself and its slots
} md_script_eval_plan_t;

// One plan per property of ir, in its order, to out_plans, against an evaluation of num_frames frames made from eval.
// And the frames which would remain to be evaluated, for all of them and their slots, to out_pending when it is not
// NULL. Allocated with alloc. eval may be NULL: nothing is taken over. False, and nothing written, if ir is not valid
// or cap is less than its number of properties (md_script_ir_property_count).
bool md_script_eval_plan(md_script_eval_plan_t* out_plans, size_t cap, md_bitset_t* out_pending, const md_script_eval_t* eval,
                         const md_script_ir_t* ir, size_t num_frames, struct md_allocator_i* alloc);

// ### SLOTS ###
// A slot is a distribution or a volume over frames of the caller's choosing, asked for under a key of its own. Two
// callers asking under the same key for the same property share the slot.
//
// A slot moves toward the frames it is asked for one frame at a time, as frames are evaluated: each frame it wants
// and does not hold is added, each it holds and no longer wants is taken out again. Asked for other frames, it keeps
// what it has in common with them, and moves on from there; nothing restarts. At any moment it holds the aggregate
// over the frames it holds, exactly.
//
// A distribution kept per frame (see md_script_eval_desc_t) is summed over the frames when it is read, and nothing is
// evaluated for its slots: they hold the frames they want which are evaluated. A volume, or a distribution not kept
// per frame, is added up as frames are evaluated, also those the evaluation has evaluated before.
//
// Every function here may be called while frames are being evaluated, from another thread. A slot asked for new
// frames needs them evaluated: the frames which remain (md_script_eval_pending_frames) include them, and an
// evaluation running with an older list of frames does not.
#define MD_SCRIPT_EVAL_MAX_SLOTS 64

typedef struct md_script_eval_slot_info_t {
    uint64_t version;       // Changes whenever what the slot holds changes, and never repeats. Never 0.
    uint32_t num_frames;    // the frames it holds
    uint32_t num_pending;   // the frames it has yet to add or take out; 0 once it holds what it is asked for
    size_t   num_values;    // the values of the property: bins of every distribution in it, or voxels
} md_script_eval_slot_info_t;

// Asks for the property ident over frames (those below the evaluation's frame count) under key. Asked again for the
// same frames, nothing changes. False if ident is not a distribution or a volume, or there is no room for another slot.
bool md_script_eval_slot_request(md_script_eval_t* eval, uint64_t key, str_t ident, const md_bitset_t* frames);

// The slot is no longer needed: what it holds is let go
void md_script_eval_slot_release(md_script_eval_t* eval, uint64_t key, str_t ident);

// What the slot holds, without its values: cheap, to find out whether it holds anything new. False if there is no
// such slot.
bool md_script_eval_slot_info(md_script_eval_slot_info_t* out_info, const md_script_eval_t* eval, uint64_t key, str_t ident);

// Copies what the slot holds: the mean of the values over its frames, and for a distribution the mean of the weights,
// num_values each (see md_script_eval_slot_info_t), as the published aggregate holds them. out_weights may be NULL,
// and is not written for a volume. With no frames, values of 0 and weights of 1. A slot added up as frames are evaluated is copied under its lock, which
// holds back adding a frame to it until the copy is done, nothing else. False if there is no such slot, or
// num_values is too small, in which case only out_info is written.
bool md_script_eval_slot_extract(md_script_eval_slot_info_t* out_info, float* out_values, float* out_weights, size_t num_values,
                                 const md_script_eval_t* eval, uint64_t key, str_t ident);

uint64_t md_script_eval_ir_fingerprint(const md_script_eval_t* eval);

void md_script_eval_free(md_script_eval_t* eval);

// Clear before evaluating (computing data) frames
void md_script_eval_clear_data(md_script_eval_t* eval);

// Compute properties
// eval             : evaluation object to hold result
// ir               : holds an IR of the script to be evaluated
// sys              : the system, whose attribute table holds the run
// run              : "run/<name>", the run whose frames are evaluated (see RUNS in md_system.h)
// frames           : the frames to evaluate. Those with nothing left to do (md_script_eval_pending_frames) are skipped.
// Safe to call from several threads at once, also with overlapping frames (a frame is added to an aggregate once):
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
//                             a volume is {x, y, z}, an array of n volumes ('sdf(...) in
//                             residue(1:n)', '{sdf(a), sdf(b)}') {n, x, y, z}, the i-th volume a
//                             slice of the first axis, visualized with subidx i
//                             (md_script_vis_eval_ref). The unit is the unit of the values.
//   script/<ident>/range      rank 0, 2 components (min, max). The range the script declared for
//                             the property, or the observed one when it declared none. Temporal:
//                             the range of the values. Distribution: the range of the bin axis.
//   script/<ident>/mean       temporal with a population only: the per frame summary over it
//   script/<ident>/variance
//   script/<ident>/extent     2 components, (min, max)
//   script/<ident>/period     periodic temporal only (MD_SCRIPT_PROPERTY_FLAG_PERIODIC): rank 0, 2 components
//                             (beg, end). A value is the same as itself plus a whole number of end - beg. The values
//                             are given in [beg, end[ where the expression gives them there (a dihedral in [-PI, PI[),
//                             and are left as the expression gives them otherwise (dihedral(...) + 1).
//   script/<ident>/seam       periodic temporal only: rank 0. Where a plot of the values is best cut, in [beg, end[:
//                             where values of their kind are least likely to be, so that they cross it as seldom as
//                             possible. Shown in [seam, seam + end - beg[, a value crossing it jumps the whole period.
//                             The procedure states it (0 for a dihedral, 270 degrees for a pseudorotation phase), and
//                             arithmetic carries it along: the difference of two angles is cut half a period from 0.
//
// The summary of a periodic property is circular, around the circle its values lie on: the mean is the circular mean,
// in [beg, end[; the variance the square of the circular standard deviation sqrt(-2 ln R), at most (period / 2)^2,
// where mean +- its root covers the whole period; the extent the shortest arc which holds every value, (min, max) with
// min in [beg, end[ and max - min the length of the arc, so max may lie past end where the arc crosses it.
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

// The number of frames the evaluation holds
size_t md_script_eval_frame_count(const md_script_eval_t* eval);
// The frames every temporal property has evaluated (see md_script_eval_property_frames for one of them).
// A copy allocated with alloc, taken under the evaluation's lock.
md_bitset_t md_script_eval_frame_mask(const md_script_eval_t* eval, struct md_allocator_i* alloc);

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

// What the value of a visualized node is, in a few lines (a selection's atoms, a contact set's pairs, a plane's centroid
// and normal, the first values of a long array), at ctx->state for a value which depends on the frame. What hovering it
// shows next to its type. subidx picks one element of an array (-1 for all), as for md_script_vis_eval_ref.
// Writes a terminated string of at most cap - 1 characters and returns its length, 0 if there is nothing to say.
size_t md_script_vis_value_text(char* buf, size_t cap, md_script_vis_ref_t ref, int subidx, const md_script_vis_ctx_t* ctx);

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

// ### ATTRIBUTES ###

// The path that attr() reads an attribute of the system's table by (attributes, the table attr is in), in the form
// completion offers it: relative to its run ("edr/potential") where that names it alone, in full
// ("run/a/edr/potential") where it does not. A view into attr->path. Empty if attr() cannot read the attribute: it
// does not vary over the frames of a run, does not hold numbers, or has more dimensions per frame than a script value
// holds (attr() given its path says which).
str_t md_script_attr_path(const struct md_attributes_t* attributes, const struct md_attribute_t* attr);

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
    bool fits;                          // It has the type that is expected at the cursor: that of the parameter the
                                        // argument is for, or a selection at the top of a filter (md_filter_complete).
                                        // The values and the names of parameters always do. Where nothing in
                                        // particular is expected, or what an item gives is not known (the variables
                                        // of a script), it does not.
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
// ranked, which is left to the editor (with fits to go by). Nothing is offered in comments and numbers, and inside a
// string only when the argument it is given to takes known values. Values from the system are only offered if sys is
// not NULL. Everything returned (items and strings) is allocated from alloc, which is meant to be a temporary arena.
md_script_completions_t md_script_complete(str_t src, int cursor, const struct md_system_t* sys, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
