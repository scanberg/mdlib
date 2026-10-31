#pragma once

// VIAL, the VIAMD Analysis Language: version 2 of the script.
//
// The language is specified in docs/vial_reference.md, and the design behind it in the VIAMD Script v2 plan. Its
// examples are the test suite (unittest/test_vial_spec.c). VIAL is built in phases; this header grows with them and
// declares only what is implemented:
//
//   phase 1   the front end: tokens with byte spans, the syntax tree with error recovery, syntax diagnostics
//   phase 2   types and checking: dtypes, shapes, units, stages, passes, the built-in table, all compile-time codes
//   phase 3   the array IR and the interpreter
//
// v1 (md_script.h) is untouched; a workspace uses one language or the other.

#include <core/md_str.h>
#include <core/md_unit.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;
struct md_system_state_t;

// A range of bytes in the source, [beg, end)
typedef struct md_vial_span_t {
    uint32_t beg;
    uint32_t end;
} md_vial_span_t;

typedef enum md_vial_severity_t {
    MD_VIAL_SEVERITY_ERROR = 0,     // E....: the script does not compile
    MD_VIAL_SEVERITY_WARNING,       // W....: it compiles
    MD_VIAL_SEVERITY_RUNTIME,       // R....: found while evaluating a frame
} md_vial_severity_t;

// One diagnostic, with the code and message of the reference's diagnostics table
typedef struct md_vial_diagnostic_t {
    char                code[6];    // "E0301"
    md_vial_severity_t  severity;
    md_vial_span_t      span;
    str_t               message;
    str_t               fix;        // text to put in place of span, the table's "Fix offered"; empty when there is none
} md_vial_diagnostic_t;

typedef enum md_vial_token_kind_t {
    MD_VIAL_TOKEN_NAME = 0,
    MD_VIAL_TOKEN_KEYWORD,
    MD_VIAL_TOKEN_INT,
    MD_VIAL_TOKEN_FLOAT,
    MD_VIAL_TOKEN_UNIT,         // the unit of a quantity literal: ns in `10 ns`
    MD_VIAL_TOKEN_STRING,
    MD_VIAL_TOKEN_OPERATOR,
    MD_VIAL_TOKEN_PUNCT,
    MD_VIAL_TOKEN_COMMENT,
    MD_VIAL_TOKEN_NEWLINE,
    MD_VIAL_TOKEN_INDENT,
    MD_VIAL_TOKEN_DEDENT,
    MD_VIAL_TOKEN_ERROR,
    MD_VIAL_TOKEN_EOF,
} md_vial_token_kind_t;

typedef struct md_vial_token_t {
    md_vial_token_kind_t kind;
    md_vial_span_t       span;
} md_vial_token_t;

typedef struct md_vial_program_t md_vial_program_t;

// A stored selection a script reads by name, as viamd keeps them: its particles, indices of the system's atoms
typedef struct md_vial_named_set_t {
    str_t           name;
    const uint32_t* atoms;
    size_t          num_atoms;
} md_vial_named_set_t;

enum {
    MD_VIAL_COMPILE_NO_PRELUDE = 1 << 0,    // compile without the prelude: for the prelude itself, and tests
};

typedef struct md_vial_compile_desc_t {
    uint32_t flags;     // MD_VIAL_COMPILE_*

    // The system the script runs on. Selectors and the predefined sets are evaluated against it when compiling, so
    // their sizes, counts and labels are known and W0601, W0602 and E0207 are given. NULL compiles without one:
    // selectors are typed with extents not known. The checker reads it only during md_vial_compile, and a program
    // compiled against it holds its topology version. Stored selections are given below.
    const struct md_system_t* sys;

    // The number of frames of the run the script will be evaluated over, or 0 when not known. at(x, frame=k) with k
    // past the last frame is E0405, and a published per-frame value over 64 MB for the run is W0404. at() is read in a
    // step before evaluation reads the other frames, and remembered; compiling never waits for a frame.
    uint64_t num_frames;

    // The run the script is evaluated over: its prefix in sys->attributes, such as "run/md". attr(path) reads
    // <run>/<path>, with its unit and shape, and the run's frame axis, <run>/time, gives the number of frames when
    // num_frames is 0. Empty for none: attr() is then typed only when called, and its paths are not checked.
    str_t run;

    // Stored selections the script may read by name: a static set each, as a predefined set is. A name the script
    // binds hides one; a name that is not valid (md_vial_name_valid) is not read. Copied when compiling.
    const md_vial_named_set_t* sets;
    size_t                     num_sets;
} md_vial_compile_desc_t;

// Compiles a script: the front end, then the checker, which types every expression, evaluates the static ones it can
// (selectors, counts, literals) against desc->sys, and gives every compile-time diagnostic of the reference. Never
// fails on broken text, which gives diagnostics; NULL only when out of memory.
// desc may be NULL. The program keeps its own copy of src. Any text is safe to compile: deep nesting is reported
// (E0104) rather than followed, so compiling needs less than 256 KB of stack, and time and memory grow linearly with
// the length of src.
md_vial_program_t* md_vial_compile(str_t src, const md_vial_compile_desc_t* desc, struct md_allocator_i* alloc);
void               md_vial_program_free(md_vial_program_t* prog);

// Every diagnostic, in source order
size_t md_vial_diagnostics(const md_vial_program_t* prog, const md_vial_diagnostic_t** out);

// Every token, comments included, in source order, for highlighting
size_t md_vial_tokens(const md_vial_program_t* prog, const md_vial_token_t** out);

typedef enum md_vial_stage_t {
    MD_VIAL_STAGE_STATIC = 0,           // once, when compiling
    MD_VIAL_STAGE_FRAME,                // once per frame
    MD_VIAL_STAGE_TRAJECTORY,           // over the frames it needs
    MD_VIAL_STAGE_TRAJECTORY_PER_FRAME, // over the frames, keeping the frame, as smooth() gives
} md_vial_stage_t;

// One published value: a top-level name that does not start with _, bound to a number, a bool or a set, of at most
// two dimensions. A set or an array of sets is a selection; a frame-stage one is a recipe, evaluated for the frame
// shown. In the order the script binds them.
enum {
    MD_VIAL_OUTPUT_TEMPORAL     = 1 << 0,   // a value per frame: frame stage, or over the trajectory keeping the frame
    MD_VIAL_OUTPUT_DISTRIBUTION = 1 << 1,   // a histogram over the frames
    MD_VIAL_OUTPUT_PERIODIC     = 1 << 2,   // values that wrap, such as angles about an axis
};

typedef struct md_vial_output_t {
    str_t           name;
    md_vial_span_t  span;           // the name, where it is bound
    char            type[128];      // as the reference prints it: float[433] Å
    md_vial_stage_t stage;
    uint32_t        pass;           // the read over the trajectory in which it is complete; 0 when static
    bool            selection;
    uint32_t        flags;          // MD_VIAL_OUTPUT_*
    // What it is computed from: the ops that compute it, their constants and types, the system's topology and the
    // attributes it reads, not its name or where it is written. An evaluation of another compilation holds what this
    // one would compute for it when the two agree (md_vial_eval_create_from). 0 when the program does not evaluate.
    uint64_t        fingerprint;
} md_vial_output_t;

size_t md_vial_outputs(const md_vial_program_t* prog, const md_vial_output_t** out);

// What hover shows of the innermost expression at a byte offset of the source, as the checker knows it: its span,
// type, stage, pass and, where known when compiling, its value. False where no expression with a type is.
typedef struct md_vial_hover_t {
    md_vial_span_t  span;
    char            type[128];
    md_vial_stage_t stage;
    uint32_t        pass;
    char            value[96];      // empty when not known when compiling
} md_vial_hover_t;

bool md_vial_hover(const md_vial_program_t* prog, uint32_t offset, md_vial_hover_t* out);

// Expressions by id: what the editor marks, hovers and draws. An id names one expression of one program, 0 none.
// The innermost expression with a type at a byte offset, as md_vial_hover finds it; and what hover shows of one.
uint32_t md_vial_expr_at(const md_vial_program_t* prog, uint32_t offset);
bool     md_vial_expr_info(const md_vial_program_t* prog, uint32_t expr, md_vial_hover_t* out);

// Every expression with a type, in source order, for the editor to mark: its span and how many marked expressions
// hold it
typedef struct md_vial_mark_t {
    uint32_t        expr;
    md_vial_span_t  span;
    uint32_t        depth;
} md_vial_mark_t;

size_t md_vial_marks(const md_vial_program_t* prog, const md_vial_mark_t** out);

// ### EDITOR ###

// The words the language reserves, every keyword of Python, for highlighting: how many there are, and up to cap of
// them into out, which may be NULL
size_t md_vial_keywords(str_t* out, size_t cap);

// Every name a script reads without binding it, for highlighting: the built-ins, the predefined sets and the prelude's
// defs, not those starting with _. How many there are, and up to cap of them into out, which may be NULL.
size_t md_vial_builtin_names(str_t* out, size_t cap);

// Whether a script, or a stored selection it reads, may bind name: an identifier that is no keyword and no built-in name
bool md_vial_name_valid(str_t name);

// The path attr() reads an attribute by, in a script evaluated over run: its name past "<run>/", empty when the
// attribute is not one of the run's
str_t md_vial_attr_path(str_t run, str_t attribute_name);

typedef enum md_vial_completion_kind_t {
    MD_VIAL_COMPLETION_KEYWORD = 0,
    MD_VIAL_COMPLETION_FUNCTION,    // a built-in, a selector or a def, to call
    MD_VIAL_COMPLETION_CONSTANT,    // a predefined set, such as protein
    MD_VIAL_COMPLETION_VARIABLE,    // a name the script binds where the cursor is
    MD_VIAL_COMPLETION_PARAMETER,   // a parameter of the call at the cursor, as name=
    MD_VIAL_COMPLETION_PROPERTY,    // after a dot: .xyz, .com, and methods
    MD_VIAL_COMPLETION_VALUE,       // inside a string: a value from the system, such as a residue name for resname=
} md_vial_completion_kind_t;

typedef struct md_vial_completion_t {
    str_t                       label;  // what to show, and match what has been typed against
    str_t                       text;   // what replaces the span
    md_vial_completion_kind_t   kind;
} md_vial_completion_t;

typedef struct md_vial_completions_t {
    md_vial_span_t          span;       // the bytes a completion replaces
    str_t                   prefix;     // what has been typed of it
    md_vial_completion_t*   items;      // those that start with prefix: parameters, names, then keywords; each once
    size_t                  count;
} md_vial_completions_t;

// What may be written at a byte offset of prog's source, from what prog knows there: the names in scope, the
// parameters of the call around the cursor, and inside a string a selector keyword's values in the system prog was
// compiled against, or attr()'s paths in its run. Everything from alloc.
bool md_vial_complete(md_vial_completions_t* out, const md_vial_program_t* prog, uint32_t cursor, struct md_allocator_i* alloc);

// ### EVALUATION ###
// A program evaluated over the frames of one run: its static part once, then sweeps over the frames in order, as many
// as its values need (see the reference's passes). The published values become attributes: "script/<name>", with the
// frame as the outer axis for a value per frame, beside "time", the evaluation's frame axis of ordinals. Selections
// are recipes and are not published as values. The interpreter is the reference of every other backend.
//
// Frames are evaluated on several threads when asked, with the same bits as on one: each frame is computed alone, and
// the trajectory operations after a sweep run on one thread in a fixed order. An evaluation made from an earlier one
// takes over what that one computed alike (md_vial_eval_create_from).

typedef struct md_vial_eval_t md_vial_eval_t;

// The caller's thread pool: calls fn(arg, i) for every i in [0, count), on up to num_threads threads at once, and
// returns when every call has returned. viamd passes its task system.
typedef void (*md_vial_parallel_for_t)(void* user, size_t count, void (*fn)(void* arg, size_t i), void* arg);

typedef struct md_vial_eval_desc_t {
    const struct md_system_t* sys;  // the system the program was compiled against
    str_t       run;                // the run whose frames are evaluated: its prefix in sys->attributes, as "run/ala"
    size_t      max_elements;       // the most elements one value may have (R0105); 0 for 2^28
    // Threads for the frames: 0 or 1 evaluates them on the calling thread. With parallel_for NULL, md_vial_eval_run
    // starts num_threads - 1 threads of its own for each sweep over the frames and works on the calling one.
    size_t                  num_threads;
    md_vial_parallel_for_t  parallel_for;
    void*                   parallel_user;
} md_vial_eval_desc_t;

// NULL when the program cannot be evaluated: it has errors, or something in it is not evaluated yet
md_vial_eval_t* md_vial_eval_create(const md_vial_program_t* prog, const md_vial_eval_desc_t* desc, struct md_allocator_i* alloc);

// The same, taking over from prev every output computed alike (md_vial_output_t.fingerprint), whatever it is called
// now: its values are published as prev published them, and nothing is evaluated for it. Only from an evaluation of
// the same run that has finished every frame; prev is not changed, may still be evaluating something else's, and may
// be freed after. Its program must still exist while this is called.
md_vial_eval_t* md_vial_eval_create_from(const md_vial_eval_t* prev, const md_vial_program_t* prog, const md_vial_eval_desc_t* desc,
                                         struct md_allocator_i* alloc);

// How each output of prog stands against prev, without evaluating anything: what md_vial_eval_create_from(prev, prog)
// would take over, and how many frames would be evaluated for the rest
typedef enum md_vial_plan_status_t {
    MD_VIAL_PLAN_NEW = 0,       // nothing to take over
    MD_VIAL_PLAN_CHANGED,       // prev has an output of its name computed otherwise
    MD_VIAL_PLAN_HELD,          // taken over from an output computed alike
} md_vial_plan_status_t;

typedef struct md_vial_plan_t {
    md_vial_plan_status_t status;
    str_t    from;              // HELD: the output it is taken over from, which is named otherwise when it was renamed
    uint32_t num_frames;        // the frames it would hold
    uint32_t num_pending;       // the frames that would be evaluated for it
} md_vial_plan_t;

// One plan per output of prog (md_vial_outputs), in its order, against an evaluation of num_frames frames. prev may be
// NULL. False, and nothing written, when cap is less than the number of outputs.
bool md_vial_eval_plan(md_vial_plan_t* out, size_t cap, const md_vial_eval_t* prev, const md_vial_program_t* prog, size_t num_frames);
void            md_vial_eval_free(md_vial_eval_t* eval);

// Evaluates every frame. False when a value failed (its R diagnostic is in md_vial_eval_diagnostics; what does not
// depend on it is still computed) or when interrupted.
bool            md_vial_eval_run(md_vial_eval_t* eval);

// From any thread while md_vial_eval_run runs: how far it is, 0 to 1, and a request to stop between frames
double          md_vial_eval_progress(const md_vial_eval_t* eval);
void            md_vial_eval_interrupt(md_vial_eval_t* eval);

// What failed while evaluating: R diagnostics, each once, with the span of the expression
size_t          md_vial_eval_diagnostics(const md_vial_eval_t* eval, const md_vial_diagnostic_t** out);

// The published values. A value per frame whose shape the type knows is there from the start of md_vial_eval_run, NaN
// (0 for ints) at the frames not computed yet, and written as frames finish, from the evaluation's threads; its
// version changes when the run ends. Everything else is there once md_vial_eval_run has finished.
const struct md_attributes_t* md_vial_eval_attributes(const md_vial_eval_t* eval);

// Which frames an output holds, from any thread while md_vial_eval_run runs: out[f] is 1 where it is computed. Every
// frame once the run has finished, and for an output taken over. Returns the number of frames; writes nothing when cap
// is less.
size_t md_vial_eval_output_frames(const md_vial_eval_t* eval, str_t name, uint8_t* out, size_t cap);

// An output over the trajectory, such as a histogram over the frames, computed again over chosen frames: those where
// frames[f] is nonzero, from the per-frame values the evaluation kept, so nothing is evaluated again. What a plot
// period shows (decided 31 October). Into out as floats, NaN where missing; returns how many values it has, and writes
// them only when cap holds them all. 0 when the run has not finished, or the output is not computed over the frames
// from per-frame values. Safe from any thread once md_vial_eval_run has returned.
size_t md_vial_eval_aggregate(const md_vial_eval_t* eval, str_t name, const uint8_t* frames, float* out, size_t cap);

// The value of the expression md_vial_hover shows at offset, at a frame of the run, as the reference prints values,
// cut at cap. It is computed for that value alone, at most 16 million entries, with a reader and an interpreter of its
// own, so it may run while md_vial_eval_run does. A value over the trajectory needs every frame: it is the
// evaluation's, once md_vial_eval_run has finished, where the evaluation computed it. False when there is none: no
// expression there, one inside a def's body, a value over the trajectory not (yet) computed, or one that fails or is
// past the budget.
bool md_vial_eval_hover(const md_vial_eval_t* eval, uint32_t offset, int64_t frame, char* buf, size_t cap);

// ### DRAWING ###
// What hover draws of an expression, as the reference's Drawing section has it: the drawing a built-in or a def
// attached with shown(), else the value by its type, a set or an array of sets as its particles, points as markers.
// In mdlib's own layout, which viamd converts to its vertices; no colours or styles.

typedef enum md_vial_prim_kind_t {
    MD_VIAL_PRIM_MARKER = 0,    // p[0]
    MD_VIAL_PRIM_LINE,          // p[0] to p[1], which is the image of the end nearest p[0], so no segment spans the box
    MD_VIAL_PRIM_ARROW,         // p[0] to p[1], the head at p[1]
    MD_VIAL_PRIM_SPHERE,        // centre p[0], radius r
    MD_VIAL_PRIM_ARC,           // the angle at p[0] from p[1] to p[2], both nearest images
    MD_VIAL_PRIM_DISK,          // centre p[0], unit normal p[1], radius r: a plane
} md_vial_prim_kind_t;

typedef struct md_vial_prim_t {
    md_vial_prim_kind_t kind;
    uint32_t            element;    // the element of the value it draws, in row-major order; 0 for a scalar
    float               p[3][3];    // Å
    float               r;
} md_vial_prim_t;

// Text at a point: a number with its unit, which viamd may format in units of its own, or a string (value NaN)
typedef struct md_vial_label_t {
    float       pos[3];
    uint32_t    element;
    str_t       text;       // always set: the number with two decimals and its unit, or the string
    double      value;
    md_unit_t   unit;
} md_vial_label_t;

// Particles highlighted: atoms[beg, end) of the drawing, ascending
typedef struct md_vial_group_t {
    uint32_t    element;
    uint32_t    beg;
    uint32_t    end;
} md_vial_group_t;

typedef struct md_vial_drawing_t {
    md_vial_prim_t*     prims;
    size_t              num_prims;
    md_vial_label_t*    labels;
    size_t              num_labels;
    md_vial_group_t*    groups;
    size_t              num_groups;
    uint32_t*           atoms;          // the groups' particles: indices of the system's atoms
    size_t              num_elements;   // of the value drawn
    size_t              num_omitted;    // primitives and labels left out past 10,000
} md_vial_drawing_t;

typedef struct md_vial_draw_desc_t {
    const md_vial_program_t*        prog;
    const struct md_system_t*       sys;    // the system prog was compiled against
    str_t                           run;    // the run: at() reads its frames, and the frame drawn when state is NULL
    int64_t                         frame;  // the frame shown: frame(), time(), and a value over the trajectory
    // The coordinates shown, which may be interpolated or recentred: values per frame are computed from them. NULL
    // reads frame from the run.
    const struct md_system_state_t* state;
    // An evaluation of prog over the run, for values over the trajectory once md_vial_eval_run has finished; may be
    // NULL, and those are then not drawn
    const md_vial_eval_t*           eval;
    int64_t                         element;    // the element of an array drawn alone; -1 for every element
} md_vial_draw_desc_t;

// What an expression (md_vial_expr_at) draws, its arrays from alloc. The expression is computed alone,
// as hover computes it, with its own reader and interpreter, so it may run while md_vial_eval_run does. False when
// there is none there, it fails, or it is over the trajectory and not computed; true with nothing in it when the
// value draws nothing, such as a number, which hover shows as text.
bool md_vial_draw(md_vial_drawing_t* out, const md_vial_draw_desc_t* desc, uint32_t expr, struct md_allocator_i* alloc);

// The value of an expression where desc draws, as the reference prints values, cut at cap: what hover shows beside the
// drawing. Returns its length, 0 when it has none there.
size_t md_vial_value_text(char* buf, size_t cap, const md_vial_draw_desc_t* desc, uint32_t expr);

// ### FILTERS ###
// An expression that picks particles, as viamd's representations, groups and queries write them: compiled alone as one
// statement, with the stored selections and the static selections of a program as names it may read.

typedef struct md_vial_filter_desc_t {
    const struct md_system_t*       sys;        // required
    // The coordinates, for an expression that depends on them, such as within(): it is then dynamic. NULL fails one.
    const struct md_system_state_t* state;
    str_t                           run;        // for attr(); may be empty
    const md_vial_named_set_t*      sets;       // stored selections it may read by name
    size_t                          num_sets;
    // A program whose published selections it may read by name, below the stored ones: a static one as it is, and a
    // recipe evaluated at state, which makes the filter dynamic. One set each; NULL for none.
    const md_vial_program_t*        prog;
    // The particles it picks within: what it picks is cut to these, and sets left empty dropped. NULL for all.
    const uint32_t*                 scope;
    size_t                          num_scope;
} md_vial_filter_desc_t;

typedef struct md_vial_filter_t {
    bool                    ok;
    bool                    dynamic;    // what it picks depends on the coordinates: filter again when they change
    // What it picks: one set, or an array of them, row-major, with its shape; bonds(...) is [n][2]
    uint8_t                 ndim;       // 0 for one set
    uint32_t                dims[2];
    size_t                  num_sets;
    uint32_t*               offsets;    // [num_sets + 1] into atoms
    uint32_t*               atoms;      // each set's particles, ascending
    // When not ok: the first error, its span in expr, empty when it is about the whole
    md_vial_diagnostic_t    error;
} md_vial_filter_t;

bool md_vial_filter(md_vial_filter_t* out, str_t expr, const md_vial_filter_desc_t* desc, struct md_allocator_i* alloc);

// What may be written at a byte offset of expr, as md_vial_complete; the span is in expr
bool md_vial_filter_complete(md_vial_completions_t* out, str_t expr, uint32_t cursor, const md_vial_filter_desc_t* desc,
                             struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
