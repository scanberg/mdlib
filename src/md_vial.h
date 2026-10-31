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

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;

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

enum {
    MD_VIAL_COMPILE_NO_PRELUDE = 1 << 0,    // compile without the prelude: for the prelude itself, and tests
};

typedef struct md_vial_compile_desc_t {
    uint32_t flags;     // MD_VIAL_COMPILE_*

    // The system the script runs on. Selectors and the predefined sets are evaluated against it when compiling, so
    // their sizes, counts and labels are known and W0601, W0602 and E0207 are given. NULL compiles without one:
    // selectors are typed with extents not known. The checker reads it only during md_vial_compile, and a program
    // compiled against it holds its topology version. Stored selections and loaded runs arrive with later steps.
    const struct md_system_t* sys;

    // The number of frames of the run the script will be evaluated over, or 0 when not known. at(x, frame=k) with k
    // past the last frame is E0405, and a published per-frame value over 64 MB for the run is W0404. at() is read in a
    // step before evaluation reads the other frames, and remembered; compiling never waits for a frame.
    uint64_t num_frames;

    // The run the script is evaluated over: its prefix in sys->attributes, such as "run/md". attr(path) reads
    // <run>/<path>, with its unit and shape, and the run's frame axis, <run>/time, gives the number of frames when
    // num_frames is 0. Empty for none: attr() is then typed only when called, and its paths are not checked.
    str_t run;
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
typedef struct md_vial_output_t {
    str_t           name;
    md_vial_span_t  span;           // the name, where it is bound
    char            type[128];      // as the reference prints it: float[433] Å
    md_vial_stage_t stage;
    uint32_t        pass;           // the read over the trajectory in which it is complete; 0 when static
    bool            selection;
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

// ### EVALUATION ###
// A program evaluated over the frames of one run: its static part once, then sweeps over the frames in order, as many
// as its values need (see the reference's passes). The published values become attributes: "script/<name>", with the
// frame as the outer axis for a value per frame, beside "time", the evaluation's frame axis of ordinals. Selections
// are recipes and are not published as values. The interpreter is the reference of every other backend.
//
// Phase 3's API: one call evaluates every frame. Evaluating a range of frames, frames in parallel, and taking over
// what an earlier evaluation computed come with viamd's integration.

typedef struct md_vial_eval_t md_vial_eval_t;

typedef struct md_vial_eval_desc_t {
    const struct md_system_t* sys;  // the system the program was compiled against
    str_t       run;                // the run whose frames are evaluated: its prefix in sys->attributes, as "run/ala"
    size_t      max_elements;       // the most elements one value may have (R0105); 0 for 2^28
} md_vial_eval_desc_t;

// NULL when the program cannot be evaluated: it has errors, or something in it is not evaluated yet
md_vial_eval_t* md_vial_eval_create(const md_vial_program_t* prog, const md_vial_eval_desc_t* desc, struct md_allocator_i* alloc);
void            md_vial_eval_free(md_vial_eval_t* eval);

// Evaluates every frame. False when a value failed (its R diagnostic is in md_vial_eval_diagnostics; what does not
// depend on it is still computed) or when interrupted.
bool            md_vial_eval_run(md_vial_eval_t* eval);

// From any thread while md_vial_eval_run runs: how far it is, 0 to 1, and a request to stop between frames
double          md_vial_eval_progress(const md_vial_eval_t* eval);
void            md_vial_eval_interrupt(md_vial_eval_t* eval);

// What failed while evaluating: R diagnostics, each once, with the span of the expression
size_t          md_vial_eval_diagnostics(const md_vial_eval_t* eval, const md_vial_diagnostic_t** out);

// The published values, complete after md_vial_eval_run
const struct md_attributes_t* md_vial_eval_attributes(const md_vial_eval_t* eval);

#ifdef __cplusplus
}
#endif
