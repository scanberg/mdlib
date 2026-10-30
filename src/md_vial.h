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
    // past the last frame is E0405, and a published per-frame value over 64 MB for the run is W0404. at() reads its
    // frame when compiling once VIAL evaluates (phase 3); until then its value is typed but not known.
    uint64_t num_frames;

    // The run the script is evaluated over: its prefix in sys->attributes, such as "run/md". attr(path) reads
    // <run>/<path>, with its unit and shape, and the run's frame axis, <run>/time, gives the number of frames when
    // num_frames is 0. Empty for none: attr() is then typed only when called, and its paths are not checked.
    str_t run;
} md_vial_compile_desc_t;

// Compiles a script. Phase 1 stops after the parser: the program holds the tokens, the syntax tree and the
// diagnostics of the front end. Never fails on broken text, which gives diagnostics; NULL only when out of memory.
// desc may be NULL. The program keeps its own copy of src. Any text is safe to compile: deep nesting is reported
// (E0104) rather than followed, so compiling needs less than 256 KB of stack, and time and memory grow linearly with
// the length of src.
md_vial_program_t* md_vial_compile(str_t src, const md_vial_compile_desc_t* desc, struct md_allocator_i* alloc);
void               md_vial_program_free(md_vial_program_t* prog);

// Every diagnostic, in source order
size_t md_vial_diagnostics(const md_vial_program_t* prog, const md_vial_diagnostic_t** out);

// Every token, comments included, in source order, for highlighting
size_t md_vial_tokens(const md_vial_program_t* prog, const md_vial_token_t** out);

#ifdef __cplusplus
}
#endif
