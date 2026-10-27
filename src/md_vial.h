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
} md_vial_diagnostic_t;

typedef enum md_vial_token_kind_t {
    MD_VIAL_TOKEN_NAME = 0,
    MD_VIAL_TOKEN_KEYWORD,
    MD_VIAL_TOKEN_INT,
    MD_VIAL_TOKEN_FLOAT,
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

#ifdef __cplusplus
}
#endif
