#pragma once

// Shared between the stages of the VIAL compiler. Not part of the public API (md_vial.h).

#include <md_vial.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>

#include <stdarg.h>

// ### DIAGNOSTICS ###
// Every stage reports through one list. A code's severity is its first letter; its message follows the wording of
// the reference's diagnostics table, filled in for the case at hand.

typedef struct vial_diags_t {
    md_vial_diagnostic_t* items;    // md_array
    str_t                 src;      // the source the spans point into
    uint32_t*             line_beg; // md_array: byte offset of each line's start, for line numbers in messages
    uint64_t*             at;       // a bit per byte of src: does a diagnostic start there
    md_allocator_i*       alloc;
} vial_diags_t;

void     vial_diags_init(vial_diags_t* d, str_t src, md_allocator_i* alloc);
void     vial_diag(vial_diags_t* d, const char code[6], md_vial_span_t span, const char* fmt, ...);
uint32_t vial_line_of(const vial_diags_t* d, uint32_t offset);   // 1-based
size_t   vial_diag_count(const vial_diags_t* d, const char code[6]);
void     vial_diags_sort(vial_diags_t* d);   // into source order, stable
bool     vial_diag_at(const vial_diags_t* d, uint32_t offset);   // does a diagnostic start at offset
void     vial_diag_set_fix(vial_diags_t* d, size_t index, str_t fix);   // copies fix

// ### TOKENS ###

typedef enum vial_kw_t {
    VIAL_KW_NONE = 0,
    // Python's keywords, all reserved so that no script changes meaning when a construct arrives
    VIAL_KW_FALSE, VIAL_KW_NONE_, VIAL_KW_TRUE, VIAL_KW_AND, VIAL_KW_AS, VIAL_KW_ASSERT, VIAL_KW_ASYNC, VIAL_KW_AWAIT,
    VIAL_KW_BREAK, VIAL_KW_CLASS, VIAL_KW_CONTINUE, VIAL_KW_DEF, VIAL_KW_DEL, VIAL_KW_ELIF, VIAL_KW_ELSE,
    VIAL_KW_EXCEPT, VIAL_KW_FINALLY, VIAL_KW_FOR, VIAL_KW_FROM, VIAL_KW_GLOBAL, VIAL_KW_IF, VIAL_KW_IMPORT,
    VIAL_KW_IN, VIAL_KW_IS, VIAL_KW_LAMBDA, VIAL_KW_NONLOCAL, VIAL_KW_NOT, VIAL_KW_OR, VIAL_KW_PASS, VIAL_KW_RAISE,
    VIAL_KW_RETURN, VIAL_KW_TRY, VIAL_KW_WHILE, VIAL_KW_WITH, VIAL_KW_YIELD,
    VIAL_KW_COUNT
} vial_kw_t;

typedef enum vial_op_t {
    VIAL_OP_NONE = 0,
    // Brackets and punctuation
    VIAL_OP_LPAREN, VIAL_OP_RPAREN, VIAL_OP_LBRACKET, VIAL_OP_RBRACKET, VIAL_OP_LBRACE, VIAL_OP_RBRACE,
    VIAL_OP_COMMA, VIAL_OP_COLON, VIAL_OP_SEMICOLON, VIAL_OP_DOT, VIAL_OP_RANGE /* .. */, VIAL_OP_ELLIPSIS,
    VIAL_OP_ARROW /* -> */, VIAL_OP_ASSIGN /* = */, VIAL_OP_WALRUS /* := */,
    // Arithmetic and bits
    VIAL_OP_PLUS, VIAL_OP_MINUS, VIAL_OP_STAR, VIAL_OP_SLASH, VIAL_OP_DSLASH, VIAL_OP_PERCENT, VIAL_OP_POW,
    VIAL_OP_AT, VIAL_OP_LSHIFT, VIAL_OP_RSHIFT, VIAL_OP_AMP, VIAL_OP_PIPE, VIAL_OP_CARET, VIAL_OP_TILDE,
    // Comparisons
    VIAL_OP_LT, VIAL_OP_GT, VIAL_OP_LE, VIAL_OP_GE, VIAL_OP_EQ, VIAL_OP_NE,
    // Augmented assignment
    VIAL_OP_PLUS_EQ, VIAL_OP_MINUS_EQ, VIAL_OP_STAR_EQ, VIAL_OP_SLASH_EQ, VIAL_OP_DSLASH_EQ, VIAL_OP_PERCENT_EQ,
    VIAL_OP_POW_EQ, VIAL_OP_AT_EQ, VIAL_OP_AMP_EQ, VIAL_OP_PIPE_EQ, VIAL_OP_CARET_EQ, VIAL_OP_LSHIFT_EQ,
    VIAL_OP_RSHIFT_EQ,
    VIAL_OP_COUNT
} vial_op_t;

enum {
    VIAL_STR_RAW    = 1 << 0,
    VIAL_STR_BYTES  = 1 << 1,
    VIAL_STR_FORMAT = 1 << 2,   // an f-string, which VIAL does not have
    VIAL_STR_TRIPLE = 1 << 3,
    VIAL_STR_UNCLOSED = 1 << 4,
    VIAL_TOK_REPORTED = 1 << 15,  // the lexer reported this token already (a closer that matched no opener)
};

typedef struct vial_token_t {
    md_vial_token_kind_t kind;
    md_vial_span_t       span;
    uint16_t             sub;       // vial_kw_t for keywords, vial_op_t for operators and punctuation
    uint16_t             flags;     // VIAL_STR_* for strings
} vial_token_t;

// Splits src into tokens, ending with MD_VIAL_TOKEN_EOF. Reports E0101, E0102 and E0103.
vial_token_t* vial_lex(str_t src, vial_diags_t* diags, md_allocator_i* alloc);

const char* vial_kw_str(vial_kw_t kw);

// ### UNITS ###
// The names that may follow a number, `10 ns`, and the strings, `1.5 'kJ/mol'`. Both fill unit when it is not NULL.

struct md_unit_t;
bool vial_unit_name(str_t name, struct md_unit_t* unit);
bool vial_unit_string(str_t text, struct md_unit_t* unit);
const char* vial_op_str(vial_op_t op);
const char* vial_token_kind_str(md_vial_token_kind_t kind);

// One line per token, "KIND beg..end text", for golden tests
str_t vial_tokens_dump(const vial_token_t* tokens, str_t src, md_allocator_i* alloc);

// The diagnostics of the prelude compiled with a program: none, unless the prelude has a mistake. For the tests.
size_t md_vial__prelude_diagnostics(const md_vial_program_t* prog, const md_vial_diagnostic_t** out);
