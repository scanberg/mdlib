#pragma once

// The VIAL syntax tree.
//
// The tree is a flat array of nodes, indexed by vial_node_id; node 0 is NONE and stands for an absent optional part
// (a slice bound left out, an if without else). A node's children are a range of ids in tree.children, in the order
// the layouts below give. Every node has the byte span of its text, so hover, diagnostics and the editor can point at
// any part of a script. Broken text still gives a tree: what failed to parse is an ERROR node, and a construct VIAL
// refuses is an UNSUPPORTED node around what it holds, so a refused class still has its body.

#include "vial_internal.h"

typedef uint32_t vial_node_id;

typedef enum vial_node_kind_t {
    VIAL_N_NONE = 0,

    // Statements                       children
    VIAL_N_MODULE,                      // statements
    VIAL_N_BLOCK,                       // statements
    VIAL_N_EXPR_STMT,                   // [expr]
    VIAL_N_ASSIGN,                      // [target, ..., value]: a = b = 1 has two targets
    VIAL_N_AUG_ASSIGN,                  // [target, value]; sub: the operator (VIAL_OP_PLUS_EQ, ...)
    VIAL_N_DEF,                         // [name, PARAMS, return annotation or NONE, BLOCK]
    VIAL_N_PARAMS,                      // [PARAM, ...]
    VIAL_N_PARAM,                       // [name, annotation or NONE, default or NONE]
    VIAL_N_RETURN,                      // [value or NONE]
    VIAL_N_IF,                          // [condition, BLOCK, else: BLOCK, IF (an elif) or NONE]
    VIAL_N_FOR,                         // [target, iterable, BLOCK]
    VIAL_N_WHILE,                       // [condition, BLOCK]
    VIAL_N_BREAK,
    VIAL_N_CONTINUE,
    VIAL_N_PASS,
    VIAL_N_ASSERT,                      // [condition, message or NONE]
    VIAL_N_UNSUPPORTED,                 // [what it holds, if anything]; sub: vial_unsupported_t
    VIAL_N_ERROR,                       // text that did not parse

    // Expressions
    VIAL_N_NAME,
    VIAL_N_INT,
    VIAL_N_FLOAT,
    VIAL_N_STRING,                      // adjacent string literals are one node
    VIAL_N_QUANTITY,                    // [number, UNIT]: 10 ns, 1.5 'kJ/mol'
    VIAL_N_UNIT,                        // a unit name, or sub 1: a unit string
    VIAL_N_CONST,                       // sub: VIAL_KW_TRUE, VIAL_KW_FALSE or VIAL_KW_NONE_
    VIAL_N_ELLIPSIS,
    VIAL_N_UNARY,                       // [operand]; sub: VIAL_OP_PLUS, VIAL_OP_MINUS or VIAL_OP_TILDE
    VIAL_N_NOT,                         // [operand]
    VIAL_N_BINARY,                      // [left, right]; sub: the operator
    VIAL_N_BOOLOP,                      // [operand, operand, ...]; sub: VIAL_KW_AND or VIAL_KW_OR
    VIAL_N_COMPARE,                     // [operand, CMPOP, operand, CMPOP, operand, ...]: a chain
    VIAL_N_CMPOP,                       // sub: vial_cmp_t
    VIAL_N_RANGE,                       // [low, high]: a..b
    VIAL_N_TERNARY,                     // [value, condition, otherwise]: value if condition else otherwise
    VIAL_N_CALL,                        // [callee, argument, ...]: an argument is an expression or a KWARG
    VIAL_N_KWARG,                       // [name, value]
    VIAL_N_ATTR,                        // [object, name]
    VIAL_N_SUBSCRIPT,                   // [object, index]: the index is an expression, a SLICE, or a TUPLE of them
    VIAL_N_SLICE,                       // [start or NONE, stop or NONE, step or NONE]
    VIAL_N_TUPLE,                       // [element, ...]
    VIAL_N_LIST,                        // [element, ...]
    VIAL_N_LISTCOMP,                    // [element, COMP_FOR, ...]
    VIAL_N_COMP_FOR,                    // [target, iterable, condition, ...]

    VIAL_N_COUNT
} vial_node_kind_t;

typedef enum vial_cmp_t {
    VIAL_CMP_LT = 1, VIAL_CMP_GT, VIAL_CMP_LE, VIAL_CMP_GE, VIAL_CMP_EQ, VIAL_CMP_NE,
    VIAL_CMP_IN, VIAL_CMP_NOT_IN, VIAL_CMP_IS, VIAL_CMP_IS_NOT,
} vial_cmp_t;

// What an UNSUPPORTED node is. Each is reported once, where it is written: E0803 for what VIAL never has, E0801 for
// what it will have, E0701 for a v1 spelling.
typedef enum vial_unsupported_t {
    VIAL_UNSUP_NONE = 0,
    // E0803
    VIAL_UNSUP_CLASS, VIAL_UNSUP_IMPORT, VIAL_UNSUP_LAMBDA, VIAL_UNSUP_EXCEPTION, VIAL_UNSUP_GENERATOR,
    VIAL_UNSUP_GLOBAL, VIAL_UNSUP_DECORATOR, VIAL_UNSUP_WITH, VIAL_UNSUP_ASYNC, VIAL_UNSUP_DICT, VIAL_UNSUP_SET,
    VIAL_UNSUP_STARRED, VIAL_UNSUP_FSTRING, VIAL_UNSUP_DEL, VIAL_UNSUP_WALRUS, VIAL_UNSUP_LOOP_ELSE,
    VIAL_UNSUP_ANNOTATION, VIAL_UNSUP_BYTES,
    // E0801
    VIAL_UNSUP_ITEM_ASSIGN,
    // E0701
    VIAL_UNSUP_V1_RANGE,
} vial_unsupported_t;

typedef struct vial_node_t {
    uint8_t         kind;           // vial_node_kind_t
    uint8_t         flags;
    uint16_t        sub;
    md_vial_span_t  span;
    uint32_t        child_beg;      // into tree.children
    uint32_t        child_count;
} vial_node_t;

typedef struct vial_tree_t {
    vial_node_t*   nodes;           // md_array; nodes[0] is NONE
    vial_node_id*  children;        // md_array
    vial_node_id   root;            // a MODULE
} vial_tree_t;

static inline const vial_node_t* vial_node(const vial_tree_t* t, vial_node_id id) { return &t->nodes[id]; }

static inline vial_node_id vial_child(const vial_tree_t* t, vial_node_id id, uint32_t i) {
    const vial_node_t* n = &t->nodes[id];
    return i < n->child_count ? t->children[n->child_beg + i] : 0;
}

// Parses tokens (from vial_lex, comments included) into a tree. Reports E0101 and E0103 for text that does not parse,
// E0701 for v1 spellings, E0801 for what is not supported yet and E0803 for what VIAL never has.
vial_tree_t vial_parse(const vial_token_t* tokens, str_t src, vial_diags_t* diags, md_allocator_i* alloc);

const char* vial_node_kind_str(vial_node_kind_t kind);
const char* vial_cmp_str(vial_cmp_t cmp);

// The tree as indented text, one node per line: "kind beg..end [sub] [text]". For golden tests.
str_t vial_tree_dump(const vial_tree_t* tree, str_t src, md_allocator_i* alloc);
