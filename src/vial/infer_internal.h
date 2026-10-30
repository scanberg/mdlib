#pragma once

// What the parts of the type checker share: infer.c (operators, indexing, statements) and infer_sets.c (selectors, the
// predefined sets, set algebra and the set functions). Later built-ins join the same way.

#include "infer.h"
#include "select.h"

typedef struct infer_t {
    const vial_tree_t*  t;
    str_t               src;
    vial_diags_t*       diags;
    const vial_names_t* n;
    vial_types_t*       out;
    vial_sys_t*         sys;        // NULL without a system
    vial_sets_t*        all;        // every particle, the scope of a script, once `not` needs it
    struct { str_t name; const vial_sets_t* sets; } predefined[8];     // evaluated once per compile
    md_allocator_i*     alloc;
} infer_t;

vial_type_t     vx_infer(infer_t* x, vial_node_id id);
void            vx_diag_fix(infer_t* x, const char code[6], vial_node_id at, const char* fix, const char* fmt, ...);
const char*     vx_tstr(const vial_type_t* t);
void            vx_set_value(infer_t* x, vial_node_id id, vial_value_t v);
vial_value_t    vx_value_of(const infer_t* x, vial_node_id id);
str_t           vx_text(const infer_t* x, vial_node_id id);

static inline vial_node_kind_t vx_kind(const infer_t* x, vial_node_id id) { return (vial_node_kind_t)x->t->nodes[id].kind; }
static inline vial_node_id vx_child(const infer_t* x, vial_node_id id, uint32_t i) { return vial_child(x->t, id, i); }
static inline uint32_t vx_num_children(const infer_t* x, vial_node_id id) { return x->t->nodes[id].child_count; }

// A description for the dimension of an array of sets; 0 when there is nothing to describe
uint32_t vx_desc(infer_t* x, uint8_t level, const vial_sets_t* sets);

// infer_sets.c

// A call to a built-in: true when it is one infer_sets.c types (selectors, len, count, labels), with its type in out
bool        vx_call_sets(infer_t* x, vial_node_id call, const vial_builtin_t* b, vial_type_t* out);

// A built-in written bare: a selector means all its units, a predefined set its particles. False for other names.
bool        vx_name_sets(infer_t* x, vial_node_id id, const vial_builtin_t* b, vial_type_t* out);

// The type and value of a set or an array of sets
vial_type_t vx_sets_type(infer_t* x, const vial_sets_t* s, vial_node_id at);

// Set algebra: the value of `a and b`, `a or b` and `not a` on known sets, or none
vial_value_t vx_sets_boolop(infer_t* x, bool is_and, const vial_value_t* a, const vial_value_t* b);
vial_value_t vx_sets_not(infer_t* x, const vial_value_t* a);

// The value of a subscript on a known set, an array of sets or an array of strings, with one index item; none when
// it is not known
vial_value_t vx_sets_subscript(infer_t* x, const vial_value_t* a, vial_node_id item);
