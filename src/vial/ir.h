#pragma once

// The array IR: what the checker lowers a script to, and what the interpreter runs.
//
// A value is a whole array, in SSA form: each is defined once, by an op or as a parameter of a region, and has the
// checker's type (dtype, shape with descriptions, unit, stage, pass). Ops are calls into runtime kernels: operators,
// indexing, built-ins. Control is structured: an if, a for, a map over dimensions and a fold own regions, and a region
// ends with a yield of what it gives back. A def is a function, one per specialization the checker made of it, until
// inlining removes them.
//
//     %3 : set[111]            = const residues(id=10..120, chain='A')      @ 5:10-5:38   static
//     %8 : float[111][3] Å     = com(%3)                                    @ 10:10-10:32 frame
//
// The text form is for golden tests and for looking at what the compiler did (vial_ir_dump).

#include "types.h"
#include "value.h"

typedef uint32_t vial_ssa_t;    // a value; 0 is none

typedef enum vial_ir_kind_t {
    VIAL_IR_NONE = 0,
    VIAL_IR_CONST,      // []; aux: the constant
    VIAL_IR_UNARY,      // [a]; sub: VIAL_OP_MINUS, VIAL_OP_PLUS or VIAL_OP_TILDE
    VIAL_IR_NOT,        // [a]: not, elementwise on bools
    VIAL_IR_BINARY,     // [a, b]; sub: the vial_op_t, NumPy's broadcasting
    VIAL_IR_FMA,        // [a, b, c]: a * b + c, rounded once
    VIAL_IR_COMPARE,    // [a, b]; sub: the vial_cmp_t
    VIAL_IR_LOGIC,      // [a, b]; sub: VIAL_KW_AND or VIAL_KW_OR, elementwise on bools
    VIAL_IR_SETOP,      // sub VIAL_KW_AND, VIAL_KW_OR: [a, b]; VIAL_KW_NOT: [scope, a]; VIAL_OP_MINUS: [a, b], a and not b
    VIAL_IR_CONVERT,    // [a]: to the result's dtype
    VIAL_IR_INDEX,      // [a, index values...]; aux: the first of its items in index_items, sub: how many
    VIAL_IR_LIST,       // [elements...]: stacked along a new first dimension
    VIAL_IR_TUPLE,      // [elements...]
    VIAL_IR_ELEM,       // [tuple]; aux: which element
    VIAL_IR_BUILTIN,    // [arguments by parameter, 0 when not given]; aux: the built-in; names: the parameters;
                        //   attr: in ints, each parameter's kind, when the check bound them
    VIAL_IR_CALL,       // [arguments]; aux: the function
    VIAL_IR_APPLY,      // [arguments]; region: the body of a nested def, which reads the names around it
    VIAL_IR_IF,         // [condition]; regions: then, else; results: what both yield
    VIAL_IR_FOR,        // [iterable, carried...]; region: body, its parameters the element and the carried values,
                        //   yielding the carried values; results: their last values
    VIAL_IR_MAP,        // [operands...]; region: one call, its parameters one element of each operand; attr: in ints,
                        //   how many leading dimensions of each operand are mapped (0: passed whole to every call);
                        //   results: the yields, stacked over the mapped dimensions broadcast
    VIAL_IR_KIND,       // [a]; sub: the kind (VX_KIND_*) a parameter takes: a set's positions, or a set of one
                        //   particle's position; a point as a group of one
    VIAL_IR_FOLD,       // [series, initial state]; region: the step, (state, value) -> state
    VIAL_IR_YIELD,      // [values...]: ends a region
    VIAL_IR_KIND_COUNT
} vial_ir_kind_t;

// VIAL_IR_KIND's sub beyond the parameter kinds: a set's particles, as an array of sets of one particle each, which a
// comprehension iterates over
#define VIAL_KIND_EACH 32

// VIAL_IR_MAP's sub: the region yields (value, keep) and the results hold only the kept values, as a comprehension
// with an if gives
#define VIAL_MAP_FILTERED 1

// One item of an index: a[i, 2:5, None, ..., mask]
typedef enum vial_index_kind_t {
    VIAL_INDEX_SCALAR = 1,  // an int: the dimension goes
    VIAL_INDEX_SLICE,       // start, stop, step: each a value or 0 when left out
    VIAL_INDEX_NEWAXIS,     // None
    VIAL_INDEX_ELLIPSIS,
    VIAL_INDEX_ARRAY,       // an int array or a bool mask; or on a set, positions of its particles
    VIAL_INDEX_SET,         // a set: a, the values of its particles along dimension c (an int) over the particles b
} vial_index_kind_t;

typedef struct vial_index_item_t {
    uint8_t     kind;       // vial_index_kind_t
    vial_ssa_t  a, b, c;    // SCALAR and ARRAY: a; SLICE: start, stop, step; SET: the set, the particles, the dimension
} vial_index_item_t;

typedef struct vial_ir_op_t {
    uint16_t        kind;       // vial_ir_kind_t
    uint16_t        sub;
    uint32_t        aux;
    uint32_t        arg_beg;    // into ir.args (and ir.names)
    uint32_t        num_args;
    vial_ssa_t      result;     // the first result; the others follow it
    uint32_t        num_results;
    uint32_t        region[2];  // owned regions, 0 when none
    uint32_t        parent;     // the region it is in
    uint32_t        attr;       // into ints, for the kinds that have attributes there
    md_vial_span_t  span;       // in the script: a def of the prelude has the span of the call that led to it
} vial_ir_op_t;

#define VIAL_IR_NO_OP UINT32_MAX

typedef struct vial_ir_value_t {
    vial_type_t     type;
    uint32_t        op;         // the op that defines it, or VIAL_IR_NO_OP for a region's parameter
    uint32_t        region;     // where it is defined
    uint32_t        index;      // which result of its op, or which parameter of its region
} vial_ir_value_t;

typedef struct vial_ir_region_t {
    uint32_t*       ops;        // md_array of op indices, in order
    vial_ssa_t      param_beg;
    uint32_t        num_params;
    uint32_t        owner;      // the op that owns it, or VIAL_IR_NO_OP for a function's body
    uint32_t        func;       // the function it is part of
} vial_ir_region_t;

typedef struct vial_ir_func_t {
    str_t           name;       // the def's name; empty for the script's top level
    uint32_t        body;       // its region, whose parameters are the def's
    md_vial_span_t  span;       // the call that made this specialization
} vial_ir_func_t;

// A published name: a top-level binding of the script
typedef struct vial_ir_output_t {
    str_t           name;
    vial_ssa_t      value;
    md_vial_span_t  span;
} vial_ir_output_t;

// What could not be lowered yet: the program checks, but cannot be evaluated
typedef struct vial_ir_failure_t {
    md_vial_span_t  span;
    char            what[96];
} vial_ir_failure_t;

typedef struct vial_ir_t {
    vial_ir_op_t*       ops;        // md_array
    vial_ir_value_t*    values;     // md_array by ssa; 0 is unused
    vial_ssa_t*         args;       // md_array
    str_t*              names;      // md_array beside args: a built-in's parameter names
    vial_ir_region_t*   regions;    // md_array; 0 is unused
    vial_ir_func_t*     funcs;      // md_array; 0 is the script's top level
    vial_val_t*         consts;     // md_array
    str_t*              const_text; // md_array beside consts: the source a constant was written as, or empty
    vial_index_item_t*  index_items;// md_array
    uint32_t*           ints;       // md_array: small integer attributes of ops
    vial_ir_output_t*   outputs;    // md_array
    vial_ir_failure_t*  failures;   // md_array
    uint32_t            insert;     // the region ops are emitted into
    const struct vial_desc_t* descs;// what dimensions run over, by the descs of the types' extents: the check's
    md_allocator_i*     alloc;
} vial_ir_t;

// An empty IR with function 0, the top level, ready to be emitted into
void vial_ir_init(vial_ir_t* ir, md_allocator_i* alloc);

// Makes op the owner of a region made before it, as its k-th
void       vial_ir_attach(vial_ir_t* ir, uint32_t op, uint32_t k, uint32_t region);

// A new region owned by op (VIAL_IR_NO_OP for a function body) in function func, with parameters of the given types
uint32_t   vial_ir_region(vial_ir_t* ir, uint32_t owner, uint32_t func, const vial_type_t* params, uint32_t num_params);

// A new function whose body has parameters of the given types; its body is not entered
uint32_t   vial_ir_func(vial_ir_t* ir, str_t name, md_vial_span_t span, const vial_type_t* params, uint32_t num_params);

// Appends an op to the insert region; its results have the given types and are numbered from the value returned
typedef struct vial_ir_emit_t {
    uint16_t            kind;
    uint16_t            sub;
    uint32_t            aux;
    const vial_ssa_t*   args;
    const str_t*        names;      // NULL, or one per argument
    uint32_t            num_args;
    const vial_type_t*  types;
    uint32_t            num_results;
    const uint32_t*     attr;       // copied into ints, or NULL
    uint32_t            num_attr;
    md_vial_span_t      span;
} vial_ir_emit_t;
uint32_t   vial_ir_emit(vial_ir_t* ir, const vial_ir_emit_t* e);

static inline vial_ssa_t vial_ir_result(const vial_ir_t* ir, uint32_t op) { return ir->ops[op].result; }

// Shorthands for one result
vial_ssa_t vial_ir_const(vial_ir_t* ir, const vial_val_t* v, const vial_type_t* t, str_t text, md_vial_span_t span);
vial_ssa_t vial_ir_op1(vial_ir_t* ir, uint16_t kind, uint16_t sub, uint32_t aux, const vial_ssa_t* args, uint32_t n,
                       const vial_type_t* t, md_vial_span_t span);
void       vial_ir_yield(vial_ir_t* ir, const vial_ssa_t* args, uint32_t n, md_vial_span_t span);

static inline const vial_type_t* vial_ir_type(const vial_ir_t* ir, vial_ssa_t v) { return &ir->values[v].type; }
static inline vial_ssa_t vial_ir_param(const vial_ir_t* ir, uint32_t region, uint32_t i) { return ir->regions[region].param_beg + i; }

const char* vial_ir_kind_str(vial_ir_kind_t kind);

// Checks that every value is defined before it is used, in its region or one enclosing it, that regions end in a
// yield of the right arity and that ops have the arguments their kind needs. Writes the first problem to err.
bool vial_ir_verify(const vial_ir_t* ir, char* err, size_t cap);

// The IR as text, functions first, one op per line with its span as line:col in src and its stage
str_t vial_ir_dump(const vial_ir_t* ir, str_t src, md_allocator_i* alloc);
