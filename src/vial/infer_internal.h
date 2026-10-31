#pragma once

// What the parts of the type checker share: infer.c (operators, indexing, statements) and infer_sets.c (selectors, the
// predefined sets, set algebra and the set functions). Later built-ins join the same way.

#include "infer.h"
#include "select.h"

// The results of defs already checked, by the def and its arguments' types
typedef struct vx_memo_t {
    md_hashmap32_t      map;        // key -> index + 1
    vial_type_t*        results;    // md_array
    vial_value_t*       values;
    uint32_t*           funcs;      // lowering: the IR function of each, 0 when none was emitted
} vx_memo_t;

// A call to a def being checked: its own types for the def's nodes and symbols, in the def's tree
typedef struct vx_frame_t {
    struct vx_frame_t*  parent;
    const vial_tree_t*  t;          // the def's tree: the script's or the prelude's
    vial_node_id        def;
    bool                prelude;    // the def is the prelude's: its diagnostics are reported in the script
    vial_node_id        site;       // prelude: the call in the script that led here, where its diagnostics go
    md_hashmap32_t      syms;       // symbol -> slot + 1
    vial_type_t*        sym_types;  // md_array by slot
    vial_value_t*       sym_values;
    vial_node_id*       sym_prov;   // prelude: the script node a parameter's argument was written as, or 0
    uint32_t*           sym_ssa;    // lowering: the symbol's value in the IR
    vial_type_t         result;     // the returned values joined
    vial_value_t        value;      // the value returned, when there is one return and it is known
    bool                has_result;
    uint32_t            depth;
    bool                lower_failed;   // lowering: something in the body was not lowered
} vx_frame_t;

typedef struct infer_t {
    const vial_tree_t*  t;
    str_t               src;
    vial_diags_t*       diags;
    const vial_names_t* n;
    vial_types_t*       out;
    vial_types_t*       top;        // the script's top level: its names' types
    vial_desc_t**       descs;      // what dimensions run over, one array for every frame: the top level's
    vial_sys_t*         sys;        // NULL without a system
    vial_sets_t*        all;        // every particle, the scope of a script, once `not` needs it
    struct { str_t name; const vial_sets_t* sets; } predefined[16];    // predefined sets and bare selectors, once per compile
    const vial_check_t* c;          // the script's and the prelude's trees and names
    vx_frame_t*         frame;      // the def being checked, or NULL at top level
    vx_memo_t*          memo;       // a def's result by its arguments' types, shared by every frame
    struct vial_ir_t*   ir;         // lowering: the IR emitted into, or NULL when nothing is emitted
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

bool vx_broadcast(const vial_type_t* a, const vial_type_t* b, vial_type_t* out);   // NumPy's, shapes only
void vx_join_stage(vial_type_t* out, const vial_type_t* a, const vial_type_t* b);  // stage, pass, missing

// The read over the trajectory in which a value is complete: 0 static, a frame value pass 1 at the earliest
static inline uint8_t vx_pass(const vial_type_t* t) { return t->stage != VIAL_STAGE_STATIC && t->pass == 0 ? 1 : t->pass; }
bool vx_extent_join(vial_extent_t a, vial_extent_t b, vial_extent_t* out);
void vx_check_pairing(infer_t* x, vial_node_id at, const vial_type_t* a, const vial_type_t* b);   // W0605

void vx_unit_join(infer_t* x, vial_node_id at, const char* verb, const vial_type_t* a, const vial_type_t* b, vial_type_t* out);
uint8_t vx_promote(const vial_type_t* a, const vial_type_t* b);    // NumPy 2's dtype of a and b together

// infer_geom.c

// A parameter of a built-in: its name, its kind (how much of an argument one call takes) and whether it must be given
typedef struct vx_param_t {
    const char* name;
    uint8_t     kind;       // VX_KIND_*
    bool        required;
} vx_param_t;

enum { VX_KIND_WHOLE = 0, VX_KIND_SET, VX_KIND_GROUP, VX_KIND_POINT, VX_KIND_POINTS, VX_KIND_FLOAT, VX_KIND_INT, VX_KIND_BOOL };

// What one call takes of an argument under a kind (core) and what it maps over (loop's dims); reports E0302, E0313
bool vx_split(infer_t* x, vial_node_id node, const vial_type_t* a, uint8_t kind, const char* fname, const char* pname,
              vial_type_t* core, vial_type_t* loop);

// The statements of a block; true when it surely returned
bool vx_check_block(infer_t* x, vial_node_id block);

// Binds a call's arguments to parameters, by position and then by keyword, and types each; E0202 and E0204
bool vx_bind(infer_t* x, vial_node_id call, const char* fname, const vx_param_t* params, uint32_t np, vial_node_id* args);

// A call to a geometry built-in, with kinds converting and mapping its arguments; false for other names
bool vx_call_geom(infer_t* x, vial_node_id call, const vial_builtin_t* b, vial_type_t* out);

// .xyz, .x, .y, .z on a set, and .com on a set or points; false when it is not one of them
bool vx_attr_geom(infer_t* x, vial_node_id attr, vial_node_id on, const vial_type_t* a, str_t name, vial_type_t* out);

// infer_defs.c

// A call to a def of the script or the prelude: its arguments bound to its parameters, split by their kinds, and its
// body checked for those types, once per distinct set of them
vial_type_t vx_call_def(infer_t* x, vial_node_id call, const vial_tree_t* t, const vial_names_t* n, str_t src, bool prelude, vial_node_id def);

// A call to a def given as a value, with the types of its first num arguments: fold's step
vial_type_t vx_call_def_types(infer_t* x, vial_node_id site, vial_node_id def, const vial_type_t* types, const vial_value_t* values, uint32_t num,
                              const uint32_t* ssas, uint32_t* out_ssa);

// shown() and the drawing primitives, until step 9 types drawings fully; false for other names
bool vx_call_draw(infer_t* x, vial_node_id call, const vial_builtin_t* b, vial_type_t* out);

// A symbol's type and value in the frames of the defs being checked; false when no frame binds it
bool vx_frame_get(const infer_t* x, uint32_t sym, vial_type_t* t, vial_value_t* v);
void vx_frame_set(infer_t* x, uint32_t sym, const vial_type_t* t, const vial_value_t* v, vial_node_id prov);
vial_node_id vx_frame_prov(const infer_t* x, uint32_t sym);

// The value a return gives, joined into the frame's result
void vx_frame_return(infer_t* x, vial_node_id at, const vial_type_t* t);

// infer_comp.c: comprehensions, loops, fold, frame() and time()
vial_type_t vx_listcomp(infer_t* x, vial_node_id id);
bool        vx_check_for(infer_t* x, vial_node_id s);
bool        vx_call_frames(infer_t* x, vial_node_id call, const vial_builtin_t* b, vial_type_t* out);

// infer_math.c: reductions, math and arrays; false for other names
bool vx_call_math(infer_t* x, vial_node_id call, const vial_builtin_t* b, vial_type_t* out);

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

// ### LOWERING (emit.c) ###
// The check emits the array IR as it goes: each checked node gets its value in the IR, a constant when the check knows
// its value. What cannot be lowered yet is recorded in the IR's failures, not reported: the program still checks.

typedef uint32_t vial_ssa_t;

// A built-in call's arguments by parameter, as vx_bind bound them
typedef struct vx_bound_t {
    const vx_param_t*   params;
    uint32_t            np;
    vial_node_id        args[12];
} vx_bound_t;

void        vx_record_bind(infer_t* x, vial_node_id call, const vx_param_t* params, uint32_t np, const vial_node_id* args);
vial_ssa_t  vx_ssa(const infer_t* x, vial_node_id id);
void        vx_set_ssa(infer_t* x, vial_node_id id, vial_ssa_t v);
void        vx_emit(infer_t* x, vial_node_id id, const vial_type_t* t);    // the node's value, once it is typed
void        vx_lower_fail(infer_t* x, vial_node_id id, const char* what);
md_vial_span_t vx_span(const infer_t* x, vial_node_id id);                  // in the script, for an op

// A symbol's value: in the frames of the defs being checked, then the top level
vial_ssa_t  vx_sym_ssa(const infer_t* x, uint32_t sym);
void        vx_frame_set_ssa(infer_t* x, uint32_t sym, vial_ssa_t v);
void        vx_sizes_ssa(infer_t* x, vial_types_t* out, size_t num_nodes, size_t num_syms);   // arrays sized for lowering

// The i-th of what a tuple target unpacks: an element of a tuple, or a position along the first dimension
vial_ssa_t  vx_unpack(infer_t* x, vial_node_id target, vial_ssa_t value, const vial_type_t* t, uint32_t i, const vial_type_t* et);

// Ends the top level: its yield, and the published names
void        vx_emit_finish(infer_t* x);

// A return in a def's body: a yield of its value (or None) that ends the body's function
void        vx_emit_return(infer_t* x, vial_node_id value, vial_node_id at);

// A call of a def's specialization fn, or of a nested def's body region: its arguments converted to the parameters' kinds, called once, or mapped over
// the dimensions the kinds leave
vial_ssa_t  vx_emit_call(infer_t* x, vial_node_id call, uint32_t fn, uint32_t body, const uint8_t* kinds, const uint32_t* argv, const vial_type_t* atypes,
                         const vial_type_t* loops, const vial_type_t* cores, uint32_t np, const vial_type_t* result);

// A comprehension's map, while its element is checked: the region, and what it maps over
typedef struct vx_comp_t {
    uint32_t    region;
    uint32_t    saved;          // the region to return to
    uint32_t    operands[8];
    uint32_t    num;
    bool        ok;
} vx_comp_t;

// Starts the map of a comprehension over it, whose elements are elem, bound to target: emits the iterable's value in
// the enclosing region and enters the region. False when the iterable is not lowered.
bool        vx_comp_begin(infer_t* x, vial_node_id it, vial_node_id target, const vial_type_t* elem, vx_comp_t* c);
// Ends it: yields value (and keep, when filters are given; 0 otherwise) and emits the map, of type t, as id's value
void        vx_comp_end(infer_t* x, vial_node_id id, vx_comp_t* c, uint32_t value, uint32_t keep, const vial_type_t* t);

// What follows a failed static assert never runs: the body ends there, and its def is not evaluated
void        vx_emit_unreachable(infer_t* x);

// What failed to check and was reported makes the def it is in not evaluated, without a lowering failure of its own
void        vx_emit_unreachable_mark(infer_t* x);
