#pragma once

// The interpreter: runs the array IR one op at a time, each a call into a runtime kernel, much like NumPy. It is the
// reference every other backend must match bit for bit, and it keeps every intermediate value, so it also serves
// hover and drawing.
//
// Static values are computed once and kept. Frame values are computed for one frame at a time, from that frame's
// positions and unit cell; regions (maps, applied bodies, folds) run their ops once per call.

#include "ir.h"
#include "select.h"

#include <md_unitcell.h>

struct md_system_t;
struct vec3_t;

// One frame of a run, as the interpreter reads it
typedef struct vial_frame_t {
    const float*    xyz;        // [num_atoms][3], Å
    size_t          num_atoms;
    md_unitcell_t   cell;
    int64_t         index;      // the frame's position in the run
    double          time;       // the run's time axis at that frame, in its unit
} vial_frame_t;

// Runtime errors: an evaluation that cannot go on, with the code of the reference's diagnostics table
typedef struct vial_rt_error_t {
    char            code[6];    // "R0103"
    char            message[160];
    md_vial_span_t  span;
} vial_rt_error_t;

typedef struct vial_interp_t {
    const vial_ir_t*            ir;
    vial_sys_t*                 sys;        // NULL without a system
    const struct md_system_t*   msys;
    md_allocator_i*             alloc;      // what lives as long as the interpreter: static values
    md_allocator_i*             temp;       // what lives for one frame: an arena, reset between frames
    vial_val_t*                 vals;       // md_array by ssa
    uint8_t*                    state;      // md_array by ssa: 0 not computed, 1 computed for this frame, 2 static
    const vial_frame_t*         frame;      // the frame being evaluated, or NULL
    vial_rt_error_t             error;
    bool                        failed;
} vial_interp_t;

bool vial_interp_init(vial_interp_t* in, const vial_ir_t* ir, vial_sys_t* sys, const struct md_system_t* msys, md_allocator_i* alloc);
void vial_interp_free(vial_interp_t* in);

// Computes the static values of the top level. False, with in->error, when one fails.
bool vial_interp_static(vial_interp_t* in);

// Computes the top level for one frame; values of earlier frames are dropped. Static values are computed first when
// they have not been. Trajectory values are not computed here.
bool vial_interp_frame(vial_interp_t* in, const vial_frame_t* frame);

// The value of v, as last computed; NULL when it was not
const vial_val_t* vial_interp_value(const vial_interp_t* in, vial_ssa_t v);

// ### KERNELS (kernels.c) ###
// A kernel takes its arguments as runtime values and writes one result. It reports through the interpreter.

typedef struct vial_kctx_t {
    vial_interp_t*      in;
    const vial_ir_op_t* op;
    const vial_type_t*  type;       // the result's type
    md_allocator_i*     alloc;      // where the result's data goes
} vial_kctx_t;

bool vial_rt_fail(vial_kctx_t* k, const char code[6], const char* fmt, ...);

bool vial_k_unary(vial_kctx_t* k, uint16_t op, const vial_val_t* a, vial_val_t* out);
bool vial_k_not(vial_kctx_t* k, const vial_val_t* a, vial_val_t* out);
bool vial_k_binary(vial_kctx_t* k, uint16_t op, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_fma(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, const vial_val_t* c, vial_val_t* out);
bool vial_k_compare(vial_kctx_t* k, uint16_t cmp, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_logic(vial_kctx_t* k, uint16_t kw, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_setop(vial_kctx_t* k, uint16_t sub, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_convert(vial_kctx_t* k, const vial_val_t* a, uint8_t dtype, vial_val_t* out);
bool vial_k_index(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* items, uint32_t num, vial_val_t* out);
bool vial_k_stack(vial_kctx_t* k, const vial_val_t* items, uint32_t n, vial_val_t* out);
bool vial_k_kind(vial_kctx_t* k, uint16_t kind, const vial_val_t* a, vial_val_t* out);

// A built-in, its arguments by parameter (NULL when not given) and the parameters' names
bool vial_k_builtin(vial_kctx_t* k, const char* name, const vial_val_t* const* args, const str_t* names, uint32_t n, vial_val_t* out);

// A value as dtype, converted elementwise; the same value when it already is
vial_val_t vial_val_as(const vial_val_t* a, uint8_t dtype, md_allocator_i* alloc);
