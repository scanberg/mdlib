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

#define VIAL_MAP_SCRATCH_DEPTH 8

typedef struct vial_interp_t {
    const vial_ir_t*            ir;
    vial_sys_t*                 sys;        // NULL without a system
    const struct md_system_t*   msys;
    float*                      masses;     // each particle's mass in msys, made the first time a kernel asks
    md_allocator_i*             alloc;      // what lives as long as the interpreter: static values
    md_allocator_i*             temp;       // what lives for one frame: an arena, reset between frames
    size_t                      max_elements;   // the most elements one value may have (R0105); 0 for the default
    vial_val_t*                 vals;       // md_array by ssa
    uint8_t*                    state;      // md_array by ssa: 0 not computed, 1 computed for this frame, 2 static
    const vial_frame_t*         frame;      // the frame being evaluated, or NULL
    vial_rt_error_t             error;      // the first error since the interpreter began
    bool                        failed;     // an op failed: what is computed from it is not computed
    vial_rt_error_t*            errors;     // md_array: every op that failed, once each, in the order they failed
    uint8_t*                    live;       // md_array by op: what the outputs need
    // at(x, frame=i) reads another frame, in a step before the frames are evaluated: the caller loads it
    bool                      (*load_frame)(void* user, int64_t index, vial_frame_t* out);
    void*                       load_user;
    // Over the trajectory (vial_interp_run): how each top-level op is evaluated, in which sweep over the frames, and
    // the per-frame values kept for what reads them as a series
    uint8_t*                    op_class;   // md_array by op: VIAL_OP_STATIC, VIAL_OP_PER_FRAME or VIAL_OP_ONCE
    uint8_t*                    op_sweep;   // md_array by op: the sweep a per-frame op runs in, after which a once op runs
    vial_val_t**                series;     // md_array by ssa: the value at each frame, or NULL when it is not kept
    double*                     times;      // [num_frames], time() of each frame
    int64_t                     num_frames;
    uint8_t                     sweep;      // the sweep being run
    vial_val_t*                 fold_state; // md_array by ssa of a fold's result: its state
    uint8_t*                    fold_has;
    // Progress and interruption, read and written from other threads: frames done of the frames to do, and a flag
    // that stops vial_interp_run between frames
    volatile uint64_t           frames_done;
    volatile uint64_t           frames_total;
    const volatile uint64_t*    stop;       // nonzero: stop, between frames or between a map's calls
    bool                        interrupted;
    // A map's calls are evaluated in a scratch arena per depth of maps in maps, reset after each call
    md_allocator_i*             map_scratch[VIAL_MAP_SCRATCH_DEPTH];
    uint32_t                    map_depth;
    bool                        no_batch;   // tests: a map's calls one by one, never batched for a built-in
} vial_interp_t;

enum { VIAL_OP_STATIC = 0, VIAL_OP_PER_FRAME = 1, VIAL_OP_ONCE = 2 };

bool vial_interp_init(vial_interp_t* in, const vial_ir_t* ir, vial_sys_t* sys, const struct md_system_t* msys, md_allocator_i* alloc);

// Computes only what the given values need from now on, as hover and at() do; init marks what the outputs need
void vial_interp_only(vial_interp_t* in, const vial_ssa_t* roots, size_t n);
void vial_interp_free(vial_interp_t* in);

// Computes the static values of the top level. An op that fails is recorded in in->errors, and what is computed from
// it is not computed; the rest is. False when something failed.
bool vial_interp_static(vial_interp_t* in);

// Computes the top level for one frame, likewise; values of earlier frames are dropped. Static values are computed
// first when they have not been. Trajectory values are not computed here.
bool vial_interp_frame(vial_interp_t* in, const vial_frame_t* frame);

// The value of v, as last computed; NULL when it was not
const vial_val_t* vial_interp_value(const vial_interp_t* in, vial_ssa_t v);

// Evaluates everything over num_frames frames, which in->load_frame gives: the static part, then sweeps over the
// frames in order, each computing the per-frame values that what was computed after the sweep before lets it, and
// after each sweep what reads whole series. False when something failed (in->errors).
bool vial_interp_run(vial_interp_t* in, int64_t num_frames);

// After vial_interp_run: the value of a per-frame value at frame f, or of a trajectory value (f ignored); NULL when it
// was not computed or not kept
const vial_val_t* vial_interp_at(const vial_interp_t* in, vial_ssa_t v, int64_t f);

// ### KERNELS (kernels.c) ###
// A kernel takes its arguments as runtime values and writes one result. It reports through the interpreter.

typedef struct vial_kctx_t {
    vial_interp_t*      in;
    const vial_ir_op_t* op;
    const vial_type_t*  type;       // the result's type
    md_allocator_i*     alloc;      // where the result's data goes
    int64_t             map_index;  // in a built-in mapped by its kinds: the call's position, for distinct=; else -1
} vial_kctx_t;

bool vial_rt_fail(vial_kctx_t* k, const char code[6], const char* fmt, ...);

// The default of vial_interp_t.max_elements: 2^28, a GB of floats
#define VIAL_MAX_ELEMENTS ((size_t)1 << 28)

// Whether a value of count elements fits the evaluation's limit; R0105 when it does not
bool vial_k_room(vial_kctx_t* k, size_t count);

// A value of dtype and shape, zeroed: false, with R0105, when it would exceed the limit
bool vial_k_make(vial_kctx_t* k, uint8_t dtype, uint8_t ndim, const uint32_t* dims, vial_val_t* out);
// The same, its data zeroed: for a kernel that writes only some elements (eye) or accumulates
bool vial_k_make_zero(vial_kctx_t* k, uint8_t dtype, uint8_t ndim, const uint32_t* dims, vial_val_t* out);

bool vial_k_unary(vial_kctx_t* k, uint16_t op, const vial_val_t* a, vial_val_t* out);
bool vial_k_not(vial_kctx_t* k, const vial_val_t* a, vial_val_t* out);
bool vial_k_binary(vial_kctx_t* k, uint16_t op, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
// One element of + - * / // % ** in VIAL_DT_FLOAT (a and b rounded to float first) or VIAL_DT_FLOAT64, exactly as
// vial_k_binary computes it. The check folds known values with it, so a value known when compiling is the value
// evaluation would give.
double vial_k_arith(uint16_t op, uint8_t dtype, double a, double b);

// A map's calls of closest_pair or farthest_pair at once, as a map whose calls are only that gives them: call c pairs
// group ia[c] of A with group ib[c] of B, each group a set's positions. out is float[count], made by the caller; its
// distances have the bits the calls one by one give, and a call with no pair is missing.
// VIAL_IR_PAIRS: red (vial_pairs_t) over axis of the distances between x [n][1] and y [m], or their test against r
bool vial_k_pairs(vial_kctx_t* k, uint32_t red, uint32_t axis, uint32_t test, const vial_val_t* x, const vial_val_t* y, const vial_val_t* r,
                  vial_val_t* out);

bool vial_k_pair_batch(vial_kctx_t* k, bool closest, const vial_sets_t* A, const vial_sets_t* B, const uint32_t* ia, const uint32_t* ib,
                       size_t count, vial_val_t* out);
bool vial_k_fma(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, const vial_val_t* c, vial_val_t* out);
bool vial_k_compare(vial_kctx_t* k, uint16_t cmp, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_logic(vial_kctx_t* k, uint16_t kw, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_setop(vial_kctx_t* k, uint16_t sub, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);
bool vial_k_convert(vial_kctx_t* k, const vial_val_t* a, uint8_t dtype, vial_val_t* out);
bool vial_k_index(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* items, uint32_t num, vial_val_t* out);
bool vial_k_stack(vial_kctx_t* k, const vial_val_t* items, uint32_t n, vial_val_t* out);
bool vial_k_kind(vial_kctx_t* k, uint16_t kind, const vial_val_t* a, vial_val_t* out);

bool vial_k_matmul(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, vial_val_t* out);

// A built-in, its arguments by parameter (NULL when not given) and the parameters' names
bool vial_k_builtin(vial_kctx_t* k, const char* name, const vial_val_t* const* args, const str_t* names, uint32_t n, vial_val_t* out);

// A value as dtype, converted elementwise; the same value when it already is
vial_val_t vial_val_as(const vial_val_t* a, uint8_t dtype, md_allocator_i* alloc);
