#pragma once

// What is added up over the frames rather than computed at one: the volumes, sdf() and density() (volume.h), and the
// radial distribution function, rdf(). Each frame adds to one running sum, from any thread and in any order, and the
// value is made from the sum once, after its sweep. The sums are of integers, or kept per frame and added in frame
// order at the end, so a value is the same bits whichever thread added which frame, as every value of the interpreter is.

#include "volume.h"

#include <math.h>

// ### ARGUMENTS ###
// An op's arguments as the accumulators take them: by parameter name when the check bound them, else by position

typedef struct vial_args_t {
    const vial_val_t* const* v;
    const md_str_t*          names;
    uint32_t                 n;
} vial_args_t;

// The argument of parameter name, else the pos-th positional one; its index in the op's arguments, -1 when not given
static inline int vial_arg_index(const vial_args_t* a, const char* name, uint32_t pos) {
    bool bound = false;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (!a->names || !a->names[i].len) continue;
        if (md_str_eq(a->names[i], md_str_from_cstr(name))) return (int)i;
        bound = true;
    }
    if (bound) return -1;
    uint32_t p = 0;
    for (uint32_t i = 0; i < a->n; ++i) {
        if (a->names && a->names[i].len) continue;
        if (p++ == pos) return (int)i;
    }
    return -1;
}

static inline const vial_val_t* vial_arg(const vial_args_t* a, const char* name, uint32_t pos) {
    const int i = vial_arg_index(a, name, pos);
    return i < 0 ? NULL : a->v[i];
}

static inline bool vial_arg_none(const vial_val_t* v) { return !v || v->dtype == VIAL_DT_NONE; }

static inline double vial_arg_num(const vial_val_t* v, size_t i) {
    switch (v->dtype) {
    case VIAL_DT_BOOL:    return ((const uint8_t*)v->data)[i];
    case VIAL_DT_INT:     return ((const int32_t*)v->data)[i];
    case VIAL_DT_FLOAT:   return ((const float*)v->data)[i];
    case VIAL_DT_FLOAT64: return ((const double*)v->data)[i];
    default:              return NAN;
    }
}

static inline bool vial_arg_is_num(const vial_val_t* v) {
    return v && v->data && (v->dtype == VIAL_DT_INT || v->dtype == VIAL_DT_FLOAT || v->dtype == VIAL_DT_FLOAT64 || v->dtype == VIAL_DT_BOOL);
}

// ### RDF ###
// rdf(a, b, cutoff, bins=None): the radial distribution function of the points of b about the points of a, normalized
// as gmx rdf does by default (-norm rdf). Each frame counts the pairs of a point of a and a point of b whose distance
// under the minimum image is in [lo, hi), bin by bin, by a spatial search, so no matrix of every pair is made; a
// particle is never paired with itself. Bin k is
//
//     g_k = n_k / (P · 4/3 π (r_{k+1}³ - r_k³)),    P = Σ over frames of (N_a N_b - N_self) / V
//
// the pairs counted over the pairs that the same points spread evenly over each frame's cell, of volume V, would put in
// the shell: 1 at every distance for points that do not see each other. A frame without a cell, or whose cell is less
// than twice the cutoff across, is left out, with a note (R0101), since the minimum image would miscount it. Without a
// frame left, every bin is missing.
//
// The counts are kept per frame while they fit in what is left of the interpreter's frame_budget, handed out in program
// order as v1 hands its distributions room, so an rdf over chosen frames is summed from them (vial_rdf_over); one that
// does not fit keeps only its sum, and is made over other frames by reading them again, as a volume is.

#define VIAL_RDF_BIN_WIDTH      0.02                    // Å, by default: GROMACS's -bin 0.002 nm
#define VIAL_RDF_FRAME_BUDGET   ((size_t)256 << 20)    // bytes the rdfs keep per frame, by default

typedef struct vial_rdf_t {
    uint32_t    bins;
    double      lo, hi;         // Å
    int64_t     num_frames;
    uint64_t*   sum;            // [bins]: the pairs counted, over the frames added
    uint32_t*   frame_counts;   // [num_frames][bins] when kept per frame, else NULL
    double*     frame_pairs;    // [num_frames]: (N_a N_b - N_self) / V of each frame added, Å⁻³; NaN for one not added
    md_mutex_t  lock;
    bool        have_lock;
} vial_rdf_t;

// The rdf's part of the accumulators below
bool vial_rdf_begin(vial_kctx_t* k, vial_rdf_t* rdf, const vial_args_t* a);
bool vial_rdf_add(vial_kctx_t* k, vial_rdf_t* rdf, const vial_args_t* a);
bool vial_rdf_end(vial_kctx_t* k, vial_rdf_t* rdf, vial_val_t* out);
bool vial_rdf_like(vial_rdf_t* out, const vial_rdf_t* rdf);
void vial_rdf_free(vial_rdf_t* rdf);

// The range of the distances an rdf bins, from its cutoff argument in Å: (0, c) for a length, (lo, hi) for lo..hi or a
// tuple. scale converts the argument's unit to Å, and scale2 that of a tuple's second element.
bool vial_rdf_range(const vial_val_t* cutoff, double scale, double scale2, double* lo, double* hi);

// The rdf of chosen frames, those where mask[f] is set, from the counts kept per frame, into out (bins floats, NaN where
// missing). False when the counts are not kept per frame.
bool vial_rdf_over(const vial_rdf_t* rdf, const uint8_t* mask, float* out);

// ### ACCUMULATORS ###

typedef enum vial_accum_kind_t {
    VIAL_ACCUM_VOLUME = 1,
    VIAL_ACCUM_RDF,
} vial_accum_kind_t;

typedef struct vial_accum_t {
    uint8_t             kind;   // vial_accum_kind_t
    union {
        vial_volume_t   vol;
        vial_rdf_t      rdf;
    };
} vial_accum_t;

// Whether a built-in is added up over the frames rather than computed at one
bool vial_accum_builtin(const char* name);

// Sets an accumulator up from its op's static arguments (k->op), as vial_volume_begin does a volume, with what it adds
// frames to zeroed. False with an R diagnostic through k.
bool vial_accum_begin(vial_kctx_t* k, vial_accum_t* acc, const char* name, const vial_val_t* const* args, const md_str_t* names, uint32_t n,
                      const vial_frame_t* first, md_allocator_i* alloc);

// Adds the frame k->in->frame holds, from the op's arguments at that frame. Any thread, any order.
bool vial_accum_add(vial_kctx_t* k, vial_accum_t* acc, const vial_val_t* const* args, const md_str_t* names, uint32_t n);

// The value over what was added. Once.
bool vial_accum_end(vial_kctx_t* k, vial_accum_t* acc, vial_val_t* out);

// A copy of acc set up alike, nothing added: for the same value over other frames
bool vial_accum_like(vial_accum_t* out, const vial_accum_t* acc, md_allocator_i* alloc);

void vial_accum_free(vial_accum_t* acc);

// The number of elements of its value, known once it is set up
size_t vial_accum_count(const vial_accum_t* acc);

static inline const vial_volume_t* vial_accum_volume(const vial_accum_t* acc) {
    return acc && acc->kind == VIAL_ACCUM_VOLUME ? &acc->vol : NULL;
}

static inline const vial_rdf_t* vial_accum_rdf(const vial_accum_t* acc) {
    return acc && acc->kind == VIAL_ACCUM_RDF ? &acc->rdf : NULL;
}
