#pragma once

// The types of VIAL, as the checker infers them and the reference prints them: a dtype, a shape whose dimensions say
// what they run over, a unit, a period, a stage and a pass, and whether the value may be missing.
//
//     float[433][3] Å              433 positions, static
//     set[?]                       a frame-stage array of sets, its length known only per frame
//     float Å | missing            a length that may be missing in some frames
//     float[4..14][3] Å            inside a mapped call: an extent that changes from call to call
//
// A type is a plain value, copied freely; tuple elements point into the checker's arena.

#include "vial_internal.h"

#include <core/md_unit.h>

#define VIAL_MAX_DIMS 6

typedef enum vial_dtype_t {
    VIAL_DT_ERROR = 0,      // a value whose expression failed: it absorbs, so what is built on it gives no diagnostic
    VIAL_DT_BOOL,
    VIAL_DT_INT,
    VIAL_DT_FLOAT,
    VIAL_DT_FLOAT64,        // time() and what is computed from it alone
    VIAL_DT_STR,
    VIAL_DT_SET,
    VIAL_DT_RANGE,
    VIAL_DT_DRAWING,
    VIAL_DT_NONE,           // the value None
    VIAL_DT_TUPLE,          // elements in elems
    VIAL_DT_DEF,            // a def, passed as an argument
    VIAL_DT_COUNT
} vial_dtype_t;

typedef enum vial_extent_kind_t {
    VIAL_EXT_STATIC = 0,    // n, known when compiling
    VIAL_EXT_DYNAMIC,       // ?, known only per frame; id is its identity, offset its related slice ?-offset
    VIAL_EXT_RANGED,        // n..m: inside a mapped call, it changes from call to call
    VIAL_EXT_UNKNOWN,       // static but not known, as without a system: printed as a letter
} vial_extent_kind_t;

typedef struct vial_extent_t {
    uint8_t  kind;          // vial_extent_kind_t
    int32_t  offset;        // VIAL_EXT_DYNAMIC: how much shorter than its identity a fixed-offset slice is
    uint32_t n;             // VIAL_EXT_STATIC: the extent; VIAL_EXT_RANGED: the least
    uint32_t m;             // VIAL_EXT_RANGED: the most
    uint32_t id;            // VIAL_EXT_DYNAMIC: the identity, the selection that produced it
    uint32_t desc;          // what the dimension runs over, an id the checker hands out; 0 when nothing is known
} vial_extent_t;

typedef enum vial_unit_kind_t {
    VIAL_UNIT_KNOWN = 0,    // unit holds it; md_unit_none() is a value without a unit
    VIAL_UNIT_OPEN,         // a bare number, or zeros() without unit=: it takes the unit it meets
    VIAL_UNIT_UNKNOWN,      // after a mismatch, or from an attribute without one: printed as ?
} vial_unit_kind_t;

typedef enum vial_stage_t {
    VIAL_STAGE_STATIC = 0,
    VIAL_STAGE_FRAME,
    VIAL_STAGE_TRAJECTORY,
} vial_stage_t;

enum {
    VIAL_TYPE_MISSING      = 1 << 0,    // may be missing: carries a mask
    VIAL_TYPE_PER_FRAME    = 1 << 1,    // a trajectory value that keeps the frame, as smooth() gives
    VIAL_TYPE_ONE_PARTICLE = 1 << 2,    // a set known to hold one particle, so it is a point
    VIAL_TYPE_EMPTY        = 1 << 3,    // a set known to be empty
    VIAL_TYPE_OVERLAP      = 1 << 4,    // an array of sets whose sets may share particles
    VIAL_TYPE_PERIODIC     = 1 << 5,    // period_beg, period_end and seam hold its period
    VIAL_TYPE_DRAWN        = 1 << 6,    // shown() attached a drawing
};

typedef struct vial_type_t {
    uint8_t         dtype;      // vial_dtype_t
    uint8_t         ndim;
    uint8_t         stage;      // vial_stage_t
    uint8_t         pass;       // the read over the trajectory in which it is complete; 0 when static
    uint8_t         unit_kind;  // vial_unit_kind_t
    uint16_t        flags;      // VIAL_TYPE_*
    vial_extent_t   dims[VIAL_MAX_DIMS];
    md_unit_t       unit;
    float           period_beg, period_end, seam;
    const struct vial_type_t* elems;    // VIAL_DT_TUPLE: its elements
    uint32_t        num_elems;
} vial_type_t;

// Scalars of each dtype, static, without a unit
vial_type_t vial_type_scalar(vial_dtype_t dtype);
vial_type_t vial_type_error(void);

static inline bool vial_type_is_error(const vial_type_t* t) { return t->dtype == VIAL_DT_ERROR; }

// Appends a dimension on the right; false when it would be the seventh (E0308)
bool vial_type_push_dim(vial_type_t* t, vial_extent_t e);
static inline vial_extent_t vial_extent_static(uint32_t n) { return (vial_extent_t){.kind = VIAL_EXT_STATIC, .n = n}; }
static inline vial_extent_t vial_extent_dynamic(uint32_t id) { return (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = id}; }

// The type as the reference prints it: dtype, extents, unit, and " | missing". Stage and pass are printed apart.
// Returns the length written, without the terminating zero; the text is cut at cap.
size_t vial_type_print(char* buf, size_t cap, const vial_type_t* t);

// The stage as the reference writes it after the type: "static", "frame", "trajectory" or "trajectory, per frame"
const char* vial_stage_str(const vial_type_t* t);

// A unit as the reference prints it: md_unit's spelling with small powers as superscripts, Å² rather than Å^2
size_t vial_unit_print(char* buf, size_t cap, md_unit_t unit);
