#pragma once

// The values the interpreter computes with, and the constants of the array IR.
//
// An array is dense and row-major, of one dtype, with up to six dimensions. Its type, unit and stage are the IR
// value's; the runtime value carries only what varies when it runs: the extents, the data and which elements are
// missing. Every array the interpreter holds is rectangular: a ragged extent exists only inside one call of a map
// region, where each call sees its own.
//
//     dtype     data                  notes
//     BOOL      uint8_t, 0 or 1
//     INT       int32_t
//     FLOAT     float
//     FLOAT64   double                time() and what is computed from it alone
//     STR       str_t                 the text is owned by the program or the evaluation's arena
//     SET       none                  sets: one set per element, in element order
//     RANGE     double[2]             a..b, both included; a scalar
//     DRAWING   vial_drawing_t*       drawings: one list of primitives per element; NULL draws nothing
//     NONE      none                  the value None
//     TUPLE     none                  elems
//
// Missing elements have their bit set in missing, which is NULL when none is; a missing element's data is zero.

#include "types.h"
#include "select.h"

#include <stdint.h>

// A drawing's primitives, as the reference's Drawing section lists them
typedef enum vial_prim_kind_t {
    VIAL_PRIM_MARKER = 0,   // p[0]
    VIAL_PRIM_LINE,         // p[0] to p[1], p[1] the nearest image of b
    VIAL_PRIM_ARROW,        // p[0] to p[1], the head at p[1]
    VIAL_PRIM_SPHERE,       // centre p[0], radius r
    VIAL_PRIM_ARC,          // the angle at p[0] from p[1] to p[2], both their nearest images
    VIAL_PRIM_DISK,         // centre p[0], normal p[1] (a unit vector), radius r
    VIAL_PRIM_LABEL,        // text at p[0]
    VIAL_PRIM_HIGHLIGHT,    // atoms
} vial_prim_kind_t;

typedef struct vial_prim_t {
    uint8_t         kind;       // vial_prim_kind_t
    bool            labelled;   // a label goes with it: line and arc with label=, and every label
    float           p[3][3];
    float           r;
    double          value;      // the label's number, with its unit; NaN when the label is text
    md_unit_t       unit;
    str_t           text;       // the label as text: the number with its unit, or the string
    const uint32_t* atoms;      // highlight: the particles, ascending
    uint32_t        num_atoms;
} vial_prim_t;

typedef struct vial_drawing_t {
    uint32_t            count;
    const vial_prim_t*  prims;
} vial_drawing_t;

// vial_val_t.flags
enum {
    VIAL_VAL_FRAMES = 1 << 0,   // the first dimension is the run's frames: a series, as trajectory operations take them
};

typedef struct vial_val_t {
    uint8_t             dtype;      // vial_dtype_t
    uint8_t             ndim;
    uint16_t            flags;
    uint32_t            dims[VIAL_MAX_DIMS];
    void*               data;
    uint64_t*           missing;    // a bit per element, or NULL
    const vial_sets_t*  sets;       // SET
    const struct vial_val_t* elems; // TUPLE
    uint32_t            num_elems;
    const vial_sets_t*  origin;     // positions and what is computed from them elementwise: the particles along their
                                    // first dimension of particles, as particles(x) gives them; NULL when not known
    const struct vial_val_t* drawn; // the drawing shown() attached, while drawing (vial_interp_t.draw); kept through
                                    // names, returns and calls, and dropped by any operation on the value
} vial_val_t;

static inline size_t vial_val_count(const vial_val_t* v) {
    size_t n = 1;
    for (uint8_t i = 0; i < v->ndim; ++i) n *= v->dims[i];
    return n;
}

// Bytes of one element's data: 0 for the dtypes that keep none
size_t vial_dtype_size(uint8_t dtype);

static inline bool vial_val_is_missing(const vial_val_t* v, size_t i) {
    return v->missing && (v->missing[i >> 6] >> (i & 63)) & 1;
}

// A value of dtype and shape, its data zeroed, from alloc
vial_val_t vial_val_make(uint8_t dtype, uint8_t ndim, const uint32_t* dims, md_allocator_i* alloc);
// The same, its data not zeroed: for a kernel that writes every element
vial_val_t vial_val_make_uninit(uint8_t dtype, uint8_t ndim, const uint32_t* dims, md_allocator_i* alloc);

// Scalars
vial_val_t vial_val_bool(bool b, md_allocator_i* alloc);
vial_val_t vial_val_int(int32_t i, md_allocator_i* alloc);
vial_val_t vial_val_float(float f, md_allocator_i* alloc);
vial_val_t vial_val_float64(double f, md_allocator_i* alloc);
vial_val_t vial_val_str(str_t s, md_allocator_i* alloc);    // copies s
vial_val_t vial_val_none(void);
vial_val_t vial_val_sets(const vial_sets_t* s);             // a set when s->single, else an array of them

// Marks element i missing, allocating the mask on first use
void vial_val_set_missing(vial_val_t* v, size_t i, md_allocator_i* alloc);

// The value as the reference writes it after ->, cut at cap: 442, 12.5, True, 'MET1/A', [1, 2, 3]. Arrays longer than
// max_items show their first items and "…". Returns the length written.
size_t vial_val_print(char* buf, size_t cap, const vial_val_t* v, uint32_t max_items);

// One drawing of the given primitives, copied into alloc; and a drawing of every primitive of a and b
const vial_drawing_t* vial_drawing_make(const vial_prim_t* prims, uint32_t count, md_allocator_i* alloc);
const vial_drawing_t* vial_drawing_concat(const vial_drawing_t* a, const vial_drawing_t* b, md_allocator_i* alloc);

// A copy of v whose data, mask, sets, elements and drawings are alloc's
vial_val_t vial_val_copy(const vial_val_t* v, md_allocator_i* alloc);

// The shortest text that reads back as the same float, as NumPy prints float32 (is32) and Python float64; a whole
// number keeps its point: 12.0
int vial_print_float(char* buf, size_t cap, double f, bool is32);

// Deep equality of dtype, shape, data and mask, bit for bit for floats
bool vial_val_equal(const vial_val_t* a, const vial_val_t* b);
