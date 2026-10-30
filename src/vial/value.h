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
//     DRAWING   vial_prim_t           drawings: one list of primitives per element
//     NONE      none                  the value None
//     TUPLE     none                  elems
//
// Missing elements have their bit set in missing, which is NULL when none is; a missing element's data is zero.

#include "types.h"
#include "select.h"

#include <stdint.h>

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

// Deep equality of dtype, shape, data and mask, bit for bit for floats
bool vial_val_equal(const vial_val_t* a, const vial_val_t* b);
