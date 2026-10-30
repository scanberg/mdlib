#pragma once

// The checker, step 3 on: the type of every expression, and the value of every static one that can be known when
// compiling. Steps add to it: literals, operators, units and indexing now; selectors, built-ins, defs and stages
// later. An expression whose type cannot be known yet is typed VIAL_DT_ERROR without a diagnostic, so nothing built
// on it warns; the harness counts those as not checked.

#include "check.h"
#include "types.h"
#include "select.h"

typedef enum vial_value_kind_t {
    VIAL_VAL_NONE = 0,      // not known when compiling
    VIAL_VAL_BOOL,
    VIAL_VAL_INT,
    VIAL_VAL_FLOAT,
    VIAL_VAL_STR,
    VIAL_VAL_NULL,          // the value None
    VIAL_VAL_SETS,          // a set or an array of sets: sets
    VIAL_VAL_STRS,          // an array of strings: strs[0..n)
    VIAL_VAL_LABELS,        // the labels of the units of sets, written only when one is taken
    VIAL_VAL_SHAPE,         // shape(a): the extents of type, with what each runs over
} vial_value_kind_t;

typedef struct vial_value_t {
    uint8_t kind;           // vial_value_kind_t
    bool    b;
    int64_t i;
    double  f;              // in the unit of its type
    str_t   s;
    const vial_sets_t*  sets;
    const str_t*        strs;
    uint32_t            n;
    const struct vial_type_t* type;
} vial_value_t;

// What a dimension runs over, for labels, set indices, by= and W0605: the sets of an array, whose units it describes
typedef struct vial_desc_t {
    uint8_t             level;  // vial_level_t
    const vial_sets_t*  sets;   // the array's sets, when known
} vial_desc_t;

typedef struct vial_types_t {
    vial_type_t*    types;      // md_array by node id
    vial_value_t*   values;     // md_array by node id
    uint8_t*        typed;      // md_array by node id: 1 when types[] holds its type
    vial_type_t*    sym_types;  // md_array by symbol
    vial_value_t*   sym_values; // md_array by symbol
    vial_desc_t*    descs;      // md_array; an extent's desc indexes it, 0 is none
} vial_types_t;

// Types the top-level statements of c's tree in order. Runs after vial_check_names.
void vial_check_types(vial_check_t* c, vial_types_t* out);

// Whether a value is one the reference writes after ->: a number, a bool, a string or None
static inline bool vial_value_is_scalar(const vial_value_t* v) { return v->kind >= VIAL_VAL_BOOL && v->kind <= VIAL_VAL_NULL; }

// The value as the reference writes it: 442, True, 'MET1/A'. Returns the length written.
size_t vial_value_print(char* buf, size_t cap, const vial_value_t* v);
