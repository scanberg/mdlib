#pragma once

// The built-in names of VIAL: functions, selectors, predefined sets, properties and kinds.

#include "vial_internal.h"

enum {
    VIAL_BUILTIN_FUNC     = 1 << 0,     // called: f(x), or as a method or property, x.f
    VIAL_BUILTIN_SELECTOR = 1 << 1,     // a selector: called, or bare for all its units
    VIAL_BUILTIN_VALUE    = 1 << 2,     // a value written bare: a predefined set, or a selector
    VIAL_BUILTIN_PROPERTY = 1 << 3,     // after a dot: .T, .xyz, .x, .y, .z, .com
    VIAL_BUILTIN_KIND     = 1 << 4,     // a kind a parameter may be annotated with
};

typedef struct vial_builtin_t {
    const char* name;
    uint32_t    flags;
} vial_builtin_t;

const vial_builtin_t* vial_builtin_find(str_t name);
size_t                vial_builtin_count(void);
const vial_builtin_t* vial_builtin_at(size_t i);

// A name v1 had: the message E0701 gives, and the name to write instead
bool        vial_builtin_v1(str_t name, const char** message, const char** fix);

// A word that VIAL spells differently, such as centre for com; NULL when there is none
const char* vial_builtin_synonym(str_t word);
