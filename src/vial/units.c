// Quantity literals: the unit names that may follow a number, as in `10 ns`, and the unit strings, as in
// `1.5 'kJ/mol'`. The names are a closed table, so whether a name after a number is a unit never depends on what
// the script binds: `ps` stays free as a variable, and `10 ps` is always ten picoseconds.

#include "vial_internal.h"

#include <core/md_common.h>
#include <core/md_unit.h>

#include <string.h>

typedef struct vial_unit_name_t {
    const char* name;   // as written after a number
    const char* unit;   // as md_unit_parse reads it
} vial_unit_name_t;

static const vial_unit_name_t vial_unit_names[] = {
    {u8"Å", u8"Å"}, {"angstrom", u8"Å"}, {"nm", "nm"}, {"pm", "pm"},
    {"fs", "fs"}, {"ps", "ps"}, {"ns", "ns"}, {u8"µs", u8"μs"}, {u8"μs", u8"μs"}, {"us", u8"μs"},
    {"deg", "deg"}, {"rad", "rad"},
    {"K", "K"}, {"e", "e"}, {"u", "u"}, {"kJ", "kJ"}, {"kcal", "kcal"},
};

bool vial_unit_name(str_t name, md_unit_t* unit) {
    for (size_t i = 0; i < ARRAY_SIZE(vial_unit_names); ++i) {
        if (str_eq(name, str_from_cstr(vial_unit_names[i].name))) {
            if (unit) md_unit_parse(unit, str_from_cstr(vial_unit_names[i].unit));
            return true;
        }
    }
    return false;
}

bool vial_unit_string(str_t text, md_unit_t* unit) {
    // md_unit spells micro with the Greek mu; accept the micro sign, which keyboards type, as well
    char buf[128];
    size_t n = 0;
    for (size_t i = 0; i < text.len; ++i) {
        if (n + 2 >= sizeof(buf)) return false;
        if ((unsigned char)text.ptr[i] == 0xC2 && i + 1 < text.len && (unsigned char)text.ptr[i + 1] == 0xB5) {
            buf[n++] = (char)0xCE; buf[n++] = (char)0xBC; ++i;
        } else {
            buf[n++] = text.ptr[i];
        }
    }
    const str_t s = str_trim((str_t){buf, n});
    if (str_empty(s)) return false;
    md_unit_t u;
    if (!md_unit_parse(&u, s)) return false;
    if (unit) *unit = u;
    return true;
}
