// Quantity literals: the unit names that may follow a number, as in `10 ns`, the units divided by them, as in
// `1.5 kJ/mol`, and the unit strings, as in `1.5 'kJ/mol'`. The names are a closed table, so whether a name after a
// number is a unit never depends on what the script binds: `ps` stays free as a variable, and `10 ps` is always ten
// picoseconds. A divisor may be any one name md_unit knows (mol, s, K), written against the slash; the checker
// refuses one the script also binds (E0503), since `10 nm/ps` divided by the variable before compound units.

#include "vial_internal.h"

#include <core/md_common.h>
#include <core/md_unit.h>

#include <string.h>

typedef struct vial_unit_name_t {
    const char* name;   // as written after a number
    const char* unit;   // as md_unit_parse reads it
} vial_unit_name_t;

// Names outside ASCII are written as their UTF-8 bytes, which every compiler keeps as they are: \xC3\x85 is Å,
// \xC2\xB5 the micro sign and \xCE\xBC the Greek mu, both read as μ
static const vial_unit_name_t vial_unit_names[] = {
    {"\xC3\x85", "\xC3\x85"}, {"angstrom", "\xC3\x85"}, {"nm", "nm"}, {"pm", "pm"},
    {"fs", "fs"}, {"ps", "ps"}, {"ns", "ns"}, {"\xC2\xB5s", "\xCE\xBCs"}, {"\xCE\xBCs", "\xCE\xBCs"}, {"us", "\xCE\xBCs"},
    {"deg", "deg"}, {"rad", "rad"},
    {"K", "K"}, {"e", "e"}, {"u", "u"}, {"kJ", "kJ"}, {"kcal", "kcal"},
};

bool vial_unit_name(str_t name, md_unit_t* unit) {
    for (size_t i = 0; i < MD_COUNTOF(vial_unit_names); ++i) {
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

bool vial_unit_divisor(str_t name, md_unit_t* unit) {
    if (str_empty(name)) return false;
    for (size_t i = 0; i < name.len; ++i) {
        const unsigned char c = (unsigned char)name.ptr[i];
        if (!(c >= 0x80 || c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (i > 0 && c >= '0' && c <= '9'))) return false;
    }
    return vial_unit_name(name, unit) || vial_unit_string(name, unit);
}

bool vial_unit_written(str_t text, md_unit_t* unit) {
    size_t e = 0;
    while (e < text.len && text.ptr[e] != '/') ++e;
    md_unit_t u;
    if (!vial_unit_name((str_t){text.ptr, e}, &u)) return false;
    while (e < text.len) {
        const size_t b = e + 1;
        e = b;
        while (e < text.len && text.ptr[e] != '/') ++e;
        md_unit_t d;
        if (!vial_unit_divisor((str_t){text.ptr + b, e - b}, &d)) return false;
        u = md_unit_div(u, d);
    }
    if (unit) *unit = u;
    return true;
}
