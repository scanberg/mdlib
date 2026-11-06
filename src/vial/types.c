#include "types.h"

#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

vial_type_t vial_type_scalar(vial_dtype_t dtype) {
    vial_type_t t;
    memset(&t, 0, sizeof(t));
    t.dtype = (uint8_t)dtype;
    t.unit = md_unit_none();
    return t;
}

vial_type_t vial_type_error(void) {
    return vial_type_scalar(VIAL_DT_ERROR);
}

bool vial_type_push_dim(vial_type_t* t, vial_extent_t e) {
    if (t->ndim >= VIAL_MAX_DIMS) return false;
    t->dims[t->ndim++] = e;
    return true;
}

static const char* dtype_name(vial_dtype_t d) {
    switch (d) {
    case VIAL_DT_ERROR:   return "error";
    case VIAL_DT_BOOL:    return "bool";
    case VIAL_DT_INT:     return "int";
    case VIAL_DT_FLOAT:   return "float";
    case VIAL_DT_FLOAT64: return "float64";
    case VIAL_DT_STR:     return "str";
    case VIAL_DT_SET:     return "set";
    case VIAL_DT_RANGE:   return "range";
    case VIAL_DT_DRAWING: return "drawing";
    case VIAL_DT_NONE:    return "None";
    case VIAL_DT_TUPLE:   return "tuple";
    case VIAL_DT_DEF:     return "def";
    default:              return "?";
    }
}

typedef struct out_t { char* buf; size_t cap; size_t len; } out_t;

static void put(out_t* o, const char* s, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (o->len + 1 < o->cap) o->buf[o->len] = s[i];
        o->len += 1;
    }
    if (o->cap) o->buf[o->len < o->cap ? o->len : o->cap - 1] = '\0';
}
static void puts_(out_t* o, const char* s) { put(o, s, strlen(s)); }
static void putf(out_t* o, const char* fmt, ...) {
    char tmp[64];
    va_list args;
    va_start(args, fmt);
    int n = vsnprintf(tmp, sizeof(tmp), fmt, args);
    va_end(args);
    if (n > 0) put(o, tmp, (size_t)n < sizeof(tmp) ? (size_t)n : sizeof(tmp) - 1);
}

size_t vial_unit_print(char* buf, size_t cap, md_unit_t unit) {
    char raw[128];
    const size_t n = md_unit_print(raw, sizeof(raw), unit);
    out_t o = {buf, cap, 0};
    if (cap) buf[0] = '\0';
    // ⁰ to ⁹ as UTF-8 bytes, which every compiler keeps as they are
    static const char* sup[] = {"\xE2\x81\xB0", "\xC2\xB9", "\xC2\xB2", "\xC2\xB3", "\xE2\x81\xB4", "\xE2\x81\xB5", "\xE2\x81\xB6", "\xE2\x81\xB7", "\xE2\x81\xB8", "\xE2\x81\xB9"};
    for (size_t i = 0; i < n; ++i) {
        // ^k for a small positive power becomes a superscript; anything else stays as md_unit spells it
        if (raw[i] == '^' && i + 1 < n && raw[i + 1] >= '2' && raw[i + 1] <= '9' && (i + 2 >= n || raw[i + 2] < '0' || raw[i + 2] > '9')) {
            puts_(&o, sup[raw[i + 1] - '0']);
            i += 1;
            continue;
        }
        put(&o, raw + i, 1);
    }
    return o.len;
}

static void print_extent(out_t* o, const vial_extent_t* e) {
    switch (e->kind) {
    case VIAL_EXT_STATIC:  putf(o, "[%u]", e->n); break;
    case VIAL_EXT_DYNAMIC: if (e->offset) putf(o, "[?-%d]", e->offset); else puts_(o, "[?]"); break;
    case VIAL_EXT_RANGED:  putf(o, "[%u..%u]", e->n, e->m); break;
    default:               puts_(o, "[n]"); break;
    }
}

static void print_type(out_t* o, const vial_type_t* t) {
    if (t->dtype == VIAL_DT_TUPLE) {
        puts_(o, "(");
        for (uint32_t i = 0; i < t->num_elems; ++i) {
            if (i) puts_(o, ", ");
            print_type(o, &t->elems[i]);
        }
        puts_(o, ")");
        return;
    }
    puts_(o, dtype_name((vial_dtype_t)t->dtype));
    for (uint8_t i = 0; i < t->ndim; ++i) print_extent(o, &t->dims[i]);
    const bool has_unit = t->dtype == VIAL_DT_FLOAT || t->dtype == VIAL_DT_FLOAT64 || t->dtype == VIAL_DT_INT;
    if (has_unit) {
        if (t->unit_kind == VIAL_UNIT_UNKNOWN) {
            puts_(o, " ?");
        } else if (t->unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(t->unit)) {
            char u[128];
            vial_unit_print(u, sizeof(u), t->unit);
            puts_(o, " ");
            puts_(o, u);
        }
    }
    if (t->flags & VIAL_TYPE_MISSING) puts_(o, " | missing");
}

size_t vial_type_print(char* buf, size_t cap, const vial_type_t* t) {
    out_t o = {buf, cap, 0};
    if (cap) buf[0] = '\0';
    print_type(&o, t);
    return o.len;
}

const char* vial_stage_str(const vial_type_t* t) {
    switch (t->stage) {
    case VIAL_STAGE_STATIC: return "static";
    case VIAL_STAGE_FRAME:  return "frame";
    default:                return (t->flags & VIAL_TYPE_PER_FRAME) ? "trajectory, per frame" : "trajectory";
    }
}
