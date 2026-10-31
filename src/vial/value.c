#include "value.h"

#include <core/md_allocator.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t vial_dtype_size(uint8_t dtype) {
    switch (dtype) {
    case VIAL_DT_BOOL:    return 1;
    case VIAL_DT_INT:     return 4;
    case VIAL_DT_FLOAT:   return 4;
    case VIAL_DT_FLOAT64: return 8;
    case VIAL_DT_STR:     return sizeof(str_t);
    case VIAL_DT_RANGE:   return 16;
    default:              return 0;
    }
}

vial_val_t vial_val_make(uint8_t dtype, uint8_t ndim, const uint32_t* dims, md_allocator_i* alloc) {
    vial_val_t v = {.dtype = dtype, .ndim = ndim};
    if (ndim) memcpy(v.dims, dims, ndim * sizeof(uint32_t));
    const size_t bytes = vial_val_count(&v) * vial_dtype_size(dtype);
    if (bytes) {
        v.data = md_alloc(alloc, bytes);
        memset(v.data, 0, bytes);
    }
    return v;
}

vial_val_t vial_val_bool(bool b, md_allocator_i* alloc) {
    vial_val_t v = vial_val_make(VIAL_DT_BOOL, 0, NULL, alloc);
    *(uint8_t*)v.data = b;
    return v;
}

vial_val_t vial_val_int(int32_t i, md_allocator_i* alloc) {
    vial_val_t v = vial_val_make(VIAL_DT_INT, 0, NULL, alloc);
    *(int32_t*)v.data = i;
    return v;
}

vial_val_t vial_val_float(float f, md_allocator_i* alloc) {
    vial_val_t v = vial_val_make(VIAL_DT_FLOAT, 0, NULL, alloc);
    *(float*)v.data = f;
    return v;
}

vial_val_t vial_val_float64(double f, md_allocator_i* alloc) {
    vial_val_t v = vial_val_make(VIAL_DT_FLOAT64, 0, NULL, alloc);
    *(double*)v.data = f;
    return v;
}

vial_val_t vial_val_str(str_t s, md_allocator_i* alloc) {
    vial_val_t v = vial_val_make(VIAL_DT_STR, 0, NULL, alloc);
    *(str_t*)v.data = str_copy(s, alloc);
    return v;
}

vial_val_t vial_val_none(void) { return (vial_val_t){.dtype = VIAL_DT_NONE}; }

vial_val_t vial_val_sets(const vial_sets_t* s) {
    vial_val_t v = {.dtype = VIAL_DT_SET, .sets = s};
    if (!s->single) { v.ndim = 1; v.dims[0] = s->count; }
    return v;
}

void vial_val_set_missing(vial_val_t* v, size_t i, md_allocator_i* alloc) {
    if (!v->missing) {
        const size_t words = (vial_val_count(v) + 63) / 64;
        v->missing = md_alloc(alloc, words * sizeof(uint64_t));
        memset(v->missing, 0, words * sizeof(uint64_t));
    }
    v->missing[i >> 6] |= 1ull << (i & 63);
}

int vial_print_float(char* buf, size_t cap, double f, bool is32) {
    if (isnan(f)) return snprintf(buf, cap, "nan");
    if (isinf(f)) return snprintf(buf, cap, f < 0 ? "-inf" : "inf");
    char tmp[40];
    const int maxp = is32 ? 9 : 17;
    int p = 1;
    for (; p <= maxp; ++p) {
        snprintf(tmp, sizeof(tmp), "%.*g", p, f);
        const double back = strtod(tmp, NULL);
        if (is32 ? (float)back == (float)f : back == f) break;
    }
    // Positional, as Python's repr, unless the exponent is below -4 or 16 and above: 40.0, not 4e+01
    if (strpbrk(tmp, "eE") && f != 0.0) {
        const int e = (int)floor(log10(fabs(f)));
        if (e >= -4 && e < 16 && e + 1 > p) snprintf(tmp, sizeof(tmp), "%.*g", e + 1, f);
    }
    // A whole number keeps a point, so a float never reads as an int: 12.0
    if (!strpbrk(tmp, ".eEn")) {
        const size_t len = strlen(tmp);
        if (len + 2 < sizeof(tmp)) { tmp[len] = '.'; tmp[len + 1] = '0'; tmp[len + 2] = 0; }
    }
    return snprintf(buf, cap, "%s", tmp);
}

static int print_elem(char* buf, size_t cap, const vial_val_t* v, size_t i) {
    if (vial_val_is_missing(v, i)) return snprintf(buf, cap, "missing");
    switch (v->dtype) {
    case VIAL_DT_BOOL:    return snprintf(buf, cap, "%s", ((const uint8_t*)v->data)[i] ? "True" : "False");
    case VIAL_DT_INT:     return snprintf(buf, cap, "%d", ((const int32_t*)v->data)[i]);
    case VIAL_DT_FLOAT:   return vial_print_float(buf, cap, ((const float*)v->data)[i], true);
    case VIAL_DT_FLOAT64: return vial_print_float(buf, cap, ((const double*)v->data)[i], false);
    case VIAL_DT_STR: {
        const str_t s = ((const str_t*)v->data)[i];
        return snprintf(buf, cap, "'%.*s'", (int)s.len, s.ptr);
    }
    case VIAL_DT_RANGE: {
        const double* r = (const double*)v->data + 2 * i;
        return snprintf(buf, cap, "%g..%g", r[0], r[1]);
    }
    case VIAL_DT_SET:     return v->sets ? snprintf(buf, cap, "set of %u", vial_sets_size(v->sets, (uint32_t)i)) : snprintf(buf, cap, "set");
    case VIAL_DT_NONE:    return snprintf(buf, cap, "None");
    case VIAL_DT_DRAWING: return snprintf(buf, cap, "drawing");
    default:              return snprintf(buf, cap, "?");
    }
}

typedef struct printer_t {
    char*  buf;
    size_t cap;
    size_t len;
    uint32_t max_items;
} printer_t;

static void put(printer_t* p, const char* s) {
    const size_t n = strlen(s);
    for (size_t i = 0; i < n && p->len + 1 < p->cap; ++i) p->buf[p->len++] = s[i];
    if (p->cap) p->buf[p->len < p->cap ? p->len : p->cap - 1] = 0;
}

static void print_dims(printer_t* p, const vial_val_t* v, uint8_t d, size_t base) {
    if (d == v->ndim) {
        char tmp[96];
        print_elem(tmp, sizeof(tmp), v, base);
        put(p, tmp);
        return;
    }
    size_t stride = 1;
    for (uint8_t k = d + 1; k < v->ndim; ++k) stride *= v->dims[k];
    put(p, "[");
    const uint32_t n = v->dims[d];
    const uint32_t shown = p->max_items && n > p->max_items ? p->max_items : n;
    for (uint32_t i = 0; i < shown; ++i) {
        if (i) put(p, ", ");
        print_dims(p, v, d + 1, base + i * stride);
    }
    if (shown < n) put(p, ", …");
    put(p, "]");
}

size_t vial_val_print(char* buf, size_t cap, const vial_val_t* v, uint32_t max_items) {
    printer_t p = {.buf = buf, .cap = cap, .max_items = max_items};
    if (cap) buf[0] = 0;
    if (v->dtype == VIAL_DT_TUPLE) {
        put(&p, "(");
        for (uint32_t i = 0; i < v->num_elems; ++i) {
            if (i) put(&p, ", ");
            char tmp[256];
            vial_val_print(tmp, sizeof(tmp), &v->elems[i], max_items);
            put(&p, tmp);
        }
        put(&p, ")");
        return p.len;
    }
    print_dims(&p, v, 0, 0);
    return p.len;
}

bool vial_val_equal(const vial_val_t* a, const vial_val_t* b) {
    if (a->dtype != b->dtype || a->ndim != b->ndim) return false;
    for (uint8_t i = 0; i < a->ndim; ++i) if (a->dims[i] != b->dims[i]) return false;
    if (a->dtype == VIAL_DT_TUPLE) {
        if (a->num_elems != b->num_elems) return false;
        for (uint32_t i = 0; i < a->num_elems; ++i) if (!vial_val_equal(&a->elems[i], &b->elems[i])) return false;
        return true;
    }
    const size_t n = vial_val_count(a);
    for (size_t i = 0; i < n; ++i) if (vial_val_is_missing(a, i) != vial_val_is_missing(b, i)) return false;
    if (a->dtype == VIAL_DT_SET) {
        if (!a->sets || !b->sets) return a->sets == b->sets;
        for (uint32_t i = 0; i < (uint32_t)n; ++i) {
            const uint32_t na = vial_sets_size(a->sets, i);
            if (na != vial_sets_size(b->sets, i)) return false;
            if (na && memcmp(vial_sets_items(a->sets, i), vial_sets_items(b->sets, i), na * sizeof(uint32_t))) return false;
        }
        return true;
    }
    if (a->dtype == VIAL_DT_STR) {
        for (size_t i = 0; i < n; ++i) if (!str_eq(((const str_t*)a->data)[i], ((const str_t*)b->data)[i])) return false;
        return true;
    }
    const size_t bytes = n * vial_dtype_size(a->dtype);
    return !bytes || memcmp(a->data, b->data, bytes) == 0;
}

static vial_sets_t* sets_copy(const vial_sets_t* s, md_allocator_i* alloc) {
    if (!s) return NULL;
    vial_sets_t* c = md_alloc(alloc, sizeof(vial_sets_t));
    *c = *s;
    const uint32_t total = s->offset[s->count];
    c->offset = md_alloc(alloc, (s->count + 1) * sizeof(uint32_t));
    memcpy(c->offset, s->offset, (s->count + 1) * sizeof(uint32_t));
    c->idx = md_alloc(alloc, (total ? total : 1) * sizeof(uint32_t));
    if (total) memcpy(c->idx, s->idx, total * sizeof(uint32_t));
    if (s->unit) {
        c->unit = md_alloc(alloc, (s->count ? s->count : 1) * sizeof(uint32_t));
        memcpy(c->unit, s->unit, s->count * sizeof(uint32_t));
    }
    c->root = NULL;
    return c;
}

vial_val_t vial_val_copy(const vial_val_t* v, md_allocator_i* alloc) {
    vial_val_t c = *v;
    const size_t n = vial_val_count(v);
    const size_t bytes = n * vial_dtype_size(v->dtype);
    if (v->data && bytes) {
        c.data = md_alloc(alloc, bytes);
        memcpy(c.data, v->data, bytes);
        if (v->dtype == VIAL_DT_STR) for (size_t i = 0; i < n; ++i) ((str_t*)c.data)[i] = str_copy(((const str_t*)v->data)[i], alloc);
    }
    if (v->missing) {
        const size_t words = (n + 63) / 64;
        c.missing = md_alloc(alloc, words * sizeof(uint64_t));
        memcpy(c.missing, v->missing, words * sizeof(uint64_t));
    }
    c.sets = sets_copy(v->sets, alloc);
    c.origin = sets_copy(v->origin, alloc);
    if (v->dtype == VIAL_DT_TUPLE && v->num_elems) {
        vial_val_t* el = md_alloc(alloc, v->num_elems * sizeof(vial_val_t));
        for (uint32_t i = 0; i < v->num_elems; ++i) el[i] = vial_val_copy(&v->elems[i], alloc);
        c.elems = el;
    }
    return c;
}
