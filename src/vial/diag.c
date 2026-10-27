#include "vial_internal.h"

#include <core/md_common.h>

#include <string.h>
#include <stdio.h>

void vial_diags_init(vial_diags_t* d, str_t src, md_allocator_i* alloc) {
    memset(d, 0, sizeof(*d));
    d->src = src;
    d->alloc = alloc;
    md_array_push(d->line_beg, 0, alloc);
    for (size_t i = 0; i < src.len; ++i) {
        if (src.ptr[i] == '\n') md_array_push(d->line_beg, (uint32_t)(i + 1), alloc);
    }
}

uint32_t vial_line_of(const vial_diags_t* d, uint32_t offset) {
    // The last line start at or before offset
    size_t lo = 0, hi = md_array_size(d->line_beg);
    while (hi - lo > 1) {
        const size_t mid = (lo + hi) / 2;
        if (d->line_beg[mid] <= offset) lo = mid; else hi = mid;
    }
    return (uint32_t)lo + 1;
}

void vial_diag(vial_diags_t* d, const char code[6], md_vial_span_t span, const char* fmt, ...) {
    char buf[512];
    va_list args;
    va_start(args, fmt);
    int len = vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    if (len < 0) len = 0;
    if ((size_t)len >= sizeof(buf)) len = (int)sizeof(buf) - 1;

    md_vial_diagnostic_t diag = {0};
    memcpy(diag.code, code, 5);
    diag.severity = code[0] == 'W' ? MD_VIAL_SEVERITY_WARNING :
                    code[0] == 'R' ? MD_VIAL_SEVERITY_RUNTIME : MD_VIAL_SEVERITY_ERROR;
    diag.span = span;
    diag.message = str_copy_cstrn(buf, (size_t)len, d->alloc);
    md_array_push(d->items, diag, d->alloc);
}

size_t vial_diag_count(const vial_diags_t* d, const char code[6]) {
    size_t n = 0;
    for (size_t i = 0; i < md_array_size(d->items); ++i) {
        if (strncmp(d->items[i].code, code, 5) == 0) n += 1;
    }
    return n;
}

void vial_diag_set_fix(vial_diags_t* d, size_t index, str_t fix) {
    if (index < md_array_size(d->items)) d->items[index].fix = str_copy(fix, d->alloc);
}

void vial_diags_sort(vial_diags_t* d) {
    // Insertion sort, stable, by where each diagnostic starts: the lexer's and the parser's interleave by position
    md_vial_diagnostic_t* a = d->items;
    const size_t n = md_array_size(a);
    for (size_t i = 1; i < n; ++i) {
        const md_vial_diagnostic_t x = a[i];
        size_t j = i;
        while (j > 0 && a[j - 1].span.beg > x.span.beg) { a[j] = a[j - 1]; --j; }
        a[j] = x;
    }
}
