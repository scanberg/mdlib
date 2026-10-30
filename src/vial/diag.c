#include "vial_internal.h"

#include <core/md_common.h>

#include <string.h>
#include <stdio.h>

void vial_diags_init(vial_diags_t* d, str_t src, md_allocator_i* alloc) {
    memset(d, 0, sizeof(*d));
    d->src = src;
    d->alloc = alloc;
    const size_t words = src.len / 64 + 1;
    d->at = md_alloc(alloc, words * sizeof(uint64_t));
    memset(d->at, 0, words * sizeof(uint64_t));
    md_array_push(d->line_beg, 0, alloc);
    for (size_t i = 0; i < src.len; ++i) {
        if (src.ptr[i] == '\n') md_array_push(d->line_beg, (uint32_t)(i + 1), alloc);
    }
}

str_t vial_excerpt(char* buf, size_t cap, str_t text) {
    // At most VIAL_EXCERPT_MAX bytes of the first line, cut on a character; "…" marks a cut
    if (cap < VIAL_EXCERPT_MAX + 4) return (str_t){buf, 0};
    size_t n = 0;
    while (n < text.len && n < VIAL_EXCERPT_MAX && text.ptr[n] != '\n' && text.ptr[n] != '\r') ++n;
    const bool cut = n < text.len;
    while (n > 0 && n < text.len && ((unsigned char)text.ptr[n] & 0xC0) == 0x80) --n;
    if (n) memcpy(buf, text.ptr, n);     // an empty string's text may have no pointer at all
    if (cut) { memcpy(buf + n, "\xE2\x80\xA6", 3); n += 3; }
    buf[n] = '\0';
    return (str_t){buf, n};
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
    if ((size_t)len >= sizeof(buf)) {
        // Cut on a character, not inside one
        len = (int)sizeof(buf) - 1;
        while (len > 0 && ((unsigned char)buf[len] & 0xC0) == 0x80) --len;
    }

    md_vial_diagnostic_t diag = {0};
    memcpy(diag.code, code, 5);
    diag.severity = code[0] == 'W' ? MD_VIAL_SEVERITY_WARNING :
                    code[0] == 'R' ? MD_VIAL_SEVERITY_RUNTIME : MD_VIAL_SEVERITY_ERROR;
    diag.span = span;
    if (span.beg <= d->src.len) d->at[span.beg / 64] |= 1ull << (span.beg % 64);
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

bool vial_diag_at(const vial_diags_t* d, uint32_t offset) {
    return offset <= d->src.len && (d->at[offset / 64] >> (offset % 64)) & 1;
}

void vial_diags_sort(vial_diags_t* d) {
    // A bottom-up merge sort, stable, by where each diagnostic starts: the lexer's and the parser's interleave by
    // position, and the lexer reports an unclosed bracket at the end, where it learns of it
    const size_t n = md_array_size(d->items);
    if (n < 2) return;
    md_vial_diagnostic_t* a = d->items;
    md_vial_diagnostic_t* b = md_alloc(d->alloc, n * sizeof(md_vial_diagnostic_t));
    for (size_t w = 1; w < n; w *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * w) {
            const size_t mid = lo + w < n ? lo + w : n;
            const size_t hi = lo + 2 * w < n ? lo + 2 * w : n;
            size_t i = lo, j = mid, k = lo;
            while (i < mid && j < hi) b[k++] = a[j].span.beg < a[i].span.beg ? a[j++] : a[i++];
            while (i < mid) b[k++] = a[i++];
            while (j < hi) b[k++] = a[j++];
        }
        md_vial_diagnostic_t* t = a; a = b; b = t;
    }
    if (a != d->items) memcpy(d->items, a, n * sizeof(md_vial_diagnostic_t));
    md_free(d->alloc, a == d->items ? b : a, n * sizeof(md_vial_diagnostic_t));
}
