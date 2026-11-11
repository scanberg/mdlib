#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_common.h>
#include <core/md_str.h>

// A short name or identifier kept in place rather than behind a pointer: an atom or residue name, a chain id. Up to 8
// characters, the bytes after them zero. With all 8 used there is no terminator, so a label is read by its length
// (md_label_len, md_label_str) and never as a C string. Its length is where its first zero byte is, found in one 64-bit
// word, and two labels are equal when their bytes are: a label is a 64-bit key (md_label_key).
typedef struct md_label_t {
    char buf[8];
#ifdef __cplusplus
    inline operator md_str_t() const;
#endif
} md_label_t;

// The same for the longer names of a force field's types ('opls_135', 'CG2R61', 'SQ4p'): up to 16 characters
typedef struct md_label16_t {
    char buf[16];
#ifdef __cplusplus
    inline operator md_str_t() const;
#endif
} md_label16_t;

// macro concatenate trick to assert that the input is a valid compile time C-string
#define MD_LABEL(cstr) {cstr""}

// The index of the first zero byte of a little endian word, 8 when it has none. The borrow of a subtraction only runs
// upwards from a zero byte, so the lowest flagged byte is the first zero byte.
static inline size_t md_word_zero_byte(uint64_t v) {
    const uint64_t z = (v - 0x0101010101010101ull) & ~v & 0x8080808080808080ull;
    return z ? (size_t)(md_ctz64(z) >> 3) : 8;
}

static inline uint64_t md_label_key(md_label_t lbl) {
    uint64_t v;
    md_memcpy(&v, lbl.buf, sizeof(v));
    return v;
}

static inline size_t md_label_len(const md_label_t* lbl) {
    return md_word_zero_byte(md_label_key(*lbl));
}

static inline size_t md_label16_len(const md_label16_t* lbl) {
    uint64_t v[2];
    md_memcpy(v, lbl->buf, sizeof(v));
    const size_t n = md_word_zero_byte(v[0]);
    return n < 8 ? n : 8 + md_word_zero_byte(v[1]);
}

// A view of the label: valid while the label is where it is
static inline md_str_t md_label_str(const md_label_t* lbl) {
    md_str_t s = {lbl->buf, md_label_len(lbl)};
    return s;
}

static inline md_str_t md_label16_str(const md_label16_t* lbl) {
    md_str_t s = {lbl->buf, md_label16_len(lbl)};
    return s;
}

#define MD_LABEL_STR(lbl) md_label_str(&(lbl))

static inline bool md_label_eq(md_label_t a, md_label_t b) {
    return md_label_key(a) == md_label_key(b);
}

static inline bool md_label_empty(md_label_t lbl) {
    return md_label_key(lbl) == 0;
}

// The first 8 characters of str
static inline md_label_t md_label_make(md_str_t str) {
    md_label_t lbl = {0};
    if (str.ptr) {
        md_memcpy(lbl.buf, str.ptr, MD_MIN(str.len, sizeof(lbl.buf)));
    }
    return lbl;
}

// The first 16 characters of str
static inline md_label16_t md_label16_make(md_str_t str) {
    md_label16_t lbl = {0};
    if (str.ptr) {
        md_memcpy(lbl.buf, str.ptr, MD_MIN(str.len, sizeof(lbl.buf)));
    }
    return lbl;
}

#ifdef __cplusplus
inline md_label_t::operator md_str_t() const { return md_label_str(this); }
inline md_label16_t::operator md_str_t() const { return md_label16_str(this); }
inline bool operator==(const md_label_t& a, const md_label_t& b) { return md_label_eq(a, b); }
inline bool operator!=(const md_label_t& a, const md_label_t& b) { return !md_label_eq(a, b); }
#endif
