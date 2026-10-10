#pragma once

#include <core/md_str.h>
#include <core/md_os.h>

#include <string.h>

// Text parsing primitives: buffered line reading, whitespace and delimiter tokenization, and numbers.
//
// NUMBERS
//
// md_parse_f64 and md_parse_i64 are the two primitives. Each reads the longest PREFIX of its input
// that is a number and returns how many characters that was, 0 when the input does not start with
// one - so the caller decides whether whatever follows is an error. Neither skips leading
// whitespace, and neither reads a byte outside the view it is given.
//
//   float:    [+-]? ( digits ( '.' digits? )? | '.' digits ) ( [eE] [+-]? digits )?
//             [+-]? ( "inf" | "infinity" | "nan" ), any case
//   integer:  [+-]? digits
//
// An exponent marker without digits after it is not part of the number: "1e" reads as 1 with
// one character consumed, the way strtod reads it.
//
// The float value is CORRECTLY ROUNDED - the same double strtod gives - and does not depend on the
// process locale: a decimal comma locale, which a GTK file dialog can switch the process into, does
// not change what "1.5" means. A plain decimal of up to ~15 digits (every fixed width MD format) is
// one exact division; up to 19 significant digits with any exponent goes through Eisel-Lemire, as
// fast_float does it; only a longer mantissa whose 19 digit truncation straddles a rounding boundary
// is handed to strtod, in a form without a decimal point. This takes IEEE arithmetic in the default
// round to nearest mode, which is why md_parse.c is compiled without fast math.
//
// parse_float, parse_int, is_float and is_int are the long standing conveniences on top:
//   parse_float / parse_int  skip leading whitespace, read a prefix, and give 0 when there is none
//   is_float / is_int        whether the WHOLE string is a number, with nothing around it; is_float
//                            accepts the finite forms only, so "inf" and "nan" are not floats here
//                            (a script identifier named nan must stay an identifier)

#ifdef __cplusplus
extern "C" {
#endif

size_t md_parse_f64(double* out_value, str_t str);

bool md_parse_is_float(str_t str);

#ifdef __cplusplus
}
#endif

// Saturates to INT64_MIN / INT64_MAX, reading every digit either way. Inline, as it has no floating
// point to keep away from fast math.
static inline size_t md_parse_i64(int64_t* out_value, str_t str) {
    ASSERT(out_value);
    if (!str.ptr) return 0;

    const char* c   = str.ptr;
    const char* end = str.ptr + str.len;

    bool negative = false;
    if (c < end && (*c == '-' || *c == '+')) {
        negative = (*c == '-');
        ++c;
    }

    const char* digits = c;
    uint64_t value = 0;

    // 18 digits cannot overflow, so they go in without a check
    const char* unchecked_end = (size_t)(end - c) > 18 ? c + 18 : end;
    while (c < unchecked_end && is_digit(*c)) {
        value = value * 10 + (unsigned)(*c - '0');
        ++c;
    }
    if (c == digits) {
        return 0;
    }

    bool overflow = false;
    for (; c < end && is_digit(*c); ++c) {
        const unsigned d = (unsigned)(*c - '0');
        if (value > UINT64_MAX / 10 || (value == UINT64_MAX / 10 && d > UINT64_MAX % 10)) {
            overflow = true;
        } else if (!overflow) {
            value = value * 10 + d;
        }
    }

    // The magnitude of INT64_MIN is one more than INT64_MAX holds
    const uint64_t limit = negative ? (UINT64_C(1) << 63) : (uint64_t)INT64_MAX;
    if (overflow || value > limit) {
        *out_value = negative ? INT64_MIN : INT64_MAX;
    } else if (negative) {
        *out_value = value == (UINT64_C(1) << 63) ? INT64_MIN : -(int64_t)value;
    } else {
        *out_value = (int64_t)value;
    }
    return (size_t)(c - str.ptr);
}

static inline double parse_float(str_t str) {
    double value = 0.0;
    md_parse_f64(&value, str_trim_beg(str));
    return value;
}

static inline int64_t parse_int(str_t str) {
    int64_t value = 0;
    md_parse_i64(&value, str_trim_beg(str));
    return value;
}

static inline bool is_float(str_t str) {
    return md_parse_is_float(str);
}

static inline bool is_int(str_t str) {
    const char* c   = str.ptr;
    const char* end = str.ptr + str.len;
    if (c < end && (*c == '-' || *c == '+')) ++c;
    if (c >= end) return false;
    while (c < end && is_digit(*c)) ++c;
    return c == end;
}

// FIXED WIDTH FIELDS
//
// md_parse_fixed_f32 reads a field of a fixed width column format whose layout the format fixes - the
// "%8.3f" coordinates of gro and pdb: right aligned in the field, leading spaces, an optional '-', at
// least one digit, '.', and exactly 'decimals' digits up to the end of the field. A field laid out that
// way of up to 8 characters is checked and converted as one 64 bit word, with no scan and no branch on
// its characters; a wider one is checked character by character and converted by md_parse_f64.
//
// A field laid out any other way - another number of decimals, an exponent, a '+', trailing spaces,
// even when it is a number - gives false with *out untouched, so the caller decides what that means:
// the field parsed again by parse_float (pdb), or the line read another way (gro, whose columns are
// only a convention some writers do not keep to).
//
// The value is the one (float)parse_float(field) gives. For a field of up to 8 characters that holds in
// a fast math build as well, where the division below may become a multiplication by an inexact
// reciprocal: the double is then off by at most an ulp or two, while a decimal of at most 7 digits
// lies at least 2^-25 / 10^6 (relative) from every float rounding boundary - and is on one exactly
// only above 2^(24 - decimals), out of reach of 7 digits. So both round to the same float.

#ifdef __cplusplus
extern "C" {
#endif

bool md_parse_fixed_f32_wide(float* out, str_t field, size_t decimals);

#ifdef __cplusplus
}
#endif

// The value of eight ASCII digits loaded little endian, the first the most significant. From
// fast_float (parse_eight_digits_unrolled), MIT licensed.
static inline uint32_t md_parse_eight_digits(uint64_t val) {
    const uint64_t mask = 0x000000FF000000FFull;
    const uint64_t mul1 = 0x000F424000000064ull;   // 100 + (1000000 << 32)
    const uint64_t mul2 = 0x0000271000000001ull;   // 1 + (10000 << 32)
    val -= 0x3030303030303030ull;
    val = (val * 10) + (val >> 8);
    val = (((val & mask) * mul1) + (((val >> 16) & mask) * mul2)) >> 32;
    return (uint32_t)val;
}

static inline bool md_parse_fixed_f32(float* out, str_t field, size_t decimals) {
    ASSERT(out);
    const size_t width = field.len;
    // A digit, the dot and the decimals at least
    if (!field.ptr || decimals == 0 || width < decimals + 2) return false;
    if (width > 8) return md_parse_fixed_f32_wide(out, field, decimals);

    // The field right aligned in a word, spaces in front, byte i the i:th character (little endian, as
    // every target). A narrower field from two loads of 4 that overlap, none of them outside the field
    // and no store to forward from.
    uint64_t w;
    if (width == 8) {
        MEMCPY(&w, field.ptr, 8);
    } else if (width >= 4) {
        uint32_t lo, hi;
        MEMCPY(&lo, field.ptr, 4);
        MEMCPY(&hi, field.ptr + width - 4, 4);
        w = ((uint64_t)hi << 32) | ((uint64_t)lo << (8 * (8 - width))) | (0x2020202020202020ull >> (8 * width));
    } else {
        w = 0x2020202020202020ull;
        for (size_t i = 0; i < width; ++i) {
            w = (w & ~(0xFFull << (8 * (8 - width + i)))) | ((uint64_t)(uint8_t)field.ptr[i] << (8 * (8 - width + i)));
        }
    }

    const unsigned dot = 7 - (unsigned)decimals;    // 1..6
    if (((w >> (8 * dot)) & 0xFF) != '.') return false;

    // The dot taken out: what is in front of it moves up a byte, and a space comes in first
    const uint64_t front = (UINT64_C(1) << (8 * dot)) - 1;
    const uint64_t back  = ~((UINT64_C(1) << (8 * (dot + 1))) - 1);
    const uint64_t r = ((w & front) << 8) | (w & back) | 0x20;

    // ASCII only, so no byte below carries into the next
    if (r & 0x8080808080808080ull) return false;
    // 0x80 in every byte that is no digit: below '0', or above '9'
    const uint64_t other = (~(r + 0x5050505050505050ull) | (r + 0x4646464646464646ull)) & 0x8080808080808080ull;
    // Those are in front of the digits, which take in the last character in front of the dot (byte
    // 'dot' now), and they are spaces but for a '-' last
    const uint64_t lead = (other >> 7) * 0xFF;
    if ((lead & (lead + 1)) != 0 || (lead >> (8 * dot)) != 0) return false;
    const uint64_t last = lead ^ (lead >> 8);
    const uint64_t sign = (r ^ 0x2020202020202020ull) & lead;
    if ((sign & ~last) != 0 || (sign != 0 && sign != (0x0D0D0D0D0D0D0D0Dull & last))) return false;     // '-' ^ ' '

    const uint64_t digits = (r & ~lead) | (0x3030303030303030ull & lead);
    static const double pow10[7] = { 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6 };
    const double v = (double)md_parse_eight_digits(digits) / pow10[decimals];
    *out = (float)(sign ? -v : v);
    return true;
}

// LINES

// Reads as many COMPLETE lines from the file as fit in buf, and leaves the file positioned at the
// start of the first line it did not return. A single line longer than cap is returned in pieces.
static inline size_t md_parse_read_lines(md_file_t file, char* buf, size_t cap) {
    if (!md_file_valid(file) || !buf || cap < 1) return 0;
    size_t len = md_file_read(file, buf, cap);
    if (len == cap) {
        size_t loc;
        const str_t str = {buf, len};
        if (str_rfind_char(&loc, str, '\n')) {
            const int64_t offset = (int64_t)loc + 1 - (int64_t)len;
            // Set file pointer to the beginning of the next line
            md_file_seek(file, offset, MD_FILE_CUR);
            len = loc + 1;
        }
    }
    return len;
}

// Line by line reading from a file through a caller owned buffer, or from a string in memory, behind
// one interface. The lines are views into the buffer or the string: a file backed line is valid
// until the next call that refills the buffer.
typedef struct md_buffered_reader_t {
    str_t   str;
    char*   buf;
    size_t  cap;
    md_file_t  file;
} md_buffered_reader_t;

static inline md_buffered_reader_t md_buffered_reader_from_file(char* buf, size_t cap, md_file_t  file) {
    ASSERT(buf);
    ASSERT(md_file_valid(file));

    md_buffered_reader_t lr = {
        .str = {0, 0},
        .buf = buf,
        .cap = cap,
        .file = file,
    };

    return lr;
}

static inline md_buffered_reader_t md_buffered_reader_from_str(str_t str) {
    md_buffered_reader_t reader = {
        .str = str,
        .cap = str.len,
        .file = {0},
    };
    return reader;
}

static inline void md_buffered_reader_ensure_lines(md_buffered_reader_t* r) {
    ASSERT(r);
    if (md_file_valid(r->file) && !r->str.len) {
        ASSERT(r->buf);
        const size_t bytes_read = md_parse_read_lines(r->file, r->buf, r->cap);
        if (bytes_read > 0) {
            r->str.ptr = r->buf;
            r->str.len = bytes_read;
        }
    }
}

static inline bool md_buffered_reader_extract_line(str_t* line, md_buffered_reader_t* r) {
    ASSERT(r);
    ASSERT(line);
    md_buffered_reader_ensure_lines(r);
    return str_extract_line(line, &r->str);
}

static inline bool md_buffered_reader_peek_line(str_t* line, md_buffered_reader_t* r) {
    ASSERT(r);
    ASSERT(line);
    md_buffered_reader_ensure_lines(r);
    return str_peek_line(line, &r->str);
}

static inline bool md_buffered_reader_skip_line(md_buffered_reader_t* r) {
    ASSERT(r);
    md_buffered_reader_ensure_lines(r);
    return str_skip_line(&r->str);
}

// Back to the first line
static inline void md_buffered_reader_reset(md_buffered_reader_t* r) {
    ASSERT(r);
    if (md_file_valid(r->file)) {
        md_file_seek(r->file, 0, MD_FILE_BEG);
        r->str.ptr = NULL;
        r->str.len = 0;
    } else {
        ASSERT(r->str.ptr && "Cannot reset uninitialized reader");
        r->str.ptr += r->str.len - r->cap;
        r->str.len = r->cap;
    }
}

// Offset of the next unread byte from the start of the file or string
static inline int64_t md_buffered_reader_tellg(const md_buffered_reader_t* r) {
    if (md_file_valid(r->file)) {
        return md_file_tell(r->file) - r->str.len;
    } else {
        return r->cap - r->str.len;
    }
}

// TOKENS

// Extracts the next token delimited by whitespace, skipping any whitespace in front of it, and
// consumes the one delimiter character after it. False when only whitespace is left.
static inline bool extract_token(str_t* tok, str_t* str) {
    ASSERT(tok);
    ASSERT(str);
    if (str_empty(*str)) return false;

    const char* end = str->ptr + str->len;
    const char* c = str->ptr;
    while (c < end && is_whitespace(*c)) ++c;
    if (c >= end) return false;

    const char* tok_beg = c;
    while (c < end && !is_whitespace(*c)) ++c;

    tok->ptr = tok_beg;
    tok->len = c - tok_beg;

    str->ptr = c < end ? c + 1 : end;
    str->len = end - str->ptr;

    return true;
}

// Up to tok_cap whitespace delimited tokens
static inline size_t extract_tokens(str_t tok_arr[], size_t tok_cap, str_t* str) {
    ASSERT(tok_arr);
    ASSERT(str);

    size_t num_tokens = 0;
    while (num_tokens < tok_cap && extract_token(&tok_arr[num_tokens], str)) {
        num_tokens += 1;
    }
    return num_tokens;
}

// Extracts the field up to the next 'delim' (or the end) and consumes the delimiter. No whitespace is
// skipped and empty fields are returned as such: "a,,b" is "a", "" and "b". A delimiter at the very
// end does NOT give a trailing empty field - "a," is the one field "a" - as the input is used up.
static inline bool extract_token_delim(str_t* tok, str_t* str, char delim) {
    ASSERT(tok);
    ASSERT(str);
    if (!str->ptr || str->len == 0) return false;

    const char* beg = str->ptr;
    const char* end = str->ptr + str->len;
    const char* c = str->ptr;
    while (c != end && *c != delim) {
        ++c;
    }
    tok->ptr = beg;
    tok->len = c - beg;

    str->ptr = c != end ? c + 1 : end;
    str->len = end - str->ptr;

    return true;
}

// Up to tok_cap fields, see extract_token_delim
static inline size_t extract_tokens_delim(str_t tok_arr[], size_t tok_cap, str_t* str, char delim) {
    ASSERT(tok_arr);
    ASSERT(str);

    size_t num_tokens = 0;
    while (num_tokens < tok_cap && extract_token_delim(&tok_arr[num_tokens], str, delim)) {
        num_tokens += 1;
    }
    return num_tokens;
}
