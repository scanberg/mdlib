#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>

#include <core/md_common.h>  // ASSERT, MIN, MAX

#ifdef __cplusplus
#define DEF_VAL(x) = x
#else
#define DEF_VAL(x)
#endif

struct md_allocator_i;

// Non owning string_view type
typedef struct str_t {
    const char* ptr;
    size_t      len;

#ifdef __cplusplus
    constexpr char operator[](size_t idx) noexcept       { return ptr[idx]; }
    constexpr char operator[](size_t idx) const noexcept { return ptr[idx]; }

    constexpr const char* beg() noexcept { return ptr; }
    constexpr const char* end() noexcept { return ptr + len; }
    constexpr operator bool() noexcept { return ptr && len > 0; }
#endif
} str_t;

// Macros to bake a string literal into a str_t.
// The empty string concatenation (cstr"") asserts that the input is a compile time string literal.
//
// STR_INIT: Braced initializer. Use when initializing a str_t (static/const variables, array elements, struct fields),
//           since compound literals are not constant expressions in C (MSVC rejects them in static initializers).
//           e.g. static const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
// STR_LIT:  Expression of type str_t. Use as function arguments, in assignments, return values, member access, etc.
//           e.g. md_attributes_find(&attr, STR_LIT("atom/position"));
#define STR_INIT(cstr) {(cstr""), (sizeof(cstr)-1)}
#ifdef __cplusplus
#define STR_LIT(cstr) (str_t STR_INIT(cstr))
#else
#define STR_LIT(cstr) ((str_t)STR_INIT(cstr))
#endif

#define STR_FMT "%.*s"
#define STR_ARG(str) (int)str.len, str.ptr

#ifdef __cplusplus
extern "C" {
#endif

// Only for ASCII character set
static inline int  char_to_digit(int c) { return c - '0'; }
static inline int  to_lower(int c)      { return ('A' <= c && c <= 'Z') ? c+32 : c; }
static inline int  to_upper(int c)      { return ('a' <= c && c <= 'z') ? c-32 : c; }
static inline bool is_digit(int c)      { return '0' <= c && c <= '9'; }
static inline bool is_alpha(int c)      { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
static inline bool is_whitespace(int c) { return c == ' ' || c == '\n' || c == '\r' || c == '\t'; }
static inline bool is_symbol(int c) {
    return (32 < c && c < 48) || (57 < c && c < 65) || (90 < c && c < 97) || (122 < c && c < 127);
}

static inline bool str_empty(str_t str) { return str.ptr == 0 || str.len == 0; }

static inline const char* str_beg(str_t str) { return str.ptr; }
static inline const char* str_end(str_t str) { return str.ptr + str.len; }

static inline size_t      str_len(str_t str) { return str.len; }
static inline const char* str_ptr(str_t str) { return str.ptr; }

static inline void str_swap(str_t* a, str_t* b) {
    str_t tmp = *a;
    *a = *b;
    *b = tmp;
}

static inline str_t str_trim_beg(str_t str) {
    const char* beg = str_beg(str);
    const char* end = str_end(str);
    while (beg < end && is_whitespace(*beg)) ++beg;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

static inline str_t str_trim_end(str_t str) {
    const char* beg = str_beg(str);
    const char* end = str_end(str);
    while (beg < end && (is_whitespace(end[-1]) || end[-1] == '\0')) --end;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

static inline str_t str_trim(str_t str) {
    const char* beg = str_beg(str);
    const char* end = str_end(str);
    while (beg < end && is_whitespace(*beg)) ++beg;
    while (beg < end && (is_whitespace(end[-1]) || end[-1] == '\0')) --end;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

// A str_t is a view: its value is the len characters at ptr, and where it points does not matter when len is 0.
// {NULL,0}, "" and an empty slice of another string are the same, empty, string - as with std::string_view.
// The comparisons below never look at ptr for an empty str, and take a NULL cstr as "" (as str_from_cstr does).
//
// They are inlined: lookups and parsers call them per element and per line.

static inline bool str_eq(str_t str_a, str_t str_b) {
    if (str_a.len != str_b.len) return false;
    if (str_a.len == 0) return true;
    if (str_a.ptr[0] != str_b.ptr[0]) return false;
    return MEMCMP(str_a.ptr, str_b.ptr, str_a.len) == 0;
}

static inline bool str_eq_ignore_case(str_t str_a, str_t str_b) {
    if (str_a.len != str_b.len) return false;
    for (size_t i = 0; i < str_a.len; ++i) {
        if (to_lower(str_a.ptr[i]) != to_lower(str_b.ptr[i])) return false;
    }
    return true;
}

// Equal in their first n characters: the prefixes of length MIN(len, n) are the same string.
// With n = prefix.len this is str_begins_with(str, prefix).
static inline bool str_eq_n(str_t str_a, str_t str_b, size_t n) {
    const size_t len = MIN(str_a.len, n);
    if (len != MIN(str_b.len, n)) return false;
    return len == 0 || MEMCMP(str_a.ptr, str_b.ptr, len) == 0;
}

static inline bool str_eq_n_ignore_case(str_t str_a, str_t str_b, size_t n) {
    const size_t len = MIN(str_a.len, n);
    if (len != MIN(str_b.len, n)) return false;
    for (size_t i = 0; i < len; ++i) {
        if (to_lower(str_a.ptr[i]) != to_lower(str_b.ptr[i])) return false;
    }
    return true;
}

static inline bool str_eq_cstr(str_t str, const char* cstr) {
    if (!cstr) return str.len == 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (cstr[i] == '\0' || cstr[i] != str.ptr[i]) return false;
    }
    return cstr[str.len] == '\0';
}

static inline bool str_eq_cstr_ignore_case(str_t str, const char* cstr) {
    if (!cstr) return str.len == 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (cstr[i] == '\0' || to_lower(str.ptr[i]) != to_lower(cstr[i])) return false;
    }
    return cstr[str.len] == '\0';
}

// str_eq_n with cstr as a str. With n = strlen(cstr) it tells whether str begins with cstr: a str shorter than
// cstr does not ("END" is not "ENDMDL").
static inline bool str_eq_cstr_n(str_t str, const char* cstr, size_t n) {
    const size_t len = MIN(str.len, n);
    if (!cstr) return len == 0;
    for (size_t i = 0; i < len; ++i) {
        if (cstr[i] == '\0' || cstr[i] != str.ptr[i]) return false;
    }
    return len == n || cstr[len] == '\0';
}

static inline bool str_eq_cstr_n_ignore_case(str_t str, const char* cstr, size_t n) {
    const size_t len = MIN(str.len, n);
    if (!cstr) return len == 0;
    for (size_t i = 0; i < len; ++i) {
        if (cstr[i] == '\0' || to_lower(str.ptr[i]) != to_lower(cstr[i])) return false;
    }
    return len == n || cstr[len] == '\0';
}

// Lexicographical order of the characters as unsigned bytes (as memcmp, strcmp and std::string_view::compare),
// with a prefix ordered before what it is a prefix of. Returns -1, 0 or 1.
// A byte loop rather than memcmp: what is ordered is mostly short (names, labels), where the call costs more than it saves.
static inline int str_cmp_lex(str_t a, str_t b) {
    const size_t len = MIN(a.len, b.len);
    for (size_t i = 0; i < len; ++i) {
        const unsigned char ca = (unsigned char)a.ptr[i];
        const unsigned char cb = (unsigned char)b.ptr[i];
        if (ca != cb) return ca < cb ? -1 : 1;
    }
    return (a.len > b.len) - (a.len < b.len);
}

static inline bool str_begins_with(str_t str, str_t prefix) {
    return str.len >= prefix.len && (prefix.len == 0 || MEMCMP(str.ptr, prefix.ptr, prefix.len) == 0);
}

static inline bool str_ends_with(str_t str, str_t suffix) {
    return str.len >= suffix.len && (suffix.len == 0 || MEMCMP(str.ptr + str.len - suffix.len, suffix.ptr, suffix.len) == 0);
}

size_t str_count_equal_chars(str_t a, str_t b);
size_t str_count_occur_char(str_t str, char character);

static inline str_t str_substr(str_t str, size_t offset, size_t length DEF_VAL(SIZE_MAX)) {
    str_t res = {0,0};
    if (offset > str.len) {
        return res;
    }
    const size_t max_len = str.len - offset;
    res.ptr = str.ptr + offset;
    res.len = length = MIN(length, max_len);
    return res;
}

static inline str_t str_from_cstr(const char* cstr) {
    str_t str = {0,0};
    if (cstr) {
        str.ptr = cstr;
        str.len = strlen(cstr);
    }
    return str;
}

static inline str_t str_from_cstrn(const char* cstr, size_t n) {
    str_t str = {0,0};
    if (cstr) {
        str.ptr = cstr;
        str.len = strnlen(cstr, n);
    }
    return str;
}

// Joins two strings by taking the beginning of the first and the end of the last
// @WARNING:
// This does no allocation of any sort, it simply creates a new string view from the two inputs which correspond to the total span
// Use at your own risk!
static inline str_t str_join(str_t first, str_t last) {
    ASSERT(first.ptr < last.ptr + last.len);
    str_t str = { first.ptr, (size_t)(last.ptr + last.len - first.ptr) };
    return str;
}

// Skip, peek and extract agree on what a line is: the text up to a '\n' or the end, without the '\r' of a CRLF,
// and a last line counts whether or not a newline ends it. Inlined, as parsers call them for every line.
static inline bool str_extract_line(str_t* out_line, str_t* in_out_str) {
    ASSERT(out_line);
    ASSERT(in_out_str);
    if (in_out_str->len == 0) return false;

    const char* beg = in_out_str->ptr;
    const char* end = beg + in_out_str->len;
    const char* nl  = (const char*)memchr(beg, '\n', in_out_str->len);
    const char* line_end = nl ? nl     : end;
    const char* next     = nl ? nl + 1 : end;
    if (nl && line_end > beg && line_end[-1] == '\r') --line_end;

    out_line->ptr = beg;
    out_line->len = (size_t)(line_end - beg);
    in_out_str->ptr = next;
    in_out_str->len = (size_t)(end - next);
    return true;
}

static inline bool str_peek_line(str_t* out_line, const str_t* in_str) {
    ASSERT(in_str);
    str_t rest = *in_str;
    return str_extract_line(out_line, &rest);
}

static inline bool str_skip_line(str_t* in_out_str) {
    ASSERT(in_out_str);
    if (in_out_str->len == 0) return false;
    const char* nl   = (const char*)memchr(in_out_str->ptr, '\n', in_out_str->len);
    const char* next = nl ? nl + 1 : in_out_str->ptr + in_out_str->len;
    in_out_str->len -= (size_t)(next - in_out_str->ptr);
    in_out_str->ptr  = next;
    return true;
}

bool str_extract_i32  (int32_t* val, str_t* in_out_str);
bool str_extract_i64  (int64_t* val, str_t* in_out_str);

bool str_extract_f32  (float*  val,  str_t* in_out_str);
bool str_extract_f64  (double* val,  str_t* in_out_str);

// Return true if the item is found, and write its offset to loc (when loc is not NULL).
// A character is matched as a byte: c is converted to char (as memchr converts it to unsigned char).
// str_find_str finds nothing when the needle is empty or longer than the haystack.
static inline bool str_find_char(size_t* loc, str_t str, int c) {
    if (str.len == 0) return false;
    const char* p = (const char*)memchr(str.ptr, c, str.len);
    if (!p) return false;
    if (loc) *loc = (size_t)(p - str.ptr);
    return true;
}

static inline bool str_rfind_char(size_t* loc, str_t str, int c) {
    const char ch = (char)c;
    for (size_t i = str.len; i > 0; --i) {
        if (str.ptr[i - 1] == ch) {
            if (loc) *loc = i - 1;
            return true;
        }
    }
    return false;
}

bool str_find_str(size_t* loc, str_t haystack, str_t needle);

// Will allocate one extra character for zero termination
str_t str_alloc(size_t len, struct md_allocator_i* alloc);
void  str_free(str_t str, struct md_allocator_i* alloc);
str_t str_copy(str_t str, struct md_allocator_i* alloc);
str_t str_copy_cstr(const char* cstr, struct md_allocator_i* alloc);
str_t str_copy_cstrn(const char* cstr, size_t len, struct md_allocator_i* alloc);

str_t str_printf(struct md_allocator_i* alloc, const char* format, ...);

#if defined(__GNUC__) || defined(__clang__)
#define MD_STR_PRINTF(fmt_idx, arg_idx) __attribute__((format(printf, fmt_idx, arg_idx)))
#else
#define MD_STR_PRINTF(fmt_idx, arg_idx)
#endif

// snprintf and vsnprintf as C99 has them, in the "C" locale whatever the process or the calling thread
// has set: the decimal point is always '.'. The process locale is not ours to rely on. A GTK file dialog
// switches it to the user's, and so may a shell extension a Windows file dialog loads, and under a decimal
// comma locale the C library writes 1.5 as "1,5" - which md_parse, as it does not follow the locale either,
// reads back as 1 and 5. Everything mdlib formats goes through these (str_printf, the writers,
// md_file_printf and the log too), so what it writes is the same text in every locale.
//
// Implemented in md_os.c: each call is handed a "C" locale of its own through the platform's API, rather
// than setlocale, which is process wide and would change what every other thread prints meanwhile.
int md_snprintf (char* buf, size_t cap, const char* format, ...) MD_STR_PRINTF(3, 4);
int md_vsnprintf(char* buf, size_t cap, const char* format, va_list args);

// This should probably be removed
str_t load_textfile(str_t path, struct md_allocator_i* alloc);

// c:/folder/file.ext -> ext
bool extract_ext(str_t* ext, str_t path);

// c:/folder/file.ext -> file.ext
bool extract_file(str_t* file, str_t path);

// c:/folder/file.ext -> c:/folder/file
bool extract_file_path_without_ext(str_t* file_path, str_t path);

// c:/folder/file.ext -> c:/folder/
bool extract_folder_path(str_t* folder_path, str_t path);

// Mutating string operations
void replace_char(char* str, size_t len, char c, char replacement);
void convert_to_lower(char* str, size_t len);
void convert_to_upper(char* str, size_t len);

// Copies the contents of a str_t into a char buffer and ensures zero termination (also for an empty str_t, which
// leaves an empty string). Truncates to fit in cap. Returns the number of characters written (excluding the zero
// termination character)
static inline size_t str_copy_to_char_buf(char* buf, size_t cap, str_t str) {
    ASSERT(buf);
    if (cap == 0) return 0;
    const size_t len = MIN(str.len, cap - 1);
    if (len) MEMCPY(buf, str.ptr, len);
    buf[len] = '\0';
    return len;
}

// Returns the Levenshtein edit distance between two strings
// https://en.wikipedia.org/wiki/Levenshtein_distance
int str_edit_distance(str_t a, str_t b);

#ifdef __cplusplus
}
#endif
