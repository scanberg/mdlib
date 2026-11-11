#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#include <core/md_common.h>  // MD_ASSERT, MD_MIN, MD_MAX

#ifdef __cplusplus
#define MD_DEF_VAL(x) = x
#else
#define MD_DEF_VAL(x)
#endif

// Non owning string_view type
typedef struct md_str_t {
    const char* ptr;
    size_t      len;

#ifdef __cplusplus
    constexpr char operator[](size_t idx) noexcept       { return ptr[idx]; }
    constexpr char operator[](size_t idx) const noexcept { return ptr[idx]; }

    constexpr const char* beg() noexcept { return ptr; }
    constexpr const char* end() noexcept { return ptr + len; }
    constexpr operator bool() noexcept { return ptr && len > 0; }
#endif
} md_str_t;

// Macros to bake a string literal into an md_str_t.
// The empty string concatenation (cstr"") asserts that the input is a compile time string literal.
//
// MD_STR_INIT: Braced initializer. Use when initializing an md_str_t (static/const variables, array elements, struct
//              fields), since compound literals are not constant expressions in C (MSVC rejects them in static initializers).
//              e.g. static const md_str_t paths[] = { MD_STR_INIT("atom/position"), MD_STR_INIT("unitcell") };
// MD_STR_LIT:  Expression of type md_str_t. Use as function arguments, in assignments, return values, member access, etc.
//              e.g. md_attributes_find(&attr, MD_STR_LIT("atom/position"));
#define MD_STR_INIT(cstr) {(cstr""), (sizeof(cstr)-1)}
#ifdef __cplusplus
#define MD_STR_LIT(cstr) (md_str_t MD_STR_INIT(cstr))
#else
#define MD_STR_LIT(cstr) ((md_str_t)MD_STR_INIT(cstr))
#endif

#define MD_STR_FMT "%.*s"
#define MD_STR_ARG(str) (int)str.len, str.ptr

#ifdef __cplusplus
extern "C" {
#endif

// Only for ASCII character set
static inline int  md_ascii_char_to_digit(int c) { return c - '0'; }
static inline int  md_ascii_to_lower(int c)      { return ('A' <= c && c <= 'Z') ? c+32 : c; }
static inline int  md_ascii_to_upper(int c)      { return ('a' <= c && c <= 'z') ? c-32 : c; }
static inline bool md_ascii_is_digit(int c)      { return '0' <= c && c <= '9'; }
static inline bool md_ascii_is_alpha(int c)      { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
static inline bool md_ascii_is_whitespace(int c) { return c == ' ' || c == '\n' || c == '\r' || c == '\t'; }
static inline bool md_ascii_is_symbol(int c) {
    return (32 < c && c < 48) || (57 < c && c < 65) || (90 < c && c < 97) || (122 < c && c < 127);
}

static inline bool md_str_empty(md_str_t str) { return str.ptr == 0 || str.len == 0; }

static inline const char* md_str_beg(md_str_t str) { return str.ptr; }
static inline const char* md_str_end(md_str_t str) { return str.ptr + str.len; }

static inline size_t      md_str_len(md_str_t str) { return str.len; }
static inline const char* md_str_ptr(md_str_t str) { return str.ptr; }

static inline void md_str_swap(md_str_t* a, md_str_t* b) {
    md_str_t tmp = *a;
    *a = *b;
    *b = tmp;
}

static inline md_str_t md_str_trim_beg(md_str_t str) {
    const char* beg = md_str_beg(str);
    const char* end = md_str_end(str);
    while (beg < end && md_ascii_is_whitespace(*beg)) ++beg;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

static inline md_str_t md_str_trim_end(md_str_t str) {
    const char* beg = md_str_beg(str);
    const char* end = md_str_end(str);
    while (beg < end && (md_ascii_is_whitespace(end[-1]) || end[-1] == '\0')) --end;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

static inline md_str_t md_str_trim(md_str_t str) {
    const char* beg = md_str_beg(str);
    const char* end = md_str_end(str);
    while (beg < end && md_ascii_is_whitespace(*beg)) ++beg;
    while (beg < end && (md_ascii_is_whitespace(end[-1]) || end[-1] == '\0')) --end;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

// A md_str_t is a view: its value is the len characters at ptr, and where it points does not matter when len is 0.
// {NULL,0}, "" and an empty slice of another string are the same, empty, string - as with std::string_view.
// The comparisons below never look at ptr for an empty str, and take a NULL cstr as "" (as md_str_from_cstr does).
//
// They are inlined: lookups and parsers call them per element and per line.

static inline bool md_str_eq(md_str_t str_a, md_str_t str_b) {
    if (str_a.len != str_b.len) return false;
    if (str_a.len == 0) return true;
    if (str_a.ptr[0] != str_b.ptr[0]) return false;
    return md_memcmp(str_a.ptr, str_b.ptr, str_a.len) == 0;
}

static inline bool md_str_eq_ignore_case(md_str_t str_a, md_str_t str_b) {
    if (str_a.len != str_b.len) return false;
    for (size_t i = 0; i < str_a.len; ++i) {
        if (md_ascii_to_lower(str_a.ptr[i]) != md_ascii_to_lower(str_b.ptr[i])) return false;
    }
    return true;
}

static inline bool md_str_eq_cstr(md_str_t str, const char* cstr) {
    if (!cstr) return str.len == 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (cstr[i] == '\0' || cstr[i] != str.ptr[i]) return false;
    }
    return cstr[str.len] == '\0';
}

static inline bool md_str_eq_cstr_ignore_case(md_str_t str, const char* cstr) {
    if (!cstr) return str.len == 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (cstr[i] == '\0' || md_ascii_to_lower(str.ptr[i]) != md_ascii_to_lower(cstr[i])) return false;
    }
    return cstr[str.len] == '\0';
}

// Lexicographical order of the characters as unsigned bytes (as memcmp, strcmp and std::string_view::compare),
// with a prefix ordered before what it is a prefix of. Returns -1, 0 or 1.
// A byte loop rather than memcmp: what is ordered is mostly short (names, labels), where the call costs more than it saves.
static inline int md_str_cmp_lex(md_str_t a, md_str_t b) {
    const size_t len = MD_MIN(a.len, b.len);
    for (size_t i = 0; i < len; ++i) {
        const unsigned char ca = (unsigned char)a.ptr[i];
        const unsigned char cb = (unsigned char)b.ptr[i];
        if (ca != cb) return ca < cb ? -1 : 1;
    }
    return (a.len > b.len) - (a.len < b.len);
}

// Prefix and suffix tests. A literal one is measured by MD_STR_LIT: md_str_begins_with(line, MD_STR_LIT("ATOM")).
static inline bool md_str_begins_with(md_str_t str, md_str_t prefix) {
    return str.len >= prefix.len && (prefix.len == 0 || md_memcmp(str.ptr, prefix.ptr, prefix.len) == 0);
}

static inline bool md_str_begins_with_ignore_case(md_str_t str, md_str_t prefix) {
    if (str.len < prefix.len) return false;
    for (size_t i = 0; i < prefix.len; ++i) {
        if (md_ascii_to_lower(str.ptr[i]) != md_ascii_to_lower(prefix.ptr[i])) return false;
    }
    return true;
}

static inline bool md_str_ends_with(md_str_t str, md_str_t suffix) {
    return str.len >= suffix.len && (suffix.len == 0 || md_memcmp(str.ptr + str.len - suffix.len, suffix.ptr, suffix.len) == 0);
}

// The number of times c occurs in str
static inline size_t md_str_count_occur_char(md_str_t str, char c) {
    size_t count = 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (str.ptr[i] == c) count += 1;
    }
    return count;
}

static inline md_str_t md_str_substr(md_str_t str, size_t offset, size_t length MD_DEF_VAL(SIZE_MAX)) {
    md_str_t res = {0,0};
    if (offset > str.len) {
        return res;
    }
    const size_t max_len = str.len - offset;
    res.ptr = str.ptr + offset;
    res.len = length = MD_MIN(length, max_len);
    return res;
}

static inline md_str_t md_str_from_cstr(const char* cstr) {
    md_str_t str = {0,0};
    if (cstr) {
        str.ptr = cstr;
        str.len = strlen(cstr);
    }
    return str;
}

static inline md_str_t md_str_from_cstrn(const char* cstr, size_t n) {
    md_str_t str = {0,0};
    if (cstr) {
        str.ptr = cstr;
        str.len = strnlen(cstr, n);
    }
    return str;
}

// Joins two strings by taking the beginning of the first and the end of the last
// A position that is not there: what the finds give when they find nothing. It is SIZE_MAX, which md_str_substr
// reads as "to the end", so md_str_substr(s, 0, md_str_find_char(s, ',')) is s up to its first ',' or all of s.
// And MD_STR_NPOS + 1 is 0: md_str_substr(s, md_str_rfind_char(s, '/') + 1, SIZE_MAX) is what follows the last
// '/', or all of s. To split in two, use md_str_cut.
#define MD_STR_NPOS SIZE_MAX

// ### CHARACTER RUNS ###
// On the len characters at chars: no terminator is needed, and none is written. (NULL, 0) is fine.
// Characters are matched as bytes.

// The position of the first (find) or the last (rfind) c, or MD_STR_NPOS
static inline size_t md_chars_find(const char* chars, size_t len, char c) {
    if (len == 0) return MD_STR_NPOS;
    const char* p = (const char*)memchr(chars, (unsigned char)c, len);
    return p ? (size_t)(p - chars) : MD_STR_NPOS;
}

static inline size_t md_chars_rfind(const char* chars, size_t len, char c) {
    // A loop: memrchr is glibc's only
    for (size_t i = len; i > 0; --i) {
        if (chars[i - 1] == c) return i - 1;
    }
    return MD_STR_NPOS;
}

// In place
static inline void md_chars_replace(char* chars, size_t len, char target, char replacement) {
    MD_ASSERT(chars || len == 0);
    for (size_t i = 0; i < len; ++i) {
        if (chars[i] == target) chars[i] = replacement;
    }
}

// md_ascii_to_lower / md_ascii_to_upper on each character
static inline void md_chars_to_lower(char* chars, size_t len) {
    MD_ASSERT(chars || len == 0);
    for (size_t i = 0; i < len; ++i) {
        chars[i] = (char)md_ascii_to_lower(chars[i]);
    }
}

static inline void md_chars_to_upper(char* chars, size_t len) {
    MD_ASSERT(chars || len == 0);
    for (size_t i = 0; i < len; ++i) {
        chars[i] = (char)md_ascii_to_upper(chars[i]);
    }
}

// ### FIND ###
// Positions are indices into str, MD_STR_NPOS when there is none.

static inline size_t md_str_find_char (md_str_t str, char c) { return md_chars_find (str.ptr, str.len, c); }
static inline size_t md_str_rfind_char(md_str_t str, char c) { return md_chars_rfind(str.ptr, str.len, c); }

// The first needle in haystack. An empty needle, or one longer than haystack, is not found.
static inline size_t md_str_find_str(md_str_t haystack, md_str_t needle) {
    if (needle.len == 0 || needle.len > haystack.len) return MD_STR_NPOS;

    // Every position where the needle still fits is a candidate. (Restarting the match after the character that broke
    // a partial match is not enough: "aab" contains "ab" starting inside the partial match "aa".)
    const char* beg  = haystack.ptr;
    const char* last = haystack.ptr + (haystack.len - needle.len);
    for (const char* c = beg; c <= last; ++c) {
        c = (const char*)memchr(c, (unsigned char)needle.ptr[0], (size_t)(last - c) + 1);
        if (!c) break;
        if (md_memcmp(c, needle.ptr, needle.len) == 0) {
            return (size_t)(c - beg);
        }
    }
    return MD_STR_NPOS;
}

static inline bool md_str_contains_char(md_str_t str, char c)               { return md_str_find_char(str, c) != MD_STR_NPOS; }
static inline bool md_str_contains_str (md_str_t haystack, md_str_t needle) { return md_str_find_str(haystack, needle) != MD_STR_NPOS; }

// Splits str at loc, a position from one of the finds, and drops the character there: before is what precedes it,
// after what follows it. md_str_cut(&key, &val, s, md_str_find_char(s, '=')) splits at the first '=', and with
// md_str_rfind_char at the last. MD_STR_NPOS (or any loc past the end) cuts nothing: before is all of str, after
// is empty, and the result is false. Either out pointer may be NULL, and either may point at str itself - but as a
// miss empties after, a loop that needs what is left once nothing more is found cuts into a variable of its own.
static inline bool md_str_cut(md_str_t* before, md_str_t* after, md_str_t str, size_t loc) {
    md_str_t b = str;
    md_str_t a = {0, 0};
    const bool found = loc < str.len;
    if (found) {
        b.len = loc;
        a.ptr = str.ptr + loc + 1;
        a.len = str.len - loc - 1;
    }
    if (before) *before = b;
    if (after)  *after  = a;
    return found;
}

// Copies the contents of a md_str_t into a char buffer and ensures zero termination (also for an empty md_str_t, which
// leaves an empty string). Truncates to fit in cap. Returns the number of characters written (excluding the zero
// termination character)
static inline size_t md_str_copy_to_char_buf(char* buf, size_t cap, md_str_t str) {
    MD_ASSERT(buf);
    if (cap == 0) return 0;
    const size_t len = MD_MIN(str.len, cap - 1);
    if (len) md_memcpy(buf, str.ptr, len);
    buf[len] = '\0';
    return len;
}

#ifdef __cplusplus
}
#endif
