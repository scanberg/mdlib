#ifndef _MD_STR_UTIL_H
#define _MD_STR_UTIL_H

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;

typedef struct str_t {
    const char* ptr;
    int64_t     len;
} str_t;

#ifdef __cplusplus
extern "C" {
#endif

#if 0
// This is just neat, used for inspiration
#define set(buffer, string) do {\
_Static_assert(sizeof(buffer)>=sizeof(string), "buffer to small");\
strcpy(buffer,string "");}\
while(0)
#endif

#define make_cstr(string) ((str_t) {.ptr = (string ""), .len = (sizeof(string)-1)})

// Only for ASCII character set
static inline int  char_to_digit(int c) { return c - '0'; }
static inline int  to_lower(int c)      { return ('A' <= c && c <= 'Z') ? c+32 : c; }
static inline int  to_upper(int c)      { return ('a' <= c && c <= 'z') ? c-32 : c; }
static inline bool is_digit(int c)      { return '0' <= c && c <= '9'; }
static inline bool is_alpha(int c)      { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
static inline bool is_whitespace(int c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }
static inline bool is_symbol(int c) {
    return (32 < c && c < 48) || (57 < c && c < 65) || (90 < c && c < 97) || (122 < c && c < 127);
}

static inline str_t trim_whitespace(str_t str) {
    const char* beg = str.ptr;
    const char* end = str.ptr + str.len;
    while (beg < end && is_whitespace(*beg)) ++beg;
    while (beg < end && (is_whitespace(end[-1]) || end[-1] == '\0')) --end;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

static inline bool compare_str(const str_t str_a, const str_t str_b) {
    if (!str_a.ptr || !str_b.ptr) return false;
    if (str_a.len != str_b.len) return false;
    for (uint64_t i = 0; i < str_a.len; ++i) {
        if (str_a.ptr[i] != str_b.ptr[i]) return false;
    }
    return true;
}

static inline str_t substr(str_t str, uint64_t offset, uint64_t length) {
    if (offset > str.len) {
        str_t res = {0};
        return res;   
    }
    if (offset + length > str.len) length = str.len - offset;
    str.ptr = str.ptr + offset;
    str.len = length;
    return str;
}

bool skip_line(str_t* in_out_str);
bool peek_line(str_t* out_line, const str_t* in_str);
bool extract_line(str_t* out_line, str_t* in_out_str);

int64_t find_char(str_t str, int c);
int64_t rfind_char(str_t str, int c);

// Make sure you ave trimmed all whitespace before using these!
// They are lean and mean
double parse_float(str_t str);
int64_t parse_int(str_t str);

// Will allocate one extra character for zero termination
//str_t make_cstr(const char* str);
str_t alloc_str(uint64_t len, struct md_allocator_i* alloc);
void  free_str(str_t str, struct md_allocator_i* alloc);
str_t copy_str(const str_t str, struct md_allocator_i* alloc);
str_t load_textfile(str_t path, struct md_allocator_i* alloc);

// c:/folder/file.ext -> ext
str_t extract_ext(str_t path);

// c:/folder/file.ext -> file.ext
str_t extract_file(str_t path);

// c:/folder/file.ext -> c:/folder/file.
str_t extract_path_without_ext(str_t path);

// c:/folder/file.ext -> c:/folder/
str_t extract_path_without_file(str_t path);

#ifdef __cplusplus
}
#endif

#endif