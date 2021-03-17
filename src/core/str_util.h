#ifndef _MD_STR_UTIL_H
#define _MD_STR_UTIL_H

#include <stdint.h>
#include <stdbool.h>

typedef struct str_t str_t;
struct str_t {
    const char* ptr;
    uint64_t    len;
};

#ifdef __cplusplus
extern "C" {
#endif

// Only for ASCII character set
inline int  char_to_digit(int c) { return c - '0'; }
inline int  to_lower(int c)      { return ('A' <= c && c <= 'Z') ? c+32 : c; }
inline int  to_upper(int c)      { return ('a' <= c && c <= 'z') ? c-32 : c; }
inline bool is_digit(int c)      { return '0' <= c && c <= '9'; }
inline bool is_alpha(int c)      { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
inline bool is_whitespace(int c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }
inline bool is_symbol(int c) {
    return (32 < c && c < 48) || (57 < c && c < 65) || (90 < c && c < 97) || (122 < c && c < 127);
}
bool compare_str(const str_t str_a, const str_t str_b);
str_t substr(str_t str, uint64_t offset, uint64_t length);
bool skip_line(str_t* in_out_str);
bool peek_line(str_t* out_line, const str_t* in_str);
bool extract_line(str_t* out_line, str_t* in_out_str);
str_t trim_whitespace(str_t str);

// Make sure you ave trimmed all whitespace before using these!
// They are lean and mean
double parse_float(str_t str);
int64_t parse_int(str_t str);

struct md_allocator_i;

// Will allocate one extra character for zero termination
str_t alloc_str(uint64_t len, struct md_allocator_i* alloc);
void  free_str(str_t str, struct md_allocator_i* alloc);
str_t copy_str(const str_t str, struct md_allocator_i* alloc);
str_t load_textfile(const char* filename, uint32_t filename_len, struct md_allocator_i* alloc);


#ifdef __cplusplus
}
#endif

#endif