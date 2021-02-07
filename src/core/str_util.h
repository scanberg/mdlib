#ifndef _MD_STR_UTIL_H
#define _MD_STR_UTIL_H

#include <stdint.h>
#include <stdbool.h>

typedef struct str_t str_t;
struct str_t {
    const char* str;
    uint64_t    len;
};

#ifdef __cplusplus
extern "C" {
#endif

inline int  char_to_digit(char c) { return c - '0'; }
inline char to_lower(char c)      { return ('A' <= c && c <= 'Z') ? c+32 : c; }
inline char to_upper(char c)      { return ('a' <= c && c <= 'z') ? c-32 : c; }
inline bool is_digit(char c)      { return '0' <= c && c <= '9'; }
inline bool is_alpha(char c)      { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
inline bool is_whitespace(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

str_t substr(str_t str, int offset, int length);
bool skip_line(str_t* in_out_str);
bool peek_line(str_t* out_line, const str_t* in_str);
bool extract_line(str_t* out_line, str_t* in_out_str);
str_t trim_whitespace(str_t str);

double parse_float(str_t str);
int64_t parse_int(str_t str);

struct md_allocator_i;

str_t load_textfile(const char* filename, uint32_t filename_len, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

#endif