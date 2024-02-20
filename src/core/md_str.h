#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

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

// Macro to bake a string literal into a str_t
#ifdef __cplusplus
#define STR_LIT(cstr) {(cstr""), (sizeof(cstr)-1)}
#else
#define STR_LIT(cstr) (str_t){(cstr""), (sizeof(cstr)-1)}
#endif

#define STR_FMT "%.*s"
#define STR_ARG(str) (int)str.len, str.ptr

#ifdef __cplusplus
extern "C" {
#endif

str_t str_from_cstr(const char* cstr);

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

void str_swap(str_t a, str_t b);

str_t str_trim(str_t str);
str_t str_trim_beg(str_t str);
str_t str_trim_end(str_t str);

bool str_eq(str_t str_a, str_t str_b);
bool str_eq_ignore_case(const str_t str_a, const str_t str_b);
bool str_eq_n(str_t str_a, str_t str_b, size_t n);
bool str_eq_n_ignore_case(const str_t str_a, const str_t str_b, size_t n);

// Lexicographical comparison
int str_cmp_lex(str_t a, str_t b);

bool str_eq_cstr(str_t str, const char* cstr);
bool str_eq_cstr_n(str_t str, const char* cstr, size_t n);
bool str_eq_cstr_ignore_case(str_t str, const char* cstr);
bool str_eq_cstr_n_ignore_case(str_t str, const char* cstr, size_t n);

size_t str_count_equal_chars(str_t a, str_t b);
size_t str_count_occur_char(str_t str, char character);

str_t str_substr(str_t str, size_t offset, size_t length DEF_VAL(SIZE_MAX));

// Joins two strings by taking the beginning of the first and the end of the last
// @WARNING:
// This does no allocation of any sort, it simply creates a new string view from the two inputs which correspond to the total span
// Use at your own risk!
str_t str_join(str_t first, str_t last);

bool str_skip_line   (str_t* in_out_str);
bool str_peek_line   (str_t* out_line, const str_t* in_str);
bool str_extract_line(str_t* out_line, str_t* in_out_str);

bool str_extract_i32  (int32_t* val, str_t* in_out_str);
bool str_extract_i64  (int64_t* val, str_t* in_out_str);

bool str_extract_f32  (float*  val,  str_t* in_out_str);
bool str_extract_f64  (double* val,  str_t* in_out_str);

// Returns the offset where the item was found or -1 if nothing was found.
bool str_find_char (size_t* loc, str_t str, int c);
bool str_rfind_char(size_t* loc, str_t str, int c);
bool str_find_str  (size_t* loc, str_t haystack, str_t needle);

bool str_begins_with(str_t str, str_t prefix);
bool str_ends_with(str_t str, str_t suffix);

// Will allocate one extra character for zero termination
str_t str_alloc(size_t len, struct md_allocator_i* alloc);
void  str_free(str_t str, struct md_allocator_i* alloc);
str_t str_copy(str_t str, struct md_allocator_i* alloc);
str_t str_copy_cstr(const char* cstr, struct md_allocator_i* alloc);
str_t str_copy_cstrn(const char* cstr, size_t len, struct md_allocator_i* alloc);

str_t str_printf(struct md_allocator_i* alloc, const char* format, ...);

// This should probably be removed
str_t load_textfile(str_t path, struct md_allocator_i* alloc);

// c:/folder/file.ext -> ext
bool extract_ext(str_t* ext, str_t path);

// c:/folder/file.ext -> file.ext
bool extract_file(str_t* file, str_t path);

// c:/folder/file.ext -> c:/folder/
bool extract_folder_path(str_t* folder_path, str_t path);

// Mutating string operations
void replace_char(char* str, size_t len, char c, char replacement);
void convert_to_lower(char* str, size_t len);
void convert_to_upper(char* str, size_t len);

// Copies the contents of a str_t into a char buffer and ensures zero termination
// Returns the number of characters written (excluding the zero termination character)
size_t str_copy_to_char_buf(char* buf, size_t cap, str_t str);

// Returns the Levenshtein edit distance between two strings
// https://en.wikipedia.org/wiki/Levenshtein_distance
int str_edit_distance(str_t a, str_t b);

#ifdef __cplusplus
}
#endif
