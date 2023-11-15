#pragma once

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
    int64_t     len;

#ifdef __cplusplus
    constexpr char operator[](int64_t idx) noexcept       { return ptr[idx]; }
    constexpr char operator[](int64_t idx) const noexcept { return ptr[idx]; }

    constexpr const char* beg() noexcept { return ptr; }
    constexpr const char* end() noexcept { return ptr + len; }
    constexpr operator bool() noexcept { return ptr && len > 0; }
#endif
} str_t;

// Macro to bake a string literal into a str_t
#ifdef __cplusplus
#define STR(cstr) {(cstr""), (sizeof(cstr)-1)}
#else
#define STR(cstr) (str_t){(cstr""), (sizeof(cstr)-1)}
#endif

#define STR_FMT(str) (int)str.len, str.ptr

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

static inline void str_swap(str_t a, str_t b) {
    str_t tmp = a;
	a = b;
	b = tmp;
}

static inline bool str_empty(str_t str) {
    return str.ptr == 0 || str.len == 0;
}

static inline const char* str_beg(str_t str) {
    return str.ptr;
}

static inline const char* str_end(str_t str) {
    return str.ptr + str.len;
}

static inline int64_t str_len(str_t str) {
    return str.len;
}

static inline const char* str_ptr(str_t str) {
    return str.ptr;
}

str_t str_trim(str_t str);

bool str_equal(str_t str_a, str_t str_b);
bool str_equal_n(str_t str_a, str_t str_b, int64_t n);
bool str_equal_ignore_case(const str_t str_a, const str_t str_b);

// Lexicographical comparison
int str_compare_lex(str_t str_a, str_t str_b);

bool str_equal_cstr(str_t str, const char* cstr);
bool str_equal_cstr_n(str_t str, const char* cstr, int64_t n);
bool str_equal_cstr_ignore_case(str_t str, const char* cstr);

int64_t str_count_equal_chars(str_t a, str_t b);
int64_t str_count_occur_char(str_t str, char character);

str_t str_substr(str_t str, int64_t offset, int64_t length DEF_VAL(-1));

bool str_skip_line(str_t* in_out_str);
bool str_peek_line(str_t* out_line, const str_t* in_str);
bool str_extract_line(str_t* out_line, str_t* in_out_str);

bool str_extract_i32  (int* val, str_t* in_out_str);
bool str_extract_i64(int64_t* val, str_t* in_out_str);

bool str_extract_f32(float* val,  str_t* in_out_str);
bool str_extract_f64(double* val,  str_t* in_out_str);

// Returns the offset where the item was found or -1 if nothing was found.
int64_t str_find_char(str_t str, int c);
int64_t str_rfind_char(str_t str, int c);
int64_t str_find_str(str_t str, str_t str_to_find);

bool str_starts_with(str_t str, str_t prefix);
bool str_ends_with(str_t str, str_t suffix);

// Will allocate one extra character for zero termination
str_t str_alloc(uint64_t len, struct md_allocator_i* alloc);
void  str_free(str_t str, struct md_allocator_i* alloc);
str_t str_copy(str_t str, struct md_allocator_i* alloc);
str_t str_copy_cstr(const char* cstr, struct md_allocator_i* alloc);
str_t str_copy_cstrn(const char* cstr, int64_t len, struct md_allocator_i* alloc);

// This should probably be removed
str_t load_textfile(str_t path, struct md_allocator_i* alloc);
str_t alloc_printf(struct md_allocator_i* alloc, const char* format, ...);

// c:/folder/file.ext -> ext
str_t extract_ext(str_t path);

// c:/folder/file.ext -> file.ext
str_t extract_file(str_t path);

// c:/folder/file.ext -> c:/folder/file
str_t extract_path_without_ext(str_t path);

// c:/folder/file.ext -> c:/folder/
str_t extract_path_without_file(str_t path);

// Converts Windows backslashes '\\' to forward slashes '/'
static inline void convert_backslashes(char* str, int64_t len) {
    ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        if (*c == '\\') *c = '/';
    }
}

static inline void convert_to_lower(char* str, int64_t len) {
    ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        *c = (char)to_lower(*c);
    }
}

static inline void convert_to_upper(char* str, int64_t len) {
    ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        *c = (char)to_upper(*c);
    }
}

// Copies the contents of a str_t into a char buffer and ensures zero termination
// Returns the number of characters written (excluding the zero termination character)
int64_t str_copy_to_char_buf(char* buf, int64_t cap, str_t str);

// Returns the Levenshtein edit distance between two strings
// https://en.wikipedia.org/wiki/Levenshtein_distance
int str_edit_distance(str_t a, str_t b);

#ifdef __cplusplus
}
#endif
