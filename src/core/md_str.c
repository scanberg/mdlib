#include "md_str.h"

#include "md_common.h"
#include "md_os.h"
#include "md_allocator.h"
#include "md_array.h"
#include "md_log.h"

#include <string.h>
#include <stdarg.h>
#include <stdio.h>

str_t str_from_cstr(const char* cstr) {
    str_t str = {cstr, strlen(cstr)};
    return str;
}

bool str_skip_line(str_t* in_out_str) {
    ASSERT(in_out_str);
    const char* c = (const char*)memchr(in_out_str->ptr, '\n', in_out_str->len);
    if (c) {
        c = c + 1;
        in_out_str->len = in_out_str->len - (c - in_out_str->ptr);
        in_out_str->ptr = c;
        return in_out_str;
    }
    return false;
}

bool peek_line(str_t* out_line, const str_t* in_str) {
    ASSERT(out_line);
    ASSERT(in_str);
    const char* beg = in_str->ptr;
    const char* end = (const char*)memchr(in_str->ptr, '\n', in_str->len);
    if (end) {
        out_line->ptr = beg;
        out_line->len = end - beg;
        return true;
    }
    return false;
}

bool str_extract_line(str_t* out_line, str_t* in_out_str) {
    ASSERT(out_line);
    ASSERT(in_out_str);
    if (in_out_str->len == 0) return false;

    const char* beg = in_out_str->ptr;
    const char* end = in_out_str->ptr + in_out_str->len;
    const char* c = (const char*)memchr(in_out_str->ptr, '\n', in_out_str->len);
    if (c) {
        end = c + 1;
        if (c > beg && c[-1] == '\r') --c;
    } else {
        c = end;
    }

    out_line->ptr = beg;
    out_line->len = MAX(0, c - beg);

    in_out_str->len = in_out_str->ptr + in_out_str->len - end;
    in_out_str->ptr = end;
    
    return true;
}

int64_t str_read_line(str_t* str, char* buf, int64_t cap) {
    ASSERT(str);
    if (!str->ptr || str->len == 0) return 0;
    const char* beg = str_beg(*str);
    const char* end = str_end(*str);
    const char* c = (const char*)memchr(str->ptr, '\n', str->len);
    if (c) {
        end = c + 1;
    }
    end = MIN(end, beg + cap);
    int64_t len = end - beg;
    if (len > 0) {
        MEMCPY(buf, beg, len - 1);
        buf[len - 1] = '\0';
    }

    str->ptr = end;
    str->len = str->len - len;
    return len - 1;
}

bool str_extract_i32(int* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_next_token(&tok, str)) return false;
    *val = (int)parse_int(str_trim(tok));
    
    return true;
}

bool str_extract_i64(int64_t* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_next_token(&tok, str)) return false;
    *val = parse_int(str_trim(tok));

    return true;
}

bool str_extract_f32(float* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_next_token(&tok, str)) return false;
    *val = (float)parse_float(str_trim(tok));

    return true;
}

bool str_extract_f64(double* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_next_token(&tok, str)) return false;
    *val = parse_float(str_trim(tok));

    return true;
}

int64_t str_find_char(str_t str, int c) {
    // @TODO: Implement support for UTF encodings here.
    // Identify the effective width used in c and search for that.
    for (int64_t i = 0; i < (int64_t)str.len; ++i) {
        if (str.ptr[i] == c) return i;
    }
    return -1;
}

int64_t str_rfind_char(str_t str, int c) {
    for (int64_t i = str.len - 1; i >= 0; --i) {
        if (str.ptr[i] == c) return i;
    }
    return -1;
}

bool str_starts_with(str_t str, str_t prefix) {
    return str.len >= prefix.len && strncmp(str.ptr, prefix.ptr, prefix.len) == 0;
}

bool str_ends_with(str_t str, str_t suffix) {
    return str.len >= suffix.len && strncmp(str.ptr + str.len - suffix.len, suffix.ptr, suffix.len) == 0;
}

str_t str_find_str(str_t haystack, str_t needle) {
    const char* h_beg = haystack.ptr;
    const char* h_end = haystack.ptr + haystack.len;
    str_t result = {0};

    if (haystack.len == 0) goto done;
    if (needle.len == 0) goto done;
    if (needle.len > haystack.len) {
        MD_LOG_ERROR("Trying to find 'needle' which is larger than supplied 'haystack'");
        goto done;
    }

    int64_t i = 0;
    const char* n_beg = 0;
    for (const char* c = h_beg; c != h_end; ++c) {
        if (*c == needle.ptr[i]) {
            if (i == 0) n_beg = c;
            ++i;
        } else {
            i = 0;
        }

        if (i == needle.len) {
            result.ptr = n_beg;
            result.len = needle.len;
            goto done;
        }
    }

done:
    return result;
}

#if MD_COMPILER_MSVC
    double __cdecl pow(double x, double y);
#   pragma intrinsic(pow)
#   define POW pow
#elif MD_COMPILER_GCC || MD_COMPILER_CLANG
#   define POW __builtin_pow
#endif

double parse_float(str_t str) {
    ASSERT(str.ptr);
    static const double pow10[32] = {
        1e+0,  1e+1,  1e+2,  1e+3,  1e+4,  1e+5,  1e+6,  1e+7,
        1e+8,  1e+9,  1e+10, 1e+11, 1e+12, 1e+13, 1e+14, 1e+15,
        1e+16, 1e+17, 1e+18, 1e+19, 1e+20, 1e+21, 1e+22, 1e+23,
        1e+24, 1e+25, 1e+26, 1e+27, 1e+28, 1e+29, 1e+30, 1e+31,
    };

    const char* c = str.ptr;
    const char* end = str.ptr + str.len;

    while (c < end && is_whitespace(*c)) ++c;

    double val = 0;
    double sign = 1;
    if (*c == '-') {
        ++c;
        sign = -1;
    }
    while (c < end && is_digit(*c)) {
        val = val * 10 + ((int)(*c) - '0');
        ++c;
    }

    if (c < end && *c == '.') {
        const char* dec = ++c;
        while (c < end && (c - dec) < (int64_t)ARRAY_SIZE(pow10) - 1 && is_digit(*c)) {
            val = val * 10 + ((int)(*c) - '0');
            ++c;
        }
        int64_t count = c - dec;
        val /= pow10[count];
    }

    if (c < end && (*c == 'e' || *c == 'E')) {
        ++c;
        int exp_sign = 1;
        if (c < end && (*c == '+' || *c == '-')) {
            exp_sign = *c == '+' ? 1 : -1;
            ++c;
        }
        int exp_val = 0;
        while (c < end && is_digit(*c)) {
            exp_val = exp_val * 10 + ((int)(*c) - '0');
            ++c;
        }
        while (exp_val) {
            int ev = MIN(exp_val, (int)ARRAY_SIZE(pow10)-1);
            val = exp_sign > 0 ? val * pow10[ev] : val / pow10[ev];
            exp_val -= ev;
        }
    }

    return sign * val;
}

int64_t parse_int(str_t str) {
    int64_t val = 0;
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;
    while (c < end && is_whitespace(*c)) ++c;
    int64_t sign = 1;
    if (*c == '-') {
        ++c;
        sign = -1;
    }
    while (c != end && is_digit(*c)) {
        val = val * 10 + ((int)(*c) - '0');
        ++c;
    }
    return sign * val;
}

str_t load_textfile(str_t filename, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    str_t result = {0,0};
    if (file) {
        int64_t file_size = md_file_size(file);
        char* mem = md_alloc(alloc, file_size + 1);
        if (mem) {
            int64_t read_size = md_file_read(file, mem, file_size);
            if (read_size == file_size) {
                mem[file_size] = '\0'; // Zero terminate as a nice guy
                result.ptr = mem;
                result.len = file_size;
            }
            else {
                MD_LOG_ERROR("Failed to read full textfile");
            }
        } else {
            MD_LOG_ERROR("Could not allocate memory for file %d", 10);
        }
        md_file_close(file);
    }
    return result;
}

str_t alloc_str(uint64_t len, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    char* mem = md_alloc(alloc, len + 1);
    MEMSET(mem, 0, len + 1);
    str_t str = {mem, len};
    return str;
}

void str_free(str_t str, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_free(alloc, (void*)str.ptr, str.len + 1);
}

str_t str_copy(const str_t str, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    str_t result = {0,0};
    if (str.ptr && str.len > 0) {
        char* data = md_alloc(alloc, str.len + 1);
        data[str.len] = '\0';
        MEMCPY(data, str.ptr, str.len);
        result.ptr = data;
        result.len = str.len;
    }
    return result;
}

str_t alloc_printf(struct md_allocator_i* alloc, const char* format, ...) {
    va_list args;

    va_start(args, format);
    int len = vsnprintf(NULL, 0, format, args);
    va_end(args);

    char* buf = md_alloc(alloc, len + 1);

    va_start(args, format);
    vsnprintf(buf, len + 1, format, args);
    va_end(args);

    return (str_t) {buf, len};
}

str_t str_concat(str_t a, str_t b, struct md_allocator_i* alloc) {
    char* buf = md_alloc(alloc, a.len + b.len + 1);
    int64_t len = a.len + b.len;
    
    if (!str_empty(a)) {
        MEMCPY(buf, a.ptr, a.len);
    }
    if (!str_empty(b)) {
        MEMCPY(buf + a.len, b.ptr, b.len);
    }
    buf[a.len + b.len] = '\0';
    return (str_t){buf, len};
}

// c:/folder/file.ext -> ext
str_t extract_ext(str_t path) {
    int64_t pos = str_rfind_char(path, '.');

    str_t res = {0,0};
    if (pos > -1) {
        pos += 1; // skip '.'
        res.ptr = path.ptr + pos;
        res.len = path.len - pos;
    }
    return res;
}

// c:/folder/file.ext -> file.ext
str_t extract_file(str_t path) {
    int64_t pos = str_rfind_char(path, '/');
    if (pos == -1) {
        pos = str_rfind_char(path, '\\');
    }

    str_t res = {0,0};
    if (pos != -1) {
        pos += 1; // skip slash or backslash
        res.ptr = path.ptr + pos;
        res.len = path.len - pos;
    }
    return res;
}

// c:/folder/file.ext -> c:/folder/file.
str_t extract_path_without_ext(str_t path) {
    const int64_t pos = str_rfind_char(path, '.');

    str_t res = {0,0};
    if (pos) {
        res.ptr = path.ptr;
        res.len = pos;
    }
    return res;
}

// c:/folder/file.ext -> c:/folder/
str_t extract_path_without_file(str_t path) {
    int64_t pos = str_rfind_char(path, '/');
    if (pos == -1) {
        pos = str_rfind_char(path, '\\');
    }

    str_t res = {0};
    if (pos != -1) {
        res.ptr = path.ptr;
        res.len = pos+1;    // include '/' or '\'
    }
    return res;
}

bool extract_next_token(str_t* tok, str_t* str) {
    ASSERT(tok);
    ASSERT(str);
    if (!str->ptr || str->len == 0) return false;

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

bool extract_next_token_delim(str_t* tok, str_t* str, char delim) {
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
