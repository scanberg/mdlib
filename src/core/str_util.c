#include "str_util.h"

#include "common.h"
#include "file.inl"
#include <md_allocator.h>
#include <md_log.h>

#include <string.h>

bool compare_str(const str_t str_a, const str_t str_b) {
    ASSERT(str_a.ptr && str_b.ptr);
    if (str_a.len != str_b.len) return false;
    for (uint64_t i = 0; i < str_a.len; ++i) {
        if (str_a.ptr[i] != str_b.ptr[i]) return false;
    }
    return true;
}

str_t substr(str_t str, uint64_t offset, uint64_t length) {
    if (offset > str.len) return (str_t){0};
    if (offset + length > str.len) length = str.len - offset;
    str.ptr = str.ptr + offset;
    str.len = length;
    return str;
}

bool skip_line(str_t* in_out_str) {
    ASSERT(in_out_str);
    const char* c = (const char*)memchr(in_out_str->ptr, '\n', in_out_str->len);
    if (c) {
        in_out_str->len = in_out_str->len - (c - in_out_str->ptr);
        in_out_str->ptr = c + 1;
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

bool extract_line(str_t* out_line, str_t* in_out_str) {
    ASSERT(out_line);
    ASSERT(in_out_str);
    const char* beg = in_out_str->ptr;
    if (skip_line(in_out_str)) {
        out_line->ptr = beg;
        out_line->len = in_out_str->ptr - beg;
    }
    return false;
}

str_t trim_whitespace(str_t str) {
    const char* beg = str.ptr;
    const char* end = str.ptr + str.len;
    while (beg != end && is_whitespace(*beg)) ++beg;
    while (beg + 1 < end && is_whitespace(end[-1])) --end;
    str.ptr = beg;
    str.len = end - beg;
    return str;
}

double parse_float(str_t str) {
    const double pow10[16] = {
        1e+0,  1e+1,  1e+2,  1e+3,
        1e+4,  1e+5,  1e+6,  1e+7,
        1e+8,  1e+9,  1e+10, 1e+11,
        1e+12, 1e+13, 1e+14, 1e+15
    };

    double val = 0;
    double sign = 1;
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;
    if (*c == '-') {
        ++c;
        sign = -1;
    }
    while (c != end && is_digit(*c)) {
        val = val * 10 + (*c - '0');
        ++c;
    }
    if (*c != '.' || c == end) return sign * val;

    ++c; // skip '.'
    const uint32_t count = (uint32_t)(end - c);
    while (c < end) {
        val = val * 10 + (*c - '0');
        ++c;
    }

    return sign * val / pow10[count];
}

int64_t parse_int(str_t str) {
    int64_t val = 0;
    const char* c = str.ptr;
    const char* end = str.ptr + str.len;
    int64_t sign = 1;
    if (*c == '-') {
        ++c;
        sign = -1;
    }
    while (c != end && is_digit(*c)) {
        val = val * 10 + (*c - '0');
        ++c;
    }
    return sign * val;
}

str_t load_textfile(const char* filename, uint32_t filename_len, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    FILE* file = md_file_open(filename, filename_len, "rb");
    str_t result = {0,0};
    if (file) {
        uint64_t file_size = md_file_size(file);
        char* mem = md_alloc(alloc, file_size + 1);
        if (mem) {
            uint64_t read_size = md_file_read(file, mem, file_size);
            md_file_close(file);

            ASSERT(read_size == file_size);
            mem[file_size] = '\0'; // Zero terminate as a nice guy

            result.ptr = mem;
            result.len = file_size;
        } else {
            md_printf(MD_LOG_TYPE_ERROR, "Could not allocate memory for file %d", 10);
        }
        md_file_close(file);
    }
    return result;
}

str_t alloc_str(uint64_t len, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    char* mem = md_alloc(alloc, len + 1);
    memset(mem, 0, len + 1);
    str_t str = {mem, len};
    return str;
}

void free_str(str_t str, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_free(alloc, (void*)str.ptr, str.len + 1);
}

str_t copy_str(const str_t str, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    str_t result = {0,0};
    if (str.ptr && str.len > 0) {
        char* data = md_alloc(alloc, str.len + 1);
        data[str.len] = '\0';
        memcpy(data, str.ptr, str.len);
        result.ptr = data;
        result.len = str.len;
    }
    return result;
}