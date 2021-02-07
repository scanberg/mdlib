#include "str_util.h"

#include "common.h"
#include "file.h"
#include <md_allocator.h>
#include <md_log.h>

#include <string.h>

str_t substr(str_t str, int offset, int length) {
    str.str = str.str + offset;
    str.len = length < (str.len - offset) ? length : (str.len - offset);
    return str;
}

bool skip_line(str_t* in_out_str) {
    ASSERT(in_out_str);
    const char* c = (const char*)memchr(in_out_str->str, '\n', in_out_str->len);
    if (c) {
        in_out_str->len = in_out_str->len - (c - in_out_str->str);
        in_out_str->str = c + 1;
        return in_out_str;
    }
    return false;
}

bool peek_line(str_t* out_line, const str_t* in_str) {
    ASSERT(out_line);
    ASSERT(in_str);
    const char* beg = in_str->str;
    const char* end = (const char*)memchr(in_str->str, '\n', in_str->len);
    if (end) {
        out_line->str = beg;
        out_line->len = end - beg;
        return true;
    }
    return false;
}

bool extract_line(str_t* out_line, str_t* in_out_str) {
    ASSERT(out_line);
    ASSERT(in_out_str);
    const char* beg = in_out_str->str;
    if (skip_line(in_out_str)) {
        out_line->str = beg;
        out_line->len = in_out_str->str - beg;
    }
    return false;
}

str_t trim_whitespace(str_t str) {
    const char* beg = str.str;
    const char* end = str.str + str.len;
    while (beg != end && is_whitespace(*beg)) ++beg;
    while (beg + 1 < end && is_whitespace(end[-1])) --end;
    str.str = beg;
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
    const char* c = str.str;
    const char* end = str.str + str.len;
    while (c != end && is_digit(*c)) {
        val = val * 10 + (*c - '0');
        ++c;
    }
    if (*c != '.') return val;

    ++c; // skip '.'
    const uint32_t count = (uint32_t)(end - c);
    while (c < end) {
        val = val * 10 + (*c - '0');
        ++c;
    }

    return val / pow10[count];
}

int64_t parse_int(str_t str) {
    int64_t val = 0;
    const char* c = str.str;
    const char* end = str.str + str.len;
    while (c != end && is_digit(*c)) {
        val = val * 10 + (*c - '0');
        ++c;
    }
    return val;
}

str_t load_textfile(const char* filename, uint32_t filename_len, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_file* file = md_file_open(filename, filename_len, "rb");
    str_t result = {0,0};
    if (file) {
        uint64_t file_size = md_file_size(file);
        char* mem = md_alloc(alloc, file_size + 1);
        if (mem) {
            uint64_t read_size = md_file_read(file, mem, file_size);
            md_file_close(file);

            ASSERT(read_size == file_size);
            mem[file_size] = '\0'; // Zero terminate as a nice guy

            result.str = mem;
            result.len = file_size;
        } else {
            md_printf(MD_LOG_TYPE_ERROR, "Could not allocate memory for file %d", 10);
        }
        md_file_close(file);
    }
    return result;
}