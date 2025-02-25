#include <core/md_str.h>

#include <core/md_common.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_parse.h>

#include <string.h>  // strlen, memchr
#include <stdarg.h>
#include <stdio.h>

bool str_eq_n(str_t str_a, str_t str_b, size_t n) {
    if (!str_a.ptr || !str_b.ptr) return false;
    if ((str_a.len < n || str_b.len < n) && str_a.len != str_b.len) return false;
    if (str_a.ptr[0] != str_b.ptr[0]) return false;
    // str_a & str_b have equal len at this point
    n = n < str_a.len ? n : str_a.len;
    return MEMCMP(str_a.ptr, str_b.ptr, n) == 0;
}


bool str_eq_ignore_case(const str_t str_a, const str_t str_b) {
    if (!str_a.ptr || !str_b.ptr) return false;
    if (str_a.len != str_b.len) return false;
    for (size_t i = 0; i < str_a.len; ++i) {
        if (to_lower(str_a.ptr[i]) != to_lower(str_b.ptr[i])) return false;
    }
    return true;
}

bool str_eq_n_ignore_case(const str_t str_a, const str_t str_b, size_t n) {
    if (!str_a.ptr || !str_b.ptr) return false;
	if ((str_a.len < n || str_b.len < n) && str_a.len != str_b.len) return false;

	// str_a & str_b have equal len.
	n = n < str_a.len ? n : str_a.len;
	for (size_t i = 0; i < n; ++i) {
    	if (to_lower(str_a.ptr[i]) != to_lower(str_b.ptr[i])) return false;
    }
	return true;
}

bool str_eq_cstr(str_t str, const char* cstr) {
    if (!str.ptr || !str.len || !cstr) return false;
    if (str.ptr[0] != cstr[0]) return false;
    return (strncmp(str.ptr, cstr, str.len) == 0) && cstr[str.len] == '\0';
}

// Compare str and cstr only up to n characters
bool str_eq_cstr_n(str_t str, const char* cstr, size_t n) {
    if (!n) return false;
    if (!str.ptr || !str.len || !cstr) return false;
    if (str.ptr[0] != cstr[0]) return false;
    n = n < str.len ? n : str.len;
    return strncmp(str.ptr, cstr, n) == 0;
}

bool str_eq_cstr_ignore_case(str_t str, const char* cstr) {
    if (!str.ptr || !str.len || !cstr) return false;
    for (size_t i = 0; i < str.len; ++i) {
        if (cstr[i] == '\0' || to_lower(str.ptr[i]) != to_lower(cstr[i])) return false;
    }
    return cstr[str.len] == '\0';
}

bool str_eq_cstr_n_ignore_case(str_t str, const char* cstr, size_t n) {
    if (!n) return false;
    if (!str.ptr || !str.len || !cstr) return false;
    n = n < str.len ? n : str.len;
    for (size_t i = 0; i < n; ++i) {
        if (cstr[i] == '\0' || to_lower(str.ptr[i]) != to_lower(cstr[i])) return false;
    }
    return true;
}

size_t str_count_equal_chars(str_t a, str_t b) {
    if (!a.ptr || a.len <= 0 || !b.ptr || b.len <= 0) return 0;
    const size_t len = MIN(a.len, b.len);
    size_t i = 0;
    for (; i < len; ++i) {
        if (a.ptr[i] != b.ptr[i]) break;
    }
    return i;
}

size_t str_count_occur_char(str_t str, char c) {
    if (!str.ptr || str.len <= 0) return 0;
    size_t count = 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (str.ptr[i] == c) count += 1;
    }
    return count;
}

str_t str_join(str_t first, str_t last) {
    ASSERT(first.ptr < last.ptr + last.len);
    str_t str = {
        str_beg(first),
        str_end(last) - str_beg(first)
    };
    return str;
}

int str_cmp_lex(str_t a, str_t b) {
    size_t len = a.len < b.len ? a.len : b.len;
    for (size_t i = 0; i < len; ++i) {
        if (a.ptr[i] < b.ptr[i]) return -1;
        if (a.ptr[i] > b.ptr[i]) return 1;
    }
    if (a.len < b.len) return -1;
    if (a.len > b.len) return 1;
    return 0;
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

bool str_peek_line(str_t* out_line, const str_t* in_str) {
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

bool str_extract_i32(int32_t* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = str_trim(tok);
    if (!is_int(tok)) return false;
    *val = (int)parse_int(tok);
    
    return true;
}

bool str_extract_i64(int64_t* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = str_trim(tok);
    if (!is_int(tok)) return false;
    *val = parse_int(tok);

    return true;
}

bool str_extract_f32(float* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = str_trim(tok);
    if (!is_float(tok)) return false;
    *val = (float)parse_float(tok);

    return true;
}

bool str_extract_f64(double* val, str_t* str) {
    ASSERT(val);
    ASSERT(str);

    str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = str_trim(tok);
    if (!is_float(tok)) return false;
    *val = parse_float(tok);

    return true;
}

bool str_find_char(size_t* loc, str_t str, int c) {
    // @TODO: Implement support for UTF encodings here.
    // Identify the effective width used in c and search for that.

    for (size_t i = 0; i < str.len; ++i) {
        if (str.ptr[i] == c) {
            if (loc) *loc = i;
            return true;
        }
    }
    return false;
}

// Insired by this
// https://github.com/Alexpux/Cygwin/blob/master/newlib/libc/string/memrchr.c#L62
// memrchr is a GNU extension
bool str_rfind_char(size_t* loc, str_t str, int c) {
    const char* ptr = str.ptr + str.len - 1;
    size_t len = str.len;

    // This should only be used until the pointer is aligned to 8-bytes
    while (true) {
        if (!len--) {
            return false;
        }
        if (*ptr == c) {
            if (loc) *loc = (size_t)(ptr - str.ptr);
            return true;
        }
        ptr--;
    }

    // @TODO: Check for character 8-characters at a time
}

bool str_find_str(size_t* loc, str_t haystack, str_t needle) {
    if (haystack.len == 0) return false;
    if (needle.len   == 0) return false;

    if (needle.len > haystack.len) {
        MD_LOG_ERROR("Trying to find 'needle' which is larger than supplied 'haystack'");
        return false;
    }

    const char* h_beg = haystack.ptr;
    const char* h_end = haystack.ptr + haystack.len;

    size_t i = 0;
    const char* n_beg = 0;
    for (const char* c = h_beg; c != h_end; ++c) {
        if (*c == needle.ptr[i]) {
            if (i == 0) n_beg = c;
            ++i;
        } else {
            i = 0;
        }

        if (i == needle.len) {
            if (loc) *loc = n_beg - h_beg;
            return true;
        }
    }

    return false;
}

bool str_begins_with(str_t str, str_t prefix) {
    return str.len >= prefix.len && MEMCMP(str.ptr, prefix.ptr, prefix.len) == 0;
}

bool str_ends_with(str_t str, str_t suffix) {
    return str.len >= suffix.len && MEMCMP(str.ptr + str.len - suffix.len, suffix.ptr, suffix.len) == 0;
}

str_t load_textfile(str_t filename, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
    str_t result = {0,0};
    if (file) {
        size_t file_size = md_file_size(file);
        char* mem = md_alloc(alloc, file_size + 1);
        if (mem) {
            size_t read_size = md_file_read(file, mem, file_size);
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

str_t str_alloc(size_t len, struct md_allocator_i* alloc) {
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

str_t str_copy(const str_t str, md_allocator_i* alloc) {
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

str_t str_copy_cstr(const char* cstr, md_allocator_i* alloc) {
    ASSERT(alloc);
    str_t result = {0,0};
    if (cstr) {
        size_t len = strlen(cstr);
        char* data = md_alloc(alloc, len + 1);
        data[len] = '\0';
        MEMCPY(data, cstr, len);
        result.ptr = data;
        result.len = len;
    }
    return result;
}

str_t str_copy_cstrn(const char* cstr, size_t len, md_allocator_i* alloc) {
    ASSERT(alloc);
    str_t result = {0,0};
    if (cstr) {
        char* data = md_alloc(alloc, len + 1);
        data[len] = '\0';
        MEMCPY(data, cstr, len);
        result.ptr = data;
        result.len = len;
    }
    return result;
}

str_t str_printf(struct md_allocator_i* alloc, const char* format, ...) {
    va_list args;

    va_start(args, format);
    int64_t len = vsnprintf(NULL, 0, format, args);
    va_end(args);

    if (len <= 0) {
        return (str_t){0};
    }

    char* buf = md_alloc(alloc, len + 1);

    va_start(args, format);
    vsnprintf(buf, len + 1, format, args);
    va_end(args);

    return (str_t) {buf, len};
}

// c:/folder/file.ext -> ext
bool extract_ext(str_t* ext, str_t path) {
    size_t loc;
    if (str_rfind_char(&loc, path, '.')) {
        if (ext) {
            loc += 1; // skip '.'
            ext->ptr = path.ptr + loc;
            ext->len = path.len - loc;
        }
        return true;
    }
    return false;
}

// c:/folder/file.ext -> file.ext
bool extract_file(str_t* file, str_t path) {
    size_t loc;
    if (!str_rfind_char(&loc, path, '/') && !str_rfind_char(&loc, path, '\\')) {
        return false;
    }

    if (file) {
        loc += 1; // skip slash or backslash
        file->ptr = path.ptr + loc;
        file->len = path.len - loc;
    }
    return true;
}

// c:/folder/file.ext -> c:/folder/file
bool extract_file_path_without_ext(str_t* file_path, str_t path) {
    size_t loc;
    if (!str_rfind_char(&loc, path, '.')) {
        return false;
    }

    if (file_path) {
        file_path->ptr = path.ptr;
        file_path->len = loc;
    }
    return true;
}

// c:/folder/file.ext -> c:/folder/
bool extract_folder_path(str_t* folder_path, str_t path) {
    size_t loc;
    if (!str_rfind_char(&loc, path, '/') && !str_rfind_char(&loc, path, '\\')) {
        return false;
    }

    if (folder_path) {
        folder_path->ptr = path.ptr;
        folder_path->len = loc+1;    // include '/' or '\'
    }
    return true;
}

void replace_char(char* str, size_t len, char target, char replacement) {
    if (!str) return;
	for (char* c = str; c != str + len; ++c) {
    	if (*c == target) *c = replacement;
    }
}

void convert_to_lower(char* str, size_t len) {
    ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        *c = (char)to_lower(*c);
    }
}

void convert_to_upper(char* str, size_t len) {
    ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        *c = (char)to_upper(*c);
    }
}

size_t str_copy_to_char_buf(char* buf, size_t cap, str_t str) {
    ASSERT(buf);
    if (cap == 0 || str_empty(str)) return 0;
    const size_t len = CLAMP(str.len, 0, cap - 1);
    MEMCPY(buf, str.ptr, len);
    buf[len] = '\0';
    return len;
}

// Ported from Martin Ettl's version available on
// https://rosettacode.org/wiki/Levenshtein_distance#C++

int str_edit_distance(str_t s1, str_t s2) {
    const int m = (int)str_len(s1);
    const int n = (int)str_len(s2);

    if (m == 0) return n;
    if (n == 0) return m;

    const int64_t bytes = ALIGN_TO(sizeof(int) * (n + 1), 16);
    int* costs = md_temp_push(bytes);

    for (int k=0; k<=n; k++) {
        // Seller's variant (=0) for fuzzy matching
        costs[k] = k;
    }

    int i = 0;
    for (const char* c1 = str_beg(s1); c1 != str_end(s1); ++c1) {
        costs[0] = i+1;
        int corner = i;
        int j = 0;

        for (const char* c2 = str_beg(s2); c2 != str_end(s2); ++c2) {
            int upper = costs[j+1];
            if (*c1 == *c2) {
                costs[j+1] = corner;
            } else {
                int t = MIN(upper, corner);
                costs[j+1] = MIN(costs[j], t) + 1;
            }

            corner = upper;
            j++;
        }
        i++;
    }

    int result = costs[n];
    md_temp_pop(bytes);

    return result;
}
