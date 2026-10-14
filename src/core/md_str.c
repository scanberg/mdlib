#include <core/md_str.h>

#include <core/md_common.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_parse.h>

#include <string.h>  // strlen, memchr
#include <stdarg.h>
#include <stdio.h>

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

bool str_find_str(size_t* loc, str_t haystack, str_t needle) {
    if (needle.len == 0 || needle.len > haystack.len) return false;

    // Every position where the needle still fits is a candidate. (Restarting the match after the character that broke
    // a partial match is not enough: "aab" contains "ab" starting inside the partial match "aa".)
    const char* h_beg = haystack.ptr;
    const char* last  = haystack.ptr + (haystack.len - needle.len);
    for (const char* c = h_beg; c <= last; ++c) {
        c = (const char*)memchr(c, needle.ptr[0], (size_t)(last - c) + 1);
        if (!c) break;
        if (MEMCMP(c, needle.ptr, needle.len) == 0) {
            if (loc) *loc = (size_t)(c - h_beg);
            return true;
        }
    }

    return false;
}

str_t load_textfile(str_t filename, struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_file_t file = {0};
    str_t result = {0};
    if (md_file_open(&file, filename, MD_FILE_READ)) {
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
        md_file_close(&file);
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
    size_t loc = path.len;
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

// Ported from Martin Ettl's version available on
// https://rosettacode.org/wiki/Levenshtein_distance#C++

int str_edit_distance(str_t s1, str_t s2) {
    const int m = (int)str_len(s1);
    const int n = (int)str_len(s2);

    if (m == 0) return n;
    if (n == 0) return m;

    md_temp_scope_t temp = md_temp_begin();
    int* costs = md_temp_alloc_array(temp, int, n + 1);

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
    md_temp_end(temp);

    return result;
}
