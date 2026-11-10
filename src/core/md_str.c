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

size_t md_str_count_equal_chars(md_str_t a, md_str_t b) {
    if (!a.ptr || a.len <= 0 || !b.ptr || b.len <= 0) return 0;
    const size_t len = MD_MIN(a.len, b.len);
    size_t i = 0;
    for (; i < len; ++i) {
        if (a.ptr[i] != b.ptr[i]) break;
    }
    return i;
}

size_t md_str_count_occur_char(md_str_t str, char c) {
    if (!str.ptr || str.len <= 0) return 0;
    size_t count = 0;
    for (size_t i = 0; i < str.len; ++i) {
        if (str.ptr[i] == c) count += 1;
    }
    return count;
}

bool md_str_extract_i32(int32_t* val, md_str_t* str) {
    MD_ASSERT(val);
    MD_ASSERT(str);

    md_str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = md_str_trim(tok);
    if (!is_int(tok)) return false;
    *val = (int)parse_int(tok);
    
    return true;
}

bool md_str_extract_i64(int64_t* val, md_str_t* str) {
    MD_ASSERT(val);
    MD_ASSERT(str);

    md_str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = md_str_trim(tok);
    if (!is_int(tok)) return false;
    *val = parse_int(tok);

    return true;
}

bool md_str_extract_f32(float* val, md_str_t* str) {
    MD_ASSERT(val);
    MD_ASSERT(str);

    md_str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = md_str_trim(tok);
    if (!is_float(tok)) return false;
    *val = (float)parse_float(tok);

    return true;
}

bool md_str_extract_f64(double* val, md_str_t* str) {
    MD_ASSERT(val);
    MD_ASSERT(str);

    md_str_t tok;
    if (!extract_token(&tok, str)) return false;
    tok = md_str_trim(tok);
    if (!is_float(tok)) return false;
    *val = parse_float(tok);

    return true;
}

bool md_str_find_str(size_t* loc, md_str_t haystack, md_str_t needle) {
    if (needle.len == 0 || needle.len > haystack.len) return false;

    // Every position where the needle still fits is a candidate. (Restarting the match after the character that broke
    // a partial match is not enough: "aab" contains "ab" starting inside the partial match "aa".)
    const char* h_beg = haystack.ptr;
    const char* last  = haystack.ptr + (haystack.len - needle.len);
    for (const char* c = h_beg; c <= last; ++c) {
        c = (const char*)memchr(c, needle.ptr[0], (size_t)(last - c) + 1);
        if (!c) break;
        if (md_memcmp(c, needle.ptr, needle.len) == 0) {
            if (loc) *loc = (size_t)(c - h_beg);
            return true;
        }
    }

    return false;
}

md_str_t md_load_textfile(md_str_t filename, struct md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    md_file_t file = {0};
    md_str_t result = {0};
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

md_str_t md_str_alloc(size_t len, struct md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    char* mem = md_alloc(alloc, len + 1);
    md_memset(mem, 0, len + 1);
    md_str_t str = {mem, len};
    return str;
}

void md_str_free(md_str_t str, struct md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    md_free(alloc, (void*)str.ptr, str.len + 1);
}

md_str_t md_str_copy(const md_str_t str, md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    md_str_t result = {0,0};
    if (str.ptr && str.len > 0) {
        char* data = md_alloc(alloc, str.len + 1);
        data[str.len] = '\0';
        md_memcpy(data, str.ptr, str.len);
        result.ptr = data;
        result.len = str.len;
    }
    return result;
}

md_str_t md_str_copy_cstr(const char* cstr, md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    md_str_t result = {0,0};
    if (cstr) {
        size_t len = strlen(cstr);
        char* data = md_alloc(alloc, len + 1);
        data[len] = '\0';
        md_memcpy(data, cstr, len);
        result.ptr = data;
        result.len = len;
    }
    return result;
}

md_str_t md_str_copy_cstrn(const char* cstr, size_t len, md_allocator_i* alloc) {
    MD_ASSERT(alloc);
    md_str_t result = {0,0};
    if (cstr) {
        char* data = md_alloc(alloc, len + 1);
        data[len] = '\0';
        md_memcpy(data, cstr, len);
        result.ptr = data;
        result.len = len;
    }
    return result;
}

md_str_t md_str_printf(struct md_allocator_i* alloc, const char* format, ...) {
    va_list args;

    va_start(args, format);
    int64_t len = md_vsnprintf(NULL, 0, format, args);
    va_end(args);

    if (len <= 0) {
        return (md_str_t){0};
    }

    char* buf = md_alloc(alloc, len + 1);

    va_start(args, format);
    md_vsnprintf(buf, len + 1, format, args);
    va_end(args);

    return (md_str_t) {buf, len};
}

// c:/folder/file.ext -> ext
bool md_extract_ext(md_str_t* ext, md_str_t path) {
    // The last '.' of the file name, not of the path: "/data/run.1/mdout" has no extension
    size_t name_beg = 0;
    for (size_t i = path.len; i > 0; --i) {
        if (path.ptr[i - 1] == '/' || path.ptr[i - 1] == '\\') {
            name_beg = i;
            break;
        }
    }
    size_t loc;
    if (md_str_rfind_char(&loc, md_str_substr(path, name_beg, SIZE_MAX), '.')) {
        if (ext) {
            loc += name_beg + 1; // skip '.'
            ext->ptr = path.ptr + loc;
            ext->len = path.len - loc;
        }
        return true;
    }
    return false;
}

// c:/folder/file.ext -> file.ext
bool md_extract_file(md_str_t* file, md_str_t path) {
    size_t loc;
    if (!md_str_rfind_char(&loc, path, '/') && !md_str_rfind_char(&loc, path, '\\')) {
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
bool md_extract_file_path_without_ext(md_str_t* file_path, md_str_t path) {
    size_t loc;
    if (!md_str_rfind_char(&loc, path, '.')) {
        return false;
    }

    if (file_path) {
        file_path->ptr = path.ptr;
        file_path->len = loc;
    }
    return true;
}

// c:/folder/file.ext -> c:/folder/
bool md_extract_folder_path(md_str_t* folder_path, md_str_t path) {
    size_t loc = path.len;
    if (!md_str_rfind_char(&loc, path, '/') && !md_str_rfind_char(&loc, path, '\\')) {
        return false;
    }

    if (folder_path) {
        folder_path->ptr = path.ptr;
        folder_path->len = loc+1;    // include '/' or '\'
    }

    return true;
}

void md_replace_char(char* str, size_t len, char target, char replacement) {
    if (!str) return;
	for (char* c = str; c != str + len; ++c) {
    	if (*c == target) *c = replacement;
    }
}

void md_convert_to_lower(char* str, size_t len) {
    MD_ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        *c = (char)to_lower(*c);
    }
}

void md_convert_to_upper(char* str, size_t len) {
    MD_ASSERT(str);
    for (char* c = str; c != str + len; ++c) {
        *c = (char)to_upper(*c);
    }
}

// Ported from Martin Ettl's version available on
// https://rosettacode.org/wiki/Levenshtein_distance#C++

int md_str_edit_distance(md_str_t s1, md_str_t s2) {
    const int m = (int)md_str_len(s1);
    const int n = (int)md_str_len(s2);

    if (m == 0) return n;
    if (n == 0) return m;

    md_temp_scope_t temp = md_temp_begin();
    int* costs = md_temp_alloc_array(temp, int, n + 1);

    for (int k=0; k<=n; k++) {
        // Seller's variant (=0) for fuzzy matching
        costs[k] = k;
    }

    int i = 0;
    for (const char* c1 = md_str_beg(s1); c1 != md_str_end(s1); ++c1) {
        costs[0] = i+1;
        int corner = i;
        int j = 0;

        for (const char* c2 = md_str_beg(s2); c2 != md_str_end(s2); ++c2) {
            int upper = costs[j+1];
            if (*c1 == *c2) {
                costs[j+1] = corner;
            } else {
                int t = MD_MIN(upper, corner);
                costs[j+1] = MD_MIN(costs[j], t) + 1;
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
