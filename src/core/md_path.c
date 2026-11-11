#include <core/md_path.h>

static inline bool is_sep(char c) {
    return c == '/' || c == '\\';
}

// Offset of the file name within path: just past the last separator, 0 if there is none
static size_t filename_beg(md_str_t path) {
    for (size_t i = path.len; i > 0; --i) {
        if (is_sep(path.ptr[i - 1])) {
            return i;
        }
    }
    return 0;
}

// Offset of the '.' which starts the extension: searched for in the file name only, and not at its first
// character, as a leading '.' marks a hidden file
static bool ext_dot(size_t* loc, md_str_t path) {
    const size_t beg = filename_beg(path);
    for (size_t i = path.len; i > beg + 1; --i) {
        if (path.ptr[i - 1] == '.') {
            *loc = i - 1;
            return true;
        }
    }
    return false;
}

md_str_t md_path_filename(md_str_t path) {
    md_str_t result = {0};
    if (!md_str_empty(path)) {
        const size_t beg = filename_beg(path);
        result.ptr = path.ptr + beg;
        result.len = path.len - beg;
    }
    return result;
}

md_str_t md_path_folder(md_str_t path) {
    md_str_t result = {0};
    const size_t beg = filename_beg(path);
    if (beg > 0) {
        result.ptr = path.ptr;
        result.len = beg;
    }
    return result;
}

md_str_t md_path_ext(md_str_t path) {
    md_str_t result = {0};
    size_t dot;
    if (ext_dot(&dot, path)) {
        result.ptr = path.ptr + dot + 1;
        result.len = path.len - dot - 1;
    }
    return result;
}

md_str_t md_path_without_ext(md_str_t path) {
    size_t dot;
    if (ext_dot(&dot, path)) {
        path.len = dot;
    }
    return path;
}

md_str_t md_path_stem(md_str_t path) {
    return md_path_without_ext(md_path_filename(path));
}

bool md_path_is_absolute(md_str_t path) {
    if (md_str_empty(path)) {
        return false;
    }
    const bool sep0 = is_sep(path.ptr[0]);
    if (path.len >= 2 && sep0 && is_sep(path.ptr[1])) {
        // UNC share
        return true;
    }
    if (path.len >= 3 && md_ascii_is_alpha(path.ptr[0]) && path.ptr[1] == ':' && is_sep(path.ptr[2])) {
        // Windows drive
        return true;
    }
    return sep0;
}
