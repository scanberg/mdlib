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

// Offset of the '.' which starts the extension. Only the file name is searched
static bool ext_dot(size_t* loc, md_str_t path) {
    const size_t beg = filename_beg(path);
    for (size_t i = path.len; i > beg; --i) {
        if (path.ptr[i - 1] == '.') {
            *loc = i - 1;
            return true;
        }
    }
    return false;
}

bool md_path_ext(md_str_t* ext, md_str_t path) {
    size_t dot;
    if (!ext_dot(&dot, path)) {
        return false;
    }
    if (ext) {
        ext->ptr = path.ptr + dot + 1;
        ext->len = path.len - dot - 1;
    }
    return true;
}

bool md_path_filename(md_str_t* filename, md_str_t path) {
    if (md_str_empty(path)) {
        return false;
    }
    const size_t beg = filename_beg(path);
    if (filename) {
        filename->ptr = path.ptr + beg;
        filename->len = path.len - beg;
    }
    return true;
}

bool md_path_without_ext(md_str_t* path_without_ext, md_str_t path) {
    size_t dot;
    if (!ext_dot(&dot, path)) {
        return false;
    }
    if (path_without_ext) {
        path_without_ext->ptr = path.ptr;
        path_without_ext->len = dot;
    }
    return true;
}

bool md_path_folder(md_str_t* folder, md_str_t path) {
    const size_t beg = filename_beg(path);
    if (beg == 0) {
        return false;
    }
    if (folder) {
        folder->ptr = path.ptr;
        folder->len = beg;
    }
    return true;
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
