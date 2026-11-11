#pragma once

#include <core/md_str.h>

// Lexical path operations: they look at the characters of a path and never touch the file system,
// so they are cheap and behave the same on every platform. Both '/' and '\' are separators, and a
// mix of the two is fine. Results are views into the supplied path.
//
// Operations which have to ask the file system (canonical and relative paths, cwd, exe, validity)
// are in md_os.h.

#ifdef __cplusplus
extern "C" {
#endif

// c:/folder/file.ext -> ext
// The '.' is looked for in the file name only: "/data/run.1/mdout" has no extension.
// False if the file name has no '.'. A leading '.' counts: ".bashrc" -> bashrc
bool md_path_ext(md_str_t* ext, md_str_t path);

// c:/folder/file.ext -> file.ext
// file.ext           -> file.ext   (no folder: the whole path is the file name)
// c:/folder/         -> ""
// False only for an empty path.
bool md_path_filename(md_str_t* filename, md_str_t path);

// c:/folder/file.ext -> c:/folder/file
// The counterpart of md_path_ext: false when md_path_ext is false, and the folder part is never cut.
bool md_path_without_ext(md_str_t* path_without_ext, md_str_t path);

// c:/folder/file.ext -> c:/folder/   (the separator is included)
// False if the path has no separator, e.g. a bare file name.
bool md_path_folder(md_str_t* folder, md_str_t path);

// Checks if a path is absolute, i.e. it does not have to be resolved against a base folder.
// Recognizes both unix roots and windows drive letters / UNC shares on every platform, so
// that a path stored on one platform can be classified on another.
bool md_path_is_absolute(md_str_t path);

#ifdef __cplusplus
}
#endif
