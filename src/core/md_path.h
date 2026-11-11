#pragma once

#include <core/md_str.h>

// Paths. Both '/' and '\' are separators, and a mix of the two is fine.
//
// Two kinds of operation: the LEXICAL ones look at the characters of a path only, so they are cheap,
// thread safe, and give the same answer on every platform for paths from any platform. Those under
// FILE SYSTEM ask the disk or the process (the cwd, what exists), so their answers depend on when and
// where they are asked; they are implemented in md_os.c.

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// ### LEXICAL ###

// The parts of a path, as views into it. Each is the part asked for, and empty exactly when the path has no
// such part: folder + filename is the path, and so is without_ext + '.' + ext whenever there is an extension.

// 'c:/folder/file.ext' -> 'file.ext'
// 'file.ext'           -> 'file.ext'   (no folder: the whole path is the file name)
// 'c:/folder/'         -> ''
md_str_t md_path_filename(md_str_t path);

// 'c:/folder/file.ext' -> 'c:/folder/'  (with its separator)
// 'file.ext'           -> ''            (a bare file name, which is in the cwd)
md_str_t md_path_folder(md_str_t path);

// 'c:/folder/file.ext' -> 'ext'
// 'file.tar.gz'        -> 'gz'
// Looked for in the file name only, and not at its first character: '/data/run.1/mdout' has none, and
// neither has '.bashrc' (a leading '.' marks a hidden file). 'file.' has an empty one, the same as none.
md_str_t md_path_ext(md_str_t path);

// 'c:/folder/file.ext' -> 'c:/folder/file'
// All of path when it has no extension.
md_str_t md_path_without_ext(md_str_t path);

// The file name without its extension
// 'c:/folder/file.tar.gz' -> 'file.tar'
// '.bashrc'               -> '.bashrc'
md_str_t md_path_stem(md_str_t path);

// Checks if a path is absolute, i.e. it does not have to be resolved against a base folder.
// Recognizes both unix roots and windows drive letters / UNC shares on every platform, so
// that a path stored on one platform can be classified on another.
bool md_path_is_absolute(md_str_t path);

// ### FILE SYSTEM ###
// These ask the disk or the process, and are implemented in md_os.c. md_path_set_cwd changes the process's
// current directory, which every thread shares.

// Gets the current working directory
// Writes into a buffer (up to buf_cap) and returns the number of characters written
size_t md_path_write_cwd(char* buf, size_t buf_cap);
bool  md_path_set_cwd(md_str_t path);

// Gets the path to the executable
// Writes into a buffer (up to buf_cap) and returns the number of characters written
size_t md_path_write_exe(char* buf, size_t buf_cap);

// Gets the path to the user directory
// Writes into a buffer (up to buf_cap) and returns the number of characters written
// On windows this would correspond to the user's home directory e.g. C:\Users\username
// On linux this would correspond to the user's home directory e.g. /home/username
// On mac this would correspond to the user's home directory e.g. /Users/username
size_t md_path_write_user_dir(char* buf, size_t buf_cap);

// Writes a canonical path to buf.
size_t md_path_write_canonical(char* buf, size_t buf_cap, md_str_t path);
md_str_t md_path_make_canonical(md_str_t path, struct md_allocator_i* alloc);

// Writes a relative path from path_from to path_to to buf.
size_t md_path_write_relative(char* buf, size_t buf_cap, md_str_t path_from, md_str_t path_to);
md_str_t md_path_make_relative(md_str_t path_from, md_str_t path_to, struct md_allocator_i* alloc);

// Checks if a path to a file or directory is valid: i.e. it points to an actual file / dir on disk.
bool md_path_is_valid(md_str_t path);
bool md_path_is_directory(md_str_t path);

#ifdef __cplusplus
}
#endif
