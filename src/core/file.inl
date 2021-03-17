#include "platform.h"
#include "compiler.h"

#if MD_COMPILER_MSVC
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif

#if MD_PLATFORM_WINDOWS
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

// This should have been file.h and file.c, but that would force us to include <stdio.h> since we don't have the definition of FILE
// And cannot forwarddeclare it since it is a define and the implementation differs for compilers. So to not violate the rules
// Of only having stdint.h and stdbool.h as includes in our header, this has been turned into a .inl ...
// It might be worthy of an exception since it is just one include... BUT HEY!!!

static FILE* md_file_open(const char* filename, uint64_t filename_len, const char* mode) {
#if MD_PLATFORM_WINDOWS
    const int file_len = (int)filename_len;
    const int mode_len = (int)strlen(mode);

    if (file_len == 0) return NULL;
    if (mode_len == 0) return NULL;

    wchar_t w_file[MAX_PATH];
    wchar_t w_mode[MAX_PATH];

    const int w_file_len = MultiByteToWideChar(CP_UTF8, 0, filename, file_len, w_file, MAX_PATH);
    if (w_file_len >= MAX_PATH) return NULL;
    w_file[w_file_len] = L'\0';

    const int w_mode_len = MultiByteToWideChar(CP_UTF8, 0, mode, mode_len, w_mode, MAX_PATH);
    if (w_mode_len >= MAX_PATH) return NULL;
    w_mode[w_mode_len] = L'\0';
    
    return _wfopen(w_file, w_mode);
#else
    return fopen(filename, mode);
#endif
}

static inline void md_file_close(FILE* file) {
    fclose(file);
}

static inline int64_t md_file_tell(FILE* file) {
#if MD_PLATFORM_WINDOWS
    return _ftelli64(file);
#else
    return ftello(file);
#endif
}

static inline bool md_file_seek(FILE* file, int64_t offset, int origin) {
#if MD_PLATFORM_WINDOWS
    return _fseeki64((FILE*)file, offset, origin) == 0;
#else
    return fseeko((FILE*)file, offset, origin) == 0;
#endif
}

static inline uint64_t md_file_size(FILE* file) {
    if (!file) return 0;
    int64_t cur = md_file_tell(file);
    md_file_seek(file, 0, SEEK_END);
    int64_t end = md_file_tell(file);
    md_file_seek(file, cur, SEEK_SET);
    return (uint32_t)end;
}

// Returns the number of successfully written/read bytes
static inline uint64_t md_file_read(FILE* file, void* ptr, uint64_t num_bytes) {
    //size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );
    return fread(ptr, 1, num_bytes, file);
}

static inline uint64_t md_file_write(FILE* file, const void* ptr, uint64_t num_bytes) {
    return fwrite(ptr, 1, num_bytes, file);
}
