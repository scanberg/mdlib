#include "md_file.h"

#include "md_platform.h"
#include "md_compiler.h"
#include "md_common.h"
#include "md_log.h"


#if MD_PLATFORM_WINDOWS
#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>

#pragma warning(disable:4063) // combined flags not valid for switch of enum
#endif

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#ifndef MD_MAX_PATH
#define MD_MAX_PATH 4096
#endif

struct md_file_o {
    FILE handle;
};

md_file_o* md_file_open(str_t filename, int flags) {
    if (!filename.ptr || !filename.len) return NULL;

#if MD_PLATFORM_WINDOWS
    wchar_t w_file[MD_MAX_PATH] = {0};
    const int w_file_len = MultiByteToWideChar(CP_UTF8, 0, filename.ptr, (int)filename.len, w_file, ARRAY_SIZE(w_file));
    if (w_file_len >= ARRAY_SIZE(w_file)) {
        md_print(MD_LOG_TYPE_ERROR, "File path exceeds stupid limit!");
        return NULL;
    }
    w_file[w_file_len] = L'\0';

    const wchar_t* w_mode = 0;
    switch (flags) {
    case MD_FILE_READ:
        w_mode = L"r";
        break;
    case (MD_FILE_READ | MD_FILE_BINARY):
        w_mode = L"rb";
        break;
    case MD_FILE_WRITE:
        w_mode = L"w";
        break;
    case (MD_FILE_WRITE | MD_FILE_BINARY):
        w_mode = L"wb";
        break;
    case MD_FILE_APPEND:
        w_mode = L"a";
        break;
    case (MD_FILE_APPEND | MD_FILE_BINARY):
        w_mode = L"ab";
        break;
    default:
        md_print(MD_LOG_TYPE_ERROR, "Invalid combination of file access flags!");
        return NULL;
    }
    
    return (md_file_o*) _wfopen(w_file, w_mode);
#else
    const char* mode = 0;
    switch (flags) {
    case MD_FILE_READ:
        mode = "r";
        break;
    case (MD_FILE_READ | MD_FILE_BINARY):
        mode = "rb";
        break;
    case MD_FILE_WRITE:
        mode = "w";
        break;
    case (MD_FILE_WRITE | MD_FILE_BINARY):
        mode = "wb";
        break;
    case MD_FILE_APPEND:
        mode = "a";
        break;
    case (MD_FILE_APPEND | MD_FILE_BINARY):
        mode = "ab";
        break;
    default:
        md_print(MD_LOG_TYPE_ERROR, "Invalid combination of file access flags!");
        return NULL;
    }

    char file[MD_MAX_PATH] = "";
    strncpy(file, filename.ptr, ARRAY_SIZE(file)-1);
    return (md_file_o*)fopen(file, mode);
#endif
}

void md_file_close(md_file_o* file) {
    fclose((FILE*)file);
}

int64_t md_file_tell(md_file_o* file) {
#if MD_PLATFORM_WINDOWS
    return _ftelli64((FILE*)file);
#else
    return ftello((FILE*)file);
#endif
}

bool md_file_seek(md_file_o* file, int64_t offset, md_file_seek_origin_t origin) {
    ASSERT(file);

    int o = 0;
    switch(origin) {
    case MD_FILE_BEG:
        o = SEEK_SET;
        break;
    case MD_FILE_CUR:
        o = SEEK_CUR;
        break;
    case MD_FILE_END:
        o = SEEK_END;
        break;
    default:
        md_print(MD_LOG_TYPE_ERROR, "Invalid seek origin value!");
        return false;
    }

#if MD_PLATFORM_WINDOWS
    return _fseeki64((FILE*)file, offset, o) == 0;
#else
    return fseeko((FILE*)file, offset, o) == 0;
#endif
}

int64_t md_file_size(md_file_o* file) {
    ASSERT(file);
    int64_t cur = md_file_tell(file);
    md_file_seek(file, 0, SEEK_END);
    int64_t end = md_file_tell(file);
    md_file_seek(file, cur, SEEK_SET);
    return end;
}

// Returns the number of successfully written/read bytes
int64_t md_file_read(md_file_o* file, void* ptr, int64_t num_bytes) {
    ASSERT(file);
    return fread(ptr, 1, num_bytes, (FILE*)file);
}

int64_t md_file_write(md_file_o* file, const void* ptr, int64_t num_bytes) {
    ASSERT(file);
    return fwrite(ptr, 1, num_bytes, (FILE*)file);
}

int64_t md_file_printf(md_file_o* file, const char* format, ...) {
    ASSERT(file);
    va_list args;
    va_start (args, format);
    int res = vfprintf((FILE*)file, format, args);
    va_end (args);
    return res;
}
