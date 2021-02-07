#include "file.h"
#include "platform.h"

#if MD_PLATFORM_WINDOWS
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <stdio.h>

struct md_file {
    void* _ptr;
};

md_file* md_file_open(const char* filename, uint32_t filename_len, const char* mode) {
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

    return (md_file*)_wfopen(w_file, w_mode);
#else
    return (md_file*)fopen(file, mode);
#endif
}

void md_file_close(md_file* file) {
    fclose((FILE*)file);
}

int64_t md_file_tell(md_file* file) {
#if MD_PLATFORM_WINDOWS
    return _ftelli64((FILE*)file);
#else
    return ftello((FILE*)file);
#endif
}

bool md_file_seek(md_file* file, int64_t offset, int origin) {
#if MD_PLATFORM_WINDOWS
    return _fseeki64((FILE*)file, offset, origin) == 0;
#else
    return fseeko((FILE*)file, offset, origin) == 0;
#endif
}

uint64_t md_file_size(md_file* file) {
    if (!file) return 0;
    int64_t cur = md_file_tell(file);
    md_file_seek(file, 0, SEEK_END);
    int64_t end = md_file_tell(file);
    md_file_seek(file, cur, SEEK_SET);
    return (uint32_t)end;
}

// Returns the number of successfully written/read bytes
uint64_t md_file_read(md_file* file, void* ptr, uint64_t num_bytes) {
    //size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );
    return fread(ptr, 1, num_bytes, (FILE*)file);
}

uint64_t md_file_write(md_file* file, const void* ptr, uint64_t num_bytes) {
    return fwrite(ptr, 1, num_bytes, (FILE*)file);
}