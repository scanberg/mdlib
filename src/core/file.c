#include "file.h"
#include "platform.h"

#if MOLD_PLATFORM_WINDOWS
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <stdio.h>

struct mold_file {
    void* _placeholder;
};

mold_file* mold_file_open(const char* filename, const char* mode) {
#if MOLD_PLATFORM_WINDOWS
    const int file_len = (int)strlen(filename);
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

    return (mold_file*)_wfopen(w_file, w_mode);
#else
    return (mold_file*)fopen(file, mode);
#endif
}

void mold_file_close(mold_file* file) {
    fclose((FILE*)file);
}

int64_t mold_file_tell(mold_file* file) {
#if MOLD_PLATFORM_WINDOWS
    return _ftelli64((FILE*)file);
#else
    return ftello((FILE*)file);
#endif
}

bool mold_file_seek(mold_file* file, int64_t offset, int origin) {
#if MOLD_PLATFORM_WINDOWS
    return _fseeki64((FILE*)file, offset, origin);
#else
    return fseeko((FILE*)file, offset, origin);
#endif
}

// Returns the number of successfully written/read bytes
uint64_t mold_file_read(mold_file* file, void* ptr, uint64_t num_bytes) {
    //size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );
    return fread(ptr, 1, num_bytes, (FILE*)file);
}

uint64_t mold_file_write(mold_file* file, const void* ptr, uint64_t num_bytes) {
    return fwrite(ptr, 1, num_bytes, (FILE*)file);
}