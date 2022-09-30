#include "md_os.h"

#include "md_common.h"
#include "md_platform.h"
#include "md_log.h"
#include "md_allocator.h"

#include <stdio.h>

#define MD_MAX_PATH 4096

#include <time.h>

#if MD_PLATFORM_WINDOWS

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <Windows.h>
#include <Shlwapi.h>
#include <direct.h>

// If we need to explicitly link against some lib
#pragma comment(lib, "User32.lib")
#pragma comment(lib, "Shlwapi.lib")

#elif MD_PLATFORM_UNIX

#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/mman.h>

#endif

#if MD_PLATFORM_OSX
#include <sys/param.h>
#endif

#include <string.h>

#define SZBUF_FROM_STR(buf, str) strncpy(buf, str.ptr, MIN(sizeof(buf)-1, str.len))

static char CWD[4096];

#if MD_PLATFORM_WINDOWS
// https://docs.microsoft.com/en-us/windows/win32/debug/retrieving-the-last-error-code
void print_windows_error() {
    LPVOID msg_buf;
    DWORD err_code = GetLastError();

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        err_code,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR)&msg_buf,
        0, NULL);

    md_printf(MD_LOG_TYPE_ERROR, "Error code (%d): %s", err_code, msg_buf);
    LocalFree(msg_buf);
}
#endif

int64_t md_os_physical_ram_in_bytes() {
    ASSERT(false);
    return 0;
}

str_t md_os_current_working_directory() {
#if MD_PLATFORM_WINDOWS
    char* val = _getcwd(CWD, sizeof(CWD));
    ASSERT(val != 0);
#elif MD_PLATFORM_UNIX
    getcwd(CWD, sizeof(CWD));
#else
    ASSERT(false);
#endif
    str_t res = { .ptr = CWD, .len = (int64_t)strnlen(CWD, sizeof(CWD)) };
    return res;
}

static inline str_t internal_fullpath(str_t path, md_allocator_i* alloc) {
    path = str_copy(path, default_temp_allocator); // Zero terminate

    char sz_buf[MD_MAX_PATH] = "";
#if MD_PLATFORM_WINDOWS
    int64_t len = (int64_t)GetFullPathName(path.ptr, sizeof(sz_buf), sz_buf, NULL);
    if (len == 0) {
        print_windows_error();
    }
#elif MD_PLATFORM_UNIX
    int64_t len = 0;
    if (realpath(path.ptr, sz_buf) != NULL) {
        len = (int64_t)strnlen(sz_buf, sizeof(sz_buf));
    }
    else {
        switch (errno) {
        case EACCES:        md_printf(MD_LOG_TYPE_ERROR, "Read or search permission was denied for a component of the path prefix."); break;
        case EINVAL:        md_printf(MD_LOG_TYPE_ERROR, "path is NULL or resolved_path is NULL."); break;
        case EIO:           md_printf(MD_LOG_TYPE_ERROR, "An I/O error occurred while reading from the filesystem."); break;
        case ELOOP:         md_printf(MD_LOG_TYPE_ERROR, "Too many symbolic links were encountered in translating the pathname."); break;
        case ENAMETOOLONG:  md_printf(MD_LOG_TYPE_ERROR, "A component of a pathname exceeded NAME_MAX characters, or an entire pathname exceeded PATH_MAX characters."); break;
        case ENOENT:        md_printf(MD_LOG_TYPE_ERROR, "The named file does not exist."); break;
        case ENOMEM:        md_printf(MD_LOG_TYPE_ERROR, "Out of memory."); break;
        case ENOTDIR:       md_printf(MD_LOG_TYPE_ERROR, "A component of the path prefix is not a directory."); break;
        default:            md_printf(MD_LOG_TYPE_ERROR, "Undefined error"); break;
        }
    }
#else
    ASSERT(false);
#endif
    str_t result = { 0 };
    if (len > 0) {
        result = str_copy((str_t) { sz_buf, len }, alloc);
    }
    return result;
}

str_t md_os_path_make_canonical(str_t path, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    path = internal_fullpath(path, alloc);

    if (path.len > 0) {
#if MD_PLATFORM_WINDOWS
        convert_backslashes((char*)path.ptr, path.len);
#endif
    }
    else {
        md_printf(MD_LOG_TYPE_ERROR, "Failed to create canonical path!");
    }

    return path;
}

str_t md_os_path_make_relative(str_t from, str_t to, struct md_allocator_i* alloc) {
    char sz_buf[MD_MAX_PATH] = "";

    str_t relative_str = { 0 };
    bool result = false;

    // Make 2 canonical paths
    from = internal_fullpath(from, default_temp_allocator);
    to = internal_fullpath(to, default_temp_allocator);

#if MD_PLATFORM_WINDOWS
    result = PathRelativePathTo(sz_buf, from.ptr, FILE_ATTRIBUTE_NORMAL, to.ptr, FILE_ATTRIBUTE_NORMAL);
    convert_backslashes(sz_buf, sizeof(sz_buf));
#elif MD_PLATFORM_UNIX

    // Find the common base
    int64_t count = str_count_equal_chars(from, to);
    result = count > 0;

    if (result) {
        // Count number of folders as N in from and add N times '../'
        str_t rfrom = substr(from, count, -1);
        str_t rto = substr(to, count, -1);

        const int64_t folder_count = str_count_char_occur(rfrom, '/');

        str_t folder_up = MAKE_STR("../");
        int len = 0;
        for (int64_t i = 0; i < folder_count; ++i) {
            len += snprintf(sz_buf + len, sizeof(sz_buf) - len, "%.*s", (int)folder_up.len, folder_up.ptr);
        }
        len += snprintf(sz_buf + len, sizeof(sz_buf) - len, "%.*s", (int)rto.len, rto.ptr);
    }
#else
    ASSERT(false);
#endif
    if (result) {
        relative_str = str_copy((str_t) { .ptr = sz_buf, .len = (int64_t)strnlen(sz_buf, sizeof(sz_buf)) }, alloc);
    }
    else {
        md_printf(MD_LOG_TYPE_ERROR, "Failed to extract relative path.");
    }
    return relative_str;
}

bool md_os_path_exists(str_t path) {
#if MD_PLATFORM_WINDOWS
    path = str_copy(path, default_temp_allocator);
    bool result = PathFileExists(path.ptr);
#elif MD_PLATFORM_UNIX
    bool result = (access(path.ptr, F_OK) == 0);
#else
    ASSERT(false);
#endif
    return result;
}

bool md_os_path_is_directory(str_t path) {
#if MD_PLATFORM_WINDOWS
    path = str_copy(path, default_temp_allocator);
    bool result = PathIsDirectory(path.ptr);
#elif MD_PLATFORM_UNIX
    bool result = false;
    struct stat s;
    if (stat(path.ptr, &s) == 0) {
        if (s.st_mode & S_IFDIR) {
            result = true;
        }
    }
#else
    ASSERT(false);
#endif
    return result;
}

timestamp_t md_os_time_current() {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return t.QuadPart;
#elif MD_PLATFORM_UNIX
    struct timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return t.tv_sec * 1000000000 + t.tv_nsec;
#else
    ASSERT(false);
#endif
}

timestamp_t md_os_time_from_milliseconds(int64_t ms) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return (ms * 3000) / frequency.QuadPart;
#elif MD_PLATFORM_UNIX
    return ms * 1.0e-6;
#else
    ASSERT(false);
#endif
}

double md_os_time_as_milliseconds(timestamp_t t) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return (double)(t * 1000) / (double)frequency.QuadPart;
#elif MD_PLATFORM_UNIX
    return t * 1.0e-6;
#else
    ASSERT(false);
#endif
}

double md_os_time_as_seconds(timestamp_t t) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    return (double)(t * 1000) / (double)frequency.QuadPart;
#elif MD_PLATFORM_UNIX
    return t * 1.0e-9;
#else
    ASSERT(false);
#endif
}

void md_os_sleep(int64_t milliseconds) {
#if MD_PLATFORM_WINDOWS
    Sleep((int32_t)milliseconds);
#elif MD_PLATFORM_UNIX
    usleep(milliseconds * 1000);
#else
    ASSERT(false);
#endif
}

// MEM
// Linux equivalent is partly taken from here
// https://forums.pcsx2.net/Thread-blog-VirtualAlloc-on-Linux

uint64_t md_os_page_size(void) {
#if MD_PLATFORM_WINDOWS
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwPageSize;
#elif MD_PLATFORM_UNIX
    return sysconf(_SC_PAGE_SIZE);
#else
    ASSERT(false);
#endif
}

void* md_os_reserve(uint64_t size) {
    uint64_t gb_snapped_size = ALIGN_TO(size, GIGABYTES(1));
#if MD_PLATFORM_WINDOWS
    return VirtualAlloc(0, gb_snapped_size, MEM_RESERVE, PAGE_NOACCESS);
#elif MD_PLATFORM_UNIX
    void* result = mmap(0, gb_snapped_size, PROT_NONE, MAP_PRIVATE | MAP_ANON, -1, 0);
    ASSERT(result != MAP_FAILED);
    return result;
#else
    ASSERT(false);
#endif
}

void md_os_release(void* ptr) {
#if MD_PLATFORM_WINDOWS
    VirtualFree(ptr, 0, MEM_RELEASE);
#elif MD_PLATFORM_UNIX
    munmap(ptr, -1);
#else
    ASSERT(false);
#endif
}

void md_os_commit(void* ptr, uint64_t size) {
    uint64_t page_snapped_size = ALIGN_TO(size, md_os_page_size());
#if MD_PLATFORM_WINDOWS
    VirtualAlloc(ptr, page_snapped_size, MEM_COMMIT, PAGE_READWRITE);
#elif MD_PLATFORM_UNIX
    int result = mprotect(ptr, page_snapped_size, PROT_READ | PROT_WRITE);
    ASSERT(result == 0);
#else
    ASSERT(false);
#endif
}

void md_os_decommit(void* ptr, uint64_t size) {
#if MD_PLATFORM_WINDOWS
    VirtualFree(ptr, size, MEM_DECOMMIT);
#elif MD_PLATFORM_UNIX
    int res;
    uint64_t page_snapped_size = ALIGN_TO(size, md_os_page_size());
    res = mprotect(ptr, page_snapped_size, PROT_NONE);
    ASSERT(res == 0);
    res = madvise(ptr, page_snapped_size, MADV_DONTNEED);
    ASSERT(res == 0);
#else
    ASSERT(false);
#endif
}