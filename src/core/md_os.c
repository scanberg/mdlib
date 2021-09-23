#include "md_os.h"

#include "md_common.h"
#include "md_platform.h"
#include "md_log.h"
#include "md_allocator.h"

#include <stdio.h>

#if MD_PLATFORM_WINDOWS

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <Windows.h>
#include "Shlwapi.h"
#include <direct.h>

// If we need to explicitly link against some lib
#pragma comment(lib, "User32.lib")
#pragma comment(lib, "Shlwapi.lib")

#define MD_MAX_PATH MAX_PATH

#elif MD_PLATFORM_UNIX
#include <time.h>
#include <unistd.h>
#include <stdlib.h>

#define MD_MAX_PATH PATH_MAX
#endif

#if MD_PLATFORM_OSX
#include <sys/param.h>
#endif

#include <string.h>

#define SZBUF_FROM_STR(buf, str) strncpy(buf, str.ptr, MIN(sizeof(buf)-1, str.len))

char CWD[1024];

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
#endif
    str_t res = {.ptr = CWD, .len = strnlen(CWD, sizeof(CWD))};
    return res;
}

str_t md_os_path_make_canonical(str_t path, struct md_allocator_i* alloc) {
    ASSERT(alloc);

    char sz_buf[MD_MAX_PATH] = "";
    char sz_path[MD_MAX_PATH] = "";

    SZBUF_FROM_STR(sz_path, path);

    str_t canonical_str = {0};
    bool result = false;

#if MD_PLATFORM_WINDOWS
    result = PathCanonicalizeA(sz_buf, sz_path);
#elif MD_PLATFORM_UNIX
    result = realpath(sz_path, sz_buf) == sz_buf;
#endif

    if (result) {
        canonical_str = copy_str((str_t) {.ptr = sz_buf, .len = (int64_t)strnlen(sz_buf, sizeof(sz_buf))}, alloc);
    } else {
        md_print(MD_LOG_TYPE_ERROR, "Failed to create canonical path!");
    }

    return canonical_str;
}

str_t md_os_path_make_relative(str_t from, str_t to, struct md_allocator_i* alloc) {
    char sz_from[MD_MAX_PATH] = "";
    char sz_to[MD_MAX_PATH] = "";
    char sz_buf[MD_MAX_PATH] = "";

    SZBUF_FROM_STR(sz_from, from);
    SZBUF_FROM_STR(sz_to, to);

    str_t relative_str = {0};
    bool result = false;

#if MD_PLATFORM_WINDOWS
    result = PathRelativePathTo(sz_buf, sz_from, FILE_ATTRIBUTE_DIRECTORY, sz_to, FILE_ATTRIBUTE_NORMAL);

#elif MD_PLATFORM_UNIX
    // Make 2 canonical paths
    str_t cfrom = md_os_path_make_canonical(from, default_temp_allocator);
    str_t cto   = md_os_path_make_canonical(to,   default_temp_allocator);

    // Find the common base
    int64_t count = str_count_equal_chars(cfrom, cto);
    result = count > 0;

    if (result) {
        // Count number of folders as N in from and add N times '../'
        str_t rfrom = substr(cfrom, count, -1);
        str_t rto = substr(cto, count, -1);

        int64_t folder_count = str_count_char_occur(rfrom, '/');

        str_t folder_up = make_cstr("../");
        int len = 0;
        for (int64_t i = 0; i < folder_count; ++i) {
            len += snprintf(sz_buf + len, sizeof(sz_buf) - len, "%.*s", (int)folder_up.len, folder_up.ptr);
        }
        len += snprintf(sz_buf + len, sizeof(sz_buf) - len, "%.*s", (int)rto.len, rto.ptr);
    }
#endif
    if (result) {
        relative_str = copy_str((str_t) {.ptr = sz_buf, .len = (int64_t)strnlen(sz_buf, sizeof(sz_buf))}, alloc);
    } else {
        md_print(MD_LOG_TYPE_ERROR, "Failed to extract relative path.");
    }
    return relative_str;
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
#endif
}

double md_os_time_delta_in_ms(timestamp_t t0, timestamp_t t1) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER start, stop, elapsed, frequency;
    QueryPerformanceFrequency(&frequency);
    start.QuadPart = t0;
    stop.QuadPart = t1;
    elapsed.QuadPart = stop.QuadPart - start.QuadPart;
    elapsed.QuadPart *= 1000;
    double ms = ((double)elapsed.QuadPart / (double)frequency.QuadPart);
    return ms;
#elif MD_PLATFORM_UNIX
    return (t1 - t0) * 1.0e-6;
#endif
}

double md_os_time_delta_in_s(timestamp_t t0, timestamp_t t1) {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER start, stop, elapsed, frequency;
    QueryPerformanceFrequency(&frequency);
    start.QuadPart = t0;
    stop.QuadPart = t1;
    elapsed.QuadPart = stop.QuadPart - start.QuadPart;
    elapsed.QuadPart *= 1000000;
    double ms = ((double)elapsed.QuadPart / (double)frequency.QuadPart);
    return ms;
#elif MD_PLATFORM_UNIX
    return (t1 - t0) * 1.0e-9;
#endif
}

void md_os_sleep(int64_t milliseconds) {
#if MD_PLATFORM_WINDOWS
    Sleep((int32_t)milliseconds);
#elif MD_PLATFORM_UNIX
    usleep(milliseconds * 1000);
#endif
}
