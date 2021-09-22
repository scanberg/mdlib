#include "md_os.h"

#include "md_common.h"
#include "md_platform.h"

#if MD_PLATFORM_WINDOWS
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include <Windows.h>
#include <direct.h>

// If we need to explicitly link against some lib
#pragma comment(lib, "User32.lib")

#elif MD_PLATFORM_UNIX
#include <time.h>
#include <unistd.h>
#endif

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

timestamp_t md_os_time_current() {
#if MD_PLATFORM_WINDOWS
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return t.QuadPart;
#elif MD_PLATFORM_LINUX || MD_PLATFORM_OSX
    timespec t;
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
    return (ts1 - ts0) * 1.0e-6;
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
    return (ts1 - ts0) * 1.0e-9;
#endif
}