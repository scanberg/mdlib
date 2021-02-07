#ifndef _MD_PLATFORM_H_
#define _MD_PLATFORM_H_

#define MD_PLATFORM_WINDOWS 0
#define MD_PLATFORM_LINUX   0
#define MD_PLATFORM_OSX     0

#if defined(_WIN32) || defined(_WIN64)
    #undef  MD_PLATFORM_WINDOWS
    #define MD_PLATFORM_WINDOWS 1
#elif __APPLE__
    #include "TargetConditionals.h"
    #ifdef  TARGET_OS_MAC
        #undef  MD_PLATFORM_OSX
        #define MD_PLATFORM_OSX 1
    #endif
#elif defined(__linux__)
    #undef  MD_PLATFORM_LINUX
    #define MD_PLATFORM_LINUX 1
#else
    #error "Platform Unknown"
#endif

#endif