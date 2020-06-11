#ifndef _MOLD_PLATFORM_H_
#define _MOLD_PLATFORM_H_

#define MOLD_PLATFORM_WINDOWS 0
#define MOLD_PLATFORM_LINUX   0
#define MOLD_PLATFORM_OSX     0

#if defined(_WIN32) || defined(_WIN64)
    #undef  MOLD_PLATFORM_WINDOWS
    #define MOLD_PLATFORM_WINDOWS 1
#elif __APPLE__
    #include "TargetConditionals.h"
    #ifdef  TARGET_OS_MAC
        #undef  MOLD_PLATFORM_OSX
        #define MOLD_PLATFORM_OSX 1
    #endif
#elif defined __linux__
    #undef  MOLD_PLATFORM_LINUX
    #define MOLD_PLATFORM_LINUX 1
#else
    #error "Platform Unknown"
#endif

#endif