#ifndef _MD_COMPILER_H_
#define _MD_COMPILER_H_

#define MD_COMPILER_GCC     0
#define MD_COMPILER_CLANG   0
#define MD_COMPILER_MSVC    0
#define MD_COMPILER_UNKNOWN 0

#define MD_DEBUG            0
#define MD_RELEASE          0

#ifdef NDEBUG
#undef MD_RELEASE
#define MD_RELEASE 1
#else
#undef MD_DEBUG
#define MD_DEBUG 1
#endif // NDEBUG

#if defined(_MSC_VER)
#undef  MD_COMPILER_MSVC
#define MD_COMPILER_MSVC 1
#elif defined(__GNUC__) && !defined(__clang__)
    #undef  MD_COMPILER_GCC
    #define MD_COMPILER_GCC 1
#elif defined(__clang__)
    #undef  MD_COMPILER_CLANG
    #define MD_COMPILER_CLANG 1
#else
    #undef  MD_COMPILER_UNKNOWN
    #define MD_COMPILER_UNKNOWN 1
#endif

#endif