#ifndef _MOLD_COMPILER_H_
#define _MOLD_COMPILER_H_

#define MOLD_COMPILER_GCC     0
#define MOLD_COMPILER_CLANG   0
#define MOLD_COMPILER_MSVC    0
#define MOLD_COMPILER_UNKNOWN 0

#if defined(_MSC_VER)
#undef  MOLD_COMPILER_MSVC
#define MOLD_COMPILER_MSVC 1
#elif defined(__GNUC__) && !defined(__clang__)
    #undef  MOLD_COMPILER_GCC
    #define MOLD_COMPILER_GCC 1
#elif defined(__clang__)
    #undef  MOLD_COMPILER_CLANG
    #define MOLD_COMPILER_CLANG 1
#else
    #undef  MOLD_COMPILER_UNKNOWN
    #define MOLD_COMPILER_UNKNOWN 1
#endif

#endif