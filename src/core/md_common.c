#include "md_common.h"

#include <stdio.h>

#if MD_COMPILER_MSVC
#include <intrin.h>
#elif MC_COMPILER_GCC || MD_COMPILER_CLANG
#else
#include <stdlib.h>
#endif

void md_assert_impl(const char* file, int line, const char* func_name, const char* expr) {
    fprintf(stderr, "Assertion \"%s\" failed at line %i in %s:%s\n", expr, line, file, func_name);
#if MD_COMPILER_MSVC
    __debugbreak();
#elif MD_COMPILER_GCC || MD_COMPILER_CLANG
    __builtin_trap();
#else
    abort();
#endif
}
