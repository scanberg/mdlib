#ifndef _MD_COMMON_H_
#define _MD_COMMON_H_

#include "md_compiler.h"

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#if MD_COMPILER_MSVC
#ifndef STATIC_ASSERT
#define STATIC_ASSERT static_assert
#endif
#else
#ifndef STATIC_ASSERT
#define STATIC_ASSERT _Static_assert
#endif
#endif
#endif

#if MD_COMPILER_MSVC
#ifndef THREAD_LOCAL
#define THREAD_LOCAL __declspec(thread)
#endif
#else
#ifndef THREAD_LOCAL
#define THREAD_LOCAL __thread
#endif
#endif

#if MD_COMPILER_MSVC
#ifndef ALIGNAS
#define ALIGNAS(x) __declspec(align(x))
#endif
#else
#ifndef ALIGNAS
#define ALIGNAS(x) __attribute__ ((aligned (x)))
#endif
#endif

// Really GCC? REALLY? DO WE REALLY NEED TO INCLUDE stddef.h for this??????
#ifndef NULL
#define NULL 0
#endif

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))
#endif

#ifndef IS_ALIGNED
#define IS_ALIGNED(ptr, alignment) (((uintptr_t)ptr % alignment) == 0)
#endif

#ifndef NEXT_ALIGNED_ADRESS
#define NEXT_ALIGNED_ADRESS(ptr, alignment) ((char*)ptr + alignment - ((uint64_t)ptr & (alignment - 1)))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef ABS
#define ABS(x) ((x) > 0 ? (x) : -(x))
#endif

#ifndef CLAMP
#define CLAMP(val, _min, _max) MIN(MAX(val, _min), _max)
#endif

// Round up division
#ifndef DIV_UP
#define DIV_UP(x, y) ((x + (y-1)) / y)
#endif

// Rounds up x to nearest multiple of y, not the weed-killer
#ifndef ROUND_UP
#define ROUND_UP(x, y) (y * DIV_UP(x,y))
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD(x) ((x)*(3.1415926535897932 / 180.0))
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG(x) ((x)*(180.0 / 3.1415926535897932))
#endif

#ifndef KILOBYTES
#define KILOBYTES(x) (x * 1024LL)
#endif

#ifndef MEGABYTES
#define MEGABYTES(x) (x * 1024LL * 1024LL)
#endif

#ifndef CONCAT
#define CONCAT_INTERNAL(x, y) x##y
#define CONCAT(x, y) CONCAT_INTERNAL(x, y)
#endif

#ifdef  __cplusplus

struct ExitScopeHelp {
    template <typename T>
    struct ExitScope {
        T lambda;
        ExitScope(T lambda) : lambda(lambda) {}
        ~ExitScope() { lambda(); }
        ExitScope& operator=(const ExitScope&) = delete;
    };

    template <typename T>
    ExitScope<T> operator+(T t) {
        return t;
    }
};

#define defer [[maybe_unused]] const auto& CONCAT(defer__, __LINE__) = ExitScopeHelp() + [&]()

#endif //  __cplusplus

#endif