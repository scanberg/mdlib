#pragma once

#include <core/md_compiler.h>

#if defined(__x86_64__) || defined(_M_X64)
#   ifndef __x86_64__
#       define __x86_64__ 1
#   endif
#endif

#if defined(__x86_64__)
#   ifndef __SSE2__
#       define __SSE2__ 1
#   endif
#endif

#ifdef __cplusplus
#   ifndef STATIC_ASSERT
#       define STATIC_ASSERT static_assert
#   endif
#else
#   if MD_COMPILER_MSVC
#       ifndef STATIC_ASSERT
#           define STATIC_ASSERT _Static_assert
#       endif
#   else
#       ifndef STATIC_ASSERT
#           define STATIC_ASSERT _Static_assert
#       endif
#   endif
#endif

#if MD_COMPILER_MSVC
#ifndef THREAD_LOCAL
#define THREAD_LOCAL __declspec(thread)
#endif
#else
#ifndef THREAD_LOCAL
#define THREAD_LOCAL _Thread_local
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

#if MD_COMPILER_MSVC
#ifndef FALLTHROUGH
#define FALLTHROUGH 
#endif
#else
#ifndef FALLTHROUGH
#define FALLTHROUGH __attribute__((fallthrough))
#endif
#endif

#ifndef RESTRICT
#define RESTRICT __restrict
#endif

#ifndef FORCE_INLINE
#if MD_COMPILER_MSVC
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE __attribute__((always_inline)) inline
#endif
#endif

#ifndef TYPEOF
#define TYPEOF(x) __typeof__(x)
#endif

// Really GCC? REALLY? DO WE REALLY NEED TO INCLUDE stddef.h for this??????
#ifndef NULL
#define NULL 0
#endif

#define DEBUG   0
#define RELEASE 0

#ifdef NDEBUG
#   undef  RELEASE
#   define RELEASE 1
#else
#   undef  DEBUG
#   define DEBUG 1
#endif // NDEBUG

#ifndef ARRAY_SIZE
#   ifdef __cplusplus
        template <typename T, unsigned long long N> constexpr unsigned long long array_size_impl(T (&)[N]) { return N; }
#       define ARRAY_SIZE(x) array_size_impl(x)
#   else
#       define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))
#   endif
#endif

#ifndef TYPE_COMPATIBLE
#define TYPE_COMPATIBLE(x, type) _Generic((x), type : 1, default : 0)
#endif

#ifndef IS_POW2
#define IS_POW2(x) (((x) & ((x) - 1)) == 0)
#endif

// This is essentially the same as ROUND_UP, but only works when alignment is a power of two
#ifndef ALIGN_TO
#define ALIGN_TO(x, alignment) (((x) + ((alignment)-1)) & (~(alignment)+1))
#endif

#ifndef IS_ALIGNED
#define IS_ALIGNED(ptr, alignment) (((uintptr_t)(ptr) % (alignment)) == 0)
#endif

#ifndef NEXT_ALIGNED_ADDRESS
#define NEXT_ALIGNED_ADDRESS(ptr, alignment) ((void*)(ALIGN_TO((uintptr_t)(ptr), (alignment))))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef CLAMP
#define CLAMP(val, _min, _max) MIN(MAX(val, _min), _max)
#endif

// Round up division
#ifndef DIV_UP
#define DIV_UP(x, y) (((x) + ((y)-1)) / (y))
#endif

// Rounds up x to nearest multiple of y. No affiliation to the weed-killer
#ifndef ROUND_UP
#define ROUND_UP(x, y) ((y) * DIV_UP(x,y))
#endif

// Rounds down x to nearest multiple of y.
#ifndef ROUND_DOWN
#define ROUND_DOWN(x, y) ((y) * ((x) / (y)))
#endif

#ifndef PI
#define PI 3.1415926535897932
#endif

#ifndef TWO_PI
#define TWO_PI (2.0 * PI)
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD(x) ((x)*(PI / 180.0))
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG(x) ((x)*(180.0 / PI))
#endif

#ifndef KILOBYTES
#define KILOBYTES(x) (x * 1024LL)
#endif

#ifndef MEGABYTES
#define MEGABYTES(x) (x * 1024LL * 1024LL)
#endif

#ifndef GIGABYTES
#define GIGABYTES(x) (x * 1024LL * 1024LL * 1024LL)
#endif

#ifndef CONCAT
#define CONCAT_INTERNAL(x, y) x##y
#define CONCAT(x, y) CONCAT_INTERNAL(x, y)
#endif

#ifndef STRINGIFY
#define STRINGIFY_VAL(s) STRINGIFY(s)
#define STRINGIFY(s) #s
#endif

// Provide declarations for common intrinsic functions which all of the supported compilers expose
#ifndef ASSERT

#   ifdef  __cplusplus
extern "C" void md_assert_impl(const char* file, int line, const char* func_name, const char* expr);
#   else
void md_assert_impl(const char* file, int line, const char* func_name, const char* expr);
#   endif

#   if DEBUG || __FORCE_ASSERTIONS__
#       define ASSERT(__e) \
            ((__e) \
                ? (void)0 \
                : md_assert_impl( \
                    __FILE__, \
                    __LINE__, \
                    __func__, \
                    #__e) \
            )
#   else
#       define ASSERT(__e)
#   endif

#endif

// Define intrinsic functions for memcpy, memset, memmove and memcmp

#if MD_COMPILER_MSVC

#ifdef  __cplusplus
extern "C" {
#endif
void *  __cdecl memcpy(void*, const void*, unsigned long long);
void *  __cdecl memset(void*, int, unsigned long long);
void *  __cdecl memmove(void*, const void*, unsigned long long);
int     __cdecl memcmp(const void*, const void*, unsigned long long);
#ifdef __cplusplus
}
#endif

#pragma intrinsic(memcpy)
#pragma intrinsic(memset)
#pragma intrinsic(memmove)
#pragma intrinsic(memcmp)

#define MEMCPY  memcpy
#define MEMSET  memset
#define MEMMOVE memmove
#define MEMCMP  memcmp

#elif MD_COMPILER_GCC || MD_COMPILER_CLANG
#define MEMCPY  __builtin_memcpy
#define MEMSET  __builtin_memset
#define MEMMOVE __builtin_memmove
#define MEMCMP  __builtin_memcmp

#endif

#ifdef  __cplusplus

#ifndef defer
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
#endif

#endif //  __cplusplus
