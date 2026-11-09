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

#ifdef  __cplusplus
#define FALLTHROUGH  [[fallthrough]]
#else
#if MD_COMPILER_MSVC
#ifndef FALLTHROUGH
#define FALLTHROUGH
#endif
#else
#ifndef FALLTHROUGH
#define FALLTHROUGH __attribute__((fallthrough))
#endif
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

// State an invariant the optimiser cannot derive on its own. Unlike ASSERT this is NOT a check:
// nothing is tested at runtime in any build, and a false expression is undefined behaviour. Use it
// only where the invariant is guaranteed by construction and the compiler's inability to see it
// costs something concrete - a spurious diagnostic, or a branch that can never be taken.
#ifndef MD_ASSUME
#   if MD_COMPILER_MSVC
#       define MD_ASSUME(expr) __assume(expr)
#   elif MD_COMPILER_GCC || MD_COMPILER_CLANG
#       define MD_ASSUME(expr) ((expr) ? (void)0 : __builtin_unreachable())
#   else
#       define MD_ASSUME(expr) ((void)0)
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

// ## Intrinsics
//
// THE RULE: what resolves without an include on every supported compiler is here; what needs a header to get
// is not (the SIMD types and their operations, md_simd.h). GCC and Clang - clang-cl too, which is MD_COMPILER_CLANG -
// have builtins for all of it. MSVC gets each intrinsic from a declaration of our own, made an intrinsic with #pragma
// intrinsic, exactly as <intrin.h> and winnt.h declare it, so that a translation unit which includes those as well
// sees the same declarations. md_intrinsics.h is what remains of the old split, kept for its includers.

// The types <stdint.h> calls uint32_t and uint64_t, without it. The same types, not merely the same sizes, so that a
// uint64_t* passes where an MD_U64* is asked for: GCC and Clang predefine them as stdint.h's own, and MSVC's
// unsigned __int64 is its unsigned long long.
#if MD_COMPILER_MSVC
#   define MD_U32 unsigned int
#   define MD_U64 unsigned __int64
#else
#   define MD_U32 __UINT32_TYPE__
#   define MD_U64 __UINT64_TYPE__
#endif

#if MD_COMPILER_MSVC

#if !defined(_M_X64) && !defined(_M_ARM64)
#   error "MSVC: x64 and ARM64 are the supported targets"
#endif

#ifdef __cplusplus
extern "C" {
#endif

    void *  __cdecl             memcpy(void*, const void*, unsigned long long);
    void *  __cdecl             memset(void*, int, unsigned long long);
    void *  __cdecl             memmove(void*, const void*, unsigned long long);
    int     __cdecl             memcmp(const void*, const void*, unsigned long long);
    unsigned short __cdecl      _byteswap_ushort (unsigned short);
    unsigned long __cdecl       _byteswap_ulong  (unsigned long);
    unsigned __int64 __cdecl    _byteswap_uint64 (unsigned __int64);
    unsigned char __cdecl       _BitScanForward  (unsigned long* index, unsigned long mask);
    unsigned char __cdecl       _BitScanReverse  (unsigned long* index, unsigned long mask);
    unsigned char __cdecl       _BitScanForward64(unsigned long* index, unsigned __int64 mask);
    unsigned char __cdecl       _BitScanReverse64(unsigned long* index, unsigned __int64 mask);
    __int64 __cdecl             _InterlockedExchangeAdd64(__int64 volatile* addend, __int64 value);
    __int64 __cdecl             _InterlockedCompareExchange64(__int64 volatile* destination, __int64 exchange, __int64 comparand);
    unsigned __int64 __cdecl    __umulh(unsigned __int64 a, unsigned __int64 b);
    void __cdecl                __debugbreak(void);
#if defined(_M_X64)
    unsigned int __cdecl        __popcnt  (unsigned int);
    unsigned __int64 __cdecl    __popcnt64(unsigned __int64);
#endif

#ifdef __cplusplus
}
#endif

#pragma intrinsic(memcpy)
#pragma intrinsic(memset)
#pragma intrinsic(memmove)
#pragma intrinsic(memcmp)
#pragma intrinsic(_byteswap_ushort)
#pragma intrinsic(_byteswap_ulong)
#pragma intrinsic(_byteswap_uint64)
#pragma intrinsic(_BitScanForward)
#pragma intrinsic(_BitScanReverse)
#pragma intrinsic(_BitScanForward64)
#pragma intrinsic(_BitScanReverse64)
#pragma intrinsic(_InterlockedExchangeAdd64)
#pragma intrinsic(_InterlockedCompareExchange64)
#pragma intrinsic(__umulh)
#pragma intrinsic(__debugbreak)
#if defined(_M_X64)
#pragma intrinsic(__popcnt)
#pragma intrinsic(__popcnt64)
#endif

#define MEMCPY  memcpy
#define MEMSET  memset
#define MEMMOVE memmove
#define MEMCMP  memcmp
#define BSWAP16 _byteswap_ushort
#define BSWAP32 _byteswap_ulong
#define BSWAP64 _byteswap_uint64
#define DEBUG_BREAK() __debugbreak()

#elif MD_COMPILER_GCC || MD_COMPILER_CLANG

#define MEMCPY   __builtin_memcpy
#define MEMSET   __builtin_memset
#define MEMMOVE  __builtin_memmove
#define MEMCMP   __builtin_memcmp
#define BSWAP16  __builtin_bswap16
#define BSWAP32  __builtin_bswap32
#define BSWAP64  __builtin_bswap64
#if defined(__x86_64__) || defined(__i386__)
#   define DEBUG_BREAK() __asm__ volatile("int3")
#else
#   define DEBUG_BREAK() __builtin_trap()
#endif

#endif

#ifndef NOINLINE
#   if MD_COMPILER_MSVC
#       define NOINLINE __declspec(noinline)
#   else
#       define NOINLINE __attribute__((noinline))
#   endif
#endif

// ### Bits
// The same results on every compiler and every CPU, for every input: clz and ctz of 0 are the width, popcnt, bsf and
// bsr are defined throughout. (MSVC's __lzcnt is not used: on a CPU without LZCNT - Sandy and Ivy Bridge, AVX or not -
// it executes as BSR and gives the index of the top bit instead of the count. _BitScan* is the same everywhere.)

#if MD_COMPILER_GCC || MD_COMPILER_CLANG

static inline unsigned int popcnt32(MD_U32 x) { return (unsigned int)__builtin_popcount(x); }
static inline unsigned int popcnt64(MD_U64 x) { return (unsigned int)__builtin_popcountll(x); }

// Count of leading (most significant) zero bits; 32 / 64 for 0
static inline unsigned int clz32(MD_U32 x) { return x ? (unsigned int)__builtin_clz(x)   : 32u; }
static inline unsigned int clz64(MD_U64 x) { return x ? (unsigned int)__builtin_clzll(x) : 64u; }

// Count of trailing (least significant) zero bits; 32 / 64 for 0. With BMI enabled this is one TZCNT, without it a BSF
// and a CMOV: the compilers know the idiom.
static inline unsigned int ctz32(MD_U32 x) { return x ? (unsigned int)__builtin_ctz(x)   : 32u; }
static inline unsigned int ctz64(MD_U64 x) { return x ? (unsigned int)__builtin_ctzll(x) : 64u; }

#elif MD_COMPILER_MSVC

#if defined(_M_X64)
static inline unsigned int popcnt32(MD_U32 x) { return __popcnt(x); }
static inline unsigned int popcnt64(MD_U64 x) { return (unsigned int)__popcnt64(x); }
#else
// ARM64: a SWAR count, which needs no declaration of _CountOneBits to be trusted
static inline unsigned int popcnt64(MD_U64 x) {
    x = x - ((x >> 1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Full;
    return (unsigned int)((x * 0x0101010101010101ull) >> 56);
}
static inline unsigned int popcnt32(MD_U32 x) { return popcnt64(x); }
#endif

static inline unsigned int clz32(MD_U32 x) { unsigned long i; return _BitScanReverse  (&i, x) ? 31u - (unsigned int)i : 32u; }
static inline unsigned int clz64(MD_U64 x) { unsigned long i; return _BitScanReverse64(&i, x) ? 63u - (unsigned int)i : 64u; }
static inline unsigned int ctz32(MD_U32 x) { unsigned long i; return _BitScanForward  (&i, x) ? (unsigned int)i : 32u; }
static inline unsigned int ctz64(MD_U64 x) { unsigned long i; return _BitScanForward64(&i, x) ? (unsigned int)i : 64u; }

#endif

// The 1-based index of the lowest (bsf) and highest (bsr) set bit, 0 for no bit set (the POSIX ffs convention)
static inline unsigned int bsf32(MD_U32 x) { return x ? ctz32(x) + 1u : 0u; }
static inline unsigned int bsf64(MD_U64 x) { return x ? ctz64(x) + 1u : 0u; }
static inline unsigned int bsr32(MD_U32 x) { return 32u - clz32(x); }
static inline unsigned int bsr64(MD_U64 x) { return 64u - clz64(x); }

// The index of the first zero nibble / byte, counted from the least significant end; the width in nibbles / bytes when
// there is none
static inline unsigned int find_first_zero_nibble16(unsigned short x) {
    const MD_U32 y = (MD_U32)((x - 0x1111u) & ~(MD_U32)x & 0x8888u);
    return ctz32(y | 0x10000u) >> 2;
}
static inline unsigned int find_first_zero_nibble32(MD_U32 x) {
    const MD_U32 y = (x - 0x11111111u) & ~x & 0x88888888u;
    return ctz32(y) >> 2;
}
static inline unsigned int find_first_zero_nibble64(MD_U64 x) {
    const MD_U64 y = (x - 0x1111111111111111ull) & ~x & 0x8888888888888888ull;
    return ctz64(y) >> 2;
}
static inline unsigned int find_first_zero_byte32(MD_U32 x) {
    const MD_U32 y = (x - 0x01010101u) & ~x & 0x80808080u;
    return ctz32(y) >> 3;
}
static inline unsigned int find_first_zero_byte64(MD_U64 x) {
    const MD_U64 y = (x - 0x0101010101010101ull) & ~x & 0x8080808080808080ull;
    return ctz64(y) >> 3;
}

// The smallest power of two >= x (x for 0 and 1)
static inline MD_U32 next_power_of_two32(MD_U32 x) { return (x < 2) ? x : (MD_U32)1 << (32 - clz32(x - 1)); }
static inline MD_U64 next_power_of_two64(MD_U64 x) { return (x < 2) ? x : (MD_U64)1 << (64 - clz64(x - 1)); }

// The 128 bit product of two 64 bit integers: the low half returned, the high half in *hi
static inline MD_U64 md_umul128(MD_U64 a, MD_U64 b, MD_U64* hi) {
#if (MD_COMPILER_GCC || MD_COMPILER_CLANG) && defined(__SIZEOF_INT128__)
    __extension__ typedef unsigned __int128 md_u128_t;    // __extension__: a GCC / Clang type, said so for -Wpedantic
    const md_u128_t r = (md_u128_t)a * b;
    *hi = (MD_U64)(r >> 64);
    return (MD_U64)r;
#elif MD_COMPILER_MSVC
    *hi = __umulh(a, b);
    return a * b;
#else
    const MD_U64 a_lo = (MD_U32)a, a_hi = a >> 32;
    const MD_U64 b_lo = (MD_U32)b, b_hi = b >> 32;
    const MD_U64 p0 = a_lo * b_lo, p1 = a_lo * b_hi, p2 = a_hi * b_lo, p3 = a_hi * b_hi;
    const MD_U64 mid = (p0 >> 32) + (MD_U32)p1 + (MD_U32)p2;
    *hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
    return (mid << 32) | (MD_U32)p0;
#endif
}

// ### Atomics
// On a naturally aligned 64-bit integer shared between threads without a lock, sequentially consistent. A fetch
// returns the value from before the operation.
#if MD_COMPILER_GCC || MD_COMPILER_CLANG
static inline MD_U64 md_atomic_fetch_add64(volatile MD_U64* ptr, MD_U64 value) {
    return __atomic_fetch_add(ptr, value, __ATOMIC_SEQ_CST);
}

static inline MD_U64 md_atomic_load64(const volatile MD_U64* ptr) {
    return __atomic_load_n(ptr, __ATOMIC_SEQ_CST);
}
#elif MD_COMPILER_MSVC
static inline MD_U64 md_atomic_fetch_add64(volatile MD_U64* ptr, MD_U64 value) {
    return (MD_U64)_InterlockedExchangeAdd64((volatile __int64*)ptr, (__int64)value);
}

static inline MD_U64 md_atomic_load64(const volatile MD_U64* ptr) {
    // A compare exchange which never exchanges: a full barrier, and a read of all 64 bits at once
    return (MD_U64)_InterlockedCompareExchange64((volatile __int64*)ptr, 0, 0);
}
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

#define ENUM_FLAGS(T) \
inline T operator|(T a, T b) { return (T)((int)a | (int)b); } \
inline T operator&(T a, T b) { return (T)((int)a & (int)b); } \
inline T operator^(T a, T b) { return (T)((int)a ^ (int)b); } \
inline T operator~(T a) { return (T)(~(int)a); } \
inline T& operator|=(T& a, T b) { return a = a | b; } \
inline T& operator&=(T& a, T b) { return a = a & b; } \
inline T& operator^=(T& a, T b) { return a = a ^ b; }
#else  //  __cplusplus
#define ENUM_FLAGS(T)
#endif //  __cplusplus
