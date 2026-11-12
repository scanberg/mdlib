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

// Everything here is prefixed: mdlib is included into other projects, where an unprefixed MIN, PI or DEBUG would
// collide with theirs. md_ for what is called like a function and behaves like one - functions, and macros that are
// only another name for one (md_memcpy, md_bswap32) - and MD_ for macros that do not: MD_MIN evaluates an argument
// twice, MD_ASSERT's expression is gone in a release build, MD_COUNTOF wants an array, MD_PI is a constant. Nothing is behind an #ifndef guard either - a guard
// silently adopts whatever definition came first - except MD_ASSERT, which a project may define before including
// this to route assertions to its own handler.

#ifdef __cplusplus
#   define MD_STATIC_ASSERT static_assert
#else
#   define MD_STATIC_ASSERT _Static_assert
#endif

#if MD_COMPILER_MSVC
#   define MD_THREAD_LOCAL  __declspec(thread)
#   define MD_ALIGNAS(x)    __declspec(align(x))
#   define MD_FORCE_INLINE  __forceinline
#   define MD_NOINLINE      __declspec(noinline)
#else
#   define MD_THREAD_LOCAL  _Thread_local
#   define MD_ALIGNAS(x)    __attribute__((aligned(x)))
#   define MD_FORCE_INLINE  __attribute__((always_inline)) inline
#   define MD_NOINLINE      __attribute__((noinline))
#endif

#ifdef __cplusplus
#define MD_ALIGNOF(T) alignof(T)
#else
#define MD_ALIGNOF(T) _Alignof(T)
#endif

#if defined(__cplusplus)
#   define MD_FALLTHROUGH [[fallthrough]]
#elif MD_COMPILER_MSVC
#   define MD_FALLTHROUGH
#else
#   define MD_FALLTHROUGH __attribute__((fallthrough))
#endif

// Has GCC and Clang check the arguments of a printf style function against its format string
#if defined(__GNUC__) || defined(__clang__)
#   define MD_PRINTF_FORMAT(fmt_idx, arg_idx) __attribute__((format(printf, fmt_idx, arg_idx)))
#else
#   define MD_PRINTF_FORMAT(fmt_idx, arg_idx)
#endif

#define MD_RESTRICT  __restrict
#define MD_TYPEOF(x) __typeof__(x)

// Really GCC? REALLY? DO WE REALLY NEED TO INCLUDE stddef.h for this??????
#ifndef NULL
#define NULL 0
#endif

// 1 in the build they name, 0 in the other: test with #if, never #ifdef
#ifdef NDEBUG
#   define MD_DEBUG   0
#   define MD_RELEASE 1
#else
#   define MD_DEBUG   1
#   define MD_RELEASE 0
#endif

// The number of elements of an array (not a pointer: in C++ that does not compile)
#ifdef __cplusplus
    template <typename T, unsigned long long N> constexpr unsigned long long md_array_size_impl(T (&)[N]) { return N; }
#   define MD_COUNTOF(x) md_array_size_impl(x)
#else
#   define MD_COUNTOF(x) (sizeof(x)/sizeof(x[0]))
#endif

// State an invariant the optimiser cannot derive on its own. Unlike MD_ASSERT this is NOT a check:
// nothing is tested at runtime in any build, and a false expression is undefined behaviour. Use it
// only where the invariant is guaranteed by construction and the compiler's inability to see it
// costs something concrete - a spurious diagnostic, or a branch that can never be taken.
#if MD_COMPILER_MSVC
#   define MD_ASSUME(expr) __assume(expr)
#elif MD_COMPILER_GCC || MD_COMPILER_CLANG
#   define MD_ASSUME(expr) ((expr) ? (void)0 : __builtin_unreachable())
#else
#   define MD_ASSUME(expr) ((void)0)
#endif

#define MD_TYPE_COMPATIBLE(x, type) _Generic((x), type : 1, default : 0)

#define MD_IS_POW2(x) (((x) & ((x) - 1)) == 0)

// This is essentially the same as MD_ROUND_UP, but only works when alignment is a power of two
#define MD_ALIGN_TO(x, alignment) (((x) + ((alignment)-1)) & (~(alignment)+1))
#define MD_IS_ALIGNED(ptr, alignment) (((uintptr_t)(ptr) % (alignment)) == 0)
#define MD_NEXT_ALIGNED_ADDRESS(ptr, alignment) ((void*)(MD_ALIGN_TO((uintptr_t)(ptr), (alignment))))

#define MD_MIN(x, y) ((x) < (y) ? (x) : (y))
#define MD_MAX(x, y) ((x) > (y) ? (x) : (y))
#define MD_CLAMP(val, _min, _max) MD_MIN(MD_MAX(val, _min), _max)

// Round up division
#define MD_DIV_UP(x, y) (((x) + ((y)-1)) / (y))

// Rounds up x to nearest multiple of y. No affiliation to the weed-killer
#define MD_ROUND_UP(x, y) ((y) * MD_DIV_UP(x,y))

// Rounds down x to nearest multiple of y.
#define MD_ROUND_DOWN(x, y) ((y) * ((x) / (y)))

#define MD_PI 3.1415926535897932
#define MD_TWO_PI (2.0 * MD_PI)
#define MD_DEG_TO_RAD(x) ((x)*(MD_PI / 180.0))
#define MD_RAD_TO_DEG(x) ((x)*(180.0 / MD_PI))

#define MD_KILOBYTES(x) ((x) * 1024LL)
#define MD_MEGABYTES(x) ((x) * 1024LL * 1024LL)
#define MD_GIGABYTES(x) ((x) * 1024LL * 1024LL * 1024LL)

#define MD_CONCAT_INTERNAL(x, y) x##y
#define MD_CONCAT(x, y) MD_CONCAT_INTERNAL(x, y)

#define MD_STRINGIFY_VAL(s) MD_STRINGIFY(s)
#define MD_STRINGIFY(s) #s

// ## Assertions
// Compiled out in release builds unless __FORCE_ASSERTIONS__ (CMake: MD_FORCE_ASSERTIONS). Define MD_ASSERT before
// including this header to replace it.
#ifndef MD_ASSERT

#   ifdef  __cplusplus
extern "C" void md_assert_impl(const char* file, int line, const char* func_name, const char* expr);
#   else
void md_assert_impl(const char* file, int line, const char* func_name, const char* expr);
#   endif

#   if MD_DEBUG || __FORCE_ASSERTIONS__
#       define MD_ASSERT(__e) \
            ((__e) \
                ? (void)0 \
                : md_assert_impl( \
                    __FILE__, \
                    __LINE__, \
                    __func__, \
                    #__e) \
            )
#   else
#       define MD_ASSERT(__e)
#   endif

#endif

// ## Intrinsics
//
// THE RULE: what resolves without an include on every supported compiler is here; what needs a header to get
// is not (the SIMD types and their operations, md_simd.h). GCC and Clang - clang-cl too, which is MD_COMPILER_CLANG -
// have builtins for all of it. MSVC gets each intrinsic from a declaration of our own, made an intrinsic with #pragma
// intrinsic, exactly as <intrin.h> and winnt.h declare it, so that a translation unit which includes those as well
// sees the same declarations.

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

#define md_memcpy  memcpy
#define md_memset  memset
#define md_memmove memmove
#define md_memcmp  memcmp
#define md_bswap16 _byteswap_ushort
#define md_bswap32 _byteswap_ulong
#define md_bswap64 _byteswap_uint64
#define md_debug_break() __debugbreak()

#elif MD_COMPILER_GCC || MD_COMPILER_CLANG

#define md_memcpy   __builtin_memcpy
#define md_memset   __builtin_memset
#define md_memmove  __builtin_memmove
#define md_memcmp   __builtin_memcmp
#define md_bswap16  __builtin_bswap16
#define md_bswap32  __builtin_bswap32
#define md_bswap64  __builtin_bswap64
#if defined(__x86_64__) || defined(__i386__)
#   define md_debug_break() __asm__ volatile("int3")
#else
#   define md_debug_break() __builtin_trap()
#endif

#endif

// ### Bits
// The same results on every compiler and every CPU, for every input: clz and ctz of 0 are the width, popcnt, bsf and
// bsr are defined throughout. (MSVC's __lzcnt is not used: on a CPU without LZCNT - Sandy and Ivy Bridge, AVX or not -
// it executes as BSR and gives the index of the top bit instead of the count. _BitScan* is the same everywhere.)

#if MD_COMPILER_GCC || MD_COMPILER_CLANG

static inline unsigned int md_popcnt32(MD_U32 x) { return (unsigned int)__builtin_popcount(x); }
static inline unsigned int md_popcnt64(MD_U64 x) { return (unsigned int)__builtin_popcountll(x); }

// Count of leading (most significant) zero bits; 32 / 64 for 0
static inline unsigned int md_clz32(MD_U32 x) { return x ? (unsigned int)__builtin_clz(x)   : 32u; }
static inline unsigned int md_clz64(MD_U64 x) { return x ? (unsigned int)__builtin_clzll(x) : 64u; }

// Count of trailing (least significant) zero bits; 32 / 64 for 0. With BMI enabled this is one TZCNT, without it a BSF
// and a CMOV: the compilers know the idiom.
static inline unsigned int md_ctz32(MD_U32 x) { return x ? (unsigned int)__builtin_ctz(x)   : 32u; }
static inline unsigned int md_ctz64(MD_U64 x) { return x ? (unsigned int)__builtin_ctzll(x) : 64u; }

#elif MD_COMPILER_MSVC

#if defined(_M_X64)
static inline unsigned int md_popcnt32(MD_U32 x) { return __popcnt(x); }
static inline unsigned int md_popcnt64(MD_U64 x) { return (unsigned int)__popcnt64(x); }
#else
// ARM64: a SWAR count, which needs no declaration of _CountOneBits to be trusted
static inline unsigned int md_popcnt64(MD_U64 x) {
    x = x - ((x >> 1) & 0x5555555555555555ull);
    x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Full;
    return (unsigned int)((x * 0x0101010101010101ull) >> 56);
}
static inline unsigned int md_popcnt32(MD_U32 x) { return md_popcnt64(x); }
#endif

static inline unsigned int md_clz32(MD_U32 x) { unsigned long i; return _BitScanReverse  (&i, x) ? 31u - (unsigned int)i : 32u; }
static inline unsigned int md_clz64(MD_U64 x) { unsigned long i; return _BitScanReverse64(&i, x) ? 63u - (unsigned int)i : 64u; }
static inline unsigned int md_ctz32(MD_U32 x) { unsigned long i; return _BitScanForward  (&i, x) ? (unsigned int)i : 32u; }
static inline unsigned int md_ctz64(MD_U64 x) { unsigned long i; return _BitScanForward64(&i, x) ? (unsigned int)i : 64u; }

#endif

// The 1-based index of the lowest (bsf) and highest (bsr) set bit, 0 for no bit set (the POSIX ffs convention)
static inline unsigned int md_bsf32(MD_U32 x) { return x ? md_ctz32(x) + 1u : 0u; }
static inline unsigned int md_bsf64(MD_U64 x) { return x ? md_ctz64(x) + 1u : 0u; }
static inline unsigned int md_bsr32(MD_U32 x) { return 32u - md_clz32(x); }
static inline unsigned int md_bsr64(MD_U64 x) { return 64u - md_clz64(x); }

// The index of the first zero nibble / byte, counted from the least significant end; the width in nibbles / bytes when
// there is none
static inline unsigned int md_find_first_zero_nibble16(unsigned short x) {
    const MD_U32 y = (MD_U32)((x - 0x1111u) & ~(MD_U32)x & 0x8888u);
    return md_ctz32(y | 0x10000u) >> 2;
}
static inline unsigned int md_find_first_zero_nibble32(MD_U32 x) {
    const MD_U32 y = (x - 0x11111111u) & ~x & 0x88888888u;
    return md_ctz32(y) >> 2;
}
static inline unsigned int md_find_first_zero_nibble64(MD_U64 x) {
    const MD_U64 y = (x - 0x1111111111111111ull) & ~x & 0x8888888888888888ull;
    return md_ctz64(y) >> 2;
}
static inline unsigned int md_find_first_zero_byte32(MD_U32 x) {
    const MD_U32 y = (x - 0x01010101u) & ~x & 0x80808080u;
    return md_ctz32(y) >> 3;
}
static inline unsigned int md_find_first_zero_byte64(MD_U64 x) {
    const MD_U64 y = (x - 0x0101010101010101ull) & ~x & 0x8080808080808080ull;
    return md_ctz64(y) >> 3;
}

// The smallest power of two >= x (x for 0 and 1)
static inline MD_U32 md_next_power_of_two32(MD_U32 x) { return (x < 2) ? x : (MD_U32)1 << (32 - md_clz32(x - 1)); }
static inline MD_U64 md_next_power_of_two64(MD_U64 x) { return (x < 2) ? x : (MD_U64)1 << (64 - md_clz64(x - 1)); }

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

// MD_DEFER { ... }; runs the block when the scope is left
struct md_defer_helper {
    template <typename T>
    struct scope {
        T lambda;
        scope(T lambda) : lambda(lambda) {}
        ~scope() { lambda(); }
        scope& operator=(const scope&) = delete;
    };

    template <typename T>
    scope<T> operator+(T t) {
        return t;
    }
};

#define MD_DEFER [[maybe_unused]] const auto& MD_CONCAT(md_defer_, __LINE__) = md_defer_helper() + [&]()

#define MD_ENUM_FLAGS(T) \
inline T operator|(T a, T b) { return (T)((int)a | (int)b); } \
inline T operator&(T a, T b) { return (T)((int)a & (int)b); } \
inline T operator^(T a, T b) { return (T)((int)a ^ (int)b); } \
inline T operator~(T a) { return (T)(~(int)a); } \
inline T& operator|=(T& a, T b) { return a = a | b; } \
inline T& operator&=(T& a, T b) { return a = a & b; } \
inline T& operator^=(T& a, T b) { return a = a ^ b; }
#else  //  __cplusplus
#define MD_ENUM_FLAGS(T)
#endif //  __cplusplus
