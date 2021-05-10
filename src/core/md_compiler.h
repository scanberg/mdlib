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

#if defined(__clang__)
#	undef  MD_COMPILER_CLANG
#	define MD_COMPILER_CLANG (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)
#elif defined(_MSC_VER)
#	undef  MD_COMPILER_MSVC
#	define MD_COMPILER_MSVC _MSC_VER
#elif defined(__GNUC__)
#	undef  MD_COMPILER_GCC
#	define MD_COMPILER_GCC (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
#	error "MD_COMPILER_* is not defined!"
#endif //

#endif