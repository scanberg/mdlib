#ifndef _MD_COMMON_H_
#define _MD_COMMON_H_

#include "compiler.h"

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#define STATIC_ASSERT static_assert
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

#define internal static

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define ABS(x) ((x) > 0 ? (x) : -(x))
#define CLAMP(val, _min, _max) MIN(MAX(val, _min), _max)

#define DEG_TO_RAD(x) ((x)*(3.14159265f / 180.0f))
#define RAD_TO_DEG(x) ((x)*(180.0f / 3.14159265f))

#define KILOBYTES(x) (x * 1024)
#define MEGABYTES(x) (x * 1024 * 1024)

#endif