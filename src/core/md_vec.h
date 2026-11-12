#pragma once

#include "md_common.h"
#include "md_compiler.h"

#include <stddef.h>

// The vector and matrix types, without the math on them (md_vec_math.h) or the SIMD that math is built on: what a
// header which only stores them includes.
//
// Named by element type and shape: md_float3_t, md_double3_t, md_float3x3_t. Functions take the name of the type they
// are of and, after an operation, that of the other operand: md_float3x3_mul_float3(M, v) is M * v, with the
// arguments in the order of the product. A function with a 1 after it takes a scalar (md_float3_mul1).
//
// Matrices are column major: col[i] (m[i] in C++) is column i, and elem[c][r] the element of column c, row r. Slang's
// m[i] is ROW i. A matrix named RxC has R rows and C columns, as in Slang: md_float4x3_t is 4 rows by 3 columns, stored as
// 3 columns of md_float4_t.
//
// The layout is fixed - these are in structs, files and GPU buffers - and the asserts at the end hold it to what it
// was when the SIMD registers were members: 16 byte alignment for the 4 wide float types (and the matrices of them)
// and md_double2_t, 32 for md_double4_t.

#if MD_COMPILER_GCC
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Wpedantic"
#elif MD_COMPILER_CLANG
#	pragma clang diagnostic push
#	pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#	pragma clang diagnostic ignored "-Wnested-anon-types"
#elif MD_COMPILER_MSVC
#	pragma warning(push)
#	pragma warning(disable: 4201)  // nonstandard extension used : nameless struct/union
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_float2_t {
    union {
        struct {
            float x, y;
        };
        float elem[2];
    };
#ifdef __cplusplus
    float& operator[](size_t i)       { return elem[i]; }
    const float& operator[](size_t i) const { return elem[i]; }
#endif
} md_float2_t;

typedef struct md_float3_t {
    union {
        struct {
            float x, y, z;
        };
        float elem[3];
    };
#ifdef __cplusplus
    float& operator[](size_t i)       { return elem[i]; }
    const float& operator[](size_t i) const { return elem[i]; }
#endif
} md_float3_t;

typedef struct MD_ALIGNAS(16) md_float4_t {
    union {
        struct {
            float x, y, z, w;
        };
        float elem[4];
    };
#ifdef __cplusplus
    float& operator[](size_t i)       { return elem[i]; }
    const float& operator[](size_t i) const { return elem[i]; }
#endif
} md_float4_t;

typedef struct MD_ALIGNAS(16) md_quat_t {
    union {
        struct {
            float x, y, z, w;
        };
        float elem[4];
    };
#ifdef __cplusplus
        float& operator[](size_t i)       { return elem[i]; }
        const float& operator[](size_t i) const { return elem[i]; }
#endif
} md_quat_t;

typedef struct md_float2x2_t {
    union {
        float elem[2][2];
        md_float2_t col[2];
    };
#ifdef __cplusplus
    md_float2_t& operator[](size_t i)       { return col[i]; }
    const md_float2_t& operator[](size_t i) const { return col[i]; }
#endif
} md_float2x2_t;

typedef struct md_float3x3_t {
    union {
        float elem[3][3];
        md_float3_t col[3];
    };
#ifdef __cplusplus
    md_float3_t& operator[](size_t i)       { return col[i]; }
    const md_float3_t& operator[](size_t i) const { return col[i]; }
#endif
} md_float3x3_t;

// ISO C does not implicitly add a qualifier through a pointer to array: float(*)[3] does NOT convert
// to const float(*)[3] the way float* converts to const float* (C11 6.5.16.1 / 6.7.6.1). So handing
// a md_float3x3_t's elem - or any float[3][3] - to a function taking const float[3][3] is a constraint
// violation that -Wpedantic reports. The cast only adds const, so it is always safe.
#define MD_AS_CONST_MAT3(m) ((const float (*)[3])(m))

typedef struct md_float4x4_t {
    union {
        float elem[4][4];
        md_float4_t col[4];
    };
#ifdef __cplusplus
    md_float4_t& operator[](size_t i)       { return col[i]; }
    const md_float4_t& operator[](size_t i) const { return col[i]; }
#endif
} md_float4x4_t;

// Braced-initialiser form of md_float4x4_ident(), for use as a C++ default member
// initialiser: 'md_float4x4_t m = MD_MAT4_IDENT_INIT;' instead of 'md_float4x4_t m = md_float4x4_ident();'.
//
// GCC 13 ICEs (gimplify.cc:774, in gimple_add_tmp_var) when a default member
// initialiser that needs a temporary -- a call returning a class type, or a copy
// from a constant -- is re-gimplified for a *second* value-initialised temporary
// of the enclosing type in the same translation unit. The front end shares one
// initialiser tree across all the use sites, and gimplifying it twice trips the
// assert. A braced initialiser constructs in place, needs no temporary, and is
// therefore immune. The value is identical.
#define MD_MAT4_IDENT_INIT {{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}}

typedef struct md_float4x3_t {
	union {
		float elem[3][4];
		md_float4_t col[3];
	};
    #ifdef __cplusplus
	md_float4_t& operator[](size_t i)       { return col[i]; }
    const md_float4_t& operator[](size_t i) const { return col[i]; }
    #endif
} md_float4x3_t;

typedef struct MD_ALIGNAS(16) md_double2_t {
    union {
        struct {
            double x, y;  
        };
        double elem[2];
    };
#ifdef __cplusplus
    double& operator[](size_t i)       { return elem[i]; }
    const double& operator[](size_t i) const { return elem[i]; }
#endif
} md_double2_t;

typedef struct md_double3_t {
    union {
        struct {
            double x, y, z;  
        };
        double elem[3];
    };
#ifdef __cplusplus
    double& operator[](size_t i)       { return elem[i]; }
    const double& operator[](size_t i) const { return elem[i]; }
#endif
} md_double3_t;

typedef struct MD_ALIGNAS(32) md_double4_t {
    union {
        struct {
            double x, y, z, w;  
        };
        double elem[4];
    };
#ifdef __cplusplus
    double& operator[](size_t i)       { return elem[i]; }
    const double& operator[](size_t i) const { return elem[i]; }
#endif
} md_double4_t;

#ifdef __cplusplus
}
#endif

#if MD_COMPILER_CLANG
#	pragma clang diagnostic pop
#elif MD_COMPILER_GCC
#	pragma GCC diagnostic pop
#elif MD_COMPILER_MSVC
#	pragma warning(pop)
#endif

MD_STATIC_ASSERT(sizeof(md_float2_t)   ==  8 && MD_ALIGNOF(md_float2_t)   ==  4, "md_float2_t layout");
MD_STATIC_ASSERT(sizeof(md_float3_t)   == 12 && MD_ALIGNOF(md_float3_t)   ==  4, "md_float3_t layout");
MD_STATIC_ASSERT(sizeof(md_float4_t)   == 16 && MD_ALIGNOF(md_float4_t)   == 16, "md_float4_t layout");
MD_STATIC_ASSERT(sizeof(md_quat_t)     == 16 && MD_ALIGNOF(md_quat_t)     == 16, "md_quat_t layout");
MD_STATIC_ASSERT(sizeof(md_float2x2_t) == 16 && MD_ALIGNOF(md_float2x2_t) ==  4, "md_float2x2_t layout");
MD_STATIC_ASSERT(sizeof(md_float3x3_t) == 36 && MD_ALIGNOF(md_float3x3_t) ==  4, "md_float3x3_t layout");
MD_STATIC_ASSERT(sizeof(md_float4x4_t) == 64 && MD_ALIGNOF(md_float4x4_t) == 16, "md_float4x4_t layout");
MD_STATIC_ASSERT(sizeof(md_float4x3_t) == 48 && MD_ALIGNOF(md_float4x3_t) == 16, "md_float4x3_t layout");
MD_STATIC_ASSERT(sizeof(md_double2_t)  == 16 && MD_ALIGNOF(md_double2_t)  == 16, "md_double2_t layout");
MD_STATIC_ASSERT(sizeof(md_double3_t)  == 24 && MD_ALIGNOF(md_double3_t)  ==  8, "md_double3_t layout");
MD_STATIC_ASSERT(sizeof(md_double4_t)  == 32 && MD_ALIGNOF(md_double4_t)  == 32, "md_double4_t layout");
