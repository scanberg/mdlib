#pragma once

#include "md_common.h"
#include "md_compiler.h"

#ifndef MD_VEC_MATH_USE_SIMD
#define MD_VEC_MATH_USE_SIMD 1
#endif

#if MD_VEC_MATH_USE_SIMD
#include "md_simd.h"
#endif

#ifndef MD_VEC_INLINE
#define MD_VEC_INLINE static FORCE_INLINE
#endif

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

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct vec2_t {
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
} vec2_t;

typedef struct vec3_t {
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
} vec3_t;

typedef struct vec4_t {
    union {
        struct {
            float x, y, z, w;
        };
        float elem[4];
#if MD_VEC_MATH_USE_SIMD
        md_128 m128;
#endif
    };
#ifdef __cplusplus
    float& operator[](size_t i)       { return elem[i]; }
    const float& operator[](size_t i) const { return elem[i]; }
#endif
} vec4_t;

typedef struct quat_t {
    union {
        struct {
            float x, y, z, w;
        };
        float elem[4];
#if MD_VEC_MATH_USE_SIMD
        md_128 m128;
#endif
    };
#ifdef __cplusplus
        float& operator[](size_t i)       { return elem[i]; }
        const float& operator[](size_t i) const { return elem[i]; }
#endif
} quat_t;

typedef struct mat2_t {
    union {
        float elem[2][2];
        vec2_t col[2];
    };
#ifdef __cplusplus
    vec2_t& operator[](size_t i)       { return col[i]; }
    const vec2_t& operator[](size_t i) const { return col[i]; }
#endif
} mat2_t;

typedef struct mat3_t {
    union {
        float elem[3][3];
        vec3_t col[3];
    };
#ifdef __cplusplus
    vec3_t& operator[](size_t i)       { return col[i]; }
    const vec3_t& operator[](size_t i) const { return col[i]; }
#endif
} mat3_t;

typedef struct mat4_t {
    union {
        float elem[4][4];
        vec4_t col[4];
    };
#ifdef __cplusplus
    vec4_t& operator[](size_t i)       { return col[i]; }
    const vec4_t& operator[](size_t i) const { return col[i]; }
#endif
} mat4_t;

typedef struct mat4x3_t {
	union {
		float elem[3][4];
		vec4_t col[3];
	};
    #ifdef __cplusplus
	vec4_t& operator[](size_t i)       { return col[i]; }
    const vec4_t& operator[](size_t i) const { return col[i]; }
    #endif
} mat4x3_t;

#if MD_COMPILER_CLANG
#	pragma clang diagnostic pop
#elif MD_COMPILER_GCC
#	pragma GCC diagnostic pop
#elif MD_COMPILER_MSVC
#	pragma warning(pop)
#endif

MD_VEC_INLINE float stepf(float edge, float x) { return (float)((x - edge) > 0); }
MD_VEC_INLINE double step(double edge, double x) { return (double)((x - edge) > 0); }

MD_VEC_INLINE float fractf(float x) { return x - floorf(x); }
MD_VEC_INLINE double fract(double x) { return x - floor(x); }

// This is the version which seems to result in the best codegen for all tested compilers (msvc, gcc, clang)
MD_VEC_INLINE int   signf(float x) { return (int)((x > 0.0f) - (x < 0.0f)); }
MD_VEC_INLINE int   sign(double x) { return (int)((x > 0.0) - (x < 0.0)); }

MD_VEC_INLINE float deperiodizef(float x, float r, float period) {
    if (period == 0.0f) return x;
    const float dx  = (x - r) / period;
    const float dxp = dx - roundf(dx);
    const float x_prim = r + dxp * period;
    return x_prim;
}

MD_VEC_INLINE double deperiodize(double x, double r, double period) {
    if (period == 0.0) return x;
    const double dx  = (x - r) / period;
    const double dxp = dx - round(dx);
    const double x_prim = r + dxp * period;
    return x_prim;
}

MD_VEC_INLINE double deperiodize2(double x, double r, double period, double r_period) {
    if (period == 0.0) return x;
    const double dx  = (x - r) * r_period;
    const double dxp = dx - round(dx);
    const double x_prim = r + dxp * period;
    return x_prim;
}

MD_VEC_INLINE float lerpf(float a, float b, float t) {
    //t = CLAMP(t, 0.0f, 1.0f);
    return a * (1.0f - t) + b * t;
}

// Collision with std::lerp in math.h on GCC (fucking hell) 
// MD_VEC_INLINE double lerp(double a, double b, double t) {
//    return a * (1.0 - t) + b * t;
//}

// Cardinal cubic spline
MD_VEC_INLINE float cubic_splinef(float p0, float p1, float p2, float p3, float t, float s) {
    const float v0 = (p2 - p0) * s;
    const float v1 = (p3 - p1) * s;
    const float t2 = t * t;
    const float t3 = t * t * t;
    return (2.0f * p1 - 2.0f * p2 + v0 + v1) * t3 + (-3.0f * p1 + 3.0f * p2 - 2.0f * v0 - v1) * t2 + v0 * t + p1;
}

// VEC2 OPERATIONS
MD_VEC_INLINE vec2_t vec2_from_vec3(vec3_t v) {
    vec2_t r = {v.x, v.y};
    return r;
}

MD_VEC_INLINE vec2_t vec2_from_vec4(vec4_t v) {
    vec2_t r = {v.x, v.y};
    return r;
}

MD_VEC_INLINE vec2_t vec2_add(vec2_t a, vec2_t b) {
    vec2_t res = {a.x + b.x, a.y + b.y};
    return res;
}

MD_VEC_INLINE vec2_t vec2_add_f(vec2_t a, float f) {
    vec2_t res = {a.x + f, a.y + f};
    return res;
}

MD_VEC_INLINE vec2_t vec2_sub(vec2_t a, vec2_t b) {
    vec2_t res = {a.x - b.x, a.y - b.y};
    return res;
}

MD_VEC_INLINE vec2_t vec2_sub_f(vec2_t a, float s) {
    vec2_t res = {a.x - s, a.y - s};
    return res;
}

MD_VEC_INLINE vec2_t vec2_mul(vec2_t a, vec2_t b) {
    vec2_t res = {a.x * b.x, a.y * b.y};
    return res;
}

MD_VEC_INLINE vec2_t vec2_mul_f(vec2_t a, float s) {
    vec2_t res = {a.x * s, a.y * s};
    return res;
}

MD_VEC_INLINE vec2_t vec2_madd(vec2_t a, vec2_t b, vec2_t c) {
    vec2_t res = {(a.x * b.x) + c.x, (a.y * b.y) + c.y};
    return res;
}

MD_VEC_INLINE vec2_t vec2_div(vec2_t a, vec2_t b) {
    vec2_t res = {a.x / b.x, a.y / b.y};
    return res;
}

MD_VEC_INLINE vec2_t vec2_div_f(vec2_t a, float s) {
    vec2_t res = {a.x / s, a.y / s};
    return res;
}

MD_VEC_INLINE float vec2_dot(vec2_t a, vec2_t b) {
    return a.x * b.x + a.y * b.y;
}

MD_VEC_INLINE float vec2_length(vec2_t v) {
    float l2 = vec2_dot(v,v);
    return sqrtf(l2);
}

MD_VEC_INLINE float vec2_dist(vec2_t a, vec2_t b) {
    vec2_t d = vec2_sub(a, b);
    return vec2_length(d);
}

static inline vec2_t vec2_normalize(vec2_t v) {
    float len = vec2_length(v);
    if (len > 1.0e-5) {
        vec2_t res = {v.x / len, v.y / len};
        return res;
    }

    vec2_t res = {0, 0};
    return res;
}

MD_VEC_INLINE vec2_t vec2_lerp(vec2_t a, vec2_t b, float t) {
    ASSERT(0 <= t && t <= 1);
    return vec2_add(vec2_mul_f(a, 1.0f - t), vec2_mul_f(b, t));
}

MD_VEC_INLINE vec2_t vec2_fract(vec2_t v) {
    v.x = fractf(v.x);
    v.y = fractf(v.y);
    return v;
}

// VEC3 OPERATIONS
MD_VEC_INLINE vec3_t vec3_zero(void) {
    vec3_t res = {0};
    return res;
}

MD_VEC_INLINE vec3_t vec3_set(float x, float y, float z) {
    vec3_t v = {x,y,z};
    return v;
}

MD_VEC_INLINE vec3_t vec3_set1(float x) {
    vec3_t v = {x,x,x};
    return v;
}

MD_VEC_INLINE vec3_t vec3_fract(vec3_t v) {
    v.x = fractf(v.x);
    v.y = fractf(v.y);
    v.z = fractf(v.z);
    return v;
}

MD_VEC_INLINE vec3_t vec3_abs(vec3_t v) {
    v.x = fabsf(v.x);
    v.y = fabsf(v.y);
    v.z = fabsf(v.z);
    return v;
}

MD_VEC_INLINE vec3_t vec3_less_than(vec3_t a, vec3_t b) {
    vec3_t v = {
        (float)(a.x < b.x),
        (float)(a.y < b.y),
        (float)(a.z < b.z),
    };
    return v;
}

MD_VEC_INLINE bool vec3_equal(vec3_t a, vec3_t b) {
    return (a.x == b.x && a.y == b.y && a.z == b.z);
}

MD_VEC_INLINE vec3_t vec3_clamp(vec3_t v, vec3_t min, vec3_t max) {
    v.x = CLAMP(v.x, min.x, max.x);
    v.y = CLAMP(v.y, min.y, max.y);
    v.z = CLAMP(v.z, min.z, max.z);
    return v;
}

MD_VEC_INLINE vec3_t vec3_clamp_f(vec3_t v, float min, float max) {
    v.x = CLAMP(v.x, min, max);
    v.y = CLAMP(v.y, min, max);
    v.z = CLAMP(v.z, min, max);
    return v;
}

MD_VEC_INLINE vec3_t vec3_from_vec2(vec2_t v, float z) {
    vec3_t r = {v.x, v.y, z};
    return r;
}

MD_VEC_INLINE vec3_t vec3_from_vec4(vec4_t v) {
    vec3_t r = {v.x, v.y, v.z};
    return r;
}

MD_VEC_INLINE vec3_t vec3_add(vec3_t a, vec3_t b) {
    vec3_t res = {a.x + b.x, a.y + b.y, a.z + b.z};
    return res;
}

MD_VEC_INLINE vec3_t vec3_add_f(vec3_t a, float f) {
    vec3_t res = {a.x + f, a.y + f, a.z + f};
    return res;
}

MD_VEC_INLINE vec3_t vec3_sub(vec3_t a, vec3_t b) {
    vec3_t res = {a.x - b.x, a.y - b.y, a.z - b.z};
    return res;
}

MD_VEC_INLINE vec3_t vec3_sub_f(vec3_t a, float f) {
    vec3_t res = {a.x - f, a.y - f, a.z - f};
    return res;
}

MD_VEC_INLINE vec3_t vec3_mul(vec3_t a, vec3_t b) {
    vec3_t res = {a.x * b.x, a.y * b.y, a.z * b.z};
    return res;
}

MD_VEC_INLINE vec3_t vec3_mul_f(vec3_t a, float f) {
    vec3_t res = {a.x * f, a.y * f, a.z * f};
    return res;
}

MD_VEC_INLINE vec3_t vec3_div(vec3_t a, vec3_t b) {
    vec3_t res = {a.x / b.x, a.y / b.y, a.z / b.z};
    return res;
}

MD_VEC_INLINE vec3_t vec3_div_f(vec3_t a, float f) {
    vec3_t res = {a.x / f, a.y / f, a.z / f};
    return res;
}

MD_VEC_INLINE vec3_t vec3_cross(vec3_t a, vec3_t b) {
    vec3_t res = {
        a.y * b.z - b.y * a.z,
        a.z * b.x - b.z * a.x,
        a.x * b.y - b.x * a.y};
    return res;
}

MD_VEC_INLINE float vec3_dot(vec3_t a, vec3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

MD_VEC_INLINE float vec3_length(vec3_t v) {
    return sqrtf(vec3_dot(v,v));
}

MD_VEC_INLINE float vec3_length_squared(vec3_t v) {
    return vec3_dot(v,v);
}

MD_VEC_INLINE float vec3_distance(vec3_t a, vec3_t b) {
    return vec3_length(vec3_sub(a, b));
}

MD_VEC_INLINE float vec3_distance_squared(vec3_t a, vec3_t b) {
    const vec3_t d = vec3_sub(a, b);
    return vec3_dot(d, d);
}

MD_VEC_INLINE vec3_t vec3_normalize(vec3_t v) {
    float len = vec3_length(v);
    if (len > 1.0e-5) {
        vec3_t res = {v.x / len, v.y / len, v.z / len};
        return res;
    }

    vec3_t res = {0, 0, 0};
    return res;
}

MD_VEC_INLINE vec3_t vec3_lerp(vec3_t a, vec3_t b, float t) {
    ASSERT(0 <= t && t <= 1);
    return vec3_add(vec3_mul_f(a, 1.0f - t), vec3_mul_f(b, t));
}

MD_VEC_INLINE vec3_t vec3_deperiodize(vec3_t v, vec3_t ref, vec3_t period) {
    vec3_t res = {
        deperiodizef(v.x, ref.x, period.x),
        deperiodizef(v.y, ref.y, period.y),
        deperiodizef(v.z, ref.z, period.z)
    };
    return res;
}

MD_VEC_INLINE vec3_t vec3_min(vec3_t a, vec3_t b) {
    vec3_t res = {MIN(a.x, b.x), MIN(a.y, b.y), MIN(a.z, b.z)};
    return res;
}

MD_VEC_INLINE vec3_t vec3_max(vec3_t a, vec3_t b) {
    vec3_t res = {MAX(a.x, b.x), MAX(a.y, b.y), MAX(a.z, b.z)};
    return res;
}

// VEC4 OPERATIONS
MD_VEC_INLINE vec4_t vec4_zero(void) {
#if MD_VEC_MATH_USE_SIMD
    vec4_t res;
    res.m128 = md_mm_setzero_ps();
    return res;
#else
    vec4_t res = {0};
    return res;
#endif
}

MD_VEC_INLINE vec4_t vec4_set(float x, float y, float z, float w) {
    vec4_t res;
#if MD_VEC_MATH_USE_SIMD
    res.m128 = md_mm_set_ps(w, z, y, x);
#else
    res.x = x;
    res.y = y;
    res.z = z;
    res.w = w;
#endif
    return res;
}

MD_VEC_INLINE vec4_t vec4_set1(float v) {
    vec4_t res;
#if MD_VEC_MATH_USE_SIMD
    res.m128 = md_mm_set1_ps(v);
#else
    res.x = v;
    res.y = v;
    res.z = v;
    res.w = v;
#endif
    return res;
}

MD_VEC_INLINE vec4_t vec4_from_float(float v) {
    return vec4_set1(v);
}

MD_VEC_INLINE vec4_t vec4_from_vec3(vec3_t v, float w) {
    vec4_t res;
#if MD_VEC_MATH_USE_SIMD
    res.m128 = md_mm_set_ps(w, v.z, v.y, v.x);
#else
    res = {v.x, v.y, v.z, w};
#endif
    return res;
}

MD_VEC_INLINE vec4_t vec4_mul(vec4_t a, vec4_t b) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_mul_ps(a.m128, b.m128);
#else
    c.x = a.x * b.x;
    c.y = a.y * b.y;
    c.z = a.z * b.z;
    c.w = a.w * b.w;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_mul_f(vec4_t a, float s) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_mul_ps(a.m128, md_mm_set1_ps(s));
#else
    c.x = a.x * s;
    c.y = a.y * s;
    c.z = a.z * s;
    c.w = a.w * s;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_div(vec4_t a, vec4_t b) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_div_ps(a.m128, b.m128);
#else
    c.x = a.x / b.x;
    c.y = a.y / b.y;
    c.z = a.z / b.z;
    c.w = a.w / b.w;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_div_f(vec4_t a, float s) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_div_ps(a.m128, md_mm_set1_ps(s));
#else
    c.x = a.x / s;
    c.y = a.y / s;
    c.z = a.z / s;
    c.w = a.w / s;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_add(vec4_t a, vec4_t b) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_add_ps(a.m128, b.m128);
#else
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    c.w = a.w + b.w;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_add_f(vec4_t a, float s) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_add_ps(a.m128, md_mm_set1_ps(s));
#else
    c.x = a.x + s;
    c.y = a.y + s;
    c.z = a.z + s;
    c.w = a.w + s;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_sub(vec4_t a, vec4_t b) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_sub_ps(a.m128, b.m128);
#else
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    c.w = a.w - b.w;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_sub_f(vec4_t a, float s) {
    vec4_t c;
#if MD_VEC_MATH_USE_SIMD
    c.m128 = md_mm_sub_ps(a.m128, md_mm_set1_ps(s));
#else
    c.x = a.x - s;
    c.y = a.y - s;
    c.z = a.z - s;
    c.w = a.w - s;
#endif
    return c;
}

MD_VEC_INLINE vec4_t vec4_fmadd(vec4_t a, vec4_t b, vec4_t c) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_fmadd_ps(a.m128, b.m128, c.m128);
#else
    r.x = a.x * b.x + c.x;
    r.y = a.y * b.y + c.y;
    r.z = a.z * b.z + c.z;
    r.w = a.w * b.w + c.w;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_fmsub(vec4_t a, vec4_t b, vec4_t c) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_fmsub_ps(a.m128, b.m128, c.m128);
#else
    r.x = a.x * b.x - c.x;
    r.y = a.y * b.y - c.y;
    r.z = a.z * b.z - c.z;
    r.w = a.w * b.w - c.w;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_sqrt(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_sqrt_ps(v.m128);
#else
    r.x = sqrtf(v.x);
    r.y = sqrtf(v.y);
    r.z = sqrtf(v.z);
    r.w = sqrtf(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_rsqrt(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_rsqrt_ps(v.m128);
    vec4_t mul = vec4_mul(vec4_mul(v, r), r);
    r = vec4_mul(vec4_mul(vec4_set1(0.5f), r), vec4_sub(vec4_set1(3.0f), mul));
#else
    r.x = 1.0f / sqrtf(v.x);
    r.y = 1.0f / sqrtf(v.y);
    r.z = 1.0f / sqrtf(v.z);
    r.w = 1.0f / sqrtf(v.w);
#endif
    return r;
}

MD_VEC_INLINE float vec4_dot(vec4_t a, vec4_t b) {
    vec4_t res = vec4_mul(a, b);
#if MD_VEC_MATH_USE_SIMD
    return md_mm_reduce_add_ps(res.m128);
#else
    return (res.x + res.y) + (res.z + res.w);
#endif
}

MD_VEC_INLINE float vec4_distance_squared(vec4_t a, vec4_t b) {
    vec4_t d = vec4_sub(a, b);
    return vec4_dot(d, d);
}

MD_VEC_INLINE float vec4_distance(vec4_t a, vec4_t b) {
    return sqrtf(vec4_distance_squared(a,b));
}

MD_VEC_INLINE float vec4_length_squared(vec4_t v) {
    return vec4_dot(v,v);
}

MD_VEC_INLINE float vec4_length(vec4_t v) {
    return sqrtf(vec4_dot(v,v));
}

MD_VEC_INLINE vec4_t vec4_normalize(vec4_t v) {
    float mag2 = vec4_dot(v,v);
    if (mag2 > 1.0e-5f) {
        v = vec4_mul_f(v, 1.0f / sqrtf(mag2));
    }
    return v;
}

MD_VEC_INLINE vec4_t vec4_lerp(vec4_t a, vec4_t b, float t) {
    return vec4_add(vec4_mul_f(a, 1.0f - t), vec4_mul_f(b, t));
}

MD_VEC_INLINE vec4_t vec4_cubic_spline(vec4_t p0, vec4_t p1, vec4_t p2, vec4_t p3, float t, float s) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    const vec4_t vt = vec4_set1(t);
    const vec4_t vs = vec4_set1(s);
    r.m128 = md_mm_cubic_spline_ps(p0.m128, p1.m128, p2.m128, p3.m128, vt.m128, vs.m128);
#else
    r.x = cubic_splinef(p0.x, p1.x, p2.x, p3.x, t, s);
    r.y = cubic_splinef(p0.y, p1.y, p2.y, p3.y, t, s);
    r.z = cubic_splinef(p0.z, p1.z, p2.z, p3.z, t, s);
    r.w = cubic_splinef(p0.w, p1.w, p2.w, p3.w, t, s);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_min(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = _mm_min_ps(a.m128, b.m128);
#else
    r.x = MIN(a.x, b.x);
    r.y = MIN(a.y, b.y);
    r.z = MIN(a.z, b.z);
    r.w = MIN(a.w, b.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_max(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = _mm_max_ps(a.m128, b.m128);
#else
    r.x = MAX(a.x, b.x);
    r.y = MAX(a.y, b.y);
    r.z = MAX(a.z, b.z);
    r.w = MAX(a.w, b.w);
#endif
    return r;
}

MD_VEC_INLINE float vec4_reduce_min(vec4_t v) {
#if MD_VEC_MATH_USE_SIMD
    return md_mm_reduce_min_ps(v.m128);
#else
    return MIN(MIN(v.x, v.y), MIN(v.z, v.w));
#endif
}

MD_VEC_INLINE float vec4_reduce_max(vec4_t v) {
#if MD_VEC_MATH_USE_SIMD
    return md_mm_reduce_max_ps(v.m128);
#else
    return MAX(MAX(v.x, v.y), MAX(v.z, v.w));
#endif
}

MD_VEC_INLINE float vec4_reduce_add(vec4_t v) {
#if MD_VEC_MATH_USE_SIMD
    return md_mm_reduce_add_ps(v.m128);
#else
    return (v.x + v.y) + (v.z + v.w);
#endif
}

MD_VEC_INLINE vec4_t vec4_clamp(vec4_t v, vec4_t min, vec4_t max) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = _mm_max_ps(_mm_min_ps(v.m128, max.m128), min.m128);
#else
    r.x = CLAMP(v.x, min.x, max.x);
    r.y = CLAMP(v.y, min.y, max.y);
    r.z = CLAMP(v.z, min.z, max.z);
    r.w = CLAMP(v.w, min.w, max.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_abs(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_abs_ps(v.m128);
#else
    r.x = ABS(v.x);
    r.y = ABS(v.y);
    r.z = ABS(v.z);
    r.w = ABS(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_fract(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_fract_ps(v.m128);
#else
    r.x = fractf(v.x);
    r.y = fractf(v.y);
    r.z = fractf(v.z);
    r.w = fractf(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_sign(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_sign_ps(v.m128);
#else
    r.x = signf(v.x);
    r.y = signf(v.y);
    r.z = signf(v.z);
    r.w = signf(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_round(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_round_ps(v.m128);
#else
    r.x = roundf(v.x);
    r.y = roundf(v.y);
    r.z = roundf(v.z);
    r.w = roundf(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_floor(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_floor_ps(v.m128);
#else
    r.x = roundf(v.x);
    r.y = roundf(v.y);
    r.z = roundf(v.z);
    r.w = roundf(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_ceil(vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_ceil_ps(v.m128);
#else
    r.x = roundf(v.x);
    r.y = roundf(v.y);
    r.z = roundf(v.z);
    r.w = roundf(v.w);
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_cmp_eq(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_cmpeq_ps(a.m128, b.m128);
#else
    r.x = (a.x == b.x) ? 1.0f : 0.0f;
    r.y = (a.y == b.y) ? 1.0f : 0.0f;
    r.z = (a.z == b.z) ? 1.0f : 0.0f;
    r.w = (a.w == b.w) ? 1.0f : 0.0f;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_cmp_lt(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_cmplt_ps(a.m128, b.m128);
#else
    r.x = (a.x < b.x) ? 1.0f : 0.0f;
    r.y = (a.y < b.y) ? 1.0f : 0.0f;
    r.z = (a.z < b.z) ? 1.0f : 0.0f;
    r.w = (a.w < b.w) ? 1.0f : 0.0f;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_cmp_le(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_cmple_ps(a.m128, b.m128);
#else
    r.x = (a.x <= b.x) ? 1.0f : 0.0f;
    r.y = (a.y <= b.y) ? 1.0f : 0.0f;
    r.z = (a.z <= b.z) ? 1.0f : 0.0f;
    r.w = (a.w <= b.w) ? 1.0f : 0.0f;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_cmp_gt(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_cmpgt_ps(a.m128, b.m128);
#else
    r.x = (a.x > b.x) ? 1.0f : 0.0f;
    r.y = (a.y > b.y) ? 1.0f : 0.0f;
    r.z = (a.z > b.z) ? 1.0f : 0.0f;
    r.w = (a.w > b.w) ? 1.0f : 0.0f;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_cmp_ge(vec4_t a, vec4_t b) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_cmpge_ps(a.m128, b.m128);
#else
    r.x = (a.x >= b.x) ? 1.0f : 0.0f;
    r.y = (a.y >= b.y) ? 1.0f : 0.0f;
    r.z = (a.z >= b.z) ? 1.0f : 0.0f;
    r.w = (a.w >= b.w) ? 1.0f : 0.0f;
#endif
    return r;
}

MD_VEC_INLINE bool vec4_equal(vec4_t a, vec4_t b) {
    bool r;
#if MD_VEC_MATH_USE_SIMD
    r = md_mm_movemask_ps(md_mm_cmpeq_ps(a.m128, b.m128)) == 0xF;
#else
    r = a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
#endif
    return r;
}

MD_VEC_INLINE vec4_t vec4_blend(vec4_t a, vec4_t b, vec4_t mask) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_blendv_ps(a.m128, b.m128, mask.m128);
#else
    r.x = (mask.x != 0.0f) ? a.x : b.x;
    r.y = (mask.y != 0.0f) ? a.y : b.y;
    r.z = (mask.z != 0.0f) ? a.z : b.z;
    r.w = (mask.w != 0.0f) ? a.w : b.w;
#endif
    return r;
}

#if MD_VEC_MATH_USE_SIMD
// We cannot invoke the SIMD version of this function directly, because the mask is not a constant expression when passed as a function argument.
#   ifdef _cplusplus
#       define vec4_blend_mask(a, b, mask) (vec4_t {.m128 = md_mm_blend_ps((a).m128, (b).m128, mask)})
#   else
#       define vec4_blend_mask(a, b, mask) ((vec4_t) {.m128 = md_mm_blend_ps((a).m128, (b).m128, mask)})
#   endif
#else
MD_VEC_INLINE vec4_t vec4_blend_mask(vec4_t a, vec4_t b, const int mask) {
    vec4_t r = {
        (mask & 1) ? a.x : b.x,
        (mask & 2) ? a.y : b.y,
        (mask & 4) ? a.z : b.z,
        (mask & 8) ? a.w : b.w,
    };
    return r;
}
#endif


MD_VEC_INLINE vec4_t vec4_deperiodize(vec4_t x, vec4_t r, vec4_t period) {
    const vec4_t dx     = vec4_div(vec4_sub(x, r), period);
    const vec4_t dxp    = vec4_sub(dx, vec4_round(dx));
    const vec4_t x_prim = vec4_add(r, vec4_mul(dxp, period));
    return vec4_blend(x_prim, x, vec4_cmp_eq(period, vec4_zero()));
}

MD_VEC_INLINE float vec4_periodic_distance_squared(vec4_t a, vec4_t b, vec4_t period) {
    const vec4_t dx = vec4_sub(a, b);
    const vec4_t d  = vec4_sub(dx, vec4_mul(vec4_round(vec4_div(dx, period)), period));
    const vec4_t dp = vec4_blend(d, dx, vec4_cmp_eq(period, vec4_zero()));
    return vec4_length_squared(dp);
}

MD_VEC_INLINE float vec4_periodic_distance(vec4_t a, vec4_t b, vec4_t period) {
    const vec4_t dx = vec4_sub(a, b);
    const vec4_t d  = vec4_sub(dx, vec4_mul(vec4_round(vec4_div(dx, period)), period));
    const vec4_t dp = vec4_blend(d, dx, vec4_cmp_eq(period, vec4_zero()));
    return vec4_length(dp);
}

MD_VEC_INLINE vec4_t vec4_from_u32(uint32_t rgba) {
    vec4_t v = {(float)((rgba >> 0) & 0xFF), (float)((rgba >> 8) & 0xFF), (float)((rgba >> 16) & 0xFF), (float)((rgba >> 24) & 0xFF)};
    return vec4_mul_f(v, 1.0f / 255.0f);
}

MD_VEC_INLINE uint32_t u32_from_vec4(vec4_t v) {
    uint32_t out;
    out  = ((uint32_t)(CLAMP(v.x, 0.0f, 1.0f) * 255.0f + 0.5f)) << 0;
    out |= ((uint32_t)(CLAMP(v.y, 0.0f, 1.0f) * 255.0f + 0.5f)) << 8;
    out |= ((uint32_t)(CLAMP(v.z, 0.0f, 1.0f) * 255.0f + 0.5f)) << 16;
    out |= ((uint32_t)(CLAMP(v.w, 0.0f, 1.0f) * 255.0f + 0.5f)) << 24;
    return out;
}

MD_VEC_INLINE vec4_t vec4_splat_x(vec4_t v) {
    vec4_t result;
    result.m128 = md_mm_splat_ps(v.m128, 0);
    return result;
}

MD_VEC_INLINE vec4_t vec4_splat_y(vec4_t v) {
    vec4_t result;
    result.m128 = md_mm_splat_ps(v.m128, 1);
    return result;
}

MD_VEC_INLINE vec4_t vec4_splat_z(vec4_t v) {
    vec4_t result;
    result.m128 = md_mm_splat_ps(v.m128, 2);
    return result;
}

MD_VEC_INLINE vec4_t vec4_splat_w(vec4_t v) {
    vec4_t result;
    result.m128 = md_mm_splat_ps(v.m128, 3);
    return result;
}

MD_VEC_INLINE void vec4_sincos(vec4_t x, vec4_t* s, vec4_t* c) {
#if MD_VEC_MATH_USE_SIMD
    md_mm_sincos_ps(x.m128, &s->m128, &c->m128);
#else
    s->x = sinf(x.x);
	s->y = sinf(x.y);
	s->z = sinf(x.z);
	s->w = sinf(x.w);
	c->x = cosf(x.x);
	c->y = cosf(x.y);
	c->z = cosf(x.z);
	c->w = cosf(x.w);
#endif
}

MD_VEC_INLINE quat_t quat_ident(void) {
	quat_t q = {0, 0, 0, 1};
	return q;
}

// quat
MD_VEC_INLINE float quat_dot(quat_t a, quat_t b) {
#if MD_VEC_MATH_USE_SIMD
    return md_mm_reduce_add_ps(md_mm_mul_ps(a.m128, b.m128));
#else
    return (a.x * b.x + a.y * b.y) + (a.z * b.z + a.w * b.w);
#endif
}

MD_VEC_INLINE quat_t quat_mul(quat_t a, quat_t b) {
    quat_t c;
#if MD_VEC_MATH_USE_SIMD
    // @NOTE: Reversed notation used here for ijkl, because it makes more sense in my reptile brain
    // shuffle from left to right (xyzw) x:0 y:1 z:2 w:3
#define MD_SHUFFLE(v,i,j,k,l) md_mm_shuffle_ps(v,v,_MM_SHUFFLE(l,k,j,i))
    __m128 t1 = _mm_mul_ps(MD_SHUFFLE(a.m128, 3,3,3,3), b.m128);
    __m128 t2 = _mm_mul_ps(MD_SHUFFLE(a.m128, 0,1,2,0), MD_SHUFFLE(b.m128, 3,3,3,0));
    __m128 t3 = _mm_mul_ps(MD_SHUFFLE(a.m128, 1,2,0,1), MD_SHUFFLE(b.m128, 2,0,1,1));
    __m128 t4 = _mm_mul_ps(MD_SHUFFLE(a.m128, 2,0,1,2), MD_SHUFFLE(b.m128, 1,2,0,2));
#undef MD_SHUFFLE
    t2 = _mm_mul_ps(t2, _mm_set_ps(-1,1,1,1));
    t3 = _mm_mul_ps(t3, _mm_set_ps(-1,1,1,1));
    c.m128 = _mm_add_ps(_mm_add_ps(t1, t2), _mm_sub_ps(t3, t4));
#else
    c.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
    c.y = a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z;
    c.z = a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x;
    c.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
#endif
    return c;
}

MD_VEC_INLINE quat_t quat_mul_f(quat_t q, float s) {
    quat_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = md_mm_mul_ps(q.m128, md_mm_set1_ps(s));
#else
    r.x = q.x * s,
    r.y = q.y * s,
    r.z = q.z * s,
    r.w = q.w * s,
#endif
    return r;
}

MD_VEC_INLINE vec3_t quat_mul_vec3(quat_t q, vec3_t v) {
    vec3_t u = {q.x, q.y, q.z};
    vec3_t t = vec3_mul_f(vec3_cross(u, v), 2.0f);
    vec3_t r = vec3_add(v, vec3_mul_f(t, q.w));
    return vec3_add(r, vec3_cross(u, t));
}

MD_VEC_INLINE quat_t quat_conj(quat_t q) {

    quat_t r = {
        -q.x,
        -q.y,
        -q.z,
         q.w,
    };
    return r;
}

MD_VEC_INLINE quat_t quat_axis_angle(vec3_t axis, float angle) {
    float half_angle = angle * 0.5f;
    float sin_angle = sinf(half_angle);
    quat_t q = {
        axis.x * sin_angle,
        axis.y * sin_angle,
        axis.z * sin_angle,
        cosf(half_angle),
    };
    return q;
}

MD_VEC_INLINE quat_t quat_normalize(quat_t q) {
    float mag2 = quat_dot(q,q);
    if (mag2 > 0.000001f) {
        q = quat_mul_f(q, 1.0f / sqrtf(mag2));
    }
    return q;
}

static inline quat_t quat_nlerp(quat_t qa, quat_t qb, float t) {
    quat_t res;

    if (quat_dot(qa, qb) < 0.0f) {
        qb = quat_conj(qb);
    }
#if MD_VEC_MATH_USE_SIMD
    res.m128 = md_mm_lerp_ps(qa.m128, qb.m128, t);
#else
    res.x = lerpf(qa.x, qb.x, t);
    res.y = lerpf(qa.y, qb.y, t);
    res.z = lerpf(qa.z, qb.z, t);
    res.w = lerpf(qa.w, qb.w, t);
#endif
    return quat_normalize(res);
}

// https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
static inline quat_t quat_slerp(quat_t qa, quat_t qb, float t) {
    quat_t qm;

    // Calculate angle between them.
    float cosHalfTheta = quat_dot(qa, qb);

    // if qa=qb or qa=-qb then theta = 0 and we can return qa
    if (fabsf(cosHalfTheta) >= 1.0) {
        qm.w = qa.w; qm.x = qa.x; qm.y = qa.y; qm.z = qa.z;
        return qm;
    }
    // Calculate temporary values.
    float halfTheta = acosf(cosHalfTheta);
    float sinHalfTheta = sqrtf(1.0f - cosHalfTheta * cosHalfTheta);
    // if theta = 180 degrees then result is not fully defined
    // we could rotate around any axis normal to qa or qb
    if (fabsf(sinHalfTheta) < 0.001f) { // fabs is floating point absolute
        qm.w = (qa.w * 0.5f + qb.w * 0.5f);
        qm.x = (qa.x * 0.5f + qb.x * 0.5f);
        qm.y = (qa.y * 0.5f + qb.y * 0.5f);
        qm.z = (qa.z * 0.5f + qb.z * 0.5f);
        return qm;
    }
    float ratioA = sinf((1 - t) * halfTheta) / sinHalfTheta;
    float ratioB = sinf(t * halfTheta) / sinHalfTheta;
    //calculate Quaternion.
    qm.w = (qa.w * ratioA + qb.w * ratioB);
    qm.x = (qa.x * ratioA + qb.x * ratioB);
    qm.y = (qa.y * ratioA + qb.y * ratioB);
    qm.z = (qa.z * ratioA + qb.z * ratioB);
    return qm;
}

static inline quat_t quat_from_mat3(mat3_t M) {
    quat_t q;
    q.w = sqrtf(MAX(0, 1 + M.elem[0][0] + M.elem[1][1] + M.elem[2][2])) * 0.5f;
    q.x = sqrtf(MAX(0, 1 + M.elem[0][0] - M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.y = sqrtf(MAX(0, 1 - M.elem[0][0] + M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.z = sqrtf(MAX(0, 1 - M.elem[0][0] - M.elem[1][1] + M.elem[2][2])) * 0.5f;

    q.x = copysignf(q.x, M.elem[2][1] - M.elem[1][2]);
    q.y = copysignf(q.y, M.elem[0][2] - M.elem[2][0]);
    q.z = copysignf(q.z, M.elem[1][0] - M.elem[0][1]);
    return q;
}

static inline quat_t quat_from_mat4(mat4_t M) {
    quat_t q;
    q.w = sqrtf(MAX(0, 1 + M.elem[0][0] + M.elem[1][1] + M.elem[2][2])) * 0.5f;
    q.x = sqrtf(MAX(0, 1 + M.elem[0][0] - M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.y = sqrtf(MAX(0, 1 - M.elem[0][0] + M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.z = sqrtf(MAX(0, 1 - M.elem[0][0] - M.elem[1][1] + M.elem[2][2])) * 0.5f;

    q.x = copysignf(q.x, M.elem[2][1] - M.elem[1][2]);
    q.y = copysignf(q.y, M.elem[0][2] - M.elem[2][0]);
    q.z = copysignf(q.z, M.elem[1][0] - M.elem[0][1]);
    return q;
}

#if MD_VEC_MATH_USE_SIMD
MD_VEC_INLINE md_128 linear_combine_4(md_128 a, md_128 B[4]) {
    md_128 res;
    res = md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x00), B[0]);
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x55), B[1]));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0xaa), B[2]));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0xff), B[3]));
    return res;
}

MD_VEC_INLINE md_128 linear_combine_3(md_128 a, md_128 B[3]) {
    md_128 res;
    res = md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x00), B[0]);
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x55), B[1]));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0xaa), B[2]));
    return res;
}

#endif

// MAT2
MD_VEC_INLINE mat2_t mat2_ident(void) {
    mat2_t M = {
        1,0,
        0,1,
    };
    return M;
}

// MAT3
MD_VEC_INLINE mat3_t mat3_ident(void) {
    mat3_t M = {
        1,0,0,
        0,1,0,
        0,0,1,
    };
    return M;
}

MD_VEC_INLINE bool mat3_equal(mat3_t A, mat3_t B) {
    return MEMCMP(A.elem, B.elem, sizeof(mat3_t)) == 0;
}

MD_VEC_INLINE mat3_t mat3_from_mat4(mat4_t M) {
    mat3_t R = {
        M.elem[0][0], M.elem[0][1], M.elem[0][2],
        M.elem[1][0], M.elem[1][1], M.elem[1][2],
        M.elem[2][0], M.elem[2][1], M.elem[2][2],
    };
    return R;
}

static inline mat3_t mat3_from_quat(quat_t q) {
    float yy2 = 2.0f * q.elem[1] * q.elem[1];
    float xy2 = 2.0f * q.elem[0] * q.elem[1];
    float xz2 = 2.0f * q.elem[0] * q.elem[2];
    float yz2 = 2.0f * q.elem[1] * q.elem[2];
    float zz2 = 2.0f * q.elem[2] * q.elem[2];
    float wz2 = 2.0f * q.elem[3] * q.elem[2];
    float wy2 = 2.0f * q.elem[3] * q.elem[1];
    float wx2 = 2.0f * q.elem[3] * q.elem[0];
    float xx2 = 2.0f * q.elem[0] * q.elem[0];

    mat3_t M;
    M.elem[0][0] = - yy2 - zz2 + 1.0f;
    M.elem[0][1] = xy2 + wz2;
    M.elem[0][2] = xz2 - wy2;
    M.elem[1][0] = xy2 - wz2;
    M.elem[1][1] = - xx2 - zz2 + 1.0f;
    M.elem[1][2] = yz2 + wx2;
    M.elem[2][0] = xz2 + wy2;
    M.elem[2][1] = yz2 - wx2;
    M.elem[2][2] = - xx2 - yy2 + 1.0f;
    return M;
}

// Create scale matrix from scalar components x, y, z which dictates the scale on each corresponding axis
MD_VEC_INLINE mat3_t mat3_scale(float x, float y, float z) {
    mat3_t M = {
        x,0,0,
        0,y,0,
        0,0,z,
    };
    return M;
}

MD_VEC_INLINE mat3_t mat3_add(mat3_t A, mat3_t B) {
    mat3_t M;
    M.col[0] = vec3_add(A.col[0], B.col[0]);
    M.col[1] = vec3_add(A.col[1], B.col[1]);
    M.col[2] = vec3_add(A.col[2], B.col[2]);
    return M;
}

MD_VEC_INLINE mat3_t mat3_sub(mat3_t A, mat3_t B) {
    mat3_t M;
    M.col[0] = vec3_sub(A.col[0], B.col[0]);
    M.col[1] = vec3_sub(A.col[1], B.col[1]);
    M.col[2] = vec3_sub(A.col[2], B.col[2]);
    return M;
}

MD_VEC_INLINE vec3_t mat3_mul_vec3(mat3_t M, vec3_t v) {
    vec3_t res;
    res.x = M.elem[0][0] * v.x + M.elem[1][0] * v.y + M.elem[2][0] * v.z;
    res.y = M.elem[0][1] * v.x + M.elem[1][1] * v.y + M.elem[2][1] * v.z;
    res.z = M.elem[0][2] * v.x + M.elem[1][2] * v.y + M.elem[2][2] * v.z;
    return res;
}

MD_VEC_INLINE mat3_t mat3_mul(mat3_t A, mat3_t B) {
    mat3_t C;
#define MULT(col, row) A.elem[0][row] * B.elem[col][0] + A.elem[1][row] * B.elem[col][1] + A.elem[2][row] * B.elem[col][2]
    C.elem[0][0] = MULT(0, 0);
    C.elem[0][1] = MULT(0, 1);
    C.elem[0][2] = MULT(0, 2);

    C.elem[1][0] = MULT(1, 0);
    C.elem[1][1] = MULT(1, 1);
    C.elem[1][2] = MULT(1, 2);

    C.elem[2][0] = MULT(2, 0);
    C.elem[2][1] = MULT(2, 1);
    C.elem[2][2] = MULT(2, 2);
#undef MULT
    return C;
}

MD_VEC_INLINE mat3_t mat3_mul_f(mat3_t M, float s) {
    mat3_t R;
    R.col[0] = vec3_mul_f(M.col[0], s);
    R.col[1] = vec3_mul_f(M.col[1], s);
    R.col[2] = vec3_mul_f(M.col[2], s);    
    return R;
}

MD_VEC_INLINE mat3_t mat3_transpose(mat3_t M) {
    mat3_t T = {
        M.elem[0][0], M.elem[1][0], M.elem[2][0],
        M.elem[0][1], M.elem[1][1], M.elem[2][1],
        M.elem[0][2], M.elem[1][2], M.elem[2][2]
    };
    return T;
}

static inline mat3_t mat3_inverse(mat3_t M) {
    const float one_over_det = 1.0f / (
        + M.elem[0][0] * (M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2])
        - M.elem[1][0] * (M.elem[0][1] * M.elem[2][2] - M.elem[2][1] * M.elem[0][2])
        + M.elem[2][0] * (M.elem[0][1] * M.elem[1][2] - M.elem[1][1] * M.elem[0][2]));

    mat3_t I;

    I.elem[0][0] = + (M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2]) * one_over_det;
    I.elem[1][0] = - (M.elem[1][0] * M.elem[2][2] - M.elem[2][0] * M.elem[1][2]) * one_over_det;
    I.elem[2][0] = + (M.elem[1][0] * M.elem[2][1] - M.elem[2][0] * M.elem[1][1]) * one_over_det;
    I.elem[0][1] = - (M.elem[0][1] * M.elem[2][2] - M.elem[2][1] * M.elem[0][2]) * one_over_det;
    I.elem[1][1] = + (M.elem[0][0] * M.elem[2][2] - M.elem[2][0] * M.elem[0][2]) * one_over_det;
    I.elem[2][1] = - (M.elem[0][0] * M.elem[2][1] - M.elem[2][0] * M.elem[0][1]) * one_over_det;
    I.elem[0][2] = + (M.elem[0][1] * M.elem[1][2] - M.elem[1][1] * M.elem[0][2]) * one_over_det;
    I.elem[1][2] = - (M.elem[0][0] * M.elem[1][2] - M.elem[1][0] * M.elem[0][2]) * one_over_det;
    I.elem[2][2] = + (M.elem[0][0] * M.elem[1][1] - M.elem[1][0] * M.elem[0][1]) * one_over_det;

    return I;
}

MD_VEC_INLINE float mat3_determinant(mat3_t M) {
    float d = M.elem[0][0] * (M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2])
            - M.elem[1][0] * (M.elem[0][1] * M.elem[2][2] - M.elem[2][1] * M.elem[0][2])
            + M.elem[2][0] * (M.elem[0][1] * M.elem[1][2] - M.elem[1][1] * M.elem[0][2]);
    return d;
}

MD_VEC_INLINE float mat3_trace(mat3_t M) {
    return M.elem[0][0] + M.elem[1][1] + M.elem[2][2];
}

MD_VEC_INLINE float mat3_abs_trace(mat3_t M) {
    return fabsf(M.elem[0][0]) + fabsf(M.elem[1][1]) + fabsf(M.elem[2][2]);
}

MD_VEC_INLINE float mat3_abs_sum(mat3_t M) {
    return
        fabsf(M.elem[0][0]) + fabsf(M.elem[0][1]) + fabsf(M.elem[0][2]) +
        fabsf(M.elem[1][0]) + fabsf(M.elem[1][1]) + fabsf(M.elem[1][2]) +
        fabsf(M.elem[2][0]) + fabsf(M.elem[2][1]) + fabsf(M.elem[2][2]);
}

MD_VEC_INLINE vec3_t mat3_diag(mat3_t M) {
    vec3_t v = {M.elem[0][0], M.elem[1][1], M.elem[2][2]};
    return v;
}

// Construct a rotation matrix from and angle and axis (Expects axis to be normalized)
MD_VEC_INLINE mat3_t mat3_angle_axis(float angle, vec3_t axis) {
    float c = cosf(angle);
	float s = sinf(angle);
	float t = 1.0f - c;

	mat3_t M = {
		t * axis.x * axis.x + c, t * axis.x * axis.y - s * axis.z, t * axis.x * axis.z + s * axis.y,
		t * axis.x * axis.y + s * axis.z, t * axis.y * axis.y + c, t * axis.y * axis.z - s * axis.x,
		t * axis.x * axis.z - s * axis.y, t * axis.y * axis.z + s * axis.x, t * axis.z * axis.z + c
	};
	return M;
}

typedef struct mat3_eigen_t {
    mat3_t vectors;
    vec3_t values;
} mat3_eigen_t;

typedef struct mat3_svd_t {
    mat3_t U;
    mat3_t V;
    vec3_t s;
} mat3_svd_t;

mat3_eigen_t mat3_eigen(mat3_t M);
mat3_svd_t   mat3_svd(mat3_t M);
mat3_t       mat3_extract_rotation(mat3_t M);
mat3_t       mat3_orthonormalize(mat3_t M);

// Computes the covariance matrix for a set of coordinates with a given center of mass.
// x,y,z / xyz: coordinates
// w:           weights (optional)
// indices:     indices into coordinates and weights (optional)
// count:       number of coordinates or indices
// mean:        mean (com if coordinates)
mat3_t mat3_covariance_matrix(const float* x, const float* y, const float* z, const float* w, const int32_t* indices, size_t count, vec3_t mean);
mat3_t mat3_covariance_matrix_vec4(const vec4_t* xyzw, const int32_t* indices, size_t count, vec3_t mean);

// Computes the cross covariance matrix for two set of coordinates with given center of mass.
// The set of points are assumed to have equal length and if w is not NULL, the same weight.
// x[2],y[2],z[2] / xyz: coordinate streams
// w[2]:                 weights (optional)
// indices[2]:           indices into coordinates and weights (optional)
// count:                number of coordinates or indices
// mean[2]:              mean (com if coordinates)
mat3_t mat3_cross_covariance_matrix(const float* const x[2], const float* const y[2], const float* const z[2], const float* const w[2], const int32_t* const indices[2], size_t count, const vec3_t mean[2]);

// Compute the cross covariance matrix for two set of coordinates with given center of mass.
// xyzw[2]: coordinate + weights
// indices: indices into coordinates and weights (optional)
// count:   number of coordinates or indices
// mean[2]: mean (com if coordinates)
mat3_t mat3_cross_covariance_matrix_vec4(const vec4_t* const xyzw[2], const int32_t* const indices[2], size_t count, const vec3_t mean[2]);

// Computes the optimal rotation matrix that minimizes the RMSD between two sets of coordinates.
// The set of points are assumed to have equal length and if w is not NULL, the same weight.
// x[2],y[2],z[2] / xyz: coordinate streams
// w[2]:                 weights (optional)
// indices[2]:           indices into coordinates and weights (optional)
// count:                number of coordinates or indices
// com[2]:               center of mass
mat3_t mat3_optimal_rotation(const float* const x[2], const float* const y[2], const float* const z[2], const float* const w[2], const int32_t* const indices[2], size_t count, const vec3_t com[2]);

// Computes the optimal rotation matrix that minimizes the RMSD between two sets of coordinates.
// xyzw[2]: coordinate + weights
// indices: indices into coordinates and weights (optional)
// count:   number of coordinates or indices
// com[2]:  center of mass
mat3_t mat3_optimal_rotation_vec4(const vec4_t* const xyzw[2], const int32_t* const indices[2], size_t count, const vec3_t com[2]);

// MAT4
MD_VEC_INLINE mat4_t mat4_ident(void) {
    mat4_t M = {
        {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        }
    };
    return M;
}

MD_VEC_INLINE mat4_t mat4_from_mat3(const mat3_t M) {
    mat4_t R = {
        {
            M.col[0].x, M.col[0].y, M.col[0].z, 0,
            M.col[1].x, M.col[1].y, M.col[1].z, 0,
            M.col[2].x, M.col[2].y, M.col[2].z, 0,
            0, 0, 0, 1,
        }
    };
    return R;
}

static inline mat4_t mat4_from_quat(quat_t q) {
    float yy2 = 2.0f * q.elem[1] * q.elem[1];
    float xy2 = 2.0f * q.elem[0] * q.elem[1];
    float xz2 = 2.0f * q.elem[0] * q.elem[2];
    float yz2 = 2.0f * q.elem[1] * q.elem[2];
    float zz2 = 2.0f * q.elem[2] * q.elem[2];
    float wz2 = 2.0f * q.elem[3] * q.elem[2];
    float wy2 = 2.0f * q.elem[3] * q.elem[1];
    float wx2 = 2.0f * q.elem[3] * q.elem[0];
    float xx2 = 2.0f * q.elem[0] * q.elem[0];

    mat4_t M;
    M.elem[0][0] = - yy2 - zz2 + 1.0f;
    M.elem[0][1] = xy2 + wz2;
    M.elem[0][2] = xz2 - wy2;
    M.elem[0][3] = 0;
    M.elem[1][0] = xy2 - wz2;
    M.elem[1][1] = - xx2 - zz2 + 1.0f;
    M.elem[1][2] = yz2 + wx2;
    M.elem[1][3] = 0;
    M.elem[2][0] = xz2 + wy2;
    M.elem[2][1] = yz2 - wx2;
    M.elem[2][2] = - xx2 - yy2 + 1.0f;
    M.elem[2][3] = 0;
    M.elem[3][0] = M.elem[3][1] = M.elem[3][2] = 0;
    M.elem[3][3] = 1;
    return M;
}

// Create mat4 scaling matrix from scalars x, y, z which dictates the corresponding scaling factor for each axis
MD_VEC_INLINE mat4_t mat4_scale(float x, float y, float z) {
    mat4_t M = {
        {
            x,0,0,0,
            0,y,0,0,
            0,0,z,0,
            0,0,0,1
        }
    };
    return M;
}

// Create mat4 scaling matrix from scalar s which dictates the corresponding scaling factor for all axes
MD_VEC_INLINE mat4_t mat4_scale_f(float s) {
    mat4_t M = {
        {
            s,0,0,0,
            0,s,0,0,
            0,0,s,0,
            0,0,0,1
        }
    };
    return M;
}

// Create mat4 scaling matrix from scalar s which dictates the corresponding scaling factor for all axes
MD_VEC_INLINE mat4_t mat4_scale_vec3(vec3_t v) {
    mat4_t M = {
        {
            v.x,0,0,0,
            0,v.y,0,0,
            0,0,v.z,0,
            0,0,0,1
        }
    };
    return M;
}

// Create mat4 translation matrix from scalars x, y, z which dictates the corresponding translation for each axis
MD_VEC_INLINE mat4_t mat4_translate(float x, float y, float z) {
    mat4_t M = {
        {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            x,y,z,1
        }
    };
    return M;
}

// Create mat4 translation matrix from scalars x, y, z which dictates the corresponding translation for each axis
MD_VEC_INLINE mat4_t mat4_translate_vec3(vec3_t t) {
    mat4_t M = {
        {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            t.x,t.y,t.z,1
        }
    };
    return M;
}

MD_VEC_INLINE mat4_t mat4_add(mat4_t A, mat4_t B) {
    mat4_t M;
    M.col[0] = vec4_add(A.col[0], B.col[0]);
    M.col[1] = vec4_add(A.col[1], B.col[1]);
    M.col[2] = vec4_add(A.col[2], B.col[2]);
    M.col[3] = vec4_add(A.col[3], B.col[3]);
    return M;
}

MD_VEC_INLINE mat4_t mat4_sub(mat4_t A, mat4_t B) {
    mat4_t M;
    M.col[0] = vec4_sub(A.col[0], B.col[0]);
    M.col[1] = vec4_sub(A.col[1], B.col[1]);
    M.col[2] = vec4_sub(A.col[2], B.col[2]);
    M.col[3] = vec4_sub(A.col[3], B.col[3]);
    return M;
}

MD_VEC_INLINE mat4_t mat4_mul(mat4_t A, mat4_t B) {
    mat4_t C;
#if MD_VEC_MATH_USE_SIMD
    C.col[0].m128 = linear_combine_4(B.col[0].m128, &A.col[0].m128);
    C.col[1].m128 = linear_combine_4(B.col[1].m128, &A.col[0].m128);
    C.col[2].m128 = linear_combine_4(B.col[2].m128, &A.col[0].m128);
    C.col[3].m128 = linear_combine_4(B.col[3].m128, &A.col[0].m128);
#else
#define MULT(col, row) \
    A.elem[0][row] * B.elem[col][0] + A.elem[1][row] * B.elem[col][1] + A.elem[2][row] * B.elem[col][2] + A.elem[3][row] * B.elem[col][3]
    C.elem[0][0] = MULT(0, 0);
    C.elem[0][1] = MULT(0, 1);
    C.elem[0][2] = MULT(0, 2);
    C.elem[0][3] = MULT(0, 3);

    C.elem[1][0] = MULT(1, 0);
    C.elem[1][1] = MULT(1, 1);
    C.elem[1][2] = MULT(1, 2);
    C.elem[1][3] = MULT(1, 3);

    C.elem[2][0] = MULT(2, 0);
    C.elem[2][1] = MULT(2, 1);
    C.elem[2][2] = MULT(2, 2);
    C.elem[2][3] = MULT(2, 3);

    C.elem[3][0] = MULT(3, 0);
    C.elem[3][1] = MULT(3, 1);
    C.elem[3][2] = MULT(3, 2);
    C.elem[3][3] = MULT(3, 3);
#undef MULT
#endif
    return C;
}

MD_VEC_INLINE vec4_t mat4_mul_vec4(mat4_t M, vec4_t v) {
    vec4_t r;
#if MD_VEC_MATH_USE_SIMD
    r.m128 = linear_combine_4(v.m128, &M.col[0].m128);
#else
    r.x = M.elem[0][0] * v.x + M.elem[1][0] * v.y + M.elem[2][0] * v.z + M.elem[3][0] * v.w;
    r.y = M.elem[0][1] * v.x + M.elem[1][1] * v.y + M.elem[2][1] * v.z + M.elem[3][1] * v.w;
    r.z = M.elem[0][2] * v.x + M.elem[1][2] * v.y + M.elem[2][2] * v.z + M.elem[3][2] * v.w;
    r.w = M.elem[0][3] * v.x + M.elem[1][3] * v.y + M.elem[2][3] * v.z + M.elem[3][3] * v.w;
    ASSERT(false);
#endif
    return r;
}

MD_VEC_INLINE vec3_t mat4_mul_vec3(mat4_t M, vec3_t v, float w) {
    vec4_t r = {v.x, v.y, v.z, w};
    r = mat4_mul_vec4(M, r);
    return vec3_from_vec4(r);
}

MD_VEC_INLINE mat4_t mat4_mul_f(mat4_t M, float s) {
    mat4_t C = {0};
#if MD_VEC_MATH_USE_SIMD
    C.col[0] = vec4_mul_f(M.col[0], s);
    C.col[1] = vec4_mul_f(M.col[1], s);
    C.col[2] = vec4_mul_f(M.col[2], s);
    C.col[3] = vec4_mul_f(M.col[3], s);
#else
    C.elem[0][0] = M.elem[0][0] * s;
    C.elem[0][1] = M.elem[0][1] * s;
    C.elem[0][2] = M.elem[0][2] * s;
    C.elem[0][3] = M.elem[0][3] * s;

    C.elem[1][0] = M.elem[1][0] * s;
    C.elem[1][1] = M.elem[1][1] * s;
    C.elem[1][2] = M.elem[1][2] * s;
    C.elem[1][3] = M.elem[1][3] * s;

    C.elem[2][0] = M.elem[2][0] * s;
    C.elem[2][1] = M.elem[2][1] * s;
    C.elem[2][2] = M.elem[2][2] * s;
    C.elem[2][3] = M.elem[2][3] * s;

    C.elem[3][0] = M.elem[3][0] * s;
    C.elem[3][1] = M.elem[3][1] * s;
    C.elem[3][2] = M.elem[3][2] * s;
    C.elem[3][3] = M.elem[3][3] * s;
#endif
    return C;
}

mat4_t mat4_inverse(mat4_t M);

MD_VEC_INLINE mat4_t mat4_transpose(mat4_t M) {
    mat4_t T;

#if MD_VEC_MATH_USE_SIMD
    T = M;
    SIMDE_MM_TRANSPOSE4_PS(T.col[0].m128, T.col[1].m128, T.col[2].m128, T.col[3].m128);
#else
    T.elem[0][0] = M.elem[0][0];
    T.elem[0][1] = M.elem[1][0];
    T.elem[0][2] = M.elem[2][0];
    T.elem[0][3] = M.elem[3][0];

    T.elem[1][0] = M.elem[0][1];
    T.elem[1][1] = M.elem[1][1];
    T.elem[1][2] = M.elem[2][1];
    T.elem[1][3] = M.elem[3][1];

    T.elem[2][0] = M.elem[0][2];
    T.elem[2][1] = M.elem[1][2];
    T.elem[2][2] = M.elem[2][2];
    T.elem[2][3] = M.elem[3][2];

    T.elem[3][0] = M.elem[0][3];
    T.elem[3][1] = M.elem[1][3];
    T.elem[3][2] = M.elem[2][3];
    T.elem[3][3] = M.elem[3][3];
#endif
    return T;
}

MD_VEC_INLINE bool mat4_equal(mat4_t A, mat4_t B) {
    return vec4_equal(A.col[0], B.col[0]) && vec4_equal(A.col[1], B.col[1]) && vec4_equal(A.col[2], B.col[2]) && vec4_equal(A.col[3], B.col[3]);
}

vec3_t mat4_unproject(vec3_t window_coords, mat4_t inv_view_proj_mat, vec4_t viewport);

mat4_t mat4_look_at(vec3_t eye, vec3_t center, vec3_t up);

mat4_t mat4_ortho(float left, float right, float bottom, float top, float near, float far);
mat4_t mat4_ortho_inv(float left, float right, float bottom, float top, float near, float far);
mat4_t mat4_ortho_2d(float left, float right, float bottom, float top);
mat4_t mat4_ortho_2d_inv(float left, float right, float bottom, float top);

mat4_t mat4_persp(float fovy, float aspect, float near, float far);
mat4_t mat4_persp_inv(float fovy, float aspect, float near, float far);

mat4_t mat4_frustum(float left, float right, float bottom, float top, float near, float far);
mat4_t mat4_frustum_inv(float left, float right, float bottom, float top, float near, float far);

// These are routines for performing the same operation on many items
void vec3_batch_translate_inplace(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, size_t count, vec3_t translation);
void vec3_batch_translate(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, size_t count, vec3_t translation);

void mat3_batch_transform_inplace(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, size_t count, mat3_t transform);
void mat3_batch_transform(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, size_t count, mat3_t transform);

void mat4_batch_transform_inplace(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, float w_comp, size_t count, mat4_t transform);
void mat4_batch_transform(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, float w_comp, size_t count, mat4_t transform);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
static inline bool operator == (vec2_t a, vec2_t b){
    return a.x == b.x && a.y == b.y;
}

static inline bool operator != (vec2_t a, vec2_t b){
    return !(a == b);
}

static inline vec2_t operator + (vec2_t a, vec2_t b) {
    vec2_t c = {a.x + b.x, a.y + b.y};
    return c;
}

static inline vec2_t operator + (vec2_t a, float s) {
    vec2_t c = {a.x + s, a.y + s};
    return c;
}

static inline vec2_t operator + (float s, vec2_t a) {
    vec2_t c = {a.x + s, a.y + s};
    return c;
}

static inline vec2_t operator - (vec2_t a, vec2_t b) {
    vec2_t c = {a.x - b.x, a.y - b.y};
    return c;
}

static inline vec2_t operator - (vec2_t a, float s) {
    vec2_t c = {a.x - s, a.y - s};
    return c;
}

static inline vec2_t operator - (const vec2_t &a) {
    vec2_t c = {-a.x, -a.y};
    return c;
}

static inline vec2_t operator - (float s, vec2_t a) {
    vec2_t c = {s - a.x, s - a.y};
    return c;
}

static inline vec2_t operator * (vec2_t a, vec2_t b) {
    vec2_t c = {a.x * b.x, a.y * b.y};
    return c;
}

static inline vec2_t operator * (vec2_t a, float s) {
    vec2_t c = {a.x * s, a.y * s};
    return c;
}

static inline vec2_t operator / (vec2_t a, vec2_t b) {
    vec2_t c = {a.x / b.x, a.y / b.y};
    return c;
}

static inline vec2_t operator / (vec2_t a, float s) {
    vec2_t c = {a.x / s, a.y / s};
    return c;
}

static inline vec2_t operator / (float s, vec2_t a) {
    vec2_t c = {s / a.x, s / a.y};
    return c;
}

// vec3
static inline bool operator == (vec3_t a, vec3_t b){
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static inline bool operator != (vec3_t a, vec3_t b){
    return !(a == b);
}

static inline vec3_t operator + (vec3_t a, vec3_t b) {
    return vec3_add(a, b);
}

static inline vec3_t operator + (vec3_t v, float s) {
    return vec3_add_f(v, s);
}

static inline vec3_t operator + (float s, vec3_t v) {
    return vec3_add_f(v, s);
}

static inline vec3_t operator += (vec3_t& a, vec3_t b) {
    a = vec3_add(a, b);
    return a;
}

static inline vec3_t operator += (vec3_t& v, float s) {
    v = vec3_add_f(v, s);
    return v;
}

static inline vec3_t operator - (vec3_t a, vec3_t b) {
    vec3_t c = {a.x - b.x, a.y - b.y, a.z - b.z};
    return c;
}

static inline vec3_t operator - (vec3_t a, float s) {
    vec3_t c = {a.x - s, a.y - s, a.z - s};
    return c;
}

static inline vec3_t operator - (float s, vec3_t a) {
    vec3_t c = {s - a.x, s - a.y, s - a.z};
    return c;
}

static inline vec3_t operator - (const vec3_t &a) {
    vec3_t c = {-a.x, -a.y, -a.z};
    return c;
}

static inline vec3_t operator -= (vec3_t& a, vec3_t b) {
    a = vec3_sub(a, b);
    return a;
}

static inline vec3_t operator -= (vec3_t& v, float s) {
    v = vec3_sub_f(v, s);
    return v;
}

static inline vec3_t operator * (vec3_t a, vec3_t b) {
    vec3_t c = {a.x * b.x, a.y * b.y, a.z * b.z};
    return c;
}

static inline vec3_t operator * (vec3_t a, float s) {
    vec3_t c = {a.x * s, a.y * s, a.z * s};
    return c;
}

static inline vec3_t operator * (float s, vec3_t a) {
    vec3_t c = {a.x * s, a.y * s, a.z * s};
    return c;
}

static inline vec3_t operator *= (vec3_t& a, vec3_t b) {
    a = vec3_mul(a, b);
    return a;
}

static inline vec3_t operator *= (vec3_t& v, float s) {
    v = vec3_mul_f(v, s);
    return v;
}

static inline vec3_t operator / (vec3_t a, vec3_t b) {
    vec3_t c = {a.x / b.x, a.y / b.y, a.z / b.z};
    return c;
}

static inline vec3_t operator / (vec3_t a, float s) {
    vec3_t c = {a.x / s, a.y / s, a.z / s};
    return c;
}

static inline vec3_t operator / (float s, vec3_t a) {
    vec3_t c = {s / a.x, s / a.y, s / a.z};
    return c;
}

static inline vec3_t operator /= (vec3_t& a, vec3_t b) {
    a = vec3_div(a, b);
    return a;
}

static inline vec3_t operator /= (vec3_t& v, float s) {
    v = vec3_div_f(v, s);
    return v;
}

// vec4_t
static inline bool operator == (vec4_t a, vec4_t b){
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

static inline bool operator != (vec4_t a, vec4_t b){
    return !(a == b);
}

static inline vec4_t& operator += (vec4_t& a, vec4_t b) {
    a = vec4_add(a,b);
    return a;
}

static inline vec4_t& operator += (vec4_t& v, float s) {
    v = vec4_add_f(v,s);
    return v;
}

static inline vec4_t operator + (vec4_t a, vec4_t b) {
    return vec4_add(a, b);
}

static inline vec4_t operator + (vec4_t v, float s) {
    return vec4_add_f(v, s);
}

static inline vec4_t operator + (float s, vec4_t v) {
    return vec4_add_f(v, s);
}

static inline vec4_t& operator -= (vec4_t& a, vec4_t b) {
    a = vec4_sub(a,b);
    return a;
}

static inline vec4_t& operator -= (vec4_t& v, float s) {
    v = vec4_sub_f(v,s);
    return v;
}

static inline vec4_t operator - (vec4_t a, vec4_t b) {
    return vec4_sub(a, b);
}

static inline vec4_t operator - (vec4_t v, float s) {
    return vec4_sub_f(v, s);
}

static inline vec4_t operator - (float s, vec4_t v) {
    return vec4_sub(vec4_from_float(s), v);
}

static inline vec4_t operator - (const vec4_t &a) {
    vec4_t c = {-a.x, -a.y, -a.z, -a.w};
    return c;
}

static inline vec4_t& operator *= (vec4_t& a, vec4_t b) {
    a = vec4_mul(a,b);
    return a;
}

static inline vec4_t& operator *= (vec4_t& v, float s) {
    v = vec4_mul_f(v,s);
    return v;
}

static inline vec4_t operator * (vec4_t a, vec4_t b) {
    return vec4_mul(a, b);
}

static inline vec4_t operator * (vec4_t v, float s) {
    return vec4_mul_f(v, s);
}

static inline vec4_t operator * (float s, vec4_t v) {
    return vec4_mul_f(v, s);
}

static inline vec4_t& operator /= (vec4_t& a, vec4_t b) {
    a = vec4_div(a,b);
    return a;
}

static inline vec4_t& operator /= (vec4_t& v, float s) {
    v = vec4_div_f(v,s);
    return v;
}

static inline vec4_t operator / (vec4_t a, vec4_t b) {
    return vec4_div(a, b);
}

static inline vec4_t operator / (vec4_t v, float s) {
    return vec4_div_f(v, s);
}

static inline vec4_t operator / (float s, vec4_t v) {
    return vec4_div(vec4_from_float(s), v);
}

// quat

static inline quat_t operator * (quat_t a, quat_t b) {
    return quat_mul(a, b);
}

static inline vec3_t operator * (quat_t q, vec3_t v) {
    return quat_mul_vec3(q, v);
}

// mat

static inline mat3_t operator * (mat3_t A, mat3_t B) {
    return mat3_mul(A, B);
}

static inline mat3_t operator * (mat3_t A, float s) {
    return mat3_mul_f(A, s);
}

static inline mat3_t operator * (float s, mat3_t A) {
    return mat3_mul_f(A, s);
}

static inline mat3_t operator + (mat3_t A, mat3_t B) {
    return mat3_add(A, B);
}

static inline mat3_t operator - (mat3_t A, mat3_t B) {
    return mat3_sub(A, B);
}

static inline bool operator == (mat3_t A, mat3_t B) {
    return A[0] == B[0] && A[1] == B[1] && A[2] == B[2];
}

static inline bool operator != (mat3_t A, mat3_t B) {
    return !(A == B);
}

static inline vec3_t operator * (mat3_t M, vec3_t v) {
    return mat3_mul_vec3(M, v);
}

static inline mat4_t operator * (mat4_t A, mat4_t B) {
    return mat4_mul(A, B);
}

static inline mat4_t operator + (mat4_t A, mat4_t B) {
    return mat4_add(A, B);
}

static inline mat4_t operator - (mat4_t A, mat4_t B) {
    return mat4_sub(A, B);
}

static inline vec4_t operator * (mat4_t M, vec4_t v) {
    return mat4_mul_vec4(M, v);
}

static inline bool operator == (mat4_t A, mat4_t B) {
    return A[0] == B[0] && A[1] == B[1] && A[2] == B[2] && A[3] == B[3];
}

static inline bool operator != (mat4_t A, mat4_t B) {
    return !(A == B);
}

template<typename T, typename V>
T lerp (T a, T b, V t) {
    return (T)((1 - t) * a + t * b);
}

template <typename T, typename V>
T cubic_spline(T p0, T p1, T p2, T p3, V t, V s = (V)1.0) {
    T v0 = (p2 - p0) * s;
    T v1 = (p3 - p1) * s;
    V t2 = t * t;
    V t3 = t * t2;
    return ((V)2.0 * p1 - (V)2.0 * p2 + v0 + v1) * t3 + (-(V)3.0 * p1 + (V)3.0 * p2 - (V)2.0 * v0 - v1) * t2 + v0 * t + p1;
}

#endif
