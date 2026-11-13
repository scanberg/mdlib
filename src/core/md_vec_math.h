#pragma once

#include "md_common.h"
#include "md_vec.h"
#include "md_simd.h"

#ifndef MD_VEC_INLINE
#define MD_VEC_INLINE static MD_FORCE_INLINE
#endif


#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

MD_VEC_INLINE float md_f32_step(float edge, float x) { return (float)((x - edge) > 0); }
MD_VEC_INLINE double md_f64_step(double edge, double x) { return (double)((x - edge) > 0); }

MD_VEC_INLINE float md_f32_fract(float x) { return x - floorf(x); }
MD_VEC_INLINE double md_f64_fract(double x) { return x - floor(x); }

// This is the version which seems to result in the best codegen for all tested compilers (msvc, gcc, clang)
MD_VEC_INLINE int   md_f32_sign(float x) { return (int)((x > 0.0f) - (x < 0.0f)); }
MD_VEC_INLINE int   md_f64_sign(double x) { return (int)((x > 0.0) - (x < 0.0)); }

MD_VEC_INLINE float md_f32_deperiodize_ortho(float x, float r, float period) {
    if (period == 0.0f) return x;
    const float dx  = (x - r) / period;
    const float dxp = dx - roundf(dx);
    const float x_prim = r + dxp * period;
    return x_prim;
}

MD_VEC_INLINE double md_f64_deperiodize_ortho(double x, double r, double period) {
    if (period == 0.0) return x;
    const double dx  = (x - r) / period;
    const double dxp = dx - round(dx);
    const double x_prim = r + dxp * period;
    return x_prim;
}

MD_VEC_INLINE double md_f64_deperiodize_ortho2(double x, double r, double period, double r_period) {
    if (period == 0.0) return x;
    const double dx  = (x - r) * r_period;
    const double dxp = MD_CLAMP(dx - round(dx), -0.5, 0.5); // As in md_float4_min_image
    const double x_prim = r + dxp * period;
    return x_prim;
}

MD_VEC_INLINE float md_f32_lerp(float a, float b, float t) {
    //t = MD_CLAMP(t, 0.0f, 1.0f);
    return a * (1.0f - t) + b * t;
}

MD_VEC_INLINE double md_f64_lerp(double a, double b, double t) {
    return a * (1.0 - t) + b * t;
}

// Cardinal cubic spline
MD_VEC_INLINE float md_f32_cubic_spline(float p0, float p1, float p2, float p3, float t, float s) {
    const float v0 = (p2 - p0) * s;
    const float v1 = (p3 - p1) * s;
    const float t2 = t * t;
    const float t3 = t * t * t;
    return (2.0f * p1 - 2.0f * p2 + v0 + v1) * t3 + (-3.0f * p1 + 3.0f * p2 - 2.0f * v0 - v1) * t2 + v0 * t + p1;
}

// VEC2 OPERATIONS
MD_VEC_INLINE md_float2_t md_float2_set(float x, float y) {
    md_float2_t v = {x, y};
    return v;
}

MD_VEC_INLINE md_float2_t md_float2_set1(float x) {
    md_float2_t v = {x, x};
    return v;
}

MD_VEC_INLINE md_float2_t md_float2_load(const float* ptr) {
    md_float2_t v;
    md_memcpy(v.elem, ptr, sizeof(md_float2_t));
    return v;
}

MD_VEC_INLINE void md_float2_store(float* ptr, md_float2_t v) {
    md_memcpy(ptr, v.elem, sizeof(md_float2_t));
}

MD_VEC_INLINE md_float2_t md_float2_from_float3(md_float3_t v) {
    md_float2_t r = {v.x, v.y};
    return r;
}

MD_VEC_INLINE md_float2_t md_float2_from_float4(md_float4_t v) {
    md_float2_t r = {v.x, v.y};
    return r;
}

MD_VEC_INLINE md_float2_t md_float2_add(md_float2_t a, md_float2_t b) {
    md_float2_t res = {a.x + b.x, a.y + b.y};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_add1(md_float2_t a, float f) {
    md_float2_t res = {a.x + f, a.y + f};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_sub(md_float2_t a, md_float2_t b) {
    md_float2_t res = {a.x - b.x, a.y - b.y};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_sub1(md_float2_t a, float s) {
    md_float2_t res = {a.x - s, a.y - s};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_mul(md_float2_t a, md_float2_t b) {
    md_float2_t res = {a.x * b.x, a.y * b.y};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_mul1(md_float2_t a, float s) {
    md_float2_t res = {a.x * s, a.y * s};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_madd(md_float2_t a, md_float2_t b, md_float2_t c) {
    md_float2_t res = {(a.x * b.x) + c.x, (a.y * b.y) + c.y};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_div(md_float2_t a, md_float2_t b) {
    md_float2_t res = {a.x / b.x, a.y / b.y};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_div1(md_float2_t a, float s) {
    md_float2_t res = {a.x / s, a.y / s};
    return res;
}

MD_VEC_INLINE float md_float2_dot(md_float2_t a, md_float2_t b) {
    return a.x * b.x + a.y * b.y;
}

MD_VEC_INLINE float md_float2_length(md_float2_t v) {
    float l2 = md_float2_dot(v,v);
    return sqrtf(l2);
}

MD_VEC_INLINE float md_float2_dist(md_float2_t a, md_float2_t b) {
    md_float2_t d = md_float2_sub(a, b);
    return md_float2_length(d);
}

static inline md_float2_t md_float2_normalize(md_float2_t v) {
    float len = md_float2_length(v);
    if (len > 1.0e-5) {
        md_float2_t res = {v.x / len, v.y / len};
        return res;
    }

    md_float2_t res = {0, 0};
    return res;
}

MD_VEC_INLINE md_float2_t md_float2_lerp(md_float2_t a, md_float2_t b, float t) {
    MD_ASSERT(0 <= t && t <= 1);
    return md_float2_add(md_float2_mul1(a, 1.0f - t), md_float2_mul1(b, t));
}

MD_VEC_INLINE md_float2_t md_float2_fract(md_float2_t v) {
    v.x = md_f32_fract(v.x);
    v.y = md_f32_fract(v.y);
    return v;
}

// VEC3 OPERATIONS
MD_VEC_INLINE md_float3_t md_float3_zero(void) {
    md_float3_t res = {0};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_set(float x, float y, float z) {
    md_float3_t v = {x,y,z};
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_set1(float x) {
    md_float3_t v = {x,x,x};
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_load(const float* v) {
    md_float3_t r;
    md_memcpy(r.elem, v, sizeof(md_float3_t));
    return r;
}

MD_VEC_INLINE void md_float3_store(float* v, md_float3_t r) {
    md_memcpy(v, r.elem, sizeof(md_float3_t));
}

MD_VEC_INLINE md_float3_t md_float3_fract(md_float3_t v) {
    v.x = md_f32_fract(v.x);
    v.y = md_f32_fract(v.y);
    v.z = md_f32_fract(v.z);
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_abs(md_float3_t v) {
    v.x = fabsf(v.x);
    v.y = fabsf(v.y);
    v.z = fabsf(v.z);
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_less_than(md_float3_t a, md_float3_t b) {
    md_float3_t v = {
        (float)(a.x < b.x),
        (float)(a.y < b.y),
        (float)(a.z < b.z),
    };
    return v;
}

MD_VEC_INLINE bool md_float3_equal(md_float3_t a, md_float3_t b) {
    return (a.x == b.x && a.y == b.y && a.z == b.z);
}

MD_VEC_INLINE md_float3_t md_float3_clamp(md_float3_t v, md_float3_t min, md_float3_t max) {
    v.x = MD_CLAMP(v.x, min.x, max.x);
    v.y = MD_CLAMP(v.y, min.y, max.y);
    v.z = MD_CLAMP(v.z, min.z, max.z);
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_clamp_f(md_float3_t v, float min, float max) {
    v.x = MD_CLAMP(v.x, min, max);
    v.y = MD_CLAMP(v.y, min, max);
    v.z = MD_CLAMP(v.z, min, max);
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_from_float2(md_float2_t v, float z) {
    md_float3_t r = {v.x, v.y, z};
    return r;
}

MD_VEC_INLINE md_float3_t md_float3_from_float4(md_float4_t v) {
    md_float3_t r = {v.x, v.y, v.z};
    return r;
}

MD_VEC_INLINE md_float3_t md_float3_add(md_float3_t a, md_float3_t b) {
    md_float3_t res = {a.x + b.x, a.y + b.y, a.z + b.z};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_add1(md_float3_t a, float f) {
    md_float3_t res = {a.x + f, a.y + f, a.z + f};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_sub(md_float3_t a, md_float3_t b) {
    md_float3_t res = {a.x - b.x, a.y - b.y, a.z - b.z};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_sub1(md_float3_t a, float f) {
    md_float3_t res = {a.x - f, a.y - f, a.z - f};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_mul(md_float3_t a, md_float3_t b) {
    md_float3_t res = {a.x * b.x, a.y * b.y, a.z * b.z};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_mul1(md_float3_t a, float f) {
    md_float3_t res = {a.x * f, a.y * f, a.z * f};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_div(md_float3_t a, md_float3_t b) {
    md_float3_t res = {a.x / b.x, a.y / b.y, a.z / b.z};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_div1(md_float3_t a, float f) {
    md_float3_t res = {a.x / f, a.y / f, a.z / f};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_cross(md_float3_t a, md_float3_t b) {
    md_float3_t res = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x};
    return res;
}

MD_VEC_INLINE float md_float3_dot(md_float3_t a, md_float3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

MD_VEC_INLINE float md_float3_length(md_float3_t v) {
    return sqrtf(md_float3_dot(v,v));
}

MD_VEC_INLINE float md_float3_length_squared(md_float3_t v) {
    return md_float3_dot(v,v);
}

MD_VEC_INLINE float md_float3_distance(md_float3_t a, md_float3_t b) {
    return md_float3_length(md_float3_sub(a, b));
}

MD_VEC_INLINE float md_float3_distance_squared(md_float3_t a, md_float3_t b) {
    const md_float3_t d = md_float3_sub(a, b);
    return md_float3_dot(d, d);
}

MD_VEC_INLINE md_float3_t md_float3_normalize(md_float3_t v) {
    float len = md_float3_length(v);
    if (len > 1.0e-5) {
        md_float3_t res = {v.x / len, v.y / len, v.z / len};
        return res;
    }

    md_float3_t res = {0, 0, 0};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_lerp(md_float3_t a, md_float3_t b, float t) {
    MD_ASSERT(0 <= t && t <= 1);
    return md_float3_add(md_float3_mul1(a, 1.0f - t), md_float3_mul1(b, t));
}

MD_VEC_INLINE md_float3_t md_float3_deperiodize_ortho(md_float3_t v, md_float3_t ref, md_float3_t period) {
    md_float3_t res = {
        md_f32_deperiodize_ortho(v.x, ref.x, period.x),
        md_f32_deperiodize_ortho(v.y, ref.y, period.y),
        md_f32_deperiodize_ortho(v.z, ref.z, period.z)
    };
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_min(md_float3_t a, md_float3_t b) {
    md_float3_t res = {MD_MIN(a.x, b.x), MD_MIN(a.y, b.y), MD_MIN(a.z, b.z)};
    return res;
}

MD_VEC_INLINE md_float3_t md_float3_max(md_float3_t a, md_float3_t b) {
    md_float3_t res = {MD_MAX(a.x, b.x), MD_MAX(a.y, b.y), MD_MAX(a.z, b.z)};
    return res;
}

MD_VEC_INLINE float md_float3_reduce_min(md_float3_t v) {
    return MD_MIN(v.x, MD_MIN(v.y, v.z));
}

MD_VEC_INLINE float md_float3_reduce_max(md_float3_t v) {
    return MD_MAX(v.x, MD_MAX(v.y, v.z));
}

MD_VEC_INLINE float md_float3_angle(md_float3_t v1, md_float3_t v2) {
    md_float3_t w = md_float3_cross(v1, v2);
    float wl = md_float3_length(w);
    float  s = md_float3_dot(v1, v2);
    return atan2f(wl, s);
}

// Each vector d is a delta vector between points, d[i] = p[i] - p[i-1]
// It is up to the caller to ensure that the resulting d1, d2 and d3 are constructed
// from points p0, p1, p2 and p4 that abide the minimum image convention
MD_VEC_INLINE float md_float3_dihedral_angle(md_float3_t d1, md_float3_t d2, md_float3_t d3) {
    md_float3_t v1 = md_float3_cross(d1, d2);
    md_float3_t v2 = md_float3_cross(d2, d3);
    float angle = md_float3_angle(v1, v2);
    float dot   = md_float3_dot(d1, v2);
    if (dot < 0) {
        angle = -angle;
    }
    return angle;
}

MD_VEC_INLINE md_float3_t md_float3_round(md_float3_t v) {
    v.x = roundf(v.x);
    v.y = roundf(v.y);
    v.z = roundf(v.z);
    return v;
}

MD_VEC_INLINE md_float3_t md_float3_min_image(md_float3_t dx, md_float3_t ext, md_float3_t inv_ext) {
    const md_float3_t d   = md_float3_mul(dx, inv_ext);
    const md_float3_t dxp = md_float3_clamp_f(md_float3_sub(d, md_float3_round(d)), -0.5f, 0.5f); // As in md_float4_min_image
    return md_float3_mul(dxp, ext);
}

// VEC4 OPERATIONS
MD_VEC_INLINE md_float4_t md_float4_zero(void) {
    md_float4_t res;
    md_mm_store_ps(res.elem, md_mm_setzero_ps());
    return res;
}

MD_VEC_INLINE md_float4_t md_float4_set(float x, float y, float z, float w) {
    md_float4_t res;
    md_mm_store_ps(res.elem, md_mm_set_ps(w, z, y, x));
    return res;
}

MD_VEC_INLINE md_float4_t md_float4_set1(float v) {
    md_float4_t res;
    md_mm_store_ps(res.elem, md_mm_set1_ps(v));
    return res;
}

MD_VEC_INLINE md_float4_t md_float4_load(const float* ptr) {
    md_float4_t v;
    md_memcpy(v.elem, ptr, sizeof(md_float4_t));
    return v;
}

MD_VEC_INLINE void md_float4_store(float* ptr, md_float4_t v) {
    md_memcpy(ptr, v.elem, sizeof(md_float4_t));
}

MD_VEC_INLINE md_float4_t md_float4_from_float(float v) {
    return md_float4_set1(v);
}

MD_VEC_INLINE md_float4_t md_float4_from_float3(md_float3_t v, float w) {
    md_float4_t res;
    md_mm_store_ps(res.elem, md_mm_set_ps(w, v.z, v.y, v.x));
    return res;
}

MD_VEC_INLINE md_float4_t md_float4_mul(md_float4_t a, md_float4_t b) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_mul_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_mul1(md_float4_t a, float s) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_mul_ps(md_mm_load_ps(a.elem), md_mm_set1_ps(s)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_div(md_float4_t a, md_float4_t b) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_div_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_div1(md_float4_t a, float s) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_div_ps(md_mm_load_ps(a.elem), md_mm_set1_ps(s)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_add(md_float4_t a, md_float4_t b) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_add_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_add1(md_float4_t a, float s) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_add_ps(md_mm_load_ps(a.elem), md_mm_set1_ps(s)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_sub(md_float4_t a, md_float4_t b) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_sub_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_sub1(md_float4_t a, float s) {
    md_float4_t c;
    md_mm_store_ps(c.elem, md_mm_sub_ps(md_mm_load_ps(a.elem), md_mm_set1_ps(s)));
    return c;
}

MD_VEC_INLINE md_float4_t md_float4_fmadd(md_float4_t a, md_float4_t b, md_float4_t c) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_fmadd_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem), md_mm_load_ps(c.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_fnmadd(md_float4_t a, md_float4_t b, md_float4_t c) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_fnmadd_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem), md_mm_load_ps(c.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_fmsub(md_float4_t a, md_float4_t b, md_float4_t c) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_fmsub_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem), md_mm_load_ps(c.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_fnmsub(md_float4_t a, md_float4_t b, md_float4_t c) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_fnmsub_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem), md_mm_load_ps(c.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_sqrt(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_sqrt_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_rsqrt(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_rsqrt_ps(md_mm_load_ps(v.elem)));
    // Perform additional Newton-Rhapson step to improve precision
    md_float4_t vr = md_float4_mul(v, r);
    md_float4_t hr = md_float4_mul(md_float4_set1(0.5f), r);
    r = md_float4_mul(hr, md_float4_fnmadd(vr, r, md_float4_set1(3.0f)));
    return r;
}

MD_VEC_INLINE float md_float4_dot(md_float4_t a, md_float4_t b) {
    md_float4_t res = md_float4_mul(a, b);
    return md_mm_reduce_add_ps(md_mm_load_ps(res.elem));
}

MD_VEC_INLINE float md_float4_distance_squared(md_float4_t a, md_float4_t b) {
    md_float4_t d = md_float4_sub(a, b);
    return md_float4_dot(d, d);
}

MD_VEC_INLINE float md_float4_distance(md_float4_t a, md_float4_t b) {
    return sqrtf(md_float4_distance_squared(a,b));
}

MD_VEC_INLINE float md_float4_length_squared(md_float4_t v) {
    return md_float4_dot(v,v);
}

MD_VEC_INLINE float md_float4_length(md_float4_t v) {
    return sqrtf(md_float4_dot(v,v));
}

MD_VEC_INLINE md_float4_t md_float4_normalize(md_float4_t v) {
    float mag2 = md_float4_dot(v,v);
    if (mag2 > 1.0e-5f) {
        v = md_float4_mul1(v, 1.0f / sqrtf(mag2));
    }
    return v;
}

// The cross product of the xyz of a and b, w 0
MD_VEC_INLINE md_float4_t md_float4_cross(md_float4_t a, md_float4_t b) {
    md_float4_t res;
    // (a * b.yzx - a.yzx * b).yzx: three shuffles where the textbook a.yzx * b.zxy - a.zxy * b.yzx takes four. Its w,
    // a.w * b.w - a.w * b.w, is not 0 in general once the multiply and the subtract fuse into one rounding: cleared.
    const md_128 a_yzx = md_mm_shuffle_ps(md_mm_load_ps(a.elem), md_mm_load_ps(a.elem), MD_SHUFFLE(3, 0, 2, 1));
    const md_128 b_yzx = md_mm_shuffle_ps(md_mm_load_ps(b.elem), md_mm_load_ps(b.elem), MD_SHUFFLE(3, 0, 2, 1));
    const md_128 c = md_mm_sub_ps(md_mm_mul_ps(md_mm_load_ps(a.elem), b_yzx), md_mm_mul_ps(a_yzx, md_mm_load_ps(b.elem)));
    md_mm_store_ps(res.elem, md_mm_blend_ps(md_mm_shuffle_ps(c, c, MD_SHUFFLE(3, 0, 2, 1)), md_mm_setzero_ps(), 0x8));
    return res;
}

MD_VEC_INLINE float md_float4_dihedral_angle(md_float4_t d1, md_float4_t d2, md_float4_t d3) {
    md_float4_t v1 = md_float4_cross(d1, d2);
    md_float4_t v2 = md_float4_cross(d2, d3);

    // normalize for stability
    v1 = md_float4_normalize(v1);
    v2 = md_float4_normalize(v2);

    float x = md_float4_dot(v1, v2);
    float y = md_float4_dot(md_float4_cross(v1, v2), md_float4_normalize(d2));

    return atan2f(y, x);
}

MD_VEC_INLINE float md_float4_angle(md_float4_t v1, md_float4_t v2) {
    md_float4_t w = md_float4_cross(v1, v2);
    float wl = md_float4_length(w);
    float s = md_float4_dot(v1, v2);
    return atan2f(wl, s);
}

MD_VEC_INLINE md_float4_t md_float4_lerp(md_float4_t a, md_float4_t b, float t) {
    return md_float4_add(md_float4_mul1(a, 1.0f - t), md_float4_mul1(b, t));
}

MD_VEC_INLINE md_float4_t md_float4_cubic_spline(md_float4_t p0, md_float4_t p1, md_float4_t p2, md_float4_t p3, float t, float s) {
    md_float4_t r;
    const md_float4_t vt = md_float4_set1(t);
    const md_float4_t vs = md_float4_set1(s);
    md_mm_store_ps(r.elem, md_mm_cubic_spline_ps(md_mm_load_ps(p0.elem), md_mm_load_ps(p1.elem), md_mm_load_ps(p2.elem), md_mm_load_ps(p3.elem), md_mm_load_ps(vt.elem), md_mm_load_ps(vs.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_min(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_min_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_max(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_max_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE float md_float4_reduce_min(md_float4_t v) {
    return md_mm_reduce_min_ps(md_mm_load_ps(v.elem));
}

MD_VEC_INLINE float md_float4_reduce_max(md_float4_t v) {
    return md_mm_reduce_max_ps(md_mm_load_ps(v.elem));
}

MD_VEC_INLINE float md_float4_reduce_add(md_float4_t v) {
    return md_mm_reduce_add_ps(md_mm_load_ps(v.elem));
}

MD_VEC_INLINE md_float4_t md_float4_clamp(md_float4_t v, md_float4_t min, md_float4_t max) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_max_ps(md_mm_min_ps(md_mm_load_ps(v.elem), md_mm_load_ps(max.elem)), md_mm_load_ps(min.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_abs(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_abs_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_fract(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_fract_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_sign(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_sign_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_neg(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_neg_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_round(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_round_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_floor(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_floor_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_ceil(md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_ceil_ps(md_mm_load_ps(v.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_cmp_eq(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_cmpeq_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_cmp_lt(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_cmplt_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_cmp_le(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_cmple_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_cmp_gt(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_cmpgt_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_cmp_ge(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_cmpge_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE bool md_float4_equal(md_float4_t a, md_float4_t b) {
    bool r;
    r = md_mm_movemask_ps(md_mm_cmpeq_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem))) == 0xF;
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_and(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_and_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

MD_VEC_INLINE md_float4_t md_float4_or(md_float4_t a, md_float4_t b) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_or_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
    return r;
}

// Construct a blend mask
MD_VEC_INLINE md_float4_t md_float4_mask(bool x, bool y, bool z, bool w) {
    float val;
    md_memset(&val, 0xFF, sizeof(val));
    md_float4_t r = {
        x ? val : 0.0f,
        y ? val : 0.0f,
        z ? val : 0.0f,
        w ? val : 0.0f,
    };
    return r;
}

// b where the mask is set and a where it is clear, for a mask whose lanes are all ones or all zeros, as
// md_float4_mask and the comparisons give them (md_mm_select_ps)
MD_VEC_INLINE md_float4_t md_float4_blend(md_float4_t a, md_float4_t b, md_float4_t mask) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_mm_select_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem), md_mm_load_ps(mask.elem)));
    return r;
}

// We cannot invoke the SIMD version of this function directly, because the mask is not a constant expression when passed as a function argument.
MD_VEC_INLINE md_float4_t md_float4__from_m128(md_128 v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, v);
    return r;
}
#   define md_float4_blend_mask(a, b, mask) md_float4__from_m128(md_mm_blend_ps(md_mm_load_ps((a).elem), md_mm_load_ps((b).elem), mask))

MD_VEC_INLINE int md_float4_move_mask(md_float4_t v) {
    return md_mm_movemask_ps(md_mm_load_ps(v.elem));
}

// dx brought within half a period of zero. The fraction of a period is clamped to [-1/2, 1/2]: a compiler that contracts
// may fuse the product into the subtraction of its own rounding and land an ulp past a half (see the deperiodize
// comment in md_simd.h).
MD_VEC_INLINE md_float4_t md_float4_min_image(md_float4_t dx, md_float4_t ext, md_float4_t inv_ext) {
    dx = md_float4_mul(dx, inv_ext);
    const md_float4_t dxp = md_float4_clamp(md_float4_sub(dx, md_float4_round(dx)), md_float4_set1(-0.5f), md_float4_set1(0.5f));
    return md_float4_mul(dxp, ext);
}

// x brought to the image nearest r, as md_mm_deperiodize2_ps: within half a period of r, in [0, ext] for r = ext / 2
MD_VEC_INLINE md_float4_t md_float4_deperiodize_ortho2(md_float4_t x, md_float4_t r, md_float4_t ext, md_float4_t inv_ext) {
    md_float4_t d = md_float4_fmsub(x, inv_ext, md_float4_mul(r, inv_ext));
    d = md_float4_sub(d, md_float4_round(d));
    const md_float4_t x_prim = md_float4_fmadd(d, ext, r);
    return md_float4_blend(x_prim, x, md_float4_cmp_eq(ext, md_float4_zero()));
}

MD_VEC_INLINE md_float4_t md_float4_deperiodize_ortho(md_float4_t x, md_float4_t r, md_float4_t ext) {
    return md_float4_deperiodize_ortho2(x, r, ext, md_float4_div(md_float4_set1(1.0f), ext));
}

MD_VEC_INLINE float md_float4_periodic_distance_squared(md_float4_t a, md_float4_t b, md_float4_t period) {
    const md_float4_t dx = md_float4_sub(a, b);
    const md_float4_t d  = md_float4_sub(dx, md_float4_mul(md_float4_round(md_float4_div(dx, period)), period));
    const md_float4_t dp = md_float4_blend(d, dx, md_float4_cmp_eq(period, md_float4_zero()));
    return md_float4_length_squared(dp);
}

MD_VEC_INLINE float md_float4_periodic_distance(md_float4_t a, md_float4_t b, md_float4_t period) {
    const md_float4_t dx = md_float4_sub(a, b);
    const md_float4_t d  = md_float4_sub(dx, md_float4_mul(md_float4_round(md_float4_div(dx, period)), period));
    const md_float4_t dp = md_float4_blend(d, dx, md_float4_cmp_eq(period, md_float4_zero()));
    return md_float4_length(dp);
}

MD_VEC_INLINE md_float4_t md_float4_from_u32(uint32_t rgba) {
    md_float4_t v = {(float)((rgba >> 0) & 0xFF), (float)((rgba >> 8) & 0xFF), (float)((rgba >> 16) & 0xFF), (float)((rgba >> 24) & 0xFF)};
    return md_float4_mul1(v, 1.0f / 255.0f);
}

MD_VEC_INLINE uint32_t md_u32_from_float4(md_float4_t v) {
    uint32_t out;
    out  = ((uint32_t)(MD_CLAMP(v.x, 0.0f, 1.0f) * 255.0f + 0.5f)) << 0;
    out |= ((uint32_t)(MD_CLAMP(v.y, 0.0f, 1.0f) * 255.0f + 0.5f)) << 8;
    out |= ((uint32_t)(MD_CLAMP(v.z, 0.0f, 1.0f) * 255.0f + 0.5f)) << 16;
    out |= ((uint32_t)(MD_CLAMP(v.w, 0.0f, 1.0f) * 255.0f + 0.5f)) << 24;
    return out;
}

MD_VEC_INLINE md_float4_t md_float4_splat_x(md_float4_t v) {
    md_float4_t result;
    md_mm_store_ps(result.elem, md_mm_splat_ps(md_mm_load_ps(v.elem), 0));
    return result;
}

MD_VEC_INLINE md_float4_t md_float4_splat_y(md_float4_t v) {
    md_float4_t result;
    md_mm_store_ps(result.elem, md_mm_splat_ps(md_mm_load_ps(v.elem), 1));
    return result;
}

MD_VEC_INLINE md_float4_t md_float4_splat_z(md_float4_t v) {
    md_float4_t result;
    md_mm_store_ps(result.elem, md_mm_splat_ps(md_mm_load_ps(v.elem), 2));
    return result;
}

MD_VEC_INLINE md_float4_t md_float4_splat_w(md_float4_t v) {
    md_float4_t result;
    md_mm_store_ps(result.elem, md_mm_splat_ps(md_mm_load_ps(v.elem), 3));
    return result;
}

// These are workhorses for linear algebra operations, e.g. matrix multiplications
MD_VEC_INLINE md_float4_t md_float4_linear_combine_4(md_float4_t a, const md_float4_t B[4]) {
    md_float4_t res;
    res = md_float4_mul(B[0], md_float4_splat_x(a));
    res = md_float4_add(res, md_float4_mul(B[1], md_float4_splat_y(a)));
    res = md_float4_add(res, md_float4_mul(B[2], md_float4_splat_z(a)));
    res = md_float4_add(res, md_float4_mul(B[3], md_float4_splat_w(a)));
    return res;
}

MD_VEC_INLINE md_float4_t md_float4_linear_combine_3(md_float4_t a, const md_float4_t B[3]) {
    md_float4_t res;
    res = md_float4_mul(B[0], md_float4_splat_x(a));
    res = md_float4_add(res, md_float4_mul(B[1], md_float4_splat_y(a)));
    res = md_float4_add(res, md_float4_mul(B[2], md_float4_splat_z(a)));
    return res;
}

MD_VEC_INLINE void md_float4_sincos(md_float4_t x, md_float4_t* s, md_float4_t* c) {
    md_128 sin_x, cos_x;
    md_mm_sincos_ps(md_mm_load_ps(x.elem), &sin_x, &cos_x);
    md_mm_store_ps(s->elem, sin_x);
    md_mm_store_ps(c->elem, cos_x);
}

MD_VEC_INLINE md_quat_t md_quat_ident(void) {
	md_quat_t q = {0, 0, 0, 1};
	return q;
}

// quat
MD_VEC_INLINE float md_quat_dot(md_quat_t a, md_quat_t b) {
    return md_mm_reduce_add_ps(md_mm_mul_ps(md_mm_load_ps(a.elem), md_mm_load_ps(b.elem)));
}

MD_VEC_INLINE md_quat_t md_quat_mul(md_quat_t a, md_quat_t b) {
    md_quat_t c;
    // @NOTE: Reversed notation used here for ijkl, because it makes more sense in my reptile brain
    // shuffle from left to right (xyzw) x:0 y:1 z:2 w:3
#define MD_PS_SHUFFLE(v,i,j,k,l) md_mm_shuffle_ps(v,v, MD_SHUFFLE(l,k,j,i))
    md_128 t1 = md_mm_mul_ps(MD_PS_SHUFFLE(md_mm_load_ps(a.elem), 3,3,3,3), md_mm_load_ps(b.elem));
    md_128 t2 = md_mm_mul_ps(MD_PS_SHUFFLE(md_mm_load_ps(a.elem), 0,1,2,0), MD_PS_SHUFFLE(md_mm_load_ps(b.elem), 3,3,3,0));
    md_128 t3 = md_mm_mul_ps(MD_PS_SHUFFLE(md_mm_load_ps(a.elem), 1,2,0,1), MD_PS_SHUFFLE(md_mm_load_ps(b.elem), 2,0,1,1));
    md_128 t4 = md_mm_mul_ps(MD_PS_SHUFFLE(md_mm_load_ps(a.elem), 2,0,1,2), MD_PS_SHUFFLE(md_mm_load_ps(b.elem), 1,2,0,2));
#undef MD_PS_SHUFFLE
    t2 = md_mm_mul_ps(t2, md_mm_set_ps(-1,1,1,1));
    t3 = md_mm_mul_ps(t3, md_mm_set_ps(-1,1,1,1));
    md_mm_store_ps(c.elem, md_mm_add_ps(md_mm_add_ps(t1, t2), md_mm_sub_ps(t3, t4)));
    return c;
}

MD_VEC_INLINE md_quat_t md_quat_mul1(md_quat_t q, float s) {
    md_quat_t r;
    md_mm_store_ps(r.elem, md_mm_mul_ps(md_mm_load_ps(q.elem), md_mm_set1_ps(s)));
    return r;
}

MD_VEC_INLINE md_float3_t md_quat_mul_float3(md_quat_t q, md_float3_t v) {
    md_float3_t u = {q.x, q.y, q.z};
    md_float3_t t = md_float3_mul1(md_float3_cross(u, v), 2.0f);
    md_float3_t r = md_float3_add(v, md_float3_mul1(t, q.w));
    return md_float3_add(r, md_float3_cross(u, t));
}

MD_VEC_INLINE md_quat_t md_quat_conj(md_quat_t q) {

    md_quat_t r = {
        -q.x,
        -q.y,
        -q.z,
         q.w,
    };
    return r;
}

MD_VEC_INLINE md_quat_t md_quat_axis_angle(md_float3_t axis, float angle) {
    float half_angle = angle * 0.5f;
    float sin_angle = sinf(half_angle);
    md_quat_t q = {
        axis.x * sin_angle,
        axis.y * sin_angle,
        axis.z * sin_angle,
        cosf(half_angle),
    };
    return q;
}

MD_VEC_INLINE md_quat_t md_quat_normalize(md_quat_t q) {
    float mag2 = md_quat_dot(q,q);
    if (mag2 > 0.000001f) {
        q = md_quat_mul1(q, 1.0f / sqrtf(mag2));
    }
    return q;
}

static inline md_quat_t md_quat_nlerp(md_quat_t qa, md_quat_t qb, float t) {
    md_quat_t res;

    if (md_quat_dot(qa, qb) < 0.0f) {
        qb = md_quat_conj(qb);
    }
    md_mm_store_ps(res.elem, md_mm_lerp_ps(md_mm_load_ps(qa.elem), md_mm_load_ps(qb.elem), t));
    return md_quat_normalize(res);
}

// https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
static inline md_quat_t md_quat_slerp(md_quat_t qa, md_quat_t qb, float t) {
    md_quat_t qm;

    // Calculate angle between them.
    float cosHalfTheta = md_quat_dot(qa, qb);

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

static inline md_quat_t md_quat_from_float3x3(md_float3x3_t M) {
    md_quat_t q;
    q.w = sqrtf(MD_MAX(0, 1 + M.elem[0][0] + M.elem[1][1] + M.elem[2][2])) * 0.5f;
    q.x = sqrtf(MD_MAX(0, 1 + M.elem[0][0] - M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.y = sqrtf(MD_MAX(0, 1 - M.elem[0][0] + M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.z = sqrtf(MD_MAX(0, 1 - M.elem[0][0] - M.elem[1][1] + M.elem[2][2])) * 0.5f;

    q.x = copysignf(q.x, M.elem[2][1] - M.elem[1][2]);
    q.y = copysignf(q.y, M.elem[0][2] - M.elem[2][0]);
    q.z = copysignf(q.z, M.elem[1][0] - M.elem[0][1]);
    return q;
}

static inline md_quat_t md_quat_from_float4x4(md_float4x4_t M) {
    md_quat_t q;
    q.w = sqrtf(MD_MAX(0, 1 + M.elem[0][0] + M.elem[1][1] + M.elem[2][2])) * 0.5f;
    q.x = sqrtf(MD_MAX(0, 1 + M.elem[0][0] - M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.y = sqrtf(MD_MAX(0, 1 - M.elem[0][0] + M.elem[1][1] - M.elem[2][2])) * 0.5f;
    q.z = sqrtf(MD_MAX(0, 1 - M.elem[0][0] - M.elem[1][1] + M.elem[2][2])) * 0.5f;

    q.x = copysignf(q.x, M.elem[2][1] - M.elem[1][2]);
    q.y = copysignf(q.y, M.elem[0][2] - M.elem[2][0]);
    q.z = copysignf(q.z, M.elem[1][0] - M.elem[0][1]);
    return q;
}

MD_VEC_INLINE md_128 md_float4__linear_combine_4(md_128 a, const md_float4_t B[4]) {
    md_128 res;
    res = md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x00), md_mm_load_ps(B[0].elem));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x55), md_mm_load_ps(B[1].elem)));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0xaa), md_mm_load_ps(B[2].elem)));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0xff), md_mm_load_ps(B[3].elem)));
    return res;
}

MD_VEC_INLINE md_128 md_float4__linear_combine_3(md_128 a, const md_float4_t B[3]) {
    md_128 res;
    res = md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x00), md_mm_load_ps(B[0].elem));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0x55), md_mm_load_ps(B[1].elem)));
    res = md_mm_add_ps(res, md_mm_mul_ps(md_mm_shuffle_ps(a, a, 0xaa), md_mm_load_ps(B[2].elem)));
    return res;
}


// MAT2
MD_VEC_INLINE md_float2x2_t md_float2x2_ident(void) {
    md_float2x2_t M = {
        1,0,
        0,1,
    };
    return M;
}

// MAT3
MD_VEC_INLINE md_float3x3_t md_float3x3_ident(void) {
    md_float3x3_t M = {
        1,0,0,
        0,1,0,
        0,0,1,
    };
    return M;
}

MD_VEC_INLINE md_float3x3_t md_float3x3_load(const float* ptr) {
    md_float3x3_t M;
    md_memcpy(M.elem, ptr, sizeof(md_float3x3_t));
    return M;
}

MD_VEC_INLINE void md_float3x3_store(float* ptr, md_float3x3_t M) {
    md_memcpy(ptr, M.elem, sizeof(md_float3x3_t));
}

MD_VEC_INLINE bool md_float3x3_equal(md_float3x3_t A, md_float3x3_t B) {
    return md_memcmp(A.elem, B.elem, sizeof(md_float3x3_t)) == 0;
}

MD_VEC_INLINE md_float3x3_t md_float3x3_from_float4x4(md_float4x4_t M) {
    md_float3x3_t R = {
        M.elem[0][0], M.elem[0][1], M.elem[0][2],
        M.elem[1][0], M.elem[1][1], M.elem[1][2],
        M.elem[2][0], M.elem[2][1], M.elem[2][2],
    };
    return R;
}

static inline md_float3x3_t md_float3x3_from_quat(md_quat_t q) {
    float yy2 = 2.0f * q.elem[1] * q.elem[1];
    float xy2 = 2.0f * q.elem[0] * q.elem[1];
    float xz2 = 2.0f * q.elem[0] * q.elem[2];
    float yz2 = 2.0f * q.elem[1] * q.elem[2];
    float zz2 = 2.0f * q.elem[2] * q.elem[2];
    float wz2 = 2.0f * q.elem[3] * q.elem[2];
    float wy2 = 2.0f * q.elem[3] * q.elem[1];
    float wx2 = 2.0f * q.elem[3] * q.elem[0];
    float xx2 = 2.0f * q.elem[0] * q.elem[0];

    md_float3x3_t M;
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
MD_VEC_INLINE md_float3x3_t md_float3x3_scale(float x, float y, float z) {
    md_float3x3_t M = {
        x,0,0,
        0,y,0,
        0,0,z,
    };
    return M;
}

MD_VEC_INLINE md_float3x3_t md_float3x3_add(md_float3x3_t A, md_float3x3_t B) {
    md_float3x3_t M;
    M.col[0] = md_float3_add(A.col[0], B.col[0]);
    M.col[1] = md_float3_add(A.col[1], B.col[1]);
    M.col[2] = md_float3_add(A.col[2], B.col[2]);
    return M;
}

MD_VEC_INLINE md_float3x3_t md_float3x3_sub(md_float3x3_t A, md_float3x3_t B) {
    md_float3x3_t M;
    M.col[0] = md_float3_sub(A.col[0], B.col[0]);
    M.col[1] = md_float3_sub(A.col[1], B.col[1]);
    M.col[2] = md_float3_sub(A.col[2], B.col[2]);
    return M;
}

MD_VEC_INLINE md_float3_t md_float3x3_mul_float3(md_float3x3_t M, md_float3_t v) {
    md_float3_t res;
    res.x = M.elem[0][0] * v.x + M.elem[1][0] * v.y + M.elem[2][0] * v.z;
    res.y = M.elem[0][1] * v.x + M.elem[1][1] * v.y + M.elem[2][1] * v.z;
    res.z = M.elem[0][2] * v.x + M.elem[1][2] * v.y + M.elem[2][2] * v.z;
    return res;
}

MD_VEC_INLINE md_float3x3_t md_float3x3_mul(md_float3x3_t A, md_float3x3_t B) {
    md_float3x3_t C;
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

MD_VEC_INLINE md_float3x3_t md_float3x3_mul1(md_float3x3_t M, float s) {
    md_float3x3_t R;
    R.col[0] = md_float3_mul1(M.col[0], s);
    R.col[1] = md_float3_mul1(M.col[1], s);
    R.col[2] = md_float3_mul1(M.col[2], s);    
    return R;
}

MD_VEC_INLINE md_float3x3_t md_float3x3_transpose(md_float3x3_t M) {
    md_float3x3_t T = {
        M.elem[0][0], M.elem[1][0], M.elem[2][0],
        M.elem[0][1], M.elem[1][1], M.elem[2][1],
        M.elem[0][2], M.elem[1][2], M.elem[2][2]
    };
    return T;
}

static inline md_float3x3_t md_float3x3_inverse(md_float3x3_t M) {
    const float one_over_det = 1.0f / (
        + M.elem[0][0] * (M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2])
        - M.elem[1][0] * (M.elem[0][1] * M.elem[2][2] - M.elem[2][1] * M.elem[0][2])
        + M.elem[2][0] * (M.elem[0][1] * M.elem[1][2] - M.elem[1][1] * M.elem[0][2]));

    md_float3x3_t I;

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

MD_VEC_INLINE float md_float3x3_determinant(md_float3x3_t M) {
    float d = M.elem[0][0] * (M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2])
            - M.elem[1][0] * (M.elem[0][1] * M.elem[2][2] - M.elem[2][1] * M.elem[0][2])
            + M.elem[2][0] * (M.elem[0][1] * M.elem[1][2] - M.elem[1][1] * M.elem[0][2]);
    return d;
}

MD_VEC_INLINE float md_float3x3_trace(md_float3x3_t M) {
    return M.elem[0][0] + M.elem[1][1] + M.elem[2][2];
}

MD_VEC_INLINE float md_float3x3_abs_trace(md_float3x3_t M) {
    return fabsf(M.elem[0][0]) + fabsf(M.elem[1][1]) + fabsf(M.elem[2][2]);
}

MD_VEC_INLINE float md_float3x3_abs_sum(md_float3x3_t M) {
    return
        fabsf(M.elem[0][0]) + fabsf(M.elem[0][1]) + fabsf(M.elem[0][2]) +
        fabsf(M.elem[1][0]) + fabsf(M.elem[1][1]) + fabsf(M.elem[1][2]) +
        fabsf(M.elem[2][0]) + fabsf(M.elem[2][1]) + fabsf(M.elem[2][2]);
}

MD_VEC_INLINE md_float3_t md_float3x3_diag(md_float3x3_t M) {
    md_float3_t v = {M.elem[0][0], M.elem[1][1], M.elem[2][2]};
    return v;
}

// Construct a rotation matrix from and angle and axis (Expects axis to be normalized)
MD_VEC_INLINE md_float3x3_t md_float3x3_angle_axis(float angle, md_float3_t axis) {
    float c = cosf(angle);
	float s = sinf(angle);
	float t = 1.0f - c;

	md_float3x3_t M = {
		t * axis.x * axis.x + c, t * axis.x * axis.y - s * axis.z, t * axis.x * axis.z + s * axis.y,
		t * axis.x * axis.y + s * axis.z, t * axis.y * axis.y + c, t * axis.y * axis.z - s * axis.x,
		t * axis.x * axis.z - s * axis.y, t * axis.y * axis.z + s * axis.x, t * axis.z * axis.z + c
	};
	return M;
}

// dvec2
MD_VEC_INLINE md_double2_t md_double2_set(double x, double y) {
    md_double2_t v = {x, y};
    return v;
}

MD_VEC_INLINE md_double2_t md_double2_set1(double s) {
    md_double2_t v = {s, s};
    return v;
}

MD_VEC_INLINE md_double2_t md_double2_load(const double* ptr) {
    md_double2_t v;
    md_memcpy(&v, ptr, sizeof(md_double2_t));
    return v;
}

MD_VEC_INLINE void md_double2_store(double* ptr, md_double2_t v) {
    md_memcpy(ptr, &v, sizeof(md_double2_t));
}

MD_VEC_INLINE double md_double2_dot(md_double2_t a, md_double2_t b) {
    return a.x * b.x + a.y * b.y;
}

MD_VEC_INLINE double md_double2_length(md_double2_t v) {
    return sqrt(md_double2_dot(v, v));
}

MD_VEC_INLINE md_double2_t md_double2_add(md_double2_t a, md_double2_t b) {
    md_double2_t r = {
        a.x + b.x,
        a.y + b.y,
    };
    return r;
}

// dvec3
MD_VEC_INLINE md_double3_t md_double3_set(double x, double y, double z) {
    md_double3_t v = {x, y, z};
    return v;
}

MD_VEC_INLINE md_double3_t md_double3_set1(double s) {
    md_double3_t v = {s, s, s};
    return v;
}

MD_VEC_INLINE md_double3_t md_double3_load(const double* ptr) {
    md_double3_t v;
    md_memcpy(&v, ptr, sizeof(md_double3_t));
    return v;
}

MD_VEC_INLINE void md_double3_store(double* ptr, md_double3_t v) {
    md_memcpy(ptr, &v, sizeof(md_double3_t));
}

MD_VEC_INLINE double md_double3_dot(md_double3_t a, md_double3_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

MD_VEC_INLINE double md_double3_length(md_double3_t v) {
    return sqrt(md_double3_dot(v, v));
}

MD_VEC_INLINE md_double3_t md_double3_add(md_double3_t a, md_double3_t b) {
    md_double3_t r = {
        a.x + b.x,
        a.y + b.y,
        a.z + b.z,
    };
    return r;
}

MD_VEC_INLINE md_double3_t md_double3_add1(md_double3_t v, double s) {
    md_double3_t r = {
        v.x + s,
        v.y + s,
        v.z + s,
    };
    return r;
}

MD_VEC_INLINE md_double3_t md_double3_sub(md_double3_t a, md_double3_t b) {
    md_double3_t r = {
        a.x - b.x,
        a.y - b.y,
        a.z - b.z,
    };
    return r;
}   

MD_VEC_INLINE md_double3_t md_double3_sub1(md_double3_t v, double s) {
    md_double3_t r = {
        v.x - s,
        v.y - s,
        v.z - s,
    };
    return r;
}

MD_VEC_INLINE md_double3_t md_double3_mul(md_double3_t a, md_double3_t b) {
    md_double3_t r = {
        a.x * b.x,
        a.y * b.y,
        a.z * b.z,
    };
    return r;
}

MD_VEC_INLINE md_double3_t md_double3_mul1(md_double3_t v, double s) {
    md_double3_t r = {
        v.x * s,
        v.y * s,
        v.z * s,
    };
    return r;
}

MD_VEC_INLINE md_double3_t md_double3_div(md_double3_t a, md_double3_t b) {
    md_double3_t r = {
        a.x / b.x,
        a.y / b.y,
        a.z / b.z,
    };
    return r;
}

MD_VEC_INLINE md_double3_t md_double3_div1(md_double3_t v, double s) {
    md_double3_t r = {
        v.x / s,
        v.y / s,
        v.z / s,
    };
    return r;
}

// dvec4
MD_VEC_INLINE md_double4_t md_double4_set(double x, double y, double z, double w) {
    md_double4_t v = {x, y, z, w};
    return v;
}

MD_VEC_INLINE md_double4_t md_double4_set1(double s) {
    md_double4_t v = {s, s, s, s};
    return v;
}

MD_VEC_INLINE md_double4_t md_double4_load(const double* ptr) {
    md_double4_t v;
    md_memcpy(&v, ptr, sizeof(md_double4_t));
    return v;
}

MD_VEC_INLINE void md_double4_store(double* ptr, md_double4_t v) {
    md_memcpy(ptr, &v, sizeof(md_double4_t));
}

MD_VEC_INLINE double md_double4_dot(md_double4_t a, md_double4_t b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

MD_VEC_INLINE double md_double4_length(md_double4_t v) {
    return sqrt(md_double4_dot(v, v));
}

// Math misc

typedef struct md_float3x3_eigen_t {
    md_float3x3_t vectors;
    md_float3_t values;
} md_float3x3_eigen_t;

typedef struct md_float3x3_svd_t {
    md_float3x3_t U;
    md_float3x3_t V;
    md_float3_t s;
} md_float3x3_svd_t;

md_float3x3_eigen_t md_float3x3_eigen(md_float3x3_t M);
md_float3x3_svd_t   md_float3x3_svd(md_float3x3_t M);
md_float3x3_t       md_float3x3_extract_rotation(md_float3x3_t M);
md_float3x3_t       md_float3x3_orthonormalize(md_float3x3_t M);

// Computes the covariance matrix for a set of coordinates with a given center of mass.
// xyz:         packed coordinates
// w:           weights (optional)
// indices:     indices into coordinates and weights (optional)
// count:       number of coordinates or indices
// mean:        mean (com if coordinates)
md_float3x3_t md_float3x3_covariance_matrix(const md_float3_t* xyz, const float* w, const int32_t* indices, size_t count, md_float3_t mean);
md_float3x3_t md_float3x3_covariance_matrix_float4(const md_float4_t* xyzw, const int32_t* indices, size_t count, md_float3_t mean);

// Computes the cross covariance matrix for two set of coordinates with given center of mass.
// The set of points are assumed to have equal length and if w is not NULL, the same weight.
// xyz[2]:               packed coordinate streams
// w[2]:                 weights (optional)
// indices[2]:           indices into coordinates and weights (optional)
// count:                number of coordinates or indices
// mean[2]:              mean (com if coordinates)
md_float3x3_t md_float3x3_cross_covariance_matrix(const md_float3_t* const xyz[2], const float* const w[2], const int32_t* const indices[2], size_t count, const md_float3_t mean[2]);

// Computes the cross covariance matrix for two set of relative coordinates (coordinates with center of mass subtracted).
md_float3x3_t md_float3x3_cross_covariance_matrix_raw_float4(const md_float4_t* const xyzw[2], size_t count);

// Compute the cross covariance matrix for two set of coordinates with given center of mass.
// xyzw[2]: coordinate + weights
// indices: indices into coordinates and weights (optional)
// count:   number of coordinates or indices
// mean[2]: mean (com if coordinates)
md_float3x3_t md_float3x3_cross_covariance_matrix_float4(const md_float4_t* const xyzw[2], const int32_t* const indices[2], size_t count, const md_float3_t mean[2]);

// Computes the optimal rotation matrix that minimizes the RMSD between two sets of coordinates.
// The set of points are assumed to have equal length and if w is not NULL, the same weight.
// xyz[2]:               packed coordinate streams
// w[2]:                 weights (optional)
// indices[2]:           indices into coordinates and weights (optional)
// count:                number of coordinates or indices
// com[2]:               center of mass
md_float3x3_t md_float3x3_optimal_rotation(const md_float3_t* const xyz[2], const float* const w[2], const int32_t* const indices[2], size_t count, const md_float3_t com[2]);

// Computes the optimal rotation matrix that minimizes the RMSD between two sets of coordinates.
// xyzw[2]: coordinate + weights
// indices: indices into coordinates and weights (optional)
// count:   number of coordinates or indices
// com[2]:  center of mass
md_float3x3_t md_float3x3_optimal_rotation_float4(const md_float4_t* const xyzw[2], const int32_t* const indices[2], size_t count, const md_float3_t com[2]);

// MAT4

MD_VEC_INLINE md_float4x4_t md_float4x4_ident(void) {
    md_float4x4_t M = {
        {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        }
    };
    return M;
}

MD_VEC_INLINE md_float4x4_t md_float4x4_load(const float* ptr) {
    md_float4x4_t M;
    md_memcpy(M.elem, ptr, sizeof(md_float4x4_t));
    return M;
}

MD_VEC_INLINE void md_float4x4_store(float* ptr, md_float4x4_t M) {
    md_memcpy(ptr, M.elem, sizeof(md_float4x4_t));
}

MD_VEC_INLINE md_float4x4_t md_float4x4_from_float3x3(const md_float3x3_t M) {
    md_float4x4_t R = {
        {
            M.col[0].x, M.col[0].y, M.col[0].z, 0,
            M.col[1].x, M.col[1].y, M.col[1].z, 0,
            M.col[2].x, M.col[2].y, M.col[2].z, 0,
            0, 0, 0, 1,
        }
    };
    return R;
}

static inline md_float4x4_t md_float4x4_from_quat(md_quat_t q) {
    float yy2 = 2.0f * q.elem[1] * q.elem[1];
    float xy2 = 2.0f * q.elem[0] * q.elem[1];
    float xz2 = 2.0f * q.elem[0] * q.elem[2];
    float yz2 = 2.0f * q.elem[1] * q.elem[2];
    float zz2 = 2.0f * q.elem[2] * q.elem[2];
    float wz2 = 2.0f * q.elem[3] * q.elem[2];
    float wy2 = 2.0f * q.elem[3] * q.elem[1];
    float wx2 = 2.0f * q.elem[3] * q.elem[0];
    float xx2 = 2.0f * q.elem[0] * q.elem[0];

    md_float4x4_t M;
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
MD_VEC_INLINE md_float4x4_t md_float4x4_scale(float x, float y, float z) {
    md_float4x4_t M = {
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
MD_VEC_INLINE md_float4x4_t md_float4x4_scale_f(float s) {
    md_float4x4_t M = {
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
MD_VEC_INLINE md_float4x4_t md_float4x4_scale_float3(md_float3_t v) {
    md_float4x4_t M = {
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
MD_VEC_INLINE md_float4x4_t md_float4x4_translate(float x, float y, float z) {
    md_float4x4_t M = {
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
MD_VEC_INLINE md_float4x4_t md_float4x4_translate_float3(md_float3_t t) {
    md_float4x4_t M = {
        {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            t.x,t.y,t.z,1
        }
    };
    return M;
}

MD_VEC_INLINE md_float4x4_t md_float4x4_add(md_float4x4_t A, md_float4x4_t B) {
    md_float4x4_t M;
    M.col[0] = md_float4_add(A.col[0], B.col[0]);
    M.col[1] = md_float4_add(A.col[1], B.col[1]);
    M.col[2] = md_float4_add(A.col[2], B.col[2]);
    M.col[3] = md_float4_add(A.col[3], B.col[3]);
    return M;
}

MD_VEC_INLINE md_float4x4_t md_float4x4_sub(md_float4x4_t A, md_float4x4_t B) {
    md_float4x4_t M;
    M.col[0] = md_float4_sub(A.col[0], B.col[0]);
    M.col[1] = md_float4_sub(A.col[1], B.col[1]);
    M.col[2] = md_float4_sub(A.col[2], B.col[2]);
    M.col[3] = md_float4_sub(A.col[3], B.col[3]);
    return M;
}

MD_VEC_INLINE md_float4x4_t md_float4x4_mul(md_float4x4_t A, md_float4x4_t B) {
    md_float4x4_t C;
    md_mm_store_ps(C.col[0].elem, md_float4__linear_combine_4(md_mm_load_ps(B.col[0].elem), A.col));
    md_mm_store_ps(C.col[1].elem, md_float4__linear_combine_4(md_mm_load_ps(B.col[1].elem), A.col));
    md_mm_store_ps(C.col[2].elem, md_float4__linear_combine_4(md_mm_load_ps(B.col[2].elem), A.col));
    md_mm_store_ps(C.col[3].elem, md_float4__linear_combine_4(md_mm_load_ps(B.col[3].elem), A.col));
    return C;
}

MD_VEC_INLINE md_float4_t md_float4x4_mul_float4(md_float4x4_t M, md_float4_t v) {
    md_float4_t r;
    md_mm_store_ps(r.elem, md_float4__linear_combine_4(md_mm_load_ps(v.elem), M.col));
    return r;
}

MD_VEC_INLINE md_float3_t md_float4x4_mul_float3(md_float4x4_t M, md_float3_t v, float w) {
    md_float4_t r = {v.x, v.y, v.z, w};
    r = md_float4x4_mul_float4(M, r);
    return md_float3_from_float4(r);
}

MD_VEC_INLINE md_float4x4_t md_float4x4_mul1(md_float4x4_t M, float s) {
    md_float4x4_t C = {0};
    C.col[0] = md_float4_mul1(M.col[0], s);
    C.col[1] = md_float4_mul1(M.col[1], s);
    C.col[2] = md_float4_mul1(M.col[2], s);
    C.col[3] = md_float4_mul1(M.col[3], s);
    return C;
}

md_float4x4_t md_float4x4_inverse(md_float4x4_t M);

MD_VEC_INLINE md_float4x4_t md_float4x4_transpose(md_float4x4_t M) {
    md_float4x4_t T;

    md_128 c0 = md_mm_load_ps(M.col[0].elem), c1 = md_mm_load_ps(M.col[1].elem), c2 = md_mm_load_ps(M.col[2].elem), c3 = md_mm_load_ps(M.col[3].elem);
    SIMDE_MM_TRANSPOSE4_PS(c0, c1, c2, c3);
    md_mm_store_ps(T.col[0].elem, c0);
    md_mm_store_ps(T.col[1].elem, c1);
    md_mm_store_ps(T.col[2].elem, c2);
    md_mm_store_ps(T.col[3].elem, c3);
    return T;
}

MD_VEC_INLINE bool md_float4x4_equal(md_float4x4_t A, md_float4x4_t B) {
    return md_float4_equal(A.col[0], B.col[0]) && md_float4_equal(A.col[1], B.col[1]) && md_float4_equal(A.col[2], B.col[2]) && md_float4_equal(A.col[3], B.col[3]);
}

// md_float4x3_t: a 3x3 matrix with its columns padded to four, for SIMD. Its product with a md_float4_t uses x, y, z
// (the three components a 4x3 matrix takes) and leaves w out.

MD_VEC_INLINE md_float4x3_t md_float4x3_from_float3x3(md_float3x3_t M) {
    md_float4x3_t R;
    R.col[0] = md_float4_from_float3(M.col[0], 0);
    R.col[1] = md_float4_from_float3(M.col[1], 0);
    R.col[2] = md_float4_from_float3(M.col[2], 0);
    return R;
}

MD_VEC_INLINE md_float4_t md_float4x3_mul_float4(md_float4x3_t M, md_float4_t v) {
    md_float4_t r = {0};
    md_mm_store_ps(r.elem, md_float4__linear_combine_3(md_mm_load_ps(v.elem), M.col));
    return r;
}


// More complex functionality

md_float3_t md_float4x4_unproject(md_float3_t window_coords, md_float4x4_t inv_view_proj_mat, md_float4_t viewport);

md_float4x4_t md_float4x4_look_at(md_float3_t eye, md_float3_t center, md_float3_t up);

// 'near'/'far' are off limits as parameter names here: windef.h defines both away to nothing
// (see the note by md_float4x4_persp in md_vec_math.c).
md_float4x4_t md_float4x4_ortho(float left, float right, float bottom, float top, float z_near, float z_far);
md_float4x4_t md_float4x4_ortho_inv(float left, float right, float bottom, float top, float z_near, float z_far);
md_float4x4_t md_float4x4_ortho_2d(float left, float right, float bottom, float top);
md_float4x4_t md_float4x4_ortho_2d_inv(float left, float right, float bottom, float top);

md_float4x4_t md_float4x4_persp(float fovy, float aspect, float z_near, float z_far);
md_float4x4_t md_float4x4_persp_inv(float fovy, float aspect, float z_near, float z_far);

md_float4x4_t md_float4x4_frustum(float left, float right, float bottom, float top, float z_near, float z_far);
md_float4x4_t md_float4x4_frustum_inv(float left, float right, float bottom, float top, float z_near, float z_far);

// These are routines for performing the same operation on many items, packed coordinates in and
// out. The out-of-place forms may be given the same array for in and out.
void md_float3_batch_translate_inplace(md_float3_t* xyz, size_t count, md_float3_t translation);
void md_float3_batch_translate(md_float3_t* out_xyz, const md_float3_t* in_xyz, size_t count, md_float3_t translation);

void md_float3x3_batch_transform_inplace(md_float3_t* xyz, size_t count, md_float3x3_t transform);
void md_float3x3_batch_transform(md_float3_t* out_xyz, const md_float3_t* in_xyz, size_t count, md_float3x3_t transform);

void md_float4x4_batch_transform_inplace(md_float3_t* xyz, float w_comp, size_t count, md_float4x4_t transform);
void md_float4x4_batch_transform(md_float3_t* out_xyz, const md_float3_t* in_xyz, float w_comp, size_t count, md_float4x4_t transform);

#ifdef __cplusplus
}

static inline bool operator == (md_float2_t a, md_float2_t b){
    return a.x == b.x && a.y == b.y;
}

static inline bool operator != (md_float2_t a, md_float2_t b){
    return !(a == b);
}

static inline md_float2_t operator + (md_float2_t a, md_float2_t b) {
    md_float2_t c = {a.x + b.x, a.y + b.y};
    return c;
}

static inline md_float2_t operator + (md_float2_t a, float s) {
    md_float2_t c = {a.x + s, a.y + s};
    return c;
}

static inline md_float2_t operator + (float s, md_float2_t a) {
    md_float2_t c = {a.x + s, a.y + s};
    return c;
}

static inline md_float2_t operator - (md_float2_t a, md_float2_t b) {
    md_float2_t c = {a.x - b.x, a.y - b.y};
    return c;
}

static inline md_float2_t operator - (md_float2_t a, float s) {
    md_float2_t c = {a.x - s, a.y - s};
    return c;
}

static inline md_float2_t operator - (const md_float2_t &a) {
    md_float2_t c = {-a.x, -a.y};
    return c;
}

static inline md_float2_t operator - (float s, md_float2_t a) {
    md_float2_t c = {s - a.x, s - a.y};
    return c;
}

static inline md_float2_t operator * (md_float2_t a, md_float2_t b) {
    md_float2_t c = {a.x * b.x, a.y * b.y};
    return c;
}

static inline md_float2_t operator * (md_float2_t a, float s) {
    md_float2_t c = {a.x * s, a.y * s};
    return c;
}

static inline md_float2_t operator / (md_float2_t a, md_float2_t b) {
    md_float2_t c = {a.x / b.x, a.y / b.y};
    return c;
}

static inline md_float2_t operator / (md_float2_t a, float s) {
    md_float2_t c = {a.x / s, a.y / s};
    return c;
}

static inline md_float2_t operator / (float s, md_float2_t a) {
    md_float2_t c = {s / a.x, s / a.y};
    return c;
}

static inline md_float2_t operator *= (md_float2_t& a, md_float2_t b) {
    a = md_float2_mul(a, b);
    return a;
}

static inline md_float2_t operator *= (md_float2_t& v, float s) {
    v = md_float2_mul1(v, s);
    return v;
}

static inline md_float2_t operator /= (md_float2_t& a, md_float2_t b) {
    a = md_float2_div(a, b);
    return a;
}

static inline md_float2_t operator /= (md_float2_t& v, float s) {
    v = md_float2_div1(v, s);
    return v;
}

static inline md_float2_t operator += (md_float2_t& a, md_float2_t b) {
    a = md_float2_add(a, b);
    return a;
}

static inline md_float2_t operator += (md_float2_t& v, float s) {
    v = md_float2_add1(v, s);
    return v;
}

static inline md_float2_t operator -= (md_float2_t& a, md_float2_t b) {
    a = md_float2_sub(a, b);
    return a;
}

static inline md_float2_t operator -= (md_float2_t& v, float s) {
    v = md_float2_sub1(v, s);
    return v;
}

// vec3
static inline bool operator == (md_float3_t a, md_float3_t b){
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static inline bool operator != (md_float3_t a, md_float3_t b){
    return !(a == b);
}

static inline md_float3_t operator + (md_float3_t a, md_float3_t b) {
    return md_float3_add(a, b);
}

static inline md_float3_t operator + (md_float3_t v, float s) {
    return md_float3_add1(v, s);
}

static inline md_float3_t operator + (float s, md_float3_t v) {
    return md_float3_add1(v, s);
}

static inline md_float3_t operator += (md_float3_t& a, md_float3_t b) {
    a = md_float3_add(a, b);
    return a;
}

static inline md_float3_t operator += (md_float3_t& v, float s) {
    v = md_float3_add1(v, s);
    return v;
}

static inline md_float3_t operator - (md_float3_t a, md_float3_t b) {
    md_float3_t c = {a.x - b.x, a.y - b.y, a.z - b.z};
    return c;
}

static inline md_float3_t operator - (md_float3_t a, float s) {
    md_float3_t c = {a.x - s, a.y - s, a.z - s};
    return c;
}

static inline md_float3_t operator - (float s, md_float3_t a) {
    md_float3_t c = {s - a.x, s - a.y, s - a.z};
    return c;
}

static inline md_float3_t operator - (const md_float3_t &a) {
    md_float3_t c = {-a.x, -a.y, -a.z};
    return c;
}

static inline md_float3_t operator -= (md_float3_t& a, md_float3_t b) {
    a = md_float3_sub(a, b);
    return a;
}

static inline md_float3_t operator -= (md_float3_t& v, float s) {
    v = md_float3_sub1(v, s);
    return v;
}

static inline md_float3_t operator * (md_float3_t a, md_float3_t b) {
    md_float3_t c = {a.x * b.x, a.y * b.y, a.z * b.z};
    return c;
}

static inline md_float3_t operator * (md_float3_t a, float s) {
    md_float3_t c = {a.x * s, a.y * s, a.z * s};
    return c;
}

static inline md_float3_t operator * (float s, md_float3_t a) {
    md_float3_t c = {a.x * s, a.y * s, a.z * s};
    return c;
}

static inline md_float3_t operator *= (md_float3_t& a, md_float3_t b) {
    a = md_float3_mul(a, b);
    return a;
}

static inline md_float3_t operator *= (md_float3_t& v, float s) {
    v = md_float3_mul1(v, s);
    return v;
}

static inline md_float3_t operator / (md_float3_t a, md_float3_t b) {
    md_float3_t c = {a.x / b.x, a.y / b.y, a.z / b.z};
    return c;
}

static inline md_float3_t operator / (md_float3_t a, float s) {
    md_float3_t c = {a.x / s, a.y / s, a.z / s};
    return c;
}

static inline md_float3_t operator / (float s, md_float3_t a) {
    md_float3_t c = {s / a.x, s / a.y, s / a.z};
    return c;
}

static inline md_float3_t operator /= (md_float3_t& a, md_float3_t b) {
    a = md_float3_div(a, b);
    return a;
}

static inline md_float3_t operator /= (md_float3_t& v, float s) {
    v = md_float3_div1(v, s);
    return v;
}

// dvec3
static inline bool operator == (md_double3_t a, md_double3_t b){
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

static inline bool operator != (md_double3_t a, md_double3_t b){
    return !(a == b);
}

static inline md_double3_t operator + (md_double3_t a, md_double3_t b) {
    return md_double3_add(a, b);
}

static inline md_double3_t operator + (md_double3_t v, double s) {
    return md_double3_add1(v, s);
}

static inline md_double3_t operator + (double s, md_double3_t v) {
    return md_double3_add1(v, s);
}

static inline md_double3_t operator += (md_double3_t& a, md_double3_t b) {
    a = md_double3_add(a, b);
    return a;
}

static inline md_double3_t operator += (md_double3_t& v, double s) {
    v = md_double3_add1(v, s);
    return v;
}

static inline md_double3_t operator - (md_double3_t a, md_double3_t b) {
    md_double3_t c = {a.x - b.x, a.y - b.y, a.z - b.z};
    return c;
}

static inline md_double3_t operator - (md_double3_t a, double s) {
    md_double3_t c = {a.x - s, a.y - s, a.z - s};
    return c;
}

static inline md_double3_t operator - (double s, md_double3_t a) {
    md_double3_t c = {s - a.x, s - a.y, s - a.z};
    return c;
}

static inline md_double3_t operator - (const md_double3_t &a) {
    md_double3_t c = {-a.x, -a.y, -a.z};
    return c;
}

static inline md_double3_t operator -= (md_double3_t& a, md_double3_t b) {
    a = md_double3_sub(a, b);
    return a;
}

static inline md_double3_t operator -= (md_double3_t& v, double s) {
    v = md_double3_sub1(v, s);
    return v;
}

static inline md_double3_t operator * (md_double3_t a, md_double3_t b) {
    md_double3_t c = {a.x * b.x, a.y * b.y, a.z * b.z};
    return c;
}

static inline md_double3_t operator * (md_double3_t a, double s) {
    md_double3_t c = {a.x * s, a.y * s, a.z * s};
    return c;
}

static inline md_double3_t operator * (double s, md_double3_t a) {
    md_double3_t c = {a.x * s, a.y * s, a.z * s};
    return c;
}

static inline md_double3_t operator *= (md_double3_t& a, md_double3_t b) {
    a = md_double3_mul(a, b);
    return a;
}

static inline md_double3_t operator *= (md_double3_t& v, double s) {
    v = md_double3_mul1(v, s);
    return v;
}

static inline md_double3_t operator / (md_double3_t a, md_double3_t b) {
    md_double3_t c = {a.x / b.x, a.y / b.y, a.z / b.z};
    return c;
}

static inline md_double3_t operator / (md_double3_t a, double s) {
    md_double3_t c = {a.x / s, a.y / s, a.z / s};
    return c;
}

static inline md_double3_t operator / (double s, md_double3_t a) {
    md_double3_t c = {s / a.x, s / a.y, s / a.z};
    return c;
}

static inline md_double3_t operator /= (md_double3_t& a, md_double3_t b) {
    a = md_double3_div(a, b);
    return a;
}

static inline md_double3_t operator /= (md_double3_t& v, double s) {
    v = md_double3_div1(v, s);
    return v;
}

// md_float4_t
static inline bool operator == (md_float4_t a, md_float4_t b){
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

static inline bool operator != (md_float4_t a, md_float4_t b){
    return !(a == b);
}

static inline md_float4_t& operator += (md_float4_t& a, md_float4_t b) {
    a = md_float4_add(a,b);
    return a;
}

static inline md_float4_t& operator += (md_float4_t& v, float s) {
    v = md_float4_add1(v,s);
    return v;
}

static inline md_float4_t operator + (md_float4_t a, md_float4_t b) {
    return md_float4_add(a, b);
}

static inline md_float4_t operator + (md_float4_t v, float s) {
    return md_float4_add1(v, s);
}

static inline md_float4_t operator + (float s, md_float4_t v) {
    return md_float4_add1(v, s);
}

static inline md_float4_t& operator -= (md_float4_t& a, md_float4_t b) {
    a = md_float4_sub(a,b);
    return a;
}

static inline md_float4_t& operator -= (md_float4_t& v, float s) {
    v = md_float4_sub1(v,s);
    return v;
}

static inline md_float4_t operator - (md_float4_t a, md_float4_t b) {
    return md_float4_sub(a, b);
}

static inline md_float4_t operator - (md_float4_t v, float s) {
    return md_float4_sub1(v, s);
}

static inline md_float4_t operator - (float s, md_float4_t v) {
    return md_float4_sub(md_float4_from_float(s), v);
}

static inline md_float4_t operator - (const md_float4_t &a) {
    md_float4_t c = {-a.x, -a.y, -a.z, -a.w};
    return c;
}

static inline md_float4_t& operator *= (md_float4_t& a, md_float4_t b) {
    a = md_float4_mul(a,b);
    return a;
}

static inline md_float4_t& operator *= (md_float4_t& v, float s) {
    v = md_float4_mul1(v,s);
    return v;
}

static inline md_float4_t operator * (md_float4_t a, md_float4_t b) {
    return md_float4_mul(a, b);
}

static inline md_float4_t operator * (md_float4_t v, float s) {
    return md_float4_mul1(v, s);
}

static inline md_float4_t operator * (float s, md_float4_t v) {
    return md_float4_mul1(v, s);
}

static inline md_float4_t& operator /= (md_float4_t& a, md_float4_t b) {
    a = md_float4_div(a,b);
    return a;
}

static inline md_float4_t& operator /= (md_float4_t& v, float s) {
    v = md_float4_div1(v,s);
    return v;
}

static inline md_float4_t operator / (md_float4_t a, md_float4_t b) {
    return md_float4_div(a, b);
}

static inline md_float4_t operator / (md_float4_t v, float s) {
    return md_float4_div1(v, s);
}

static inline md_float4_t operator / (float s, md_float4_t v) {
    return md_float4_div(md_float4_from_float(s), v);
}

// quat

static inline md_quat_t operator * (md_quat_t a, md_quat_t b) {
    return md_quat_mul(a, b);
}

static inline md_float3_t operator * (md_quat_t q, md_float3_t v) {
    return md_quat_mul_float3(q, v);
}

// mat

static inline md_float3x3_t operator * (md_float3x3_t A, md_float3x3_t B) {
    return md_float3x3_mul(A, B);
}

static inline md_float3x3_t operator * (md_float3x3_t A, float s) {
    return md_float3x3_mul1(A, s);
}

static inline md_float3x3_t operator * (float s, md_float3x3_t A) {
    return md_float3x3_mul1(A, s);
}

static inline md_float3x3_t operator + (md_float3x3_t A, md_float3x3_t B) {
    return md_float3x3_add(A, B);
}

static inline md_float3x3_t operator - (md_float3x3_t A, md_float3x3_t B) {
    return md_float3x3_sub(A, B);
}

static inline bool operator == (md_float3x3_t A, md_float3x3_t B) {
    return A[0] == B[0] && A[1] == B[1] && A[2] == B[2];
}

static inline bool operator != (md_float3x3_t A, md_float3x3_t B) {
    return !(A == B);
}

static inline md_float3_t operator * (md_float3x3_t M, md_float3_t v) {
    return md_float3x3_mul_float3(M, v);
}

static inline md_float4x4_t operator * (md_float4x4_t A, md_float4x4_t B) {
    return md_float4x4_mul(A, B);
}

static inline md_float4x4_t operator + (md_float4x4_t A, md_float4x4_t B) {
    return md_float4x4_add(A, B);
}

static inline md_float4x4_t operator - (md_float4x4_t A, md_float4x4_t B) {
    return md_float4x4_sub(A, B);
}

static inline md_float4_t operator * (md_float4x4_t M, md_float4_t v) {
    return md_float4x4_mul_float4(M, v);
}

static inline bool operator == (md_float4x4_t A, md_float4x4_t B) {
    return A[0] == B[0] && A[1] == B[1] && A[2] == B[2] && A[3] == B[3];
}

static inline bool operator != (md_float4x4_t A, md_float4x4_t B) {
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
