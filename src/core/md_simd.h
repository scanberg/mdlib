#ifndef _MD_SIMD_H_
#define _MD_SIMD_H_

#include "md_intrinsics.h"

#ifdef __AVX__
    #define md_simd_width 8
    #define md_simd_typef __m256
    #define md_simd_loadf  md_simd_load_f256
    #define md_simd_storef md_simd_store_f256

    #define md_simd_set1f md_simd_set1_f256
    #define md_simd_setf  md_simd_set_f256
    #define md_simd_zerof md_simd_zero_f256

    #define md_simd_addf md_simd_add_f256
    #define md_simd_subf md_simd_sub_f256
    #define md_simd_mulf md_simd_mul_f256
    #define md_simd_divf md_simd_div_f256

    #define md_simd_andf     md_simd_and_f256
    #define md_simd_and_notf md_simd_and_not_f256
    #define md_simd_orf      md_simd_or_f256
    #define md_simd_xorf     md_simd_xor_f256

    #define md_simd_cmp_gtf  md_simd_cmp_gt_f256
    #define md_simd_cmp_gef  md_simd_cmp_ge_f256
    #define md_simd_cmp_ltf  md_simd_cmp_lt_f256
    #define md_simd_cmp_lef  md_simd_cmp_le_f256
    #define md_simd_cmp_eqf  md_simd_cmp_eq_f256
    #define md_simd_cmp_neqf md_simd_cmp_neq_f256

    #define md_simd_absf    md_simd_abs_f256
    #define md_simd_minf    md_simd_min_f256
    #define md_simd_maxf    md_simd_max_f256
    #define md_simd_signf   md_simd_sign_f256

    #define md_simd_horizontal_minf md_simd_horizontal_min_f256
    #define md_simd_horizontal_maxf md_simd_horizontal_max_f256
    #define md_simd_horizontal_addf md_simd_horizontal_add_f256

    #define md_simd_stepf         md_simd_step_f256
    #define md_simd_lerpf         md_simd_lerp_f256
    #define md_simd_cubic_splinef md_simd_cubic_spline_f256
#else
    #define md_simd_width 4
    #define md_simd_typef __m128
    #define md_simd_loadf  md_simd_load_f128
    #define md_simd_storef md_simd_store_f128

    #define md_simd_set1f md_simd_set1_f128
    #define md_simd_setf  md_simd_set_f128
    #define md_simd_zerof md_simd_zero_f128

    #define md_simd_addf md_simd_add_f128
    #define md_simd_subf md_simd_sub_f128
    #define md_simd_mulf md_simd_mul_f128
    #define md_simd_divf md_simd_div_f128

    #define md_simd_andf     md_simd_and_f128
    #define md_simd_and_notf md_simd_and_not_f128
    #define md_simd_orf      md_simd_or_f128
    #define md_simd_xorf     md_simd_xor_f128

    #define md_simd_cmp_gtf  md_simd_cmp_gt_f128
    #define md_simd_cmp_gef  md_simd_cmp_ge_f128
    #define md_simd_cmp_ltf  md_simd_cmp_lt_f128
    #define md_simd_cmp_lef  md_simd_cmp_le_f128
    #define md_simd_cmp_eqf  md_simd_cmp_eq_f128
    #define md_simd_cmp_neqf md_simd_cmp_neq_f128

    #define md_simd_absf    md_simd_abs_f128
    #define md_simd_minf    md_simd_min_f128
    #define md_simd_maxf    md_simd_max_f128
    #define md_simd_signf   md_simd_sign_f128

    #define md_simd_horizontal_minf md_simd_horizontal_min_f128
    #define md_simd_horizontal_maxf md_simd_horizontal_max_f128
    #define md_simd_horizontal_addf md_simd_horizontal_add_f128

    #define md_simd_stepf         md_simd_step_f128
    #define md_simd_lerpf         md_simd_lerp_f128
    #define md_simd_cubic_splinef md_simd_cubic_spline_f128
#endif

// 128-bit wide
// float operations

static inline __m128 md_simd_set1_f128(float x) { return _mm_set1_ps(x); }
static inline __m128 md_simd_set_f128(float x, float y, float z, float w) { return _mm_set_ps(w, z, y, x); }

static inline __m128 md_simd_zero_f128() { return _mm_setzero_ps(); }

static inline __m128 md_simd_load_f128(const float* addr) { return _mm_loadu_ps(addr); }
static inline __m128 md_simd_load_aligned_f128(const float* addr) {
    ASSERT(IS_ALIGNED(addr, 16));
    return _mm_load_ps(addr);
}

static inline void md_simd_store_f128(float* addr, __m128 v) { _mm_storeu_ps(addr, v); }
static inline void md_simd_store_aligned_f128(float* addr, __m128 v) {
    ASSERT(IS_ALIGNED(addr, 16));
    _mm_store_ps(addr, v);
}

static inline __m128 md_simd_add_f128(__m128 a, __m128 b) { return _mm_add_ps(a, b); }
static inline __m128 md_simd_sub_f128(__m128 a, __m128 b) { return _mm_sub_ps(a, b); }
static inline __m128 md_simd_mul_f128(__m128 a, __m128 b) { return _mm_mul_ps(a, b); }
static inline __m128 md_simd_div_f128(__m128 a, __m128 b) { return _mm_div_ps(a, b); }

static inline __m128 md_simd_and_f128    (__m128 a, __m128 b) { return _mm_and_ps(a, b); }
static inline __m128 md_simd_and_not_f128(__m128 a, __m128 b) { return _mm_andnot_ps(a, b); }
static inline __m128 md_simd_or_f128     (__m128 a, __m128 b) { return _mm_or_ps(a, b); }
static inline __m128 md_simd_xor_f128    (__m128 a, __m128 b) { return _mm_xor_ps(a, b); }

static inline __m128 md_simd_cmp_gt_f128 (__m128 a, __m128 b) { return _mm_cmpgt_ps(a, b); }
static inline __m128 md_simd_cmp_ge_f128 (__m128 a, __m128 b) { return _mm_cmpge_ps(a, b); }
static inline __m128 md_simd_cmp_lt_f128 (__m128 a, __m128 b) { return _mm_cmplt_ps(a, b); }
static inline __m128 md_simd_cmp_le_f128 (__m128 a, __m128 b) { return _mm_cmple_ps(a, b); }
static inline __m128 md_simd_cmp_eq_f128 (__m128 a, __m128 b) { return _mm_cmpeq_ps(a, b); }
static inline __m128 md_simd_cmp_neq_f128(__m128 a, __m128 b) { return _mm_cmpneq_ps(a, b); }

static inline __m128 md_simd_abs_f128(__m128 a) { return _mm_and_ps(a, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
static inline __m128 md_simd_sign_f128(__m128 x) {
    __m128 sgn = _mm_and_ps(x, _mm_set1_ps(-0.0f));
    __m128 res = _mm_xor_ps(sgn, _mm_set1_ps(1.0f));
    return res;
}

static inline __m128 md_simd_min_f128(__m128 a, __m128 b) { return _mm_min_ps(a, b); }
static inline __m128 md_simd_max_f128(__m128 a, __m128 b) { return _mm_max_ps(a, b); }

static inline float md_simd_horizontal_min_f128(__m128 x) {
    __m128 min1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
    __m128 min2 = _mm_min_ps(x, min1);
    __m128 min3 = _mm_shuffle_ps(min2, min2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128 min4 = _mm_min_ps(min2, min3);
    return _mm_cvtss_f32(min4);
}

static inline float md_simd_horizontal_max_f128(__m128 x) {
    __m128 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
    __m128 max2 = _mm_max_ps(x, max1);
    __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128 max4 = _mm_max_ps(max2, max3);
    return _mm_cvtss_f32(max4);
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
static inline float md_simd_horizontal_add_f128(__m128 x) {
    __m128 shuf = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));
    __m128 sums = _mm_add_ps(x, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

static inline __m128 md_simd_step_f128(__m128 edge, __m128 x) {
    __m128 cmp = _mm_cmpge_ps(x, edge);
    return _mm_and_ps(cmp, _mm_set1_ps(1.f));
}

static inline __m128 md_simd_lerp_f128(__m128 a, __m128 b, float t) {
    __m128 one_minus_alpha = _mm_set1_ps(1.0f - t);
    __m128 alpha = _mm_set1_ps(t);
    return _mm_add_ps(_mm_mul_ps(a, one_minus_alpha), _mm_mul_ps(b, alpha));
}

static inline __m128 md_simd_cubic_spline_f128(__m128 p0, __m128 p1, __m128 p2, __m128 p3, float s, float tension) {
    __m128 vt = _mm_set1_ps(tension);
    __m128 s1 = _mm_set1_ps(s);
    __m128 s2 = _mm_mul_ps(s1, s1);
    __m128 s3 = _mm_mul_ps(s2, s1);
    __m128 v0 = _mm_mul_ps(_mm_sub_ps(p2, p0), vt);
    __m128 v1 = _mm_mul_ps(_mm_sub_ps(p3, p1), vt);
    __m128 x0 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(2), _mm_sub_ps(p1, p2)), _mm_add_ps(v0, v1));
    __m128 x1 = _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(3), _mm_sub_ps(p2, p1)), _mm_add_ps(_mm_mul_ps(_mm_set1_ps(2), v0), v1));
    __m128 r0 = _mm_add_ps(_mm_mul_ps(x0, s3), _mm_mul_ps(x1, s2));
    __m128 r1 = _mm_add_ps(_mm_mul_ps(v0, s1), p1);
    return _mm_add_ps(r0, r1);
}

// int operations (fill this in when the need comes)
static inline __m128i md_simd_set1_i128(int x) { return _mm_set1_epi32(x); }
static inline __m128i md_simd_set_i128(int x, int y, int z, int w) { return _mm_set_epi32(w, z, y, x); }

static inline __m128i md_simd_zero_i128() { return _mm_setzero_si128(); }

// 256-bit wide
// float operations


static inline __m256 md_simd_set1_f256(float x) { return _mm256_set1_ps(x); }
static inline __m256 md_simd_set_f256(float x0, float y0, float z0, float w0, float x1, float y1, float z1, float w1) { return _mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0); }

static inline __m256 md_simd_zero_f256() { return _mm256_setzero_ps(); }

static inline __m256 md_simd_load_f256(const float* addr) { return _mm256_loadu_ps(addr); }
static inline __m256 md_simd_load_aligned_f256(const float* addr) {
    ASSERT(IS_ALIGNED(addr, 16));
    return _mm256_load_ps(addr);
}

static inline void md_simd_store_f256(float* addr, __m256 v) { _mm256_storeu_ps(addr, v); }
static inline void md_simd_store_aligned_f256(float* addr, __m256 v) {
    ASSERT(IS_ALIGNED(addr, 16));
    _mm256_store_ps(addr, v);
}

static inline __m256 md_simd_add_f256(__m256 a, __m256 b) { return _mm256_add_ps(a, b); }
static inline __m256 md_simd_sub_f256(__m256 a, __m256 b) { return _mm256_sub_ps(a, b); }
static inline __m256 md_simd_mul_f256(__m256 a, __m256 b) { return _mm256_mul_ps(a, b); }
static inline __m256 md_simd_div_f256(__m256 a, __m256 b) { return _mm256_div_ps(a, b); }

static inline __m256 md_simd_and_f256    (__m256 a, __m256 b) { return _mm256_and_ps(a, b); }
static inline __m256 md_simd_and_not_f256(__m256 a, __m256 b) { return _mm256_andnot_ps(a, b); }
static inline __m256 md_simd_or_f256     (__m256 a, __m256 b) { return _mm256_or_ps(a, b); }
static inline __m256 md_simd_xor_f256    (__m256 a, __m256 b) { return _mm256_xor_ps(a, b); }

static inline __m256 md_simd_cmp_gt_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_GT_OQ); }
static inline __m256 md_simd_cmp_ge_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_GE_OQ); }
static inline __m256 md_simd_cmp_lt_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }
static inline __m256 md_simd_cmp_le_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_LE_OQ); }
static inline __m256 md_simd_cmp_eq_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }
static inline __m256 md_simd_cmp_neq_f256(__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ); }

static inline __m256 md_simd_abs_f256(__m256 a) { return _mm256_and_ps(a, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
static inline __m256 md_simd_sign_f256(__m256 x) {
    __m256 sgn = _mm256_and_ps(x, _mm256_set1_ps(-0.0f));
    __m256 res = _mm256_xor_ps(sgn, _mm256_set1_ps(1.0f));
    return res;
}

static inline __m256 md_simd_min_f256(__m256 a, __m256 b) { return _mm256_min_ps(a, b); }
static inline __m256 md_simd_max_f256(__m256 a, __m256 b) { return _mm256_max_ps(a, b); }

static inline float md_simd_horizontal_min_f256(__m256 x) {
    __m128 lo = _mm256_castps256_ps128(x);
    __m128 hi = _mm256_extractf128_ps(x, 0x1);
    __m128 min_val = _mm_min_ps(lo, hi);
    return md_simd_horizontal_min_f128(min_val);
}

static inline float md_simd_horizontal_max_f256(__m256 x) {
    __m128 lo = _mm256_castps256_ps128(x);
    __m128 hi = _mm256_extractf128_ps(x, 0x1);
    __m128 max_val = _mm_max_ps(lo, hi);
    return md_simd_horizontal_max_f128(max_val);
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
static inline float md_simd_horizontal_add_f256(__m256 x) {
    __m128 vlow = _mm256_castps256_ps128(x);
    __m128 vhigh = _mm256_extractf128_ps(x, 1); // high 128
    vlow = _mm_add_ps(vlow, vhigh);     // add the low 128
    __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(vlow, shuf);
    shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

static inline __m256 md_simd_step_f256(__m256 edge, __m256 x) {
    __m256 cmp = _mm256_cmp_ps(x, edge, _CMP_GE_OQ);
    return _mm256_and_ps(cmp, _mm256_set1_ps(1.f));
}

static inline __m256 md_simd_lerp_f256(__m256 a, __m256 b, float t) {
    __m256 one_minus_alpha = _mm256_set1_ps(1.0f - t);
    __m256 alpha = _mm256_set1_ps(t);
    return _mm256_add_ps(_mm256_mul_ps(a, one_minus_alpha), _mm256_mul_ps(b, alpha));
}

static inline __m256 md_simd_cubic_spline_f256(__m256 p0, __m256 p1, __m256 p2, __m256 p3, float s, float tension) {
    __m256 vt = _mm256_set1_ps(tension);
    __m256 s1 = _mm256_set1_ps(s);
    __m256 s2 = _mm256_mul_ps(s1, s1);
    __m256 s3 = _mm256_mul_ps(s2, s1);
    __m256 v0 = _mm256_mul_ps(_mm256_sub_ps(p2, p0), vt);
    __m256 v1 = _mm256_mul_ps(_mm256_sub_ps(p3, p1), vt);
    __m256 x0 = _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(2), _mm256_sub_ps(p1, p2)), _mm256_add_ps(v0, v1));
    __m256 x1 = _mm256_sub_ps(_mm256_mul_ps(_mm256_set1_ps(3), _mm256_sub_ps(p2, p1)), _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(2), v0), v1));
    __m256 r0 = _mm256_add_ps(_mm256_mul_ps(x0, s3), _mm256_mul_ps(x1, s2));
    __m256 r1 = _mm256_add_ps(_mm256_mul_ps(v0, s1), p1);
    return _mm256_add_ps(r0, r1);
}

#endif
