#pragma once

#include "md_intrinsics.h"

#define md_simd_f128_t __m128
#define md_simd_f256_t __m256
#define md_simd_i128_t __m128i
#define md_simd_i256_t __m256i

#define MD_SIMD_INLINE static FORCE_INLINE

#ifdef __AVX__
// Float
    #define md_simd_widthf 8
    #define md_simd_typef  md_simd_f256_t
    #define md_simd_loadf  md_simd_load_f256
    #define md_simd_storef md_simd_store_f256

    #define md_simd_set1f md_simd_set1_f256
    #define md_simd_setf  md_simd_set_f256
    #define md_simd_zerof md_simd_zero_f256

    #define MD_SIMD_CONSTF 

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
    #define md_simd_signcopyf   md_simd_signcopy_f256
    #define md_simd_fractf  md_simd_fract_f256
    #define md_simd_roundf  md_simd_round_f256
    #define md_simd_blendf  md_simd_blend_f256
    #define md_simd_sqrtf   md_simd_sqrtf_f256

    #define md_simd_deperiodizef md_simd_deperiodize_f256

    #define md_simd_horizontal_minf md_simd_horizontal_min_f256
    #define md_simd_horizontal_maxf md_simd_horizontal_max_f256
    #define md_simd_horizontal_addf md_simd_horizontal_add_f256

    #define md_simd_convert_i_to_f md_simd_convert_to_f256_i256
    #define md_simd_convert_f_to_i md_simd_convert_to_i256_f256

    #define md_simd_cast_i_to_f md_simd_cast_f256_i256
    #define md_simd_cast_f_to_i md_simd_cast_i256_f256

// Int
    #define md_simd_widthi 8
    #define md_simd_typei   md_simd_i256_t
    #define md_simd_set1i   md_simd_set1_i256
    #define md_simd_zeroi   md_simd_zero_i256
    #define md_simd_ori     md_simd_or_i256
    #define md_simd_andi    md_simd_and_i256
    #define md_simd_andnoti md_simd_andnot_i256
    #define md_simd_xori    md_simd_xor_i256
    #define md_simd_noti    md_simd_not_i256
    #define md_simd_addi    md_simd_add_i256
    #define md_simd_subi    md_simd_sub_i256
    #define md_simd_cmp_eqi md_simd_cmp_eq_i256
    #define md_simd_shift_lefti md_simd_shift_left_i256

#else
    #define md_simd_widthf 4
    #define md_simd_typef  md_simd_f128_t
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
    #define md_simd_signcopyf   md_simd_signcopy_f128
    #define md_simd_fractf  md_simd_fract_f128
    #define md_simd_roundf  md_simd_round_f128
    #define md_simd_blendf  md_simd_blend_f128

    #define md_simd_horizontal_minf md_simd_horizontal_min_f128
    #define md_simd_horizontal_maxf md_simd_horizontal_max_f128
    #define md_simd_horizontal_addf md_simd_horizontal_add_f128

    #define md_simd_convert_i_to_f md_simd_convert_to_f128_i128
    #define md_simd_convert_f_to_i md_simd_convert_to_i128_f128

    #define md_simd_cast_i_to_f md_simd_cast_f128_i128
    #define md_simd_cast_f_to_i md_simd_cast_i128_f128

// Int
    #define md_simd_widthi 4
    #define md_simd_typei   md_simd_i128_t
    #define md_simd_set1i   md_simd_set1_i128
    #define md_simd_zeroi   md_simd_zero_i128
    #define md_simd_ori     md_simd_or_i128
    #define md_simd_andi    md_simd_and_i128
    #define md_simd_andnoti md_simd_andnot_i128
    #define md_simd_xori    md_simd_xor_i128
    #define md_simd_noti    md_simd_not_i128
    #define md_simd_addi    md_simd_add_i128
    #define md_simd_subi    md_simd_sub_i128
    #define md_simd_cmp_eqi md_simd_cmp_eq_i128
    #define md_simd_shift_lefti md_simd_shift_left_i128

#endif

// 128-bit wide
// float operations

MD_SIMD_INLINE __m128 md_simd_set1_f128(float x) { return _mm_set1_ps(x); }
MD_SIMD_INLINE __m128 md_simd_set_f128(float x, float y, float z, float w) { return _mm_set_ps(w, z, y, x); }

MD_SIMD_INLINE __m128 md_simd_zero_f128() { return _mm_setzero_ps(); }

MD_SIMD_INLINE __m128 md_simd_load_f128(const float* addr) { return _mm_loadu_ps(addr); }
MD_SIMD_INLINE __m128 md_simd_load_aligned_f128(const float* addr) {
    ASSERT(IS_ALIGNED(addr, 16));
    return _mm_load_ps(addr);
}

MD_SIMD_INLINE void md_simd_store_f128(float* addr, __m128 v) { _mm_storeu_ps(addr, v); }
MD_SIMD_INLINE void md_simd_store_aligned_f128(float* addr, __m128 v) {
    ASSERT(IS_ALIGNED(addr, 16));
    _mm_store_ps(addr, v);
}

MD_SIMD_INLINE __m128 md_simd_add_f128(__m128 a, __m128 b) { return _mm_add_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_sub_f128(__m128 a, __m128 b) { return _mm_sub_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_mul_f128(__m128 a, __m128 b) { return _mm_mul_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_div_f128(__m128 a, __m128 b) { return _mm_div_ps(a, b); }

MD_SIMD_INLINE __m128 md_simd_and_f128    (__m128 a, __m128 b) { return _mm_and_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_and_not_f128(__m128 a, __m128 b) { return _mm_andnot_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_or_f128     (__m128 a, __m128 b) { return _mm_or_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_xor_f128    (__m128 a, __m128 b) { return _mm_xor_ps(a, b); }

MD_SIMD_INLINE __m128 md_simd_cmp_gt_f128 (__m128 a, __m128 b) { return _mm_cmpgt_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_cmp_ge_f128 (__m128 a, __m128 b) { return _mm_cmpge_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_cmp_lt_f128 (__m128 a, __m128 b) { return _mm_cmplt_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_cmp_le_f128 (__m128 a, __m128 b) { return _mm_cmple_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_cmp_eq_f128 (__m128 a, __m128 b) { return _mm_cmpeq_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_cmp_neq_f128(__m128 a, __m128 b) { return _mm_cmpneq_ps(a, b); }

MD_SIMD_INLINE __m128 md_simd_abs_f128(__m128 a) { return _mm_and_ps(a, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
MD_SIMD_INLINE __m128 md_simd_sign_f128(__m128 x) {
    __m128 sgn = _mm_and_ps(x, _mm_set1_ps(-0.0f));
    __m128 res = _mm_xor_ps(sgn, _mm_set1_ps(1.0f));
    return res;
}

MD_SIMD_INLINE __m128 md_simd_copysign_f128(__m128 mag, __m128 sign) {
    __m128 const mask = _mm_set1_ps(-0.0f);
    return _mm_or_ps(_mm_and_ps(mask, mag), _mm_andnot_ps(mask, sign));
}

MD_SIMD_INLINE __m128 md_simd_fract_f128(__m128 x) {
    return _mm_sub_ps(x, _mm_round_ps(x, _MM_FROUND_TRUNC | _MM_FROUND_NO_EXC));
}

MD_SIMD_INLINE __m128 md_simd_round_f128(__m128 x) {
    return _mm_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}

MD_SIMD_INLINE __m128 md_simd_blend_f128(__m128 a, __m128 b, __m128 mask) {
    return _mm_blendv_ps(a, b, mask);
}

MD_SIMD_INLINE __m128 md_simd_sqrt_f128(__m128 x) {
    return _mm_sqrt_ps(x);
}
MD_SIMD_INLINE __m128 md_simd_min_f128(__m128 a, __m128 b) { return _mm_min_ps(a, b); }
MD_SIMD_INLINE __m128 md_simd_max_f128(__m128 a, __m128 b) { return _mm_max_ps(a, b); }

MD_SIMD_INLINE float md_simd_horizontal_min_f128(__m128 x) {
    __m128 min1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
    __m128 min2 = _mm_min_ps(x, min1);
    __m128 min3 = _mm_shuffle_ps(min2, min2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128 min4 = _mm_min_ps(min2, min3);
    return _mm_cvtss_f32(min4);
}

MD_SIMD_INLINE float md_simd_horizontal_max_f128(__m128 x) {
    __m128 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
    __m128 max2 = _mm_max_ps(x, max1);
    __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
    __m128 max4 = _mm_max_ps(max2, max3);
    return _mm_cvtss_f32(max4);
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
MD_SIMD_INLINE float md_simd_horizontal_add_f128(__m128 x) {
    __m128 shuf = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));
    __m128 sums = _mm_add_ps(x, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

MD_SIMD_INLINE __m128 md_simd_cubic_spline_f128(__m128 p0, __m128 p1, __m128 p2, __m128 p3, __m128 t, __m128 tension) {
    const __m128 vt = tension;
    const __m128 t1 = t;
    const __m128 t2 = _mm_mul_ps(t, t);
    const __m128 t3 = _mm_mul_ps(t2, t);
    const __m128 v0 = _mm_mul_ps(_mm_sub_ps(p2, p0), vt);
    const __m128 v1 = _mm_mul_ps(_mm_sub_ps(p3, p1), vt);
    const __m128 x0 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(2), _mm_sub_ps(p1, p2)), _mm_add_ps(v0, v1));
    const __m128 x1 = _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(3), _mm_sub_ps(p2, p1)), _mm_add_ps(_mm_mul_ps(_mm_set1_ps(2), v0), v1));
    const __m128 r0 = _mm_add_ps(_mm_mul_ps(x0, t3), _mm_mul_ps(x1, t2));
    const __m128 r1 = _mm_add_ps(_mm_mul_ps(v0, t1), p1);
    return _mm_add_ps(r0, r1);
}

// Cast
MD_SIMD_INLINE __m128  md_simd_cast_f128_i128(__m128i x) { return _mm_castsi128_ps(x); }
MD_SIMD_INLINE __m128i md_simd_cast_i128_f128(__m128  x) { return _mm_castps_si128(x); }

// Conversion
MD_SIMD_INLINE __m128  md_simd_convert_to_f128_i128(__m128i x)    { return _mm_cvtepi32_ps(x); }
MD_SIMD_INLINE __m128i md_simd_convert_to_i128_f128(__m128  x)    { return _mm_cvtps_epi32(x); }

// Int operations
MD_SIMD_INLINE __m128i md_simd_set1_i128(int x) { return _mm_set1_epi32(x); }
MD_SIMD_INLINE __m128i md_simd_set_i128(int x, int y, int z, int w) { return _mm_set_epi32(w, z, y, x); }
MD_SIMD_INLINE __m128i md_simd_zero_i128() { return _mm_setzero_si128(); }

MD_SIMD_INLINE __m128i md_simd_or_i128(__m128i a, __m128i b)     { return _mm_or_si128(a, b); }
MD_SIMD_INLINE __m128i md_simd_and_i128(__m128i a, __m128i b)    { return _mm_and_si128(a, b); }
MD_SIMD_INLINE __m128i md_simd_andnot_i128(__m128i a, __m128i b) { return _mm_andnot_si128(a, b); }
MD_SIMD_INLINE __m128i md_simd_xor_i128(__m128i a, __m128i b)    { return _mm_xor_si128(a, b); }
MD_SIMD_INLINE __m128i md_simd_not_i128(__m128i x)               { return _mm_andnot_si128(x, _mm_set1_epi64x(-1)); }

MD_SIMD_INLINE __m128i md_simd_add_i128(__m128i a, __m128i b)    { return _mm_add_epi32(a, b); }
MD_SIMD_INLINE __m128i md_simd_sub_i128(__m128i a, __m128i b)    { return _mm_sub_epi32(a, b); }

MD_SIMD_INLINE __m128i md_simd_cmp_eq_i128(__m128i a, __m128i b) { return _mm_cmpeq_epi32(a, b); }

MD_SIMD_INLINE __m128i md_simd_shift_left_i128(__m128i a, int count) { return _mm_slli_epi32(a, count); }

#ifdef __AVX__
// 256-bit wide
// float operations

MD_SIMD_INLINE __m256 md_simd_set1_f256(float x) { return _mm256_set1_ps(x); }
MD_SIMD_INLINE __m256 md_simd_set_f256(float x0, float y0, float z0, float w0, float x1, float y1, float z1, float w1) { return _mm256_set_ps(w1, z1, y1, x1, w0, z0, y0, x0); }

MD_SIMD_INLINE __m256 md_simd_zero_f256() { return _mm256_setzero_ps(); }

MD_SIMD_INLINE __m256 md_simd_load_f256(const float* addr) { return _mm256_loadu_ps(addr); }
MD_SIMD_INLINE __m256 md_simd_load_aligned_f256(const float* addr) {
    ASSERT(IS_ALIGNED(addr, 16));
    return _mm256_load_ps(addr);
}

MD_SIMD_INLINE void md_simd_store_f256(float* addr, __m256 v) { _mm256_storeu_ps(addr, v); }
MD_SIMD_INLINE void md_simd_store_aligned_f256(float* addr, __m256 v) {
    ASSERT(IS_ALIGNED(addr, 16));
    _mm256_store_ps(addr, v);
}

MD_SIMD_INLINE __m256 md_simd_add_f256(__m256 a, __m256 b) { return _mm256_add_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_sub_f256(__m256 a, __m256 b) { return _mm256_sub_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_mul_f256(__m256 a, __m256 b) { return _mm256_mul_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_div_f256(__m256 a, __m256 b) { return _mm256_div_ps(a, b); }

MD_SIMD_INLINE __m256 md_simd_and_f256    (__m256 a, __m256 b) { return _mm256_and_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_and_not_f256(__m256 a, __m256 b) { return _mm256_andnot_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_or_f256     (__m256 a, __m256 b) { return _mm256_or_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_xor_f256    (__m256 a, __m256 b) { return _mm256_xor_ps(a, b); }

MD_SIMD_INLINE __m256 md_simd_cmp_gt_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_GT_OQ); }
MD_SIMD_INLINE __m256 md_simd_cmp_ge_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_GE_OQ); }
MD_SIMD_INLINE __m256 md_simd_cmp_lt_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }
MD_SIMD_INLINE __m256 md_simd_cmp_le_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_LE_OQ); }
MD_SIMD_INLINE __m256 md_simd_cmp_eq_f256 (__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }
MD_SIMD_INLINE __m256 md_simd_cmp_neq_f256(__m256 a, __m256 b) { return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ); }

MD_SIMD_INLINE __m256 md_simd_abs_f256(__m256 a) { return _mm256_and_ps(a, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF))); }

// @NOTE: If 0.0f is given as input, it will be mapped to 1.0f
MD_SIMD_INLINE __m256 md_simd_sign_f256(__m256 x) {
    __m256 sgn = _mm256_and_ps(x, _mm256_set1_ps(-0.0f));
    __m256 res = _mm256_xor_ps(sgn, _mm256_set1_ps(1.0f));
    return res;
}

MD_SIMD_INLINE __m256 md_simd_copysign_f256(__m256 mag, __m256 sign) {
    __m256 const mask = _mm256_set1_ps(-0.0f);
    return _mm256_or_ps(_mm256_and_ps(mask, mag), _mm256_andnot_ps(mask, sign));
}

MD_SIMD_INLINE __m256 md_simd_fract_f256(__m256 x) {
    return _mm256_sub_ps(x, _mm256_round_ps(x, _MM_FROUND_TRUNC | _MM_FROUND_NO_EXC));
}

MD_SIMD_INLINE __m256 md_simd_round_f256(__m256 x) {
    return _mm256_round_ps(x, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
}

MD_SIMD_INLINE __m256 md_simd_blend_f256(__m256 a, __m256 b, __m256 mask) {
    return _mm256_blendv_ps(a, b, mask);
}

MD_SIMD_INLINE __m256 md_simd_min_f256(__m256 a, __m256 b) { return _mm256_min_ps(a, b); }
MD_SIMD_INLINE __m256 md_simd_max_f256(__m256 a, __m256 b) { return _mm256_max_ps(a, b); }

MD_SIMD_INLINE float md_simd_horizontal_min_f256(__m256 x) {
    __m128 lo = _mm256_castps256_ps128(x);
    __m128 hi = _mm256_extractf128_ps(x, 0x1);
    __m128 min_val = _mm_min_ps(lo, hi);
    return md_simd_horizontal_min_f128(min_val);
}

MD_SIMD_INLINE float md_simd_horizontal_max_f256(__m256 x) {
    __m128 lo = _mm256_castps256_ps128(x);
    __m128 hi = _mm256_extractf128_ps(x, 0x1);
    __m128 max_val = _mm_max_ps(lo, hi);
    return md_simd_horizontal_max_f128(max_val);
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
MD_SIMD_INLINE float md_simd_horizontal_add_f256(__m256 x) {
    __m128 vlow = _mm256_castps256_ps128(x);
    __m128 vhigh = _mm256_extractf128_ps(x, 1); // high 128
    vlow = _mm_add_ps(vlow, vhigh);             // add the low 128
    __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(vlow, shuf);
    shuf = _mm_movehl_ps(shuf, sums);           // high half -> low half
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

// Cast
MD_SIMD_INLINE __m256  md_simd_cast_f256_i256(__m256i x) { return _mm256_castsi256_ps(x); }
MD_SIMD_INLINE __m256i md_simd_cast_i256_f256(__m256  x) { return _mm256_castps_si256(x); }

// Conversions
MD_SIMD_INLINE __m256  md_simd_convert_to_f256_i256(__m256i x)    { return _mm256_cvtepi32_ps(x); }
MD_SIMD_INLINE __m256i md_simd_convert_to_i256_f256(__m256  x)    { return _mm256_cvtps_epi32(x); }

// Int operations

MD_SIMD_INLINE __m256i md_simd_set1_i256(int x) { return _mm256_set1_epi32(x); }
MD_SIMD_INLINE __m256i md_simd_set_i256(int x0, int y0, int z0, int w0, int x1, int y1, int z1, int w1) { return _mm256_set_epi32(w1, z1, y1, x1, w0, z0, y0, x0); }
MD_SIMD_INLINE __m256i md_simd_zero_i256() { return _mm256_setzero_si256(); }


MD_SIMD_INLINE __m128i md_simd_extract_lo_i256(__m256i x) { return _mm256_castsi256_si128(x); }
MD_SIMD_INLINE __m128i md_simd_extract_hi_i256(__m256i x) { return _mm256_extractf128_si256(x, 1); }

#ifdef __AVX2__
MD_SIMD_INLINE __m256i md_simd_or_i256(__m256i a, __m256i b)     { return _mm256_or_si256(a, b); }
MD_SIMD_INLINE __m256i md_simd_and_i256(__m256i a, __m256i b)    { return _mm256_and_si256(a, b); }
MD_SIMD_INLINE __m256i md_simd_andnot_i256(__m256i a, __m256i b) { return _mm256_andnot_si256(a, b); }
MD_SIMD_INLINE __m256i md_simd_xor_i256(__m256i a, __m256i b)    { return _mm256_xor_si256(a, b); }
MD_SIMD_INLINE __m256i md_simd_not_i256(__m256i x)               { return _mm256_andnot_si256(x, _mm256_set1_epi64x(-1)); }

MD_SIMD_INLINE __m256i md_simd_add_i256(__m256i a, __m256i b)    { return _mm256_add_epi32(a, b); }
MD_SIMD_INLINE __m256i md_simd_sub_i256(__m256i a, __m256i b)    { return _mm256_sub_epi32(a, b); }

MD_SIMD_INLINE __m256i md_simd_cmp_eq_i256(__m256i a, __m256i b) { return _mm256_cmpeq_epi32(a, b); }

MD_SIMD_INLINE __m256i md_simd_shift_left_i256(__m256i a, int count) { return _mm256_slli_epi32(a, count); }


#else
// Fallback in case AVX2 is not supported
MD_SIMD_INLINE __m256i md_simd_or_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_or_i128(a_hi, b_hi), md_simd_or_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_and_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_and_i128(a_hi, b_hi), md_simd_and_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_andnot_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_andnot_i128(a_hi, b_hi), md_simd_andnot_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_xor_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_xor_i128(a_hi, b_hi), md_simd_xor_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_not_i256(__m256i x) {
    __m128i x_lo = md_simd_extract_lo_i256(x);
    __m128i x_hi = md_simd_extract_hi_i256(x);
    return _mm256_set_m128i(md_simd_not_i128(x_hi), md_simd_not_i128(x_lo));
}

MD_SIMD_INLINE __m256i md_simd_add_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_add_i128(a_hi, b_hi), md_simd_add_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_sub_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_sub_i128(a_hi, b_hi), md_simd_sub_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_cmp_eq_i256(__m256i a, __m256i b) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    __m128i b_lo = md_simd_extract_lo_i256(b);
    __m128i b_hi = md_simd_extract_hi_i256(b);
    return _mm256_set_m128i(md_simd_cmp_eq_i128(a_hi, b_hi), md_simd_cmp_eq_i128(a_lo, b_lo));
}

MD_SIMD_INLINE __m256i md_simd_shift_left_i256(__m256i a, int count) {
    __m128i a_lo = md_simd_extract_lo_i256(a);
    __m128i a_hi = md_simd_extract_hi_i256(a);
    return _mm256_set_m128i(md_simd_shift_left_i128(a_hi, count), md_simd_shift_left_i128(a_lo, count));
}

#endif



#endif

// Common higher level operation which are not specific for any instruction set

MD_SIMD_INLINE md_simd_typef md_simd_stepf(md_simd_typef edge, md_simd_typef x) {
    const md_simd_typef cmp = md_simd_cmp_gef(x, edge);
    return md_simd_andf(cmp, md_simd_set1f(1.f));
}

MD_SIMD_INLINE md_simd_typef md_simd_lerpf(md_simd_typef a, md_simd_typef b, float t) {
    const md_simd_typef one_minus_alpha = md_simd_set1f(1.0f - t);
    const md_simd_typef alpha = md_simd_set1f(t);
    return md_simd_addf(md_simd_mulf(a, one_minus_alpha), md_simd_mulf(b, alpha));
}

MD_SIMD_INLINE md_simd_typef md_simd_cubic_splinef(md_simd_typef p0, md_simd_typef p1, md_simd_typef p2, md_simd_typef p3, float s, float tension) {
    const md_simd_typef vt = md_simd_set1f(tension);
    const md_simd_typef s1 = md_simd_set1f(s);
    const md_simd_typef s2 = md_simd_mulf(s1, s1);
    const md_simd_typef s3 = md_simd_mulf(s2, s1);
    const md_simd_typef v0 = md_simd_mulf(md_simd_subf(p2, p0), vt);
    const md_simd_typef v1 = md_simd_mulf(md_simd_subf(p3, p1), vt);
    const md_simd_typef x0 = md_simd_addf(md_simd_mulf(md_simd_set1f(2), md_simd_subf(p1, p2)), md_simd_addf(v0, v1));
    const md_simd_typef x1 = md_simd_subf(md_simd_mulf(md_simd_set1f(3), md_simd_subf(p2, p1)), md_simd_addf(md_simd_mulf(md_simd_set1f(2), v0), v1));
    const md_simd_typef r0 = md_simd_addf(md_simd_mulf(x0, s3), md_simd_mulf(x1, s2));
    const md_simd_typef r1 = md_simd_addf(md_simd_mulf(v0, s1), p1);
    return md_simd_addf(r0, r1);
}

MD_SIMD_INLINE void md_simd_sincosf(md_simd_typef x, md_simd_typef *s, md_simd_typef *c) {
    md_simd_typef xmm1, xmm2, xmm3, sign_bit_sin, y;
    md_simd_typei imm0, imm2, imm4;

    sign_bit_sin = x;
    /* take the absolute value */
    x = md_simd_absf(x);
    /* extract the sign bit (upper one) */
    sign_bit_sin = md_simd_andf(x, md_simd_set1f(-0.0f));

    /* scale by 4/Pi */
    y = md_simd_mulf(x, md_simd_set1f(1.27323954473516)); // 4 / M_PI

    /* store the integer part of y in imm2 */
    imm2 = md_simd_convert_f_to_i(y);

    /* j=(j+1) & (~1) (see the cephes sources) */
    
    imm2 = md_simd_addi(imm2, md_simd_set1i(1));
    imm2 = md_simd_andi(imm2, md_simd_set1i(~1));

    y = md_simd_convert_i_to_f(imm2);
    imm4 = imm2;

    /* get the swap sign flag for the sine */
    imm0 = md_simd_andi(imm2, md_simd_set1i(4));
    imm0 = md_simd_shift_lefti(imm0, 29);
    //md_simd_typef swap_sign_bit_sin = _mm256_castsi256_ps(imm0);

    /* get the polynom selection mask for the sine*/
    imm2 = md_simd_andi(imm2, md_simd_set1i(2));
    imm2 = md_simd_cmp_eqi(imm2, md_simd_zeroi());
    //md_simd_typef poly_mask = _mm256_castsi256_ps(imm2);

    md_simd_typef swap_sign_bit_sin = md_simd_cast_i_to_f(imm0);
    md_simd_typef poly_mask = md_simd_cast_i_to_f(imm2);

    /* The magic pass: "Extended precision modular arithmetic" 
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = md_simd_set1f(-0.78515625);
    xmm2 = md_simd_set1f(-2.4187564849853515625e-4);
    xmm3 = md_simd_set1f(-3.77489497744594108e-8);
    xmm1 = md_simd_mulf(y, xmm1);
    xmm2 = md_simd_mulf(y, xmm2);
    xmm3 = md_simd_mulf(y, xmm3);
    x = md_simd_addf(x, xmm1);
    x = md_simd_addf(x, xmm2);
    x = md_simd_addf(x, xmm3);

    imm4 = md_simd_subi(imm4, md_simd_set1i(2));
    imm4 = md_simd_andnoti(imm4, md_simd_set1i(4));
    imm4 = md_simd_shift_lefti(imm4, 29);

    md_simd_typef sign_bit_cos = md_simd_convert_i_to_f(imm4);

    sign_bit_sin = md_simd_xorf(sign_bit_sin, swap_sign_bit_sin);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    md_simd_typef z = md_simd_mulf(x,x);
    y = md_simd_set1f(2.443315711809948E-005);

    y = md_simd_mulf(y, z);
    y = md_simd_addf(y, md_simd_set1f(-1.388731625493765E-003));
    y = md_simd_mulf(y, z);
    y = md_simd_addf(y, md_simd_set1f(4.166664568298827E-002));
    y = md_simd_mulf(y, z);
    y = md_simd_mulf(y, z);
    md_simd_typef tmp = md_simd_mulf(z, md_simd_set1f(0.5));
    y = md_simd_subf(y, tmp);
    y = md_simd_addf(y, md_simd_set1f(1));

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

    md_simd_typef y2 = md_simd_set1f(-1.9515295891E-4);
    y2 = md_simd_mulf(y2, z);
    y2 = md_simd_addf(y2, md_simd_set1f(8.3321608736E-3));
    y2 = md_simd_mulf(y2, z);
    y2 = md_simd_addf(y2, md_simd_set1f(-1.6666654611E-1));
    y2 = md_simd_mulf(y2, z);
    y2 = md_simd_mulf(y2, x);
    y2 = md_simd_addf(y2, x);

    /* select the correct result from the two polynoms */  
    xmm3 = poly_mask;
    md_simd_typef ysin2 = md_simd_andf(xmm3, y2);
    md_simd_typef ysin1 = md_simd_and_notf(xmm3, y);
    y2 = md_simd_subf(y2,ysin2);
    y = md_simd_subf(y, ysin1);

    xmm1 = md_simd_addf(ysin1,ysin2);
    xmm2 = md_simd_addf(y,y2);

    /* update the sign */
    *s = md_simd_xorf(xmm1, sign_bit_sin);
    *c = md_simd_xorf(xmm2, sign_bit_cos);
}
