/**
 * This is an attempt of a generic SIMD interface.
 * It is not complete, but it is a start.
 * The idea is to support SSE2, AVX (AVX2) and NEON.
 * Using defined types (based on hardware) and generics on top of operations
 * 
 * TODO: Cleanup
**/

// If the instruction set supports AVX2 it should support FMA as well.
// However, we cannot rely on that (in case of emulation or whatever may gimp your instruction set)
// so we check for it in cmake and define if it is available

#pragma once

#include <core/md_common.h>
#include <core/md_intrinsics.h>

#if defined(__x86_64__)
// x86 integer intrinsics does not have distinct types
// so we define distinct types for each in order to enable C11 Generics
typedef struct md_i32x4_t {
    __m128i m128i;
} md_i32x4_t;

typedef struct md_i64x2_t {
    __m128i m128i;
} md_i64x2_t;

typedef struct md_i32x8_t {
    __m256i m256i;
} md_i32x8_t;

typedef struct md_i64x4_t {
    __m256i m256i;
} md_i64x4_t;

#   define md_f32x4_t __m128
#   define md_f64x2_t __m128d
#   define md_f32x8_t __m256
#   define md_f64x4_t __m256d
#   if defined(__AVX__)
#       define md_simd_f32_t    md_f32x8_t
#       define md_simd_f64_t    md_f64x4_t
#       define md_simd_i32_t    md_i32x8_t
#       define md_simd_i64_t    md_i64x4_t
#       define md_simd_width_f32 8
#       define md_simd_width_i32 8
#       define md_simd_width_f64 4
#       define md_simd_width_i64 4
#   else
#       define md_simd_f32_t    md_f32x4_t
#       define md_simd_f64_t    md_f64x2_t
#       define md_simd_i32_t    md_i32x4_t
#       define md_simd_i64_t    md_i64x2_t
#       define md_simd_width_f32 4
#       define md_simd_width_i32 4
#       define md_simd_width_f64 2
#       define md_simd_width_i64 2
#   endif
#elif defined(__ARM_NEON__)
#   error "NEON not supported yet"
#else
#   error "Unsupported platform"
#endif

#define MD_SIMD_INLINE static FORCE_INLINE

#ifdef __x86_64__

// Integers have to be explicitly handled since we use distinct types
MD_SIMD_INLINE void md_simd_store_i32x4(int* x, md_i32x4_t val) { _mm_storeu_si128((__m128i*)x, val.m128i); }
MD_SIMD_INLINE void md_simd_store_i64x2(int64_t* x, md_i64x2_t val) { _mm_storeu_si128((__m128i*)x, val.m128i); }
MD_SIMD_INLINE void md_simd_store_i32x8(int* x, md_i32x8_t val) { _mm256_storeu_si256((__m256i*)x, val.m256i); }
MD_SIMD_INLINE void md_simd_store_i64x4(int64_t* x, md_i64x4_t val) { _mm256_storeu_si256((__m256i*)x, val.m256i); }

MD_SIMD_INLINE md_i32x4_t md_simd_load_i32x4(const int* x) {
    md_i32x4_t val = { _mm_loadu_si128((const __m128i*)x) };
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_load_i64x2(const int64_t* x) {
    md_i64x2_t val = {_mm_loadu_si128((const __m128i*)x)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_load_i32x8(const int* x) {
    md_i32x8_t val = {_mm256_loadu_si256((const __m256i*)x)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_load_i64x4(const int64_t* x) {
    md_i64x4_t val = {_mm256_loadu_si256((const __m256i*)x)};
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_set1_i32x4(int x) {
    md_i32x4_t val = {_mm_set1_epi32(x)};
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_set1_i64x2(int64_t x) {
    md_i64x2_t val = { _mm_set1_epi64x(x) };
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_set1_i32x8(int x) {
    md_i32x8_t val = {_mm256_set1_epi32(x)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_set1_i64x4(int64_t x) {
    md_i64x4_t val = {_mm256_set1_epi64x(x)};
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_set_i32x4(int x, int y, int z, int w) {
    md_i32x4_t val = {_mm_set_epi32(w, z, y, x)};
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_set_i64x2(int64_t x, int64_t y) {
    md_i64x2_t val = {_mm_set_epi64x(y, x)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_set_i32x8(int x, int y, int z, int w, int q, int r, int s, int t) {
    md_i32x8_t val = {_mm256_set_epi32(t, s, r, q, w, z, y, x)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_set_i64x4(int64_t x, int64_t y, int64_t z, int64_t w) {
    md_i64x4_t val = {_mm256_set_epi64x(w, z, y, x)};
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_zero_i32x4() {
    md_i32x4_t val = {_mm_setzero_si128()};
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_zero_i64x2() {
    md_i64x2_t val = {_mm_setzero_si128()};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_zero_i32x8() {
    md_i32x8_t val = {_mm256_setzero_si256()};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_zero_i64x4() {
    md_i64x4_t val = {_mm256_setzero_si256()};
    return val;
}

#define MD_SIMD_EXTRACT_HI_SI128(X) _mm256_extractf128_si256(X, 1)
#define MD_SIMD_EXTRACT_LO_SI128(X) _mm256_castsi256_si128(X)
#define MD_SIMD_DOUBLE_PUMP1_SI128(X, OP) _mm256_set_m128i(OP(MD_SIMD_EXTRACT_HI_SI128(X)), OP(MD_SIMD_EXTRACT_LO_SI128(X)))
#define MD_SIMD_DOUBLE_PUMP1_ARGS_SI128(X, OP, ...) _mm256_set_m128i(OP(MD_SIMD_EXTRACT_HI_SI128(X),##__VA_ARGS__), OP(MD_SIMD_EXTRACT_LO_SI128(X),##__VA_ARGS__))
#define MD_SIMD_DOUBLE_PUMP2_SI128(A, B, OP) _mm256_set_m128i(OP(MD_SIMD_EXTRACT_HI_SI128(A), MD_SIMD_EXTRACT_HI_SI128(B)), OP(MD_SIMD_EXTRACT_LO_SI128(A), MD_SIMD_EXTRACT_LO_SI128(B)))

MD_SIMD_INLINE md_i32x4_t md_simd_and_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_and_si128(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_and_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_and_si256(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_and_si128)};
#endif
    return val;
}
    
MD_SIMD_INLINE md_i64x2_t md_simd_and_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_and_si128(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_and_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_and_si256(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_and_si128)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_or_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_or_si128(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_or_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_or_si256(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_or_si128)};
#endif
	return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_or_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_or_si128(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_or_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_or_si256(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_or_si128)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_xor_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_xor_si128(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_xor_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_xor_si256(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_xor_si128)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_xor_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_xor_si128(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_xor_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_xor_si256(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_xor_si128)};
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_not_i32x4(md_i32x4_t a) {
    return md_simd_xor_i32x4(a, md_simd_set1_i32x4(-1));
}

MD_SIMD_INLINE md_i32x8_t md_simd_not_i32x8(md_i32x8_t a) {
	return md_simd_xor_i32x8(a, md_simd_set1_i32x8(-1));
}

MD_SIMD_INLINE md_i64x2_t md_simd_not_i64x2(md_i64x2_t a) {
    return md_simd_xor_i64x2(a, md_simd_set1_i64x2(-1));
}

MD_SIMD_INLINE md_i64x4_t md_simd_not_i64x4(md_i64x4_t a) {
    return md_simd_xor_i64x4(a, md_simd_set1_i64x4(-1));
}

MD_SIMD_INLINE md_i32x4_t md_simd_and_not_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_andnot_si128(b.m128i, a.m128i)};
	return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_and_not_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_andnot_si256(b.m256i, a.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(b.m256i, a.m256i, _mm_andnot_si128)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_and_not_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_andnot_si128(b.m128i, a.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_and_not_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_andnot_si256(b.m256i, a.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(b.m256i, a.m256i, _mm_andnot_si128)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_add_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_add_epi32(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_add_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_add_epi32(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_add_epi32)};
#endif
	return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_add_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_add_epi64(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_add_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_add_epi64(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_add_epi64)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_sub_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_sub_epi32(a.m128i, b.m128i)};
	return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_sub_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_sub_epi32(a.m256i, b.m256i)};
#else
        md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_sub_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_sub_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_sub_epi64(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_sub_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_sub_epi64(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_sub_epi64)};
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_abs_i32x4(md_i32x4_t a) {
    md_i32x4_t val = {_mm_abs_epi32(a.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_abs_i32x8(md_i32x8_t a) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_abs_epi32(a.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP1_SI128(a.m256i, _mm_abs_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_abs_i64x2(md_i64x2_t a) {
    md_i64x2_t val = {_mm_abs_epi64(a.m128i)};
	return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_abs_i64x4(md_i64x4_t a) {
    md_i64x4_t val = {_mm256_abs_epi64(a.m256i)};
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_min_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_min_epi32(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_min_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_min_epi32(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_min_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_min_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_min_epi64(a.m128i, b.m128i)};
	return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_min_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_min_epi64(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_min_epi64)};
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_max_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_max_epi32(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_max_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_max_epi32(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_max_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_max_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_max_epi64(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_max_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_max_epi64(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_max_epi64)};
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_cmp_gt_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_cmpgt_epi32(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_cmp_gt_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_cmpgt_epi32(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmpgt_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_cmp_gt_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_cmpgt_epi64(a.m128i, b.m128i)};
	return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_cmp_gt_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_cmpgt_epi64(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmpgt_epi64)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_cmp_eq_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_cmpeq_epi32(a.m128i, b.m128i)};
	return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_cmp_eq_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_cmpeq_epi32(a.m256i, b.m256i)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmpeq_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_cmp_eq_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = {_mm_cmpeq_epi64(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_cmp_eq_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_cmpeq_epi64(a.m256i, b.m256i)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmpeq_epi64)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_cmp_neq_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = md_simd_not_i32x4(md_simd_cmp_eq_i32x4(a, b));
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_cmp_neq_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = md_simd_not_i32x8(md_simd_cmp_eq_i32x8(a, b));
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmpeq_epi32)};
    val = md_simd_not_i32x8(val);
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_cmp_neq_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = md_simd_not_i64x2(md_simd_cmp_eq_i64x2(a, b));
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_cmp_neq_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = md_simd_not_i64x4(md_simd_cmp_eq_i64x4(a, b));
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmpeq_epi64)};
    val = md_simd_not_i64x4(val);
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_cmp_lt_i32x4(md_i32x4_t a, md_i32x4_t b) {
    md_i32x4_t val = {_mm_cmplt_epi32(a.m128i, b.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_cmp_lt_i32x8(md_i32x8_t a, md_i32x8_t b) {
#ifdef __AVX2__
    md_i32x8_t val = md_simd_cmp_gt_i32x8(b, a);
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(a.m256i, b.m256i, _mm_cmplt_epi32)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_cmp_lt_i64x2(md_i64x2_t a, md_i64x2_t b) {
    md_i64x2_t val = md_simd_not_i64x2(md_simd_cmp_gt_i64x2(a, b));
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_cmp_lt_i64x4(md_i64x4_t a, md_i64x4_t b) {
#ifdef __AVX2__
    md_i64x4_t val = md_simd_cmp_gt_i64x4(b, a);
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP2_SI128(b.m256i, a.m256i, _mm_cmpgt_epi64)};
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_shift_left_i32x4(md_i32x4_t a, int32_t b) {
    md_i32x4_t val = {_mm_slli_epi32(a.m128i, b)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_shift_left_i32x8(md_i32x8_t a, int32_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_slli_epi32(a.m256i, b)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP1_ARGS_SI128(a.m256i, _mm_slli_epi32, b)};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_shift_left_i64x2(md_i64x2_t a, int32_t b) {
    md_i64x2_t val = {_mm_slli_epi64(a.m128i, b)};
	return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_shift_left_i64x4(md_i64x4_t a, int32_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_slli_epi64(a.m256i, b)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP1_ARGS_SI128(a.m256i, _mm_slli_epi64, b)};
#endif
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_shift_right_i32x4(md_i32x4_t a, int32_t b) {
    md_i32x4_t val = {_mm_srli_epi32(a.m128i, b)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_shift_right_i32x8(md_i32x8_t a, int32_t b) {
#ifdef __AVX2__
    md_i32x8_t val = {_mm256_srli_epi32(a.m256i, b)};
#else
    md_i32x8_t val = {MD_SIMD_DOUBLE_PUMP1_ARGS_SI128(a.m256i, _mm_srli_epi32, b)};
#endif
	return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_shift_right_i64x2(md_i64x2_t a, int32_t b) {
    md_i64x2_t val = {_mm_srli_epi64(a.m128i, b)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_shift_right_i64x4(md_i64x4_t a, int32_t b) {
#ifdef __AVX2__
    md_i64x4_t val = {_mm256_srli_epi64(a.m256i, b)};
#else
    md_i64x4_t val = {MD_SIMD_DOUBLE_PUMP1_ARGS_SI128(a.m256i, _mm_srli_epi64, b)};
#endif
	return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_blend_i32x4(md_i32x4_t a, md_i32x4_t b, md_i32x4_t mask) {
    md_i32x4_t val = {_mm_blendv_epi8(a.m128i, b.m128i, mask.m128i)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_blend_i32x8(md_i32x8_t a, md_i32x8_t b, md_i32x8_t mask) {
#if __AVX2__
    md_i32x8_t val = {_mm256_blendv_epi8(a.m256i, b.m256i, mask.m256i)};
#else
    md_i32x8_t val = {_mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(a.m256i), _mm256_castsi256_ps(b.m256i), _mm256_castsi256_ps(mask.m256i)))};
#endif
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_blend_i64x2(md_i64x2_t a, md_i64x2_t b, md_i64x2_t mask) {
    md_i64x2_t val = {_mm_blendv_epi8(a.m128i, b.m128i, mask.m128i)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_blend_i64x4(md_i64x4_t a, md_i64x4_t b, md_i64x4_t mask) {
#if __AVX2__
    md_i64x4_t val = {_mm256_blendv_epi8(a.m256i, b.m256i, mask.m256i)};
#else
    md_i64x4_t val = {_mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(a.m256i), _mm256_castsi256_pd(b.m256i), _mm256_castsi256_pd(mask.m256i)))};
#endif
    return val;
}

#define md_simd_load_f64x2 _mm_loadu_pd
#define md_simd_load_f64x4 _mm256_loadu_pd
#define md_simd_load_f32x4 _mm_loadu_ps
#define md_simd_load_f32x8 _mm256_loadu_ps

#define md_simd_store_f32x4 _mm_storeu_ps
#define md_simd_store_f64x2 _mm_storeu_pd
#define md_simd_store_f32x8 _mm256_storeu_ps
#define md_simd_store_f64x4 _mm256_storeu_pd

#define md_simd_set1_f32x4 _mm_set1_ps
#define md_simd_set1_f64x2 _mm_set1_pd
#define md_simd_set1_f32x8 _mm256_set1_ps
#define md_simd_set1_f64x4 _mm256_set1_pd

// @NOTE: Arguments are reversed compared to default intel intrinsics
#define md_simd_set_f32x4(X,Y,Z,W)          _mm_set_ps(W,Z,Y,X)
#define md_simd_set_f64x2(X,Y)              _mm_set_pd(Y,X)
#define md_simd_set_f32x8(X,Y,Z,W,Q,R,S,T)  _mm256_set_ps(T,S,R,Q,W,Z,Y,X)
#define md_simd_set_f64x4(X,Y,Z,W)          _mm256_set_pd(W,Z,Y,X)

#define md_simd_zero_f32x4 _mm_setzero_ps
#define md_simd_zero_f64x2 _mm_setzero_pd
#define md_simd_zero_f32x8 _mm256_setzero_ps
#define md_simd_zero_f64x4 _mm256_setzero_pd

#define md_simd_and_f32x4 _mm_and_ps
#define md_simd_and_f64x2 _mm_and_pd
#define md_simd_and_f32x8 _mm256_and_ps
#define md_simd_and_f64x4 _mm256_and_pd

#define md_simd_or_f32x4 _mm_or_ps
#define md_simd_or_f64x2 _mm_or_pd
#define md_simd_or_f32x8 _mm256_or_ps
#define md_simd_or_f64x4 _mm256_or_pd

#define md_simd_xor_f32x4 _mm_xor_ps
#define md_simd_xor_f64x2 _mm_xor_pd
#define md_simd_xor_f32x8 _mm256_xor_ps
#define md_simd_xor_f64x4 _mm256_xor_pd

MD_SIMD_INLINE md_f32x4_t md_simd_gather_f32x4(const float* base, const int* indices) {
#if __AVX2__
    __m128i i32x4 = _mm_lddqu_si128((const __m128i*)indices);
    return _mm_i32gather_ps(base, i32x4, 4);
#else
#   if MD_COMPILER_MSVC
    return _mm_set_ps(
        base[indices[0]],
        base[indices[1]],
        base[indices[2]],
        base[indices[3]]
    );
#   elif MD_COMPILER_GCC || MD_COMPILER_CLANG
    return _mm_set_ps(
        base[indices[0]],
        base[indices[1]],
        base[indices[2]],
        base[indices[3]]
    );
#   else
#       error "Unsupported compiler :<"
#   endif
#endif
}

MD_SIMD_INLINE md_f64x2_t md_simd_gather_f64x2(const double* base, const int* indices) {
#if __AVX2__
    __m128i i32x4 = _mm_lddqu_si128((const __m128i*)indices);
    return _mm_i32gather_pd(base, i32x4, 8);
#else
#   if MD_COMPILER_MSVC
    return _mm_set_pd(
        base[indices[0]],
        base[indices[1]]
    );
#   elif MD_COMPILER_GCC || MD_COMPILER_CLANG
    return _mm_set_pd(
        base[indices[0]],
        base[indices[1]]
    );
#   else
#       error "Unsupported compiler :<"
#   endif
#endif
}

MD_SIMD_INLINE md_f32x8_t md_simd_gather_f32x8(const float* base, const int* indices) {
#if __AVX2__
    __m256i i32x8 = _mm256_lddqu_si256((const __m256i*)indices);
    return _mm256_i32gather_ps(base, i32x8, 4);
#else
#   if MD_COMPILER_MSVC
    return _mm256_set_ps(
        base[indices[0]],
        base[indices[1]],
        base[indices[2]],
        base[indices[3]],
        base[indices[4]],
        base[indices[5]],
        base[indices[6]],
        base[indices[7]]
    );
#   elif MD_COMPILER_GCC || MD_COMPILER_CLANG
    return _mm256_set_ps(
        base[indices[0]],
        base[indices[1]],
        base[indices[2]],
        base[indices[3]],
        base[indices[4]],
        base[indices[5]],
        base[indices[6]],
        base[indices[7]]
    );
#   else
#       error "Unsupported compiler :<"
#   endif
#endif
}

MD_SIMD_INLINE md_f64x4_t md_simd_gather_f64x4(const double* base, const int* indices) {
#if __AVX2__
    __m128i i32x4 = _mm_lddqu_si128((const __m128i*)indices);
    return _mm256_i32gather_pd(base, i32x4, 8);
#else
#   if MD_COMPILER_MSVC
    return _mm256_set_pd(
        base[indices[0]],
        base[indices[1]],
        base[indices[2]],
        base[indices[3]]
    );
#   elif MD_COMPILER_GCC || MD_COMPILER_CLANG
    return _mm256_set_pd(
        base[indices[0]],
        base[indices[1]],
        base[indices[2]],
        base[indices[3]]
    );
#   else
#       error "Unsupported compiler :<"
#   endif
#endif
}

#if md_simd_width_f32 == 4
#define md_simd_gather_f32 md_simd_gather_f32x4
#elif md_simd_width_f32 == 8
#define md_simd_gather_f32 md_simd_gather_f32x8
#endif

#if md_simd_width_f64 == 2
#define md_simd_gather_f64 md_simd_gather_f64x2
#elif md_simd_width_f64 == 4
#define md_simd_gather_f64 md_simd_gather_f64x4
#endif

#define LOAD_STRIDED_128(ptr, offset, stride) _mm_load_ps((const float*)((const char*)ptr + offset * stride_in_bytes))

MD_SIMD_INLINE void md_simd_unpack_xyz_f32x4(md_f32x4_t* out_x, md_f32x4_t* out_y, md_f32x4_t* out_z, const float* in_xyz, size_t stride_in_bytes) {
    __m128  r0, r1, r2, r3;
    __m128  t0, t1, t2, t3;

    r0 = LOAD_STRIDED_128(in_xyz, 0, stride_in_bytes);
    r1 = LOAD_STRIDED_128(in_xyz, 1, stride_in_bytes);
    r2 = LOAD_STRIDED_128(in_xyz, 2, stride_in_bytes);
    r3 = LOAD_STRIDED_128(in_xyz, 3, stride_in_bytes);

    t0 = _mm_unpacklo_ps(r0,r1); // xxyy xxyy
    t1 = _mm_unpackhi_ps(r0,r1); // zzww zzww
    t2 = _mm_unpacklo_ps(r2,r3); // xxyy xxyy
    t3 = _mm_unpackhi_ps(r2,r3); // zzww zzww

    *out_x = _mm_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));  // xxxx xxxx
    *out_y = _mm_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));  // yyyy yyyy
    *out_z = _mm_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));  // zzzz zzzz
}

MD_SIMD_INLINE void md_simd_unpack_xyz_f32x8(md_f32x8_t* out_x, md_f32x8_t* out_y, md_f32x8_t* out_z, const float* in_xyz, size_t stride_in_bytes) {
    __m256 r0, r1, r2, r3;
    __m256 t0, t1, t2, t3;

    // @TODO: Try and implement this using 256-bit loads and then shuffle all the way.
    // It's a tradeoff between the number of loads issued and the port pressure on the generally few ports that are used for shuffle instructions.
    r0 = _mm256_insertf128_ps(_mm256_castps128_ps256(LOAD_STRIDED_128(in_xyz, 0, stride_in_bytes)), LOAD_STRIDED_128(in_xyz, 4, stride_in_bytes), 1);
    r1 = _mm256_insertf128_ps(_mm256_castps128_ps256(LOAD_STRIDED_128(in_xyz, 1, stride_in_bytes)), LOAD_STRIDED_128(in_xyz, 5, stride_in_bytes), 1);
    r2 = _mm256_insertf128_ps(_mm256_castps128_ps256(LOAD_STRIDED_128(in_xyz, 2, stride_in_bytes)), LOAD_STRIDED_128(in_xyz, 6, stride_in_bytes), 1);
    r3 = _mm256_insertf128_ps(_mm256_castps128_ps256(LOAD_STRIDED_128(in_xyz, 3, stride_in_bytes)), LOAD_STRIDED_128(in_xyz, 7, stride_in_bytes), 1);

    t0 = _mm256_unpacklo_ps(r0,r1); // xxyy xxyy
    t1 = _mm256_unpackhi_ps(r0,r1); // zzww zzww
    t2 = _mm256_unpacklo_ps(r2,r3); // xxyy xxyy
    t3 = _mm256_unpackhi_ps(r2,r3); // zzww zzww

    *out_x = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));  // xxxx xxxx
    *out_y = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));  // yyyy yyyy
    *out_z = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));  // zzzz zzzz
}

#undef LOAD_STRIDED_128

MD_SIMD_INLINE md_f32x4_t md_simd_and_not_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_andnot_ps(b, a); }
MD_SIMD_INLINE md_f32x8_t md_simd_and_not_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_andnot_ps(b, a); }
MD_SIMD_INLINE md_f64x2_t md_simd_and_not_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_andnot_pd(b, a); }
MD_SIMD_INLINE md_f64x4_t md_simd_and_not_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_andnot_pd(b, a); }

MD_SIMD_INLINE md_f32x4_t md_simd_abs_f32x4(md_f32x4_t a) { return _mm_and_ps(a, _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF))); }
MD_SIMD_INLINE md_f32x8_t md_simd_abs_f32x8(md_f32x8_t a) { return _mm256_and_ps(a, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF))); }
MD_SIMD_INLINE md_f64x2_t md_simd_abs_f64x2(md_f64x2_t a) { return _mm_and_pd(a, _mm_castsi128_pd(_mm_set1_epi64x(0x7FFFFFFFFFFFFFFF))); }
MD_SIMD_INLINE md_f64x4_t md_simd_abs_f64x4(md_f64x4_t a) { return _mm256_and_pd(a, _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF))); }

#define md_simd_add_f32x4 _mm_add_ps
#define md_simd_add_f64x2 _mm_add_pd
#define md_simd_add_f32x8 _mm256_add_ps
#define md_simd_add_f64x4 _mm256_add_pd

#define md_simd_sub_f32x4 _mm_sub_ps
#define md_simd_sub_f64x2 _mm_sub_pd
#define md_simd_sub_f32x8 _mm256_sub_ps
#define md_simd_sub_f64x4 _mm256_sub_pd

#define md_simd_min_f64x2 _mm_min_pd   
#define md_simd_min_f64x4 _mm256_min_pd
#define md_simd_min_f32x4 _mm_min_ps   
#define md_simd_min_f32x8 _mm256_min_ps

#define md_simd_max_f64x2 _mm_max_pd   
#define md_simd_max_f64x4 _mm256_max_pd
#define md_simd_max_f32x4 _mm_max_ps   
#define md_simd_max_f32x8 _mm256_max_ps

MD_SIMD_INLINE md_f32x4_t md_simd_cmp_gt_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_cmpgt_ps  (a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_cmp_gt_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_cmp_ps (a, b, _CMP_GT_OQ); }
MD_SIMD_INLINE md_f64x2_t md_simd_cmp_gt_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_cmpgt_pd  (a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_cmp_gt_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_cmp_pd (a, b, _CMP_GT_OQ); }

MD_SIMD_INLINE md_f32x4_t md_simd_cmp_ge_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_cmpge_ps  (a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_cmp_ge_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_cmp_ps (a, b, _CMP_GE_OQ); }
MD_SIMD_INLINE md_f64x2_t md_simd_cmp_ge_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_cmpge_pd  (a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_cmp_ge_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_cmp_pd (a, b, _CMP_GE_OQ); }

MD_SIMD_INLINE md_f32x4_t md_simd_cmp_lt_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_cmplt_ps  (a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_cmp_lt_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_cmp_ps (a, b, _CMP_LT_OQ); }
MD_SIMD_INLINE md_f64x2_t md_simd_cmp_lt_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_cmplt_pd  (a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_cmp_lt_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_cmp_pd (a, b, _CMP_LT_OQ); }

MD_SIMD_INLINE md_f32x4_t md_simd_cmp_le_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_cmple_ps  (a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_cmp_le_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_cmp_ps (a, b, _CMP_LE_OQ); }
MD_SIMD_INLINE md_f64x2_t md_simd_cmp_le_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_cmple_pd  (a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_cmp_le_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_cmp_pd (a, b, _CMP_LE_OQ); }

MD_SIMD_INLINE md_f32x4_t md_simd_cmp_eq_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_cmpeq_ps  (a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_cmp_eq_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_cmp_ps (a, b, _CMP_EQ_OQ); }
MD_SIMD_INLINE md_f64x2_t md_simd_cmp_eq_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_cmpeq_pd  (a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_cmp_eq_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_cmp_pd (a, b, _CMP_EQ_OQ); }

MD_SIMD_INLINE md_f32x4_t md_simd_cmp_ne_f32x4(md_f32x4_t a, md_f32x4_t b) { return _mm_cmpneq_ps (a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_cmp_ne_f32x8(md_f32x8_t a, md_f32x8_t b) { return _mm256_cmp_ps (a, b, _CMP_NEQ_OQ); }
MD_SIMD_INLINE md_f64x2_t md_simd_cmp_ne_f64x2(md_f64x2_t a, md_f64x2_t b) { return _mm_cmpneq_pd (a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_cmp_ne_f64x4(md_f64x4_t a, md_f64x4_t b) { return _mm256_cmp_pd (a, b, _CMP_NEQ_OQ); }

// FLOAT SPECIFIC
// for mul, there are operations defined for int as well, but the semantics differ: it multiplies the lanes as expected
// but the results are stored in two lanes rather than storing the truncated value in each lane.

#define md_simd_mul_f32x4 _mm_mul_ps
#define md_simd_mul_f32x8 _mm256_mul_ps
#define md_simd_mul_f64x2 _mm_mul_pd
#define md_simd_mul_f64x4 _mm256_mul_pd

#define md_simd_div_f32x4 _mm_div_ps
#define md_simd_div_f32x8 _mm256_div_ps
#define md_simd_div_f64x2 _mm_div_pd
#define md_simd_div_f64x4 _mm256_div_pd

#ifdef __FMA__
#define md_simd_fmadd_f32x4 _mm_fmadd_ps
#define md_simd_fmadd_f32x8 _mm256_fmadd_ps
#define md_simd_fmadd_f64x2 _mm_fmadd_pd
#define md_simd_fmadd_f64x4 _mm256_fmadd_pd
#else
MD_SIMD_INLINE md_f32x4_t md_simd_fmadd_f32x4(md_f32x4_t a, md_f32x4_t b, md_f32x4_t c) { return _mm_add_ps(_mm_mul_ps(a,b), c); }
MD_SIMD_INLINE md_f32x8_t md_simd_fmadd_f32x8(md_f32x8_t a, md_f32x8_t b, md_f32x8_t c) { return _mm256_add_ps(_mm256_mul_ps(a,b), c); }
MD_SIMD_INLINE md_f64x2_t md_simd_fmadd_f64x2(md_f64x2_t a, md_f64x2_t b, md_f64x2_t c) { return _mm_add_pd(_mm_mul_pd(a,b), c); }
MD_SIMD_INLINE md_f64x4_t md_simd_fmadd_f64x4(md_f64x4_t a, md_f64x4_t b, md_f64x4_t c) { return _mm256_add_pd(_mm256_mul_pd(a,b), c); }
#endif

MD_SIMD_INLINE md_f32x4_t md_simd_round_f32x4(md_f32x4_t a) { return _mm_round_ps     (a, _MM_FROUND_NO_EXC | _MM_FROUND_TO_NEAREST_INT); }
MD_SIMD_INLINE md_f32x8_t md_simd_round_f32x8(md_f32x8_t a) { return _mm256_round_ps  (a, _MM_FROUND_NO_EXC | _MM_FROUND_TO_NEAREST_INT); }
MD_SIMD_INLINE md_f64x2_t md_simd_round_f64x2(md_f64x2_t a) { return _mm_round_pd     (a, _MM_FROUND_NO_EXC | _MM_FROUND_TO_NEAREST_INT); }
MD_SIMD_INLINE md_f64x4_t md_simd_round_f64x4(md_f64x4_t a) { return _mm256_round_pd  (a, _MM_FROUND_NO_EXC | _MM_FROUND_TO_NEAREST_INT); }

MD_SIMD_INLINE md_f32x4_t md_simd_floor_f32x4(md_f32x4_t a) { return _mm_round_ps     (a, _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR); }
MD_SIMD_INLINE md_f32x8_t md_simd_floor_f32x8(md_f32x8_t a) { return _mm256_round_ps  (a, _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR); }
MD_SIMD_INLINE md_f64x2_t md_simd_floor_f64x2(md_f64x2_t a) { return _mm_round_pd     (a, _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR); }
MD_SIMD_INLINE md_f64x4_t md_simd_floor_f64x4(md_f64x4_t a) { return _mm256_round_pd  (a, _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR); }

MD_SIMD_INLINE md_f32x4_t md_simd_ceil_f32x4(md_f32x4_t a) { return _mm_round_ps      (a, _MM_FROUND_NO_EXC | _MM_FROUND_CEIL); }
MD_SIMD_INLINE md_f32x8_t md_simd_ceil_f32x8(md_f32x8_t a) { return _mm256_round_ps   (a, _MM_FROUND_NO_EXC | _MM_FROUND_CEIL); }
MD_SIMD_INLINE md_f64x2_t md_simd_ceil_f64x2(md_f64x2_t a) { return _mm_round_pd      (a, _MM_FROUND_NO_EXC | _MM_FROUND_CEIL); }
MD_SIMD_INLINE md_f64x4_t md_simd_ceil_f64x4(md_f64x4_t a) { return _mm256_round_pd   (a, _MM_FROUND_NO_EXC | _MM_FROUND_CEIL); }

MD_SIMD_INLINE md_f32x4_t md_simd_fract_f32x4(md_f32x4_t a) { return _mm_sub_ps(a,    _mm_round_ps(a,     _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR)); }
MD_SIMD_INLINE md_f32x8_t md_simd_fract_f32x8(md_f32x8_t a) { return _mm256_sub_ps(a, _mm256_round_ps(a,  _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR)); }
MD_SIMD_INLINE md_f64x2_t md_simd_fract_f64x2(md_f64x2_t a) { return _mm_sub_pd(a,    _mm_round_pd(a,     _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR)); }
MD_SIMD_INLINE md_f64x4_t md_simd_fract_f64x4(md_f64x4_t a) { return _mm256_sub_pd(a, _mm256_round_pd(a,  _MM_FROUND_NO_EXC | _MM_FROUND_FLOOR)); }

MD_SIMD_INLINE md_f32x4_t md_simd_sign_f32x4(md_f32x4_t a) { return _mm_xor_ps(   _mm_and_ps(a,    _mm_set1_ps(-0.0f)),    _mm_set1_ps(1.0f));    }
MD_SIMD_INLINE md_f32x8_t md_simd_sign_f32x8(md_f32x8_t a) { return _mm256_xor_ps(_mm256_and_ps(a, _mm256_set1_ps(-0.0f)), _mm256_set1_ps(1.0f)); }
MD_SIMD_INLINE md_f64x2_t md_simd_sign_f64x2(md_f64x2_t a) { return _mm_xor_pd(   _mm_and_pd(a,    _mm_set1_pd(-0.0)),     _mm_set1_pd(1.0));     }
MD_SIMD_INLINE md_f64x4_t md_simd_sign_f64x4(md_f64x4_t a) { return _mm256_xor_pd(_mm256_and_pd(a, _mm256_set1_pd(-0.0)),  _mm256_set1_pd(1.0));  }

#define md_simd_blend_f32x4 _mm_blendv_ps
#define md_simd_blend_f64x2 _mm_blendv_pd
#define md_simd_blend_f32x8 _mm256_blendv_ps
#define md_simd_blend_f64x4 _mm256_blendv_pd

#define md_simd_blend_mask_f32x4 _mm_blend_ps
#define md_simd_blend_mask_f64x2 _mm_blend_pd
#define md_simd_blend_mask_f32x8 _mm256_blend_ps
#define md_simd_blend_mask_f64x4 _mm256_blend_pd

#define md_simd_sqrt_f32x4 _mm_sqrt_ps
#define md_simd_sqrt_f64x2 _mm_sqrt_pd
#define md_simd_sqrt_f32x8 _mm256_sqrt_ps
#define md_simd_sqrt_f64x4 _mm256_sqrt_pd

MD_SIMD_INLINE float md_simd_hmin_f32x4(__m128 x) {
    __m128 a = _mm_min_ps(x, _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2)));
    __m128 b = _mm_min_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 0, 0, 1)));
    return _mm_cvtss_f32(_mm_min_ps(a,b));
}

MD_SIMD_INLINE float md_simd_hmax_f32x4(__m128 x) {
    __m128 a = _mm_max_ps(x, _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2)));
    __m128 b = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 0, 0, 1)));
    return _mm_cvtss_f32(_mm_max_ps(a,b));
}

MD_SIMD_INLINE double md_simd_hmin_f64x4(__m256d x) {
    __m256d a = _mm256_min_pd(x, _mm256_shuffle_pd(x, x, _MM_SHUFFLE(0, 0, 3, 2)));
    __m256d b = _mm256_min_pd(a, _mm256_shuffle_pd(a, a, _MM_SHUFFLE(0, 0, 0, 1)));
    return _mm256_cvtsd_f64(_mm256_min_pd(a,b));
}

MD_SIMD_INLINE double md_simd_hmax_f64x4(__m256d x) {
    __m256d a = _mm256_max_pd(x, _mm256_shuffle_pd(x, x, _MM_SHUFFLE(0, 0, 3, 2)));
    __m256d b = _mm256_max_pd(a, _mm256_shuffle_pd(a, a, _MM_SHUFFLE(0, 0, 0, 1)));
    return _mm256_cvtsd_f64(_mm256_max_pd(a,b));
}

// From here https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx
MD_SIMD_INLINE double md_simd_hsum_f64x4(__m256d x) {
    __m128d vlow   = _mm256_castpd256_pd128(x);   // low  128
    __m128d vhigh  = _mm256_extractf128_pd(x, 1); // high 128
    vlow           = _mm_add_pd(vlow, vhigh);     // reduce down to 128
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return           _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
MD_SIMD_INLINE float md_simd_hsum_f32x4(__m128 x) {
    //__m128 shuf = _mm_movehdup_ps(x);        // broadcast elements 3,1 to 2,0 (this instruction is SSE3 and we avoid it by using shuffle instead)
    __m128 shuf = _mm_shuffle_ps(x, x, _MM_SHUFFLE(3, 3, 1, 1));
    __m128 sums = _mm_add_ps(x, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

MD_SIMD_INLINE double md_simd_hsum_f64x2(md_f64x2_t x) { return _mm_cvtsd_f64(_mm_add_pd(x, _mm_shuffle_pd(x, x, _MM_SHUFFLE(0, 0, 0, 1)))); }
MD_SIMD_INLINE double md_simd_hmin_f64x2(md_f64x2_t x) { return _mm_cvtsd_f64(_mm_min_pd(x, _mm_shuffle_pd(x, x, _MM_SHUFFLE(0, 0, 0, 1)))); }
MD_SIMD_INLINE double md_simd_hmax_f64x2(md_f64x2_t x) { return _mm_cvtsd_f64(_mm_max_pd(x, _mm_shuffle_pd(x, x, _MM_SHUFFLE(0, 0, 0, 1)))); }

MD_SIMD_INLINE float md_simd_hsum_f32x8(md_f32x8_t x) { return md_simd_hsum_f32x4(_mm_add_ps(_mm256_castps256_ps128(x), _mm256_extractf128_ps(x, 0x1))); }
MD_SIMD_INLINE float md_simd_hmin_f32x8(md_f32x8_t x) { return md_simd_hmin_f32x4(_mm_min_ps(_mm256_castps256_ps128(x), _mm256_extractf128_ps(x, 0x1))); }
MD_SIMD_INLINE float md_simd_hmax_f32x8(md_f32x8_t x) { return md_simd_hmax_f32x4(_mm_max_ps(_mm256_castps256_ps128(x), _mm256_extractf128_ps(x, 0x1))); }

#define md_simd_movemask_f32x4 _mm_movemask_ps
#define md_simd_movemask_f32x8 _mm256_movemask_ps
#define md_simd_movemask_f64x2 _mm_movemask_pd
#define md_simd_movemask_f64x4 _mm256_movemask_pd

// CAST AND CONVERSIONS BETWEEN CORRESPONDING TYPES OF FLOAT AND INT
// NAMING CONVENTION GIVES THE SOURCE OPERAND TYPE

MD_SIMD_INLINE md_i32x4_t md_simd_cast_f32x4(md_f32x4_t a) {
    md_i32x4_t val = {_mm_castps_si128(a)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_cast_f32x8(md_f32x8_t a) {
    md_i32x8_t val = {_mm256_castps_si256(a)};
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_cast_f64x2(md_f64x2_t a) {
    md_i64x2_t val = {_mm_castpd_si128(a)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_cast_f64x4(md_f64x4_t a) {
    md_i64x4_t val = {_mm256_castpd_si256(a)};
    return val;
}

MD_SIMD_INLINE md_i32x4_t md_simd_convert_f32x4(md_f32x4_t a) {
    md_i32x4_t val = {_mm_cvtps_epi32(a)};
    return val;
}

MD_SIMD_INLINE md_i32x8_t md_simd_convert_f32x8(md_f32x8_t a) {
    md_i32x8_t val = {_mm256_cvtps_epi32(a)};
    return val;
}

MD_SIMD_INLINE md_i64x2_t md_simd_convert_f64x2(md_f64x2_t a) {
    md_i64x2_t val = {_mm_cvtpd_epi64(a)};
    return val;
}

MD_SIMD_INLINE md_i64x4_t md_simd_convert_f64x4(md_f64x4_t a) {
    md_i64x4_t val = {_mm256_cvtpd_epi64(a)};
    return val;
}

#define md_simd_cast_i32x4(X) _mm_castsi128_ps(X.m128i)
#define md_simd_cast_i32x8(X) _mm256_castsi256_ps(X.m256i)
#define md_simd_cast_i64x2(X) _mm_castsi128_pd(X.m128i)
#define md_simd_cast_i64x4(X) _mm256_castsi256_pd(X.m256i)

#define md_simd_convert_i32x4(X) _mm_cvtepi32_ps(X.m128i)
#define md_simd_convert_i32x8(X) _mm256_cvtepi32_ps(X.m256i)
#define md_simd_convert_i64x2(X) _mm_cvtepi64_pd(X.m128i)
#define md_simd_convert_i64x4(X) _mm256_cvtepi64_pd(X.m256i)

#endif

MD_SIMD_INLINE md_f32x4_t md_simd_deperiodize_f32x4(md_f32x4_t x, md_f32x4_t r, md_f32x4_t p) {
    md_f32x4_t d  = md_simd_sub_f32x4(x, r);
    md_f32x4_t dx = md_simd_div_f32x4(d, p);
    dx = md_simd_sub_f32x4(dx, md_simd_round_f32x4(dx));
    md_f32x4_t x_prim = md_simd_add_f32x4(r, md_simd_mul_f32x4(dx, p));
    return md_simd_blend_f32x4(x_prim, x, md_simd_cmp_eq_f32x4(p, md_simd_zero_f32x4()));
}

MD_SIMD_INLINE md_f32x8_t md_simd_deperiodize_f32x8(md_f32x8_t x, md_f32x8_t r, md_f32x8_t p) {
    md_f32x8_t d  = md_simd_sub_f32x8(x, r);
    md_f32x8_t dx = md_simd_div_f32x8(d, p);
    dx = md_simd_sub_f32x8(dx, md_simd_round_f32x8(dx));
    md_f32x8_t x_prim = md_simd_add_f32x8(r, md_simd_mul_f32x8(dx, p));
    return md_simd_blend_f32x8(x_prim, x, md_simd_cmp_eq_f32x8(p, md_simd_zero_f32x8()));
}

MD_SIMD_INLINE md_f64x2_t md_simd_deperiodize_f64x2(md_f64x2_t x, md_f64x2_t r, md_f64x2_t p) {
    md_f64x2_t d  = md_simd_sub_f64x2(x, r);
    md_f64x2_t dx = md_simd_div_f64x2(d, p);
    dx = md_simd_sub_f64x2(dx, md_simd_round_f64x2(dx));
    md_f64x2_t x_prim = md_simd_add_f64x2(r, md_simd_mul_f64x2(dx, p));
    return md_simd_blend_f64x2(x_prim, x, md_simd_cmp_eq_f64x2(p, md_simd_zero_f64x2()));
}

MD_SIMD_INLINE md_f64x4_t md_simd_deperiodize_f64x4(md_f64x4_t x, md_f64x4_t r, md_f64x4_t p) {
    md_f64x4_t d  = md_simd_sub_f64x4(x, r);
    md_f64x4_t dx = md_simd_div_f64x4(d, p);
    dx = md_simd_sub_f64x4(dx, md_simd_round_f64x4(dx));
    md_f64x4_t x_prim = md_simd_add_f64x4(r, md_simd_mul_f64x4(dx, p));
    return md_simd_blend_f64x4(x_prim, x, md_simd_cmp_eq_f64x4(p, md_simd_zero_f64x4()));
}

MD_SIMD_INLINE md_f32x4_t md_simd_minimum_image_f32x4(md_f32x4_t dx, md_f32x4_t p, md_f32x4_t rp) {
    return md_simd_sub_f32x4(dx, md_simd_mul_f32x4(p, md_simd_round_f32x4(md_simd_mul_f32x4(dx, rp))));
}

MD_SIMD_INLINE md_f32x8_t md_simd_minimum_image_f32x8(md_f32x8_t dx, md_f32x8_t p, md_f32x8_t rp) {
    return md_simd_sub_f32x8(dx, md_simd_mul_f32x8(p, md_simd_round_f32x8(md_simd_mul_f32x8(dx, rp))));
}

MD_SIMD_INLINE md_f64x2_t md_simd_minimum_image_f64x2(md_f64x2_t dx, md_f64x2_t p, md_f64x2_t rp) {
    return md_simd_sub_f64x2(dx, md_simd_mul_f64x2(p, md_simd_round_f64x2(md_simd_mul_f64x2(dx, rp))));
}

MD_SIMD_INLINE md_f64x4_t md_simd_minimum_image_f64x4(md_f64x4_t dx, md_f64x4_t p, md_f64x4_t rp) {
    return md_simd_sub_f64x4(dx, md_simd_mul_f64x4(p, md_simd_round_f64x4(md_simd_mul_f64x4(dx, rp))));
}

MD_SIMD_INLINE md_f32x4_t md_simd_cubic_spline_f32x4(md_f32x4_t p0, md_f32x4_t p1, md_f32x4_t p2, md_f32x4_t p3, md_f32x4_t t, md_f32x4_t s) {
    const md_f32x4_t t1 = t;
    const md_f32x4_t t2 = md_simd_mul_f32x4(t, t);
    const md_f32x4_t t3 = md_simd_mul_f32x4(t2, t);
    const md_f32x4_t v0 = md_simd_mul_f32x4(md_simd_sub_f32x4(p2, p0), s);
    const md_f32x4_t v1 = md_simd_mul_f32x4(md_simd_sub_f32x4(p3, p1), s);
    const md_f32x4_t x0 = md_simd_add_f32x4(md_simd_mul_f32x4(md_simd_set1_f32x4(2), md_simd_sub_f32x4(p1, p2)), md_simd_add_f32x4(v0, v1));
    const md_f32x4_t x1 = md_simd_sub_f32x4(md_simd_mul_f32x4(md_simd_set1_f32x4(3), md_simd_sub_f32x4(p2, p1)), md_simd_add_f32x4(md_simd_mul_f32x4(md_simd_set1_f32x4(2), v0), v1));
    const md_f32x4_t r0 = md_simd_add_f32x4(md_simd_mul_f32x4(x0, t3), md_simd_mul_f32x4(x1, t2));
    const md_f32x4_t r1 = md_simd_add_f32x4(md_simd_mul_f32x4(v0, t1), p1);
    return md_simd_add_f32x4(r0, r1);
}

MD_SIMD_INLINE md_f32x8_t md_simd_cubic_spline_f32x8(md_f32x8_t p0, md_f32x8_t p1, md_f32x8_t p2, md_f32x8_t p3, md_f32x8_t t, md_f32x8_t s) {
    const md_f32x8_t t2 = md_simd_mul_f32x8(t, t);
    const md_f32x8_t t3 = md_simd_mul_f32x8(t2, t);
    const md_f32x8_t v0 = md_simd_mul_f32x8(md_simd_sub_f32x8(p2, p0), s);
    const md_f32x8_t v1 = md_simd_mul_f32x8(md_simd_sub_f32x8(p3, p1), s);
    const md_f32x8_t x0 = md_simd_add_f32x8(md_simd_mul_f32x8(md_simd_set1_f32x8(2), md_simd_sub_f32x8(p1, p2)), md_simd_add_f32x8(v0, v1));
    const md_f32x8_t x1 = md_simd_sub_f32x8(md_simd_mul_f32x8(md_simd_set1_f32x8(3), md_simd_sub_f32x8(p2, p1)), md_simd_add_f32x8(md_simd_mul_f32x8(md_simd_set1_f32x8(2), v0), v1));
    const md_f32x8_t r0 = md_simd_add_f32x8(md_simd_mul_f32x8(x0, t3), md_simd_mul_f32x8(x1, t2));
    const md_f32x8_t r1 = md_simd_add_f32x8(md_simd_mul_f32x8(v0, t), p1);
    return md_simd_add_f32x8(r0, r1);
}

MD_SIMD_INLINE md_f64x2_t md_simd_cubic_spline_f64x2(md_f64x2_t p0, md_f64x2_t p1, md_f64x2_t p2, md_f64x2_t p3, md_f64x2_t t, md_f64x2_t s) {
    const md_f64x2_t t2 = md_simd_mul_f64x2(t, t);
    const md_f64x2_t t3 = md_simd_mul_f64x2(t2, t);
    const md_f64x2_t v0 = md_simd_mul_f64x2(md_simd_sub_f64x2(p2, p0), s);
    const md_f64x2_t v1 = md_simd_mul_f64x2(md_simd_sub_f64x2(p3, p1), s);
    const md_f64x2_t x0 = md_simd_add_f64x2(md_simd_mul_f64x2(md_simd_set1_f64x2(2), md_simd_sub_f64x2(p1, p2)), md_simd_add_f64x2(v0, v1));
    const md_f64x2_t x1 = md_simd_sub_f64x2(md_simd_mul_f64x2(md_simd_set1_f64x2(3), md_simd_sub_f64x2(p2, p1)), md_simd_add_f64x2(md_simd_mul_f64x2(md_simd_set1_f64x2(2), v0), v1));
    const md_f64x2_t r0 = md_simd_add_f64x2(md_simd_mul_f64x2(x0, t3), md_simd_mul_f64x2(x1, t2));
    const md_f64x2_t r1 = md_simd_add_f64x2(md_simd_mul_f64x2(v0, t), p1);
    return md_simd_add_f64x2(r0, r1);
}

MD_SIMD_INLINE md_f64x4_t md_simd_cubic_spline_f64x4(md_f64x4_t p0, md_f64x4_t p1, md_f64x4_t p2, md_f64x4_t p3, md_f64x4_t t, md_f64x4_t s) {
    const md_f64x4_t t2 = md_simd_mul_f64x4(t, t);
    const md_f64x4_t t3 = md_simd_mul_f64x4(t2, t);
    const md_f64x4_t v0 = md_simd_mul_f64x4(md_simd_sub_f64x4(p2, p0), s);
    const md_f64x4_t v1 = md_simd_mul_f64x4(md_simd_sub_f64x4(p3, p1), s);
    const md_f64x4_t x0 = md_simd_add_f64x4(md_simd_mul_f64x4(md_simd_set1_f64x4(2), md_simd_sub_f64x4(p1, p2)), md_simd_add_f64x4(v0, v1));
    const md_f64x4_t x1 = md_simd_sub_f64x4(md_simd_mul_f64x4(md_simd_set1_f64x4(3), md_simd_sub_f64x4(p2, p1)), md_simd_add_f64x4(md_simd_mul_f64x4(md_simd_set1_f64x4(2), v0), v1));
    const md_f64x4_t r0 = md_simd_add_f64x4(md_simd_mul_f64x4(x0, t3), md_simd_mul_f64x4(x1, t2));
    const md_f64x4_t r1 = md_simd_add_f64x4(md_simd_mul_f64x4(v0, t), p1);
    return md_simd_add_f64x4(r0, r1);
}

MD_SIMD_INLINE md_f32x4_t md_simd_step_f32x4(md_f32x4_t edge, md_f32x4_t x) {
    return md_simd_and_f32x4(md_simd_cmp_ge_f32x4(x, edge), md_simd_set1_f32x4(1.f));
}

MD_SIMD_INLINE md_f32x8_t md_simd_step_f32x8(md_f32x8_t edge, md_f32x8_t x) {
    return md_simd_and_f32x8(md_simd_cmp_ge_f32x8(x, edge), md_simd_set1_f32x8(1.f));
}

MD_SIMD_INLINE md_f64x2_t md_simd_step_f64x2(md_f64x2_t edge, md_f64x2_t x) {
    return md_simd_and_f64x2(md_simd_cmp_ge_f64x2(x, edge), md_simd_set1_f64x2(1.f));
}

MD_SIMD_INLINE md_f64x4_t md_simd_step_f64x4(md_f64x4_t edge, md_f64x4_t x) {
    return md_simd_and_f64x4(md_simd_cmp_ge_f64x4(x, edge), md_simd_set1_f64x4(1.f));
}

MD_SIMD_INLINE md_f32x4_t md_simd_lerp_f32x4(md_f32x4_t a, md_f32x4_t b, float t) {
    return md_simd_add_f32x4(md_simd_mul_f32x4(a, md_simd_set1_f32x4(1.0f - t)), md_simd_mul_f32x4(b, md_simd_set1_f32x4(t)));
}

MD_SIMD_INLINE md_f32x8_t md_simd_lerp_f32x8(md_f32x8_t a, md_f32x8_t b, float t) {
    return md_simd_add_f32x8(md_simd_mul_f32x8(a, md_simd_set1_f32x8(1.0f - t)), md_simd_mul_f32x8(b, md_simd_set1_f32x8(t)));
}

MD_SIMD_INLINE md_f64x2_t md_simd_lerp_f64x2(md_f64x2_t a, md_f64x2_t b, float t) {
    return md_simd_add_f64x2(md_simd_mul_f64x2(a, md_simd_set1_f64x2(1.0f - t)), md_simd_mul_f64x2(b, md_simd_set1_f64x2(t)));
}

MD_SIMD_INLINE md_f64x4_t md_simd_lerp_f64x4(md_f64x4_t a, md_f64x4_t b, float t) {
    return md_simd_add_f64x4(md_simd_mul_f64x4(a, md_simd_set1_f64x4(1.0f - t)), md_simd_mul_f64x4(b, md_simd_set1_f64x4(t)));
}

#ifdef __cplusplus
// C++ function overload
MD_SIMD_INLINE md_f32x4_t md_simd_add(md_f32x4_t a, md_f32x4_t b) { return md_simd_add_f32x4(a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_add(md_f32x8_t a, md_f32x8_t b) { return md_simd_add_f32x8(a, b); }
MD_SIMD_INLINE md_f64x2_t md_simd_add(md_f64x2_t a, md_f64x2_t b) { return md_simd_add_f64x2(a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_add(md_f64x4_t a, md_f64x4_t b) { return md_simd_add_f64x4(a, b); }

MD_SIMD_INLINE md_f32x4_t md_simd_sub(md_f32x4_t a, md_f32x4_t b) { return md_simd_sub_f32x4(a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_sub(md_f32x8_t a, md_f32x8_t b) { return md_simd_sub_f32x8(a, b); }
MD_SIMD_INLINE md_f64x2_t md_simd_sub(md_f64x2_t a, md_f64x2_t b) { return md_simd_sub_f64x2(a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_sub(md_f64x4_t a, md_f64x4_t b) { return md_simd_sub_f64x4(a, b); }

MD_SIMD_INLINE md_f32x4_t md_simd_mul(md_f32x4_t a, md_f32x4_t b) { return md_simd_mul_f32x4(a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_mul(md_f32x8_t a, md_f32x8_t b) { return md_simd_mul_f32x8(a, b); }
MD_SIMD_INLINE md_f64x2_t md_simd_mul(md_f64x2_t a, md_f64x2_t b) { return md_simd_mul_f64x2(a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_mul(md_f64x4_t a, md_f64x4_t b) { return md_simd_mul_f64x4(a, b); }

MD_SIMD_INLINE md_f32x4_t md_simd_div(md_f32x4_t a, md_f32x4_t b) { return md_simd_div_f32x4(a, b); }
MD_SIMD_INLINE md_f32x8_t md_simd_div(md_f32x8_t a, md_f32x8_t b) { return md_simd_div_f32x8(a, b); }
MD_SIMD_INLINE md_f64x2_t md_simd_div(md_f64x2_t a, md_f64x2_t b) { return md_simd_div_f64x2(a, b); }
MD_SIMD_INLINE md_f64x4_t md_simd_div(md_f64x4_t a, md_f64x4_t b) { return md_simd_div_f64x4(a, b); }

MD_SIMD_INLINE md_i32x4_t md_simd_cast(md_f32x4_t v) { return md_simd_cast_f32x4(v); }
MD_SIMD_INLINE md_i32x8_t md_simd_cast(md_f32x8_t v) { return md_simd_cast_f32x8(v); }
MD_SIMD_INLINE md_i64x2_t md_simd_cast(md_f64x2_t v) { return md_simd_cast_f64x2(v); }
MD_SIMD_INLINE md_i64x4_t md_simd_cast(md_f64x4_t v) { return md_simd_cast_f64x4(v); }

MD_SIMD_INLINE md_f32x4_t md_simd_cast(md_i32x4_t v) { return md_simd_cast_i32x4(v); }
MD_SIMD_INLINE md_f32x8_t md_simd_cast(md_i32x8_t v) { return md_simd_cast_i32x8(v); }
MD_SIMD_INLINE md_f64x2_t md_simd_cast(md_i64x2_t v) { return md_simd_cast_i64x2(v); }
MD_SIMD_INLINE md_f64x4_t md_simd_cast(md_i64x4_t v) { return md_simd_cast_i64x4(v); }

MD_SIMD_INLINE md_i32x4_t md_simd_convert(md_f32x4_t v) { return md_simd_convert_f32x4(v); }
MD_SIMD_INLINE md_i32x8_t md_simd_convert(md_f32x8_t v) { return md_simd_convert_f32x8(v); }
MD_SIMD_INLINE md_i64x2_t md_simd_convert(md_f64x2_t v) { return md_simd_convert_f64x2(v); }
MD_SIMD_INLINE md_i64x4_t md_simd_convert(md_f64x4_t v) { return md_simd_convert_f64x4(v); }

MD_SIMD_INLINE md_f32x4_t md_simd_convert(md_i32x4_t v) { return md_simd_convert_i32x4(v); }
MD_SIMD_INLINE md_f32x8_t md_simd_convert(md_i32x8_t v) { return md_simd_convert_i32x8(v); }
MD_SIMD_INLINE md_f64x2_t md_simd_convert(md_i64x2_t v) { return md_simd_convert_i64x2(v); }
MD_SIMD_INLINE md_f64x4_t md_simd_convert(md_i64x4_t v) { return md_simd_convert_i64x4(v); }

// @TODO: Complete this

#else
// C11 generics

#define md_simd_store(a, b) _Generic((b),   \
            md_i32x4_t : md_simd_store_i32x4, \
            md_i32x8_t : md_simd_store_i32x8, \
            md_i64x2_t : md_simd_store_i64x2, \
            md_i64x4_t : md_simd_store_i64x4, \
            const md_i32x4_t : md_simd_store_i32x4, \
            const md_i32x8_t : md_simd_store_i32x8, \
            const md_i64x2_t : md_simd_store_i64x2, \
            const md_i64x4_t : md_simd_store_i64x4, \
            md_f32x4_t : md_simd_store_f32x4, \
            md_f32x8_t : md_simd_store_f32x8, \
            md_f64x2_t : md_simd_store_f64x2, \
            md_f64x4_t : md_simd_store_f64x4, \
            const md_f32x4_t : md_simd_store_f32x4, \
            const md_f32x8_t : md_simd_store_f32x8, \
            const md_f64x2_t : md_simd_store_f64x2, \
            const md_f64x4_t : md_simd_store_f64x4)(a,b)

#define md_simd_and(a, b) _Generic((a),     \
            md_i32x4_t : md_simd_and_i32x4, \
            md_i32x8_t : md_simd_and_i32x8, \
            md_i64x2_t : md_simd_and_i64x2, \
            md_i64x4_t : md_simd_and_i64x4, \
            const md_i32x4_t : md_simd_and_i32x4, \
            const md_i32x8_t : md_simd_and_i32x8, \
            const md_i64x2_t : md_simd_and_i64x2, \
            const md_i64x4_t : md_simd_and_i64x4, \
            md_f32x4_t : md_simd_and_f32x4, \
            md_f32x8_t : md_simd_and_f32x8, \
            md_f64x2_t : md_simd_and_f64x2, \
            md_f64x4_t : md_simd_and_f64x4, \
            const md_f32x4_t : md_simd_and_f32x4, \
            const md_f32x8_t : md_simd_and_f32x8, \
            const md_f64x2_t : md_simd_and_f64x2, \
            const md_f64x4_t : md_simd_and_f64x4)(a,b)

#define md_simd_or(a, b) _Generic((a),      \
            md_i32x4_t : md_simd_or_i32x4,  \
            md_i32x8_t : md_simd_or_i32x8,  \
            md_i64x2_t : md_simd_or_i64x2,  \
            md_i64x4_t : md_simd_or_i64x4,  \
            const md_i32x4_t : md_simd_or_i32x4,  \
            const md_i32x8_t : md_simd_or_i32x8,  \
            const md_i64x2_t : md_simd_or_i64x2,  \
            const md_i64x4_t : md_simd_or_i64x4,  \
            md_f32x4_t : md_simd_or_f32x4,  \
            md_f32x8_t : md_simd_or_f32x8,  \
            md_f64x2_t : md_simd_or_f64x2,  \
            md_f64x4_t : md_simd_or_f64x4,  \
            const md_f32x4_t : md_simd_or_f32x4,  \
            const md_f32x8_t : md_simd_or_f32x8,  \
            const md_f64x2_t : md_simd_or_f64x2,  \
            const md_f64x4_t : md_simd_or_f64x4)(a,b)

#define md_simd_xor(a, b) _Generic((a),     \
            md_i32x4_t : md_simd_xor_i32x4, \
            md_i32x8_t : md_simd_xor_i32x8, \
            md_i64x2_t : md_simd_xor_i64x2, \
            md_i64x4_t : md_simd_xor_i64x4, \
            const md_i32x4_t : md_simd_xor_i32x4, \
            const md_i32x8_t : md_simd_xor_i32x8, \
            const md_i64x2_t : md_simd_xor_i64x2, \
            const md_i64x4_t : md_simd_xor_i64x4, \
            md_f32x4_t : md_simd_xor_f32x4, \
            md_f32x8_t : md_simd_xor_f32x8, \
            md_f64x2_t : md_simd_xor_f64x2, \
            md_f64x4_t : md_simd_xor_f64x4, \
            const md_f32x4_t : md_simd_xor_f32x4, \
            const md_f32x8_t : md_simd_xor_f32x8, \
            const md_f64x2_t : md_simd_xor_f64x2, \
            const md_f64x4_t : md_simd_xor_f64x4)(a,b)

#define md_simd_and_not(a, b) _Generic((a),     \
            md_i32x4_t : md_simd_and_not_i32x4, \
            md_i32x8_t : md_simd_and_not_i32x8, \
            md_i64x2_t : md_simd_and_not_i64x2, \
            md_i64x4_t : md_simd_and_not_i64x4, \
            const md_i32x4_t : md_simd_and_not_i32x4, \
            const md_i32x8_t : md_simd_and_not_i32x8, \
            const md_i64x2_t : md_simd_and_not_i64x2, \
            const md_i64x4_t : md_simd_and_not_i64x4, \
            md_f32x4_t : md_simd_and_not_f32x4, \
            md_f32x8_t : md_simd_and_not_f32x8, \
            md_f64x2_t : md_simd_and_not_f64x2, \
            md_f64x4_t : md_simd_and_not_f64x4, \
            const md_f32x4_t : md_simd_and_not_f32x4, \
            const md_f32x8_t : md_simd_and_not_f32x8, \
            const md_f64x2_t : md_simd_and_not_f64x2, \
            const md_f64x4_t : md_simd_and_not_f64x4)(a,b)

#define md_simd_not(x) _Generic((x),        \
            md_i32x4_t : md_simd_not_i32x4, \
            md_i32x8_t : md_simd_not_i32x8, \
            md_i64x2_t : md_simd_not_i64x2, \
            md_i64x4_t : md_simd_not_i64x4, \
            const md_i32x4_t : md_simd_not_i32x4, \
            const md_i32x8_t : md_simd_not_i32x8, \
            const md_i64x2_t : md_simd_not_i64x2, \
            const md_i64x4_t : md_simd_not_i64x4)(x)

#define md_simd_add(a, b) _Generic((a),     \
            md_f32x4_t : md_simd_add_f32x4, \
            md_f32x8_t : md_simd_add_f32x8, \
            md_f64x2_t : md_simd_add_f64x2, \
            md_f64x4_t : md_simd_add_f64x4, \
            const md_f32x4_t : md_simd_add_f32x4, \
            const md_f32x8_t : md_simd_add_f32x8, \
            const md_f64x2_t : md_simd_add_f64x2, \
            const md_f64x4_t : md_simd_add_f64x4, \
            md_i32x4_t : md_simd_add_i32x4, \
            md_i32x8_t : md_simd_add_i32x8, \
            md_i64x2_t : md_simd_add_i64x2, \
            md_i64x4_t : md_simd_add_i64x4, \
            const md_i32x4_t : md_simd_add_i32x4, \
            const md_i32x8_t : md_simd_add_i32x8, \
            const md_i64x2_t : md_simd_add_i64x2, \
            const md_i64x4_t : md_simd_add_i64x4)(a, b)
            
#define md_simd_sub(a, b) _Generic((a),     \
            md_f32x4_t : md_simd_sub_f32x4, \
            md_f32x8_t : md_simd_sub_f32x8, \
            md_f64x2_t : md_simd_sub_f64x2, \
            md_f64x4_t : md_simd_sub_f64x4, \
            const md_f32x4_t : md_simd_sub_f32x4, \
            const md_f32x8_t : md_simd_sub_f32x8, \
            const md_f64x2_t : md_simd_sub_f64x2, \
            const md_f64x4_t : md_simd_sub_f64x4, \
            md_i32x4_t : md_simd_sub_i32x4, \
            md_i32x8_t : md_simd_sub_i32x8, \
            md_i64x2_t : md_simd_sub_i64x2, \
            md_i64x4_t : md_simd_sub_i64x4, \
            const md_i32x4_t : md_simd_sub_i32x4, \
            const md_i32x8_t : md_simd_sub_i32x8, \
            const md_i64x2_t : md_simd_sub_i64x2, \
            const md_i64x4_t : md_simd_sub_i64x4)(a, b)

#define md_simd_mul(a, b) _Generic((a),     \
            md_f32x4_t : md_simd_mul_f32x4, \
            md_f32x8_t : md_simd_mul_f32x8, \
            md_f64x2_t : md_simd_mul_f64x2, \
            md_f64x4_t : md_simd_mul_f64x4, \
            const md_f32x4_t : md_simd_mul_f32x4, \
            const md_f32x8_t : md_simd_mul_f32x8, \
            const md_f64x2_t : md_simd_mul_f64x2, \
            const md_f64x4_t : md_simd_mul_f64x4)(a, b)

#define md_simd_div(a, b) _Generic((a),     \
            md_f32x4_t : md_simd_div_f32x4, \
            md_f32x8_t : md_simd_div_f32x8, \
            md_f64x2_t : md_simd_div_f64x2, \
            md_f64x4_t : md_simd_div_f64x4, \
            const md_f32x4_t : md_simd_div_f32x4, \
            const md_f32x8_t : md_simd_div_f32x8, \
            const md_f64x2_t : md_simd_div_f64x2, \
            const md_f64x4_t : md_simd_div_f64x4)(a, b)

#define md_simd_fmadd(a, b, c) _Generic((a),       \
            md_f32x4_t : md_simd_fmadd_f32x4,   \
            md_f32x8_t : md_simd_fmadd_f32x8,   \
            md_f64x2_t : md_simd_fmadd_f64x2,   \
            md_f64x4_t : md_simd_fmadd_f64x4,   \
            const md_f32x4_t : md_simd_fmadd_f32x4,   \
            const md_f32x8_t : md_simd_fmadd_f32x8,   \
            const md_f64x2_t : md_simd_fmadd_f64x2,   \
            const md_f64x4_t : md_simd_fmadd_f64x4)(a, b, c)

#define md_simd_abs(x) _Generic((x),        \
            md_f32x4_t : md_simd_abs_f32x4, \
            md_f32x8_t : md_simd_abs_f32x8, \
            md_f64x2_t : md_simd_abs_f64x2, \
            md_f64x4_t : md_simd_abs_f64x4, \
            const md_f32x4_t : md_simd_abs_f32x4, \
            const md_f32x8_t : md_simd_abs_f32x8, \
            const md_f64x2_t : md_simd_abs_f64x2, \
            const md_f64x4_t : md_simd_abs_f64x4, \
            md_i32x4_t : md_simd_abs_i32x4, \
            md_i32x8_t : md_simd_abs_i32x8, \
            md_i64x2_t : md_simd_abs_i64x2, \
            md_i64x4_t : md_simd_abs_i64x4, \
            const md_i32x4_t : md_simd_abs_i32x4, \
            const md_i32x8_t : md_simd_abs_i32x8, \
            const md_i64x2_t : md_simd_abs_i64x2, \
            const md_i64x4_t : md_simd_abs_i64x4)(x)

#define md_simd_min(a,b) _Generic((a),      \
            md_f32x4_t : md_simd_min_f32x4, \
            md_f32x8_t : md_simd_min_f32x8, \
            md_f64x2_t : md_simd_min_f64x2, \
            md_f64x4_t : md_simd_min_f64x4, \
            const md_f32x4_t : md_simd_min_f32x4, \
            const md_f32x8_t : md_simd_min_f32x8, \
            const md_f64x2_t : md_simd_min_f64x2, \
            const md_f64x4_t : md_simd_min_f64x4, \
            md_i32x4_t : md_simd_min_i32x4, \
            md_i32x8_t : md_simd_min_i32x8, \
            md_i64x2_t : md_simd_min_i64x2, \
            md_i64x4_t : md_simd_min_i64x4, \
            const md_i32x4_t : md_simd_min_i32x4, \
            const md_i32x8_t : md_simd_min_i32x8, \
            const md_i64x2_t : md_simd_min_i64x2, \
            const md_i64x4_t : md_simd_min_i64x4)(a,b)

#define md_simd_max(a,b) _Generic((a),      \
            md_f32x4_t : md_simd_max_f32x4, \
            md_f32x8_t : md_simd_max_f32x8, \
            md_f64x2_t : md_simd_max_f64x2, \
            md_f64x4_t : md_simd_max_f64x4, \
            const md_f32x4_t : md_simd_max_f32x4, \
            const md_f32x8_t : md_simd_max_f32x8, \
            const md_f64x2_t : md_simd_max_f64x2, \
            const md_f64x4_t : md_simd_max_f64x4, \
            md_i32x4_t : md_simd_max_i32x4, \
            md_i32x8_t : md_simd_max_i32x8, \
            md_i64x2_t : md_simd_max_i64x2, \
            md_i64x4_t : md_simd_max_i64x4, \
            const md_i32x4_t : md_simd_max_i32x4, \
            const md_i32x8_t : md_simd_max_i32x8, \
            const md_i64x2_t : md_simd_max_i64x2, \
            const md_i64x4_t : md_simd_max_i64x4)(a,b)

#define md_simd_cmp_gt(a,b) _Generic((a),       \
            md_f32x4_t : md_simd_cmp_gt_f32x4,  \
            md_f32x8_t : md_simd_cmp_gt_f32x8,  \
            md_f64x2_t : md_simd_cmp_gt_f64x2,  \
            md_f64x4_t : md_simd_cmp_gt_f64x4,  \
            const md_f32x4_t : md_simd_cmp_gt_f32x4,  \
            const md_f32x8_t : md_simd_cmp_gt_f32x8,  \
            const md_f64x2_t : md_simd_cmp_gt_f64x2,  \
            const md_f64x4_t : md_simd_cmp_gt_f64x4,  \
            md_i32x4_t : md_simd_cmp_gt_i32x4,  \
            md_i32x8_t : md_simd_cmp_gt_i32x8,  \
            md_i64x2_t : md_simd_cmp_gt_i64x2,  \
            md_i64x4_t : md_simd_cmp_gt_i64x4,  \
            const md_i32x4_t : md_simd_cmp_gt_i32x4,  \
            const md_i32x8_t : md_simd_cmp_gt_i32x8,  \
            const md_i64x2_t : md_simd_cmp_gt_i64x2,  \
            const md_i64x4_t : md_simd_cmp_gt_i64x4)(a,b)

#define md_simd_cmp_lt(a,b) _Generic((a),         \
            md_f32x4_t : md_simd_cmp_lt_f32x4,  \
            md_f32x8_t : md_simd_cmp_lt_f32x8,  \
            md_f64x2_t : md_simd_cmp_lt_f64x2,  \
            md_f64x4_t : md_simd_cmp_lt_f64x4,  \
            const md_f32x4_t : md_simd_cmp_lt_f32x4,  \
            const md_f32x8_t : md_simd_cmp_lt_f32x8,  \
            const md_f64x2_t : md_simd_cmp_lt_f64x2,  \
            const md_f64x4_t : md_simd_cmp_lt_f64x4,  \
            md_i32x4_t : md_simd_cmp_lt_i32x4,  \
            md_i32x8_t : md_simd_cmp_lt_i32x8,  \
            md_i64x2_t : md_simd_cmp_lt_i64x2,  \
            md_i64x4_t : md_simd_cmp_lt_i64x4,  \
            const md_i32x4_t : md_simd_cmp_lt_i32x4,  \
            const md_i32x8_t : md_simd_cmp_lt_i32x8,  \
            const md_i64x2_t : md_simd_cmp_lt_i64x2,  \
            const md_i64x4_t : md_simd_cmp_lt_i64x4)(a,b)

#define md_simd_cmp_eq(a,b) _Generic((a),       \
            md_f32x4_t : md_simd_cmp_eq_f32x4,  \
            md_f32x8_t : md_simd_cmp_eq_f32x8,  \
            md_f64x2_t : md_simd_cmp_eq_f64x2,  \
            md_f64x4_t : md_simd_cmp_eq_f64x4,  \
            const md_f32x4_t : md_simd_cmp_eq_f32x4,  \
            const md_f32x8_t : md_simd_cmp_eq_f32x8,  \
            const md_f64x2_t : md_simd_cmp_eq_f64x2,  \
            const md_f64x4_t : md_simd_cmp_eq_f64x4,  \
            md_i32x4_t : md_simd_cmp_eq_i32x4,  \
            md_i32x8_t : md_simd_cmp_eq_i32x8,  \
            md_i64x2_t : md_simd_cmp_eq_i64x2,  \
            md_i64x4_t : md_simd_cmp_eq_i64x4,  \
            const md_i32x4_t : md_simd_cmp_eq_i32x4,  \
            const md_i32x8_t : md_simd_cmp_eq_i32x8,  \
            const md_i64x2_t : md_simd_cmp_eq_i64x2,  \
            const md_i64x4_t : md_simd_cmp_eq_i64x4)(a,b)

#define md_simd_cmp_ne(a,b) _Generic((a),       \
            md_f32x4_t : md_simd_cmp_ne_f32x4,  \
            md_f32x8_t : md_simd_cmp_ne_f32x8,  \
            md_f64x2_t : md_simd_cmp_ne_f64x2,  \
            md_f64x4_t : md_simd_cmp_ne_f64x4,  \
            const md_f32x4_t : md_simd_cmp_ne_f32x4,  \
            const md_f32x8_t : md_simd_cmp_ne_f32x8,  \
            const md_f64x2_t : md_simd_cmp_ne_f64x2,  \
            const md_f64x4_t : md_simd_cmp_ne_f64x4,  \
            md_i32x4_t : md_simd_cmp_ne_i32x4,  \
            md_i32x8_t : md_simd_cmp_ne_i32x8,  \
            md_i64x2_t : md_simd_cmp_ne_i64x2,  \
            md_i64x4_t : md_simd_cmp_ne_i64x4,  \
            const md_i32x4_t : md_simd_cmp_ne_i32x4,  \
            const md_i32x8_t : md_simd_cmp_ne_i32x8,  \
            const md_i64x2_t : md_simd_cmp_ne_i64x2,  \
            const md_i64x4_t : md_simd_cmp_ne_i64x4)(a,b)

#define md_simd_cmp_le(a,b) _Generic((a),       \
            md_f32x4_t : md_simd_cmp_le_f32x4,  \
            md_f32x8_t : md_simd_cmp_le_f32x8,  \
            md_f64x2_t : md_simd_cmp_le_f64x2,  \
            md_f64x4_t : md_simd_cmp_le_f64x4,  \
            const md_f32x4_t : md_simd_cmp_le_f32x4,  \
            const md_f32x8_t : md_simd_cmp_le_f32x8,  \
            const md_f64x2_t : md_simd_cmp_le_f64x2,  \
            const md_f64x4_t : md_simd_cmp_le_f64x4)(a,b)

#define md_simd_cmp_le(a,b) _Generic((a),       \
            md_f32x4_t : md_simd_cmp_le_f32x4,  \
            md_f32x8_t : md_simd_cmp_le_f32x8,  \
            md_f64x2_t : md_simd_cmp_le_f64x2,  \
            md_f64x4_t : md_simd_cmp_le_f64x4,  \
            const md_f32x4_t : md_simd_cmp_le_f32x4,  \
            const md_f32x8_t : md_simd_cmp_le_f32x8,  \
            const md_f64x2_t : md_simd_cmp_le_f64x2,  \
            const md_f64x4_t : md_simd_cmp_le_f64x4)(a,b)

#define md_simd_round(x) _Generic((x),          \
            md_f32x4_t : md_simd_round_f32x4,   \
            md_f32x8_t : md_simd_round_f32x8,   \
            md_f64x2_t : md_simd_round_f64x2,   \
            md_f64x4_t : md_simd_round_f64x4,   \
            const md_f32x4_t : md_simd_round_f32x4,   \
            const md_f32x8_t : md_simd_round_f32x8,   \
            const md_f64x2_t : md_simd_round_f64x2,   \
            const md_f64x4_t : md_simd_round_f64x4)(x)

#define md_simd_floor(x) _Generic((x),          \
            md_f32x4_t : md_simd_floor_f32x4,   \
            md_f32x8_t : md_simd_floor_f32x8,   \
            md_f64x2_t : md_simd_floor_f64x2,   \
            md_f64x4_t : md_simd_floor_f64x4,   \
            const md_f32x4_t : md_simd_floor_f32x4,   \
            const md_f32x8_t : md_simd_floor_f32x8,   \
            const md_f64x2_t : md_simd_floor_f64x2,   \
            const md_f64x4_t : md_simd_floor_f64x4)(x)

#define md_simd_ceil(x) _Generic((x),           \
            md_f32x4_t : md_simd_ceil_f32x4,    \
            md_f32x8_t : md_simd_ceil_f32x8,    \
            md_f64x2_t : md_simd_ceil_f64x2,    \
            md_f64x4_t : md_simd_ceil_f64x4,    \
            const md_f32x4_t : md_simd_ceil_f32x4,    \
            const md_f32x8_t : md_simd_ceil_f32x8,    \
            const md_f64x2_t : md_simd_ceil_f64x2,    \
            const md_f64x4_t : md_simd_ceil_f64x4)(x)

#define md_simd_fract(x) _Generic((x),          \
            md_f32x4_t : md_simd_fract_f32x4,   \
            md_f32x8_t : md_simd_fract_f32x8,   \
            md_f64x2_t : md_simd_fract_f64x2,   \
            md_f64x4_t : md_simd_fract_f64x4,   \
            const md_f32x4_t : md_simd_fract_f32x4,   \
            const md_f32x8_t : md_simd_fract_f32x8,   \
            const md_f64x2_t : md_simd_fract_f64x2,   \
            const md_f64x4_t : md_simd_fract_f64x4)(x)

#define md_simd_sign(x) _Generic((x),           \
            md_f32x4_t : md_simd_sign_f32x4,    \
            md_f32x8_t : md_simd_sign_f32x8,    \
            md_f64x2_t : md_simd_sign_f64x2,    \
            md_f64x4_t : md_simd_sign_f64x4)(x)

#define md_simd_sqrt(x) _Generic((x),           \
            md_f32x4_t : md_simd_sqrt_f32x4,    \
            md_f32x8_t : md_simd_sqrt_f32x8,    \
            md_f64x2_t : md_simd_sqrt_f64x2,    \
            md_f64x4_t : md_simd_sqrt_f64x4,    \
            const md_f32x4_t : md_simd_sqrt_f32x4,    \
            const md_f32x8_t : md_simd_sqrt_f32x8,    \
            const md_f64x2_t : md_simd_sqrt_f64x2,    \
            const md_f64x4_t : md_simd_sqrt_f64x4)(x)

#define md_simd_blend(a,b,mask) _Generic((a),   \
            md_f32x4_t : md_simd_blend_f32x4,   \
            md_f32x8_t : md_simd_blend_f32x8,   \
            md_f64x2_t : md_simd_blend_f64x2,   \
            md_f64x4_t : md_simd_blend_f64x4,   \
            const md_f32x4_t : md_simd_blend_f32x4,   \
            const md_f32x8_t : md_simd_blend_f32x8,   \
            const md_f64x2_t : md_simd_blend_f64x2,   \
            const md_f64x4_t : md_simd_blend_f64x4,   \
            md_i32x4_t : md_simd_blend_i32x4,   \
            md_i32x8_t : md_simd_blend_i32x8,   \
            md_i64x2_t : md_simd_blend_i64x2,   \
            md_i64x4_t : md_simd_blend_i64x4,   \
            const md_i32x4_t : md_simd_blend_i32x4,   \
            const md_i32x8_t : md_simd_blend_i32x8,   \
            const md_i64x2_t : md_simd_blend_i64x2,   \
            const md_i64x4_t : md_simd_blend_i64x4)(a,b,mask)

#define md_simd_movemask(a) _Generic((a),       \
            md_f32x4_t : md_simd_movemask_f32x4,\
            md_f32x8_t : md_simd_movemask_f32x8,\
            md_f64x2_t : md_simd_movemask_f64x2,\
            md_f64x4_t : md_simd_movemask_f64x4,\
            const md_f32x4_t : md_simd_movemask_f32x4,\
            const md_f32x8_t : md_simd_movemask_f32x8,\
            const md_f64x2_t : md_simd_movemask_f64x2,\
            const md_f64x4_t : md_simd_movemask_f64x4)(a)

#define md_simd_hmin(x) _Generic((x),           \
            md_f32x4_t : md_simd_hmin_f32x4,    \
            md_f32x8_t : md_simd_hmin_f32x8,    \
            md_f64x2_t : md_simd_hmin_f64x2,    \
            md_f64x4_t : md_simd_hmin_f64x4,    \
            const md_f32x4_t : md_simd_hmin_f32x4,    \
            const md_f32x8_t : md_simd_hmin_f32x8,    \
            const md_f64x2_t : md_simd_hmin_f64x2,    \
            const md_f64x4_t : md_simd_hmin_f64x4)(x)

#define md_simd_hmax(x) _Generic((x),           \
            md_f32x4_t : md_simd_hmax_f32x4,    \
            md_f32x8_t : md_simd_hmax_f32x8,    \
            md_f64x2_t : md_simd_hmax_f64x2,    \
            md_f64x4_t : md_simd_hmax_f64x4,    \
            const md_f32x4_t : md_simd_hmax_f32x4,    \
            const md_f32x8_t : md_simd_hmax_f32x8,    \
            const md_f64x2_t : md_simd_hmax_f64x2,    \
            const md_f64x4_t : md_simd_hmax_f64x4)(x)

#define md_simd_hsum(x) _Generic((x),           \
            md_f32x4_t : md_simd_hsum_f32x4,    \
            md_f32x8_t : md_simd_hsum_f32x8,    \
            md_f64x2_t : md_simd_hsum_f64x2,    \
            md_f64x4_t : md_simd_hsum_f64x4,    \
            const md_f32x4_t : md_simd_hsum_f32x4,    \
            const md_f32x8_t : md_simd_hsum_f32x8,    \
            const md_f64x2_t : md_simd_hsum_f64x2,    \
            const md_f64x4_t : md_simd_hsum_f64x4)(x)

#define md_simd_deperiodize(x, r, p) _Generic((x),  \
            md_f32x4_t : md_simd_deperiodize_f32x4, \
            md_f32x8_t : md_simd_deperiodize_f32x8, \
            md_f64x2_t : md_simd_deperiodize_f64x2, \
            md_f64x4_t : md_simd_deperiodize_f64x4, \
            const md_f32x4_t : md_simd_deperiodize_f32x4, \
            const md_f32x8_t : md_simd_deperiodize_f32x8, \
            const md_f64x2_t : md_simd_deperiodize_f64x2, \
            const md_f64x4_t : md_simd_deperiodize_f64x4)(x, r, p)

#define md_simd_minimum_image(dx, p, rp) _Generic((dx),  \
            md_f32x4_t : md_simd_minimum_image_f32x4, \
            md_f32x8_t : md_simd_minimum_image_f32x8, \
            md_f64x2_t : md_simd_minimum_image_f64x2, \
            md_f64x4_t : md_simd_minimum_image_f64x4, \
            const md_f32x4_t : md_simd_minimum_image_f32x4, \
            const md_f32x8_t : md_simd_minimum_image_f32x8, \
            const md_f64x2_t : md_simd_minimum_image_f64x2, \
            const md_f64x4_t : md_simd_minimum_image_f64x4)(dx, p, rp)

#define md_simd_step(edge, x) _Generic((x),  \
            md_f32x4_t : md_simd_step_f32x4, \
            md_f32x8_t : md_simd_step_f32x8, \
            md_f64x2_t : md_simd_step_f64x2, \
            md_f64x4_t : md_simd_step_f64x4, \
            const md_f32x4_t : md_simd_step_f32x4, \
            const md_f32x8_t : md_simd_step_f32x8, \
            const md_f64x2_t : md_simd_step_f64x2, \
            const md_f64x4_t : md_simd_step_f64x4)(edge, x)

#define md_simd_lerp(a, b, t) _Generic((a),  \
            md_f32x4_t : md_simd_lerp_f32x4, \
            md_f32x8_t : md_simd_lerp_f32x8, \
            md_f64x2_t : md_simd_lerp_f64x2, \
            md_f64x4_t : md_simd_lerp_f64x4, \
            const md_f32x4_t : md_simd_lerp_f32x4, \
            const md_f32x8_t : md_simd_lerp_f32x8, \
            const md_f64x2_t : md_simd_lerp_f64x2, \
            const md_f64x4_t : md_simd_lerp_f64x4)(a, b, t)

#define md_simd_cubic_spline(p0, p1, p2, p3, t, s) _Generic((p0),  \
            md_f32x4_t : md_simd_cubic_spline_f32x4, \
            md_f32x8_t : md_simd_cubic_spline_f32x8, \
            md_f64x2_t : md_simd_cubic_spline_f64x2, \
            md_f64x4_t : md_simd_cubic_spline_f64x4, \
            const md_f32x4_t : md_simd_cubic_spline_f32x4, \
            const md_f32x8_t : md_simd_cubic_spline_f32x8, \
            const md_f64x2_t : md_simd_cubic_spline_f64x2, \
            const md_f64x4_t : md_simd_cubic_spline_f64x4)(p0, p1, p2, p3, t, s)

#define md_simd_unpack_xyz(x, y, z, stream, stride) _Generic(x,  \
            md_f32x4_t* : md_simd_unpack_xyz_f32x4, \
            md_f32x8_t* : md_simd_unpack_xyz_f32x8)(x, y, z, stream, stride)

#define md_simd_shift_left(x, i) _Generic((x),      \
            md_i32x4_t : md_simd_shift_left_i32x4,  \
            md_i32x8_t : md_simd_shift_left_i32x8,  \
            md_i64x2_t : md_simd_shift_left_i64x2,  \
            md_i64x4_t : md_simd_shift_left_i64x4,  \
            const md_i32x4_t : md_simd_shift_left_i32x4,  \
            const md_i32x8_t : md_simd_shift_left_i32x8,  \
            const md_i64x2_t : md_simd_shift_left_i64x2,  \
            const md_i64x4_t : md_simd_shift_left_i64x4)(x, i)

#define md_simd_shift_right(x, i) _Generic((x),     \
            md_i32x4_t : md_simd_shift_right_i32x4, \
            md_i32x8_t : md_simd_shift_right_i32x8, \
            md_i64x2_t : md_simd_shift_right_i64x2, \
            md_i64x4_t : md_simd_shift_right_i64x4, \
            const md_i32x4_t : md_simd_shift_right_i32x4, \
            const md_i32x8_t : md_simd_shift_right_i32x8, \
            const md_i64x2_t : md_simd_shift_right_i64x2, \
            const md_i64x4_t : md_simd_shift_right_i64x4)(x, i)

#endif

#if md_simd_width_f32 == 8

// Float
#define md_simd_load_f32    md_simd_load_f32x8
#define md_simd_load_f64    md_simd_load_f64x4

#define md_simd_store_f32   md_simd_store_f32x8
#define md_simd_store_f64   md_simd_store_f64x4

#define md_simd_set1_f32    md_simd_set1_f32x8
#define md_simd_set1_f64    md_simd_set1_f64x4

#define md_simd_set_f32     md_simd_set_f32x8
#define md_simd_set_f64     md_simd_set_f64x4

#define md_simd_zero_f32    md_simd_zero_f32x8
#define md_simd_zero_f64    md_simd_zero_f64x4

#define md_simd_convert_f32 md_simd_convert_f32x8
#define md_simd_convert_f64 md_simd_convert_f64x4

#define md_simd_cast_f32 md_simd_cast_f32x8
#define md_simd_cast_f64 md_simd_cast_f64x4

// Int
#define md_simd_load_i32    md_simd_load_i32x8
#define md_simd_load_i64    md_simd_load_i64x4

#define md_simd_store_i32   md_simd_store_i32x8
#define md_simd_store_i64   md_simd_store_i64x4

#define md_simd_set1_i32    md_simd_set1_i32x8
#define md_simd_set1_i64    md_simd_set1_i64x4

#define md_simd_set_i32     md_simd_set_i32x8
#define md_simd_set_i64     md_simd_set_i64x4

#define md_simd_zero_i32    md_simd_zero_i32x8
#define md_simd_zero_i64    md_simd_zero_i64x4

#define md_simd_convert_i32 md_simd_convert_i32x8
#define md_simd_convert_i64 md_simd_convert_i64x4

#define md_simd_cast_i32    md_simd_cast_i32x8
#define md_simd_cast_i64    md_simd_cast_i64x4

#elif md_simd_width_f32 == 4
// Float
#define md_simd_load_f32    md_simd_load_f32x4
#define md_simd_load_f64    md_simd_load_f64x2

#define md_simd_store_f32   md_simd_store_f32x4
#define md_simd_store_f64   md_simd_store_f64x2

#define md_simd_set1_f32    md_simd_set1_f32x4
#define md_simd_set1_f64    md_simd_set1_f64x2

#define md_simd_set_f32     md_simd_set_f32x4
#define md_simd_set_f64     md_simd_set_f64x2

#define md_simd_zero_f32    md_simd_zero_f32x4
#define md_simd_zero_f64    md_simd_zero_f64x2

#define md_simd_convert_f32 md_simd_convert_f32x4
#define md_simd_convert_f64 md_simd_convert_f64x2

#define md_simd_cast_f32 md_simd_cast_f32x4
#define md_simd_cast_f64 md_simd_cast_f64x2

// Int
#define md_simd_load_i32    md_simd_load_i32x4
#define md_simd_load_i64    md_simd_load_i64x2

#define md_simd_store_i32   md_simd_store_i32x4
#define md_simd_store_i64   md_simd_store_i64x2

#define md_simd_set1_i32    md_simd_set1_i32x4
#define md_simd_set1_i64    md_simd_set1_i64x2

#define md_simd_set_i32     md_simd_set_i32x4
#define md_simd_set_i64     md_simd_set_i64x2

#define md_simd_zero_i32    md_simd_zero_i32x4
#define md_simd_zero_i64    md_simd_zero_i64x2

#define md_simd_convert_i32 md_simd_convert_i32x4
#define md_simd_convert_i64 md_simd_convert_i64x2

#define md_simd_cast_i32    md_simd_cast_i32x4
#define md_simd_cast_i64    md_simd_cast_i64x2

#endif

