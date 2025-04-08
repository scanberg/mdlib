/*

SIMDE handles most of the portability issues for us.
It automatically maps the intrinsics written in x86-64 to the corresponding ARM NEON intrinsics.
It also emulates any lower instruction set, e.g. SSE2 when we are compiling for AVX2, by double-pumping the intrinsics.
We have chosen a sufficiently high baseline, AVX2 so that we can use the same code for both x86-64 and ARM.
In the future, when the support for AVX512 matures, or it is superseeded by something else, we may transition to that.

*/

#pragma once

#include <core/md_common.h>

#if defined(__FMA__)
#define SIMDE_X86_FMA_NATIVE
#endif

#if defined(__AVX512F__)
#include <immintrin.h>
#endif

#if MD_COMPILER_GCC
// Disable warning for _Float16 type on gcc
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

#include <simde/x86/avx2.h>
#include <simde/x86/fma.h>

#if MD_COMPILER_GCC
#pragma GCC diagnostic pop
#endif

#define MD_SIMD_INLINE static FORCE_INLINE

#define md_128  simde__m128
#define md_128d simde__m128d
#define md_128i simde__m128i
#define md_256  simde__m256
#define md_256d simde__m256d
#define md_256i simde__m256i

#if defined(__AVX512F__)
#define md_512  __m512
#define md_512d __m512d
#define md_512i __m512i
// @TODO: Add mask types
#endif

// BASE FLOAT OPERATIONS

#define md_mm_loadu_ps simde_mm_loadu_ps
#define md_mm_loadu_pd simde_mm_loadu_pd
#define md_mm256_loadu_ps simde_mm256_loadu_ps
#define md_mm256_loadu_pd simde_mm256_loadu_pd

#define md_mm_load_ps simde_mm_load_ps
#define md_mm_load_pd simde_mm_load_pd
#define md_mm256_load_ps simde_mm256_load_ps
#define md_mm256_load_pd simde_mm256_load_pd

#define md_mm_i32gather_ps simde_mm_i32gather_ps
#define md_mm_i32gather_pd simde_mm_i32gather_pd
#define md_mm256_i32gather_ps simde_mm256_i32gather_ps
#define md_mm256_i32gather_pd simde_mm256_i32gather_pd

#define md_mm_mask_i32gather_ps simde_mm_mask_i32gather_ps
#define md_mm_mask_i32gather_pd simde_mm_mask_i32gather_pd
#define md_mm256_mask_i32gather_ps simde_mm256_mask_i32gather_ps
#define md_mm256_mask_i32gather_pd simde_mm256_mask_i32gather_pd

#define md_mm_storeu_ps simde_mm_storeu_ps
#define md_mm_storeu_pd simde_mm_storeu_pd
#define md_mm256_storeu_ps simde_mm256_storeu_ps
#define md_mm256_storeu_pd simde_mm256_storeu_pd

#define md_mm_store_ps simde_mm_store_ps
#define md_mm_store_pd simde_mm_store_pd
#define md_mm256_store_ps simde_mm256_store_ps
#define md_mm256_store_pd simde_mm256_store_pd

#define md_mm_set1_ps simde_mm_set1_ps
#define md_mm_set1_pd simde_mm_set1_pd
#define md_mm256_set1_ps simde_mm256_set1_ps
#define md_mm256_set1_pd simde_mm256_set1_pd

#define md_mm_set_ps simde_mm_set_ps
#define md_mm_set_pd simde_mm_set_pd
#define md_mm256_set_ps simde_mm256_set_ps
#define md_mm256_set_pd simde_mm256_set_pd

#define md_mm_setzero_ps simde_mm_setzero_ps
#define md_mm_setzero_pd simde_mm_setzero_pd
#define md_mm256_setzero_ps simde_mm256_setzero_ps
#define md_mm256_setzero_pd simde_mm256_setzero_pd

#define md_mm_and_ps simde_mm_and_ps
#define md_mm_and_pd simde_mm_and_pd
#define md_mm256_and_ps simde_mm256_and_ps
#define md_mm256_and_pd simde_mm256_and_pd

#define md_mm_or_ps simde_mm_or_ps
#define md_mm_or_pd simde_mm_or_pd
#define md_mm256_or_ps simde_mm256_or_ps
#define md_mm256_or_pd simde_mm256_or_pd

#define md_mm_xor_ps simde_mm_xor_ps
#define md_mm_xor_pd simde_mm_xor_pd
#define md_mm256_xor_ps simde_mm256_xor_ps
#define md_mm256_xor_pd simde_mm256_xor_pd

#define md_mm_andnot_ps simde_mm_andnot_ps
#define md_mm_andnot_pd simde_mm_andnot_pd
#define md_mm256_andnot_ps simde_mm256_andnot_ps
#define md_mm256_andnot_pd simde_mm256_andnot_pd

#define md_mm_cmp_ps simde_mm_cmp_ps
#define md_mm_cmp_pd simde_mm_cmp_pd
#define md_mm256_cmp_ps simde_mm_cmp_ps
#define md_mm256_cmp_pd simde_mm_cmp_pd

#define md_mm_cmpeq_ps simde_mm_cmpeq_ps
#define md_mm_cmpeq_pd simde_mm_cmpeq_pd
#define md_mm256_cmpeq_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_EQ_OQ)
#define md_mm256_cmpeq_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_EQ_OQ)

#define md_mm_cmpneq_ps simde_mm_cmpneq_ps
#define md_mm_cmpneq_pd simde_mm_cmpneq_pd
#define md_mm256_cmpneq_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_NEQ_OQ)
#define md_mm256_cmpneq_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_NEQ_OQ)

#define md_mm_cmplt_ps simde_mm_cmplt_ps
#define md_mm_cmplt_pd simde_mm_cmplt_pd
#define md_mm256_cmplt_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_LT_OQ)
#define md_mm256_cmplt_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_LT_OQ)

#define md_mm_cmple_ps simde_mm_cmple_ps
#define md_mm_cmple_pd simde_mm_cmple_pd
#define md_mm256_cmple_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_LE_OQ)
#define md_mm256_cmple_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_LE_OQ)

#define md_mm_cmpgt_ps simde_mm_cmpgt_ps
#define md_mm_cmpgt_pd simde_mm_cmpgt_pd
#define md_mm256_cmpgt_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_GT_OQ)
#define md_mm256_cmpgt_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_GT_OQ)

#define md_mm_cmpge_ps simde_mm_cmpge_ps
#define md_mm_cmpge_pd simde_mm_cmpge_pd
#define md_mm256_cmpge_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_GE_OQ)
#define md_mm256_cmpge_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_GE_OQ)

#define md_mm_round_ps(x) simde_mm_round_ps(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC)
#define md_mm_round_pd(x) simde_mm_round_pd(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC)
#define md_mm256_round_ps(x) simde_mm256_round_ps(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC)
#define md_mm256_round_pd(x) simde_mm256_round_pd(x, SIMDE_MM_FROUND_TO_NEAREST_INT | SIMDE_MM_FROUND_NO_EXC)

#define md_mm_floor_ps simde_mm_floor_ps
#define md_mm_floor_pd simde_mm_floor_pd
#define md_mm256_floor_ps simde_mm256_floor_ps
#define md_mm256_floor_pd simde_mm256_floor_pd

#define md_mm_ceil_ps simde_mm_ceil_ps
#define md_mm_ceil_pd simde_mm_ceil_pd
#define md_mm256_ceil_ps simde_mm256_ceil_ps
#define md_mm256_ceil_pd simde_mm256_ceil_pd

#define md_mm_add_ps simde_mm_add_ps
#define md_mm_add_pd simde_mm_add_pd
#define md_mm256_add_ps simde_mm256_add_ps
#define md_mm256_add_pd simde_mm256_add_pd

#define md_mm_add_ss simde_mm_add_ss
#define md_mm_add_sd simde_mm_add_sd
#define md_mm256_add_ss simde_mm256_add_ss
#define md_mm256_add_sd simde_mm256_add_sd

#define md_mm_sub_ps simde_mm_sub_ps
#define md_mm_sub_pd simde_mm_sub_pd
#define md_mm256_sub_ps simde_mm256_sub_ps
#define md_mm256_sub_pd simde_mm256_sub_pd

#define md_mm_mul_ps simde_mm_mul_ps
#define md_mm_mul_pd simde_mm_mul_pd
#define md_mm256_mul_ps simde_mm256_mul_ps
#define md_mm256_mul_pd simde_mm256_mul_pd

#define md_mm_div_ps simde_mm_div_ps
#define md_mm_div_pd simde_mm_div_pd
#define md_mm256_div_ps simde_mm256_div_ps
#define md_mm256_div_pd simde_mm256_div_pd

#define md_mm_fmadd_ps simde_mm_fmadd_ps
#define md_mm_fmadd_pd simde_mm_fmadd_pd
#define md_mm256_fmadd_ps simde_mm256_fmadd_ps
#define md_mm256_fmadd_pd simde_mm256_fmadd_pd

#define md_mm_fnmadd_ps simde_mm_fnmadd_ps
#define md_mm_fnmadd_pd simde_mm_fnmadd_pd
#define md_mm256_fnmadd_ps simde_mm256_fnmadd_ps
#define md_mm256_fnmadd_pd simde_mm256_fnmadd_pd

#define md_mm_fmsub_ps simde_mm_fmsub_ps
#define md_mm_fmsub_pd simde_mm_fmsub_pd
#define md_mm256_fmsub_ps simde_mm256_fmsub_ps
#define md_mm256_fmsub_pd simde_mm256_fmsub_pd

#define md_mm_fnmsub_ps simde_mm_fnmsub_ps
#define md_mm_fnmsub_pd simde_mm_fnmsub_pd
#define md_mm256_fnmsub_ps simde_mm256_fnmsub_ps
#define md_mm256_fnmsub_pd simde_mm256_fnmsub_pd

#define md_mm_blendv_ps simde_mm_blendv_ps
#define md_mm_blendv_pd simde_mm_blendv_pd
#define md_mm256_blendv_ps simde_mm256_blendv_ps
#define md_mm256_blendv_pd simde_mm256_blendv_pd

#define MD_SIMD_BLEND_MASK(x,y,z,w) (((w) << 3) | ((z) << 2) | ((y) << 1) | (x))

#define md_mm_blend_ps simde_mm_blend_ps
#define md_mm_blend_pd simde_mm_blend_pd
#define md_mm256_blend_ps simde_mm256_blend_ps
#define md_mm256_blend_pd simde_mm256_blend_pd

#define md_mm_sqrt_ps simde_mm_sqrt_ps
#define md_mm_sqrt_pd simde_mm_sqrt_pd
#define md_mm256_sqrt_ps simde_mm256_sqrt_ps
#define md_mm256_sqrt_pd simde_mm256_sqrt_pd

#define md_mm_rsqrt_ps simde_mm_rsqrt_ps
#define md_mm_rsqrt_pd simde_mm_rsqrt_pd
#define md_mm256_rsqrt_ps simde_mm256_rsqrt_ps
#define md_mm256_rsqrt_pd simde_mm256_rsqrt_pd

#define md_mm_min_ps simde_mm_min_ps   
#define md_mm_min_pd simde_mm_min_pd   
#define md_mm256_min_ps simde_mm256_min_ps
#define md_mm256_min_pd simde_mm256_min_pd

#define md_mm_max_ps simde_mm_max_ps   
#define md_mm_max_pd simde_mm_max_pd   
#define md_mm256_max_ps simde_mm256_max_ps
#define md_mm256_max_pd simde_mm256_max_pd

#define md_mm_shuffle_ps simde_mm_shuffle_ps
#define md_mm_shuffle_pd simde_mm_shuffle_pd
#define md_mm256_shuffle_ps simde_mm256_shuffle_ps
#define md_mm256_shuffle_pd simde_mm256_shuffle_pd

#define md_mm_movehl_ps simde_mm_movehl_ps
#define md_mm_movemask_ps simde_mm_movemask_ps
#define md_mm_movemask_pd simde_mm_movemask_pd
#define md_mm256_movemask_ps simde_mm256_movemask_ps
#define md_mm256_movemask_pd simde_mm256_movemask_pd

// BASE INT OPERATIONS

#define md_mm_loadu_si128 simde_mm_loadu_si128
#define md_mm_loadu_epi32 simde_mm_loadu_si32
#define md_mm_loadu_epi64 simde_mm_loadu_si64
#define md_mm256_loadu_si256 simde_mm256_loadu_si256
#define md_mm256_loadu_epi32 simde_mm256_loadu_epi32
#define md_mm256_loadu_epi64 simde_mm256_loadu_epi64

#define md_mm_load_si128 simde_mm_load_si128
#define md_mm_load_epi32 simde_mm_load_si128
#define md_mm_load_epi64 simde_mm_load_si128
#define md_mm256_load_si256 simde_mm256_load_si256
#define md_mm256_load_epi32 simde_mm256_load_si256
#define md_mm256_load_epi64 simde_mm256_load_si256

#define md_mm_i32gather_epi32 simde_mm_i32gather_epi32
#define md_mm_i32gather_epi64 simde_mm_i32gather_epi64
#define md_mm256_i32gather_epi32 simde_mm256_i32gather_epi32
#define md_mm256_i32gather_epi64 simde_mm256_i32gather_epi64

#define md_mm_storeu_si128 simde_mm_storeu_si128
#define md_mm_storeu_epi32 simde_mm_storeu_si128
#define md_mm_storeu_epi64 simde_mm_storeu_si128
#define md_mm256_storeu_si256 simde_mm256_storeu_si256
#define md_mm256_storeu_epi32 simde_mm256_storeu_si256
#define md_mm256_storeu_epi64 simde_mm256_storeu_si256

#define md_mm_store_si128 simde_mm_store_si128
#define md_mm_store_epi32 simde_mm_store_si128
#define md_mm_store_epi64 simde_mm_store_si128
#define md_mm256_store_si256 simde_mm256_store_si256
#define md_mm256_store_epi32 simde_mm256_store_si256
#define md_mm256_store_epi64 simde_mm256_store_si256

#define md_mm_set1_epi32 simde_mm_set1_epi32
#define md_mm_set1_epi64 simde_mm_set1_epi64x
#define md_mm256_set1_epi32 simde_mm256_set1_epi32
#define md_mm256_set1_epi64 simde_mm256_set1_epi64x

#define md_mm_set_epi32 simde_mm_set_epi32
#define md_mm_set_epi64 simde_mm_set_epi64x
#define md_mm256_set_epi32 simde_mm256_set_epi32
#define md_mm256_set_epi64 simde_mm256_set_epi64x

#define md_mm_setzero_si128 simde_mm_setzero_si128
#define md_mm256_setzero_si256 simde_mm256_setzero_si256

#define md_mm_and_si128 simde_mm_and_si128
#define md_mm_and_epi32 simde_mm_and_si128
#define md_mm_and_epi64 simde_mm_and_si128
#define md_mm256_and_si256 simde_mm256_and_si256
#define md_mm256_and_epi32 simde_mm256_and_si256
#define md_mm256_and_epi64 simde_mm256_and_si256

#define md_mm_or_si128 simde_mm_or_si128
#define md_mm256_or_si256 simde_mm256_or_si256

#define md_mm_xor_si128 simde_mm_xor_si128
#define md_mm256_xor_si256 simde_mm256_xor_si256

#define md_mm_andnot_si128 simde_mm_andnot_si128
#define md_mm256_andnot_si256 simde_mm256_andnot_si256

#define md_mm_blendv_epi8 simde_mm_blendv_epi8
#define md_mm256_blendv_epi8 simde_mm256_blendv_epi8

#define md_mm_not_si128(a) md_mm_xor_si128(md_mm_set1_epi32(0xFFFFFFFF), a)
#define md_mm256_not_si256(a) md_mm256_xor_si256(md_mm256_set1_epi32(0xFFFFFFFF), a)

#define md_mm_add_epi32 simde_mm_add_epi32
#define md_mm_add_epi64 simde_mm_add_epi64
#define md_mm256_add_epi32 simde_mm256_add_epi32
#define md_mm256_add_epi64 simde_mm256_add_epi64

#define md_mm_sub_epi32 simde_mm_sub_epi32
#define md_mm_sub_epi64 simde_mm_sub_epi64
#define md_mm256_sub_epi32 simde_mm256_sub_epi32
#define md_mm256_sub_epi64 simde_mm256_sub_epi64

#define md_mm_slli_epi32 simde_mm_slli_epi32
#define md_mm_slli_epi64 simde_mm_slli_epi64
#define md_mm256_slli_epi32 simde_mm256_slli_epi32
#define md_mm256_slli_epi64 simde_mm256_slli_epi64

#define md_mm_srli_epi32 simde_mm_srli_epi32
#define md_mm_srli_epi64 simde_mm_srli_epi64
#define md_mm256_srli_epi32 simde_mm256_srli_epi32
#define md_mm256_srli_epi64 simde_mm256_srli_epi64

#define md_mm_cmpeq_epi32 simde_mm_cmpeq_epi32
#define md_mm_cmpeq_epi64 simde_mm_cmpeq_epi64
#define md_mm256_cmpeq_epi32 simde_mm256_cmpeq_epi32
#define md_mm256_cmpeq_epi64 simde_mm256_cmpeq_epi64

#define md_mm_cmpneq_epi32(a,b) md_mm_not_si128(md_mm_cmpeq_epi32(a,b))
#define md_mm_cmpneq_epi64(a,b) md_mm_not_si128(md_mm_cmpeq_epi64(a,b))
#define md_mm256_cmpneq_epi32(a,b) md_mm256_not_si256(md_mm256_cmpeq_epi32(a,b))
#define md_mm256_cmpneq_epi64(a,b) md_mm256_not_si256(md_mm256_cmpeq_epi64(a,b))

#define md_mm_cmpgt_epi32 simde_mm_cmpgt_epi32
#define md_mm_cmpgt_epi64 simde_mm_cmpgt_epi64
#define md_mm256_cmpgt_epi32 simde_mm256_cmpgt_epi32
#define md_mm256_cmpgt_epi64 simde_mm256_cmpgt_epi64

#define md_mm_cmplt_epi32 simde_mm_cmplt_epi32

#define md_mm_movemask_epi8 simde_mm_movemask_epi8
#define md_mm256_movemask_epi8 simde_mm256_movemask_epi8

// CASTS AND CONVERSIONS

#define md_mm_castps_si128 simde_mm_castps_si128
#define md_mm_castpd_si128 simde_mm_castpd_si128
#define md_mm256_castps_si256 simde_mm256_castps_si256
#define md_mm256_castpd_si256 simde_mm256_castpd_si256

#define md_mm_castsi128_ps simde_mm_castsi128_ps
#define md_mm_castsi128_pd simde_mm_castsi128_pd
#define md_mm256_castsi256_ps simde_mm256_castsi256_ps
#define md_mm256_castsi256_pd simde_mm256_castsi256_pd

#define md_mm_castps_pd simde_mm_castps_pd
#define md_mm_castpd_ps simde_mm_castpd_ps
#define md_mm256_castps_pd simde_mm256_castps_pd
#define md_mm256_castpd_ps simde_mm256_castpd_ps

#define md_mm_castsi128_pd simde_mm_castsi128_pd
#define md_mm_castpd_si128 simde_mm_castpd_si128
#define md_mm256_castsi256_pd simde_mm256_castsi256_pd
#define md_mm256_castpd_si256 simde_mm256_castpd_si256

#define md_mm_cvtepi32_ps simde_mm_cvtepi32_ps
#define md_mm256_cvtepi32_ps simde_mm256_cvtepi32_ps

#define md_mm_cvtps_epi32 simde_mm_cvtps_epi32
#define md_mm_cvttps_epi32 simde_mm_cvttps_epi32
#define md_mm256_cvtps_epi32 simde_mm256_cvtps_epi32
#define md_mm256_cvttps_epi32 simde_mm256_cvttps_epi32
#define md_mm256_cvtepu8_epi32 simde_mm256_cvtepu8_epi32

#define md_mm_cvtss_f32 simde_mm_cvtss_f32
#define md_mm_cvtsd_f64 simde_mm_cvtsd_f64
#define md_mm256_cvtss_f32 simde_mm256_cvtss_f32
#define md_mm256_cvtsd_f64 simde_mm256_cvtsd_f64

// HIGH LEVEL OPERATIONS

#define MD_LOAD_STRIDED_128(ptr, offset, stride) simde_mm_loadu_ps((const float*)((const char*)ptr + offset * stride_in_bytes))

MD_SIMD_INLINE void md_mm_unpack_xyz_ps(md_128* out_x, md_128* out_y, md_128* out_z, const float* in_xyz, size_t stride_in_bytes) {
    md_128  r0, r1, r2, r3;
    md_128  t0, t1, t2, t3;

    r0 = MD_LOAD_STRIDED_128(in_xyz, 0, stride_in_bytes);
    r1 = MD_LOAD_STRIDED_128(in_xyz, 1, stride_in_bytes);
    r2 = MD_LOAD_STRIDED_128(in_xyz, 2, stride_in_bytes);
    r3 = MD_LOAD_STRIDED_128(in_xyz, 3, stride_in_bytes);

    t0 = simde_mm_unpacklo_ps(r0, r1); // xxyy xxyy
    t1 = simde_mm_unpackhi_ps(r0, r1); // zzww zzww
    t2 = simde_mm_unpacklo_ps(r2, r3); // xxyy xxyy
    t3 = simde_mm_unpackhi_ps(r2, r3); // zzww zzww

    *out_x = simde_mm_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(1,0,1,0));  // xxxx xxxx
    *out_y = simde_mm_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(3,2,3,2));  // yyyy yyyy
    *out_z = simde_mm_shuffle_ps(t1, t3, SIMDE_MM_SHUFFLE(1,0,1,0));  // zzzz zzzz
}

MD_SIMD_INLINE void md_mm256_unpack_xyz_ps(md_256* out_x, md_256* out_y, md_256* out_z, const float* in_xyz, size_t stride_in_bytes) {
    md_256 r0, r1, r2, r3;
    md_256 t0, t1, t2, t3;

    // @TODO: Try and implement this using 256-bit loads and then shuffle all the way.
    // It's a tradeoff between the number of loads issued and the port pressure on the generally few ports that are used for shuffle instructions.
    r0 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyz, 0, stride_in_bytes)), MD_LOAD_STRIDED_128(in_xyz, 4, stride_in_bytes), 1);
    r1 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyz, 1, stride_in_bytes)), MD_LOAD_STRIDED_128(in_xyz, 5, stride_in_bytes), 1);
    r2 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyz, 2, stride_in_bytes)), MD_LOAD_STRIDED_128(in_xyz, 6, stride_in_bytes), 1);
    r3 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyz, 3, stride_in_bytes)), MD_LOAD_STRIDED_128(in_xyz, 7, stride_in_bytes), 1);

    t0 = simde_mm256_unpacklo_ps(r0, r1); // xxyy xxyy
    t1 = simde_mm256_unpackhi_ps(r0, r1); // zzww zzww
    t2 = simde_mm256_unpacklo_ps(r2, r3); // xxyy xxyy
    t3 = simde_mm256_unpackhi_ps(r2, r3); // zzww zzww

    *out_x = simde_mm256_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(1,0,1,0));  // xxxx xxxx
    *out_y = simde_mm256_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(3,2,3,2));  // yyyy yyyy
    *out_z = simde_mm256_shuffle_ps(t1, t3, SIMDE_MM_SHUFFLE(1,0,1,0));  // zzzz zzzz
}

#undef MD_LOAD_STRIDED_128

MD_SIMD_INLINE md_128  md_mm_abs_ps(md_128 a)     { return simde_mm_and_ps(a, simde_mm_castsi128_ps( simde_mm_set1_epi32(0x7FFFFFFF))); }
MD_SIMD_INLINE md_128d md_mm_abs_pd(md_128d a)    { return simde_mm_and_pd(a, simde_mm_castsi128_pd( simde_mm_set1_epi64x(0x7FFFFFFFFFFFFFFF))); }
MD_SIMD_INLINE md_256  md_mm256_abs_ps(md_256 a)  { return simde_mm256_and_ps(a, simde_mm256_castsi256_ps(simde_mm256_set1_epi32(0x7FFFFFFF))); }
MD_SIMD_INLINE md_256d md_mm256_abs_pd(md_256d a) { return simde_mm256_and_pd(a, simde_mm256_castsi256_pd(simde_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF))); }

// This naming convention clashes a bit with the intel intrinsics,
// The operation broadcasts a single component into all components of the vector.
// Maybe swizzle is a better name?
#define md_mm_splat_ps(v, i) simde_mm_shuffle_ps(v,v, SIMDE_MM_SHUFFLE(i, i, i, i))
#define md_mm_splat_pd(v, i) simde_mm_shuffle_pd(v,v, SIMDE_MM_SHUFFLE2(i, i))

MD_SIMD_INLINE md_128  md_mm_fract_ps(md_128 a)  { return simde_mm_sub_ps(a, simde_mm_round_ps(a, SIMDE_MM_FROUND_NO_EXC | SIMDE_MM_FROUND_FLOOR)); }
MD_SIMD_INLINE md_128d md_mm_fract_pd(md_128d a) { return simde_mm_sub_pd(a, simde_mm_round_pd(a, SIMDE_MM_FROUND_NO_EXC | SIMDE_MM_FROUND_FLOOR)); }
MD_SIMD_INLINE md_256  md_mm256_fract_ps(md_256 a)  { return simde_mm256_sub_ps(a, simde_mm256_round_ps(a,  SIMDE_MM_FROUND_NO_EXC | SIMDE_MM_FROUND_FLOOR)); }
MD_SIMD_INLINE md_256d md_mm256_fract_pd(md_256d a) { return simde_mm256_sub_pd(a, simde_mm256_round_pd(a,  SIMDE_MM_FROUND_NO_EXC | SIMDE_MM_FROUND_FLOOR)); }

MD_SIMD_INLINE md_128  md_mm_sign_ps(md_128 a)  { return simde_mm_xor_ps(   simde_mm_and_ps(a,  simde_mm_set1_ps(-0.0)),    simde_mm_set1_ps(1.0));    }
MD_SIMD_INLINE md_128d md_mm_sign_pd(md_128d a) { return simde_mm_xor_pd(   simde_mm_and_pd(a,  simde_mm_set1_pd(-0.0)),    simde_mm_set1_pd(1.0));     }
MD_SIMD_INLINE md_256  md_mm256_sign_ps(md_256 a)  { return simde_mm256_xor_ps(simde_mm256_and_ps(a, simde_mm256_set1_ps(-0.0)), simde_mm256_set1_ps(1.0)); }
MD_SIMD_INLINE md_256d md_mm256_sign_pd(md_256d a) { return simde_mm256_xor_pd(simde_mm256_and_pd(a, simde_mm256_set1_pd(-0.0)), simde_mm256_set1_pd(1.0));  }

MD_SIMD_INLINE float md_mm_reduce_min_ps(md_128 x) {
    md_128 a = simde_mm_min_ps(x, simde_mm_shuffle_ps(x, x, SIMDE_MM_SHUFFLE(0, 0, 3, 2)));
    md_128 b = simde_mm_min_ps(a, simde_mm_shuffle_ps(a, a, SIMDE_MM_SHUFFLE(0, 0, 0, 1)));
    return simde_mm_cvtss_f32(simde_mm_min_ps(a,b));
}

MD_SIMD_INLINE float md_mm_reduce_max_ps(md_128 x) {
    md_128 a = simde_mm_max_ps(x, simde_mm_shuffle_ps(x, x, SIMDE_MM_SHUFFLE(0, 0, 3, 2)));
    md_128 b = simde_mm_max_ps(a, simde_mm_shuffle_ps(a, a, SIMDE_MM_SHUFFLE(0, 0, 0, 1)));
    return simde_mm_cvtss_f32(simde_mm_max_ps(a,b));
}

MD_SIMD_INLINE double md_mm256_reduce_min_pd(md_256d x) {
    md_256d a = simde_mm256_min_pd(x, simde_mm256_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 3, 2)));
    md_256d b = simde_mm256_min_pd(a, simde_mm256_shuffle_pd(a, a, SIMDE_MM_SHUFFLE(0, 0, 0, 1)));
    return simde_mm256_cvtsd_f64(simde_mm256_min_pd(a,b));
}

MD_SIMD_INLINE double md_mm256_reduce_max_pd(md_256d x) {
    md_256d a = simde_mm256_max_pd(x, simde_mm256_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 3, 2)));
    md_256d b = simde_mm256_max_pd(a, simde_mm256_shuffle_pd(a, a, SIMDE_MM_SHUFFLE(0, 0, 0, 1)));
    return simde_mm256_cvtsd_f64(simde_mm256_max_pd(a,b));
}

// From here https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx
MD_SIMD_INLINE double md_mm256_reduce_add_pd(md_256d x) {
    md_128d vlow   = simde_mm256_castpd256_pd128(x);   // low  128
    md_128d vhigh  = simde_mm256_extractf128_pd(x, 1); // high 128
    vlow           = simde_mm_add_pd(vlow, vhigh);     // reduce down to 128
    md_128d high64 = simde_mm_unpackhi_pd(vlow, vlow);
    return           simde_mm_cvtsd_f64(simde_mm_add_sd(vlow, high64));  // reduce to scalar
}

// https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
MD_SIMD_INLINE float md_mm_reduce_add_ps(md_128 x) {
    //__m128 shuf = simde_mm_movehdup_ps(x);        // broadcast elements 3,1 to 2,0 (this instruction is SSE3 and we avoid it by using shuffle instead)
    md_128 shuf = md_mm_shuffle_ps(x, x, SIMDE_MM_SHUFFLE(3, 3, 1, 1));
    md_128 sums = md_mm_add_ps(x, shuf);
    shuf        = md_mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = md_mm_add_ps(sums, shuf);
    return        md_mm_cvtss_f32(sums);
}

MD_SIMD_INLINE double md_mm_reduce_add_pd(md_128d x) { return simde_mm_cvtsd_f64(simde_mm_add_pd(x, md_mm_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 0, 1)))); }
MD_SIMD_INLINE double md_mm_reduce_min_pd(md_128d x) { return simde_mm_cvtsd_f64(simde_mm_min_pd(x, md_mm_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 0, 1)))); }
MD_SIMD_INLINE double md_mm_reduce_max_pd(md_128d x) { return simde_mm_cvtsd_f64(simde_mm_max_pd(x, md_mm_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 0, 1)))); }

MD_SIMD_INLINE float md_mm256_reduce_add_ps(md_256 x) { return md_mm_reduce_add_ps(md_mm_add_ps(simde_mm256_castps256_ps128(x), simde_mm256_extractf128_ps(x, 0x1))); }
MD_SIMD_INLINE float md_mm256_reduce_min_ps(md_256 x) { return md_mm_reduce_min_ps(md_mm_min_ps(simde_mm256_castps256_ps128(x), simde_mm256_extractf128_ps(x, 0x1))); }
MD_SIMD_INLINE float md_mm256_reduce_max_ps(md_256 x) { return md_mm_reduce_max_ps(md_mm_max_ps(simde_mm256_castps256_ps128(x), simde_mm256_extractf128_ps(x, 0x1))); }

#if defined(__AVX512F__)
#define md_mm512_reduce_add_ps _mm512_reduce_add_ps
#define md_mm512_reduce_min_ps _mm512_reduce_min_ps
#define md_mm512_reduce_max_ps _mm512_reduce_max_ps
#endif

// Deperiodize is the operation which calculates closest image of x in relation to r,
// where r is the reference point and p is the period.
// It comes in two flavors, one where only the period is given and one where the reciprocal period (rp) is also given.
// The latter is more efficient if the reciprocal period is already computed.
// The operation is used in the periodic boundary condition calculations.

MD_SIMD_INLINE md_128 md_mm_deperiodize_ps(md_128 x, md_128 r, md_128 p) {
    md_128 d  = md_mm_sub_ps(x, r);
    md_128 dx = md_mm_div_ps(d, p);
    dx = md_mm_sub_ps(dx, md_mm_round_ps(dx));
    md_128 x_prim = md_mm_add_ps(r, md_mm_mul_ps(dx, p));
    return md_mm_blendv_ps(x_prim, x, md_mm_cmpeq_ps(p, md_mm_setzero_ps()));
}

MD_SIMD_INLINE md_128 md_mm_deperiodize2_ps(md_128 x, md_128 r, md_128 p, md_128 rp) {
	md_128 d  = md_mm_sub_ps(x, r);
	md_128 dx = md_mm_mul_ps(d, rp);
	dx = md_mm_sub_ps(dx, md_mm_round_ps(dx));
	md_128 x_prim = md_mm_add_ps(r, md_mm_mul_ps(dx, p));
	return md_mm_blendv_ps(x_prim, x, md_mm_cmpeq_ps(p, md_mm_setzero_ps()));
}

MD_SIMD_INLINE md_128d md_mm_deperiodize_pd(md_128d x, md_128d r, md_128d p) {
    md_128d d  = md_mm_sub_pd(x, r);
    md_128d dx = md_mm_div_pd(d, p);
    dx = md_mm_sub_pd(dx, md_mm_round_pd(dx));
    md_128d x_prim = md_mm_add_pd(r, md_mm_mul_pd(dx, p));
    return md_mm_blendv_pd(x_prim, x, md_mm_cmpeq_pd(p, md_mm_setzero_pd()));
}

MD_SIMD_INLINE md_128d md_mm_deperiodize2_pd(md_128d x, md_128d r, md_128d p, md_128d rp) {
	md_128d d  = md_mm_sub_pd(x, r);
	md_128d dx = md_mm_mul_pd(d, rp);
	dx = md_mm_sub_pd(dx, md_mm_round_pd(dx));
	md_128d x_prim = md_mm_add_pd(r, md_mm_mul_pd(dx, p));
	return md_mm_blendv_pd(x_prim, x, md_mm_cmpeq_pd(p, md_mm_setzero_pd()));
}

MD_SIMD_INLINE md_256 md_mm256_deperiodize_ps(md_256 x, md_256 r, md_256 p) {
    md_256 d  = md_mm256_sub_ps(x, r);
    md_256 dx = md_mm256_div_ps(d, p);
    dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
    md_256 x_prim = md_mm256_add_ps(r, md_mm256_mul_ps(dx, p));
    return md_mm256_blendv_ps(x_prim, x, md_mm256_cmpeq_ps(p, md_mm256_setzero_ps()));
}

MD_SIMD_INLINE md_256 md_mm256_deperiodize2_ps(md_256 x, md_256 r, md_256 p, md_256 rp) {
	md_256 d  = md_mm256_sub_ps(x, r);
	md_256 dx = md_mm256_mul_ps(d, rp);
	dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
	md_256 x_prim = md_mm256_add_ps(r, md_mm256_mul_ps(dx, p));
	return md_mm256_blendv_ps(x_prim, x, md_mm256_cmpeq_ps(p, md_mm256_setzero_ps()));
}

MD_SIMD_INLINE md_256d md_mm256_deperiodize_pd(md_256d x, md_256d r, md_256d p) {
    md_256d d  = md_mm256_sub_pd(x, r);
    md_256d dx = md_mm256_div_pd(d, p);
    dx = md_mm256_sub_pd(dx, md_mm256_round_pd(dx));
    md_256d x_prim = md_mm256_add_pd(r, md_mm256_mul_pd(dx, p));
    return md_mm256_blendv_pd(x_prim, x, md_mm256_cmpeq_pd(p, md_mm256_setzero_pd()));
}

MD_SIMD_INLINE md_256d md_mm256_deperiodize2_pd(md_256d x, md_256d r, md_256d p, md_256d rp) {
	md_256d d  = md_mm256_sub_pd(x, r);
	md_256d dx = md_mm256_mul_pd(d, rp);
	dx = md_mm256_sub_pd(dx, md_mm256_round_pd(dx));
	md_256d x_prim = md_mm256_add_pd(r, md_mm256_mul_pd(dx, p));
	return md_mm256_blendv_pd(x_prim, x, md_mm256_cmpeq_pd(p, md_mm256_setzero_pd()));
}

#if defined(__AVX512F__)
MD_SIMD_INLINE md_512 md_mm512_deperiodize_ps(md_512 x, md_512 r, md_512 p) {
    md_512 d  = _mm512_sub_ps(x, r);
    md_512 dx = _mm512_div_ps(d, p);
    dx = _mm512_sub_ps(dx, _mm512_roundscale_ps(dx, _MM_FROUND_TO_NEAREST_INT));
    md_512 x_prim = _mm512_add_ps(r, _mm512_mul_ps(dx, p));
    __mmask16 mask = _mm512_cmpeq_ps_mask(p, _mm512_setzero_ps());
    return _mm512_mask_blend_ps(mask, x_prim, x);
}

MD_SIMD_INLINE md_512 md_mm512_deperiodize2_ps(md_512 x, md_512 r, md_512 p, md_512 rp) {
	md_512 d  = _mm512_sub_ps(x, r);
	md_512 dx = _mm512_mul_ps(d, rp);
	dx = _mm512_sub_ps(dx, _mm512_roundscale_ps(dx, _MM_FROUND_TO_NEAREST_INT));
	md_512 x_prim = _mm512_add_ps(r, _mm512_mul_ps(dx, p));
	__mmask16 mask = _mm512_cmpeq_ps_mask(p, _mm512_setzero_ps());
	return _mm512_mask_blend_ps(mask, x_prim, x);
}

MD_SIMD_INLINE md_512d md_mm512_deperiodize_pd(md_512d x, md_512d r, md_512d p) {
    md_512d d  = _mm512_sub_pd(x, r);
    md_512d dx = _mm512_div_pd(d, p);
    dx = _mm512_sub_pd(dx, _mm512_roundscale_pd(dx, _MM_FROUND_TO_NEAREST_INT ));
    md_512d x_prim = _mm512_add_pd(r, _mm512_mul_pd(dx, p));
    __mmask8 mask = _mm512_cmpeq_pd_mask(p, _mm512_setzero_pd());
    return _mm512_mask_blend_pd(mask , x_prim, x);
}

MD_SIMD_INLINE md_512d md_mm512_deperiodize2_pd(md_512d x, md_512d r, md_512d p, md_512d rp) {
	md_512d d  = _mm512_sub_pd(x, r);
	md_512d dx = _mm512_mul_pd(d, rp);
	dx = _mm512_sub_pd(dx, _mm512_roundscale_pd(dx, _MM_FROUND_TO_NEAREST_INT));
	md_512d x_prim = _mm512_add_pd(r, _mm512_mul_pd(dx, p));
	__mmask8 mask = _mm512_cmpeq_pd_mask(p, _mm512_setzero_pd());
	return _mm512_mask_blend_pd(mask, x_prim, x);
}

#endif

MD_SIMD_INLINE md_128 md_mm_minimage_ps(md_128 dx, md_128 p, md_128 rp) {
    return md_mm_sub_ps(dx, md_mm_mul_ps(p, md_mm_round_ps(md_mm_mul_ps(dx, rp))));
}

MD_SIMD_INLINE md_128d md_mm_minimage_pd(md_128d dx, md_128d p, md_128d rp) {
    return md_mm_sub_pd(dx, md_mm_mul_pd(p, md_mm_round_pd(md_mm_mul_pd(dx, rp))));
}

MD_SIMD_INLINE md_256 md_mm256_minimage_ps(md_256 dx, md_256 p, md_256 rp) {
    return md_mm256_sub_ps(dx, md_mm256_mul_ps(p, md_mm256_round_ps(md_mm256_mul_ps(dx, rp))));
}

MD_SIMD_INLINE md_256d md_mm256_minimage_pd(md_256d dx, md_256d p, md_256d rp) {
    return md_mm256_sub_pd(dx, md_mm256_mul_pd(p, md_mm256_round_pd(md_mm256_mul_pd(dx, rp))));
}

MD_SIMD_INLINE md_128 md_mm_cubic_spline_ps(md_128 p0, md_128 p1, md_128 p2, md_128 p3, md_128 t, md_128 s) {
    const md_128 t1 = t;
    const md_128 t2 = md_mm_mul_ps(t, t);
    const md_128 t3 = md_mm_mul_ps(t2, t);
    const md_128 v0 = md_mm_mul_ps(md_mm_sub_ps(p2, p0), s);
    const md_128 v1 = md_mm_mul_ps(md_mm_sub_ps(p3, p1), s);
    const md_128 x0 = md_mm_add_ps(md_mm_mul_ps(md_mm_set1_ps(2), md_mm_sub_ps(p1, p2)), md_mm_add_ps(v0, v1));
    const md_128 x1 = md_mm_sub_ps(md_mm_mul_ps(md_mm_set1_ps(3), md_mm_sub_ps(p2, p1)), md_mm_add_ps(md_mm_mul_ps(md_mm_set1_ps(2), v0), v1));
    const md_128 r0 = md_mm_add_ps(md_mm_mul_ps(x0, t3), md_mm_mul_ps(x1, t2));
    const md_128 r1 = md_mm_add_ps(md_mm_mul_ps(v0, t1), p1);
    return md_mm_add_ps(r0, r1);
}

MD_SIMD_INLINE md_256 md_mm256_cubic_spline_ps(md_256 p0, md_256 p1, md_256 p2, md_256 p3, md_256 t, md_256 s) {
    const md_256 t2 = md_mm256_mul_ps(t, t);
    const md_256 t3 = md_mm256_mul_ps(t2, t);
    const md_256 v0 = md_mm256_mul_ps(md_mm256_sub_ps(p2, p0), s);
    const md_256 v1 = md_mm256_mul_ps(md_mm256_sub_ps(p3, p1), s);
    const md_256 x0 = md_mm256_add_ps(md_mm256_mul_ps(md_mm256_set1_ps(2), md_mm256_sub_ps(p1, p2)), md_mm256_add_ps(v0, v1));
    const md_256 x1 = md_mm256_sub_ps(md_mm256_mul_ps(md_mm256_set1_ps(3), md_mm256_sub_ps(p2, p1)), md_mm256_add_ps(md_mm256_mul_ps(md_mm256_set1_ps(2), v0), v1));
    const md_256 r0 = md_mm256_add_ps(md_mm256_mul_ps(x0, t3), md_mm256_mul_ps(x1, t2));
    const md_256 r1 = md_mm256_add_ps(md_mm256_mul_ps(v0, t), p1);
    return md_mm256_add_ps(r0, r1);
}

MD_SIMD_INLINE md_128d md_mm_cubic_spline_pd(md_128d p0, md_128d p1, md_128d p2, md_128d p3, md_128d t, md_128d s) {
    const md_128d t2 = md_mm_mul_pd(t, t);
    const md_128d t3 = md_mm_mul_pd(t2, t);
    const md_128d v0 = md_mm_mul_pd(md_mm_sub_pd(p2, p0), s);
    const md_128d v1 = md_mm_mul_pd(md_mm_sub_pd(p3, p1), s);
    const md_128d x0 = md_mm_add_pd(md_mm_mul_pd(md_mm_set1_pd(2), md_mm_sub_pd(p1, p2)), md_mm_add_pd(v0, v1));
    const md_128d x1 = md_mm_sub_pd(md_mm_mul_pd(md_mm_set1_pd(3), md_mm_sub_pd(p2, p1)), md_mm_add_pd(md_mm_mul_pd(md_mm_set1_pd(2), v0), v1));
    const md_128d r0 = md_mm_add_pd(md_mm_mul_pd(x0, t3), md_mm_mul_pd(x1, t2));
    const md_128d r1 = md_mm_add_pd(md_mm_mul_pd(v0, t), p1);
    return md_mm_add_pd(r0, r1);
}

MD_SIMD_INLINE md_256d md_mm256_cubic_spline_pd(md_256d p0, md_256d p1, md_256d p2, md_256d p3, md_256d t, md_256d s) {
    const md_256d t2 = md_mm256_mul_pd(t, t);
    const md_256d t3 = md_mm256_mul_pd(t2, t);
    const md_256d v0 = md_mm256_mul_pd(md_mm256_sub_pd(p2, p0), s);
    const md_256d v1 = md_mm256_mul_pd(md_mm256_sub_pd(p3, p1), s);
    const md_256d x0 = md_mm256_add_pd(md_mm256_mul_pd(md_mm256_set1_pd(2), md_mm256_sub_pd(p1, p2)), md_mm256_add_pd(v0, v1));
    const md_256d x1 = md_mm256_sub_pd(md_mm256_mul_pd(md_mm256_set1_pd(3), md_mm256_sub_pd(p2, p1)), md_mm256_add_pd(md_mm256_mul_pd(md_mm256_set1_pd(2), v0), v1));
    const md_256d r0 = md_mm256_add_pd(md_mm256_mul_pd(x0, t3), md_mm256_mul_pd(x1, t2));
    const md_256d r1 = md_mm256_add_pd(md_mm256_mul_pd(v0, t), p1);
    return md_mm256_add_pd(r0, r1);
}

MD_SIMD_INLINE md_128 md_mm_step_ps(md_128 edge, md_128 x) {
    return md_mm_and_ps(md_mm_cmpge_ps(x, edge), md_mm_set1_ps(1.f));
}

MD_SIMD_INLINE md_128d md_mm_step_pd(md_128d edge, md_128d x) {
    return md_mm_and_pd(md_mm_cmpge_pd(x, edge), md_mm_set1_pd(1.f));
}

MD_SIMD_INLINE md_256 md_mm256_step_ps(md_256 edge, md_256 x) {
    return md_mm256_and_ps(md_mm256_cmpge_ps(x, edge), md_mm256_set1_ps(1.f));
}

MD_SIMD_INLINE md_256d md_mm256_step_pd(md_256d edge, md_256d x) {
    return md_mm256_and_pd(md_mm256_cmpge_pd(x, edge), md_mm256_set1_pd(1.f));
}

MD_SIMD_INLINE md_128 md_mm_lerp_ps(md_128 a, md_128 b, float t) {
    return md_mm_add_ps(md_mm_mul_ps(a, md_mm_set1_ps(1.0f - t)), md_mm_mul_ps(b, md_mm_set1_ps(t)));
}

MD_SIMD_INLINE md_256 md_mm256_lerp_ps(md_256 a, md_256 b, float t) {
    return md_mm256_add_ps(md_mm256_mul_ps(a, md_mm256_set1_ps(1.0f - t)), md_mm256_mul_ps(b, md_mm256_set1_ps(t)));
}

MD_SIMD_INLINE md_128d md_mm_lerp_pd(md_128d a, md_128d b, float t) {
    return md_mm_add_pd(md_mm_mul_pd(a, md_mm_set1_pd(1.0f - t)), md_mm_mul_pd(b, md_mm_set1_pd(t)));
}

MD_SIMD_INLINE md_256d md_mm256_lerp_pd(md_256d a, md_256d b, float t) {
    return md_mm256_add_pd(md_mm256_mul_pd(a, md_mm256_set1_pd(1.0f - t)), md_mm256_mul_pd(b, md_mm256_set1_pd(t)));
}

#if 0
// These are versions from ispc, but they are not as accurate, nor as fast as the cephes versions below
MD_SIMD_INLINE void md_mm_sincos_ps(md_128 in_x, md_128* out_sin, md_128* out_cos) {
    const md_128 two_over_pi = md_mm_set1_ps(0.6366197466850280761718);
    const md_128 scaled = md_mm_mul_ps(in_x, two_over_pi);
    const md_128 k_real = md_mm_floor_ps(scaled);
    const md_128i k = md_mm_cvtps_epi32(k_real);

    const md_128i k_mod4 = md_mm_and_si128(k, md_mm_set1_epi32(3));

    // These can probably be simplified
    const md_128 cos_usecos = simde_mm_castsi128_ps(
        md_mm_or_si128(
            md_mm_cmpeq_epi32(k_mod4, md_mm_setzero_si128()),
            md_mm_cmpeq_epi32(k_mod4, md_mm_set1_epi32(2))
        ));

    const md_128 sin_usecos = md_mm_castsi128_ps(
        md_mm_or_si128(
            md_mm_cmpeq_epi32(k_mod4, md_mm_set1_epi32(1)),
            md_mm_cmpeq_epi32(k_mod4, md_mm_set1_epi32(3))
        ));

    const md_128 sin_sign_bit = md_mm_castsi128_ps(
        md_mm_slli_epi32(
            md_mm_cmpgt_epi32(k_mod4, md_mm_set1_epi32(1)),
            31)
        );

    const md_128 cos_sign_bit = simde_mm_castsi128_ps(
        md_mm_slli_epi32(
            md_mm_or_si128(
                md_mm_cmpeq_epi32(k_mod4, simde_mm_set1_epi32(1)),
                md_mm_cmpeq_epi32(k_mod4, simde_mm_set1_epi32(2))
            ),
            31)
        );

    const md_128 one     = md_mm_set1_ps(1.0f);
    const md_128 sin_c2  = md_mm_set1_ps(-0.16666667163372039794921875);
    const md_128 sin_c4  = md_mm_set1_ps(8.333347737789154052734375e-3);
    const md_128 sin_c6  = md_mm_set1_ps(-1.9842604524455964565277099609375e-4);
    const md_128 sin_c8  = md_mm_set1_ps(2.760012648650445044040679931640625e-6);
    const md_128 sin_c10 = md_mm_set1_ps(-2.50293279435709337121807038784027099609375e-8);

    const md_128 cos_c2  = md_mm_set1_ps(-0.5);
    const md_128 cos_c4  = md_mm_set1_ps(4.166664183139801025390625e-2);
    const md_128 cos_c6  = md_mm_set1_ps(-1.388833043165504932403564453125e-3);
    const md_128 cos_c8  = md_mm_set1_ps(2.47562347794882953166961669921875e-5);
    const md_128 cos_c10 = md_mm_set1_ps(-2.59630184018533327616751194000244140625e-7);

    md_128 x = in_x;

    // FMA-enhanced Cody-Waite style reduction
    x = md_mm_fmadd_ps(k_real, md_mm_set1_ps(-0x1.921fb0p+00f), x);
    x = md_mm_fmadd_ps(k_real, md_mm_set1_ps(-0x1.5110b4p-22f), x);
    x = md_mm_fmadd_ps(k_real, md_mm_set1_ps(-0x1.846988p-48f), x);

    const md_128 x2 = md_mm_mul_ps(x, x);

    md_128 sin_res, cos_res;
    sin_res = md_mm_fmadd_ps(x2, sin_c10, sin_c8);
    cos_res = md_mm_fmadd_ps(x2, cos_c10, cos_c8);

    sin_res = md_mm_fmadd_ps(x2, sin_res, sin_c6);
    cos_res = md_mm_fmadd_ps(x2, cos_res, cos_c6);

    sin_res = md_mm_fmadd_ps(x2, sin_res, sin_c4);
    cos_res = md_mm_fmadd_ps(x2, cos_res, cos_c4);

    sin_res = md_mm_fmadd_ps(x2, sin_res, sin_c2);
    cos_res = md_mm_fmadd_ps(x2, cos_res, cos_c2);

    sin_res = md_mm_fmadd_ps(x2, sin_res, one);
    cos_res = md_mm_fmadd_ps(x2, cos_res, one);

    sin_res = md_mm_mul_ps(x, sin_res);

    *out_sin = md_mm_xor_ps(_mm_blendv_ps(sin_res, cos_res, sin_usecos), sin_sign_bit);
    *out_cos = md_mm_xor_ps(_mm_blendv_ps(sin_res, cos_res, cos_usecos), cos_sign_bit);
}

// We cannot reserve the name mm256_sincos_ps because it is already defined in immintrin.h
MD_SIMD_INLINE void md_mm256_sincos_ps(md_256 in_x, md_256* out_sin, md_256* out_cos) {
    const md_256 two_over_pi = md_mm256_set1_ps(0.6366197466850280761718);
    const md_256 scaled = md_mm256_mul_ps(in_x, two_over_pi);
    const md_256 k_real = md_mm256_floor_ps(scaled);
    const md_256i k = md_mm256_cvtps_epi32(k_real);

    const md_256i k_mod4 = md_mm256_and_si256(k, md_mm256_set1_epi32(3));

    // These can probably be simplified
    const md_256 cos_usecos = md_mm256_castsi256_ps(
        md_mm256_or_si256(
            md_mm256_cmpeq_epi32(k_mod4, md_mm256_setzero_si256()),
            md_mm256_cmpeq_epi32(k_mod4, md_mm256_set1_epi32(2))
        ));

    const md_256 sin_usecos = md_mm256_castsi256_ps(
        md_mm256_or_si256(
            md_mm256_cmpeq_epi32(k_mod4, md_mm256_set1_epi32(1)),
            md_mm256_cmpeq_epi32(k_mod4, md_mm256_set1_epi32(3))
        ));

    const md_256 sin_sign_bit = md_mm256_castsi256_ps(
        md_mm256_slli_epi32(
            md_mm256_cmpgt_epi32(k_mod4, md_mm256_set1_epi32(1)),
            31)
    );

    const md_256 cos_sign_bit = md_mm256_castsi256_ps(
        md_mm256_slli_epi32(
            md_mm256_or_si256(
                md_mm256_cmpeq_epi32(k_mod4, md_mm256_set1_epi32(1)),
                md_mm256_cmpeq_epi32(k_mod4, md_mm256_set1_epi32(2))
            ),
            31)
    );

    const md_256 one     = md_mm256_set1_ps(1.0f);
    const md_256 sin_c2  = md_mm256_set1_ps(-0.16666667163372039794921875);
    const md_256 sin_c4  = md_mm256_set1_ps(8.333347737789154052734375e-3);
    const md_256 sin_c6  = md_mm256_set1_ps(-1.9842604524455964565277099609375e-4);
    const md_256 sin_c8  = md_mm256_set1_ps(2.760012648650445044040679931640625e-6);
    const md_256 sin_c10 = md_mm256_set1_ps(-2.50293279435709337121807038784027099609375e-8);

    const md_256 cos_c2  = md_mm256_set1_ps(-0.5);
    const md_256 cos_c4  = md_mm256_set1_ps(4.166664183139801025390625e-2);
    const md_256 cos_c6  = md_mm256_set1_ps(-1.388833043165504932403564453125e-3);
    const md_256 cos_c8  = md_mm256_set1_ps(2.47562347794882953166961669921875e-5);
    const md_256 cos_c10 = md_mm256_set1_ps(-2.59630184018533327616751194000244140625e-7);

    md_256 x = in_x;

    // FMA-enhanced Cody-Waite style reduction
    x = md_mm256_fmadd_ps(k_real, md_mm256_set1_ps(-0x1.921fb0p+00f), x);
    x = md_mm256_fmadd_ps(k_real, md_mm256_set1_ps(-0x1.5110b4p-22f), x);
    x = md_mm256_fmadd_ps(k_real, md_mm256_set1_ps(-0x1.846988p-48f), x);

    const md_256 x2 = md_mm256_mul_ps(x, x);

    md_256 sin_res, cos_res;
    sin_res = md_mm256_fmadd_ps(x2, sin_c10, sin_c8);
    cos_res = md_mm256_fmadd_ps(x2, cos_c10, cos_c8);

    sin_res = md_mm256_fmadd_ps(x2, sin_res, sin_c6);
    cos_res = md_mm256_fmadd_ps(x2, cos_res, cos_c6);

    sin_res = md_mm256_fmadd_ps(x2, sin_res, sin_c4);
    cos_res = md_mm256_fmadd_ps(x2, cos_res, cos_c4);

    sin_res = md_mm256_fmadd_ps(x2, sin_res, sin_c2);
    cos_res = md_mm256_fmadd_ps(x2, cos_res, cos_c2);

    sin_res = md_mm256_fmadd_ps(x2, sin_res, one);
    cos_res = md_mm256_fmadd_ps(x2, cos_res, one);

    sin_res = md_mm256_mul_ps(x, sin_res);

    *out_sin = md_mm256_xor_ps(_mm256_blendv_ps(sin_res, cos_res, sin_usecos), sin_sign_bit);
    *out_cos = md_mm256_xor_ps(_mm256_blendv_ps(sin_res, cos_res, cos_usecos), cos_sign_bit);
}
#endif

#if defined(_MSC_VER)
#pragma float_control(precise, on, push)
#endif

// These are based on the SSEmathfun implementations by Julien Pommier.
// http://gruntthepeon.free.fr/ssemath/
// Who's implementations which are based on the cephes library
// http://www.netlib.org/cephes/

MD_SIMD_INLINE void md_mm_sincos_ps(md_128 xx, md_128* s, md_128* c) {
    md_128  xmm1, xmm2, sign_bit_sin, x, y, y2;
    md_128i imm0, imm2, imm3;
    
    // sign bit
    sign_bit_sin = md_mm_and_ps(xx, md_mm_castsi128_ps( md_mm_set1_epi32(0x80000000)));

    // abs
    x = md_mm_abs_ps(xx);

    /* scale by 4/Pi */
    y = md_mm_mul_ps(x, md_mm_set1_ps(1.27323954473516f));

    /* store integer part of y */
    imm2 = md_mm_cvttps_epi32(y);

    /* j=(j+1) & (~1) (see the cephes sources) */
    imm2 = md_mm_add_epi32(imm2, md_mm_set1_epi32(1));
    imm2 = md_mm_and_epi32(imm2, md_mm_set1_epi32(~1));

    y = md_mm_cvtepi32_ps(imm2);
    imm3 = imm2;

    /* get the swap sign flag for the sine */
    imm0 = md_mm_and_epi32(imm2, md_mm_set1_epi32(4));
    imm0 = md_mm_slli_epi32(imm0, 29);

    /* get polynom selection mask for the sine */
    imm2 = md_mm_and_epi32(imm2, md_mm_set1_epi32(2));
    imm2 = md_mm_cmpeq_epi32(imm2, md_mm_setzero_si128());

    const md_128 swap_sign_bit_sin = md_mm_castsi128_ps(imm0);
    const md_128 poly_mask         = md_mm_castsi128_ps(imm2);

    imm3 = md_mm_sub_epi32(imm3, md_mm_set1_epi32(2));
    imm3 = md_mm_andnot_si128(imm3, md_mm_set1_epi32(4));
    imm3 = md_mm_slli_epi32(imm3, 29);

    const md_128 sign_bit_cos = md_mm_castsi128_ps(imm3);
    sign_bit_sin = md_mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);

    /* The magic pass: "Extended precision modular arithmetic" 
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    x = md_mm_fmadd_ps(y, md_mm_set1_ps(-0.78515625f), x);
    x = md_mm_fmadd_ps(y, md_mm_set1_ps(-2.4187564849853515625e-4f), x);
    x = md_mm_fmadd_ps(y, md_mm_set1_ps(-3.77489497744594108e-8f), x);

    const md_128 x2 = md_mm_mul_ps(x,x);
    const md_128 x3 = md_mm_mul_ps(x2,x);
    const md_128 x4 = md_mm_mul_ps(x2,x2);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    const md_128 c0 = md_mm_set1_ps(2.443315711809948E-005f);
    const md_128 c1 = md_mm_set1_ps(-1.388731625493765E-003f);
    const md_128 c2 = md_mm_set1_ps(4.166664568298827E-002f);

    y = md_mm_fmadd_ps(x2, md_mm_fmadd_ps(x2, c0, c1), c2);
    //y = md_mm_mul_ps(y, x4);
    //y = md_mm_sub_ps(y, md_mm_mul_ps(x2, md_mm_set1_ps(0.5)));
    y = md_mm_fmadd_ps(x2, md_mm_set1_ps(-0.5), md_mm_mul_ps(y, x4));
    y = md_mm_add_ps(y, md_mm_set1_ps(1.0));

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
    const md_128 s0 = md_mm_set1_ps(-1.9515295891E-4f);
    const md_128 s1 = md_mm_set1_ps(8.3321608736E-3f);
    const md_128 s2 = md_mm_set1_ps(-1.6666654611E-1f);
    y2 = md_mm_fmadd_ps(x2, md_mm_fmadd_ps(x2, s0, s1), s2);
    y2 = md_mm_fmadd_ps(y2, x3, x);

    /* select the correct result from the two polynoms */  
    md_128 ysin2 = md_mm_and_ps(poly_mask, y2);
    md_128 ysin1 = md_mm_andnot_ps(poly_mask, y);

    y2 = md_mm_sub_ps(y2,ysin2);
    y  = md_mm_sub_ps(y, ysin1);

    xmm1 = md_mm_add_ps(ysin1,ysin2);
    xmm2 = md_mm_add_ps(y,y2);

    /* update the sign */
    *s = md_mm_xor_ps(xmm1, sign_bit_sin);
    *c = md_mm_xor_ps(xmm2, sign_bit_cos);
}

MD_SIMD_INLINE void md_mm256_sincos_ps(md_256 x, md_256* s, md_256* c) {
    md_256  xmm1, xmm2, sign_bit_sin, y, y2;
    md_256i imm0, imm2, imm3;

    // sign bit
    sign_bit_sin = md_mm256_and_ps(x, md_mm256_castsi256_ps( md_mm256_set1_epi32(0x80000000)));

    // abs
    x = md_mm256_abs_ps(x);

    /* scale by 4/Pi */
    y = md_mm256_mul_ps(x, md_mm256_set1_ps(1.27323954473516f));

    /* store integer part of y */
    imm2 = md_mm256_cvttps_epi32(y);

    /* j=(j+1) & (~1) (see the cephes sources) */
    imm2 = md_mm256_add_epi32(imm2, md_mm256_set1_epi32(1));
    imm2 = md_mm256_and_epi32(imm2, md_mm256_set1_epi32(~1));

    y = md_mm256_cvtepi32_ps(imm2);
    imm3 = imm2;

    /* get the swap sign flag for the sine */
    imm0 = md_mm256_and_epi32(imm2, md_mm256_set1_epi32(4));
    imm0 = md_mm256_slli_epi32(imm0, 29);

    /* get polynom selection mask for the sine */
    imm2 = md_mm256_and_epi32(imm2, md_mm256_set1_epi32(2));
    imm2 = md_mm256_cmpeq_epi32(imm2, md_mm256_setzero_si256());

    const md_256 swap_sign_bit_sin = md_mm256_castsi256_ps(imm0);
    const md_256 poly_mask         = md_mm256_castsi256_ps(imm2);

    imm3 = md_mm256_sub_epi32(imm3, md_mm256_set1_epi32(2));
    imm3 = md_mm256_andnot_si256(imm3, md_mm256_set1_epi32(4));
    imm3 = md_mm256_slli_epi32(imm3, 29);

    const md_256 sign_bit_cos = md_mm256_castsi256_ps(imm3);
    sign_bit_sin = md_mm256_xor_ps(sign_bit_sin, swap_sign_bit_sin);

    /* The magic pass: "Extended precision modular arithmetic" 
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    x = md_mm256_fmadd_ps(y, md_mm256_set1_ps(-0.78515625f), x);
    x = md_mm256_fmadd_ps(y, md_mm256_set1_ps(-2.4187564849853515625E-4f), x);
    x = md_mm256_fmadd_ps(y, md_mm256_set1_ps(-3.77489470793079817668E-8f), x);

    const md_256 x2 = md_mm256_mul_ps(x,x);
    const md_256 x3 = md_mm256_mul_ps(x2,x);
    const md_256 x4 = md_mm256_mul_ps(x2,x2);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    const md_256 c0 = md_mm256_set1_ps(2.443315711809948E-5f);
    const md_256 c1 = md_mm256_set1_ps(-1.388731625493765E-3f);
    const md_256 c2 = md_mm256_set1_ps(4.166664568298827E-2f);

    y = md_mm256_fmadd_ps(x2, md_mm256_fmadd_ps(x2, c0, c1), c2);
    y = md_mm256_fmadd_ps(x2, md_mm256_set1_ps(-0.5), md_mm256_mul_ps(y, x4));
    y = md_mm256_add_ps(y, md_mm256_set1_ps(1.0));

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
    const md_256 s0 = md_mm256_set1_ps(-1.9515295891E-4f);
    const md_256 s1 = md_mm256_set1_ps(8.3321608736E-3f);
    const md_256 s2 = md_mm256_set1_ps(-1.6666654611E-1f);
    y2 = md_mm256_fmadd_ps(x2, md_mm256_fmadd_ps(x2, s0, s1), s2);
    y2 = md_mm256_fmadd_ps(y2, x3, x);

    /* select the correct result from the two polynoms */  
    md_256 ysin2 = md_mm256_and_ps(poly_mask, y2);
    md_256 ysin1 = md_mm256_andnot_ps(poly_mask, y);

    y2 = md_mm256_sub_ps(y2,ysin2);
    y  = md_mm256_sub_ps(y, ysin1);

    xmm1 = md_mm256_add_ps(ysin1,ysin2);
    xmm2 = md_mm256_add_ps(y,y2);

    /* update the sign */
    *s = md_mm256_xor_ps(xmm1, sign_bit_sin);
    *c = md_mm256_xor_ps(xmm2, sign_bit_cos);
}

#if defined(__AVX512F__) && defined(__AVX512DQ__)
MD_SIMD_INLINE void md_mm512_sincos_ps(__m512 x, __m512* s, __m512* c) {
    __m512  xmm1, xmm2, sign_bit_sin, y, y2;
    __m512i imm0, imm2, imm3;

    // sign bit
    sign_bit_sin = _mm512_and_ps(x, _mm512_castsi512_ps( _mm512_set1_epi32(0x80000000)));

    // abs
    x = _mm512_abs_ps(x);

    /* scale by 4/Pi */
    y = _mm512_mul_ps(x, _mm512_set1_ps(1.27323954473516f));

    /* store integer part of y */
    imm2 = _mm512_cvttps_epi32(y);

    /* j=(j+1) & (~1) (see the cephes sources) */
    imm2 = _mm512_add_epi32(imm2, _mm512_set1_epi32(1));
    imm2 = _mm512_and_epi32(imm2, _mm512_set1_epi32(~1));

    y = _mm512_cvtepi32_ps(imm2);
    imm3 = imm2;

    /* get the swap sign flag for the sine */
    imm0 = _mm512_and_epi32(imm2, _mm512_set1_epi32(4));
    imm0 = _mm512_slli_epi32(imm0, 29);

    /* get polynom selection mask for the sine */
    imm2 = _mm512_and_epi32(imm2, _mm512_set1_epi32(2));
    
    const __mmask16 mask = _mm512_cmpeq_epi32_mask(imm2, _mm512_setzero_si512());
    imm2 = _mm512_maskz_mov_epi32(mask, _mm512_set1_epi32(0xffffffff));

    const __m512 swap_sign_bit_sin = _mm512_castsi512_ps(imm0);
    const __m512 poly_mask         = _mm512_castsi512_ps(imm2);

    imm3 = _mm512_sub_epi32(imm3, _mm512_set1_epi32(2));
    imm3 = _mm512_andnot_epi32(imm3, _mm512_set1_epi32(4));
    imm3 = _mm512_slli_epi32(imm3, 29);

    const __m512 sign_bit_cos = _mm512_castsi512_ps(imm3);
    sign_bit_sin = _mm512_xor_ps(sign_bit_sin, swap_sign_bit_sin);

    /* The magic pass: "Extended precision modular arithmetic" 
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    x = _mm512_fmadd_ps(y, _mm512_set1_ps(-0.78515625f), x);
    x = _mm512_fmadd_ps(y, _mm512_set1_ps(-2.4187564849853515625E-4f), x);
    x = _mm512_fmadd_ps(y, _mm512_set1_ps(-3.77489470793079817668E-8f), x);

    const __m512 x2 = _mm512_mul_ps(x,x);
    const __m512 x3 = _mm512_mul_ps(x2,x);
    const __m512 x4 = _mm512_mul_ps(x2,x2);

    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    const __m512 c0 = _mm512_set1_ps(2.443315711809948E-5f);
    const __m512 c1 = _mm512_set1_ps(-1.388731625493765E-3f);
    const __m512 c2 = _mm512_set1_ps(4.166664568298827E-2f);

    y = _mm512_fmadd_ps(x2, _mm512_fmadd_ps(x2, c0, c1), c2);
    y = _mm512_fmadd_ps(x2, _mm512_set1_ps(-0.5), _mm512_mul_ps(y, x4));
    y = _mm512_add_ps(y, _mm512_set1_ps(1.0));

    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
    const __m512 s0 = _mm512_set1_ps(-1.9515295891E-4f);
    const __m512 s1 = _mm512_set1_ps(8.3321608736E-3f);
    const __m512 s2 = _mm512_set1_ps(-1.6666654611E-1f);
    y2 = _mm512_fmadd_ps(x2, _mm512_fmadd_ps(x2, s0, s1), s2);
    y2 = _mm512_fmadd_ps(y2, x3, x);

    /* select the correct result from the two polynoms */  
    __m512 ysin2 = _mm512_and_ps(poly_mask, y2);
    __m512 ysin1 = _mm512_andnot_ps(poly_mask, y);

    y2 = _mm512_sub_ps(y2,ysin2);
    y  = _mm512_sub_ps(y, ysin1);

    xmm1 = _mm512_add_ps(ysin1,ysin2);
    xmm2 = _mm512_add_ps(y,y2);

    /* update the sign */
    *s = _mm512_xor_ps(xmm1, sign_bit_sin);
    *c = _mm512_xor_ps(xmm2, sign_bit_cos);
}
#endif

MD_SIMD_INLINE md_128 md_mm_exp_ps(md_128 x) {
    md_128 tmp = md_mm_setzero_ps(), fx;
    md_128i emm0;
    md_128 one = md_mm_set1_ps(1.0f);

    x = md_mm_min_ps(x, md_mm_set1_ps( 88.3762626647949f));
    x = md_mm_max_ps(x, md_mm_set1_ps(-88.3762626647949f));

    /* express exp(x) as exp(g + n*log(2)) */
    fx = md_mm_mul_ps(x, md_mm_set1_ps(1.44269504088896341f));
    fx = md_mm_add_ps(fx, md_mm_set1_ps(0.5f));

    /* floor */
    emm0 = md_mm_cvttps_epi32(fx);
    tmp  = md_mm_cvtepi32_ps(emm0);

    /* if greater, substract 1 */
    md_128 mask = md_mm_cmpgt_ps(tmp, fx);    
    mask = md_mm_and_ps(mask, one);
    fx = md_mm_sub_ps(tmp, mask);

    tmp = md_mm_mul_ps(fx, md_mm_set1_ps(0.693359375f));
    md_128 z = md_mm_mul_ps(fx, md_mm_set1_ps(-2.12194440e-4f));
    x = md_mm_sub_ps(x, tmp);
    x = md_mm_sub_ps(x, z);

    z = md_mm_mul_ps(x,x);

    md_128 y = md_mm_set1_ps(1.9875691500E-4f);
    y = md_mm_fmadd_ps(y, x, md_mm_set1_ps(1.3981999507E-3f));
    y = md_mm_fmadd_ps(y, x, md_mm_set1_ps(8.3334519073E-3f));
    y = md_mm_fmadd_ps(y, x, md_mm_set1_ps(4.1665795894E-2f));
    y = md_mm_fmadd_ps(y, x, md_mm_set1_ps(1.6666665459E-1f));
    y = md_mm_fmadd_ps(y, x, md_mm_set1_ps(5.0000001201E-1f));
    y = md_mm_fmadd_ps(y, z, md_mm_add_ps(x, one));

    /* build 2^n */
    emm0 = md_mm_cvttps_epi32(fx);
    emm0 = md_mm_add_epi32(emm0, md_mm_set1_epi32(0x7f));
    emm0 = md_mm_slli_epi32(emm0, 23);
    md_128 pow2n = md_mm_castsi128_ps(emm0);

    y = md_mm_mul_ps(y, pow2n);
    return y;
}

MD_SIMD_INLINE md_256 md_mm256_exp_ps(md_256 x) {
    md_256 tmp = md_mm256_setzero_ps(), fx;
    md_256i emm0;
    md_256 one = md_mm256_set1_ps(1.0f);

    x = md_mm256_min_ps(x, md_mm256_set1_ps( 88.3762626647949f));
    x = md_mm256_max_ps(x, md_mm256_set1_ps(-88.3762626647949f));

    /* express exp(x) as exp(g + n*log(2)) */
    fx = md_mm256_mul_ps(x, md_mm256_set1_ps(1.44269504088896341f));
    fx = md_mm256_add_ps(fx, md_mm256_set1_ps(0.5f));

    /* floor */
    emm0 = md_mm256_cvttps_epi32(fx);
    tmp  = md_mm256_cvtepi32_ps(emm0);

    /* if greater, substract 1 */
    md_256 mask = md_mm256_cmpgt_ps(tmp, fx);    
    mask = md_mm256_and_ps(mask, one);
    fx = md_mm256_sub_ps(tmp, mask);

    tmp = md_mm256_mul_ps(fx, md_mm256_set1_ps(0.693359375f));
    md_256 z = md_mm256_mul_ps(fx, md_mm256_set1_ps(-2.12194440e-4f));
    x = md_mm256_sub_ps(x, tmp);
    x = md_mm256_sub_ps(x, z);

    z = md_mm256_mul_ps(x,x);

    md_256 y = md_mm256_set1_ps(1.9875691500E-4f);
    y = md_mm256_fmadd_ps(y, x, md_mm256_set1_ps(1.3981999507E-3f));
    y = md_mm256_fmadd_ps(y, x, md_mm256_set1_ps(8.3334519073E-3f));
    y = md_mm256_fmadd_ps(y, x, md_mm256_set1_ps(4.1665795894E-2f));
    y = md_mm256_fmadd_ps(y, x, md_mm256_set1_ps(1.6666665459E-1f));
    y = md_mm256_fmadd_ps(y, x, md_mm256_set1_ps(5.0000001201E-1f));
    y = md_mm256_fmadd_ps(y, z, md_mm256_add_ps(x, one));

    /* build 2^n */
    emm0 = md_mm256_cvttps_epi32(fx);
    emm0 = md_mm256_add_epi32(emm0, md_mm256_set1_epi32(0x7f));
    emm0 = md_mm256_slli_epi32(emm0, 23);
    md_256 pow2n = md_mm256_castsi256_ps(emm0);

    y = md_mm256_mul_ps(y, pow2n);
    return y;
}

#ifdef __AVX512F__
MD_SIMD_INLINE __m512 md_mm512_exp_ps(__m512 x) {
    __m512 tmp = _mm512_setzero_ps(), fx;
    __m512i emm0;
    __m512 one = _mm512_set1_ps(1.0f);

    x = _mm512_min_ps(x, _mm512_set1_ps( 88.3762626647949f));
    x = _mm512_max_ps(x, _mm512_set1_ps(-88.3762626647949f));

    /* express exp(x) as exp(g + n*log(2)) */
    fx = _mm512_mul_ps(x, _mm512_set1_ps(1.44269504088896341f));
    fx = _mm512_add_ps(fx, _mm512_set1_ps(0.5f));

    /* floor */
    emm0 = _mm512_cvttps_epi32(fx);
    tmp  = _mm512_cvtepi32_ps(emm0);

    /* if greater, substract 1 */
    __mmask16 mask = _mm512_cmp_ps_mask(tmp, fx, _CMP_GT_OQ);
    fx = _mm512_sub_ps(tmp, _mm512_maskz_mov_ps(mask, _mm512_set1_ps(1.0f)));

    tmp = _mm512_mul_ps(fx, _mm512_set1_ps(0.693359375f));
    __m512 z = _mm512_mul_ps(fx, _mm512_set1_ps(-2.12194440e-4f));
    x = _mm512_sub_ps(x, tmp);
    x = _mm512_sub_ps(x, z);

    z = _mm512_mul_ps(x,x);

    __m512 y = _mm512_set1_ps(1.9875691500E-4f);
    y = _mm512_fmadd_ps(y, x, _mm512_set1_ps(1.3981999507E-3f));
    y = _mm512_fmadd_ps(y, x, _mm512_set1_ps(8.3334519073E-3f));
    y = _mm512_fmadd_ps(y, x, _mm512_set1_ps(4.1665795894E-2f));
    y = _mm512_fmadd_ps(y, x, _mm512_set1_ps(1.6666665459E-1f));
    y = _mm512_fmadd_ps(y, x, _mm512_set1_ps(5.0000001201E-1f));
    y = _mm512_fmadd_ps(y, z, _mm512_add_ps(x, one));

    /* build 2^n */
    emm0 = _mm512_cvttps_epi32(fx);
    emm0 = _mm512_add_epi32(emm0, _mm512_set1_epi32(0x7f));
    emm0 = _mm512_slli_epi32(emm0, 23);
    __m512 pow2n = _mm512_castsi512_ps(emm0);

    y = _mm512_mul_ps(y, pow2n);
    return y;
}
#endif

#if defined(_MSC_VER)
#pragma float_control(pop)
#endif

#if 0
#ifdef __cplusplus
// C++ function overload
MD_SIMD_INLINE md_128 md_mm_add(md_128 a, md_128 b) { return md_simd_add_ps(a, b); }
MD_SIMD_INLINE md_256 md_mm256_add(md_256 a, md_256 b) { return md_simd_add_ps(a, b); }
MD_SIMD_INLINE md_128d md_simd_add(md_128d a, md_128d b) { return md_simd_add_pd(a, b); }
MD_SIMD_INLINE md_256d md_simd_add(md_256d a, md_256d b) { return md_simd_add_pd(a, b); }

MD_SIMD_INLINE md_128 md_mm_sub(md_128 a, md_128 b) { return md_simd_sub_ps(a, b); }
MD_SIMD_INLINE md_256 md_mm256_sub(md_256 a, md_256 b) { return md_simd_sub_ps(a, b); }
MD_SIMD_INLINE md_128d md_simd_sub(md_128d a, md_128d b) { return md_simd_sub_pd(a, b); }
MD_SIMD_INLINE md_256d md_simd_sub(md_256d a, md_256d b) { return md_simd_sub_pd(a, b); }

MD_SIMD_INLINE md_128 md_mm_mul(md_128 a, md_128 b) { return md_simd_mul_ps(a, b); }
MD_SIMD_INLINE md_256 md_mm256_mul(md_256 a, md_256 b) { return md_simd_mul_ps(a, b); }
MD_SIMD_INLINE md_128d md_simd_mul(md_128d a, md_128d b) { return md_simd_mul_pd(a, b); }
MD_SIMD_INLINE md_256d md_simd_mul(md_256d a, md_256d b) { return md_simd_mul_pd(a, b); }

MD_SIMD_INLINE md_128 md_mm_div(md_128 a, md_128 b) { return md_simd_div_ps(a, b); }
MD_SIMD_INLINE md_256 md_mm256_div(md_256 a, md_256 b) { return md_simd_div_ps(a, b); }
MD_SIMD_INLINE md_128d md_simd_div(md_128d a, md_128d b) { return md_simd_div_pd(a, b); }
MD_SIMD_INLINE md_256d md_simd_div(md_256d a, md_256d b) { return md_simd_div_pd(a, b); }

MD_SIMD_INLINE md_128i md_simd_cast(md_128 v) { return md_simd_cast_ps(v); }
MD_SIMD_INLINE md_256i md_simd_cast(md_256 v) { return md_simd_cast_ps(v); }
MD_SIMD_INLINE md_128i md_simd_cast(md_128d v) { return md_simd_cast_pd(v); }
MD_SIMD_INLINE md_256i md_simd_cast(md_256d v) { return md_simd_cast_pd(v); }

MD_SIMD_INLINE md_128 md_mm_cast(md_128i v) { return md_simd_cast_i32x4(v); }
MD_SIMD_INLINE md_256 md_mm256_cast(md_256i v) { return md_simd_cast_i32x8(v); }
MD_SIMD_INLINE md_128d md_simd_cast(md_128i v) { return md_simd_cast_i64x2(v); }
MD_SIMD_INLINE md_256d md_simd_cast(md_256i v) { return md_simd_cast_i64x4(v); }

MD_SIMD_INLINE md_128i md_simd_convert(md_128 v) { return md_simd_convert_ps(v); }
MD_SIMD_INLINE md_256i md_simd_convert(md_256 v) { return md_simd_convert_ps(v); }
MD_SIMD_INLINE md_128i md_simd_convert(md_128d v) { return md_simd_convert_pd(v); }
MD_SIMD_INLINE md_256i md_simd_convert(md_256d v) { return md_simd_convert_pd(v); }

MD_SIMD_INLINE md_128 md_mm_convert(md_128i v) { return md_simd_convert_i32x4(v); }
MD_SIMD_INLINE md_256 md_mm256_convert(md_256i v) { return md_simd_convert_i32x8(v); }
MD_SIMD_INLINE md_128d md_simd_convert(md_128i v) { return md_simd_convert_i64x2(v); }
MD_SIMD_INLINE md_256d md_simd_convert(md_256i v) { return md_simd_convert_i64x4(v); }

// @TODO: Complete this

#else
// C11 generics

#define md_simd_store(a, b) _Generic((b),   \
            md_128i : md_simd_store_i32x4, \
            md_256i : md_simd_store_i32x8, \
            md_128i : md_simd_store_i64x2, \
            md_256i : md_simd_store_i64x4, \
            const md_128i : md_simd_store_i32x4, \
            const md_256i : md_simd_store_i32x8, \
            const md_128i : md_simd_store_i64x2, \
            const md_256i : md_simd_store_i64x4, \
            md_128 : md_simd_store_ps, \
            md_256 : md_simd_store_ps, \
            md_128d : md_simd_store_pd, \
            md_256d : md_simd_store_pd, \
            const md_128 : md_simd_store_ps, \
            const md_256 : md_simd_store_ps, \
            const md_128d : md_simd_store_pd, \
            const md_256d : md_simd_store_pd)(a,b)

#define md_simd_and(a, b) _Generic((a),     \
            md_128i : md_simd_and_i32x4, \
            md_256i : md_simd_and_i32x8, \
            md_128i : md_simd_and_i64x2, \
            md_256i : md_simd_and_i64x4, \
            const md_128i : md_simd_and_i32x4, \
            const md_256i : md_simd_and_i32x8, \
            const md_128i : md_simd_and_i64x2, \
            const md_256i : md_simd_and_i64x4, \
            md_128 : md_simd_and_ps, \
            md_256 : md_simd_and_ps, \
            md_128d : md_simd_and_pd, \
            md_256d : md_simd_and_pd, \
            const md_128 : md_simd_and_ps, \
            const md_256 : md_simd_and_ps, \
            const md_128d : md_simd_and_pd, \
            const md_256d : md_simd_and_pd)(a,b)

#define md_simd_or(a, b) _Generic((a),      \
            md_128i : md_simd_or_i32x4,  \
            md_256i : md_simd_or_i32x8,  \
            md_128i : md_simd_or_i64x2,  \
            md_256i : md_simd_or_i64x4,  \
            const md_128i : md_simd_or_i32x4,  \
            const md_256i : md_simd_or_i32x8,  \
            const md_128i : md_simd_or_i64x2,  \
            const md_256i : md_simd_or_i64x4,  \
            md_128 : md_simd_or_ps,  \
            md_256 : md_simd_or_ps,  \
            md_128d : md_simd_or_pd,  \
            md_256d : md_simd_or_pd,  \
            const md_128 : md_simd_or_ps,  \
            const md_256 : md_simd_or_ps,  \
            const md_128d : md_simd_or_pd,  \
            const md_256d : md_simd_or_pd)(a,b)

#define md_simd_xor(a, b) _Generic((a),     \
            md_128i : md_simd_xor_i32x4, \
            md_256i : md_simd_xor_i32x8, \
            md_128i : md_simd_xor_i64x2, \
            md_256i : md_simd_xor_i64x4, \
            const md_128i : md_simd_xor_i32x4, \
            const md_256i : md_simd_xor_i32x8, \
            const md_128i : md_simd_xor_i64x2, \
            const md_256i : md_simd_xor_i64x4, \
            md_128 : md_simd_xor_ps, \
            md_256 : md_simd_xor_ps, \
            md_128d : md_simd_xor_pd, \
            md_256d : md_simd_xor_pd, \
            const md_128 : md_simd_xor_ps, \
            const md_256 : md_simd_xor_ps, \
            const md_128d : md_simd_xor_pd, \
            const md_256d : md_simd_xor_pd)(a,b)

#define md_simd_and_not(a, b) _Generic((a),     \
            md_128i : md_simd_and_not_i32x4, \
            md_256i : md_simd_and_not_i32x8, \
            md_128i : md_simd_and_not_i64x2, \
            md_256i : md_simd_and_not_i64x4, \
            const md_128i : md_simd_and_not_i32x4, \
            const md_256i : md_simd_and_not_i32x8, \
            const md_128i : md_simd_and_not_i64x2, \
            const md_256i : md_simd_and_not_i64x4, \
            md_128 : md_simd_and_not_ps, \
            md_256 : md_simd_and_not_ps, \
            md_128d : md_simd_and_not_pd, \
            md_256d : md_simd_and_not_pd, \
            const md_128 : md_simd_and_not_ps, \
            const md_256 : md_simd_and_not_ps, \
            const md_128d : md_simd_and_not_pd, \
            const md_256d : md_simd_and_not_pd)(a,b)

#define md_simd_not(x) _Generic((x),        \
            md_128i : md_simd_not_i32x4, \
            md_256i : md_simd_not_i32x8, \
            md_128i : md_simd_not_i64x2, \
            md_256i : md_simd_not_i64x4, \
            const md_128i : md_simd_not_i32x4, \
            const md_256i : md_simd_not_i32x8, \
            const md_128i : md_simd_not_i64x2, \
            const md_256i : md_simd_not_i64x4)(x)

#define md_simd_add(a, b) _Generic((a),     \
            md_128 : md_simd_add_ps, \
            md_256 : md_simd_add_ps, \
            md_128d : md_simd_add_pd, \
            md_256d : md_simd_add_pd, \
            const md_128 : md_simd_add_ps, \
            const md_256 : md_simd_add_ps, \
            const md_128d : md_simd_add_pd, \
            const md_256d : md_simd_add_pd, \
            md_128i : md_simd_add_i32x4, \
            md_256i : md_simd_add_i32x8, \
            md_128i : md_simd_add_i64x2, \
            md_256i : md_simd_add_i64x4, \
            const md_128i : md_simd_add_i32x4, \
            const md_256i : md_simd_add_i32x8, \
            const md_128i : md_simd_add_i64x2, \
            const md_256i : md_simd_add_i64x4)(a, b)
            
#define md_simd_sub(a, b) _Generic((a),     \
            md_128 : md_simd_sub_ps, \
            md_256 : md_simd_sub_ps, \
            md_128d : md_simd_sub_pd, \
            md_256d : md_simd_sub_pd, \
            const md_128 : md_simd_sub_ps, \
            const md_256 : md_simd_sub_ps, \
            const md_128d : md_simd_sub_pd, \
            const md_256d : md_simd_sub_pd, \
            md_128i : md_simd_sub_i32x4, \
            md_256i : md_simd_sub_i32x8, \
            md_128i : md_simd_sub_i64x2, \
            md_256i : md_simd_sub_i64x4, \
            const md_128i : md_simd_sub_i32x4, \
            const md_256i : md_simd_sub_i32x8, \
            const md_128i : md_simd_sub_i64x2, \
            const md_256i : md_simd_sub_i64x4)(a, b)

#define md_simd_mul(a, b) _Generic((a),     \
            md_128 : md_simd_mul_ps, \
            md_256 : md_simd_mul_ps, \
            md_128d : md_simd_mul_pd, \
            md_256d : md_simd_mul_pd, \
            const md_128 : md_simd_mul_ps, \
            const md_256 : md_simd_mul_ps, \
            const md_128d : md_simd_mul_pd, \
            const md_256d : md_simd_mul_pd)(a, b)

#define md_simd_div(a, b) _Generic((a),     \
            md_128 : md_simd_div_ps, \
            md_256 : md_simd_div_ps, \
            md_128d : md_simd_div_pd, \
            md_256d : md_simd_div_pd, \
            const md_128 : md_simd_div_ps, \
            const md_256 : md_simd_div_ps, \
            const md_128d : md_simd_div_pd, \
            const md_256d : md_simd_div_pd)(a, b)

#define md_simd_fmadd(a, b, c) _Generic((a),       \
            md_128 : md_simd_fmadd_ps,   \
            md_256 : md_simd_fmadd_ps,   \
            md_128d : md_simd_fmadd_pd,   \
            md_256d : md_simd_fmadd_pd,   \
            const md_128 : md_simd_fmadd_ps,   \
            const md_256 : md_simd_fmadd_ps,   \
            const md_128d : md_simd_fmadd_pd,   \
            const md_256d : md_simd_fmadd_pd)(a, b, c)

#define md_simd_abs(x) _Generic((x),        \
            md_128 : md_simd_abs_ps, \
            md_256 : md_simd_abs_ps, \
            md_128d : md_simd_abs_pd, \
            md_256d : md_simd_abs_pd, \
            const md_128 : md_simd_abs_ps, \
            const md_256 : md_simd_abs_ps, \
            const md_128d : md_simd_abs_pd, \
            const md_256d : md_simd_abs_pd, \
            md_128i : md_simd_abs_i32x4, \
            md_256i : md_simd_abs_i32x8, \
            md_128i : md_simd_abs_i64x2, \
            md_256i : md_simd_abs_i64x4, \
            const md_128i : md_simd_abs_i32x4, \
            const md_256i : md_simd_abs_i32x8, \
            const md_128i : md_simd_abs_i64x2, \
            const md_256i : md_simd_abs_i64x4)(x)

#define md_simd_min(a,b) _Generic((a),      \
            md_128 : md_simd_min_ps, \
            md_256 : md_simd_min_ps, \
            md_128d : md_simd_min_pd, \
            md_256d : md_simd_min_pd, \
            const md_128 : md_simd_min_ps, \
            const md_256 : md_simd_min_ps, \
            const md_128d : md_simd_min_pd, \
            const md_256d : md_simd_min_pd, \
            md_128i : md_simd_min_i32x4, \
            md_256i : md_simd_min_i32x8, \
            md_128i : md_simd_min_i64x2, \
            md_256i : md_simd_min_i64x4, \
            const md_128i : md_simd_min_i32x4, \
            const md_256i : md_simd_min_i32x8, \
            const md_128i : md_simd_min_i64x2, \
            const md_256i : md_simd_min_i64x4)(a,b)

#define md_simd_max(a,b) _Generic((a),      \
            md_128 : md_simd_max_ps, \
            md_256 : md_simd_max_ps, \
            md_128d : md_simd_max_pd, \
            md_256d : md_simd_max_pd, \
            const md_128 : md_simd_max_ps, \
            const md_256 : md_simd_max_ps, \
            const md_128d : md_simd_max_pd, \
            const md_256d : md_simd_max_pd, \
            md_128i : md_simd_max_i32x4, \
            md_256i : md_simd_max_i32x8, \
            md_128i : md_simd_max_i64x2, \
            md_256i : md_simd_max_i64x4, \
            const md_128i : md_simd_max_i32x4, \
            const md_256i : md_simd_max_i32x8, \
            const md_128i : md_simd_max_i64x2, \
            const md_256i : md_simd_max_i64x4)(a,b)

#define md_simd_cmp_gt(a,b) _Generic((a),       \
            md_128 : md_simd_cmp_gt_ps,  \
            md_256 : md_simd_cmp_gt_ps,  \
            md_128d : md_simd_cmp_gt_pd,  \
            md_256d : md_simd_cmp_gt_pd,  \
            const md_128 : md_simd_cmp_gt_ps,  \
            const md_256 : md_simd_cmp_gt_ps,  \
            const md_128d : md_simd_cmp_gt_pd,  \
            const md_256d : md_simd_cmp_gt_pd,  \
            md_128i : md_simd_cmp_gt_i32x4,  \
            md_256i : md_simd_cmp_gt_i32x8,  \
            md_128i : md_simd_cmp_gt_i64x2,  \
            md_256i : md_simd_cmp_gt_i64x4,  \
            const md_128i : md_simd_cmp_gt_i32x4,  \
            const md_256i : md_simd_cmp_gt_i32x8,  \
            const md_128i : md_simd_cmp_gt_i64x2,  \
            const md_256i : md_simd_cmp_gt_i64x4)(a,b)

#define md_simd_cmp_lt(a,b) _Generic((a),         \
            md_128 : md_simd_cmp_lt_ps,  \
            md_256 : md_simd_cmp_lt_ps,  \
            md_128d : md_simd_cmp_lt_pd,  \
            md_256d : md_simd_cmp_lt_pd,  \
            const md_128 : md_simd_cmp_lt_ps,  \
            const md_256 : md_simd_cmp_lt_ps,  \
            const md_128d : md_simd_cmp_lt_pd,  \
            const md_256d : md_simd_cmp_lt_pd,  \
            md_128i : md_simd_cmp_lt_i32x4,  \
            md_256i : md_simd_cmp_lt_i32x8,  \
            md_128i : md_simd_cmp_lt_i64x2,  \
            md_256i : md_simd_cmp_lt_i64x4,  \
            const md_128i : md_simd_cmp_lt_i32x4,  \
            const md_256i : md_simd_cmp_lt_i32x8,  \
            const md_128i : md_simd_cmp_lt_i64x2,  \
            const md_256i : md_simd_cmp_lt_i64x4)(a,b)

#define md_simd_cmp_eq(a,b) _Generic((a),       \
            md_128 : md_simd_cmp_eq_ps,  \
            md_256 : md_simd_cmp_eq_ps,  \
            md_128d : md_simd_cmp_eq_pd,  \
            md_256d : md_simd_cmp_eq_pd,  \
            const md_128 : md_simd_cmp_eq_ps,  \
            const md_256 : md_simd_cmp_eq_ps,  \
            const md_128d : md_simd_cmp_eq_pd,  \
            const md_256d : md_simd_cmp_eq_pd,  \
            md_128i : md_simd_cmp_eq_i32x4,  \
            md_256i : md_simd_cmp_eq_i32x8,  \
            md_128i : md_simd_cmp_eq_i64x2,  \
            md_256i : md_simd_cmp_eq_i64x4,  \
            const md_128i : md_simd_cmp_eq_i32x4,  \
            const md_256i : md_simd_cmp_eq_i32x8,  \
            const md_128i : md_simd_cmp_eq_i64x2,  \
            const md_256i : md_simd_cmp_eq_i64x4)(a,b)

#define md_simd_cmp_ne(a,b) _Generic((a),       \
            md_128 : md_simd_cmp_ne_ps,  \
            md_256 : md_simd_cmp_ne_ps,  \
            md_128d : md_simd_cmp_ne_pd,  \
            md_256d : md_simd_cmp_ne_pd,  \
            const md_128 : md_simd_cmp_ne_ps,  \
            const md_256 : md_simd_cmp_ne_ps,  \
            const md_128d : md_simd_cmp_ne_pd,  \
            const md_256d : md_simd_cmp_ne_pd,  \
            md_128i : md_simd_cmp_ne_i32x4,  \
            md_256i : md_simd_cmp_ne_i32x8,  \
            md_128i : md_simd_cmp_ne_i64x2,  \
            md_256i : md_simd_cmp_ne_i64x4,  \
            const md_128i : md_simd_cmp_ne_i32x4,  \
            const md_256i : md_simd_cmp_ne_i32x8,  \
            const md_128i : md_simd_cmp_ne_i64x2,  \
            const md_256i : md_simd_cmp_ne_i64x4)(a,b)

#define md_simd_cmp_le(a,b) _Generic((a),       \
            md_128 : md_simd_cmp_le_ps,  \
            md_256 : md_simd_cmp_le_ps,  \
            md_128d : md_simd_cmp_le_pd,  \
            md_256d : md_simd_cmp_le_pd,  \
            const md_128 : md_simd_cmp_le_ps,  \
            const md_256 : md_simd_cmp_le_ps,  \
            const md_128d : md_simd_cmp_le_pd,  \
            const md_256d : md_simd_cmp_le_pd)(a,b)

#define md_simd_cmp_le(a,b) _Generic((a),       \
            md_128 : md_simd_cmp_le_ps,  \
            md_256 : md_simd_cmp_le_ps,  \
            md_128d : md_simd_cmp_le_pd,  \
            md_256d : md_simd_cmp_le_pd,  \
            const md_128 : md_simd_cmp_le_ps,  \
            const md_256 : md_simd_cmp_le_ps,  \
            const md_128d : md_simd_cmp_le_pd,  \
            const md_256d : md_simd_cmp_le_pd)(a,b)

#define md_simd_round(x) _Generic((x),          \
            md_128 : md_simd_round_ps,   \
            md_256 : md_simd_round_ps,   \
            md_128d : md_simd_round_pd,   \
            md_256d : md_simd_round_pd,   \
            const md_128 : md_simd_round_ps,   \
            const md_256 : md_simd_round_ps,   \
            const md_128d : md_simd_round_pd,   \
            const md_256d : md_simd_round_pd)(x)

#define md_simd_floor(x) _Generic((x),          \
            md_128 : md_simd_floor_ps,   \
            md_256 : md_simd_floor_ps,   \
            md_128d : md_simd_floor_pd,   \
            md_256d : md_simd_floor_pd,   \
            const md_128 : md_simd_floor_ps,   \
            const md_256 : md_simd_floor_ps,   \
            const md_128d : md_simd_floor_pd,   \
            const md_256d : md_simd_floor_pd)(x)

#define md_simd_ceil(x) _Generic((x),           \
            md_128 : md_simd_ceil_ps,    \
            md_256 : md_simd_ceil_ps,    \
            md_128d : md_simd_ceil_pd,    \
            md_256d : md_simd_ceil_pd,    \
            const md_128 : md_simd_ceil_ps,    \
            const md_256 : md_simd_ceil_ps,    \
            const md_128d : md_simd_ceil_pd,    \
            const md_256d : md_simd_ceil_pd)(x)

#define md_simd_fract(x) _Generic((x),          \
            md_128 : md_simd_fract_ps,   \
            md_256 : md_simd_fract_ps,   \
            md_128d : md_simd_fract_pd,   \
            md_256d : md_simd_fract_pd,   \
            const md_128 : md_simd_fract_ps,   \
            const md_256 : md_simd_fract_ps,   \
            const md_128d : md_simd_fract_pd,   \
            const md_256d : md_simd_fract_pd)(x)

#define md_simd_sign(x) _Generic((x),           \
            md_128 : md_simd_sign_ps,    \
            md_256 : md_simd_sign_ps,    \
            md_128d : md_simd_sign_pd,    \
            md_256d : md_simd_sign_pd)(x)

#define md_simd_sqrt(x) _Generic((x),           \
            md_128 : md_simd_sqrt_ps,    \
            md_256 : md_simd_sqrt_ps,    \
            md_128d : md_simd_sqrt_pd,    \
            md_256d : md_simd_sqrt_pd,    \
            const md_128 : md_simd_sqrt_ps,    \
            const md_256 : md_simd_sqrt_ps,    \
            const md_128d : md_simd_sqrt_pd,    \
            const md_256d : md_simd_sqrt_pd)(x)

#define md_simd_blend(a,b,mask) _Generic((a),   \
            md_128 : md_simd_blend_ps,   \
            md_256 : md_simd_blend_ps,   \
            md_128d : md_simd_blend_pd,   \
            md_256d : md_simd_blend_pd,   \
            const md_128 : md_simd_blend_ps,   \
            const md_256 : md_simd_blend_ps,   \
            const md_128d : md_simd_blend_pd,   \
            const md_256d : md_simd_blend_pd,   \
            md_128i : md_simd_blend_i32x4,   \
            md_256i : md_simd_blend_i32x8,   \
            md_128i : md_simd_blend_i64x2,   \
            md_256i : md_simd_blend_i64x4,   \
            const md_128i : md_simd_blend_i32x4,   \
            const md_256i : md_simd_blend_i32x8,   \
            const md_128i : md_simd_blend_i64x2,   \
            const md_256i : md_simd_blend_i64x4)(a,b,mask)

#define md_simd_movemask(a) _Generic((a),       \
            md_128 : md_simd_movemask_ps,\
            md_256 : md_simd_movemask_ps,\
            md_128d : md_simd_movemask_pd,\
            md_256d : md_simd_movemask_pd,\
            const md_128 : md_simd_movemask_ps,\
            const md_256 : md_simd_movemask_ps,\
            const md_128d : md_simd_movemask_pd,\
            const md_256d : md_simd_movemask_pd)(a)

#define md_simd_hmin(x) _Generic((x),           \
            md_128 : md_simd_reduce_min_ps,    \
            md_256 : md_simd_reduce_min_ps,    \
            md_128d : md_simd_reduce_min_pd,    \
            md_256d : md_simd_reduce_min_pd,    \
            const md_128 : md_simd_reduce_min_ps,    \
            const md_256 : md_simd_reduce_min_ps,    \
            const md_128d : md_simd_reduce_min_pd,    \
            const md_256d : md_simd_reduce_min_pd)(x)

#define md_simd_hmax(x) _Generic((x),           \
            md_128 : md_simd_reduce_max_ps,    \
            md_256 : md_simd_reduce_max_ps,    \
            md_128d : md_simd_reduce_max_pd,    \
            md_256d : md_simd_reduce_max_pd,    \
            const md_128 : md_simd_reduce_max_ps,    \
            const md_256 : md_simd_reduce_max_ps,    \
            const md_128d : md_simd_reduce_max_pd,    \
            const md_256d : md_simd_reduce_max_pd)(x)

#define md_simd_hsum(x) _Generic((x),           \
            md_128 : md_simd_reduce_add_ps,    \
            md_256 : md_simd_reduce_add_ps,    \
            md_128d : md_simd_reduce_add_pd,    \
            md_256d : md_simd_reduce_add_pd,    \
            const md_128 : md_simd_reduce_add_ps,    \
            const md_256 : md_simd_reduce_add_ps,    \
            const md_128d : md_simd_reduce_add_pd,    \
            const md_256d : md_simd_reduce_add_pd)(x)

#define md_simd_deperiodize(x, r, p) _Generic((x),  \
            md_128 : md_simd_deperiodize_ps, \
            md_256 : md_simd_deperiodize_ps, \
            md_128d : md_simd_deperiodize_pd, \
            md_256d : md_simd_deperiodize_pd, \
            const md_128 : md_simd_deperiodize_ps, \
            const md_256 : md_simd_deperiodize_ps, \
            const md_128d : md_simd_deperiodize_pd, \
            const md_256d : md_simd_deperiodize_pd)(x, r, p)

#define md_simd_minimum_image(dx, p, rp) _Generic((dx),  \
            md_128 : md_simd_minimum_image_ps, \
            md_256 : md_simd_minimum_image_ps, \
            md_128d : md_simd_minimum_image_pd, \
            md_256d : md_simd_minimum_image_pd, \
            const md_128 : md_simd_minimum_image_ps, \
            const md_256 : md_simd_minimum_image_ps, \
            const md_128d : md_simd_minimum_image_pd, \
            const md_256d : md_simd_minimum_image_pd)(dx, p, rp)

#define md_simd_step(edge, x) _Generic((x),  \
            md_128 : md_simd_step_ps, \
            md_256 : md_simd_step_ps, \
            md_128d : md_simd_step_pd, \
            md_256d : md_simd_step_pd, \
            const md_128 : md_simd_step_ps, \
            const md_256 : md_simd_step_ps, \
            const md_128d : md_simd_step_pd, \
            const md_256d : md_simd_step_pd)(edge, x)

#define md_simd_lerp(a, b, t) _Generic((a),  \
            md_128 : md_simd_lerp_ps, \
            md_256 : md_simd_lerp_ps, \
            md_128d : md_simd_lerp_pd, \
            md_256d : md_simd_lerp_pd, \
            const md_128 : md_simd_lerp_ps, \
            const md_256 : md_simd_lerp_ps, \
            const md_128d : md_simd_lerp_pd, \
            const md_256d : md_simd_lerp_pd)(a, b, t)

#define md_simd_cubic_spline(p0, p1, p2, p3, t, s) _Generic((p0),  \
            md_128 : md_simd_cubic_spline_ps, \
            md_256 : md_simd_cubic_spline_ps, \
            md_128d : md_simd_cubic_spline_pd, \
            md_256d : md_simd_cubic_spline_pd, \
            const md_128 : md_simd_cubic_spline_ps, \
            const md_256 : md_simd_cubic_spline_ps, \
            const md_128d : md_simd_cubic_spline_pd, \
            const md_256d : md_simd_cubic_spline_pd)(p0, p1, p2, p3, t, s)

#define md_simd_unpack_xyz(x, y, z, stream, stride) _Generic(x,  \
            md_128* : md_mm_unpack_xyz_ps, \
            md_256* : md_mm256_unpack_xyz_ps)(x, y, z, stream, stride)

#define md_simd_sincos(x, s, c) _Generic(x,  \
            md_128 : md_simd_sincos_ps, \
            md_256 : md_simd_sincos_ps)(x, s, c)

#define md_simd_shift_left(x, i) _Generic((x),      \
            md_128i : md_simd_shift_left_i32x4,  \
            md_256i : md_simd_shift_left_i32x8,  \
            md_128i : md_simd_shift_left_i64x2,  \
            md_256i : md_simd_shift_left_i64x4,  \
            const md_128i : md_simd_shift_left_i32x4,  \
            const md_256i : md_simd_shift_left_i32x8,  \
            const md_128i : md_simd_shift_left_i64x2,  \
            const md_256i : md_simd_shift_left_i64x4)(x, i)

#define md_simd_shift_right(x, i) _Generic((x),     \
            md_128i : md_simd_shift_right_i32x4, \
            md_256i : md_simd_shift_right_i32x8, \
            md_128i : md_simd_shift_right_i64x2, \
            md_256i : md_simd_shift_right_i64x4, \
            const md_128i : md_simd_shift_right_i32x4, \
            const md_256i : md_simd_shift_right_i32x8, \
            const md_128i : md_simd_shift_right_i64x2, \
            const md_256i : md_simd_shift_right_i64x4)(x, i)

#endif

#if md_simd_width_f32 == 8

// Float
#define md_mm256_loadu_ps    md_simd_load_ps
#define md_simd_load_f64    md_simd_load_pd

#define md_simd_store_f32   md_simd_store_ps
#define md_simd_store_f64   md_simd_store_pd

#define md_simd_set1_f32    md_simd_set1_ps
#define md_simd_set1_f64    md_simd_set1_pd

#define md_simd_set_f32     md_simd_set_ps
#define md_simd_set_f64     md_simd_set_pd

#define md_simd_zero_f32    md_simd_zero_ps
#define md_simd_zero_f64    md_simd_zero_pd

#define md_simd_convert_f32 md_simd_convert_ps
#define md_simd_convert_f64 md_simd_convert_pd

#define md_simd_cast_f32 md_simd_cast_ps
#define md_simd_cast_f64 md_simd_cast_pd

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
#define md_mm256_loadu_ps    md_simd_load_ps
#define md_simd_load_f64    md_simd_load_pd

#define md_simd_store_f32   md_simd_store_ps
#define md_simd_store_f64   md_simd_store_pd

#define md_simd_set1_f32    md_simd_set1_ps
#define md_simd_set1_f64    md_simd_set1_pd

#define md_simd_set_f32     md_simd_set_ps
#define md_simd_set_f64     md_simd_set_pd

#define md_simd_zero_f32    md_simd_zero_ps
#define md_simd_zero_f64    md_simd_zero_pd

#define md_simd_convert_f32 md_simd_convert_ps
#define md_simd_convert_f64 md_simd_convert_pd

#define md_simd_cast_f32 md_simd_cast_ps
#define md_simd_cast_f64 md_simd_cast_pd

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

#endif
