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

#define MD_SHUFFLE SIMDE_MM_SHUFFLE

// NATIVE PATHS
//
// simde gives the x86 semantics everywhere. Where NEON has no instruction with those semantics it
// composes one: a horizontal sum from shuffles, a min that returns the second operand for a NaN, a
// truncating conversion that gives INT_MIN for a NaN or an out of range lane. md_simd takes the
// NEON instruction instead where the x86 corner is not something mdlib relies on, and documents the
// difference where it does so. Every other target, x86 included, goes through simde as before.
//
// On AArch64 (MD_SIMD_NEON):
//   md_*_min_*, md_*_max_*   fminnm and fmaxnm: a NaN in either operand gives the other operand, where
//                            x86 gives the second operand whichever is the NaN.
//   md_*_reduce_*            across the vector in one instruction; min and max ignore NaN lanes.
//   md_*_cvttps_epi32        out of range lanes saturate and a NaN gives 0, where x86 gives INT_MIN.
//   md_*_cvtps_epi32         rounds to nearest even, which is the x86 default rounding mode.
//   md_*_rsqrt_ps            one Newton-Raphson step on the 8 bit estimate of NEON, so it is at least as
//                            accurate as the 12 bit estimate of x86 (which callers refine further).
//   md_mm256_movemask_ps     the two halves of the 128 bit form.
//   md_*_fnmadd_*, fmsub_*,  fused, as on x86 (simde multiplies and subtracts in places, which the compiler
//   fnmsub_*                 may or may not fuse: MSVC does not).
#if defined(SIMDE_ARM_NEON_A64V8_NATIVE)
#define MD_SIMD_NEON 1
#endif

#if defined(MD_SIMD_NEON)
#define MD_NEON_F32(x)  simde__m128_to_neon_f32(x)
#define MD_NEON_F64(x)  simde__m128d_to_neon_f64(x)
#define MD_NEON_I32(x)  simde__m128i_to_neon_i32(x)
#define MD_FROM_F32(x)  simde__m128_from_neon_f32(x)
#define MD_FROM_F64(x)  simde__m128d_from_neon_f64(x)
#define MD_FROM_I32(x)  simde__m128i_from_neon_i32(x)

// An eight (four) lane operation as two four (two) lane ones: simde holds the 256 bit types as two 128 bit
// halves on NEON, so the split and the join cost nothing. Each argument is used twice: arguments are names.
#define MD_NEON_256_PS_1(op, a)     simde_mm256_set_m128(op(simde_mm256_extractf128_ps(a, 1)), op(simde_mm256_castps256_ps128(a)))
#define MD_NEON_256_PS_2(op, a, b)  simde_mm256_set_m128(op(simde_mm256_extractf128_ps(a, 1), simde_mm256_extractf128_ps(b, 1)), op(simde_mm256_castps256_ps128(a), simde_mm256_castps256_ps128(b)))
#define MD_NEON_256_PD_2(op, a, b)  simde_mm256_set_m128d(op(simde_mm256_extractf128_pd(a, 1), simde_mm256_extractf128_pd(b, 1)), op(simde_mm256_castpd256_pd128(a), simde_mm256_castpd256_pd128(b)))
#define MD_NEON_256_PS_EPI32(op, a) simde_mm256_set_m128i(op(simde_mm256_extractf128_ps(a, 1)), op(simde_mm256_castps256_ps128(a)))
#endif

static const uint32_t sse_compress_lut[16] = {
    0xFFFFFFFF,
    0xFFFFFF00,
    0xFFFFFF01,
    0xFFFF0100,
    0xFFFFFF02,
    0xFFFF0200,
    0xFFFF0201,
    0xFF020100,
    0xFFFFFF03,
    0xFFFF0300,
    0xFFFF0301,
    0xFF030100,
    0xFFFF0302,
    0xFF030200,
    0xFF030201,
    0x03020100,
};

static const uint64_t avx_compress_lut[256] = {
    0x0000000000000000ULL,   0x0000000000000000ULL,   0x0000000000000001ULL,   0x0000000000000100ULL,
    0x0000000000000002ULL,   0x0000000000000200ULL,   0x0000000000000201ULL,   0x0000000000020100ULL,
    0x0000000000000003ULL,   0x0000000000000300ULL,   0x0000000000000301ULL,   0x0000000000030100ULL,
    0x0000000000000302ULL,   0x0000000000030200ULL,   0x0000000000030201ULL,   0x0000000003020100ULL,
    0x0000000000000004ULL,   0x0000000000000400ULL,   0x0000000000000401ULL,   0x0000000000040100ULL,
    0x0000000000000402ULL,   0x0000000000040200ULL,   0x0000000000040201ULL,   0x0000000004020100ULL,
    0x0000000000000403ULL,   0x0000000000040300ULL,   0x0000000000040301ULL,   0x0000000004030100ULL,
    0x0000000000040302ULL,   0x0000000004030200ULL,   0x0000000004030201ULL,   0x0000000403020100ULL,
    0x0000000000000005ULL,   0x0000000000000500ULL,   0x0000000000000501ULL,   0x0000000000050100ULL,
    0x0000000000000502ULL,   0x0000000000050200ULL,   0x0000000000050201ULL,   0x0000000005020100ULL,
    0x0000000000000503ULL,   0x0000000000050300ULL,   0x0000000000050301ULL,   0x0000000005030100ULL,
    0x0000000000050302ULL,   0x0000000005030200ULL,   0x0000000005030201ULL,   0x0000000503020100ULL,
    0x0000000000000504ULL,   0x0000000000050400ULL,   0x0000000000050401ULL,   0x0000000005040100ULL,
    0x0000000000050402ULL,   0x0000000005040200ULL,   0x0000000005040201ULL,   0x0000000504020100ULL,
    0x0000000000050403ULL,   0x0000000005040300ULL,   0x0000000005040301ULL,   0x0000000504030100ULL,
    0x0000000005040302ULL,   0x0000000504030200ULL,   0x0000000504030201ULL,   0x0000050403020100ULL,
    0x0000000000000006ULL,   0x0000000000000600ULL,   0x0000000000000601ULL,   0x0000000000060100ULL,
    0x0000000000000602ULL,   0x0000000000060200ULL,   0x0000000000060201ULL,   0x0000000006020100ULL,
    0x0000000000000603ULL,   0x0000000000060300ULL,   0x0000000000060301ULL,   0x0000000006030100ULL,
    0x0000000000060302ULL,   0x0000000006030200ULL,   0x0000000006030201ULL,   0x0000000603020100ULL,
    0x0000000000000604ULL,   0x0000000000060400ULL,   0x0000000000060401ULL,   0x0000000006040100ULL,
    0x0000000000060402ULL,   0x0000000006040200ULL,   0x0000000006040201ULL,   0x0000000604020100ULL,
    0x0000000000060403ULL,   0x0000000006040300ULL,   0x0000000006040301ULL,   0x0000000604030100ULL,
    0x0000000006040302ULL,   0x0000000604030200ULL,   0x0000000604030201ULL,   0x0000060403020100ULL,
    0x0000000000000605ULL,   0x0000000000060500ULL,   0x0000000000060501ULL,   0x0000000006050100ULL,
    0x0000000000060502ULL,   0x0000000006050200ULL,   0x0000000006050201ULL,   0x0000000605020100ULL,
    0x0000000000060503ULL,   0x0000000006050300ULL,   0x0000000006050301ULL,   0x0000000605030100ULL,
    0x0000000006050302ULL,   0x0000000605030200ULL,   0x0000000605030201ULL,   0x0000060503020100ULL,
    0x0000000000060504ULL,   0x0000000006050400ULL,   0x0000000006050401ULL,   0x0000000605040100ULL,
    0x0000000006050402ULL,   0x0000000605040200ULL,   0x0000000605040201ULL,   0x0000060504020100ULL,
    0x0000000006050403ULL,   0x0000000605040300ULL,   0x0000000605040301ULL,   0x0000060504030100ULL,
    0x0000000605040302ULL,   0x0000060504030200ULL,   0x0000060504030201ULL,   0x0006050403020100ULL,
    0x0000000000000007ULL,   0x0000000000000700ULL,   0x0000000000000701ULL,   0x0000000000070100ULL,
    0x0000000000000702ULL,   0x0000000000070200ULL,   0x0000000000070201ULL,   0x0000000007020100ULL,
    0x0000000000000703ULL,   0x0000000000070300ULL,   0x0000000000070301ULL,   0x0000000007030100ULL,
    0x0000000000070302ULL,   0x0000000007030200ULL,   0x0000000007030201ULL,   0x0000000703020100ULL,
    0x0000000000000704ULL,   0x0000000000070400ULL,   0x0000000000070401ULL,   0x0000000007040100ULL,
    0x0000000000070402ULL,   0x0000000007040200ULL,   0x0000000007040201ULL,   0x0000000704020100ULL,
    0x0000000000070403ULL,   0x0000000007040300ULL,   0x0000000007040301ULL,   0x0000000704030100ULL,
    0x0000000007040302ULL,   0x0000000704030200ULL,   0x0000000704030201ULL,   0x0000070403020100ULL,
    0x0000000000000705ULL,   0x0000000000070500ULL,   0x0000000000070501ULL,   0x0000000007050100ULL,
    0x0000000000070502ULL,   0x0000000007050200ULL,   0x0000000007050201ULL,   0x0000000705020100ULL,
    0x0000000000070503ULL,   0x0000000007050300ULL,   0x0000000007050301ULL,   0x0000000705030100ULL,
    0x0000000007050302ULL,   0x0000000705030200ULL,   0x0000000705030201ULL,   0x0000070503020100ULL,
    0x0000000000070504ULL,   0x0000000007050400ULL,   0x0000000007050401ULL,   0x0000000705040100ULL,
    0x0000000007050402ULL,   0x0000000705040200ULL,   0x0000000705040201ULL,   0x0000070504020100ULL,
    0x0000000007050403ULL,   0x0000000705040300ULL,   0x0000000705040301ULL,   0x0000070504030100ULL,
    0x0000000705040302ULL,   0x0000070504030200ULL,   0x0000070504030201ULL,   0x0007050403020100ULL,
    0x0000000000000706ULL,   0x0000000000070600ULL,   0x0000000000070601ULL,   0x0000000007060100ULL,
    0x0000000000070602ULL,   0x0000000007060200ULL,   0x0000000007060201ULL,   0x0000000706020100ULL,
    0x0000000000070603ULL,   0x0000000007060300ULL,   0x0000000007060301ULL,   0x0000000706030100ULL,
    0x0000000007060302ULL,   0x0000000706030200ULL,   0x0000000706030201ULL,   0x0000070603020100ULL,
    0x0000000000070604ULL,   0x0000000007060400ULL,   0x0000000007060401ULL,   0x0000000706040100ULL,
    0x0000000007060402ULL,   0x0000000706040200ULL,   0x0000000706040201ULL,   0x0000070604020100ULL,
    0x0000000007060403ULL,   0x0000000706040300ULL,   0x0000000706040301ULL,   0x0000070604030100ULL,
    0x0000000706040302ULL,   0x0000070604030200ULL,   0x0000070604030201ULL,   0x0007060403020100ULL,
    0x0000000000070605ULL,   0x0000000007060500ULL,   0x0000000007060501ULL,   0x0000000706050100ULL,
    0x0000000007060502ULL,   0x0000000706050200ULL,   0x0000000706050201ULL,   0x0000070605020100ULL,
    0x0000000007060503ULL,   0x0000000706050300ULL,   0x0000000706050301ULL,   0x0000070605030100ULL,
    0x0000000706050302ULL,   0x0000070605030200ULL,   0x0000070605030201ULL,   0x0007060503020100ULL,
    0x0000000007060504ULL,   0x0000000706050400ULL,   0x0000000706050401ULL,   0x0000070605040100ULL,
    0x0000000706050402ULL,   0x0000070605040200ULL,   0x0000070605040201ULL,   0x0007060504020100ULL,
    0x0000000706050403ULL,   0x0000070605040300ULL,   0x0000070605040301ULL,   0x0007060504030100ULL,
    0x0000070605040302ULL,   0x0007060504030200ULL,   0x0007060504030201ULL,   0x0706050403020100ULL,
};

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

// Non-temporal: for output written once and not read soon. Weakly ordered on x86, so a writer ends
// with md_mm_sfence before the values are handed to another thread.
#define md_mm_stream_pd simde_mm_stream_pd
#define md_mm_sfence simde_mm_sfence

#define md_mm_maskstore_ps simde_mm_maskstore_ps
#define md_mm_maskstore_pd simde_mm_maskstore_pd
#define md_mm256_maskstore_ps simde_mm256_maskstore_ps
#define md_mm256_maskstore_pd simde_mm256_maskstore_pd

#define md_mm_set1_ps simde_mm_set1_ps
#define md_mm_set1_pd simde_mm_set1_pd
#define md_mm256_set1_ps simde_mm256_set1_ps
#define md_mm256_setr_ps simde_mm256_setr_ps
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
#define md_mm256_cmp_ps simde_mm256_cmp_ps
#define md_mm256_cmp_pd simde_mm256_cmp_pd

#define md_mm_cmpeq_ps simde_mm_cmpeq_ps
#define md_mm_cmpeq_pd simde_mm_cmpeq_pd
#define md_mm256_cmpeq_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_EQ_OQ)
#define md_mm256_cmpeq_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_EQ_OQ)

#define md_mm_cmpneq_ps simde_mm_cmpneq_ps
#define md_mm_cmpneq_pd simde_mm_cmpneq_pd
#define md_mm256_cmpneq_ps(a,b) simde_mm256_cmp_ps(a,b, SIMDE_CMP_NEQ_UQ)    // As _mm_cmpneq_ps: true for a NaN
#define md_mm256_cmpneq_pd(a,b) simde_mm256_cmp_pd(a,b, SIMDE_CMP_NEQ_UQ)

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

#if defined(MD_SIMD_NEON)
// simde composes these from a multiply and a subtract on NEON, which the compiler may or may not fuse: one rounding
// here as on x86. -(a * b) + c, a * b - c, and -(a * b) - c. For fnmadd simde has the fused instruction at four lanes
// but builds the eight lane float form from scalar code, which clang contracts and MSVC does not.
MD_SIMD_INLINE md_128  md_mm_fnmadd_ps(md_128 a, md_128 b, md_128 c)    { return MD_FROM_F32(vfmsq_f32(MD_NEON_F32(c), MD_NEON_F32(a), MD_NEON_F32(b))); }
MD_SIMD_INLINE md_128d md_mm_fnmadd_pd(md_128d a, md_128d b, md_128d c) { return MD_FROM_F64(vfmsq_f64(MD_NEON_F64(c), MD_NEON_F64(a), MD_NEON_F64(b))); }
MD_SIMD_INLINE md_256  md_mm256_fnmadd_ps(md_256 a, md_256 b, md_256 c) {
    return simde_mm256_set_m128(md_mm_fnmadd_ps(simde_mm256_extractf128_ps(a, 1), simde_mm256_extractf128_ps(b, 1), simde_mm256_extractf128_ps(c, 1)),
                                md_mm_fnmadd_ps(simde_mm256_castps256_ps128(a), simde_mm256_castps256_ps128(b), simde_mm256_castps256_ps128(c)));
}
MD_SIMD_INLINE md_256d md_mm256_fnmadd_pd(md_256d a, md_256d b, md_256d c) {
    return simde_mm256_set_m128d(md_mm_fnmadd_pd(simde_mm256_extractf128_pd(a, 1), simde_mm256_extractf128_pd(b, 1), simde_mm256_extractf128_pd(c, 1)),
                                 md_mm_fnmadd_pd(simde_mm256_castpd256_pd128(a), simde_mm256_castpd256_pd128(b), simde_mm256_castpd256_pd128(c)));
}
MD_SIMD_INLINE md_128  md_mm_fmsub_ps(md_128 a, md_128 b, md_128 c)     { return MD_FROM_F32(vfmaq_f32(vnegq_f32(MD_NEON_F32(c)), MD_NEON_F32(a), MD_NEON_F32(b))); }
MD_SIMD_INLINE md_128d md_mm_fmsub_pd(md_128d a, md_128d b, md_128d c)  { return MD_FROM_F64(vfmaq_f64(vnegq_f64(MD_NEON_F64(c)), MD_NEON_F64(a), MD_NEON_F64(b))); }
MD_SIMD_INLINE md_128  md_mm_fnmsub_ps(md_128 a, md_128 b, md_128 c)    { return MD_FROM_F32(vfmsq_f32(vnegq_f32(MD_NEON_F32(c)), MD_NEON_F32(a), MD_NEON_F32(b))); }
MD_SIMD_INLINE md_128d md_mm_fnmsub_pd(md_128d a, md_128d b, md_128d c) { return MD_FROM_F64(vfmsq_f64(vnegq_f64(MD_NEON_F64(c)), MD_NEON_F64(a), MD_NEON_F64(b))); }
MD_SIMD_INLINE md_256  md_mm256_fmsub_ps(md_256 a, md_256 b, md_256 c) {
    return simde_mm256_set_m128(md_mm_fmsub_ps(simde_mm256_extractf128_ps(a, 1), simde_mm256_extractf128_ps(b, 1), simde_mm256_extractf128_ps(c, 1)),
                                md_mm_fmsub_ps(simde_mm256_castps256_ps128(a), simde_mm256_castps256_ps128(b), simde_mm256_castps256_ps128(c)));
}
MD_SIMD_INLINE md_256d md_mm256_fmsub_pd(md_256d a, md_256d b, md_256d c) {
    return simde_mm256_set_m128d(md_mm_fmsub_pd(simde_mm256_extractf128_pd(a, 1), simde_mm256_extractf128_pd(b, 1), simde_mm256_extractf128_pd(c, 1)),
                                 md_mm_fmsub_pd(simde_mm256_castpd256_pd128(a), simde_mm256_castpd256_pd128(b), simde_mm256_castpd256_pd128(c)));
}
MD_SIMD_INLINE md_256  md_mm256_fnmsub_ps(md_256 a, md_256 b, md_256 c) {
    return simde_mm256_set_m128(md_mm_fnmsub_ps(simde_mm256_extractf128_ps(a, 1), simde_mm256_extractf128_ps(b, 1), simde_mm256_extractf128_ps(c, 1)),
                                md_mm_fnmsub_ps(simde_mm256_castps256_ps128(a), simde_mm256_castps256_ps128(b), simde_mm256_castps256_ps128(c)));
}
MD_SIMD_INLINE md_256d md_mm256_fnmsub_pd(md_256d a, md_256d b, md_256d c) {
    return simde_mm256_set_m128d(md_mm_fnmsub_pd(simde_mm256_extractf128_pd(a, 1), simde_mm256_extractf128_pd(b, 1), simde_mm256_extractf128_pd(c, 1)),
                                 md_mm_fnmsub_pd(simde_mm256_castpd256_pd128(a), simde_mm256_castpd256_pd128(b), simde_mm256_castpd256_pd128(c)));
}
#else
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
#endif

#define md_mm_blendv_ps simde_mm_blendv_ps
#define md_mm_blendv_pd simde_mm_blendv_pd
#define md_mm256_blendv_ps simde_mm256_blendv_ps
#define md_mm256_blendv_pd simde_mm256_blendv_pd

// a where the mask is clear and b where it is set, for a mask whose lanes are all ones or all zeros, as a
// comparison gives them. blendv decides by the sign bit alone, one instruction on x86 and two on NEON (the
// sign spread over the lane, then the select); select is one on both.
#if defined(MD_SIMD_NEON)
MD_SIMD_INLINE md_128  md_mm_select_ps(md_128 a, md_128 b, md_128 mask) {
    return MD_FROM_F32(vbslq_f32(vreinterpretq_u32_f32(MD_NEON_F32(mask)), MD_NEON_F32(b), MD_NEON_F32(a)));
}
MD_SIMD_INLINE md_128d md_mm_select_pd(md_128d a, md_128d b, md_128d mask) {
    return MD_FROM_F64(vbslq_f64(vreinterpretq_u64_f64(MD_NEON_F64(mask)), MD_NEON_F64(b), MD_NEON_F64(a)));
}
MD_SIMD_INLINE md_256  md_mm256_select_ps(md_256 a, md_256 b, md_256 mask) {
    return simde_mm256_set_m128(md_mm_select_ps(simde_mm256_extractf128_ps(a, 1), simde_mm256_extractf128_ps(b, 1), simde_mm256_extractf128_ps(mask, 1)),
                                md_mm_select_ps(simde_mm256_castps256_ps128(a), simde_mm256_castps256_ps128(b), simde_mm256_castps256_ps128(mask)));
}
MD_SIMD_INLINE md_256d md_mm256_select_pd(md_256d a, md_256d b, md_256d mask) {
    return simde_mm256_set_m128d(md_mm_select_pd(simde_mm256_extractf128_pd(a, 1), simde_mm256_extractf128_pd(b, 1), simde_mm256_extractf128_pd(mask, 1)),
                                 md_mm_select_pd(simde_mm256_castpd256_pd128(a), simde_mm256_castpd256_pd128(b), simde_mm256_castpd256_pd128(mask)));
}
#else
#define md_mm_select_ps    simde_mm_blendv_ps
#define md_mm_select_pd    simde_mm_blendv_pd
#define md_mm256_select_ps simde_mm256_blendv_ps
#define md_mm256_select_pd simde_mm256_blendv_pd
#endif

#define MD_SIMD_BLEND_MASK(x,y,z,w) (((w) << 3) | ((z) << 2) | ((y) << 1) | (x))

#define md_mm_blend_ps simde_mm_blend_ps
#define md_mm_blend_pd simde_mm_blend_pd
#define md_mm256_blend_ps simde_mm256_blend_ps
#define md_mm256_blend_pd simde_mm256_blend_pd

#define md_mm_sqrt_ps simde_mm_sqrt_ps
#define md_mm_sqrt_pd simde_mm_sqrt_pd
#define md_mm256_sqrt_ps simde_mm256_sqrt_ps
#define md_mm256_sqrt_pd simde_mm256_sqrt_pd

#if defined(MD_SIMD_NEON)
MD_SIMD_INLINE md_128 md_mm_rsqrt_ps(md_128 a) {
    const float32x4_t x = MD_NEON_F32(a);
    const float32x4_t e = vrsqrteq_f32(x);
    return MD_FROM_F32(vmulq_f32(e, vrsqrtsq_f32(vmulq_f32(x, e), e)));
}
MD_SIMD_INLINE md_256 md_mm256_rsqrt_ps(md_256 a) { return MD_NEON_256_PS_1(md_mm_rsqrt_ps, a); }
#else
#define md_mm_rsqrt_ps simde_mm_rsqrt_ps
#define md_mm256_rsqrt_ps simde_mm256_rsqrt_ps
#endif

#if defined(MD_SIMD_NEON)
MD_SIMD_INLINE md_128  md_mm_min_ps(md_128 a, md_128 b)   { return MD_FROM_F32(vminnmq_f32(MD_NEON_F32(a), MD_NEON_F32(b))); }
MD_SIMD_INLINE md_128d md_mm_min_pd(md_128d a, md_128d b) { return MD_FROM_F64(vminnmq_f64(MD_NEON_F64(a), MD_NEON_F64(b))); }
MD_SIMD_INLINE md_128  md_mm_max_ps(md_128 a, md_128 b)   { return MD_FROM_F32(vmaxnmq_f32(MD_NEON_F32(a), MD_NEON_F32(b))); }
MD_SIMD_INLINE md_128d md_mm_max_pd(md_128d a, md_128d b) { return MD_FROM_F64(vmaxnmq_f64(MD_NEON_F64(a), MD_NEON_F64(b))); }
MD_SIMD_INLINE md_256  md_mm256_min_ps(md_256 a, md_256 b)   { return MD_NEON_256_PS_2(md_mm_min_ps, a, b); }
MD_SIMD_INLINE md_256d md_mm256_min_pd(md_256d a, md_256d b) { return MD_NEON_256_PD_2(md_mm_min_pd, a, b); }
MD_SIMD_INLINE md_256  md_mm256_max_ps(md_256 a, md_256 b)   { return MD_NEON_256_PS_2(md_mm_max_ps, a, b); }
MD_SIMD_INLINE md_256d md_mm256_max_pd(md_256d a, md_256d b) { return MD_NEON_256_PD_2(md_mm_max_pd, a, b); }
#else
#define md_mm_min_ps simde_mm_min_ps
#define md_mm_min_pd simde_mm_min_pd
#define md_mm256_min_ps simde_mm256_min_ps
#define md_mm256_min_pd simde_mm256_min_pd

#define md_mm_max_ps simde_mm_max_ps
#define md_mm_max_pd simde_mm_max_pd
#define md_mm256_max_ps simde_mm256_max_ps
#define md_mm256_max_pd simde_mm256_max_pd
#endif

#define md_mm_shuffle_ps simde_mm_shuffle_ps
#define md_mm_shuffle_pd simde_mm_shuffle_pd
#define md_mm256_shuffle_ps simde_mm256_shuffle_ps
#define md_mm256_shuffle_pd simde_mm256_shuffle_pd

#define md_mm_movehl_ps simde_mm_movehl_ps
#define md_mm_movemask_ps simde_mm_movemask_ps
#define md_mm_movemask_pd simde_mm_movemask_pd
#define md_mm256_movemask_pd simde_mm256_movemask_pd
#if defined(MD_SIMD_NEON)
// The sign bits of the eight lanes, the low half in bits 0-3 and the high in 4-7: each lane's bit shifted to
// its place and the lanes summed, both halves in one sum
MD_SIMD_INLINE int md_mm256_movemask_ps(md_256 a) {
    static const int32_t shift_lo[4] = { 0, 1, 2, 3 };
    static const int32_t shift_hi[4] = { 4, 5, 6, 7 };
    const uint32x4_t lo = vshrq_n_u32(vreinterpretq_u32_f32(MD_NEON_F32(simde_mm256_castps256_ps128(a))), 31);
    const uint32x4_t hi = vshrq_n_u32(vreinterpretq_u32_f32(MD_NEON_F32(simde_mm256_extractf128_ps(a, 1))), 31);
    return (int)vaddvq_u32(vaddq_u32(vshlq_u32(lo, vld1q_s32(shift_lo)), vshlq_u32(hi, vld1q_s32(shift_hi))));
}
#else
#define md_mm256_movemask_ps simde_mm256_movemask_ps
#endif

// BASE INT OPERATIONS

#define md_mm_loadu_si128 simde_mm_loadu_si128
#define md_mm_loadu_epi32 simde_mm_loadu_epi32
#define md_mm_loadu_epi64 simde_mm_loadu_epi64
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

#define md_mm_set1_epi8  simde_mm_set1_epi8
#define md_mm_set1_epi16 simde_mm_set1_epi16
#define md_mm_set1_epi32 simde_mm_set1_epi32
#define md_mm_set1_epi64 simde_mm_set1_epi64x
#define md_mm256_set1_epi32 simde_mm256_set1_epi32
#define md_mm256_set1_epi64 simde_mm256_set1_epi64x

#define md_mm_set_epi8  simde_mm_set_epi8
#define md_mm_set_epi16 simde_mm_set_epi16
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

#define md_mm_cmpeq_epi8  simde_mm_cmpeq_epi8
#define md_mm_cmpeq_epi16 simde_mm_cmpeq_epi16
#define md_mm_cmpeq_epi32 simde_mm_cmpeq_epi32
#define md_mm_cmpeq_epi64 simde_mm_cmpeq_epi64

#define md_mm256_cmpeq_epi8  simde_mm256_cmpeq_epi8
#define md_mm256_cmpeq_epi16 simde_mm256_cmpeq_epi16
#define md_mm256_cmpeq_epi32 simde_mm256_cmpeq_epi32
#define md_mm256_cmpeq_epi64 simde_mm256_cmpeq_epi64

#define md_mm_cmpneq_epi32(a,b) md_mm_not_si128(md_mm_cmpeq_epi32(a,b))
#define md_mm_cmpneq_epi64(a,b) md_mm_not_si128(md_mm_cmpeq_epi64(a,b))
#define md_mm256_cmpneq_epi32(a,b) md_mm256_not_si256(md_mm256_cmpeq_epi32(a,b))
#define md_mm256_cmpneq_epi64(a,b) md_mm256_not_si256(md_mm256_cmpeq_epi64(a,b))

#define md_mm_cmpgt_epi8  simde_mm_cmpgt_epi8
#define md_mm_cmpgt_epi16 simde_mm_cmpgt_epi16
#define md_mm_cmpgt_epi32 simde_mm_cmpgt_epi32
#define md_mm_cmpgt_epi64 simde_mm_cmpgt_epi64

#define md_mm256_cmpgt_epi8  simde_mm256_cmpgt_epi8
#define md_mm256_cmpgt_epi16 simde_mm256_cmpgt_epi16
#define md_mm256_cmpgt_epi32 simde_mm256_cmpgt_epi32
#define md_mm256_cmpgt_epi64 simde_mm256_cmpgt_epi64

#define md_mm_cmplt_epi8  simde_mm_cmplt_epi8
#define md_mm_cmplt_epi16 simde_mm_cmplt_epi16
#define md_mm_cmplt_epi32 simde_mm_cmplt_epi32
#define md_mm_cmplt_epi64(a,b) simde_mm_cmpgt_epi64(b,a)

#define md_mm256_cmplt_epi8(a,b)  simde_mm256_cmpgt_epi8(b,a)
#define md_mm256_cmplt_epi16(a,b) simde_mm256_cmpgt_epi16(b,a)
#define md_mm256_cmplt_epi32(a,b) simde_mm256_cmpgt_epi32(b,a)
#define md_mm256_cmplt_epi64(a,b) simde_mm256_cmpgt_epi64(b,a)

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

#define md_mm_cvtepi32_ps simde_mm_cvtepi32_ps
#define md_mm256_cvtepi32_ps simde_mm256_cvtepi32_ps

#if defined(MD_SIMD_NEON)
MD_SIMD_INLINE md_128i md_mm_cvtps_epi32(md_128 a)  { return MD_FROM_I32(vcvtnq_s32_f32(MD_NEON_F32(a))); }
MD_SIMD_INLINE md_128i md_mm_cvttps_epi32(md_128 a) { return MD_FROM_I32(vcvtq_s32_f32(MD_NEON_F32(a))); }
MD_SIMD_INLINE md_256i md_mm256_cvtps_epi32(md_256 a)  { return MD_NEON_256_PS_EPI32(md_mm_cvtps_epi32, a); }
MD_SIMD_INLINE md_256i md_mm256_cvttps_epi32(md_256 a) { return MD_NEON_256_PS_EPI32(md_mm_cvttps_epi32, a); }
#else
#define md_mm_cvtps_epi32 simde_mm_cvtps_epi32
#define md_mm_cvttps_epi32 simde_mm_cvttps_epi32
#define md_mm256_cvtps_epi32 simde_mm256_cvtps_epi32
#define md_mm256_cvttps_epi32 simde_mm256_cvttps_epi32
#endif
#define md_mm_cvtepu8_epi32 simde_mm_cvtepu8_epi32
#define md_mm_cvtsi64_si128 simde_mm_cvtsi64_si128
#define md_mm256_cvtepu8_epi32 simde_mm256_cvtepu8_epi32

#define md_mm_cvtss_f32 simde_mm_cvtss_f32
#define md_mm_cvtsd_f64 simde_mm_cvtsd_f64
#define md_mm256_cvtss_f32 simde_mm256_cvtss_f32
#define md_mm256_cvtsd_f64 simde_mm256_cvtsd_f64

// PERMUTE
#define md_mm256_permutevar8x32_epi32 simde_mm256_permutevar8x32_epi32
#define md_mm256_permutevar8x32_ps simde_mm256_permutevar8x32_ps

// HIGH LEVEL OPERATIONS

#define MD_LOAD_STRIDED_128(ptr, offset, stride) simde_mm_loadu_ps((const float*)((const char*)(ptr) + (size_t)(offset) * (stride)))

// Four packed xyz triplets - twelve contiguous floats, x0 y0 z0 x1 y1 ... z3 - split into x, y and z.
// It reads exactly the twelve floats, so the last four atoms of an array are as safe to load as any
// others. The strided variant below reads 16 bytes per atom and so needs the array padded past its
// end.
//
// On NEON this is one instruction, a de-interleaving load. The x86 path below - three loads and
// five shuffles - is native on x86 and what simde would otherwise translate everywhere else.
MD_SIMD_INLINE void md_mm_load_xyz_packed_ps(md_128* out_x, md_128* out_y, md_128* out_z, const float* in_xyz) {
#if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
    const float32x4x3_t v = vld3q_f32(in_xyz);
    *out_x = simde__m128_from_neon_f32(v.val[0]);
    *out_y = simde__m128_from_neon_f32(v.val[1]);
    *out_z = simde__m128_from_neon_f32(v.val[2]);
#else
    const md_128 v0 = simde_mm_loadu_ps(in_xyz + 0);  // x0 y0 z0 x1
    const md_128 v1 = simde_mm_loadu_ps(in_xyz + 4);  // y1 z1 x2 y2
    const md_128 v2 = simde_mm_loadu_ps(in_xyz + 8);  // z2 x3 y3 z3
    const md_128 t0 = simde_mm_shuffle_ps(v1, v2, SIMDE_MM_SHUFFLE(2, 1, 3, 2));  // x2 y2 x3 y3
    const md_128 t1 = simde_mm_shuffle_ps(v0, v1, SIMDE_MM_SHUFFLE(1, 0, 2, 1));  // y0 z0 y1 z1
    *out_x = simde_mm_shuffle_ps(v0, t0, SIMDE_MM_SHUFFLE(2, 0, 3, 0));           // x0 x1 x2 x3
    *out_y = simde_mm_shuffle_ps(t1, t0, SIMDE_MM_SHUFFLE(3, 1, 2, 0));           // y0 y1 y2 y3
    *out_z = simde_mm_shuffle_ps(t1, v2, SIMDE_MM_SHUFFLE(3, 0, 3, 1));           // z0 z1 z2 z3
#endif
}

// The inverse: four x, y and z back into packed triplets. Writes exactly the twelve floats.
MD_SIMD_INLINE void md_mm_store_xyz_packed_ps(float* out_xyz, md_128 x, md_128 y, md_128 z) {
#if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
    float32x4x3_t v;
    v.val[0] = simde__m128_to_neon_f32(x);
    v.val[1] = simde__m128_to_neon_f32(y);
    v.val[2] = simde__m128_to_neon_f32(z);
    vst3q_f32(out_xyz, v);
#else
    const md_128 a = simde_mm_shuffle_ps(x, y, SIMDE_MM_SHUFFLE(2, 0, 2, 0));  // x0 x2 y0 y2
    const md_128 b = simde_mm_shuffle_ps(y, z, SIMDE_MM_SHUFFLE(3, 1, 3, 1));  // y1 y3 z1 z3
    const md_128 c = simde_mm_shuffle_ps(z, x, SIMDE_MM_SHUFFLE(3, 1, 2, 0));  // z0 z2 x1 x3
    simde_mm_storeu_ps(out_xyz + 0, simde_mm_shuffle_ps(a, c, SIMDE_MM_SHUFFLE(2, 0, 2, 0)));  // x0 y0 z0 x1
    simde_mm_storeu_ps(out_xyz + 4, simde_mm_shuffle_ps(b, a, SIMDE_MM_SHUFFLE(3, 1, 2, 0)));  // y1 z1 x2 y2
    simde_mm_storeu_ps(out_xyz + 8, simde_mm_shuffle_ps(c, b, SIMDE_MM_SHUFFLE(3, 1, 3, 1)));  // z2 x3 y3 z3
#endif
}

// Eight packed triplets, 24 floats. On x86 each 128 bit lane does the four atom split above at once,
// atoms 0-3 in the low lane and 4-7 in the high; on NEON it is two de-interleaving loads.
MD_SIMD_INLINE void md_mm256_load_xyz_packed_ps(md_256* out_x, md_256* out_y, md_256* out_z, const float* in_xyz) {
#if defined(SIMDE_ARM_NEON_A32V7_NATIVE)
    md_128 x0, y0, z0, x1, y1, z1;
    md_mm_load_xyz_packed_ps(&x0, &y0, &z0, in_xyz);
    md_mm_load_xyz_packed_ps(&x1, &y1, &z1, in_xyz + 12);
    *out_x = simde_mm256_set_m128(x1, x0);
    *out_y = simde_mm256_set_m128(y1, y0);
    *out_z = simde_mm256_set_m128(z1, z0);
#else
    const md_256 v0 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(simde_mm_loadu_ps(in_xyz + 0)), simde_mm_loadu_ps(in_xyz + 12), 1);
    const md_256 v1 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(simde_mm_loadu_ps(in_xyz + 4)), simde_mm_loadu_ps(in_xyz + 16), 1);
    const md_256 v2 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(simde_mm_loadu_ps(in_xyz + 8)), simde_mm_loadu_ps(in_xyz + 20), 1);
    const md_256 t0 = simde_mm256_shuffle_ps(v1, v2, SIMDE_MM_SHUFFLE(2, 1, 3, 2));
    const md_256 t1 = simde_mm256_shuffle_ps(v0, v1, SIMDE_MM_SHUFFLE(1, 0, 2, 1));
    *out_x = simde_mm256_shuffle_ps(v0, t0, SIMDE_MM_SHUFFLE(2, 0, 3, 0));
    *out_y = simde_mm256_shuffle_ps(t1, t0, SIMDE_MM_SHUFFLE(3, 1, 2, 0));
    *out_z = simde_mm256_shuffle_ps(t1, v2, SIMDE_MM_SHUFFLE(3, 0, 3, 1));
#endif
}

MD_SIMD_INLINE void md_mm256_store_xyz_packed_ps(float* out_xyz, md_256 x, md_256 y, md_256 z) {
    md_mm_store_xyz_packed_ps(out_xyz,      simde_mm256_castps256_ps128(x), simde_mm256_castps256_ps128(y), simde_mm256_castps256_ps128(z));
    md_mm_store_xyz_packed_ps(out_xyz + 12, simde_mm256_extractf128_ps(x, 1), simde_mm256_extractf128_ps(y, 1), simde_mm256_extractf128_ps(z, 1));
}

// x, y and z of the atoms idx names in packed triplets: the same three gathers as from separate
// arrays, at three times the index.
MD_SIMD_INLINE void md_mm_i32gather_xyz_ps(md_128* out_x, md_128* out_y, md_128* out_z, const float* in_xyz, md_128i idx) {
    const md_128i i3 = simde_mm_add_epi32(idx, simde_mm_add_epi32(idx, idx));
    *out_x = simde_mm_i32gather_ps(in_xyz + 0, i3, 4);
    *out_y = simde_mm_i32gather_ps(in_xyz + 1, i3, 4);
    *out_z = simde_mm_i32gather_ps(in_xyz + 2, i3, 4);
}

MD_SIMD_INLINE void md_mm256_i32gather_xyz_ps(md_256* out_x, md_256* out_y, md_256* out_z, const float* in_xyz, md_256i idx) {
    const md_256i i3 = simde_mm256_add_epi32(idx, simde_mm256_add_epi32(idx, idx));
    *out_x = simde_mm256_i32gather_ps(in_xyz + 0, i3, 4);
    *out_y = simde_mm256_i32gather_ps(in_xyz + 1, i3, 4);
    *out_z = simde_mm256_i32gather_ps(in_xyz + 2, i3, 4);
}

#if defined(__AVX512F__)
MD_SIMD_INLINE __m512 md_mm512_from_m256_ps(__m256 lo, __m256 hi) {
    return _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_castps256_ps512(lo)), _mm256_castps_pd(hi), 1));
}

// Sixteen packed triplets, 48 floats
MD_SIMD_INLINE void md_mm512_load_xyz_packed_ps(__m512* out_x, __m512* out_y, __m512* out_z, const float* in_xyz) {
    md_256 x0, y0, z0, x1, y1, z1;
    md_mm256_load_xyz_packed_ps(&x0, &y0, &z0, in_xyz);
    md_mm256_load_xyz_packed_ps(&x1, &y1, &z1, in_xyz + 24);
    *out_x = md_mm512_from_m256_ps(x0, x1);
    *out_y = md_mm512_from_m256_ps(y0, y1);
    *out_z = md_mm512_from_m256_ps(z0, z1);
}

MD_SIMD_INLINE void md_mm512_i32gather_xyz_ps(__m512* out_x, __m512* out_y, __m512* out_z, const float* in_xyz, __m512i idx) {
    const __m512i i3 = _mm512_add_epi32(idx, _mm512_add_epi32(idx, idx));
    *out_x = _mm512_i32gather_ps(i3, in_xyz + 0, 4);
    *out_y = _mm512_i32gather_ps(i3, in_xyz + 1, 4);
    *out_z = _mm512_i32gather_ps(i3, in_xyz + 2, 4);
}
#endif

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

MD_SIMD_INLINE void md_mm_unpack_xyzw_ps(md_128* out_x, md_128* out_y, md_128* out_z, md_128* out_w, const float* in_xyzw, size_t stride_in_bytes) {
    md_128 r0, r1, r2, r3;
    md_128 t0, t1, t2, t3;

    r0 = MD_LOAD_STRIDED_128(in_xyzw, 0, stride_in_bytes);
    r1 = MD_LOAD_STRIDED_128(in_xyzw, 1, stride_in_bytes);
    r2 = MD_LOAD_STRIDED_128(in_xyzw, 2, stride_in_bytes);
    r3 = MD_LOAD_STRIDED_128(in_xyzw, 3, stride_in_bytes);

    t0 = simde_mm_unpacklo_ps(r0, r1);  // xxyy xxyy
    t1 = simde_mm_unpackhi_ps(r0, r1);  // zzww zzww
    t2 = simde_mm_unpacklo_ps(r2, r3);  // xxyy xxyy
    t3 = simde_mm_unpackhi_ps(r2, r3);  // zzww zzww

    *out_x = simde_mm_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(1, 0, 1, 0));  // xxxx xxxx
    *out_y = simde_mm_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(3, 2, 3, 2));  // yyyy yyyy
    *out_z = simde_mm_shuffle_ps(t1, t3, SIMDE_MM_SHUFFLE(1, 0, 1, 0));  // zzzz zzzz
    *out_w = simde_mm_shuffle_ps(t1, t3, SIMDE_MM_SHUFFLE(3, 2, 3, 2));  // wwww wwww
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

MD_SIMD_INLINE void md_mm256_unpack_xyzw_ps(md_256* out_x, md_256* out_y, md_256* out_z, md_256* out_w, const float* in_xyzw, size_t stride_in_bytes) {
    md_256 r0, r1, r2, r3;
    md_256 t0, t1, t2, t3;

    // @TODO: Try and implement this using 256-bit loads and then shuffle all the way.
    // It's a tradeoff between the number of loads issued and the port pressure on the generally few ports that are used for shuffle instructions.
    r0 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyzw, 0, stride_in_bytes)),
                                   MD_LOAD_STRIDED_128(in_xyzw, 4, stride_in_bytes), 1);
    r1 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyzw, 1, stride_in_bytes)),
                                   MD_LOAD_STRIDED_128(in_xyzw, 5, stride_in_bytes), 1);
    r2 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyzw, 2, stride_in_bytes)),
                                   MD_LOAD_STRIDED_128(in_xyzw, 6, stride_in_bytes), 1);
    r3 = simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(MD_LOAD_STRIDED_128(in_xyzw, 3, stride_in_bytes)),
                                   MD_LOAD_STRIDED_128(in_xyzw, 7, stride_in_bytes), 1);

    t0 = simde_mm256_unpacklo_ps(r0, r1);  // xxyy xxyy
    t1 = simde_mm256_unpackhi_ps(r0, r1);  // zzww zzww
    t2 = simde_mm256_unpacklo_ps(r2, r3);  // xxyy xxyy
    t3 = simde_mm256_unpackhi_ps(r2, r3);  // zzww zzww

    *out_x = simde_mm256_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(1, 0, 1, 0));  // xxxx xxxx
    *out_y = simde_mm256_shuffle_ps(t0, t2, SIMDE_MM_SHUFFLE(3, 2, 3, 2));  // yyyy yyyy
    *out_z = simde_mm256_shuffle_ps(t1, t3, SIMDE_MM_SHUFFLE(1, 0, 1, 0));  // zzzz zzzz
    *out_w = simde_mm256_shuffle_ps(t1, t3, SIMDE_MM_SHUFFLE(3, 2, 3, 2));  // wwww wwww
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

// -a: the sign bit flipped
MD_SIMD_INLINE md_128   md_mm_neg_ps(md_128 a) { return simde_mm_xor_ps(a, simde_mm_set1_ps(-0.0f)); }
MD_SIMD_INLINE md_128d  md_mm_neg_pd(md_128d a) { return simde_mm_xor_pd(a, simde_mm_set1_pd(-0.0)); }
MD_SIMD_INLINE md_256   md_mm256_neg_ps(md_256 a) { return simde_mm256_xor_ps(a, simde_mm256_set1_ps(-0.0f)); }
MD_SIMD_INLINE md_256d  md_mm256_neg_pd(md_256d a) { return simde_mm256_xor_pd(a, simde_mm256_set1_pd(-0.0)); }


// HORIZONTAL REDUCTIONS
//
// The sums add in the same order on every target, ((x0 + x1) + (x2 + x3)) for four lanes and the two halves
// first for eight, so a sum is the same bit for bit on x86 and on NEON.

MD_SIMD_INLINE float md_mm_reduce_add_ps(md_128 x) {
#if defined(MD_SIMD_NEON)
    return vaddvq_f32(MD_NEON_F32(x));
#else
    // https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
    md_128 shuf = md_mm_shuffle_ps(x, x, SIMDE_MM_SHUFFLE(3, 3, 1, 1));    // x1 x1 x3 x3 (movshdup, which compilers pick when SSE3 is there)
    md_128 sums = md_mm_add_ps(x, shuf);
    shuf        = md_mm_movehl_ps(shuf, sums);                              // high half -> low half
    sums        = md_mm_add_ps(sums, shuf);
    return        md_mm_cvtss_f32(sums);
#endif
}

MD_SIMD_INLINE float md_mm_reduce_min_ps(md_128 x) {
#if defined(MD_SIMD_NEON)
    return vminnmvq_f32(MD_NEON_F32(x));
#else
    const md_128 a = simde_mm_min_ps(x, simde_mm_shuffle_ps(x, x, SIMDE_MM_SHUFFLE(0, 0, 3, 2)));
    const md_128 b = simde_mm_min_ps(a, simde_mm_shuffle_ps(a, a, SIMDE_MM_SHUFFLE(0, 0, 0, 1)));
    return simde_mm_cvtss_f32(b);
#endif
}

MD_SIMD_INLINE float md_mm_reduce_max_ps(md_128 x) {
#if defined(MD_SIMD_NEON)
    return vmaxnmvq_f32(MD_NEON_F32(x));
#else
    const md_128 a = simde_mm_max_ps(x, simde_mm_shuffle_ps(x, x, SIMDE_MM_SHUFFLE(0, 0, 3, 2)));
    const md_128 b = simde_mm_max_ps(a, simde_mm_shuffle_ps(a, a, SIMDE_MM_SHUFFLE(0, 0, 0, 1)));
    return simde_mm_cvtss_f32(b);
#endif
}

#if defined(MD_SIMD_NEON)
MD_SIMD_INLINE double md_mm_reduce_add_pd(md_128d x) { return vaddvq_f64(MD_NEON_F64(x)); }
MD_SIMD_INLINE double md_mm_reduce_min_pd(md_128d x) { return vminnmvq_f64(MD_NEON_F64(x)); }
MD_SIMD_INLINE double md_mm_reduce_max_pd(md_128d x) { return vmaxnmvq_f64(MD_NEON_F64(x)); }
#else
MD_SIMD_INLINE double md_mm_reduce_add_pd(md_128d x) { return simde_mm_cvtsd_f64(simde_mm_add_pd(x, md_mm_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 0, 1)))); }
MD_SIMD_INLINE double md_mm_reduce_min_pd(md_128d x) { return simde_mm_cvtsd_f64(simde_mm_min_pd(x, md_mm_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 0, 1)))); }
MD_SIMD_INLINE double md_mm_reduce_max_pd(md_128d x) { return simde_mm_cvtsd_f64(simde_mm_max_pd(x, md_mm_shuffle_pd(x, x, SIMDE_MM_SHUFFLE(0, 0, 0, 1)))); }
#endif

// Eight (four) lanes: the two halves combined, then reduced as four (two)
MD_SIMD_INLINE float  md_mm256_reduce_add_ps(md_256 x)  { return md_mm_reduce_add_ps(md_mm_add_ps(simde_mm256_castps256_ps128(x), simde_mm256_extractf128_ps(x, 1))); }
MD_SIMD_INLINE float  md_mm256_reduce_min_ps(md_256 x)  { return md_mm_reduce_min_ps(md_mm_min_ps(simde_mm256_castps256_ps128(x), simde_mm256_extractf128_ps(x, 1))); }
MD_SIMD_INLINE float  md_mm256_reduce_max_ps(md_256 x)  { return md_mm_reduce_max_ps(md_mm_max_ps(simde_mm256_castps256_ps128(x), simde_mm256_extractf128_ps(x, 1))); }
MD_SIMD_INLINE double md_mm256_reduce_add_pd(md_256d x) { return md_mm_reduce_add_pd(md_mm_add_pd(simde_mm256_castpd256_pd128(x), simde_mm256_extractf128_pd(x, 1))); }
MD_SIMD_INLINE double md_mm256_reduce_min_pd(md_256d x) { return md_mm_reduce_min_pd(md_mm_min_pd(simde_mm256_castpd256_pd128(x), simde_mm256_extractf128_pd(x, 1))); }
MD_SIMD_INLINE double md_mm256_reduce_max_pd(md_256d x) { return md_mm_reduce_max_pd(md_mm_max_pd(simde_mm256_castpd256_pd128(x), simde_mm256_extractf128_pd(x, 1))); }

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
//
// The reciprocal flavor forms the number of periods d = x/p - r/p with a fused multiply-subtract. With a plain product
// (x - r) * rp, a compiler that contracts may fuse it into d - round(d) even though round(d) uses it too (Clang does
// for Apple CPUs): round(d) then sees the rounded product and the subtraction the exact one, and the fraction can land
// an ulp past a half. With d from a fused operation there is nothing left to fuse, the fraction is within [-1/2, 1/2],
// and the image within half a period of r: in [0, p] for r = p / 2.

MD_SIMD_INLINE md_128 md_mm_deperiodize_ps(md_128 x, md_128 r, md_128 p) {
    md_128 d  = md_mm_sub_ps(x, r);
    md_128 dx = md_mm_div_ps(d, p);
    dx = md_mm_sub_ps(dx, md_mm_round_ps(dx));
    md_128 x_prim = md_mm_add_ps(r, md_mm_mul_ps(dx, p));
    return md_mm_select_ps(x_prim, x, md_mm_cmpeq_ps(p, md_mm_setzero_ps()));
}

MD_SIMD_INLINE md_128 md_mm_deperiodize2_ps(md_128 x, md_128 r, md_128 p, md_128 rp) {
	md_128 dx = md_mm_fmsub_ps(x, rp, md_mm_mul_ps(r, rp));
	dx = md_mm_sub_ps(dx, md_mm_round_ps(dx));
	md_128 x_prim = md_mm_fmadd_ps(dx, p, r);
	return md_mm_select_ps(x_prim, x, md_mm_cmpeq_ps(p, md_mm_setzero_ps()));
}

MD_SIMD_INLINE md_128d md_mm_deperiodize_pd(md_128d x, md_128d r, md_128d p) {
    md_128d d  = md_mm_sub_pd(x, r);
    md_128d dx = md_mm_div_pd(d, p);
    dx = md_mm_sub_pd(dx, md_mm_round_pd(dx));
    md_128d x_prim = md_mm_add_pd(r, md_mm_mul_pd(dx, p));
    return md_mm_select_pd(x_prim, x, md_mm_cmpeq_pd(p, md_mm_setzero_pd()));
}

MD_SIMD_INLINE md_128d md_mm_deperiodize2_pd(md_128d x, md_128d r, md_128d p, md_128d rp) {
	md_128d dx = md_mm_fmsub_pd(x, rp, md_mm_mul_pd(r, rp));
	dx = md_mm_sub_pd(dx, md_mm_round_pd(dx));
	md_128d x_prim = md_mm_fmadd_pd(dx, p, r);
	return md_mm_select_pd(x_prim, x, md_mm_cmpeq_pd(p, md_mm_setzero_pd()));
}

MD_SIMD_INLINE md_256 md_mm256_deperiodize_ps(md_256 x, md_256 r, md_256 p) {
    md_256 d  = md_mm256_sub_ps(x, r);
    md_256 dx = md_mm256_div_ps(d, p);
    dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
    md_256 x_prim = md_mm256_add_ps(r, md_mm256_mul_ps(dx, p));
    return md_mm256_select_ps(x_prim, x, md_mm256_cmpeq_ps(p, md_mm256_setzero_ps()));
}

MD_SIMD_INLINE md_256 md_mm256_deperiodize2_ps(md_256 x, md_256 r, md_256 p, md_256 rp) {
	md_256 dx = md_mm256_fmsub_ps(x, rp, md_mm256_mul_ps(r, rp));
	dx = md_mm256_sub_ps(dx, md_mm256_round_ps(dx));
	md_256 x_prim = md_mm256_fmadd_ps(dx, p, r);
	return md_mm256_select_ps(x_prim, x, md_mm256_cmpeq_ps(p, md_mm256_setzero_ps()));
}

MD_SIMD_INLINE md_256d md_mm256_deperiodize_pd(md_256d x, md_256d r, md_256d p) {
    md_256d d  = md_mm256_sub_pd(x, r);
    md_256d dx = md_mm256_div_pd(d, p);
    dx = md_mm256_sub_pd(dx, md_mm256_round_pd(dx));
    md_256d x_prim = md_mm256_add_pd(r, md_mm256_mul_pd(dx, p));
    return md_mm256_select_pd(x_prim, x, md_mm256_cmpeq_pd(p, md_mm256_setzero_pd()));
}

MD_SIMD_INLINE md_256d md_mm256_deperiodize2_pd(md_256d x, md_256d r, md_256d p, md_256d rp) {
	md_256d dx = md_mm256_fmsub_pd(x, rp, md_mm256_mul_pd(r, rp));
	dx = md_mm256_sub_pd(dx, md_mm256_round_pd(dx));
	md_256d x_prim = md_mm256_fmadd_pd(dx, p, r);
	return md_mm256_select_pd(x_prim, x, md_mm256_cmpeq_pd(p, md_mm256_setzero_pd()));
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
	md_512 dx = _mm512_fmsub_ps(x, rp, _mm512_mul_ps(r, rp));
	dx = _mm512_sub_ps(dx, _mm512_roundscale_ps(dx, _MM_FROUND_TO_NEAREST_INT));
	md_512 x_prim = _mm512_fmadd_ps(dx, p, r);
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
	md_512d dx = _mm512_fmsub_pd(x, rp, _mm512_mul_pd(r, rp));
	dx = _mm512_sub_pd(dx, _mm512_roundscale_pd(dx, _MM_FROUND_TO_NEAREST_INT));
	md_512d x_prim = _mm512_fmadd_pd(dx, p, r);
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

MD_SIMD_INLINE md_128d md_mm_lerp_pd(md_128d a, md_128d b, double t) {
    return md_mm_add_pd(md_mm_mul_pd(a, md_mm_set1_pd(1.0 - t)), md_mm_mul_pd(b, md_mm_set1_pd(t)));
}

MD_SIMD_INLINE md_256d md_mm256_lerp_pd(md_256d a, md_256d b, double t) {
    return md_mm256_add_pd(md_mm256_mul_pd(a, md_mm256_set1_pd(1.0 - t)), md_mm256_mul_pd(b, md_mm256_set1_pd(t)));
}

// Constructs a permutation vector for compressing a vector based on a mask. The mask is a bitfield where each set bit indicates that the corresponding element in the input vector should be kept.
// This is typically the mask returned from movemask instructions after a comparison, and the resulting permutation vector can be used with a permute instruction to compress the vector to only the elements where the mask is set.

MD_SIMD_INLINE md_128i md_mm_compression_mask_4x32(int mask) {
    uint32_t key = sse_compress_lut[mask];
    return simde_mm_cvtepu8_epi32(simde_mm_cvtsi32_si128(key));
}

MD_SIMD_INLINE md_256i md_mm256_compression_mask_8x32(int mask) {
    uint64_t key = avx_compress_lut[mask];
    return simde_mm256_cvtepu8_epi32(simde_mm_cvtsi64_si128((long long)key));
}

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

    /* select the correct result from the two polynoms: where poly_mask is set the sine is y2 and the cosine y */
    xmm1 = md_mm_select_ps(y, y2, poly_mask);
    xmm2 = md_mm_select_ps(y2, y, poly_mask);

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
    x = md_mm256_fmadd_ps(y, md_mm256_set1_ps(-3.77489497744594108e-8f), x);

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

    /* select the correct result from the two polynoms: where poly_mask is set the sine is y2 and the cosine y */
    xmm1 = md_mm256_select_ps(y, y2, poly_mask);
    xmm2 = md_mm256_select_ps(y2, y, poly_mask);

    /* update the sign */
    *s = md_mm256_xor_ps(xmm1, sign_bit_sin);
    *c = md_mm256_xor_ps(xmm2, sign_bit_cos);
}

// NATIVE WIDTH
//
// md_xv is the widest float vector the target executes natively: eight lanes with AVX, four
// otherwise. simde provides the eight lane type everywhere, but on NEON it is two four lane halves
// behind a layer the compiler does not see through, and a loop over it runs several times slower
// than the same loop over four lanes. Kernels that have no reason to prefer eight lanes use this.
#if defined(__AVX__)
typedef md_256 md_xv;
#define MD_XV_WIDTH 8
#define md_xv_loadu_ps              md_mm256_loadu_ps
#define md_xv_storeu_ps             md_mm256_storeu_ps
#define md_xv_set1_ps               md_mm256_set1_ps
#define md_xv_setzero_ps            md_mm256_setzero_ps
#define md_xv_add_ps                md_mm256_add_ps
#define md_xv_sub_ps                md_mm256_sub_ps
#define md_xv_mul_ps                md_mm256_mul_ps
#define md_xv_fmadd_ps              md_mm256_fmadd_ps
#define md_xv_min_ps                md_mm256_min_ps
#define md_xv_max_ps                md_mm256_max_ps
#define md_xv_round_ps              md_mm256_round_ps
#define md_xv_blendv_ps             md_mm256_blendv_ps
#define md_xv_cmpeq_ps              md_mm256_cmpeq_ps
#define md_xv_sincos_ps             md_mm256_sincos_ps
#define md_xv_cubic_spline_ps       md_mm256_cubic_spline_ps
#define md_xv_reduce_add_ps         md_mm256_reduce_add_ps
#define md_xv_load_xyz_packed_ps    md_mm256_load_xyz_packed_ps
#define md_xv_store_xyz_packed_ps   md_mm256_store_xyz_packed_ps
#define md_xv_load_ps               md_mm256_load_ps
#define md_xv_store_ps              md_mm256_store_ps
#define md_xv_div_ps                md_mm256_div_ps
#define md_xv_sqrt_ps               md_mm256_sqrt_ps
#define md_xv_rsqrt_ps              md_mm256_rsqrt_ps
#define md_xv_abs_ps                md_mm256_abs_ps
#define md_xv_neg_ps                md_mm256_neg_ps
#define md_xv_floor_ps              md_mm256_floor_ps
#define md_xv_ceil_ps               md_mm256_ceil_ps
#define md_xv_fract_ps              md_mm256_fract_ps
#define md_xv_fnmadd_ps             md_mm256_fnmadd_ps
#define md_xv_fmsub_ps              md_mm256_fmsub_ps
#define md_xv_and_ps                md_mm256_and_ps
#define md_xv_or_ps                 md_mm256_or_ps
#define md_xv_xor_ps                md_mm256_xor_ps
#define md_xv_andnot_ps             md_mm256_andnot_ps
#define md_xv_cmplt_ps              md_mm256_cmplt_ps
#define md_xv_cmple_ps              md_mm256_cmple_ps
#define md_xv_cmpgt_ps              md_mm256_cmpgt_ps
#define md_xv_cmpge_ps              md_mm256_cmpge_ps
#define md_xv_cmpneq_ps             md_mm256_cmpneq_ps
#define md_xv_select_ps             md_mm256_select_ps
#define md_xv_movemask_ps           md_mm256_movemask_ps
#define md_xv_reduce_min_ps         md_mm256_reduce_min_ps
#define md_xv_reduce_max_ps         md_mm256_reduce_max_ps
#define md_xv_exp_ps                md_mm256_exp_ps
#define md_xv_minimage_ps           md_mm256_minimage_ps
#else
typedef md_128 md_xv;
#define MD_XV_WIDTH 4
#define md_xv_loadu_ps              md_mm_loadu_ps
#define md_xv_storeu_ps             md_mm_storeu_ps
#define md_xv_set1_ps               md_mm_set1_ps
#define md_xv_setzero_ps            md_mm_setzero_ps
#define md_xv_add_ps                md_mm_add_ps
#define md_xv_sub_ps                md_mm_sub_ps
#define md_xv_mul_ps                md_mm_mul_ps
#define md_xv_fmadd_ps              md_mm_fmadd_ps
#define md_xv_min_ps                md_mm_min_ps
#define md_xv_max_ps                md_mm_max_ps
#define md_xv_round_ps              md_mm_round_ps
#define md_xv_blendv_ps             md_mm_blendv_ps
#define md_xv_cmpeq_ps              md_mm_cmpeq_ps
#define md_xv_sincos_ps             md_mm_sincos_ps
#define md_xv_cubic_spline_ps       md_mm_cubic_spline_ps
#define md_xv_reduce_add_ps         md_mm_reduce_add_ps
#define md_xv_load_xyz_packed_ps    md_mm_load_xyz_packed_ps
#define md_xv_store_xyz_packed_ps   md_mm_store_xyz_packed_ps
#define md_xv_load_ps               md_mm_load_ps
#define md_xv_store_ps              md_mm_store_ps
#define md_xv_div_ps                md_mm_div_ps
#define md_xv_sqrt_ps               md_mm_sqrt_ps
#define md_xv_rsqrt_ps              md_mm_rsqrt_ps
#define md_xv_abs_ps                md_mm_abs_ps
#define md_xv_neg_ps                md_mm_neg_ps
#define md_xv_floor_ps              md_mm_floor_ps
#define md_xv_ceil_ps               md_mm_ceil_ps
#define md_xv_fract_ps              md_mm_fract_ps
#define md_xv_fnmadd_ps             md_mm_fnmadd_ps
#define md_xv_fmsub_ps              md_mm_fmsub_ps
#define md_xv_and_ps                md_mm_and_ps
#define md_xv_or_ps                 md_mm_or_ps
#define md_xv_xor_ps                md_mm_xor_ps
#define md_xv_andnot_ps             md_mm_andnot_ps
#define md_xv_cmplt_ps              md_mm_cmplt_ps
#define md_xv_cmple_ps              md_mm_cmple_ps
#define md_xv_cmpgt_ps              md_mm_cmpgt_ps
#define md_xv_cmpge_ps              md_mm_cmpge_ps
#define md_xv_cmpneq_ps             md_mm_cmpneq_ps
#define md_xv_select_ps             md_mm_select_ps
#define md_xv_movemask_ps           md_mm_movemask_ps
#define md_xv_reduce_min_ps         md_mm_reduce_min_ps
#define md_xv_reduce_max_ps         md_mm_reduce_max_ps
#define md_xv_exp_ps                md_mm_exp_ps
#define md_xv_minimage_ps           md_mm_minimage_ps
#endif

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
    
    const __mmask16 poly_mask = _mm512_cmpeq_epi32_mask(imm2, _mm512_setzero_si512());

    const __m512 swap_sign_bit_sin = _mm512_castsi512_ps(imm0);

    imm3 = _mm512_sub_epi32(imm3, _mm512_set1_epi32(2));
    imm3 = _mm512_andnot_epi32(imm3, _mm512_set1_epi32(4));
    imm3 = _mm512_slli_epi32(imm3, 29);

    const __m512 sign_bit_cos = _mm512_castsi512_ps(imm3);
    sign_bit_sin = _mm512_xor_ps(sign_bit_sin, swap_sign_bit_sin);

    /* The magic pass: "Extended precision modular arithmetic" 
    x = ((x - y * DP1) - y * DP2) - y * DP3; */
    x = _mm512_fmadd_ps(y, _mm512_set1_ps(-0.78515625f), x);
    x = _mm512_fmadd_ps(y, _mm512_set1_ps(-2.4187564849853515625E-4f), x);
    x = _mm512_fmadd_ps(y, _mm512_set1_ps(-3.77489497744594108e-8f), x);

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

    /* select the correct result from the two polynoms: where poly_mask is set the sine is y2 and the cosine y */
    xmm1 = _mm512_mask_blend_ps(poly_mask, y, y2);
    xmm2 = _mm512_mask_blend_ps(poly_mask, y2, y);

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

    fx = md_mm_floor_ps(fx);

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

    fx = md_mm256_floor_ps(fx);

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

    fx = _mm512_roundscale_ps(fx, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

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
