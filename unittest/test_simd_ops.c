#include "utest.h"

#include <core/md_common.h>
#include <core/md_simd.h>

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

// THE OPERATIONS OF md_simd, one by one, against a scalar reference, at every width they come in: 4 and 8 lanes of
// float, 2 and 4 of double, 4 and 8 of int32, 2 and 4 of int64, and the native width (md_xv). Every name md_simd
// defines is used here, so a name that maps to the wrong simde function (or to none) fails to compile or fails a test
// rather than waiting for its first user. md_mm512_* is tested where AVX-512 is enabled.
//
// The inputs are random (a fixed sequence) with the special values tested on their own: NaN semantics differ between
// targets by design (see md_simd.h) and the scalar references are compiled with -ffast-math on some targets.

typedef union { md_128  v; float   f[4]; uint32_t u[4]; } u128;
typedef union { md_256  v; float   f[8]; uint32_t u[8]; } u256;
typedef union { md_128d v; double  d[2]; uint64_t u[2]; } u128d;
typedef union { md_256d v; double  d[4]; uint64_t u[4]; } u256d;
typedef union { md_128i v; int32_t i[4]; int64_t l[2]; int16_t s[8];  int8_t b[16]; uint32_t u[4]; } u128i;
typedef union { md_256i v; int32_t i[8]; int64_t l[4]; int16_t s[16]; int8_t b[32]; uint32_t u[8]; } u256i;

#define NUM_CASES 256

static uint64_t rng;
static uint32_t rnd_u32(void) {
    rng = rng * 6364136223846793005ull + 1442695040888963407ull;
    return (uint32_t)(rng >> 32);
}
static uint64_t rnd_u64(void) { return ((uint64_t)rnd_u32() << 32) | rnd_u32(); }

// Never zero, |x| in [2^-8, 2^9), every mantissa bit random
static float rnd_f32(void) {
    const uint32_t r = rnd_u32();
    const uint32_t bits = (r & 0x80000000u) | ((uint32_t)(127 - 8 + (int)((r >> 23) % 17)) << 23) | (rnd_u32() & 0x7FFFFFu);
    float f;
    memcpy(&f, &bits, 4);
    return f;
}
static double rnd_f64(void) {
    const uint64_t r = rnd_u64();
    const uint64_t bits = (r & 0x8000000000000000ull) | ((uint64_t)(1023 - 8 + (int)((r >> 52) % 17)) << 52) | (rnd_u64() & 0xFFFFFFFFFFFFFull);
    double d;
    memcpy(&d, &bits, 8);
    return d;
}
static void fill_f32(float* x, int n) { for (int i = 0; i < n; ++i) x[i] = rnd_f32(); }
static void fill_f64(double* x, int n) { for (int i = 0; i < n; ++i) x[i] = rnd_f64(); }

static uint32_t f32_bits(float f)  { uint32_t u; memcpy(&u, &f, 4); return u; }
static uint64_t f64_bits(double d) { uint64_t u; memcpy(&u, &d, 8); return u; }
static bool f32_nan(float f)  { return (f32_bits(f) & 0x7FFFFFFFu) > 0x7F800000u; }
static bool f64_nan(double d) { return (f64_bits(d) & 0x7FFFFFFFFFFFFFFFull) > 0x7FF0000000000000ull; }
// Bit for bit, any NaN being the same as any other
static bool f32_same(float a, float b)   { return f32_bits(a) == f32_bits(b) || (f32_nan(a) && f32_nan(b)); }
static bool f64_same(double a, double b) { return f64_bits(a) == f64_bits(b) || (f64_nan(a) && f64_nan(b)); }
static float  f32_from_bits(uint32_t u) { float f; memcpy(&f, &u, 4); return f; }
static uint32_t f32_ulps(float a, float b) {
    const int32_t x = (int32_t)f32_bits(a), y = (int32_t)f32_bits(b);
    const int64_t ox = x < 0 ? (int64_t)INT32_MIN - x : x, oy = y < 0 ? (int64_t)INT32_MIN - y : y;   // ordered
    return (uint32_t)(ox > oy ? ox - oy : oy - ox);
}

// GCC with -ffast-math computes a vector division of floats as an approximate reciprocal and a Newton-Raphson step:
// within 2 ulp, but a / a is not 1. mdlib does not build with it; division is exact unless a build adds it.
#if defined(__FAST_MATH__)
#define DIV_ULPS 2
#else
#define DIV_ULPS 0
#endif
static double f64_from_bits(uint64_t u) { double d; memcpy(&d, &u, 8); return d; }

#define FAIL_F32(name, lane, e, a) printf("%s, lane %d: expected %.9g (0x%08x), got %.9g (0x%08x)\n", name, (int)(lane), (double)(e), f32_bits(e), (double)(a), f32_bits(a))
#define FAIL_F64(name, lane, e, a) printf("%s, lane %d: expected %.17g, got %.17g\n", name, (int)(lane), (e), (a))

#define CHECK_F32(name, lane, expected, actual) do { const float  e_ = (expected), a_ = (actual); if (!f32_same(e_, a_)) { FAIL_F32(name, lane, e_, a_); fails += 1; } } while (0)
#define CHECK_F64(name, lane, expected, actual) do { const double e_ = (expected), a_ = (actual); if (!f64_same(e_, a_)) { FAIL_F64(name, lane, e_, a_); fails += 1; } } while (0)
#define CHECK_NEAR_F32(name, lane, expected, actual, rel) do { const double e_ = (expected), a_ = (actual); if (!(fabs(e_ - a_) <= (rel) * (fabs(e_) + 1e-30))) { FAIL_F32(name, lane, (float)e_, (float)a_); fails += 1; } } while (0)
#define CHECK_NEAR_F64(name, lane, expected, actual, rel) do { const double e_ = (expected), a_ = (actual); if (!(fabs(e_ - a_) <= (rel) * (fabs(e_) + 1e-300))) { FAIL_F64(name, lane, e_, a_); fails += 1; } } while (0)
#define CHECK_ULPS_F32(name, lane, expected, actual, ulps) do { const float e_ = (expected), a_ = (actual); if (f32_ulps(e_, a_) > (ulps)) { FAIL_F32(name, lane, e_, a_); fails += 1; } } while (0)
#define CHECK_I64(name, lane, expected, actual) do { const int64_t e_ = (int64_t)(expected), a_ = (int64_t)(actual); if (e_ != a_) { printf("%s, lane %d: expected %lld (0x%llx), got %lld (0x%llx)\n", name, (int)(lane), (long long)e_, (unsigned long long)e_, (long long)a_, (unsigned long long)a_); fails += 1; } } while (0)

// Float references: one rounding of the exact result. A sum, difference, product, quotient or square root of two
// floats is exact in double up to the rounding to float, so these hold whatever the flags of this file.
static float ref_add(float a, float b) { return (float)((double)a + (double)b); }
static float ref_sub(float a, float b) { return (float)((double)a - (double)b); }
static float ref_mul(float a, float b) { return (float)((double)a * (double)b); }
static float ref_div(float a, float b) { return (float)((double)a / (double)b); }
static float ref_min(float a, float b) { return a < b ? a : b; }
static float ref_max(float a, float b) { return a > b ? a : b; }
static float ref_and(float a, float b)    { return f32_from_bits(f32_bits(a) & f32_bits(b)); }
static float ref_or(float a, float b)     { return f32_from_bits(f32_bits(a) | f32_bits(b)); }
static float ref_xor(float a, float b)    { return f32_from_bits(f32_bits(a) ^ f32_bits(b)); }
static float ref_andnot(float a, float b) { return f32_from_bits(~f32_bits(a) & f32_bits(b)); }
static float ref_cmpeq(float a, float b)  { return f32_from_bits(a == b ? ~0u : 0u); }
static float ref_cmpneq(float a, float b) { return f32_from_bits(a != b ? ~0u : 0u); }
static float ref_cmplt(float a, float b)  { return f32_from_bits(a <  b ? ~0u : 0u); }
static float ref_cmple(float a, float b)  { return f32_from_bits(a <= b ? ~0u : 0u); }
static float ref_cmpgt(float a, float b)  { return f32_from_bits(a >  b ? ~0u : 0u); }
static float ref_cmpge(float a, float b)  { return f32_from_bits(a >= b ? ~0u : 0u); }
static float ref_step(float edge, float x) { return x >= edge ? 1.0f : 0.0f; }
static float ref_sqrt(float a)  { return (float)sqrt(fabs((double)a)); }
static float ref_abs(float a)   { return f32_from_bits(f32_bits(a) & 0x7FFFFFFFu); }
static float ref_neg(float a)   { return f32_from_bits(f32_bits(a) ^ 0x80000000u); }
static float ref_sign(float a)  { return (f32_bits(a) & 0x80000000u) ? -1.0f : 1.0f; }
static float ref_floor(float a) { return (float)floor((double)a); }
static float ref_ceil(float a)  { return (float)ceil((double)a); }
static float ref_round(float a) { return (float)rint((double)a); }
static float ref_fract(float a) { return (float)((double)a - floor((double)a)); }
static float ref_fmadd(float a, float b, float c)  { return fmaf(a, b, c); }
static float ref_fnmadd(float a, float b, float c) { return fmaf(-a, b, c); }
static float ref_fmsub(float a, float b, float c)  { return fmaf(a, b, -c); }
static float ref_fnmsub(float a, float b, float c) { return fmaf(-a, b, -c); }

static double ref_add_d(double a, double b) { return a + b; }
static double ref_sub_d(double a, double b) { return a - b; }
static double ref_mul_d(double a, double b) { return a * b; }
static double ref_div_d(double a, double b) { return a / b; }
static double ref_min_d(double a, double b) { return a < b ? a : b; }
static double ref_max_d(double a, double b) { return a > b ? a : b; }
static double ref_and_d(double a, double b)    { return f64_from_bits(f64_bits(a) & f64_bits(b)); }
static double ref_or_d(double a, double b)     { return f64_from_bits(f64_bits(a) | f64_bits(b)); }
static double ref_xor_d(double a, double b)    { return f64_from_bits(f64_bits(a) ^ f64_bits(b)); }
static double ref_andnot_d(double a, double b) { return f64_from_bits(~f64_bits(a) & f64_bits(b)); }
static double ref_cmpeq_d(double a, double b)  { return f64_from_bits(a == b ? ~0ull : 0ull); }
static double ref_cmpneq_d(double a, double b) { return f64_from_bits(a != b ? ~0ull : 0ull); }
static double ref_cmplt_d(double a, double b)  { return f64_from_bits(a <  b ? ~0ull : 0ull); }
static double ref_cmple_d(double a, double b)  { return f64_from_bits(a <= b ? ~0ull : 0ull); }
static double ref_cmpgt_d(double a, double b)  { return f64_from_bits(a >  b ? ~0ull : 0ull); }
static double ref_cmpge_d(double a, double b)  { return f64_from_bits(a >= b ? ~0ull : 0ull); }
static double ref_step_d(double edge, double x) { return x >= edge ? 1.0 : 0.0; }
static double ref_sqrt_d(double a)  { return sqrt(fabs(a)); }
static double ref_abs_d(double a)   { return f64_from_bits(f64_bits(a) & 0x7FFFFFFFFFFFFFFFull); }
static double ref_neg_d(double a)   { return f64_from_bits(f64_bits(a) ^ 0x8000000000000000ull); }
static double ref_sign_d(double a)  { return (f64_bits(a) & 0x8000000000000000ull) ? -1.0 : 1.0; }
static double ref_floor_d(double a) { return floor(a); }
static double ref_ceil_d(double a)  { return ceil(a); }
static double ref_round_d(double a) { return rint(a); }
static double ref_fract_d(double a) { return a - floor(a); }
static double ref_fmadd_d(double a, double b, double c)  { return fma(a, b, c); }
static double ref_fnmadd_d(double a, double b, double c) { return fma(-a, b, c); }
static double ref_fmsub_d(double a, double b, double c)  { return fma(a, b, -c); }
static double ref_fnmsub_d(double a, double b, double c) { return fma(-a, b, -c); }

// A float operation at four and at eight lanes
#define TEST_PS1(OP4, OP8, REF, ARG) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        float a[8]; fill_f32(a, 8); \
        for (int i = 0; i < 8; ++i) a[i] = ARG(a[i]); \
        u128 r4 = { .v = OP4(md_mm_loadu_ps(a)) }; \
        u256 r8 = { .v = OP8(md_mm256_loadu_ps(a)) }; \
        for (int i = 0; i < 8; ++i) { \
            if (i < 4) CHECK_F32(#OP4, i, REF(a[i]), r4.f[i]); \
            CHECK_F32(#OP8, i, REF(a[i]), r8.f[i]); \
        } \
    }

#define TEST_PS2(OP4, OP8, REF) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        float a[8], b[8]; fill_f32(a, 8); fill_f32(b, 8); \
        if (c & 1) for (int i = 0; i < 8; i += 2) b[i] = a[i];  /* equal lanes for the comparisons */ \
        u128 r4 = { .v = OP4(md_mm_loadu_ps(a), md_mm_loadu_ps(b)) }; \
        u256 r8 = { .v = OP8(md_mm256_loadu_ps(a), md_mm256_loadu_ps(b)) }; \
        for (int i = 0; i < 8; ++i) { \
            if (i < 4) CHECK_F32(#OP4, i, REF(a[i], b[i]), r4.f[i]); \
            CHECK_F32(#OP8, i, REF(a[i], b[i]), r8.f[i]); \
        } \
    }

// Within a number of ulp of the reference
#define TEST_PS2_ULPS(OP4, OP8, REF, ULPS) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        float a[8], b[8]; fill_f32(a, 8); fill_f32(b, 8); \
        if (c & 1) for (int i = 0; i < 8; i += 2) b[i] = a[i]; \
        u128 r4 = { .v = OP4(md_mm_loadu_ps(a), md_mm_loadu_ps(b)) }; \
        u256 r8 = { .v = OP8(md_mm256_loadu_ps(a), md_mm256_loadu_ps(b)) }; \
        for (int i = 0; i < 8; ++i) { \
            if (i < 4) CHECK_ULPS_F32(#OP4, i, REF(a[i], b[i]), r4.f[i], ULPS); \
            CHECK_ULPS_F32(#OP8, i, REF(a[i], b[i]), r8.f[i], ULPS); \
        } \
    }

#define TEST_PS3(OP4, OP8, REF) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        float a[8], b[8], d[8]; fill_f32(a, 8); fill_f32(b, 8); fill_f32(d, 8); \
        u128 r4 = { .v = OP4(md_mm_loadu_ps(a), md_mm_loadu_ps(b), md_mm_loadu_ps(d)) }; \
        u256 r8 = { .v = OP8(md_mm256_loadu_ps(a), md_mm256_loadu_ps(b), md_mm256_loadu_ps(d)) }; \
        for (int i = 0; i < 8; ++i) { \
            if (i < 4) CHECK_F32(#OP4, i, REF(a[i], b[i], d[i]), r4.f[i]); \
            CHECK_F32(#OP8, i, REF(a[i], b[i], d[i]), r8.f[i]); \
        } \
    }

#define TEST_PD1(OP2, OP4, REF, ARG) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        double a[4]; fill_f64(a, 4); \
        for (int i = 0; i < 4; ++i) a[i] = ARG(a[i]); \
        u128d r2 = { .v = OP2(md_mm_loadu_pd(a)) }; \
        u256d r4 = { .v = OP4(md_mm256_loadu_pd(a)) }; \
        for (int i = 0; i < 4; ++i) { \
            if (i < 2) CHECK_F64(#OP2, i, REF(a[i]), r2.d[i]); \
            CHECK_F64(#OP4, i, REF(a[i]), r4.d[i]); \
        } \
    }

#define TEST_PD2(OP2, OP4, REF) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        double a[4], b[4]; fill_f64(a, 4); fill_f64(b, 4); \
        if (c & 1) for (int i = 0; i < 4; i += 2) b[i] = a[i]; \
        u128d r2 = { .v = OP2(md_mm_loadu_pd(a), md_mm_loadu_pd(b)) }; \
        u256d r4 = { .v = OP4(md_mm256_loadu_pd(a), md_mm256_loadu_pd(b)) }; \
        for (int i = 0; i < 4; ++i) { \
            if (i < 2) CHECK_F64(#OP2, i, REF(a[i], b[i]), r2.d[i]); \
            CHECK_F64(#OP4, i, REF(a[i], b[i]), r4.d[i]); \
        } \
    }

#define TEST_PD3(OP2, OP4, REF) \
    for (int c = 0; c < NUM_CASES; ++c) { \
        double a[4], b[4], d[4]; fill_f64(a, 4); fill_f64(b, 4); fill_f64(d, 4); \
        u128d r2 = { .v = OP2(md_mm_loadu_pd(a), md_mm_loadu_pd(b), md_mm_loadu_pd(d)) }; \
        u256d r4 = { .v = OP4(md_mm256_loadu_pd(a), md_mm256_loadu_pd(b), md_mm256_loadu_pd(d)) }; \
        for (int i = 0; i < 4; ++i) { \
            if (i < 2) CHECK_F64(#OP2, i, REF(a[i], b[i], d[i]), r2.d[i]); \
            CHECK_F64(#OP4, i, REF(a[i], b[i], d[i]), r4.d[i]); \
        } \
    }

#define AS_IS(x) (x)
#define POSITIVE(x) ((float)fabs(x))
#define POSITIVE_D(x) fabs(x)
// Up to 2^20 with fractions, where floor, ceil and round have something to do
#define SPREAD(x) ((x) * 4096.0f)
#define SPREAD_D(x) ((x) * 4096.0)
// Halves, for the ties of round
#define HALVES(x) (floorf((x) * 8.0f) * 0.5f)
#define HALVES_D(x) (floor((x) * 8.0) * 0.5)

// ARITHMETIC

UTEST(simd_ops, arithmetic_ps) {
    rng = 1;
    int fails = 0;
    TEST_PS2(md_mm_add_ps, md_mm256_add_ps, ref_add);
    TEST_PS2(md_mm_sub_ps, md_mm256_sub_ps, ref_sub);
    TEST_PS2(md_mm_mul_ps, md_mm256_mul_ps, ref_mul);
    TEST_PS2_ULPS(md_mm_div_ps, md_mm256_div_ps, ref_div, DIV_ULPS);
    TEST_PS2(md_mm_min_ps, md_mm256_min_ps, ref_min);
    TEST_PS2(md_mm_max_ps, md_mm256_max_ps, ref_max);
    TEST_PS3(md_mm_fmadd_ps,  md_mm256_fmadd_ps,  ref_fmadd);
    TEST_PS3(md_mm_fnmadd_ps, md_mm256_fnmadd_ps, ref_fnmadd);
    TEST_PS3(md_mm_fmsub_ps,  md_mm256_fmsub_ps,  ref_fmsub);
    TEST_PS3(md_mm_fnmsub_ps, md_mm256_fnmsub_ps, ref_fnmsub);
    TEST_PS1(md_mm_sqrt_ps, md_mm256_sqrt_ps, ref_sqrt, POSITIVE);
    TEST_PS1(md_mm_abs_ps,  md_mm256_abs_ps,  ref_abs,  AS_IS);
    TEST_PS1(md_mm_neg_ps,  md_mm256_neg_ps,  ref_neg,  AS_IS);
    TEST_PS1(md_mm_sign_ps, md_mm256_sign_ps, ref_sign, AS_IS);
    TEST_PS1(md_mm_floor_ps, md_mm256_floor_ps, ref_floor, SPREAD);
    TEST_PS1(md_mm_ceil_ps,  md_mm256_ceil_ps,  ref_ceil,  SPREAD);
    TEST_PS1(md_mm_round_ps, md_mm256_round_ps, ref_round, SPREAD);
    TEST_PS1(md_mm_round_ps, md_mm256_round_ps, ref_round, HALVES);
    TEST_PS1(md_mm_fract_ps, md_mm256_fract_ps, ref_fract, SPREAD);
    TEST_PS2(md_mm_step_ps, md_mm256_step_ps, ref_step);
    EXPECT_EQ(0, fails);

    // The single lane forms
    u128 s = { .v = md_mm_add_ss(md_mm_set_ps(4, 3, 2, 1), md_mm_set_ps(40, 30, 20, 10)) };
    EXPECT_EQ(11.0f, s.f[0]);
    EXPECT_EQ(2.0f, s.f[1]);
    u128d sd = { .v = md_mm_add_sd(md_mm_set_pd(2, 1), md_mm_set_pd(20, 10)) };
    EXPECT_EQ(11.0, sd.d[0]);
    EXPECT_EQ(2.0, sd.d[1]);
}

UTEST(simd_ops, arithmetic_pd) {
    rng = 2;
    int fails = 0;
    TEST_PD2(md_mm_add_pd, md_mm256_add_pd, ref_add_d);
    TEST_PD2(md_mm_sub_pd, md_mm256_sub_pd, ref_sub_d);
    TEST_PD2(md_mm_mul_pd, md_mm256_mul_pd, ref_mul_d);
    TEST_PD2(md_mm_div_pd, md_mm256_div_pd, ref_div_d);
    TEST_PD2(md_mm_min_pd, md_mm256_min_pd, ref_min_d);
    TEST_PD2(md_mm_max_pd, md_mm256_max_pd, ref_max_d);
    TEST_PD3(md_mm_fmadd_pd,  md_mm256_fmadd_pd,  ref_fmadd_d);
    TEST_PD3(md_mm_fnmadd_pd, md_mm256_fnmadd_pd, ref_fnmadd_d);
    TEST_PD3(md_mm_fmsub_pd,  md_mm256_fmsub_pd,  ref_fmsub_d);
    TEST_PD3(md_mm_fnmsub_pd, md_mm256_fnmsub_pd, ref_fnmsub_d);
    TEST_PD1(md_mm_sqrt_pd, md_mm256_sqrt_pd, ref_sqrt_d, POSITIVE_D);
    TEST_PD1(md_mm_abs_pd,  md_mm256_abs_pd,  ref_abs_d,  AS_IS);
    TEST_PD1(md_mm_neg_pd,  md_mm256_neg_pd,  ref_neg_d,  AS_IS);
    TEST_PD1(md_mm_sign_pd, md_mm256_sign_pd, ref_sign_d, AS_IS);
    TEST_PD1(md_mm_floor_pd, md_mm256_floor_pd, ref_floor_d, SPREAD_D);
    TEST_PD1(md_mm_ceil_pd,  md_mm256_ceil_pd,  ref_ceil_d,  SPREAD_D);
    TEST_PD1(md_mm_round_pd, md_mm256_round_pd, ref_round_d, SPREAD_D);
    TEST_PD1(md_mm_round_pd, md_mm256_round_pd, ref_round_d, HALVES_D);
    TEST_PD1(md_mm_fract_pd, md_mm256_fract_pd, ref_fract_d, SPREAD_D);
    TEST_PD2(md_mm_step_pd, md_mm256_step_pd, ref_step_d);
    EXPECT_EQ(0, fails);
}

// The approximate reciprocal square root: within 1.5 * 2^-12 of the exact value, the x86 bound, on every target
UTEST(simd_ops, rsqrt) {
    rng = 3;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        float a[8];
        fill_f32(a, 8);
        for (int i = 0; i < 8; ++i) a[i] = fabsf(a[i]);
        u128 r4 = { .v = md_mm_rsqrt_ps(md_mm_loadu_ps(a)) };
        u256 r8 = { .v = md_mm256_rsqrt_ps(md_mm256_loadu_ps(a)) };
        for (int i = 0; i < 8; ++i) {
            const double e = 1.0 / sqrt((double)a[i]);
            if (i < 4) CHECK_NEAR_F32("md_mm_rsqrt_ps", i, e, r4.f[i], 1.5 / 4096.0);
            CHECK_NEAR_F32("md_mm256_rsqrt_ps", i, e, r8.f[i], 1.5 / 4096.0);
        }
    }
    EXPECT_EQ(0, fails);
}

// LOGIC AND COMPARISON

UTEST(simd_ops, logic_and_compare_ps) {
    rng = 4;
    int fails = 0;
    TEST_PS2(md_mm_and_ps,    md_mm256_and_ps,    ref_and);
    TEST_PS2(md_mm_or_ps,     md_mm256_or_ps,     ref_or);
    TEST_PS2(md_mm_xor_ps,    md_mm256_xor_ps,    ref_xor);
    TEST_PS2(md_mm_andnot_ps, md_mm256_andnot_ps, ref_andnot);
    TEST_PS2(md_mm_cmpeq_ps,  md_mm256_cmpeq_ps,  ref_cmpeq);
    TEST_PS2(md_mm_cmpneq_ps, md_mm256_cmpneq_ps, ref_cmpneq);
    TEST_PS2(md_mm_cmplt_ps,  md_mm256_cmplt_ps,  ref_cmplt);
    TEST_PS2(md_mm_cmple_ps,  md_mm256_cmple_ps,  ref_cmple);
    TEST_PS2(md_mm_cmpgt_ps,  md_mm256_cmpgt_ps,  ref_cmpgt);
    TEST_PS2(md_mm_cmpge_ps,  md_mm256_cmpge_ps,  ref_cmpge);
#define CMP_LT_PS(a, b)    md_mm_cmp_ps(a, b, SIMDE_CMP_LT_OQ)
#define CMP_LT_PS256(a, b) md_mm256_cmp_ps(a, b, SIMDE_CMP_LT_OQ)
#define CMP_GE_PS(a, b)    md_mm_cmp_ps(a, b, SIMDE_CMP_GE_OQ)
#define CMP_GE_PS256(a, b) md_mm256_cmp_ps(a, b, SIMDE_CMP_GE_OQ)
    TEST_PS2(CMP_LT_PS, CMP_LT_PS256, ref_cmplt);
    TEST_PS2(CMP_GE_PS, CMP_GE_PS256, ref_cmpge);
    EXPECT_EQ(0, fails);
}

UTEST(simd_ops, logic_and_compare_pd) {
    rng = 5;
    int fails = 0;
    TEST_PD2(md_mm_and_pd,    md_mm256_and_pd,    ref_and_d);
    TEST_PD2(md_mm_or_pd,     md_mm256_or_pd,     ref_or_d);
    TEST_PD2(md_mm_xor_pd,    md_mm256_xor_pd,    ref_xor_d);
    TEST_PD2(md_mm_andnot_pd, md_mm256_andnot_pd, ref_andnot_d);
    TEST_PD2(md_mm_cmpeq_pd,  md_mm256_cmpeq_pd,  ref_cmpeq_d);
    TEST_PD2(md_mm_cmpneq_pd, md_mm256_cmpneq_pd, ref_cmpneq_d);
    TEST_PD2(md_mm_cmplt_pd,  md_mm256_cmplt_pd,  ref_cmplt_d);
    TEST_PD2(md_mm_cmple_pd,  md_mm256_cmple_pd,  ref_cmple_d);
    TEST_PD2(md_mm_cmpgt_pd,  md_mm256_cmpgt_pd,  ref_cmpgt_d);
    TEST_PD2(md_mm_cmpge_pd,  md_mm256_cmpge_pd,  ref_cmpge_d);
#define CMP_LE_PD(a, b)    md_mm_cmp_pd(a, b, SIMDE_CMP_LE_OQ)
#define CMP_LE_PD256(a, b) md_mm256_cmp_pd(a, b, SIMDE_CMP_LE_OQ)
    TEST_PD2(CMP_LE_PD, CMP_LE_PD256, ref_cmple_d);
    EXPECT_EQ(0, fails);
}

// NaN: what is the same on every target. A NaN in the first operand of min and max gives the second (what
// min(x, bound) clamps by); comparisons with a NaN are false, except not equal. Skipped where the build assumes there
// are no NaN (-ffast-math), as the compiler may then fold these any way it likes.
#if !defined(__FINITE_MATH_ONLY__) || __FINITE_MATH_ONLY__ == 0
UTEST(simd_ops, nan) {
    const float qnan = f32_from_bits(0x7FC00000u);
    const md_128 n4 = md_mm_set1_ps(qnan), x4 = md_mm_set1_ps(2.0f);
    const md_256 n8 = md_mm256_set1_ps(qnan), x8 = md_mm256_set1_ps(2.0f);
    u128 r;
    u256 r8;
    r.v = md_mm_min_ps(n4, x4);       EXPECT_EQ(2.0f, r.f[0]);
    r.v = md_mm_max_ps(n4, x4);       EXPECT_EQ(2.0f, r.f[0]);
    r8.v = md_mm256_min_ps(n8, x8);   EXPECT_EQ(2.0f, r8.f[7]);
    r8.v = md_mm256_max_ps(n8, x8);   EXPECT_EQ(2.0f, r8.f[7]);
    r.v = md_mm_cmplt_ps(n4, x4);     EXPECT_EQ(0u, r.u[0]);
    r.v = md_mm_cmpge_ps(n4, x4);     EXPECT_EQ(0u, r.u[0]);
    r.v = md_mm_cmpeq_ps(n4, n4);     EXPECT_EQ(0u, r.u[0]);
    r.v = md_mm_cmpneq_ps(n4, x4);    EXPECT_EQ(~0u, r.u[0]);
    r8.v = md_mm256_cmpneq_ps(n8, x8); EXPECT_EQ(~0u, r8.u[5]);
    const double dnan = f64_from_bits(0x7FF8000000000000ull);
    u128d d;
    d.v = md_mm_min_pd(md_mm_set1_pd(dnan), md_mm_set1_pd(2.0)); EXPECT_EQ(2.0, d.d[1]);
    d.v = md_mm_max_pd(md_mm_set1_pd(dnan), md_mm_set1_pd(2.0)); EXPECT_EQ(2.0, d.d[1]);
}
#endif

// BLENDS

UTEST(simd_ops, blend_and_select) {
    rng = 6;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        float a[8], b[8], m[8];
        double ad[4], bd[4], md[4];
        fill_f32(a, 8); fill_f32(b, 8); fill_f32(m, 8);
        fill_f64(ad, 4); fill_f64(bd, 4); fill_f64(md, 4);
        // blendv takes the sign bit of any value; select a comparison mask
        u128 r4 = { .v = md_mm_blendv_ps(md_mm_loadu_ps(a), md_mm_loadu_ps(b), md_mm_loadu_ps(m)) };
        u256 r8 = { .v = md_mm256_blendv_ps(md_mm256_loadu_ps(a), md_mm256_loadu_ps(b), md_mm256_loadu_ps(m)) };
        u128 s4 = { .v = md_mm_select_ps(md_mm_loadu_ps(a), md_mm_loadu_ps(b), md_mm_cmplt_ps(md_mm_loadu_ps(m), md_mm_setzero_ps())) };
        u256 s8 = { .v = md_mm256_select_ps(md_mm256_loadu_ps(a), md_mm256_loadu_ps(b), md_mm256_cmplt_ps(md_mm256_loadu_ps(m), md_mm256_setzero_ps())) };
        for (int i = 0; i < 8; ++i) {
            const float e = m[i] < 0 ? b[i] : a[i];
            if (i < 4) { CHECK_F32("md_mm_blendv_ps", i, e, r4.f[i]); CHECK_F32("md_mm_select_ps", i, e, s4.f[i]); }
            CHECK_F32("md_mm256_blendv_ps", i, e, r8.f[i]);
            CHECK_F32("md_mm256_select_ps", i, e, s8.f[i]);
        }
        u128d r2 = { .v = md_mm_blendv_pd(md_mm_loadu_pd(ad), md_mm_loadu_pd(bd), md_mm_loadu_pd(md)) };
        u256d r4d = { .v = md_mm256_blendv_pd(md_mm256_loadu_pd(ad), md_mm256_loadu_pd(bd), md_mm256_loadu_pd(md)) };
        u128d s2 = { .v = md_mm_select_pd(md_mm_loadu_pd(ad), md_mm_loadu_pd(bd), md_mm_cmplt_pd(md_mm_loadu_pd(md), md_mm_setzero_pd())) };
        u256d s4d = { .v = md_mm256_select_pd(md_mm256_loadu_pd(ad), md_mm256_loadu_pd(bd), md_mm256_cmplt_pd(md_mm256_loadu_pd(md), md_mm256_setzero_pd())) };
        for (int i = 0; i < 4; ++i) {
            const double e = md[i] < 0 ? bd[i] : ad[i];
            if (i < 2) { CHECK_F64("md_mm_blendv_pd", i, e, r2.d[i]); CHECK_F64("md_mm_select_pd", i, e, s2.d[i]); }
            CHECK_F64("md_mm256_blendv_pd", i, e, r4d.d[i]);
            CHECK_F64("md_mm256_select_pd", i, e, s4d.d[i]);
        }
    }
    // The immediate blends: b where the bit is set
    float a[8], b[8];
    double ad[4], bd[4];
    fill_f32(a, 8); fill_f32(b, 8); fill_f64(ad, 4); fill_f64(bd, 4);
    u128 i4 = { .v = md_mm_blend_ps(md_mm_loadu_ps(a), md_mm_loadu_ps(b), MD_SIMD_BLEND_MASK(1, 0, 1, 1)) };
    u256 i8 = { .v = md_mm256_blend_ps(md_mm256_loadu_ps(a), md_mm256_loadu_ps(b), 0xA6) };
    u128d i2 = { .v = md_mm_blend_pd(md_mm_loadu_pd(ad), md_mm_loadu_pd(bd), 2) };
    u256d i4d = { .v = md_mm256_blend_pd(md_mm256_loadu_pd(ad), md_mm256_loadu_pd(bd), 9) };
    for (int i = 0; i < 8; ++i) {
        if (i < 4) CHECK_F32("md_mm_blend_ps", i, (0xD >> i) & 1 ? b[i] : a[i], i4.f[i]);
        CHECK_F32("md_mm256_blend_ps", i, (0xA6 >> i) & 1 ? b[i] : a[i], i8.f[i]);
        if (i < 2) CHECK_F64("md_mm_blend_pd", i, (2 >> i) & 1 ? bd[i] : ad[i], i2.d[i]);
        if (i < 4) CHECK_F64("md_mm256_blend_pd", i, (9 >> i) & 1 ? bd[i] : ad[i], i4d.d[i]);
    }
    EXPECT_EQ(0, fails);
}

// INTEGERS

UTEST(simd_ops, integer) {
    rng = 7;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        int32_t a[8], b[8];
        int64_t al[4], bl[4];
        for (int i = 0; i < 8; ++i) {
            // Small values half the time, so that comparisons see equal and near lanes
            a[i] = (c & 2) ? (int32_t)(rnd_u32() % 7) - 3 : (int32_t)rnd_u32();
            b[i] = (c & 2) ? (int32_t)(rnd_u32() % 7) - 3 : (int32_t)rnd_u32();
        }
        for (int i = 0; i < 4; ++i) {
            al[i] = (c & 2) ? (int64_t)(rnd_u32() % 7) - 3 : (int64_t)rnd_u64();
            bl[i] = (c & 2) ? (int64_t)(rnd_u32() % 7) - 3 : (int64_t)rnd_u64();
        }
        const md_128i a4 = md_mm_loadu_si128(a), b4 = md_mm_loadu_si128(b);
        const md_256i a8 = md_mm256_loadu_si256(a), b8 = md_mm256_loadu_si256(b);
        const md_128i al2 = md_mm_loadu_si128(al), bl2 = md_mm_loadu_si128(bl);
        const md_256i al4 = md_mm256_loadu_si256(al), bl4 = md_mm256_loadu_si256(bl);

#define I32_CASE(OP4, OP8, EXPR) do { \
            u128i r4 = { .v = OP4 }; u256i r8 = { .v = OP8 }; \
            for (int i = 0; i < 8; ++i) { const uint32_t x = (uint32_t)a[i], y = (uint32_t)b[i]; (void)x; (void)y; \
                if (i < 4) CHECK_I64(#OP4, i, (int32_t)(EXPR), r4.i[i]); CHECK_I64(#OP8, i, (int32_t)(EXPR), r8.i[i]); } } while (0)
#define I64_CASE(OP2, OP4, EXPR) do { \
            u128i r2 = { .v = OP2 }; u256i r4 = { .v = OP4 }; \
            for (int i = 0; i < 4; ++i) { const uint64_t x = (uint64_t)al[i], y = (uint64_t)bl[i]; (void)x; (void)y; \
                if (i < 2) CHECK_I64(#OP2, i, (int64_t)(EXPR), r2.l[i]); CHECK_I64(#OP4, i, (int64_t)(EXPR), r4.l[i]); } } while (0)

        I32_CASE(md_mm_add_epi32(a4, b4), md_mm256_add_epi32(a8, b8), x + y);
        I32_CASE(md_mm_sub_epi32(a4, b4), md_mm256_sub_epi32(a8, b8), x - y);
        I32_CASE(md_mm_slli_epi32(a4, 5), md_mm256_slli_epi32(a8, 5), x << 5);
        I32_CASE(md_mm_srli_epi32(a4, 7), md_mm256_srli_epi32(a8, 7), x >> 7);
        I32_CASE(md_mm_and_si128(a4, b4), md_mm256_and_si256(a8, b8), x & y);
        I32_CASE(md_mm_and_epi32(a4, b4), md_mm256_and_epi32(a8, b8), x & y);
        I32_CASE(md_mm_or_si128(a4, b4),  md_mm256_or_si256(a8, b8),  x | y);
        I32_CASE(md_mm_xor_si128(a4, b4), md_mm256_xor_si256(a8, b8), x ^ y);
        I32_CASE(md_mm_andnot_si128(a4, b4), md_mm256_andnot_si256(a8, b8), ~x & y);
        I32_CASE(md_mm_not_si128(a4), md_mm256_not_si256(a8), ~x);
        I32_CASE(md_mm_cmpeq_epi32(a4, b4),  md_mm256_cmpeq_epi32(a8, b8),  a[i] == b[i] ? -1 : 0);
        I32_CASE(md_mm_cmpneq_epi32(a4, b4), md_mm256_cmpneq_epi32(a8, b8), a[i] != b[i] ? -1 : 0);
        I32_CASE(md_mm_cmpgt_epi32(a4, b4),  md_mm256_cmpgt_epi32(a8, b8),  a[i] >  b[i] ? -1 : 0);
        I32_CASE(md_mm_cmplt_epi32(a4, b4),  md_mm256_cmplt_epi32(a8, b8),  a[i] <  b[i] ? -1 : 0);

        I64_CASE(md_mm_add_epi64(al2, bl2), md_mm256_add_epi64(al4, bl4), x + y);
        I64_CASE(md_mm_sub_epi64(al2, bl2), md_mm256_sub_epi64(al4, bl4), x - y);
        I64_CASE(md_mm_slli_epi64(al2, 11), md_mm256_slli_epi64(al4, 11), x << 11);
        I64_CASE(md_mm_srli_epi64(al2, 13), md_mm256_srli_epi64(al4, 13), x >> 13);
        I64_CASE(md_mm_and_epi64(al2, bl2), md_mm256_and_epi64(al4, bl4), x & y);
        I64_CASE(md_mm_cmpeq_epi64(al2, bl2),  md_mm256_cmpeq_epi64(al4, bl4),  al[i] == bl[i] ? -1 : 0);
        I64_CASE(md_mm_cmpneq_epi64(al2, bl2), md_mm256_cmpneq_epi64(al4, bl4), al[i] != bl[i] ? -1 : 0);
        I64_CASE(md_mm_cmpgt_epi64(al2, bl2),  md_mm256_cmpgt_epi64(al4, bl4),  al[i] >  bl[i] ? -1 : 0);
        I64_CASE(md_mm_cmplt_epi64(al2, bl2),  md_mm256_cmplt_epi64(al4, bl4),  al[i] <  bl[i] ? -1 : 0);

        // The 8 and 16 bit comparisons, lane by lane over the same bytes
        u128i x4 = { .v = a4 }, y4 = { .v = b4 };
        u256i x8 = { .v = a8 }, y8 = { .v = b8 };
        u128i e8 = { .v = md_mm_cmpeq_epi8(a4, b4) },  g8 = { .v = md_mm_cmpgt_epi8(a4, b4) },  l8 = { .v = md_mm_cmplt_epi8(a4, b4) };
        u128i e16 = { .v = md_mm_cmpeq_epi16(a4, b4) }, g16 = { .v = md_mm_cmpgt_epi16(a4, b4) }, l16 = { .v = md_mm_cmplt_epi16(a4, b4) };
        u256i E8 = { .v = md_mm256_cmpeq_epi8(a8, b8) },  G8 = { .v = md_mm256_cmpgt_epi8(a8, b8) },  L8 = { .v = md_mm256_cmplt_epi8(a8, b8) };
        u256i E16 = { .v = md_mm256_cmpeq_epi16(a8, b8) }, G16 = { .v = md_mm256_cmpgt_epi16(a8, b8) }, L16 = { .v = md_mm256_cmplt_epi16(a8, b8) };
        for (int i = 0; i < 32; ++i) {
            if (i < 16) {
                CHECK_I64("md_mm_cmpeq_epi8", i, x4.b[i] == y4.b[i] ? -1 : 0, e8.b[i]);
                CHECK_I64("md_mm_cmpgt_epi8", i, x4.b[i] >  y4.b[i] ? -1 : 0, g8.b[i]);
                CHECK_I64("md_mm_cmplt_epi8", i, x4.b[i] <  y4.b[i] ? -1 : 0, l8.b[i]);
            }
            CHECK_I64("md_mm256_cmpeq_epi8", i, x8.b[i] == y8.b[i] ? -1 : 0, E8.b[i]);
            CHECK_I64("md_mm256_cmpgt_epi8", i, x8.b[i] >  y8.b[i] ? -1 : 0, G8.b[i]);
            CHECK_I64("md_mm256_cmplt_epi8", i, x8.b[i] <  y8.b[i] ? -1 : 0, L8.b[i]);
            if (i < 8) {
                CHECK_I64("md_mm_cmpeq_epi16", i, x4.s[i] == y4.s[i] ? -1 : 0, e16.s[i]);
                CHECK_I64("md_mm_cmpgt_epi16", i, x4.s[i] >  y4.s[i] ? -1 : 0, g16.s[i]);
                CHECK_I64("md_mm_cmplt_epi16", i, x4.s[i] <  y4.s[i] ? -1 : 0, l16.s[i]);
            }
            if (i < 16) {
                CHECK_I64("md_mm256_cmpeq_epi16", i, x8.s[i] == y8.s[i] ? -1 : 0, E16.s[i]);
                CHECK_I64("md_mm256_cmpgt_epi16", i, x8.s[i] >  y8.s[i] ? -1 : 0, G16.s[i]);
                CHECK_I64("md_mm256_cmplt_epi16", i, x8.s[i] <  y8.s[i] ? -1 : 0, L16.s[i]);
            }
        }

        // blendv on bytes: the sign bit of each byte of the mask
        u128i bv4 = { .v = md_mm_blendv_epi8(a4, b4, md_mm_loadu_si128(bl)) };
        u256i bv8 = { .v = md_mm256_blendv_epi8(a8, b8, md_mm256_loadu_si256(bl)) };
        u256i m8 = { .v = md_mm256_loadu_si256(bl) };
        for (int i = 0; i < 32; ++i) {
            if (i < 16) CHECK_I64("md_mm_blendv_epi8", i, m8.b[i] < 0 ? y4.b[i] : x4.b[i], bv4.b[i]);
            CHECK_I64("md_mm256_blendv_epi8", i, m8.b[i] < 0 ? y8.b[i] : x8.b[i], bv8.b[i]);
        }

        // movemask on bytes
        int mm4 = 0, mm8 = 0;
        for (int i = 0; i < 32; ++i) {
            if (i < 16 && x4.b[i] < 0) mm4 |= 1 << i;
            if (x8.b[i] < 0) mm8 |= (int)(1u << i);
        }
        CHECK_I64("md_mm_movemask_epi8", 0, mm4, md_mm_movemask_epi8(a4));
        CHECK_I64("md_mm256_movemask_epi8", 0, mm8, md_mm256_movemask_epi8(a8));
    }
    EXPECT_EQ(0, fails);
}

// LOADS, STORES, SETS AND GATHERS

UTEST(simd_ops, load_store_set) {
    int fails = 0;
    ALIGNAS(32) float   f[8];
    ALIGNAS(32) double  d[4];
    ALIGNAS(32) int32_t i32[8];
    ALIGNAS(32) int64_t i64[4];
    for (int i = 0; i < 8; ++i) { f[i] = (float)(i + 1) * 1.5f; i32[i] = -3 * i + 7; }
    for (int i = 0; i < 4; ++i) { d[i] = (double)(i + 1) * -2.25; i64[i] = ((int64_t)1 << 40) * (i + 1) - i; }

    ALIGNAS(32) float   fo[8];
    ALIGNAS(32) double  dout[4];
    ALIGNAS(32) int32_t io[8];
    ALIGNAS(32) int64_t lo[4];

#define ROUND_TRIP(store, load, out, in, n) do { memset(out, 0, sizeof(out)); store(out, load(in)); if (memcmp(out, in, n) != 0) { printf("%s/%s\n", #store, #load); fails += 1; } } while (0)
    ROUND_TRIP(md_mm_storeu_ps,    md_mm_loadu_ps,    fo, f, 16);
    ROUND_TRIP(md_mm_store_ps,     md_mm_load_ps,     fo, f, 16);
    ROUND_TRIP(md_mm256_storeu_ps, md_mm256_loadu_ps, fo, f, 32);
    ROUND_TRIP(md_mm256_store_ps,  md_mm256_load_ps,  fo, f, 32);
    ROUND_TRIP(md_mm_storeu_pd,    md_mm_loadu_pd,    dout, d, 16);
    ROUND_TRIP(md_mm_store_pd,     md_mm_load_pd,     dout, d, 16);
    ROUND_TRIP(md_mm256_storeu_pd, md_mm256_loadu_pd, dout, d, 32);
    ROUND_TRIP(md_mm256_store_pd,  md_mm256_load_pd,  dout, d, 32);
#define STOREU_SI128(p, v) md_mm_storeu_si128((md_128i*)(p), v)
#define STORE_SI128(p, v)  md_mm_store_si128((md_128i*)(p), v)
#define LOAD_SI128(p)      md_mm_load_si128((const md_128i*)(p))
#define STOREU_SI256(p, v) md_mm256_storeu_si256((md_256i*)(p), v)
#define STORE_SI256(p, v)  md_mm256_store_si256((md_256i*)(p), v)
#define LOAD_SI256(p)      md_mm256_load_si256((const md_256i*)(p))
#define STOREU_EPI32(p, v)    md_mm_storeu_epi32((md_128i*)(p), v)
#define STORE_EPI32(p, v)     md_mm_store_epi32((md_128i*)(p), v)
#define LOAD_EPI32(p)         md_mm_load_epi32((const md_128i*)(p))
#define STOREU_EPI64(p, v)    md_mm_storeu_epi64((md_128i*)(p), v)
#define STORE_EPI64(p, v)     md_mm_store_epi64((md_128i*)(p), v)
#define LOAD_EPI64(p)         md_mm_load_epi64((const md_128i*)(p))
#define STOREU256_EPI32(p, v) md_mm256_storeu_epi32((md_256i*)(p), v)
#define STORE256_EPI32(p, v)  md_mm256_store_epi32((md_256i*)(p), v)
#define LOAD256_EPI32(p)      md_mm256_load_epi32((const md_256i*)(p))
#define STOREU256_EPI64(p, v) md_mm256_storeu_epi64((md_256i*)(p), v)
#define STORE256_EPI64(p, v)  md_mm256_store_epi64((md_256i*)(p), v)
#define LOAD256_EPI64(p)      md_mm256_load_epi64((const md_256i*)(p))
    ROUND_TRIP(STOREU_SI128,    md_mm_loadu_si128,    io, i32, 16);
    ROUND_TRIP(STORE_SI128,     LOAD_SI128,           io, i32, 16);
    ROUND_TRIP(STOREU_EPI32,    md_mm_loadu_epi32,    io, i32, 16);
    ROUND_TRIP(STORE_EPI32,     LOAD_EPI32,           io, i32, 16);
    ROUND_TRIP(STOREU_EPI64,    md_mm_loadu_epi64,    lo, i64, 16);
    ROUND_TRIP(STORE_EPI64,     LOAD_EPI64,           lo, i64, 16);
    ROUND_TRIP(STOREU_SI256,    md_mm256_loadu_si256, io, i32, 32);
    ROUND_TRIP(STORE_SI256,     LOAD_SI256,           io, i32, 32);
    ROUND_TRIP(STOREU256_EPI32, md_mm256_loadu_epi32, io, i32, 32);
    ROUND_TRIP(STORE256_EPI32,  LOAD256_EPI32,        io, i32, 32);
    ROUND_TRIP(STOREU256_EPI64, md_mm256_loadu_epi64, lo, i64, 32);
    ROUND_TRIP(STORE256_EPI64,  LOAD256_EPI64,        lo, i64, 32);

    // set takes the lanes high to low, setr low to high
    u128 s4 = { .v = md_mm_set_ps(4, 3, 2, 1) };
    u256 s8 = { .v = md_mm256_set_ps(8, 7, 6, 5, 4, 3, 2, 1) };
    u256 r8 = { .v = md_mm256_setr_ps(1, 2, 3, 4, 5, 6, 7, 8) };
    u128d s2 = { .v = md_mm_set_pd(2, 1) };
    u256d s4d = { .v = md_mm256_set_pd(4, 3, 2, 1) };
    u128i si4 = { .v = md_mm_set_epi32(4, 3, 2, 1) };
    u256i si8 = { .v = md_mm256_set_epi32(8, 7, 6, 5, 4, 3, 2, 1) };
    u128i sl2 = { .v = md_mm_set_epi64(2, 1) };
    u256i sl4 = { .v = md_mm256_set_epi64(4, 3, 2, 1) };
    u128i ss8 = { .v = md_mm_set_epi16(8, 7, 6, 5, 4, 3, 2, 1) };
    u128i sb16 = { .v = md_mm_set_epi8(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1) };
    for (int i = 0; i < 16; ++i) {
        if (i < 4) { CHECK_F32("md_mm_set_ps", i, (float)(i + 1), s4.f[i]); CHECK_I64("md_mm_set_epi32", i, i + 1, si4.i[i]); }
        if (i < 8) {
            CHECK_F32("md_mm256_set_ps", i, (float)(i + 1), s8.f[i]);
            CHECK_F32("md_mm256_setr_ps", i, (float)(i + 1), r8.f[i]);
            CHECK_I64("md_mm256_set_epi32", i, i + 1, si8.i[i]);
            CHECK_I64("md_mm_set_epi16", i, i + 1, ss8.s[i]);
        }
        if (i < 2) { CHECK_F64("md_mm_set_pd", i, (double)(i + 1), s2.d[i]); CHECK_I64("md_mm_set_epi64", i, i + 1, sl2.l[i]); }
        if (i < 4) { CHECK_F64("md_mm256_set_pd", i, (double)(i + 1), s4d.d[i]); CHECK_I64("md_mm256_set_epi64", i, i + 1, sl4.l[i]); }
        CHECK_I64("md_mm_set_epi8", i, i + 1, sb16.b[i]);
    }

    // set1 and setzero
    u128 o4 = { .v = md_mm_set1_ps(2.5f) };
    u256 o8 = { .v = md_mm256_set1_ps(2.5f) };
    u128d o2 = { .v = md_mm_set1_pd(-1.25) };
    u256d o4d = { .v = md_mm256_set1_pd(-1.25) };
    u128i b16 = { .v = md_mm_set1_epi8(-5) }, h8 = { .v = md_mm_set1_epi16(-300) }, w4 = { .v = md_mm_set1_epi32(-70000) }, q2 = { .v = md_mm_set1_epi64(-((int64_t)1 << 35)) };
    u256i w8 = { .v = md_mm256_set1_epi32(-70000) }, q4 = { .v = md_mm256_set1_epi64(-((int64_t)1 << 35)) };
    u128 z4 = { .v = md_mm_setzero_ps() };
    u256 z8 = { .v = md_mm256_setzero_ps() };
    u128d z2 = { .v = md_mm_setzero_pd() };
    u256d z4d = { .v = md_mm256_setzero_pd() };
    u128i zi = { .v = md_mm_setzero_si128() };
    u256i zi8 = { .v = md_mm256_setzero_si256() };
    for (int i = 0; i < 16; ++i) {
        if (i < 4) { CHECK_F32("md_mm_set1_ps", i, 2.5f, o4.f[i]); CHECK_F32("md_mm_setzero_ps", i, 0.0f, z4.f[i]); CHECK_I64("md_mm_set1_epi32", i, -70000, w4.i[i]); CHECK_I64("md_mm_setzero_si128", i, 0, zi.i[i]); }
        if (i < 8) { CHECK_F32("md_mm256_set1_ps", i, 2.5f, o8.f[i]); CHECK_F32("md_mm256_setzero_ps", i, 0.0f, z8.f[i]); CHECK_I64("md_mm256_set1_epi32", i, -70000, w8.i[i]); CHECK_I64("md_mm256_setzero_si256", i, 0, zi8.i[i]); CHECK_I64("md_mm_set1_epi16", i, -300, h8.s[i]); }
        if (i < 2) { CHECK_F64("md_mm_set1_pd", i, -1.25, o2.d[i]); CHECK_F64("md_mm_setzero_pd", i, 0.0, z2.d[i]); CHECK_I64("md_mm_set1_epi64", i, -((int64_t)1 << 35), q2.l[i]); }
        if (i < 4) { CHECK_F64("md_mm256_set1_pd", i, -1.25, o4d.d[i]); CHECK_F64("md_mm256_setzero_pd", i, 0.0, z4d.d[i]); CHECK_I64("md_mm256_set1_epi64", i, -((int64_t)1 << 35), q4.l[i]); }
        CHECK_I64("md_mm_set1_epi8", i, -5, b16.b[i]);
    }

    // maskstore: only where the mask is set
    for (int i = 0; i < 8; ++i) fo[i] = -1.0f;
    for (int i = 0; i < 4; ++i) dout[i] = -1.0;
    md_mm_maskstore_ps(fo, md_mm_set_epi32(-1, 0, -1, 0), md_mm_loadu_ps(f));
    md_mm_maskstore_pd(dout, md_mm_set_epi64(-1, 0), md_mm_loadu_pd(d));
    md_mm256_maskstore_pd(dout, md_mm256_set_epi64(0, -1, 0, 0), md_mm256_loadu_pd(d));
    CHECK_F32("md_mm_maskstore_ps", 0, -1.0f, fo[0]);
    CHECK_F32("md_mm_maskstore_ps", 1, f[1], fo[1]);
    CHECK_F32("md_mm_maskstore_ps", 2, -1.0f, fo[2]);
    CHECK_F32("md_mm_maskstore_ps", 3, f[3], fo[3]);
    CHECK_F32("md_mm_maskstore_ps", 4, -1.0f, fo[4]);
    for (int i = 0; i < 8; ++i) fo[i] = -1.0f;
    md_mm256_maskstore_ps(fo, md_mm256_set_epi32(0, -1, 0, 0, 0, 0, 0, -1), md_mm256_loadu_ps(f));
    for (int i = 0; i < 8; ++i) CHECK_F32("md_mm256_maskstore_ps", i, (i == 0 || i == 6) ? f[i] : -1.0f, fo[i]);
    CHECK_F64("md_mm_maskstore_pd", 0, -1.0, dout[0]);
    CHECK_F64("md_mm_maskstore_pd", 1, d[1], dout[1]);
    CHECK_F64("md_mm256_maskstore_pd", 2, d[2], dout[2]);
    CHECK_F64("md_mm256_maskstore_pd", 3, -1.0, dout[3]);
    EXPECT_EQ(0, fails);
}

UTEST(simd_ops, gather) {
    int fails = 0;
    float f[64];
    double d[64];
    int32_t i32[64];
    int64_t i64[64];
    for (int i = 0; i < 64; ++i) { f[i] = (float)i * 0.5f; d[i] = (double)i * -0.25; i32[i] = i * 3 - 50; i64[i] = ((int64_t)i << 33) - i; }
    const int32_t idx[8] = { 63, 0, 17, 5, 5, 42, 1, 30 };
    const md_128i ix4 = md_mm_loadu_si128(idx);
    const md_256i ix8 = md_mm256_loadu_si256(idx);
    u128 gf4 = { .v = md_mm_i32gather_ps((const float*)f, ix4, 4) };
    u256 gf8 = { .v = md_mm256_i32gather_ps((const float*)f, ix8, 4) };
    u128d gd2 = { .v = md_mm_i32gather_pd((const double*)d, ix4, 8) };
    u256d gd4 = { .v = md_mm256_i32gather_pd((const double*)d, ix4, 8) };
    u128i gi4 = { .v = md_mm_i32gather_epi32((const int32_t*)i32, ix4, 4) };
    u256i gi8 = { .v = md_mm256_i32gather_epi32((const int32_t*)i32, ix8, 4) };
    u128i gl2 = { .v = md_mm_i32gather_epi64((const int64_t*)i64, ix4, 8) };
    u256i gl4 = { .v = md_mm256_i32gather_epi64((const int64_t*)i64, ix4, 8) };
    // The masked gathers keep src where the mask is clear
    u128 mf4 = { .v = md_mm_mask_i32gather_ps(md_mm_set1_ps(-9.0f), (const float*)f, ix4, md_mm_castsi128_ps(md_mm_set_epi32(-1, 0, -1, 0)), 4) };
    u256 mf8 = { .v = md_mm256_mask_i32gather_ps(md_mm256_set1_ps(-9.0f), (const float*)f, ix8, md_mm256_castsi256_ps(md_mm256_set_epi32(-1, 0, 0, -1, -1, 0, -1, 0)), 4) };
    u128d md2 = { .v = md_mm_mask_i32gather_pd(md_mm_set1_pd(-9.0), (const double*)d, ix4, md_mm_castsi128_pd(md_mm_set_epi64(-1, 0)), 8) };
    u256d md4 = { .v = md_mm256_mask_i32gather_pd(md_mm256_set1_pd(-9.0), (const double*)d, ix4, md_mm256_castsi256_pd(md_mm256_set_epi64(0, -1, -1, 0)), 8) };
    for (int i = 0; i < 8; ++i) {
        if (i < 4) {
            CHECK_F32("md_mm_i32gather_ps", i, f[idx[i]], gf4.f[i]);
            CHECK_I64("md_mm_i32gather_epi32", i, i32[idx[i]], gi4.i[i]);
            CHECK_F64("md_mm256_i32gather_pd", i, d[idx[i]], gd4.d[i]);
            CHECK_I64("md_mm256_i32gather_epi64", i, i64[idx[i]], gl4.l[i]);
            CHECK_F32("md_mm_mask_i32gather_ps", i, (0xA >> i) & 1 ? f[idx[i]] : -9.0f, mf4.f[i]);
            CHECK_F64("md_mm256_mask_i32gather_pd", i, (0x6 >> i) & 1 ? d[idx[i]] : -9.0, md4.d[i]);
        }
        if (i < 2) {
            CHECK_F64("md_mm_i32gather_pd", i, d[idx[i]], gd2.d[i]);
            CHECK_I64("md_mm_i32gather_epi64", i, i64[idx[i]], gl2.l[i]);
            CHECK_F64("md_mm_mask_i32gather_pd", i, (0x2 >> i) & 1 ? d[idx[i]] : -9.0, md2.d[i]);
        }
        CHECK_F32("md_mm256_i32gather_ps", i, f[idx[i]], gf8.f[i]);
        CHECK_I64("md_mm256_i32gather_epi32", i, i32[idx[i]], gi8.i[i]);
        CHECK_F32("md_mm256_mask_i32gather_ps", i, (0x9A >> i) & 1 ? f[idx[i]] : -9.0f, mf8.f[i]);
    }
    EXPECT_EQ(0, fails);
}

// CASTS AND CONVERSIONS

UTEST(simd_ops, casts) {
    rng = 8;
    int fails = 0;
    float f[8];
    fill_f32(f, 8);
    const md_128 v4 = md_mm_loadu_ps(f);
    const md_256 v8 = md_mm256_loadu_ps(f);
    // Bit casts keep the bits
    u128 a = { .v = md_mm_castsi128_ps(md_mm_castps_si128(v4)) };
    u128 b = { .v = md_mm_castpd_ps(md_mm_castps_pd(v4)) };
    u128 c = { .v = md_mm_castpd_ps(md_mm_castsi128_pd(md_mm_castpd_si128(md_mm_castps_pd(v4)))) };
    u256 A = { .v = md_mm256_castsi256_ps(md_mm256_castps_si256(v8)) };
    u256 B = { .v = md_mm256_castpd_ps(md_mm256_castps_pd(v8)) };
    u256 C = { .v = md_mm256_castpd_ps(md_mm256_castsi256_pd(md_mm256_castpd_si256(md_mm256_castps_pd(v8)))) };
    for (int i = 0; i < 8; ++i) {
        if (i < 4) { CHECK_F32("md_mm_cast*_ps/si128", i, f[i], a.f[i]); CHECK_F32("md_mm_cast*_ps/pd", i, f[i], b.f[i]); CHECK_F32("md_mm_cast*_pd/si128", i, f[i], c.f[i]); }
        CHECK_F32("md_mm256_cast*_ps/si256", i, f[i], A.f[i]);
        CHECK_F32("md_mm256_cast*_ps/pd", i, f[i], B.f[i]);
        CHECK_F32("md_mm256_cast*_pd/si256", i, f[i], C.f[i]);
    }
    CHECK_F32("md_mm_cvtss_f32", 0, f[0], md_mm_cvtss_f32(v4));
    CHECK_F32("md_mm256_cvtss_f32", 0, f[0], md_mm256_cvtss_f32(v8));
    CHECK_F64("md_mm_cvtsd_f64", 0, -3.5, md_mm_cvtsd_f64(md_mm_set_pd(7.0, -3.5)));
    CHECK_F64("md_mm256_cvtsd_f64", 0, -3.5, md_mm256_cvtsd_f64(md_mm256_set_pd(1.0, 2.0, 7.0, -3.5)));
    EXPECT_EQ(0, fails);
}

UTEST(simd_ops, convert) {
    rng = 9;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        float f[8];
        int32_t n[8];
        fill_f32(f, 8);
        for (int i = 0; i < 8; ++i) {
            f[i] = (c & 1) ? HALVES(f[i]) : f[i] * 1.0e6f;   // ties, and values up to 5e8
            n[i] = (int32_t)rnd_u32() >> (rnd_u32() % 31);
        }
        u128i t4 = { .v = md_mm_cvttps_epi32(md_mm_loadu_ps(f)) };
        u256i t8 = { .v = md_mm256_cvttps_epi32(md_mm256_loadu_ps(f)) };
        u128i r4 = { .v = md_mm_cvtps_epi32(md_mm_loadu_ps(f)) };
        u256i r8 = { .v = md_mm256_cvtps_epi32(md_mm256_loadu_ps(f)) };
        u128 i4 = { .v = md_mm_cvtepi32_ps(md_mm_loadu_si128(n)) };
        u256 i8 = { .v = md_mm256_cvtepi32_ps(md_mm256_loadu_si256(n)) };
        for (int i = 0; i < 8; ++i) {
            const int32_t trunc = (int32_t)f[i];
            const int32_t round = (int32_t)rint((double)f[i]);
            const float back = (float)n[i];
            if (i < 4) { CHECK_I64("md_mm_cvttps_epi32", i, trunc, t4.i[i]); CHECK_I64("md_mm_cvtps_epi32", i, round, r4.i[i]); CHECK_F32("md_mm_cvtepi32_ps", i, back, i4.f[i]); }
            CHECK_I64("md_mm256_cvttps_epi32", i, trunc, t8.i[i]);
            CHECK_I64("md_mm256_cvtps_epi32", i, round, r8.i[i]);
            CHECK_F32("md_mm256_cvtepi32_ps", i, back, i8.f[i]);
        }
    }
    // Bytes to int32 lanes, zero extended
    u128i e4 = { .v = md_mm_cvtepu8_epi32(md_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (char)200, 3, (char)255, 1)) };
    u256i e8 = { .v = md_mm256_cvtepu8_epi32(md_mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 9, 8, 7, 6, (char)200, 3, (char)255, 1)) };
    const int32_t ext[8] = { 1, 255, 3, 200, 6, 7, 8, 9 };
    for (int i = 0; i < 8; ++i) {
        if (i < 4) CHECK_I64("md_mm_cvtepu8_epi32", i, ext[i], e4.i[i]);
        CHECK_I64("md_mm256_cvtepu8_epi32", i, ext[i], e8.i[i]);
    }
    u128i q = { .v = md_mm_cvtsi64_si128(-((int64_t)1 << 40)) };
    CHECK_I64("md_mm_cvtsi64_si128", 0, -((int64_t)1 << 40), q.l[0]);
    CHECK_I64("md_mm_cvtsi64_si128", 1, 0, q.l[1]);
    EXPECT_EQ(0, fails);
}

// SHUFFLES AND MASKS

UTEST(simd_ops, shuffle_and_permute) {
    int fails = 0;
    const md_128 a = md_mm_set_ps(3, 2, 1, 0), b = md_mm_set_ps(7, 6, 5, 4);
    const md_256 A = md_mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0), B = md_mm256_set_ps(15, 14, 13, 12, 11, 10, 9, 8);
    // shuffle: lanes 0 and 1 from a, 2 and 3 from b, in each 128 bit half
    u128 s = { .v = md_mm_shuffle_ps(a, b, MD_SHUFFLE(1, 3, 0, 2)) };
    u256 S = { .v = md_mm256_shuffle_ps(A, B, MD_SHUFFLE(1, 3, 0, 2)) };
    const float es[4] = { 2, 0, 7, 5 };
    const float eS[8] = { 2, 0, 11, 9, 6, 4, 15, 13 };
    for (int i = 0; i < 8; ++i) {
        if (i < 4) CHECK_F32("md_mm_shuffle_ps", i, es[i], s.f[i]);
        CHECK_F32("md_mm256_shuffle_ps", i, eS[i], S.f[i]);
    }
    u128d sd = { .v = md_mm_shuffle_pd(md_mm_set_pd(1, 0), md_mm_set_pd(3, 2), 1) };
    u256d SD = { .v = md_mm256_shuffle_pd(md_mm256_set_pd(3, 2, 1, 0), md_mm256_set_pd(7, 6, 5, 4), 0x6) };
    CHECK_F64("md_mm_shuffle_pd", 0, 1.0, sd.d[0]);
    CHECK_F64("md_mm_shuffle_pd", 1, 2.0, sd.d[1]);
    const double eSD[4] = { 0, 5, 3, 6 };
    for (int i = 0; i < 4; ++i) CHECK_F64("md_mm256_shuffle_pd", i, eSD[i], SD.d[i]);
    // splat
    for (int k = 0; k < 4; ++k) {
        u128 sp;
        switch (k) {
        case 0: sp.v = md_mm_splat_ps(a, 0); break;
        case 1: sp.v = md_mm_splat_ps(a, 1); break;
        case 2: sp.v = md_mm_splat_ps(a, 2); break;
        default: sp.v = md_mm_splat_ps(a, 3); break;
        }
        for (int i = 0; i < 4; ++i) CHECK_F32("md_mm_splat_ps", i, (float)k, sp.f[i]);
    }
    u128d sp0 = { .v = md_mm_splat_pd(md_mm_set_pd(1, 0), 0) }, sp1 = { .v = md_mm_splat_pd(md_mm_set_pd(1, 0), 1) };
    CHECK_F64("md_mm_splat_pd", 1, 0.0, sp0.d[1]);
    CHECK_F64("md_mm_splat_pd", 0, 1.0, sp1.d[0]);
    // movehl: the high half of b, then the high half of a
    u128 hl = { .v = md_mm_movehl_ps(a, b) };
    const float ehl[4] = { 6, 7, 2, 3 };
    for (int i = 0; i < 4; ++i) CHECK_F32("md_mm_movehl_ps", i, ehl[i], hl.f[i]);
    // permutevar8x32: out[i] = in[idx[i]]
    const int32_t perm[8] = { 7, 0, 3, 3, 6, 1, 5, 2 };
    u256 pv = { .v = md_mm256_permutevar8x32_ps(A, md_mm256_loadu_si256(perm)) };
    u256i pi = { .v = md_mm256_permutevar8x32_epi32(md_mm256_set_epi32(70, 60, 50, 40, 30, 20, 10, 0), md_mm256_loadu_si256(perm)) };
    for (int i = 0; i < 8; ++i) {
        CHECK_F32("md_mm256_permutevar8x32_ps", i, (float)perm[i], pv.f[i]);
        CHECK_I64("md_mm256_permutevar8x32_epi32", i, perm[i] * 10, pi.i[i]);
    }
    EXPECT_EQ(0, fails);
}

UTEST(simd_ops, movemask_and_compression) {
    int fails = 0;
    for (int m = 0; m < 256; ++m) {
        float f[8];
        double d[4];
        for (int i = 0; i < 8; ++i) f[i] = (m >> i) & 1 ? -1.0f - (float)i : 1.0f + (float)i;
        for (int i = 0; i < 4; ++i) d[i] = (m >> i) & 1 ? -2.0 : 2.0;
        CHECK_I64("md_mm_movemask_ps", m, m & 15, md_mm_movemask_ps(md_mm_loadu_ps(f)));
        CHECK_I64("md_mm256_movemask_ps", m, m, md_mm256_movemask_ps(md_mm256_loadu_ps(f)));
        CHECK_I64("md_mm_movemask_pd", m, m & 3, md_mm_movemask_pd(md_mm_loadu_pd(d)));
        CHECK_I64("md_mm256_movemask_pd", m, m & 15, md_mm256_movemask_pd(md_mm256_loadu_pd(d)));

        // The compression masks list the set lanes in order
        u256i c8 = { .v = md_mm256_compression_mask_8x32(m) };
        int k = 0;
        for (int i = 0; i < 8; ++i) if ((m >> i) & 1) { CHECK_I64("md_mm256_compression_mask_8x32", k, i, c8.i[k]); k += 1; }
        if (m < 16) {
            u128i c4 = { .v = md_mm_compression_mask_4x32(m) };
            k = 0;
            for (int i = 0; i < 4; ++i) if ((m >> i) & 1) { CHECK_I64("md_mm_compression_mask_4x32", k, i, c4.i[k]); k += 1; }
        }
    }
    EXPECT_EQ(0, fails);
}

// HIGHER LEVEL OPERATIONS

static float ref_lerp(float a, float b, float t) { return a * (1.0f - t) + b * t; }
static double ref_cubic_spline(double p0, double p1, double p2, double p3, double t, double s) {
    const double v0 = (p2 - p0) * s;
    const double v1 = (p3 - p1) * s;
    const double x0 = 2 * (p1 - p2) + v0 + v1;
    const double x1 = 3 * (p2 - p1) - 2 * v0 - v1;
    return x0 * t * t * t + x1 * t * t + v0 * t + p1;
}

UTEST(simd_ops, interpolation) {
    rng = 10;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        float p[4][8], t[8];
        double pd[4][4], td[4];
        for (int k = 0; k < 4; ++k) { fill_f32(p[k], 8); fill_f64(pd[k], 4); }
        for (int i = 0; i < 8; ++i) t[i] = (float)(rnd_u32() >> 8) / 16777216.0f;
        for (int i = 0; i < 4; ++i) td[i] = (double)(rnd_u32() >> 8) / 16777216.0;
        const float tt = t[0];
        u128 l4 = { .v = md_mm_lerp_ps(md_mm_loadu_ps(p[0]), md_mm_loadu_ps(p[1]), tt) };
        u256 l8 = { .v = md_mm256_lerp_ps(md_mm256_loadu_ps(p[0]), md_mm256_loadu_ps(p[1]), tt) };
        u128d l2 = { .v = md_mm_lerp_pd(md_mm_loadu_pd(pd[0]), md_mm_loadu_pd(pd[1]), td[0]) };
        u256d l4d = { .v = md_mm256_lerp_pd(md_mm256_loadu_pd(pd[0]), md_mm256_loadu_pd(pd[1]), td[0]) };
        u128 c4 = { .v = md_mm_cubic_spline_ps(md_mm_loadu_ps(p[0]), md_mm_loadu_ps(p[1]), md_mm_loadu_ps(p[2]), md_mm_loadu_ps(p[3]), md_mm_loadu_ps(t), md_mm_set1_ps(0.5f)) };
        u256 c8 = { .v = md_mm256_cubic_spline_ps(md_mm256_loadu_ps(p[0]), md_mm256_loadu_ps(p[1]), md_mm256_loadu_ps(p[2]), md_mm256_loadu_ps(p[3]), md_mm256_loadu_ps(t), md_mm256_set1_ps(0.5f)) };
        u128d c2 = { .v = md_mm_cubic_spline_pd(md_mm_loadu_pd(pd[0]), md_mm_loadu_pd(pd[1]), md_mm_loadu_pd(pd[2]), md_mm_loadu_pd(pd[3]), md_mm_loadu_pd(td), md_mm_set1_pd(0.5)) };
        u256d c4d = { .v = md_mm256_cubic_spline_pd(md_mm256_loadu_pd(pd[0]), md_mm256_loadu_pd(pd[1]), md_mm256_loadu_pd(pd[2]), md_mm256_loadu_pd(pd[3]), md_mm256_loadu_pd(td), md_mm256_set1_pd(0.5)) };
        for (int i = 0; i < 8; ++i) {
            // The tolerances are those of the cancellation in the formulas, relative to the size of the points
            const double scale = fabs(p[0][i]) + fabs(p[1][i]) + fabs(p[2][i]) + fabs(p[3][i]);
            const double el = ref_lerp(p[0][i], p[1][i], tt);
            const double ec = ref_cubic_spline(p[0][i], p[1][i], p[2][i], p[3][i], t[i], 0.5);
            if (i < 4) {
                if (fabs(l4.f[i] - el) > 1e-6 * scale) { FAIL_F32("md_mm_lerp_ps", i, (float)el, l4.f[i]); fails += 1; }
                if (fabs(c4.f[i] - ec) > 4e-6 * scale) { FAIL_F32("md_mm_cubic_spline_ps", i, (float)ec, c4.f[i]); fails += 1; }
            }
            if (fabs(l8.f[i] - el) > 1e-6 * scale) { FAIL_F32("md_mm256_lerp_ps", i, (float)el, l8.f[i]); fails += 1; }
            if (fabs(c8.f[i] - ec) > 4e-6 * scale) { FAIL_F32("md_mm256_cubic_spline_ps", i, (float)ec, c8.f[i]); fails += 1; }
        }
        for (int i = 0; i < 4; ++i) {
            const double scale = fabs(pd[0][i]) + fabs(pd[1][i]) + fabs(pd[2][i]) + fabs(pd[3][i]);
            const double el = pd[0][i] * (1.0 - td[0]) + pd[1][i] * td[0];
            const double ec = ref_cubic_spline(pd[0][i], pd[1][i], pd[2][i], pd[3][i], td[i], 0.5);
            if (i < 2) {
                if (fabs(l2.d[i] - el) > 4e-15 * scale) { FAIL_F64("md_mm_lerp_pd", i, el, l2.d[i]); fails += 1; }
                if (fabs(c2.d[i] - ec) > 1e-14 * scale) { FAIL_F64("md_mm_cubic_spline_pd", i, ec, c2.d[i]); fails += 1; }
            }
            if (fabs(l4d.d[i] - el) > 4e-15 * scale) { FAIL_F64("md_mm256_lerp_pd", i, el, l4d.d[i]); fails += 1; }
            if (fabs(c4d.d[i] - ec) > 1e-14 * scale) { FAIL_F64("md_mm256_cubic_spline_pd", i, ec, c4d.d[i]); fails += 1; }
        }
    }
    EXPECT_EQ(0, fails);
}

// The image of x nearest r for the period p (p = 0: no periodicity), and the minimum image of a difference
UTEST(simd_ops, periodic) {
    rng = 11;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        float x[8], r[8], p[8], rp[8];
        double xd[4], rd[4], pdd[4], rpd[4];
        for (int i = 0; i < 8; ++i) {
            p[i]  = (i % 3 == 2) ? 0.0f : 5.0f + 50.0f * (float)(rnd_u32() >> 8) / 16777216.0f;
            rp[i] = p[i] != 0.0f ? 1.0f / p[i] : 0.0f;
            r[i]  = 100.0f * ((float)(rnd_u32() >> 8) / 16777216.0f - 0.5f);
            x[i]  = r[i] + 300.0f * ((float)(rnd_u32() >> 8) / 16777216.0f - 0.5f);
        }
        for (int i = 0; i < 4; ++i) { xd[i] = x[i]; rd[i] = r[i]; pdd[i] = p[i]; rpd[i] = p[i] != 0.0f ? 1.0 / p[i] : 0.0; }
        const md_128 X4 = md_mm_loadu_ps(x), R4 = md_mm_loadu_ps(r), P4 = md_mm_loadu_ps(p), RP4 = md_mm_loadu_ps(rp);
        const md_256 X8 = md_mm256_loadu_ps(x), R8 = md_mm256_loadu_ps(r), P8 = md_mm256_loadu_ps(p), RP8 = md_mm256_loadu_ps(rp);
        u128 a4 = { .v = md_mm_deperiodize_ps(X4, R4, P4) }, b4 = { .v = md_mm_deperiodize2_ps(X4, R4, P4, RP4) };
        u256 a8 = { .v = md_mm256_deperiodize_ps(X8, R8, P8) }, b8 = { .v = md_mm256_deperiodize2_ps(X8, R8, P8, RP8) };
        u128 m4 = { .v = md_mm_minimage_ps(md_mm_sub_ps(X4, R4), P4, RP4) };
        u256 m8 = { .v = md_mm256_minimage_ps(md_mm256_sub_ps(X8, R8), P8, RP8) };
        u128d a2 = { .v = md_mm_deperiodize_pd(md_mm_loadu_pd(xd), md_mm_loadu_pd(rd), md_mm_loadu_pd(pdd)) };
        u128d b2 = { .v = md_mm_deperiodize2_pd(md_mm_loadu_pd(xd), md_mm_loadu_pd(rd), md_mm_loadu_pd(pdd), md_mm_loadu_pd(rpd)) };
        u256d a4d = { .v = md_mm256_deperiodize_pd(md_mm256_loadu_pd(xd), md_mm256_loadu_pd(rd), md_mm256_loadu_pd(pdd)) };
        u256d b4d = { .v = md_mm256_deperiodize2_pd(md_mm256_loadu_pd(xd), md_mm256_loadu_pd(rd), md_mm256_loadu_pd(pdd), md_mm256_loadu_pd(rpd)) };
        u128d m2 = { .v = md_mm_minimage_pd(md_mm_sub_pd(md_mm_loadu_pd(xd), md_mm_loadu_pd(rd)), md_mm_loadu_pd(pdd), md_mm_loadu_pd(rpd)) };
        u256d m4d = { .v = md_mm256_minimage_pd(md_mm256_sub_pd(md_mm256_loadu_pd(xd), md_mm256_loadu_pd(rd)), md_mm256_loadu_pd(pdd), md_mm256_loadu_pd(rpd)) };
        for (int i = 0; i < 8; ++i) {
            const double d = (double)x[i] - (double)r[i];
            const double img = p[i] != 0.0f ? d - p[i] * rint(d / p[i]) : d;
            const double e = p[i] != 0.0f ? r[i] + img : x[i];
            // |img| <= p / 2; a lane exactly half a period away may round either way
            if (fabs(fabs(img) - 0.5 * p[i]) < 1e-3) continue;
            // The fraction of a period is taken of d / p: it carries the rounding of a number as large as |d| / p, p times
            const double tol = 4e-7 * (fabs(x[i]) + fabs(r[i]) + p[i]) + 1e-6;
#define NEAR_LANE(name, val, ref) if (fabs((val) - (ref)) > tol) { FAIL_F32(name, i, (float)(ref), (float)(val)); fails += 1; }
            if (i < 4) { NEAR_LANE("md_mm_deperiodize_ps", a4.f[i], e); NEAR_LANE("md_mm_deperiodize2_ps", b4.f[i], e); NEAR_LANE("md_mm_minimage_ps", m4.f[i], img); }
            NEAR_LANE("md_mm256_deperiodize_ps", a8.f[i], e);
            NEAR_LANE("md_mm256_deperiodize2_ps", b8.f[i], e);
            NEAR_LANE("md_mm256_minimage_ps", m8.f[i], img);
            if (i < 4) {
                if (i < 2) { NEAR_LANE("md_mm_deperiodize_pd", a2.d[i], e); NEAR_LANE("md_mm_deperiodize2_pd", b2.d[i], e); NEAR_LANE("md_mm_minimage_pd", m2.d[i], img); }
                NEAR_LANE("md_mm256_deperiodize_pd", a4d.d[i], e);
                NEAR_LANE("md_mm256_deperiodize2_pd", b4d.d[i], e);
                NEAR_LANE("md_mm256_minimage_pd", m4d.d[i], img);
            }
            // No period: x as it is, bit for bit
            if (p[i] == 0.0f) {
                if (i < 4) CHECK_F32("md_mm_deperiodize_ps", i, x[i], a4.f[i]);
                CHECK_F32("md_mm256_deperiodize2_ps", i, x[i], b8.f[i]);
            }
        }
    }
    EXPECT_EQ(0, fails);
}

UTEST(simd_ops, unpack_strided) {
    int fails = 0;
    // Records of five floats (20 bytes): x, y, z, w and one more
    float rec[8 * 5 + 4];
    for (int i = 0; i < 8 * 5 + 4; ++i) rec[i] = (float)i;
    md_128 x4, y4, z4, w4;
    md_256 x8, y8, z8, w8;
    md_mm_unpack_xyz_ps(&x4, &y4, &z4, rec, 5 * sizeof(float));
    u128 ux = { .v = x4 }, uy = { .v = y4 }, uz = { .v = z4 };
    md_mm_unpack_xyzw_ps(&x4, &y4, &z4, &w4, rec, 5 * sizeof(float));
    u128 uw = { .v = w4 };
    md_mm256_unpack_xyz_ps(&x8, &y8, &z8, rec, 5 * sizeof(float));
    u256 Ux = { .v = x8 }, Uy = { .v = y8 }, Uz = { .v = z8 };
    md_mm256_unpack_xyzw_ps(&x8, &y8, &z8, &w8, rec, 5 * sizeof(float));
    u256 Uw = { .v = w8 };
    for (int i = 0; i < 8; ++i) {
        if (i < 4) {
            CHECK_F32("md_mm_unpack_xyz_ps x", i, rec[5 * i + 0], ux.f[i]);
            CHECK_F32("md_mm_unpack_xyz_ps y", i, rec[5 * i + 1], uy.f[i]);
            CHECK_F32("md_mm_unpack_xyz_ps z", i, rec[5 * i + 2], uz.f[i]);
            CHECK_F32("md_mm_unpack_xyzw_ps w", i, rec[5 * i + 3], uw.f[i]);
        }
        CHECK_F32("md_mm256_unpack_xyz_ps x", i, rec[5 * i + 0], Ux.f[i]);
        CHECK_F32("md_mm256_unpack_xyz_ps y", i, rec[5 * i + 1], Uy.f[i]);
        CHECK_F32("md_mm256_unpack_xyz_ps z", i, rec[5 * i + 2], Uz.f[i]);
        CHECK_F32("md_mm256_unpack_xyzw_ps w", i, rec[5 * i + 3], Uw.f[i]);
    }
    EXPECT_EQ(0, fails);
}

// THE NATIVE WIDTH: every md_xv operation against the scalar reference, over MD_XV_WIDTH lanes

UTEST(simd_ops, native_width) {
    rng = 12;
    int fails = 0;
    enum { W = MD_XV_WIDTH };
    for (int c = 0; c < 64; ++c) {
        ALIGNAS(32) float a[8], b[8], d[8], o[8];
        fill_f32(a, 8); fill_f32(b, 8); fill_f32(d, 8);
        const md_xv A = md_xv_loadu_ps(a), B = md_xv_loadu_ps(b), D = md_xv_loadu_ps(d);
#define XV_CHECK(name, expr, ref) do { md_xv_storeu_ps(o, expr); for (int i = 0; i < W; ++i) CHECK_F32(name, i, ref, o[i]); } while (0)
        XV_CHECK("md_xv_add_ps", md_xv_add_ps(A, B), ref_add(a[i], b[i]));
        XV_CHECK("md_xv_sub_ps", md_xv_sub_ps(A, B), ref_sub(a[i], b[i]));
        XV_CHECK("md_xv_mul_ps", md_xv_mul_ps(A, B), ref_mul(a[i], b[i]));
        md_xv_storeu_ps(o, md_xv_div_ps(A, B));
        for (int i = 0; i < W; ++i) CHECK_ULPS_F32("md_xv_div_ps", i, ref_div(a[i], b[i]), o[i], DIV_ULPS);
        XV_CHECK("md_xv_min_ps", md_xv_min_ps(A, B), ref_min(a[i], b[i]));
        XV_CHECK("md_xv_max_ps", md_xv_max_ps(A, B), ref_max(a[i], b[i]));
        XV_CHECK("md_xv_fmadd_ps",  md_xv_fmadd_ps(A, B, D),  ref_fmadd(a[i], b[i], d[i]));
        XV_CHECK("md_xv_fnmadd_ps", md_xv_fnmadd_ps(A, B, D), ref_fnmadd(a[i], b[i], d[i]));
        XV_CHECK("md_xv_fmsub_ps",  md_xv_fmsub_ps(A, B, D),  ref_fmsub(a[i], b[i], d[i]));
        XV_CHECK("md_xv_sqrt_ps",  md_xv_sqrt_ps(md_xv_abs_ps(A)), ref_sqrt(a[i]));
        XV_CHECK("md_xv_abs_ps",   md_xv_abs_ps(A),  ref_abs(a[i]));
        XV_CHECK("md_xv_neg_ps",   md_xv_neg_ps(A),  ref_neg(a[i]));
        XV_CHECK("md_xv_floor_ps", md_xv_floor_ps(md_xv_mul_ps(A, md_xv_set1_ps(4096.0f))), ref_floor(SPREAD(a[i])));
        XV_CHECK("md_xv_ceil_ps",  md_xv_ceil_ps(md_xv_mul_ps(A, md_xv_set1_ps(4096.0f))),  ref_ceil(SPREAD(a[i])));
        XV_CHECK("md_xv_round_ps", md_xv_round_ps(md_xv_mul_ps(A, md_xv_set1_ps(4096.0f))), ref_round(SPREAD(a[i])));
        XV_CHECK("md_xv_fract_ps", md_xv_fract_ps(md_xv_mul_ps(A, md_xv_set1_ps(4096.0f))), ref_fract(SPREAD(a[i])));
        XV_CHECK("md_xv_and_ps",    md_xv_and_ps(A, B),    ref_and(a[i], b[i]));
        XV_CHECK("md_xv_or_ps",     md_xv_or_ps(A, B),     ref_or(a[i], b[i]));
        XV_CHECK("md_xv_xor_ps",    md_xv_xor_ps(A, B),    ref_xor(a[i], b[i]));
        XV_CHECK("md_xv_andnot_ps", md_xv_andnot_ps(A, B), ref_andnot(a[i], b[i]));
        XV_CHECK("md_xv_cmpeq_ps",  md_xv_cmpeq_ps(A, A),  ref_cmpeq(a[i], a[i]));
        XV_CHECK("md_xv_cmpneq_ps", md_xv_cmpneq_ps(A, B), ref_cmpneq(a[i], b[i]));
        XV_CHECK("md_xv_cmplt_ps",  md_xv_cmplt_ps(A, B),  ref_cmplt(a[i], b[i]));
        XV_CHECK("md_xv_cmple_ps",  md_xv_cmple_ps(A, B),  ref_cmple(a[i], b[i]));
        XV_CHECK("md_xv_cmpgt_ps",  md_xv_cmpgt_ps(A, B),  ref_cmpgt(a[i], b[i]));
        XV_CHECK("md_xv_cmpge_ps",  md_xv_cmpge_ps(A, B),  ref_cmpge(a[i], b[i]));
        XV_CHECK("md_xv_blendv_ps", md_xv_blendv_ps(A, B, D), d[i] < 0 ? b[i] : a[i]);
        XV_CHECK("md_xv_select_ps", md_xv_select_ps(A, B, md_xv_cmplt_ps(D, md_xv_setzero_ps())), d[i] < 0 ? b[i] : a[i]);
        XV_CHECK("md_xv_set1_ps",   md_xv_set1_ps(a[0]), a[0]);
        XV_CHECK("md_xv_setzero_ps", md_xv_setzero_ps(), 0.0f);
        md_xv_store_ps(o, A);
        for (int i = 0; i < W; ++i) CHECK_F32("md_xv_store_ps", i, a[i], o[i]);
        XV_CHECK("md_xv_load_ps", md_xv_load_ps(a), a[i]);

        int mask = 0;
        float sum_ref = 0, mn = a[0], mx = a[0];
        for (int i = 0; i < W; ++i) { mask |= (a[i] < 0) << i; mn = MIN(mn, a[i]); mx = MAX(mx, a[i]); }
        CHECK_I64("md_xv_movemask_ps", 0, mask, md_xv_movemask_ps(A));
        CHECK_F32("md_xv_reduce_min_ps", 0, mn, md_xv_reduce_min_ps(A));
        CHECK_F32("md_xv_reduce_max_ps", 0, mx, md_xv_reduce_max_ps(A));
        // The sum in its order, ((x0 + x1) + (x2 + x3)), with the two halves added first for eight lanes. One
        // rounding per statement through volatile, which -ffast-math may not reassociate.
        {
            volatile float h[4], p0, p1;
            for (int i = 0; i < 4; ++i) h[i] = W == 8 ? a[i] + a[i + 4] : a[i];
            p0 = h[0] + h[1];
            p1 = h[2] + h[3];
            sum_ref = p0 + p1;
        }
        CHECK_F32("md_xv_reduce_add_ps", 0, sum_ref, md_xv_reduce_add_ps(A));

        md_xv_storeu_ps(o, md_xv_rsqrt_ps(md_xv_abs_ps(A)));
        for (int i = 0; i < W; ++i) CHECK_NEAR_F32("md_xv_rsqrt_ps", i, 1.0 / sqrt(fabs((double)a[i])), o[i], 1.5 / 4096.0);
        float xin[8];
        md_xv_storeu_ps(xin, md_xv_mul_ps(A, md_xv_set1_ps(0.1f)));
        md_xv_storeu_ps(o, md_xv_exp_ps(md_xv_loadu_ps(xin)));
        for (int i = 0; i < W; ++i) CHECK_NEAR_F32("md_xv_exp_ps", i, exp((double)xin[i]), o[i], 3e-7);
        md_xv S, C;
        md_xv_sincos_ps(A, &S, &C);
        md_xv_storeu_ps(o, S);
        for (int i = 0; i < W; ++i) if (fabs(o[i] - sin((double)a[i])) > 2e-7) { FAIL_F32("md_xv_sincos_ps sin", i, (float)sin(a[i]), o[i]); fails += 1; }
        md_xv_storeu_ps(o, md_xv_minimage_ps(A, md_xv_set1_ps(4.0f), md_xv_set1_ps(0.25f)));
        for (int i = 0; i < W; ++i) CHECK_F32("md_xv_minimage_ps", i, (float)((double)a[i] - 4.0 * rint((double)a[i] * 0.25)), o[i]);
        md_xv_storeu_ps(o, md_xv_cubic_spline_ps(A, B, D, A, md_xv_set1_ps(0.0f), md_xv_set1_ps(0.5f)));
        for (int i = 0; i < W; ++i) CHECK_F32("md_xv_cubic_spline_ps", i, b[i], o[i]);   // at t = 0, p1
    }
    // The packed triplets at the native width
    float xyz[24], back[24];
    for (int i = 0; i < 24; ++i) xyz[i] = (float)i;
    md_xv X, Y, Z;
    md_xv_load_xyz_packed_ps(&X, &Y, &Z, xyz);
    md_xv_store_xyz_packed_ps(back, X, Y, Z);
    for (int i = 0; i < 3 * W; ++i) CHECK_F32("md_xv_load/store_xyz_packed_ps", i, xyz[i], back[i]);
    EXPECT_EQ(0, fails);
}

// AVX-512, where it is enabled: against the 256 bit forms
#if defined(__AVX512F__) && defined(__AVX512DQ__)
UTEST(simd_ops, avx512) {
    rng = 13;
    int fails = 0;
    for (int c = 0; c < NUM_CASES; ++c) {
        float a[16], p[16], r[16];
        fill_f32(a, 16);
        for (int i = 0; i < 16; ++i) { p[i] = (i % 4 == 3) ? 0.0f : 10.0f; r[i] = 0.5f * a[(i + 3) & 15]; }
        const __m512 A = _mm512_loadu_ps(a);
        const md_256 lo = md_mm256_loadu_ps(a), hi = md_mm256_loadu_ps(a + 8);
        float o[16], e[16];
        // Built from two halves
        _mm512_storeu_ps(o, md_mm512_from_m256_ps(lo, hi));
        for (int i = 0; i < 16; ++i) CHECK_F32("md_mm512_from_m256_ps", i, a[i], o[i]);
        // exp and sincos agree with the 256 bit forms bit for bit
        _mm512_storeu_ps(o, md_mm512_exp_ps(A));
        md_mm256_storeu_ps(e, md_mm256_exp_ps(lo)); md_mm256_storeu_ps(e + 8, md_mm256_exp_ps(hi));
        for (int i = 0; i < 16; ++i) CHECK_F32("md_mm512_exp_ps", i, e[i], o[i]);
        __m512 S, C;
        md_256 s0, c0, s1, c1;
        md_mm512_sincos_ps(A, &S, &C);
        md_mm256_sincos_ps(lo, &s0, &c0);
        md_mm256_sincos_ps(hi, &s1, &c1);
        _mm512_storeu_ps(o, S); md_mm256_storeu_ps(e, s0); md_mm256_storeu_ps(e + 8, s1);
        for (int i = 0; i < 16; ++i) CHECK_F32("md_mm512_sincos_ps sin", i, e[i], o[i]);
        _mm512_storeu_ps(o, C); md_mm256_storeu_ps(e, c0); md_mm256_storeu_ps(e + 8, c1);
        for (int i = 0; i < 16; ++i) CHECK_F32("md_mm512_sincos_ps cos", i, e[i], o[i]);
        // deperiodize
        const __m512 P = _mm512_loadu_ps(p), R = _mm512_loadu_ps(r);
        __m512 rp = _mm512_div_ps(_mm512_set1_ps(1.0f), _mm512_max_ps(P, _mm512_set1_ps(1.0f)));
        md_mm256_storeu_ps(e, md_mm256_deperiodize_ps(lo, md_mm256_loadu_ps(r), md_mm256_loadu_ps(p)));
        md_mm256_storeu_ps(e + 8, md_mm256_deperiodize_ps(hi, md_mm256_loadu_ps(r + 8), md_mm256_loadu_ps(p + 8)));
        for (int k = 0; k < 2; ++k) {
            _mm512_storeu_ps(o, k ? md_mm512_deperiodize2_ps(A, R, P, rp) : md_mm512_deperiodize_ps(A, R, P));
            for (int i = 0; i < 16; ++i) {
                if (p[i] == 0.0f) CHECK_F32("md_mm512_deperiodize(2)_ps", i, a[i], o[i]);
                else if (fabs(o[i] - e[i]) > 4e-7 * (fabs(a[i]) + fabs(r[i]) + p[i]) + 1e-6) { FAIL_F32("md_mm512_deperiodize(2)_ps", i, e[i], o[i]); fails += 1; }
            }
        }
        double ad[8], pd[8], rd[8], od[8];
        for (int i = 0; i < 8; ++i) { ad[i] = a[i]; pd[i] = p[i]; rd[i] = r[i]; }
        _mm512_storeu_pd(od, md_mm512_deperiodize_pd(_mm512_loadu_pd(ad), _mm512_loadu_pd(rd), _mm512_loadu_pd(pd)));
        u256d e0 = { .v = md_mm256_deperiodize_pd(md_mm256_loadu_pd(ad), md_mm256_loadu_pd(rd), md_mm256_loadu_pd(pd)) };
        for (int i = 0; i < 4; ++i) CHECK_NEAR_F64("md_mm512_deperiodize_pd", i, e0.d[i], od[i], 1e-12);
        _mm512_storeu_pd(od, md_mm512_deperiodize2_pd(_mm512_loadu_pd(ad), _mm512_loadu_pd(rd), _mm512_loadu_pd(pd), _mm512_set1_pd(0.1)));
        for (int i = 0; i < 8; ++i) if (pd[i] == 0.0) CHECK_F64("md_mm512_deperiodize2_pd", i, ad[i], od[i]);
        // reductions
        double sum = 0, scale = 0;
        float mn = a[0], mx = a[0];
        for (int i = 0; i < 16; ++i) { mn = MIN(mn, a[i]); mx = MAX(mx, a[i]); sum += a[i]; scale += fabs(a[i]); }
        if (fabs(md_mm512_reduce_add_ps(A) - sum) > 1e-6 * scale) { FAIL_F32("md_mm512_reduce_add_ps", 0, (float)sum, md_mm512_reduce_add_ps(A)); fails += 1; }
        CHECK_F32("md_mm512_reduce_min_ps", 0, mn, md_mm512_reduce_min_ps(A));
        CHECK_F32("md_mm512_reduce_max_ps", 0, mx, md_mm512_reduce_max_ps(A));
    }
    // packed triplets and their gather
    float xyz[60];
    for (int i = 0; i < 60; ++i) xyz[i] = (float)i;
    __m512 X, Y, Z;
    float ox[16], oz[16];
    md_mm512_load_xyz_packed_ps(&X, &Y, &Z, xyz);
    _mm512_storeu_ps(ox, X); _mm512_storeu_ps(oz, Z);
    for (int i = 0; i < 16; ++i) { CHECK_F32("md_mm512_load_xyz_packed_ps x", i, (float)(3 * i), ox[i]); CHECK_F32("md_mm512_load_xyz_packed_ps z", i, (float)(3 * i + 2), oz[i]); }
    int32_t idx[16];
    for (int i = 0; i < 16; ++i) idx[i] = (i * 7) % 20;
    md_mm512_i32gather_xyz_ps(&X, &Y, &Z, xyz, _mm512_loadu_si512(idx));
    _mm512_storeu_ps(ox, Y);
    for (int i = 0; i < 16; ++i) CHECK_F32("md_mm512_i32gather_xyz_ps y", i, (float)(3 * idx[i] + 1), ox[i]);
    EXPECT_EQ(0, fails);
}
#endif
