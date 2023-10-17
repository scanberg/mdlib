#include "utest.h"

#ifdef __AVX__
// Reference implementation
#include "avx_mathfun.h"
#endif

#include <core/md_common.h>
#include <core/md_simd.h>
#include <math.h>

#include <stdbool.h>

typedef union {
	md_simd_f32_t vec;
	float val[md_simd_width_f32];
} vec_t;

static inline bool validate_eqf(md_simd_f32_t x, float ref[]) {
	vec_t v  = { .vec = x };
	for (int i = 0; i < md_simd_width_f32; ++i) {
		if (v.val[i] != ref[i]) return false;
	}
	return true;
}

static inline bool validate_close1f(md_simd_f32_t x, float ref, float eps) {
	vec_t v  = { .vec = x };
	for (int i = 0; i < md_simd_width_f32; ++i) {
		if (fabsf(v.val[i] - ref) > eps) return false;
	}
	return true;
}

static inline bool validate_closef(md_simd_f32_t x, float ref[], float eps) {
	vec_t v  = { .vec = x };
	for (int i = 0; i < md_simd_width_f32; ++i) {
		if (fabsf(v.val[i] - ref[i]) > eps) return false;
	}
	return true;
}

UTEST(simd, sin_cos) {
    md_f32x8_t sin, cos;
    md_simd_sincos_f32x8(md_simd_set1_f32(1.0f), &sin, &cos);
    EXPECT_TRUE(validate_close1f(sin, sinf(1.0f), 0.0001f));
    EXPECT_TRUE(validate_close1f(cos, cosf(1.0f), 0.0001f));
}

UTEST(simd, hsum) {
#if md_simd_width_f32 >= 4
    {
    md_f32x4_t x = md_simd_set_f32x4(1.0f, 2.0f, 3.0f, 4.0f);
    float sum = md_simd_hsum(x);
    EXPECT_NEAR(sum, 10.0f, 0.0001f);
    }
#endif
#if md_simd_width_f32 >= 8
    {
    md_f32x8_t x = md_simd_set_f32x8(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
    float sum = md_simd_hsum(x);
    EXPECT_NEAR(sum, 36.0f, 0.0001f);
    }
#endif
}

UTEST(simd, hmax) {
#if md_simd_width_f32 >= 4
    {
    md_f32x4_t x = md_simd_set_f32x4(1.0f, 2.0f, 3.0f, 4.0f);
    float max = md_simd_hmax(x);
    EXPECT_NEAR(max, 4.0f, 0.0001f);
    }
#endif
#if md_simd_width_f32 >= 8
    {
    md_f32x8_t x = md_simd_set_f32x8(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
    float max = md_simd_hmax(x);
    EXPECT_NEAR(max, 8.0f, 0.0001f);
    }
#endif
}

UTEST(simd, hmin) {
#if md_simd_width_f32 >= 4
    {
    md_f32x4_t x = md_simd_set_f32x4(1.0f, 2.0f, -3.0f, 4.0f);
    float min = md_simd_hmin(x);
    EXPECT_NEAR(min, -3.0f, 0.0001f);
    }
#endif
#if md_simd_width_f32 >= 8
    {
    md_f32x8_t x = md_simd_set_f32x8(1.0f, 2.0f, 3.0f, -4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
    float min = md_simd_hmin(x);
    EXPECT_NEAR(min, -4.0f, 0.0001f);
    }
#endif
}

