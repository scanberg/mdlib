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
    md_128 vec;
	float val[4];
} v4_t;

typedef union {
    md_256 vec;
    float val[8];
} v8_t;

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1985, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

static float FOPI = 1.27323954473516;

static float DP1 = 0.78515625;
static float DP2 = 2.4187564849853515625e-4;
static float DP3 = 3.77489497744594108e-8;

static float sincof[] = {
    -1.9515295891E-4,
    8.3321608736E-3,
    -1.6666654611E-1
};
static float coscof[] = {
    2.443315711809948E-005,
    -1.388731625493765E-003,
    4.166664568298827E-002
};

float cephes_sinf( float xx )
{
    float *p;
    float x, y, z;
    register unsigned long j;
    register int sign;

    if (xx == 2.0f) {
        while(0) {};
    }

    sign = 1;
    x = xx;
    if( xx < 0 )
    {
        sign = -1;
        x = -xx;
    }
    j = FOPI * x; /* integer part of x/(PI/4) */
    y = j;
    /* map zeros to origin */
    if( j & 1 )
    {
        j += 1;
        y += 1.0;
    }
    j &= 7; /* octant modulo 360 degrees */
    /* reflect in x axis */
    if( j > 3)
    {
        sign = -sign;
        j -= 4;
    }

    /* Extended precision modular arithmetic */
    x = ((x - y * DP1) - y * DP2) - y * DP3;

    z = x * x;
    if( (j==1) || (j==2) )
    {
        p = coscof;
        y = *p++;
        y = y * z + *p++;
        y = y * z + *p++;
        y *= z * z;
        y -= 0.5 * z;
        y += 1.0;
    }
    else
    {
        p = sincof;
        y = *p++;
        y = y * z + *p++;
        y = y * z + *p++;
        y *= z * x;
        y += x;
    }
    if(sign < 0)
        y = -y;
    return( y);
}

// https://marc-b-reynolds.github.io/math/2019/04/24/ULPDiff.html

// get the bit pattern of 'x'
static uint32_t f32_to_bits(float x)
{
    uint32_t u; memcpy(&u, &x, 4); return u;
}

// okay...the name's goofy, but you get what I mean
static uint32_t u32_abs(uint32_t x)
{
    return (int32_t)x >= 0 ? x : -x;
}

static uint32_t f32_ulp_dist(float a, float b)
{
    uint32_t ua = f32_to_bits(a);
    uint32_t ub = f32_to_bits(b);
    uint32_t s  = ub^ua;

    if ((int32_t)s >= 0)
        return u32_abs(ua-ub);

    return ua+ub+0x80000000;
}

#if defined(_MSC_VER)
#pragma float_control(precise, on, push)
#endif

UTEST(simd, sin_cos) {

    const float values[] = {
        0.00100000005f,
        0.00999999978f,
        0.100000001f,
        1.0f,
        2.0f,
        3.0f,
        4.0f,
        10.f,
        100.f,
        1000.f,
		8192.f,
        8193.f,
        10001.f,
        10002.f,
        100001.f,
        1000000.f,
        10000000.f,
        100000000.f,
    };

    // Reference sin values computed by wolfram alpha
    const double reference_sin[] = {
        0.0009999998833333166666673015874223985902802826153752216985433025,
        0.0099998331141776645906340811925307282622351416106754492775577823,
        0.0998334176418323175349232560169359110982006816705633483355436365,
        0.8414709848078965066525023216302989996225630607983710656727517099,
        0.9092974268256816953960198659117448427022549714478902683789730115,
        0.1411200080598672221007448028081102798469332642522655841518826412,
        -0.756802495307928251372639094511829094135912887336472571485416773,
        -0.544021110889369813404747661851377281683643012916223891574184012,
        -0.506365641109758793656557610459785432065032721290657323443392473,
        0.8268795405320025602558874291092181412127249678477883209081232758,
        -0.956173152843146286228747525773835843477602162161292598199027227,
        -0.270238329228266638196168321085249325579720412657590400914146482,
        -0.966335274441843598152339541522975821062664960879442783733899435,
        -0.738611965157047291838492797313355681872029717023198046517078829,
        -0.821617964837152459458049701499136588670334311523615525610326265,
        -0.349993502171292952117652486780771469061406605328716273857059054,
        0.4205477931907824912985065897409455951671752475308045898687660412,
        0.9316390271097260080275166536120429704729018385275364343082838951,
    };
    
    for (int i = 0; i < (int)ARRAY_SIZE(values); ++i) {
        const float x = values[i];
        const float ref_sin = sin(x);
        const float ref_cos = cos(x);
        printf("Testing for value: %f\n", x);

        {
            float sin = sinf(x);
            const uint32_t ulp_diff_libc  = f32_ulp_dist(ref_sin, sin);
            const double diff_libc        = fabs((double)ref_sin - (double)sin);
            const uint32_t ulp_diff_ref   = f32_ulp_dist(sin, (float)reference_sin[i]);
            const double diff_ref         = fabs((double)sin - reference_sin[i]);
            printf("libc    ulp diff libc, f32 diff libc, ulp diff ref, f32 diff ref: %10d %10.10f %10d %10.10f\n", ulp_diff_libc, diff_libc, ulp_diff_ref, diff_ref);
        }

        {
            float sin = cephes_sinf(x);
            const uint32_t ulp_diff_libc  = f32_ulp_dist(ref_sin, sin);
            const double diff_libc        = fabs((double)ref_sin - (double)sin);
            const uint32_t ulp_diff_ref   = f32_ulp_dist(sin, (float)reference_sin[i]);
            const double diff_ref         = fabs((double)sin - reference_sin[i]);
            printf("cephes  ulp diff libc, f32 diff libc, ulp diff ref, f32 diff ref: %10d %10.10f %10d %10.10f\n", ulp_diff_libc, diff_libc, ulp_diff_ref, diff_ref);
        }

        {
            v4_t vsin, vcos;
            md_mm_sincos_ps(md_mm_set1_ps(x), &vsin.vec, &vcos.vec);
            float sin = vsin.val[0];
            const uint32_t ulp_diff_libc  = f32_ulp_dist(ref_sin, sin);
            const double diff_libc        = fabs((double)ref_sin - (double)sin);
            const uint32_t ulp_diff_ref   = f32_ulp_dist(sin, (float)reference_sin[i]);
            const double diff_ref         = fabs((double)sin - reference_sin[i]);
            printf("sse     ulp diff libc, f32 diff libc, ulp diff ref, f32 diff ref: %10d %10.10f %10d %10.10f\n", ulp_diff_libc, diff_libc, ulp_diff_ref, diff_ref);
        }
    }
}

#if defined(_MSC_VER)
#pragma float_control(pop)
#endif

UTEST(simd, hsum) {
    {
        md_128 x = md_mm_set_ps(1.0f, 2.0f, 3.0f, 4.0f);
        float sum = md_mm_reduce_add_ps(x);
        EXPECT_NEAR(sum, 10.0f, 0.0001f);
    }
    {
        md_256 x = md_mm256_set_ps(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
        float sum = md_mm256_reduce_add_ps(x);
        EXPECT_NEAR(sum, 36.0f, 0.0001f);
    }
}

UTEST(simd, hmax) {
    {
        md_128 x = md_mm_set_ps(1.0f, 2.0f, 3.0f, 4.0f);
        float max = md_mm_reduce_max_ps(x);
        EXPECT_NEAR(max, 4.0f, 0.0001f);
    }
    {
        md_256 x = md_mm256_set_ps(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
        float max = md_mm256_reduce_max_ps(x);
        EXPECT_NEAR(max, 8.0f, 0.0001f);
    }
}

UTEST(simd, hmin) {
    {
        md_128 x = md_mm_set_ps(1.0f, 2.0f, -3.0f, 4.0f);
        float min = md_mm_reduce_min_ps(x);
        EXPECT_NEAR(min, -3.0f, 0.0001f);
    }

    {
        md_256 x = md_mm256_set_ps(1.0f, 2.0f, 3.0f, -4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
        float min = md_mm256_reduce_min_ps(x);
        EXPECT_NEAR(min, -4.0f, 0.0001f);
    }
}

