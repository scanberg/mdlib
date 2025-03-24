#include "ubench.h"

#include <core/md_os.h>
#include <core/md_arena_allocator.h>
#include <md_gro.h>
#include <md_util.h>

#include <stdlib.h>

static const vec3_t pbc_ext = {100, 100, 100};
static const int size = 1024 * 1024;

#define SETUP_XYZW \
	float* x = md_alloc(alloc, size * sizeof(float)); \
	float* y = md_alloc(alloc, size * sizeof(float)); \
	float* z = md_alloc(alloc, size * sizeof(float)); \
	float* w = md_alloc(alloc, size * sizeof(float)); \
	srand(0); \
	for (int i = 0; i < size; ++i) { \
		x[i] = (float)rand() / RAND_MAX * pbc_ext.x; \
		y[i] = (float)rand() / RAND_MAX * pbc_ext.y; \
		z[i] = (float)rand() / RAND_MAX * pbc_ext.z; \
		w[i] = (float)rand() / RAND_MAX * 1.0f + 1.0f; \
	}

#define SETUP_XYZW_VEC4 \
	vec4_t* xyzw = md_alloc(alloc, size * sizeof(vec4_t)); \
	srand(0); \
	for (int i = 0; i < size; ++i) { \
		xyzw[i] = vec4_set( \
			(float)rand() / RAND_MAX * pbc_ext.x, \
			(float)rand() / RAND_MAX * pbc_ext.y, \
			(float)rand() / RAND_MAX * pbc_ext.z, \
			(float)rand() / RAND_MAX * 1.0f + 1.0f); \
	}

#define SETUP_INDICES \
	int* indices = md_alloc(alloc, size * sizeof(int)); \
	for (int i = 0; i < size; ++i) { \
		indices[i] = i; \
	}

#if defined(_MSC_VER)
#pragma float_control(precise, on, push)
#elif defined (MD_COMPILER_GCC)
#pragma GCC push_options
#pragma GCC optimize ("no-fast-math")
#endif
static double com_trig_ref(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const double scl = TWO_PI / x_max;

    double acc_c = 0.0;
    double acc_s = 0.0;
    double acc_w = 0.0;
    for (int64_t i = 0; i < count; ++i) {
        double theta = in_x[i] * scl;
        double w = in_w[i];
        acc_c += w * cos(theta);
        acc_s += w * sin(theta);
        acc_w += w;
    }

    const double y = acc_s / acc_w;
    const double x = acc_c / acc_w;
    const double r2 = x*x + y*y;

    double theta_prim = PI;
    if (r2 > 1.0e-8) {
        theta_prim += atan2(-y, -x);
    }

    return (theta_prim / TWO_PI) * x_max;
}

static double com_min_image_ref(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const double ext = x_max;
    const double r_ext = 1.0 / ext;
    double acc_x = in_x[0];
    double acc_w = in_w[0];
    
    for (int64_t i = 1; i < count; ++i) {
        double r = acc_x / acc_w;
        double x = in_x[i];
        double w = in_w[i];
        x = deperiodize_ortho2(x, r, ext, r_ext);
        acc_x += x * w;
        acc_w += w;
    }

    return (acc_x / acc_w);
}
#if defined(_MSC_VER)
#pragma float_control(pop)
#elif defined (MD_COMPILER_GCC)
#pragma GCC pop_options
#endif

static double com_min_image_sse(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const md_128 v_ext = md_mm_set1_ps(x_max);
    const md_128 v_r_ext = md_mm_set1_ps(1.0f / x_max);
    md_128 v_acc_x = md_mm_loadu_ps(in_x);
    md_128 v_acc_w = md_mm_loadu_ps(in_w);
    for (int64_t i = 4; i < count; i += 4) {
        md_128 r = md_mm_div_ps(v_acc_x, v_acc_w);
        md_128 x = md_mm_loadu_ps(in_x + i);
        md_128 w = md_mm_loadu_ps(in_w + i);
        x = md_mm_deperiodize2_ps(x, r, v_ext, v_r_ext);
        v_acc_x = md_mm_add_ps(v_acc_x, md_mm_mul_ps(x, w));
        v_acc_w = md_mm_add_ps(v_acc_w, w);
    }

    double acc_x = md_mm_reduce_add_ps(v_acc_x);
    double acc_w = md_mm_reduce_add_ps(v_acc_w);
    return (acc_x / acc_w);
}

static double com_min_image_avx2(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const md_256 v_ext = md_mm256_set1_ps(x_max);
    const md_256 v_r_ext = md_mm256_set1_ps(1.0f / x_max);
    md_256 v_acc_x = md_mm256_loadu_ps(in_x);
    md_256 v_acc_w = md_mm256_loadu_ps(in_w);
    for (int64_t i = 8; i < count; i += 8) {
        md_256 x = md_mm256_loadu_ps(in_x + i);
        md_256 w = md_mm256_loadu_ps(in_w + i);
        md_256 r = md_mm256_div_ps(v_acc_x, v_acc_w);
        x = md_mm256_deperiodize2_ps(x, r, v_ext, v_r_ext);
        v_acc_x = md_mm256_add_ps(v_acc_x, md_mm256_mul_ps(x, w));
        v_acc_w = md_mm256_add_ps(v_acc_w, w);
    }

    double acc_x = md_mm256_reduce_add_ps(v_acc_x);
    double acc_w = md_mm256_reduce_add_ps(v_acc_w);
    return (acc_x / acc_w);
}

#if defined(__AVX512F__)
static double com_min_image_avx512(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const md_512 v_ext = _mm512_set1_ps(x_max);
    const md_512 v_r_ext = _mm512_set1_ps(1.0f / x_max);
    md_512 v_acc_x = _mm512_loadu_ps(in_x);
    md_512 v_acc_w = _mm512_loadu_ps(in_w);
    for (int64_t i = 16; i < count; i += 16) {
        md_512 r = _mm512_div_ps(v_acc_x, v_acc_w);
        md_512 x = _mm512_loadu_ps(in_x + i);
        md_512 w = _mm512_loadu_ps(in_w + i);
        x = md_mm512_deperiodize2_ps(x, r, v_ext, v_r_ext);
        v_acc_x = _mm512_add_ps(v_acc_x, _mm512_mul_ps(x, w));
        v_acc_w = _mm512_add_ps(v_acc_w, w);
    }

    double acc_x = _mm512_reduce_add_ps(v_acc_x);
    double acc_w = _mm512_reduce_add_ps(v_acc_w);
    return (acc_x / acc_w);
}
#endif

static double com_trig_sse(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const md_128 v_scl = md_mm_set1_ps((float)(TWO_PI / x_max));
    md_128 v_acc_c = md_mm_setzero_ps();
    md_128 v_acc_s = md_mm_setzero_ps();
    md_128 v_acc_w = md_mm_setzero_ps();

    for (int64_t i = 0; i < count; i += 4) {
        md_128 v_x = md_mm_loadu_ps(in_x + i);
        md_128 v_w = md_mm_loadu_ps(in_w + i);
        md_128 v_theta = md_mm_mul_ps(v_x, v_scl);
        md_128 v_c, v_s;
        md_mm_sincos_ps(v_theta, &v_s, &v_c);
        v_acc_s = md_mm_add_ps(v_acc_s, md_mm_mul_ps(v_s, v_w));
        v_acc_c = md_mm_add_ps(v_acc_c, md_mm_mul_ps(v_c, v_w));
        v_acc_w = md_mm_add_ps(v_acc_w, v_w);
    }

    double acc_s = md_mm_reduce_add_ps(v_acc_s);
    double acc_c = md_mm_reduce_add_ps(v_acc_c);
    double acc_w = md_mm_reduce_add_ps(v_acc_w);

    const double y = acc_s / acc_w;
    const double x = acc_c / acc_w;
    const double r2 = x*x + y*y;

    double theta_prim = PI;
    if (r2 > 1.0e-8) {
        theta_prim += atan2(-y, -x);
    }

    return (theta_prim / TWO_PI) * x_max;
}

static double com_trig_avx2(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const md_256 v_scl = md_mm256_set1_ps((float)(TWO_PI / x_max));
    md_256 v_acc_c = md_mm256_setzero_ps();
    md_256 v_acc_s = md_mm256_setzero_ps();
    md_256 v_acc_w = md_mm256_setzero_ps();

    for (int64_t i = 0; i < count; i += 8) {
        md_256 v_x = md_mm256_loadu_ps(in_x + i);
        md_256 v_w = md_mm256_loadu_ps(in_w + i);
        md_256 v_theta = md_mm256_mul_ps(v_x, v_scl);
        md_256 v_c, v_s;
        md_mm256_sincos_ps(v_theta, &v_s, &v_c);
        v_acc_s = md_mm256_add_ps(v_acc_s, md_mm256_mul_ps(v_s, v_w));
        v_acc_c = md_mm256_add_ps(v_acc_c, md_mm256_mul_ps(v_c, v_w));
        v_acc_w = md_mm256_add_ps(v_acc_w, v_w);
    }

    double acc_s = md_mm256_reduce_add_ps(v_acc_s);
    double acc_c = md_mm256_reduce_add_ps(v_acc_c);
    double acc_w = md_mm256_reduce_add_ps(v_acc_w);

    const double y = acc_s / acc_w;
    const double x = acc_c / acc_w;
    const double r2 = x*x + y*y;

    double theta_prim = PI;
    if (r2 > 1.0e-8) {
        theta_prim += atan2(-y, -x);
    }

    return (theta_prim / TWO_PI) * x_max;
}

#if defined(__AVX512F__)
static double com_trig_avx512(const float* in_x, const float* in_w, int64_t count, float x_max) {
    const double scl = TWO_PI / x_max;
    const md_512 v_scl = _mm512_set1_ps((float)(TWO_PI / x_max));
    md_512 v_acc_c = _mm512_setzero_ps();
    md_512 v_acc_s = _mm512_setzero_ps();
    md_512 v_acc_w = _mm512_setzero_ps();

    for (int64_t i = 0; i < count; i += 16) {
        md_512 v_x = _mm512_loadu_ps(in_x + i);
        md_512 v_w = _mm512_loadu_ps(in_w + i);
        md_512 v_theta = _mm512_mul_ps(v_x, v_scl);
        md_512 v_c, v_s;
        md_mm512_sincos_ps(v_theta, &v_s, &v_c);
        v_acc_s = _mm512_add_ps(v_acc_s, _mm512_mul_ps(v_s, v_w));
        v_acc_c = _mm512_add_ps(v_acc_c, _mm512_mul_ps(v_c, v_w));
        v_acc_w = _mm512_add_ps(v_acc_w, v_w);
    }

    double acc_s = _mm512_reduce_add_ps(v_acc_s);
    double acc_c = _mm512_reduce_add_ps(v_acc_c);
    double acc_w = _mm512_reduce_add_ps(v_acc_w);

    const double y = acc_s / acc_w;
    const double x = acc_c / acc_w;
    const double r2 = x*x + y*y;

    double theta_prim = PI;
    if (r2 > 1.0e-8) {
        theta_prim += atan2(-y, -x);
    }

    return (theta_prim / TWO_PI) * x_max;
}
#endif

static double touch_mem(const float* in_x, const float* in_w, int64_t count, float x_len) {
    (void)x_len;
#if defined(__AVX512F__)
    md_512 v_acc_x = _mm512_setzero_ps();
    md_512 v_acc_w = _mm512_setzero_ps();
    for (int64_t i = 0; i < count; i += 16) {
        md_512 x = _mm512_loadu_ps(in_x + i);
        md_512 w = _mm512_loadu_ps(in_w + i);
        v_acc_x = _mm512_add_ps(v_acc_x, x);
        v_acc_w = _mm512_add_ps(v_acc_w, w);
}
    return _mm512_cvtss_f32(v_acc_x) + _mm512_cvtss_f32(v_acc_w);
#elif defined(__AVX2__)
    md_256 v_acc_x = md_mm256_setzero_ps();
    md_256 v_acc_w = md_mm256_setzero_ps();
    for (int64_t i = 0; i < count; i += 8) {
        md_256 x = md_mm256_loadu_ps(in_x + i);
        md_256 w = md_mm256_loadu_ps(in_w + i);
        v_acc_x = md_mm256_add_ps(v_acc_x, x);
        v_acc_w = md_mm256_add_ps(v_acc_w, w);
    }
    return md_mm256_cvtss_f32(v_acc_x) + md_mm256_cvtss_f32(v_acc_w);
#elif defined(__SSE2__)
    md_128 v_acc_x = md_mm_setzero_ps();
    md_128 v_acc_w = md_mm_setzero_ps();
    for (int64_t i = 0; i < count; i += 4) {
        md_128 x = md_mm_loadu_ps(in_x + i);
        md_128 w = md_mm_loadu_ps(in_w + i);
	    v_acc_x = md_mm_add_ps(v_acc_x, x);
	    v_acc_w = md_mm_add_ps(v_acc_w, w);
    }
    return md_mm_cvtss_f32(v_acc_x) + md_mm_cvtss_f32(v_acc_w);
#endif
}

typedef double (*com_func)(const float* in_x, const float* in_w, int64_t count, float x_max);

typedef struct com_func_entry_t {
    com_func    func;
	const char* name;
} com_func_entry_t;

void bench_com(void) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    SETUP_XYZW

    com_func_entry_t tests[] = {
        { touch_mem,            "Speed of light"},
        { com_min_image_ref,    "com_min_image_ref" },
        { com_trig_ref,         "com_trig_ref" },
#if defined (__SSE2__)
        { com_min_image_sse,    "com_min_image_sse" },
        { com_trig_sse,         "com_trig_sse" },
#endif
#if defined (__AVX2__)
        { com_min_image_avx2,   "com_min_image_avx2" },
        { com_trig_avx2,        "com_trig_avx2" },
#endif
#if defined (__AVX512F__)
        { com_min_image_avx512, "com_min_image_avx512" },
        { com_trig_avx512,      "com_trig_avx512" },
#endif
    };

    const int64_t num_iter = 100;

    printf("%-20s%25s%25s%25s\n", "name", "min iter time (ms)", "min elem time (ns)", "throughput (MB/s)");
    for (int i = 0; i < (int)ARRAY_SIZE(tests); ++i) {
        int64_t min_time = INT64_MAX;
        double acc = 0.0;
        for (int j = 0; j < num_iter; ++j) {
            int64_t start = ubench_ns();
            double res_x = tests[i].func(x, w, size, pbc_ext.x);
            double res_y = tests[i].func(y, w, size, pbc_ext.y);
            double res_z = tests[i].func(z, w, size, pbc_ext.z);
            int64_t end = ubench_ns();
            acc += res_x + res_y + res_z;
            min_time = MIN(min_time, end - start);
        }
        double GBps = (double)(size * sizeof(float) * 4) / (double)min_time;
        double MBps = GBps * 1024.0;
        double ms = (double)min_time / 1000000.0;
        double ns = (double)min_time / (double)(size * 3);
        printf("%-20s%25.2f%25.2f%25.2f\n", tests[i].name, ms, ns, MBps);
    }

    md_arena_allocator_destroy(alloc);
}
