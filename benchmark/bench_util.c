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
#endif
static double ref_com_trig(const float* in_x, const float* in_w, int64_t count, float x_max) {
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
    if (r2 > 1.0e-15) {
        theta_prim += atan2(-y, -x);
    }

    return (theta_prim / TWO_PI) * x_max;
}
#if defined(_MSC_VER)
#pragma float_control(precise, on, pop)
#endif

UBENCH_EX(util, com_speed_of_light) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));

    SETUP_XYZW

    md_256 vx = md_mm256_setzero_ps();
    md_256 vy = md_mm256_setzero_ps();
    md_256 vz = md_mm256_setzero_ps();
    md_256 vw = md_mm256_setzero_ps();

    UBENCH_DO_BENCHMARK() {
        for (int64_t i = 0; i < size; i += 8) {
            vx = md_mm256_add_ps(vx, md_mm256_load_ps(x + i));
            vy = md_mm256_add_ps(vy, md_mm256_load_ps(y + i));
            vz = md_mm256_add_ps(vz, md_mm256_load_ps(z + i));
            vw = md_mm256_add_ps(vw, md_mm256_load_ps(w + i));
        }
    }

    printf("com: %.3f %.3f %.3f %.3f\n", md_mm256_cvtss_f32(vx), md_mm256_cvtss_f32(vy), md_mm256_cvtss_f32(vz), md_mm256_cvtss_f32(vw));

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_ref) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));

    SETUP_XYZW

    vec3_t com = {0};
    UBENCH_DO_BENCHMARK() {
        com.x = (float)ref_com_trig(x, w, size, pbc_ext.x);
        com.y = (float)ref_com_trig(y, w, size, pbc_ext.y);
        com.z = (float)ref_com_trig(z, w, size, pbc_ext.z);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_soa) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));

    SETUP_XYZW

    vec3_t com = {0};
    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_ortho(x, y, z, w, 0, size, pbc_ext);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_soa_indices) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    
    SETUP_XYZW
    SETUP_INDICES

    vec3_t com = {0};
    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_ortho(x, y, z, w, indices, size, pbc_ext);        
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_vec4) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));

    SETUP_XYZW_VEC4

    vec3_t com = {0};
    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_vec4_ortho(xyzw, 0, size, pbc_ext);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_vec4_indices) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));

    SETUP_XYZW_VEC4
    SETUP_INDICES

    vec3_t com = {0};
    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_vec4_ortho(xyzw, indices, size, pbc_ext);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}