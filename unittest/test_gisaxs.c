#include "utest.h"

#include <md_gisaxs.h>
#include <core/md_fft.h>
#include <core/md_allocator.h>
#include <core/md_common.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static float frand(uint32_t* state) {
    // xorshift32
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return (float)(x & 0xFFFFFF) / (float)0x1000000;
}

UTEST(gisaxs, fft_2d_matches_dft) {
    const int nx = 64;
    const int ny = 48;
    ASSERT_EQ(md_fft_valid_size(nx, true), nx);

    md_fft_2d_t* fft = md_fft_2d_create(nx, ny);
    ASSERT_TRUE(fft != NULL);

    const int nkx = nx / 2 + 1;
    float* in  = (float*)md_fft_alloc(sizeof(float) * nx * ny);
    float* out = (float*)md_fft_alloc(sizeof(float) * nkx * ny * 2);
    float* scratch = (float*)md_fft_alloc(sizeof(float) * md_fft_2d_scratch_size(fft));

    uint32_t seed = 1234567;
    for (int i = 0; i < nx * ny; ++i) in[i] = frand(&seed) - 0.5f;

    md_fft_2d_r2c(fft, in, out, 0, scratch);

    double max_err = 0.0;
    for (int ky = 0; ky < ny; ++ky) {
        for (int kx = 0; kx < nkx; ++kx) {
            double re = 0.0, im = 0.0;
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    const double ph = -2.0 * PI * ((double)kx * x / nx + (double)ky * y / ny);
                    re += in[y * nx + x] * cos(ph);
                    im += in[y * nx + x] * sin(ph);
                }
            }
            const float* v = out + 2 * (ky * nkx + kx);
            max_err = MAX(max_err, fabs(re - v[0]));
            max_err = MAX(max_err, fabs(im - v[1]));
        }
    }
    EXPECT_LT(max_err, 1.0e-3);

    md_fft_free(in);
    md_fft_free(out);
    md_fft_free(scratch);
    md_fft_2d_destroy(fft);
}

typedef struct test_system_t {
    size_t count;
    float* x;
    float* y;
    float* z;
    float* w;
    float* s;
} test_system_t;

static test_system_t make_system(size_t count, double bx, double by, double z0, double z1, bool two_classes) {
    test_system_t sys = {0};
    sys.count = count;
    sys.x = (float*)malloc(sizeof(float) * count);
    sys.y = (float*)malloc(sizeof(float) * count);
    sys.z = (float*)malloc(sizeof(float) * count);
    sys.w = (float*)malloc(sizeof(float) * count);
    sys.s = (float*)malloc(sizeof(float) * count);
    uint32_t seed = 42;
    for (size_t i = 0; i < count; ++i) {
        sys.x[i] = (float)(frand(&seed) * bx);
        sys.y[i] = (float)(frand(&seed) * by);
        sys.z[i] = (float)(z0 + frand(&seed) * (z1 - z0));
        sys.w[i] = 1.0f + (float)(i % 3);
        sys.s[i] = (two_classes && (i % 2)) ? 5.0f : 3.0f;
    }
    return sys;
}

static void free_system(test_system_t* sys) {
    free(sys->x); free(sys->y); free(sys->z); free(sys->w); free(sys->s);
}

UTEST(gisaxs, born_matches_direct_sum) {
    const double bx = 120.0, by = 100.0;
    test_system_t sys = make_system(60, bx, by, 10.0, 70.0, true);

    md_gisaxs_input_t input = {
        .count = sys.count,
        .x = sys.x, .y = sys.y, .z = sys.z,
        .weight = sys.w,
        .sigma = sys.s,
        .box_x = bx, .box_y = by,
    };
    md_gisaxs_params_t params = {
        .q_par_max = 0.4,
        .q_z_max = 0.4,
        .oversampling = 2.0,
    };

    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);

    md_gisaxs_info_t info;
    md_gisaxs_get_info(ctx, &info);
    EXPECT_EQ(info.num_classes, (size_t)2);

    ASSERT_TRUE(md_gisaxs_compute(ctx));

    const double qz[] = {0.0, 0.05, -0.12, 0.2, 0.33, 0.4};
    const size_t num_qz = ARRAY_SIZE(qz);
    const size_t R = md_gisaxs_num_rings(ctx);
    float*  I   = (float*)malloc(sizeof(float) * num_qz * R);
    double* ref = (double*)malloc(sizeof(double) * num_qz * R);

    md_gisaxs_model_t model = { .dwba = false };
    md_gisaxs_evaluate(ctx, &model, qz, num_qz, I);
    md_gisaxs_reference_born(ctx, &input, qz, num_qz, ref);

    double max_ref = 0.0;
    for (size_t i = 0; i < num_qz * R; ++i) max_ref = MAX(max_ref, ref[i]);

    double max_rel = 0.0;
    for (size_t i = 0; i < num_qz * R; ++i) {
        const double err = fabs(I[i] - ref[i]) / (ref[i] + 1.0e-4 * max_ref);
        max_rel = MAX(max_rel, err);
    }
    printf("born vs direct sum: max relative error %.3e\n", max_rel);
    EXPECT_LT(max_rel, 2.0e-2);

    free(I);
    free(ref);
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}

UTEST(gisaxs, dwba_without_contrast_equals_born) {
    const double bx = 150.0, by = 150.0;
    test_system_t sys = make_system(80, bx, by, 5.0, 60.0, false);

    md_gisaxs_input_t input = {
        .count = sys.count,
        .x = sys.x, .y = sys.y, .z = sys.z,
        .weight = sys.w,
        .sigma_uniform = 4.0f,
        .box_x = bx, .box_y = by,
    };
    md_gisaxs_params_t params = { .q_par_max = 0.3, .q_z_max = 0.3 };
    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(md_gisaxs_compute(ctx));

    const double lambda = 1.0;
    const double alpha_i = 0.3 * PI / 180.0;
    const double k0 = 2.0 * PI / lambda;
    const double p = k0 * sin(alpha_i);
    double qz[16];
    for (int i = 0; i < 16; ++i) qz[i] = p + 0.28 * i / 15.0;

    const size_t R = md_gisaxs_num_rings(ctx);
    float* I_dwba = (float*)malloc(sizeof(float) * 16 * R);
    float* I_born = (float*)malloc(sizeof(float) * 16 * R);

    md_gisaxs_model_t model = {
        .wavelength = lambda,
        .alpha_i = alpha_i,
        .dwba = true,
        .z_substrate = 0.0,
        .sld_ambient = 1.0e-5,
        .sld_substrate = 1.0e-5,
    };
    md_gisaxs_evaluate(ctx, &model, qz, 16, I_dwba);
    model.dwba = false;
    md_gisaxs_evaluate(ctx, &model, qz, 16, I_born);

    double max_rel = 0.0;
    for (size_t i = 0; i < 16 * R; ++i) {
        max_rel = MAX(max_rel, fabs(I_dwba[i] - I_born[i]) / (I_born[i] + 1e-30));
    }
    printf("dwba (no contrast) vs born: max relative error %.3e\n", max_rel);
    EXPECT_LT(max_rel, 1.0e-3);

    free(I_dwba);
    free(I_born);
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}

UTEST(gisaxs, dwba_yoneda_peak) {
    // Particles on a silicon substrate. The exit angle dependence of the diffuse scattering at fixed q_par should
    // be dominated by the transmission function |1 + r_f|^2 which peaks at the critical angle (Yoneda peak).
    const double bx = 400.0, by = 400.0;
    test_system_t sys = make_system(200, bx, by, 10.0, 30.0, false);
    md_gisaxs_input_t input = {
        .count = sys.count,
        .x = sys.x, .y = sys.y, .z = sys.z,
        .weight = sys.w,
        .sigma_uniform = 5.0f,
        .box_x = bx, .box_y = by,
    };
    md_gisaxs_params_t params = { .q_par_max = 0.05, .q_z_max = 0.12 };
    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(md_gisaxs_compute(ctx));

    const double lambda = 1.0332;   // 12 keV
    const double k0 = 2.0 * PI / lambda;
    const double sld_si = MD_GISAXS_R_E * 0.6991;
    const double alpha_c = sqrt(4.0 * PI * sld_si) / k0;   // sin(alpha_c) ~ alpha_c
    const double alpha_i = 2.0 * alpha_c;

    md_gisaxs_model_t model = {
        .wavelength = lambda,
        .alpha_i = alpha_i,
        .dwba = true,
        .z_substrate = 0.0,
        .sld_ambient = 0.0,
        .sld_substrate = sld_si,
        .sld_substrate_abs = 2.0 * PI * 1.0e-7 / (lambda * lambda),
    };

    enum { N = 200 };
    double qz[N];
    double af[N];
    for (int i = 0; i < N; ++i) {
        af[i] = 3.0 * alpha_c * (i + 1) / N;
        qz[i] = k0 * (sin(alpha_i) + sin(af[i]));
    }
    const size_t R = md_gisaxs_num_rings(ctx);
    float* I = (float*)malloc(sizeof(float) * N * R);
    md_gisaxs_evaluate(ctx, &model, qz, N, I);

    // Sum over rings for robustness
    int best = 0;
    double best_val = -1.0;
    for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (size_t r = 0; r < R; ++r) s += I[i * R + r];
        if (s > best_val) { best_val = s; best = i; }
    }
    const double ratio = af[best] / alpha_c;
    printf("yoneda peak at alpha_f / alpha_c = %.3f\n", ratio);
    EXPECT_GT(ratio, 0.85);
    EXPECT_LT(ratio, 1.25);

    free(I);
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}

// Independent brute force DWBA for particles in the ambient above a single interface (full plane average over +-q_par)
typedef struct { double re, im; } tc_t;
static tc_t tc_mul(tc_t a, tc_t b) { tc_t r = {a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re}; return r; }
static tc_t tc_add(tc_t a, tc_t b) { tc_t r = {a.re + b.re, a.im + b.im}; return r; }
static tc_t tc_div(tc_t a, tc_t b) { double d = b.re*b.re + b.im*b.im; tc_t r = {(a.re*b.re + a.im*b.im)/d, (a.im*b.re - a.re*b.im)/d}; return r; }
static tc_t tc_sqrt(tc_t a) {
    double m = sqrt(a.re*a.re + a.im*a.im);
    tc_t r = {sqrt(0.5*(m + a.re)), sqrt(0.5*(m - a.re))};
    if (a.im < 0) r.im = -r.im;
    return r;
}
static tc_t fresnel(double kz, double dsld_re, double sld_abs) {
    tc_t kz2 = {kz*kz - 4.0*PI*dsld_re, 4.0*PI*sld_abs};
    tc_t ks = tc_sqrt(kz2);
    tc_t num = {kz - ks.re, -ks.im};
    tc_t den = {kz + ks.re, ks.im};
    return tc_div(num, den);
}

UTEST(gisaxs, dwba_matches_direct_sum) {
    const double bx = 110.0, by = 90.0;
    test_system_t sys = make_system(40, bx, by, 8.0, 50.0, true);
    md_gisaxs_input_t input = {
        .count = sys.count,
        .x = sys.x, .y = sys.y, .z = sys.z,
        .weight = sys.w,
        .sigma = sys.s,
        .box_x = bx, .box_y = by,
    };
    md_gisaxs_params_t params = { .q_par_max = 0.35, .q_z_max = 0.35 };
    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(md_gisaxs_compute(ctx));

    const double lambda = 1.0332;
    const double k0 = 2.0 * PI / lambda;
    const double alpha_i = 0.25 * PI / 180.0;
    const double p = k0 * sin(alpha_i);
    const double z_sub = 2.0;
    const double sld_amb = 0.0;
    const double sld_sub = MD_GISAXS_R_E * 0.6991;
    const double sld_abs = 2.0 * PI * 1.7e-7 / (lambda * lambda);

    md_gisaxs_model_t model = {
        .wavelength = lambda, .alpha_i = alpha_i, .dwba = true, .z_substrate = z_sub,
        .sld_ambient = sld_amb, .sld_substrate = sld_sub, .sld_substrate_abs = sld_abs,
    };

    enum { N = 12 };
    double qz[N];
    for (int i = 0; i < N; ++i) qz[i] = p + 0.003 + 0.3 * i / (N - 1);

    const size_t R = md_gisaxs_num_rings(ctx);
    float* I = (float*)malloc(sizeof(float) * N * R);
    md_gisaxs_evaluate(ctx, &model, qz, N, I);

    // Reference, reusing the ring point definition from the Born reference by evaluating each term separately is not
    // possible (cross terms), so the ring points are reconstructed here from the same rules.
    md_gisaxs_info_t info;
    md_gisaxs_get_info(ctx, &info);
    const double dqx = 2.0 * PI / bx, dqy = 2.0 * PI / by;
    const double dq = MAX(dqx, dqy);
    double* ref = (double*)calloc(N * R, sizeof(double));
    double* cnt = (double*)calloc(R, sizeof(double));
    const int kxm = (int)floor(params.q_par_max / dqx), kym = (int)floor(params.q_par_max / dqy);
    for (int kx = -kxm; kx <= kxm; ++kx) {
        for (int ky = -kym; ky <= kym; ++ky) {
            if (kx == 0 && ky == 0) continue;
            const double qx = kx * dqx, qy = ky * dqy;
            const double qp = sqrt(qx*qx + qy*qy);
            if (qp > params.q_par_max) continue;
            long r = (long)floor(qp / dq + 0.5);
            if (r < 1) r = 1;
            r -= 1;
            if (r >= (long)R) continue;
            cnt[r] += 1.0;
            for (int iq = 0; iq < N; ++iq) {
                const double q = qz[iq] - p;
                tc_t ri = fresnel(p, sld_sub - sld_amb, sld_abs);
                tc_t rf = fresnel(q, sld_sub - sld_amb, sld_abs);
                const double Qs[4] = {p + q, q - p, p - q, -(p + q)};
                tc_t coef[4] = {{1,0}, ri, rf, tc_mul(ri, rf)};
                tc_t F = {0,0};
                for (int t = 0; t < 4; ++t) {
                    tc_t Ft = {0,0};
                    for (size_t j = 0; j < sys.count; ++j) {
                        const double s = sys.s[j];
                        const double amp = sys.w[j] * exp(-0.5 * (qp*qp + Qs[t]*Qs[t]) * s * s);
                        const double ph = -(qx * sys.x[j] + qy * sys.y[j] + Qs[t] * (sys.z[j] - z_sub));
                        tc_t e = {amp * cos(ph), amp * sin(ph)};
                        Ft = tc_add(Ft, e);
                    }
                    F = tc_add(F, tc_mul(coef[t], Ft));
                }
                ref[iq * R + r] += F.re*F.re + F.im*F.im;
            }
        }
    }
    double max_ref = 0.0;
    for (size_t i = 0; i < N * R; ++i) {
        const size_t r = i % R;
        ref[i] = cnt[r] > 0 ? ref[i] / cnt[r] / (bx * by) : 0.0;
        max_ref = MAX(max_ref, ref[i]);
    }
    double max_rel = 0.0;
    for (size_t i = 0; i < N * R; ++i) {
        max_rel = MAX(max_rel, fabs(I[i] - ref[i]) / (ref[i] + 1.0e-4 * max_ref));
    }
    printf("dwba vs direct sum: max relative error %.3e\n", max_rel);
    EXPECT_LT(max_rel, 2.0e-2);

    free(ref);
    free(cnt);
    free(I);
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}

UTEST(gisaxs, graded_without_film_contrast_equals_simple) {
    const double bx = 200.0, by = 200.0;
    test_system_t sys = make_system(150, bx, by, 10.0, 120.0, false);
    md_gisaxs_input_t input = {
        .count = sys.count, .x = sys.x, .y = sys.y, .z = sys.z, .weight = sys.w,
        .sigma_uniform = 4.0f, .box_x = bx, .box_y = by,
    };
    md_gisaxs_params_t params = { .q_par_max = 0.2, .q_z_max = 0.2 };
    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(md_gisaxs_compute(ctx));

    const double lambda = 1.0;
    const double k0 = 2.0 * PI / lambda;
    const double alpha_i = 0.3 * PI / 180.0;
    enum { N = 24 };
    double qz[N];
    for (int i = 0; i < N; ++i) qz[i] = k0 * sin(alpha_i) + 0.19 * i / (N - 1);

    md_gisaxs_model_t model = {
        .wavelength = lambda, .alpha_i = alpha_i, .dwba = true, .z_substrate = 5.0,
        .sld_ambient = 0.0, .sld_substrate = 2.0e-5, .sld_substrate_abs = 3.0e-7,
        .substrate_roughness = 3.0,
    };
    const size_t R = md_gisaxs_num_rings(ctx);
    float* a = (float*)malloc(sizeof(float) * N * R);
    float* b = (float*)malloc(sizeof(float) * N * R);
    md_gisaxs_evaluate(ctx, &model, qz, N, a);
    model.graded = true;
    model.profile_sld_scale = 0.0;
    model.profile_abs_scale = 0.0;
    md_gisaxs_evaluate(ctx, &model, qz, N, b);

    double max_rel = 0.0;
    for (size_t i = 0; i < N * R; ++i) {
        max_rel = MAX(max_rel, fabs(a[i] - b[i]) / (fabs(a[i]) + 1e-30));
    }
    printf("graded (no film contrast) vs simple: max relative error %.3e\n", max_rel);
    // Same physics, different phase bookkeeping; differences are float rounding in c^H S c
    EXPECT_LT(max_rel, 1.0e-3);

    free(a); free(b);
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}

UTEST(gisaxs, reflectivity_single_interface) {
    // Fresnel reflectivity of a bare substrate: total reflection below the critical angle, analytic above
    const double bx = 100.0, by = 100.0;
    test_system_t sys = make_system(10, bx, by, 10.0, 20.0, false);
    md_gisaxs_input_t input = { .count = sys.count, .x = sys.x, .y = sys.y, .z = sys.z, .weight = sys.w,
        .sigma_uniform = 3.0f, .box_x = bx, .box_y = by };
    md_gisaxs_params_t params = { .q_par_max = 0.2, .q_z_max = 0.2 };
    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);

    const double sld = 2.0e-5;
    const double qc = 4.0 * sqrt(PI * sld);
    md_gisaxs_model_t model = { .wavelength = 1.0, .dwba = true, .z_substrate = 0.0, .sld_substrate = sld };
    double qz[3] = {0.5 * qc, 2.0 * qc, 5.0 * qc};
    double Rf[3];
    md_gisaxs_reflectivity(ctx, &model, qz, 3, Rf);
    EXPECT_NEAR(Rf[0], 1.0, 1e-6);
    for (int i = 1; i < 3; ++i) {
        const double a = qz[i], b = sqrt(qz[i]*qz[i] - qc*qc);
        const double ref = ((a - b) / (a + b)) * ((a - b) / (a + b));
        EXPECT_NEAR(Rf[i] / ref, 1.0, 1e-6);
    }
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}

UTEST(gisaxs, graded_film_yoneda_peak) {
    // A dense particle film without substrate contrast: the diffuse scattering should show a Yoneda peak at the
    // critical angle of the film itself, which only the graded reference medium can produce.
    const double bx = 300.0, by = 300.0;
    const size_t count = 20000;
    const double z0 = 10.0, z1 = 410.0;
    test_system_t sys = make_system(count, bx, by, z0, z1, false);
    const double rho_film = 0.5;    // e/Å^3
    const float w = (float)(rho_film * bx * by * (z1 - z0) / count);
    for (size_t i = 0; i < count; ++i) sys.w[i] = w;

    md_gisaxs_input_t input = { .count = sys.count, .x = sys.x, .y = sys.y, .z = sys.z, .weight = sys.w,
        .sigma_uniform = 3.0f, .box_x = bx, .box_y = by };
    md_gisaxs_params_t params = { .q_par_max = 0.1, .q_z_max = 0.1 };
    md_gisaxs_t* ctx = md_gisaxs_create(&input, &params, md_get_heap_allocator());
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(md_gisaxs_compute(ctx));

    const double lambda = 1.0;
    const double k0 = 2.0 * PI / lambda;
    const double sld_film = MD_GISAXS_R_E * rho_film;
    const double alpha_c = asin(sqrt(4.0 * PI * sld_film) / k0);
    const double alpha_i = 3.0 * alpha_c;

    md_gisaxs_model_t model = {
        .wavelength = lambda, .alpha_i = alpha_i, .dwba = true, .graded = true, .z_substrate = z0 - 5.0,
        .sld_ambient = 0.0, .sld_substrate = 0.0,
        .profile_sld_scale = MD_GISAXS_R_E, .profile_abs_scale = 0.0,
    };
    enum { N = 200 };
    double qz[N], af[N];
    for (int i = 0; i < N; ++i) {
        af[i] = 3.0 * alpha_c * (i + 1) / N;
        qz[i] = k0 * (sin(alpha_i) + sin(af[i]));
    }
    const size_t R = md_gisaxs_num_rings(ctx);
    float* I = (float*)malloc(sizeof(float) * N * R);
    md_gisaxs_evaluate(ctx, &model, qz, N, I);
    int best = 0; double best_val = -1.0;
    for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (size_t r = 0; r < R; ++r) s += I[i * R + r];
        if (s > best_val) { best_val = s; best = i; }
    }
    const double ratio = af[best] / alpha_c;
    printf("graded film yoneda peak at alpha_f / alpha_c(film) = %.3f\n", ratio);
    EXPECT_GT(ratio, 0.85);
    EXPECT_LT(ratio, 1.25);

    // Kiessig fringes of the film: spacing of reflectivity minima ~ 2 pi / thickness
    {
        enum { NR = 4000 };
        double* qr = (double*)malloc(sizeof(double) * NR);
        double* Rr = (double*)malloc(sizeof(double) * NR);
        const double qc = 4.0 * sqrt(PI * sld_film);
        for (int i = 0; i < NR; ++i) qr[i] = 1.5 * qc + 0.1 * i / (NR - 1);
        md_gisaxs_reflectivity(ctx, &model, qr, NR, Rr);
        double first = -1, last = -1; int nmin = 0;
        for (int i = 1; i < NR - 1; ++i) {
            if (Rr[i] < Rr[i-1] && Rr[i] < Rr[i+1]) {
                if (first < 0) first = qr[i];
                last = qr[i];
                nmin++;
            }
        }
        ASSERT_GT(nmin, 3);
        const double spacing = (last - first) / (nmin - 1);
        const double expected = 2.0 * PI / (z1 - z0);
        printf("kiessig spacing %.5f, expected %.5f\n", spacing, expected);
        EXPECT_NEAR(spacing / expected, 1.0, 0.1);
        free(qr); free(Rr);
    }

    // Profile sanity: film density inside the slab
    const double* prof = md_gisaxs_slice_profile(ctx);
    const double* sz = md_gisaxs_slice_z(ctx);
    double sum = 0.0; int n = 0;
    for (size_t k = 0; k < md_gisaxs_num_slices(ctx); ++k) {
        if (sz[k] > z0 + 20 && sz[k] < z1 - 20) { sum += prof[k]; ++n; }
    }
    EXPECT_NEAR(sum / n, rho_film, 0.05 * rho_film);

    free(I);
    md_gisaxs_destroy(ctx);
    free_system(&sys);
}
