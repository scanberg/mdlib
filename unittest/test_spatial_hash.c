#include "utest.h"

#include <core/md_coord_stream.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_spatial_acc.h>
#include <core/md_str.h>
#include <core/md_intrinsics.h>
#include <core/md_os.h>
#include <core/md_hash.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_lammps.h>
#include <md_system.h>
#include <md_util.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <float.h>

typedef struct {
    uint32_t i, j;
    double d2;
} dist_pair_t;

static inline double rnd_rng(double min, double max);

typedef struct spatial_acc_data_t {
    md_array(dist_pair_t)* pairs;
    md_allocator_i* alloc;
} spatial_acc_data_t;

static void spatial_acc_neighbor_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx;
    (void)j_idx;
    uint32_t* count = (uint32_t*)user_param;

    const md_256 v_r2 = md_mm256_set1_ps(25.0f);  // 5.0^2

    const size_t vec_count = num_pairs & ~(size_t)7;
    for (size_t k = 0; k < vec_count; k += 8) {
        md_256 v_d2 = md_mm256_loadu_ps(ij_dist2 + k);
        md_256 v_mask = md_mm256_cmplt_ps(v_d2, v_r2);
        *count += popcnt32(md_mm256_movemask_ps(v_mask));
    }

    for (size_t k = vec_count; k < num_pairs; ++k) {
        *count += (ij_dist2[k] < 25.0f);
    }
}

static void spatial_acc_cutoff_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    spatial_acc_data_t* data = (spatial_acc_data_t*)user_param;
    for (size_t i = 0; i < num_pairs; i++) {
        dist_pair_t pair = {
            .i = i_idx[i],
            .j = j_idx[i],
            .d2 = ij_dist2[i]
        };
        md_array_push(*data->pairs, pair, data->alloc);
    }
}

static void spatial_acc_pair_count_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {  
    (void)i_idx;
    (void)j_idx;
    (void)ij_dist2;
    uint32_t* count = (uint32_t*)user_param;
    *count += (uint32_t)num_pairs;
}

static void spatial_acc_point_count_callback(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user_param) {
    (void)idx;
    (void)x;
    (void)y;
    (void)z;
    uint32_t* count = (uint32_t*)user_param;
    *count += (uint32_t)num_points;
}

// Points within a sphere, from the AABB query: the box around the sphere, then the distance to the image of the centre
// which the query works in. The library has no sphere query of its own.
typedef struct sphere_filter_t {
    md_spatial_acc_point_callback_t callback;
    void* user_param;
    float c[3];
    float r2;
} sphere_filter_t;

static void sphere_filter_callback(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user_param) {
    const sphere_filter_t* f = (const sphere_filter_t*)user_param;
    uint32_t ki[64];
    float kx[64], ky[64], kz[64];
    size_t n = 0;
    for (size_t k = 0; k < num_points; ++k) {
        const float dx = x[k] - f->c[0];
        const float dy = y[k] - f->c[1];
        const float dz = z[k] - f->c[2];
        if (dx * dx + dy * dy + dz * dz <= f->r2) {
            ki[n] = idx[k]; kx[n] = x[k]; ky[n] = y[k]; kz[n] = z[k];
            if (++n == 64) {
                f->callback(ki, kx, ky, kz, n, f->user_param);
                n = 0;
            }
        }
    }
    if (n) f->callback(ki, kx, ky, kz, n, f->user_param);
}

static void query_points_in_sphere(const md_spatial_acc_t* acc, const double center[3], double radius, md_spatial_acc_point_callback_t callback, void* user_param) {
    double qc[3];
    md_spatial_acc_aabb_query_center(qc, acc, center);
    const float r2 = (float)(radius * radius);
    sphere_filter_t f = { callback, user_param, { (float)qc[0], (float)qc[1], (float)qc[2] }, nextafterf(r2, r2 + 1.0f) };
    const double r3[3] = { radius, radius, radius };
    md_spatial_acc_for_each_point_in_aabb(acc, center, r3, sphere_filter_callback, &f);
}

typedef struct spatial_acc_point_collect_t {
    md_array(uint32_t) idx;
    md_allocator_i* alloc;
} spatial_acc_point_collect_t;

static void spatial_acc_point_collect_callback(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user_param) {
    (void)x;
    (void)y;
    (void)z;
    spatial_acc_point_collect_t* data = (spatial_acc_point_collect_t*)user_param;
    for (size_t i = 0; i < num_points; ++i) {
        md_array_push(data->idx, idx[i], data->alloc);
    }
}

static int cmp_u32_asc(const void* a, const void* b) {
    const uint32_t va = *(const uint32_t*)a;
    const uint32_t vb = *(const uint32_t*)b;
    return (va > vb) - (va < vb);
}

#define EXPECT_U32_SET_EQ(got_arr, exp_ptr, exp_count) do { \
    const size_t got_count__ = md_array_size(got_arr); \
    EXPECT_EQ((size_t)(exp_count), got_count__); \
    if (got_count__ == (size_t)(exp_count)) { \
        qsort((got_arr), got_count__, sizeof(uint32_t), cmp_u32_asc); \
        for (size_t i__ = 1; i__ < got_count__; ++i__) { \
            EXPECT_NE((got_arr)[i__ - 1], (got_arr)[i__]); \
        } \
        for (size_t i__ = 0; i__ < got_count__; ++i__) { \
            EXPECT_EQ((exp_ptr)[i__], (got_arr)[i__]); \
        } \
    } \
} while(0)

#define EXPECT_U32_MDARRAY_SET_EQ(got_arr, exp_arr) do { \
    const size_t got_count__ = md_array_size(got_arr); \
    const size_t exp_count__ = md_array_size(exp_arr); \
    EXPECT_EQ(exp_count__, got_count__); \
    if (got_count__ == exp_count__) { \
        qsort((got_arr), got_count__, sizeof(uint32_t), cmp_u32_asc); \
        qsort((exp_arr), exp_count__, sizeof(uint32_t), cmp_u32_asc); \
        for (size_t i__ = 1; i__ < got_count__; ++i__) { \
            EXPECT_NE((got_arr)[i__ - 1], (got_arr)[i__]); \
        } \
        for (size_t i__ = 0; i__ < got_count__; ++i__) { \
            EXPECT_EQ((exp_arr)[i__], (got_arr)[i__]); \
        } \
    } \
} while(0)

static inline bool point_in_aabb_cart(const double p[3], const double c[3], const double r[3]) {
    const double eps = 1.0e-6;
    return (fabs(p[0] - c[0]) <= r[0] + eps) && (fabs(p[1] - c[1]) <= r[1] + eps) && (fabs(p[2] - c[2]) <= r[2] + eps);
}

static inline double wrap_mic_ortho(double d, double L) {
    // Wrap into [-L/2, L/2] using nearest-integer convention
    return d - round(d / L) * L;
}

static inline void dmat3_mul(double out[3][3], const double a[3][3], const double b[3][3]) {
    // Matrices are indexed as [col][row] (column vectors are stored in the first index)
    // out = a * b
    // out[col][row] = sum_k a[k][row] * b[col][k]
    for (int col = 0; col < 3; ++col) {
        for (int row = 0; row < 3; ++row) {
            out[col][row] =
                a[0][row] * b[col][0] +
                a[1][row] * b[col][1] +
                a[2][row] * b[col][2];
        }
    }
}

static inline void dmat3_mul_vec3(double out[3], const double m[3][3], const double v[3]) {
    // Matrices are indexed as [col][row]
    // out = m * v
    for (int row = 0; row < 3; ++row) {
        out[row] =
            m[0][row] * v[0] +
            m[1][row] * v[1] +
            m[2][row] * v[2];
    }
}

// Convert cartesian coordinates to fractional coordinates using the inverse of the unit cell matrix
static inline void cart_to_fract(double out_s[3], const double in_x[3], const double I[3][3]) {
	dmat3_mul_vec3(out_s, I, in_x);
}

static inline void fract_to_cart(double out_x[3], const double in_s[3], const double A[3][3]) {
    dmat3_mul_vec3(out_x, A, in_s);
}

static inline double distance_ref_mic27(const double G[3][3], const double s0[3], const double s1[3]) {
    double ds0[3] = { s1[0] - s0[0], s1[1] - s0[1], s1[2] - s0[2] };

    // Start from the rounded guess
    ds0[0] -= round(ds0[0]);
    ds0[1] -= round(ds0[1]);
    ds0[2] -= round(ds0[2]);

    double best = DBL_MAX;
    for (int ix = -1; ix <= 1; ++ix)
    for (int iy = -1; iy <= 1; ++iy)
    for (int iz = -1; iz <= 1; ++iz) {
        double d[3] = { ds0[0] + ix, ds0[1] + iy, ds0[2] + iz };
        double d2 = 0.0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                d2 += G[a][b] * d[a] * d[b];
        if (d2 < best) best = d2;
    }
    return best;
}

UTEST(spatial_hash, small_periodic) {
    float x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    float y[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    md_unitcell_t cell = md_unitcell_from_extent(10, 0, 0);

    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, 10);
    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 10.0, &cell, 0);
    
    uint32_t count = 0;
    double p0[3] = {5, 0, 0};
    query_points_in_sphere(&acc, p0, 1.5f, spatial_acc_point_count_callback, &count);
    EXPECT_EQ(3, count);

    count = 0;
    double p1[3] = {8.5f, 0, 0};
    query_points_in_sphere(&acc, p1, 3, spatial_acc_point_count_callback, &count);
    EXPECT_EQ(6, count);

    md_spatial_acc_free(&acc);
}

UTEST(spatial_hash, aabb_periodic_ortho) {
    // Ortho periodic unit cell. Query an AABB that crosses x=0 seam.
    // Expected: points close to x=0 and x=L are both reported.
    float x[] = {0.10f, 9.90f, 9.80f, 5.00f};
    float y[] = {5.00f, 5.00f, 5.00f, 5.00f};
    float z[] = {5.00f, 5.00f, 5.00f, 5.00f};

    md_unitcell_t cell = md_unitcell_from_extent(10.0, 10.0, 10.0);
    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, ARRAY_SIZE(x));

    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);

    spatial_acc_point_collect_t data = {
        .idx = NULL,
        .alloc = md_get_heap_allocator(),
    };

    const double aabb_cen[3] = {0.20, 5.00, 5.00};
    const double aabb_rad[3] = {0.35, 0.20, 0.20};
    md_spatial_acc_for_each_point_in_aabb(&acc, aabb_cen, aabb_rad, spatial_acc_point_collect_callback, &data);

    const uint32_t exp[] = {0, 1};
    EXPECT_U32_SET_EQ(data.idx, exp, ARRAY_SIZE(exp));

    md_array_free(data.idx, data.alloc);
    md_spatial_acc_free(&acc);
}

UTEST(spatial_hash, aabb_periodic_triclinic) {
    // Triclinic periodic unit cell. Query an AABB that crosses the periodic seam
    // in fractional Y (which corresponds to a slanted shift in cartesian space).
    const double A[3][3] = {
        {10.0, 0.0, 0.0},
        {3.0,  9.0, 0.0},
        {2.0,  1.0, 8.0},
    };

    md_unitcell_t cell = md_unitcell_from_matrix_double(A);

    double s0[3] = {0.50, 0.02, 0.50};
    double s1[3] = {0.50, 0.98, 0.50};
    double s2[3] = {0.50, 0.50, 0.50};

    double x0[3];
    double x1[3];
    double x2[3];
    fract_to_cart(x0, s0, A);
    fract_to_cart(x1, s1, A);
    fract_to_cart(x2, s2, A);

    float x[3] = {(float)x0[0], (float)x1[0], (float)x2[0]};
    float y[3] = {(float)x0[1], (float)x1[1], (float)x2[1]};
    float z[3] = {(float)x0[2], (float)x1[2], (float)x2[2]};

    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, 3);
    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);

    spatial_acc_point_collect_t data = {
        .idx = NULL,
        .alloc = md_get_heap_allocator(),
    };

    // AABB around x0. s1 is far in cartesian, but its periodic image (shifted by -b)
    // is close to x0 and must be reported.
    const double aabb_cen[3] = {x0[0], x0[1], x0[2]};
    const double aabb_rad[3] = {0.25, 0.50, 0.25};
    md_spatial_acc_for_each_point_in_aabb(&acc, aabb_cen, aabb_rad, spatial_acc_point_collect_callback, &data);

    const uint32_t exp[] = {0, 1};
    EXPECT_U32_SET_EQ(data.idx, exp, ARRAY_SIZE(exp));

    md_array_free(data.idx, data.alloc);
    md_spatial_acc_free(&acc);
}

UTEST(spatial_hash, aabb_periodic_ortho_randomized_reference) {
    // Robust randomized test for periodic ortho AABB queries.
    // Compares `md_spatial_acc_for_each_point_in_aabb` against a brute-force MIC reference.
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    const double Lx = 50.0;
    const double Ly = 60.0;
    const double Lz = 70.0;
    md_unitcell_t cell = md_unitcell_from_extent(Lx, Ly, Lz);

    const size_t N = 4096;
    float* x = (float*)md_temp_alloc(temp, N * sizeof(float));
    float* y = (float*)md_temp_alloc(temp, N * sizeof(float));
    float* z = (float*)md_temp_alloc(temp, N * sizeof(float));

    srand(1337);
    for (size_t i = 0; i < N; ++i) {
        x[i] = (float)rnd_rng(0.0, Lx);
        y[i] = (float)rnd_rng(0.0, Ly);
        z[i] = (float)rnd_rng(0.0, Lz);
    }

    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, N);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);

    spatial_acc_point_collect_t got = { .idx = NULL, .alloc = alloc };
    uint8_t* seen = (uint8_t*)md_temp_alloc(temp, N);
    ASSERT_TRUE(seen);

    const int iters = 250;
    for (int iter = 0; iter < iters; ++iter) {
        md_array_shrink(got.idx, 0);
        memset(seen, 0, N);

        const double cen[3] = { rnd_rng(0.0, Lx), rnd_rng(0.0, Ly), rnd_rng(0.0, Lz) };
        // Keep radii < 0.5 box length to ensure MIC reference is sufficient.
        const double rad[3] = { rnd_rng(0.0, 0.45 * Lx), rnd_rng(0.0, 0.45 * Ly), rnd_rng(0.0, 0.45 * Lz) };

        const double eps = MAX(1.0e-6, 128.0 * (double)FLT_EPSILON * (Lx + Ly + Lz + rad[0] + rad[1] + rad[2] + 1.0));

        md_spatial_acc_for_each_point_in_aabb(&acc, cen, rad, spatial_acc_point_collect_callback, &got);

        for (size_t k = 0; k < md_array_size(got.idx); ++k) {
            const uint32_t idx = got.idx[k];
            EXPECT_LT(idx, (uint32_t)N);
            EXPECT_EQ(0, seen[idx]);
            seen[idx] = 1;
        }

        // Reference classification with slack:
        //   margin < -eps => definitely inside => must be reported
        //   margin > +eps => definitely outside => must NOT be reported
        //   otherwise ambiguous near boundary => accept either
        for (uint32_t i = 0; i < (uint32_t)N; ++i) {
            const double dx = wrap_mic_ortho((double)x[i] - cen[0], Lx);
            const double dy = wrap_mic_ortho((double)y[i] - cen[1], Ly);
            const double dz = wrap_mic_ortho((double)z[i] - cen[2], Lz);
            const double mx = fabs(dx) - rad[0];
            const double my = fabs(dy) - rad[1];
            const double mz = fabs(dz) - rad[2];
            const double margin = MAX(mx, MAX(my, mz));

            if (margin < -eps) {
                EXPECT_EQ(1, seen[i]);
            } else if (margin > eps) {
                EXPECT_EQ(0, seen[i]);
            }
        }
    }

    md_array_free(got.idx, got.alloc);
    md_spatial_acc_free(&acc);
    md_temp_end(temp);
}

UTEST(spatial_hash, aabb_periodic_triclinic_randomized_reference) {
    // Robust randomized test for periodic triclinic AABB queries.
    // Reference uses brute-force over 27 periodic images (sufficient for chosen small radii).
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    const double A[3][3] = {
        {30.0,  0.0,  0.0},
        {10.0, 28.0,  0.0},
        {-5.0,  7.0, 22.0},
    };
    md_unitcell_t cell = md_unitcell_from_matrix_double(A);

    const double a_vec[3] = {A[0][0], A[0][1], A[0][2]};
    const double b_vec[3] = {A[1][0], A[1][1], A[1][2]};
    const double c_vec[3] = {A[2][0], A[2][1], A[2][2]};

    const size_t N = 2048;
    float* x = (float*)md_temp_alloc(temp, N * sizeof(float));
    float* y = (float*)md_temp_alloc(temp, N * sizeof(float));
    float* z = (float*)md_temp_alloc(temp, N * sizeof(float));

    srand(7331);
    for (size_t i = 0; i < N; ++i) {
        const double s[3] = { rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0) };
        double p[3];
        fract_to_cart(p, s, A);
        x[i] = (float)p[0];
        y[i] = (float)p[1];
        z[i] = (float)p[2];
    }

    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, N);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);

    spatial_acc_point_collect_t got = { .idx = NULL, .alloc = alloc };
    uint8_t* seen = (uint8_t*)md_temp_alloc(temp, N);
    ASSERT_TRUE(seen);

    const int iters = 200;
    for (int iter = 0; iter < iters; ++iter) {
        md_array_shrink(got.idx, 0);
        memset(seen, 0, N);

        const double sc[3] = { rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0) };
        double cen[3];
        fract_to_cart(cen, sc, A);

        // Keep radii small enough that checking 27 images is sufficient.
        const double rad[3] = { rnd_rng(0.0, 6.0), rnd_rng(0.0, 6.0), rnd_rng(0.0, 6.0) };

        const double eps = MAX(1.0e-6, 256.0 * (double)FLT_EPSILON * (fabs(cen[0]) + fabs(cen[1]) + fabs(cen[2]) + rad[0] + rad[1] + rad[2] + 1.0));

        md_spatial_acc_for_each_point_in_aabb(&acc, cen, rad, spatial_acc_point_collect_callback, &got);

        for (size_t k = 0; k < md_array_size(got.idx); ++k) {
            const uint32_t idx = got.idx[k];
            EXPECT_LT(idx, (uint32_t)N);
            EXPECT_EQ(0, seen[idx]);
            seen[idx] = 1;
        }

        // Reference classification with slack using best (minimum) margin over 27 periodic images.
        // margin < -eps => definitely inside => must be reported
        // margin > +eps => definitely outside => must NOT be reported
        for (uint32_t i = 0; i < (uint32_t)N; ++i) {
            const double p0[3] = { (double)x[i], (double)y[i], (double)z[i] };
            double best_margin = DBL_MAX;
            for (int ia = -1; ia <= 1; ++ia) {
                for (int ib = -1; ib <= 1; ++ib) {
                    for (int ic = -1; ic <= 1; ++ic) {
                        const double shift[3] = {
                            ia * a_vec[0] + ib * b_vec[0] + ic * c_vec[0],
                            ia * a_vec[1] + ib * b_vec[1] + ic * c_vec[1],
                            ia * a_vec[2] + ib * b_vec[2] + ic * c_vec[2],
                        };
                        const double p[3] = { p0[0] + shift[0], p0[1] + shift[1], p0[2] + shift[2] };
                        const double mx = fabs(p[0] - cen[0]) - rad[0];
                        const double my = fabs(p[1] - cen[1]) - rad[1];
                        const double mz = fabs(p[2] - cen[2]) - rad[2];
                        const double margin = MAX(mx, MAX(my, mz));
                        best_margin = MIN(best_margin, margin);
                    }
                }
            }

            if (best_margin < -eps) {
                EXPECT_EQ(1, seen[i]);
            } else if (best_margin > eps) {
                EXPECT_EQ(0, seen[i]);
            }
        }
    }

    md_array_free(got.idx, got.alloc);
    md_spatial_acc_free(&acc);
    md_temp_end(temp);
}

static size_t do_brute_force_double(const float* in_x, const float* in_y, const float* in_z, size_t num_points, double cutoff, const double G[3][3], const double I[3][3], md_array(dist_pair_t)* pairs, md_allocator_i* alloc) {
    size_t count = 0;
    const double r2 = cutoff * cutoff;

    for (size_t i = 0; i < num_points - 1; ++i) {
        // Fractional coords of i
        double xi[3] = { in_x[i], in_y[i], in_z[i] };
        double si[3];
        cart_to_fract(si, xi, I);
        for (size_t j = i + 1; j < num_points; ++j) {
            double xj[3] = { in_x[j], in_y[j], in_z[j] };
            double sj[3];
            cart_to_fract(sj, xj, I);
            double d2 = distance_ref_mic27(G, si, sj);
            if (d2 < r2) {
                if (pairs) {
                    dist_pair_t pair = { .i = (uint32_t)i, .j = (uint32_t)j, .d2 = d2 };
                    md_array_push(*pairs, pair, alloc);
                }
                count += 1;
            }
        }
    }
    return count;
}

static inline double rnd_rng(double min, double max) {
    double r = ((double)rand() / (double)RAND_MAX);
    return r * (max - min) + min;
}

// The point queries hand their coordinates to the callback, and those have to be usable: cartesian, an
// image of the point that was put in, and inside the queried region around the image of the centre the
// query works in (md_spatial_acc_aabb_query_center). Tests which only look at the indices cannot tell
// fractional coordinates from cartesian ones, which the triclinic AABB query used to hand out, and the
// triclinic sphere query in its last batch.
typedef struct spatial_acc_point_coord_collect_t {
    md_array(uint32_t) idx;
    md_array(vec3_t)   xyz;
    md_allocator_i*    alloc;
} spatial_acc_point_coord_collect_t;

static void spatial_acc_point_coord_collect_callback(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user_param) {
    spatial_acc_point_coord_collect_t* data = (spatial_acc_point_coord_collect_t*)user_param;
    for (size_t i = 0; i < num_points; ++i) {
        md_array_push(data->idx, idx[i], data->alloc);
        md_array_push(data->xyz, vec3_set(x[i], y[i], z[i]), data->alloc);
    }
}

// Largest deviation of 'd' from the nearest lattice vector, as a length
static double lattice_residual(const double A[3][3], const double I[3][3], const double d[3]) {
    double s[3];
    dmat3_mul_vec3(s, I, d);
    for (int k = 0; k < 3; ++k) s[k] -= round(s[k]);
    double r[3];
    dmat3_mul_vec3(r, A, s);
    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

static void point_query_coordinates(int* utest_result, const double A[3][3]) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_unitcell_t cell = md_unitcell_from_matrix_double(A);
    double I[3][3];
    md_unitcell_I_extract_double(I, &cell);

    // Enough points that a query overflows the staging buffer, so the in loop flushes run as well as
    // the tail flush
    const size_t N = 20000;
    float* x = (float*)md_temp_alloc(temp, N * sizeof(float));
    float* y = (float*)md_temp_alloc(temp, N * sizeof(float));
    float* z = (float*)md_temp_alloc(temp, N * sizeof(float));
    srand(4242);
    for (size_t i = 0; i < N; ++i) {
        const double s[3] = { rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0) };
        double p[3];
        fract_to_cart(p, s, A);
        x[i] = (float)p[0]; y[i] = (float)p[1]; z[i] = (float)p[2];
    }

    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, N);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 4.0, &cell, 0);

    spatial_acc_point_coord_collect_t got = { .alloc = alloc };
    size_t total = 0;

    for (int iter = 0; iter < 60; ++iter) {
        // Centres inside the cell and well outside it
        const double sc[3] = { rnd_rng(-2.0, 3.0), rnd_rng(-2.0, 3.0), rnd_rng(-2.0, 3.0) };
        double cen[3];
        fract_to_cart(cen, sc, A);
        const double rad = rnd_rng(2.0, 8.0);

        double qc[3];
        md_spatial_acc_aabb_query_center(qc, &acc, cen);

        // The query image is the centre moved by a lattice vector
        const double dc[3] = { qc[0] - cen[0], qc[1] - cen[1], qc[2] - cen[2] };
        EXPECT_LT(lattice_residual(A, I, dc), 1.0e-3);

        for (int kind = 0; kind < 2; ++kind) {
            md_array_shrink(got.idx, 0);
            md_array_shrink(got.xyz, 0);
            if (kind == 0) {
                const double r3[3] = { rad, rad, rad };
                md_spatial_acc_for_each_point_in_aabb(&acc, cen, r3, spatial_acc_point_coord_collect_callback, &got);
            } else {
                query_points_in_sphere(&acc, cen, rad, spatial_acc_point_coord_collect_callback, &got);
            }

            size_t bad_image = 0, bad_region = 0;
            for (size_t k = 0; k < md_array_size(got.idx); ++k) {
                const uint32_t i = got.idx[k];
                const vec3_t p = got.xyz[k];
                const double d_in[3] = { p.x - x[i], p.y - y[i], p.z - z[i] };
                if (lattice_residual(A, I, d_in) > 1.0e-3) bad_image += 1;

                const double d[3] = { p.x - qc[0], p.y - qc[1], p.z - qc[2] };
                const double eps = 1.0e-3;
                if (kind == 0) {
                    if (fabs(d[0]) > rad + eps || fabs(d[1]) > rad + eps || fabs(d[2]) > rad + eps) bad_region += 1;
                } else {
                    if (sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) > rad + eps) bad_region += 1;
                }
            }
            if (bad_image || bad_region) {
                printf("  %s query, centre (%.2f %.2f %.2f) r %.2f: %zu of %zu not an image of their input, %zu outside the region\n",
                    kind == 0 ? "aabb" : "sphere", cen[0], cen[1], cen[2], rad, bad_image, md_array_size(got.idx), bad_region);
            }
            EXPECT_EQ((size_t)0, bad_image);
            EXPECT_EQ((size_t)0, bad_region);
            total += md_array_size(got.idx);
        }
    }
    EXPECT_GT(total, (size_t)0);

    md_temp_end(temp);
}

UTEST(spatial_hash, point_query_coordinates_ortho) {
    const double A[3][3] = {
        {31.0,  0.0,  0.0},
        { 0.0, 27.0,  0.0},
        { 0.0,  0.0, 24.0},
    };
    point_query_coordinates(utest_result, A);
}

UTEST(spatial_hash, point_query_coordinates_triclinic) {
    const double A[3][3] = {
        {30.0,  0.0,  0.0},
        {10.0, 28.0,  0.0},
        {-5.0,  7.0, 22.0},
    };
    point_query_coordinates(utest_result, A);
}

UTEST(spatial_hash, n2) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));

    {

#define TEST_COUNT 2048
        float x[TEST_COUNT];
        float y[TEST_COUNT];
        float z[TEST_COUNT];

        const double A[3][3] = {
            {40.0,  0.0,   0.0},
            {10.0, 50.0,   0.0},
            {-20.0, -10.0,  60.0}
		};

        srand(0);

        md_unitcell_t test_cell = md_unitcell_from_matrix_double(A);
        double G[3][3];
        double I[3][3];
        md_unitcell_G_extract_double(G, &test_cell);
        md_unitcell_I_extract_double(I, &test_cell);

        double X[3][3];
		// Ensure that A*I = Identity
		dmat3_mul(X, A, I);

		EXPECT_NEAR(X[0][0], 1.0, 1.0e-16);
        EXPECT_NEAR(X[0][1], 0.0, 1.0e-16);
        EXPECT_NEAR(X[0][2], 0.0, 1.0e-16);

        EXPECT_NEAR(X[1][0], 0.0, 1.0e-16);
        EXPECT_NEAR(X[1][1], 1.0, 1.0e-16);
        EXPECT_NEAR(X[1][2], 0.0, 1.0e-16);

        EXPECT_NEAR(X[2][0], 0.0, 1.0e-16);
        EXPECT_NEAR(X[2][1], 0.0, 1.0e-16);
        EXPECT_NEAR(X[2][2], 1.0, 1.0e-16);

        // Generate random points
        for (size_t i = 0; i < TEST_COUNT; ++i) {
            const double cx = rnd_rng(0.0, 100.0);
            const double cy = rnd_rng(0.0, 100.0);
            const double cz = rnd_rng(0.0, 100.0);

            x[i] = cx;
            y[i] = cy;
            z[i] = cz;
        }

        md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, TEST_COUNT);
        md_spatial_acc_t sa = { .alloc = alloc };
		md_spatial_acc_init(&sa, &stream, 6.0, &test_cell, 0);

        EXPECT_EQ(sa.G00, (float)G[0][0]);
        EXPECT_EQ(sa.G11, (float)G[1][1]);
        EXPECT_EQ(sa.G22, (float)G[2][2]);

        EXPECT_EQ(sa.H01, (float)(2.0 * G[0][1]));
        EXPECT_EQ(sa.H02, (float)(2.0 * G[0][2]));
        EXPECT_EQ(sa.H12, (float)(2.0 * G[1][2]));

        EXPECT_EQ(sa.I[0][0], (float)I[0][0]);
        EXPECT_EQ(sa.I[0][1], (float)I[0][1]);
        EXPECT_EQ(sa.I[0][2], (float)I[0][2]);

        EXPECT_EQ(sa.I[1][0], (float)I[1][0]);
        EXPECT_EQ(sa.I[1][1], (float)I[1][1]);
        EXPECT_EQ(sa.I[1][2], (float)I[1][2]);

        EXPECT_EQ(sa.I[2][0], (float)I[2][0]);
        EXPECT_EQ(sa.I[2][1], (float)I[2][1]);
        EXPECT_EQ(sa.I[2][2], (float)I[2][2]);

        md_array(dist_pair_t) sa_pairs = NULL;
        md_array(dist_pair_t) bf_pairs = NULL;

        spatial_acc_data_t usr_data = {
            .pairs = &sa_pairs,
            .alloc = alloc
        };

        for (double rad = 3.0; rad <= 6.0; rad += 0.5) {
            md_array_shrink(bf_pairs, 0);
            md_array_shrink(sa_pairs, 0);

            size_t bf_count = do_brute_force_double(x, y, z, TEST_COUNT, rad, G, I, &bf_pairs, alloc);
            md_spatial_acc_for_each_internal_pair_within_cutoff(&sa, rad, spatial_acc_cutoff_callback, &usr_data);

            size_t sa_count = md_array_size(sa_pairs);
            //EXPECT_EQ(bf_count, sa_count);
            if (bf_count != sa_count) {
                printf("wierd expected: %zu, but got: %zu, cutoff: %f\n", bf_count, sa_count, rad);
                md_hashmap32_t map_sa = { .allocator = alloc };
                md_hashmap32_t map_bf = { .allocator = alloc };

                // populate reference hashmap with brute force pairs
                for (size_t i = 0; i < bf_count; ++i) {
                    uint64_t key = ((uint64_t)bf_pairs[i].i << 32) | bf_pairs[i].j;
                    md_hashmap_add(&map_bf, key, 1);
                }

                for (size_t idx = 0; idx < sa_count; ++idx) {
                    uint32_t i = MIN(sa_pairs[idx].i, sa_pairs[idx].j);
                    uint32_t j = MAX(sa_pairs[idx].i, sa_pairs[idx].j);
                    uint64_t key = ((uint64_t)i << 32) | j;
                    md_hashmap_add(&map_sa, key, 1);

                    if (!md_hashmap_get(&map_bf, key)) {
                        printf("SA only pair: %u %u %f\n", i, j, sqrt(sa_pairs[idx].d2));
                        double xi[3] = { x[i], y[i], z[i] };
                        double xj[3] = { x[j], y[j], z[j] };
                        double si[3], sj[3];
                        cart_to_fract(si, xi, I);
                        cart_to_fract(sj, xj, I);
                        double dist_ref = sqrt(distance_ref_mic27(G, si, sj));
                        printf("Reference distance: %f\n", dist_ref);
                    }
                }

                for (size_t idx = 0; idx < bf_count; ++idx) {
                    uint32_t i = bf_pairs[idx].i;
                    uint32_t j = bf_pairs[idx].j;
                    uint64_t key = ((uint64_t)i << 32) | j;
                    if (!md_hashmap_get(&map_sa, key)) {
                        printf("BF only pair: %u %u %f\n", i, j, sqrt(bf_pairs[idx].d2));
                    }
                }
            }
        }
#undef TEST_COUNT
    }

#if 1
    md_unitcell_t cell = sys_state.unitcell;
    // Pairs closer than 5 A. The coordinates sit on a 0.01 A grid, so a few pairs are exactly 5 A apart and fall on
    // either side of the strict test in the callback depending on rounding; hence the tolerance.
    const size_t expected_count = 3711875;

    double G[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &cell);
    md_unitcell_I_extract_double(I, &cell);

    md_tick_t start, end;
    uint32_t count = 0;

    // Spatial acc implementation
    start = md_tick_now();
    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)sys_state.xyz, sizeof(vec3_t), NULL, sys.atom.count);
    md_spatial_acc_t acc = {.alloc = alloc};
    md_spatial_acc_init(&acc, &stream, 5.0, &cell, 0);
    md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, 5.0, spatial_acc_neighbor_callback, &count);
	//md_spatial_acc_for_each_pair_within_cutoff(&acc, 5.0, spatial_acc_cutoff_callback, &count);
    end = md_tick_now();
    size_t sa_count = count;
    const size_t internal_count = sa_count;
    //end = md_tick_now();
    //printf("Spatial acc cell neighborhood: %f ms\n", md_tick_to_milliseconds(end - start));
    EXPECT_NEAR(expected_count, sa_count, 5);
    if (sa_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, sa_count);
    }

    start = md_tick_now();
    count = 0;
    md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &stream, 5.0, spatial_acc_neighbor_callback, &count, 0);
	end = md_tick_now();
    sa_count = count;
	//printf("Spatial acc external query: %f ms\n", md_tick_to_milliseconds(end - start));
    // Every point against the structure: each pair from both sides, and every point with itself. The two queries
    // compute the same distances, so this holds exactly whatever the rounding.
	size_t ext_expected_count = internal_count * 2 + sys.atom.count;
    EXPECT_EQ(ext_expected_count, sa_count);

#if 0
    // This is so slow that we don't want to run it by default, but it can be useful for validating the reference implementation
    // Brute force
    start = md_tick_now();
    size_t bf_count = do_brute_force_double(sys_state.xyz, sys.atom.count, 5.0, G, I, NULL, NULL);
    end = md_tick_now();
    printf("Brute force: %f ms\n", md_tick_to_milliseconds(end - start));
    EXPECT_EQ(expected_count, bf_count);
    if (bf_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, bf_count);
    }
#endif
#endif

    md_temp_end(temp);
}

struct spatial_hash {
    md_allocator_i* arena;
};

UTEST_F_SETUP(spatial_hash) {
    utest_fixture->arena = md_vm_arena_create(GIGABYTES(4));
}

UTEST_F_TEARDOWN(spatial_hash) {
    md_vm_arena_destroy(utest_fixture->arena);
}

static inline float rnd() {
    return rand() / (float)RAND_MAX;
}

UTEST_F(spatial_hash, test_correctness_centered) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));

    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)sys_state.xyz, sizeof(vec3_t), NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys_state.unitcell, 0);
    
    srand(31);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys_state.unitcell);
    md_unitcell_A_extract_double(A, &sys_state.unitcell);
    md_unitcell_I_extract_double(I, &sys_state.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        double radius = rnd() * 20;

        int ref_count = 0;
        // NOT 'rad2': something in the ARM64-only include chain (md_intrinsics.h pulls in
        // <intrin.h>, which on ARM64 drags in the NEON/intrinsic headers) defines rad2 as an
        // object-like macro expanding to an integer constant, so the declaration becomes
        // 'const double <literal> = ...'. x64 never sees that header, which is why this only
        // breaks on the arm64 build.
        const double radius_sq = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys_state.xyz[i].x, sys_state.xyz[i].y, sys_state.xyz[i].z };
            double si[3];
            cart_to_fract(si, xi, I);

            if (distance_ref_mic27(G, s0, si) < radius_sq) {
                ref_count += 1;
            }
        }

        uint32_t sa_count = 0;
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        md_coord_stream_t ext_stream = md_coord_stream_from_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }

        sa_count = 0;
		query_points_in_sphere(&acc, x0, radius, spatial_acc_point_count_callback, &sa_count);
		EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}

UTEST_F(spatial_hash, test_correctness_ala) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));

    srand(31);
    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)sys_state.xyz, sizeof(vec3_t), NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys_state.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys_state.unitcell);
    md_unitcell_A_extract_double(A, &sys_state.unitcell);
    md_unitcell_I_extract_double(I, &sys_state.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        double radius = rnd() * 20.0;

        int ref_count = 0;
        const double radius_sq = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys_state.xyz[i].x, sys_state.xyz[i].y, sys_state.xyz[i].z };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < radius_sq) {
                ref_count += 1;
            }
        }

        uint32_t sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_from_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        query_points_in_sphere(&acc, x0, (float)radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }
    }
}

UTEST_F(spatial_hash, test_correctness_water) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro")));

    srand(31);

    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)sys_state.xyz, sizeof(vec3_t), NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys_state.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys_state.unitcell);
    md_unitcell_A_extract_double(A, &sys_state.unitcell);
    md_unitcell_I_extract_double(I, &sys_state.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        double radius = rnd() * 20.0;

        int ref_count = 0;
        const double radius_sq = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys_state.xyz[i].x, sys_state.xyz[i].y, sys_state.xyz[i].z };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < radius_sq) {
                ref_count += 1;
            }
        }

        uint32_t sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_from_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        query_points_in_sphere(&acc, x0, radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}

UTEST_F(spatial_hash, test_correctness_water_ethane_triclinic) {
    md_allocator_i* alloc = utest_fixture->arena;

    const char** atom_formats = md_lammps_atom_format_strings();
    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Triclinic_Init.data"), atom_format));

    srand(31);

    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)sys_state.xyz, sizeof(vec3_t), NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys_state.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys_state.unitcell);
    md_unitcell_A_extract_double(A, &sys_state.unitcell);
    md_unitcell_I_extract_double(I, &sys_state.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        double radius = rnd() * 10.0;

        uint32_t ref_count = 0;
        uint32_t sa_count = 0;

#if 0
		// Do N^2 test as well to validate the reference implementation
        ref_count = do_brute_force_double(sys_state.xyz, sys.atom.count, radius, G, I, NULL, NULL);
		md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, (float)radius, spatial_acc_pair_count_callback, &sa_count);
        EXPECT_NEAR(ref_count, sa_count, 2);
#endif

        ref_count = 0;
        const double radius_sq = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys_state.xyz[i].x, sys_state.xyz[i].y, sys_state.xyz[i].z };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < radius_sq) {
                ref_count += 1;
            }
        }

        sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_from_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        query_points_in_sphere(&acc, x0, (float)radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}


UTEST_F(spatial_hash, npt_triclinic) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/npt.gro")));

    srand(31);

    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)sys_state.xyz, sizeof(vec3_t), NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys_state.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys_state.unitcell);
    md_unitcell_A_extract_double(A, &sys_state.unitcell);
    md_unitcell_I_extract_double(I, &sys_state.unitcell);

	EXPECT_EQ((float)G[0][0], acc.G00);
	EXPECT_EQ((float)G[1][1], acc.G11);
	EXPECT_EQ((float)G[2][2], acc.G22);

	EXPECT_EQ((float)A[0][0], acc.A[0][0]);
    EXPECT_EQ((float)A[0][1], acc.A[0][1]);
    EXPECT_EQ((float)A[0][2], acc.A[0][2]);

    EXPECT_EQ((float)A[1][0], acc.A[1][0]);
    EXPECT_EQ((float)A[1][1], acc.A[1][1]);
    EXPECT_EQ((float)A[1][2], acc.A[1][2]);

    EXPECT_EQ((float)A[2][0], acc.A[2][0]);
    EXPECT_EQ((float)A[2][1], acc.A[2][1]);
    EXPECT_EQ((float)A[2][2], acc.A[2][2]);

	EXPECT_EQ((float)I[0][0], acc.I[0][0]);
	EXPECT_EQ((float)I[0][1], acc.I[0][1]);
	EXPECT_EQ((float)I[0][2], acc.I[0][2]);

	EXPECT_EQ((float)I[1][0], acc.I[1][0]);
	EXPECT_EQ((float)I[1][1], acc.I[1][1]);
	EXPECT_EQ((float)I[1][2], acc.I[1][2]);

	EXPECT_EQ((float)I[2][0], acc.I[2][0]);
	EXPECT_EQ((float)I[2][1], acc.I[2][1]);
	EXPECT_EQ((float)I[2][2], acc.I[2][2]);

    double X[3][3];
    dmat3_mul(X, A, I);

	EXPECT_NEAR(X[0][0], 1.0, 1.0e-15);
	EXPECT_NEAR(X[0][1], 0.0, 1.0e-15);
	EXPECT_NEAR(X[0][2], 0.0, 1.0e-15);
	EXPECT_NEAR(X[1][0], 0.0, 1.0e-15);
	EXPECT_NEAR(X[1][1], 1.0, 1.0e-15);
	EXPECT_NEAR(X[1][2], 0.0, 1.0e-15);
	EXPECT_NEAR(X[2][0], 0.0, 1.0e-15);
	EXPECT_NEAR(X[2][1], 0.0, 1.0e-15);
	EXPECT_NEAR(X[2][2], 1.0, 1.0e-15);

    {
        int i0 = 835;
        int i1 = 6160;

        double x0[3] = { sys_state.xyz[i0].x, sys_state.xyz[i0].y, sys_state.xyz[i0].z };
		double x1[3] = { sys_state.xyz[i1].x, sys_state.xyz[i1].y, sys_state.xyz[i1].z };

        double s0[3];
        double s1[3];

        cart_to_fract(s0, x0, I);
        cart_to_fract(s1, x1, I);

        double d2_ref = distance_ref_mic27(G, s0, s1);

		vec3_t dx = { x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2] };
        md_util_min_image_vec3(&dx, 1, &sys_state.unitcell);
        double d2 = dx.x * dx.x + dx.y * dx.y + dx.z * dx.z;

		EXPECT_NEAR(d2_ref, d2, 1.0e-4);
    }

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        double radius = rnd() * 10.0;

        uint32_t ref_count = 0;
        uint32_t sa_count = 0;

#if 0
        // Do N^2 test as well to validate the reference implementation
        ref_count = do_brute_force_double(sys_state.xyz, sys.atom.count, radius, G, I, NULL, NULL);
        md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, (float)radius, spatial_acc_pair_count_callback, &sa_count);
        EXPECT_NEAR(ref_count, sa_count, 2);
#endif

        ref_count = 0;
        const double radius_sq = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys_state.xyz[i].x, sys_state.xyz[i].y, sys_state.xyz[i].z };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < radius_sq) {
                ref_count += 1;
            }
        }

        sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_from_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        query_points_in_sphere(&acc, x0, radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}

// ### CELL INDEX ###
// Pairs against a brute force reference, in boxes and point sets which stress the sparse cell index: clusters in a
// large box, a box of micrometres, subsets, partial periodicity, triclinic cells.

typedef struct cell_pairs_t {
    uint64_t* keys;
    size_t count;
    size_t cap;
} cell_pairs_t;

static void cell_pairs_push(cell_pairs_t* p, uint32_t i, uint32_t j) {
    if (p->count == p->cap) {
        p->cap  = p->cap ? p->cap * 2 : 1024;
        p->keys = (uint64_t*)realloc(p->keys, p->cap * sizeof(uint64_t));
    }
    p->keys[p->count++] = ((uint64_t)MIN(i, j) << 32) | MAX(i, j);
}

static void cell_pairs_callback(const uint32_t* i, const uint32_t* j, const float* d2, size_t n, void* user) {
    (void)d2;
    for (size_t k = 0; k < n; ++k) cell_pairs_push((cell_pairs_t*)user, i[k], j[k]);
}

static int cell_cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

static cell_pairs_t cell_collect_pairs(const md_spatial_acc_t* acc, double cutoff) {
    cell_pairs_t p = {0};
    md_spatial_acc_for_each_internal_pair_within_cutoff(acc, cutoff, cell_pairs_callback, &p);
    if (p.count) qsort(p.keys, p.count, sizeof(uint64_t), cell_cmp_u64);
    return p;
}

static bool cell_has_pair(const cell_pairs_t* p, uint64_t key) {
    return p->count && bsearch(&key, p->keys, p->count, sizeof(uint64_t), cell_cmp_u64) != NULL;
}

// Squared minimum image distance, in double, along the periodic axes of cell (NULL: none). The cutoffs used here stay
// below half the perpendicular width of the cell, so the nearest of the 27 images is the minimum image.
static double cell_ref_d2(const float* a, const float* b, const md_unitcell_t* cell) {
    double d[3] = { (double)b[0] - a[0], (double)b[1] - a[1], (double)b[2] - a[2] };
    if (!cell) return d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    double A[3][3], I[3][3];
    md_unitcell_A_extract_double(A, cell);
    md_unitcell_I_extract_double(I, cell);
    const int pbc[3] = { (cell->flags & MD_UNITCELL_PBC_X) != 0, (cell->flags & MD_UNITCELL_PBC_Y) != 0, (cell->flags & MD_UNITCELL_PBC_Z) != 0 };
    double s[3];
    for (int r = 0; r < 3; ++r) s[r] = I[0][r] * d[0] + I[1][r] * d[1] + I[2][r] * d[2];
    for (int k = 0; k < 3; ++k) if (pbc[k]) s[k] -= round(s[k]);
    double best = DBL_MAX;
    for (int iz = -1; iz <= 1; ++iz) for (int iy = -1; iy <= 1; ++iy) for (int ix = -1; ix <= 1; ++ix) {
        const int o[3] = { ix, iy, iz };
        if ((!pbc[0] && ix) || (!pbc[1] && iy) || (!pbc[2] && iz)) continue;
        double t[3] = { s[0] + o[0], s[1] + o[1], s[2] + o[2] };
        double c[3];
        for (int r = 0; r < 3; ++r) c[r] = A[0][r] * t[0] + A[1][r] * t[1] + A[2][r] * t[2];
        best = MIN(best, c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    }
    return best;
}

// The pairs of the structure against a brute force search over the stream. Pairs within 1e-3 A of the cutoff may go
// either way: the structure works in single precision.
static bool cell_matches_reference(const md_coord_stream_t* stream, const md_spatial_acc_t* acc, double cutoff, const md_unitcell_t* cell, bool stream_idx) {
    cell_pairs_t got = cell_collect_pairs(acc, cutoff);
    const double r_in  = (cutoff - 1.0e-3) * (cutoff - 1.0e-3);
    const double r_out = (cutoff + 1.0e-3) * (cutoff + 1.0e-3);
    size_t missing = 0, extra = 0, found = 0;
    for (size_t a = 0; a < stream->count; ++a) {
        const size_t ia = stream->idx ? (size_t)stream->idx[a] : a;
        const float pa[3] = { stream->soa.x[ia], stream->soa.y[ia], stream->soa.z[ia] };
        for (size_t b = a + 1; b < stream->count; ++b) {
            const size_t ib = stream->idx ? (size_t)stream->idx[b] : b;
            const float pb[3] = { stream->soa.x[ib], stream->soa.y[ib], stream->soa.z[ib] };
            const double d2 = cell_ref_d2(pa, pb, cell);
            if (d2 >= r_out) continue;
            const uint32_t ka = stream_idx ? (uint32_t)ia : (uint32_t)a;
            const uint32_t kb = stream_idx ? (uint32_t)ib : (uint32_t)b;
            const uint64_t key = ((uint64_t)MIN(ka, kb) << 32) | MAX(ka, kb);
            const bool has = cell_has_pair(&got, key);
            if (has) found += 1;
            else if (d2 < r_in) missing += 1;
        }
    }
    // Everything reported has to be within the cutoff (up to the tolerance), and reported once
    extra = got.count - found;
    for (size_t k = 1; k < got.count; ++k) if (got.keys[k] == got.keys[k - 1]) extra += 1;
    if (missing || extra) {
        printf("reference: %zu pairs missing, %zu reported which should not be (of %zu)\n", missing, extra, got.count);
    }
    free(got.keys);
    return missing == 0 && extra == 0;
}

UTEST(spatial_acc_cells, matches_reference_randomized) {
    srand(31337);
    const size_t max_count = 3000;
    float* x   = (float*)malloc(max_count * sizeof(float));
    float* y   = (float*)malloc(max_count * sizeof(float));
    float* z   = (float*)malloc(max_count * sizeof(float));
    int*   idx = (int*)malloc(max_count * sizeof(int));

    for (int iter = 0; iter < 60; ++iter) {
        const size_t count = 1 + (size_t)rnd_rng(0, (double)(max_count - 1));
        const double cutoff = rnd_rng(1.5, 9.5);
        // At least three cutoffs across, so the minimum image is the only image within the cutoff
        const double ext[3] = { rnd_rng(30, 200), rnd_rng(30, 200), rnd_rng(30, 400) };
        // Mostly cells of the cutoff, sometimes down to half of it (two cells of reach), sometimes chosen by the structure
        const int cell_choice = (int)rnd_rng(0, 3);
        const double cell_ext = cell_choice == 0 ? cutoff : cell_choice == 1 ? cutoff * rnd_rng(0.51, 1.0) : 0.0;
        const bool clustered = rnd_rng(0, 1) < 0.5;
        const bool tri = rnd_rng(0, 1) < 0.25;
        const double sk[3] = { tri ? rnd_rng(-0.2, 0.2) * ext[1] : 0, tri ? rnd_rng(-0.2, 0.2) * ext[2] : 0, tri ? rnd_rng(-0.2, 0.2) * ext[2] : 0 };
        for (size_t i = 0; i < count; ++i) {
            double f[3];
            if (clustered) {
                const int c = (int)rnd_rng(0, 4);
                f[0] = 0.1 + 0.2 * c + rnd_rng(0, 0.05);
                f[1] = rnd_rng(0, 1);
                f[2] = 0.95 + rnd_rng(0, 0.1);   // Across the boundary in z
            } else {
                f[0] = rnd_rng(-0.1, 1.1);       // Partly outside the cell
                f[1] = rnd_rng(0, 1);
                f[2] = rnd_rng(0, 1);
            }
            x[i] = (float)(f[0] * ext[0] + f[1] * sk[0] + f[2] * sk[1]);
            y[i] = (float)(f[1] * ext[1] + f[2] * sk[2]);
            z[i] = (float)(f[2] * ext[2]);
        }

        md_unitcell_t cell = md_unitcell_from_basis_parameters(ext[0], ext[1], ext[2], sk[0], sk[1], sk[2]);
        const int pbc = tri ? 7 : (int)rnd_rng(0, 8);
        if (!(pbc & 1)) cell.flags &= ~MD_UNITCELL_PBC_X;
        if (!(pbc & 2)) cell.flags &= ~MD_UNITCELL_PBC_Y;
        if (!(pbc & 4)) cell.flags &= ~MD_UNITCELL_PBC_Z;
        const md_unitcell_t* cell_ptr = (!tri && rnd_rng(0, 1) < 0.15) ? NULL : &cell;

        // A subset through an index, reported by stream index or by position
        const bool subset = rnd_rng(0, 1) < 0.3;
        size_t num = count;
        if (subset) {
            num = 0;
            for (size_t i = 0; i < count; ++i) {
                if (rnd_rng(0, 1) < 0.4) idx[num++] = (int)i;
            }
            if (num == 0) idx[num++] = 0;
        }
        const bool stream_idx = subset && rnd_rng(0, 1) < 0.5;

        md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, subset ? idx : NULL, num);
        md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
        md_spatial_acc_desc_t desc = { .coords = &stream, .cell_ext = cell_ext, .cutoff = cutoff, .unitcell = cell_ptr,
                                       .flags = stream_idx ? MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX : 0 };
        md_spatial_acc_init_desc(&acc, &desc);
        const bool ok = cell_matches_reference(&stream, &acc, cutoff, cell_ptr, stream_idx);
        if (!ok) {
            printf("iteration %d: %zu points, extent %.1f %.1f %.1f, %s, pbc %d, cutoff %.2f, cell %.2f\n", iter, num, ext[0], ext[1], ext[2], tri ? "triclinic" : "orthorhombic", cell_ptr ? pbc : -1, cutoff, cell_ext);
        }
        EXPECT_TRUE(ok);
        md_spatial_acc_free(&acc);
    }

    free(x);
    free(y);
    free(z);
    free(idx);
}

// A handful of clusters in a box of 3000 A with 5 A cells: as a dense grid 216 million cells, of which the index
// stores the few thousand occupied ones
UTEST(spatial_acc_cells, clustered_large_box) {
    srand(4711);
    const size_t count = 6000;
    float* x = (float*)malloc(count * sizeof(float));
    float* y = (float*)malloc(count * sizeof(float));
    float* z = (float*)malloc(count * sizeof(float));
    for (size_t i = 0; i < count; ++i) {
        const int c = (int)rnd_rng(0, 5);
        x[i] = (float)(100 + c * 600 + rnd_rng(0, 30));
        y[i] = (float)(2900 + rnd_rng(0, 200));   // Across the periodic boundary in y
        z[i] = (float)(c * 37 + rnd_rng(0, 30));
    }
    md_unitcell_t cell = md_unitcell_from_extent(3000, 3000, 3000);
    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, count);

    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 5.0, &cell, 0);
    EXPECT_EQ(600u, acc.cell_dim[0]);
    EXPECT_LE(acc.num_cells, count);
    EXPECT_TRUE(cell_matches_reference(&stream, &acc, 5.0, &cell, false));
    md_spatial_acc_free(&acc);

    free(x);
    free(y);
    free(z);
}

// Far more cells than a dense grid could hold: 2 x 2 x 8 um with 3 A cells
UTEST(spatial_acc_cells, micrometre_box) {
    srand(1234);
    const size_t count = 10000;
    float* x = (float*)malloc(count * sizeof(float));
    float* y = (float*)malloc(count * sizeof(float));
    float* z = (float*)malloc(count * sizeof(float));
    // Pairs of points 2 A apart, scattered over the box, plus the image of one across each periodic boundary
    for (size_t i = 0; i < count; i += 2) {
        x[i] = (float)rnd_rng(0, 20000);
        y[i] = (float)rnd_rng(0, 20000);
        z[i] = (float)rnd_rng(0, 80000);
        x[i + 1] = x[i] + 2.0f;
        y[i + 1] = y[i];
        z[i + 1] = z[i];
    }
    x[0] = 19999.0f; x[1] = 1.0f;       // 2 A apart through the boundary in x
    md_unitcell_t cell = md_unitcell_from_extent(20000, 20000, 80000);
    md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, count);

    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);
    EXPECT_EQ(6666u, acc.cell_dim[0]);
    EXPECT_EQ(26666u, acc.cell_dim[2]);
    EXPECT_LE(acc.num_cells, count);

    cell_pairs_t p = cell_collect_pairs(&acc, 2.5);
    // Every pair, and the odd chance encounter between two pairs
    EXPECT_GE(p.count, count / 2);
    EXPECT_TRUE(cell_has_pair(&p, ((uint64_t)0 << 32) | 1));
    free(p.keys);
    md_spatial_acc_free(&acc);

    free(x);
    free(y);
    free(z);
}

// The cell extent from the cutoff: the cutoff where the points are dense, larger (at most 4x) where they are sparse,
// and the same pairs either way
UTEST(spatial_acc_cells, cell_extent_from_cutoff) {
    srand(2718);
    const size_t count = 20000;
    float* x = (float*)malloc(count * sizeof(float));
    float* y = (float*)malloc(count * sizeof(float));
    float* z = (float*)malloc(count * sizeof(float));
    const double cutoff = 5.0;

    // Dense: ~0.1 points per A^3, about 12 per cell of the cutoff. Sparse: ~1 point per cell of the cutoff.
    const double extents[2] = { 58.0, 136.0 };
    for (int pass = 0; pass < 2; ++pass) {
        const double ext = extents[pass];
        for (size_t i = 0; i < count; ++i) {
            x[i] = (float)rnd_rng(0, ext);
            y[i] = (float)rnd_rng(0, ext);
            z[i] = (float)rnd_rng(0, ext);
        }
        md_unitcell_t cell = md_unitcell_from_extent(ext, ext, ext);
        md_coord_stream_t stream = md_coord_stream_from_soa(x, y, z, NULL, count);

        md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
        md_spatial_acc_desc_t desc = { .coords = &stream, .cutoff = cutoff, .unitcell = &cell };
        md_spatial_acc_init_desc(&acc, &desc);
        const double cell_ext = ext / acc.cell_dim[0];
        if (pass == 0) {
            EXPECT_LT(cell_ext, cutoff * 1.1);
        } else {
            EXPECT_GT(cell_ext, cutoff * 1.5);
        }
        EXPECT_GE(cell_ext, cutoff);
        EXPECT_LE(cell_ext, cutoff * 4.0 * 1.01);
        EXPECT_TRUE(cell_matches_reference(&stream, &acc, cutoff, &cell, false));
        md_spatial_acc_free(&acc);
    }

    free(x);
    free(y);
    free(z);
}
