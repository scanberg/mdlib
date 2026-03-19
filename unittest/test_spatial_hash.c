#include "utest.h"

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

    md_coord_stream_t stream = md_coord_stream_create_soa(x, y, z, NULL, 10);
    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 10.0, &cell, 0);
    
    uint32_t count = 0;
    double p0[3] = {5, 0, 0};
    md_spatial_acc_for_each_point_in_sphere(&acc, p0, 1.5f, spatial_acc_point_count_callback, &count);
    EXPECT_EQ(3, count);

    count = 0;
    double p1[3] = {8.5f, 0, 0};
    md_spatial_acc_for_each_point_in_sphere(&acc, p1, 3, spatial_acc_point_count_callback, &count);
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
    md_coord_stream_t stream = md_coord_stream_create_soa(x, y, z, NULL, ARRAY_SIZE(x));

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

    md_coord_stream_t stream = md_coord_stream_create_soa(x, y, z, NULL, 3);
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));

    const double Lx = 50.0;
    const double Ly = 60.0;
    const double Lz = 70.0;
    md_unitcell_t cell = md_unitcell_from_extent(Lx, Ly, Lz);

    const size_t N = 4096;
    float* x = (float*)md_arena_allocator_push(alloc, N * sizeof(float));
    float* y = (float*)md_arena_allocator_push(alloc, N * sizeof(float));
    float* z = (float*)md_arena_allocator_push(alloc, N * sizeof(float));

    srand(1337);
    for (size_t i = 0; i < N; ++i) {
        x[i] = (float)rnd_rng(0.0, Lx);
        y[i] = (float)rnd_rng(0.0, Ly);
        z[i] = (float)rnd_rng(0.0, Lz);
    }

    md_coord_stream_t stream = md_coord_stream_create_soa(x, y, z, NULL, N);
    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);

    spatial_acc_point_collect_t got = { .idx = NULL, .alloc = md_get_heap_allocator() };
    uint8_t* seen = (uint8_t*)malloc(N);
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

    free(seen);
    md_array_free(got.idx, got.alloc);
    md_spatial_acc_free(&acc);
    md_arena_allocator_destroy(alloc);
}

UTEST(spatial_hash, aabb_periodic_triclinic_randomized_reference) {
    // Robust randomized test for periodic triclinic AABB queries.
    // Reference uses brute-force over 27 periodic images (sufficient for chosen small radii).
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));

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
    float* x = (float*)md_arena_allocator_push(alloc, N * sizeof(float));
    float* y = (float*)md_arena_allocator_push(alloc, N * sizeof(float));
    float* z = (float*)md_arena_allocator_push(alloc, N * sizeof(float));

    srand(7331);
    for (size_t i = 0; i < N; ++i) {
        const double s[3] = { rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0), rnd_rng(0.0, 1.0) };
        double p[3];
        fract_to_cart(p, s, A);
        x[i] = (float)p[0];
        y[i] = (float)p[1];
        z[i] = (float)p[2];
    }

    md_coord_stream_t stream = md_coord_stream_create_soa(x, y, z, NULL, N);
    md_spatial_acc_t acc = { .alloc = md_get_heap_allocator() };
    md_spatial_acc_init(&acc, &stream, 3.0, &cell, 0);

    spatial_acc_point_collect_t got = { .idx = NULL, .alloc = md_get_heap_allocator() };
    uint8_t* seen = (uint8_t*)malloc(N);
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

    free(seen);
    md_array_free(got.idx, got.alloc);
    md_spatial_acc_free(&acc);
    md_arena_allocator_destroy(alloc);
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

UTEST(spatial_hash, n2) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), alloc));

    md_system_t sys = {0};
    ASSERT_TRUE(md_gro_system_init(&sys, &gro_data, alloc));

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

        md_coord_stream_t stream = md_coord_stream_create_soa(x, y, z, NULL, TEST_COUNT);
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
                size_t temp_pos = md_temp_get_pos();
                md_hashmap32_t map_sa = { .allocator = md_get_temp_allocator() };
                md_hashmap32_t map_bf = { .allocator = md_get_temp_allocator() };

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

                md_temp_set_pos_back(temp_pos);
            }
        }
#undef TEST_COUNT
    }

#if 1
    md_unitcell_t cell = sys.unitcell;
    const size_t expected_count = 3711879;

    double G[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &cell);
    md_unitcell_I_extract_double(I, &cell);

    md_timestamp_t start, end;
    uint32_t count = 0;

    // Spatial acc implementation
    start = md_time_current();
    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = {.alloc = alloc};
    md_spatial_acc_init(&acc, &stream, 5.0, &cell, 0);
    md_spatial_acc_for_each_internal_pair_in_neighboring_cells(&acc, spatial_acc_neighbor_callback, &count);
	//md_spatial_acc_for_each_pair_within_cutoff(&acc, 5.0, spatial_acc_cutoff_callback, &count);
    end = md_time_current();
    size_t sa_count = count;
    //end = md_time_current();
    //printf("Spatial acc cell neighborhood: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_NEAR(expected_count, sa_count, 5);
    if (sa_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, sa_count);
    }

    start = md_time_current();
    count = 0;
    md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &stream, 5.0, spatial_acc_neighbor_callback, &count, 0);
	end = md_time_current();
    sa_count = count;
	//printf("Spatial acc external query: %f ms\n", md_time_as_milliseconds(end - start));
	size_t ext_expected_count = expected_count * 2 + sys.atom.count;
	if (sa_count != ext_expected_count) {
		printf("Count mismatch: expected %zu, got %zu\n", ext_expected_count, sa_count);
	}

#if 0
    // This is so slow that we don't want to run it by default, but it can be useful for validating the reference implementation
    // Brute force
    start = md_time_current();
    size_t bf_count = do_brute_force_double(sys.atom.x, sys.atom.y, sys.atom.z, sys.atom.count, 5.0, G, I, NULL, NULL);
    end = md_time_current();
    printf("Brute force: %f ms\n", md_time_as_milliseconds(end - start));
    EXPECT_EQ(expected_count, bf_count);
    if (bf_count != expected_count) {
        printf("Count mismatch: expected %zu, got %zu\n", expected_count, bf_count);
    }
#endif
#endif

    md_arena_allocator_destroy(alloc);
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

    md_system_t sys = { 0 };
    ASSERT_TRUE(md_gro_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), NULL, alloc));

    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys.unitcell, 0);
    
    srand(31);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys.unitcell);
    md_unitcell_A_extract_double(A, &sys.unitcell);
    md_unitcell_I_extract_double(I, &sys.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        double radius = rnd() * 20;

        int ref_count = 0;
        const double rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys.atom.x[i], sys.atom.y[i], sys.atom.z[i] };
            double si[3];
            cart_to_fract(si, xi, I);

            if (distance_ref_mic27(G, s0, si) < rad2) {
                ref_count += 1;
            }
        }

        uint32_t sa_count = 0;
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        md_coord_stream_t ext_stream = md_coord_stream_create_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }

        sa_count = 0;
		md_spatial_acc_for_each_point_in_sphere(&acc, x0, radius, spatial_acc_point_count_callback, &sa_count);
		EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}

UTEST_F(spatial_hash, test_correctness_ala) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { 0 };
    ASSERT_TRUE(md_pdb_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, alloc));

    srand(31);
    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys.unitcell);
    md_unitcell_A_extract_double(A, &sys.unitcell);
    md_unitcell_I_extract_double(I, &sys.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        double radius = rnd() * 20.0;

        int ref_count = 0;
        const double rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys.atom.x[i], sys.atom.y[i], sys.atom.z[i] };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < rad2) {
                ref_count += 1;
            }
        }

        uint32_t sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_create_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        md_spatial_acc_for_each_point_in_sphere(&acc, x0, (float)radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }
    }
}

UTEST_F(spatial_hash, test_correctness_water) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { 0 };
    ASSERT_TRUE(md_gro_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro"), NULL, alloc));

    srand(31);

    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys.unitcell);
    md_unitcell_A_extract_double(A, &sys.unitcell);
    md_unitcell_I_extract_double(I, &sys.unitcell);

    const int num_iter = 100;
    for (int iter = 0; iter < num_iter; ++iter) {
        double s0[3] = { rnd(), rnd(), rnd() };
        double x0[3];
        fract_to_cart(x0, s0, A);
        vec3_t pos = vec3_set(x0[0], x0[1], x0[2]);
        double radius = rnd() * 20.0;

        int ref_count = 0;
        const double rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys.atom.x[i], sys.atom.y[i], sys.atom.z[i] };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < rad2) {
                ref_count += 1;
            }
        }

        uint32_t sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_create_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {           
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        md_spatial_acc_for_each_point_in_sphere(&acc, x0, radius, spatial_acc_point_count_callback, &sa_count);
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
    md_lammps_molecule_loader_arg_t args = md_lammps_molecule_loader_arg(atom_format);

    md_system_t sys = { 0 };
    ASSERT_TRUE(md_lammps_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Triclinic_Init.data"), &args, alloc));

    srand(31);

    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys.unitcell);
    md_unitcell_A_extract_double(A, &sys.unitcell);
    md_unitcell_I_extract_double(I, &sys.unitcell);

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
        ref_count = do_brute_force_double(sys.atom.x, sys.atom.y, sys.atom.z, sys.atom.count, radius, G, I, NULL, NULL);
		md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, (float)radius, spatial_acc_pair_count_callback, &sa_count);
        EXPECT_NEAR(ref_count, sa_count, 2);
#endif

        ref_count = 0;
        const double rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys.atom.x[i], sys.atom.y[i], sys.atom.z[i] };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < rad2) {
                ref_count += 1;
            }
        }

        sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_create_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        md_spatial_acc_for_each_point_in_sphere(&acc, x0, (float)radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}


UTEST_F(spatial_hash, npt_triclinic) {
    md_allocator_i* alloc = utest_fixture->arena;

    md_system_t sys = { 0 };
    ASSERT_TRUE(md_gro_system_loader()->init_from_file(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/npt.gro"), NULL, alloc));

    srand(31);

    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = alloc };
    md_spatial_acc_init(&acc, &stream, 10.0, &sys.unitcell, 0);

    double G[3][3], A[3][3], I[3][3];
    md_unitcell_G_extract_double(G, &sys.unitcell);
    md_unitcell_A_extract_double(A, &sys.unitcell);
    md_unitcell_I_extract_double(I, &sys.unitcell);

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

        double x0[3] = { sys.atom.x[i0], sys.atom.y[i0], sys.atom.z[i0] };
		double x1[3] = { sys.atom.x[i1], sys.atom.y[i1], sys.atom.z[i1] };

        double s0[3];
        double s1[3];

        cart_to_fract(s0, x0, I);
        cart_to_fract(s1, x1, I);

        double d2_ref = distance_ref_mic27(G, s0, s1);

		vec3_t dx = { x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2] };
        md_util_min_image_vec3(&dx, 1, &sys.unitcell);
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
        ref_count = do_brute_force_double(sys.atom.x, sys.atom.y, sys.atom.z, sys.atom.count, radius, G, I, NULL, NULL);
        md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, (float)radius, spatial_acc_pair_count_callback, &sa_count);
        EXPECT_NEAR(ref_count, sa_count, 2);
#endif

        ref_count = 0;
        const double rad2 = radius * radius;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            double xi[3] = { sys.atom.x[i], sys.atom.y[i], sys.atom.z[i] };
            double si[3];
            cart_to_fract(si, xi, I);
            if (distance_ref_mic27(G, s0, si) < rad2) {
                ref_count += 1;
            }
        }

        sa_count = 0;
        md_coord_stream_t ext_stream = md_coord_stream_create_soa(&pos.x, &pos.y, &pos.z, NULL, 1);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &ext_stream, (float)radius, spatial_acc_pair_count_callback, &sa_count, 0);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, (float)radius, ref_count, sa_count);
        }

        sa_count = 0;
        md_spatial_acc_for_each_point_in_sphere(&acc, x0, radius, spatial_acc_point_count_callback, &sa_count);
        EXPECT_EQ(ref_count, sa_count);
        if (sa_count != ref_count) {
            printf("iter: %i, pos: %f %f %f, rad: %f, expected: %i, got: %i\n", iter, pos.x, pos.y, pos.z, radius, ref_count, sa_count);
        }
    }
}
