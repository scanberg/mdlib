#include "utest.h"

#include <core/md_sort.h>
#include <core/md_unitcell.h>
#include <core/md_vec_math.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

#include <string.h>
#include <stdio.h>

// splitmix64: the same keys on every platform (unlike rand), with all 64 bits random
static uint64_t radix_test_rand(uint64_t* state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// Sizes around the insertion sort of small arrays (64), the 256 buckets of a pass, and one past the 16 bits of an index
static const size_t radix_test_sizes[] = { 0, 1, 2, 3, 16, 63, 64, 65, 255, 256, 257, 1000, 65537 };

#define RADIX_TEST_MAX_SIZE 65537

// The distributions of keys: the passes the sorts skip (an odd number of passes left leaves the keys in the scratch),
// long runs of equal keys for the stability of the indices, keys in order (left as they are) and in order but the last,
// and the two ends of the range
enum { RADIX_TEST_DIST_COUNT = 11 };

static const char* radix_test_dist_name[RADIX_TEST_DIST_COUNT] = {
    "random", "few distinct", "all equal", "ascending", "descending", "top digit only", "three low digits",
    "packed index pairs", "first and last digit", "both ends of the range", "ascending but the last",
};

static uint64_t radix_test_key64(int dist, size_t i, size_t n, uint64_t* rng) {
    const uint64_t r = radix_test_rand(rng);
    switch (dist) {
    case 0: return r;
    case 1: return r % 7;
    case 2: return 0xDEADBEEFCAFEF00DULL;
    case 3: return (uint64_t)i;
    case 4: return UINT64_MAX - (uint64_t)i;
    case 5: return (r << 56) | 0x0000123456789ABCULL;
    case 6: return (r & 0xFFFFFF) | 0xAB00000000000000ULL;
    case 7: return ((r % 5000) << 32) | ((r >> 32) % 5000);
    case 8: return (r & 0xFF000000000000FFULL) | 0x0011223344556600ULL;
    case 9: return (r & 1) ? UINT64_MAX - (r >> 1) % 1000 : (r >> 1) % 1000;
    case 10: return i + 1 == n ? 0 : (uint64_t)i + 1;
    default: return 0;
    }
}

static uint32_t radix_test_key32(int dist, size_t i, size_t n, uint64_t* rng) {
    const uint64_t r = radix_test_rand(rng);
    switch (dist) {
    case 0: return (uint32_t)r;
    case 1: return (uint32_t)(r % 7);
    case 2: return 0xCAFEF00Du;
    case 3: return (uint32_t)i;
    case 4: return UINT32_MAX - (uint32_t)i;
    case 5: return (uint32_t)(r << 24) | 0x123456u;
    case 6: return (uint32_t)(r & 0xFFFFFF) | 0xAB000000u;
    case 7: return (uint32_t)(((r % 5000) << 16) | ((r >> 32) % 5000));
    case 8: return (uint32_t)(r & 0xFF0000FFu) | 0x00112200u;
    case 9: return (r & 1) ? UINT32_MAX - (uint32_t)((r >> 1) % 1000) : (uint32_t)((r >> 1) % 1000);
    case 10: return i + 1 == n ? 0 : (uint32_t)i + 1;
    default: return 0;
    }
}

typedef struct radix_test_ref_t {
    uint64_t key;
    uint32_t idx;
} radix_test_ref_t;

static int radix_test_ref_cmp(const void* a, const void* b) {
    const radix_test_ref_t* x = (const radix_test_ref_t*)a;
    const radix_test_ref_t* y = (const radix_test_ref_t*)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    if (x->idx != y->idx) return x->idx < y->idx ? -1 : 1;
    return 0;
}

// The reference, from qsort: the keys in order, and the indices of a stable sort (equal keys by index)
static void radix_test_reference(uint64_t* out_sorted, uint32_t* out_idx, radix_test_ref_t* scratch, const uint64_t* keys, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        scratch[i] = (radix_test_ref_t){ keys[i], (uint32_t)i };
    }
    qsort(scratch, n, sizeof(radix_test_ref_t), radix_test_ref_cmp);
    for (size_t i = 0; i < n; ++i) {
        out_sorted[i] = scratch[i].key;
        out_idx[i] = scratch[i].idx;
    }
}

#define RADIX_TEST_GUARD 0xA5A5A5A5u

// Every size and distribution against the reference: the keys sorted in place, and the indices exactly those of a stable
// sort. A guard past the end of each output catches a write beyond count.
UTEST(sort, radix_sort_u32) {
    md_temp_scope_t temp = md_temp_begin();
    const size_t cap = RADIX_TEST_MAX_SIZE + 1;
    uint32_t* keys       = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t* sorted     = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t* idx        = md_temp_alloc_array(temp, uint32_t, cap);
    uint64_t* keys_wide  = md_temp_alloc_array(temp, uint64_t, cap);
    uint64_t* ref_sorted = md_temp_alloc_array(temp, uint64_t, cap);
    uint32_t* ref_idx    = md_temp_alloc_array(temp, uint32_t, cap);
    radix_test_ref_t* scratch = md_temp_alloc_array(temp, radix_test_ref_t, cap);

    uint64_t rng = 0x5EED32;
    for (size_t s = 0; s < MD_COUNTOF(radix_test_sizes); ++s) {
        const size_t n = radix_test_sizes[s];
        for (int d = 0; d < RADIX_TEST_DIST_COUNT; ++d) {
            char msg[128];
            snprintf(msg, sizeof(msg), "%s, n = %zu", radix_test_dist_name[d], n);

            for (size_t i = 0; i < n; ++i) {
                keys[i] = radix_test_key32(d, i, n, &rng);
                keys_wide[i] = keys[i];
            }
            radix_test_reference(ref_sorted, ref_idx, scratch, keys_wide, n);

            memcpy(sorted, keys, n * sizeof(uint32_t));
            sorted[n] = RADIX_TEST_GUARD;
            md_sort_radix_u32(sorted, n);
            bool same = true;
            for (size_t i = 0; i < n; ++i) {
                same &= sorted[i] == (uint32_t)ref_sorted[i];
            }
            EXPECT_TRUE_MSG(same, msg);
            EXPECT_EQ(RADIX_TEST_GUARD, sorted[n]);

            idx[n] = RADIX_TEST_GUARD;
            md_argsort_radix_u32(idx, keys, n);
            EXPECT_TRUE_MSG(n == 0 || memcmp(idx, ref_idx, n * sizeof(uint32_t)) == 0, msg);
            EXPECT_EQ(RADIX_TEST_GUARD, idx[n]);
        }
    }

    md_temp_end(temp);
}

UTEST(sort, radix_sort_u64) {
    md_temp_scope_t temp = md_temp_begin();
    const size_t cap = RADIX_TEST_MAX_SIZE + 1;
    uint64_t* keys       = md_temp_alloc_array(temp, uint64_t, cap);
    uint64_t* sorted     = md_temp_alloc_array(temp, uint64_t, cap);
    uint32_t* idx        = md_temp_alloc_array(temp, uint32_t, cap);
    uint64_t* ref_sorted = md_temp_alloc_array(temp, uint64_t, cap);
    uint32_t* ref_idx    = md_temp_alloc_array(temp, uint32_t, cap);
    radix_test_ref_t* scratch = md_temp_alloc_array(temp, radix_test_ref_t, cap);

    uint64_t rng = 0x5EED64;
    for (size_t s = 0; s < MD_COUNTOF(radix_test_sizes); ++s) {
        const size_t n = radix_test_sizes[s];
        for (int d = 0; d < RADIX_TEST_DIST_COUNT; ++d) {
            char msg[128];
            snprintf(msg, sizeof(msg), "%s, n = %zu", radix_test_dist_name[d], n);

            for (size_t i = 0; i < n; ++i) {
                keys[i] = radix_test_key64(d, i, n, &rng);
            }
            radix_test_reference(ref_sorted, ref_idx, scratch, keys, n);

            memcpy(sorted, keys, n * sizeof(uint64_t));
            sorted[n] = RADIX_TEST_GUARD;
            md_sort_radix_u64(sorted, n);
            EXPECT_TRUE_MSG(n == 0 || memcmp(sorted, ref_sorted, n * sizeof(uint64_t)) == 0, msg);
            EXPECT_EQ((uint64_t)RADIX_TEST_GUARD, sorted[n]);

            idx[n] = RADIX_TEST_GUARD;
            md_argsort_radix_u64(idx, keys, n);
            EXPECT_TRUE_MSG(n == 0 || memcmp(idx, ref_idx, n * sizeof(uint32_t)) == 0, msg);
            EXPECT_EQ(RADIX_TEST_GUARD, idx[n]);
        }
    }

    md_temp_end(temp);
}

// Timings on large arrays of random keys, each result checked for order
UTEST(sort, radix_sort) {
    uint32_t arr[] = { 1, 278, 128312745, 4, 5, 0, 12382, 26, 12, 14, 7 };
    size_t len = MD_COUNTOF(arr);

    uint32_t idx[16];

    md_argsort_radix_u32(idx, arr, len);
    for (size_t i = 0; i < len - 1; ++i) {
        EXPECT_LE(arr[idx[i]], arr[idx[i+1]]);
    }

    md_sort_radix_u32(arr, len);

    for (size_t i = 0; i < len - 1; ++i) {
    	EXPECT_LE(arr[i], arr[i+1]);
    }

    md_temp_scope_t temp = md_temp_begin();

    size_t N = 10000000;
    uint32_t* values   = md_temp_alloc_array(temp, uint32_t, N);
    uint32_t* indices  = md_temp_alloc_array(temp, uint32_t, N);

    for (size_t i = 0; i < N; ++i) {
        values[i] = (uint32_t)rand() * rand();
    }

    md_tick_t t0, t1;
    bool ordered;

    t0 = md_tick_now();
    md_argsort_radix_u32(indices, values, N);
    t1 = md_tick_now();

    printf("Time for radix argsort (u32): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values[indices[i - 1]] <= values[indices[i]];
    EXPECT_TRUE(ordered);

    t0 = md_tick_now();
    md_sort_radix_u32(values, N);
    t1 = md_tick_now();

    printf("Time for radix sort (u32): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values[i - 1] <= values[i];
    EXPECT_TRUE(ordered);

    uint64_t* values64 = md_temp_alloc_array(temp, uint64_t, N);
    uint64_t rng = 0x5EED;
    for (size_t i = 0; i < N; ++i) {
        values64[i] = radix_test_rand(&rng);
    }

    t0 = md_tick_now();
    md_argsort_radix_u64(indices, values64, N);
    t1 = md_tick_now();

    printf("Time for radix argsort (u64): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values64[indices[i - 1]] <= values64[indices[i]];
    EXPECT_TRUE(ordered);

    t0 = md_tick_now();
    md_sort_radix_u64(values64, N);
    t1 = md_tick_now();

    printf("Time for radix sort (u64): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values64[i - 1] <= values64[i];
    EXPECT_TRUE(ordered);

    md_temp_end(temp);
}

// The order is a permutation of the atoms, and atoms that share a position end up together.
UTEST(sort, sort_spatial) {
    enum { N = 64 };
    md_float3_t xyz[N];
    for (int i = 0; i < N; ++i) {
        const int cluster = i % 4;
        xyz[i] = md_float3_set(cluster * 20.0f, cluster * 20.0f, 0.0f);
    }
    uint32_t order[N];
    md_sort_spatial(order, (const float*)xyz, sizeof(md_float3_t), N);
    bool seen[N] = {0};
    for (int i = 0; i < N; ++i) {
        ASSERT_LT(order[i], (uint32_t)N);
        EXPECT_FALSE(seen[order[i]]);
        seen[order[i]] = true;
    }
    int changes = 0;
    for (int i = 1; i < N; ++i) {
        changes += (order[i] % 4) != (order[i - 1] % 4);
    }
    EXPECT_EQ(3, changes);
}
