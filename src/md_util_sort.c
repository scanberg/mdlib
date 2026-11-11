#include <md_util_sort.h>
#include <md_util_geom.h>

#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>
#include <core/md_simd.h>

#include <string.h>
#include <float.h>

// ### RADIX SORT ###

// Least significant digit first, 8 bits a digit. The histograms of every digit come from one read of the keys, and a pass
// whose digit all keys share would leave the order as it is, so it is skipped: keys with constant upper bytes (indices,
// pairs of them packed) take only the passes they need. Up to RADIX_SMALL_COUNT keys an insertion sort is faster than
// the histograms, and keys already in order are left as they are.
//
// The histograms have a row per bucket and a column per digit, so the prefix sums of every digit are one vector sum per
// row. The argsorts move each key together with its index, packed in one value, so a pass reads them in order and writes
// one stream per bucket, rather than looking up the key of each index.

#define RADIX_SMALL_COUNT 64

static void insertion_sort_u32(uint32_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const uint32_t k = data[i];
        size_t j = i;
        while (j > 0 && data[j - 1] > k) {
            data[j] = data[j - 1];
            --j;
        }
        data[j] = k;
    }
}

static void insertion_sort_u64(uint64_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const uint64_t k = data[i];
        size_t j = i;
        while (j > 0 && data[j - 1] > k) {
            data[j] = data[j - 1];
            --j;
        }
        data[j] = k;
    }
}

// The key above the index: the values are unique, in the order of a stable sort
static void insertion_argsort_u32(uint32_t* out_indices, const uint32_t* keys, size_t count) {
    MD_ASSERT(count <= RADIX_SMALL_COUNT);
    uint64_t v[RADIX_SMALL_COUNT];
    for (size_t i = 0; i < count; ++i) {
        v[i] = ((uint64_t)keys[i] << 32) | i;
    }
    insertion_sort_u64(v, count);
    for (size_t i = 0; i < count; ++i) {
        out_indices[i] = (uint32_t)v[i];
    }
}

static void insertion_argsort_u64(uint32_t* out_indices, const uint64_t* keys, size_t count) {
    MD_ASSERT(count <= RADIX_SMALL_COUNT);
    uint64_t k[RADIX_SMALL_COUNT];
    for (size_t i = 0; i < count; ++i) {
        const uint64_t key = keys[i];
        size_t j = i;
        // Only past greater keys: equal keys keep their order
        while (j > 0 && k[j - 1] > key) {
            k[j] = k[j - 1];
            out_indices[j] = out_indices[j - 1];
            --j;
        }
        k[j] = key;
        out_indices[j] = (uint32_t)i;
    }
}

static bool is_sorted_u32(const uint32_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        if (data[i] < data[i - 1]) return false;
    }
    return true;
}

static bool is_sorted_u64(const uint64_t* data, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        if (data[i] < data[i - 1]) return false;
    }
    return true;
}

static void radix_hist_u32(uint32_t hist[256][4], const uint32_t* keys, size_t count) {
    md_memset(hist, 0, sizeof(uint32_t) * 256 * 4);
    for (size_t i = 0; i < count; ++i) {
        const uint32_t k = keys[i];
        hist[(k >>  0) & 255][0]++;
        hist[(k >>  8) & 255][1]++;
        hist[(k >> 16) & 255][2]++;
        hist[(k >> 24) & 255][3]++;
    }
}

static void radix_hist_u64(uint32_t hist[256][8], const uint64_t* keys, size_t count) {
    md_memset(hist, 0, sizeof(uint32_t) * 256 * 8);
    for (size_t i = 0; i < count; ++i) {
        const uint64_t k = keys[i];
        hist[(k >>  0) & 255][0]++;
        hist[(k >>  8) & 255][1]++;
        hist[(k >> 16) & 255][2]++;
        hist[(k >> 24) & 255][3]++;
        hist[(k >> 32) & 255][4]++;
        hist[(k >> 40) & 255][5]++;
        hist[(k >> 48) & 255][6]++;
        hist[(k >> 56) & 255][7]++;
    }
}

// The digits whose pass reorders anything, then the counts turned into the offsets of the buckets. Any key tells whether
// a digit is shared by all: then its bucket holds every key.
static int radix_offsets_u32(int passes[4], uint32_t hist[256][4], uint32_t any_key, size_t count) {
    int num_passes = 0;
    for (int p = 0; p < 4; ++p) {
        if (hist[(any_key >> (p * 8)) & 255][p] != count) passes[num_passes++] = p;
    }
    md_128i sum = md_mm_setzero_si128();
    for (size_t i = 0; i < 256; ++i) {
        const md_128i c = md_mm_loadu_si128(hist[i]);
        md_mm_storeu_si128(hist[i], sum);
        sum = md_mm_add_epi32(sum, c);
    }
    return num_passes;
}

static int radix_offsets_u64(int passes[8], uint32_t hist[256][8], uint64_t any_key, size_t count) {
    int num_passes = 0;
    for (int p = 0; p < 8; ++p) {
        if (hist[(any_key >> (p * 8)) & 255][p] != count) passes[num_passes++] = p;
    }
    md_128i sum_lo = md_mm_setzero_si128();
    md_128i sum_hi = md_mm_setzero_si128();
    for (size_t i = 0; i < 256; ++i) {
        const md_128i lo = md_mm_loadu_si128(hist[i] + 0);
        const md_128i hi = md_mm_loadu_si128(hist[i] + 4);
        md_mm_storeu_si128(hist[i] + 0, sum_lo);
        md_mm_storeu_si128(hist[i] + 4, sum_hi);
        sum_lo = md_mm_add_epi32(sum_lo, lo);
        sum_hi = md_mm_add_epi32(sum_hi, hi);
    }
    return num_passes;
}

// The offsets of one digit, out of its column into a table of their own, for the pass to keep close
static inline void radix_column(uint32_t offset[256], const uint32_t* hist, int num_digits, int digit) {
    for (int b = 0; b < 256; ++b) {
        offset[b] = hist[b * num_digits + digit];
    }
}

static inline void radix_pass_u32(uint32_t* dst, const uint32_t* src, size_t count, uint32_t hist[256][4], int digit) {
    const int shift = digit * 8;
    for (size_t i = 0; i < count; ++i) {
        const uint32_t k = src[i];
        dst[hist[(k >> shift) & 255][digit]++] = k;
    }
}

static void sort_radix_u32(uint32_t* data, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][4];
    radix_hist_u32(hist, data, count);
    int passes[4];
    const int num_passes = radix_offsets_u32(passes, hist, data[0], count);
    if (num_passes == 0) return;

    uint32_t* src = data;
    uint32_t* dst = md_temp_alloc_array(temp, uint32_t, count);
    for (int p = 0; p < num_passes; ++p) {
        // The digit a constant in each: these loops are the time the sort takes
        switch (passes[p]) {
        case 0: radix_pass_u32(dst, src, count, hist, 0); break;
        case 1: radix_pass_u32(dst, src, count, hist, 1); break;
        case 2: radix_pass_u32(dst, src, count, hist, 2); break;
        case 3: radix_pass_u32(dst, src, count, hist, 3); break;
        default: MD_ASSERT(false); break;
        }
        uint32_t* tmp = src;
        src = dst;
        dst = tmp;
    }

    // An odd number of passes leaves the keys in the scratch
    if (src != data) {
        md_memcpy(data, src, count * sizeof(uint32_t));
    }
}

static void sort_radix_u64(uint64_t* data, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][8];
    radix_hist_u64(hist, data, count);
    int passes[8];
    const int num_passes = radix_offsets_u64(passes, hist, data[0], count);
    if (num_passes == 0) return;

    uint64_t* src = data;
    uint64_t* dst = md_temp_alloc_array(temp, uint64_t, count);
    uint32_t offset[256];
    for (int p = 0; p < num_passes; ++p) {
        const int shift = passes[p] * 8;
        radix_column(offset, hist[0], 8, passes[p]);
        for (size_t i = 0; i < count; ++i) {
            const uint64_t k = src[i];
            dst[offset[(k >> shift) & 255]++] = k;
        }
        uint64_t* tmp = src;
        src = dst;
        dst = tmp;
    }

    if (src != data) {
        md_memcpy(data, src, count * sizeof(uint64_t));
    }
}

// Each key with its index below it in one 64-bit value: the passes sort on the upper half, and the last keeps the index
static void argsort_radix_u32(uint32_t* out_indices, const uint32_t* keys, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][4];
    radix_hist_u32(hist, keys, count);
    int passes[4];
    const int num_passes = radix_offsets_u32(passes, hist, keys[0], count);
    uint32_t offset[256];

    if (num_passes == 0) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }
    if (num_passes == 1) {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 4, passes[0]);
        for (size_t i = 0; i < count; ++i) {
            out_indices[offset[(keys[i] >> shift) & 255]++] = (uint32_t)i;
        }
        return;
    }

    uint64_t* buf[2] = {
        md_temp_alloc_array(temp, uint64_t, count),
        num_passes > 2 ? md_temp_alloc_array(temp, uint64_t, count) : NULL,
    };

    {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 4, passes[0]);
        uint64_t* dst = buf[0];
        for (size_t i = 0; i < count; ++i) {
            const uint32_t k = keys[i];
            dst[offset[(k >> shift) & 255]++] = ((uint64_t)k << 32) | i;
        }
    }
    for (int p = 1; p < num_passes - 1; ++p) {
        const int shift = 32 + passes[p] * 8;
        radix_column(offset, hist[0], 4, passes[p]);
        const uint64_t* src = buf[(p - 1) & 1];
        uint64_t* dst = buf[p & 1];
        for (size_t i = 0; i < count; ++i) {
            const uint64_t v = src[i];
            dst[offset[(v >> shift) & 255]++] = v;
        }
    }
    {
        const int shift = 32 + passes[num_passes - 1] * 8;
        radix_column(offset, hist[0], 4, passes[num_passes - 1]);
        const uint64_t* src = buf[(num_passes - 2) & 1];
        for (size_t i = 0; i < count; ++i) {
            const uint64_t v = src[i];
            out_indices[offset[(v >> shift) & 255]++] = (uint32_t)v;
        }
    }
}

// A 64-bit key leaves no room for its index: the two go together as a record of 16 bytes
typedef struct radix_record_u64_t {
    uint64_t key;
    uint64_t idx;
} radix_record_u64_t;

static void argsort_radix_u64(uint32_t* out_indices, const uint64_t* keys, size_t count, md_temp_scope_t temp) {
    uint32_t hist[256][8];
    radix_hist_u64(hist, keys, count);
    int passes[8];
    const int num_passes = radix_offsets_u64(passes, hist, keys[0], count);
    uint32_t offset[256];

    if (num_passes == 0) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }
    if (num_passes == 1) {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 8, passes[0]);
        for (size_t i = 0; i < count; ++i) {
            out_indices[offset[(keys[i] >> shift) & 255]++] = (uint32_t)i;
        }
        return;
    }

    radix_record_u64_t* buf[2] = {
        md_temp_alloc_array(temp, radix_record_u64_t, count),
        num_passes > 2 ? md_temp_alloc_array(temp, radix_record_u64_t, count) : NULL,
    };

    {
        const int shift = passes[0] * 8;
        radix_column(offset, hist[0], 8, passes[0]);
        radix_record_u64_t* dst = buf[0];
        for (size_t i = 0; i < count; ++i) {
            const uint64_t k = keys[i];
            dst[offset[(k >> shift) & 255]++] = (radix_record_u64_t){ k, i };
        }
    }
    for (int p = 1; p < num_passes - 1; ++p) {
        const int shift = passes[p] * 8;
        radix_column(offset, hist[0], 8, passes[p]);
        const radix_record_u64_t* src = buf[(p - 1) & 1];
        radix_record_u64_t* dst = buf[p & 1];
        for (size_t i = 0; i < count; ++i) {
            const radix_record_u64_t r = src[i];
            dst[offset[(r.key >> shift) & 255]++] = r;
        }
    }
    {
        const int shift = passes[num_passes - 1] * 8;
        radix_column(offset, hist[0], 8, passes[num_passes - 1]);
        const radix_record_u64_t* src = buf[(num_passes - 2) & 1];
        for (size_t i = 0; i < count; ++i) {
            out_indices[offset[(src[i].key >> shift) & 255]++] = (uint32_t)src[i].idx;
        }
    }
}

/*
This is blatantly stolen and modified from 'Arseny Kapoulkines' goldnugget 'Mesh optimizer'
https://github.com/zeux/meshoptimizer/

MIT License
Copyright (c) 2016-2022 Arseny Kapoulkine
*/

// "Insert" two 0 bits after each of the 10 low bits of x
static inline uint32_t part1By2(uint32_t x) {
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x << 8))  & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x << 4))  & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x << 2))  & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return x;
}

static void computeOrder(uint32_t* result, const float* vx, const float* vy, const float* vz, size_t count, size_t byte_stride) {
    float minv[3] = {+FLT_MAX, +FLT_MAX, +FLT_MAX};
    float maxv[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    byte_stride = MD_MAX(byte_stride, sizeof(float));
    for (size_t i = 0; i < count; ++i) {
        float x = *(float*)((const char*)vx + byte_stride * i);
        float y = *(float*)((const char*)vy + byte_stride * i);
        float z = *(float*)((const char*)vz + byte_stride * i);

        minv[0] = MD_MIN(minv[0], x);
        maxv[0] = MD_MAX(maxv[0], x);

        minv[1] = MD_MIN(minv[1], y);
        maxv[1] = MD_MAX(maxv[1], y);

        minv[2] = MD_MIN(minv[2], z);
        maxv[2] = MD_MAX(maxv[2], z);
    }

    float extent = 0.f;

    extent = (maxv[0] - minv[0]) < extent ? extent : (maxv[0] - minv[0]);
    extent = (maxv[1] - minv[1]) < extent ? extent : (maxv[1] - minv[1]);
    extent = (maxv[2] - minv[2]) < extent ? extent : (maxv[2] - minv[2]);

    float scale = extent == 0 ? 0.f : 1.f / extent;

    // generate Morton order based on the position inside a unit cube
    for (size_t i = 0; i < count; ++i) {
        const float px = *(const float*)((const char*)vx + byte_stride * i);
        const float py = *(const float*)((const char*)vy + byte_stride * i);
        const float pz = *(const float*)((const char*)vz + byte_stride * i);
        int x = (int)((px - minv[0]) * scale * 1023.f + 0.5f);
        int y = (int)((py - minv[1]) * scale * 1023.f + 0.5f);
        int z = (int)((pz - minv[2]) * scale * 1023.f + 0.5f);

        result[i] = part1By2(x) | (part1By2(y) << 1) | (part1By2(z) << 2);
    }
}

static void computeHistogram10(uint32_t hist[1024][3], const uint32_t* data, size_t count) {
    // compute 3 10-bit histograms in parallel
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = data[i];

        hist[(id >> 0)  & 1023][0]++;
        hist[(id >> 10) & 1023][1]++;
        hist[(id >> 20) & 1023][2]++;
    }

    uint32_t sumx = 0, sumy = 0, sumz = 0;

    // replace histogram data with prefix histogram sums in-place
    for (int i = 0; i < 1024; ++i) {
        uint32_t hx = hist[i][0], hy = hist[i][1], hz = hist[i][2];

        hist[i][0] = sumx;
        hist[i][1] = sumy;
        hist[i][2] = sumz;

        sumx += hx;
        sumy += hy;
        sumz += hz;
    }

    MD_ASSERT(sumx == count && sumy == count && sumz == count);
}

static void radixPass10(uint32_t* destination, const uint32_t* source, const uint32_t* keys, size_t count, uint32_t hist[1024][3], int pass) {
    int bitoff = pass * 10;
    for (size_t i = 0; i < count; ++i) {
        uint32_t id = (keys[source[i]] >> bitoff) & 1023;
        destination[hist[id][pass]++] = source[i];
    }
}

// Packed coordinates: the bounds come from the packed AABB, which reads them as they lie, and the
// keys from each atom in turn.
void md_util_sort_spatial(uint32_t* source, const md_f32x3_t* xyz, size_t count) {
    if (!source || !xyz || count == 0) return;

    md_temp_scope_t temp_scope = md_temp_begin();

    float minv[3], maxv[3];
    md_util_aabb_compute(minv, maxv, xyz, NULL, NULL, count);
    const float extent = MD_MAX(maxv[0] - minv[0], MD_MAX(maxv[1] - minv[1], maxv[2] - minv[2]));
    const float scale = extent == 0 ? 0.f : 1023.f / extent;

    // Morton order based on the position inside a unit cube
    uint32_t* keys = md_temp_alloc_array(temp_scope, uint32_t, count);
    for (size_t i = 0; i < count; ++i) {
        const int x = (int)((xyz[i].x - minv[0]) * scale + 0.5f);
        const int y = (int)((xyz[i].y - minv[1]) * scale + 0.5f);
        const int z = (int)((xyz[i].z - minv[2]) * scale + 0.5f);
        keys[i] = part1By2(x) | (part1By2(y) << 1) | (part1By2(z) << 2);
    }

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram10(hist, keys, count);

    uint32_t* scratch = md_temp_alloc_array(temp_scope, uint32_t, count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass 10-bit radix sort computes the resulting order into scratch
    radixPass10(source, scratch, keys, count, hist, 0);
    radixPass10(scratch, source, keys, count, hist, 1);
    radixPass10(source, scratch, keys, count, hist, 2);

    md_temp_end(temp_scope);
}

void md_util_sort_spatial_xyz(uint32_t* source, const float* xyz, size_t stride_in_bytes, size_t count) {
    if (!source || !xyz || count <= 0) return;

    md_temp_scope_t temp_scope = md_temp_begin();

    uint32_t* keys = md_temp_alloc_array(temp_scope, uint32_t, count);
    const float* base = (const float*)xyz;
    const float* x = base + 0;
    const float* y = base + 1;
    const float* z = base + 2;
    computeOrder(keys, x, y, z, count, stride_in_bytes);

    // Important to zero the data here, since we increment when computing the histogram
    uint32_t hist[1024][3] = {0};
    computeHistogram10(hist, keys, count);

    uint32_t* scratch = md_temp_alloc_array(temp_scope, uint32_t, count);
    for (uint32_t i = 0; i < (uint32_t)count; ++i) {
        scratch[i] = i;
    }

    // 3-pass radix sort computes the resulting order into scratch
    radixPass10(source, scratch, keys, count, hist, 0);
    radixPass10(scratch, source, keys, count, hist, 1);
    radixPass10(source, scratch, keys, count, hist, 2);

    md_temp_end(temp_scope);
}

void md_util_sort_radix_u32(uint32_t* data, size_t count) {
    if (!data || count <= 1) return;
    MD_ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_sort_u32(data, count);
        return;
    }
    if (is_sorted_u32(data, count)) return;

    md_temp_scope_t temp_scope = md_temp_begin();
    sort_radix_u32(data, count, temp_scope);
    md_temp_end(temp_scope);
}

void md_util_sort_radix_u64(uint64_t* data, size_t count) {
    if (!data || count <= 1) return;
    MD_ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_sort_u64(data, count);
        return;
    }
    if (is_sorted_u64(data, count)) return;

    md_temp_scope_t temp_scope = md_temp_begin();
    sort_radix_u64(data, count, temp_scope);
    md_temp_end(temp_scope);
}

void md_util_argsort_radix_u32(uint32_t* out_indices, const uint32_t* keys, size_t count) {
    if (count == 0) return;
    MD_ASSERT(out_indices);
    MD_ASSERT(keys);
    MD_ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_argsort_u32(out_indices, keys, count);
        return;
    }
    if (is_sorted_u32(keys, count)) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    argsort_radix_u32(out_indices, keys, count, temp_scope);
    md_temp_end(temp_scope);
}

void md_util_argsort_radix_u64(uint32_t* out_indices, const uint64_t* keys, size_t count) {
    if (count == 0) return;
    MD_ASSERT(out_indices);
    MD_ASSERT(keys);
    MD_ASSERT(count <= UINT32_MAX);

    if (count <= RADIX_SMALL_COUNT) {
        insertion_argsort_u64(out_indices, keys, count);
        return;
    }
    if (is_sorted_u64(keys, count)) {
        for (size_t i = 0; i < count; ++i) out_indices[i] = (uint32_t)i;
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    argsort_radix_u64(out_indices, keys, count, temp_scope);
    md_temp_end(temp_scope);
}
