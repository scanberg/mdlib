#include <md_geom.h>

#include <md_system.h>
#include "md_system_fields.h"

#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_simd.h>

#include <math.h>
#include <string.h>
#include <float.h>

// PACKED XYZ, COMPONENT WISE
//
// Eight packed atoms are 24 floats, and float p of them is component p % 3. An operation that treats
// x, y and z alike but apart - a min image in an orthorhombic cell, a blend, a sum - works on those
// floats as they are loaded, with its per axis constants laid out in the same pattern, and never
// splits them into x, y and z. The 24 floats are XYZ_GROUP_VECS native vectors (md_xv): three of
// eight lanes with AVX, six of four on NEON.
#define XYZ_GROUP_VECS (24 / MD_XV_WIDTH)

// Per axis constants in the lane pattern of a group
static inline void xyz_lane_pattern(md_xv out[XYZ_GROUP_VECS], float x, float y, float z) {
    const float v[3] = { x, y, z };
    float p[24];
    for (int k = 0; k < 24; ++k) p[k] = v[k % 3];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) out[j] = md_xv_loadu_ps(p + j * MD_XV_WIDTH);
}

// Eight per atom values (weights) spread to the lanes of their atoms' three components
static inline void xyz_lane_spread(md_xv out[XYZ_GROUP_VECS], const float* w8) {
#if MD_XV_WIDTH == 8
    const md_256 w = md_mm256_loadu_ps(w8);
    out[0] = md_mm256_permutevar8x32_ps(w, simde_mm256_setr_epi32(0, 0, 0, 1, 1, 1, 2, 2));
    out[1] = md_mm256_permutevar8x32_ps(w, simde_mm256_setr_epi32(2, 3, 3, 3, 4, 4, 4, 5));
    out[2] = md_mm256_permutevar8x32_ps(w, simde_mm256_setr_epi32(5, 5, 6, 6, 6, 7, 7, 7));
#else
    for (int h = 0; h < 2; ++h) {
        const md_128 w = md_mm_loadu_ps(w8 + 4 * h);
        out[3 * h + 0] = md_mm_shuffle_ps(w, w, SIMDE_MM_SHUFFLE(1, 0, 0, 0));   // w0 w0 w0 w1
        out[3 * h + 1] = md_mm_shuffle_ps(w, w, SIMDE_MM_SHUFFLE(2, 2, 1, 1));   // w1 w1 w2 w2
        out[3 * h + 2] = md_mm_shuffle_ps(w, w, SIMDE_MM_SHUFFLE(3, 3, 3, 2));   // w2 w3 w3 w3
    }
#endif
}

// A group's lanes summed per component
static inline void xyz_lane_reduce_add(double out[3], const md_xv v[XYZ_GROUP_VECS]) {
    float t[24];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) md_xv_storeu_ps(t + j * MD_XV_WIDTH, v[j]);
    for (int p = 0; p < 24; ++p) out[p % 3] += t[p];
}

// x brought to the image nearest r, lane by lane; a lane with no period (zero) is left alone. As
// md_mm256_deperiodize2_ps: d from a fused multiply-subtract, see the deperiodize comment in md_simd.h.
static inline md_xv xv_deperiodize_ortho(md_xv x, md_xv r, md_xv period, md_xv inv_period) {
    md_xv d = md_xv_fmsub_ps(x, inv_period, md_xv_mul_ps(r, inv_period));
    d = md_xv_sub_ps(d, md_xv_round_ps(d));
    const md_xv p = md_xv_fmadd_ps(d, period, r);
    return md_xv_select_ps(p, x, md_xv_cmpeq_ps(period, md_xv_setzero_ps()));
}

// Atoms to step over before packed coordinates at p are aligned to the given boundary, so every load
// of a group stays within its cache line; from a 16 byte boundary, which is where an md_array's
// elements begin, half of the 32 byte loads would straddle two. Zero if p is not float aligned.
static inline size_t xyz_align_peel_to(const void* p, size_t count, size_t alignment) {
    const size_t a = (size_t)((uintptr_t)p & (alignment - 1));
    if (a & 3) return 0;
    for (size_t k = 0; k < alignment / 4; ++k) {
        if (((a + 12 * k) & (alignment - 1)) == 0) return MD_MIN(k, count);
    }
    return 0;
}

static inline size_t xyz_align_peel(const void* p, size_t count) {
    return xyz_align_peel_to(p, count, 32);
}

// Blatantly stolen from MDAnalysis project
// MDAnalysis --- https://www.mdanalysis.org
// Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
// https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/include/calc_distances.h

static inline void min_image_triclinic(float dx[3], const float box[3][3]) {
    /*
    * Minimum image convention for triclinic systems, modelled after domain.cpp
    * in LAMMPS.
    * Assumes that there is a maximum separation of 1 box length (enforced in
    * dist functions by moving all particles to inside the box before
    * calculating separations).
    * Assumes box having zero values for box[1], box[2] and box[5]:
    *   /  a_x   0    0   \                 /  0    1    2  \
    *   |  b_x  b_y   0   |       indices:  |  3    4    5  |
    *   \  c_x  c_y  c_z  /                 \  6    7    8  /
    */

    // The search below only looks one lattice vector in each direction, so first bring dx into the
    // zone it covers. Without this a separation of more than one box - any point of an unwrapped
    // trajectory which has diffused out of the cell - silently gets a wrong image.
    // Reduce along c, then b, then a: the box is lower triangular, so each step leaves the
    // components already reduced untouched.
    for (int i = 2; i >= 0; --i) {
        if (box[i][i] != 0.0f) {
            const float n = roundf(dx[i] / box[i][i]);
            if (n != 0.0f) {
                for (int j = i; j >= 0; --j) {
                    dx[j] -= n * box[i][j];
                }
            }
        }
    }

    double dx_min[3] = {0.0, 0.0, 0.0};
    double dsq_min = FLT_MAX;
    double dsq;
    double rx;
    double ry[2];
    double rz[3];
    int ix, iy, iz;
    for (ix = -1; ix < 2; ++ix) {
        rx = dx[0] + box[0][0] * ix;
        for (iy = -1; iy < 2; ++iy) {
            ry[0] = rx + box[1][0] * iy;
            ry[1] = dx[1] + box[1][1] * iy;
            for (iz = -1; iz < 2; ++iz) {
                rz[0] = ry[0] + box[2][0] * iz;
                rz[1] = ry[1] + box[2][1] * iz;
                rz[2] = dx[2] + box[2][2] * iz;
                dsq = rz[0] * rz[0] + rz[1] * rz[1] + rz[2] * rz[2];
                if (dsq < dsq_min) {
                    dsq_min = dsq;
                    dx_min[0] = rz[0];
                    dx_min[1] = rz[1];
                    dx_min[2] = rz[2];
                }
            }
        }
    }
    dx[0] = (float)dx_min[0];
    dx[1] = (float)dx_min[1];
    dx[2] = (float)dx_min[2];
}

static inline void simd_min_image_triclinic(md_256 dx[3], const float box[3][3]) {
    // Into the zone the one lattice vector search covers first - see min_image_triclinic
    for (int i = 2; i >= 0; --i) {
        if (box[i][i] != 0.0f) {
            const md_256 n = md_mm256_round_ps(md_mm256_mul_ps(dx[i], md_mm256_set1_ps(1.0f / box[i][i])));
            for (int j = i; j >= 0; --j) {
                dx[j] = md_mm256_sub_ps(dx[j], md_mm256_mul_ps(n, md_mm256_set1_ps(box[i][j])));
            }
        }
    }

    md_256 dx_min[3] = {md_mm256_set1_ps(0), md_mm256_set1_ps(0), md_mm256_set1_ps(0)};
    md_256 dsq_min   = md_mm256_set1_ps(FLT_MAX);
    md_256 dsq;
    md_256 rx;
    md_256 ry[2];
    md_256 rz[3];
    int ix, iy, iz;
    for (ix = -1; ix < 2; ++ix) {
        md_256 fix = md_mm256_set1_ps((float)ix);
        rx = md_mm256_add_ps(dx[0], md_mm256_mul_ps(md_mm256_set1_ps(box[0][0]), fix));
        for (iy = -1; iy < 2; ++iy) {
            md_256 fiy = md_mm256_set1_ps((float)iy);
            ry[0] = md_mm256_add_ps(rx,    md_mm256_mul_ps(md_mm256_set1_ps(box[1][0]), fiy));
            ry[1] = md_mm256_add_ps(dx[1], md_mm256_mul_ps(md_mm256_set1_ps(box[1][1]), fiy));
            for (iz = -1; iz < 2; ++iz) {
                md_256 fiz = md_mm256_set1_ps((float)iz);
                rz[0] = md_mm256_add_ps(ry[0], md_mm256_mul_ps(md_mm256_set1_ps(box[2][0]), fiz));
                rz[1] = md_mm256_add_ps(ry[1], md_mm256_mul_ps(md_mm256_set1_ps(box[2][1]), fiz));
                rz[2] = md_mm256_add_ps(dx[2], md_mm256_mul_ps(md_mm256_set1_ps(box[2][2]), fiz));
                dsq = md_mm256_add_ps(md_mm256_add_ps(md_mm256_mul_ps(rz[0], rz[0]), md_mm256_mul_ps(rz[1], rz[1])), md_mm256_mul_ps(rz[2], rz[2]));
                md_256 mask = md_mm256_cmplt_ps(dsq, dsq_min);
                dsq_min = md_mm256_select_ps(dsq_min, dsq, mask);
                dx_min[0] = md_mm256_select_ps(dx_min[0], rz[0], mask);
                dx_min[1] = md_mm256_select_ps(dx_min[1], rz[1], mask);
                dx_min[2] = md_mm256_select_ps(dx_min[2], rz[2], mask);
            }
        }
    }
    dx[0] = dx_min[0];
    dx[1] = dx_min[1];
    dx[2] = dx_min[2];
}

static inline void deperiodize_triclinic(float x[3], const float rx[3], const float box[3][3]) {
    float dx[3];

    dx[0] = x[0] - rx[0];
    dx[1] = x[1] - rx[1];
    dx[2] = x[2] - rx[2];

    min_image_triclinic(dx, box);

    x[0] = rx[0] + dx[0];
    x[1] = rx[1] + dx[1];
    x[2] = rx[2] + dx[2];
}

static inline void simd_deperiodize_triclinic(md_256 x[3], const md_256 rx[3], const float box[3][3]) {
    md_256 dx[3];

    dx[0] = md_mm256_sub_ps(x[0], rx[0]);
    dx[1] = md_mm256_sub_ps(x[1], rx[1]);
    dx[2] = md_mm256_sub_ps(x[2], rx[2]);

    simd_min_image_triclinic(dx, box);

    x[0] = md_mm256_add_ps(rx[0], dx[0]);
    x[1] = md_mm256_add_ps(rx[1], dx[1]);
    x[2] = md_mm256_add_ps(rx[2], dx[2]);
}

void md_geom_aabb(float out_min[3], float out_max[3], const float* xyz, size_t stride, const int32_t* idx, size_t count) {
    MD_ASSERT(out_min);
    MD_ASSERT(out_max);
    MD_ASSERT(xyz || count == 0);

    if (stride == 0) stride = 3 * sizeof(float);
    md_float4_t aabb_min = md_float4_set1( FLT_MAX);
    md_float4_t aabb_max = md_float4_set1(-FLT_MAX);
    size_t i = 0;

    if (stride == 3 * sizeof(float)) {
        if (!idx) {
            // A minimum and a maximum need no split: float p of a group of eight atoms is component p % 3
            // Scalar up to 32 byte alignment (below)
            const size_t head = xyz_align_peel(xyz, count);
            for (; i < head; ++i) {
                const md_float4_t c = md_float4_set(xyz[i * 3 + 0], xyz[i * 3 + 1], xyz[i * 3 + 2], 0);
                aabb_min = md_float4_min(aabb_min, c);
                aabb_max = md_float4_max(aabb_max, c);
            }
            md_xv v_min[XYZ_GROUP_VECS], v_max[XYZ_GROUP_VECS];
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                v_min[j] = md_xv_set1_ps(+FLT_MAX);
                v_max[j] = md_xv_set1_ps(-FLT_MAX);
            }
            for (; i + 8 <= count; i += 8) {
                for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                    const md_xv v = md_xv_loadu_ps(xyz + i * 3 + j * MD_XV_WIDTH);
                    v_min[j] = md_xv_min_ps(v_min[j], v);
                    v_max[j] = md_xv_max_ps(v_max[j], v);
                }
            }
            float lo[24], hi[24];
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                md_xv_storeu_ps(lo + j * MD_XV_WIDTH, v_min[j]);
                md_xv_storeu_ps(hi + j * MD_XV_WIDTH, v_max[j]);
            }
            for (int p = 0; p < 24; ++p) {
                aabb_min.elem[p % 3] = MD_MIN(aabb_min.elem[p % 3], lo[p]);
                aabb_max.elem[p % 3] = MD_MAX(aabb_max.elem[p % 3], hi[p]);
            }
        } else {
            md_256 vx_min = md_mm256_set1_ps(+FLT_MAX);
            md_256 vy_min = md_mm256_set1_ps(+FLT_MAX);
            md_256 vz_min = md_mm256_set1_ps(+FLT_MAX);
            md_256 vx_max = md_mm256_set1_ps(-FLT_MAX);
            md_256 vy_max = md_mm256_set1_ps(-FLT_MAX);
            md_256 vz_max = md_mm256_set1_ps(-FLT_MAX);
            for (; i + 8 <= count; i += 8) {
                md_256 x, y, z;
                md_mm256_i32gather_xyz_ps(&x, &y, &z, xyz, md_mm256_loadu_si256(idx + i));
                vx_min = md_mm256_min_ps(vx_min, x);
                vy_min = md_mm256_min_ps(vy_min, y);
                vz_min = md_mm256_min_ps(vz_min, z);
                vx_max = md_mm256_max_ps(vx_max, x);
                vy_max = md_mm256_max_ps(vy_max, y);
                vz_max = md_mm256_max_ps(vz_max, z);
            }
            aabb_min = md_float4_min(aabb_min, md_float4_set(md_mm256_reduce_min_ps(vx_min), md_mm256_reduce_min_ps(vy_min), md_mm256_reduce_min_ps(vz_min), 0));
            aabb_max = md_float4_max(aabb_max, md_float4_set(md_mm256_reduce_max_ps(vx_max), md_mm256_reduce_max_ps(vy_max), md_mm256_reduce_max_ps(vz_max), 0));
        }
    }

    // The remainder, and any other stride
    for (; i < count; ++i) {
        const size_t j = idx ? (size_t)idx[i] : i;
        const float* p = (const float*)((const char*)xyz + j * stride);
        const md_float4_t c = md_float4_set(p[0], p[1], p[2], 0);
        aabb_min = md_float4_min(aabb_min, c);
        aabb_max = md_float4_max(aabb_max, c);
    }

    md_memcpy(out_min, &aabb_min, sizeof(float) * 3);
    md_memcpy(out_max, &aabb_max, sizeof(float) * 3);
}


#define TRIG_ATAN2_R2_THRESHOLD 1.0e-8

// Contiguous atoms without an index: the sums run on the packed floats, the weights spread to the
// lanes of their atoms' three components. Groups of eight atoms from the first, and the last few as
// one more group padded with zeros, so how the sums group follows the atoms' indices and count alone,
// never where in memory they lie: the same positions give the same bits from any buffer.
static void com_packed(float out_com[3], const md_float3_t* in_xyz, const float* in_w, size_t count) {
    const float* f = (const float*)in_xyz;
    double acc[3] = {0};
    double acc_w = 0;

    md_xv s[XYZ_GROUP_VECS];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) s[j] = md_xv_setzero_ps();
    md_xv sw = md_xv_setzero_ps();
    size_t i = 0;
    // The last count % 8 atoms padded to a group: zero positions and, weighted, zero weights
    float tail_f[24] = {0};
    float tail_w[8]  = {0};
    const size_t body = MD_ROUND_DOWN(count, 8);
    if (body < count) {
        md_memcpy(tail_f, f + body * 3, (count - body) * 3 * sizeof(float));
        if (in_w) md_memcpy(tail_w, in_w + body, (count - body) * sizeof(float));
    }
    for (; i < count; i += 8) {
        const float* g = i < body ? f + i * 3 : tail_f;
        if (in_w) {
            const float* gw = i < body ? in_w + i : tail_w;
            md_xv w[XYZ_GROUP_VECS];
            xyz_lane_spread(w, gw);
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                s[j] = md_xv_fmadd_ps(md_xv_loadu_ps(g + j * MD_XV_WIDTH), w[j], s[j]);
            }
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(gw));
#if MD_XV_WIDTH == 4
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(gw + 4));
#endif
        } else {
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                s[j] = md_xv_add_ps(s[j], md_xv_loadu_ps(g + j * MD_XV_WIDTH));
            }
        }
    }
    acc_w = in_w ? md_xv_reduce_add_ps(sw) : (double)count;
    xyz_lane_reduce_add(acc, s);

    out_com[0] = (float)(acc[0] / acc_w);
    out_com[1] = (float)(acc[1] / acc_w);
    out_com[2] = (float)(acc[2] / acc_w);
}

// The periodic mean in an orthorhombic cell, contiguous atoms without an index. The angle of each
// coordinate is its own axis scaled (M is diagonal), so the sines and cosines are taken on the
// packed floats as they lie and summed per lane; the same angles and the same epilogue as _com_pbc.
static void com_pbc_ortho_packed(float out_com[3], const md_float3_t* in_xyz, const float* in_w, size_t count, const float M[3][3], const float I[3][3]) {
    const float* f = (const float*)in_xyz;
    double acc_c[3] = {0};
    double acc_s[3] = {0};
    double acc_w = 0;

    // From the first atom, sixteen or eight at a time, and the last few as one more group of eight padded with zero
    // weights: every atom's sine and cosine come from the same vector sincos, and how the sums group follows the
    // atoms' indices and count alone, never where in memory they lie
    size_t i = 0;

#if defined(__AVX512F__) && defined(__AVX512DQ__)
    {
        // Sixteen atoms, 48 floats: lane k of the j:th vector holds component (16j + k) % 3
        float   pat[48];
        int32_t spread[48];
        for (int p = 0; p < 48; ++p) {
            pat[p]    = M[p % 3][p % 3];
            spread[p] = p / 3;
        }
        __m512  m5[3], c5[3], s5[3];
        __m512i w5_idx[3];
        for (int j = 0; j < 3; ++j) {
            m5[j] = _mm512_loadu_ps(pat + 16 * j);
            w5_idx[j] = _mm512_loadu_si512(spread + 16 * j);
            c5[j] = _mm512_setzero_ps();
            s5[j] = _mm512_setzero_ps();
        }
        __m512 sw5 = _mm512_setzero_ps();
        for (; i + 16 <= count; i += 16) {
            __m512 w16 = _mm512_setzero_ps();
            if (in_w) {
                w16 = _mm512_loadu_ps(in_w + i);
                sw5 = _mm512_add_ps(sw5, w16);
            }
            for (int j = 0; j < 3; ++j) {
                const __m512 t = _mm512_mul_ps(_mm512_loadu_ps(f + i * 3 + j * 16), m5[j]);
                __m512 s, c;
                md_mm512_sincos_ps(t, &s, &c);
                if (in_w) {
                    const __m512 w = _mm512_permutexvar_ps(w5_idx[j], w16);
                    c5[j] = _mm512_fmadd_ps(c, w, c5[j]);
                    s5[j] = _mm512_fmadd_ps(s, w, s5[j]);
                } else {
                    c5[j] = _mm512_add_ps(c5[j], c);
                    s5[j] = _mm512_add_ps(s5[j], s);
                }
            }
        }
        float tc[48], ts[48];
        for (int j = 0; j < 3; ++j) {
            _mm512_storeu_ps(tc + 16 * j, c5[j]);
            _mm512_storeu_ps(ts + 16 * j, s5[j]);
        }
        for (int p = 0; p < 48; ++p) {
            acc_c[p % 3] += tc[p];
            acc_s[p % 3] += ts[p];
        }
        if (in_w) acc_w += _mm512_reduce_add_ps(sw5);
    }
#endif

    md_xv m[XYZ_GROUP_VECS];
    xyz_lane_pattern(m, M[0][0], M[1][1], M[2][2]);
    md_xv vc[XYZ_GROUP_VECS];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) vc[j] = md_xv_setzero_ps();
    md_xv vs[XYZ_GROUP_VECS];
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) vs[j] = md_xv_setzero_ps();
    md_xv sw = md_xv_setzero_ps();
    for (; i + 8 <= count; i += 8) {
        md_xv w[XYZ_GROUP_VECS];
        if (in_w) {
            xyz_lane_spread(w, in_w + i);
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(in_w + i));
#if MD_XV_WIDTH == 4
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(in_w + i + 4));
#endif
        }
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            const md_xv t = md_xv_mul_ps(md_xv_loadu_ps(f + i * 3 + j * MD_XV_WIDTH), m[j]);
            md_xv s, c;
            md_xv_sincos_ps(t, &s, &c);
            if (in_w) {
                vc[j] = md_xv_fmadd_ps(c, w[j], vc[j]);
                vs[j] = md_xv_fmadd_ps(s, w[j], vs[j]);
            } else {
                vc[j] = md_xv_add_ps(vc[j], c);
                vs[j] = md_xv_add_ps(vs[j], s);
            }
        }
    }
    if (i < count) {
        // The rest: a weight of 1 (or its own) for each atom and 0 for the padding, so the padding's cos 0 = 1
        // adds nothing; c * 1 + acc is c + acc, the same as the unweighted lanes
        const size_t rest = count - i;
        float tf[24] = {0};
        float tw[8]  = {0};
        md_memcpy(tf, f + i * 3, rest * 3 * sizeof(float));
        for (size_t q = 0; q < rest; ++q) tw[q] = in_w ? in_w[i + q] : 1.0f;
        md_xv w[XYZ_GROUP_VECS];
        xyz_lane_spread(w, tw);
        if (in_w) {
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(tw));
#if MD_XV_WIDTH == 4
            sw = md_xv_add_ps(sw, md_xv_loadu_ps(tw + 4));
#endif
        }
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            const md_xv t = md_xv_mul_ps(md_xv_loadu_ps(tf + j * MD_XV_WIDTH), m[j]);
            md_xv s, c;
            md_xv_sincos_ps(t, &s, &c);
            vc[j] = md_xv_fmadd_ps(c, w[j], vc[j]);
            vs[j] = md_xv_fmadd_ps(s, w[j], vs[j]);
        }
        i = count;
    }
    acc_w += in_w ? md_xv_reduce_add_ps(sw) : (double)count;
    xyz_lane_reduce_add(acc_c, vc);
    xyz_lane_reduce_add(acc_s, vs);

    const double inv_w = 1.0 / acc_w;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = MD_PI;
        const double x = acc_c[j] * inv_w;
        const double y = acc_s[j] * inv_w;
        if (x * x + y * y > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

static void com(float* out_com, const md_float3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count) {
    if (!in_idx) {
        com_packed(out_com, in_xyz, in_w, count);
        return;
    }
    size_t i = 0;
    double acc_x = 0;
    double acc_y = 0;
    double acc_z = 0;
    double acc_w = 0;
#if defined(__AVX512F__)
    __m512 vx = _mm512_setzero_ps();
    __m512 vy = _mm512_setzero_ps();
    __m512 vz = _mm512_setzero_ps();
    __m512 vw = _mm512_setzero_ps();

    const size_t simd_count = MD_ROUND_DOWN(count, 16);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 x, y, z;
                md_mm512_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);
                __m512 w    = _mm512_i32gather_ps(idx, in_w, 4);

                vx = _mm512_add_ps(vx, _mm512_mul_ps(x, w));
                vy = _mm512_add_ps(vy, _mm512_mul_ps(y, w));
                vz = _mm512_add_ps(vz, _mm512_mul_ps(z, w));
                vw = _mm512_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512i idx = _mm512_loadu_si512(in_idx + i);
                __m512 x, y, z;
                md_mm512_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);

                vx = _mm512_add_ps(vx, x);
                vy = _mm512_add_ps(vy, y);
                vz = _mm512_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 16) {
                __m512 x, y, z;
                md_mm512_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));
                __m512 w = _mm512_loadu_ps(in_w + i);

                vx = _mm512_add_ps(vx, _mm512_mul_ps(x, w));
                vy = _mm512_add_ps(vy, _mm512_mul_ps(y, w));
                vz = _mm512_add_ps(vz, _mm512_mul_ps(z, w));
                vw = _mm512_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 16) {
                __m512 x, y, z;
                md_mm512_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

                vx = _mm512_add_ps(vx, x);
                vy = _mm512_add_ps(vy, y);
                vz = _mm512_add_ps(vz, z);
            }
        }
    }

    acc_x = _mm512_reduce_add_ps(vx);
    acc_y = _mm512_reduce_add_ps(vy);
    acc_z = _mm512_reduce_add_ps(vz);
    acc_w = _mm512_reduce_add_ps(vw);
#elif defined (__AVX__)
    md_256 vx = md_mm256_setzero_ps();
    md_256 vy = md_mm256_setzero_ps();
    md_256 vz = md_mm256_setzero_ps();
    md_256 vw = md_mm256_setzero_ps();

    const size_t simd_count = MD_ROUND_DOWN(count, 8);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 x, y, z;
                md_mm256_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);
                md_256 w    = md_mm256_i32gather_ps(in_w, idx, 4);

                vx = md_mm256_add_ps(vx, md_mm256_mul_ps(x, w));
                vy = md_mm256_add_ps(vy, md_mm256_mul_ps(y, w));
                vz = md_mm256_add_ps(vz, md_mm256_mul_ps(z, w));
                vw = md_mm256_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256i idx = md_mm256_loadu_si256(in_idx + i);
                md_256 x, y, z;
                md_mm256_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);

                vx = md_mm256_add_ps(vx, x);
                vy = md_mm256_add_ps(vy, y);
                vz = md_mm256_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 8) {
                md_256 x, y, z;
                md_mm256_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));
                md_256 w = md_mm256_loadu_ps(in_w + i);

                vx = md_mm256_add_ps(vx, md_mm256_mul_ps(x, w));
                vy = md_mm256_add_ps(vy, md_mm256_mul_ps(y, w));
                vz = md_mm256_add_ps(vz, md_mm256_mul_ps(z, w));
                vw = md_mm256_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 8) {
                md_256 x, y, z;
                md_mm256_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

                vx = md_mm256_add_ps(vx, x);
                vy = md_mm256_add_ps(vy, y);
                vz = md_mm256_add_ps(vz, z);
            }
        }
    }

    acc_x = md_mm256_reduce_add_ps(vx);
    acc_y = md_mm256_reduce_add_ps(vy);
    acc_z = md_mm256_reduce_add_ps(vz);
    acc_w = md_mm256_reduce_add_ps(vw);

#else
    md_128 vx = md_mm_setzero_ps();
    md_128 vy = md_mm_setzero_ps();
    md_128 vz = md_mm_setzero_ps();
    md_128 vw = md_mm_setzero_ps();

    const size_t simd_count = MD_ROUND_DOWN(count, 4);

    if (in_idx) {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 x, y, z;
                md_mm_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);
                md_128 w    = md_mm_i32gather_ps(in_w, idx, 4);

                vx = md_mm_add_ps(vx, md_mm_mul_ps(x, w));
                vy = md_mm_add_ps(vy, md_mm_mul_ps(y, w));
                vz = md_mm_add_ps(vz, md_mm_mul_ps(z, w));
                vw = md_mm_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128i idx = md_mm_loadu_si128(in_idx + i);
                md_128 x, y, z;
                md_mm_i32gather_xyz_ps(&x, &y, &z, (const float*)in_xyz, idx);

                vx = md_mm_add_ps(vx, x);
                vy = md_mm_add_ps(vy, y);
                vz = md_mm_add_ps(vz, z);
            }
        }
    } else {
        if (in_w) {
            for (; i < simd_count; i += 4) {
                md_128 x, y, z;
                md_mm_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));
                md_128 w = md_mm_loadu_ps(in_w + i);

                vx = md_mm_add_ps(vx, md_mm_mul_ps(x, w));
                vy = md_mm_add_ps(vy, md_mm_mul_ps(y, w));
                vz = md_mm_add_ps(vz, md_mm_mul_ps(z, w));
                vw = md_mm_add_ps(vw, w);
            }
        } else {
            for (; i < simd_count; i += 4) {
                md_128 x, y, z;
                md_mm_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

                vx = md_mm_add_ps(vx, x);
                vy = md_mm_add_ps(vy, y);
                vz = md_mm_add_ps(vz, z);
            }
        }
    }

    acc_x = md_mm_reduce_add_ps(vx);
    acc_y = md_mm_reduce_add_ps(vy);
    acc_z = md_mm_reduce_add_ps(vz);
    acc_w = md_mm_reduce_add_ps(vw);
#endif

    if (in_idx) {
        if (in_w) {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                float w = in_w[idx];
                acc_x += in_xyz[idx].x * w;
                acc_y += in_xyz[idx].y * w;
                acc_z += in_xyz[idx].z * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                int64_t idx = in_idx[i];
                acc_x += in_xyz[idx].x;
                acc_y += in_xyz[idx].y;
                acc_z += in_xyz[idx].z;
            }
            acc_w = (double)count;
        }
    } else {
        if (in_w) {
            for (; i < count; ++i) {
                float w = in_w[i];
                acc_x += in_xyz[i].x * w;
                acc_y += in_xyz[i].y * w;
                acc_z += in_xyz[i].z * w;
                acc_w += w;
            }
        } else {
            for (; i < count; ++i) {
                acc_x += in_xyz[i].x;
                acc_y += in_xyz[i].y;
                acc_z += in_xyz[i].z;
            }
            acc_w = (double)count;
        }
    }

    out_com[0] = (float)(acc_x / acc_w);
    out_com[1] = (float)(acc_y / acc_w);
    out_com[2] = (float)(acc_z / acc_w);
}

// Internal versions for COM computation supporting triclinic and ortho PBC
static void _com_pbc_w(float out_com[3], const md_float3_t* in_xyz, const float* in_w, size_t count, const float M[3][3], const float I[3][3]) {
    MD_ASSERT(out_com);
    MD_ASSERT(in_xyz);
    MD_ASSERT(in_w);
    MD_ASSERT(M);
    MD_ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};
    double acc_w = 0;

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    md_512 v_acc_w = _mm512_setzero_ps();
    const size_t simd_count = MD_ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512 v_x, v_y, v_z;
        md_mm512_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        md_512 v_w = _mm512_loadu_ps(in_w + i);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = _mm512_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = _mm512_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = _mm512_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = _mm512_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = _mm512_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = _mm512_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
    acc_w    = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    md_256 v_acc_w = md_mm256_setzero_ps();
    const size_t simd_count = MD_ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256 v_x, v_y, v_z;
        md_mm256_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        md_256 v_w = md_mm256_loadu_ps(in_w + i);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm256_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm256_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm256_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm256_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm256_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm256_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm256_reduce_add_ps(v_acc_w);
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = MD_ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128 v_x, v_y, v_z;
        md_mm_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        md_128 v_w = md_mm_loadu_ps(in_w + i);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm_reduce_add_ps(v_acc_w);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        double x = in_xyz[i].x;
        double y = in_xyz[i].y;
        double z = in_xyz[i].z;
        double w = in_w[i];
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += w * cos(tx);
        acc_s[0] += w * sin(tx);
        acc_c[1] += w * cos(ty);
        acc_s[1] += w * sin(ty);
        acc_c[2] += w * cos(tz);
        acc_s[2] += w * sin(tz);
        acc_w    += w;
    }

    const double inv_w = 1.0 / acc_w;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = MD_PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

// Internal versions for COM computation supporting triclinic and ortho PBC
static void _com_pbc(float out_com[3], const md_float3_t* in_xyz, size_t count, const float M[3][3], const float I[3][3]) {
    MD_ASSERT(out_com);
    MD_ASSERT(in_xyz);
    MD_ASSERT(M);
    MD_ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    const size_t simd_count = MD_ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512 v_x, v_y, v_z;
        md_mm512_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = _mm512_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = _mm512_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = _mm512_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = _mm512_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = _mm512_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    const size_t simd_count = MD_ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256 v_x, v_y, v_z;
        md_mm256_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm256_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm256_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm256_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm256_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm256_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    const size_t simd_count = MD_ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128 v_x, v_y, v_z;
        md_mm_load_xyz_packed_ps(&v_x, &v_y, &v_z, (const float*)(in_xyz + i));
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        double x = in_xyz[i].x;
        double y = in_xyz[i].y;
        double z = in_xyz[i].z;
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += cos(tx);
        acc_s[0] += sin(tx);
        acc_c[1] += cos(ty);
        acc_s[1] += sin(ty);
        acc_c[2] += cos(tz);
        acc_s[2] += sin(tz);
    }

    const double inv_w = 1.0 / (double)count;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = MD_PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}


static void _com_pbc_i(float out_com[3], const md_float3_t* in_xyz, const int32_t* in_idx, size_t count, const float M[3][3], const float I[3][3]) {
    MD_ASSERT(out_com);
    MD_ASSERT(in_xyz);
    MD_ASSERT(M);
    MD_ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    const size_t simd_count = MD_ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512i idx = _mm512_loadu_si512(in_idx + i);
        md_512 v_x, v_y, v_z;
        md_mm512_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = _mm512_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = _mm512_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = _mm512_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = _mm512_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = _mm512_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    const size_t simd_count = MD_ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256i idx = md_mm256_loadu_si256(in_idx + i);
        md_256 v_x, v_y, v_z;
        md_mm256_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm256_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm256_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm256_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm256_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm256_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    const size_t simd_count = MD_ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128i idx = md_mm_loadu_si128(in_idx + i);
        md_128 v_x, v_y, v_z;
        md_mm_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_add_ps(v_acc_s[0], v_sx);
        v_acc_c[0] = md_mm_add_ps(v_acc_c[0], v_cx);
        v_acc_s[1] = md_mm_add_ps(v_acc_s[1], v_sy);
        v_acc_c[1] = md_mm_add_ps(v_acc_c[1], v_cy);
        v_acc_s[2] = md_mm_add_ps(v_acc_s[2], v_sz);
        v_acc_c[2] = md_mm_add_ps(v_acc_c[2], v_cz);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        int32_t idx = in_idx[i];
        double x = in_xyz[idx].x;
        double y = in_xyz[idx].y;
        double z = in_xyz[idx].z;
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += cos(tx);
        acc_s[0] += sin(tx);
        acc_c[1] += cos(ty);
        acc_s[1] += sin(ty);
        acc_c[2] += cos(tz);
        acc_s[2] += sin(tz);
    }

    const double inv_w = 1.0 / (double)count;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = MD_PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

static void _com_pbc_iw(float out_com[3], const md_float3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const float M[3][3], const float I[3][3]) {
    MD_ASSERT(out_com);
    MD_ASSERT(in_xyz);
    MD_ASSERT(in_w);
    MD_ASSERT(M);
    MD_ASSERT(I);

    size_t i = 0;
    double acc_c[3] = {0};
    double acc_s[3] = {0};
    double acc_w = 0;

#if defined(__AVX512F__)
    md_512 v_acc_c[3] = { 0 };
    md_512 v_acc_s[3] = { 0 };
    md_512 v_acc_w = _mm512_setzero_ps();
    const size_t simd_count = MD_ROUND_DOWN(count, 16);
    for (; i < simd_count; i += 16) {
        // Load
        md_512i idx = _mm512_loadu_si512(in_idx + i);
        md_512 v_x, v_y, v_z;
        md_mm512_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        md_512 v_w  = _mm512_i32gather_ps(idx, in_w, 4);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_512 v_tx = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][0]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][0]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][0]))));
        md_512 v_ty = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][1]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][1]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][1]))));
        md_512 v_tz = _mm512_fmadd_ps(v_x, _mm512_set1_ps(M[0][2]), _mm512_fmadd_ps(v_y, _mm512_set1_ps(M[1][2]), _mm512_mul_ps(v_z, _mm512_set1_ps(M[2][2]))));
        // Compute sin cos
        md_512 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm512_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm512_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm512_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = _mm512_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = _mm512_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = _mm512_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = _mm512_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = _mm512_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = _mm512_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = _mm512_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = _mm512_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = _mm512_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = _mm512_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = _mm512_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = _mm512_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = _mm512_reduce_add_ps(v_acc_c[2]);
    acc_w    = _mm512_reduce_add_ps(v_acc_w);
#elif defined(__AVX2__)
    md_256 v_acc_c[3] = { 0 };
    md_256 v_acc_s[3] = { 0 };
    md_256 v_acc_w = md_mm256_setzero_ps();
    const size_t simd_count = MD_ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        // Load
        md_256i idx = md_mm256_loadu_si256(in_idx + i);
        md_256 v_x, v_y, v_z;
        md_mm256_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        md_256 v_w  = md_mm256_i32gather_ps(in_w, idx, 4);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_256 v_tx = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][0]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][0]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][0]))));
        md_256 v_ty = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][1]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][1]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][1]))));
        md_256 v_tz = md_mm256_fmadd_ps(v_x, md_mm256_set1_ps(M[0][2]), md_mm256_fmadd_ps(v_y, md_mm256_set1_ps(M[1][2]), md_mm256_mul_ps(v_z, md_mm256_set1_ps(M[2][2]))));
        // Compute sin cos
        md_256 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm256_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm256_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm256_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm256_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm256_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm256_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm256_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm256_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm256_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm256_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm256_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm256_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm256_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm256_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm256_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm256_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm256_reduce_add_ps(v_acc_w);
#else
    md_128 v_acc_c[3] = { 0 };
    md_128 v_acc_s[3] = { 0 };
    md_128 v_acc_w = md_mm_setzero_ps();
    const size_t simd_count = MD_ROUND_DOWN(count, 4);
    for (; i < simd_count; i += 4) {
        // Load
        md_128i idx = md_mm_loadu_si128(in_idx + i);
        md_128 v_x, v_y, v_z;
        md_mm_i32gather_xyz_ps(&v_x, &v_y, &v_z, (const float*)in_xyz, idx);
        md_128 v_w  = md_mm_i32gather_ps(in_w, idx, 4);
        // Compute thetas: theta = 2*pi * Ai * p, the fractional coordinate carried into angle space.
        // Every component of theta draws on all three Cartesian components, because a triclinic Ai is
        // not diagonal. Resolving an axis against itself alone is only correct for an ortho cell.
        md_128 v_tx = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][0]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][0]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][0]))));
        md_128 v_ty = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][1]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][1]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][1]))));
        md_128 v_tz = md_mm_fmadd_ps(v_x, md_mm_set1_ps(M[0][2]), md_mm_fmadd_ps(v_y, md_mm_set1_ps(M[1][2]), md_mm_mul_ps(v_z, md_mm_set1_ps(M[2][2]))));
        // Compute sin cos
        md_128 v_cx, v_cy, v_cz, v_sx, v_sy, v_sz;
        md_mm_sincos_ps(v_tx, &v_sx, &v_cx);
        md_mm_sincos_ps(v_ty, &v_sy, &v_cy);
        md_mm_sincos_ps(v_tz, &v_sz, &v_cz);
        // Accumulate
        v_acc_s[0] = md_mm_fmadd_ps(v_sx, v_w, v_acc_s[0]);
        v_acc_c[0] = md_mm_fmadd_ps(v_cx, v_w, v_acc_c[0]);
        v_acc_s[1] = md_mm_fmadd_ps(v_sy, v_w, v_acc_s[1]);
        v_acc_c[1] = md_mm_fmadd_ps(v_cy, v_w, v_acc_c[1]);
        v_acc_s[2] = md_mm_fmadd_ps(v_sz, v_w, v_acc_s[2]);
        v_acc_c[2] = md_mm_fmadd_ps(v_cz, v_w, v_acc_c[2]);
        v_acc_w    = md_mm_add_ps(v_acc_w, v_w);
    }
    // Reduce
    acc_s[0] = md_mm_reduce_add_ps(v_acc_s[0]);
    acc_c[0] = md_mm_reduce_add_ps(v_acc_c[0]);
    acc_s[1] = md_mm_reduce_add_ps(v_acc_s[1]);
    acc_c[1] = md_mm_reduce_add_ps(v_acc_c[1]);
    acc_s[2] = md_mm_reduce_add_ps(v_acc_s[2]);
    acc_c[2] = md_mm_reduce_add_ps(v_acc_c[2]);
    acc_w    = md_mm_reduce_add_ps(v_acc_w);
#endif
    // Scalar remainder
    for (; i < count; ++i) {
        int32_t idx = in_idx[i];
        double x = in_xyz[idx].x;
        double y = in_xyz[idx].y;
        double z = in_xyz[idx].z;
        double w = in_w[idx];
        double tx = x * M[0][0] + y * M[1][0] + z * M[2][0];
        double ty = x * M[0][1] + y * M[1][1] + z * M[2][1];
        double tz = x * M[0][2] + y * M[1][2] + z * M[2][2];
        acc_c[0] += w * cos(tx);
        acc_s[0] += w * sin(tx);
        acc_c[1] += w * cos(ty);
        acc_s[1] += w * sin(ty);
        acc_c[2] += w * cos(tz);
        acc_s[2] += w * sin(tz);
        acc_w    += w;
    }

    const double inv_w = 1.0 / acc_w;
    double theta[3];
    for (int j = 0; j < 3; ++j) {
        double t = MD_PI;
        double x = acc_c[j] * inv_w;
        double y = acc_s[j] * inv_w;
        double r2 = x * x + y * y;
        if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
            t += atan2(-y, -x);
        }
        theta[j] = t;
    }

    // The three mean angles ARE a fractional coordinate (times two pi). Carrying it back out to
    // Cartesian is a full matrix vector product - component j draws on all three angles, not only
    // on its own. I is column major, so I[col][row].
    for (int j = 0; j < 3; ++j) {
        out_com[j] = (float)(theta[0] * I[0][j] + theta[1] * I[1][j] + theta[2] * I[2][j]);
    }
}

// This is uses the trigonometric com algorithm presented in: INSERT REF HERE
static void com_pbc(float out_com[3], const md_float3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {

    md_float3x3_t A = { 0 };
    md_unitcell_A_extract_float(A.elem, cell);

    md_float3x3_t Ai = { 0 };
    md_unitcell_I_extract_float(Ai.elem, cell);

    // Should transform cartesian coordinates to fractional coordinates and scale by 2 pi
    md_float3x3_t M = md_float3x3_mul(md_float3x3_scale(MD_TWO_PI, MD_TWO_PI, MD_TWO_PI), Ai);
    // Inverse transform
    md_float3x3_t I = md_float3x3_mul(A, md_float3x3_scale(1.0/MD_TWO_PI, 1.0/MD_TWO_PI, 1.0/MD_TWO_PI));

    if (!in_idx && md_unitcell_is_orthorhombic(cell)) {
        com_pbc_ortho_packed(out_com, in_xyz, in_w, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        return;
    }

    // Here we select the correct internal version based on the available inputs
    if (in_w) {
        if (in_idx) {
            _com_pbc_iw(out_com, in_xyz, in_w, in_idx, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        } else {
            _com_pbc_w(out_com, in_xyz, in_w, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        }
    } else {
        if (in_idx) {
            _com_pbc_i(out_com, in_xyz, in_idx, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        } else {
            _com_pbc(out_com, in_xyz, count, (const float (*)[3])M.elem, (const float (*)[3])I.elem);
        }
    }
}

// Point i of a strided set as (x, y, z, w), w one without weights
static inline md_float4_t com_load_xyzw(const float* xyz, size_t xyz_stride, const float* w, size_t w_stride, size_t i) {
    const float* p = (const float*)((const char*)xyz + i * xyz_stride);
    return md_float4_set(p[0], p[1], p[2], w ? *(const float*)((const char*)w + i * w_stride) : 1.0f);
}

// Any strides: a point at a time
static void com_strided(float out_com[3], const float* xyz, size_t xyz_stride, const float* w, size_t w_stride, const int32_t* in_idx, size_t count) {
    md_float4_t acc = {0,0,0,0};
    for (size_t i = 0; i < count; ++i) {
        const size_t idx = in_idx ? (size_t)in_idx[i] : i;
        md_float4_t xyzw = com_load_xyzw(xyz, xyz_stride, w, w_stride, idx);
        md_float4_t www1 = md_float4_blend_mask(md_float4_splat_w(xyzw), md_float4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
        acc = md_float4_add(acc, md_float4_mul(xyzw, www1));
    }

    acc = md_float4_div(acc, md_float4_splat_w(acc));
    md_memcpy(out_com, &acc, sizeof(float) * 3);
}

static void com_pbc_strided(float* out_com, const float* xyz, size_t xyz_stride, const float* w, size_t w_stride, const int32_t* in_idx, size_t count, const md_unitcell_t* cell) {
    MD_ASSERT(out_com);
    MD_ASSERT(xyz);
    MD_ASSERT(cell);

    md_float4_t acc_c = {0};
    md_float4_t acc_s = {0};
    md_float4_t acc_xyzw = {0};

    uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        md_float3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
        const md_float4_t scl = md_float4_div(md_float4_set1(MD_TWO_PI), md_float4_from_float3(ext, MD_TWO_PI));

        if (in_idx) {
            for (size_t i = 0; i < count; ++i) {
                int32_t idx  = in_idx[i];
                md_float4_t xyzw  = com_load_xyzw(xyz, xyz_stride, w, w_stride, (size_t)idx);
                md_float4_t www1  = md_float4_blend_mask(md_float4_splat_w(xyzw), md_float4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                md_float4_t theta = md_float4_mul(xyzw, scl);
                md_float4_t s,c;
                md_float4_sincos(theta, &s, &c);
                acc_s = md_float4_add(acc_s, md_float4_mul(s, www1));
                acc_c = md_float4_add(acc_c, md_float4_mul(c, www1));
                acc_xyzw = md_float4_add(acc_xyzw, md_float4_mul(xyzw, www1));
            }
        } else {
            for (size_t i = 0; i < count; ++i) {
                md_float4_t xyzw  = com_load_xyzw(xyz, xyz_stride, w, w_stride, i);
                md_float4_t www1  = md_float4_blend_mask(md_float4_splat_w(xyzw), md_float4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                md_float4_t theta = md_float4_mul(xyzw, scl);
                md_float4_t s,c;
                md_float4_sincos(theta, &s, &c);
                acc_s = md_float4_add(acc_s, md_float4_mul(s, www1));
                acc_c = md_float4_add(acc_c, md_float4_mul(c, www1));
                acc_xyzw = md_float4_add(acc_xyzw, md_float4_mul(xyzw, www1));
            }
        }

        const float w = acc_xyzw.w;
        for (int i = 0; i < 3; ++i) {
            const double y = acc_s.elem[i] / w;
            const double x = acc_c.elem[i] / w;
            const double r2 = x*x + y*y;
            double theta_prim = MD_PI;
            if (r2 > 1.0e-15) {
                theta_prim += atan2(-y, -x);
            }
            out_com[i] = (float)((theta_prim / MD_TWO_PI) * ext.elem[i]);
        }
    } else {
        MD_ASSERT(flags & MD_UNITCELL_TRICLINIC);
        md_float3x3_t Ai = { 0 };
        md_unitcell_I_extract_float(Ai.elem, cell);

        md_float3x3_t A = { 0 };
        md_unitcell_A_extract_float(A.elem, cell);

        // Cartesian -> fractional, scaled into angle space, and the inverse which takes the mean
        // angle back out to a Cartesian position. A triclinic basis mixes the components, so both
        // directions have to be full matrix products: resolving an axis on its own is only valid
        // when the basis is diagonal, which is why the orthorhombic branch above can get away with
        // per component scalars and this one cannot.
        const md_float4x3_t M = md_float4x3_from_float3x3(md_float3x3_mul(md_float3x3_scale(MD_TWO_PI, MD_TWO_PI, MD_TWO_PI), Ai));
        const md_float3x3_t   I = md_float3x3_mul(A, md_float3x3_scale(1.0f / MD_TWO_PI, 1.0f / MD_TWO_PI, 1.0f / MD_TWO_PI));

        if (in_idx) {
            for (size_t i = 0; i < count; ++i) {
                int32_t idx  = in_idx[i];
                md_float4_t xyzw  = com_load_xyzw(xyz, xyz_stride, w, w_stride, (size_t)idx);
                md_float4_t www1  = md_float4_blend_mask(md_float4_splat_w(xyzw), md_float4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                md_float4_t theta = md_float4x3_mul_float4(M, xyzw);
                md_float4_t s,c;
                md_float4_sincos(theta, &s, &c);
                acc_s = md_float4_add(acc_s, md_float4_mul(s, www1));
                acc_c = md_float4_add(acc_c, md_float4_mul(c, www1));
                acc_xyzw = md_float4_add(acc_xyzw, md_float4_mul(xyzw, www1));
            }
        } else {
            for (size_t i = 0; i < count; ++i) {
                md_float4_t xyzw  = com_load_xyzw(xyz, xyz_stride, w, w_stride, i);
                md_float4_t www1  = md_float4_blend_mask(md_float4_splat_w(xyzw), md_float4_set1(1.0f), MD_SIMD_BLEND_MASK(0,0,0,1));
                md_float4_t theta = md_float4x3_mul_float4(M, xyzw);
                md_float4_t s,c;
                md_float4_sincos(theta, &s, &c);
                acc_s = md_float4_add(acc_s, md_float4_mul(s, www1));
                acc_c = md_float4_add(acc_c, md_float4_mul(c, www1));
                acc_xyzw = md_float4_add(acc_xyzw, md_float4_mul(xyzw, www1));
            }
        }

        md_float3_t theta_vec = { 0 };
        for (int i = 0; i < 3; ++i) {
            const double y = acc_s.elem[i] / acc_xyzw.w;
            const double x = acc_c.elem[i] / acc_xyzw.w;
            const double r2 = x*x + y*y;
            double theta_prim = MD_PI;
            if (r2 > TRIG_ATAN2_R2_THRESHOLD) {
                theta_prim += atan2(-y, -x);
            }
            theta_vec.elem[i] = (float)theta_prim;
        }

        // The three mean angles ARE a fractional coordinate (times two pi). Carry it back to
        // Cartesian as one vector - component i of the result draws on all three angles.
        const md_float3_t com = md_float3x3_mul_float3(I, theta_vec);
        out_com[0] = com.x;
        out_com[1] = com.y;
        out_com[2] = com.z;
    }
}

md_float3_t md_geom_com(const float* xyz, size_t xyz_stride, const float* w, size_t w_stride, const int32_t* idx, size_t count, const md_unitcell_t* cell) {
    MD_ASSERT(xyz || count == 0);

    md_float3_t com_xyz = {0};
    if (count == 0) {
        return com_xyz;
    }
    if (xyz_stride == 0) xyz_stride = 3 * sizeof(float);
    if (w_stride == 0)   w_stride   = sizeof(float);

    const bool periodic = cell && (cell->flags & (MD_UNITCELL_ORTHO | MD_UNITCELL_TRICLINIC));
    if (cell && cell->flags != MD_UNITCELL_NONE && !periodic) {
        MD_LOG_ERROR("Invalid unit cell flags: %d", cell->flags);
        return com_xyz;
    }

    // Packed coordinates and weights: the vectorized kernels
    if (xyz_stride == 3 * sizeof(float) && (!w || w_stride == sizeof(float))) {
        if (periodic) {
            com_pbc(com_xyz.elem, (const md_float3_t*)xyz, w, idx, count, cell);
        } else {
            com(com_xyz.elem, (const md_float3_t*)xyz, w, idx, count);
        }
    } else {
        if (periodic) {
            com_pbc_strided(com_xyz.elem, xyz, xyz_stride, w, w_stride, idx, count, cell);
        } else {
            com_strided(com_xyz.elem, xyz, xyz_stride, w, w_stride, idx, count);
        }
    }
    return com_xyz;
}

#if MD_COMPILER_MSVC
#   pragma warning( push )
#   pragma warning( disable : 4244 )
#endif

void md_geom_distance_array(float* out_dist, const md_float3_t* coord_a, size_t len_a, const md_float3_t* coord_b, size_t len_b, const md_unitcell_t* cell) {
    if (cell->flags == 0) {
        for (size_t i = 0; i < len_a; ++i) {
            for (size_t j = 0; j < len_b; ++j) {
                out_dist[i * len_b + j] = md_float3_distance(coord_a[i], coord_b[j]);
            }
        }
    }
    else if (cell->flags & MD_UNITCELL_ORTHO) {
		md_float4_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
        for (size_t i = 0; i < len_a; ++i) {
            md_float4_t a = md_float4_from_float3(coord_a[i], 0);
            for (size_t j = 0; j < len_b; ++j) {
                md_float4_t b = md_float4_from_float3(coord_b[j], 0);
                out_dist[i * len_b + j] = md_float4_periodic_distance(a, b, ext);
            }
        }
    } else if (cell->flags & MD_UNITCELL_TRICLINIC) {
        float A[3][3] = {0};
        md_unitcell_A_extract_float(A, cell);
        // We make the assumption that we are not beyond 1 cell unit in distance
        for (size_t i = 0; i < len_a; ++i) {
            for (size_t j = 0; j < len_b; ++j) {
                md_float3_t dx = md_float3_sub(coord_a[i], coord_b[j]);
                min_image_triclinic(dx.elem, MD_AS_CONST_MAT3(A));
                out_dist[i * len_b + j] = md_float3_length(dx);
            }
        }
    }
}

// ### MINIMUM AND MAXIMUM DISTANCE BETWEEN SETS ###
//
// For each group of points of a, the nearest (farthest) point of b under the minimum image convention.
//
// There is no spatial structure. The points of a group are taken eight at a time, a cluster in one SIMD register, and
// b is streamed once per cluster, each point of b tested against the sphere which bounds the cluster: a point of b can
// only improve on the best pair so far if it lies within best + R of the centre (beyond best - R for the farthest
// pair). All but a handful are rejected by that one test, and those are tested against every point of the cluster.
// The best pair is shared by the clusters of a group and starts from the pair of the previous group, which is usually
// close, so the test is tight from the first point. A grid over b would make the cost of a cluster sublinear in b,
// but building one costs more than the whole stream until a has a few hundred clusters, and the common case is a
// small a (a ligand, a residue, a QM region) against a large b.
//
// Minimum image. A displacement is first reduced into the brick |d_x| <= x/2, |d_y| <= y/2, |d_z| <= z/2: along c,
// then b, then a, which the lower triangular basis allows without disturbing the components already reduced. In an
// orthorhombic cell that is the minimum image. In a triclinic one it is the minimum image whenever that is shorter
// than R_SAFE = min(x, y, z) / 2: the brick is a fundamental domain of the lattice, the reduction returns the image
// inside it, and a vector shorter than R_SAFE lies inside it. So a reduced length below R_SAFE is exact, and one at or
// above it only says that the true length lies between R_SAFE and the reduced one. Where that does not decide the
// comparison at hand, the minimum image is searched for (dist_search_d2). That takes a group about as far from b as
// half the cell, or the farthest pair.
//
// The search is bounded rather than assumed: md_geom_min_image looks one lattice vector either way of the reduced
// image, which is the minimum image for a cell of a sensible shape but not for every cell. The bound is two at most for
// the usual cells (a reduced triclinic box, a rhombic dodecahedron, a truncated octahedron), and only a degenerate cell
// would reach DIST_MAX_SEARCH.

enum {
    DIST_PBC_NONE      = 0,
    DIST_PBC_ORTHO     = 1,
    DIST_PBC_TRICLINIC = 2,
};

// The search goes no further than this many lattice vectors along an axis. A cell which needs more is degenerate.
#define DIST_MAX_SEARCH 8

typedef struct dist_cell_t {
    // The basis, lower triangular: a = (ax, 0, 0), b = (bx, by, 0), c = (cx, cy, cz). Zero along an axis without a
    // period, as is the inverse of its diagonal (ix, iy, iz), which leaves such an axis alone in the reduction.
    md_256 ax, bx, by, cx, cy, cz;
    md_256 ix, iy, iz;
    float  basis[3][3];     // [vector][residue]
    float  dual[3];         // Lengths of the rows of the inverse basis, zero along an axis without a period
    float  rsafe2;          // A reduced length squared below this is the minimum image
    float  eps;             // Slack for the rounding of the coordinates in the sphere test
    int    mode;
} dist_cell_t;

static void dist_cell_init(dist_cell_t* c, const md_unitcell_t* cell, float coord_scale) {
    md_memset(c, 0, sizeof(*c));
    const uint32_t flags = cell ? (uint32_t)cell->flags : 0;
    const bool pbc[3] = {
        (flags & MD_UNITCELL_PBC_X) != 0,
        (flags & MD_UNITCELL_PBC_Y) != 0,
        (flags & MD_UNITCELL_PBC_Z) != 0,
    };

    double A[3][3] = {{0}};
    if (cell) md_unitcell_A_extract_double(A, cell);
    for (int i = 0; i < 3; ++i) {
        if (!pbc[i] || A[i][i] == 0.0) {
            A[i][0] = A[i][1] = A[i][2] = 0.0;
        }
    }

    const bool any = A[0][0] != 0.0 || A[1][1] != 0.0 || A[2][2] != 0.0;
    c->mode = !any ? DIST_PBC_NONE : (flags & MD_UNITCELL_TRICLINIC) ? DIST_PBC_TRICLINIC : DIST_PBC_ORTHO;
    if (c->mode == DIST_PBC_ORTHO) {
        A[1][0] = A[2][0] = A[2][1] = 0.0;
    }

    double extent = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            c->basis[i][j] = (float)A[i][j];
            extent += fabs(A[i][j]);
        }
    }

    c->ax = md_mm256_set1_ps((float)A[0][0]);
    c->bx = md_mm256_set1_ps((float)A[1][0]);
    c->by = md_mm256_set1_ps((float)A[1][1]);
    c->cx = md_mm256_set1_ps((float)A[2][0]);
    c->cy = md_mm256_set1_ps((float)A[2][1]);
    c->cz = md_mm256_set1_ps((float)A[2][2]);
    c->ix = md_mm256_set1_ps(A[0][0] != 0.0 ? (float)(1.0 / A[0][0]) : 0.0f);
    c->iy = md_mm256_set1_ps(A[1][1] != 0.0 ? (float)(1.0 / A[1][1]) : 0.0f);
    c->iz = md_mm256_set1_ps(A[2][2] != 0.0 ? (float)(1.0 / A[2][2]) : 0.0f);

    // Several ulps of the largest magnitude a difference is taken between
    c->eps = (float)((coord_scale + extent) * 4.0e-6 + 1.0e-5);

    // The rows of the inverse of the basis, completed by a unit vector along each axis without a period. A lattice
    // vector A n has n_i = r_i . A n, so |n_i| <= |r_i| |A n|.
    double M[3][3];
    md_memcpy(M, A, sizeof(M));
    for (int i = 0; i < 3; ++i) {
        if (M[i][i] == 0.0) {
            M[i][0] = M[i][1] = M[i][2] = 0.0;
            M[i][i] = 1.0;
        }
    }
    const double x = M[0][0], xy = M[1][0], y = M[1][1], xz = M[2][0], yz = M[2][1], z = M[2][2];
    const double r0[3] = { 1.0 / x, -xy / (x * y), (xy * yz - xz * y) / (x * y * z) };
    const double r1[3] = { 0.0, 1.0 / y, -yz / (y * z) };
    const double r2[3] = { 0.0, 0.0, 1.0 / z };
    c->dual[0] = A[0][0] != 0.0 ? (float)sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]) : 0.0f;
    c->dual[1] = A[1][1] != 0.0 ? (float)sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]) : 0.0f;
    c->dual[2] = A[2][2] != 0.0 ? (float)sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]) : 0.0f;

    if (c->mode == DIST_PBC_TRICLINIC) {
        double rsafe = DBL_MAX;
        for (int i = 0; i < 3; ++i) {
            if (A[i][i] != 0.0) rsafe = MD_MIN(rsafe, 0.5 * fabs(A[i][i]));
        }
        rsafe = MD_MAX(0.0, rsafe - c->eps);
        c->rsafe2 = (float)(rsafe * rsafe);
    } else {
        c->rsafe2 = FLT_MAX;
    }
}

// Into the brick, along c, then b, then a
static MD_FORCE_INLINE void dist_reduce(md_256* dx, md_256* dy, md_256* dz, const dist_cell_t* c, const int mode) {
    if (mode == DIST_PBC_ORTHO) {
        *dx = md_mm256_fnmadd_ps(md_mm256_round_ps(md_mm256_mul_ps(*dx, c->ix)), c->ax, *dx);
        *dy = md_mm256_fnmadd_ps(md_mm256_round_ps(md_mm256_mul_ps(*dy, c->iy)), c->by, *dy);
        *dz = md_mm256_fnmadd_ps(md_mm256_round_ps(md_mm256_mul_ps(*dz, c->iz)), c->cz, *dz);
    } else if (mode == DIST_PBC_TRICLINIC) {
        const md_256 nz = md_mm256_round_ps(md_mm256_mul_ps(*dz, c->iz));
        *dx = md_mm256_fnmadd_ps(nz, c->cx, *dx);
        *dy = md_mm256_fnmadd_ps(nz, c->cy, *dy);
        *dz = md_mm256_fnmadd_ps(nz, c->cz, *dz);
        const md_256 ny = md_mm256_round_ps(md_mm256_mul_ps(*dy, c->iy));
        *dx = md_mm256_fnmadd_ps(ny, c->bx, *dx);
        *dy = md_mm256_fnmadd_ps(ny, c->by, *dy);
        const md_256 nx = md_mm256_round_ps(md_mm256_mul_ps(*dx, c->ix));
        *dx = md_mm256_fnmadd_ps(nx, c->ax, *dx);
    }
}

static MD_FORCE_INLINE md_256 dist_length2(md_256 dx, md_256 dy, md_256 dz) {
    return md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
}

// The minimum image length squared of reduced displacements, by search. A lattice vector A n which shortens a reduced
// d has |A n| <= |d| + |d + A n| <= 2 |d|, so |n_i| <= 2 |d| |r_i|. In a cell of a sensible shape that is one vector
// either way along each axis, or two along a short one.
static md_256 dist_search_d2(md_256 dx, md_256 dy, md_256 dz, const dist_cell_t* c) {
    md_256 best = dist_length2(dx, dy, dz);
    const float len = sqrtf(md_mm256_reduce_max_ps(best));
    int k[3];
    for (int i = 0; i < 3; ++i) {
        const float n = floorf(2.0f * len * c->dual[i] * 1.0001f);
        k[i] = (int)MD_CLAMP(n, 0.0f, (float)DIST_MAX_SEARCH);
    }
    for (int nz = -k[2]; nz <= k[2]; ++nz) {
        for (int ny = -k[1]; ny <= k[1]; ++ny) {
            for (int nx = -k[0]; nx <= k[0]; ++nx) {
                if (nx == 0 && ny == 0 && nz == 0) continue;
                const float sx = nx * c->basis[0][0] + ny * c->basis[1][0] + nz * c->basis[2][0];
                const float sy = ny * c->basis[1][1] + nz * c->basis[2][1];
                const float sz = nz * c->basis[2][2];
                const md_256 vx = md_mm256_add_ps(dx, md_mm256_set1_ps(sx));
                const md_256 vy = md_mm256_add_ps(dy, md_mm256_set1_ps(sy));
                const md_256 vz = md_mm256_add_ps(dz, md_mm256_set1_ps(sz));
                best = md_mm256_min_ps(best, dist_length2(vx, vy, vz));
            }
        }
    }
    return best;
}

// Up to eight points of a group, as one image of them: each moved next to the first. Only the size of the bounding
// sphere depends on which image, the distances do not. Lanes past the end repeat the first point.
typedef struct dist_cluster_t {
    md_256  x, y, z;
    float   cen[3];
    float   rad;
    int     count;      // The lanes which are points of the group
    int64_t idx[8];
} dist_cluster_t;

static MD_FORCE_INLINE void dist_cluster_init(dist_cluster_t* cl, const md_float3_t* a, size_t beg, size_t n, const dist_cell_t* c, const int mode) {
    const md_float3_t p0 = a[beg];
    float x[8], y[8], z[8];
    for (size_t i = 0; i < 8; ++i) {
        const size_t k = beg + (i < n ? i : 0);
        x[i] = a[k].x - p0.x;
        y[i] = a[k].y - p0.y;
        z[i] = a[k].z - p0.z;
        cl->idx[i] = (int64_t)k;
    }
    cl->count = (int)n;
    md_256 vx = md_mm256_loadu_ps(x);
    md_256 vy = md_mm256_loadu_ps(y);
    md_256 vz = md_mm256_loadu_ps(z);
    dist_reduce(&vx, &vy, &vz, c, mode);
    cl->x = md_mm256_add_ps(vx, md_mm256_set1_ps(p0.x));
    cl->y = md_mm256_add_ps(vy, md_mm256_set1_ps(p0.y));
    cl->z = md_mm256_add_ps(vz, md_mm256_set1_ps(p0.z));
    md_mm256_storeu_ps(x, cl->x);
    md_mm256_storeu_ps(y, cl->y);
    md_mm256_storeu_ps(z, cl->z);

    float lo[3] = { x[0], y[0], z[0] };
    float hi[3] = { x[0], y[0], z[0] };
    for (int i = 1; i < 8; ++i) {
        lo[0] = MD_MIN(lo[0], x[i]); hi[0] = MD_MAX(hi[0], x[i]);
        lo[1] = MD_MIN(lo[1], y[i]); hi[1] = MD_MAX(hi[1], y[i]);
        lo[2] = MD_MIN(lo[2], z[i]); hi[2] = MD_MAX(hi[2], z[i]);
    }
    cl->cen[0] = 0.5f * (lo[0] + hi[0]);
    cl->cen[1] = 0.5f * (lo[1] + hi[1]);
    cl->cen[2] = 0.5f * (lo[2] + hi[2]);
    float r2 = 0.0f;
    for (int i = 0; i < 8; ++i) {
        const float dx = x[i] - cl->cen[0], dy = y[i] - cl->cen[1], dz = z[i] - cl->cen[2];
        r2 = MD_MAX(r2, dx * dx + dy * dy + dz * dz);
    }
    cl->rad = sqrtf(r2);
}

// The points of b, one array per axis, padded to a multiple of eight
typedef struct dist_set_t {
    const float* x;
    const float* y;
    const float* z;
    size_t n;
} dist_set_t;

// The pairs near the extreme: those whose distance is within tol of it, as the kernel computes distances. Kept while
// they are near the best so far, which only ever gets more extreme, so every pair near the final one is among them.
typedef struct dist_near_t {
    float    tol;
    float    bound2;     // dist_near_bound2 of the best so far
    int64_t* a;
    int64_t* b;
    float*   d2;
    size_t   cap;
    size_t   count;
    bool     overflow;   // A near pair did not fit: the ones kept are not all of them
} dist_near_t;

// What a squared distance must be within to be near best2: at most (best + tol)^2, or for the farthest at least
// (best - tol)^2, where -1 takes every pair. Never short of best2 itself, whatever the rounding.
static MD_FORCE_INLINE float dist_near_bound2(float best2, float tol, const bool largest) {
    if (largest) {
        if (best2 < 0.0f) return -1.0f;
        const float t = sqrtf(best2) - tol;
        return t > 0.0f ? MD_MIN(t * t, best2) : -1.0f;
    }
    if (best2 >= FLT_MAX) return FLT_MAX;
    const float t = sqrtf(best2) + tol;
    return MD_MAX(t * t, best2);
}

static void dist_near_prune(dist_near_t* nr, const bool largest) {
    const float e = nr->bound2;
    size_t w = 0;
    for (size_t r = 0; r < nr->count; ++r) {
        if (largest ? nr->d2[r] >= e : nr->d2[r] <= e) {
            nr->a[w] = nr->a[r];
            nr->b[w] = nr->b[r];
            nr->d2[w] = nr->d2[r];
            ++w;
        }
    }
    nr->count = w;
}

// The lanes of mask, near pairs with point j of b; lane i is point first + i of a
static void dist_near_record(dist_near_t* nr, int64_t first, const float v[8], size_t j, int mask, const bool largest) {
    while (mask) {
        const int lane = md_ctz32((uint32_t)mask);
        mask &= mask - 1;
        if (nr->count == nr->cap) dist_near_prune(nr, largest);
        if (nr->count == nr->cap) {
            nr->overflow = true;
            continue;
        }
        nr->a[nr->count] = first + lane;
        nr->b[nr->count] = (int64_t)j;
        nr->d2[nr->count] = v[lane];
        nr->count += 1;
    }
}

// The squared distances from the points of the cluster to point j of b, exact wherever that can change the
// comparison with best2. If one improves on best2 it becomes the best pair and true is returned. Given nr, they are
// exact wherever they can be near, and the near ones are recorded unless record is false.
static MD_FORCE_INLINE bool dist_cluster_test(float* best2, int64_t* ia, int64_t* ib, const dist_cluster_t* cl, const dist_set_t* b, size_t j, const dist_cell_t* c, const int mode, const bool largest,
                                           dist_near_t* nr, bool record) {
    md_256 dx = md_mm256_sub_ps(cl->x, md_mm256_set1_ps(b->x[j]));
    md_256 dy = md_mm256_sub_ps(cl->y, md_mm256_set1_ps(b->y[j]));
    md_256 dz = md_mm256_sub_ps(cl->z, md_mm256_set1_ps(b->z[j]));
    dist_reduce(&dx, &dy, &dz, c, mode);
    md_256 d2 = dist_length2(dx, dy, dz);
    if (mode == DIST_PBC_TRICLINIC) {
        // At or beyond R_SAFE a reduced length is an upper bound only: the true one lies in [R_SAFE, reduced]
        const md_256 far = md_mm256_cmpge_ps(d2, md_mm256_set1_ps(c->rsafe2));
        const int far_mask = md_mm256_movemask_ps(far);
        if (far_mask) {
            // Equal to the bound matters too: a pair at it is near
            const float bound = nr ? nr->bound2 : *best2;
            const bool matters = largest ? (md_mm256_movemask_ps(md_mm256_and_ps(far, md_mm256_cmpge_ps(d2, md_mm256_set1_ps(bound)))) != 0) : (c->rsafe2 <= bound);
            if (matters) {
                d2 = dist_search_d2(dx, dy, dz, c);
            }
        }
    }
    const float m = largest ? md_mm256_reduce_max_ps(d2) : md_mm256_reduce_min_ps(d2);
    // With near pairs, one comparison for the usual lanes, as without: not near the best is not better than it either
    if (nr && (largest ? !(m >= nr->bound2) : !(m <= nr->bound2))) return false;
    const bool better = largest ? m > *best2 : m < *best2;
    if (better) {
        const int lane = md_ctz32((uint32_t)md_mm256_movemask_ps(md_mm256_cmpeq_ps(d2, md_mm256_set1_ps(m))));
        *best2 = m;
        *ia = cl->idx[lane];
        *ib = (int64_t)j;
        if (nr) nr->bound2 = dist_near_bound2(m, nr->tol, largest);
    }
    if (nr && record) {
        const md_256 e = md_mm256_set1_ps(nr->bound2);
        const int mask = md_mm256_movemask_ps(largest ? md_mm256_cmpge_ps(d2, e) : md_mm256_cmple_ps(d2, e)) & ((1 << cl->count) - 1);
        if (mask) {
            float v[8];
            md_mm256_storeu_ps(v, d2);
            dist_near_record(nr, cl->idx[0], v, j, mask, largest);
        }
    }
    return better;
}

// The bound of the sphere test, squared: a point of b at a reduced distance from the centre of the cluster at or beyond
// it (for the farthest pair: at or within it) cannot improve on best2. Negative for the farthest pair while every
// point can.
static MD_FORCE_INLINE float dist_threshold2(float best2, const dist_cluster_t* cl, const dist_cell_t* c, const bool largest) {
    if (largest) {
        if (best2 < 0.0f) return -1.0f;
        const float t = sqrtf(best2) - cl->rad - c->eps;
        return t > 0.0f ? t * t : -1.0f;
    }
    if (best2 >= FLT_MAX) return FLT_MAX;
    const float t = sqrtf(best2) + cl->rad + c->eps;
    return t * t;
}

// Given nr, which takes one group, the culling leaves room for its tolerance and the pairs near the extreme are
// recorded in it
static MD_FORCE_INLINE void dist_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const md_float3_t* a, const size_t* a_offsets, size_t num_groups,
                                     const dist_set_t* b, const dist_cell_t* c, const int mode, const bool largest, dist_near_t* nr) {
    int64_t prev = -1;
    // In locals: what the near pairs' stores could alias as far as the compiler knows, it reloads every iteration.
    // The cell too, whose address goes nowhere the stores do.
    const dist_cell_t cell_local = *c;
    c = &cell_local;
    const float* bx = b->x;
    const float* by = b->y;
    const float* bz = b->z;
    const size_t bn = b->n;
    for (size_t g = 0; g < num_groups; ++g) {
        const size_t beg = a_offsets[g];
        const size_t end = a_offsets[g + 1];
        float   best2 = largest ? -1.0f : FLT_MAX;
        int64_t ia = -1;
        int64_t ib = -1;
        bool done = false;
        if (nr) nr->bound2 = dist_near_bound2(best2, nr->tol, largest);

        for (size_t i = beg; i < end && !done; i += 8) {
            dist_cluster_t cl;
            dist_cluster_init(&cl, a, i, MD_MIN(8, end - i), c, mode);

            // The partner of the previous group: a real pair, usually close to the best one. Not recorded: the scan
            // below meets it again if it is near.
            if (prev >= 0) {
                dist_cluster_test(&best2, &ia, &ib, &cl, b, (size_t)prev, c, mode, largest, nr, false);
            }

            const md_256 cx = md_mm256_set1_ps(cl.cen[0]);
            const md_256 cy = md_mm256_set1_ps(cl.cen[1]);
            const md_256 cz = md_mm256_set1_ps(cl.cen[2]);
            float  thr2 = dist_threshold2(nr ? nr->bound2 : best2, &cl, c, largest);
            md_256 vthr = md_mm256_set1_ps(thr2);
            const md_256 vrsafe2 = md_mm256_set1_ps(c->rsafe2);

            for (size_t j = 0; j < bn && !done; j += 8) {
                md_256 dx = md_mm256_sub_ps(cx, md_mm256_loadu_ps(bx + j));
                md_256 dy = md_mm256_sub_ps(cy, md_mm256_loadu_ps(by + j));
                md_256 dz = md_mm256_sub_ps(cz, md_mm256_loadu_ps(bz + j));
                dist_reduce(&dx, &dy, &dz, c, mode);
                const md_256 d2 = dist_length2(dx, dy, dz);
                int mask = md_mm256_movemask_ps(largest ? md_mm256_cmpgt_ps(d2, vthr) : md_mm256_cmplt_ps(d2, vthr));
                if (mode == DIST_PBC_TRICLINIC && !largest && thr2 > c->rsafe2) {
                    // The sphere reaches past R_SAFE: a reduced distance there may overstate the true one
                    const int far = md_mm256_movemask_ps(md_mm256_cmpge_ps(d2, vrsafe2)) & ~mask;
                    if (far) {
                        mask |= md_mm256_movemask_ps(md_mm256_cmplt_ps(dist_search_d2(dx, dy, dz, c), vthr));
                    }
                }
                if (j + 8 > bn) {
                    mask &= (1 << (bn - j)) - 1;
                }
                while (mask) {
                    const size_t k = j + md_ctz32((uint32_t)mask);
                    mask &= mask - 1;
                    if (dist_cluster_test(&best2, &ia, &ib, &cl, b, k, c, mode, largest, nr, true)) {
                        if (!largest && best2 <= 0.0f && !nr) {
                            // Nothing is nearer than a shared point
                            done = true;
                            break;
                        }
                        thr2 = dist_threshold2(nr ? nr->bound2 : best2, &cl, c, largest);
                        vthr = md_mm256_set1_ps(thr2);
                    }
                }
            }
        }

        out_dist[g] = ib >= 0 ? sqrtf(MD_MAX(best2, 0.0f)) : 0.0f;
        if (out_idx_a) out_idx_a[g] = ia;
        if (out_idx_b) out_idx_b[g] = ib;
        if (ib >= 0) prev = ib;
        if (nr) dist_near_prune(nr, largest);
    }
}

// One function per kind of cell and extreme, so each loop is compiled for it; the plain ones without near pairs
static void dist_groups_none_min(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { (void)nr; dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_NONE, false, NULL); }
static void dist_groups_ortho_min(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { (void)nr; dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_ORTHO, false, NULL); }
static void dist_groups_tri_min(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { (void)nr; dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_TRICLINIC, false, NULL); }
static void dist_groups_none_max(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { (void)nr; dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_NONE, true, NULL); }
static void dist_groups_ortho_max(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { (void)nr; dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_ORTHO, true, NULL); }
static void dist_groups_tri_max(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { (void)nr; dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_TRICLINIC, true, NULL); }
static void dist_near_none_min(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_NONE, false, nr); }
static void dist_near_ortho_min(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_ORTHO, false, nr); }
static void dist_near_tri_min(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_TRICLINIC, false, nr); }
static void dist_near_none_max(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_NONE, true, nr); }
static void dist_near_ortho_max(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_ORTHO, true, nr); }
static void dist_near_tri_max(float* d, int64_t* ia, int64_t* ib, const md_float3_t* a, const size_t* o, size_t n, const dist_set_t* b, const dist_cell_t* c, dist_near_t* nr) { dist_groups(d, ia, ib, a, o, n, b, c, DIST_PBC_TRICLINIC, true, nr); }

static void distance_extent_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const md_float3_t* coord_a, const size_t* a_offsets, size_t num_groups, const md_float3_t* coord_b, size_t num_b, const md_unitcell_t* cell, bool largest, dist_near_t* nr) {
    MD_ASSERT(out_dist);
    MD_ASSERT(a_offsets || num_groups == 0);
    if (num_groups == 0) return;

    if (num_b == 0) {
        for (size_t g = 0; g < num_groups; ++g) {
            out_dist[g] = 0.0f;
            if (out_idx_a) out_idx_a[g] = -1;
            if (out_idx_b) out_idx_b[g] = -1;
        }
        return;
    }
    MD_ASSERT(coord_b);

    md_temp_scope_t temp = md_temp_begin();

    // b, one array per axis. The lanes past the end repeat the last point and are masked off.
    const size_t cap = MD_ALIGN_TO(num_b, 8);
    float* bx = md_temp_alloc_array(temp, float, cap);
    float* by = md_temp_alloc_array(temp, float, cap);
    float* bz = md_temp_alloc_array(temp, float, cap);
    float scale = 0.0f;
    for (size_t i = 0; i < cap; ++i) {
        const md_float3_t p = coord_b[MD_MIN(i, num_b - 1)];
        bx[i] = p.x;
        by[i] = p.y;
        bz[i] = p.z;
        scale = MD_MAX(scale, MD_MAX(fabsf(p.x), MD_MAX(fabsf(p.y), fabsf(p.z))));
    }
    for (size_t g = 0; g < num_groups; ++g) {
        for (size_t i = a_offsets[g]; i < a_offsets[g + 1]; ++i) {
            const md_float3_t p = coord_a[i];
            scale = MD_MAX(scale, MD_MAX(fabsf(p.x), MD_MAX(fabsf(p.y), fabsf(p.z))));
        }
    }

    dist_cell_t c;
    dist_cell_init(&c, cell, scale);
    const dist_set_t b = { bx, by, bz, num_b };

    typedef void (*dist_fn)(float*, int64_t*, int64_t*, const md_float3_t*, const size_t*, size_t, const dist_set_t*, const dist_cell_t*, dist_near_t*);
    static const dist_fn fns[2][2][3] = {
        { { dist_groups_none_min, dist_groups_ortho_min, dist_groups_tri_min }, { dist_groups_none_max, dist_groups_ortho_max, dist_groups_tri_max } },
        { { dist_near_none_min,   dist_near_ortho_min,   dist_near_tri_min   }, { dist_near_none_max,   dist_near_ortho_max,   dist_near_tri_max   } },
    };
    MD_ASSERT(c.mode >= DIST_PBC_NONE && c.mode <= DIST_PBC_TRICLINIC);
    fns[nr != NULL][largest][c.mode](out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, &b, &c, nr);

    md_temp_end(temp);
}

void md_geom_min_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const md_float3_t* coord_a, const size_t* a_offsets, size_t num_groups, const md_float3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    distance_extent_groups(out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, coord_b, num_b, cell, false, NULL);
}

void md_geom_max_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const md_float3_t* coord_a, const size_t* a_offsets, size_t num_groups, const md_float3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    distance_extent_groups(out_dist, out_idx_a, out_idx_b, coord_a, a_offsets, num_groups, coord_b, num_b, cell, true, NULL);
}

static size_t distance_extent_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const md_float3_t* coord_a, size_t num_a, const md_float3_t* coord_b, size_t num_b, const md_unitcell_t* cell, bool largest) {
    if (num_a == 0 || num_b == 0) return 0;
    MD_ASSERT(out_idx_a && out_idx_b);
    md_temp_scope_t temp = md_temp_begin();
    dist_near_t nr = {
        .tol = MD_MAX(tol, 0.0f),
        .a = out_idx_a,
        .b = out_idx_b,
        .d2 = md_temp_alloc_array(temp, float, MD_MAX(cap, 1)),
        .cap = cap,
    };
    const size_t offsets[2] = { 0, num_a };
    float   dist = 0.0f;
    int64_t ia = -1;
    int64_t ib = -1;
    distance_extent_groups(&dist, &ia, &ib, coord_a, offsets, 1, coord_b, num_b, cell, largest, &nr);
    md_temp_end(temp);
    return nr.overflow ? cap + 1 : nr.count;
}

size_t md_geom_min_distance_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const md_float3_t* coord_a, size_t num_a, const md_float3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    return distance_extent_near(out_idx_a, out_idx_b, cap, tol, coord_a, num_a, coord_b, num_b, cell, false);
}

size_t md_geom_max_distance_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const md_float3_t* coord_a, size_t num_a, const md_float3_t* coord_b, size_t num_b, const md_unitcell_t* cell) {
    return distance_extent_near(out_idx_a, out_idx_b, cap, tol, coord_a, num_a, coord_b, num_b, cell, true);
}

#if MD_COMPILER_MSVC
#   pragma warning( pop )
#endif

static inline void min_image_ortho(float dx[3], const float ext[3], const float inv_ext[3]) {
    dx[0] -= ext[0] * nearbyintf(dx[0] * inv_ext[0]);
    dx[1] -= ext[1] * nearbyintf(dx[1] * inv_ext[1]);
    dx[2] -= ext[2] * nearbyintf(dx[2] * inv_ext[2]);
}

// Reciprocal extent for min_image_ortho.
// A zero extent means that axis is not periodic; a zero reciprocal makes the minimum image a no-op
// along it rather than producing a NaN.
static inline void min_image_inv_ext(float out_inv_ext[3], const float ext[3]) {
    out_inv_ext[0] = ext[0] != 0.0f ? 1.0f / ext[0] : 0.0f;
    out_inv_ext[1] = ext[1] != 0.0f ? 1.0f / ext[1] : 0.0f;
    out_inv_ext[2] = ext[2] != 0.0f ? 1.0f / ext[2] : 0.0f;
}

void md_geom_min_image(float* dx, size_t stride, size_t count, const md_unitcell_t* cell) {
    MD_ASSERT(dx || count == 0);
    if (!cell) return;
    if (stride == 0) stride = 3 * sizeof(float);
    char* base = (char*)dx;
    const uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        float ext[3] = { 0 };
        float inv_ext[3];
        md_unitcell_diag_extract_float(ext, cell);
        min_image_inv_ext(inv_ext, ext);
        for (size_t i = 0; i < count; ++i) {
            min_image_ortho((float*)(base + i * stride), ext, inv_ext);
        }
    } else if (flags & MD_UNITCELL_TRICLINIC) {
        float A[3][3] = { 0 };
        md_unitcell_A_extract_float(A, cell);
        for (size_t i = 0; i < count; ++i) {
            min_image_triclinic((float*)(base + i * stride), MD_AS_CONST_MAT3(A));
        }
    }
}

// Wraps into [0, ext) along each axis with a period; an axis without one (ext zero) is left alone.
// Packed without an index the axes are independent and it runs on the floats directly.
static void pbc_ortho(float* xyz, size_t stride, const int32_t* indices, size_t count, md_float3_t box_ext) {
    MD_ASSERT(xyz);

    const md_float4_t ext = md_float4_from_float3(box_ext, 0);
    const md_float4_t ref = md_float4_mul1(ext, 0.5f);

    size_t i = 0;
    if (!indices && stride == 3 * sizeof(float)) {
        md_xv v_ext[XYZ_GROUP_VECS], v_inv[XYZ_GROUP_VECS], v_ref[XYZ_GROUP_VECS], v_refi[XYZ_GROUP_VECS], v_keep[XYZ_GROUP_VECS];
        xyz_lane_pattern(v_ext, ext.x, ext.y, ext.z);
        xyz_lane_pattern(v_inv, ext.x > 0 ? 1.0f / ext.x : 0.0f, ext.y > 0 ? 1.0f / ext.y : 0.0f, ext.z > 0 ? 1.0f / ext.z : 0.0f);
        xyz_lane_pattern(v_ref, ref.x, ref.y, ref.z);
        for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
            v_keep[j] = md_xv_cmpeq_ps(v_ext[j], md_xv_setzero_ps());
            v_refi[j] = md_xv_mul_ps(v_ref[j], v_inv[j]);
        }

        // Scalar up to alignment (below), then eight atoms at a time
        const size_t head = xyz_align_peel(xyz, count);
        for (; i < head; ++i) {
            float* p = xyz + i * 3;
            const md_float4_t pos = md_float4_deperiodize_ortho(md_float4_set(p[0], p[1], p[2], 0), ref, ext);
            p[0] = pos.x; p[1] = pos.y; p[2] = pos.z;
        }
        float* f = xyz;
        for (; i + 8 <= count; i += 8) {
            for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
                const md_xv v = md_xv_loadu_ps(f + i * 3 + j * MD_XV_WIDTH);
                // Into [0, ext]: d from a fused multiply-subtract, see the deperiodize comment in md_simd.h
                md_xv d = md_xv_fmsub_ps(v, v_inv[j], v_refi[j]);
                d = md_xv_sub_ps(d, md_xv_round_ps(d));
                const md_xv p = md_xv_fmadd_ps(d, v_ext[j], v_ref[j]);
                md_xv_storeu_ps(f + i * 3 + j * MD_XV_WIDTH, md_xv_select_ps(p, v, v_keep[j]));
            }
        }
    }
    for (; i < count; ++i) {
        const size_t idx = indices ? (size_t)indices[i] : i;
        float* p = (float*)((char*)xyz + idx * stride);
        const md_float4_t pos = md_float4_deperiodize_ortho(md_float4_set(p[0], p[1], p[2], 0), ref, ext);
        p[0] = pos.x; p[1] = pos.y; p[2] = pos.z;
    }
}

// Into the cell: to fractional coordinates, the fraction of each, and back.
static void pbc_triclinic(float* xyz, size_t stride, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    MD_ASSERT(xyz);

    md_float3x3_t A = { 0 };
    md_unitcell_A_extract_float(A.elem, cell);

    md_float3x3_t I = { 0 };
    md_unitcell_I_extract_float(I.elem, cell);

    size_t i = 0;
    if (!indices && stride == 3 * sizeof(float)) {
        md_256 vA[3][3], vI[3][3];
        for (int c = 0; c < 3; ++c) {
            for (int r = 0; r < 3; ++r) {
                vA[c][r] = md_mm256_set1_ps(A.elem[c][r]);
                vI[c][r] = md_mm256_set1_ps(I.elem[c][r]);
            }
        }
        for (; i + 8 <= count; i += 8) {
            md_256 x, y, z;
            md_mm256_load_xyz_packed_ps(&x, &y, &z, xyz + i * 3);
            md_256 f[3];
            for (int r = 0; r < 3; ++r) {
                f[r] = md_mm256_fmadd_ps(vI[0][r], x, md_mm256_fmadd_ps(vI[1][r], y, md_mm256_mul_ps(vI[2][r], z)));
                f[r] = md_mm256_fract_ps(f[r]);
            }
            md_256 c[3];
            for (int r = 0; r < 3; ++r) {
                c[r] = md_mm256_fmadd_ps(vA[0][r], f[0], md_mm256_fmadd_ps(vA[1][r], f[1], md_mm256_mul_ps(vA[2][r], f[2])));
            }
            md_mm256_store_xyz_packed_ps(xyz + i * 3, c[0], c[1], c[2]);
        }
    }

    const md_float4x3_t I4x3 = md_float4x3_from_float3x3(I);
    const md_float4x3_t A4x3 = md_float4x3_from_float3x3(A);
    for (; i < count; ++i) {
        const size_t idx = indices ? (size_t)indices[i] : i;
        float* p = (float*)((char*)xyz + idx * stride);
        md_float4_t c = md_float4x3_mul_float4(I4x3, md_float4_set(p[0], p[1], p[2], 0));
        c = md_float4_fract(c);
        c = md_float4x3_mul_float4(A4x3, c);
        p[0] = c.x; p[1] = c.y; p[2] = c.z;
    }
}

bool md_geom_pbc(float* xyz, size_t stride, const int32_t* indices, size_t count, const md_unitcell_t* cell) {
    if (!xyz) {
        MD_LOG_ERROR("Missing required input: xyz");
        return false;
    }

    if (!cell) {
        MD_LOG_ERROR("Missing unitcell");
        return false;
    }

    if (stride == 0) stride = 3 * sizeof(float);
    uint32_t flags = md_unitcell_flags(cell);
    if (flags & MD_UNITCELL_ORTHO) {
        md_float3_t ext = { 0 };
        md_unitcell_diag_extract_float(ext.elem, cell);
        pbc_ortho(xyz, stride, indices, count, ext);
        return true;
    } else if (flags & MD_UNITCELL_TRICLINIC) {
        pbc_triclinic(xyz, stride, indices, count, cell);
        return true;
    }

    // The unit cell is not initialized or is simply not periodic
    MD_LOG_ERROR("Unrecognized unit_cell type");
    return false;
}

void md_geom_unwrap_molecule(md_system_state_t* state, const md_molecule_t* mol) {
    MD_ASSERT(state);
    MD_ASSERT(mol);

    if (mol->count == 0) {
        return;
    }

    const bool is_ortho     = md_unitcell_is_orthorhombic(&state->unitcell);
    const bool is_triclinic = md_unitcell_is_triclinic(&state->unitcell);
    if (!is_ortho && !is_triclinic) {
        // Nothing to unwrap against
        return;
    }

    float ext[3] = { 0 };
    float inv_ext[3] = { 0 };
    float box[3][3] = { 0 };

    if (is_ortho) {
        md_unitcell_diag_extract_float(ext, &state->unitcell);
        min_image_inv_ext(inv_ext, ext);
    } else {
        md_unitcell_A_extract_float(box, &state->unitcell);
    }

    // md_molecule_data_t stores each molecule in BFS order, so a single forward pass places every
    // atom after the one it was reached from. No queue, no visited set and no bond lookups: the
    // traversal was already computed by md_molecule_infer.
    // The seed is its own parent, so it deperiodizes against itself and stays put - no branch.
    for (size_t i = 0; i < mol->count; ++i) {
        const int32_t a = mol->atom_idx[i];
        const int32_t p = mol->parent_idx[i];

        const md_float3_t x = state->xyz[a];
        const md_float3_t r = state->xyz[p];
	    md_float3_t dx = {x.x - r.x, x.y - r.y, x.z - r.z};

        if (is_ortho) {
            min_image_ortho(dx.elem, ext, inv_ext);
        } else {
            min_image_triclinic(dx.elem, MD_AS_CONST_MAT3(box));
        }

        // The atom is placed relative to its PARENT, not relative to where it happened to sit:
        // dx is already the minimum image of (x - r), so the unwrapped position is r + dx.
        // Adding dx to x instead yields 2x - r, which is wrong even when no image correction applies.
        state->xyz[a] = md_float3_add(r, dx);
    }
}

void md_geom_unwrap_system(md_system_state_t* state, const md_system_t* sys) {
    MD_ASSERT(state);
    MD_ASSERT(sys);

    for (size_t i = 0; i < sys->molecule.count; ++i) {
        const md_molecule_t mol = md_system_molecule(sys, i);
        md_geom_unwrap_molecule(state, &mol);
    }
}

// PERIODIC IMAGE SELECTION
//
// Everything below shares one idea: the estimator picks the periodic image of each point, choosing
// whichever image minimises the estimator's own objective. That is strictly better than deciding the
// images up front and estimating afterwards, because a preprocessing pass has to commit using some
// criterion unrelated to what is about to be measured.
//
// The payoff is stability. Every image choice has a margin - how far the deciding quantity sits from
// the half cell flip boundary - and the margin is what determines whether a point flips images
// between two nearly identical frames. Deciding by distance to a centre gives an atom at radius r a
// margin of only (L/2 - r), and it is worst for the outer atoms which dominate a rotational fit.
// Deciding by fit residual gives a margin of nearly L/2, because minimising the residual is the
// whole point. The method that is most accurate here is also the one that is hardest to flip.

// Decoded cell, hoisted out of the inner loops.
typedef struct pbc_image_t {
    bool   periodic;
    bool   ortho;
    md_float4_t ext;   // orthorhombic extent, w = 0 so the weight lane is carried through untouched
    md_float3x3_t A;     // triclinic basis
    md_float3x3_t I;     // triclinic inverse basis, Cartesian -> fractional
} pbc_image_t;

static pbc_image_t pbc_image_init(const md_unitcell_t* cell) {
    pbc_image_t im = {0};
    if (!cell) {
        return im;
    }
    if (md_unitcell_is_orthorhombic(cell)) {
        md_float3_t ext = {0};
        md_unitcell_diag_extract_float(ext.elem, cell);
        im.periodic = true;
        im.ortho    = true;
        im.ext      = md_float4_from_float3(ext, 0);
    } else if (md_unitcell_is_triclinic(cell)) {
        md_unitcell_A_extract_float(im.A.elem, cell);
        md_unitcell_I_extract_float(im.I.elem, cell);
        im.periodic = true;
    }
    return im;
}

// The image of 'pos' nearest to 'ref'. The weight in w is carried through untouched.
static inline md_float4_t pbc_image_nearest(const pbc_image_t* im, md_float4_t pos, md_float3_t ref) {
    if (!im->periodic) {
        return pos;
    }
    if (im->ortho) {
        return md_float4_deperiodize_ortho(pos, md_float4_from_float3(ref, 0), im->ext);
    }
    md_float4_t out = pos;
    deperiodize_triclinic(out.elem, ref.elem, im->A.elem);
    return out;
}

// The lattice vector nearest to 'delta'. An axis with no extent is not periodic and contributes
// nothing.
static inline md_float3_t pbc_image_lattice_round(const pbc_image_t* im, md_float3_t delta) {
    if (!im->periodic) {
        return md_float3_zero();
    }
    if (im->ortho) {
        const md_float3_t ext = md_float3_from_float4(im->ext);
        return md_float3_set(ext.x != 0.0f ? roundf(delta.x / ext.x) * ext.x : 0.0f,
                        ext.y != 0.0f ? roundf(delta.y / ext.y) * ext.y : 0.0f,
                        ext.z != 0.0f ? roundf(delta.z / ext.z) * ext.z : 0.0f);
    }
    const md_float3_t f = md_float3x3_mul_float3(im->I, delta);
    return md_float3x3_mul_float3(im->A, md_float3_set(roundf(f.x), roundf(f.y), roundf(f.z)));
}

// KEEPING THE SET IN ITS OWN IMAGE
//
// The alternations below settle on images that are consistent RELATIVE to a centre, and that centre
// is seeded from the circular mean - whose angle is taken in [0, 2pi), so it always lands in the
// reference cell no matter where the input lives. Left at that, the whole set is quietly translated
// into image (0,0,0) as a side effect of being made whole.
//
// That is invisible to a caller that only reads relative geometry (an rmsd, a covariance, a
// rotation) and wrong for one that mixes the result back in with coordinates it did not hand over:
// a transform applied to the full system, an AABB queried against raw positions, a plane offset
// drawn over the atoms it was fitted to. Those callers get a result a whole cell vector off, and
// with a rotation in play, off by R times a cell vector, which is not even an image.
//
// So the contract is: make the set mutually consistent and otherwise leave it where it came in.
// The shift back is the lattice vector closest to the displacement of the weighted mean, i.e. the
// one that moves the input least. A set that genuinely straddles a boundary has no image to
// preserve - the two candidates are one lattice vector apart and equally good - and a caller that
// needs the set inside the cell wraps afterwards, explicitly.
static void pbc_image_restore(const pbc_image_t* im, md_float4_t* xyzw, size_t count, md_float3_t input_com, md_float3_t* in_out_com) {
    const md_float3_t off = pbc_image_lattice_round(im, md_float3_sub(input_com, *in_out_com));
    if (md_float3_dot(off, off) == 0.0f) {
        return;
    }
    const md_float4_t off4 = md_float4_from_float3(off, 0);
    for (size_t i = 0; i < count; ++i) {
        xyzw[i] = md_float4_add(xyzw[i], off4);
    }
    *in_out_com = md_float3_add(*in_out_com, off);
}

// Squared distance below which two placements count as the same image. Generous: a genuine image
// change moves a point by a whole cell vector, so anything smaller is float noise.
#define PBC_IMAGE_EPS_SQ 1.0e-6f

// Alternation is a fixed point iteration and settles almost immediately - the assignment is usually
// stable on the second pass. The cap only bounds pathological input.
#define PBC_MAX_ITER 8

bool md_geom_deperiodize_self(md_float4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell, md_float3_t* out_com) {
    if (!in_out_xyzw) {
        MD_LOG_ERROR("Missing required input: in_out_xyzw");
        return false;
    }

    if (count == 0) {
        if (out_com) *out_com = md_float3_zero();
        return true;
    }

    const pbc_image_t im = pbc_image_init(cell);

    if (!im.periodic) {
        if (out_com) *out_com = md_geom_com((const float*)in_out_xyzw, sizeof(md_float4_t), (const float*)in_out_xyzw + 3, sizeof(md_float4_t), NULL, count, NULL);
        return true;
    }

    // Where the input sits, before anything is moved. Only the image it implies is used - see
    // pbc_image_restore.
    const md_float3_t input_com = md_geom_com((const float*)in_out_xyzw, sizeof(md_float4_t), (const float*)in_out_xyzw + 3, sizeof(md_float4_t), NULL, count, NULL);

    // Seed with the circular mean. It needs no image assignment to compute, so unlike an arbitrary
    // reference point it cannot be biased by whichever images the coordinates happened to arrive in.
    md_float3_t com = md_geom_com((const float*)in_out_xyzw, sizeof(md_float4_t), (const float*)in_out_xyzw + 3, sizeof(md_float4_t), NULL, count, cell);

    for (int iter = 0; iter < PBC_MAX_ITER; ++iter) {
        bool stable = true;
        for (size_t i = 0; i < count; ++i) {
            const md_float4_t before = in_out_xyzw[i];
            const md_float4_t after  = pbc_image_nearest(&im, before, com);
            if (md_float3_distance_squared(md_float3_from_float4(before), md_float3_from_float4(after)) > PBC_IMAGE_EPS_SQ) {
                stable = false;
            }
            in_out_xyzw[i] = after;
        }

        // Recompute the centre from the placed points. Both halves decrease the same objective.
        com = md_geom_com((const float*)in_out_xyzw, sizeof(md_float4_t), (const float*)in_out_xyzw + 3, sizeof(md_float4_t), NULL, count, NULL);

        if (stable) {
            break;
        }
    }

    // Whole, and back in the image it arrived in.
    pbc_image_restore(&im, in_out_xyzw, count, input_com, &com);

    if (out_com) *out_com = com;
    return true;
}

float md_geom_optimal_rotation_pbc(md_float3x3_t* out_rot, md_float3_t* out_com, md_float4_t* out_trg_xyzw,
                                             const md_float4_t* ref_xyzw, md_float3_t ref_com,
                                             const md_float4_t* trg_xyzw, size_t count,
                                             const md_unitcell_t* cell, int max_iter, float tol) {
    MD_ASSERT(ref_xyzw);
    MD_ASSERT(trg_xyzw);

    if (out_rot) *out_rot = md_float3x3_ident();
    if (out_com) *out_com = md_float3_zero();

    if (count == 0) {
        return 0.0f;
    }

    const pbc_image_t im = pbc_image_init(cell);

    md_temp_scope_t temp = md_temp_begin();

    // The caller can have the placed points if it wants them (rmsd does); otherwise they are scratch.
    md_float4_t* cur = out_trg_xyzw ? out_trg_xyzw : md_temp_alloc_array(temp, md_float4_t, count);

    // Where the target sits before anything is moved. Read now: 'cur' is allowed to alias trg_xyzw
    // and the loop below writes through it.
    const md_float3_t input_com = md_geom_com((const float*)trg_xyzw, sizeof(md_float4_t), (const float*)trg_xyzw + 3, sizeof(md_float4_t), NULL, count, NULL);

    // Image free seed, safe to compute before any assignment exists.
    md_float3_t com = md_geom_com((const float*)trg_xyzw, sizeof(md_float4_t), (const float*)trg_xyzw + 3, sizeof(md_float4_t), NULL, count, im.periodic ? cell : NULL);
    md_float3x3_t rot = md_float3x3_ident();

    // Alternating minimisation of  sum_k w_k |R (q_k + A n_k - c) - p_k|^2  over R, c and n_k.
    // Each block has a closed form optimum given the other, so the objective decreases monotonically
    // and the discrete part reaches a genuine fixed point rather than a tolerance.
    bool stable = false;
    for (int iter = 0; iter < max_iter && !stable; ++iter) {
        // rot maps the target frame onto the reference frame, so its transpose carries the reference
        // template back out into the target's frame, which is where the prediction has to live.
        const md_float3x3_t rot_inv = md_float3x3_transpose(rot);

        stable = (iter > 0);
        for (size_t i = 0; i < count; ++i) {
            const md_float3_t p    = md_float3_sub(md_float3_from_float4(ref_xyzw[i]), ref_com);
            const md_float3_t pred = md_float3_add(com, md_float3x3_mul_float3(rot_inv, p));
            const md_float4_t placed = pbc_image_nearest(&im, trg_xyzw[i], pred);

            if (iter > 0 && md_float3_distance_squared(md_float3_from_float4(cur[i]), md_float3_from_float4(placed)) > tol) {
                stable = false;
            }
            cur[i] = placed;
        }

        // Refit. Both are the exact optimum for this assignment.
        com = md_geom_com((const float*)cur, sizeof(md_float4_t), (const float*)cur + 3, sizeof(md_float4_t), NULL, count, NULL);

        const md_float4_t* const sets[2] = { ref_xyzw, cur };
        const md_float3_t        coms[2] = { ref_com,  com };
        rot = md_float3x3_extract_rotation(md_float3x3_cross_covariance_matrix_float4(sets, NULL, count, coms));
    }

    // Back in the image the target arrived in. The rotation and the residual below are measured
    // relative to 'com', so both are untouched by this.
    pbc_image_restore(&im, cur, count, input_com, &com);

    // The largest per point residual is the distance to an ambiguous image choice. Small compared to
    // half a cell means the assignment above is not close to flipping.
    float max_dist_sq = 0.0f;
    for (size_t i = 0; i < count; ++i) {
        const md_float3_t p = md_float3_sub(md_float3_from_float4(ref_xyzw[i]), ref_com);
        const md_float3_t q = md_float3x3_mul_float3(rot, md_float3_sub(md_float3_from_float4(cur[i]), com));
        const float  d = md_float3_distance_squared(p, q);
        max_dist_sq = MD_MAX(max_dist_sq, d);
    }

    md_temp_end(temp);

    if (out_rot) *out_rot = rot;
    if (out_com) *out_com = com;

    return sqrtf(max_dist_sq);
}

md_float3x3_t md_geom_optimal_rotation_rel(const md_float4_t* ref_rel_xyzw, const md_float4_t* trg_rel_xyzw, size_t count) {
    MD_ASSERT(ref_rel_xyzw);
    MD_ASSERT(trg_rel_xyzw);

    if (count == 0) {
        return md_float3x3_ident();
    }

    double A[3][3] = { 0 };
    double w_sum = 0.0;

    for (size_t i = 0; i < count; ++i) {
		const md_float4_t p = ref_rel_xyzw[i];
		const md_float4_t q = trg_rel_xyzw[i];

        // The question here is how to combine the weights.
        // For now we just take the average.
        // This should be equivalent to the other case where the weights for both sets are equal.
        const float w = (p.w + q.w) * 0.5f;

        A[0][0] += w * p.x * q.x;
        A[0][1] += w * p.x * q.y;
        A[0][2] += w * p.x * q.z;
        A[1][0] += w * p.y * q.x;
        A[1][1] += w * p.y * q.y;
        A[1][2] += w * p.y * q.z;
        A[2][0] += w * p.z * q.x;
        A[2][1] += w * p.z * q.y;
        A[2][2] += w * p.z * q.z;
        w_sum += w;
    }

	if (w_sum == 0.0) {
		return md_float3x3_ident();
	}

    md_float3x3_t C;
	C.elem[0][0] = (float)(A[0][0] / w_sum);
    C.elem[0][1] = (float)(A[0][1] / w_sum);
	C.elem[0][2] = (float)(A[0][2] / w_sum);
	C.elem[1][0] = (float)(A[1][0] / w_sum);
	C.elem[1][1] = (float)(A[1][1] / w_sum);
	C.elem[1][2] = (float)(A[1][2] / w_sum);
	C.elem[2][0] = (float)(A[2][0] / w_sum);
	C.elem[2][1] = (float)(A[2][1] / w_sum);
	C.elem[2][2] = (float)(A[2][2] / w_sum);

    return md_float3x3_extract_rotation(C);
}

double md_geom_rmsd(const md_float4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const md_float3_t in_com[2]) {
    const md_float3x3_t R = md_float3x3_optimal_rotation_float4(in_xyzw, in_idx, count, in_com);
    double d_sum = 0;
    double w_sum = 0;
    if (in_idx) {
        for (size_t i = 0; i < count; ++i) {
            int32_t idx0 = in_idx[0][i];
            int32_t idx1 = in_idx[1][i];
            md_float3_t u = {in_xyzw[0][idx0].x - in_com[0].x, in_xyzw[0][idx0].y - in_com[0].y, in_xyzw[0][idx0].z - in_com[0].z};
            md_float3_t v = {in_xyzw[1][idx1].x - in_com[1].x, in_xyzw[1][idx1].y - in_com[1].y, in_xyzw[1][idx1].z - in_com[1].z};
            md_float3_t vp = md_float3x3_mul_float3(R, v);
            md_float3_t d = md_float3_sub(u, vp);
            float w = (in_xyzw[0][idx0].w + in_xyzw[1][idx1].w) * 0.5f;
            d_sum += w * md_float3_dot(d, d);
            w_sum += w;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            md_float3_t u = {in_xyzw[0][i].x - in_com[0].x, in_xyzw[0][i].y - in_com[0].y, in_xyzw[0][i].z - in_com[0].z};
            md_float3_t v = {in_xyzw[1][i].x - in_com[1].x, in_xyzw[1][i].y - in_com[1].y, in_xyzw[1][i].z - in_com[1].z};
            md_float3_t vp = md_float3x3_mul_float3(R, v);
            md_float3_t d = md_float3_sub(u, vp);
            float w = (in_xyzw[0][i].w + in_xyzw[1][i].w) * 0.5f;
            d_sum += w * md_float3_dot(d, d);
            w_sum += w;
        }
    }

    return sqrt(d_sum / w_sum);
}

md_float3_t md_geom_shape_weights(const md_float3x3_t* covariance_matrix) {
    MD_ASSERT(covariance_matrix);
    md_float3x3_eigen_t eigen = md_float3x3_eigen(*covariance_matrix);
    float scl = 1.0f / (eigen.values.x + eigen.values.y + eigen.values.z);
    md_float3_t weights = {(eigen.values.x - eigen.values.y) * scl, 2.0f * (eigen.values.y - eigen.values.z) * scl, 3.0f * eigen.values.z * scl};
    return weights;
}

// One group of eight atoms, packed floats as they lie. The second frame is brought to the minimum
// image of the first when the cell is periodic.
static inline void lerp_xyz8(float* dst, const float* a, const float* b, const md_xv ext[XYZ_GROUP_VECS], const md_xv inv[XYZ_GROUP_VECS], bool periodic, float t) {
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
        const md_xv x0 = md_xv_loadu_ps(a + j * MD_XV_WIDTH);
        md_xv x1 = md_xv_loadu_ps(b + j * MD_XV_WIDTH);
        if (periodic) x1 = xv_deperiodize_ortho(x1, x0, ext[j], inv[j]);
        md_xv_storeu_ps(dst + j * MD_XV_WIDTH, md_xv_fmadd_ps(md_xv_sub_ps(x1, x0), md_xv_set1_ps(t), x0));
    }
}

// One group of eight atoms over four frames, weighted (see md_geom_interpolate_cubic_weighted). The weights are applied
// in order, so a weight of exactly one with the others zero gives that frame back exactly: an interpolating scheme
// lands on its frames.
static inline void cubic_xyz8(float* dst, const float* const src[4], const md_xv ext[XYZ_GROUP_VECS], const md_xv inv[XYZ_GROUP_VECS], bool periodic, const md_xv w[4]) {
    for (int j = 0; j < XYZ_GROUP_VECS; ++j) {
        md_xv x0 = md_xv_loadu_ps(src[0] + j * MD_XV_WIDTH);
        const md_xv x1 = md_xv_loadu_ps(src[1] + j * MD_XV_WIDTH);
        md_xv x2 = md_xv_loadu_ps(src[2] + j * MD_XV_WIDTH);
        md_xv x3 = md_xv_loadu_ps(src[3] + j * MD_XV_WIDTH);
        if (periodic) {
            x0 = xv_deperiodize_ortho(x0, x1, ext[j], inv[j]);
            x2 = xv_deperiodize_ortho(x2, x1, ext[j], inv[j]);
            x3 = xv_deperiodize_ortho(x3, x2, ext[j], inv[j]);
        }
        md_xv r = md_xv_mul_ps(w[0], x0);
        r = md_xv_fmadd_ps(w[1], x1, r);
        r = md_xv_fmadd_ps(w[2], x2, r);
        r = md_xv_fmadd_ps(w[3], x3, r);
        md_xv_storeu_ps(dst + j * MD_XV_WIDTH, r);
    }
}

// Linear interpolation of packed coordinates. Without a cell or in an orthorhombic one the axes are
// independent and it runs on the packed floats directly, eight atoms at a time, with a few leading
// atoms first so the loads are aligned and a partial group last, each through a small staging copy;
// any count works. A triclinic cell couples the axes: that path splits eight atoms into x, y and z
// and relies on the arrays being padded to a multiple of eight atoms.
bool md_geom_interpolate_linear(md_float3_t* out_xyz, const md_float3_t* const in_xyz[2], size_t count, const md_unitcell_t* cell, float t) {
    MD_ASSERT(out_xyz);
    MD_ASSERT(in_xyz);
    MD_ASSERT(in_xyz[0]);
    MD_ASSERT(in_xyz[1]);

    t = MD_CLAMP(t, 0.0f, 1.0f);

    const bool ortho = cell && md_unitcell_is_orthorhombic(cell);
    const bool tricl = cell && md_unitcell_is_triclinic(cell);

    md_float3x3_t A = { 0 };
    md_float3x3_t I = { 0 };
    if (ortho | tricl) {
        md_unitcell_A_extract_float(A.elem, cell);
        md_unitcell_I_extract_float(I.elem, cell);
    }

    const float* src0 = (const float*)in_xyz[0];
    const float* src1 = (const float*)in_xyz[1];
    float* dst = (float*)out_xyz;

    if (!tricl) {
        md_xv ext[XYZ_GROUP_VECS], inv[XYZ_GROUP_VECS];
        xyz_lane_pattern(ext, A.elem[0][0], A.elem[1][1], A.elem[2][2]);
        xyz_lane_pattern(inv, I.elem[0][0], I.elem[1][1], I.elem[2][2]);

        size_t i = 0;
        const size_t head = xyz_align_peel(src0, count);
        while (i < count) {
            const size_t n = (i == 0 && head) ? head : MD_MIN(8, count - i);
            if (n == 8) {
                lerp_xyz8(dst + i * 3, src0 + i * 3, src1 + i * 3, ext, inv, ortho, t);
            } else {
                float a[24] = {0}, b[24] = {0}, r[24];
                md_memcpy(a, src0 + i * 3, n * sizeof(md_float3_t));
                md_memcpy(b, src1 + i * 3, n * sizeof(md_float3_t));
                lerp_xyz8(r, a, b, ext, inv, ortho, t);
                md_memcpy(dst + i * 3, r, n * sizeof(md_float3_t));
            }
            i += n;
        }
        return true;
    }

    for (size_t i = 0; i < count; i += 8) {
        md_256 x0[3], x1[3];
        md_mm256_load_xyz_packed_ps(&x0[0], &x0[1], &x0[2], src0 + i * 3);
        md_mm256_load_xyz_packed_ps(&x1[0], &x1[1], &x1[2], src1 + i * 3);

        simd_deperiodize_triclinic(x1, x0, MD_AS_CONST_MAT3(A.elem));

        const md_256 x = md_mm256_lerp_ps(x0[0], x1[0], t);
        const md_256 y = md_mm256_lerp_ps(x0[1], x1[1], t);
        const md_256 z = md_mm256_lerp_ps(x0[2], x1[2], t);

        md_mm256_store_xyz_packed_ps(dst + i * 3, x, y, z);
    }

    return true;
}

void md_geom_cubic_spline_weights(float w[4], float t, float s) {
    MD_ASSERT(w);
    // The Hermite form, p(t) = h00 p1 + h10 m1 + h01 p2 + h11 m2, with the tangents m1 = s (p2 - p0) and m2 = s (p3 - p1)
    const float t2 = t * t;
    const float t3 = t2 * t;
    const float h00 =  2.0f * t3 - 3.0f * t2 + 1.0f;
    const float h10 =         t3 - 2.0f * t2 + t;
    const float h01 = -2.0f * t3 + 3.0f * t2;
    const float h11 =         t3 -        t2;
    w[0] = -s * h10;
    w[1] = h00 - s * h11;
    w[2] = h01 + s * h10;
    w[3] = s * h11;
}

void md_geom_cubic_bspline_weights(float w[4], float t) {
    MD_ASSERT(w);
    const float t2 = t * t;
    const float t3 = t2 * t;
    const float u  = 1.0f - t;
    w[0] = u * u * u * (1.0f / 6.0f);
    w[1] = ( 3.0f * t3 - 6.0f * t2 + 4.0f) * (1.0f / 6.0f);
    w[2] = (-3.0f * t3 + 3.0f * t2 + 3.0f * t + 1.0f) * (1.0f / 6.0f);
    w[3] = t3 * (1.0f / 6.0f);
}

// As md_geom_interpolate_linear: packed floats directly without a cell or in an orthorhombic one,
// split into x, y and z (and padded) in a triclinic one.
bool md_geom_interpolate_cubic_weighted(md_float3_t* out_xyz, const md_float3_t* const in_xyz[4], size_t count, const md_unitcell_t* cell, const float w[4]) {
    MD_ASSERT(out_xyz);
    MD_ASSERT(in_xyz);
    MD_ASSERT(in_xyz[0]);
    MD_ASSERT(in_xyz[1]);
    MD_ASSERT(in_xyz[2]);
    MD_ASSERT(in_xyz[3]);
    MD_ASSERT(w);

    const bool ortho = cell && md_unitcell_is_orthorhombic(cell);
    const bool tricl = cell && md_unitcell_is_triclinic(cell);
    md_float3x3_t A = { 0 };
    md_float3x3_t I = { 0 };
    if (ortho | tricl) {
        md_unitcell_A_extract_float(A.elem, cell);
        md_unitcell_I_extract_float(I.elem, cell);
    }

    const md_256 vw[4] = { md_mm256_set1_ps(w[0]), md_mm256_set1_ps(w[1]), md_mm256_set1_ps(w[2]), md_mm256_set1_ps(w[3]) };
    const md_xv  xw[4] = { md_xv_set1_ps(w[0]),    md_xv_set1_ps(w[1]),    md_xv_set1_ps(w[2]),    md_xv_set1_ps(w[3]) };

    const float* src[4] = { (const float*)in_xyz[0], (const float*)in_xyz[1], (const float*)in_xyz[2], (const float*)in_xyz[3] };
    float* dst = (float*)out_xyz;

    if (!tricl) {
        md_xv ext[XYZ_GROUP_VECS], inv[XYZ_GROUP_VECS];
        xyz_lane_pattern(ext, A.elem[0][0], A.elem[1][1], A.elem[2][2]);
        xyz_lane_pattern(inv, I.elem[0][0], I.elem[1][1], I.elem[2][2]);

        size_t i = 0;
        const size_t head = xyz_align_peel(src[1], count);
        while (i < count) {
            const size_t n = (i == 0 && head) ? head : MD_MIN(8, count - i);
            if (n == 8) {
                const float* const g[4] = { src[0] + i * 3, src[1] + i * 3, src[2] + i * 3, src[3] + i * 3 };
                cubic_xyz8(dst + i * 3, g, ext, inv, ortho, xw);
            } else {
                float stage[4][24] = {0}, r[24];
                for (int k = 0; k < 4; ++k) {
                    md_memcpy(stage[k], src[k] + i * 3, n * sizeof(md_float3_t));
                }
                const float* const g[4] = { stage[0], stage[1], stage[2], stage[3] };
                cubic_xyz8(r, g, ext, inv, ortho, xw);
                md_memcpy(dst + i * 3, r, n * sizeof(md_float3_t));
            }
            i += n;
        }
        return true;
    }

    for (size_t i = 0; i < count; i += 8) {
        md_256 x0[3], x1[3], x2[3], x3[3];
        md_mm256_load_xyz_packed_ps(&x0[0], &x0[1], &x0[2], src[0] + i * 3);
        md_mm256_load_xyz_packed_ps(&x1[0], &x1[1], &x1[2], src[1] + i * 3);
        md_mm256_load_xyz_packed_ps(&x2[0], &x2[1], &x2[2], src[2] + i * 3);
        md_mm256_load_xyz_packed_ps(&x3[0], &x3[1], &x3[2], src[3] + i * 3);

        simd_deperiodize_triclinic(x0, x1, MD_AS_CONST_MAT3(A.elem));
        simd_deperiodize_triclinic(x2, x1, MD_AS_CONST_MAT3(A.elem));
        simd_deperiodize_triclinic(x3, x2, MD_AS_CONST_MAT3(A.elem));

        md_256 r[3];
        for (int k = 0; k < 3; ++k) {
            r[k] = md_mm256_mul_ps(vw[0], x0[k]);
            r[k] = md_mm256_fmadd_ps(vw[1], x1[k], r[k]);
            r[k] = md_mm256_fmadd_ps(vw[2], x2[k], r[k]);
            r[k] = md_mm256_fmadd_ps(vw[3], x3[k], r[k]);
        }

        md_mm256_store_xyz_packed_ps(dst + i * 3, r[0], r[1], r[2]);
    }

    return true;
}
