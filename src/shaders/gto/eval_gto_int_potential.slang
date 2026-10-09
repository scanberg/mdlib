// eval_gto_int_potential.slang
//
// The electrostatic potential of an md_gto_int charge distribution on a 3D grid:
//
//   V(C) = sum_k q_k / |C - A_k|  +  sum_g sum_{t+u+v<=L_g} h^g_tuv R_tuv(p_g, P_g - C)
//
// One entry point per Hermite order (l0..l8) plus a resolve pass. The host sorts the gaussians by
// order and launches each order's entry point over one slab of the grid and a run of that order's
// gaussians, so every kernel contains exactly one unrolled recursion (gto_int_hermite.slang) and
// is register-sized for it: the rare high orders do not cost the common low ones occupancy.
//
// One thread per voxel. Every thread of a dispatch walks the same gaussians in the same order, so
// the gaussian and coefficient loads are uniform and there is no divergence beyond the Boys
// function's two branches.
//
// SUMMATION. Each term is computed in float and added to a 64 bit fixed point accumulator with 32
// fractional bits, kept per voxel in `accum` between dispatches. The terms cancel against the
// nuclei to a few percent of either, and a float sum over tens of thousands of them is off by
// ~5e-4 au on C60; integer addition loses nothing, is independent of how the work was split, and
// cannot be reassociated by a fast-math compiler (Metal's default). resolve adds the point charges
// the same way, converts once, and applies the md_gto_op_t.

#include "../md_gpu.slang"
#include "gto_common.slang"
#include "gto_int_hermite.slang"

#define WG_X 4
#define WG_Y 4
#define WG_Z 4

// Boys function table, see gto_int_boys_table() in md_gto_int.c: F_n(k * DT) for n < NN.
#define GTO_INT_BOYS_DT     0.1f
#define GTO_INT_BOYS_INV_DT 10.0f
#define GTO_INT_BOYS_TMAX   30.0f
#define GTO_INT_BOYS_NN     16

struct RootArgs {
    float4x4 index_to_world;   // grid index (sample offset included) -> world, bohr
    uint4    box_off;          // first voxel of the slab, xyz
    uint4    box_dim;          // voxels in the slab, xyz
    uint     first;            // first gaussian of this dispatch
    uint     count;            // gaussians in this dispatch, all of one order
    uint     coeff_first;      // index of the first coefficient of 'first'
    uint     num_points;       // resolve: point charges
    uint     operation;        // resolve: md_gto_op_t
    uint     _pad0;
    uint     _pad1;
    uint     _pad2;
    float4*  gaussians;        // (P.xyz, p)
    float*   coeffs;           // h', md_gto_int.h order, contiguous per order
    float*   boys;             // [rows][GTO_INT_BOYS_NN]
    float4*  points;           // (xyz, q)
    uint2*   accum;            // [slab voxels], fixed point (lo, hi)
    DescriptorHandle<RWTexture3D<float>> out_tex;
};
struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;
#define args (*md_root.args_ptr)

// ---------------------------------------------------------------------------
// 64 bit fixed point, 32 fractional bits, as (lo, hi) with hi two's complement.
// ---------------------------------------------------------------------------
uint2 gto_int_fixed_add(uint2 acc, float v) {
    // Magnitude into (integer, fraction); a - floor(a) is exact for a >= 0, and the fraction
    // times 2^32 stays below 2^32 (at most 2^32 - 2^8), so round() cannot overflow.
    const float a  = min(abs(v), 2147483520.0f);
    const float fl = floor(a);
    const uint  ih = uint(fl);
    const uint  il = uint(round((a - fl) * 4294967296.0f));
    uint2 r;
    if (v >= 0.0f) {
        r.x = acc.x + il;
        r.y = acc.y + ih + ((r.x < acc.x) ? 1u : 0u);
    } else {
        r.x = acc.x - il;
        r.y = acc.y - ih - ((acc.x < il) ? 1u : 0u);
    }
    return r;
}

float gto_int_fixed_to_float(uint2 acc) {
    return float(int(acc.y)) + float(acc.x) * 2.3283064365386963e-10f;
}

// ---------------------------------------------------------------------------
// s[n] = (-2p)^n F_n(T), n = 0..L
// ---------------------------------------------------------------------------
// Below TMAX: 5 term Taylor expansion of F_L around the nearest tabulated T (|dT| <= 0.05, error
// below 3e-9 relative), then downward recursion, which is stable. Above: the asymptotic F_0 and
// upward recursion, stable there for every order this kernel handles.
void gto_int_boys_scaled<let L : int>(out float s[GTO_INT_NS], float T, float p) {
    float F[GTO_INT_NS];
    const float e = exp(-T);
    if (T < GTO_INT_BOYS_TMAX) {
        const int   k = int(T * GTO_INT_BOYS_INV_DT + 0.5f);
        const float d = float(k) * GTO_INT_BOYS_DT - T;
        const uint  row = uint(k) * GTO_INT_BOYS_NN + uint(L);
        float f = args.boys[row + 4u] * (1.0f / 24.0f);
        f = args.boys[row + 3u] * (1.0f / 6.0f) + d * f;
        f = args.boys[row + 2u] * 0.5f + d * f;
        f = args.boys[row + 1u] + d * f;
        f = args.boys[row + 0u] + d * f;
        F[L] = f;
        [ForceUnroll] for (int n = L - 1; n >= 0; --n) {
            F[n] = (2.0f * T * F[n + 1] + e) * (1.0f / float(2 * n + 1));
        }
    } else {
        const float inv2T = 0.5f / T;
        F[0] = 0.5f * sqrt(3.14159265358979f / T);
        [ForceUnroll] for (int n = 0; n < L; ++n) {
            F[n + 1] = (float(2 * n + 1) * F[n] - e) * inv2T;
        }
    }
    const float m2p = -2.0f * p;
    float scale = 1.0f;
    [ForceUnroll] for (int n = 0; n <= L; ++n) {
        s[n] = scale * F[n];
        scale *= m2p;
    }
    [ForceUnroll] for (int n = L + 1; n < GTO_INT_NS; ++n) s[n] = 0.0f;
}

float gto_int_hermite<let L : int>(float3 X, float s[GTO_INT_NS], uint o) {
    if (L == 0) return gto_int_hermite_L0(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 1) return gto_int_hermite_L1(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 2) return gto_int_hermite_L2(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 3) return gto_int_hermite_L3(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 4) return gto_int_hermite_L4(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 5) return gto_int_hermite_L5(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 6) return gto_int_hermite_L6(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 7) return gto_int_hermite_L7(X.x, X.y, X.z, s, args.coeffs, o);
    return gto_int_hermite_L8(X.x, X.y, X.z, s, args.coeffs, o);
}

float3 gto_int_voxel_position(uint3 vox) {
    return mul(args.index_to_world, float4(float3(vox), 1.0f)).xyz;
}

uint gto_int_slab_index(uint3 tid) {
    return (tid.z * args.box_dim.y + tid.y) * args.box_dim.x + tid.x;
}

void gto_int_accumulate<let L : int>(uint3 tid) {
    if (any(tid >= args.box_dim.xyz)) return;
    const float3 C   = gto_int_voxel_position(args.box_off.xyz + tid);
    const uint   lin = gto_int_slab_index(tid);
    const uint   nh  = uint(((L + 1) * (L + 2) * (L + 3)) / 6);

    uint2 acc = args.accum[lin];
    for (uint i = 0u; i < args.count; ++i) {
        const float4 g = args.gaussians[args.first + i];
        const float3 X = g.xyz - C;
        float s[GTO_INT_NS];
        gto_int_boys_scaled<L>(s, g.w * dot(X, X), g.w);
        acc = gto_int_fixed_add(acc, gto_int_hermite<L>(X, s, args.coeff_first + i * nh));
    }
    args.accum[lin] = acc;
}

[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l0(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<0>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l1(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<1>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l2(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<2>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l3(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<3>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l4(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<4>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l5(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<5>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l6(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<6>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l7(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<7>(tid); }
[shader("compute")] [numthreads(WG_X, WG_Y, WG_Z)] void l8(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate<8>(tid); }

// Point charges, conversion, and the md_gto_op_t into the output texture.
[shader("compute")]
[numthreads(WG_X, WG_Y, WG_Z)]
void resolve(uint3 tid : SV_DispatchThreadID) {
    if (any(tid >= args.box_dim.xyz)) return;
    const uint3  vox = args.box_off.xyz + tid;
    const float3 C   = gto_int_voxel_position(vox);
    uint2 acc = args.accum[gto_int_slab_index(tid)];
    for (uint k = 0u; k < args.num_points; ++k) {
        const float4 q = args.points[k];
        const float3 d = C - q.xyz;
        const float r2 = dot(d, d);
        if (r2 > 1.0e-20f) acc = gto_int_fixed_add(acc, q.w * rsqrt(r2));
    }
    const float V = gto_int_fixed_to_float(acc);
    const float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
    args.out_tex[vox] = gto_apply_op(V, prev, args.operation);
}
