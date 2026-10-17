// eval_gto_int_potential_xyz.slang
//
// The electrostatic potential of an md_gto_int charge distribution at arbitrary points, read from a
// device buffer: the GPU counterpart of md_gto_int_potential_xyz (md_gto_int_gpu_potential_xyz_launch).
// The same evaluation as eval_gto_int_potential.slang on a grid - one entry point per Hermite order
// (l0..l8) over a run of that order's gaussians, then a resolve pass for the point charges - with one
// thread per point instead of per voxel, and the result written to a float buffer.
//
// The points come in the absolute frame (bohr), as float; the distribution is stored relative to
// args.ref, and the points are moved into that frame here. A point carries the precision of its
// float coordinates, as on the CPU, whose md_gto_int_potential_xyz takes float positions too.
//
// SUMMATION as in eval_gto_int_potential.slang: each term in float, summed in 64 bit fixed point.

#include "../md_gpu.slang"
#include "gto_common.slang"
#include "gto_int_hermite.slang"

#define WG_X 64

#define GTO_INT_BOYS_DT     0.1f
#define GTO_INT_BOYS_INV_DT 10.0f
#define GTO_INT_BOYS_TMAX   30.0f
#define GTO_INT_BOYS_NN     16

struct RootArgs {
    float4   ref;              // xyz: the frame the distribution is stored in, bohr
    uint     first;            // first gaussian of this dispatch
    uint     count;            // gaussians in this dispatch, all of one order
    uint     coeff_first;      // index of the first coefficient of 'first'
    uint     num_points;       // resolve: point charges
    uint     xyz_first;        // first point of this batch
    uint     xyz_count;        // points in this batch
    uint     xyz_stride;       // floats between points
    uint     _pad0;
    float*   xyz;              // the points, bohr: xyz[i * xyz_stride + 0..2]
    float4*  gaussians;        // (P.xyz, p), relative to ref
    float*   coeffs;           // h', md_gto_int.h order, contiguous per order
    float*   boys;             // [rows][GTO_INT_BOYS_NN]
    float4*  points;           // (xyz, q), relative to ref
    uint2*   accum;            // [batch points], fixed point (lo, hi)
    float*   out;              // [points], the potential, hartree / e
    uint2    _pad_end;         // the struct is 16 byte aligned (ref): its size, said out loud for the C mirror
};
struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;
#define args (*md_root.args_ptr)

#include "gto_int_potential_common.slang"

float3 gto_int_point(uint i) {
    const uint o = (args.xyz_first + i) * args.xyz_stride;
    return float3(args.xyz[o + 0u], args.xyz[o + 1u], args.xyz[o + 2u]) - args.ref.xyz;
}

void gto_int_accumulate_xyz<let L : int>(uint i) {
    if (i >= args.xyz_count) return;
    const float3 C  = gto_int_point(i);
    const uint   nh = uint(((L + 1) * (L + 2) * (L + 3)) / 6);

    uint2 acc = args.accum[i];
    for (uint k = 0u; k < args.count; ++k) {
        const float4 g = args.gaussians[args.first + k];
        const float3 X = g.xyz - C;
        float s[GTO_INT_NS];
        gto_int_boys_scaled<L>(s, g.w * dot(X, X), g.w);
        acc = gto_int_fixed_add(acc, gto_int_hermite<L>(X, s, args.coeff_first + k * nh));
    }
    args.accum[i] = acc;
}

[shader("compute")] [numthreads(WG_X, 1, 1)] void l0(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<0>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l1(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<1>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l2(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<2>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l3(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<3>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l4(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<4>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l5(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<5>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l6(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<6>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l7(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<7>(tid.x); }
[shader("compute")] [numthreads(WG_X, 1, 1)] void l8(uint3 tid : SV_DispatchThreadID) { gto_int_accumulate_xyz<8>(tid.x); }

// Point charges, conversion, and out
[shader("compute")]
[numthreads(WG_X, 1, 1)]
void resolve(uint3 tid : SV_DispatchThreadID) {
    const uint i = tid.x;
    if (i >= args.xyz_count) return;
    const float3 C = gto_int_point(i);
    uint2 acc = args.accum[i];
    for (uint k = 0u; k < args.num_points; ++k) {
        const float4 q = args.points[k];
        const float3 d = C - q.xyz;
        const float r2 = dot(d, d);
        if (r2 > 1.0e-20f) acc = gto_int_fixed_add(acc, q.w * rsqrt(r2));
    }
    args.out[args.xyz_first + i] = gto_int_fixed_to_float(acc);
}
