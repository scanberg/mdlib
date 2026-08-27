// gpu_test.slang — kernels used by unittest/test_gpu.c.
//
// Each entry point is compiled separately. They deliberately cover the parts
// of md_gpu that cannot be tested from the host side alone: device pointers
// as kernel arguments, program-order dependencies between launches, bindless
// textures, and indirect launches driven by device-side counts.

#include "../../src/shaders/md_gpu.slang"

// ---------------------------------------------------------------------------
// fill: dst[i] = base + i
// ---------------------------------------------------------------------------
struct FillArgs {
    uint   n;
    uint   base;
    uint   _pad0;
    uint   _pad1;
    uint*  dst;
};

struct FillRoot { FillArgs* args; };
[[vk::push_constant]] ConstantBuffer<FillRoot> fill_root;

[shader("compute")]
[numthreads(64, 1, 1)]
void fill(uint3 tid : SV_DispatchThreadID) {
    FillArgs a = *fill_root.args;
    if (tid.x >= a.n) return;
    a.dst[tid.x] = a.base + tid.x;
}

// ---------------------------------------------------------------------------
// scale_add: dst[i] = src[i] * mul + add
// Chained after fill to prove program order gives the barrier for free.
// ---------------------------------------------------------------------------
struct ScaleArgs {
    uint   n;
    uint   mul;
    uint   add;
    uint   _pad;
    uint*  src;
    uint*  dst;
};

struct ScaleRoot { ScaleArgs* args; };
[[vk::push_constant]] ConstantBuffer<ScaleRoot> scale_root;

[shader("compute")]
[numthreads(64, 1, 1)]
void scale_add(uint3 tid : SV_DispatchThreadID) {
    ScaleArgs a = *scale_root.args;
    if (tid.x >= a.n) return;
    a.dst[tid.x] = a.src[tid.x] * a.mul + a.add;
}

// ---------------------------------------------------------------------------
// sum_reduce: single group, atomically accumulates src into *out.
// Used to produce a device-side count for the indirect-launch test.
// ---------------------------------------------------------------------------
struct SumArgs {
    uint   n;
    uint   _pad0;
    uint   _pad1;
    uint   _pad2;
    uint*  src;
    uint*  out_val;
};

struct SumRoot { SumArgs* args; };
[[vk::push_constant]] ConstantBuffer<SumRoot> sum_root;

[shader("compute")]
[numthreads(64, 1, 1)]
void sum_reduce(uint3 tid : SV_DispatchThreadID) {
    SumArgs a = *sum_root.args;
    if (tid.x >= a.n) return;
    uint prev;
    InterlockedAdd(a.out_val[0], a.src[tid.x], prev);
}

// ---------------------------------------------------------------------------
// tex_write: writes a ramp into a texture reached through a DescriptorHandle.
// Note there is no global array and no binding: the handle is a typed value in
// the argument struct.
// ---------------------------------------------------------------------------
struct TexArgs {
    uint   dim_x;
    uint   dim_y;
    uint   dim_z;
    float  scale;
    DescriptorHandle<RWTexture3D<float>> tex;
};

struct TexRoot { TexArgs* args; };
[[vk::push_constant]] ConstantBuffer<TexRoot> tex_root;

[shader("compute")]
[numthreads(4, 4, 4)]
void tex_write(uint3 tid : SV_DispatchThreadID) {
    TexArgs a = *tex_root.args;
    uint3 dim = uint3(a.dim_x, a.dim_y, a.dim_z);
    if (any(tid >= dim)) return;
    uint linear = tid.z * dim.x * dim.y + tid.y * dim.x + tid.x;
    a.tex[tid] = float(linear) * a.scale;
}

// ---------------------------------------------------------------------------
// tex_read: copies a bindless texture back into a buffer.
// ---------------------------------------------------------------------------
struct TexReadArgs {
    uint    dim_x;
    uint    dim_y;
    uint    dim_z;
    uint    _pad;
    DescriptorHandle<RWTexture3D<float>> tex;
    float*  dst;
};

struct TexReadRoot { TexReadArgs* args; };
[[vk::push_constant]] ConstantBuffer<TexReadRoot> tex_read_root;

[shader("compute")]
[numthreads(4, 4, 4)]
void tex_read(uint3 tid : SV_DispatchThreadID) {
    TexReadArgs a = *tex_read_root.args;
    uint3 dim = uint3(a.dim_x, a.dim_y, a.dim_z);
    if (any(tid >= dim)) return;
    uint linear = tid.z * dim.x * dim.y + tid.y * dim.x + tid.x;
    a.dst[linear] = a.tex[tid];
}

// ---------------------------------------------------------------------------
// layout_probe: writes every argument field back to a buffer.
//
// md_gpu.h's md_gpu_float4 / md_gpu_uint2 / ... exist to make one C struct match
// both SPIR-V's scalar layout and MSL's vector layout, and nothing in the suite
// exercised them. Reading the fields back turns a layout mismatch into "field X
// came through as Y" rather than a wrong result three launches later.
//
// The layout is chosen to be a real test of the ABI rule rather than a trivial
// one: `pair` sits at offset 24, which the two targets agree on only because 24
// is a multiple of 8 -- drop `dim_y` and SPIR-V would put it at 20 while MSL
// kept it at 24, which is exactly what check_gpu_arg_layout.py rejects. The
// three scalar dimensions stand in for the 3-vector the ABI does not allow, and
// the trailing pointer proves the struct's total footprint agrees too.
// ---------------------------------------------------------------------------
struct LayoutProbeArgs {
    float4  v4;         // @0
    uint    dim_x;      // @16
    uint    dim_y;      // @20
    uint2   pair;       // @24  -- only agrees because 24 % 8 == 0
    uint    dim_z;      // @32
    float   scale;      // @36
    uint*   dst;        // @40
};

struct LayoutProbeRoot { LayoutProbeArgs* args; };
[[vk::push_constant]] ConstantBuffer<LayoutProbeRoot> layout_probe_root;

[shader("compute")]
[numthreads(1, 1, 1)]
void layout_probe(uint3 tid : SV_DispatchThreadID) {
    LayoutProbeArgs a = *layout_probe_root.args;
    if (tid.x != 0) return;
    a.dst[0] = a.dim_x;
    a.dst[1] = a.dim_y;
    a.dst[2] = a.dim_z;
    a.dst[3] = asuint(a.scale);
    a.dst[4] = a.pair.x;
    a.dst[5] = a.pair.y;
    a.dst[6] = asuint(a.v4.x);
    a.dst[7] = asuint(a.v4.y);
    a.dst[8] = asuint(a.v4.z);
    a.dst[9] = asuint(a.v4.w);
}

// ---------------------------------------------------------------------------
// tex_probe: separates "did the argument struct arrive intact" from "did the
// texture handle reach the shader".
//
// A DescriptorHandle sits *inside* the struct, and the two backends reach it by
// completely different routes -- a heap index on Vulkan, a resource id on
// Metal. If a backend's toolchain lifts the handle out of the struct into its
// own binding slot, every member after it shifts by 8 bytes and the C mirror is
// silently wrong. `marker` sits after the handle for exactly that reason: it is
// a plain uint, so a shifted read is a wrong value rather than a wild pointer
// dereference, and the test reports instead of crashing.
//
// dst deliberately sits *before* the handle so it survives such a shift and the
// shader can always report what it saw.
// ---------------------------------------------------------------------------
struct TexProbeArgs {
    uint    n;
    uint    _pad0;
    uint    _pad1;
    uint    _pad2;
    uint*   dst;
    DescriptorHandle<RWTexture3D<float>> tex;
    uint    marker;
    uint    _pad3;
};

struct TexProbeRoot { TexProbeArgs* args; };
[[vk::push_constant]] ConstantBuffer<TexProbeRoot> tex_probe_root;

[shader("compute")]
[numthreads(1, 1, 1)]
void tex_probe(uint3 tid : SV_DispatchThreadID) {
    TexProbeArgs a = *tex_probe_root.args;
    if (tid.x != 0) return;
    a.dst[0] = a.marker;      // struct layout survived the handle
    a.dst[1] = a.n;
    a.tex[uint3(0, 0, 0)] = 1.0f;   // the handle actually resolved
    a.tex[uint3(1, 0, 0)] = 2.0f;
}

// ---------------------------------------------------------------------------
// bump: dst[i] += delta. Launched indirectly, so the number of groups comes
// from a device-side count written by an earlier launch.
// ---------------------------------------------------------------------------
struct BumpArgs {
    uint   n;
    uint   delta;
    uint   _pad0;
    uint   _pad1;
    uint*  dst;
};

struct BumpRoot { BumpArgs* args; };
[[vk::push_constant]] ConstantBuffer<BumpRoot> bump_root;

[shader("compute")]
[numthreads(64, 1, 1)]
void bump(uint3 tid : SV_DispatchThreadID) {
    BumpArgs a = *bump_root.args;
    if (tid.x >= a.n) return;
    a.dst[tid.x] = a.dst[tid.x] + a.delta;
}
