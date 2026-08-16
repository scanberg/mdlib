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
