// gpu_volume.slang -- kernels for examples/gpu_volume.c
//
// A miniature of the shape viamd actually runs: evaluate a scalar field into a
// 3D volume, compact the voxels above a threshold into a list with a
// device-side counter, then process exactly that many with an indirect launch.
//
// Vectors in argument structs are fine -- mirror them on the host with
// md_gpu_uint4 / md_gpu_float4 / ..., which carry the alignment the ABI
// requires. The one restriction is that a 3-vector may not appear in an
// argument struct (12 bytes in SPIR-V, 16 in MSL), so a volume dimension
// travels as a uint4 and is read back as .xyz.
// tools/check_gpu_arg_layout.py enforces both at build time.

#include "../../src/shaders/md_gpu.slang"

// ---------------------------------------------------------------------------
// 1. eval_field -- write an analytic scalar field into a 3D volume.
//    Stands in for md_gto's density evaluation.
// ---------------------------------------------------------------------------
struct EvalArgs {
    uint4 dim;         // xyz used
    float freq;
    uint  _pad;        // handles need an 8-byte offset, same as a vec2
    DescriptorHandle<RWTexture3D<float>> vol;
};
struct EvalRoot { EvalArgs* args; };
[[vk::push_constant]] ConstantBuffer<EvalRoot> eval_root;

[shader("compute")]
[numthreads(4, 4, 4)]
void eval_field(uint3 tid : SV_DispatchThreadID) {
    EvalArgs a = *eval_root.args;
    uint3 dim = a.dim.xyz;
    if (any(tid >= dim)) return;

    float3 p = (float3(tid) + 0.5f) / float3(dim) * 2.0f - 1.0f;
    float  r = length(p);
    a.vol[tid] = exp(-a.freq * r * r);
}

// ---------------------------------------------------------------------------
// 2. compact_above -- append the linear index of every voxel above `threshold`
//    to `indices`, bumping a device-side counter. The counter never touches
//    the CPU; it drives the indirect launch below.
// ---------------------------------------------------------------------------
struct CompactArgs {
    uint4  dim;        // xyz used
    float  threshold;
    uint   _pad0;      // handles need an 8-byte offset, same as a vec2
    DescriptorHandle<RWTexture3D<float>> vol;
    uint*  count;      // single uint32
    uint*  indices;    // capacity assumed sufficient
    uint   capacity;
    uint   _pad1;
};
struct CompactRoot { CompactArgs* args; };
[[vk::push_constant]] ConstantBuffer<CompactRoot> compact_root;

[shader("compute")]
[numthreads(4, 4, 4)]
void compact_above(uint3 tid : SV_DispatchThreadID) {
    CompactArgs a = *compact_root.args;
    uint3 dim = a.dim.xyz;
    if (any(tid >= dim)) return;
    if (a.vol[tid] <= a.threshold) return;

    uint slot;
    InterlockedAdd(a.count[0], 1u, slot);
    if (slot < a.capacity) {
        a.indices[slot] = tid.z * dim.x * dim.y + tid.y * dim.x + tid.x;
    }
}

// ---------------------------------------------------------------------------
// 3. gather -- launched indirectly over exactly the number of compacted
//    voxels. Reads the volume back through the same bindless handle.
// ---------------------------------------------------------------------------
struct GatherArgs {
    uint4  dim;        // xyz used
    DescriptorHandle<RWTexture3D<float>> vol;
    uint*  count;
    uint*  indices;
    float* values;
};
struct GatherRoot { GatherArgs* args; };
[[vk::push_constant]] ConstantBuffer<GatherRoot> gather_root;

[shader("compute")]
[numthreads(64, 1, 1)]
void gather(uint3 tid : SV_DispatchThreadID) {
    GatherArgs a = *gather_root.args;
    if (tid.x >= a.count[0]) return;

    uint  li = a.indices[tid.x];
    uint3 c  = uint3(li % a.dim.x, (li / a.dim.x) % a.dim.y, li / (a.dim.x * a.dim.y));
    a.values[tid.x] = a.vol[c];
}
