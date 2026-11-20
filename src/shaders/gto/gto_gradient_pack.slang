// gto_gradient_pack.slang
//
// md_gto_gpu_gradient_pack: the direction of a gradient (RGBA32F, xyz) packed into 32 bits as RGB10A2 is, x in the low
// bits (GL_UNSIGNED_INT_2_10_10_10_REV): round((c / 2 + 1 / 2) 1023) per component, alpha 3 where the gradient is
// not 0, else 0 (and the components 512).

#include "../md_gpu.slang"

struct RootArgs {
    uint4 grid_dim;
    DescriptorHandle<RWTexture3D<float4>> src;
    DescriptorHandle<RWTexture3D<uint>>   dst;
};

struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;

#define args (*md_root.args_ptr)

[shader("compute")]
[numthreads(8, 8, 4)]
void main(uint3 id : SV_DispatchThreadID) {
    if (!all(id < args.grid_dim.xyz)) return;
    const float3 g  = args.src[id].xyz;
    const float  m  = max(abs(g.x), max(abs(g.y), abs(g.z)));     // scaled by it first: |g|^2 neither under- nor overflows
    const bool   ok = m > 0.0f && m < 3.4e38f;                      // finite and not 0
    const float3 n  = ok ? normalize(g / m) : float3(0.0f);
    const uint3  q  = uint3(clamp(round((n * 0.5f + 0.5f) * 1023.0f), 0.0f, 1023.0f));
    args.dst[id] = q.x | (q.y << 10u) | (q.z << 20u) | ((ok ? 3u : 0u) << 30u);
}
