// gpu_imgui_test.slang -- the ImGui renderer's shaders, as an application
// would write them against md_gpu.h. Used by the ImGui test in test_gpu.c.

#include "../../src/shaders/md_gpu.slang"

// Scalars on purpose: {float2, float2, uint} would be 20 bytes in SPIR-V and
// 24 in MSL, and indexing through the pointer would drift. This is 20 bytes on
// both targets, the same as ImDrawVert.
struct Vert { float x, y, u, v; uint col; };

struct Args {
    float2 scale;
    float2 translate;
    Vert*  verts;
    DescriptorHandle<Texture2D>    tex;
    DescriptorHandle<SamplerState> smp;
};
MD_SHADER_ARGS(Args);

struct VOut {
    float4 pos : SV_Position;
    float2 uv  : UV;
    float4 col : COL;
};

[shader("vertex")]
VOut imgui_vs(uint vid : SV_VertexID) {       // vid = the index value
    Args a = MD_ARGS;
    Vert v = a.verts[vid];
    VOut o;
    o.pos = float4(float2(v.x, v.y) * a.scale + a.translate, 0, 1);
    o.uv  = float2(v.u, v.v);
    o.col = unpackUnorm4x8ToFloat(v.col);
    return o;
}

[shader("fragment")]
float4 imgui_fs(VOut i) : SV_Target0 {
    Args a = MD_ARGS;
    return i.col * a.tex.Sample(a.smp, i.uv);
}
