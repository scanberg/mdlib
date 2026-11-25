// md_post_ui.slang -- md_post_ui's (md_post.h): the triangles of a UI's draw
// lists over the target, in order, one draw per command. The vertices are
// pulled from the list's (20 bytes each, read as words so that the layout is
// the caller's to the byte: ImDrawVert's), each index offset by the command's
// vertex offset through the address. The colour times the texture's sample,
// straight alpha, blended over the target (SRC_ALPHA, ONE_MINUS_SRC_ALPHA;
// alpha ONE, ONE_MINUS_SRC_ALPHA). Layout rules of md_gpu.slang. Mirrored in
// md_post.c.

#include "../md_gpu.slang"

struct Args {
    float4 xf;                  // pixels to clip: x * xf.x + xf.z, y * xf.y + xf.w
    uint*  vertices;            // the command's first vertex: 5 words each
    DescriptorHandle<Texture2D<float4>> tex;
    DescriptorHandle<SamplerState>      smp;
};
MD_SHADER_ARGS(Args);

struct UV {
    float4 pos : SV_Position;
    float2 uv  : TEXCOORD0;
    float4 col : COLOR0;
};

[shader("vertex")]
UV ui_vs(uint vid : SV_VertexID) {
    Args a = MD_ARGS;
    uint* v = a.vertices + vid * 5u;
    float2 p = asfloat(uint2(v[0], v[1]));
    uint   c = v[4];
    UV o;
    o.pos = float4(p.x * a.xf.x + a.xf.z, p.y * a.xf.y + a.xf.w, 0.0, 1.0);
    o.uv  = asfloat(uint2(v[2], v[3]));
    o.col = float4(float(c & 0xFFu), float((c >> 8) & 0xFFu), float((c >> 16) & 0xFFu), float(c >> 24)) / 255.0;
    return o;
}

[shader("fragment")]
float4 ui_fs(UV i) : SV_Target0 {
    Args a = MD_ARGS;
    return i.col * a.tex.Sample(a.smp, i.uv);
}
